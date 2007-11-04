/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "libgretl.h"
#include "gretl_string_table.h"
#include "libset.h"
#include "csvdata.h"

#include <ctype.h>
#include <glib.h>

#define QUOTE      '\''
#define CSVSTRLEN  72
#define NON_NUMERIC 1.0e99

enum {
    CSV_HAVEDATA = 1 << 0,
    CSV_GOTDELIM = 1 << 1,
    CSV_GOTTAB   = 1 << 2,
    CSV_BLANK1   = 1 << 3,
    CSV_OBS1     = 1 << 4,
    CSV_TRAIL    = 1 << 5,
    CSV_AUTONAME = 1 << 6,
    CSV_NONNUM   = 1 << 7
};

typedef struct csvdata_ csvdata;

struct csvdata_ {
    int flags;
    char delim;
    int markerpd;
    int real_n;
    double **Z;
    DATAINFO *dinfo;
    int ncols, nrows;
    char str[CSVSTRLEN];
    char skipstr[8];
    int *codelist;
    char *descrip;
    gretl_string_table *st;
};

#define csv_has_trailing_comma(c) (c->flags & CSV_TRAIL)
#define csv_has_obs_column(c)     (c->flags & CSV_OBS1)
#define csv_has_blank_column(c)   (c->flags & CSV_BLANK1)
#define csv_got_tab(c)            (c->flags & CSV_GOTTAB)
#define csv_got_delim(c)          (c->flags & CSV_GOTDELIM)
#define csv_autoname(c)           (c->flags & CSV_AUTONAME)
#define csv_skip_column(c)        (c->flags & (CSV_OBS1 | CSV_BLANK1))
#define csv_have_data(c)          (c->flags & CSV_HAVEDATA)
#define csv_force_nonnum(c)       (c->flags & CSV_NONNUM)

#define csv_set_trailing_comma(c)   (c->flags |= CSV_TRAIL)
#define csv_unset_trailing_comma(c) (c->flags &= ~CSV_TRAIL)
#define csv_set_obs_column(c)       (c->flags |= CSV_OBS1)
#define csv_set_blank_column(c)     (c->flags |= CSV_BLANK1)
#define csv_set_got_tab(c)          (c->flags |= CSV_GOTTAB)
#define csv_set_got_delim(c)        (c->flags |= CSV_GOTDELIM)
#define csv_set_autoname(c)         (c->flags |= CSV_AUTONAME)
#define csv_set_force_nonnum(c)     (c->flags |= CSV_NONNUM)

#define csv_skipping(c)        (*c->skipstr != '\0')
#define csv_has_non_numeric(c) (c->st != NULL)

static int 
time_series_label_check (DATAINFO *pdinfo, char *skipstr, PRN *prn);

static void csvdata_free (csvdata *c)
{
    if (c->descrip != NULL) {
	free(c->descrip);
    }

    if (c->st != NULL) {
	gretl_string_table_destroy(c->st);
    }

    if (c->codelist != NULL) {
	free(c->codelist);
    }

    destroy_dataset(c->Z, c->dinfo);

    free(c);
}

static csvdata *csvdata_new (double **Z, const DATAINFO *pdinfo,
			     gretlopt opt)
{
    csvdata *c = malloc(sizeof *c);

    if (c == NULL) {
	return NULL;
    }

    c->flags = (opt & OPT_C)? CSV_NONNUM : 0;
    c->delim = '\t';
    c->markerpd = -1;
    c->real_n = 0;
    c->Z = NULL;
    c->dinfo = NULL;
    c->ncols = 0;
    c->nrows = 0;
    *c->str = '\0';
    *c->skipstr = '\0';
    c->codelist = NULL;
    c->descrip = NULL;
    c->st = NULL;

    c->dinfo = datainfo_new();
    if (c->dinfo == NULL) {
	free(c);
	c = NULL;
    } else {
	if (pdinfo != NULL) {
	    c->delim = c->dinfo->delim = pdinfo->delim;
	}
	c->dinfo->delim = c->delim;
	if (Z != NULL) {
	    c->flags |= CSV_HAVEDATA;
	}
    }

    return c;
}

static int add_obs_marker (DATAINFO *pdinfo, int n)
{
    char **S;

    S = realloc(pdinfo->S, n * sizeof *S);
    if (S == NULL) {
	return 1;
    }

    pdinfo->S = S;

    pdinfo->S[n-1] = malloc(OBSLEN);
    if (pdinfo->S[n-1] == NULL) {
	return 1;
    }

    strcpy(pdinfo->S[n-1], "NA");

    return 0;
}

static int add_single_obs (double ***pZ, DATAINFO *pdinfo)
{
    double *x;
    int i, err = 0;

    for (i=0; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    x = realloc((*pZ)[i], (pdinfo->n + 1) * sizeof *x);
	    if (x != NULL) {
		(*pZ)[i] = x;
	    } else {
		return 1;
	    }
	}
    }

    pdinfo->n += 1;

    (*pZ)[0][pdinfo->n - 1] = 1.0;

    for (i=1; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    (*pZ)[i][pdinfo->n - 1] = NADBL;
	}
    }

    if (pdinfo->S != NULL) {
	err = add_obs_marker(pdinfo, pdinfo->n);
    }

    return err;
}

static int pad_weekly_data (double ***pZ, DATAINFO *pdinfo, int add)
{
    int oldn = pdinfo->n;
    int ttarg, offset = 0, skip = 0;
    int i, s, t, tc, err;

    err = dataset_add_observations(add, pZ, pdinfo, OPT_A); 

    if (!err) {
	for (t=0; t<oldn; t++) {
	    tc = calendar_obs_number(pdinfo->S[t], pdinfo) - offset;
	    if (tc != t) {
		skip = tc - t;
		fprintf(stderr, "Gap of size %d at original t = %d\n", skip, t);
		offset += skip;
		ttarg = oldn - 1 + offset;
		for (s=0; s<oldn-t+skip; s++) {
		    for (i=1; i<pdinfo->v; i++) {
			if (s < oldn - t) {
			    if (s == 0 || s == oldn-t-1) {
				fprintf(stderr, "shifting obs %d to obs %d\n",
					ttarg-skip, ttarg);
			    }
			    (*pZ)[i][ttarg] = (*pZ)[i][ttarg - skip];
			} else {
			    fprintf(stderr, "inserting NA at obs %d\n", ttarg);
			    (*pZ)[i][ttarg] = NADBL;
			}
		    }
		    ttarg--;
		}
	    }
	}
    }

    return err;
}

/* FIXME the following needs to be made more flexible */

static int csv_weekly_data (double ***pZ, DATAINFO *pdinfo)
{
    char *lbl2 = pdinfo->S[pdinfo->n - 1];
    int ret = 1;
    int misscount = 0;
    int t, tc;

    for (t=0; t<pdinfo->n; t++) {
	tc = calendar_obs_number(pdinfo->S[t], pdinfo) - misscount;
	if (tc != t) {
	    misscount += tc - t;
	}
    }

    if (misscount > 0) {
	double missfrac = (double) misscount / pdinfo->n;

	fprintf(stderr, "nobs = %d, misscount = %d (%.2f%%)\n", 
		pdinfo->n, misscount, 100.0 * missfrac);
	if (missfrac > 0.05) {
	    ret = 0;
	} else {
	    int Tc = calendar_obs_number(lbl2, pdinfo) + 1;
	    int altmiss = Tc - pdinfo->n;

	    fprintf(stderr, "check: Tc = %d, missing = %d\n", Tc, altmiss);
	    if (altmiss != misscount) {
		ret = 0;
	    } else {
		int err;

		fprintf(stderr, "OK, consistent\n");
		err = pad_weekly_data(pZ, pdinfo, misscount);
		if (err) ret = 0;
	    } 
	} 
    }
    
    return ret;
}

#define DAY_DEBUG 0

static int check_daily_dates (DATAINFO *pdinfo, int *pd, PRN *prn)
{
    int T = pdinfo->n;
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[T - 1];
    int fulln = 0, n, t, nbak;
    int oldpd = pdinfo->pd;
    double oldsd0 = pdinfo->sd0;
    long ed1, ed2;
    int nmiss = 0, err = 0;

    *pd = 0;
    
    ed1 = get_epoch_day(lbl1);
    if (ed1 < 0) {
	err = 1;
    }

#if DAY_DEBUG    
    fprintf(stderr, "S[0] = '%s', ed1 = %ld\n", lbl1, ed1);
#endif

    pdinfo->pd = guess_daily_pd(pdinfo);
    pdinfo->structure = TIME_SERIES;

#if DAY_DEBUG    
    fprintf(stderr, "guessed daily pd = %d\n", pdinfo->pd);
#endif

    if (!err) {
	ed2 = get_epoch_day(lbl2);
	if (ed2 <= ed1) {
	    err = 1;
	} else {
	    pdinfo->sd0 = ed1;
	}
    }

    if (!err) {
	int n1 = calendar_obs_number(lbl1, pdinfo);
	int n2 = calendar_obs_number(lbl2, pdinfo);

	fulln = n2 - n1 + 1;
	if (T > fulln) {
	    err = 1;
	} else {
	    nmiss = fulln - T;
	    pprintf(prn, "Observations: %d; days in sample: %d\n", 
		    T, fulln);
	    if (nmiss > 300 * T) {
		pprintf(prn, "Probably annual data\n");
		*pd = 1;
	    } else if (nmiss > 50 * T) {
		pprintf(prn, "Probably quarterly data\n");
		*pd = 4;
	    } else if (nmiss > 20 * T) {
		pprintf(prn, "Probably monthly data\n");
		*pd = 12;
	    } else if (nmiss > 5 * T) {
		pprintf(prn, "Probably weekly data\n");
		*pd = pdinfo->pd = 52;
	    } else {
		pprintf(prn, "Missing daily observations: %d\n", nmiss);
	    }
	}
    }

    nbak = 0;
    for (t=0; t<pdinfo->n && !err; t++) {
	n = calendar_obs_number(pdinfo->S[t], pdinfo);
	if (n < t) {
	    pprintf(prn, "Daily dates error at t = %d:\n"
		    "  calendar_obs_number() for '%s' = %d but t = %d\n", 
		    t, pdinfo->S[t], n, t);
	    err = 1;
	} else if (n > fulln - 1) {
	    pprintf(prn, "Error: date '%s' out of bounds\n", pdinfo->S[t]);
	    err = 1;
	} else if (nbak > 0 && n == nbak) {
	    pprintf(prn, "Error: date '%s' is repeated\n", pdinfo->S[t]);
	    err = 1;
	}
	nbak = n;
    }

    if (err) {
	pdinfo->pd = oldpd;
	pdinfo->sd0 = oldsd0;
	pdinfo->structure = CROSS_SECTION;
    } else {
	strcpy(pdinfo->stobs, lbl1);
	strcpy(pdinfo->endobs, lbl2);
	pdinfo->t2 = pdinfo->n - 1;
	if (nmiss > 0 && *pd == 0) {
	    pdinfo->markers = DAILY_DATE_STRINGS;
	}
    }

#if DAY_DEBUG
    fprintf(stderr, "check_daily_dates: pd = %d, err = %d\n", 
	    pdinfo->pd, err);
#endif

    return (err)? -1 : pdinfo->pd;
}

/* There's a special case (ugh!) where observation strings are
   given as in monthly data, but the frequency is in fact
   quarterly, as in:

   1947.06
   1947.09
   1947.12
   1948.03 

   we'll make a brave attempt to handle this.
*/ 

#define fakequarter(m) (m==3 || m==6 || m==9 || m==12) 

static int complete_qm_labels (DATAINFO *pdinfo, char *skipstr,
			       int *ppd, const char *fmt,
			       PRN *prn)
{
    char bad[16], skip[8];
    int t, y, p, Ey, Ep;
    int pmin = 1;
    int pd, pd0;
    int ret = 1;

    pd = pd0 = *ppd;

 restart:

    if (sscanf(pdinfo->S[0], fmt, &y, &p) != 2) {
	return 0;
    }

    for (t=1; t<pdinfo->n; t++) {
	Ey = (p == pd)? y + 1 : y;
	Ep = (p == pd)? pmin : p + pmin;
	if (sscanf(pdinfo->S[t], fmt, &y, &p) != 2) {
	    ret = 0;
	} else if (Ep == 1 && pd == pd0 && p == pd + 1 
		   && skipstr != NULL) {
	    *skip = *bad = '\0';
	    strncat(skip, pdinfo->S[t] + 4, 7); 
	    strncat(bad, pdinfo->S[t], 15); 
	    pd = pd0 + 1;
	    goto restart;
	} else if (p == Ep + 2 && pmin == 1 && fakequarter(p)) {
	    *bad = '\0';
	    strncat(bad, pdinfo->S[t], 15); 
	    pmin = 3;
	    goto restart;
	} else if (y != Ey || p != Ep) {
	    ret = 0;
	}
	if (!ret) {
	    pprintf(prn, "   %s: not a consistent date\n", 
		    pdinfo->S[t]);
	    break;
	}
    }

    if (ret) {
	if (pmin == 3) {
	    pprintf(prn, "   \"%s\": quarterly data pretending to be monthly?\n", 
		    bad);
	    *ppd = 4;
	} else if (pd == pd0 + 1) {
	    pprintf(prn, "   \"%s\": BLS-type nonsense? Trying again\n", 
		    bad);
	    strcpy(skipstr, skip);
	}
    }

    return ret;
}

static int complete_year_labels (const DATAINFO *pdinfo)
{
    int t, yr, yrbak = atoi(pdinfo->S[0]);
    int ret = 1;

    for (t=1; t<pdinfo->n; t++) {
	yr = atoi(pdinfo->S[t]);
	if (yr != yrbak + 1) {
	    ret = 0;
	    break;
	}
	yrbak = yr;
    }

    return ret;
}

static int compress_daily (DATAINFO *pdinfo, int pd)
{
    int t, yr, mon, day;

    for (t=0; t<pdinfo->n; t++) {
	sscanf(pdinfo->S[t], "%d/%d/%d", &yr, &mon, &day);
	if (pd == 1) {
	    sprintf(pdinfo->S[t], "%d", yr);
	} else if (pd == 12) {
	    sprintf(pdinfo->S[t], "%d:%02d", yr, mon);
	} else if (pd == 4) {
	    sprintf(pdinfo->S[t], "%d:%d", yr, mon / 3 + (mon % 3 != 0));
	} 
    }

    return 0;
}

enum date_orders {
    YYYYMMDD = 1,
    MMDDYYYY,
    DDMMYYYY
};

int get_date_order (int f0, int fn) 
{
    if (f0 > 31 || fn > 31) {
	/* first field must be year */
	return YYYYMMDD;
    } else if (f0 > 12 || fn > 12) {
	/* first field must be day */
	return DDMMYYYY;
    } else {
	return MMDDYYYY;
    }
}

static int transform_daily_dates (DATAINFO *pdinfo, int dorder)
{
    int t, yr, mon, day;
    char sep1[2], sep2[2];
    int sret, err = 0;

    for (t=0; t<pdinfo->n && !err; t++) {
	if (dorder == YYYYMMDD) {
	    sret = sscanf(pdinfo->S[t], "%d%1[/-]%d%1[/-]%d", 
			  &yr, sep1, &mon, sep2, &day);
	} else if (dorder == DDMMYYYY) {
	    sret = sscanf(pdinfo->S[t], "%d%1[/-]%d%1[/-]%d", 
			  &day, sep1, &mon, sep2, &yr);
	} else {
	    sret = sscanf(pdinfo->S[t], "%d%1[/-]%d%1[/-]%d", 
			  &mon, sep1, &day, sep2, &yr);
	}
	if (sret == 5) {
	    sprintf(pdinfo->S[t], "%02d/%02d/%02d", yr, mon, day);
	} else {
	    err = 1;
	}
    }

    return err;
}

static int 
csv_daily_date_check (double ***pZ, DATAINFO *pdinfo, char *skipstr, 
		      PRN *prn)
{
    int d1[3], d2[3];
    char sep1[2], sep2[2];
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[pdinfo->n - 1];

    if (pZ == NULL) {
	return E_DATA;
    }

    if (sscanf(lbl1, "%d%1[/-]%d%1[/-]%d", &d1[0], sep1, &d1[1], sep2, &d1[2]) == 5 &&
	sscanf(lbl2, "%d%1[/-]%d%1[/-]%d", &d2[0], sep1, &d2[1], sep2, &d2[2]) == 5) {
	int yr1, mon1, day1;
	int yr2, mon2, day2;
	int dorder = get_date_order(d1[0], d2[0]);
	int pd, ret = 0;

	if (dorder == YYYYMMDD) {
	    pputs(prn, "Trying date order YYYYMMDD\n");
	    yr1 = d1[0];
	    mon1 = d1[1];
	    day1 = d1[2];
	    yr2 = d2[0];
	    mon2 = d2[1];
	    day2 = d2[2];
	} else if (dorder == DDMMYYYY) {
	    pputs(prn, "Trying date order DDMMYYYY\n");
	    day1 = d1[0];
	    mon1 = d1[1];
	    yr1 = d1[2];
	    day2 = d2[0];
	    mon2 = d2[1];
	    yr2 = d2[2];
	} else {
	    pputs(prn, "Trying date order MMDDYYYY\n");
	    mon1 = d1[0];
	    day1 = d1[1];
	    yr1 = d1[2];
	    mon2 = d2[0];
	    day2 = d2[1];
	    yr2 = d2[2];
	}		
	    
	if (yr2 >= yr1 && 
	    mon1 > 0 && mon1 < 13 &&
	    mon2 > 0 && mon2 < 13 && 
	    day1 > 0 && day1 < 32 &&
	    day2 > 0 && day2 < 32) {
	    /* looks promising for calendar dates */
	    if (dorder != YYYYMMDD || *sep1 != '/' || *sep2 != '/') {
		if (transform_daily_dates(pdinfo, dorder)) {
		    return -1;
		}
	    }
	    pprintf(prn, "? %s - %s\n", lbl1, lbl2);
	    ret = check_daily_dates(pdinfo, &pd, prn);
	    if (ret >= 0 && pd > 0) {
		if (pd == 52) {
		    if (csv_weekly_data(pZ, pdinfo)) {
			ret = 52;
		    } else {
			ret = -1;
		    }
		} else {
		    compress_daily(pdinfo, pd);
		    ret = time_series_label_check(pdinfo, 
						  skipstr, 
						  prn);
		}
	    } 
	    return ret;
	}
    } 

    return -1;
}

static int pd_from_date_label (const char *lbl, char *year, char *subp,
			       char *format, PRN *prn)
{
    const char *subchars = ".:QqMmPp";
    int len = strlen(lbl);
    int try, pd = -1;

    *year = '\0';
    strncat(year, lbl, 4);
    try = atoi(year);

    if (try > 0 && try < 3000) {
	pprintf(prn, M_("   %s: probably a year... "), year);
    } else {
	pprintf(prn, M_("   %s: probably not a year\n"), year);
    }

    if (len == 5) {
	pputs(prn, M_("   but I can't make sense of the extra bit\n"));
    } else if (len == 4) {
	pputs(prn, M_("and just a year\n"));
	pd = 1;
    } else {
	char sep = lbl[4];
	char sub[3], *s = NULL;
	int p;

	if (strchr(subchars, sep)) {
	    *sub = '\0';
	    strncat(sub, lbl + 5, 2);
	    s = sub;
	    if (len == 6 || (len == 7 && (sep == 'q' || sep == 'Q'))) {
		if (len == 7) s++;
		p = atoi(s);
		if (p > 0 && p < 5) {
		    pprintf(prn, M_("quarter %s?\n"), s);
		    pd = 4;
		} else {
		    pprintf(prn, "quarter %d: not possible\n", p);
		}
	    } else if (len == 7) {
		p = atoi(s);
		if (p > 0 && p < 13) {
		    pprintf(prn, M_("month %s?\n"), s);
		    pd = 12;
		} else {
		    pprintf(prn, "month %d: not possible\n", p);
		}
	    }
	    strcpy(subp, s);
	    if (format != NULL && (pd == 4 || pd == 12)) {
		sprintf(format, "%%d%c%%d", sep);
	    }
	}
    }

    return pd;
}

static int 
time_series_label_check (DATAINFO *pdinfo, char *skipstr, PRN *prn)
{
    char year[5], sub[3];
    char format[8] = {0};
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[pdinfo->n - 1];
    int pd = -1;

    pd = pd_from_date_label(lbl1, year, sub, format, prn);

    if (pd == 1) {
	if (complete_year_labels(pdinfo)) {
	    pdinfo->pd = pd;
	    strcpy(pdinfo->stobs, year);
	    pdinfo->sd0 = atof(pdinfo->stobs);
	    strcpy(pdinfo->endobs, lbl2);
	} else {
	    pputs(prn, M_("   but the dates are not complete and consistent\n"));
	    pd = -1;
	}
    } else if (pd == 4 || pd == 12) {
	int savepd = pd;

	if (complete_qm_labels(pdinfo, skipstr, &pd, format, prn)) {
	    pdinfo->pd = pd;
	    if (savepd == 12 && pd == 4) {
		int s = atoi(sub) / 3;

		sprintf(pdinfo->stobs, "%s:%d", year, s);
	    } else {
		sprintf(pdinfo->stobs, "%s:%s", year, sub);
	    }
	    pdinfo->sd0 = obs_str_to_double(pdinfo->stobs);
	    ntodate(pdinfo->endobs, pdinfo->n - 1, pdinfo);
	} else {
	    pputs(prn, M_("   but the dates are not complete and consistent\n"));
	    pd = -1;
	}
    }

    return pd;
}

/* attempt to parse csv row labels as dates.  Return -1 if this
   doesn't work out, or 0 if the labels seem to be just integer
   observation numbers, else return the inferred data frequency 
*/

int test_markers_for_dates (double ***pZ, DATAINFO *pdinfo, 
			    char *skipstr, PRN *prn)
{
    char endobs[OBSLEN];
    int n = pdinfo->n;
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[n - 1];
    int len1 = strlen(lbl1);

    if (skipstr != NULL && *skipstr != '\0') {
	return time_series_label_check(pdinfo, skipstr, prn);
    }

    pprintf(prn, M_("   first row label \"%s\", last label \"%s\"\n"), 
	    lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", n);
    if (!strcmp(lbl1, "1") && !strcmp(lbl2, endobs)) {
	return 0;
    }

    /* labels are of different lengths? */
    if (len1 != strlen(lbl2)) {
	pputs(prn, M_("   label strings can't be consistent dates\n"));
	return -1;
    }

    pputs(prn, M_("trying to parse row labels as dates...\n"));

    if (len1 == 8 || len1 == 10) {
	/* daily data? */
	return csv_daily_date_check(pZ, pdinfo, skipstr, prn);
    } else if (len1 >= 4) {
	/* annual, quarterly, monthly? */
	if (isdigit((unsigned char) lbl1[0]) &&
	    isdigit((unsigned char) lbl1[1]) &&
	    isdigit((unsigned char) lbl1[2]) && 
	    isdigit((unsigned char) lbl1[3])) {
	    return time_series_label_check(pdinfo, skipstr, prn);
	} else {
	    pputs(prn, M_("   definitely not a four-digit year\n"));
	}
    }

    return -1;
}

/* The function below checks for the maximum line length in the given
   file.  It also checks for extraneous binary data (the file is 
   supposed to be plain text), and checks whether the 'delim'
   character is present in the file, on a non-comment line (where
   a comment line is one that starts with '#').  

   In addition, we check whether the file has a trailing comma on every
   line.
*/

static int csv_max_line_length (FILE *fp, csvdata *cdata, PRN *prn)
{
    int c, c1, cbak = 0, cc = 0;
    int comment = 0, maxlen = 0;

    csv_set_trailing_comma(cdata);

    while ((c = fgetc(fp)) != EOF) {
	if (c == 0x0d) {
	    /* CR */
	    c1 = fgetc(fp);
	    if (c1 == EOF) {
		break;
	    } else if (c1 == 0x0a) {
		/* CR + LF -> LF */
		c = c1;
	    } else {
		/* Mac-style: CR not followed by LF */
		c = 0x0a;
		ungetc(c1, fp);
	    }
	}
	if (c == 0x0a) {
	    if (cc > maxlen) {
		maxlen = cc;
	    }
	    cc = 0;
	    if (cbak != 0 && cbak != ',') {
		csv_unset_trailing_comma(cdata);
	    }
	    continue;
	}
	cbak = c;
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ)) {
	    pprintf(prn, M_("Binary data (%d) encountered: this is not a valid "
			   "text file\n"), c);
	    return -1;
	}
	if (cc == 0) {
	    comment = (c == '#');
	}
	if (!comment) {
	    if (c == '\t') {
		csv_set_got_tab(cdata);
	    }
	    if (c == cdata->delim) {
		csv_set_got_delim(cdata);
	    }
	}
	cc++;
    }

    if (maxlen == 0) {
	pprintf(prn, M_("Data file is empty\n"));
    } else if (csv_has_trailing_comma(cdata)) {
	pprintf(prn, M_("Data file has trailing commas\n"));
    }

    if (maxlen > 0) {
	/* allow for newline and null terminator */
	maxlen += 3;
    }

    return maxlen;
}

static int count_csv_fields (const char *line, char delim)
{
    int cbak, nf = 0;
    const char *p = line;

    if (*p == delim && *p == ' ') p++;

    while (*p) {
	if (*p == delim) nf++;
	cbak = *p;
	p++;
	/* Problem: (when) should trailing delimiter be read as implicit "NA"? */
	if (*p == '\0' && cbak == delim && cbak != ',') {
	    nf--;
	}
    }

    return nf + 1;
}

static void remove_quoted_commas (char *s)
{
    int inquote = 0;

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	}
	if (inquote && *s == ',') {
	    *s = ' ';
	}
	s++;
    }
}

static void compress_csv_line (char *line, csvdata *c)
{
    int n = strlen(line);
    char *p = line + n - 1;

    if (*p == '\n') {
	*p = '\0';
	p--;
    }

    if (*p == '\r') *p = '\0';

    if (c->delim == ',') {
	remove_quoted_commas(line);
    }

    if (c->delim != ' ') {
	delchar(' ', line);
    } else {
	compress_spaces(line);
    }

    delchar('"', line);

    if (csv_has_trailing_comma(c)) {
	/* chop trailing comma */
	n = strlen(line);
	if (n > 0) {
	    line[n-1] = '\0';
	}
    }
}

int import_obs_label (const char *s)
{
    char tmp[32];

    *tmp = '\0';
    strncat(tmp, s, 31);
    lower(tmp);

    return (!strcmp(s, "obs") ||
	    !strcmp(s, "date") || 
	    !strcmp(s, "year") || 
	    !strcmp(s, "period"));    
}

static void check_first_field (const char *line, csvdata *c, PRN *prn)
{
    if (c->delim != ' ' && *line == c->delim) {
	csv_set_blank_column(c);
    } else {
	char field1[16];
	int i = 0;

	if (c->delim == ' ' && *line == ' ') line++;

	while (*line && i < 15) {
	    if (*line == c->delim) break;
	    field1[i++] = *line++;
	}
	field1[i] = '\0';
	iso_to_ascii(field1);

	pprintf(prn, M_("   first field: '%s'\n"), field1);

	if (import_obs_label(field1)) {
	    pputs(prn, M_("   seems to be observation label\n"));
	    csv_set_obs_column(c);
	}
    }
}

static int csv_missval (const char *str, int i, int t, PRN *prn)
{
    int miss = 0;

    if (*str == '\0') {
	if (t < 100) {
	    pprintf(prn, M_("   the cell for variable %d, obs %d "
			    "is empty: treating as missing value\n"), 
		    i, t);
	}
	miss = 1;
    }

    if (import_na_string(str)) {
	if (t < 100) {
	    pprintf(prn, M_("   warning: missing value for variable "
			    "%d, obs %d\n"), i, t);
	}
	miss = 1;
    }

    return miss;
}

static int non_numeric_check (csvdata *c, PRN *prn)
{
    int *list = NULL;
    int i, j, t, nn = 0;
    int err = 0;

    for (i=1; i<c->dinfo->v; i++) {
	if (is_codevar(c->dinfo->varname[i])) {
	    nn++;
	} else {
	    for (t=0; t<c->dinfo->n; t++) {
		if (c->Z[i][t] == NON_NUMERIC) {
		    nn++;
		    break;
		}
	    }
	}
    }

    if (nn > 0) {
	list = gretl_list_new(nn);
	if (list == NULL) {
	    err = E_ALLOC;
	}
    }

    if (list != NULL) {
	j = 1;
	for (i=1; i<c->dinfo->v; i++) {
	    if (is_codevar(c->dinfo->varname[i])) {
		list[j++] = i;
	    } else {
		for (t=0; t<c->dinfo->n; t++) {
		    if (c->Z[i][t] == NON_NUMERIC) {
			list[j++] = i;
			break;
		    }
		}
	    }
	}

	printlist(list, "non-numeric vars list");

	for (i=1; i<=list[0]; i++) {
	    int v = list[i];
	    int cv = is_codevar(c->dinfo->varname[v]);

	    set_var_discrete(c->dinfo, v, 1);

	    if (!cv) {
		double nnfrac;
		int nnon = 0;
		int nok = 0;
		int tn = 0;

		for (t=0; t<c->dinfo->n; t++) {
		    if (c->Z[v][t] == NON_NUMERIC) {
			if (!tn) tn = t + 1;
			nnon++;
		    } else if (!na(c->Z[v][t])) {
			nok++;
		    }
		}

		nnfrac = (nok == 0)? 1.0 : (double) nnon / nok;
		pprintf(prn, "variable %d (%s): non-numeric values = %d "
			"(%.2f percent)\n", v, c->dinfo->varname[v], 
			nnon, 100 * nnfrac);
		if (!csv_force_nonnum(c) && (nnon < 2 || nnfrac < 0.01)) {
		    pprintf(prn, M_("ERROR: variable %d (%s), observation %d, "
				    "non-numeric value\n"), 
			    v, c->dinfo->varname[v], tn);
		    err = E_DATA;
		}
	    }
	}

	if (!err) {
	    pputs(prn, "allocating string table\n");
	    c->st = string_table_new_from_cols_list(list);
	    if (c->st == NULL) {
		err = E_ALLOC;
		free(list);
	    } else {
		c->codelist = list;
	    }
	}
    }

    return err;
}

static int process_csv_obs (csvdata *c, int i, int t, 
			    gretlopt opt, PRN *prn)
{
    int ix, err = 0;

    if (c->st != NULL) {
	if (in_gretl_list(c->codelist, i) && !na(c->Z[i][t])) {
	    ix = gretl_string_table_index(c->st, c->str, i, 0, prn);
	    if (ix > 0) {
		c->Z[i][t] = (double) ix;
	    } else {
		err = E_DATA;
	    }
	}
    } else if (csv_missval(c->str, i, t+1, prn)) {
	c->Z[i][t] = NADBL;
    } else if (check_atof(c->str)) {
	c->Z[i][t] = NON_NUMERIC;
    } else {
	c->Z[i][t] = atof(c->str);
    }

    return err;
}

/* wrapper for fgets(), designed to handle any sort of line
   termination (unix, DOS, Mac or an unholy mixture)
*/

static char *csv_fgets (char *s, int n, FILE *fp)
{
    int i, c1, c = 0;

    for (i=0; i<n-1 && c!=0x0a; i++) {
	c = fgetc(fp);
	if (c == EOF) {
	    if (i == 0) {
		/* signal end of read */
		return NULL;
	    } else {
		break;
	    }
	} else if (c == 0x0d) {
	    /* CR: convert to LF and peek at next char: if it's
	       LF swallow it, otherwise put it back */
	    c = 0x0a;
	    c1 = fgetc(fp);
	    if (c1 != 0x0a) {
		ungetc(c1, fp);
	    }
	}
	s[i] = c;
    }

    s[i] = '\0';

    return s;
}

/* pick up any comments following the data block in a CSV file */

static char *get_csv_descrip (char *line, int n, FILE *fp)
{
    char *desc = NULL;

    while (csv_fgets(line, n, fp)) {
	tailstrip(line);
	if (desc == NULL) {
	    desc = malloc(strlen(line) + 2);
	    if (desc == NULL) {
		return NULL;
	    }
	    sprintf(desc, "%s\n", line);
	} else {
	    char *tmp;

	    tmp = realloc(desc, strlen(desc) + strlen(line) + 2);
	    if (tmp == NULL) {
		free(desc);
		return NULL;
	    }
	    desc = tmp;
	    strcat(desc, line);
	    strcat(desc, "\n");
	}
    }

    if (desc != NULL && string_is_blank(desc)) {
	free(desc);
	desc = NULL;
    }

    return desc;
}

static const char *
csv_msg = N_("\nPlease note:\n"
	     "- The first row of the CSV file should contain the "
	     "names of the variables.\n"
	     "- The first column may optionally contain date "
	     "strings or other 'markers':\n  in that case its row 1 entry "
	     "should be blank, or should say 'obs' or 'date'.\n"
	     "- The remainder of the file must be a rectangular "
	     "array of data.\n");

static int csv_fields_check (char *line, int maxlen, FILE *fp,
			     csvdata *c, PRN *prn)
{
    int gotdata = 0;
    int chkcols = 0;
    int err = 0;

    c->ncols = c->nrows = 0;

    while (csv_fgets(line, maxlen, fp) && !err) {

	/* skip comment lines */
	if (*line == '#') {
	    continue;
	}

	/* skip blank lines -- but finish if the blank comes after data */
	if (string_is_blank(line)) {
	    if (gotdata) {
		if (!csv_have_data(c)) {
		    c->descrip = get_csv_descrip(line, maxlen, fp);
		}
		break;
	    } else {
		continue;
	    }
	}
	
	c->nrows += 1;
	compress_csv_line(line, c);

	if (!gotdata) {
	    /* scrutinize first "real" line */
	    check_first_field(line, c, prn);
	    gotdata = 1;
	} 

	chkcols = count_csv_fields(line, c->delim);
	if (c->ncols == 0) {
	    c->ncols = chkcols;
	    pprintf(prn, M_("   number of columns = %d\n"), c->ncols);	    
	} else if (chkcols != c->ncols) {
	    pprintf(prn, M_("   ...but row %d has %d fields: aborting\n"),
		    c->nrows, chkcols);
	    pputs(prn, M_(csv_msg));
	    err = E_DATA;
	}
    }

    return err;
}

static int 
csv_reconfigure_for_markers (double ***pZ, DATAINFO *pdinfo)
{
    if (dataset_allocate_obs_markers(pdinfo)) {
	return 1;
    }

    return dataset_drop_last_variables(1, pZ, pdinfo);
}

#define starts_number(c) (isdigit((unsigned char) c) || c == '-' || \
                          c == '+' || c == '.')

#define obs_labels_no_varnames(o,c,n)  (!o && c->v > 3 && n == c->v - 2)

static int csv_varname_scan (csvdata *c, char *line, int maxlen, FILE *fp,
			     PRN *prn, PRN *mprn)
{
    char *p;
    int obscol = csv_has_obs_column(c);
    int i, k, numcount;
    int err = 0;

    pputs(mprn, M_("scanning for variable names...\n"));

    while (csv_fgets(line, maxlen, fp)) {
	if (*line == '#' || string_is_blank(line)) {
	    continue;
	} else {
	    break;
	}
    }

    compress_csv_line(line, c);   

    p = line;
    if (c->delim == ' ' && *p == ' ') p++;
    iso_to_ascii(p);
    pprintf(mprn, M_("   line: %s\n"), p);
    
    numcount = 0;

    for (k=0; k<c->ncols && !err; k++) {
	int nv = 0;

	i = 0;
	while (*p && *p != c->delim) {
	    if (i < CSVSTRLEN - 1) {
		c->str[i++] = *p;
	    }
	    p++;
	}
	if (*p == c->delim) p++;

	c->str[i] = 0;

	if (k == 0 && (csv_skip_column(c))) {
	    ;
	} else {
	    nv = (csv_skip_column(c))? k : k + 1;

	    if (*c->str == '\0') {
		pprintf(prn, M_("   variable name %d is missing: aborting\n"), nv);
		pputs(prn, M_(csv_msg));
		err = E_DATA;
	    } else {
		c->dinfo->varname[nv][0] = 0;
		strncat(c->dinfo->varname[nv], c->str, VNAMELEN - 1);
		if (starts_number(*c->str)) {
		    numcount++;
		} else {
		    iso_to_ascii(c->dinfo->varname[nv]);
		    if (check_varname(c->dinfo->varname[nv])) {
			errmsg(1, prn);
			err = E_DATA;
		    }
		}
	    }
	}
	if (nv == c->dinfo->v - 1) break;
    }

    if (err) {
	return err;
    }

    if (numcount == c->dinfo->v - 1 || 
	obs_labels_no_varnames(obscol, c->dinfo, numcount)) {
	pputs(prn, M_("it seems there are no variable names\n"));
	/* then we undercounted the observations by one */
	if (add_single_obs(&c->Z, c->dinfo)) {
	    err = E_ALLOC;
	} else {
	    csv_set_autoname(c);
	    rewind(fp);
	    if (obs_labels_no_varnames(obscol, c->dinfo, numcount)) {
		if (csv_reconfigure_for_markers(&c->Z, c->dinfo)) {
		    err = E_ALLOC;
		} else {
		    csv_set_obs_column(c);
		}
	    }
	}
    } else if (numcount > 0) {
	for (i=1; i<c->dinfo->v; i++) {
	    if (check_varname(c->dinfo->varname[i])) {
		errmsg(1, prn);
		break;
	    }
	}
	err = E_DATA;
    }

    return err;
}

static int 
real_read_labels_and_data (csvdata *c, char *line, int maxlen, 
			   FILE *fp, gretlopt opt, PRN *prn)
{
    char *p;
    int i, k, nv, t = 0;
    int err = 0;

    c->real_n = c->dinfo->n;

    while (csv_fgets(line, maxlen, fp) && !err) {

	if (*line == '#' || string_is_blank(line)) {
	    continue;
	}

	if (csv_skipping(c) && strstr(line, c->skipstr)) {
	    c->real_n -= 1;
	    continue;
	}

	compress_csv_line(line, c);
	p = line;
	if (c->delim == ' ' && *p == ' ') p++;

	for (k=0; k<c->ncols; k++) {
	    i = 0;
	    while (*p && *p != c->delim) {
		if (i < CSVSTRLEN - 1) {
		    c->str[i++] = *p;
		} else {
		    pprintf(prn, M_("warning: truncating data at row %d, column %d\n"),
			    t+1, k+1);
		}
		p++;
	    }
	    if (*p == c->delim) {
		p++;
	    }
	    c->str[i] = '\0';
	    if (k == 0 && csv_skip_column(c) && c->dinfo->S != NULL) {
		char *S = c->str;

		c->dinfo->S[t][0] = 0;
		if (*S == '"' || *S == '\'') {
		    S++;
		}
		strncat(c->dinfo->S[t], S, OBSLEN - 1);
		iso_to_ascii(c->dinfo->S[t]);
	    } else {
		nv = (csv_skip_column(c))? k : k + 1;
		err = process_csv_obs(c, nv, t, opt, prn);
		if (err) {
		    break;
		}
	    }
	}

	if (++t == c->dinfo->n) {
	    break;
	}
    }

    if (!err && c->real_n < c->dinfo->n) {
	int drop = c->dinfo->n - c->real_n;

	err = dataset_drop_observations(drop, &c->Z, c->dinfo);
    }

    return err;
}

static int csv_read_data (csvdata *c, char *line, int maxlen,
			  FILE *fp, gretlopt opt, PRN *prn,
			  PRN *mprn)
{
    int err;

    pputs(mprn, M_("scanning for row labels and data...\n"));
    err = real_read_labels_and_data(c, line, maxlen, fp, opt, prn);

    if (!err && csv_skip_column(c)) {
	c->markerpd = test_markers_for_dates(&c->Z, c->dinfo, c->skipstr, 
					     prn);
    }

    return err;
}

static void csv_parsing_header (const char *fname, PRN *prn)
{
    if (!g_utf8_validate(fname, -1, NULL)) {
	gchar *trfname = g_locale_to_utf8(fname, -1, NULL, NULL, NULL);

	pprintf(prn, "%s %s...\n", M_("parsing"), trfname);
	g_free(trfname);
    } else {
	pprintf(prn, "%s %s...\n", M_("parsing"), fname);
    }
}

/**
 * import_csv:
 * @pZ: pointer to data set.
 * @ppdinfo: pointer to data information struct.
 * @fname: name of CSV file.
 * @opt: not yet documented.
 * @prn: gretl printing struct (can be NULL).
 * 
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int import_csv (double ***pZ, DATAINFO **ppdinfo, 
		const char *fname, gretlopt opt, PRN *prn)
{
    csvdata *c = NULL;
    int popit = 0;
    int i, maxlen;
    FILE *fp = NULL;
    PRN *mprn = NULL;
    char *line = NULL;
    long datapos;
    int err = 0;

#ifdef ENABLE_NLS
    if (prn != NULL) {
	check_for_console(prn);
    }
#endif

    if (gretl_messages_on()) {
	mprn = prn;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto csv_bailout;
    }

    c = csvdata_new(*pZ, *ppdinfo, opt);
    if (c == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }	

    if (mprn != NULL) {
	csv_parsing_header(fname, mprn);
    }

    /* get line length, also check for binary data, etc. */
    maxlen = csv_max_line_length(fp, c, prn);    
    if (maxlen <= 0) {
	err = E_DATA;
	goto csv_bailout;
    } 

    if (!csv_got_delim(c)) {
	/* set default delimiter */
	if (csv_got_tab(c)) {
	    c->delim = c->dinfo->delim = '\t';
	} else {
	    c->delim = c->dinfo->delim = ' ';
	}
    }

    pprintf(mprn, M_("using delimiter '%c'\n"), c->delim);
    pprintf(mprn, M_("   longest line: %d characters\n"), maxlen - 1);

    if (csv_has_trailing_comma(c) && c->delim != ',') {
	csv_unset_trailing_comma(c);
    }

    /* buffer to hold lines */
    line = malloc(maxlen);
    if (line == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }  
    
    rewind(fp);

    /* read lines, check for consistency in number of fields */
    err = csv_fields_check(line, maxlen, fp, c, prn);
    if (err) {
	goto csv_bailout;
    }

    c->dinfo->n = c->nrows - 1; /* allow for var headings */
    c->dinfo->v = (csv_skip_column(c))? c->ncols : c->ncols + 1;

    pprintf(mprn, M_("   number of variables: %d\n"), c->dinfo->v - 1);
    pprintf(mprn, M_("   number of non-blank lines: %d\n"), c->nrows);

    if (c->dinfo->n == 0) {
	pputs(prn, M_("Invalid data file\n"));
	err = E_DATA;
	goto csv_bailout;
    }

    /* initialize datainfo and Z */
    if (start_new_Z(&c->Z, c->dinfo, 0)) {
	err = E_ALLOC;
	goto csv_bailout;
    }

    if (csv_skip_column(c)) {
	if (dataset_allocate_obs_markers(c->dinfo)) {
	    err = E_ALLOC;
	    goto csv_bailout;
	}
    }

    /* second pass */

    rewind(fp);

    err = csv_varname_scan(c, line, maxlen, fp, prn, mprn);
    if (err) {
	goto csv_bailout;
    }

    if (*ppdinfo != NULL && (*ppdinfo)->decpoint != ',') {
	gretl_push_c_numeric_locale();
	popit = 1;
    }

    datapos = ftell(fp);

    err = csv_read_data(c, line, maxlen, fp, opt, prn, mprn);

    if (!err && csv_skipping(c)) {
	/* try again */
	fseek(fp, datapos, SEEK_SET);
	err = csv_read_data(c, line, maxlen, fp, opt, prn, NULL);
    }

    if (!err) {
	err = non_numeric_check(c, prn);
	if (!err && csv_has_non_numeric(c)) {
	    /* try once more */
	    fseek(fp, datapos, SEEK_SET);
	    err = csv_read_data(c, line, maxlen, fp, opt, prn, NULL);
	}
    }	

    if (popit) {
	gretl_pop_c_numeric_locale();
    }

    if (err) {
	goto csv_bailout;
    }

    c->dinfo->t1 = 0;
    c->dinfo->t2 = c->dinfo->n - 1;

    if (c->markerpd > 0) {
	pputs(mprn, M_("taking date information from row labels\n\n"));
	if (csv_skipping(c)) {
	    pprintf(prn, "WARNING: Check your data! gretl has stripped out "
		    "what appear to be\nextraneous lines in a %s dataset: " 
		    "this may not be right.\n\n",
		    (c->dinfo->pd == 4)? "quarterly" : "monthly");
	}
    } else {
	pputs(mprn, M_("treating these as undated data\n\n"));
	dataset_obs_info_default(c->dinfo);
    }

    if (c->dinfo->pd != 1 || strcmp(c->dinfo->stobs, "1")) { 
        c->dinfo->structure = TIME_SERIES;
    }

    if (c->st != NULL) {
	gretl_string_table_print(c->st, c->dinfo, fname, prn);
    }

    /* If there were observation labels and they were not interpretable
       as dates, and they weren't simply "1, 2, 3, ...", then they 
       should probably be preserved; otherwise discard them. 
    */
    if (c->dinfo->S != NULL && c->markerpd >= 0 && 
	c->dinfo->markers != DAILY_DATE_STRINGS) {
	dataset_destroy_obs_markers(c->dinfo);
    }

    if (csv_autoname(c)) {
	/* no variable names were found */
	for (i=1; i<c->dinfo->v; i++) {
	    sprintf(c->dinfo->varname[i], "v%d", i);
	}
    } else if (fix_varname_duplicates(c->dinfo)) {
	pputs(prn, M_("warning: some variable names were duplicated\n"));
    }

    if (*pZ == NULL) {
	/* no dataset currently in place */
	*pZ = c->Z;
	c->Z = NULL;
	if (*ppdinfo != NULL) {
	    free(*ppdinfo);
	}
	if (c->descrip != NULL) {
	    c->dinfo->descrip = c->descrip;
	    c->descrip = NULL;
	}
	*ppdinfo = c->dinfo;
	c->dinfo = NULL;
    } else {
	err = merge_data(pZ, *ppdinfo, c->Z, c->dinfo, prn);
	if (err) {
	    /* merge_data() frees the dataset */
	    c->Z = NULL;
	    c->dinfo = NULL;
	    goto csv_bailout;
	}
    }

 csv_bailout:

    if (fp != NULL) {
	fclose(fp);
    }

    free(line);
    csvdata_free(c);

    if (err == E_ALLOC) {
	pputs(prn, M_("Out of memory\n"));
    }    

#ifdef ENABLE_NLS
    console_off();
#endif

    return err;
}
