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
#include <errno.h>
#include <glib.h>

#define QUOTE      '\''
#define CSVSTRLEN  72
#define NON_NUMERIC 1.0e99

enum {
    CSV_HAVEDATA = 1 << 0,
    CSV_GOTDELIM = 1 << 1,
    CSV_GOTTAB   = 1 << 2,
    CSV_GOTSEMI  = 1 << 3,
    CSV_BLANK1   = 1 << 4,
    CSV_OBS1     = 1 << 5,
    CSV_TRAIL    = 1 << 6,
    CSV_AUTONAME = 1 << 7,
    CSV_NONNUM   = 1 << 8,
    CSV_REVERSED = 1 << 9
};

typedef struct csvdata_ csvdata;

struct csvdata_ {
    int flags;
    char delim;
    int markerpd;
    int maxlen;
    int real_n;
    char *line;
    double **Z;
    DATAINFO *dinfo;
    int ncols, nrows;
    char str[CSVSTRLEN];
    char skipstr[8];
    int *codelist;
    char *descrip;
    gretl_string_table *st;
    int *cols_list;
    int *width_list;
};

#define csv_has_trailing_comma(c) (c->flags & CSV_TRAIL)
#define csv_has_obs_column(c)     (c->flags & CSV_OBS1)
#define csv_has_blank_column(c)   (c->flags & CSV_BLANK1)
#define csv_got_tab(c)            (c->flags & CSV_GOTTAB)
#define csv_got_semi(c)           (c->flags & CSV_GOTSEMI)
#define csv_got_delim(c)          (c->flags & CSV_GOTDELIM)
#define csv_autoname(c)           (c->flags & CSV_AUTONAME)
#define csv_skip_column(c)        (c->flags & (CSV_OBS1 | CSV_BLANK1))
#define csv_have_data(c)          (c->flags & CSV_HAVEDATA)
#define csv_force_nonnum(c)       (c->flags & CSV_NONNUM)
#define csv_data_reversed(c)      (c->flags & CSV_REVERSED)

#define csv_set_trailing_comma(c)   (c->flags |= CSV_TRAIL)
#define csv_unset_trailing_comma(c) (c->flags &= ~CSV_TRAIL)
#define csv_set_obs_column(c)       (c->flags |= CSV_OBS1)
#define csv_set_blank_column(c)     (c->flags |= CSV_BLANK1)
#define csv_set_got_tab(c)          (c->flags |= CSV_GOTTAB)
#define csv_set_got_semi(c)         (c->flags |= CSV_GOTSEMI)
#define csv_set_got_delim(c)        (c->flags |= CSV_GOTDELIM)
#define csv_set_autoname(c)         (c->flags |= CSV_AUTONAME)
#define csv_set_force_nonnum(c)     (c->flags |= CSV_NONNUM)
#define csv_set_data_reversed(c)    (c->flags |= CSV_REVERSED)

#define csv_skipping(c)        (*c->skipstr != '\0')
#define csv_has_non_numeric(c) (c->st != NULL)

static int 
time_series_label_check (DATAINFO *pdinfo, int reversed, char *skipstr, PRN *prn);

static void csvdata_free (csvdata *c)
{
    if (c == NULL) {
	return;
    }

    if (c->descrip != NULL) {
	free(c->descrip);
    }

    if (c->st != NULL) {
	gretl_string_table_destroy(c->st);
    }

    if (c->codelist != NULL) {
	free(c->codelist);
    }

    if (c->line != NULL) {
	free(c->line);
    }  

    if (c->cols_list != NULL) {
	free(c->cols_list);
	free(c->width_list);
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
    c->maxlen = 0;
    c->real_n = 0;
    c->line = NULL;
    c->Z = NULL;
    c->dinfo = NULL;
    c->ncols = 0;
    c->nrows = 0;
    *c->str = '\0';
    *c->skipstr = '\0';
    c->codelist = NULL;
    c->descrip = NULL;
    c->st = NULL;
    c->cols_list = NULL;
    c->width_list = NULL;

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

static int csvdata_add_cols_list (csvdata *c, const char *s)
{
    int *list, *clist = NULL, *wlist = NULL;
    int i, n, m, err = 0;

    list = gretl_list_from_string(s, &err);

    if (!err) {
	n = list[0];
	if (n == 0 || n % 2 != 0) {
	    err = E_DATA;
	} else {
	    m = n / 2;
	    clist = gretl_list_new(m);
	    wlist = gretl_list_new(m);
	    if (clist == NULL || wlist == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	int j = 1;

	for (i=1; i<=n; i+=2, j++) {
	    clist[j] = list[i];
	    wlist[j] = list[i+1];
	}

	/* clist = column start list: must be a set of increasing
	   positive integers; and wlist = respective column widths,
	   must all be positive 
	*/
	for (i=1; i<=m; i++) {
	    if (wlist[i] <= 0 || clist[i] <= 0 || 
		(i > 1 && clist[i] <= clist[i-1])) {
		err = E_DATA;
		break;
	    } else if (wlist[i] >= CSVSTRLEN) {
		fprintf(stderr, "Warning: field %d too wide (%d), truncating\n", 
			i, wlist[i]);
		wlist[i] = CSVSTRLEN - 1;
	    }
	}
    }

    free(list);

    if (!err) {
	c->cols_list = clist;
	c->width_list = wlist;
    } else {
	free(clist);
	free(wlist);
	if (err == E_DATA) {
	    gretl_errmsg_set(_("Invalid column specification"));
	}
    }

    return err;
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
	x = realloc((*pZ)[i], (pdinfo->n + 1) * sizeof *x);
	if (x != NULL) {
	    (*pZ)[i] = x;
	} else {
	    return 1;
	}
    }

    pdinfo->n += 1;

    (*pZ)[0][pdinfo->n - 1] = 1.0;

    for (i=1; i<pdinfo->v; i++) {
	(*pZ)[i][pdinfo->n - 1] = NADBL;
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
	    } else if (pZ != NULL) {
		int err;

		fprintf(stderr, "OK, consistent\n");
		err = pad_weekly_data(pZ, pdinfo, misscount);
		if (err) ret = 0;
	    } 
	} 
    }
    
    return ret;
}

#define DAY_DEBUG 1

static int check_daily_dates (DATAINFO *pdinfo, int *pd, int *reversed, PRN *prn)
{
    int T = pdinfo->n;
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[T - 1];
    int fulln = 0, n, t, nbak;
    int alt_pd = 0;
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
	if (ed2 < 0) {
	    err = 1;
	} else if (ed2 < ed1) {
#if DAY_DEBUG    
	    fprintf(stderr, "check_daily_dates: data are reversed?\n");
#endif
	    pdinfo->sd0 = ed2;
	    *reversed = 1;
	} else {
	    pdinfo->sd0 = ed1;
	}
    }

 recompute:

    alt_pd = 0;
    nbak = 0;

    if (!err) {
	int n1 = calendar_obs_number((*reversed)? lbl2 : lbl1, pdinfo);
	int n2 = calendar_obs_number((*reversed)? lbl1 : lbl2, pdinfo);

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
	    } else if (nmiss > 3 * T) {
		pprintf(prn, "Probably weekly data\n");
		*pd = pdinfo->pd = 52;
	    } else {
		pprintf(prn, "Missing daily observations: %d\n", nmiss);
	    }
	}
    }

    nbak = 0;

    for (t=0; t<pdinfo->n && !err; t++) {
	int wd, s = (*reversed)? (pdinfo->n - 1 - t) : t;

	wd = get_day_of_week(pdinfo->S[s]);

	if (pdinfo->pd == 5 && (wd == 6 || wd == 0)) {
	    /* Got Sat or Sun, can't be 5-day daily? */
	    alt_pd = (wd == 6)? 6 : 7;
	    pprintf(prn, "Found a Saturday (%s): re-trying with pd = %d\n", 
		    pdinfo->S[s], alt_pd);
	    break;
	} else if (pdinfo->pd == 6 && wd == 0) {
	    /* Got Sun, can't be 6-day daily? */
	    alt_pd = 7;
	    pprintf(prn, "Found a Sunday (%s): re-trying with pd = %d\n", 
		    pdinfo->S[s], alt_pd);
	    break;
	}
	    
	n = calendar_obs_number(pdinfo->S[s], pdinfo);
	if (n < t) {
	    pprintf(prn, "Daily dates error at t = %d:\n"
		    "  calendar_obs_number() for '%s' = %d but t = %d\n", 
		    t, pdinfo->S[t], n, t);
	    err = 1;
	} else if (n > fulln - 1) {
	    pprintf(prn, "Error: date '%s' out of bounds\n", pdinfo->S[s]);
	    err = 1;
	} else if (nbak > 0 && n == nbak) {
	    pprintf(prn, "Error: date '%s' is repeated\n", pdinfo->S[s]);
	    err = 1;
	}
	nbak = n;
    }

    if (alt_pd > 0) {
	pdinfo->pd = alt_pd;
	goto recompute;
    }

    if (err) {
	pdinfo->pd = oldpd;
	pdinfo->sd0 = oldsd0;
	pdinfo->structure = CROSS_SECTION;
    } else {
	strcpy(pdinfo->stobs, (*reversed)? lbl2 : lbl1);
	strcpy(pdinfo->endobs, (*reversed)? lbl1 : lbl2);
	pdinfo->t2 = pdinfo->n - 1;
	if (nmiss > 0 && *pd == 0) {
	    pdinfo->markers = DAILY_DATE_STRINGS;
	}
    }

#if DAY_DEBUG
    fprintf(stderr, "check_daily_dates: pd = %d, reversed = %d, err = %d\n", 
	    pdinfo->pd, *reversed, err);
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

static int complete_qm_labels (DATAINFO *pdinfo, int reversed,
			       char *skipstr, int *ppd, 
			       const char *fmt,
			       PRN *prn)
{
    char bad[16], skip[8];
    int t, s, yr, per, Ey, Ep;
    int pmin = 1;
    int pd, pd0;
    int ret = 1;

    pd = pd0 = *ppd;

 restart:

    s = (reversed)? (pdinfo->n - 1) : 0;
    if (sscanf(pdinfo->S[s], fmt, &yr, &per) != 2) {
	return 0;
    }

    for (t=1; t<pdinfo->n; t++) {
	s = (reversed)? (pdinfo->n - 1 - t) : t;
	Ey = (per == pd)? yr + 1 : yr;
	Ep = (per == pd)? pmin : per + pmin;
	if (sscanf(pdinfo->S[s], fmt, &yr, &per) != 2) {
	    ret = 0;
	} else if (Ep == 1 && pd == pd0 && per == pd + 1 
		   && skipstr != NULL) {
	    *skip = *bad = '\0';
	    strncat(skip, pdinfo->S[s] + 4, 7); 
	    strncat(bad, pdinfo->S[s], 15); 
	    pd = pd0 + 1;
	    goto restart;
	} else if (per == Ep + 2 && pmin == 1 && fakequarter(per)) {
	    *bad = '\0';
	    strncat(bad, pdinfo->S[s], 15); 
	    pmin = 3;
	    goto restart;
	} else if (yr != Ey || per != Ep) {
	    ret = 0;
	}
	if (!ret) {
	    pprintf(prn, "   %s: not a consistent date\n", 
		    pdinfo->S[s]);
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

static int complete_year_labels (const DATAINFO *pdinfo, int reversed)
{
    int s, t, yr, yrbak;
    int ret = 1;

    s = (reversed)? (pdinfo->n - 1) : 0;
    yrbak = atoi(pdinfo->S[s]);

    for (t=1; t<pdinfo->n; t++) {
	s = (reversed)? (pdinfo->n - 1 - t) : t;
	yr = atoi(pdinfo->S[s]);
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
	/* could be wrong here */
	return MMDDYYYY;
    }
}

static void retransform_daily_dates (DATAINFO *pdinfo)
{
    int t, y, m, d;

    /* we apparently guessed wrongly at MMDDYYYY, so
       put the dates back as they were for another try,
       at DDMMYYYY.
    */

    for (t=0; t<pdinfo->n; t++) {
	sscanf(pdinfo->S[t], "%d/%d/%d", &y, &d, &m);
	sprintf(pdinfo->S[t], "%d/%d/%d", d, m, y);
    }
}

static int transform_daily_dates (DATAINFO *pdinfo, int dorder)
{
    char *label;
    int t, yr, mon, day;
    char s1, s2;
    int n, err = 0;

    for (t=0; t<pdinfo->n && !err; t++) {
	label = pdinfo->S[t];
	if (dorder == YYYYMMDD) {
	    n = sscanf(label, "%d%c%d%c%d", &yr, &s1, &mon, &s2, &day);
	} else if (dorder == DDMMYYYY) {
	    n = sscanf(label, "%d%c%d%c%d", &day, &s1, &mon, &s2, &yr);
	} else {
	    n = sscanf(label, "%d%c%d%c%d", &mon, &s1, &day, &s2, &yr);
	}
	if (n == 5) {
	    sprintf(label, "%02d/%02d/%02d", yr, mon, day);
	} else {
	    err = 1;
	}
    }

    return err;
}

void reverse_data (double **Z, DATAINFO *pdinfo, PRN *prn)
{
    char tmp[OBSLEN];
    double x;
    int T = pdinfo->n / 2;
    int i, t, s;

    pprintf(prn, M_("reversing the data!\n"));

    for (t=0; t<T; t++) {
	s = pdinfo->n - 1 - t;
	for (i=1; i<pdinfo->v; i++) {
	    x = Z[i][t];
	    Z[i][t] = Z[i][s];
	    Z[i][s] = x;
	}
	if (pdinfo->S != NULL) {
	    strcpy(tmp, pdinfo->S[t]);
	    strcpy(pdinfo->S[t], pdinfo->S[s]);
	    strcpy(pdinfo->S[s], tmp);
	}
    }
}

static int 
csv_daily_date_check (double ***pZ, DATAINFO *pdinfo, int *reversed,
		      char *skipstr, PRN *prn)
{
    int d1[3], d2[3];
    char s1, s2;
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[pdinfo->n - 1];
    int dorder = 0;

    if (sscanf(lbl1, "%d%c%d%c%d", &d1[0], &s1, &d1[1], &s2, &d1[2]) == 5 &&
	sscanf(lbl2, "%d%c%d%c%d", &d2[0], &s1, &d2[1], &s2, &d2[2]) == 5 &&
	s1 == s2 && ispunct(s1)) {
	int yr1, mon1, day1;
	int yr2, mon2, day2;
	int pd, ret = 0;

	dorder = get_date_order(d1[0], d2[0]);

    tryagain:

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
	    
	if (mon1 > 0 && mon1 < 13 &&
	    mon2 > 0 && mon2 < 13 && 
	    day1 > 0 && day1 < 32 &&
	    day2 > 0 && day2 < 32) {
	    /* looks promising for calendar dates */
	    if (dorder != YYYYMMDD || s1 != '/' || s2 != '/') {
		if (transform_daily_dates(pdinfo, dorder)) {
		    return -1;
		}
	    }
	    pprintf(prn, "Could be %s - %s\n", lbl1, lbl2);
	    ret = check_daily_dates(pdinfo, &pd, reversed, prn);
	    if (ret >= 0 && pd > 0) {
		if (pd == 52) {
		    if (csv_weekly_data(pZ, pdinfo)) {
			ret = 52;
		    } else if (dorder == MMDDYYYY) {
			/* maybe we guessed wrong */
			retransform_daily_dates(pdinfo);
			dorder = DDMMYYYY;
			goto tryagain;
		    } else {
			ret = -1;
		    }
		} else {
		    compress_daily(pdinfo, pd);
		    ret = time_series_label_check(pdinfo, 
						  *reversed,
						  skipstr, 
						  prn);
		}
	    } 
	    return ret;
	}
    } else {
	pprintf(prn, "'%s' and '%s': couldn't get dates\n", lbl1, lbl2);
    }

    return -1;
}

static int pd_from_date_label (const char *lbl, char *year, char *subp,
			       char *format, PRN *prn)
{
    const char *subchars = ".:QqMmPp-";
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

static int time_series_label_check (DATAINFO *pdinfo, int reversed,
				    char *skipstr, PRN *prn)
{
    char year[5], sub[3];
    char format[8] = {0};
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[pdinfo->n - 1];
    int pd = -1;

    pd = pd_from_date_label((reversed)? lbl2 : lbl1, year, sub, 
			    format, prn);

    if (pd == 1) {
	if (complete_year_labels(pdinfo, reversed)) {
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

	if (complete_qm_labels(pdinfo, reversed, skipstr, &pd, format, prn)) {
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

static int dates_maybe_reversed (const char *s1, const char *s2,
				 PRN *prn)
{
    char d1[5], d2[5];
    int ret = 0;

    *d1 = *d2 = '\0';

    strncat(d1, s1, 4);
    strncat(d2, s2, 4);

    ret = atoi(d1) > atoi(d2);

    if (ret) {
	pputs(prn, M_("   dates are reversed?\n"));
    }
    
    return ret;
}

/* e.g. "M1 1957", "M12 2009", with spaces removed */

static int fix_IFS_data_labels (DATAINFO *pdinfo)
{
    char *s1 = pdinfo->S[0];
    char *s2 = pdinfo->S[pdinfo->n - 1];
    int ret = 0;

    if ((*s1 == 'M' || *s1 == 'Q') && *s2 == *s1) {
	const char *dig = "0123456789";
	int n1 = strlen(s1);
	int n2 = strlen(s2);

	if ((n1 == 6 || n1 == 7) && (n2 == 6 || n2 == 7) &&
	    strspn(s1 + 1, dig) == n1 - 1 &&
	    strspn(s2 + 1, dig) == n2 - 1) {
	    char sp[3], tmp[8], *s;
	    int pmax = (*s1 == 'M')? 12 : 4;
	    int y, p, pbak = 0;
	    int i, n, doit = 1;

	    for (i=0; i<pdinfo->n; i++) {
		s = pdinfo->S[i];
		if (*s != *s1) {
		    doit = 0;
		    break;
		}
		n = strlen(s);
		if (n != 6 && n != 7) {
		    doit = 0;
		    break;
		}
		if (strspn(s + 1, dig) != n - 1) {
		    doit = 0;
		    break;
		}
		y = atoi(s + n - 4);
		*sp = '\0';
		strncat(sp, s + 1, n - 5);
		p = atoi(sp);
		if (y < 1800 || y > 2500 || p <= 0 || p > pmax) {
		    doit = 0;
		    break;
		} 
		if (i > 0 && p != pbak + 1 && p != 1) {
		    doit = 0;
		    break;
		}		    
		pbak = p;
	    }

	    if (doit) {
		for (i=0; i<pdinfo->n; i++) {
		    s = pdinfo->S[i];
		    n = strlen(s);
		    y = atoi(s + n - 4);
		    *sp = '\0';
		    strncat(sp, s + 1, n - 5);
		    p = atoi(sp);
		    if (pmax == 12) {
			sprintf(tmp, "%d:%02d", y, p);
		    } else {
			sprintf(tmp, "%d:%d", y, p);
		    }
		    if (strlen(tmp) > strlen(s)) {
			free(s);
			pdinfo->S[i] = gretl_strdup(tmp);
		    } else {
			strcpy(s, tmp);
		    }		    
		}
		ret = 1;
	    }
	}
    }

    return ret;
}

/* Attempt to parse csv row labels as dates.  Return -1 if this
   doesn't work out, or 0 if the labels seem to be just integer
   observation numbers, else return the inferred data frequency.  The
   pZ argument may be NULL; in that case we will not attempt to pad
   weekly data.
*/

int test_markers_for_dates (double ***pZ, DATAINFO *pdinfo, 
			    int *reversed, char *skipstr, 
			    PRN *prn)
{
    char endobs[OBSLEN];
    int n = pdinfo->n;
    char *lbl1 = pdinfo->S[0];
    char *lbl2 = pdinfo->S[n - 1];
    int len1 = strlen(lbl1);
    int pd = -1;

    if (skipstr != NULL && *skipstr != '\0') {
	return time_series_label_check(pdinfo, *reversed, skipstr, prn);
    }

    pprintf(prn, M_("   first row label \"%s\", last label \"%s\"\n"), 
	    lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", n);
    if (!strcmp(lbl1, "1") && !strcmp(lbl2, endobs)) {
	return 0;
    }

    if (fix_IFS_data_labels(pdinfo)) {
	lbl1 = pdinfo->S[0];
	lbl2 = pdinfo->S[n - 1];
	len1 = strlen(lbl1);
    }

    /* labels are of different lengths? */
    if (len1 != strlen(lbl2)) {
	pputs(prn, M_("   label strings can't be consistent dates\n"));
	return -1;
    }

    pputs(prn, M_("trying to parse row labels as dates...\n"));

    if (len1 == 8 || len1 == 10) {
	/* daily data? */
	pd = csv_daily_date_check(pZ, pdinfo, reversed, skipstr, prn);
    } else if (len1 >= 4) {
	/* annual, quarterly, monthly? */
	if (isdigit((unsigned char) lbl1[0]) &&
	    isdigit((unsigned char) lbl1[1]) &&
	    isdigit((unsigned char) lbl1[2]) && 
	    isdigit((unsigned char) lbl1[3])) {
	    *reversed = dates_maybe_reversed(lbl1, lbl2, prn);
	    pd = time_series_label_check(pdinfo, *reversed, skipstr, prn);
	} else {
	    pputs(prn, M_("   definitely not a four-digit year\n"));
	}
    }

    if (pd <= 0 && *reversed) {
	/* give up the "reversed" notion if we didn't get
	   a workable time-series interpretation */
	*reversed = 0;
    }

    return pd;
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
	    if (c == ';') {
		csv_set_got_semi(cdata);
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

#define nonspace_delim(d) (d != ',' && d != ';')

static int count_csv_fields (const char *s, char delim)
{
    int cbak, nf = 0;

    if (*s == delim && *s == ' ') {
	s++;
    }

    while (*s) {
	if (*s == delim) {
	    nf++;
	}
	cbak = *s;
	s++;
	/* Problem: (when) should a trailing delimiter be read as an
	   implicit NA?  For now we'll so treat it if the delimiter
	   is not white space.
	*/
	if (*s == '\0' && cbak == delim && nonspace_delim(delim)) {
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

static void compress_csv_line (csvdata *c)
{
    int n = strlen(c->line);
    char *p = c->line + n - 1;

    if (*p == '\n') {
	*p = '\0';
	p--;
    }

    if (*p == '\r') *p = '\0';

    if (c->delim == ',') {
	remove_quoted_commas(c->line);
    }

    if (c->delim != ' ') {
	delchar(' ', c->line);
    } else {
	compress_spaces(c->line);
    }

    delchar('"', c->line);

    if (csv_has_trailing_comma(c)) {
	/* chop trailing comma */
	n = strlen(c->line);
	if (n > 0) {
	    c->line[n-1] = '\0';
	}
    }
}

int import_obs_label (const char *s)
{
    char tmp[OBSLEN];

    if (s == NULL) {
	return 1;
    }

    if (*s == '"' || *s == '\'') s++;

    if (*s == '\0') {
	return 1;
    }

    if (strlen(s) > OBSLEN - 1) {
	return 0;
    }

    *tmp = '\0';
    strncat(tmp, s, OBSLEN - 1);
    lower(tmp);

    return (!strcmp(tmp, "obs") ||
	    !strcmp(tmp, "date") || 
	    !strcmp(tmp, "year") || 
	    !strcmp(tmp, "period"));    
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

	    series_set_flag(c->dinfo, v, VAR_DISCRETE);

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

		nnfrac = (nok == 0)? 1.0 : (double) nnon / (nnon + nok);
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

/* remedial function: we may be trying to read a "CSV" file
   that uses ',' as decimal separator, on a platform where
   '.' is the decimal separator
*/

static double try_comma_atof (const char *orig)
{
    char *test, s[32];
    double x;

    if (strlen(orig) > 31) {
	return NON_NUMERIC;
    }

    strcpy(s, orig);
    charsub(s, ',', '.');

    errno = 0;
    x = strtod(s, &test);

    if (*test == '\0' && errno != ERANGE) {
	return x;
    } else {
	return NON_NUMERIC;
    }
}

static int process_csv_obs (csvdata *c, int i, int t, PRN *prn)
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
	if (c->delim != ',' && get_local_decpoint() != ',' &&
	    strchr(c->str, ',') != NULL) {
	    c->Z[i][t] = try_comma_atof(c->str);
	} else {
	    c->Z[i][t] = NON_NUMERIC;
	}
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

static int csv_fields_check (FILE *fp, csvdata *c, PRN *prn)
{
    int gotdata = 0;
    int chkcols = 0;
    int err = 0;

    c->ncols = c->nrows = 0;

    while (csv_fgets(c->line, c->maxlen, fp) && !err) {

	/* skip comment lines */
	if (*c->line == '#') {
	    continue;
	}

	/* skip blank lines -- but finish if the blank comes after data */
	if (string_is_blank(c->line)) {
	    if (gotdata) {
		if (!csv_have_data(c)) {
		    c->descrip = get_csv_descrip(c->line, c->maxlen, fp);
		}
		break;
	    } else {
		continue;
	    }
	}
	
	c->nrows += 1;

	if (c->cols_list != NULL) {
	    tailstrip(c->line);
	    gotdata = 1;
	    chkcols = strlen(c->line);
	    if (chkcols < c->cols_list[c->cols_list[0]]) {
		gretl_errmsg_set(_("Invalid column specification"));
		err = E_DATA;
		break;
	    } else {
		continue;
	    }
	}

	compress_csv_line(c);

	if (!gotdata) {
	    /* scrutinize first "real" line */
	    check_first_field(c->line, c, prn);
	    gotdata = 1;
	} 

	chkcols = count_csv_fields(c->line, c->delim);
	if (c->ncols == 0) {
	    c->ncols = chkcols;
	    pprintf(prn, M_("   number of columns = %d\n"), c->ncols);	    
	} else if (chkcols != c->ncols) {
	    pprintf(prn, M_("   ...but row %d has %d fields: aborting\n"),
		    c->nrows, chkcols);
	    err = E_DATA;
	}
    }

    if (!err && c->cols_list != NULL) {
	c->ncols = c->cols_list[0];
    }

    return err;
}

static void strip_illegals (char *s)
{
    int i;

    if (*s == '\0') return;

    for (i=1; s[i]!='\0'; i++) {
	if (!isalnum(s[i])) {
	    s[i] = '_';
	}
    }
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

static int csv_varname_scan (csvdata *c, FILE *fp, PRN *prn, PRN *mprn)
{
    char *p;
    int obscol = csv_has_obs_column(c);
    int i, k, numcount;
    int err = 0;

    pputs(mprn, M_("scanning for variable names...\n"));

    while (csv_fgets(c->line, c->maxlen, fp)) {
	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	} else {
	    break;
	}
    }

    compress_csv_line(c);   

    p = c->line;
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
		    strip_illegals(c->dinfo->varname[nv]);
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

/* read numerical data when we've been given a fixed column-reading
   specification */

static int fixed_format_read (csvdata *c, FILE *fp, PRN *prn)
{
    char *p;
    int i, k, n, m, t = 0;
    int err = 0;

    c->real_n = c->dinfo->n;

    while (csv_fgets(c->line, c->maxlen, fp) && !err) {

	tailstrip(c->line);

	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	}

	m = strlen(c->line);

	for (i=1; i<=c->ncols; i++) {
	    k = c->cols_list[i];
	    n = c->width_list[i];
	    if (k + n - 1 > m) {
		/* attempting to read out of bounds */
		fprintf(stderr, "row %d, column %d: start=%d, width=%d, "
			"but line length = %d\n", t+1, i, k, n, m);
		err = E_DATA;
		break;
	    }		
	    p = c->line + k - 1;
	    *c->str = '\0';
	    strncat(c->str, p, n);
	    if (csv_missval(c->str, i, t+1, prn)) {
		c->Z[i][t] = NADBL;
	    } else if (check_atof(c->str)) {
		if (c->delim != ',' && get_local_decpoint() != ',' &&
		    strchr(c->str, ',') != NULL) {
		    c->Z[i][t] = try_comma_atof(c->str);
		} else {
		    err = E_DATA;
		}
	    } else {
		c->Z[i][t] = atof(c->str);
	    }
	}

	if (++t == c->dinfo->n) {
	    break;
	}
    }

    if (err == E_DATA) {
	gretl_errmsg_set(_("Invalid column specification"));
    }

    return err;
}

static int 
real_read_labels_and_data (csvdata *c, FILE *fp, PRN *prn)
{
    char *p;
    int i, k, nv, t = 0;
    int err = 0;

    c->real_n = c->dinfo->n;

    while (csv_fgets(c->line, c->maxlen, fp) && !err) {

	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	}

	if (csv_skipping(c) && strstr(c->line, c->skipstr)) {
	    c->real_n -= 1;
	    continue;
	}

	compress_csv_line(c);
	p = c->line;
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
		err = process_csv_obs(c, nv, t, prn);
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

static int csv_read_data (csvdata *c, FILE *fp, PRN *prn, PRN *mprn)
{
    int reversed = csv_data_reversed(c);
    int err;

    pputs(mprn, M_("scanning for row labels and data...\n"));
    err = real_read_labels_and_data(c, fp, prn);

    if (!err && csv_skip_column(c)) {
	c->markerpd = test_markers_for_dates(&c->Z, c->dinfo, &reversed,
					     c->skipstr, prn);
	if (reversed) {
	    csv_set_data_reversed(c);
	}
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
 * @fname: name of CSV file.
 * @pZ: pointer to data set.
 * @pdinfo: pointer to data information struct.
 * @opt: use %OPT_C to force interpretation of data colums containing
 * strings as coded values and not errors; for use of %OPT_T see
 * the help for "append".
 * @prn: gretl printing struct (or %NULL).
 * 
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_csv (const char *fname, double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn)
{
    csvdata *c = NULL;
    const char *cols = NULL;
    int popit = 0;
    FILE *fp = NULL;
    PRN *mprn = NULL;
    int newdata = (*pZ == NULL);
    char save_delim = pdinfo->delim;
    int fixed_format = 0;
    long datapos;
    int i, err = 0;

    if (opt & OPT_Q) {
	/* quiet */
	prn = NULL;
    }

    if (opt & OPT_F) {
	/* fixed format: should have --cols=XXX specification */
	cols = get_optval_string(OPEN, OPT_F);
    }

    if (prn != NULL) {
	check_for_console(prn);
    }

    if (gretl_messages_on()) {
	mprn = prn;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, M_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto csv_bailout;
    }

    c = csvdata_new(*pZ, pdinfo, opt);
    if (c == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }

    if (cols != NULL && *cols != '\0') {
	fixed_format = 1;
	pprintf(mprn, M_("using fixed column format\n"));
	err = csvdata_add_cols_list(c, cols);
	if (err) {
	    goto csv_bailout;
	}
    }

    if (mprn != NULL) {
	csv_parsing_header(fname, mprn);
    }

    /* get line length, also check for binary data, etc. */
    c->maxlen = csv_max_line_length(fp, c, prn);    
    if (c->maxlen <= 0) {
	err = E_DATA;
	goto csv_bailout;
    } 

    if (!fixed_format && !csv_got_delim(c)) {
	/* set default delimiter */
	if (csv_got_tab(c)) {
	    c->delim = c->dinfo->delim = '\t';
	} else if (csv_got_semi(c)) {
	    c->delim = c->dinfo->delim = ';';
	} else {
	    c->delim = c->dinfo->delim = ' ';
	}
    }

 alt_delim:

    if (!fixed_format) {
	pprintf(mprn, M_("using delimiter '%c'\n"), c->delim);
    }

    pprintf(mprn, M_("   longest line: %d characters\n"), c->maxlen - 1);

    if (csv_has_trailing_comma(c) && c->delim != ',') {
	csv_unset_trailing_comma(c);
    }

    /* buffer to hold lines */
    c->line = malloc(c->maxlen);
    if (c->line == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }  
    
    rewind(fp);

    /* read lines, check for consistency in number of fields */
    err = csv_fields_check(fp, c, mprn);
    if (err && !fixed_format) {
	if (c->delim != ';' && csv_got_semi(c)) {
	    c->delim = c->dinfo->delim = ';';
	    err = 0;
	    goto alt_delim;
	}
	pputs(prn, M_(csv_msg));
	goto csv_bailout;
    }

    if (fixed_format) {
	c->dinfo->n = c->nrows;
	c->dinfo->v = c->ncols + 1;
    } else {
	c->dinfo->n = c->nrows - 1; /* allow for var headings */
	c->dinfo->v = (csv_skip_column(c))? c->ncols : c->ncols + 1;
    }

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

    if (fixed_format) {
	err = fixed_format_read(c, fp, prn);
	if (err) {
	    goto csv_bailout;
	} else {
	    csv_set_autoname(c);
	    goto csv_continue;
	}
    }

    err = csv_varname_scan(c, fp, prn, mprn);
    if (err) {
	goto csv_bailout;
    }

    if (pdinfo != NULL && pdinfo->decpoint != ',') {
	gretl_push_c_numeric_locale();
	popit = 1;
    }

    datapos = ftell(fp);

    err = csv_read_data(c, fp, prn, mprn);

    if (!err && csv_skipping(c)) {
	/* try again */
	fseek(fp, datapos, SEEK_SET);
	err = csv_read_data(c, fp, prn, NULL);
    }

    if (!err) {
	err = non_numeric_check(c, prn);
	if (!err && csv_has_non_numeric(c)) {
	    /* try once more */
	    fseek(fp, datapos, SEEK_SET);
	    err = csv_read_data(c, fp, prn, NULL);
	}
    }	

    if (popit) {
	gretl_pop_c_numeric_locale();
    }

    if (err) {
	goto csv_bailout;
    }

    if (csv_data_reversed(c)) {
	reverse_data(c->Z, c->dinfo, mprn);
    }

 csv_continue:

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

    err = merge_or_replace_data(pZ, pdinfo, &c->Z, &c->dinfo, opt, prn);

    if (!err && newdata && c->descrip != NULL) {
	pdinfo->descrip = c->descrip;
	c->descrip = NULL;
    }

    if (!err) {
	dataset_add_import_info(pdinfo, fname, GRETL_CSV);
    }

 csv_bailout:

    if (fp != NULL) {
	fclose(fp);
    }

    csvdata_free(c);

    if (err == E_ALLOC) {
	pputs(prn, M_("Out of memory\n"));
    }    

    console_off();

    pdinfo->delim = save_delim;

    return err;
}
