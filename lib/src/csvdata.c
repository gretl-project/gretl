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
#include "usermat.h"
#include "genparse.h"
#include "csvdata.h"

#include <ctype.h>
#include <errno.h>
#include <glib.h>

#define CDEBUG 0

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
    CSV_REVERSED = 1 << 8,
    CSV_DOTSUB   = 1 << 9
};

enum {
    JOIN_KEY = 1,
    JOIN_VAL,
    JOIN_LHS,
    JOIN_RHS,
    JOIN_KEY2,
    JOIN_AUX
};

#define JOIN_MAXCOL 6

typedef struct csvjoin_ csvjoin;
typedef struct csvdata_ csvdata;

struct csvjoin_ {
    const char *colnames[JOIN_MAXCOL];
    char colnums[JOIN_MAXCOL];
    csvdata *c;
};    

struct csvdata_ {
    int flags;
    char delim;
    char decpoint;
    int markerpd;
    int maxlen;
    int real_n;
    char *line;
    DATASET *dset;
    int ncols, nrows;
    char str[CSVSTRLEN];
    char skipstr[8];
    int *codelist;
    char *descrip;
    gretl_string_table *st;
    int *cols_list;
    int *width_list;
    const gretl_matrix *rowmask;
    int masklen;
    csvjoin *jspec; /* info used for "join" command */
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
#define csv_data_reversed(c)      (c->flags & CSV_REVERSED)
#define csv_do_dotsub(c)          (c->flags & CSV_DOTSUB)

#define csv_set_trailing_comma(c)   (c->flags |= CSV_TRAIL)
#define csv_unset_trailing_comma(c) (c->flags &= ~CSV_TRAIL)
#define csv_set_obs_column(c)       (c->flags |= CSV_OBS1)
#define csv_set_blank_column(c)     (c->flags |= CSV_BLANK1)
#define csv_set_got_tab(c)          (c->flags |= CSV_GOTTAB)
#define csv_set_got_semi(c)         (c->flags |= CSV_GOTSEMI)
#define csv_set_got_delim(c)        (c->flags |= CSV_GOTDELIM)
#define csv_set_autoname(c)         (c->flags |= CSV_AUTONAME)
#define csv_set_data_reversed(c)    (c->flags |= CSV_REVERSED)
#define csv_set_dotsub(c)           (c->flags |= CSV_DOTSUB)

#define csv_skip_bad(c)        (*c->skipstr != '\0')
#define csv_has_non_numeric(c) (c->st != NULL)

#define fixed_format(c) (c->cols_list != NULL && c->width_list != NULL)
#define cols_subset(c) (c->cols_list != NULL && c->width_list == NULL)
#define rows_subset(c) (c->rowmask != NULL)

#define joining(c) (c->jspec != NULL)

static int
time_series_label_check (DATASET *dset, int reversed, char *skipstr, PRN *prn);

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

    destroy_dataset(c->dset);

    free(c);
}

static csvdata *csvdata_new (DATASET *dset)
{
    csvdata *c = malloc(sizeof *c);

    if (c == NULL) {
	return NULL;
    }

    c->flags = 0;
    c->delim = '\t';
    c->markerpd = -1;
    c->maxlen = 0;
    c->real_n = 0;
    c->line = NULL;
    c->dset = NULL;
    c->ncols = 0;
    c->nrows = 0;
    *c->str = '\0';
    *c->skipstr = '\0';
    c->codelist = NULL;
    c->descrip = NULL;
    c->st = NULL;
    c->cols_list = NULL;
    c->width_list = NULL;
    c->rowmask = NULL;
    c->masklen = 0;

    c->jspec = NULL;

    c->dset = datainfo_new();
    if (c->dset == NULL) {
	free(c);
	c = NULL;
    } else {
	c->delim = get_data_export_delimiter();
	c->decpoint = get_data_export_decpoint();
	if (dset->Z != NULL) {
	    c->flags |= CSV_HAVEDATA;
	}
    }

    return c;
}

static int *cols_list_from_matrix (const char *s, int *err)
{
    gretl_matrix *m = get_matrix_by_name(s);
    int i, n = gretl_vector_get_length(m);
    int *list = NULL;

    if (n == 0) {
	*err = E_DATA;
    } else {
	list = gretl_list_new(n);
	if (list == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		list[i+1] = gretl_vector_get(m, i);
	    }
	}
    }

    return list;
}

/* The interpretation of the "cols" specification depends on
   @opt: if this includes OPT_L then it should provide a 1-based 
   list of columns to be read; but if @opt includes OPT_F it 
   should provide a fixed-format spec, consisting of pairs 
   (start column, width).
*/

static int csvdata_add_cols_list (csvdata *c, const char *s,
				  gretlopt opt)
{
    int delimited = (opt & OPT_L);
    int *list, *clist = NULL, *wlist = NULL;
    int i, n, m = 0;
    int err = 0;

    if (gretl_is_matrix(s)) {
	list = cols_list_from_matrix(s, &err);
    } else {
	list = gretl_list_from_string(s, &err);
    }

    if (!err) {
	n = list[0];
	if (n == 0) {
	    err = E_DATA;
	} else if (delimited) {
	    m = n;
	    clist = list;
	} else {
	    /* fixed format: we need two lists */
	    if (n % 2 != 0) {
		err = E_DATA;
	    } else {
		m = n / 2;
		clist = gretl_list_new(m);
		wlist = gretl_list_new(m);
		if (clist == NULL || wlist == NULL) {
		    err = E_ALLOC;
		} else {
		    int j = 1;
	
		    for (i=1; i<=n; i+=2, j++) {
			clist[j] = list[i];
			wlist[j] = list[i+1];
		    }
		}
	    }
	}
    }

    /* clist = column (start) list: must be a set of increasing
       positive integers; and wlist = respective column widths,
       must all be positive, if present
    */

    for (i=1; i<=m && !err; i++) {
	if (clist[i] <= 0 || (i > 1 && clist[i] <= clist[i-1])) {
	    err = E_DATA;
	} else if (wlist != NULL && wlist[i] <= 0) {
	    err = E_DATA;
	} else if (wlist != NULL && wlist[i] >= CSVSTRLEN) {
	    fprintf(stderr, "Warning: field %d too wide (%d), truncating\n", 
		    i, wlist[i]);
	    wlist[i] = CSVSTRLEN - 1;
	}
    }

    if (list != clist) {
	free(list);
    }

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

static int csvdata_add_row_mask (csvdata *c, const char *s)
{
    int err = 0;

    c->rowmask = get_matrix_by_name(s);
    if (c->rowmask == NULL) {
	gretl_errmsg_sprintf(_("'%s': no such matrix"), s);
	err = E_DATA;
    } else {
	c->masklen = gretl_vector_get_length(c->rowmask);
	if (c->masklen == 0) {
	    err = E_NONCONF;
	}
    }

    return err;
}

static int n_from_row_mask (csvdata *c)
{
    int i, n = 0;

    for (i=0; i<c->masklen && i<=c->nrows; i++) {
	if (gretl_vector_get(c->rowmask, i) != 0) {
	    n++;
	}
    }

    return n;
}

static int add_obs_marker (DATASET *dset, int n)
{
    char **S = realloc(dset->S, n * sizeof *S);
    int err = 0;

    if (S == NULL) {
	err = E_ALLOC;
    } else {
	dset->S = S;
	dset->S[n-1] = malloc(OBSLEN);
	if (dset->S[n-1] == NULL) {
	    err = E_ALLOC;
	} else {
	    strcpy(dset->S[n-1], "NA");
	}
    }

    return err;
}

static int add_single_obs (DATASET *dset)
{
    double *x;
    int i, err = 0;

    for (i=0; i<dset->v && !err; i++) {
	x = realloc(dset->Z[i], (dset->n + 1) * sizeof *x);
	if (x != NULL) {
	    dset->Z[i] = x;
	} else {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	dset->n += 1;
	dset->Z[0][dset->n - 1] = 1.0;
	for (i=1; i<dset->v; i++) {
	    dset->Z[i][dset->n - 1] = NADBL;
	}
	if (dset->S != NULL) {
	    err = add_obs_marker(dset, dset->n);
	}
    }

    return err;
}

static int pad_weekly_data (DATASET *dset, int add)
{
    int oldn = dset->n;
    int ttarg, offset = 0, skip = 0;
    int i, s, t, tc, err;

    err = dataset_add_observations(add, dset, OPT_A); 

    if (!err) {
	for (t=0; t<oldn; t++) {
	    tc = calendar_obs_number(dset->S[t], dset) - offset;
	    if (tc != t) {
		skip = tc - t;
		fprintf(stderr, "Gap of size %d at original t = %d\n", skip, t);
		offset += skip;
		ttarg = oldn - 1 + offset;
		for (s=0; s<oldn-t+skip; s++) {
		    for (i=1; i<dset->v; i++) {
			if (s < oldn - t) {
			    if (s == 0 || s == oldn-t-1) {
				fprintf(stderr, "shifting obs %d to obs %d\n",
					ttarg-skip, ttarg);
			    }
			    dset->Z[i][ttarg] = dset->Z[i][ttarg - skip];
			} else {
			    fprintf(stderr, "inserting NA at obs %d\n", ttarg);
			    dset->Z[i][ttarg] = NADBL;
			}
		    }
		    ttarg--;
		}
	    }
	}
    }

    return err;
}

/* FIXME the following needs to be made more flexible? */

static int csv_weekly_data (DATASET *dset)
{
    char *lbl2 = dset->S[dset->n - 1];
    int ret = 1;
    int misscount = 0;
    int t, tc;

    for (t=0; t<dset->n; t++) {
	tc = calendar_obs_number(dset->S[t], dset) - misscount;
	if (tc != t) {
	    misscount += tc - t;
	}
    }

    if (misscount > 0) {
	double missfrac = (double) misscount / dset->n;

	fprintf(stderr, "nobs = %d, misscount = %d (%.2f%%)\n", 
		dset->n, misscount, 100.0 * missfrac);
	if (missfrac > 0.05) {
	    ret = 0;
	} else {
	    int Tc = calendar_obs_number(lbl2, dset) + 1;
	    int altmiss = Tc - dset->n;

	    fprintf(stderr, "check: Tc = %d, missing = %d\n", Tc, altmiss);
	    if (altmiss != misscount) {
		ret = 0;
	    } else if (dset->Z != NULL) {
		int err;

		fprintf(stderr, "OK, consistent\n");
		err = pad_weekly_data(dset, misscount);
		if (err) ret = 0;
	    } 
	} 
    }
    
    return ret;
}

#define DAY_DEBUG 1

static int check_daily_dates (DATASET *dset, int *pd, 
			      int *reversed, PRN *prn)
{
    int T = dset->n;
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[T - 1];
    int fulln = 0, n, t, nbak;
    int alt_pd = 0;
    int oldpd = dset->pd;
    double oldsd0 = dset->sd0;
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

    dset->pd = guess_daily_pd(dset);
    dset->structure = TIME_SERIES;

#if DAY_DEBUG    
    fprintf(stderr, "guessed daily pd = %d\n", dset->pd);
#endif

    if (!err) {
	ed2 = get_epoch_day(lbl2);
	if (ed2 < 0) {
	    err = 1;
	} else if (ed2 < ed1) {
#if DAY_DEBUG    
	    fprintf(stderr, "check_daily_dates: data are reversed?\n");
#endif
	    dset->sd0 = ed2;
	    *reversed = 1;
	} else {
	    dset->sd0 = ed1;
	}
    }

 recompute:

    alt_pd = 0;
    nbak = 0;

    if (!err) {
	int n1 = calendar_obs_number((*reversed)? lbl2 : lbl1, dset);
	int n2 = calendar_obs_number((*reversed)? lbl1 : lbl2, dset);

	fulln = n2 - n1 + 1;

	if (T > fulln) {
	    err = 1;
	} else {
	    nmiss = fulln - T;
	    pprintf(prn, A_("Observations: %d; days in sample: %d\n"), 
		    T, fulln);
	    if (nmiss > 300 * T) {
		pprintf(prn, A_("Probably annual data\n"));
		*pd = 1;
	    } else if (nmiss > 50 * T) {
		pprintf(prn, A_("Probably quarterly data\n"));
		*pd = 4;
	    } else if (nmiss > 20 * T) {
		pprintf(prn, A_("Probably monthly data\n"));
		*pd = 12;
	    } else if (nmiss > 3 * T) {
		pprintf(prn, A_("Probably weekly data\n"));
		*pd = dset->pd = 52;
	    } else {
		pprintf(prn, A_("Missing daily observations: %d\n"), nmiss);
	    }
	}
    }

    nbak = 0;

    for (t=0; t<dset->n && !err; t++) {
	int wd, s = (*reversed)? (dset->n - 1 - t) : t;

	wd = get_day_of_week(dset->S[s]);

	if (dset->pd == 5 && (wd == 6 || wd == 0)) {
	    /* Got Sat or Sun, can't be 5-day daily? */
	    alt_pd = (wd == 6)? 6 : 7;
	    pprintf(prn, "Found a Saturday (%s): re-trying with pd = %d\n", 
		    dset->S[s], alt_pd);
	    break;
	} else if (dset->pd == 6 && wd == 0) {
	    /* Got Sun, can't be 6-day daily? */
	    alt_pd = 7;
	    pprintf(prn, "Found a Sunday (%s): re-trying with pd = %d\n", 
		    dset->S[s], alt_pd);
	    break;
	}
	    
	n = calendar_obs_number(dset->S[s], dset);
	if (n < t) {
	    pprintf(prn, "Daily dates error at t = %d:\n"
		    "  calendar_obs_number() for '%s' = %d but t = %d\n", 
		    t, dset->S[t], n, t);
	    err = 1;
	} else if (n > fulln - 1) {
	    pprintf(prn, "Error: date '%s' out of bounds\n", dset->S[s]);
	    err = 1;
	} else if (nbak > 0 && n == nbak) {
	    pprintf(prn, "Error: date '%s' is repeated\n", dset->S[s]);
	    err = 1;
	}
	nbak = n;
    }

    if (alt_pd > 0) {
	dset->pd = alt_pd;
	goto recompute;
    }

    if (err) {
	dset->pd = oldpd;
	dset->sd0 = oldsd0;
	dset->structure = CROSS_SECTION;
    } else {
	strcpy(dset->stobs, (*reversed)? lbl2 : lbl1);
	strcpy(dset->endobs, (*reversed)? lbl1 : lbl2);
	dset->t2 = dset->n - 1;
	if (nmiss > 0 && *pd == 0) {
	    dset->markers = DAILY_DATE_STRINGS;
	}
    }

#if DAY_DEBUG
    fprintf(stderr, "check_daily_dates: pd = %d, reversed = %d, err = %d\n", 
	    dset->pd, *reversed, err);
#endif

    return (err)? -1 : dset->pd;
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

static int complete_qm_labels (DATASET *dset, int reversed,
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

    s = (reversed)? (dset->n - 1) : 0;
    if (sscanf(dset->S[s], fmt, &yr, &per) != 2) {
	return 0;
    }

    for (t=1; t<dset->n; t++) {
	s = (reversed)? (dset->n - 1 - t) : t;
	Ey = (per == pd)? yr + 1 : yr;
	Ep = (per == pd)? pmin : per + pmin;
	if (sscanf(dset->S[s], fmt, &yr, &per) != 2) {
	    ret = 0;
	} else if (Ep == 1 && pd == pd0 && per == pd + 1 
		   && skipstr != NULL) {
	    *skip = *bad = '\0';
	    strncat(skip, dset->S[s] + 4, 7); 
	    strncat(bad, dset->S[s], 15); 
	    pd = pd0 + 1;
	    goto restart;
	} else if (per == Ep + 2 && pmin == 1 && fakequarter(per)) {
	    *bad = '\0';
	    strncat(bad, dset->S[s], 15); 
	    pmin = 3;
	    goto restart;
	} else if (yr != Ey || per != Ep) {
	    ret = 0;
	}
	if (!ret) {
	    pprintf(prn, "   %s: not a consistent date\n", 
		    dset->S[s]);
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

static int complete_year_labels (const DATASET *dset, int reversed)
{
    int s, t, yr, yrbak;
    int ret = 1;

    s = (reversed)? (dset->n - 1) : 0;
    yrbak = atoi(dset->S[s]);

    for (t=1; t<dset->n; t++) {
	s = (reversed)? (dset->n - 1 - t) : t;
	yr = atoi(dset->S[s]);
	if (yr != yrbak + 1) {
	    ret = 0;
	    break;
	}
	yrbak = yr;
    }

    return ret;
}

static int compress_daily (DATASET *dset, int pd)
{
    int t, yr, mon, day;

    for (t=0; t<dset->n; t++) {
	sscanf(dset->S[t], "%d/%d/%d", &yr, &mon, &day);
	if (pd == 1) {
	    sprintf(dset->S[t], "%d", yr);
	} else if (pd == 12) {
	    sprintf(dset->S[t], "%d:%02d", yr, mon);
	} else if (pd == 4) {
	    sprintf(dset->S[t], "%d:%d", yr, mon / 3 + (mon % 3 != 0));
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

static void retransform_daily_dates (DATASET *dset)
{
    int t, y, m, d;

    /* we apparently guessed wrongly at MMDDYYYY, so
       put the dates back as they were for another try,
       at DDMMYYYY.
    */

    for (t=0; t<dset->n; t++) {
	sscanf(dset->S[t], "%d/%d/%d", &y, &d, &m);
	sprintf(dset->S[t], "%d/%d/%d", d, m, y);
    }
}

static int transform_daily_dates (DATASET *dset, int dorder)
{
    char *label;
    int t, yr, mon, day;
    char s1, s2;
    int n, err = 0;

    for (t=0; t<dset->n && !err; t++) {
	label = dset->S[t];
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

void reverse_data (DATASET *dset, PRN *prn)
{
    char tmp[OBSLEN];
    double x;
    int T = dset->n / 2;
    int i, t, s;

    pprintf(prn, A_("reversing the data!\n"));

    for (t=0; t<T; t++) {
	s = dset->n - 1 - t;
	for (i=1; i<dset->v; i++) {
	    x = dset->Z[i][t];
	    dset->Z[i][t] = dset->Z[i][s];
	    dset->Z[i][s] = x;
	}
	if (dset->S != NULL) {
	    strcpy(tmp, dset->S[t]);
	    strcpy(dset->S[t], dset->S[s]);
	    strcpy(dset->S[s], tmp);
	}
    }
}

static int csv_daily_date_check (DATASET *dset, int *reversed,
				 char *skipstr, PRN *prn)
{
    int d1[3], d2[3];
    char s1, s2;
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[dset->n - 1];
    int dorder = 0;

    if (sscanf(lbl1, "%d%c%d%c%d", &d1[0], &s1, &d1[1], &s2, &d1[2]) == 5 &&
	sscanf(lbl2, "%d%c%d%c%d", &d2[0], &s1, &d2[1], &s2, &d2[2]) == 5 &&
	s1 == s2 && ispunct(s1)) {
	int mon1, day1;
	int mon2, day2;
	int pd, ret = 0;

	dorder = get_date_order(d1[0], d2[0]);

    tryagain:

	if (dorder == YYYYMMDD) {
	    pputs(prn, A_("Trying date order YYYYMMDD\n"));
	    mon1 = d1[1];
	    day1 = d1[2];
	    mon2 = d2[1];
	    day2 = d2[2];
	} else if (dorder == DDMMYYYY) {
	    pputs(prn, A_("Trying date order DDMMYYYY\n"));
	    day1 = d1[0];
	    mon1 = d1[1];
	    day2 = d2[0];
	    mon2 = d2[1];
	} else {
	    pputs(prn, A_("Trying date order MMDDYYYY\n"));
	    mon1 = d1[0];
	    day1 = d1[1];
	    mon2 = d2[0];
	    day2 = d2[1];
	}		
	    
	if (mon1 > 0 && mon1 < 13 &&
	    mon2 > 0 && mon2 < 13 && 
	    day1 > 0 && day1 < 32 &&
	    day2 > 0 && day2 < 32) {
	    /* looks promising for calendar dates */
	    if (dorder != YYYYMMDD || s1 != '/' || s2 != '/') {
		if (transform_daily_dates(dset, dorder)) {
		    return -1;
		}
	    }
	    pprintf(prn, A_("Could be %s - %s\n"), lbl1, lbl2);
	    ret = check_daily_dates(dset, &pd, reversed, prn);
	    if (ret >= 0 && pd > 0) {
		if (pd == 52) {
		    if (csv_weekly_data(dset)) {
			ret = 52;
		    } else if (dorder == MMDDYYYY) {
			/* maybe we guessed wrong */
			retransform_daily_dates(dset);
			dorder = DDMMYYYY;
			goto tryagain;
		    } else {
			ret = -1;
		    }
		} else {
		    compress_daily(dset, pd);
		    ret = time_series_label_check(dset, 
						  *reversed,
						  skipstr, 
						  prn);
		}
	    } 
	    return ret;
	}
    } else {
	pprintf(prn, A_("'%s' and '%s': couldn't get dates\n"), lbl1, lbl2);
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
	pprintf(prn, A_("   %s: probably a year... "), year);
    } else {
	pprintf(prn, A_("   %s: probably not a year\n"), year);
    }

    if (len == 5) {
	pputs(prn, A_("   but I can't make sense of the extra bit\n"));
    } else if (len == 4) {
	pputs(prn, A_("and just a year\n"));
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
		    pprintf(prn, A_("quarter %s?\n"), s);
		    pd = 4;
		} else {
		    pprintf(prn, "quarter %d: not possible\n", p);
		}
	    } else if (len == 7) {
		p = atoi(s);
		if (p > 0 && p < 13) {
		    pprintf(prn, A_("month %s?\n"), s);
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

static int time_series_label_check (DATASET *dset, int reversed,
				    char *skipstr, PRN *prn)
{
    char year[5], sub[3];
    char format[8] = {0};
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[dset->n - 1];
    int pd = -1;

    pd = pd_from_date_label((reversed)? lbl2 : lbl1, year, sub, 
			    format, prn);

    if (pd == 1) {
	if (complete_year_labels(dset, reversed)) {
	    dset->pd = pd;
	    strcpy(dset->stobs, year);
	    dset->sd0 = atof(dset->stobs);
	    strcpy(dset->endobs, lbl2);
	} else {
	    pputs(prn, A_("   but the dates are not complete and consistent\n"));
	    pd = -1;
	}
    } else if (pd == 4 || pd == 12) {
	int savepd = pd;

	if (complete_qm_labels(dset, reversed, skipstr, &pd, format, prn)) {
	    dset->pd = pd;
	    if (savepd == 12 && pd == 4) {
		int s = atoi(sub) / 3;

		sprintf(dset->stobs, "%s:%d", year, s);
	    } else {
		sprintf(dset->stobs, "%s:%s", year, sub);
	    }
	    dset->sd0 = obs_str_to_double(dset->stobs);
	    ntodate(dset->endobs, dset->n - 1, dset);
	} else {
	    pputs(prn, A_("   but the dates are not complete and consistent\n"));
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
	pputs(prn, A_("   dates are reversed?\n"));
    }
    
    return ret;
}

/* e.g. "M1 1957", "M12 2009", with spaces removed */

static int fix_IFS_data_labels (DATASET *dset)
{
    char *s1 = dset->S[0];
    char *s2 = dset->S[dset->n - 1];
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

	    for (i=0; i<dset->n; i++) {
		s = dset->S[i];
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
		for (i=0; i<dset->n; i++) {
		    s = dset->S[i];
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
			dset->S[i] = gretl_strdup(tmp);
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

int test_markers_for_dates (DATASET *dset, int *reversed, 
			    char *skipstr, PRN *prn)
{
    char endobs[OBSLEN];
    int n = dset->n;
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[n - 1];
    int len1 = strlen(lbl1);
    int pd = -1;

    if (skipstr != NULL && *skipstr != '\0') {
	return time_series_label_check(dset, *reversed, skipstr, prn);
    }

    pprintf(prn, A_("   first row label \"%s\", last label \"%s\"\n"), 
	    lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", n);
    if (!strcmp(lbl1, "1") && !strcmp(lbl2, endobs)) {
	return 0;
    }

    if (fix_IFS_data_labels(dset)) {
	lbl1 = dset->S[0];
	lbl2 = dset->S[n - 1];
	len1 = strlen(lbl1);
    }

    /* labels are of different lengths? */
    if (len1 != strlen(lbl2)) {
	return -1;
    }

    pputs(prn, A_("trying to parse row labels as dates...\n"));

    if (len1 == 8 || len1 == 10) {
	/* daily data? */
	pd = csv_daily_date_check(dset, reversed, skipstr, prn);
    } else if (len1 >= 4) {
	/* annual, quarterly, monthly? */
	if (isdigit((unsigned char) lbl1[0]) &&
	    isdigit((unsigned char) lbl1[1]) &&
	    isdigit((unsigned char) lbl1[2]) && 
	    isdigit((unsigned char) lbl1[3])) {
	    *reversed = dates_maybe_reversed(lbl1, lbl2, prn);
	    pd = time_series_label_check(dset, *reversed, skipstr, prn);
	} else {
	    pputs(prn, A_("   definitely not a four-digit year\n"));
	}
    }

    if (pd <= 0 && *reversed) {
	/* give up the "reversed" notion if we didn't get
	   a workable time-series interpretation */
	*reversed = 0;
    }

    return pd;
}

static int utf8_ok (FILE *fp, int pos)
{
    long mark = ftell(fp);
    int len = pos + 9;
    char *test = malloc(len + 1);
    int i, ret = 0;

    fseek(fp, mark - pos - 1, SEEK_SET);

    for (i=0; i<len; i++) {
	test[i] = fgetc(fp);
    }
    test[i] = '\0';

    if (g_utf8_validate(test, -1, NULL)) {
	ret = 1;
    } 

    free(test);

    fseek(fp, mark, SEEK_SET);

    return ret;
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
    int lines = 0;

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
	    lines++;
	    continue;
	}
	cbak = c;
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ) && !utf8_ok(fp, cc)) {
	    pprintf(prn, A_("Binary data (%d) encountered (line %d:%d): "
			    "this is not a valid text file\n"), 
		    c, lines, cc);
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
	pprintf(prn, A_("Data file is empty\n"));
    } else if (csv_has_trailing_comma(cdata)) {
	pprintf(prn, A_("Data file has trailing commas\n"));
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

static void purge_quoted_commas (char *s)
{
    int inquote = 0;

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	} else if (inquote && *s == ',') {
	    *s = ' ';
	}
	s++;
    }
}

static void purge_unquoted_spaces (char *s)
{
    int inquote = 0;

    while (*s) {
	if (*s == '"') {
	    inquote = !inquote;
	} else if (!inquote && *s == ' ') {
	    shift_string_left(s, 1);
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
	purge_quoted_commas(c->line);
    }

    if (c->delim != ' ') {
	/* 2012-07-14: was gretl_delchar(' ', c->line); */
	purge_unquoted_spaces(c->line);
    } else {
	compress_spaces(c->line);
    }

    gretl_delchar('"', c->line);

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
    gretl_lower(tmp);

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

	pprintf(prn, A_("   first field: '%s'\n"), field1);

	if (import_obs_label(field1)) {
	    pputs(prn, A_("   seems to be observation label\n"));
	    csv_set_obs_column(c);
	}
    }
}

static int csv_missval (const char *str, int i, int t, 
			int *miss_shown, PRN *prn)
{
    int miss = 0;

    if (*str == '\0') {
	if (t < 80 || *miss_shown < i) {
	    pprintf(prn, A_("   the cell for variable %d, obs %d "
			    "is empty: treating as missing value\n"), 
		    i, t);
	    *miss_shown += 1;
	}
	miss = 1;
    }

    if (import_na_string(str)) {
	if (t < 80 || *miss_shown < i) {
	    pprintf(prn, A_("   warning: missing value for variable "
			    "%d, obs %d\n"), i, t);
	    *miss_shown += 1;
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

    for (i=1; i<c->dset->v; i++) {
	for (t=0; t<c->dset->n; t++) {
	    if (c->dset->Z[i][t] == NON_NUMERIC) {
		nn++;
		break;
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
	for (i=1; i<c->dset->v; i++) {
	    for (t=0; t<c->dset->n; t++) {
		if (c->dset->Z[i][t] == NON_NUMERIC) {
		    list[j++] = i;
		    break;
		}
	    }
	}

#if CDEBUG > 1
	printlist(list, "non-numeric vars list");
#endif

	for (i=1; i<=list[0]; i++) {
	    double nnfrac;
	    int nnon = 0;
	    int nok = 0;
	    int tn = 0;
	    int v = list[i];

	    series_set_flag(c->dset, v, VAR_DISCRETE);

	    for (t=0; t<c->dset->n; t++) {
		if (c->dset->Z[v][t] == NON_NUMERIC) {
		    if (!tn) tn = t + 1;
		    nnon++;
		} else if (!na(c->dset->Z[v][t])) {
		    nok++;
		}
	    }

	    nnfrac = (nok == 0)? 1.0 : (double) nnon / (nnon + nok);
	    pprintf(prn, "variable %d (%s): non-numeric values = %d "
		    "(%.2f percent)\n", v, c->dset->varname[v], 
		    nnon, 100 * nnfrac);
	    if (nnon < 2 || nnfrac < 0.01) {
		pprintf(prn, A_("ERROR: variable %d (%s), observation %d, "
				"non-numeric value\n"), 
			v, c->dset->varname[v], tn);
		err = E_DATA;
	    }
	}

	if (!err) {
	    pputs(prn, "allocating string table\n");
	    c->st = gretl_string_table_new(list);
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

static double csv_atof (csvdata *c, const char *s)
{
    double x = NON_NUMERIC;
    char *test;

    errno = 0;

    if (c->decpoint == '.' || !csv_do_dotsub(c) || strchr(s, ',') == NULL) {
	/* we should currently be set to the correct locale,
	   or there's no problematic decimal point in @s
	*/
	x = strtod(s, &test);
	if (*test == '\0' && errno == 0) {
	    return x;
	} else {
	    x = NON_NUMERIC;
	}
    } else if (csv_do_dotsub(c) && strlen(s) <= 31) {
	/* substitute dot for comma */
	char tmp[32];

	strcpy(tmp, s);
	gretl_charsub(tmp, ',', '.');

	x = strtod(tmp, &test);
	if (*test == '\0' && errno == 0) {
	    return x;
	} else {
	    x = NON_NUMERIC;
	}
    }

    if (c->decpoint == '.' && strchr(s, ',') != NULL && strlen(s) <= 31) {
	/* try remediation for decimal comma? */
	char tmp[32];

	strcpy(tmp, s);
	gretl_charsub(tmp, ',', '.');
	errno = 0;

	x = strtod(tmp, &test);
	if (*test != '\0' || errno != 0) {
	    x = NON_NUMERIC;
	} 
    }

    return x;
}

static int process_csv_obs (csvdata *c, int i, int t, int *miss_shown,
			    PRN *prn)
{
    int ix, err = 0;

    if (c->st != NULL) {
	if (in_gretl_list(c->codelist, i) && !na(c->dset->Z[i][t])) {
	    ix = gretl_string_table_index(c->st, c->str, i, 0, prn);
	    if (ix > 0) {
		c->dset->Z[i][t] = (double) ix;
	    } else {
		err = E_DATA;
	    }
	}
    } else if (csv_missval(c->str, i, t+1, miss_shown, prn)) {
	c->dset->Z[i][t] = NADBL;
    } else {
	c->dset->Z[i][t] = csv_atof(c, gretl_chopstr(c->str));
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
    size_t len;

    while (csv_fgets(line, n, fp)) {
	tailstrip(line);
	if (desc == NULL) {
	    len = strlen(line) + 3;
	    desc = malloc(len);
	    if (desc == NULL) {
		return NULL;
	    }
	    sprintf(desc, "%s\n", line);
	} else {
	    char *tmp;

	    len = strlen(desc) + strlen(line) + 3;
	    tmp = realloc(desc, len);
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

	if (fixed_format(c)) {
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
	    pprintf(prn, A_("   number of columns = %d\n"), c->ncols);	    
	} else if (chkcols != c->ncols) {
	    pprintf(prn, A_("   ...but row %d has %d fields: aborting\n"),
		    c->nrows, chkcols);
	    err = E_DATA;
	} else if (cols_subset(c)) {
	    int datacols = csv_skip_column(c) ? (c->ncols - 1) : c->ncols;

	    if (c->cols_list[c->cols_list[0]] > datacols) {
		gretl_errmsg_set(_("Invalid column specification"));
		err = E_DATA;
	    }
	}		
    }

    if (!err && fixed_format(c)) {
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

static int csv_reconfigure_for_markers (DATASET *dset)
{
    int err = dataset_allocate_obs_markers(dset);

    if (!err) {
	err = dataset_drop_last_variables(1, dset);
    }

    return err;
}

static int skip_data_column (csvdata *c, int k)
{
    int col = csv_skip_column(c) ? k : k + 1;

    if (!in_gretl_list(c->cols_list, col)) {
	return 1;
    } else {
	return 0;
    }
}

static int update_join_cols_list (csvdata *c, int k)
{
    int *test = gretl_list_append_term(&c->cols_list, k+1);
    int err = 0;

    if (test == NULL) {
	err = E_ALLOC;
    }

#if CDEBUG
    printlist(c->cols_list, "c->cols_list for join");
#endif

    return err;
}

static int handle_join_varname (csvdata *c, int k, int *pj)
{
    int i, j = *pj;

    for (i=0; i<JOIN_MAXCOL; i++) {
#if CDEBUG
	fprintf(stderr, "i = %d; %s <-> %s\n", i, c->str, c->jspec->colnames[i]);
#endif
	if (c->jspec->colnames[i] != NULL &&
	    !strcmp(c->str, c->jspec->colnames[i])) {
	    c->dset->varname[j][0] = '\0';
	    strncat(c->dset->varname[j], c->str, VNAMELEN - 1);
	    update_join_cols_list(c, k);
	    c->jspec->colnums[j-1] = i + 1;
	    *pj += 1;
	    break;
	}
    }

    return 0;
}

#define starts_number(c) (isdigit((unsigned char) c) || c == '-' || \
                          c == '+' || c == '.')

#define obs_labels_no_varnames(o,c,n)  (!o && c->v > 3 && n == c->v - 2)

static int csv_varname_scan (csvdata *c, FILE *fp, PRN *prn, PRN *mprn)
{
    char *p;
    int obscol = csv_has_obs_column(c);
    int i, j, k, numcount;
    int err = 0;

    pputs(mprn, A_("scanning for variable names...\n"));

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
    if (strlen(p) > 118) {
	pprintf(mprn, A_("   line: %.115s...\n"), p);
    } else {
	pprintf(mprn, A_("   line: %s\n"), p);
    }
    
    numcount = 0;
    j = 1; /* for the constant */

    for (k=0; k<c->ncols && !err; k++) {
	i = 0;
	while (*p && *p != c->delim) {
	    if (i < CSVSTRLEN - 1) {
		c->str[i++] = *p;
	    }
	    p++;
	}
	c->str[i] = '\0';
	if (*p == c->delim) p++;

	if (k == 0 && (csv_skip_column(c))) {
	    ; /* no-op */
	} else if (!joining(c) && cols_subset(c) && skip_data_column(c, k)) {
	    ; /* no-op */
	} else {
	    if (*c->str == '\0') {
		pprintf(prn, A_("   variable name %d is missing: aborting\n"), j);
		pputs(prn, A_(csv_msg));
		err = E_DATA;
	    } else if (joining(c)) {
		handle_join_varname(c, k, &j);
	    } else {
		c->dset->varname[j][0] = '\0';
		strncat(c->dset->varname[j], c->str, VNAMELEN - 1);
		if (starts_number(*c->str)) {
		    numcount++;
		} else {
		    iso_to_ascii(c->dset->varname[j]);
		    strip_illegals(c->dset->varname[j]);
		    if (check_varname(c->dset->varname[j])) {
			errmsg(1, prn);
			err = E_DATA;
		    }
		}
		j++;
	    }
	}
	if (j == c->dset->v) {
#if CDEBUG
	    fprintf(stderr, "breaking on j = %d (k = %d)\n", j, k);
#endif	    
	    break;
	}
    }

    if (!err && joining(c) && c->cols_list == NULL) {
	/* no relevant columns were found */
	err = E_UNKVAR;
    }

    if (err) {
	return err;
    }

    if (numcount == c->dset->v - 1 || 
	obs_labels_no_varnames(obscol, c->dset, numcount)) {
	pputs(prn, A_("it seems there are no variable names\n"));
	/* then we undercounted the observations by one? */
	if (!rows_subset(c)) {
	    err = add_single_obs(c->dset);
	}
	if (!err) {
	    csv_set_autoname(c);
	    rewind(fp);
	    if (obs_labels_no_varnames(obscol, c->dset, numcount)) {
		err = csv_reconfigure_for_markers(c->dset);
		if (!err) {
		    csv_set_obs_column(c);
		}
	    }
	}
    } else if (numcount > 0) {
	for (i=1; i<c->dset->v; i++) {
	    if (check_varname(c->dset->varname[i])) {
		errmsg(1, prn);
		break;
	    }
	}
	err = E_DATA;
    }

    return err;
}

static int row_not_wanted (csvdata *c, int t)
{
    if (c->rowmask != NULL) {
	if (t >= c->masklen) {
	    return 1;
	} else if (gretl_vector_get(c->rowmask, t) == 0) {
	    return 1;
	}
    }

    return 0;
}

/* read numerical data when we've been given a fixed column-reading
   specification */

static int fixed_format_read (csvdata *c, FILE *fp, PRN *prn)
{
    char *p;
    int miss_shown = 0;
    int t = 0, s = 0;
    int i, k, n, m;
    int err = 0;

    c->real_n = c->dset->n;

    while (csv_fgets(c->line, c->maxlen, fp) && !err) {

	tailstrip(c->line);

	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	}

	if (row_not_wanted(c, s)) {
	    s++;
	    continue;
	}

	m = strlen(c->line);

	for (i=1; i<=c->ncols && !err; i++) {
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
	    if (csv_missval(c->str, i, t+1, &miss_shown, prn)) {
		c->dset->Z[i][t] = NADBL;
	    } else {
		c->dset->Z[i][t] = csv_atof(c, c->str);
		if (c->dset->Z[i][t] == NON_NUMERIC) {
		    gretl_errmsg_sprintf(_("At row %d, column %d:\n"), t+1, k);
		    gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!"),
					 c->str);
		    err = E_DATA;
		}
	    }
	}

	s++;
	if (++t == c->dset->n) {
	    break;
	}
    }

    if (err == E_DATA) {
	gretl_errmsg_set(_("Invalid column specification"));
    }

    return err;
}

#define XML1_OK(u) ((u>=0x0020 && u<=0xD7FF) || \
		    (u>=0xE000 && u<=0xFFFD))

/* check that observation label contains only 
   valid UTF-8, and moreover that every character
   is valid in XML 1.0
*/

static int label_is_valid (gchar *s)
{
    if (!g_utf8_validate(s, -1, NULL)) {
	return 0;
    } else {
	int i, n = g_utf8_strlen(s, -1);
	gunichar u;

	for (i=0; i<n; i++) {
	    u = g_utf8_get_char(s);
	    if (!XML1_OK(u)) {
		return 0;
	    }
	    s = g_utf8_next_char(s);
	}
    }

    return 1;
}

static int 
real_read_labels_and_data (csvdata *c, FILE *fp, PRN *prn)
{
    char *p;
    int miss_shown = 0;
    int t = 0, s = 0;
    int i, j, k;
    int err = 0;
    int truncated = 0;

    c->real_n = c->dset->n;

    while (csv_fgets(c->line, c->maxlen, fp) && !err) {

	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	}

	if (*c->skipstr != '\0' && strstr(c->line, c->skipstr)) {
	    c->real_n -= 1;
	    continue;
	}

	if (row_not_wanted(c, s)) {
	    s++;
	    continue;
	}

	compress_csv_line(c);
	p = c->line;
	if (c->delim == ' ' && *p == ' ') p++;

	j = 1;
	for (k=0; k<c->ncols && !err; k++) {
	    i = 0;
	    while (*p && *p != c->delim) {
		if (i < CSVSTRLEN - 1) {
		    c->str[i++] = *p;
		} else {
		    truncated++;
		}
		p++;
	    }
	    if (*p == c->delim) {
		p++;
	    }
	    c->str[i] = '\0';
	    if (k == 0 && csv_skip_column(c) && c->dset->S != NULL) {
		char *S = c->str;

		c->dset->S[t][0] = 0;
		if (*S == '"' || *S == '\'') {
		    S++;
		}
		strncat(c->dset->S[t], S, OBSLEN - 1);
		if (!label_is_valid((gchar *) c->dset->S[t])) {
		    iso_to_ascii(c->dset->S[t]);
		} 
	    } else if (cols_subset(c) && skip_data_column(c, k)) {
		; /* no-op */
	    } else {
		err = process_csv_obs(c, j++, t, &miss_shown, prn);
	    }
	}

	s++;
	if (++t == c->dset->n) {
	    break;
	}
    }

    if (truncated) {
	pprintf(prn, A_("warning: %d labels were truncated.\n"), truncated);
    }

    if (!err && c->real_n < c->dset->n) {
	int drop = c->dset->n - c->real_n;

	err = dataset_drop_observations(drop, c->dset);
    }

    return err;
}

static int csv_read_data (csvdata *c, FILE *fp, PRN *prn, PRN *mprn)
{
    int reversed = csv_data_reversed(c);
    int err;

    pputs(mprn, A_("scanning for row labels and data...\n"));
    err = real_read_labels_and_data(c, fp, prn);

    if (!err && csv_skip_column(c) && !rows_subset(c)) {
	c->markerpd = test_markers_for_dates(c->dset, &reversed,
					     c->skipstr, prn);
	if (reversed) {
	    csv_set_data_reversed(c);
	}
    }

    return err;
}

static void print_csv_parsing_header (const char *fname, PRN *prn)
{
    if (!g_utf8_validate(fname, -1, NULL)) {
	gchar *trfname = g_locale_to_utf8(fname, -1, NULL, NULL, NULL);

	pprintf(prn, "%s %s...\n", A_("parsing"), trfname);
	g_free(trfname);
    } else {
	pprintf(prn, "%s %s...\n", A_("parsing"), fname);
    }
}

static void csv_set_dataset_dimensions (csvdata *c)
{
    if (rows_subset(c)) {
	c->dset->n = n_from_row_mask(c);
    }

    if (fixed_format(c)) {
	if (c->dset->n == 0) {
	    c->dset->n = c->nrows;
	}
	c->dset->v = c->ncols + 1;
    } else {
	if (c->dset->n == 0) {
	    c->dset->n = c->nrows - 1; /* allow for varnames row */
	}
	if (joining(c)) {
	    int i;

	    c->dset->v = 1;
	    for (i=0; i<JOIN_MAXCOL; i++) {
		if (c->jspec->colnames[i] != NULL) {
		    c->dset->v += 1;
		}
	    }
	} else if (cols_subset(c)) {
	    c->dset->v = c->cols_list[0] + 1;
	} else {
	    c->dset->v = csv_skip_column(c) ? c->ncols : c->ncols + 1;
	}
    }

#if CDEBUG
    if (joining(c)) {
	fprintf(stderr, "csv dataset dimensions: v=%d, n=%d\n",
		c->dset->v, c->dset->n);
    }
#endif
}

/*
 * real_import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @cols: column specification.
 * @rows: row specification.
 * @keyspec: join key specification.
 * @cptr: optional location to grab CSV data struct.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; for use of OPT_T see
 * the help for "append".
 * @prn: gretl printing struct (or NULL).
 * 
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

static int real_import_csv (const char *fname, 
			    DATASET *dset, 
			    const char *cols, 
			    const char *rows,
			    csvjoin *join,
			    gretlopt opt, 
			    PRN *prn)
{
    csvdata *c = NULL;
    FILE *fp = NULL;
    PRN *mprn = NULL;
    int newdata = (dset->Z == NULL);
    int popit = 0;
    long datapos;
    int i, err = 0;

    if (opt & OPT_Q) {
	/* quiet */
	prn = NULL;
    }

    if (prn != NULL) {
	set_alt_gettext_mode(prn);
    }

    if (gretl_messages_on()) {
	mprn = prn;
    }

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	pprintf(prn, A_("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto csv_bailout;
    }

    c = csvdata_new(dset);
    if (c == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }

    if (cols != NULL) {
	err = csvdata_add_cols_list(c, cols, opt);
	if (err) {
	    goto csv_bailout;
	} else if (fixed_format(c)) {
	    pprintf(mprn, A_("using fixed column format\n"));
	}
    }

    if (rows != NULL) {
	err = csvdata_add_row_mask(c, rows);
	if (err) {
	    goto csv_bailout;
	} 
    }

    if (join != NULL) {
	c->jspec = join;
    }

    if (mprn != NULL) {
	print_csv_parsing_header(fname, mprn);
    }

    /* get line length, also check for binary data, etc. */
    c->maxlen = csv_max_line_length(fp, c, prn);    
    if (c->maxlen <= 0) {
	err = E_DATA;
	goto csv_bailout;
    } 

    if (!fixed_format(c) && !csv_got_delim(c)) {
	/* set default delimiter */
	if (csv_got_tab(c)) {
	    c->delim = '\t';
	} else if (csv_got_semi(c)) {
	    c->delim = ';';
	} else {
	    c->delim = ' ';
	}
    }

    /* buffer to hold lines */
    c->line = malloc(c->maxlen);
    if (c->line == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }  

 alt_delim:

    if (!fixed_format(c)) {
	pprintf(mprn, A_("using delimiter '%c'\n"), c->delim);
    }

    pprintf(mprn, A_("   longest line: %d characters\n"), c->maxlen - 1);

    if (csv_has_trailing_comma(c) && c->delim != ',') {
	csv_unset_trailing_comma(c);
    }

    rewind(fp);

    /* read lines, check for consistency in number of fields */
    err = csv_fields_check(fp, c, mprn);
    if (err && !fixed_format(c)) {
	if (c->delim != ';' && csv_got_semi(c)) {
	    c->delim = ';';
	    err = 0;
	    goto alt_delim;
	}
	pputs(prn, A_(csv_msg));
	goto csv_bailout;
    }

    csv_set_dataset_dimensions(c);

    pprintf(mprn, A_("   number of variables: %d\n"), c->dset->v - 1);
    pprintf(mprn, A_("   number of non-blank lines: %d\n"), c->nrows);

    if (c->dset->n == 0) {
	pputs(prn, A_("Invalid data file\n"));
	err = E_DATA;
	goto csv_bailout;
    }

    /* initialize CSV dataset */
    err = start_new_Z(c->dset, 0);
    if (!err && csv_skip_column(c)) {
	err = dataset_allocate_obs_markers(c->dset);
    }

    if (err) {
	goto csv_bailout;
    }

    /* second pass */

    rewind(fp);

    if (fixed_format(c)) {
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

    if (c->decpoint == '.' && get_local_decpoint() == ',') {
	gretl_push_c_numeric_locale();
	popit = 1;
    } else if (c->decpoint == ',' && get_local_decpoint() == '.') {
	csv_set_dotsub(c);
    }

    datapos = ftell(fp);

    err = csv_read_data(c, fp, prn, mprn);

    if (!err && csv_skip_bad(c)) {
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
	reverse_data(c->dset, mprn);
    }

 csv_continue:

    c->dset->t1 = 0;
    c->dset->t2 = c->dset->n - 1;

    if (c->markerpd > 0) {
	pputs(mprn, A_("taking date information from row labels\n\n"));
	if (csv_skip_bad(c)) {
	    pprintf(prn, "WARNING: Check your data! gretl has stripped out "
		    "what appear to be\nextraneous lines in a %s dataset: " 
		    "this may not be right.\n\n",
		    (c->dset->pd == 4)? "quarterly" : "monthly");
	}
    } else {
	pputs(mprn, A_("treating these as undated data\n\n"));
	dataset_obs_info_default(c->dset);
    }

    if (c->dset->pd != 1 || strcmp(c->dset->stobs, "1")) { 
        c->dset->structure = TIME_SERIES;
    }

    if (c->st != NULL) {
	if (joining(c)) {
	    gretl_string_table_save(c->st, c->dset);
	} else {
	    gretl_string_table_print(c->st, c->dset, fname, prn);
	}
    }

    /* If there were observation labels and they were not interpretable
       as dates, and they weren't simply "1, 2, 3, ...", then they 
       should probably be preserved; otherwise discard them. 
    */
    if (c->dset->S != NULL && c->markerpd >= 0 && 
	c->dset->markers != DAILY_DATE_STRINGS) {
	dataset_destroy_obs_markers(c->dset);
    }

    if (csv_autoname(c)) {
	/* no variable names were found */
	for (i=1; i<c->dset->v; i++) {
	    sprintf(c->dset->varname[i], "v%d", i);
	}
    } else if (fix_varname_duplicates(c->dset)) {
	pputs(prn, A_("warning: some variable names were duplicated\n"));
    }

    if (!joining(c)) {
	/* not doing a special "join" operation */
	err = merge_or_replace_data(dset, &c->dset, opt, prn);

	if (!err && newdata && c->descrip != NULL) {
	    dset->descrip = c->descrip;
	    c->descrip = NULL;
	}

	if (!err) {
	    dataset_add_import_info(dset, fname, GRETL_CSV);
	}
    }

 csv_bailout:

    if (fp != NULL) {
	fclose(fp);
    }

    if (!err && c->jspec != NULL) {
	c->jspec->c = c;
    } else {
	csvdata_free(c);
    }

    if (err == E_ALLOC) {
	pputs(prn, A_("Out of memory\n"));
    }    

    return err;
}

/**
 * import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; for use of OPT_T see
 * the help for "append".
 * @prn: gretl printing struct (or NULL).
 * 
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 * 
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_csv (const char *fname, DATASET *dset, 
		gretlopt opt, PRN *prn)
{
    const char *cols = NULL;
    const char *rows = NULL;
    int ci, err;

    err = incompatible_options(opt, OPT_F | OPT_L);
    if (err) {
	/* --cols and --fixed-cols */
	return err;
    }

    ci = (dset != NULL && dset->v > 0)? APPEND : OPEN;

    if (opt & OPT_F) {
	/* we should have a "--fixed-cols=XXX" specification */
	cols = get_optval_string(ci, OPT_F);
	if (cols == NULL || *cols == '\0') {
	    return E_PARSE;
	}
    } else if (opt & OPT_L) {
	/* should have a "--cols=XXX" specification */
	cols = get_optval_string(ci, OPT_L);
	if (cols == NULL || *cols == '\0') {
	    return E_PARSE;
	}
    }	

    if (opt & OPT_M) {
	/* we should have a "--rowmask=XXX" specification */
	rows = get_optval_string(ci, OPT_M);
	if (rows == NULL || *rows == '\0') {
	    return E_PARSE;
	}
    }	

    return real_import_csv(fname, dset, cols, rows, 
			   NULL, opt, prn);
}

/* below: apparatus to implement the "join" command */

struct jr_row_ {
    int n_keys;  /* number of keys (needed for qsort callback) */
    int keyval;  /* primary key value */
    int keyval2; /* seconary key value, if applicable */
    double val;  /* data value */
    double aux;  /* auxiliary value */
};

typedef struct jr_row_ jr_row;

struct joiner_ {
    int n_rows;     /* number of rows in data table */
    int n_keys;     /* number of keys used (1 or 2) */
    int n_unique;   /* number of unique keys found, primary key */
    jr_row *rows;   /* array of table rows */
    int *keys;      /* array of unique (primary) key values as integers */
    int *key_freq;  /* counts of occurrences of (primary) key values */
    int str_keys;   /* flag for string comparison of primary keys */
    int str_keys2;  /* flag for string comparison of secondary keys */
    const int *l_keyno; /* for string comparison: list of ikey IDs in lhs dset */
    const int *r_keyno; /* for string comparison: list of okey IDs in rhs dset */
    AggrType aggr;    /* aggregation method for 1:n joining */
    int seqval;       /* aux. sequence number for aggregation */
    int auxcol;       /* aux. data column for aggregation */
    int valcol;       /* column of RHS dataset holding payload */
    DATASET *l_dset;  /* the main dataset */
    DATASET *r_dset;  /* the temporary CSV dataset */
};

typedef struct joiner_ joiner;

struct jr_filter_ {
    char *lhname;
    double lhval;
    int lhcol;
    char *rhname;
    double rhval;
    int rhcol;
    int op;
    int is_string;
};

typedef struct jr_filter_ jr_filter;

static void jr_filter_destroy (jr_filter *f)
{
    if (f != NULL) {
	free(f->lhname);
	free(f->rhname);
	free(f);
    }
}

static void joiner_destroy (joiner *jr)
{
    if (jr != NULL) {
	free(jr->rows);
	free(jr->keys);
	free(jr->key_freq);
	free(jr);
    }
}

static int join_row_wanted (DATASET *dset, int i,
			    jr_filter *filter, 
			    int *err)
{
    int ret = 0;

    if (filter == NULL) {
	return 1;
    }

    if (filter->is_string) {
	const char *x;
	const char *y;
	size_t slen;

	if (filter->lhcol) {
	    x = series_get_string_val(dset, filter->lhcol, i);
	} else {
	    x = filter->lhname;
	}

	if (filter->rhcol) {
	    y = series_get_string_val(dset, filter->rhcol, i);
	} else {
	    y = filter->rhname;
	}

	if (filter->op == B_EQ) {
	    ret = (strcmp(x, y) == 0);
	} else if (filter->op == B_GT) {
	    slen = strlen(y);
	    ret = (strlen(x) > slen) && (strncmp(x, y, slen) == 0);
	} else if (filter->op == B_LT) {
	    slen = strlen(x);
	    ret = (strlen(y) > slen) && (strncmp(x, y, slen) == 0);
	} else if (filter->op == B_NEQ) {
	    ret = (strcmp(x, y) != 0);
	} else {
	    *err = E_PARSE;
	}
    } else {
	double x = filter->lhcol ? dset->Z[filter->lhcol][i] : filter->lhval;
	double y = filter->rhcol ? dset->Z[filter->rhcol][i] : filter->rhval;

	if (filter->op == B_EQ) {
	    ret = x == y;
	} else if (filter->op == B_GT) {
	    ret = x > y;
	} else if (filter->op == B_LT) {
	    ret = x < y;
	} else if (filter->op == B_GTE) {
	    ret = x >= y;
	} else if (filter->op == B_LTE) {
	    ret = x <= y;
	} else if (filter->op == B_NEQ) {
	    ret = x != y;
	} else {
	    *err = E_PARSE;
	}
    }

#if CDEBUG
    fprintf(stderr, "join filter: %s row %d\n",
	    ret ? "keeping" : "discarding", i);
#endif

    return ret;
}

static joiner *joiner_new (csvjoin *jspec, 
			   DATASET *l_dset,
			   jr_filter *filter,
			   AggrType aggr,
			   int seqval,
			   int *err)
{
    joiner *jr = NULL;
    DATASET *r_dset = jspec->c->dset;
    int keycol = 0, key2col = 0, valcol = 0;
    int lhcol = 0, rhcol = 0, auxcol = 0;
    int i, nrows = 0;

    for (i=0; i<JOIN_MAXCOL; i++) {
	if (jspec->colnums[i] == JOIN_KEY) {
	    keycol = i+1;
	} else if (jspec->colnums[i] == JOIN_VAL) {
	    valcol = i+1;
	} else if (jspec->colnums[i] == JOIN_LHS) {
	    lhcol = i+1;
	} else if (jspec->colnums[i] == JOIN_RHS) {
	    rhcol = i+1;
	} else if (jspec->colnums[i] == JOIN_KEY2) {
	    key2col = i+1;
	} else if (jspec->colnums[i] == JOIN_AUX) {
	    auxcol = i+1;
	}
    }

    if (filter != NULL) {
	if (filter->lhname != NULL && lhcol == 0) {
	    /* a required filter column is missing */
	    fprintf(stderr, "join: filter column '%s' was not found\n", 
		    filter->lhname);
	    *err = E_DATA;
	} else {
	    filter->lhcol = lhcol;
	    filter->is_string = series_has_string_table(r_dset, lhcol);
	    if (filter->rhname != NULL) {
		if (rhcol > 0) {
		    filter->rhcol = rhcol;
		} else if (!filter->is_string) {
		    fprintf(stderr, "join: filter column '%s' was not found\n", 
			    filter->rhname);
		    *err = E_DATA;
		}
	    }
	}
    }

    if (!*err) {
	jr = malloc(sizeof *jr);
	if (jr == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NULL;
    }

    if (filter != NULL) {
	/* count the filtered rows */
	for (i=0; i<r_dset->n; i++) {
	    if (join_row_wanted(r_dset, i, filter, err)) {
		nrows++;
	    }
	}
    } else {
	nrows = r_dset->n;
    }

    jr->keys = NULL;
    jr->key_freq = NULL;
    jr->l_keyno = NULL;
    jr->r_keyno = NULL;

    jr->rows = malloc(nrows * sizeof *jr->rows);
    if (jr->rows == NULL) {
	*err = E_ALLOC;
	joiner_destroy(jr);
	jr = NULL;
    } else {
	jr->n_rows = nrows;
	jr->n_unique = 0;
	jr->aggr = aggr;
	jr->seqval = seqval;
	jr->auxcol = auxcol;
	jr->valcol = valcol;
	jr->l_dset = l_dset;
	jr->r_dset = r_dset;
    }

    if (jr != NULL) {
	/* transcribe the rows we want */
	int j = 0;

	for (i=0; i<r_dset->n; i++) {
	    if (join_row_wanted(r_dset, i, filter, err)) {
		/* the keys */
		if (keycol > 0) {
		    jr->rows[j].keyval = (int) r_dset->Z[keycol][i];
		    if (key2col > 0) {
			/* double key */
			jr->rows[j].n_keys = 2;
			jr->rows[j].keyval2 = (int) r_dset->Z[key2col][i];
		    } else {
			/* single key */
			jr->rows[j].n_keys = 1;
			jr->rows[j].keyval2 = 0;
		    }
		} else {
		    /* no keys have been specified */
		    jr->rows[j].n_keys = 0;
		    jr->rows[j].keyval = 0;
		    jr->rows[j].keyval2 = 0;
		}
		/* the "payload" data */
		if (valcol > 0) {
		    jr->rows[j].val = r_dset->Z[valcol][i];
		} else {
		    jr->rows[j].val = 0;
		}
		/* the auxiliary data */
		if (auxcol > 0) {
		    jr->rows[j].aux = r_dset->Z[auxcol][i];
		} else {
		    jr->rows[j].aux = 0;
		}
		j++;
	    }
	}
    }

    return jr;
}

static int compare_jr_rows (const void *a, const void *b)
{
    const jr_row *ra = a;
    const jr_row *rb = b;
    int ret;

    ret = (ra->keyval > rb->keyval) - (ra->keyval < rb->keyval);

    if (ret == 0 && ra->n_keys > 1) {
	ret = (ra->keyval2 > rb->keyval2) - (ra->keyval2 < rb->keyval2);
    }
    
    return ret;
}

static int joiner_sort (joiner *jr)
{
    int i, err = 0;

    qsort(jr->rows, jr->n_rows, sizeof *jr->rows, compare_jr_rows);

    jr->n_unique = 1;
    for (i=1; i<jr->n_rows; i++) {
	if (jr->rows[i].keyval != jr->rows[i-1].keyval) {
	    jr->n_unique += 1;
	}
    }

    jr->keys = malloc(jr->n_unique * sizeof *jr->keys);
    jr->key_freq = malloc(jr->n_unique * sizeof *jr->key_freq);

    if (jr->keys == NULL || jr->key_freq == NULL) {
	err = E_ALLOC;
    } else {
	int j = 0, nj = 1;

	for (i=0; i<jr->n_unique; i++) {
	    jr->key_freq[i] = 0;
	}

	jr->keys[0] = jr->rows[0].keyval;

	for (i=1; i<jr->n_rows; i++) {
	    if (jr->rows[i].keyval != jr->rows[i-1].keyval) {
		jr->keys[j] = jr->rows[i-1].keyval;
		jr->key_freq[j] = nj;
		nj = 1;
		j++;
	    } else {
		nj++;
	    }
	}

	jr->keys[j] = jr->rows[i-1].keyval;
	jr->key_freq[j] = nj;
    }

    return err;
}

#if CDEBUG > 1

static void joiner_print (joiner *jr)
{
    int i;

    fprintf(stderr, "\njoiner: n_rows = %d\n", jr->n_rows);
    for (i=0; i<jr->n_rows; i++) {
	fprintf(stderr, " row %d: keyval=%d, val=%g\n", i, jr->rows[i].keyval,
		jr->rows[i].val);
    }

    if (jr->keys != NULL) {
	fprintf(stderr, " n_unique = %d\n", jr->n_unique);
	for (i=0; i<jr->n_unique; i++) {
	    fprintf(stderr,"  key value %d : count = %d\n", 
		    jr->keys[i], jr->key_freq[i]);
	}
    }
}

#endif

static double aggr_retval (int key, const char *lstr,
			   const char **rlabels,
			   int key2, const char *lstr2,
			   const char **rlabels2,
			   joiner *jr, double *xmatch,
			   double *auxmatch, int *err)
{
    const char *rstr = NULL;
    double x, y, xa;
    int pos, imax;
    int i, n1, n;

#if CDEBUG
    if (lstr != NULL) {
	fprintf(stderr, "  left-hand (primary) key string = '%s'\n", lstr);	
    }
    if (lstr2 != NULL) {
	fprintf(stderr, "  left-hand secondary key string = '%s'\n", lstr2);	
    }    
#endif

    /* is the (primary) key present in the freq rectangle? */

    pos = -1;
    if (jr->str_keys) {
	for (i=0; i<jr->n_unique && pos<0; i++) {
	    rstr = rlabels[jr->keys[i] - 1];
	    if (!strcmp(lstr, rstr)) {
		/* got a match on strings */
		pos = i;
	    }
	}
    } else {
	for (i=0; i<jr->n_unique && pos<0; i++) {
	    if (key == jr->keys[i]) {
		/* got a numerical match */
		pos = i;
	    }
	}
    }

#if CDEBUG
    fprintf(stderr, "  pos on right, among unique keys = %d\n", pos);
#endif

    if (pos < 0) {
	/* (primary) key not found */
	return (jr->aggr == AGGR_COUNT)? 0 : NADBL;
    }

    n1 = jr->key_freq[pos];

    if (jr->n_keys == 1) {
	if (jr->aggr == AGGR_COUNT) {
	    /* simple, we're done */
	    return n1;
	} else if (jr->aggr == AGGR_SEQ && jr->seqval > n1) {
	    /* out of bounds sequence index */
	    return NADBL;
	} else if (n1 > 1 && jr->aggr == AGGR_NONE) {
	    *err = E_DATA;
	    gretl_errmsg_set(_("You need to specify an aggregation "
			       "method for a 1:n join"));
	    return NADBL;
	}
    }

    pos = -1;
    
    /* At this point jr is already sorted: we now need to
       find the starting row of the target (primary) key 
       value in the rectangle proper.
    */

    if (jr->str_keys) {
	/* matching by string values */
	for (i=0; i<jr->n_rows && pos<0; i++) {
	    rstr = rlabels[jr->rows[i].keyval - 1];
	    if (!strcmp(lstr, rstr)) {
		pos = i;
	    }
	}
    } else {
	/* matching numerically */
	for (i=0; i<jr->n_rows && pos<0; i++) {
	    if (key == jr->rows[i].keyval) {
		pos = i;
	    }
	}
    }

    if (pos < 0) {
	/* "can't happen" */
	return NADBL;
    }

    imax = pos + n1;
    n = 0; /* will hold count of non-NA matches */

    /* If we also have a secondary key, we need to find how
       many instances of the secondary match fall under the
       primary match. While we're at it fill out the array
       @xmatch with non-missing values from the relevant
       rows.
    */

    if (jr->n_keys > 1) {
	/* note: @totcount ignores the OK/NA distinction */
	int match, totcount = 0;

	for (i=pos; i<imax; i++) {
	    if (jr->str_keys2) {
		rstr = rlabels2[jr->rows[i].keyval2 - 1];
		match = (strcmp(lstr2, rstr) == 0);
	    } else {
		match = (key2 == jr->rows[i].keyval2);
	    }
	    if (match) {
		totcount++;
		x = jr->rows[i].val;
		if (jr->auxcol) {
		    xa = jr->rows[i].aux;
		    if (!na(xa)) {
			auxmatch[n] = xa;
			xmatch[n++] = x;
		    }
		} else if (!na(x)) {
		    xmatch[n++] = x;
		}
	    }		
	}

	if (jr->aggr == AGGR_COUNT) {
	    return totcount;
	} else if (totcount > 1 && jr->aggr == AGGR_NONE) {
	    *err = E_DATA;
	    gretl_errmsg_set(_("You need to specify an aggregation "
			       "method for a 1:n join"));
	    return NADBL;
	}
    } else {
	/* just one key */
	for (i=pos; i<imax; i++) {
	    x = jr->rows[i].val;
	    if (jr->auxcol) {
		xa = jr->rows[i].aux;
		if (!na(xa)) {
		    auxmatch[n] = xa;
		    xmatch[n++] = x;
		}
	    } else if (!na(x)) {
		xmatch[n++] = x;
	    }
	}
    }

    x = NADBL;

    if (n == 0) {
	; /* all obs. are NAs */
    } else if (jr->aggr == AGGR_NONE) {
	x = xmatch[0];
    } else if (jr->aggr == AGGR_SEQ) {
	i = jr->seqval - 1;
	if (i >= 0 && i < n) {
	    x = xmatch[i];
	}
    } else if (jr->aggr == AGGR_MAX) {
	if (jr->auxcol) {
	    /* using the max of the auxiliary var */
	    xa = auxmatch[0];
	}
	x = xmatch[0];
	for (i=1; i<n; i++) {
	    if (jr->auxcol) {
		y = auxmatch[i];
		if (y > xa) {
		    x = xmatch[i];
		}
	    } else {
		y = xmatch[i];
		if (y > x) {
		    x = y;
		}
	    }
	}
    } else if (jr->aggr == AGGR_MIN) {
	if (jr->auxcol) {
	    /* using the min of the auxiliary var */
	    xa = auxmatch[0];
	}
	x = xmatch[0];
	for (i=1; i<n; i++) {
	    if (jr->auxcol) {
		y = auxmatch[i];
		if (y < xa) {
		    x = xmatch[i];
		}
	    } else {	    
		y = xmatch[i];
		if (y < x) {
		    x = y;
		}
	    }
	}
    } else if (jr->aggr == AGGR_SUM || jr->aggr == AGGR_AVG) {
	x = 0.0;
	for (i=0; i<n; i++) {
	    x += xmatch[i];
	}
	if (jr->aggr == AGGR_AVG) {
	    x /= n;
	}
    }

    return x;
}

static int key_from_double (double x, int *err)
{
    if (xna(x) || fabs(x) > INT_MAX) {
	*err = E_INVARG;
	return -1;
    } else {
	return (int) trunc(x);
    }
}

/* Handle the case where (a) the "aggregated" value from the right,
   @rz, is actually the coding of a string value, and (b) the LHS
   series is pre-existing and already has a string table attached. The
   RHS coding must be made consistent with that on the left. We reach
   this function only if we've verified that there are string tables
   on both sides, and if @rz is not NA.
*/

static double maybe_adjust_string_code (series_table *rst,
					series_table *lst,
					double rz, int *err)
{
    const char *rstr = series_table_get_string(rst, rz);
    double lz = series_table_get_value(lst, rstr);

    if (!na(lz)) {
	/* use the LHS encoding */
	rz = lz;
    } else {
	/* we need to append to the LHS string table */
	int n = series_table_add_string(lst, rstr);

	if (n < 0) {
	    *err = E_ALLOC;
	} else {
	    rz = n;
	}
    }

    return rz;
}

static int aggregate_data (joiner *jr, const int *ikeyvars, int newvar)
{
    series_table *rst = NULL;
    series_table *lst = NULL;
    const char **llabels = NULL;
    const char **rlabels = NULL;
    const char **llabels2 = NULL;
    const char **rlabels2 = NULL;
    const char *keystr = NULL;
    const char *key2str = NULL;
    DATASET *dset = jr->l_dset;
    double *xmatch = NULL;
    double *auxmatch = NULL;
    double z, z2 = 0.0;
    int strcheck = 0;
    int i, nmax, key, key2 = 0;
    int err = 0;

#if CDEBUG
    fputs("\naggregate data:\n", stderr);
#endif

    if (jr->str_keys) {
	/* matching on primary key by strings */
	llabels = series_get_string_vals(jr->l_dset, jr->l_keyno[1], NULL);
	rlabels = series_get_string_vals(jr->r_dset, jr->r_keyno[1], NULL);
    }

    if (jr->str_keys2) {
	/* matching on secondary key by strings */
	llabels2 = series_get_string_vals(jr->l_dset, jr->l_keyno[2], NULL);
	rlabels2 = series_get_string_vals(jr->r_dset, jr->r_keyno[2], NULL);
    }

    /* find the greatest (primary) key frequency */
    nmax = 0;
    for (i=0; i<jr->n_unique; i++) {
	if (jr->key_freq[i] > nmax) {
	    nmax = jr->key_freq[i];
	}
    }    

    if (nmax > 0) {
	xmatch = malloc(nmax * sizeof *xmatch);
	if (xmatch == NULL) {
	    return E_ALLOC;
	}
	if (jr->auxcol) {
	    auxmatch = malloc(nmax * sizeof *auxmatch);
	    if (auxmatch == NULL) {
		free(xmatch);
		return E_ALLOC;
	    }
	}	    
    }

    if (jr->valcol > 0) {
	rst = series_get_string_table(jr->r_dset, jr->valcol);
	lst = series_get_string_table(jr->l_dset, newvar);
	strcheck = (rst != NULL && lst != NULL);
    }

    for (i=dset->t1; i<=dset->t2 && !err; i++) {
	z = dset->Z[ikeyvars[1]][i];
	if (jr->n_keys == 2) {
	    z2 = dset->Z[ikeyvars[2]][i];
	}
#if CDEBUG
	if (jr->n_keys == 2) {
	    fprintf(stderr, " left-hand key values = (%g,%g)\n", z, z2);
	} else {
	    fprintf(stderr, " left-hand key value = %g\n", z);
	}
#endif
	if (xna(z) || xna(z2)) {
	    dset->Z[newvar][i] = NADBL;
	} else {
	    key = key_from_double(z, &err);
	    if (!err && llabels != NULL) {
		keystr = llabels[key-1];
	    }
	    if (!err && jr->n_keys == 2) {
		key2 = key_from_double(z2, &err);
		if (!err && llabels2 != NULL) {
		    key2str = llabels2[key2-1];
		}
	    }
	    if (!err) {
		z = aggr_retval(key, keystr, rlabels,
				key2, key2str, rlabels2,
				jr, xmatch, auxmatch,
				&err);
	    }
#if CDEBUG
	    fprintf(stderr, " aggr_retval gives %g (err = %d)\n", z, err);
#endif
	    if (!err && strcheck && !na(z)) {
		z = maybe_adjust_string_code(rst, lst, z, &err);
	    }
	    if (!err) {
		dset->Z[newvar][i] = z;
	    }
	}
    }

    free(xmatch);
    free(auxmatch);

    return err;
}

/* for use when no keys are given */

static int join_fetch_data (joiner *jr, int newvar)
{
    series_table *rst = NULL;
    series_table *lst = NULL;
    DATASET *dset = jr->l_dset;
    double z;
    int strcheck = 0;
    int i, err = 0;

    if (jr->valcol > 0) {
	rst = series_get_string_table(jr->r_dset, jr->valcol);
	lst = series_get_string_table(dset, newvar);
	strcheck = (rst != NULL && lst != NULL);
    }

    for (i=0; i<jr->n_rows && !err; i++) {
	z = jr->rows[i].val;
	if (strcheck && !na(z)) {
	    z = maybe_adjust_string_code(rst, lst, z, &err);
	}
	if (!err) {
	    dset->Z[newvar][dset->t1 + i] = z;
	}
    }

    return err;
}

/* parse a string of the form <lhs> <op> <rhs> */

static jr_filter *make_join_filter (const char *s,
				    int *err)
{
    jr_filter *filter = NULL;
    const char *opchars = "=!><";
    char *lhs = NULL, *opstr = NULL, *rhs = NULL;
    size_t nop, len = strlen(s);
    size_t nlhs = strcspn(s, opchars);
    int op = 0;

    if (nlhs == len) {
	*err = E_PARSE;
    } else {
	lhs = gretl_strndup(s, nlhs);
	g_strstrip(lhs);
	nop = strspn(s + nlhs, opchars);
	opstr = gretl_strndup(s + nlhs, nop);
	if (nlhs + nop == len) {
	    *err = E_PARSE;
	} else {
	    rhs = gretl_strdup(s + nlhs + nop);
	    g_strstrip(rhs);
	}

#if CDEBUG
	fprintf(stderr,"lhs = '%s'\n", lhs);
	fprintf(stderr,"op = '%s'\n", opstr);
	fprintf(stderr,"rhs = '%s'\n", rhs);
#endif
    }

    if (!*err) {
	if (!strcmp(opstr, "==")) {
	    op = B_EQ;
	} else if (!strcmp(opstr, "<")) {
	    op = B_LT;
	} else if (!strcmp(opstr, ">")) {
	    op = B_GT;
	} else if (!strcmp(opstr, "<=")) {
	    op = B_LTE;
	} else if (!strcmp(opstr, ">=")) {
	    op = B_GTE;
	} else if (!strcmp(opstr, "!=")) {
	    op = B_NEQ;
	} else {
	    *err = E_PARSE;
	}
    }

    if (!*err) {
	filter = malloc(sizeof *filter);
	if (filter == NULL) {
	    *err = E_ALLOC;
	} else {
	    filter->lhname = NULL;
	    filter->rhname = NULL;
	    filter->lhval = NADBL;
	    filter->rhval = NADBL;
	    filter->lhcol = 0;
	    filter->rhcol = 0;
	}
    }

    if (!*err) {
	filter->op = op;
	if (numeric_string(lhs)) {
	    filter->lhval = dot_atof(lhs);
	} else {
	    filter->lhname = lhs;
	    lhs = NULL;
	}
	if (numeric_string(rhs)) {
	    filter->rhval = dot_atof(rhs);
	} else {
	    filter->rhname = rhs;
	    rhs = NULL;
	}
    }

    free(lhs);
    free(opstr);
    free(rhs);

    return filter;
}

/* get the series ID for the left-hand side variable, which
   may be an existing series or a new one
*/

static int get_target_varnum (const char *vname,
			      DATASET *dset,
			      int *err)
{
    int i, targ = -1;

    *err = dataset_add_series(1, dset);

    if (!*err) {
	targ = dset->v - 1;
	strcpy(dset->varname[targ], vname);
	for (i=0; i<dset->n; i++) {
	    dset->Z[targ][i] = NADBL;
	}
    }

    return targ;
}

/* Parse either one column name or two comma-separated names out of
   @s. If @s contains a comma, we accept a zero-length name on either
   the left or the right -- but not both -- as indicating that we
   should use the corresponding inner key name.
*/

static int process_outer_key (const char *s, int n_keys, 
			      char *name1, char *name2)
{
    int n_okeys = 0;
    int err = 0;

    if (strchr(s, ',') == NULL) {
	/* just one outer key */
	strncat(name1, s, CSVSTRLEN - 1);
	n_okeys = 1;
    } else {
	/* two comma-separated keys */
	int n2, n1 = strcspn(s, ",");

	if (n1 >= CSVSTRLEN) {
	    err = E_PARSE;
	} else {
	    strncat(name1, s, n1);
	    s += n1 + 1;
	    n2 = strlen(s);
	    if (n2 >= CSVSTRLEN) {
		err = E_PARSE;
	    } else {
		strncat(name2, s, n2);
	    }
	}

	if (!err && n1 == 0 && n2 == 0) {
	    /* both fields empty */
	    err = E_PARSE;
	}

	if (!err) {
	    n_okeys = 2;
	}
    }

    if (!err && n_okeys != n_keys) {
	err = E_PARSE;
    }

    return err;
}

static int numerical_aggr (joiner *jr, int aggr)
{
    if (aggr == AGGR_SUM || aggr == AGGR_AVG) {
	return 1;
    } else if (aggr == AGGR_MIN || aggr == AGGR_MAX) {
	return jr->auxcol == 0;
    } else {
	return 0;
    }
}

#define lr_mismatch(l,r) ((l > 0 && r == 0) || (r > 0 && l == 0))

/* Run some checks pertaining to the nature of the payload
   (string-valued vs numeric) in relation to the aggregation
   method specified, and the nature of the existing left-hand
   series, if any.
*/

static int join_data_type_check (joiner *jr, int targvar,
				 int aggr)
{
    int lstr = -1, rstr = -1;
    int err = 0;

    if (targvar > 0) {
	/* there's an existing LHS series */
	lstr = series_has_string_table(jr->l_dset, targvar);
	if (lstr && aggr == AGGR_COUNT) {
	    /* count values can't be mixed with strings */
	    err = E_TYPES;
	}
    }

    if (!err && jr->valcol > 0) {
	/* there's a payload variable on the right */
	rstr = series_has_string_table(jr->r_dset, jr->valcol);
	if (rstr && numerical_aggr(jr, aggr)) {
	    /* if the RHS series is string-valued, numerical
	       aggregation methods are not meaningful
	    */
	    err = E_TYPES;
	}
    }

    if (!err && lr_mismatch(lstr, rstr)) {
	/* one of (L, R) is string-valued, but not the other */
	err = E_TYPES;
    }    

    if (!err && jr->auxcol > 0) {
	/* we're using an aux. column for min/max aggregation:
	   that variable cannot be string-valued 
	*/
	if (series_has_string_table(jr->r_dset, jr->auxcol)) {
	    err = E_TYPES;
	}
    }

    return err;
}

int join_from_csv (const char *fname,
		   const char *varname,
		   DATASET *dset, 
		   const int *ikeyvars,
		   const char *okey,
		   const char *filtstr,
		   const char *data,
		   AggrType aggr,
		   int seqval,
		   const char *auxname,
		   gretlopt opt,
		   PRN *prn)
{
    csvjoin jspec = {0};
    joiner *jr = NULL;
    jr_filter *filter = NULL;
    int okeyvars[3] = {0, -1, -1};
    char okeyname1[CSVSTRLEN] = {0};
    char okeyname2[CSVSTRLEN] = {0};
    int targvar, orig_v = dset->v;
    int str_keys = 0;
    int str_keys2 = 0;
    int n_keys = 0;
    int i, err = 0;

    targvar = current_series_index(dset, varname);
    if (targvar == 0) {
	/* can't modify const */
	return E_DATA;
    }    

    if (ikeyvars != NULL) {
	n_keys = ikeyvars[0];
    }

#if CDEBUG
    fputs("*** join_from_csv:\n", stderr);
    fprintf(stderr, " filename = '%s'\n", fname);
    fprintf(stderr, " target series name = '%s'\n", varname);
    if (n_keys > 0) {
	fprintf(stderr, " inner key series = %d\n", ikeyvars[1]);
	if (n_keys == 2) {
	    fprintf(stderr, " second inner key series = %d\n", ikeyvars[2]);
	}
    }
    if (okey != NULL) {
	fprintf(stderr, " outer key = '%s'\n", okey);
    } else if (n_keys > 0) {
	fprintf(stderr, " outer key = '%s' (from inner key)\n", 
		dset->varname[ikeyvars[1]]);
	if (n_keys == 2) {
	    fprintf(stderr, " second outer key = '%s' (from inner)\n", 
		    dset->varname[ikeyvars[2]]);
	}
    }
    if (filtstr != NULL) {
	fprintf(stderr, " filter = '%s'\n", filtstr);
    }    
    if (data != NULL) {
	fprintf(stderr, " source data series = '%s'\n", data);
    } else {
	fprintf(stderr, " source data series = '%s' (from inner varname)\n", 
		varname);
    }
    fprintf(stderr, " aggregation = %d\n", aggr);
    if (auxname != NULL) {
	fprintf(stderr, " aggr auxiliary col = '%s'\n", auxname);
    }
#endif

    if (filtstr != NULL) {
	filter = make_join_filter(filtstr, &err);
	if (err) {
	    fprintf(stderr, "join: error %d processing row filter\n", err);
	}
    }

    if (!err && okey != NULL) {
	err = process_outer_key(okey, n_keys, okeyname1, okeyname2);
	if (err) {
	    fprintf(stderr, "join: error %d processing outer key(s)\n", err);
	} 
    }

    if (!err) {
	/* handle the "outer" key column, if any */
	if (*okeyname1 != '\0') {
	    jspec.colnames[0] = okeyname1;
	} else if (n_keys > 0) {
	    jspec.colnames[0] = dset->varname[ikeyvars[1]];
	}

	/* the data or "payload" column */
	if (aggr != AGGR_COUNT) {
	    if (data != NULL) {
		jspec.colnames[1] = data;
	    } else {
		jspec.colnames[1] = varname;
	    }
	}

	/* handle filter columns, if applicable */
	if (filter != NULL) {
	    jspec.colnames[2] = filter->lhname;
	    jspec.colnames[3] = filter->rhname;
	}

	/* the second outer key, if present */
	if (*okeyname2 != '\0') {
	    jspec.colnames[4] = okeyname2;
	} else if (n_keys > 1) {
	    jspec.colnames[4] = dset->varname[ikeyvars[2]];
	}

	/* the auxiliary var for aggregation, if present */
	if (auxname != NULL) {
	    jspec.colnames[5] = auxname;
	}

	err = real_import_csv(fname, dset, NULL, NULL,
			      &jspec, opt, prn);
	if (err) {
	    fprintf(stderr, "join: error %d from real_import_csv\n", err);
	}
    }

    if (!err && aggr != AGGR_COUNT) {
	/* run some sanity tests on the payload */
	int valcol = -1;

	for (i=0; i<JOIN_MAXCOL && valcol<0; i++) {
	    if (jspec.colnums[i] == JOIN_VAL) {
		valcol = i+1;
	    }
	}

	if (valcol < 0) {
	    if (data != NULL) {
		fprintf(stderr, "join: data column '%s' was not found\n", 
			jspec.colnames[1]);
		err = E_UNKVAR;
	    }
	} else if (aggr != AGGR_NONE && aggr != AGGR_SEQ) {
	    if (series_has_string_table(jspec.c->dset, valcol)) {
		/* maybe this should just be a warning? */
		fprintf(stderr, "'%s' is a string variable: aggregation type "
			"is not applicable\n", jspec.colnames[1]);
		err = E_TYPES;
	    }
	}
    }

    if (!err && jspec.colnames[0] != NULL) {
	/* check that outer key was found in the right-hand-side
	   file, and is conformable to the inner key 
	*/
	for (i=0; i<JOIN_MAXCOL && okeyvars[0]<2; i++) {
	    if (jspec.colnums[i] == JOIN_KEY) {
		okeyvars[0] += 1;
		okeyvars[1] = i + 1;
	    } else if (jspec.colnums[i] == JOIN_KEY2) {
		okeyvars[0] += 1;
		okeyvars[2] = i + 1;
	    }		
	}

	if ((jspec.colnames[0] != NULL && okeyvars[1] < 0) ||
	    (jspec.colnames[4] != NULL && okeyvars[2] < 0)) {
	    fprintf(stderr, "join: error finding outer key columns\n");
	    err = E_DATA;
	} else {
	    int lstr = series_has_string_table(dset, ikeyvars[1]);
	    int rstr = series_has_string_table(jspec.c->dset, okeyvars[1]);

	    if (lstr != rstr) {
		fprintf(stderr, "key 1: numeric/string mismatch\n");
		err = E_TYPES; 
	    } else if (lstr) {
		str_keys = 1;
	    }

	    if (!err && okeyvars[2] > 0) {
		lstr = series_has_string_table(dset, ikeyvars[2]);
		rstr = series_has_string_table(jspec.c->dset, okeyvars[2]);

		if (lstr != rstr) {
		    fprintf(stderr, "key 2: numeric/string mismatch\n");
		    err = E_TYPES; 
		} else if (lstr) {
		    str_keys2 = 1;
		}
	    }		
	}
    }

    if (!err) {
#if CDEBUG > 1
	pprintf(prn, "Data extracted from %s:\n", fname);
	printdata(NULL, NULL, jspec.c->dset, OPT_O, prn);
#endif
	jr = joiner_new(&jspec, dset, filter, aggr, seqval, &err);
	if (err) {
	    fprintf(stderr, "join: error %d from joiner_new()\n", err);
	}
    }

    if (!err) {
	/* initial check for mash-up of string data and numeric */
	err = join_data_type_check(jr, targvar, aggr);
    }

    if (!err) {
	jr->n_keys = n_keys;
	jr->str_keys = str_keys;
	jr->str_keys2 = str_keys2;
	if (jr->str_keys || jr->str_keys2) {
	    jr->l_keyno = ikeyvars;
	    jr->r_keyno = okeyvars;
	}
	if (jr->n_keys > 0) {
	    err = joiner_sort(jr);
	}	
#if CDEBUG > 1
	joiner_print(jr);
#endif
    }

    if (!err && jr->n_keys == 0 && jr->n_rows != sample_size(jr->l_dset)) {
	gretl_errmsg_set(_("Series length does not match the dataset"));
	err = E_DATA;
    }

    if (!err && targvar < 0) {
	targvar = get_target_varnum(varname, dset, &err);
	if (err) {
	    fprintf(stderr, "join: error %d from get_target_varnum()\n", err);
	}	
    }

    if (!err) {
	if (jr->n_keys == 0) {
	    err = join_fetch_data(jr, targvar);
	} else {
	    err = aggregate_data(jr, ikeyvars, targvar);
	}
	if (err) {
	    fprintf(stderr, "join: error %d from aggregate_data()\n", err);
	}
    }

    if (!err && dset->v > orig_v && jr->valcol > 0) {
	/* we got a newly added payload series */
	if (series_has_string_table(jr->r_dset, jr->valcol)) {
	    /* let the new series grab the RHS string table */
	    steal_string_table(jr->l_dset, targvar, jr->r_dset, jr->valcol); 
	}
    }

    if (err) {
	dataset_drop_last_variables(dset->v - orig_v, dset);
    }

    csvdata_free(jspec.c);
    joiner_destroy(jr);
    jr_filter_destroy(filter);

    return err;
}
