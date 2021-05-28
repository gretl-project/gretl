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
#include "uservar.h"
#include "genparse.h"
#include "gretl_xml.h"
#include "gretl_midas.h"
#include "matrix_extra.h"
#include "gretl_www.h"
#include "csvdata.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <errno.h>

#define CDEBUG 0    /* CSV reading in general */
#define AGGDEBUG 0  /* aggregation in "join" */
#define TDEBUG 0    /* handling of time keys in "join" */

#define CSVSTRLEN 128

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
    CSV_DOTSUB   = 1 << 9,
    CSV_ALLCOLS  = 1 << 10,
    CSV_BOM      = 1 << 11,
    CSV_VERBOSE  = 1 << 12,
    CSV_THOUSEP  = 1 << 13,
    CSV_NOHEADER = 1 << 14,
    CSV_QUOTES   = 1 << 15,
    CSV_AS_MAT   = 1 << 16
};

enum {
    JOIN_KEY,
    JOIN_F1,
    JOIN_F2,
    JOIN_F3,
    JOIN_KEY2,
    JOIN_AUX,
    JOIN_TARG
};

typedef struct joinspec_ joinspec;
typedef struct csvprobe_ csvprobe;
typedef struct csvdata_ csvdata;

struct joinspec_ {
    int ncols;
    const char **colnames;
    const char *mdsbase;
    int *colnums;
    int *timecols;
    csvdata *c;
    DATASET *dset;
    int wildcard;
    int auto_midas;
    int midas_pd;
    char **wildnames;
    char **mdsnames;
    char **tmpnames;
    int n_tmp;
};

struct csvprobe_ {
    DATASET *dset; /* more info might be wanted */
};

struct csvdata_ {
    int flags;
    char delim;
    char decpoint;
    char thousep;
    char qchar;
    int markerpd;
    int maxlinelen;
    int real_n;
    char *line;
    DATASET *dset;
    int ncols, nrows;
    long datapos;
    char str[CSVSTRLEN];
    char skipstr[8];
    int *codelist;
    char *descrip;
    const char *user_na;
    gretl_string_table *st;
    int *cols_list;
    int *width_list;
    const gretl_matrix *rowmask;
    int masklen;
    joinspec *jspec; /* info used for "join" command */
    csvprobe *probe; /* used in connection with "join" */
};

#define csv_has_trailing_comma(c) (c->flags & CSV_TRAIL)
#define csv_has_obs_column(c)     (c->flags & CSV_OBS1)
#define csv_has_blank_column(c)   (c->flags & CSV_BLANK1)
#define csv_got_tab(c)            (c->flags & CSV_GOTTAB)
#define csv_got_semi(c)           (c->flags & CSV_GOTSEMI)
#define csv_got_delim(c)          (c->flags & CSV_GOTDELIM)
#define csv_autoname(c)           (c->flags & CSV_AUTONAME)
#define csv_skip_col_1(c)         (c->flags & (CSV_OBS1 | CSV_BLANK1))
#define csv_have_data(c)          (c->flags & CSV_HAVEDATA)
#define csv_data_reversed(c)      (c->flags & CSV_REVERSED)
#define csv_do_dotsub(c)          (c->flags & CSV_DOTSUB)
#define csv_all_cols(c)           (c->flags & CSV_ALLCOLS)
#define csv_has_bom(c)            (c->flags & CSV_BOM)
#define csv_is_verbose(c)         (c->flags & CSV_VERBOSE)
#define csv_scrub_thousep(c)      (c->flags & CSV_THOUSEP)
#define csv_no_header(c)          (c->flags & CSV_NOHEADER)
#define csv_keep_quotes(c)        (c->flags & CSV_QUOTES)
#define csv_as_matrix(c)          (c->flags & CSV_AS_MAT)

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
#define csv_set_all_cols(c)         (c->flags |= CSV_ALLCOLS)
#define csv_set_has_bom(c)          (c->flags |= CSV_BOM)
#define csv_set_verbose(c)          (c->flags |= CSV_VERBOSE)
#define csv_set_scrub_thousep(c)    (c->flags |= CSV_THOUSEP)
#define csv_set_no_header(c)        (c->flags |= CSV_NOHEADER)
#define csv_unset_keep_quotes(c)    (c->flags &= ~CSV_QUOTES)
#define csv_set_as_matrix(c)        (c->flags |= CSV_AS_MAT)

#define csv_skip_bad(c)        (*c->skipstr != '\0')
#define csv_has_non_numeric(c) (c->st != NULL)

#define fixed_format(c) (c->cols_list != NULL && c->width_list != NULL)
#define cols_subset(c) (c->cols_list != NULL && c->width_list == NULL)
#define rows_subset(c) (c->rowmask != NULL)

#define joining(c) (c->jspec != NULL)
#define probing(c) (c->probe != NULL)

#define is_wildstr(s) (strchr(s, '*') || strchr(s, '?'))

static int
time_series_label_check (DATASET *dset, int reversed, char *skipstr,
			 int convert_pd, PRN *prn);

static int format_uses_quarterly (char *fmt);

/* file-scope global */
static char import_na[8];

struct time_mapper {
    int ncols;         /* number of "timeconv" columns */
    char **colnames;   /* array of outer-dataset column names */
    char *tname;       /* the name of the "tkey", if among colnames, or NULL */
    char **fmt;        /* array of up to two time-format strings, or NULL */
    char m_means_q[2]; /* array of "monthly means quarterly" flags */
};

/* file-scope global */
struct time_mapper tconv_map;

enum {
    TCONV_FMT = 0,
    TKEY_FMT = 1
};

#define no_formats(map) (map.fmt == NULL)
#define no_tkey_format(map) (map.tname == NULL)
#define has_tconv_format(map) (map.fmt[TCONV_FMT] != NULL)
#define is_tkey_variable(name, map) (strcmp(name, map.tname) == 0)

static void timeconv_map_set (int ncols, char **colnames,
			      char *tname, char **fmt)
{
    tconv_map.ncols = ncols;
    tconv_map.colnames = colnames;
    tconv_map.tname = tname;
    tconv_map.fmt = fmt;

    if (fmt != NULL) {
	if (fmt[TCONV_FMT] != NULL) {
	    tconv_map.m_means_q[TCONV_FMT] =
		format_uses_quarterly(fmt[TCONV_FMT]);
	}
	if (fmt[TKEY_FMT] != NULL) {
	    tconv_map.m_means_q[TKEY_FMT] =
		format_uses_quarterly(fmt[TKEY_FMT]);
	}
    }
}

static void timeconv_map_init (void)
{
    timeconv_map_set(0, NULL, NULL, NULL);
}

static void timeconv_map_destroy (void)
{
    if (tconv_map.colnames != NULL) {
	strings_array_free(tconv_map.colnames, tconv_map.ncols);
    }
    if (tconv_map.fmt != NULL) {
	strings_array_free(tconv_map.fmt, 2);
    }
    timeconv_map_init();
}

static int timecol_get_format (const DATASET *dset, int v,
			       char **pfmt, int *q)
{
    if (no_formats(tconv_map)) {
	return 0;
    } else if (no_tkey_format(tconv_map)) {
	/* get the common "tconvert" format */
	*pfmt = tconv_map.fmt[TCONV_FMT];
	*q = tconv_map.m_means_q[TCONV_FMT];
	return 1;
    } else if (is_tkey_variable(dset->varname[v], tconv_map)) {
	/* get the tkey-specific format */
	*pfmt = tconv_map.fmt[TKEY_FMT];
	*q = tconv_map.m_means_q[TKEY_FMT];
	return 1;
    } else if (has_tconv_format(tconv_map)) {
	/* get the other one */
	*pfmt = tconv_map.fmt[TCONV_FMT];
	*q = tconv_map.m_means_q[TCONV_FMT];
	return 1;
    }

    return 0;
}

static int column_is_timecol (const char *colname)
{
    int i, n = tconv_map.ncols;

    for (i=0; i<n; i++) {
	if (!strcmp(colname, tconv_map.colnames[i])) {
	    return 1;
	}
    }

    return 0;
}

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

    c->flags = CSV_QUOTES;
    c->delim = '\t';
    c->thousep = 0;
    c->qchar = 0;
    c->markerpd = -1;
    c->maxlinelen = 0;
    c->real_n = 0;
    c->line = NULL;
    c->dset = NULL;
    c->ncols = 0;
    c->nrows = 0;
    c->datapos = 0;
    *c->str = '\0';
    *c->skipstr = '\0';
    c->codelist = NULL;
    c->descrip = NULL;
    c->user_na = NULL;
    c->st = NULL;
    c->cols_list = NULL;
    c->width_list = NULL;
    c->rowmask = NULL;
    c->masklen = 0;

    if (strcmp(import_na, "default")) {
	c->user_na = import_na;
    }

    c->jspec = NULL;
    c->probe = NULL;

    c->dset = datainfo_new();

    if (c->dset == NULL) {
	free(c);
	c = NULL;
    } else {
	c->delim = get_data_export_delimiter();
	c->decpoint = get_data_export_decpoint();
	if (dset != NULL && dset->Z != NULL) {
	    c->flags |= CSV_HAVEDATA;
	}
#if CDEBUG
	fprintf(stderr, "csvdata_new: c->delim = '%c', c->decpoint = '%c'\n",
		c->delim, c->decpoint);
#endif
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

    if (get_matrix_by_name(s)) {
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

    err = dataset_add_observations(dset, add, OPT_A);

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
    guint32 ed1, ed2;
    int nmiss = 0, err = 0;

    *pd = 0;

    ed1 = get_epoch_day(lbl1);
    ed2 = get_epoch_day(lbl2);
    if (ed1 <= 0 || ed2 <= 0) {
	err = 1;
    }

#if DAY_DEBUG
    fprintf(stderr, "check_daily_dates: '%s' -> %d, '%s' -> %d\n",
	    lbl1, (int) ed1, lbl2, (int) ed2);
#endif

    dset->pd = guess_daily_pd(dset);
    dset->structure = TIME_SERIES;

#if DAY_DEBUG
    fprintf(stderr, "guessed at daily pd = %d\n", dset->pd);
#endif

    if (!err) {
	if (ed2 < ed1) {
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
	guint32 n1 = (*reversed)? ed2 : ed1;
	guint32 n2 = (*reversed)? ed1 : ed2;

	fulln = n2 - n1 + 1;

	if (T > fulln) {
	    err = 1;
	} else {
	    nmiss = fulln - T;
	    pprintf(prn, _("Observations: %d; days in sample: %d\n"),
		    T, fulln);
	    if (nmiss > 300 * T) {
		pprintf(prn, _("Probably annual data\n"));
		*pd = 1;
	    } else if (nmiss > 50 * T) {
		pprintf(prn, _("Probably quarterly data\n"));
		*pd = 4;
	    } else if (nmiss > 20 * T) {
		pprintf(prn, _("Probably monthly data\n"));
		*pd = 12;
	    } else if (nmiss > 3 * T) {
		pprintf(prn, _("Probably weekly data\n"));
		*pd = dset->pd = 52;
	    } else {
		pprintf(prn, _("Missing daily rows: %d\n"), nmiss);
	    }
	}
    }

    nbak = 0;

    for (t=0; t<dset->n && !err; t++) {
	int wd, s = (*reversed)? (dset->n - 1 - t) : t;

	wd = weekday_from_date(dset->S[s]);

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
		    t, dset->S[s], n, t);
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
    fprintf(stderr, "check_daily_dates: daily pd = %d, reversed = %d, err = %d\n",
	    dset->pd, *reversed, err);
#endif

    return (err)? -1 : dset->pd;
}

/* convert from daily date label to a lower frequency --
   annual, monthly or quarterly -- if @pd indicates this
   is required
*/

static void convert_daily_label (char *targ, const char *src,
				 int pd)
{
    int y, m, d;

    sscanf(src, YMD_READ_FMT, &y, &m, &d);

    if (pd == 1) {
	sprintf(targ, "%d", y);
    } else if (pd == 12) {
	sprintf(targ, "%d:%02d", y, m);
    } else if (pd == 4) {
	sprintf(targ, "%d:%d", y, m / 3 + (m % 3 != 0));
    }
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

static int consistent_qm_labels (DATASET *dset, int reversed,
				 int convert_pd, char *skipstr,
				 int *ppd, const char *fmt,
				 int *extra_zero, PRN *prn)
{
    char bad[16], skip[8];
    char label[OBSLEN];
    int Ey; /* expected year */
    int Ep; /* expected sub-period */
    int t, s, yr, per;
    int pmin = 1;
    int pd, pd0;
    int ret = 1;

    pd = pd0 = *ppd;

 restart:

    s = reversed ? (dset->n - 1) : 0;

    if (convert_pd) {
	convert_daily_label(label, dset->S[s], pd);
    } else {
	strcpy(label, dset->S[s]);
    }

    if (sscanf(label, fmt, &yr, &per) != 2) {
	return 0;
    }

    for (t=1; t<dset->n; t++) {
	s = (reversed)? (dset->n - 1 - t) : t;
	Ey = (per == pd)? yr + 1 : yr;
	Ep = (per == pd)? pmin : per + pmin;

	if (convert_pd) {
	    convert_daily_label(label, dset->S[s], pd);
	} else {
	    strcpy(label, dset->S[s]);
	}

	if (sscanf(label, fmt, &yr, &per) != 2) {
	    ret = 0;
	} else if (Ep == 1 && pd == pd0 && per == pd + 1
		   && skipstr != NULL) {
	    *skip = *bad = '\0';
	    strncat(skip, label + 4, 7);
	    strncat(bad, label, OBSLEN-1);
	    pd = pd0 + 1;
	    goto restart;
	} else if (per == Ep + 2 && pmin == 1 && fakequarter(per)) {
	    *bad = '\0';
	    strncat(bad, label, OBSLEN-1);
	    pmin = 3;
	    goto restart;
	} else if (pd == 12 && Ep == 5 && per == 1 && yr == Ey + 1) {
	    /* apparently monthly but really quarterly? */
	    pprintf(prn, "   \"%s\": quarterly date with spurious zero?\n", label);
	    *extra_zero = 1;
	    *ppd = pd0 = pd = 4;
	    goto restart;
	} else if (yr != Ey || per != Ep) {
	    ret = 0;
	}

	if (!ret) {
	    pprintf(prn, "   %s: not a consistent date\n", label);
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

static int consistent_year_labels (const DATASET *dset,
				   int reversed,
				   int convert_pd)
{
    char label[OBSLEN];
    int s, t, yr, yprev;
    int ret = 1;

    s = (reversed)? (dset->n - 1) : 0;
    yprev = atoi(dset->S[s]);

    for (t=1; t<dset->n; t++) {
	s = reversed ? (dset->n - 1 - t) : t;
	if (convert_pd) {
	    convert_daily_label(label, dset->S[s], 1);
	    yr = atoi(label);
	} else {
	    yr = atoi(dset->S[s]);
	}
	if (yr != yprev + 1) {
	    ret = 0;
	    break;
	}
	yprev = yr;
    }

    return ret;
}

/* check for all 1s in first column of dates: this may
   indicate start-of-period dates, day first */

static int all_day_ones (DATASET *dset)
{
    int t;

    for (t=1; t<dset->n; t++) {
	if (atoi(dset->S[t]) != 1) {
	    return 0;
	} else if (t > 31) {
	    /* "1" can't mean January */
	    return 1;
	}
    }

    return 0;
}

enum date_orders {
    YYYYMMDD = 1,
    MMDDYYYY,
    DDMMYYYY
};

static int get_date_order (int f0, int fn, DATASET *dset)
{
    if (f0 > 31 || fn > 31) {
	/* first field must be year */
	return YYYYMMDD;
    } else if (f0 > 12 || fn > 12) {
	/* first field must be day */
	return DDMMYYYY;
    } else if (f0 == 1 && fn == 1 && all_day_ones(dset)) {
	/* start-of-period dates, day first? */
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
	sscanf(dset->S[t], YMD_READ_FMT, &y, &d, &m);
	sprintf(dset->S[t], YMD_WRITE_FMT, d, m, y);
    }
}

static int transform_daily_dates (DATASET *dset, int dorder,
				  char sep)
{
    char *label, fmt[16];
    int t, yr, mon, day;
    int n, err = 0;

    if (sep > 0) {
	sprintf(fmt, "%%d%c%%d%c%%d", sep, sep);
    } else {
	strcpy(fmt, "%4d%2d%2d");
    }

    for (t=0; t<dset->n && !err; t++) {
	label = dset->S[t];
	if (dorder == YYYYMMDD) {
	    n = sscanf(label, fmt, &yr, &mon, &day);
	} else if (dorder == DDMMYYYY) {
	    n = sscanf(label, fmt, &day, &mon, &yr);
	} else {
	    n = sscanf(label, fmt, &mon, &day, &yr);
	}
	if (n == 3) {
	    sprintf(label, YMD_WRITE_Y2_FMT, yr, mon, day);
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

    pprintf(prn, _("reversing the data!\n"));

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
    char s1 = 0, s2 = 0;
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[dset->n - 1];
    int dorder = 0;

    if ((sscanf(lbl1, "%d%c%d%c%d", &d1[0], &s1, &d1[1], &s2, &d1[2]) == 5 &&
	 sscanf(lbl2, "%d%c%d%c%d", &d2[0], &s1, &d2[1], &s2, &d2[2]) == 5 &&
	 s1 == s2 && ispunct(s1)) ||
	(sscanf(lbl1, "%4d%2d%2d", &d1[0], &d1[1], &d1[2]) == 3 &&
	 sscanf(lbl2, "%4d%2d%2d", &d2[0], &d2[1], &d2[2]) == 3)) {
	int mon1, day1;
	int mon2, day2;
	int pd, ret = 0;

	dorder = get_date_order(d1[0], d2[0], dset);

    tryagain:

	if (dorder == YYYYMMDD) {
	    pputs(prn, _("Trying date order YYYYMMDD\n"));
	    mon1 = d1[1];
	    day1 = d1[2];
	    mon2 = d2[1];
	    day2 = d2[2];
	} else if (dorder == DDMMYYYY) {
	    pputs(prn, _("Trying date order DDMMYYYY\n"));
	    day1 = d1[0];
	    mon1 = d1[1];
	    day2 = d2[0];
	    mon2 = d2[1];
	} else {
	    pputs(prn, _("Trying date order MMDDYYYY\n"));
	    mon1 = d1[0];
	    day1 = d1[1];
	    mon2 = d2[0];
	    day2 = d2[1];
	}

	if (mon1 > 0 && mon1 < 13 &&
	    mon2 > 0 && mon2 < 13 &&
	    day1 > 0 && day1 < 32 &&
	    day2 > 0 && day2 < 32) {
	    /* looks promising for calendar dates, but check
	       further if we don't have the canonical order
	       or separator
	    */
	    if (dorder != YYYYMMDD || s1 != '-') {
		if (transform_daily_dates(dset, dorder, s1)) {
		    return -1;
		}
		s1 = '-';
	    }
	    pprintf(prn, _("Could be %s - %s\n"), lbl1, lbl2);
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
		    int convert_pd = 0;

		    if (pd == 1 || pd == 4 || pd == 12) {
			convert_pd = pd;
		    }
		    ret = time_series_label_check(dset,
						  *reversed,
						  skipstr,
						  convert_pd,
						  prn);
		    if (ret < 0 && dorder == MMDDYYYY) {
			retransform_daily_dates(dset);
			dorder = DDMMYYYY;
			goto tryagain;
		    }
		}
	    }
	    return ret;
	}
    } else {
	pprintf(prn, _("'%s' and '%s': couldn't get dates\n"), lbl1, lbl2);
    }

    return -1;
}

static int pd_from_date_label (const char *lbl, char *year, char *subp,
			       char *format, PRN *prn)
{
    const char *subchars = ".:QqMmPp-";
    int len = strlen(lbl);
    int try, pd = -1;

    strncat(year, lbl, 4);
    try = atoi(year);

    if (try > 0 && try < 3000) {
	pprintf(prn, _("   %s: probably a year... "), year);
    } else {
	pprintf(prn, _("   %s: probably not a year\n"), year);
    }

    if (len == 5) {
	pputs(prn, _("   but I can't make sense of the extra bit\n"));
    } else if (len == 4) {
	pputs(prn, _("and just a year\n"));
	pd = 1;
    } else {
	char sep = lbl[4];
	char sub[3], *s = NULL;
	int dashQ = 0;
	int p;

	if (strchr(subchars, sep)) {
	    *sub = '\0';
	    strncat(sub, lbl + 5, 2);
	    s = sub;
	    if (len == 6 || (len == 7 && (sep == 'q' || sep == 'Q'))) {
		if (len == 7) s++;
		p = atoi(s);
		if (p > 0 && p < 5) {
		    pprintf(prn, _("quarter %s?\n"), s);
		    pd = 4;
		} else {
		    pprintf(prn, "quarter %d: not possible\n", p);
		}
	    } else if (len == 7) {
		if (*s == 'Q') {
		    /* YYYY-Qn? This is supported by SDMX */
		    dashQ = 1;
		    s++;
		}
		p = atoi(s);
		if (dashQ) {
		    if (p > 0 && p < 5) {
			pprintf(prn, _("quarter %d?\n"), p);
			pd = 4;
		    } else {
			pprintf(prn, "quarter %d: not possible\n", p);
		    }
		} else {
		    if (p > 0 && p < 13) {
			pprintf(prn, _("month %s?\n"), s);
			pd = 12;
		    } else {
			pprintf(prn, "month %d: not possible\n", p);
		    }
		}
	    }
	    strcpy(subp, s);
	    if (format != NULL && (pd == 4 || pd == 12)) {
		if (dashQ) {
		    sprintf(format, "%%d%cQ%%d", sep);
		} else {
		    sprintf(format, "%%d%c%%d", sep);
		}
	    }
	}
    }

    return pd;
}

static int time_series_label_check (DATASET *dset, int reversed,
				    char *skipstr, int convert_pd,
				    PRN *prn)
{
    char year[5], sub[3];
    char format[8] = {0};
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[dset->n - 1];
    char *label;
    int pd = -1;

    *year = *sub = '\0';
    label = reversed ? lbl2 : lbl1;

    if (convert_pd) {
	char altobs[OBSLEN];

	convert_daily_label(altobs, label, convert_pd);
	pd = pd_from_date_label(altobs, year, sub, format, prn);
    } else {
	pd = pd_from_date_label(label, year, sub, format, prn);
    }

    if (pd == 1) {
	if (consistent_year_labels(dset, reversed, convert_pd)) {
	    dset->pd = pd;
	    strcpy(dset->stobs, year);
	    dset->sd0 = atof(dset->stobs);
	    strcpy(dset->endobs, lbl2);
	    dset->structure = TIME_SERIES;
	} else {
	    pputs(prn, _("   but the dates are not complete and consistent\n"));
	    pd = -1;
	}
    } else if (pd == 4 || pd == 12) {
	int savepd = pd;
	int extra_zero = 0;

	if (consistent_qm_labels(dset, reversed, convert_pd,
				 skipstr, &pd, format,
				 &extra_zero, prn)) {
	    dset->pd = pd;
	    if (savepd == 12 && pd == 4) {
		/* we switched the interpretation from
		   monthly to quarterly */
		int s;

		if (extra_zero) {
		    /* e.g. 1960Q1 written as 1960:01 */
		    s = atoi(sub + 1);
		} else {
		    /* e.g. 1960Q1 written as 1960:03 */
		    s = atoi(sub) / 3;
		}
		sprintf(dset->stobs, "%s:%d", year, s);
	    } else {
		sprintf(dset->stobs, "%s:%s", year, sub);
	    }
	    dset->sd0 = obs_str_to_double(dset->stobs);
	    ntolabel(dset->endobs, dset->n - 1, dset);
	} else {
	    pputs(prn, _("   but the dates are not complete and consistent\n"));
	    pd = -1;
	}
    }

    return pd;
}

static int dates_maybe_reversed (const char *s1,
				 const char *s2,
				 PRN *prn)
{
    char d1[5], d2[5];
    int ret = 0;

    *d1 = *d2 = '\0';

    strncat(d1, s1, 4);
    strncat(d2, s2, 4);

    ret = atoi(d1) > atoi(d2);

    if (ret) {
	pputs(prn, _("   dates are reversed?\n"));
    }

    return ret;
}

/* e.g. "M1 1957", "M12 2009" */

static int fix_IFS_data_labels (DATASET *dset)
{
    char *s1 = dset->S[0];
    char *s2 = dset->S[dset->n - 1];
    int ret = 0;

    if ((*s1 == 'M' || *s1 == 'Q') && *s2 == *s1) {
	int n1 = strlen(s1);
	int n2 = strlen(s2);

	if ((n1 == 7 || n1 == 8) && (n2 == 7 || n2 == 8) &&
	    isdigit(s1[1]) && isdigit(s2[1])) {
	    int pmax = (*s1 == 'M')? 12 : 4;
	    char c, tmp[8], *s;
	    int y, p, pbak = 0;
	    int i, n, doit = 1;

	    for (i=0; i<dset->n; i++) {
		s = dset->S[i];
		n = strlen(s);
		if (n != 7 && n != 8) {
		    doit = 0;
		    break;
		}
		n = sscanf(s, "%c%d %d", &c, &p, &y);
		if (n != 3 || c != *s1) {
		    doit = 0;
		    break;
		}
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
		    sscanf(s, "%c%d %d", &c, &p, &y);
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

static int month_number (char *s)
{
    const char *mo[] = {
	"jan", "feb", "mar", "apr",
	"may", "jun", "jul", "aug",
	"sep", "oct", "nov", "dec"
    };
    int i;

    gretl_lower(s);

    for (i=0; i<12; i++) {
	if (!strcmp(s, mo[i])) {
	    return i+1;
	}
    }

    return 0;
}

/* e.g. "Jan-1980", for monthly or quarterly data */

static int fix_mon_year_labels (DATASET *dset)
{
    char *s1 = dset->S[0];
    char *s2 = dset->S[dset->n - 1];
    char m1[4] = {0};
    char m2[4] = {0};
    int yr1 = 0, yr2 = 0;
    int ret = 0;

    if (strlen(s1) == 8 && strlen(s2) == 8 &&
	s1[3] == '-' && s2[3] == '-') {
	yr1 = atoi(s1 + 4);
	yr2 = atoi(s2 + 4);
	strncat(m1, s1, 3);
	strncat(m2, s2, 3);
    }

    if (yr1 > 999 && yr1 < 3000 && yr2 > 999 && yr2 < 3000 &&
	month_number(m1) && month_number(m2)) {
	int i, p, pbak = 0;
	int dt, pd = 0;
	char *s;

	for (i=0; i<dset->n; i++) {
	    s = dset->S[i];
	    if (strlen(s) != 8 || s[3] != '-') {
		pd = 0;
		break;
	    }
	    yr1 = atoi(s + 4);
	    *m1 = '\0';
	    strncat(m1, s, 3);
	    if (yr1 < 1000 || yr1 >= 3000 ||
		(p = month_number(m1)) < 1) {
		pd = 0;
		break;
	    }
	    if (i > 0) {
		dt = p - pbak;
		if (dt != 1 && dt != 3 && p != 1) {
		    pd = 0;
		    break;
		}
		if (pd == 0 && dt > 0) {
		    pd = (dt == 1)? 12 : 4;
		}
	    }
	    pbak = p;
	}

	if (pd > 0) {
	    for (i=0; i<dset->n; i++) {
		s = dset->S[i];
		yr1 = atoi(s + 4);
		*m1 = '\0';
		strncat(m1, s, 3);
		p = month_number(m1);
		if (pd == 12) {
		    sprintf(dset->S[i], "%d:%02d", yr1, p);
		} else {
		    sprintf(dset->S[i], "%d:%g", yr1, ceil((3+p)/4.0));
		}
	    }
	    ret = 1;
	}
    }

    return ret;
}

/* Attempt to parse CSV row labels as dates.  Return -1 if this
   doesn't work out, or 0 if the labels seem to be just integer
   observation numbers, else return the inferred data frequency.
*/

int test_markers_for_dates (DATASET *dset, int *reversed,
			    char *skipstr, PRN *prn)
{
    char endobs[OBSLEN];
    int n = dset->n;
    char *lbl1 = dset->S[0];
    char *lbl2 = dset->S[n - 1];
    int len1 = strlen(lbl1);
    int len2 = strlen(lbl2);
    int pd = -1;

    if (skipstr != NULL && *skipstr != '\0') {
	return time_series_label_check(dset, *reversed, skipstr, 0, prn);
    }

    pprintf(prn, _("   first row label \"%s\", last label \"%s\"\n"),
	    lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs, "%d", n);
    if (!strcmp(lbl1, "1") && !strcmp(lbl2, endobs)) {
	return 0;
    }

    if (fix_IFS_data_labels(dset) || fix_mon_year_labels(dset)) {
	lbl1 = dset->S[0];
	lbl2 = dset->S[n - 1];
	len1 = strlen(lbl1);
    }

    /* labels are of different lengths? */
    if (len1 != len2) {
	if (abs(len1 - len2) > 1) {
	    return -1;
	} else if (len2 > len1) {
	    len1 = len2;
	}
    }

    pputs(prn, _("trying to parse row labels as dates...\n"));

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
	    pd = time_series_label_check(dset, *reversed, skipstr, 0, prn);
	} else {
	    pputs(prn, _("   definitely not a four-digit year\n"));
	}
    }

    if (pd <= 0 && *reversed) {
	/* give up the "reversed" notion if we didn't get
	   a workable time-series interpretation */
	*reversed = 0;
    }

    return pd;
}

static int utf8_ok (gzFile fp, int pos)
{
    long mark = gztell(fp);
    int len = pos + 9;
    char *test = malloc(len + 1);
    int i, ret = 0;

    gzseek(fp, mark - pos - 1, SEEK_SET);

    for (i=0; i<len; i++) {
	test[i] = gzgetc(fp);
    }
    test[i] = '\0';

    if (g_utf8_validate(test, -1, NULL)) {
	ret = 1;
    } else {
	GError *gerr = NULL;
	gsize wrote = 0;
	gchar *tr;

	/* try for iso-8859? */
	tr = g_convert(test, -1, "UTF-8", "ISO-8859-15",
		       NULL, &wrote, &gerr);
	if (gerr != NULL) {
	    g_error_free(gerr);
	} else {
	    g_free(tr);
	    ret = 1;
	}
    }

    free(test);

    gzseek(fp, mark, SEEK_SET);

    return ret;
}

enum {
    UTF_8 = 1,
    UTF_16,
    UTF_32
};

/* If we got a UTF-16 or UTF-32 BOM, try recoding to
   UTF-8 before parsing data. We write the recoded text
   to a temporary file in the user's "dotdir" (and
   then delete that file once we're done).
*/

static int csv_recode_input (gzFile *fpp,
			     const char *fname,
			     gchar **pfname,
			     int ucode,
			     PRN *prn)
{
    const gchar *from_set =
	(ucode == UTF_32)? "UTF-32" : "UTF-16";
    gchar *altname = NULL;
    int err = 0;

    /* the current stream is not useable as is,
       so shut it down
    */
    gzclose(*fpp);
    *fpp = NULL;

    /* we'll recode to a temp file in dotdir */
    altname = g_strdup_printf("%srecode_tmp.u8", gretl_dotdir());

    err = gretl_recode_file(fname, altname,
			    from_set, "UTF-8",
			    prn);

    if (!err) {
	/* try reattaching the stream */
	*fpp = gretl_gzopen(altname, "rb");
	if (*fpp == NULL) {
	    gretl_remove(altname);
	    err = E_FOPEN;
	} else {
	    pputs(prn, "switched to recoded input\n");
	    *pfname = altname;
	    altname = NULL;
	}
    }

    g_free(altname);

    return err;
}

/* Check the first 4 bytes of "CSV" input for a Byte Order
   Mark. If we find the UTF-8 BOM (typically written by
   Microsoft tools), simply record the fact so that we can
   skip it on reading. But if we find a BOM indicating a
   16-bit or 32-bit unicode encoding, flag this by returning
   a non-zero @ucode value; in that case we'll attempt a
   full recording of the input (via GLib) before we start
   reading data.
*/

static int csv_unicode_check (gzFile fp, csvdata *c, PRN *prn)
{
    unsigned char b[4];
    int n = gzread(fp, b, 4);
    int ucode = 0;

    if (n == 4) {
	if (b[0] == 0xEF && b[1] == 0xBB && b[2] == 0xBF) {
	    pputs(prn, "got UTF-8 BOM\n");
	    ucode = UTF_8;
	} else if (b[0] == 0xFE && b[1] == 0xFF) {
	    pputs(prn, "got UTF-16BE, will try recoding\n");
	    ucode = UTF_16;
	} else if (b[0] == 0xFF && b[1] == 0xFE) {
	    if (b[2] == 0 && b[3] == 0) {
		pputs(prn, "got UTF-32LE, will try recoding\n");
		ucode = UTF_32;
	    } else {
		pputs(prn, "got UTF-16LE, will try recoding\n");
		ucode = UTF_16;
	    }
	} else if (b[0] == 0 && b[1] == 0 &&
		   b[0] == 0xFE && b[1] == 0xFF) {
	    pputs(prn, "got UTF-32BE, will try recoding\n");
	    ucode = UTF_32;
	}
    }

    if (ucode == UTF_8) {
	csv_set_has_bom(c);
	gzseek(fp, 3, SEEK_SET);
	ucode = 0;
    } else {
	gzrewind(fp);
    }

    return ucode;
}

/* The function below checks for the maximum line length in the given
   file.  It also checks for extraneous binary data (the file is
   supposed to be plain text), and checks whether the 'delim'
   character is present in the file, on a non-comment line (where
   a comment line is one that starts with '#').

   In addition, we check whether the file has a trailing comma on every
   line, and for the numbers of double- and single-quote characters
   to try to determine which, if either, is used to indicate quoted
   fields in the input.
*/

static int csv_max_line_length (gzFile fp, csvdata *cdata, PRN *prn)
{
    int c, c1, cbak = 0, cc = 0;
    int comment = 0, maxlinelen = 0;
    int max_ldquo = 0, max_lsquo = 0;
    int min_ldquo = 0, min_lsquo = 0;
    int ldquo = 0, lsquo = 0;
    int ndquo = 0, nsquo = 0;
    int crlf = 0, lines = 0;

    csv_set_trailing_comma(cdata); /* just provisionally */

    while ((c = gzgetc(fp)) != EOF) {
	if (c == 0x0d) {
	    /* CR */
	    c1 = gzgetc(fp);
	    if (c1 == EOF) {
		break;
	    } else if (c1 == 0x0a) {
		/* CR + LF -> LF */
		crlf = 1;
		c = c1;
	    } else {
		/* Mac-style: CR not followed by LF */
		c = 0x0a;
		gzungetc(c1, fp);
	    }
	}
	if (c == 0x0a) {
	    if (cc > maxlinelen) {
		maxlinelen = cc;
	    }
	    cc = 0;
	    if (cbak != 0 && cbak != ',') {
		csv_unset_trailing_comma(cdata);
	    }
	    lines++;
	    if (ldquo > max_ldquo) {
		max_ldquo = ldquo;
	    } else if (ldquo > 0 && ldquo < max_ldquo) {
		min_ldquo = ldquo;
	    }
	    if (lsquo > max_lsquo) {
		max_lsquo = lsquo;
	    } else if (lsquo > 0 && lsquo < max_lsquo) {
		min_lsquo = lsquo;
	    }
	    ldquo = lsquo = 0;
	    continue;
	}
	cbak = c;
	if (!isspace((unsigned char) c) && !isprint((unsigned char) c) &&
	    !(c == CTRLZ) && !utf8_ok(fp, cc)) {
	    pprintf(prn, _("Binary data (%d) encountered (line %d:%d): "
			   "this is not a valid text file\n"),
		    c, lines + 1, cc + 1);
	    return -1;
	}
	if (cc == 0) {
	    comment = (c == '#');
	}
	if (!comment) {
	    if (c == '\t') {
		/* let's ignore trailing tabs in this heuristic */
		c1 = gzgetc(fp);
		if (c1 != 0x0d && c1 != 0x0a) {
		    csv_set_got_tab(cdata);
		}
		gzungetc(c1, fp);
	    }
	    if (c == ';') {
		csv_set_got_semi(cdata);
	    }
	    if (c == cdata->delim) {
		csv_set_got_delim(cdata);
	    } else if (c == '"') {
		ldquo++;
		ndquo++;
	    } else if (c == '\'') {
		lsquo++;
		nsquo++;
	    }
	}
	cc++;
    }

    if (maxlinelen == 0) {
	pputs(prn, _("Data file is empty\n"));
    } else if (csv_has_trailing_comma(cdata)) {
	pputs(prn, _("Data file has trailing commas\n"));
    }

    if (ndquo > 0 || nsquo > 0) {
	/* candidates for quotation character? */
	int cands[2] = {0};

	if (ndquo > 0) {
	    pprintf(prn, _("Found %d double-quotes, max %d per line\n"),
		    ndquo, max_ldquo);
	}
	if (nsquo > 0) {
	    pprintf(prn, _("Found %d single-quotes, max %d per line\n"),
		    nsquo, max_lsquo);
	}
	if (max_ldquo > 0 && max_ldquo % 2 == 0) {
	    /* double-quote is a candidate? */
	    if (min_ldquo > 0 && min_ldquo % 2) {
		; /* nope */
	    } else {
		cands[0] = 1;
	    }
	}
	if (max_lsquo > 0 && max_lsquo % 2 == 0) {
	    /* single-quote is a candidate? */
	    if (min_lsquo > 0 && min_lsquo % 2) {
		; /* nope */
	    } else {
		cands[1] = 1;
	    }
	}
	if (cands[0] && cands[1]) {
	    /* hmm, rule one out: prefer the more numerous */
	    if (nsquo > ndquo) {
		cands[0] = 0;
	    } else {
		cands[1] = 0;
	    }
	}
	if (cands[0]) {
	    pputs(prn, _("Assuming double-quote is the relevant "
			 "quotation character\n"));
	    cdata->qchar = '"';
	} else if (cands[1]) {
	    pputs(prn, _("Assuming single-quote is the relevant "
			 "quotation character\n"));
	    cdata->qchar = '\'';
	}
    }

    if (maxlinelen > 0) {
	/* allow for newline and null terminator */
	maxlinelen += 2 + crlf;
    }

    return maxlinelen;
}

#define nonspace_delim(d) (d != ',' && d != ';' && d != '\t')

static int count_csv_fields (csvdata *c)
{
    const char *s = c->line;
    int inquote = 0;
    int cbak, nf = 0;

    if (*s == c->delim && *s == ' ') {
	s++;
    }

    while (*s) {
	if (csv_keep_quotes(c) && *s == c->qchar) {
	    inquote = !inquote;
	} else if (!inquote && *s == c->delim) {
	    nf++;
	}
	cbak = *s;
	s++;
	/* Problem: (when) should a trailing delimiter be read as an
	   implicit NA?  For now we'll so treat it if the delimiter
	   is not plain space.
	*/
	if (*s == '\0' && cbak == c->delim && nonspace_delim(c->delim)) {
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

static void compress_csv_line (csvdata *c, int nospace)
{
    int n = strlen(c->line);
    char *p = c->line + n - 1;

    if (*p == 0x0a) {
	*p = '\0';
	p--;
    }

    if (*p == 0x0d) {
	*p = '\0';
    }

    if (!csv_keep_quotes(c) && c->delim == ',') {
	purge_quoted_commas(c->line);
    }

    if (c->delim != ' ') {
	if (nospace) {
	    purge_unquoted_spaces(c->line);
	}
    } else {
	compress_spaces(c->line);
    }

    if (!csv_keep_quotes(c)) {
        gretl_delchar('"', c->line);
    }

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
    char tmp[VNAMELEN];

    if (s == NULL) {
	return 1;
    }

    if (!strcmp(s, "\"\"") || !strcmp(s, "''")) {
	return 1;
    }

    if (*s == '"' || *s == '\'') s++;

    if (*s == '\0') {
	return 1;
    }

    if (strlen(s) > VNAMELEN - 1) {
	return 0;
    }

    *tmp = '\0';
    strncat(tmp, s, VNAMELEN - 1);
    gretl_lower(tmp);

    return (!strcmp(tmp, "obs") ||
	    !strcmp(tmp, "date") ||
	    !strcmp(tmp, "year") ||
	    !strcmp(tmp, "period") ||
	    !strcmp(tmp, "observation") ||
	    !strcmp(tmp, "observation_date"));
}

static int join_wants_col_zero (csvdata *c, const char *s)
{
    const char *colname;
    int i;

    if (*s == '\0') {
	return 0;
    }

    for (i=0; i<c->jspec->ncols; i++) {
	colname = c->jspec->colnames[i];
	if (colname != NULL && !strcmp(s, colname)) {
	    return 1;
	}
    }

    return 0;
}

static void check_first_field (const char *line, csvdata *c, PRN *prn)
{
    const char *s;

 tryagain:
    s = line;

    if (c->delim != ' ' && *s == c->delim) {
	csv_set_blank_column(c);
    } else {
	char field1[OBSLEN];
	int i = 0;

	if (c->delim == ' ' && *s == ' ') {
	    s++;
	}

	while (*s && i < sizeof field1) {
	    if (*s == c->delim) {
		break;
	    } else if (*s == '\t') {
		/* presence of a tab must indicate tab-separation? */
		c->delim = '\t';
		goto tryagain;
	    }
	    field1[i++] = *s++;
	}

	field1[i] = '\0';
	iso_to_ascii(field1);

	if (joining(c) && join_wants_col_zero(c, field1)) {
	    return;
	} else if (csv_all_cols(c)) {
	    /* open/append wants all columns as data */
	    return;
	}

	pprintf(prn, _("   first field: '%s'\n"), field1);

	if (import_obs_label(field1)) {
	    pputs(prn, _("   seems to be observation label\n"));
	    csv_set_obs_column(c);
	}
    }
}

void import_na_init (void)
{
    const char *s = get_csv_na_read_string();

    strcpy(import_na, s);
}

/* Returns 1 if the string @s should be counted representing
   an NA or missing value, 0 otherwise. If there is a user-set
   "csv_read_na" value this is used for comparison, otherwise
   a set of default values is consulted.
*/

int import_na_string (const char *s)
{
    if (*import_na != '\0' && strcmp(import_na, "default")) {
	/* the user has set a specific "NA" string, so
	   respect it */
	return !strcmp(s, import_na);
    } else {
	/* consult a list of common representations of NA */
	const char *defaults[] = {
	    "NA",
	    "N.A.",
	    "n.a.",
	    "na",
	    "n/a",
	    "N/A",
	    "#N/A",
	    "NaN",
	    ".NaN",
	    ".",
	    "..",
	    "-999",
	    "-9999",
	    "-",
	    NULL
	};
	int i;

	for (i=0; defaults[i] != NULL; i++) {
	    if (!strcmp(s, defaults[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

static int csv_missval (const char *str, int i, int t,
			int *miss_shown, PRN *prn)
{
    int miss = 0;

    if (*str == '\0' || !strcmp(str, "\"\"")) {
	/* 2021-03-03: let '""' indicate missing */
	if (miss_shown != NULL) {
	    if (t < 80 || *miss_shown < i) {
		pprintf(prn, _("   the cell for variable %d, obs %d "
			       "is empty: treating as missing value\n"),
			i, t);
		*miss_shown += 1;
	    }
	}
	miss = 1;
    }

    if (import_na_string(str)) {
	if (miss_shown != NULL) {
	    if (t < 80 || *miss_shown < i) {
		pprintf(prn, _("   warning: missing value for variable "
			       "%d, obs %d\n"), i, t);
		*miss_shown += 1;
	    }
	}
	miss = 1;
    }

    return miss;
}

/* In the case where we think we've found thousands
   separators in numerical input, provisionally mark
   all "non-numeric" values as NAs; we do this prior
   to a second pass through the data.
*/

static void revise_non_numeric_values (csvdata *c)
{
    int i, t;

    for (i=1; i<c->dset->v; i++) {
	for (t=0; t<c->dset->n; t++) {
	    if (c->dset->Z[i][t] == NON_NUMERIC) {
		c->dset->Z[i][t] = NADBL;
	    }
	}
    }
}

int non_numeric_check (DATASET *dset, int **plist,
		       gretl_string_table **pst,
		       PRN *prn)
{
    int *list = NULL;
    int i, j, t, nn = 0;
    int err = 0;

#if CDEBUG > 1
    fprintf(stderr, "non_numeric_check: testing %d series, pst = %p\n",
	    dset->v - 1, (void *) pst);
#endif

    if (pst == NULL) {
	/* not interested in string-valued series/columns */
	for (i=1; i<dset->v; i++) {
	    for (t=0; t<dset->n; t++) {
		if (dset->Z[i][t] == NON_NUMERIC) {
		    dset->Z[i][t] = NADBL;
		}
	    }
	}
	return 0;
    }

    for (i=1; i<dset->v; i++) {
	for (t=0; t<dset->n; t++) {
	    if (dset->Z[i][t] == NON_NUMERIC) {
		nn++;
		break;
	    }
	}
    }

#if CDEBUG > 1
    fprintf(stderr, " found %d candidate series\n", nn);
#endif

    if (nn == 0) {
	return 0; /* nothing to be done */
    }

    list = gretl_list_new(nn);
    if (list == NULL) {
	return E_ALLOC;
    }

    j = 1;
    for (i=1; i<dset->v; i++) {
	for (t=0; t<dset->n; t++) {
	    if (dset->Z[i][t] == NON_NUMERIC) {
		list[j++] = i;
		break;
	    }
	}
    }

#if CDEBUG > 1
    printlist(list, "non-numeric vars list");
#endif

    for (i=1; i<=list[0]; i++) {
	/* check each member of @list */
	double nnfrac;
	int nnon = 0;
	int tnon = -1;
	int nok = 0;
	int v = list[i];

	series_set_flag(dset, v, VAR_DISCRETE);

	for (t=0; t<dset->n; t++) {
	    if (dset->Z[v][t] == NON_NUMERIC) {
		if (tnon < 0) {
		    /* record the first non-numeric obs */
		    tnon = t + 1;
		}
		nnon++;
	    } else if (!na(dset->Z[v][t])) {
		nok++;
	    }
	}

	nnfrac = (nok == 0)? 1.0 : nnon / (double) (nnon + nok);
	pprintf(prn, _("variable %d (%s): non-numeric values = %d "
		       "(%.2f percent)\n"), v, dset->varname[v],
		nnon, 100 * nnfrac);
	if ((nnon < 2 && dset->n > 2) || nnfrac < 0.05) {
	    /* if we got just a few non-numeric values, we'll assume
	       that the data file is broken
	    */
	    pprintf(prn, _("ERROR: variable %d (%s), observation %d, "
			   "expected numeric value\n"),
		    v, dset->varname[v], tnon);
	    err = E_DATA;
	    break;
	}
    }

    if (!err) {
	pputs(prn, _("allocating string table\n"));
	*pst = gretl_string_table_new(list);
	if (*pst == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	free(list);
    } else {
	*plist = list;
    }

    return err;
}

static int csv_non_numeric_check (csvdata *c, PRN *prn)
{
    gretl_string_table *st = NULL;
    int *nlist = NULL;
    int err = 0;

    if (csv_as_matrix(c)) {
	err = non_numeric_check(c->dset, &nlist, NULL, prn);
    } else {
	err = non_numeric_check(c->dset, &nlist, &st, prn);
    }

    if (!err) {
	c->codelist = nlist;
	c->st = st;
    }

    return err;
}

/* Handle the case in "join" where the user specified some time
   columns for conversion to numeric and also gave a specific format
   for the conversion.
*/

static double special_time_val (const char *s, const char *fmt,
				int m_means_q)
{
    struct tm t = {0};
    char *test;

    test = strptime(s, fmt, &t);

    if (test == NULL || *test != '\0') {
	/* conversion didn't work right */
	return NADBL;
    } else {
	int y, m, d;

	y = t.tm_year + 1900;
	m = t.tm_mon + 1;
	d = t.tm_mday;

	if (m_means_q) {
	    /* convert to 1st month of quarter */
	    if (m == 2) m = 4;
	    else if (m == 3) m = 7;
	    else if (m == 4) m = 10;
	    else if (m != 1) {
		return NADBL;
	    }
	}

	if (d == 0) d = 1;

	return 10000*y + 100*m + d;
    }
}

static int char_count (char c, const char *s)
{
    int n = 0;

    while (*s) {
	if (*s == c) n++;
	s++;
    }

    return n;
}

/* Follow-up check for the case where we think we might
   have found a thousands separator: each occurrence of
   the putative separator must be followed by exactly 3
   digits: we set c->thousep to an invalid value if this
   is not the case.
*/

static void validate_thousep (csvdata *c, const char *s)
{
    int nd;

    while (*s) {
	if (*s == c->thousep) {
	    nd = 0;
	    s++;
	    while (*s) {
		if (isdigit(*s)) {
		    nd++;
		    s++;
		} else {
		    break;
		}
	    }
	    if (nd != 3) {
		/* nope! */
#if CDEBUG
		fprintf(stderr, "validate_thousep: no: '%c' is followed by %d digits\n",
			c->thousep, nd);
#endif
		c->thousep = -1;
		break;
	    }
	} else {
	    s++;
	}
    }
}

/* Initial heuristic for detecting a thousands separator,
   where the string @s has been determined to contain
   nothing but digits, dot and comma (allowing for a leading
   minus).

   1) If the string contains both comma and dot, whichever
   character appears to the left cannot be the decimal
   separator and may be a thousands separator.

   2) If more than one comma appears in the string, comma
   cannot be the decimal character and might be a thousands
   separator; mutatis mutandis for dot.
*/

static void test_for_thousands_sep (csvdata *c, const char *s)
{
    const char *p1 = strrchr(s, '.');
    const char *p2 = strrchr(s, ',');
    char thousep = 0;

    if (p1 != NULL && p2 != NULL) {
	thousep = (p2 - p1 > 0)? '.' : ',';
    } else if (p1 != NULL && char_count('.', s) > 0) {
	thousep = '.';
    } else if (p2 != NULL && char_count(',', s) > 0) {
	thousep = ',';
    }

    if (c->thousep > 0) {
	if (thousep != 0 && thousep != c->thousep) {
	    /* no consistent interpretation exists */
	    c->thousep = -1; /* invalid */
	}
    } else if (thousep != 0) {
	/* we have a candidate for testing */
	char *test, tmp[CSVSTRLEN];

	strcpy(tmp, s);
	gretl_delchar(thousep, tmp);
	if (thousep == '.' && get_local_decpoint() == '.') {
	    gretl_charsub(tmp, ',', '.');
	}
	errno = 0;
	strtod(tmp, &test);
	if (*test == '\0' && errno == 0) {
	    c->thousep = thousep;
	}
    }

    if (c->thousep && thousep != 0) {
	validate_thousep(c, s);
    }
}

static int all_digits_and_seps (const char *s)
{
    const char *test = "0123456789.,";

    if (*s == '-') s++;

    return strspn(s, test) == strlen(s);
}

static double eval_non_numeric (csvdata *c, int i, const char *s)
{
    double x = NON_NUMERIC;

    if (series_get_flags(c->dset, i) & VAR_TIMECOL) {
	char *fmt = NULL;
	int mq = 0;

	if (timecol_get_format(c->dset, i, &fmt, &mq)) {
	    /* the user gave a specific format for this */
	    x = special_time_val(s, fmt, mq);
	} else {
	    /* default: ISO 8601 extended */
	    int y, m, d, n;

	    n = sscanf(s, "%d-%d-%d", &y, &m, &d);
	    if (n == 3) {
		x = 10000*y + 100*m + d;
	    } else {
		x = NADBL;
	    }
	}
    } else if (c->thousep >= 0 && !csv_scrub_thousep(c)) {
	/* Here we consider the possibility although @s does not
	   validate as numeric according to the C library, it is by
	   intent numeric but includes one or more thousands
	   separators.

	   The condition c->thousep >= 0 requires that we haven't
	   already ruled out this interpretation due to inconsistency,
	   and !csv_scrub_thousep(c) requires that we're not on a
	   second pass through the data.
	*/
	if (all_digits_and_seps(s)) {
	    test_for_thousands_sep(c, s);
	}
    }

    return x;
}

static int converted_ok (const char *s, char *test, double x)
{
    if (*test != '\0') {
	if (errno) perror(s);
	return 0; /* definitely not OK */
    } else if (errno == ERANGE && fabs(x) > 0 && fabs(x) < 0.001) {
	return 1; /* subnormal, but we'll let that pass */
    } else if (errno) {
	perror(s);
	return 0;
    } else {
	return 1;
    }
}

static char *csv_unquote (char *s)
{
    if (s[0] == '"') {
	int i, n = strlen(s);

	if (n > 1 && s[n-1] == '"') {
	    for (i=0; i<n-2; i++) {
		s[i] = s[i+1];
	    }
	    s[i] = '\0';
	}
    }
    return s;
}

static double csv_atof (csvdata *c, int i)
{
    char tmp[CSVSTRLEN], clean[CSVSTRLEN];
    double x = NON_NUMERIC;
    const char *s = c->str;
    char *test;

    if (csv_scrub_thousep(c) && strchr(s, c->thousep) &&
	all_digits_and_seps(s)) {
	/* second pass through the data: pre-process fields
	   that we reckon include thousands separators
	*/
	strcpy(clean, s);
	gretl_delchar(c->thousep, clean);
	s = clean;
    }

    if (c->decpoint == '.' || !csv_do_dotsub(c) || strchr(s, ',') == NULL) {
	/* either we're currently set to the correct locale,
	   or there's no problematic decimal point in @s
	*/
	errno = 0;
	x = strtod(s, &test);
	if (converted_ok(s, test, x)) {
	    return x; /* handled */
	}
    } else if (csv_do_dotsub(c)) {
	/* in C numeric locale: substitute dot for comma */
	strcpy(tmp, s);
	gretl_charsub(tmp, ',', '.');
	errno = 0;
	x = strtod(tmp, &test);
	if (converted_ok(s, test, x)) {
	    return x; /* handled */
	}
    }

    if (c->decpoint == '.' && strchr(s, ',') != NULL) {
	/* try remediation for decimal comma? */
	strcpy(tmp, s);
	gretl_charsub(tmp, ',', '.');
	errno = 0;
	x = strtod(tmp, &test);
	if (converted_ok(s, test, x)) {
	    return x; /* handled */
	}
    }

    /* fallback */
    /* revised 2020-02-13 to use csv_unquote */
    return eval_non_numeric(c, i, csv_unquote(c->str));
}

static int process_csv_obs (csvdata *c, int i, int t, int *miss_shown,
			    PRN *prn)
{
    int err = 0;

    if (c->st != NULL) {
	/* second round, handling string-valued variables */
	if (in_gretl_list(c->codelist, i)) {
	    double zit = c->dset->Z[i][t];
	    int ix;

	    if (na(zit) && *c->str != '\0' && c->user_na == NULL) {
		/* by default (no user_na) only blanks count as NAs */
		zit = NON_NUMERIC;
	    }
	    if (!na(zit)) {
		ix = gretl_string_table_index(c->st, c->str, i, 0, prn);
		if (ix > 0) {
		    c->dset->Z[i][t] = (double) ix;
		} else {
		    err = E_DATA;
		}
	    }
	}
    } else if (csv_missval(c->str, i, t+1, miss_shown, prn)) {
	c->dset->Z[i][t] = NADBL;
    } else {
	gretl_strstrip(c->str);
	c->dset->Z[i][t] = csv_atof(c, i);
    }

    return err;
}

/* Emulation of fgets(), designed to handle any sort of line
   termination (unix, DOS, Mac or even an unholy mixture).
   Line-endings are converted to LF (0x0a).
*/

static char *csv_fgets (csvdata *cdata, gzFile fp)
{
    char *s = cdata->line;
    int n = cdata->maxlinelen;
    int i, c1, c = 0;

    for (i=0; i<n-1 && c!=0x0a; i++) {
	c = gzgetc(fp);
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
	    c1 = gzgetc(fp);
	    if (c1 != 0x0a) {
		gzungetc(c1, fp);
	    }
	}
	s[i] = c;
    }

    s[i] = '\0';

    return s;
}

/* pick up any comments following the data block in a CSV file */

static char *get_csv_descrip (csvdata *c, gzFile fp)
{
    char *line = c->line;
    char *desc = NULL;
    size_t llen, totlen;

    while (csv_fgets(c, fp)) {
	tailstrip(line);
	llen = strlen(line);
	if (desc == NULL) {
	    totlen = llen + 4;
	    desc = malloc(totlen);
	    if (desc == NULL) {
		return NULL;
	    }
	    sprintf(desc, "%s\n", line);
	} else {
	    char *tmp;

	    totlen = strlen(desc) + llen + 4;
	    tmp = realloc(desc, totlen);
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

/* Here we check whether we get a consistent reading on
   the number of fields per line in the CSV file
*/

static int csv_fields_check (gzFile fp, csvdata *c, PRN *prn)
{
    int gotdata = 0;
    int chkcols = 0;
    int err = 0;

    c->ncols = c->nrows = 0;

    if (csv_has_bom(c)) {
	gzseek(fp, 3, SEEK_SET);
    }

    while (csv_fgets(c, fp) && !err) {

	/* skip comment lines */
	if (*c->line == '#') {
	    continue;
	}

	/* skip blank lines -- but finish if the blank comes after data */
	if (string_is_blank(c->line)) {
	    if (gotdata) {
		if (!csv_have_data(c)) {
		    c->descrip = get_csv_descrip(c, fp);
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

	compress_csv_line(c, 1);

	if (!gotdata) {
	    /* scrutinize the first "real" line */
	    check_first_field(c->line, c, prn);
	    gotdata = 1;
	}

	chkcols = count_csv_fields(c);
	if (c->ncols == 0) {
	    c->ncols = chkcols;
	    pprintf(prn, _("   number of columns = %d\n"), c->ncols);
	} else if (chkcols != c->ncols) {
	    pprintf(prn, _("   ...but row %d has %d fields: aborting\n"),
		    c->nrows, chkcols);
	    err = E_DATA;
	} else if (cols_subset(c)) {
	    int datacols = csv_skip_col_1(c) ? (c->ncols - 1) : c->ncols;

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
    char name[VNAMELEN] = {0};
    int i, j = 0;

    for (i=0; s[i] != '\0'; i++) {
	if (isalnum(s[i]) || s[i] == '_') {
	    name[j++] = s[i];
	}
    }

    name[j] = '\0';
    strcpy(s, name);
}

static int intercept_nan_as_name (const char *s)
{
    if (strlen(s) == 3) {
	char screen[4];

	strcpy(screen, s);
	gretl_lower(screen);
	if (!strcmp(screen, "nan")) {
	    return 1;
	}
    }

    return 0;
}

static int csv_is_numeric (const char *s, csvdata *c)
{
    int ret = 0;

    if (c->decpoint == '.') {
	ret = numeric_string(s);
    } else {
	/* decimal comma in force */
	char *tmp = gretl_strdup(s);

	gretl_charsub(tmp, ',', '.');
	ret = numeric_string(tmp);
	free(tmp);
    }

    return ret;
}

static int process_csv_varname (csvdata *c, int j, int *numcount,
				PRN *prn)
{
    char *vname = c->dset->varname[j];
    char *src = c->str;
    int err = 0;

    *vname = '\0';

    if (intercept_nan_as_name(src)) {
	gretl_errmsg_sprintf(_("If '%s' is intended as the name of a variable, "
			       "please change it --\nstrings of this sort usually "
			       "mean 'not a number'."), src);
	err = E_DATA;
    } else if (*src == '\0') {
	fprintf(stderr, "variable name %d is missing\n", j);
	sprintf(vname, "v%d", j);
    } else if (csv_is_numeric(src, c)) {
	*numcount += 1;
    } else {
	const char *s = src;

	while (*s && !isalpha(*s)) s++;
	if (*s == '\0') {
	    fprintf(stderr, "variable name %d (%s) is garbage\n", j, src);
	    sprintf(vname, "v%d", j);
	} else {
	    strncat(vname, s, VNAMELEN - 1);
	}
	iso_to_ascii(vname);
	strip_illegals(vname);
	if (check_varname(vname)) {
	    errmsg(1, prn);
	    err = E_DATA;
	}
    }

    return err;
}

static int csv_reconfigure_for_markers (DATASET *dset)
{
    int err = dataset_allocate_obs_markers(dset);

    if (!err) {
	err = dataset_drop_last_variables(dset, 1);
    }

    return err;
}

static int skip_data_column (csvdata *c, int k)
{
    int col = csv_skip_col_1(c) ? k : k + 1;

    if (!in_gretl_list(c->cols_list, col)) {
	return 1;
    } else {
	return 0;
    }
}

static int update_join_cols_list (csvdata *c, int k)
{
    int *test;
    int err = 0;

    test = gretl_list_append_term(&c->cols_list, k);
    if (test == NULL) {
	err = E_ALLOC;
    }

#if CDEBUG
    printlist(c->cols_list, "c->cols_list for join");
#endif

    return err;
}

/* handle_join_varname: the index @k contains the column number
   relative to the entire CSV file, while @pj points to j, the column
   number relative to the reduced dataset that will be constructed by
   selection of columns from the file.

   Here we're examining a column heading read from file (c->str) to
   see whether it matches any of the column-names required for an
   ongoing join operation (held in c->jspec->colnames). If so, we
   write the index j into the appropriate slot in c->jspec->colnums
   (which starts off filled with zeros), so the joiner will know where
   to find the required data. (The j value is bound to be at least 1
   since column 0 is reserved to the constant.)

   In some cases a given named column may perform more than one role in
   a join operation -- for example, it may serve as an element in a
   filter and also as the auxiliary variable in an "aggregation"
   method. To allow for this we don't stop scanning at the first match
   of c->str with a required column name.

   The call to update_join_cols_list() uses the index @k to record the
   overall column position of "wanted data", for use by the CSV
   reader.
*/

static int handle_join_varname (csvdata *c, int k, int *pj)
{
    const char *colname;
    char okname[VNAMELEN];
    int matched = 0;
    int i, j = *pj;

    if (!csv_skip_col_1(c)) {
	k++;
    }

    if (csv_no_header(c)) {
	sprintf(okname, "col%d", k);
    } else {
	/* convert to valid gretl identifier */
	gretl_normalize_varname(okname, c->str, 0, k);
    }

#if CDEBUG
    fprintf(stderr, "handle_join_varname: looking at '%s' (%s)\n", c->str, okname);
#endif

    for (i=0; i<c->jspec->ncols; i++) {
	/* find "wanted name" i */
	colname = c->jspec->colnames[i];
	if (colname == NULL || c->jspec->colnums[i] > 0) {
	    /* name not wanted, or already found */
	    continue;
	}
	if (!strcmp(okname, colname)) {
#if CDEBUG
	    fprintf(stderr, " target %d matched at CSV col %d, j=%d\n", i, k, j);
#endif
	    c->jspec->colnums[i] = j;
	    if (!matched) {
		matched = 1;
		strcpy(c->dset->varname[j], okname);
		update_join_cols_list(c, k);
		*pj += 1;
		if (in_gretl_list(c->jspec->timecols, i)) {
		    series_set_flag(c->dset, j, VAR_TIMECOL);
		}
	    }
	}
    }

    return 0;
}

#define starts_number(c) (isdigit((unsigned char) c) || c == '-' ||	\
                          c == '+' || c == '.')

#define obs_labels_no_varnames(o,c,n)  (!o && c->v > 3 && n == c->v - 2)

static int csv_varname_scan (csvdata *c, gzFile fp, PRN *prn, PRN *mprn)
{
    char *p;
    int obscol = csv_has_obs_column(c);
    int i, j, k, numcount;
    int err = 0;

    if (!csv_no_header(c)) {
	pputs(mprn, _("scanning for variable names...\n"));
    }

    if (csv_has_bom(c)) {
	gzseek(fp, 3, SEEK_SET);
    }

    while (csv_fgets(c, fp)) {
	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	} else {
	    break;
	}
    }

    c->datapos = gztell(fp);

    compress_csv_line(c, 1);

    p = c->line;
    if (c->delim == ' ' && *p == ' ') p++;
    iso_to_ascii(p);

    if (strlen(p) > 118) {
	pprintf(mprn, _("   line: %.115s...\n"), p);
    } else {
	pprintf(mprn, _("   line: %s\n"), p);
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

	if (k == 0 && csv_skip_col_1(c)) {
	    ; /* no-op */
	} else if (!joining(c) && cols_subset(c) && skip_data_column(c, k)) {
	    ; /* no-op */
	} else {
	    if (joining(c)) {
		handle_join_varname(c, k, &j);
	    } else if (probing(c) && csv_no_header(c)) {
		sprintf(c->dset->varname[j], "col%d", j);
		j++;
	    } else {
		err = process_csv_varname(c, j, &numcount, prn);
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
	gretl_errmsg_set("No relevant columns were found");
	err = E_UNKVAR;
    }

    if (err) {
	return err;
    }

    if (csv_no_header(c) || numcount == c->dset->v - 1 ||
	obs_labels_no_varnames(obscol, c->dset, numcount)) {
	if (!csv_no_header(c)) {
	    pputs(prn, _("it seems there are no variable names\n"));
	    /* then we undercounted the observations by one? */
	    if (!rows_subset(c)) {
		err = add_single_obs(c->dset);
	    }
	}
	if (!err) {
	    /* set up to handle the "no varnames" case */
	    csv_set_autoname(c);
	    c->datapos = csv_has_bom(c) ? 3 : 0;
	    if (!csv_all_cols(c)) {
		if (obs_labels_no_varnames(obscol, c->dset, numcount)) {
		    err = csv_reconfigure_for_markers(c->dset);
		    if (!err) {
			csv_set_obs_column(c);
		    }
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
	fprintf(stderr, "numcount = %d\n", numcount);
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

static int fixed_format_read (csvdata *c, gzFile fp, PRN *prn)
{
    char *p;
    int miss_shown = 0;
    int *missp = NULL;
    int t = 0, s = 0;
    int i, k, n, m;
    int err = 0;

    c->real_n = c->dset->n;

    if (csv_has_bom(c)) {
	gzseek(fp, 3, SEEK_SET);
    }

    if (csv_is_verbose(c)) {
	missp = &miss_shown;
    }

    while (csv_fgets(c, fp) && !err) {
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
	    /* Added 2016-11-16: allow trailing blanks in a field
	       of specified width. This is required for handling
	       US CPS data.
	    */
	    tailstrip(c->str);
	    if (csv_missval(c->str, i, t+1, missp, prn)) {
		c->dset->Z[i][t] = NADBL;
	    } else {
		c->dset->Z[i][t] = csv_atof(c, i);
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

/* Check that an observation label contains only
   valid UTF-8, and moreover that every character
   is valid in XML 1.0. If not, try recoding from
   ISO 8859.
*/

static int maybe_fix_csv_string (gchar *s)
{
    int err = 0;

    if (!g_utf8_validate(s, -1, NULL)) {
	GError *gerr = NULL;
	gsize wrote = 0;
	gchar *tr;

	/* try for iso-8859? */
	tr = g_convert(s, -1, "UTF-8", "ISO-8859-15",
		       NULL, &wrote, &gerr);
	if (gerr != NULL) {
	    gretl_errmsg_set(gerr->message);
	    g_error_free(gerr);
	    err = E_DATA;
	} else {
	    *s = '\0';
	    gretl_utf8_strncat(s, tr, CSVSTRLEN-1);
	    g_free(tr);
	}
    }

    if (!err) {
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

    return err;
}

static void transcribe_obs_label (csvdata *c, int t)
{
    char *s = c->str;
    char c0 = *s;
    int n = strlen(s);

    /* skip a leading quote, and unquote fully
       if a matching trailing quote is found
    */

    if (c0 == '"' || c0 == '\'') {
	if (s[n-1] == c0) {
	    s[n-1] = '\0';
	    n--;
	}
	s++;
	n--;
	/* and once more, with feeling... */
	if (s[0] == '\'') {
	    s++;
	    n--;
	}
    }

    if (n > OBSLEN - 1) {
	n = OBSLEN - 1;
    }

    c->dset->S[t][0] = '\0';
    gretl_utf8_strncat(c->dset->S[t], s, n);
}

static int real_read_labels_and_data (csvdata *c, gzFile fp, PRN *prn)
{
    char *p;
    int miss_shown = 0;
    int *missp = NULL;
    int truncated = 0;
    int t = 0, s = 0;
    int i, j, k;
    int err = 0;

    if (csv_is_verbose(c)) {
	missp = &miss_shown;
    }

    c->real_n = c->dset->n;

    while (csv_fgets(c, fp) && !err) {
	int inquote = 0;

	if (*c->line == '#' || string_is_blank(c->line)) {
	    continue;
	} else if (*c->skipstr != '\0' && strstr(c->line, c->skipstr)) {
	    c->real_n -= 1;
	    continue;
	} else if (row_not_wanted(c, s)) {
	    s++;
	    continue;
	}

	compress_csv_line(c, 0);
	p = c->line;

	if (c->delim == ' ') {
	    if (*p == ' ') p++;
	} else {
	    p += strspn(p, " ");
	}

	j = 1;
	for (k=0; k<c->ncols && !err; k++) {
	    i = 0;
	    while (*p) {
		if (csv_keep_quotes(c) && *p == c->qchar) {
		    inquote = !inquote;
		} else if (!inquote && *p == c->delim) {
		    break;
		}
		if (i < CSVSTRLEN - 1) {
		    c->str[i++] = *p;
		} else {
		    truncated++;
		}
		p++;
	    }
	    c->str[i] = '\0';
	    err = maybe_fix_csv_string(c->str);
	    if (!err) {
		if (k == 0 && csv_skip_col_1(c) && c->dset->S != NULL) {
		    transcribe_obs_label(c, t);
		} else if (cols_subset(c) && skip_data_column(c, k)) {
		    ; /* no-op */
		} else {
		    err = process_csv_obs(c, j++, t, missp, prn);
		}
	    }
	    if (!err) {
		/* prep for next column */
		if (*p == c->delim) {
		    p++;
		}
		if (c->delim != ' ') {
		    p += strspn(p, " ");
		}
	    }
	}

	s++;
	if (++t == c->dset->n) {
	    break;
	}
    }

    if (truncated) {
	pprintf(prn, _("warning: %d labels were truncated.\n"), truncated);
    }

    if (!err && c->real_n < c->dset->n) {
	int drop = c->dset->n - c->real_n;

	err = dataset_drop_observations(c->dset, drop);
    }

    return err;
}

/* When reading a CSV file, should we attempt to parse observation
   strings as dates (and impose time-series structure on the data
   if this is successful)? In general, yes, but maybe not if we're
   reading the data in the context of a "join" operation, since
   in this case automatic detection may collide with time-key
   information supplied by the user. Current status: we'll skip
   the auto-dating stuff when joining unless (a) it's a MIDAS
   join (mixed frequencies) and the user has _not_ supplied any
   time key specification.
*/

static int csv_skip_dates (csvdata *c)
{
    if (c->jspec != NULL) {
	/* with --aggr=spread (MIDAS) we'll need dates info,
	   unless the user have a time key spec
	*/
	return c->jspec->auto_midas == 0;
    } else {
	return 0;
    }
}

static int csv_read_data (csvdata *c, gzFile fp, PRN *prn, PRN *mprn)
{
    int reversed = csv_data_reversed(c);
    int err;

    if (mprn != NULL) {
	if (csv_all_cols(c)) {
	    pputs(mprn, _("scanning for data...\n"));
	} else {
	    pputs(mprn, _("scanning for row labels and data...\n"));
	}
    }

    gzseek(fp, c->datapos, SEEK_SET);

    err = real_read_labels_and_data(c, fp, prn);

    if (!err && csv_skip_col_1(c) && !rows_subset(c) && !csv_skip_dates(c)) {
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

	pprintf(prn, "%s %s...\n", _("parsing"), trfname);
	g_free(trfname);
    } else {
	pprintf(prn, "%s %s...\n", _("parsing"), fname);
    }
}

static int join_unique_columns (csvdata *c)
{
    const char **cnames = c->jspec->colnames;
    char *counted;
    int i, j, ncols = 0;

    counted = calloc(c->jspec->ncols, 1);

    for (i=0; i<c->jspec->ncols; i++) {
	if (cnames[i] != NULL && counted[i] == 0) {
	    counted[i] = 1;
	    /* mark any duplicates as counted too */
	    for (j=i+1; j<c->jspec->ncols; j++) {
		if (cnames[j] != NULL && !strcmp(cnames[i], cnames[j])) {
		    counted[j] = 1;
		}
	    }
#if CDEBUG
	    fprintf(stderr, "join_unique_columns: '%s'\n", cnames[i]);
#endif
	    ncols++;
	}
    }

    free(counted);

    return ncols;
}

static int csv_set_dataset_dimensions (csvdata *c)
{
    int err = 0;

    c->dset->v = 0;

    if (rows_subset(c)) {
	c->dset->n = n_from_row_mask(c);
    }

    if (fixed_format(c)) {
	if (c->dset->n == 0) {
	    c->dset->n = c->nrows;
	}
	c->dset->v = c->ncols + 1;
    } else {
	int cols_wanted, cols_present;

	if (c->dset->n == 0) {
	    if (csv_no_header(c)) {
		c->dset->n = c->nrows;
	    } else {
		/* allow for varnames row */
		c->dset->n = c->nrows - 1;
	    }
	}

	cols_present = csv_skip_col_1(c) ? (c->ncols - 1) : c->ncols;

	if (joining(c)) {
	    cols_wanted = join_unique_columns(c);
	} else if (cols_subset(c)) {
	    cols_wanted = c->cols_list[0];
	} else {
	    cols_wanted = cols_present;
	}

	if (cols_wanted > cols_present) {
	    gretl_errmsg_set(_("Invalid column specification"));
	    err = E_DATA;
	} else {
	    /* allow for the constant */
	    c->dset->v = cols_wanted + 1;
	}
    }

    if (probing(c)) {
	/* don't allocate tons of space for data that
	   we won't read right now */
	c->dset->n = 1;
    }

#if CDEBUG
    if (joining(c)) {
	fprintf(stderr, "csv dataset dimensions: v=%d, n=%d\n",
		c->dset->v, c->dset->n);
    }
#endif

    return err;
}

/*
 * real_import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @cols: column specification.
 * @rows: row specification.
 * @join: specification pertaining to "join" command.
 * @probe: also pertains to "join" (via GUI).
 * @pm: location of matrix to accept the data or NULL.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; use OPT_H
 * to indicate absence of a header row; use OPT_A to indicate that
 * all columns should be read as data series (i.e. do not try to
 * interpret the first column as observation labels); for use of
 * OPT_T see the help text for the "append" command.
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
			    joinspec *join,
			    csvprobe *probe,
			    gretl_matrix **pm,
			    gretlopt opt,
			    PRN *prn)
{
    csvdata *c = NULL;
    gzFile fp = NULL;
    PRN *mprn = NULL;
    gchar *altname = NULL;
    int recode = 0;
    int popit = 0;
    int i, err = 0;

    import_na_init();

    if (gretl_messages_on()) {
	mprn = prn;
    }

    fp = gretl_gzopen(fname, "rb");
    if (fp == NULL) {
	pprintf(prn, _("Couldn't open %s\n"), fname);
	err = E_FOPEN;
	goto csv_bailout;
    }

    c = csvdata_new(dset);
    if (c == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }

    recode = csv_unicode_check(fp, c, prn);
    if (recode) {
	err = csv_recode_input(&fp, fname, &altname, recode, prn);
	if (err) {
	    goto csv_bailout;
	}
    }

    if (cols != NULL) {
	err = csvdata_add_cols_list(c, cols, opt);
	if (err) {
	    goto csv_bailout;
	} else if (fixed_format(c)) {
	    pprintf(mprn, _("using fixed column format\n"));
	}
    }

    if (rows != NULL) {
	err = csvdata_add_row_mask(c, rows);
	if (err) {
	    goto csv_bailout;
	}
    }

    if (opt & OPT_H) {
	csv_set_no_header(c);
    }

    if (join != NULL) {
	c->jspec = join;
	c->flags |= CSV_HAVEDATA;
    } else if (probe != NULL) {
	c->probe = probe;
	c->flags |= CSV_HAVEDATA;
    } else {
	if (pm != NULL) {
	    csv_set_as_matrix(c);
	}
        if (opt & OPT_A) {
	    csv_set_all_cols(c);
        }
        if (opt & OPT_V) {
            csv_set_verbose(c);
        }
    }

    if (opt & OPT_I) {
	csv_unset_keep_quotes(c);
    }

    if (mprn != NULL) {
	print_csv_parsing_header(fname, mprn);
    }

    /* get line length, also check for binary data, etc. */
    c->maxlinelen = csv_max_line_length(fp, c, prn);
    if (c->maxlinelen <= 0) {
	err = E_DATA;
	goto csv_bailout;
    }

    if (csv_as_matrix(c) && csv_got_semi(c)) {
	if (c->delim == ',' && csv_got_delim(c)) {
	    c->decpoint = ',';
	}
	c->delim = ';';
    } else if (!fixed_format(c) && !csv_got_delim(c)) {
	/* set default delimiter */
	if (csv_got_tab(c)) {
	    c->delim = '\t';
	} else if (csv_got_semi(c)) {
	    c->delim = ';';
	} else {
	    c->delim = ' ';
	}
    }

#if CDEBUG
    fprintf(stderr, "fixed_format? %s; got_delim (%c)? %s; got_tab? %s; ",
	    fixed_format(c) ? "yes" : "no", c->delim,
	    csv_got_delim(c) ? "yes" : "no",
	    csv_got_tab(c)? "yes" : "no");
    fprintf(stderr, "decpoint '%c'\n", c->decpoint);
#endif

    /* buffer to hold lines */
    c->line = malloc(c->maxlinelen);
    if (c->line == NULL) {
	err = E_ALLOC;
	goto csv_bailout;
    }

 alt_delim:

    if (mprn != NULL) {
	if (!fixed_format(c)) {
	    pprintf(mprn, _("using delimiter '%c'\n"), c->delim);
	}
	pprintf(mprn, _("   longest line: %d characters\n"), c->maxlinelen - 1);
    }

    if (csv_has_trailing_comma(c) && c->delim != ',') {
	csv_unset_trailing_comma(c);
    }

    gzrewind(fp);

    /* read lines, check for consistency in number of fields */
    err = csv_fields_check(fp, c, mprn);
    if (err && !fixed_format(c)) {
	if (c->delim != ';' && csv_got_semi(c)) {
	    c->delim = ';';
	    err = 0;
	    goto alt_delim;
	}
	pputs(prn, _(csv_msg));
	goto csv_bailout;
    }

    err = csv_set_dataset_dimensions(c);
    if (err) {
	err = E_DATA;
	goto csv_bailout;
    }

    pprintf(mprn, _("   number of variables: %d\n"), c->dset->v - 1);
    pprintf(mprn, _("   number of non-blank lines: %d\n"), c->nrows);

    if (c->dset->n == 0) {
	pputs(prn, _("Invalid data file\n"));
	err = E_DATA;
	goto csv_bailout;
    }

    /* initialize CSV dataset */
    err = start_new_Z(c->dset, 0);
    if (!err && csv_skip_col_1(c)) {
	err = dataset_allocate_obs_markers(c->dset);
    }

    if (err) {
	goto csv_bailout;
    }

    /* second pass */

    gzrewind(fp);

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
    if (err || probing(c)) {
	goto csv_bailout;
    }

    if (c->decpoint == '.' && get_local_decpoint() == ',') {
	/* we're in a locale that uses decimal comma:
	   switch to the C locale */
	gretl_push_c_numeric_locale();
	popit = 1;
    } else if (c->decpoint == ',' && get_local_decpoint() == '.') {
	/* dotsub: define this if we're in a '.' locale and
	   we've figured that the decimal character is ',' in
	   the file we're reading
	*/
	csv_set_dotsub(c);
    }

    err = csv_read_data(c, fp, prn, mprn);

    if (!err) {
	/* try again, under certain conditions */
	if (csv_skip_bad(c)) {
	    err = csv_read_data(c, fp, prn, NULL);
	} else if (c->thousep > 0) {
	    pprintf(mprn, _("WARNING: it seems '%c' is being used "
			    "as thousands separator\n"), c->thousep);
	    c->decpoint = (c->thousep == '.')? ',' : '.';
	    if (c->decpoint == ',') {
		if (get_local_decpoint() == '.') {
		    csv_set_dotsub(c);
		} else if (popit) {
		    gretl_pop_c_numeric_locale();
		    popit = 0;
		}
	    }
	    revise_non_numeric_values(c);
	    csv_set_scrub_thousep(c);
	    err = csv_read_data(c, fp, prn, NULL);
	}
    }

    if (!err && !probing(c)) {
	err = csv_non_numeric_check(c, prn);
	if (!err && csv_has_non_numeric(c)) {
	    /* try once more */
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
	pputs(mprn, _("taking date information from row labels\n\n"));
	if (csv_skip_bad(c)) {
	    pprintf(prn, "WARNING: Check your data! gretl has stripped out "
		    "what appear to be\nextraneous lines in a %s dataset: "
		    "this may not be right.\n\n",
		    (c->dset->pd == 4)? "quarterly" : "monthly");
	}
    } else {
	pputs(mprn, _("treating these as undated data\n\n"));
	dataset_obs_info_default(c->dset);
    }

    if (c->dset->pd != 1 || strcmp(c->dset->stobs, "1")) {
        c->dset->structure = TIME_SERIES;
    }

    if (c->st != NULL) {
	err = gretl_string_table_validate(c->st, OPT_NONE);
	if (err) {
	    pputs(prn, _("Failed to interpret the data as numeric\n"));
	    goto csv_bailout;
	} else if (joining(c)) {
	    gretl_string_table_save(c->st, c->dset);
	} else {
	    gretl_string_table_print(c->st, c->dset, fname, prn);
	}
    }

    if (csv_as_matrix(c)) {
	/* FIXME placement of this */
	if (csv_autoname(c)) {
	    strings_array_free(c->dset->varname, c->dset->v);
	    c->dset->varname = NULL;
	}
	*pm = gretl_matrix_data_subset(NULL, c->dset, -1, -1,
				       M_MISSING_OK, &err);
	goto csv_bailout;
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
    } else {
#if CDEBUG
	int ii;

	for (ii=0; ii<c->dset->v; ii++) {
	    fprintf(stderr, " c->dset->varname[%d] = '%s'\n", ii, c->dset->varname[ii]);
	}
#endif
	if (fix_varname_duplicates(c->dset)) {
	    pputs(prn, _("warning: some variable names were duplicated\n"));
	}
    }

    if (!joining(c) && !probing(c)) {
	int newdata = (dset->Z == NULL);

	/* not doing a special "join" operation */
	err = merge_or_replace_data(dset, &c->dset, get_merge_opts(opt), prn);

	if (!err && newdata && c->descrip != NULL) {
	    dset->descrip = c->descrip;
	    c->descrip = NULL;
	}

	if (!err && newdata) {
	    dataset_add_import_info(dset, fname, GRETL_CSV);
	}
    }

 csv_bailout:

    if (fp != NULL) {
	gzclose(fp);
    }

    if (!err && c->jspec != NULL) {
	c->jspec->c = c;
    } else if (!err && c->probe != NULL) {
	c->probe->dset = c->dset;
	c->dset = NULL;
	csvdata_free(c);
    } else {
	csvdata_free(c);
    }

    if (altname != NULL) {
	gretl_remove(altname);
	g_free(altname);
    }

    if (err == E_ALLOC) {
	pputs(prn, _("Out of memory\n"));
    }

    return err;
}

/**
 * import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; for use of
 * OPT_T see the help for "append".
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
			   NULL, NULL, NULL, opt, prn);
}

gretl_matrix *import_csv_as_matrix (const char *fname, int *err)
{
#if CDEBUG
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
#else
    PRN *prn = NULL;
#endif
    gretl_matrix *m = NULL;
    char csvname[MAXLEN] = {0};
    gretlopt opt = OPT_A; /* --all-cols */
    int http = 0;

    *err = try_http(fname, csvname, &http);

    if (!*err && http) {
	*err = real_import_csv(csvname, NULL, NULL, NULL,
			       NULL, NULL, &m, opt, prn);
    } else if (!*err) {
	char fullname[FILENAME_MAX];

	strcpy(fullname, fname);
	gretl_maybe_prepend_dir(fullname);
	*err = real_import_csv(fullname, NULL, NULL, NULL,
			       NULL, NULL, &m, opt, prn);
    }

    gretl_print_destroy(prn);

    return m;
}

static int probe_varnames_check (DATASET *dset, gretlopt opt,
				 int *rerun)
{
    int missnames = 0;
    int i, err = 0;

    for (i=1; i<dset->v; i++) {
	if (dset->varname[i][0] == '\0') {
	    missnames = 1;
	    break;
	}
    }

    if (missnames) {
	if (opt & OPT_H) {
	    gretl_errmsg_set("Couldn't find all variable names");
	    err = E_DATA;
	} else {
	    *rerun = 1;
	}
    }

    return err;
}

/**
 * probe_csv:
 * @fname: name of CSV file.
 * @varnames: location to receive variable names.
 * @nvars: location to receive number of variables (columns).
 * @opt: on input, may contain any extra options to pass to
 * real_import_csv(); on return, OPT_H (indicating that the
 * CSV file has no header) may be added if it seems to be
 * required (no header).
 *
 * Open a Comma-Separated Values data file and read enough to
 * determine the variable names.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int probe_csv (const char *fname, char ***varnames,
	       int *nvars, gretlopt *opt)
{
    csvprobe probe = {0};
    int err;

    err = real_import_csv(fname, NULL, NULL, NULL, NULL,
			  &probe, NULL, *opt, NULL);

    if (!err) {
	int rerun = 0;

	err = probe_varnames_check(probe.dset, *opt, &rerun);

	if (err || rerun) {
	    destroy_dataset(probe.dset);
	    probe.dset = NULL;
	}

	if (!err && rerun) {
	    /* try again with --no-header flag */
	    *opt |= OPT_H;
	    err = real_import_csv(fname, NULL, NULL, NULL, NULL,
				  &probe, NULL, *opt, NULL);
	}

	if (!err) {
	    /* steal the varname array */
	    *varnames = probe.dset->varname;
	    *nvars = probe.dset->v;
	    probe.dset->varname = NULL;
	}

	destroy_dataset(probe.dset);
    }

    return err;
}

int csv_open_needs_matrix (gretlopt opt)
{
    int ret = 0;

    if (opt & OPT_M) {
	/* --rowmask=matrix */
	ret = 1;
    } else if (opt & OPT_F) {
	/* --fixed-cols=whatever */
	const char *s = get_optval_string(OPEN, OPT_F);

	ret = get_matrix_by_name(s) != NULL;
    }

    return ret;
}

typedef double keynum;

/* below: apparatus to implement the "join" command */

struct jr_row_ {
    int n_keys;     /* number of keys (needed for qsort callback) */
    keynum keyval;  /* primary key value */
    keynum keyval2; /* secondary key value, if applicable */
    int micro;      /* high-frequency "key", if any */
    int dset_row;   /* associated row in the RHS or outer dataset */
    double aux;     /* auxiliary value */
};

typedef struct jr_row_ jr_row;

struct obskey_ {
    char *timefmt; /* time format, as in strptime */
    int keycol;    /* the column holding the outer time-key */
    int m_means_q; /* "monthly means quarterly" */
    int numdates;  /* flag for conversion from numeric to string */
    int native;    /* native time-series info */
};

typedef struct obskey_ obskey;

struct joiner_ {
    int n_rows;     /* number of rows in data table */
    int n_keys;     /* number of keys used (0, 1 or 2) */
    int n_unique;   /* number of unique primary key values on right */
    jr_row *rows;   /* array of table rows */
    keynum *keys;   /* array of unique (primary) key values as 64-bit ints */
    int *key_freq;  /* counts of occurrences of (primary) key values */
    int *key_row;   /* record of starting row in joiner table for primary keys */
    int *str_keys;  /* flags for string comparison of key(s) */
    const int *l_keyno; /* list of key columns in left-hand dataset */
    const int *r_keyno; /* list of key columns in right-hand dataset */
    AggrType aggr;      /* aggregation method for 1:n joining */
    int seqval;         /* sequence number for aggregation */
    int auxcol;         /* auxiliary data column for aggregation */
    int midas_m;        /* midas frequency ratio */
    int midas_pd;       /* frequency of outer dataset */
    obskey *auto_keys;  /* struct to hold info on obs-based key(s) */
    DATASET *l_dset;    /* the left-hand or inner dataset */
    DATASET *r_dset;    /* the right-hand or outer temporary dataset */
};

typedef struct joiner_ joiner;

struct jr_filter_ {
    const char *expr;  /* expression to be run through "genr" */
    const double *val; /* (series) result of evaluating @expr */
    char *vname1;      /* first right-hand variable name */
    char *vname2;      /* second right-hand variable name */
    char *vname3;      /* third right-hand variable name */
};

typedef struct jr_filter_ jr_filter;

static int expand_jspec (joinspec *jspec, int addvars);

static void jr_filter_destroy (jr_filter *f)
{
    if (f != NULL) {
	free(f->vname1);
	free(f->vname2);
	free(f->vname3);
	free(f);
    }
}

static void joiner_destroy (joiner *jr)
{
    if (jr != NULL) {
	free(jr->rows);
	free(jr->keys);
	free(jr->key_freq);
	free(jr->key_row);
	free(jr);
    }
}

static joiner *joiner_new (int nrows)
{
    joiner *jr = malloc(sizeof *jr);

    if (jr != NULL) {
	jr->rows = calloc(nrows, sizeof *jr->rows);
	if (jr->rows == NULL) {
	    free(jr);
	    jr = NULL;
	}
    }

    if (jr != NULL) {
	jr->n_rows = nrows;
	jr->n_unique = 0;
	jr->keys = NULL;
	jr->key_freq = NULL;
	jr->key_row = NULL;
	jr->l_keyno = NULL;
	jr->r_keyno = NULL;
    }

    return jr;
}

static int real_set_outer_auto_keys (joiner *jr, const char *s,
				     int j, struct tm *tp)
{
    DATASET *l_dset = jr->l_dset;
    int err = 0;

    if (calendar_data(l_dset)) {
	int y, m, d, eday;

	y = tp->tm_year + 1900;
	m = tp->tm_mon + 1;
	d = tp->tm_mday;
	eday = epoch_day_from_ymd(y, m, d);
	if (eday < 0) {
	    if (s != NULL) {
		gretl_errmsg_sprintf("'%s' is not a valid date", s);
	    }
	    err = E_DATA;
	} else if (jr->n_keys == 2) {
	    /* use year and month */
	    jr->rows[j].n_keys = 2;
	    jr->rows[j].keyval = y;
	    jr->rows[j].keyval2 = m;
	    jr->rows[j].micro = 0;
	} else {
	    /* use epoch day */
	    jr->rows[j].n_keys = 1;
	    jr->rows[j].keyval = eday;
	    jr->rows[j].keyval2 = 0;
	    jr->rows[j].micro = 0;
	}
    } else {
	int lpd = l_dset->pd;
	int major = tp->tm_year + 1900;
	int minor = tp->tm_mon + 1;
	int micro = 0;

	if (jr->auto_keys->m_means_q) {
	    /* using the gretl-specific "%q" conversion */
	    if (minor > 4) {
		gretl_errmsg_sprintf("'%s' is not a valid date", s);
		err = E_DATA;
	    }
	} else if (lpd == 4) {
	    /* map from month on right to quarter on left, but
	       preserve the month info in case we need it
	    */
	    micro = minor;
	    minor = (int) ceil(minor / 3.0);
	}
	if (!err && micro == 0) {
	    micro = tp->tm_mday;
	}
	if (!err) {
	    jr->rows[j].n_keys = 2;
	    jr->rows[j].keyval = major;
	    jr->rows[j].keyval2 = minor;
	    jr->rows[j].micro = micro;
	}
    }

    return err;
}

static int set_time_format (obskey *auto_keys, const char *fmt)
{
    if (auto_keys->timefmt != NULL) {
	free(auto_keys->timefmt);
    }
    auto_keys->timefmt = gretl_strdup(fmt);
    return auto_keys->timefmt == NULL ? E_ALLOC : 0;
}

/* convert a numerical value to string for use with strptime */

static int numdate_to_string (char *targ, double x)
{
    if (na(x)) {
	return E_MISSDATA;
    } else {
	sprintf(targ, "%.16g", x);
	return 0;
    }
}

/* Parse a string from row @i of the outer dataset and set the
   key(s) on row @j of the joiner struct. The indices @i and @j may
   not be equal if a filter is being used. Note: we don't come
   here if the outer time-key column is subject to "tconvert"
   treatment; in that case we use read_iso_basic instead.
*/

static int read_outer_auto_keys (joiner *jr, int j, int i)
{
    char *tfmt = jr->auto_keys->timefmt;
    int numdates = jr->auto_keys->numdates;
    int tcol = jr->auto_keys->keycol;
    int pd = jr->l_dset->pd;
    struct tm t = {0};
    const char *s = NULL;
    char *test = NULL;
    char sconv[32];
    int s_src = 0;
    int err = 0;

    if (tcol >= 0) {
	/* using a specified column */
	if (numdates) {
	    /* column is numeric, conversion needed */
	    numdate_to_string(sconv, jr->r_dset->Z[tcol][i]);
	    s = sconv;
	    s_src = 1;
	} else {
	    /* column is string-valued, OK */
	    s = series_get_string_for_obs(jr->r_dset, tcol, i);
	    s_src = 2;
	}
    } else if (jr->auto_keys->native) {
	/* using native time-series info on right */
	ntolabel_8601(sconv, i, jr->r_dset);
	s = sconv;
	s_src = 1;
    } else {
	/* using first-column observation strings */
	s = jr->r_dset->S[i];
	s_src = 3;
    }

    if (s != NULL) {
	/* note: with strptime, a NULL return means that an error
	   occurred while a non-NULL and non-empty return string
	   means a trailing portion of the input was not
	   processed.
	*/
	test = strptime(s, tfmt, &t);
    }

    if (test == NULL || *test != '\0') {
	err = E_DATA;
	if (j == 0 && test != NULL && (pd == 12 || pd == 4 || pd == 1)) {
	    /* If we're looking at the first row of the filtered data,
	       allow for the possibility that we got "excess
	       precision", i.e. a daily date string when the left-hand
	       dataset is monthly, quarterly or annual.
	    */
	    char *chk = strptime(s, "%Y-%m-%d", &t);

	    if (chk != NULL && *chk == '\0') {
		set_time_format(jr->auto_keys, "%Y-%m-%d");
		err = 0; /* we might be OK, cancel the error for now */
	    }
	}
	if (err) {
	    gretl_errmsg_sprintf("'%s' does not match the format '%s'", s, tfmt);
	    fprintf(stderr, "time-format match error in read_outer_auto_keys:\n"
		    " remainder = '%s' (source = %s)\n", test ? test : "null",
		    s_src < 3 ? "specified time column" : "first-column strings");
	}
    }

    if (!err) {
	err = real_set_outer_auto_keys(jr, s, j, &t);
    }

    return err;
}

static int read_iso_basic (joiner *jr, int j, int i)
{
    int tcol = jr->auto_keys->keycol;
    double x;
    int err = 0;

    x = jr->r_dset->Z[tcol][i];

    if (na(x)) {
	err = E_MISSDATA;
    } else {
	int y = (int) floor(x / 10000);
	int m = (int) floor((x - 10000*y) / 100);
	int d = (int) (x - 10000*y - 100*m);
	guint32 ed = epoch_day_from_ymd(y, m, d);

	if (ed <= 0) {
	    gretl_errmsg_sprintf("'%.8g' is not a valid date", x);
	    err = E_DATA;
	} else if (calendar_data(jr->l_dset)) {
	    /* note: no need to go via struct tm */
	    jr->rows[j].n_keys = 1;
	    jr->rows[j].keyval = ed;
	    jr->rows[j].keyval2 = 0;
	    jr->rows[j].micro = 0;
	} else {
	    struct tm t = {0};

	    t.tm_year = y - 1900;
	    t.tm_mon = m - 1;
	    t.tm_mday = d;
	    err = real_set_outer_auto_keys(jr, NULL, j, &t);
	}
    }

    return err;
}

/* Evaluate the filter expression provided by the user, and if it
   works OK count the number of rows on which the filter returns
   non-zero.  Flag an error if the filter gives NA on any row, since
   it is then indeterminate.
*/

static int evaluate_filter (jr_filter *filter, DATASET *r_dset,
			    int *nrows)
{
    char *line;
    int i, err = 0;

    line = gretl_strdup_printf("filtered__=%s", filter->expr);
    if (line == NULL) {
	err = E_ALLOC;
    } else {
	err = generate(line, r_dset, GRETL_TYPE_SERIES,
		       OPT_P | OPT_Q, NULL);
    }

    if (!err) {
	int v = r_dset->v - 1;

	filter->val = r_dset->Z[v];
	*nrows = 0;

#if CDEBUG > 1
	fprintf(stderr, "filter genr: '%s':\n", line);
	for (i=0; i<r_dset->n; i++) {
	    fprintf(stderr, " %d: %g\n", i+1, filter->val[i]);
	}
#endif
	for (i=0; i<r_dset->n; i++) {
	    if (na(filter->val[i])) {
		gretl_errmsg_sprintf("join filter: indeterminate "
				     "value on row %d", i+1);
		err = E_MISSDATA;
		break;
	    } else if (filter->val[i] != 0.0) {
		*nrows += 1;
	    }
	}
    }

    free(line);

    return err;
}

static keynum dtoll (double x, int *err)
{
    if (na(x)) {
	*err = E_DATA;
	return -1;
    } else {
	return x;
    }
}

static keynum dtoll_full (double x, int key, int row, int *err)
{
    if (na(x)) {
	if (key == 2) {
	    gretl_errmsg_sprintf("%s: invalid secondary outer key value on row %d",
				 "join", row);
	} else {
	    gretl_errmsg_sprintf("%s: invalid (primary) outer key value on row %d",
				 "join", row);
	}
	*err = E_DATA;
	return -1;
    } else {
	return x;
    }
}

/* Determine whether or not row @i of the outer data satisfies the
   filter criterion; return 1 if the condition is met, 0 otherwise.
*/

static int join_row_wanted (jr_filter *filter, int i)
{
    int ret = filter->val[i] != 0;

#if CDEBUG > 2
    fprintf(stderr, "join filter: %s row %d\n",
	    ret ? "keeping" : "discarding", i);
#endif

    return ret;
}

static DATASET *outer_dataset (joinspec *jspec)
{
    if (jspec->c != NULL) {
	return jspec->c->dset;
    } else {
	return jspec->dset;
    }
}

#define using_auto_keys(j) (j->auto_keys->timefmt != NULL)

static joiner *build_joiner (joinspec *jspec,
			     DATASET *l_dset,
			     jr_filter *filter,
			     AggrType aggr,
			     int seqval,
			     obskey *auto_keys,
			     int n_keys,
			     int *err)
{
    joiner *jr = NULL;
    DATASET *r_dset = outer_dataset(jspec);
    int keycol  = jspec->colnums[JOIN_KEY];
    int valcol  = jspec->colnums[JOIN_TARG];
    int key2col = jspec->colnums[JOIN_KEY2];
    int auxcol  = jspec->colnums[JOIN_AUX];
    int i, nrows = r_dset->n;

#if CDEBUG
    fprintf(stderr, "joiner column numbers:\n"
	    "KEY %d, VAL %d, F1 %d, F2 %d, F3 %d, KEY2 %d, AUX %d\n",
	    keycol, valcol, jspec->colnums[JOIN_F1],
	    jspec->colnums[JOIN_F2], jspec->colnums[JOIN_F3],
	    key2col, auxcol);
#endif

    if (filter != NULL) {
	*err = evaluate_filter(filter, r_dset, &nrows);
	if (*err) {
	    return NULL;
	} else if (nrows == 0) {
	    gretl_warnmsg_set(_("No matching data after filtering"));
	    return NULL;
	}
    }

#if CDEBUG
    fprintf(stderr, "after filtering: dset->n = %d, nrows = %d\n",
	    r_dset->n, nrows);
#endif

    jr = joiner_new(nrows);

    if (jr == NULL) {
	*err = E_ALLOC;
    } else {
	double **Z = r_dset->Z;
	int use_iso_basic = 0;
	int j = 0;

	jr->aggr = aggr;
	jr->seqval = seqval;
	jr->auxcol = auxcol;
	jr->l_dset = l_dset;
	jr->r_dset = r_dset;
	jr->auto_keys = auto_keys;
	jr->n_keys = n_keys;
	jr->midas_m = 0;

	if (using_auto_keys(jr)) {
	    /* check for the case where the outer time-key
	       column is in the "tconvert" set: if so we
	       know it will be in YYYYMMDD format and we'll
	       give it special treatment
	    */
	    int tcol = jr->auto_keys->keycol;

	    if (tcol > 0 && jr->auto_keys->numdates) {
		if (column_is_timecol(jr->r_dset->varname[tcol])) {
		    use_iso_basic = 1;
		}
	    }
	}

#if CDEBUG
	fprintf(stderr, "use_iso_basic = %d\n", use_iso_basic);
#endif

	/* Now transcribe the data we want: we're pulling from the
	   full outer dataset and writing into the array of joiner row
	   structs. At this point we're applying the join filter (if
	   any) but are not doing any matching by key to the inner
	   dataset.
	*/

	for (i=0; i<r_dset->n && !*err; i++) {
	    if (filter != NULL && !join_row_wanted(filter, i)) {
		continue;
	    }
	    /* the keys */
	    if (use_iso_basic) {
		*err = read_iso_basic(jr, j, i);
	    } else if (using_auto_keys(jr)) {
		*err = read_outer_auto_keys(jr, j, i);
	    } else if (keycol > 0) {
		jr->rows[j].keyval = dtoll_full(Z[keycol][i], 1, i+1, err);
		if (!*err && key2col > 0) {
		    /* double key */
		    jr->rows[j].n_keys = 2;
		    jr->rows[j].keyval2 = dtoll_full(Z[key2col][i], 2, i+1, err);
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
	    /* "payload" data: record the dataset row */
	    jr->rows[j].dset_row = valcol > 0 ? i : -1;
	    /* the auxiliary data */
	    jr->rows[j].aux = auxcol > 0 ? Z[auxcol][i] : 0;
	    j++;
	}
    }

    return jr;
}

/* qsort callback for sorting rows of the joiner struct */

static int compare_jr_rows (const void *a, const void *b)
{
    const jr_row *ra = a;
    const jr_row *rb = b;
    int ret;

    ret = (ra->keyval > rb->keyval) - (ra->keyval < rb->keyval);

    if (ret == 0 && ra->n_keys > 1) {
	ret = (ra->keyval2 > rb->keyval2) - (ra->keyval2 < rb->keyval2);
    }

    if (ret == 0) {
	/* ensure stable sort */
	ret = a - b > 0 ? 1 : -1;
    }

    return ret;
}

/* Sort the rows of the joiner struct, by either one or two keys, then
   figure out how many unique (primary) key values we have and
   construct (a) an array of frequency of occurrence of these values
   and (b) an array which records the first row of the joiner on
   which each of these values is found.
*/

static int joiner_sort (joiner *jr)
{
    int matches = jr->n_rows;
    int i, err = 0;

    /* If there are string keys, we begin by mapping from the string
       indices on the right -- held in the keyval and/or keyval2
       members of the each joiner row -- to the indices for the same
       strings on the left. This enables us to avoid doing string
       comparisons when running aggr_value() later; we can just
       compare the indices of the strings. In addition, if on any
       given row we get no match for the right-hand key string on the
       left (signalled by a strmap value of -1) we can exploit this
       information by shuffling such rows to the end of the joiner
       rectangle and ignoring them when aggregating.
    */

    if (jr->str_keys[0] || jr->str_keys[1]) {
	series_table *stl, *str;
	int *strmap;
	int k, kmin, kmax, lkeyval, rkeyval;

	kmin = jr->str_keys[0] ? 1 : 2;
	kmax = jr->str_keys[1] ? 2 : 1;

	for (k=kmin; k<=kmax; k++) {
	    stl = series_get_string_table(jr->l_dset, jr->l_keyno[k]);
	    str = series_get_string_table(jr->r_dset, jr->r_keyno[k]);
	    strmap = series_table_map(str, stl);

	    if (strmap == NULL) {
		err = E_ALLOC;
		break;
	    }

	    for (i=0; i<jr->n_rows; i++) {
		if (k == 1) {
		    rkeyval = jr->rows[i].keyval;
		} else if (jr->rows[i].keyval == INT_MAX) {
		    continue;
		} else {
		    rkeyval = jr->rows[i].keyval2;
		}
		lkeyval = strmap[rkeyval];
#if CDEBUG > 1
		fprintf(stderr, "k = %d, row %d, keyval: %d -> %d\n", k, i, rkeyval, lkeyval);
#endif
		if (lkeyval > 0) {
		    if (k == 1) {
			jr->rows[i].keyval = lkeyval;
		    } else {
			jr->rows[i].keyval2 = lkeyval;
		    }
		} else {
		    /* arrange for qsort to move row to end */
		    jr->rows[i].keyval = G_MAXDOUBLE;
		    matches--;
		}
	    }

	    free(strmap);
	}
    }

    if (err) {
	return err;
    }

    qsort(jr->rows, jr->n_rows, sizeof *jr->rows, compare_jr_rows);

    if (matches < jr->n_rows) {
	jr->n_rows = matches;
    }

    jr->n_unique = 1;
    for (i=1; i<jr->n_rows; i++) {
	if (jr->rows[i].keyval != jr->rows[i-1].keyval) {
	    jr->n_unique += 1;
	}
    }

    jr->keys = malloc(jr->n_unique * sizeof *jr->keys);
    jr->key_freq = malloc(jr->n_unique * sizeof *jr->key_freq);
    jr->key_row = malloc(jr->n_unique * sizeof *jr->key_row);

    if (jr->keys == NULL || jr->key_freq == NULL || jr->key_row == NULL) {
	err = E_ALLOC;
    } else {
	int j = 0, nj = 1;

	for (i=0; i<jr->n_unique; i++) {
	    jr->key_freq[i] = 0;
	}

	jr->keys[0] = jr->rows[0].keyval;
	jr->key_row[0] = 0;

	for (i=1; i<jr->n_rows; i++) {
	    if (jr->rows[i].keyval != jr->rows[i-1].keyval) {
		/* finalize info for key j */
		jr->keys[j] = jr->rows[i-1].keyval;
		jr->key_freq[j] = nj;
		/* and initialize for next key */
		nj = 1;
		if (j < jr->n_unique - 1) {
		    jr->key_row[j+1] = i;
		}
		j++;
	    } else {
		nj++;
	    }
	}

	/* make sure the last row is right */
	jr->keys[j] = jr->rows[i-1].keyval;
	jr->key_freq[j] = nj;
    }

    return err;
}

#if CDEBUG > 1

static void joiner_print (joiner *jr)
{
    char **labels = NULL;
    jr_row *row;
    int i;

    if (jr->str_keys[0]) {
	labels = series_get_string_vals(jr->l_dset, jr->l_keyno[1], NULL, 0);
    }

    fprintf(stderr, "\njoiner: n_rows = %d\n", jr->n_rows);
    for (i=0; i<jr->n_rows; i++) {
	row = &jr->rows[i];
	if (row->n_keys > 1) {
	    fprintf(stderr, " row %d: keyvals=(%g,%g)\n",
		    i, row->keyval, row->keyval2);
	} else {
	    int k = lrint(row->keyval) - 1;

	    if (jr->str_keys[0] && row->keyval >= 0) {
		fprintf(stderr, " row %d: keyval=%g (%s)\n",
			i, row->keyval, labels[k]);
	    } else {
		fprintf(stderr, " row %d: keyval=%g\n",
			i, row->keyval);
	    }
	}
    }

    if (jr->keys != NULL) {
	fprintf(stderr, " for primary key: n_unique = %d\n", jr->n_unique);
	for (i=0; i<jr->n_unique; i++) {
	    fprintf(stderr,"  key value %g: count = %d\n",
		    jr->keys[i], jr->key_freq[i]);
	}
    }
}

# if CDEBUG > 2
static void print_outer_dataset (const DATASET *dset, const char *fname)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

    pprintf(prn, "Data extracted from %s:\n", fname);
    printdata(NULL, NULL, dset, OPT_O, prn);
    gretl_print_destroy(prn);
}
# endif

#endif

static int seqval_out_of_bounds (joiner *jr, int seqmax)
{
    if (jr->seqval < 0) {
	/* counting down from last match */
	return -jr->seqval > seqmax;
    } else {
	/* counting up from first match */
	return jr->seqval > seqmax;
    }
}

/* Do a binary search for the left-hand key value @targ in the sorted
   array of unique right-hand key values, @vals; return the position
   among @vals at which @targ matches, or -1 for no match.
*/

static int binsearch (keynum targ, const keynum *vals, int n, int offset)
{
    int m = n/2;

    if (fabs((targ) - (vals[m])) < 1.0e-7) {
	return m + offset;
    } else if (targ < vals[0] || targ > vals[n-1]) {
	return -1;
    } else if (targ < vals[m]) {
	return binsearch(targ, vals, m, offset);
    } else {
	return binsearch(targ, vals + m, n - m, offset + m);
    }
}

/* In some cases we can figure out what aggr_value() should return
   just based on the number of matches, @n, and the characteristics
   of the joiner. If so, write the value into @x and return 1; if
   not, return 0.
*/

static int aggr_val_determined (joiner *jr, int n, double *x, int *err)
{
    if (jr->aggr == AGGR_COUNT) {
	/* just return the number of matches */
	*x = n;
	return 1;
    } else if (jr->aggr == AGGR_SEQ && seqval_out_of_bounds(jr, n)) {
	/* out of bounds sequence index: return NA */
	*x = NADBL;
	return 1;
    } else if (n > 1 && jr->aggr == AGGR_NONE) {
	/* fail */
#if AGGDEBUG
	fprintf(stderr, "aggr_val_determined(): got n=%d\n", n);
#endif
	*err = E_DATA;
	gretl_errmsg_set(_("You need to specify an aggregation "
			   "method for a 1:n join"));
	*x = NADBL;
	return 1;
    } else {
	/* not enough information so far */
	return 0;
    }
}

/* get month-day index from @dset time-series info */

static int midas_day_index (int t, DATASET *dset)
{
    char obs[OBSLEN];
    int y, m, d, idx = -1;

    ntolabel(obs, t, dset);
    if (sscanf(obs, YMD_READ_FMT, &y, &m, &d) == 3) {
	idx = month_day_index(y, m, d, dset->pd);
    }

    return idx;
}

#define midas_daily(j) (j->midas_m > 20)

#define min_max_cond(x,y,a) ((a==AGGR_MAX && x>y) || (a==AGGR_MIN && x<y))

/* aggr_value: here we're working on a given row of the left-hand
   dataset. The values @key and (if applicable) @key2 are the
   left-hand keys for this row. We count the key-matches on the
   right and apply an aggregation procedure if the user specified
   one. We return the value that should be entered for the imported
   series on this row.

   Note: @xmatch and @auxmatch are workspace arrays allocated by
   the caller.
*/

static double aggr_value (joiner *jr,
			  keynum key1,
			  keynum key2,
			  int v, int revseq,
			  double *xmatch,
			  double *auxmatch,
			  int *nomatch,
			  int *err)
{
    double x, xa;
    int imin, imax, pos;
    int i, n, ntotal;

    /* find the position of the inner (primary) key in the
       array of unique outer key values */
    pos = binsearch(key1, jr->keys, jr->n_unique, 0);

#if AGGDEBUG
    if (pos < 0) {
	fprintf(stderr, " key1 = %g: no match\n", key1);
    } else {
	fprintf(stderr, " key1 = %g: matched at position %d\n", key1, pos);
    }
#endif

    if (pos < 0) {
	/* (primary) inner key value not found */
	*nomatch = 1;
	return jr->aggr == AGGR_COUNT ? 0 : NADBL;
    }

    /* how many matches at @pos? (must be at least 1 unless
       something very bad has happened)
    */
    n = jr->key_freq[pos];

#if AGGDEBUG
    fprintf(stderr, "  number of primary matches = %d (n_keys=%d)\n",
	    n, jr->n_keys);
#endif

    if (jr->n_keys == 1) {
	/* if there's just a single key, we can figure some
	   cases out already */
	if (aggr_val_determined(jr, n, &x, err)) {
	    return x;
	}
    }

    if (jr->key_row[pos] < 0) {
	/* "can't happen" */
	return NADBL;
    }

    /* set the range of rows for reading from the joiner rectangle */
    imin = jr->key_row[pos];
    imax = imin + n;

#if AGGDEBUG
    fprintf(stderr, "  aggregation row range: %d to %d\n", imin+1, imax);
#endif

    if (jr->aggr == AGGR_MIDAS) {
	/* special case: MIDAS "spreading" */
	int daily = dated_daily_data(jr->r_dset);
	int gotit = 0;

	x = NADBL;

	for (i=imin; i<imax && !gotit; i++) {
	    /* loop across primary key matches */
	    jr_row *r = &jr->rows[i];

	    if (jr->n_keys == 1 || key2 == r->keyval2) {
		/* got secondary key match */
		int sub, t = r->dset_row;
#if AGGDEBUG
		fprintf(stderr, "  i=%d: 2-key match: %d,%d (revseq=%d)\n",
			i, (int) key1, (int) key2, revseq);
#endif
		if (daily) {
		    /* outer dataset has known daily structure */
		    sub = midas_day_index(t, jr->r_dset);
		    gotit = sub == revseq;
		} else if (midas_daily(jr) && r->micro > 0) {
		    /* "other" daily data: r->micro holds day */
		    sub = month_day_index((int) key1, (int) key2,
					  r->micro, jr->midas_pd);
		    gotit = sub == revseq;
		} else {
		    if (r->micro > 0) {
			/* if present, this is derived from the outer
			   time-key specification
			*/
			sub = r->micro;
		    } else {
			date_maj_min(t, jr->r_dset, NULL, &sub);
		    }
		    gotit = (sub - 1) % jr->midas_m + 1 == revseq;
		}
		if (gotit) {
		    x = jr->r_dset->Z[v][t];
		}
	    }
	}

	/* and we're done */
	return x;
    }

    /* We now fill out the array @xmatch with non-missing values
       from the matching outer rows. If we have a secondary key
       we screen for matches on that as we go.
    */

    n = 0;      /* will now hold count of non-NA matches */
    ntotal = 0; /* will ignore the OK/NA distinction */

    for (i=imin; i<imax; i++) {
	jr_row *r = &jr->rows[i];

	if (jr->n_keys == 1 || key2 == r->keyval2) {
	    ntotal++;
	    x = jr->r_dset->Z[v][r->dset_row];
	    if (jr->auxcol) {
		xa = r->aux;
		if (!na(x) && na(xa)) {
		    /* we can't know the min/max of the aux var */
		    *err = E_MISSDATA;
		    return NADBL;
		}
		if (!na(xa)) {
		    auxmatch[n] = xa;
		    xmatch[n++] = x;
		}
	    } else if (!na(x)) {
		xmatch[n++] = x;
	    }
	}
    }

    if (jr->n_keys > 1) {
	/* we've already checked this for the 1-key case */
	if (aggr_val_determined(jr, n, &x, err)) {
	    return x;
	}
    }

    x = NADBL;

    if (n == 0) {
	; /* all matched observations are NA */
    } else if (jr->aggr == AGGR_NONE) {
	x = xmatch[0];
    } else if (jr->aggr == AGGR_SEQ) {
	int sval = jr->seqval;

	i = sval < 0 ? n + sval : sval - 1;
	if (i >= 0 && i < n) {
	    x = xmatch[i];
	}
    } else if (jr->aggr == AGGR_MAX || jr->aggr == AGGR_MIN) {
	if (jr->auxcol) {
	    /* using the max/min of an auxiliary var */
	    int idx = 0;

	    x = auxmatch[0];
	    for (i=1; i<n; i++) {
		if (min_max_cond(auxmatch[i], x, jr->aggr)) {
		    x = auxmatch[i];
		    idx = i;
		}
	    }
	    x = xmatch[idx];
	} else {
	    /* max/min of the actual data */
	    x = xmatch[0];
	    for (i=1; i<n; i++) {
		if (min_max_cond(xmatch[i], x, jr->aggr)) {
		    x = xmatch[i];
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

/* Handle the case where (a) the value from the right, @rz, is
   actually the coding of a string value, and (b) the LHS series is
   pre-existing and already has a string table attached. The RHS
   coding must be made consistent with that on the left. We reach this
   function only if we've verified that there are string tables on
   both sides, and that @rz is not NA.
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

static int get_inner_key_values (joiner *jr, int i,
				 const int *ikeyvars,
				 keynum *pk1, keynum *pk2,
				 int *missing)
{
    DATASET *l_dset = jr->l_dset;
    int err = 0;

    *pk1 = *pk2 = 0;

    if (using_auto_keys(jr)) {
	/* using the LHS dataset obs info */
	char obs[OBSLEN];

	ntolabel(obs, i, l_dset);
	if (calendar_data(l_dset)) {
	    guint32 ed = get_epoch_day(obs);

	    if (jr->n_keys == 2) {
		/* inner daily, outer monthly */
		int y, m, d;

		ymd_bits_from_epoch_day(ed, &y, &m, &d);
		*pk1 = y;
		*pk2 = m;
	    } else {
		*pk1 = ed;
	    }
	} else {
	    /* monthly or quarterly (FIXME any others?) */
	    *pk1 = atoi(obs);
	    *pk2 = atoi(obs + 5);
	}
    } else {
	/* using regular LHS key series */
	double dk1, dk2 = 0;
	keynum k1 = 0, k2 = 0;

	dk1 = l_dset->Z[ikeyvars[1]][i];
	if (jr->n_keys == 2) {
	    dk2 = l_dset->Z[ikeyvars[2]][i];
	}
	if (na(dk1) || na(dk2)) {
	    *missing = 1;
	} else {
	    k1 = dtoll(dk1, &err);
	    if (!err && jr->n_keys == 2) {
		k2 = dtoll(dk2, &err);
	    }
	}
	if (!err && !*missing) {
	    *pk1 = k1;
	    *pk2 = k2;
	}
    }

    return err;
}

/* Returns the 0-based column index in the outer dataset associated
   with the i^th target variable for joining. Note that the arg @i is
   1-based, being a position within a gretl list.
*/

static int outer_series_index (joinspec *jspec, int i)
{
    if (i >= 1) {
	return jspec->colnums[JOIN_TARG + i - 1];
    } else {
	return -1;
    }
}

static int aggregate_data (joiner *jr, const int *ikeyvars,
			   const int *targvars, joinspec *jspec,
			   int orig_v, int *modified)
{
    series_table *rst = NULL;
    series_table *lst = NULL;
    DATASET *dset = jr->l_dset;
    double *xmatch = NULL;
    double *auxmatch = NULL;
    keynum key, key2 = 0;
    int revseq = 0;
    int i, t, nmax;
    int err = 0;

    /* find the greatest (primary) key frequency */
    nmax = 0;
    for (i=0; i<jr->n_unique; i++) {
	if (jr->key_freq[i] > nmax) {
	    nmax = jr->key_freq[i];
	}
    }

#if AGGDEBUG
    fprintf(stderr, "\naggregate data: max primary matches = %d\n", nmax);
#endif

    if (jr->aggr == AGGR_MIDAS) {
	/* reverse sequence number for MIDAS join */
	revseq = targvars[0];
	jr->midas_m = revseq;
#if AGGDEBUG
	fprintf(stderr, "midas m = %d\n", jr->midas_m);
#endif
    } else if (nmax > 0) {
	/* allocate workspace for aggregation */
	int nx = (jr->auxcol > 0)? 2 * nmax : nmax;

	xmatch = malloc(nx * sizeof *xmatch);
	if (xmatch == NULL) {
	    return E_ALLOC;
	}
	if (jr->auxcol) {
	    auxmatch = xmatch + nmax;
	}
    }

    for (i=1; i<=targvars[0]; i++) {
	/* loop across the series to be added/modified */
	int rv, lv = targvars[i];
	int strcheck = 0;

	if (jr->aggr == AGGR_MIDAS) {
	    rv = jspec->colnums[JOIN_TARG];
	    jr->midas_pd = jspec->midas_pd;
	} else {
	    rv = outer_series_index(jspec, i);
	    jr->midas_pd = 0;
	}

	if (rv > 0) {
	    /* check for the case where both the target variable on the
	       left and the series to be imported are string-valued
	    */
	    rst = series_get_string_table(jr->r_dset, rv);
	    lst = series_get_string_table(jr->l_dset, lv);
	    strcheck = (rst != NULL && lst != NULL);
	}

	/* run through the rows in the current sample range of the
	   left-hand dataset, pick up the value of the inner key(s), and
	   call aggr_value() to determine the value that should be
	   imported from the right
	*/

	for (t=dset->t1; t<=dset->t2 && !err; t++) {
	    int missing = 0;
	    int nomatch = 0;
	    double z;

	    err = get_inner_key_values(jr, t, ikeyvars, &key, &key2, &missing);

	    if (err) {
		break;
	    } else if (missing) {
		dset->Z[lv][t] = NADBL;
		continue;
	    }

	    z = aggr_value(jr, key, key2, rv, revseq, xmatch, auxmatch,
			   &nomatch, &err);
#if AGGDEBUG
	    if (na(z)) {
		fprintf(stderr, " aggr_value: got NA (keys=%g,%g, err=%d)\n",
			key, key2, err);
	    } else {
		fprintf(stderr, " aggr_value: got %.12g (keys=%g,%g, err=%d)\n",
			z, key, key2, err);
	    }
#endif
	    if (!err && strcheck && !na(z)) {
		z = maybe_adjust_string_code(rst, lst, z, &err);
	    }
	    if (!err) {
		if (lv >= orig_v) {
		    /* @lv is a newly added series */
		    dset->Z[lv][t] = z;
		} else if (z != dset->Z[lv][t]) {
		    if (nomatch && !na(dset->Z[lv][t])) {
			; /* leave existing data alone (?) */
		    } else {
			dset->Z[lv][t] = z;
			*modified += 1;
		    }
		}
	    }
	}

	if (!err && jr->aggr == AGGR_MIDAS) {
	    /* set MIDAS-specific info on series @lv */
	    char label[MAXLABEL];

	    series_set_midas_period(dset, lv, revseq);
	    sprintf(label, "%s in sub-period %d",
		    jr->r_dset->varname[rv], revseq);
	    series_record_label(dset, lv, label);
	    series_set_midas_freq(dset, lv, jr->r_dset->pd);
	    if (i == 1) {
		series_set_midas_anchor(dset, lv);
	    }
	}

	revseq--;
    }

    free(xmatch);

    return err;
}

/* Simple transcription: we come here only if there are no keys, and
   we've verified that the number of rows on the right is no greater
   than the number of rows in the current sample range on the left.
*/

static int join_transcribe_data (joiner *jr, int lv, int newvar,
				 joinspec *jspec, int *modified)
{
    series_table *rst = NULL;
    series_table *lst = NULL;
    DATASET *dset = jr->l_dset;
    double zi;
    int strcheck = 0;
    int i, t, rv;
    int err = 0;

    rv = outer_series_index(jspec, 1);
    if (rv < 0) {
	return E_DATA;
    }

    rst = series_get_string_table(jr->r_dset, rv);
    lst = series_get_string_table(jr->l_dset, lv);
    strcheck = (rst != NULL && lst != NULL);

    for (i=0; i<jr->n_rows && !err; i++) {
	jr_row *r = &jr->rows[i];

	zi = jr->r_dset->Z[rv][r->dset_row];
	if (strcheck && !na(zi)) {
	    zi = maybe_adjust_string_code(rst, lst, zi, &err);
	}
	if (!err) {
	    t = dset->t1 + i;
	    if (newvar) {
		dset->Z[lv][t] = zi;
	    } else if (zi != dset->Z[lv][t]) {
		dset->Z[lv][t] = zi;
		*modified += 1;
	    }
	}
    }

    return err;
}

#include "tsjoin.c"

static int join_transcribe_multi_data (DATASET *l_dset,
				       DATASET *r_dset,
				       int *targlist,
				       int orig_v,
				       joinspec *jspec,
				       ts_joiner *tjr,
				       int *modified)
{
    series_table *rst = NULL;
    series_table *lst = NULL;
    int i, s, t, lv, rv;
    int t1, t2, m = 0;
    int strcheck, newvar;
    double xit;
    int err = 0;

    if (tjr == NULL) {
	t1 = l_dset->t1;
	t2 = l_dset->t2;
    } else {
	t1 = tjr->t1;
	t2 = tjr->t2;
    }

    for (i=1; i<=targlist[0] && !err; i++) {
	lv = targlist[i];
	rv = outer_series_index(jspec, i);
	if (rv < 0) {
	    gretl_errmsg_sprintf("join: '%s' not matched", l_dset->varname[lv]);
	    err = E_DATA;
	} else {
	    newvar = lv >= orig_v;
	    if (newvar) {
		strcheck = 0;
	    } else {
		rst = series_get_string_table(r_dset, rv);
		lst = series_get_string_table(l_dset, lv);
		strcheck = (rst != NULL && lst != NULL);
	    }
	    if (tjr != NULL) {
		m = tjr->rminor;
		s = tjr->rt1;
	    } else {
		s = 0;
	    }
	    for (t=t1; t<=t2 && !err; t++) {
		xit = r_dset->Z[rv][s];
		if (strcheck && !na(xit)) {
		    xit = maybe_adjust_string_code(rst, lst, xit, &err);
		}
		if (newvar) {
		    l_dset->Z[lv][t] = xit;
		} else if (xit != l_dset->Z[lv][t]) {
		    l_dset->Z[lv][t] = xit;
		    *modified += 1;
		}
		if (tjr != NULL) {
		    m = tj_continue(tjr, m, &s);
		} else {
		    s++;
		}
	    }
	}
    }

    return err;
}

static int join_simple_range_check (DATASET *l_dset,
				    DATASET *r_dset,
				    int *targlist)
{
    int err = 0;

    if (r_dset->n != sample_size(l_dset)) {
	gretl_errmsg_set("join: the observation ranges don't match");
	err = E_DATA;
    } else if (r_dset->v - 1 < targlist[0]) {
	gretl_errmsg_set("join: series missing on the right");
	err = E_DATA;
    }

    return err;
}

static jr_filter *join_filter_new (int *err)
{
    jr_filter *filter = malloc(sizeof *filter);

    if (filter == NULL) {
	*err = E_ALLOC;
    } else {
	filter->expr = NULL;
	filter->val = NULL;
	filter->vname1 = NULL;
	filter->vname2 = NULL;
	filter->vname3 = NULL;
    }

    return filter;
}

#if CDEBUG

static void print_filter_vnames (jr_filter *f)
{
    if (f == NULL) return;

    if (f->vname1 != NULL) {
	fprintf(stderr, "filter varname 1 (target %d, JOIN_F1): %s\n",
		JOIN_F1, f->vname1);
    }
    if (f->vname2 != NULL) {
	fprintf(stderr, "filter varname 2 (target %d, JOIN_F2): %s\n",
		JOIN_F2, f->vname2);
    }
    if (f->vname3 != NULL) {
	fprintf(stderr, "filter varname 3 (target %d, JOIN_F3): %s\n",
		JOIN_F3, f->vname3);
    }
}

#endif

/* Allocate the filter struct and crawl along the filter expression,
   @s, looking for up to three names of right-hand columns. We need to
   record these names so we can be sure to read the associated column
   data, or else the filter won't work. We could increase the maximum
   number of column-names to store but probably 3 is enough.

   The heuristic for column-name detection is that we find a portion
   of @s which is legal as a gretl identifier but which is not
   enclosed in quotes, is not directly followed by '(', and is not
   the name of a scalar or string variable on the left.
*/

static jr_filter *make_join_filter (const char *s, int *err)
{
    jr_filter *filter = join_filter_new(err);

    if (filter != NULL) {
	char test[VNAMELEN];
	int n, ngot = 0;

	filter->expr = s;

	while (*s && ngot < 3) {
	    if (*s == '"') {
		/* skip double-quoted stuff */
		s = strchr(s + 1, '"');
		if (s == NULL) {
		    *err = E_PARSE;
		    break;
		} else {
		    s++;
		    if (*s == '\0') {
			break;
		    }
		}
	    }
	    n = gretl_namechar_spn(s);
	    if (n > 0) {
		if (n < VNAMELEN && s[n] != '(') {
		    *test = '\0';
		    strncat(test, s, n);
		    if (!gretl_is_scalar(test) && !gretl_is_string(test)) {
			if (ngot == 0) {
			    filter->vname1 = gretl_strdup(test);
			} else if (ngot == 1) {
			    filter->vname2 = gretl_strdup(test);
			} else if (ngot == 2) {
			    filter->vname3 = gretl_strdup(test);
			}
			ngot++;
		    }
		}
		s += n;
	    } else {
		s++;
	    }
	}
    }

#if CDEBUG
    print_filter_vnames(filter);
#endif

    return filter;
}

/* Add series to hold the join data.  We come here only if the
   target series is/are not already present in the left-hand
   dataset.
*/

static int add_target_series (const char **vnames,
			      DATASET *dset,
			      int *targvars,
			      int n_add)
{
    int i, v = dset->v;
    int err;

    err = dataset_add_NA_series(dset, n_add);

    if (!err) {
	for (i=1; i<=targvars[0]; i++) {
	    if (targvars[i] < 0) {
		strcpy(dset->varname[v], vnames[i-1]);
		targvars[i] = v++;
	    }
	}
    }

    return err;
}

/* Parse either one or two elements (time-key column name, time format
   string) out of @s. If both elements are present they must be
   comma-separated; if only the second element is present it must be
   preceded by a comma.
*/

static int process_time_key (const char *s, char *tkeyname,
			     char *tkeyfmt)
{
    int err = 0;

    if (*s == ',') {
	/* only time-key format supplied */
	strncat(tkeyfmt, s + 1, 31);
    } else {
	const char *p = strchr(s, ',');

	if (p == NULL) {
	    /* only column name given */
	    strncat(tkeyname, s, VNAMELEN - 1);
	} else {
	    int n = p - s;

	    n = (n > VNAMELEN - 1)? VNAMELEN - 1 : n;
	    strncat(tkeyname, s, n);
	    strncat(tkeyfmt, p + 1, 31);
	}
    }

#if TDEBUG
    fprintf(stderr, "time key: name='%s', format='%s'\n",
	    tkeyname, tkeyfmt);
#endif

    if (*tkeyname == '\0' && *tkeyfmt == '\0') {
	err = E_DATA;
    }

    return err;
}

/* Parse either one column-name or two comma-separated names out of
   @s. If @s contains a comma, we accept a zero-length name on either
   the left or the right -- but not both -- as indicating that we
   should use the corresponding inner key name.
*/

static int process_outer_key (const char *s, int n_keys,
			      char *name1, char *name2,
			      gretlopt opt)
{
    int n_okeys = 0;
    int err = 0;

    if (strchr(s, ',') == NULL) {
	/* just one outer key */
	strncat(name1, s, VNAMELEN - 1);
	n_okeys = 1;
    } else {
	/* two comma-separated keys */
	int n2 = 0, n1 = strcspn(s, ",");

	if (n1 >= VNAMELEN) {
	    err = E_PARSE;
	} else {
	    strncat(name1, s, n1);
	    s += n1 + 1;
	    n2 = strlen(s);
	    if (n2 >= VNAMELEN) {
		err = E_PARSE;
	    } else {
		strcat(name2, s);
	    }
	}

	if (!err && n1 == 0 && n2 == 0) {
	    /* both fields empty: wrong */
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

/* Handle the case where the user gave a "%q" or "%Q" conversion
   specifier (which we take to mean quarter). We convert this to %m
   for use with strptime(), but record that the fact that "month means
   quarter".
*/

static int format_uses_quarterly (char *fmt)
{
    char *s = fmt;
    int i, ret = 0;

    for (i=0; s[i]; i++) {
	if (s[i] == '%' &&
	    (s[i+1] == 'q' || s[i+1] == 'Q') &&
	    (i == 0 || s[i-1] != '%')) {
	    s[i+1] = 'm';
	    ret = 1;
	}
    }

    return ret;
}

static int check_for_quarterly_format (obskey *auto_keys, int pd)
{
    char *s = auto_keys->timefmt;
    int i, err = 0;

#if CDEBUG
    fprintf(stderr, "check_for_quarterly_format: '%s'\n", s);
#endif

    for (i=0; s[i]; i++) {
	if (s[i] == '%' &&
	    (s[i+1] == 'q' || s[i+1] == 'Q') &&
	    (i == 0 || s[i-1] != '%')) {
	    if (pd == 4 || pd == 1) {
		s[i+1] = 'm';
		auto_keys->m_means_q = 1;
	    } else {
		err = E_DATA;
		gretl_errmsg_sprintf("The '%c' format is not applicable "
				     "for data with frequency %d",
				     s[i+1], pd);
	    }
	    break;
	}
    }

    return err;
}

/* for use in determining optimal auto keys */
#define daily(d) (d->pd >= 5 && d->pd <= 7)

/* time-series data on the left, and no explicit keys supplied */

static int auto_keys_check (const DATASET *l_dset,
			    const DATASET *r_dset,
			    gretlopt opt,
			    const char *tkeyfmt,
			    obskey *auto_keys,
			    int *n_keys,
			    int *do_tsjoin)
{
    int lpd = l_dset->pd;
    int rpd = 0;
    int err = 0;

    if (!dataset_is_time_series(l_dset)) {
	/* On the left we need a time-series dataset */
	err = E_DATA;
	goto bailout;
    }

    if (dataset_is_time_series(r_dset)) {
	rpd = r_dset->pd;
	if (use_tsjoin(l_dset, r_dset)) {
	    *do_tsjoin = 1;
	    return 0;
	}
	auto_keys->native = 1;
    } else if (r_dset->S == NULL && auto_keys->keycol < 0) {
	/* On the right, we need either obs strings or a specified
	   time column
	*/
	err = E_DATA;
	goto bailout;
    }

    if (*tkeyfmt != '\0') {
	/* the user supplied a time-format spec */
	err = set_time_format(auto_keys, tkeyfmt);
	if (!err) {
	    err = check_for_quarterly_format(auto_keys, lpd);
	}
	if (!err) {
	    if (annual_data(l_dset)) {
		*n_keys = 1;
	    } else if (calendar_data(l_dset)) {
		*n_keys = 1;
	    } else {
		*n_keys = 2;
	    }
	}
    } else {
	/* default formats */
	if (calendar_data(l_dset)) {
	    err = set_time_format(auto_keys, "%Y-%m-%d");
	    if (daily(l_dset) && rpd == 12) {
		*n_keys = 2; /* use year and month */
	    } else if (!err) {
		*n_keys = 1; /* use epoch day */
	    }
	} else if (lpd == 12) {
	    err = set_time_format(auto_keys, "%Y-%m");
	    if (!err) {
		*n_keys = 2;
	    }
	} else if (annual_data(l_dset)) {
	    err = set_time_format(auto_keys, "%Y");
	    if (!err) {
		*n_keys = 1;
	    }
	} else if (lpd == 4) {
	    /* try "excess precision" ISO daily? */
	    err = set_time_format(auto_keys, "%Y-%m-%d");
	    if (!err) {
		*n_keys = 2;
	    }
	} else {
	    err = E_PDWRONG;
	}
    }

 bailout:

    /* we should flag an error here only if the user
       explicitly requested use of this apparatus,
       by giving --tkey=<whatever> (OPT_K)
    */
    if (err && !(opt & OPT_K)) {
	err = 0;
    }

    return err;
}

static int make_time_formats_array (char const **fmts, char ***pS)
{
    char **S = strings_array_new(2);
    int i, err = 0;

    if (S == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<2 && !err; i++) {
	    if (fmts[i] != NULL) {
		S[i] = gretl_strdup(fmts[i]);
		if (S[i] == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    }

    if (err) {
	strings_array_free(S, 2);
    } else {
	*pS = S;
    }

    return err;
}

/* Crawl along the string containing comma-separated names of columns
   that fall under the "tconvert" option. For each name, look for a
   match among the names of columns selected from the CSV file for
   some role in the join operation (data, key or whatever). It is
   not considered an error if there's no match for a given "tconvert"
   name; in that case the column will be ignored.

   If @tconvfmt is non-NULL, this should hold a common format for
   date conversion.

   If @keyfmt is not an empty string, it holds a specific format
   for the --tkey column. In that case we should check for tkey
   among the tconvert columns; if it's present we need to add its
   format to the timeconv_map apparatus (as an override to the
   common format, or an addendum if no common format is given).
*/

static int process_tconvert_info (joinspec *jspec,
				  const char *tconvcols,
				  const char *tconvfmt,
				  const char *keyfmt)
{
    int *list = NULL;
    char **names = NULL;
    char **fmts = NULL;
    const char *colname;
    const char *tkeyfmt = NULL;
    char *tkeyname = NULL;
    int nnames = 0;
    int i, j, err = 0;

    names = gretl_string_split(tconvcols, &nnames, " ,");
    if (names == NULL) {
	err = E_ALLOC;
    }

    /* match the names we got against the set of "wanted" columns */

    for (i=0; i<nnames && !err; i++) {
	for (j=0; j<jspec->ncols; j++) {
	    colname = jspec->colnames[j];
	    if (colname != NULL && !strcmp(names[i], colname)) {
		gretl_list_append_term(&list, j);
		if (list == NULL) {
		    err = E_ALLOC;
		} else if (*keyfmt != '\0' && j == JOIN_KEY) {
		    /* we've got the time-key variable here,
		       and a format has been given for it
		    */
		    tkeyfmt = keyfmt;
		    tkeyname = names[i];
		}
		break;
	    }
	}
    }

    /* allocate and record the time-format info, if any */
    if (!err && list != NULL && (tconvfmt != NULL || tkeyfmt != NULL)) {
	char const *tmp[2] = {tconvfmt, tkeyfmt};

	err = make_time_formats_array(tmp, &fmts);
    }

#if CDEBUG
    printlist(list, "timeconv list");
#endif

    if (!err && list != NULL) {
	timeconv_map_set(nnames, names, tkeyname, fmts);
    }

    if (err) {
	/* clean up if need be */
	strings_array_free(names, nnames);
	free(list);
    } else {
	jspec->timecols = list;
    }

    return err;
}

static void obskey_init (obskey *keys)
{
    keys->timefmt = NULL;
    keys->keycol = -1;
    keys->m_means_q = 0;
    keys->numdates = 0;
    keys->native = 0;
}

#define lr_mismatch(l,r) ((l > 0 && r == 0) || (r > 0 && l == 0))

/* Run a check pertaining to the nature of the "payload"
   (string-valued vs numeric) in relation to the aggregation
   method specified and the nature of the existing left-hand
   series, if any.
*/

static int join_data_type_check (joinspec *jspec,
				 const DATASET *l_dset,
				 int *targvars,
				 AggrType aggr)
{
    DATASET *r_dset = outer_dataset(jspec);
    int lstr, rstr;
    int i, vl, vr;
    int err = 0;

    for (i=1; i<=targvars[0]; i++) {
	lstr = rstr = -1;
	vl = targvars[i];
	vr = outer_series_index(jspec, i);
	if (vl > 0) {
	    /* there's an existing LHS series */
	    lstr = is_string_valued(l_dset, vl);
	    if (lstr && aggr == AGGR_COUNT) {
		/* count values can't be mixed with strings */
		err = E_TYPES;
	    }
	}
	if (!err && vr > 0) {
	    /* there's a payload variable on the right */
	    rstr = is_string_valued(r_dset, vr);
	}
	if (!err && lr_mismatch(lstr, rstr)) {
	    /* one of (L, R) is string-valued, but not the other */
	    err = E_TYPES;
	}
    }

    return err;
}

static int aggregation_type_check (joinspec *jspec, AggrType aggr)
{
    int err = 0;

    if (aggr == AGGR_NONE || aggr == AGGR_COUNT || aggr == AGGR_SEQ) {
	; /* no problem */
    } else {
	/* the aggregation method requires numerical data: flag
	   an error if we got strings instead
	*/
	const DATASET *dset = outer_dataset(jspec);
	int aggcol = 0;

	if (jspec->colnums[JOIN_AUX] > 0) {
	    aggcol = jspec->colnums[JOIN_AUX];
	} else if (jspec->colnums[JOIN_TARG] > 0) {
	    aggcol = jspec->colnums[JOIN_TARG];
	}

	if (aggcol > 0 && is_string_valued(dset, aggcol)) {
	    gretl_errmsg_sprintf("'%s' is a string variable: aggregation type "
				 "is not applicable", dset->varname[aggcol]);
	    err = E_TYPES;
	}
    }

    return err;
}

static int check_for_missing_columns (joinspec *jspec)
{
    const char *name;
    int i;

    /* Note: it's possible, though we hope unlikely, that our
       heuristic for extracting column names from the filter
       expression (if present) gave a false positive. In that case the
       name at position JOIN_F* might be spuriously "missing": not
       found, but not really wanted. To guard against this eventuality
       we'll skip the check for the JOIN_F* columns here. If either
       one is really missing, that will show up before long, when the
       filter is evaluated.
    */

    for (i=0; i<jspec->ncols; i++) {
	if (i == JOIN_F1 || i == JOIN_F2 || i == JOIN_F3) {
	    continue;
	}
	name = jspec->colnames[i];
	if (name != NULL && jspec->colnums[i] == 0) {
	    gretl_errmsg_sprintf(_("%s: column '%s' was not found"), "join", name);
	    return E_DATA;
	}
#if CDEBUG
	if (name != NULL) {
	    fprintf(stderr, "colname '%s' -> colnum %d\n", name, jspec->colnums[i]);
	}
#endif
    }

    return 0;
}

/* Handle the situation where for one or other of the keys
   there's a type mismatch: string on left but not on right, or
   vice versa.
*/

static int key_types_error (int lstr, int rstr)
{
    if (lstr) {
	gretl_errmsg_sprintf(_("%s: string key on left but numeric on right"), "join");
    } else {
	gretl_errmsg_sprintf(_("%s: string key on right but numeric on left"), "join");
    }

    return E_TYPES;
}

static int set_up_outer_keys (joinspec *jspec, const DATASET *l_dset,
			      gretlopt opt, const int *ikeyvars,
			      int *okeyvars, obskey *auto_keys,
			      int *str_keys)
{
    const DATASET *r_dset = outer_dataset(jspec);
    int lstr, rstr;
    int err = 0;

    if (jspec->colnums[JOIN_KEY] > 0) {
	okeyvars[0] += 1;
	okeyvars[1] = jspec->colnums[JOIN_KEY];
    }

    if (jspec->colnums[JOIN_KEY2] > 0) {
	okeyvars[0] += 1;
	okeyvars[2] = jspec->colnums[JOIN_KEY2];
    }

    if (opt & OPT_K) {
	/* time key on right */
	rstr = is_string_valued(r_dset, okeyvars[1]);
	auto_keys->keycol = okeyvars[1];
	if (!rstr) {
	    /* flag the need to convert to string later */
	    auto_keys->numdates = 1;
	}
    } else {
	/* regular key(s) on right */
	lstr = is_string_valued(l_dset, ikeyvars[1]);
	rstr = is_string_valued(r_dset, okeyvars[1]);

	if (lstr != rstr) {
	    err = key_types_error(lstr, rstr);
	} else if (lstr) {
	    str_keys[0] = 1;
	}

	if (!err && okeyvars[2] > 0) {
	    lstr = is_string_valued(l_dset, ikeyvars[2]);
	    rstr = is_string_valued(r_dset, okeyvars[2]);
	    if (lstr != rstr) {
		err = key_types_error(lstr, rstr);
	    } else if (lstr) {
		str_keys[1] = 1;
	    }
	}
    }

    return err;
}

static int first_available_index (joinspec *jspec)
{
    int i;

    for (i=JOIN_TARG; i<jspec->ncols; i++) {
	if (jspec->colnums[i] == 0) {
	    return i;
	}
    }

    return -1;
}

static int determine_gdt_matches (const char *fname,
				  joinspec *jspec,
				  int **plist,
				  int *addvars,
				  int *omm,
				  PRN *prn)
{
    char **vnames = NULL;
    int nv = 0;
    int err = 0;

    err = gretl_read_gdt_varnames(fname, &vnames, &nv);

    if (!err) {
	GPatternSpec *pspec;
	int *vlist = NULL;
	char **S = NULL;
	int i, j, ns = 0;

	/* form array of unique wanted identifiers */
	for (i=0; i<jspec->ncols && !err; i++) {
	    if (jspec->colnames[i] != NULL) {
		err = strings_array_add_uniq(&S, &ns, jspec->colnames[i], NULL);
	    }
	}

	for (i=0; i<ns && !err; i++) {
	    int match = 0;

	    if (!strcmp(S[i], "$obsmajor")) {
		omm[0] = 1;
		continue;
	    } else if (!strcmp(S[i], "$obsminor")) {
		omm[1] = 1;
		continue;
	    }

	    if (prn != NULL) {
		pprintf(prn, "checking for '%s'\n", S[i]);
	    }

	    if (is_wildstr(S[i])) {
		pspec = g_pattern_spec_new(S[i]);
		for (j=1; j<nv; j++) {
		    if (g_pattern_match_string(pspec, vnames[j])) {
			match++;
			if (!in_gretl_list(vlist, j)) {
			    gretl_list_append_term(&vlist, j);
			}
		    }
		}
		g_pattern_spec_free(pspec);
	    } else {
		for (j=1; j<nv; j++) {
		    if (!strcmp(S[i], vnames[j])) {
			match = 1;
			if (!in_gretl_list(vlist, j)) {
			    gretl_list_append_term(&vlist, j);
			}
			break;
		    }
		}
		if (!match) {
		    err = E_DATA;
		    gretl_errmsg_sprintf("'%s': not found", S[i]);
		}
	    }
	    if (prn != NULL) {
		pprintf(prn, " found %d match(es)\n", match);
	    }
	}

	if (!err && (vlist == NULL || vlist[0] == 0)) {
	    gretl_errmsg_set("No matching data were found");
	    err = E_DATA;
	}

	if (!err) {
	    *addvars = vlist[0] - ns;
	    *plist = vlist;
	} else {
	    free(vlist);
	}

	strings_array_free(S, ns);
	strings_array_free(vnames, nv);
    }

    return err;
}

static int rhs_add_obsmajmin (int *omm, DATASET *dset)
{
    int err = dataset_add_series(dset, omm[0] + omm[1]);

    if (!err) {
	int t, maj, min, vmaj = 0, vmin = 0;
	int *pmaj = NULL, *pmin = NULL;

	if (omm[0]) {
	    pmaj = &maj; vmaj = dset->v - 1 - omm[1];
	    strcpy(dset->varname[vmaj], "$obsmajor");
	}
	if (omm[1]) {
	    pmin = &min; vmin = dset->v - 1;
	    strcpy(dset->varname[vmin], "$obsminor");
	}
	for (t=0; t<dset->n; t++) {
	    date_maj_min(t, dset, pmaj, pmin);
	    if (vmaj) {
		dset->Z[vmaj][t] = maj;
	    }
	    if (vmin) {
		dset->Z[vmin][t] = min;
	    }
	}
    }

    return err;
}

static int join_import_gdt (const char *fname,
			    joinspec *jspec,
			    gretlopt opt,
			    PRN *prn)
{
    const char *cname;
    int *vlist = NULL;
    int orig_ncols = jspec->ncols;
    int i, vi, addvars = 0;
    int omm[2] = {0};
    int err = 0;

    err = determine_gdt_matches(fname, jspec, &vlist, &addvars,
				omm, prn);

    if (!err) {
	jspec->dset = datainfo_new();
	if (jspec->dset == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_read_gdt_subset(fname, jspec->dset, vlist, opt);
    }

    if (!err && addvars > 0) {
	/* we have some extra vars due to wildcard expansion */
	err = expand_jspec(jspec, addvars);
    }

    if (!err && (omm[0] || omm[1])) {
	/* we need to add $obsmajor and/or $obsminor on the right */
	err = rhs_add_obsmajmin(omm, jspec->dset);
    }

    if (!err) {
	/* match up the imported series with their roles */
	for (i=0; i<orig_ncols && !err; i++) {
	    cname = jspec->colnames[i];
	    if (cname != NULL && !is_wildstr(cname)) {
		vi = current_series_index(jspec->dset, cname);
		if (vi < 0) {
		    err = E_DATA;
		} else {
		    jspec->colnums[i] = vi;
		}
	    }
	}
    }

    if (!err && addvars > 0) {
	/* register any extra imported series */
	int j, pos, idx;

	for (i=1; i<jspec->dset->v; i++) {
	    pos = 0;
	    for (j=0; j<jspec->ncols; j++) {
		if (jspec->colnums[j] == i) {
		    /* already registered */
		    pos = i;
		    break;
		}
	    }
	    if (pos == 0) {
		idx = first_available_index(jspec);
		if (idx < 0) {
		    err = E_DATA;
		} else {
		    jspec->colnums[idx] = i;
		    jspec->colnames[idx] = jspec->dset->varname[i];
		}
	    }
	}
    }

    free(vlist);

    return err;
}

/* add/replace an entry in @jspec's colnames array,
   and record the fact that @src is now owned by
   @jspec so we can free it when we're finished
   and hence avoid leaking memory
*/

static int jspec_push_tmpname (joinspec *jspec,
			       int pos,
			       char *src)
{
    jspec->colnames[pos] = src;
    return strings_array_donate(&jspec->tmpnames, &jspec->n_tmp, src);
}

static int determine_csv_matches (const char *fname,
				  joinspec *jspec,
				  PRN *prn)
{
    gretlopt opt = OPT_NONE; /* FIXME? */
    char **vnames = NULL;
    int nv = 0;
    int err = 0;

    err = probe_csv(fname, &vnames, &nv, &opt);

    if (!err) {
	GPatternSpec *pspec;
	int nmatch = 0;
	int i, err = 0;

	pspec = g_pattern_spec_new(jspec->colnames[JOIN_TARG]);

	/* first determine the number of matches to @pspec */
	for (i=0; i<nv; i++) {
	    if (g_pattern_match_string(pspec, vnames[i])) {
		nmatch++;
	    }
	}

	if (nmatch > 1) {
	    /* we have some extra vars due to wildcard expansion */
	    err = expand_jspec(jspec, nmatch - 1);
	}

	if (!err && nmatch > 0) {
	    int j = JOIN_TARG;

	    for (i=0; i<nv; i++) {
		if (g_pattern_match_string(pspec, vnames[i])) {
		    jspec_push_tmpname(jspec, j++, vnames[i]);
		    vnames[i] = NULL;
		}
	    }
	}

	g_pattern_spec_free(pspec);
	strings_array_free(vnames, nv);
    }

    return err;
}

static int join_import_csv (const char *fname,
			    joinspec *jspec,
			    gretlopt opt,
			    PRN *prn)
{
    int err = 0;

    if (jspec->wildcard) {
	err = determine_csv_matches(fname, jspec, prn);
	if (err) {
	    pputs(prn, "join_import_csv: failed in matching varnames\n");
	}
    }

    if (!err) {
	err = real_import_csv(fname, NULL, NULL, NULL, jspec,
			      NULL, NULL, opt, prn);
	if (0 && !err) {
	    /* question, 2021-01-09: this is zeroed out: why? */
	    DATASET *dset = jspec->c->dset;
	    int pd, reversed = 0;

	    fprintf(stderr, "join_import_csv: n=%d, v=%d, pd=%d, markers=%d\n",
		    dset->n, dset->v, dset->pd, dset->markers);
	    pd = test_markers_for_dates(dset, &reversed, NULL, prn);
	    fprintf(stderr, "pd from markers: %d\n", pd);
	}
    }

    return err;
}

static int join_range_check (joiner *jr, DATASET *dset, AggrType aggr)
{
    int T = sample_size(dset);
    int err = 0;

    if (jr->n_rows < T) {
	/* can we handle jr->n_rows < T? */
	if (aggr != AGGR_NONE) {
	    err = E_DATA; /* No */
	} else if (dset->v == 1) {
	    ; /* only const, OK */
	} else if (dset->v == 2 && !strcmp(dset->varname[1], "index")) {
	    ; /* "nulldata" default, OK */
	} else {
	    err = E_DATA; /* Not OK */
	}
    } else if (jr->n_rows != T) {
	gretl_errmsg_set(_("Series length does not match the dataset"));
	err = E_DATA;
    }

    if (err) {
	gretl_errmsg_set(_("Series length does not match the dataset"));
    }

    return err;
}

static int *get_series_indices (const char **vnames,
				int nvars,
				DATASET *dset,
				int *n_add,
				int *any_wild,
				int *err)
{
    int *ret = gretl_list_new(nvars);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	int i, v;

	for (i=0; i<nvars && !*err; i++) {
	    if (is_wildstr(vnames[i])) {
		*any_wild = 1;
		continue;
	    }
	    v = current_series_index(dset, vnames[i]);
	    if (v == 0) {
		*err = E_DATA;
	    } else {
		ret[i+1] = v;
		if (v < 0) {
		    *n_add += 1;
		}
	    }
	}
    }

    if (!*err && nvars > 1 && *any_wild) {
	/* wildcard spec must be singleton varname spec */
	gretl_errmsg_set(_("Invalid join specification"));
	*err = E_DATA;
    }

    if (*err) {
	free(ret);
	ret = NULL;
    }

    return ret;
}

/* figure the number of series to import */

static int jspec_n_vars (joinspec *jspec)
{
    return jspec->ncols - JOIN_TARG;
}

/* we come here if we have determined that the
   import series specification includes a wildcard
   ('*' or '?')
*/

static int *revise_series_indices (joinspec *jspec,
				   DATASET *dset,
				   int *n_add,
				   int *err)
{
    int nvars = jspec_n_vars(jspec);
    int *ret = gretl_list_new(nvars);

    ret = gretl_list_new(nvars);

    if (ret == NULL) {
	*err = E_ALLOC;
    } else {
	jspec->wildnames = strings_array_new(nvars);
	if (jspec->wildnames == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	const char *cname;
	int i, v;

	/* zero the count of added vars */
	*n_add = 0;

	for (i=0; i<nvars && !*err; i++) {
	    cname = jspec->colnames[JOIN_TARG + i];
	    v = current_series_index(dset, cname);
	    if (v == 0) {
		*err = E_DATA;
	    } else {
		ret[i+1] = v;
		if (v < 0) {
		    /* not a current series */
		    if (gretl_type_from_name(cname, NULL)) {
			*err = E_TYPES;
		    } else {
			*n_add += 1;
		    }
		}
		if (!*err) {
		    jspec->wildnames[i] = gretl_strdup(cname);
		}
	    }
	}
    }

    if (*err) {
	free(ret);
	ret = NULL;
    }

    return ret;
}

static int set_up_jspec (joinspec *jspec,
			 const char **vnames,
			 int nvars,
			 gretlopt opt,
			 int any_wild,
			 AggrType aggr,
			 int midas_pd)
{
    int i, j, ncols = JOIN_TARG + nvars;

    jspec->colnames = malloc(ncols * sizeof *jspec->colnames);
    jspec->colnums = malloc(ncols * sizeof *jspec->colnums);

    if (jspec->colnames == NULL || jspec->colnums == NULL) {
	return E_ALLOC;
    }

    jspec->ncols = ncols;
    jspec->wildcard = any_wild;
    jspec->wildnames = NULL;
    jspec->tmpnames = NULL;
    jspec->n_tmp = 0;
    jspec->mdsbase = NULL;
    jspec->mdsnames = NULL;
    jspec->auto_midas = 0;
    jspec->midas_pd = 0;

    if (aggr == AGGR_MIDAS) {
	jspec->mdsbase = vnames[0];
	jspec->auto_midas = 1; /* provisional! */
	jspec->midas_pd = midas_pd;
	for (i=0; i<ncols; i++) {
	    jspec->colnames[i] = NULL;
	    jspec->colnums[i] = 0;
	}
    } else {
	j = 1;
	for (i=0; i<ncols; i++) {
	    if (i > JOIN_TARG) {
		jspec->colnames[i] = vnames[j++];
	    } else {
		jspec->colnames[i] = NULL;
	    }
	    jspec->colnums[i] = 0;
	}
    }

    return 0;
}

static int expand_jspec (joinspec *jspec, int addvars)
{
    int i, ncols = jspec->ncols + addvars;
    char **colnames;
    int *colnums;

    colnames = realloc(jspec->colnames, ncols * sizeof *colnames);
    colnums = realloc(jspec->colnums, ncols * sizeof *colnums);

    if (colnames == NULL || colnums == NULL) {
	return E_ALLOC;
    }

    jspec->colnames = (const char **) colnames;
    jspec->colnums = colnums;

    for (i=jspec->ncols; i<ncols; i++) {
	jspec->colnames[i] = NULL;
	jspec->colnums[i] = 0;
    }

    jspec->ncols = ncols;

    return 0;
}

static void clear_jspec (joinspec *jspec, joiner *jr)
{
    free(jspec->colnames);
    free(jspec->colnums);

    if (jspec->timecols != NULL) {
	free(jspec->timecols);
    }

    if (jspec->c != NULL) {
	csvdata_free(jspec->c);
    } else if (jspec->dset != NULL) {
	destroy_dataset(jspec->dset);
    }

    if (jspec->wildnames != NULL) {
	strings_array_free(jspec->wildnames, jspec_n_vars(jspec));
    }

    if (jr != NULL && jspec->mdsnames != NULL) {
	strings_array_free(jspec->mdsnames, jr->midas_m);
    }

    if (jspec->tmpnames != NULL) {
	strings_array_free(jspec->tmpnames, jspec->n_tmp);
    }
}

static int *midas_revise_jspec (joinspec *jspec,
				DATASET *dset,
				int *n_add,
				int *err)
{
    DATASET *rdset = outer_dataset(jspec);
    int m, rpd = 0, nvars = 0;
    int *ret = NULL;

    if (jspec->midas_pd == 0) {
	/* outer pd not specified on input: so take
	   it from the discovered dataset
	*/
	rpd = jspec->midas_pd = rdset->pd;
    }

    m = midas_m_from_pd(dset, jspec->midas_pd);

    if (m == 0) {
	gretl_errmsg_sprintf("frequency %d in import data: \"spread\" will "
			     "not work", rpd);
	*err = E_PDWRONG;
	return NULL;
    } else {
	nvars = m;
    }

    ret = gretl_list_new(nvars);
    jspec->mdsnames = strings_array_new(nvars);

    if (ret == NULL || jspec->mdsnames == NULL) {
	*err = E_ALLOC;
    }

    if (!*err) {
	char *cname, tmp[VNAMELEN];
	int i, v, extlen;
	char mc;

	/* zero the count of added vars */
	*n_add = 0;

	/* create base for naming vars */
	mc = rpd == 12 ? 'm' : rpd == 4 ? 'q' : 'd';
	extlen = m < 10 ? 3 : 4;
	*tmp = '\0';
	strncat(tmp, jspec->mdsbase, VNAMELEN - 1);
	gretl_trunc(tmp, VNAMELEN - extlen - 1);

	for (i=0; i<nvars && !*err; i++) {
	    cname = gretl_strdup_printf("%s_%c%d", tmp, mc, nvars - i);
	    if (cname == NULL) {
		*err = E_ALLOC;
		break;
	    }
	    v = current_series_index(dset, cname);
	    ret[i+1] = v;
	    if (v < 0) {
		/* not a current series */
		if (gretl_type_from_name(cname, NULL)) {
		    *err = E_TYPES;
		} else {
		    *n_add += 1;
		}
	    }
	    if (*err) {
		free(cname);
	    } else {
		jspec->mdsnames[i] = cname;
	    }
	}
    }

    return ret;
}

static void maybe_transfer_string_table (DATASET *l_dset,
					 DATASET *r_dset,
					 joinspec *jspec,
					 int *targvars,
					 int orig_v)
{
    int i, lv, rv;

    for (i=1; i<=targvars[0]; i++) {
	lv = targvars[i];
	if (lv >= orig_v) {
	    /* it's a new series */
	    rv = outer_series_index(jspec, i);
	    if (rv > 0 && is_string_valued(r_dset, rv)) {
		/* let the new series grab the RHS string table */
		steal_string_table(l_dset, lv, r_dset, rv);
	    }
	}
    }
}

static int initial_midas_check (int nvars, int any_wild, int pd,
				DATASET *dset)
{
    int err;

    if (pd != 0 && pd != 12 && pd != 4 && pd != 5 && pd != 6 && pd != 7) {
	/* unacceptable outer data frequency */
	err = E_PDWRONG;
    } else if (nvars == 1 && (annual_data(dset) || quarterly_or_monthly(dset))) {
	/* might be OK, if no wildcard */
	err = any_wild ? E_DATA : 0;
    } else {
	err = E_DATA;
    }

    if (err) {
	gretl_errmsg_set(_("Invalid join specification"));
    }

    return err;
}

static int has_native_suffix (const char *fname)
{
    return has_suffix(fname, ".gdt") || has_suffix(fname, ".gdtb");
}

/**
 * gretl_join_data:
 * @fname: name of data file.
 * @vnames: name(s) of variables to create or modify.
 * @nvars: the number of elements in @vnames.
 * @dset: pointer to dataset.
 * @ikeyvars: list of 1 or 2 "inner" key variables, or NULL.
 * @okey: string specifying "outer" key(s) or NULL.
 * @filtstr: string specifying filter, or NULL.
 * @srcname: name of variable to import at source, or NULL.
 * @aggr: aggregation method specifier.
 * @seqval: 1-based sequence number for aggregation, or 0.
 * @auxname: name of auxiliary column for max or min aggregation,
 * or NULL.
 * @tconvstr: string specifying date columns for conversion, or NULL.
 * @tconvfmt: string giving format(s) for "timeconv" columns, or NULL.
 * @midas_pd: hint regarding pd for --aggr=spread (?).
 * @opt: may contain OPT_V for verbose operation, OPT_H to assume
 * no header row.
 * @prn: gretl printing struct (or NULL).
 *
 * Opens a delimited text data file or gdt file and carries out a
 * "join" operation to pull data into the current working dataset.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_join_data (const char *fname,
		     const char **vnames,
		     int nvars,
		     DATASET *dset,
		     const int *ikeyvars,
		     const char *okey,
		     const char *filtstr,
		     const char *srcname,
		     AggrType aggr,
		     int seqval,
		     const char *auxname,
		     const char *tconvstr,
		     const char *tconvfmt,
		     int midas_pd,
		     gretlopt opt,
		     PRN *prn)
{
    DATASET *outer_dset = NULL;
    joinspec jspec = {0};
    joiner *jr = NULL;
    jr_filter *filter = NULL;
    const char *varname;
    int okeyvars[3] = {0, 0, 0};
    char okeyname1[VNAMELEN] = {0};
    char okeyname2[VNAMELEN] = {0};
    char tkeyfmt[32] = {0};
    obskey auto_keys;
    int *targvars = NULL;
    int do_tsjoin = 0;
    int orig_v = dset->v;
    int add_v = 0;
    int modified = 0;
    int any_wild = 0;
    int verbose = (opt & OPT_V);
    int str_keys[2] = {0};
    int n_keys = 0;
    int err = 0;

    /** Step 0: preliminaries **/

    if (vnames == NULL || nvars < 1) {
	return E_DATA;
    }

    varname = vnames[0];

    targvars = get_series_indices(vnames, nvars, dset, &add_v,
				  &any_wild, &err);

    if (!err && srcname != NULL) {
	/* If we have a spec for the original name of a series
	   to import (@srcname), we cannot accept more than one
	   target name (in @vnames), nor can we accept any
	   wildcard specification -- unless we're doing a
	   "MIDAS join", in which case a single wildcard spec
	   is OK for the target series.
	*/
	if ((nvars > 1 || any_wild) && aggr != AGGR_MIDAS) {
	    gretl_errmsg_set(_("Invalid join specification"));
	    err = E_DATA;
	}
    }

    if (!err && aggr == AGGR_MIDAS) {
	err = initial_midas_check(nvars, any_wild, midas_pd, dset);
    }

    if (err) {
	return err;
    }

    err = set_up_jspec(&jspec, vnames, nvars, opt, any_wild,
		       aggr, midas_pd);
    if (err) {
	return err;
    }

    if (ikeyvars != NULL) {
	n_keys = ikeyvars[0];
    }

    obskey_init(&auto_keys);
    timeconv_map_init();

#if CDEBUG
    fputs("*** gretl_join_data:\n", stderr);
    fprintf(stderr, " filename = '%s'\n", fname);
    if (nvars > 1) {
	int i;

	fputs(" target series names:\n", stderr);
	for (i=0; i<nvars; i++) {
	    fprintf(stderr, "  '%s'\n", vnames[i]);
	}
    } else {
	fprintf(stderr, " target series name = '%s'\n", varname);
    }
    if (n_keys > 0) {
	fprintf(stderr, " inner key series ID = %d\n", ikeyvars[1]);
	if (n_keys == 2) {
	    fprintf(stderr, " second inner key series ID = %d\n", ikeyvars[2]);
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
    if (srcname != NULL) {
	fprintf(stderr, " source data series = '%s'\n", srcname);
    } else if (aggr != AGGR_COUNT) {
	fprintf(stderr, " source data series: assuming '%s' (from inner varname)\n",
		varname);
    }
    fprintf(stderr, " aggregation method = %d\n", aggr);
    if (auxname != NULL) {
	fprintf(stderr, " aggr auxiliary column = '%s'\n", auxname);
    }
    if (tconvstr != NULL) {
	fprintf(stderr, " timeconv = '%s'\n", tconvstr);
    }
    if (tconvfmt != NULL) {
	fprintf(stderr, " tconvfmt = '%s'\n", tconvfmt);
    }
#endif

    /* Step 1: process the arguments we got with regard to filtering
       and keys: extract the names of the columns that are required
       from the outer datafile, checking for errors as we go.
    */

    if (filtstr != NULL) {
	filter = make_join_filter(filtstr, &err);
    }

    if (!err && okey != NULL) {
	if (opt & OPT_K) {
	    /* cancel automatic MIDAS flag if present */
	    jspec.auto_midas = 0;
	    err = process_time_key(okey, okeyname1, tkeyfmt);
	} else {
	    err = process_outer_key(okey, n_keys, okeyname1, okeyname2, opt);
	}
    }

    /* Step 2: set up the array of required column names,
       jspec.colnames. This is an array of const *char, pointers
       to strings that "live" elsewhere. We leave any unneeded
       elements as NULL.
    */

    if (!err) {
	/* handle the primary outer key column, if any */
	if (*okeyname1 != '\0') {
	    jspec.colnames[JOIN_KEY] = okeyname1;
	} else if (n_keys > 0) {
	    jspec.colnames[JOIN_KEY] = dset->varname[ikeyvars[1]];
	}

	/* and the secondary outer key, if any */
	if (*okeyname2 != '\0') {
	    jspec.colnames[JOIN_KEY2] = okeyname2;
	} else if (n_keys > 1) {
	    jspec.colnames[JOIN_KEY2] = dset->varname[ikeyvars[2]];
	}

	/* the data or "payload" column */
	if (aggr != AGGR_COUNT) {
	    if (srcname != NULL) {
		jspec.colnames[JOIN_TARG] = srcname;
	    } else {
		jspec.colnames[JOIN_TARG] = varname;
	    }
	}

	/* the filter columns, if applicable */
	if (filter != NULL) {
	    jspec.colnames[JOIN_F1] = filter->vname1;
	    jspec.colnames[JOIN_F2] = filter->vname2;
	    jspec.colnames[JOIN_F3] = filter->vname3;
	}

	/* the auxiliary var for aggregation, if present */
	if (auxname != NULL) {
	    jspec.colnames[JOIN_AUX] = auxname;
	}
    }

    /* Step 3: handle the tconvert and tconv-fmt options */

    if (!err && tconvstr != NULL) {
	err = process_tconvert_info(&jspec, tconvstr, tconvfmt, tkeyfmt);
    }

    /* Step 4: read data from the outer file; check we got all the
       required columns; check that nothing is screwed up type-wise
    */

    if (!err) {
	PRN *vprn = verbose ? prn : NULL;

	if (has_native_suffix(fname)) {
	    gretlopt gdt_opt = OPT_NONE;

	    if (dataset_is_time_series(dset)) {
		/* import obs markers: may be needed */
		gdt_opt = OPT_M;
	    }
	    err = join_import_gdt(fname, &jspec, gdt_opt, vprn);
	} else {
	    err = join_import_csv(fname, &jspec, opt, vprn);
	}
	if (!err) {
	    outer_dset = outer_dataset(&jspec);
	}
#if CDEBUG > 2
	if (!err) {
	    print_outer_dataset(outer_dset, fname);
	}
#endif
	if (!err) {
	    err = check_for_missing_columns(&jspec);
	}
	if (!err) {
	    err = aggregation_type_check(&jspec, aggr);
	}
	if (!err) {
	    err = join_data_type_check(&jspec, dset, targvars, aggr);
	}
    }

    if (!err && verbose) {
	int i;

	pprintf(prn, _("Outer dataset: read %d columns and %d rows\n"),
		outer_dset->v - 1, outer_dset->n);
	for (i=1; i<outer_dset->v; i++) {
	    pprintf(prn, " col %d: '%s'\n", i, outer_dset->varname[i]);
	}
    }

    /* Step 5: set up keys and check for conformability errors */

    if (!err && jspec.colnames[JOIN_KEY] != NULL) {
	err = set_up_outer_keys(&jspec, dset, opt, ikeyvars, okeyvars,
				&auto_keys, str_keys);
    }

    if (!err && n_keys == 0 && dataset_is_time_series(dset)) {
	err = auto_keys_check(dset, outer_dset, opt, tkeyfmt,
			      &auto_keys, &n_keys, &do_tsjoin);
	if (do_tsjoin) {
	    goto transcribe;
	}
    }

    if (!err && n_keys == 0 && filter == NULL && aggr == 0 &&
	auto_keys.timefmt == NULL) {
	/* the simple case: no need to build joiner struct */
	err = join_simple_range_check(dset, outer_dset, targvars);
	goto transcribe;
    }

    /* Step 6: build the joiner struct from the outer dataset,
       applying a filter if one is specified
    */

    if (!err) {
	jr = build_joiner(&jspec, dset, filter, aggr, seqval,
			  &auto_keys, n_keys, &err);
	if (!err && jr == NULL) {
	    /* no matching data to join */
	    goto bailout;
	}
    }

    if (!err && filter != NULL && verbose) {
	pprintf(prn, "Filter: %d rows were selected\n", jr->n_rows);
    }

    /* Step 7: transcribe more info and sort the "joiner" struct */

    if (!err) {
	jr->n_keys = n_keys;
	jr->str_keys = str_keys;
	jr->l_keyno = ikeyvars;
	jr->r_keyno = okeyvars;
	if (jr->n_keys > 0) {
	    err = joiner_sort(jr);
	}
#if CDEBUG > 1
	if (!err) joiner_print(jr);
#endif
    }

    /* Step 8: another check now the joiner struct is ready */

    if (!err && jr != NULL && jr->n_keys == 0) {
	err = join_range_check(jr, dset, aggr);
    }

 transcribe:

    /* step 9: revise information on the series to be imported
       if we came across a wildcard specification or a case
       of MIDAS importation
    */

    if (!err && (jspec.wildcard || aggr == AGGR_MIDAS)) {
	free(targvars);
	if (aggr == AGGR_MIDAS) {
	    targvars = midas_revise_jspec(&jspec, dset, &add_v, &err);
	} else {
	    targvars = revise_series_indices(&jspec, dset, &add_v, &err);
	}
    }

    /* Step 10: transcribe or aggregate the data */

    if (!err && add_v > 0) {
	/* we need to add one or more new series on the left */
	if (jspec.wildnames != NULL) {
	    err = add_target_series((const char **) jspec.wildnames,
				    dset, targvars, add_v);
	} else if (jspec.mdsnames != NULL) {
	    err = add_target_series((const char **) jspec.mdsnames,
				    dset, targvars, add_v);
	} else {
	    err = add_target_series(vnames, dset, targvars, add_v);
	}
    }

    if (!err) {
	if (jr == NULL && do_tsjoin) {
	    ts_joiner tjr = {0};

	    fill_ts_joiner(dset, outer_dset, &tjr);
	    err = join_transcribe_multi_data(dset, outer_dset, targvars,
					     orig_v, &jspec, &tjr,
					     &modified);
	} else if (jr == NULL) {
	    err = join_transcribe_multi_data(dset, outer_dset, targvars,
					     orig_v, &jspec, NULL,
					     &modified);
	} else if (jr->n_keys == 0) {
	    err = join_transcribe_data(jr, targvars[1], add_v,
				       &jspec, &modified);
	} else {
	    err = aggregate_data(jr, ikeyvars, targvars, &jspec,
				 orig_v, &modified);
	    /* complete the job for MIDAS daily import */
	    if (!err && aggr == AGGR_MIDAS && midas_daily(jr)) {
		postprocess_daily_data(dset, targvars);
	    }
	}
    }

#if CDEBUG
    fprintf(stderr, "join: add_v = %d, modified = %d\n",
	    add_v, modified);
#endif

    if (!err && add_v > 0 && jspec.colnums[JOIN_TARG] > 0) {
	/* we added one or more new series */
	if (aggr != AGGR_MIDAS) {
	    maybe_transfer_string_table(dset, outer_dset, &jspec,
					targvars, orig_v);
	}
    }

    if (err) {
	dataset_drop_last_variables(dset, dset->v - orig_v);
    } else {
	if (add_v || modified) {
	    set_dataset_is_changed(dset, 1);
	}
	if (gretl_messages_on()) {
	    if (add_v) {
		pputs(prn, _("Data appended OK\n"));
	    } else if (modified) {
		pputs(prn, _("Data modified OK\n"));
	    } else {
		pputs(prn, _("No changes were made to the dataset\n"));
	    }
	}
    }

 bailout:

    /* null out file-scope "timeconv" globals */
    timeconv_map_destroy();

    if (auto_keys.timefmt != NULL) {
	free(auto_keys.timefmt);
    }

    clear_jspec(&jspec, jr);
    joiner_destroy(jr);
    jr_filter_destroy(filter);
    free(targvars);

    return err;
}
