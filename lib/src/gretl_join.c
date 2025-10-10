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
#include "gretl_xml.h"
#include "gretl_midas.h"
#include "uservar.h"
#include "csvdata.h"
#include "join_priv.h"
#include "gretl_join.h"

#ifdef WIN32
# include "gretl_win32.h" /* for strptime() */
#endif

#define AGGDEBUG 0  /* aggregation in "join" */
#define TDEBUG 0    /* handling of time keys in "join" */
#define JDEBUG 0    /* joining in general */

enum {
    JOIN_KEY,
    JOIN_F1,
    JOIN_F2,
    JOIN_F3,
    JOIN_KEY2,
    JOIN_AUX,
    JOIN_TARG
};

enum {
    TCONV_FMT = 0,
    TKEY_FMT = 1
};

typedef double keynum;

struct jr_row_ {
    int n_keys;     /* number of keys (needed for qsort callback) */
    keynum keyval;  /* primary key value */
    keynum keyval2; /* secondary key value, if present */
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
    keynum *keys;   /* array of unique (primary) key values as doubles */
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

struct time_mapper {
    int ncols;         /* number of "timeconv" columns */
    char **colnames;   /* array of outer-dataset column names */
    char *tname;       /* the name of the "tkey", if among colnames, or NULL */
    char **fmt;        /* array of up to two time-format strings, or NULL */
    char m_means_q[2]; /* array of "monthly means quarterly" flags */
};

struct jr_matcher_ {
    keynum *k1;    /* first key value, per observation */
    keynum *k2;    /* second key value per obs, or NULL */
    int *pos;      /* position of match in outer key array */
};

typedef struct jr_matcher_ jr_matcher;

#define KEYMISS -999

#define is_wildstr(s) (strchr(s, '*') || strchr(s, '?'))

/* file-scope global */
struct time_mapper tconv_map;

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
                gretl_errmsg_sprintf(_("'%s' is not a valid date"), s);
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
                gretl_errmsg_sprintf(_("'%s' is not a valid date"), s);
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

static void timeconv_map_set (int ncols, char **colnames,
                              char *tname, char **fmt)
{
    tconv_map.ncols = ncols;
    tconv_map.colnames = colnames;
    tconv_map.tname = tname;
    tconv_map.fmt = fmt;

#if TDEBUG
    if (fmt != NULL) {
	fprintf(stderr, "timeconv_map_set: tname '%s' tconv_fmt '%s' tkey_fmt '%s'\n",
	    tname, fmt[TCONV_FMT], fmt[TKEY_FMT]);
    }
#endif

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
    struct tm t = {0,0,0,1,0,0,0,0,-1};
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
           occurred, while a non-NULL and non-empty return string
           means a trailing portion of the input was not
           processed.
        */
#if TDEBUG
	if (i < 6) {
	    fprintf(stderr, "i=%d, src=%d, right-hand pd=%d, tmft='%s', s='%s'\n",
		    i, s_src, jr->r_dset->pd, tfmt, s);
	}
#endif
        test = strptime(s, tfmt, &t);
	if (test == NULL && s_src == 3 && j == 0) {
	    if (strchr(tfmt, '-') && strlen(s) == 4 &&
		jr->r_dset->pd == 1 && integer_string(s)) {
		/* annual data from CSV? */
		set_time_format(jr->auto_keys, "%Y");
		goto finish;
	    }
	}
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
            gretl_errmsg_sprintf(_("'%s' does not match the format '%s'"), s, tfmt);
            fprintf(stderr, "time-format match error in read_outer_auto_keys:\n"
                    " remainder = '%s' (source = %s)\n", test ? test : "null",
                    s_src < 3 ? "specified time column" : "first-column strings");
        }
    }

 finish:

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
            gretl_errmsg_sprintf(_("'%.8g' is not a valid date"), x);
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

#if JDEBUG > 1
        fprintf(stderr, "filter genr: '%s':\n", line);
        for (i=0; i<r_dset->n; i++) {
            fprintf(stderr, " %d: %g\n", i+1, filter->val[i]);
        }
#endif
        for (i=0; i<r_dset->n; i++) {
            if (na(filter->val[i])) {
                gretl_errmsg_sprintf(_("join filter: indeterminate "
                                     "value on row %d"), i+1);
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
            gretl_errmsg_sprintf(_("%s: invalid secondary outer key value on row %d"),
                                 "join", row);
        } else {
            gretl_errmsg_sprintf(_("%s: invalid (primary) outer key value on row %d"),
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

#if JDEBUG > 2
    fprintf(stderr, "join filter: %s row %d\n",
            ret ? "keeping" : "discarding", i);
#endif

    return ret;
}

static DATASET *outer_dataset (joinspec *jspec)
{
    if (jspec->c != NULL) {
        return csvdata_get_dataset(jspec->c);
    } else {
        return jspec->dset;
    }
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

#if JDEBUG
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

#if JDEBUG
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

#if JDEBUG
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
#if JDEBUG > 1
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

#if JDEBUG > 1

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

# if JDEBUG > 2
static void print_outer_dataset (const DATASET *dset, const char *fname)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

    pprintf(prn, "Data extracted from %s:\n", fname);
    printdata(NULL, NULL, (DATASET *) dset, OPT_O, prn);
    gretl_print_destroy(prn);
}
# endif

#endif /* JDEBUG */

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

#define KEYMISS -999

static keynum matcher_get_k2 (jr_matcher *matcher, int i)
{
    return matcher->k2 == NULL ? 0 : matcher->k2[i];
}

static void matcher_set_k2 (jr_matcher *matcher, int i,
			    keynum val)
{
    if (matcher->k2 != NULL) {
	matcher->k2[i] = val;
    }
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
			  jr_matcher *matcher,
			  int s, int v,
                          int revseq,
                          double *xmatch,
                          double *auxmatch,
                          int *nomatch,
                          int *err)
{
    keynum key1 = matcher->k1[s];
    keynum key2 = matcher_get_k2(matcher, s);
    int pos = matcher->pos[s];
    double x, xa;
    int imin, imax;
    int i, n, ntotal;

#if AGGDEBUG
    fprintf(stderr, " key1 = %g, key2 = %g\n", key1, key2);
    fprintf(stderr, " key1 matched at position %d\n", pos);
#endif

    /* how many matches at @pos? */
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
    } /* end of special MIDAS case */

    /* We now fill out the array @xmatch with non-missing values
       from the matching outer rows. If we have a secondary key
       we screen for matches on that as we go.
    */

    n = 0;      /* will hold count of non-NA matches */
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
        /* handle case of no match on secondary key */
        if (ntotal == 0 && jr->aggr != AGGR_COUNT) {
            *nomatch = 1;
            return NADBL;
        }
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

static void handle_midas_setup (joiner *jr, int i, int lv, int rv,
                                int revseq)
{
    char label[MAXLABEL];

    series_set_midas_period(jr->l_dset, lv, revseq);
    sprintf(label, "%s in sub-period %d",
            jr->r_dset->varname[rv], revseq);
    series_record_label(jr->l_dset, lv, label);
    series_set_midas_freq(jr->l_dset, lv, jr->r_dset->pd);
    if (i == 1) {
        series_set_midas_anchor(jr->l_dset, lv);
    }
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

static void jr_matcher_free (jr_matcher *matcher)
{
    free(matcher->k1);
    free(matcher->k2);
    free(matcher->pos);
}

static int jr_matcher_init (jr_matcher *matcher, int nobs,
                            int nkeys)
{
    int i, err = 0;

    matcher->k1 = malloc(nobs * sizeof *matcher->k1);
    matcher->pos = malloc(nobs * sizeof *matcher->pos);
    matcher->k2 = NULL;

    if (matcher->k1 == NULL || matcher->pos == NULL) {
        err = E_ALLOC;
    } else if (nkeys == 2) {
        matcher->k2 = malloc(nobs * sizeof *matcher->k2);
        if (matcher->k2 == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        jr_matcher_free(matcher);
        return err;
    }

    for (i=0; i<nobs; i++) {
        matcher->k1[i] = 0;
        matcher->pos[i] = 0;
        if (matcher->k2 != NULL) {
            matcher->k2[i] = 0;
        }
    }

    return 0;
}

static int get_all_inner_key_values (joiner *jr,
                                     const int *ikeyvars,
                                     jr_matcher *matcher)
{
    DATASET *dset = jr->l_dset;
    int i, j, n = sample_size(dset);
    int err;

    err = jr_matcher_init(matcher, n, jr->n_keys);
    if (err) {
        return err;
    }

    for (i=dset->t1, j=0; i<=dset->t2 && !err; i++, j++) {
        if (using_auto_keys(jr)) {
            /* using the LHS dataset obs info */
            char obs[OBSLEN];

            ntolabel(obs, i, dset);
            if (calendar_data(dset)) {
                guint32 ed = get_epoch_day(obs);

                if (jr->n_keys == 2) {
                    /* inner daily, outer monthly */
                    int y, m, d;

                    ymd_bits_from_epoch_day(ed, &y, &m, &d);
                    matcher->k1[j] = y;
                    matcher->k2[j] = m;
                } else {
                    matcher->k1[j] = ed;
                }
            } else {
                /* monthly or quarterly (FIXME any others?) */
                matcher->k1[j] = atoi(obs);
		matcher_set_k2(matcher, j, atoi(obs + 5));
            }
        } else {
            /* using regular LHS key series */
            double dk1, dk2 = 0;
            keynum k1 = 0, k2 = 0;

            dk1 = dset->Z[ikeyvars[1]][i];
            if (jr->n_keys == 2) {
                dk2 = dset->Z[ikeyvars[2]][i];
            }
            if (na(dk1) || na(dk2)) {
                matcher->pos[j] = KEYMISS;
            } else {
                k1 = dtoll(dk1, &err);
                if (!err && jr->n_keys == 2) {
                    k2 = dtoll(dk2, &err);
                }
            }
            if (!err && matcher->pos[j] != KEYMISS) {
                matcher->k1[j] = k1;
		matcher_set_k2(matcher, j, k2);
            }
        }

        if (!err && matcher->pos[j] != KEYMISS) {
            /* look up position in outer keys array */
            matcher->pos[j] = binsearch(matcher->k1[j], jr->keys,
                                        jr->n_unique, 0);
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
    jr_matcher matcher = {0};
    series_table *rst = NULL;
    series_table *lst = NULL;
    DATASET *dset = jr->l_dset;
    double *xmatch = NULL;
    double *auxmatch = NULL;
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

    err = get_all_inner_key_values(jr, ikeyvars, &matcher);

    for (i=1; i<=targvars[0] && !err; i++) {
        /* loop across the series to be added/modified */
        int s, rv, lv = targvars[i];
        int strcheck = 0;

#if AGGDEBUG
        fprintf(stderr, "\nworking on series %d\n", i);
#endif

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

        for (t=dset->t1, s=0; t<=dset->t2 && !err; t++, s++) {
	    int nomatch = 0;
            double zt;

#if AGGDEBUG
            fprintf(stderr, " working on LHS obs %d (v=%d, value %g), s=%d\n",
                    t, lv, dset->Z[lv][t], s);
#endif
            if (matcher.pos[s] == KEYMISS) {
                dset->Z[lv][t] = NADBL;
                continue;
            } else if (matcher.pos[s] < 0) {
		nomatch = 1;
		zt = (jr->aggr == AGGR_COUNT)? 0 : NADBL;
                continue;
	    } else {
		zt = aggr_value(jr, &matcher, s, rv, revseq, xmatch,
				auxmatch, &nomatch, &err);
	    }
#if AGGDEBUG
            if (na(zt)) {
                fprintf(stderr, " aggregate_data: got NA (keys=%g,%g, err=%d)\n",
                        matcher.k1[s], matcher_get_k2(&matcher, s), err);
            } else {
                fprintf(stderr, " aggregate_data: got %.12g (keys=%g,%g, err=%d)\n",
                        zt, matcher.k1[s], matcher_get_k2(&matcher, s), err);
            }
#endif
            if (!err && strcheck && !na(zt)) {
                zt = maybe_adjust_string_code(rst, lst, zt, &err);
            }
            if (!err) {
                if (lv >= orig_v) {
                    /* @lv is a newly added series */
                    dset->Z[lv][t] = zt;
                } else if (zt != dset->Z[lv][t]) {
                    if (nomatch && !na(dset->Z[lv][t])) {
                        ; /* leave existing data alone */
                    } else {
                        dset->Z[lv][t] = zt;
                        *modified += 1;
                    }
                }
            }
        }

        if (!err && jr->aggr == AGGR_MIDAS) {
            handle_midas_setup(jr, i, lv, rv, revseq);
        }

        revseq--;
    }

    free(xmatch);
    jr_matcher_free(&matcher);

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
            gretl_errmsg_sprintf(_("join: '%s' not matched"), l_dset->varname[lv]);
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
        gretl_errmsg_set(_("join: the observation ranges don't match"));
        err = E_DATA;
    } else if (r_dset->v - 1 < targlist[0]) {
        gretl_errmsg_set(_("join: series missing on the right"));
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

#if JDEBUG

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

#if JDEBUG
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

static int check_for_quarterly_format (obskey *auto_keys, int pd)
{
    char *s = auto_keys->timefmt;
    int i, err = 0;

#if TDEBUG
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
                gretl_errmsg_sprintf(_("The '%c' format is not applicable "
                                     "for data with frequency %d"),
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

#if JDEBUG
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
            gretl_errmsg_sprintf(_("'%s' is a string variable: aggregation type "
                                 "is not applicable"), dset->varname[aggcol]);
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
#if JDEBUG
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
                pprintf(prn, _("checking for '%s'\n"), S[i]);
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
                    gretl_errmsg_sprintf(_("'%s': not found"), S[i]);
                }
            }
            if (prn != NULL) {
                pprintf(prn, _(" found %d match(es)\n"), match);
            }
        }

        if (!err && (vlist == NULL || vlist[0] == 0)) {
            gretl_errmsg_set(_("No matching data were found"));
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

/* add/replace an entry in @jspec's colnames array, and record the
   fact that @src is now owned by @jspec so we can free it when we're
   finished and hence avoid leaking memory
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
			    DATASET *ldset,
                            PRN *prn)
{
    int err = 0;

    if (jspec->wildcard) {
        err = determine_csv_matches(fname, jspec, prn);
        if (err) {
            pputs(prn, _("join_import_csv: failed in matching varnames\n"));
        }
    }

    if (!err) {
        err = real_import_csv(fname, NULL, NULL, NULL, jspec,
                              NULL, NULL, OPT_NONE, prn);
    }

    if (!err && !jspec->user_tkey && dataset_is_time_series(ldset)) {
	/* If we have time series data on the left, check on
	   the imported CSV side for a clear time series
	   interpretation -- unless the caller has specified an
	   outer key.
	*/
	DATASET *rdset = csvdata_get_dataset(jspec->c);

	if (rdset->S != NULL) {
	    int pd, reversed = 0;

	    pd = test_markers_for_dates(rdset, &reversed, NULL, prn);
	    if (pd > 0 && !reversed) {
		rdset->structure = TIME_SERIES;
	    }
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
        gretl_errmsg_sprintf(_("frequency %d in import data: \"spread\" will "
                             "not work"), rpd);
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

/* For newly added series: import labels and/or string tables
   from the outer dataset, if applicable.
*/

static void maybe_import_strings (DATASET *l_dset,
				  DATASET *r_dset,
				  joinspec *jspec,
				  int *targvars,
				  int orig_v)
{
    const char *label;
    int i, lv, rv;

    for (i=1; i<=targvars[0]; i++) {
	if (targvars[i] < orig_v) {
	    /* not a new series, skip it */
	    continue;
	}
	rv = outer_series_index(jspec, i);
	if (rv > 0) {
	    lv = targvars[i];
	    label = series_get_label(r_dset, rv);
	    if (label != NULL && *label != '\0') {
		/* transcribe descriptive label */
		series_set_label(l_dset, lv, label);
	    }
	    if (is_string_valued(r_dset, rv)) {
		/* grab the RHS string table */
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

    if (okey != NULL || tconvstr != NULL || tconvfmt != NULL) {
	jspec.user_tkey = 1;
    }

#if JDEBUG
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
            err = join_import_csv(fname, &jspec, dset, vprn);
        }
        if (!err) {
            outer_dset = outer_dataset(&jspec);
        }
#if JDEBUG > 2
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
        pprintf(prn, _("Filter: %d rows were selected\n"), jr->n_rows);
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
#if JDEBUG > 1
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

#if JDEBUG
    fprintf(stderr, "join: add_v = %d, modified = %d\n",
            add_v, modified);
#endif

    if (!err && add_v > 0 && jspec.colnums[JOIN_TARG] > 0) {
        /* we added one or more new series */
        if (aggr != AGGR_MIDAS) {
            maybe_import_strings(dset, outer_dset, &jspec,
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

/* called in csvdata.c */

int timecol_get_format (const DATASET *dset, int v,
			char **pfmt, int *q)
{
    if (tconv_map.fmt == NULL) {
	/* no formats present */
        return 0;
    } else if (tconv_map.tname == NULL) {
        /* get the common "tconvert" format */
        *pfmt = tconv_map.fmt[TCONV_FMT];
        *q = tconv_map.m_means_q[TCONV_FMT];
        return 1;
    } else if (!strcmp(dset->varname[v], tconv_map.tname)) {
        /* get the tkey-specific format */
        *pfmt = tconv_map.fmt[TKEY_FMT];
        *q = tconv_map.m_means_q[TKEY_FMT];
        return 1;
    } else if (tconv_map.fmt[TCONV_FMT] != NULL) {
        /* get the other one */
        *pfmt = tconv_map.fmt[TCONV_FMT];
        *q = tconv_map.m_means_q[TCONV_FMT];
        return 1;
    }

    return 0;
}
