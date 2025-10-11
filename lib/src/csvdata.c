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
#include "gretl_xml.h"
#include "gretl_midas.h"
#include "matrix_extra.h"
#include "gretl_www.h"
#include "gretl_join.h"
#include "join_priv.h"
#include "csvdata.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <errno.h>

#define CDEBUG 0    /* CSV reading in general */
#define COMPDEBUG 0 /* line compression */

#define CSVSTRLEN 128

typedef enum {
    CSV_HAVEDATA = 1 << 0,
    CSV_GOTDELIM = 1 << 1,
    CSV_GOTTAB   = 1 << 2,
    CSV_GOTSEMI  = 1 << 3,
    CSV_GOTPIPE  = 1 << 4,
    CSV_BLANK1   = 1 << 5,
    CSV_OBS1     = 1 << 6,
    CSV_TRAIL    = 1 << 7,
    CSV_AUTONAME = 1 << 8,
    CSV_REVERSED = 1 << 9,
    CSV_DOTSUB   = 1 << 10,
    CSV_ALLCOLS  = 1 << 11,
    CSV_BOM      = 1 << 12,
    CSV_VERBOSE  = 1 << 13,
    CSV_THOUSEP  = 1 << 14,
    CSV_NOHEADER = 1 << 15,
    CSV_QUOTES   = 1 << 16,
    CSV_AS_MAT   = 1 << 17
} csvflags;

struct csvprobe_ {
    DATASET *dset; /* more info might be wanted */
};

typedef struct csvprobe_ csvprobe;

struct csvdata_ {
    csvflags flags;
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
#define csv_got_pipe(c)           (c->flags & CSV_GOTPIPE)
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
#define csv_set_got_pipe(c)         (c->flags |= CSV_GOTPIPE)
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

static int
time_series_label_check (DATASET *dset, int reversed, char *skipstr,
                         int convert_pd, PRN *prn);

/* file-scope global */
static char import_na[8];


/* for use in gretl_join.c */

DATASET *csvdata_get_dataset (csvdata *c)
{
    return c->dset;
}

/* shared with gretl_join.c */

void csvdata_free (csvdata *c)
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
        fprintf(stderr, "\ncsvdata_new: c->delim = '%c', c->decpoint = '%c'\n",
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
            tc = calendar_obs_number(dset->S[t], dset, 0) - offset;
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
        tc = calendar_obs_number(dset->S[t], dset, 0) - misscount;
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
            int Tc = calendar_obs_number(lbl2, dset, 0) + 1;
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

        if (dset->pd == 5 && (wd == 6 || wd == 7)) {
            /* Got Sat or Sun, can't be 5-day daily? */
            alt_pd = wd;
            pprintf(prn, _("Found a Saturday (%s): re-trying with pd = %d\n"),
                    dset->S[s], alt_pd);
            break;
        } else if (dset->pd == 6 && wd == 7) {
            /* Got Sun, can't be 6-day daily? */
            alt_pd = wd;
            pprintf(prn, _("Found a Sunday (%s): re-trying with pd = %d\n"),
                    dset->S[s], alt_pd);
            break;
        }

        n = calendar_obs_number(dset->S[s], dset, 0);
        if (n < t) {
            pprintf(prn, _("Daily dates error at t = %d:\n"
                    "  calendar_obs_number() for '%s' = %d but t = %d\n"),
                    t, dset->S[s], n, t);
            err = 1;
        } else if (n > fulln - 1) {
            pprintf(prn, _("Error: date '%s' out of bounds\n"), dset->S[s]);
            err = 1;
        } else if (nbak > 0 && n == nbak) {
            pprintf(prn, _("Error: date '%s' is repeated\n"), dset->S[s]);
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
            pprintf(prn, _("   \"%s\": quarterly date with spurious zero?\n"), label);
            *extra_zero = 1;
            *ppd = pd0 = pd = 4;
            goto restart;
        } else if (yr != Ey || per != Ep) {
            ret = 0;
        }

        if (!ret) {
            pprintf(prn, "   %s: %s\n", label, _("not a consistent date"));
            break;
        }
    }

    if (ret) {
        if (pmin == 3) {
            pprintf(prn, "   \"%s\": %s\n", bad, _("quarterly data pretending to be monthly?"));
            *ppd = 4;
        } else if (pd == pd0 + 1) {
            pprintf(prn, "   \"%s\": %s\n", bad, _("BLS-type nonsense? Trying again"));
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
                    pprintf(prn, _("quarter %d: not possible\n"), p);
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
                        pprintf(prn, _("quarter %d: not possible\n"), p);
                    }
                } else {
                    if (p > 0 && p < 13) {
                        pprintf(prn, _("month %s?\n"), s);
                        pd = 12;
                    } else {
                        pprintf(prn, _("month %d: not possible\n"), p);
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
    char *lbl1, *lbl2;
    int len1, len2;
    int pd = -1;
    int n = dset->n;

    if (dset->S == NULL) {
	return pd;
    }

    lbl1 = dset->S[0];
    lbl2 = dset->S[n - 1];
    len1 = strlen(lbl1);
    len2 = strlen(lbl2);

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

/* If we got a UTF-16 or UTF-32 BOM, try recoding to UTF-8 before
   parsing data. We write the recoded text to a temporary file in the
   user's "dotdir" (and then delete that file once we're done).
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

/* Check the first 4 bytes of "CSV" input for a Byte Order Mark. If we
   find the UTF-8 BOM (typically written by Microsoft tools), simply
   record the fact so that we can skip it on reading. But if we find a
   BOM indicating a 16-bit or 32-bit unicode encoding, flag this by
   returning a non-zero @ucode value; in that case we'll attempt a full
   recording of the input (via GLib) before we start reading data.
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

static int report_quotation_info (int *qtotal, int *nqmax,
				  char qchar, PRN *prn)
{
    int err = 0;

    if (qtotal[0] > 0) {
	pprintf(prn, _("Found %d double-quotes, max %d per line\n"),
		qtotal[0], nqmax[0]);
    }
    if (qtotal[1] > 0) {
	pprintf(prn, _("Found %d single-quotes, max %d per line\n"),
		qtotal[1], nqmax[1]);
    }
    if (qchar == '"') {
	pputs(prn, _("Assuming double-quote is the relevant "
		     "quotation character\n"));
    } else if (qchar == '\'') {
	pputs(prn, _("Assuming single-quote is the relevant "
		     "quotation character\n"));
    } else {
        gretl_warnmsg_set(_("unmatched quotation marks in CSV file"));
	pputs(prn, _("Quotation broken: marks not paired\n"));
    }

    return err;
}

static int handle_CR (int *crlf, gzFile fp)
{
    int c, c1 = gzgetc(fp);

    if (c1 == EOF) {
        return c1;
    } else if (c1 == 0x0a) {
        /* CR + LF -> LF */
        *crlf = 1;
        c = c1;
    } else {
        /* old Mac-style: CR not followed by LF */
        c = 0x0a;
        gzungetc(c1, fp);
    }

    return c;
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

   We return the maximum line length, or -1 on error.
*/

static int csv_max_line_length (gzFile fp, csvdata *cdata, PRN *prn)
{
    int c, c1, cbak = 0, cc = 0;
    int comment = 0, maxlinelen = 0;
    int qcand[2] = {1, 1};
    int nq[2] = {0};
    int nqmax[2] = {0};
    int qtotal[2] = {0};
    int dquoted = 0;
    int crlf = 0, lines = 0;
    int i, err = 0;

    csv_set_trailing_comma(cdata); /* just provisionally */

    while ((c = gzgetc(fp)) != EOF) {
        if (c == 0x0d) {
            c = handle_CR(&crlf, fp);
            if (c == EOF) {
                break;
            }
        }
        if (c == 0x0a) {
	    /* we reached the end of a line */
            if (cc > maxlinelen) {
                maxlinelen = cc;
            }
            cc = 0;
            if (cbak != 0 && cbak != ',') {
                csv_unset_trailing_comma(cdata);
            }
            lines++;
	    for (i=0; i<2; i++) {
                /* quotation--mark info, double and single */
		if (nq[i] % 2) {
                    /* unpaired: not a valid candidate */
		    qcand[i] = 0;
		}
		if (nq[i] > nqmax[i]) {
		    nqmax[i] = nq[i];
		}
		qtotal[i] += nq[i];
		nq[i] = 0; /* reset count */
	    }
            dquoted = 0;
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
	    /* line starts with hash mark */
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
            } else if (c == '|') {
		csv_set_got_pipe(cdata);
	    }
            if (c == cdata->delim && !dquoted) {
                csv_set_got_delim(cdata);
            } else if (c == '"') {
		dquoted = !dquoted;
		nq[0] += 1;
            } else if (c == '\'') {
		nq[1] += 1;
            }
        }
        cc++;
    }

    if (maxlinelen == 0) {
        pputs(prn, _("Data file is empty\n"));
    } else if (csv_has_trailing_comma(cdata)) {
        pputs(prn, _("Data file has trailing commas\n"));
    }

    for (i=0; i<2; i++) {
	if (qcand[i] && qtotal[i] == 0) {
	    qcand[i] = 0;
	}
    }

    if (qtotal[0] || qtotal[1]) {
        if (qcand[0] && qcand[1]) {
            /* quotes: double and single are both paired,
               maybe prefer the more numerous as @qchar?
            */
            cdata->qchar = qtotal[1] > qtotal[0] ? '\'' : '"';
        } else if (qcand[0]) {
	    cdata->qchar = '"';
	} else if (qcand[1]) {
	    cdata->qchar = '\'';
	}
	err = report_quotation_info(qtotal, nqmax, cdata->qchar, prn);
    }

    if (!err && maxlinelen > 0) {
        /* allow for newline and null terminator */
        maxlinelen += 2 + crlf;
    }

    return err ? -1 : maxlinelen;
}

#define nonspace_delim(d) (d != ',' && d != ';' && d != '\t' && d != '|')

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

static void purge_quoted_delimiters (csvdata *c)
{
    char *s = c->line;
    int inquote = 0;

    while (*s) {
        if (*s == '"') {
            inquote = !inquote;
        } else if (inquote && *s == c->delim) {
            *s = (c->delim == ' ')? '_' : ' ';
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

#if COMPDEBUG
    fprintf(stderr, "HERE compress: keep_quotes=%d, qchar=%c\n",
            csv_keep_quotes(c) ? 1 : 0, c->qchar ? c->qchar : '0');
    fprintf(stderr, " before:\n  '%s'\n", c->line);
#endif

    if (!csv_keep_quotes(c)) {
        purge_quoted_delimiters(c);
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

#if COMPDEBUG
    fprintf(stderr, " after:\n  '%s'\n", c->line);
#endif


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
        char field1[VNAMELEN];
	int quoted = 0;
        int i = 0;

        if (c->delim == ' ' && *s == ' ') {
            s++;
        }

	if (*s == c->qchar) {
	    quoted = 1;
	}

        while (*s && i < sizeof field1) {
            if (!quoted && *s == c->delim) {
                break;
            } else if (!quoted && *s == '\t') {
                /* presence of a tab must indicate tab-separation? */
                c->delim = '\t';
                goto tryagain;
            } else if (i > 0 && *s == c->qchar) {
		quoted = 0;
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

/* Returns 1 if the string @s should be counted as representing a
   missing value, 0 otherwise. If there is a user-set "csv_read_na"
   value we consult it, otherwise we consult a set of default values.
*/

int import_na_string (const char *s)
{
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

    if (*import_na != '\0' && strcmp(import_na, "default")) {
        /* the user has set a specific "NA" string */
	return !strcmp(s, import_na);
    }

    /* Should we continue in the following way if "csv_read_na"
       is user-specified? */
    for (i=0; defaults[i] != NULL; i++) {
	if (!strcmp(s, defaults[i])) {
	    return 1;
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

/* In the case where we think we've found thousands separators in
   numerical input, provisionally mark all "non-numeric" values as NAs;
   we do this prior to a second pass through the data.
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
    struct tm t = {0,0,0,1,0,0,0,0,-1};
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

/* Here we check whether we get a consistent reading on the number of
   fields per line in the CSV file.
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
        } else if (chkcols > c->ncols) {
            pprintf(prn, _("   ...but row %d has %d fields: aborting\n"),
                    c->nrows, chkcols);
            err = E_DATA;
        } else if (chkcols < c->ncols) {
            pprintf(prn, _("Warning: row %d has %d fields: %d entries set to missing\n"),
                    c->nrows, chkcols, c->ncols - chkcols);
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

    if (!strcmp(s, ".")) {
	/* treat as NA? */
	ret = 1;
    } else if (c->decpoint == '.') {
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
        if (gretl_reserved_word(vname)) {
            /* try a fix for this */
            int n = strlen(vname);

            if (n < VNAMELEN-1) {
                strcat(vname, "_");
            } else {
                vname[n-1] = '_';
            }
        }
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

#define starts_number(c) (isdigit((unsigned char) c) || c == '-' ||     \
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
	int quoted = (*p == c->qchar);

        i = 0;
        while (*p && (*p != c->delim || quoted)) {
            if (i < CSVSTRLEN - 1) {
                c->str[i++] = *p;
            }
	    if (i > 1 && *p == c->qchar) {
		quoted = 0;
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
        gretl_errmsg_set(_("No relevant columns were found"));
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
        pprintf(prn, _("warning: %d strings were truncated.\n"), truncated);
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
            } else if (c->nrows == 1) {
		csv_set_no_header(c);
		c->dset->n = 1;
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
 * the current work space. Shared with gretl_join.c.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int real_import_csv (const char *fname,
		     DATASET *dset,
		     const char *cols,
		     const char *rows,
		     joinspec *join,
		     void *probe,
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
	} else if (csv_got_pipe(c)) {
	    c->delim = '|';
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
        } else if (c->delim != '|' && csv_got_pipe(c)) {
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

    if (c->dset->v == 1 || c->dset->n == 0) {
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
            pprintf(prn, _("WARNING: Check your data! gretl has stripped out "
                    "what appear to be\nextraneous lines in a %s dataset: "
                    "this may not be right.\n\n"),
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
            gretl_string_table_finalize(c->st, c->dset);
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
            gretl_errmsg_set(_("Couldn't find all variable names"));
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

/* Read the first @n_lines lines of @fname and write the
   content to @prn.
*/

int peek_at_csv (const char *fname, int n_lines, PRN *prn)
{
    csvdata *c = NULL;
    gzFile fp = NULL;
    gchar *altname = NULL;
    int recode;
    int err = 0;

    fp = gretl_gzopen(fname, "rb");
    if (fp == NULL) {
        return E_FOPEN;
    }

    c = csvdata_new(NULL);
    if (c == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        recode = csv_unicode_check(fp, c, prn);
        if (recode) {
            err = csv_recode_input(&fp, fname, &altname, recode, prn);
        }
    }

    if (!err) {
        c->maxlinelen = csv_max_line_length(fp, c, prn);
        if (c->maxlinelen <= 0) {
            err = E_DATA;
        }
    }

    if (!err) {
        c->line = malloc(c->maxlinelen);
        if (c->line == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        char *s;
        int i;

        gzrewind(fp);
        gretl_print_reset_buffer(prn);

        for (i=0; i<n_lines; i++) {
            s = gzgets(fp, c->line, c->maxlinelen);
            if (s == Z_NULL) {
                break;
            } else {
                pputs(prn, s);
            }
        }
    }

    gzclose(fp);
    csvdata_free(c);

    if (altname != NULL) {
        gretl_remove(altname);
        g_free(altname);
    }

    return err;
}
