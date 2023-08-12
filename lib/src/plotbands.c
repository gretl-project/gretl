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

/* plotbands.c for gretl */

#include "libgretl.h"
#include "usermat.h"
#include "uservar.h"
#include "plot_priv.h"

typedef enum {
    BAND_LINE,
    BAND_FILL,
    BAND_DASH,
    BAND_BARS,
    BAND_STEP
} BandStyle;

typedef struct {
    int center;      /* ID of center-of-band series */
    int width;       /* ID of width-of-band series */
    double factor;   /* factor by which to multiply width */
    int bdummy;      /* flag for drawing rectangles */
    BandStyle style; /* see enumeration above */
    char rgb[10];    /* specific color, if wanted */
} band_info;

static band_info *band_info_new (void)
{
    band_info *bi = malloc(sizeof *bi);

    if (bi != NULL) {
	bi->center = -1;
	bi->width = -1;
	bi->factor = 1.0;
	bi->bdummy = 0;
	bi->style = BAND_LINE;
	bi->rgb[0] = '\0';
    }

    return bi;
}

static void print_pm_filledcurve_line (band_info *bi,
                                       const char *title,
                                       FILE *fp)
{
    char cstr[10];

    if (bi->rgb[0] != '\0') {
        *cstr = '\0';
        strncat(cstr, bi->rgb, 9);
    } else {
        print_rgb_hash(cstr, get_shadecolor());
    }

    if (title == NULL) {
        fprintf(fp, "'-' using 1:($2-%g*$3):($2+%g*$3) "
                "notitle lc rgb \"%s\" w filledcurve, \\\n",
                bi->factor, bi->factor, cstr);
    } else {
        fprintf(fp, "'-' using 1:($2-%g*$3):($2+%g*$3) "
                "title '%s' lc rgb \"%s\" w filledcurve, \\\n",
                bi->factor, bi->factor, title, cstr);
    }
}

static void print_user_pm_data (const double *x,
                                const double *c,
                                const double *w,
                                int t1, int t2,
                                FILE *fp)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (na(c[t]) || na(w[t])) {
            fprintf(fp, "%.10g %s %s\n", x[t], GPNA, GPNA);
        } else {
            fprintf(fp, "%.10g %.10g %.10g\n", x[t], c[t], w[t]);
        }
    }

    fputs("e\n", fp);
}

/* Handle the special case where we get to the band-plot code
   from a "plot" block in which the data to be plotted (and
   hence also the band specification) are given in matrix
   form. By this point the plot-data have been converted to
   (temporary) DATASET form; here we retrieve the band-spec
   matrix, check it for conformability, and stick the two
   extra columns onto the dataset (borrowing pointers into
   the matrix content).
*/

static int process_band_matrix (const int *list,
                                DATASET *dset,
                                band_info *bi,
                                int **plist)
{
    const char *s = get_optval_string(PLOT, OPT_N);
    gretl_matrix *m = NULL;
    gchar **S;
    int i = 0;
    int err = 0;

    if (s == NULL) {
        return E_INVARG;
    }

    S = g_strsplit(s, ",", -1);

    while (S != NULL && S[i] != NULL && !err) {
        if (i == 0) {
            m = get_matrix_by_name(S[i]);
            if (m == NULL || m->cols != 2 || m->rows != dset->n) {
                /* missing or non-conformable */
                err = invalid_field_error(S[i]);
            } else {
                /* the last two series in expanded dataset */
                bi->center = dset->v;
                bi->width = dset->v + 1;
            }
        } else if (i == 1) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi->factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi->factor = gretl_scalar_get_value(S[i], &err);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    g_strfreev(S);

    if (!err && (bi->factor < 0 || na(bi->factor))) {
        err = E_INVARG;
    }

    if (!err) {
        /* enlarge the dset->Z array */
        int newv = dset->v + 2;
        double **tmp = realloc(dset->Z, newv * sizeof *tmp);

        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            /* note: we don't need varnames here */
            dset->Z = tmp;
            dset->Z[dset->v] = m->val;
            dset->Z[dset->v+1] = m->val + m->rows;
            dset->v += 2;
        }
    }

    if (!err) {
        *plist = gretl_list_copy(list);
        gretl_list_append_term(plist, bi->center);
        gretl_list_append_term(plist, bi->width);
    }

    return err;
}

/* Handle the band plus-minus option for all cases apart
   from the special one handled just above. Here we require
   two comma-separated series identifiers for center and
   width.
*/

static int parse_band_pm_option (const int *list,
                                 const DATASET *dset,
                                 gretlopt opt,
                                 band_info *bi,
                                 int **plist)
{
    int pci = get_effective_plot_ci();
    const char *s = get_optval_string(pci, OPT_N);
    gchar **S;
    int cpos = 0, wpos = 0;
    int v, pos, i = 0;
    int err = 0;

    if (s == NULL) {
        return E_INVARG;
    }

    if (strchr(s, ',') == NULL) {
        /* a single field: try for "recession bars" */
        v = current_series_index(dset, s);
        if (v >= 0 && v < dset->v) {
            if (gretl_isdummy(dset->t1, dset->t2, dset->Z[v])) {
                bi->bdummy = v;
            } else {
                err = E_INVARG;
                fprintf(stderr, "%s: not a dummy variable\n",
                        dset->varname[v]);
            }
        } else {
            err = E_INVARG;
        }
        return err;
    }

    /* at this point, can't be a recession-style band */
    S = g_strsplit(s, ",", -1);

    while (S != NULL && S[i] != NULL && !err) {
        if (i < 2) {
            /* specs for the "center" and "width" series: required */
            if (opt & OPT_X) {
                /* special for matrix-derived dataset */
                v = (i == 0)? dset->v - 2 : dset->v - 1;
            } else if (integer_string(S[i])) {
                /* var ID number? */
                v = atoi(S[i]);
            } else {
                /* varname? */
                v = current_series_index(dset, S[i]);
            }
            if (v >= 0 && v < dset->v) {
                pos = in_gretl_list(list, v);
                if (i == 0) {
                    bi->center = v;
                    cpos = pos;
                } else {
                    bi->width = v;
                    wpos = pos;
                }
            } else {
                err = invalid_field_error(S[i]);
            }
        } else if (i == 2) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi->factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi->factor = gretl_scalar_get_value(S[i], &err);
            } else {
                /* FIXME support a named vector? */
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    g_strfreev(S);

#if 0
    fprintf(stderr, "pm err = %d\n", err);
    fprintf(stderr, "pm center = %d (pos %d)\n", bi->center, cpos);
    fprintf(stderr, "pm width = %d (pos %d)\n", bi->width, wpos);
    fprintf(stderr, "pm factor = %g\n", bi->factor);
#endif

    if (!err) {
        if (bi->center < 0 || bi->width < 0 ||
            bi->factor < 0 || na(bi->factor)) {
            err = E_INVARG;
        }
    }

    if (!err && (cpos == 0 || wpos == 0)) {
        /* stick the "extra" series into *plist so we
           can check all series for NAs
        */
        *plist = gretl_list_copy(list);
        if (cpos == 0) {
            gretl_list_append_term(plist, bi->center);
        }
        if (wpos == 0) {
            gretl_list_append_term(plist, bi->width);
        }
    }

    return err;
}

/* We're looking here for any one of three patterns:

   <style>
   <style>,<color>
   <color>

   where <style> should be "fill", "dash" or "line" (the
   default) and <color> should be a hex string such as
   "#00ff00" or "0x00ff00".
*/

static int parse_band_style_option (band_info *bi)
{
    int pci = get_effective_plot_ci();
    const char *s = get_optval_string(pci, OPT_J);
    int err = 0;

    if (s != NULL) {
        const char *p = strchr(s, ',');

        if (bi->bdummy && *s != ',') {
            /* must be just a color */
            err = parse_gnuplot_color(s, bi->rgb);
        } else if (*s == ',') {
            /* skipping field 1, going straight to color */
            err = parse_gnuplot_color(s + 1, bi->rgb);
        } else if (p == NULL) {
            /* just got field 1, style spec */
            if (!strcmp(s, "fill")) {
                bi->style = BAND_FILL;
            } else if (!strcmp(s, "dash")) {
                bi->style = BAND_DASH;
            } else if (!strcmp(s, "line")) {
                bi->style = BAND_LINE;
            } else if (!strcmp(s, "bars")) {
                bi->style = BAND_BARS;
            } else if (!strcmp(s, "step")) {
                bi->style = BAND_STEP;
            } else {
                err = invalid_field_error(s);
            }
        } else {
            /* embedded comma: style + color */
            if (strlen(s) < 8) {
                err = invalid_field_error(s);
            } else if (!strncmp(s, "fill,", 5)) {
                bi->style = BAND_FILL;
            } else if (!strncmp(s, "dash,", 5)) {
                bi->style = BAND_DASH;
            } else if (!strncmp(s, "line,", 5)) {
                bi->style = BAND_LINE;
            } else if (!strncmp(s, "bars,", 5)) {
                bi->style = BAND_BARS;
            } else if (!strncmp(s, "step,", 5)) {
                bi->style = BAND_STEP;
            } else {
                err = invalid_field_error(s);
            }
            if (!err) {
                err = parse_gnuplot_color(s + 5, bi->rgb);
            }
        }
    }

    return err;
}

/* write "recession bars" as gnuplot rectangle objects, using
   the dummy variable @d for on/off information
*/

static int write_rectangles (gnuplot_info *gi,
                             char *rgb,
                             const double *d,
                             int t1, int t2,
                             DATASET *dset,
                             FILE *fp)
{
    char stobs[16], endobs[16];
    int bar_on = 0, obj = 1;
    int t, err = 0;

    if (gi->x == NULL) {
        return E_DATA;
    }

    if (*rgb == '\0') {
        strcpy(rgb, "#dddddd");
    }

    *stobs = *endobs = '\0';

    for (t=t1; t<=t2; t++) {
        if (na(d[t])) {
            err = E_MISSDATA;
            break;
        }
        if (bar_on && d[t] == 0) {
            /* finalize a bar */
            sprintf(endobs, "%g", gi->x[t]);
            fprintf(fp, "set object %d rectangle from %s, graph 0 to %s, graph 1 back "
                    "fillstyle solid 0.5 noborder fc rgb \"%s\"\n",
                    obj++, stobs, endobs, rgb);
            bar_on = 0;
        } else if (!bar_on && d[t] != 0) {
            /* start a bar */
            sprintf(stobs, "%g", gi->x[t]);
            bar_on = 1;
        }
    }

    if (bar_on) {
        /* terminate an unfinished bar */
        sprintf(endobs, "%g", gi->x[t2]);
        fprintf(fp, "set object rectangle from %s, graph 0 to %s, graph 1 back "
                "fillstyle solid 0.5 noborder fc rgb \"%s\"\n",
                stobs, endobs, rgb);
    }

    return err;
}

static int band_straddles_zero (const double *c,
                                const double *w,
                                double factor,
                                int t1, int t2)
{
    int t, lt0 = 0, gt0 = 0;
    double b1, b2;

    for (t=t1; t<=t2; t++) {
        b1 = c[t] - w[t] * factor;
        b2 = c[t] + w[t] * factor;
        if (b1 < 0 || b2 < 0) {
            lt0 = 1;
        }
        if (b1 > 0 || b2 > 0) {
            gt0 = 1;
        }
        if (lt0 && gt0) {
            return 1;
        }
    }

    return 0;
}

int plot_with_band (BPMode mode, gnuplot_info *gi,
                    const char *literal,
                    DATASET *dset,
                    gretlopt opt)
{
    band_info *bi = NULL;
    FILE *fp = NULL;
    const double *x = NULL;
    const double *y = NULL;
    const double *c = NULL;
    const double *w = NULL;
    const double *d = NULL;
    char yname[MAXDISP];
    char xname[MAXDISP];
    char wspec[16] = {0};
    int *biglist = NULL;
    int show_zero = 0;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int i, n_yvars = 0;
    int err = 0;

    bi = band_info_new();
    if (bi == NULL) {
	return E_ALLOC;
    }

    if (mode == BP_BLOCKMAT) {
        /* Coming from a "plot" block in matrix mode: in this case the
           band should be given in the form of a named matrix with
           two columns holding center and width, respectively.
        */
        err = process_band_matrix(gi->list, dset, bi, &biglist);
    } else {
        err = parse_band_pm_option(gi->list, dset, opt, bi, &biglist);
    }

    if (!err && (opt & OPT_J)) {
        err = parse_band_style_option(bi);
    }

    if (!err) {
        if (biglist != NULL) {
            err = graph_list_adjust_sample(biglist, gi, dset, 1);
        } else {
            err = graph_list_adjust_sample(gi->list, gi, dset, 1);
        }
        if (!err) {
            t1 = gi->t1;
            t2 = gi->t2;
        }
    }

    free(biglist);

    if (err) {
	free(bi);
        return err;
    }

    if (gi->flags & (GPT_TS | GPT_IDX)) {
        x = gi->x;
        *xname = '\0';
        if (gi->flags & GPT_TS) {
            gi->flags |= GPT_LETTERBOX;
        }
    } else {
        int xno = gi->list[gi->list[0]];

        x = dset->Z[xno];
        strcpy(xname, series_get_graph_name(dset, xno));
    }

    n_yvars = gi->list[0] - 1;

    fp = open_plot_input_file(PLOT_BAND, gi->flags, &err);
    if (err) {
	free(bi);
        return err;
    }

    /* assemble the data we'll need */
    if (bi->bdummy) {
        d = dset->Z[bi->bdummy];
    } else {
        c = dset->Z[bi->center];
        w = dset->Z[bi->width];
        show_zero = band_straddles_zero(c, w, bi->factor, t1, t2);
    }

    if (gi->flags & GPT_TS) {
        PRN *prn = gretl_print_new_with_stream(fp);

        make_time_tics(gi, dset, 0, NULL, prn);
        gretl_print_detach_stream(prn);
        gretl_print_destroy(prn);
    }

    if (n_yvars == 1) {
        fputs("set nokey\n", fp);
        strcpy(yname, series_get_graph_name(dset, gi->list[1]));
        fprintf(fp, "set ylabel \"%s\"\n", yname);
    }
    if (*xname != '\0') {
        fprintf(fp, "set xlabel \"%s\"\n", xname);
    }
    if (show_zero && bi->style != BAND_FILL) {
        fputs("set xzeroaxis\n", fp);
    }

    gretl_push_c_numeric_locale();

    if (gi->x != NULL) {
        /* FIXME case of gi->x == NULL? */
        print_x_range(gi, fp);
    }

    print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fp);

    if (bi->bdummy) {
        /* write out the rectangles */
        write_rectangles(gi, bi->rgb, d, t1, t2, dset, fp);
    }

    if (bi->bdummy) {
        int oddman = 0;

        if (!(opt & OPT_Y)) {
            check_for_yscale(gi, (const double **) dset->Z, &oddman);
            if (gi->flags & GPT_Y2AXIS) {
                fputs("set ytics nomirror\n", fp);
                fputs("set y2tics\n", fp);
            }
        }

        fputs("plot \\\n", fp);

        /* plot the actual data */
        for (i=1; i<=n_yvars; i++) {
            const char *iname = series_get_graph_name(dset, gi->list[i]);

            set_plot_withstr(gi, i, wspec);
            if (gi->flags & GPT_Y2AXIS) {
                fprintf(fp, "'-' using 1:2 axes %s title \"%s (%s)\" %s lt %d",
                        (i == oddman)? "x1y2" : "x1y1", iname,
                        (i == oddman)? _("right") : _("left"),
                        wspec, i);
            } else {
                fprintf(fp, "'-' using 1:2 title \"%s\" %s lt %d", iname, wspec, i);
            }
            if (i < n_yvars) {
                fputs(", \\\n", fp);
            } else {
                fputc('\n', fp);
            }
        }
        /* and write the data block */
        for (i=0; i<n_yvars; i++) {
            y = dset->Z[gi->list[i+1]];
            print_user_y_data(x, y, t1, t2, fp);
        }
        goto finish;
    }

    fputs("plot \\\n", fp);

    if (bi->style == BAND_FILL) {
        /* plot the confidence band first, so the other lines
           come out on top */
        print_pm_filledcurve_line(bi, NULL, fp);
        if (show_zero) {
            fputs("0 notitle w lines lt 0, \\\n", fp);
        }
        /* plot the non-band data */
        for (i=1; i<=n_yvars; i++) {
            const char *iname = series_get_graph_name(dset, gi->list[i]);

            set_plot_withstr(gi, i, wspec);
            fprintf(fp, "'-' using 1:2 title '%s' %s lt %d", iname, wspec, i);
            if (i == n_yvars) {
                fputc('\n', fp);
            } else {
                fputs(", \\\n", fp);
            }
        }
    } else {
        char lspec[24], dspec[8];

        *lspec = *dspec = '\0';

        /* plot the non-band data first */
        for (i=1; i<=n_yvars; i++) {
            const char *iname = series_get_graph_name(dset, gi->list[i]);

            set_plot_withstr(gi, i, wspec);
            fprintf(fp, "'-' using 1:2 title '%s' %s lt %d, \\\n", iname, wspec, i);
        }
        if (bi->rgb[0] != '\0') {
            sprintf(lspec, "lc rgb \"%s\"", bi->rgb);
        } else {
            sprintf(lspec, "lt %d", n_yvars + 1);
        }
        if (bi->style == BAND_DASH) {
            strcpy(dspec, " dt 2");
        }
        /* then the band */
        if (bi->style == BAND_BARS) {
            fprintf(fp, "'-' using 1:2:(%g*$3) notitle w errorbars %s%s\n",
                    bi->factor, lspec, dspec);
        } else {
            char *wstr = bi->style == BAND_STEP ? "steps" : "lines";

            fprintf(fp, "'-' using 1:($2-%g*$3) notitle w %s %s%s, \\\n",
                    bi->factor, wstr, lspec, dspec);
            fprintf(fp, "'-' using 1:($2+%g*$3) notitle w %s %s%s\n",
                    bi->factor, wstr, lspec, dspec);
        }
    }

    /* write out the inline data, the order depending on whether
       or not we're using fill style for the band
    */

    if (bi->style == BAND_FILL) {
        print_user_pm_data(x, c, w, t1, t2, fp);
        for (i=0; i<n_yvars; i++) {
            y = dset->Z[gi->list[i+1]];
            print_user_y_data(x, y, t1, t2, fp);
        }
    } else {
        for (i=0; i<n_yvars; i++) {
            y = dset->Z[gi->list[i+1]];
            print_user_y_data(x, y, t1, t2, fp);
        }
        print_user_pm_data(x, c, w, t1, t2, fp);
        if (bi->style != BAND_BARS) {
            print_user_pm_data(x, c, w, t1, t2, fp);
        }
    }

 finish:

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    clear_gpinfo(gi);
    free(bi);

    if (mode == BP_BLOCKMAT) {
        /* hide the two extra dataset columns
           representing the band */
        dset->v -= 2;
    }

    return err;
}
