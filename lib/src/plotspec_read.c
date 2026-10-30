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

/* Reading of gnuplot command files, as written by libgretl, into a
   GPT_SPEC struct, as needed by the GUI plot editor. The inverse
   operation (writing a script given a GPT_SPEC) is in plotspec.c.
*/

#include "libgretl.h"
#include "plotspec.h"
#include "libset.h"

#include <errno.h>

#define GPDEBUG 0
#define LT_DEBUG 0

/* the following needs more testing */
#define HANDLE_HEREDATA 1

/* max number of observations for which we'll read data markers */
#define MAX_MARKERS 250

typedef struct linestyle_ linestyle;

struct linestyle_ {
    char lc[20]; /* color */
    float lw;    /* line width */
    int dt;      /* dash type */
    int pt;      /* point type */
};

static int is_plot_line (const char *s)
{
    return !strncmp(s, "plot", 4) || !strncmp(s, "splot", 5);
}

/* chop trailing comma, if present; return 1 if comma chopped,
   zero otherwise */

static int chop_comma (char *str)
{
    size_t i, n = strlen(str);

    for (i=n-1; i>0; i--) {
        if (isspace((unsigned char) str[i])) {
            continue;
        }
        if (str[i] == ',') {
            str[i] = 0;
            return 1;
        } else {
            break;
        }
    }

    return 0;
}

/* returns non-zero on obtaining a marker */

static int get_gpt_marker (const char *line, char *label,
                           const char *format)
{
    const char *p = strchr(line, '#');

#if GPDEBUG > 1
    fprintf(stderr, "get_gpt_marker, p='%s'\n", p);
#endif

    *label = '\0';

    if (p != NULL) {
        sscanf(p + 2, format, label);
#if GPDEBUG > 1
        fprintf(stderr, "read marker: '%s'\n", label);
#endif
    }

    return *label != '\0';
}


#define gp_missing(s) (s[0] == '?' || !strcmp(s, "NaN"))

static int get_gpt_data (GPT_SPEC *spec,
                         int *do_markers,
                         const char *buf)
{
    char s[MAXLEN];
    char *got = NULL;
    double *x[5] = {NULL};
    char test[5][32];
    char obsfmt[12] = {0};
    int started_data_lines = 0;
    int i, j, t, imin = 0;
    int err = 0;

    spec->okobs = spec->nobs;

    gretl_push_c_numeric_locale();

    /* first handle "shaded bars" info, if present */

    for (i=0; i<spec->nbars && !err; i++) {
        double y1, y2, dx[2];

        /* start date */
        got = bufgets(s, sizeof s, buf);
        if (got == NULL ||
            sscanf(s, "%lf %lf %lf", &dx[0], &y1, &y2) != 3) {
            err = 1;
            break;
        }
        /* stop date */
        got = bufgets(s, sizeof s, buf);
        if (got == NULL ||
            sscanf(s, "%lf %lf %lf", &dx[1], &y1, &y2) != 3) {
            err = 1;
            break;
        }
        /* trailing 'e' */
        got = bufgets(s, sizeof s, buf);
        if (got) {
            plotspec_set_bar_info(spec, i, dx[0], dx[1]);
        } else {
            err = 1;
        }
    }

    if (*do_markers) {
        sprintf(obsfmt, "%%%d[^\r\n]", OBSLEN - 1);
    }

    /* then get the regular plot data */

    for (i=0; i<spec->n_lines && !err; i++) {
        int ncols = gp_line_data_columns(spec, i);
        int okobs = spec->nobs;
        int offset = 1;

        if (ncols == 0) {
            if (i == 0) {
                imin = 1;
            }
            continue;
        }

#if GPDEBUG > 1
        fprintf(stderr, "reading plot data, line %d, ncols %d\n", i, ncols);
#endif

        if (!started_data_lines) {
            offset = 0;
            x[0] = spec->data->val;
            x[1] = x[0] + spec->nobs;
            started_data_lines = 1;
        } else if (spec->code == PLOT_FACTORIZED) {
	    offset = 0;
	}

        x[2] = x[1] + spec->nobs;
        x[3] = x[2] + spec->nobs;
        x[4] = x[3] + spec->nobs;

        for (t=0; t<spec->nobs && !err; t++) {
            int nf, missing = 0;

            got = bufgets(s, sizeof s, buf);
            if (got == NULL) {
                err = 1;
                break;
            }

            nf = 0;
            if (ncols == 5) {
                nf = sscanf(s, "%31s %31s %31s %31s %31s", test[0], test[1], test[2],
                            test[3], test[4]);
            } else if (ncols == 4) {
                nf = sscanf(s, "%31s %31s %31s %31s", test[0], test[1], test[2], test[3]);
            } else if (ncols == 3) {
                nf = sscanf(s, "%31s %31s %31s", test[0], test[1], test[2]);
            } else if (ncols == 2) {
                nf = sscanf(s, "%31s %31s", test[0], test[1]);
            }
            if (nf != ncols) {
                err = 1;
            }

            for (j=offset; j<nf && !err; j++) {
                if (gp_missing(test[j])) {
                    x[j][t] = NADBL;
                    missing++;
                } else {
                    x[j][t] = atof(test[j]);
                }
            }
            if (missing) {
                okobs--;
            }

            if (i <= imin && *do_markers) {
                *do_markers = get_gpt_marker(s, spec->markers[t], obsfmt);
                if (spec->code == PLOT_FACTORIZED && imin == 0) {
		    /* FIXME ? */
                    imin = 1;
                }
            }
        }

        if (okobs < spec->okobs) {
            spec->okobs = okobs;
        }

        /* trailer line for data block */
        bufgets(s, sizeof s, buf);

	if (spec->code == PLOT_FACTORIZED) {
	    /* shift 'x'and 'y' write locations */
	    x[0] += ncols * spec->nobs;
	    x[1] += ncols * spec->nobs;
	} else {
	    /* shift 'y' write location only */
	    x[1] += (ncols - 1) * spec->nobs;
	}
    }

#if GPDEBUG > 1
    gretl_matrix_print(spec->data, "spec->data");
#endif

    gretl_pop_c_numeric_locale();

    return err;
}

static int get_gpt_heredata (GPT_SPEC *spec,
                             int *do_markers,
                             long barpos,
                             long datapos,
                             const char *buf)
{
    gretl_matrix *m = spec->data;
    char line[MAXLEN], test[32];
    char obsfmt[12] = {0};
    char *s, *got = NULL;
    double xij;
    int i, j;
    int err = 0;

    spec->okobs = spec->nobs;

    gretl_push_c_numeric_locale();

    /* first handle "shaded bars" info, if present */

    if (barpos > 0) {
        bufseek(buf, barpos);
        for (i=0; i<spec->nbars && !err; i++) {
            double y1, y2, dx[2];

            /* start date */
            got = bufgets(line, sizeof line, buf);
            if (got == NULL ||
                sscanf(line, "%lf %lf %lf", &dx[0], &y1, &y2) != 3) {
                err = 1;
                break;
            }
            /* stop date */
            got = bufgets(line, sizeof line, buf);
            if (got == NULL ||
                sscanf(line, "%lf %lf %lf", &dx[1], &y1, &y2) != 3) {
                err = 1;
                break;
            }
            plotspec_set_bar_info(spec, i, dx[0], dx[1]);
        }
    }

    if (*do_markers) {
        sprintf(obsfmt, "%%%d[^\r\n]", OBSLEN - 1);
    }

    /* then get the regular plot data */

    bufseek(buf, datapos);

    for (i=0; i<m->rows && !err; i++) {
        got = bufgets(line, sizeof line, buf);
        if (got == NULL) {
            err = 1;
        }
        s = line;
        for (j=0; j<m->cols && !err; j++) {
            s += strspn(s, " ");
            sscanf(s, "%31s", test);
            if (gp_missing(test)) {
                xij = NADBL;
            } else {
                xij = atof(test);
            }
            gretl_matrix_set(m, i, j, xij);
            s += strlen(test);
            if (*do_markers) {
                *do_markers = get_gpt_marker(s, spec->markers[i], obsfmt);
            }
        }
    }

    gretl_pop_c_numeric_locale();

#if GPDEBUG > 1
    gretl_matrix_print(m, "gp heredata");
#endif

    return err;
}

static int line_starts_heredata (const char *line,
                                 char **eod)
{
    int ret = 0;

    if (*line == '$') {
        const char *p = strstr(line, "<<");

        if (p != NULL) {
            int len;

            p += 2;
            p += strspn(p, " ");
            len = gretl_namechar_spn(p);
            if (len > 0) {
                *eod = gretl_strndup(p, len);
                ret = 1;
            }
        }
    }

    return ret;
}

/* This recognizes the limitations of the GPT_LABEL struct (see
   lib/src/plotspec.h). In principle that struct could be rejigged to
   hold information on "level" (front/back), point type and point size.
*/

static int unhandled_label_spec (const char *s)
{
    if (strstr(s, " front") || strstr(s, " back")) {
        return 1;
    } else if (strstr(s, " point") || strstr(s, " pt") ||
               strstr(s, "ps")) {
        return 1;
    } else {
        return 0;
    }
}

/* read a gnuplot source line specifying a text label */

static int parse_label_line (GPT_SPEC *spec, const char *line,
                             int *unhandled)
{
    const char *p, *s;
    char *text = NULL;
    double x, y;
    int nc, just = GP_JUST_LEFT;
    int err = 0;

    /* Examples:
       (a) set label "this is a label" at 1998.26,937.557 left front
       (b) set label 'foobar' at 1500,350 left
       (c) set label "" at 0.0127,0.6308 front center point pt 8 ps 2
    */

    /* find first double or single quote */
    p = strchr(line, '"');
    if (p == NULL) {
        p = strchr(line, '\'');
    }

    if (p == NULL) {
        /* no label text found */
        return 1;
    }

    /* note: text shouldn't be NULL but it may be empty */
    text = gretl_quoted_string_strdup(p, &s);
    if (text == NULL) {
        return 1;
    }

    /* get the position */
    p = strstr(s, "at");
    if (p == NULL) {
        err = E_DATA;
    } else {
        p += 2;
        gretl_push_c_numeric_locale();
        nc = sscanf(p, "%lf,%lf", &x, &y);
        gretl_pop_c_numeric_locale();
        if (nc != 2) {
            err = E_DATA;
        }
    }

    /* We need to watch out for case like (c) above, where the
       specification exceeds what we're currently able to pack into
       a GPT_LABEL struct in @spec.
    */
    if (!err && unhandled != NULL && unhandled_label_spec(p)) {
        *unhandled = 1;
        return 1;
    }

    if (!err) {
        /* justification */
        if (strstr(p, "right")) {
            just = GP_JUST_RIGHT;
        } else if (strstr(p, "center")) {
            just = GP_JUST_CENTER;
        }
    }

    if (!err) {
        err = plotspec_add_label(spec);
        if (!err) {
            int i = spec->n_labels - 1;

            strncat(spec->labels[i].text, text, PLOT_LABEL_TEXT_LEN);
            spec->labels[i].pos[0] = x;
            spec->labels[i].pos[1] = y;
            spec->labels[i].just = just;
        }
    }

    free(text);

    return err;
}

/* read a gnuplot source line specifying an arrow */

static int parse_arrow_line (GPT_SPEC *spec, const char *line,
                             int *unhandled)
{
    double x0, y0, x1, y1;
    char head[12] = {0};
    int n, err = 0;

    /* Examples:
       (a) set arrow from 2000,150 to 2000,500 nohead [ lt 0 ]
       (b) set arrow 1 from 0.1, graph 0 to 0.1, graph 1
       (c) set arrow 1 as 1
    */

    gretl_push_c_numeric_locale();
    n = sscanf(line, "set arrow from %lf,%lf to %lf,%lf %s",
               &x0, &y0, &x1, &y1, head);
    gretl_pop_c_numeric_locale();

    if (n < 5) {
        if (unhandled != NULL) {
            *unhandled = 1;
        }
        err = E_DATA;
    } else {
        err = plotspec_add_arrow(spec);
    }

    if (!err) {
        int flags = 0;

        if (strcmp(head, "nohead")) {
            flags |= GP_ARROW_HEAD;
        }
        if (strstr(line, "lt 0")) {
            flags |= GP_ARROW_DOTS;
        }
        n = spec->n_arrows - 1;
        spec->arrows[n].x0 = x0;
        spec->arrows[n].y0 = y0;
        spec->arrows[n].x1 = x1;
        spec->arrows[n].y1 = y1;
        spec->arrows[n].flags = flags;
    }

    return err;
}

static int
read_plotspec_range (const char *obj, const char *s, GPT_SPEC *spec)
{
    double r0, r1;
    int i = 0, err = 0;

    if (!strcmp(obj, "xrange")) {
        i = GP_X_RANGE;
    } else if (!strcmp(obj, "yrange")) {
        i = GP_Y_RANGE;
    } else if (!strcmp(obj, "y2range")) {
        i = GP_Y2_RANGE;
    } else if (!strcmp(obj, "trange")) {
        i = GP_T_RANGE;
    } else if (!strcmp(obj, "x2range")) {
        i = GP_X2_RANGE;
    } else if (!strcmp(obj, "cbrange")) {
        i = GP_CB_RANGE;
    } else {
        err = 1;
    }

    if (!err) {
        gretl_push_c_numeric_locale();
        if (!strcmp(s, "[*:*]")) {
            r0 = r1 = NADBL;
        } else if (strstr(s, "[*:")) {
            r0 = NADBL;
            r1 = atof(s + 3);
        } else if (strstr(s, ":*]")) {
            r0 = atof(s + 1);
            r1 = NADBL;
        } else if (sscanf(s, "[%lf:%lf]", &r0, &r1) != 2) {
            err = 1;
        }
        gretl_pop_c_numeric_locale();
    }

    if (!err) {
        spec->range[i][0] = r0;
        spec->range[i][1] = r1;
    }

    return err;
}

static int read_plot_logscale (const char *s, GPT_SPEC *spec)
{
    char axis[3] = {0};
    double base = 0;
    int i, n, err = 0;

    n = sscanf(s, "%2s %lf", axis, &base);

    if (n < 1 || (n == 2 && base < 1.1)) {
        err = 1;
    } else {
        if (n == 1) {
            base = 10.0;
        }
        if (!strcmp(axis, "x")) {
            i = 0;
        } else if (!strcmp(axis, "y")) {
            i = 1;
        } else if (!strcmp(axis, "y2")) {
            i = 2;
        } else {
            err = 1;
        }
    }

    if (!err) {
        spec->logbase[i] = base;
    }

    return err;
}

static int read_plot_format (const char *s, GPT_SPEC *spec,
                             int timefmt)
{
    char fmt[16];
    char axis = '\0';
    int n, err = 0;

    if (timefmt) {
        if (strcmp(s, "%s") != 0) {
            gretl_errmsg_set("Invalid 'timefmt' setting");
            err = 1;
        } else {
            spec->flags |= GPT_TIMEFMT;
        }
    } else {
        n = sscanf(s, "%c \"%15[^\"]", &axis, fmt);
        err = n < 2;
    }

    if (!err) {
        if (axis == 'x') {
            *spec->xfmt = '\0';
            strncat(spec->xfmt, fmt, 15);
        } else if (axis == 'y') {
            *spec->yfmt = '\0';
            strncat(spec->yfmt, fmt, 15);
        }
    }

    return err;
}

static int catch_value (char *targ, const char *src, int maxlen)
{
    int i, n, q = 0;
    int ret = 0;

    src += strspn(src, " \t\r\n");
    if (*src == '\'' || *src == '"') {
        q = *src;
        src++;
    }

    *targ = '\0';

    if (*src != '\0') {
        strncat(targ, src, maxlen - 1);
        n = strlen(targ);

        for (i=n-1; i>=0; i--) {
            if (isspace((unsigned char) targ[i])) {
                targ[i] = '\0';
            } else {
                break;
            }
        }
        if (targ[i] == '\'' || targ[i] == '"') {
            if (i == 0 && targ[i] == q) {
                ret = 1; /* accept empty string as "value" */
            }
            targ[i] = '\0';
        }
    }

    return ret || (*targ != '\0');
}

static void fix_old_roots_plot (GPT_SPEC *spec)
{
    strings_array_free(spec->literal, spec->n_literal);
    spec->literal = NULL;
    spec->n_literal = 0;

    spec->flags |= (GPT_POLAR | GPT_XZEROAXIS | GPT_YZEROAXIS);
    spec->border = 0;
    spec->keyspec = GP_KEY_NONE;
    strcpy(spec->xtics, "none");
    strcpy(spec->ytics, "none");
}

static int unhandled_gp_line_error (const char *s)
{
    int n = strlen(s);

    if (string_is_blank(s)) {
        return 0;
    }

    if (s[n-1] == '\n') {
        char *tmp = gretl_strndup(s, n-1);

        fprintf(stderr, "unhandled gnuplot command: '%s'\n", tmp);
        free(tmp);
    } else {
        fprintf(stderr, "unhandled gnuplot command: '%s'\n", s);
    }

    return 1;
}

static int parse_gp_unset_line (GPT_SPEC *spec, const char *s)
{
    char key[16] = {0};
    int err = 0;

    if (sscanf(s + 6, "%15s", key) != 1) {
        err = 1;
    } else if (!strcmp(key, "border")) {
        spec->border = 0;
    } else if (!strcmp(key, "key")) {
        spec->keyspec = GP_KEY_NONE;
    } else if (!strcmp(key, "xtics")) {
        strcpy(spec->xtics, "none");
    } else if (!strcmp(key, "ytics")) {
        strcpy(spec->ytics, "none");
    } else {
        /* unrecognized "unset" parameter */
        err = 1;
    }

    if (err) {
        unhandled_gp_line_error(s);
    }

    return err;
}

/* here we're appending a second specification (e.g. rotation) to a
   prior tics spec (e.g. "nomirror")
*/

static void xtics_concat (char *targ, const char *src, size_t sz)
{
    int n = strlen(targ);

    if (targ[n-1] != ' ') {
        gchar *tmp = g_strdup_printf(" %s", src);

        strncat(targ, tmp, sz - n);
        g_free(tmp);
    } else {
        strncat(targ, src, sz - n);
    }
}

static void read_xtics_setting (GPT_SPEC *spec,
                                const char *key,
                                const char *val)
{
    if (!strcmp(key, "x2tics")) {
        spec->x2ticstr = gretl_strdup(val);
    } else if (*val == '(') {
        /* the 'special' form of xtics */
        int len = strlen(val);

        if (val[len-1] == '\\') {
            fprintf(stderr, "got continuing xtics\n");
        }
        spec->xticstr = gretl_strdup(val);
    } else {
        /* regular xtics specification: this might come in
           two parts if rotation is specified
        */
        size_t sz = sizeof(spec->xtics) - 1;

        if (*spec->xtics == '\0') {
            strncat(spec->xtics, val, sz);
        } else if (!strcmp(spec->xtics, "none")) {
            *spec->xtics = '\0';
            strncat(spec->xtics, val, sz);
        } else {
            /* concatenate */
            xtics_concat(spec->xtics, val, sz);
        }
    }
}

/* Try to accept variant RGB specifications besides the one that's
   standard in gretl plot files, namely "#RRGGBB" in hex. The candidate
   string -- which might be a color name of up to 17 characters -- is in
   the @cstr argument; on success it's overwritten by a standard gnuplot
   hex string.
*/

static int verify_rgb (char *cstr)
{
    char *s = cstr;
    char delim = s[0];
    int err = 0;

    if (delim == '"' || delim == '\'') {
        /* valid input will be quoted */
        const char *p;

        s++;
        p = strchr(s, delim);
        if (p != NULL) {
            gretlRGB rgb;
            char test[18] = {0};
            int len = p - s;

            if (len >= 3 && len <= 17) {
                strncat(test, s, len);
                rgb = numeric_color_from_string(test, &err);
                if (!err) {
                    print_rgb_hash(cstr, rgb);
                }
            } else {
                /* too short or too long */
                err = E_DATA;
            }
        } else {
            /* unmatched delimiter */
            err = E_DATA;
        }
    } else if (strlen(cstr) > 19) {
	/* Unquoted input must fit the canonical "#RRGGBB" form (or
	   similar short token) that callers' fixed-size destination
	   buffers assume; anything longer is not a color we accept.
	*/
	err = E_DATA;
    }

    return err;
}

static int grab_rgb_spec (char *targ, size_t targlen, const char *src)
{
    char tmp[20];
    int err = 0;

    if (sscanf(src, "%19s", tmp)) {
        /* validate and convert to hex if needed */
        err = verify_rgb(tmp);
        if (!err && strlen(tmp) < targlen) {
            strcpy(targ, tmp);
        }
    } else {
        err = 1;
    }

    return err;
}

/* e.g. set linetype 1 pt 6 lc rgb "#1B9E77" # dark teal
        set linetype 2 pt 7 lc rgb "#D95F02" # dark orange
*/

static int parse_linetype (const char *s, linestyle *styles)
{
    const char *p;
    int n, i, idx = 0;
    int err = 0;

    /* read the 1-based index of the linetype */
    n = sscanf(s, " %d", &idx);

    if (n != 1 || idx <= 0 || idx > N_GP_LINETYPES) {
        /* get out on error */
        return 1;
    } else {
        /* convert index to zero-based */
        i = idx - 1;
    }

    if (!err && (p = strstr(s, " lc ")) != NULL) {
        /* check for a color specification (20 bytes
           allows for quoted named colors)
        */
        p += 4;
        if (!strncmp(p, "rgb ", 4)) {
            p += 4;
        }
        p += strspn(p, " ");
        err = grab_rgb_spec(styles[i].lc, sizeof styles[i].lc, p);
    }
    if (!err && (p = strstr(s, " lw ")) != NULL) {
        /* check line-width specification */
        float lw = 0.0;

        if (sscanf(p + 4, "%f", &lw)) {
            styles[i].lw = lw;
        } else {
            err = 1;
        }
    }
    if (!err && (p = strstr(s, " dt ")) != NULL) {
        /* check dash-type specification */
        int dt = 0;

        if (sscanf(p + 4, "%d", &dt)) {
            styles[i].dt = dt;
        } else {
            err = 1;
        }
    }
    if (!err && (p = strstr(s, " pt ")) != NULL) {
        /* check point-type specification */
        int pt = 0;

        if (sscanf(p + 4, "%d", &pt)) {
            styles[i].pt = pt;
        } else {
            err = 1;
        }
    }

#if LT_DEBUG
    fprintf(stderr, "styles[%d]: idx %d lc '%s' pt %d lw %g dt %d\n",
            i, idx, styles[i].lc, styles[i].pt, (double) styles[i].lw,
            styles[i].dt);
#endif

    return err;
}

static int parse_gp_set_line (GPT_SPEC *spec,
                              const char *s,
                              linestyle *styles,
                              int *unhandled)
{
    char *p, key[16] = {0};
    char val[MAXLEN] = {0};
    char *extra = NULL;
    int lt_pos = 0;
    int err = 0;

    if (!strncmp(s, "set linetype", 12)) {
        lt_pos = 12;
    } else if (!strncmp(s, "set style fill solid", 20)) {
        /* special with multi-word key */
        strncat(val, s + 20, 8);
        spec->fillfrac = (float) atof(val);
        val[0] = '\0';
    }

    if (lt_pos > 0) {
        err = parse_linetype(s + lt_pos, styles);
        if (err && unhandled != NULL) {
            *unhandled = 1;
        }
        return err;
    }

    if (unhandled != NULL) {
        /* we're peeking at "literal" lines */
        if (strstr(s, "set style data histogram")) {
            spec->code = PLOT_BAR;
            *unhandled = 1;
            return 0;
        } else if (strstr(s, "histogram rowstacked")) {
            spec->code = PLOT_STACKED_BAR;
            *unhandled = 1;
            return 0;
        }
    }

    if (sscanf(s + 4, "%11s", key) != 1) {
        return unhandled_gp_line_error(s);
    }

    if ((p = strchr(key, '[')) != NULL) {
        /* try fixing this */
        *p = '\0';
    }

#if GPDEBUG
    fprintf(stderr, "parse_gp_set_line: key = '%s'\n", key);
#endif

    if (!strcmp(key, "term") ||
        !strcmp(key, "output") ||
        !strcmp(key, "encoding") ||
        !strcmp(key, "datafile")) {
        /* we ignore these */
        return 0;
    }

    /* first, settings that don't require a parameter */

    if (!strcmp(key, "y2tics")) {
        spec->flags |= GPT_Y2AXIS;
        return 0;
    } else if (!strcmp(key, "parametric")) {
        spec->flags |= GPT_PARAMETRIC;
        return 0;
    } else if (!strcmp(key, "polar")) {
        spec->flags |= GPT_POLAR;
        return 0;
    } else if (!strcmp(key, "mono") ||
               !strcmp(key, "monochrome")) {
        spec->flags |= GPT_MONO;
        return 0;
    } else if (!strcmp(key, "xzeroaxis")) {
        spec->flags |= GPT_XZEROAXIS;
        return 0;
    } else if (!strcmp(key, "yzeroaxis")) {
        spec->flags |= GPT_YZEROAXIS;
        return 0;
    } else if (!strcmp(key, "noxtics")) {
        strcpy(spec->xtics, "none");
        return 0;
    } else if (!strcmp(key, "noytics")) {
        strcpy(spec->ytics, "none");
        return 0;
    } else if (!strcmp(key, "nokey")) {
        spec->keyspec = GP_KEY_NONE;
        return 0;
    } else if (!strcmp(key, "label")) {
        if (parse_label_line(spec, s, unhandled) && unhandled != NULL) {
            /* nothing was stored (e.g. position given in "graph"
               coordinates), so the line must not be treated as
               handled */
            *unhandled = 1;
        }
        return 0;
    } else if (!strcmp(key, "arrow")) {
        parse_arrow_line(spec, s, unhandled);
        return 0;
    }

    /* grid lines: parameter is optional */
    if (!strcmp(key, "grid")) {
        if (catch_value(val, s + 4 + strlen(key), MAXLEN)) {
            if (!strcmp(val, "ytics")) {
                spec->flags |= GPT_GRID_Y;
            } else if (!strcmp(val, "xtics")) {
                spec->flags |= GPT_GRID_X;
            }
        } else {
            spec->flags |= GPT_GRID_Y;
            spec->flags |= GPT_GRID_X;
        }
    }

    /* now catch value for settings that need a parameter */

    if (!catch_value(val, s + 4 + strlen(key), MAXLEN)) {
        return 0;
    }

#if GPDEBUG
    fprintf(stderr, " value = '%s'\n", val);
#endif

    if (strstr(key, "range")) {
        if (read_plotspec_range(key, val, spec)) {
            return unhandled_gp_line_error(s);
        }
    } else if (!strcmp(key, "logscale")) {
        if (read_plot_logscale(val, spec)) {
            return unhandled_gp_line_error(s);
        }
    } else if (!strcmp(key, "format")) {
        if (read_plot_format(val, spec, 0)) {
            return unhandled_gp_line_error(s);
        }
    } else if (!strcmp(key, "timefmt")) {
        if (read_plot_format(val, spec, 1)) {
            return unhandled_gp_line_error(s);
        }
    } else if (!strcmp(key, "title")) {
        spec->titles[0] = g_strdup(val);
    } else if (!strcmp(key, "xlabel")) {
        spec->titles[1] = g_strdup(val);
        *spec->xvarname = '\0';
        strncat(spec->xvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "ylabel")) {
        spec->titles[2] = g_strdup(val);
        *spec->yvarname = '\0';
        strncat(spec->yvarname, val, MAXDISP-1);
    } else if (!strcmp(key, "y2label")) {
        spec->titles[3] = g_strdup(val);
    } else if (!strcmp(key, "x2label")) {
        spec->titles[4] = g_strdup(val);
    } else if (!strcmp(key, "key")) {
        spec->keyspec = gp_keyspec_from_string(val);
    } else if (!strcmp(key, "xtics") || !strcmp(key, "x2tics")) {
        read_xtics_setting(spec, key, val);
    } else if (!strcmp(key, "mxtics")) {
        *spec->mxtics = '\0';
        strncat(spec->mxtics, val, sizeof(spec->mxtics) - 1);
    } else if (!strcmp(key, "ytics")) {
        *spec->ytics = '\0';
        strncat(spec->ytics, val, sizeof(spec->ytics) - 1);
    } else if (!strcmp(key, "border")) {
        spec->border = atoi(val);
        if ((extra = strstr(val, "lc rgb \"")) != NULL) {
            spec->border_lc[0] = '\0';
            strncat(spec->border_lc, extra + 8, 7);
        }
    } else if (!strcmp(key, "bmargin")) {
        spec->bmargin = atoi(val);
    } else if (!strcmp(key, "boxwidth")) {
        spec->boxwidth = (float) atof(val);
        if (strstr(s, "absolute")) {
            spec->boxwidth = -spec->boxwidth;
        }
    } else if (!strcmp(key, "samples")) {
        spec->samples = atoi(val);
    } else if (!strcmp(key, "xdata")) {
        if (!strncmp(val, "time", 4)) {
            spec->flags |= GPT_TIMEFMT;
        }
    } else if (unhandled != NULL) {
        *unhandled = 1;
    }

    return 0;
}

/* allocate markers for identifying particular data points */

static int plotspec_allocate_markers (GPT_SPEC *spec)
{
    spec->markers = strings_array_new_with_length(spec->nobs,
                                                  OBSLEN);
    if (spec->markers == NULL) {
        spec->n_markers = 0;
        return E_ALLOC;
    } else {
        spec->n_markers = spec->nobs;
        return 0;
    }
}

static void plotspec_destroy_markers (GPT_SPEC *spec)
{
    strings_array_free(spec->markers, spec->n_markers);
    spec->markers = NULL;
    spec->n_markers = 0;
}

static int get_boxplot_auxdata (GPT_SPEC *spec, char *line,
                                const char *buf)
{
    double x, y;
    int i, n, r = 0, c = 0;
    int err = 0;

    n = sscanf(line + 9, "%d %d", &r, &c);

    if (n == 2 && r > 0 && c == 2) {
        spec->auxdata = gretl_matrix_alloc(r, c);
        if (spec->auxdata != NULL) {
            gretl_push_c_numeric_locale();
            if (spec->heredata) {
                /* skip a line: "$aux << EOA" */
                if (bufgets(line, MAXLEN - 1, buf) == NULL) {
                    err = E_DATA;
                }
            }
            for (i=0; i<r && !err; i++) {
                if (bufgets(line, MAXLEN - 1, buf) == NULL) {
                    err = E_DATA;
                } else if (sscanf(line, "%lf %lf", &x, &y) != 2) {
                    err = E_DATA;
                } else {
                    gretl_matrix_set(spec->auxdata, i, 0, x);
                    gretl_matrix_set(spec->auxdata, i, 1, y);
                }
            }
            gretl_pop_c_numeric_locale();
            if (err) {
                gretl_matrix_free(spec->auxdata);
                spec->auxdata = NULL;
            }
        }
    }

    return err;
}

/* Determine the number of data points in a plot. While we're at it,
   determine the type of plot, check whether there are any data-point
   markers along with the data, and see if there are are definitions of
   time-series vertical bars.
*/

static void get_plot_nobs (GPT_SPEC *spec,
                           PlotReadInfo *info,
                           const char *buf,
                           int *do_markers,
                           long *barpos,
                           long *datapos)
{
    int nobs = 0, started = -1;
    int startmin = 1;
    long auxpos = 0;
    char line[MAXLEN], test[12];
    char *eod = NULL;
    char *p = NULL;

    spec->code = PLOT_REGULAR;
    *do_markers = 0;

    while (bufgets(line, MAXLEN - 1, buf)) {
        if (*line == '#' && spec->code == PLOT_REGULAR) {
            tailstrip(line);
            spec->code = plot_type_from_string(line);
            if (spec->code == PLOT_GEOMAP) {
                int w, h;

                if (sscanf(line + 9, "%d %d", &w, &h) == 2) {
                    info->pixel_width = w;
                    info->pixel_height = h;
                }
                break;
            }
        }
        if (sscanf(line, "# n_bars = %d", &spec->nbars) == 1) {
            /* for old-style data representation */
            startmin += spec->nbars;
            continue;
        }
        if (spec->nbars > 0 && !strncmp(line, "$bars <", 7)) {
            *barpos = buftell(buf);
        }
        if (!strncmp(line, "# auxdata", 9)) {
            auxpos = buftell(buf);
        }
        if (started < 0 && line_starts_heredata(line, &eod)) {
            /* "heredoc" method of handling plot data */
            spec->heredata = 1;
            *datapos = buftell(buf);
            continue;
        }
        if (eod != NULL) {
            if (!strncmp(line, eod, strlen(eod))) {
                free(eod);
                eod = NULL;
                break;
            } else if (*do_markers == 0 && (p = strchr(line, '#')) != NULL) {
                if (sscanf(p + 1, "%8s", test) == 1) {
                    *do_markers = 1;
                }
            }
            nobs++;
        } else {
            if (is_plot_line(line)) {
                started = 0;
            }
            if (started == 0 && strchr(line, '\\') == NULL) {
                started = 1;
                continue;
            }
            if (started > 0 && started < startmin) {
                if (*line == 'e') {
                    started++;
                    continue;
                }
            }
            if (started == startmin) {
                if (*line == 'e' || !strncmp(line, "set ", 4)) {
                    /* end of data, or onto "set print" for bounds */
                    break;
                } else if (*do_markers == 0 && (p = strchr(line, '#')) != NULL) {
                    if (sscanf(p + 1, "%11s", test) == 1) {
                        *do_markers = 1;
                    }
                }
                nobs++;
            }
        }
    }

    spec->nobs = nobs;

    if (spec->code == PLOT_BOXPLOTS) {
        /* In the case of boxplots the gnuplot file may contain extra
           (outlier) data, which have to be read into an auxiliary
           matrix.
        */
        if (auxpos > 0) {
            /* we already went past it */
            bufseek(buf, auxpos);
            buf_back_lines(buf, 1);
            bufgets(line, MAXLEN - 1, buf);
        } else if (!spec->heredata) {
            /* try reading further */
            while (bufgets(line, MAXLEN - 1, buf)) {
                if (!strncmp(line, "# auxdata", 9)) {
                    auxpos = buftell(buf);
                    break;
                }
            }
        }
        if (auxpos > 0) {
            get_boxplot_auxdata(spec, line, buf);
        }
    }

#if GPDEBUG
    fprintf(stderr, "get_plot_nobs: %d (factorized? %d)\n",
	    spec->nobs, spec->code == PLOT_FACTORIZED);
#endif
}

static int grab_fit_coeffs (GPT_SPEC *spec, const char *s)
{
    gretl_matrix *b = NULL;
    int k, f = spec->fit;
    int n = 0, err = 0;

    k = (f == PLOT_FIT_CUBIC)? 4 : (f == PLOT_FIT_QUADRATIC)? 3 : 2;

    b = gretl_column_vector_alloc(k);
    if (b == NULL) {
        return E_ALLOC;
    }

    if (!strncmp(s, "exp(", 4)) {
        s += 4;
    }

    gretl_push_c_numeric_locale();

    if (k == 2) {
        n = sscanf(s, "%lf + %lf", &b->val[0], &b->val[1]);
    } else if (k == 3) {
        n = sscanf(s, "%lf + %lf*x + %lf", &b->val[0],
                   &b->val[1], &b->val[2]);
    } else if (k == 4) {
        n = sscanf(s, "%lf + %lf*x + %lf*x**2 + %lf", &b->val[0],
                   &b->val[1], &b->val[2], &b->val[3]);
    }

    gretl_pop_c_numeric_locale();

    if (n != k) {
        err = E_DATA;
        gretl_matrix_free(b);
        spec->flags &= ~GPT_AUTO_FIT;
    } else if (f == PLOT_FIT_OLS) {
        spec->b_ols = b;
    } else if (f == PLOT_FIT_INVERSE) {
        spec->b_inv = b;
    } else if (f == PLOT_FIT_LOGLIN) {
        spec->b_log = b;
    } else if (f == PLOT_FIT_QUADRATIC) {
        spec->b_quad = b;
    } else if (f == PLOT_FIT_CUBIC) {
        spec->b_cub = b;
    } else if (f == PLOT_FIT_LINLOG) {
        spec->b_linlog = b;
    }

    return err;
}

/* scan the stuff after "title '" or 'title "' */

static void grab_line_title (GPT_LINE *line, const char *src)
{
    char tmp[128];

    *tmp = '\0';

    if (*src == '\'') {
        sscanf(src + 1, "%127[^']'", tmp);
    } else {
        sscanf(src + 1, "%127[^\"]\"", tmp);
    }

    if (*tmp != '\0') {
        line->title = g_strdup(tmp);
    }
}

/* Examples:

   using 1:2
   using 1:n:m:p
   using 1:xtic(...)
   using 1:2:xtic(""):ytic("")
   using 1:($2+x)

   we need to extract the data column numbers, and if we find special
   stuff inline record the 'using' string as a whole
*/

static int process_using_spec (const char **ps,
                               GPT_SPEC *spec,
                               int i)
{
    GPT_LINE *line = &spec->lines[i];
    const char *f, *p, *s = *ps;
    int *cols = NULL;
    int inparen = 0;
    int anyparen = 0;
    int ticspecs = 0;
    int n_uniq = 0;
    int k, n = 0;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "process_using_spec: s = '%s'\n", s);
#endif

    if (isdigit(*s)) {
        k = atoi(s);
        gretl_list_append_term(&cols, k);
        n_uniq++;
    }

    p = s;

    while (*s) {
        if (*s == '(') {
            anyparen = 1;
            inparen++;
        } else if (*s == ')') {
            inparen--;
        } else if (!inparen) {
            if (*s == ':') {
                f = s + 1;
                if (*f == 'x' || *f == 'y') {
                    ticspecs++;
                } else {
                    k = (*f == '$')? atoi(f + 1) : atoi(f);
                    if (k > 0) {
                        if (!in_gretl_list(cols, k)) {
                            n_uniq++;
                        }
                        gretl_list_append_term(&cols, k);
                    }
                }
            } else if (isspace(*s) || *s == ',') {
                break;
            }
        } else if (*s == '$' && isdigit(*(s+1))) {
            /* inparen: a tic specification may make reference
               to a data column */
            k = atoi(s + 1);
            if (k > 0) {
                if (!in_gretl_list(cols, k)) {
                    n_uniq++;
                }
                gretl_list_append_term(&cols, k);
            }
        }
        s++;
        n++;
    }

    /* pass back pointer to remainder of line, if any */
    *ps = s;

    line->style = GP_STYLE_AUTO;

    if (n > 0 && (anyparen || ticspecs > 0)) {
        /* record the 'using' string as is */
        char *ustr = gretl_strndup(p, n);
#if GPDEBUG
        fprintf(stderr, "saving ustr = '%s'\n", ustr);
#endif
        line->ustr = ustr;
    }

    if (cols != NULL) {
        line->ncols = n_uniq;
        if (cols[0] == 5 && n_uniq == 2) {
            /* boxplot special */
            line->flags |= GP_LINE_BOXDATA;
        }
        if (spec->heredata) {
            gretl_list_merge_list(&spec->heredata_cols, cols);
            spec->n_datacols = spec->heredata_cols[0];
        } else if (spec->n_datacols > 0 &&
		   in_gretl_list(cols, 1) &&
		   spec->code != PLOT_FACTORIZED) {
            /* col 1 should already be counted */
            spec->n_datacols += line->ncols - 1;
        } else {
            spec->n_datacols += line->ncols;
        }
#if GPDEBUG
        fprintf(stderr, "number of unique columns %d\n", n_uniq);
        printlist(cols, "cols list");
        fprintf(stderr, "spec->n_datacols now = %d\n", spec->n_datacols);
#endif
        if (line->ustr == NULL) {
            line->mcols = cols;
        } else {
            free(cols);
        }
    }

    return err;
}

/* Parse the "using..." portion of plot specification for a
   given plot element. The full form is like:

     using <cols> axes <axes> title <title> w XX lt XX lw XX
     using <cols> title <title> lc rgb <rgb> w <style>
*/

static int parse_gp_line_line (const char *s, GPT_SPEC *spec,
                               int auto_linewidth)
{
    GPT_LINE *line;
    char tmp[16];
    const char *p;
    int i, err;

#if GPDEBUG
    fprintf(stderr, "parse_gp_line_line, starting, heredata = %d\n",
            spec->heredata);
#endif

    err = plotspec_add_line(spec);
    if (err) {
        return err;
    }

    i = spec->n_lines - 1;
    line = &spec->lines[i];

    if (!strncmp(s, "plot ", 5)) {
        s += 5;
    } else if (!strncmp(s, "splot ", 6)) {
        s += 6;
    }
    s += strspn(s, " ");

    if ((p = strstr(s, " using "))) {
        /* data column spec */
        p += 7;
        err = process_using_spec(&p, spec, i);
        s = p; /* pointer to remainder of line */
    } else if (*s == '\'' || *s == '"') {
        /* name of data file, without 'using': implicitly
           using cols 1 and 2
        */
        if (*(s+1) != '-') {
            fprintf(stderr, "plotting datafile, not supported\n");
        } else if (spec->heredata) {
            int ltmp[] = {2, 1, 3};

            line->ncols = 2;
            gretl_list_merge_list(&spec->heredata_cols, ltmp);
            spec->n_datacols = spec->heredata_cols[0];
        } else {
            line->ncols = 2;
            spec->n_datacols += 2;
        }
    } else {
        /* Absence of "using" should mean that the line plots
           a formula, not a set of data columns. So get the
           formula: it runs up to "title" or "notitle".
        */
        p = strstr(s, " title");
        if (p == NULL) {
            p = strstr(s, " notitle");
        }
        if (p != NULL) {
            line->formula = g_strndup(s, p - s);
            if (i == 1 && spec->flags & GPT_AUTO_FIT) {
                grab_fit_coeffs(spec, line->formula);
            }
        } else {
            /* title must be implicit? */
            gchar *tmp = g_strstrip(g_strdup(s));

            line->formula = tmp;
            line->title = g_strdup(tmp);
        }
    }

    if (strstr(s, "axes x1y2")) {
        line->yaxis = 2;
    }
    if ((p = strstr(s, " title "))) {
        grab_line_title(line, p + 7);
    }
    if ((p = strstr(s, " w "))) {
        sscanf(p + 3, "%15[^, ]", tmp);
        line->style = gp_style_index_from_name(tmp);
    } else if ((p = strstr(s, " with "))) {
        sscanf(p + 6, "%15[^, ]", tmp);
        line->style = gp_style_index_from_name(tmp);
    }
    if ((p = strstr(s, " lt "))) {
        sscanf(p + 4, "%d", &line->type);
    } else if ((p = strstr(s, " lc rgb "))) {
        grab_rgb_spec(line->rgb, sizeof line->rgb, p + 8);
    }
    if ((p = strstr(s, " ps "))) {
        sscanf(p + 4, "%f", &line->pscale);
    }
    if ((p = strstr(s, " pt "))) {
        sscanf(p + 4, "%d", &line->ptype);
    }
    if ((p = strstr(s, " dt "))) {
        sscanf(p + 4, "%d", &line->dtype);
    }
    if (!auto_linewidth && (p = strstr(s, " lw "))) {
        sscanf(p + 4, "%f", &line->width);
    }
    if ((p = strstr(s, " whiskerbars "))) {
        double ww;

        sscanf(p + 13, "%lf", &ww);
        line->whiskwidth = (float) ww;
    }

    if (line->ncols == 0 && line->formula == NULL) {
        /* got neither data column spec nor formula */
        err = 1;
    }

#if GPDEBUG
    fprintf(stderr, "parse_gp_line_line, returning %d\n", err);
#endif

    return err;
}

/* We got a special comment supposedly indicating the name and ID number
   of the X or Y variable in an OLS fitted line. Here we check this info
   for validity: @vname should be the name of a bona fide series
   variable, and its ID number should match the given @v.  We return 1
   if this works out OK, otherwise 0.
*/

static int plot_ols_var_ok (const DATASET *dset,
                            const char *vname, int v)
{
    int vcheck;

    if (dset == NULL) {
        return 0;
    }

    vcheck = current_series_index(dset, vname);

    return vcheck == v;
}

static void maybe_set_add_fit_ok (GPT_SPEC *spec,
                                  const DATASET *dset)
{
    if (spec->n_lines == 2 && spec->fit != PLOT_FIT_NONE) {
        ; /* already OK */
    } else if (spec->data != NULL &&
               spec->code == PLOT_REGULAR &&
               spec->n_lines == 1 &&
               spec->lines[0].ncols == 2) {
        if (spec->flags & GPT_TS) {
            spec->fit = (dset != NULL && dataset_is_time_series(dset)) ?
                PLOT_FIT_NONE : PLOT_FIT_NA;
        } else {
            spec->fit = PLOT_FIT_NONE;
        }
    } else {
        spec->fit = PLOT_FIT_NA;
    }
}

static int plot_get_data_and_markers (GPT_SPEC *spec,
                                      const char *buf,
                                      int *do_markers,
                                      long barpos,
                                      long datapos)
{
    int err = 0;

    if (spec->nobs == 0 || spec->n_datacols == 0) {
        /* nothing to be done */
        return 0;
    }

#if GPDEBUG
    fprintf(stderr, "plot_get_data, allocating: nobs=%d, n_datacols=%d",
            spec->nobs, spec->n_datacols);
    if (spec->code == PLOT_FACTORIZED) {
	fputs(" (factorized)\n", stderr);
    } else {
	fputc('\n', stderr);
    }
#endif

    /* allocate for the plot data... */
    spec->data = gretl_matrix_alloc(spec->nobs, spec->n_datacols);
    if (spec->data == NULL) {
        err = E_ALLOC;
    }

    /* and markers if any */
    if (!err && *do_markers) {
        err = plotspec_allocate_markers(spec);
    }

    /* and time-series bars, if any */
    if (!err && spec->nbars > 0) {
        err = plotspec_allocate_bars(spec);
    }

    /* read the data (and perhaps markers) from the plot file */
    if (!err) {
        if (spec->heredata) {
            err = get_gpt_heredata(spec, do_markers,
                                   barpos, datapos, buf);
        } else {
            err = get_gpt_data(spec, do_markers, buf);
        }
    } else {
        gretl_matrix_free(spec->data);
        spec->data = NULL;
    }

    if (spec->markers != NULL && *do_markers == 0) {
        /* something must have gone wrong: clean up */
        plotspec_destroy_markers(spec);
    }

#if GPDEBUG
    fprintf(stderr, "plot_get_data_and_markers:\n"
            " spec->data = %p, spec->markers = %p, spec->n_markers = %d, err = %d\n",
            spec->data, (void *) spec->markers, spec->n_markers, err);
#endif

    return err;
}

/* Parse special comment to get the 0-based index numbers of any
   user-defined lines that have been added to the plot.  The comment
   looks like:

   # %d user-defined lines: %d %d ...
*/

static int *get_user_lines_list (const char *s)
{
    int *list = NULL;
    char id[6];
    int i, n = 0;

    sscanf(s, "# %d", &n);

    if (n > 0) {
        s = strchr(s, ':');
        if (s != NULL) {
            list = gretl_list_new(n);
            if (list != NULL) {
                s++;
                for (i=0; i < n && *s != '\0'; i++) {
                    s++;
                    if (sscanf(s, "%5[^ ]", id)) {
                        list[i+1] = atoi(id);
                        s += strlen(id);
                    }
                }
            }
        }
    }

    return list;
}

static FitType recognize_fit_string (const char *s)
{
    if (strstr(s, "OLS")) {
        return PLOT_FIT_OLS;
    } else if (strstr(s, "quadratic")) {
        return PLOT_FIT_QUADRATIC;
    } else if (strstr(s, "cubic")) {
        return PLOT_FIT_CUBIC;
    } else if (strstr(s, "inverse")) {
        return PLOT_FIT_INVERSE;
    } else if (strstr(s, "loess")) {
        return PLOT_FIT_LOESS;
    } else if (strstr(s, "semilog")) {
        return PLOT_FIT_LOGLIN;
    } else if (strstr(s, "linlog")) {
        return PLOT_FIT_LINLOG;
    } else {
        return PLOT_FIT_NONE;
    }
}

static void check_for_plot_size (GPT_SPEC *spec, gchar *buf)
{
    char line[128];
    int i = 0;

    while (bufgets(line, sizeof line, buf) && i < 6) {
        if (!strncmp(line, "# multiple ", 11)) {
            if (strstr(line, "extra-large")) {
                spec->flags |= GPT_XXL;
                break;
            } else if (strstr(line, "large")) {
                spec->flags |= GPT_XL;
                break;
            }
        } else if (strstr(line, "extra-wide")) {
            spec->flags |= GPT_XW;
            break;
        }
        i++;
    }
}

static void check_for_plot_time_data (GPT_SPEC *spec, gchar *buf)
{
    char line[128];

    while (bufgets(line, sizeof line, buf)) {
        if (!strncmp(line, "set timefmt ", 12)) {
            spec->flags |= GPT_TIMEFMT;
        } else if (*line == '#' && strstr(line, "letterbox")) {
            spec->flags |= GPT_LETTERBOX;
        } else if (!strncmp(line, "plot", 4)) {
            break;
        }
    }
}

static void linestyle_init (linestyle *ls)
{
    ls->lc[0] = '\0';
    ls->lw = 0;
    ls->dt = 0;
    ls->pt = 0;
}

static int push_z_row (gretl_matrix *z, int i, int n, char *line)
{
    char *p = line;
    double x;
    int j, err = 0;

    errno = 0;

    for (j=0; j<n && *p; j++) {
        if (*p == ' ') {
            p++;
        }
        if (*p == '?') {
            x = NADBL;
            p += 2;
        } else if (!strncmp(p, "NaN", 3)) {
            x = NADBL;
            p += 4;
        } else {
            x = strtod(p, &p);
            if (errno) {
                err = 1;
                break;
            }
        }
        gretl_matrix_set(z, i, j, x);
    }

    return err;
}

static void get_heatmap_matrix (GPT_SPEC *spec, gchar *buf,
                                char *line, size_t len,
                                long datapos)
{
    gretl_matrix *z;
    int n = spec->nobs;
    int i, err = 0;

    z = gretl_matrix_alloc(n, n);
    if (z == NULL) {
        return;
    }

    bufseek(buf, datapos);

    gretl_push_c_numeric_locale();

    for (i=0; i<n && !err; i++) {
        bufgets(line, len, buf);
        err = push_z_row(z, i, n, line);
    }

    gretl_pop_c_numeric_locale();

    if (err) {
        gretl_matrix_free(z);
    } else {
        spec->auxdata = z;
    }
}

/* Here we're seeing if we should "promote" literal plot lines to
   supplement or override "set" variables that are recognized by gretl
   and form part of @spec. If we don't do this the GUI plot editor may
   be broken. However, at the same time we have to be careful not to
   break literal lines by mistakenly supposing they can be parlayed into
   elements of @spec.
*/

static void maybe_promote_literal_lines (GPT_SPEC *spec,
                                         linestyle *styles)
{
    const char *line;
    int *rmlines = NULL;
    int unhandled;
    int i, err = 0;

    for (i=0; i<spec->n_literal; i++) {
        line = spec->literal[i];
        if (!strncmp(line, "set ", 4)) {
            unhandled = 0;
            gretl_push_c_numeric_locale();
            err = parse_gp_set_line(spec, line, styles, &unhandled);
            gretl_pop_c_numeric_locale();
            if (!err && !unhandled) {
                gretl_list_append_term(&rmlines, i);
            }
        } else if (!strncmp(line, "unset ", 6)) {
            err = parse_gp_unset_line(spec, line);
            if (!err) {
                gretl_list_append_term(&rmlines, i);
            }
        }
    }

    if (rmlines != NULL) {
        /* at least one "literal" line has been promoted */
        if (rmlines[0] == spec->n_literal) {
            strings_array_free(spec->literal, spec->n_literal);
            spec->literal = NULL;
            spec->n_literal = 0;
        } else {
            char **S = spec->literal;
            int ns, orig_ns = spec->n_literal;

            ns = orig_ns - rmlines[0];
            spec->literal = strings_array_new(ns);

            if (spec->literal == NULL) {
                /* unlikely, but... */
                spec->literal = S;
            } else {
                /* reduced version of the "literal" array */
                int j = 0;

                for (i=0; i<orig_ns; i++) {
                    if (in_gretl_list(rmlines, i)) {
                        free(S[i]);
                    } else {
                        spec->literal[j++] = S[i];
                    }
                    S[i] = NULL;
                }
                strings_array_free(S, orig_ns);
                spec->n_literal = ns;
            }
        }
    }
}

static void plot_get_ols_info (const char *line,
                               const DATASET *dset,
                               int *reglist)
{
    char fmt[16], vname[VNAMELEN];
    int v;

    sprintf(fmt, "'%%%d[^\\']' (%%d)", VNAMELEN - 1);

    if (sscanf(line + 6, fmt, vname, &v) == 2) {
        if (line[2] == 'X') {
            if (plot_ols_var_ok(dset, vname, v)) {
                reglist[3] = v;
            }
        } else {
            /* 'Y' */
            if (reglist[3] > 0 && plot_ols_var_ok(dset, vname, v)) {
                reglist[0] = 3;
                reglist[1] = v;
            }
        }
    }
}

#define plot_needs_obs(c) (c != PLOT_ELLIPSE && \
                           c != PLOT_PROB_DIST && \
                           c != PLOT_CURVE && \
                           c != PLOT_USER)

/* Read plotspec struct from gnuplot command file.  This is not a
   general parser for gnuplot files; it is designed specifically for
   files auto-generated by gretl and even then there are some sorts of
   plot that are not fully handled.
*/

int plotspec_read_file (GPT_SPEC *spec, const DATASET *dset,
                        PlotReadInfo *info)
{
    linestyle styles[N_GP_LINETYPES];
    int do_markers = 0;
    int auto_linewidth = 0;
    int reglist[4] = {0};
    int *uservec = NULL;
    long barpos = 0;
    long datapos = 0;
    char gpline[MAXLEN];
    GError *gerr = NULL;
    gchar *buf = NULL;
    char *got = NULL;
#if HANDLE_HEREDATA
    char *eod = NULL;
#endif
    int ignore = 0;
    int i, done;
    int err = 0;

#if GPDEBUG
    fprintf(stderr, "plotspec_read_file: spec=%p\n",
            (void *) spec);
#endif

    /* get the content of the plot file */
    if (!g_file_get_contents(spec->fname, &buf, NULL, &gerr)) {
        if (gerr != NULL) {
            gretl_errmsg_set(gerr->message);
            g_error_free(gerr);
        }
        return E_FOPEN;
    }

    spec->n_datacols = 0;
    bufgets_init(buf);

    /* get the number of data-points, plot type, and check for
       observation markers and time-series bars */
    get_plot_nobs(spec, info, buf, &do_markers, &barpos, &datapos);

    if (spec->code == PLOT_GEOMAP) {
        /* FIXME */
        goto bailout;
    }

    if (spec->nobs == 0 && plot_needs_obs(spec->code)) {
        /* failed reading plot data */
#if GPDEBUG
        fprintf(stderr, " got spec->nobs = 0\n");
#endif
        err = 1;
        goto bailout;
    }

    if (do_markers && spec->nobs > MAX_MARKERS) {
        do_markers = 0;
    }

    if (cant_edit(spec->code)) {
        fprintf(stderr, "plotspec_read_file: plot is not editable\n");
        if (maybe_big_multiplot(spec->code)) {
            buf_rewind(buf);
            check_for_plot_size(spec, buf);
        } else if (spec->code == PLOT_BAND) {
            buf_rewind(buf);
            check_for_plot_time_data(spec, buf);
        } else if (spec->code == PLOT_HEATMAP) {
            buf_rewind(buf);
            get_heatmap_matrix(spec, buf, gpline, sizeof gpline, datapos);
        }
        goto bailout;
    } else {
        buf_rewind(buf);
        check_for_plot_size(spec, buf);
    }

    for (i=0; i<N_GP_LINETYPES; i++) {
        linestyle_init(&styles[i]);
    }

    buf_rewind(buf);

    /* get the preamble and "set" lines */

    while ((got = bufgets(gpline, sizeof gpline, buf)) && !err) {
#if GPDEBUG
        tailstrip(gpline);
        fprintf(stderr, "gpline: '%s'\n", gpline);
#endif
        if (!strncmp(gpline, "plot ", 5)) {
            /* we're done with 'set' */
            break;
        }
        if (ignore) {
            if (!strncmp(gpline, "# end inline", 12)) {
                ignore = 0;
            }
        } else if (!strncmp(gpline, "# start inline", 14)) {
            ignore = 1;
        } else if (!strncmp(gpline, "# timeseries", 12)) {
            if (sscanf(gpline, "# timeseries %d", &spec->pd)) {
                info->pd = spec->pd;
            }
            spec->flags |= GPT_TS;
            if (strstr(gpline, "letterbox")) {
                spec->flags |= GPT_LETTERBOX;
            }
        } else if (!strncmp(gpline, "# multiple timeseries", 21)) {
            if (sscanf(gpline, "# multiple timeseries %d", &spec->pd)) {
                info->pd = spec->pd;
            }
            spec->flags |= GPT_TS;
        } else if (!strncmp(gpline, "# scale = ", 10)) {
            double x = dot_atof(gpline + 10);

            if (x >= 0.5 && x <= 2.0) {
                spec->scale = x;
            }
        } else if (!strncmp(gpline, "# auto linewidth", 16)) {
            auto_linewidth = 1;
        } else if (!strncmp(gpline, "# fontspec: ", 12)) {
            free(spec->fontstr);
            spec->fontstr = gretl_strdup(gpline + 12);
            tailstrip(spec->fontstr);
        } else if (!strncmp(gpline, "# boxplots", 10)) {
            ; /* OK */
        } else if (!strncmp(gpline, "# X = ", 6) ||
                   !strncmp(gpline, "# Y = ", 6)) {
            plot_get_ols_info(gpline, dset, reglist);
        } else if (sscanf(gpline, "# literal lines = %d", &spec->n_literal)) {
            spec->literal = strings_array_new(spec->n_literal);
            if (spec->literal == NULL) {
                err = E_ALLOC;
            } else {
                for (i=0; i<spec->n_literal; i++) {
                    if (!bufgets(gpline, MAXLEN - 1, buf)) {
                        gretl_errmsg_set(_("Plot file is corrupted"));
                        err = E_DATA;
                    } else {
                        g_strstrip(gpline);
                        spec->literal[i] = gretl_strdup(gpline);
                    }
                }
            }
        } else if (sscanf(gpline, "# xtics lines = %d", &spec->n_xtics)) {
            spec->multi_xtics = strings_array_new(spec->n_xtics);
            if (spec->multi_xtics == NULL) {
                err = E_ALLOC;
            } else {
                for (i=0; i<spec->n_xtics; i++) {
                    if (!bufgets(gpline, MAXLEN - 1, buf)) {
                        gretl_errmsg_set(_("Plot file is corrupted"));
                        err = E_DATA;
                    } else {
                        g_strstrip(gpline);
                        spec->multi_xtics[i] = gretl_strdup(gpline);
                    }
                }
            }
        } else if (!strncmp(gpline, "# start literal lines", 21) &&
                   spec->literal == NULL) {
            for (i=0; ; i++) {
                if (!bufgets(gpline, MAXLEN - 1, buf)) {
                    gretl_errmsg_set(_("Plot file is corrupted"));
                    err = E_DATA;
                    break;
                } else if (!strncmp(gpline, "# end literal lines", 19)) {
                    break;
                } else {
                    top_n_tail(gpline, 0, NULL);
                    strings_array_add(&spec->literal, &spec->n_literal,
                                      gpline);
                }
            }
        } else if (strstr(gpline, "automatic fit")) {
            spec->flags |= GPT_AUTO_FIT;
            spec->fit = recognize_fit_string(gpline);
        } else  if (strstr(gpline, "user-defined")) {
            uservec = get_user_lines_list(gpline);
        } else if (strstr(gpline, "printing data labels")) {
            spec->flags |= GPT_PRINT_MARKERS;
        } else if (!strncmp(gpline, "# ", 2)) {
            ; /* ignore unknown comment lines */
        } else if (!strncmp(gpline, "set ", 4)) {
            gretl_push_c_numeric_locale();
            err = parse_gp_set_line(spec, gpline, styles, NULL);
            gretl_pop_c_numeric_locale();
        } else if (!strncmp(gpline, "unset ", 6)) {
            err = parse_gp_unset_line(spec, gpline);
        } else {
            /* done reading "set" lines? */
            break;
        }
    }

    if (!err && got == NULL) {
        err = E_DATA;
    }

    if (err) {
        goto bailout;
    }

    for (i=0; i<4; i++) {
        if (!string_is_blank(spec->titles[i])) {
            gretl_delchar('"', spec->titles[i]);
        }
    }

#if HANDLE_HEREDATA
    /* Hmm, what's this doing here? (2020-07-13) */
    if (line_starts_heredata(gpline, &eod)) {
        while (bufgets(gpline, MAXLEN - 1, buf) != NULL) {
            if (!strncmp(gpline, eod, strlen(eod))) {
                bufgets(gpline, MAXLEN - 1, buf);
                break;
            }
        }
        free(eod);
    }
#endif

    /* then get the "plot" lines */
    if (strncmp(gpline, "plot ", 5) ||
        (strlen(gpline) < 10 && bufgets(gpline, MAXLEN - 1, buf) == NULL)) {
        err = unhandled_gp_line_error(gpline);
        goto bailout;
    }

    /* if there are any time-series bars, skip past them */
    for (i=0; i<spec->nbars; i++) {
        if ((got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
            err = E_DATA;
        }
    }

    done = 0;

    while (!err) {
        top_n_tail(gpline, 0, NULL);
        if (!chop_comma(gpline)) {
            /* line did not end with comma -> no continuation of
               the plot command */
            done = 1;
        }
        err = parse_gp_line_line(gpline, spec, auto_linewidth);
        if (err || done || (got = bufgets(gpline, MAXLEN - 1, buf)) == NULL) {
            break;
        }
    }

    if (!err && got == NULL) {
        err = E_DATA;
    }

    if (!err && spec->literal != NULL) {
        maybe_promote_literal_lines(spec, styles);
    }

    /* transcribe styles info into lines for use in the
       GUI plot-editor */

    for (i=0; i<spec->n_lines && !err; i++) {
        GPT_LINE *line = &spec->lines[i];
        int idx = line->type; /* this will be 1-based */

        if (idx == LT_AUTO) {
            /* automatic sequential line-types */
            idx = i + 1;
        }
        if (uservec != NULL && in_gretl_list(uservec, i)) {
            line->flags |= GP_LINE_USER;
        }
        if (idx > 0 && idx < N_GP_LINETYPES) {
            int j = idx - 1;

            if (line->rgb[0] == '\0') {
                strcpy(line->rgb, styles[j].lc);
            }
            if (line->width == 1 && styles[j].lw > 0) {
                line->width = styles[j].lw;
            }
            if (line->dtype == 0 && styles[j].dt > 0) {
                line->dtype = styles[j].dt;
            }
            if (line->ptype == 0 && styles[j].pt > 0) {
                line->ptype = styles[j].pt;
            }
        }
#if LT_DEBUG
        fprintf(stderr, "spec->lines[%d]: type %d idx %d rgb '%s' pt %d lw %g\n",
                i, line->type, idx, line->rgb, line->ptype, (double) line->width);
#endif
        if (spec->auxdata != NULL && i == spec->n_lines - 1) {
            /* the last "line" doesn't use the regular
               data mechanism */
            line->flags = GP_LINE_AUXDATA;
            continue;
        }
    }

#if GPDEBUG
    fprintf(stderr, "plotspec, after line transcription, err=%d\n", err);
#endif

    if (!err) {
        err = plot_get_data_and_markers(spec, buf, &do_markers,
                                        barpos, datapos);
    }

    if (!err && reglist[0] > 0) {
        spec->reglist = gretl_list_copy(reglist);
    }

    if (!err && spec->fit == PLOT_FIT_NONE) {
        maybe_set_add_fit_ok(spec, dset);
    }

    if (!err && spec->code == PLOT_ROOTS && spec->n_literal == 8) {
        /* backward compatibility */
        fprintf(stderr, "old roots plot: fixing\n");
        fix_old_roots_plot(spec);
    }

 bailout:

    if (err) {
        fprintf(stderr, "plotspec_read_file: err = %d\n", err);
    }

    bufgets_finalize(buf);
    g_free(buf);
    free(uservec);

    return err;
}
