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

/* plotbands.c for gretl : support multiple "bands" in a given plot */

#include "libgretl.h"
#include "usermat.h"
#include "uservar.h"
#include "gretl_typemap.h"
#include "plot_priv.h"

#define PB_DEBUG 0
#define SUPPORT_LEGACY 1

typedef enum {
    BAND_LINE,
    BAND_FILL,
    BAND_DASH,
    BAND_BARS,
    BAND_STEP
} BandStyle;

static const char *style_specs[] = {
    "line", "fill", "dash", "bars", "step"
};

typedef struct {
    int center;      /* series ID of center-of-band data */
    int width;       /* series ID of width-of-band data */
    double factor;   /* factor by which to multiply width */
    int bdummy;      /* flag for drawing rectangles */
    BandStyle style; /* see enumeration above */
    char rgb[10];    /* specific color, if wanted */
    const char *title;
} band_info;

#if SUPPORT_LEGACY
static band_info **legacy_get_band_info (const char *spec,
                                         int matrix_mode,
					 gnuplot_info *gi,
					 DATASET *dset,
					 gretlopt opt,
					 int *err);
#endif

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
	bi->title = NULL;
    }

    return bi;
}

static void free_bbi (band_info **bbi, int n)
{
    if (bbi != NULL) {
	int i;

	for (i=0; i<n; i++) {
	    free(bbi[i]);
	}
	free(bbi);
    }
}

static BandStyle style_from_string (const char *s, int *err)
{
    int i, n = G_N_ELEMENTS(style_specs);

    for (i=0; i<n; i++) {
        if (!strcmp(s, style_specs[i])) {
            return (BandStyle) i;
        }
    }

    *err = invalid_field_error(s);
    return (BandStyle) 0;
}

/* handle the band-as-matrix case */

static int process_band_matrix (const gretl_matrix *m,
				band_info *bi,
				gnuplot_info *gi,
				DATASET *dset)
{
    int err = 0;

    if (m == NULL || m->cols != 2 || m->rows != dset->n) {
	/* matrix missing or non-conformable */
	err = E_NONCONF;
    } else {
        /* enlarge the dset->Z array and stick bi->center and
           bi->width onto the end of it
        */
	err = matrix_dataset_expand_Z(dset, 2);
        if (!err) {
            bi->center = dset->v - 2;
            bi->width = dset->v - 1;
            dset->Z[bi->center] = m->val;
            dset->Z[bi->width] = m->val + m->rows;
        }
    }

    if (!err) {
        gretl_list_append_term(&gi->list, bi->center);
	bi->center = gi->list[0] + 1;
        gretl_list_append_term(&gi->list, bi->width);
	bi->width = gi->list[0] + 1;
    }

    return err;
}

/* check the validity of a series that's supposed to be a dummy
   variable (to produce something resembling recession bars)
*/

static int check_assign_bdummy (band_info *bi, int vnum,
				const DATASET *dset)
{
    int err = 0;

    if (gretl_isdummy(dset->t1, dset->t2, dset->Z[vnum])) {
	bi->bdummy = vnum;
    } else {
	err = E_INVARG;
	fprintf(stderr, "%s: not a dummy variable\n",
		dset->varname[vnum]);
    }

    return err;
}

static void do_center_or_width (band_info *bi,
				gnuplot_info *gi,
				int i, int v)
{
    int pos = in_gretl_list(gi->list, v);

    if (pos == 0) {
	/* not already in list */
	gretl_list_append_term(&gi->list, v);
	pos = gi->list[0];
    }

    if (i == 0) {
	bi->center = pos + 1;
    } else {
	bi->width = pos + 1;
    }
}

static int has_bogus_key (gretl_bundle *b)
{
    const char *known_keys[] = {
	"center", "width", "dummy", "style", "color",
	"factor", "bandmat", "title", NULL
    };
    char **S = NULL;
    int i, j, ns = 0;
    int ok = 0;

    S = gretl_bundle_get_keys_raw(b, &ns);

    for (i=0; i<ns; i++) {
	ok = 0;
	for (j=0; known_keys[j] != NULL; j++) {
	    if (!strcmp(S[i], known_keys[j])) {
		ok = 1;
		break;
	    }
	}
	if (!ok) {
	    gretl_errmsg_sprintf("Unknown key '%s'", S[i]);
	    break;
	}
    }

    strings_array_free(S, ns);

    return !ok;
}

static int get_input_color (gretl_bundle *b, char *rgb)
{
    GretlType t = gretl_bundle_get_member_type(b, "color", NULL);
    guint32 u;
    int err = 0;

    if (t == GRETL_TYPE_STRING) {
	const char *s = gretl_bundle_get_string(b, "color", &err);

	if (!err) {
	    if (!strcmp(s, "default")) {
		/* ensure @rgb is blank */
		rgb[0] = '\0';
		return 0;
	    } else {
		u = numeric_color_from_string(s, &err);
	    }
	}
    } else if (gretl_is_scalar_type(t)) {
	u = gretl_bundle_get_unsigned(b, "color", &err);
    } else {
	err = E_TYPES;
    }

    if (!err) {
	sprintf(rgb, "#%x", u);
    }

    return err;
}

static const gretl_matrix *retrieve_bandmat (gretl_bundle *b,
                                             int *err)
{
    gretl_matrix *m = NULL;
    GretlType t;
    void *ptr;

    ptr = gretl_bundle_get_data(b, "bandmat", &t, NULL, err);

    if (!*err) {
        if (t == GRETL_TYPE_MATRIX) {
            m = (gretl_matrix *) ptr;
        } else if (t == GRETL_TYPE_STRING) {
            m = get_matrix_by_name((const char *) ptr);
        } else {
            *err = E_TYPES;
        }
    }

    return m;
}

static band_info *band_info_from_bundle (int matrix_mode,
					 gretl_bundle *b,
					 gnuplot_info *gi,
                                         DATASET *dset,
                                         int *err)
{
    band_info *bi = band_info_new();
    const char *S[5] = {NULL};
    const char *strkeys[] = {
        "center", "width", "dummy", "style", "title"
    };
    char colspec[32];
    int v, i, imin = 0;

    if (has_bogus_key(b)) {
	*err = E_DATA;
	return NULL;
    }

    /* Allow for the "matrix_mode" case where the band (center, width)
       is represented by a two-column matrix. This arises when
       "gnuplot" or "plot" is fed data in matrix form.
    */
    if (matrix_mode) {
        const gretl_matrix *m = retrieve_bandmat(b, err);

	if (!*err) {
	    *err = process_band_matrix(m, bi, gi, dset);
	}
	if (!*err) {
	    imin = 3;
	} else {
	    return NULL;
	}
    }

    *colspec = '\0';

    for (i=imin; i<5 && !*err; i++) {
        if (gretl_bundle_has_key(b, strkeys[i])) {
            S[i] = gretl_bundle_get_string(b, strkeys[i], err);
        }
    }
    if (!*err && gretl_bundle_has_key(b, "color")) {
	*err = get_input_color(b, bi->rgb);
    }
    if (!*err && gretl_bundle_has_key(b, "factor")) {
        bi->factor = gretl_bundle_get_scalar(b, "factor", err);
    }

    /* try cashing out the string members */
    for (i=imin; i<5 && !*err; i++) {
        if (S[i] == NULL) {
            continue;
        }
        if (i < 3) {
            v = current_series_index(dset, S[i]);
            if (v < 0) {
                *err = invalid_field_error(S[i]);
            } else if (i < 2) {
		do_center_or_width(bi, gi, i, v);
	    } else {
		/* dummy? */
		*err = check_assign_bdummy(bi, v, dset);
	    }
        } else if (i == 3) {
            bi->style = style_from_string(S[i], err);
        } else if (i == 4) {
	    bi->title = S[i];
	}
    }

    if (!*err) {
        /* further checks on validity of input */
        if (bi->bdummy > 0) {
            /* not compatible with width and center */
            if (bi->center > 0 || bi->width > 0) {
                *err = E_BADOPT;
            }
        } else if (bi->center <= 0) {
            /* without bdummy, center is required */
            *err = E_ARGS;
        }
        if (!*err && (bi->factor <= 0 || na(bi->factor))) {
            *err = E_INVARG;
        }
    }

    if (*err) {
        free(bi);
        bi = NULL;
    }

    return bi;
}

static int band_matrices_differ (const gretl_matrix *a,
                                 const gretl_matrix *b)
{
    const void *p0 = a;
    const void *p1 = b;

    if (p1 - p0 == 0) {
        return 0;
    } else if (a->rows * a->cols != b->rows * b->cols) {
        return 1;
    } else {
        int i, n = a->rows * a->cols;

        for (i=0; i<n; i++) {
            if (a->val[i] != b->val[i]) {
                return 1;
            }
        }
    }

    return 0;
}

/* Here we're trying to determine a plausible z-order for two bands
   that both employ the "fill" style, to avoid having a narrower band
   obscured by a wider one. The idea is that if a multiplicative
   "factor" is present in both cases, and the two bands have the same
   center and width, they can be ordered by the factor. Otherwise we
   treat them as unordered.
*/

static int swap_bands_order (gretl_array *a, int *err)
{
    gretl_bundle *b0 = gretl_array_get_data(a, 0);
    gretl_bundle *b1 = gretl_array_get_data(a, 1);
    const char *s0;
    const char *s1;
    double f0, f1;

    f0 = gretl_bundle_get_scalar(b0, "factor", NULL);
    f1 = gretl_bundle_get_scalar(b1, "factor", NULL);
    if (na(f0) || na(f1)) {
        /* we won't try to order the two */
        return 0;
    }

    s0 = gretl_bundle_get_string(b0, "style", NULL);
    s1 = gretl_bundle_get_string(b1, "style", NULL);
    if (s0 == NULL || s1 == NULL) {
        return 0;
    } else if (strcmp(s0, "fill") || strcmp(s1, "fill")) {
        /* condition ordering on "fill" for both bands? */
        return 0;
    }

    if (gretl_bundle_has_key(b0, "bandmat") &&
        gretl_bundle_has_key(b1, "bandmat")) {
        const gretl_matrix *m0 = retrieve_bandmat(b0, err);
        const gretl_matrix *m1 = retrieve_bandmat(b1, err);

        if (m0 == NULL || m1 == NULL ||
            band_matrices_differ(m0, m1)) {
            return 0;
        }
    } else if (gretl_bundle_has_key(b0, "center") &&
               gretl_bundle_has_key(b1, "center") &&
               gretl_bundle_has_key(b0, "width") &&
               gretl_bundle_has_key(b1, "width")) {
        s0 = gretl_bundle_get_string(b0, "center", NULL);
        s1 = gretl_bundle_get_string(b1, "center", NULL);
        if (s0 == NULL || s1 == NULL || strcmp(s0, s1)) {
            /* we don't have the same centers */
            return 0;
        }
        s0 = gretl_bundle_get_string(b0, "width", NULL);
        s1 = gretl_bundle_get_string(b1, "width", NULL);
        if (s0 == NULL || s1 == NULL || strcmp(s0, s1)) {
            /* we don't have the same widths */
            return 0;
        }
    }

    return f1 > f0;
}

/* Here we're supposing we got --bands=@aname. We try to retrieve an
   array of bundles, and if that succeeds we build and return an array
   of band_info structs.
*/

static band_info **get_band_info_array (int matrix_mode,
					int *n_bands,
					gnuplot_info *gi,
                                        DATASET *dset,
                                        int *err)
{
    const char *s = get_optval_string(GNUPLOT, OPT_a);
    gretl_array *a = get_array_by_name(s);
    band_info **pbi = NULL;
    int n = 0;
    int swap = 0;

    if ((n = gretl_array_get_length(a)) < 1) {
	fprintf(stderr, "get_band_info_array: array has no content\n");
        *err = E_DATA;
        return NULL;
    } else if (gretl_array_get_type(a) != GRETL_TYPE_BUNDLES) {
	fprintf(stderr, "get_band_info_array: array has wrong type\n");
        *err = E_TYPES;
        return NULL;
    }

    if (n == 2) {
        swap = swap_bands_order(a, err);
    }

    pbi = calloc(n, sizeof *pbi);
    if (pbi == NULL) {
        *err = E_ALLOC;
    } else {
	gretl_bundle *b;
        int i, j;

        for (i=0; i<n && !*err; i++) {
            b = gretl_array_get_data(a, swap ? 1 - i : i);
            pbi[i] = band_info_from_bundle(matrix_mode, b, gi, dset, err);
        }
        if (*err) {
            for (j=0; j<=i; j++) {
                free(pbi[j]);
            }
            free(pbi);
            pbi = NULL;
        } else {
            *n_bands = n;
        }
    }

    return pbi;
}

static const char *get_lw_string (void)
{
    if (get_default_png_scale() > 1.0) {
	return " lw 2";
    } else {
	return "";
    }
}

static void gp_newline (int contd, FILE *fp)
{
    if (contd) {
	fputs(", \\\n", fp);
    } else {
	fputc('\n', fp);
    }
}

static void band_title (band_info *bi, FILE *fp)
{
    if (bi->title != NULL) {
	fprintf(fp, "title '%s' ", bi->title);
    } else {
	fputs("notitle ", fp);
    }
}

static void print_pm_filledcurve (band_info *bi, int contd,
				  FILE *fp)
{
    double f = bi->factor;
    int c = bi->center;
    int w = bi->width;
    char cstr[10];

    if (bi->rgb[0] != '\0') {
        strcpy(cstr, bi->rgb);
    } else {
        print_rgb_hash(cstr, get_shadecolor());
    }

    if (bi->factor == 1.0) {
	fprintf(fp, "$data using 1:($%d-$%d):($%d+$%d) ", c, w, c, w);
    } else {
	fprintf(fp, "$data using 1:($%d-%g*$%d):($%d+%g*$%d) ",
		c, f, w, c, f, w);
    }
    band_title(bi, fp);
    fprintf(fp, "lc rgb \"%s\" w filledcurve", cstr);
    gp_newline(contd, fp);
}

/* for plain lines, dashed lines, or steps */

static void print_pm_lines (band_info *bi, int n_yvars,
			    const char *lw, int contd,
			    FILE *fp)
{
    char *wstr = bi->style == BAND_STEP ? "steps" : "lines";
    char lspec[24], dspec[8] = {0};
    double f = bi->factor;
    int c = bi->center;
    int w = bi->width;

    if (bi->rgb[0] != '\0') {
	sprintf(lspec, "lc rgb \"%s\"", bi->rgb);
    } else {
	/* default to gray */
	strcpy(lspec, "lc rgb \"#bbbbbb\"");
    }
    if (bi->style == BAND_DASH) {
	strcpy(dspec, " dt 2");
    }
    /* lower line */
    fprintf(fp, "$data using 1:($%d-%g*$%d) ", c, f, w);
    band_title(bi, fp);
    fprintf(fp, "w %s %s%s%s, \\\n", wstr, lspec, dspec, lw);
    /* upper line */
    fprintf(fp, "$data using 1:($%d+%g*$%d) notitle w %s %s%s%s",
	    c, f, w, wstr, lspec, dspec, lw);
    gp_newline(contd, fp);
}

/* for errorbars only */

static void print_pm_bars (band_info *bi, int n_yvars,
			   int contd, FILE *fp)
{
    char lspec[24];

    /* note: "pt 7" is a solid circle */

    if (bi->rgb[0] != '\0') {
	sprintf(lspec, "lc rgb \"%s\" pt 7", bi->rgb);
    } else {
	sprintf(lspec, "lt %d pt 7", n_yvars + 1);
    }
    fprintf(fp, "$data using 1:%d:(%g*$%d) ", bi->center, bi->factor, bi->width);
    band_title(bi, fp);
    fprintf(fp, "w errorbars %s", lspec);
    gp_newline(contd, fp);
}

/* Print the data in the form of a gnuplot "heredoc": this should be
   relatively efficient, avoiding duplication of the x-axis variable
   and the band center.
*/

static void print_data_block (gnuplot_info *gi,
			      const double *x,
			      const DATASET *dset,
			      int *show_zero,
			      FILE *fp)
{
    double yt;
    int gt_zero = 0;
    int lt_zero = 0;
    int n = gi->list[0];
    int t, i;

    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (t=gi->t1; t<=gi->t2; t++) {
	fprintf(fp, "%g", x[t]);
	for (i=0; i<n; i++) {
	    yt = dset->Z[gi->list[i+1]][t];
	    if (na(yt)) {
		fprintf(fp, " %s", GPNA);
	    } else {
		fprintf(fp, " %.10g", yt);
	    }
	    if (yt > 0) {
		gt_zero++;
	    } else if (yt < 0) {
		lt_zero++;
	    }
	}
	fputc('\n', fp);
    }
    fputs("EOD\n", fp);
    fputs("# end inline data\n", fp);

    if (gt_zero && lt_zero) {
	*show_zero = 1;
    }
}

static band_info **get_single_band_info (int matrix_mode,
					 gnuplot_info *gi,
					 DATASET *dset,
					 gretlopt opt,
					 int *err)
{
    const char *spec = get_optval_string(GNUPLOT, OPT_N);
    gretl_bundle *b = NULL;
    band_info **bi = NULL;

    if (spec == NULL || *spec == '\0') {
        *err = E_INVARG;
	return NULL;
    }

    b = get_bundle_by_name(spec);
    if (b != NULL) {
        bi = malloc(sizeof *bi);
        bi[0] = band_info_from_bundle(matrix_mode, b, gi, dset, err);
        if (*err) {
            free(bi);
            bi = NULL;
        }
        return bi; /* done */
    }

#if SUPPORT_LEGACY
    bi = legacy_get_band_info(spec, matrix_mode, gi, dset, opt, err);
    gretl_warnmsg_set(_("Deprecated --band syntax: please use the current syntax\n"
                        "as described in the help for the 'gnuplot' command"));
#else
    *err = E_INVARG;
#endif

    return bi;
}

/* Write full-height bars as gnuplot rectangle objects, using
   the dummy variable bi->bdummy for on/off information.
*/

static int write_rectangles (band_info *bi,
			     gnuplot_info *gi,
			     const double *x,
                             DATASET *dset,
                             FILE *fp)
{
    double *d = dset->Z[bi->bdummy];
    char stobs[16], endobs[16];
    int bar_on = 0, obj = 1;
    int t, err = 0;

    if (bi->rgb[0] == '\0') {
        strcpy(bi->rgb, "#dddddd");
    }

    *stobs = *endobs = '\0';

    for (t=gi->t1; t<=gi->t2; t++) {
        if (na(d[t])) {
            err = E_MISSDATA;
            break;
        }
        if (bar_on && d[t] == 0) {
            /* finalize a bar */
            sprintf(endobs, "%g", x[t]);
            fprintf(fp, "set object %d rectangle from %s, graph 0 to %s, graph 1 back "
                    "fillstyle solid 0.5 noborder fc rgb \"%s\"\n",
                    obj++, stobs, endobs, bi->rgb);
            bar_on = 0;
        } else if (!bar_on && d[t] != 0) {
            /* start a bar */
            sprintf(stobs, "%g", x[t]);
            bar_on = 1;
        }
    }

    if (bar_on) {
        /* terminate an unfinished bar */
        sprintf(endobs, "%g", x[gi->t2]);
        fprintf(fp, "set object rectangle from %s, graph 0 to %s, graph 1 back "
                "fillstyle solid 0.5 noborder fc rgb \"%s\"\n",
                stobs, endobs, bi->rgb);
    }

    return err;
}

/* Below: the bars in question don't have to indicate recessions,
   but they're the same sort of thing: full-height bars indicating
   the presence of absence of some binary time-varying attribute.
   This function handles the case where the only "band" instance
   is a case of such bars.
*/

static void recession_bars_plot (band_info *bi,
				 gnuplot_info *gi,
				 int n_yvars,
				 const double *x,
				 DATASET *dset,
				 gretlopt opt,
				 FILE *fp)
{
    const double *y;
    char wspec[16] = {0};
    int oddman = 0;
    int i;

    /* write out the rectangles as objects */
    write_rectangles(bi, gi, x, dset, fp);

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
	const char *iname = plotname(dset, gi->list[i], 1);

	set_plot_withstr(gi, i, wspec);
	if (gi->flags & GPT_Y2AXIS) {
	    fprintf(fp, "'-' using 1:2 axes %s title \"%s (%s)\" %s lt %d",
		    (i == oddman)? "x1y2" : "x1y1", iname,
		    (i == oddman)? _("right") : _("left"),
		    wspec, i);
	} else {
	    fprintf(fp, "'-' using 1:2 title \"%s\" %s lt %d", iname, wspec, i);
	}
	gp_newline(i < n_yvars, fp);
    }

    /* and write the data block */
    for (i=0; i<n_yvars; i++) {
	y = dset->Z[gi->list[i+1]];
	print_user_y_data(x, y, gi->t1, gi->t2, fp);
    }
}

/* The last entry in gi->list pertains to the x-axis variable:
   either the index of a "genuine" x-var or a dummy placeholder
   of 0. Either way we want to drop this term from the list,
   and return a pointer to the actual x-axis data.
*/

static const double *gi_get_xdata (gnuplot_info *gi,
				   char *xname,
				   DATASET *dset)
{
    int xpos = gi->list[0];
    const double *x;

    if (gi->x != NULL) {
	/* the "plotx" (dummy) case */
	x = gi->x;
	*xname = '\0';
        if (gi->flags & GPT_TS) {
            gi->flags |= GPT_LETTERBOX;
        }
    } else {
	int xno = gi->list[xpos];

	x = dset->Z[xno];
	strcpy(xname, plotname(dset, xno, 1));
    }

    gretl_list_delete_at_pos(gi->list, xpos);

    return x;
}

/* If we have a single data series to plot and no band is
   titled, omit the key and put the name of the series on
   the y axis.
*/

static int maybe_suppress_key (band_info **bb, int n_bands,
			       char *yname, const char *src,
			       FILE *fp)
{
    int i, have_titles = 0;

    for (i=0; i<n_bands; i++) {
	if (bb[i]->title != NULL) {
	    have_titles = 1;
	    break;
	}
    }
    if (!have_titles) {
	strcpy(yname, src);
	fprintf(fp, "set ylabel \"%s\"\n", yname);
	return 1;
    } else {
	return 0;
    }
}

int plot_with_band (BPMode mode,
		    gnuplot_info *gi,
                    const char *literal,
                    DATASET *dset,
                    gretlopt opt)
{
    band_info **bbi = NULL;
    band_info *bi = NULL;
    FILE *fp = NULL;
    const double *x = NULL;
    const char *lwstr;
    char yname[MAXDISP];
    char xname[MAXDISP];
    char wspec[16] = {0};
    int show_zero = 0;
    int no_key = 0;
    int i, j, n_yvars = 0;
    int matrix_mode;
    int n_bands = 1;
    int err = 0;

    matrix_mode = (opt & OPT_X) || mode == BP_BLOCKMAT;

#if PB_DEBUG
    fprintf(stderr, "\nplot_with_band: matrix_mode = %d\n", matrix_mode);
    printlist(gi->list, "gi->list");
#endif

    /* subtract 1 for x */
    n_yvars = gi->list[0] - 1;

    x = gi_get_xdata(gi, xname, dset);

#if PB_DEBUG
    printlist(gi->list, "gi->list, rev 1");
#endif

    if (opt & OPT_a) {
	/* --bands=@array */
	bbi = get_band_info_array(matrix_mode, &n_bands, gi, dset, &err);
    } else {
	/* --band=<whatever> */
	bbi = get_single_band_info(matrix_mode, gi, dset, opt, &err);
    }
    if (err) {
	return err;
    }

#if PB_DEBUG
    printlist(gi->list, "gi->list, rev 2");
#endif

    err = graph_list_adjust_sample(gi, dset, 1);
    if (err) {
        free_bbi(bbi, n_bands);
        return err;
    }

    fp = open_plot_input_file(PLOT_BAND, gi->flags, &err);
    if (err) {
        free_bbi(bbi, n_bands);
        return err;
    }

    if (gi->flags & GPT_TS) {
        PRN *prn = gretl_print_new_with_stream(fp);

        make_time_tics(gi, dset, 0, NULL, prn);
        gretl_print_detach_stream(prn);
        gretl_print_destroy(prn);
    }

    if (n_yvars == 1) {
	const char *s = plotname(dset, gi->list[1], 1);

	no_key = maybe_suppress_key(bbi, n_bands, yname, s, fp);
    }

    if (no_key) {
	fputs("set nokey\n", fp);
    } else {
	fputs("set key left top\n", fp);
    }

    if (*xname != '\0') {
        fprintf(fp, "set xlabel \"%s\"\n", xname);
    }

    /* for adjusting line width if wanted */
    lwstr = get_lw_string();

    gretl_push_c_numeric_locale();

    if (gi->x != NULL) {
        print_x_range(gi, fp);
    }

    print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fp);

    if (show_zero && bi->style != BAND_FILL) {
	/* in the "fill" case xzeroaxis won't be visible */
        fputs("set xzeroaxis\n", fp);
    }

    if (bbi[0]->bdummy) {
	/* special case */
	if (n_bands == 1) {
	    recession_bars_plot(bbi[0], gi, n_yvars,
				x, dset, opt, fp);
	    goto finish;
	} else {
	    err = write_rectangles(bbi[0], gi, x, dset, fp);
	}
    }

    print_data_block(gi, x, dset, &show_zero, fp);

    fputs("plot \\\n", fp);

    /* bands first */
    for (j=0; j<n_bands; j++) {
	bi = bbi[j];
	if (bi->bdummy) {
	    continue;
	} else if (bi->style == BAND_FILL) {
	    print_pm_filledcurve(bi, 1, fp);
	    if (show_zero) {
		/* custom xzeroaxis on top of the fill */
		fputs("0 notitle w lines lc rgb \"#888888\" dt 2, \\\n", fp);
	    }
	} else if (bi->style == BAND_BARS) {
	    print_pm_bars(bi, n_yvars, 1, fp);
	} else {
	    print_pm_lines(bi, n_yvars, lwstr, 1, fp);
	}
    }
    /* then the non-band data */
    for (i=1; i<=n_yvars; i++) {
	const char *iname = plotname(dset, gi->list[i], 1);

	set_plot_withstr(gi, i, wspec);
	fprintf(fp, "$data using 1:%d title '%s' %s lt %d%s", i+1,
		iname, wspec, i, lwstr);
	gp_newline(i < n_yvars, fp);
    }

 finish:

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    clear_gpinfo(gi);
    free_bbi(bbi, n_bands);

    if (matrix_mode) {
        /* hide the two extra dataset columns
           representing the band */
        dset->v -= 2;
    }

    return err;
}

#ifdef SUPPORT_LEGACY
# include "pb_legacy.c"
#endif
