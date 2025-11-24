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
#include "libset.h"

typedef struct ARGB_ {
    guint8 a; /* alpha channel [0,255] */
    guint8 r; /* red channel [0,255] */
    guint8 g; /* green channel [0,255] */
    guint8 b; /* blue channel [0,255] */
} ARGB;

union channels {
    guint32 u32;
    guint8 u8[4];
};

static int show_mixed_colors (gretlRGB c1, gretlRGB c2,
			      guint32 *mix, const double *f,
			      int n);

/* rgb color mixing: (1-f) * c1 + f * c2 */

static ARGB argb_mix (ARGB *c1, ARGB *c2, double f)
{
    ARGB ret = {0};
    double x;

    x = (1-f) * c1->a + f * c2->a;
    ret.a = (guint8) nearbyint(x);
    x = (1-f) * c1->r + f * c2->r;
    ret.r = (guint8) nearbyint(x);
    x = (1-f) * c1->g + f * c2->g;
    ret.g = (guint8) nearbyint(x);
    x = (1-f) * c1->b + f * c2->b;
    ret.b = (guint8) nearbyint(x);

    return ret;
}

static ARGB ARGB_from_guint32 (guint32 u)
{
    union channels ch;
    ARGB out;

    ch.u32 = u;
    out.a = ch.u8[3];
    out.r = ch.u8[2];
    out.g = ch.u8[1];
    out.b = ch.u8[0];

    return out;
}

void decompose_argb (guint32 u,
		     guint8 *a,
		     guint8 *r,
		     guint8 *g,
		     guint8 *b)
{
    union channels ch;

    ch.u32 = u;
    *a = ch.u8[3];
    *r = ch.u8[2];
    *g = ch.u8[1];
    *b = ch.u8[0];
}

gretl_array *colormix_array (gretlRGB c0, gretlRGB c1,
			     const double *f, int nf,
			     int do_plot, int *err)
{
    gretl_array *ret = NULL;
    guint32 *uvec = NULL;
    char color[12] = {0};
    ARGB argb0, argb1, mix;
    guint32 u;
    int i;

    ret = gretl_array_new(GRETL_TYPE_STRINGS, nf, err);
    if (*err) {
	return NULL;
    }

    if (do_plot) {
	uvec = malloc(nf * sizeof *uvec);
    }

    argb0 = ARGB_from_guint32(c0);
    argb1 = ARGB_from_guint32(c1);

    for (i=0; i<nf; i++) {
	mix = argb_mix(&argb0, &argb1, f[i]);
	if (mix.a > 0) {
	    u = (mix.a << 24) | (mix.r << 16) | (mix.g << 8) | mix.b;
	    sprintf(color, "0x%08x", u);
	} else {
	    u = (mix.r << 16) | (mix.g << 8) | mix.b;
	    sprintf(color, "0x%06x", u);
	}
	gretl_array_set_string(ret, i, color, 1);
	if (uvec != NULL) {
	    uvec[i] = u;
	}
    }

    if (*err == 0 && uvec != NULL) {
	show_mixed_colors(c0, c1, uvec, f, nf);
    }

    free(uvec);

    return ret;
}

static int show_mixed_colors (gretlRGB c0, gretlRGB c1,
			      guint32 *mix, const double *f,
			      int n)
{
    FILE *fp;
    char cs[3][10];
    int i, j;
    int err = 0;

    set_optval_string(GNUPLOT, OPT_U, "display");
    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	return err;
    }

    print_rgb_hash(cs[0], c0);
    print_rgb_hash(cs[1], c1);
    print_rgb_hash(cs[2], mix[0]);

    fprintf(fp, "set title \"mix %s and %s via RGB\\n", cs[0], cs[1]);
    if (n == 1) {
	fprintf(fp, "f = %g\"\n", f[0]);
	fprintf(fp, "set xtics (\"%s\" 1, \"%s\" 2, \"%s\" 3)\n",
		cs[0], cs[1], cs[2]);
    } else {
	fprintf(fp, "f = ");
	for (i=0; i<n; i++) {
	    fprintf(fp, "%g", f[i]);
	    fputs(i < n-1 ? ", " : "\"\n", fp);
	}
	fputs("set xtics (", fp);
	for (i=0; i<n; i++) {
	    print_rgb_hash(cs[2], mix[i]);
	    fprintf(fp, "\"%s\" %d", cs[2], i+1);
	    fputs(i < n-1 ? ", " : ")\n", fp);
	}
    }
    fputs("set border 0\n", fp);
    fputs("set nokey\n", fp);
    fputs("set noytics\n", fp);
    fputs("set boxwidth 0.9 relative\n", fp);
    fputs("set style fill solid 1.0\n", fp);
    if (n == 1) {
	/* show c0, c1, mix */
	fprintf(fp, "set linetype 1 lc rgb \"%s\"\n", cs[0]);
	fprintf(fp, "set linetype 2 lc rgb \"%s\"\n", cs[1]);
	fprintf(fp, "set linetype 3 lc rgb \"%s\"\n", cs[2]);
	n = 3;
    } else {
	/* just show the mixes */
	for (i=0; i<n; i++) {
	    print_rgb_hash(cs[2], mix[i]);
	    fprintf(fp, "set linetype %d lc rgb \"%s\"\n", i+1, cs[2]);
	}
    }
    fputs("$data << EOD\n", fp);
    for (i=0; i<n; i++) {
	fprintf(fp, "%d ", i+1);
	for (j=0; j<n; j++) {
	    fprintf(fp, "%d", j == i ? 1 : 0);
	    fputc(j < n-1 ? ' ' : '\n', fp);
	}
    }
    fputs("EOD\n", fp);
    fputs("plot \\\n", fp);
    for (i=0; i<n; i++) {
	fprintf(fp, " '$data' using 1:%d w boxes", i+2);
	fputs(i < n-1 ? ", \\\n" : "\n", fp);
    }

    err = finalize_plot_input_file(fp);
    if (!err && gretl_in_gui_mode()) {
	manufacture_gui_callback(GNUPLOT);
    }

    return err;
}

