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

typedef struct RGB_ {
    guint8 r; /* red channel [0,255] */
    guint8 g; /* green channel [0,255] */
    guint8 b; /* blue channel [0,255] */
} RGB;

union channels {
    guint32 u32;
    guint8 u8[4];
};

static int show_mixed_colors (gretlRGB c1, gretlRGB c2,
			      guint32 *mix, const double *f,
			      int n);

/* rgb color mixing: (1-f) * c1 + f * c2 */

static RGB rgb_mix (RGB *c1, RGB *c2, double f)
{
    RGB ret = {0};
    double x;

    x = (1-f) * c1->r + f * c2->r;
    ret.r = (guint8) nearbyint(x);
    x = (1-f) * c1->g + f * c2->g;
    ret.g = (guint8) nearbyint(x);
    x = (1-f) * c1->b + f * c2->b;
    ret.b = (guint8) nearbyint(x);

    return ret;
}

static RGB RGB_from_guint32 (guint32 u)
{
    union channels ch;
    RGB out;

    ch.u32 = u;
    out.r = ch.u8[2];
    out.g = ch.u8[1];
    out.b = ch.u8[0];

    return out;
}

gretl_array *colormix_array (gretlRGB c1, gretlRGB c2,
			     const double *f, int nf,
			     int do_plot, int *err)
{
    gretl_array *ret = NULL;
    guint32 *uvec = NULL;
    char color[9] = {0};
    RGB rgb1, rgb2, mix;
    guint32 u;
    int i;

    ret = gretl_array_new(GRETL_TYPE_STRINGS, nf, err);
    if (*err) {
	return NULL;
    }

    if (do_plot) {
	uvec = malloc(nf * sizeof *uvec);
    }

    rgb1 = RGB_from_guint32(c1);
    rgb2 = RGB_from_guint32(c2);

    for (i=0; i<nf; i++) {
	mix = rgb_mix(&rgb1, &rgb2, f[i]);
	u = (mix.r << 16) | (mix.g << 8) | mix.b;
	sprintf(color, "0x%06x", u);
	gretl_array_set_string(ret, i, color, 1);
	if (do_plot) {
	    uvec[i] = u;
	}
    }

    if (*err == 0 && do_plot) {
	show_mixed_colors(c1, c2, uvec, f, nf);
    }

    free(uvec);

    return ret;
}

static int show_mixed_colors (gretlRGB c1, gretlRGB c2,
			      guint32 *mix, const double *f,
			      int n)
{
    FILE *fp;
    int i, j;
    int err = 0;

    set_optval_string(GNUPLOT, OPT_U, "display");
    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "set title \"mix #%06x and #%06x via RGB\\n", c1, c2);
    if (n == 1) {
	fprintf(fp, "f = %g\"\n", f[0]);
	fprintf(fp, "set xtics (\"#%06x\" 1, \"#%06x\" 2, \"#%06x\" 3)\n",
		c1, c2, mix[0]);
    } else {
	fprintf(fp, "f = ");
	for (i=0; i<n; i++) {
	    fprintf(fp, "%g", f[i]);
	    fputs(i < n-1 ? ", " : "\"\n", fp);
	}
	fputs("set xtics (", fp);
	for (i=0; i<n; i++) {
	    fprintf(fp, "\"#%06x\" %d", mix[i], i+1);
	    fputs(i < n-1 ? ", " : ")\n", fp);
	}
    }
    fputs("set border 0\n", fp);
    fputs("set nokey\n", fp);
    fputs("set noytics\n", fp);
    fputs("set boxwidth 0.9 relative\n", fp);
    fputs("set style fill solid 1.0\n", fp);
    if (n == 1) {
	/* show c1, c2, mix */
	fprintf(fp, "set linetype 1 lc rgb \"#%06x\"\n", c1);
	fprintf(fp, "set linetype 2 lc rgb \"#%06x\"\n", c2);
	fprintf(fp, "set linetype 3 lc rgb \"#%06x\"\n", mix[0]);
	n = 3;
    } else {
	/* just show the mixes */
	for (i=0; i<n; i++) {
	    fprintf(fp, "set linetype %d lc rgb \"#%06x\"\n", i+1, mix[i]);
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
