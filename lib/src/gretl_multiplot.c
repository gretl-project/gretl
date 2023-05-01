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

/* Support for the built-in "gridplot" command for producing
   multiple plots.
*/

#include "libgretl.h"
#include "gretl_multiplot.h"

typedef struct {
    int row;
    int col;
    gchar *buf;
} mplot_element;

typedef struct {
    gretlopt flag;
    int *target;
    int min;
    int max;
} mplot_option;

static GArray *multiplot;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

static void set_multiplot_defaults (void)
{
    mp_fontsize = 10;
    mp_width = 800;
    mp_height = 600;
    mp_rows = 0;
    mp_cols = 0;
}

static void mplot_element_clear (mplot_element *element)
{
    g_free(element->buf);
}

int gretl_multiplot_active (void)
{
    return multiplot != NULL && mp_collecting;
}

void gretl_multiplot_destroy (void)
{
    if (multiplot != NULL) {
        g_array_unref(multiplot);
        multiplot = NULL;
    }
}

static const mplot_option mp_options[] = {
    { OPT_F, &mp_fontsize, 4, 24 },
    { OPT_W, &mp_width,    300, 2048 },
    { OPT_H, &mp_height,   300, 2048 },
    { OPT_R, &mp_rows,     1, 10 },
    { OPT_C, &mp_cols,     1, 10 }
};

static int set_multiplot_sizes (gretlopt opt)
{
    const mplot_option *mpo;
    int i, k, err = 0;

    for (i=0; i<G_N_ELEMENTS(mp_options) && !err; i++) {
	mpo = &mp_options[i];
	if (opt & mpo->flag) {
	    k = get_optval_int(GRIDPLOT, mpo->flag, &err);
	    if (!err && (k < mpo->min || k > mpo->max)) {
		gretl_errmsg_set("gridplot: bad option value");
		err = E_INVARG;
	    }
	    if (!err) {
		*mpo->target = k;
	    }
	}
    }

    return err;
}

int gretl_multiplot_start (gretlopt opt)
{
    int err = 0;

    if (multiplot == NULL) {
	if (opt == OPT_NONE) {
	    set_multiplot_defaults();
	} else {
	    err = set_multiplot_sizes(opt);
	}
	if (!err) {
	    multiplot = g_array_new(FALSE, FALSE, sizeof(mplot_element));
	    g_array_set_clear_func(multiplot, (GDestroyNotify) mplot_element_clear);
	    mp_collecting = 1;
	}
    } else {
        gretl_errmsg_set("multplot: cannot be nested");
        err = E_DATA;
    }

    return err;
}

int gretl_multiplot_add_plot (int row, int col, gchar *buf)
{
    if (multiplot != NULL && buf != NULL) {
        mplot_element element = {row, col, buf};

        g_array_append_val(multiplot, element);
        return 0;
    } else {
	gretl_errmsg_set("gretl_multiplot_add_plot: failed");
        return E_DATA;
    }
}

static int multiplot_set_grid (int n, int *pr, int *pc)
{
    int nr = *pr;
    int nc = *pc;
    int err = 0;

    if (nr == 0 && nc == 0) {
	/* fully automatic grid */
	*pr = ceil(sqrt((double) n));
	*pc = ceil((double) n / *pr);
    } else if (nr == 0) {
	/* automatic rows */
	*pr = ceil((double) n / nc);
    } else if (nc == 0) {
	/* automatic cols */
	*pc = ceil((double) n / nr);
    } else if (nr * nc < n) {
	gretl_errmsg_sprintf("Specified grid (%d by %d) is too small "
			     "for %d sub-plots", nr, nc, n);
	err = E_INVARG;
    } else if (nr * nc > n) {
	int ar = ceil(sqrt((double) n));
	int ac = ceil((double) n / ar);

	if (nr * nc > ar * ac) {
	    gretl_errmsg_sprintf("Specified grid (%d by %d) is too big "
				 "for %d sub-plots", nr, nc, n);
	    err = E_INVARG;
	}
    }

    return err;
}

int gretl_multiplot_finalize (gretlopt opt)
{
    int np, err = 0;
    int rows = 0;
    int cols = 0;

    if (multiplot == NULL) {
	gretl_errmsg_set("end multiplot: multiplot not started");
	return E_DATA;
    }

    mp_collecting = 0;
    np = multiplot->len;
    if (np > 0) {
	err = multiplot_set_grid(np, &rows, &cols);
    }

    if (np > 0 && !err) {
	mplot_element *element;
	FILE *fp = NULL;
	int i;

	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
	fp = open_plot_input_file(PLOT_USER_MULTI, 0, &err);

	if (!err) {
	    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", rows, cols);
	    gretl_push_c_numeric_locale();
	    for (i=0; i<multiplot->len; i++) {
		fprintf(fp, "# multiplot: subplot %d\n", i);
		if (i > 0) {
		    fputs("reset\n", fp);
		}
		element = &g_array_index(multiplot, mplot_element, i);
		fputs(element->buf, fp);
	    }
	    gretl_pop_c_numeric_locale();
	    fputs("unset multiplot\n", fp);
	    err = finalize_plot_input_file(fp);
	}
    }

    gretl_multiplot_destroy();

    return err;
}
