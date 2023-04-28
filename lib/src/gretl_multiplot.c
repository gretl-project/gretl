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
// #include "gretl_multiplot.h"

typedef struct {
    int row;
    int col;
    gchar *buf;
} mplot_element;

static GArray *multiplot;
static int mp_sizes[3] = {10, 800, 600};
static int mp_collecting;

static void set_multiplot_defaults (void)
{
    mp_sizes[0] = 10;  /* font size */
    mp_sizes[1] = 800; /* width */
    mp_sizes[2] = 600; /* height */
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

static int set_multiplot_sizes (gretlopt opt)
{
    int k, err = 0;

    if (opt & OPT_F) {
	k = get_optval_int(MULTIPLT, OPT_F, &err);
	if (!err && (k < 6 || k > 20)) {
	    gretl_errmsg_set("multiplt: bad fontsize value");
	    err = E_INVARG;
	}
	if (!err) {
	    mp_sizes[0] = k;
	}
    }
    if (!err && (opt & OPT_W)) {
	k = get_optval_int(MULTIPLT, OPT_W, &err);
	if (!err && (k < 400 || k > 2048)) {
	    gretl_errmsg_set("multiplt: bad width value");
	    err = E_INVARG;
	}
	if (!err) {
	    mp_sizes[1] = k;
	}
    }
    if (!err && (opt & OPT_H)) {
	k = get_optval_int(MULTIPLT, OPT_H, &err);
	if (!err && (k < 300 || k > 2048)) {
	    gretl_errmsg_set("multiplt: bad width value");
	    err = E_INVARG;
	}
	if (!err) {
	    mp_sizes[2] = k;
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
	    fprintf(stderr, "gretl_multiplot_start: OK\n");
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
	fprintf(stderr, "gretl_multiplot_add_plot: OK\n");
        return 0;
    } else {
	fprintf(stderr, "gretl_multiplot_add_plot: failed\n");
        return E_DATA;
    }
}

int gretl_multiplot_finalize (gretlopt opt)
{
    mp_collecting = 0;

    if (multiplot != NULL) {
	const char *s = get_optval_string(MULTIPLT, OPT_U);
	mplot_element *element;
	int i;

	if (s != NULL) {
	    fprintf(stderr, "HERE gretl_multiplot_finalize, output='%s'\n", s);
	}

	/* FIXME reference mp_sizes */

        fprintf(stderr, "HERE gretl_multiplot_finalize, OK\n");
	for (i=0; i<multiplot->len; i++) {
	    element = &g_array_index(multiplot, mplot_element, i);
	    fprintf(stderr, "plot %d: r %d, c %d, bufstart '%.32s'\n",
		    i, element->row, element->col, element->buf);
	}
        gretl_multiplot_destroy();
        return 0;
    } else {
        fprintf(stderr, "HERE gretl_multiplot_finalize, no plots\n");
        return E_DATA;
    }
}
