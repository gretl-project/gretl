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

static void mplot_element_clear (mplot_element *element)
{
    g_free(element->buf);
}

void gretl_multiplot_destroy (void)
{
    if (multiplot != NULL) {
        g_array_unref(multiplot);
        multiplot = NULL;
    }
}

int gretl_multiplot_start (char *param /* unused at present */)
{
    if (multiplot == NULL) {
        multiplot = g_array_new(FALSE, FALSE, sizeof(mplot_element));
        g_array_set_clear_func(multiplot, (GDestroyNotify) mplot_element_clear);
        fprintf(stderr, "gretl_multiplot_start: OK\n");
        return 0;
    } else {
        gretl_errmsg_set("multplot: cannot be nested");
        return E_DATA;
    }
}

int gretl_multiplot_add_plot (int row, int col, gchar *buf)
{
    if (multiplot != NULL) {
        mplot_element element = {row, col, buf};

        g_array_append_val(multiplot, element);
        return 0;
    } else {
        return E_DATA;
    }
}

int gretl_multiplot_get_plot (int i, int *prow, int *pcol,
                              gchar **pbuf)
{
    if (multiplot != NULL && i < multiplot->len) {
        mplot_element *element;

        element = &g_array_index(multiplot, mplot_element, i);
        *prow = element->row;
        *pcol = element->col;
        *pbuf = element->buf;
        return 0;
    } else {
        return E_DATA;
    }    
}

int gretl_multiplot_finalize (gretlopt opt)
{
    if (multiplot != NULL) {
        fprintf(stderr, "HERE gretl_multiplot_finalize, OK\n");
        gretl_multiplot_destroy();
        return 0;
    } else {
        fprintf(stderr, "HERE gretl_multiplot_finalize, no plots\n");
        return E_DATA;
    }    
}



