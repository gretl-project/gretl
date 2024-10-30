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

/* Sample program to estimate an ARMA model via libgretl */

#include <gretl/libgretl.h>

int arma_estimate (DATASET *dset, PRN *prn)
{
    MODEL *model;
    int *list;
    int err;

    model = gretl_model_new();
    list = gretl_list_new(5);

    list[1] = 1;        /* AR order */
    list[2] = 0;        /* order of integration */
    list[3] = 1;        /* MA order */
    list[4] = LISTSEP;  /* separator */
    list[5] = 1;        /* position of dependent variable in dataset */

    *model = arma(list, NULL, dset, OPT_NONE, prn);
    err = model->errcode;

    if (err) {
        errmsg(err, prn);
    } else {
        printmodel(model, dset, OPT_NONE, prn);
    }

    gretl_model_free(model);
    free(list);

    return err;
}

int main (void)
{
    DATASET *dset;
    PRN *prn;
    int err;

    libgretl_init();

    dset = datainfo_new();
    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL);

    /* Give the full path to the gretl datafile unless it's
       in the current working directory. Note that PREFIX is
       defined in the Makefile
    */
#if defined(_WIN32)
	err = gretl_read_native_data(PREFIX "\\data\\ramanathan\\data9-7.gdt", dset);
#else
	err = gretl_read_native_data(PREFIX "/share/gretl/data/ramanathan/data9-7.gdt", dset);
#endif

    if (err) {
        errmsg(err, prn);
        exit(EXIT_FAILURE);
    }

    err = arma_estimate(dset, prn);

    destroy_dataset(dset);
    gretl_print_destroy(prn);

    libgretl_cleanup();

    return 0;
}
