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

/* Sample program to estimate a model via nonlinear least squares */

#include <gretl/libgretl.h>

int nls_estimate (DATASET *dset, PRN *prn)
{
    MODEL *model = NULL;
    nlspec *spec;
    const char *nls_function;
    char *pnames[] = {"alpha", "beta"}; /* parameter names */
    double parms[] = {-20.0, 3.0};      /* and initial values */
    int err = 0;

    spec = nlspec_new(NLS, dset);
    if (spec == NULL) {
        return E_ALLOC;
    }

    /* specify the function: we're assuming the dataset contains
       variables named x and y
    */
    nls_function = "y = alpha + beta*x1 + (1/beta)*x2";

    err = nlspec_set_regression_function(spec, nls_function, dset);

    if (!err) {
        err = nlspec_add_param_list(spec, 2, parms, pnames);
    }

    if (!err) {
        model = gretl_model_new();
        *model = model_from_nlspec(spec, dset, OPT_NONE, prn);
        err = model->errcode;
    }

    if (!err) {
        printmodel(model, dset, OPT_NONE, prn);
    }

    gretl_model_free(model);
    nlspec_destroy(spec);

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

    /* note: PREFIX is defined in the Makefile */
#if defined(_WIN32)
	err = gretl_read_native_data(PREFIX "\\data\\misc\\ects_nls.gdt", dset);
#else
	err = gretl_read_native_data(PREFIX "/share/gretl/data/misc/ects_nls.gdt",
				 dset);
#endif

    if (err) {
        errmsg(err, prn);
        exit(EXIT_FAILURE);
    }

    err = nls_estimate(dset, prn);
    if (err) {
        errmsg(err, prn);
    }

    destroy_dataset(dset);
    gretl_print_destroy(prn);

    libgretl_cleanup();

    return 0;
}
