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

/* Trivial sample client program for libgretl */

#include <gretl/libgretl.h>

void noalloc (void)
{
    printf("Couldn't allocate memory.\n");
    exit(EXIT_FAILURE);
}

int main (void)
{

    DATASET *dataset;       /* pointer to dataset struct */
    int *list;              /* list of regressors etc. */
    MODEL *model;           /* pointer to model struct */
    PRN *prn;               /* pointer to struct for printing */
    int model_count = 0;    /* keep a tally of models estimated */
    int i;                  /* index variable */

    /* the first data series to load */
    double z1[] = {
	199.9, 228, 235, 285, 239, 293, 285,
	365, 295, 290, 385, 505, 425, 415
    };
    /* and the second series */
    double z2[] = {
	1065, 1254, 1300, 1577, 1600, 1750, 1800,
	1870, 1935, 1948, 2254, 2600, 2800, 3000
    };

    /* basic initialization of library */
    libgretl_init();

    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL); /* simple printing */

    /* allocate the dataset struct: the parameters are the number of
       variables (here 3, allowing for the constant in position 0),
       the number of observations on each variable (here 14), and
       a 0/1 flag indicating whether we want to supply "case marker"
       strings for the observations (here we don't).
    */
    dataset = create_new_dataset(3, 14, 0);
    if (dataset == NULL) noalloc();

    /* copy in the names of the variables (starting at [1]
       because [0] refers to the constant) */
    strcpy(dataset->varname[1], "price");
    strcpy(dataset->varname[2], "sqft");

    /* Fill in the data array, starting at variable 1. Note that
       this array may be a superset of the data actually used in
       the regression equation. Note that dataset->n records the
       number of observations.
    */

    for (i=0; i<dataset->n; i++) {
        dset_set_data(dataset, 1, i, z1[i]);
        dset_set_data(dataset, 2, i, z2[i]);
    }

    /* Set up the "list", which is fed to the regression function.
       The first element of list represents the length of the list
       vector itself, counting from zero.  The second entry is the ID
       number of the dependent variable (i.e. its place in the data
       set Z) counting from one (zero being reserved for the
       constant).  The third entry (and there can be more) is the ID
       number of the first independent variable.
    */

    list = gretl_list_new(3); /* number of terms will be 3 */
    list[1] = 1;   /* the dependent variable is the one with ID# 1 */
    list[2] = 0;   /* we include a constant (ID# 0) */
    list[3] = 2;   /* the independent variable has ID# 2 */

    /* Now we call the lsq function from libgretl to get least squares
       estimates and associated statistics. */
    model = gretl_model_new();
    if (model == NULL) noalloc();
    *model = lsq(list,     /* regressand and regressors */
		 dataset,  /* the dataset */
		 OLS,      /* use Ordinary Least Squares */
		 OPT_NONE  /* no special options */
		 );

    /* Handle case where lsq bombed */
    if (model->errcode) {
        printf("model->errcode: %d\n", model->errcode);
        printf("error message: %s\n", gretl_errmsg_get());
        return 1;
    }

    /* Otherwise give this model an ID number for reference */
    model->ID = ++model_count;

    /* and print the regression results */
    printmodel(model, dataset, OPT_NONE, prn);

    /* memory management check -- try explicitly freeing all allocated
       memory */
    gretl_model_free(model);
    free(list);
    destroy_dataset(dataset);
    gretl_print_destroy(prn);

    libgretl_cleanup();

    return 0;
}
