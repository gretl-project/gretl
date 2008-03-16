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

#define GRETLCLI /* don't include libxml headers */

#include <gretl/libgretl.h>   

void noalloc(void)
{
    printf("Couldn't allocate memory.\n");
    exit(EXIT_FAILURE);
}

int main (void)
{

    DATAINFO *datainfo;         /* data information struct */
    double **Z;                 /* the data array */
    int *list;                  /* list of regressors etc. */
    MODEL *model;               /* pointer to model struct */
    PRN *prn;                   /* pointer to struct for printing */
    int model_count = 0;        /* keep a tally of models estimated */

    /* basic initialization of library */
    libgretl_init();

    logo(); /* print version info and session time */
    prn = gretl_print_new(GRETL_PRINT_STDOUT, NULL); /* simple printing */

    /* create the datainfo struct and data matrix -- pass in pointer
       to data array; specify the number of variables (allowing one
       for the constant term), and the number of observations on each
       variable (here 3 variables, 14 observations).  The last
       parameter is a 0/1 flag indicating whether we want to supply
       "case marker" strings for the observations: here we don't
    */
    datainfo = create_new_dataset(&Z, 3, 14, 0);
    if (datainfo == NULL) noalloc();

    /* copy in the names of the variables (starting at [1]
       because [0] refers to the constant) */
    strcpy(datainfo->varname[1], "price");
    strcpy(datainfo->varname[2], "sqft");

    /* Fill in the dataset, Z.  Note that Z may be a superset of the 
       data actually used in the regression equation.

       The elements of Z[0] are automatically reserved for the
       constant (via create_new_dataset()), so we fill in data values
       starting from column 1.  
    */

    Z[1][0] = 199.9;	Z[2][0]  = 1065;
    Z[1][1]  = 228;	Z[2][1]  = 1254;
    Z[1][2]  = 235;	Z[2][2]  = 1300;
    Z[1][3]  = 285;	Z[2][3]  = 1577;
    Z[1][4]  = 239;	Z[2][4]  = 1600;
    Z[1][5]  = 293;	Z[2][5]  = 1750;
    Z[1][6]  = 285;	Z[2][6]  = 1800;
    Z[1][7]  = 365;	Z[2][7]  = 1870;
    Z[1][8]  = 295;	Z[2][8]  = 1935;
    Z[1][9]  = 290;	Z[2][9]  = 1948;
    Z[1][10] = 385;	Z[2][10] = 2254;
    Z[1][11] = 505;	Z[2][11] = 2600;
    Z[1][12] = 425;	Z[2][12] = 2800;
    Z[1][13] = 415;	Z[2][13] = 3000;

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
		 &Z,       /* data matrix */
		 datainfo, /* data information */
		 OLS,      /* use Ordinary Least Squares */
		 OPT_NONE  /* no special options */
		 );

    /* Handle case where lsq bombed */
    if (model->errcode) {
        printf("model->errcode: %d\n", model->errcode);
        printf("error message: %s\n", gretl_errmsg_get());
        return 1;
    }

    /* Otherwise give this model an ID number for reference... */
    ++model_count;
    model->ID = model_count;

    /* ...and print info from the regression. */
    printmodel(model, datainfo, OPT_NONE, prn);

    /* memory management check -- try explicitly freeing all allocated
       memory */
    gretl_model_free(model);
    free(list);
    destroy_dataset(Z, datainfo); 
    gretl_print_destroy(prn);

    libgretl_cleanup();

    return 0;

}

