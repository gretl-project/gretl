/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* Trivial sample client program for libesl */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gretl/libgretl.h>   

void noalloc(void)
{
    printf("Couldn't allocate memory.\n");
    exit(EXIT_FAILURE);
}

int main (void)
{

    DATAINFO *datainfo;         /* data information struct */
    double *Z;                  /* the data set */
    int *list;                  /* list of regressors etc. */
    MODEL *model;               /* pointer to model struct */
    print_t prn;                /* struct for printing */
    int model_count = 0;        /* keep a tally of models estimated */

    logo(); /* print version info and session time */
    prn.fp = stdout; /* simple printing */

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

       Elements 0 to n-1 in Z are automatically reserved for the
       constant (via create_new_dataset()), so we fill in data values
       starting from position n.  
       
       The first column below gives 14 (that is, n) values for the
       first variable in the data set, the second column supplies
       14 observations on the second variable...
    */

    Z[14] = 199.9;	Z[28] = 1065;
    Z[15] = 228;	Z[29] = 1254;
    Z[16] = 235;	Z[30] = 1300;
    Z[17] = 285;	Z[31] = 1577;
    Z[18] = 239;	Z[32] = 1600;
    Z[19] = 293;	Z[33] = 1750;
    Z[20] = 285;	Z[34] = 1800;
    Z[21] = 365;	Z[35] = 1870;
    Z[22] = 295;	Z[36] = 1935;
    Z[23] = 290;	Z[37] = 1948;
    Z[24] = 385;	Z[38] = 2254;
    Z[25] = 505;	Z[39] = 2600;
    Z[26] = 425;	Z[40] = 2800;
    Z[27] = 415;	Z[41] = 3000;

    /* Set up the "list", which is fed to the regression function.
       The first element of list represents the length of the list
       vector itself, counting from zero.  The second entry is the ID
       number of the dependent variable (i.e. its place in the data
       set Z) counting from one (zero being reserved for the
       constant).  The third entry (and there can be more) is the ID
       number of the first independent variable.  "list" should be
       malloc'ed: it will be realloc'ed by libesl.  
    */
    list = malloc(4 * sizeof *list);
    if (list == NULL) noalloc(); 
    list[0] = 3;   /* three variables follow */
    list[1] = 1;   /* the dependent variable is the one with ID# 1 */
    list[2] = 0;   /* we include a constant (ID# 0) */
    list[3] = 2;   /* the independent variable has ID# 2 */

    /* Now we call the lsq function from libesl to get least squares 
       estimates and associated statistics. */
    model = gretl_model_new();
    if (model == NULL) noalloc();
    *model = lsq(list, Z, datainfo, OLS, 1, 0.0);

    /* Handle case where lsq bombed */
    if (model->errcode) {
        printf("model->errcode: %d\n", model->errcode);
        printf("model->errmsg: %s\n", model->errmsg);
        return 1;
    }

    /* Otherwise give this model an ID number for reference... */
    ++model_count;
    model->ID = model_count;

    /* ...and print info from the regression. */
    printmodel(model, datainfo, &prn);

    /* memory management check -- try explicitly freeing all allocated
       memory */
    free(Z); 
    free_model(model);
    free(list);
    free_datainfo(datainfo);

    return 0;

}

