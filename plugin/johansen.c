/*
 *  Copyright (c) by Allin Cottrell
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

#include "libgretl.h"
#include "gretl_matrix.h"

/* Critical values for Johansen's likelihood ratio test of the
   null hypothesis of h cointegrating relations against the
   alternative of no restrictions.
*/

const double johansen_pvals_a[15][3] = {
    /* .1      .05      .01  = p-value */ 
    /* case 1 (no constant included in auxiliary regressions) */
    {  2.86,   3.84,   6.51}, /* g = 1 (# of random walks) */
    { 10.47,  12.53,  16.31}, /* g = 2 */
    { 21.63,  24.31,  29.75}, /* g = 3 */
    { 36.58,  39.89,  45.58}, /* g = 4 */
    { 55.44,  59.46,  55.52}, /* g = 5 */
    /* case 2 (no deterministic trends) */
    { 6.691,  8.083, 11.576}, /* g = 1 */
    {15.583, 17.844, 21.962}, /* g = 2 */
    {28.436, 31.256, 37.291}, /* g = 3 */
    {45.248, 48.419, 55.551}, /* g = 4 */
    {65.956, 69.977, 77.911}, /* g = 5 */
    /* case 3 (one or more vars has deterministic trend) */
    { 2.816,  3.962,  6.936}, /* g = 1 */
    {13.338, 15.197, 19.310}, /* g = 2 */
    {26.791, 29.509, 35.397}, /* g = 3 */
    {43.964, 47.181, 53.792}, /* g = 4 */
    {65.063, 68.905, 76.955}  /* g = 5 */
};

/* Critical values for Johansen's likelihood ratio test of the
   null hypothesis of h cointegrating relations against the
   alternative of h + 1 relations.
*/

const double johansen_pvals_b[15][3] = {
    /* .1      .05      .01  = p-value */ 
    /* case 1 (no constant included in auxiliary regressions) */
    {  2.86,   3.84,   6.51}, /* g = 1 (# of random walks) */
    {  9.52,  11.44,  15.69}, /* g = 2 */
    { 15.59,  17.89,  22.99}, /* g = 3 */
    { 21.58,  23.80,  28.82}, /* g = 4 */
    { 27.62,  30.04,  35.17}, /* g = 5 */
    /* case 2 (no deterministic trends) */
    { 6.691,  8.083, 11.576}, /* g = 1 */
    {12.783, 14.595, 18.782}, /* g = 2 */
    {18.959, 21.279, 26.154}, /* g = 3 */
    {24.917, 27.341, 32.616}, /* g = 4 */
    {30.818, 33.262, 38.858}, /* g = 5 */
    /* case 3 (one or more vars has deterministic trend) */
    { 2.816,  3.962,  6.936}, /* g = 1 */
    {12.099, 14.036, 17.936}, /* g = 2 */
    {18.697, 20.778, 25.521}, /* g = 3 */
    {24.712, 27.169, 31.943}, /* g = 4 */
    {30.774, 33.178, 38.341}  /* g = 5 */
};

static int inverse_compare_doubles (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da < *db) - (*da > *db);
}

int johansen_eigenvals (const double **X, const double **Y, const double **Z, 
			int k, double *evals)
{
    gretl_matrix *Suu, *Svv, *Suv;
    gretl_matrix *Inv, *TmpL, *TmpR, *M;
    double *eigvals;
    int err = 0;

    Suu = gretl_matrix_from_2d_array(X, k, k);
    Svv = gretl_matrix_from_2d_array(Y, k, k);
    Suv = gretl_matrix_from_2d_array(Z, k, k);

    Inv = gretl_matrix_alloc(k, k);
    TmpL = gretl_matrix_alloc(k, k);
    TmpR = gretl_matrix_alloc(k, k);
    M = gretl_matrix_alloc(k, k);

    /* calculate Suu^{-1} Suv */
    gretl_invert_general_matrix(Suu);
    gretl_matmult(Suu, Suv, TmpR);

    /* calculate Svv^{-1} Suv' */
    gretl_invert_general_matrix(Svv);
    gretl_matmult_mod(Svv, GRETL_MOD_NONE,
		      Suv, GRETL_MOD_TRANSPOSE, 
		      TmpL);

    gretl_matmult(TmpL, TmpR, M);

    eigvals = gretl_general_matrix_eigenvals(M);

    if (eigvals != NULL) {
	int i;

	for (i=0; i<k; i++) {
	    evals[i] = eigvals[i];
	}
	free(eigvals);
    }

    qsort(evals, k, sizeof *evals, inverse_compare_doubles);

    /* free stuff */
    gretl_matrix_free(Svv);
    gretl_matrix_free(Suu);
    gretl_matrix_free(Suv);

    gretl_matrix_free(Inv);
    gretl_matrix_free(TmpL);
    gretl_matrix_free(TmpR);
    gretl_matrix_free(M);

    return err;
}
