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

static const int rank_sum_lower[39][3] = {
  /* (nA, nB) 0.01, 0.05, 0.1 */
    { -1, 11, 13 },
    { 10, 12, 14 },
    { 11, 13, 15 },
    { 11, 14, 16 },
    { 12, 15, 17 },
    { 13, 16, 19 },
    { 13, 17, 20 },
    { 14, 18, 21 },
    { 15, 19, 22 },
    { 16, 19, 20 },
    { 17, 20, 22 },
    { 18, 21, 23 },
    { 19, 23, 25 },
    { 20, 24, 27 },
    { 21, 26, 28 },
    { 22, 27, 30 },
    { 23, 28, 32 },
    { 24, 28, 30 },
    { 25, 29, 32 },
    { 27, 31, 34 },
    { 28, 33, 36 },
    { 29, 35, 38 },
    { 30, 37, 40 },
    { 32, 38, 42 },
    { 34, 39, 41 },
    { 35, 41, 44 },
    { 37, 43, 46 },
    { 39, 45, 49 },
    { 40, 47, 51 },
    { 42, 49, 54 },
    { 45, 51, 55 },
    { 47, 54, 58 },
    { 49, 56, 60 },
    { 51, 59, 63 },
    { 53, 62, 66 },
    { 59, 66, 70 },
    { 61, 69, 73 },
    { 63, 72, 76 },
    { 66, 75, 80 }
};

static const int rank_sum_upper[39][3] = {
    /* (nA, nB) 0.1, 0.05, 0.01 */
    { 23, 25, -1 },
    { 26, 28, 30 },
    { 29, 31, 33 },
    { 32, 34, 37 },
    { 35, 37, 40 },
    { 37, 40, 43 },
    { 40, 43, 47 },
    { 43, 46, 50 },
    { 46, 49, 53 },
    { 35, 36, 39 },
    { 38, 40, 43 },
    { 42, 44, 47 },
    { 45, 47, 51 },
    { 48, 51, 55 },
    { 52, 54, 59 },
    { 55, 58, 63 },
    { 58, 62, 67 },
    { 48, 50, 54 },
    { 52, 55, 59 },
    { 56, 59, 63 },
    { 60, 63, 68 },
    { 64, 67, 73 },
    { 68, 71, 78 },
    { 72, 76, 82 },
    { 64, 66, 71 },
    { 68, 71, 77 },
    { 73, 76, 82 },
    { 77, 81, 87 },
    { 82, 86, 93 },
    { 86, 91, 98 },
    { 81, 85, 91 },
    { 86, 90, 97 },
    { 92, 96, 103 },
    { 97, 101, 109 },
    { 102, 106, 115 },
    { 101, 105, 112 },
    { 107, 111, 119 },
    { 113, 117, 126 },
    { 118, 123, 132 }
};

/* D-W lookup apparatus: thanks to Marcin Blazejowski and Tadeusz Kufel
   See also http://www.stanford.edu/~clint/bench/dwcrit.htm 
*/

/* dw_row: returns the row of the data table on which to find the 
   appropriate DW values.

   This table runs from n = 6 to n = 2000:
   - all values are represented from 6 to 200 
   - n goes by 10s from 200 to 500
   - n goes by 50s from 500 to 2000
 */

static int dw_row (int *pn)
{
    int rem, pos, row = 254;
    int n = *pn;

    if (n <= 200) {
	row = n - 6;
    } else if (n <= 500) {
	pos = (n - 200) / 10 + 194;
	rem = n % 10;
	row = (rem > 5)? pos + 1 : pos;
	n = (n/10) * 10 + ((rem > 5)? 10 : 0);
    } else if (n < 2000) {
	pos = (n - 500) / 50 + 194 + 30;
	rem = n % 50;
	row = (rem > 25)? pos + 1 : pos;
	n = (n/50) * 50 + ((rem > 25)? 50 : 0);
    }

    *pn = n;

    return row;
}

int dw_lookup (int n, int k, gretl_matrix **pm)
{
    gzFile fz;
    char datfile[FILENAME_MAX];
    double dl = 0, du = 0;
    int dn = n, dk = k;
    char buf[14];
    int r, c;
    int err = 0;

    if (n < 6) {
	gretl_errmsg_set("DW: n must be at least 6");
	return E_DATA;
    }
	
    /* Open data file */
#ifdef WIN32
    sprintf(datfile, "%splugins\\data\\dwdata.gz", gretl_lib_path());
#else
    sprintf(datfile, "%sdata/dwdata.gz", gretl_lib_path());
#endif
    fz = gretl_gzopen(datfile, "rb");
    if (fz == NULL) {
	gretl_errmsg_set("Couldn't open D-W table");
	return E_FOPEN;
    }

    if (dk > 20) dk = 20;
    if (dn > 2000) dn = 2000;

    r = dw_row(&dn);
    c = 14 * (dk - 1);

    gzseek(fz, (z_off_t) (r * 280 + c), SEEK_SET);
    gzgets(fz, buf, sizeof buf);

    gretl_push_c_numeric_locale();
    sscanf(buf, "%lf %lf", &dl, &du);
    gretl_pop_c_numeric_locale();

    gzclose(fz);
		
    if (dl == 0 || du == 0) {
	gretl_errmsg_sprintf("No critical values available for n=%d and k=%d\n", n, k);
	err = E_DATA;
    } else {
	gretl_vector *v = gretl_vector_alloc(4);

	if (v == NULL) {
	    err = E_ALLOC;
	} else {
	    /* fill vector with dl, du, and effective n, k */
	    gretl_vector_set(v, 0, dl);
	    gretl_vector_set(v, 1, du);
	    gretl_vector_set(v, 2, (double) dn);
	    gretl_vector_set(v, 3, (double) dk);
	    *pm = v;
	}
    }
	
    return err;
}

static int rank_table_row (int na, int nb)
{
    int step = 9, ret = 0;

    if (na < 4 || na > 9 || nb < na || nb > 12) {
	return -1;
    }

    nb -= na;
    na -= 4;

    while (na-- > 0) {
	ret += step--;
    }

    while (nb-- > 0) {
	ret++;
    }

    return ret;
}

void rank_sum_lookup (int na, int nb, PRN *prn)
{
    int i = rank_table_row(na, nb);

    if (i < 0) {
	return;
    }

    pprintf(prn, "\n%s:\n", _("Critical values"));

    if (i > 0) {
	pprintf(prn, "  %s: %2d%% %d, %2d%% %d, %2d%% %d\n", _("lower tail"),
		1,  rank_sum_lower[i][0], 
		5,  rank_sum_lower[i][1], 
		10, rank_sum_lower[i][2]);
	pprintf(prn, "  %s: %2d%% %d, %2d%% %d, %2d%% %d\n", _("upper tail"),
		10, rank_sum_upper[i][0], 
		5,  rank_sum_upper[i][1], 
		1,  rank_sum_upper[i][2]);
    } else {
	pprintf(prn, "  %s: %2d%% %d, %2d%% %d\n", _("lower tail"),
		5,  rank_sum_lower[i][1], 
		10, rank_sum_lower[i][2]);
	pprintf(prn, "  %s: %2d%% %d, %2d%% %d\n",_("upper tail"),
		10, rank_sum_upper[i][1], 
		5,  rank_sum_upper[i][2]);
    }	
}

/* 
   Stock and Yogo, 2003, Table 1.
   Critical values for weak instrument test based on TSLS relative bias

   cols 1 - 4: n = 1; b = {0.05, 0.10, 0.20, 0.30}
   cols 5 - 8: n = 2; b = {0.05, 0.10, 0.20, 0.30}
   cols 9 -12: n = 3; b = {0.05, 0.10, 0.20, 0.30}

   where b is desired maximal bias relative to OLS. 
   Rows: K2 values, 3 to 30. 
*/

double tsls_bias_vals[][12] = {
    { 13.91,  9.08, 6.46, 5.39,     0,     0,    0,    0,     0,     0,    0,    0 },
    { 16.85, 10.27, 6.71, 5.34, 11.04,  7.56, 5.57, 4.73,     0,     0,    0,    0 },
    { 18.37, 10.83, 6.77, 5.25, 13.97,  8.78, 5.91, 4.79,  9.53,  6.61, 4.99, 4.30 },
    { 19.28, 11.12, 6.76, 5.15, 15.72,  9.48, 6.08, 4.78, 12.20,  7.77, 5.35, 4.40 },
    { 19.86, 11.29, 6.73, 5.07, 16.88,  9.92, 6.16, 4.76, 13.95,  8.50, 5.56, 4.44 },
    { 20.25, 11.39, 6.69, 4.99, 17.70, 10.22, 6.20, 4.73, 15.18,  9.01, 5.69, 4.46 },
    { 20.53, 11.46, 6.65, 4.92, 18.30, 10.43, 6.22, 4.69, 16.10,  9.37, 5.78, 4.46 },
    { 20.74, 11.49, 6.61, 4.86, 18.76, 10.58, 6.23, 4.66, 16.80,  9.64, 5.83, 4.45 },
    { 20.90, 11.51, 6.56, 4.80, 19.12, 10.69, 6.23, 4.62, 17.35,  9.85, 5.87, 4.44 },
    { 21.01, 11.52, 6.53, 4.75, 19.40, 10.78, 6.22, 4.59, 17.80, 10.01, 5.90, 4.42 },
    { 21.10, 11.52, 6.49, 4.71, 19.64, 10.84, 6.21, 4.56, 18.17, 10.14, 5.92, 4.41 },
    { 21.18, 11.52, 6.45, 4.67, 19.83, 10.89, 6.20, 4.53, 18.47, 10.25, 5.93, 4.39 },
    { 21.23, 11.51, 6.42, 4.63, 19.98, 10.93, 6.19, 4.50, 18.73, 10.33, 5.94, 4.37 },
    { 21.28, 11.50, 6.39, 4.59, 20.12, 10.96, 6.17, 4.48, 18.94, 10.41, 5.94, 4.36 },
    { 21.31, 11.49, 6.36, 4.56, 20.23, 10.99, 6.16, 4.45, 19.13, 10.47, 5.94, 4.34 },
    { 21.34, 11.48, 6.33, 4.53, 20.33, 11.00, 6.14, 4.43, 19.29, 10.52, 5.94, 4.32 },
    { 21.36, 11.46, 6.31, 4.51, 20.41, 11.02, 6.13, 4.41, 19.44, 10.56, 5.94, 4.31 },
    { 21.38, 11.45, 6.28, 4.48, 20.48, 11.03, 6.11, 4.39, 19.56, 10.60, 5.93, 4.29 },
    { 21.39, 11.44, 6.26, 4.46, 20.54, 11.04, 6.10, 4.37, 19.67, 10.63, 5.93, 4.28 },
    { 21.40, 11.42, 6.24, 4.43, 20.60, 11.05, 6.08, 4.35, 19.77, 10.65, 5.92, 4.27 },
    { 21.41, 11.41, 6.22, 4.41, 20.65, 11.05, 6.07, 4.33, 19.86, 10.68, 5.92, 4.25 },
    { 21.41, 11.40, 6.20, 4.39, 20.69, 11.05, 6.06, 4.32, 19.94, 10.70, 5.91, 4.24 },
    { 21.42, 11.38, 6.18, 4.37, 20.73, 11.06, 6.05, 4.30, 20.01, 10.71, 5.90, 4.23 },
    { 21.42, 11.37, 6.16, 4.35, 20.76, 11.06, 6.03, 4.29, 20.07, 10.73, 5.90, 4.21 },
    { 21.42, 11.36, 6.14, 4.34, 20.79, 11.06, 6.02, 4.27, 20.13, 10.74, 5.89, 4.20 },
    { 21.42, 11.34, 6.13, 4.32, 20.82, 11.05, 6.01, 4.26, 20.18, 10.75, 5.88, 4.19 },
    { 21.42, 11.33, 6.11, 4.31, 20.84, 11.05, 6.00, 4.24, 20.23, 10.76, 5.88, 4.18 },
    { 21.42, 11.32, 6.09, 4.29, 20.86, 11.05, 5.99, 4.23, 20.27, 10.77, 5.87, 4.17 }
};

/* 
   Stock and Yogo, 2003, Table 2.
   Critical values for weak instrument test based on TSLS size,
   for nominal significance level 5%.

   cols 1-4: n = 1; r = {0.10, 0.15, 0.20, 0.25}
   cols 5-8: n = 2; r = {0.10, 0.15, 0.20, 0.25}

   where r is desired maximal size. Rows: K2 values, 1 to 30. 
*/

double tsls_size_vals[][8] = {
    { 16.38,  8.96,  6.66,  5.53,     0,     0,     0,     0 },
    { 19.93, 11.59,  8.75,  7.25,  7.03,  4.58,  3.95,  3.63 },
    { 22.30, 12.83,  9.54,  7.80, 13.43,  8.18,  6.40,  5.45 },
    { 24.58, 13.96, 10.26,  8.31, 16.87,  9.93,  7.54,  6.28 },
    { 26.87, 15.09, 10.98,  8.84, 19.45, 11.22,  8.38,  6.89 },
    { 29.18, 16.23, 11.72,  9.38, 21.68, 12.33,  9.10,  7.42 },
    { 31.50, 17.38, 12.48,  9.93, 23.72, 13.34,  9.77,  7.91 },
    { 33.84, 18.54, 13.24, 10.50, 25.64, 14.31, 10.41,  8.39 },
    { 36.19, 19.71, 14.01, 11.07, 27.51, 15.24, 11.03,  8.85 },
    { 38.54, 20.88, 14.78, 11.65, 29.32, 16.16, 11.65,  9.31 },
    { 40.90, 22.06, 15.56, 12.23, 31.11, 17.06, 12.25,  9.77 },
    { 43.27, 23.24, 16.35, 12.82, 32.88, 17.95, 12.86, 10.22 },
    { 45.64, 24.42, 17.14, 13.41, 34.62, 18.84, 13.45, 10.68 },
    { 48.01, 25.61, 17.93, 14.00, 36.36, 19.72, 14.05, 11.13 },
    { 50.39, 26.80, 18.72, 14.60, 38.08, 20.60, 14.65, 11.58 },
    { 52.77, 27.99, 19.51, 15.19, 39.80, 21.48, 15.24, 12.03 },
    { 55.15, 29.19, 20.31, 15.79, 41.51, 22.35, 15.83, 12.49 },
    { 57.53, 30.38, 21.10, 16.39, 43.22, 23.22, 16.42, 12.94 },
    { 59.92, 31.58, 21.90, 16.99, 44.92, 24.09, 17.02, 13.39 },
    { 62.30, 32.77, 22.70, 17.60, 46.62, 24.96, 17.61, 13.84 },
    { 64.69, 33.97, 23.50, 18.20, 48.31, 25.82, 18.20, 14.29 },
    { 67.07, 35.17, 24.30, 18.80, 50.01, 26.69, 18.79, 14.74 },
    { 69.46, 36.37, 25.10, 19.41, 51.70, 27.56, 19.38, 15.19 },
    { 71.85, 37.57, 25.90, 20.01, 53.39, 28.42, 19.97, 15.64 },
    { 74.24, 38.77, 26.71, 20.61, 55.07, 29.29, 20.56, 16.10 },
    { 76.62, 39.97, 27.51, 21.22, 56.76, 30.15, 21.15, 16.55 },
    { 79.01, 41.17, 28.31, 21.83, 58.45, 31.02, 21.74, 17.00 },
    { 81.40, 42.37, 29.12, 22.43, 60.13, 31.88, 22.33, 17.45 },
    { 83.79, 43.57, 29.92, 23.04, 61.82, 32.74, 22.92, 17.90 },
    { 86.17, 44.78, 30.72, 23.65, 63.51, 33.61, 23.51, 18.35 }
};

/* 
   Stock and Yogo, 2003, Table 3.
   Critical values for weak instrument test based on LIML size,
   for nominal significance level 5%.

   cols 1-4: n = 1; r = {0.10, 0.15, 0.20, 0.25}
   cols 5-8: n = 2; r = {0.10, 0.15, 0.20, 0.25}

   where r is desired maximal size. Rows: K2 values, 1 to 30. 
*/

double liml_size_vals[][8] = {
    { 16.38, 8.96, 6.66, 5.53,    0,    0,    0,    0 },
    {  8.68, 5.33, 4.42, 3.92, 7.03, 4.58, 3.95, 3.63 },
    {  6.46, 4.36, 3.69, 3.32, 5.44, 3.81, 3.32, 3.09 },
    {  5.44, 3.87, 3.30, 2.98, 4.72, 3.39, 2.99, 2.79 },
    {  4.84, 3.56, 3.05, 2.77, 4.32, 3.13, 2.78, 2.60 },
    {  4.45, 3.34, 2.87, 2.61, 4.06, 2.95, 2.63, 2.46 },
    {  4.18, 3.18, 2.73, 2.49, 3.90, 2.83, 2.52, 2.35 },
    {  3.97, 3.04, 2.63, 2.39, 3.78, 2.73, 2.43, 2.27 },
    {  3.81, 2.93, 2.54, 2.32, 3.70, 2.66, 2.36, 2.20 },
    {  3.68, 2.84, 2.46, 2.25, 3.64, 2.60, 2.30, 2.14 },
    {  3.58, 2.76, 2.40, 2.19, 3.60, 2.55, 2.25, 2.09 },
    {  3.50, 2.69, 2.34, 2.14, 3.58, 2.52, 2.21, 2.05 },
    {  3.42, 2.63, 2.29, 2.10, 3.56, 2.48, 2.17, 2.02 },
    {  3.36, 2.57, 2.25, 2.06, 3.55, 2.46, 2.14, 1.99 },
    {  3.31, 2.52, 2.21, 2.03, 3.54, 2.44, 2.11, 1.96 },
    {  3.27, 2.48, 2.18, 2.00, 3.55, 2.42, 2.09, 1.93 },
    {  3.24, 2.44, 2.14, 1.97, 3.55, 2.41, 2.07, 1.91 },
    {  3.20, 2.41, 2.11, 1.94, 3.56, 2.40, 2.05, 1.89 },
    {  3.18, 2.37, 2.09, 1.92, 3.57, 2.39, 2.03, 1.87 },
    {  3.21, 2.34, 2.06, 1.90, 3.58, 2.38, 2.02, 1.86 },
    {  3.39, 2.32, 2.04, 1.88, 3.59, 2.38, 2.01, 1.84 },
    {  3.57, 2.29, 2.02, 1.86, 3.60, 2.37, 1.99, 1.83 },
    {  3.68, 2.27, 2.00, 1.84, 3.62, 2.37, 1.98, 1.81 },
    {  3.75, 2.25, 1.98, 1.83, 3.64, 2.37, 1.98, 1.80 },
    {  3.79, 2.24, 1.96, 1.81, 3.65, 2.37, 1.97, 1.79 },
    {  3.82, 2.22, 1.95, 1.80, 3.67, 2.38, 1.96, 1.78 },
    {  3.85, 2.21, 1.93, 1.78, 3.74, 2.38, 1.96, 1.77 },
    {  3.86, 2.20, 1.92, 1.77, 3.87, 2.38, 1.95, 1.77 },
    {  3.87, 2.19, 1.90, 1.76, 4.02, 2.39, 1.95, 1.76 },
    {  3.88, 2.18, 1.89, 1.75, 4.12, 2.39, 1.95, 1.75 }
};

double sy_bvals[] = {0.05, 0.10, 0.20, 0.30}; /* maximal bias */
double sy_rvals[] = {0.10, 0.15, 0.20, 0.25}; /* maximal size */

/* 'which' codes: 

    1 = TSLS relative bias critical values
    2 = TSLS test size critical values
    3 = LIML test size critical values
*/

gretl_matrix *stock_yogo_lookup (int n, int K2, int which)
{
    gretl_matrix *v;
    const double *valrow;
    int nmax = (which == 1)? 3 : 2;
    int K2min = (which == 1)? 3 : 1;
    int i, c;

    if (n < 1 || n > nmax || K2 < K2min || K2 > 30 || K2 < n) {
	/* can't do it */
	return NULL;
    }

    v = gretl_matrix_alloc(2, 4);
    if (v == NULL) {
	return NULL;
    }

    if (which == 1) {
	valrow = tsls_bias_vals[K2 - 3];
	c = (n == 1)? 0 : (n == 2)? 4 : 8;
    } else if (which == 2) {
	valrow = tsls_size_vals[K2 - 1];
	c = (n == 1)? 0 : 4;
    } else {
	valrow = liml_size_vals[K2 - 1];
	c = (n == 1)? 0 : 4;
    }	

    for (i=0; i<4; i++) {
	if (which == 1) {
	    gretl_matrix_set(v, 0, i, sy_bvals[i]);
	} else {
	    gretl_matrix_set(v, 0, i, sy_rvals[i]);
	}
	gretl_matrix_set(v, 1, i, valrow[c+i]);
    }

    return v;
}
