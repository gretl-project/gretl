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
   Stock and Yogo, 2003, Table 2.
   Critical Values for the Weak Instrument Test Based on TSLS Size,
   for nominal significance level 5%

   cols 1-4: n = 1; r = 0.10, 0.15, 0.20, 0.25
   cols 5-8: n = 2; r = 0.10, 0.15, 0.20, 0.25

   where r is desired maximal size. left margin: K2 
*/

double syogo_critvals[][8] = {
    /*  1 */ { 16.38,  8.96,  6.66,  5.53,     0,     0,     0,     0 },
    /*  2 */ { 19.93, 11.59,  8.75,  7.25,  7.03,  4.58,  3.95,  3.63 },
    /*  3 */ { 22.30, 12.83,  9.54,  7.80, 13.43,  8.18,  6.40,  5.45 },
    /*  4 */ { 24.58, 13.96, 10.26,  8.31, 16.87,  9.93,  7.54,  6.28 },
    /*  5 */ { 26.87, 15.09, 10.98,  8.84, 19.45, 11.22,  8.38,  6.89 },
    /*  6 */ { 29.18, 16.23, 11.72,  9.38, 21.68, 12.33,  9.10,  7.42 },
    /*  7 */ { 31.50, 17.38, 12.48,  9.93, 23.72, 13.34,  9.77,  7.91 },
    /*  8 */ { 33.84, 18.54, 13.24, 10.50, 25.64, 14.31, 10.41,  8.39 },
    /*  9 */ { 36.19, 19.71, 14.01, 11.07, 27.51, 15.24, 11.03,  8.85 },
    /* 10 */ { 38.54, 20.88, 14.78, 11.65, 29.32, 16.16, 11.65,  9.31 },
    /* 11 */ { 40.90, 22.06, 15.56, 12.23, 31.11, 17.06, 12.25,  9.77 },
    /* 12 */ { 43.27, 23.24, 16.35, 12.82, 32.88, 17.95, 12.86, 10.22 },
    /* 13 */ { 45.64, 24.42, 17.14, 13.41, 34.62, 18.84, 13.45, 10.68 },
    /* 14 */ { 48.01, 25.61, 17.93, 14.00, 36.36, 19.72, 14.05, 11.13 },
    /* 15 */ { 50.39, 26.80, 18.72, 14.60, 38.08, 20.60, 14.65, 11.58 },
    /* 16 */ { 52.77, 27.99, 19.51, 15.19, 39.80, 21.48, 15.24, 12.03 },
    /* 17 */ { 55.15, 29.19, 20.31, 15.79, 41.51, 22.35, 15.83, 12.49 },
    /* 18 */ { 57.53, 30.38, 21.10, 16.39, 43.22, 23.22, 16.42, 12.94 },
    /* 19 */ { 59.92, 31.58, 21.90, 16.99, 44.92, 24.09, 17.02, 13.39 },
    /* 20 */ { 62.30, 32.77, 22.70, 17.60, 46.62, 24.96, 17.61, 13.84 },
    /* 21 */ { 64.69, 33.97, 23.50, 18.20, 48.31, 25.82, 18.20, 14.29 },
    /* 22 */ { 67.07, 35.17, 24.30, 18.80, 50.01, 26.69, 18.79, 14.74 },
    /* 23 */ { 69.46, 36.37, 25.10, 19.41, 51.70, 27.56, 19.38, 15.19 },
    /* 24 */ { 71.85, 37.57, 25.90, 20.01, 53.39, 28.42, 19.97, 15.64 },
    /* 25 */ { 74.24, 38.77, 26.71, 20.61, 55.07, 29.29, 20.56, 16.10 },
    /* 26 */ { 76.62, 39.97, 27.51, 21.22, 56.76, 30.15, 21.15, 16.55 },
    /* 27 */ { 79.01, 41.17, 28.31, 21.83, 58.45, 31.02, 21.74, 17.00 },
    /* 28 */ { 81.40, 42.37, 29.12, 22.43, 60.13, 31.88, 22.33, 17.45 },
    /* 29 */ { 83.79, 43.57, 29.92, 23.04, 61.82, 32.74, 22.92, 17.90 },
    /* 30 */ { 86.17, 44.78, 30.72, 23.65, 63.51, 33.61, 23.51, 18.35 }
};

gretl_matrix *stock_yogo_lookup (int n, int K2)
{
    gretl_matrix *v;
    const double *valrow;
    int i, c;

    if (n < 1 || n > 2 || K2 < 1 || K2 > 30 || (n == 2 && K2 == 1)) {
	/* can't do it */
	return NULL;
    }

    v = gretl_vector_alloc(4);
    if (v == NULL) {
	return NULL;
    }

    valrow = syogo_critvals[K2 - 1];
    c = (n == 1)? 0 : 4;

    for (i=0; i<4; i++) {
	v->val[i] = valrow[c+i];
    }

    return v;
}
