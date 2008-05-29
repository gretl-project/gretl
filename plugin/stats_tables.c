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
