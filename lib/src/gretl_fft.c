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
#include "gretl_fft.h"
#include <fftw3.h>

static int fft_allocate (double **px, gretl_matrix **pm,
			 fftw_complex **pc, int r, int c)
{
    *pm = gretl_matrix_alloc(r, c);
    if (*pm == NULL) {
	return E_ALLOC;
    }

    *px = fftw_malloc(r * sizeof **px);
    if (*px == NULL) {
	gretl_matrix_free(*pm);
	return E_ALLOC;
    }

    *pc = fftw_malloc((r/2 + 1 + r % 2) * sizeof **pc);
    if (*pc == NULL) {
	gretl_matrix_free(*pm);
	free(*px);
	return E_ALLOC;
    }

    return 0;
}

/**
 * gretl_matrix_fft:
 * @y: input matrix.
 * @err: location to receive error code.
 *
 * Add description here.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err)
{
    gretl_matrix *ft = NULL;
    fftw_plan p = NULL;
    double *tmp = NULL;
    fftw_complex *out;
    int r = gretl_matrix_rows(y);
    int c, m, odd, cr, ci;
    int i, j;

    if (r < 2) {
	*err = E_DATA;
	return NULL;
    }

    c = gretl_matrix_cols(y);
    m = r / 2;
    odd = r % 2;
    cr = 0;
    ci = 1;

    *err = fft_allocate(&tmp, &ft, &out, r, 2 * c);
    if (*err) {
	return NULL;
    }

    for (j=0; j<c; j++) {

	for (i=0; i<r; i++) {
	    tmp[i] = gretl_matrix_get(y, i, j);
	}
		
	if (j == 0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_r2c_1d(r, tmp, out, FFTW_ESTIMATE);
	}

	fftw_execute(p);

	for (i=0; i<=m+odd; i++) {
	    gretl_matrix_set(ft, i, cr, out[i][0]);
	    gretl_matrix_set(ft, i, ci, out[i][1]);
	}

	for (i=m; i>0; i--) {
	    gretl_matrix_set(ft, r-i, cr,  out[i][0]);
	    gretl_matrix_set(ft, r-i, ci, -out[i][1]);
	}

	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    fftw_free(tmp);

    return ft;
}

/**
 * gretl_matrix_ffti:
 * @y: input matrix.
 * @err: location to receive error code.
 *
 * Add description here.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_ffti (const gretl_matrix *y, int *err)
{
    gretl_matrix *ft = NULL;
    fftw_plan p = NULL;
    double *tmp = NULL;
    fftw_complex *in;
    int c, r = gretl_matrix_rows(y);
    int m, odd, cr, ci;
    int i, j;

    if (r < 2) {
	*err = E_DATA;
	return NULL;
    }

    c = gretl_matrix_cols(y) / 2;
    m = r / 2;
    odd = r % 2;

    if (c == 0) {
	*err = E_NONCONF;
	return NULL;
    }

    *err = fft_allocate(&tmp, &ft, &in, r, c);
    if (*err) {
	return NULL;
    }

    cr = 0;
    ci = 1;
    
    for (j=0; j<c; j++) {

	for (i=0; i<=m+odd; i++) {
	    in[i][0] = gretl_matrix_get(y, i, cr);
	    in[i][1] = gretl_matrix_get(y, i, ci);
	}

	if (j == 0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_c2r_1d(r, in, tmp, FFTW_ESTIMATE);
	}

	fftw_execute(p);

	for (i=0; i<r; i++) {
	    gretl_matrix_set(ft, i, j, tmp[i] / r);
	}

	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(tmp);

    return ft;
}
