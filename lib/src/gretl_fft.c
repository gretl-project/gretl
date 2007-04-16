/*
 *  Copyright (c) by Allin Cottrell and Riccardo "Jack" Lucchetti
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
#include "gretl_fft.h"
#include <fftw3.h>

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
    int i, j;

    int r = gretl_matrix_rows(y);
    int c = gretl_matrix_cols(y);
    int is_odd = (r % 2);
    int m = r/2;
    int cr = 0;
    int ci = 1;

    tmp = malloc(r * sizeof *tmp);
    ft = gretl_matrix_alloc(r, 2 * c);
    out = fftw_malloc(r * sizeof *out);

    if (tmp == NULL || ft == NULL || out == NULL) {
	free(tmp);
	gretl_matrix_free(ft);
	fftw_free(out);
	*err = E_ALLOC;
	return NULL;
    }
    
    /* now do the FFT proper */

    for (j=0; j<c; j++) {

	for (i=0; i<r; i++) {
	    tmp[i] = gretl_matrix_get(y, i, j);
	}
		
	if (j==0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_r2c_1d(r, tmp, out, FFTW_ESTIMATE);
	}

	fftw_execute(p);

	for (i=0; i<m+1+is_odd; i++) {
	    gretl_matrix_set(ft, i, cr, out[i][0]);
	    gretl_matrix_set(ft, i, ci, out[i][1]);
	}

	for (i=m; i>0; i--) {
	    gretl_matrix_set(ft, r-i, cr,  out[i][0]);
	    gretl_matrix_set(ft, r-i, ci, -out[i][1]);
	}

	cr +=2;
	ci +=2;
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    free(tmp);

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
    int i, j;

    int r = gretl_matrix_rows(y);
    int c = gretl_matrix_cols(y) / 2;
    int is_odd = (r % 2);
    int m = r/2;
    int cr = 0;
    int ci = 1;

    tmp = malloc(r * sizeof *tmp);
    ft = gretl_matrix_alloc(r, c);
    in = fftw_malloc(r * sizeof *in);

    if (tmp == NULL || ft == NULL || in == NULL) {
	free(tmp);
	gretl_matrix_free(ft);
	fftw_free(in);
	*err = E_ALLOC;
	return NULL;
    }
    
    /* now do the inverse FFT proper */

    for (j=0; j<c; j++) {

	for (i=0; i<m+1+is_odd; i++) {
	    in[i][0] = gretl_matrix_get(y, i, cr);
	    in[i][1] = gretl_matrix_get(y, i, ci);
	}

	if (j == 0) {
	    /* make the plan just once */
	    p = fftw_plan_dft_c2r_1d(r, in, tmp, FFTW_ESTIMATE);
	}

	fftw_execute(p);

	for (i=0; i<r; i++) {
	    gretl_matrix_set(ft, i, j, tmp[i]/r);
	}

	cr += 2;
	ci += 2;
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    free(tmp);

    return ft;
}
