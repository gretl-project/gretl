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

/*
  odd_sort: creates the vector y = [x_0, x_2, ldots , x_{N-2}, x_1,
  x_3, ldots, x_{N-1}]' by rearranging the N = 2T elements of x =
  [x_0, x_1, ldots, x_{T-1}, x_T, x_{T+1}, ldots, x_{N-1}]'.
*/

static void odd_sort (int ncap, double *y)
{
    int tcap, t, j, k; 
    double x;
    
    tcap = ncap / 2;
    
    for (j = 1; j <= tcap - 1; j++) { 
	k = j;
	if (j > 1)
	    do {
		k = (2 * k) % (2 * tcap - 1);
            } while (k > j);
	if (k == j) { 
	    x = y[j];
	    t = j;
	    do {
		k = t;
		t = (2 * k) % (2 * tcap - 1);
		y[k] = y[t];
            } while (t != j);
	    y[k] = x;
	}
    }
} 

/*
  bit_reverse: maps the number j = j_1 + j_2 2 + cdots + j_g 2^{g-1}
  into t = j_1 2^{g-1} + j_2 2^{g-2} + cdots + j_g. This is used in
  the process of unscrambling the transformed vector in the final
  stage of the calculation of a Base-2 FFT.
*/

static int bit_reverse (int j, int g)
{
    int i, t = 0, r = j;
    
    for (i=0; i<g; i++) {
	t = t * 2 + r % 2;
	r = r / 2;
    }

    return t;
} 

/*
  base_2_fft: calculates the discrete Fourier transform of a
  complex-valued vector of T = 2^g elements. The real and the
  imaginary components of the vector are stored in a single array in
  the order [y^re_0, y^re_1, ldots , y^re_{T-1}, y^im_0, y^im_1,
  ldots , y^im_{T-1}] such that y^re_t and y^im_t are coded as
  y[t] and y[tcap + t] respectively.

  In order to conform to the usual definition of the Fourier
  transform, the elements of the vector y must be multiplied by
  the factor T^{-1}, either at the start of the procedure or on
  its completion.
*/

static void base_2_fft (double *y, int tcap, int g)
{
    int a, t, i, j, k, l;
    int c = 0, P = tcap, Q = 0;
    double theta, sine, cosine; 
    double yR, yI, W = M_2PI / tcap; 
    
    for (l = 1; l <= g; l++) { 
	a = 0;
	t = 0;
	if (l == 1) {
	    Q = 1;
	} else {
	    Q *= 2;
	}
	P = P / 2;
       
	for (i = 0; i <= Q - 1; i++) { 
	    for (k = 0; k <= P - 1; k++) { 
		t = a + c;
		yR = y[t];
		yI = y[t + tcap];
		y[t] +=  y[t + P];
		y[t + tcap] +=  y[t + P + tcap];
		y[t + P] = yR - y[t + P];
		y[t + P + tcap] = yI - y[t + P + tcap];
		
		if (l < g) { 
		    theta = W * ((c * Q) % tcap);
		    cosine = cos(theta);
		    sine = sin(theta);
		    yR = y[t + P];
		    yI = y[t + P + tcap];
		    y[t + P] = yR * cosine + yI * sine;
		    y[t + P + tcap] = yI * cosine - yR * sine;
		} 
		c++;
            } 
	    c = 0;
	    a += 2 * P;
	} 
    }
      
    for (j = 1; j <= tcap - 2; j++) { 
	t = bit_reverse(j, g);
	if (t > j) { 
	    yR = y[t];
	    yI = y[t + tcap];
	    y[t] = y[j];
	    y[t + tcap] = y[j + tcap];
	    y[j] = yR;
	    y[j + tcap] = yI;
	} 
    } 
} 

/*
  compact_real_fft: finds the discrete Fourier transform of a
  real-valued matrix y = [y_0, y_1, ldots, y_{T-1}, y_T, y_{T+1},
  ldots , y_{N-1}]' of N = 2T elements, where T = 2^(g+1). The
  real and imaginary parts are stored in two adjacent columns of the
  matrix ft.

  In order to conform to the usual definition of the Fourier
  transform, the elements of the vector y must be multiplied by
  the factor N^{-1}, either at the start of the procedure or on
  its completion.
*/

int compact_real_fft (const gretl_matrix *y, int g, gretl_matrix *ft)
{
    int ncap = gretl_matrix_rows(y); 
    int c = gretl_matrix_cols(y);
    int t, j, tcap; 
    double tmp, xReven, xRodd, xIeven, xIodd; 
    double theta, incr, sine, cosine; 
    double *x;
    int ffcr = 0;
    int ffci = 1;

    x = malloc(ncap * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    tcap = ncap / 2;
    incr = M_PI / tcap;
    
    for (j=0; j<c; j++) {

	for (t=0; t<ncap; t++) {
	    x[t] = gretl_matrix_get(y, t, j);
	}

	odd_sort(ncap, x);
	base_2_fft(x, tcap, g-1);
	theta = 0;

	for (t = 1; t <= tcap / 2; t++) {
	    theta += incr;
	    cosine = cos(theta);
	    sine = sin(theta);

	    xReven = (x[t] + x[tcap - t]) / 2;
	    xRodd = (x[t] - x[tcap - t]) / 2;
	    xIeven = (x[t + tcap] + x[ncap - t]) / 2;
	    xIodd = (x[t + tcap] - x[ncap - t]) / 2;

	    tmp = xReven + cosine * xIeven - sine * xRodd;
	    gretl_matrix_set(ft, t, ffcr, tmp);
	    gretl_matrix_set(ft, ncap-t, ffcr, tmp);

	    tmp = xReven - cosine * xIeven + sine * xRodd;
	    gretl_matrix_set(ft, tcap-t, ffcr, tmp);
	    gretl_matrix_set(ft, tcap+t, ffcr, tmp);

	    tmp = xIodd - cosine * xRodd - sine * xIeven;
	    gretl_matrix_set(ft, t, ffci, tmp);
	    gretl_matrix_set(ft, ncap-t, ffci, -tmp);
	
	    tmp  = -xIodd - cosine * xRodd - sine * xIeven;
	    gretl_matrix_set(ft, tcap-t, ffci, tmp);
	    gretl_matrix_set(ft, tcap+t, ffci, -tmp);
	}
      
	gretl_matrix_set(ft, 0, ffcr, x[0] + x[tcap]);
	gretl_matrix_set(ft, 0, ffci, 0);
	gretl_matrix_set(ft, tcap, ffcr, x[0] - x[tcap]);

	ffcr += 2;
	ffci += 2;
    }

    free(x);

    return 0;
} 

/**
 * gretl_matrix_fft:
 * @y: n x 2 input matrix (?)
 * @err: location to receive error code.
 *
 * Add description here.
 *
 * Returns: the generated matrix, or %NULL on failure.
 */

gretl_matrix *gretl_matrix_fft (const gretl_matrix *y, int *err)
{
    gretl_matrix *tmp = NULL;
    gretl_matrix *ft = NULL;
    int r = gretl_matrix_rows(y);
    int c = gretl_matrix_cols(y);
    int i, j;
    int N = 1;
    int g = 0;

    /* determine suitable power of 2 for padding */
    while (N < r) {
	g++;
	N <<= 1;
    }

    tmp = gretl_zero_matrix_new(N, c);
    ft = gretl_matrix_alloc(N, 2*c);

    if (tmp == NULL || ft == NULL) {
	gretl_matrix_free(tmp);
	gretl_matrix_free(ft);
	*err = E_ALLOC;
	return NULL;
    }
    
    /* now fill in the top of tmp */
    for (i=0; i<r; i++) {
	for (j=0; j<c; j++) {
	    gretl_matrix_set(tmp, i, j, gretl_matrix_get(y, i, j));
	}
    }

    *err = compact_real_fft(tmp, g, ft);
    if (*err) {
	gretl_matrix_free(ft);
	ft = NULL;
    }

    gretl_matrix_free(tmp);

    return ft;
}
