/* 
   lbfgsb.c -- translated from lbfgsb.f (plus a couple of functions
   from linpack.f) via f2c, then subject to minimal clean-up by
   Allin Cottrell, 2015-07-26. The original source is version 3.0
   of L-BFGS-B, downloaded from

   http://users.iems.northwestern.edu/~nocedal/lbfgsb.html

   and subject to the "New BSD License" or "Modified BSD License,"
   as follows:

   Copyright (c) <year>, <copyright holder>
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "gretl_f2c.h"
#include "clapack_double.h"
#include "minpack.h"

static double c_b9 = 0.;
static int c__1 = 1;
static int c__11 = 11;
static double c_b280 = .001;
static double c_b281 = .9;
static double c_b282 = .1;

static int dpofa_(double *a, int *lda, int *n, int *info)
{
    int a_dim1, a_offset, i__1, i__2, i__3;

    int j, k;
    double s, t;
    int jm1;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

   i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k - 1;
	    t = a[k + j * a_dim1] - ddot_(&i__3, &a[k * a_dim1 + 1], &c__1, &
					  a[j * a_dim1 + 1], &c__1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
	}
    L20:
	s = a[j + j * a_dim1] - s;
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
    }
    *info = 0;
 L40:
    return 0;
} /* dpofa_ */

static int dtrsl_(double *t, int *ldt, int *n, 
		  double *b, int *job, int *info)
{
    int t_dim1, t_offset, i__1, i__2;

    int j, jj, case__;
    double temp;

    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
	if (t[*info + *info * t_dim1] == 0.) {
	    goto L150;
	}
    }
    *info = 0;

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

    /* solve t*x=b for t lower triangular */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = -b[j - 1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	b[j] /= t[j + j * t_dim1];
    }
L40:
    goto L140;

    /* solve t*x=b for t upper triangular. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	daxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
    }
L70:
    goto L140;

    /* solve trans(t)*x=b for t lower triangular. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = jj - 1;
	b[j] -= ddot_(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	b[j] /= t[j + j * t_dim1];
    }
L100:
    goto L140;

    /* solve trans(t)*x=b for t upper triangular. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	b[j] -= ddot_(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
    }
L130:
L140:
L150:
    return 0;
} /* dtrsl_ */

static int active_(int *n, double *l, double *u, int *nbd,
		   double *x, int *iwhere, logical *prjctd,
		   logical *cnstnd, logical *boxed)
{
    int i__1;
    int i__;

    --iwhere;
    --x;
    --nbd;
    --u;
    --l;

    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] > 0) {
	    if (nbd[i__] <= 2 && x[i__] <= l[i__]) {
		if (x[i__] < l[i__]) {
		    *prjctd = 1;
		    x[i__] = l[i__];
		}
	    } else if (nbd[i__] >= 2 && x[i__] >= u[i__]) {
		if (x[i__] > u[i__]) {
		    *prjctd = 1;
		    x[i__] = u[i__];
		}
	    }
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] != 2) {
	    *boxed = 0;
	}
	if (nbd[i__] == 0) {
	    iwhere[i__] = -1;
	} else {
	    *cnstnd = 1;
	    if (nbd[i__] == 2 && u[i__] - l[i__] <= 0.) {
		iwhere[i__] = 3;
	    } else {
		iwhere[i__] = 0;
	    }
	}
    }

    return 0;
} /* active_ */

static int bmv_(int *m, double *sy, double *wt, int *col,
		double *v, double *p, int *info)
{
    int sy_dim1, sy_offset, wt_dim1, wt_offset, i__1, i__2;

    int i__, k, i2;
    double sum;

    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    --p;
    --v;

    if (*col == 0) {
	return 0;
    }
    p[*col + 1] = v[*col + 1];
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i2 = *col + i__;
	sum = 0.;
	i__2 = i__ - 1;
	for (k = 1; k <= i__2; ++k) {
	    sum += sy[i__ + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
	}
	p[i2] = v[i2] + sum;
   }
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = v[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
   }
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = -p[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *col;
	for (k = i__ + 1; k <= i__2; ++k) {
	    sum += sy[k + i__ * sy_dim1] * p[*col + k] / sy[i__ + i__ * 
		    sy_dim1];
	}
	p[i__] += sum;
    }
    return 0;
} /* bmv_ */

static int cmprlb_(int *n, int *m, double *x, 
		   double *g, double *ws, double *wy, double *sy, 
		   double *wt, double *z__, double *r__, double *wa, 
		   int *index, double *theta, int *col, int *head, 
		   int *nfree, logical *cnstnd, int *info)
{
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	wt_dim1, wt_offset, i__1, i__2;

    int i__, j, k;
    double a1, a2;
    int pointr;

    --index;
    --r__;
    --z__;
    --g;
    --x;
    --wa;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    if (! (*cnstnd) && *col > 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__] = -g[i__];
	}
    } else {
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    r__[i__] = -(*theta) * (z__[k] - x[k]) - g[k];
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(*m << 1) + 1], &wa[
		1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
	}
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    a1 = wa[j];
	    a2 = *theta * wa[*col + j];
	    i__2 = *nfree;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		k = index[i__];
		r__[i__] = r__[i__] + wy[k + pointr * wy_dim1] * a1 + ws[k + 
			pointr * ws_dim1] * a2;
	    }
	    pointr = pointr % *m + 1;
	}
    }
    return 0;
} /* cmprlb_ */

static int errclb_(int *n, int *m, double *factr, 
		   double *l, double *u, int *nbd, char *task,
		   int *info, int *k)
{
    int i__1;

    int i__;

    --nbd;
    --u;
    --l;

    if (*n <= 0) {
	strcpy(task, "ERROR: N .LE. 0");
    }
    if (*m <= 0) {
	strcpy(task, "ERROR: M .LE. 0");
    }
    if (*factr < 0.) {
	strcpy(task, "ERROR: FACTR .LT. 0");
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] < 0 || nbd[i__] > 3) {
	    strcpy(task, "ERROR: INVALID NBD");
	    *info = -6;
	    *k = i__;
	}
	if (nbd[i__] == 2) {
	    if (l[i__] > u[i__]) {
		strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
		*info = -7;
		*k = i__;
	    }
	}
    }
    return 0;
} /* errclb_ */

static int formk_(int *n, int *nsub, int *ind, int *nenter,
		  int *ileave, int *indx2, int *iupdat, logical *
		  updatd, double *wn, double *wn1, int *m, double *ws, 
		  double *wy, double *sy, double *theta, int *col, 
		  int *head, int *info)
{
    int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset, 
	    wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2, i__3;

    int i__, k, k1, m2, is, js, iy, jy, is1, js1, col2, dend, pend;
    int upcl;
    double temp1, temp2, temp3, temp4;
    int ipntr, jpntr, dbegin, pbegin;

    --indx2;
    --ind;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wn1_dim1 = 2 * *m;
    wn1_offset = 1 + wn1_dim1;
    wn1 -= wn1_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;

    if (*updatd) {
	if (*iupdat > *m) {
	    i__1 = *m - 1;
	    for (jy = 1; jy <= i__1; ++jy) {
		js = *m + jy;
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[jy + 1 + (jy + 1) * wn1_dim1], &c__1,
		       &wn1[jy + jy * wn1_dim1], &c__1);
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[js + 1 + (js + 1) * wn1_dim1], &c__1,
		       &wn1[js + js * wn1_dim1], &c__1);
		i__2 = *m - 1;
		dcopy_(&i__2, &wn1[*m + 2 + (jy + 1) * wn1_dim1], &c__1,
		       &wn1[*m + 1 + jy * wn1_dim1], &c__1);
	    }
	}
	pbegin = 1;
	pend = *nsub;
	dbegin = *nsub + 1;
	dend = *n;
	iy = *col;
	is = *m + *col;
	ipntr = *head + *col - 1;
	if (ipntr > *m) {
	    ipntr -= *m;
	}
	jpntr = *head;
	i__1 = *col;
	for (jy = 1; jy <= i__1; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    i__2 = dend;
	    for (k = dbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    wn1[iy + jy * wn1_dim1] = temp1;
	    wn1[is + js * wn1_dim1] = temp2;
	    wn1[is + jy * wn1_dim1] = temp3;
	    jpntr = jpntr % *m + 1;
	}
	jy = *col;
	jpntr = *head + *col - 1;
	if (jpntr > *m) {
	    jpntr -= *m;
	}
	ipntr = *head;
	i__1 = *col;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = *m + i__;
	    temp3 = 0.;
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    ipntr = ipntr % *m + 1;
	    wn1[is + jy * wn1_dim1] = temp3;
	}
	upcl = *col - 1;
    } else {
	upcl = *col;
    }
    ipntr = *head;
    i__1 = upcl;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *m + iy;
	jpntr = *head;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
	    }
	    wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
	    wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
	    jpntr = jpntr % *m + 1;
	}
	ipntr = ipntr % *m + 1;
    }
    ipntr = *head;
    i__1 = *m + upcl;
    for (is = *m + 1; is <= i__1; ++is) {
	jpntr = *head;
	i__2 = upcl;
	for (jy = 1; jy <= i__2; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
	    }
	    if (is <= jy + *m) {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] + temp1 - 
			temp3;
	    } else {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] - temp1 + 
			temp3;
	    }
	    jpntr = jpntr % *m + 1;
	}
	ipntr = ipntr % *m + 1;
    }
    m2 = *m << 1;
    i__1 = *col;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *col + iy;
	is1 = *m + iy;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *col + jy;
	    js1 = *m + jy;
	    wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
	    wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
	}
	i__2 = iy - 1;
	for (jy = 1; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
	}
	i__2 = *col;
	for (jy = iy; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
	}
	wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
    }
    dpofa_(&wn[wn_offset], &m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }
    col2 = *col << 1;
    i__1 = col2;
    for (js = *col + 1; js <= i__1; ++js) {
	dtrsl_(&wn[wn_offset], &m2, col, &wn[js * wn_dim1 + 1], &c__11, info);
    }
    i__1 = col2;
    for (is = *col + 1; is <= i__1; ++is) {
	i__2 = col2;
	for (js = is; js <= i__2; ++js) {
	    wn[is + js * wn_dim1] += ddot_(col, &wn[is * wn_dim1 + 1], &c__1, 
		    &wn[js * wn_dim1 + 1], &c__1);
	}
    }
    dpofa_(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }
    return 0;
} /* formk_ */

static int formt_(int *m, double *wt, double *sy, 
		  double *ss, int *col, double *theta,
		  int *info)
{
    int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1, 
	i__2, i__3;

    int i__, j, k, k1;
    double ddum;

    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;

    i__1 = *col;
    for (j = 1; j <= i__1; ++j) {
	wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
    }
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *col;
	for (j = i__; j <= i__2; ++j) {
	    k1 = min(i__,j) - 1;
	    ddum = 0.;
	    i__3 = k1;
	    for (k = 1; k <= i__3; ++k) {
		ddum += sy[i__ + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k + 
			k * sy_dim1];
	    }
	    wt[i__ + j * wt_dim1] = ddum + *theta * ss[i__ + j * ss_dim1];
	}
    }
    dpofa_(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt_ */

static int freev_(int *n, int *nfree, int *index, 
		  int *nenter, int *ileave, int *indx2, int *iwhere, 
		  logical *wrk, logical *updatd, logical *cnstnd,
		  int *iter)
{
    int i__1;

    int i__, k, iact;

    --iwhere;
    --indx2;
    --index;

    *nenter = 0;
    *ileave = *n + 1;
    if (*iter > 0 && *cnstnd) {
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] > 0) {
		--(*ileave);
		indx2[*ileave] = k;
	    }
	}
	i__1 = *n;
	for (i__ = *nfree + 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] <= 0) {
		++(*nenter);
		indx2[*nenter] = k;
	    }
	}
    }
    *wrk = *ileave < *n + 1 || *nenter > 0 || *updatd;
    *nfree = 0;
    iact = *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iwhere[i__] <= 0) {
	    ++(*nfree);
	    index[*nfree] = i__;
	} else {
	    --iact;
	    index[iact] = i__;
	}
    }

    return 0;
} /* freev_ */

static int hpsolb_(int *n, double *t, int *iorder, 
		   int *iheap)
{
    int i__1;

    int i__, j, k;
    double out, ddum;
    int indxin, indxou;

    --iorder;
    --t;

    if (*iheap == 0) {
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
	    i__ = k;
L10:
	    if (i__ > 1) {
		j = i__ / 2;
		if (ddum < t[j]) {
		    t[i__] = t[j];
		    iorder[i__] = iorder[j];
		    i__ = j;
		    goto L10;
		}
	    }
	    t[i__] = ddum;
	    iorder[i__] = indxin;
	}
    }
    if (*n > 1) {
	i__ = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[*n];
	indxin = iorder[*n];
L30:
	j = i__ + i__;
	if (j <= *n - 1) {
	    if (t[j + 1] < t[j]) {
		++j;
	    }
	    if (t[j] < ddum) {
		t[i__] = t[j];
		iorder[i__] = iorder[j];
		i__ = j;
		goto L30;
	    }
	}
	t[i__] = ddum;
	iorder[i__] = indxin;
	t[*n] = out;
	iorder[*n] = indxou;
    }
    return 0;
} /* hpsolb_ */

static int matupd_(int *n, int *m, double *ws, 
		   double *wy, double *sy, double *ss, double *d__, 
		   double *r__, int *itail, int *iupdat, int *col, 
		   int *head, double *theta, double *rr, double *dr, 
		   double *stp, double *dtd)
{
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    ss_dim1, ss_offset, i__1, i__2;

    int j;
    int pointr;

    --r__;
    --d__;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;

    if (*iupdat <= *m) {
	*col = *iupdat;
	*itail = (*head + *iupdat - 2) % *m + 1;
    } else {
	*itail = *itail % *m + 1;
	*head = *head % *m + 1;
    }
    dcopy_(n, &d__[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
    dcopy_(n, &r__[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
    *theta = *rr / *dr;
    if (*iupdat > *m) {
	i__1 = *col - 1;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1, &ss[j * ss_dim1 + 1]
		    , &c__1);
	    i__2 = *col - j;
	    dcopy_(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1, &sy[j + j * 
		    sy_dim1], &c__1);
	}
    }
    pointr = *head;
    i__1 = *col - 1;
    for (j = 1; j <= i__1; ++j) {
	sy[*col + j * sy_dim1] = ddot_(n, &d__[1], &c__1, &wy[pointr * 
		wy_dim1 + 1], &c__1);
	ss[j + *col * ss_dim1] = ddot_(n, &ws[pointr * ws_dim1 + 1], &c__1, &
		d__[1], &c__1);
	pointr = pointr % *m + 1;
    }
    if (*stp == 1.) {
	ss[*col + *col * ss_dim1] = *dtd;
    } else {
	ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
    }
    sy[*col + *col * sy_dim1] = *dr;
    return 0;
} /* matupd_ */

static int projgr_(int *n, double *l, double *u, int *nbd,
		   double *x, double *g, double *sbgnrm)
{
    int i__1;
    double d__1, d__2;

    int i__;
    double gi;

    --g;
    --x;
    --nbd;
    --u;
    --l;

    *sbgnrm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = g[i__];
	if (nbd[i__] != 0) {
	    if (gi < 0.) {
		if (nbd[i__] >= 2) {
		    d__1 = x[i__] - u[i__];
		    gi = max(d__1,gi);
		}
	    } else {
		if (nbd[i__] <= 2) {
		    d__1 = x[i__] - l[i__];
		    gi = min(d__1,gi);
		}
	    }
	}
	d__1 = *sbgnrm, d__2 = fabs(gi);
	*sbgnrm = max(d__1,d__2);
    }
    return 0;
} /* projgr_ */

static int subsm_(int *n, int *m, int *nsub, int *ind,
		  double *l, double *u, int *nbd, double *x, 
		  double *d__, double *xp, double *ws, double *wy, 
		  double *theta, double *xx, double *gg, int *col, 
		  int *head, int *iword, double *wv, double *wn, 
		  int *info)
{
    int ws_dim1, ws_offset, wy_dim1, wy_offset, wn_dim1, wn_offset,
	i__1, i__2;
    double d__1, d__2;

    int i__, j, k, m2;
    double dk;
    int js, jy;
    double xk;
    int ibd, col2;
    double dd_p__, temp1, temp2, alpha;
    int pointr;

    --gg;
    --xx;
    --xp;
    --d__;
    --x;
    --nbd;
    --u;
    --l;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    --wv;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    --ind;

    if (*nsub <= 0) {
	return 0;
    }
    pointr = *head;
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 0.;
	temp2 = 0.;
	i__2 = *nsub;
	for (j = 1; j <= i__2; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * wy_dim1] * d__[j];
	    temp2 += ws[k + pointr * ws_dim1] * d__[j];
	}
	wv[i__] = temp1;
	wv[*col + i__] = *theta * temp2;
	pointr = pointr % *m + 1;
    }
    m2 = *m << 1;
    col2 = *col << 1;
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wv[i__] = -wv[i__];
    }
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
    pointr = *head;
    i__1 = *col;
    for (jy = 1; jy <= i__1; ++jy) {
	js = *col + jy;
	i__2 = *nsub;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = ind[i__];
	    d__[i__] = d__[i__] + wy[k + pointr * wy_dim1] * wv[jy] / *theta 
		    + ws[k + pointr * ws_dim1] * wv[js];
	}
	pointr = pointr % *m + 1;
    }
    d__1 = 1. / *theta;
    dscal_(nsub, &d__1, &d__[1], &c__1);
    *iword = 0;
    dcopy_(n, &x[1], &c__1, &xp[1], &c__1);
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	dk = d__[i__];
	xk = x[k];
	if (nbd[k] != 0) {
	    if (nbd[k] == 1) {
		/* lower bounds only */
		d__1 = l[k], d__2 = xk + dk;
		x[k] = max(d__1,d__2);
		if (x[k] == l[k]) {
		    *iword = 1;
		}
	    } else {
		if (nbd[k] == 2) {
		    /* upper and lower bounds */
		    d__1 = l[k], d__2 = xk + dk;
		    xk = max(d__1,d__2);
		    d__1 = u[k];
		    x[k] = min(d__1,xk);
		    if (x[k] == l[k] || x[k] == u[k]) {
			*iword = 1;
		    }
		} else {
		    if (nbd[k] == 3) {
			/* upper bounds only */
			d__1 = u[k], d__2 = xk + dk;
			x[k] = min(d__1,d__2);
			if (x[k] == u[k]) {
			    *iword = 1;
			}
		    }
		}
	    }
	} else {
	    /* free variables */
	    x[k] = xk + dk;
	}
    }
    if (*iword == 0) {
	goto L911;
    }
    dd_p__ = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dd_p__ += (x[i__] - xx[i__]) * gg[i__];
    }
    if (dd_p__ > 0.) {
	dcopy_(n, &xp[1], &c__1, &x[1], &c__1);
    } else {
	goto L911;
    }
    alpha = 1.;
    temp1 = alpha;
    ibd = 0;
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	dk = d__[i__];
	if (nbd[k] != 0) {
	    if (dk < 0. && nbd[k] <= 2) {
		temp2 = l[k] - x[k];
		if (temp2 >= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha < temp2) {
		    temp1 = temp2 / dk;
		}
	    } else if (dk > 0. && nbd[k] >= 2) {
		temp2 = u[k] - x[k];
		if (temp2 <= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha > temp2) {
		    temp1 = temp2 / dk;
		}
	    }
	    if (temp1 < alpha) {
		alpha = temp1;
		ibd = i__;
	    }
	}
    }
    if (alpha < 1.) {
	dk = d__[ibd];
	k = ind[ibd];
	if (dk > 0.) {
	    x[k] = u[k];
	    d__[ibd] = 0.;
	} else if (dk < 0.) {
	    x[k] = l[k];
	    d__[ibd] = 0.;
	}
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	x[k] += alpha * d__[i__];
    }
L911:

    return 0;
} /* subsm_ */

static int dcstep_(double *stx, double *fx, double *dx, 
		   double *sty, double *fy, double *dy, double *stp, 
		   double *fp, double *dp, logical *brackt,
		   double *stpmin, double *stpmax)
{
    double d__1, d__2, d__3;

    double p, q, r__, s, sgnd, stpc, stpf, stpq, gamma, theta;

    sgnd = *dp * (*dx / fabs(*dx));
    if (*fp > *fx) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d__1 = fabs(theta), d__2 = fabs(*dx),
	    d__1 = max(d__1,d__2), d__2 = fabs(*dp);
	s = max(d__1,d__2);
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + *dp;
	r__ = p / q;
	stpc = *stx + r__ * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp 
		- *stx);
	if ((d__1 = stpc - *stx, fabs(d__1)) < (d__2 = stpq - *stx, fabs(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	*brackt = 1;
    } else if (sgnd < 0.) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = max(d__1,d__2),
	    d__2 = fabs(*dp);
	s = max(d__1,d__2);
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma - *dp + gamma + *dx;
	r__ = p / q;
	stpc = *stp + r__ * (*stx - *stp);
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if ((d__1 = stpc - *stp, fabs(d__1)) > (d__2 = stpq - *stp, fabs(d__2))) {
	    stpf = stpc;
	} else {
	    stpf = stpq;
	}
	*brackt = 1;
    } else if (fabs(*dp) < fabs(*dx)) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
	d__1 = fabs(theta), d__2 = fabs(*dx), d__1 = max(d__1,d__2),
	    d__2 = fabs(*dp);
	s = max(d__1,d__2);
	d__3 = theta / s;
	d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
	gamma = s * sqrt((max(d__1,d__2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma + (*dx - *dp) + gamma;
	r__ = p / q;
	if (r__ < 0. && gamma != 0.) {
	    stpc = *stp + r__ * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = *stpmax;
	} else {
	    stpc = *stpmin;
	}
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if (*brackt) {
	    if ((d__1 = stpc - *stp, fabs(d__1)) < (d__2 = stpq - *stp, fabs(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    if (*stp > *stx) {
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = min(d__1,stpf);
	    } else {
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = max(d__1,stpf);
	    }
	} else {
	    if ((d__1 = stpc - *stp, fabs(d__1)) >
		(d__2 = stpq - *stp, fabs(d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min(*stpmax,stpf);
	    stpf = max(*stpmin,stpf);
	}
    } else {
	if (*brackt) {
	    theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
	    d__1 = fabs(theta), d__2 = fabs(*dy),
		d__1 = max(d__1,d__2), d__2 = fabs(*dp);
	    s = max(d__1,d__2);
	    d__1 = theta / s;
	    gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - *dp + theta;
	    q = gamma - *dp + gamma + *dy;
	    r__ = p / q;
	    stpc = *stp + r__ * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = *stpmax;
	} else {
	    stpf = *stpmin;
	}
    }
    if (*fp > *fx) {
	*sty = *stp;
	*fy = *fp;
	*dy = *dp;
    } else {
	if (sgnd < 0.) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = *fp;
	*dx = *dp;
    }
    *stp = stpf;
    return 0;
} /* dcstep_ */

static int dcsrch_(double *f, double *g, double *stp, 
		   double *ftol, double *gtol, double *xtol,
		   double *stpmin, double *stpmax, char *task,
		   int *isave, double *dsave)
{
    double d__1;

    double fm, gm, fx, fy, gx, gy, fxm, fym, gxm, gym, stx, sty;
    int stage;
    double finit, ginit, width, ftest, gtest, stmin, stmax, width1;
    logical brackt;

    --dsave;
    --isave;

    if (strncmp(task, "START", 5) == 0) {
	if (*stp < *stpmin) {
	    strcpy(task, "ERROR: STP .LT. STPMIN");
	}
	if (*stp > *stpmax) {
	    strcpy(task, "ERROR: STP .GT. STPMAX");
	}
	if (*g >= 0.) {
	    strcpy(task, "ERROR: INITIAL G .GE. ZERO");
	}
	if (*ftol < 0.) {
	    strcpy(task, "ERROR: FTOL .LT. ZERO");
	}
	if (*gtol < 0.) {
	    strcpy(task, "ERROR: GTOL .LT. ZERO");
	}
	if (*xtol < 0.) {
	    strcpy(task, "ERROR: XTOL .LT. ZERO");
	}
	if (*stpmin < 0.) {
	    strcpy(task, "ERROR: STPMIN .LT. ZERO");
	}
	if (*stpmax < *stpmin) {
	    strcpy(task, "ERROR: STPMAX .LT. STPMIN");
	}
	if (strncmp(task, "ERROR", 5) == 0) {
	    return 0;
	}
	brackt = 0;
	stage = 1;
	finit = *f;
	ginit = *g;
	gtest = *ftol * ginit;
	width = *stpmax - *stpmin;
	width1 = width / .5;
	stx = 0.;
	fx = finit;
	gx = ginit;
	sty = 0.;
	fy = finit;
	gy = ginit;
	stmin = 0.;
	stmax = *stp + *stp * 4.;
	strcpy(task, "FG");
	goto L1000;
    } else {
	if (isave[1] == 1) {
	    brackt = 1;
	} else {
	    brackt = 0;
	}
	stage = isave[2];
	ginit = dsave[1];
	gtest = dsave[2];
	gx = dsave[3];
	gy = dsave[4];
	finit = dsave[5];
	fx = dsave[6];
	fy = dsave[7];
	stx = dsave[8];
	sty = dsave[9];
	stmin = dsave[10];
	stmax = dsave[11];
	width = dsave[12];
	width1 = dsave[13];
    }
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
	stage = 2;
    }
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
	strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	strcpy(task, "WARNING: XTOL TEST SATISFIED");
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	strcpy(task, "WARNING: STP = STPMAX");
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	strcpy(task, "WARNING: STP = STPMIN");
    }
    if (*f <= ftest && fabs(*g) <= *gtol * (-ginit)) {
	strcpy(task, "CONVERGENCE");
    }
    if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0) {
	goto L1000;
    }
    if (stage == 1 && *f <= fx && *f > ftest) {
	fm = *f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = *g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
		stmin, &stmax);
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
		stmax);
    }
    if (brackt) {
	if ((d__1 = sty - stx, fabs(d__1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d__1 = sty - stx, fabs(d__1));
    }
    if (brackt) {
	stmin = min(stx,sty);
	stmax = max(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }
    *stp = max(*stp,*stpmin);
    *stp = min(*stp,*stpmax);
    if ((brackt && (*stp <= stmin || *stp >= stmax)) ||
	(brackt && stmax - stmin <= *xtol * stmax)) {
	*stp = stx;
    }
    strcpy(task, "FG");
L1000:
    if (brackt) {
	isave[1] = 1;
    } else {
	isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
} /* dcsrch_ */

static int lnsrlb_(int *n, double *l, double *u, 
		   int *nbd, double *x, double *f, double *fold, 
		   double *gd, double *gdold, double *g, double *d__, 
		   double *r__, double *t, double *z__, double *stp, 
		   double *dnorm, double *dtd, double *xstep, double *stpmx,
		   int *iter, int *ifun, int *iback, int *nfgv, 
		   int *info, char *task, logical *boxed, logical *cnstnd,
		   char *csave, int *isave, double *dsave)
{
    int i__1;
    double d__1;

    int i__;
    double a1, a2;

    --z__;
    --t;
    --r__;
    --d__;
    --g;
    --x;
    --nbd;
    --u;
    --l;
    --isave;
    --dsave;

    if (strncmp(task, "FG_LN", 5) == 0) {
	goto L556;
    }
    *dtd = ddot_(n, &d__[1], &c__1, &d__[1], &c__1);
    *dnorm = sqrt(*dtd);
    *stpmx = 1e10;
    if (*cnstnd) {
	if (*iter == 0) {
	    *stpmx = 1.;
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a1 = d__[i__];
		if (nbd[i__] != 0) {
		    if (a1 < 0. && nbd[i__] <= 2) {
			a2 = l[i__] - x[i__];
			if (a2 >= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx < a2) {
			    *stpmx = a2 / a1;
			}
		    } else if (a1 > 0. && nbd[i__] >= 2) {
			a2 = u[i__] - x[i__];
			if (a2 <= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx > a2) {
			    *stpmx = a2 / a1;
			}
		    }
		}
	    }
	}
    }
    if (*iter == 0 && ! (*boxed)) {
	d__1 = 1. / *dnorm;
	*stp = min(d__1,*stpmx);
    } else {
	*stp = 1.;
    }
    dcopy_(n, &x[1], &c__1, &t[1], &c__1);
    dcopy_(n, &g[1], &c__1, &r__[1], &c__1);
    *fold = *f;
    *ifun = 0;
    *iback = 0;
    strcpy(csave, "START");
L556:
    *gd = ddot_(n, &g[1], &c__1, &d__[1], &c__1);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
	    *info = -4;
	    return 0;
	}
    }
    dcsrch_(f, gd, stp, &c_b280, &c_b281, &c_b282, &c_b9, stpmx, csave, &
	    isave[1], &dsave[1]);
    *xstep = *stp * *dnorm;
    if (strncmp(csave, "CONV", 4) != 0 &&
	strncmp(csave, "WARN", 4) != 0) {
	strcpy(task, "FG_LNSRCH");
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, &z__[1], &c__1, &x[1], &c__1);
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = *stp * d__[i__] + t[i__];
	    }
	}
    } else {
	strcpy(task, "NEW_X");
    }
    return 0;
} /* lnsrlb_ */

static void timer_(double *ttime)
{
    clock_t tim = clock();

    *ttime = (double) tim / CLOCKS_PER_SEC;
}

static int cauchy_(int *n, double *x, double *l, 
		   double *u, int *nbd, double *g, int *iorder,
		   int *iwhere, double *t, double *d__, double *xcp, int *m, 
		   double *wy, double *ws, double *sy, double *wt, 
		   double *theta, int *col, int *head, double *p, 
		   double *c__, double *wbp, double *v, int *nseg, 
		   double *sbgnrm, int *info, double *epsmch)
{
    int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset, 
	wt_dim1, wt_offset, i__1, i__2;
    double d__1;

    int i__, j;
    double f1, f2, dt, tj, tj0;
    double tu = 0.0, tl = 0.0;
    int ibp;
    double dtm;
    double wmc, wmp, wmw;
    int col2;
    double dibp;
    int iter;
    double zibp, tsum, dibp2;
    logical bnded;
    double neggi;
    int nfree;
    double bkmin;
    int nleft;
    double f2_org__;
    int nbreak, ibkmin;
    int pointr;
    logical xlower, xupper;

    --xcp;
    --d__;
    --t;
    --iwhere;
    --iorder;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --v;
    --wbp;
    --c__;
    --p;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;

    if (*sbgnrm <= 0.) {
	dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
	return 0;
    }
    bnded = 1;
    nfree = *n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = *col << 1;
    f1 = 0.;
    i__1 = col2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neggi = -g[i__];
	if (iwhere[i__] != 3 && iwhere[i__] != -1) {
	    if (nbd[i__] <= 2) {
		tl = x[i__] - l[i__];
	    }
	    if (nbd[i__] >= 2) {
		tu = u[i__] - x[i__];
	    }
	    xlower = nbd[i__] <= 2 && tl <= 0.;
	    xupper = nbd[i__] >= 2 && tu <= 0.;
	    iwhere[i__] = 0;
	    if (xlower) {
		if (neggi <= 0.) {
		    iwhere[i__] = 1;
		}
	    } else if (xupper) {
		if (neggi >= 0.) {
		    iwhere[i__] = 2;
		}
	    } else {
		if (fabs(neggi) <= 0.) {
		    iwhere[i__] = -3;
		}
	    }
	}
	pointr = *head;
	if (iwhere[i__] != 0 && iwhere[i__] != -1) {
	    d__[i__] = 0.;
	} else {
	    d__[i__] = neggi;
	    f1 -= neggi * neggi;
	    i__2 = *col;
	    for (j = 1; j <= i__2; ++j) {
		p[j] += wy[i__ + pointr * wy_dim1] * neggi;
		p[*col + j] += ws[i__ + pointr * ws_dim1] * neggi;
		pointr = pointr % *m + 1;
	    }
	    if (nbd[i__] <= 2 && nbd[i__] != 0 && neggi < 0.) {
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tl / (-neggi);
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else if (nbd[i__] >= 2 && neggi > 0.) {
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tu / neggi;
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else {
		--nfree;
		iorder[nfree] = i__;
		if (fabs(neggi) > 0.) {
		    bnded = 0;
		}
	    }
	}
    }
    if (*theta != 1.) {
	dscal_(col, theta, &p[*col + 1], &c__1);
    }
    dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
    if (nbreak == 0 && nfree == *n + 1) {
	return 0;
    }
    i__1 = col2;
    for (j = 1; j <= i__1; ++j) {
	c__[j] = 0.;
    }
    f2 = -(*theta) * f1;
    f2_org__ = f2;
    if (*col > 0) {
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(&col2, &v[1], &c__1, &p[1], &c__1);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nseg = 1;
    if (nbreak == 0) {
	goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;
L777:
    tj0 = tj;
    if (iter == 1) {
	tj = bkmin;
	ibp = iorder[ibkmin];
    } else {
	if (iter == 2) {
	    if (ibkmin != nbreak) {
		t[ibkmin] = t[nbreak];
		iorder[ibkmin] = iorder[nbreak];
	    }
	}
	i__1 = iter - 2;
	hpsolb_(&nleft, &t[1], &iorder[1], &i__1);
	tj = t[nleft];
	ibp = iorder[nleft];
    }
    dt = tj - tj0;
    if (dtm < dt) {
	goto L888;
    }
    tsum += dt;
    --nleft;
    ++iter;
    dibp = d__[ibp];
    d__[ibp] = 0.;
    if (dibp > 0.) {
	zibp = u[ibp] - x[ibp];
	xcp[ibp] = u[ibp];
	iwhere[ibp] = 2;
    } else {
	zibp = l[ibp] - x[ibp];
	xcp[ibp] = l[ibp];
	iwhere[ibp] = 1;
    }
    if (nleft == 0 && nbreak == *n) {
	dtm = dt;
	goto L999;
    }
    ++(*nseg);
    d__1 = dibp;
    dibp2 = d__1 * d__1;
    f1 = f1 + dt * f2 + dibp2 - *theta * dibp * zibp;
    f2 -= *theta * dibp2;
    if (*col > 0) {
	daxpy_(&col2, &dt, &p[1], &c__1, &c__[1], &c__1);
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    wbp[j] = wy[ibp + pointr * wy_dim1];
	    wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
	    pointr = pointr % *m + 1;
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(&col2, &c__[1], &c__1, &v[1], &c__1);
	wmp = ddot_(&col2, &p[1], &c__1, &v[1], &c__1);
	wmw = ddot_(&col2, &wbp[1], &c__1, &v[1], &c__1);
	d__1 = -dibp;
	daxpy_(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }
    d__1 = *epsmch * f2_org__;
    f2 = max(d__1,f2);
    if (nleft > 0) {
	dtm = -f1 / f2;
	goto L777;
    } else if (bnded) {
	f1 = 0.;
	f2 = 0.;
	dtm = 0.;
    } else {
	dtm = -f1 / f2;
    }
L888:
    if (dtm <= 0.) {
	dtm = 0.;
    }
    tsum += dtm;
    daxpy_(n, &tsum, &d__[1], &c__1, &xcp[1], &c__1);
L999:
    if (*col > 0) {
	daxpy_(&col2, &dtm, &p[1], &c__1, &c__[1], &c__1);
    }

    return 0;
} /* cauchy_ */

static int mainlb_(int *n, int *m, double *x, 
		   double *l, double *u, int *nbd, double *f, double *g,
		   double *factr, double *pgtol, double *ws, double *wy,
		   double *sy, double *ss, double *wt, double *wn, 
		   double *snd, double *z__, double *r__, double *d__, 
		   double *t, double *xp, double *wa, int *index, 
		   int *iwhere, int *indx2, char *task,
		   char *csave, logical *lsave, int *isave, double *dsave)
{
    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	ss_dim1, ss_offset, wt_dim1, wt_offset, wn_dim1, wn_offset, 
	snd_dim1, snd_offset, i__1;
    double d__1, d__2;

    int i__, k;
    double gd, dr, rr, dtd;
    int col;
    double tol;
    logical wrk;
    double stp, cpu1, cpu2;
    int head;
    double fold;
    int nact;
    double ddum;
    int info, nseg;
    int nfgv, ifun, iter;
    char word[4];
    double time1, time2;
    int iback;
    double gdold;
    int nfree;
    logical boxed;
    int itail;
    double theta;
    double dnorm;
    int nskip, iword;
    double xstep, stpmx;
    int ileave;
    double cachyt;
    int itfile;
    double epsmch;
    logical updatd;
    double sbtime;
    logical prjctd;
    int iupdat;
    double sbgnrm;
    logical cnstnd;
    int nenter;
    double lnscht;
    int nintol;

    --indx2;
    --iwhere;
    --index;
    --xp;
    --t;
    --d__;
    --r__;
    --z__;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --wa;
    snd_dim1 = 2 * *m;
    snd_offset = 1 + snd_dim1;
    snd -= snd_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1;
    wn -= wn_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1;
    wt -= wt_offset;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1;
    ws -= ws_offset;
    --lsave;
    --isave;
    --dsave;

    if (strncmp(task, "START", 5) == 0) {
	epsmch = 2.22e-16;
	timer_(&time1);
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	iback = 0;
	itail = 0;
	iword = 0;
	nact = 0;
	ileave = 0;
	nenter = 0;
	fold = 0.;
	dnorm = 0.;
	cpu1 = 0.;
	gd = 0.;
	stpmx = 0.;
	sbgnrm = 0.;
	stp = 0.;
	gdold = 0.;
	dtd = 0.;
	iter = 0;
	nfgv = 0;
	nseg = 0;
	nintol = 0;
	nskip = 0;
	nfree = *n;
	ifun = 0;
	tol = *factr * epsmch;
	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;
	strcpy(word, "---");
	info = 0;
	itfile = 8;
	errclb_(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
	if (strncmp(task, "ERROR", 5) == 0) {
	    return 0;
	}
	active_(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], &prjctd, 
		&cnstnd, &boxed);
    } else {
	prjctd = lsave[1];
	cnstnd = lsave[2];
	boxed = lsave[3];
	updatd = lsave[4];
	nintol = isave[1];
	itfile = isave[3];
	iback = isave[4];
	nskip = isave[5];
	head = isave[6];
	col = isave[7];
	itail = isave[8];
	iter = isave[9];
	iupdat = isave[10];
	nseg = isave[12];
	nfgv = isave[13];
	info = isave[14];
	ifun = isave[15];
	iword = isave[16];
	nfree = isave[17];
	nact = isave[18];
	ileave = isave[19];
	nenter = isave[20];
	theta = dsave[1];
	fold = dsave[2];
	tol = dsave[3];
	dnorm = dsave[4];
	epsmch = dsave[5];
	cpu1 = dsave[6];
	cachyt = dsave[7];
	sbtime = dsave[8];
	lnscht = dsave[9];
	time1 = dsave[10];
	gd = dsave[11];
	stpmx = dsave[12];
	sbgnrm = dsave[13];
	stp = dsave[14];
	gdold = dsave[15];
	dtd = dsave[16];
	if (strncmp(task, "FG_LN", 5) == 0) {
	    goto L666;
	}
	if (strncmp(task, "NEW_X", 5) == 0) {
	    goto L777;
	}
	if (strncmp(task, "FG_ST", 5) == 0) {
	    goto L111;
	}
	if (strncmp(task, "STOP", 4) == 0) {
	    if (strncmp(task + 6, "CPU", 3) == 0) {
		dcopy_(n, &t[1], &c__1, &x[1], &c__1);
		dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
		*f = fold;
	    }
	    goto L999;
	}
    }
    strcpy(task, "FG_START");
    goto L1000;
L111:
    nfgv = 1;
    projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
    if (sbgnrm <= *pgtol) {
	strcpy(task, "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL");
	goto L999;
    }
L222:
    iword = -1;
    if (! cnstnd && col > 0) {
	dcopy_(n, &x[1], &c__1, &z__[1], &c__1);
	wrk = updatd;
	nseg = 0;
	goto L333;
    }
    timer_(&cpu1);
    cauchy_(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1],
	    &t[1], &d__[1], &z__[1], m, &wy[wy_offset], &ws[ws_offset],
	    &sy[sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1],
	    &wa[(*m << 1) + 1], &wa[(*m << 2) + 1], &wa[*m * 6 + 1], &nseg,
	    &sbgnrm, &info, &epsmch);
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	cachyt = cachyt + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    cachyt = cachyt + cpu2 - cpu1;
    nintol += nseg;
    freev_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iwhere[1],
	   &wrk, &updatd, &cnstnd, &iter);
    nact = *n - nfree;
L333:
    if (nfree == 0 || col == 0) {
	goto L555;
    }
    timer_(&cpu1);
    if (wrk) {
	formk_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iupdat,
	       &updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset],
	       &wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
    }
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
    cmprlb_(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset],
	    &wt[wt_offset], &z__[1], &r__[1], &wa[1], &index[1], &theta, &col,
	    &head, &nfree, &cnstnd, &info);
    if (info != 0) {
	goto L444;
    }
    subsm_(n, m, &nfree, &index[1], &l[1], &u[1], &nbd[1], &z__[1], &r__[1],
	   &xp[1], &ws[ws_offset], &wy[wy_offset], &theta, &x[1], &g[1], &col,
	   &head, &iword, &wa[1], &wn[wn_offset], &info);
L444:
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    sbtime = sbtime + cpu2 - cpu1;
L555:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = z__[i__] - x[i__];
    }
    timer_(&cpu1);
L666:
    lnsrlb_(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1],
	    &d__[1], &r__[1], &t[1], &z__[1], &stp, &dnorm, &dtd, &xstep,
	    &stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed,
	    &cnstnd, csave, &isave[22], &dsave[17]);
    if (info != 0 || iback >= 20) {
	dcopy_(n, &t[1], &c__1, &x[1], &c__1);
	dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
	*f = fold;
	if (col == 0) {
	    if (info == 0) {
		info = -9;
		--nfgv;
		--ifun;
		--iback;
	    }
	    strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
	    ++iter;
	    goto L999;
	} else {
	    if (info == 0) {
		--nfgv;
	    }
	    info = 0;
	    col = 0;
	    head = 1;
	    theta = 1.;
	    iupdat = 0;
	    updatd = 0;
	    strcpy(task, "RESTART_FROM_LNSRCH");
	    timer_(&cpu2);
	    lnscht = lnscht + cpu2 - cpu1;
	    goto L222;
	}
    } else if (strncmp(task, "FG_LN", 5) == 0) {
	goto L1000;
    } else {
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
	projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
	goto L1000;
    }
L777:
    if (sbgnrm <= *pgtol) {
	strcpy(task, "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL");
	goto L999;
    }
    d__1 = fabs(fold), d__2 = fabs(*f), d__1 = max(d__1, d__2);
    ddum = max(d__1,1.);
    if (fold - *f <= tol * ddum) {
	strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH");
	if (iback >= 10) {
	    info = -5;
	}
	goto L999;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = g[i__] - r__[i__];
    }
    rr = ddot_(n, &r__[1], &c__1, &r__[1], &c__1);
    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, &stp, &d__[1], &c__1);
	ddum = -gdold * stp;
    }
    if (dr <= epsmch * ddum) {
	++nskip;
	updatd = 0;
	goto L888;
    }
    updatd = 1;
    ++iupdat;
    matupd_(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset],
	    &ss[ss_offset], &d__[1], &r__[1], &itail, &iupdat, &col,
	    &head, &theta, &rr, &dr, &stp, &dtd);
    formt_(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col,
	   &theta, &info);
    if (info != 0) {
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = 0;
	goto L222;
    }
L888:
    goto L222;
L999:
    timer_(&time2);
#if 0    
    fprintf(stderr, "time = %g\n", time2 - time1);
#endif    
L1000:
    lsave[1] = prjctd;
    lsave[2] = cnstnd;
    lsave[3] = boxed;
    lsave[4] = updatd;
    isave[1] = nintol;
    isave[3] = itfile;
    isave[4] = iback;
    isave[5] = nskip;
    isave[6] = head;
    isave[7] = col;
    isave[8] = itail;
    isave[9] = iter;
    isave[10] = iupdat;
    isave[12] = nseg;
    isave[13] = nfgv;
    isave[14] = info;
    isave[15] = ifun;
    isave[16] = iword;
    isave[17] = nfree;
    isave[18] = nact;
    isave[19] = ileave;
    isave[20] = nenter;
    dsave[1] = theta;
    dsave[2] = fold;
    dsave[3] = tol;
    dsave[4] = dnorm;
    dsave[5] = epsmch;
    dsave[6] = cpu1;
    dsave[7] = cachyt;
    dsave[8] = sbtime;
    dsave[9] = lnscht;
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx;
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;
    return 0;
} /* mainlb_ */

int setulb_(int *n, int *m, double *x, 
	    double *l, double *u, int *nbd, double *f, double 
	    *g, double *factr, double *pgtol, double *wa, int *iwa,
	    char *task, char *csave, logical *lsave, 
	    int *isave, double *dsave)
{
    int i__1;

    int ld, lr, lt, lz, lwa, lwn, lss, lxp, lws, lwt, lsy, lwy, lsnd;

    --iwa;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --wa;
    --lsave;
    --isave;
    --dsave;

    if (strncmp(task, "START", 5) == 0) {
	isave[1] = *m * *n;
	i__1 = *m;
	isave[2] = i__1 * i__1;
	i__1 = *m;
	isave[3] = i__1 * i__1 << 2;
	isave[4] = 1;
	isave[5] = isave[4] + isave[1];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[2];
	isave[8] = isave[7] + isave[2];
	isave[9] = isave[8] + isave[2];
	isave[10] = isave[9] + isave[3];
	isave[11] = isave[10] + isave[3];
	isave[12] = isave[11] + *n;
	isave[13] = isave[12] + *n;
	isave[14] = isave[13] + *n;
	isave[15] = isave[14] + *n;
	isave[16] = isave[15] + *n;
    }
    
    lws = isave[4];
    lwy = isave[5];
    lsy = isave[6];
    lss = isave[7];
    lwt = isave[8];
    lwn = isave[9];
    lsnd = isave[10];
    lz = isave[11];
    lr = isave[12];
    ld = isave[13];
    lt = isave[14];
    lxp = isave[15];
    lwa = isave[16];
    
    mainlb_(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], factr, pgtol,
	    &wa[lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lwt], &wa[lwn], &wa[lsnd],
	    &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lxp], &wa[lwa], &iwa[1], 
	    &iwa[*n + 1], &iwa[(*n << 1) + 1], task, csave, &lsave[1],
	    &isave[22], &dsave[1]);
    return 0;
} /* setulb_ */
