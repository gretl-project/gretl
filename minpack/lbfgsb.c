/* 
   lbfgsb.c -- translated from lbfgsb.f (plus a couple of functions
   from linpack.f) via f2c, then subject to minimal clean-up by
   Allin Cottrell, 2015-07-26. Rewritten in idiomatic C (0-based
   indexing, callbacks in place of reverse communication, no BLAS),
   2026-10. The original source is version 3.0 of L-BFGS-B,
   downloaded from

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

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "minpack.h"

/* Basic vector operations (formerly calls to the BLAS) */

static inline double dot (int n, const double *x, const double *y)
{
    double s = 0.0;
    int i;

    for (i=0; i<n; i++) {
        s += x[i] * y[i];
    }

    return s;
}

/* y += a * x */

static inline void axpy (int n, double a, const double *x, double *y)
{
    int i;

    for (i=0; i<n; i++) {
        y[i] += a * x[i];
    }
}

/* x *= a */

static inline void scal (int n, double a, double *x)
{
    int i;

    for (i=0; i<n; i++) {
        x[i] *= a;
    }
}

static inline void copy (int n, const double *x, double *y)
{
    if (n > 0) {
        memcpy(y, x, n * sizeof *y);
    }
}

/* These mirror the f2c max and min macros, as opposed to fmax() and
   fmin(), which differ in their treatment of NaN.
*/

static inline double dmax (double a, double b)
{
    return a >= b ? a : b;
}

static inline double dmin (double a, double b)
{
    return a <= b ? a : b;
}

/* project x onto the feasible set */

static void project (int n, const double *l, const double *u,
                     const int *nbd, double *x)
{
    int i;

    for (i=0; i<n; i++) {
        if ((nbd[i] == 1 || nbd[i] == 2) && x[i] < l[i]) {
            x[i] = l[i];
        } else if ((nbd[i] == 2 || nbd[i] == 3) && x[i] > u[i]) {
            x[i] = u[i];
        }
    }
}

/* are f and all elements of g finite? */

static int fg_finite (int n, double f, const double *g)
{
    int i;

    if (!isfinite(f)) {
        return 0;
    }
    for (i=0; i<n; i++) {
        if (!isfinite(g[i])) {
            return 0;
        }
    }

    return 1;
}

/* Cholesky factorization of the symmetric positive definite n x n
   matrix @a (column-major, leading dimension @lda), as per LINPACK's
   dpofa: on exit the upper triangle of @a holds R such that A = R'R.
   Only the upper triangle is referenced. Returns 0 on success, or
   (j + 1) if the leading minor of order j + 1 is not positive
   definite.
*/

static int dpofa (double *a, int lda, int n)
{
    int j, k;

    for (j=0; j<n; j++) {
        double *aj = a + j*lda;
        double s = 0.0;

        for (k=0; k<j; k++) {
            double t = aj[k] - dot(k, a + k*lda, aj);

            t /= a[k + k*lda];
            aj[k] = t;
            s += t * t;
        }
        s = aj[j] - s;
        if (s <= 0.0) {
            return j + 1;
        }
        aj[j] = sqrt(s);
    }

    return 0;
}

/* The two cases of LINPACK's dtrsl that L-BFGS-B uses: solve
   R x = b (rsolve) or R' x = b (rtsolve) in place, where R is the
   n x n upper triangular matrix @r with leading dimension @ldr.
   Both return 0 on success, or (j + 1) if R(j,j) is zero.
*/

static int check_rdiag (const double *r, int ldr, int n)
{
    int j;

    for (j=0; j<n; j++) {
        if (r[j + j*ldr] == 0.0) {
            return j + 1;
        }
    }

    return 0;
}

static int rsolve (const double *r, int ldr, int n, double *b)
{
    int j, info;

    if (n == 0 || (info = check_rdiag(r, ldr, n)) != 0) {
        return n == 0 ? 0 : info;
    }
    b[n-1] /= r[(n-1) + (n-1)*ldr];
    for (j=n-2; j>=0; j--) {
        axpy(j + 1, -b[j+1], r + (j+1)*ldr, b);
        b[j] /= r[j + j*ldr];
    }

    return 0;
}

static int rtsolve (const double *r, int ldr, int n, double *b)
{
    int j, info;

    if (n == 0 || (info = check_rdiag(r, ldr, n)) != 0) {
        return n == 0 ? 0 : info;
    }
    b[0] /= r[0];
    for (j=1; j<n; j++) {
        b[j] -= dot(j, r + j*ldr, b);
        b[j] /= r[j + j*ldr];
    }

    return 0;
}

/* Initialize iwhere and project the initial x onto the feasible set.
   On return, *prjctd indicates whether x was changed by projection,
   *cnstnd whether any variable is bounded, and *boxed whether all
   variables are bounded on both sides.
*/

static void active (int n, const double *l, const double *u,
                    const int *nbd, double *x, int *iwhere,
                    int *prjctd, int *cnstnd, int *boxed)
{
    int i;

    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;

    for (i=0; i<n; i++) {
        if (nbd[i] > 0) {
            if (nbd[i] <= 2 && x[i] <= l[i]) {
                if (x[i] < l[i]) {
                    *prjctd = 1;
                    x[i] = l[i];
                }
            } else if (nbd[i] >= 2 && x[i] >= u[i]) {
                if (x[i] > u[i]) {
                    *prjctd = 1;
                    x[i] = u[i];
                }
            }
        }
    }

    for (i=0; i<n; i++) {
        if (nbd[i] != 2) {
            *boxed = 0;
        }
        if (nbd[i] == 0) {
            /* this variable is always free */
            iwhere[i] = -1;
        } else {
            *cnstnd = 1;
            if (nbd[i] == 2 && u[i] - l[i] <= 0.0) {
                /* this variable is always fixed */
                iwhere[i] = 3;
            } else {
                iwhere[i] = 0;
            }
        }
    }
}

/* Compute the product of the 2m x 2m middle matrix of the compact
   L-BFGS formula with the 2col vector @v, writing the result into @p.
   @sy holds S'Y and @wt the Cholesky factor of (theta*S'S + LD^(-1)L'),
   both m x m. Returns 0 on success, non-zero if a triangular solve
   fails.
*/

static int bmv (int m, const double *sy, const double *wt, int col,
                const double *v, double *p)
{
    double sum;
    int i, k, info;

    if (col == 0) {
        return 0;
    }

    /* solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
             [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ]
       by solving J p2 = v2 + L D^(-1) v1 */
    p[col] = v[col];
    for (i=1; i<col; i++) {
        sum = 0.0;
        for (k=0; k<i; k++) {
            sum += sy[i + k*m] * v[k] / sy[k + k*m];
        }
        p[col + i] = v[col + i] + sum;
    }
    info = rtsolve(wt, m, col, p + col);
    if (info != 0) {
        return info;
    }

    /* solve D^(1/2) p1 = v1 */
    for (i=0; i<col; i++) {
        p[i] = v[i] / sqrt(sy[i + i*m]);
    }

    /* solve [ -D^(1/2)   D^(-1/2)*L' ] [ p1 ] = [ p1 ]
             [  0         J'          ] [ p2 ]   [ p2 ]
       by solving J' p2 = p2 and then D^(1/2) p1 = -p1 + D^(-1/2) L' p2 */
    info = rsolve(wt, m, col, p + col);
    if (info != 0) {
        return info;
    }
    for (i=0; i<col; i++) {
        p[i] = -p[i] / sqrt(sy[i + i*m]);
    }
    for (i=0; i<col; i++) {
        sum = 0.0;
        for (k=i+1; k<col; k++) {
            sum += sy[k + i*m] * p[col + k] / sy[i + i*m];
        }
        p[i] += sum;
    }

    return 0;
}

/* Compute r = -Z'B(xcp - xk) - Z'g, using wa[2m..4m-1] from cauchy().
   The ring buffers @ws and @wy (n x m) hold S and Y, with the oldest
   correction in column @head. Returns 0 on success or -8 on failure.
*/

static int cmprlb (int n, int m, const double *x, const double *g,
                   const double *ws, const double *wy, const double *sy,
                   const double *wt, const double *z, double *r,
                   double *wa, const int *index, double theta, int col,
                   int head, int nfree, int cnstnd)
{
    int i, j, k, pointr;

    if (!cnstnd && col > 0) {
        for (i=0; i<n; i++) {
            r[i] = -g[i];
        }
    } else {
        for (i=0; i<nfree; i++) {
            k = index[i];
            r[i] = -theta * (z[k] - x[k]) - g[k];
        }
        if (bmv(m, sy, wt, col, wa + 2*m, wa) != 0) {
            return -8;
        }
        pointr = head;
        for (j=0; j<col; j++) {
            const double *wyj = wy + pointr*n;
            const double *wsj = ws + pointr*n;
            double a1 = wa[j];
            double a2 = theta * wa[col + j];

            for (i=0; i<nfree; i++) {
                k = index[i];
                r[i] = r[i] + wyj[k] * a1 + wsj[k] * a2;
            }
            pointr = (pointr + 1) % m;
        }
    }

    return 0;
}

/* Check the input arguments for errors; if any are found, write a
   message into @task (if there are several errors, the last one
   found is reported).
*/

static void errclb (int n, int m, double factr, const double *l,
                    const double *u, const int *nbd, char *task)
{
    int i;

    if (n <= 0) {
        strcpy(task, "ERROR: N .LE. 0");
    }
    if (m <= 0) {
        strcpy(task, "ERROR: M .LE. 0");
    }
    if (factr < 0.0) {
        strcpy(task, "ERROR: FACTR .LT. 0");
    }
    for (i=0; i<n; i++) {
        if (nbd[i] < 0 || nbd[i] > 3) {
            strcpy(task, "ERROR: INVALID NBD");
        }
        if (nbd[i] == 2 && l[i] > u[i]) {
            strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
        }
    }
}

/* Form the LEL' factorization of the indefinite matrix

     K = [-D - Y'ZZ'Y/theta     L_a' - R_z'  ]
         [L_a - R_z             theta*S'AA'S ]

   where E = [-I 0; 0 I], storing the factors in the upper triangle
   of @wn (2m x 2m). The matrix @wn1 (2m x 2m) holds the inner
   products of S and Y restricted to the free and active variables,
   and is updated incrementally from one iteration to the next using
   the entering (indx2[0..nenter-1]) and leaving (indx2[ileave..n-1])
   variables. @ind holds the free variables in its first @nsub
   elements, the active ones in the rest. Returns 0 on success, or -1
   or -2 if the first or second Cholesky factorization fails.
*/

#define WN(i,j) wn[(i) + (j)*m2]
#define WN1(i,j) wn1[(i) + (j)*m2]

static int formk (int n, int nsub, const int *ind, int nenter,
                  int ileave, const int *indx2, int iupdat, int updatd,
                  double *wn, double *wn1, int m, const double *ws,
                  const double *wy, const double *sy, double theta,
                  int col, int head)
{
    const double *wyi, *wyj, *wsi, *wsj;
    double temp1, temp2, temp3, temp4;
    int i, k, k1, is, js, iy, jy, is1, js1;
    int ipntr, jpntr, upcl;
    int m2 = 2*m, col2 = 2*col;

    /* Form the lower triangular part of
         WN1 = [Y' ZZ'Y   L_a'+R_z']
               [L_a+R_z   S'AA'S   ]
       where L_a is the strictly lower triangular part of S'AA'Y
       and R_z is the upper triangular part of S'ZZ'Y.
    */

    if (updatd) {
        if (iupdat > m) {
            /* shift old part of WN1 */
            for (jy=0; jy<m-1; jy++) {
                js = m + jy;
                copy(m - jy - 1, &WN1(jy+1, jy+1), &WN1(jy, jy));
                copy(m - jy - 1, &WN1(js+1, js+1), &WN1(js, js));
                copy(m - 1, &WN1(m+1, jy+1), &WN1(m, jy));
            }
        }

        /* put new rows in blocks (1,1), (2,1) and (2,2): free
           variables are ind[0..nsub-1], active ones ind[nsub..n-1] */
        iy = col - 1;
        is = m + col - 1;
        ipntr = head + col - 1;
        if (ipntr >= m) {
            ipntr -= m;
        }
        wyi = wy + ipntr*n;
        wsi = ws + ipntr*n;
        jpntr = head;
        for (jy=0; jy<col; jy++) {
            js = m + jy;
            wyj = wy + jpntr*n;
            wsj = ws + jpntr*n;
            temp1 = temp2 = temp3 = 0.0;
            for (k=0; k<nsub; k++) {
                k1 = ind[k];
                temp1 += wyi[k1] * wyj[k1];
            }
            for (k=nsub; k<n; k++) {
                k1 = ind[k];
                temp2 += wsi[k1] * wsj[k1];
                temp3 += wsi[k1] * wyj[k1];
            }
            WN1(iy, jy) = temp1;
            WN1(is, js) = temp2;
            WN1(is, jy) = temp3;
            jpntr = (jpntr + 1) % m;
        }

        /* put new column in block (2,1) */
        jy = col - 1;
        jpntr = head + col - 1;
        if (jpntr >= m) {
            jpntr -= m;
        }
        wyj = wy + jpntr*n;
        ipntr = head;
        for (i=0; i<col; i++) {
            is = m + i;
            wsi = ws + ipntr*n;
            temp3 = 0.0;
            for (k=0; k<nsub; k++) {
                k1 = ind[k];
                temp3 += wsi[k1] * wyj[k1];
            }
            ipntr = (ipntr + 1) % m;
            WN1(is, jy) = temp3;
        }
        upcl = col - 1;
    } else {
        upcl = col;
    }

    /* modify the old parts in blocks (1,1) and (2,2) due to
       changes in the set of free variables */
    ipntr = head;
    for (iy=0; iy<upcl; iy++) {
        is = m + iy;
        wyi = wy + ipntr*n;
        wsi = ws + ipntr*n;
        jpntr = head;
        for (jy=0; jy<=iy; jy++) {
            js = m + jy;
            wyj = wy + jpntr*n;
            wsj = ws + jpntr*n;
            temp1 = temp2 = temp3 = temp4 = 0.0;
            for (k=0; k<nenter; k++) {
                k1 = indx2[k];
                temp1 += wyi[k1] * wyj[k1];
                temp2 += wsi[k1] * wsj[k1];
            }
            for (k=ileave; k<n; k++) {
                k1 = indx2[k];
                temp3 += wyi[k1] * wyj[k1];
                temp4 += wsi[k1] * wsj[k1];
            }
            WN1(iy, jy) = WN1(iy, jy) + temp1 - temp3;
            WN1(is, js) = WN1(is, js) - temp2 + temp4;
            jpntr = (jpntr + 1) % m;
        }
        ipntr = (ipntr + 1) % m;
    }

    /* modify the old parts in block (2,1) */
    ipntr = head;
    for (is=m; is<m+upcl; is++) {
        wsi = ws + ipntr*n;
        jpntr = head;
        for (jy=0; jy<upcl; jy++) {
            wyj = wy + jpntr*n;
            temp1 = temp3 = 0.0;
            for (k=0; k<nenter; k++) {
                k1 = indx2[k];
                temp1 += wsi[k1] * wyj[k1];
            }
            for (k=ileave; k<n; k++) {
                k1 = indx2[k];
                temp3 += wsi[k1] * wyj[k1];
            }
            if (is <= jy + m) {
                WN1(is, jy) = WN1(is, jy) + temp1 - temp3;
            } else {
                WN1(is, jy) = WN1(is, jy) - temp1 + temp3;
            }
            jpntr = (jpntr + 1) % m;
        }
        ipntr = (ipntr + 1) % m;
    }

    /* Form the upper triangle of WN = [D+Y'ZZ'Y/theta   -L_a'+R_z' ]
                                       [-L_a +R_z        S'AA'S*theta]
    */
    for (iy=0; iy<col; iy++) {
        is = col + iy;
        is1 = m + iy;
        for (jy=0; jy<=iy; jy++) {
            js = col + jy;
            js1 = m + jy;
            WN(jy, iy) = WN1(iy, jy) / theta;
            WN(js, is) = WN1(is1, js1) * theta;
        }
        for (jy=0; jy<iy; jy++) {
            WN(jy, is) = -WN1(is1, jy);
        }
        for (jy=iy; jy<col; jy++) {
            WN(jy, is) = WN1(is1, jy);
        }
        WN(iy, iy) += sy[iy + iy*m];
    }

    /* Form the upper triangle of WN = [  LL'            L^-1(-L_a'+R_z')]
                                       [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
       first Cholesky factor (1,1) block of wn to get LL'
       with L' stored in the upper triangle of wn */
    if (dpofa(wn, m2, col) != 0) {
        return -1;
    }
    /* then form L^-1(-L_a'+R_z') in the (1,2) block */
    for (js=col; js<col2; js++) {
        rtsolve(wn, m2, col, &WN(0, js));
    }

    /* Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
       upper triangle of (2,2) block of wn */
    for (is=col; is<col2; is++) {
        for (js=is; js<col2; js++) {
            WN(is, js) += dot(col, &WN(0, is), &WN(0, js));
        }
    }

    /* Cholesky factorization of (2,2) block of wn */
    if (dpofa(&WN(col, col), m2, col) != 0) {
        return -2;
    }

    return 0;
}

#undef WN
#undef WN1

/* Form the upper half of the positive definite symmetric matrix
   T = theta*S'S + L D^(-1) L', store it in the upper triangle of
   @wt, then Cholesky-factorize it in place. Returns 0 on success
   or -3 if T is not positive definite.
*/

static int formt (int m, double *wt, const double *sy,
                  const double *ss, int col, double theta)
{
    double ddum;
    int i, j, k;

    for (j=0; j<col; j++) {
        wt[j*m] = theta * ss[j*m];
    }
    for (i=1; i<col; i++) {
        for (j=i; j<col; j++) {
            ddum = 0.0;
            for (k=0; k<i; k++) {
                ddum += sy[i + k*m] * sy[j + k*m] / sy[k + k*m];
            }
            wt[i + j*m] = ddum + theta * ss[i + j*m];
        }
    }

    return dpofa(wt, m, col) != 0 ? -3 : 0;
}

/* Count the entering and leaving variables for the subspace
   minimization at iteration @iter, and set up the index sets of
   free and active variables at the generalized Cauchy point.

   On return, index[0..nfree-1] holds the free variables and
   index[nfree..n-1] the active ones; indx2[0..nenter-1] holds the
   variables that have entered the free set since the last
   iteration, and indx2[ileave..n-1] those that have left it.
   Returns non-zero if the matrix K must be recomputed by formk().
*/

static int freev (int n, int *nfree, int *index, int *nenter,
                  int *ileave, int *indx2, const int *iwhere,
                  int updatd, int cnstnd, int iter)
{
    int i, k, iact;

    *nenter = 0;
    *ileave = n;

    if (iter > 0 && cnstnd) {
        for (i=0; i<*nfree; i++) {
            k = index[i];
            if (iwhere[k] > 0) {
                indx2[--(*ileave)] = k;
            }
        }
        for (i=*nfree; i<n; i++) {
            k = index[i];
            if (iwhere[k] <= 0) {
                indx2[(*nenter)++] = k;
            }
        }
    }

    *nfree = 0;
    iact = n;
    for (i=0; i<n; i++) {
        if (iwhere[i] <= 0) {
            index[(*nfree)++] = i;
        } else {
            index[--iact] = i;
        }
    }

    return *ileave < n || *nenter > 0 || updatd;
}

/* Heap operations on the breakpoints @t (with companion indices
   @iorder) of length n. If @iheap is zero, first arrange t into a
   heap (t[0] the least element). Then move the least element to
   t[n-1] and restore the heap property on t[0..n-2].
*/

static void hpsolb (int n, double *t, int *iorder, int iheap)
{
    double out, ddum;
    int i, j, k;
    int indxin, indxou;

    if (iheap == 0) {
        for (k=1; k<n; k++) {
            ddum = t[k];
            indxin = iorder[k];
            /* sift t[k] up toward the root */
            i = k;
            while (i > 0) {
                j = (i - 1) / 2;
                if (ddum < t[j]) {
                    t[i] = t[j];
                    iorder[i] = iorder[j];
                    i = j;
                } else {
                    break;
                }
            }
            t[i] = ddum;
            iorder[i] = indxin;
        }
    }

    if (n > 1) {
        /* remove the least element t[0] and sift t[n-1] down
           from the root into the heap t[0..n-2] */
        i = 0;
        out = t[0];
        indxou = iorder[0];
        ddum = t[n-1];
        indxin = iorder[n-1];
        while (1) {
            j = 2*i + 1;
            if (j >= n - 1) {
                break;
            }
            if (t[j+1] < t[j]) {
                j++;
            }
            if (t[j] < ddum) {
                t[i] = t[j];
                iorder[i] = iorder[j];
                i = j;
            } else {
                break;
            }
        }
        t[i] = ddum;
        iorder[i] = indxin;
        /* put the least member in t[n-1] */
        t[n-1] = out;
        iorder[n-1] = indxou;
    }
}

/* Update the limited-memory matrices WS and WY (ring buffers of the
   last m correction pairs, oldest in column @head) and form the
   middle matrices SY = S'Y and SS = S'S, given the new correction
   pair (d, r) and the step @stp.
*/

static void matupd (int n, int m, double *ws, double *wy, double *sy,
                    double *ss, const double *d, const double *r,
                    int *itail, int iupdat, int *col, int *head,
                    double *theta, double rr, double dr, double stp,
                    double dtd)
{
    int j, last, pointr;

    /* set pointers for matrices WS and WY */
    if (iupdat <= m) {
        *col = iupdat;
        *itail = (*head + iupdat - 1) % m;
    } else {
        *itail = (*itail + 1) % m;
        *head = (*head + 1) % m;
    }

    /* update matrices WS and WY */
    copy(n, d, ws + *itail*n);
    copy(n, r, wy + *itail*n);

    /* set theta = yy/ys */
    *theta = rr / dr;

    /* Form the middle matrix in B: update the upper triangle of SS
       and the lower triangle of SY */
    if (iupdat > m) {
        /* move old information */
        for (j=0; j<*col-1; j++) {
            copy(j + 1, &ss[1 + (j+1)*m], &ss[j*m]);
            copy(*col - j - 1, &sy[(j+1) + (j+1)*m], &sy[j + j*m]);
        }
    }

    /* add new information: the last row of SY and the last column
       of SS */
    last = *col - 1;
    pointr = *head;
    for (j=0; j<last; j++) {
        sy[last + j*m] = dot(n, d, wy + pointr*n);
        ss[j + last*m] = dot(n, ws + pointr*n, d);
        pointr = (pointr + 1) % m;
    }
    if (stp == 1.0) {
        ss[last + last*m] = dtd;
    } else {
        ss[last + last*m] = stp * stp * dtd;
    }
    sy[last + last*m] = dr;
}

/* Return the infinity norm of the projected gradient */

static double projgr (int n, const double *l, const double *u,
                      const int *nbd, const double *x, const double *g)
{
    double gi, sbgnrm = 0.0;
    int i;

    for (i=0; i<n; i++) {
        gi = g[i];
        if (nbd[i] != 0) {
            if (gi < 0.0) {
                if (nbd[i] >= 2) {
                    gi = dmax(x[i] - u[i], gi);
                }
            } else if (nbd[i] <= 2) {
                gi = dmin(x[i] - l[i], gi);
            }
        }
        sbgnrm = dmax(sbgnrm, fabs(gi));
    }

    return sbgnrm;
}

/* Subspace minimization: given xcp (passed in as @x) and the reduced
   gradient @d (length nsub), compute the minimizer of the quadratic
   model on the subspace of free variables ind[0..nsub-1], writing it
   into @x, then project it onto the feasible set (falling back on
   a backtrack toward xcp if the projected point is not a descent
   direction relative to @xx with gradient @gg). Returns 0 on
   success or non-zero if a triangular solve failed.
*/

static int subsm (int n, int m, int nsub, const int *ind,
                  const double *l, const double *u, const int *nbd,
                  double *x, double *d, double *xp, const double *ws,
                  const double *wy, double theta, const double *xx,
                  const double *gg, int col, int head,
                  double *wv, const double *wn)
{
    double alpha, dk, xk, dd_p, temp1, temp2;
    int i, j, k, js, jy, ibd, pointr, info;
    int m2 = 2*m, col2 = 2*col;
    int iword;

    if (nsub <= 0) {
        return 0;
    }

    /* compute wv = W'Zd */
    pointr = head;
    for (i=0; i<col; i++) {
        const double *wyi = wy + pointr*n;
        const double *wsi = ws + pointr*n;

        temp1 = temp2 = 0.0;
        for (j=0; j<nsub; j++) {
            k = ind[j];
            temp1 += wyi[k] * d[j];
            temp2 += wsi[k] * d[j];
        }
        wv[i] = temp1;
        wv[col + i] = theta * temp2;
        pointr = (pointr + 1) % m;
    }

    /* compute wv := K^(-1) wv */
    info = rtsolve(wn, m2, col2, wv);
    if (info != 0) {
        return info;
    }
    for (i=0; i<col; i++) {
        wv[i] = -wv[i];
    }
    info = rsolve(wn, m2, col2, wv);
    if (info != 0) {
        return info;
    }

    /* compute d = (1/theta)d + (1/theta**2)Z'W wv */
    pointr = head;
    for (jy=0; jy<col; jy++) {
        const double *wyj = wy + pointr*n;
        const double *wsj = ws + pointr*n;

        js = col + jy;
        for (i=0; i<nsub; i++) {
            k = ind[i];
            d[i] = d[i] + wyj[k] * wv[jy] / theta + wsj[k] * wv[js];
        }
        pointr = (pointr + 1) % m;
    }
    scal(nsub, 1.0 / theta, d);

    /* try the projection of d */
    iword = 0;
    copy(n, x, xp);
    for (i=0; i<nsub; i++) {
        k = ind[i];
        dk = d[i];
        xk = x[k];
        if (nbd[k] == 0) {
            /* free variable */
            x[k] = xk + dk;
        } else if (nbd[k] == 1) {
            /* lower bound only */
            x[k] = dmax(l[k], xk + dk);
            if (x[k] == l[k]) {
                iword = 1;
            }
        } else if (nbd[k] == 2) {
            /* upper and lower bounds */
            xk = dmax(l[k], xk + dk);
            x[k] = dmin(u[k], xk);
            if (x[k] == l[k] || x[k] == u[k]) {
                iword = 1;
            }
        } else if (nbd[k] == 3) {
            /* upper bound only */
            x[k] = dmin(u[k], xk + dk);
            if (x[k] == u[k]) {
                iword = 1;
            }
        }
    }

    if (iword == 0) {
        return 0;
    }

    /* check sign of the directional derivative */
    dd_p = 0.0;
    for (i=0; i<n; i++) {
        dd_p += (x[i] - xx[i]) * gg[i];
    }
    if (dd_p <= 0.0) {
        return 0;
    }

    /* not a descent direction: restore x and backtrack along d
       to the feasible region */
    copy(n, xp, x);
    alpha = 1.0;
    temp1 = alpha;
    ibd = 0;
    for (i=0; i<nsub; i++) {
        k = ind[i];
        dk = d[i];
        if (nbd[k] != 0) {
            if (dk < 0.0 && nbd[k] <= 2) {
                temp2 = l[k] - x[k];
                if (temp2 >= 0.0) {
                    temp1 = 0.0;
                } else if (dk * alpha < temp2) {
                    temp1 = temp2 / dk;
                }
            } else if (dk > 0.0 && nbd[k] >= 2) {
                temp2 = u[k] - x[k];
                if (temp2 <= 0.0) {
                    temp1 = 0.0;
                } else if (dk * alpha > temp2) {
                    temp1 = temp2 / dk;
                }
            }
            if (temp1 < alpha) {
                alpha = temp1;
                ibd = i;
            }
        }
    }

    if (alpha < 1.0) {
        dk = d[ibd];
        k = ind[ibd];
        if (dk > 0.0) {
            x[k] = u[k];
            d[ibd] = 0.0;
        } else if (dk < 0.0) {
            x[k] = l[k];
            d[ibd] = 0.0;
        }
    }
    for (i=0; i<nsub; i++) {
        k = ind[i];
        x[k] += alpha * d[i];
    }

    return 0;
}

/* Compute a safeguarded step for a line search and update the
   interval that contains a step satisfying a sufficient decrease and
   a curvature condition (More' and Thuente). On input, (stx, fx, dx)
   is the step with the least function value so far, (sty, fy, dy)
   the other endpoint of the interval of uncertainty, and stp the
   current step, with function value fp and derivative dp; all are
   updated on output, as is the flag *brackt which records whether
   a minimizer has been bracketed. The new step is confined to
   [stpmin, stpmax].
*/

static void dcstep (double *pstx, double *pfx, double *pdx,
                    double *psty, double *pfy, double *pdy,
                    double *pstp, double fp, double dp, int *brackt,
                    double stpmin, double stpmax)
{
    double stx = *pstx, fx = *pfx, dx = *pdx;
    double sty = *psty, fy = *pfy, dy = *pdy;
    double stp = *pstp;
    double p, q, r, s, ts, sgnd, stpc, stpf, stpq, gamma, theta;

    sgnd = dp * (dx / fabs(dx));

    if (fp > fx) {
        /* case 1: a higher function value; the minimum is bracketed */
        theta = (fx - fp) * 3. / (stp - stx) + dx + dp;
        s = dmax(dmax(fabs(theta), fabs(dx)), fabs(dp));
        ts = theta / s;
        gamma = s * sqrt(ts * ts - dx / s * (dp / s));
        if (stp < stx) {
            gamma = -gamma;
        }
        p = gamma - dx + theta;
        q = gamma - dx + gamma + dp;
        r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2. * (stp - stx);
        if (fabs(stpc - stx) < fabs(stpq - stx)) {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.;
        }
        *brackt = 1;
    } else if (sgnd < 0.0) {
        /* case 2: lower function value and derivatives of opposite
           sign; the minimum is bracketed */
        theta = (fx - fp) * 3. / (stp - stx) + dx + dp;
        s = dmax(dmax(fabs(theta), fabs(dx)), fabs(dp));
        ts = theta / s;
        gamma = s * sqrt(ts * ts - dx / s * (dp / s));
        if (stp > stx) {
            gamma = -gamma;
        }
        p = gamma - dp + theta;
        q = gamma - dp + gamma + dx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + dp / (dp - dx) * (stx - stp);
        if (fabs(stpc - stp) > fabs(stpq - stp)) {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = 1;
    } else if (fabs(dp) < fabs(dx)) {
        /* case 3: lower function value, derivatives of the same sign,
           and the magnitude of the derivative decreases */
        theta = (fx - fp) * 3. / (stp - stx) + dx + dp;
        s = dmax(dmax(fabs(theta), fabs(dx)), fabs(dp));
        ts = theta / s;
        gamma = s * sqrt(dmax(0., ts * ts - dx / s * (dp / s)));
        if (stp > stx) {
            gamma = -gamma;
        }
        p = gamma - dp + theta;
        q = gamma + (dx - dp) + gamma;
        r = p / q;
        if (r < 0.0 && gamma != 0.0) {
            stpc = stp + r * (stx - stp);
        } else if (stp > stx) {
            stpc = stpmax;
        } else {
            stpc = stpmin;
        }
        stpq = stp + dp / (dp - dx) * (stx - stp);
        if (*brackt) {
            /* a minimizer has been bracketed */
            if (fabs(stpc - stp) < fabs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            if (stp > stx) {
                stpf = dmin(stp + (sty - stp) * .66, stpf);
            } else {
                stpf = dmax(stp + (sty - stp) * .66, stpf);
            }
        } else {
            if (fabs(stpc - stp) > fabs(stpq - stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            stpf = dmin(stpmax, stpf);
            stpf = dmax(stpmin, stpf);
        }
    } else {
        /* case 4: lower function value, derivatives of the same sign,
           and the magnitude of the derivative does not decrease */
        if (*brackt) {
            theta = (fp - fy) * 3. / (sty - stp) + dy + dp;
            s = dmax(dmax(fabs(theta), fabs(dy)), fabs(dp));
            ts = theta / s;
            gamma = s * sqrt(ts * ts - dy / s * (dp / s));
            if (stp > sty) {
                gamma = -gamma;
            }
            p = gamma - dp + theta;
            q = gamma - dp + gamma + dy;
            r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        } else if (stp > stx) {
            stpf = stpmax;
        } else {
            stpf = stpmin;
        }
    }

    /* update the interval which contains a minimizer */
    if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if (sgnd < 0.0) {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }

    *pstx = stx; *pfx = fx; *pdx = dx;
    *psty = sty; *pfy = fy; *pdy = dy;
    *pstp = stpf;
}

/* State of the line search, carried over between calls to dcsrch()
   and from lnsrlb() to the main loop */

typedef struct {
    double stp;    /* current step length */
    double fold;   /* f at the start of the line search */
    double gd;     /* directional derivative g'd at stp */
    double gdold;  /* directional derivative at the start */
    double dtd;    /* squared norm of the search direction d */
    /* internal state of dcsrch() */
    int brackt, stage;
    double ginit, gtest, gx, gy, finit, fx, fy;
    double stx, sty, stmin, stmax, width, width1;
} line_search;

/* dcsrch() tasks */
enum { LS_START, LS_FG, LS_CONV, LS_WARN, LS_ERROR };

/* Find a step that satisfies a sufficient decrease condition and a
   curvature condition (More' and Thuente, "Line search algorithms
   with guaranteed sufficient decrease", ACM TOMS 20, 1994).

   This is a reverse-communication routine: it is first called with
   task LS_START, with @f and @g the function value and derivative at
   step zero and *stp an initial estimate of a satisfactory step. It
   then returns LS_FG with a new trial step in *stp, and should be
   called again (passing that task back) with @f and @g evaluated at
   the trial step; this is repeated until it returns LS_CONV (a
   satisfactory step was found) or LS_WARN (no further progress is
   possible; *stp is the best step found). LS_ERROR indicates invalid
   input on a call with LS_START.
*/

static int dcsrch (double f, double g, double *stp, double ftol,
                   double gtol, double xtol, double stpmin,
                   double stpmax, int task, line_search *s)
{
    double ftest, fm, fxm, fym, gm, gxm, gym;
    int ret = LS_FG;

    if (task == LS_START) {
        /* check the input arguments for errors */
        if (*stp < stpmin || *stp > stpmax || g >= 0.0 || ftol < 0.0 ||
            gtol < 0.0 || xtol < 0.0 || stpmin < 0.0 || stpmax < stpmin) {
            return LS_ERROR;
        }
        /* initialize local variables: the variables stx, fx, gx
           contain the values of the step, function, and derivative
           at the best step; sty, fy, gy those at the other endpoint
           of the interval of uncertainty; stmin, stmax are the
           minimum and maximum steps for the current interval */
        s->brackt = 0;
        s->stage = 1;
        s->finit = f;
        s->ginit = g;
        s->gtest = ftol * s->ginit;
        s->width = stpmax - stpmin;
        s->width1 = s->width / .5;
        s->stx = 0.0;
        s->fx = s->finit;
        s->gx = s->ginit;
        s->sty = 0.0;
        s->fy = s->finit;
        s->gy = s->ginit;
        s->stmin = 0.0;
        s->stmax = *stp + *stp * 4.;
        return LS_FG;
    }

    /* if psi(stp) <= 0 and f'(stp) >= 0 for some step, the algorithm
       enters the second stage */
    ftest = s->finit + *stp * s->gtest;
    if (s->stage == 1 && f <= ftest && g >= 0.0) {
        s->stage = 2;
    }

    /* test for warnings */
    if (s->brackt && (*stp <= s->stmin || *stp >= s->stmax)) {
        /* rounding errors prevent progress */
        ret = LS_WARN;
    }
    if (s->brackt && s->stmax - s->stmin <= xtol * s->stmax) {
        /* xtol test satisfied */
        ret = LS_WARN;
    }
    if (*stp == stpmax && f <= ftest && g <= s->gtest) {
        /* stp = stpmax */
        ret = LS_WARN;
    }
    if (*stp == stpmin && (f > ftest || g >= s->gtest)) {
        /* stp = stpmin */
        ret = LS_WARN;
    }

    /* test for convergence */
    if (f <= ftest && fabs(g) <= gtol * (-s->ginit)) {
        ret = LS_CONV;
    }

    if (ret != LS_FG) {
        return ret;
    }

    /* A modified function is used to predict the step during the
       first stage if a lower function value has been obtained but
       the decrease is not sufficient */
    if (s->stage == 1 && f <= s->fx && f > ftest) {
        /* define the modified function and derivative values */
        fm = f - *stp * s->gtest;
        fxm = s->fx - s->stx * s->gtest;
        fym = s->fy - s->sty * s->gtest;
        gm = g - s->gtest;
        gxm = s->gx - s->gtest;
        gym = s->gy - s->gtest;
        /* compute the new step and update the interval */
        dcstep(&s->stx, &fxm, &gxm, &s->sty, &fym, &gym, stp, fm, gm,
               &s->brackt, s->stmin, s->stmax);
        /* reset the function and derivative values for f */
        s->fx = fxm + s->stx * s->gtest;
        s->fy = fym + s->sty * s->gtest;
        s->gx = gxm + s->gtest;
        s->gy = gym + s->gtest;
    } else {
        dcstep(&s->stx, &s->fx, &s->gx, &s->sty, &s->fy, &s->gy, stp,
               f, g, &s->brackt, s->stmin, s->stmax);
    }

    /* decide if a bisection step is needed */
    if (s->brackt) {
        if (fabs(s->sty - s->stx) >= s->width1 * .66) {
            *stp = s->stx + (s->sty - s->stx) * .5;
        }
        s->width1 = s->width;
        s->width = fabs(s->sty - s->stx);
    }

    /* set the minimum and maximum steps allowed for stp */
    if (s->brackt) {
        s->stmin = dmin(s->stx, s->sty);
        s->stmax = dmax(s->stx, s->sty);
    } else {
        s->stmin = *stp + (*stp - s->stx) * 1.1;
        s->stmax = *stp + (*stp - s->stx) * 4.;
    }

    /* force the step to be within the bounds stpmax and stpmin */
    *stp = dmax(*stp, stpmin);
    *stp = dmin(*stp, stpmax);

    /* if further progress is not possible, let stp be the best
       point obtained during the search */
    if ((s->brackt && (*stp <= s->stmin || *stp >= s->stmax)) ||
        (s->brackt && s->stmax - s->stmin <= xtol * s->stmax)) {
        *stp = s->stx;
    }

    return LS_FG;
}

/* lnsrlb() return values */
enum { LNSRCH_OK, LNSRCH_ASCENT, LNSRCH_MAXBACK, LNSRCH_ABORT };

/* Line search along the direction @d from the current iterate x
   (with function value *f and gradient @g), using dcsrch(); the
   point z = x + d is the subspace minimizer. On successful return
   (LNSRCH_OK), x, f and g hold the new iterate, @t and @r hold the
   previous x and g, and @ls holds the step length and the
   directional derivatives needed for the BFGS update. If f or g is
   not finite at a trial point the step is halved. Returns
   LNSRCH_ASCENT if d is not a descent direction, LNSRCH_MAXBACK if
   the search needed too many function evaluations, or LNSRCH_ABORT
   if @fgfunc requested a stop.
*/

static int lnsrlb (int n, const double *l, const double *u,
                   const int *nbd, double *x, double *f, double *g,
                   const double *d, const double *z, double *t,
                   double *r, int iter, int boxed, int cnstnd,
                   int *nfgv, line_search *ls, LBFGSB_FG fgfunc,
                   void *data)
{
    const double ftol = .001, gtol = .9, xtol = .1;
    double a1, a2, dnorm, stpmx;
    int i, ifun = 0;
    int task = LS_START;

    ls->dtd = dot(n, d, d);
    dnorm = sqrt(ls->dtd);

    /* determine the maximum step length */
    stpmx = 1e10;
    if (cnstnd) {
        if (iter == 0) {
            stpmx = 1.0;
        } else {
            for (i=0; i<n; i++) {
                a1 = d[i];
                if (nbd[i] != 0) {
                    if (a1 < 0.0 && nbd[i] <= 2) {
                        a2 = l[i] - x[i];
                        if (a2 >= 0.0) {
                            stpmx = 0.0;
                        } else if (a1 * stpmx < a2) {
                            stpmx = a2 / a1;
                        }
                    } else if (a1 > 0.0 && nbd[i] >= 2) {
                        a2 = u[i] - x[i];
                        if (a2 <= 0.0) {
                            stpmx = 0.0;
                        } else if (a1 * stpmx > a2) {
                            stpmx = a2 / a1;
                        }
                    }
                }
            }
        }
    }

    if (iter == 0 && !boxed) {
        ls->stp = dmin(1.0 / dnorm, stpmx);
    } else {
        ls->stp = 1.0;
    }

    copy(n, x, t);
    copy(n, g, r);
    ls->fold = *f;

    while (1) {
        ls->gd = dot(n, g, d);
        if (ifun == 0) {
            ls->gdold = ls->gd;
            if (ls->gd >= 0.0) {
                /* the directional derivative is >= 0:
                   the line search is impossible */
                return LNSRCH_ASCENT;
            }
        }

        task = dcsrch(*f, ls->gd, &ls->stp, ftol, gtol, xtol, 0.0,
                      stpmx, task, ls);
        if (task == LS_CONV || task == LS_WARN) {
            return LNSRCH_OK;
        }

        /* evaluate f and g at a new trial point, unless we've
           already backtracked too many times */
        while (1) {
            if (ifun >= 20) {
                return LNSRCH_MAXBACK;
            }
            ifun++;
            (*nfgv)++;
            if (ls->stp == 1.0) {
                copy(n, z, x);
            } else {
                for (i=0; i<n; i++) {
                    x[i] = ls->stp * d[i] + t[i];
                }
            }
            /* Not in the original: rounding error can put x a few ulps
               outside of the bounds, which matters if the function is
               undefined there */
            project(n, l, u, nbd, x);
            if (fgfunc(x, f, g, data)) {
                return LNSRCH_ABORT;
            }
            if (fg_finite(n, *f, g)) {
                break;
            }
            /* Not in the original: f or g is not finite at the trial
               point (dcsrch() would then extrapolate, typically to
               the same point). Halve the step and restart the line
               search, with the maximum step set to the new one so we
               never go back toward the failed point.
            */
            stpmx = 0.5 * ls->stp;
            ls->stp = stpmx;
            dcsrch(ls->fold, ls->gdold, &ls->stp, ftol, gtol, xtol, 0.0,
                   stpmx, LS_START, ls);
            task = LS_FG;
        }
    }
}

/* Compute the generalized Cauchy point (GCP) xcp: the first local
   minimizer of the quadratic model along the projected gradient
   path x(t) = P(x - t g). On exit @iwhere records the status of
   each variable at the GCP and @c holds W'(xcp - x). The arrays
   @iorder and @t (length n) are workspace for the breakpoints; @p,
   @c, @wbp and @v are vectors of length 2m. Returns 0 on success or
   non-zero if bmv() fails.
*/

static int cauchy (int n, const double *x, const double *l,
                   const double *u, const int *nbd, const double *g,
                   int *iorder, int *iwhere, double *t, double *d,
                   double *xcp, int m, const double *wy,
                   const double *ws, const double *sy, const double *wt,
                   double theta, int col, int head, double *p, double *c,
                   double *wbp, double *v, double sbgnrm,
                   double epsmch)
{
    double f1, f2, f2_org, dt, dtm, tj, tj0, tsum, bkmin;
    double dibp, dibp2, zibp, neggi, wmc, wmp, wmw;
    double tu = 0.0, tl = 0.0;
    int i, j, ibp = 0, iter, nfree, nbreak, ibkmin, nleft, pointr;
    int col2 = 2*col;
    int bnded, xlower, xupper, allfixed = 0;
    int info;

    /* check the status of the variables, reset iwhere(i) if
       necessary; compute the Cauchy direction d and the breakpoints
       t; initialize the derivative f1 and the vector p = W'd */

    if (sbgnrm <= 0.0) {
        /* x is the GCP */
        copy(n, x, xcp);
        return 0;
    }

    bnded = 1;
    nfree = n;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.0;
    f1 = 0.0;

    for (i=0; i<col2; i++) {
        p[i] = 0.0;
    }

    /* iwhere[i] = 0 if x[i] is free and has bounds
                 = 1 if x[i] is at its lower bound
                 = 2 if x[i] is at its upper bound
                 = 3 if x[i] is always fixed (l[i] == u[i])
                 = -1 if x[i] is always free (no bounds)
                 = -3 if x[i] is free but its gradient is zero
    */

    for (i=0; i<n; i++) {
        neggi = -g[i];
        if (iwhere[i] != 3 && iwhere[i] != -1) {
            /* if x[i] is not a constant and has bounds,
               compute the difference between x[i] and its bounds */
            if (nbd[i] <= 2) {
                tl = x[i] - l[i];
            }
            if (nbd[i] >= 2) {
                tu = u[i] - x[i];
            }
            /* if a variable is close enough to a bound
               we treat it as at bound */
            xlower = nbd[i] <= 2 && tl <= 0.0;
            xupper = nbd[i] >= 2 && tu <= 0.0;
            /* reset iwhere[i] */
            iwhere[i] = 0;
            if (xlower) {
                if (neggi <= 0.0) {
                    iwhere[i] = 1;
                }
            } else if (xupper) {
                if (neggi >= 0.0) {
                    iwhere[i] = 2;
                }
            } else if (fabs(neggi) <= 0.0) {
                iwhere[i] = -3;
            }
        }
        pointr = head;
        if (iwhere[i] != 0 && iwhere[i] != -1) {
            d[i] = 0.0;
        } else {
            d[i] = neggi;
            f1 -= neggi * neggi;
            /* calculate p := p - W'e_i * (g_i) */
            for (j=0; j<col; j++) {
                p[j] += wy[i + pointr*n] * neggi;
                p[col + j] += ws[i + pointr*n] * neggi;
                pointr = (pointr + 1) % m;
            }
            if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.0) {
                /* x[i] + d[i] is bounded; compute t[i] */
                iorder[nbreak] = i;
                t[nbreak] = tl / (-neggi);
                if (nbreak == 0 || t[nbreak] < bkmin) {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
                nbreak++;
            } else if (nbd[i] >= 2 && neggi > 0.0) {
                /* x[i] + d[i] is bounded; compute t[i] */
                iorder[nbreak] = i;
                t[nbreak] = tu / neggi;
                if (nbreak == 0 || t[nbreak] < bkmin) {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
                nbreak++;
            } else {
                /* x[i] + d[i] is not bounded */
                iorder[--nfree] = i;
                if (fabs(neggi) > 0.0) {
                    bnded = 0;
                }
            }
        }
    }

    /* the indices of the nonzero components of d are now stored
       in iorder[0..nbreak-1] and iorder[nfree..n-1];
       the smallest of the nbreak breakpoints is in t[ibkmin] = bkmin */

    if (theta != 1.0) {
        /* complete the initialization of p for theta not= one */
        scal(col, theta, p + col);
    }

    /* initialize GCP xcp = x */
    copy(n, x, xcp);

    if (nbreak == 0 && nfree == n) {
        /* d is a zero vector, return with the initial xcp as GCP */
        return 0;
    }

    /* initialize c = W'(xcp - x) = 0 */
    for (j=0; j<col2; j++) {
        c[j] = 0.0;
    }

    /* initialize derivative f2 */
    f2 = -theta * f1;
    f2_org = f2;
    if (col > 0) {
        info = bmv(m, sy, wt, col, p, v);
        if (info != 0) {
            return info;
        }
        f2 -= dot(col2, v, p);
    }
    dtm = -f1 / f2;
    tsum = 0.0;

    if (nbreak > 0) {
        nleft = nbreak;
        iter = 1;
        tj = 0.0;

        /* loop over the breakpoints, looking for the GCP */
        while (1) {
            /* find the next smallest breakpoint, tj */
            tj0 = tj;
            if (iter == 1) {
                /* since we already have the smallest breakpoint we need
                   not do heapsort yet; often only one breakpoint is used
                   and the cost of heapsort is avoided */
                tj = bkmin;
                ibp = iorder[ibkmin];
            } else {
                if (iter == 2) {
                    /* replace the already used smallest breakpoint with
                       the last one, before the heapsort call */
                    if (ibkmin != nbreak - 1) {
                        t[ibkmin] = t[nbreak - 1];
                        iorder[ibkmin] = iorder[nbreak - 1];
                    }
                }
                /* update the heap structure of breakpoints (if iter = 2,
                   initialize heap) */
                hpsolb(nleft, t, iorder, iter - 2);
                tj = t[nleft - 1];
                ibp = iorder[nleft - 1];
            }
            dt = tj - tj0;

            /* if a minimizer is within this interval, locate the GCP
               and return */
            if (dtm < dt) {
                break;
            }

            /* otherwise fix one variable and reset the corresponding
               component of d to zero */
            tsum += dt;
            nleft--;
            iter++;
            dibp = d[ibp];
            d[ibp] = 0.0;
            if (dibp > 0.0) {
                zibp = u[ibp] - x[ibp];
                xcp[ibp] = u[ibp];
                iwhere[ibp] = 2;
            } else {
                zibp = l[ibp] - x[ibp];
                xcp[ibp] = l[ibp];
                iwhere[ibp] = 1;
            }
            if (nleft == 0 && nbreak == n) {
                /* all n variables are fixed, return with xcp as GCP */
                dtm = dt;
                allfixed = 1;
                break;
            }

            /* update the derivative information */
            dibp2 = dibp * dibp;

            /* update f1 and f2; temporarily set f1 and f2 for col = 0 */
            f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
            f2 -= theta * dibp2;

            if (col > 0) {
                /* update c = c + dt*p */
                axpy(col2, dt, p, c);
                /* choose wbp, the row of W corresponding to the
                   breakpoint encountered */
                pointr = head;
                for (j=0; j<col; j++) {
                    wbp[j] = wy[ibp + pointr*n];
                    wbp[col + j] = theta * ws[ibp + pointr*n];
                    pointr = (pointr + 1) % m;
                }
                /* compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)' */
                info = bmv(m, sy, wt, col, wbp, v);
                if (info != 0) {
                    return info;
                }
                wmc = dot(col2, c, v);
                wmp = dot(col2, p, v);
                wmw = dot(col2, wbp, v);
                /* update p = p - dibp*wbp */
                axpy(col2, -dibp, wbp, p);
                /* complete updating f1 and f2 while col > 0 */
                f1 += dibp * wmc;
                f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
            }

            f2 = dmax(epsmch * f2_org, f2);
            if (nleft > 0) {
                /* repeat the loop for unsearched intervals */
                dtm = -f1 / f2;
            } else {
                if (bnded) {
                    f1 = 0.0;
                    f2 = 0.0;
                    dtm = 0.0;
                } else {
                    dtm = -f1 / f2;
                }
                break;
            }
        }
    }

    if (!allfixed) {
        if (dtm <= 0.0) {
            dtm = 0.0;
        }
        tsum += dtm;
        /* move free variables (i.e., the ones w/o breakpoints) and
           the variables whose breakpoints haven't been reached */
        axpy(n, tsum, d, xcp);
    }

    /* update c = c + dtm*p = W'(x^c - x), which will be used in
       computing r = Z'(B(x^c - x) + g) */
    if (col > 0) {
        axpy(col2, dtm, p, c);
    }

    return 0;
}

/* State of the limited-memory BFGS matrix */

typedef struct {
    int col;      /* number of correction pairs currently stored */
    int head;     /* ring-buffer position of the oldest pair */
    int itail;    /* ring-buffer position of the newest pair */
    int iupdat;   /* number of updates since the last reset */
    int updatd;   /* was the matrix updated at the last iteration? */
    double theta; /* scaling factor */
} bfgs_memory;

/* discard the stored correction pairs */

static void reset_memory (bfgs_memory *bm)
{
    bm->col = 0;
    bm->head = 0;
    bm->theta = 1.0;
    bm->iupdat = 0;
    bm->updatd = 0;
}

/**
 * lbfgsb:
 * @n: number of variables.
 * @m: number of correction pairs to store (3 <= m <= 20 recommended).
 * @x: on input, the starting point; on output, the final iterate.
 * @l: array of lower bounds.
 * @u: array of upper bounds.
 * @nbd: array of bound types: for each i, 0 means x[i] is unbounded,
 * 1 means it has only a lower bound l[i], 2 means both lower and
 * upper bounds, and 3 means only an upper bound u[i].
 * @factr: iteration stops when the relative reduction in f is at
 * most factr times the machine precision.
 * @pgtol: iteration stops when the infinity norm of the projected
 * gradient is at most pgtol.
 * @f: on output, the function value at @x.
 * @g: array of length @n; on output, the gradient at @x.
 * @fgfunc: callback to compute the function value and gradient.
 * @iterfunc: optional callback (may be NULL), called at the starting
 * point and after each iteration.
 * @data: pointer passed to @fgfunc and @iterfunc.
 * @iter: location to receive the number of iterations, or NULL.
 * @nfg: location to receive the number of function evaluations,
 * or NULL.
 * @msg: buffer of at least 60 bytes to receive a message describing
 * the outcome, or NULL.
 *
 * Minimizes a function subject to simple bounds, using the L-BFGS-B
 * algorithm of Byrd, Lu, Nocedal and Zhu (version 3.0, with the
 * modifications described by Morales and Nocedal, 2011). If the
 * line search fails, @x, @f and @g hold the last successful iterate.
 * Unlike the original, this version tolerates non-finite values of
 * f or g at trial points (e.g. where the function is undefined),
 * backtracking toward the current iterate.
 *
 * Returns: LBFGSB_CONVERGED, LBFGSB_ABNORMAL (the line search failed
 * to make progress), LBFGSB_STOPPED (a callback requested a stop),
 * LBFGSB_BADINPUT, LBFGSB_NOMEM, or LBFGSB_NONFINITE (f or g is
 * not finite at the starting point).
 */

int lbfgsb (int n, int m, double *x, const double *l, const double *u,
            const int *nbd, double factr, double pgtol, double *f,
            double *g, LBFGSB_FG fgfunc, LBFGSB_ITER iterfunc,
            void *data, int *iter, int *nfg, char *msg)
{
    /* the original code's approximation to machine epsilon */
    const double epsmch = 2.22e-16;
    double *wspace = NULL;
    int *ispace = NULL;
    double *ws, *wy, *sy, *ss, *wt, *wn, *snd;
    double *z, *r, *d, *t, *xp, *wa;
    int *index, *iwhere, *indx2;
    bfgs_memory bm = {0};
    line_search ls = {0};
    double tol, sbgnrm, ddum, rr, dr;
    int prjctd, cnstnd, boxed, wrk;
    int nfree = n, nenter = 0, ileave = 0;
    int niter = 0, nfgv = 0;
    int i, info, ret;
    char task[60] = "";

    errclb(n, m, factr, l, u, nbd, task);
    if (task[0] != '\0') {
        ret = LBFGSB_BADINPUT;
        goto finish;
    }

    wspace = calloc((size_t) (2*m + 5)*n + 11*m*m + 8*m, sizeof *wspace);
    ispace = calloc(3*n, sizeof *ispace);
    if (wspace == NULL || ispace == NULL) {
        strcpy(task, "ERROR: OUT OF MEMORY");
        ret = LBFGSB_NOMEM;
        goto finish;
    }

    ws = wspace;           /* n x m: S, the corrections in x */
    wy = ws + m*n;         /* n x m: Y, the corrections in g */
    sy = wy + m*n;         /* m x m: S'Y */
    ss = sy + m*m;         /* m x m: S'S */
    wt = ss + m*m;         /* m x m: Cholesky factor of theta*S'S + LD^(-1)L' */
    wn = wt + m*m;         /* 2m x 2m: factorization of the matrix K */
    snd = wn + 4*m*m;      /* 2m x 2m: inner products for K */
    z = snd + 4*m*m;       /* the GCP, then the subspace minimizer */
    r = z + n;             /* reduced gradient; previous g */
    d = r + n;             /* search direction */
    t = d + n;             /* breakpoints; previous x */
    xp = t + n;            /* workspace for subsm() */
    wa = xp + n;           /* 8m workspace */

    index = ispace;        /* free and active variables */
    iwhere = index + n;    /* status of the variables */
    indx2 = iwhere + n;    /* entering and leaving variables */

    reset_memory(&bm);
    tol = factr * epsmch;

    /* project x onto the feasible set and initialize iwhere */
    active(n, l, u, nbd, x, iwhere, &prjctd, &cnstnd, &boxed);

    /* evaluate f and g at the starting point */
    if (fgfunc(x, f, g, data)) {
        ret = LBFGSB_STOPPED;
        goto stopped;
    }
    nfgv = 1;
    if (!fg_finite(n, *f, g)) {
        strcpy(task, "ERROR: F OR G NOT FINITE AT STARTING POINT");
        ret = LBFGSB_NONFINITE;
        goto finish;
    }
    if (iterfunc != NULL && iterfunc(0, nfgv, *f, x, g, 0.0, data)) {
        ret = LBFGSB_STOPPED;
        goto stopped;
    }

    /* is the initial x stationary? */
    sbgnrm = projgr(n, l, u, nbd, x, g);
    if (sbgnrm <= pgtol) {
        strcpy(task, "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL");
        ret = LBFGSB_CONVERGED;
        goto finish;
    }

    while (1) {
        info = 0;

        if (!cnstnd && bm.col > 0) {
            /* unconstrained and not the first iteration: skip the
               search for the generalized Cauchy point */
            copy(n, x, z);
            wrk = bm.updatd;
        } else {
            /* compute the generalized Cauchy point (GCP) */
            info = cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m,
                          wy, ws, sy, wt, bm.theta, bm.col, bm.head,
                          wa, wa + 2*m, wa + 4*m, wa + 6*m, sbgnrm,
                          epsmch);
            if (info != 0) {
                /* singular triangular system detected: refresh the
                   memory and restart the iteration */
                reset_memory(&bm);
                continue;
            }
            /* count the entering and leaving variables, and find the
               index set of free and active variables at the GCP */
            wrk = freev(n, &nfree, index, &nenter, &ileave, indx2,
                        iwhere, bm.updatd, cnstnd, niter);
        }

        /* subspace minimization, if there are free variables at the
           GCP and the BFGS matrix is not the identity */
        if (nfree != 0 && bm.col != 0) {
            if (wrk) {
                /* form the LEL' factorization of the indefinite
                   matrix K */
                info = formk(n, nfree, index, nenter, ileave, indx2,
                             bm.iupdat, bm.updatd, wn, snd, m, ws, wy, sy,
                             bm.theta, bm.col, bm.head);
            }
            if (info == 0) {
                /* compute r = -Z'B(xcp - xk) - Z'g */
                info = cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index,
                              bm.theta, bm.col, bm.head, nfree, cnstnd);
            }
            if (info == 0) {
                /* compute the subspace minimizer in z */
                info = subsm(n, m, nfree, index, l, u, nbd, z, r, xp, ws,
                             wy, bm.theta, x, g, bm.col, bm.head, wa, wn);
            }
            if (info != 0) {
                /* nonpositive definiteness in a Cholesky factorization
                   or a singular triangular system: refresh the memory
                   and restart the iteration */
                reset_memory(&bm);
                continue;
            }
        }

        /* line search along the direction d = z - x */
        for (i=0; i<n; i++) {
            d[i] = z[i] - x[i];
        }
        info = lnsrlb(n, l, u, nbd, x, f, g, d, z, t, r, niter, boxed,
                      cnstnd, &nfgv, &ls, fgfunc, data);
        if (info == LNSRCH_ABORT) {
            ret = LBFGSB_STOPPED;
            goto stopped;
        } else if (info != LNSRCH_OK) {
            /* restore the previous iterate */
            copy(n, t, x);
            copy(n, r, g);
            *f = ls.fold;
            if (bm.col == 0) {
                /* we were already using the steepest descent
                   direction: give up */
                strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
                niter++;
                ret = LBFGSB_ABNORMAL;
                goto finish;
            } else {
                /* refresh the memory and restart the iteration */
                reset_memory(&bm);
                continue;
            }
        }

        /* we have a new iterate */
        niter++;
        sbgnrm = projgr(n, l, u, nbd, x, g);
        if (iterfunc != NULL &&
            iterfunc(niter, nfgv, *f, x, g, ls.stp, data)) {
            ret = LBFGSB_STOPPED;
            goto stopped;
        }

        /* test for termination */
        if (sbgnrm <= pgtol) {
            strcpy(task, "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL");
            ret = LBFGSB_CONVERGED;
            goto finish;
        }
        ddum = dmax(dmax(fabs(ls.fold), fabs(*f)), 1.0);
        if (ls.fold - *f <= tol * ddum) {
            strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH");
            ret = LBFGSB_CONVERGED;
            goto finish;
        }

        /* compute d = newx - oldx (as stp * d), r = newg - oldg,
           rr = y'y and dr = y's */
        for (i=0; i<n; i++) {
            r[i] = g[i] - r[i];
        }
        rr = dot(n, r, r);
        if (ls.stp == 1.0) {
            dr = ls.gd - ls.gdold;
            ddum = -ls.gdold;
        } else {
            dr = (ls.gd - ls.gdold) * ls.stp;
            scal(n, ls.stp, d);
            ddum = -ls.gdold * ls.stp;
        }

        if (dr <= epsmch * ddum) {
            /* skip the L-BFGS update */
            bm.updatd = 0;
            continue;
        }

        /* update the L-BFGS matrix */
        bm.updatd = 1;
        bm.iupdat++;
        matupd(n, m, ws, wy, sy, ss, d, r, &bm.itail, bm.iupdat, &bm.col,
               &bm.head, &bm.theta, rr, dr, ls.stp, ls.dtd);

        /* form the upper half of T = theta*SS + L*D^(-1)*L' and
           Cholesky-factorize it */
        if (formt(m, wt, sy, ss, bm.col, bm.theta) != 0) {
            reset_memory(&bm);
        }
    }

 stopped:
    strcpy(task, "STOP: REQUESTED BY CALLBACK");

 finish:
    if (msg != NULL) {
        strcpy(msg, task);
    }
    if (iter != NULL) {
        *iter = niter;
    }
    if (nfg != NULL) {
        *nfg = nfgv;
    }
    free(wspace);
    free(ispace);

    return ret;
}
