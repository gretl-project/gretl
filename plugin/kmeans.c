/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2017 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/* Code for the K means problem, based on John Burkardt's C rendition of
   Applied Statistics Algorithm 136, originally coded in FORTRAN by
   J. A. Hartigan and M. A. Wong, Journal of the Royal Statistical
   Society Series C: Applied Statistics, Volume 28, Issue 1, March 1979,
   Pages 100–108, https://doi.org/10.2307/2346830 .  Adapted for gretl
   by Jack Lucchetti, May 2026.
*/

#include "libgretl.h"
#include "libset.h"

static int init_centers (int *nc, gretl_matrix *c, int *ic1,
			 const gretl_matrix *a)
{
    int i, j, l;
    int k = c->rows;
    int m = a->rows;
    int n = a->cols;
    int err = 0;
    double temp;

    for (l=0; l<k; l++) {
	nc[l] = 0;
    }

    gretl_matrix_zero(c);

    for (i=0; i<m; i++) {
	l = ic1[i] - 1;
	nc[l]++;
	for (j=0; j<n; j++) {
	    temp = gretl_matrix_get(c, l, j) + gretl_matrix_get(a, i, j);
	    gretl_matrix_set(c, l, j, temp);
	}
    }

    /* Check to see if there is any empty cluster at this stage */
    for (l = 0; l < k; l++)  {
	if (nc[l] == 0)  {
	    err = 1;
	    break;
	}
    }

    return err;
}


static double compute_ith_distance (int i, const gretl_matrix *a,
				    const gretl_matrix *c, int l1)
{
    double dif, dist = 0.0;
    int j;

    for (j=0; j<a->cols; j++) {
	dif = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, l1-1, j);
	dist += dif * dif;
    }

    return dist;
}

static void update_an (gretl_matrix *an, int l1, int l2, double al1, double al2)
{
    double tmp = (al1 > 2.0) ? (al1 - 1.0) / (al1 - 2.0) : 1.0e100;

    gretl_matrix_set(an, l1-1, 1, (al1 - 1.0) / al1);
    gretl_matrix_set(an, l1-1, 0, tmp);
    gretl_matrix_set(an, l2-1, 0, (al2 + 1.0) / al2);
    gretl_matrix_set(an, l2-1, 1, (al2 + 1.0) / (al2 + 2.0));
}

/******************************************************************************/

static void optra (const gretl_matrix *a, gretl_matrix *c, int *ic1,
		   int *ic2, int *nc, gretl_matrix *an, int ncp[],
		   gretl_vector *d, int itran[], int live[], int *indx)

/******************************************************************************/
/*
  Purpose:

    OPTRA carries out the optimal transfer stage.

  Discussion:

    This is the optimal transfer stage.

    Each point is re-allocated, if necessary, to the cluster that
    will induce a maximum reduction in the within-cluster sum of
    squares.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Author:

    Original FORTRAN77 version by John Hartigan, Manchek Wong.
    Original C version by John Burkardt.
    Adapted for libgretl by Jack Lucchetti.

  Reference:

    John Hartigan, Manchek Wong,
    Algorithm AS 136:
    A K-Means Clustering Algorithm,
    Applied Statistics,
    Volume 28, Number 1, 1979, pages 100-108.

  Parameters:

    Input, double A(M,N), the points.
    Input/output, double C(K,N), the cluster centers.
    Input/output, int IC1(M), the cluster to which each point is assigned.
    Input/output, int IC2(M), used to store the cluster which each point
    	is most likely to be transferred to at each step.

    Input/output, int NC(K), the number of points in each cluster.
    Input/output, double AN(K,2).
    Input/output, int NCP(K).
    Input/output, double D(M).
    Input/output, int ITRAN(K).
    Input/output, int LIVE(K).
    Input/output, int *INDX, the number of steps since a transfer took place.
*/
{
    int m = a->rows;
    int n = a->cols;
    int k = c->rows;
    double aij;
    double al1;
    double al2;
    double alt;
    double alw;
    double da;
    double dc;
    double de;
    int i, j, l;
    int l1, l2, ll;
    double r2;
    double rr;
    double tmp;

    /* If cluster L is updated in the last quick-transfer stage, it
       belongs to the live set throughout this stage.  Otherwise, at
       each step, it is not in the live set if it has not been updated
       in the last M optimal transfer steps.
    */
    for (l=0; l<k; l++) {
	if (itran[l] == 1) {
	    live[l] = m + 1;
	}
    }

    for (i=0; i<m; i++) {
	*indx = *indx + 1;
	l1 = ic1[i];
	l2 = ic2[i];
	ll = l2;

	/* If point I is the only member of cluster L1, no transfer */

	if (nc[l1-1] > 1) {
	    /* If L1 has not yet been updated in this stage, no need to
	       re-compute D(I).
	    */
	    if (ncp[l1-1] != 0) {
		de = compute_ith_distance(i, a, c, l1);
		gretl_vector_set(d, i, de * gretl_matrix_get(an, l1-1, 0));
	    }

	    /* Find the cluster with minimum R2 */
	    da = compute_ith_distance(i, a, c, l2);
	    r2 = da * gretl_matrix_get(an, l2-1, 1);

	    for (l=1; l<=k; l++) {
		/* If LIVE(L1) <= I, then L1 is not in the live set.  If
		   this is true, we only need to consider clusters that
		   are in the live set for possible transfer of point I.
		   Otherwise, we need to consider all possible clusters.
		*/
		if ((i < live[l1-1]-1 || i < live[l2-1]-1) && l != l1 && l != ll) {
		    rr = r2 / gretl_matrix_get(an, l-1, 1);
		    dc = compute_ith_distance(i, a, c, l);
		    if (dc < rr)  {
			r2 = dc * gretl_matrix_get(an, l-1, 1);
			l2 = l;
		    }
		}
	    }

	    /* If no transfer is necessary, L2 is the new IC2(I) */
	    if (gretl_vector_get(d, i) <= r2) {
		ic2[i] = l2;
	    } else {
		/* Update cluster centers, LIVE, NCP and AN for clusters
		   L1 and L2, and update IC1(I) and IC2(I).
		*/
		*indx = 0;
		live[l1-1] = m + i;
		live[l2-1] = m + i;
		ncp[l1-1] = i + 1;
		ncp[l2-1] = i + 1;
		al1 = (double) ( nc[l1-1] );
		alw = al1 - 1.0;
		al2 = (double) ( nc[l2-1] );
		alt = al2 + 1.0;
		for (j=0; j<n; j++) {
		    aij = gretl_matrix_get(a, i, j);
		    tmp = (gretl_matrix_get(c, l1-1, j) * al1 - aij ) / alw;
		    gretl_matrix_set(c, l1-1, j, tmp);
		    tmp = (gretl_matrix_get(c, l2-1, j) * al2 + aij ) / alt;
		    gretl_matrix_set(c, l2-1, j, tmp);
		}

		nc[l1-1] = nc[l1-1] - 1;
		nc[l2-1] = nc[l2-1] + 1;
		update_an(an, l1, l2, al1, al2);
		ic1[i] = l2;
		ic2[i] = l1;
	    }
	}

	if (*indx == m) {
	    return;
	}
    }

    /* ITRAN(L) = 0 before entering QTRAN. Also, LIVE(L) has to be
       decreased by M before re-entering OPTRA.
    */
    for (l=1; l<=k; l++) {
	itran[l-1] = 0;
	live[l-1] = live[l-1] - m;
    }
}

/******************************************************************************/

static void qtran (const gretl_matrix *a, gretl_matrix *c, int *ic1,
		   int *ic2, int *nc, gretl_matrix *an, int ncp[],
		   gretl_vector *d, int itran[], int *indx)

/******************************************************************************/
/*
  Purpose:

    QTRAN carries out the quick transfer stage.

  Discussion:

    This is the quick transfer stage.

    IC1(I) is the cluster which point I belongs to.
    IC2(I) is the cluster which point I is most likely to be
    transferred to.

    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
    reduce within-cluster sum of squares.  The cluster centers are
    updated after each step.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Author:

    Original FORTRAN77 version by John Hartigan, Manchek Wong.
    Original C version by John Burkardt.
    Adapted for libgretl by Jack Lucchetti.

  Reference:

    John Hartigan, Manchek Wong,
    Algorithm AS 136:
    A K-Means Clustering Algorithm,
    Applied Statistics,
    Volume 28, Number 1, 1979, pages 100-108.

  Parameters:

    Input, double A(M,N), the points.
    Input/output, double C(K,N), the cluster centers.

    Input/output, int IC1(M), the cluster to which each
    point is assigned.

    Input/output, int IC2(M), used to store the cluster
    which each point is most likely to be transferred to at each step.

    Input/output, int NC(K), the number of points in
    each cluster.

    Input/output, double AN(K, 2).

    Input/output, int NCP(K).

    Input/output, double D(M).

    Input/output, int ITRAN(K).

    Input/output, int INDX, counts the number of steps
    since the last transfer.
*/
{
    int m = a->rows;
    int n = a->cols;
    double aij;
    double al1;
    double al2;
    double alt;
    double alw;
    double di;
    int i;
    int icoun;
    int istep;
    int done = 0;
    int j;
    int l1;
    int l2;
    double r2;
    double tmp;

    /* In the optimal transfer stage, NCP(L) indicates the step at which
       cluster L is last updated.  In the quick transfer stage, NCP(L)
       is equal to the step at which cluster L is last updated plus M.
    */
    icoun = 0;
    istep = 0;

    while (!done) {
	for (i=0; i<m; i++) {
	    icoun++;
	    istep++;
	    l1 = ic1[i];
	    l2 = ic2[i];

	    /* If point I is the only member of cluster L1, no transfer */
	    if (nc[l1-1] > 1) {
		/* If NCP(L1) < ISTEP, no need to re-compute distance
		   from point I to cluster L1.  Note that if cluster L1
		   is last updated exactly M steps ago, we still need to
		   compute the distance from point I to cluster L1.
		*/
		if (istep <= ncp[l1-1]) {
		    di = compute_ith_distance(i, a, c, l1);
		    gretl_vector_set(d, i, di * gretl_matrix_get(an, l1-1, 0));
		}
		/* If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will
		   be no transfer of point I at this step.
		*/
		if (istep < ncp[l1-1] || istep < ncp[l2-1]) {
		    r2 = gretl_vector_get(d, i) / gretl_matrix_get(an, l2-1, 1);
		    di = compute_ith_distance(i, a, c, l2);

		    /* Update cluster centers, NCP, NC, ITRAN, AN1 and
		      AN2 for clusters L1 and L2.  Also update IC1(I)
		      and IC2(I).  Note that if any updating occurs in
		      this stage, INDX is set back to 0.
		    */
		    if (di < r2) {
			icoun = 0;
			*indx = 0;
			itran[l1-1] = 1;
			itran[l2-1] = 1;
			ncp[l1-1] = istep + m;
			ncp[l2-1] = istep + m;
			al1 = (double) (nc[l1-1]);
			alw = al1 - 1.0;
			al2 = (double) (nc[l2-1]);
			alt = al2 + 1.0;
			for (j=0; j<n; j++) {
			    aij = gretl_matrix_get(a, i, j);
			    tmp = (gretl_matrix_get(c, l1-1, j) * al1 - aij) / alw;
			    gretl_matrix_set(c, l1-1, j, tmp);
			    tmp = (gretl_matrix_get(c, l2-1, j) * al2 + aij) / alt;
			    gretl_matrix_set(c, l2-1, j, tmp);
			}
			nc[l1-1] = nc[l1-1] - 1;
			nc[l2-1] = nc[l2-1] + 1;
			update_an(an, l1, l2, al1, al2);
			ic1[i] = l2;
			ic2[i] = l1;
		    }
		}
	    }

	    /* If no re-allocation took place in the last M steps,
	       we're finished.
	    */
	    done = icoun == m;
	    if (done) {
		break;
	    }
	}
    }
}

/******************************************************************************/

gretl_matrix *kmeans (const gretl_matrix *a, int k,
		      gretl_matrix **clustinfo,
		      int *err)

/******************************************************************************/
/*
  Purpose:

    kmeans carries out the K-means algorithm.

  Discussion:

    This routine attempts to divide M points in N-dimensional space into
    K clusters so that the within cluster sum of squares is minimized.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Author:

    Original FORTRAN77 version by John Hartigan, Manchek Wong.
    Original C version by John Burkardt.
    Adapted for libgretl by Jack Lucchetti.

  Reference:

    John Hartigan, Manchek Wong,
    Algorithm AS 136:
    A K-Means Clustering Algorithm,
    Applied Statistics,
    Volume 28, Number 1, 1979, pages 100-108.

  Parameters:

    Input, double A(M,N), the points.
    Input/output, double C(K,N), the cluster centers.
    Output, int IC1(M), the cluster to which each point is assigned.
    Output, double CLUSTINFO(K,2): the number of points in each cluster (col. 1),
    the within-cluster sum of squares of each cluster (col. 2).

  Returs: error code.
    0, no error was detected.
    1, at least one cluster is empty after the initial assignment.  A better
       set of initial cluster centers is needed.
    E_NOCONV, the allowed maximum number off iterations was exceeded.
    E_NONCONF, K is less than or equal to 1, or greater than or equal to M.
*/
{
    gretl_vector *d;
    gretl_vector *clustid;
    gretl_matrix *an;
    gretl_matrix *c;
    double aa, da, db, dc;
    double dt[2];
    int i;
    int *ic1;
    int *ic2;
    int ij;
    int il;
    int indx;
    int *itran;
    int j;
    int l;
    int *live;
    int *nc;
    int *ncp;
    double temp;
    int maxiter = 128;
    double HUGE = libset_get_double(CONV_HUGE);
    int m = a->rows;
    int n = a->cols;

    if (k <= 1 || m <= k) {
	*err = E_NONCONF;
	return NULL;
    }

    ic1 = malloc(m * sizeof *ic2);
    ic2 = malloc(m * sizeof *ic2);
    nc = malloc(k * sizeof *nc);
    ncp = malloc(k * sizeof *ncp);
    itran = malloc(k * sizeof *itran);
    live = malloc(k * sizeof *live);

    d = gretl_vector_alloc(m);
    an = gretl_matrix_alloc(k, 2);
    c = gretl_matrix_alloc(k, n);
    clustid = gretl_column_vector_alloc(m);

    if (ic1 == NULL || ic2 == NULL || nc == NULL || ncp == NULL ||
	itran == NULL || live == NULL || d == NULL || an == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    /* Initialize the cluster centers. Here, we arbitrarily make the
       first K data points cluster centers.
    */
    for (i=0; i<k; i++)  {
	for (j=0; j<n; j++)  {
	    gretl_matrix_set(c, i, j, gretl_matrix_get(a, i, j));
	}
    }

    /* For each point I, find its two closest centers, IC1(I) and
       IC2(I).  Assign the point to IC1(I).
    */
    for (i=0; i<m; i++) {
	ic1[i] = 1;
	ic2[i] = 2;

	for (il=0; il<2; il++) {
	    dt[il] = 0.0;
	    for (j=0; j<n; j++) {
		da = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, il, j);
		dt[il] = dt[il] + da * da;
	    }
	}

	if (dt[1] < dt[0]) {
	    ic1[i] = 2;
	    ic2[i] = 1;
	    temp = dt[0];
	    dt[0] = dt[1];
	    dt[1] = temp;
	}

	for (l=2; l<k; l++) {
	    db = 0.0;
	    for (j=0; j<n; j++) {
		dc = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, l, j);
		db = db + dc * dc;
	    }

	    if (db < dt[1]) {
		if (dt[0] <= db) {
		    dt[1] = db;
		    ic2[i] = l+1;
		} else {
		    dt[1] = dt[0];
		    ic2[i] = ic1[i];
		    dt[0] = db;
		    ic1[i] = l+1;
		}
	    }
	}
    }

    /* Update cluster centers to be the average of points contained
       within them.
    */
    *err = init_centers(nc, c, ic1, a);
    if (*err) {
	/* there's at least one empty cluster */
	goto bailout;
    }

    for (l=0; l<k; l++)  {
	/* compute centroids as averages */
	aa = (double) (nc[l]);
	for (j=0; j<n; j++) {
	    temp = gretl_matrix_get(c, l, j) / aa;
	    gretl_matrix_set(c, l, j, temp);
	}
	/* Initialize AN, ITRAN and NCP.

	   AN(L,1) = NC(L) / (NC(L) - 1)
	   AN(L,2) = NC(L) / (NC(L) + 1)
	   ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
	   = 0 otherwise

	   In the optimal-transfer stage, NCP(L) stores the step at
	   which cluster L is last updated.

	   In the quick-transfer stage, NCP(L) stores the step at which
	   cluster L is last updated plus M.
	*/
	temp = (1.0 < aa) ? aa / (aa - 1.0) : HUGE;
	gretl_matrix_set(an, l, 0, temp);
	gretl_matrix_set(an, l, 1, aa / (aa + 1.0));
	itran[l] = 1;
	ncp[l] = -1;
    }

    indx = 0;
    *err = E_NOCONV;

    for (ij=0; ij<maxiter; ij++)  {
	/* In this stage, there is only one pass through the data.  Each
	   point is re-allocated, if necessary, to the cluster that will
	   induce the maximum reduction in within-cluster sum of
	   squares.
	*/
	optra(a, c, ic1, ic2, nc, an, ncp, d, itran, live, &indx);

	/* Stop if no transfer took place in the last M optimal transfer
	   steps.
	*/
	if (indx == m) {
	    *err = 0;
	    break;
	}

	/* Each point is tested in turn to see if it should be
	   re-allocated to the cluster to which it is most likely to be
	   transferred, IC2(I), from its present cluster, IC1(I).  Loop
	   through the data until no further change is to take place.
	*/
	qtran(a, c, ic1, ic2, nc, an, ncp, d, itran, &indx);

	/* If there are only two clusters, there is no need to re-enter
	   the optimal transfer stage.
	*/
	if (k == 2) {
	    *err = 0;
	    break;
	}

	/* NCP has to be set to 0 before entering OPTRA */
	for (l=0; l<k; l++) {
	    ncp[l] = 0;
	}
    }

    if (*err) {
	goto bailout;
    }

    if (clustinfo != NULL) {
	/* Compute the within-cluster sum of squares for each cluster */
	gretl_matrix *cinfo = gretl_zero_matrix_new(k, 2);

	gretl_matrix_zero(c);
	for (i=0; i<m; i++) {
	    l = ic1[i];
	    for (j=0; j<n; j++)  {
		temp = gretl_matrix_get(c, l-1, j) + gretl_matrix_get(a, i, j);
		gretl_matrix_set(c, l-1, j, temp);
	    }
	}
	for (j=0; j<n; j++) {
	    for (l=0; l<k; l++) {
		temp = gretl_matrix_get(c, l, j) / (double) (nc[l]);
		gretl_matrix_set(c, l, j, temp);
	    }
	    for (i=0; i<m; i++) {
		l = ic1[i];
		da = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, l-1, j);
		temp = gretl_matrix_get(cinfo, l-1, 1) + da * da;
		gretl_matrix_set(cinfo, l-1, 1, temp);
	    }
	}
	for (i=0; i<k; i++) {
	    gretl_matrix_set(cinfo, i, 0, nc[i]);
	}

	if (*clustinfo != NULL) {
	    /* or we could insist on a correctly sized matrix on input */
	    gretl_matrix_free(*clustinfo);
	}
	*clustinfo = cinfo;
    }

    /* fill the output vector */
    for (i=0; i<m; i++) {
	gretl_vector_set(clustid, i, ic1[i]);
    }

 bailout:

    gretl_matrix_free(an);
    gretl_matrix_free(c);
    gretl_vector_free(d);

    free(ic1);
    free(ic2);
    free(nc);
    free(ncp);
    free(itran);
    free(live);

    return clustid;
}
