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

/* Code for the K-means problem, based on John Burkardt's C rendition of
   Applied Statistics Algorithm 136 under the GNU GPL; see
   https://people.math.sc.edu/Burkardt/cpp_src/asa136/asa136.html

   This algorithm was originally coded in FORTRAN by J. A. Hartigan and
   M. A. Wong. See Journal of the Royal Statistical Society Series C:
   Applied Statistics, Volume 28, Issue 1, March 1979, Pages 100–108.
   See https://doi.org/10.2307/2346830

   Adapted for gretl by Jack Lucchetti, May 2026.
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
    for (l=0; l<k; l++)  {
	if (nc[l] == 0)  {
	    err = 1;
	    break;
	}
    }

    return err;
}

/* Get the vector of per-cluster SSTs */

static gretl_matrix *compute_sst_full (const gretl_matrix *a,
				       const gretl_matrix *c,
				       const int *ic1)
{
    gretl_vector *sstvec;
    double d;
    int i, j, l;

    sstvec = gretl_zero_matrix_new(c->rows, 1);

    for (i=0; i<a->rows; i++) {
	l = ic1[i] - 1;
	for (j=0; j<a->cols; j++) {
	    d = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, l, j);
	    sstvec->val[l] += d * d;
	}
    }

    return sstvec;
}

#if 0 /* not yet: for use with iteration */

/* Just get the sum of per-cluster SSTs */

static double compute_sst (const gretl_matrix *a,
			   const gretl_matrix *c,
			   const int *ic1)
{
    double sst = 0;
    double d;
    int i, j, l;

    for (i=0; i<a->rows; i++) {
	l = ic1[i] - 1;
	for (j=0; j<a->cols; j++) {
	    d = gretl_matrix_get(a, i, j) - gretl_matrix_get(c, l, j);
	    sst += d * d;
	}
    }

    return sst;
}

static void use_k_random_candidates (const gretl_matrix *a,
				     gretl_matrix *c)
{
    int k = c->rows;
    int n = c->cols;
    int m = a->rows;
    gretl_vector *v = gretl_vector_alloc(k);
    int i, j, r;

    /* fill @v with @k draws from 1,2,...,@m without replacement */
    fill_permutation_vector(v, m);

    for (i=0; i<k; i++)  {
	r = v->val[i];
	for (j=0; j<n; j++)  {
	    gretl_matrix_set(c, r, j, gretl_matrix_get(a, r, j));
	}
    }

    gretl_matrix_free(v);
}

#endif /* not yet */

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

/* optra() carries out the optimal transfer stage: each point is
   reallocated, if necessary, to the cluster that will induce the
   maximum reduction in the within-cluster sum of squares. This
   involves a single pass through the data.

   @a (m x n): the data points
   @c (k x n): the cluster centers
   @ic1 (m): the cluster to which each point is assigned
   @ic2 (m): the cluster to which each point is most likely
     to be transferred to at each step
   @nc (k): the number of points in each cluster
   @an (k x 2)
   @ncp (k)
   @d (m)
   @itran (k)
   @live (k)
   @indx: the number of steps since a transfer took place

*/

static void optra (const gretl_matrix *a, gretl_matrix *c, int *ic1,
		   int *ic2, int *nc, gretl_matrix *an, int ncp[],
		   gretl_vector *d, int itran[], int live[], int *indx)
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

    /* If cluster l is updated in the last quick-transfer stage, it
       belongs to the live set throughout this stage.  Otherwise, at
       each step, it is not in the live set if it has not been updated
       in the last m optimal transfer steps.
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

	/* If point i is the only member of cluster l1, no transfer */

	if (nc[l1-1] > 1) {
	    /* If l1 has not yet been updated in this stage, no need to
	       re-compute d[i].
	    */
	    if (ncp[l1-1] != 0) {
		de = compute_ith_distance(i, a, c, l1);
		gretl_vector_set(d, i, de * gretl_matrix_get(an, l1-1, 0));
	    }

	    /* Find the cluster with minimum r2 */
	    da = compute_ith_distance(i, a, c, l2);
	    r2 = da * gretl_matrix_get(an, l2-1, 1);

	    for (l=1; l<=k; l++) {
		/* If live[l1] <= i, then l1 is not in the live set.  If
		   this is true, we only need to consider clusters that
		   are in the live set for possible transfer of point i.
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

	    /* If no transfer is necessary, l2 is the new ic2[i] */
	    if (gretl_vector_get(d, i) <= r2) {
		ic2[i] = l2;
	    } else {
		/* Update cluster centers, line, ncp and an for clusters
		   l1 and l2, and update ic1[i] and ic2[i].
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

    /* itran(l) = 0 before entering qtran(). Also, live(l) has to be
       decreased by m before re-entering optra().
    */
    for (l=1; l<=k; l++) {
	itran[l-1] = 0;
	live[l-1] = live[l-1] - m;
    }
}

/* qtran() carries out the quick transfer stage.

   @ic1 and @ic2 record the cluster to which point i belongs and the
   cluster to which it is most likely to be transferred, respectively.

   For each point i, ic1[i] and ic2[i] are switched, if necessary, to
   reduce the within-cluster sum of squares.  The cluster centers are
   updated after each step. This involves as many passes through the
   data as are required.

   @a (m x n): the data points
   @c (k x n): the cluster centers
   @ic1 (m): the cluster to which each point is assigned
   @ic2 (m): the cluster to which each point is most likely
     to be transferred to at each step
   @nc (k): the number of points in each cluster
   @an (k x 2)
   @ncp (k)
   @d (m)
   @itran (k)
   @indx: the number of steps since a transfer took place
*/

static void qtran (const gretl_matrix *a, gretl_matrix *c, int *ic1,
		   int *ic2, int *nc, gretl_matrix *an, int ncp[],
		   gretl_vector *d, int itran[], int *indx)
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

    /* In the optimal transfer stage, ncp[l] indicates the step at which
       cluster l is last updated.  In the quick transfer stage, ncp[l]
       is equal to the step at which cluster l is last updated plus m.
    */
    icoun = 0;
    istep = 0;

    while (!done) {
	for (i=0; i<m; i++) {
	    icoun++;
	    istep++;
	    l1 = ic1[i];
	    l2 = ic2[i];

	    /* If point i is the only member of cluster l1, no transfer */
	    if (nc[l1-1] > 1) {
		/* If ncp[l1] < istep, no need to recompute distance
		   from point i to cluster l1.  Note that if cluster l1
		   is last updated exactly m steps ago, we still need to
		   compute the distance from point i to cluster l1.
		*/
		if (istep <= ncp[l1-1]) {
		    di = compute_ith_distance(i, a, c, l1);
		    gretl_vector_set(d, i, di * gretl_matrix_get(an, l1-1, 0));
		}
		/* If ncp[l1] <= istep and ncp[l2] <= istep, there will
		   be no transfer of point i at this step.
		*/
		if (istep < ncp[l1-1] || istep < ncp[l2-1]) {
		    r2 = gretl_vector_get(d, i) / gretl_matrix_get(an, l2-1, 1);
		    di = compute_ith_distance(i, a, c, l2);
		    /* Update cluster centers, ncp, nc, itran, and an
		      for clusters l1 and l2.  Also update ic1[i] and
		      ic2[i].  If any updating occurs in this stage,
		      *indx is set back to 0.
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
	    /* If no re-allocation took place in the last m steps,
	       we're finished.
	    */
	    done = icoun == m;
	    if (done) {
		break;
	    }
	}
    }
}

static void find_nearest_neighbors (const gretl_matrix *a,
				    const gretl_matrix *c,
				    int *ic1, int *ic2,
				    double *dt)
{
    double da, db, dc;
    double temp;
    int i, j, l, il;

    for (i=0; i<a->rows; i++) {
	ic1[i] = 1;
	ic2[i] = 2;
	for (il=0; il<2; il++) {
	    dt[il] = 0.0;
	    for (j=0; j<a->cols; j++) {
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
	for (l=2; l<c->rows; l++) {
	    db = 0.0;
	    for (j=0; j<a->cols; j++) {
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
}

gretl_matrix *get_clustinfo_target (gretl_matrix **clustinfo,
				    int k, int n, int *reuse)
{
    gretl_matrix *cinfo = NULL;

    if (*clustinfo != NULL) {
	cinfo = *clustinfo;
	if (cinfo->rows == k && cinfo->cols == n+2) {
	    *reuse = 1;
	} else {
	    gretl_matrix_free(*clustinfo);
	    cinfo = NULL;
	}
    }
    if (cinfo == NULL) {
	cinfo = gretl_zero_matrix_new(k, n+2);
    }

    return cinfo;
}

static void add_clustinfo_colnames (const gretl_matrix *a,
				    gretl_matrix *cinfo)
{
    int n = cinfo->cols;
    char **S = strings_array_new(n);

    if (S != NULL) {
	const char **Sa = gretl_matrix_get_colnames(a);
	char tmp[12];
	int i;

	for (i=0; i<n; i++) {
	    if (i == 0) {
		S[i] = gretl_strdup("nc");
	    } else if (i == n-1) {
		S[i] = gretl_strdup("SST");
	    } else if (Sa != NULL) {
		S[i] = gretl_strdup(Sa[i-1]);
	    } else {
		sprintf(tmp, "a%d", i);
		S[i] = gretl_strdup(tmp);
	    }
	}
	gretl_matrix_set_colnames(cinfo, S);
    }
}

int init(const gretl_matrix *a, gretl_matrix *c,
	 int *ic1, int *ic2, gretl_matrix *an,
	 int *nc, int *itran, int *ncp)
{
    int i, j;
    int k = c->rows;
    int n = a->cols;
    int err = 0;

    double aa;
    double dt[2];
#ifdef WIN32	
    double HUGE_win32 = libset_get_double(CONV_HUGE);
#else
	double HUGE = libset_get_double(CONV_HUGE);
#endif

    /* Initialize the cluster centers. Here, we arbitrarily make the
       first k data points cluster centers.
    */
    for (i=0; i<k; i++)  {
	for (j=0; j<n; j++)  {
	    gretl_matrix_set(c, i, j, gretl_matrix_get(a, i, j));
	}
    }

    /** At this point the basic set-up is complete. This is the point to
	which we'd return if iterating over random initial selection of
	the cluster centers.
    **/

    /* For each point i, find its two closest centers, ic1[i] and
       ic2[i].  Assign the point to ic1[i].
    */
    find_nearest_neighbors(a, c, ic1, ic2, dt);

    err = init_centers(nc, c, ic1, a);
    if (err) {
	/* there's at least one empty cluster */
	return err;
    }

    double temp;
    
    for (i=0; i<k; i++)  {
	/* compute centroids per cluster */
	aa = (double) (nc[i]);
	for (j=0; j<n; j++) {
	    temp = gretl_matrix_get(c, i, j) / aa;
	    gretl_matrix_set(c, i, j, temp);
	}
	/* initialize an, itran, ncp */
#ifdef WIN32	
	temp = (aa > 1.0) ? aa / (aa - 1.0) : HUGE_win32;
#else
	temp = (aa > 1.0) ? aa / (aa - 1.0) : HUGE;
#endif
	gretl_matrix_set(an, i, 0, temp);
	gretl_matrix_set(an, i, 1, aa / (aa + 1.0));
	itran[i] = 1;
	ncp[i] = -1;
    }

    return err;
}

	 
/* kmeans() carries out the K-means algorithm.

   @a (m x n): the data points
   @k: the assumed number of clusters
   @clustinfo: pointer to get information on the clusters
   @err: location to receive error code
*/

gretl_matrix *kmeans (const gretl_matrix *a, int k,
		      gretl_matrix **clustinfo,
		      int *err)
{
    gretl_vector *d;
    gretl_vector *clustid;
    gretl_matrix *an;
    gretl_matrix *c;
    int i, j, l;
    int *ic1;
    int *ic2;
    int indx;
    int *itran;
    int *live;
    int *nc;
    int *ncp;
    double temp;
    int maxiter = 128;
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

    *err = init(a, c, ic1, ic2, an, nc, itran, ncp);

    if (*err) {
	goto bailout;
    }
    
    indx = 0;
    *err = E_NOCONV;

    for (i=0; i<maxiter; i++)  {
	optra(a, c, ic1, ic2, nc, an, ncp, d, itran, live, &indx);
	if (indx == m) {
	    /* No optimal transfer in the last m steps: OK, stop */
	    *err = 0;
	    break;
	}
	qtran(a, c, ic1, ic2, nc, an, ncp, d, itran, &indx);
	if (k == 2) {
	    /* only two clusters: no need to repeat optra() */
	    *err = 0;
	    break;
	}
	/* reset ncp to 0 before re-entering optra() */
	for (l=0; l<k; l++) {
	    ncp[l] = 0;
	}
    }

    /* If we were iterating over random initial selection of centers,
       this is presumably the point at which we'd zip back to the top
       of the iteration (possibly wiping put a non-zero *err if it
       could have been bad luck?).
    */

    if (*err) {
	goto bailout;
    }

#if 0 /* not yet */
    double SST = compute_sst(a, c, ic1);
    fprintf(stderr, "HERE SST = %g\n", SST);
#endif    

    if (clustinfo != NULL) {
	gretl_matrix *cinfo;
	gretl_vector *sst;
	int reuse = 0;

	cinfo = get_clustinfo_target(clustinfo, k, n, &reuse);
	sst = compute_sst_full(a, c, ic1);

	for (i=0; i<k; i++) {
	    /* count of points in cluster i (first col) */
	    gretl_matrix_set(cinfo, i, 0, (double) nc[i]);
	    /* per-cluster centroids (middle cols) */
	    for (j=0; j<n; j++) {
		temp = gretl_matrix_get(c, i, j);
		gretl_matrix_set(cinfo, i, j+1, temp);
	    }
	    /* per-cluster SST (last col) */
	    gretl_matrix_set(cinfo, i, n+1, sst->val[i]);
	}

	gretl_matrix_free(sst);
	add_clustinfo_colnames(a, cinfo);
	if (!reuse) {
	    *clustinfo = cinfo;
	}
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
