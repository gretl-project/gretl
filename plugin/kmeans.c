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

   Adapted for gretl by Allin Cottrell and Jack Lucchetti, May 2026.
*/

#include "libgretl.h"
#include "libset.h"
#include "matrix_extra.h"

typedef enum InitFlag_ {
    INIT_HW,
    INIT_FAST,
    INIT_USER,
    INIT_RAND,
    INIT_FIN
} InitFlag;

typedef struct hw_info_ {
    const gretl_matrix *a; /* data points, m x n */
    gretl_matrix *c;       /* cluster centroids, k x n */
    gretl_matrix *cmin;    /* for random sequence */
    double *ameans;        /* column means of @a */
    int *ic1;              /* best centroid per point */
    int *ic2;              /* next best centroid per point */
    int m;                 /* number of observations */
    int n;                 /* dimension of observations */
    int k;                 /* number of clusters */
} hw_info;

static double *get_ameans (const gretl_matrix *a)
{
    double *ameans = malloc(a->cols * sizeof *ameans);
    int i, j;

    for (j=0; j<a->cols; j++) {
	ameans[j] = 0.0;
	for (i=0; i<a->rows; i++) {
	    ameans[j] += gretl_matrix_get(a, i, j);
	}
	ameans[j] /= a->rows;
    }

    return ameans;
}

static int build_hw_info (hw_info *hw,
			  const gretl_matrix *a,
			  const gretl_matrix *c0,
			  int k,
			  int rand_starts,
			  InitFlag *iflag)
{
    int err = 0;

    hw->m = a->rows;
    hw->n = a->cols;
    hw->k = k;
    hw->a = a;
    if (c0 != NULL) {
	hw->c = gretl_matrix_copy(c0);
	*iflag = INIT_USER; /* override */
    } else {
	hw->c = gretl_matrix_alloc(hw->k, hw->n);
    }
    hw->ic1 = malloc(2 * hw->m * sizeof *hw->ic1);
    hw->ic2 = hw->ic1 + hw->m;
    if (rand_starts > 0) {
	/* recorder for "best so far" */
	hw->cmin = gretl_matrix_alloc(hw->k, hw->n);
    } else {
	hw->cmin = NULL;
    }
    hw->ameans = get_ameans(a);

    return err;
}

static void destroy_hw_info (hw_info *hw)
{
    gretl_matrix_free(hw->c);
    gretl_matrix_free(hw->cmin);
    free(hw->ameans);
    free(hw->ic1);
}

static int init_centers (hw_info *hw, int *nc)
{
    int i, j, l;
    int err = 0;
    double tmp;

    for (l=0; l<hw->k; l++) {
	nc[l] = 0;
    }

    gretl_matrix_zero(hw->c);

    for (i=0; i<hw->m; i++) {
	l = hw->ic1[i] - 1;
	nc[l]++;
	for (j=0; j<hw->n; j++) {
	    tmp = gretl_matrix_get(hw->c, l, j) + gretl_matrix_get(hw->a, i, j);
	    gretl_matrix_set(hw->c, l, j, tmp);
	}
    }

    /* Check to see if there is any empty cluster at this stage */
    for (l=0; l<hw->k; l++)  {
	if (nc[l] == 0)  {
	    err = E_ZERO;
	    gretl_errmsg_sprintf(_("Cluster %d is empty; stopping."), l+1);
	    break;
	}
    }

    return err;
}

/* Get the k-vector of per-cluster SSTs */

static double *compute_sst_full (hw_info *hw)
{
    double *sst;
    double d;
    int i, j, l;

    sst = malloc(hw->k * sizeof *sst);
    for (l=0; l<hw->k; l++) {
	sst[l] = 0.0;
    }

    for (i=0; i<hw->m; i++) {
	l = hw->ic1[i] - 1;
	for (j=0; j<hw->n; j++) {
	    d = gretl_matrix_get(hw->a, i, j) - gretl_matrix_get(hw->c, l, j);
	    sst[l] += d * d;
	}
    }

    return sst;
}

/* Just get the sum of per-cluster SSTs */

static double compute_sst (hw_info *hw)
{
    double sst = 0;
    double d;
    int i, j, l;

    for (i=0; i<hw->m; i++) {
	l = hw->ic1[i] - 1;
	for (j=0; j<hw->n; j++) {
	    d = gretl_matrix_get(hw->a, i, j) - gretl_matrix_get(hw->c, l, j);
	    sst += d * d;
	}
    }

    return sst;
}

/* The initialization of @c suggested in the last paragraph of
   Hartigan and Wong (1979): sort the data points by euclidean
   distance from the global centroid, and select k evenly spaced
   points from the sorted array.
*/

static int hartigan_wong_init (hw_info *hw)
{
    gretl_matrix *dmat;
    gretl_matrix *s;
    double d, d2i, cmean;
    int i, j, l, r;
    int step;
    int err = 0;

    dmat = gretl_zero_matrix_new(hw->m, 2);

    for (j=0; j<hw->n; j++) {
	cmean = hw->ameans[j];
	for (i=0; i<hw->m; i++) {
	    d = gretl_matrix_get(hw->a, i, j) - cmean;
	    d2i = gretl_matrix_get(dmat, i, 0) + d * d;
	    gretl_matrix_set(dmat, i, 0, d2i);
	    if (j == 0) {
		gretl_matrix_set(dmat, i, 1, (double) i);
	    }
	}
    }

    s = gretl_matrix_sort_by_column(dmat, 0, &err);
    step = hw->m / hw->k;

    for (l=0; l<hw->k; l++) {
	r = gretl_matrix_get(s, l * step, 1);
	for (j=0; j<hw->n; j++) {
	    d = gretl_matrix_get(hw->a, r, j);
	    gretl_matrix_set(hw->c, l, j, d);
	}
    }

    gretl_matrix_free(dmat);
    gretl_matrix_free(s);

    return err;
}

static double global_sst (hw_info *hw)
{
    double sst = 0;
    double cmean, d;
    int i, j;

    for (j=0; j<hw->n; j++) {
	cmean = hw->ameans[j];
	for (i=0; i<hw->m; i++) {
	    d = gretl_matrix_get(hw->a, i, j) - cmean;
	    sst += d * d;
	}
    }

    return sst;
}

static void get_k_random_candidates (hw_info *hw)
{
    gretl_vector *v = gretl_vector_alloc(hw->k);
    double tmp;
    int i, j, r;

    /* fill @v with @k draws from 1,2,...,@m without replacement */
    fill_permutation_vector(v, hw->m);

    for (i=0; i<hw->k; i++)  {
	r = (int) v->val[i] - 1;
	for (j=0; j<hw->n; j++)  {
	    tmp = gretl_matrix_get(hw->a, r, j);
	    gretl_matrix_set(hw->c, i, j, tmp);
	}
    }

    gretl_matrix_free(v);
}

static double compute_ith_distance (hw_info *hw,
				    int i, int l1)

{
    double d, dist = 0.0;
    int j;

    for (j=0; j<hw->n; j++) {
	d = gretl_matrix_get(hw->a, i, j) - gretl_matrix_get(hw->c, l1-1, j);
	dist += d * d;
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

   @hw->a (m x n): the data points
   @hw->c (k x n): the cluster centers
   @hw->ic1 (m): the cluster to which each point is assigned
   @hw->ic2 (m): the cluster to which each point is most likely
     to be transferred to at each step
   @nc (k): the number of points in each cluster
   @an (k x 2)
   @ncp (k)
   @d (m)
   @itran (k)
   @live (k)
   @indx: the number of steps since a transfer took place
*/

static void optra (hw_info *hw, int *nc, gretl_matrix *an,
		   int ncp[], gretl_vector *d, int itran[],
		   int live[], int *indx)
{
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
    for (l=0; l<hw->k; l++) {
	if (itran[l] == 1) {
	    live[l] = hw->m + 1;
	}
    }

    for (i=0; i<hw->m; i++) {
	*indx = *indx + 1;
	l1 = hw->ic1[i];
	l2 = hw->ic2[i];
	ll = l2;

	/* If point i is the only member of cluster l1, no transfer */

	if (nc[l1-1] > 1) {
	    /* If l1 has not yet been updated in this stage, no need to
	       re-compute d[i].
	    */
	    if (ncp[l1-1] != 0) {
		de = compute_ith_distance(hw, i, l1);
		gretl_vector_set(d, i, de * gretl_matrix_get(an, l1-1, 0));
	    }

	    /* Find the cluster with minimum r2 */
	    da = compute_ith_distance(hw, i, l2);
	    r2 = da * gretl_matrix_get(an, l2-1, 1);

	    for (l=1; l<=hw->k; l++) {
		/* If live[l1] <= i, then l1 is not in the live set.  If
		   this is true, we only need to consider clusters that
		   are in the live set for possible transfer of point i.
		   Otherwise, we need to consider all possible clusters.
		*/
		if ((i < live[l1-1]-1 || i < live[l2-1]-1) && l != l1 && l != ll) {
		    rr = r2 / gretl_matrix_get(an, l-1, 1);
		    dc = compute_ith_distance(hw, i, l);
		    if (dc < rr)  {
			r2 = dc * gretl_matrix_get(an, l-1, 1);
			l2 = l;
		    }
		}
	    }

	    /* If no transfer is necessary, l2 is the new ic2[i] */
	    if (gretl_vector_get(d, i) <= r2) {
		hw->ic2[i] = l2;
	    } else {
		/* Update cluster centers, line, ncp and an for clusters
		   l1 and l2, and update ic1[i] and ic2[i].
		*/
		*indx = 0;
		live[l1-1] = hw->m + i;
		live[l2-1] = hw->m + i;
		ncp[l1-1] = i + 1;
		ncp[l2-1] = i + 1;
		al1 = (double) (nc[l1-1]);
		alw = al1 - 1.0;
		al2 = (double) (nc[l2-1]);
		alt = al2 + 1.0;
		for (j=0; j<hw->n; j++) {
		    aij = gretl_matrix_get(hw->a, i, j);
		    tmp = (gretl_matrix_get(hw->c, l1-1, j) * al1 - aij) / alw;
		    gretl_matrix_set(hw->c, l1-1, j, tmp);
		    tmp = (gretl_matrix_get(hw->c, l2-1, j) * al2 + aij) / alt;
		    gretl_matrix_set(hw->c, l2-1, j, tmp);
		}
		nc[l1-1] = nc[l1-1] - 1;
		nc[l2-1] = nc[l2-1] + 1;
		update_an(an, l1, l2, al1, al2);
		hw->ic1[i] = l2;
		hw->ic2[i] = l1;
	    }
	}

	if (*indx == hw->m) {
	    return;
	}
    }

    /* itran(l) = 0 before entering qtran(). Also, live(l) has to be
       decreased by m before re-entering optra().
    */
    for (l=1; l<=hw->k; l++) {
	itran[l-1] = 0;
	live[l-1] = live[l-1] - hw->m;
    }
}

/* qtran() carries out the quick transfer stage.

   @ic1 and @ic2 record the cluster to which point i belongs and the
   cluster to which it is most likely to be transferred, respectively.

   For each point i, ic1[i] and ic2[i] are switched, if necessary, to
   reduce the within-cluster sum of squares.  The cluster centers are
   updated after each step. This involves as many passes through the
   data as are required.

   @hw->a (m x n): the data points
   @hw->c (k x n): the cluster centers
   @hw->ic1 (m): the cluster to which each point is assigned
   @hw->ic2 (m): the cluster to which each point is most likely
     to be transferred to at each step
   @nc (k): the number of points in each cluster
   @an (k x 2)
   @ncp (k)
   @d (m)
   @itran (k)
   @indx: the number of steps since a transfer took place
*/

static void qtran (hw_info *hw, int *nc, gretl_matrix *an,
		   int ncp[], gretl_vector *d, int itran[],
		   int *indx)
{
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
	for (i=0; i<hw->m; i++) {
	    icoun++;
	    istep++;
	    l1 = hw->ic1[i];
	    l2 = hw->ic2[i];

	    /* If point i is the only member of cluster l1, no transfer */
	    if (nc[l1-1] > 1) {
		/* If ncp[l1] < istep, no need to recompute distance
		   from point i to cluster l1.  Note that if cluster l1
		   is last updated exactly m steps ago, we still need to
		   compute the distance from point i to cluster l1.
		*/
		if (istep <= ncp[l1-1]) {
		    di = compute_ith_distance(hw, i, l1);
		    gretl_vector_set(d, i, di * gretl_matrix_get(an, l1-1, 0));
		}
		/* If ncp[l1] <= istep and ncp[l2] <= istep, there will
		   be no transfer of point i at this step.
		*/
		if (istep < ncp[l1-1] || istep < ncp[l2-1]) {
		    r2 = gretl_vector_get(d, i) / gretl_matrix_get(an, l2-1, 1);
		    di = compute_ith_distance(hw, i, l2);
		    /* Update cluster centers, ncp, nc, itran, and an
		      for clusters l1 and l2.  Also update ic1[i] and
		      ic2[i].  If any updating occurs in this stage,
		      *indx is reset to 0.
		    */
		    if (di < r2) {
			icoun = 0;
			*indx = 0;
			itran[l1-1] = 1;
			itran[l2-1] = 1;
			ncp[l1-1] = istep + hw->m;
			ncp[l2-1] = istep + hw->m;
			al1 = (double) (nc[l1-1]);
			alw = al1 - 1.0;
			al2 = (double) (nc[l2-1]);
			alt = al2 + 1.0;
			for (j=0; j<hw->n; j++) {
			    aij = gretl_matrix_get(hw->a, i, j);
			    tmp = (gretl_matrix_get(hw->c, l1-1, j) * al1 - aij) / alw;
			    gretl_matrix_set(hw->c, l1-1, j, tmp);
			    tmp = (gretl_matrix_get(hw->c, l2-1, j) * al2 + aij) / alt;
			    gretl_matrix_set(hw->c, l2-1, j, tmp);
			}
			nc[l1-1] = nc[l1-1] - 1;
			nc[l2-1] = nc[l2-1] + 1;
			update_an(an, l1, l2, al1, al2);
			hw->ic1[i] = l2;
			hw->ic2[i] = l1;
		    }
		}
	    }
	    /* If no reallocation took place in the last m steps,
	       we're finished.
	    */
	    done = icoun == hw->m;
	    if (done) {
		break;
	    }
	}
    }
}

static void find_nearest_neighbors (hw_info *hw)
{
    double da, db, dc;
    double tmp;
    double dt[2];
    int i, j, l, il;

    for (i=0; i<hw->m; i++) {
	hw->ic1[i] = 1;
	hw->ic2[i] = 2;
	for (il=0; il<2; il++) {
	    dt[il] = 0.0;
	    for (j=0; j<hw->n; j++) {
		tmp = gretl_matrix_get(hw->c, il, j);
		da = gretl_matrix_get(hw->a, i, j) - tmp;
		dt[il] = dt[il] + da * da;
	    }
	}
	if (dt[1] < dt[0]) {
	    hw->ic1[i] = 2;
	    hw->ic2[i] = 1;
	    tmp = dt[0];
	    dt[0] = dt[1];
	    dt[1] = tmp;
	}
	for (l=2; l<hw->k; l++) {
	    db = 0.0;
	    for (j=0; j<hw->n; j++) {
		tmp = gretl_matrix_get(hw->c, l, j);
		dc = gretl_matrix_get(hw->a, i, j) - tmp;
		db = db + dc * dc;
	    }
	    if (db < dt[1]) {
		if (dt[0] <= db) {
		    dt[1] = db;
		    hw->ic2[i] = l+1;
		} else {
		    dt[1] = dt[0];
		    hw->ic2[i] = hw->ic1[i];
		    dt[0] = db;
		    hw->ic1[i] = l+1;
		}
	    }
	}
    }
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

static int kmeans_init (hw_info *hw,
			gretl_matrix *an,
			int *nc, int *itran,
			int *ncp, InitFlag iflag)
{
    int i, j;
    double aa, tmp;
    double huge = libset_get_double(CONV_HUGE);
    int err = 0;

    if (iflag == INIT_FAST) {
	/* Make the first k data points the cluster centers */
	for (i=0; i<hw->k; i++) {
	    for (j=0; j<hw->n; j++) {
		tmp = gretl_matrix_get(hw->a, i, j);
		gretl_matrix_set(hw->c, i, j, tmp);
	    }
	}
    } else if (iflag == INIT_RAND) {
	get_k_random_candidates(hw);
    } else if (iflag == INIT_HW) {
	hartigan_wong_init(hw);
    } else {
	; /* INIT_USER or INIT_FIN: use the incoming @c as is */
    }

    /* For each point i, find its two closest centers, ic1[i] and
       ic2[i].  Assign the point to ic1[i].
    */
    find_nearest_neighbors(hw);

    err = init_centers(hw, nc);
    if (err) {
	/* there's at least one empty cluster */
	return err;
    }

    for (i=0; i<hw->k; i++)  {
	/* compute centroids per cluster */
	aa = (double) (nc[i]);
	for (j=0; j<hw->n; j++) {
	    tmp = gretl_matrix_get(hw->c, i, j) / aa;
	    gretl_matrix_set(hw->c, i, j, tmp);
	}
	/* initialize an, itran, ncp */
	tmp = (aa > 1.0) ? aa / (aa - 1.0) : huge;
	gretl_matrix_set(an, i, 0, tmp);
	gretl_matrix_set(an, i, 1, aa / (aa + 1.0));
	itran[i] = 1;
	ncp[i] = -1;
    }

    return err;
}

static int check_opts (gretl_bundle *b,
		       int *rand_starts,
		       int *verbosity,
		       InitFlag *iflag)
{
    int err = 0;

    if (gretl_bundle_has_key(b, "rand_starts")) {
	*rand_starts = gretl_bundle_get_int(b, "rand_starts", &err);
    }
    if (gretl_bundle_has_key(b, "verbosity")) {
	*verbosity = gretl_bundle_get_int(b, "verbosity", &err);
    }
    if (gretl_bundle_has_key(b, "init")) {
	int itype = gretl_bundle_get_int(b, "init", &err);

	if (itype == 2) {
	    *iflag = INIT_FAST;
	} else if (itype != 1) {
	    err = E_INVARG;
	}
    }

    return err;
}

/* kmeans() carries out the K-means algorithm.

   @a (m x n): the data points
   @k: the assumed number of clusters (or 0)
   @c0: initial specification of clusters (or NULL)
   @clustinfo: pointer to get information on the clusters
   @err: location to receive error code
*/

gretl_bundle *kmeans (const gretl_matrix *a,
		      int k, const gretl_matrix *c0,
		      const gretl_bundle *opts,
		      PRN *prn, int *err)
{
    gretl_bundle *ret = NULL;
    hw_info hw = {0};
    gretl_vector *d = NULL;
    gretl_vector *clustid = NULL;
    gretl_matrix *an = NULL;
    double tmp;
    double SSTmin, SST;
    int *iwork;
    int *nc;
    int *ncp;
    int *itran;
    int *live;
    int indx;
    int i, j, l;
    int maxiter = 128;
    int m = a->rows;
    int n = a->cols;
    int ri = 0;
    InitFlag iflag = INIT_HW;
    int rand_starts = 0;
    int verbosity = 0;

    /* initial checks */
    if (c0 != NULL) {
	k = c0->rows;
	if (c0->cols != n) {
	    *err = E_NONCONF;
	}
    }
    if (!*err && (k <= 1 || m <= k)) {
	*err = E_NONCONF;
    }
    if (!*err && opts != NULL) {
	*err = check_opts((gretl_bundle *) opts, &rand_starts,
			  &verbosity, &iflag);
    }
    if (*err) {
	return NULL;
    }

    build_hw_info(&hw, a, c0, k, rand_starts, &iflag);

    if (verbosity) {
	pprintf(prn, "_kmeans: m=%d, n=%d, k=%d, initial centers %s\n",
		m, n, k, c0 == NULL ? "automatic" : "user-specified");
	if (c0 == NULL) {
	    pprintf(prn, "automatic method: %s\n", iflag == INIT_FAST ?
		    "simple, fast" : "Hartigan-Wong");
	}
	pprintf(prn, "%d random starts requested\n", rand_starts);
    }

    /* integer-valued workspace */
    iwork = malloc(4 * k * sizeof *iwork);
    nc = iwork;
    ncp = nc + k;
    itran = ncp + k;
    live = itran + k;

    d = gretl_vector_alloc(m);
    an = gretl_matrix_alloc(k, 2);
    clustid = gretl_column_vector_alloc(m);

    if (iwork == NULL || d == NULL || an == NULL || clustid == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    SSTmin = libset_get_double(CONV_HUGE);

 start_outer_loop:

    *err = kmeans_init(&hw, an, nc, itran, ncp, iflag);
    if (*err) {
	goto bailout;
    }

    indx = 0;
    *err = E_NOCONV;

    for (i=0; i<maxiter; i++)  {
	optra(&hw, nc, an, ncp, d, itran, live, &indx);
	if (indx == m) {
	    /* No optimal transfer in the last m steps: OK, stop */
	    *err = 0;
	    break;
	}
	qtran(&hw, nc, an, ncp, d, itran, &indx);
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

    if (*err) {
	goto bailout;
    }

    if (iflag != INIT_FIN) {
	if (verbosity > 2) {
	    pprintf(prn, "  point transfers converged in %d iterations\n", i+1);
	}
	if (rand_starts > 0 && ri <= rand_starts) {
	    SST = compute_sst(&hw);
	    if (verbosity > 1) {
		if (ri == 0) {
		    pprintf(prn, "initial SST = %g\n", SST);
		} else {
		    pprintf(prn, "start %d: SST = %g\n", ri, SST);
		}
	    }
	    if (SST < SSTmin) {
		SSTmin = SST;
		gretl_matrix_copy_values(hw.cmin, hw.c);
	    }
	    iflag = INIT_RAND;
	    ri++;
	    goto start_outer_loop;
	}
    }

    if (iflag == INIT_RAND) {
	if (verbosity > 1) {
	    pprintf(prn, "Minimum SST = %g\n", SSTmin);
	}
	gretl_matrix_copy_values(hw.c, hw.cmin);
	iflag = INIT_FIN;
	goto start_outer_loop;
    }

    ret = gretl_bundle_new();

    /* build the clustinfo matrix */
    {
	gretl_matrix *cinfo = gretl_zero_matrix_new(k, n+2);
	double *sst = compute_sst_full(&hw);

	for (i=0; i<k; i++) {
	    /* count of points in cluster i (first col) */
	    gretl_matrix_set(cinfo, i, 0, (double) nc[i]);
	    /* per-cluster centroids (middle cols) */
	    for (j=0; j<n; j++) {
		tmp = gretl_matrix_get(hw.c, i, j);
		gretl_matrix_set(cinfo, i, j+1, tmp);
	    }
	    /* per-cluster SST (last col) */
	    gretl_matrix_set(cinfo, i, n+1, sst[i]);
	}

	free(sst);
	add_clustinfo_colnames(a, cinfo);
	if (verbosity) {
	    gretl_matrix_print_to_prn(cinfo, "Cluster information:", prn);
	}
	gretl_bundle_donate_data(ret, "clustinfo", cinfo,
				 GRETL_TYPE_MATRIX, 0);
    }

    /* fill the output vector */
    for (i=0; i<m; i++) {
	gretl_vector_set(clustid, i, hw.ic1[i]);
    }
    gretl_bundle_donate_data(ret, "clustid", clustid,
			     GRETL_TYPE_MATRIX, 0);
    /* for reference, add the global SST */
    gretl_bundle_set_scalar(ret, "global_SST", global_sst(&hw));

 bailout:

    destroy_hw_info(&hw);
    gretl_matrix_free(an);
    gretl_vector_free(d);
    free(iwork);

    return ret;
}
