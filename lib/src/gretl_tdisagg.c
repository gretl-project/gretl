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

/* support for temporal disaggregation functionality in gretl */

#include "libgretl.h"
#include "gretl_midas.h"
#include "gretl_tdisagg.h"

/* Convert from @x or @list (whichever is non-NULL) to matrix,
   respecting the current sample range. If @cfac > 1 we take
   only every cfac-th row from the data source.
*/

gretl_matrix *tdisagg_matrix_from_series (const double *x,
					  int xnum,
					  const int *list,
					  const DATASET *dset,
					  int cfac, int *err)
{
    gretl_matrix *m = NULL;
    char **S = NULL;
    int k = list != NULL ? list[0] : 1;
    int T = dset->t2 - dset->t1 + 1;
    int i, s, t;
    int serr = 0;

    if (cfac > 1) {
	int rem = T % cfac;

	T = T / cfac + (rem > 0);
    }

    m = gretl_matrix_alloc(T, k);
    if (m == NULL) {
	return NULL;
    }

    if (list != NULL) {
	for (i=0; i<k; i++) {
	    x = dset->Z[list[i+1]];
	    s = dset->t1;
	    for (t=0; t<T; t++) {
		gretl_matrix_set(m, t, i, x[s]);
		s += cfac;
	    }
	}
	S = gretl_list_get_names_array(list, dset, &serr);
    } else {
	s = dset->t1;
	for (t=0; t<T; t++) {
	    m->val[t] = x[s];
	    s += cfac;
	}
	if (xnum > 0) {
	    S = strings_array_new(1);
	    if (S != NULL) {
		S[0] = gretl_strdup(dset->varname[xnum]);
	    }
	}
    }

    gretl_matrix_set_t1(m, dset->t1);
    gretl_matrix_set_t2(m, dset->t2);

    if (S != NULL) {
	gretl_matrix_set_colnames(m, S);
    }

    return m;
}

/**
 * matrix_tdisagg:
 * @Y: T x g: holds the original data to be disaggregated, series
 * in columns.
 * @X: (optionally) holds covariates of Y at the higher frequency.
 * @s: the expansion factor: e.g. 3 for quarterly to monthly,
 * 4 for annual to quarterly, 12 for annual to monthly.
 * @b: bundle containing optional arguments, or NULL.
 * @r: bundle for retrieving details, or NULL.
 * @err: location to receive error code.
 *
 * Temporal disaggregation, via regression or Denton-type method.
 *
 * Returns: matrix containing the disaggregated series, or
 * NULL on failure.
 */

gretl_matrix *matrix_tdisagg (const gretl_matrix *Y,
			      const gretl_matrix *X,
			      int s, gretl_bundle *b,
			      gretl_bundle *res,
			      DATASET *dset,
			      PRN *prn, int *err)
{
    gretl_matrix *(*tdisagg) (const gretl_matrix *,
			      const gretl_matrix *,
			      int, gretl_bundle *,
			      gretl_bundle *,
			      DATASET *,
			      PRN *, int *);
    gretl_matrix *ret = NULL;

    if (s <= 1) {
	*err = E_INVARG;
	return NULL;
    }

    if (X != NULL) {
	double xs = X->rows / (double) Y->rows;

	if (xs < s) {
	    *err = E_INVARG;
	    return NULL;
	} else if (X->is_complex) {
	    *err = E_CMPLX;
	    return NULL;
	}
    }

    tdisagg = get_plugin_function("time_disaggregate");

    if (tdisagg == NULL) {
	*err = E_FOPEN;
    } else {
	ret = (*tdisagg) (Y, X, s, b, res, dset, prn, err);
    }

    return ret;
}

/* For backward compatibility only: support the old chowlin()
   function, with "avg" as default aggregation type.
*/

gretl_matrix *matrix_chowlin (const gretl_matrix *Y,
			      const gretl_matrix *X,
			      int s, int *err)
{
    gretl_bundle *b = gretl_bundle_new();
    gretl_matrix *m = NULL;

    gretl_bundle_set_string(b, "aggtype", "avg");
    m = matrix_tdisagg(Y, X, s, b, NULL, NULL, NULL, err);
    gretl_bundle_destroy(b);

    return m;
}

/* Called in the context of tdisagg driver, when we're trying
   to determine if the target y (series or list) needs
   compressing. This will be the case if y just repeats
   low-frequency values.
*/

static int tdisagg_get_y_compression (int ynum, int xconv,
                                      int s, GretlType targ,
				      DATASET *dset)
{
    if (ynum > 0 && series_get_orig_pd(dset, ynum)) {
        return s;
    } else if (targ == GRETL_TYPE_SERIES) {
        return s;
    } else if (dset->pd == 4 && s == 4) {
        return s;
    } else if (dset->pd == 12 && s == 12) {
        return s;
    } else if (xconv == 1) {
        /* X was given as a series or list */
        return s;
    } else {
        return 1;
    }
}

/* tdisagg: advance the sample start if the first
   high-frequency X observation is not aligned to
   the first sub-period.
*/

static int maybe_advance_t1 (int t1, const DATASET *dset)
{
    int subper = 0;

    date_maj_min(t1, dset, NULL, &subper);
    if (subper > 1) {
        t1 += dset->pd - subper + 1;
    }
    return t1;
}

/* tdisagg: when both Y and X are dataset objects, try to
   restrict the sample ranges appropriately and ensure
   alignment at the start of the data.
*/

static int tdisagg_get_start_stop (int ynum, const int *ylist,
                                   int xnum, const int *xlist,
                                   const DATASET *dset,
                                   int cfac, int xmidas,
                                   int *start, int *ystop,
                                   int *xstop)
{
    int yvars[2] = {1, ynum};
    int xvars[2] = {1, xnum};
    int t1 = dset->t1;
    int t2 = dset->t2;
    int err = 0;

    if ((ynum == 0 && ylist == NULL) ||
        (xnum == 0 && xlist == NULL)) {
        /* can't do this */
        return 0;
    }

    if (ylist == NULL) {
        ylist = yvars;
    }
    if (xlist == NULL) {
        xlist = xvars;
    }

    err = list_adjust_sample(xlist, &t1, &t2, dset, NULL);

    if (!err && !xmidas) {
        t1 = maybe_advance_t1(t1, dset);
    }

    if (!err) {
        int yt1 = t1, yt2 = t2;
        int nmiss = 0;

        if (cfac > 1) {
            err = list_adjust_sample(ylist, &yt1, &yt2, dset, &nmiss);
        } else {
            err = list_adjust_sample(ylist, &yt1, &yt2, dset, NULL);
        }
        if (!err) {
            if (yt1 > t1) {
                t1 = yt1;
                if (!xmidas) {
                    t1 = maybe_advance_t1(t1, dset);
                }
            }
            *start = t1;
            *ystop = yt2;
            *xstop = t2;
        }
    }

    return err;
}

/* tdisagg: when Y is a dataset object and no stochastic
   X is given, try to restrict the sample range for Y
   appropriately.
*/

static int tdisagg_get_y_start_stop (int ynum, const int *ylist,
                                     const DATASET *dset, int cfac,
                                     int *t1, int *t2)
{
    int yvars[2] = {1, ynum};
    int err = 0;

    if (ynum == 0 && ylist == NULL) {
        /* can't do this */
        return 0;
    } else if (ylist == NULL) {
        ylist = yvars;
    }

    if (cfac > 1) {
        int nmiss = 0;

        err = list_adjust_sample(ylist, t1, t2, dset, &nmiss);
    } else {
        err = list_adjust_sample(ylist, t1, t2, dset, NULL);
    }

    return err;
}

int tdisagg_data_to_matrix (int yconv, int ynum,
			    const int *ylist, const double *yval,
			    int xconv, int xnum,
			    const int *xlist, const double *xval,
			    int xmidas, int fac, GretlType targ, 
			    gretl_matrix **X, gretl_matrix **Y,
			    DATASET *dset)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int yt2 = 0, xt2 = 0;
    int cfac = 1;
    int err = 0;

    if (yconv) {
	cfac = tdisagg_get_y_compression(ynum, xconv, fac, targ, dset);
    }
    if (yconv && (xconv || xmidas)) {
	err = tdisagg_get_start_stop(ynum, ylist, xnum, xlist,
				     dset, cfac, xmidas, &t1,
				     &yt2, &xt2);
	if (!err) {
	    dset->t1 = t1;
	}
    } else if (yconv) {
	err = tdisagg_get_y_start_stop(ynum, ylist, dset,
				       cfac, &t1, &t2);
	if (!err) {
	    dset->t1 = t1;
	    dset->t2 = t2;
	}
    }
    if (!err && yconv) {
	if (yt2 > 0) {
	    dset->t2 = yt2;
	}
	*Y = tdisagg_matrix_from_series(yval, ynum, ylist, dset,
				       cfac, &err);
    }
    if (!err && xconv) {
	if (xt2 > 0) {
	    dset->t2 = xt2;
	}
	*X = tdisagg_matrix_from_series(xval, xnum, xlist, dset,
				       1, &err);
    } else if (!err && xmidas) {
	if (xt2 > 0) {
	    dset->t2 = xt2;
	}
	*X = midas_list_to_vector(xlist, dset, &err);
    }

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}
