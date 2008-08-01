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

/* subsample.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_panel.h"
#include "cmd_private.h"
#include "dbread.h"
#include "gretl_scalar.h"

#define SUBDEBUG 0

/*
  The purpose of the static pointers below: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting cases
  rather than just moving the starting or ending points of the data
  range -- we create a new sub-dataset, and we need to keep the full
  dataset around so that it can be restored later.  The pointers fullZ
  and fullinfo are used to record the addresses of the full data
  matrix and DATAINFO struct respectively.

  In addition, peerinfo keeps track of the location of the DATAINFO
  struct associated with the back-up full dataset; by means of this,
  we can know when to free the full dataset and when not to (for
  instance, if we're freeing an auxiliary dataset).
*/

static double **fullZ;
static DATAINFO *fullinfo;
static DATAINFO *peerinfo;

#define SUBMASK_SENTINEL 127

static int smpl_get_int (const char *s, double ***pZ, DATAINFO *pdinfo,
			 int *err);

/* accessors for full dataset, when sub-sampled */

double ***fetch_full_Z (void)
{
    return &fullZ;
}

void reset_full_Z (double ***pZ)
{
#if SUBDEBUG
    fprintf(stderr, "reset_full_Z -> %p\n", (void *) *pZ);
#endif    
    fullZ = *pZ;
}

DATAINFO *fetch_full_datainfo (void)
{
    return fullinfo;
}

static int full_data_length (const DATAINFO *pdinfo)
{
    int n = 0;

    if (fullinfo != NULL) {
	n = fullinfo->n;
    } else if (pdinfo != NULL) {
	n = pdinfo->n;
    } 

    return n;
}

static int get_submask_length (const char *s)
{
    int n = 1;

    if (s == RESAMPLED) {
	n = 0;
    } else {
	while (*s != SUBMASK_SENTINEL) {
	    n++;
	    s++;
	}
    }

    return n;
}

void free_subsample_mask (char *s)
{
    if (s != RESAMPLED && s != NULL) {
	free(s);
    }
}

char *copy_subsample_mask (const char *src)
{
    char *ret = NULL;

    if (src == RESAMPLED) {
	ret = RESAMPLED;
    } else if (src != NULL) {
	int n = get_submask_length(src);

	ret = malloc(n * sizeof *ret);
	if (ret != NULL) {
	    memcpy(ret, src, n);
	}
    }

    return ret;
}

char *copy_datainfo_submask (const DATAINFO *pdinfo)
{
    char *mask = NULL;

    if (complex_subsampled()) {
	mask = copy_subsample_mask(pdinfo->submask);
    }

    return mask;
}

int write_model_submask (const MODEL *pmod, FILE *fp)
{
    int ret = 0;

    if (pmod->submask == RESAMPLED) {
	fputs("<submask length=\"0\"></submask>\n", fp);
	ret = 1;
    } else if (pmod->submask != NULL) {
	int i, n = get_submask_length(pmod->submask);

	fprintf(fp, "<submask length=\"%d\">", n);
	for (i=0; i<n; i++) {
	    fprintf(fp, "%d ", (int) pmod->submask[i]);
	}
	fputs("</submask>\n", fp);
	ret = 1;
    }

    return ret;
}

int write_datainfo_submask (const DATAINFO *pdinfo, FILE *fp)
{
    int ret = 0;

    if (pdinfo->submask == RESAMPLED) {
	unsigned int seed = get_resampling_seed();

	fprintf(fp, "<resample seed=\"%u\" n=\"%d\"/>\n", seed, pdinfo->n);
	ret = 1;
    } else if (complex_subsampled()) {
	int i, n = get_submask_length(pdinfo->submask);

	fprintf(fp, "<submask length=\"%d\">", n);
	for (i=0; i<n; i++) {
	    fprintf(fp, "%d ", (int) pdinfo->submask[i]);
	}
	fputs("</submask>\n", fp);
	ret = 1;
    }

    return ret;
}

int submask_cmp (const char *m1, const char *m2)
{
    if (m1 == RESAMPLED || m2 == RESAMPLED) {
	return m1 == RESAMPLED && m2 == RESAMPLED;
    }

    while (*m1 != SUBMASK_SENTINEL && *m2 != SUBMASK_SENTINEL) {
	if (*m1 != *m2) {
	    return 1;
	}
	m1++;
	m2++;
    }

    return 0;
}

static char *make_submask (int n)
{
    char *mask = calloc(n + 1, 1);
    
    if (mask != NULL) {
	mask[n] = SUBMASK_SENTINEL;
    }

    return mask;
}

void maybe_free_full_dataset (const DATAINFO *pdinfo)
{
    if (pdinfo == peerinfo) {
	if (fullZ != NULL) {
	    free_Z(fullZ, fullinfo);
	    fullZ = NULL;
	}
	if (fullinfo != NULL) {
	    clear_datainfo(fullinfo, CLEAR_SUBSAMPLE);
	    free(fullinfo);
	    fullinfo = NULL;
	}
	peerinfo = NULL;
    }
}

/* we do this on "restore full sample" */

static void 
relink_to_full_dataset (double ***pZ, DATAINFO *pdinfo)
{
    *pZ = fullZ;
    *pdinfo = *fullinfo;

#if SUBDEBUG
    fprintf(stderr, "relink_to_full_dataset: fullZ = %p, fullinfo = %p\n",
	    (void *) fullZ, (void *) fullinfo);
#endif

    fullZ = NULL;
    free(fullinfo);
    fullinfo = NULL;
    peerinfo = NULL;
}

/* sync malloced elements of the fullinfo struct that might
   have been moved via realloc
*/

static void sync_dataset_elements (const DATAINFO *pdinfo)
{
    fullinfo->varname = pdinfo->varname;
    fullinfo->varinfo = pdinfo->varinfo;
    fullinfo->descrip = pdinfo->descrip;
}

/* attach_subsample_to_model:
 * @pmod: model to which subsample should be attached.
 * @pdinfo: pointer to current dataset info.
 *
 * If the dataset is currently subsampled, record the subsample
 * information with the model so that it can be retrieved later.
 * 
 * Returns: 0 if the recording is not needed, or on success; non-zero
 * error code failure.
 */

int attach_subsample_to_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    int err = 0;

#if SUBDEBUG
    fprintf(stderr, "attach_subsample_to_model: fullZ = %p\n",
	    (void *) fullZ);
#endif

    if (fullZ != NULL) {
	/* sync, in case anything has moved */
	sync_dataset_elements(pdinfo);

	if (pmod->submask != NULL) {
	    free_subsample_mask(pmod->submask);
	}

	pmod->submask = copy_subsample_mask(pdinfo->submask);
	if (pmod->submask == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int resample_update_dataset (double ***RZ, DATAINFO *pdinfo)
{
    int newv = fullinfo->v;
    double **newZ = NULL;
    int drop = 0;
    int i, j, err = 0;

    /* carry back values of pre-existing scalars */

    for (i=1; i<fullinfo->v && i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    fullZ[i][0] = (*RZ)[i][0];
	} 
    }

    /* figure how many scalars and how many series have been added */

    for (i=fullinfo->v; i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    newv++;
	} else {
	    drop++;
	}
    }

    /* if new series have been added, drop them all */

    if (drop > 0) {
	int *list = gretl_list_new(drop);

	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    j = 1;
	    for (i=fullinfo->v; i<pdinfo->v; i++) {
		if (var_is_series(pdinfo, i)) {
		    list[j++] = i;
		}
	    }	    
	}
	
	if (list != NULL) {
	    err = dataset_drop_listed_variables(list, RZ, pdinfo, NULL, NULL);
	    free(list);
	}
    }

    /* if new scalars have been added, add them to fullZ */

    if (!err && newv > fullinfo->v) {
	newZ = realloc(fullZ, newv * sizeof *fullZ);
	if (newZ == NULL) {
	    newv = fullinfo->v;
	    err = E_ALLOC;
	} else {
	    fullZ = newZ;
	    j = fullinfo->v;
	    for (i=fullinfo->v; i<pdinfo->v && !err; i++) {
		if (var_is_scalar(pdinfo, i)) {
		    fullZ[j] = malloc(sizeof **newZ);
		    if (fullZ[j] == NULL) {
			err = E_ALLOC;
		    } else {
			fullZ[j][0] = (*RZ)[i][0];
			j++;
		    }
		}
	    }
	}
    }

    if (!err) {
	fullinfo->v = newv;
    }

    return err;
}

/* Apparatus for updating full dataset when restoring full sample
   after sub-sampling.  
*/

static void
update_full_data_values (const double **subZ, const DATAINFO *pdinfo)
{
    int i, s, t;

#if SUBDEBUG
    fprintf(stderr, "update_full_data_values: fullZ=%p, subZ=%p, pdinfo=%p\n",
	    (void *) fullZ, (void *) subZ, (void *) pdinfo);
#endif

    for (i=1; i<fullinfo->v && i<pdinfo->v; i++) {
	s = 0;
	if (var_is_scalar(pdinfo, i)) {
	    fullZ[i][0] = subZ[i][0];
	} else {
	    for (t=0; t<fullinfo->n; t++) {
		if (pdinfo->submask[t] == 1) {
		    fullZ[i][t] = subZ[i][s++];
		} else if (pdinfo->submask[t] == 'p') {
		    /* skip panel padding (?) */
		    s++;
		}
	    }
	}
    }
}

static int update_case_markers (const DATAINFO *pdinfo)
{
    int err = 0;

    if (pdinfo->markers && !fullinfo->markers) {
	dataset_allocate_obs_markers(fullinfo);
	if (fullinfo->S == NULL) {
	    err = 1;
	} else {
	    int t, subt = 0;

	    for (t=0; t<fullinfo->n; t++) {
		if (pdinfo->submask[t]) {
		    strcpy(fullinfo->S[t], pdinfo->S[subt++]);
		} else {
		    sprintf(fullinfo->S[t], "%d", t + 1);
		}
	    }
	}
    }	
	
    return err;
}

static int add_new_vars_to_full (const double **Z, DATAINFO *pdinfo)
{
    int V1 = pdinfo->v;
    int V0 = fullinfo->v;
    int N = fullinfo->n;
    double **newZ = NULL;
    int i, t, s;
    int err = 0;

    if (V1 <= V0) {
	return 0;
    }

    if (pdinfo->submask == NULL) {
	return E_NOMERGE;
    }

#if SUBDEBUG
    fprintf(stderr, "add_new_vars_to_full:\n");
    fprintf(stderr, " V1 = pdinfo->v = %d; V0 = fullinfo->v = %d\n",
	    V1, V0);
    fprintf(stderr, " Z = %p, fullZ = %p\n", (void *) Z, (void *) fullZ);
#endif

    /* allocate expanded data array */
    newZ = realloc(fullZ, V1 * sizeof *fullZ);

    if (newZ == NULL) {
	return E_ALLOC;
    } 

    fullZ = newZ;

    for (i=V0; i<pdinfo->v && !err; i++) {
	if (var_is_series(pdinfo, i)) {
	    fullZ[i] = malloc(N * sizeof **newZ);
	} else {
#if SUBDEBUG
	    fprintf(stderr, " new var %d (%s) is a scalar\n", i, 
		    pdinfo->varname[i]);
#endif
	    fullZ[i] = malloc(sizeof **newZ);
	}
	if (fullZ[i] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "After expansion, fullZ = %p (%d vars)\n", (void *) fullZ,
	    pdinfo->v);
#endif

    for (i=V0; i<pdinfo->v; i++) {
	if (var_is_scalar(pdinfo, i)) {
	    fullZ[i][0] = Z[i][0]; 
	} else {
	    s = 0;
	    for (t=0; t<N; t++) {
		fullZ[i][t] = (pdinfo->submask[t])? Z[i][s++] : NADBL;
	    }
	}
    }

    fullinfo->v = V1;

    return 0;
}

static char *make_current_sample_mask (DATAINFO *pdinfo)
{
    int n = full_data_length(pdinfo);
    char *currmask = NULL;
    int s, t;

    if (pdinfo->submask == NULL) {
	/* no pre-existing mask: not sub-sampled */
	currmask = make_submask(n);
	if (currmask != NULL) {
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		currmask[t] = 1;
	    }
	}
    } else {
	/* there's a pre-existing mask */
	currmask = copy_subsample_mask(pdinfo->submask);
	if (currmask != NULL) {
	    s = -1;
	    for (t=0; t<n; t++) {
		if (pdinfo->submask[t]) s++;
		if (s < pdinfo->t1 || s > pdinfo->t2) {
		    currmask[t] = 0;
		} 
	    }
	}
    }

    return currmask;
}

/* restore_full_sample: 
 * @pZ: pointer to data array.  
 * @pdinfo: pointer to dataset information.
 * @state: structure representing program execution state,
 * or %NULL. 
 *
 * Restore the full data range, undoing any sub-sampling that was
 * previously performed (either by shifting the endpoints of the
 * sample range or by selecting cases on some criterion or other).
 * However, if @state is not %NULL, and if it contains a submask,
 * the "full" range is relative to that mask.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restore_full_sample (double ***pZ, DATAINFO *pdinfo, ExecState *state)
{
    int err = 0;

    gretl_error_clear();

    if (!complex_subsampled()) {
	if (pdinfo->t1 != 0 || pdinfo->t2 != pdinfo->n - 1) {
	    pdinfo->t1 = 0;
	    pdinfo->t2 = pdinfo->n - 1;
#if SUBDEBUG
	    fprintf(stderr, "restore_full_sample: reset t1 and t2\n");
#endif
	}
	return 0;
    }

#if SUBDEBUG
    fprintf(stderr, "\nrestore_full_sample: pZ=%p, *pZ=%p, pdinfo=%p\n",
	    (void *) pZ, (void *) *pZ, (void *) pdinfo);
#endif

    /* beyond this point we are doing a non-trivial restoration
       of a stored full dataset, which has previously been
       subsampled, e.g., by some boolean criterion 
    */

    /* reattach malloc'd elements of subsampled dataset,
       which may have moved */
    sync_dataset_elements(pdinfo);

    if (pdinfo->submask == RESAMPLED) {
	/* not regular subsample but resampled dataset */
	resample_update_dataset(pZ, pdinfo);
    } else {
	/* update values for pre-existing series, which may have been
	   modified via genr etc */
	update_full_data_values((const double **) *pZ, pdinfo);

	/* if case markers were added when subsampled, carry them back */
	update_case_markers(pdinfo);

	/* delete any newly added hidden vars */
	err = dataset_destroy_hidden_variables(pZ, pdinfo, fullinfo->v);

	/* in case any new vars were added when subsampled, try to merge
	   them into the full dataset */
	if (!err) {
	    err = add_new_vars_to_full((const double **) *pZ, pdinfo);
	}
    }

    if (err == E_NOMERGE) {
        sprintf(gretl_errmsg, 
		_("Missing sub-sample information; can't merge data\n"));
    }

    /* destroy sub-sampled data array */
#if SUBDEBUG
    fprintf(stderr, "restore_full_sample: freeing Z at %p\n", (void *) *pZ);
#endif
    free_Z(*pZ, pdinfo);

#if SUBDEBUG
    fprintf(stderr, "restore_full_sample: clearing pdinfo at %p\n", 
	    (void *) pdinfo);
#endif
    clear_datainfo(pdinfo, CLEAR_SUBSAMPLE);

    relink_to_full_dataset(pZ, pdinfo);

    if (!err && state != NULL && state->submask != NULL) {
	/* "full" really means, relative to the original state */
	err = restrict_sample_from_mask(state->submask, pZ, pdinfo);
    }

    return err;
}

static int overlay_masks (char *targ, const char *src, int n)
{
    int i, sn = 0;
    
    for (i=0; i<n; i++) {
	if (targ[i] == 1 && src[i] == 1) {
	    targ[i] = 1;
	    sn++;
	} else {
	    targ[i] = 0;
	}
    }
	
    return sn;
}

static int 
make_missing_mask (const int *list, const double **Z, const DATAINFO *pdinfo,
		   char *mask)
{
    int i, vi, t;

    if (list != NULL) {
	/* check specified list of variables */
	for (t=0; t<pdinfo->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (var_is_series(pdinfo, vi) && na(Z[vi][t])) {
		    mask[t] = 0;
		    break;
		}
	    }
	}
    } else {	
	/* check all variables */
	for (t=0; t<pdinfo->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<pdinfo->v; i++) {
		if (var_is_series(pdinfo, i) && na(Z[i][t])) {
		    mask[t] = 0;
		    break;
		}
	    }
	}
    }

    return 0;
} 

static int copy_dummy_to_mask (char *mask, const double *x, int n)
{
    int t, err = 0;

    for (t=0; t<n && !err; t++) {
	if (x[t] == 1.0) {
	    mask[t] = 1;
	} else if (!na(x[t]) && x[t] != 0.0) { /* NA? */
	    err = 1;
	}
    }

    return err;
}

static int mask_from_temp_dummy (const char *s, double ***pZ, DATAINFO *pdinfo, 
				 char *mask)
{
    char formula[MAXLEN];
    double *x;
    int err = 0;

    *formula = '\0';
    strncat(formula, s, MAXLEN - 1);

    x = generate_series(formula, pZ, pdinfo, &err);

    if (!err) {
	err = copy_dummy_to_mask(mask, x, pdinfo->n);
	if (err) {
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), "mask");
	}
    }

    free(x);

    return err;
}

static int mask_from_dummy (const char *s, const double **Z, const DATAINFO *pdinfo,
			    char *mask)
{
    char dname[VNAMELEN] = {0};
    int dnum, err = 0;

    sscanf(s, "%15s", dname);
    dnum = varindex(pdinfo, dname);

    if (dnum == pdinfo->v) {
	sprintf(gretl_errmsg, _("Variable '%s' not defined"), dname);
	err = 1;
    } else {
	err = copy_dummy_to_mask(mask, Z[dnum], pdinfo->n);
	if (err) {
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), dname);
	}
    }

    return err;
}

/* how many observations are selected by the given 
   subsample mask? */

static int 
count_selected_cases (const char *mask, const DATAINFO *pdinfo)
{
    int i, n = 0;

    for (i=0; i<pdinfo->n; i++) {
	if (mask[i]) {
	    n++;
	}
    }

    return n;
}

/* panel: how many distinct cross-sectional units are included 
   in the masked subset of observations? */

static int 
count_panel_units (const char *mask, const DATAINFO *pdinfo)
{
    int u, ubak = -1;
    int i, n = 0;

    for (i=0; i<pdinfo->n; i++) {
	if (mask[i]) {
	    u = pdinfo->paninfo->unit[i];
	    if (u != ubak) {
		n++;
		ubak = u;
	    }
	}
    }

#if SUBDEBUG
    fprintf(stderr, "count_panel_units: got n = %d\n", n);
#endif

    return n;
}

/* construct mask for taking random sub-sample from dataset */

static int make_random_mask (const char *s, const char *oldmask, 
			     double ***pZ, DATAINFO *pdinfo,
			     char *mask)
{
    unsigned u;
    int i, subn;
    int cases = 0;
    int err = 0;

    subn = smpl_get_int(s, pZ, pdinfo, &err);

    if (subn <= 0 || subn >= pdinfo->n) {
	err = 1;
    } else if (oldmask != NULL) {
	int oldn = 0;

	for (i=0; i<pdinfo->n; i++) {
	    if (oldmask[i]) {
		oldn++;
	    }
	}
	if (subn >= oldn) {
	    err = 1;
	}
    }	

    if (err) {
	sprintf(gretl_errmsg, _("Invalid number of cases %d"), subn);
	return err;
    }	

    for (i=0; i<pdinfo->n; i++) {
	mask[i] = 0;
    }

    for (i=0; cases != subn; i++) {
	u = gretl_rand_int_max(pdinfo->n);
	if (oldmask == NULL || oldmask[u]) {
	    mask[u] = 1;
	}
	if (i >= subn - 1) {
	    cases = count_selected_cases(mask, pdinfo);
	}
    }

    return err;
}

int backup_full_dataset (double **Z, DATAINFO *pdinfo)
{
    fullZ = Z;

    if (fullinfo == NULL) {
	fullinfo = malloc(sizeof *fullinfo);
	if (fullinfo == NULL) {
	    return E_ALLOC;
	}
    }

    if (pdinfo != NULL) {
	*fullinfo = *pdinfo;
	peerinfo = pdinfo;
    } 

#if SUBDEBUG
    fprintf(stderr, "backup_full_dataset: fullZ = %p, fullinfo = %p\n",
	    (void *) fullZ, (void *) fullinfo);
#endif

    return 0;
}

int complex_subsampled (void)
{
    return (fullZ != NULL);
}

int get_full_length_n (void)
{
    int n = 0;

    if (fullinfo != NULL) {
	n = fullinfo->n;
    }

    return n;
}

/* When sub-sampling on some boolean criterion, check to see if we can
   meet the criterion by simply adjusting the endpoints of the sample
   range: life will be simpler if that is so.
*/

static int mask_contiguous (const char *mask,
			    const DATAINFO *pdinfo,
			    int *pt1, int *pt2)
{
    int t, t1 = 0, t2 = pdinfo->n - 1;
    int contig = 1;

    for (t=0; t<pdinfo->n; t++) {
	if (mask[t] == 0) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=pdinfo->n - 1; t>=0; t--) {
	if (mask[t] == 0) {
	    t2--;
	} else {
	    break;
	}
    }

    for (t=t1; t<=t2; t++) {
	if (mask[t] == 0) {
	    /* there's a hole inside the range */
	    contig = 0;
	    break;
	}
    }

    if (contig && dataset_is_panel(pdinfo)) {
	int n = t2 - t1 + 1;

	/* sample must leave a whole number of panel units; moreover,
	   to retain "panelness" this number must be greater than 1 
	*/
	if (t1 % pdinfo->pd != 0 || n % pdinfo->pd != 0) {
	    contig = 0;
	} else if (n == pdinfo->pd) {
	    contig = 0;
	}
    }

    if (contig) {
	*pt1 = t1;
	*pt2 = t2;
    }

    return contig;
}

static void 
copy_data_to_subsample (double **subZ, DATAINFO *subinfo,
			const double **Z, const DATAINFO *pdinfo,
			const char *mask)
{
    int i, t, st;

#if SUBDEBUG
    fprintf(stderr, "copy_data_to_subsample: subinfo = %p, pdinfo = %p\n",
	    (void *) subinfo, (void *) pdinfo);
#endif

    /* copy data values */
    for (i=1; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    st = 0;
	    for (t=0; t<pdinfo->n; t++) {
		if (mask == NULL) {
		    subZ[i][st++] = Z[i][t];
		} else if (mask[t] == 1) {
		    subZ[i][st++] = Z[i][t];
		} else if (mask[t] == 'p') {
		    /* panel padding */
		    subZ[i][st++] = NADBL;
		}
	    }
	} else {
	    subZ[i][0] = Z[i][0];
	}
    }

    /* copy observation markers, if any */
    if (pdinfo->markers && subinfo->markers) {
	st = 0;
	for (t=0; t<pdinfo->n; t++) {
	    if (mask == NULL || mask[t] == 1 || mask[t] == 'p') {
		strcpy(subinfo->S[st++], pdinfo->S[t]);
	    }
	}
    }

    strcpy(subinfo->stobs, "1");
    sprintf(subinfo->endobs, "%d", subinfo->n);
}

int get_restriction_mode (gretlopt opt)
{
    int mode = SUBSAMPLE_UNKNOWN;

    if (opt & OPT_M) {
	mode = SUBSAMPLE_DROP_MISSING;
    } else if (opt & OPT_R) {
	mode = SUBSAMPLE_BOOLEAN;
    } else if (opt & OPT_N) {
	mode = SUBSAMPLE_RANDOM;
    } else if (opt & OPT_O) {
	mode = SUBSAMPLE_USE_DUMMY;
    } 

    return mode;
}

static int 
make_panel_submask (char *mask, const DATAINFO *pdinfo, int *err)
{
    char *umask = NULL;
    char *pmask = NULL;
    int i, np = 0;

    umask = calloc(pdinfo->paninfo->nunits, 1);
    if (umask == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    pmask = calloc(pdinfo->paninfo->Tmax, 1);
    if (pmask == NULL) {
	free(umask);
	*err = E_ALLOC;
	return 0;
    }

    for (i=0; i<pdinfo->n; i++) {
	if (mask[i]) {
	    umask[pdinfo->paninfo->unit[i]] = 1;
	    pmask[pdinfo->paninfo->period[i]] = 1;
	}
    }

    for (i=0; i<pdinfo->n; i++) {
	if (!mask[i]) {
	    if (umask[pdinfo->paninfo->unit[i]] &&
		pmask[pdinfo->paninfo->period[i]]) {
		mask[i] = 'p'; /* mark as padding row */
		np++;
	    }
	}
    }

#if SUBDEBUG
    fprintf(stderr, "make_panel_submask: number of padding rows = %d\n",
	    np);
#endif

    free(umask);
    free(pmask);

    return np;
}

#define needs_string_arg(m) (m == SUBSAMPLE_RANDOM || \
	                     m == SUBSAMPLE_USE_DUMMY || \
                             m == SUBSAMPLE_BOOLEAN)

static int 
make_restriction_mask (int mode, const char *s, const int *list, 
		       double ***pZ, DATAINFO *pdinfo,
		       PRN *prn, const char *oldmask, char **pmask)
{
    char *mask = NULL;
    int sn = 0, err = 0;

    if (needs_string_arg(mode) && (s == NULL || *s == '\0')) {
	return E_ARGS;
    }

    mask = make_submask(pdinfo->n);
    if (mask == NULL) {
	return E_ALLOC;
    }

    /* construct subsample mask in one of several possible ways */
    if (mode == SUBSAMPLE_DROP_MISSING) {   
	err = make_missing_mask(list, (const double **) *pZ, pdinfo, mask);
    } else if (mode == SUBSAMPLE_RANDOM) {
	err = make_random_mask(s, oldmask, pZ, pdinfo, mask);
    } else if (mode == SUBSAMPLE_USE_DUMMY) {
	err = mask_from_dummy(s, (const double **) *pZ, pdinfo, mask);
    } else if (mode == SUBSAMPLE_BOOLEAN) {
	err = mask_from_temp_dummy(s, pZ, pdinfo, mask);
    } else {
	strcpy(gretl_errmsg, _("Sub-sample command failed mysteriously"));
	err = 1;
    }

    /* exit now on unrecoverable error */
    if (err) {
	free(mask);
	return err;
    }

    /* cumulate sample restrictions, if appropriate */
    if (oldmask != NULL && mode != SUBSAMPLE_RANDOM) {
	sn = overlay_masks(mask, oldmask, pdinfo->n);
    } else {
	sn = count_selected_cases(mask, pdinfo);
    }

    /* does this policy lead to an empty sample, or no change in the
       sample, perchance? */

    if (sn == 0) {
	strcpy(gretl_errmsg, _("No observations would be left!"));
	err = 1;
    } else if (sn == pdinfo->n) {
	/* not really an error, just a no-op */
	if (gretl_messages_on()) {
	    pputs(prn, _("No observations were dropped!"));
	    pputc(prn, '\n');
	}
	free(mask);
	mask = NULL;
    }

    if (err) {
	free(mask);
    } else {
	*pmask = mask;
    }

    return err;
}

/* This is also used in session.c, when re-establishing a previously
   sub-sampled data state on reopening a saved session, and also on
   exit from a user function when the dataset was sub-sampled on
   entry to the function.
*/

int 
restrict_sample_from_mask (char *mask, double ***pZ, DATAINFO *pdinfo)
{
    double **subZ = NULL;
    DATAINFO *subinfo;
    int err = 0;

    if (mask == RESAMPLED) {
	fprintf(stderr, "restrict_sample_from_mask: got RESAMPLED!\n");
	return E_DATA;
    }

    subinfo = datainfo_new();
    if (subinfo == NULL) {
	return E_ALLOC;
    }

    subinfo->n = count_selected_cases(mask, pdinfo);
    subinfo->v = pdinfo->v;

#if SUBDEBUG
    fprintf(stderr, "restrict_sample_from_mask: new subinfo=%p, %d obs in subsample\n",
	    (void *) subinfo, subinfo->n);
#endif

    if (dataset_is_panel(pdinfo)) {
	/* are we able to reconstitute a panel? */
	int npad, nunits = count_panel_units(mask, pdinfo);

	if (nunits > 1 && subinfo->n > nunits) {
	    /* add padding rows if need be */
	    npad = make_panel_submask(mask, pdinfo, &err);
	    if (err) {
		free(subinfo);
		return err;
	    }
	    subinfo->structure = STACKED_TIME_SERIES;
	    subinfo->n += npad;
	    subinfo->pd = subinfo->n / nunits;
	    /* note: panel indices are added below */
	} else if (nunits == 1 && subinfo->n == pdinfo->pd) {
	    /* time series for single panel unit */
	    subinfo->structure = SPECIAL_TIME_SERIES;
	}
    }

    /* set up the sub-sampled datainfo */
    if (start_new_Z(&subZ, subinfo, 1)) {
	free(subinfo);
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "incoming *pZ at %p, new subZ is at %p\n", 
	    (void *) *pZ, (void *) subZ);
#endif

    /* link (don't copy) varnames and descriptions, since these are
       not dependent on the series length */
    subinfo->varname = pdinfo->varname;
    subinfo->varinfo = pdinfo->varinfo;
    subinfo->descrip = pdinfo->descrip;

    if (subinfo->structure == STACKED_TIME_SERIES) {
#if SUBDEBUG
	fprintf(stderr, "panel subset: nT = %d, pd = %d\n", 
		subinfo->n, subinfo->pd);
#endif
	err = dataset_add_default_panel_indices(subinfo);
	if (err) {
	    free_Z(subZ, subinfo);
	    free(subinfo);
	    return E_ALLOC;
	}
    }

    /* set up case markers? */
    if (pdinfo->markers) {
	err = dataset_allocate_obs_markers(subinfo);
	if (err) {
	    free_Z(subZ, subinfo);
	    free(subinfo);
	    return E_ALLOC;
	}
    }

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(subZ, subinfo, (const double **) *pZ, 
			   pdinfo, mask);

    err = backup_full_dataset(*pZ, pdinfo);

    subinfo->submask = copy_subsample_mask(mask);

    /* switch pointers */
    *pZ = subZ;
    *pdinfo = *subinfo;
    free(subinfo);

    return err;
}

/* restrict_sample: 
 * @line: command line (or %NULL).  
 * @pZ: pointer to original data array.  
 * @pdinfo: data info pointer. 
 * @state: structure representing program state (or %NULL).
 * @list: list of variables in case of OPT_M (or %NULL).  
 * @opt: option flags.
 * @prn: printing apparatus.
 *
 * Sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values (OPT_M); or using as a mask a
 * specified dummy variable (OPT_O); or masking with a specified
 * boolean condition (OPT_R); or selecting at random (OPT_N).
 *
 * In case OPT_M a @list of variables may be supplied; in cases
 * OPT_O, OPT_R and OPT_N, @line must contain specifics.
 *
 * In case OPT_P is included, the restriction will rePlace any
 * existing sample restriction, otherwise the resulting restriction
 * will be the logical product of the new restriction and any
 * existing restriction.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample (const char *line, const int *list,
		     double ***pZ, DATAINFO *pdinfo, 
		     ExecState *state, gretlopt opt, 
		     PRN *prn)
{
    char *oldmask = NULL;
    char *mask = NULL;
    int free_oldmask = 0;
    int mode, err = 0;

    gretl_error_clear();

#if SUBDEBUG
    fprintf(stderr, "\nrestrict_sample: '%s'\n", line);
    fprintf(stderr, " pdinfo=%p, pZ=%p, *pZ=%p, state=%p\n", 
	    (void *) pdinfo, (void *) pZ, (void *) *pZ,
	    (void *) state);
#endif

    if (line != NULL) {
	/* skip the command word */
	line += strcspn(line, " ");
	line += strspn(line, " ");
    }

    mode = get_restriction_mode(opt);
    if (mode == SUBSAMPLE_UNKNOWN) {
	strcpy(gretl_errmsg, "Unrecognized sample command");
	return 1;
    }

    if (!(opt & OPT_P)) {
	/* not replacing but cumulating any existing restrictions */
	oldmask = make_current_sample_mask(pdinfo);
	if (oldmask == NULL) {
	    return E_ALLOC;
	}
	free_oldmask = 1;
    } else if (state != NULL && state->submask != NULL) {
	/* subsampling within a function, with incoming
	   restriction recorded in state
	*/
	oldmask = state->submask;
    }

    /* restore the full data range, for housekeeping purposes */
    err = restore_full_sample(pZ, pdinfo, NULL);
    if (err) {
	return err;
    }

    if (!err) {
	err = make_restriction_mask(mode, line, list, pZ, pdinfo, 
				    prn, oldmask, &mask);
    }

    if (!err && mask != NULL) {
	int t1 = 0, t2 = 0;
	int contig = 0;

	if (mode != SUBSAMPLE_RANDOM &&
	    (dataset_is_time_series(pdinfo) ||
	     dataset_is_panel(pdinfo))) {
	    contig = mask_contiguous(mask, pdinfo, &t1, &t2);
	}

	if (contig) {
	    pdinfo->t1 = t1;
	    pdinfo->t2 = t2;
	} else {
	    err = restrict_sample_from_mask(mask, pZ, pdinfo);
	}
    }

    free(mask);

    if (free_oldmask) {
	free(oldmask);
    }

    return err;
}

enum {
    SMPL_T1,
    SMPL_T2
};

static int panel_round (const DATAINFO *pdinfo, int t, int code)
{
    if (code == SMPL_T1) {
	while ((t + 1) % pdinfo->paninfo->Tmax != 1) {
	    t++;
	}
    } else {
	while ((t + 1) % pdinfo->paninfo->Tmax != 0) {
	    t--;
	}
    }

    return t;
}

static int smpl_get_int (const char *s, double ***pZ, DATAINFO *pdinfo,
			 int *err)
{
    int k = 0;

    if (integer_string(s)) {
	k = atoi(s);
    } else if (gretl_is_scalar(s)) {
	k = gretl_scalar_get_value(s);
    } else {
	int v = varindex(pdinfo, s);

	if (v < pdinfo->v && var_is_scalar(pdinfo, v)) {
	    k = (int) (*pZ)[v][0];
	} else {
	    k = (int) generate_scalar(s, pZ, pdinfo, err);
	}
    }

    return k;
}

static int get_sample_limit (char *s, double ***pZ, DATAINFO *pdinfo,
			     int code)
{
    int ret = -1;
    int err = 0;

    if (*s == '-' || *s == '+') {
	/* increment/decrement form */
	int incr = smpl_get_int(s + 1, pZ, pdinfo, &err);

	if (!err) {
	    if (*s == '-') {
		incr = -incr;
	    }
	    if (dataset_is_panel(pdinfo)) {
		incr *= pdinfo->paninfo->Tmax;
	    }
	    if (code == SMPL_T1) {
		ret = pdinfo->t1 + incr;
	    } else {
		ret = pdinfo->t2 + incr;
	    }
	}
    } else {
	/* absolute form */
	ret = get_t_from_obs_string(s, (const double **) *pZ, pdinfo);
	if (ret < 0) {
	    ret = smpl_get_int(s, pZ, pdinfo, &err);
	    /* convert to base 0 */
	    if (!err) ret--;
	}
	if (ret >= 0 && dataset_is_panel(pdinfo)) {
	    ret = panel_round(pdinfo, ret, code);
	}
    }

    return ret;
}

int set_sample (const char *line, double ***pZ, DATAINFO *pdinfo)
{
    int nf, new_t1 = pdinfo->t1, new_t2 = pdinfo->t2;
    char newstart[OBSLEN], newstop[OBSLEN];

    gretl_error_clear();

    line += strcspn(line, " ");
    line += strspn(line, " ");

    nf = count_fields(line);

#if SUBDEBUG
    fprintf(stderr, "set_sample: line='%s', nf=%d, pdinfo=%p\n", 
	    line, nf, (void *) pdinfo);
    if (pdinfo != NULL) {
	fprintf(stderr, "pdinfo->v = %d, pdinfo->n = %d\n",
		pdinfo->v, pdinfo->n);
    }
#endif

    if (nf == 2 && pdinfo->n == 0) {
	/* database special */
	return db_set_sample(line, pdinfo);
    }

    if (nf == 0) {
	/* no-op, just print the current sample */
	return 0;
    }
	
    if (nf == 1) {
	if (sscanf(line, "%10s", newstart) != 1) {
	    strcpy(gretl_errmsg, _("error reading smpl line"));
	    return 1;
	} else {
	    new_t1 = get_sample_limit(newstart, pZ, pdinfo, SMPL_T1);
	    if (new_t1 < 0 || new_t1 >= pdinfo->n) {
		strcpy(gretl_errmsg, _("error in new starting obs"));
		return 1;
	    }
	    pdinfo->t1 = new_t1;
	    return 0;
	}
    }

    /* now we're looking at nf = 2 (2 fields) case */

    if (sscanf(line, "%10s %10s", newstart, newstop) != 2) {
	strcpy(gretl_errmsg, _("error reading smpl line"));
	return 1;
    }

    if (strcmp(newstart, ";")) {
	new_t1 = get_sample_limit(newstart, pZ, pdinfo, SMPL_T1);
	if (new_t1 < 0 || new_t1 >= pdinfo->n) {
	    strcpy(gretl_errmsg, _("error in new starting obs"));
	    return 1;
	}	
    }

    if (strcmp(newstop, ";")) {
	new_t2 = get_sample_limit(newstop, pZ, pdinfo, SMPL_T2);
	if (new_t2 < 0 || new_t2 >= pdinfo->n) {
	    strcpy(gretl_errmsg, _("error in new ending obs"));
	    return 1;
	}
    }

    if (new_t1 < 0 || new_t1 > new_t2) {
	strcpy(gretl_errmsg, _("Invalid null sample"));
	return 1;
    }

    pdinfo->t1 = new_t1;
    pdinfo->t2 = new_t2;

    return 0;
}

int count_missing_values (double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int i, t, tmiss;
    int missval = 0, missobs = 0, totvals = 0, oldmiss = 0;
    int *missvec;

    missvec = malloc(pdinfo->v * sizeof missvec);
    if (missvec != NULL) {
	for (i=0; i<pdinfo->v; i++) missvec[i] = 0;
    }

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	tmiss = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (var_is_hidden(pdinfo, i) || var_is_scalar(pdinfo, i)) {
		continue;
	    }
	    if (na((*pZ)[i][t])) {
		if (missvec[i] == 0) {
		    missvec[0] += 1;
		}
		missvec[i] += 1;
		missval++;
	    }
	    totvals++;
	}

	tmiss = missval - oldmiss;

	if (tmiss) {
	    missobs++;

	    if (pdinfo->markers) { 
		pprintf(prn, "%8s %4d %s\n", pdinfo->S[t], tmiss,
			_("missing values"));
	    } else {
		char tmp[OBSLEN];

		ntodate(tmp, t, pdinfo);
		pprintf(prn, "%8s %4d %s\n", tmp, tmiss,
			_("missing values"));
	    }
	}
	oldmiss = missval;
    }

    pprintf(prn, _("\nNumber of observations (rows) with missing data "
	    "values = %d (%.2f%%)\n"), missobs, 
	    (100.0 * missobs / (pdinfo->t2 - pdinfo->t1 + 1)));
    pprintf(prn, _("Total number of missing data values = %d (%.2f%% "
	    "of total data values)\n"), missval, 
	    (100.0 * missval / totvals));
    if (missvec[0] > 0) {
	pputc(prn, '\n');
	for (i=1; i<pdinfo->v; i++) {
	    if (missvec[i] > 0) {
		pprintf(prn, "%8s: %d %s\n", pdinfo->varname[i], 
			missvec[i], _("missing values"));
	    }
	}
    }

    free(missvec);

    return missval;
}

static void copy_series_info (DATAINFO *dest, const DATAINFO *src)
{
    int i;

    for (i=1; i<src->v; i++) {
	strcpy(dest->varname[i], src->varname[i]);
	if (src->varinfo != NULL) {
	    copy_varinfo(dest->varinfo[i], src->varinfo[i]);
	}
    }
}

/* Below: For a model that was estimated using a data sample which
   differs from the current sample, reconstitute the dataset
   appropriate to the model.  If the model used a subsample, we use
   the special dummy variable (submask) saved with the model to select
   observations from the full dataset.  Another possibility is that
   the model used the full dataset, which has by now been subsampled.

   We attach this dataset to the model, so that it can be used for,
   e.g., residual or fitted versus actual plots.  The attached data
   set will be destroyed when the model itself is destroyed.
*/

int add_dataset_to_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    int sn = 0;
    double **modZ = NULL;
    DATAINFO *modinfo = NULL;
    char *mask = NULL;

    if (fullZ == NULL || fullinfo == NULL) {
#if SUBDEBUG
	fputs("add_subsampled_dataset_to_model: got NULL full dataset\n",
	      stderr);
#endif
	return 1;
    }

    if (pmod->dataset != NULL) {
#if SUBDEBUG
	fputs("add_subsampled_dataset_to_model: job already done\n",
	      stderr);
#endif	
	return 0;
    }

    /* sync malloced elements that may have moved */
    sync_dataset_elements(pdinfo);

    if (pmod->submask == NULL) {
	/* no subsample info: model was estimated on full dataset,
	   so we reconstruct the full dataset */
	sn = fullinfo->n;
#if SUBDEBUG
	fprintf(stderr, "pmod->submask = NULL, set sn = %d\n", sn);
#endif
    } else {
	/* model estimated on subsample, which has to be reconstructed */
	int t;

	mask = calloc(fullinfo->n, 1);
	if (mask == NULL) {
	    return E_ALLOC;
	}
#if SUBDEBUG
	fprintf(stderr, "setting mask from pmod->submask: fullinfo->n = %d\n", 
		fullinfo->n);
#endif
	for (t=0; t<fullinfo->n; t++) {
	    if (pmod->submask[t] > 0) {
		mask[t] = 1;
		sn++;
	    }
	}
	if (sn == 0) {
	    free(mask);
	    return 1;
	}
    }

    pmod->dataset = malloc(sizeof *pmod->dataset);
    if (pmod->dataset == NULL) {
	free(mask);
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "pmod->dataset allocated at %p\n", 
	    (void *) pmod->dataset);
#endif

    /* initial allocation of sub-sampled dataset */
    modinfo = create_new_dataset(&modZ, fullinfo->v, sn,
				 fullinfo->markers);
    if (modinfo == NULL) {
	free(mask);
	free(pmod->dataset);
	pmod->dataset = NULL;
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "dataset created, copying series info\n");
#endif

    /* copy across info on series */
    copy_series_info(modinfo, fullinfo);

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(modZ, modinfo,
			   (const double **) fullZ, fullinfo,
			   mask);

#if SUBDEBUG
    fputs("data copied to subsampled dataset\n", stderr);
#endif

    /* dataset characteristics such as pd: if we're rebuilding the
       _full_ dataset copy these across; but if we're reconstructing a
       subsampled dataset these features from the full dataset will be
       wrong, and we stay with the simple defaults
    */
    if (pmod->submask == NULL) {
	modinfo->pd = fullinfo->pd;
	modinfo->sd0 = fullinfo->sd0;
	strcpy(modinfo->stobs, fullinfo->stobs);
	strcpy(modinfo->endobs, fullinfo->endobs);
	modinfo->structure = fullinfo->structure;
    }

    pmod->dataset->Z = modZ;
    pmod->dataset->dinfo = modinfo;

    free(mask);

#if SUBDEBUG
    fputs("add_subsampled_dataset_to_model: success\n", stderr);
#endif

    return 0;
}

void free_model_dataset (MODEL *pmod)
{
    if (pmod->dataset != NULL) {
#if SUBDEBUG
	fprintf(stderr, "Deep freeing model->dataset\n");
#endif
	destroy_dataset(pmod->dataset->Z, pmod->dataset->dinfo);
	free(pmod->dataset);
	pmod->dataset = NULL;
    }
}
