/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* subsample.c for gretl */

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"

#define SUBDEBUG 1

/*
  The purpose of the static pointers below: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting cases
  rather than just moving the starting or ending points of the data
  range -- we create a new sub-dataset, and we need to keep the full
  dataset around so that it can be restored later.  The pointers fullZ
  and fullinfo are used to record the addresses of the full data
  matrix and DATAINFO struct respectively.

  Another issue arises: if the user replaces or clears a dataset while
  it is subsampled, we want to free the associated full dataset also.
  The "peerinfo" pointer is used to facilitate this.  On subsampling,
  when fullZ and fullinfo are assigned to, peerinfo is pointed at the
  associated subsampled DATAINFO struct.  Then, on freeing the
  subsampled dataset, we check whether its DATAINFO address matches
  peerinfo: if so, we free up fullZ and fullinfo.

*/

static double **fullZ;
static DATAINFO *fullinfo;
static const DATAINFO *peerinfo;

static int parent_subdum_index (const DATAINFO *pdinfo);

#define SUBMASK_SENTINEL 127

static int get_submask_length (const char *s)
{
    int n = 1;

    while (*s != SUBMASK_SENTINEL) {
	n++;
	s++;
    }

    return n;
}

char *copy_subsample_mask (const char *src)
{
    int n = get_submask_length(src);
    char *ret = NULL;

    if (src != NULL) {
	ret = malloc(n * sizeof *ret);
	if (ret != NULL) {
	    memcpy(ret, src, n);
	}
    }

    return ret;
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

static void relink_full_dataset (double ***pZ, DATAINFO **ppdinfo)
{
    *pZ = fullZ;
    *ppdinfo = fullinfo;

    fullZ = NULL;
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
    fullinfo->vector = pdinfo->vector;
    fullinfo->descrip = pdinfo->descrip;      
}

static int 
attach_subsample_to_datainfo (DATAINFO *subinfo, const double **Z, 
			      const DATAINFO *pdinfo)
{
    int s, t;

    subinfo->submask = malloc((pdinfo->n + 1) * sizeof *subinfo->submask);
    if (subinfo->submask == NULL) {
	return E_ALLOC;
    }

    s = varindex(pdinfo, "subdum");

    if (s == pdinfo->v) {
	s = parent_subdum_index(pdinfo);
    }	

    if (s == pdinfo->v) { 
	/* safety measure: should be impossible */
	fprintf(stderr, "mystery failure in attach_subsample_to_dataset\n");
	return 1;   
    } 

    for (t=0; t<pdinfo->n; t++) {
	subinfo->submask[t] = Z[s][t];
    }

    subinfo->submask[pdinfo->n] = SUBMASK_SENTINEL;

    return 0;
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

    if (fullZ != NULL) {
	/* sync, in case anything has moved */
	sync_dataset_elements(pdinfo);

	pmod->submask = copy_subsample_mask(pdinfo->submask);
	if (pmod->submask == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* apparatus for updating full dataset when restoring full sample
   after sub-sampling */

static void
update_full_data_values (const double **subZ, const DATAINFO *pdinfo)
{
    int i, t;

    for (i=1; i<fullinfo->v; i++) {
	int subt = 0;

	if (!pdinfo->vector[i]) {
	    fullZ[i][0] = subZ[i][0];
	} else {
	    for (t=0; t<fullinfo->n; t++) {
		if (pdinfo->submask[t]) {
		    fullZ[i][t] = subZ[i][subt++];
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
    int i, t, subt;
    int newvars = pdinfo->v - fullinfo->v;
    int n = fullinfo->n;
    double **newZ = NULL;
    int err = 0;

    if (newvars <= 0) {
	return 0;
    }

    if (pdinfo->submask == NULL) {
	return E_NOMERGE;
    }

    /* allocate expanded data array */
    newZ = realloc(fullZ, pdinfo->v * sizeof *fullZ);
    if (newZ != NULL) {
	for (i=0; i<newvars; i++) {
	    if (pdinfo->vector[fullinfo->v+i]) {
		newZ[fullinfo->v+i] = malloc(n * sizeof **newZ);
	    } else {
		newZ[fullinfo->v+i] = malloc(sizeof **newZ);
	    }
	    if (newZ[fullinfo->v+i] == NULL) {
		err = 1;
		break;
	    }
	}
    } else {
	err = 1;
    }

    if (err) {
	return E_ALLOC;
    }
    
    fullZ = newZ;

    for (i=fullinfo->v; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) {
	   fullZ[i][0] = Z[i][0]; 
	}
    }

    subt = 0;

    for (t=0; t<n; t++) {
	if (pdinfo->submask[t]) {
	    for (i=fullinfo->v; i<pdinfo->v; i++) {
		if (pdinfo->vector[i]) {
		    fullZ[i][t] = Z[i][subt];
		}
	    }
	    subt++;
	} else {
	    for (i=fullinfo->v; i<pdinfo->v; i++) { 
		if (pdinfo->vector[i]) {
		    fullZ[i][t] = NADBL;
		}
	    }
	}
    }

    fullinfo->v = pdinfo->v;

    return 0;
}

static int parent_subdum_index (const DATAINFO *pdinfo)
{
    int i, d = gretl_function_stack_depth();
    int v = pdinfo->v;

    fprintf(stderr, "parent_subdum_index: fn depth = %d\n", d);

    if (d > 0) {
	for (i=1; i<pdinfo->v; i++) {
	    if (STACK_LEVEL(pdinfo, i) == d - 1 &&
		!strcmp(pdinfo->varname[i], "subdum")) {
		v = i;
		break;
	    }
	}
    }

    return v;
}

static int maybe_add_subdum (double ***pZ, DATAINFO *pdinfo, int *snum)
{
    int v = varindex(pdinfo, "subdum");

    if (v == pdinfo->v) {
	v = parent_subdum_index(pdinfo);
    }

    if (v == pdinfo->v) {
	/* variable doesn't exist: create it */
	if (dataset_add_series(1, pZ, pdinfo)) {
	    return 1;
	}
	strcpy(pdinfo->varname[v], "subdum");
	strcpy(VARLABEL(pdinfo, v), _("automatic sub-sampling dummy"));
    } 

    *snum = v;

    return 0;
}

static int make_smpl_mask (double ***pZ, DATAINFO *pdinfo)
{
    int old_vmax = pdinfo->v;
    int v, t;
    int err = 0;

    /* add "subdum" variable if it's not already present */
    if (maybe_add_subdum(pZ, pdinfo, &v)) {
	err = 1;
    } else if (v == old_vmax) {
	/* no pre-existing mask: write it from scratch */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[v][t] = (t < pdinfo->t1 || t > pdinfo->t2)? 0.0 : 1.0;
	}
    } else {
	/* there was a pre-existing mask: adjust it */
	for (t=0; t<pdinfo->n; t++) {
	    if (t < pdinfo->t1 || t > pdinfo->t2) {
		(*pZ)[v][t] = 0.0;
	    }
	}	
    }

    return err;
}

/* If the third parameter, "cumulate", is non-zero, we're restoring
   the full data set just as a preliminary to applying a sample
   restriction that will be the logical product of some new
   restriction and any existing restriction.
*/

static int 
real_restore_full_sample (double ***pZ, DATAINFO **ppdinfo, int cumulate)
{
    int i, t, err = 0;

    *gretl_errmsg = '\0';

#if SUBDEBUG
    fprintf(stderr, "real_restore_full_sample: pZ=%p, ppdinfo=%p cum=%d\n",
	    (void *) pZ, (void *) ppdinfo, cumulate);
    fprintf(stderr, "*pZ=%p, *ppdinfo=%p\n",
	    (void *) *pZ, (void *) *ppdinfo);
#endif

    if (*pZ != NULL && cumulate) {  
	/* cumulating, not replacing restrictions: we need to
	   ensure that an appropriate mask is in place */
	err = make_smpl_mask(pZ, *ppdinfo);
	if (err) {
	    return err;
	}
    }

    /* simple case: merely resetting the sample starting
       and ending points to full range */
    if (!complex_subsampled()) {
	(*ppdinfo)->t1 = 0;
	(*ppdinfo)->t2 = (*ppdinfo)->n - 1;
#if SUBDEBUG
	fprintf(stderr, " reset t1 and t2; all done\n");
#endif
	return 0;
    }

    /* beyond this point we are doing a non-trivial restoration
       of a stored full dataset, which has previously been
       subsampled, e.g., by some boolean criterion 
    */

    /* reattach malloc'd elements of subsampled dataset,
       which may have moved */
    sync_dataset_elements(*ppdinfo);

    /* update values for pre-existing series, which may have been
       modified via genr */
    update_full_data_values((const double **) *pZ, *ppdinfo);

    /* if case markers were added when subsampled, carry them back */
    update_case_markers(*ppdinfo);

    /* in case any new vars were added when subsampled, try to merge
       them into the full dataset */
    err = add_new_vars_to_full((const double **) *pZ, *ppdinfo);
    if (err == E_ALLOC) {
        sprintf(gretl_errmsg, _("Out of memory expanding data set\n"));
    } else if (err == E_NOMERGE) {
        sprintf(gretl_errmsg, 
		_("Missing sub-sample information; can't merge data\n"));
    }

    /* If we plan to destroy/replace any existing sample restrictions,
       zero out the "subdum" dummy variable. */
    if (!cumulate) {
	i = varindex(fullinfo, "subdum");
	if (i < fullinfo->v) {
	    for (t=0; t<fullinfo->n; t++) {
		fullZ[i][t] = 0.0;
	    }
	}
    }

    free_Z(*pZ, *ppdinfo);
    clear_datainfo(*ppdinfo, CLEAR_SUBSAMPLE);
    free(*ppdinfo);

    relink_full_dataset(pZ, ppdinfo);

    return 0;
}

/* restore_full_sample: 
 * @pZ: pointer to data array.  
 * @ppdinfo: address of data info pointer. 
 *
 * Restore the full data range, undoing any sub-sampling that was
 * previously performed (either by shifting the endpoints of the
 * sample range or by selecting cases on some criterion or other).
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restore_full_sample (double ***pZ, DATAINFO **ppdinfo)
{
    return real_restore_full_sample(pZ, ppdinfo, 0);
}

static double *get_old_mask (double **Z, const DATAINFO *pdinfo)
{
    double *mask = NULL;
    int v = varindex(pdinfo, "subdum");

    fprintf(stderr, "get_old_mask: v = %d\n", v);

    if (v == pdinfo->v) {
	fprintf(stderr, "get_old_mask: trying parent_subdum_index\n");
	v = parent_subdum_index(pdinfo);
    }

    if (v < pdinfo->v && gretl_isdummy(0, pdinfo->n - 1, Z[v])) {
	mask = Z[v];
    }    

    return mask;
}

static int overlay_masks (char *targ, const double *src, int n)
{
    int i, sn = 0;
    
    for (i=0; i<n; i++) {
	if (targ[i] == 1 && src[i] == 1.0) {
	    targ[i] = 1;
	    sn++;
	} else {
	    targ[i] = 0;
	}
    }
	
    return sn;
}

static int 
make_missing_mask (const double **Z, const DATAINFO *pdinfo,
		   const int *list, char *mask)
{
    int i, t, v, sn = 0;
    int lmax;

    if (list != NULL) {
	/* examine only a specified list of variables */
	lmax = list[0];
    } else {
	/* default: examine all variables */
	lmax = pdinfo->v - 1;
    }

    for (t=0; t<pdinfo->n; t++) {
	mask[t] = 1;
	for (i=1; i<=lmax; i++) {
	    if (list != NULL) {
		v = list[i];
	    } else {
		v = i;
	    }
	    if (pdinfo->vector[v] && na(Z[v][t])) {
		mask[t] = 0;
		break;
	    }
	}
	if (mask[t] == 1) sn++;
    }

    return sn;
} 

static void copy_to_mask (char *mask, const double *x, int n)
{
    int t;

    for (t=0; t<n; t++) {
	mask[t] = (x[t] != 0.0 && !na(x[t]));
    }
}

static int sn_from_tmp_dummy (double ***pZ, DATAINFO *pdinfo, 
			      const char *line, char *mask)
{
    char formula[MAXLEN];
    int dnum, isdum, dlist[2];
    int err;

    /* + 4 to skip the command word "smpl" */
    sprintf(formula, "__tmpmsk=%s", line + 4);

    err = generate(formula, pZ, pdinfo, OPT_P, NULL);

    if (err) {
	return -1;
    }

    dnum = varindex(pdinfo, "tmpmsk");
    isdum = gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[dnum]);
    copy_to_mask(mask, (*pZ)[dnum], pdinfo->n);

    dlist[0] = 1;
    dlist[1] = dnum;
    dataset_drop_listed_variables(dlist, pZ, pdinfo, NULL);

    return isdum;
}

static int sn_from_dummy (const double **Z, const DATAINFO *pdinfo,
			  const char *dname, char *mask)
{
    int dnum = varindex(pdinfo, dname);
    int isdum;

    if (dnum == pdinfo->v) {
	sprintf(gretl_errmsg, _("Variable '%s' not defined"), dname);
	return -1;
    } 

    isdum = gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[dnum]);

    copy_to_mask(mask, Z[dnum], pdinfo->n);

    return isdum;
} 

static int count_selected_cases (const char *x, int n)
{
    int i, count = 0;

    for (i=0; i<n; i++) {
	if (x[i] != 0) {
	    count++;
	}
    }

    return count;
}

static int make_random_mask (double *oldmask, int fulln, int subn, 
			     char *mask)
{
    int i, cases = 0, err = 0;
    unsigned u;

    if (subn <= 0 || subn >= fulln) {
	err = 1;
    } else if (oldmask != NULL) {
	int oldn = 0;

	for (i=0; i<fulln; i++) {
	    if (oldmask[i] != 0.0) oldn++;
	}

	if (subn >= oldn) err = 1;
    }	

    if (err) {
	sprintf(gretl_errmsg, _("Invalid number of cases %d"), subn);
	return 0;
    }	

    for (i=0; i<fulln; i++) {
	mask[i] = 0;
    }

    for (i=0; (cases != subn); i++) {
	u = gretl_rand_int_max(fulln);
	if (oldmask == NULL || oldmask[u] == 1.0) {
	    mask[u] = 1;
	}
	if (i >= subn - 1) {
	    cases = count_selected_cases(mask, fulln);
	}
    }

    return cases;
}

enum {
    SUBSAMPLE_UNKNOWN,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN,
    SUBSAMPLE_RANDOM
} subsample_options;

static void backup_full_dataset (double ***pZ, DATAINFO **ppdinfo,
				 const DATAINFO *newinfo)
{
    fullZ = *pZ;
    fullinfo = *ppdinfo;
    peerinfo = newinfo;
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

static int block_retained (const char *s)
{
    int nrows = 0;
    int ret = 0;

#if SUBDEBUG
    fprintf(stderr, "block_retained: s = '%s'\n", s);
#endif

    while (*s) {
	if (*s == '1') {
	    nrows++;
	}
	if (nrows > 1) {
	    ret = 1;
	    break;
	}
	s++;
    }

    return ret;
}

static char *make_panel_signature (const DATAINFO *pdinfo,
				   const char *mask,
				   int mt1, int mt2,
				   int neg)
{
    char *sig;
    size_t sz = pdinfo->n + pdinfo->n / pdinfo->pd;
    int t;

    sig = calloc(sz, 1);

#if SUBDEBUG
    for (t=0; t<pdinfo->n; t++) {
	if (t >= mt1 && t <= mt2) {
	    fprintf(stderr, "mask[%d]=%c\n", t-mt1, mask[t-mt1]);
	}
    }
#endif

    if (sig != NULL) {
	int i = 0;

	for (t=0; t<pdinfo->n; t++) {
	    if (t > 0 && t % pdinfo->pd == 0) {
		i++;
	    } 
	    if (t < mt1 || t > mt2) {
		sig[i++] = '0';
	    } else {
		if (neg) {
		    /* using model missmask: '1' for missing vs '0' */
		    sig[i++] = (mask[t - mt1] != '1')? '1' : '0';
		} else {
		    /* using sample missmask: 0 for missing vs 1 */
		    sig[i++] = (mask[t - mt1] != 0)? '1' : '0';
		}
	    }
	}
    }

    return sig;
}

static int 
real_mask_leaves_balanced_panel (const char *mask,
				 const DATAINFO *pdinfo,
				 int *new_blocks,
				 int mt1, int mt2, int neg)
{
    int i;
    int ok_blocks = 0, unbal = 0;
    int nblocks = pdinfo->n / pdinfo->pd;
    char *sig = make_panel_signature(pdinfo, mask, mt1, mt2, neg);
    char *sig_0 = NULL;

    /* FIXME: allow for possibility that starting obs is not 1:1 */

    if (sig == NULL) {
	return 0;
    }

    for (i=0; i<nblocks; i++) {
	char *sig_n = sig + i * (pdinfo->pd + 1);

#if SUBDEBUG
	fprintf(stderr, "block %d: sig_n = '%s'\n", i, sig_n);
#endif

	if (block_retained(sig_n)) {
	    if (ok_blocks == 0) {
		sig_0 = sig_n;
	    } else {
		if (strcmp(sig_n, sig_0)) {
		    unbal = 1;
		}
	    }
	    if (!unbal) {
		ok_blocks++;
	    }
	} 
	if (unbal) {
	    break;
	}
    }

    free(sig);

    if (new_blocks != NULL) {
	*new_blocks = ok_blocks;
    }

    return (ok_blocks > 0 && unbal == 0);
}

int model_mask_leaves_balanced_panel (const MODEL *pmod,
				      const DATAINFO *pdinfo)
{
    return real_mask_leaves_balanced_panel(pmod->missmask, 
					   pdinfo, NULL,
					   pmod->t1, pmod->t2,
					   1);
}

/* when subsampling cases from a time series, if the resulting dataset
   is still a time series without internal "holes" then automatically
   set a time series interpretation of the reduced dataset
*/

static void
maybe_reconstitute_time_series (const DATAINFO *pdinfo,
				const char *mask,
				DATAINFO *subinfo)
{
    int missing = 1, switches = 0;
    int t, t1 = 0, ts_ok = 1;

    for (t=0; t<pdinfo->n; t++) {
	if (missing && mask[t] == 1) {
	    t1 = t;
	    switches++;
	    missing = 0;
	} else if (!missing && mask[t] == 0) {
	    switches++;
	    missing = 1;
	}
	if (switches > 2) {
	    ts_ok = 0;
	    break;
	}
    }

    if (ts_ok) {
	char line[32];
	char stobs[OBSLEN];

	ntodate(stobs, t1, pdinfo);
	sprintf(line, "setobs %d %s", pdinfo->pd, stobs);
	set_obs(line, subinfo, OPT_NONE);
    } 
}

/* when subsampling cases from a panel dataset, if the resulting
   dataset is still a balanced panel then automatically set the
   interpretation of the reduced dataset as a panel
*/

static void
maybe_reconstitute_panel (const DATAINFO *pdinfo,
			  const char *mask,
			  DATAINFO *subinfo)
{
    int bal, ok_blocks = 0;

    bal = real_mask_leaves_balanced_panel(mask, pdinfo, &ok_blocks,
					  0, pdinfo->n - 1, 0);

    if (bal && ok_blocks > 1) {
	char line[32];
	int pd = subinfo->n / ok_blocks;
	gretlopt opt = OPT_C;

	sprintf(line, "setobs %d 1.1", pd);
	if (pdinfo->structure == STACKED_TIME_SERIES) {
	    opt = OPT_S;
	} 
	set_obs(line, subinfo, opt);
    } 
}

static void 
copy_data_to_subsample (double **subZ, DATAINFO *subinfo,
			const double **Z, const DATAINFO *pdinfo,
			const char *mask)
{
    int i, t, st;

    /* copy data values */
    for (i=1; i<pdinfo->v; i++) {
	if (pdinfo->vector[i]) {
	    st = 0;
	    for (t=0; t<pdinfo->n; t++) {
		if (mask == NULL || mask[t] == 1) {
		    subZ[i][st++] = Z[i][t];
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
	   if (mask == NULL || mask[t] == 1) {
	       strcpy(subinfo->S[st++], pdinfo->S[t]);
	   }
       }
   }

    strcpy(subinfo->stobs, "1");
    sprintf(subinfo->endobs, "%d", subinfo->n);
}

/* restrict_sample: 
 * @line: command line (or %NULL).  
 * @pZ: pointer to original data array.  
 * @ppdinfo: address of original data info pointer. 
 * @list: list of variables in case of OPT_M (or %NULL).  
 * @opt: option flags.
 * @prn: printing apparatus.
 *
 * Sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values (OPT_M); or using as a mask a
 * specified dummy variable (OPT_O); or masking with a specified
 * boolean condition (OPT_T); or selecting at random (OPT_N).
 *
 * In case OPT_M a @list of variables may be supplied; in cases
 * OPT_O, OPT_T and OPT_N, @line must contain specifics.
 *
 * In case OPT_R is included, the restriction will Replace any
 * existing sample restriction, otherwise the resulting restriction
 * will be the logical product of the new restriction and any
 * existing restriction.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample (const char *line, 
		     double ***pZ, DATAINFO **ppdinfo, 
		     const int *list, gretlopt opt,
		     PRN *prn)
{
    int mode = SUBSAMPLE_UNKNOWN;
    char dname[VNAMELEN] = {0};
    double *oldmask = NULL;
    double **subZ = NULL;
    DATAINFO *subinfo = NULL;
    char *mask = NULL;
    int cumulate = 1;
    int subnum = 0;
    int t, sn = 0;
    int err = 0;

    *gretl_errmsg = '\0';

    if (opt & OPT_M) {
	mode = SUBSAMPLE_DROP_MISSING;
    } else if (opt & OPT_T) {
	mode = SUBSAMPLE_BOOLEAN;
    } else if (opt & OPT_N) {
	mode = SUBSAMPLE_RANDOM;
    } else if (opt & OPT_O) {
	if (line != NULL && sscanf(line, "%*s %s", dname)) {
	    mode = SUBSAMPLE_USE_DUMMY;
	} else {
	    mode = SUBSAMPLE_DROP_MISSING;
	}
    } else {
	strcpy(gretl_errmsg, "Unrecognized sample command");
	return 1;
    }

    if (opt & OPT_R) {
	/* Replace any existing restrictions, don't cumulate. */
	cumulate = 0;
    }

    /* The first step is to restore the full data range for
       housekeeping purposes. */
    err = real_restore_full_sample(pZ, ppdinfo, cumulate);
    if (err) {
	return err;
    }

    /* allocate temporary (char) mask series */
    mask = calloc((*ppdinfo)->n, 1);
    if (mask == NULL) {
	return E_ALLOC;
    }

    /* if there's no currently existing sub-sample mask, 
       oldmask will be NULL (which is OK) */
    oldmask = get_old_mask(*pZ, *ppdinfo);
#if SUBDEBUG
    fprintf(stderr, "restrict_sample: oldmask = %p\n", (void *) oldmask);
#endif

    /* construct subsample mask in one of several possible ways */

    if (mode == SUBSAMPLE_DROP_MISSING) {   
	sn = make_missing_mask((const double **) *pZ, *ppdinfo, list, mask);
    } else if (mode == SUBSAMPLE_RANDOM) {
	sn = make_random_mask(oldmask, (*ppdinfo)->n, atoi(line + 4), mask);
	if (sn == 0) {
	    err = 1;
	}
    } else if (mode == SUBSAMPLE_USE_DUMMY || mode == SUBSAMPLE_BOOLEAN) {
	if (mode == SUBSAMPLE_USE_DUMMY) {
	    sn = sn_from_dummy((const double **) *pZ, *ppdinfo, dname, mask);
	} else {
	    sn = sn_from_tmp_dummy(pZ, *ppdinfo, line, mask);
#if SUBDEBUG
	    fprintf(stderr, "SUBSAMPLE_BOOLEAN: sn = %d\n", sn);
#endif
	}
	if (sn < 0) {
	    err = 1;
	}
    } else {
	/* impossible */
	strcpy(gretl_errmsg, _("Sub-sample command failed mysteriously"));
	err = 1;
    }

    /* exit now on unrecoverable error */
    if (err) {
	free(mask);
	return 1;
    }

    /* cumulate sample restrictions, if appropriate */
    if (oldmask != NULL && mode != SUBSAMPLE_RANDOM) {
	sn = overlay_masks(mask, oldmask, (*ppdinfo)->n);
    }

    /* does this policy lead to an empty sample, or no change
       in the sample, perchance? */
    if (sn == 0) {
	if (mode == SUBSAMPLE_USE_DUMMY) {
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), dname);
	} else if (mode == SUBSAMPLE_DROP_MISSING) {
	    strcpy(gretl_errmsg, _("No observations would be left!"));
	} else if (mode == SUBSAMPLE_BOOLEAN) { 
	    if (mask[(*ppdinfo)->t1] == 0) {
		strcpy(gretl_errmsg, _("No observations would be left!"));
	    } else {
		/* this is not really an error, just a no-op */
		if (gretl_messages_on()) {
		    pputs(prn, _("No observations were dropped!"));
		    pputc(prn, '\n');
		}
		free(mask);
		return 0;
	    }
	}
	err = 1;
    } else if (sn == (*ppdinfo)->n) {
	strcpy(gretl_errmsg, _("No observations were dropped!"));
	err = 1;
    }

    /* exit now on nonsense subsample */
    if (err) {
	free(mask);
	return 1;
    }

    /* now we start making the new dataset */

    /* allocate new datainfo */
    subinfo = datainfo_new();
    if (subinfo == NULL) {
	free(mask);
	return E_ALLOC;
    }

    /* create "hidden" dummy in full dataset to record sub-sample, 
       if need be */
    if (maybe_add_subdum(pZ, *ppdinfo, &subnum)) {
	free(subinfo);
	free(mask);
	return E_ALLOC;
    } 

    /* write the new mask into the "subdum" variable */
    for (t=0; t<(*ppdinfo)->n; t++) {
	(*pZ)[subnum][t] = mask[t];
    }

    /* set up the sub-sampled datainfo */
    subinfo->n = sn;
    subinfo->v = (*ppdinfo)->v;
    if (start_new_Z(&subZ, subinfo, 1)) {
	free(subinfo);
	free(mask);
	return E_ALLOC;
    }

    /* link (don't copy) varnames and descriptions, 
       these are not dependent on the series length */
    subinfo->varname = (*ppdinfo)->varname;
    subinfo->varinfo = (*ppdinfo)->varinfo;
    subinfo->descrip = (*ppdinfo)->descrip;
    subinfo->vector = (*ppdinfo)->vector;

    /* set up case markers? */
    if ((*ppdinfo)->markers) {
	err = dataset_allocate_obs_markers(subinfo);
	if (err) {
	    free_Z(subZ, subinfo);
	    free(subinfo);
	    free(mask);
	    return E_ALLOC;
	}
	subinfo->markers = (*ppdinfo)->markers;
    }

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(subZ, subinfo, 
			   (const double **) *pZ, *ppdinfo,
			   mask);

    attach_subsample_to_datainfo(subinfo, 
				 (const double **) *pZ, 
				 *ppdinfo);

    if (mode == SUBSAMPLE_USE_DUMMY || 
	mode == SUBSAMPLE_BOOLEAN ||
	mode == SUBSAMPLE_DROP_MISSING) {
	if (dataset_is_panel(*ppdinfo)) {
	    maybe_reconstitute_panel(*ppdinfo, mask, subinfo);
	} else if (dataset_is_time_series(*ppdinfo)) {
	    maybe_reconstitute_time_series(*ppdinfo, mask, subinfo);
	}
    }

    /* save state */
    backup_full_dataset(pZ, ppdinfo, subinfo);

    /* and switch pointers */
    *pZ = subZ;
    *ppdinfo = subinfo;

    free(mask);

    return 0;
}

enum {
    SMPL_T1,
    SMPL_T2
};

static int 
get_sample_limit (char *s, const double **Z, DATAINFO *pdinfo,
		  int code)
{
    int v, ret = -1;

    if (*s == '-' || *s == '+') {
	/* incremental form */
	int incr = 0;

	if (isdigit((unsigned char) s[1])) {
	    incr = atoi(s);
	} else {
	    v = varindex(pdinfo, s + 1);
	    if (v < pdinfo->v) {
		incr = (int) Z[v][0];
		if (*s == '-') {
		    incr = -incr;
		}
	    }
	}
	if (code == SMPL_T1) {
	    ret = pdinfo->t1 + incr;
	} else {
	    ret = pdinfo->t2 + incr;
	}
    } else {
	/* absolute form */
	ret = get_t_from_obs_string(s, Z, pdinfo);
    }

    return ret;
}

int set_sample (const char *line, const double **Z, DATAINFO *pdinfo)
{
    int nf, new_t1 = pdinfo->t1, new_t2 = pdinfo->t2;
    char cmd[5], newstart[OBSLEN], newstop[OBSLEN];

    *gretl_errmsg = '\0';

    nf = count_fields(line);

#if SUBDEBUG
    fprintf(stderr, "set_sample: line='%s', nf=%d, pdinfo=%p\n", 
	    line, nf, (void *) pdinfo);
    if (pdinfo != NULL) {
	fprintf(stderr, "pdinfo->v = %d, pdinfo->n = %d\n",
		pdinfo->v, pdinfo->n);
    }
#endif

    if (nf == 3 && pdinfo->n == 0) {
	/* database special */
	return db_set_sample(line, pdinfo);
    }

    if (nf == 1) return 0;
	
    if (nf == 2) {
	if (sscanf(line, "%4s %10s", cmd, newstart) != 2) {
	    strcpy(gretl_errmsg, _("error reading smpl line"));
	    return 1;
	} else {
	    new_t1 = get_sample_limit(newstart, Z, pdinfo, SMPL_T1);
	    if (new_t1 < 0 || new_t1 >= pdinfo->n) {
		strcpy(gretl_errmsg, _("error in new starting obs"));
		return 1;
	    }
	    pdinfo->t1 = new_t1;
	    return 0;
	}
    }

    /* now we're looking at nf = 3 (3 fields) case */

    if (sscanf(line, "%4s %10s %10s", cmd, newstart, newstop) != 3) {
	strcpy(gretl_errmsg, _("error reading smpl line"));
	return 1;
    }

    if (strcmp(newstart, ";")) {
	new_t1 = get_sample_limit(newstart, Z, pdinfo, SMPL_T1);
	if (new_t1 < 0 || new_t1 >= pdinfo->n) {
	    strcpy(gretl_errmsg, _("error in new starting obs"));
	    return 1;
	}	
    }

    if (strcmp(newstop, ";")) {
	new_t2 = get_sample_limit(newstop, Z, pdinfo, SMPL_T2);
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
	    if (is_hidden_variable(i, pdinfo) || !pdinfo->vector[i]) {
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

static void copy_varinfo (VARINFO *dest, const VARINFO *src)
{
    strcpy(dest->label, src->label);
    strcpy(dest->display_name, src->display_name);
    dest->compact_method = src->compact_method;
}

static void copy_series_info (DATAINFO *dest, const DATAINFO *src)
{
    int i;

    for (i=1; i<src->v; i++) {
	strcpy(dest->varname[i], src->varname[i]);
	if (src->varinfo != NULL) {
	    copy_varinfo(dest->varinfo[i], src->varinfo[i]);
	}
	dest->vector[i] = src->vector[i];	
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
