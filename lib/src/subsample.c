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
#include "gretl_private.h"

#undef SUBDEBUG

/* private to libgretl */

/*
  The purpose of these static pointers: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting
  cases rather than just moving the starting or ending points of
  the data range -- we create a new sub-dataset, and we need to
  keep the full dataset around so that it can be restored later.
  The pointers fullZ and fullinfo are used to record the addresses
  of the full data matrix and DATAINFO struct respectively.

  Another issue arises: if the user replaces or clears a dataset
  while it is subsampled, we want to free the associated full
  dataset also.  The peerinfo pointer is used to facilitate
  this.  On subsampling, when fullZ and fullinfo are assigned
  to, peerinfo is pointed at the associated subsampled
  DATAINFO struct.  Then, on freeing the subsampled dataset,
  we check whether its DATAINFO address matches peerinfo: if
  so, we free up fullZ and fullinfo.
*/

static double **fullZ;
static DATAINFO *fullinfo;
static DATAINFO *peerinfo;

/* .......................................................... */

char *copy_subdum (const char *src, int n)
{
    char *ret;

    if (n == 0 || src == NULL) return NULL;

    ret = malloc(n * sizeof *ret);
    if (ret == NULL) return NULL;

    memcpy(ret, src, n);

    return ret;
}

/* .......................................................... */

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

/* .......................................................... */

static int 
attach_subsample_to_dataset (DATAINFO *subinfo, double ***pZ, 
			     const DATAINFO *pdinfo)
{
    int i, t, n = pdinfo->n;

    /* no subsample currently in force */
    if (pZ == NULL) return 0;

    subinfo->subdum = malloc(n * sizeof *subinfo->subdum);
    if (subinfo->subdum == NULL) return E_ALLOC;

    i = varindex(pdinfo, "subdum");
    if (i == pdinfo->v) { /* safety measure: should be impossible */
	fprintf(stderr, "mystery failure in attach_subsample_to_dataset\n");
	return 1;   
    } 

    for (t=0; t<n; t++) {
	subinfo->subdum[t] = (*pZ)[i][t];
    }

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
	/* sync in case of any changes */
	fullinfo->varname = pdinfo->varname;
	fullinfo->varinfo = pdinfo->varinfo;

	pmod->subdum = copy_subdum(pdinfo->subdum, fullinfo->n);
	if (pmod->subdum == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* allocate_case_markers:
 * @n: number of observations.
 *
 * Allocate storage for a set of @n case markers or observation
 * labels.
 * 
 * Returns: pointer to storage, or %NULL if allocation fails.
 */

char **allocate_case_markers (int n)
{
    char **S;
    int t;

    S = malloc(n * sizeof *S);
    if (S == NULL) return NULL;

    for (t=0; t<n; t++) {
	S[t] = malloc(OBSLEN);
	if (S[t] == NULL) {
	    int j;

	    for (j=0; j<t; j++) free(S[j]);
	    free(S);
	    return NULL;
	}
	S[t][0] = '\0';
    }

    return S;
}

/* .......................................................... */

static void prep_subdinfo (DATAINFO *dinfo, int markers, int n)
{
    dinfo->sd0 = 1.;
    dinfo->pd = 1;
    dinfo->time_series = 0;
    strcpy(dinfo->stobs, "1");
    sprintf(dinfo->endobs, "%d", n);
    if (markers) {
	dinfo->markers = 1;
    } else {
	dinfo->markers = 0;
    }
}

/* .......................................................... */

static double *get_old_mask (double **Z, const DATAINFO *pdinfo)
{
    double *mask = NULL;
    int v = varindex(pdinfo, "subdum");

    if (v < pdinfo->v && isdummy(Z[v], 0, pdinfo->n - 1)) {
	mask = Z[v];
    }    

    return mask;
}

static int dummy_with_missing (const double *x, int t1, int t2)
{
    int t, m = 0;
    double xx;

    for (t=t1; t<=t2; t++) {
	xx = x[t];
	if (!na(xx) && floatneq(xx, 0.0) && floatneq(xx, 1.0)) { 
	    return 0;
	}
	if (floateq(xx, 1.0)) m++;
    }

    if (m < t2 - t1 + 1) return m;

    return 0;
} 

static int overlay_masks (double *targ, const double *src, int n)
{
    int i, sn = 0;
    
    for (i=0; i<n; i++) {
	if (targ[i] == 1.0 && src[i] == 1.0) {
	    targ[i] = 1.0;
	    sn++;
	} else {
	    targ[i] = 0.0;
	}
    }
	
    return sn;
}

static int make_boolean_mask (double ***pZ, DATAINFO *pdinfo, 
			      const char *line,
			      double *oldmask, int *dnum)
{
    char formula[MAXLEN];
    int t, subv = 0;

    if (oldmask != NULL) {
	/* copy across the old mask */
	subv = varindex(pdinfo, "subdum");
	for (t=0; t<pdinfo->n; t++) {
	    oldmask[t] = (*pZ)[subv][t];
	}
    }

    /* + 4 to skip the command word "smpl" */
    sprintf(formula, "__subdum=%s", line + 4);

    if (generate(pZ, pdinfo, formula, NULL))
	return -1;

    if (subv == 0) {
	subv = varindex(pdinfo, "subdum");
    }

    *dnum = subv;

    return isdummy((*pZ)[subv], pdinfo->t1, pdinfo->t2);
}

static int 
make_missing_mask (const double **Z, const DATAINFO *pdinfo,
		   const int *list, double *mask)
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
	mask[t] = 1.0;
	for (i=1; i<=lmax; i++) {
	    if (list != NULL) {
		v = list[i];
	    } else {
		v = i;
	    }
	    if (pdinfo->vector[v] && na(Z[v][t])) {
		mask[t] = 0.;
		break;
	    }
	}
	if (mask[t] == 1.0) sn++;
    }

    return sn;
} 

static int sn_from_dummy (const double **Z, const DATAINFO *pdinfo,
			  const char *dname, int *dnum)
{
    *dnum = varindex(pdinfo, dname);

    if (*dnum == pdinfo->v) {
	sprintf(gretl_errmsg, _("Variable '%s' not defined"), dname);
	return -1;
    } 

    return isdummy(Z[*dnum], pdinfo->t1, pdinfo->t2);
} 

static int count_selected_cases (double *x, int n)
{
    int i, count = 0;

    for (i=0; i<n; i++) {
	if (x[i] > 0.0) {
	    count++;
	}
    }

    return count;
}

static int make_random_mask (double *mask, double *oldmask, 
			     int fulln, int subn)
{
    int i, cases = 0, err = 0;
    unsigned u;

    if (subn <= 0 || subn >= fulln) {
	err = 1;
    } else if (oldmask != NULL) {
	int oldn = count_selected_cases(oldmask, fulln);

	if (subn >= oldn) err = 1;
    }	

    if (err) {
	sprintf(gretl_errmsg, _("Invalid number of cases %d"), subn);
	return 0;
    }	

    for (i=0; i<fulln; i++) {
	mask[i] = 0.0;
    }

    for (i=0; (cases != subn); i++) {
	u = gretl_rand_int_max(fulln);
	if (oldmask == NULL || oldmask[u] == 1.0) {
	    mask[u] = 1.0;
	}
	if (i >= subn - 1) {
	    cases = count_selected_cases(mask, fulln);
	}
    }

    return cases;
}

static int maybe_add_subdum (double ***pZ, DATAINFO *pdinfo, int *subnum)
{
    int v = varindex(pdinfo, "subdum");

    if (v == pdinfo->v) {
	/* variable doesn't exist: create it */
	if (dataset_add_vars(1, pZ, pdinfo)) {
	    return 1;
	}
	strcpy(pdinfo->varname[v], "subdum");
	strcpy(VARLABEL(pdinfo, v), _("automatic sub-sampling dummy"));
    } 

    *subnum = v;

    return 0;
}

enum {
    SUBSAMPLE_UNKNOWN,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN,
    SUBSAMPLE_RANDOM
} subsample_options;

static void backup_full_dataset (double ***pZ, DATAINFO **ppdinfo,
				 DATAINFO *newinfo)
{
    fullZ = *pZ;
    fullinfo = *ppdinfo;
    peerinfo = newinfo;
}

static void relink_full_dataset (double ***pZ, DATAINFO **ppdinfo)
{
    *pZ = fullZ;
    *ppdinfo = fullinfo;

    fullZ = NULL;
    fullinfo = NULL;
    peerinfo = NULL;
}

int complex_subsampled (void)
{
    if (fullZ == NULL) return 0;
    else return 1;
}

int get_full_length_n (void)
{
    if (fullinfo != NULL) {
	return fullinfo->n;
    } else {
	return 0;
    }
}

/* restrict_sample: 
 * @line: command line (or %NULL).  
 * @pZ: pointer to original data array.  
 * @ppdinfo: address of original data info pointer. 
 * @list: list of variables in case of OPT_M (or %NULL).  
 * @oflag: option flag.
 *
 * sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values (OPT_M); or using as a mask a
 * specified dummy variable (OPT_O); or masking with a specified
 * boolean condition (OPT_R); or selecting at random (OPT_N).
 *
 * In case OPT_M a @list of variables may be supplied; in cases
 * OPT_O, OPT_R and OPT_N, @line must contain specifics.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample (const char *line, 
		     double ***pZ, DATAINFO **ppdinfo, 
		     const int *list, gretlopt oflag)
{
    double xx, *tmpdum = NULL;
    char **S = NULL, dname[VNAMELEN] = {0};
    int subnum = 0;
    int i, t, st, sn = 0, n = (*ppdinfo)->n;
    int opt = SUBSAMPLE_UNKNOWN;
    double *mask = NULL, *oldmask = NULL;
    double **subZ = NULL;
    DATAINFO *subinfo = NULL;

    *gretl_errmsg = '\0';

    if (oflag & OPT_M) {
	opt = SUBSAMPLE_DROP_MISSING;
    } else if (oflag & OPT_R) {
	opt = SUBSAMPLE_BOOLEAN;
    } else if (oflag & OPT_N) {
	opt = SUBSAMPLE_RANDOM;
    } else if (oflag & OPT_O) {
	if (line != NULL && sscanf(line, "%*s %s", dname)) {
	    opt = SUBSAMPLE_USE_DUMMY;
	} else {
	    opt = SUBSAMPLE_DROP_MISSING;
	}
    } else {
	strcpy(gretl_errmsg, "Unrecognized sample command");
	return 1;
    }

    /* if no existing sub-sample mask, oldmask will be NULL */
    oldmask = get_old_mask(*pZ, *ppdinfo);

    if (opt == SUBSAMPLE_DROP_MISSING || opt == SUBSAMPLE_RANDOM ||
	(opt == SUBSAMPLE_BOOLEAN && oldmask != NULL)) {
	tmpdum = malloc(n * sizeof *tmpdum);
	if (tmpdum == NULL) return E_ALLOC;
	mask = tmpdum;
    }

    if (opt == SUBSAMPLE_DROP_MISSING) {   
	sn = make_missing_mask((const double **) *pZ, *ppdinfo, list, tmpdum);
    }  

    else if (opt == SUBSAMPLE_RANDOM) {
	sn = make_random_mask(tmpdum, oldmask, n, atoi(line + 4));
	if (sn == 0) {
	    free(tmpdum);
	    return 1;
	}
    }

    else if (opt == SUBSAMPLE_USE_DUMMY || opt == SUBSAMPLE_BOOLEAN) {
	int dnum;

	if (opt == SUBSAMPLE_USE_DUMMY) {
	    sn = sn_from_dummy((const double **) *pZ, *ppdinfo, dname, &dnum);
	    mask = (*pZ)[dnum];
	} else {
	    sn = make_boolean_mask(pZ, *ppdinfo, line, tmpdum, &dnum);
	    mask = (*pZ)[dnum];
	    oldmask = tmpdum;
	}
	if (sn < 0) {
	    return 1;
	}
    } 

    else {
	/* impossible */
	strcpy(gretl_errmsg, _("Sub-sample command failed mysteriously"));
	return 1;
    }

    /* cumulate restrictions, if appropriate */
    if (oldmask != NULL && opt != SUBSAMPLE_RANDOM) {
	sn = overlay_masks(mask, oldmask, n);
    }

    if (sn == 0) { /* "not a dummy variable"? */
	if (dummy_with_missing(mask, (*ppdinfo)->t1, (*ppdinfo)->t2)) {
	    strcpy(gretl_errmsg, _("Missing values found when applying criterion"));
	    return 1;
	}
    }

    /* does this policy lead to an empty sample, or no change
       in the sample, perchance? */
    if (sn == 0) {
	if (opt == SUBSAMPLE_USE_DUMMY) {
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), dname);
	} else if (opt == SUBSAMPLE_DROP_MISSING) {
	    strcpy(gretl_errmsg, _("No observations would be left!"));
	} else if (opt == SUBSAMPLE_BOOLEAN) { 
	    if (mask[(*ppdinfo)->t1] == 0) {
		strcpy(gretl_errmsg, _("No observations would be left!"));
	    } else {
		strcpy(gretl_errmsg, _("No observations were dropped!"));
	    }
	}
	return 1;
    }

    if (sn == n) {
	strcpy(gretl_errmsg, _("No observations were dropped!"));
	return 1;
    }

    /* allocate new datainfo */
    subinfo = malloc(sizeof *subinfo);
    if (subinfo == NULL) {
	return E_ALLOC;
    }

    /* create "hidden" dummy to record sub-sample, if need be */
    if (maybe_add_subdum(pZ, *ppdinfo, &subnum)) {
	free(tmpdum);
	free(subinfo);
	return E_ALLOC;
    } 

    /* write the new mask into the "subdum" variable */
    for (t=0; t<n; t++) {
	(*pZ)[subnum][t] = mask[t];
    }

    /* set up the sub-sampled datainfo */
    subinfo->n = sn;
    subinfo->v = (*ppdinfo)->v;
    if (start_new_Z(&subZ, subinfo, 1)) {
	free(tmpdum);
	free(subinfo);
	return E_ALLOC;
    }

    /* link varnames and descriptions (not dependent on series length) */
    subinfo->varname = (*ppdinfo)->varname;
    subinfo->varinfo = (*ppdinfo)->varinfo;
    subinfo->descrip = (*ppdinfo)->descrip;
    subinfo->vector = (*ppdinfo)->vector;

    /* case markers */
    if ((*ppdinfo)->markers) {
	S = allocate_case_markers(sn);
	if (S == NULL) {
	    free_Z(subZ, subinfo);
	    free(tmpdum);
	    free(subinfo);
	    return E_ALLOC;
	}
    }

    /* copy across data and case markers, if any */

    for (i=1; i<(*ppdinfo)->v; i++) {
	if (!(*ppdinfo)->vector[i]) {
	    /* copy any scalars */
	    subZ[i][0] = (*pZ)[i][0];
	}
    }

    st = 0;
    for (t=0; t<n; t++) {
	xx = mask[t];
	if (xx == 1.) {
	    for (i=1; i<(*ppdinfo)->v; i++) {
		if ((*ppdinfo)->vector[i]) {
		    subZ[i][st] = (*pZ)[i][t];
		}
	    } if ((*ppdinfo)->markers) {
		strcpy(S[st], (*ppdinfo)->S[t]);
	    }
	    st++;
	}
    }

    prep_subdinfo(subinfo, (*ppdinfo)->markers, sn);
    if ((*ppdinfo)->markers) {
	subinfo->S = S;
    }

    attach_subsample_to_dataset(subinfo, pZ, *ppdinfo);

    if (tmpdum != NULL) free(tmpdum);

    /* save state */
    backup_full_dataset(pZ, ppdinfo, subinfo);

    /* and switch pointers */
    *pZ = subZ;
    *ppdinfo = subinfo;

    return 0;
}

/* .......................................................... */

static int get_sample_increment (const char *str)
{
    if ((*str == '-' || *str == '+') && isdigit((unsigned char) str[1])) {
	return atoi(str);
    }

    return 0;
}

/* .......................................................... */

int set_sample (const char *line, DATAINFO *pdinfo)
{
    int nf, new_t1 = pdinfo->t1, new_t2 = pdinfo->t2;
    char cmd[5], newstart[OBSLEN], newstop[OBSLEN];

    *gretl_errmsg = '\0';

    nf = count_fields(line);

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
	    new_t1 = get_sample_increment(newstart);
	    if (new_t1) {
		new_t1 = pdinfo->t1 + new_t1;
		if (new_t1 < 0) {
		   strcpy(gretl_errmsg, _("Observation number out of bounds"));
		} 
	    } else {
		new_t1 = dateton(newstart, pdinfo);
	    }
	    if (*gretl_errmsg) return 1;
	    if (new_t1 < 0 || new_t1 > pdinfo->n) {
		strcpy(gretl_errmsg, _("error in new starting obs"));
		return 1;
	    }
	    pdinfo->t1 = new_t1;
	    return 0;
	}
    }

    if (sscanf(line, "%4s %10s %10s", cmd, newstart, newstop) != 3) {
	strcpy(gretl_errmsg, _("error reading smpl line"));
	return 1;
    }

    if (strcmp(newstart, ";")) {
	new_t1 = get_sample_increment(newstart);
	if (new_t1) {
	    new_t1 = pdinfo->t1 + new_t1;
	    if (new_t1 < 0) {
		strcpy(gretl_errmsg, _("Observation number out of bounds"));
	    }	
	} else {
	    new_t1 = dateton(newstart, pdinfo);
	}
	if (*gretl_errmsg) {
	    return 1;
	}
    }

    if (strcmp(newstop, ";")) {
	new_t2 = get_sample_increment(newstop);
	if (new_t2) {
	    new_t2 = pdinfo->t2 + new_t2;
	} else {
	    new_t2 = dateton(newstop, pdinfo);
	}
	if (*gretl_errmsg) return 1;
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

/* ........................................................... */

static int datamerge (double ***pZ, DATAINFO *pdinfo)
{
    int i, t, subt;
    int newvars = pdinfo->v - fullinfo->v;
    int n = fullinfo->n;
    double **newZ = NULL;
    int err = 0;

    if (newvars <= 0) return 0;

    if (pdinfo->subdum == NULL) {
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
    } else err = 1;

    if (err) return E_ALLOC;
    
    fullZ = newZ;

    for (i=fullinfo->v; i<pdinfo->v; i++) {
	if (!pdinfo->vector[i]) {
	   fullZ[i][0] = (*pZ)[i][0]; 
	}
    }

    subt = 0;
    for (t=0; t<n; t++) {
	if (pdinfo->subdum[t]) {
	    for (i=fullinfo->v; i<pdinfo->v; i++) {
		if (pdinfo->vector[i]) {
		    fullZ[i][t] = (*pZ)[i][subt];
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

/* ........................................................... */

static int make_smpl_mask (double ***pZ, DATAINFO *pdinfo)
{
    int old_vmax = pdinfo->v;
    int v, t;

    if (maybe_add_subdum(pZ, pdinfo, &v))
	return 1;

    if (v == old_vmax) {
	/* no pre-existing mask */
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[v][t] = (t < pdinfo->t1 || t > pdinfo->t2)? 0.0 : 1.0;
	}
    } else {
	/* there was a pre-existing mask */
	for (t=0; t<pdinfo->n; t++) {
	    if (t < pdinfo->t1 || t > pdinfo->t2) {
		(*pZ)[v][t] = 0.0;
	    }
	}	
    }

    return 0;
}

/* ........................................................... */

/* OPT_C -> replaCe sample restrictions */

int restore_full_sample (double ***pZ, DATAINFO **ppdinfo, gretlopt opt)
{
    int i, t, n, err = 0;

    *gretl_errmsg = '\0';

#ifdef SUBDEBUG
    fprintf(stderr, "restore_full_sample: pZ=%p, ppdinfo=%p opt=%lu\n",
	    (void *) pZ, (void *) ppdinfo, opt);
    fprintf(stderr, "*pZ=%p, *ppdinfo=%p\n",
	    (void *) *pZ, (void *) *ppdinfo);
#endif

    if (*pZ != NULL && !(opt & OPT_C)) {  
	/* FIXME */
	/* cumulating, not replacing restrictions */
	err = make_smpl_mask(pZ, *ppdinfo);
    }

    if (err) return err;

    if (!complex_subsampled()) {
	(*ppdinfo)->t1 = 0;
	(*ppdinfo)->t2 = (*ppdinfo)->n - 1;
	return 0;
    }

    /* set n to full series length */
    n = fullinfo->n;

    /* in case any new vars added, try to merge them in */
    err = datamerge(pZ, *ppdinfo);
    if (err == E_ALLOC) {
        sprintf(gretl_errmsg, _("Out of memory expanding data set\n"));
    } else if (err == E_NOMERGE) {
        sprintf(gretl_errmsg, 
		_("Missing sub-sample information; can't merge data\n"));
    }

    /* reattach malloc'd elements, which may have moved */
    fullinfo->varname = (*ppdinfo)->varname;
    fullinfo->varinfo = (*ppdinfo)->varinfo;
    fullinfo->vector = (*ppdinfo)->vector;
    fullinfo->descrip = (*ppdinfo)->descrip;  

    /* zero out the "subdum" dummy variable, if replacing
       sample restrictions */
    if (opt & OPT_C) {
	i = varindex(fullinfo, "subdum");
	if (i < fullinfo->v) {
	    for (t=0; t<n; t++) {
		fullZ[i][t] = 0.0;
	    }
	}
    }

    /* copy any scalars, which may have been modified */
    for (i=1; i<(*ppdinfo)->v; i++) {
	if (!(*ppdinfo)->vector[i]) {
	    fullZ[i][0] = (*pZ)[i][0];
	}
    }

    free_Z(*pZ, *ppdinfo);
    clear_datainfo(*ppdinfo, CLEAR_SUBSAMPLE);
    free(*ppdinfo);

    relink_full_dataset(pZ, ppdinfo);

    return 0;
}

/* ........................................................... */

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
	    if (hidden_var(i, pdinfo) || !pdinfo->vector[i]) continue;
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
		pprintf(prn, "%*s: %d %s\n", VNAMELEN, pdinfo->varname[i], 
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
	copy_varinfo(dest->varinfo[i], src->varinfo[i]);
	dest->vector[i] = src->vector[i];	
    }
}

/* Below: For a model that was estimated using a subsample which
   differs from the current sample, reconstitute the dataset
   appropriate to the model.  Use the special dummy variable
   saved with the model to select observations from the
   full dataset.  Attach this dataset to the model, so that
   it can be used for, e.g., residual or fitted versus actual
   plots.  The attached data set will be destroyed when the
   model itself is destroyed.
*/

int add_subsampled_dataset_to_model (MODEL *pmod)
{
    int i, t, st, sn = 0;
    double **modZ = NULL;
    char **S = NULL;

    if (fullZ == NULL || fullinfo == NULL) {
#ifdef SUBDEBUG
	fputs("add_subsampled_dataset_to_model: got NULL full dataset\n",
	      stderr);
#endif
	return 1;
    }

    if (pmod->dataset != NULL) {
#ifdef SUBDEBUG
	fputs("add_subsampled_dataset_to_model: job already done\n",
	      stderr);
#endif	
	return 0;
    }

    if (pmod->subdum == NULL) {
	/* no subsample info: model was estimated on full dataset */
	sn = fullinfo->n;
    } else {
	for (t=0; t<fullinfo->n; t++) {
	    if (pmod->subdum[t] > 0) sn++;
	}
	if (sn == 0) return 1;
    }

    pmod->dataset = malloc(sizeof *pmod->dataset);
    if (pmod->dataset == NULL) return E_ALLOC;

    pmod->dataset->dinfo = malloc(sizeof *pmod->dataset->dinfo);
    if (pmod->dataset->dinfo == NULL) {
	free(pmod->dataset);
	pmod->dataset = NULL;
	return E_ALLOC;	
    }

    /* set up the sub-sampled datainfo */
    pmod->dataset->dinfo->n = sn;
    pmod->dataset->dinfo->v = fullinfo->v;
    if (start_new_Z(&modZ, pmod->dataset->dinfo, 0)) {
	free(pmod->dataset->dinfo);
	free(pmod->dataset);
	pmod->dataset = NULL;
	return E_ALLOC;
    }

    /* copy across info on series */
    copy_series_info(pmod->dataset->dinfo, fullinfo);

    /* case markers */
    if (fullinfo->markers) {
	S = allocate_case_markers(sn);
	if (S == NULL) {
	    free_Z(modZ, pmod->dataset->dinfo);
	    free(pmod->dataset->dinfo);
	    free(pmod->dataset);
	    pmod->dataset = NULL;
	    return E_ALLOC;
	}
    }

    /* copy across data and case markers, if any */

    for (i=1; i<fullinfo->v; i++) {
	if (!fullinfo->vector[i]) {
	    /* copy any scalars */
	    modZ[i][0] = fullZ[i][0];
	}
    }

    st = 0;
    for (t=0; t<fullinfo->n; t++) {
	if (sn == fullinfo->n || pmod->subdum[t] == 1.) {
	    for (i=1; i<fullinfo->v; i++) {
		if (fullinfo->vector[i]) {
		    modZ[i][st] = fullZ[i][t];
		}
	    } if (fullinfo->markers) {
		strcpy(S[st], fullinfo->S[t]);
	    }
	    st++;
	}
    }

    pmod->dataset->Z = modZ;

    prep_subdinfo(pmod->dataset->dinfo, fullinfo->markers, sn);
    if (fullinfo->markers) {
	pmod->dataset->dinfo->S = S;
    }

#ifdef SUBDEBUG
    fputs("add_subsampled_dataset_to_model: success\n", stderr);
#endif

    return 0;
}

void free_model_dataset (MODEL *pmod)
{
    if (pmod->dataset != NULL) {
#ifdef SUBDEBUG
	fprintf(stderr, "Deep freeing model->dataset\n");
#endif
	free_Z(pmod->dataset->Z, pmod->dataset->dinfo);
	clear_datainfo(pmod->dataset->dinfo, CLEAR_FULL);
	free(pmod->dataset->dinfo);
	free(pmod->dataset);
	pmod->dataset = NULL;
    }
}
