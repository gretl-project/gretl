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

#define SUBDEBUG

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

static int 
attach_subsample_to_dataset (DATAINFO *pdinfo, double ***fullZ, 
			     const DATAINFO *fullinfo)
{
    int i, t, n = fullinfo->n;

    /* no subsample currently in force */
    if (fullZ == NULL) return 0;

    pdinfo->subdum = malloc(n * sizeof *pdinfo->subdum);
    if (pdinfo->subdum == NULL) return E_ALLOC;

    i = varindex(fullinfo, "subdum");
    if (i == fullinfo->v) { /* safety measure: should be impossible */
	fprintf(stderr, "mystery failure in attach_subsample_to_dataset\n");
	return 1;   
    } 

    for (t=0; t<n; t++) {
	pdinfo->subdum[t] = (*fullZ)[i][t];
    }

    return 0;
}

int attach_subsample_to_model (MODEL *pmod, const DATAINFO *pdinfo, int n)
     /* if the data set is currently subsampled, record the
	subsample info in pmod->subdum */
{
    /* no subsample currently in force */
    if (pdinfo == NULL || pdinfo->subdum == NULL) {
	return 0;
    }

    pmod->subdum = copy_subdum(pdinfo->subdum, n);
    if (pmod->subdum == NULL) return E_ALLOC;

    return 0;
}

/* .......................................................... */

int allocate_case_markers (char ***S, int n)
{
    int t;

    *S = malloc(n * sizeof **S);
    if (*S == NULL) {
	return E_ALLOC;
    }
    for (t=0; t<n; t++) {
	(*S)[t] = malloc(OBSLEN);
	if ((*S)[t] == NULL) {
	    free(*S);
	    return E_ALLOC;
	}
    }
    return 0;
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

static int make_boolean_mask (double ***pZ, DATAINFO *pdinfo, const char *line,
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

static int make_random_mask (double *mask, double *oldmask, int fulln, int subn)
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

enum {
    SUBSAMPLE_UNKNOWN,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN,
    SUBSAMPLE_RANDOM
} subsample_options;

/* restrict_sample: 
 * sub-sample the data set, based on the criterion of skipping all
 * observations with missing data values; or using as a mask a
 * specified dummy variable; or masking with a specified boolean
 * condition; or selecting at random.
*/

int restrict_sample (const char *line, 
		     double ***oldZ, double ***newZ,
		     DATAINFO *oldinfo, DATAINFO *newinfo,
		     const int *list, gretlopt oflag)
{
    double xx, *tmpdum = NULL;
    char **S = NULL, dname[VNAMELEN] = {0};
    int subnum = 0;
    int i, t, st, sn = 0, n = oldinfo->n;
    int opt = SUBSAMPLE_UNKNOWN;
    double *mask = NULL, *oldmask = NULL;

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
    oldmask = get_old_mask(*oldZ, oldinfo);

    if (opt == SUBSAMPLE_DROP_MISSING || opt == SUBSAMPLE_RANDOM ||
	(opt == SUBSAMPLE_BOOLEAN && oldmask != NULL)) {
	tmpdum = malloc(n * sizeof *tmpdum);
	if (tmpdum == NULL) return E_ALLOC;
	mask = tmpdum;
    }

    if (opt == SUBSAMPLE_DROP_MISSING) {   
	sn = make_missing_mask((const double **) *oldZ, oldinfo, list, tmpdum);
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
	    sn = sn_from_dummy((const double **) *oldZ, oldinfo, dname, &dnum);
	    mask = (*oldZ)[dnum];
	} else {
	    sn = make_boolean_mask(oldZ, oldinfo, line, tmpdum, &dnum);
	    mask = (*oldZ)[dnum];
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
	if (dummy_with_missing(mask, oldinfo->t1, oldinfo->t2)) {
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
	    if (mask[oldinfo->t1] == 0) {
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

    /* create "hidden" dummy to record sub-sample, if need be */
    subnum = varindex(oldinfo, "subdum");
    if (subnum == oldinfo->v) {
	/* variable doesn't exist: create it */
	if (dataset_add_vars(1, oldZ, oldinfo)) {
	    free(tmpdum);
	    return E_ALLOC;
	}
	strcpy(oldinfo->varname[subnum], "subdum");
	strcpy(VARLABEL(oldinfo, subnum), _("automatic sub-sampling dummy"));
    } 

    /* write the new mask into the "subdum" variable */
    for (t=0; t<n; t++) {
	(*oldZ)[subnum][t] = mask[t];
    }

    /* set up the sub-sampled datainfo */
    newinfo->n = sn;
    newinfo->v = oldinfo->v;
    if (start_new_Z(newZ, newinfo, 1)) {
	free(tmpdum);
	return E_ALLOC;
    }

    /* link varnames and descriptions (not dependent on series length) */
    newinfo->varname = oldinfo->varname;
    newinfo->varinfo = oldinfo->varinfo;
    newinfo->descrip = oldinfo->descrip;
    newinfo->vector = oldinfo->vector;

    /* case markers */
    if (oldinfo->markers && allocate_case_markers(&S, sn)) {
	free_Z(*newZ, newinfo);
	free(tmpdum);
	return E_ALLOC;
    }

    /* copy across data and case markers, if any */

    for (i=1; i<oldinfo->v; i++) {
	if (!oldinfo->vector[i]) {
	    /* copy any scalars */
	    (*newZ)[i][0] = (*oldZ)[i][0];
	}
    }

    st = 0;
    for (t=0; t<n; t++) {
	xx = mask[t];
	if (xx == 1.) {
	    for (i=1; i<oldinfo->v; i++) {
		if (oldinfo->vector[i]) {
		    (*newZ)[i][st] = (*oldZ)[i][t];
		}
	    } if (oldinfo->markers) {
		strcpy(S[st], oldinfo->S[t]);
	    }
	    st++;
	}
    }

    prep_subdinfo(newinfo, oldinfo->markers, sn);
    if (oldinfo->markers) newinfo->S = S;

    attach_subsample_to_dataset(newinfo, oldZ, oldinfo);

    if (tmpdum != NULL) free(tmpdum);

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
    int nf, new_t1 = 0, new_t2 = 0;
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

static int datamerge (double ***fullZ, DATAINFO *fullinfo,
		      double ***subZ, DATAINFO *subinfo)
{
    int i, t, subt;
    int newvars = subinfo->v - fullinfo->v;
    int n = fullinfo->n;
    double **newZ = NULL;
    int err = 0;

    if (newvars <= 0) return 0;

    if (subinfo->subdum == NULL) {
	return E_NOMERGE;
    }

    /* allocate expanded data array */
    newZ = realloc(*fullZ, subinfo->v * sizeof **fullZ);
    if (newZ != NULL) {
	for (i=0; i<newvars; i++) {
	    if (subinfo->vector[fullinfo->v+i]) {
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
    
    *fullZ = newZ;

    for (i=fullinfo->v; i<subinfo->v; i++) {
	if (!subinfo->vector[i]) {
	   (*fullZ)[i][0] = (*subZ)[i][0]; 
	}
    }

    subt = 0;
    for (t=0; t<n; t++) {
	if (subinfo->subdum[t]) {
	    for (i=fullinfo->v; i<subinfo->v; i++) {
		if (subinfo->vector[i]) {
		    (*fullZ)[i][t] = (*subZ)[i][subt];
		}
	    }
	    subt++;
	} else {
	    for (i=fullinfo->v; i<subinfo->v; i++) { 
		if (subinfo->vector[i]) {
		    (*fullZ)[i][t] = NADBL;
		}
	    }
	}
    }

    fullinfo->v = subinfo->v;

    return 0;
}

/* ........................................................... */

int restore_full_sample (double ***subZ, double ***fullZ, double ***Z,
			 DATAINFO **subinfo, DATAINFO **fullinfo,
			 DATAINFO **datainfo, gretlopt opt)
{
    int i, t, n, err = 0;

    *gretl_errmsg = '\0';

    /* simple case: merely a change of start or end of sample */
    if (*subZ == NULL) {
        (*datainfo)->t1 = 0;
        (*datainfo)->t2 = (*datainfo)->n - 1;
        return 0;
    }

    if (fullinfo == NULL || *fullinfo == NULL) return 1;

    /* reset n to full series length */
    n = (*fullinfo)->n;

    /* in case any new vars added, try to merge them in */
    err = datamerge(fullZ, *fullinfo, Z, *subinfo);
    if (err == E_ALLOC)
        sprintf(gretl_errmsg, _("Out of memory expanding data set\n"));
    if (err == E_NOMERGE)
        sprintf(gretl_errmsg, 
		_("Missing sub-sample information; can't merge data\n"));

    /* reattach the malloc'd elements, which might have moved */
    (*fullinfo)->varname = (*subinfo)->varname;
    (*fullinfo)->varinfo = (*subinfo)->varinfo;
    (*fullinfo)->vector = (*subinfo)->vector;
    (*fullinfo)->descrip = (*subinfo)->descrip;  

    /* zero out the "subdum" dummy variable, if not cumulating
       sample restrictions */
    if (!(opt & OPT_C)) {
	i = varindex(*fullinfo, "subdum");
	if (i < (*fullinfo)->v) {
	    for (t=0; t<n; t++) {
		(*fullZ)[i][t] = 0.0;
	    }
	}
    }

    /* copy any scalars, which may have been modified */
    for (i=1; i<(*subinfo)->v; i++) {
	if (!(*subinfo)->vector[i]) {
	    (*fullZ)[i][0] = (*Z)[i][0];
	}
    }

    /* reorganize pointers for data set */
    *subZ = *Z;
    *Z = *fullZ;
    free_Z(*subZ, *subinfo); 
    *subZ = NULL;
    *fullZ = NULL;

    /* and data info struct */
    *subinfo = *datainfo;
    *datainfo = *fullinfo;
    clear_datainfo(*subinfo, CLEAR_SUBSAMPLE);

    free(*subinfo);
    *subinfo = NULL;
    *fullinfo = NULL;
    
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

int add_subsampled_dataset_to_model (MODEL *pmod, 
				     const double **fullZ, 
				     const DATAINFO *fullinfo)
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
    if (fullinfo->markers && allocate_case_markers(&S, sn)) {
	free_Z(modZ, pmod->dataset->dinfo);
	free(pmod->dataset->dinfo);
	free(pmod->dataset);
	pmod->dataset = NULL;
	return E_ALLOC;
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
