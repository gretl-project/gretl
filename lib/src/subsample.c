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
#include "internal.h"

#define SUBDEBUG

/* .......................................................... */

int attach_subsample_to_model (MODEL *pmod, double ***fullZ, 
			       const DATAINFO *fullinfo)
     /* if the data set is currently subsampled, record the
	subsample info in pmod->subdum */
{
    int i, t, n = fullinfo->n;

    /* no subsample currently in force */
    if (fullZ == NULL) return 0;

    pmod->subdum = malloc(n * sizeof *pmod->subdum);
    if (pmod->subdum == NULL) return E_ALLOC;

    i = varindex(fullinfo, "subdum");
    if (i == fullinfo->v) { /* safety measure: should be impossible */
	fprintf(stderr, I_("mystery failure in attach_subsample_to_model\n"));
	return 1;   
    } 

    for (t=0; t<n; t++) {
	pmod->subdum[t] = (*fullZ)[i][t];
    }
    
    return 0;
}

/* .......................................................... */

static int subsampled (double **Z, const DATAINFO *pdinfo, 
		       const int subnum)
     /* Is the data set currently "sub-sampled" via selection of 
	cases?  Use this func _only_ if subnum tests < pdinfo->v,
        i.e. only if the dummy variable exists. 
     */
{
    int t, n = pdinfo->n;

    for (t=0; t<n; t++) {
	if (floatneq(Z[subnum][t], 0.0)) return 1;
    }
    return 0;
}

/* .......................................................... */

int model_sample_issue (const MODEL *pmod, MODELSPEC *spec, 
			double **Z, const DATAINFO *pdinfo)
     /* check a model (or modelspec) against the data info to see if 
	it may have been estimated on a different (subsampled) data 
	set from the current one */
{
    int i, n = pdinfo->n;
    double *subdum;

    if (pmod == NULL && spec == NULL) return 0;

    /* if no sub-sampling has been done, we're OK */
    if ((i = varindex(pdinfo, "subdum")) == pdinfo->v) return 0;

    if (pmod != NULL) subdum = pmod->subdum;
    else subdum = spec->subdum;

    /* case: model has no sub-sampling info recorded */
    if (subdum == NULL) {
	/* if data set is not currently sub-sampled, we're OK */
	if (!subsampled(Z, pdinfo, i)) {
	    return 0;
	} else {
	    strcpy(gretl_errmsg, _("dataset is subsampled, model is not\n"));
	    return 1;
	}
    }

    /* case: model has sub-sampling info recorded */
    if (!subsampled(Z, pdinfo, i)) {
	strcpy(gretl_errmsg, _("model is subsampled, dataset is not\n"));
	return 1;
    } else { 
	/* do the subsamples (model and current data set) agree? */
	if (vars_identical(Z[i], subdum, n)) {
	    return 0;
	} else {
	    strcpy(gretl_errmsg, _("model and dataset subsamples not the same\n"));
	    return 1;
	}
    }

    /* not reached */
    return 1;
}

/* .......................................................... */

int allocate_case_markers (char ***S, int n)
{
    int t;

    *S = malloc(n * sizeof(char *));
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
    if (markers) dinfo->markers = 1;
    else dinfo->markers = 0;
    strcpy(dinfo->stobs, "1");
    sprintf(dinfo->endobs, "%d", n);
}

/* .......................................................... */

static int dummy_with_missing (double *x, int t1, int t2)
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

/* .......................................................... */

enum {
    SUBSAMPLE_UNKNOWN,
    SUBSAMPLE_DROP_MISSING,
    SUBSAMPLE_USE_DUMMY,
    SUBSAMPLE_BOOLEAN
} subsample_options;

int restrict_sample (const char *line, 
		     double ***oldZ, double ***newZ,
		     DATAINFO *oldinfo, DATAINFO *newinfo,
		     unsigned long oflag)
     /* sub-sample the data set, based on the criterion of skipping
	all observations with missing data values; or using as a
	mask a specified dummy variable; or masking with a specified
	boolean condition */
{
    double xx, *dum = NULL;
    char **S = NULL, dumv[VNAMELEN];
    int subnum = 0, dumnum = 0;
    int i, t, st, sn, n = oldinfo->n;
    int opt = SUBSAMPLE_UNKNOWN;

    *gretl_errmsg = '\0';
    *dumv = '\0';

    if (oflag == OPT_O && 
	(line == NULL || sscanf(line, "%*s %s", dumv) <= 0)) {
	opt = SUBSAMPLE_DROP_MISSING;
    } else if (oflag == OPT_O) {
	opt = SUBSAMPLE_USE_DUMMY;
    } else if (oflag == OPT_R) {
	opt = SUBSAMPLE_BOOLEAN;
    } else {
	strcpy(gretl_errmsg, "Unrecognized sample command");
	return 1;
    }

    if (opt == SUBSAMPLE_DROP_MISSING) { 
	dum = malloc(n * sizeof *dum);
	if (dum == NULL) return E_ALLOC;
	sn = 0;
	for (t=0; t<n; t++) {
	    dum[t] = 1.0;
	    for (i=1; i<oldinfo->v; i++) {
		if (oldinfo->vector[i] && na((*oldZ)[i][t])) {
		    dum[t] = 0.;
		    break;
		}
	    }
	    if (floateq(dum[t], 1.0)) sn++;
	}
    } 

    else if (opt == SUBSAMPLE_USE_DUMMY) { 
	dumnum = varindex(oldinfo, dumv);
	if (dumnum == oldinfo->v) {
	    sprintf(gretl_errmsg, _("Variable '%s' not defined"), dumv);
	    return 1;
	} 
	sn = isdummy((*oldZ)[dumnum], oldinfo->t1, oldinfo->t2);
    } 

    else if (opt == SUBSAMPLE_BOOLEAN) { 
	char formula[MAXLEN];
	int err;

	/* + 4 below to omit the word "smpl" */
	sprintf(formula, "subdum=%s", line + 4);
	err = generate(oldZ, oldinfo, formula, 0, NULL, 1);
	if (err) return err;
	subnum = varindex(oldinfo, "subdum");
	dumnum = subnum;
	sn = isdummy((*oldZ)[subnum], oldinfo->t1, oldinfo->t2);
    } 
    else {
	/* impossible */
	strcpy(gretl_errmsg, _("Sub-sample command failed mysteriously"));
	return 1;
    }

    if (sn == 0) { /* "not a dummy variable" */
	if (dummy_with_missing((*oldZ)[subnum], oldinfo->t1, oldinfo->t2)) {
	    strcpy(gretl_errmsg, _("Missing values found when applying criterion"));
	    return 1;
	}
    }

    /* does this policy lead to an empty sample, or no change
       in the sample, perchance? */
    if (sn == 0) {
	if (opt == SUBSAMPLE_USE_DUMMY) {
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), dumv);
	} else if (opt == SUBSAMPLE_DROP_MISSING) {
	    strcpy(gretl_errmsg, _("No observations would be left!"));
	} else { 
	    /* case of boolean expression */
	    if ((*oldZ)[subnum][oldinfo->t1] == 0) {
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

    /* create or reuse "hidden" dummy to record sub-sample */
    subnum = varindex(oldinfo, "subdum");
    if (subnum == oldinfo->v) {
	if (dataset_add_vars(1, oldZ, oldinfo)) return E_ALLOC;
	strcpy(oldinfo->varname[subnum], "subdum");
	strcpy(VARLABEL(oldinfo, subnum), _("automatic sub-sampling dummy"));
    }

    for (t=0; t<n; t++) {
	if (opt == SUBSAMPLE_DROP_MISSING) {
	    (*oldZ)[subnum][t] = dum[t];
	} else if (opt == SUBSAMPLE_USE_DUMMY) {
	    /* possibility of missing values here? */
	    (*oldZ)[subnum][t] = (*oldZ)[dumnum][t];
	}
    }

    /* set up the sub-sampled datainfo */
    newinfo->n = sn;
    newinfo->v = oldinfo->v;
    if (start_new_Z(newZ, newinfo, 1)) {
	if (dum != NULL) free(dum);
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
	free(dum);
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
	xx = (opt == SUBSAMPLE_DROP_MISSING)? dum[t] : (*oldZ)[dumnum][t];
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

    if (dum != NULL) free(dum);

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

    nf = _count_fields(line);

    if (nf == 3 && pdinfo->n == 0) {
	/* database special */
	return db_set_sample(line, pdinfo);
    }

    if (nf == 1) return 0;
	
    if (nf == 2) {
	if (sscanf(line, "%4s %8s", cmd, newstart) != 2) {
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

    if (sscanf(line, "%4s %8s %8s", cmd, newstart, newstop) != 3) {
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
    int i, t, dumn, subt;
    int newvars = subinfo->v - fullinfo->v;
    int n = fullinfo->n;
    double **newZ = NULL;
    int err = 0;

    if (newvars <= 0) return 0;

    dumn = varindex(subinfo, "subdum"); 
    if (dumn == subinfo->v) return E_NOMERGE;

    /* allocate expanded data array */
    newZ = realloc(*fullZ, subinfo->v * sizeof **fullZ);
    if (newZ != NULL) {
	for (i=0; i<newvars; i++) {
	    if (subinfo->vector[fullinfo->v+i])
		newZ[fullinfo->v+i] = malloc(n * sizeof **newZ);
	    else
		newZ[fullinfo->v+i] = malloc(sizeof **newZ);
	    if (newZ[fullinfo->v+i] == NULL) {
		err = 1;
		break;
	    }
	}
    } else err = 1;

    if (err) return E_ALLOC;
    else *fullZ = newZ;

    for (i=fullinfo->v; i<subinfo->v; i++) 
	if (!subinfo->vector[i])
	   (*fullZ)[i][0] = (*subZ)[i][0]; 

    subt = 0;
    for (t=0; t<n; t++) {
	if ((*fullZ)[dumn][t] == 1.0) {
	    for (i=fullinfo->v; i<subinfo->v; i++) {
		if (subinfo->vector[i]) 
		    (*fullZ)[i][t] = (*subZ)[i][subt];
	    }
	    subt++;
	} else {
	    for (i=fullinfo->v; i<subinfo->v; i++) { 
		if (subinfo->vector[i]) 
		    (*fullZ)[i][t] = NADBL;
	    }
	}
    }

    fullinfo->v = subinfo->v;

    return 0;
}

/* ........................................................... */

int restore_full_sample (double ***subZ, double ***fullZ, double ***Z,
			 DATAINFO **subinfo, DATAINFO **fullinfo,
			 DATAINFO **datainfo)
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

    /* zero out the "subdum" dummy variable */
    i = varindex(*fullinfo, "subdum");
    if (i < (*fullinfo)->v) {
        for (t=0; t<n; t++) {
	    (*fullZ)[i][t] = 0.0;
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
