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

/* .......................................................... */

int attach_subsample_to_model (MODEL *pmod, double ***fullZ, 
			       const DATAINFO *fullinfo)
     /* if the data set is currently subsampled, record the
	subsample info in pmod->subdum */
{
    int i, t, n = fullinfo->n;

    /* no subsample currently in force */
    if (fullZ == NULL) return 0;

    pmod->subdum = malloc(n * sizeof(double));
    if (pmod->subdum == NULL) return E_ALLOC;

    i = varindex(fullinfo, "subdum");
    if (i == fullinfo->v) { /* safety measure: should be impossible */
	fprintf(stderr, _("mystery failure in attach_subsample_to_model\n"));
	return 1;   
    } 

    for (t=0; t<n; t++)
	pmod->subdum[t] = (*fullZ)[i][t];
    
    return 0;
}

/* .......................................................... */

static int subsampled (double **Z, const DATAINFO *pdinfo, 
		       const int subnum)
     /* Is the data set currently "sub-sampled" via selection of 
	cases?  Use this func _only_ if subnum tests < pdinfo->v */
{
    int t, n = pdinfo->n;

    for (t=0; t<n; t++)
	if (floatneq(Z[subnum][t], 0.0)) return 1;
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
    extern int _identical (const double *x, const double *y, const int n);

    if (pmod == NULL && spec == NULL) return 0;

    /* if no sub-sampling has been done, we're OK */
    if ((i = varindex(pdinfo, "subdum")) == pdinfo->v) return 0;

    if (pmod != NULL) subdum = pmod->subdum;
    else subdum = spec->subdum;

    /* case: model has no sub-sampling info recorded */
    if (subdum == NULL) {
	/* if data set is not currently sub-sampled, we're OK */
	if (!subsampled(Z, pdinfo, i)) return 0;
	/* data set is subsampled, model is not: problem */
	else {
	    fprintf(stderr, _("dataset is subsampled, model is not\n"));
	    return 1;
	}
    }

    /* case: model has sub-sampling info recorded */
    if (!subsampled(Z, pdinfo, i)) {
	/* data set not sub-sampled: problem */
	fprintf(stderr, _("model is subsampled, dataset is not\n"));
	return 1;
    } else { /* do the subsamples (model and current data set) agree? */
	if (_identical(Z[i], subdum, n))
	    return 0;
	else {
	    fprintf(stderr, _("model and dataset subsamples not the same\n"));
	    return 1;
	}
    }

    /* can't be reached */
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
	(*S)[t] = malloc(9);
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
    dinfo->extra = 0;
    if (markers) dinfo->markers = 1;
    else dinfo->markers = 0;
    strcpy(dinfo->stobs, "1");
    sprintf(dinfo->endobs, "%d", n);
}

/* .......................................................... */

int set_sample_dummy (const char *line, 
		      double ***oldZ, double ***newZ,
		      DATAINFO *oldinfo, DATAINFO *newinfo,
		      const int opt)
     /* sub-sample the data set, based on the criterion of skipping
	all observations with missing data values; or using as a
	mask a specified dummy variable;, or masking with a specified
	boolean condition */
{
    double xx, *dum = NULL;
    char **S = NULL, dumv[9];
    int missobs = 0, subnum = 0, dumnum = 0;
    int i, t, st, sn, n = oldinfo->n;

    gretl_errmsg[0] = '\0';

    dumv[0] = '\0';
    if (opt == OPT_O && 
	(line == NULL || sscanf(line, "%*s %s", dumv) <= 0))
	missobs = 1; 

    if (missobs) { /* construct missing obs dummy on the fly */
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
    else if (opt == OPT_O) { /* the name of a dummy var was given */ 
	dumnum = varindex(oldinfo, dumv);
	if (dumnum == oldinfo->v) {
	    sprintf(gretl_errmsg, _("Variable '%s' not defined"), dumv);
	    return 1;
	} 
	sn = isdummy(dumnum, oldinfo->t1, oldinfo->t2, *oldZ);
    } 
    else if (opt == OPT_R) { /* construct dummy from boolean */
	char formula[MAXLEN];
	int err;

	/* + 4 below to omit the word "smpl" */
	sprintf(formula, "subdum=%s", line + 4);
	err = generate(oldZ, oldinfo, formula, 0, NULL, 1);
	if (err) return err;
	subnum = varindex(oldinfo, "subdum");
	dumnum = subnum;
	sn = isdummy(subnum, oldinfo->t1, oldinfo->t2, *oldZ);
    } else {
	/* impossible */
	strcpy(gretl_errmsg, _("Sub-sample command failed mysteriously"));
	return 1;
    }

    /* does this policy lead to an empty sample, or no change
       in the sample, perchance? */
    if (sn == 0) {
	if (opt == OPT_O && !missobs)
	    sprintf(gretl_errmsg, _("'%s' is not a dummy variable"), dumv);
	else if (missobs)
	    strcpy(gretl_errmsg, _("No observations would be left!"));
	else { /* case of boolean expression */
	    if ((*oldZ)[subnum][oldinfo->t1] == 0)
		strcpy(gretl_errmsg, _("No observations would be left!"));
	    else
		strcpy(gretl_errmsg, _("No observations were dropped!"));
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
	strcpy(oldinfo->label[subnum], _("automatic sub-sampling dummy"));
    }

    for (t=0; t<n; t++) {
	if (missobs) 
	    (*oldZ)[subnum][t] = dum[t];
	else if (opt == OPT_O)
	    /* ?possibility of missing values here? */
	    (*oldZ)[subnum][t] = (*oldZ)[dumnum][t];
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
    newinfo->label = oldinfo->label;
    newinfo->descrip = oldinfo->descrip;
    newinfo->vector = oldinfo->vector;

    /* case markers */
    if (oldinfo->markers && allocate_case_markers(&S, sn)) {
	free_Z(*newZ, newinfo);
	free(dum);
	return E_ALLOC;
    }

    /* copy across data and case markers, if any */
    for (i=1; i<oldinfo->v; i++) 
	if (!oldinfo->vector[i])
	    (*newZ)[i][0] = (*oldZ)[i][0];
    st = 0;
    for (t=0; t<n; t++) {
	xx = (missobs)? dum[t] : (*oldZ)[dumnum][t];
	if (xx == 1.) {
	    for (i=1; i<oldinfo->v; i++) {
		if (oldinfo->vector[i]) 
		    (*newZ)[i][st] = (*oldZ)[i][t];
	    } if (oldinfo->markers) 
		strcpy(S[st], oldinfo->S[t]);
	    st++;
	}
    }

    prep_subdinfo(newinfo, oldinfo->markers, sn);
    if (oldinfo->markers) newinfo->S = S;

    if (dum != NULL) free(dum);

    return 0;
}

/* .......................................................... */

int set_sample (const char *line, DATAINFO *pdinfo)
{
    int nf, new_t1 = 0, new_t2 = 0;
    char cmd[5], newstart[9], newstop[9];

    gretl_errmsg[0] = '\0';

    nf = _count_fields(line);

    if (nf == 1) return 0;
	
    if (nf == 2) {
	if (sscanf(line, "%s %s", cmd, newstart) != 2) {
	    sprintf(gretl_errmsg, _("error reading smpl line"));
	    return 1;
	} else {
	    new_t1 = dateton(newstart, pdinfo);
	    if (new_t1 < 0 || strlen(gretl_errmsg)) return 1;
	    if (new_t1 > pdinfo->n) {
		sprintf(gretl_errmsg, _("error in new starting obs"));
		return 1;
	    }
	    pdinfo->t1 = new_t1;
	    return 0;
	}
    }
    if (sscanf(line, "%s %s %s", cmd, newstart, newstop) != 3) {
	sprintf(gretl_errmsg, _("error reading smpl line"));
	return 1;
    }
    if (strcmp(newstart, ";")) {
	new_t1 = dateton(newstart, pdinfo);
	if (new_t1 < 0 || strlen(gretl_errmsg)) {
	    return 1;
	}
    }
    if (strcmp(newstop, ";")) {
	new_t2 = dateton(newstop, pdinfo);
	if (strlen(gretl_errmsg)) return 1;
	if (new_t2 >= pdinfo->n) {
	    sprintf(gretl_errmsg, _("error in new ending obs"));
	    return 1;
	}
    }

    if (new_t1 > new_t2) {
	sprintf(gretl_errmsg, _("Invalid null sample"));
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

    gretl_errmsg[0] = '\0';

    /* simple case: merely a change of start or end of sample */
    if (*subZ == NULL) {
        (*datainfo)->t1 = 0;
        (*datainfo)->t2 = (*datainfo)->n - 1;
        return 0;
    }

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
    (*fullinfo)->label = (*subinfo)->label;
    (*fullinfo)->vector = (*subinfo)->vector;
    (*fullinfo)->descrip = (*subinfo)->descrip;    

    /* zero out the "subdum" dummy variable */
    i = varindex(*fullinfo, "subdum");
    if (i < (*fullinfo)->v)
        for (t=0; t<n; t++) (*fullZ)[i][t] = 0.;

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
    int i, v, t;
    int missval = 0, missobs = 0, oldmiss = 0, tmiss;
    int year = 0, yearmiss = 0, totvals = 0, yearbak = 0;

    v = varindex(pdinfo, "year");
    if (v == pdinfo->v) v = varindex(pdinfo, "YEAR");
    if (v == pdinfo->v) v = 0;
    else yearbak = (int) (*pZ)[v][pdinfo->t1];

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	tmiss = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (hidden_var(i, pdinfo)) continue;
	    if (na((*pZ)[i][t])) missval++;
	    totvals++;
	}
	if ((tmiss = missval - oldmiss)) missobs++;
	if (v) {
	    year = (int) (*pZ)[v][t];
	    if (year != yearbak) {
		pprintf(prn, _("%d: %4d missing data values\n"), 
			yearbak, yearmiss);
		yearmiss = tmiss;
		yearbak = year;
	    } else
		yearmiss += tmiss;
	}
	oldmiss = missval;
    }
    if (v) 
	pprintf(prn, _("%d: %4d missing data values\n"), 
		year, yearmiss);
    
    pprintf(prn, _("\nNumber of observations (rows) with missing data "
	    "values = %d (%d%%)\n"), missobs, 
	    (int) (100.0 * missobs / (pdinfo->t2 - pdinfo->t1 + 1)));
    pprintf(prn, _("Total number of missing data values = %d (%d%% "
	    "of total data values)\n"), missval, 
	    (int) (100.0 * missval / totvals));
    return missval;
}
