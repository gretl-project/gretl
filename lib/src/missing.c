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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* mechanisms for handling missing observations */

#include "libgretl.h"
#include "gretl_private.h"

typedef struct {
    int misscount;
    char *missvec;
} MISSOBS;

int model_missval_count (const MODEL *pmod)
{
    int mc = 0;

    if (pmod->missmask != NULL) {
	int t;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (pmod->missmask[t - pmod->t1]) mc++;
	}
    }

    return mc;
}

/* The first set of functions here intended for use with daily data,
   where "missing observations" are likely to be non-existent
   observations (e.g. days on which trading did not take place due to
   holidays).
*/

/* reconstitute full-length series for residuals and fitted values,
   when the model was estimated using a data set from which missing
   observations had been temporarily purged.
*/

static int reorganize_uhat_yhat (MODEL *pmod) 
{
    MISSOBS *mobs = (MISSOBS *) pmod->data;
    double *tmp;
    int t, g;

    tmp = malloc(pmod->nobs * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    /* first do uhat */
    for (t=0; t<pmod->nobs; t++) {
	tmp[t] = pmod->uhat[t + pmod->t1];
    }

    g = 0;
    for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	if (mobs->missvec[t - pmod->t1]) {
	    pmod->uhat[t] = NADBL;
	} else {
	    pmod->uhat[t] = tmp[g++];
	}
    }

    /* then yhat */
    for (t=0; t<pmod->nobs; t++) {
	tmp[t] = pmod->yhat[t + pmod->t1];
    }

    g = 0;
    for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	if (mobs->missvec[t - pmod->t1]) {
	    pmod->yhat[t] = NADBL;
	} else {
	    pmod->yhat[t] = tmp[g++];
	}
    }

    free(tmp);

    return 0;
}

/* reverse the effect of "repack_missing" (below): put the
   observations that were reshuffled to the end of the sample range
   back into their proper places. 
*/

int undo_daily_repack (MODEL *pmod, double **Z, 
		       const DATAINFO *pdinfo)
{
    int i, j, t, m, g;
    double *tmpmiss = NULL, *tmpgood = NULL;
    /* note: pmod->t2 is still "shrunk" at this point */
    MISSOBS *mobs;
    int err = 0;

    if (!gretl_model_get_int(pmod, "daily_repack")) {
	return 1;
    }
    
    if (pmod->data == NULL) {
	return E_DATA;
    }
			   
    mobs = (MISSOBS *) pmod->data;

    tmpmiss = malloc(mobs->misscount * sizeof *tmpmiss);
    if (tmpmiss == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    tmpgood = malloc(pmod->nobs * sizeof *tmpgood);
    if (tmpgood == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (j=1; j<=pmod->list[0]; j++) {
	i = pmod->list[j];

	if (i == 0 || i == LISTSEP) {
	    continue;
	}

	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	     tmpgood[g++] = Z[i][t];
	}
	for (t=pmod->t2 + 1; t<=pmod->t2 + mobs->misscount; t++) {
	    tmpmiss[m++] = Z[i][t];
	}
	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	    if (mobs->missvec[t - pmod->t1]) {
		Z[i][t] = tmpmiss[m++];
	    } else {
		Z[i][t] = tmpgood[g++];
	    }
	}
    }

 bailout:

    free(tmpmiss);
    free(tmpgood);

    if (!err) {
	err = reorganize_uhat_yhat(pmod); 
    }  

    /* undo temporary shrinkage of pmod->t2 */
    pmod->t2 += mobs->misscount;

    /* trash the daily missing obs stuff */
    free(mobs->missvec);
    free(mobs);
    pmod->data = NULL;

    pmod->errcode = err;

    return err;
}

/* reshuffle missing observations from within a sample range to
   the end of the range (preserving relative order)
*/

static int 
repack_missing (const MODEL *pmod, double **Z, const DATAINFO *pdinfo, 
		const char *missvec, int misscount)
{
    int i, j, t, m, g;
    double *tmpmiss, *tmpgood;
    int modn = pmod->t2 - pmod->t1 + 1;

    tmpmiss = malloc(misscount * sizeof *tmpmiss);
    if (tmpmiss == NULL) {
	return 1;
    }

    tmpgood = malloc((modn - misscount) * sizeof *tmpgood);
    if (tmpgood == NULL) {
	free(tmpmiss);
	return 1;
    }

    for (j=1; j<=pmod->list[0]; j++) {
	i = pmod->list[j];

	if (i == 0 || i == LISTSEP) {
	    continue;
	}

	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (missvec[t - pmod->t1]) {
		tmpmiss[m++] = Z[i][t];
	    } else {
		tmpgood[g++] = Z[i][t];
	    }
	}
	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2 - misscount; t++) {
	     Z[i][t] = tmpgood[g++];
	}
	for (t=pmod->t2 + 1 - misscount; t<=pmod->t2; t++) {
	    Z[i][t] = tmpmiss[m++];
	}
    }

    free(tmpmiss);
    free(tmpgood);

    return 0;
}

/* construct a mask to code for missing observations
   within sample range for a model */

static char *model_missmask (const int *list, int t1, int t2,
			     const double **Z, int dwt,
			     int *misscount)
{
    char *mask;
    double xx;
    int i, j, t;

    mask = calloc(t2 - t1 + 1, sizeof *mask);

    if (mask == NULL) {
	return NULL;
    }

    if (misscount != NULL) {
	*misscount = 0;
    }

    for (t=t1; t<=t2; t++) {
	for (j=1; j<=list[0]; j++) {
	    i = list[j];
	    if (i == 0 || i == LISTSEP) {
		continue;
	    }
	    xx = Z[i][t];
	    if (dwt > 0) {
		xx *= Z[dwt][t];
	    }
	    if (na(xx)) {
		/* FIXME dwt case and nobs?? */
		mask[t - t1] = 1;
		if (misscount != NULL) {
		    *misscount += 1;
		}
		break;
	    }
	}
    }

    return mask;
}

/* with daily data, work around the missing obs by temporarily
   reorganizing the dataset 
*/

int repack_missing_daily_obs (MODEL *pmod, double **Z, 
			      const DATAINFO *pdinfo)
{
    char *missvec;
    int misscount;
    MISSOBS *mobs;
    int err = 0;

    missvec = model_missmask(pmod->list, pmod->t1, pmod->t2,
			     (const double **) Z, 0, &misscount);
    if (missvec == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    mobs = malloc(sizeof *mobs);
    if (mobs == NULL) {
	free(missvec);
	pmod->errcode = E_ALLOC;
	return 1;
    }
	
    err = repack_missing(pmod, Z, pdinfo, missvec, misscount);
    if (err) {
	pmod->errcode = E_ALLOC;
	free(missvec);
    } else {
	gretl_model_set_int(pmod, "daily_repack", 1);
	/* have to be sure to undo this after estimation */
	pmod->t2 -= misscount;
	mobs->missvec = missvec;
	mobs->misscount = misscount;
	pmod->data = mobs;
    }

    return err;
}

/* Drop first/last observations from sample if missing obs encountered.
   Also check for missing vals within the remaining sample: in case
   missing values are encountered there, either (a) construct a mask
   for them (if misst == NULL), or (b) flag an error.
*/

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 const double **Z, int *misst)
{
    int i, t, dwt = 0, t1min = *t1, t2max = *t2;
    int missobs, ret = 0;
    double xx;

    if (pmod != NULL && gretl_model_get_int(pmod, "wt_dummy")) {
	/* we have a weight variable which is a 0/1 dummy */
	dwt = pmod->nwt;
    }

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
	missobs = 0;
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) continue;
	    xx = Z[list[i]][t];
	    if (dwt) {
		xx *= Z[dwt][t];
	    }
	    if (na(xx)) {
		missobs = 1;
		break;
	    }
	}
	if (missobs) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
	missobs = 0;
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) continue;
	    xx = Z[list[i]][t];
	    if (dwt) {
		xx *= Z[dwt][t];
	    }
	    if (na(xx)) {
		missobs = 1;
		break;
	    }
	}
	if (missobs) {
	    t2max--;
	} else {
	    break;
	}	
    }

    /* check for missing values within remaining range */
    if (misst != NULL) {
	for (t=t1min; t<=t2max; t++) {
	    for (i=1; i<=list[0]; i++) {
		if (list[i] == LISTSEP) continue;
		xx = Z[list[i]][t];
		if (dwt) {
		    xx *= Z[dwt][t];
		}
		if (na(xx)) {
		    *misst = t + 1;
		    ret = list[i];
		    break;
		}
	    }
	    if (ret) {
		break;
	    }
	}     
    }

    /* construct a mask for missing values within remaining range? 
       Note: we do this only if misst == NULL 
    */
    else if (pmod != NULL) {
	missobs = 0;
	for (t=t1min; t<=t2max; t++) {
	    for (i=1; i<=list[0]; i++) {
		if (list[i] == LISTSEP) continue;
		xx = Z[list[i]][t];
		if (dwt) {
		    xx *= Z[dwt][t];
		}
		if (na(xx)) {
		    missobs++;
		    break;
		}
	    }
	}
	if (missobs > 0) {
	    /* FIXME: special treatment if no valid obs left? */
	    pmod->missmask = model_missmask(list, t1min, t2max, 
					    Z, dwt, NULL);
	}
    }    

    *t1 = t1min; 
    *t2 = t2max;

    return ret;
}


/* ........................................................... */

static int real_setmiss (double missval, int varno, 
			 double **Z, DATAINFO *pdinfo) 
{
    int i, t, count = 0;
    int start = 1, end = pdinfo->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	for (t=0; t<pdinfo->n; t++) {
	    if (Z[i][t] == missval) {
		Z[i][t] = NADBL;
		count++;
	    }
	}	
    }

    return count;
}

/**
 * set_miss:
 * @list: list of variables to process.
 * @param: string with specification of value to treat as missing.
 * @Z: data matrix.
 * @pdinfo: pointer to data information struct.
 * @PRN: pointer to printing struct.
 * 
 * Set to "missing" each observation of each series in list that
 * has the specified value, as in @param.
 *
 */

void set_miss (const int *list, const char *param, double **Z,
	       DATAINFO *pdinfo, PRN *prn)
{
    double missval;
    int i, count;

    missval = atof(param);

    if (list[0] == 0) {
	count = real_setmiss(missval, 0, Z, pdinfo);
	if (count) { 
	    pprintf(prn, _("Set %d values to \"missing\"\n"), count);
	} else {
	    pputs(prn, _("Didn't find any matching observations\n"));
	}
	return;
    }

    for (i=1; i<=list[0]; i++) {
	if (!pdinfo->vector[list[i]]) {
	    pprintf(prn, _("The variable %s is a scalar\n"), 
		    pdinfo->varname[list[i]]);
	    continue;
	}
	count = real_setmiss(missval, list[i], Z, pdinfo);
	if (count) { 
	    pprintf(prn, _("%s: set %d observations to \"missing\"\n"), 
		    pdinfo->varname[list[i]], count);
	} else { 
	    pprintf(prn, _("%s: Didn't find any matching observations\n"),
		    pdinfo->varname[list[i]]);
	}
    }
}
