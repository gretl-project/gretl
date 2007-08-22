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

/* mechanisms for handling missing observations */

#include "libgretl.h"
#include "missing_private.h"

typedef struct MISSOBS_ MISSOBS;

struct MISSOBS_ {
    int misscount;
    char *missvec;
};

static char *model_missmask (const int *list, int t1, int t2,
			     int n, const double **Z, int dwt,
			     int *misscount);

/* The first set of functions here intended for use with daily data,
   where "missing observations" are likely to be non-existent
   observations (e.g. days on which trading did not take place due to
   holidays).
*/

/* reconstitute full-length series for residuals and fitted values,
   when the model was estimated using a data set from which missing
   observations had been temporarily purged.
*/

static int reorganize_uhat_yhat (MODEL *pmod, MISSOBS *mobs) 
{
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
	if (mobs->missvec[t] == '1') {
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
	if (mobs->missvec[t] == '1') {
	    pmod->yhat[t] = NADBL;
	} else {
	    pmod->yhat[t] = tmp[g++];
	}
    }

    free(tmp);

    return 0;
}

/**
 * undo_daily_repack:
 * @pmod: pointer to model.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Reverses the effect of repack_missing_daily_obs(), hence
 * restoring a daily dataset to its correct order.  The
 * model pointer, @pmod, must be the same as that given 
 * previously to repack_missing_daily_obs().
 *
 * Returns: 0 on success, non-zero error code on error.
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

    mobs = (MISSOBS *) gretl_model_get_data(pmod, "missobs");
    if (mobs == NULL) {
	return E_DATA;
    }

    /* take charge of MISSOBS pointer */
    gretl_model_detach_data_item(pmod, "missobs");
			   
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
	    if (mobs->missvec[t] == '1') {
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
	err = reorganize_uhat_yhat(pmod, mobs); 
    }  

    /* undo temporary shrinkage of pmod->t2 */
    pmod->t2 += mobs->misscount;

    /* Trash the daily missing obs stuff. TODO: might save this info
       for future use (e.g. in forecasting) */
    free(mobs->missvec);
    free(mobs);

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
	    if (missvec[t] == '1') {
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

/**
 * repack_missing_daily_obs:
 * @pmod: pointer to model.
 * @Z: data array.
 * @pdinfo: information on dataset.
 *
 * Reorganizes a daily dataset, moving any missing observations
 * to the end of the array.  Permits estimation of a model
 * using daily data, where one wants to treat missing values
 * as not really "missing" but rather non-existent (e.g. values
 * of financial market variables on dates of trading holidays).
 * The original dataset can (and generally should) be restored
 * using undo_daily_repack() after model estimation.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int repack_missing_daily_obs (MODEL *pmod, double **Z, 
			      const DATAINFO *pdinfo)
{
    char *missvec;
    int misscount;
    MISSOBS *mobs;
    int err = 0;

    missvec = model_missmask(pmod->list, pmod->t1, pmod->t2, 
			     pdinfo->n, (const double **) Z, 
			     0, &misscount);
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
	/* Be sure to undo this after estimation! */
	pmod->t2 -= misscount;
	mobs->missvec = missvec;
	mobs->misscount = misscount;
	err = gretl_model_set_data(pmod, "missobs", mobs, 
				   MODEL_DATA_STRUCT,
				   sizeof *mobs);
    }

    return err;
}

#define MASKDEBUG 0

/**
 * model_missval_count:
 * @pmod: pointer to model.
 *
 * Returns: a count of the missing values within the sample
 * range over which @pmod was estimated.
 */

int model_missval_count (const MODEL *pmod)
{
    int mc = 0;

    if (pmod->missmask != NULL) {
	int t;

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (pmod->missmask[t] == '1') {
		mc++;
	    }
	}
    }

    return mc;
}

/* Construct a mask to code for missing observations within the sample
   range for a model.  The mask is an array of char with '0's 
   for OK observations and '1's for missing obs, terminated with
   a NUL byte.
*/

static char *model_missmask (const int *list, int t1, int t2,
			     int n, const double **Z, int dwt,
			     int *misscount)
{
    char *mask;
    double xx;
    int i, li, t;

    mask = malloc(n + 1);
    if (mask == NULL) {
	return NULL;
    }

    memset(mask, '0', n);
    mask[n] = 0;

#if MASKDEBUG
    fprintf(stderr, "model_missmask: using series length %d\n", n);
#endif

    if (misscount != NULL) {
	*misscount = 0;
    }

    for (t=t1; t<=t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    li = list[i];
	    if (li == 0 || li == LISTSEP) {
		continue;
	    }
	    xx = Z[li][t];
	    if (dwt > 0) {
		/* dummy weight variable */
		xx *= Z[dwt][t];
	    }
	    if (na(xx)) {
#if MASKDEBUG
		fprintf(stderr, "model_missmask: NA at list[%d] (%d), obs %d\n",
			i, list[i], t);
#endif
		/* FIXME dwt case and nobs?? */
		mask[t] = '1';
		if (misscount != NULL) {
		    *misscount += 1;
		}
		break;
	    }
	}
    }

    return mask;
}

/**
 * adjust_t1t2: 
 * @pmod: pointer to model, or %NULL.
 * @list: list of variables to be tested for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @n: full length of data array.
 * @Z: data array.
 * @misst: location to receive the first observation with a
 *         missing value inside the sample range, or %NULL.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2, if missing values are 
 * found among the variables given in @list.  Also checks for missing 
 * values within the adjusted sample range.  If missing values are 
 * encountered there, either (a) flag an error (if @misst != %NULL), 
 * or (b) if @pmod != %NULL, construct a "missing mask" and attach it 
 * to @pmod.
 * 
 * Returns: if @misst is not %NULL, either the ID number of the
 * variable for which a missing value is first found inside the
 * adjusted sample range or 0 if there is no such variable.  If
 * @pmod is not %NULL, returns 1 if there is an error creating
 * the missing obs mask, otherwise 0.  If both @misst and @pmod 
 * are %NULL, always returns 0.
 */

int adjust_t1t2 (MODEL *pmod, const int *list, int *t1, int *t2, 
		 int n, const double **Z, int *misst)
{
    int i, t, dwt = 0, t1min = *t1, t2max = *t2;
    int vi, missobs, ret = 0;
    int move_ends = 1;
    double xx;

    if (pmod != NULL && gretl_model_get_int(pmod, "wt_dummy")) {
	/* we have a weight variable which is a 0/1 dummy */
	dwt = pmod->nwt;
    }

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
	missobs = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == 0 || vi == LISTSEP) {
		continue;
	    }
	    xx = Z[vi][t];
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
	    vi = list[i];
	    if (vi == 0 || vi == LISTSEP) {
		continue;
	    }
	    xx = Z[vi][t];
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

    if (misst != NULL) {
	/* check for missing values within remaining range and
	   flag an error in case any are found */
	for (t=t1min; t<=t2max; t++) {
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (vi == 0 || vi == LISTSEP) {
		    continue;
		}
		xx = Z[vi][t];
		if (dwt) {
		    xx *= Z[dwt][t];
		}
		if (na(xx)) {
		    /* identify first missing obs and var */
		    *misst = t + 1;
		    ret = vi;
		    break;
		}
	    }
	    if (ret) {
		break;
	    }
	}     
    } else if (pmod != NULL) {
	/* construct a mask for missing values within remaining range? 
	   Note: we do this only if misst == NULL */
	missobs = 0;
	for (t=t1min; t<=t2max; t++) {
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (vi == 0 || vi == LISTSEP) {
		    continue;
		}
		xx = Z[vi][t];
		if (dwt) {
		    xx *= Z[dwt][t];
		}
		if (na(xx)) {
		    missobs++;
		    break;
		}
	    }
	}
	
	if (missobs == t2max - t1min + 1) {
	    /* no valid observations */
	    pmod->errcode = E_MISSDATA;
	    ret = 1;
	} else if (missobs > 0) {
#if 1
	    pmod->missmask = model_missmask(list, *t1, *t2, 
					    n, Z, dwt, NULL);
	    move_ends = 0;
#else	
	    pmod->missmask = model_missmask(list, t1min, t2max, 
					    n, Z, dwt, NULL);
#endif
	    if (pmod->missmask == NULL) {
		pmod->errcode = E_ALLOC;
		ret = 1;
	    }
	}
    } 

    if (move_ends) {
	*t1 = t1min; 
	*t2 = t2max;
    }

    return ret;
}

/* drop first/last observations from sample if missing obs 
   encountered -- also check for missing vals within the
   remaining sample: return non-zero if there are such.
   Adjust the t1 and t2 members of pdinfo if need be.
*/

int list_adjust_t1t2 (const int *list, const double **Z, DATAINFO *pdinfo)
{
    int err, misst = 0;

    err = adjust_t1t2(NULL, list, &pdinfo->t1, &pdinfo->t2, 
		      pdinfo->n, Z, &misst);
    if (err) {
	err = E_MISSDATA;
    }
    
    return err;
}

/* drop first/last observations from sample if missing obs 
   encountered -- also check for missing vals within the
   remaining sample */

int array_adjust_t1t2 (const double *x, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;

    for (t=t1min; t<t2max; t++) {
	if (na(x[t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(x[t])) t2max--;
	else break;
    }

    for (t=t1min; t<=t2max; t++) {
	if (na(x[t])) {
	    return t;
	}
    }

    *t1 = t1min; *t2 = t2max;

    return 0;
}

/**
 * varlist_adjust_sample: 
 * @list: list of variables to be tested for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @Z: data array.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2, if missing values are 
 * found among the variables given in @list at the start or end of
 * the range.  
 *
 * If you want to check for missing values inside the sample
 * range, use check_for_missing_obs() instead.
 * 
 * Returns: 1 if an adjustment was made, otherwise 0.
 */

int varlist_adjust_sample (const int *list, int *t1, int *t2, 
			   const double **Z)
{
    int oldt1 = *t1, oldt2 = *t2;
    int ret = 0;

    adjust_t1t2(NULL, list, t1, t2, 0, Z, NULL);

    if (*t1 != oldt1 || *t2 != oldt2) {
	ret = 1;
    }

    return ret;
}

/**
 * check_for_missing_obs: 
 * @list: list of variables to be tested for missing values.
 * @t1: on entry, intial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @Z: data array.
 * @misst: return location for index of the first missing observation
 *         inside the (possibly reduced) sample range, or %NULL.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2, if missing values are 
 * found among the variables given in @list.  Then checks for any
 * missing values within the adjusted range.  If such are found,
 * the return will be non-zero (see below).  In addition, if
 * @misst is non-%NULL it will receive the index number of the
 * observation where the first such missing value was found.
 * 
 * If you don't care about missing values inside the sample range,
 * use the simpler varlist_adjust_sample().
 *
 * Returns: the (non-zero) ID number of the first variable for 
 * which a missing value is first found inside the adjusted sample 
 * range or 0 if there is no such variable.
 */

int check_for_missing_obs (const int *list, int *t1, int *t2,
			   const double **Z, int *misst)
{
    int missv = 0;

    if (misst != NULL) {
	missv = adjust_t1t2(NULL, list, t1, t2, 0, Z, misst);
    } else {
	int tmp;

	missv = adjust_t1t2(NULL, list, t1, t2, 0, Z, &tmp);
    }

    return missv;
}

/* For handling the "omit" command, applied to a model that has
   missing values within the sample range.  The model as re-estimated
   with a reduced set of regressors must use the sample sample range
   as the original model, or else the comparison of the original and
   reduced models will be invalid.
*/

/* copy of a given model's missmask, or NULL */

static char *refmask;

/* Set the "reference" mask based on a target model */

void set_reference_missmask_from_model (const MODEL *pmod)
{
#if MASKDEBUG
    fprintf(stderr, "set_reference_missmask: using model = %p\n", 
	    (void *) pmod);
#endif
    if (pmod != NULL) {
	refmask = gretl_strdup(pmod->missmask);
    } 
}

int copy_to_reference_missmask (const char *mask)
{
    if (refmask != NULL) {
	free(refmask);
    }

    refmask = gretl_strdup(mask);

#if MASKDEBUG
    fprintf(stderr, "copy_to_reference_missmask: refmask = %p\n", 
	    (void *) refmask);
#endif

    return (refmask == NULL)? E_ALLOC : 0;
}

int set_reference_missmask_from_list (const int *list,
				      const double **Z,
				      const DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    char *mask;
    int nmiss = 0;
    int err = 0;

    mask = model_missmask(list, pdinfo->t1, pdinfo->t2, 
			  pdinfo->n, Z, 0, &nmiss);

    if (nmiss == T) {
	return E_DATA;
    } else if (nmiss == 0 && mask != NULL) {
	free(mask);
	mask = NULL;
    } else if (mask == NULL) {
	err = E_ALLOC;
    }

#if MASKDEBUG
    fprintf(stderr, "set_reference_missmask_from_list\n");
    printlist(list, "list");
    fprintf(stderr, "nmiss = %d, mask = %p\n", nmiss, (void *) mask);
#endif

    if (!err) {
	free(refmask);
	refmask = mask;
    }

    return err;
}

/* attach the reference missing obs mask to a specified model */

int apply_reference_missmask (MODEL *pmod)
{
    if (refmask != NULL) {
#if MASKDEBUG
	fprintf(stderr, "applying reference mask to model = %p\n", 
		(void *) pmod);
#endif
	pmod->missmask = refmask;
	refmask = NULL;
    }

    return 0;
}

int reference_missmask_present (void)
{
    return (refmask != NULL);
}

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
 * @list: list of variables to process, or an empty list or %NULL
 *        to process all variables.
 * @param: string representation of the numerical value to treat as missing.
 * @Z: data matrix.
 * @pdinfo: pointer to data information struct.
 * @prn: pointer to printing struct.
 * 
 * Set to "missing" each observation of each series in @list that
 * has the value represented by @param.
 *
 * Returns: 1 if at least one observation was set as missing,
 * otherwise 0.
 */

int set_miss (const int *list, const char *param, double **Z,
	      DATAINFO *pdinfo, PRN *prn)
{
    double missval = atof(param);
    int i, count, ret = 0;

    if (list == NULL || list[0] == 0) {
	count = real_setmiss(missval, 0, Z, pdinfo);
	if (count) { 
	    pprintf(prn, _("Set %d values to \"missing\"\n"), count);
	    ret = 1;
	} else {
	    pputs(prn, _("Didn't find any matching observations\n"));
	}
    } else {
	for (i=1; i<=list[0]; i++) {
	    if (var_is_scalar(pdinfo, list[i])) {
		pprintf(prn, _("The variable %s is a scalar\n"), 
			pdinfo->varname[list[i]]);
		continue;
	    }
	    count = real_setmiss(missval, list[i], Z, pdinfo);
	    if (count) { 
		pprintf(prn, _("%s: set %d observations to \"missing\"\n"), 
			pdinfo->varname[list[i]], count);
		ret = 1;
	    } else { 
		pprintf(prn, _("%s: Didn't find any matching observations\n"),
			pdinfo->varname[list[i]]);
	    }
	}
    }

    return ret;
}
