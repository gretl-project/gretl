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
 * @dset: dataset struct.
 *
 * Reverses the effect of repack_missing_daily_obs(), hence
 * restoring a daily dataset to its correct order.  The
 * model pointer, @pmod, must be the same as that given 
 * previously to repack_missing_daily_obs().
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int undo_daily_repack (MODEL *pmod, DATASET *dset)
{
    int i, j, t, m, g;
    double *tmpmiss = NULL, *tmpgood = NULL;
    /* note: pmod->t2 is still "shrunk" at this point */
    MISSOBS *mobs;
    int err = 0;

    if (!gretl_model_get_int(pmod, "daily_repack")) {
	return E_DATA;
    }

    mobs = (MISSOBS *) gretl_model_get_data(pmod, "missobs");
    if (mobs == NULL) {
	return E_DATA;
    }

    /* take charge of MISSOBS pointer */
    gretl_model_detach_data_item(pmod, "missobs");
			   
    tmpmiss = malloc((mobs->misscount + pmod->nobs) * sizeof *tmpmiss);
    if (tmpmiss == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    tmpgood = tmpmiss + mobs->misscount;

    for (j=1; j<=pmod->list[0]; j++) {
	i = pmod->list[j];

	if (i == 0 || i == LISTSEP) {
	    continue;
	}

	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	     tmpgood[g++] = dset->Z[i][t];
	}
	for (t=pmod->t2 + 1; t<=pmod->t2 + mobs->misscount; t++) {
	    tmpmiss[m++] = dset->Z[i][t];
	}
	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2 + mobs->misscount; t++) {
	    if (mobs->missvec[t] == '1') {
		dset->Z[i][t] = tmpmiss[m++];
	    } else {
		dset->Z[i][t] = tmpgood[g++];
	    }
	}
    }

 bailout:

    free(tmpmiss);

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
repack_missing (const MODEL *pmod, DATASET *dset, 
		const char *missvec, int misscount)
{
    int i, v, t, m, g;
    double *tmpmiss, *tmpgood;
    int modn = pmod->t2 - pmod->t1 + 1;

    tmpmiss = malloc(modn * sizeof *tmpmiss);
    if (tmpmiss == NULL) {
	return E_ALLOC;
    }

    tmpgood = tmpmiss + misscount;

    for (i=1; i<=pmod->list[0]; i++) {
	v = pmod->list[i];

	if (v == 0 || v == LISTSEP) {
	    continue;
	}

	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (missvec[t] == '1') {
		tmpmiss[m++] = dset->Z[v][t];
	    } else {
		tmpgood[g++] = dset->Z[v][t];
	    }
	}

	m = g = 0;
	for (t=pmod->t1; t<=pmod->t2 - misscount; t++) {
	     dset->Z[v][t] = tmpgood[g++];
	}
	for (t=pmod->t2 + 1 - misscount; t<=pmod->t2; t++) {
	    dset->Z[v][t] = tmpmiss[m++];
	}
    }

    free(tmpmiss);

    return 0;
}

/**
 * repack_missing_daily_obs:
 * @pmod: pointer to model.
 * @dset: dataset struct.
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

int repack_missing_daily_obs (MODEL *pmod, DATASET *dset)
{
    char *missvec;
    int misscount;
    MISSOBS *mobs;
    int err = 0;

    missvec = model_missmask(pmod->list, pmod->t1, pmod->t2, 
			     dset->n, (const double **) dset->Z, 
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
	
    err = repack_missing(pmod, dset, missvec, misscount);

    if (err) {
	pmod->errcode = err;
	free(missvec);
    } else {
	gretl_model_set_int(pmod, "daily_repack", 1);
	/* Be sure to undo this after estimation! */
	pmod->t2 -= misscount;
	mobs->missvec = missvec;
	mobs->misscount = misscount;
	err = gretl_model_set_data(pmod, "missobs", mobs, 
				   GRETL_TYPE_STRUCT,
				   sizeof *mobs);
    }

    return err;
}

#define MASKDEBUG 0

/**
 * model_has_missing_obs:
 * @pmod: pointer to model.
 *
 * Returns: 1 if there are missing observations in the
 * model's sample range, otherwise 0.
 */

int model_has_missing_obs (const MODEL *pmod)
{
    int sample = pmod->t2 - pmod->t1 + 1;

    return (pmod->nobs < sample);
}

static int really_missing (int v, int t, const double **Z, int d)
{
    if (d > 0 && Z[d][t] == 0) {
	/* the obs is dummied out */
	return 0;
    } else {
	return na(Z[v][t]);
    } 
}

/* Construct a mask to code for missing observations within the sample
   range for a model.  The mask is an array of char with '0's 
   for OK observations and '1's for missing obs, terminated with
   a NUL byte. This format allows copying of the mask using
   strcpy.
*/

static char *model_missmask (const int *list, int t1, int t2,
			     int n, const double **Z, int dwt,
			     int *misscount)
{
    char *mask;
    int i, vi, t;

    mask = malloc(n + 1);
    if (mask == NULL) {
	return NULL;
    }

    memset(mask, '0', n);
    mask[n] = 0; /* note NUL-termination */

#if MASKDEBUG
    fprintf(stderr, "model_missmask: using series length %d\n", n);
#endif

    if (misscount != NULL) {
	*misscount = 0;
    }

    for (t=t1; t<=t2; t++) {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (really_missing(vi, t, Z, dwt)) {
#if MASKDEBUG > 1
		    fprintf(stderr, "model_missmask: NA at list[%d] (%d), "
			    "obs %d\n", i, vi, t);
#endif
		    /* FIXME dwt case and nobs? */
		    mask[t] = '1';
		    if (misscount != NULL) {
			*misscount += 1;
		    }
		    break;
		}
	    }
	}
    }

    return mask;
}

/**
 * model_adjust_sample: 
 * @pmod: pointer to gretl model.
 * @n: full length of data array.
 * @Z: data array.
 * @misst: location to receive the first observation with a
 *         missing value inside the sample range, or NULL.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in the t1 and t2 members of
 * @pmod, if missing values are found among the variables given 
 * in the list member of @pmod.
 *
 * Also checks for missing values within the adjusted sample range.  
 * If such values are encountered, either flag an error (if 
 * @misst != %NULL), or construct a "missing mask" and attach it 
 * to @pmod.
 *
 * Returns: if @misst is not NULL, either the ID number of the
 * variable for which a missing value is first found inside the
 * adjusted sample range or 0 if there is no such variable.  If
 * @misst is NULL, returns 1 if there is an error creating
 * the missing obs mask, otherwise 0.
 */

int model_adjust_sample (MODEL *pmod, int n, const double **Z, 
			 int *misst)
{
    int i, t, dwt = 0, t1min = pmod->t1, t2max = pmod->t2;
    int vi, missobs, ret = 0;
    int move_ends = 1;

    if (gretl_model_get_int(pmod, "wt_dummy")) {
	/* we have a weight variable which is a 0/1 dummy */
	dwt = pmod->nwt;
    }

    /* advance start of sample range to skip missing obs */
    for (t=t1min; t<t2max; t++) {
	missobs = 0;
	for (i=1; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (really_missing(vi, t, Z, dwt)) {
		    missobs = 1;
		    break;
		}
	    }
	}
	if (missobs) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip missing obs */
    for (t=t2max; t>t1min; t--) {
	missobs = 0;
	for (i=1; i<=pmod->list[0]; i++) {
	    vi = pmod->list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (really_missing(vi, t, Z, dwt)) {
		    missobs = 1;
		    break;
		}		
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
	for (t=t1min; t<=t2max && !ret; t++) {
	    for (i=1; i<=pmod->list[0]; i++) {
		vi = pmod->list[i];
		if (vi > 0 && vi != LISTSEP) {
		    if (really_missing(vi, t, Z, dwt)) {
			/* identify first missing obs and var */
			*misst = t + 1;
			ret = vi;
			break;
		    }		    
		}
	    }
	}     
    } else {
	/* construct a mask for missing values within remaining range? 
	   we do this only if misst == NULL */
	missobs = 0;
	for (t=t1min; t<=t2max; t++) {
	    for (i=1; i<=pmod->list[0]; i++) {
		vi = pmod->list[i];
		if (vi > 0 || vi != LISTSEP) {
		    if (really_missing(vi, t, Z, dwt)) {
			missobs++;
			break;
		    }
		}
	    }
	}
	
	if (missobs == t2max - t1min + 1) {
	    /* no valid observations */
	    pmod->errcode = E_MISSDATA;
	    ret = 1;
	} else if (missobs > 0) {
	    pmod->missmask = model_missmask(pmod->list, pmod->t1, pmod->t2, 
					    n, Z, dwt, NULL);
	    move_ends = 0;
	    if (pmod->missmask == NULL) {
		pmod->errcode = E_ALLOC;
		ret = 1;
	    }
	}
    } 

    if (move_ends) {
	pmod->t1 = t1min; 
	pmod->t2 = t2max;
    }

#if MASKDEBUG
    if (pmod->missmask != NULL) {
	fprintf(stderr, "model at %p: now has mask at %p\n",
		(void *) pmod, (void *) pmod->missmask);
    }
#endif

    return ret;
}

/**
 * first_missing_index:
 * @x: array to be checked for missing values.
 * @t1: start of range to check.
 * @t2: end of range to check.
 *
 * Returns: the index of the first missing observation in @x 
 * over the sample range @t1 to @t2, or -1 if there is no
 * such observation.
 */

int first_missing_index (const double *x, int t1, int t2)
{
    int t;

    for (t=t1; t<=t2; t++) {
	if (na(x[t])) {
	    return t;
	}
    }

    return -1;
}

/**
 * series_adjust_sample: 
 * @x: series to be checked for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 *
 * Adjusts @t1 and @t2 so as to drop any leading or trailing 
 * missing observations. 
 *
 * Returns: E_MISSDATA if interior missing values were found
 * within the (possibly adjusted) sample range, otherwise 0.
 */

int series_adjust_sample (const double *x, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;
    int err = 0;

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
	    err = E_MISSDATA;
	    break;
	}
    }

    *t1 = t1min; 
    *t2 = t2max;

    return err;
}

/**
 * list_adjust_sample: 
 * @list: list of variables to be tested for missing values.
 * @t1: on entry, initial start of sample range; on exit,
 *      start of sample range adjusted for missing values.
 * @t2: on entry, initial end of sample range; on exit, end
 *      of sample range adjusted for missing values.
 * @dset: dataset struct.
 * @nmiss: location to receive number of missing values within
 * (possibly adjusted) sample range.
 *
 * Drops leading or trailing observations from the sample range
 * initially given by the values in @t1 and @t2 if missing values
 * are found for any of the variables given in @list.
 *
 * If @nmiss is non-NULL it receives the number of missing values
 * inside the (possibly reduced) sample range, otherwise it is
 * considered an error if there are any such missing values.
 *
 * Returns: 0 on success or %E_MISSDATA or error.
 */

int list_adjust_sample (const int *list, int *t1, int *t2, 
			const DATASET *dset, int *nmiss)
{
    int i, t, t1min = *t1, t2max = *t2;
    int vi, missing, err = 0;

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
	missing = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    missing = 1;
		    break;
		}
	    }
	}
	if (missing) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
	missing = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    missing = 1;
		    break;
		}
	    }
	}
	if (missing) {
	    t2max--;
	} else {
	    break;
	}	
    }

    if (nmiss != NULL) {
	*nmiss = 0;
    }

    /* check for missing values within remaining range */
    for (t=t1min; t<=t2max && !err; t++) {
	missing = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && vi != LISTSEP) {
		if (na(dset->Z[vi][t])) {
		    if (nmiss == NULL) {
			err = E_MISSDATA;
		    } else {
			*nmiss += 1;
		    }
		    break;
		}
	    }
	}
    } 

    *t1 = t1min; 
    *t2 = t2max;

    return err;
}

/* For handling the "omit" command, applied to a model that has
   missing values within the sample range.  The model as re-estimated
   with a reduced set of regressors must use the sample sample range
   as the original model, or else the comparison of the original and
   reduced models will be invalid.
*/

/* copy of a given model's missmask, or NULL */

static char *refmask;

static char *build_refmask_from_model (const MODEL *pmod)
{
    int t, n = pmod->full_n;
    char *mask = malloc(n + 1);

    if (mask != NULL) {
	memset(mask, '0', n);
	mask[n] = 0; /* NUL-terminate */
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (na(pmod->uhat[t])) {
		mask[t] = '1';
	    }
	}
    }

    return mask;
}

/* Set the "reference" mask based on a given model */

void set_reference_missmask_from_model (const MODEL *pmod)
{
    if (pmod != NULL) {
	if (refmask != NULL) {
	    free(refmask);
	    refmask = NULL;
	}
	if (pmod->missmask != NULL) {
	    refmask = gretl_strdup(pmod->missmask);
	} else if (model_has_missing_obs(pmod)) {
	    refmask = build_refmask_from_model(pmod);
	}
    } 
}

int copy_to_reference_missmask (const char *mask)
{
    if (refmask != NULL) {
	free(refmask);
	refmask = NULL;
    }

    refmask = gretl_strdup(mask);

#if MASKDEBUG
    fprintf(stderr, "copy_to_reference_missmask: refmask = %p\n", 
	    (void *) refmask);
#endif

    return (refmask == NULL)? E_ALLOC : 0;
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
	refmask = NULL; /* note: apply it only once */
    }

    return 0;
}

int reference_missmask_present (void)
{
    return (refmask != NULL);
}

static int real_setmiss (double missval, int varno, 
			 DATASET *dset) 
{
    int i, t, count = 0;
    int start = 1, end = dset->v;

    if (varno) {
	start = varno;
	end = varno + 1;
    }

    for (i=start; i<end; i++) {
	for (t=0; t<dset->n; t++) {
	    if (dset->Z[i][t] == missval) {
		dset->Z[i][t] = NADBL;
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
 * @dset: dataset struct.
 * @prn: pointer to printing struct.
 * 
 * Set to "missing" each observation of each series in @list that
 * has the value represented by @param.
 *
 * Returns: 1 if at least one observation was set as missing,
 * otherwise 0.
 */

int set_miss (const int *list, const char *param, 
	      DATASET *dset, PRN *prn)
{
    double missval = atof(param);
    int i, vi, count, ret = 0;

    if (list == NULL || list[0] == 0) {
	count = real_setmiss(missval, 0, dset);
	if (count) { 
	    pprintf(prn, _("Set %d values to \"missing\"\n"), count);
	    ret = 1;
	} else {
	    pputs(prn, _("Didn't find any matching observations\n"));
	}
    } else {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    count = real_setmiss(missval, vi, dset);
	    if (count) { 
		pprintf(prn, _("%s: set %d observations to \"missing\"\n"), 
			dset->varname[vi], count);
		ret = 1;
	    } else { 
		pprintf(prn, _("%s: Didn't find any matching observations\n"),
			dset->varname[vi]);
	    }
	}
    }

    return ret;
}

/**
 * missing_obs_fraction:
 * @dset: dataset struct.
 * 
 * Returns: the fraction of the observations in @Z for which
 * all variables have missing values (empty rows).
 */

double missing_obs_fraction (const DATASET *dset)
{
    int missrow, totmiss = 0;
    int i, t;

    if (dset->n == 0) {
	return 0.0;
    }

    for (t=0; t<dset->n; t++) {
	missrow = 1;
	for (i=1; i<dset->v; i++) {
	    if (!na(dset->Z[i][t])) {
		missrow = 0;
		break;
	    }
	}
	totmiss += missrow;
    }

    return (double) totmiss / dset->n;
}

/**
 * any_missing_user_values:
 * @dset: dataset struct.
 * 
 * Returns: 1 if there are missing values for any non-hidden
 * variables within the current sample range, otherwise 0.
 */

int any_missing_user_values (const DATASET *dset)
{
    int i, t;

    if (dset->n == 0) {
	return 0;
    }

    for (i=1; i<dset->v; i++) {
	if (!series_is_hidden(dset, i)) {
	    for (t=dset->t1; t<=dset->t2; t++) {
		if (na(dset->Z[i][t])) {
		    return 1;
		}
	    }
	}
    }

    return 0;
}
