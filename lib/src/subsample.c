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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_panel.h"
#include "cmd_private.h"
#include "dbread.h"
#include "uservar.h"
#include "gretl_xml.h"

#define SUBDEBUG 0
#define FULLDEBUG 0

/*
  The purpose of the static pointers below: When the user subsamples
  the current dataset in a non-trivial way -- i.e., by selecting cases
  rather than just moving the starting or ending points of the data
  range -- we create a new sub-dataset, and we need to keep the full
  dataset around so that it can be restored later.  The pointer
  @fullset is used to record the address of the full dataset.

  In addition, @peerset keeps track of the location of the DATASET
  struct associated with the backed-up full dataset; by means of this,
  we can know when to free the full dataset and when not to (for
  instance, if we're freeing an auxiliary dataset).
*/

static DATASET *fullset;
static DATASET *peerset;

#define SUBMASK_SENTINEL 127

static int smpl_get_int (const char *s, DATASET *dset, int *err);

/* accessors for full dataset, when sub-sampled */

DATASET *fetch_full_dataset (void)
{
    return fullset;
}

static int full_data_length (const DATASET *dset)
{
    int n = 0;

    if (fullset != NULL) {
	n = fullset->n;
    } else if (dset != NULL) {
	n = dset->n;
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

char *copy_subsample_mask (const char *src, int *err)
{
    char *ret = NULL;

    if (src == RESAMPLED) {
	ret = RESAMPLED;
    } else if (src != NULL) {
	int n = get_submask_length(src);

	ret = malloc(n * sizeof *ret);
	if (ret != NULL) {
	    memcpy(ret, src, n);
	} else {
	    *err = E_ALLOC;
	}
    }

    return ret;
}

char *copy_datainfo_submask (const DATASET *dset, int *err)
{
    char *mask = NULL;

    if (complex_subsampled()) {
	mask = copy_subsample_mask(dset->submask, err);
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

int write_datainfo_submask (const DATASET *dset, FILE *fp)
{
    int ret = 0;

    if (dset->submask == RESAMPLED) {
	unsigned int seed = get_resampling_seed();

	fprintf(fp, "<resample seed=\"%u\" n=\"%d\"/>\n", seed, dset->n);
	ret = 1;
    } else if (complex_subsampled()) {
	int i, n = get_submask_length(dset->submask);

	fprintf(fp, "<submask length=\"%d\">", n);
	for (i=0; i<n; i++) {
	    fprintf(fp, "%d ", (int) dset->submask[i]);
	}
	fputs("</submask>\n", fp);

	if (dset->restriction != NULL) {
	    gretl_xml_put_tagged_string("restriction", dset->restriction, fp);
	}

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

/* all values apart from the sentinel are initialized to zero; once
   the mask is used, 1s will indicate included observations and
   0s will indicate excluded observations
*/

static char *make_submask (int n)
{
    char *mask = calloc(n + 1, 1);
    
    if (mask != NULL) {
	mask[n] = SUBMASK_SENTINEL;
    }

    return mask;
}

void set_dataset_resampled (DATASET *dset)
{
    dset->submask = RESAMPLED;
}

int dataset_is_resampled (const DATASET *dset)
{
    return (dset != NULL && dset->submask == RESAMPLED);
}

void maybe_free_full_dataset (const DATASET *dset)
{
    if (dset == peerset) {
	if (fullset != NULL) {
#if SUBDEBUG
	    fprintf(stderr, "maybe_free_full_dataset: freeing fullset at %p (Z at %p)\n",
		    (void *) fullset, (void *) fullset->Z);
#endif
	    if (fullset->Z != NULL) {
		free_Z(fullset);
	    }
	    clear_datainfo(fullset, CLEAR_SUBSAMPLE);
	    free(fullset);
	    fullset = NULL;
	}
	peerset = NULL;
    }
}

/* we do this on "restore full sample" */

static void relink_to_full_dataset (DATASET *dset)
{
#if SUBDEBUG
    fprintf(stderr, "relink_to_full_dataset: fullset = %p (freeing and nulling)\n",
	    (void *) fullset);
    fprintf(stderr, "fullset: v = %d, n = %d\n", fullset->v, fullset->n);
#endif

    *dset = *fullset;
    free(fullset);
    fullset = NULL;
    peerset = NULL;
}

/* sync malloced elements of the fullset struct that might
   have been moved via realloc
*/

static void sync_datainfo_members (const DATASET *dset)
{
    if (fullset->v > dset->v) {
	int i;

#if FULLDEBUG
	fprintf(stderr, "*** sync_datainfo_members: fullset->v = %d but dset->v = %d\n",
		fullset->v, dset->v);
	fprintf(stderr, " deleting the last %d element(s) of fullZ\n", 
		fullset->v - dset->v);
#endif
	for (i=dset->v; i<fullset->v; i++) {
	    free(fullset->Z[i]);
	    fullset->Z[i] = NULL;
	}
	fullset->v = dset->v;
    }

    fullset->varname = dset->varname;
    fullset->varinfo = dset->varinfo;
    fullset->descrip = dset->descrip;
}

/* attach_subsample_to_model:
 * @pmod: model to which subsample should be attached.
 * @dset: pointer to current dataset info.
 *
 * If the dataset is currently subsampled, record the subsample
 * information with the model so that it can be retrieved later.
 * 
 * Returns: 0 if the recording is not needed, or on success; non-zero
 * error code failure.
 */

int attach_subsample_to_model (MODEL *pmod, const DATASET *dset)
{
    int err = 0;

#if SUBDEBUG
    fprintf(stderr, "attach_subsample_to_model: fullset = %p\n",
	    (void *) fullset);
#endif

    if (fullset != NULL) {
	/* sync, in case anything has moved */
	sync_datainfo_members(dset);

	if (pmod->submask != NULL) {
	    free_subsample_mask(pmod->submask);
	}

	pmod->submask = copy_subsample_mask(dset->submask, &err);
    }

    return err;
}

/* If series have been added to a resampled dataset, we can't
   bring these back to the "full" dataset, which may have a
   longer or shorter series length, and from which there is
   no definite mapping by row. So we just delete them. In
   this function we destroy their varnames and varinfo
   structures; the numerical arrays get deleted later.
*/

static int resample_sync_dataset (DATASET *dset)
{
    if (dset->v > fullset->v) {
	char **varname;
	VARINFO **varinfo;
	int i, nv = fullset->v;

	for (i=fullset->v; i<dset->v; i++) {
	    free(dset->varname[i]);
	    free(dset->varinfo[i]);
	}

	varname = realloc(dset->varname, nv * sizeof *varname);
	if (varname == NULL) {
	    return E_ALLOC;
	}
	dset->varname = varname;

	varinfo = realloc(dset->varinfo, nv * sizeof *varinfo);
	if (varinfo == NULL) {
	    return E_ALLOC;
	}
	dset->varinfo = varinfo;
    }

    /* sync */
    fullset->varname = dset->varname;
    fullset->varinfo = dset->varinfo;
    fullset->descrip = dset->descrip;

    return 0;
}

/* Apparatus for updating full dataset when restoring full sample
   after sub-sampling.  
*/

static void
update_full_data_values (const DATASET *dset)
{
    int i, s, t;

#if SUBDEBUG
    fprintf(stderr, "update_full_data_values: fullset->Z=%p, dset->Z=%p, dset=%p\n",
	    (void *) fullset->Z, (void *) dset->Z, (void *) dset);
#endif

    for (i=1; i<fullset->v && i<dset->v; i++) {
	s = 0;
	for (t=0; t<fullset->n; t++) {
	    if (dset->submask[t] == 1) {
		fullset->Z[i][t] = dset->Z[i][s++];
	    } else if (dset->submask[t] == 'p') {
		/* skip panel padding (?) */
		s++;
	    }
	}
    }
}

static int update_case_markers (const DATASET *dset)
{
    int err = 0;

    if (dset->markers && !fullset->markers) {
	dataset_allocate_obs_markers(fullset);
	if (fullset->S == NULL) {
	    err = 1;
	} else {
	    int t, subt = 0;

	    for (t=0; t<fullset->n; t++) {
		if (dset->submask[t]) {
		    strcpy(fullset->S[t], dset->S[subt++]);
		} else {
		    sprintf(fullset->S[t], "%d", t + 1);
		}
	    }
	}
    }	
	
    return err;
}

static int add_new_vars_to_full (DATASET *dset)
{
    int V1 = dset->v;
    int V0 = fullset->v;
    int N = fullset->n;
    double **newZ = NULL;
    int i, t, s;
    int err = 0;

    if (V1 <= V0) {
	return 0;
    }

    if (dset->submask == NULL) {
	return E_NOMERGE;
    }

#if SUBDEBUG
    fprintf(stderr, "add_new_vars_to_full:\n");
    fprintf(stderr, " V1 = dset->v = %d; V0 = fullset->v = %d\n",
	    V1, V0);
    fprintf(stderr, " dset->Z = %p, fullset->Z = %p\n", (void *) dset->Z, 
	    (void *) fullset->Z);
#endif

    /* allocate expanded data array */
    newZ = realloc(fullset->Z, V1 * sizeof *fullset->Z);

    if (newZ == NULL) {
	return E_ALLOC;
    } 

    fullset->Z = newZ;

    for (i=V0; i<dset->v && !err; i++) {
#if FULLDEBUG
	fprintf(stderr, "adding to full: var %d (%s, level %d)\n",
		i, dset->varname[i], series_get_stack_level(dset, i));
#endif
	fullset->Z[i] = malloc(N * sizeof **fullset->Z);
	if (fullset->Z[i] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (err) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "After expansion, fullset->Z = %p (%d vars)\n", 
	    (void *) fullset->Z, dset->v);
#endif

    for (i=V0; i<dset->v; i++) {
	s = 0;
	for (t=0; t<N; t++) {
	    fullset->Z[i][t] = (dset->submask[t])? 
		dset->Z[i][s++] : NADBL;
	}
    }

    fullset->v = V1;

    return 0;
}

static int sync_data_to_full (DATASET *dset)
{
    int err;

    /* update values for pre-existing series, which may have been
       modified via "genr" etc. */
    update_full_data_values(dset);

    /* if case markers were added when subsampled, carry them back */
    update_case_markers(dset);

    /* delete any newly added hidden vars */
    err = dataset_destroy_hidden_variables(dset, fullset->v);

    /* in case any new vars were added when subsampled, try to merge
       them into the full dataset */
    if (!err) {
	err = add_new_vars_to_full(dset);
    }

    return err;
}

static char *make_current_sample_mask (DATASET *dset)
{
    int n = full_data_length(dset);
    char *currmask = NULL;
    int s, t, err = 0;

    if (dset->submask == NULL) {
	/* no pre-existing mask: not sub-sampled */
	currmask = make_submask(n);
	if (currmask != NULL) {
	    for (t=dset->t1; t<=dset->t2; t++) {
		currmask[t] = 1;
	    }
	}
    } else {
	/* there's a pre-existing mask */
	currmask = copy_subsample_mask(dset->submask, &err);
	if (currmask != NULL) {
	    s = -1;
	    for (t=0; t<n; t++) {
		if (dset->submask[t]) s++;
		if (s < dset->t1 || s > dset->t2) {
		    currmask[t] = 0;
		} 
	    }
	}
    }

    return currmask;
}

/* Deal with the case where sampling has been done simply by
   moving the sample-range endpoints
*/

static int restore_full_easy (DATASET *dset, ExecState *state)
{
    int t1min, t2max;

    if (state == NULL) {
	/* not inside a function */
	t1min = 0;
	t2max = dset->n - 1;
    } else {
	/* don't go outside the bounds set on entry to
	   a function */
	sample_range_get_extrema(dset, &t1min, &t2max);
    }

    if (dset->t1 != t1min || dset->t2 != t2max) {
	dset->t1 = t1min;
	dset->t2 = t2max;
#if SUBDEBUG
	fprintf(stderr, "restore_full_sample: just set t1=%d and t2=%d\n",
		t1min, t2max);
#endif
    }

    return 0;
}

/* restore_full_sample: 
 * @dset: dataset struct.
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

int restore_full_sample (DATASET *dset, ExecState *state)
{
    int err = 0;

    if (dset == NULL) {
	return E_NODATA;
    } else if (!complex_subsampled()) {
	return restore_full_easy(dset, state);
    }

#if FULLDEBUG || SUBDEBUG
    fprintf(stderr, "\nrestore_full_sample: dset=%p, state=%p, fullset=%p\n", 
	    (void *) dset, (void *) state, (void *) fullset);
#endif

    /* Beyond this point we are doing a non-trivial restoration
       of a stored "full" dataset which has previously been
       subsampled, e.g., by some boolean criterion.
    */

    if (dset->submask == RESAMPLED) {
	err = resample_sync_dataset(dset);
    } else {
	if (dset->padmask != NULL) {
	    fprintf(stderr, "restore_full_sample: first undo panel padding\n");
	    err = undo_panel_padding(dset);
	}
	if (!err) {
	    sync_datainfo_members(dset);
	    err = sync_data_to_full(dset);
	}
    }

    if (err == E_NOMERGE) {
        gretl_errmsg_set(_("Missing sub-sample information; can't merge data\n"));
    }

    if (err) {
	return err;
    }

    /* destroy sub-sampled data array */
#if SUBDEBUG
    fprintf(stderr, "restore_full_sample: freeing sub-sampled Z at %p (v = %d, n = %d)\n"
	    " and clearing dset at %p\n", (void *) dset->Z, dset->v, dset->n,
	    (void *) dset);
#endif
    free_Z(dset);
    clear_datainfo(dset, CLEAR_SUBSAMPLE);

    relink_to_full_dataset(dset);

    if (state != NULL) {
	/* in this case restoring the "full" sample really means, relative 
	   to the original state 
	*/
	if (state->submask != NULL) {
	    err = restrict_sample_from_mask(state->submask, dset,
					    OPT_NONE);
	} else {
	    int t1min, t2max;

	    sample_range_get_extrema(dset, &t1min, &t2max);
	    if (dset->t1 < t1min) {
		dset->t1 = t1min;
	    }
	    if (dset->t2 > t2max) {
		dset->t2 = t2max;
	    }
	}
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
make_missing_mask (const int *list, const DATASET *dset, char *mask)
{
    int i, vi, t;

    if (list != NULL && list[0] > 0) {
	/* check specified list of variables */
	for (t=0; t<dset->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<=list[0]; i++) {
		vi = list[i];
		if (na(dset->Z[vi][t])) {
		    mask[t] = 0;
		    break;
		}
	    }
	}
    } else {	
	/* check all variables */
	for (t=0; t<dset->n; t++) {
	    mask[t] = 1;
	    for (i=1; i<dset->v; i++) {
		if (!series_is_hidden(dset, i) && na(dset->Z[i][t])) {
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

static int mask_from_temp_dummy (const char *s, DATASET *dset, 
				 char *mask, PRN *prn)
{
    char formula[MAXLINE];
    double *x;
    int err = 0;

    *formula = '\0';
    strncat(formula, s, MAXLINE - 1);

    x = generate_series(formula, dset, prn, &err);

    if (!err) {
	err = copy_dummy_to_mask(mask, x, dset->n);
	if (err) {
	    gretl_errmsg_sprintf(_("'%s' is not a dummy variable"), "mask");
	}
    }

    free(x);

    return err;
}

static int mask_from_dummy (const char *s, const DATASET *dset,
			    char *mask)
{
    char dname[VNAMELEN] = {0};
    int dnum, err = 0;

    gretl_scan_varname(s, dname);
    dnum = series_index(dset, dname);

    if (dnum == dset->v) {
	gretl_errmsg_sprintf(_("Variable '%s' not defined"), dname);
	err = 1;
    } else {
	err = copy_dummy_to_mask(mask, dset->Z[dnum], dset->n);
	if (err) {
	    gretl_errmsg_sprintf(_("'%s' is not a dummy variable"), dname);
	}
    }

    return err;
}

/* how many observations are selected by the given 
   subsample mask? */

static int 
count_selected_cases (const char *mask, const DATASET *dset)
{
    int i, n = 0;

    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    n++;
	}
    }

    return n;
}

/* panel: how many distinct cross-sectional units are included 
   in the masked subset of observations? */

static int 
count_panel_units (const char *mask, const DATASET *dset)
{
    int u, ubak = -1;
    int i, n = 0;

    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    u = i / dset->pd;
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

/* construct mask for taking random sub-sample from dataset:
   we're selecting 'subn' cases without replacement, and
   there may or may not be an existing mask in place that
   has to be respected.
*/

static int make_random_mask (const char *s, const char *oldmask, 
			     DATASET *dset, char *mask)
{
    unsigned u;
    int targ, avail, oldn = dset->n;
    int i, subn, cases, rejn;
    int err = 0;

    /* how many cases are requested? */
    subn = smpl_get_int(s, dset, &err);

    if (subn <= 0 || subn >= dset->n) {
	err = 1;
    } else if (oldmask != NULL) {
	/* dataset is already sub-sampled */
	oldn = 0;
	for (i=0; i<dset->n; i++) {
	    if (oldmask[i]) {
		oldn++;
	    }
	}
	if (subn >= oldn) {
	    err = 1;
	}
    }	

    if (err) {
	gretl_errmsg_sprintf(_("Invalid number of cases %d"), subn);
	return err;
    }	

    /* Which is smaller: the number of cases to be selected or the
       complement, the number to be discarded?  For the sake of
       efficiency we'll go for the smaller value.
    */
    rejn = oldn - subn;
    if (rejn < subn) {
	/* select rejn observations to discard */
	targ = rejn;
	avail = 1;
    } else {
	/* select subn observations to include */
	targ = subn;
	avail = 0;
    }

    for (i=0; i<dset->n; i++) {
	if (oldmask != NULL && oldmask[i] == 0) {
	    /* obs is already excluded, hence not selectable */
	    mask[i] = -1;
	} else {
	    mask[i] = avail;
	}
    }

    cases = 0;

    while (cases < targ) {
	u = gretl_rand_int_max(dset->n);
	if (mask[u] == avail) {
	    /* obs is available, not yet selected */
	    mask[u] = !avail;
	    cases++;
	}
    }

    if (oldmask != NULL) {
	/* undo the 'already excluded' coding above */
	for (i=0; i<dset->n; i++) {
	    if (mask[i] == -1) {
		mask[i] = 0;
	    }
	}
    }

    return err;
}

int backup_full_dataset (DATASET *dset)
{
#if SUBDEBUG
    int newfull = 0;
#endif

    if (fullset == NULL) {
	fullset = malloc(sizeof *fullset);
	if (fullset == NULL) {
	    return E_ALLOC;
	}
#if SUBDEBUG
	newfull = 1;
#endif
    }

    if (dset != NULL) {
	*fullset = *dset;
	peerset = dset;
    } 

#if SUBDEBUG
    fprintf(stderr, "backup_full_dataset: fullset = %p (%s)\n",
	    (void *) fullset, newfull ? "new" : "old");
#endif

    return 0;
}

int complex_subsampled (void)
{
    return (fullset != NULL && fullset->Z != NULL);
}

int get_full_length_n (void)
{
    return (fullset != NULL) ? fullset->n : 0;
}

/* When sub-sampling on some boolean criterion, check to see if we can
   meet the criterion by simply adjusting the endpoints of the sample
   range: life will be simpler if that is so.
*/

static int mask_contiguous (const char *mask,
			    const DATASET *dset,
			    int *pt1, int *pt2)
{
    int t, t1 = 0, t2 = dset->n - 1;
    int contig = 1;

    for (t=0; t<dset->n; t++) {
	if (mask[t] == 0) {
	    t1++;
	} else {
	    break;
	}
    }

    for (t=dset->n - 1; t>=0; t--) {
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

    if (contig && dataset_is_panel(dset)) {
	int n = t2 - t1 + 1;

	/* sample must leave a whole number of panel units; moreover,
	   to retain "panelness" this number must be greater than 1 
	*/
	if (t1 % dset->pd != 0 || n % dset->pd != 0) {
	    contig = 0;
	} else if (n == dset->pd) {
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
copy_data_to_subsample (DATASET *subset, const DATASET *dset,
			int maxv, const char *mask)
{
    int i, t, s;

#if SUBDEBUG
    fprintf(stderr, "copy_data_to_subsample: subset = %p, dset = %p\n",
	    (void *) subset, (void *) dset);
#endif

    /* copy data values */
    for (i=1; i<maxv; i++) {
	s = 0;
	for (t=0; t<dset->n; t++) {
	    if (mask == NULL) {
		subset->Z[i][s++] = dset->Z[i][t];
	    } else if (mask[t] == 1) {
		subset->Z[i][s++] = dset->Z[i][t];
	    } else if (mask[t] == 'p') {
		/* panel padding */
		subset->Z[i][s++] = NADBL;
	    }
	}
    }

    /* copy observation markers, if any */
    if (dset->markers && subset->markers) {
	s = 0;
	for (t=0; t<dset->n; t++) {
	    if (mask == NULL || mask[t] == 1 || mask[t] == 'p') {
		strcpy(subset->S[s++], dset->S[t]);
	    }
	}
    }

    /* copy panel time info? */
    if (dataset_is_panel(subset)) {
	if (subset->pd == dset->pd) {
	    subset->panel_pd = dset->panel_pd;
	    subset->panel_sd0 = dset->panel_sd0;
	}
    }

    strcpy(subset->stobs, "1");
    sprintf(subset->endobs, "%d", subset->n);
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

/* Copy list of panel periods from Ti to T0; return 1 if members of
   the list Ti are consecutive, 0 if they are not. 
*/

static int copy_periods_list (int *T0, const int *Ti)
{
    int j, ret = 1;

    for (j=0; j<=Ti[0] && ret; j++) {
	if (j > 1 && Ti[j] != Ti[j-1] + 1) {
	    ret = 0;
	} else {
	    T0[j] = Ti[j];
	}
    }

    return ret;
}

/* When sub-sampling panel data on some boolean criterion: see if the
   exclusion of certain rows leaves a balanced panel.  Note that the
   requirement is not simply that each unit has the same number of
   temporal observations -- they must have the _same_ observations,
   and in addition the observations must be temporally contiguous.  
*/

static int mask_leaves_balanced_panel (char *mask, const DATASET *dset)
{
    int T = dset->pd;
    int *T0, *Ti;
    int i, u, ubak = -1;
    int ret = 1;

    T0 = gretl_list_new(T);
    Ti = gretl_list_new(T);

    if (T0 == NULL || Ti == NULL) {
	free(T0);
	free(Ti);
	return 0;
    }

    T0[0] = Ti[0] = 0;

    for (i=0; i<dset->n && ret; i++) {
	if (mask[i]) {
	    u = i / T;
	    if (u != ubak) {
		if (Ti[0] > 0) {
		    if (T0[0] == 0) {
			/* we haven't made the reference list, T0, yet */
			ret = copy_periods_list(T0, Ti);
		    } else if (gretl_list_cmp(Ti, T0)) {
			/* the current list differs from the reference one */
			ret = 0;
		    }
		}
		Ti[0] = 1;
		Ti[1] = i % T;
	    } else {
		Ti[0] += 1;
		Ti[Ti[0]] = i % T;
	    }
	    ubak = u;
	}
    }

    free(T0);
    free(Ti);

    return ret;
}

static int 
make_panel_submask (char *mask, const DATASET *dset, int *err)
{
    int T = dset->pd;
    int N = dset->n / T;
    char *umask, *pmask;
    int i, np = 0;

    umask = calloc(N + T, 1);
    if (umask == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    pmask = umask + N;
 
    for (i=0; i<dset->n; i++) {
	if (mask[i]) {
	    umask[i / T] = 1;
	    pmask[i % T] = 1;
	}
    }

    for (i=0; i<dset->n; i++) {
	if (!mask[i]) {
	    if (umask[i / T] && pmask[i % T]) {
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

    return np;
}

#define needs_string_arg(m) (m == SUBSAMPLE_RANDOM || \
	                     m == SUBSAMPLE_USE_DUMMY || \
                             m == SUBSAMPLE_BOOLEAN)

static int 
make_restriction_mask (int mode, const char *s,
		       const int *list, DATASET *dset,
		       const char *oldmask, char **pmask,
		       PRN *prn)
{
    char *mask = NULL;
    int sn = 0, err = 0;

    if (needs_string_arg(mode) && (s == NULL || *s == '\0')) {
	return E_ARGS;
    }

    mask = make_submask(dset->n);
    if (mask == NULL) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "make_restriction_mask: oldmask = %p\n", (void *) oldmask);
#endif

    /* construct subsample mask in one of several possible ways */

    if (mode == SUBSAMPLE_DROP_MISSING) {   
	err = make_missing_mask(list, dset, mask);
    } else if (mode == SUBSAMPLE_RANDOM) {
	err = make_random_mask(s, oldmask, dset, mask);
    } else if (mode == SUBSAMPLE_USE_DUMMY) {
	err = mask_from_dummy(s, dset, mask);
    } else if (mode == SUBSAMPLE_BOOLEAN) {
	err = mask_from_temp_dummy(s, dset, mask, prn);
    } else {
	gretl_errmsg_set(_("Sub-sample command failed mysteriously"));
	err = 1;
    }

    if (err) {
	/* exit now on unrecoverable error */
	free(mask);
	return err;
    }

    /* cumulate sample restrictions, if appropriate */
    if (oldmask != NULL && mode != SUBSAMPLE_RANDOM) {
	sn = overlay_masks(mask, oldmask, dset->n);
    } else {
	sn = count_selected_cases(mask, dset);
    }

    /* does this policy lead to an empty sample, or no change in the
       sample, perchance? */

    if (sn == 0) {
	gretl_errmsg_set(_("No observations would be left!"));
	err = 1;
    } else if (sn == dset->n) {
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

static void finalize_panel_subset (DATASET *subset)
{
    int pdp = subset->pd;
    int den = 10.0;

    while ((pdp = pdp / 10)) {
	den *= 10;
    }

    subset->sd0 = 1.0 + 1.0 / den;
    ntodate(subset->stobs, 0, subset); 
    ntodate(subset->endobs, subset->n - 1, subset);
}    

/* This is also used elsewhere: 

   in gui2/session.c, when re-establishing a previously sub-sampled data 
   state on reopening a saved session

   in gretl_func.c, on exit from a user function when the dataset was 
   sub-sampled on entry to the function (and we need to re-establish the
   original sub-sample)
*/

int 
restrict_sample_from_mask (char *mask, DATASET *dset, gretlopt opt)
{
    DATASET *subset;
    int err = 0;

    if (mask == RESAMPLED) {
	fprintf(stderr, "restrict_sample_from_mask: got RESAMPLED!\n");
	return E_DATA;
    }

    if ((opt & OPT_B) && !dataset_is_panel(dset)) {
	gretl_errmsg_set(_("--balanced option is invalid: the (full) "
			   "dataset is not a panel"));
	return E_BADOPT;
    }

    subset = datainfo_new();
    if (subset == NULL) {
	return E_ALLOC;
    }

    subset->n = count_selected_cases(mask, dset);
    subset->v = dset->v;

#if SUBDEBUG
    fprintf(stderr, "restrict_sample_from_mask:\n new subset=%p, "
	    "%d obs in subsample vs %d in full dataset\n",
	    (void *) subset, subset->n, dset->n);
#endif

    if (dataset_is_panel(dset)) {
	/* are we able to reconstitute a panel? */
	int nunits = count_panel_units(mask, dset);
	int ok = 0, npad = 0;

	if (nunits > 1 && subset->n > nunits) {
	    /* there's some possibility of doing so */
	    if (opt & OPT_B) {
		/* add padding rows? only if this was requested */
		npad = make_panel_submask(mask, dset, &err);
		if (err) {
		    free(subset);
		    return err;
		}
		ok = 1;
	    } else {		
		ok = mask_leaves_balanced_panel(mask, dset);
	    } 
	    if (ok) {
		subset->structure = STACKED_TIME_SERIES;
		subset->n += npad;
		subset->pd = subset->n / nunits;
		finalize_panel_subset(subset);
	    }
	} else if (nunits == 1 && subset->n == dset->pd) {
	    /* time series for single panel unit */
	    subset->structure = SPECIAL_TIME_SERIES;
	}
    }

    /* set up the sub-sampled datainfo */
    err = start_new_Z(subset, OPT_R);
    if (err) { 
	free(subset);
	return err;
    }

#if SUBDEBUG
    fprintf(stderr, "started new Z for subset (v=%d, n=%d, Z=%p)\n", 
	    subset->v, subset->n, (void *) subset->Z);
#endif

    /* link (don't copy) varnames and descriptions, since these are
       not dependent on the series length */
    subset->varname = dset->varname;
    subset->varinfo = dset->varinfo;
    subset->descrip = dset->descrip;

    /* set up case markers? */
    if (dset->markers) {
	err = dataset_allocate_obs_markers(subset);
	if (err) {
	    free_Z(subset);
	    free(subset);
	    return E_ALLOC;
	}
    }

    /* copy across data (and case markers, if any) */
    copy_data_to_subsample(subset, dset, dset->v, mask);

    err = backup_full_dataset(dset);

    subset->submask = copy_subsample_mask(mask, &err);

    /* switch pointers */
    *dset = *subset;
    free(subset);

    return err;
}

/* Below: we do this "precompute" thing if the dataset is already
   subsampled and the user wants to compound the restriction with a
   new one of the form "obs=x" or "obs!=x".  The "obs" references may
   get out of whack if we restore the full dataset first, as we
   usually do.  For example, say the spec is "obs!=50" to exclude
   observation 50: presumably the user means to exclude the 50th
   observation in the current, subsampled dataset, which may not be
   the same as the 50th observation in the full dataset.
*/

static char *precompute_mask (const char *s, const char *oldmask,
			      DATASET *dset, PRN *prn, int *err)
{
    char *tmp = make_submask(dset->n);
    char *mask = NULL;
    int i, t;

#if SUBDEBUG
    fprintf(stderr, "restrict_sample: precomputing new mask\n");
#endif

    if (tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* fill out mask relative to current, restricted dataset */
    *err = mask_from_temp_dummy(s, dset, tmp, prn);

    if (!*err) {
	/* make blank full-length mask */
	mask = make_submask(fullset->n);
	if (mask == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	/* map the new mask onto full data range */
	i = 0;
	for (t=0; t<fullset->n; t++) {
	    if (oldmask[t]) {
		mask[t] = tmp[i++];
	    }
	}
    }

    free(tmp);

    return mask;
}

/* Intended for time series data: trim any missing values
   at the start and end of the current sample range, then
   check the remaining range for missing values and flag
   an error if any are found.
*/

static int set_contiguous_sample (const int *list,
				  DATASET *dset,
				  gretlopt opt)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int err = 0;

    if (opt & OPT_P) {
	/* can't combine this with the "replace" option */
	return E_BADOPT;
    }

    if (list != NULL && list[0] > 0) {
	err = list_adjust_sample(list, &dset->t1, &dset->t2, dset, NULL);
    } else {
	int *biglist = NULL;
	int nvars = 0;

	biglist = full_var_list(dset, &nvars);
	if (nvars == 0) {
	    ; /* no-op */
	} else if (biglist == NULL) {
	    err = E_ALLOC;
	} else {
	    err = list_adjust_sample(biglist, &dset->t1, &dset->t2, 
				     dset, NULL);
	    free(biglist);
	}
    } 

    if (err) {
	dset->t1 = save_t1;
	dset->t2 = save_t2;
    }

    return err;
}

#if 1 /* let's be conservative here */

# define restriction_uses_obs(s) (strstr(s, "obs") != NULL)

#else

static int restriction_uses_obs (const char *s)
{
    return gretl_namechar_spn(s) == 3 && !strncmp(s, "obs", 3);
}

#endif

/* Make a string representing the sample restriction. This is
   for reporting purposes. Note that in some cases the
   incoming @restr string may be NULL, for instance if the
   no-missing option is chosen via the gretl GUI.
*/

static int make_restriction_string (DATASET *dset, char *old, 
				    const char *restr, int mode)
{
    char *s = NULL;
    int n = 0, err = 0;

    if (dset->restriction != NULL) {
	free(dset->restriction);
	dset->restriction = NULL;
    }

    if (old != NULL) {
	n = strlen(old);
    }

    if (mode == SUBSAMPLE_RANDOM) {
	n += strlen("random");
    } else if (mode == SUBSAMPLE_DROP_MISSING) {
	n += strlen("no-missing");
    } else if (restr != NULL) {
	n += strlen(restr);
    }

    if (n > 0) {
	s = malloc(n + 5);
	if (s == NULL) {
	    err = E_ALLOC;
	} 
    }

    if (s != NULL) {
	*s = '\0';
	if (old != NULL) {
	    strcpy(s, old);
	    strcat(s, " && ");
	}
	if (mode == SUBSAMPLE_RANDOM) {
	    strcat(s, "random");
	} else if (mode == SUBSAMPLE_DROP_MISSING) {
	    strcat(s, "no-missing");
	} else if (restr != NULL) {
	    strcat(s, restr);
	}
	dset->restriction = s;
    }

    if (old != NULL) {
	free(old);
    }

    return err;
}

/* restrict_sample: 
 * @line: command line (or %NULL).  
 * @dset: dataset struct.
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
 * In case OPT_M or OPT_C a @list of variables may be supplied; in 
 * cases OPT_O, OPT_R and OPT_N, @line must contain specifics.
 *
 * In case OPT_P is included, the restriction will rePlace any
 * existing sample restriction, otherwise the resulting restriction
 * will be the logical product of the new restriction and any
 * existing restriction.
 *
 * In case the original dataset was a panel and OPT_B was given,
 * we'll pad with missing values if necessary, to try to reconstitute 
 * a balanced panel.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int restrict_sample (const char *line, const int *list,
		     DATASET *dset, ExecState *state, 
		     gretlopt opt, PRN *prn)
{
    char *oldrestr = NULL;
    char *oldmask = NULL;
    char *mask = NULL;
    int free_oldmask = 0;
    int mode, err = 0;
    
    if (dset == NULL || dset->Z == NULL) {
	return E_NODATA;
    }

    gretl_error_clear();

#if FULLDEBUG || SUBDEBUG
    fprintf(stderr, "\nrestrict_sample: '%s'\n", line);
    fprintf(stderr, " dset=%p, state=%p, fullset=%p\n", (void *) dset, 
	    (void *) state, (void *) fullset);
#endif

    if (line != NULL) {
	/* skip the command word */
	line += strcspn(line, " ");
	line += strspn(line, " ");
    }

    if (opt & OPT_C) {
	return set_contiguous_sample(list, dset, opt);
    }	

    mode = get_restriction_mode(opt);

    if (mode == SUBSAMPLE_UNKNOWN) {
	gretl_errmsg_set("Unrecognized sample command");
	return 1;
    } 

    if (!(opt & OPT_P)) {
	/* not replacing but cumulating any existing restrictions */
	oldmask = make_current_sample_mask(dset);
	if (oldmask == NULL) {
	    return E_ALLOC;
	}
	free_oldmask = 1;
	if (dset->restriction != NULL) {
	    oldrestr = gretl_strdup(dset->restriction);
	}
    } else if (state != NULL && state->submask != NULL) {
	/* subsampling within a function, with incoming
	   restriction recorded in state
	*/
	oldmask = state->submask;
    }

    if (mode == SUBSAMPLE_BOOLEAN && fullset != NULL && 
	oldmask != NULL && restriction_uses_obs(line)) {
	mask = precompute_mask(line, oldmask, dset, prn, &err);
    }

    /* restore the full data range, for housekeeping purposes */
    if (!err) {
	err = restore_full_sample(dset, NULL);
    }

    if (err) {
	return err;
    }

    if (mask == NULL) {
	/* not already handled by "precompute" above */
	err = make_restriction_mask(mode, line, list, dset, 
				    oldmask, &mask, prn);
    }

    if (!err && mask != NULL) {
	int t1 = 0, t2 = 0;
	int contig = 0;

	if (mode != SUBSAMPLE_RANDOM) {
	    contig = mask_contiguous(mask, dset, &t1, &t2);
	}

	if (contig) {
	    /* The specified subsample consists of contiguous
	       observations, so we'll just adjust the range, avoiding
	       the overhead of creating a parallel dataset.
	    */
#if SUBDEBUG
	    fprintf(stderr, "restrict sample: got contiguous range\n");
#endif
	    dset->t1 = t1;
	    dset->t2 = t2;
	} else {
#if SUBDEBUG
	    fprintf(stderr, "restrict sample: using mask\n");
#endif
	    err = restrict_sample_from_mask(mask, dset, opt);
	}
    }

    free(mask);

    if (free_oldmask) {
	free(oldmask);
    }

    make_restriction_string(dset, oldrestr, line, mode);

    return err;
}

enum {
    SMPL_T1,
    SMPL_T2
};

static int panel_round (const DATASET *dset, int t, int code)
{
    if (code == SMPL_T1) {
	while ((t + 1) % dset->pd != 1) {
	    t++;
	}
    } else {
	while ((t + 1) % dset->pd != 0) {
	    t--;
	}
    }

    return t;
}

static int smpl_get_int (const char *s, DATASET *dset, int *err)
{
    int k;

    if (integer_string(s)) {
	k = atoi(s);
    } else if (gretl_is_scalar(s)) {
	k = gretl_scalar_get_value(s, NULL);
    } else {
	k = (int) generate_scalar(s, dset, err);
    }

    return k;
}

static int get_sample_limit (char *s, DATASET *dset, int code)
{
    int ret = -1;
    int err = 0;

#if SUBDEBUG
    fprintf(stderr, "get_sample_limit: s = '%s'\n", s);
#endif

    if (*s == '-' || *s == '+') {
	/* increment/decrement form */
	int incr = smpl_get_int(s + 1, dset, &err);

	if (!err) {
	    if (*s == '-') {
		incr = -incr;
	    }
	    if (dataset_is_panel(dset)) {
		incr *= dset->pd;
	    }
	    if (code == SMPL_T1) {
		ret = dset->t1 + incr;
	    } else {
		ret = dset->t2 + incr;
	    }
	}
    } else {
	/* absolute form */
	ret = get_t_from_obs_string(s, dset);
	if (ret < 0) {
	    gretl_error_clear();
	    ret = smpl_get_int(s, dset, &err);
	    /* convert to base 0 */
	    if (!err) ret--;
	}
	if (ret >= 0 && dataset_is_panel(dset)) {
	    ret = panel_round(dset, ret, code);
	}
    }

    return ret;
}

/* Catch the case where we're in a function and attempting to
   move t1 or t2 out of the range established on entry to the
   function: we don't want to carry forward any error message
   that implies t was out of the full data range.
*/

static void maybe_clear_range_error (int t, DATASET *dset)
{
    if (t >= 0 && t < dset->n) {
	gretl_error_clear();
    }
}

int set_sample (const char *line, DATASET *dset)
{
    int nf, new_t1 = dset->t1, new_t2 = dset->t2;
    int tmin = 0, tmax = 0;
    char newstart[64];
    char newstop[64];

    if (dset == NULL) {
	return E_NODATA;
    }

    gretl_error_clear();

    line += strcspn(line, " ");
    line += strspn(line, " ");

    nf = count_fields(line, NULL);

#if SUBDEBUG
    fprintf(stderr, "set_sample: line='%s', nf=%d, dset=%p\n", 
	    line, nf, (void *) dset);
    if (dset != NULL) {
	fprintf(stderr, "dset->v = %d, dset->n = %d, pd = %d\n",
		dset->v, dset->n, dset->pd);
    }
#endif

    if (nf == 2 && dset->n == 0) {
	/* database special */
	return db_set_sample(line, dset);
    }

    if (nf == 0) {
	/* no-op, just print the current sample */
	return 0;
    }

    sample_range_get_extrema(dset, &tmin, &tmax);

#if SUBDEBUG
    fprintf(stderr, "sample extrema: lo = %d, hi = %d\n", tmin, tmax);
#endif
	
    if (nf == 1) {
	if (sscanf(line, "%63s", newstart) != 1) {
	    gretl_errmsg_set(_("error reading smpl line"));
	    return 1;
	} else {
	    new_t1 = get_sample_limit(newstart, dset, SMPL_T1);
	    if (new_t1 < tmin || new_t1 > tmax) {
		maybe_clear_range_error(new_t1, dset);
		gretl_errmsg_set(_("error in new starting obs"));
		return 1;
	    }
	    dset->t1 = new_t1;
	    return 0;
	}
    }

    /* now we're looking at nf = 2 (2 fields) case */

    if (sscanf(line, "%63s %63s", newstart, newstop) != 2) {
	gretl_errmsg_set(_("error reading smpl line"));
	return 1;
    }

    if (strcmp(newstart, ";")) {
	new_t1 = get_sample_limit(newstart, dset, SMPL_T1);
	if (new_t1 < tmin || new_t1 > tmax) {
	    maybe_clear_range_error(new_t1, dset);
	    gretl_errmsg_set(_("error in new starting obs"));
	    return 1;
	}	
    }

    if (strcmp(newstop, ";")) {
	new_t2 = get_sample_limit(newstop, dset, SMPL_T2);
	if (new_t2 < tmin || new_t2 > tmax) {
	    maybe_clear_range_error(new_t2, dset);
	    gretl_errmsg_set(_("error in new ending obs"));
	    return 1;
	}
    }

    if (new_t1 < tmin || new_t1 > new_t2) {
	gretl_error_clear();
	gretl_errmsg_set(_("Invalid null sample"));
	return 1;
    }

    dset->t1 = new_t1;
    dset->t2 = new_t2;

    return 0;
}

/**
 * count_missing_values:
 * @dset: dataset struct.
 * @opt: use %OPT_V for verbose operation, %OPT_A to 
 * examine all data.
 * @prn: printing struct.
 * @err: location to receive error code.
 *
 * Prints a count of missing values (if any) in the current
 * dataset over the currently defined sample range (or the 
 * entire data range if %OPT_A is given). If %OPT_V is given 
 * this includes a count of missing values at each observation; 
 * otherwise it just includes global and per-variable counts.
 *
 * Returns: 0 if no missing values are found (or on error),
 * otherwise the total number of missing values.
 */

int count_missing_values (const DATASET *dset, gretlopt opt, 
			  PRN *prn, int *err)
{
    int missval = 0, missobs = 0, totvals = 0, oldmiss = 0;
    int T, t1, t2;
    int *missvec;
    double missfrac;
    int i, t, tmiss;

    if (opt & OPT_A) {
	t1 = 0;
	t2 = dset->n - 1;
    } else {
	t1 = dset->t1;
	t2 = dset->t2;
    }

    T = t2 - t1 + 1;

    missvec = malloc(dset->v * sizeof missvec);

    if (missvec == NULL) {
	*err = E_ALLOC;
	return 0;
    }

    for (i=0; i<dset->v; i++) {
	missvec[i] = 0;
    }

    for (t=t1; t<=t2; t++) {
	tmiss = 0;
	for (i=1; i<dset->v; i++) {
	    if (series_is_hidden(dset, i)) {
		continue;
	    }
	    if (na(dset->Z[i][t])) {
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

	    if (opt & OPT_V) {
		/* verbose print by observation */
		if (dset->markers) { 
		    pprintf(prn, "%8s %4d %s\n", dset->S[t], tmiss,
			    _("missing values"));
		} else {
		    char tmp[OBSLEN];

		    ntodate(tmp, t, dset);
		    pprintf(prn, "%8s %4d %s\n", tmp, tmiss,
			    _("missing values"));
		}
	    }
	}
	oldmiss = missval;
    }

    missfrac = 100.0 * (double) missobs / T;

    pprintf(prn, _("\nNumber of observations (rows) with missing data "
		   "values = %d (%.2f%%)\n"), missobs, missfrac);

    pprintf(prn, _("Total number of missing data values = %d (%.2f%% "
	    "of total data values)\n"), missval, 
	    (100.0 * (double) missval / totvals));

    if (missvec[0] > 0) {
	pputc(prn, '\n');
	for (i=1; i<dset->v; i++) {
	    if (missvec[i] > 0) {
		missfrac = 100.0 * (double) missvec[i] / T;
		pprintf(prn, "%8s: %d %s (%.2f%%); %d %s (%.2f%%)\n", 
			dset->varname[i], missvec[i], _("missing values"),
			missfrac, T - missvec[i], _("valid values"),
			100.0 - missfrac);
	    }
	}
    }

    free(missvec);

    return missval;
}

static void copy_series_info (DATASET *dest, const DATASET *src,
			      int maxv)
{
    int i;

    for (i=1; i<maxv; i++) {
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

   By default we create a dataset containing series up to the
   highest numbered series associated with the model (regressand
   and regressors).  If OPT_F is given we include all series; if
   OPT_G is given we include only the constant.
*/

int add_dataset_to_model (MODEL *pmod, const DATASET *dset,
			  gretlopt opt)
{
    const DATASET *srcset;
    char *mask = NULL;
    int maxv, sn = 0;

    if (pmod->dataset != NULL) {
	/* FIXME? */
	destroy_dataset(pmod->dataset);
	pmod->dataset = NULL;
    }

#if SUBDEBUG
    fprintf(stderr, "add_dataset_to_model: fullset=%p, pmod->submask=%p\n",
	    (void *) fullset, (void *) pmod->submask);
#endif

    if (fullset != NULL) {
	sync_datainfo_members(dset);
	srcset = fullset;
    } else {
	srcset = dset;
    }

    if (pmod->submask == NULL) {
	/* no subsample info: pmod was estimated on the full dataset,
	   so we'll reconstruct the full dataset */
	sn = srcset->n;
    } else {
	/* pmod was estimated on a subsample, which has to 
	   be reconstructed */
	int t;

	mask = calloc(srcset->n, 1);
	if (mask == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<srcset->n; t++) {
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

    if (opt & OPT_F) {
	maxv = srcset->v;
    } else if (opt & OPT_G) {
	maxv = 1;
    } else {
	maxv = highest_numbered_var_in_model(pmod, dset) + 1;
    }

    /* allocate auxiliary dataset */
    pmod->dataset = create_auxiliary_dataset(maxv, sn, 0);
    if (pmod->dataset == NULL) {
	return E_ALLOC;
    }

#if SUBDEBUG
    fprintf(stderr, "pmod->dataset allocated at %p\n", 
	    (void *) pmod->dataset);
#endif

    /* copy across info on series */
    copy_series_info(pmod->dataset, srcset, maxv);

    /* copy across data */
    copy_data_to_subsample(pmod->dataset, srcset, maxv, mask);

    /* dataset characteristics such as pd: if we're rebuilding the
       full dataset copy these across; but if we're reconstructing a
       subsampled dataset these features from the full dataset will be
       wrong, and we stay with the simple defaults
    */
    if (pmod->submask == NULL) {
	pmod->dataset->pd = srcset->pd;
	pmod->dataset->sd0 = srcset->sd0;
	strcpy(pmod->dataset->stobs, srcset->stobs);
	strcpy(pmod->dataset->endobs, srcset->endobs);
	pmod->dataset->structure = srcset->structure;
    }

    free(mask);

#if SUBDEBUG
    fputs("add_subsampled_dataset_to_model: success\n", stderr);
#endif

    return 0;
}

static int submask_match (const char *s1, const char *s2, int n)
{
    int t;

    if (s1 == RESAMPLED || s2 == RESAMPLED) {
	return s1 == RESAMPLED && s2 == RESAMPLED;
    }

    for (t=0; t<n; t++) {
	if (s1[t] != s2[t]) return 0;
    }

    return 1;
}

/* check the subsample mask from a model against datainfo to 
   see if it may have been estimated on a different
   (subsampled) data set from the current one
*/

int model_sample_problem (MODEL *pmod, const DATASET *dset)
{
    int n = dset->n;
    int ret = 1;

    if (pmod->submask == NULL) {
	/* the model has no sub-sampling info recorded */
	if (dset->submask == NULL) {
	    /* data set is not sub-sampled either, OK */
	    ret = 0;
	} else {
	    fputs(I_("dataset is subsampled, model is not\n"), stderr);
	    gretl_errmsg_set(_("dataset is subsampled, model is not\n"));
	    ret = 1;
	}
    } else {
	/* the model does have sub-sampling info recorded */
	if (dset->submask == NULL) {
	    fputs(I_("model is subsampled, dataset is not\n"), stderr);
	    gretl_errmsg_set(_("model is subsampled, dataset is not\n"));
	    ret = 1;
	} else if (submask_match(dset->submask, pmod->submask, n)) {
	    /* the subsamples (model and current data set) agree, OK */
	    ret = 0;
	} else {
	    /* the subsamples differ */
	    gretl_errmsg_set(_("model and dataset subsamples not the same\n"));
	    ret = 1;
	}
    }

    return ret;
}

static void dataset_type_string (char *str, const DATASET *dset)
{
    if (dataset_is_time_series(dset)) {
	strcpy(str, _("time series"));
    } else if (dataset_is_panel(dset)) {
        strcpy(str, _("panel"));
    } else {
        strcpy(str, _("undated"));
    }
}

static void pd_string (char *str, const DATASET *dset)
{
    if (custom_time_series(dset)) {
	strcpy(str, _("special"));
    } else {
	switch (dset->pd) {
	case 1:
	    strcpy(str, _("annual")); break;
	case 4:
	    strcpy(str, _("quarterly")); break;
	case 12:
	    strcpy(str, _("monthly")); break;
	case 24:
	    strcpy(str, _("hourly")); break;
	case 52:
	    strcpy(str, _("weekly")); break;
	case 5:
	case 6:
	case 7:
	    strcpy(str, _("daily")); break;
	case 10:
	    strcpy(str, _("decennial")); break;
	default:
	    strcpy(str, _("unknown")); break;
	}
    }
}

void print_sample_obs (const DATASET *dset, PRN *prn)
{
    char d1[OBSLEN], d2[OBSLEN];

    ntodate(d1, dset->t1, dset);
    ntodate(d2, dset->t2, dset);

    pprintf(prn, "%s: %s - %s", _("Current sample"), d1, d2);
    pprintf(prn, " (n = %d)\n", dset->t2 - dset->t1 + 1);
}

void print_sample_status (const DATASET *dset, PRN *prn)
{
    char tmp[128];

    if (complex_subsampled()) {
	pprintf(prn, "%s\n\n", _("Full dataset"));
	dataset_type_string(tmp, fullset);
	pprintf(prn, "%s: %s\n", _("Type"), tmp);
	if (dataset_is_time_series(fullset)) {
	    pd_string(tmp, fullset);
	    pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
	} else if (dataset_is_panel(fullset)) {
	    int nu = fullset->n / fullset->pd;

	    pprintf(prn, "%s: %d\n", _("Number of cross-sectional units"), nu);
	    pprintf(prn, "%s: %d\n", _("Number of time periods"), fullset->pd);
	}
	pprintf(prn, "%s: %s - %s (n = %d)\n", _("Range"), 
		fullset->stobs, fullset->endobs, fullset->n);

	pprintf(prn, "\n%s\n", _("Subsampled data"));
	if (dset->restriction != NULL) {
	    pprintf(prn, "(%s: %s)\n\n", _("restriction"), dset->restriction);
	} else {
	    pputc(prn, '\n');
	}
    }	

    dataset_type_string(tmp, dset);
    pprintf(prn, "%s: %s\n", _("Type"), tmp);
    if (dataset_is_time_series(dset)) {
	pd_string(tmp, dset);
	pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
    } else if (dataset_is_panel(dset)) {
	int nu = dset->n / dset->pd;

	pprintf(prn, "%s: %d\n", _("Number of cross-sectional units"), nu);
	pprintf(prn, "%s: %d\n", _("Number of time periods"), dset->pd);
    }	
    pprintf(prn, "%s: %s - %s (n = %d)\n", _("Full range"), 
	    dset->stobs, dset->endobs, dset->n);
    print_sample_obs(dset, prn); 
}

/**
 * data_report:
 * @dset: data information struct.
 * @fname: filename for current datafile.
 * @prn: gretl printing struct.
 * 
 * Write out a summary of the content of the current data set.
 * 
 * Returns: 0 on successful completion, non-zero on error.
 * 
 */

int data_report (const DATASET *dset, const char *fname, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN], tmp[MAXLEN];
    char tstr[48];
    int i;

    ntodate(startdate, 0, dset);
    ntodate(enddate, dset->n - 1, dset);

    sprintf(tmp, _("Data file %s\nas of"), 
	    (*fname != '\0')? fname : _("(unsaved)"));

    print_time(tstr);
    pprintf(prn, "%s %s\n\n", tmp, tstr);

    if (dset->descrip != NULL && *dset->descrip != '\0') {
	pprintf(prn, "%s:\n\n", _("Description"));
	pputs(prn, dset->descrip);
	pputs(prn, "\n\n");
    }

    dataset_type_string(tmp, dset);
    pprintf(prn, "%s: %s\n", _("Type of data"), tmp);
    
    if (dataset_is_time_series(dset)) {
	pd_string(tmp, dset);
	pprintf(prn, "%s: %s\n", _("Frequency"), tmp);
    }	

    pprintf(prn, "%s: %s - %s (n = %d)\n\n", _("Range"),
	    startdate, enddate, dset->n);

    pprintf(prn, "%s:\n\n", _("Listing of variables"));

    for (i=1; i<dset->v; i++) {
	pprintf(prn, "%*s  %s\n", VNAMELEN, dset->varname[i], 
		series_get_label(dset, i));
    }

    return 0;
}


