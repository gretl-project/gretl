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

/* Mechanism for recording and retrieving the preferred set of lags
   for a given variable.  We save the preferences that are set, and
   destroy the stacked preferences when we switch data sets.
*/

#define LDEBUG 0

#include "gretl.h"
#include "selector.h"
#include "lagpref.h"

enum {
    LAGS_NONE,    /* no lags specified for variable */
    LAGS_MINMAX,  /* min and max lags given (consecutive) */
    LAGS_LIST,    /* list of lags specific given */
    LAGS_TMP      /* provisional list when working from model */
} SpecType;

typedef struct lagpref_ lagpref;

struct lagpref_ {
    int v;
    char context;
    char spectype;
    union lspec {
	int lminmax[2];
	int *laglist;
    } lspec;
};

static lagpref **lprefs;
static int n_prefs;

void destroy_lag_preferences (void)
{
    int i;
    
    for (i=0; i<n_prefs; i++) {
	if (lprefs[i]->spectype == LAGS_LIST) {
	    free(lprefs[i]->lspec.laglist);
	}
	free(lprefs[i]);
    } 

    free(lprefs);
    lprefs = NULL;
    n_prefs = 0;
}

static int add_lpref_to_stack (lagpref *lpref)
{
    lagpref **tmp;
    int err = 0;

    tmp = realloc(lprefs, (n_prefs + 1) * sizeof *lprefs);
    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	lprefs = tmp;
	lprefs[n_prefs] = lpref;
	n_prefs++;
    }

    return err;
}

static lagpref *lpref_new (int v, char context, int type)
{
    lagpref *lpref = malloc(sizeof *lpref);
    
    if (lpref == NULL) {
	return NULL;
    }

    lpref->v = v;
    lpref->context = context;

#if LDEBUG > 1
    fprintf(stderr, "lpref_new: added lpref for var %d\n", v);
#endif

    if (type == LAGS_LIST) {
	lpref->spectype = type;
	lpref->lspec.laglist = NULL;
    } else if (context == LAG_Y_X || context == LAG_Y_W) {
	/* by default: one lag, but not enabled */
	lpref->spectype = LAGS_MINMAX;
	lpref->lspec.lminmax[0] = 1;
	lpref->lspec.lminmax[1] = 1;
    } else {
	lpref->spectype = LAGS_NONE;
    }

    return lpref;
}

#if LDEBUG > 1
static void print_lpref (lagpref *lpref)
{
    if (lpref == NULL) {
	fprintf(stderr, "No lag reference recorded\n");
	return;
    }

    fprintf(stderr, "lagpref at %p: v = %d, context = %d\n",
	    (void *) lpref, lpref->v, lpref->context);

    if (lpref->spectype == LAGS_NONE) {
	fprintf(stderr, "type == LAGS_NONE (%d)\n", LAGS_NONE);
    } else if (lpref->spectype == LAGS_MINMAX) {
	fprintf(stderr, "type == LAGS_MINMAX (%d), min=%d, max=%d\n",
		LAGS_MINMAX, lpref->lspec.lminmax[0], lpref->lspec.lminmax[1]);
    } else if (lpref->spectype == LAGS_LIST) {
	printlist(lpref->lspec.laglist, "type == LAGS_LIST");
    } 
}
#endif

/* Modify the lag preferences for a given variable, based either on
   selections made via the lags dialog box or on the specification of
   a pre-existing model.  Return 1 if the preferences in question were
   in fact changed, otherwise 0.
*/

static int 
modify_lpref (lagpref *lpref, char spectype, int lmin, int lmax, int *laglist)
{
    int mod = 1;

    if (spectype == LAGS_TMP) {
	/* special for compiling lag preferences based on a saved
	   model specification */
	gretl_list_append_term(&lpref->lspec.laglist, lmin);
#if LDEBUG > 1
	fprintf(stderr, "modify_lpref: added lag %d\n", lmin);
#endif
	return 1;
    }

    /* if we got something presented as a list of specific lags, but 
       actually it's a consecutive list, convert to "minmax" form 
    */
    if (spectype == LAGS_LIST) {
	gretl_list_sort(laglist);
	if (gretl_list_is_consecutive(laglist)) {
	    lmin = laglist[1];
	    lmax = laglist[laglist[0]];
	    free(laglist);
	    spectype = LAGS_MINMAX;
	}
    }

    if (spectype == lpref->spectype) {
	if (spectype == LAGS_LIST) {
	    if (!gretl_list_cmp(laglist, lpref->lspec.laglist)) {
		free(laglist);
		mod = 0;
	    }
	} else if (spectype == LAGS_MINMAX) {
	    if (lmin == lpref->lspec.lminmax[0] &&
		lmax == lpref->lspec.lminmax[1]) {
		mod = 0;
	    }
	} else if (spectype == LAGS_NONE) {
	    mod = 0;
	} 
    }

    if (mod == 0) {
	/* no-op: no (effective) change made */
	return mod;
    }

    /* beyond this point we're making changes to the saved
       lag preferences */

    if (lpref->spectype == LAGS_LIST) {
	free(lpref->lspec.laglist);
	lpref->lspec.laglist = NULL;
    }    

    if (spectype == LAGS_MINMAX) {
	lpref->lspec.lminmax[0] = lmin;
	lpref->lspec.lminmax[1] = lmax;
    } else if (spectype == LAGS_LIST) {
	lpref->lspec.laglist = laglist;
    } else if (spectype == LAGS_NONE) {
	lpref->lspec.lminmax[0] = 0;
	lpref->lspec.lminmax[1] = 0;
    }

    lpref->spectype = spectype;

#if LDEBUG > 1
    fprintf(stderr, "modify_lpref: type=%d, lmin=%d, lmax=%d list=%p\n", 
	    spectype, lmin, lmax, (void *) laglist);
    print_lpref(lpref);
#endif    

    return mod;
}

static lagpref *get_saved_lpref (int v, char context)
{
    lagpref *lpref = NULL;
    int i;

    for (i=0; i<n_prefs; i++) {
	if (lprefs[i]->v == v && lprefs[i]->context == context) {
	    lpref = lprefs[i];
	    break;
	}
    }

#if LDEBUG > 1
    fprintf(stderr, "get_saved_lpref: v=%d, context=%d\n", v, context);
    print_lpref(lpref);
#endif    

    return lpref;
}

const int *get_VAR_lags_list (void)
{
    lagpref *lp = get_saved_lpref(VDEFLT, LAG_Y_V);

    if (lp != NULL && lp->spectype == LAGS_LIST) {
	return lp->lspec.laglist;
    } else {
	return NULL;
    }
}

void set_VAR_max_lag (int lmax)
{
    lagpref *lp = get_saved_lpref(VDEFLT, LAG_Y_V);

    if (lp != NULL) {
	lp->lspec.lminmax[1] = lmax;
    }
}

/* determine if a variable in a listbox is just a "dummy" lag
   entry or not */

int is_lag_dummy (int v, int lag, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int ynum = selector_get_depvar_number(NULL);
    int ret = 0;

    if (v > 0 && v == ynum) {
	ret = 1;
    } else if (lpref != NULL) {
	if (lpref->spectype == LAGS_LIST &&
	    lag > lpref->lspec.laglist[1]) {
	    ret = 1;
	} else if (lpref->spectype == LAGS_MINMAX &&
		   lag > lpref->lspec.lminmax[0]) {
	    ret = 1;
	}
    }

    return ret;
}

static void maybe_destroy_depvar_lags (lagpref *lpref, char context)
{
    if (lpref != NULL) {
	if (lpref->spectype == LAGS_MINMAX &&
	    lpref->lspec.lminmax[0] == 0 &&
	    lpref->lspec.lminmax[1] == 0) {
	    lpref->lspec.lminmax[0] = 1;
	    lpref->lspec.lminmax[1] = 1;
	    enable_lags_for_context(context, FALSE);
	}
    } else {
	enable_lags_for_context(context, FALSE);
    }	    
}

/* called when a specific lag, e.g. foo(-3), is removed from
   a list of selected variables with the mouse */

int remove_specific_lag (int v, int lag, char context)
{
    int ynum = selector_get_depvar_number(NULL);
    lagpref *lpref;
    int lmin, lmax;
    int err = 0;

    if (v == ynum) {
	if (context == LAG_X) {
	    context = LAG_Y_X;
	} else if (context == LAG_W) {
	    context = LAG_Y_W;
	}
    }

    lpref = get_saved_lpref(v, context);

    if (lpref == NULL) {
	err = 1;
    } else {
	if (lpref->spectype == LAGS_LIST) {
	    int pos = in_gretl_list(lpref->lspec.laglist, lag);

	    if (pos == 0) {
		err = 1;
	    } else {
		gretl_list_delete_at_pos(lpref->lspec.laglist, pos);
		if (lpref->lspec.laglist[0] == 0) {
		    /* no values left */
		    modify_lpref(lpref, LAGS_NONE, 0, 0, NULL);
		} else if (gretl_list_is_consecutive(lpref->lspec.laglist)) {
		    /* convert lag spec to minmax form */
		    int l0 = lpref->lspec.laglist[0];

		    lmin = lpref->lspec.laglist[1];
		    lmax = lpref->lspec.laglist[l0];
		    modify_lpref(lpref, LAGS_MINMAX, lmin, lmax, NULL);
		}
	    }
	} else if (lpref->spectype == LAGS_MINMAX) {
	    lmin = lpref->lspec.lminmax[0];
	    lmax = lpref->lspec.lminmax[1];
	    if (lag < lmin || lag > lmax) {
		err = 1;
	    } else if (lag == lmin && lmin == lmax) {
		lpref->lspec.lminmax[0] = 0;
		lpref->lspec.lminmax[1] = 0;
	    } else if (lag == lmin) {
		lpref->lspec.lminmax[0] += 1;
	    } else if (lag == lmax) {
		lpref->lspec.lminmax[1] -= 1;
	    } else {
		/* convert lag spec to list form */
		int *llist;

		llist = gretl_list_new(lmax - lmin);
		if (llist == NULL) {
		    err = 1;
		} else {
		    int i, j = 1;

		    for (i=lmin; i<=lmax; i++) {
			if (i != lag) {
			    llist[j++] = i;
			}
		    }
		    modify_lpref(lpref, LAGS_LIST, 0, 0, llist);
		}
	    }
	} else {
	    err = 1;
	}
    }

    /* special handling of dependent var lags */
    if (context == LAG_Y_X || context == LAG_Y_W) {
	maybe_destroy_depvar_lags(lpref, context);
    }

    return err;
}

static lagpref *lpref_add (int v, char context, int type)
{
    lagpref *lpref = lpref_new(v, context, type);

    if (lpref != NULL) {
	if (add_lpref_to_stack(lpref)) {
	    free(lpref);
	    lpref = NULL;
	}
    }

    return lpref;
}

int set_lag_prefs_from_list (int v, int *llist, char context,
			     int *changed)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int mod, err = 0;

    *changed = 0;

    if (lpref == NULL) {
	*changed = 1;
	lpref = lpref_add(v, context, LAGS_NONE);
    }

    if (lpref == NULL) {
	err = E_ALLOC;
    } else {
	mod = modify_lpref(lpref, LAGS_LIST, 0, 0, llist);
	if (!*changed && mod) {
	    *changed = 1;
	}
    }

    return err;
}

static int minmax_defaults (int lmin, int lmax, char context)
{
    if ((context == LAG_X || context == LAG_W) &&
	lmin == 0 && lmax == 0) {
	return 1;
    }

    if ((context == LAG_Y_X || context == LAG_Y_W) &&
	lmin == 1 && lmax == 1) {
	return 1;
    }

    return 0;
}

int set_lag_prefs_from_minmax (int v, int lmin, int lmax,
			       char context, int *changed)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int mod, err = 0;

    *changed = 0;

    if (lpref == NULL && minmax_defaults(lmin, lmax, context)) {
	return 0;
    }

    if (lpref == NULL) {
	*changed = 1;
	lpref = lpref_add(v, context, LAGS_NONE);
    }

    if (lpref == NULL) {
	err = E_ALLOC;
    } else {
	mod = modify_lpref(lpref, LAGS_MINMAX, lmin, lmax, NULL);
	if (mod) {
	    *changed = 1;
	}	
    } 

    return err;
}

/* push a specific lag onto the list of lags for a variable */

static int set_lag_pref_from_lag (int v, int lag, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int err = 0;

    if (lpref == NULL) {
	lpref = lpref_add(v, context, LAGS_LIST);
    } 

    if (lpref == NULL) {
	err = E_ALLOC;
    } else {
	modify_lpref(lpref, LAGS_TMP, lag, 0, NULL);
	if (lpref->lspec.laglist == NULL) {
	    err = E_ALLOC;
	}
    } 

    return err;
}

static int parent_in_list (const int *list, int p)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == p || list[i] == -p) {
	    return 1;
	}
    }

    return 0;
}

/* get lag preferences out of a list of regressors or instruments
   from a pre-existing model */

static int set_lag_prefs_from_xlist (int *list, int dv, char cbase,
				     int *pnset)
{
    char context;
    int i, vi, lag, pv;
    int nset = 0;
    int err = 0;

    /* Search the list for recognizable lags and if we find any,
       consolidate the lag information by variable
    */

    for (i=1; i<=list[0] && !err; i++) {
	vi = list[i];
	if (vi == 0 || vi == LISTSEP) {
	    continue;
	}
	pv = 0;
	lag = series_get_lag(dataset, vi);
	if (lag) {
	    pv = series_get_parent_id(dataset, vi);
	}	
	if (pv > 0) {
	    context = (pv == dv)? cbase + 1 : cbase;
	    err = set_lag_pref_from_lag(pv, lag, context);
	    if (!err) {
		if (pv != dv && !parent_in_list(list, pv)) {
		    /* convert to "ghost" list entry for parent variable */
		    list[i] = -pv;
		} else {
		    /* delete lagged var from list to avoid duplication */
		    gretl_list_delete_at_pos(list, i--);
		}
		enable_lags_for_context(context, TRUE);
		nset++;
	    } 
	}
    }

    /* Now check whether any lagged vars were also present in
       contemporaneous form.  Note that recognizable lags have been
       removed by now, and also that we don't need to bother about the
       dependent variable.  Also note that negative variable numbers
       mark positions in the list where a "ghost" entry has been
       inserted for a parent variable that is present only in lagged
       form: at this point we restore the correct variable number.
    */

    if (!err && nset > 0) {
	lagpref *lpref;

	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == LISTSEP) {
		continue;
	    }
	    if (vi < 0) {
		/* fix "ghosted" contemporaneous var */
		list[i] = -vi;
	    } else {
		/* genuine contemporaneous var */
		lpref = get_saved_lpref(vi, cbase);
		if (lpref != NULL) {
		    /* insert lag 0 */
		    set_lag_pref_from_lag(vi, 0, cbase);
		}
	    }
	}
    }

    if (pnset != NULL) {
	*pnset = nset;
    }

    return err;
}

/* convert any singleton or consecutive lists of lags (produced
   piecewise from a model's lists) to min/max form 
*/

static void clean_up_model_prefs (void)
{
    int *list;
    int i, lmin, lmax;

    for (i=0; i<n_prefs; i++) {
	list = gretl_list_sort(lprefs[i]->lspec.laglist);
	lmin = list[1];
	lmax = list[list[0]];
	if (lmax - lmin == list[0] - 1) {
	    /* singleton or consecutive */
	    modify_lpref(lprefs[i], LAGS_MINMAX, lmin, lmax, NULL);
	} 
    }

}

/* Read the lists of regressors and instruments from a saved model,
   parse the lags info in these lists and convert to saved
   "lagpref" form.
*/

int set_lag_prefs_from_model (int dv, int *xlist, int *zlist)
{
    int n, nset = 0;
    int err = 0;

    /* start with a clean slate */
    destroy_lag_preferences();
    enable_lags_for_context(LAG_Y_X, FALSE);
    enable_lags_for_context(LAG_Y_W, FALSE);

    if (xlist != NULL) {
	/* regressors */
	err = set_lag_prefs_from_xlist(xlist, dv, LAG_X, &n);
	nset += n;
    }

    if (!err && zlist != NULL) {
	/* instruments */
	err = set_lag_prefs_from_xlist(zlist, dv, LAG_W, &n);
	nset += n;
    }

    if (nset > 0) {
	clean_up_model_prefs();
    }

    return err;
}

int set_lag_prefs_from_VAR (const int *lags, int *xlist)
{
    int err = 0;    

    /* start with a clean slate */
    destroy_lag_preferences();

    /* lags of endogenous variables */
    if (lags != NULL) {
	int *list = gretl_list_copy(lags);
	int ch;

	if (list == NULL) {
	    err = E_ALLOC;
	} else {
	    /* set directly from a given lag list */
	    err = set_lag_prefs_from_list(VDEFLT, list, LAG_Y_V, &ch);
	}
    }

    /* lags among exogenous vars */
    if (!err && xlist != NULL) {
	int nset = 0;

	err = set_lag_prefs_from_xlist(xlist, -1, LAG_X, &nset);
	if (nset > 0) {
	    clean_up_model_prefs();
	}
    }

    return err;
}

void set_null_lagpref (int v, char context, int *changed)
{
    lagpref *lpref = get_saved_lpref(v, context);

    if (lpref != NULL) {
	*changed = modify_lpref(lpref, LAGS_NONE, 0, 0, NULL);
    } else {
	*changed = 0;
    }
}

void
get_lag_preference (int v, int *lmin, int *lmax, const int **laglist,
		    char context, selector *sr)
{
    lagpref *lpref = get_saved_lpref(v, context);

    *lmin = *lmax = 0;
    *laglist = NULL;

    if (context == LAG_Y_V && lpref == NULL) {
	*lmin = 1;
	if (sr != NULL) {
	    *lmax = selector_get_VAR_order(sr);
	}
	return;
    }

    if (!lags_enabled_for_context(context)) {
	return;
    }

    if ((context == LAG_Y_X || context == LAG_Y_W) && lpref == NULL) {
	*lmin = *lmax = 1;
	return;
    }  
    
    if (lpref == NULL || v >= dataset->v) {
	return;
    }

    if (lpref->spectype == LAGS_LIST) {
	*laglist = lpref->lspec.laglist;
    } else if (lpref->spectype == LAGS_MINMAX) {
	*lmin = lpref->lspec.lminmax[0];
	*lmax = lpref->lspec.lminmax[1];
    } 
}

int *get_lag_pref_as_list (int v, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int *list = NULL;

#if LDEBUG
    fprintf(stderr, "get_lag_pref_as_list: var = %d, context = %d\n", v, context);
#endif

    if (!lags_enabled_for_context(context)) {
	return NULL;
    }

    if ((context == LAG_Y_X || context == LAG_Y_W) && lpref == NULL) {
	/* the default: a single lag */
	list = gretl_list_new(1);
	if (list != NULL) {
	    list[1] = 1;
	}
	return list;
    }

    if (lpref != NULL) { 
	if (lpref->spectype == LAGS_LIST) {
	    list = gretl_list_copy(lpref->lspec.laglist);
	} else if (lpref->spectype == LAGS_MINMAX) {
	    list = gretl_consecutive_list_new(lpref->lspec.lminmax[0],
					      lpref->lspec.lminmax[1]);
	}
    }

#if LDEBUG
    printlist(list, "returned list");
#endif

    return list;
}
