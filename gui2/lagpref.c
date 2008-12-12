
/* Mechanism for recording and retrieving the preferred set of lags
   for a given variable.  We save the preferences that are set, and
   destroy the stacked preferences when we switch data sets.
*/

#define LDEBUG 0

enum {
    LAGS_NONE,    /* no lags specified */
    LAGS_MINMAX,  /* min and max lags given (consecutive) */
    LAGS_LIST,    /* list of lags specific given */
    LAGS_TMP      /* provisional list when working from model */
} SpecType;

enum {
    LAG_X = 1,    /* lags set for regular variable context */
    LAG_Y_X,      /* lags for dependent variable */
    LAG_W,        /* lags set for variable as instrument */
    LAG_Y_W       /* lags for dependent var as instrument */
} LagContext;

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

static void destroy_lag_preferences (void)
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

/* return 1 if actually changed, else 0 */

static int 
modify_lpref (lagpref *lpref, char spectype, int lmin, int lmax, int *laglist)
{
    int mod = 1;

    if (spectype == LAGS_TMP) {
	gretl_list_append_term(&lpref->lspec.laglist, lmin);
	return 1;
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
#if LDEBUG
	fprintf(stderr, "modify_lpref: no changes made\n");
#endif
	return mod;
    }

#if LDEBUG == 1
    fprintf(stderr, "modify_lpref: lmin=%d, lmax=%d, laglist=%p\n",
	    lmin, lmax, (void *) laglist);
#endif

    if (lpref->spectype == LAGS_LIST) {
	free(lpref->lspec.laglist);
    }    

    if (spectype == LAGS_MINMAX) {
	lpref->lspec.lminmax[0] = lmin;
	lpref->lspec.lminmax[1] = lmax;
    } else if (spectype == LAGS_LIST) {
	lpref->lspec.laglist = gretl_list_sort(laglist);
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

/* determine if a variable in a listbox is just a "dummy" lag
   entry or not */

static int is_lag_dummy (int v, int lag, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int ynum = selector_get_depvar_number(open_selector);
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
	    if (context == LAG_Y_X) {
		y_x_lags_enabled = 0;
	    } else {
		y_w_lags_enabled = 0;
	    }
	}
    } else {
	if (context == LAG_Y_X) {
	    y_x_lags_enabled = 0;
	} else {
	    y_w_lags_enabled = 0;
	}
    }	    
}

/* called when a specific lag, e.g. foo(-3), is removed from
   a list of selected variables with the mouse */

static int remove_specific_lag (int v, int lag, char context)
{
    int ynum = selector_get_depvar_number(open_selector);
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
	    int pos = gretl_list_position(lag, lpref->lspec.laglist);

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

static int set_lag_prefs_from_list (int v, int *llist, char context,
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

static int set_lag_prefs_from_minmax (int v, int lmin, int lmax,
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

/* Read the lists of regressors and instruments from a saved model,
   parse the lags info in these lists and convert to saved
   "lagpref" form.
*/

static int set_lag_prefs_from_model (int dv, int *xlist, int *zlist)
{
    lagpref *lpref;
    int *list;
    char cbase, context;
    int i, j, vi, lag, pv;
    int err = 0, nset = 0;

    /* start with a clean slate */
    destroy_lag_preferences();

    for (j=0; j<2 && !err; j++) {
	list = (j == 0)? xlist : zlist;
	if (list != NULL) {
	    cbase = (j == 0)? LAG_X : LAG_W;
	    for (i=1; i<=list[0] && !err; i++) {
		vi = list[i];
		lag = is_standard_lag(vi, datainfo, &pv);
		if (lag != 0) {
		    context = (pv == dv)? cbase + 1 : cbase;
		    err = set_lag_pref_from_lag(pv, lag, context);
		    if (!err) {
			/* delete lagvar from list to avoid duplication */
			gretl_list_delete_at_pos(list, i--);
			if (context == LAG_Y_X) {
			    y_x_lags_enabled = 1;
			} else if (context == LAG_Y_W) {
			    y_w_lags_enabled = 1;
			}
			nset++;
		    } 
		}
	    }
	}
    }

    /* Now check whether any lagged vars were also present in
       contemporaneous form.  Note that recognizable lags have been
       removed by now, and also that we don't need to bother about the
       dependent variable.
    */

    if (nset > 0) {
	for (j=0; j<2; j++) {
	    list = (j == 0)? xlist : zlist;
	    if (list != NULL) {
		context = (j == 0)? LAG_X : LAG_W;
		for (i=1; i<=list[0]; i++) {
		    vi = list[i];
		    lpref = get_saved_lpref(vi, context);
		    if (lpref != NULL) {
			set_lag_pref_from_lag(vi, 0, context);
		    }
		}
	    }
	}
    }

    /* Finally, convert any singleton or consecutive lists of lags
       to min/max form 
    */

    if (nset > 0) {
	int lmin, lmax;

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

    return nset;
}

static void set_null_lagpref (int v, char context, int *changed)
{
    lagpref *lpref = get_saved_lpref(v, context);

    if (lpref != NULL) {
	*changed = modify_lpref(lpref, LAGS_NONE, 0, 0, NULL);
    } else {
	*changed = 0;
    }
}

static void
get_lag_preference (int v, int *lmin, int *lmax, const int **laglist,
		    char context)
{
    lagpref *lpref = get_saved_lpref(v, context);

    *lmin = *lmax = 0;
    *laglist = NULL;

    if ((context == LAG_Y_X && !y_x_lags_enabled) ||
	(context == LAG_Y_W && !y_w_lags_enabled)) {
	return;
    }  

    if ((context == LAG_Y_X || context == LAG_Y_W) && lpref == NULL) {
	*lmin = *lmax = 1;
	return;
    }  
    
    if (lpref == NULL || v >= datainfo->v) {
	return;
    }

    if (lpref->spectype == LAGS_LIST) {
	*laglist = lpref->lspec.laglist;
    } else if (lpref->spectype == LAGS_MINMAX) {
	*lmin = lpref->lspec.lminmax[0];
	*lmax = lpref->lspec.lminmax[1];
    } 
}

static int *get_lag_pref_as_list (int v, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int *list = NULL;

    if ((context == LAG_Y_X && !y_x_lags_enabled) ||
	(context == LAG_Y_W && !y_w_lags_enabled)) {
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

    return list;
}
