
/* Mechanism for recording and retrieving the preferred set of lags
   for a given variable.  We save the preferences that are set, and
   destroy the stacked preferences when we switch data sets.
*/

enum {
    LAGS_NONE,    /* no lags specified */
    LAGS_MINMAX,  /* min and max lags given (consecutive) */
    LAGS_LIST     /* list of lags specific given */
} SpecType;

enum {
    LAG_X,        /* lags set for regular variable context */
    LAG_INSTR     /* lags set for the variable as instrument */
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

static lagpref *lpref_new (int v, char context)
{
    lagpref *lpref = malloc(sizeof *lpref);
    
    if (lpref == NULL) {
	return NULL;
    }

    lpref->v = v;
    lpref->spectype = LAGS_NONE;
    lpref->context = context;

    return lpref;
}

static int 
modify_lpref (lagpref *lpref, char spectype, int lmin, int lmax, int *laglist)
{
    if (lpref->spectype == LAGS_LIST) {
	free(lpref->lspec.laglist);
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

    return 0;
}

static lagpref *get_saved_lpref (int v, char context)
{
    int i;

    for (i=0; i<n_prefs; i++) {
	if (lprefs[i]->v == v && lprefs[i]->context == context) {
	    return lprefs[i];
	}
    }

    return NULL;
}

static int set_lag_prefs_from_list (int v, int *llist, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int err = 0;

    if (lpref == NULL) {
	lpref = lpref_new(v, context);
	if (lpref == NULL) {
	    err = E_ALLOC;
	} else {
	    err = add_lpref_to_stack(lpref);
	}
    }

    if (!err) {
	modify_lpref(lpref, LAGS_LIST, 0, 0, llist);
    }

    return err;
}

static int set_lag_prefs_from_minmax (int v, int lmin, int lmax,
				      char context)
{
    lagpref *lpref = get_saved_lpref(v, context);
    int err = 0;

    if (lpref == NULL) {
	lpref = lpref_new(v, context);
	if (lpref == NULL) {
	    err = E_ALLOC;
	} else {
	    err = add_lpref_to_stack(lpref);
	}
    }

    if (!err) {
	modify_lpref(lpref, LAGS_MINMAX, lmin, lmax, NULL);
    }    

    return err;
}

static void set_null_lagpref (int v, char context)
{
    lagpref *lpref = get_saved_lpref(v, context);

    if (lpref != NULL) {
	modify_lpref(lpref, LAGS_NONE, 0, 0, NULL);
    }
}

static void
get_lag_preference (int v, int *lmin, int *lmax, const int **laglist,
		    char context)
{
    lagpref *lpref = get_saved_lpref(v, context);

    *lmin = *lmax = 0;
    *laglist = NULL;
    
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
