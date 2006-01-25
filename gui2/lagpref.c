
/* Mechanism for recording and retrieving the preferred set of lags
   for a given variable.  We save the preferences that are set, and
   destroy the stacked preferences when we switch data sets.
*/

enum {
    LAGS_NONE,
    LAGS_MINMAX,
    LAGS_LIST
};

typedef struct lagpref_ lagpref;

struct lagpref_ {
    int v;
    char ltype;
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

static lagpref *lpref_new (int v)
{
    lagpref *lpref = malloc(sizeof *lpref);
    
    if (lpref == NULL) {
	return NULL;
    }

    lpref->v = v;
    lpref->ltype = LAGS_NONE;

    return lpref;
}

static int 
modify_lpref (lagpref *lpref, char ltype, int lmin, int lmax, int *laglist)
{
    if (lpref->ltype == LAGS_LIST) {
	free(lpref->lspec.laglist);
    }    

    if (ltype == LAGS_MINMAX) {
	lpref->lspec.lminmax[0] = lmin;
	lpref->lspec.lminmax[1] = lmax;
    } else if (ltype == LAGS_LIST) {
	lpref->lspec.laglist = laglist;
    } else if (ltype == LAGS_NONE) {
	lpref->lspec.lminmax[0] = 0;
	lpref->lspec.lminmax[1] = 0;
    }

    lpref->ltype = ltype;

    return 0;
}

static lagpref *get_lpref_by_varnum (int v)
{
    int i;

    for (i=0; i<n_prefs; i++) {
	if (lprefs[i]->v == v) {
	    return lprefs[i];
	}
    }

    return NULL;
}

static int set_lag_prefs_from_list (int v, int *llist)
{
    lagpref *lpref = get_lpref_by_varnum(v);
    int err = 0;

    if (lpref == NULL) {
	lpref = lpref_new(v);
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

static int set_lag_prefs_from_minmax (int v, int lmin, int lmax)
{
    lagpref *lpref = get_lpref_by_varnum(v);
    int err = 0;

    if (lpref == NULL) {
	lpref = lpref_new(v);
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

static void set_null_lagpref (int v)
{
    lagpref *lpref = get_lpref_by_varnum(v);

    if (lpref != NULL) {
	modify_lpref(lpref, LAGS_NONE, 0, 0, NULL);
    }
}

static void
get_lag_preference (int v, int *lmin, int *lmax, int **laglist)
{
    lagpref *lpref = get_lpref_by_varnum(v);

    *lmin = *lmax = 0;
    *laglist = NULL;
    
    if (lpref == NULL || v >= datainfo->v) {
	return;
    }

    if (lpref->ltype == LAGS_LIST) {
	*laglist = lpref->lspec.laglist;
    } else if (lpref->ltype == LAGS_MINMAX) {
	*lmin = lpref->lspec.lminmax[0];
	*lmax = lpref->lspec.lminmax[1];
    } 
}

static void destroy_lag_preferences (void)
{
    int i;
    
    for (i=0; i<n_prefs; i++) {
	if (lprefs[i]->ltype == LAGS_LIST) {
	    free(lprefs[i]->lspec.laglist);
	}
	free(lprefs[i]);
    } 

    free(lprefs);
    lprefs = NULL;
    n_prefs = 0;
}
