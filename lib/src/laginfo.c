/* for inclusion in interact.c: provides a mechanism for handling
   the fallout from the auto-generation of lags when an estimation
   command using the syntax of "foo(-1)" in the regression list,
   or "foo(-1 to -4)", etc.
*/

struct Laginfo_ {
    int *reflist;    /* list of distinct var for which we'll generate lags */
    int **lag_lists; /* list of lags to be generated, per var */
    int *srclist;    /* "shadow" of command list, with "source" IDs in place
			of lag vars IDs */
};

#define LLDEBUG 0

static Laginfo *list_lag_info_new (void)
{
    Laginfo *linfo = malloc(sizeof *linfo);

    if (linfo != NULL) {
	linfo->reflist = NULL;
	linfo->srclist = NULL;
	linfo->lag_lists = NULL;
    }

    return linfo;
}

static void list_lag_info_destroy (Laginfo *linfo)
{
    if (linfo != NULL) {
	if (linfo->reflist != NULL) {
	    int i;

	    for (i=0; i<linfo->reflist[0]; i++) {
		free(linfo->lag_lists[i]);
	    }
	    free(linfo->lag_lists);
	    free(linfo->reflist);
	}
	free(linfo->srclist);
	free(linfo);
    }
}

static void cmd_lag_info_destroy (CMD *cmd)
{
    list_lag_info_destroy(cmd->linfo);
    cmd->linfo = NULL;
}

/* unlike gretl_list_separator_position(), we start the check
   at position 1 here */

static int laglist_sep_pos (const int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    return i;
	}
    }

    return 0;
}

/* retrieve list of lags for variable v */

static const int *get_lag_list_by_varnum (int v, const Laginfo *linfo)
{
    int *list = NULL;
    int i;

    if (linfo != NULL && linfo->reflist != NULL) {
	for (i=1; i<=linfo->reflist[0]; i++) {
	    if (linfo->reflist[i] == v) {
		list = linfo->lag_lists[i-1];
		break;
	    }
	}
    }

    return list;
}

/* add a specific lag to the record of lags for a given
   variable; insert a list separator first, if needed 
*/

static int 
add_lag_to_laglist (int srcv, int lag, int spos, Laginfo *linfo)
{
    int lspos, j, i = -1;
    int err = 0;

    for (j=1; j<=linfo->reflist[0]; j++) {
	if (linfo->reflist[j] == srcv) {
	    i = j - 1;
	    break;
	}
    }

    if (i < 0) {
	/* "can't happen" */
	return E_DATA;
    }

#if LLDEBUG
    fprintf(stderr, "add_lag_to_laglist: srcv=%d, i=%d, lag=%d, spos=%d, starting\n",
	    srcv, i, lag, spos);
    printlist(linfo->lag_lists[i], NULL);
#endif

    lspos = laglist_sep_pos(linfo->lag_lists[i]);

    if (spos && !lspos) {
	linfo->lag_lists[i] = gretl_list_append_term(&linfo->lag_lists[i], LISTSEP);
	if (linfo->lag_lists[i] == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	linfo->lag_lists[i] = gretl_list_append_term(&linfo->lag_lists[i], lag);
	if (linfo->lag_lists[i] == NULL) {
	    err = E_ALLOC;
	}
    }	

#if LLDEBUG
    printlist(linfo->lag_lists[i], "after adding");
#endif

    return err;
}

/* add a new lag list for a variable that doesn't already have one */

static int laginfo_add_lags_list (int n, Laginfo *linfo)
{
    int **llists;
    int err = 0;

    llists = realloc(linfo->lag_lists, n * sizeof *llists);

    if (llists != NULL) {
#if LLDEBUG
	fprintf(stderr, " realloced lag_lists, size %d\n", n);
#endif
	linfo->lag_lists = llists;
	linfo->lag_lists[n - 1] = gretl_null_list();
	if (linfo->lag_lists[n - 1] == NULL) {
	    err = E_ALLOC;
	}
    } else {
	err = E_ALLOC;
    }

    return err;
}

/* expand the master list that keeps track of which variables 
   have lag lists */

static int laginfo_expand_reflist (int n, int v, Laginfo *linfo)
{
    int *reflist;
    int l0, err = 0;

    if (linfo->reflist == NULL) {
	l0 = 1;
    } else {
	l0 = linfo->reflist[0] + 1;
    }

#if LLDEBUG
    fprintf(stderr, " now reallocing %p, %d elements\n", 
	    (void *) linfo->reflist, n);
#endif

    reflist = realloc(linfo->reflist, n * sizeof *reflist);

    if (reflist != NULL) {
	linfo->reflist = reflist;
	linfo->reflist[0] = l0;
	linfo->reflist[l0] = v;
    } else {
	err = E_ALLOC;
    }

    return err;
}

/* expand linfo->srclist, which shadows the main command list, keeping
   a record of the "source" variable for any lagged terms
*/

static int expand_srclist (int srcv, int lpos, Laginfo *linfo,
			   const int *cmdlist)
{
    int i, n = lpos + 1;
    int err = 0;

#if LLDEBUG
    printlist(linfo->srclist, "srclist, before expansion");
#endif

    if (linfo->srclist == NULL) {
	linfo->srclist = gretl_list_new(lpos);
	if (linfo->srclist == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=1; i<lpos; i++) {
		linfo->srclist[i] = cmdlist[i];
	    }
	}
    } else {
	int m = linfo->srclist[0];
	int *tmp = realloc(linfo->srclist, n * sizeof *tmp);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    linfo->srclist = tmp;
	    linfo->srclist[0] = lpos;
	    for (i=m+1; i<lpos; i++) {
		linfo->srclist[i] = cmdlist[i];
	    }
	}
    }

    if (!err) {
	linfo->srclist[lpos] = srcv;
    }

#if LLDEBUG
    printlist(linfo->srclist, "srclist, after");
#endif

    return err;
}

/* Add lag info: 'srcv' is the "source" variable; 'lag'
   is the lag order; 'lagv' is the ID of the created variable;
   and 'lpos' is the current position in the command list.
*/

static int list_lag_info_add (int srcv, int lag, int lagv, int lpos, CMD *cmd)
{
    int spos = gretl_list_separator_position(cmd->list);
    int add_to_reflist = 0;
    int nl, err = 0;

#if LLDEBUG
    fprintf(stderr, "*** list_lag_info_add: srcv=%d, lag=%d -> lagv=%d, lpos = %d\n",
	    srcv, lag, lagv, lpos);
    printlist(cmd->list, "cmd->list");
#endif

    if (cmd->linfo == NULL) {
	cmd->linfo = list_lag_info_new();
	if (cmd->linfo == NULL) {
	    return E_ALLOC;
	}
	add_to_reflist = 1;
    }

    if (!add_to_reflist) {
	if (get_lag_list_by_varnum(srcv, cmd->linfo) == NULL) {
	    add_to_reflist = 1;
	}
    }

    if (cmd->linfo->reflist == NULL) {
	nl = 1;
    } else {
	nl = cmd->linfo->reflist[0] + 1;
    }

    if (add_to_reflist) {
	err = laginfo_add_lags_list(nl, cmd->linfo);
	if (!err) {
	    err = laginfo_expand_reflist(nl + 1, srcv, cmd->linfo);
	}
    }

    if (!err) {
	err = expand_srclist(srcv, lpos, cmd->linfo, cmd->list);
    }

    if (!err) {
#if LLDEBUG
	printlist(cmd->linfo->reflist, "reflist");
#endif
	err = add_lag_to_laglist(srcv, lag, spos, cmd->linfo);
    }

    return err;
}

static int 
var_lags_contiguous (const int *laglist, int lstart, int lmax)
{
    int i, ret = 1;

    /* actual lags */
    for (i=lstart+1; i<=lmax; i++) {
	if (laglist[i] != laglist[i-1] + 1) {
	    ret = 0;
	    break;
	}
    }

    if (ret == 0) {
	/* check for contiguous leads? */
	int test = 1;

	for (i=lstart+1; i<=lmax; i++) {
	    if (laglist[i] != laglist[i-1] - 1) {
		test = 0;
		break;
	    }
	}
	if (test == 1) {
	    ret = 1;
	}
    }	

    return ret;    
}

static const char *lag_sign_str (int lag)
{
    if (lag > 0) return "-";
    if (lag < 0) return "+";
    else return "";
}

static void 
get_lstart_lmax (const int *llist, int gotsep,
		 int *lstart, int *lmax)
{
    int spos = laglist_sep_pos(llist);

    if (spos == 0) {
	*lstart = 1;
	*lmax = llist[0];
    } else if (!gotsep) {
	if (spos == 1) {
	    *lstart = 2;
	    *lmax = 0;
	} else {
	    *lstart = 1;
	    *lmax = spos - 1;
	}
    } else {
	*lstart = spos + 1;
	*lmax = llist[0];
    }
}

/* returns number of bytes printed */

static int print_var_lags (const int *laglist, int gotsep,
			   PRN *prn)
{
    char tmp[32];
    int lag, lsign;
    int lstart, lmax;
    int i, n, ret = 0;

    get_lstart_lmax(laglist, gotsep, &lstart, &lmax);
    n = lmax - lstart + 1;

    if (n < 1) {
	/* no actual lags, this sublist */
	return 0;
    } else if (n == 1) {
	/* just one lag */
	lsign = laglist[lstart];
	lag = abs(laglist[lstart]);
	sprintf(tmp, "(%s%d)", lag_sign_str(lsign), lag);
	ret += pputs(prn, tmp);	
    } else if (var_lags_contiguous(laglist, lstart, lmax)) {
	/* first lag */
	lsign = laglist[lstart];
	lag = abs(laglist[lstart]);
	sprintf(tmp, "(%s%d to ", lag_sign_str(lsign), lag);
	ret += pputs(prn, tmp);	
	/* last lag */ 
	lsign = laglist[lmax];
	lag = abs(laglist[lmax]);
	sprintf(tmp, "%s%d)", lag_sign_str(lsign), lag);
	ret += pputs(prn, tmp);	
    } else {
	pputc(prn, '(');
	ret++;
	for (i=lstart; i<=lmax; i++) {
	    lsign = laglist[i];
	    lag = fabs(laglist[i]);
	    sprintf(tmp, "%s%d", lag_sign_str(lsign), lag);
	    ret += pputs(prn, tmp);
	    if (i < lmax) {
		ret += pputs(prn, ", ");
	    } else {
		pputc(prn, ')');
		ret++;
	    }
	}
    }

    return ret;
}

static int 
print_lags_by_varnum (int v, const Laginfo *linfo, 
		      const DATASET *dset, 
		      int gotsep, PRN *prn)
{
    const int *laglist = NULL;
    int ret = 0;

#if LLDEBUG
    fprintf(stderr, "print_lags_by_varnum: v = %d, gotsep = %d\n",
	    v, gotsep);
#endif

    laglist = get_lag_list_by_varnum(v, linfo);
    if (laglist != NULL) {
	pputc(prn, ' ');
	ret = 1 + pputs(prn, dset->varname[v]);
	ret += print_var_lags(laglist, gotsep, prn);
    } 

    return ret;
}

/* below: FIXME case of i == 1 ?? (auto-lagged dependent var) */

static int 
is_auto_generated_lag (int i, const int *cmdlist, const Laginfo *linfo)
{
    int v = cmdlist[i];

    if (linfo == NULL || linfo->srclist == NULL) {
	return 0;
    }

    if (i > linfo->srclist[0]) {
	return 0;
    }

    if (!in_gretl_list(linfo->srclist, v)) {
	return 1;
    }    

    if (in_gretl_list(linfo->reflist, v)) {
	/* v could be the unlagged "parent" for lagged followers:
	   to get this right we need to take sublists into account
	*/
	int spos = gretl_list_separator_position(cmdlist);

	if (spos > 0) {
	    int lmin = (i < spos)? 2 : spos + 1;
	    int lmax = (i < spos)? spos - 1 : linfo->srclist[0];
	    int j;

	    if (lmax > linfo->srclist[0]) {
		lmax = linfo->srclist[0];
	    }

	    for (j=lmin; j<=lmax; j++) {
		if (linfo->srclist[j] == v) {
		    return 1;
		}
	    }
	} else {
	    return 1;
	}
    }

    return 0;
}

/* Determine whether the var that appears at position i in
   cmd->list is the first lag of any variable that appears in
   the command list: if so, it will be used as an 'anchor' for echoing
   other lags of the same variable.  Complication: 'first' really
   means 'first within a given sublist' if the command list contains
   sublists.  (Conversely, if a lag var is *not* 'first' in this
   sense it should not be echoed separately, since it will have
   benn handled already.)
*/

static int 
is_first_lag (int i, const int *cmdlist, int sep, const Laginfo *linfo, int *src)
{
    int spos = gretl_list_separator_position(cmdlist);
    int srcv = linfo->srclist[i];
    int lmin, lmax, j;

    if (spos) {
	lmin = (sep)? spos + 1 : 2;
	lmax = (sep)? linfo->srclist[0] : spos - 1;
    } else {
	lmin = 2;
	lmax = linfo->srclist[0];
    }

#if LLDEBUG
    fprintf(stderr, "is_first_lag: i=%d, sep=%d, spos=%d, lmin=%d, lmax=%d, srcv=%d\n",
	    i, sep, spos, lmin, lmax, srcv);
#endif

    for (j=lmin; j<=lmax; j++) {
	if (linfo->srclist[j] == srcv) {
	    if (j < i) {
		return 0;
	    } else if (j == i) {
		if (src != NULL) {
		    *src = srcv;
		}
		return 1;
	    }
	}
    }

    return 0;
}
