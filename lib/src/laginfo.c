/* for inclusion in interact.c: provides a mechanism for handling
   the fallout from the auto-generation of lags when an estimation
   command using the syntax of "foo(-1)" in the regression list,
   or "foo(-1 to -4)", etc.
*/

struct Laginfo_ {
    int *reflist;
    int *genlist;
    int **lag_lists;
};

static Laginfo *list_lag_info_new (void)
{
    Laginfo *linfo = malloc(sizeof *linfo);

    if (linfo != NULL) {
	linfo->reflist = NULL;
	linfo->genlist = NULL;
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
	free(linfo->genlist);
	free(linfo);
    }
}

static void cmd_lag_info_destroy (CMD *cmd)
{
    list_lag_info_destroy(cmd->linfo);
    cmd->linfo = NULL;
}

static const int *
get_lag_list_by_varnum (int v, const Laginfo *linfo)
{
    int *list = NULL;

    if (linfo != NULL && linfo->reflist != NULL) {
	int i;
	
	for (i=1; i<=linfo->reflist[0]; i++) {
	    if (linfo->reflist[i] == v) {
		list = linfo->lag_lists[i-1];
		break;
	    }
	}
    }

    return list;
}

static int add_lagv_to_genlist (int lv, Laginfo *linfo)
{
    int *genlist;
    int n, err = 0;

    if (linfo->genlist != NULL) {
	n = linfo->genlist[0] + 2;
    } else {
	n = 2;
    }

    genlist = realloc(linfo->genlist, n * sizeof *genlist);
    if (genlist != NULL) {
	linfo->genlist = genlist;
	genlist[n - 1] = lv;
	genlist[0] = n - 1;
    } else {
	err = 1;
    }

    return err;
}

static int add_lag_to_laglist (int i, int lag, Laginfo *linfo)
{
    int err = 0;

    if (linfo->lag_lists[i] == NULL) {
	linfo->lag_lists[i] = gretl_list_new(1);
	if (linfo->lag_lists[i] != NULL) {
	    linfo->lag_lists[i][1] = lag;
	} else {
	    err = 1;
	}
    } else {
	int n = linfo->lag_lists[i][0] + 2;
	int *laglist;

	laglist = realloc(linfo->lag_lists[i], n * sizeof *laglist);
	if (laglist != NULL) {
	    linfo->lag_lists[i] = laglist;
	    laglist[n - 1] = lag;
	    laglist[0] = n - 1;
	} else {
	    err = 1;
	}
    }

    return err;
}

static int add_to_list_lag_info (int v, int lag, int lagv, CMD *cmd)
{
    int *reflist = NULL;
    int **llists = NULL;
    int add_to_reflist = 0;
    int nl, llnum = 0, err = 0;

    fprintf(stderr, "add_to_list_lag_info: v=%d, lag=%d, lagv=%d\n",
	    v, lag, lagv);

    if (cmd->linfo == NULL) {
	cmd->linfo = list_lag_info_new();
	add_to_reflist = 1;
    }

    if (cmd->linfo == NULL) {
	return 1;
    }

    if (!add_to_reflist) {
	if (get_lag_list_by_varnum(v, cmd->linfo) == NULL) {
	    /* there's no list already started for this var */
	    add_to_reflist = 1;
	}
    }

    if (cmd->linfo->reflist == NULL) {
	nl = 1;
    } else {
	nl = cmd->linfo->reflist[0] + 1;
    }

    if (add_to_reflist) {
	fprintf(stderr, "reallocing %p, %d elements\n", 
		(void *) cmd->linfo->reflist, nl + 1);
	reflist = realloc(cmd->linfo->reflist, (nl + 1) * sizeof *reflist);
	if (reflist != NULL) {
	    cmd->linfo->reflist = reflist;
	    fprintf(stderr, "set cmd->linfo->reflist = %p\n",
		    (void *) cmd->linfo->reflist);
	    cmd->linfo->reflist[0] = 0;
	} else {
	    err = 1;
	}
 
	if (!err) {
	    llists = realloc(cmd->linfo->lag_lists, nl * sizeof *llists);
	    if (llists != NULL) {
		fprintf(stderr, " realloced lag_lists, size %d\n", nl);
		cmd->linfo->lag_lists = llists;
		cmd->linfo->lag_lists[nl - 1] = NULL;
	    } else {
		err = 1;
	    }
	}
    }

    if (!err) {
	err = add_lagv_to_genlist(lagv, cmd->linfo);
	fprintf(stderr, " add_lagv_to_genlist: err = %d\n", err);
    }

    if (!err) {
	if (cmd->linfo->reflist[0] > 0) {
	    llnum = cmd->linfo->reflist[0] - 1;
	}
	fprintf(stderr, " llnum = %d\n", llnum);
    }    

    if (!err) {
	fprintf(stderr, " doing add_lag_to_laglist(%d, %d, %p)\n", 
		llnum, lag, (void *) cmd->linfo);
	err = add_lag_to_laglist(llnum, lag, cmd->linfo);
    }

    if (!err) {
	if (add_to_reflist) {
	    cmd->linfo->reflist[0] += 1;
	    cmd->linfo->reflist[nl] = v;
	}
    }

    return err;
}

static int var_lags_contiguous (const int *laglist)
{
    int i, ret = 1;

    for (i=2; i<=laglist[0]; i++) {
	if (laglist[i] != laglist[i-1] + 1) {
	    ret = 0;
	}
    }

    return ret;    
}

/* returns number of bytes printed */

static int print_var_lags (const int *laglist, PRN *prn)
{
    char tmp[32];
    int lmax = laglist[0];
    int ret = 0;
    
    /* FIXME leads */

    if (lmax == 1) {
	sprintf(tmp, "(-%d)", laglist[1]);
	ret += pputs(prn, tmp);	
    } else if (var_lags_contiguous(laglist)) {
	sprintf(tmp, "(-%d to -%d)", laglist[1], laglist[lmax]);
	ret += pputs(prn, tmp);	
    } else {
	int i;

	pputc(prn, '(');
	ret++;
	for (i=1; i<=lmax; i++) {
	    sprintf(tmp, "%d", laglist[i]);
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
print_lags_by_varnum (int v, const Laginfo *linfo, int cli, 
		      const DATAINFO *pdinfo, PRN *prn)
{
    PRN *myprn = NULL;
    const int *laglist = NULL;
    int ret = 0;

    if (cli) {
	myprn = gretl_print_new(GRETL_PRINT_STDOUT);
    } else {
	myprn = prn;
    }

    laglist = get_lag_list_by_varnum(v, linfo);
    if (laglist != NULL) {
	pputc(myprn, ' ');
	ret = 1 + pputs(myprn, pdinfo->varname[v]);
	ret += print_var_lags(laglist, myprn);
    } 

    if (cli) {
	gretl_print_destroy(myprn);
    }

    return ret;
}

static int is_auto_generated_lag (int v, const Laginfo *linfo)
{
    int i, ret = 0;

    if (linfo != NULL && linfo->genlist != NULL) {
	for (i=1; i<=linfo->genlist[0]; i++) {
	    if (v == linfo->genlist[i]) {
		ret = i;
		break;
	    }
	}
    }

    return ret;
}

static int is_first_lag (int v, const Laginfo *linfo, int *src)
{
    int i, j, k = 0, ret = 0;

    for (i=0; i<linfo->reflist[0]; i++) {
	for (j=0; j<linfo->lag_lists[i][0]; j++) {
	    k++;
	    if (k == v) {
		if (j == 0) {
		    ret = 1;
		    *src = linfo->reflist[i+1];
		}
		break;
	    }
	}
    }

    return ret;
}



	    

    
