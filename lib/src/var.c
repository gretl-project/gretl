/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/* var.c - vector autoregressions */  

#include "libgretl.h"   
#include "internal.h"

struct var_resids {
    int *levels_list;
    double **uhat;
    int m;
    int t1, t2;
};

enum {
    VAR_PRINT_MODELS = 1 << 0,
    VAR_DO_FTESTS    = 1 << 1,
    VAR_PRINT_PAUSE  = 1 << 2
} var_flags;

/* ......................................................  */

static int gettrend (double ***pZ, DATAINFO *pdinfo)
{
    int index;
    int t, n = pdinfo->n, v = pdinfo->v;

    index = varindex(pdinfo, "time");
    if (index < v) return index;
    
    if (dataset_add_vars(1, pZ, pdinfo)) return 999;

    for (t=0; t<n; t++) (*pZ)[v][t] = (double) t+1;
    strcpy(pdinfo->varname[v], "time");
    strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
	    
    return index;
}

/* ...................................................................  */

static int diffvarnum (int index, const DATAINFO *pdinfo)
     /* Given an "ordinary" variable name, construct the name of the
	corresponding first difference and find its ID number */
{
    char diffname[16], s[16];
    
    strcpy(s, pdinfo->varname[index]);
    _esl_trunc(s, 6);
    strcpy(diffname, "d_");
    strcat(diffname, s);
    return varindex(pdinfo, diffname);
}

/* ......................................................  */

static int diffgenr (int iv, double ***pZ, DATAINFO *pdinfo)
{
    char word[32];
    char s[32];
    int t, t1, n = pdinfo->n, v = pdinfo->v;
    double x0, x1;

    strcpy(word, pdinfo->varname[iv]);
    _esl_trunc(word, 6);
    strcpy(s, "d_");
    strcat(s, word);

    /* "s" should now contain the new variable name --
     check whether it already exists: if so, get out */
    if (varindex(pdinfo, s) < v) return 0;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->time_series == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}	
	x0 = (*pZ)[iv][t];
	x1 = (*pZ)[iv][t-1];
	if (na(x0) || na(x1)) 
	    (*pZ)[v][t] = NADBL;
	else				      
	    (*pZ)[v][t] = x0 - x1;
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(VARLABEL(pdinfo, v), _("%s = first difference of %s"),
	    pdinfo->varname[v], pdinfo->varname[iv]);
	    
    return 0;
}

/* ......................................................  */

static int ldiffgenr (int iv, double ***pZ, DATAINFO *pdinfo)
{
    char word[32];
    char s[32];
    int t, t1, n = pdinfo->n, v = pdinfo->v;
    double x0, x1;

    strcpy(word, pdinfo->varname[iv]);
    _esl_trunc(word, 5);
    strcpy(s, "ld_");
    strcat(s, word);

    /* "s" should now contain the new variable name --
     check whether it already exists: if so, get out */
    if (varindex(pdinfo, s) < v) return 0;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->time_series == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}
	x0 = (*pZ)[iv][t];
	x1 = (*pZ)[iv][t-1];
	if (na(x0) || na(x1) || x0/x1 < 0.) 
	    (*pZ)[v][t] = NADBL;
	else				      
	    (*pZ)[v][t] = log(x0/x1);
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(VARLABEL(pdinfo, v), _("%s = log difference of %s"),
	    pdinfo->varname[v], pdinfo->varname[iv]);
	    
    return 0;
}

/**
 * list_diffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate first-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_diffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++)
	if (diffgenr(list[i], pZ, pdinfo)) return 1;
    return 0;
}

/**
 * list_ldiffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate log-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_ldiffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++)
	if (ldiffgenr(list[i], pZ, pdinfo)) return 1;
    return 0;
}

/**
 * _lagvarnum:
 * @iv: ID number of the variable.
 * @lag: Desired lag length.
 * @pdinfo: data information struct.
 *
 * Given an "ordinary" variable, construct the name of the
 * corresponding lagged variable and find its ID number.
 *
 * Returns: the ID number of the lagged variable.
 *
 */

int _lagvarnum (int iv, int lag, const DATAINFO *pdinfo)
{
    char lagname[16], s[4];
    
    strcpy(lagname, pdinfo->varname[iv]);
    if (pdinfo->pd >=10) _esl_trunc(lagname, 5);
    else _esl_trunc(lagname, 6);
    sprintf(s, "_%d", lag);
    strcat(lagname, s);
    return varindex(pdinfo, lagname);
}

/* ...................................................................  */

static void reset_list (int *list1, int *list2)
{
    int i;
    
    for (i=2; i<=list1[0]; i++) list1[i] = list2[i];
}

/* ...................................................................  */

static int get_listlen (const int *varlist, int order, double **Z, 
			const DATAINFO *pdinfo)
     /* parse varlist (for a VAR) and determine how long the augmented 
	list will be, once all the appropriate lag terms are inserted */
{
    int i, v = 1;
    
    for (i=1; i<=varlist[0]; i++) {
	if (strcmp(pdinfo->varname[varlist[i]], "time") == 0 ||
	    strcmp(pdinfo->varname[varlist[i]], "const") == 0 ||
	    isdummy(Z[varlist[i]], pdinfo->t1, pdinfo->t2)) v++;
	else v += order;
    }
    return v;
}

static int real_var (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
		     PRN *prn, struct var_resids *resids, char flags)
{
    /* construct the respective lists by adding the appropriate
       number of lags ("order") to the variables in list 

       Say the list is "x_1 const time x_2 x_3", and the order is 2.
       Then the first list should be

       x_1 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

       the second:

       x_2 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

       and so on.

       Run the regressions and print the results.
    */

    int i, j, index, l, listlen, end, neqns = 0;
    int *varlist = NULL, *depvars = NULL, *shortlist = NULL;
    int t1, t2, oldt1, oldt2, dfd = 0;
    int missv = 0, misst = 0;
    double essu = 0.0, F = 0.0;
    MODEL var_model;

    _init_model(&var_model, pdinfo);

    if (resids == NULL && order < 1) {
	fprintf(stderr, _("Not much point in a zero-order \"VAR\" surely?\n"));
	return 1;
    }

    /* how long will our list have to be? */
    listlen = get_listlen(list, order, *pZ, pdinfo);

    varlist = malloc((listlen + 1) * sizeof *varlist);
    depvars = malloc((listlen + 1) * sizeof *depvars);
    shortlist = malloc(listlen * sizeof *shortlist);
    if (varlist == NULL || depvars == NULL || shortlist == NULL)
	return E_ALLOC;

    varlist[0] = listlen;
    index = 2; /* skip beyond the counter and the dep var */
    end = listlen;

    /* now fill out the list */
    for (i=1; i<=list[0]; i++) {
	/* if we're looking at a dummy variable (or time trend) just include 
	   it unmodified -- at the end of the list */
	if (strcmp(pdinfo->varname[list[i]], "time") == 0 || 
	    strcmp(pdinfo->varname[list[i]], "const") == 0 || 
	    isdummy((*pZ)[list[i]], pdinfo->t1, pdinfo->t2)) { 
	    varlist[end] = list[i];
	    end--;
	    continue;	    
	}
	/* otherwise it's a "real" variable and we replace it with
	   <order> lags of itself */
	if (varindex(pdinfo, pdinfo->varname[list[i]]) 
	    < pdinfo->v) {
	    depvars[neqns] = list[i];
	    neqns++;
	    for (l=1; l<=order; l++) {
		_laggenr(list[i], l, 1, pZ, pdinfo);
		/* note: the lagvar may already exist */
		varlist[index] = _lagvarnum(list[i], l, pdinfo); 
		index++;
	    }
	}
    }

    /* sort out sample range */
    t1 = oldt1 = pdinfo->t1;
    t2 = oldt2 = pdinfo->t2;
    varlist[1] = depvars[0];
    if ((missv = _adjust_t1t2(NULL, varlist, &t1, &t2, *pZ, &misst))) {
	free(varlist);
	free(depvars);
	free(shortlist);
	return 1;
    }
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    /* run and print out the several regressions */
    if (flags & VAR_PRINT_MODELS) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), order);
    }
    shortlist[0] = listlen - order;

    /* apparatus for saving the residuals */
    if (resids != NULL) {
	resids->m = 2 * neqns;
	resids->uhat = malloc(2 * neqns * sizeof(double *));
	if (resids->uhat == NULL) return E_ALLOC;
    }

    for (i=0; i<neqns; i++) {
	varlist[1] = depvars[i];
	/* run an OLS regression for the current dep var */
	var_model = lsq(varlist, pZ, pdinfo, VAR, 0, 0.0);
	var_model.aux = VAR;
	/* save the residuals if required */
	if (resids != NULL) {
	    resids->t1 = var_model.t1;
	    resids->t2 = var_model.t2;
	    resids->uhat[i] = var_model.uhat;
	    var_model.uhat = NULL;
	} 
	if (flags & VAR_PRINT_MODELS) {
	    printmodel(&var_model, pdinfo, prn);
	}
	if (flags & VAR_DO_FTESTS) {
	    /* keep some results for hypothesis testing */
	    essu = var_model.ess;
	    dfd = var_model.dfd;	    
	}
	clear_model(&var_model, pdinfo);
	if (resids != NULL) {
	    /* estimate equations for Johansen test too */
	    varlist[1] = resids->levels_list[i + 1]; 
	    var_model = lsq(varlist, pZ, pdinfo, VAR, 0, 0.0);
	    if (flags & VAR_PRINT_MODELS) {
		var_model.aux = VAR;
		printmodel(&var_model, pdinfo, prn);
	    }
	    resids->uhat[i + neqns] = var_model.uhat;
	    var_model.uhat = NULL;
	    clear_model(&var_model, pdinfo);
	}
	if (flags & VAR_DO_FTESTS) {
	    /* now build truncated lists for hyp. tests */
	    shortlist[1] = varlist[1];
	    pputs(prn, _("\nF-tests of zero restrictions:\n\n"));
	    for (j=0; j<neqns; j++) {
		reset_list(shortlist, varlist);
		for (l=1; l<=order; l++) {
		    index = l + 1 + j * order;
		    if (index > shortlist[0]) break;
		    shortlist[index] = varlist[index+order];
		}
		end = 0;
		for (l=shortlist[0]; l>index; l--) {
		    shortlist[l] = varlist[varlist[0]-end];
		    end++;
		}
		pprintf(prn, _("All lags of %-8s "), 
			pdinfo->varname[depvars[j]]);
		/*  printlist(shortlist); */
		var_model = lsq(shortlist, pZ, pdinfo, VAR, 0, 0.0);
		F = ((var_model.ess - essu)/order)/(essu/dfd);
		clear_model(&var_model, pdinfo);
		pprintf(prn, "F(%d, %d) = %f, ", order, dfd, F);
		pprintf(prn, _("p-value %f\n"), fdist(F, order, dfd));
	    }
	    if (order > 1) {
		pprintf(prn, _("All vars, lag %-6d "), order);
		reset_list(shortlist, varlist);
		index = 2;
		for (j=1; j<=neqns*(order); j++) {
		    if (j % order) {
			shortlist[index] = varlist[j+1];
			index++;
		    }
		}
		end = 0;
		for (l=shortlist[0]; l>=index; l--) {
		    shortlist[l] = varlist[varlist[0]-end];
		    end++;
		}
		/*  printlist(shortlist); */
		var_model = lsq(shortlist, pZ, pdinfo, VAR, 0, 0.0);
		F = ((var_model.ess - essu)/neqns)/(essu/dfd);
		clear_model(&var_model, pdinfo);
		pprintf(prn, "F(%d, %d) = %f, ", neqns, dfd, F);
		pprintf(prn, _("p-value %f\n"), fdist(F, neqns, dfd)); 
	    }
	    pputs(prn, "\n");
	    if (flags & VAR_PRINT_PAUSE) page_break(0, NULL, 0);
	}
    }
    pputs(prn, "\n");

    free(varlist);
    free(shortlist);
    free(depvars);

    /* reset sample range to what it was before */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;
    return 0;
}


/**
 * var:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @pause: if = 1, pause after showing each model.
 * @prn: gretl printing struct.
 *
 * Estimate a vector auto-regression (VAR) and print the results.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int var (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
	 int pause, PRN *prn)
{
    char flags = VAR_PRINT_MODELS | VAR_DO_FTESTS;

    if (pause) {
	flags |= VAR_PRINT_PAUSE;
    }
    return real_var(order, list, pZ, pdinfo, prn, NULL, flags);
}

/**
 * coint:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Test for cointegration.  
 *
 * Returns: 0 on successful completion.
 *
 */

int coint (int order, const LIST list, double ***pZ, 
	   DATAINFO *pdinfo, PRN *prn)
     /* FIXME - need proper error checking here */
{
    int i, t, n, nv, l0 = list[0];
    MODEL coint_model;
    int *cointlist = NULL;

    _init_model(&coint_model, pdinfo);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (i > 1) pputs(prn, "\n");
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	adf_test(order, list[i], pZ, pdinfo, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (_hasconst(list) == 0) {
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) return E_ALLOC;
	for (i=0; i<=l0; i++) cointlist[i] = list[i];
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else 
	copylist(&cointlist, list);

    pputs(prn, "\n");
    pprintf(prn, _("Step %d: cointegration\n"), l0 + 1);
    
    coint_model = lsq(cointlist, pZ, pdinfo, OLS, 1, 0.0); 
    coint_model.aux = AUX_COINT;
    printmodel(&coint_model, pdinfo, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
    nv = pdinfo->v - 1;
    for (t=0; t<coint_model.t1; t++)
	(*pZ)[nv][t] = NADBL;
    for (t = coint_model.t1; t<=coint_model.t2; t++)
	(*pZ)[nv][t] = coint_model.uhat[t];
    for (t=coint_model.t2 + 1; t<n; t++) 
	(*pZ)[nv][t] = NADBL;
    strcpy(pdinfo->varname[nv], "uhat");

    /* Run ADF test on these residuals */
    pputs(prn, "\n");
    adf_test(order, pdinfo->v - 1, pZ, pdinfo, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
	    "(a) The unit-root hypothesis is not rejected for the individual"
	    " variables.\n(b) The unit-root hypothesis is rejected for the "
	    "residuals (uhat) from the \n    cointegrating regression.\n"
	    "\n(Note that significance levels for the D-W and F statistics here "
	    "cannot be \nread from the usual statistical tables.)\n"));

    /* clean up and get out */
    clear_model(&coint_model, pdinfo);
    free(cointlist);
    dataset_drop_vars(1, pZ, pdinfo);
    return 0;
}

/**
 * adf_test:
 * @order: lag order for the test.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller test for 
 * a unit root.
 *
 * Returns: 0 on successful completion, non-zero on error.
 *
 */

int adf_test (int order, int varno, double ***pZ,
	      DATAINFO *pdinfo, PRN *prn)
{
    int i, l, T, k, row, orig_nvars = pdinfo->v;
    int *adflist;
    int *shortlist;
    MODEL adf_model;
    double essu, F, DFt;
    char pval[40];

                                 /* 1%    2.5%    5%    10% */
    double t_crit_vals[6][4] = {{ -3.75, -3.33, -3.00, -2.62 },  /* T=25 */
				{ -3.58, -3.22, -2.93, -2.60 },  /* T=50 */
				{ -3.51, -3.17, -2.89, -2.58 },  /* T=100 */
				{ -3.46, -3.14, -2.88, -2.57 },  /* T=250 */
				{ -3.44, -3.13, -2.87, -2.57 },  /* T=500 */
				{ -3.43, -3.12, -2.86, -2.57 }}; /* inf */

                                /* 1%   2.5%   5%    10% */
    double f_crit_vals[6][4] = {{ 5.91, 7.24, 8.65, 10.61 },  /* T = 25 */
 			        { 5.61, 6.73, 7.81,  9.31 },  /* T = 50 */
			        { 5.47, 6.49, 7.44,  8.73 },  /* T = 100 */
			        { 5.39, 6.34, 7.25,  8.43 },  /* T = 250 */
			        { 5.36, 6.30, 7.20,  8.34 },  /* T = 500 */
			        { 5.34, 6.25, 7.16,  8.27 }}; /* inf */
    

    if (varno == 0) return E_DATA;

    _init_model(&adf_model, pdinfo);
    k = 3 + order;
    adflist = malloc((5 + order) * sizeof *adflist);
    shortlist = malloc(k * sizeof *shortlist);
    if (adflist == NULL || shortlist == NULL) return E_ALLOC;

    i = pdinfo->t1;
    pdinfo->t1 = 0;
    diffgenr(varno, pZ, pdinfo);
    _laggenr(varno, 1, 1, pZ, pdinfo);
    pdinfo->t1 = i;

    adflist[1] = diffvarnum(varno, pdinfo);

    /* do the more familiar Dickey-Fuller t-test first */
    adflist[0] = 3;
    adflist[2] = 0;
    adflist[3] = _lagvarnum(varno, 1, pdinfo);

    adf_model = lsq(adflist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode) {
	return adf_model.errcode;
    }

    DFt = adf_model.coeff[1] / adf_model.sderr[1];
    T = adf_model.nobs;

    row = (T > 500)? 5 : 
	(T > 450)? 4 : 
	(T > 240)? 3 : 
	(T > 90)? 2 : 
	(T > 40)? 1 : 
	(T > 24)? 0 : -1;

    if (row < 0) {
	sprintf(pval, _("significance level unknown"));
    } else {
	if (DFt < t_crit_vals[row][0])
	    sprintf(pval, _("significant at the 1 percent level"));
	else if (DFt < t_crit_vals[row][1])
	    sprintf(pval, _("significant at the 2.5 percent level"));
	else if (DFt < t_crit_vals[row][2])
	    sprintf(pval, _("significant at the 5 percent level"));
	else if (DFt < t_crit_vals[row][3])
	    sprintf(pval, _("significant at the 10 percent level"));
	else
	    sprintf(pval, _("not significant at the 10 percent level"));
    }
    
    pprintf(prn, _("\nDickey-Fuller test with constant\n\n"
	    "   model: (1 - L)%s = m + g * %s(-1) + e\n"
	    "   unit-root null hypothesis: g = 0\n"
	    "   estimated value of g: %f\n"
	    "   test statistic: t = %f, with sample size %d\n"
	    "   %s\n"),
	    pdinfo->varname[varno], pdinfo->varname[varno],
	    adf_model.coeff[1], DFt, adf_model.nobs, pval);
    clear_model(&adf_model, pdinfo);

    /* then do ADF test using F-statistic */
    adflist[0] = 4 + order;
    adflist[3] = _lagvarnum(varno, 1, pdinfo);

    for (l=1; l<=order; l++) {
	_laggenr(adflist[1], l, 1, pZ, pdinfo);
	/* note: the lagvar may already exist */
	adflist[l+3] = _lagvarnum(adflist[1], l, pdinfo); 
    }

    adflist[adflist[0]] = 0;
    if ((adflist[2] = gettrend(pZ, pdinfo)) == 999) {
	free(adflist);
	free(shortlist);
	return E_ALLOC;
    }
    /*  printlist(adflist); */
    adf_model = lsq(adflist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode)
	return adf_model.errcode;
    adf_model.aux = AUX_ADF;
    printmodel(&adf_model, pdinfo, prn);
    essu = adf_model.ess;
    T = adf_model.nobs;
    clear_model(&adf_model, pdinfo);

    shortlist[0] = adflist[0] - 2;
    shortlist[1] = adflist[1];
    for (i=0; i<=order; i++) 
	shortlist[2+i] = adflist[4+i];
    /*  printlist(shortlist); */
    adf_model = lsq(shortlist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode)
	return adf_model.errcode;	
    F = (adf_model.ess - essu) * (T - k)/(2 * essu);
    clear_model(&adf_model, pdinfo);

    row = -1;
    if (T > 500) row = 5;
    else if (T > 250) row = 4;
    else if (T > 100) row = 3;
    else if (T > 50) row = 2;
    else if (T > 25) row = 1;
    else if (T > 23) row = 0;

    if (row == -1) strcpy(pval, _("unknown pvalue"));
    else {
	if (F > f_crit_vals[row][3]) strcpy(pval, _("pvalue < .01"));
	else if (F > f_crit_vals[row][2]) strcpy(pval, _(".025 > pvalue > .01"));
	else if (F > f_crit_vals[row][1]) strcpy(pval, _(".05 > pvalue > .025"));
	else if (F > f_crit_vals[row][0]) strcpy(pval, _(".10 > pvalue > .05"));
	else strcpy(pval, _("pvalue > .10"));
    }

    pprintf(prn, _("Augmented Dickey-Fuller test on %s:\n   F(2, %d) = %f, "
	   "with %s\n"), pdinfo->varname[varno], T - k, F, pval);
    pprintf(prn, _("The null hypothesis is that %s has a unit root, i.e. "
	    "the parameters on\nthe time trend and %s are both zero.\n"),
	    pdinfo->varname[varno], pdinfo->varname[adflist[3]]);

    free(adflist);
    free(shortlist);
    dataset_drop_vars(pdinfo->v - orig_nvars, pZ, pdinfo);

    return 0;
}

/* ....................................................... */

#ifdef notyet

int ma_model (LIST list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int t, v = pdinfo->v, err = 0;
    int malist[4], iv = list[2];
    double a, aopt, essmin, diff;
    int step, t0 = pdinfo->t1, T = pdinfo->t2;
    MODEL mamod;

    if (list[0] != 2) {
	pputs(prn, "mvavg: takes a list of two variables\n");
	return 1;
    }
    
    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
    strcpy(pdinfo->varname[v], "Z_t");

    malist[0] = 3;
    malist[1] = list[1]; /* original dependent variable */
    malist[2] = v;       /* new var: moving average of indep var */
    malist[3] = 0;

    _init_model(&mamod, pdinfo);

    a = aopt = 0.0;
    essmin = 0.0;
    diff = 0.01;
    for (step=1; step<=100; step++) {
	a += diff;
	if (a > 0.995) break;
	(*pZ)[v][t0] = (*pZ)[iv][t0] / (1 - a);
	for (t=t0+1; t<T; t++) { 
	    (*pZ)[v][t] = (*pZ)[iv][t] + a * (*pZ)[v][t-1];
	    /*  printf("newvars[%d] %g %g\n", t, 
		(*pZ)[v*n + t], (*pZ)[(v+1)*n + t]); */
	}
	clear_model(&mamod, pdinfo);
	mamod = lsq(malist, pZ, pdinfo, OLS, 0, 0.0);
	if ((err = mamod.errcode)) {
	    clear_model(&mamod, pdinfo);
	    return err;
	}	
	if (step == 1) {
	    pputs(prn, "\n ADJ       ESS      ADJ       ESS      "
		    "ADJ       ESS      ADJ       ESS     \n");
	}
	pprintf(prn, "%5.2f %10.4g", a, mamod.ess);
	if (step%4 == 0) pputs(prn, "\n");
	else _bufspace(3, prn);
	if (step == 1 || mamod.ess < essmin) {
	    essmin = mamod.ess;
	    aopt = a;
	}
    }

    pprintf(prn, "\n\nESS is minimum for adj = %.2f\n\n", aopt);
    a = aopt;
    (*pZ)[v][t0] = (*pZ)[iv][t0] / (1 - a);
    for (t=t0+1; t<T; t++) { 
	(*pZ)[v][t] = (*pZ)[iv][t] + a * (*pZ)[v][t-1];
    }
    mamod = lsq(malist, pZ, pdinfo, OLS, 1, 0.0);
    printmodel(&mamod, pdinfo, prn);

    pputs(prn, "\nEstimates of original parameters:\n");
    pprintf(prn, "constant: %.4g\n", mamod.coeff[0]);
    pprintf(prn, "slope:    %.4g\n", mamod.coeff[1] / (1 - a));
    pprintf(prn, "adaptive coefficient: %.2f\n", a);
	   
    clear_model(&mamod, pdinfo);

    dataset_drop_vars(1, pZ, pdinfo);

    return 0;
}

#endif /* notyet */

static int 
has_time_trend (LIST varlist, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    int origv = pdinfo->v;
    int tlist[4];
    int trends = 0;
    MODEL tmod;

    _init_model(&tmod, pdinfo);

    tlist[0] = 3;
    tlist[2] = 0;

    for (i=1; i<=varlist[0]; i++) {
	double tstat;
	int v, vl;

	v = varlist[i];
	if (v == 0) continue;
	if (diffgenr(v, pZ, pdinfo)) {
	    trends = -1;
	    break;
	}
	vl = _lagvarnum(v, 1, pdinfo);
	tlist[1] = v;
	tlist[3] = vl;
	tmod = lsq(tlist, pZ, pdinfo, OLS, 0, 0.0);
	if (tmod.errcode) {
	    trends = -1;
	    clear_model(&tmod, pdinfo);
	    break;
	}
	tstat = tmod.coeff[0] / tmod.sderr[0];
	if (tprob(tstat, tmod.dfd) < 0.05) {
	    trends = 1;
	}
	clear_model(&tmod, pdinfo);
	if (trends) break;
    }

    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

    return trends;
}

static int allocate_sigmas (double ***X, double ***Y, double ***Z, int k)
{
    int i, j;
    double **Suu, **Svv, **Suv;

    Suu = malloc(k * sizeof *Suu);
    Svv = malloc(k * sizeof *Svv);
    Suv = malloc(k * sizeof *Suv);

    if (Suu == NULL || Svv == NULL || Suv == NULL) return 1;

    for (i=0; i<k; i++) {
	Suu[i] = malloc(k * sizeof **Suu);
	Svv[i] = malloc(k * sizeof **Svv);
	Suv[i] = malloc(k * sizeof **Suv);
	if (Suu[i] == NULL || Svv[i] == NULL || Suv[i] == NULL) {
	    free(Suu);
	    free(Svv);
	    free(Suv);
	    return 1;
	}
	for (j=0; j<k; j++) {
	    Suu[i][j] = 0.0;
	    Svv[i][j] = 0.0;
	    Suv[i][j] = 0.0;
	}
    }

    *X = Suu;
    *Y = Svv;
    *Z = Suv;

    return 0;
}

static void free_sigmas (double **X, double **Y, double **Z, int k)
{
    int i;

    for (i=0; i<k; i++) {
	if (X != NULL) free(X[i]);
	if (Y != NULL) free(Y[i]);
	if (Z != NULL) free(Z[i]);
    }

    free(X);
    free(Y);
    free(Z);
}

static void scatter_product (const double **u, const double **v, 
			     double **X, int T, int k)
{
    int i, j, t;

    for (t=0; t<T; t++) {
	for (i=0; i<k; i++) {
	    for (j=0; j<k; j++) {
		X[i][j] += u[i][t] * v[j][t];
	    }
	}
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    X[i][j] /= (double) T;
	}
    }
}

static void 
print_sigmas (const double **X, const double **Y, const double **Z, 
	      int k, PRN *prn)
{
    int i, j, l;
    const double **P = NULL;

    pprintf(prn, "\n%s\n\n", _("Sample variance-covariance matrices for residuals"));

    for (l=0; l<3; l++) {
	if (l == 0) {
	    P = X;
	    pprintf(prn, " %s\n\n", _("VAR system in first differences"));
	}
	else if (l == 1) {
	    P = Y;
	    pprintf(prn, " %s\n\n", _("System with levels as dependent variable"));
	}
	else {
	    P = Z;
	    pprintf(prn, " %s\n\n", _("Cross-products"));
	}
	for (i=0; i<k; i++) {
	    for (j=0; j<k; j++) {
		pprintf(prn, "%#12.6g", P[i][j]);
	    }
	    pputs(prn, "\n");
	}
	pputs(prn, "\n");
    }
}

static int 
johansen_complete (const double **X, const double **Y, const double **Z,
		   int k, int T, int trends, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (const double **, const double **, const double **,
		     int, int, int, PRN *);
    int err = 0;

    *gretl_errmsg = 0;
    
    if (open_plugin("johansen", &handle)) {
        err = 1;
        strcpy(gretl_errmsg, _("Couldn't load plugin function\n"));
        goto system_bailout;
    }

    johansen = get_plugin_function("johansen_eigenvals", handle);
    if (johansen == NULL) {
        err = 1;
        strcpy(gretl_errmsg, _("Couldn't load plugin function\n"));
        goto system_bailout;
    }
        
    err = (* johansen) (X, Y, Z, k, T, trends, prn);
    
 system_bailout:
    if (handle != NULL) {
        close_plugin(handle);
    }

    return err;
}

int johansen_test (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
		   int verbose, PRN *prn)
{
    PRN *varprn = NULL;
    struct var_resids resids;
    char flags;
    int err = 0;
    int i, j;
    int orig_t1 = pdinfo->t1;
    int orig_v = pdinfo->v;
    int *varlist;
    int hasconst = 0;
    int trends = 0;

    /* we're assuming that the list we are fed is in levels */
    resids.levels_list = malloc((1 + list[0]) * sizeof *list);
    if (resids.levels_list == NULL) return E_ALLOC;
    resids.levels_list[0] = list[0];

    varlist = malloc((2 + list[0]) * sizeof *list);
    if (varlist == NULL) return E_ALLOC;
    varlist[0] = list[0];

    j = 1;
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    resids.levels_list[0] -= 1;
	    hasconst = 1;
	    continue;
	}
	_laggenr(list[i], 1, 1, pZ, pdinfo);
	resids.levels_list[j++] = _lagvarnum(list[i], 1, pdinfo);
    }

    /* now get differences and put them into list */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	diffgenr(list[i], pZ, pdinfo);
	varlist[i] = diffvarnum(list[i], pdinfo);
    }

    if (!hasconst) {
	varlist[0] += 1;
	varlist[varlist[0]] = 0;
    }

    if (verbose) {
	flags = VAR_PRINT_MODELS;
	varprn = prn;
    } else {
	flags = 0;
	varprn = gretl_print_new(GRETL_PRINT_NULL, NULL);
    }

    /* FIXME? */
    pdinfo->t1 += (order + 1);
    /* Check Hamilton: what if order for test = 1? */
    err = real_var(order - 1, varlist, pZ, pdinfo, varprn, &resids, flags); 
    
    if (!verbose) {
	gretl_print_destroy(varprn);
    }

    if (!err) {
	int k = resids.m / 2;
	int T = resids.t2 - resids.t1 + 1;
	double **Suu, **Svv, **Suv;
	double **u = NULL, **v = NULL;
	char stobs[9], endobs[9];

	if (allocate_sigmas(&Suu, &Svv, &Suv, k)) {
	    err = E_ALLOC;
	    goto johansen_bailout;
	}

	u = malloc(k * sizeof *u);
	v = malloc(k * sizeof *v);

	if (u == NULL || v == NULL) {
	    err = E_ALLOC;
	    goto johansen_bailout;
	}

	for (i=0; i<k; i++) {
	    u[i] = &(resids.uhat[i][resids.t1]);
	    v[i] = &(resids.uhat[i + k][resids.t1]);
	}

	scatter_product((const double **) u, (const double **) u, Suu, T, k);
	scatter_product((const double **) v, (const double **) v, Svv, T, k);
	scatter_product((const double **) u, (const double **) v, Suv, T, k);

	pprintf(prn, "%s:\n", _("Johansen test"));
	pprintf(prn, "%s = %d\n", _("Number of equations"), k);
	pprintf(prn, "%s: %s - %s (T = %d)\n", _("Estimation period"),
		ntodate(stobs, resids.t1, pdinfo), 
		ntodate(endobs, resids.t2, pdinfo), T);

	if (verbose) {
	    print_sigmas((const double **) Suu, 
			 (const double **) Svv, 
			 (const double **) Suv, k, prn);
	}

#ifdef JOHANSEN_DEBUG
	for (i=0; i<resids.m; i++) {
	    int t;

	    pprintf(prn, "Residuals from VAR model %d\n", i);
	    for (t=resids.t1; t<=resids.t2; t++) {
		ntodate(date, t, pdinfo);
		pprintf(prn, "%8s %#.*g\n", ntodate(stobs, t, pdinfo), 
			GRETL_DIGITS, resids.uhat[i][t]);
	    }
	}
#endif

	trends = has_time_trend(list, pZ, pdinfo);
	if (trends == -1) {
	    pprintf(prn, "%s\n", _("Error checking for time trends"));
	    goto johansen_bailout;
	}

	/* now get johansen plugin to finish the job */
	err = johansen_complete((const double **) Suu, 
				(const double **) Svv, 
				(const double **) Suv, k, T, trends, prn);

    johansen_bailout:
	
	for (i=0; i<resids.m; i++) {
	    free(resids.uhat[i]);
	}
	free(resids.uhat);

	free_sigmas(Suu, Svv, Suv, k);
	free(u);
	free(v);
    } 

    free(resids.levels_list);
    free(varlist);

    pdinfo->t1 = orig_t1;

    dataset_drop_vars(pdinfo->v - orig_v, pZ, pdinfo);

    return err;
}

