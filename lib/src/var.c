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

/* ......................................................  */

static int gettrend (double **pZ, DATAINFO *pdinfo)
{
    int index;
    int t, n = pdinfo->n, v = pdinfo->v;

    index = varindex(pdinfo, "time");
    if (index < v) return index;
    
    if (_grow_Z(1, pZ, pdinfo)) return 999;

    for (t=0; t<n; t++) (*pZ)[n*v + t] = (double) t+1;
    strcpy(pdinfo->varname[v], "time");
    strcpy(pdinfo->label[v], "time trend variable");
	    
    return index;
}

/* ...................................................................  */

static int diffvarnum (const int index, const DATAINFO *pdinfo)
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

static int diffgenr (const int iv, double **pZ, DATAINFO *pdinfo)
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

    if (_grow_Z(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[n*v + t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	x0 = (*pZ)[iv*n + t];
	x1 = (*pZ)[iv*n + t-1];
	if (na(x0) || na(x1)) 
	    (*pZ)[n*v + t] = NADBL;
	else				      
	    (*pZ)[n*v + t] = x0 - x1;
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(pdinfo->label[v], "%s = first difference of %s",
	    pdinfo->varname[v], pdinfo->varname[iv]);
	    
    return 0;
}

/* ......................................................  */

static int ldiffgenr (const int iv, double **pZ, DATAINFO *pdinfo)
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

    if (_grow_Z(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[n*v + t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	x0 = (*pZ)[iv*n + t];
	x1 = (*pZ)[iv*n + t-1];
	if (na(x0) || na(x1) || x0/x1 < 0.) 
	    (*pZ)[n*v + t] = NADBL;
	else				      
	    (*pZ)[n*v + t] = log(x0/x1);
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(pdinfo->label[v], "%s = log difference of %s",
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

int list_diffgenr (const LIST list, double **pZ, DATAINFO *pdinfo)
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

int list_ldiffgenr (const LIST list, double **pZ, DATAINFO *pdinfo)
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

int _lagvarnum (const int iv, const int lag, const DATAINFO *pdinfo)
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

static int get_listlen (const int *varlist, const int order, double *Z, 
			const DATAINFO *pdinfo)
     /* parse varlist (for a VAR) and determine how long the augmented 
	list will be, once all the appropriate lag terms are inserted */
{
    int i, v = 1;
    
    for (i=1; i<=varlist[0]; i++) {
	if (strcmp(pdinfo->varname[varlist[i]], "time") == 0 ||
	    strcmp(pdinfo->varname[varlist[i]], "const") == 0 ||
	    isdummy(varlist[i], pdinfo->t1, pdinfo->t2, Z, 
		    pdinfo->n)) v++;
	else v += order;
    }
    return v;
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

int var (const int order, const LIST list, double **pZ, DATAINFO *pdinfo,
	 const int pause, PRN *prn)
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
    int *varlist, *depvars, *shortlist;
    int t1, t2, oldt1, oldt2, dfd;
    int missv = 0, misst = 0;
    double essu, F;
    MODEL var_model;

    _init_model(&var_model);

    if (order < 1) {
	fprintf(stderr, "Not much point in a zero-order \"VAR\" surely?\n");
	return 1;
    }

    /* how long will our list have to be? */
    listlen = get_listlen(list, order, *pZ, pdinfo);

    varlist = malloc((listlen + 1) * sizeof(int));
    depvars = malloc((listlen + 1) * sizeof(int));
    shortlist = malloc(listlen * sizeof(int));
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
	    isdummy(list[i], pdinfo->t1, pdinfo->t2, *pZ, 
		    pdinfo->n)) {
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
    if ((missv = _adjust_t1t2(NULL, varlist, &t1, &t2, *pZ, pdinfo->n, 
			      &misst))) {
	free(varlist);
	free(depvars);
	free(shortlist);
	return 1;
    }
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    /* run and print out the several regressions */
    pprintf(prn, "\nVAR system, lag order %d\n\n", order);
    shortlist[0] = listlen - order;
    
    for (i=0; i<neqns; i++) {
	varlist[1] = depvars[i];
	/* run an OLS regression for the current dep var */
	var_model = lsq(varlist, pZ, pdinfo, VAR, 0, 0.0);
	var_model.aux = VAR;
	printmodel(&var_model, pdinfo, prn);
	/* keep some results for hypothesis testing */
	essu = var_model.ess;
	dfd = var_model.dfd;
	clear_model(&var_model, NULL, NULL);
	/* now build truncated lists for hyp. tests */
	shortlist[1] = varlist[1];
	pprintf(prn, "\nF-tests of zero restrictions:\n\n");
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
	    pprintf(prn, "All lags of %-8s ", 
		   pdinfo->varname[depvars[j]]);
	    /*  printlist(shortlist); */
	    var_model = lsq(shortlist, pZ, pdinfo, VAR, 0, 0.0);
	    F = ((var_model.ess - essu)/order)/(essu/dfd);
	    clear_model(&var_model, NULL, NULL);
	    pprintf(prn, "F(%d, %d) = %f, ", order, dfd, F);
	    pprintf(prn, "p-value %f\n", fdist(F, order, dfd));
	}
	if (order > 1) {
	    pprintf(prn, "All vars, lag %-6d ", order);
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
	    clear_model(&var_model, NULL, NULL);
	    pprintf(prn, "F(%d, %d) = %f, ", neqns, dfd, F);
	    pprintf(prn, "p-value %f\n", fdist(F, neqns, dfd)); 
	}
	pprintf(prn, "\n");
	if (pause) page_break(0, NULL, 0);
    }
    pprintf(prn, "\n");

    free(varlist);
    free(shortlist);
    free(depvars);

    /* reset sample range to what it was before */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;
    return 0;
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

int coint (const int order, const LIST list, double **pZ, 
	   DATAINFO *pdinfo, PRN *prn)
     /* FIXME - need proper error checking here */
{
    int i, t, n, nv, l0 = list[0];
    MODEL coint_model;
    int *cointlist;

    _init_model(&coint_model);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	pprintf(prn, "\n");
	adf_test(order, list[i], pZ, pdinfo, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (_hasconst(list) == 0) {
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) return E_ALLOC;
	for (i=0; i<=l0; i++) cointlist[i] = list[i];
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else copylist(&cointlist, list);
    
    coint_model = lsq(cointlist, pZ, pdinfo, OLS, 1, 0.0); 
    coint_model.aux = AUX_COINT;
    printmodel(&coint_model, pdinfo, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (_grow_Z(1, pZ, pdinfo)) return E_ALLOC;
    nv = pdinfo->v - 1;
    for (t=0; t<coint_model.t1; t++)
	(*pZ)[n*nv + t] = NADBL;
    for (t = coint_model.t1; t<=coint_model.t2; t++)
	(*pZ)[n*nv + t] = coint_model.uhat[t];
    for (t=coint_model.t2 + 1; t<n; t++) 
	(*pZ)[n*nv + t] = NADBL;
    strcpy(pdinfo->varname[nv], "uhat");

    /* Run ADF test on these residuals */
    pprintf(prn, "\n");
    adf_test(order, pdinfo->v - 1, pZ, pdinfo, prn);

    pprintf(prn, "\nThere is evidence for a cointegrating relationship if:\n"
	    "(a) The unit-root hypothesis is not rejected for the individual"
	    " variables.\n(b) The unit-root hypothesis is rejected for the "
	    "residuals (uhat) from the \n    cointegrating regression.\n"
	    "\n(Note that significance levels for the D-W and F statistics here "
	    "cannot be \nread from the usual statistical tables.)\n");

    /* clean up and get out */
    clear_model(&coint_model, NULL, NULL);
    free(cointlist);
    _shrink_Z(1, pZ, pdinfo);
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
 * Carries out and prints the results of the Augmented Dickey-Fuller test for a unit root.
 *
 * Returns: 0 on successful completion.
 *
 */

int adf_test (const int order, const int varno, double **pZ,
	      DATAINFO *pdinfo, PRN *prn)
{
    int i, l, T, k, row, orig_nvars = pdinfo->v;
    int *adflist;
    int *shortlist;
    MODEL adf_model;
    double essu, F, DFt;
    char pval[40];

                                /* 99%    97.5%  95%    90%    10%    5%    2.5%  1% */
    double t_crit_vals[6][8] = {{-3.75, -3.33, -3.00, -2.62, -0.37,  0.00, 0.34, 0.72}, /* T=25 */
				{-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, 0.29, 0.66}, /* T=50 */
				{-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, 0.26, 0.63}, /* T=100 */
				{-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, 0.24, 0.62}, /* T=250 */
				{-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, 0.24, 0.61}, /* T=500 */
				{-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, 0.23, 0.60}}; /* T>500 */

                            /* .100  .050  .025  .010 */
    double crit_vals[6][4] = {{5.91, 7.24, 8.65, 10.61}, /* T = 25 */
			      {5.61, 6.73, 7.81, 9.31},  /* T = 50 */
			      {5.47, 6.49, 7.44, 8.73},  /* T = 100 */
			      {5.39, 6.34, 7.25, 8.43},  /* T = 250 */
			      {5.36, 6.30, 7.20, 8.34},  /* T = 500 */
			      {5.34, 6.25, 7.16, 8.27}}; /* infinity */
    

    _init_model(&adf_model);
    k = 3 + order;
    adflist = malloc((5 + order) * sizeof(int));
    shortlist = malloc(k * sizeof(int));
    if (adflist == NULL || shortlist == NULL) return E_ALLOC;

    i = pdinfo->t1;
    pdinfo->t1 = 0;
    diffgenr(varno, pZ, pdinfo);
    _laggenr(varno, 1, 1, pZ, pdinfo);
    pdinfo->t1 = i;

    adflist[1] = diffvarnum(varno, pdinfo);

    /* do the more familiar Dickey-Fuller t-test first */
    adflist[0] = 3;
    adflist[2] = _lagvarnum(varno, 1, pdinfo);
    adflist[3] = 0;
    adf_model = lsq(adflist, pZ, pdinfo, OLS, 0, 0.0);
    DFt = adf_model.coeff[1] / adf_model.sderr[1];
    T = adf_model.nobs;
    row = (T > 500)? 5 : (T > 450)? 4 : (T > 240)? 3 : (T > 90)? 2 : 
	(T > 40)? 1 : (T > 24)? 0 : -1;
    if (row < 0) {
	sprintf(pval, "significance level unknown");
    } else {
	if (DFt < t_crit_vals[row][0] || DFt > t_crit_vals[row][7])
	    sprintf(pval, "significant at the 1 percent level");
	else if (DFt < t_crit_vals[row][1] || DFt > t_crit_vals[row][6])
	    sprintf(pval, "significant at the 2.5 percent level");
	else if (DFt < t_crit_vals[row][2] || DFt > t_crit_vals[row][5])
	    sprintf(pval, "significant at the 5 percent level");
	else if (DFt < t_crit_vals[row][3] || DFt > t_crit_vals[row][4])
	    sprintf(pval, "significant at the 10 percent level");
	else
	    sprintf(pval, "not significant at the 10 percent level");
    }
    
    pprintf(prn, "\nDickey-Fuller test with constant\n\n"
	    "   model: (1 - L)%s = m + g * %s(-1) + e\n"
	    "   unit-root null hypothesis: g = 0\n"
	    "   estimated value of g: %f\n"
	    "   test statistic: t = %f, with sample size %d\n"
	    "   %s\n",
	    pdinfo->varname[varno], pdinfo->varname[varno],
	    adf_model.coeff[1], DFt, adf_model.nobs, pval);
    clear_model(&adf_model, NULL, NULL);

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
    adf_model.aux = AUX_ADF;
    printmodel(&adf_model, pdinfo, prn);
    essu = adf_model.ess;
    T = adf_model.nobs;
    clear_model(&adf_model, NULL, NULL);

    shortlist[0] = adflist[0] - 2;
    shortlist[1] = adflist[1];
    for (i=0; i<=order; i++) 
	shortlist[2+i] = adflist[4+i];
    /*  printlist(shortlist); */
    adf_model = lsq(shortlist, pZ, pdinfo, OLS, 0, 0.0);
    F = (adf_model.ess - essu) * (T - k)/(2 * essu);
    clear_model(&adf_model, NULL, NULL);

    row = -1;
    if (T > 500) row = 5;
    else if (T > 250) row = 4;
    else if (T > 100) row = 3;
    else if (T > 50) row = 2;
    else if (T > 25) row = 1;
    else if (T > 23) row = 0;

    if (row == -1) strcpy(pval, "unknown pvalue");
    else {
	if (F > crit_vals[row][3]) strcpy(pval, "pvalue < .01");
	else if (F > crit_vals[row][2]) strcpy(pval, ".025 > pvalue > .01");
	else if (F > crit_vals[row][1]) strcpy(pval, ".05 > pvalue > .025");
	else if (F > crit_vals[row][0]) strcpy(pval, ".10 > pvalue > .05");
	else strcpy(pval, "pvalue > .10");
    }

    pprintf(prn, "Augmented Dickey-Fuller test on %s:\n   F(2, %d) = %f, "
	   "with %s\n", pdinfo->varname[varno], T - k, F, pval);
    pprintf(prn, "The null hypothesis is that %s has a unit root, i.e. "
	    "the parameters on\nthe time trend and %s are both zero.\n",
	    pdinfo->varname[varno], pdinfo->varname[adflist[3]]);

    free(adflist);
    free(shortlist);
    _shrink_Z(pdinfo->v - orig_nvars, pZ, pdinfo);
    return 0;
}

/* ....................................................... */

int ma_model (LIST list, double **pZ, DATAINFO *pdinfo, PRN *prn)
{
    int t, n = pdinfo->n, v = pdinfo->v, err = 0;
    int malist[4], iv = list[2];
    double a, aopt, essmin, diff;
    int step, t0 = pdinfo->t1, T = pdinfo->t2;
    MODEL mamod;

    if (list[0] != 2) {
	pprintf(prn, "mvavg: takes a list of two variables\n");
	return 1;
    }
    
    if (_grow_Z(1, pZ, pdinfo)) return E_ALLOC;
    strcpy(pdinfo->varname[v], "Z_t");

    malist[0] = 3;
    malist[1] = list[1]; /* original dependent variable */
    malist[2] = v;       /* new var: moving average of indep var */
    malist[3] = 0;

    _init_model(&mamod);

    a = aopt = 0.0;
    essmin = 0.0;
    diff = 0.01;
    for (step=1; step<=100; step++) {
	a += diff;
	if (a > 0.995) break;
	(*pZ)[v*n + t0] = (*pZ)[iv*n + t0] / (1 - a);
	for (t=t0+1; t<T; t++) { 
	    (*pZ)[v*n + t] = (*pZ)[iv*n + t] + a * (*pZ)[v*n + t-1];
	    /*  printf("newvars[%d] %g %g\n", t, 
		(*pZ)[v*n + t], (*pZ)[(v+1)*n + t]); */
	}
	clear_model(&mamod, NULL, NULL);
	mamod = lsq(malist, pZ, pdinfo, OLS, 0, 0.0);
	if ((err = mamod.errcode)) {
	    clear_model(&mamod, NULL, NULL);
	    return err;
	}	
	if (step == 1) {
	    pprintf(prn, "\n ADJ       ESS      ADJ       ESS      "
		    "ADJ       ESS      ADJ       ESS     \n");
	}
	pprintf(prn, "%5.2f %10.4g", a, mamod.ess);
	if (step%4 == 0) pprintf(prn, "\n");
	else _bufspace(3, prn);
	if (step == 1 || mamod.ess < essmin) {
	    essmin = mamod.ess;
	    aopt = a;
	}
    }

    pprintf(prn, "\n\nESS is minimum for adj = %.2f\n\n", aopt);
    a = aopt;
    (*pZ)[v*n + t0] = (*pZ)[iv*n + t0] / (1 - a);
    for (t=t0+1; t<T; t++) { 
	(*pZ)[v*n + t] = (*pZ)[iv*n + t] + a * (*pZ)[v*n + t-1];
    }
    mamod = lsq(malist, pZ, pdinfo, OLS, 1, 0.0);
    printmodel(&mamod, pdinfo, prn);

    pprintf(prn, "\nEstimates of original parameters:\n");
    pprintf(prn, "constant: %.4g\n", mamod.coeff[2]);
    pprintf(prn, "slope:    %.4g\n", mamod.coeff[1] / (1 - a));
    pprintf(prn, "adaptive coefficient: %.2f\n", a);
	   
    clear_model(&mamod, NULL, NULL);

    _shrink_Z(1, pZ, pdinfo);

    return 0;
}
