/*
 *  Copyright (c) 2002 by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* mp_ols.c - gretl least squares with multiple precision (GMP) */

#include "libgretl.h"
#include <float.h>
#include <gmp.h>

#ifdef HAVE_MPFR
#include <mpfr.h>
#endif

/* #define MP_DEBUG 1 */

#define DEFAULT_GRETL_MP_BITS 256   /* min. bits of precision for GMP */

static mpf_t MPF_ONE;
static mpf_t MPF_ZERO;
static mpf_t MPF_MINUS_ONE;

/* struct to hold model results */
typedef struct {
    int ID;                      /* ID number for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by (translated) ID number */
    int *varlist;                /* counterpart of list, using "real" ID #s */
    const int *polylist;         /* list of polynomial powers */
    int ifc;                     /* = 1 if the equation includes a constant,
				    else = 0 */
    mpf_t *coeff;                /* array of coefficient estimates */
    mpf_t *sderr;                /* array of estimated std. errors */
    mpf_t *xpx;
    mpf_t ess, tss;              /* Error and Total Sums of Squares */
    mpf_t sigma;                 /* Standard error of regression */
    mpf_t rsq, adjrsq;           /* Unadjusted and adjusted R^2 */     
    mpf_t fstt;                  /* F-statistic */
    int errcode;                 /* Error code in case of failure */
    int polyvar;                 /* number of the variable to be raised to
				    specified powers, if any */
} MPMODEL;

typedef struct {
    mpf_t *xpx;
    mpf_t *xpy;
    int ivalue;
    int nv;
    int errcode;
} MPXPXXPY;	

typedef struct {
    MPXPXXPY xpxxpy;
    mpf_t *coeff;
    mpf_t rss;
    int errcode;
} MPCHOLBETA;

static void set_gretl_mp_bits (void);
#ifdef HAVE_MPFR
static void set_gretl_mpfr_bits (void);
#endif
static MPXPXXPY mp_xpxxpy_func (const int *list, int n, mpf_t **mpZ);
static void mp_regress (MPMODEL *pmod, MPXPXXPY xpxxpy, mpf_t **mpZ, int n,
			char *errbuf);
static MPCHOLBETA mp_cholbeta (MPXPXXPY xpxxpy);
static void mp_diaginv (MPXPXXPY xpxxpy, mpf_t *diag);
static int mp_rearrange (int *list);

static void mpf_constants_init (void)
{
    mpf_init_set_d (MPF_ONE, 1.0);
    mpf_init_set_d (MPF_ZERO, 0.0);
    mpf_init_set_d (MPF_MINUS_ONE, -1.0);
}

static void mpf_constants_clear (void)
{
    mpf_clear (MPF_ONE);
    mpf_clear (MPF_ZERO);
    mpf_clear (MPF_MINUS_ONE);
}

static int data_problems (const int *list, double **Z, DATAINFO *pdinfo,
			  char *errbuf)
{
    /* reject (a) if any values are missing, (b) if any vars are all
       zero */
    int i, t, allzero;

    for (i=1; i<=list[0]; i++) {
	/* no need to check the constant */
	if (list[i] == 0) continue; 
	allzero = 1;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	    if (na(Z[list[i]][t])) {
		sprintf(errbuf, _("Missing observations for variable '%s'"), 
			pdinfo->varname[list[i]]);	
		return 1;
	    }
	    if (!floateq(Z[list[i]][t], 0.0)) allzero = 0;
	}
	if (allzero) {
	    sprintf(errbuf, _("Variable '%s' is all zeros"), 
		    pdinfo->varname[list[i]]);
	    return 1;
	}
    }
    return 0;
}

/* Given a data set containing doubles, build a data set using
   GMP's mpf_t floating point type.  

   For "ordinary" data we either (a) simply initialize the mpf_t
   straight from the corresponding double, or (b) attempt a clever
   trick, namely print the double to a string using a precision that
   was recorded at the time the original data were read, then set the
   mpf_t from that string.  This is designed to avoid the transmission
   to the mpf_t of garbage lying beyond DBL_DIGITS into the double.
   The trick is applicable only for data read from some original
   source, e.g. the NIST data files; also, it works only if the
   original data had a precision of not more than DBL_DIGITS.  

   Besides converting ordinary data, this function is also used to
   generate powers of x in the case of a polynomial regression of y on
   x, as found in several of the NIST examples.  Accuracy of the
   regression results may suffer if the generation of successive
   powers of x is done in regular double precision, so we do it here.
*/

static mpf_t **make_mpZ (MPMODEL *pmod, double **Z, DATAINFO *pdinfo)
{
    int i, s, t, n = pmod->t2 - pmod->t1 + 1;
    int l0 = pmod->list[0];
    int npoly = (pmod->polylist == NULL)? 0 : pmod->polylist[0];
    int mp_poly_pos = 0;
    int dropc = 0, nvars = 0, listpt;
    mpf_t **mpZ = NULL;
    unsigned char **digits = (unsigned char **) pdinfo->data;
    char numstr[64];

    if (n <= 0) return NULL;

    pmod->varlist = malloc ((l0 + 1) * sizeof(int));
    if (pmod->varlist == NULL) return NULL;
    pmod->varlist[0] = l0;

    mpZ = malloc(l0 * sizeof *mpZ);
    if (mpZ == NULL) return NULL;

    if (pmod->ifc) { /* create the constant vector */
	mpZ[0] = malloc(n * sizeof **mpZ);
	s = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    mpf_init_set_d (mpZ[0][s++], 1.0);
	}
	nvars++;
    }

    /* constant will have been moved to the end of the list by now */
    if (npoly && pmod->ifc) dropc = 1;

    /* process the ordinary data */
    for (i=0; i<l0-npoly-dropc; i++) {
	if (pmod->list[i+1] == 0) {
	    continue;
	} else {
	    mpZ[nvars] = malloc(n * sizeof **mpZ);
	    if (mpZ[nvars] == NULL) {
		return NULL;
	    }

	    /* record position in mpZ of the var to be raised 
	       to various powers, if applicable */
	    if (pmod->list[i+1] == pmod->polyvar) {
#ifdef MP_DEBUG
		printf("var to be raised to powers: it's "
		       "at position %d in the regression list,\n"
		       "and at slot %d in mpZ\n", i + 1, nvars);
#endif
		mp_poly_pos = nvars;
	    }

	    s = 0; /* obs. index for mpZ */
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		/* do trick with strings? */
		if (digits != NULL && digits[i+1] != NULL) {
		    sprintf(numstr, "%.*g", digits[i+1][t], 
			    Z[pmod->list[i+1]][t]);
#ifdef MP_DEBUG
		    printf("setting mpZ[%d][%d] from '%s'\n", 
			   nvars, s, numstr);
#endif
		    mpf_init_set_str (mpZ[nvars][s], numstr, 10);
		} else { /* or do straight conversion */
		    mpf_init_set_d (mpZ[nvars][s], 
				    Z[pmod->list[i+1]][t]);
		}
		s++;
	    }
	    pmod->varlist[i+1] = pmod->list[i+1];
	    pmod->list[i+1] = nvars++;
	}
    } /* end processing ordinary data */

    listpt = i + 1;

    /* generate polynomial data (if applicable) */
    for (i=0; i<npoly; i++) {  
	mpZ[nvars] = malloc(n * sizeof **mpZ);
	if (mpZ[nvars] == NULL) {
	    return NULL;
	}
	s = 0;
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    unsigned long pwr = pmod->polylist[i+1];

#ifdef MP_DEBUG
	    printf("generating mpZ[%d][%d] from mpZ[%d][%d],\n"
		   "using power %lu taken from polylist[%d]\n", 
		   nvars, s, mp_poly_pos, s, pwr, i+1);
#endif
	    mpf_init(mpZ[nvars][s]);
	    mpf_pow_ui(mpZ[nvars][s],          /* target */
		       mpZ[mp_poly_pos][s],    /* source */ 
		       pwr);                   /* power */
	    s++;
	}
	pmod->varlist[i+listpt] = pmod->polyvar;
	pmod->list[i+listpt] = nvars++;
    }

    /* ensure that the constant is properly registered in varlist */
    for (i=0; i<=pmod->list[0]; i++) {
	if (pmod->list[i] == 0) pmod->varlist[i] = 0;
    }

    return mpZ;
}

static void free_mpZ (mpf_t **mpZ, int v, int n)
{
    int i, t;

    for (i=0; i<v; i++) {
	for (t=0; t<n; t++) mpf_clear (mpZ[i][t]);
	free(mpZ[i]);
    }
    free(mpZ);
}

static void mp_model_free (MPMODEL *pmod)
{
    int i, l0 = pmod->list[0];

    free (pmod->list);
    free (pmod->varlist);

    if (pmod->coeff != NULL) {
	for (i=0; i<=pmod->ncoeff; i++) 
	    mpf_clear (pmod->coeff[i]);
	free (pmod->coeff);
    }

    if (pmod->sderr != NULL) {    
	for (i=0; i<=pmod->ncoeff; i++) 
	    mpf_clear (pmod->sderr[i]);
	free (pmod->sderr);
    }

    if (pmod->xpx != NULL) {
	for (i=0; i<=(l0-1)*l0/2; i++) 
	    mpf_clear (pmod->xpx[i]);
	free (pmod->xpx);
    }

    mpf_clear (pmod->ess);
    mpf_clear (pmod->tss);
    mpf_clear (pmod->sigma);
    mpf_clear (pmod->rsq);
    mpf_clear (pmod->adjrsq);
    mpf_clear (pmod->fstt); 
}

static void mp_model_init (MPMODEL *pmod, DATAINFO *pdinfo)
{
    pmod->ID = 0;
    pmod->t1 = pdinfo->t1;
    pmod->t2 = pdinfo->t2;
    pmod->list = NULL;
    pmod->varlist = NULL;
    pmod->polylist = NULL; /* don't free, the caller does that */
    pmod->ifc = 1;
    pmod->coeff = NULL;
    pmod->sderr = NULL;
    pmod->xpx = NULL;
    mpf_init (pmod->ess);
    mpf_init (pmod->tss);
    mpf_init (pmod->sigma);
    mpf_init (pmod->rsq);
    mpf_init (pmod->adjrsq);
    mpf_init (pmod->fstt); 
    pmod->errcode = 0;
    pmod->polyvar = 0;
}

static void get_mp_varname (const MPMODEL *pmod, const DATAINFO *pdinfo,
			    int c, char *vname)
{
    int realv = (pmod->polyvar)? 
	pmod->list[0] - pmod->polylist[0] : 0;

    if (pmod->polyvar && c >= realv 
	&& !(pmod->ifc && c == pmod->list[0])
	&& !(pmod->ifc == 0 && c == realv)) {
	int pwrpos = c - realv + pmod->ifc;

	sprintf(vname, "%s^%d", pdinfo->varname[pmod->polyvar],
		pmod->polylist[pwrpos]);
    } else {
	strcpy(vname, pdinfo->varname[pmod->varlist[c]]);
    } 
}

static void print_mp_coeff (const MPMODEL *pmod, const DATAINFO *pdinfo,
			    int c, PRN *prn)
{
    char vname[12];
    double xx = mpf_get_d (pmod->coeff[c-1]);
    double yy = mpf_get_d (pmod->sderr[c-1]);

    /* a bit of a fiddle getting the variable names to come
       out right */

    get_mp_varname (pmod, pdinfo, c, vname);

    pprintf(prn, " %3d) %8s ", pmod->varlist[c], vname);

    gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    gretl_print_fullwidth_double(yy, GRETL_MP_DIGITS, prn); 

    pputs(prn, "\n");
}

static void other_stats (const MPMODEL *pmod, PRN *prn)
{
    double xx;
    char fstr[16];
    int len = 24;

    if (doing_nls()) len = 36;
    
    xx = mpf_get_d (pmod->sigma);
    if (xx < DBL_EPSILON) xx = 0.0;
    pprintf(prn, "%-*s", len, _("Standard error"));
    gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    pputs(prn, "\n");

    xx = mpf_get_d (pmod->ess);
    if (xx < DBL_EPSILON) xx = 0.0;
    pprintf(prn, "%-*s", len, _("Error Sum of Squares"));
    gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    pputs(prn, "\n");

    xx = mpf_get_d (pmod->rsq);
    pprintf(prn, "%-*s", len, _("Unadjusted R-squared"));
    gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    pputs(prn, "\n");

    xx = mpf_get_d (pmod->adjrsq);
    pprintf(prn, "%-*s", len, _("Adjusted R-squared"));
    gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    pputs(prn, "\n");

    xx = mpf_get_d (pmod->fstt);
    sprintf(fstr, "F(%d, %d)", pmod->dfn, pmod->dfd);
    pprintf(prn, "%-*s", len, fstr);
    if (na(xx)) {
	pprintf(prn, "            %s", _("undefined"));
    } else {
	gretl_print_fullwidth_double(xx, GRETL_MP_DIGITS, prn);
    }
    pputs(prn, "\n");
}

static int print_mp_ols (const MPMODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    int i, ncoeff;
    char startdate[9], enddate[9];
    int t1 = pmod->t1, t2 = pmod->t2;

    ncoeff = pmod->list[0];
    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, _("Multiple-precision OLS estimates using "
		   "the %d observations %s-%s\n"),
	    pmod->nobs, startdate, enddate);
    pprintf(prn, "%s: %s\n\n", _("Dependent variable"),
	    pdinfo->varname[pmod->varlist[1]]);

    pputs(prn, _("      VARIABLE         COEFFICIENT          "
		 "        STD. ERROR\n"));

    if (pmod->ifc) {
	print_mp_coeff(pmod, pdinfo, ncoeff, prn);
	ncoeff--;
    }
    for (i=2; i<=ncoeff; i++) {
	print_mp_coeff(pmod, pdinfo, i, prn);
    }
    pputs(prn, "\n");

    other_stats (pmod, prn);

    return 0;
}

/**
 * mp_vector_raise_to_power:
 * @srcvec: source vector (doubles)
 * @targvec: vector to be filled in with results
 * @n: length of vector
 * @power: integer power to which elements of @srcvec should
 * be raised, using multiple precision arithmetic.
 *
 * Returns: 0 on success, error code on failure.
 */

int mp_vector_raise_to_power (const double *srcvec, double *targvec,
			      int n, unsigned pwr)
{
    int t;
    mpf_t src, targ;

    set_gretl_mp_bits();

    mpf_init (src);
    mpf_init (targ);

    for (t=0; t<n; t++) {
	if (na(srcvec[t])) {
	    targvec[t] = NADBL;
	    continue;
	}
	mpf_set_d (src, srcvec[t]);
	mpf_pow_ui(targ, src, (unsigned long) pwr);
	targvec[t] = mpf_get_d (targ);
    }

    mpf_clear (src);
    mpf_clear (targ);

    return 0;
}

#ifdef HAVE_MPFR
/**
 * mp_vector_ln:
 * @srcvec: source vector (doubles)
 * @targvec: vector to be filled in with results
 * @n: length of vector
 *
 * Returns: 0 on success, error code on failure.
 */

int mp_vector_ln (const double *srcvec, double *targvec, int n)
{
    int t;
    mpfr_t src, targ;

    set_gretl_mpfr_bits();

    mpfr_init (src);
    mpfr_init (targ);

    for (t=0; t<n; t++) {
	if (na(srcvec[t])) {
	    targvec[t] = NADBL;
	    continue;
	}
	mpfr_set_d (src, srcvec[t], GMP_RNDN);
	mpfr_log(targ, src, GMP_RNDN);
	targvec[t] = mpfr_get_d (targ, GMP_RNDN);
    }

    mpfr_clear (src);
    mpfr_clear (targ);

    return 0;
}
#endif

static int poly_check (MPMODEL *pmod, const int *list)
{
    int i;

    /* check that all powers are > 1 */

    for (i=1; i<=pmod->polylist[0]; i++) {
	if (pmod->polylist[i] < 2) return 1;
    }

    /* take the rightmost var in the regression list (other than 
       the constant) as the one to be raised to various powers */

    for (i=list[0]; i>1; i--) {
	if (list[i] != 0) {
	    pmod->polyvar = list[i];
	    break;
	}
    }

    if (pmod->polyvar == 0) return 1;
    return 0;
}

static int poly_copy_list (int **targ, const int *list, const int *poly)
{
    int i;

    *targ = malloc((list[0] + poly[0] + 1) * sizeof **targ);
    
    if (*targ == NULL) return 1;

    (*targ)[0] = list[0] + poly[0];

    for (i=1; i<=list[0]; i++) {
	(*targ)[i] = list[i];
    }

    for (i=1; i<=poly[0]; i++) {
	(*targ)[list[0] + i] = list[0] + i - 1; 
    }  

    return 0;
}

static void set_gretl_mp_bits (void)
{
    char *user_bits = getenv("GRETL_MP_BITS");
    
    if (user_bits != NULL) {
	mpf_set_default_prec ((unsigned long) atoi(user_bits));
    } else {
	mpf_set_default_prec ((unsigned long) DEFAULT_GRETL_MP_BITS);
    }
}

#ifdef HAVE_MPFR
static void set_gretl_mpfr_bits (void)
{
    char *user_bits = getenv("GRETL_MP_BITS");
    
    if (user_bits != NULL) {
	mpfr_set_default_prec ((unsigned long) atoi(user_bits));
    } else {
	mpfr_set_default_prec ((unsigned long) DEFAULT_GRETL_MP_BITS);
    }
} 
#endif

static int copy_mp_results (MPMODEL *pmod, DATAINFO *pdinfo,
			    mp_results *results)
{
    int i, err = 0;

    for (i=0; i<pmod->ncoeff; i++) {
	results->coeff[i] = mpf_get_d (pmod->coeff[i]);
	results->sderr[i] = mpf_get_d (pmod->sderr[i]);
    }

    if (pmod->ifc || (!pmod->ifc && pmod->list[0] == 2)) {
	results->coeff[0] = mpf_get_d (pmod->coeff[pmod->ncoeff]);
	results->sderr[0] = mpf_get_d (pmod->sderr[pmod->ncoeff]);
    } 

    results->sigma = mpf_get_d (pmod->sigma);
    results->ess = mpf_get_d (pmod->ess);
    results->rsq = mpf_get_d (pmod->rsq);
    results->fstt = mpf_get_d (pmod->fstt);

    if (results->varnames != NULL) { /* will use results for printing */
	int ncoeff = pmod->list[0];

	if (pmod->ifc) {
	    get_mp_varname(pmod, pdinfo, ncoeff, results->varnames[ncoeff]);
	    ncoeff--;
	}
	for (i=2; i<=ncoeff; i++) {
	    get_mp_varname(pmod, pdinfo, i, results->varnames[i]);
	}

	strcpy(results->varnames[0], pdinfo->varname[pmod->varlist[1]]);

	results->t1 = pmod->t1;
	results->t2 = pmod->t2;
	results->ifc = pmod->ifc;
	results->dfn = pmod->dfn;
	results->dfd = pmod->dfd;
	results->adjrsq = mpf_get_d (pmod->adjrsq);
	if (copylist(&results->varlist, pmod->varlist)) {
	    err = 1;
	}
    }

    return err;
}

/**
 * mplsq:
 * @list: dependent variable plus list of regressors.
 * @polylist: list of polynomial terms (or NULL).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @prn: where to print results.
 * @errbuf: where to print any error message.
 * @results: pointer to struct to hold results (or NULL).
 *
 * Computes multiple-precision OLS estimates of the model 
 * specified by @list.  If @results is NULL, prints the estimates
 * to @prn, otherwise puts them into @results.
 * 
 * Returns: 0 on success, error code on failure.
 */

int mplsq (const int *list, const int *polylist,
	   double ***pZ, DATAINFO *pdinfo, PRN *prn,
	   char *errbuf, mp_results *results) 
{
    int l0, i, ret = 0;
    mpf_t **mpZ = NULL;
    MPXPXXPY xpxxpy;
    MPMODEL model;

    *errbuf = 0;

    if (list == NULL || pZ == NULL || *pZ == NULL || pdinfo == NULL ||
	list[0] == 1 || pdinfo->v == 1) return E_DATA;

    set_gretl_mp_bits();

    mp_model_init (&model, pdinfo);

    /* preserve a copy of the list supplied, for future reference */
    if (polylist == NULL) {
	copylist(&(model.list), list);
    } else {
	poly_copy_list(&(model.list), list, polylist);
    }

    if (model.list == NULL) return E_ALLOC;

    model.polylist = polylist; /* attached for convenience */

    if (polylist != NULL && poly_check(&model, list)) {
	mp_model_free(&model);
	return E_DATA;
    }

    /* check for missing obs in sample */
    if (data_problems(list, *pZ, pdinfo, errbuf)) {
	mp_model_free(&model);
	return E_DATA;
    }

    /* see if the regressor list contains a constant */
    model.ifc = mp_rearrange(model.list);

    /* construct multiple-precision data matrix */
    mpZ = make_mpZ(&model, *pZ, pdinfo);

    if (mpZ == NULL) {
	mp_model_free(&model);
	return E_ALLOC;
    }

    mpf_constants_init();

    l0 = model.list[0];  
    model.ncoeff = l0 - 1; 
    model.nobs = model.t2 - model.t1 + 1;

    /* check degrees of freedom */
    if (model.nobs < model.ncoeff) { 
        sprintf(errbuf, _("No. of obs (%d) is less than no. "
			  "of parameters (%d)"), model.nobs, model.ncoeff);
	mp_model_free(&model);
	free_mpZ(mpZ, l0, model.nobs);
	mpf_constants_clear();
        return E_DF; 
    }

    /* calculate regression results */
    xpxxpy = mp_xpxxpy_func(model.list, model.nobs, mpZ);
    mpf_set (model.tss, xpxxpy.xpy[l0]);

    mp_regress(&model, xpxxpy, mpZ, model.nobs, errbuf);

    for (i=0; i<=l0; i++) mpf_clear (xpxxpy.xpy[i]);
    free(xpxxpy.xpy);
    xpxxpy.xpy = NULL;

    ret = model.errcode;

    if (ret == 0) {
	if (results == NULL) {
	    print_mp_ols (&model, pdinfo, prn);
	} else {
	    copy_mp_results (&model, pdinfo, results);
	}
    }

    /* free all the mpf stuff */
    free_mpZ(mpZ, l0, model.nobs);
    mp_model_free(&model);
    mpf_constants_clear();

    return ret;
}

/* .......................................................... */

static MPXPXXPY mp_xpxxpy_func (const int *list, int n, mpf_t **mpZ)
{
    int i, j, li, lj, m, l0 = list[0], yno = list[1], t;
    mpf_t xx, yy, z1, z2, tmp;
    MPXPXXPY xpxxpy;

    i = l0 - 1;
    m = i * (i + 1) / 2;

    if ((xpxxpy.xpy = malloc((l0 + 1) * sizeof(mpf_t))) == NULL ||
	(xpxxpy.xpx = malloc((m + 1) * sizeof(mpf_t))) == NULL) {
        xpxxpy.errcode = E_ALLOC;
        return xpxxpy;
    }

    for (i=0; i<=l0; i++) mpf_init(xpxxpy.xpy[i]);
    for (i=0; i<=m; i++) mpf_init(xpxxpy.xpx[i]);

    mpf_init(xx);
    mpf_init(yy);
    mpf_init(z1);
    mpf_init(z2);
    mpf_init(tmp);

    xpxxpy.nv = l0 - 1;

    for (t=0; t<n; t++) {
        mpf_set (xx, mpZ[yno][t]);
	mpf_add (xpxxpy.xpy[0], xpxxpy.xpy[0], xx);
	mpf_mul (yy, xx, xx);
	mpf_add (xpxxpy.xpy[l0], xpxxpy.xpy[l0], yy);
    }
    if (mpf_sgn (xpxxpy.xpy[l0]) == 0) {
         xpxxpy.ivalue = yno; 
         return xpxxpy; 
    }    
    m = 0;

    for (i=2; i<=l0; i++) {
        li = list[i];
        for (j=i; j<=l0; j++) {
            lj = list[j];
            mpf_set (xx, MPF_ZERO);
            for (t=0; t<n; t++) {
		mpf_mul (tmp, mpZ[li][t], mpZ[lj][t]);
		mpf_add (xx, xx, tmp);
	    }
            if ((mpf_sgn (xx) == 0) && li == lj)  {
                xpxxpy.ivalue = li;
                return xpxxpy;  
            }
            mpf_set (xpxxpy.xpx[++m], xx);
        }
        mpf_set (xx, MPF_ZERO);
        for (t=0; t<n; t++) {
	    mpf_mul (tmp, mpZ[yno][t], mpZ[li][t]);
	    mpf_add (xx, xx, tmp);
	}
        mpf_set (xpxxpy.xpy[i-1], xx);
    }
    xpxxpy.ivalue = 0;

    mpf_clear(xx);
    mpf_clear(yy);
    mpf_clear(z1);
    mpf_clear(z2);
    mpf_clear(tmp);

    return xpxxpy; 
}

/* .......................................................... */

static void mp_regress (MPMODEL *pmod, MPXPXXPY xpxxpy, mpf_t **mpZ, int n,
			char *errbuf)
{
    int i, v, nobs, nv, yno;
    mpf_t *diag, ysum, ypy, zz, rss, tss;
    mpf_t den, sgmasq, tmp;
    double ess;
    MPCHOLBETA cb;

    nv = xpxxpy.nv;
    yno = pmod->list[1];

    if ((pmod->sderr = malloc((nv + 1) * sizeof(mpf_t))) == NULL) {
        pmod->errcode = E_ALLOC;
        return;
    }

    for (i=0; i<nv+1; i++) mpf_init (pmod->sderr[i]);

    mpf_init (den);
    mpf_init (sgmasq);
    mpf_init (ysum);
    mpf_init (ypy);
    mpf_init (zz);
    mpf_init (rss);
    mpf_init (tss);
    mpf_init (tmp);

    nobs = pmod->nobs;
    pmod->ncoeff = nv;
    pmod->dfd = nobs - nv;

    if (pmod->dfd < 0) { 
       pmod->errcode = E_DF; 
       return; 
    }

    pmod->dfn = nv - pmod->ifc;
    mpf_set (ysum, xpxxpy.xpy[0]);
    mpf_set (ypy, xpxxpy.xpy[nv + 1]);
    if (mpf_sgn(ypy) == 0) { 
        pmod->errcode = E_YPY;
        return; 
    }

    mpf_mul (zz, ysum, ysum);
    mpf_set_d (tmp, (double) nobs);
    mpf_div (zz, zz, tmp);
    mpf_sub (tss, ypy, zz);
    if (mpf_sgn(tss) < 0) { 
        pmod->errcode = E_TSS; 
        return; 
    }

    /* Choleski-decompose X'X and find the coefficients */
    cb = mp_cholbeta(xpxxpy);
    pmod->coeff = cb.coeff;
    pmod->xpx = cb.xpxxpy.xpx;

    if (cb.errcode) {
        pmod->errcode = E_ALLOC;
        return;
    } 

    mpf_set (rss, cb.rss);
    mpf_clear (cb.rss);
    if (mpf_cmp(rss, MPF_MINUS_ONE) == 0) { 
        pmod->errcode = E_SINGULAR;
        return; 
    }

    mpf_sub (pmod->ess, ypy, rss);
    ess = mpf_get_d(pmod->ess);
    if (fabs(ess) < DBL_EPSILON) {
	mpf_set (pmod->ess, MPF_ZERO);
    }
    if (mpf_sgn(pmod->ess) < 0) { 
	sprintf(errbuf, _("Error sum of squares is not >= 0"));
        return; 
    }

    if (pmod->dfd == 0) {
	mpf_set (pmod->sigma, MPF_ZERO);
	mpf_set_d (pmod->adjrsq, NADBL);
    } else {
	mpf_set_d (tmp, (double) pmod->dfd);
	mpf_div (sgmasq, pmod->ess, tmp);
	mpf_sqrt (pmod->sigma, sgmasq);
	mpf_mul (den, tss, tmp);
    }

    if (mpf_sgn (tss) <= 0) {
	mpf_set_d (pmod->rsq, NADBL);
	mpf_set_d (pmod->adjrsq, NADBL);
	pmod->errcode = E_TSS;
	return;
    }       

    if (pmod->errcode) {
	fprintf(stderr, "mp_ols: pmod->errcode = %d\n", pmod->errcode);
	return;
    }

    mpf_div (tmp, pmod->ess, tss);
    mpf_sub (pmod->rsq, MPF_ONE, tmp);

    if (pmod->dfd > 0) {
	mpf_set_d (tmp, (double) (nobs - 1));
	mpf_div (tmp, tmp, den);
	mpf_mul (tmp, tmp, pmod->ess);
	mpf_sub (pmod->adjrsq, MPF_ONE, tmp);
	if (!pmod->ifc) { 
	    mpf_t df;

	    mpf_div (tmp, pmod->ess, ypy);
	    mpf_sub (pmod->rsq, MPF_ONE, tmp);
	    mpf_sub (tmp, MPF_ONE, pmod->rsq);
	    mpf_init_set_d (df, (double) (nobs - 1));
	    mpf_mul (tmp, tmp, df);
	    mpf_set_d (df, (double) pmod->dfd);
	    mpf_div (tmp, tmp, df);
	    mpf_sub (pmod->adjrsq, MPF_ONE, tmp);
	    mpf_clear (df);
	}
    }

    if (pmod->ifc && nv == 1) {
        mpf_set (zz, MPF_ZERO);
        pmod->dfn = 1;
    }

    if (mpf_sgn(sgmasq) != 1 || pmod->dfd == 0) {
	mpf_set_d (pmod->fstt, NADBL);
    } else { 
	mpf_set_d (tmp, (double) pmod->ifc);
	mpf_mul (tmp, zz, tmp);
	mpf_sub (pmod->fstt, rss, tmp);
	mpf_div (pmod->fstt, pmod->fstt, sgmasq);
	mpf_set_d (tmp, (double) pmod->dfn);
	mpf_div (pmod->fstt, pmod->fstt, tmp);
    }

    diag = malloc((nv + 1) * sizeof *diag); 
    if (diag == NULL) {
	pmod->errcode = E_ALLOC;
	return;
    }

    for (i=0; i<nv+1; i++) mpf_init (diag[i]);

    mp_diaginv(xpxxpy, diag);

    for (v=1; v<=nv; v++) { 
	mpf_sqrt (zz, diag[v]);
	mpf_mul (pmod->sderr[v], pmod->sigma, zz);
    }

    for (i=0; i<nv+1; i++) mpf_clear (diag[i]);
    free(diag); 

    mpf_clear (den);
    mpf_clear (sgmasq);
    mpf_clear (ysum);
    mpf_clear (ypy);
    mpf_clear (zz);
    mpf_clear (rss);
    mpf_clear (tss);
    mpf_clear (tmp);
    
    return;  
}

/* .......................................................... */

static MPCHOLBETA mp_cholbeta (MPXPXXPY xpxxpy)
{
    int nm1, i, j, k, kk, l, jm1, nv;
    mpf_t e, d, d1, test, xx, tmp;
    MPCHOLBETA cb;

    nv = xpxxpy.nv; 
    cb.errcode = 0; 
    mpf_init (cb.rss);

    if ((cb.coeff = malloc((nv + 1) * sizeof(mpf_t))) == NULL) {
        cb.errcode = E_ALLOC;
        return cb;
    }
    for (j=0; j<nv+1; j++) mpf_init (cb.coeff[j]);

    mpf_init (e);
    mpf_init (d);
    mpf_init (d1);
    mpf_init (test);
    mpf_init (xx);
    mpf_init (tmp);

    cb.xpxxpy = xpxxpy;

    nm1 = nv - 1;
    mpf_sqrt (tmp, xpxxpy.xpx[1]);
    mpf_div (e, MPF_ONE, tmp);
    mpf_set (xpxxpy.xpx[1], e);
    mpf_mul(xpxxpy.xpy[1], xpxxpy.xpy[1], e);
    for (i=2; i<=nv; i++) { 
	mpf_mul (xpxxpy.xpx[i], xpxxpy.xpx[i], e);
    }
    kk = nv + 1;

    for (j=2; j<=nv; j++) {
	/* diagonal elements */
	mpf_set (d, MPF_ZERO);
	mpf_set (d1,MPF_ZERO);
        k = j;
        jm1 = j - 1;
        for (l=1; l<=jm1; l++) {
	    mpf_set (xx, xpxxpy.xpx[k]);
	    mpf_mul (tmp, xx, xpxxpy.xpy[l]);
	    mpf_add (d1, d1, tmp);
	    mpf_mul (tmp, xx, xx);
	    mpf_add (d, d, tmp);
            k += nv-l;
        }
	mpf_sub (test, xpxxpy.xpx[kk], d);
        if (mpf_sgn(test) != 1) {
           mpf_set (cb.rss, MPF_MINUS_ONE); 
           return cb;
        }   
	mpf_sqrt (tmp, test);
	mpf_div (e, MPF_ONE, tmp);
	mpf_set (xpxxpy.xpx[kk], e);
	mpf_sub (tmp, xpxxpy.xpy[j], d1);
	mpf_mul (xpxxpy.xpy[j], tmp, e);

        /* off-diagonal elements */
        for (i=j+1; i<=nv; i++) {
            kk++;
            mpf_set (d, MPF_ZERO);
            k = j;
            for (l=1; l<=jm1; l++) {
		mpf_mul (tmp, xpxxpy.xpx[k], xpxxpy.xpx[k-j+i]);
		mpf_add (d, d, tmp);
                k += nv - l;
            }
	    mpf_sub (tmp, xpxxpy.xpx[kk], d);
	    mpf_mul (xpxxpy.xpx[kk], tmp, e);
        }
        kk++;
    }
    kk--;
    mpf_set (d, MPF_ZERO);
    for(j=1; j<=nv; j++) {
	mpf_mul (tmp, xpxxpy.xpy[j], xpxxpy.xpy[j]);
	mpf_add (d, d, tmp);
    }
    mpf_set (cb.rss, d);
    mpf_mul (cb.coeff[nv], xpxxpy.xpy[nv], xpxxpy.xpx[kk]);
    for(j=nm1; j>=1; j--) {
	mpf_set (d, xpxxpy.xpy[j]);
        for (i=nv; i>=j+1; i--) {
            kk--;
	    mpf_mul (tmp, cb.coeff[i], xpxxpy.xpx[kk]);
	    mpf_sub (d, d, tmp);
        }
        kk--;
	mpf_mul (cb.coeff[j], d, xpxxpy.xpx[kk]);
    }   

    mpf_clear (e);
    mpf_clear (d);
    mpf_clear (d1);
    mpf_clear (test);
    mpf_clear (xx);
    mpf_clear (tmp);
 
    return cb; 
}

/* ...............................................................    */

static void mp_diaginv (MPXPXXPY xpxxpy, mpf_t *diag)
{
    int kk = 1, l, m, nstop, k, i, j, nv;
    mpf_t d, e, tmp;

    mpf_init (d);
    mpf_init (e);
    mpf_init (tmp);

    nv = xpxxpy.nv;
    nstop = nv * (nv+1)/2;
    for (l=1; l<=nv-1; l++) {
	mpf_set (d, xpxxpy.xpx[kk]);
	mpf_set (xpxxpy.xpy[l], d);
	mpf_mul (e, d, d);
        m = 0;
        if (l > 1) 
	    for (j=1; j<=l-1; j++) m += nv - j;
        for (i=l+1; i<=nv; i++) {
	    mpf_set (d, MPF_ZERO); 
            k = i + m;
            for (j=l; j<=i-1; j++) {
		mpf_mul (tmp, xpxxpy.xpy[j], xpxxpy.xpx[k]);
		mpf_add (d, d, tmp);
                k += nv - j;
            }
	    mpf_mul (d, d, xpxxpy.xpx[k]);
	    mpf_mul (d, d, MPF_MINUS_ONE);
	    mpf_set (xpxxpy.xpy[i], d);
	    mpf_mul (tmp, d, d);
	    mpf_add (e, e, tmp);
        }
        kk += nv + 1 - l;
        mpf_set (diag[l], e);
    }

    mpf_mul (diag[nv], xpxxpy.xpx[nstop], xpxxpy.xpx[nstop]);

    mpf_clear (d);
    mpf_clear (e);
    mpf_clear (tmp);
}

/* .........................................................   */

static int mp_rearrange (int *list)
/* checks a list for a constant term (ID # 0), and if present, 
   move it to the last position.  Return 1 if constant found,
   else 0.
*/
{
    int lo = list[0], v;

    for (v=2; v<=lo; v++) {
        if (list[v] == 0)  {
            list_exclude(v, list);
            list[0] = lo;
            list[lo] = 0;
            return 1;
        }
    }
    return 0;
}
