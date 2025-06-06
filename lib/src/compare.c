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

#include "libgretl.h"
#include "libset.h"
#include "gretl_drivers.h"
#include "gretl_panel.h"
#include "var.h"
#include "system.h"
#include "missing_private.h"
#include "matrix_extra.h"
#include "plotspec.h"
#include "tsls.h"
#include "uservar.h"

/**
 * SECTION:compare
 * @short_description: diagnostic and specification tests for models
 * @title: Model tests
 * @include: libgretl.h
 *
 * Included here are several tests for "pathologies" of the error
 * term in regression models, as well as specification tests
 * covering nonlinearity and the omission or addition of variables.
 */

#define WDEBUG 0

struct COMPARE {
    int ci;        /* command index: ADD or OMIT */
    int model_id;  /* ID of model under test */
    int model_ci;  /* estimator code for the model under test */
    int dfn;       /* numerator degrees of freedom */
    int dfd;       /* denominator degrees of freedom */
    double test;   /* test statistic */
    double pval;   /* p-value for test statistic */
    int stat;      /* code identifying the test statistic */
    double F;      /* F test statistic */
    double X2;     /* Wald chi-square test statistic */
    double LR;     /* likelihood ratio test statistic */
    int score;     /* number of info stats showing improvement */
    int robust;    /* = 1 when robust vcv is in use, else 0 */
    const int *testvars; /* list of variables added or omitted */
};

/* Given a list of variables, check them against the independent
   variables included in a model, and construct a mask with 1s in
   positions where there is a match, 0s otherwise.  If the test list
   is NULL, match all variables except the constant.
*/

static char *
mask_from_test_list (const int *list, const MODEL *pmod, int *err)
{
    char *mask;
    int off1 = 2, off2 = 0;
    int cmax = pmod->ncoeff;
    int nmask = 0;
    int i, j;

    mask = calloc(pmod->ncoeff, 1);
    if (mask == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (pmod->ci == DPANEL) {
	/* find correct offset into independent vars in list */
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == LISTSEP) {
		off1 = i + 2;
	    }
	}
	off2 = pmod->list[1];
    } else if (pmod->ci == NEGBIN) {
	cmax--;
    }

    for (i=0; i<cmax; i++) {
	if (list != NULL) {
	    for (j=1; j<=list[0]; j++) {
		if (pmod->list[i + off1] == list[j]) {
#if WDEBUG
		    fprintf(stderr, "matched var %d at pmod->list[%d]: "
			    "set mask[%d] = 1\n", list[j], i + off1, i + off2);
		    printlist(list, "test list");
		    printlist(pmod->list, "pmod->list");
#endif
		    mask[i + off2] = 1;
		    nmask++;
		}
	    }
	} else if (pmod->list[i + off1] != 0) {
	    mask[i + off2] = 1;
	}
    }

    if (list != NULL && nmask != list[0]) {
	fprintf(stderr, "mask from list: list[0] = %d but nmask = %d\n",
		list[0], nmask);
	*err = E_DATA;
    }

#if WDEBUG
    for (i=0; i<pmod->ncoeff; i++) {
	fprintf(stderr, "mask[%d] = %d\n", i, mask[i]);
    }
#endif

    return mask;
}

/* Wald (chi-square and/or F) test for a set of zero restrictions on
   the parameters of a given model, based on the covariance matrix of
   the unrestricted model. Suitable for use where the original model
   is estimated by FGLS or IV.  Note that if list is NULL, we do an
   automatic test, for the significance of all vars but the constant.
*/

static int
wald_test (const int *list, MODEL *pmod, double *chisq, double *F)
{
    char *mask = NULL;
    gretl_matrix *C = NULL;
    gretl_vector *b = NULL;
    double wX = NADBL;
    double wF = NADBL;
    int err = 0;

    mask = mask_from_test_list(list, pmod, &err);
    if (err) {
	free(mask);
	return err;
    }

    if (!err) {
	C = gretl_vcv_matrix_from_model(pmod, mask, &err);
    }

    if (!err) {
	b = gretl_coeff_vector_from_model(pmod, mask, &err);
    }

#if WDEBUG
    gretl_matrix_print(C, "VCV");
    gretl_matrix_print(b, "coeff");
#endif

    if (!err) {
	/* use "left division" to explicit inversion of @C */
	gretl_matrix *tmp, *mx;

	mx = gretl_matrix_alloc(1, 1);
	tmp = gretl_matrix_divide(C, b, GRETL_MOD_NONE, &err);
	if (tmp != NULL) {
	    gretl_matrix_multiply_mod(b, GRETL_MOD_TRANSPOSE,
				      tmp, GRETL_MOD_NONE,
				      mx, GRETL_MOD_NONE);
	    wX = mx->val[0];
	    gretl_matrix_free(tmp);
	    gretl_matrix_free(mx);
	}
    }

#if WDEBUG
    fprintf(stderr, "wX (quadratic form) = %g\n", wX);
#endif

    if (!err) {
	if (wX < 0) {
	    wF = wX = NADBL;
	} else {
	    wF = wX / gretl_vector_get_length(b);
	}
    }

    if (!err) {
	if (chisq != NULL) {
	    *chisq = wX;
	}
	if (F != NULL) {
	    *F = wF;
	}
    }

#if WDEBUG
    fprintf(stderr, "Wald test: F = %g, Chi^2 = %g\n", wF, wX);
#endif

    free(mask);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    return err;
}

/**
 * wald_omit_F:
 * @list: list of variables to omit, or NULL.
 * @pmod: model to be tested.
 *
 * Simple form of Wald F-test for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated F-value, or #NADBL on failure.
 */

double wald_omit_F (const int *list, MODEL *pmod)
{
    double F = NADBL;

    wald_test(list, pmod, NULL, &F);
    return F;
}

/**
 * wald_omit_chisq:
 * @list: list of variables to omit, or NULL.
 * @pmod: model to be tested.
 *
 * Simple form of Wald chi-square for omission of variables.  If @list
 * is non-NULL, do the test for the omission of the variables in
 * @list from the model @pmod.  Otherwise test for omission of
 * all variables in @pmod except for the constant.
 *
 * Returns: Calculated chi-square value, or #NADBL on failure.
 */

double wald_omit_chisq (const int *list, MODEL *pmod)
{
    double X = NADBL;

    wald_test(list, pmod, &X, NULL);
    return X;
}

static int add_diffvars_to_test (ModelTest *test, const int *list,
				 const DATASET *dset)
{
    char *vnames;
    int i, len = 0;
    int err = 0;

    for (i=1; i<=list[0]; i++) {
	len += strlen(dset->varname[list[i]]) + 1;
    }

    vnames = malloc(len);

    if (vnames == NULL) {
	err = 1;
    } else {
	*vnames = '\0';
	for (i=1; i<=list[0]; i++) {
	    strcat(vnames, dset->varname[list[i]]);
	    if (i < list[0]) {
		strcat(vnames, " ");
	    }
	}
	model_test_set_allocated_param(test, vnames);
    }

    return err;
}

void print_add_omit_null (const int *list, const DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    int nl = list == NULL ? 0 : list[0];

    if (nl == 1 && !(opt & OPT_S)) {
	/* single parameter, not in system context */
	pputs(prn, "\n  ");
	pprintf(prn, _("Null hypothesis: the regression parameter is zero for %s"),
		dset->varname[list[1]]);
	pputc(prn, '\n');
    } else if (nl == 0) {
	/* VAR omit specials */
	if ((opt & OPT_E) && (opt & OPT_T)) {
	    pprintf(prn, "\n  %s: %s\n", _("Null hypothesis"),
		    _("no seasonal effects or trend"));
	} else if (opt & OPT_E) {
	    pprintf(prn, "\n  %s: %s\n", _("Null hypothesis"),
		    _("no seasonal effects"));
	} else if (opt & OPT_T) {
	    pprintf(prn, "\n  %s: %s\n", _("Null hypothesis"),
		    _("no trend"));
	}
    } else {
	/* all other cases */
	const char *vname;
	int i, nc = 0;

	pputs(prn, _("\n  Null hypothesis: the regression parameters are "
		     "zero for the variables\n"));

	pputs(prn, "    ");
	for (i=1; i<=list[0]; i++) {
	    vname = dset->varname[list[i]];
	    nc += strlen(vname) + 2;
	    pprintf(prn, "%s", vname);
	    if (i < list[0]) {
		if (nc > 60) {
		    pputs(prn, ",\n    ");
		    nc = 0;
		} else {
		    pputs(prn, ", ");
		}
	    }
	}
	pputc(prn, '\n');
	if (opt & OPT_E) {
	    /* seasonals */
	    pprintf(prn, "    %s\n", _("seasonal dummies"));
	}
	if (opt & OPT_T) {
	    /* trend */
	    pputs(prn, "    time\n");
	}
    }
}

/* add/omit: are we printing a revised model following the
   test results? */

static int printing_revised_model (gretlopt opt)
{
    if (opt & OPT_Q) {
	/* --quiet: no */
	return 0;
    } else if (opt & (OPT_L | OPT_W)) {
	/* --lm or --wald test variants: no */
	return 0;
    } else {
	/* then yes */
	return 1;
    }
}

/* print the results from an add or omit test */

static void print_compare (struct COMPARE *cmp,
			   const DATASET *dset,
			   gretlopt opt,
			   PRN *prn)
{
    if (opt & OPT_Q) {
	if (cmp->model_id > 0) {
	    pprintf(prn, _("Test on Model %d:"), cmp->model_id);
	    pputc(prn, '\n');
	}
    } else if (cmp->model_id > 0) {
	if (opt & OPT_A) {
	    /* --auto-omit: add vertical space */
	    pputc(prn, '\n');
	}
	pprintf(prn, _("Test on Model %d:"), cmp->model_id);
	pputc(prn, '\n');
    }

    print_add_omit_null(cmp->testvars, dset, OPT_NONE, prn);

    if (cmp->stat == GRETL_STAT_WALD_CHISQ) {
	pprintf(prn, "  %s: %s(%d) = %g, %s %g\n",  _("Wald test"),
		_("Chi-square"), cmp->dfn, cmp->test,
		_("p-value"), cmp->pval);
	if (!na(cmp->LR)) {
	    /* auxiliary LR test? */
	    double pval = chisq_cdf_comp(cmp->dfn, cmp->LR);

	    pprintf(prn, "  (%s: %s(%d) = %g, %s %g)\n",  _("LR test"),
		    _("Chi-square"), cmp->dfn, cmp->LR,
		    _("p-value"), pval);
	} else if (!na(cmp->F)) {
	    /* alternate F-form for Wald? */
	    double pval = snedecor_cdf_comp(cmp->dfn, cmp->dfd, cmp->F);

	    pprintf(prn, "  (%s: F(%d, %d) = %g, %s %g)\n",
		    _("F-form"), cmp->dfn, cmp->dfd, cmp->F,
		    _("p-value"), pval);
	}
    } else if (cmp->stat == GRETL_STAT_LR) {
	/* note: unused */
	pprintf(prn, "  %s: %s(%d) = %g, %s %g\n",  _("LR test"),
		_("Chi-square"), cmp->dfn, cmp->test,
		_("p-value"), cmp->pval);
    } else if (cmp->stat == GRETL_STAT_F) {
	pprintf(prn, "  %s: %s(%d, %d) = %g, %s %g\n", _("Test statistic"),
		(cmp->robust)? _("Robust F") : "F",
		cmp->dfn, cmp->dfd, cmp->test,
		_("p-value"), cmp->pval);
    } else if (cmp->stat == GRETL_STAT_LM) {
	pprintf(prn, "  %s: %s(%d) = %g, %s %g\n",  _("LM test"),
		_("Chi-square"), cmp->dfn, cmp->test,
		_("p-value"), cmp->pval);
    }

    if (cmp->score >= 0) {
	pputs(prn, "  ");
	if (cmp->ci == ADD) {
	    pprintf(prn, _("Adding variables improved %d of %d information "
			   "criteria.\n"), cmp->score, C_MAX);
	} else {
	    pprintf(prn, _("Omitting variables improved %d of %d information "
			   "criteria.\n"), cmp->score, C_MAX);
	}
    }

    if (!printing_revised_model(opt)) {
	pputc(prn, '\n');
    }
}

/* try attaching the add or omit test to @pmod */

static void add_omit_attach (struct COMPARE *cmp, MODEL *pmod,
			     const DATASET *dset)
{
    int tcode = (cmp->ci == OMIT)? GRETL_TEST_OMIT : GRETL_TEST_ADD;
    ModelTest *mtest = model_test_new(tcode);

    if (mtest != NULL) {
	model_test_set_teststat(mtest, cmp->stat);
	model_test_set_dfn(mtest, cmp->dfn);
	if (cmp->stat == GRETL_STAT_F) {
	    model_test_set_dfd(mtest, cmp->dfd);
	}
	model_test_set_value(mtest, cmp->test);
	model_test_set_pvalue(mtest, cmp->pval);
	add_diffvars_to_test(mtest, cmp->testvars, dset);
	maybe_add_test_to_model(pmod, mtest);
    }
}

/* add a record of improvement/degeneration of the
   model selection criteria, for non-quiet printing:
   @pmodA is the original model, and @pmodB the
   restricted (OMIT) or augmented (ADD) variant
*/

static void maybe_add_info_score (struct COMPARE *cmp,
				  const MODEL *pmodA,
				  const MODEL *pmodB)
{
    if (pmodA != NULL && pmodB != NULL) {
	int i;

	cmp->score = 0;
	for (i=0; i<C_MAX; i++) {
	    if (na(pmodB->criterion[i]) || na(pmodA->criterion[i])) {
		cmp->score = -1;
		break;
	    } else if (pmodB->criterion[i] < pmodA->criterion[i]) {
		cmp->score++;
	    }
	}
    }
}

static void compare_init (struct COMPARE *cmp,
			  MODEL *pmod, int ci,
			  const int *testvars)
{
    cmp->ci = ci;
    cmp->stat = 0;
    cmp->test = cmp->pval = NADBL;
    cmp->F = cmp->X2 = cmp->LR = NADBL;
    cmp->score = -1;
    cmp->robust = 0;
    cmp->dfn = testvars[0];
    cmp->testvars = testvars;
    cmp->model_ci = pmod->ci;
    cmp->model_id = pmod->ID;
}

static int add_or_omit_compare (MODEL *pmodA, MODEL *pmodB,
				const int *testvars, const DATASET *dset,
				int ci, gretlopt opt, PRN *prn)
{
    struct COMPARE cmp;
    MODEL *umod, *rmod;
    int iv_special = 0;
    int err = 0;

    compare_init(&cmp, pmodA, ci, testvars);

    if (ci == OMIT) {
	/* the unrestricted model is the original one, 'A' */
	umod = pmodA;
	rmod = pmodB;
    } else {
	/* add: the unrestricted model is the new one, 'B' */
	umod = pmodB;
	rmod = pmodA;
    }

    if (opt & OPT_B) {
	/* ivreg special: A and B have different
	   instruments, not just different regressors
	*/
	iv_special = 1;
	cmp.model_id = -1;
    } else if (opt & OPT_L) {
	cmp.model_id = -1;
    } else if (umod != NULL && rmod != NULL) {
	cmp.dfn = umod->ncoeff - rmod->ncoeff;
    }

    if (pmodA->opt & OPT_R) {
	cmp.robust = 1;
    }

    if (pmodA->ci == DPANEL) {
	/* FIXME: plus some other cases? */
	opt |= OPT_X;
    }

    cmp.dfd = umod->dfd;

    if (opt & OPT_L) {
	/* "add" with LM option: the auxiliary regression
	   results are in pmodB */
	cmp.stat = GRETL_STAT_LM;
	cmp.test = pmodB->nobs * pmodB->rsq;
	cmp.pval = chisq_cdf_comp(cmp.dfn, cmp.test);
    } else if (cmp.model_ci == OLS && !cmp.robust &&
	       rmod != NULL && umod != NULL &&
	       rmod->nobs == umod->nobs &&
	       !(opt & OPT_X)) {
	/* plain OLS with both the restricted and unrestricted
	   models available: base F-test on sums of squared
	   residuals
	*/
	cmp.stat = GRETL_STAT_F;
	cmp.test = ((rmod->ess - umod->ess) / umod->ess) *
	    cmp.dfd / cmp.dfn;
	cmp.pval = snedecor_cdf_comp(cmp.dfn, cmp.dfd, cmp.test);
    } else if (opt & OPT_X) {
	/* chi-square form of Wald test is requested */
	cmp.stat = GRETL_STAT_WALD_CHISQ;
	err = wald_test(testvars, umod, &cmp.test, &cmp.F);
	if (!err) {
	    cmp.pval = chisq_cdf_comp(cmp.dfn, cmp.test);
	}
    } else {
	/* F-form of Wald test */
	cmp.stat = GRETL_STAT_F;
	err = wald_test(testvars, umod, &cmp.X2, &cmp.test);
	if (!err) {
	    cmp.pval = snedecor_cdf_comp(cmp.dfn, cmp.dfd, cmp.test);
	}
    }

    /* auxiliary LR test? */
    if (umod != NULL && rmod != NULL && !(opt & OPT_L)) {
	if (!na(umod->lnL) && !na(rmod->lnL)) {
	    cmp.LR = 2.0 * (umod->lnL - rmod->lnL);
	}
    }

    if (!err && (na(cmp.test) || na(cmp.pval))) {
	err = E_DATA;
    }

    if (!err) {
	record_test_result(cmp.test, cmp.pval);
	if (opt & OPT_S) {
	    /* attach test to model */
	    add_omit_attach(&cmp, pmodA, dset);
	}
	if (!(opt & OPT_I)) {
	    /* not --silent */
	    if (!(opt & OPT_Q) && cmp.stat != GRETL_STAT_LM &&
		!iv_special) {
		maybe_add_info_score(&cmp, pmodA, pmodB);
	    }
	    print_compare(&cmp, dset, opt, prn);
	}
    }

    return err;
}

static int get_extra_var (const MODEL *pmod)
{
    if (COUNT_MODEL(pmod->ci)) {
	return gretl_model_get_int(pmod, "offset_var");
    } else if (pmod->ci == DURATION) {
	return gretl_model_get_int(pmod, "cens_var");
    } else {
	return 0;
    }
}

/* reconstitute full varlist for WLS, AR, count and
   duation models */

static int *
full_model_list (const MODEL *pmod, const int *inlist)
{
    int *flist = NULL;

    if (pmod->ci == AR) {
	/* cobble together arlist and @inlist */
	flist = gretl_lists_join_with_separator(pmod->arinfo->arlist,
						inlist);
    } else {
	int i, len = inlist[0];

	if (pmod->ci == WLS) {
	    /* prepend the weight variable */
	    len += 1;
	} else if (COUNT_MODEL(pmod->ci) || pmod->ci == DURATION) {
	    /* append list separator and offset or censoring var */
	    len += 2;
	}

	flist = gretl_list_new(len);
	if (flist == NULL) {
	    return NULL;
	}

	if (pmod->ci == WLS) {
	    flist[1] = pmod->nwt;
	    for (i=1; i<=inlist[0]; i++) {
		flist[i+1] = inlist[i];
	    }
	} else if (COUNT_MODEL(pmod->ci) || pmod->ci == DURATION) {
	    int extra = get_extra_var(pmod);

	    for (i=1; i<=inlist[0]; i++) {
		flist[i] = inlist[i];
	    }
	    flist[flist[0]-1] = LISTSEP;
	    flist[flist[0]] = extra;
	}
    }

    return flist;
}

static gretlopt retrieve_dpanel_opts (const MODEL *pmod)
{
    gretlopt opt = OPT_NONE;

    if (pmod->opt & OPT_D) {
	/* --dpdstyle */
	opt |= OPT_D;
    }

    if (pmod->opt & OPT_L) {
	/* --system (include levels) */
	opt |= OPT_L;
    }

    if (gretl_model_get_int(pmod, "asy")) {
	opt |= OPT_A;
    }

    if (gretl_model_get_int(pmod, "step") == 2) {
	opt |=OPT_T;
    }

    return opt;
}

static int obs_diff_ok (const MODEL *m_old, const MODEL *m_new)
{
    int tdiff, ndiff = m_new->nobs - m_old->nobs;

    if (m_old->ci == AR1) {
	return 0;
    }

    if (ndiff > 0) {
	tdiff = (m_new->t2 - m_new->t1) - (m_old->t2 - m_old->t1);
	if (ndiff == tdiff) {
	    return 1;
	}
    }

    return 0;
}

/* This function is used for "add" and "omit", when we are estimating
   an augmented or reduced version of the original model. It's also
   used in the special case of calculation of the p-value for the
   Durbin-Watson statistic, which requires re-estimation of the
   original specification. In these cases we need to ensure
   comparability, which means we have to retrieve any relevant options
   from the original model and re-apply them.
*/

static MODEL replicate_estimator (const MODEL *orig, int *list,
				  DATASET *dset, gretlopt myopt,
				  PRN *prn)
{
    MODEL rep;
    const char *param = NULL;
    const int *laglist = NULL;
    const char *cname = NULL;
    int *full_list = NULL;
    char altparm[32] = {0};
    int mc = get_model_count();
    int repci = orig->ci;
    int cv, order = 0;
    int first = 1;

    gretl_model_init(&rep, dset);

    /* recreate options and auxiliary vars, if required */

    transcribe_option_flags(&myopt, orig->opt, OPT_D | OPT_J | OPT_R);

    cname = gretl_model_get_cluster_vname(orig);
    if (cname != NULL) {
	/* FIXME second cluster-var name? */
	cv = current_series_index(dset, cname);
	if (cv > 0) {
	    myopt |= OPT_C;
	    if (orig->ci == OLS && gretl_model_get_int(orig, "pooled")) {
		/* transcribe the cluster spec to "panel", since
		   we'll be calling panel_model() below
		*/
		set_optval_string(PANEL, OPT_C, cname);
	    } else {
		set_optval_string(orig->ci, OPT_C, cname);
	    }
	}
    }

    if (orig->ci == AR1) {
	if (orig->opt & OPT_H) {
	    myopt |= OPT_H;
	    if (gretl_model_get_int(orig, "no-corc")) {
		myopt |= OPT_B;
	    }
	} else if (orig->opt & OPT_P) {
	    myopt |= OPT_P;
	}
    } else if (orig->ci == WLS || orig->ci == AR || get_extra_var(orig)) {
	full_list = full_model_list(orig, list);
	if (full_list == NULL) {
	    rep.errcode = E_ALLOC;
	} else {
	    list = full_list;
	}
    } else if (orig->ci == HSK) {
	if (gretl_model_get_int(orig, "no-squares")) {
	    myopt |= OPT_N;
	}
    } else if (orig->ci == DPANEL) {
	param = gretl_model_get_data(orig, "istr");
	laglist = gretl_model_get_list(orig, "ylags");
	myopt |= retrieve_dpanel_opts(orig);
    } else if (orig->ci == ARCH) {
	order = gretl_model_get_int(orig, "arch_order");
    } else if (orig->ci == LOGIT || orig->ci == PROBIT) {
	if (orig->opt & OPT_P) {
	    /* p-values, not slopes, selected */
	    myopt |= OPT_P;
	} else if (gretl_model_get_int(orig, "ordered")) {
	    myopt |= OPT_D;
	} else if (gretl_model_get_int(orig, "multinom")) {
	    myopt |= OPT_M;
	} else if (orig->ci == PROBIT && (orig->opt & OPT_E)) {
	    /* random effects */
	    int qp = gretl_model_get_int(orig, "quadpoints");

	    myopt |= (OPT_E | OPT_G);
	    set_optval_double(PROBIT, OPT_G, qp);
	}
    } else if (orig->ci == PANEL) {
	if (gretl_model_get_int(orig, "pooled")) {
	    /* pooled OLS */
	    myopt |= OPT_P;
	} else if (orig->opt & OPT_U) {
	    /* random effects */
	    myopt |= OPT_U;
	} else if (orig->opt & OPT_H) {
	    /* unit weights */
	    myopt |= OPT_H;
	    if (gretl_model_get_int(orig, "iters")) {
		myopt |= OPT_I;
	    }
	}
    } else if (orig->ci == LAD && gretl_model_get_int(orig, "rq")) {
	double x;

	x = gretl_model_get_double(orig, "tau");
	sprintf(altparm, "%g", x);
	if (gretl_model_get_int(orig, "rq_nid")) {
	    myopt |= OPT_R;
	}
	x = gretl_model_get_double(orig, "rq_alpha");
	if (!na(x)) {
	    myopt |= OPT_I;
	    set_optval_double(QUANTREG, OPT_I, x);
	}
    } else if (orig->ci == TOBIT) {
	double x;

	x = gretl_model_get_double(orig, "llimit");
	if (!na(x)) {
	    myopt |= OPT_L;
	    set_optval_double(TOBIT, OPT_L, x);
	}
	x = gretl_model_get_double(orig, "rlimit");
	if (!na(x)) {
	    myopt |= OPT_M;
	    set_optval_double(TOBIT, OPT_M, x);
	}
    } else if (orig->ci == IVREG) {
	transcribe_option_flags(&myopt, orig->opt, OPT_L | OPT_G);
    }

    if (rep.errcode) {
	goto bailout;
    }

 try_again:

    switch (orig->ci) {

    case AR:
	rep = ar_model(list, dset, myopt, NULL);
	break;
    case AR1:
	rep = ar1_model(list, dset, myopt, NULL);
	break;
    case DPANEL:
	rep = dpd_model(list, laglist, param, dset, myopt, prn);
	break;
    case ARCH:
	rep = arch_model(list, order, dset, myopt);
	break;
    case LOGIT:
    case PROBIT:
	rep = logit_probit(list, dset, orig->ci, myopt, NULL);
	break;
    case TOBIT:
	rep = tobit_driver(list, dset, myopt, NULL);
	break;
    case LAD:
	if (gretl_model_get_int(orig, "rq")) {
	    rep = quantreg_driver(altparm, list, dset, myopt, NULL);
	} else {
	    rep = lad_model(list, dset, myopt);
	}
	break;
    case POISSON:
    case NEGBIN:
	rep = count_model(list, orig->ci, dset, myopt, NULL);
	break;
    case DURATION:
	rep = duration_model(list, dset, myopt, NULL);
	break;
    case HECKIT:
	rep = heckit_model(list, dset, myopt, NULL);
	break;
    case IVREG:
	rep = ivreg(list, dset, myopt);
	break;
    case LOGISTIC:
	{
	    double lmax = gretl_model_get_double(orig, "lmax");

	    rep = logistic_model(list, lmax, dset, myopt);
	}
	break;
    case PANEL:
	rep = panel_model(list, dset, myopt, prn);
	break;
    case HSK:
	rep = hsk_model(list, dset, myopt);
	break;
    default:
	/* handles OLS, WLS, etc. */
	if (gretl_model_get_int(orig, "pooled")) {
	    myopt |= OPT_P;
	    rep = panel_model(list, dset, myopt, prn);
	} else {
	    rep = lsq(list, dset, repci, myopt);
	}
	break;
    }

#if 0
    fprintf(stderr, "replicate_estimator:\n"
	    " orig: t1=%d, t2=%d, nobs = %d\n"
	    " rep:  t1=%d, t2=%d, nobs = %d\n",
	    orig->t1, orig->t2, orig->nobs,
	    rep.t1, rep.t2, rep.nobs);
#endif

    /* check that we got the same sample as the original */
    if (!rep.errcode && rep.nobs != orig->nobs) {
	if (first && obs_diff_ok(orig, &rep)) {
	    dset->t1 = orig->t1;
	    dset->t2 = orig->t2;
	    clear_model(&rep);
	    first = 0;
	    goto try_again;
	} else {
	    fprintf(stderr, "Original obs = %d but new = %d\n", orig->nobs, rep.nobs);
	    rep.errcode = E_DATA;
	}
    }

 bailout:

    free(full_list);

    /* if the model count went up for an aux regression,
       bring it back down */
    if ((myopt & OPT_A) && get_model_count() > mc) {
	model_count_minus(NULL);
    }

#if 0
    printmodel(&rep, dset, OPT_NONE, prn);
#endif

    return rep;
}

static int add_residual_to_dataset (MODEL *pmod, DATASET *dset)
{
    int err = 0;

    if (dataset_add_series(dset, 1)) {
	err = E_ALLOC;
    } else {
	int t, v = dset->v - 1;

	for (t=0; t<dset->n; t++) {
	    dset->Z[v][t] = pmod->uhat[t];
	}

	strcpy(dset->varname[v], "uhat");
	series_set_label(dset, v, _("residual"));
    }

    return err;
}

static void nonlin_test_header (int code, PRN *prn)
{
    pputc(prn, '\n');
    if (code == AUX_SQ) {
	pputs(prn, _("Non-linearity test (squared terms)"));
    } else {
	pputs(prn, _("Non-linearity test (log terms)"));
    }
    pputs(prn, "\n\n");
}

static int
real_nonlinearity_test (MODEL *pmod, int *list,
			DATASET *dset, int aux_code,
			gretlopt opt, PRN *prn)
{
    MODEL aux;
    int err;

    err = add_residual_to_dataset(pmod, dset);
    if (err) {
	return err;
    }

    /* replace the dependent var */
    list[1] = dset->v - 1;

    aux = lsq(list, dset, OLS, OPT_A);
    if (aux.errcode) {
	err = aux.errcode;
	fprintf(stderr, "auxiliary regression failed\n");
    } else {
	double pval, trsq = aux.rsq * aux.nobs;
	int df = aux.ncoeff - pmod->ncoeff;

	if (df <= 0) {
	    /* couldn't add any regressors */
	    err = E_SINGULAR;
	    goto bailout;
	}

	pval = chisq_cdf_comp(df, trsq);
	aux.aux = aux_code;

	if (!(opt & OPT_I)) {
	    /* OPT_I = --silent */
	    if (opt & OPT_Q) {
		nonlin_test_header(aux_code, prn);
	    } else {
		printmodel(&aux, dset, OPT_NONE, prn);
		pputc(prn, '\n');
	    }
	    pprintf(prn, "  %s: TR^2 = %g,\n  ", _("Test statistic"), trsq);
	    pprintf(prn, "%s = P(%s(%d) > %g) = %g\n\n",
		    _("with p-value"), _("Chi-square"), df, trsq, pval);
	}

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = model_test_new((aux_code == AUX_SQ)?
				  GRETL_TEST_SQUARES : GRETL_TEST_LOGS);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LM);
		model_test_set_dfn(test, df);
		model_test_set_value(test, trsq);
		model_test_set_pvalue(test, chisq_cdf_comp(df, trsq));
		maybe_add_test_to_model(pmod, test);
	    }
	}

	record_test_result(trsq, pval);
    }

 bailout:

    clear_model(&aux);

    return err;
}

/**
 * nonlinearity_test:
 * @pmod: pointer to original model.
 * @dset: dataset struct.
 * @aux: AUX_SQ for squares or AUX_LOG for logs
 * @opt: if contains OPT_S, save test results to model; if
 * contains OPT_I, run silently.
 * @prn: gretl printing struct.
 *
 * Run an auxiliary regression to test @pmod for nonlinearity,
 * via the addition of either squares or logs of the original
 * indepdendent variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int nonlinearity_test (MODEL *pmod, DATASET *dset, ModelAuxCode aux,
		       gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *tmplist = NULL;
    const int orig_nvar = dset->v;
    int err = 0;

    if (!command_ok_for_model(ADD, 0, pmod)) {
	return E_NOTIMP;
    }

    if (pmod->ci == LOGISTIC || pmod->ci == LAD) {
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod, dset);
    if (err) {
	return err;
    }

    /* re-impose the sample that was in force when the original model
       was estimated */
    impose_model_smpl(pmod, dset);

    /* add squares or logs */
    tmplist = augment_regression_list(pmod->list, aux, dset, &err);
    if (err) {
	return err;
    } else if (tmplist[0] == pmod->list[0]) {
	/* no vars were added */
	if (aux == AUX_SQ) {
	    fprintf(stderr, "gretl: generation of squares failed\n");
	    err = E_SQUARES;
	} else if (aux == AUX_LOG) {
	    fprintf(stderr, "gretl: generation of logs failed\n");
	    err = E_LOGS;
	}
    }

    if (!err) {
	err = real_nonlinearity_test(pmod, tmplist, dset, aux,
				     opt, prn);
    }

    /* trash any extra variables generated (squares, logs) */
    dataset_drop_last_variables(dset, dset->v - orig_nvar);

    /* put back into dset what was there on input */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(tmplist);

    return err;
}

static int add_vars_missing (const MODEL *pmod,
			     const int *addvars,
			     const DATASET *dset)
{
    int i, vi, t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t) ||
	    (pmod->yhat != NULL && na(pmod->yhat[t]))) {
	    continue;
	}
	for (i=1; i<=addvars[0]; i++) {
	    vi = addvars[i];
	    if (na(dset->Z[vi][t])) {
		fprintf(stderr, "add: obs %d OK in model but missing for series %s\n",
			t+1, dset->varname[vi]);
		return E_MISSDATA;
	    }
	}
    }

    return 0;
}

static MODEL LM_add_test (MODEL *pmod, DATASET *dset, int *list,
			  gretlopt opt, PRN *prn)
{
    MODEL aux;
    int err;

    err = add_residual_to_dataset(pmod, dset);
    if (err) {
	gretl_model_init(&aux, dset);
	aux.errcode = err;
	return aux;
    }

    list[1] = dset->v - 1;

    aux = lsq(list, dset, OLS, OPT_A | OPT_Z);

    if (aux.errcode) {
	fprintf(stderr, "auxiliary regression failed\n");
    } else {
	int df = aux.ncoeff - pmod->ncoeff;

	if (df <= 0) {
	    /* couldn't add any regressors */
	    aux.errcode = E_SINGULAR;
	} else {
	    if (!(opt & (OPT_Q | OPT_I))) {
		aux.aux = AUX_ADD;
		printmodel(&aux, dset, OPT_S, prn);
	    }
	}
    }

    return aux;
}

/**
 * add_test_full:
 * @orig: pointer to original model.
 * @pmod: pointer to receive augmented model.
 * @addvars: list of variables to add to original model.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the new model, OPT_O to print its covariance matrix,
 * OPT_I for silent operation.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after adding the specified
 * variables, and records a joint test on the additional
 * variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_test_full (MODEL *orig, MODEL *pmod, const int *addvars,
		   DATASET *dset, gretlopt opt, PRN *prn)
{
    MODEL umod;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *biglist = NULL;
    const int orig_nvar = dset->v;
    int n_add = 0;
    int err = 0;

    if (orig == NULL || orig->list == NULL || addvars == NULL) {
	return E_DATA;
    }

    n_add = addvars[0];
    if (n_add == 0) {
	return E_NOADD;
    }

    if (incompatible_options(opt, OPT_A | OPT_L)) {
        /* --auto and --lm cannot be combined */
        return E_BADOPT;
    }

    if (!command_ok_for_model(ADD, opt, orig)) {
	return E_NOTIMP;
    }

    if ((opt & OPT_L) && pmod != NULL) {
	/* --lm option incompatible with "full" variant of add test */
	return E_BADOPT;
    }

    if (exact_fit_check(orig, prn)) {
	return 0;
    }

    /* check for changes in original list members */
    err = list_members_replaced(orig, dset);
    if (err) {
	return err;
    }

    /* check for NAs in add list relative to model */
    err = add_vars_missing(orig, addvars, dset);
    if (err) {
	return err;
    }

    if (!(opt & OPT_A)) {
        /* create augmented regression list, unless --auto (stepwise) */
        if (orig->ci == IVREG) {
            biglist = ivreg_list_add(orig->list, addvars, opt, &err);
        } else if (orig->ci == DPANEL) {
            biglist = panel_list_add(orig, addvars, &err);
        } else {
            biglist = gretl_list_add(orig->list, addvars, &err);
        }
    }

    if (err) {
	return err;
    }

    /* impose the sample range in force when the
       original model was estimated
    */
    impose_model_smpl(orig, dset);

    if (opt & OPT_A) {
        /* Do stepwise augmentation of the original model */
        MODEL (*stepwise_add) (MODEL *, const int *, DATASET *,
                               gretlopt, PRN *);

        stepwise_add = get_plugin_function("stepwise_add");
        if (stepwise_add == NULL) {
            err = E_FOPEN;
        } else {
            umod = stepwise_add(orig, addvars, dset, opt, prn);
        }
    } else if (opt & OPT_L) {
	/* run an LM test */
	umod = LM_add_test(orig, dset, biglist, opt, prn);
    } else {
	/* Run augmented regression, matching the original estimation
	   method; use OPT_Z to suppress the elimination of perfectly
	   collinear variables.
	*/
	gretlopt ropt = OPT_Z;

	if (opt & (OPT_Q | OPT_I)) {
	    /* not printing @umod, so pass along OPT_Q */
	    ropt |= OPT_Q;
	}
	umod = replicate_estimator(orig, biglist, dset, ropt, prn);
    }

    if (umod.errcode) {
	err = umod.errcode;
	errmsg(err, prn);
    }

    if (!(opt & OPT_A)) {
        if (!err && umod.ncoeff - orig->ncoeff != n_add) {
            gretl_errmsg_sprintf(_("Failed to add %d variable(s)"), n_add);
            err = E_DATA;
        }
        if (!err) {
            err = add_or_omit_compare(orig, &umod, addvars,
                                      dset, ADD, opt, prn);
        }
    }

    if (err || pmod == NULL) {
	clear_model(&umod);
    } else {
	*pmod = umod;
    }

    /* put dset back as it was on input */
    dataset_drop_last_variables(dset, dset->v - orig_nvar);
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(biglist);

    return err;
}

/**
 * add_test:
 * @pmod: pointer to model to be tested.
 * @addvars: list of variables to test.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the auxiliary model, OPT_I to suppress all printing
 * of results.
 * @prn: gretl printing struct.
 *
 * Performs an LM test on @pmod for the null hypothesis
 * that the @addvars variables do not contribute
 * significant explanatory power.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_test (MODEL *pmod, const int *addvars,
	      DATASET *dset, gretlopt opt, PRN *prn)
{
    return add_test_full(pmod, NULL, addvars, dset, opt, prn);
}

static int wald_omit_test (const int *list, MODEL *pmod,
			   const DATASET *dset, gretlopt opt,
			   PRN *prn)
{
    int *test = NULL;
    int err = 0;

    /* test validity of omissions */

    if (pmod->ci == IVREG) {
	test = ivreg_list_omit(pmod->list, list, opt, &err);
    } else if (pmod->ci == PANEL || pmod->ci == DPANEL) {
	test = panel_list_omit(pmod, list, &err);
    } else {
	test = gretl_list_omit(pmod->list, list, 2, &err);
    }

    if (!err) {
	free(test);
	err = add_or_omit_compare(pmod, NULL, list, dset, OMIT,
				  opt, prn);
    }

    return err;
}

/* Check whether coefficient @i corresponds to a variable
   that is removable from the model: this is the case if either
   (a) the @cands list is empty, or (b) coefficient @i is the
   coefficient on one of the variables in @cands.
*/

static int coeff_is_removable (const int *cands, const MODEL *pmod,
			       DATASET *dset, int i)
{
    int ret = 1;

    if (cands != NULL && cands[0] > 0) {
	const char *vname;
	int j, pj;

	ret = 0; /* reverse the presumption */

	for (j=1; j<=cands[0]; j++) {
	    vname = dset->varname[cands[j]];
	    pj = gretl_model_get_param_number(pmod, dset, vname);
	    if (pj == i) {
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

/* Determine if @pmod contains a variable with p-value
   greater than some cutoff alpha_max; and if so, remove
   this variable from @list.

   If the list @cands is non-empty then confine the search
   to candidate variables in that list.

   Returns 1 if a variable was dropped, else 0.
*/

static int auto_drop_var (const MODEL *pmod,
			  int *list, const int *cands,
			  DATASET *dset, double alpha_max,
			  int starting, gretlopt opt,
			  PRN *prn)
{
    double tstat, pv = 0.0, tmin = 4.0;
    int imin, imax = pmod->ncoeff;
    int i, k = -1;
    int ret = 0;

    if ((pmod->ci == LOGIT || pmod->ci == PROBIT) &&
	gretl_model_get_int(pmod, "ordered")) {
	/* FIXME other problematic cases? */
	imax = pmod->list[0] - 1;
    }

    /* if the constant is the sole regressors, allow it
       to be dropped */
    imin = (pmod->ncoeff == 1)? 0 : pmod->ifc;

    for (i=imin; i<imax; i++) {
	if (coeff_is_removable(cands, pmod, dset, i)) {
	    tstat = fabs(pmod->coeff[i] / pmod->sderr[i]);
	    if (tstat < tmin) {
		tmin = tstat;
		k = i;
	    }
	}
    }

    if (k >= 0) {
	pv = coeff_pval(pmod->ci, tmin, pmod->dfd);
    }

    if (pv > alpha_max) {
	char pname[VNAMELEN];
	int err;

	if (starting && !(opt & OPT_I)) {
	    /* not --silent */
	    pputc(prn, '\n');
	    pprintf(prn, _("Sequential elimination using two-sided alpha = %.2f"),
		    alpha_max);
	    pputs(prn, "\n\n");
	}

	gretl_model_get_param_name(pmod, dset, k, pname);
	err = gretl_list_delete_at_pos(list, k + 2);
	if (!err) {
	    if (!(opt & OPT_I)) {
		pprintf(prn, _(" Dropping %-16s (p-value %.3f)\n"), pname, pv);
	    }
	    ret = 1;
	}
    }

    return ret;
}

static void list_copy_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

/* run a loop in which the least significant variable is dropped
   from the regression list, provided its p-value exceeds some
   specified cutoff.  FIXME this probably still needs work for
   estimators other than OLS.  If @omitlist is non-empty the
   routine is confined to members of the list.
*/

static MODEL auto_omit (MODEL *orig, const int *omitlist,
			DATASET *dset, gretlopt opt,
			PRN *prn)
{
    MODEL omod;
    double amax;
    int *tmplist;
    int allgone = 0;
    int i, drop;
    int err = 0;

    gretl_model_init(&omod, dset);

    tmplist = gretl_list_copy(orig->list);
    if (tmplist == NULL) {
	omod.errcode = E_ALLOC;
	return omod;
    }

    amax = get_optval_double(OMIT, OPT_A, &err);
    if (err) {
	omod.errcode = err;
	return omod;
    }

    if (na(amax) || amax <= 0.0 || amax >= 1.0) {
	amax = 0.10;
    }

    drop = auto_drop_var(orig, tmplist, omitlist, dset, amax,
			 1, opt, prn);
    if (!drop) {
	/* nothing was dropped: set a "benign" error code */
	err = omod.errcode = E_NOOMIT;
    }

    for (i=0; drop > 0 && !allgone; i++) {
	if (i > 0) {
	    set_reference_missmask_from_model(orig);
	}
	omod = replicate_estimator(orig, tmplist, dset, OPT_A, prn);
	err = omod.errcode;
	if (err) {
	    fprintf(stderr, "auto_omit: error %d from replicate_estimator\n",
		    err);
	    drop = 0; /* will break */
	} else {
	    list_copy_values(tmplist, omod.list);
	    drop = auto_drop_var(&omod, tmplist, omitlist, dset,
				 amax, 0, opt, prn);
	    if (drop && omod.ncoeff == 1) {
		allgone = 1; /* will break */
	    }
	    clear_model(&omod);
	}
    }

    if (allgone) {
	pputc(prn, '\n');
	pprintf(prn, _("No coefficient has a p-value less than %g"), amax);
	pputc(prn, '\n');
    } else if (!err) {
	/* re-estimate the final model without the auxiliary flag */
	gretlopt ropt = OPT_NONE;

	if (opt & (OPT_Q | OPT_I)) {
	    /* not printing @omod, so pass along OPT_Q */
	    ropt |= OPT_Q;
	}
	set_reference_missmask_from_model(orig);
	omod = replicate_estimator(orig, tmplist, dset, ropt, prn);
    }

    free(tmplist);

    return omod;
}

/* create reduced list for "omit" test on model, based on
   the list of variables to be dropped, @omitvars
*/

static int make_short_list (MODEL *orig, const int *omitvars,
			    gretlopt opt, int **plist)
{
    int *list = NULL;
    int err = 0;

    if (omitvars == NULL || omitvars[0] == 0) {
	return E_PARSE;
    }

    if (orig->ci == IVREG) {
	list = ivreg_list_omit(orig->list, omitvars, opt, &err);
    } else if (orig->ci == PANEL || orig->ci == DPANEL) {
	list = panel_list_omit(orig, omitvars, &err);
    } else {
	list = gretl_list_omit(orig->list, omitvars, 2, &err);
    }

    if (list != NULL && list[0] == 1) {
	/* only the dependent variable would be left */
	err = E_NOVARS;
    }

    *plist = list;

    return err;
}

static int omit_options_inconsistent (MODEL *pmod, gretlopt opt)
{
    if (opt & OPT_B) {
	/* 2sls: omitting variable as instrument */
	if (opt & (OPT_W | OPT_A)) {
	    /* can't use Wald method on original VCV,
	       and can't do --auto */
	    return 1;
	} else if (pmod->ci != IVREG) {
	    return 1;
	}
    }

    if ((opt & OPT_A) && (opt & OPT_W)) {
	/* --auto and --wald options incompatible */
	return 1;
    }

    return 0;
}

static int omit_test_precheck (MODEL *pmod, gretlopt opt)
{
    int err = 0;

    if (pmod == NULL || pmod->list == NULL) {
	err = E_DATA;
    } else if (!command_ok_for_model(OMIT, 0, pmod)) {
	err = E_NOTIMP;
    } else if (omit_options_inconsistent(pmod, opt)) {
	err = E_BADOPT;
    }

    return err;
}

/**
 * omit_test_full:
 * @orig: pointer to original model.
 * @pmod: pointer to receive new model, with vars omitted.
 * @omitvars: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of the new model, OPT_O to print its covariance matrix,
 * OPT_I for silent operation; for OPT_A, see below.
 * @prn: gretl printing struct.
 *
 * Re-estimate a given model after removing the variables
 * specified in @omitvars.  Or if OPT_A is given, proceed
 * sequentially, at each step dropping the least significant
 * variable provided its p-value is above a certain threshold
 * (currently 0.10, two-sided).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test_full (MODEL *orig, MODEL *pmod, const int *omitvars,
		    DATASET *dset, gretlopt opt, PRN *prn)
{
    MODEL rmod;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *tmplist = NULL;
    int err;

    err = omit_test_precheck(orig, opt);
    if (err) {
	return err;
    }

    /* check that vars to omit have not been redefined */
    if ((err = list_members_replaced(orig, dset))) {
	return err;
    }

    if (!(opt & OPT_A)) {
	/* not doing auto-omit */
	err = make_short_list(orig, omitvars, opt, &tmplist);
	if (err) {
	    free(tmplist);
	    return err;
	}
    }

    /* impose the sample range used for the original model */
    impose_model_smpl(orig, dset);

    /* set the mask for missing obs within the sample range, based
       on the original model */
    set_reference_missmask_from_model(orig);

    if (opt & OPT_A) {
#if 1 /* not yet */
        /* Do stepwise reduction of the original model */
        MODEL (*stepwise_omit) (MODEL *, const int *, DATASET *,
                               gretlopt, PRN *);

        stepwise_omit = get_plugin_function("stepwise_omit");
        if (stepwise_omit == NULL) {
            err = E_FOPEN;
        } else {
            rmod = stepwise_omit(orig, omitvars, dset, opt, prn);
        }
#else
	rmod = auto_omit(orig, omitvars, dset, opt, prn);
#endif
    } else {
	gretlopt ropt = OPT_NONE;

	if (opt & (OPT_Q | OPT_I)) {
	    /* not printing @rmod, so pass along OPT_Q */
	    ropt |= OPT_Q;
	}
	rmod = replicate_estimator(orig, tmplist, dset, ropt, prn);
    }

    err = rmod.errcode;

    if (err) {
	if (err == E_NOOMIT && (opt & OPT_I)) {
	    ; /* --silent: keep quiet */
	} else {
	    errmsg(err, prn);
	}
    } else {
	MODEL *newmod = NULL;
	int minpos = 2;
	int *omitlist = NULL;

	if (orig->ci == DPANEL) {
	    /* skip AR spec, separator, and dep var */
	    minpos = 4;
	}

	if (rmod.list == NULL) {
	    /* handle the case where sequential elimination
	       ended up dropping all regressors */
	    int *minlist = gretl_list_new(1);

	    minlist[1] = orig->list[1];
	    omitlist = gretl_list_diff_new(orig->list, minlist, minpos);
	    free(minlist);
	} else {
	    newmod = &rmod;
	    omitlist = gretl_list_diff_new(orig->list, rmod.list, minpos);
	}
	if (omitlist != NULL) {
	    err = add_or_omit_compare(orig, newmod, omitlist,
				      dset, OMIT, opt, prn);
	    free(omitlist);
	}
    }

    if (err || pmod == NULL) {
	clear_model(&rmod);
    } else {
	*pmod = rmod;
    }

    /* put back into dset what was there on input */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(tmplist);

    return err;
}

/**
 * omit_test:
 * @pmod: pointer to model to be tested.
 * @omitvars: list of variables to test.
 * @dset: dataset struct.
 * @opt: can contain OPT_Q (quiet) to suppress printing
 * of results.
 * @prn: gretl printing struct.
 *
 * Performs a Wald test on @pmod for the null hypothesis
 * that the @omitvars variables do not contribute
 * explanatory power.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int omit_test (MODEL *pmod, const int *omitvars,
	       DATASET *dset, gretlopt opt, PRN *prn)
{
    int err = omit_test_precheck(pmod, opt);

    if (!err) {
	err = wald_omit_test(omitvars, pmod, dset, opt, prn);
    }

    return err;
}

/**
 * get_DW_pvalue_for_model:
 * @pmod: model to be tested.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Computes the p-value for the Durbin-Watson statistic for the
 * given model, using the Imhof method.
 *
 * Returns: the p-value, or #NADBL on error.
 */

double get_DW_pvalue_for_model (MODEL *pmod, DATASET *dset,
				int *err)
{
    MODEL dwmod;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *list = NULL;
    double pv = NADBL;

    /* maybe this has already been done? */
    pv = gretl_model_get_double(pmod, "dw_pval");
    if (!na(pv)) {
	return pv;
    }

    if (pmod->ci == PANEL && panel_DW_pval_ok(pmod)) {
	/* Use Bhargava et al approximation */
	return BFN_panel_DW_pvalue(pmod, dset, err);
    }

    if (dset == NULL || dset->Z == NULL) {
	*err = E_NODATA;
    } else if (pmod == NULL || pmod->list == NULL) {
	*err = E_DATA;
    } else if ((pmod->ci != OLS && pmod->ci != PANEL) ||
	       na(pmod->dw) || model_has_missing_obs(pmod)) {
	*err = E_BADSTAT;
    } else {
	/* check that relevant vars have not been redefined */
	*err = list_members_replaced(pmod, dset);
    }

    if (!*err) {
	list = gretl_list_copy(pmod->list);
	if (list == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (*err) {
	return NADBL;
    }

    gretl_model_init(&dwmod, dset);

    /* impose the sample range used for the original model */
    impose_model_smpl(pmod, dset);

    dwmod = replicate_estimator(pmod, list, dset, OPT_A | OPT_I, NULL);
    *err = dwmod.errcode;

    if (!*err) {
	pv = gretl_model_get_double(&dwmod, "dw_pval");
	if (!na(pv)) {
	    /* record result on the incoming model */
	    gretl_model_set_double(pmod, "dw_pval", pv);
	}
    }

    /* put back into dset what was there on input */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    clear_model(&dwmod);
    free(list);

    return pv;
}

/**
 * reset_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if contains %OPT_S, save test results to model. %OPT_Q
 * suppresses the printout of the auxiliary regression. %OPT_U and
 * %OPT_C stand for "squares only" and "cubes only", respectively.
 * %OPT_I produces silent operation.
 * @prn: gretl printing struct.
 *
 * Carries out and prints Ramsey's RESET test for model specification.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int reset_test (MODEL *pmod, DATASET *dset,
		gretlopt opt, PRN *prn)
{
    int *newlist = NULL;
    MODEL aux;
    double RF;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int i, t, orig_v = dset->v;
    int addv, use_square, use_cube;
    int robust = 0;
    const char *mode;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if ((opt & OPT_R) || (pmod->opt & OPT_R)) {
	robust = 1;
    }

    err = incompatible_options(opt, OPT_C | OPT_U);

    if (err) {
	return err;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    use_square = !(opt & OPT_C); /* not cubes-only */
    use_cube = !(opt & OPT_U);   /* not squares-only */

    gretl_model_init(&aux, dset);

    if (opt & OPT_U) {
	addv = 1;
	mode = N_("squares only");
    } else if (opt & OPT_C) {
	addv = 1;
	mode = N_("cubes only");
    } else {
	addv = 2;
	mode = N_("squares and cubes");
    }

    impose_model_smpl(pmod, dset);

    if (pmod->ncoeff + addv >= dset->t2 - dset->t1) {
	/* can't run aux regression */
	err = E_DF;
    }

    if (!err) {
	newlist = gretl_list_new(pmod->list[0] + addv);
	if (newlist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	for (i=1; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_series(dset, addv)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add scaled yhat^2 and/or yhat^3 to data set */
	int vs = orig_v;
	int vc = (opt & OPT_C)? orig_v : orig_v + 1;
	int k = pmod->list[0] + 1;
        double s;

        s = gretl_stddev(pmod->t1, pmod->t2, pmod->yhat);

	for (t = pmod->t1; t<=pmod->t2; t++) {
	    double x = pmod->yhat[t] / s;

	    if (use_square) {
		dset->Z[vs][t] = x * x;
	    }
	    if (use_cube) {
		dset->Z[vc][t] = x * x * x;
	    }
	}

	if (use_square) {
	    strcpy(dset->varname[vs], "yhat^2");
	    newlist[k++] = vs;
	}
	if (use_cube) {
	    strcpy(dset->varname[vc], "yhat^3");
	    newlist[k++] = vc;
	}
    }

    if (!err) {
	gretlopt auxopt = OPT_A;

	if (robust) {
	    auxopt |= OPT_R;
	}
	aux = lsq(newlist, dset, OLS, auxopt);
	err = aux.errcode;
	if (err) {
	    errmsg(aux.errcode, prn);
	}
    }

    if (!err) {
	int silent = (opt & OPT_I);
	double pval;

	aux.aux = AUX_RESET;

	if (!silent) {
	    if (!(opt & OPT_Q)) {
		printmodel(&aux, dset, OPT_NONE, prn);
	    } else {
		if (!(opt & OPT_G)) {
		    /* GUI special; see gui2/library.c */
		    pputc(prn, '\n');
		}
		if (robust) {
		    pputs(prn, _("Robust RESET test for specification"));
		} else {
		    pputs(prn, _("RESET test for specification"));
		}
		pprintf(prn, " (%s)\n", _(mode));
	    }
	}

	if (robust) {
	    /* FIXME wald_omit_F(omitlist, pmod); */
	    int *omitlist = gretl_list_diff_new(aux.list, pmod->list, 2);

	    RF = wald_omit_F(omitlist, &aux);
	    free(omitlist);
	} else {
	    RF = ((pmod->ess - aux.ess) / addv) / (aux.ess / aux.dfd);
	}
	pval = snedecor_cdf_comp(addv, aux.dfd, RF);

	if (!silent) {
	    const char *h0 = get_h0_string_for_test(GRETL_TEST_RESET);

	    pprintf(prn, "%s: %s\n", _("Null hypothesis"), _(h0));
	    pprintf(prn, "%s: F = %f,\n", _("Test statistic"), RF);
	    pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"),
		    addv, aux.dfd, RF, pval);
	    pputc(prn, '\n');
	}

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_RESET);
	    gretlopt topt = OPT_NONE;

	    if (test != NULL) {
		if (robust) {
		    topt = OPT_R;
		}
		if (opt & OPT_U) {
		    topt = OPT_U;
		} else if (opt & OPT_C) {
		    topt = OPT_C;
		}
		model_test_set_teststat(test, GRETL_STAT_RESET);
		model_test_set_dfn(test, addv);
		model_test_set_dfd(test, aux.dfd);
		model_test_set_value(test, RF);
		model_test_set_pvalue(test, pval);
		model_test_set_opt(test, topt);
		maybe_add_test_to_model(pmod, test);
	    }
	}

	record_test_result(RF, pval);
    }

    free(newlist);
    dataset_drop_last_variables(dset, addv);
    clear_model(&aux);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

static void bg_test_header (int order, PRN *prn, int ivreg)
{
    if (ivreg) {
	pprintf(prn, "\n%s ", _("Godfrey (1994) test for"));
    } else {
	pprintf(prn, "\n%s ", _("Breusch-Godfrey test for"));
    }

    if (order > 1) {
	pprintf(prn, "%s %d\n", _("autocorrelation up to order"),
		order);
    } else {
	pprintf(prn, "%s\n", _("first-order autocorrelation"));
    }

    pputc(prn, '\n');
}

static double ivreg_autocorr_wald_stat (MODEL *aux, int order, int *err)
{
    gretl_vector *b = gretl_vector_alloc(order);
    gretl_matrix *V1 = gretl_matrix_alloc(order, order);
    gretl_vector *WT = gretl_vector_alloc(1);
    gretl_matrix *V0 = NULL;
    double x = NADBL;
    int i, j, ki, kj;

    if (b == NULL || V1 == NULL || WT == NULL) {
	*err = E_ALLOC;
    } else {
	V0 = gretl_model_get_matrix(aux, M_VCV, err);
    }

    if (!*err) {
	ki = aux->ncoeff - order;

	for (i=0; i<order; i++) {
	    x = aux->coeff[ki];
	    gretl_vector_set(b, i, x);
	    kj = ki;
	    for (j=i; j<order; j++) {
		x = gretl_matrix_get(V0, ki, kj);
		gretl_matrix_set(V1, i, j, x);
		gretl_matrix_set(V1, j, i, x);
		kj++;
	    }
	    ki++;
	}
	*err = gretl_invert_symmetric_matrix(V1);
    }

    if (!*err) {
	gretl_matrix_qform(b, GRETL_MOD_NONE, V1, WT, GRETL_MOD_NONE);
	x = gretl_vector_get(WT, 0) / order;
    }

    gretl_vector_free(WT);
    gretl_vector_free(b);
    gretl_matrix_free(V1);
    gretl_matrix_free(V0);

    return x;
}

/**
 * ivreg_autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose.
 * @prn: gretl printing struct.
 *
 * Tests the given IV model for autocorrelation of order equal
 * to the specified value, or equal to the frequency of the data if
 * the supplied @order is zero, as per Godfrey (1994), "Testing for
 * Serial Correlation by Variable Addition in Dynamic Models Estimated
 * by Instrumental Variables", RES. Note that none of the
 * asymptotically equivalent tests given on page 553 is used
 * here. Instead, we estimate the model augmented with lags and then
 * perform a Wald-type test. The resulting chi-square statistic is
 * divided by its degrees of freedom as a finite-sample adjustment and
 * compared to an F distribution.
 *
 * Returns: 0 on successful completion, error code on error.
 */

static int ivreg_autocorr_test (MODEL *pmod, int order,
				DATASET *dset, gretlopt opt,
				PRN *prn)
{
    int smpl_t1 = dset->t1;
    int smpl_t2 = dset->t2;
    int v = dset->v;
    int *addlist = NULL;
    int *testlist = NULL;
    double x = 0, pval = 1;
    MODEL aux;
    int i, t;
    int err = 0;

    if (dataset_is_panel(dset)) {
	return E_NOTIMP;
    }

    if (model_has_missing_obs(pmod)) {
	return E_MISSDATA;
    }

    /* impose original sample range */
    impose_model_smpl(pmod, dset);

    gretl_model_init(&aux, dset);

    if (order <= 0) {
	order = dset->pd;
    }

    if (pmod->ncoeff + order >= dset->t2 - dset->t1) {
	return E_DF;
    }

    addlist = gretl_list_new(order);

    if (addlist == NULL) {
	err = E_ALLOC;
    } else {
	err = add_residual_to_dataset(pmod, dset);
    }

    if (!err) {
	/* add lags of residual */
	for (i=1; i<=order && !err; i++) {
	    int lnum;

	    lnum = laggenr(v, i, dset);

	    if (lnum < 0) {
		gretl_errmsg_set(_("lagging uhat failed"));
		err = E_LAGS;
	    } else {
		/* set the first entries to 0 for compatibility with Godfrey (1994)
		   and PcGive (not perfect) */
		for (t=smpl_t1; t<smpl_t1+i; t++) {
		    dset->Z[lnum][t] = 0;
		}
		addlist[i] = lnum;
	    }
	}
	if (!err) {
	    /* compose augmented regression list */
	    testlist = ivreg_list_add(pmod->list, addlist, OPT_B, &err);
	}
    }

    if (!err) {
	gretlopt ivopt = OPT_A;

	transcribe_option_flags(&ivopt, pmod->opt,
				OPT_L | OPT_G | OPT_R);
	aux = ivreg(testlist, dset, ivopt);
	err = aux.errcode;
    }

    if (!err) {
	x = ivreg_autocorr_wald_stat(&aux, order, &err);
    }

    if (!err) {
	aux.aux = AUX_AR;
	gretl_model_set_int(&aux, "BG_order", order);
	pval = snedecor_cdf_comp(order, aux.nobs - pmod->ncoeff - order, x);

	if (opt & OPT_Q) {
	    bg_test_header(order, prn, 1);
	} else {
	    printmodel(&aux, dset, OPT_S, prn);
	}

	pputc(prn, '\n');
	pprintf(prn, "%s: Pseudo-LMF = %f,\n", _("Test statistic"), x);
	pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"),
		order, aux.nobs - pmod->ncoeff, x, pval);
	pputc(prn, '\n');
	record_test_result(x / order, pval);

	if (opt & OPT_S) {
	    ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, aux.nobs - pmod->ncoeff);
		model_test_set_order(test, order);
		model_test_set_value(test, x);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }
	}
    }

    free(addlist);
    free(testlist);

    dataset_drop_last_variables(dset, dset->v - v);
    clear_model(&aux);

    /* reset sample as it was */
    dset->t1 = smpl_t1;
    dset->t2 = smpl_t2;

    return err;
}

static int lb_autocorr_test (MODEL *pmod, int order,
			     gretlopt opt, PRN *prn)
{
    double lb, pval = NADBL;
    int df, k;
    int err = 0;

    k = arma_model_get_n_arma_coeffs(pmod);
    df = order - k;

    if (df <= 0) {
	gretl_errmsg_set(_("Insufficient degrees of freedom for test"));
	return E_DATA;
    }

    lb = ljung_box(order, pmod->t1, pmod->t2, pmod->uhat, &err);

    if (!err) {
	pval = chisq_cdf_comp(df, lb);
	if (na(pval)) {
	    err = E_DATA;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error calculating Ljung-Box statistic"));
    } else {
	pputc(prn, '\n');
	pprintf(prn, _("Test for autocorrelation up to order %d"),
		order);
	pputs(prn, "\n\n");
	pprintf(prn, "Ljung-Box Q' = %g,\n", lb);
	pprintf(prn, "%s = P(%s(%d) > %g) = %#.4g\n", _("with p-value"),
		_("Chi-square"), df, lb, chisq_cdf_comp(df, lb));
	pputc(prn, '\n');
	record_test_result(lb, pval);
    }

    if (!err && (opt & OPT_S)) {
	ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_LB_CHISQ);
	    model_test_set_dfn(test, df);
	    model_test_set_order(test, order);
	    model_test_set_value(test, lb);
	    model_test_set_pvalue(test, pval);
	    maybe_add_test_to_model(pmod, test);
	}
    }

    return err;
}

/**
 * regular_autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose; if OPT_I, be silent.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero. Prints TR^2 and LMF test statistics.
 *
 * Returns: 0 on successful completion, error code on error.
 */

static int regular_autocorr_test (MODEL *pmod, int order, DATASET *dset,
				  gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *newlist = NULL;
    MODEL aux;
    double RSSxe, RSSx = pmod->ess;
    int i, t, n = dset->n, v = dset->v;
    double trsq, LMF, lb, pval = 1.0;
    int err = 0;

    gretl_model_init(&aux, dset);

    if (order <= 0) {
	order = dset->pd;
    }

    if (pmod->ncoeff + order >= pmod->t2 - pmod->t1) {
	return E_DF;
    }

    newlist = gretl_list_new(pmod->list[0] + order);

    if (newlist == NULL) {
	err = E_ALLOC;
    } else {
	newlist[0] = pmod->list[0] + order;
	for (i=2; i<=pmod->list[0]; i++) {
	    newlist[i] = pmod->list[i];
	}
	if (dataset_add_series(dset, 1 + order)) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	/* add uhat to data set: substitute zeros for
	   pre-sample values */
	for (t=0; t<n; t++) {
	    if (t < pmod->t1) {
		dset->Z[v][t] = 0.0;
	    } else {
		dset->Z[v][t] = pmod->uhat[t];
	    }
	}
	strcpy(dset->varname[v], "uhat");
	series_set_label(dset, v, _("residual"));
	/* then order lags of same */
	for (i=1; i<=order; i++) {
	    int s, lv = v + i;
	    double ul;

	    sprintf(dset->varname[lv], "uhat_%d", i);
	    newlist[pmod->list[0] + i] = lv;
	    for (t=0; t<dset->n; t++) {
		s = t - i;
		if (s < 0) {
		    dset->Z[lv][t] = 0.0;
		} else {
		    ul = dset->Z[v][s];
		    dset->Z[lv][t] = (na(ul))? 0.0 : ul;
		}
	    }
	}
    }

    /* LMF apparatus: see Kiviet, Review of Economic Studies,
       53/2, 1986, equation (5), p. 245.
    */

    if (!err) {
	/* regression on [X~E], using original sample */
	impose_model_smpl(pmod, dset);
	newlist[1] = v;
	aux = lsq(newlist, dset, OLS, OPT_A);
	err = aux.errcode;
	if (err) {
	    errmsg(err, prn);
	} else {
	    RSSxe = aux.ess;
	}
    }

    if (!err) {
	int dfd = aux.nobs - pmod->ncoeff - order;
	int lberr;

	aux.aux = AUX_AR;
	gretl_model_set_int(&aux, "BG_order", order);
	trsq = aux.rsq * aux.nobs;
	LMF = ((RSSx - RSSxe) / RSSxe) * dfd / order;
	pval = snedecor_cdf_comp(order, dfd, LMF);

	if (pmod->aux == AUX_VAR || (opt & OPT_I)) {
	    ; /* don't print anything here */
	} else {
	    if (opt & OPT_Q) {
		bg_test_header(order, prn, 0);
	    } else {
		printmodel(&aux, dset, OPT_NONE, prn);
		pputc(prn, '\n');
	    }
	    pprintf(prn, "%s: LMF = %f,\n", _("Test statistic"), LMF);
	    pprintf(prn, "%s = P(F(%d,%d) > %g) = %.3g\n", _("with p-value"),
		    order, aux.nobs - pmod->ncoeff - order, LMF, pval);

	    pprintf(prn, "\n%s: TR^2 = %f,\n",
		    _("Alternative statistic"), trsq);
	    pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n\n", _("with p-value"),
		    _("Chi-square"), order, trsq, chisq_cdf_comp(order, trsq));

	    lb = ljung_box(order, pmod->t1, pmod->t2, dset->Z[v], &lberr);
	    if (!na(lb)) {
		pprintf(prn, "Ljung-Box Q' = %g,\n", lb);
		pprintf(prn, "%s = P(%s(%d) > %g) = %.3g\n", _("with p-value"),
			_("Chi-square"), order, lb, chisq_cdf_comp(order, lb));
	    }
	    pputc(prn, '\n');
	}

	record_test_result(LMF, pval);

	if (opt & OPT_S) {
	    /* save the test onto @pmod */
	    ModelTest *test = model_test_new(GRETL_TEST_AUTOCORR);

	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_LMF);
		model_test_set_dfn(test, order);
		model_test_set_dfd(test, aux.nobs - pmod->ncoeff - order);
		model_test_set_order(test, order);
		model_test_set_value(test, LMF);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }
	}
    }

    free(newlist);
    dataset_drop_last_variables(dset, dset->v - v);
    clear_model(&aux);

    /* reset sample as it was */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/**
 * autocorr_test:
 * @pmod: pointer to model to be tested.
 * @order: lag order for test.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose; if OPT_I, be silent.
 * @prn: gretl printing struct.
 *
 * Tests the given model for autocorrelation of order equal to
 * the specified value, or equal to the frequency of the data if
 * the supplied @order is zero.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int autocorr_test (MODEL *pmod, int order, DATASET *dset,
		   gretlopt opt, PRN *prn)
{
    if (pmod->ci == IVREG) {
	return ivreg_autocorr_test(pmod, order, dset, opt, prn);
    } else if (pmod->ci == ARMA) {
	return lb_autocorr_test(pmod, order, opt, prn);
    } else if (gretl_is_regular_panel_model(pmod)) {
	if (dset->pd < 3) {
	    /* time series not long enough */
	    return E_NOTIMP;
	} else {
	    return panel_autocorr_test(pmod, dset, opt, prn);
	}
    }

    if (pmod->ci != OLS && pmod->ci != VAR) {
	return E_NOTIMP;
    } else if (model_has_missing_obs(pmod)) {
	return E_MISSDATA;
    }

    return regular_autocorr_test(pmod, order, dset, opt, prn);
}

static int chow_active (int split, const double *x, int t)
{
    if (x != NULL) {
	return x[t] == 1.0;
    } else {
	return (t >= split);
    }
}

static int *get_break_list (const MODEL *pmod, int ci, int *err)
{
    const char *lname = get_optval_string(ci, OPT_L);
    int *list = NULL;

    if (lname == NULL) {
	*err = E_DATA;
    } else {
	list = get_list_by_name(lname);
    }

    if (list != NULL) {
	if (list[0] == 0) {
	    *err = E_DATA;
	} else {
	    int i;

	    for (i=1; i<=list[0] && !*err; i++) {
		if (!in_gretl_list(pmod->list, list[i])) {
		    *err = E_DATA;
		} else if (list[i] == 0) {
		    *err = E_DATA;
		}
	    }
	}
    }

    return list;
}

/* compose list of variables to be used for the Chow test and add
   them to the data set */

static int *
make_chow_list (const MODEL *pmod, DATASET *dset,
		int split, int dumv, int ci,
		gretlopt opt, int *err)
{
    int *chowlist = NULL;
    int *brklist = NULL;
    int l0 = pmod->list[0];
    int ninter = 0, newvars = 0;
    int havedum = (dumv > 0);
    int i, t, v = dset->v;

    if (havedum && in_gretl_list(pmod->list, dumv)) {
	gretl_errmsg_sprintf(_("The model already contains %s"),
			     dset->varname[dumv]);
	*err = E_DATA;
	return NULL;
    }

    if (opt & OPT_L) {
	/* --limit-to */
	brklist = get_break_list(pmod, ci, err);
	if (*err) {
	    return NULL;
	} else {
	    ninter = brklist[0];
	}
    } else {
	/* number of interaction terms */
	ninter = pmod->ncoeff - pmod->ifc;
    }

    newvars = ninter + 1 - havedum;

    if (dataset_add_series(dset, newvars)) {
	*err = E_ALLOC;
    } else {
	chowlist = gretl_list_new(pmod->list[0] + ninter + 1);
	if (chowlist == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	const double *cdum = NULL;

	for (i=1; i<=l0; i++) {
	    chowlist[i] = pmod->list[i];
	}

	if (dumv > 0) {
	    /* we have a user-supplied dummy var */
	    cdum = dset->Z[dumv];
	} else {
	    /* generate the split variable */
	    for (t=0; t<dset->n; t++) {
		dset->Z[v][t] = (double) (t >= split);
	    }
	    strcpy(dset->varname[v], "splitdum");
	    series_set_label(dset, v, _("dummy variable for Chow test"));
	}

	chowlist[l0 + 1] = (dumv > 0)? dumv : v;

	/* and the interaction terms */
	for (i=0; i<ninter; i++) {
	    int pv, sv = v + i + 1 - havedum;

	    if (brklist != NULL) {
		pv = brklist[i+1];
	    } else {
		pv = pmod->list[i + 2 + pmod->ifc];
	    }

	    for (t=0; t<dset->n; t++) {
		if (model_missing(pmod, t)) {
		    dset->Z[sv][t] = NADBL;
		} else if (chow_active(split, cdum, t)) {
		    dset->Z[sv][t] = dset->Z[pv][t];
		} else {
		    dset->Z[sv][t] = 0.0;
		}
	    }

	    if (havedum) {
		sprintf(dset->varname[sv], "%.2s_", dset->varname[dumv]);
	    } else {
		strcpy(dset->varname[sv], "sd_");
	    }
	    strncat(dset->varname[sv], dset->varname[pv], VNAMELEN - 4);
	    chowlist[l0 + 2 + i] = sv;
	}
    }

    return chowlist;
}

static void write_plot_x_range (const double *x, int t1, int t2,
				FILE *fp)
{
    double xmin0, xmin, xmax0, xmax;
    double xrange;

    gretl_minmax(t1, t2, x, &xmin0, &xmax0);
    xrange = xmax0 - xmin0;
    xmin = xmin0 - xrange * .025;
    if (xmin0 >= 0.0 && xmin < 0.0) {
	xmin = 0.0;
    }
    xmax = xmax0 + xrange * .025;
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
}

/* try to ensure that the critical value, @cv, is not squashed
   against the bottom or top of the plot
*/

static void write_plot_y_range (const double *y, int T,
				double cv, FILE *fp)
{
    double ymin, ymax, yrange;
    double f = 0.025;

    gretl_minmax(0, T, y, &ymin, &ymax);
    if (cv < ymin) {
	ymin = cv;
	f = 0.05;
    } else if (cv > ymax) {
	ymax = cv;
	f = 0.05;
    }
    yrange = ymax - ymin;
    ymin -= f * yrange;
    ymax += f * yrange;

    fprintf(fp, "set yrange [%.10g:%.10g]\n", ymin, ymax);
}

static int QLR_graph (const double *testvec, int t1, int t2,
		      const DATASET *dset, int df, int robust)
{
    const double *x = gretl_plotx(dset, OPT_NONE);
    const char *titles[] = {
	N_("Chow F-test for break"),
	N_("Robust Wald test for break")
    };
    const char *title;
    double (*qlr_critval) (int);
    double critval = NADBL;
    GptFlags flags = 0;
    FILE *fp;
    int t;
    int err = 0;

    set_effective_plot_ci(QLRTEST);

    if (dataset_is_time_series(dset)) {
	flags = GPT_TS | GPT_LETTERBOX;
    }

    fp = open_plot_input_file(PLOT_REGULAR, flags, &err);
    if (err) {
	return err;
    }

    qlr_critval = get_plugin_function("qlr_critval_15_05");
    if (qlr_critval != NULL) {
	critval = qlr_critval(df);
    }

    print_keypos_string(GP_KEY_LEFT | GP_KEY_TOP, fp);

    if (robust) {
	title = titles[1];
    } else {
	title = titles[0];
	if (!na(critval)) {
	    critval /= df;
	}
    }

    gretl_push_c_numeric_locale();

    write_plot_x_range(x, t1, t2, fp);

    if (!na(critval)) {
	write_plot_y_range(testvec, t2 - t1, critval, fp);
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 title \"%s\" w lines , \\\n"
		"%g title \"%s\" w lines dt 2\n",
		_(title), critval, _("5% QLR critical value"));
    } else {
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 title \"%s\" w lines\n", _(title));
    }

    for (t=t1; t<=t2; t++) {
	fprintf(fp, "%g %g\n", x[t], testvec[t-t1]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    reset_effective_plot_ci();

    return err;
}

static void save_QLR_test (MODEL *pmod, const char *datestr,
			   double X2, double pval, int df)
{
    ModelTest *test = model_test_new(GRETL_TEST_QLR);

    if (test != NULL) {
	model_test_set_teststat(test, GRETL_STAT_SUP_WALD);
	model_test_set_param(test, datestr);
	model_test_set_value(test, X2);
	model_test_set_pvalue(test, pval);
	model_test_set_dfn(test, df);
	maybe_add_test_to_model(pmod, test);
    }
}

/* for internal use by the "qlrtest" command */

static double get_QLR_pval (double test, int df,
			    int k1, int k2,
			    MODEL *pmod)
{
    double (*qlr_asy_pvalue) (double, int, double);
    double pi_1, pi_2, lam0;
    double pval;

    qlr_asy_pvalue = get_plugin_function("qlr_asy_pvalue");

    if (qlr_asy_pvalue == NULL) {
	return NADBL;
    }

    pi_1 = (k1 - pmod->t1 + 1) / (double) pmod->nobs;
    pi_2 = (k2 - pmod->t1 + 1) / (double) pmod->nobs;
    lam0 = (pi_2*(1.0 - pi_1)) / (pi_1*(1.0 - pi_2));

    pval = (*qlr_asy_pvalue) (test, df, lam0);

#if 0
    fprintf(stderr, "k1 = %d, k2 = %d, pi_1 = %g pi_2 = %g\n",
	    k1, k2, pi_1, pi_2);
    fprintf(stderr, "lambda0 = %g, pi0 = %g, test = %g, pval = %g\n",
	    lam0, 1.0/(1 + sqrt(lam0)), test, pval);
#endif

    return pval;
}

/* for use by the qlrpval() function */

double QLR_pval (double X2, int df, double p1, double p2)
{
    double (*qlr_asy_pvalue) (double, int, double);
    double lamda, pval;

    if (X2 < 0 || df <= 0 || p1 <= 0 || p2 <= p1 || p2 >= 1) {
	return NADBL;
    }

    qlr_asy_pvalue = get_plugin_function("qlr_asy_pvalue");

    if (qlr_asy_pvalue == NULL) {
	return NADBL;
    }

    lamda = (p2*(1.0 - p1)) / (p1*(1.0 - p2));
    pval = (*qlr_asy_pvalue) (X2, df, lamda);

    return pval;
}

static void QLR_print_result (double test,
			      double pval,
			      const char *datestr,
			      int dfn, int dfd,
			      int robust,
			      PRN *prn)
{
    pputs(prn, _("Quandt likelihood ratio test for structural break at an "
		 "unknown point,\nwith 15 percent trimming"));
    pputs(prn, ":\n");
    if (robust) {
	pprintf(prn, _("The maximum Wald test = %g occurs "
		       "at observation %s"), test, datestr);
    } else {
	pprintf(prn, _("The maximum F(%d, %d) = %g occurs "
		       "at observation %s"), dfn, dfd, test, datestr);
	test *= dfn;
    }
    pputc(prn, '\n');

    if (!na(pval)) {
	pprintf(prn, _("Asymptotic p-value = %.6g for "
		       "chi-square(%d) = %g"),
		pval, dfn, test);
	pputc(prn, '\n');
    }
    pputc(prn, '\n');
}

static void save_chow_test (MODEL *pmod, const char *chowstr,
			    double test, double pval,
			    int dfn, int dfd,
			    gretlopt opt)
{
    int ttype = (opt & OPT_D)? GRETL_TEST_CHOWDUM : GRETL_TEST_CHOW;
    ModelTest *mt = model_test_new(ttype);

    if (mt != NULL) {
	if (dfd == 0) {
	    model_test_set_teststat(mt, GRETL_STAT_WALD_CHISQ);
	} else {
	    model_test_set_teststat(mt, GRETL_STAT_F);
	}
	model_test_set_param(mt, chowstr);
	model_test_set_value(mt, test);
	model_test_set_pvalue(mt, pval);
	model_test_set_dfn(mt, dfn);
	model_test_set_dfd(mt, dfd);
	maybe_add_test_to_model(pmod, mt);
    }
}

/*
 * real_chow_test:
 * @chowparm: sample breakpoint; or ID number of dummy
 * variable if given OPT_D; or ignored if given OPT_T.
 * @pmod: pointer to model to be tested.
 * @dset: dataset information.
 * @ci: either CHOW or QLRTEST.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_D included, do the Chow test based on a given dummy
 * variable; if OPT_T included, do the QLR test. If OPT_M,
 * just save the test on the model, don't set $test and
 * $pvalue.
 * @prn: gretl printing struct.
 *
 * Returns: 0 on successful completion, error code on error.
 */

static int real_chow_test (int chowparm, MODEL *pmod, DATASET *dset,
			   int ci, gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int *chowlist = NULL;
    int origv = dset->v;
    MODEL chow_mod;
    int QLR = (opt & OPT_T);
    int dumv = 0, split = 0, smax = 0;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    /* temporarily impose the sample that was in force when the
       original model was estimated */
    impose_model_smpl(pmod, dset);

    gretl_model_init(&chow_mod, dset);

    if (QLR) {
	/* "15 percent trimming": exactly how this should be
	   defined is debatable, but the following agrees
	   with the strucchange package for R, when the
	   Fstats() function is given "from=0.15, to=0.85".
	*/
	split = pmod->t1 + floor(0.15 * pmod->nobs) - 1;
	smax = pmod->t1 + floor(0.85 * pmod->nobs) - 1;
    } else if (opt & OPT_D) {
	/* Chow, using a predefined dummy */
	dumv = chowparm;
    } else {
	/* Chow, using break observation */
	smax = split = chowparm;
    }

    if (!err) {
	chowlist = make_chow_list(pmod, dset, split, dumv, ci,
				  opt, &err);
    }

    if (err) {
	goto bailout;
    }

    if (QLR) {
	/* Quandt likelihood ratio */
	int robust = (pmod->opt & OPT_R);
	gretlopt lsqopt = OPT_A;
	double test, testmax = 0.0;
	double *testvec = NULL;
	int *testlist = NULL;
	int dfn = 0, dfd = 0;
	int nextra = dset->v - origv;
	int tmax = 0;
	int i, t;

	set_effective_plot_ci(QLRTEST);
	if (gnuplot_graph_wanted(PLOT_REGULAR, opt, NULL)) {
	    testvec = malloc((smax - split + 1) * sizeof *testvec);
	}
	reset_effective_plot_ci();

	if (robust) {
	    lsqopt |= OPT_R;
	    testlist = gretl_list_diff_new(chowlist, pmod->list, 2);
	    if (testlist == NULL) {
		err = E_ALLOC;
	    }
	}

	for (t=split; t<=smax && !err; t++) {
	    chow_mod = lsq(chowlist, dset, OLS, lsqopt);
	    if (chow_mod.errcode) {
		err = chow_mod.errcode;
		errmsg(err, prn);
		break;
	    }
	    dfn = chow_mod.ncoeff - pmod->ncoeff;
	    if (robust) {
		test = wald_omit_chisq(testlist, &chow_mod);
#if 0
		fprintf(stderr, "%d %g\n", t+1, test);
#endif
	    } else {
		dfd = chow_mod.dfd;
		test = (pmod->ess - chow_mod.ess) * dfd / (chow_mod.ess * dfn);
#if 0
		fprintf(stderr, "%d %g\n", t+1, test*dfn);
#endif
	    }
	    if (test > testmax) {
		tmax = t;
		testmax = test;
	    }
	    if (testvec != NULL) {
		testvec[t - split] = test;
	    }
#if 0
	    fprintf(stderr, "split at t=%d: F(%d,%d)=%g (X2=%g)\n", t,
		    dfn, dfd, F, test*dfn);
	    fprintf(stderr, " pmod->ess = %g, chow_mod.ess = %g\n",
		    pmod->ess, chow_mod.ess);
#endif
	    clear_model(&chow_mod);
	    for (i=0; i<nextra; i++) {
		dset->Z[origv+i][t] = 0.0;
	    }
	}

	if (!err) {
	    char datestr[OBSLEN];
	    double pval, X2;

	    X2 = robust ? testmax : dfn * testmax;
	    pval = get_QLR_pval(X2, dfn, split, smax, pmod);
	    ntolabel(datestr, tmax, dset);

	    if (!(opt & OPT_Q)) {
		QLR_print_result(testmax, pval, datestr, dfn, dfd,
				 robust, prn);
	    }
	    if (!(opt & OPT_M)) {
		record_QLR_test_result(X2, pval, tmax + 1);
	    }
	    if (opt & (OPT_S | OPT_M)) {
		save_QLR_test(pmod, datestr, X2, pval, dfn);
	    }
	    if (testvec != NULL) {
		QLR_graph(testvec, split, smax, dset, dfn, robust);
	    }
	}

	free(testvec);
	free(testlist);
    } else {
	/* regular (or robust) Chow test */
	int robust = (pmod->opt & OPT_R);
	gretlopt lsqopt = OPT_A;
	int *testlist = NULL;

	if (robust) {
	    lsqopt |= OPT_R;
	}
	chow_mod = lsq(chowlist, dset, OLS, lsqopt);

	if (chow_mod.errcode) {
	    err = chow_mod.errcode;
	    errmsg(err, prn);
	} else if (chow_mod.ncoeff <= pmod->ncoeff) {
	    err = chow_mod.errcode = E_DATA;
	    errmsg(err, prn);
	} else {
	    int dfd = (robust)? 0 : chow_mod.dfd;
	    int dfn = chow_mod.ncoeff - pmod->ncoeff;
	    double test, pval = NADBL;

	    if (!(opt & OPT_Q)) {
		chow_mod.aux = AUX_CHOW;
		printmodel(&chow_mod, dset, OPT_NONE, prn);
	    }

	    if (robust) {
		testlist = gretl_list_diff_new(chow_mod.list, pmod->list, 2);
		if (testlist == NULL) {
		    err = E_ALLOC;
		    goto bailout;
		} else {
		    test = wald_omit_chisq(testlist, &chow_mod);
		}
		if (!na(test)) {
		    pval = chisq_cdf_comp(dfn, test);
		}
	    } else {
		test = (pmod->ess - chow_mod.ess) * dfd / (chow_mod.ess * dfn);
		if (!na(test)) {
		    pval = snedecor_cdf_comp(dfn, dfd, test);
		}
	    }

	    if (na(test)) {
		pputs(prn, _("Couldn't compute Chow test statistic\n"));
	    } else if (na(pval)) {
		pprintf(prn, _("Couldn't compute Chow test p-value (test = %g)\n"), test);
	    }

	    if (!na(test) && !na(pval)) {
		char chowstr[VNAMELEN];

		if (opt & OPT_Q) {
		    pputc(prn, '\n');
		}
		if (opt & OPT_D) {
		    strcpy(chowstr, dset->varname[chowparm]);
		    pprintf(prn, _("Chow test for structural difference with respect to %s"),
			    chowstr);
		} else {
		    ntolabel(chowstr, chowparm, dset);
		    pprintf(prn, _("Chow test for structural break at observation %s"),
			    chowstr);
		}
		pputc(prn, '\n');

		if (robust) {
		    pprintf(prn, "  %s(%d) = %g %s %.4f\n", _("Chi-square"),
			    dfn, test, _("with p-value"), pval);
		    pprintf(prn, "  %s: F(%d, %d) = %g %s %.4f\n\n", _("F-form"),
			    dfn, chow_mod.dfd, test / dfn, _("with p-value"),
			    snedecor_cdf_comp(dfn, chow_mod.dfd, test / dfn));
		} else {
		    pprintf(prn, "  F(%d, %d) = %g %s %.4f\n\n",
			    dfn, dfd, test, _("with p-value"), pval);
		}

		if (opt & OPT_S) {
		    save_chow_test(pmod, chowstr, test, pval, dfn, dfd, opt);
		}

		record_test_result(test, pval);
	    }
	}
	clear_model(&chow_mod);
	free(testlist);
    }

 bailout:

    /* clean up extra variables */
    dataset_drop_last_variables(dset, dset->v - origv);
    free(chowlist);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/**
 * chow_test:
 * @splitobs: the 0-based observation number at which to split
 * the sample.
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test)
 * using the sample break-point given by @splitobs
 * and prints the results to @prn. (The second portion of
 * the sample runs from observation @splitobs to the
 * end of the original sample range.)
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test (int splitobs, MODEL *pmod, DATASET *dset,
	       gretlopt opt, PRN *prn)
{
    int err = 0;

    if (splitobs <= 0 || splitobs >= dset->n) {
	gretl_errmsg_set(_("Invalid sample split for Chow test"));
	err = E_DATA;
    } else {
	err = real_chow_test(splitobs, pmod, dset, CHOW, opt, prn);
    }

    return err;
}

/**
 * chow_test_from_dummy:
 * @splitvar: the ID number of a dummy variable that should
 * be used to divide the sample.
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability (Chow test),
 * using the dummy variable specified by @splitvar to divide the
 * original sample range of @pmod into two portions.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test_from_dummy (int splitvar, MODEL *pmod, DATASET *dset,
			  gretlopt opt, PRN *prn)
{
    int err = 0;

    if (splitvar <= 0 || splitvar >= dset->v) {
	gretl_errmsg_set(_("Invalid sample split for Chow test"));
	err = E_DATA;
    } else {
	err = real_chow_test(splitvar, pmod, dset, CHOW,
			     opt | OPT_D, prn);
    }

    return err;
}

/**
 * QLR_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for structural stability via the
 * Quandt Likelihood Ratio test for a structural break at
 * an unknown point in the sample range.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int QLR_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    return real_chow_test(0, pmod, dset, QLRTEST,
			  opt | OPT_T, prn);
}

/* compute v'Mv, for symmetric M */

static double vprime_M_v (double *v, double *M, int n)
{
    int i, j, jmin, jmax, k;
    double xx, val = 0.0;

    k = jmin = 0;
    for (i=0; i<n; i++) {
	xx = 0.0;
	for (j=jmin; j<n; j++) {
	    xx += v[j] * M[k++];
	}
	val += xx * v[i];
	jmin++;
    }

    jmax = 1;
    for (i=1; i<n; i++) {
	k = i;
	xx = 0.0;
	for (j=0; j<jmax; j++) {
	    xx += v[j] * M[k];
	    k += n - j - 1;
	}
	val += xx * v[i];
	jmax++;
    }

    return val;
}

/* Compute scaled recursive residuals and cumulate \bar{w} via a
   sequence of OLS regressions. See William Greene, Econometric
   Analysis 5e, pp. 135-136.
*/

static int cusum_compute (MODEL *pmod, double *cresid, int T, int k,
			  double *wbar, DATASET *dset)
{
    MODEL cmod;
    gretlopt opt = OPT_X | OPT_A;
    double yhat, den, *xt;
    int n = T - k;
    int i, j, t;
    int err = 0;

    /* vector of regressors at t */
    xt = malloc(k * sizeof *xt);
    if (xt == NULL) {
	return E_ALLOC;
    }

    /* set minimal sample based on model to be tested */
    dset->t1 = pmod->t1;
    dset->t2 = pmod->t1 + k - 1;

    for (j=0; j<n && !err; j++) {
	cmod = lsq(pmod->list, dset, OLS, opt);
	if (cmod.errcode) {
	    err = cmod.errcode;
	} else {
	    /* compute one step ahead prediction error: \hat{y}_t
	       based on estimates through t - 1
	    */
	    t = dset->t2 + 1;
	    yhat = 0.0;
	    for (i=0; i<cmod.ncoeff; i++) {
		xt[i] = dset->Z[cmod.list[i+2]][t];
		yhat += cmod.coeff[i] * xt[i];
	    }
	    cresid[j] = dset->Z[pmod->list[1]][t] - yhat;
	    cmod.ci = CUSUM;
	    err = makevcv(&cmod, 1.0); /* (X'X)^{-1} */
	    if (!err) {
		/* compute scaled residual: divide by
		   \sqrt{1 + x_t'(X_{t-1}'X_{t-1})^{-1} x_t}
		*/
		den = vprime_M_v(xt, cmod.vcv, cmod.ncoeff);
		cresid[j] /= sqrt(1.0 + den);
		*wbar += cresid[j];
		/* extend the sample by one observation */
		dset->t2 += 1;
	    }
	}
	clear_model(&cmod);
    }

    free(xt);

    return err;
}

#define okfreq(p) (p == 1 || p == 4 || p == 12 || p == 24 || p == 52)

static int cusum_do_graph (double a, double b, const double *W,
			   int t1, int k, int m,
			   DATASET *dset, gretlopt opt)
{
    FILE *fp = NULL;
    GptFlags flags = 0;
    const double *obs = NULL;
    double frac = 1.0;
    double x0 = 0.0;
    int j, t, err = 0;

    if (dataset_is_time_series(dset)) {
	flags = GPT_TS | GPT_LETTERBOX;
    }

    fp = open_plot_input_file(PLOT_CUSUM, flags, &err);
    if (err) {
	return err;
    }

    if (dataset_is_time_series(dset) && okfreq(dset->pd)) {
	b *= dset->pd;
	frac /= dset->pd;
        obs = gretl_plotx(dset, OPT_NONE);
	if (obs != NULL) {
	    x0 = obs[t1 + k];
	}
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xlabel '%s'\n", _("Observation"));
    fputs("set nokey\n", fp);

    if (opt & OPT_R) {
	fprintf(fp, "set title '%s'\n",
		/* xgettext:no-c-format */
		_("CUSUMSQ plot with 95% confidence band"));
	fprintf(fp, "plot \\\n%g*(x-%g) title '' w dots lt 2, \\\n", b, x0 - frac);
	fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n", -a, b, x0 - frac);
	fprintf(fp, "%g+%g*(x-%g) title '' w lines lt 2, \\\n", a, b, x0 - frac);
    } else {
	fputs("set xzeroaxis\n", fp);
	fprintf(fp, "set title '%s'\n",
		/* xgettext:no-c-format */
		_("CUSUM plot with 95% confidence band"));
	fprintf(fp, "plot \\\n%g+%g*(x-%g) title '' w lines lt 2, \\\n", a, b, x0);
	fprintf(fp, "%g-%g*(x-%g) title '' w lines lt 2, \\\n", -a, b, x0);
    }

    fputs("'-' using 1:2 w linespoints lt 1\n", fp);

    for (j=0; j<m; j++) {
	t = t1 + k + j;
	if (obs != NULL) {
	    fprintf(fp, "%g %g\n", obs[t], W[j]);
	} else {
	    fprintf(fp, "%d %g\n", t, W[j]);
	}
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static void cusum_harvey_collier (double wbar, double sigma, int m,
				  MODEL *pmod, gretlopt opt,
				  PRN *prn)
{
    double hct, pval;

    hct = (sqrt((double) m) * wbar) / sigma;
    pval = student_pvalue_2(m - 1, hct);
    pprintf(prn, _("\nHarvey-Collier t(%d) = %g with p-value %.4g\n\n"),
	    m - 1, hct, pval);

    if (opt & OPT_S) {
	ModelTest *test = model_test_new(GRETL_TEST_CUSUM);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_HARVEY_COLLIER);
	    model_test_set_dfn(test, m - 1);
	    model_test_set_value(test, hct);
	    model_test_set_pvalue(test, pval);
	    maybe_add_test_to_model(pmod, test);
	}
    }

    record_test_result(hct, pval);
}

/**
 * cusum_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include %OPT_S, save results of test to model.
 * @prn: gretl printing struct.
 *
 * Tests the given model for parameter stability via the CUSUM test,
 * or if @opt includes %OPT_R, via the CUSUMSQ test; %OPT_Q makes
 * the test quiet; %OPT_U governs the associated plot, if wanted.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int cusum_test (MODEL *pmod, DATASET *dset,
		gretlopt opt, PRN *prn)
{
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    char cumdate[OBSLEN];
    double wbar = 0.0;
    double *cresid = NULL, *W = NULL;
    int quiet = opt & OPT_Q;
    int m, i, j;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    if (model_has_missing_obs(pmod)) {
	return E_MISSDATA;
    }

    /* the number of forecasts */
    m = T - k;

    cresid = malloc(m * sizeof *cresid);
    W = malloc(m * sizeof *W);

    if (cresid == NULL || W == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = cusum_compute(pmod, cresid, T, k, &wbar, dset);
	if (err) {
	    errmsg(err, prn);
	}
    }

    if (!err) {
	double a, b, den = 0.0, sigma = 0.0;
	int sig;

	if (opt & OPT_R) {
	    double a1, a2, a3;
	    double n = 0.5 * m - 1;

	    pprintf(prn, "\n%s\n\n", _("CUSUMSQ test for stability of parameters"));

	    for (j=0; j<m; j++) {
		den += cresid[j] * cresid[j];
	    }

	    /* see Edgerton and Wells, Oxford Bulletin of Economics and
	       Statistics, 56, 1994, pp. 355-365 */
	    a1 = 1.358015 / sqrt(n);
	    a2 = -0.6701218 / n;
	    a3 = -0.8858694 / pow(n, 1.5);

	    /* 0.5 * total height of band */
	    a = a1 + a2 + a3;
	    /* slope of expectation wrt time */
	    b = 1.0 / m;
	    if (!quiet) {
		pputs(prn, _("Cumulated sum of squared residuals"));
	    }
	} else {
	    wbar /= T - k;
	    pprintf(prn, "\n%s\n\n", _("CUSUM test for stability of parameters"));
	    pprintf(prn, _("mean of scaled residuals = %g\n"), wbar);

	    for (j=0; j<m; j++) {
		sigma += (cresid[j] - wbar) * (cresid[j] - wbar);
	    }
	    sigma /= T - k - 1;
	    sigma = sqrt(sigma);
	    pprintf(prn, _("sigmahat                 = %g\n\n"), sigma);

	    /* height of confidence band for first prediction */
	    a = 0.948 * sqrt((double) m);
	    /* slope of confidence band limit wrt time */
	    b = 2.0 * a / m;
	    if (!quiet) {
		pputs(prn, _("Cumulated sum of scaled residuals"));
	    }
	}

	pputc(prn, '\n');
	pputs(prn, /* xgettext:no-c-format */
	      _("('*' indicates a value outside of 95% confidence band)"));
	pputs(prn, "\n\n");

	for (j=0; j<m; j++) {
	    W[j] = 0.0;
	    if (opt & OPT_R) {
		for (i=0; i<=j; i++) {
		    W[j] += cresid[i] * cresid[i] / den;
		}
		sig = fabs(W[j] - (j+1) / (double) m) > a;
	    } else {
		for (i=0; i<=j; i++) {
		    W[j] += cresid[i];
		}
		W[j] /= sigma;
		sig = fabs(W[j]) > a + j * b;
	    }
	    if (!quiet) {
		ntolabel(cumdate, pmod->t1 + k + j, dset);
		pprintf(prn, " %s %9.3f %s\n", cumdate, W[j], sig? "*" : "");
	    }
	}

	if (!(opt & OPT_R)) {
	    cusum_harvey_collier(wbar, sigma, m, pmod, opt, prn);
	}

	/* plot with 95% confidence band if wanted */
	if (gnuplot_graph_wanted(PLOT_CUSUM, opt, NULL)) {
	    err = cusum_do_graph(a, b, W, pmod->t1, k, m, dset, opt);
	}
    }

    /* restore original sample range */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(cresid);
    free(W);

    return err;
}

/**
 * comfac_test:
 * @pmod: pointer to original model.
 * @dset: dataset struct.
 * @opt: if contains %OPT_S, save test results to model.
 * @prn: gretl printing struct.
 *
 * If @pmod was estimated via an AR(1) estimator, run an
 * auxiliary regression to test the implied common-factor
 * restriction.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int comfac_test (MODEL *pmod, DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    MODEL cmod;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int v = dset->v;
    int *biglist = NULL;
    int clearit = 0;
    int nadd, i, k, t;
    int err;

    if (pmod->ci != AR1 || (pmod->opt & OPT_P)) {
	/* can't handle Prais-Winsten? */
	return E_NOTIMP;
    }

    /* check for changes in original list members */
    err = list_members_replaced(pmod, dset);
    if (err) {
	return err;
    }

    biglist = gretl_list_copy(pmod->list);
    if (biglist == NULL) {
	return E_ALLOC;
    }

    nadd = 1 + pmod->ncoeff - pmod->ifc;

    err = dataset_add_series(dset, nadd);
    if (err) {
	free(biglist);
	return err;
    }

    /* add lags of the dependent variable and all regressors: some of
       these may be redundant but we'll just let the redundant terms be
       eliminated automatically via gretl's collinearity checking.
    */

    k = v;
    for (i=1; i<=pmod->list[0]; i++) {
	int src, lag, parent;

	src = pmod->list[i];
	if (src == 0) {
	    continue;
	}
	for (t=0; t<dset->n; t++) {
	    if (t == 0 || na(dset->Z[src][t-1])) {
		dset->Z[k][t] = NADBL;
	    } else {
		dset->Z[k][t] = dset->Z[src][t-1];
	    }
	}
	biglist = gretl_list_append_term(&biglist, k);
	if (biglist == NULL) {
	    err = E_ALLOC;
	    break;
	}
	lag = series_get_lag(dset, src);
	parent = series_get_parent_id(dset, src);
	if (lag > 0 && parent > 0) {
	    char tmp[16];

	    sprintf(tmp, "_%d", lag + 1);
	    strcpy(dset->varname[k], dset->varname[parent]);
	    gretl_trunc(dset->varname[k], 15 - strlen(tmp));
	    strcat(dset->varname[k], tmp);
	} else {
	    strcpy(dset->varname[k], dset->varname[src]);
	    gretl_trunc(dset->varname[k], 13);
	    strcat(dset->varname[k], "_1");
	}
	k++;
    }

    if (!err) {
	/* re-impose the sample that was in force when the original model
	   was estimated */
	impose_model_smpl(pmod, dset);
	cmod = lsq(biglist, dset, OLS, OPT_A);
	clearit = 1;
	err = cmod.errcode;
    }

    if (!err) {
	if (cmod.nobs != pmod->nobs || cmod.ess > pmod->ess || cmod.dfd >= pmod->dfd) {
	    /* something has gone wrong */
	    err = E_DATA;
	}
    }

    if (!err) {
	/* construct an F-test based on the SSR from the original
	   AR(1) model and the SSR from the unrestricted model, cmod.
	*/
	int dfd = cmod.dfd;
	int dfn = pmod->dfd - dfd - 1; /* account for rho */
	double SSRr = pmod->ess;
	double SSRu = cmod.ess;
	double Ftest = ((SSRr - SSRu)/dfn) / (SSRu/dfd);
	double pval = snedecor_cdf_comp(dfn, dfd, Ftest);

	if (!(opt & OPT_I)) {
	    if (!(opt & OPT_Q)) {
		cmod.aux = AUX_COMFAC;
		printmodel(&cmod, dset, OPT_S, prn);
		pputc(prn, '\n');
	    }

	    pputs(prn, _("Test of common factor restriction"));
	    pputs(prn, "\n\n");

	    pprintf(prn, "  %s: %s(%d, %d) = %g, ", _("Test statistic"),
		    "F", dfn, dfd, Ftest);
	    pprintf(prn, _("with p-value = %g\n"), pval);
	    pputc(prn, '\n');
	}

	if (opt & OPT_S) {
	    ModelTest *test;

	    test = model_test_new(GRETL_TEST_COMFAC);
	    if (test != NULL) {
		model_test_set_teststat(test, GRETL_STAT_F);
		model_test_set_dfn(test, dfn);
		model_test_set_dfd(test, dfd);
		model_test_set_value(test, Ftest);
		model_test_set_pvalue(test, pval);
		maybe_add_test_to_model(pmod, test);
	    }
	}

	record_test_result(Ftest, pval);
    }

    if (clearit) {
	clear_model(&cmod);
    }

    /* delete the added variables and restore the original
       sample range */

    dataset_drop_last_variables(dset, nadd);
    free(biglist);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}

/**
 * panel_specification_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Tests the given pooled model for fixed and random effects.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int panel_specification_test (MODEL *pmod, DATASET *dset,
			      gretlopt opt, PRN *prn)
{
    if (pmod->ci != OLS || !dataset_is_panel(dset)) {
	pputs(prn, _("This test is only relevant for pooled models\n"));
	return 1;
    }

    if (pmod->ifc == 0) {
	pputs(prn, _("This test requires that the model contains a constant\n"));
	return 1;
    }

    if (!balanced_panel(dset)) { /* ?? */
	pputs(prn, _("Sorry, can't do this test on an unbalanced panel.\n"
		"You need to have the same number of observations\n"
		"for each cross-sectional unit"));
	return 1;
    } else {
	panel_diagnostics(pmod, dset, opt, prn);
    }

    return 0;
}

/* wrapper for backward compatibility */

int panel_hausman_test (MODEL *pmod, DATASET *dset,
			gretlopt opt, PRN *prn)
{
    return panel_specification_test(pmod, dset, opt, prn);
}

/**
 * add_leverage_values_to_dataset:
 * @dset: dataset struct.
 * @m: matrix containing leverage values.
 * @opt: OPT_O = overwrite series on save.
 * @flags: may include SAVE_LEVERAGE, SAVE_INFLUENCE,
 * and/or SAVE_DFFITS.
 *
 * Adds to the working dataset one or more series calculated by
 * the gretl test for leverage/influence of data points.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int add_leverage_values_to_dataset (DATASET *dset, gretl_matrix *m,
				    gretlopt opt, int flags)
{
    int overwrite = (opt & OPT_O);
    int t1, t2;
    int addvars = 0;

    if (flags & SAVE_LEVERAGE) addvars++;
    if (flags & SAVE_INFLUENCE) addvars++;
    if (flags & SAVE_DFFITS) addvars++;

    if (addvars == 0) {
	return 0;
    }

    if (dataset_add_series(dset, addvars)) {
	return E_ALLOC;
    }

    t1 = gretl_matrix_get_t1(m);
    t2 = t1 + gretl_matrix_rows(m);

    /* add leverage? */
    if (flags & SAVE_LEVERAGE) {
	int t, v = dset->v - addvars;
	int j = 0;

	for (t=0; t<dset->n; t++) {
	    if (t < t1 || t >= t2) {
		dset->Z[v][t] = NADBL;
	    } else {
		dset->Z[v][t] = gretl_matrix_get(m, j++, 0);
	    }
	}
	strcpy(dset->varname[v], "lever");
	if (!overwrite) {
	    make_varname_unique(dset->varname[v], v, dset);
	}
	series_set_label(dset, v, "leverage values");
    }

    /* add influence? */
    if (flags & SAVE_INFLUENCE) {
	int t, v = dset->v - (addvars - 1);
	int j = 0;

	for (t=0; t<dset->n; t++) {
	    if (t < t1 || t >= t2) {
		dset->Z[v][t] = NADBL;
	    } else {
		dset->Z[v][t] = gretl_matrix_get(m, j++, 1);
	    }
	}
	strcpy(dset->varname[v], "influ");
	if (!overwrite) {
	    make_varname_unique(dset->varname[v], v, dset);
	}
	series_set_label(dset, v, "influence values");
    }

    /* add DFFITS? */
    if (flags & SAVE_DFFITS) {
	int t, v = dset->v - (addvars - 2);
	int j = 0;

	for (t=0; t<dset->n; t++) {
	    double s, h;

	    if (t < t1 || t >= t2) {
		dset->Z[v][t] = NADBL;
	    } else {
		/* s = studentized residuals */
		h = gretl_matrix_get(m, j, 0);
		s = gretl_matrix_get(m, j, 2);
		if (na(h) || na(s)) {
		    dset->Z[v][t] = NADBL;
		} else {
		    dset->Z[v][t] = s * sqrt(h / (1.0 - h));
		}
		j++;
	    }
	}
	strcpy(dset->varname[v], "dffits");
	if (!overwrite) {
	    make_varname_unique(dset->varname[v], v, dset);
	}
	series_set_label(dset, v, "DFFITS values");
    }

    set_dataset_is_changed(dset, 1);

    return 0;
}

/**
 * leverage_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if OPT_S, add calculated series to data set; operate
 * silently if OPT_Q. If includes OPT_O, do not adjust save
 * names (overwrite any existing series of the same names).
 * @prn: gretl printing struct.
 *
 * Tests the data used in the given model for points with
 * high leverage and influence on the estimates.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int leverage_test (MODEL *pmod, DATASET *dset,
		   gretlopt opt, PRN *prn)
{
    gretl_matrix *(*model_leverage) (const MODEL *, const DATASET *,
				     gretlopt, PRN *, int *);
    gretl_matrix *m;
    int err = 0;

    if (pmod->ci != OLS) {
	return E_OLSONLY;
    }

    model_leverage = get_plugin_function("model_leverage");
    if (model_leverage == NULL) {
	return 1;
    }

    m = (*model_leverage)(pmod, dset, opt, prn, &err);

    if (!err) {
	/* set the $results accessor */
	set_last_result_data(m, GRETL_TYPE_MATRIX);
	if (opt & OPT_S) {
	    /* and respond to the --save option */
	    err = add_leverage_values_to_dataset(dset, m, opt,
						 SAVE_LEVERAGE |
						 SAVE_INFLUENCE |
						 SAVE_DFFITS);
	}
    }

    return err;
}

/**
 * vif_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: may contain OPT_Q for quiet operation.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Variance Inflation Factors for
 * the independent variables in the given model.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int vif_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    int (*compute_vifs) (MODEL *, DATASET *, gretlopt, PRN *);

    gretl_error_clear();

    compute_vifs = get_plugin_function("compute_vifs");
    if (compute_vifs == NULL) {
	return 1;
    }

    return (*compute_vifs)(pmod, dset, opt, prn);
}

/**
 * bkw_test:
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: may contain OPT_Q for quiet operation.
 * @prn: gretl printing struct.
 *
 * Calculates and displays the Belsley-Kuh-Welsch collinearity
 * diagnostics for @pmod.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int bkw_test (MODEL *pmod, DATASET *dset,
	      gretlopt opt, PRN *prn)
{
    int (*compute_bkw) (MODEL *, DATASET *, gretlopt, PRN *);

    gretl_error_clear();

    compute_bkw = get_plugin_function("compute_bkw");
    if (compute_bkw == NULL) {
	return 1;
    }

    return (*compute_bkw)(pmod, dset, opt, prn);
}
