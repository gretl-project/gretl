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

/* Miscellaneous functions to bridge between gretl commands and
   the corresponding libgretl functions. This glue allows a cleaner
   interface for the latter, hiding some parsing of strings
   coming from the command line.
*/

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "gretl_panel.h"
#include "usermat.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "boxplots.h"
#include "gretl_string_table.h"
#include "gretl_drivers.h"

/*
 * model_test_driver:
 * @order: lag order for --autocorr and --arch.
 * @dset: dataset struct.
 * @opt: controls which test(s) will be performed; OPT_Q
 * gives less verbose results, OPT_I gives silent operation.
 * @prn: gretl printing struct.
 *
 * Performs some subset of gretl's "modtest" tests on the
 * model last estimated, and prints the results to @prn.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int model_test_driver (int order, DATASET *dset,
		       gretlopt opt, PRN *prn)
{
    GretlObjType type;
    gretlopt testopt = OPT_NONE;
    void *ptr;
    int k = 0;
    int err = 0;

    if (opt == OPT_NONE || opt == OPT_Q || opt == OPT_I) {
	/* note: OPT_Q and OPT_I are just quiet and silent respectively */
	pprintf(prn, _("modtest: no options selected\n"));
	return 0;
    }

    err = incompatible_options(opt, OPT_A | OPT_H | OPT_L | OPT_S |
			       OPT_N | OPT_P | OPT_W | OPT_X | OPT_D);
    if (err) {
	return err;
    }

    ptr = get_last_model(&type);
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN && exact_fit_check(ptr, prn)) {
	return 0;
    }

    if (opt & (OPT_A | OPT_H)) {
	/* autocorrelation and arch: lag order */
	k = order > 0 ? order : dset->pd;
    }

    /* transcribe the quietness flags */
    if (opt & OPT_I) {
	testopt = OPT_I | OPT_Q;
    } else if (opt & OPT_Q) {
	testopt = OPT_Q;
    }

    /* non-linearity (squares) */
    if (!err && (opt & OPT_S)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, dset, AUX_SQ,
				    testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* non-linearity (logs) */
    if (!err && (opt & OPT_L)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, dset, AUX_LOG,
				    testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* heteroskedasticity (White or Breusch-Pagan) */
    if (!err && (opt & (OPT_W | OPT_X | OPT_B))) {
	if (type == GRETL_OBJ_EQN) {
	    transcribe_option_flags(&testopt, opt, OPT_B | OPT_X);
	    if ((opt & OPT_B) && (opt & OPT_R)) {
		testopt |= OPT_R;
	    }
	    err = whites_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* autocorrelation */
    if (!err && (opt & OPT_A)) {
	if (type == GRETL_OBJ_EQN) {
	    err = autocorr_test(ptr, k, dset, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    if (opt & OPT_U) {
		/* --univariate */
		testopt |= OPT_U;
	    }
	    err = gretl_VAR_autocorrelation_test(ptr, k, dset,
						 testopt, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_autocorrelation_test(ptr, k, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* ARCH */
    if (!err && (opt & OPT_H)) {
	if (type == GRETL_OBJ_EQN) {
	    err = arch_test(ptr, k, dset, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    if (opt & OPT_U) {
		/* --univariate */
		testopt |= OPT_U;
	    }
	    err = gretl_VAR_arch_test(ptr, k, dset,
				      testopt, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_arch_test(ptr, k, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* normality of residual */
    if (!err && (opt & OPT_N)) {
	return last_model_test_uhat(dset, testopt, prn);
    }

    /* groupwise heteroskedasticity */
    if (!err && (opt & OPT_P)) {
	if (type == GRETL_OBJ_EQN) {
	    err = groupwise_hetero_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* common factor restriction */
    if (!err && (opt & OPT_C)) {
	if (type == GRETL_OBJ_EQN) {
	    err = comfac_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    /* cross-sectional dependence */
    if (!err && (opt & OPT_D)) {
	if (type == GRETL_OBJ_EQN) {
	    err = panel_xdepend_test(ptr, dset, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
	return err;
    }

    return err;
}

static int get_chow_dummy (const char *s, const DATASET *dset,
			   int *err)
{
    int v = current_series_index(dset, s);

    if (v < 0) {
	*err = E_UNKVAR;
    } else if (!gretl_isdummy(dset->t1, dset->t2, dset->Z[v])) {
	*err = E_DATA;
    }

    return v;
}

/*
 * chow_test_driver:
 * @param: parameter (observation or name of dummy)
 * @pmod: pointer to model to be tested.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_D included, do the Chow test based on a given dummy
 * variable.
 * @prn: gretl printing struct.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test_driver (const char *param, MODEL *pmod, DATASET *dset,
		      gretlopt opt, PRN *prn)
{
    int chowparm = 0;
    int err = 0;

    if (param == NULL || *param == '\0') {
	return E_DATA;
    }

    if (opt & OPT_D) {
	chowparm = get_chow_dummy(param, dset, &err);
    } else {
	chowparm = dateton(param, dset);
    }

    if (!err) {
	if (opt & OPT_D) {
	    err = chow_test_from_dummy(chowparm, pmod, dset, opt, prn);
	} else {
	    err = chow_test(chowparm, pmod, dset, opt, prn);
	}
    }

    return err;
}

/* The @param here may contain a scalar or a matrix: in either case,
   convert to a list of lag-orders before handing off to the real
   Levin-Lin-Chu code.
*/

int llc_test_driver (const char *param, const int *list,
		     DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_matrix *m = NULL;
    int *plist = NULL;
    int p0 = -1;
    int err = 0;

    if (param == NULL) {
	err = E_DATA;
    } else if (*param == '{') {
	m = generate_matrix(param, dset, &err);
	if (!err) {
	    plist = gretl_auxlist_from_vector(m, &err);
	}
	gretl_matrix_free(m);
    } else if (get_matrix_by_name(param)) {
	m = get_matrix_by_name(param);
	plist = gretl_auxlist_from_vector(m, &err);
    } else if (integer_string(param)) {
	p0 = atoi(param);
    } else if (gretl_is_scalar(param)) {
	p0 = gretl_scalar_get_value(param, NULL);
    } else {
	err = E_DATA;
    }

    if (!err) {
	if (plist != NULL) {
	    err = levin_lin_test(list[1], plist, dset, opt, prn);
	    free(plist);
	} else {
	    int tmplist[2] = {1, p0};

	    err = levin_lin_test(list[1], tmplist, dset, opt, prn);
	}
    }

    return err;
}

static void bds_print (const gretl_matrix *m,
		       const char *vname,
		       int order, double eps,
		       int c1, int boot,
		       double *detail, PRN *prn)
{
    double z, pv;
    int i;

    /* header */
    pputc(prn, '\n');
    pprintf(prn, _("BDS test for %s, maximum order %d"), vname, order);
    pputc(prn, '\n');
    pputs(prn, _("H0: the series is linear/IID"));
    pputc(prn, '\n');
    if (boot > 0) {
	pputs(prn, _("Bootstrapped p-values in []"));
    } else {
	pputs(prn, _("Asymptotic p-values in []"));
    }
    pputs(prn, "\n\n");

    /* test statistics */
    for (i=0; i<order-1; i++) {
	z = gretl_matrix_get(m, 0, i);
	pv = gretl_matrix_get(m, 1, i);
	pputs(prn, "  ");
	pprintf(prn, _("test order %d: z = %.3f [%.3f]"), i+2, z, pv);
	pputc(prn, '\n');
    }
    pputc(prn, '\n');

    /* trailer */
    if (c1) {
	pputs(prn, _("Distance criterion based on first-order correlation"));
    } else {
	pprintf(prn, _("Distance criterion based on sd(%s)"), vname);
    }
    pputc(prn, '\n');
    pprintf(prn, _("eps = %g, first-order correlation %.3f"),
	    detail[0], detail[1]);
    pputs(prn, "\n\n");
}

static int get_vector_x (const double **px, int *n,
			 const char **pvname)
{
    const char *mname = get_optval_string(BDS, OPT_X);
    int err = 0;

    if (mname != NULL) {
	gretl_matrix *m = get_matrix_by_name(mname);

	if (gretl_is_null_matrix(m)) {
	    err = E_INVARG;
	} else {
	    *n = gretl_vector_get_length(m);
	    if (*n == 0) {
		err = E_INVARG;
	    } else {
		*px = m->val;
		*pvname = mname;
	    }
	}
    } else {
	err = E_INVARG;
    }

    return err;
}

int bds_test_driver (int order, int *list, DATASET *dset,
		     gretlopt opt, PRN *prn)
{
    gretl_matrix *res = NULL;
    const double *x = NULL;
    const char *vname = NULL;
    double detail[2] = {0};
    double eps = 0.7;
    int t1 = dset->t1;
    int t2 = dset->t2;
    int boot = -1;
    int ctarget = 1;
    int n, v = 0;
    int err = 0;

    if (list == NULL) {
	err = get_vector_x(&x, &n, &vname);
	if (err) {
	    return err;
	}
    } else {
	v = list[1];
	x = dset->Z[v];
	vname = dset->varname[v];
    }

    if (order < 2) {
	err = E_INVARG;
    } else {
	err = series_adjust_sample(x, &t1, &t2);
    }

    if (!err) {
	err = incompatible_options(opt, OPT_S | OPT_C);
    }

    /* Note: @ctarget = 1 means that the value of @eps we pass to the
       bdstest() function is actually the target first-order
       correlation (with a default value of 0.7 that can be adjusted
       via the --corr1 option to the command). The alternative is that
       the user can specify @eps as a multiple of the std dev of the
       input series, in which case we set @ctarget to 0.
    */

    if (!err) {
	if (opt & OPT_S) {
	    /* eps as multiple of std dev of @x */
	    eps = get_optval_double(BDS, OPT_S, &err);
	    if (!err && eps <= 0) {
		err = E_INVARG;
	    }
	    ctarget = 0;
	} else if (opt & OPT_C) {
	    /* eps as target first-order correlation */
	    eps = get_optval_double(BDS, OPT_C, &err);
	    if (!err && (eps < 0.1 || eps > 0.9)) {
		err = E_INVARG;
	    }
	}
    }

    if (!err && (opt & OPT_B)) {
	boot = get_optval_int(BDS, OPT_B, &err);
	if (boot < 0) {
	    err = E_INVARG;
	}
    }

    if (!err) {
	gretl_matrix *(*bdstest) (const double *, int, int, double,
				  int, int, double *, int *);

	bdstest = get_plugin_function("bdstest");
	if (bdstest == NULL) {
	    err = E_FOPEN;
	} else {
	    n = t2 - t1 + 1;
	    if (boot < 0) {
		/* auto selection */
		boot = n < 600;
	    }
	    res = bdstest(x + t1, n, order, eps, ctarget, boot, detail, &err);
	}
    }

    if (res != NULL) {
	if (!(opt & OPT_Q)) {
	    bds_print(res, vname, order, eps, ctarget, boot, detail, prn);
	}
	set_last_result_data(res, GRETL_TYPE_MATRIX);
    }

    return err;
}

/* parse the tau vector out of @param before calling the
   "real" quantreg function
*/

MODEL quantreg_driver (const char *param, const int *list,
		       DATASET *dset, gretlopt opt, PRN *prn)
{
    gretl_vector *tau;
    MODEL mod;
    int err = 0;

    tau = generate_matrix(param, dset, &err);

    if (!err && gretl_vector_get_length(tau) == 0) {
	err = E_DATA;
    }

    if (err) {
	gretl_model_init(&mod, dset);
	mod.errcode = err;
    } else {
	mod = quantreg(tau, list, dset, opt, prn);
    }

    gretl_matrix_free(tau);

    return mod;
}

static int handle_binary_strval (DATASET *dset, int v)
{
    series_table *st = series_get_string_table(dset, v);
    int t, ns = series_table_get_n_strings(st);

    if (ns == 2) {
	/* convert temporarily to 0/1 */
	for (t=dset->t1; t<=dset->t2; t++) {
	    dset->Z[v][t] -= 1.0;
	}
	return 1;
    }

    return 0;
}

static void restore_binary_strval (MODEL *pmod, DATASET *dset, int v)
{
    const char *yname = dset->varname[v];
    const char *valstr;
    int n, t;

    /* restore the original numeric codes */
    for (t=dset->t1; t<=dset->t2; t++) {
	dset->Z[v][t] += 1.0;
    }

    /* and construct a suitable depvar name */
    valstr = series_get_string_for_value(dset, v, 2.0);
    n = strlen(yname) + strlen(valstr) + 5;
    pmod->depvar = calloc(n, 1);
    sprintf(pmod->depvar, "%s==\"%s\"", yname, valstr);
}

/* wrapper for the various sorts of logit and probit models
   that gretl supports
*/

MODEL logit_probit (int *list, DATASET *dset, int ci,
		    gretlopt opt, PRN *prn)
{
    int yv = list[1];
    int restore = 0;
    MODEL ret;

    if (ci == LOGIT && (opt & OPT_M)) {
	return multinomial_logit(list, dset, opt, prn);
    } else if (ci == PROBIT && (opt & OPT_E)) {
	return reprobit_model(list, dset, opt, prn);
    }

    if (is_string_valued(dset, yv)) {
	restore = handle_binary_strval(dset, yv);
    }

    if (gretl_isdummy(dset->t1, dset->t2, dset->Z[yv])) {
        ret = binary_model(ci, list, dset, opt, prn);
    } else {
        ret = ordered_estimate(ci, list, dset, opt, prn);
    }

    if (restore) {
	restore_binary_strval(&ret, dset, yv);
    }

    return ret;
}

/* parse out optional "ymax=..." parameter before calling the real
   logistic model function
*/

MODEL logistic_driver (const int *list, DATASET *dset,
		       gretlopt opt)
{
    double lmax = NADBL;

    if (opt & OPT_M) {
	int err = 0;

	lmax = get_optval_double(LOGISTIC, OPT_M, &err);
	if (err) {
	    MODEL mdl;

	    gretl_model_init(&mdl, dset);
	    mdl.errcode = err;
	    return mdl;
	}
    }

    return logistic_model(list, lmax, dset, opt);
}

/* assemble the left and right limits for tobit using gretl's
   option apparatus before calling the real tobit function
*/

MODEL tobit_driver (const int *list, DATASET *dset,
		    gretlopt opt, PRN *prn)
{
    MODEL model;
    double llim = -1.0e300;
    double rlim = NADBL;
    int err = 0;

    if (opt & OPT_L) {
	/* we should have an explicit lower limit */
	llim = get_optval_double(TOBIT, OPT_L, &err);
	if (!err && na(llim)) {
	    err = E_INVARG;
	}
    }

    if (!err && (opt & OPT_M)) {
	/* we should have an explicit upper limit */
	rlim = get_optval_double(TOBIT, OPT_M, &err);
	if (!err && (na(rlim) || rlim <= llim)) {
	    err = E_INVARG;
	}
    }

    if (err) {
	gretl_model_init(&model, dset);
	model.errcode = err;
	return model;
    }

    if (!(opt & (OPT_L | OPT_M))) {
	/* the default: left-censoring at zero */
	llim = 0;
    }

    return tobit_model(list, llim, rlim, dset, opt, prn);
}

static gretl_array *strings_array_from_string (const char *s,
					       int n, int *err)
{
    gretl_array *names = NULL;
    const char *sep = ",";
    char *tmp;
    int i;

    if (s == NULL) {
	*err = E_DATA;
	return NULL;
    }

    /* copy the incoming string @s before applying strtok */
    tmp = gretl_strdup(s);
    if (tmp == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    names = gretl_array_new(GRETL_TYPE_STRINGS, n, err);
    if (*err) {
	free(tmp);
	return NULL;
    }

    if (strchr(s, ',') == NULL) {
	sep = " ";
    }

    for (i=0; i<n && !*err; i++) {
	char *name = strtok((i == 0)? tmp : NULL, sep);

	if (name == NULL || *name == '\0') {
	    gretl_errmsg_sprintf(_("modprint: expected %d names"), n);
	    *err = E_DATA;
	} else {
	    while (isspace(*name)) {
		name++;
	    }
	    gretl_array_set_element(names, i, name,
				    GRETL_TYPE_STRING, 1);
	}
    }

    free(tmp);

    if (*err) {
	gretl_array_destroy(names);
	names = NULL;
    }

    return names;
}

/*
 * do_modprint:
 * @line: command line.
 * @opt: may contain %OPT_O for specifying output, and if
 * TeX output is called for then %OPT_C calls for
 * a complete LaTeX document.
 * @prn: gretl printer.
 *
 * Prints to @prn the coefficient table and optional additional statistics
 * for a model estimated "by hand". Mainly useful for user-written functions.
 *
 * The string @line must contain, in order: (1) the name of a k x 2 matrix
 * containing k coefficients and k associated standard errors and (2) the
 * name of a string variable containing at least k comma- or space-
 * separated names for the coefficients (or a string literal on that
 * pattern).
 *
 * Optionally, @line may contain a third element, the name of a vector
 * containing p additional statistics.  In that case element (2) should
 * contain k + p names, the additional p names to be associated with the
 * additional statistics.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int do_modprint (const char *mname, const char *names,
		 gretlopt opt, PRN *prn)
{
    gretl_matrix *coef_se = NULL;
    gretl_matrix *addstats = NULL;
    gretl_array *parnames = NULL;
    const char *parstr = NULL;
    const char **rnames = NULL;
    int free_coef_se = 0;
    int free_parnames = 0;
    int nnames = 0;
    int ncoef = 0;
    int err = 0;

    if (mname == NULL) {
	return E_ARGS;
    }

    /* k x 2 matrix: coeffs and standard errors */
    coef_se = get_matrix_by_name(mname);
    if (coef_se == NULL) {
	gretl_errmsg_set(_("modprint: expected the name of a matrix"));
	return E_INVARG;
    } else if (gretl_matrix_cols(coef_se) != 2) {
	gretl_errmsg_set(_("modprint: the first matrix argument must have 2 columns"));
	return E_INVARG;
    }

    nnames = ncoef = coef_se->rows;

    if (names != NULL) {
	/* names for coeffs: string literal, string variable,
	   or array of strings */
	if (opt & OPT_L) {
	    /* treat as string _L_iteral */
	    parstr = names;
	} else {
	    parstr = get_string_by_name(names);
	    if (parstr == NULL) {
		parnames = get_array_by_name(names);
		if (parnames == NULL ||
		    gretl_array_get_type(parnames) != GRETL_TYPE_STRINGS) {
		    err = E_TYPES;
		}
	    }
	}
    } else {
	rnames = gretl_matrix_get_rownames(coef_se);
	if (rnames == NULL) {
	    return E_ARGS;
	}
    }

    if (!err && (opt & OPT_A)) {
	/* optional third field: extra statistics */
	const char *aname = get_optval_string(MODPRINT, OPT_A);

	if (aname != NULL) {
	    addstats = get_matrix_by_name(aname);
	    if (addstats == NULL) {
		err = E_TYPES;
	    } else {
		nnames += gretl_vector_get_length(addstats);
	    }
	}
    }

    if (!err && nnames > ncoef && rnames != NULL) {
	/* for now reject this case */
	err = E_ARGS;
    }

    if (!err) {
	if (rnames != NULL) {
	    /* use matrix rownames: convert to array */
	    parnames = gretl_array_from_strings((char **) rnames, ncoef, 1, &err);
	    free_parnames = 1;
	} else if (parnames == NULL) {
	    /* we need to construct the strings array */
	    parnames = strings_array_from_string(parstr, nnames, &err);
	    free_parnames = 1;
	} else if (gretl_array_get_length(parnames) < nnames) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	PrnFormat fmt = GRETL_FORMAT_TXT;
	char fname[FILENAME_MAX];

	*fname = '\0';

	if (opt & OPT_O) {
	    /* try for --output=filename, and if found let
	       the suffix determine the output type
	    */
	    const char *s = get_optval_string(MODPRINT, OPT_O);

	    if (s != NULL && *s != '\0') {
		strcpy(fname, s);
		if (has_suffix(fname, ".tex")) {
		    fmt = GRETL_FORMAT_TEX;
		    if (opt & OPT_C) {
			fmt |= GRETL_FORMAT_DOC;
		    }
		} else if (has_suffix(fname, ".rtf")) {
		    fmt = GRETL_FORMAT_RTF;
		} else if (has_suffix(fname, ".csv")) {
		    fmt = GRETL_FORMAT_CSV;
		}
	    }
	}

	if (*fname != '\0') {
	    PRN *myprn;

	    gretl_maybe_switch_dir(fname);
	    myprn = gretl_print_new_with_filename(fname, &err);
	    if (!err) {
		gretl_print_set_format(myprn, fmt);
		err = print_model_from_matrices(coef_se, addstats,
						parnames, 0, OPT_NONE,
						myprn);
		gretl_print_destroy(myprn);
	    }
	} else {
	    gretl_print_set_format(prn, fmt);
	    err = print_model_from_matrices(coef_se, addstats,
					    parnames, 0, OPT_NONE,
					    prn);
	}
    }

    if (free_coef_se) {
	gretl_matrix_free(coef_se);
    }
    if (free_parnames) {
	gretl_array_destroy(parnames);
    }

    return err;
}

static void maybe_add_tsinfo (DATASET *mdset,
			      const gretl_matrix *m,
			      const DATASET *dset)
{
    int ts_ok = 0;

    if (dset != NULL && dataset_is_time_series(dset)) {
	int t1, t2;

	if (m->rows == dset->n) {
	    t1 = 0;
	    t2 = dset->n - 1;
	    ts_ok = 1;
	} else {
	    /* subsampled? */
	    t1 = gretl_matrix_get_t1(m);
	    t2 = gretl_matrix_get_t2(m);
	    if (t2 > t1 && t2 < dset->n) {
		ts_ok = 1;
	    }
	}
	if (ts_ok) {
	    mdset->pd = dset->pd;
	    mdset->structure = dset->structure;
	    ntolabel(mdset->stobs, t1, dset);
	    ntolabel(mdset->endobs, t2, dset);
	    mdset->sd0 = get_date_x(mdset->pd, mdset->stobs);
	}
    }
}

int matrix_command_driver (int ci,
			   const int *list,
			   const char *param,
			   const DATASET *dset,
			   gretlopt opt,
			   PRN *prn)
{
    gretl_matrix *m = NULL;
    DATASET *mdset = NULL;
    int *collist = NULL;
    const char *mname;
    int err = 0;

    mname = get_optval_string(ci, OPT_X);
    if (mname != NULL) {
	m = get_matrix_by_name(mname);
    }
    if (gretl_is_null_matrix(m)) {
	return E_DATA;
    }

    if (ci == SCATTERS) {
	/* note: this is a special case, for now */
	return matrix_multi_plots(m, list, dset, opt);
    } else if (list != NULL && list[0] == 0) {
	/* use all columns of the matrix */
	mdset = gretl_dataset_from_matrix(m, NULL, OPT_B, &err);
    } else if (list != NULL && list[0] == 1 && ci == SUMMARY) {
	/* summary stats for a single specified column */
	mdset = gretl_dataset_from_matrix(m, list, OPT_B | OPT_N, &err);
    } else {
	/* note that a NULL list is OK here */
	mdset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
    }

    if (!err) {
	dataset_set_matrix_name(mdset, mname);
	collist = gretl_consecutive_list_new(1, mdset->v - 1);
	if (collist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err && ci == GNUPLOT) {
	if (opt & OPT_T) {
	    maybe_add_tsinfo(mdset, m, dset);
	}
	err = gnuplot(collist, param, mdset, opt);
    } else if (!err) {
	opt &= ~OPT_X;
	if (ci == BXPLOT) {
	    err = boxplots(collist, param, mdset, opt);
	} else if (ci == SUMMARY) {
	    err = list_summary(collist, 0, mdset, opt, prn);
	} else if (ci == CORR) {
	    err = gretl_corrmx(collist, mdset, opt, prn);
	} else {
	    err = E_DATA;
	}
    }

    destroy_dataset(mdset);
    free(collist);

    return err;
}

int matrix_freq_driver (const int *list,
			gretlopt opt,
			PRN *prn)
{
    gretl_matrix *m = NULL;
    DATASET *mdset = NULL;
    const char *mname;
    int err = 0;

    if (list != NULL && list[0] != 1) {
	return E_DATA;
    }

    mname = get_optval_string(FREQ, OPT_X);
    if (mname != NULL) {
	m = get_matrix_by_name(mname);
    }

    if (gretl_is_null_matrix(m)) {
	err = E_DATA;
    } else {
	if (list == NULL) {
	    /* this is OK if m is a column vector */
	    if (m->cols == 1) {
		int mlist[2] = {1, 1};

		mdset = gretl_dataset_from_matrix(m, mlist, OPT_B, &err);
	    } else {
		err = E_ARGS;
	    }
	} else {
	    mdset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
	}
    }

    if (!err) {
	err = freqdist(1, mdset, opt, prn);
    }

    destroy_dataset(mdset);

    return err;
}

int list_summary_driver (const int *list, const DATASET *dset,
			 gretlopt opt, PRN *prn)
{
    int wtvar = 0;
    int err = 0;

    if (opt & OPT_W) {
	const char *wname = get_optval_string(SUMMARY, OPT_W);

	if (wname == NULL) {
	    err = E_DATA;
	} else {
	    wtvar = current_series_index(dset, wname);
	    if (wtvar < 0) {
		err = E_UNKVAR;
	    }
	}
    }

    if (!err) {
	err = list_summary(list, wtvar, dset, opt, prn);
    }

    return err;
}
