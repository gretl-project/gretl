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

/* mechanisms for defining and handling systems of equations */

#include "libgretl.h"
#include "system.h"
#include "objstack.h"
#include "gretl_xml.h"

enum {
    OP_PLUS,
    OP_MINUS
} identity_ops;

enum {
    ENDOG_LIST,
    INSTR_LIST
} aux_list_types;

enum {
    SYS_NO_TEST,
    SYS_LR_TEST,
    SYS_F_TEST,
    SYS_TEST_NOTIMP
} system_test_types;

struct id_atom_ {
    int op;         /* operator (plus or miinus) */
    int varnum;     /* ID number of variable to right of operator */
};

struct identity_ {
    int n_atoms;    /* number of "atomic" elements in identity */
    int depvar;     /* LHS variable in indentity */
    id_atom *atoms; /* pointer to RHS "atoms" */
};

const char *gretl_system_method_strings[] = {
    "sur",
    "3sls",
    "fiml",
    "liml",
    "ols",
    "tsls",
    "wls",
    NULL
};

const char *gretl_system_short_strings[] = {
    N_("SUR"),
    N_("3SLS"),
    N_("FIML"),
    N_("LIML"),
    N_("OLS"),
    N_("TSLS"),
    N_("WLS"),
    NULL
};

const char *gretl_system_long_strings[] = {
    N_("Seemingly Unrelated Regressions"),
    N_("Three-Stage Least Squares"),
    N_("Full Information Maximum Likelihood"),
    N_("Limited Information Maximum Likelihood"),
    N_("Ordinary Least Squares"),
    N_("Two-Stage Least Squares"),
    N_("Weighted Least Squares"),
    NULL
};

const char *nosystem = N_("No system of equations has been defined");
const char *badsystem = N_("Unrecognized equation system type");
const char *toofew = N_("An equation system must have at least two equations");

static void destroy_ident (identity *ident);
static int make_instrument_list (gretl_equation_system *sys);

static void 
print_system_equation (const int *list, const DATAINFO *pdinfo, 
		       PRN *prn)
{
    int i, v;

    pputs(prn, "equation");

    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v == LISTSEP) {
	    pputs(prn, " ;");
	} else if (v > 0 && v < pdinfo->v) {
	    pprintf(prn, " %s", pdinfo->varname[v]);
	} else {
	    pprintf(prn, " %d", v);
	}
    }

    pputc(prn, '\n');
}

static void 
print_system_identity (const identity *ident, const DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    int i;

    if (opt & OPT_H) {
	pprintf(prn, "Identity: %s = %s ", 
		pdinfo->varname[ident->depvar],
		pdinfo->varname[ident->atoms[0].varnum]);
    } else {
	pprintf(prn, "identity %s = %s ", 
		pdinfo->varname[ident->depvar],
		pdinfo->varname[ident->atoms[0].varnum]);
    }

    for (i=1; i<ident->n_atoms; i++) {
	pprintf(prn, "%c %s ", (ident->atoms[i].op == OP_PLUS)? '+' : '-',
		pdinfo->varname[ident->atoms[i].varnum]);
    }

    pputc(prn, '\n');
}

/**
 * print_equation_system_info:
 * @sys: gretl equation system.
 * @pdinfo: dataset information.
 * @opt: use %OPT_H for printing in a form designed to appear
 * in the header when results of estimation are printed.
 * @prn: printing struct.
 * 
 * Prints details of an equation system to @prn.
 */

void 
print_equation_system_info (const gretl_equation_system *sys, 
			    const DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn)
{
    int i;

    
    if ((opt & OPT_H) && sys->name != NULL) {
	pprintf(prn, "Equation system %s\n", sys->name);
    }

    if (!(opt & OPT_H)) {
	for (i=0; i<sys->n_equations; i++) {
	    print_system_equation(sys->lists[i], pdinfo, prn);
	}    
    }

    for (i=0; i<sys->n_identities; i++) {
	print_system_identity(sys->idents[i], pdinfo, opt, prn);
    }

    if (sys->endog_vars != NULL) {
	pputs(prn, (opt & OPT_H)? "Endogenous variables:" : "endog");
	for (i=1; i<=sys->endog_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->endog_vars[i]]);
	}
	pputc(prn, '\n');
    }

    if (sys->instr_vars != NULL) {
	pputs(prn, (opt & OPT_H)? "Exogenous variables:" : "instr");
	for (i=1; i<=sys->instr_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->instr_vars[i]]);
	}
	pputc(prn, '\n');
    }
}

int gretl_system_method_from_string (const char *s)
{
    int i = 0;

    while (gretl_system_method_strings[i] != NULL) {
	if (!strcmp(s, gretl_system_method_strings[i])) {
	    return i;
	}
	i++;
    }

    return i;
}

const char *system_method_full_string (int method)
{
    if (method >= SYS_SUR && method < SYS_MAX) {
	return gretl_system_long_strings[method];
    } else {
	return NULL;
    }
}

const char *system_method_short_string (int method)
{
    if (method >= SYS_SUR && method < SYS_MAX) {
	return gretl_system_method_strings[method];
    } else {
	return NULL;
    }
}

static gretl_equation_system *
gretl_equation_system_new (int method, const char *name)
{
    gretl_equation_system *sys;

    if (method < 0 && name == NULL) {
	return NULL;
    }

    sys = malloc(sizeof *sys);
    if (sys == NULL) return NULL;

    if (name != NULL) {
	sys->name = gretl_strdup(name);
	if (sys->name == NULL) {
	    free(sys);
	    return NULL;
	}
    } else {
	sys->name = NULL;
    }

    sys->refcount = 0;
    sys->method = method;

    sys->t1 = sys->t2 = 0;

    sys->n_equations = 0;
    sys->n_identities = 0;

    sys->R = NULL;
    sys->q = NULL;

    sys->n_obs = 0;
    sys->iters = 0;
    sys->flags = 0;

    sys->ll = sys->llu = NADBL;
    sys->X2 = NADBL;
    sys->ess = NADBL;
    sys->diag = 0.0;
    sys->bdiff = 0.0;

    sys->b = NULL;
    sys->vcv = NULL;
    sys->sigma = NULL;
    sys->uhat = NULL;

    sys->lists = NULL;
    sys->endog_vars = NULL;
    sys->instr_vars = NULL;
    sys->idents = NULL;

    sys->models = NULL;

    return sys;
}

static void system_clear_restrictions (gretl_equation_system *sys)
{
    if (sys->R != NULL) {
	free(sys->R);
	sys->R = NULL;
    }

    if (sys->q != NULL) {
	free(sys->q);
	sys->q = NULL;
    }

    sys->flags &= ~GRETL_SYS_RESTRICT;
}

static void system_clear_results (gretl_equation_system *sys)
{
    sys->iters = 0;

    sys->t1 = sys->t2 = 0;

    sys->ll = NADBL;
    sys->llu = NADBL;
    sys->X2 = NADBL;
    sys->ess = NADBL;

    sys->diag = 0.0;
    sys->bdiff = 0.0;

    if (sys->b != NULL) {
	gretl_matrix_free(sys->b);
	sys->b = NULL;
    }

    if (sys->vcv != NULL) {
	gretl_matrix_free(sys->vcv);
	sys->vcv = NULL;
    }

    if (sys->sigma != NULL) {
	gretl_matrix_free(sys->sigma);
	sys->sigma = NULL;
    }

    if (sys->uhat != NULL) {
	gretl_matrix_free(sys->uhat);
	sys->uhat = NULL;
    }
}

void gretl_equation_system_destroy (gretl_equation_system *sys)
{
    int i;

    if (sys == NULL || sys->lists == NULL) {
	return;
    }

    sys->refcount -= 1;
    if (sys->refcount > 0) {
	return;
    }

    for (i=0; i<sys->n_equations; i++) {
	free(sys->lists[i]);
    }
    free(sys->lists);
    sys->lists = NULL;

    for (i=0; i<sys->n_identities; i++) {
	destroy_ident(sys->idents[i]);
    }
    free(sys->idents);

    free(sys->endog_vars);
    free(sys->instr_vars);

    free(sys->name);

    if (sys->R != NULL) {
	gretl_matrix_free(sys->R);
    }

    if (sys->q != NULL) {
	gretl_matrix_free(sys->q);
    }

    system_clear_results(sys);

    free(sys);
}

static void sur_rearrange_lists (gretl_equation_system *sys,
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    if (sys->method == SYS_SUR) {
	int i;

	for (i=0; i<sys->n_equations; i++) {
	    reglist_check_for_const(sys->lists[i], Z, pdinfo);
	}
    }
}

/**
 * gretl_equation_system_append:
 * @sys: initialized equation system.
 * @list: list containing dependent variable and regressors.
 * 
 * Adds an equation (as represented by @list) to @sys.
 * 
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int gretl_equation_system_append (gretl_equation_system *sys, 
				  const int *list)
{
    int i, neq;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    neq = sys->n_equations;

    sys->lists = realloc(sys->lists, (neq + 1) * sizeof *sys->lists);
    if (sys->lists == NULL) return E_ALLOC;

    sys->lists[neq] = gretl_list_new(list[0]);

    if (sys->lists[neq] == NULL) {
	gretl_equation_system_destroy(sys);
	return E_ALLOC;
    }

    for (i=0; i<=list[0]; i++) {
	sys->lists[neq][i] = list[i];
    }

    sys->n_equations += 1;

    return 0;
}

/* retrieve the name -- possible quoted with embedded spaces -- for
   an equation system */

char *get_system_name_from_line (const char *s)
{
    char *name = NULL;
    const char *p = strstr(s, " name");
    int pchars = 0;

    if (p != NULL) {
	p += 5;
    } else {
	p = strstr(s, "estimate ");
	if (p == NULL) {
	    p = strstr(s, "restrict ");
	}
	if (p != NULL) {
	    p += 9;
	}	
    }

    if (p == NULL) {
	return NULL;
    }

    s = p;
    while (isspace((unsigned char) *s) || *s == '=') {
	s++;
    }

    if (*s == '"') {
	if (*(s + 1) != '\0') s++;
	p = s;
	while (*p && *p != '"') {
	    if (!isspace((unsigned char) *p)) pchars++;
	    p++;
	}
	if (*p != '"') {
	    /* no closing quote */
	    pchars = 0;
	}
    } else {
	p = s;
	while (*p && !isspace((unsigned char) *p)) {
	    pchars++;
	    p++;
	}
    }

    if (pchars > 0) {
	name = gretl_strndup(s, p - s);
    }

    return name;
}

/* parse a system estimation method out of a command line */

static int get_estimation_method_from_line (const char *s)
{
    const char *p = strstr(s, " method");
    char mstr[9];
    int method = -1;

    if (p != NULL) {
	p += 7;
    } else {
	p = strstr(s, " type");
	if (p != NULL) {
	    p += 5;
	}
    }

    if (p == NULL) {
	return -1;
    }

    while (isspace((unsigned char) *p) || *p == '=') {
	p++;
    }

    if (sscanf(p, "%8s", mstr) == 1) {
	lower(mstr);
	method = gretl_system_method_from_string(mstr);
    }

    return method;
}

static void 
system_set_flags (gretl_equation_system *sys, const char *s, gretlopt opt)
{
    if (opt & OPT_T) {
	sys->flags |= GRETL_SYS_ITERATE;
    }

    s = strstr(s, " save");

    if (s != NULL) {
	s += 5;
	if (*s != ' ' && *s != '=') {
	    return;
	}
	if (strstr(s, "resids") || strstr(s, "uhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_UHAT;
	}
	if (strstr(s, "fitted") || strstr(s, "yhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_YHAT;
	}
    }
}

/**
 * system_start:
 * @line: command line.
 * @opt: may include %OPT_T for iterative estimation, if the
 * estimation method supports this.
 * 
 * Start compiling an equation system: @line must specify a "method" 
 * (estimation method) and/or a name for the system.  If a method is
 * given, the system will be estimated as soon as its definition is 
 * complete.  If a name is given, the system definition is saved on a 
 * stack, and it can subsequently be estimated via various methods
 * If both a name and an estimation method are given, the system 
 * is both estimated and saved.
 * 
 * Returns: pointer to a new equation system, or %NULL on error.
 */

gretl_equation_system *system_start (const char *line, gretlopt opt)
{
    gretl_equation_system *sys = NULL;
    char *sysname = NULL;
    int method;

    method = get_estimation_method_from_line(line);

    if (method == SYS_MAX) {
	/* invalid method was given */
	strcpy(gretl_errmsg, _(badsystem));
	return NULL;
    }

    sysname = get_system_name_from_line(line);

    if (method < 0 && sysname == NULL) {
	/* neither a method nor a name was specified */
	strcpy(gretl_errmsg, _(badsystem));
	return NULL;
    }

    sys = gretl_equation_system_new(method, sysname);
    if (sys == NULL) {
	return NULL;
    }

    system_set_flags(sys, line, opt);

    if (sysname != NULL) {
	free(sysname);
    }

    return sys;
}

static int system_get_dfu (const gretl_equation_system *sys)
{
    int dfu = sys->n_obs * sys->n_equations;
    int i;

    for (i=0; i<sys->n_equations; i++) {
	dfu -= sys->lists[i][0] - 1;
    }

    return dfu;
}

/* Asymptotic F-test, as in Greene:

     (1/J) * (Rb-q)' * [R*Var(b)*R']^{-1} * (Rb-q) 
*/

static void 
system_print_F_test (const gretl_equation_system *sys,
		     const gretl_matrix *b, const gretl_matrix *vcv,
		     PRN *prn)
{
    const gretl_matrix *R = sys->R;
    const gretl_matrix *q = sys->q;
    int Rrows = gretl_matrix_rows(R);
    int dfu = system_get_dfu(sys);
    int dfn = gretl_matrix_rows(R);

    gretl_matrix *Rbq = NULL;
    gretl_matrix *RvR = NULL;

    double F = NADBL;
    int err = 0;

    if (R == NULL || q == NULL || b == NULL || vcv == NULL) {
	pputs(prn, "Missing data in F test\n");
	goto bailout;
    }

    Rbq = gretl_matrix_alloc(Rrows, 1);
    RvR = gretl_matrix_alloc(Rrows, Rrows);

    if (Rbq == NULL || RvR == NULL) {
	pputs(prn, "Out of memory in F test\n");
	goto bailout;
    }

    gretl_matrix_multiply(R, b, Rbq);
    gretl_matrix_subtract_from(Rbq, q);

    gretl_matrix_qform(R, GRETL_MOD_NONE, vcv,
		       RvR, GRETL_MOD_NONE);

    err = gretl_invert_symmetric_matrix(RvR);
    if (err) {
	pputs(prn, "Matrix inversion failed in F test\n");
	goto bailout;
    }

    F = gretl_scalar_qform(Rbq, RvR, &err);
    if (err) {
	pputs(prn, "Matrix multiplication failed in F test\n");
	goto bailout;
    }

    F /= dfn;

    pprintf(prn, "%s:\n", _("F test for the specified restrictions"));
    pprintf(prn, "  F(%d,%d) = %g %s %g\n", dfn, dfu, F,
	    _("with p-value"), snedecor_cdf_comp(F, dfn, dfu));
    pputc(prn, '\n');    

 bailout:
    
    gretl_matrix_free(Rbq);
    gretl_matrix_free(RvR);
}

static void 
system_print_LR_test (const gretl_equation_system *sys,
		      double llu, PRN *prn)
{
    double llr = sys->ll;
    int df = gretl_matrix_rows(sys->R);
    double X2;

    if (na(llr) || na(llu) || llr == 0.0 || llu == 0.0 || df <= 0) {
	fputs("bad or missing data in system LR test\n", stderr);
	return;
    }

    X2 = 2.0 * (llu - llr);

    pprintf(prn, "%s:\n", _("LR test for the specified restrictions"));
    pprintf(prn, "  %s = %g\n", _("Restricted log-likelihood"), llr);
    pprintf(prn, "  %s = %g\n", _("Unrestricted log-likelihood"), llu);
    pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
	    df, X2, _("with p-value"), chisq_cdf_comp(X2, df));
    pputc(prn, '\n');    
}

static int sys_test_type (gretl_equation_system *sys)
{
    int ret = SYS_NO_TEST;

    if (system_n_restrictions(sys) > 0) {
	if (sys->method == SYS_SUR || sys->method == SYS_WLS) {
	    if (sys->flags & GRETL_SYS_ITERATE) {
		ret = SYS_LR_TEST;
	    } else {
		ret = SYS_F_TEST;
	    }
	} else if (sys->method == SYS_OLS ||
		   sys->method == TSLS ||
		   sys->method == SYS_3SLS) {
	    ret = SYS_F_TEST;
	} else if (sys->method == SYS_LIML) {
	    /* experimental */
	    ret = SYS_F_TEST;
	} else if (sys->method == SYS_FIML) {
	    ret = SYS_LR_TEST;
	} 
    }

    return ret;
}

static int estimate_with_test (gretl_equation_system *sys, 
			       double ***pZ, DATAINFO *pdinfo, 
			       gretlopt opt, int stest, 
			       int (*system_est)(), PRN *prn)
{
    gretl_matrix *vcv = NULL;
    gretl_matrix *b = NULL;
    double llu = 0.0;
    int err = 0;

    /* estimate the unrestricted system first */

    sys->flags &= ~GRETL_SYS_RESTRICT;

    if (stest == SYS_F_TEST) {
	/* save unrestricted coeffs and vcv */
	sys->flags |= GRETL_SYS_SAVE_VCV;
    }

    err = (* system_est) (sys, pZ, pdinfo, opt | OPT_Q, prn);

    if (stest == SYS_F_TEST) {
	sys->flags ^= GRETL_SYS_SAVE_VCV;
    }

    sys->flags ^= GRETL_SYS_RESTRICT;

    if (err) {
	goto bailout;
    }

    /* grab the data from unrestricted estimation */

    if (stest == SYS_LR_TEST) {
	llu = sys->ll;
    } else if (stest == SYS_F_TEST) {
	b = sys->b;
	sys->b = NULL;
	vcv = sys->vcv;
	sys->vcv = NULL;
    }

    /* now estimate the restricted system */

    system_clear_results(sys);
    err = (* system_est) (sys, pZ, pdinfo, opt, prn);

    if (!err) {
	if (stest == SYS_LR_TEST) {
	    system_print_LR_test(sys, llu, prn);
	} else if (stest == SYS_F_TEST) {
	    system_print_F_test(sys, b, vcv, prn);
	}
    }

 bailout:

    if (b != NULL) gretl_matrix_free(b);
    if (vcv != NULL) gretl_matrix_free(vcv);

    return err;
}

static void 
adjust_sys_flags_for_method (gretl_equation_system *sys, int method)
{
    char oldflags = sys->flags;

    sys->flags = 0;

    if (oldflags & GRETL_SYS_ITERATE) {
	/* the iterate option is only available for WLS, SUR or 3SLS */
	if (sys->method == SYS_WLS || sys->method == SYS_SUR ||
	    sys->method == SYS_3SLS) {
	    sys->flags |= GRETL_SYS_ITERATE;
	}
    }

    /* by default, apply a df correction for single-equation methods */
    if (sys->method == SYS_OLS || sys->method == SYS_WLS ||
	sys->method == SYS_TSLS || sys->method == SYS_LIML) {
	if (oldflags & GRETL_SYSTEM_DFCORR) {
	    sys->flags |= GRETL_SYSTEM_DFCORR;
	}
    } 

    if (oldflags & GRETL_SYS_VCV_GEOMEAN) {
	sys->flags |= GRETL_SYS_VCV_GEOMEAN;
    }    
}

static void 
set_sys_flags_from_opt (gretl_equation_system *sys, gretlopt opt)
{
    char oldflags = sys->flags;

    sys->flags = 0;

    if (opt & OPT_T) {
	/* the iterate option is only available for WLS, SUR or 3SLS */
	if (sys->method == SYS_WLS || sys->method == SYS_SUR ||
	    sys->method == SYS_3SLS) {
	    sys->flags |= GRETL_SYS_ITERATE;
	}
    }

    /* by default, apply a df correction for single-equation methods */
    if (sys->method == SYS_OLS || sys->method == SYS_WLS ||
	sys->method == SYS_TSLS || sys->method == SYS_LIML) {
	if (!(opt & OPT_N)) {
	    sys->flags |= GRETL_SYSTEM_DFCORR;
	}
    } 

    if (opt & OPT_M) {
	sys->flags |= GRETL_SYS_VCV_GEOMEAN;
    } 

    if (oldflags & GRETL_SYSTEM_SAVE_UHAT) {
	sys->flags |= GRETL_SYSTEM_SAVE_UHAT;
    }

    if (oldflags & GRETL_SYSTEM_SAVE_YHAT) {
	sys->flags |= GRETL_SYSTEM_SAVE_YHAT;
    }
}

/**
 * gretl_equation_system_estimate:
 * @sys: pre-defined equation system.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: may include %OPT_Q for relatively quiet operation.
 * @prn: printing struct.
 * 
 * Estimate a pre-defined equation system and print the results
 * to @prn.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int 
gretl_equation_system_estimate (gretl_equation_system *sys, 
				double ***pZ, DATAINFO *pdinfo, 
				gretlopt opt, PRN *prn)
{
    int err = 0;
    void *handle = NULL;
    int (*system_est) (gretl_equation_system *, 
		       double ***, DATAINFO *, gretlopt, PRN *);
    int stest = 0;

    gretl_error_clear();

    if (opt == OPT_UNSET) {
	opt = OPT_NONE;
	adjust_sys_flags_for_method(sys, sys->method);
    } else {
	set_sys_flags_from_opt(sys, opt);
    }

    /* if restrictions are in place, reset the restricted flag */
    if (sys->R != NULL && sys->q != NULL) {
	sys->flags |= GRETL_SYS_RESTRICT;
    } 

    /* in case we're re-estimating */
    system_clear_results(sys);

    err = make_instrument_list(sys);
    if (err) goto system_bailout;
    
    if (sys->method == SYS_SUR) {
	sur_rearrange_lists(sys, (const double **) *pZ, pdinfo);
    }

    system_est = get_plugin_function("system_estimate", &handle);

    if (system_est == NULL) {
	err = 1;
	goto system_bailout;
    }

    stest = sys_test_type(sys);

    if (stest == SYS_TEST_NOTIMP) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
	err = 1;
    } else if (stest != SYS_NO_TEST) {
	err = estimate_with_test(sys, pZ, pdinfo, opt, stest, 
				 system_est, prn);
    } else {
	err = (* system_est) (sys, pZ, pdinfo, opt, prn);
    }
    
 system_bailout:

    if (handle != NULL) {
	close_plugin(handle);
    }

    if (!err) {
	set_as_last_model(sys, GRETL_OBJ_SYS);
    } 

    return err;
}

/**
 * gretl_equation_system_finalize:
 * @sys: pre-defined equation system.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @prn: printing struct.
 * 
 * Finalize an equation system, e.g. in response to "end system".
 * If the system has a specified name, we save it on a stack of 
 * defined systems.  If it has a specified estimation method, 
 * we go ahead and estimate it.  If it has both, we do both.
 * 
 * Returns: 0 on success, non-zero on error.  If the system is
 * mal-formed, it is destroyed.
 */

int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    gretlopt opt = OPT_NONE;
    int err = 0;

    gretl_error_clear();

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    if (sys->n_equations < 2) {
	strcpy(gretl_errmsg, _(toofew));
	gretl_equation_system_destroy(sys);
	return 1;
    }

    if (sys->method >= SYS_MAX) {
	strcpy(gretl_errmsg, _(badsystem));
	gretl_equation_system_destroy(sys);
	return 1;
    }    

    if (sys->name != NULL) {
	/* save the system for subsequent estimation */
	err = gretl_stack_object_as(sys, GRETL_OBJ_SYS, sys->name);
    }

    if (!err && sys->method >= 0) {
	err = gretl_equation_system_estimate(sys, pZ, pdinfo, opt, prn);
    }

    return err;
}


/* Implement the "estimate" command, which must give the name of a pre-defined
   equation system and an estimation method, as in:

             estimate "Klein Model 1" method=FIML 
*/

int estimate_named_system (const char *line, double ***pZ, DATAINFO *pdinfo, 
			   gretlopt opt, PRN *prn)
{
    gretl_equation_system *sys;
    char *sysname;
    int method;

    if (strlen(line) < 12) {
	strcpy(gretl_errmsg, "estimate: no system name was provided");
	return 1;
    }

    sysname = get_system_name_from_line(line);
    if (sysname == NULL) {
	strcpy(gretl_errmsg, "estimate: no system name was provided");
	return 1;
    }

    sys = get_equation_system_by_name(sysname);
    if (sys == NULL) {
	sprintf(gretl_errmsg, _("'%s': unrecognized name"), sysname);
	free(sysname);
	return 1;
    }

    free(sysname);

    method = get_estimation_method_from_line(line);
    if (method < 0 || method >= SYS_MAX) {
	method = sys->method;
    }

    if (method < 0 || method >= SYS_MAX) {
	strcpy(gretl_errmsg, "estimate: no valid method was specified");
	return 1;
    }

    sys->method = method;

    return gretl_equation_system_estimate(sys, pZ, pdinfo, opt, prn);
}

static int get_real_list_length (const int *list)
{
    int i, len = list[0];

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    len = i - 1;
	    break;
	}
    }

    return len;
}

int system_max_indep_vars (const gretl_equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->n_equations; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	if (nvi > nv) nv = nvi;
    }

    return nv;
}

int system_n_indep_vars (const gretl_equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->n_equations; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	nv += nvi;
    }

    return nv;
}

const char *gretl_system_short_string (const MODEL *pmod)
{
    int i = gretl_model_get_int(pmod, "method");

    return gretl_system_short_strings[i];
}

void gretl_system_set_name (gretl_equation_system *sys, const char *name)
{
    if (sys->name != NULL && !strcmp(sys->name, name)) {
	return;
    }

    if (sys->name != NULL) {
	free(sys->name);
    }

    sys->name = gretl_strdup(name);
}

int system_adjust_t1t2 (gretl_equation_system *sys,
			int *t1, int *t2, const double **Z)
{
    int i, err = 0;

    for (i=0; i<sys->n_equations && !err; i++) {
	err = check_for_missing_obs(sys->lists[i], t1, t2, Z, NULL);
    }

    if (!err) {
	sys->t1 = *t1;
	sys->t2 = *t2;
	sys->n_obs = *t2 - *t1 + 1;
    }

    return err;
}

int *compose_tsls_list (gretl_equation_system *sys, int i)
{
    int *list;
    int j, k1, k2;

    if (i >= sys->n_equations) {
	return NULL;
    }

    if (sys->instr_vars == NULL && make_instrument_list(sys)) {
	return NULL;
    }

    k1 = sys->lists[i][0];
    k2 = sys->instr_vars[0];

    list = malloc((k1 + k2 + 2) * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    list[0] = k1 + k2 + 1;

    for (j=1; j<=list[0]; j++) {
	if (j <= k1) {
	    list[j] = sys->lists[i][j];
	} else if (j == k1 + 1) {
	    list[j] = LISTSEP;
	} else {
	    list[j] = sys->instr_vars[j - (k1 + 1)];
	}
    }

    return list;
}

int system_normality_test (const gretl_equation_system *sys, PRN *prn)
{
    int err = 0;

    if (sys->uhat == NULL || sys->sigma == NULL) {
	err = 1;
    } else {
	err = gretl_system_normality_test(sys->uhat, sys->sigma, prn);
    }

    return err;
}

void make_system_data_info (gretl_equation_system *sys, int eqn, 
			    DATAINFO *pdinfo, int v, int code)
{
    if (code == GRETL_SYSTEM_SAVE_UHAT) {
	sprintf(pdinfo->varname[v], "uhat_s%02d", eqn);
	if (sys->method == SYS_SUR) {
	    sprintf(VARLABEL(pdinfo, v), _("SUR residual, equation %d"), eqn);
	} else if (sys->method == SYS_3SLS) {
	    sprintf(VARLABEL(pdinfo, v), _("3SLS residual, equation %d"), eqn);
	} else {
	    sprintf(VARLABEL(pdinfo, v), _("system residual, equation %d"), eqn);
	}
    } else if (code == GRETL_SYSTEM_SAVE_YHAT) {
	sprintf(pdinfo->varname[v], "yhat_s%02d", eqn);
	if (sys->method == SYS_SUR) {
	    sprintf(VARLABEL(pdinfo, v), _("SUR fitted value, equation %d"), eqn);
	} else if (sys->method == SYS_3SLS) {
	    sprintf(VARLABEL(pdinfo, v), _("3SLS fitted value, equation %d"), eqn);
	} else {
	    sprintf(VARLABEL(pdinfo, v), _("system fitted value, equation %d"), eqn);
	}
    }
}

int gretl_system_add_resids_to_dataset (gretl_equation_system *sys, int eqnum,
					double ***pZ, DATAINFO *pdinfo)
{
    int v, t;

    if (sys->uhat == NULL) {
	return E_DATA;
    }

    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    v = pdinfo->v - 1;

    for (t=0; t<pdinfo->n; t++) {
	if (t < sys->t1 || t > sys->t2) {
	    (*pZ)[v][t] = NADBL;
	} else {
	    (*pZ)[v][t] = gretl_matrix_get(sys->uhat, t - sys->t1, eqnum);
	}
    }

    make_system_data_info(sys, eqnum + 1, pdinfo, v, GRETL_SYSTEM_SAVE_UHAT);

    return 0;
}

/* simple accessor functions */

const char *system_get_full_string (const gretl_equation_system *sys)
{
    if (sys->flags & GRETL_SYS_ITERATE) {
	static char sysstr[64];

	sprintf(sysstr, _("iterated %s"), gretl_system_long_strings[sys->method]);
	return sysstr;
    } else {
	return gretl_system_long_strings[sys->method];
    }
}

int system_save_uhat (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYSTEM_SAVE_UHAT;
}

int system_save_yhat (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYSTEM_SAVE_YHAT;
}

int system_save_vcv (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYS_SAVE_VCV;
}

int system_want_df_corr (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYSTEM_DFCORR;
}

int system_n_restrictions (const gretl_equation_system *sys)
{
    int nr = 0;

    if (sys->R != NULL && (sys->flags & GRETL_SYS_RESTRICT)) {
	nr = gretl_matrix_rows(sys->R);
    }

    return nr;
}

int *system_get_list (const gretl_equation_system *sys, int i)
{
    if (i >= sys->n_equations) return NULL;

    return sys->lists[i];
}

int system_get_depvar (const gretl_equation_system *sys, int i)
{
    if (i >= sys->n_equations) return 0;

    return sys->lists[i][1];
}

int system_get_method (const gretl_equation_system *sys)
{
    return sys->method;
}

int *system_get_endog_vars (const gretl_equation_system *sys)
{
    return sys->endog_vars;
}

int *system_get_instr_vars (const gretl_equation_system *sys)
{
    return sys->instr_vars;
}

void system_attach_coeffs (gretl_equation_system *sys, gretl_matrix *b)
{
    if (sys->b != NULL) {
	gretl_matrix_free(sys->b);
    }

    sys->b = b;
}

void system_attach_vcv (gretl_equation_system *sys, gretl_matrix *vcv)
{
    if (sys->vcv != NULL) {
	gretl_matrix_free(sys->vcv);
    }

    sys->vcv = vcv;
}

void system_attach_sigma (gretl_equation_system *sys, gretl_matrix *sigma)
{
    if (sys->sigma != NULL) {
	gretl_matrix_free(sys->sigma);
    }

    sys->sigma = sigma;
}

void system_attach_uhat (gretl_equation_system *sys, gretl_matrix *uhat)
{
    if (sys->uhat != NULL) {
	gretl_matrix_free(sys->uhat);
    }

    sys->uhat = uhat;
}

MODEL *system_get_model (const gretl_equation_system *sys, int i)
{
    if (sys->models == NULL || i >= sys->n_equations) {
	return NULL; 
    }   

    return sys->models[i];
}

/* for case of applying df correction to cross-equation 
   residual covariance matrix */

int system_vcv_geomean (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYS_VCV_GEOMEAN;
}

static int sys_eqn_indep_coeffs (const gretl_equation_system *sys, int eq)
{
    int nc;

    if (eq >= sys->n_equations) {
	nc = 0;
    } else if (sys->R == NULL) {
	nc = sys->lists[eq][0] - 1;
    } else {
	int nr = gretl_matrix_rows(sys->R);
	int cols = gretl_matrix_cols(sys->R);
	int cmax, cmin = 0;
	int single, i, j;

	nc = sys->lists[eq][0] - 1;

	for (i=0; i<eq; i++) {
	    cmin += sys->lists[i][0] - 1;
	}

	cmax = cmin + nc;

	for (i=0; i<nr; i++) {
	    single = 1;
	    for (j=0; j<cols; j++) {
		if (j >= cmin && j < cmax) {
		    if (gretl_matrix_get(sys->R, i, j) != 0.0) {
			single = 0;
			break;
		    }
		}
	    }
	    if (single) {
		nc--;
	    }
	}
    }

    return nc;
}

double 
system_vcv_denom (const gretl_equation_system *sys, int i, int j)
{
    double den = sys->n_obs;

    if ((sys->flags & GRETL_SYS_VCV_GEOMEAN) &&
	i < sys->n_equations && j < sys->n_equations) {
	int ki = sys_eqn_indep_coeffs(sys, i);

	if (j == i) {
	    den = sys->n_obs - ki;
	} else {
	    int kj = sys_eqn_indep_coeffs(sys, j);

	    den = (sys->n_obs - ki) * (sys->n_obs - kj);
	    den = sqrt(den);
	}
    }

    return den;
}

/* for system over-identification test */

int system_get_overid_df (const gretl_equation_system *sys)
{
    int gl, i, k = 0;

    gl = sys->n_equations * sys->instr_vars[0];

    for (i=0; i<sys->n_equations; i++) {
	k += sys->lists[i][0] - 1;
    }

    return gl - k;
}

/* dealing with identities (FIML, LIML) */

int rhs_var_in_identity (const gretl_equation_system *sys, int lhsvar,
			 int rhsvar)
{
    const identity *ident;
    int i, j;

    for (i=0; i<sys->n_identities; i++) {
	ident = sys->idents[i];
	if (ident->depvar == lhsvar) {
	    for (j=0; j<ident->n_atoms; j++) {
		if (ident->atoms[j].varnum == rhsvar) {
		    return (ident->atoms[j].op == OP_PLUS)? 1 : -1;
		}
	    }
	}
    }

    return 0;
}

static void destroy_ident (identity *ident)
{
    free(ident->atoms);
    free(ident);
}

static identity *ident_new (int nv)
{
    identity *ident;

    ident = malloc(sizeof *ident);
    if (ident == NULL) return NULL;

    ident->n_atoms = nv;
    ident->atoms = malloc(nv * sizeof *ident->atoms);
    if (ident->atoms == NULL) {
	free(ident);
	ident = NULL;
    }

    return ident;
}

static identity *
parse_identity (const char *str, const DATAINFO *pdinfo, int *err)
{
    identity *ident;
    const char *p;
    char f1[24], f2[16];
    char op, vname1[VNAMELEN], vname2[VNAMELEN];
    int i, nv;

    sprintf(f1, "%%%ds = %%%d[^+ -]", VNAMELEN - 1, VNAMELEN - 1);
    sprintf(f2, "%%c %%%d[^+ -]", VNAMELEN - 1);

    if (sscanf(str, f1, vname1, vname2) != 2) {
	*err = E_PARSE;
	return NULL;
    }

    p = str;
    nv = 1;
    while (*p) {
	if (*p == '+' || *p == '-') nv++;
	p++;
    }

    ident = ident_new(nv);
    if (ident == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    ident->depvar = varindex(pdinfo, vname1);
    if (ident->depvar == pdinfo->v) {
	destroy_ident(ident);
	*err = E_UNKVAR;
	return NULL;
    }

    ident->atoms[0].op = OP_PLUS;
    ident->atoms[0].varnum = varindex(pdinfo, vname2);
    if (ident->atoms[0].varnum == pdinfo->v) {
	destroy_ident(ident);
	*err = E_UNKVAR;
	return NULL;
    }

    p = str;
    for (i=1; i<nv && !*err; i++) {
	p += strcspn(p, "+-");
	sscanf(p, f2, &op, vname1);
	if (op == '+') op = OP_PLUS;
	else if (op == '-') op = OP_MINUS;
	else *err = E_PARSE;
	if (!*err) {
	    ident->atoms[i].op = op;
	    ident->atoms[i].varnum = varindex(pdinfo, vname1);
	    if (ident->atoms[i].varnum == pdinfo->v) {
		*err = E_UNKVAR;
	    }
	}
	p++;
    }

    if (*err) {
	destroy_ident(ident);
	ident = NULL;
    }
       
    return ident;
}

static int 
add_identity_to_sys (gretl_equation_system *sys, const char *line,
		     const DATAINFO *pdinfo)
{
    identity **pident;
    identity *ident;
    int ni = sys->n_identities;
    int err = 0;

    ident = parse_identity(line, pdinfo, &err);
    if (ident == NULL) return err;

    /* connect the identity to the equation system */
    pident = realloc(sys->idents, (ni + 1) * sizeof *sys->idents);
    if (pident == NULL) {
	destroy_ident(ident);
	return E_ALLOC;
    }

    sys->idents = pident;
    sys->idents[ni] = ident;
    sys->n_identities += 1;

    return 0;
}

static int
add_aux_list_to_sys (gretl_equation_system *sys, const char *line,
		     const DATAINFO *pdinfo, int which)
{
    const char *p;
    char vname[VNAMELEN];
    int *list;
    int i, j, v, nf, len, cplen;
    int err = 0;

    if (which == ENDOG_LIST) {
	if (sys->endog_vars != NULL) {
	    strcpy(gretl_errmsg, "Only one list of endogenous variables may be given");
	    return 1;
	} 
    } else if (which == INSTR_LIST) {
	if (sys->instr_vars != NULL) {
	    strcpy(gretl_errmsg, "Only one list of instruments may be given");
	    return 1;
	} else if (0 && sys->method != SYS_3SLS && sys->method != SYS_TSLS) {
	    /* FIXME ?? */
	    strcpy(gretl_errmsg, "Instruments may be specified only "
		   "for 3SLS or TSLS");
	    return 1;
	}
    } else {
	return 1;
    }

    nf = count_fields(line);
    if (nf < 1) {
	return 1;
    }

    list = gretl_list_new(nf);
    if (list == NULL) {
	return E_ALLOC;
    }

    p = line;
    j = 1;
    for (i=1; i<=nf && !err; i++) {
	while (isspace(*p)) p++;
	*vname = '\0';
	cplen = len = strcspn(p, " \t\n");
	if (cplen > VNAMELEN - 1) {
	    cplen = VNAMELEN - 1;
	}
	strncat(vname, p, cplen);
	if (isdigit(*vname)) {
	    v = atoi(vname);
	} else {
	    v = varindex(pdinfo, vname);
	}
	if (v < 0 || v >= pdinfo->v) {
	    /* maybe a named list */
	    int *sublist = get_list_by_name(vname);

	    if (sublist == NULL) {
		sprintf(gretl_errmsg, _("Undefined variable '%s'."), vname);
		err = 1;
	    } else {
		err = gretl_list_insert_list_minus(&list, sublist, j);
		if (!err) {
		    j += sublist[0];
		}
	    }
	} else {
	    list[j++] = v;
	}
	p += len;
    }
	
    if (err) {
	free(list);
	return err;
    }

    if (which == ENDOG_LIST) {
	sys->endog_vars = list;
    } else {
	sys->instr_vars = list;
    }

    return 0;
}

/**
 * system_parse_line:
 * @sys: initialized equation system.
 * @line: command line.
 * @pdinfo: dataset information.
 * 
 * Modifies @sys according to the command supplied in @line,
 * which must start with "identity" (and supply an identity
 * to be added to the system, or "endog" (and supply a list
 * of endogenous variables), or "instr" (and supply a list of
 * instrumental variables).
 * 
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int 
system_parse_line (gretl_equation_system *sys, const char *line,
		   const DATAINFO *pdinfo)
{
    int err = 0;

    gretl_error_clear();

    if (strncmp(line, "identity", 8) == 0) {
	err = add_identity_to_sys(sys, line + 8, pdinfo);
    } else if (strncmp(line, "endog", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, ENDOG_LIST);
    } else if (strncmp(line, "instr", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, INSTR_LIST);
    } else {
	err = E_PARSE;
    }

    if (err) {
	gretl_equation_system_destroy(sys);
    }

    return err;
}

static int sys_in_list (const int *list, int k)
{
    int i;

    if (list == NULL) {
	return 0;
    }

    for (i=1; i<=list[0]; i++) {
	if (list[i] < 0) break;
	if (k == list[i]) return 1;
    }

    return 0;
}

static int make_instrument_list (gretl_equation_system *sys)
{
    int *ilist, *elist = sys->endog_vars;
    int i, j, k, nexo, maxnexo = 0;

    if (sys->instr_vars != NULL) {
	/* job is already done */
	return 0;
    }

    if (sys->method != SYS_SUR && 
	sys->method != SYS_OLS && 
	sys->method != SYS_WLS &&
	elist == NULL) {
	/* no list of endog vars: can't proceed */
	return 1;
    }

    /* First pass: get a count of the max possible number of
       exogenous variables (probably an over-estimate due to
       double-counting).
    */

    for (i=0; i<sys->n_equations; i++) {
	const int *slist = sys->lists[i];

	for (j=2; j<=slist[0]; j++) {
	    if (!sys_in_list(elist, slist[j])) {
		maxnexo++;
	    }
	}
    }

    for (i=0; i<sys->n_identities; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    if (!sys_in_list(elist, ident->atoms[j].varnum)) {
		maxnexo++;
	    }
	}
    }

    ilist = malloc((maxnexo + 1) * sizeof *ilist);
    if (ilist == NULL) {
	return E_ALLOC;
    }

    ilist[0] = maxnexo;
    for (i=1; i<=maxnexo; i++) {
	ilist[i] = -1;
    }

    /* Form list of exogenous variables, drawing on both the
       stochastic equations and the identities, if any. 
    */

    nexo = 0;

    for (i=0; i<sys->n_equations; i++) {
	const int *slist = sys->lists[i];

	for (j=2; j<=slist[0]; j++) {
	    k = slist[j];
	    if (!sys_in_list(elist, k) && 
		!sys_in_list(ilist, k)) {
		ilist[++nexo] = k;
	    }
	}
    } 

    for (i=0; i<sys->n_identities; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    k = ident->atoms[j].varnum;
	    if (!sys_in_list(elist, k) &&
		!sys_in_list(ilist, k)) {
		ilist[++nexo] = k;
	    }
	}
    }

    /* Trim the list (remove effect of double-counting) */

    if (nexo < maxnexo) {
	int *plist = realloc(ilist, (nexo + 1) * sizeof *plist);

	if (plist != NULL) {
	    ilist = plist;
	}
	ilist[0] = nexo;
    }

    sys->instr_vars = ilist;

    return 0;
}

void 
system_set_restriction_matrices (gretl_equation_system *sys,
				 gretl_matrix *R, gretl_matrix *q)
{
    system_clear_restrictions(sys);

    sys->R = R;
    sys->q = q;

    sys->flags |= GRETL_SYS_RESTRICT;
}

double *
gretl_equation_system_get_series (const gretl_equation_system *sys, 
				  const DATAINFO *pdinfo,
				  int idx, const char *key, int *err)
{
    double *x = NULL;
    const char *msel;
    int t, col = 0;

    if (sys == NULL || idx != M_UHAT) { /* FIXME generalize this */
	*err = E_BADSTAT;
	return NULL;
    }

    msel = strchr(key, '[');
    if (msel == NULL || sscanf(msel, "[,%d]", &col) != 1) {
	*err = E_PARSE;
    } else if (col <= 0 || col > sys->n_equations) {
	*err = E_DATA;
    }

    if (!*err) {
	x = malloc(pdinfo->n * sizeof *x);
	if (x == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	if (sys->uhat == NULL) {
	    *err = E_DATA;
	    free(x);
	    x = NULL;
	} else {
	    int s = 0;

	    col--;
	    for (t=0; t<pdinfo->n; t++) {
		if (t < sys->t1 || t > sys->t2) {
		    x[t] = NADBL;
		} else {
		    x[t] = gretl_matrix_get(sys->uhat, s++, col);
		}
	    }
	}
    }

    return x;
}

gretl_matrix *
gretl_equation_system_get_matrix (const gretl_equation_system *sys, int idx, 
				  int *err)
{
    gretl_matrix *M = NULL;

    if (sys == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    switch (idx) {  
    case M_COEFF:
	if (sys->b == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(sys->b);
	}
	break;
    case M_UHAT:
	M = gretl_matrix_copy(sys->uhat);
	break;
    case M_VCV:
	if (sys->vcv == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(sys->vcv);
	}
	break;
    case M_MSIGMA:
	M = gretl_matrix_copy(sys->sigma);
	break;
    default:
	*err = E_BADSTAT;
	break;
    }

    if (M == NULL && !*err) {
	*err = E_ALLOC;
    }

    return M;
}

int highest_numbered_var_in_system (const gretl_equation_system *sys, 
				    const DATAINFO *pdinfo)
{
    int i, j, v, vmax = 0;

    for (i=0; i<sys->n_equations; i++) {
	for (j=1; j<=sys->lists[i][0]; j++) {
	    v = sys->lists[i][j];
	    if (v == LISTSEP || v >= pdinfo->v) {
		/* temporary variables, already gone? */
		continue;
	    }
	    if (v > vmax) {
		vmax = v;
	    }
	}
    }

    return vmax;
}

static identity *
sys_retrieve_identity (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    identity *ident;
    xmlNodePtr cur;
    int n_atoms, depvar;
    int i, got = 0;

    got += gretl_xml_get_prop_as_int(node, "n_atoms", &n_atoms);
    got += gretl_xml_get_prop_as_int(node, "depvar", &depvar);
    if (got < 2) {
	*err = E_DATA;
	return NULL;
    }

    ident = ident_new(n_atoms);
    if (ident == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    ident->depvar = depvar;

    cur = node->xmlChildrenNode;

    i = 0;
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "id_atom")) {
	    got = gretl_xml_get_prop_as_int(cur, "op", &ident->atoms[i].op);
	    got += gretl_xml_get_prop_as_int(cur, "varnum", &ident->atoms[i].varnum);
	    if (got < 2) {
		*err = E_DATA;
	    } else {
		i++;
	    }
	}
	cur = cur->next;
    }

    if (!*err && i != n_atoms) {
	*err = E_DATA;
    }

    if (*err) {
	destroy_ident(ident);
	ident = NULL;
    }

    return ident;
}

gretl_equation_system *
gretl_system_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    gretl_equation_system *sys;
    xmlNodePtr cur;
    char *name;
    int method = 0;
    int i, j, got = 0;

    got += gretl_xml_get_prop_as_string(node, "name", &name);
    got += gretl_xml_get_prop_as_int(node, "method", &method);

    if (got < 2) {
	*err = E_DATA;
	return NULL;
    }     

    sys = gretl_equation_system_new(method, name);
    if (sys == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    got = 0;
    got += gretl_xml_get_prop_as_int(node, "n_equations", &sys->n_equations);
    got += gretl_xml_get_prop_as_int(node, "n_identities", &sys->n_identities);
    got += gretl_xml_get_prop_as_char(node, "flags", &sys->flags);

    if (got < 3) {
	*err = E_DATA;
	goto bailout;
    } 

    sys->lists = malloc(sys->n_equations * sizeof sys->lists);
    if (sys->lists == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (sys->n_identities > 0) {
	sys->idents = malloc(sys->n_identities * sizeof sys->idents);
	if (sys->idents == NULL) {
	    *err = E_ALLOC;
	    goto bailout;
	}
    }

    cur = node->xmlChildrenNode;

    i = j = 0;
    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "eqnlist")) {
	    sys->lists[i++] = gretl_xml_node_get_list(cur, doc, err); 
	} else if (!xmlStrcmp(cur->name, (XUC) "endog_vars")) {
	    sys->endog_vars = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "instr_vars")) {
	    sys->instr_vars = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "identity")) {
	    sys->idents[j++] = sys_retrieve_identity(cur, doc, err); 
	} else if (!xmlStrcmp(cur->name, (XUC) "R")) {
	    sys->R = gretl_xml_get_matrix(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "q")) {
	    sys->q = gretl_xml_get_matrix(cur, doc, err);
	}
	cur = cur->next;
    } 

    if (!*err && (i != sys->n_equations || j != sys->n_identities)) {
	*err = E_DATA;
	gretl_equation_system_destroy(sys);
	sys = NULL;
    }

 bailout:

    return sys;
}

static void xml_print_identity (identity *ident, FILE *fp)
{
    int i;

    fprintf(fp, "<identity n_atoms=\"%d\" depvar=\"%d\">\n",
	    ident->n_atoms, ident->depvar);
    for (i=0; i<ident->n_atoms; i++) {
	fprintf(fp, " <id_atom op=\"%d\" varnum=\"%d\"/>\n",
		ident->atoms[i].op, ident->atoms[i].varnum);
    }
    fputs("</identity>\n", fp);
}

int gretl_system_serialize (gretl_equation_system *sys, SavedObjectFlags flags,
			    FILE *fp)
{
    int i, err = 0;

    fprintf(fp, "<gretl-equation-system name=\"%s\" saveflags=\"%d\" method=\"%d\" ",  
	    (sys->name != NULL)? sys->name : "none", flags, sys->method);

    fprintf(fp, "n_equations=\"%d\" n_identities=\"%d\" flags=\"%d\">\n",
	    sys->n_equations, sys->n_identities, (int) sys->flags);

    for (i=0; i<sys->n_equations; i++) {
	gretl_xml_put_tagged_list("eqnlist", sys->lists[i], fp);
    }

    gretl_xml_put_tagged_list("endog_vars", sys->endog_vars, fp);
    gretl_xml_put_tagged_list("instr_vars", sys->instr_vars, fp);

    for (i=0; i<sys->n_identities; i++) {
	xml_print_identity(sys->idents[i], fp);
    }    

    gretl_xml_put_matrix(sys->R, "R", fp);
    gretl_xml_put_matrix(sys->q, "q", fp);

    fputs("</gretl-equation-system>\n", fp);

    return err;
}

void gretl_system_set_save_flag (gretl_equation_system *sys)
{
    sys->flags |= GRETL_SYS_SAVEIT;
}

void gretl_system_unset_save_flag (gretl_equation_system *sys)
{
    sys->flags &= ~GRETL_SYS_SAVEIT;
}

int gretl_system_save_flag_set (gretl_equation_system *sys)
{
    if (sys == NULL) {
	return 0;
    } else if (sys->flags & GRETL_SYS_SAVEIT) {
	return 1;
    } else {
	return 0;
    }
}
