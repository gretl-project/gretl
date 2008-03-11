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

#define SYSDEBUG 0

enum {
    OP_PLUS,
    OP_MINUS
} identity_ops;

enum {
    ENDOG_LIST,
    INSTR_LIST
} aux_list_types;

enum {
    SYS_TEST_NONE,
    SYS_TEST_LR,
    SYS_TEST_F,
    SYS_TEST_NOTIMP
} system_test_types;

struct id_atom_ {
    int op;         /* operator (plus or minus) */
    int varnum;     /* ID number of variable to right of operator */
};

struct identity_ {
    int n_atoms;    /* number of "atomic" elements in identity */
    int depvar;     /* LHS variable in indentity */
    id_atom *atoms; /* pointer to RHS "atoms" */
};

struct predet_ {
    int id;         /* ID number of lag variable */
    int src;        /* ID number of "parent" */
    int lag;        /* lag order */
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
static int make_instrument_list (equation_system *sys,
				 const DATAINFO *pdinfo);

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

static int instrument_is_predet (const equation_system *sys, int v)
{
    int i;

    for (i=0; i<sys->n_predet; i++) {
	if (v == sys->pre_vars[i].id) {
	    return 1;
	}
    }

    return 0;
}

static int sys_max_predet_lag (const equation_system *sys)
{
    int i, m = 0;

    for (i=0; i<sys->n_predet; i++) {
	if (sys->pre_vars[i].lag > m) {
	    m = sys->pre_vars[i].lag;
	}
    }

    return m;
}

static int get_predet_parent (const equation_system *sys, int v, int *lag)
{
    int i;

    for (i=0; i<sys->n_predet; i++) {
	if (v == sys->pre_vars[i].id) {
	    *lag = sys->pre_vars[i].lag;
	    return sys->pre_vars[i].src;
	}
    }

    *lag = 0;

    return -1;
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
print_equation_system_info (const equation_system *sys, 
			    const DATAINFO *pdinfo, 
			    gretlopt opt, PRN *prn)
{
    int header = (opt & OPT_H);
    int i, vi, lag;
    
    if (header && sys->name != NULL) {
	pprintf(prn, "Equation system %s\n", sys->name);
    }

    if (!header) {
	for (i=0; i<sys->neqns; i++) {
	    print_system_equation(sys->lists[i], pdinfo, prn);
	}    
    }

    for (i=0; i<sys->nidents; i++) {
	print_system_identity(sys->idents[i], pdinfo, opt, prn);
    }

    if (sys->ylist != NULL) {
	pputs(prn, (header)? "Endogenous variables:" : "endog");
	for (i=1; i<=sys->ylist[0]; i++) {
	    vi = sys->ylist[i];
	    pprintf(prn, " %s", pdinfo->varname[vi]);
	}
	pputc(prn, '\n');
    }

    if (header) {
	if (sys->pre_vars != NULL) {
	    pputs(prn, "Predetermined variables:");
	    for (i=0; i<sys->n_predet; i++) {
		vi = sys->pre_vars[i].src;
		lag = sys->pre_vars[i].lag;
		pprintf(prn, " %s(-%d)", pdinfo->varname[vi], lag);
	    }
	    pputc(prn, '\n');
	}    
	if (sys->ilist != NULL && sys->ilist[0] > sys->n_predet) {
	    pputs(prn, "Exogenous variables:");
	    for (i=1; i<=sys->ilist[0]; i++) {
		vi = sys->ilist[i];
		if (!instrument_is_predet(sys, vi)) {
		    pprintf(prn, " %s", pdinfo->varname[vi]);
		}
	    }
	    pputc(prn, '\n');
	}
    } else if (sys->ilist != NULL) {
	pputs(prn, "instr");
	for (i=1; i<=sys->ilist[0]; i++) {
	    vi = sys->ilist[i];
	    pprintf(prn, " %s", pdinfo->varname[vi]);
	}
	pputc(prn, '\n');
    }	
}

int system_method_from_string (const char *s)
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
    if (method >= SYS_METHOD_SUR && method < SYS_METHOD_MAX) {
	return gretl_system_long_strings[method];
    } else {
	return NULL;
    }
}

const char *system_method_short_string (int method)
{
    if (method >= SYS_METHOD_SUR && method < SYS_METHOD_MAX) {
	return gretl_system_method_strings[method];
    } else {
	return NULL;
    }
}

static equation_system *
equation_system_new (int method, const char *name)
{
    equation_system *sys;

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
    sys->T = sys->df = 0;

    sys->neqns = 0;
    sys->nidents = 0;

    sys->R = NULL;
    sys->q = NULL;

    sys->iters = 0;
    sys->flags = 0;
    sys->n_predet = 0;
    sys->order = 0;

    sys->ll = sys->llu = sys->ldet = NADBL;
    sys->X2 = NADBL;
    sys->ess = NADBL;
    sys->diag = 0.0;
    sys->bdiff = 0.0;

    sys->b = NULL;
    sys->vcv = NULL;
    sys->S = NULL;
    sys->E = NULL;
    sys->yhat = NULL;

    sys->Gamma = NULL;
    sys->B = NULL;
    sys->A = NULL;
    sys->Sr = NULL;
    sys->F = NULL;

    sys->lists = NULL;
    sys->ylist = NULL;
    sys->ilist = NULL;
    sys->xlist = NULL;
    sys->pre_vars = NULL;
    sys->idents = NULL;

    sys->models = NULL;

    return sys;
}

static void system_clear_restrictions (equation_system *sys)
{
    if (sys->R != NULL) {
	free(sys->R);
	sys->R = NULL;
    }

    if (sys->q != NULL) {
	free(sys->q);
	sys->q = NULL;
    }

    sys->flags &= ~SYSTEM_RESTRICT;
}

static void system_clear_results (equation_system *sys)
{
    sys->iters = 0;

    sys->t1 = sys->t2 = 0;
    sys->T = sys->df = 0;

    sys->ll = NADBL;
    sys->llu = NADBL;
    sys->ldet = NADBL;
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

    if (sys->S != NULL) {
	gretl_matrix_free(sys->S);
	sys->S = NULL;
    }

    if (sys->E != NULL) {
	gretl_matrix_free(sys->E);
	sys->E = NULL;
    }

    if (sys->yhat != NULL) {
	gretl_matrix_free(sys->yhat);
	sys->yhat = NULL;
    }

    if (sys->Gamma != NULL) {
	gretl_matrix_free(sys->Gamma);
	sys->Gamma = NULL;
    }

    if (sys->B != NULL) {
	gretl_matrix_free(sys->B);
	sys->B = NULL;
    }

    if (sys->A != NULL) {
	gretl_matrix_free(sys->A);
	sys->A = NULL;
    }

    if (sys->Sr != NULL) {
	gretl_matrix_free(sys->Sr);
	sys->Sr = NULL;
    }    

    if (sys->F != NULL) {
	gretl_matrix_free(sys->F);
	sys->F = NULL;
    }
}

void equation_system_destroy (equation_system *sys)
{
    int i;

    if (sys == NULL || sys->lists == NULL) {
	return;
    }

    sys->refcount -= 1;
    if (sys->refcount > 0) {
	return;
    }

    for (i=0; i<sys->neqns; i++) {
	free(sys->lists[i]);
    }

    free(sys->lists);
    sys->lists = NULL;

    for (i=0; i<sys->nidents; i++) {
	destroy_ident(sys->idents[i]);
    }

    free(sys->idents);

    free(sys->ylist);
    free(sys->ilist);
    free(sys->xlist);
    free(sys->pre_vars);

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

static void sur_rearrange_lists (equation_system *sys,
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    if (sys->method == SYS_METHOD_SUR) {
	int i;

	for (i=0; i<sys->neqns; i++) {
	    reglist_check_for_const(sys->lists[i], Z, pdinfo);
	}
    }
}

/**
 * equation_system_append:
 * @sys: initialized equation system.
 * @list: list containing dependent variable and regressors.
 * 
 * Adds an equation (as represented by @list) to @sys.
 * 
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int equation_system_append (equation_system *sys, 
			    const int *list)
{
    int n;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return E_DATA;
    }

    n = sys->neqns;

    sys->lists = realloc(sys->lists, (n + 1) * sizeof *sys->lists);
    if (sys->lists == NULL) {
	return E_ALLOC;
    }

    sys->lists[n] = gretl_list_copy(list);

    if (sys->lists[n] == NULL) {
	equation_system_destroy(sys);
	return E_ALLOC;
    }

    sys->neqns += 1;

    return 0;
}

/* retrieve the name -- possible quoted with embedded spaces -- for
   an equation system */

char *get_system_name_from_line (const char *s, int context)
{
    char *tests[] = {
	" name",
	"estimate ",
	"restrict "
    };
    const char *p = NULL;
    char *name = NULL;
    int pchars = 0;

    if (context < 0 || context > 3) {
	return NULL;
    }

    p = strstr(s, tests[context]);

    if (context == SYSNAME_NEW && p == NULL) {
	char savename[MAXSAVENAME];

	gretl_cmd_get_savename(savename);
	if (*savename != '\0') {
	    return gretl_strdup(savename);
	} 
    } 

    if (p == NULL) {
	return NULL;
    }

    p += strlen(tests[context]);

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
	method = system_method_from_string(mstr);
    }

    return method;
}

static void 
system_set_flags (equation_system *sys, const char *s, gretlopt opt)
{
    if (opt & OPT_T) {
	sys->flags |= SYSTEM_ITERATE;
    }

    s = strstr(s, " save");

    if (s != NULL) {
	s += 5;
	if (*s != ' ' && *s != '=') {
	    return;
	}
	if (strstr(s, "resids") || strstr(s, "uhat")) {
	    sys->flags |= SYSTEM_SAVE_UHAT;
	}
	if (strstr(s, "fitted") || strstr(s, "yhat")) {
	    sys->flags |= SYSTEM_SAVE_YHAT;
	}
    }
}

/**
 * equation_system_start:
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

equation_system *equation_system_start (const char *line, 
					gretlopt opt)
{
    equation_system *sys = NULL;
    char *sysname = NULL;
    int method;

    method = get_estimation_method_from_line(line);

    if (method == SYS_METHOD_MAX) {
	/* invalid method was given */
	strcpy(gretl_errmsg, _(badsystem));
	return NULL;
    }

    sysname = get_system_name_from_line(line, SYSNAME_NEW);

    if (method < 0 && sysname == NULL) {
	/* neither a method nor a name was specified */
	strcpy(gretl_errmsg, _(badsystem));
	return NULL;
    }

    sys = equation_system_new(method, sysname);
    if (sys == NULL) {
	return NULL;
    }

    system_set_flags(sys, line, opt);

    if (sysname != NULL) {
	free(sysname);
    }

    return sys;
}

static int system_get_dfu (const equation_system *sys)
{
    int dfu = sys->T * sys->neqns;
    int i;

    for (i=0; i<sys->neqns; i++) {
	dfu -= sys->lists[i][0] - 1;
    }

    return dfu;
}

/* Asymptotic F-test, as in Greene:

   (1/J) * (Rb-q)' * [R*Var(b)*R']^{-1} * (Rb-q) 
*/

static void 
system_print_F_test (const equation_system *sys,
		     const gretl_matrix *b, 
		     const gretl_matrix *vcv,
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
system_print_LR_test (const equation_system *sys,
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

static int sys_test_type (equation_system *sys)
{
    int ret = SYS_TEST_NONE;

    if (system_n_restrictions(sys) > 0) {
	if (sys->method == SYS_METHOD_SUR || 
	    sys->method == SYS_METHOD_WLS) {
	    if (sys->flags & SYSTEM_ITERATE) {
		ret = SYS_TEST_LR;
	    } else {
		ret = SYS_TEST_F;
	    }
	} else if (sys->method == SYS_METHOD_OLS ||
		   sys->method == SYS_METHOD_TSLS ||
		   sys->method == SYS_METHOD_3SLS) {
	    ret = SYS_TEST_F;
	} else if (sys->method == SYS_METHOD_LIML) {
	    /* experimental */
	    ret = SYS_TEST_F;
	} else if (sys->method == SYS_METHOD_FIML) {
	    ret = SYS_TEST_LR;
	} 
    }

    return ret;
}

/* Handle the case where sys->b and sys->vcv have been augmented in
   order to test some restrictions: we now have to trim these matrices
   down to their final size. The vector 'b' contains unrestricted
   coefficient estimates, so its length can be used to gauge the
   true number of coeffs.
*/

static int shrink_b_and_vcv (const gretl_matrix *b,
			     equation_system *sys)
{
    int nc = gretl_vector_get_length(b);
    gretl_matrix *V;
    double x;
    int i, j;

    if (sys->vcv->rows == nc) {
	/* no-op (shouldn't happen) */
	return 0;
    }

    V = gretl_matrix_alloc(nc, nc);
    if (V == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_reuse(sys->b, nc, 1);

    for (i=0; i<nc; i++) {
	for (j=0; j<nc; j++) {
	    x = gretl_matrix_get(sys->vcv, i, j);
	    gretl_matrix_set(V, i, j, x);
	}
    }

    gretl_matrix_free(sys->vcv);
    sys->vcv = V;

    return 0;
}

static int estimate_with_test (equation_system *sys, 
			       double ***pZ, DATAINFO *pdinfo, 
			       gretlopt opt, int stest, 
			       int (*system_est)(), PRN *prn)
{
    gretl_matrix *vcv = NULL;
    gretl_matrix *b = NULL;
    double llu = 0.0;
    int err = 0;

    /* estimate the unrestricted system first */

    sys->flags &= ~SYSTEM_RESTRICT;

    err = (* system_est) (sys, pZ, pdinfo, opt | OPT_Q, prn);

    sys->flags ^= SYSTEM_RESTRICT;

    if (err) {
	goto bailout;
    }

    /* grab the data from unrestricted estimation */

    b = sys->b;
    sys->b = NULL;

    if (stest == SYS_TEST_LR) {
	llu = sys->ll;
    } else if (stest == SYS_TEST_F) {
	vcv = sys->vcv;
	sys->vcv = NULL;
    }

    /* now estimate the restricted system */

    system_clear_results(sys);
    err = (* system_est) (sys, pZ, pdinfo, opt, prn);

    if (!err) {
	if (stest == SYS_TEST_LR) {
	    system_print_LR_test(sys, llu, prn);
	} else if (stest == SYS_TEST_F) {
	    system_print_F_test(sys, b, vcv, prn);
	}
	err = shrink_b_and_vcv(b, sys);
    }

 bailout:

    if (b != NULL) gretl_matrix_free(b);
    if (vcv != NULL) gretl_matrix_free(vcv);

    return err;
}

static void 
adjust_sys_flags_for_method (equation_system *sys, int method)
{
    char oldflags = sys->flags;

    sys->flags = 0;

    /* the iterate option is only available for WLS, SUR or 3SLS */

    if (oldflags & SYSTEM_ITERATE) {
	if (sys->method == SYS_METHOD_WLS || 
	    sys->method == SYS_METHOD_SUR ||
	    sys->method == SYS_METHOD_3SLS) {
	    sys->flags |= SYSTEM_ITERATE;
	}
    }

    /* by default, apply a df correction for single-equation methods */

    if (sys->method == SYS_METHOD_OLS || 
	sys->method == SYS_METHOD_WLS ||
	sys->method == SYS_METHOD_TSLS || 
	sys->method == SYS_METHOD_LIML) {
	if (oldflags & SYSTEM_DFCORR) {
	    sys->flags |= SYSTEM_DFCORR;
	}
    } 

    /* carry forward the GEOMEAN flag */

    if (oldflags & SYSTEM_VCV_GEOMEAN) {
	sys->flags |= SYSTEM_VCV_GEOMEAN;
    }    
}

static void 
set_sys_flags_from_opt (equation_system *sys, gretlopt opt)
{
    char oldflags = sys->flags;

    sys->flags = 0;

    /* the iterate option is available for WLS, SUR or 3SLS */

    if (opt & OPT_T) {
	if (sys->method == SYS_METHOD_WLS || 
	    sys->method == SYS_METHOD_SUR ||
	    sys->method == SYS_METHOD_3SLS) {
	    sys->flags |= SYSTEM_ITERATE;
	}
    }

    /* by default, apply a df correction for single-equation methods */

    if (sys->method == SYS_METHOD_OLS || 
	sys->method == SYS_METHOD_WLS ||
	sys->method == SYS_METHOD_TSLS || 
	sys->method == SYS_METHOD_LIML) {
	if (!(opt & OPT_N)) {
	    sys->flags |= SYSTEM_DFCORR;
	}
    } 

    if (opt & OPT_M) {
	sys->flags |= SYSTEM_VCV_GEOMEAN;
    } 

    if (oldflags & SYSTEM_SAVE_UHAT) {
	sys->flags |= SYSTEM_SAVE_UHAT;
    }

    if (oldflags & SYSTEM_SAVE_YHAT) {
	sys->flags |= SYSTEM_SAVE_YHAT;
    }
}

/**
 * equation_system_estimate:
 * @sys: pre-defined equation system.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: may include %OPT_V for more verbose operation.
 * @prn: printing struct.
 * 
 * Estimate a pre-defined equation system and print the results
 * to @prn.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int 
equation_system_estimate (equation_system *sys, 
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, PRN *prn)
{
    void *handle = NULL;
    int (*system_est) (equation_system *, 
		       double ***, DATAINFO *, gretlopt, PRN *);
    int stest = 0;
    int err = 0;

    gretl_error_clear();

    if (opt == OPT_UNSET) {
	opt = OPT_NONE;
	adjust_sys_flags_for_method(sys, sys->method);
    } else {
	set_sys_flags_from_opt(sys, opt);
    }

    /* if restrictions are in place, reset the restricted flag */
    if (sys->R != NULL && sys->q != NULL) {
	sys->flags |= SYSTEM_RESTRICT;
    } 

    /* in case we're re-estimating */
    system_clear_results(sys);

    err = make_instrument_list(sys, pdinfo);
    if (err) {
	goto system_bailout;
    }

    if (sys->method == SYS_METHOD_SUR) {
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
    } else if (stest != SYS_TEST_NONE) {
	err = estimate_with_test(sys, pZ, pdinfo, opt, stest, 
				 system_est, prn);
    } else {
	err = (*system_est) (sys, pZ, pdinfo, opt, prn);
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
 * equation_system_finalize:
 * @sys: pre-defined equation system.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: may include %OPT_V for verbose operaton.
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

int equation_system_finalize (equation_system *sys, 
			      double ***pZ, DATAINFO *pdinfo,
			      gretlopt opt, PRN *prn)
{
    int err = 0;

    gretl_error_clear();

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    if (sys->neqns < 2) {
	strcpy(gretl_errmsg, _(toofew));
	equation_system_destroy(sys);
	return 1;
    }

    if (sys->method >= SYS_METHOD_MAX) {
	strcpy(gretl_errmsg, _(badsystem));
	equation_system_destroy(sys);
	return 1;
    }    

    if (sys->name != NULL) {
	/* save the system for subsequent estimation */
	err = gretl_stack_object_as(sys, GRETL_OBJ_SYS, sys->name);
    }

    if (!err && sys->method >= 0) {
	err = equation_system_estimate(sys, pZ, pdinfo, opt, prn);
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
    equation_system *sys;
    char *sysname;
    int method;

    if (strlen(line) < 12) {
	strcpy(gretl_errmsg, "estimate: no system name was provided");
	return 1;
    }

    sysname = get_system_name_from_line(line, SYSNAME_EST);
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
    if (method < 0 || method >= SYS_METHOD_MAX) {
	method = sys->method;
    }

    if (method < 0 || method >= SYS_METHOD_MAX) {
	strcpy(gretl_errmsg, "estimate: no valid method was specified");
	return 1;
    }

    sys->method = method;

    return equation_system_estimate(sys, pZ, pdinfo, opt, prn);
}

/* effective list length, allowance made for TSLS-style lists */

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

/* below: used in formulating restriction matrices for
   an equation system */

int system_get_list_length (const equation_system *sys, int i)
{
    if (i >= 0 && i < sys->neqns) {
	return get_real_list_length(sys->lists[i]);
    } else {
	return 0;
    }
}

int system_max_indep_vars (const equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->neqns; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	if (nvi > nv) nv = nvi;
    }

    return nv;
}

int system_n_indep_vars (const equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->neqns; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	nv += nvi;
    }

    return nv;
}

const char *system_short_string (const MODEL *pmod)
{
    int i = gretl_model_get_int(pmod, "method");

    return gretl_system_short_strings[i];
}

void equation_system_set_name (equation_system *sys, const char *name)
{
    if (sys->name != NULL && !strcmp(sys->name, name)) {
	return;
    }

    if (sys->name != NULL) {
	free(sys->name);
    }

    sys->name = gretl_strdup(name);
}

int system_adjust_t1t2 (equation_system *sys,
			int *t1, int *t2, const double **Z)
{
    int i, err = 0;

    for (i=0; i<sys->neqns && !err; i++) {
	err = check_for_missing_obs(sys->lists[i], t1, t2, Z, NULL);
    }

    if (!err) {
	sys->t1 = *t1;
	sys->t2 = *t2;
	sys->T = *t2 - *t1 + 1;
    }

    return err;
}

int *compose_tsls_list (equation_system *sys, int i)
{
    int *list;
    int j, k1, k2;

    if (i >= sys->neqns) {
	return NULL;
    }

    if (sys->ilist == NULL && make_instrument_list(sys, NULL)) {
	return NULL;
    }

    k1 = sys->lists[i][0];
    k2 = sys->ilist[0];

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
	    list[j] = sys->ilist[j - (k1 + 1)];
	}
    }

    return list;
}

int system_normality_test (const equation_system *sys, PRN *prn)
{
    int err = 0;

    if (sys->E == NULL || sys->S == NULL) {
	err = 1;
    } else {
	err = multivariate_normality_test(sys->E, 
					  sys->S, 
					  prn);
    }

    return err;
}

static void
add_system_var_info (equation_system *sys, int i, 
		     DATAINFO *pdinfo, int v, int code)
{
    char *label = VARLABEL(pdinfo, v);

    if (code == SYSTEM_SAVE_UHAT) {
	sprintf(pdinfo->varname[v], "uhat_s%02d", i);
	if (sys->method == SYS_METHOD_SUR) {
	    sprintf(label, _("SUR residual, equation %d"), i);
	} else if (sys->method == SYS_METHOD_3SLS) {
	    sprintf(label, _("3SLS residual, equation %d"), i);
	} else {
	    sprintf(label, _("system residual, equation %d"), i);
	}
    } else if (code == SYSTEM_SAVE_YHAT) {
	sprintf(pdinfo->varname[v], "yhat_s%02d", i);
	if (sys->method == SYS_METHOD_SUR) {
	    sprintf(label, _("SUR fitted value, equation %d"), i);
	} else if (sys->method == SYS_METHOD_3SLS) {
	    sprintf(label, _("3SLS fitted value, equation %d"), i);
	} else {
	    sprintf(label, _("system fitted value, equation %d"), i);
	}
    }
}

int 
system_add_resids_to_dataset (equation_system *sys, 
			      int eqnum, double ***pZ, 
			      DATAINFO *pdinfo)
{
    int v, t;

    if (sys->E == NULL) {
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
	    (*pZ)[v][t] = gretl_matrix_get(sys->E, t - sys->t1, eqnum);
	}
    }

    add_system_var_info(sys, eqnum + 1, pdinfo, v, SYSTEM_SAVE_UHAT);

    return 0;
}

/* simple accessor functions */

const char *system_get_full_string (const equation_system *sys)
{
    if (sys->flags & SYSTEM_ITERATE) {
	static char sysstr[64];

	sprintf(sysstr, _("iterated %s"), gretl_system_long_strings[sys->method]);
	return sysstr;
    } else {
	return gretl_system_long_strings[sys->method];
    }
}

int system_want_df_corr (const equation_system *sys)
{
    return sys->flags & SYSTEM_DFCORR;
}

int system_n_restrictions (const equation_system *sys)
{
    int nr = 0;

    if (sys->R != NULL && (sys->flags & SYSTEM_RESTRICT)) {
	nr = gretl_matrix_rows(sys->R);
    }

    return nr;
}

int *system_get_list (const equation_system *sys, int i)
{
    if (i >= sys->neqns) return NULL;

    return sys->lists[i];
}

int system_get_depvar (const equation_system *sys, int i)
{
    if (i >= sys->neqns) return 0;

    return sys->lists[i][1];
}

int system_get_method (const equation_system *sys)
{
    return sys->method;
}

int *system_get_endog_vars (const equation_system *sys)
{
    return sys->ylist;
}

int *system_get_instr_vars (const equation_system *sys)
{
    return sys->ilist;
}

void system_attach_coeffs (equation_system *sys, gretl_matrix *b)
{
    if (sys->b != NULL) {
	gretl_matrix_free(sys->b);
    }

    sys->b = b;
}

void system_attach_vcv (equation_system *sys, gretl_matrix *vcv)
{
    if (sys->vcv != NULL) {
	gretl_matrix_free(sys->vcv);
    }

    sys->vcv = vcv;
    gretl_matrix_xtr_symmetric(sys->vcv);
}

void system_attach_sigma (equation_system *sys, gretl_matrix *S)
{
    if (sys->S != NULL) {
	gretl_matrix_free(sys->S);
    }

    sys->S = S;
}

void system_attach_uhat (equation_system *sys, gretl_matrix *E)
{
    if (sys->E != NULL) {
	gretl_matrix_free(sys->E);
    }

    sys->E = E;
}

MODEL *system_get_model (const equation_system *sys, int i)
{
    if (sys->models == NULL || i >= sys->neqns) {
	return NULL; 
    }   

    return sys->models[i];
}

/* for case of applying df correction to cross-equation 
   residual covariance matrix */

int system_vcv_geomean (const equation_system *sys)
{
    return sys->flags & SYSTEM_VCV_GEOMEAN;
}

static int sys_eqn_indep_coeffs (const equation_system *sys, int eq)
{
    int nc;

    if (eq >= sys->neqns) {
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
system_vcv_denom (const equation_system *sys, int i, int j)
{
    double den = sys->T;

    if ((sys->flags & SYSTEM_VCV_GEOMEAN) &&
	i < sys->neqns && j < sys->neqns) {
	int ki = sys_eqn_indep_coeffs(sys, i);

	if (j == i) {
	    den = sys->T - ki;
	} else {
	    int kj = sys_eqn_indep_coeffs(sys, j);

	    den = (sys->T - ki) * (sys->T - kj);
	    den = sqrt(den);
	}
    }

    return den;
}

/* for system over-identification test */

int system_get_overid_df (const equation_system *sys)
{
    int gl, i, k = 0;

    gl = sys->neqns * sys->ilist[0];

    for (i=0; i<sys->neqns; i++) {
	k += sys->lists[i][0] - 1;
    }

    return gl - k;
}

/* dealing with identities (FIML, LIML) */

int rhs_var_in_identity (const equation_system *sys, int lhsvar,
			 int rhsvar)
{
    const identity *ident;
    int i, j;

    for (i=0; i<sys->nidents; i++) {
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

/* add information regarding a predetermined regressor */

static int 
add_predet_to_sys (equation_system *sys, const DATAINFO *pdinfo,
		   int id, int src, int lag)
{
    int i, n = sys->n_predet;
    predet *pre;

    if (id < 0 || src < 0 || id >= pdinfo->v || src >= pdinfo->v) {
	/* something screwy */
	return E_DATA;
    }

    for (i=0; i<n; i++) {
	if (sys->pre_vars[i].id == id) {
	    /* already present */
	    return 0;
	}
    }

    pre = realloc(sys->pre_vars, (n + 1) * sizeof *pre);
    if (pre == NULL) {
	return E_ALLOC;
    }

    sys->pre_vars = pre;
    sys->pre_vars[n].id = id;
    sys->pre_vars[n].src = src;
    sys->pre_vars[n].lag = lag;
    sys->n_predet += 1;

    return 0;
}

static int 
add_identity_to_sys (equation_system *sys, const char *line,
		     const DATAINFO *pdinfo)
{
    identity **pident;
    identity *ident;
    int ni = sys->nidents;
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
    sys->nidents += 1;

    return 0;
}

static int
add_aux_list_to_sys (equation_system *sys, const char *line,
		     const DATAINFO *pdinfo, int which)
{
    const char *p;
    char vname[VNAMELEN];
    int *list;
    int i, j, v, nf, len, cplen;
    int err = 0;

    if (which == ENDOG_LIST) {
	if (sys->ylist != NULL) {
	    strcpy(gretl_errmsg, "Only one list of endogenous variables may be given");
	    return 1;
	} 
    } else if (which == INSTR_LIST) {
	if (sys->ilist != NULL) {
	    strcpy(gretl_errmsg, "Only one list of instruments may be given");
	    return 1;
	} else if (0 && sys->method != SYS_METHOD_3SLS && 
		   sys->method != SYS_METHOD_TSLS) {
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
	sys->ylist = list;
    } else {
	sys->ilist = list;
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
system_parse_line (equation_system *sys, const char *line,
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
	equation_system_destroy(sys);
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

/* Handle the case where we were not given an explicit list 
   of endogenous variables, but rather the individual equations
   were given in TSLS style: "y X ; Z".  This is supported
   for the TSLS and 3SLS methods.
*/

static int make_tsls_style_ylist (equation_system *sys)
{
    int *test, *ylist = NULL;
    int i, vi, sep = 0;
    int err = 0;

    for (i=0; i<sys->neqns; i++) {
	if (gretl_list_has_separator(sys->lists[i])) {
	    sep++;
	}
    }

    if (sep > 0) {
	ylist = gretl_null_list();
	if (ylist == NULL) {
	    err = E_ALLOC;
	}
	for (i=0; i<sys->neqns && !err; i++) {
	    vi = sys->lists[i][1];
	    if (!in_gretl_list(ylist, vi)) {
		test = gretl_list_append_term(&ylist, vi);
		if (test == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    }

    if (!err) {
	sys->ylist = ylist;
    } else {
	free(ylist);
    }

    return err;
}

static int infer_ylist_from_insts (equation_system *sys)
{
    const int *Z = sys->ilist;
    const int *slist;
    int *ylist = NULL;
    int i, j, k;
    int err = 0;

    for (i=0; i<sys->neqns && !err; i++) {
	slist = sys->lists[i];
	if (gretl_list_has_separator(slist)) {
	    continue;
	}
	for (j=1; j<=slist[0] && !err; j++) {
	    k = slist[j];
	    if (!in_gretl_list(Z, k)) {
		gretl_list_append_term(&ylist, k);
		if (ylist == NULL) {
		    err = E_ALLOC;
		}
	    }
	}
    }

    if (!err) {
	sys->ylist = ylist;
    } 

    return err;
}

static int predet_check (equation_system *sys, const DATAINFO *pdinfo)
{
    const int *ylist = sys->ylist;
    const int *ilist = sys->ilist;
    const char *vname;
    int id = 0, src = 0, lag = 0;
    int i, err = 0;

    if (ylist == NULL) {
	/* FIXME SUR */
	return 0;
    }

    for (i=1; i<=ylist[0] && !err; i++) {
	id = ylist[i];
	lag = pdinfo->varinfo[id]->lag;
	if (lag > 0) {
	    vname = pdinfo->varinfo[id]->parent;
	    src = varindex(pdinfo, vname);
	    err = add_predet_to_sys(sys, pdinfo, id, src, lag);
	}
    }

    for (i=1; i<=ilist[0] && !err; i++) {
	id = ilist[i];
	lag = pdinfo->varinfo[id]->lag;
	if (lag > 0) {
	    vname = pdinfo->varinfo[id]->parent;
	    src = varindex(pdinfo, vname);
	    if (in_gretl_list(ylist, src)) {
		err = add_predet_to_sys(sys, pdinfo, id, src, lag);
	    }
	}
    }

    return err;
}

static int sys_max_nexo (equation_system *sys)
{
    const int *slist;
    int i, j, vj, gotsep;
    int n = 0;

    for (i=0; i<sys->neqns; i++) {
	slist = sys->lists[i];
	gotsep = 0;
	for (j=2; j<=slist[0]; j++) {
	    vj = slist[j];
	    if (gotsep) {
		n++;
	    } else if (vj == LISTSEP) {
		gotsep = 1;
	    } else if (!sys_in_list(sys->ylist, vj)) {
		n++;
	    }
	}
    }

    for (i=0; i<sys->nidents; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    if (!sys_in_list(sys->ylist, ident->atoms[j].varnum)) {
		n++;
	    }
	}
    }

    return n;
}

#define system_needs_endog_list(s) (s->method != SYS_METHOD_SUR && \
				    s->method != SYS_METHOD_OLS && \
				    s->method != SYS_METHOD_WLS)

static int 
make_instrument_list (equation_system *sys, const DATAINFO *pdinfo)
{
    int *ilist, *ylist;
    int nexo, maxnexo;
    int i, j, vj;
    int err = 0;

    if (system_needs_endog_list(sys) && sys->ylist == NULL) {
	/* first pass: handle 3SLS? */
	err = make_tsls_style_ylist(sys);
    }

    if (!err && system_needs_endog_list(sys) && sys->ylist == NULL) {
	if (sys->ilist != NULL) {
	    err = infer_ylist_from_insts(sys);
	}
    }

    if (!err && system_needs_endog_list(sys) && sys->ylist == NULL) {
	gretl_errmsg_set(_("No list of endogenous variables was given"));
	err = E_DATA;
    }

    if (err) {
	return err;
    }

    if (sys->ilist != NULL) {
	/* job is already done? */
	if (pdinfo != NULL && sys->n_predet == 0) {
	    err = predet_check(sys, pdinfo);
	}
	return err;
    }

    ylist = sys->ylist;

    /* Get a count of the max possible number of exogenous variables
       (probably an over-estimate due to double-counting).
    */

    maxnexo = sys_max_nexo(sys);

    ilist = gretl_list_new(maxnexo);
    if (ilist == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=maxnexo; i++) {
	ilist[i] = -1;
    }

    /* Form list of exogenous variables, drawing on both the
       stochastic equations and the identities, if any. 
    */

    nexo = 0;

    for (i=0; i<sys->neqns; i++) {
	const int *slist = sys->lists[i];

	for (j=2; j<=slist[0]; j++) {
	    vj = slist[j];
	    if (vj != LISTSEP && 
		!sys_in_list(ylist, vj) && 
		!sys_in_list(ilist, vj)) {
		ilist[++nexo] = vj;
	    }
	}
    } 

    for (i=0; i<sys->nidents; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    vj = ident->atoms[j].varnum;
	    if (!sys_in_list(ylist, vj) &&
		!sys_in_list(ilist, vj)) {
		ilist[++nexo] = vj;
	    }
	}
    }

    /* Trim the list (remove effect of double-counting) */
    if (nexo < maxnexo) {
	ilist[0] = nexo;
    }

    sys->ilist = ilist;

    /* now check for predetermined regressors? */
    if (pdinfo != NULL) {
	err = predet_check(sys, pdinfo);
    }

    return err;
}

void 
system_set_restriction_matrices (equation_system *sys,
				 gretl_matrix *R, gretl_matrix *q)
{
    system_clear_restrictions(sys);

    sys->R = R;
    sys->q = q;

    sys->flags |= SYSTEM_RESTRICT;
}

double *
equation_system_get_series (const equation_system *sys, 
			    const DATAINFO *pdinfo,
			    int idx, const char *key, int *err)
{
    const gretl_matrix *M = NULL;
    double *x = NULL;
    const char *msel;
    int t, col = 0;

    if (sys == NULL || (idx != M_UHAT && idx != M_YHAT)) { 
	*err = E_BADSTAT;
	return NULL;
    }

    msel = strchr(key, '[');
    if (msel == NULL || sscanf(msel, "[,%d]", &col) != 1) {
	*err = E_PARSE;
    } else if (col <= 0 || col > sys->neqns) {
	*err = E_DATA;
    }

    if (!*err) {
	M = (idx == M_UHAT)? sys->E : sys->yhat;
	if (M == NULL) {
	    *err = E_DATA;
	}
    }

    if (!*err) {
	x = malloc(pdinfo->n * sizeof *x);
	if (x == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	int s = 0;

	col--; /* switch to 0-based */
	for (t=0; t<pdinfo->n; t++) {
	    if (t < sys->t1 || t > sys->t2) {
		x[t] = NADBL;
	    } else {
		x[t] = gretl_matrix_get(M, s++, col);
	    }
	}
    }

    return x;
}

static int system_add_yhat_matrix (equation_system *sys)
{
    double x;
    int avc, nc = 0;
    int i, s, t;

    sys->yhat = gretl_matrix_alloc(sys->T, sys->neqns);
    if (sys->yhat == NULL) {
	return E_ALLOC;
    }

    sys->yhat->t1 = sys->t1;
    sys->yhat->t2 = sys->t2;

    for (i=0; i<sys->neqns; i++) {
	s = 0;
	for (t=sys->t1; t<=sys->t2; t++) {
	    x = sys->models[i]->yhat[t];
	    gretl_matrix_set(sys->yhat, s++, i, x);
	}
	nc += sys->models[i]->ncoeff;
    }

    avc = ceil((double) nc / sys->neqns);
    sys->df = sys->T - avc;

    return 0;
}

/* retrieve a copy of a specified matrix from an equation
   system */

gretl_matrix *
equation_system_get_matrix (const equation_system *sys, int idx, 
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
	M = gretl_matrix_copy(sys->E);
	break;
    case M_YHAT:
	M = gretl_matrix_copy(sys->yhat);
	break;
    case M_VCV:
	if (sys->vcv == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(sys->vcv);
	}
	break;
    case M_SIGMA:
	M = gretl_matrix_copy(sys->S);
	break;
    case M_SYSGAM:
	M = gretl_matrix_copy(sys->Gamma);
	break;
    case M_SYSA:
	M = gretl_matrix_copy(sys->A);
	break;
    case M_SYSB:
	M = gretl_matrix_copy(sys->B);
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

int highest_numbered_var_in_system (const equation_system *sys, 
				    const DATAINFO *pdinfo)
{
    int i, j, v, vmax = 0;

    for (i=0; i<sys->neqns; i++) {
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

equation_system *
equation_system_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    equation_system *sys;
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

    sys = equation_system_new(method, name);
    if (sys == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    got = 0;
    got += gretl_xml_get_prop_as_int(node, "n_equations", &sys->neqns);
    got += gretl_xml_get_prop_as_int(node, "nidents", &sys->nidents);
    got += gretl_xml_get_prop_as_char(node, "flags", &sys->flags);

    if (got < 3) {
	*err = E_DATA;
	goto bailout;
    } 

    sys->lists = malloc(sys->neqns * sizeof sys->lists);
    if (sys->lists == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    if (sys->nidents > 0) {
	sys->idents = malloc(sys->nidents * sizeof sys->idents);
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
	    sys->ylist = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "instr_vars")) {
	    sys->ilist = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "exog_vars")) {
	    sys->xlist = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "identity")) {
	    sys->idents[j++] = sys_retrieve_identity(cur, doc, err); 
	} else if (!xmlStrcmp(cur->name, (XUC) "R")) {
	    sys->R = gretl_xml_get_matrix(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "q")) {
	    sys->q = gretl_xml_get_matrix(cur, doc, err);
	}
	cur = cur->next;
    } 

    if (!*err && (i != sys->neqns || j != sys->nidents)) {
	*err = E_DATA;
	equation_system_destroy(sys);
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

int equation_system_serialize (equation_system *sys, 
			       SavedObjectFlags flags,
			       FILE *fp)
{
    int i, err = 0;

    fprintf(fp, "<gretl-equation-system name=\"%s\" saveflags=\"%d\" method=\"%d\" ",  
	    (sys->name != NULL)? sys->name : "none", flags, sys->method);

    fprintf(fp, "n_equations=\"%d\" nidents=\"%d\" flags=\"%d\">\n",
	    sys->neqns, sys->nidents, (int) sys->flags);

    for (i=0; i<sys->neqns; i++) {
	gretl_xml_put_tagged_list("eqnlist", sys->lists[i], fp);
    }

    gretl_xml_put_tagged_list("endog_vars", sys->ylist, fp);
    gretl_xml_put_tagged_list("instr_vars", sys->ilist, fp);
    gretl_xml_put_tagged_list("exog_vars", sys->xlist, fp);

    for (i=0; i<sys->nidents; i++) {
	xml_print_identity(sys->idents[i], fp);
    }    

    gretl_xml_put_matrix(sys->R, "R", fp);
    gretl_xml_put_matrix(sys->q, "q", fp);

    fputs("</gretl-equation-system>\n", fp);

    return err;
}

void system_set_save_flag (equation_system *sys)
{
    sys->flags |= SYSTEM_SAVEIT;
}

void system_unset_save_flag (equation_system *sys)
{
    sys->flags &= ~SYSTEM_SAVEIT;
}

int system_save_flag_is_set (equation_system *sys)
{
    if (sys == NULL) {
	return 0;
    } else if (sys->flags & SYSTEM_SAVEIT) {
	return 1;
    } else {
	return 0;
    }
}

static int sur_ols_diag (equation_system *sys)
{
    double s2, ls2sum = 0.0;
    int i, err = 0;

    for (i=0; i<sys->neqns; i++) {
	s2 = gretl_model_get_double(sys->models[i], "ols_sigma_squared");
	if (na(s2)) {
	    err = 1;
	    break;
	}
	ls2sum += log(s2);
    }

    if (!err) {
	sys->diag = ls2sum;
    }

    return err;
}

static void
print_system_overid_test (const equation_system *sys,
			  PRN *prn)
{
    int df = system_get_overid_df(sys);

    if (sys->method == SYS_METHOD_FIML && df > 0) {
	double X2;

	if (na(sys->ll) || na(sys->llu) || 
	    sys->ll == 0.0 || sys->llu == 0.0) {
	    return;
	}

	X2 = 2.0 * (sys->llu - sys->ll);

	pprintf(prn, "%s:\n", _("LR over-identification test"));
	pprintf(prn, "  %s = %g\n", _("Restricted log-likelihood"), sys->ll);
	pprintf(prn, "  %s = %g\n", _("Unrestricted log-likelihood"), sys->llu);
	pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
		df, X2, _("with p-value"), chisq_cdf_comp(X2, df));
	pputc(prn, '\n');
    } else if ((sys->method == SYS_METHOD_3SLS || 
		sys->method == SYS_METHOD_SUR) && df > 0) {
	if (na(sys->X2) || sys->X2 <= 0.0) {
	    pputs(prn, _("Warning: the Hansen-Sargan over-identification test "
		  "failed.\nThis probably indicates that the estimation "
		  "problem is ill-conditioned.\n"));
	    return;
	}

	pprintf(prn, "%s:\n", _("Hansen-Sargan over-identification test"));
	pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
		df, sys->X2, _("with p-value"), chisq_cdf_comp(sys->X2, df));
	pputc(prn, '\n');
    }
}

int system_print_VCV (const equation_system *sys, PRN *prn)
{
    int k, df;

    if (sys->S == NULL) {
	return E_DATA;
    }

    k = sys->S->rows;
    df = k * (k - 1) / 2;

    print_contemp_covariance_matrix(sys->S, sys->ldet, prn);

    if (sys->method == SYS_METHOD_SUR && sys->iters > 0) {
	if (!na(sys->ldet) && sys->diag != 0.0) {
	    double lr = sys->T * (sys->diag - sys->ldet);

	    pprintf(prn, "%s:\n", _("LR test for diagonal covariance matrix"));
	    pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
		    df, lr, _("with p-value"), chisq_cdf_comp(lr, df));
	}
    } else {
	double lm = sys->diag;

	if (lm > 0) {
	    pprintf(prn, "%s:\n", 
		    _("Breusch-Pagan test for diagonal covariance matrix"));
	    pprintf(prn, "  %s(%d) = %g %s %g\n", _("Chi-square"),
		    df, lm, _("with p-value"), chisq_cdf_comp(lm, df));
	}
    }

    pputc(prn, '\n');

    return 0;
}

static int *make_ylist_from_models (equation_system *sys)
{
    int *ylist;
    int i, vi, pos;

    ylist = gretl_list_new(sys->neqns);
    if (ylist == NULL) {
	return NULL;
    }

    for (i=0; i<sys->neqns; i++) {
	vi = sys->models[i]->list[1];
	ylist[i+1] = vi;
	pos = in_gretl_list(sys->ilist, vi);
	if (pos > 0) {
	    gretl_list_delete_at_pos(sys->ilist, pos);
	}
    }

    return ylist;
}

enum {
    ENDOG,
    EXOG,
    PREDET
};

static int categorize_variable (int vnum, const equation_system *sys,
				const int *xlist, int *col,
				int *lag)
{
    int pos, ret = -1;

    *lag = 0;

    pos = in_gretl_list(sys->ylist, vnum);
    if (pos > 0) {
	ret = ENDOG;
    } else {
	pos = in_gretl_list(xlist, vnum);
	if (pos > 0) {
	    ret = EXOG;
	} else {
	    int v = get_predet_parent(sys, vnum, lag);

	    pos = in_gretl_list(sys->ylist, v);
	    if (pos > 0) {
		ret = PREDET;
	    }
	}
    }  

    if (pos > 0) {
	*col = pos - 1;
    }

    return ret;
}

/* Assemble list of exogenous variables for a system that
   does not already have one (e.g. SUR?).  We go through
   all the model lists and proceed by elimination: the
   variable is exogenous if it's not marked as endogenous
   and neither is it marked as predetermined (i.e. is a
   lag of an endogenous variable).
*/

static int *make_exogenous_list (equation_system *sys)
{
    const MODEL *pmod;
    const int *mlist;
    int *xlist, *test;
    int i, j, vj;

    xlist = gretl_null_list();
    if (xlist == NULL) {
	return NULL;
    }
    
    for (i=0; i<sys->neqns; i++) {
	pmod = sys->models[i];
	mlist = pmod->list;
	for (j=1; j<=mlist[0] && mlist[j]!=LISTSEP; j++) {
	    vj = mlist[j];
	    if (!in_gretl_list(sys->ylist, vj) &&
		!in_gretl_list(xlist, vj) &&
		!instrument_is_predet(sys, vj)) {
		test = gretl_list_append_term(&xlist, vj);
		if (test == NULL) {
		    free(xlist);
		    return NULL;
		}
	    }
	}
    }

    return xlist;
}

/* reduced-form error covariance matrix */

static int sys_add_RF_covariance_matrix (equation_system *sys)
{
    gretl_matrix *G = NULL, *S = NULL;
    int n = sys->neqns + sys->nidents;
    int err = 0;

    if (!gretl_is_identity_matrix(sys->Gamma)) {
	G = gretl_matrix_copy(sys->Gamma);
	if (G == NULL) {
	    return E_ALLOC;
	}
    } 

    if (sys->nidents > 0) {
	S = gretl_zero_matrix_new(n, n);
	if (S == NULL) {
	    gretl_matrix_free(G);
	    return E_ALLOC;
	}
	gretl_matrix_inscribe_matrix(S, sys->S, 0, 0, GRETL_MOD_NONE);
    }

    sys->Sr = gretl_matrix_alloc(n, n);
    if (sys->Sr == NULL) {
	gretl_matrix_free(G);
	gretl_matrix_free(S);
	return E_ALLOC;
    }

    if (G != NULL) {
	gretl_SVD_invert_matrix(G);
	gretl_matrix_qform(G, GRETL_MOD_NONE,
			   (S != NULL)? S : sys->S,
			   sys->Sr, GRETL_MOD_NONE);
	gretl_matrix_free(G);
    } else {
	gretl_matrix_copy_values(sys->Sr, (S != NULL)? S : sys->S);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->Sr, "G^{-1} S G^{-1}'");
#endif

    if (S != NULL) {
	gretl_matrix_free(S);
    }

    return err;
}

static int sys_add_structural_form (equation_system *sys,
				    const DATAINFO *pdinfo)
{
    const int *ylist = sys->ylist;
    const int *xlist = sys->xlist;
    int ne = sys->neqns;
    int ni = sys->nidents;
    int n = ne + ni;
    double x = 0.0;
    int type, col;
    int i, j, k, vj, lag;
    int err = 0;

    /* get all the required lists in order */

    if (ylist == NULL) {
	/* SUR and possibly some others? */
	sys->ylist = make_ylist_from_models(sys);
	if (sys->ylist == NULL) {
	    return E_ALLOC;
	}
	ylist = sys->ylist;
    }

    if (xlist == NULL) {
	sys->xlist = make_exogenous_list(sys);
	if (sys->xlist == NULL) {
	    return E_ALLOC;
	}
	xlist = sys->xlist;
    }

#if SYSDEBUG
    printlist(ylist, "endogenous vars");
    printlist(xlist, "exogenous vars");
#endif

    sys->order = sys_max_predet_lag(sys);

    /* allocate coefficient matrices */

    sys->Gamma = gretl_zero_matrix_new(n, n);
    if (sys->Gamma == NULL) {
	return E_ALLOC;
    }

    if (sys->order > 0) {
	sys->A = gretl_zero_matrix_new(n, sys->order * ylist[0]);
	if (sys->A == NULL) {
	    return E_ALLOC;
	}
    }

    if (xlist[0] > 0) {
	sys->B = gretl_zero_matrix_new(n, xlist[0]);
	if (sys->B == NULL) {
	    return E_ALLOC;
	}
    }

    /* process stochastic equations */

    k = 0;
    for (i=0; i<sys->neqns && !err; i++) {
	const MODEL *pmod = sys->models[i];
	const int *mlist = pmod->list;

	for (j=1; j<=mlist[0] && mlist[j]!=LISTSEP; j++) {
	    vj = mlist[j];
	    type = categorize_variable(vj, sys, xlist, &col, &lag);
	    x = (j > 1)? pmod->coeff[j-2] : 1.0;
	    if (type == ENDOG) {
		if (j == 1) {
		    gretl_matrix_set(sys->Gamma, i, col, 1.0);
		} else {
		    gretl_matrix_set(sys->Gamma, i, col, -x);
		}
	    } else if (type == EXOG) {
		gretl_matrix_set(sys->B, i, col, x);
	    } else if (type == PREDET) {
		col += n * (lag - 1);
		gretl_matrix_set(sys->A, i, col, x);
	    } else {
		fprintf(stderr, "sys_add_structural_form: i=%d, j=%d, vj=%d, type=%d\n",
			i, j, vj, type);
		printlist(mlist, "model list");
		err = E_DATA;
	    }
	}
    }

    /* process identities */

    for (i=0; i<sys->nidents && !err; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<=ident->n_atoms; j++) {
	    if (j == 0) {
		vj = ident->depvar;
		x = 1.0;
	    } else {
		vj = ident->atoms[j-1].varnum;
		x = (ident->atoms[j-1].op)? -1.0 : 1.0;
	    }
	    type = categorize_variable(vj, sys, xlist, &col, &lag);
	    if (type == ENDOG) {
		if (j == 0) {
		    gretl_matrix_set(sys->Gamma, ne+i, col, 1.0);
		} else {
		    gretl_matrix_set(sys->Gamma, ne+i, col, -x);
		}
	    } else if (type == EXOG) {
		gretl_matrix_set(sys->B, ne+i, col, x);
	    } else if (type == PREDET) {
		col += n * (lag - 1);
		gretl_matrix_set(sys->A, ne+i, col, x);
	    } else {
		err = E_DATA;
	    }
	}
    }

    if (!err) {
	err = sys_add_RF_covariance_matrix(sys);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->Gamma, "sys->Gamma");

    if (sys->A != NULL) {
	gretl_matrix_print(sys->A, "sys->A");
    } else {
	fputs("No lagged endogenous variables used as instruments\n", stderr);
    }

    if (sys->B != NULL) {
	gretl_matrix_print(sys->B, "sys->B");
    } else {
	fputs("No truly exogenous variables are present\n", stderr);
    }
#endif

    return err;
}

static gretl_matrix *sys_companion_matrix (equation_system *sys)
{
    gretl_matrix *C;
    int m = sys->A->rows;
    int n = sys->A->cols;

    if (m == n) {
	return sys->A;
    }

    C = gretl_zero_matrix_new(n, n);

    if (C != NULL) {
	gretl_matrix_inscribe_matrix(C, sys->A, 0, 0,
				     GRETL_MOD_NONE);
	gretl_matrix_inscribe_I(C, m, 0, n - m);
    }

#if SYSDEBUG   
    gretl_matrix_print(C, "system companion matrix");
#endif

    return C;
}

static gretl_matrix *
sys_get_fcast_se (equation_system *sys, int periods)
{
    int k = sys->A->cols;
    int n = sys->neqns + sys->nidents;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *V0 = NULL, *Vt = NULL;
    gretl_matrix *C = NULL, *se = NULL;
    double vti;
    int i, t, err = 0;

    if (periods <= 0) {
	fprintf(stderr, "Invalid number of periods\n");
	return NULL;
    }

    C = sys_companion_matrix(sys);
    if (C == NULL) {
	return NULL;
    }

    se = gretl_zero_matrix_new(periods, n);
    if (se == NULL) {
	if (C != sys->A) {
	    gretl_matrix_free(C);
	}
	return NULL;
    }
    
    Vt = gretl_matrix_alloc(k, k);
    V0 = gretl_zero_matrix_new(k, k);
    Tmp = gretl_matrix_alloc(k, k);

    if (Vt == NULL || V0 == NULL || Tmp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (t=0; t<periods; t++) {
	if (t == 0) {
	    /* initial variance */
	    gretl_matrix_inscribe_matrix(V0, sys->Sr, 0, 0, GRETL_MOD_NONE);
	    gretl_matrix_copy_values(Vt, V0);
	} else {
	    /* calculate further variances */
	    gretl_matrix_copy_values(Tmp, Vt);
	    gretl_matrix_qform(C, GRETL_MOD_NONE,
			       Tmp, Vt, GRETL_MOD_NONE);
	    gretl_matrix_add_to(Vt, V0);
	}

	for (i=0; i<n; i++) {
	    vti = gretl_matrix_get(Vt, i, i);
	    gretl_matrix_set(se, t, i, sqrt(vti));
	}
    }

 bailout:

    gretl_matrix_free(V0);
    gretl_matrix_free(Vt);
    gretl_matrix_free(Tmp);
    if (C != sys->A) {
	gretl_matrix_free(C);
    }
    
    if (err) {
	gretl_matrix_free(se);
	se = NULL;
    }

    return se;
}

static int sys_add_fcast_variance (equation_system *sys, gretl_matrix *F,
				   int n_static)
{
    gretl_matrix *se = NULL;
    double ftj, vti;
    int n = sys->neqns + sys->nidents;
    int k = F->rows - n_static;
    int i, j, s;
    int err = 0;

    if (k > 0) {
	se = sys_get_fcast_se(sys, k);
	if (se == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    for (j=0; j<n; j++) {
	for (s=0; s<F->rows; s++) {
	    ftj = gretl_matrix_get(F, s, j);
	    if (na(ftj)) {
		gretl_matrix_set(F, s, n + j, NADBL);
	    } else {
		i = s - n_static;
		if (i < 0) {
		    if (j < sys->neqns) {
			vti = sqrt(gretl_matrix_get(sys->Sr, j, j));
		    } else {
			/* LHS of identity */
			vti = 0.0;
		    }
		} else {
		    vti = gretl_matrix_get(se, i, j);
		}
		gretl_matrix_set(F, s, n + j, vti);
	    }
	}
    }

    if (se != NULL) {
	gretl_matrix_free(se);
    }

 bailout:

    if (err) {
	for (i=0; i<n; i++) {
	    for (s=0; s<F->rows; s++) {
		gretl_matrix_set(F, s, n + i, NADBL);
	    }
	}
    }

    return err;
}

static int sys_add_forecast (equation_system *sys,
			     int t1, int t2,
			     const double **Z, const DATAINFO *pdinfo,
			     gretlopt opt)
{
    gretl_matrix *G = NULL;
    gretl_matrix *yl = NULL, *x = NULL;
    gretl_matrix *y = NULL, *yh = NULL;
    const int *ylist;
    const int *ilist;
    const int *xlist;
    int n = sys->neqns + sys->nidents;
    double xit, xitd;
    int tdyn, type, col, T, ncols;
    int i, vi, s, t, lag;
    int err = 0;

    if (sys->Gamma == NULL) {
	return E_DATA;
    }

    ylist = sys->ylist;
    ilist = sys->ilist;
    xlist = sys->xlist;

#if SYSDEBUG
    printlist(ylist, "ylist");
    printlist(ilist, "ilist");
    printlist(xlist, "xlist");
    fprintf(stderr, "sys->order = %d\n", sys->order);
#endif

    if (!gretl_is_identity_matrix(sys->Gamma)) {
	G = gretl_matrix_copy(sys->Gamma);
	if (G == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_SVD_invert_matrix(G);
	}
	if (err) {
	    gretl_matrix_free(G);
	    return err;
	}
    }

    T = t2 - t1 + 1;
    ncols = 2 * n;

    if (opt & OPT_S) {
	tdyn = t2 + 1;
    } else if (opt & OPT_D) {
	tdyn = t1;
    } else {
	tdyn = sys->t2 + 1;
    }

    y = gretl_matrix_alloc(n, 1);
    yh = gretl_matrix_alloc(n, 1);
    sys->F = gretl_zero_matrix_new(T, ncols);

    if (y == NULL || yh == NULL || sys->F == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (sys->order > 0) {
	yl = gretl_matrix_alloc(sys->order * ylist[0], 1);
	if (yl == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    if (xlist[0] > 0) {
	x = gretl_matrix_alloc(xlist[0], 1);
	if (x == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    for (t=t1, s=0; t<=t2; t++, s++) {
	int miss = 0;

	/* lags of endogenous vars */
	if (sys->order > 0) {
	    gretl_matrix_zero(yl);
	    for (i=1; i<=ilist[0] && !miss; i++) {
		vi = ilist[i];
		type = categorize_variable(vi, sys, xlist, &col, &lag);
		if (type == PREDET) {
		    xitd = NADBL;
		    if (t < tdyn || s - lag < 0) {
			/* pre-forecast value */
			xit = Z[vi][t];
		    } else {
			/* prior forecast value preferred */
			if (s - lag >= 0) {
			    xitd = xit = Z[vi][t];
			} 
			xit = gretl_matrix_get(sys->F, s - lag, col);
		    }
		    col += n * (lag - 1);
		    if (!na(xit)) {
			yl->val[col] = xit;
		    } else if (!na(xitd)) {
			yl->val[col] = xitd;
		    } else {
			miss = 1;
		    } 
		}
	    }
	    if (!miss) {
		gretl_matrix_multiply(sys->A, yl, y);
	    }
	} else {
	    gretl_matrix_zero(y);
	}

#if SYSDEBUG > 1
	gretl_matrix_print(yl, "yl");
#endif

	/* exogenous vars */
	if (xlist[0] > 0 && !miss) {
	    for (i=1; i<=xlist[0] && !miss; i++) {
		xit = Z[xlist[i]][t];
		if (na(xit)) {
		    miss = 1;
		} else {
		    x->val[i-1] = xit;
		}
	    }
	    if (!miss) {
		gretl_matrix_multiply_mod(sys->B, GRETL_MOD_NONE,
					  x, GRETL_MOD_NONE,
					  y, GRETL_MOD_CUMULATE);
	    }
	}

	if (miss) {
	    for (i=0; i<n; i++) {
		gretl_matrix_set(sys->F, s, i, NADBL);
	    }
	} else {
	    /* multiply by Gamma^{-1} */
	    if (G != NULL) {
		gretl_matrix_multiply(G, y, yh);
	    } else {
		gretl_matrix_multiply(sys->Gamma, y, yh);
	    }
	    for (i=0; i<n; i++) {
		gretl_matrix_set(sys->F, s, i, yh->val[i]);
	    }
	}
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(yh);
    gretl_matrix_free(yl);
    gretl_matrix_free(x);
    gretl_matrix_free(G);

    if (err) {
	gretl_matrix_free(sys->F);
	sys->F = NULL;
    } else {
	gretl_matrix_set_t1(sys->F, t1);
	gretl_matrix_set_t2(sys->F, t2);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->F, "sys->F, forecasts only");
#endif

    if (!err) {
	sys_add_fcast_variance(sys, sys->F, tdyn - t1);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->F, "sys->F, with variance");
#endif

    return err;
}

const gretl_matrix *
system_get_forecast_matrix (equation_system *sys, int t1, int t2,
			    const double **Z, DATAINFO *pdinfo, 
			    gretlopt opt, int *err)
{
    if (sys->F != NULL) {
	gretl_matrix_free(sys->F);
	sys->F = NULL;
    }
	
    *err = sys_add_forecast(sys, t1, t2, Z, pdinfo, opt);

    return sys->F;
}

/* respond to "save=..." in system command: add residuals and/or
   fitted values to dataset as series
*/

static int sys_save_uhat_yhat (equation_system *sys, double ***pZ,
			       DATAINFO *pdinfo)
{
    int v = pdinfo->v;
    int addvars = 0;
    int err;

    if (sys->flags & SYSTEM_SAVE_UHAT) {
	addvars += sys->neqns;
    } 

    if (sys->flags & SYSTEM_SAVE_YHAT) {
	addvars += sys->neqns;
    }

    if (addvars == 0) {
	return 0;
    }

    err = dataset_add_series(addvars, pZ, pdinfo);

    if (!err) {
	const MODEL *pmod;
	int i, t;

	for (i=0; i<sys->neqns; i++) {
	    pmod = sys->models[i];

	    if (sys->flags & SYSTEM_SAVE_UHAT) {
		for (t=0; t<pdinfo->n; t++) {
		    if (t < pmod->t1 || t > pmod->t2) {
			(*pZ)[v][t] = NADBL;
		    } else {
			(*pZ)[v][t] = pmod->uhat[t];
		    }
		}
		add_system_var_info(sys, i + 1, pdinfo, v++, SYSTEM_SAVE_UHAT);
	    }

	    if (sys->flags & SYSTEM_SAVE_YHAT) {
		for (t=0; t<pdinfo->n; t++) {
		    if (t < pmod->t1 || t > pmod->t2) {
			(*pZ)[v][t] = NADBL;
		    } else {
			(*pZ)[v][t] = pmod->yhat[t];
		    }
		}
		add_system_var_info(sys, i + 1, pdinfo, v++, SYSTEM_SAVE_YHAT);
	    }
	}
    }
    
    return err;
}

int gretl_system_print (const equation_system *sys, const DATAINFO *pdinfo, 
			gretlopt opt, PRN *prn)
{
    int nr = system_n_restrictions(sys);
    int i;

    if (sys->models == NULL) {
	return E_DATA;
    }

    pputc(prn, '\n');

    if (sys->name != NULL) {
	pprintf(prn, "%s, %s\n", _("Equation system"), sys->name);
	pprintf(prn, "%s: %s\n", _("Estimator"), 
		system_get_full_string(sys));
    } else {
	pprintf(prn, "%s, %s\n", _("Equation system"),
		system_get_full_string(sys));
    }

    if (sys->iters > 0) {
	pprintf(prn, _("Convergence achieved after %d iterations\n"), sys->iters);
	if (sys->method == SYS_METHOD_SUR || 
	    sys->method == SYS_METHOD_FIML) {
	    pprintf(prn, "%s = %g\n", _("Log-likelihood"), sys->ll);
	}
    }

    pputc(prn, '\n');

    for (i=0; i<sys->neqns; i++) {
	printmodel(sys->models[i], pdinfo, OPT_NONE, prn);
    }

    system_print_VCV(sys, prn);

    if (nr == 0 && (sys->method == SYS_METHOD_FIML || 
		    sys->method == SYS_METHOD_3SLS || 
		    sys->method == SYS_METHOD_SUR)) {
	print_system_overid_test(sys, prn);
    }

    return 0;
}

int 
system_save_and_print_results (equation_system *sys,
			       double ***pZ, DATAINFO *pdinfo,
			       gretlopt opt, PRN *prn)
{
    int nr = system_n_restrictions(sys);
    int err = 0;

    if (sys->E != NULL) {
	sys->E->t1 = sys->t1;
	sys->E->t2 = sys->t2;
    }

    if (sys->iters > 0 && sys->method == SYS_METHOD_SUR && nr == 0) {
	sur_ols_diag(sys);
    }

    sys->ldet = gretl_vcv_log_determinant(sys->S);

    err = system_add_yhat_matrix(sys);

    if (!err) {
	err = sys_add_structural_form(sys, pdinfo);
    } 

    if (!(opt & OPT_Q)) {
	gretl_system_print(sys, pdinfo, opt, prn);
    }

    if (!err) {
	if (sys->flags & (SYSTEM_SAVE_UHAT | SYSTEM_SAVE_YHAT)) {
	    err = sys_save_uhat_yhat(sys, pZ, pdinfo);
	}
    }

    return err;
}

