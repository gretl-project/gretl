/*
 *  Copyright (c) by Allin Cottrell
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

/* mechanisms for defining and handling systems of equations */

#include "libgretl.h"
#include "system.h"
#include "system_private.h"

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

static void destroy_ident (identity *pident);
static int make_instrument_list (gretl_equation_system *sys);


/* ------------------------------------------------------------ */

static gretl_equation_system **system_stack;
static int n_systems;

gretl_equation_system *
get_equation_system_by_name (const char *sysname, int *snum)
{
    int i;

    for (i=0; i<n_systems; i++) {
	if (!strcmp(sysname, system_stack[i]->name)) {
	    if (snum != NULL) {
		*snum = i;
	    }
	    return system_stack[i];
	}
    }

    return NULL;
}

static int stack_system (gretl_equation_system *sys, PRN *prn)
{
    gretl_equation_system *orig;
    int snum;

    /* only feasible for named systems */
    if (sys == NULL || sys->name == NULL) {
	return 1;
    }

    orig = get_equation_system_by_name(sys->name, &snum);

    if (orig != NULL) {
	/* replace existing system of same name */
	gretl_equation_system_destroy(orig);
	system_stack[snum] = sys;
	pprintf(prn, "Replaced equation system '%s'\n", sys->name);
    } else {
	gretl_equation_system **sstack;

	sstack = realloc(system_stack, (n_systems + 1) * sizeof *sstack);
	if (sstack == NULL) {
	    return E_ALLOC;
	}
	system_stack = sstack;
	system_stack[n_systems++] = sys;
	pprintf(prn, "Added equation system '%s'\n", sys->name);
    }

    return 0;
}

void gretl_equation_systems_cleanup (void)
{
    int i;

    for (i=0; i<n_systems; i++) {
	gretl_equation_system_destroy(system_stack[i]);
    }

    free(system_stack);
    system_stack = NULL;
    n_systems = 0;
}

/* ------------------------------------------------------------ */

static void 
print_system_identity (const identity *pident, const DATAINFO *pdinfo, 
		       PRN *prn)
{
    int i;

    pprintf(prn, "Identity: %s = %s ", 
	    pdinfo->varname[pident->depvar],
	    pdinfo->varname[pident->atoms[0].varnum]);

    for (i=1; i<pident->n_atoms; i++) {
	pprintf(prn, "%c %s ", (pident->atoms[i].op == OP_PLUS)? '+' : '-',
		pdinfo->varname[pident->atoms[i].varnum]);
    }

    pputc(prn, '\n');
}

void 
print_equation_system_info (const gretl_equation_system *sys, 
			    const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    if (sys->name != NULL) {
	pprintf(prn, "Equation system %s\n", sys->name);
    }

    for (i=0; i<sys->n_identities; i++) {
	print_system_identity(sys->idents[i], pdinfo, prn);
    }

    if (sys->endog_vars != NULL) {
	pputs(prn, "Endogenous variables:");
	for (i=1; i<=sys->endog_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->endog_vars[i]]);
	}
	pputc(prn, '\n');
    }

    if (sys->instr_vars != NULL) {
	pputs(prn, "Exogenous variables:");
	for (i=1; i<=sys->instr_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->instr_vars[i]]);
	}
	pputc(prn, '\n');
    }

}

static int gretl_system_method_from_string (const char *str)
{
    int i = 0;

    while (gretl_system_method_strings[i] != NULL) {
	if (!strcmp(str, gretl_system_method_strings[i])) {
	    return i;
	}
	i++;
    }

    return i;
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

    sys->method = method;

    sys->t1 = sys->t2 = 0;

    sys->n_equations = 0;
    sys->n_identities = 0;

    sys->R = NULL;
    sys->q = NULL;

    sys->n_obs = 0;
    sys->iters = 0;
    sys->flags = 0;

    sys->ll = sys->llu = 0.0;
    sys->X2 = 0.0;
    sys->ess = 0.0;
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

    sys->ll = 0.0;
    sys->llu = 0.0;
    sys->X2 = 0.0;
    sys->ess = 0.0;
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

static void gretl_equation_system_clear (gretl_equation_system *sys)
{
    if (sys == NULL || sys->lists == NULL) return;

    sys->flags = 0;
    sys->method = -1;

    /* if restrictions in place, reset the restricted flag */
    if (sys->R != NULL && sys->q != NULL) {
	sys->flags |= GRETL_SYS_RESTRICT;
    }

    system_clear_results(sys);
}

void gretl_equation_system_destroy (gretl_equation_system *sys)
{
    int i;

    if (sys == NULL || sys->lists == NULL) return;

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

static void sur_rearrange_lists (gretl_equation_system *sys)
{
    if (sys->method == SYS_SUR) {
	int i;

	for (i=0; i<sys->n_equations; i++) {
	    rearrange_list(sys->lists[i]);
	}
    }
}

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

    sys->lists[neq] = malloc((list[0] + 1) * sizeof *list);
    if (sys->lists[neq] == NULL) {
	for (i=0; i<neq; i++) {
	    free(sys->lists[i]);
	}
	free(sys->lists);
	sys->lists = NULL;
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
    const char *p;
    int pchars = 0;

    while (isspace((unsigned char) *s)) s++;

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

static int get_estimation_method (const char *s)
{
    char mstr[9];
    int method = -1;

    while (isspace((unsigned char) *s)) s++;

    if (sscanf(s, "%8s", mstr) == 1) {
	lower(mstr);
	method = gretl_system_method_from_string(mstr);
    }

    return method;
}

static char *system_start_get_name (const char *s)
{
    char *sysname = NULL;
    const char *p = strstr(s, "system name=");

    if (p != NULL) {
	sysname = get_system_name_from_line(p + 12);
    }

    return sysname;
}

static int system_start_get_method (const char *s)
{
    int method = -1;
    int offset = 14;
    const char *p;

    p = strstr(s, "system method=");

    if (p == NULL) {
	/* backward compatibility */
	p = strstr(s, "system type=");
	offset = 12;
    }

    if (p != NULL) {
	method = get_estimation_method(p + offset);
    }

    return method;
}

static int named_system_get_method (const char *s)
{
    int method = -1;
    const char *p = strstr(s, "method=");

    if (p != NULL) {
	method = get_estimation_method(p + 7);
    }

    return method;
}

/* Start compiling an equation system in response to gretl's "system"
   command: the command must specify either a "method" (estimation
   method) or a name for the system.  In the former case (method
   given, but no name), the system will be estimated as soon as its
   definition is complete, then it will be destroyed.  In the latter
   case the system definition is saved on a stack, and it can be
   estimated via various methods (the "estimate" command).
*/

gretl_equation_system *system_start (const char *line)
{
    gretl_equation_system *sys = NULL;
    char *sysname = NULL;
    int method;

    method = system_start_get_method(line);

    if (method == SYS_MAX) {
	/* invalid method was given */
	strcpy(gretl_errmsg, _(badsystem));
	return NULL;
    }

    if (method < 0) {
	/* no method was specified: look for a name */
	sysname = system_start_get_name(line);
	if (sysname == NULL) {
	    strcpy(gretl_errmsg, _(badsystem));
	    return NULL;
	}
    }

    sys = gretl_equation_system_new(method, sysname);
    if (sys == NULL) {
	return NULL;
    }

    if (strstr(line, "save=")) {
	if (strstr(line, "resids") || strstr(line, "uhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_UHAT;
	}
	if (strstr(line, "fitted") || strstr(line, "yhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_YHAT;
	}
    }

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
    int vcols = gretl_matrix_cols(vcv);
    int dfu = system_get_dfu(sys);
    int dfn = gretl_matrix_rows(R);

    gretl_matrix *Rbq = NULL;
    gretl_matrix *Rv = NULL;
    gretl_matrix *RvR = NULL;

    double F = NADBL;
    int err = 0;

    if (R == NULL || q == NULL || b == NULL || vcv == NULL) {
	pputs(prn, "Missing data in F test\n");
	goto bailout;
    }

    Rbq = gretl_matrix_alloc(Rrows, 1);
    Rv = gretl_matrix_alloc(Rrows, vcols);
    RvR = gretl_matrix_alloc(Rrows, Rrows);

    if (Rbq == NULL || Rv == NULL || RvR == NULL) {
	pputs(prn, "Out of memory in F test\n");
	goto bailout;
    }

    gretl_matrix_multiply(R, b, Rbq);
    gretl_matrix_subtract_from(Rbq, q);

    gretl_matrix_multiply(R, vcv, Rv);
    gretl_matrix_multiply_mod(Rv, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      RvR);

    err = gretl_invert_symmetric_matrix(RvR);
    if (err) {
	pputs(prn, "Matrix inversion failed in F test\n");
	goto bailout;
    }

    F = gretl_scalar_b_prime_X_b(Rbq, RvR, &err);
    if (err) {
	pputs(prn, "Matrix multiplication failed in F test\n");
	goto bailout;
    }

    F /= dfn;

    pprintf(prn, "%s:\n", _("F test for the specified restrictions"));
    pprintf(prn, "  F(%d,%d) = %g %s %g\n", dfn, dfu, F,
	    _("with p-value"), fdist(F, dfn, dfu));
    pputc(prn, '\n');    

 bailout:
    
    gretl_matrix_free(Rbq);
    gretl_matrix_free(Rv);
    gretl_matrix_free(RvR);
}

static void 
system_print_LR_test (const gretl_equation_system *sys,
		      double llu, PRN *prn)
{
    double llr = system_get_ll(sys);
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
	    df, X2, _("with p-value"), chisq(X2, df));
    pputc(prn, '\n');    
}

static int sys_test_type (gretl_equation_system *sys, gretlopt opt)
{
    int ret = SYS_NO_TEST;

    if (system_n_restrictions(sys) > 0) {
	if (sys->method == SYS_SUR || sys->method == SYS_WLS) {
	    if (opt & OPT_T) {
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

/* driver function for the routines in the "sysest" plugin */

static int 
gretl_equation_system_estimate (gretl_equation_system *sys, 
				double ***pZ, DATAINFO *pdinfo, 
				gretlopt opt, PRN *prn)
{
    int err = 0;
    void *handle = NULL;
    int (*system_est) (gretl_equation_system *, 
		       double ***, DATAINFO *, gretlopt, PRN *);
    int stest = 0;

    *gretl_errmsg = 0;

    err = make_instrument_list(sys);
    if (err) goto system_bailout;
    
    if (sys->method == SYS_SUR) {
	sur_rearrange_lists(sys);
    }

    system_est = get_plugin_function("system_estimate", &handle);

    if (system_est == NULL) {
	err = 1;
	goto system_bailout;
    }

    stest = sys_test_type(sys, opt);

    if (stest == SYS_TEST_NOTIMP) {
	pprintf(prn, _("Sorry, command not available for this estimator"));
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

    if (sys->name == NULL) {
	/* discard the system */
	gretl_equation_system_destroy(sys);
    } else {
	/* retain the system for possible re-estimation */
	gretl_equation_system_clear(sys);
    }

    return err;
}

/* Finalize an equation system in response to "end system".  If the
   system has no name but has a method specified, we go ahead and
   estimate it; otherwise we save it on a stack of defined systems.
*/

int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    gretlopt opt = OPT_NONE;
    *gretl_errmsg = 0;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    if (sys->n_equations < 2) {
	strcpy(gretl_errmsg, _(toofew));
	gretl_equation_system_destroy(sys);
	return 1;
    }

    if (sys->name != NULL) {
	/* save the system for subsequent estimation */
	return stack_system(sys, prn);
    }

    if (sys->method < 0 || sys->method >= SYS_MAX) {
	strcpy(gretl_errmsg, _(badsystem));
	gretl_equation_system_destroy(sys);
	return 1;
    }

    return gretl_equation_system_estimate(sys, pZ, pdinfo, opt, prn);
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

    sysname = get_system_name_from_line(line + 9);
    if (sysname == NULL) {
	strcpy(gretl_errmsg, "estimate: no system name was provided");
	return 1;
    }

    sys = get_equation_system_by_name(sysname, NULL);
    if (sys == NULL) {
	sprintf(gretl_errmsg, "'%s': unrecognized name", sysname);
	free(sysname);
	return 1;
    }

    free(sysname);

    method = named_system_get_method(line);
    if (method < 0 || method >= SYS_MAX) {
	strcpy(gretl_errmsg, "estimate: no valid method was specified");
	return 1;
    }

    sys->method = method;

    if (opt & OPT_T) {
	/* the iterate option is available for WLS, SUR or 3SLS */
	if (method != SYS_WLS && method != SYS_SUR && method != SYS_3SLS) {
	    opt ^= OPT_T;
	}
    }

    /* by default, apply a df correction for single-equation methods */
    if (method == SYS_OLS || method == SYS_WLS ||
	method == SYS_TSLS || method == SYS_LIML) {
	if (!(opt & OPT_N)) {
	    sys->flags |= GRETL_SYSTEM_DFCORR;
	}
    } 

    if (opt & OPT_M) {
	sys->flags |= GRETL_SYS_VCV_GEOMEAN;
    }    

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

const char *gretl_system_get_name (const gretl_equation_system *sys)
{
    return sys->name;
}

int system_adjust_t1t2 (gretl_equation_system *sys,
			int *t1, int *t2, const double **Z)
{
    int i, misst, err = 0;

    for (i=0; i<sys->n_equations && !err; i++) {
	err = adjust_t1t2(NULL, sys->lists[i], t1, t2, Z, &misst);
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

void system_set_n_obs (gretl_equation_system *sys, int n)
{
    sys->n_obs = n;
}

void system_set_iters (gretl_equation_system *sys, int n)
{
    sys->iters = n;
}

/* simple accessor functions */

const char *system_get_full_string (const gretl_equation_system *sys,
				    gretlopt opt)
{
    if (opt & OPT_T) {
	static char sysstr[64];

	sprintf(sysstr, _("iterated %s"), gretl_system_long_strings[sys->method]);
	return sysstr;
    } else {
	return gretl_system_long_strings[sys->method];
    }
}

const gretl_matrix *system_get_R_matrix (const gretl_equation_system *sys)
{
    if (sys->flags & GRETL_SYS_RESTRICT) {
	return sys->R;
    } else {
	return NULL;
    }
}

const gretl_matrix *system_get_q_matrix (const gretl_equation_system *sys)
{
    return sys->q;
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

int system_n_equations (const gretl_equation_system *sys)
{
    return sys->n_equations;
}

int system_n_indentities (const gretl_equation_system *sys)
{
    return sys->n_identities;
}

int system_n_restrictions (const gretl_equation_system *sys)
{
    int nr = 0;

    if (sys->R != NULL && (sys->flags & GRETL_SYS_RESTRICT)) {
	nr = gretl_matrix_rows(sys->R);
    }

    return nr;
}

int system_n_obs (const gretl_equation_system *sys)
{
    return sys->n_obs;
}

int system_iters (const gretl_equation_system *sys)
{
    return sys->iters;
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

gretl_matrix *system_get_uhat (const gretl_equation_system *sys)
{
    return sys->uhat;
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

gretl_matrix *system_get_sigma (const gretl_equation_system *sys)
{
    return sys->sigma;
}



MODEL *system_get_model (const gretl_equation_system *sys, int i)
{
    if (sys->models == NULL || i >= sys->n_equations) {
	return NULL; 
    }   

    return sys->models[i];
}

double system_get_ll (const gretl_equation_system *sys)
{
    return sys->ll;
}

double system_get_llu (const gretl_equation_system *sys)
{
    return sys->llu;
}

double system_get_X2 (const gretl_equation_system *sys)
{
    return sys->X2;
}

double system_get_diag_stat (const gretl_equation_system *sys)
{
    return sys->diag;
}

void system_set_ll (gretl_equation_system *sys, double ll)
{
    sys->ll = ll;
}

void system_set_llu (gretl_equation_system *sys, double llu)
{
    sys->llu = llu;
}

void system_set_X2 (gretl_equation_system *sys, double X2)
{
    sys->X2 = X2;
}

void system_set_ess (gretl_equation_system *sys, double ess)
{
    sys->ess = ess;
}

void system_set_diag_stat (gretl_equation_system *sys, double s)
{
    sys->diag = s;
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

static void destroy_ident (identity *pident)
{
    free(pident->atoms);
    free(pident);
}

static identity *ident_new (int nv)
{
    identity *pident;

    pident = malloc(sizeof *pident);
    if (pident == NULL) return NULL;

    pident->n_atoms = nv;
    pident->atoms = malloc(nv * sizeof *pident->atoms);
    if (pident->atoms == NULL) {
	free(pident);
	pident = NULL;
    }

    return pident;
}

static identity *
parse_identity (const char *str, const DATAINFO *pdinfo, int *err)
{
    identity *pident;
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

    pident = ident_new(nv);
    if (pident == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    pident->depvar = varindex(pdinfo, vname1);
    if (pident->depvar == pdinfo->v) {
	destroy_ident(pident);
	*err = E_UNKVAR;
	return NULL;
    }

    pident->atoms[0].op = OP_PLUS;
    pident->atoms[0].varnum = varindex(pdinfo, vname2);
    if (pident->atoms[0].varnum == pdinfo->v) {
	destroy_ident(pident);
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
	    pident->atoms[i].op = op;
	    pident->atoms[i].varnum = varindex(pdinfo, vname1);
	    if (pident->atoms[i].varnum == pdinfo->v) {
		*err = E_UNKVAR;
	    }
	}
	p++;
    }

    if (*err) {
	destroy_ident(pident);
	pident = NULL;
    }
       
    return pident;
}

static int 
add_identity_to_sys (gretl_equation_system *sys, const char *line,
		     const DATAINFO *pdinfo)
{
    identity **ppident;
    identity *pident;
    int ni = sys->n_identities;
    int err = 0;

    pident = parse_identity(line, pdinfo, &err);
    if (pident == NULL) return err;

    /* connect the identity to the equation system */
    ppident = realloc(sys->idents, (ni + 1) * sizeof *sys->idents);
    if (ppident == NULL) {
	destroy_ident(pident);
	return E_ALLOC;
    }

    sys->idents = ppident;
    sys->idents[ni] = pident;
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
    int i, v, nf, len, cplen;
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
	} else if (sys->method != SYS_3SLS && sys->method != SYS_TSLS) {
	    strcpy(gretl_errmsg, "Instruments may be specified only "
		   "for 3SLS or TSLS");
	    return 1;
	}
    } else {
	return 1;
    }

    nf = count_fields(line);
    if (nf < 1) return 1;

    list = malloc((nf + 1) * sizeof *list);
    if (list == NULL) return E_ALLOC;

    list[0] = nf;
    
    p = line;
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
	    sprintf(gretl_errmsg, "Undefined variable '%s'.", vname);
	    err = 1;
	} else {
	    list[i] = v;
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

int 
system_parse_line (gretl_equation_system *sys, const char *line,
		   const DATAINFO *pdinfo)
{
    int err = 1;

    *gretl_errmsg = '\0';

    if (strncmp(line, "identity", 8) == 0) {
	err = add_identity_to_sys(sys, line + 8, pdinfo);
    } else if (strncmp(line, "endog", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, ENDOG_LIST);
    } else if (strncmp(line, "instr", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, INSTR_LIST);
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



