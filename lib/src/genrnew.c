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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* generate.c for gretl */

/* #define GENR_DEBUG */

#include "libgretl.h"
#include "internal.h"

#include <errno.h>

#define ATOMLEN 32  /* length of auxiliary string in genr atom */

typedef struct _GENERATE GENERATE;
typedef struct _genatom genatom;

struct _GENERATE {
    int err;
    double *xvec;
    int varnum;
    char varname[VNAMELEN];
    char label[MAXLABEL];
    int scalar; 
    int tmpv;
    double **tmpZ;
    DATAINFO *pdinfo;
    double ***pZ;
    MODEL *pmod;
};

struct _genatom {
    int level;
    int scalar;
    int varnum;
    int tmpvar;
    int lag;
    double val;
    int func;
    char op;
    char str[ATOMLEN];
};

static double calc_xy (double x, double y, char op, int t);

static void genrfree (GENERATE *genr);
static void get_lag (int v, int lag, double *lagvec, double **Z, 
		     const DATAINFO *pdinfo);
static double genr_cov (const char *str, double ***pZ,
			const DATAINFO *pdinfo);
static double genr_corr (const char *str, double ***pZ,
			 const DATAINFO *pdinfo);
static double genr_vcv (const char *str, const DATAINFO *pdinfo, 
			MODEL *pmod);
static int genr_mpow (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
static int obs_num (const char *s, const DATAINFO *pdinfo);

#ifdef HAVE_MPFR
static int genr_mlog (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
#endif

static void genr_msg (GENERATE *genr, int oldv);
static void varerror (const char *s);
static int ismatch (int lv, const int *list);
static int genrtime (double ***pZ, DATAINFO *pdinfo, GENERATE *genr, int time);
static int add_new_var (double ***pZ, DATAINFO *pdinfo, GENERATE *genr);

static int math_tokenize (char *s, GENERATE *genr, int level);
static double get_obs_value (const char *s, double **Z, 
			     const DATAINFO *pdinfo);
static int model_variable_index (const char *s);
static int dataset_var_index (const char *s);

static double *get_model_series (double **Z, const DATAINFO *pdinfo,
				 const MODEL *pmod, int v);
static double *get_random_series (DATAINFO *pdinfo, int fn);
static double *get_mp_series (const char *s, GENERATE *genr,
			      int fn, int *err);
static double *get_tmp_series (double **Z, const DATAINFO *pdinfo, 
			       int fn, int varnum, int *err);

static double get_tnum (const DATAINFO *pdinfo, int t);
static double evaluate_statistic (double *z, GENERATE *genr, int fn);
static double get_model_data_element (const char *s, GENERATE *genr,
				      MODEL *pmod, int idx);
static double get_model_statistic (const MODEL *pmod, int idx, int *err);
static double get_dataset_statistic (DATAINFO *pdinfo, int idx);
static double evaluate_math_function (double arg, int fn, int *err);
static double evaluate_bivariate_statistic (const char *s, 
					    GENERATE *genr, 
					    int idx, int *err);
static double evaluate_pvalue (const char *s, double **Z,
			       const DATAINFO *pdinfo, int *err);


enum transformations {
    T_LOG = 1, 
    T_EXP, 
    T_SIN, 
    T_COS,
    T_TAN,
    T_ATAN,
    T_DIFF,
    T_LDIFF, 
    T_MEAN, 
    T_SD, 
    T_MIN,
    T_MAX,
    T_SORT, 
    T_INT, 
    T_LN, 
    T_COEFF,
    T_ABS, 
    T_RHO, 
    T_SQRT, 
    T_SUM, 
    T_NOBS,
    T_NORMAL, 
    T_UNIFORM, 
    T_STDERR,
    T_CUM, 
    T_MISSING,
    T_MISSZERO,
    T_CORR,
    T_VCV,
    T_VAR,
    T_COV,
    T_MEDIAN,
    T_ZEROMISS,
    T_PVALUE,
    T_MPOW,
#ifdef HAVE_MPFR
    T_MLOG,
#endif
    T_SST
};

enum retrieve {
    R_ESS = 1,
    R_NOBS,
    R_T,
    R_RSQ,
    R_SIGMA,
    R_DF,
    R_LNL,
    R_TRSQ,
    R_PD
};

enum composites {
    NEQ = 21,
    GEQ,
    LEQ
};

enum {
    STACK_PUSH,
    STACK_POP,
    STACK_GET,
    STACK_RESET,
    STACK_DESTROY
};

static char *math[] = {
    "log", 
    "exp", 
    "sin", 
    "cos", 
    "tan",
    "atan",
    "diff",
    "ldiff", 
    "mean", 
    "sd", 
    "min",
    "max",
    "sort", 
    "int", 
    "ln", 
    "coeff",
    "abs", 
    "rho", 
    "sqrt", 
    "sum", 
    "nobs",
    "normal", 
    "uniform", 
    "stderr",
    "cum", 
    "missing",
    "misszero",
    "corr",
    "vcv",
    "var",
    "cov",
    "median",
    "zeromiss",
    "pvalue",
    "mpow",
#ifdef HAVE_MPFR
    "mlog",
#endif
    "sst",
    NULL
};

#define LEVELS 7

#define STANDARD_MATH(f) (f == T_LOG || f == T_LN || f == T_EXP || \
                          f == T_SIN || f == T_COS || f == T_TAN || \
                          f == T_ATAN || f == T_INT || f == T_ABS || \
                          f == T_SQRT)

#define SCALAR_SCOPE(t) (t == T_MEAN || t == T_SD || t == T_SUM || \
                         t == T_CORR || t == T_COV || \
                         t == T_VAR || t == T_MEDIAN || t == T_MIN || \
                         t == T_SST || t == T_MAX || t == T_NOBS)

#define MODEL_DATA_ELEMENT(f) (f == T_COEFF || f == T_STDERR || \
                               f == T_RHO || f == T_VCV)

#ifdef HAVE_MPFR
# define MP_MATH(f) (f == T_MPOW || f == T_MLOG)
#else
# define MP_MATH(f) (f == T_MPOW)
#endif

#define MAXTERMS  64
#define TOKLEN   128

/* ...................................................... */

static genatom *make_atom (int scalar, int varnum,
			   int lag, double val,
			   int func, char op, char *str,
			   int level)
{
    genatom *atom;

    atom = malloc(sizeof *atom);
    if (atom == NULL) return NULL;

    atom->scalar = scalar;
    atom->varnum = varnum;
    atom->tmpvar = -1;
    atom->lag = lag;
    atom->val = val;
    atom->func = func;
    atom->op = op;
    atom->level = level;

    *atom->str = '\0';
    strncat(atom->str, str, ATOMLEN - 1);

    return atom;
}

static int print_atom (genatom *atom)
{
    if (atom == NULL) return 1;

    printf(" atom->level = %d\n", atom->level);

    printf(" atom->scalar = %d\n", atom->scalar);

    if (atom->varnum >= 0) {
	printf(" atom->varnum = %d\n", atom->varnum);
    }
    if (atom->tmpvar >= 0) {
	printf(" atom->tmpvar = %d\n", atom->tmpvar);
    }
    if (atom->lag > 0) {
	printf(" atom->lag = %d\n", atom->lag);
    }
    if (atom->scalar) {
	printf(" atom->val = %g\n", atom->val);
    }
    if (atom->func > 0) {
	printf(" atom->func = %d\n", atom->func);
    }

    if (*atom->str) {
	printf(" atom->str = '%s'\n", atom->str);
    }    

    if (isprint(atom->op)) {
	printf(" atom->op = '%c'\n", atom->op);
    } else {
	printf(" atom->op = %d\n", atom->op);
    }

    return 0;
}

static int get_lagvar (const char *s, int *lag, GENERATE *genr)
{
    char vname[VNAMELEN];
    char format[16];
    int m, v = 0;

    sprintf(format, "%%%d[^(](%%d)", VNAMELEN - 1);

    if (sscanf(s, format, vname, &m) == 2) {
	v = varindex(genr->pdinfo, vname);
	if (v == genr->pdinfo->v) {
	    v = m = 0;
	} else if (genr->pdinfo->vector[v] == 0) {
	    sprintf(gretl_errmsg, _("Variable %s is a scalar; "
				    "can't do lags/leads"), 
		    genr->pdinfo->varname[v]);
	    genr->err = E_DATA;
	    v = m = 0;
	}
    }

    if (lag != NULL) *lag = -m;

    return v;
}

static int get_arg_string (char *str, const char *s)
{
    const char *p = strchr(s, '(');
    int n = strlen(p) - 2;
    int err = 0;

    *str = '\0';

    if (n < 0 || n > ATOMLEN - 1) err = 1;
    else strncat(str, p + 1, n);

    return err;
}

static genatom *parse_token (const char *s, char op,
			     GENERATE *genr, int level)
{
    int v = -1;
    int scalar = 0, lag = 0, func = 0;
    double val = 0.0;
    char str[ATOMLEN] = {0};

    if (isalpha(*s)) {
	if (!strchr(s, '(') && !strchr(s, '[')) {
	    if (!strcmp(s, "pi")) {
		val = M_PI;
		scalar = 1;
	    }
	    v = varindex(genr->pdinfo, s);
	    if (v == genr->pdinfo->v) { 
		fprintf(stderr, "unrecognized varname '%s'\n", s);
		genr->err = E_UNKVAR;
	    } else {
		fprintf(stderr, "recognized var '%s' (#%d)\n", s, v);
	    }

	    if (v == INDEXNUM) { /* internal index variable */
		int k = genr_scalar_index(0, 0);

		val = k;
		scalar = 1;
	    }
	    /* handle scalar vars here */
	    if (!genr->pdinfo->vector[v]) {
		val = (*genr->pZ)[v][0];
		scalar = 1;
	    }
	}
	else if (strchr(s, '(')) {
	    /* try for lag variable first */
	    v = get_lagvar(s, &lag, genr);
	    if (v != 0) {
		fprintf(stderr, "recognized var #%d lag %d\n", v, lag);
	    } else {
		/* try for a function */
		v = -1;
		func = math_word_index(s);
		if (func) {
		    fprintf(stderr, "recognized function #%d\n", func);
		    if (MP_MATH(func) || func == T_PVALUE ||
			MODEL_DATA_ELEMENT(func)) {
			genr->err = get_arg_string(str, s);
		    } 
		    if (func == T_PVALUE) {
			if (!genr->err) {
			    val = evaluate_pvalue(str, *genr->pZ,
						  genr->pdinfo,
						  &genr->err);
			    scalar = 1;
			}
		    }
		    else if (MODEL_DATA_ELEMENT(func)) {
			val = get_model_data_element(str, genr,
						     genr->pmod, 
						     func);
			scalar = 1;
		    }
		    if (SCALAR_SCOPE(func)) { 
			fprintf(stderr, "function returns scalar\n");
			/* do stuff */
		    }
		} else {
		    /* not a variable of function: dead end? */
		    genr->err = E_SYNTAX; 
		}
	    } 
	}
	else if (strchr(s, '[')) {
	    val = get_obs_value(s, *genr->pZ, genr->pdinfo);
	    if (val == NADBL) {
		genr->err = E_SYNTAX; 
	    } else {
		scalar = 1;
	    }
	}
	else if (*s == '"') {
	    /* observation label: basis for dummy series */
	    val = obs_num(s, genr->pdinfo);
	}
    }

    else if (*s == '$') {
	int i = model_variable_index(s);

	if (i > 0) {
	    fprintf(stderr, "recognized '%s' as model variable, index #%d\n", 
		    s, i);
	    val = get_model_statistic(genr->pmod, i, &genr->err);
	    scalar = 1;
	} else {
	    i = dataset_var_index(s);
	    if (i > 0) {
		fprintf(stderr, "recognized '%s' as dataset var, index #%d\n", 
			s, i);
		val = get_dataset_statistic(genr->pdinfo, i);
		scalar = 1;
	    } else {
		genr->err = E_UNKVAR;
	    }
	}
    }

    else if (_isnumber(s)) {
	val = dot_atof(s);
	scalar = 1;
    } else {
	genr->err = E_SYNTAX;
    }

    if (genr->err) return NULL;

    return make_atom(scalar, v, lag, val, func, op, str, level);
}

#define STACKSIZE 32

static double calc_stack (double val, int op)
{
    static double valstack[STACKSIZE]; /* ? */
    static int nvals;
    int i;
    double x = 0.0;

    if (op == STACK_PUSH) {
	if (nvals == STACKSIZE - 1) {
	    fprintf(stderr, "stack depth exceeded\n");
	    return x;
	} else {
	    for (i=nvals; i>0; i--) {
		valstack[i] = valstack[i-1];
	    }
	    valstack[0] = val;
	    nvals++;
	}
    }
    else if (op == STACK_POP) {
	x = valstack[0];
	for (i=0; i<nvals; i++) {
	    valstack[i] = valstack[i+1];
	}
	if (nvals > 0) nvals--;
    }
    else if (op == STACK_RESET) {
	for (i=0; i<STACKSIZE; i++) {
	    valstack[i] = 0.0;
	}
	nvals = 0;
    }

    return x;
}

static void calc_push (double x)
{
    calc_stack(x, STACK_PUSH);
}

static double calc_pop (void)
{
    return calc_stack(0., STACK_POP);
}

static void reset_calc_stack (void)
{
    calc_stack(0., STACK_RESET);
}

static genatom *atom_stack (genatom *atom, int op, int i)
{
    static genatom **atoms;
    static int n_atoms;
    static int n_popped;
    genatom *ret = NULL;
    int j;

    if (op == STACK_PUSH && atom != NULL) {
	atoms = realloc(atoms, (n_atoms + 1) * sizeof *atoms);
	if (atoms != NULL) {
	    atoms[n_atoms++] = atom;
	    ret = atom;
	}
    } else if (op == STACK_POP) {
	if (n_popped < n_atoms) {
	    ret = atoms[n_popped++];
	}
    } else if (op == STACK_RESET) {
	n_popped = 0;
    } else if (op == STACK_DESTROY) {
	for (j=0; j<n_atoms; j++) free(atoms[j]);
	free(atoms);
	atoms = NULL;
	n_atoms = 0;
	n_popped = 0;
    }

    return ret;
}

static char *paren_term_stack (char *term, int op, int i)
{
    static char **terms;
    static int n_terms;
    char *ret = NULL;
    int j;

    if (op == STACK_PUSH) {
	terms = realloc(terms, (n_terms + 1) * sizeof *terms);
	if (terms != NULL) {
	    terms[n_terms++] = term;
	    ret = term;
	}
    } else if (op == STACK_GET) {
	if (i >= 0 && i < n_terms) {
	    ret = terms[i];
	}
    } else if (op == STACK_DESTROY) {
	for (j=0; j<n_terms; j++) free(terms[j]);
	free(terms);
	terms = NULL;
	n_terms = 0;
    }

    return ret;
}

static int push_paren_term (char *term)
{
    return (paren_term_stack(term, STACK_PUSH, 0) == NULL);
}

static char *get_paren_term (int i)
{
    return paren_term_stack(NULL, STACK_GET, i);
}

static void destroy_paren_term_stack (void)
{
    paren_term_stack(NULL, STACK_DESTROY, 0);
}

static int push_atom (genatom *atom)
{
    return (atom_stack(atom, STACK_PUSH, 0) == NULL);
}

static genatom *pop_atom (void)
{
    return atom_stack(NULL, STACK_POP, 0);
}

static void reset_atom_stack (void)
{
    atom_stack(NULL, STACK_RESET, 0);
}

static void destroy_atom_stack (void)
{
    atom_stack(NULL, STACK_DESTROY, 0);
}

static double eval_atom (genatom *atom, GENERATE *genr, int t, 
			 double a)
{
    double x = NADBL;

    /* case of constants, scalar variables and specific 
       observation from a series */
    if (atom->scalar) {
	x = atom->val;
	fprintf(stderr, "eval_atom: got scalar = %g\n", x);
    } 

    /* case of temporary variables */
    if (atom->tmpvar >= 0) {
	x = genr->tmpZ[atom->tmpvar][t];
	fprintf(stderr, "eval_atom: got temp obs = %g\n", x);
    }

    /* trend/index variables */
    else if (atom->varnum == TNUM) {
	x = get_tnum(genr->pdinfo, t);
    }

    /* regular variables and lagged terms */
    else if (atom->varnum >= 0) {
	if (!atom->lag) {
	    x = (*genr->pZ)[atom->varnum][t];
	    fprintf(stderr, "eval_atom: got data obs (var %d) = %g\n", 
		    atom->varnum, x);
	} else {
	    int tl = t - atom->lag;

	    if (tl >= 0 && tl < genr->pdinfo->n) {
		x = (*genr->pZ)[atom->varnum][tl];
	    }
	    fprintf(stderr, "eval_atom: got lagged data obs (var %d, "
		    "lag %d)= %g\n", atom->varnum, atom->lag, x);
	}
    }

    /* functions */
    else if (atom->func) {
	if (STANDARD_MATH(atom->func)) {
	    x = evaluate_math_function(a, atom->func, &genr->err);
	    fprintf(stderr, "evaluated func %d arg %g, "
		    "giving %g\n", atom->func, a, x);
	}
    }

    return x;
}

static int genr_add_var (GENERATE *genr, double *x)
{
    double **gZ;    
    int v = genr->tmpv; 
  
    gZ = realloc(genr->tmpZ, (v + 1) * sizeof *gZ);
    if (gZ == NULL) {
	genr->err = E_ALLOC;
	return 1;
    }

    gZ[v] = x;
    genr->tmpZ = gZ;
    genr->tmpv += 1;

    return 0;
}

static int add_random_series_to_genr (GENERATE *genr, genatom *atom)
{
    double *x;

    x = get_random_series(genr->pdinfo, atom->func);
    if (x == NULL) return 1;

    if (genr_add_var(genr, x)) {
	free(x);
	return 1;
    }

    atom->tmpvar = genr->tmpv - 1;

    return 0;
}

static int add_model_series_to_genr (GENERATE *genr, genatom *atom)
{
    double *x;

    x = get_model_series(*genr->pZ, (const DATAINFO *) genr->pdinfo, 
			 genr->pmod, atom->varnum);
    if (x == NULL) return 1;

    if (genr_add_var(genr, x)) {
	free(x);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;

    return genr->err;
}

static int add_mp_series_to_genr (GENERATE *genr, genatom *atom)
{
    double *x;

    x = get_mp_series(atom->str, genr, atom->func, &genr->err);
    if (x == NULL) return 1;

    if (genr_add_var(genr, x)) {
	free(x);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;

    return genr->err;
}

static int add_tmp_series_to_genr (GENERATE *genr, genatom *atom)
{
    double *x;

    /* diff, sort, cum : FIXME compund args */
    x = get_tmp_series(*genr->pZ, (const DATAINFO *) genr->pdinfo, 
		       atom->func, atom->varnum, &genr->err);
    if (x == NULL) return 1;

    if (genr_add_var(genr, x)) {
	free(x);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;

    return genr->err;
}

/* ...................................................... */

static int evaluate_genr (GENERATE *genr)
{
    int t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;
    genatom *atom;
    int err = 0;

    reset_atom_stack();
    while (!err && (atom = pop_atom())) {
	if (atom->varnum == UHATNUM || atom->varnum == YHATNUM) {
	    err = add_model_series_to_genr(genr, atom);
	}
	else if (atom->func == T_UNIFORM || atom->func == T_NORMAL) {
	    err = add_random_series_to_genr(genr, atom);
	}
	else if (atom->func == T_DIFF || atom->func == T_LDIFF ||
		 atom->func == T_CUM || atom->func == T_SORT) {
	    /* FIXME: compound args */
	    err = add_tmp_series_to_genr(genr, atom);
	}
	else if (MP_MATH(atom->func)) {
	    err = add_mp_series_to_genr(genr, atom);
	}
    }

    if (err) {
	genr->err = err;
	return err;
    }

    genr->scalar = 1;
    reset_atom_stack();
    while ((atom = pop_atom())) {
	if (!atom->scalar) genr->scalar = 0;
    }

    if (genr->scalar) t2 = t1;

    for (t=t1; t<=t2; t++) {
	double xbak = 0.0, x = 0.0;
	int level = 0;

	reset_atom_stack();
	while ((atom = pop_atom())) {
	    double y = eval_atom(atom, genr, t, x);

	    if (err) break;
	    if (y == NADBL) {
		x = NADBL;
	    } else {
		if (atom->level < level) {
		    x = calc_pop();
		    fprintf(stderr, "popped %g\n", x);
		}
		if (t == t1) { 
		    fprintf(stderr, "calc_xy: in: x=%g, y=%g, ", x, y);
		    if (isprint(atom->op)) fprintf(stderr, "op='%c'\n", atom->op);
		    else fprintf(stderr, "op=%d\n", atom->op);
		}
		x = calc_xy(x, y, atom->op, t);
		if (t == t1) fprintf(stderr, "calc_xy: out: x=%g\n", x); 
	    }
	    if (atom->level > level) {
		calc_push(xbak);
		fprintf(stderr, "pushed %g\n", xbak);
	    }
	    level = atom->level;
	    xbak = x;
	}
	if (err) break;
	genr->xvec[t] = x;
    }

    genr->err = err;
    return err;
}

/* ...................................................... */

static void catch_double_symbols (char *s)
{
    while (*s) {
	if (*s == '!' && *(s+1) == '=') {
	    *s = NEQ;
	    *(++s) = ' ';
	}
	else if (*s == '>' && *(s+1) == '=') {
	    *s = GEQ;
	    *(++s) = ' ';
	}
	else if (*s == '<' && *(s+1) == '=') {
	    *s = LEQ;
	    *(++s) = ' ';
	}
	else if (*s == '*' && *(s+1) == '*') {
	    *s = '^';
	    *(++s) = ' ';
	}
	s++;
    }
}

/* ...................................................... */

static int op_level (int c)
{
    if (c == '^' || c == '!') 
	return 1;
    if (c == '*' || c == '/' || c == '%') 
	return 2;
    if (c == '+' || c == '-') 
	return 3;
    if (c == '>' || c == '<' || c == GEQ || c == LEQ) 
	return 4;
    if (c == '=' || c == NEQ) 
	return 5;
    if (c == '&') 
	return 6;
    if (c == '|') 
	return 7;
    return 0;
}

static int token_is_atomic (const char *s, GENERATE *genr)
{
    int count = 0;
    const char *p = s;

    /* treat lag variable as atom */
    if (get_lagvar(s, NULL, genr)) return 1;

    while (*p++) {
	/* token is non-atomic if it contains an operator,
	   or parentheses */
	if (op_level(*p) || *p == '(') count++;
	if (count > 0) break;
    }

    return (count == 0);
}

static int string_arg_function (const char *s)
{
    if (!strncmp(s, "coeff", 5) ||
	!strncmp(s, "stderr", 6) ||
	!strncmp(s, "rho", 3) ||
	!strncmp(s, "vcv", 3) ||
	!strncmp(s, "pvalue", 6) ||
	!strncmp(s, "mpow", 4) ||
	!strncmp(s, "mlog", 4)) {
	return 1;
    }

    return 0;
}

static int token_is_function (char *s, GENERATE *genr, int level)
{
    int wlen = 0;
    int ret;
    const char *p = s;

    while (*p) {
	if (isalpha(*p)) wlen++;
	else break;
	p++;
    }

    ret = (wlen > 0 && *p == '(' && p[strlen(p)-1] == ')'); 

    if (ret) {
	fprintf(stderr, "token is function...\n");
	if (string_arg_function(s)) {
	    return ret;
	} else {
	    char subtok[TOKLEN];

	    strcpy(subtok, strchr(s, '(') + 1);
	    subtok[strlen(subtok)-1] = '\0';
	    strcpy(strchr(s, '(') + 1, "#)");
	    if (*subtok) {
		math_tokenize(subtok, genr, ++level);
	    }
	}
    }

    return ret;
}	

static int stack_op_and_token (const char *s, GENERATE *genr, int level)
{
    char tok[TOKLEN];
    int wrapped = 0;
    char op = 0;
    size_t n;

    fprintf(stderr, "stack_op_and_token: looking at '%s'\n", s);

    if (op_level(*s)) {
	/* leading character is an operator */
	op = *s;
	s++;
    } 

    if (isprint(op)) fprintf(stderr, "op = '%c'\n", op);
    else fprintf(stderr, "op = %d\n", op);

    if (*s == '#') {
	char *term = get_paren_term(atoi(s + 1));

	if (term == NULL) genr->err = 1;
	else {
	    if (op) sprintf(tok, "%c", op);
	    else *tok = '\0';
	    strcat(tok, term);
	}
    } else {
	if (*s == '(') {
	    s++;
	    wrapped = 1;
	}
	strcpy(tok, s);
    }

    if (!genr->err) {
	n = strlen(tok);
	if (n > 0 && wrapped && tok[n-1] == ')') tok[n-1] = '\0';
	
	if (!token_is_atomic(tok, genr) &&
	    !token_is_function(tok, genr, level)) {
	    fprintf(stderr, "token is not atomic...\n");
	    genr->err = math_tokenize(tok, genr, ++level);
	} else {
	    genatom *atom;

	    if (*tok == '\0') fprintf(stderr, "got a blank token\n");
	    fprintf(stderr, "token = '%s'\n", tok);
	    atom = parse_token(tok, op, genr, level);
	    if (atom == NULL) {
		genr->err = 1;
		fprintf(stderr, "stack_op_and_token: parse_token() failed\n");
	    } else {
		genr->err = push_atom(atom);
		if (!genr->err) fprintf(stderr, "pushed atom OK\n");
	    }
	}
    }

    return genr->err;
}

static int strip_wrapper_parens (char *s)
{
    int n = strlen(s);
    int strip = 0;

    if (*s == '(' && s[n-1] == ')') {
	int i, pcount = 1;

	strip = 1;
	for (i=1; i<n-1; i++) {
	    if (s[i] == '(') pcount++;
	    else if (s[i] == ')') pcount--;
	    if (pcount == 0) { 
		/* opening paren matched before end */
		strip = 0;
		break;
	    }
	}
	if (strip) s[n-1] = '\0';
    }    

    return strip;
}

static int math_tokenize (char *s, GENERATE *genr, int level)
{
    char tok[TOKLEN];
    const char *q, *p;
    int inparen = 0;

    if (strip_wrapper_parens(s)) s++;

    q = p = s;

    while (*p) {
	if (*p == '(') inparen++;
	else if (*p == ')') inparen--;
	if (inparen < 0) {
	    fprintf(stderr, "error: inparen < 0: '%s'\n", s);
	    return E_UNBAL;
	}

	if (!inparen && op_level(*p)) {
	    if (p - q > TOKLEN - 1) {
		fprintf(stderr, "error: token too long: '%s'\n", q);
		return 1;
	    }
	    if (p - q > 0) {
		*tok = '\0';
		strncat(tok, q, p - q);
		stack_op_and_token(tok, genr, level);
		q = p;
	    }
	}
	p++;
    }

    if (*q) {
	if (strlen(q) > TOKLEN - 1) {
	    fprintf(stderr, "error: remainder too long: '%s'\n", q);
	    return 1;
	}
	strcpy(tok, q);
	genr->err = stack_op_and_token(tok, genr, level);
    }

    return genr->err;
}

/* ...................................................... */

static int count_ops (char *s, int *opcount)
{
    int level, maxlevel = 0;

    while (*s++) {
	level = op_level(*s);
	opcount[level] += 1;
	if (level > maxlevel) maxlevel = level;
    }

    return maxlevel;
}

/* ...................................................... */

static char *cutstr (const char *s, int l, int r)
{
    char *cut = malloc(r - l);

    if (cut == NULL) return NULL;

    *cut = '\0';
    strncat(cut, s + l + 1, r - l - 1);
    return cut;
}

static int pterm_sub (char *orig, int termnum, size_t repl)
{
    char *tmp, termstr[8];

    tmp = malloc(strlen(orig) - repl + 1);
    if (tmp == NULL) return 1;

    strcpy(tmp, orig + repl);
    sprintf(termstr, "#%d", termnum);
    strcpy(orig, termstr);
    strcat(orig, tmp);
    free(tmp);

    return 0;
}

/* ...................................................... */

int insert_paren (char *s, int pos, char lr)
{
    static int lpos;
    static int termnum;
    int i, rpos, n = strlen(s);
    char *cut = NULL;

    if (n + 1 >= MAXLEN) return 1;
    /* move material right */
    for (i=n+1; i>=pos+1; i--) s[i] = s[i - 1];

    if (lr == 'L') {
	lpos = pos + 1;
	s[lpos] = '(';
    } else {
	rpos = pos + 1;
	s[rpos] = ')';
	cut = cutstr(s, lpos, rpos);
	if (cut == NULL) return 1;
	fprintf(stderr, "Made atom: '%s'\n", cut);
	push_paren_term(cut);
	pterm_sub(s + lpos, termnum++, rpos - lpos + 1);
    }

    return 0;
}

/* ...................................................... */

static int paren_state (char c, int *state, char lr)
{
    int s = *state;

    if (c == '(') {
	if (lr == 'L') {
	    if (s > 0) s--;
	} else s++;
    }
    else if (c == ')') {
	if (lr == 'R') {
	    if (s > 0) s--;
	} else s++;
    }

    *state = s;
    return s;
}

/* ...................................................... */

static int parenthesize (char *str)
{
    int i, k, oppos, n;
    int level1 = 0, level2;  
    int priority, start, lpins, inparens;
    int rpar, pbak, maxlevel;
    int opcount[LEVELS + 1];

    for (i=0; i<=LEVELS; i++) opcount[i] = 0;
    maxlevel = count_ops(str, opcount);

    priority = 1;
    k = 0;
    oppos = 0;
    while (priority <= maxlevel) {
	n = strlen(str);
	if (opcount[priority] == 0) {
	    priority++;
	    continue;
	}
	start = oppos + 1;
	oppos = 0;
	lpins = 0;
	for (i=start; i<n; i++) {
	    if ((level1 = op_level(str[i])) == priority) {
		oppos = i;
		break;
	    }
	}
	if (oppos == 0) {
	    k = 0; /* added May 2003 */
	    priority++;
	    continue;
	}

	/* work to left of operator... */
	inparens = 0;
	pbak = 0; 
	for (i=oppos; i>=0; i--) {
	    if (str[i] == '(') pbak++; 
	    else if (str[i] == ')') pbak--;
	    paren_state(str[i], &inparens, 'L');
	    if (inparens) continue;
	    level2 = op_level(str[i]);
	    if (level2 > level1) {
		if (!pbak) {
		    if (insert_paren(str, i, 'L')) return 1;
		    n++;
		    lpins = 1;
		    oppos++;
		}
		break;
	    }
	}
	if (lpins == 0) {
	    continue;
	}

	/* ...and to right of operator */
	inparens = 0;
	rpar = 0;
	for (i=oppos; i<n; i++) {
	    paren_state(str[i], &inparens, 'R');
	    if (inparens) continue;
	    level2 = op_level(str[i]);
	    if (str[i] == '(') rpar--;
	    if (str[i] == ')') rpar++;
	    if (level2 > level1 || i == n - 1 || 
		(str[i] == ')' && rpar == 1)) {
		if (insert_paren(str, (i == n - 1)? i: i - 1, 'R'))
		    return 1;
		n++;
		break;
	    }
	}
	k++;
	if (k == opcount[priority]) {
	    k = 0;
	    oppos = 0;
	    priority++;
	}
    }
    return 0;
}

/* ...................................................... */

int vars_identical (const double *x, const double *y, int n)
     /* check whether two vars are identical or not */
{
    int t;

    for (t=0; t<n; t++) {
	if (floatneq(x[t], y[t])) {
	    return 0;
	}
    }
    return 1;
}

/* ........................................................  */

static void otheruse (const char *str1, const char *str2)
{
    sprintf(gretl_errmsg, _("'%s' refers to a %s and may not be used as a "
			    "variable name"), str1, str2); 
}

/* .......................................................... */

/* not static because used in nls.c */

int math_word_index (const char *s)
{
    char word[9];
    const char *p;
    int i;

    *word = 0;

    p = strchr(s, '(');
    if (p != NULL && p - s <= 8) {
	strncat(word, s, p - s);
    } else {
	strncat(word, s, 8);
    }

    for (i=0; math[i] != NULL; i++) {
	if (!strcmp(word, math[i])) {
	    return i + 1;
	}
    }
    return 0;
}

/* .......................................................... */

int _reserved (const char *str)
{
    const char *resword[] = {"uhat", "yhat",
			     "const", "CONST", "pi",
			     "coeff", "stderr", "rho",
			     "mean", "median", "var", "cov", "vcv", "sd",
			     "full", "subdum", 
			     "t", "annual", "qtrs", "months", "hrs", 
			     "i", "obs", 
			     NULL};
    int i = 0;

    while (resword[i] != NULL) {
	if (strcmp(str, resword[i]) == 0) {
	    switch (i) {
	    case 0: 
		otheruse(str, _("residual vector"));
		break;
	    case 1: 
		otheruse(str, _("fitted values vector"));
		break;
	    case 2: case 3: case 4:
		otheruse(str, _("constant"));
		break;
	    case 5:
		otheruse(str, _("regr. coeff."));
		break;
	    case 6:
		otheruse(str, _("standard error"));
		break;
	    case 7:
		otheruse(str, _("autocorr. coeff."));
		break;
	    case 8: case 9: case 10: case 11: case 12: case 13:
		otheruse(str, _("stats function"));
		break;
	    case 14: case 15:
		otheruse(str, _("sampling concept"));
		break;
	    case 16: case 17: case 18: case 19: case 20:
		otheruse(str, _("plotting variable"));
		break;
	    case 21: case 22: 
		otheruse(str, _("internal variable"));
		break;
	    default:
		otheruse(str, _("math function"));
		break;
	    }
	    return 1;
	}
	i++; 
    } 

    if (math_word_index(str)) {
	otheruse(str, _("math function"));
	return 1;
    }
 
    return 0;
}

/* .........................................................    */

static int getword (char *word, char *str, char c, unsigned char oflag)
     /* Scans string str for char c, gets word to the left of it as
	"word" and deletes word from str.
	Returns number of chars deleted, or -1 if no occurrence of c, 
	or 0 if reserved word is used 
     */
{
    int i;

    i = haschar(c, str);
    if (i == -1) return -1;

    *word = '\0';
    strncat(word, str, i);
    _delete(str, 0, ++i);

    if (oflag && !strcmp(word, "subdum")) ;
    else if (_reserved(word)) i = 0;

    return i;
}

/* ........................................................... */

static void fix_decimal_commas (char *str)
{
    char *p = str;

    if (p == NULL || *p == 0) return;
    p++;
    
    while (*p && *(p + 1)) {
	if (*p == ',' && isdigit(*(p - 1)) && isdigit(*(p + 1)))
	    *p = '.';
	p++;
    }
}

/* ........................................................... */

static void get_genr_formula (char *formula, const char *line)
{
    if (line == NULL || *line == 0) return;

    /* skip over " genr " */
    while (isspace((unsigned char) *line)) line++;
    if (!strncmp(line, "genr", 4)) {
	line += 4;
	while (isspace((unsigned char) *line)) line++;
    }
    *formula = '\0';
    strncat(formula, line, MAXLEN - 10);
}

/**
 * genr_scalar_index:
 * @opt: If opt = 1, set the value of the (static) index, using
 * the value of @put.  If opt = 2, increment the static index by
 * the value of @put.
 *
 * Reads the value of a static index variable (after setting or
 * incrementing the index if @opt is non-zero).
 * 
 * Returns: the new value of the index.
 */

int genr_scalar_index (int opt, int put)
{
    /* opt = 1, set index (using "put")
       opt = 2, increment index value
       Refers to an "internal" variable named "i",
       available in genr commands, and with ID number 1001
    */
    static int i;

    if (opt == 1) i = put;
    else if (opt == 2) i += put;
    return i;
}

static int gentoler (const char *s)
{
    int ret = 0;
    double x;

    if (_isnumber(s)) {
	x = dot_atof(s);
	set_nls_toler(x);
	sprintf(gretl_msg, _("Set tolerance to %g"), x);

    } else {
	strcpy(gretl_errmsg, _("The setting for \"toler\" must be numeric"));
	ret = 1;
    }

    return ret;
}

static void make_genr_label (int replmsg, char *genrs, 
			     int model_count,
			     GENERATE *genr)
{
    if (replmsg) {
	sprintf(genr->label, _("Replaced after model %d: "), 
		model_count);
    }	
    if (strlen(genrs) > MAXLABEL - 1) {
	strncat(genr->label, genrs, MAXLABEL - 4);
	strcat(genr->label, "...");
    } else {
	strncat(genr->label, genrs, MAXLABEL - 1);
    }
}

/**
 * generate:
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @line: command line for parsing.
 * @model_count: count of models estimated so far.
 * @pmod: pointer to a model, or NULL.
 * @oflag: option flag (relates to generation of dummy variables).
 *
 * Generates a new variable, usually via some transformation of
 * existing variables, or by retrieving an internal variable associated
 * with the estimation of a model (@pmod).
 * 
 * Returns: 0 on success, integer error code on error.
 */

int generate (double ***pZ, DATAINFO *pdinfo, 
	      const char *line, int model_count, 
	      MODEL *pmod, unsigned char oflag)
{
    int i, v = 0;
    char s[MAXLEN], genrs[MAXLEN];
    char newvar[32];
    int oldv = pdinfo->v;
    GENERATE genr;
    genatom *atom;

    *gretl_errmsg = '\0';

    genr.err = 0;
    genr.scalar = 1;
    *genr.label = '\0';
    genr.pdinfo = pdinfo;
    genr.pZ = pZ;
    genr.pmod = pmod;
    genr.tmpZ = NULL;
    genr.tmpv = 0;

    *s = *genrs = '\0';
    get_genr_formula(s, line);
    delchar('\n', s);
    strcpy(genrs, s);
    catch_double_symbols(s);
    delchar(' ', s);

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint())
	fix_decimal_commas(s);
#endif

#ifdef GENR_DEBUG
    fprintf(stderr, "\n*** starting genr, s='%s'\n", s);
#endif
 
    if (strcmp(s, "dummy") == 0) {
	genr.err = dummy(pZ, pdinfo);
	if (!genr.err)
	    strcpy(gretl_msg, _("Periodic dummy variables generated.\n"));
	return genr.err;
    }
    else if (strcmp(s, "paneldum") == 0) {
	genr.err = paneldum(pZ, pdinfo, oflag);
	if (!genr.err)
	    strcpy(gretl_msg, _("Panel dummy variables generated.\n"));
	return genr.err;
    }
    else if (strcmp(s, "index") == 0) { 
	genr.err = genrtime(pZ, pdinfo, &genr, 0);
	if (!genr.err) genr_msg(&genr, oldv);
	return genr.err;
    }
    else if (strcmp(s, "time") == 0) {
	genr.err = genrtime(pZ, pdinfo, &genr, 1);
	if (!genr.err) genr_msg(&genr, oldv);
	return genr.err;
    }
    else if (strncmp(s, "toler=", 6) == 0) {
	genr.err = gentoler(s + 6);
	return genr.err;
    }

    *newvar = '\0';
    
    /* get equation newvar = s, where s is expression */
    i = getword(newvar, s, '=', oflag);

    if (i > 0) {
	if (*newvar == '\0') {
	    genr.err = E_NOVAR;
	    goto genr_return;
	}
	_esl_trunc(newvar, VNAMELEN - 1);
	if (!isalpha((unsigned char) *newvar) &&
	    strncmp(newvar, "$nls", 4)) {
	    genr.err = E_NOTALPH;
	    goto genr_return;
	}
	v = varindex(pdinfo, newvar);
	if (v == 0) { 
	    genr.err = E_CONST;
	    goto genr_return;
	}
	if (lastchar('=', s)) {
	    genr.err = E_EQN;
	    goto genr_return;
	}
    } else {
	genr.err = E_NOEQ;
	goto genr_return;
    }

    /* basic memory allocation */
    if ((genr.xvec = malloc(pdinfo->n * sizeof *genr.xvec)) == NULL) {
	genr.err = E_ALLOC;
	goto genr_return;
    }

    for (i=0; i<pdinfo->n; i++) {
	genr.xvec[i] = 0.0;
    }

    /* impose operator hierarchy */
    if (parenthesize(s)) { 
	fprintf(stderr, "genr: parenthesize failed\n");
	genr.err = E_ALLOC;
	goto genr_return;
    }
#ifdef GENR_DEBUG
    fprintf(stderr, "after parenthesize: s='%s'\n", s);
#endif

    genr.err = math_tokenize(s, &genr, 0);

    destroy_paren_term_stack();

    i = 0;
    while ((atom = pop_atom())) {
	printf("*** atom %d ***\n", i++);
	print_atom(atom);
    }

    if (!genr.err) {
	evaluate_genr(&genr);
    }

    destroy_atom_stack();
    reset_calc_stack();

 genr_return:
    if (genr.err) {
	genrfree(&genr);
    } else {
	strcpy(genr.varname, newvar);
	genr.varnum = v;
	genr_msg(&genr, oldv);
	make_genr_label(v < oldv && !oflag && model_count > 0,
			genrs, model_count, &genr);
	genr.err = add_new_var(pZ, pdinfo, &genr); 
    }

    return genr.err;
}

/* ........................................................... */
    
static int add_new_var (double ***pZ, DATAINFO *pdinfo, GENERATE *genr)
{
    int t, n = pdinfo->n, v = genr->varnum;
    int modify = 0, old_scalar = 0;
    double xx;

    /* is the new variable an addition to data set? */
    if (v >= pdinfo->v) {
	if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
	strcpy(pdinfo->varname[v], genr->varname);
    } else {
	modify = 1;
	if (!pdinfo->vector[v]) old_scalar = 1;
    }

    strcpy(VARLABEL(pdinfo, v), genr->label);
    pdinfo->vector[v] = !genr->scalar;
    xx = genr->xvec[pdinfo->t1];

    if (genr->scalar) {
	strcat(VARLABEL(pdinfo, v), _(" (scalar)"));
	(*pZ)[v] = realloc((*pZ)[v], sizeof ***pZ);
	(*pZ)[v][0] = xx;
    } else {
	if (old_scalar) {
	    (*pZ)[v] = realloc((*pZ)[v], pdinfo->n * sizeof ***pZ);
	    if ((*pZ)[v] == NULL) return E_ALLOC;
	}
	if (!modify) {
	    for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
	}
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	    (*pZ)[v][t] = genr->xvec[t];
    }

    genrfree(genr);

    return 0;
}

/* ............................................................ */

static double calc_xy (double x, double y, char op, int t) 
{
    long int ny;
    double xx, yy;

    switch (op) {
    case '\0':
	x = y;
	break;
    case '+':
	x += y;
	break;
    case '|':
	if (floatneq(x, 0.) || floatneq(y, 0.))
	    x = 1.0;
	else
	    x = 0.0;
	break;
    case '-':
	x -= y;
	break;
    case '*':
	x *= y;
	break;
    case '&':
	x *= y;
	if (x != 0.) x = 1.0;
	break;
    case '%':
	x = (double) ((int) x % (int) y);
	break;
    case '/':
	xx = y;
	if (floateq(xx, 0.0)) {  
	    sprintf(gretl_errmsg, _("Zero denominator for obs %d"), t+1);
	    return 1;
	}
	x /= xx;
	break;
    case '^':
	xx = x;
	yy = y;
	ny = (long) yy;
	if ((floateq(xx, 0.0) && yy <= 0.0) || 
	    (xx < 0.0 && (double) ny != yy)) {
	    sprintf(gretl_errmsg, 
		    _("Invalid power function args for obs. %d"
		      "\nbase value = %f, exponent = %f"), t, xx, yy);
	    return 1;
	}
	if (floateq(xx, 0.0)) x = 0.0;
	else x = pow(xx, yy);
	break;
    case '<':
	if (x < y) x = 1.0;
	else x = 0.0;
	break;
    case '>':
	if (x > y) x = 1.0;
	else x = 0.0;
	break;
    case '=':
	if (floateq(x, y)) x = 1.0;
	else x = 0.0;
	break;
    case NEQ: /* not equals */
	if (floateq(x, y)) x = 0.0;
	else x = 1.0;
	break;
    case GEQ: /* greater than or equal */
	if (floateq(x, y)) x = 1.0;
	else if (x > y) x = 1.0;
	else x = 0.0;
	break;
    case LEQ: /* less than or equal */
	if (floateq(x, y)) x = 1.0;
	else if (x < y) x = 1.0;
	else x = 0.0;
	break;
    case '!':
	if (floatneq(y, 0.0)) x = 0.0;
	else x = 1.0;
	break;
    } 

    return x;
}

/* ........................................................  */

int panel_unit_first_obs (int t, const DATAINFO *pdinfo)
{
    char *p, obs[OBSLEN];

    ntodate(obs, t, pdinfo);
    p = strchr(obs, ':');
    if (p != NULL && atoi(p + 1) == 1) return 1;
    return 0;
}

/* below: math functions taking scalar arg and returning scalar */

static double evaluate_math_function (double arg, int fn, int *err)
{
    double x = NADBL;

    if (na(arg)) return x;

    *err = 0;

    switch (fn) {

    case T_LOG:
    case T_LN:
	if (arg <= 0.0) {
	    *err = E_LOGS;
	} else {
	    x = log(arg);
	}
	break;	
    case T_EXP:
	if (exp(arg) == HUGE_VAL) {
	    fprintf(stderr, "genr: exponent = %g\n", arg);
	    *err = E_HIGH;
	} else {
	    x = exp(arg);
	}
	break;
    case T_SIN:
	x = sin(arg);
	break;
    case T_COS:
	x = cos(arg);
	break;
    case T_TAN:
	x = tan(arg);
	break;
    case T_ATAN:
	x = atan(arg);
	break;
    case T_INT:
	x = (double) (int) arg;
	break;
    case T_ABS:
	x = (arg < 0.0)? -arg : arg;
	break;
    case T_SQRT:
	if (arg < 0.0) {
	    *err = E_SQRT;
	} else {
	    x = sqrt(arg);
	}
	break;
    default:
	break;

    }

    return x;
}

static double *get_mp_series (const char *s, GENERATE *genr,
			      int fn, int *err)
{
    double *x;

    x = malloc(genr->pdinfo->n * sizeof *x);
    if (x == NULL) return NULL;

    if (fn == T_MPOW) {
	*err = genr_mpow(s, x, *genr->pZ, genr->pdinfo);
	if (*err) *err = E_INVARG;
    }
#ifdef HAVE_MPFR
    else if (fn == T_MLOG) {
	*err = genr_mlog(s, x, *genr->pZ, genr->pdinfo);
	if (*err) *err = E_INVARG;
    }
#endif

    return x;
}

static double *get_random_series (DATAINFO *pdinfo, int fn)
{
    double *x;

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) return NULL;

    if (fn == T_NORMAL) {
	gretl_normal_dist(x, pdinfo->t1, pdinfo->t2);
    }   
    else if (fn == T_UNIFORM) {
	gretl_uniform_dist(x, pdinfo->t1, pdinfo->t2);
    }

    return x;
}

/* below: functions taking series (z) as input and returning a scalar statistic */
/* thinks: should these be rolling calculations? */

static double evaluate_statistic (double *z, GENERATE *genr, int fn)
{
    double x = NADBL;
    double *tmp = NULL;
    int i, t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;

    if (fn == T_NOBS) {
	/* doesn't need a series allocated */
	i = 0;
	for (t=t1; t<=t2; t++) {
	    if (!na(z[t])) i++;
	}
	return (double) i;
    }	

    tmp = malloc((t2 - t1 + 1) * sizeof *tmp);
    if (tmp == NULL) {
	genr->err = E_ALLOC;
	return x;
    }

    i = -1;
    for (t=t1; t<=t2; t++) {
	x = z[t];
	if (na(x)) continue;
	tmp[++i] = x;
    }

    if (fn == T_MEAN) {
	x = _esl_mean(0, i, tmp);
    }
    else if (fn == T_SUM) {
	x = _esl_mean(0, i, tmp);
	x *= (i + 1);
    }
    else if (fn == T_SD) {
	x = _esl_stddev(0, i, tmp);
    }
    else if (fn == T_VAR) {
	x = _esl_variance(0, i, tmp);
    }
    else if (fn == T_SST) {
	x = _esl_sst(0, i, tmp);
    }
    else if (fn == T_MEDIAN) {
	x = gretl_median(tmp, i + 1);
    }
    else if (fn == T_MIN || fn == T_MAX) {
	double min, max;

	_minmax(0, i, tmp, &min, &max);
	x = (fn == T_MIN)? min : max;
    }

    free(tmp);
    
    return x;
}

static double evaluate_pvalue (const char *s, double **Z,
			       const DATAINFO *pdinfo, int *err)
{
    double x;

    x = batch_pvalue(s, Z, pdinfo, NULL);
    if (na(x) || x == -1.0) *err = E_INVARG;
    return x;
}

static double 
evaluate_bivariate_statistic (const char *s, GENERATE *genr, 
			      int idx, int *err)
{
    double x = NADBL;

    if (idx == T_COV) {
	x = genr_cov(s, genr->pZ, genr->pdinfo);
	if (na(x)) *err = E_INVARG;
    }
    else if (idx == T_CORR) {
	x = genr_corr(s, genr->pZ, genr->pdinfo);
	if (na(x)) *err = E_INVARG;
    }

    return x;
}

/* below: functions pertaining to "missing" status of data values */

static double evaluate_missval_func (double arg, int fn)
{
    double x = NADBL;

    if (fn == T_MISSING) {
	/* check whether obs is missing or not */
	x = (na(arg))? 1.0 : 0.0;
    }
    else if (fn == T_MISSZERO) {
	/* change missing obs to zero */
	x = (na(arg))? 0.0 : arg;
    }
    else if (fn == T_ZEROMISS) {
	/* change zero to missing obs */
	x = (floateq(arg, 0.0))? NADBL : arg;
    }

    return x;
}

static double *get_tmp_series (double **Z, const DATAINFO *pdinfo, 
			       int fn, int varnum, int *err)
{
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2; 
    double *x, *mvec;
    double xx, yy;

    x = malloc(pdinfo->n * sizeof *x); 
    if (x == NULL) return NULL;

    mvec = Z[varnum];

    if (fn == T_DIFF) {
	for (t=t1+1; t<=t2; t++) {
	    if (pdinfo->time_series == STACKED_TIME_SERIES &&
		panel_unit_first_obs(t, pdinfo)) {
		x[t] = NADBL;
		continue;
	    }
	    xx = mvec[t];
	    if (pdinfo->time_series == STACKED_CROSS_SECTION) {
		yy = (t - pdinfo->pd >= 0)? mvec[t-pdinfo->pd] : NADBL;
	    } else {
		yy = mvec[t-1];
	    }
	    x[t] = (na(xx) || na(yy))? NADBL : xx - yy;
	}
	x[t1] = NADBL;
    }

    else if (fn == T_LDIFF) {
	for (t=t1+1; t<=t2; t++) {
	    if (pdinfo->time_series == STACKED_TIME_SERIES &&
		panel_unit_first_obs(t, pdinfo)) {
		x[t] = NADBL;
		continue;
	    }
	    xx = mvec[t];
	    if (pdinfo->time_series == STACKED_CROSS_SECTION) {
		yy = (t - pdinfo->pd >= 0)? mvec[t-pdinfo->pd] : NADBL;
	    } else {
		yy = mvec[t-1];
	    }
	    if (na(xx) || na(yy)) {
		x[t] = NADBL;
		continue;
	    }   
	    else if (xx <= 0.0 || yy <= 0.0) *err = E_LOGS;
	    else x[t] = log(xx) - log(yy);
	}
	x[t1] = NADBL;
    }

    else if (fn == T_CUM) {
	x[t1] = (na(mvec[t1])) ? 0.0 : mvec[t1];
	for (t=t1; t<=t2; t++) {
	    if (na(mvec[t])) x[t] = x[t-1];
	    else x[t] = x[t-1] + mvec[t];
	}
    }

    else if (fn == T_SORT) {
	double *tmp = malloc((t2 - t1 + 1) * sizeof *tmp);
	int i;

	if (tmp == NULL) {
	    free(x);
	    return NULL;
	}

	i = -1;
	for (t=t1; t<=t2; t++) {
	    if (na(mvec[t])) continue;
	    tmp[++i] = mvec[t];
	}

	qsort(tmp, i + 1, sizeof *tmp, _compare_doubles);

	for (t=t1; t<=t2; t++) {
	    x[t] = tmp[t-t1];
	}
	free(tmp);
    }

    return x;
}

/* ...........................................................*/

static int check_modelstat (const MODEL *pmod, int idx)
{
    if (pmod == NULL || pmod->list == NULL) {
	switch (idx) {
	case R_T:
	    strcpy(gretl_errmsg, 
		   _("No $T (number of obs for model) value is available"));
	    return 1;
	case R_ESS:
	    strcpy(gretl_errmsg, 
		   _("No $ess (error sum of squares) value is available"));
	    return 1;
	case R_RSQ:
	    strcpy(gretl_errmsg, 
		   _("No $rsq (R-squared) value is available"));
	    return 1;
	case R_TRSQ:
	    strcpy(gretl_errmsg, 
		   _("No $trsq (T*R-squared) value is available"));
	    return 1;
	case R_DF:
	    strcpy(gretl_errmsg, 
		   _("No $df (degrees of freedom) value is available"));
	    return 1;
	case R_SIGMA:
	    strcpy(gretl_errmsg, 
		   _("No $sigma (std. err. of model) value is available"));
	    return 1;
	case R_LNL:
	    strcpy(gretl_errmsg, 
		   _("No $lnl (log-likelihood) value is available"));
	    return 1;
	default:
	    return 0;
	}
    }

    if (pmod != NULL && pmod->ci != LOGIT && pmod->ci != PROBIT &&
	idx == R_LNL) {
	strcpy(gretl_errmsg, 
	       _("$lnl (log-likelihood) is not available for the last model"));
	return 1;
    }

    return 0;
}

static double 
get_model_data_element (const char *s, GENERATE *genr,
			MODEL *pmod, int idx)
{
    int lv, vi;
    double x = NADBL;

    if (pmod == NULL) return x;

    if (idx == T_RHO) {
	if (!(_isnumber(s))) {
	    genr->err = E_INVARG;
	}
	else if (dot_atof(s) == 1 && 
		 (pmod->ci == CORC || pmod->ci == HILU)) {
	    x = pmod->rho_in;
	}
	else if (pmod->ci != AR && dot_atof(s) == 1) {
	    x = pmod->rho;
	}
	else if (pmod->arinfo == NULL || 
		 pmod->arinfo->arlist == NULL || 
		 pmod->arinfo->rho == NULL) {
	    genr->err = E_INVARG;
	}
	else if (!(vi = ismatch(atoi(s), pmod->arinfo->arlist))) {
	    genr->err = E_INVARG;
	} else {
	    x = pmod->arinfo->rho[vi];
	}
    }

    else if (idx == T_VCV) {
	x = genr_vcv(s, genr->pdinfo, pmod);
	if (na(x)) genr->err = E_INVARG;
    }

    else if (idx == T_COEFF || idx == T_STDERR) {
	if (pmod == NULL || pmod->list == NULL) {
	    genr->err = E_INVARG;
	    return NADBL;
	}

	lv = _isnumber(s)? atoi(s) : varindex(genr->pdinfo, s);
	vi = ismatch(lv, pmod->list);

	if (vi == 1) vi = 0;
	if (!vi) {
	    genr->err = E_INVARG;
	    return NADBL;
	}

	if (idx == T_COEFF && pmod->coeff != NULL) { 
	    x = pmod->coeff[vi-2];
	} else if (pmod->sderr != NULL) {
	    x = pmod->sderr[vi-2];
	} else {
	    genr->err = E_INVARG;
	}
    } 

    return x;
}

/* retrieve scalar statistic from model */

static double 
get_model_statistic (const MODEL *pmod, int idx, int *err)
{
    double x = NADBL;

    if (pmod == NULL) {
	*err = 1;
    }
    else if (check_modelstat(pmod, idx)) {
	*err = 1;
    }
    else if ((pmod->ci == LOGIT || pmod->ci == PROBIT) &&
	(idx == R_RSQ || idx == R_ESS || idx == R_SIGMA || 
	 idx == R_TRSQ)) {
	*err = E_BADSTAT;
    }
    else if (pmod->ci == LAD && (idx == R_RSQ || idx == R_TRSQ)) {
	*err = E_BADSTAT;
    }

    if (*err) return x;

    switch (idx) {  
    case R_ESS:
	x = pmod->ess;
	break;
    case R_RSQ:
	x = pmod->rsq;
	break;
    case R_LNL:
	x = pmod->lnL;
	break;
    case R_SIGMA:
	if (pmod->nwt) x = pmod->sigma_wt;
	else x = pmod->sigma;
	break;
    case R_TRSQ:
	x = pmod->nobs * pmod->rsq;
	break;
    case R_DF:
	x = (double) pmod->dfd;
	break;
    case R_T:
	x = (double) pmod->nobs;
	break;	
    }

    return x;
}

/* ...........................................................*/

static double get_dataset_statistic (DATAINFO *pdinfo, int idx)
{
    double x = NADBL;

    if (pdinfo == NULL) return x;

    switch (idx) {
    case R_NOBS:
	x = (double) (pdinfo->t2 - pdinfo->t1 + 1);
	break;
    case R_PD:
	x = (double) pdinfo->pd;
	break;
    default:
	break;
    }

    return x;
}

/* ...........................................................*/

static double get_obsnum_boolean (double val, int t)
{
    if ((int) val == t + 1) return 1.0; /* FIXME? */
    else return 0.0;
}

/* ...........................................................*/

static double get_obs_value (const char *s, double **Z, 
			     const DATAINFO *pdinfo)
{
    char vname[VNAMELEN], obs[OBSLEN];

    if (sscanf(s, "%8[^[][%8[^]]]", vname, obs) != 2) {
	return NADBL;
    } else {
	int i = varindex(pdinfo, vname);
	int t = dateton(obs, pdinfo);

	if (i < pdinfo->v && pdinfo->vector[i] && 
	    t >= 0 && t < pdinfo->n) {
	    return Z[i][t];
	}
    }
    return NADBL;
}

/* ...........................................................*/

static double *get_model_series (double **Z, const DATAINFO *pdinfo,
				 const MODEL *pmod, int v)
{
    int t, t2, n = pdinfo->n;
    double *x;

    if (pmod == NULL || (v != UHATNUM && v != YHATNUM))
	return NULL;

    if (pmod->t2 - pmod->t1 + 1 > n ||
	model_sample_issue(pmod, NULL, Z, pdinfo)) {
	strcpy(gretl_errmsg, (v == UHATNUM)? 
	       _("Can't retrieve uhat: data set has changed") :
	       _("Can't retrieve yhat: data set has changed"));
	return NULL;
    }   

    if ((v == UHATNUM && pmod->uhat == NULL) ||
	(v == YHATNUM && pmod->yhat == NULL)) {
	return NULL;
    }

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) return NULL;

    if (pmod->data != NULL) {
	t2 = pmod->t2 + get_misscount(pmod);
    } else {
	t2 = pmod->t2;
    }    

    for (t=0; t<pmod->t1; t++) x[t] = NADBL;
    for (t=pmod->t1; t<=t2; t++) {
	x[t] = (v == UHATNUM)? pmod->uhat[t] : pmod->yhat[t];
    }
    for (t=t2+1; t<n; t++) x[t] = NADBL;
	    
	    
    return x;
}

/* ...........................................................*/

static double get_tnum (const DATAINFO *pdinfo, int t)
{
    if (pdinfo->time_series && pdinfo->pd == 1) {
	/* annual data: let 't' be the year */ 
	return pdinfo->sd0 + t;
    } else {
	/* let 't' be the 1-based observation number */
	return (double) (t + 1);
    }
}

/* ..................................................................*/

static void get_lag (int v, int lag, double *lagvec, double **Z, 
		     const DATAINFO *pdinfo)
{
    register int t;
    int t1, lt;

    t1 = (lag > pdinfo->t1)? lag : pdinfo->t1;

    for (t=0; t<pdinfo->n; t++) lagvec[t] = NADBL;

    /* stacked X-section needs rather special handling */
    if (pdinfo->time_series == STACKED_CROSS_SECTION) {
	for (t=t1; t<=pdinfo->t2; t++) { 
	    lt = t - lag * pdinfo->pd;
	    if (lt < 0 || lt >= pdinfo->n) continue;
	    lagvec[t] = Z[v][lt];
	}
    }
    else if (dated_daily_data(pdinfo)) {
	for (t=t1; t<=pdinfo->t2; t++) {
	    lt = t - lag;
	    while (lt >= 0 && na(Z[v][lt])) lt--;
	    lagvec[t] = Z[v][lt];
	}
    } 
    else { /* the "standard" time-series case */
	for (t=t1; t<=pdinfo->t2; t++) {
	    lt = t - lag;
	    if (lt < 0 || lt >= pdinfo->n) continue;
	    lagvec[t] = Z[v][lt];
	}
    }

    /* post-process missing panel values */
    if (pdinfo->time_series == STACKED_TIME_SERIES) {
	char *p, obs[OBSLEN];
	int j;

	for (t=t1; t<=pdinfo->t2; t++) {
	    ntodate(obs, t, pdinfo);
	    p = strchr(obs, ':');
	    j = atoi(p + 1);
	    if (j <= lag) lagvec[t] = NADBL;
	}
    }
}


/* ......................................................   */

static int obs_num (const char *s, const DATAINFO *pdinfo)
{
    char test[OBSLEN];
    size_t n;
    int t;

    *test = 0;
    strncat(test, (*s == '"')? s + 1 : s, OBSLEN - 1);
    n = strlen(test);
    if (test[n-1] == '"') test[n-1] = '\0';

    if (pdinfo->markers && pdinfo->S != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    if (!strcmp(test, pdinfo->S[t])) return t + 1;
	}
    }

    if (pdinfo->time_series == TIME_SERIES) {
	t = dateton(test, pdinfo);
	if (t >= 0) return t + 1;
    }

    return 0;
}

/* ......................................................   */

static int dataset_var_index (const char *s)
{
    char test[9];

    *test = '\0';
    strncat(test, s, 8);
    lower(test);

    if (!strcmp(test, "$nobs")) 
	return R_NOBS;
    if (!strcmp(test, "$pd")) 
	return R_PD;

    return 0;
}

static int model_variable_index (const char *s)
{
    char test[9];

    *test = '\0';
    strncat(test, s, 8);
    lower(test);

    if (!strcmp(test, "$ess"))  
	return R_ESS;
    if (!strcmp(test, "$t")) 
	return R_T;
    if (!strcmp(test, "$rsq"))  
	return R_RSQ;
    if (!strcmp(test, "$sigma"))  
	return R_SIGMA;
    if (!strcmp(test, "$df"))   
	return R_DF;
    if (!strcmp(test, "$lnl"))   
	return R_LNL;
    if (!strcmp(test, "$nrsq") || 
	!strcmp(test, "$trsq")) 
	return R_TRSQ;

    return 0;
}

/* ........................................................  */

static void get_month_name (char *mname, int m)
{
    struct tm mt;

    mt.tm_sec = 0;
    mt.tm_min = 0;
    mt.tm_hour = 0;
    mt.tm_mday = 1;
    mt.tm_mon = m - 1;
    mt.tm_year = 100;

    strftime(mname, 7, "%b", &mt);
    *mname = tolower(*mname);
}

/**
 * dummy:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Adds to the data set a set of periodic (usually seasonal)
 * dummy variables.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int dummy (double ***pZ, DATAINFO *pdinfo)
{
    static char word[16];
    int vi, t, yy, pp, mm;
    int nvar = pdinfo->v;
    int ndummies = pdinfo->pd;
    double xx;

    if (ndummies == 1) return E_PDWRONG;
    if (dataset_add_vars(ndummies, pZ, pdinfo)) return E_ALLOC;

    mm = (pdinfo->pd < 10)? 10 : 100;
    for (vi=1; vi<=ndummies; vi++) {
	if (pdinfo->pd == 4 && pdinfo->time_series == TIME_SERIES) {
	    sprintf(word, "dq%d", vi);
	    sprintf(VARLABEL(pdinfo, nvar+vi-1), _("= 1 if quarter = %d, "
						   "0 otherwise"), vi);
	} 
	else if (pdinfo->pd == 12 && pdinfo->time_series == TIME_SERIES) {
	    char mname[8];

	    get_month_name(mname, vi);
	    sprintf(word, "d%s", mname);
	    sprintf(VARLABEL(pdinfo, nvar+vi-1), _("= 1 if month is %s, "
						   "0 otherwise"), mname);
	} else {
	    sprintf(word, "dummy_%d", vi);
	    sprintf(VARLABEL(pdinfo, nvar+vi-1), _("%s = 1 if period is %d, "
						   "0 otherwise"), word, vi);
	}
	strcpy(pdinfo->varname[nvar+vi-1], word);

	for (t=0; t<pdinfo->n; t++) {
	    xx = date(t, pdinfo->pd, pdinfo->sd0);
	    yy = (int) xx;
	    pp = (int) (mm * (xx - yy) + 0.5);
	    (*pZ)[nvar+vi-1][t] = (pp == vi)? 1.0 : 0.0;
	}
    }
    return 0;
}

/**
 * paneldum:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: 0 for stacked time-series, 1 for stacked cross-sections.
 *
 * Adds to the data set a set of panel data dummy variables (for
 * both unit and period).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int paneldum (double ***pZ, DATAINFO *pdinfo, unsigned char opt)
     /* creates panel data dummies (unit and period) 
	opt = 0 for stacked time-series, 
	non-zero for stacked cross-section
     */
{
    static char word[16];
    int vi, t, yy, pp, mm;
    int nvar = pdinfo->v;
    int ndum, nudum, ntdum;
    double xx;

    ntdum = pdinfo->pd;
    if (ntdum == 1) return E_PDWRONG;

    nudum = pdinfo->n / pdinfo->pd;
    if (nudum == 1) return E_PDWRONG;

    ndum = ntdum + nudum;
    if (dataset_add_vars(ndum, pZ, pdinfo)) return E_ALLOC;

    /* first generate the frequency-based dummies */
    mm = (pdinfo->pd < 10)? 10 : 100;
    for (vi=1; vi<=ntdum; vi++) {
	if (opt) sprintf(word, "du_%d", vi);
	else sprintf(word, "dt_%d", vi);
	strcpy(pdinfo->varname[nvar+vi-1], word);
	sprintf(VARLABEL(pdinfo, nvar+vi-1), _("%s = 1 if %s is %d, "
					       "0 otherwise"), word, 
		(opt)? _("unit"): _("period"), vi);
	for (t=0; t<pdinfo->n; t++) {
	    xx = date(t, pdinfo->pd, pdinfo->sd0);
	    yy = (int) xx;
	    pp = (int) (mm*(xx - yy) + 0.5);
	    (*pZ)[nvar+vi-1][t] = (pp == vi)? 1.0 : 0.0;
	}
    }

    /* and then the block-based ones */
    for (vi=1; vi<=nudum; vi++) {
	if (opt) sprintf(word, "dt_%d", vi);
	else sprintf(word, "du_%d", vi);
	strcpy(pdinfo->varname[nvar+ntdum+vi-1], word);
	sprintf(VARLABEL(pdinfo, nvar+ntdum+vi-1), _("%s = 1 if %s is %d, "
						     "0 otherwise"), word, 
		(opt)? _("period"): _("unit"), vi);
	for (t=0; t<pdinfo->n; t++) 
	    (*pZ)[nvar+ntdum+vi-1][t] = 0.0;
	for (t=(vi-1)*pdinfo->pd; t<vi*pdinfo->pd; t++) 
	    (*pZ)[nvar+ntdum+vi-1][t] = 1.0;
    }

    return 0;
}

/* ........................................................  */

static int genrtime (double ***pZ, DATAINFO *pdinfo, GENERATE *genr,
		     int time)
     /* create time trend variable */
{
    int i, t, n = pdinfo->n, v = pdinfo->v;

    if (time) i = varindex(pdinfo, "time");
    else i = varindex(pdinfo, "index");

    if (i == v) {
	if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
    }
    if (time) {
	strcpy(genr->varname, "time");
	strcpy(pdinfo->varname[i], "time");
	strcpy(VARLABEL(pdinfo, i), _("time trend variable"));
    } else {
	strcpy(genr->varname, "index");
	strcpy(pdinfo->varname[i], "index");
	strcpy(VARLABEL(pdinfo, i), _("data index variable"));
    }

    for (t=0; t<n; t++) (*pZ)[i][t] = (double) (t + 1);

    genr->varnum = i;
    genr->scalar = 0;

    return 0;
}

/**
 * plotvar:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @period: string to identify periodicity: "annual", "qtrs",
 * "months", "time" or "index".
 *
 * Adds to the data set a special dummy variable for use in plotting.
 *
 * Returns: the ID number of the variable (> 0) or -1 on failure
 */

int plotvar (double ***pZ, DATAINFO *pdinfo, const char *period)
{
    int t, vi, y1, n = pdinfo->n;
    float rm;

    vi = varindex(pdinfo, period);
    if (vi < pdinfo->v) return vi;

    if (dataset_add_vars(1, pZ, pdinfo)) return -1;

    strcpy(pdinfo->varname[vi], period);

    y1 = (int) pdinfo->sd0;
    rm = pdinfo->sd0 - y1;

    switch(period[0]) {
    case 'a':
	strcpy(VARLABEL(pdinfo, vi), _("annual plotting variable")); 
	for (t=0; t<n; t++) 
	    (*pZ)[vi][t] = (double) (t + atoi(pdinfo->stobs));
	break;
    case 'q':
	strcpy(VARLABEL(pdinfo, vi), _("quarterly plotting variable"));
	(*pZ)[vi][0] = y1 + (10.0 * rm - 1.0)/4.0;
	for (t=1; t<n; t++) 
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + .25;
	break;
    case 'm':
	strcpy(VARLABEL(pdinfo, vi), _("monthly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0)/12.0;
	for (t=1; t<n; t++) 
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0/12.0);
	break;
    case 'h':
	strcpy(VARLABEL(pdinfo, vi), _("hourly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0)/24.0;
	for (t=1; t<n; t++) 
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0/24.0);
	break; 
    case 'i':
	strcpy(VARLABEL(pdinfo, vi), _("index variable"));
	for (t=0; t<n; t++) (*pZ)[vi][t] = (double) (t + 1);
	break;
    case 't':
	strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	for (t=0; t<n; t++) (*pZ)[vi][t] = (double) (t + 1);
	break;
    default:
	break;
    }

    return vi;
}

/* ......................................................  */

/* laggenr: create Z[iv][t-lag] if this variable does not
   already exist.  

   Return the ID number of the lag var, or -1 on error.
*/

int newlag; /* library global */

int laggenr (int parent, int lag, int opt, double ***pZ, 
	     DATAINFO *pdinfo)
{
    char word[32];
    char s[32];
    int lno;
    double *lx;

    /* can't do lags of a scalar */
    if (!pdinfo->vector[parent]) return -1;

    lx = malloc(pdinfo->n * sizeof *lx);
    if (lx == NULL) return -1;

    strcpy(s, pdinfo->varname[parent]);
    if (pdinfo->pd >=10) _esl_trunc(s, 5);
    else _esl_trunc(s, 6);
    sprintf(word, "_%d", lag);
    strcat(s, word);
    lno = varindex(pdinfo, s);

    /* put the lag values into array lx */
    get_lag(parent, lag, lx, *pZ, pdinfo);

    newlag = 1;

    if (lno < pdinfo->v) {
	/* a variable of this name already exists */
	if (vars_identical(lx, (*pZ)[lno], pdinfo->n)) {
	    /* and it is just what we want */
	    free(lx);
	    newlag = 0;
	} else {
	    /* but the values are wrong: swap them */
	    free((*pZ)[lno]);
	    (*pZ)[lno] = lx;
	}
    } else {
	/* no var of this name, working from scratch */
	dataset_add_allocated_var(lx, pZ, pdinfo);
	strcpy(pdinfo->varname[lno], s);
	if (opt) { 
	    sprintf(VARLABEL(pdinfo, lno), "%s = %s(-%d)", s, 
		    pdinfo->varname[parent], lag);
	}
    }

    return lno;
}

/**
 * varlist:
 * @pdinfo: data information struct.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the variables currently defined.
 */

void varlist (const DATAINFO *pdinfo, PRN *prn)
{
    int i, n = pdinfo->v;

    pprintf(prn, _("Listing %d variables:\n"), n);
    for (i=0; i<n; i++) {
	pprintf(prn, "%3d) %-10s", i, pdinfo->varname[i]);
	if ((i+1) % 5 == 0) 
	    pputc(prn, '\n');
    }
    if (n % 5) pputc(prn, '\n');
}

/**
 * varindex:
 * @pdinfo: data information struct.
 * @varname: name of variable to test.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name.
 *
 */

int varindex (const DATAINFO *pdinfo, const char *varname)
{
    int i;

    if (varname == NULL) return pdinfo->v;

    if (!strcmp(varname, "uhat")) return UHATNUM; 
    if (!strcmp(varname, "yhat")) return YHATNUM; 
    if (!strcmp(varname, "i"))    return INDEXNUM;
    if (!strcmp(varname, "t"))    return TNUM;
    if (!strcmp(varname, "obs"))  return TNUM;

    if (!strcmp(varname, "const") || !strcmp(varname, "CONST"))
	return 0;

    for (i=0; i<pdinfo->v; i++) { 
	if (!strcmp(pdinfo->varname[i], varname)) { 
	    return i;
	}
    }

    return pdinfo->v;
}

/* ........................................................ */

static void genrfree (GENERATE *genr)
{
    int i;

    if (genr == NULL) return;

    if (genr->tmpv > 0) {
	for (i=0; i<genr->tmpv; i++) free(genr->tmpZ[i]);
	free(genr->tmpZ);
    }

    free(genr->xvec);
}

/**
 * logs:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set the natural logs of the
 * variables given in @list.
 *
 * Returns: the number of variables generated, or -1 on failure.
 */

int logs (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    register int i;
    int j, t, v, nvar = pdinfo->v, n = pdinfo->n;
    int check, le_zero;
    int l0 = list[0];
    double xx;
    char s[32];

    if (dataset_add_vars(l0, pZ, pdinfo)) return -1;

    j = 0;
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	if (v == 0) continue; /* dont try to take log of constant */
	/* and don't try to take the log of a dummy variable */
	if (isdummy((*pZ)[v], pdinfo->t1, pdinfo->t2))
	    continue;
	if (v < nvar)  { 
	    le_zero = 0;
	    for (t=0; t<n; t++) (*pZ)[nvar+j][t] = NADBL;
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		xx = (pdinfo->vector[v])? (*pZ)[v][t] : (*pZ)[v][0];
		if (xx <= 0.0) {
		    (*pZ)[nvar+j][t] = NADBL;
		    if (!na(xx)) {
			sprintf(gretl_errmsg, 
				_("Log error: Variable '%s', obs %d,"
				  " value = %g\n"), pdinfo->varname[v],
				t+1, xx);
			le_zero = 1;
		    }
		}
		else (*pZ)[nvar+j][t] = log(xx); 
	    }
	    if (le_zero) continue;
	    strcpy(s, "l_");
	    strcat(s, pdinfo->varname[v]);
	    _esl_trunc(s, 8);
	    strcpy(pdinfo->varname[nvar+j], s);
	    strcat(s, _(" = log of "));
	    strcat(s, pdinfo->varname[v]);
	    strcpy(VARLABEL(pdinfo, nvar+j), s);
	    check = varindex(pdinfo, pdinfo->varname[j]);
	    if (check < nvar) {
		if (pdinfo->vector[check]) {
		    if (vars_identical((*pZ)[check], (*pZ)[nvar+j], n)) {
			j--;
		    }
		}
	    } 
	} else varerror(s);
	j++;
    }

    /* shrink Z if warranted (not all vars logged) */
    if (j < l0) dataset_drop_vars(l0 - j, pZ, pdinfo);

    if (j == 0) j = -1;

    return j;
}

/**
 * lags:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set lagged values of the 
 * variables given in @list (up to the frequency of the data).
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int lags (const LIST list, double ***pZ, DATAINFO *pdinfo)
     /* generates lag variables for each var in list */
{
    int check, l, v, lv;
    int maxlag = pdinfo->pd;

    /* play safe with panel data */
    if (dataset_is_panel(pdinfo)) maxlag = 1;
    
    for (v=1; v<=list[0]; v++) {
	lv = list[v];
	if (lv == 0 || !pdinfo->vector[lv]) continue;
	for (l=1; l<=maxlag; l++) {
	    check = laggenr(lv, l, 1, pZ, pdinfo);
	    if (check < 0) return 1;
	}
    }

    return 0;
}

/**
 * xpxgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: If = 0, only squares are generated, if non-zero, both
 * squares and cross-products are generated.
 * @nodup: If non-zero, variables will not be created if they
 * are already present in the data set.
 *
 * Generates and adds to the data set squares and (if @opt is non-zero) 
 * cross-products of the variables given in @list.
 *
 * Returns: The number of variables generated, or -1 on error.
 */

int xpxgenr (const LIST list, double ***pZ, DATAINFO *pdinfo, 
	     int opt, int nodup)
{
    int check, i, j, t, li, lj, l0 = list[0];
    int maxterms, terms, n = pdinfo->n, v = pdinfo->v;
    double zi, zj;
    char s[12], s1[VNAMELEN];

    /* maximum number of terms if none are "bad" */
    if (opt) maxterms = (l0*l0 + l0)/2;
    else maxterms = l0;

    if (dataset_add_vars(maxterms, pZ, pdinfo)) return -1;

    terms = 0;
    for (i=1; i<=l0; i++) {
	li = list[i];
	if (!isdummy((*pZ)[li], 0, n-1)) {
	    for (t=0; t<n; t++) (*pZ)[v+terms][t] = NADBL;
	    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		zi = (*pZ)[li][t];
		if (na(zi)) (*pZ)[v+terms][t] = NADBL;
		else (*pZ)[v+terms][t] = zi * zi;
	    }
	    if (_iszero(0, n-1, (*pZ)[v+terms])) continue; 
	    /*
	      prefix varname by sq, truncate if too long and save under 
	      new varname; new label is "varname = oldname squared"
	    */
	    strcpy(s, "sq_");
	    strcat(s, pdinfo->varname[li]);
	    _esl_trunc(s, 8);
	    strcpy(pdinfo->varname[v+terms], s);
	    /* check if an identical variable exists? */
	    if (nodup) {
		check = varindex(pdinfo, pdinfo->varname[(v+terms)]);
		if (check < v) {
		    if (vars_identical((*pZ)[check], (*pZ)[v+terms], n)) 
			continue;
		}
	    }
	    sprintf(VARLABEL(pdinfo, v+terms), _("%s = %s squared"), s,
		    pdinfo->varname[li]);  
	    terms++;
	}
	/* also do cross-products if wanted */
	if (opt) {
	    for (j=i+1; j<=l0; j++) {
		lj = list[j];
		for (t=0; t<n; t++) (*pZ)[v+terms][t] = NADBL;
		for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		    zi = (*pZ)[li][t];
		    zj = (*pZ)[lj][t];
		    if (na(zi) || na(zj)) 
			(*pZ)[v+terms][t] = NADBL;
		    else (*pZ)[v+terms][t] = zi*zj;
		}
		if (_iszero(0, n-1, (*pZ)[v+terms])) continue;
		/*
		  trunc varname i and varname j if needed and cat them.
		  save as newvarname.  Also make label.
		*/
		strcpy(s, pdinfo->varname[li]);
		_esl_trunc(s, 3);
		strcat(s, "_");
		strcpy(s1, pdinfo->varname[lj]);
		_esl_trunc(s1, 4);
		strcat(s, s1);
		strcpy(pdinfo->varname[v+terms], s);
		sprintf(VARLABEL(pdinfo, v+terms), _("%s = %s times %s"),
			s, pdinfo->varname[li], pdinfo->varname[lj]);
		terms++;
	    }
	}
    }

    if (terms < maxterms) 
	dataset_drop_vars(maxterms - terms, pZ, pdinfo);

    return terms;
}

/**
 * rhodiff:
 * @param: please see the gretl help on rhodiff() for syntax.
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generates and adds to the data set rho-differenced versions
 * of the variables given in @list.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int rhodiff (char *param, const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i, j, maxlag, p, t, t1, nv, v = pdinfo->v, n = pdinfo->n;
    char s[64], parmbit[VNAMELEN];
    double xx, *rhot;

#ifdef GENR_DEBUG
    fprintf(stderr, "rhodiff: param = '%s'\n", param);
#endif
    maxlag = _count_fields(param);
    rhot = malloc(maxlag * sizeof *rhot);
    if (rhot == NULL) return E_ALLOC;
    if (maxlag > pdinfo->t1) t1 = maxlag;
    else t1 = pdinfo->t1;

#ifdef GENR_DEBUG
    fprintf(stderr, "rhodiff: maxlag = %d, t1 = %d\n", maxlag, t1);
#endif

    /* parse "param" string */
    j = strlen(param);
    p = 0;
    for (i=0; i<j; i++) {
	if ((i == 0 || param[i] == ' ') && i < (j - 1)) {
	    sscanf(param + i + (i? 1: 0), "%8s", parmbit); 
#ifdef GENR_DEBUG
	    fprintf(stderr, "rhodiff: parmbit = '%s'\n", parmbit);
#endif
	    if (isalpha((unsigned char) parmbit[0])) {
		nv = varindex(pdinfo, parmbit);
		if (nv == v) {
		    free(rhot);
		    return E_UNKVAR;
		}
		rhot[p] = get_xvalue(nv, *pZ, pdinfo);
	    } else {
		rhot[p] = dot_atof(parmbit);
	    }
	    p++;
	}
    }

    if (dataset_add_vars(list[0], pZ, pdinfo)) return E_ALLOC;

    for (i=1; i<=list[0]; i++) {
	j = list[i];
#ifdef GENR_DEBUG
	fprintf(stderr, "rhodiff: doing list[%d] = %d\n", i, list[i]);
#endif
	/* make name and label */
	strcpy(s, pdinfo->varname[j]);
	_esl_trunc(s, 7);
	strcat(s, "#");
	strcpy(pdinfo->varname[v+i-1], s);
	sprintf(VARLABEL(pdinfo, v+i-1), _("%s = rho-differenced %s"), 
		pdinfo->varname[v+i-1], pdinfo->varname[j]);
	/* fill out values */
	for (t=0; t<n; t++) (*pZ)[v+i-1][t] = NADBL;
	for (t=t1; t<=pdinfo->t2; t++) {
	    xx = (*pZ)[j][t];
	    if (na(xx)) {
		(*pZ)[v+i-1][t] = NADBL;
		continue;
	    }
	    for (p=0; p<maxlag; p++) {
		if (na((*pZ)[j][t-p-1])) {
		    xx = NADBL;
		    break;
		}
		else xx -= rhot[p] * (*pZ)[j][t-p-1];
	    }
	    (*pZ)[v+i-1][t] = xx;
	}
    }

    free(rhot);

    return 0;
}

/* ...................................................... */

static int genr_mpow (const char *str, double *xvec, double **Z, 
		      DATAINFO *pdinfo)
{
    int err, v;
    unsigned pwr;
    char vname[VNAMELEN];
    void *handle = NULL;
    int (*mp_raise) (const double *, double *, int, unsigned);
    
    if (sscanf(str, "%[^,],%u", vname, &pwr) != 2) {
	return 1;
    }

    v = varindex(pdinfo, vname);
    if (v >= pdinfo->v) {
	return 1;
    } 

    if (open_plugin("mp_ols", &handle)) {
	fprintf(stderr, I_("Couldn't access GMP plugin\n"));
	return 1;
    }

    mp_raise = 
	get_plugin_function("mp_vector_raise_to_power", handle);

    if (mp_raise == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	close_plugin(handle);
	return 1;
    }

    err = mp_raise(Z[v], xvec, pdinfo->n, pwr);

    close_plugin(handle);
    
    return err;
}

#ifdef HAVE_MPFR

static int genr_mlog (const char *str, double *xvec, double **Z, 
		      DATAINFO *pdinfo)
{
    int err, v;
    char vname[VNAMELEN];
    void *handle = NULL;
    int (*mp_log) (const double *, double *, int);
    
    if (sscanf(str, "%8s", vname) != 1) {
	return 1;
    }

    v = varindex(pdinfo, vname);
    if (v >= pdinfo->v) {
	return 1;
    } 

    if (open_plugin("mp_ols", &handle)) {
	fprintf(stderr, I_("Couldn't access GMP plugin\n"));
	return 1;
    }

    mp_log = get_plugin_function("mp_vector_ln", handle);

    if (mp_log == NULL) {
	fprintf(stderr, I_("Couldn't load plugin function\n"));
	close_plugin(handle);
	return 1;
    }

    err = mp_log(Z[v], xvec, pdinfo->n);

    close_plugin(handle);
    
    return err;
}

#endif /* HAVE_MPFR */

/* ...................................................... */

static double genr_cov (const char *str, double ***pZ, 
			const DATAINFO *pdinfo)
{
    int i, n, p, v1, v2;
    char v1str[VNAMELEN], v2str[VNAMELEN];

    n = strlen(str);
    if (n > 17) return NADBL;

    p = haschar(',', str);
    if (p < 0 || p > 8) return NADBL;

    /* get first var name */
    for (i=0; i<p; i++) v1str[i] = str[i];
    v1str[p] = '\0';

    /* get second var name */
    n = n - p - 1;
    for (i=0; i<n; i++) v2str[i] = str[p+1+i];
    v2str[i] = '\0';

    /* and look up the two */
    v1 = varindex(pdinfo, v1str);
    v2 = varindex(pdinfo, v2str);
    if (v1 >= pdinfo->v || v2 >= pdinfo->v)
	return NADBL;

    return _covar(pdinfo->t2 - pdinfo->t1 + 1,
		  &(*pZ)[v1][pdinfo->t1], 
		  &(*pZ)[v2][pdinfo->t1]);
}

/* ...................................................... */

static double genr_corr (const char *str, double ***pZ, 
			 const DATAINFO *pdinfo)
{
    int i, n, p, v1, v2;
    char v1str[VNAMELEN], v2str[VNAMELEN];

    n = strlen(str);
    if (n > 17) return NADBL;

    p = haschar(',', str);
    if (p < 0 || p > 8) return NADBL;

    /* get first var name */
    for (i=0; i<p; i++) v1str[i] = str[i];
    v1str[p] = '\0';

    /* get second var name */
    n = n - p - 1;
    for (i=0; i<n; i++) v2str[i] = str[p+1+i];
    v2str[i] = '\0';

    /* and look up the two */
    v1 = varindex(pdinfo, v1str);
    v2 = varindex(pdinfo, v2str);
    if (v1 >= pdinfo->v || v2 >= pdinfo->v)
	return NADBL;

    return _corr(pdinfo->t2 - pdinfo->t1 + 1,
		 &(*pZ)[v1][pdinfo->t1], &(*pZ)[v2][pdinfo->t1]);
}

/* ...................................................... */

static int get_nls_param_number (const MODEL *pmod, 
				 const char *vname)
{
    int i;

    if (pmod->params == NULL) return 0;

    for (i=0; i<=pmod->ncoeff; i++) {
	if (!strcmp(vname, pmod->params[i])) return i + 1;
    }

    return 0;
}

/* ...................................................... */

static double genr_vcv (const char *str, const DATAINFO *pdinfo, 
			MODEL *pmod)
{
    int i, j, k, n, nv, p, v1, v2, v1l, v2l;
    char v1str[VNAMELEN], v2str[VNAMELEN];

    if (pmod == NULL || pmod->list == NULL) return NADBL;

    n = strlen(str);
    if (n > 17) return NADBL;

    p = haschar(',', str);
    if (p < 0 || p > 8) return NADBL;

    /* get first var name */
    for (i=0; i<p; i++) v1str[i] = str[i];
    v1str[p] = '\0';

    /* get second var name */
    n = n - p - 1;
    for (i=0; i<n; i++) v2str[i] = str[p+1+i];
    v2str[i] = '\0';

    /* are they valid? */
    v1 = varindex(pdinfo, v1str);
    v2 = varindex(pdinfo, v2str);
    if (v1 >= pdinfo->v || v2 >= pdinfo->v) return NADBL;

    /* check model list */
    if (pmod->ci == NLS) {
	v1l = get_nls_param_number(pmod, v1str);
	v2l = get_nls_param_number(pmod, v2str);
    } else {
	v1l = ismatch(v1, pmod->list);
	v2l = ismatch(v2, pmod->list);
    }
    if (!v1l || !v2l) return NADBL;

    /* make model vcv matrix if need be */
    if (pmod->vcv == NULL && makevcv(pmod)) return NADBL;

    /* now find the right entry */
    nv = pmod->list[0];
    if (v1l > v2l) {
	k = v1l;
	v1l = v2l;
	v2l = k;
    }
    k = 0;
    for (i=2; i<=nv; i++) {
	for (j=2; j<=nv; j++) {
	    if (j < i) continue;
	    if (i == v1l && j == v2l) return pmod->vcv[k];
	    k++;
	}
    }

    return NADBL;
}

/* ...................................................... */

static void genr_msg (GENERATE *genr, int oldv)
{
    sprintf(gretl_msg, "%s %s %s (ID %d)", 
	    (genr->varnum < oldv)? _("Replaced") : _("Generated"), 
	    (genr->scalar)? _("scalar") : _("vector"),
	    genr->varname, genr->varnum);
    if (genr->scalar) {
	char numstr[24];

	sprintf(numstr, " = %g", genr->xvec[0]);
	strcat(gretl_msg, numstr);
    }
}

/* ......................................................  */

static int ismatch (int lv, const int *list)
{
    int i;

    for (i=list[0]; i>=1; i--) {
	if (lv == list[i]) return i;
    }
    return 0;
}

/* .......................................................... */

static void varerror (const char *s)
     /* print error message for variable not in name list */
{
    sprintf(gretl_errmsg, _("Undefined variable name '%s'"), s);

    if (!strcmp(s, "const")) {
	sprintf(gretl_errmsg, _("const cannot be used to store values"));
    } else if (!strcmp(s, "uhat")) {
	sprintf(gretl_errmsg,
		_("uhat can be used only in genr.  First use the command: "
		  "genr newname = uhat"));
    } else if (*s == '$') {
	sprintf(gretl_errmsg, _("Reserved var. names starting with "
				"$ can be used only in genr.\nFirst use the "
				"command:  genr newname = %s"), s);
    }
}

/* .......................................................... */

/* ensure there's no gap betweeen a minus sign and the term
   it qualifies, in a "sim" command */

static char *close_minus (char *str)
{
    char *q, *p = str;

    while ((p = strchr(p, '-'))) {
	q = ++p;
	while (*q == ' ') q++;
	if (q != p) {
	    memmove(p, q, strlen(q) + 1);
	}
    }

    return str;
}

/* construct a descriptive label for a variable that is
   modified via the "sim" command */

static char *make_sim_label (char *label, const char *vname,
			     char **parm, int nparm)
{
    int k, neg, started = 0;
    char term[32];

    sprintf(label, "%s(t)=", vname);

    for (k=0; k<nparm; k++) {
	if (isdigit(*parm[k]) && dot_atof(parm[k]) == 0.0) {
	    continue;
	}
	if (k == 0) {
	    strcpy(term, parm[k]);
	} else {
	    neg = (*parm[k] == '-');
	    sprintf(term, "%s%s*%s(t-%d)", (neg)? "-" : 
		    (started)? "+" : "",
		    parm[k] + neg, vname, k);
	}
	if (strlen(label) + strlen(term) >= MAXLABEL - 4) {
	    if (strlen(label) < MAXLABEL - 4) {
		strcat(label, "...");
	    }
	    break;
	} else {
	    strcat(label, term);
	}
	started = 1;
    }
	
    return label;
}

int simulate (char *cmd, double ***pZ, DATAINFO *pdinfo)
     /* implements the "sim" command */
{
    int f, i, t, t1, t2, tstart, m, nv = 0, pv;
    char varname[VNAMELEN], parm[16], tmpstr[MAXLEN];
    char *isconst = NULL, **toks = NULL;
    double xx, yy, *a = NULL;
    int vtok = 0, err = 0;

    *gretl_errmsg = '\0';

    close_minus(cmd);

    f = _count_fields(cmd);
    m = f - 2; /* default: allow for command word varname */

    a = malloc(m * sizeof *a);
    isconst = malloc(m * sizeof *isconst);
    toks = malloc((f - 1) * sizeof *toks);

    if (a == NULL || isconst == NULL || toks == NULL) {
	err = E_ALLOC;
	goto sim_bailout;
    }

    for (i=0; i<m; i++) isconst[i] = 1;

    *tmpstr = 0;
    strncat(tmpstr, cmd, MAXLEN - 1);
    
    strtok(tmpstr, " "); /* discard the "sim" command word */
    for (i=0; i<f-1; i++) {
	toks[i] = strtok(NULL, " ");
    }

    /* allow for implicit starting and ending dates */
    if (isalpha(*toks[0])) {
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
    } else {
	m -= 2;
	vtok = 2;
	/* try getting valid obs from stobs and endobs */
	t1 = dateton(toks[0], pdinfo);
	t2 = dateton(toks[1], pdinfo);
	if (*gretl_errmsg || t1 < 0 || t1 > t2 || t2 >= pdinfo->n) {

	    if (t1 < 0 || t2 >= pdinfo->n) {
		strcpy(gretl_errmsg, _("Observation number out of bounds"));
	    } else if (t1 > t2 ) {
		strcpy(gretl_errmsg, _("Invalid null sample"));
	    }

	    err = 1;
	    goto sim_bailout;
	}
    }

    /* name of var to simulate */
    *varname = 0;
    strncat(varname, toks[vtok], 8);
    nv = varindex(pdinfo, varname);

    if (nv > 0 && nv < pdinfo->v && pdinfo->vector[nv] == 0) {
	sprintf(gretl_errmsg, _("variable %s is a scalar"), 
		pdinfo->varname[nv]);
	err = 1;
	goto sim_bailout;
    }
		
    if (nv == 0 || nv >= pdinfo->v) {
	sprintf(gretl_errmsg, (nv)? _("For 'sim', the variable must already "
				      "exist") :
		_("You can't use the constant for this purpose"));
	err = 1;
	goto sim_bailout;
    }

    /* get the parameter terms */
    for (i=0; i<m; i++) {
	int neg = 0;
	const char *p = parm;

	*parm = '\0';
	strncat(parm, toks[i + vtok + 1], sizeof parm - 1);

	if (*parm == '-') {
	    neg = 1;
	    p++;
	}

	if (isalpha((unsigned char) *p)) {
	    pv = varindex(pdinfo, p);
	    if (pv == 0 || pv >= pdinfo->v) {
		sprintf(gretl_errmsg, _("Bad varname '%s' in sim"), p);
		err = 1;
		goto sim_bailout;
	    } else {
		isconst[i] = !pdinfo->vector[pv];
		a[i] = (isconst[i])? (*pZ)[pv][0] : (double) pv;
	    }
	} else {
	    a[i] = dot_atof(p);
	}

	if (neg) a[i] = -a[i];
    }

    tstart = t1;
    if (tstart < m - 1) tstart = m - 1;

    for (t=tstart; t<=t2; t++) {
	xx = 0.;
	for (i=0; i<m; i++) {
	    if (isconst[i]) {
		if (i == 0) xx += a[i];
		else xx += a[i] * (*pZ)[nv][t-i];
	    } else {
		int neg = 0;

		pv = (int) a[i];
		if (pv < 0) {
		    neg = 1;
		    pv = -pv;
		}
		yy = (*pZ)[pv][t];
		if (na(yy)) {
		    xx = NADBL;
		    break;
		}
		if (neg) yy = -yy;
		if (i == 0) xx += yy;
		else xx += yy * (*pZ)[nv][t-i];
	    }
	}
	(*pZ)[nv][t] = xx;
    }

#if 0 /* this buggers up Jack's OPG loop regressions */
    for (t=t1; t<tstart; t++) {
	(*pZ)[nv][t] = NADBL;
    }
#endif

 sim_bailout:

    if (!err && nv > 0) {
	sprintf(gretl_msg, "%s %s %s (ID %d)", 
		_("Replaced"), _("vector"), pdinfo->varname[nv], nv);
	make_sim_label(VARLABEL(pdinfo, nv), pdinfo->varname[nv],
		       toks + vtok + 1, m);
    }

    free(a);
    free(isconst);
    free(toks);

    return err;
}

/* .......................................................... */

int _multiply (char *s, int *list, char *sfx, double ***pZ,
	       DATAINFO *pdinfo)
{
    int i, t, v = 0, nv, n = pdinfo->n, lv, l0 = list[0];
    int slen;
    double m = 0;
    char tmp[VNAMELEN];

    /* parse s */
    if (isdigit((unsigned char) *s)) m = dot_atof(s);
    else {
	v = varindex(pdinfo, s);
	if (v == pdinfo->v) return E_UNKVAR; 
    }

    if (dataset_add_vars(l0, pZ, pdinfo)) return E_ALLOC;
    slen = strlen(sfx);

    /* fill out values */
    for (i=1; i<=l0; i++) {
	nv = pdinfo->v - l0 - 1 + i;
	lv = list[i];
	for (t=0; t<n; t++) (*pZ)[nv][t] = NADBL;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (na((*pZ)[lv][t])) {
		(*pZ)[nv][t] = NADBL;
		continue;
	    }
	    if (v) {
		double yy = 
		    (pdinfo->vector[v])? 
		    (*pZ)[v][t] : (*pZ)[v][0];

		if (na(yy)) (*pZ)[nv][t] = NADBL;
		else (*pZ)[nv][t] = yy * (*pZ)[lv][t];
	    } else 
		(*pZ)[nv][t] = m * (*pZ)[lv][t];
	}
	/* do names and labels */
	strcpy(tmp, pdinfo->varname[lv]);
	_esl_trunc(tmp, 8 - slen);
	strcat(tmp, sfx);
	strcpy(pdinfo->varname[nv], tmp);
	if (v) 
	    sprintf(VARLABEL(pdinfo, nv), "%s = %s * %s",
		    pdinfo->varname[nv], pdinfo->varname[v], 
		    pdinfo->varname[lv]); 
	else 
	    sprintf(VARLABEL(pdinfo, nv), "%s = %g * %s",
		    pdinfo->varname[nv], m, pdinfo->varname[lv]); 
    }

    return 0;
}

/**
 * genr_fit_resid:
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @code: GENR_RESID or GENR_FITTED or GENR_RESID2.
 * @undo: if non-zero, don't bother labeling the variables
 *
 * Adds residuals or fitted values or squared residuals from a
 * given model to the data set.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int genr_fit_resid (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int code, int undo)
{
    char vname[VNAMELEN], vlabel[MAXLABEL];
    int i, n, t, t1 = pmod->t1, t2 = pmod->t2;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    i = pdinfo->v - 1;
    n = pdinfo->n;

    if (pmod->data != NULL) t2 += get_misscount(pmod);

    for (t=0; t<t1; t++) (*pZ)[i][t] = NADBL;
    for (t=t2+1; t<n; t++) (*pZ)[i][t] = NADBL;

    if (code == GENR_RESID) { /* residuals */
	sprintf(vname, "uhat%d", pmod->ID);
	sprintf(vlabel, _("residual from model %d"), pmod->ID);
	for (t=t1; t<=t2; t++) 
	    (*pZ)[i][t] = pmod->uhat[t];
    }
    else if (code == GENR_FITTED) { /* fitted values */
	sprintf(vname, "yhat%d", pmod->ID);
	sprintf(vlabel, _("fitted value from model %d"), pmod->ID);
	for (t=t1; t<=t2; t++) 
	    (*pZ)[i][t] = pmod->yhat[t];
    }
    else if (code == GENR_RESID2) { /* squared residuals */
	sprintf(vname, "usq%d", pmod->ID);
	sprintf(vlabel, _("squared residual from model %d"), pmod->ID);
	for (t=t1; t<=t2; t++) 
	    (*pZ)[i][t] = pmod->uhat[t] * pmod->uhat[t];
    }
    strcpy(pdinfo->varname[i], vname);

    if (!undo) 
	strcpy(VARLABEL(pdinfo, i), vlabel);

    return 0;
}


    
