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

#include "libgretl.h"
#include "gretl_func.h"
#include "loop_private.h"
#include "genstack.h"
#include "genrfuncs.h"
#include "libset.h"
#include "modelspec.h"
#include "objstack.h"
#include "usermat.h"

#include "../../cephes/libprob.h"

#include <time.h>
#include <errno.h>

enum {
    HIGHNUM = 5000,
    TNUM,
    INDEXNUM
} genr_numbers;

enum {
    E_MATVAR = E_MAX + 1, /* matrix found where regular var was expected */
    E_SERIES              /* required scalar but got series */
} genr_errors;

static double calc_xy (double x, double y, char op, int t, int *err);

static int genr_mpow (const char *str, double *xvec, double **Z, 
		      DATAINFO *pdinfo);

#ifdef HAVE_MPFR
static int genr_mlog (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
#endif

static void eval_msg (const GENERATOR *genr);
static void compose_genr_msg (const GENERATOR *genr, int oldv);

static void genr_deposit_result (GENERATOR *genr);
static int genr_write_var (GENERATOR *genr, double ***pZ);
static int genr_retrieve_result (double *x, gretl_matrix **m);

static int math_tokenize (char *s, GENERATOR *genr, int level);
static int get_obs_value (const char *s, const double **Z, 
			  const DATAINFO *pdinfo, int *v, int *obsnum);

static int op_level (int c);

static double *get_random_series (DATAINFO *pdinfo, int fn);
static double *get_mp_series (const char *s, GENERATOR *genr,
			      int fn, int *err);
static double *get_tmp_series (double *x, GENERATOR *genr,
			       int fn, double param);
static int genr_add_temp_var (double *x, GENERATOR *genr, genatom *atom);

static double get_tnum (const DATAINFO *pdinfo, int t);
static double evaluate_statistic (double *z, GENERATOR *genr, int fn);
static double get_dataset_statistic (DATAINFO *pdinfo, int idx);
static double get_test_stat_value (char *label, int idx);
static double evaluate_math_function (double arg, int fn, int *err);
static double evaluate_missval_func (double arg, int fn);
static double evaluate_bivariate_statistic (const char *s, 
					    GENERATOR *genr, 
					    int fn);
static double evaluate_pvalue (const char *s, const double **Z,
			       const DATAINFO *pdinfo, int *err);
static double evaluate_critval (const char *s, const double **Z,
				const DATAINFO *pdinfo, int *err);
static void free_genr_S (GENERATOR *genr);

enum retrieve {
    R_NOBS = 1,  /* number of observations in current sample range */
    R_NVARS,     /* number of variables in dataset (including the constant) */
    R_PD,        /* periodicity of dataset */
    R_TEST_STAT, /* test statistic from last explicit test performed */
    R_TEST_PVAL  /* p-value from last explicit test performed */
};

struct genr_func {
    int fnum;
    const char *fword;
};

struct genr_func funcs[] = {
    { T_LOG,      "log" }, 
    { T_EXP,      "exp" }, 
    { T_SIN,      "sin" }, 
    { T_COS,      "cos" }, 
    { T_TAN,      "tan" },
    { T_ATAN,     "atan" },
    { T_DIFF,     "diff" },
    { T_LDIFF,    "ldiff" }, 
    { T_SDIFF,    "sdiff" },
    { T_MEAN,     "mean" }, 
    { T_SD,       "sd" }, 
    { T_MIN,      "min" },
    { T_MAX,      "max" },
    { T_SORT,     "sort" }, 
    { T_DSORT,    "dsort" }, 
    { T_INT,      "int" }, 
    { T_ABS,      "abs" }, 
    { T_SQRT,     "sqrt" }, 
    { T_SUM,      "sum" }, 
    { T_NOBS,     "nobs" },
    { T_T1,       "firstobs" },
    { T_T2,       "lastobs" },
    { T_NORMAL,   "normal" }, 
    { T_UNIFORM,  "uniform" }, 
    { T_CUM,      "cum" }, 
    { T_MISSING,  "missing" },
    { T_OK,       "ok" },        /* opposite of missing */
    { T_MISSZERO, "misszero" },
    { T_CORR,     "corr" },
    { T_VAR,      "var" },
    { T_LRVAR,    "lrvar" },
    { T_SST,      "sst" },
    { T_COV,      "cov" },
    { T_MEDIAN,   "median" },
    { T_GINI,     "gini" },
    { T_ZEROMISS, "zeromiss" },
    { T_PVALUE,   "pvalue" },
    { T_CRIT,     "critical" },
    { T_OBSNUM,   "obsnum" },
    { T_MPOW,     "mpow" },
    { T_DNORM,    "dnorm" },
    { T_CNORM,    "cnorm" },
    { T_QNORM,    "qnorm" },
    { T_GAMMA,    "gamma" },
    { T_LNGAMMA,  "lngamma" },
    { T_RESAMPLE, "resample" },
    { T_HPFILT,   "hpfilt" },    /* Hodrick-Prescott filter */
    { T_BKFILT,   "bkfilt" },    /* Baxter-King filter */
    { T_FRACDIFF, "fracdiff" },  /* fractional difference */
    { T_VARNUM,   "varnum" },    /* variable's ID number from its name */
    { T_SERIES,   "isseries" },
    { T_ISLIST,   "islist" },
    { T_NELEM,    "nelem" },
    { T_DET,      "det" },
    { T_INV,      "inv" },
    { T_CHOL,     "cholesky" },
    { T_CDMEAN,   "cdemean" },
    { T_DIAG,     "diag" },
    { T_QR,       "qrdecomp" },
    { T_EIGSYM,   "eigensym" },
    { T_EIGGEN,   "eigengen" },
    { T_LDET,     "ldet" },
    { T_TRACE,    "tr" },
    { T_SUMR,     "sumr" },
    { T_SUMC,     "sumc" },
    { T_1NORM,    "onenorm" },
    { T_RCOND,    "rcond" },
    { T_ROWS,     "rows" },
    { T_COLS,     "cols" },
    { T_TRANSP,   "transp" },
    { T_IMAT,     "I" },
    { T_ZEROS,    "zeros" },
    { T_ONES,     "ones" },
#ifdef HAVE_MPFR
    { T_MLOG,     "mlog" },
#endif
    { T_IDENTITY, "ident" },
    { 0, NULL }
};

#define LEVELS 7

#define STANDARD_MATH(f) (f == T_LOG || f == T_EXP || \
                          f == T_SIN || f == T_COS || f == T_TAN || \
                          f == T_ATAN || f == T_INT || f == T_ABS || \
                          f == T_DNORM || f == T_CNORM || f == T_QNORM || \
                          f == T_SQRT || f == T_GAMMA || f == T_LNGAMMA)

#define UNIVARIATE_STAT(t) (t == T_MEAN || t == T_SD || t == T_SUM || \
                            t == T_VAR || t == T_MEDIAN || t == T_MIN || \
                            t == T_SST || t == T_MAX || t == T_NOBS || \
                            t == T_T1 || t == T_T2 || t == T_VARNUM || \
                            t == T_SERIES || t == T_ISLIST || t == T_NELEM || \
                            t == T_GINI || t == T_LRVAR)

#define BIVARIATE_STAT(t) (t == T_CORR || t == T_COV)

#define MISSVAL_FUNC(f) (f == T_MISSING || f == T_OK || \
                         f == T_MISSZERO || f == T_ZEROMISS)

#define MATRIX_SCALAR_FUNC(f) (f == T_DET || f == T_LDET || f == T_TRACE || \
			       f == T_ROWS || f == T_COLS || f == T_1NORM || \
                               f == T_RCOND)

#define MATRIX_FILL_FUNC(f) (f == T_IMAT || f == T_ZEROS || f == T_ONES)

#define MATRIX_MATRIX_FUNC(f) (f == T_TRANSP || f == T_DIAG || f == T_CDMEAN || \
                               f == T_INV || f == T_CHOL || f == T_QR || \
                               f == T_SUMR || f == T_SUMC || f == T_EIGSYM || \
			       f == T_EIGGEN)

#define MULTI_MATRIX_FUNC(f) (f == T_QR || f == T_EIGSYM || f == T_EIGGEN)

#define RAND_FUNC(f) (f == T_UNIFORM || f == T_NORMAL)

#define needs_arg(f) (f > 0 && f != T_UNIFORM && f != T_NORMAL)

#ifdef HAVE_MPFR
# define MP_MATH(f) (f == T_MPOW || f == T_MLOG)
#else
# define MP_MATH(f) (f == T_MPOW)
#endif

#define MAXTERMS  64
#define TOKLEN  1024
#define ARGLEN   128

#define genr_is_matrix(g) ((g)->flags & GENR_MATRIX)
#define genr_set_matrix(g) ((g)->flags |= GENR_MATRIX)
#define genr_is_submatrix(g) (((g)->flags & GENR_MATRIX) && *g->label)

#define genr_is_scalar(g) ((g)->flags & GENR_SCALAR)
#define genr_is_series(g) (!((g)->flags & GENR_SCALAR))
#define genr_set_scalar(g) ((g)->flags |= GENR_SCALAR)
#define genr_unset_scalar(g) ((g)->flags &= ~GENR_SCALAR)

#define genr_force_series(g) ((g)->flags |= GENR_FORCE_SERIES)
#define genr_forcing_series(g) ((g)->flags & GENR_FORCE_SERIES)

#define genr_set_require_scalar(g) ((g)->flags |= GENR_NEED_SCALAR)
#define genr_require_scalar(g) ((g)->flags & GENR_NEED_SCALAR)

#define genr_doing_save(g) ((g)->flags & GENR_SAVE)
#define genr_set_save(g) ((g)->flags |= GENR_SAVE)
#define genr_unset_save(g) ((g)->flags &= ~GENR_SAVE)

#define genr_is_private(g) ((g)->flags & GENR_PRIVATE)
#define genr_set_private(g) ((g)->flags |= GENR_PRIVATE)

#define genr_warn(g) ((g)->flags & GENR_WARN)
#define genr_set_warn(g) ((g)->flags |= GENR_WARN)

#define genr_simple_sort(g) ((g)->flags & GENR_SIMPLE_SORT)
#define genr_set_simple_sort(g) ((g)->flags |= GENR_SIMPLE_SORT)

#define genr_depositing_result(g) ((g)->flags & GENR_DEPOSIT)
#define genr_set_depositing(g) ((g)->flags |= GENR_DEPOSIT)

#define genr_single_obs(g) ((g)->obs >= 0)

#define set_atom_is_series_from_matrix(a) (strcpy(a->str, "matrix_series"))
#define atom_is_series_from_matrix(a) (!strcmp(a->str, "matrix_series"))

static int genr_matrix_init (GENERATOR *genr, gretlopt opt)
{
    int i, err = 0;

    genr->mstack = NULL;
    genr->nmats = 0;

    if (opt & OPT_M) {
	genr->mstack = malloc(MATSTACK_SIZE * sizeof *genr->mstack);
	if (genr->mstack == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<MATSTACK_SIZE; i++) {
		genr->mstack[i] = NULL;
	    }
	}
    }

    return err;
}

static GENERATOR *genr_new (double ***pZ, DATAINFO *pdinfo, gretlopt opt,
			    PRN *prn)
{
    GENERATOR *genr;

    genr = malloc(sizeof *genr);
    if (genr == NULL) {
	return NULL;
    }

    genr->err = 0;
    genr->done = 0;

    if (opt & OPT_M) {
	genr->flags = GENR_SAVE | GENR_MATRIX;
    } else {
	genr->flags = GENR_SAVE | GENR_SCALAR;
    }

    genr->prn = prn;

    genr->xvec = NULL;
    genr->varnum = 0;
    genr->obs = -1;

    *genr->varname = '\0';
    *genr->label = '\0';
    *genr->orig_s = '\0';

    genr->tmpv = 0;
    genr->tmpZ = NULL;

    genr->aset = NULL;
    genr->S = NULL;

    if (opt & OPT_P) {
	genr_set_private(genr);
    }

    if (opt & OPT_D) {
	genr_set_depositing(genr);
    }

    /* hook up pointers to "outside world" */
    genr->pdinfo = pdinfo;
    genr->pZ = pZ;

    reset_calc_stack(genr);

    if (genr_matrix_init(genr, opt)) {
	free(genr);
	genr = NULL;
    }

    return genr;
}

static genatom *atom_new (int level, char op, atomset *aset)
{
    genatom *atom = malloc(sizeof *atom);

    if (atom == NULL) {
	return NULL;
    }

    atom->atype = ATOM_SERIES;
    atom->varnum = -1;
    atom->varobs = -1;
    atom->tmpvar = -1;
    atom->lag = 0;
    atom->val = 0.0;
    atom->M = NULL;
    atom->func = 0;

    atom->op = op;
    atom->level = level;

    *atom->str = '\0';
    atom->parent = NULL;
    atom->aset = aset;
    atom->popped = 0;

    return atom;
}

#if GENR_DEBUG
static int print_atom (genatom *atom)
{
    if (atom == NULL) return 1;

    fprintf(stderr, " atom->level = %d\n", atom->level);

    fprintf(stderr, " atom->atype = %d\n", atom->atype);

    if (atom->varnum != 0) {
	fprintf(stderr, " atom->varnum = %d\n", atom->varnum);
    }
    if (atom->varobs >= 0) {
	fprintf(stderr, " atom->varobs = %d\n", atom->varobs);
    }
    if (atom->tmpvar >= 0) {
	fprintf(stderr, " atom->tmpvar = %d\n", atom->tmpvar);
    }
    if (atom->lag > 0) {
	fprintf(stderr, " atom->lag = %d\n", atom->lag);
    }
    if (atom->atype & ATOM_SCALAR) {
	fprintf(stderr, " atom->val = %g\n", atom->val);
    }
    if (atom->M != NULL) {
	fprintf(stderr, " atom->M = %p (int = %d)\n", (void *) atom->M,
		gretl_matrix_get_int(atom->M));
    }
    if (atom->func > 0) {
	const char *fword = get_genr_func_word(atom->func);

	if (fword != NULL) {
	    fprintf(stderr, " atom->func = %d (%s)\n", 
		    atom->func, get_genr_func_word(atom->func));
	} else {
	    fprintf(stderr, " atom->func = %d (UNKNOWN!)\n", 
		    atom->func);
	}
    }

    if (*atom->str) {
	fprintf(stderr, " atom->str = '%s'\n", atom->str);
    }    

    if (isprint(atom->op)) {
	fprintf(stderr, " atom->op = '%c'\n", atom->op);
    } else {
	fprintf(stderr, " atom->op = %d\n", atom->op);
    }

    return 0;
}
#endif

static int get_lagvar (const char *s, int *lag, GENERATOR *genr)
{
    static char format[16] = {0};
    const char *p;
    char vname[VNAMELEN];
    char ls[8];
    int m = 0, v = 0;

    /* an atomic lagvar must end with the first right paren */
    if ((p = strchr(s, ')')) != NULL && *(p+1) != 0) {
	return 0;
    }

    if (*format == 0) {
	sprintf(format, "%%%d[^(](%%7[^)])", VNAMELEN - 1);
    }

    DPRINTF(("get_lagvar: looking at '%s'\n", s));

    if (sscanf(s, format, vname, ls) == 2 && sscanf(ls, "%d", &m) == 1) {
	v = varindex(genr->pdinfo, vname);
	DPRINTF(("get_lagvar: vname = '%s' (v == %d)\n", vname, v));
	if (v >= genr->pdinfo->v) {
	    DPRINTF(("get_lagvar: rejecting '%s'\n", s));
	    v = m = 0;
	} else if (v < genr->pdinfo->v && genr->pdinfo->vector[v] == 0) {
	    sprintf(gretl_errmsg, _("Variable %s is a scalar; "
				    "can't do lags/leads"), 
		    genr->pdinfo->varname[v]);
	    genr->err = E_DATA;
	    v = m = 0;
	} 
    }

    if (lag == NULL) {
	/* just testing for atomicity */
	if (v > 0) {
	    char *test;

	    strtol(ls, &test, 10);
	    if (*test != 0) {
		/* string ls contains superfluous stuff */
		v = 0;
	    }
	}
    } else {
	*lag = -m;
    }

    DPRINTF(("get_lagvar: v = %d, m = %d\n", v, m));

    return v;
}

/* also used in do_printf function */

int get_generated_value (const char *rhs, double *val,
			 double ***pZ, DATAINFO *pdinfo,
			 int t)
{
    char genline[MAXLINE];
    int v, err = 0;

    if (numeric_string(rhs)) {
	*val = dot_atof(rhs);
	return 0;
    }

    if (t < pdinfo->n) {
	v = varindex(pdinfo, rhs);
	if (v < pdinfo->v) {
	    *val = (*pZ)[v][t];
	    return 0;
	}
    }

    sprintf(genline, "genr x___=%s", rhs);
#if GENR_DEBUG
    fprintf(stderr, "get_generated_value: trying '%s' (t=%d)\n", genline, t);
#endif
    err = generate(genline, pZ, pdinfo, OPT_P | OPT_D, NULL);
    if (!err) {
	err = genr_retrieve_result(val, NULL);
    }

    return err;
}

static int get_generated_matrix (const char *rhs, gretl_matrix **M,
				 double ***pZ, DATAINFO *pdinfo)
{
    char genline[MAXLINE];
    int err = 0;

    sprintf(genline, "genr M___=%s", rhs);

#if GENR_DEBUG
    fprintf(stderr, "get_generated_matrix: trying '%s'\n", genline);
#endif

    err = generate(genline, pZ, pdinfo, OPT_P | OPT_M | OPT_D, NULL);
    if (!err) {
	err = genr_retrieve_result(NULL, M);
    }

    return err;
}

/* get the values of the scalar arguments to a genr function */

static int evaluate_genr_function_args (char *s, GENERATOR *genr)
{
    char *tmp = gretl_strdup(s);
    char *st, *arg;
    double vals[3];
    int i, nf;
    int err = 0;

    if (tmp == NULL) {
	return 1;
    }

    st = strtok(tmp, ",");
    if (st == NULL) {
	free(tmp);
	return 1;
    }

    nf = 0;
    for (i=0; i<3 && !err; i++) {
	arg = strtok(NULL, ",");
	if (arg == NULL) break;
	if (numeric_string(arg)) {
	    vals[i] = dot_atof(arg);
	} else {
	    err = get_generated_value(arg, &vals[i], genr->pZ, 
				      genr->pdinfo, 0);
	}
	nf++;
    }

    if (!err) {
	char numstr[32];

	*s = '\0';
	strcat(s, st);
	gretl_push_c_numeric_locale();
	for (i=0; i<nf; i++) {
	    sprintf(numstr, ",%.15g", vals[i]);
	    strcat(s, numstr);
	}
	gretl_pop_c_numeric_locale();
    }

    free(tmp);
	
    return err;
}

static int 
set_atom_arg_string (const char *s, GENERATOR *genr, genatom *atom)
{
    const char *p = strchr(s, '(');
    int n = strlen(p) - 2;
    int err = 0;

    DPRINTF(("set_atom_arg_string: looking at '%s'\n", s));

    if (n < 0 || n > ARGLEN - 1) {
	err = 1;
    } else {
	strncat(atom->str, p + 1, n);
    }

    if (atom->func == T_PVALUE || atom->func == T_CRIT) {
	err = evaluate_genr_function_args(atom->str, genr);
    } else if (atom->func == T_FRACDIFF) {
	char vname[9];
	double param;

	err = sscanf(atom->str, "%15[^,],%lf", vname, &param) != 2;
    }

    DPRINTF(("get_arg_string: got '%s'\n", atom->str));

    return err;
}

static gretl_matrix *genr_get_matrix_by_name (const char *s, GENERATOR *genr)
{
    gretl_matrix *m = atom_stack_get_matrix(genr, s);

    if (m != NULL) {
	DPRINTF(("recognized atom matrix '%s'\n", s));
    } else {
	m = get_matrix_by_name(s, genr->pdinfo);
	if (m != NULL) {
	    DPRINTF(("recognized matrix '%s'\n", s));
	} else {
	    m = get_matrix_transpose_by_name(s, genr->pdinfo);
	    if (m != NULL) {
		DPRINTF(("recognized transposed matrix '%s'\n", s));
	    }
	}
    }

    return m;
}

static double 
genr_get_matrix_scalar (const char *s, GENERATOR *genr, int func)
{
    const gretl_matrix *m = NULL;
    gretl_matrix *g = NULL;
    char mstr[128];
    double x = NADBL;

    s = strchr(s, '(');
    
    if (s == NULL || *(s+1) == 0) {
	goto bailout;
    }

    /* simple argument to matrix function? */
    sscanf(s+1, "%31[^)]", mstr);
    m = genr_get_matrix_by_name(mstr, genr);

    if (m == NULL) {
	/* compound argument to matrix function? */
	int len = strlen(s+1) - 1;
	    
	if (len > 0 && len < 128) {
	    *mstr = 0;
	    strncat(mstr, s+1, len);
	    genr->err = get_generated_matrix(mstr, &g, genr->pZ, genr->pdinfo);
	    m = g;
	}
    }

    if (m != NULL) {
	if (func == T_DET) {
	    x = user_matrix_get_determinant(m, &genr->err);
	} else if (func == T_LDET) {
	    x = user_matrix_get_log_determinant(m, &genr->err);
	} else if (func == T_TRACE) {
	    x = gretl_matrix_trace(m, &genr->err);
	} else if (func == T_ROWS) {
	    x = gretl_matrix_rows(m);
	} else if (func == T_COLS) {
	    x = gretl_matrix_cols(m);
	} else if (func == T_1NORM) {
	    x = gretl_matrix_one_norm(m);
	} else if (func == T_RCOND) {
	    x = gretl_symmetric_matrix_rcond(m, &genr->err);
	}
    } 

    if (g != NULL) {
	gretl_matrix_free(g);
    }

 bailout:

    if (na(x)) {
	sprintf(gretl_errmsg, "%s%s: error in function argument", 
		get_genr_func_word(func), s);
	if (genr->err == 0) {
	    genr->err = 1;
	}
    }

    return x;
}

static void undefined_var_message (const char *s)
{
    sprintf(gretl_errmsg, _("Undefined variable name '%s' in genr"), s);
}

static void
get_var_from_matrix (const char *s, GENERATOR *genr, genatom *atom)
{
    double *x = NULL;
    int len = 0;

    /* We have to be careful with this: we'll try conversion from
       matrix to variable (scalar or series) only if the command
       has explicitly typed the desired result as scalar or
       series.
    */
    if (!genr_require_scalar(genr) && !genr_forcing_series(genr)) {
	undefined_var_message(s);
	genr->err = E_UNKVAR;
    }

    /* generating a scalar: we'll only accept 1 x 1 matrices */
    if (genr_require_scalar(genr)) {
	len = 1;
    }

    genr->err = named_matrix_get_variable(s, genr->pZ, genr->pdinfo, 
					  &x, &len);

    if (!genr->err) {
	if (len == 1) {
	    atom->atype = ATOM_SCALAR;
	    atom->val = x[0];
	    free(x);
	} else {
	    genr_add_temp_var(x, genr, atom);
	    set_atom_is_series_from_matrix(atom);
	}
    } else if (genr_require_scalar(genr) && genr->err != E_UNKVAR) {
	/* Flag the possibility that the RHS is a compound
	   matrix expression that evaluates to a 1 x 1 matrix,
	   which we can then convert to scalar.
	*/
	genr->err = E_MATVAR;
    }
}

static void atom_set_matrix (genatom *atom, gretl_matrix *M, const char *s)
{
    atom->M = M;
    atom->atype = ATOM_MATRIX;
    set_matrix_on_atom(atom->M);
    if (s != NULL) {
	strncat(atom->str, s, ATOMLEN - 1); 
    }
}

static void
atom_get_variable_or_constant (const char *s, GENERATOR *genr,
			       genatom *atom)
{
    gretl_matrix *M = NULL;
    int v;

    if (!strcmp(s, "pi")) {
	atom->val = M_PI;
	atom->atype = ATOM_SCALAR;
    } else if (!strcmp(s, "NA")) {
	atom->val = NADBL;
	atom->atype = ATOM_SCALAR;
    } else if (genr_is_matrix(genr)) {
	M = genr_get_matrix_by_name(s, genr);
	if (M != NULL) {
	    atom_set_matrix(atom, M, s);
	} else {	
	    /* try for a scalar? */
	    v = varindex(genr->pdinfo, s);
	    if (v < genr->pdinfo->v && !genr->pdinfo->vector[v]) {
		atom->val = (*genr->pZ)[v][0];
		atom->atype = ATOM_SCALAR;
	    } else {
		sprintf(gretl_errmsg, _("'%s': not a matrix or scalar"), s);
		genr->err = E_UNKVAR;
	    }
	} 
    } else {
	v = varindex(genr->pdinfo, s);

	if (v == genr->pdinfo->v) { 
	    if (get_list_by_name(s) != NULL) {
		/* name of saved list, not single variable */
		strncat(atom->str, s, ATOMLEN - 1);
	    } else {
		get_var_from_matrix(s, genr, atom);
		if (genr->err && genr->err != E_MATVAR) {
		    undefined_var_message(s);
		    genr->err = E_UNKVAR;
		}
	    }
	} else {
	    DPRINTF(("recognized var '%s' (#%d)\n", s, v));
	    atom->varnum = v;
	}

	if (v == INDEXNUM) { 
	    int k = loop_scalar_read(*s);

	    atom->val = k;
	    atom->atype = ATOM_SCALAR;
	} else if (v >= 0 && v < genr->pdinfo->v && !genr->pdinfo->vector[v]) {
	    /* regular scalar variables are handled here */
	    atom->varnum = v;
	    atom->varobs = 0;
	    atom->atype = ATOM_SCALAR;
	}
    }
}

static int get_matrix_dim (const char *s, GENERATOR *genr)
{
    return gretl_int_from_string(s, (const double **) *genr->pZ,
				 genr->pdinfo, &genr->err);
}

static void
matrix_gen_function (const char *s, GENERATOR *genr, genatom *atom)
{
    gretl_matrix *M = NULL;
    char rstr[VNAMELEN], cstr[VNAMELEN];
    int r = 0, c = 0;
    int nf = 0;

    s = strchr(s, '(') + 1;

    if (sscanf(s, "%15[^,],%15[^)]", rstr, cstr) == 2) {
	r = get_matrix_dim(rstr, genr);
	c = get_matrix_dim(cstr, genr);
	if (r <= 0 || c <= 0) {
	    genr->err = E_SYNTAX;
	} else {
	    nf = 2;
	}
    } else if (sscanf(s, "%15[^)]", rstr)) {
	r = get_matrix_dim(rstr, genr);
	if (r <= 0) {
	    genr->err = E_SYNTAX;
	} else {
	    nf = 1;
	}
    }

    if (genr->err) {
	return;
    }

    if ((atom->func == T_IMAT && nf != 1) || 
	(atom->func != T_IMAT && nf != 2)) {
	genr->err = 1;
	sprintf(gretl_errmsg, "wrong number of arguments for function %s()\n",
		get_genr_func_word(atom->func));
    } else {
	if (atom->func == T_IMAT) {
	    M = gretl_identity_matrix_new(r);
	} else if (atom->func == T_ZEROS) {
	    M = gretl_zero_matrix_new(r, c);
	} else if (atom->func == T_ONES) {
	    M = gretl_unit_matrix_new(r, c);
	} else if ((atom->func == T_UNIFORM || atom->func == T_NORMAL)) {
	    M = gretl_matrix_alloc(r, c);
	    if (M != NULL) {
		gretl_matrix_random_fill(M, atom->func);
	    }
	}
	if (M == NULL) {
	    genr->err = E_ALLOC;
	} else {
	    atom_set_matrix(atom, M, s);
	} 
    }
}

static void
atom_get_function_data (const char *s, GENERATOR *genr, genatom *atom)
{
    DPRINTF(("recognized function #%d (%s)\n", atom->func, 
	     get_genr_func_word(atom->func)));

    if (MP_MATH(atom->func) || atom->func == T_PVALUE || 
	atom->func == T_CRIT || atom->func == T_FRACDIFF ||
	BIVARIATE_STAT(atom->func)) {
	genr->err = set_atom_arg_string(s, genr, atom);
    } else if (MULTI_MATRIX_FUNC(atom->func) && genr_is_matrix(genr)) {
	genr->err = set_atom_arg_string(s, genr, atom);
    }

    if (atom->func == T_PVALUE) {
	if (!genr->err) {
	    atom->val = evaluate_pvalue(atom->str, (const double **) *genr->pZ,
					genr->pdinfo, &genr->err);
	    atom->atype = ATOM_SCALAR;
	    atom->func = 0;
	}
    } else if (atom->func == T_CRIT) {
	if (!genr->err) {
	    atom->val = evaluate_critval(atom->str, (const double **) *genr->pZ,
					 genr->pdinfo, &genr->err);
	    atom->atype = ATOM_SCALAR;
	    atom->func = 0;
	}
    } else if (atom->func == T_VARNUM || atom->func == T_SERIES ||
	       atom->func == T_ISLIST || atom->func == T_NELEM) {
	atom->atype = ATOM_SCALAR;
    } else if (BIVARIATE_STAT(atom->func)) {
	atom->val = evaluate_bivariate_statistic(atom->str, genr, atom->func);
	atom->atype = ATOM_SCALAR;
	atom->func = 0;
    } else if (MATRIX_SCALAR_FUNC(atom->func) && !genr_is_matrix(genr)) {
	atom->val = genr_get_matrix_scalar(s, genr, atom->func);
	atom->atype = ATOM_SCALAR;
	atom->func = 0;
    } else if (RAND_FUNC(atom->func) && genr_is_matrix(genr)) {
	matrix_gen_function(s, genr, atom);
    } else if (MATRIX_FILL_FUNC(atom->func)) {
	if (genr_is_matrix(genr)) {
	    matrix_gen_function(s, genr, atom);
	} else if (genr_require_scalar(genr)) {
	    genr->err = E_MATVAR; /* allow second pass in matrix mode */
	} else {
	    genr->err = E_TYPES;
	}
    } else if (MATRIX_MATRIX_FUNC(atom->func) && !genr_is_matrix(genr)) {
	if (genr_require_scalar(genr)) {
	    genr->err = E_MATVAR; /* allow second pass */
	} else {
	    genr->err = E_TYPES;
	}
    }
}

static int dataset_var_index (const char *s)
{
    char test[VNAMELEN];
    int ret = 0;

    *test = '\0';
    strncat(test, s, VNAMELEN - 1);
    lower(test);

    if (!strcmp(test, "$nobs")) {
	ret = R_NOBS;
    } else if (!strcmp(test, "$pd")) {
	ret = R_PD;
    } else if (!strcmp(test, "$nvars")) {
	ret = R_NVARS;
    }

    return ret;
}

static int test_stat_index (const char *s)
{
    char test[VNAMELEN];
    int ret = 0;

    *test = '\0';
    strncat(test, s, VNAMELEN - 1);
    lower(test);

    if (!strcmp(test, "$pvalue")) { 
	ret = R_TEST_PVAL;
    } else if (!strcmp(test, "$test")) {
	ret = R_TEST_STAT;
    }

    return ret;
}

static void 
atom_get_dollar_var (const char *s, GENERATOR *genr, genatom *atom) 
{
    int i, found = 0;

    DPRINTF(("atom_get_dollar_var: s='%s'\n", s));

    if (*s == '$') {
	if ((i = dataset_var_index(s)) > 0) {
	    DPRINTF(("recognized '%s' as dataset var, index #%d\n", s, i));
	    atom->val = get_dataset_statistic(genr->pdinfo, i);
	    atom->atype = ATOM_SCALAR;
	    found = 1;
	} else if ((i = test_stat_index(s)) > 0) {
	    DPRINTF(("recognized '%s' as test-related var, index #%d\n", s, i));
	    atom->val = get_test_stat_value(genr->label, i);
	    atom->atype = ATOM_SCALAR;
	    found = 1;
	}
    }

    if (!found) {
	const char *key = s;
	char *oname = NULL;
	int idx = 0;

	if (*s != '$') {
	    key = strstr(s, ".$");
	    oname = gretl_strndup(s, key - s);
	    key++;
	}

	DPRINTF(("dollar var: key='%s'\n", key));

	idx = gretl_model_data_index(key);

	if (genr_is_matrix(genr)) {
	    gretl_matrix *M = NULL;

	    if (model_data_is_scalar(idx)) {
		double x = saved_object_get_scalar(oname, key, &genr->err);

		if (!genr->err) {
		    M = gretl_matrix_from_scalar(x);
		} 
		if (M == NULL) {
		    genr->err = E_ALLOC;
		} else {
		    atom_set_matrix(atom, M, s);
		    found = 1;
		}
	    } else if (model_data_is_matrix(idx)) {
		M = saved_object_get_matrix(oname, key, 
					    genr->pZ, genr->pdinfo,
					    &genr->err);
		if (!genr->err) {
		    atom_set_matrix(atom, M, s);
		    found = 1;
		}
	    }
	} else {
	    /* not generating a matrix */
	    if (model_data_is_scalar(idx)) {
		atom->val = saved_object_get_scalar(oname, key, &genr->err);
		atom->atype = ATOM_SCALAR;
		found = 1;
	    } else if (model_data_is_scalar_element(idx)) {
		atom->val = saved_object_get_scalar_element(oname, key, genr->pdinfo,
							    &genr->err);
		if (!genr->err) {
		    atom->atype = ATOM_SCALAR;
		    found = 1;
		}
	    } else if (model_data_is_series(idx)) {
		atom->atype = ATOM_SERIES | ATOM_MODELDAT;
		strncat(atom->str, s, ATOMLEN - 1);
		found = 1;
	    }
	}

	free(oname);
    }

    if (!found && !genr->err) {
	genr->err = E_UNKVAR;
    }
}

#define looks_like_dollar_token(s) (*s == '$' || \
                                    (isalpha(*s) && strstr(s, ".$")))

#define looks_like_var_or_func(s) (isalpha(*s) || *s == '\'')

#define looks_like_plain_var(s) (!strchr(s, '(') && !strchr(s, '[') && *s != '\'')

/* Parse a "genr" token and construct a corresponding "atom".  In the
   case of scalar tokens we could go ahead and attach a specific value
   to the atom, but we defer this so far as possible.  We're currently
   "compiling" a genr expression, and we may wish to "run" it more
   than once, with differing values of the inputs -- so we prefer to
   record the info required to obtain the scalar value, where
   possible.
*/

static genatom *
token_make_atom (const char *s, char op, GENERATOR *genr, int level)
{
    genatom *atom = NULL;
    int v, lag;

    DPRINTF(("token_make_atom: looking at '%s'\n", s));

    atom = atom_new(level, op, genr->aset);
    if (atom == NULL) {
	genr->err = E_ALLOC;
	return NULL;
    }
    
    if (looks_like_dollar_token(s)) {
	/* dollarized (internal) variable? */
	atom_get_dollar_var(s, genr, atom);
    } else if (looks_like_var_or_func(s)) {
	if (looks_like_plain_var(s)) {
	    /* ordinary variable, list or matrix? */
	    atom_get_variable_or_constant(s, genr, atom);
	} else if (strchr(s, '(')) {
	    /* function or lagged variable? */
	    if ((atom->func = genr_function_from_string(s))) {
		atom_get_function_data(s, genr, atom);
	    } else {
		/* not a function: try a lagged variable */
		v = get_lagvar(s, &lag, genr);
		if (v > 0) {
		    DPRINTF(("recognized var #%d lag %d\n", v, lag));
		    atom->varnum = v;
		    atom->lag = lag;
		} else {
		    /* not a variable or function: dead end? */
		    DPRINTF(("dead end in token_make_atom, s='%s'\n", s));
		    genr->err = E_SYNTAX; 
		}
	    }
	} else if (strchr(s, '[')) {
	    /* specific observation from a series, or slice of matrix? */
	    gretl_matrix *M = NULL;
	    int v = -1, varobs = -1;
	    int err;

	    err = get_obs_value(s, (const double **) *genr->pZ, genr->pdinfo,
				&v, &varobs);
	    if (!err) {
		if (genr_is_matrix(genr)) {
		    M = gretl_matrix_from_scalar((*genr->pZ)[v][varobs]);
		    atom_set_matrix(atom, M, s);
		} else {
		    atom->varnum = v;
		    atom->varobs = varobs;
		    atom->atype = ATOM_SCALAR;
		}
	    } else if (genr_is_matrix(genr)) {
		M = user_matrix_get_slice(s, genr->pZ, genr->pdinfo, 
					  &genr->err);
		if (M == NULL) {
		    if (genr->err == E_UNKVAR) {
			undefined_var_message(s);
		    }
		} else {
		    atom_set_matrix(atom, M, s);
		}
	    } else {
		get_var_from_matrix(s, genr, atom);
		if (genr->err) {
		    DPRINTF(("dead end at get_obs_value, s='%s'\n", s));
		}
	    }
	}
    } else if (numeric_string(s)) {
	/* plain number */
	atom->val = dot_atof(s);
	atom->atype = ATOM_SCALAR;
    } else if (*s == '"') {
	/* observation label? (basis for dummy series) */
	atom->val = get_observation_number(s, genr->pdinfo);
	if (atom->val > 0) {
	    atom->func = T_OBSNUM;
	} else{
	    genr->err = E_SYNTAX;
	}
    } else if (strchr(s, ':')) {
	/* time-series observation? (as in "obs > 1986:4") */
	atom->val = get_observation_number(s, genr->pdinfo);
	if (atom->val > 0) {
	    atom->atype = ATOM_SCALAR;
	} else {
	    genr->err = E_SYNTAX;
	}
    } else {
	DPRINTF(("dead end in token_make_atom, s='%s'\n", s));
	genr->err = E_SYNTAX;
    }

    if (genr->err) {
	free(atom);
	atom = NULL;
    }

    return atom;
}

static double get_lag_at_obs (int v, int tmp, int lag, 
			      const GENERATOR *genr, int t)
{
    int lt;
    double **Z;
    double x = NADBL;

    if (tmp) {
	Z = genr->tmpZ;
    } else {
	Z = *genr->pZ;
    }

    /* stacked X-section needs rather special handling */
    if (genr->pdinfo->structure == STACKED_CROSS_SECTION) {
	lt = t - lag * genr->pdinfo->pd;
	if (lt >= 0 && lt < genr->pdinfo->n) {
	    x = Z[v][lt];
	}
    } else if (dated_daily_data(genr->pdinfo)) {
	lt = t - lag;
	if (lt >= 0 && lt < genr->pdinfo->n) {
	    while (lt >= 0 && na(Z[v][lt])) {
		lt--;
	    }
	    x = Z[v][lt];
	}
    } else { /* the "standard" time-series case */
	lt = t - lag;
	if (lt >= 0 && lt < genr->pdinfo->n) {
	    x = Z[v][lt];
	}
    }

    /* post-process missing panel values */
    if (genr->pdinfo->structure == STACKED_TIME_SERIES) {
	char *p, obs[OBSLEN];
	int j;

	ntodate(obs, t, genr->pdinfo);
	p = strchr(obs, ':');
	j = atoi(p + 1);
	if (j <= lag) x = NADBL;
    }

    return x;
}

static double 
eval_atom (genatom *atom, GENERATOR *genr, int t, double a)
{
    double x = NADBL;

    if (atom->atype == ATOM_SCALAR) {
	/* constant, scalar variable, or specific obs from a series */
	if (atom->varnum >= 0 && atom->varnum < HIGHNUM) {
	    x = (*genr->pZ)[atom->varnum][atom->varobs];
	} else {
	    x = atom->val;
	}
	DPRINTF(("eval_atom: got scalar = %g\n", x));
    } else if (atom->tmpvar >= 0) {
	/* temporary variable */
	if (!atom->lag) {
	    x = genr->tmpZ[atom->tmpvar][t];
	    DPRINTF(("eval_atom: got temp obs = %g\n", x));
	} else {
	    x = get_lag_at_obs(atom->tmpvar, 1, atom->lag, genr, t);
	    DPRINTF(("eval_atom: got lagged temp obs (tmpvar %d, "
		    "lag %d) = %g\n", atom->tmpvar, atom->lag, x));
	}
    } else if (atom->varnum == TNUM) {
	/* trend/index variable */
	x = get_tnum(genr->pdinfo, t);
    } else if (atom->varnum >= 0) {
	/* regular variable or lagged term */
	if (!atom->lag) {
	    x = (*genr->pZ)[atom->varnum][t];
	    DPRINTF(("eval_atom: got data obs (var %d) = %g\n", 
		    atom->varnum, x));
	} else {
	    x = get_lag_at_obs(atom->varnum, 0, atom->lag, genr, t);
	    DPRINTF(("eval_atom: got lagged data obs (var %d, "
		    "lag %d) = %g\n", atom->varnum, atom->lag, x));
	}
    } else if (atom->func) {
	/* a function */
	if (STANDARD_MATH(atom->func)) {
	    x = evaluate_math_function(a, atom->func, &genr->err);
	    DPRINTF(("evaluated math func %d: %g -> %g\n", 
		     atom->func, a, x));
	} else if (MISSVAL_FUNC(atom->func)) {
	    x = evaluate_missval_func(a, atom->func);
	    DPRINTF(("evaluated missval func %d: %g -> %g\n", 
		     atom->func, a, x));
	} else if (atom->func == T_OBSNUM) {
	    x = (t + 1 == atom->val)? atom->val : 0.0;
	    DPRINTF(("evaluated obsnum at t=%d, returning %g\n",
		     t, x));
	} else if (atom->func == T_IDENTITY) {
	    DPRINTF(("identity func: passed along %g\n", a));
	    x = a;
	}
    } else if (atom->varnum < 0) {
	/* named list: not acceptable in this context */
	DPRINTF(("eval_atom: got named list: error\n"));
	if (atom->str != NULL) {
	    undefined_var_message(atom->str);
	}
	genr->err = E_UNKVAR;
    }

    return x;
}

static gretl_matrix *eval_matrix_atom (genatom *atom, GENERATOR *genr,
				       gretl_matrix *M)
{
    gretl_matrix *R = NULL;
    double x;

    if (atom->M != NULL) {
	R = atom->M;
    } else if (atom->func) {
	if (atom->func == T_DET) {
	    x = user_matrix_get_determinant(M, &genr->err);
	    R = gretl_matrix_from_scalar(x);
	} else if (atom->func == T_LDET) {
	    x = user_matrix_get_log_determinant(M, &genr->err);
	    R = gretl_matrix_from_scalar(x);
	} else if (atom->func == T_INV) {
	    R = user_matrix_get_inverse(M);
	} else if (atom->func == T_CHOL) {
	    R = user_matrix_cholesky_decomp(M);
	} else if (atom->func == T_SUMC) {
	    R = gretl_matrix_column_sum(M);
	} else if (atom->func == T_SUMR) {
	    R = gretl_matrix_row_sum(M);
	} else if (atom->func == T_CDMEAN) {
	    R = user_matrix_column_demean(M);
	} else if (atom->func == T_QR) {
	    R = user_matrix_QR_decomp(atom->str, genr->pZ, genr->pdinfo, 
				      genr->prn, &genr->err);
	} else if (atom->func == T_EIGSYM) {
	    R = user_matrix_eigen_analysis(atom->str, genr->pZ, genr->pdinfo, 
					   genr->prn, &genr->err, 1);
	} else if (atom->func == T_EIGGEN) {
	    R = user_matrix_eigen_analysis(atom->str, genr->pZ, genr->pdinfo, 
					   genr->prn, &genr->err, 0);
	} else if (atom->func == T_DIAG) {
	    R = gretl_matrix_get_diagonal(M, &genr->err);
	} else if (atom->func == T_TRANSP) {
	    R = gretl_matrix_copy_transpose(M);
	} else if (atom->func == T_SORT) {
	    R = user_matrix_get_sorted_vector(M, SORT_ASCENDING, &genr->err);
	} else if (atom->func == T_DSORT) {
	    R = user_matrix_get_sorted_vector(M, SORT_DESCENDING, &genr->err);
	} else if (atom->func >= T_NONE && atom->func < T_MATHMAX) {
	    R = user_matrix_get_transformation(M, atom->func);
	} else if (atom->func == T_IDENTITY) {
	    DPRINTF(("identity func: passing along M (%p) as R\n", (void *) M));
	    R = M;
	} else if (M != NULL) {
	    if (atom->func == T_ROWS) {
		x = gretl_matrix_rows(M);
		R = gretl_matrix_from_scalar(x);
	    } else if (atom->func == T_COLS) {
		x = gretl_matrix_cols(M);
		R = gretl_matrix_from_scalar(x);
	    } else {
		genr->err = E_TYPES;
		fprintf(stderr, "type error in eval_matrix_atom()\n");
	    }
	} else {
	    genr->err = E_TYPES;
	    fprintf(stderr, "type error in eval_matrix_atom()\n");
	}
    } else if (atom->atype == ATOM_SCALAR && !na(atom->val)) {
	R = gretl_matrix_from_scalar(atom->val);
    }

    if (R == NULL && !genr->err) {
	genr->err = 1;
    }

    return R;
}

/* Add a pre-allocated series to genr's temporary dataset: this
   function also takes charge of freeing x if the addition fails,
   and of setting atom's tmpvar number */

static int genr_add_temp_var (double *x, GENERATOR *genr, genatom *atom)
{
    int v = genr->tmpv;
    double **gZ;    

    gZ = realloc(genr->tmpZ, (v + 1) * sizeof *gZ);

    if (gZ == NULL) {
	free(x);
	genr->err = E_ALLOC;
    } else {
	atom->tmpvar = v;

	gZ[v] = x;
	genr->tmpZ = gZ;
	genr->tmpv += 1;
    }

    return genr->err;
}

static int add_random_series_to_genr (GENERATOR *genr, genatom *atom)
{
    double *x;

    x = get_random_series(genr->pdinfo, atom->func);
    if (x == NULL) return 1;

    return genr_add_temp_var(x, genr, atom);
}

static int add_model_series_to_genr (GENERATOR *genr, genatom *atom)
{
    char *s = atom->str;
    char *oname = NULL;
    char *key;
    double *x;
    int err = 0;

    if (*s == '$') {
	key = s;
    } else if ((key = strstr(s, ".$")) != NULL) {
	oname = gretl_strndup(s, key - s);
	key++;
    } else {
	return 1;
    }

    DPRINTF(("add_model_series_to_genr: oname=%s, key='%s'\n",
	     (oname)? oname : "none", key));

    x = saved_object_get_series(oname, key, genr->pdinfo, &err);
    DPRINTF((" x=%p, err=%d\n", (void *) x, err)); 

    free(oname);

    if (!err) {
	err = genr_add_temp_var(x, genr, atom);
    }

    return err;
}

static int add_mp_series_to_genr (GENERATOR *genr, genatom *atom)
{
    double *x;

    x = get_mp_series(atom->str, genr, atom->func, &genr->err);
    if (x == NULL) return 1;

    return genr_add_temp_var(x, genr, atom);
}

static double *eval_compound_arg (GENERATOR *genr,
				  genatom *this_atom)
{
    int t, t1, t2;
    genatom *atom;
    double *xtmp;

    /* changed 2005/11/28 */
#if 1
    t1 = 0; t2 = genr->pdinfo->n - 1;
#else
    t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;
#endif

    DPRINTF(("eval_compound arg, atom->func = %d\n", this_atom->func));

    if (peek_child_atom(this_atom) == NULL) {
	DPRINTF((" error: no child found\n"));
	genr->err = E_SYNTAX;
	return NULL;
    }

    xtmp = malloc(genr->pdinfo->n * sizeof *xtmp);
    if (xtmp == NULL) {
	genr->err = E_ALLOC;
	return NULL;
    }

    for (t=t1; t<=t2; t++) {
	double xbak = 0.0, x = 0.0;
	int level = 0;

	reset_atom_stack(genr);

	while ((atom = pop_child_atom(this_atom))) {
	    double y = eval_atom(atom, genr, t, x);
	    int err;

	    if (genr->err) {
		break;
	    }

	    if (atom->level < level) {
		x = calc_pop(genr);
	    }

	    x = calc_xy(x, y, atom->op, t, &err);

	    if (err) {
		genr_set_warn(genr);
	    }

	    if (atom->level > level) { 
		genr->err = calc_push(xbak, genr);
	    }

	    level = atom->level;
	    xbak = x;
	}

	if (genr->err) {
	    break;
	}

	reset_calc_stack(genr);
	xtmp[t] = x;
    }

    return xtmp;
}

static double *
get_target_fracdiff_series (GENERATOR *genr, genatom *atom,
			    double *param)
{
    char vname[VNAMELEN];
    double *x = NULL;
    int v, t;

    if (*atom->str == '\0') {
	genr->err = 1;
	return NULL;
    }

#if GENR_DEBUG
    fprintf(stderr, "fracdiff_series: atom->str = '%s'\n", atom->str);
#endif

    if (sscanf(atom->str, "%15[^,],%lf", vname, param) != 2) {
	genr->err = 1;
	return NULL;
    }

    if (fabs(*param) > 1.0) {
	genr->err = 1;
	return NULL;
    }

    v = varindex(genr->pdinfo, vname);
    if (v == 0 || v >= genr->pdinfo->v) {
	genr->err = 1;
	return NULL;
    }

    x = malloc(genr->pdinfo->n * sizeof *x);
    if (x != NULL) {
	for (t=0; t<genr->pdinfo->n; t++) {
	    x[t] = (*genr->pZ)[v][t];
	}
    }

#if GENR_DEBUG
    fprintf(stderr, "fracdiff_series: v = %d, frac = %g, x = %p\n", 
	    v, *param, (void *) x);
#endif

    return x;
}

static int add_tmp_series_to_genr (GENERATOR *genr, genatom *atom)
{
    double param = 0.0;
    double *x = NULL;
    double *y = NULL;

    atom_stack_set_parentage(genr);

    if (atom->func == T_FRACDIFF) {
	x = get_target_fracdiff_series(genr, atom, &param);
    } else {
	/* evaluate possibly compound arg */
	x = eval_compound_arg(genr, atom);
    }

    if (x == NULL) {
	return genr->err;
    }

    y = get_tmp_series(x, genr, atom->func, param);

    free(x);

    if (y == NULL) {
	return genr->err;
    }

    genr_add_temp_var(y, genr, atom);

    if (!genr->err) {
	atom_eat_children(atom);
    }

    return genr->err;
}

static double genr_get_child_status (GENERATOR *genr, genatom *atom)
{
    double val = NADBL;
    genatom *child = pop_child_atom(atom);

    if (child != NULL) {
	if (atom->func == T_VARNUM || atom->func == T_SERIES) {
	    if (child->varnum >= 0 && child->varnum < genr->pdinfo->v) {
		val = (atom->func == T_VARNUM)? child->varnum :
		    genr->pdinfo->vector[child->varnum];
	    } 
	} else if (atom->func == T_ISLIST || atom->func == T_NELEM) {
	    int *list = get_list_by_name(child->str);

	    if (atom->func == T_ISLIST) {
		val = (list != NULL);
	    } else {
		val = (list == NULL)? NADBL : list[0];
	    }
	} 
    }

    if (pop_child_atom(atom)) {
	/* should have only one child */
	genr->err = E_SYNTAX;
    }

    DPRINTF(("genr_get_object_status: child=%p, got %g\n", 
	     (void *) child, val));

    return val;
}

static int add_statistic_to_genr (GENERATOR *genr, genatom *atom)
{
    double val;

    atom_stack_set_parentage(genr);

    if (atom->func == T_VARNUM || atom->func == T_SERIES ||
	atom->func == T_ISLIST || atom->func == T_NELEM) {
	val = genr_get_child_status(genr, atom);
    } else {
	double *x = eval_compound_arg(genr, atom);

	if (x == NULL) {
	    return genr->err;
	}

	val = evaluate_statistic(x, genr, atom->func);
	free(x);
    }

    if (genr->err) {
	return genr->err;
    }

#if GENR_DEBUG
    fprintf(stderr, "add_statistic_to_genr:\n atom->func = %d (%s), val = %g, "
	    "now atom->scalar = 1\n", atom->func, get_genr_func_word(atom->func),
	    val);
#endif

    atom->val = val;
    atom->atype = ATOM_SCALAR;
    atom_eat_children(atom);
    atom->func = 0;

    return genr->err;
}

static int maybe_free_genr_matrix (gretl_matrix *M, const char *s)
{
    if (M == NULL || is_user_matrix(M) || matrix_is_on_atom(M)) {
	return 0;
    }

    MPRINTF(("tmp matrix '%s' at %p: not user or atom matrix, freeing\n", 
	     s, (void *) M));
    gretl_matrix_free(M);

    return 1;
}

static gretl_matrix *scalar_calc_AB (gretl_matrix *A, gretl_matrix *B,
				     char op, int *err)
{
    gretl_matrix *C = NULL;
    double x, y, z;

    x = gretl_vector_get(A, 0);
    y = gretl_vector_get(B, 0);
    z = calc_xy(x, y, op, 0, err);

    MPRINTF(("scalar_calc_AB: x=%g, y=%g, z=%g\n", x, y, z));

    if (!*err) {
	if (na(z)) {
	    *err = 1;
	} else {
	    C = gretl_matrix_from_scalar(z);
	    if (C == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    return C;
}

static int evaluate_matrix_genr (GENERATOR *genr)
{
    gretl_matrix *A = NULL;
    gretl_matrix *Abak = NULL;
    genatom *atom;
    int level = 0, npush = 0, npop = 0;
#if GENR_DEBUG
    int i = 0;
#endif

    reset_atom_stack(genr);

    MPRINTF(("\n*** starting evaluate_matrix_genr on\n"
	     "'%s'\n", genr->orig_s));

    while ((atom = pop_atom(genr)) && !genr->err) {
	gretl_matrix *Apop = NULL;
	int freeB = 1;

	gretl_matrix *B = eval_matrix_atom(atom, genr, A);
	if (B == A) {
	    freeB = 0;
	}

	DPRINTF(("\natom %d, atom->level %d\n", 
		 i++, atom->level));
	MPRINTF(("\nInitial A = %p; eval_matrix_atom gave B = %p\n", 
		 (void *) A, (void *) B)); 
#if GENR_DEBUG
	gretl_matrix_print(B, "eval_matrix_atom: got B =");
#endif

	if (genr->err) break;

	if (atom->level < level) { 
	    A = matrix_calc_pop(genr);
	    Apop = A;
	    npop++;
	}

	if (gretl_matrix_is_scalar(A) && gretl_matrix_is_scalar(B)) {
	    A = scalar_calc_AB(A, B, atom->op, &genr->err);
	} else {
	    A = matrix_calc_AB(A, B, atom->op, &genr->err);
	}

	MPRINTF(("after calc_AB, A = %p\n", (void *) A)); 
	if (A == B) {
	    freeB = 0;
	}

	if (genr->err) break;

	if (atom->level > level) {
	    int pad = atom->level - level - 1;

	    genr->err = matrix_calc_push(Abak, genr);
	    npush++;
	    while (pad--) {
		genr->err = matrix_calc_push(NULL, genr);
		DPRINTF(("pushed NULL at level %d\n", level));
		npush++;
	    }
	} else if (Abak != A && maybe_free_genr_matrix(Abak, "Abak")) {
	    /* not pushing matrix Abak: free it, if it's not "spoken for" */
	    Abak = NULL;
	}

	/* matrix 'B': free it if there's no further call on it */
	if (freeB && maybe_free_genr_matrix(B, "B")) {
	    B = NULL;
	}

	/* matrix 'Apop': free it if there's no further call on it */
	maybe_free_genr_matrix(Apop, "Apop");

	level = atom->level;

	MPRINTF(("old Abak = %p, setting Abak = %p\n", (void *) Abak, (void *) A));
	MPRINTF(("old Abak: user_matrix? %s; atom_matrix? %s\n",
		 is_user_matrix(Abak)? "yes" : "no",
		 matrix_is_on_atom(Abak)? "yes" : "no"));
	Abak = A;
    }

    if (!genr->err && npop > npush) {
	/* excess pushes are harmless? */
	fprintf(stderr, "genr error: npush = %d, npop = %d\n",
		npush, npop);
	genr->err = 1;
    }

    reset_matrix_calc_stack(genr);

    if (!genr->err) {
	MPRINTF(("done: pushing A = %p onto stack\n", (void *) A)); 
	genr->err = matrix_calc_push(A, genr);
    }

    return genr->err;
}

static int evaluate_genr (GENERATOR *genr)
{
    int t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;
    int m = 0, tstart = t1;
    int mseries = 0;
    int n_atoms = 0;
    genatom *atom;

    reset_atom_stack(genr);

    DPRINTF(("Doing evaluate_genr\n"));

    while ((atom = pop_atom(genr)) && !genr->err) {

	DPRINTF((" looking at atom %d (atype = %s)\n", n_atoms, 
		 ATOM_TYPE_STRING(atom->atype)));
	n_atoms++;

	if (mseries && n_atoms > 0) {
	    genr->err = E_TYPES;
	    fprintf(stderr, "type error in evaluate_genr()\n");
	    break;
	} 

	if (atom->varnum == genr->varnum && atom->lag > m) {
	    m = atom->lag;
	}

	if (atom->atype & ATOM_MODELDAT) {
	    genr->err = add_model_series_to_genr(genr, atom);
	} else if (atom->func == T_UNIFORM || atom->func == T_NORMAL) {
	    genr->err = add_random_series_to_genr(genr, atom);
	} else if (atom->func == T_DIFF || atom->func == T_LDIFF ||
		   atom->func == T_SDIFF || atom->func == T_CUM ||
		   atom->func == T_SORT || atom->func == T_DSORT ||
		   atom->func == T_RESAMPLE || atom->func == T_HPFILT ||
		   atom->func == T_BKFILT || atom->func == T_FRACDIFF) {
	    atom_stack_bookmark(genr);
	    genr->err = add_tmp_series_to_genr(genr, atom);
	    atom_stack_resume(genr);
	} else if (UNIVARIATE_STAT(atom->func)) {
	    atom_stack_bookmark(genr);
	    genr->err = add_statistic_to_genr(genr, atom);
	    atom_stack_resume(genr);
	} else if (MP_MATH(atom->func)) {
	    genr->err = add_mp_series_to_genr(genr, atom);
	} else if (needs_arg(atom->func) && !arg_atom_available(atom)) {
	    genr->err = E_SYNTAX; /* FIXME? */
	    fprintf(stderr, "%s(): needs an argument\n", 
		    get_genr_func_word(atom->func));
	} else if (atom_is_series_from_matrix(atom)) {
	    if (n_atoms > 1) {
		genr->err = E_TYPES;
	    } else {
		mseries = 1;
	    }
	}
    }

    DPRINTF(("evaluate_genr: n_atoms = %d\n", n_atoms));

    if (genr->err) {
	return genr->err;
    }

    if (n_atoms == 2) {
	reset_atom_stack(genr);
	atom = pop_atom(genr);
	if (atom->func == T_SORT || atom->func == T_DSORT) {
	    genr_set_simple_sort(genr);
	}
    }

    if (atom_stack_check_for_scalar(genr)) {
	genr_set_scalar(genr);
    } else {
	genr_unset_scalar(genr);
    }

    DPRINTF(("evaluate_genr: check for scalar, result = %d\n", 
	     genr_is_scalar(genr)));

    if (genr_is_scalar(genr)) {
	t2 = tstart;
    } else if (tstart < m) {
	/* autoregressive genr */
	tstart = m;
	for (t=t1; t<tstart; t++) {
	    /* not quite sure here: this replicates "sim" behavior */
	    genr->xvec[t] = (*genr->pZ)[genr->varnum][t];
	}
    }

    for (t=tstart; t<=t2; t++) {
	double xbak = 0.0, x = 0.0;
	int level = 0, npush = 0, npop = 0;

	reset_atom_stack(genr);

	while ((atom = pop_atom(genr))) {
	    double y = eval_atom(atom, genr, t, x);
	    int err;

	    if (genr->err) {
		break;
	    }

	    if (atom->level < level) { 
		x = calc_pop(genr);
		npop++;
		DPRINTF(("popped %g\n", x));
	    }

	    x = calc_xy(x, y, atom->op, t, &err);

	    if (err) {
		genr_set_warn(genr);
	    }

	    if (atom->level > level) {
		int pad = atom->level - level - 1;

		genr->err = calc_push(xbak, genr);
		npush++;
		DPRINTF(("pushed %g at level %d\n", xbak, level));
		while (pad--) {
		    genr->err = calc_push(0, genr);
		    DPRINTF(("pushed 0 at level %d\n", level));
		    npush++;
		}
	    }

	    level = atom->level;
	    xbak = x;
	}

	if (!genr->err && npop > npush) {
	    /* excess pushes are harmless? */
	    fprintf(stderr, "genr error: npush = %d, npop = %d\n",
		    npush, npop);
	    genr->err = 1;
	}

	reset_calc_stack(genr);

	if (genr->err) {
	    break;
	}

	genr->xvec[t] = x;
	if (m > 0 && !na(x)) {
	    /* autoregressive genr */
	    (*genr->pZ)[genr->varnum][t] = x;
	}
    }

    return genr->err;
}

/* insert right paren at the end of the numeric portion of
   the given string.
*/

static int insert_right_paren (char *s)
{
    char *rt;
    char *ins;

    errno = 0;

    strtod(s, &rt);
    
    if (errno == ERANGE) {
	/* bad fp value: give up, and pass the buck */
	return 0;
    }

    if (!strcmp(s, rt)) {
	/* no fp conversion: crawl forward to end of atom? */
	while (*s && !op_level(*s)) s++;
	ins = s;
    }  else if (*rt != '\0') {
	/* partial fp conversion: put paren at end of numeric part */
	ins = rt;
    } else {
	/* full fp conversion */
	ins = s + strlen(s);
    }

    s = ins;

    /* gotcha: next is '^' operator, with precedence > unary op? */
    if (*(s + 1) == '^' || (*(s + 1) == '*' && *(s + 2) == '*')) {
	while (*s && !op_level(*s)) s++;
    }

    if (*s == '\0') {
	*s = ')';
	*(s + 1) = '\0';
    } else {
	char *p = s;

	memmove(s + 1, s, strlen(s) + 1);
	*p = ')';
    }
	
    return 0;
}

/* insert a "ghost" zero before a unary operator, and if need
   be insert a pair of parentheses around the value
*/

static int insert_ghost_zero (char *start, char *p, int *np)
{
    int n = strlen(start);
    int pn = strlen(p);
    char *ins = p;
    int err = 0;

    if (p - start > 0 && *(p - 1) != '(') {
	*np = 1;
    } else {
	*np = 0;
    }

    /* do we have space to make the insertion? */
    if (n + (2 * (*np)) >= MAXLINE) {
	return 1;
    }

    /* move material right */
    memmove(p + 1 + *np, p, pn + 1);

    p = ins;

    if (*np) {
	*p = '(';
	*(p + 1) = '0';
	err = insert_right_paren(p + 2);
    } else {
	*p = '0';
    }

    return err;
}

/* check if there's a genr function name to the left of point:
   have to take care in case there's a variable name which 
   includes a substring equal to a function word.
*/

#define VNCHAR(c) (isalnum(c) || c == '_')

static int function_word_to_left (const char *s, int n)
{
    int i, ret = 0;

    DPRINTF(("function_word_to_left: s='%s', n=%d\n", s, n));

    for (i=1; i<=n; i++) {
	const char *p = s + n - i;

	if (!isalpha(*p)) {
	    break;
	} else if ((i == n || !VNCHAR(*(p-1))) && 
		   genr_function_from_string(p)) {
	    ret = 1;
	    break;
	}
    }

    return ret;
}

/* Given: there is a '-' or '+' at *p.  N.B. This needs care, because
   plus/minus may be part of a lag variable specification.
*/

static int unary_op_context (char *start, char *p)
{
    int pos = p - start;
    int ret = 0;

    DPRINTF(("unary_op_context: start='%s', *p='%s', pos=%d\n", 
	     start, p, pos)); 

    /* plus/minus at very start of formula, or directly preceded
       by another operator: must be plain unary */
    if (pos == 0 || op_level(*(p-1))) {
	ret = 1;
    }

    /* plus/minus preceded _only_ by left paren: unary */
    else if (pos == 1 && *start == '(') {
	ret = 1;
    }

    /* plus/minus preceded by left paren, preceded by operator
       or function: again unary */
    else if (pos >= 2 && *(p-1) == '(') {
	if (op_level(*(p-2))) {
	    ret = 1;
	} else if (function_word_to_left(start, pos - 1)) {
	    ret = 1;
	}
    }

    DPRINTF(("unary_op_context: ret=%d\n", ret));

    return ret;
}

static int catch_special_operators (GENERATOR *genr, char *s)
{
    char *p = s;
    int err = 0;
    int lshift;

    while (*s) {
	lshift = 0;

	if (*s == '!' && *(s+1) == '=') {
	    *s = OP_NEQ;
	    lshift = 1;
	} else if (*s == '>' && *(s+1) == '=') {
	    *s = OP_GTE;
	    lshift = 1;
	} else if (*s == '<' && *(s+1) == '=') {
	    *s = OP_LTE;
	    lshift = 1;
	} else if (*s == '.' && *(s+1) == '*') {
	    *s = OP_DOTMULT;
	    lshift = 1;
	} else if (*s == '.' && *(s+1) == '/') {
	    *s = OP_DOTDIV;
	    lshift = 1;
	} else if (*s == '.' && *(s+1) == '^') {
	    *s = OP_DOTPOW;
	    lshift = 1;
	} else if (*s == '*' && *(s+1) == '*') {
	    if (genr_is_matrix(genr)) {
		*s = OP_KRON;
	    } else {
		*s = '^';
	    }
	    lshift = 1;
	} else if ((*s == '-' || *s == '+') && unary_op_context(p, s)) {
	    int np; /* need to insert parentheses? */

	    err = insert_ghost_zero(p, s, &np);
	    s += 1 + np;
	} else if (*s == '&' && *(s+1) == '&') {
	    lshift = 1;
	} else if (*s == '|' && *(s+1) == '|') {
	    lshift = 1;
	} else if (*s == '=' && *(s+1) == '=') {
	    lshift = 1;
	}

	if (lshift) {
	    memmove(s + 1, s + 2, 1 + strlen(s + 2));
	}

	s++;
    }

    return err;
}

static int op_level (int c)
{
    if (c == '^' || c == '!' || c == OP_DOTPOW) 
	return 1;
    if (c == '*' || c == '/' || c == '%') 
	return 2;
    if (c == OP_DOTMULT || c == OP_DOTDIV || c == OP_KRON) 
	return 2;
    if (c == '+' || c == '-' || c == '~') 
	return 3;
    if (c == '>' || c == '<' || c == OP_GTE || c == OP_LTE) 
	return 4;
    if (c == '=' || c == OP_NEQ) 
	return 5;
    if (c == '&') 
	return 6;
    if (c == '|') 
	return 7;

    return 0;
}

static int string_arg_function_word (const char *s, GENERATOR *genr)
{
    if (!strncmp(s, "vcv", 3) ||
	!strncmp(s, "corr", 4) ||
	!strncmp(s, "cov", 3) ||
	!strncmp(s, "pvalue", 6) ||
	!strncmp(s, "critical", 8) ||
	!strncmp(s, "fracdiff", 8) ||
	!strncmp(s, "mpow", 4) ||
	!strncmp(s, "mlog", 4) ||
	!strncmp(s, "zeros", 5) ||
	!strncmp(s, "ones", 4) ||
	!strncmp(s, "I", 1)) {
	return 1;
    }

    if (genr_is_matrix(genr)) {
	if (!strncmp(s, "uniform", 7) ||
	    !strncmp(s, "normal", 6) ||
	    !strncmp(s, "qrdecomp", 8) ||
	    !strncmp(s, "eigensym", 8) ||
	    !strncmp(s, "eigengen", 8)) {
	    return 1;
	}
    }

    return 0;
}

static int matrix_scalar_function_word (const char *s)
{
    if (!strncmp(s, "det", 3) ||
	!strncmp(s, "ldet", 4) ||
	!strncmp(s, "rows", 4) ||
	!strncmp(s, "cols", 4) ||
	!strncmp(s, "onenorm", 7) ||
	!strncmp(s, "rcond", 5) ||
	!strncmp(s, "tr", 2)) {
	return 1;
    }

    return 0;
}

static int is_varname_or_model_data_selector (const char *s)
{
    int ret = 0;

    if (*s != '$') {
	s += gretl_varchar_spn(s);
	if (*s == '\0') {
	    /* got a possible varname */
	    return 1;
	}
	if (!strncmp(s, ".$", 2)) {
	    s++;
	} else if (*s == '[') {
	    /* matrix slice? */
	    if (s[strlen(s) - 1] == ']') {
		return 1;
	    } else {
		return 0;
	    }
	} else {
	    return 0;
	}
    }

    if (!strncmp(s, "$coeff(", 7) ||
	!strncmp(s, "$stderr(", 8) ||
	!strncmp(s, "$vcv(", 5) ||
	!strncmp(s, "$rho(", 5)) {
	int rp = 0;

	s = strchr(s, '(') + 1;
	while (*s) {
	    if (*s == ')') rp++;
	    s++;
	}
	if (rp == 1 && *(s-1) == ')') {
	    ret = 1;
	} 
    }

    return ret;
}

static int token_is_atomic (const char *s, GENERATOR *genr)
{
    int count = 0;
    int atomic = 0;

    DPRINTF(("token_is_atomic: looking at '%s'\n", s));

    if (numeric_string(s)) {
	/* number in scientific notation */
	atomic = 1;
    } else if (isdigit((unsigned char) *s)) {
	/* not numeric string but begins with digit: can't be atomic */
	atomic = 0;
    } else if (get_lagvar(s, NULL, genr)) {
	/* treat lag variable as atom */
	atomic = 1;
    } else if (is_varname_or_model_data_selector(s)) {
	/* treat, e.g., "$coeff(x)" as atom */
	atomic = 1;
    } else {
	while (*s) {
	    /* token is non-atomic if it contains an operator
	       or parentheses, otherwise atomic */
	    if (op_level(*s) || *s == '(') {
		count++;
	    }
	    if (count > 0) {
		break;
	    }
	    s++;
	}
	if (count == 0) {
	    atomic = 1;
	} 
    }

    DPRINTF(("token_is_atomic: returning %d\n", atomic));

    return atomic;
}

static int token_is_function (char *s, GENERATOR *genr, int level)
{
    int wlen = 0;
    int ret;
    const char *p = s;

    DPRINTF(("token_is_function: looking at '%s'\n", s));

    if (*s == '$') {
	/* can look like a function, but is not */
	return 0;
    }

    while (*p) {
	if (!op_level(*p) && *p != '(') {
	    wlen++; /* might be a problem */
	} else {
	    break;
	}
	p++;
    }

    DPRINTF((" after processing, looking at '%s'\n", p));

    ret = (*p == '(' && p[strlen(p) - 1] == ')'); 

    if (ret) {
	DPRINTF((" yes, token is function...\n"));
	/* FIXME compound matrix arg? */
	if (string_arg_function_word(s, genr) || 
	    (matrix_scalar_function_word(s) && !genr_is_matrix(genr))) {
	    return ret;
	} else {
	    char subtok[TOKLEN];

	    strcpy(subtok, strchr(s, '(') + 1);
	    subtok[strlen(subtok) - 1] = '\0';

	    if (wlen == 0) {
		strcpy(s, "ident(#)");
	    } else {
		strcpy(strchr(s, '(') + 1, "#)");
	    }

	    if (*subtok) {
		math_tokenize(subtok, genr, ++level);
	    }
	}
    }

    return ret;
}

static int contains_no_operator (const char *s)
{
    while (*s++) {
	if (op_level(*s)) {
	    return 0;
	}
    }

    return 1;
}

#define TOKENIZE_LEVEL_MAX 256	

static int stack_op_and_token (char *s, GENERATOR *genr, int level)
{
    char tok[TOKLEN];
    int wrapped = 0;
    char op = 0, last = 0;
    size_t n = strlen(s);

    *tok = '\0';

    DPRINTF(("stack_op_and_token, level %d: looking at '%s'\n", level, s));

    if (level > TOKENIZE_LEVEL_MAX) {
	genr->err = E_PARSE;
	return genr->err;
    }

    if (n > 0) last = s[n-1];

    if (op_level(*s)) {
	/* leading character is an operator */
	op = *s;
	s++;
    } 

#if GENR_DEBUG
    if (isprint(op)) fprintf(stderr, "op = '%c'\n", op);
    else fprintf(stderr, "op = %d\n", op);
#endif

    if (*s == '(') {
	if (last == ')') {
	    wrapped = 1;
	    DPRINTF(("term is wrapped in parens\n"));
	    if (contains_no_operator(s)) {
		/* unwrap the term */
		s++;
		s[strlen(s)-1] = '\0';
	    }
	} else {
	    genr->err = E_UNBAL;
	}
    }

    if (!genr->err) {
	strcpy(tok, s);
	if (!token_is_atomic(tok, genr) &&
	    !token_is_function(tok, genr, level)) {
	    DPRINTF(("token is not atomic...\n"));
	    genr->err = math_tokenize(tok, genr, ++level);
	} else {
	    genatom *atom;

#if GENR_DEBUG
	    if (*tok == '\0') fprintf(stderr, "genr: got a blank token\n");
	    else fprintf(stderr, "token = '%s'\n", tok);
#endif
	    atom = token_make_atom(tok, op, genr, level);
	    if (atom != NULL) {
		genr->err = push_atom(atom);
#if GENR_DEBUG
		if (!genr->err) fprintf(stderr, "pushed atom OK\n");
#endif
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
	    if (s[i] == '(') {
		pcount++;
	    } else if (s[i] == ')') {
		pcount--;
	    }
	    if (pcount == 0) { 
		/* opening paren matched before end */
		strip = 0;
		break;
	    }
	}
	if (strip) {
	    s[n-1] = '\0';
	}
    }    

    return strip;
}

static int math_tokenize (char *s, GENERATOR *genr, int level)
{
    static char prev[TOKLEN];
    static int oldlevel;
    char tok[TOKLEN];
    const char *q, *p;
    int inparen = 0;
    int strip = strip_wrapper_parens(s);

    if (strip) {
	DPRINTF(("math_tokenize: stripped wrapper parens\n"));
	s++;
    }

    DPRINTF(("math_tokenize, level %d: looking at '%s'\n", level, s));

    q = p = s;

    if (level > TOKENIZE_LEVEL_MAX) {
	genr->err = E_PARSE;
	return genr->err;
    }

    if (token_is_atomic(s, genr)) {
	DPRINTF(("math_tokenize: string is atomic token\n"));
	goto atomic_case; 
    } else if (!strcmp(prev, s) && level == oldlevel + 1) {
	DPRINTF(("math_tokenize: going round in circles\n"));
#if 0
	/* This needs thought, work: generating an error here causes
	   failure in some instances where, if we let computation
	   proceed, we get OK results!  Specifically, this breaks the
	   computation of derivatives in some of the NIST nls test
	   cases (big, ugly formulae, e.g. Hahn1).  AC 2005/11/04.
	*/
	genr->err = E_PARSE;
	return genr->err;
#endif
    }

    *prev = 0;
    strncat(prev, s, TOKLEN - 1);

    while (*p && !genr->err) {
	DPRINTF(("math_tokenize: inner loop '%s'\n", p));

	/* 2006/3/16: added '[' and ']' here */
	if (*p == '(' || *p == '[') inparen++;
	else if (*p == ')' || *p == ']') inparen--;

	if (inparen < 0) {
	    DPRINTF(("error: inparen < 0: '%s'\n", s));
	    return E_UNBAL;
	}

	if (!inparen && op_level(*p)) {
	    if (p - q > TOKLEN - 1) {
		fprintf(stderr, "genr error: token too long: '%s'\n", q);
		return 1;
	    }
	    /* found a break point: ship out the left-hand term */
	    if (p - q > 0) {
		*tok = '\0';
		strncat(tok, q, p - q);
		oldlevel = level;
		stack_op_and_token(tok, genr, level);
		q = p;
	    }
	    /* ... and peek at the right-hand term */
	    if (numeric_string(p)) {
		DPRINTF(("got a constant numeric token\n"));
		goto atomic_case;
	    }
	}
	p++;
    }

    if (inparen != 0) {
	fprintf(stderr, "error: inparen != 0: '%s'\n", s);
	return E_UNBAL;
    }    

 atomic_case:
	
    if (*q && !genr->err) {
	if (strlen(q) > TOKLEN - 1) {
	    fprintf(stderr, "genr error: remainder too long: '%s'\n", q);
	    return 1;
	}
	strcpy(tok, q);
	oldlevel = level;
	genr->err = stack_op_and_token(tok, genr, level);
    }

    return genr->err;
}

static int count_ops (char *s, int *opcount)
{
    int level, maxlev = 0;

    while (*s++) {
	level = op_level(*s);
	opcount[level] += 1;
	if (level > maxlev) {
	    maxlev = level;
	}
    }

    return maxlev;
}

int insert_paren (char *s, int pos, char lr)
{
    static int lpos;
    int i, rpos, n = strlen(s);

    if (n + 1 >= MAXLINE) {
	return 1;
    }

    /* move material right */
    for (i=n+1; i>=pos+1; i--) {
	s[i] = s[i - 1];
    }

    if (lr == 'L') {
	lpos = pos + 1;
	s[lpos] = '(';
    } else {
	rpos = pos + 1;
	s[rpos] = ')';
    }

    return 0;
}

static int paren_state (char c, int *state, char lr)
{
    int s = *state;

    /* modification 2006/3/16: make '[' and ']' count as
       parens: this may cause trouble, though it solves
       an immediate problem with matrix selections
    */

    if (c == '(' || c == '[') {
	if (lr == 'L') {
	    if (s > 0) {
		s--;
	    }
	} else {
	    s++;
	}
    } else if (c == ')' || c == ']') {
	if (lr == 'R') {
	    if (s > 0) {
		s--;
	    }
	} else {
	    s++;
	}
    }

    *state = s;

    return s;
}

static int parenthesize (char *str)
{
    int i, k, oppos, n;
    int level1 = 0, level2;  
    int priority, start, lpins, inparens;
    int rpar, pbak, maxlevel;
    int opcount[LEVELS + 1] = {0};

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
	for (i=oppos; i>0; i--) { /* was ; i>-0; */
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
		if (insert_paren(str, (i == n - 1 && str[i] != ')')? 
				 i: i - 1, 'R')) {
		    return 1;
		}
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

const char *get_genr_func_word (int fnum)
{
    int i;

    for (i=0; funcs[i].fnum != 0; i++) {
	if (fnum == funcs[i].fnum) {
	    return funcs[i].fword;
	}
    }

    return NULL;
}

int genr_function_from_string (const char *s)
{
    char word[VNAMELEN];
    const char *p;
    int i;

    *word = 0;

    p = strchr(s, '(');
    if (p != NULL && p - s <= 8) {
	strncat(word, s, p - s);
    } else {
	strncat(word, s, 8);
    }

    if (word[0] == '\'' && word[1] == 0) {
	return T_TRANSP;
    }

    for (i=0; funcs[i].fnum != 0; i++) {
	if (!strcmp(word, funcs[i].fword)) {
	    return funcs[i].fnum;
	}
    }

    /* aliases */
    if (!strcmp(word, "ln")) {
	return T_LOG;
    }

    return 0;
}

/**
 * gretl_reserved_word:
 * @str: string to be tested.
 *
 * Returns non-zero if @str is a reserved word that cannot figure as the
 * name of a user-defined variable, otherwise 0.
 */

int gretl_reserved_word (const char *str)
{
    const char *uses[] = {
	N_("constant"),
	N_("plotting variable"),
	N_("internal variable"),
	N_("math function")
    };
    int ret = 0;

    if (!strcmp(str, "const") ||
	!strcmp(str, "CONST") ||
	!strcmp(str, "pi") ||
	!strcmp(str, "NA") ||
	!strcmp(str, "null")) {
	ret = 1;
    } else if (!strcmp(str, "t") ||
	       !strcmp(str, "annual") ||
	       !strcmp(str, "qtrs") ||
	       !strcmp(str, "months") ||
	       !strcmp(str, "hours")) {
	ret = 2;
    } else if (!strcmp(str, "obs") ||
	       !strcmp(str, "series")) {
	ret = 3;
    } else if (genr_function_from_string(str)) {
	ret = 4;
    }

    if (ret > 0) {
	sprintf(gretl_errmsg, _("'%s' refers to a %s and may not be used as a "
			    "variable name"), str, _(uses[ret - 1])); 
    }
 
    return ret;
}

/* allow stuff like "genr foo += 3.0", as abbreviation for
   "genr foo = foo + 3.0"
*/

static int 
expand_operator_abbrev (char *s, const char *lhs, char op)
{
    int llen = strlen(lhs);
    int add_parens = 0;
    int i, err = 0;

    if (*s != '(') {
	add_parens = 1;
    }

    /* do we have space to make the insertion? */
    if (strlen(s) + llen + 2 + 2 * add_parens >= MAXLINE) {
	err = 1;
    } else {
	memmove(s + llen + 1 + add_parens, s, strlen(s) + 1);
	for (i=0; i<llen; i++) {
	    s[i] = lhs[i];
	}
	s[i] = op;
	if (add_parens) {
	    s[i+1] = '(';
	    strcat(s, ")");
	}
    }

    return err;
}

static void excise_bracketed_term (char *s)
{
    char *p, *q;

    if ((p = strchr(s, '[')) && (q = strchr(p, ']'))) {
	memmove(p, q + 1, strlen(q));
    }
}

static int split_genr_formula (char *s, GENERATOR *genr)
{
    char *p;
    int len, err = 0;

    if (genr->obs >= 0 || *genr->label != '\0') {
	/* obs number or matrix slice on left */
	excise_bracketed_term(s);
    }    

    *genr->lhs = '\0';

    if ((p = strchr(s, '=')) != NULL) {
	char op = 0;

	*p = '\0';

	if (*(p + 1) == '\0') {
	    /* got "equals", but no rhs */
	    err = E_SYNTAX;
	} else {
	    len = strlen(s);
	    if (len > 1) {
		int l = op_level(s[len - 1]);
		
		if (l == 2 || l == 3) {
		    op = s[len - 1];
		    s[len - 1] = '\0';
		}
	    }

	    /* should we warn if lhs name is truncated? */
	    strncat(genr->lhs, s, VNAMELEN - 1);

	    if (gretl_reserved_word(genr->lhs)) {
		err = 1;
	    } else {
		p++;
		memmove(s, p, strlen(p) + 1);
	    }
	}

	if (!err && op != 0) {
	    err = expand_operator_abbrev(s, genr->lhs, op);
	}
    } else {
	/* no '=' sign: maybe "foo++" or "baz--"? */
	len = strlen(s);
	if (len < VNAMELEN + 1 && gretl_varchar_spn(s) == len - 2) {
	    if (strstr(s, "++")) {
		strncat(genr->lhs, s, len - 2);
		sprintf(s, "%s+1", genr->lhs);
	    } else if (strstr(s, "--")) {
		strncat(genr->lhs, s, len - 2);
		sprintf(s, "%s-1", genr->lhs);
	    }
	}
    }

    return err;
}

/* for daily dates in brackets: replace '/' with ':', and
   while we're at it, check for unbalanced '['
*/

static int fix_obs_in_brackets (char *s)
{
    int err = 0;

    if (s != NULL && *s != 0) {
	while ((s = strchr(s, '['))) {
	    s++;
	    if (strchr(s, ']') == NULL) {
		err = E_SYNTAX;
		break;
	    }
	    while (*s != ']') {
		if (*s == '/') *s = ':';
		s++;
	    }
	}
    }

    return err;
}

/* copy src to targ (up to len), eliminating all spaces */

static void copy_compress (char *targ, const char *src, int len)
{
    int j = 0;

    while (*src && j < len) {
	if (*src != ' ') {
	    targ[j++] = *src;
	}
	src++;
    }

    targ[j] = '\0';
}

static void genr_set_submatrix_string (GENERATOR *genr, const char *s)
{
    sprintf(genr->label, "[%s]", s);
}

static void get_genr_formula (char *formula, const char *line,
			      GENERATOR *genr)
{
    char name[32], tag[32];
    char first[VNAMELEN];

    if (string_is_blank(line)) {
	return;
    }

    /* look at first 'word' in line */
    sscanf(line, "%15s", first);

    if (!strcmp(first, "genr")) {
	line += 4;
    } else if (!strcmp(first, "eval")) {
	genr_unset_save(genr);
	line += 4;
    } else if (!strcmp(first, "series")) {
	genr_force_series(genr);
	line += 6;
    } else if (!strcmp(first, "scalar")) {
	genr_set_require_scalar(genr);
	line += 6;
    } 

    while (isspace((unsigned char) *line)) {
	line++;
    } 

    /* allow for generating a single value in a series, or
       for generating a submatrix */
    if (sscanf(line, "%31[^[ =][%31[^]]", name, tag) == 2) {
	if (genr_is_matrix(genr)) {
	    genr_set_submatrix_string(genr, tag);
	} else {
	    genr->obs = get_t_from_obs_string(tag, (const double **) *genr->pZ, 
					      genr->pdinfo);
	    if (genr->obs < 0) {
		/* observation string not found: implicit matrix? */
		if (get_matrix_by_name(name, genr->pdinfo) != NULL) {
		    genr_set_submatrix_string(genr, tag);
		} else {
		    genr->err = 1;
		    return;
		}
	    }
	}
    }

    *formula = '\0';

    copy_compress(formula, line, MAXLINE - 10);
}

static int gentoler (const char *s, PRN *prn)
{
    int err = 0;

    if (numeric_string(s)) {
	double x = dot_atof(s);

	err = set_nls_toler(x);
	if (!err) {
	    pprintf(prn, _("Set tolerance to %g"), x);
	    pputc(prn, '\n');
	}
    } else {
	strcpy(gretl_errmsg, _("The setting for \"toler\" must be numeric"));
	err = 1;
    }

    return err;
}

static void 
make_genr_varname (GENERATOR *genr, const char *vname)
{
    if (!strncmp(vname, "__", 2)) {
	strcpy(genr->varname, vname + 2);
    } else {
	strcpy(genr->varname, vname);
    }
}

/* create the label that will form the description of
   the generated variable */

static void write_genr_label (GENERATOR *genr, int oldv)
{
    if (*genr->label == '\0') {
	int mc, len = 0;

	if (genr->varnum < oldv) {
	    mc = get_model_count();
	    if (mc > 0) {
		sprintf(genr->label, _("Replaced after model %d: "), mc);
		len = strlen(genr->label);
	    }
	}	

	if (strlen(genr->orig_s) > MAXLABEL - 1 - len) {
	    strncat(genr->label, genr->orig_s, MAXLABEL - 4 - len);
	    strcat(genr->label, "...");
	} else {
	    strncat(genr->label, genr->orig_s, MAXLABEL - 1);
	}
    }

    strcpy(VARLABEL(genr->pdinfo, genr->varnum), genr->label);
}

/* allocate temporary series used as workspace */

static int genr_add_xvec (GENERATOR *genr)
{
    int i, n = genr->pdinfo->n;

    genr->xvec = malloc(n * sizeof *genr->xvec);

    if (genr->xvec == NULL) {
	genr->err = E_ALLOC;
    } else {
	for (i=0; i<n; i++) {
	    genr->xvec[i] = 0.0;
	}
    }

    return genr->err;
}

#define genr_special_func(s) (strcmp(s, "dummy") == 0 || \
                              strcmp(s, "paneldum") == 0 || \
                              strcmp(s, "unitdum") == 0 || \
                              strcmp(s, "time") == 0 || \
                              strcmp(s, "index") == 0 || \
                              strcmp(s, "unit") == 0 || \
                              strncmp(s, "toler=", 6) == 0)

/* special uses of genr which are not of the form "lhs = rhs" */

static int genr_handle_special (const char *s, GENERATOR *genr, 
				double ***pZ, DATAINFO *pdinfo)
{
    int orig_v = pdinfo->v;
    int do_message = 0;
    int err = 0;

    if (!strcmp(s, "dummy")) {
	int di0 = dummy(pZ, pdinfo, 0);

	if (di0 == 0) {
	    err = 1;
	} else if (di0 == orig_v) {
	    pputs(genr->prn, _("Periodic dummy variables generated.\n"));
	} else {
	    pputs(genr->prn, _("Periodic dummy variables already present.\n"));
	}
    } else if (!strcmp(s, "paneldum")) {
	err = paneldum(pZ, pdinfo);
	if (!err) {
	    pputs(genr->prn, _("Panel dummy variables generated.\n"));
	}
    } else if (strcmp(s, "unitdum") == 0) {
	err = panel_unit_dummies(pZ, pdinfo);
	if (!err) {
	    pputs(genr->prn, _("Panel dummy variables generated.\n"));
	}
    } else if (!strncmp(s, "toler=", 6)) {
	err = gentoler(s + 6, genr->prn);
    } else if (!strcmp(s, "time")) {
	err = genrtime(pZ, pdinfo, 1);
	do_message = 1;
    } else if (!strcmp(s, "index")) {
	err = genrtime(pZ, pdinfo, 0);
	do_message = 1;
    } else if (!strcmp(s, "unit")) {
	err = genrunit(pZ, pdinfo);
	do_message = 1;
    } 

    if (!err && do_message) {
	strcpy(genr->varname, s);
	genr->varnum = varindex(pdinfo, s);
	genr_unset_scalar(genr);
	compose_genr_msg(genr, orig_v);
    }	    

    return err;
}

static void maybe_set_matrix_genr (GENERATOR *genr)
{
    if (!genr_is_matrix(genr) && 
	get_matrix_by_name(genr->lhs, genr->pdinfo) != NULL) {
	genr->flags = GENR_SAVE | GENR_MATRIX;
	genr->err = genr_matrix_init(genr, OPT_M);
    }
}

static void genr_type_check (GENERATOR *genr)
{
    int v = genr->varnum;

    if (genr->pdinfo->vector[v]) {
	/* assigning to existing series variable */
	if (genr_require_scalar(genr)) {
	    sprintf(gretl_errmsg, _("'%s' is the name of a data series"), 
		    genr->pdinfo->varname[v]);
	    genr->err = E_DATA;
	}
    }
}

GENERATOR *
genr_compile (const char *line, double ***pZ, DATAINFO *pdinfo, gretlopt opt,
	      PRN *prn)
{
#if GENR_DEBUG
    genatom *atom;
#endif
    GENERATOR *genr;
    char s[MAXLINE];

    *gretl_errmsg = *s = '\0';

    genr = genr_new(pZ, pdinfo, opt, prn);
    if (genr == NULL) {
	return NULL;
    }

    /* grab the expression, skipping the command word 
       and compressing spaces */
    get_genr_formula(s, line, genr);

    if (*s == '\0') {
	DPRINTF(("\n*** get_genr_formula failed: line='%s'\n", line));
	genr->err = E_EQN;
	return genr;
    }

    /* record the full genr expression */
    strcpy(genr->orig_s, s);

    DPRINTF(("\n*** starting genr, s='%s'\n", s));

    /* special cases which are not of the form "lhs=rhs" */
    if (genr_special_func(s)) {
	genr->err = genr_handle_special(s, genr, pZ, pdinfo);
	genr->done = 1;
	return genr;
    }

    /* split into lhs = rhs */
    genr->err = split_genr_formula(s, genr);
    if (genr->err) {
	return genr;
    }

    DPRINTF(("after split, genr->lhs='%s', s='%s'\n", genr->lhs, s));

    if (*genr->lhs != '\0') {
	if (strncmp(genr->lhs, "$nl", 3) && 
	    strncmp(genr->lhs, "__", 2) && 
	    check_varname(genr->lhs)) {
	    genr->err = E_SYNTAX;
	}
	genr->varnum = varindex(pdinfo, genr->lhs);
	if (genr->varnum == pdinfo->v) {
	    maybe_set_matrix_genr(genr);
	} else {
	    genr_type_check(genr);
	}
    } else {
	/* no "lhs=" bit */
	if (!(genr_doing_save(genr))) {
	    strcpy(genr->lhs, "$eval");
	} else {
	    genr->err = E_SYNTAX;
	}
    }

    if (genr->err) {
	return genr;
    }    

    /* process any daily dates in brackets */
    genr->err = fix_obs_in_brackets(s);

    /* special case of generating a single observation */
    if (!genr->err && genr->obs >= 0) {
	if (genr->varnum >= pdinfo->v) {
	    genr->err = E_UNKVAR;
	} else if (!pdinfo->vector[genr->varnum]) {
	    genr->err = E_DATA;
	}
    }

    if (genr->err) {
	return genr;
    }

    /* special case of stacking a group of series */
    if (!strncmp(s, "stack(", 6)) {
	genr->err = dataset_stack_variables(pZ, pdinfo, genr->lhs, 
					    s, genr->prn);
	genr->done = 1;
	return genr;
    }

    /* special case of generating observation labels */
    if (!strcmp(genr->lhs, "markers")) {
	genr->err = generate_obs_markers(pZ, pdinfo, s);
	genr->done = 1;
	return genr;
    }

    /* pre-process special operators */
    if ((genr->err = catch_special_operators(genr, s))) {
	return genr;
    }

    DPRINTF(("after catch_special_operators: s='%s'\n", s));

    /* basic memory allocation */
    if (!genr_is_matrix(genr)) {
	if ((genr->err = genr_add_xvec(genr))) {
	    return genr;
	}
    } else {
	/* matrix genr: handle tranpose of compound matrices */
	genr->err = reposition_transpose_symbol(s);
	if (genr->err) {
	    return genr;
	}
    }	

    /* impose operator hierarchy */
    if (parenthesize(s)) {
	fprintf(stderr, "genr: parenthesize failed\n");
	genr->err = E_ALLOC;
	return genr;
    }

    DPRINTF(("after parenthesize: s='%s'\n", s));

    if (!genr->err) {
	genr->err = attach_atomset(genr);
    }

    if (!genr->err) {
	genr->err = math_tokenize(s, genr, 0);
    }

#if GENR_DEBUG
    int i = 0;
    while ((atom = pop_atom(genr))) {
	fprintf(stderr, "*** atom %d ***\n", i++);
	print_atom(atom);
    }
#endif

#if GENR_DEBUG
    fprintf(stderr, "\ngenerate: err = %d\n\n", genr->err);
#endif

    return genr;
}

void destroy_genr (GENERATOR *genr)
{
    if (genr == NULL) {
	return;
    }
    
    if (genr_is_matrix(genr) && genr->err) {
	const gretl_matrix *R = get_matrix_by_name(genr->lhs, genr->pdinfo);

	/* make sure we don't delete the matrix we were trying to
	   assign to, if it pre-existed this failed call */
	if (R != NULL) {
	    atom_stack_nullify_matrix(R, genr);
	}
    }

    DPRINTF(("destroy_genr: freeing atom stack\n"));
    destroy_atom_stack(genr);

    DPRINTF(("genrfree: freeing %d vars\n", genr->tmpv));
    if (genr->tmpv > 0 && genr->tmpZ != NULL) {
	int i;

	for (i=0; i<genr->tmpv; i++) {
	    free(genr->tmpZ[i]);
	}
	free(genr->tmpZ);
    }

    if (genr->xvec != NULL) {
	free(genr->xvec);
    }

    if (genr->mstack != NULL) {
	free(genr->mstack);
    }

    free(genr);
}

static int genr_write_matrix (GENERATOR *genr)
{
    gretl_matrix *M = matrix_calc_pop(genr);
    gretl_matrix *Mptr = M;
    gretl_matrix *R = NULL;
    const char *mask = genr->label;
    int err = 0;

    if (M == NULL) {
	err = 1;
    } else {
	err = add_or_replace_user_matrix(M, genr->varname, mask, &R,
					 genr->pZ, genr->pdinfo);

	if (R != NULL) {
	    /* avoid double-freeing replaced matrix R */
	    atom_stack_nullify_matrix(R, genr);
	}

	if (genr_is_submatrix(genr)) { 
	    /* and avoid double-freeing temporary matrix M */
	    atom_stack_nullify_matrix(Mptr, genr);
	}
	    
	if (!err) {
	    if (R != NULL) {
		if (*genr->label != '\0') { 
		    DPRINTF(("genr_write_matrix: rewrote matrix '%s' at %p "
			     "using submatrix spec %s\n", genr->varname,
			     (void *) R, genr->label));
		} else {
		    DPRINTF(("genr_write_matrix: replaced matrix '%s' at %p "
			     "with new matrix at %p\n", genr->varname,
			     (void *) R, (void *) M));
		}
	    } else {
		DPRINTF(("genr_write_matrix: added matrix %p as '%s'\n", 
			 (void *) M, genr->varname));
	    }
	} else {
	    DPRINTF(("genr_write_matrix: err = %d\n", err));
	}
    }

    return err;
}

static void finalize_genr (GENERATOR *genr, int oldv)
{
    if (!genr_doing_save(genr)) {
	eval_msg(genr);
    } else {
	make_genr_varname(genr, genr->lhs);
	if (genr_depositing_result(genr)) {
	    genr_deposit_result(genr);
	} else if (genr_is_matrix(genr)) {
	    genr->err = genr_write_matrix(genr);
	} else {
	    genr->err = genr_write_var(genr, genr->pZ);
	    if (!genr->err && !genr_is_private(genr)) {
		if (!genr_single_obs(genr)) {
		    write_genr_label(genr, oldv);
		}
		compose_genr_msg(genr, oldv);
	    }
	}
    } 
}

int execute_genr (GENERATOR *genr, int oldv)
{
    if (genr->err) {
	/* reset the error in case we're recomputing a compiled
	   generator */
	genr->err = 0;
    }

    if (genr_is_matrix(genr)) {
	evaluate_matrix_genr(genr);
    } else {
	evaluate_genr(genr);
    }

    if (!genr->err) {
	finalize_genr(genr, oldv);
    }

    return genr->err;
}

/**
 * genr_get_err:
 * @genr: pointer to variable-generator struct.
 *
 * Returns: the error code on @genr, or %E_ALLOC if @genr is
 * %NULL (the code will be 0 on success).
 */

int genr_get_err (const GENERATOR *genr)
{
    if (genr == NULL) {
	return E_ALLOC;
    } else {
	return genr->err;
    }
}

/* below: mechanism for depositing and retrieving a scalar
   or matrix result (bypassing the usual mechanism whereby
   genr adds its result to the dataset, or to the stack
   of named matrices) 
*/

static double generated_x = NADBL;
static gretl_matrix *generated_m;

static void genr_deposit_result (GENERATOR *genr)
{
    if (genr_is_matrix(genr)) {
	if (genr->err) {
	    generated_m = NULL;
	} else {
	    generated_m = matrix_calc_pop(genr);
	    atom_stack_nullify_matrix(generated_m, genr);
	}
	generated_x = NADBL;
    } else {
	if (genr->err) {
	    generated_x = NADBL;
	} else {
	    generated_x = genr->xvec[genr->pdinfo->t1];
	}
	generated_m = NULL;
    }	
}

static int genr_retrieve_result (double *x, gretl_matrix **m)
{
    int err = 0;

    if (x != NULL) {
	if (na(generated_x)) {
	    err = 1;
	}
	*x = generated_x;
	generated_x = NADBL;
    } 

    if (m != NULL) {
	if (generated_m == NULL) {
	    err = 1;
	}
	*m = generated_m;
	generated_m = NULL;
    }

    return err;
}

/* For generating scalars only (so far): we hit an error
   on the first try, but it could be because the RHS is
   a compound matrix-generating expression that yields a
   1 x 1 matrix on evaluation via the matrix track.  Try
   again on this assumption.
*/

static int genr_try_matrix_expression (GENERATOR *genr, int oldv)
{
    char *p = strchr(genr->orig_s, '=');
    gretl_matrix *g = NULL;
    int err = 0;

    /* reset */
    genr->err = 0;

    if (p == NULL) {
	err = 1;
    }

    if (!err) {
	err = get_generated_matrix(p + 1, &g, genr->pZ, genr->pdinfo);
    }

    if (g != NULL) {
	if (gretl_matrix_is_scalar(g)) {
	    genr->xvec[0] = gretl_vector_get(g, 0);
	    genr_set_scalar(genr);
	} else {
	    err = 1;
	}
	gretl_matrix_free(g);
    }

    if (!err) {
	finalize_genr(genr, oldv);
    }

    if (err && !genr->err) {
	genr->err = err;
    }

    return genr->err;
}
    
/**
 * generate:
 * @line: command line (formula for generating variable).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: option flags: %OPT_P means we're generating a "private"
 * variable; %OPT_M means the result should be a matrix.
 * @prn: printing struct.
 *
 * Generates a new variable, usually via some transformation of
 * existing variables, or by retrieving an internal variable associated
 * with the estimation of the last model.
 * 
 * Returns: 0 on success, integer error code on error.
 */

int generate (const char *line, double ***pZ, DATAINFO *pdinfo, 
	      gretlopt opt, PRN *prn)
{
    GENERATOR *genr;
    int oldv = pdinfo->v;
    int err = 0;

    genr = genr_compile(line, pZ, pdinfo, opt, prn);
    err = genr_get_err(genr);

    if (!genr->done) {
	if (!err) {
	    err = execute_genr(genr, oldv);
	} else if (err == E_MATVAR) {
	    err = genr_try_matrix_expression(genr, oldv);
	}
    }

    if (genr != NULL) {
	destroy_genr(genr);
    }

    return err;
}

/* simple accessor, used in nls.c */

int genr_get_varnum (const GENERATOR *genr)
{
    return genr->varnum;
}

static int genr_write_var (GENERATOR *genr, double ***pZ)
{
    DATAINFO *pdinfo = genr->pdinfo;
    double xt = genr->xvec[pdinfo->t1];
    int t, v = genr->varnum;
    int modifying = 0;
    int err = 0;

    /* check that any request for a scalar can be honored,
       abort if not */
    if (genr_require_scalar(genr) && genr_is_series(genr)) {
	strcpy(gretl_errmsg, _("Specified a scalar, but a series "
			       "was produced."));
	return 1;
    }	

    /* take into account any forcing the user has attempted */
    if (genr_forcing_series(genr) && genr_is_scalar(genr)) {
	genr_unset_scalar(genr);
	for (t=0; t<pdinfo->n; t++) {
	    genr->xvec[t] = xt;
	}
    } 

    if (v >= genr->pdinfo->v) {
	/* the generated var is a new addition to the data set */
	if (genr_is_series(genr)) {
	    err = dataset_add_series(1, pZ, pdinfo);
	} else {
	    err = dataset_add_scalar(pZ, pdinfo);
	}
	if (!err) {
	    strcpy(pdinfo->varname[v], genr->varname);
	}
    } else {
	/* we're modifying an existing variable */
	modifying = 1;
    }

    /* generally we allow scalar -> series expansion for existing
       vars, but not if the dataset is subsampled */
    if (modifying && complex_subsampled() &&
	genr_is_series(genr) && !pdinfo->vector[v]) {
	strcpy(gretl_errmsg, _("You cannot turn a scalar into a vector "
			       "when the dataset is subsampled."));
	err = 1;
    }

#if GENR_DEBUG
    if (err) {
	fprintf(stderr, "genr_write_var: err = %d\n", err);
    } else {
	fprintf(stderr, "genr_write_var: adding %s '%s' (#%d, %s)\n",
		(pdinfo->vector[v])? "vector" : "scalar",
		pdinfo->varname[v], v,
		(modifying)? "replaced" : "newly created");
    }
#endif

    if (!err) {
	if (genr->obs >= 0) {
	    /* we're replacing a single observation */
	    if (genr_is_scalar(genr)) {
		(*pZ)[v][genr->obs] = xt;
	    } else {
		(*pZ)[v][genr->obs] = genr->xvec[genr->obs];
	    }
	} else if (genr_is_scalar(genr)) {
	    if (pdinfo->vector[v]) {
		/* we never allow series -> scalar conversion for
		   existing vars, so expand the result */
		for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
		    (*pZ)[v][t] = xt;
		}
	    } else {
		/* just transcribe scalar */
		(*pZ)[v][0] = xt;
		if (!genr_is_private(genr)) {
		    /* don't write labels for "private" vars */
		    strcat(VARLABEL(pdinfo, v), _(" (scalar)"));
		}
	    }
	} else {
	    /* we generated, and will now transcribe, a series result */
	    if (modifying && !pdinfo->vector[v]) {
		err = dataset_scalar_to_vector(v, pZ, pdinfo);
	    }

	    if (!modifying) {
		/* initialize all vals to missing */
		for (t=0; t<pdinfo->n; t++) {
		    (*pZ)[v][t] = NADBL;
		}
	    }

	    if (!err) {
		/* transcribe the generated values */
		for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
		    (*pZ)[v][t] = genr->xvec[t];
		}
	    }

	    if (!err && genr->S != NULL) {
		/* handle sorted observation markers */
		if (genr_simple_sort(genr)) {
		    set_sorted_markers(pdinfo, v, genr->S);
		    genr->S = NULL;
		} else {
		    free_genr_S(genr);
		}
	    }
	}
    } /* !err */

    return err;
}

static double calc_xy (double x, double y, char op, int t, int *err) 
{
    long int ny;
    double xx, yy;

    *err = 0;

#if GENR_DEBUG
    fprintf(stderr, "calc_xy: in: x=%g, y=%g, ", x, y);
    if (isprint(op)) fprintf(stderr, "op='%c'\n", op);
    else fprintf(stderr, "op=%d\n", op);
#endif

    /* special case: 0.0 * anything (including even NA) = 0.0 */
    if (op == '*' && (x == 0.0 || y == 0.0)) {
	return 0.0;
    }

    /* otherwise, NA propagates to the result */
    if (op && (na(x) || na(y))) {
	return NADBL;
    }

    switch (op) {
    case '\0':
	x = y;
	break;
    case '+':
	x += y;
	break;
    case '|':
	if (floatneq(x, 0.) || floatneq(y, 0.)) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case '-':
	x -= y;
	break;
    case '*':
	x *= y;
	break;
    case '&':
	x *= y;
	if (x != 0.) {
	    x = 1.0;
	}
	break;
    case '%':
	x = (double) ((int) x % (int) y);
	break;
    case '/':
	if (floateq(y, 0.0)) { 
	    sprintf(gretl_errmsg, _("Zero denominator for obs %d"), t + 1);
	    x = NADBL;
	    *err = 1;
	} else {
	    x /= y;
	}
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
	    x = NADBL;
	    *err = 1;
	} else if (floateq(xx, 0.0)) {
	    x = 0.0;
	} else {
	    x = pow(xx, yy);
	}
	break;
    case '<':
	if (x < y) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case '>':
	if (x > y) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case '=':
	if (floateq(x, y)) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case OP_NEQ: /* not equals */
	if (floateq(x, y)) {
	    x = 0.0;
	} else {
	    x = 1.0;
	}
	break;
    case OP_GTE: /* greater than or equal */
	if (floateq(x, y)) {
	    x = 1.0;
	} else if (x > y) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case OP_LTE: /* less than or equal */
	if (floateq(x, y)) {
	    x = 1.0;
	} else if (x < y) {
	    x = 1.0;
	} else {
	    x = 0.0;
	}
	break;
    case '!':
	if (floatneq(y, 0.0)) {
	    x = 0.0;
	} else {
	    x = 1.0;
	}
	break;
    default:
	fprintf(stderr, "unsupported operator\n");
	*err = 1;
	break;
    } 

    DPRINTF(("calc_xy: out: x=%g\n", x));

    return x;
}

/* below: math functions taking scalar arg and returning scalar */

static double evaluate_math_function (double arg, int fn, int *err)
{
    double x = NADBL;

    if (na(arg)) {
	return x;
    }

    *err = 0;

    switch (fn) {

    case T_LOG:
	if (arg <= 0.0) {
	    fprintf(stderr, "genr: log arg = %g\n", arg);
	    *err = E_LOGS;
	} else {
	    x = log(arg);
	}
	break;	
    case T_EXP:
	if (exp(arg) == HUGE_VAL) {
	    fprintf(stderr, "genr: excessive exponent = %g\n", arg);
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
    case T_CNORM:
	x = normal_cdf(arg);
	break;
    case T_DNORM:
	x = normal_pdf(arg);
	break;
    case T_QNORM:
	x = ndtri(arg);
	break;
    case T_GAMMA:
	x = cephes_gamma(arg);
	break;
    case T_LNGAMMA:
	x = cephes_lgamma(arg);
	break;
    default:
	break;
    }

    return x;
}

static double *get_mp_series (const char *s, GENERATOR *genr,
			      int fn, int *err)
{
    double *x = malloc(genr->pdinfo->n * sizeof *x);

    if (x == NULL) {
	return NULL;
    }

    if (fn == T_MPOW) {
	*err = genr_mpow(s, x, *genr->pZ, genr->pdinfo);
    }
#ifdef HAVE_MPFR
    else if (fn == T_MLOG) {
	*err = genr_mlog(s, x, *genr->pZ, genr->pdinfo);
    }
#endif

    if (*err) {
	*err = E_INVARG;
    }

    return x;
}

static double *get_random_series (DATAINFO *pdinfo, int fn)
{
    double *x = malloc(pdinfo->n * sizeof *x);

    if (x == NULL) {
	return NULL;
    }

    if (fn == T_NORMAL) {
	gretl_normal_dist(x, pdinfo->t1, pdinfo->t2);
    } else if (fn == T_UNIFORM) {
	gretl_uniform_dist(x, pdinfo->t1, pdinfo->t2);
    }

    return x;
}

/* below: functions taking series (z) as input and returning a scalar
   statistic.  The whole series must be evaluated before these stats
   can be calculated.
*/

static double evaluate_statistic (double *z, GENERATOR *genr, int fn)
{
    double x = NADBL;
    double *tmp = NULL;
    int i, t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;

    if (fn == T_NOBS || fn == T_T1 || fn == T_T2) {
	/* special: they don't need a series allocated */
	if (fn == T_NOBS) {
	    i = 0;
	    for (t=t1; t<=t2; t++) {
		if (!na(z[t])) i++;
	    }
	    return (double) i;
	} else if (fn == T_T1) {
	    for (t=0; t<genr->pdinfo->n; t++) {
		if (!na(z[t])) {
		    break;
		}
	    }
	    return (double) (t + 1);
	} else {
	    for (t=genr->pdinfo->n - 1; t>=0; t--) {
		if (!na(z[t])) {
		    break;
		}
	    }
	    return (double) (t + 1);
	}
    }

    tmp = malloc((t2 - t1 + 1) * sizeof *tmp);
    if (tmp == NULL) {
	genr->err = E_ALLOC;
	return x;
    }

    i = -1;
    for (t=t1; t<=t2; t++) {
	if (!na(z[t])) {
	    tmp[++i] = z[t];
	}
    }

    if (fn == T_MEAN) {
	x = gretl_mean(0, i, tmp);
    } else if (fn == T_SUM) {
	x = gretl_mean(0, i, tmp);
	x *= (i + 1);
    } else if (fn == T_SD) {
	x = gretl_stddev(0, i, tmp);
    } else if (fn == T_VAR) {
	x = gretl_variance(0, i, tmp);
    } else if (fn == T_LRVAR) {
	int pd = genr->pdinfo->pd;

	x = gretl_long_run_variance(0, i, tmp, 2 * pd);
    } else if (fn == T_SST) {
	x = gretl_sst(0, i, tmp);
    } else if (fn == T_MEDIAN) {
	x = gretl_median(0, i, tmp);
    } else if (fn == T_GINI) {
	x = gretl_gini(0, i, tmp);
    } else if (fn == T_MIN || fn == T_MAX) {
	double min, max;

	gretl_minmax(0, i, tmp, &min, &max);
	x = (fn == T_MIN)? min : max;
    } 

    free(tmp);
    
    return x;
}

static double evaluate_pvalue (const char *s, const double **Z,
			       const DATAINFO *pdinfo, int *err)
{
    double x = batch_pvalue(s, Z, pdinfo, NULL, OPT_G);

    if (na(x)) {
	*err = E_INVARG;
    }

    return x;
}

static double evaluate_critval (const char *s, const double **Z,
				const DATAINFO *pdinfo, int *err)
{
    double x = genr_get_critical(s, Z, pdinfo);

    if (na(x) || x == -1.0) {
	*err = E_INVARG;
    }

    return x;
}

static double 
evaluate_bivariate_statistic (const char *s, GENERATOR *genr, 
			      int fn)
{
    double x = NADBL;

    if (fn == T_COV || fn == T_CORR) {
	x = genr_cov_corr(s, genr->pZ, genr->pdinfo, fn);
    }

    if (na(x)) {
	genr->err = E_INVARG;
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
    } else if (fn == T_OK) {
	/* check whether obs is present or not */
	x = (na(arg))? 0.0 : 1.0;
    } else if (fn == T_MISSZERO) {
	/* change missing obs to zero */
	x = (na(arg))? 0.0 : arg;
    } else if (fn == T_ZEROMISS) {
	/* change zero to missing obs */
	x = (floateq(arg, 0.0))? NADBL : arg;
    }

    return x;
}

/* apparatus for saving sorted case markers */

struct val_mark {
    double x;
    char mark[OBSLEN];
};

static int compare_vms (const void *a, const void *b)
{
    const struct val_mark *va = (const struct val_mark *) a;
    const struct val_mark *vb = (const struct val_mark *) b;
     
    return (va->x > vb->x) - (va->x < vb->x);
}

static int inverse_compare_vms (const void *a, const void *b)
{
    const struct val_mark *va = (const struct val_mark *) a;
    const struct val_mark *vb = (const struct val_mark *) b;
     
    return (vb->x > va->x) - (vb->x < va->x);
}

static void free_genr_S (GENERATOR *genr)
{
    if (genr->S != NULL) {
	int i, n = genr->pdinfo->n;

	for (i=0; i<n; i++) {
	    if (genr->S[i] != NULL) {
		free(genr->S[i]);
	    }
	}
	free(genr->S);
	genr->S = NULL;
    }
}

static int allocate_genr_S (GENERATOR *genr)
{
    int i, n = genr->pdinfo->n;
    int err = 0;

    genr->S = malloc(n * sizeof *genr->S);

    if (genr->S == NULL) {
	err = 1;
    } else {
	for (i=0; i<n; i++) {
	    genr->S[i] = NULL;
	}
	for (i=0; i<n; i++) {
	    genr->S[i] = malloc(OBSLEN);
	    if (genr->S[i] == NULL) {
		err = 1;
		free_genr_S(genr);
		break;
	    }
	    genr->S[i][0] = '\0';
	}
    }

    return err;
}

/* do a simple sort if there are no case markers in the dataset,
   but if there are case markers, keep a record of the sorted
   markers and attach it to the newly generated variable
*/

static int 
sort_series (const double *mvec, double *x, GENERATOR *genr, int fn)
{
    DATAINFO *pdinfo = genr->pdinfo;
    double *tmp = NULL;
    struct val_mark *vm = NULL;
    int markers = 0;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int i, t;

    if (pdinfo->S != NULL && !complex_subsampled()) {
	markers = 1;
    }

    if (markers) {
	allocate_genr_S(genr);
	if (genr->S == NULL) {
	    markers = 0;
	}
    }

    if (markers) {
	vm = malloc(T * sizeof *vm);
	if (vm == NULL) {
	    free_genr_S(genr);
	    return 1;
	}
    } else {
	tmp = malloc(T * sizeof *tmp);
	if (tmp == NULL) {
	    return 1;
	}
    }

    i = -1;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(mvec[t])) {
	    continue;
	}
	++i;
	if (markers) {
	    vm[i].x = mvec[t];
	    strcpy(vm[i].mark, pdinfo->S[t]);
	} else {
	    tmp[i] = mvec[t];
	}
    }

    if (fn == T_DSORT) {
	if (markers) {
	    qsort(vm, i + 1, sizeof *vm, inverse_compare_vms);
	} else {
	    qsort(tmp, i + 1, sizeof *tmp, gretl_inverse_compare_doubles);
	}
    } else {
	if (markers) {
	    qsort(vm, i + 1, sizeof *vm, compare_vms);
	} else {
	    qsort(tmp, i + 1, sizeof *tmp, gretl_compare_doubles);
	}
    }

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(mvec[t])) {
	    x[t] = NADBL;
	    continue;
	} else if (markers) {
	    x[t] = vm[i].x;
	    strcpy(genr->S[t], vm[i].mark);
	} else {
	    x[t] = tmp[i];
	}
	i++;
    }

    if (tmp != NULL) {
	free(tmp);
    }

    if (vm != NULL) {
	free(vm);
    }

    return 0;
}

/* below: create a temporary series after evaluating the full-
   length argument. 
*/

static double *get_tmp_series (double *mvec, GENERATOR *genr, 
			       int fn, double param)
{
    DATAINFO *pdinfo = genr->pdinfo;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2; 
    double *x;
    double xx, yy;

#if GENR_DEBUG
    fprintf(stderr, "*** Doing get_tmp_series, fn = %d ***\n", fn);
#endif

    if (fn == T_SDIFF && !dataset_is_seasonal(pdinfo)) {
	genr->err = E_PDWRONG;
	return NULL;
    }

    x = malloc(pdinfo->n * sizeof *x); 
    if (x == NULL) {
	return NULL;
    }

    if (fn == T_DIFF || fn == T_LDIFF || fn == T_SDIFF) {
	int t0 = (fn == T_SDIFF)? pdinfo->pd : 1;

	if (t1 < t0) {
	    for (t=t1; t<t0; t++) {
		x[t] = NADBL;
	    }
	    t1 = t0;
	}
	
	for (t=t1; t<=t2; t++) {
	    /* get "later" value */
	    if (pdinfo->structure == STACKED_TIME_SERIES &&
		panel_unit_first_obs(t, pdinfo)) {
		x[t] = NADBL;
		continue;
	    }
	    xx = mvec[t];

	    /* get "earlier" value */
	    if (pdinfo->structure == STACKED_CROSS_SECTION) {
		yy = (t - pdinfo->pd >= 0)? mvec[t-pdinfo->pd] : NADBL;
	    } else {
		yy = mvec[t - t0];
	    }

	    if (na(xx) || na(yy)) {
		x[t] = NADBL;
		continue;
	    }

	    /* perform the differencing */
	    if (fn == T_DIFF || fn == T_SDIFF) {
		x[t] = xx - yy;
	    } else {
		/* log difference */
		if (xx <= 0.0 || yy <= 0.0) {
		    genr->err = E_LOGS;
		    x[t] = NADBL;
		} else {
		    x[t] = log(xx) - log(yy);
		}
	    }
	}
    } else if (fn == T_CUM) {
	int t0 = t1;

	for (t=t1; t<=t2; t++) {
	    if (na(mvec[t])) {
		x[t0++] = NADBL;
	    } else {
		break;
	    }
	}
	if (t0 < t2) {
	    x[t0] = (na(mvec[t0])) ? 0.0 : mvec[t0];
	    for (t=t0+1; t<=t2; t++) {
		if (na(mvec[t])) {
		    x[t] = x[t-1];
		} else {
		    x[t] = x[t-1] + mvec[t];
		}
	    }
	}
    } else if (fn == T_SORT || fn == T_DSORT) {
	genr->err = sort_series(mvec, x, genr, fn);
	if (genr->err) {
	    free(x);
	    x = NULL;
	}
    } else if (fn == T_RESAMPLE) {
	int i, n, rt1 = t1, rt2 = t2;
	double *tmp = NULL;

	array_adjust_t1t2(mvec, &rt1, &rt2);

	n = rt2 - rt1 + 1;
	if (n <= 1) {
	    genr->err = E_DATA;
	    free(x);
	    return NULL;
	}

	tmp = malloc(n * sizeof *tmp);

	if (tmp == NULL) {
	    free(x);
	    return NULL;
	}

	for (t=t1; t<=t2; t++) {
	    if (t < rt1 || t > rt2) {
		x[t] = NADBL;
	    }
	}

	/* generate uniform random series */
	gretl_uniform_dist(tmp, 0, n - 1);

	/* sample from source series based on indices */
	for (t=rt1; t<=rt2; t++) {
	    i = rt1 + n * tmp[t-rt1];
	    if (i > rt2) i = rt2;
	    x[t] = mvec[i];
	}

	free(tmp);
    } else if (fn == T_HPFILT) { 
	genr->err = hp_filter(mvec, x, pdinfo);	
    } else if (fn == T_BKFILT) { 
	genr->err = bkbp_filter(mvec, x, pdinfo);	
    } else if (fn == T_FRACDIFF) {
	genr->err = get_fracdiff(mvec, x, param, pdinfo);
    }

    return x;
}

static double get_dataset_statistic (DATAINFO *pdinfo, int idx)
{
    double x = NADBL;

    if (pdinfo == NULL) return x;

    if (idx == R_NOBS) {
	x = (double) (pdinfo->t2 - pdinfo->t1 + 1);
    } else if (idx == R_PD) {
	x = (double) pdinfo->pd;
    } else if (idx == R_NVARS) {
	x = (double) pdinfo->v;
    }

    return x;
}

static double get_test_stat_value (char *label, int idx)
{
    double x = NADBL;

    if (idx == R_TEST_PVAL) {
	x = get_last_pvalue(label);
    } else if (idx == R_TEST_STAT) {
	x = get_last_test_statistic(label);
    }

    return x;
}

static int 
get_obs_value (const char *s, const double **Z, const DATAINFO *pdinfo, 
	       int *v, int *obsnum)
{
    char vname[VNAMELEN], obs[OBSLEN];
    int err = 0;

    if (sscanf(s, "%15[^[][%10[^]]]", vname, obs) == 2) {
	int t, i = varindex(pdinfo, vname);

	if (i < pdinfo->v && pdinfo->vector[i]) {
	    t = get_t_from_obs_string(obs, Z, pdinfo);
	    if (t >= 0 && t < pdinfo->n) {
		*v = i;
		*obsnum = t;
	    } else {
		err = 1;
	    }
	} else {
	    err = 1;
	}
    } else {
	err = 1;
    }

    return err;
}

static double get_tnum (const DATAINFO *pdinfo, int t)
{
    if (pdinfo->structure == TIME_SERIES && pdinfo->pd == 1) {
	/* annual data: let 't' be the year */ 
	return pdinfo->sd0 + t;
    } else {
	/* let 't' be the 1-based observation number */
	return (double) (t + 1);
    }
}

#define GEN_LEVEL_DEBUG 0

/**
 * varindex:
 * @pdinfo: data information struct.
 * @varname: name of variable to test.
 *
 * Returns: the ID number of the variable whose name is given,
 * or the next available ID number if there is no variable of
 * that name.
 */

int varindex (const DATAINFO *pdinfo, const char *varname)
{
    const char *check = varname;
    int fsd = 0;
    int i, ret = pdinfo->v;

    if (varname == NULL) {
	return ret;
    }

    if (isdigit(*check)) {
	return ret;
    } 

    while (*check == '_') {
	check++;
    }

    if (!strcmp(check, "t") || !strcmp(check, "obs")) {
	return TNUM;
    } else if (!strcmp(check, "const") || !strcmp(check, "CONST")) {
	return 0;
    }

    if (*(check + 1) == 0 && is_active_index_loop_char(*check)) {
	/* single-char loop index variable: 'i', 'j' or such */
	return INDEXNUM;
    }

    if (gretl_executing_function()) {
	fsd = gretl_function_stack_depth();
    }

#if GEN_LEVEL_DEBUG
    fprintf(stderr, "varindex for '%s': fsd=%d\n", check, fsd);
#endif

    if (fsd > 0) {
	/* inside a function: see only vars at that level */
	for (i=1; i<pdinfo->v; i++) { 
	    if (STACK_LEVEL(pdinfo, i) == fsd && 
		!strcmp(pdinfo->varname[i], check)) {
		ret = i;
		break;
	    }
	}
    } else {
	/* see all vars */
	for (i=1; i<pdinfo->v; i++) { 
	    if (!strcmp(pdinfo->varname[i], check)) { 
		ret = i;
		break;
	    }
	}
    }

#if GEN_LEVEL_DEBUG
    fprintf(stderr, "varindex for '%s': returning %d\n", check, ret);
#endif 

    return ret;
}

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

    mp_raise = get_plugin_function("mp_vector_raise_to_power", &handle);
    if (mp_raise == NULL) {
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
    
    if (sscanf(str, "%15s", vname) != 1) {
	return 1;
    }

    v = varindex(pdinfo, vname);
    if (v >= pdinfo->v) {
	return 1;
    } 

    mp_log = get_plugin_function("mp_vector_ln", &handle);
    if (mp_log == NULL) {
	return 1;
    }

    err = mp_log(Z[v], xvec, pdinfo->n);

    close_plugin(handle);
    
    return err;
}

#endif /* HAVE_MPFR */

static void eval_msg (const GENERATOR *genr)
{
    double x = genr->xvec[genr->pdinfo->t1];

    if (na(x)) {
	pputs(genr->prn, " NA\n");
    } else {
	pprintf(genr->prn, " %g\n", x);
    }
}

static void compose_genr_msg (const GENERATOR *genr, int oldv)
{
    double x;
    int scalar = genr_is_scalar(genr);
    int mutant = 0;

    if (genr->prn == NULL || !gretl_messages_on()) {
	return;
    }

    /* no message for special internal vars */
    if (!strcmp(genr->varname, "argv") ||
	!strncmp(genr->varname, "$nl", 3) ||
	!strcmp(genr->varname, "tmpmsk")) {
	return;
    }

    if (genr->varnum < oldv) {
	if (genr->pdinfo->vector[genr->varnum]) {
	    scalar = 0;
	} else if (!scalar) {
	    mutant = 1;
	}
    }

    pprintf(genr->prn, "%s %s %s (ID %d)", 
	    (genr->obs >= 0)? _("Modified") :
	    (genr->varnum < oldv)? _("Replaced") : _("Generated"), 
	    (mutant)? _("variable") :
	    (scalar)? _("scalar") : _("vector"),
	    genr->varname, genr->varnum);

    if (scalar) {
	x = genr->xvec[genr->pdinfo->t1];
	if (na(x)) {
	    pputs(genr->prn, " = NA");
	} else {
	    pprintf(genr->prn, " = %g", x);
	}
    }

    if (genr_warn(genr)) {
	pputc(genr->prn, '\n');
	pputs(genr->prn, gretl_errmsg);
    }

    pputc(genr->prn, '\n');
}
