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
#include "libset.h"
#include "modelspec.h"
#include "objstack.h"

#include "../../cephes/libprob.h"

#include <time.h>
#include <errno.h>

enum {
    HNUM = 2000,
    UHATNUM,
    YHATNUM,
    TNUM,
    OBSBOOLNUM,
    INDEXNUM
} genr_numbers;

static double calc_xy (double x, double y, char op, int t, int *err);

static void genrfree (GENERATE *genr);
static double genr_cov_corr (const char *str, double ***pZ,
			     const DATAINFO *pdinfo, int fn);
static double genr_vcv (const char *str, const DATAINFO *pdinfo, 
			MODEL *pmod);
static int genr_mpow (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
static int obs_num (const char *s, const DATAINFO *pdinfo);

#ifdef HAVE_MPFR
static int genr_mlog (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
#endif

static void eval_msg (const GENERATE *genr);
static void compose_genr_msg (const GENERATE *genr, int oldv);
static int listpos (int v, const int *list);
static int genr_write_var (double ***pZ, DATAINFO *pdinfo, 
			   GENERATE *genr);

static int math_tokenize (char *s, GENERATE *genr, int level);
static double get_obs_value (const char *s, const double **Z, 
			     const DATAINFO *pdinfo, int *err);
static int model_vector_index (const char *s);
static int dataset_var_index (const char *s);
static int test_stat_index (const char *s);
static int op_level (int c);

static double *get_model_series (const DATAINFO *pdinfo,
				 const MODEL *pmod, int v);
static double *get_random_series (DATAINFO *pdinfo, int fn);
static double *get_mp_series (const char *s, GENERATE *genr,
			      int fn, int *err);
static double *get_tmp_series (double *x, GENERATE *genr,
			       int fn, double param);

static double get_tnum (const DATAINFO *pdinfo, int t);
static double evaluate_statistic (double *z, GENERATE *genr, int fn);
static double get_model_data_element (const char *s, GENERATE *genr,
				      MODEL *pmod, int idx);
static double get_dataset_statistic (DATAINFO *pdinfo, int idx);
static double get_test_stat_value (char *label, int idx);
static double evaluate_math_function (double arg, int fn, int *err);
static double evaluate_missval_func (double arg, int fn);
static double evaluate_bivariate_statistic (const char *s, 
					    GENERATE *genr, 
					    int fn);
static double evaluate_pvalue (const char *s, const double **Z,
			       const DATAINFO *pdinfo, int *err);
static double evaluate_critval (const char *s, const double **Z,
				const DATAINFO *pdinfo, int *err);
static void free_genr_S (GENERATE *genr);

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
    { T_MEAN,     "mean" }, 
    { T_SD,       "sd" }, 
    { T_MIN,      "min" },
    { T_MAX,      "max" },
    { T_SORT,     "sort" }, 
    { T_INT,      "int" }, 
    { T_LN,       "ln" }, 
    { T_COEFF,    "coeff" },
    { T_ABS,      "abs" }, 
    { T_RHO,      "rho" }, 
    { T_SQRT,     "sqrt" }, 
    { T_SUM,      "sum" }, 
    { T_NOBS,     "nobs" },
    { T_T1,       "firstobs" },
    { T_T2,       "lastobs" },
    { T_NORMAL,   "normal" }, 
    { T_UNIFORM,  "uniform" }, 
    { T_STDERR,   "stderr" },
    { T_CUM,      "cum" }, 
    { T_MISSING,  "missing" },
    { T_OK,       "ok" },        /* opposite of missing */
    { T_MISSZERO, "misszero" },
    { T_CORR,     "corr" },
    { T_VCV,      "vcv" },
    { T_VAR,      "var" },
    { T_SST,      "sst" },
    { T_COV,      "cov" },
    { T_MEDIAN,   "median" },
    { T_ZEROMISS, "zeromiss" },
    { T_PVALUE,   "pvalue" },
    { T_CRIT,     "critical" },
    { T_OBSNUM,   "obsnum" },
    { T_MPOW,     "mpow" },
    { T_DNORM,    "dnorm" },
    { T_CNORM,    "cnorm" },
    { T_GAMMA,    "gamma" },
    { T_LNGAMMA,  "lngamma" },
    { T_RESAMPLE, "resample" },
    { T_HPFILT,   "hpfilt" },    /* Hodrick-Prescott filter */
    { T_BKFILT,   "bkfilt" },    /* Baxter-King filter */
    { T_FRACDIFF, "fracdiff" },  /* fractional difference */
    { T_VARNUM,   "varnum" },    /* variable's ID number from its name */
    { T_VECTOR,   "isvector" },
    { T_ISLIST,   "islist" },
    { T_NELEM,    "nelem" },
#ifdef HAVE_MPFR
    { T_MLOG,     "mlog" },
#endif
    { T_IDENTITY, "ident" },
    { 0, NULL }
};

#define LEVELS 7

#define STANDARD_MATH(f) (f == T_LOG || f == T_LN || f == T_EXP || \
                          f == T_SIN || f == T_COS || f == T_TAN || \
                          f == T_ATAN || f == T_INT || f == T_ABS || \
                          f == T_DNORM || f == T_CNORM || f == T_SQRT || \
                          f == T_GAMMA || f == T_LNGAMMA)

#define UNIVARIATE_STAT(t) (t == T_MEAN || t == T_SD || t == T_SUM || \
                            t == T_VAR || t == T_MEDIAN || t == T_MIN || \
                            t == T_SST || t == T_MAX || t == T_NOBS || \
                            t == T_T1 || t == T_T2 || t == T_VARNUM || \
                            t == T_VECTOR || t == T_ISLIST || t == T_NELEM)

#define BIVARIATE_STAT(t) (t == T_CORR || t == T_COV)

#define MODEL_DATA_ELEMENT(f) (f == T_COEFF || f == T_STDERR || \
                               f == T_RHO || f == T_VCV)

#define MISSVAL_FUNC(f) (f == T_MISSING || f == T_OK || \
                         f == T_MISSZERO || f == T_ZEROMISS)

#define MODEL_VAR_INDEX(v) (v == HNUM || v == UHATNUM || v == YHATNUM)

#ifdef HAVE_MPFR
# define MP_MATH(f) (f == T_MPOW || f == T_MLOG)
#else
# define MP_MATH(f) (f == T_MPOW)
#endif

#if GENR_DEBUG
static const char *get_func_word (int fnum);
#endif

#define MAXTERMS  64
#define TOKLEN   128
#define ARGLEN    64

#define genr_is_scalar(g) ((g)->flags & GENR_SCALAR)
#define genr_is_vector(g) (!((g)->flags & GENR_SCALAR))
#define genr_set_scalar(g) ((g)->flags |= GENR_SCALAR)
#define genr_unset_scalar(g) ((g)->flags &= ~GENR_SCALAR)

#define genr_force_vector(g) ((g)->flags |= GENR_FORCE_VECTOR)
#define genr_forcing_vector(g) ((g)->flags & GENR_FORCE_VECTOR)

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

#define genr_single_obs(g) ((g)->obs >= 0)

static void genr_init (GENERATE *genr, double ***pZ, DATAINFO *pdinfo,
		       MODEL *pmod, gretlopt opt)
{
    genr->err = 0;
    genr->flags = GENR_SAVE | GENR_SCALAR;
    genr->xvec = NULL;
    genr->varnum = 0;
    genr->obs = -1;
    *genr->varname = '\0';
    *genr->label = '\0';
    genr->tmpv = 0;
    genr->tmpZ = NULL;
    genr->pdinfo = pdinfo;
    genr->pZ = pZ;
    genr->pmod = pmod;
    genr->aset = NULL;
    genr->S = NULL;

    if (opt & OPT_P) {
	genr_set_private(genr);
    }

    reset_calc_stack(genr);
}

static genatom *make_atom (int scalar, int varnum,
			   int lag, double val,
			   int func, char op, char *str,
			   int level, atomset *aset)
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

    fprintf(stderr, " atom->scalar = %d\n", atom->scalar);

    if (atom->varnum >= 0) {
	fprintf(stderr, " atom->varnum = %d\n", atom->varnum);
    }
    if (atom->tmpvar >= 0) {
	fprintf(stderr, " atom->tmpvar = %d\n", atom->tmpvar);
    }
    if (atom->lag > 0) {
	fprintf(stderr, " atom->lag = %d\n", atom->lag);
    }
    if (atom->scalar) {
	fprintf(stderr, " atom->val = %g\n", atom->val);
    }
    if (atom->func > 0) {
	const char *fword = get_func_word(atom->func);

	if (fword != NULL) {
	    fprintf(stderr, " atom->func = %d (%s)\n", 
		    atom->func, get_func_word(atom->func));
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

static int get_lagvar (const char *s, int *lag, GENERATE *genr)
{
    static char format[16] = {0};
    char vname[USER_VLEN];
    int m = 0, v = 0;

    if (*format == 0) {
	sprintf(format, "%%%d[^(](%%d)", USER_VLEN - 1);
    }

    if (sscanf(s, format, vname, &m) == 2) {
	v = varindex(genr->pdinfo, vname);
	DPRINTF(("get_lagvar: looking at '%s' (v == %d)\n", vname, v));
	if (v >= genr->pdinfo->v && !MODEL_VAR_INDEX(v)) {
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
	    char test[USER_VLEN + 8];

	    sprintf(test, "%s(%d)", vname, m);
	    if (strcmp(test, s)) {
		/* string s contains extra stuff */
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

int get_generated_value (const char *argv, double *val,
			 double ***pZ, DATAINFO *pdinfo,
			 MODEL *pmod, int t)
{
    char genline[MAXLINE];
    int err = 0;

    sprintf(genline, "genr argv=%s", argv);
#if GENR_DEBUG
    fprintf(stderr, "get_generated_value: trying '%s'\n", genline);
#endif
    err = generate(genline, pZ, pdinfo, pmod, OPT_P);
    if (!err) {
	int v = pdinfo->v - 1;

	if (pdinfo->vector[v]) {
	    *val = (*pZ)[v][0];
	} else {
	    *val = (*pZ)[v][t];
	}
	err = dataset_drop_last_variables(1, pZ, pdinfo);
    }

    return err;
}

/* get the values of the scalar arguments to a genr function */

static int evaluate_genr_function_args (char *s, GENERATE *genr)
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
				      genr->pdinfo, genr->pmod, 
				      0);
	}
	nf++;
    }

    if (!err) {
	char numstr[24];

	*s = '\0';
	strcat(s, st);
	for (i=0; i<nf; i++) {
	    sprintf(numstr, ",%.10g", vals[i]);
	    strcat(s, numstr);
	}
    }

    free(tmp);
	
    return err;
}

static int get_arg_string (char *str, const char *s, int func, GENERATE *genr)
{
    const char *p = strchr(s, '(');
    int n = strlen(p) - 2;
    int err = 0;

    *str = '\0';

    DPRINTF(("get_arg_string: looking at '%s'\n", s));

    if (n < 0 || n > ARGLEN - 1) {
	err = 1;
    } else {
	strncat(str, p + 1, n);
    }

    if (func == T_PVALUE || func == T_CRIT) {
	err = evaluate_genr_function_args(str, genr);
    } else if (func == T_FRACDIFF) {
	char vname[9];
	double param;

	err = sscanf(str, "%8[^,],%lf", vname, &param) != 2;
    }

    DPRINTF(("get_arg_string: got '%s'\n", str));

    return err;
}

static int catch_saved_object_scalar (const char *s, double *x)
{
    const char *p = strstr(s, ".$");
    int len = strcspn(s, ".");
    char *test;
    int err = 0, ret = 0;

    test = gretl_strndup(s, len);

    if (test != NULL) {
	*x = saved_object_get_value(test, p + 1, &err);
	if (!err) ret = 1;
	free(test);
    }

    return ret;
}

static genatom *parse_token (const char *s, char op,
			     GENERATE *genr, int level)
{
    int v = -1;
    int scalar = 0, lag = 0, func = 0;
    double val = 0.0;
    char str[ARGLEN] = {0};

    DPRINTF(("parse_token: looking at '%s'\n", s));

    if (isalpha((unsigned char) *s)) {
	if (!strchr(s, '(') && !strchr(s, '[')) {
	    if (!strcmp(s, "pi")) {
		val = M_PI;
		scalar = 1;
	    } else if (!strcmp(s, "NA")) {
		val = NADBL;
		scalar = 1;
	    } else if (strstr(s, ".$")) {
		if (catch_saved_object_scalar(s, &val)) {
		    scalar = 1;
		} else {
		    sprintf(gretl_errmsg, _("Undefined variable name '%s' in genr"), s);
		    genr->err = E_UNKVAR;
		}
	    } else {
		v = varindex(genr->pdinfo, s);

		if (v == genr->pdinfo->v) { 
		    if (get_list_by_name(s) != NULL) {
			/* name of saved list, not single variable */
			v = -1;
			strncat(str, s, ARGLEN - 1);
		    } else {
			sprintf(gretl_errmsg, _("Undefined variable name '%s' in genr"), s);
			genr->err = E_UNKVAR;
		    }
		} else {
		    DPRINTF(("recognized var '%s' (#%d)\n", s, v));
		}

		if (v == INDEXNUM) { 
		    int k = loop_scalar_read(*s);

		    val = k;
		    scalar = 1;
		} else if (v >= 0 && v < genr->pdinfo->v && !genr->pdinfo->vector[v]) {
		    /* handle regular scalar variables here */
		    val = (*genr->pZ)[v][0];
		    scalar = 1;
		}
	    }
	} else if (strchr(s, '(')) {
	    /* try for a function first */
	    lag = 0;
	    v = -1;

	    func = genr_function_from_string(s);

	    if (func) {
		DPRINTF(("recognized function #%d (%s)\n", func, 
			 get_func_word(func)));
		if (MP_MATH(func) || func == T_PVALUE || 
		    func == T_CRIT || func == T_FRACDIFF ||
		    BIVARIATE_STAT(func) || MODEL_DATA_ELEMENT(func)) {
		    genr->err = get_arg_string(str, s, func, genr);
		} 
		if (func == T_PVALUE) {
		    if (!genr->err) {
			val = evaluate_pvalue(str, (const double **) *genr->pZ,
					      genr->pdinfo,
					      &genr->err);
			scalar = 1;
		    }
		} else if (func == T_CRIT) {
		    if (!genr->err) {
			val = evaluate_critval(str, (const double **) *genr->pZ,
					       genr->pdinfo,
					       &genr->err);
			scalar = 1;
		    }
		} else if (func == T_VARNUM || func == T_VECTOR ||
			   func == T_ISLIST || func == T_NELEM) {
		    scalar = 1;
		} else if (MODEL_DATA_ELEMENT(func)) {
		    val = get_model_data_element(str, genr,
						 genr->pmod, 
						 func);
		    scalar = 1;
		} else if (BIVARIATE_STAT(func)) {
		    val = evaluate_bivariate_statistic(str, genr, func);
		    scalar = 1;
		}
	    } else {
		/* not a function: try a lagged variable */
		v = get_lagvar(s, &lag, genr);
		if (v > 0) {
		    DPRINTF(("recognized var #%d lag %d\n", v, lag));
		} else {
		    /* not a variable or function: dead end? */
		    DPRINTF(("dead end in parse_token, s='%s'\n", s));
		    genr->err = E_SYNTAX; 
		}
	    }
	} else if (strchr(s, '[')) {
	    int err;

	    val = get_obs_value(s, (const double **) *genr->pZ, genr->pdinfo,
				&err);
	    if (err) {
		DPRINTF(("dead end at get_obs_value, s='%s'\n", s));
		genr->err = E_SYNTAX; 
	    } else {
		scalar = 1;
	    }
	}
    }

    /* by now, first char is not alpha */

    else if (*s == '$') {
	int i;

	if ((i = get_lagvar(s, &lag, genr)) > 0) {
	    DPRINTF(("recognized var #%d lag %d\n", i, lag));
	    v = i;
	} else if ((i = gretl_model_stat_index(s)) > 0) {
	    DPRINTF(("recognized '%s' as model variable, index #%d\n", 
		     s, i));
	    val = gretl_model_get_scalar(genr->pmod, i, &genr->err);
	    scalar = 1;
	} else if ((i = model_vector_index(s)) > 0) { 
	    DPRINTF(("recognized '%s' as model vector, index #%d\n", 
		     s, i));
	    v = i;
	} else if ((i = dataset_var_index(s)) > 0) {
	    DPRINTF(("recognized '%s' as dataset var, index #%d\n", 
		     s, i));
	    val = get_dataset_statistic(genr->pdinfo, i);
	    scalar = 1;
	} else if ((i = test_stat_index(s)) > 0) {
	    DPRINTF(("recognized '%s' as test-related var, index #%d\n", 
		     s, i));
	    val = get_test_stat_value(genr->label, i);
	    scalar = 1;
	} else {
	    genr->err = E_UNKVAR;
	}
    } else if (numeric_string(s)) {
	val = dot_atof(s);
	scalar = 1;
    } else if (*s == '"') {
	/* observation label? (basis for dummy series) */
	val = obs_num(s, genr->pdinfo);
	if (val > 0) {
	    func = T_OBSNUM;
	} else{
	    genr->err = E_SYNTAX;
	}
    } else if (strchr(s, ':')) {
	/* time-series observation? */
	val = obs_num(s, genr->pdinfo);
	if (val > 0) {
	    scalar = 1;
	} else {
	    genr->err = E_SYNTAX;
	}
    } else {
	DPRINTF(("dead end in parse_token, s='%s'\n", s));
	genr->err = E_SYNTAX;
    }

    if (genr->err) return NULL;

    return make_atom(scalar, v, lag, val, func, op, str, 
		     level, genr->aset);
}

static double get_lag_at_obs (int v, int tmp, int lag, 
			      const GENERATE *genr, int t)
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

static double eval_atom (genatom *atom, GENERATE *genr, int t, 
			 double a)
{
    double x = NADBL;

    /* constant, scalar variable, or specific obs from a series */
    if (atom->scalar) {
	x = atom->val;
	DPRINTF(("eval_atom: got scalar = %g\n", x));
    } 

    /* temporary variable */
    else if (atom->tmpvar >= 0) {
	if (!atom->lag) {
	    x = genr->tmpZ[atom->tmpvar][t];
	    DPRINTF(("eval_atom: got temp obs = %g\n", x));
	} else {
	    x = get_lag_at_obs(atom->tmpvar, 1, atom->lag, genr, t);
	    DPRINTF(("eval_atom: got lagged temp obs (tmpvar %d, "
		    "lag %d) = %g\n", atom->tmpvar, atom->lag, x));
	}
    }

    /* trend/index variable */
    else if (atom->varnum == TNUM) {
	x = get_tnum(genr->pdinfo, t);
    }

    /* regular variable or lagged term */
    else if (atom->varnum >= 0) {
	if (!atom->lag) {
	    x = (*genr->pZ)[atom->varnum][t];
	    DPRINTF(("eval_atom: got data obs (var %d) = %g\n", 
		    atom->varnum, x));
	} else {
	    x = get_lag_at_obs(atom->varnum, 0, atom->lag, genr, t);
	    DPRINTF(("eval_atom: got lagged data obs (var %d, "
		    "lag %d) = %g\n", atom->varnum, atom->lag, x));
	}
    }

    /* function */
    else if (atom->func) {
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
    }

    return x;
}

static int genr_add_temp_var (GENERATE *genr, double *x)
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

    if (genr_add_temp_var(genr, x)) {
	free(x);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;

    return 0;
}

static int add_model_series_to_genr (GENERATE *genr, genatom *atom)
{
    double *x;

    x = get_model_series((const DATAINFO *) genr->pdinfo, 
			 genr->pmod, atom->varnum);
    if (x == NULL) return 1;

    if (genr_add_temp_var(genr, x)) {
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

    if (genr_add_temp_var(genr, x)) {
	free(x);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;

    return genr->err;
}

static double *eval_compound_arg (GENERATE *genr,
				  genatom *this_atom)
{
    int t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;
    genatom *atom;
    double *xtmp;

    if (peek_child_atom(this_atom) == NULL) {
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

	    if (0 && y == NADBL) { /* watch out? */
		x = NADBL;
	    } else {
		if (atom->level < level) {
		    x = calc_pop(genr);
		}
		x = calc_xy(x, y, atom->op, t, &err);
		if (err) {
		    genr_set_warn(genr);
		}
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
get_target_fracdiff_series (GENERATE *genr, genatom *atom,
			    double *param)
{
    char vname[9];
    double *x = NULL;
    int v, t;

    if (*atom->str == '\0') {
	genr->err = 1;
	return NULL;
    }

#if GENR_DEBUG
    fprintf(stderr, "fracdiff_series: atom->str = '%s'\n", atom->str);
#endif

    if (sscanf(atom->str, "%8[^,],%lf", vname, param) != 2) {
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

static int add_tmp_series_to_genr (GENERATE *genr, genatom *atom)
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

    if (genr_add_temp_var(genr, y)) {
	free(y);
	return E_ALLOC;
    }

    atom->tmpvar = genr->tmpv - 1;
    atom_eat_children(atom);

    return genr->err;
}

static double genr_get_child_status (GENERATE *genr, genatom *atom)
{
    double val = NADBL;
    genatom *child = pop_child_atom(atom);

    if (child != NULL) {
	if (atom->func == T_VARNUM || atom->func == T_VECTOR) {
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

static int add_statistic_to_genr (GENERATE *genr, genatom *atom)
{
    double val;

    atom_stack_set_parentage(genr);

    if (atom->func == T_VARNUM || atom->func == T_VECTOR ||
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
	    "now atom->scalar = 1\n", atom->func, get_func_word(atom->func),
	    val);
#endif

    atom->val = val;
    atom->scalar = 1;
    atom_eat_children(atom);

    return genr->err;
}

/* ...................................................... */

static int evaluate_genr (GENERATE *genr)
{
    int t, t1 = genr->pdinfo->t1, t2 = genr->pdinfo->t2;
    int m = 0, tstart = t1;
    int n_atoms = 0;
    genatom *atom;

    reset_atom_stack(genr);

    DPRINTF(("Doing evaluate_genr\n"));

    while (!genr->err && (atom = pop_atom(genr))) {
	n_atoms++;

	DPRINTF((" looking at atom %d\n", n_atoms));

	if (atom->varnum == genr->varnum && atom->lag > m) {
	    m = atom->lag;
	}

	if (MODEL_VAR_INDEX(atom->varnum)) {
	    genr->err = add_model_series_to_genr(genr, atom);
	} else if (atom->func == T_UNIFORM || atom->func == T_NORMAL) {
	    genr->err = add_random_series_to_genr(genr, atom);
	} else if (atom->func == T_DIFF || atom->func == T_LDIFF ||
		 atom->func == T_CUM || atom->func == T_SORT ||
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
	}	
    }

    DPRINTF(("evaluate_genr: n_atoms = %d\n", n_atoms));

    if (genr->err) {
	return genr->err;
    }

    if (n_atoms == 2) {
	reset_atom_stack(genr);
	atom = pop_atom(genr);
	if (atom->func == T_SORT) {
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

	    if (0 && y == NADBL) { /* watch out? */
		x = NADBL;
	    } else {
		if (atom->level < level) { 
		    x = calc_pop(genr);
		    npop++;
		    DPRINTF(("popped %g\n", x));
		}
		x = calc_xy(x, y, atom->op, t, &err);
		if (err) {
		    genr_set_warn(genr);
		}
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

static int catch_special_operators (char *s)
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
	} else if (*s == '*' && *(s+1) == '*') {
	    *s = '^';
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
    if (c == '^' || c == '!') 
	return 1;
    if (c == '*' || c == '/' || c == '%') 
	return 2;
    if (c == '+' || c == '-') 
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

static int string_arg_function (const char *s)
{
    if (!strncmp(s, "coeff", 5) ||
	!strncmp(s, "stderr", 6) ||
	!strncmp(s, "rho", 3) ||
	!strncmp(s, "vcv", 3) ||
	!strncmp(s, "corr", 4) ||
	!strncmp(s, "cov", 3) ||
	!strncmp(s, "pvalue", 6) ||
	!strncmp(s, "critical", 8) ||
	!strncmp(s, "fracdiff", 8) ||
	!strncmp(s, "mpow", 4) ||
	!strncmp(s, "mlog", 4)) {
	return 1;
    }

    return 0;
}

static int token_is_atomic (const char *s, GENERATE *genr)
{
    int count = 0;

    DPRINTF(("token_is_atomic: looking at '%s'\n", s));

    /* number in scientific notation */
    if (numeric_string(s)) {
	return 1;
    }

    /* if not numeric string but begins with digit, can't
       be atomic */
    if (isdigit((unsigned char) *s)) {
	return 0;
    }

    /* treat lag variable as atom */
    if (get_lagvar(s, NULL, genr)) {
	return 1;
    }

    while (*s) {
	/* token is non-atomic if it contains an operator,
	   or parentheses */
	if (op_level(*s) || *s == '(') {
	    count++;
	}
	if (count > 0) {
	    break;
	}
	s++;
    }

    return (count == 0);
}

static int token_is_function (char *s, GENERATE *genr, int level)
{
    int wlen = 0;
    int ret;
    const char *p = s;

    while (*p) {
	if (!op_level(*p) && *p != '(') {
	    wlen++; /* might be a problem */
	} else {
	    break;
	}
	p++;
    }

    ret = (*p == '(' && p[strlen(p) - 1] == ')'); 

    if (ret) {
	DPRINTF(("token is function...\n"));
	if (string_arg_function(s)) {
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

static int stack_op_and_token (char *s, GENERATE *genr, int level)
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
	    atom = parse_token(tok, op, genr, level);
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
	genr->err = E_PARSE;
	return genr->err;
    }

    *prev = 0;
    strncat(prev, s, TOKLEN - 1);

    while (*p) {
	DPRINTF(("math_tokenize: inner loop '%s'\n", p));

	if (*p == '(') inparen++;
	else if (*p == ')') inparen--;

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
	
    if (*q) {
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

    if (c == '(') {
	if (lr == 'L') {
	    if (s > 0) {
		s--;
	    }
	} else {
	    s++;
	}
    } else if (c == ')') {
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

#if GENR_DEBUG
static const char *get_func_word (int fnum)
{
    int i;

    for (i=0; funcs[i].fnum != 0; i++) {
	if (fnum == funcs[i].fnum) {
	    return funcs[i].fword;
	}
    }

    return NULL;
}
#endif

int genr_function_from_string (const char *s)
{
    char word[USER_VLEN];
    const char *p;
    int i;

    *word = 0;

    p = strchr(s, '(');
    if (p != NULL && p - s <= 8) {
	strncat(word, s, p - s);
    } else {
	strncat(word, s, 8);
    }

    for (i=0; funcs[i].fnum != 0; i++) {
	if (!strcmp(word, funcs[i].fword)) {
	    return funcs[i].fnum;
	}
    }

    return 0;
}

static void otheruse (const char *s1, const char *s2)
{
    sprintf(gretl_errmsg, _("'%s' refers to a %s and may not be used as a "
			    "variable name"), s1, s2); 
}

/**
 * gretl_reserved_word:
 * @str: string to be tested.
 *
 * Returns 1 if @str is a reserved word that cannot figure as the
 * name of a user-defined variable, otherwise 0.
 */

int gretl_reserved_word (const char *str)
{
    const char *resword[] = {
	"uhat", "yhat",
	"const", "CONST", "pi", "NA",
	"coeff", "stderr", "rho",
	/* stats functions */
	"mean", "median", "var", "cov", "vcv", "sd",
	/* internal sampling vars */
	"full", "subdum", 
	/* plotting vars */
	"t", "annual", "qtrs", "months", "hrs", 
	/* other internal vars */
	"i", "obs", "series",
	NULL
    };
    int i, ret = 0;

    for (i=0; resword[i] != NULL && !ret; i++) {
	if (!strcmp(str, resword[i])) {
	    if (i == 0) {
		otheruse(str, _("residual vector"));
	    } else if (i == 1) {
		otheruse(str, _("fitted values vector"));
	    } else if (i >= 2 && i <= 5) {
		otheruse(str, _("constant"));
	    } else if (i == 6) {
		otheruse(str, _("regr. coeff."));
	    } else if (i == 7) {
		otheruse(str, _("standard error"));
	    } else if (i == 8) {
		otheruse(str, _("autocorr. coeff."));
	    } else if (i >= 9 && i <= 14) {
		otheruse(str, _("stats function"));
	    } else if (i == 15 || i == 16) {
		otheruse(str, _("sampling concept"));
	    } else if (i >= 17 && i <= 21) {
		otheruse(str, _("plotting variable"));
	    } else if (i >= 22 && i <= 24) {
		otheruse(str, _("internal variable"));
	    } else {
		otheruse(str, _("math function"));
	    }
	    ret = 1;
	}
    } 

    if (!ret && genr_function_from_string(str)) {
	otheruse(str, _("math function"));
	ret = 1;
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
    int i, err = 0;

    /* do we have space to make the insertion? */
    if (strlen(s) + llen + 2 >= MAXLINE) {
	err = 1;
    } else {
	memmove(s + llen + 1, s, strlen(s) + 1);
	for (i=0; i<llen; i++) {
	    s[i] = lhs[i];
	}
	s[i] = op;
    }

    return err;
}

static void excise_obs (char *s)
{
    char *p, *q;

    if ((p = strchr(s, '[')) && (q = strchr(p, ']'))) {
	memmove(p, q + 1, strlen(q));
    }
}

static int split_genr_formula (char *lhs, char *s, int obs)
{
    char *p;
    int err = 0;

    if (obs >= 0) {
	excise_obs(s);
    }    

    *lhs = '\0';

    if ((p = strchr(s, '=')) != NULL) {
	char op = 0;

	*p = '\0';

	if (*(p + 1) == '\0') {
	    /* got "equals", but no rhs */
	    err = E_SYNTAX;
	} else {
	    int len = strlen(s);

	    if (len > 1) {
		int l = op_level(s[len - 1]);
		
		if (l == 2 || l == 3) {
		    op = s[len - 1];
		    s[len - 1] = '\0';
		}
	    }

	    /* should we warn if lhs name is truncated? */
	    strncat(lhs, s, USER_VLEN - 1);

	    if (gretl_reserved_word(lhs)) {
		err = 1;
	    } else {
		p++;
		memmove(s, p, strlen(p) + 1);
	    }
	}

	if (!err && op != 0) {
	    err = expand_operator_abbrev(s, lhs, op);
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

/* standardize on '.' for decimal point character with a
   genr formula
*/

static void fix_decimal_commas (char *str)
{
    char *p = str;

    if (p != NULL && *p != 0) {
	p++;
	while (*p && *(p + 1)) {
	    if (*p == ',' && isdigit(*(p - 1)) && isdigit(*(p + 1)))
		*p = '.';
	    p++;
	}
    }
}

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

#define OBS_DEBUG 0

static int plain_obs_number (const char *obs, const DATAINFO *pdinfo)
{
    char *test;
    int t = -1;

#if OBS_DEBUG
    fprintf(stderr, "plain_obs_number: looking at '%s'\n", obs);
#endif

    errno = 0;

    strtol(obs, &test, 10);

    if (*test != '\0' || !strcmp(obs, test) || errno == ERANGE) {
	fprintf(stderr, "plain_obs_number: failed on '%s'\n", obs);
    } else {
	t = atoi(obs) - 1; /* convert from 1-based to 0-based */ 
    } 

    if (t >= pdinfo->n) {
	t = -1;
    }

    return t;
}

static void fix_calendar_date (char *s)
{
    while (*s) {
	if (*s == ':') *s = '/';
	s++;
    }
}

/* Given what looks like an observation number or date within "[" and
   "]", try to determine the observation number.  This is quite tricky
   since we try to handle both dates and plain observation numbers
   (and in addition, variables representing the latter); and we may
   have to deal with the translation from 1-based indexing in user
   space to 0-based indexing for internal purposes.
*/

int get_t_from_obs_string (char *s, const double **Z, 
			   const DATAINFO *pdinfo)
{
    int t;

#if OBS_DEBUG
    fprintf(stderr, "\nget_t_from_obs_string: s ='%s'\n", s);
#endif

    if (calendar_data(pdinfo)) {
	fix_calendar_date(s);
    } 

    t = dateton(s, pdinfo);

#if OBS_DEBUG
    fprintf(stderr, " dateton gives t = %d\n", t);
#endif

    if (t < 0) {
	if (isdigit((unsigned char) *s)) {
	    t = plain_obs_number(s, pdinfo);
#if OBS_DEBUG
	    fprintf(stderr, " plain_obs_number gives t = %d\n", t);
#endif
	} else {
	    int v = varindex(pdinfo, s);

	    if (v < pdinfo->v) {
		t = (int) Z[v][0];
#if OBS_DEBUG
		fprintf(stderr, " based on var %d: t = %d\n", v, t);
#endif
		if (t >= pdinfo->n) {
		    char try[16];

		    sprintf(try, "%d", t);
		    t = dateton(try, pdinfo);
#if OBS_DEBUG
		    fprintf(stderr, " revised via dateton: t = %d\n", t);
#endif
		} else {
		    /* convert to 0-based */
		    t--;
		}
	    }
	}
    }

    if (t < 0) {
	sprintf(gretl_errmsg, _("Observation number out of bounds"));
    }

#if OBS_DEBUG
    fprintf(stderr, " return value: t = %d\n", t);
#endif

    return t;
}

static void get_genr_formula (char *formula, const char *line,
			      GENERATE *genr)
{
    char vname[USER_VLEN], obs[OBSLEN];
    char first[9];

    if (string_is_blank(line)) {
	return;
    }

    /* look at first 'word' in line */
    sscanf(line, "%8s", first);

    if (!strcmp(first, "genr")) {
	line += 4;
    } else if (!strcmp(first, "eval")) {
	genr_unset_save(genr);
	line += 4;
    } else if (!strcmp(first, "series")) {
	genr_force_vector(genr);
	line += 6;
    } else if (!strcmp(first, "scalar")) {
	genr_set_require_scalar(genr);
	line += 6;
    } 

    while (isspace((unsigned char) *line)) {
	line++;
    } 

    /* allow for generating a single value in a series */
    if (sscanf(line, "%8[^[ =][%10[^]]", vname, obs) == 2) {
	genr->obs = get_t_from_obs_string(obs, (const double **) *genr->pZ, 
					  genr->pdinfo);
	if (genr->obs < 0) {
	    genr->err = 1;
	    return;
	}
    }

    *formula = '\0';

    copy_compress(formula, line, MAXLINE - 10);
}

static int gentoler (const char *s)
{
    int err = 0;

    if (numeric_string(s)) {
	double x = dot_atof(s);

	err = set_nls_toler(x);
	if (!err) {
	    sprintf(gretl_msg, _("Set tolerance to %g"), x);
	}
    } else {
	strcpy(gretl_errmsg, _("The setting for \"toler\" must be numeric"));
	err = 1;
    }

    return err;
}

static void 
make_genr_varname (GENERATE *genr, const char *vname)
{
    if (!strncmp(vname, "__", 2)) {
	strcpy(genr->varname, vname + 2);
    } else {
	strcpy(genr->varname, vname);
    }
}

/* substitute something more informative in the genr label
   when the user has called for "$pvalue" or "$test"
*/

static void substitute_in_genrs (char *genrs, char *src)
{
    const char *targ[] = {
	"$pvalue",
	"$test",
	NULL
    };
    int i, slen = 0;
    char *p;

    for (i=0; targ[i] != NULL; i++) {
	if ((p = strstr(genrs, targ[i])) != NULL) {
	    slen = strlen(targ[i]);
	    break;
	}
    }

    if (slen > 0) {
	int srclen = strlen(src);
	
	if (strlen(genrs) + srclen < MAXLINE) {
	    int tail = strlen(p) + 1;

	    *p = ' ';
	    memmove(p + srclen, p, tail);
	    memcpy(p, src, srclen);
	}
    }
}

static void 
write_genr_label (GENERATE *genr, char *genrs, int oldv)
{
    char tmp[64] = {0};
    int llen = 0;

    if (*genr->label != '\0') {
	sprintf(tmp, "%.63s", genr->label);
    }

    if (genr->varnum < oldv) {
	int mc = get_model_count();

	if (mc > 0) {
	    sprintf(genr->label, _("Replaced after model %d: "), mc);
	    llen = 48;
	}
    }	

    if (*tmp != '\0') {
	*genr->label = '\0';
	substitute_in_genrs(genrs, tmp);
    }

    if (strlen(genrs) > MAXLABEL - 1 - llen) {
	strncat(genr->label, genrs, MAXLABEL - 4 - llen);
	strcat(genr->label, "...");
    } else {
	strncat(genr->label, genrs, MAXLABEL - 1);
    }

    strcpy(VARLABEL(genr->pdinfo, genr->varnum), genr->label);    
}

/**
 * generate:
 * @line: command line (formula for generating variable).
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @pmod: pointer to a model, or NULL if not needed.
 * @opt: option flags (for future development).
 *
 * Generates a new variable, usually via some transformation of
 * existing variables, or by retrieving an internal variable associated
 * with the estimation of a model (@pmod).
 * 
 * Returns: 0 on success, integer error code on error.
 */

int generate (const char *line, double ***pZ, DATAINFO *pdinfo, 
	      MODEL *pmod, gretlopt opt)
{
    int i;
    char s[MAXLINE], genrs[MAXLINE];
    char newvar[USER_VLEN];
    int oldv = pdinfo->v;
    GENERATE genr;
#if GENR_DEBUG
    genatom *atom;
#endif

    *gretl_errmsg = *s = *genrs = '\0';
    *gretl_msg = '\0';

    genr_init(&genr, pZ, pdinfo, pmod, opt);

    /* grab the expression, skipping the command word 
       and compressing spaces */
    get_genr_formula(s, line, &genr);

    if (*s == '\0') {
	return E_EQN;
    }

    /* record the full genr expression */
    strcpy(genrs, s);

    DPRINTF(("\n*** starting genr, s='%s'\n", s));

#ifdef ENABLE_NLS
    if (',' == get_local_decpoint())
	fix_decimal_commas(s);
#endif

    /* special cases which are not of the form "lhs=rhs" */
    if (strcmp(s, "dummy") == 0) {
	int di0, orig_v = pdinfo->v;

	di0 = dummy(pZ, pdinfo, 0);
	if (di0 == 0) {
	    genr.err = 1;
	} else if (di0 == orig_v) {
	    strcpy(gretl_msg, _("Periodic dummy variables generated.\n"));
	} else {
	    strcpy(gretl_msg, _("Periodic dummy variables already present.\n"));
	}
	return genr.err;
    } else if (strcmp(s, "paneldum") == 0) {
	genr.err = paneldum(pZ, pdinfo);
	if (!genr.err) {
	    strcpy(gretl_msg, _("Panel dummy variables generated.\n"));
	}
	return genr.err;
    } else if (strcmp(s, "unitdum") == 0) {
	genr.err = panel_unit_dummies(pZ, pdinfo);
	if (!genr.err) {
	    strcpy(gretl_msg, _("Panel dummy variables generated.\n"));
	}
	return genr.err;
    } else if (!strcmp(s, "time") || !strcmp(s, "index")) {
	int tm = !strcmp(s, "time");

	genr.err = genrtime(pZ, pdinfo, tm);
	if (!genr.err) {
	    strcpy(genr.varname, s);
	    genr.varnum = varindex(pdinfo, s);
	    genr_unset_scalar(&genr);
	    compose_genr_msg(&genr, oldv);
	}
	return genr.err;
    } else if (!strcmp(s, "unit")) {
	genr.err = genrunit(pZ, pdinfo);
	if (!genr.err) {
	    strcpy(genr.varname, s);
	    genr.varnum = varindex(pdinfo, s);
	    genr_unset_scalar(&genr);
	    compose_genr_msg(&genr, oldv);
	}
	return genr.err;	
    } else if (strncmp(s, "toler=", 6) == 0) {
	genr.err = gentoler(s + 6);
	return genr.err;
    }

    /* split into lhs = rhs */
    if ((genr.err = split_genr_formula(newvar, s, genr.obs))) {
	return genr.err;
    }
    
    DPRINTF(("after split, newvar='%s', s='%s'\n", newvar, s));

    if (*newvar != '\0') {
	if (strncmp(newvar, "$nl", 3) && 
	    strncmp(newvar, "__", 2) && 
	    check_varname(newvar)) {
	    genr.err = E_SYNTAX;
	    goto genr_return;
	}
	genr.varnum = varindex(pdinfo, newvar);
    } else {
	/* no "lhs=" bit */
	if (!(genr_doing_save(&genr))) {
	    strcpy(newvar, "$eval");
	} else {
	    genr.err = E_SYNTAX;
	    goto genr_return;
	}
    }

    /* process any daily dates in brackets */
    if ((genr.err = fix_obs_in_brackets(s))) {
	return genr.err;
    }

    /* special case of generating a single observation */
    if (genr.obs >= 0) {
	if (genr.varnum >= pdinfo->v) {
	    return E_UNKVAR;
	}
	if (!pdinfo->vector[genr.varnum]) {
	    return E_DATA;
	}
    }

    /* special case of stacking a group of series */
    if (!strncmp(s, "stack(", 6)) {
	genr.err = dataset_stack_variables(pZ, pdinfo, newvar, s);
	return genr.err;
    }

    /* special case of generating observation labels */
    if (!strcmp(newvar, "markers")) {
	genr.err = generate_obs_markers(pZ, pdinfo, s);
	return genr.err;
    }

    /* pre-process special operators */
    if ((genr.err = catch_special_operators(s))) {
	return genr.err;
    }

    DPRINTF(("after catch_special_operators: s='%s'\n", s));

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

    DPRINTF(("after parenthesize: s='%s'\n", s));

    genr.err = attach_atomset(&genr);

    if (!genr.err) {
	genr.err = math_tokenize(s, &genr, 0);
    }

#if GENR_DEBUG
    i = 0;
    while ((atom = pop_atom(&genr))) {
	fprintf(stderr, "*** atom %d ***\n", i++);
	print_atom(atom);
    }
#endif

    if (!genr.err) {
	evaluate_genr(&genr);
    }

    destroy_atom_stack(&genr);
    reset_calc_stack(&genr);

 genr_return:

    if (!genr.err) {
	if (!genr_doing_save(&genr)) {
	    eval_msg(&genr);
	} else {
	    make_genr_varname(&genr, newvar);
	    genr.err = genr_write_var(pZ, pdinfo, &genr);
	    if (!genr.err && !genr_is_private(&genr)) {
		if (!genr_single_obs(&genr)) {
		    write_genr_label(&genr, genrs, oldv);
		}
		compose_genr_msg(&genr, oldv);
	    }
	} 
    }

    genrfree(&genr);

#if GENR_DEBUG
    fprintf(stderr, "\ngenerate: err = %d\n\n", genr.err);
#endif

    return genr.err;
}

static int genr_write_var (double ***pZ, DATAINFO *pdinfo, GENERATE *genr)
{
    double xt = genr->xvec[pdinfo->t1];
    int t, v = genr->varnum;
    int modifying = 0;
    int err = 0;

    /* check that any request for a scalar can be honored,
       abort if not */
    if (genr_require_scalar(genr) && genr_is_vector(genr)) {
	strcpy(gretl_errmsg, _("Specified a scalar, but a vector "
			       "was produced."));
	return 1;
    }	

    /* take into account any forcing the user has attempted */
    if (genr_forcing_vector(genr) && genr_is_scalar(genr)) {
	genr_unset_scalar(genr);
	for (t=0; t<pdinfo->n; t++) {
	    genr->xvec[t] = xt;
	}
    } 

    if (v >= pdinfo->v) {
	/* the generated var is an addition to the data set */
	if (genr_is_vector(genr)) {
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

    /* generally we allow scalar -> vector expansion for existing
       vars, but not if the dataset is subsampled */
    if (modifying && complex_subsampled() &&
	genr_is_vector(genr) && !pdinfo->vector[v]) {
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
		/* we never allow vector -> scalar conversion for
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
	    /* we generated, and will now transcribe, a vector result */
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
    } 

    DPRINTF(("calc_xy: out: x=%g\n", x));

    return x;
}

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

    if (na(arg)) {
	return x;
    }

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

static double hp_lambda (const DATAINFO *pdinfo)
{
    double la = get_hp_lambda();

    if (na(la)) {
	la = 100 * pdinfo->pd * pdinfo->pd;
    }

    return la;
}

/* drop first/last observations from sample if missing obs 
   encountered -- also check for missing vals within the
   remaining sample */

static int array_adjust_t1t2 (const double *x, int *t1, int *t2)
{
    int t, t1min = *t1, t2max = *t2;

    for (t=t1min; t<t2max; t++) {
	if (na(x[t])) t1min++;
	else break;
    }

    for (t=t2max; t>t1min; t--) {
	if (na(x[t])) t2max--;
	else break;
    }

    for (t=t1min; t<=t2max; t++) {
	if (na(x[t])) return t;
    }

    *t1 = t1min; *t2 = t2max;

    return 0;
}

/*
  Hodrick-Prescott filter: adapted from the original FORTRAN code
  by E. Prescott. Very few changes.

  Parameters:
  x: vector of original data
  hp: pointer to a T-vector, returns Hodrick-Prescott "cycle"
*/

static int hp_filter (const double *x, double *hp, const DATAINFO *pdinfo)
{
    int i, t, T, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int err = 0;
    double v00 = 1.0, v11 = 1.0, v01 = 0.0;
    double det, tmp0, tmp1;
    double lambda;

    double **V = NULL;
    double m[2], tmp[2];

    int tb;
    double e0, e1, b00, b01, b11;

    for (t=t1; t<=t2; t++) {
	hp[t] = NADBL;
    }

    err = array_adjust_t1t2(x, &t1, &t2);
    if (err) {
	err = E_DATA;
	goto bailout;
    }

    T = t2 - t1 + 1;
    if (T < 4) {
	err = E_DATA;
	goto bailout;
    }

    lambda = hp_lambda(pdinfo);

    V = malloc(4 * sizeof *V);
    if (V == NULL) return E_ALLOC;

    for (i=0; i<4; i++) {
	V[i] = malloc(T * sizeof **V);
	if (V[i] == NULL) {
	    int j;
	    
	    for (j=0; j<i; j++) {
		free(V[j]);
	    }
	    free(V);
	    return E_ALLOC;
	}
    }

    /* adjust starting points */
    x += t1;
    hp += t1;

    /* covariance matrices for each obs */

    for (t=2; t<T; t++) {
	tmp0 = v00;
	tmp1 = v01;
	v00 = 1.0 / lambda + 4.0 * (tmp0 - tmp1) + v11;
	v01 = 2.0 * tmp0 - tmp1;
	v11 = tmp0;

	det = v00 * v11 - v01 * v01;

	V[0][t] =  v11 / det;
	V[1][t] = -v01 / det;
	V[2][t] =  v00 / det;

	tmp0 = v00 + 1.0;
	tmp1 = v00;
      
	v00 -= v00 * v00 / tmp0;
	v11 -= v01 * v01 / tmp0;
	v01 -= (tmp1 / tmp0) * v01;

    }

    m[0] = x[0];
    m[1] = x[1];

    /* forward pass */
    for (t=2; t<T; t++) {
	tmp[0] = m[1];
	m[1] = 2.0 * m[1] - m[0];
	m[0] = tmp[0];

	V[3][t-1] = V[0][t] * m[1] + V[1][t] * m[0];
	hp[t-1]   = V[1][t] * m[1] + V[2][t] * m[0];
	  
	det = V[0][t] * V[2][t] - V[1][t] * V[1][t];
	  
	v00 =  V[2][t] / det;
	v01 = -V[1][t] / det;
	  
	tmp[1] = (x[t] - m[1]) / (v00 + 1.0);
	m[1] += v00 * tmp[1];
	m[0] += v01 * tmp[1];
    }

    V[3][T-2] = m[0];
    V[3][T-1] = m[1];
    m[0] = x[T-2];
    m[1] = x[T-1];

    /* backward pass */
    for (t=T-3; t>=0; t--) {
	t1 = t+1;
	tb = T - t - 1;
      
	tmp[0] = m[0];
	m[0] = 2.0 * m[0] - m[1];
	m[1] = tmp[0];

	if (t > 1) {
	    /* combine info for y < i with info for y > i */
	    e0 = V[2][tb] * m[1] + V[1][tb] * m[0] + V[3][t];
	    e1 = V[1][tb] * m[1] + V[0][tb] * m[0] + hp[t];
	    b00 = V[2][tb] + V[0][t1];
	    b01 = V[1][tb] + V[1][t1];
	    b11 = V[0][tb] + V[2][t1];
	      
	    det = b00 * b11 - b01 * b01;
	      
	    V[3][t] = (b00 * e1 - b01 * e0) / det;
	}
	  
	det = V[0][tb] * V[2][tb] - V[1][tb] * V[1][tb];
	v00 =  V[2][tb] / det;
	v01 = -V[1][tb] / det;

	tmp[1] = (x[t] - m[0]) / (v00 + 1.0);
	m[1] += v01 * tmp[1];
	m[0] += v00 * tmp[1];
    }

    V[3][0] = m[0];
    V[3][1] = m[1];

    for (t=0; t<T; t++) {
	hp[t] = x[t] - V[3][t];
    }

 bailout:

    for (i=0; i<4; i++) {
	free(V[i]);
    }
    free(V);

    return err;
}

/*
  Baxter & King bandpass filter

  Parameters:
  y: vector of original data
  bk: pointer to a T-vector, returns filtered series
*/

static int bkbp_filter (const double *y, double *bk, const DATAINFO *pdinfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int periods[2];

    double omubar, omlbar;
    double avg_a;
    double *a;

    int i, k, t;
    int err = 0;

    /*
      periods[0] and periods[1]: threshold periodicities for business cycle
      k: order of the approximation
    */

    /* get user settings if available (or the defaults) */
    get_bkbp_periods(periods);
    k = get_bkbp_k();

#if BK_DEBUG
    fprintf(stderr, "lower limit = %d, upper limit = %d, \n", 
	    periods[0], periods[1]);
#endif

    if (periods[0] >= periods[1]) {
	strcpy(gretl_errmsg, "Error in Baxter-King frequencies");
	return 1;
    }

    err = array_adjust_t1t2(y, &t1, &t2);
    if (err) {
	return err;
    } 

    if (2 * k >= t2 - t1 + 1) {
	strcpy(gretl_errmsg, "Insufficient observations");
	return E_DATA;
    }

    a = malloc((k + 1) * sizeof *a);
    if (a == NULL) {
	return E_ALLOC;
    }
    
    omubar = 2.0 * M_PI / periods[0];
    omlbar = 2.0 * M_PI / periods[1];
    
    /* first we compute the coefficients */

    avg_a = a[0] = (omubar - omlbar) / M_PI;

    for (i=1; i<=k; i++) {
	a[i] = (sin(i * omubar) - sin(i * omlbar)) / (i * M_PI);
	avg_a += 2 * a[i];
    }

    avg_a /= (2 * k + 1);

    for (i=0; i<=k; i++) {
	a[i] -= avg_a;
#if BK_DEBUG
	fprintf(stderr, "a[%d] = %#9.6g\n", i, a[i]);
#endif
    }

    /* now we filter the series, skipping the first
       and last k observations */

    for (t=0; t<pdinfo->n; t++) {
	if (t < t1 + k || t >= t2 - k) {
	    bk[t] = NADBL;
	} else {
	    bk[t] = a[0] * y[t];
	    for (i=1; i<=k; i++) {
		bk[t] += a[i] * (y[t-i] + y[t+i]);
	    }
	}
    }

    free(a);

    return err;
}

static int get_fracdiff (const double *y, double *diffvec, double d,
			 const DATAINFO *pdinfo)
{
    int dd, t, T;
    const double TOL = 1.0E-07;
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    double phi = -d;
    int err;

    DPRINTF(("Doing get_fracdiff, with d = %g\n", d));

    err = array_adjust_t1t2(y, &t1, &t2);
    if (err) {
	return E_DATA;
    } 

    T = t2 - t1 + 1;

    for (t=0; t<pdinfo->n; t++) {
	if (t >= t1 && t <= t2) {
	    diffvec[t] = y[t];
	} else {
	    diffvec[t] = NADBL;
	}
    }   

    for (dd=1; dd<=T && fabs(phi)>TOL; dd++) {
	for (t=t1+dd; t<=t2; t++) {
	    diffvec[t] += phi * y[t - dd];
	}
	phi *= (dd - d)/(dd + 1);
    }

    return 0;
}

static double *get_mp_series (const char *s, GENERATE *genr,
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

static double evaluate_statistic (double *z, GENERATE *genr, int fn)
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
    } else if (fn == T_SST) {
	x = gretl_sst(0, i, tmp);
    } else if (fn == T_MEDIAN) {
	x = gretl_median(0, i, tmp);
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
    double x = batch_pvalue(s, Z, pdinfo, NULL);

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
evaluate_bivariate_statistic (const char *s, GENERATE *genr, 
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

static void free_genr_S (GENERATE *genr)
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

static int allocate_genr_S (GENERATE *genr)
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
sort_series (const double *mvec, double *x, GENERATE *genr)
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

    if (markers) {
	qsort(vm, i + 1, sizeof *vm, compare_vms);
    } else {
	qsort(tmp, i + 1, sizeof *tmp, gretl_compare_doubles);
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

static double *get_tmp_series (double *mvec, GENERATE *genr, 
			       int fn, double param)
{
    DATAINFO *pdinfo = genr->pdinfo;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2; 
    double *x;
    double xx, yy;

#if GENR_DEBUG
    fprintf(stderr, "*** Doing get_tmp_series, fn = %d ***\n", fn);
#endif

    x = malloc(pdinfo->n * sizeof *x); 
    if (x == NULL) {
	return NULL;
    }

    if (fn == T_DIFF || fn == T_LDIFF) {
	for (t=t1+1; t<=t2; t++) {
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
		yy = mvec[t-1];
	    }

	    if (na(xx) || na(yy)) {
		x[t] = NADBL;
		continue;
	    }

	    /* perform the differencing */
	    if (fn == T_DIFF) {
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
	x[t1] = NADBL;
    }

    else if (fn == T_CUM) {
	x[t1] = (na(mvec[t1])) ? 0.0 : mvec[t1];
	for (t=t1+1; t<=t2; t++) {
	    if (na(mvec[t])) x[t] = x[t-1];
	    else x[t] = x[t-1] + mvec[t];
	}
    }

    else if (fn == T_SORT) {
	genr->err = sort_series(mvec, x, genr);
	if (genr->err) {
	    free(x);
	    x = NULL;
	}
    }

    else if (fn == T_RESAMPLE) {
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
    }

    else if (fn == T_HPFILT) { 
	genr->err = hp_filter(mvec, x, pdinfo);	
    }

    else if (fn == T_BKFILT) { 
	genr->err = bkbp_filter(mvec, x, pdinfo);	
    }

    else if (fn == T_FRACDIFF) {
	genr->err = get_fracdiff(mvec, x, param, pdinfo);
    }

    return x;
}

static int arma_model_stat_pos (const char *s, const MODEL *pmod)
{
    int p = -1;

    if (numeric_string(s)) {
	p = atoi(s) - 1;
	if (p >= pmod->ncoeff) p = -1;
	return p;
    } else if (pmod->params != NULL) {
	int i;

	for (i=1; i<=pmod->ncoeff; i++) {
	    if (!strcmp(s, pmod->params[i])) {
		p = i - 1;
		break;
	    }
	}
    }

    return p;
}

#define AR1_MODEL(c) (c == CORC || c == HILU || c == PWE)

static double 
get_model_data_element (const char *s, GENERATE *genr,
			MODEL *pmod, int idx)
{
    int lv, vi = 0;
    double x = NADBL;

    DPRINTF(("get_model_data_element: looking at '%s'\n", s));

    if (pmod == NULL) return x;

    if (idx == T_RHO) {
	if (!(numeric_string(s))) {
	    genr->err = E_INVARG;
	} else if (dot_atof(s) == 1 && AR1_MODEL(pmod->ci)) {
	    x = gretl_model_get_double(pmod, "rho_in");
	} else if (pmod->ci != AR && dot_atof(s) == 1) {
	    x = pmod->rho;
	} else if (pmod->arinfo == NULL || 
		   pmod->arinfo->arlist == NULL || 
		   pmod->arinfo->rho == NULL) {
	    genr->err = E_INVARG;
	} else if (!(vi = listpos(atoi(s), pmod->arinfo->arlist))) {
	    genr->err = E_INVARG;
	} else {
	    x = pmod->arinfo->rho[vi-1];
	}
    }

    else if (idx == T_VCV) {
	x = genr_vcv(s, genr->pdinfo, pmod);
	if (na(x)) {
	    genr->err = E_INVARG;
	}
    }

    else if (idx == T_COEFF || idx == T_STDERR) {
	if (pmod == NULL || pmod->list == NULL) {
	    genr->err = E_INVARG;
	} else if (pmod->ci == ARMA) {
	    vi = arma_model_stat_pos(s, pmod);
	    if (vi < 0) {
		genr->err = E_INVARG;
	    }
	} else {
	    lv = numeric_string(s)? atoi(s) : varindex(genr->pdinfo, s);
	    vi = listpos(lv, pmod->list);

	    if (vi < 2) {
		genr->err = E_INVARG;
	    } else {
		vi -= 2;
	    }
	}

	if (!genr->err) {
	    if (idx == T_COEFF && pmod->coeff != NULL) { 
		x = pmod->coeff[vi];
	    } else if (pmod->sderr != NULL) {
		x = pmod->sderr[vi];
	    } else {
		genr->err = E_INVARG;
	    }
	}
    } 

    if (genr->err) {
	gretl_errno = genr->err;
    }

    DPRINTF(("get_model_data_element: err = %d\n", genr->err));

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

static double get_obs_value (const char *s, const double **Z, 
			     const DATAINFO *pdinfo, int *err)
{
    char vname[USER_VLEN], obs[OBSLEN];
    double val = NADBL;

    *err = 0;

    if (sscanf(s, "%8[^[][%10[^]]]", vname, obs) == 2) {
	int t, i = varindex(pdinfo, vname);

	if (i < pdinfo->v && pdinfo->vector[i]) {
	    t = get_t_from_obs_string(obs, Z, pdinfo);
	    if (t >= 0 && t < pdinfo->n) {
		val = Z[i][t];
	    } else {
		*err = 1;
	    }
	} else {
	    *err = 1;
	}
    } else {
	*err = 1;
    }

    return val;
}

static double *get_model_series (const DATAINFO *pdinfo,
				 const MODEL *pmod, int v)
{
    double *x, *garch_h = NULL;
    int t;

    if (pmod == NULL || !MODEL_VAR_INDEX(v)) {
	return NULL;
    }

    if (pmod->t2 - pmod->t1 + 1 > pdinfo->n || 
	model_sample_issue(pmod, NULL, 0, pdinfo)) {
	strcpy(gretl_errmsg, 
	       (v == UHATNUM)? 
	       _("Can't retrieve uhat: data set has changed") :
	       (v == YHATNUM)?
	       _("Can't retrieve yhat: data set has changed") :
	       _("Can't retrieve ht: data set has changed"));
	return NULL;
    }   

    if ((v == UHATNUM && pmod->uhat == NULL) ||
	(v == YHATNUM && pmod->yhat == NULL)) {
	return NULL;
    }

    if (v == HNUM) {
	garch_h = gretl_model_get_data(pmod, "garch_h");
	if (garch_h == NULL) {
	    strcpy(gretl_errmsg, _("Can't retrieve error variance"));
	    return NULL;
	}
    }

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) {
	return NULL;
    }

    for (t=0; t<pdinfo->n; t++) {
	if (t < pmod->t1 || t > pmod->t2) {
	    x[t] = NADBL;
	} else {
	    if (v == UHATNUM) {
		x[t] = pmod->uhat[t];
	    } else if (v == YHATNUM) {
		x[t] = pmod->yhat[t];
	    } else if (v == HNUM) {
		x[t] = garch_h[t];
	    }
	}
    }
	    
    return x;
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

static int obs_num (const char *s, const DATAINFO *pdinfo)
{
    char test[OBSLEN];
    size_t n;
    int t;

    *test = 0;
    strncat(test, (*s == '"')? s + 1 : s, OBSLEN - 1);

    n = strlen(test);
    if (test[n-1] == '"') {
	test[n-1] = '\0';
    }

    if (pdinfo->markers && pdinfo->S != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    if (!strcmp(test, pdinfo->S[t])) return t + 1;
	}
	if (calendar_data(pdinfo)) {
	    charsub(test, ':', '/');
	    for (t=0; t<pdinfo->n; t++) {
		if (!strcmp(test, pdinfo->S[t]) ||
		    !strcmp(test, pdinfo->S[t] + 2)) {
		    return t + 1;
		}
	    }
	}
    }

    if (pdinfo->structure == TIME_SERIES) {
	t = dateton(test, pdinfo);
	if (t >= 0) return t + 1;
    }

    if (calendar_data(pdinfo)) {
	char datestr[OBSLEN];

	charsub(test, ':', '/');
	for (t=0; t<pdinfo->n; t++) {
	    calendar_date_string(datestr, t, pdinfo);
	    if (!strcmp(test, datestr) ||
		!strcmp(test, datestr + 2)) {
		return t + 1;
	    }
	}
    }

    return 0;
}

static int dataset_var_index (const char *s)
{
    char test[USER_VLEN];

    *test = '\0';
    strncat(test, s, USER_VLEN - 1);
    lower(test);

    if (!strcmp(test, "$nobs")) {
	return R_NOBS;
    }
    if (!strcmp(test, "$pd")) {
	return R_PD;
    }
    if (!strcmp(test, "$nvars")) {
	return R_NVARS;
    }

    return 0;
}

static int test_stat_index (const char *s)
{
    char test[USER_VLEN];

    *test = '\0';
    strncat(test, s, USER_VLEN - 1);
    lower(test);

    if (!strcmp(test, "$pvalue"))  
	return R_TEST_PVAL;
    if (!strcmp(test, "$test")) 
	return R_TEST_STAT;

    return 0;
}

static int model_vector_index (const char *s)
{
    char test[USER_VLEN];

    *test = '\0';
    strncat(test, s, USER_VLEN - 1);
    lower(test);

    if (!strcmp(test, "$uhat"))  
	return UHATNUM;
    if (!strcmp(test, "$yhat")) 
	return YHATNUM;
    if (!strcmp(test, "$h"))
	return HNUM;

    return 0;
}

static int panel_x_offset (const DATAINFO *pdinfo, int *bad)
{
    char *p = strchr(pdinfo->stobs, ':');
    int offset = 0;

    if (p == NULL) {
	p = strchr(pdinfo->stobs, '.');
    }

    if (p == NULL) {
	*bad = 1;
    } else {
	offset = atoi(p + 1) - 1;
    }

    return offset;
}

static void 
make_x_panel_dummy (double *x, const DATAINFO *pdinfo, int i)
{
    int t, offset, bad = 0;
    int dmin, dmax;

    offset = panel_x_offset(pdinfo, &bad);

    dmin = (i - 1) * pdinfo->pd;
    dmax = i * pdinfo->pd - offset;

    if (i > 1) dmin -= offset;

    for (t=0; t<pdinfo->n; t++) {
	if (bad) {
	    x[t] = NADBL;
	} else if (t >= dmin && t < dmax) {
	    x[t] = 1.0;
	} else {
	    x[t] = 0.0;
	}
    }
}

static int n_new_dummies (const DATAINFO *pdinfo,
			  int nunits, int nperiods)
{
    char dname[VNAMELEN];
    int i, nnew = nunits + nperiods;

    for (i=0; i<nunits; i++) {
	sprintf(dname, "du_%d", i + 1);
	if (varindex(pdinfo, dname) < pdinfo->v) {
	    nnew--;
	}
    }

    for (i=0; i<nperiods; i++) {
	sprintf(dname, "dt_%d", i + 1);
	if (varindex(pdinfo, dname) < pdinfo->v) {
	    nnew--;
	}
    }

    return nnew;
}

static void 
make_dummy_name_and_label (int vi, const DATAINFO *pdinfo, int center,
			   char *vname, char *vlabel)
{
    if (center > 0) {
	sprintf(vname, "S%d", vi);
	strcpy(vlabel, "centered periodic dummy");
    } else if (center < 0) {
	sprintf(vname, "S%d", vi);
	strcpy(vlabel, "uncentered periodic dummy");
    } else if (pdinfo->pd == 4 && pdinfo->structure == TIME_SERIES) {
	sprintf(vname, "dq%d", vi);
	sprintf(vlabel, 
		_("= 1 if quarter = %d, 0 otherwise"), vi);
    } else if (pdinfo->pd == 12 && pdinfo->structure == TIME_SERIES) {
	char mname[8];

	get_month_name(mname, vi);
	sprintf(vname, "d%s", mname);
	sprintf(vlabel, _("= 1 if month is %s, 0 otherwise"), mname);
    } else {
	char dumstr[8] = "dummy_";
	char numstr[8];
	int len;

	sprintf(numstr, "%d", vi);
	len = strlen(numstr);
	dumstr[8 - len] = '\0';
	sprintf(vname, "%s%d", dumstr, vi);
	sprintf(vlabel, _("%s = 1 if period is %d, 0 otherwise"), vname, vi);
    }
}

/**
 * dummy:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @center: if greater than zero subtract the population mean from
 * each of the generated dummies; if less than zero, do not
 * subtract the mean but generate dummies with labels on the
 * same pattern as centered dummies (for internal use in VECMs).
 * Usually this argument is set to zero.
 *
 * Adds to the data set (if these variables are not already 
 * present) a set of periodic (usually seasonal) dummy variables.
 *
 * Returns: the ID number of the first dummy variable on success,
 * or 0 on error.
 */

int dummy (double ***pZ, DATAINFO *pdinfo, int center)
{
    char vname[USER_VLEN];
    char vlabel[MAXLABEL];
    int vi, t, yy, pp, mm;
    int ndums, nnew = 0;
    int di, di0 = pdinfo->v;
    double xx, dx;

    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	ndums = pdinfo->n / pdinfo->pd;
	if (pdinfo->n % pdinfo->pd) {
	    ndums++;
	}
    } else {
	ndums = pdinfo->pd;
    }

    if (ndums == 1 || ndums > 99999) {
	strcpy(gretl_errmsg, _("This command won't work with the current periodicity"));
	return 0;
    }

    for (vi=0; vi<ndums; vi++) {
	make_dummy_name_and_label(vi + 1, pdinfo, center, vname, vlabel);
	di = varindex(pdinfo, vname);
	if (di >= pdinfo->v || strcmp(vlabel, VARLABEL(pdinfo, di))) {
	    nnew++;
	} else if (vi == 0) {
	    di0 = di;
	} else if (di != di0 + vi) {
	    /* dummies not consecutive: problem */
	    di0 = pdinfo->v;
	    nnew = ndums;
	    break;
	}
    }

    if (nnew == 0) {
	/* all dummies already present */
	return di0;
    } else if (pZ == NULL) {
	return -1;
    }

    if (dataset_add_series(ndums, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Out of memory error"));
	return 0;
    }

    pp = pdinfo->pd;
    mm = 10;
    while ((pp = pp / 10)) {
	mm *= 10;
    }

    for (vi=1, di = di0; vi<=ndums; vi++, di++) {
	make_dummy_name_and_label(vi, pdinfo, center, vname, vlabel);
	strcpy(pdinfo->varname[di], vname);
	strcpy(VARLABEL(pdinfo, di), vlabel);

	if (pdinfo->structure == STACKED_CROSS_SECTION) {
	    make_x_panel_dummy((*pZ)[di], pdinfo, vi);
	} else {
	    for (t=0; t<pdinfo->n; t++) {
		xx = date(t, pdinfo->pd, pdinfo->sd0);
		if (dataset_is_daily(pdinfo)) { /* FIXME weekly? */
		    xx += .1;
		}
		yy = (int) xx;
		pp = (int) (mm * (xx - yy) + 0.5);
		dx = (pp == vi)? 1.0 : 0.0;
		(*pZ)[di][t] = dx;
	    }
	}
    }

    if (center > 0) {
	double cx = 1.0 / pdinfo->pd;
	int vimax = di0 + pdinfo->pd - 1;

	for (vi=di0; vi<=vimax; vi++) {
	    for (t=0; t<pdinfo->n; t++) {
		(*pZ)[vi][t] -= cx;
	    }
	}	
    }

    return di0;
}

/* if both == 1, generate both unit and period dummies,
   otherwise just generate unit dummies
*/

static int real_paneldum (double ***pZ, DATAINFO *pdinfo,
			  int both)
{
    char vname[16];
    int vi, t, yy, pp, mm;
    int xsect, orig_v = pdinfo->v;
    int ndum, nnew, n_blockdum = 0, n_freqdum = 0;
    int newvnum, offset, bad = 0;
    double xx;

    xsect = (pdinfo->structure == STACKED_CROSS_SECTION);

    /* in case xsect, block dummies are per time-period,
       frequency dummies are for units;
       in case of stacked time series, block dummies are
       per cross-sectional unit, frequency ones are for periods
    */

    if (both || xsect) {
	n_freqdum = pdinfo->pd;
	if (n_freqdum == 1) {
	    return E_PDWRONG;
	}
    }

    if (both || !xsect) {
	n_blockdum = pdinfo->n / pdinfo->pd;
	if (pdinfo->n % pdinfo->pd) {
	    n_blockdum++;
	}
	if (n_blockdum == 1) {
	    return E_PDWRONG;
	}
    }

    ndum = n_freqdum + n_blockdum;

    nnew = n_new_dummies(pdinfo, 
			 (xsect)? n_freqdum : n_blockdum,
			 (xsect)? n_blockdum : n_freqdum);

    if (dataset_add_series(nnew, pZ, pdinfo)) {
	return E_ALLOC;
    }

    pp = pdinfo->pd;
    mm = 10;
    while ((pp = pp / 10)) {
	mm *= 10;
    }

    newvnum = orig_v;

    /* first generate the frequency-based dummies */
    for (vi=1; vi<=n_freqdum; vi++) {
	int dnum;

	if (xsect) {
	    sprintf(vname, "du_%d", vi);
	} else {
	    sprintf(vname, "dt_%d", vi);
	}

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		(xsect)? _("unit"): _("period"), vi);

	for (t=0; t<pdinfo->n; t++) {
	    xx = date(t, pdinfo->pd, pdinfo->sd0);
	    yy = (int) xx;
	    pp = (int) (mm * (xx - yy) + 0.5);
	    (*pZ)[dnum][t] = (pp == vi)? 1.0 : 0.0;
	}
    }

    offset = panel_x_offset(pdinfo, &bad);

    /* and then the block-based ones */
    for (vi=1; vi<=n_blockdum; vi++) {
	int dmin = (vi-1) * pdinfo->pd;
	int dmax = vi * pdinfo->pd - offset;
	int dnum;

	if (vi > 1) dmin -= offset;

	if (xsect) {
	    sprintf(vname, "dt_%d", vi);
	} else {
	    sprintf(vname, "du_%d", vi);
	}

	dnum = varindex(pdinfo, vname);
	if (dnum >= orig_v) {
	    dnum = newvnum++;
	}	

	strcpy(pdinfo->varname[dnum], vname);
	sprintf(VARLABEL(pdinfo, dnum), 
		_("%s = 1 if %s is %d, 0 otherwise"), vname, 
		(xsect)? _("period"): _("unit"), vi);

	for (t=0; t<pdinfo->n; t++) {
	    if (bad) {
		(*pZ)[dnum][t] = NADBL;
	    } else if (t >= dmin && t < dmax) {
		(*pZ)[dnum][t] = 1.0;
	    } else {
		(*pZ)[dnum][t] = 0.0;
	    }
	}
    }

    return 0;
}

/**
 * panel_unit_dummies:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Adds to the data set a set of dummy variables corresponding
 * to the cross-sectional units in a panel.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int panel_unit_dummies (double ***pZ, DATAINFO *pdinfo)
{
    return real_paneldum(pZ, pdinfo, 0);
}

/**
 * paneldum:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Adds to the data set a set of panel data dummy variables (for
 * both unit and period).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int paneldum (double ***pZ, DATAINFO *pdinfo)
{
    return real_paneldum(pZ, pdinfo, 1);
}

/**
 * genrunit:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * (For panel data only) adds to the data set an index variable 
 * that uniquely identifies the cross-sectional units.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int genrunit (double ***pZ, DATAINFO *pdinfo)
{
    int xt = 0;
    int i, t;

    if (pdinfo->structure != STACKED_TIME_SERIES &&
	pdinfo->structure != STACKED_CROSS_SECTION) {
	strcpy(gretl_errmsg, "'genr unit' can be used only with "
	       "panel data");
	return 1;
    }

    i = varindex(pdinfo, "unit");

    if (i == pdinfo->v) {
	if (dataset_add_series(1, pZ, pdinfo)) return E_ALLOC;
    }

    strcpy(pdinfo->varname[i], "unit");
    strcpy(VARLABEL(pdinfo, i), _("cross-sectional unit index"));

    if (pdinfo->structure == STACKED_CROSS_SECTION) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt = 1;
	    }
	    (*pZ)[i][t] = (double) xt++;
	}
    } else {
	/* stacked time series */
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt++;
	    }
	    (*pZ)[i][t] = (double) xt;
	}
    }    

    return 0;
}

/* make special time variable for panel data */

static void 
make_panel_time_var (double *x, const DATAINFO *pdinfo)
{
    int t, xt = 0;

    if (pdinfo->structure == STACKED_TIME_SERIES) {
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt = 1;
	    }
	    x[t] = (double) xt++;
	}
    } else {
	/* stacked cross-sections */
	for (t=0; t<pdinfo->n; t++) {
	    if (t % pdinfo->pd == 0) {
		xt++;
	    }
	    x[t] = (double) xt;
	}
    }
}

/* create obs index or time trend variable */

int genrtime (double ***pZ, DATAINFO *pdinfo, int tm)
{
    int i, t;

    i = varindex(pdinfo, (tm)? "time" : "index");

    if (i == pdinfo->v && dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    if (tm) {
	strcpy(pdinfo->varname[i], "time");
	strcpy(VARLABEL(pdinfo, i), _("time trend variable"));
    } else {
	strcpy(pdinfo->varname[i], "index");
	strcpy(VARLABEL(pdinfo, i), _("data index variable"));
    }
    
    if (tm && 
	(pdinfo->structure == STACKED_TIME_SERIES ||
	 pdinfo->structure == STACKED_CROSS_SECTION)) {
	make_panel_time_var((*pZ)[i], pdinfo);
    } else {
	for (t=0; t<pdinfo->n; t++) {
	    (*pZ)[i][t] = (double) (t + 1);
	}
    }

    return 0;
}

static int plotvar_is_full_size (int v, int n, const double *x)
{
    int t, ret = 1;

    for (t=0; t<n; t++) {
	if (na(x[t])) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * plotvar:
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @period: string to identify periodicity: "annual", "qtrs",
 * "months", "decdate" (calendar data), "time" or "index".
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

    if (vi < pdinfo->v) {
	if (plotvar_is_full_size(vi, pdinfo->n, (*pZ)[vi])) {
	    return vi;
	} 
    } else if (dataset_add_series(1, pZ, pdinfo)) {
	return -1;
    }

    strcpy(pdinfo->varname[vi], period);

    y1 = (int) pdinfo->sd0;
    rm = pdinfo->sd0 - y1;

    switch (period[0]) {
    case 'a':
	strcpy(VARLABEL(pdinfo, vi), _("annual plotting variable")); 
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + atoi(pdinfo->stobs));
	}
	break;
    case 'q':
	strcpy(VARLABEL(pdinfo, vi), _("quarterly plotting variable"));
	(*pZ)[vi][0] = y1 + (10.0 * rm - 1.0) / 4.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + .25;
	}
	break;
    case 'm':
	strcpy(VARLABEL(pdinfo, vi), _("monthly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0) / 12.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0 / 12.0);
	}
	break;
    case 'h':
	strcpy(VARLABEL(pdinfo, vi), _("hourly plotting variable"));
	(*pZ)[vi][0] = y1 + (100.0 * rm - 1.0) / 24.0;
	for (t=1; t<n; t++) {
	    (*pZ)[vi][t] = (*pZ)[vi][t-1] + (1.0 / 24.0);
	}
	break;
    case 'd':
	if ((dated_daily_data(pdinfo) && pdinfo->n > 365) ||
	    (dated_weekly_data(pdinfo) && pdinfo->n > 52)) {
	    strcpy(VARLABEL(pdinfo, vi), _("daily plotting variable"));
	    for (t=0; t<n; t++) {
		if (pdinfo->S != NULL) {
		    (*pZ)[vi][t] = get_dec_date(pdinfo->S[t]);
		} else {
		    char datestr[OBSLEN];
		    
		    calendar_date_string(datestr, t, pdinfo);
		    (*pZ)[vi][t] = get_dec_date(datestr);
		}
	    } 
	} else if (dataset_is_decennial(pdinfo)) {
	    strcpy(pdinfo->varname[vi], "time");
	    strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	    for (t=0; t<n; t++) {
		(*pZ)[vi][t] = pdinfo->sd0 + 10 * t;
	    }	    
	} else {
	    strcpy(pdinfo->varname[vi], "time");
	    strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	    for (t=0; t<n; t++) {
		(*pZ)[vi][t] = (double) (t + 1);
	    }
	}
	break; 
    case 'i':
	strcpy(VARLABEL(pdinfo, vi), _("index variable"));
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + 1);
	}
	break;
    case 't':
	strcpy(VARLABEL(pdinfo, vi), _("time trend variable"));
	for (t=0; t<n; t++) {
	    (*pZ)[vi][t] = (double) (t + 1);
	}
	break;
    default:
	break;
    }

    return vi;
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
    int level = 0;
    int i, j, n;

    if (gretl_executing_function()) {
	level = gretl_function_stack_depth();

	n = 0;
	for (i=0; i<pdinfo->v; i++) {
	    if (STACK_LEVEL(pdinfo, i) == level) {
		n++;
	    }
	}
    } else {
	n = pdinfo->v;
    }

    pprintf(prn, _("Listing %d variables:\n"), n);

    j = 1;
    for (i=0; i<pdinfo->v; i++) {
	if (level > 0 && STACK_LEVEL(pdinfo, i) != level) {
	    continue;
	}
	pprintf(prn, "%3d) %-10s", i, pdinfo->varname[i]);
	if (j % 5 == 0) {
	    pputc(prn, '\n');
	}
	j++;
    }

    if (n % 5) pputc(prn, '\n');

    pputc(prn, '\n');
}

/**
 * maybe_list_vars:
 * @pdinfo: data information struct.
 * @prn: gretl printing struct
 *
 * Prints a list of the names of the variables currently defined,
 * unless gretl messaging is turned off.
 */

void maybe_list_vars (const DATAINFO *pdinfo, PRN *prn)
{
    if (gretl_messages_on()) {
	varlist(pdinfo, prn);
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

    while (*check == '_') {
	check++;
    }

    if (!strcmp(check, "uhat") || !strcmp(check, "$uhat")) {
	return UHATNUM; 
    } else if (!strcmp(check, "yhat") || !strcmp(check, "$yhat")) {
	return YHATNUM; 
    } else if (!strcmp(check, "$h")) {
	return HNUM; 
    } else if (!strcmp(check, "t") || !strcmp(check, "obs")) {
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

static void genrfree (GENERATE *genr)
{
    int i;

    if (genr == NULL) return;

    DPRINTF(("genrfree: freeing %d vars\n", genr->tmpv));

    if (genr->tmpv > 0) {
	for (i=0; i<genr->tmpv; i++) {
	    free(genr->tmpZ[i]);
	}
	free(genr->tmpZ);
    }

    free(genr->xvec);
}

static int genr_mpow (const char *str, double *xvec, double **Z, 
		      DATAINFO *pdinfo)
{
    int err, v;
    unsigned pwr;
    char vname[USER_VLEN];
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
    char vname[USER_VLEN];
    void *handle = NULL;
    int (*mp_log) (const double *, double *, int);
    
    if (sscanf(str, "%8s", vname) != 1) {
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

static int varnames_from_arg (const char *s, char *v1str, char *v2str)
{
    int i, p, n = strlen(s);

    if (n > 17) {
	return 1;
    }

    p = haschar(',', s);
    if (p < 0 || p > 8) {
	return 1;
    }

    /* get first var name */
    for (i=0; i<p; i++) {
	v1str[i] = s[i];
    }
    v1str[p] = '\0';

    /* get second var name */
    n = n - p - 1;
    for (i=0; i<n; i++) {
	v2str[i] = s[p+1+i];
    }
    v2str[i] = '\0';

    return 0;
}

static double genr_cov_corr (const char *str, double ***pZ, 
			     const DATAINFO *pdinfo, int fn)
{
    char v1str[USER_VLEN], v2str[USER_VLEN];
    int v1, v2;
    double ret = NADBL;

    if (varnames_from_arg(str, v1str, v2str)) {
	return NADBL;
    }

    v1 = varindex(pdinfo, v1str);
    v2 = varindex(pdinfo, v2str);

    if (v1 >= pdinfo->v || v2 >= pdinfo->v) {
	return NADBL;
    }

    if (fn == T_COV) {
	ret = gretl_covar(pdinfo->t1, pdinfo->t2, (*pZ)[v1], (*pZ)[v2]);
    } else if (fn == T_CORR) {
	ret = gretl_corr(pdinfo->t1, pdinfo->t2, (*pZ)[v1], (*pZ)[v2], NULL);
    }

    return ret;
}

static int get_model_param_number (const MODEL *pmod, 
				   const char *vname)
{
    int i, ret = 0;

    if (pmod->params == NULL) {
	return 0;
    }

    for (i=0; i<=pmod->ncoeff; i++) {
	if (!strcmp(vname, pmod->params[i])) {
	    ret = i + 1;
	    break;
	}
    }

    return ret;
}

static double genr_vcv (const char *str, const DATAINFO *pdinfo, 
			MODEL *pmod)
{
    int v1 = 0, v2 = 0;
    int i, j, k, v1l, v2l;
    char v1str[USER_VLEN], v2str[USER_VLEN];
    int gotit;
    double ret = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	return NADBL;
    }

    if (varnames_from_arg(str, v1str, v2str)) {
	return NADBL;
    }

    /* are the varnames valid? */
    if (pmod->ci != NLS && pmod->ci != ARMA) {
	v1 = varindex(pdinfo, v1str);
	v2 = varindex(pdinfo, v2str);
	if (v1 >= pdinfo->v || v2 >= pdinfo->v) {
	    return NADBL;
	}
    }

    /* check model list */
    if (pmod->ci == NLS || pmod->ci == ARMA) {
	v1l = get_model_param_number(pmod, v1str);
	v2l = get_model_param_number(pmod, v2str);
    } else {
	v1l = listpos(v1, pmod->list);
	v2l = listpos(v2, pmod->list);
    }

    if (v1l == 0 || v2l == 0) {
	return NADBL;
    }

    v1l -= 2;
    v2l -= 2;

    /* make model vcv matrix if need be */
    if (pmod->vcv == NULL && makevcv(pmod)) {
	return NADBL;
    }

    /* now find the right entry */
    if (v1l > v2l) {
	k = v1l;
	v1l = v2l;
	v2l = k;
    }

    gotit = 0;
    k = 0;
    for (i=0; i<pmod->ncoeff && !gotit; i++) {
	for (j=0; j<pmod->ncoeff; j++) {
	    if (j < i) {
		continue;
	    }
	    if (i == v1l && j == v2l) {
		ret = pmod->vcv[k];
		gotit = 1;
		break;
	    }
	    k++;
	}
    }

    return ret;
}

static void eval_msg (const GENERATE *genr)
{
    double x = genr->xvec[genr->pdinfo->t1];

    if (na(x)) {
	strcpy(gretl_msg, " NA");
    } else {
	sprintf(gretl_msg, " %g", x);
    }
}

static void compose_genr_msg (const GENERATE *genr, int oldv)
{
    double x;
    int scalar = genr_is_scalar(genr);
    int mutant = 0;

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

    sprintf(gretl_msg, "%s %s %s (ID %d)", 
	    (genr->obs >= 0)? _("Modified") :
	    (genr->varnum < oldv)? _("Replaced") : _("Generated"), 
	    (mutant)? _("variable") :
	    (scalar)? _("scalar") : _("vector"),
	    genr->varname, genr->varnum);

    if (scalar) {
	char numstr[24];

	x = genr->xvec[genr->pdinfo->t1];
	if (na(x)) {
	    strcpy(numstr, " = NA");
	} else {
	    sprintf(numstr, " = %g", x);
	}
	strcat(gretl_msg, numstr);
    }

    if (genr_warn(genr)) {
	strcat(gretl_msg, "\n");
	strcat(gretl_msg, gretl_errmsg);
	*gretl_errmsg = '\0';
    }
}

static int listpos (int v, const int *list)
{
    int i, lmax = list[0];

    /* handle special TSLS list */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    lmax = i - 1;
	    break;
	}
    }
	    
    for (i=lmax; i>=1; i--) {
	if (v == list[i]) return i;
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

int genr_fit_resid (const MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		    int code, int undo)
{
    char vname[USER_VLEN], vlabel[MAXLABEL];
    int i, t;
    double *h = NULL;

    if (code == GENR_H) {
	h = gretl_model_get_data(pmod, "garch_h");
	if (h == NULL) return E_DATA;
    }

    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }

    i = pdinfo->v - 1;

    for (t=0; t<pdinfo->n; t++) {
	(*pZ)[i][t] = NADBL;
    }

    if (code == GENR_RESID) {
	sprintf(vname, "uhat%d", pmod->ID);
	sprintf(vlabel, _("residual from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = pmod->uhat[t];
	}
    } else if (code == GENR_FITTED) {
	sprintf(vname, "yhat%d", pmod->ID);
	sprintf(vlabel, _("fitted value from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = pmod->yhat[t];
	}
    } else if (code == GENR_RESID2) { 
	/* squared residuals */
	sprintf(vname, "usq%d", pmod->ID);
	sprintf(vlabel, _("squared residual from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (na(pmod->uhat[t])) {
		(*pZ)[i][t] = NADBL;
	    } else {
		(*pZ)[i][t] = pmod->uhat[t] * pmod->uhat[t];
	    }
	}
    } else if (code == GENR_H) { 
	/* garch variance */
	sprintf(vname, "h%d", pmod->ID);
	sprintf(vlabel, _("fitted variance from model %d"), pmod->ID);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    (*pZ)[i][t] = h[t];
	}
    }

    strcpy(pdinfo->varname[i], vname);

    if (!undo) {
	strcpy(VARLABEL(pdinfo, i), vlabel);
    }

    return 0;
}


    
