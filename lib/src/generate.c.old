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

typedef struct _GENERATE GENERATE;

struct _GENERATE {
    double *xvec;
    int varnum;
    char varname[VNAMELEN];
    char label[MAXLABEL];
    int scalar; 
    double scalval;
};

static int cstack (double *xstack, double *xvec, char op, 
		   const DATAINFO *pdinfo, int expand);
static int domath (double *xvec, const double *mvec, int nt, 
		   const DATAINFO *pdinfo, int *scalar);
static int evalexp (char *ss, int nt, double *mvec, double *xvec, 
		    double **Z, const DATAINFO *pdinfo, 
		    const MODEL *pmod, GENERATE *genr);
static void getvar (char *str, char *word, char *c);
static int getxvec (char *s, double *xvec, 
		    double **Z, const DATAINFO *pdinfo, 
		    const MODEL *pmod, int *scalar);
static int scanb (const char *ss, char *word);
static int strtype (char *ss, const DATAINFO *pdinfo);
static int whichtrans (const char *ss);
static int createvar (double *xvec, char *snew, char *sleft, 
		      char *sright, int ssnum, double ***pZ, 
		      DATAINFO *pdinfo, int scalar);
static void genrfree (double ***pZ, DATAINFO *pdinfo, GENERATE *pgenr,
		      double *mstack, double *mvec, int nv);
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
static int obs_num (const char *ss, const DATAINFO *pdinfo);
static int has_lhs_lags (const char *str, const char *vname);
static char *genr_to_sim (const char *s, int nlags, const char *vname,
			  const DATAINFO *pdinfo);

#ifdef HAVE_MPFR
static int genr_mlog (const char *str, double *xvec, double **pZ, 
		      DATAINFO *pdinfo);
#endif

static void genr_msg (GENERATE *pgenr, int nv);
static int ismatch (int lv, const int *list);
static void varerror (const char *ss);
static int genrtime (double ***pZ, DATAINFO *pdinfo, GENERATE *genr, int time);
static int add_new_var (DATAINFO *pdinfo, double ***pZ, GENERATE *genr);

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
    R_PD,
    R_NUMERIC,
    R_MATH,
    R_VARNAME,
    R_VAROBS,
    R_OBSNUM,
    R_UNKNOWN
};

enum composites {
    NEQ = 21,
    GEQ,
    LEQ
};

enum set_or_get {
    OBS_SET,
    OBS_GET
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

static char operators[] = {
    '+', '-', '|',
    '*', '/', '%', '&',
    '^', '<', '>', '=', '!', NEQ, GEQ, LEQ, 0
};

#define LEVELS 7

#define SCALAR_SCOPE(t) (t == T_MEAN || t == T_SD || t == T_SUM || \
                         t == T_CORR || t == T_COV || \
                         t == T_VAR || t == T_MEDIAN || t == T_MIN || \
                         t == T_SST || t == T_MAX || t == T_NOBS)

#define MAXTERMS 64
#define _VSMALL 1.0e-14

/* ...................................................... */

static int is_operator (char c)
{
    int i;

    for (i=0; operators[i] != 0; i++) 
	if (c == operators[i]) return 1;
    return 0;
}

/* ...................................................... */

static void catch_double_symbols (char *str)
{
    int i, n = strlen(str);

    for (i=1; i<n; i++) {
	if (str[i] == '=') {
	    if (str[i-1] == '!') {
		str[i-1] = NEQ;
		str[i] = ' ';
	    }
	    else if (str[i-1] == '>') {
		str[i-1] = GEQ;	
		str[i] = ' ';
	    }
	    else if (str[i-1] == '<') {
		str[i-1] = LEQ;	
		str[i] = ' ';
	    }	    
	}
	/* recognize "**" as well as "^" for power */
	else if (str[i] == '*') {
	    if (str[i-1] == '*') {
		str[i-1] = '^';
		str[i] = ' ';
	    }
	}
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

static int insert_paren (char *s, int pos, char lr)
{
    int i, n = strlen(s);

    if (n + 1 >= MAXLEN) return 1;
    for (i=n+1; i>=pos+1; i--) s[i] = s[i - 1];
    if (lr == 'L') s[pos + 1] = '(';
    else s[pos + 1] = ')';

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
    int i, k, oppos, n = strlen(str);
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
    register int t;

    for (t=0; t<n; t++) 
	if (floatneq(x[t], y[t])) 
	    return 0;
    return 1;
}

/* ........................................................  */

static void otheruse (const char *str1, const char *str2)
{
    sprintf(gretl_errmsg, _("'%s' refers to a %s and may not be used as a "
			    "variable name"), str1, str2); 
}

/* .......................................................... */

int get_function (const char *s)
{
    int i;

    for (i=0; ; i++) {
	if (math[i] == NULL) break;
        if (strcmp(s, math[i]) == 0) return 1;
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

    if (get_function(str)) {
	otheruse(str, _("math function"));
	return 1;
    }
 
    return 0;
}

/* .......................................................... */

static const char *set_or_get_obs_marker (const char *s, int opt)
{
    static char obsstr[OBSLEN];

    if (opt == OBS_SET) {
	size_t n;

	*obsstr = 0;
	strncat(obsstr, s + 1, 8);
	n = strlen(obsstr);
	if (obsstr[n - 1] == '"') obsstr[n - 1] = 0;
#ifdef GENR_DEBUG
	fprintf(stderr, "set_or_get_obs_marker: obsstr='%s'\n", obsstr);
#endif
	return NULL;
    } else {
	return obsstr;
    }
}

static void make_obs_dummy (double *x, const DATAINFO *pdinfo)
{
    const char *obs;
    int t, gotit = 0;

    if (pdinfo->S == NULL) return;

    obs = set_or_get_obs_marker(NULL, OBS_GET);

    for (t=0; t<pdinfo->n; t++) {
	if (gotit) x[t] = 0.0;
	else {
	    if (!strcmp(pdinfo->S[t], obs)) {
		x[t] = t + 1.0;
		gotit = 1;
	    } else {
		x[t] = 0.0;
	    }
	}
#ifdef GENR_DEBUG
	fprintf(stderr, "obs_dummy: x[%d] = %g\n", t, x[t]);
#endif
    }
}

/* .......................................................... */

static void copy (const char *str, int indx, 
		  int count, char *dest)
     /* copies count chars from indx in str to dest */
{
    int i;

    *dest = '\0';
    for (i=0; i<count; ++i) dest[i] = str[indx + i];
    dest[count] = '\0';
}

/* .........................................................    */

static int getword (char c, char *str, char *word, unsigned char oflag)

     /* Scans string str for char c, gets word to the left of it as
	"word" and deletes word from str.
	Returns number of chars deleted, or -1 if no occurrence of c, 
	or 0 if reserved word is used 
     */
{
    int i;

    *word = '\0';

    i = haschar(c, str);
    if (i == -1) return -1;

    copy(str, 0, i, word);
    _delete(str, 0, i + 1);

    /* special case for auto sub-sampling dummy */
    if (oflag && strcmp(word, "subdum") == 0)
	return i + 1;

    if (_reserved(word)) return 0;

    return i + 1;
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

static int handle_type2 (int type2, char *word, char *sexpr, 
			 GENERATE *genr, double *mvec,
			 double ***pZ, DATAINFO *pdinfo, 
			 MODEL *pmod) 
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2, n = pdinfo->n;
    int i, vi, nt, lv;
    double xx;
    int err = 0;

    switch (type2) {

    case R_VARNAME:    
	if (!(_isnumber(sexpr)))  {
	    err = E_NOTINTG;
	} else {
	    vi = varindex(pdinfo, word);
	    if (!pdinfo->vector[vi]) {
		sprintf(gretl_errmsg, _("Variable %s is a scalar; "
					"can't do lags/leads"), 
			pdinfo->varname[vi]);
		err = 1;
	    } else {
		genr->scalar = 0;
		get_lag(vi, -atoi(sexpr), mvec, *pZ, pdinfo);
		for (i=t1; i<=t2; i++) {
		    genr->xvec[i] = mvec[i];
		}
	    }
	}
	break;

    case R_MATH:  
	nt = whichtrans(word);
#ifdef GENR_DEBUG
	fprintf(stderr, "R_MATH: nt=%d\n", nt);
#endif  
	if (nt == T_RHO) {
	    if (!(_isnumber(sexpr))) {
		return E_INVARG;
	    }
	    if (dot_atof(sexpr) == 1 && (pmod->ci == CORC ||
					 pmod->ci == HILU)) {
		for (i=t1; i<=t2; i++)
		    genr->xvec[i] = pmod->rho_in;
	    }
	    else if (pmod->ci != AR && dot_atof(sexpr) == 1) {
		for (i=t1; i<=t2; i++)
		    genr->xvec[i] = pmod->rho;
	    }
	    else if (pmod->arinfo == NULL || 
		pmod->arinfo->arlist == NULL || pmod->arinfo->rho == NULL) {
		err = E_INVARG;
	    }
	    else if (!(vi = ismatch(atoi(sexpr), pmod->arinfo->arlist))) {
		err = E_INVARG;
	    }
	    else {
		for (i=0; i<n; i++) 
		    genr->xvec[i] = pmod->arinfo->rho[vi];
	    }
	}
	else if (nt == T_NORMAL) {
	    genr->scalar = 0;
	    gretl_normal_dist(genr->xvec, t1, t2);
	}   
	else if (nt == T_UNIFORM) {
	    genr->scalar = 0;
	    gretl_uniform_dist(genr->xvec, t1, t2);
	}
	else if (nt == T_COV) {
	    xx = genr_cov(sexpr, pZ, pdinfo);
	    if (na(xx)) {
		err = E_INVARG;
	    } else {
		for (i=0; i<n; i++) 
		    genr->xvec[i] = xx;
	    }
	}
	else if (nt == T_CORR) {
	    xx = genr_corr(sexpr, pZ, pdinfo);
	    if (na(xx)) {
		err = E_INVARG;
	    } else {
		for (i=0; i<n; i++) 
		    genr->xvec[i] = xx;
	    }
	}
	else if (nt == T_VCV) {
	    xx = genr_vcv(sexpr, pdinfo, pmod);
	    if (na(xx)) {
		err = E_INVARG;
	    } else {
		for (i=0; i<n; i++) 
		    genr->xvec[i] = xx;
	    }
	}
	else if (nt == T_PVALUE) {
	    xx = batch_pvalue(sexpr, *pZ, pdinfo, NULL);
	    if (na(xx) || xx == -1.0) {
		err = E_INVARG;
	    } else {
		for (i=0; i<n; i++) 
		    genr->xvec[i] = xx;
	    }
	}
	else if (nt == T_MPOW) {
	    genr->scalar = 0;
	    err = genr_mpow(sexpr, genr->xvec, *pZ, pdinfo);
	    if (err) err = E_INVARG;
	}
#ifdef HAVE_MPFR
	else if (nt == T_MLOG) {
	    genr->scalar = 0;
	    err = genr_mlog(sexpr, genr->xvec, *pZ, pdinfo);
	    if (err) err = E_INVARG;
	}
#endif
	else if (nt == T_COEFF || nt == T_STDERR) {
	    if (pmod == NULL || pmod->list == NULL) {
		return E_INVARG;
	    }
	    lv = _isnumber(sexpr)? atoi(sexpr) : 
		varindex(pdinfo, sexpr);
	    vi = ismatch(lv, pmod->list);
	    if (vi == 1) vi = 0;
	    if (!vi) {
		return E_INVARG;
	    }
	    if (nt == T_COEFF && pmod->coeff != NULL) { 
		for (i=0; i<n; i++) {
		    genr->xvec[i] = pmod->coeff[vi-2];
		}
#ifdef GENR_DEBUG
		fprintf(stderr, "got coeff=%g\n", pmod->coeff[vi-2]);
#endif
	    } else if (pmod->sderr != NULL) {
		for (i=0; i<n; i++) {
		    genr->xvec[i] = pmod->sderr[vi-2];
		}
	    } else {
		return E_INVARG;
	    }
	} 
	else {
	    err = evalexp(sexpr, nt, mvec, genr->xvec, 
			 *pZ, pdinfo, pmod, genr);
	    if (err) {  
		err = E_IGNONZERO;
	    } else {
		for (i=t1; i<=t2; i++) mvec[i] = genr->xvec[i];
		err = domath(genr->xvec, mvec, nt, pdinfo, &genr->scalar);
	    }
	}
	break;

    case R_UNKNOWN: 
	err = E_CASEU;
	break;

    default:
	if (*word != '\0') { 
	    sprintf(gretl_errmsg, 
		    _("%s is not a variable or function"), word);
	}
	err = E_UNSPEC;
	break;

    }  /* end of switch on type2 */

    return err;
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
    int nt, v, ls, nv1, type2, nvtmp = 0;
    int t1 = pdinfo->t1, t2 = pdinfo->t2, n = pdinfo->n;
    char s[MAXLEN], sright[MAXLEN], sleft[MAXLEN];
    char snew[MAXLEN], word[16], s1[MAXLEN];
    char newvar[32], genrs[MAXLEN];
    int err = 0;
    char op0, op1;
    register int i;
    double *mstack = NULL, *mvec = NULL;
    int nv = pdinfo->v;
    GENERATE genr;

    /* mstack cumulates value of expression
       genr.xvec cumulates value of expression inside ()
       mvec gets values for each variable or function
    */

    genr.scalar = 1;
    genr.scalval = NADBL;
    *gretl_errmsg = '\0';
    *genr.label = '\0';

    *s = *genrs = *snew = *sleft = '\0';
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
	err = dummy(pZ, pdinfo);
	if (!err)
	    strcpy(gretl_msg, _("Periodic dummy variables generated.\n"));
	return err;
    }
    else if (strcmp(s, "paneldum") == 0) {
	err = paneldum(pZ, pdinfo, oflag);
	if (!err)
	    strcpy(gretl_msg, _("Panel dummy variables generated.\n"));
	return err;
    }
    else if (strcmp(s, "index") == 0) { 
	err = genrtime(pZ, pdinfo, &genr, 0);
	if (!err) genr_msg(&genr, nv);
	return err;
    }
    else if (strcmp(s, "time") == 0) {
	err = genrtime(pZ, pdinfo, &genr, 1);
	if (!err) genr_msg(&genr, nv);
	return err;
    }
    else if (strncmp(s, "toler=", 6) == 0) {
	err = gentoler(s + 6);
	return err;
    }

    *newvar = '\0';
    op0 = '\0';
    
    /* get equation newvar = s, where s is expression */
    i = getword('=', s, newvar, oflag);

    if (i > 0) {
	if (*newvar == '\0') {
	    err = E_NOVAR;
	    goto genr_return;
	}
	_esl_trunc(newvar, VNAMELEN - 1);
	if (!isalpha((unsigned char) *newvar) &&
	    strncmp(newvar, "$nls", 4)) {
	    err = E_NOTALPH;
	    goto genr_return;
	}
	v = varindex(pdinfo, newvar);
	if (v == 0) { 
	    err = E_CONST;
	    goto genr_return;
	}
	if (lastchar('=', s)) {
	    err = E_EQN;
	    goto genr_return;
	}
    } else {
	err = E_NOEQ;
	goto genr_return;
    }

    /* trap recursive formulae that need to be run though "simulate" */
    if ((nt = has_lhs_lags(s, newvar))) {
	char *simstr = genr_to_sim(s, nt, newvar, pdinfo);

#ifdef GENR_DEBUG
	fprintf(stderr, "genr: has_lhs_lags = %d\n", nt);
#endif
	
	if (simstr != NULL) {
	    int simerr;

	    simerr = simulate(simstr, pZ, pdinfo);
	    free(simstr);
	    return simerr;
	}
    }

    /* basic memory allocation */
    if ((genr.xvec = malloc(n * sizeof *genr.xvec)) == NULL) {
	err = E_ALLOC;
	goto genr_return;
    }

    if ((mstack = malloc(n * sizeof *mstack)) == NULL) {
	err = E_ALLOC;
	goto genr_return;
    } 
    if ((mvec = malloc(n * sizeof *mvec)) == NULL) {
	err = E_ALLOC;
	goto genr_return;
    } 

    for (i=0; i<n; i++) {
	genr.xvec[i] = mstack[i] = mvec[i] = 0.0;
    }

    /* deal with leading (unary) minus */
    if (*s == '-') {
	strcpy(s1, "0");
	strcat(s1, s);
	strcpy(s, s1);
    }

    /* impose operator hierarchy */
    if (parenthesize(s)) { 
	fprintf(stderr, "genr: parenthesize failed\n");
	err = E_ALLOC;
	goto genr_return;
    }
#ifdef GENR_DEBUG
    fprintf(stderr, "after parenthesize: s='%s'\n", s);
#endif

    while ((ls = strlen(s)) > 0) {
	char *indx1, *indx2;

#ifdef GENR_DEBUG
	fprintf(stderr, "s='%s', zeroing mvec, xvec\n", s);
	if (_isnumber(s)) {
	    fprintf(stderr, "Got valid floating point number, '%s'\n", s);
	}
#endif
	for (i=0; i<pdinfo->n; i++) mvec[i] = genr.xvec[i] = 0.0;

	indx1 = strrchr(s, '('); /* point to last '(' */
	if (indx1 == NULL) { /* no left parenthesis */
	    indx2 = strchr(s, ')');
	    if (indx2 != NULL) {
		err = E_UNBAL;
		goto genr_return;
	    }
            getvar(s, s1, &op1);

#ifdef GENR_DEBUG
	    if (isprint((unsigned char) op1)) 
		fprintf(stderr, "after getvar: s='%s', s1='%s', op1='%c'\n",
			s, s1, op1);
	    else
		fprintf(stderr, "after getvar: s='%s', s1='%s', op1=%d\n",
			s, s1, op1);
#endif 

	    if (is_operator(op1) && strlen(s) == 0) {
		err = E_SYNTAX;
		goto genr_return;
	    } 
	    else if (op1 == '\0' || is_operator(op1)) {
		err = getxvec(s1, genr.xvec, *pZ, pdinfo, pmod, &genr.scalar);
		if (err == E_BADSTAT) {
		    err = E_BADSTAT;
		    goto genr_return;
		}
		if (err != 0) {
		    err = E_UNKVAR;
		    goto genr_return;
		}
	    } else {
		err = E_BADOP;
		goto genr_return;
            }
#ifdef GENR_DEBUG
	    fprintf(stderr, "scalar=%d, genr.xvec[1] = %f\n", 
		    genr.scalar, genr.xvec[1]);
#endif
            if (cstack(mstack, genr.xvec, op0, pdinfo, genr.scalar)) {
		err = E_UNSPEC;
		goto genr_return;
	    }
            op0 = op1;
            if (strlen(s) == 0) {
		/* add or replace transformed variable */
                if (v < nv && !oflag && model_count > 0) 
		    sprintf(genr.label, _("Replaced after model %d: "), 
			    model_count);
		if (strlen(genrs) > MAXLABEL - 1) {
		    strncat(genr.label, genrs, MAXLABEL - 4);
		    strcat(genr.label, "...");
		} else {
		    strncat(genr.label, genrs, MAXLABEL - 1);
		}
                for (i=t1; i<=t2; i++) genr.xvec[i] = mstack[i];
                strcpy(genr.varname, newvar);
		genr.varnum = v;
		if (genr.scalar) genr.scalval = genr.xvec[pdinfo->t1];
		genr_msg(&genr, nv);
		goto genr_return;
            }
        } else { /* indx1 != NULL: left paren was found */
	    int iw, nright1, nright2;
	    int nleft1, nleft2;
	    char sexpr[MAXLEN];

            nright1 = strlen(indx1);    /* no. of characters to right of ( */
            nleft1 = ls - nright1;      /* no. of characters before ( */
	    *sleft = '\0';
	    strncat(sleft, s, nleft1);  /* string to left of (  */
            /* calculate equation inside parenthesis */
            strcpy(sright, indx1);         /* string to right of ( */
            indx2 = strchr(sright, ')');  /* point to first ) */
            if (indx2 == NULL) {
                err = E_UNBAL;
		goto genr_return;
            }
            nright2 = strlen(indx2);        /* no chars at end of string */
            nleft2 = nright1 - nright2 - 1; /* no of character inside the (),
					       including */
            indx1++;
            strcpy(sright, indx1);
	    *sexpr = '\0';
	    strncat(sexpr, sright, nleft2);   /* sexpr is expr inside ()  */
            iw = scanb(sleft, word);  /* scan backwards for word in
					 front of ( */

#ifdef GENR_DEBUG
	    fprintf(stderr, "genr: scanb gave word = '%s'\n", word);
#endif

	    if (++nvtmp > MAXTERMS) {
		err = E_NEST;
		goto genr_return;   
	    }

	    nv1 = nv + nvtmp;	    

            if (iw == 0) {
		/* there is an operator in front of (  */
                err = evalexp(sexpr, 0, mvec, genr.xvec, 
			     *pZ, pdinfo, pmod, &genr);
                if (err) {
		    err = E_IGNONZERO;
		    goto genr_return;  
                }
            } else {
		/* there is a math function or lag/lead */
                type2 = strtype(word, pdinfo);
		err = handle_type2(type2, word, sexpr, &genr, mvec,
				   pZ, pdinfo, pmod);
		if (err) goto genr_return; 

		*(sleft + nleft1 - strlen(word)) = '\0';
            } /* end of if (iw == 0) */

	    memmove(sright, indx2, strlen(indx2) + 1);
	    /* create temp var */
	    err = createvar(genr.xvec, snew, sleft, sright, 
			   nv + nvtmp, pZ, pdinfo, genr.scalar);
	    if (err) {
		err = E_UNSPEC; 
		goto genr_return; 
	    }
	    strcpy(s, snew);	    

        }  /* end of if (indx1 == '\0') loop */
    }  /* end of while loop */

 genr_return:
    if (err) {
	genrfree(pZ, pdinfo, &genr, mstack, mvec, nv);
    } else {
	genrfree(pZ, pdinfo, NULL, mstack, mvec, nv);
	err = add_new_var(pdinfo, pZ, &genr);
    }

    return err;
}

/* ........................................................... */
    
static int add_new_var (DATAINFO *pdinfo, double ***pZ, GENERATE *genr)
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

    if (genr->xvec != NULL) free(genr->xvec);

    return 0;
}

/* ........................................................  */

static void expand_vec (double *xx, const DATAINFO *pdinfo)
{
    int t, j;

    for (t=0; t<pdinfo->n; t++) {
#ifdef GENR_DEBUG
	fprintf(stderr, "   expand_vec: vec[%d]=%g\n", t, xx[t]);
#endif	
	if (na(xx[t])) continue;
	else {
	    for (j=0; j<=t; j++) xx[j] = xx[t];
	    break;
	}
    }
}

/* ............................................................ */

static int cstack (double *mstack, double *xvec, char op, 
		   const DATAINFO *pdinfo, int expand)
     /* calculate stack vector */
{
    register int t;
    long int ny;
    double xx, yy, *st2;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;

    st2 = malloc(pdinfo->n * sizeof *st2);
    if (st2 == NULL) return E_ALLOC;

#ifdef GENR_DEBUG
    fprintf(stderr, "in cstack() with op='%c', expand=%d\n", op, expand);
#endif

    if (expand) {
	for (t=0; t<t1; t++) mstack[t] = NADBL;
	for (t=t2+1; t<pdinfo->n; t++) mstack[t] = NADBL;
    }

    for (t=t1; t<=t2; t++) st2[t] = mstack[t];

    switch (op) {
    case '\0':
	for (t=t1; t<=t2; t++) mstack[t] = xvec[t];
	break;
    case '+':
	for (t=t1; t<=t2; t++) mstack[t] += xvec[t];
	break;
    case '|':
	for (t=t1; t<=t2; t++) {
	    mstack[t] = mstack[t] + xvec[t];
	    if (floatneq(mstack[t], 0.)) mstack[t] = 1.0;
	}
	break;
    case '-':
	for (t=t1; t<=t2; t++) mstack[t] -= xvec[t];
	break;
    case '*':
	for (t=t1; t<=t2; t++) mstack[t] *= xvec[t];
	break;
    case '&':
	for (t=t1; t<=t2; t++) {
	    mstack[t] = mstack[t] * xvec[t];
	    if (mstack[t] != 0.) mstack[t] = 1.0;
	}
	break;
    case '%':
	for (t=t1; t<=t2; t++) 
	    mstack[t] = (double) ((int) mstack[t] % (int) xvec[t]);
	break;
    case '/':
	for (t=t1; t<=t2; t++)  {
	    xx = xvec[t];
	    if (floateq(xx, 0.0)) {  
		sprintf(gretl_errmsg, _("Zero denominator for obs %d"), t+1);
		free(st2);
		return 1;
	    }
	    mstack[t] /= xx;
	}
	break;
    case '^':
	for (t=t1; t<=t2; t++) {
	    xx = mstack[t];
	    yy = xvec[t];
	    ny = (long) yy;
	    if ((floateq(xx, 0.0) && yy <= 0.0) || 
		(xx < 0.0 && (double) ny != yy)) {
		sprintf(gretl_errmsg, 
			_("Invalid power function args for obs. %d"
			  "\nbase value = %f, exponent = %f"), t, xx, yy);
		free(st2);
		return 1;
	    }
	    if (floateq(xx, 0.0)) mstack[t] = 0.0;
	    else mstack[t] = pow(xx, yy);
	}
	break;
    case '<':
	for (t=t1; t<=t2; t++) 
            if (mstack[t] < xvec[t]) mstack[t] = 1.0;
            else mstack[t] = 0.0;
	break;
    case '>':
	for (t=t1; t<=t2; t++) { 
#ifdef GENR_DEBUG
	    fprintf(stderr, "comparing mstack[%d]=%g with xvec[%d]=%g\n",
		    t, mstack[t], t, xvec[t]);
#endif
            if (mstack[t] > xvec[t]) mstack[t] = 1.0;
            else mstack[t] = 0.0;
	}
	break;
    case '=':
	for (t=t1; t<=t2; t++) 
            if (floateq(mstack[t], xvec[t])) mstack[t] = 1.0;
            else mstack[t] = 0.0;
	break;
    case NEQ: /* not equals */
	for (t=t1; t<=t2; t++) 
            if (floateq(mstack[t], xvec[t])) mstack[t] = 0.0;
            else mstack[t] = 1.0;
	break;
    case GEQ: /* greater than or equal */
	for (t=t1; t<=t2; t++) {
            if (floateq(mstack[t], xvec[t])) mstack[t] = 1.0;
            else if (mstack[t] > xvec[t]) mstack[t] = 1.0;
	    else mstack[t] = 0.0;
	}
	break;
    case LEQ: /* less than or equal */
	for (t=t1; t<=t2; t++) {
            if (floateq(mstack[t], xvec[t])) mstack[t] = 1.0;
	    else if (mstack[t] < xvec[t]) mstack[t] = 1.0;
            else mstack[t] = 0.0;
	}
	break;
    case '!':
	for (t=t1; t<=t2; t++)
	    if (floatneq(xvec[t], 0.0)) mstack[t] = 0.0;
	    else mstack[t] = 1.0;
	break;
    } /* end of operator switch */

    for (t=t1; t<=t2; t++) 
        if (na(xvec[t]) || na(st2[t])) mstack[t] = NADBL;

    if (expand) {
#ifdef GENR_DEBUG
	fprintf(stderr, "cstack: calling expand_vec() on mstack\n");
#endif
	expand_vec(mstack, pdinfo);
    }

    free(st2);
    return 0;
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

/* ........................................................  */

static int domath (double *xvec, const double *mvec, int nt,
		   const DATAINFO *pdinfo, int *scalar)
     /* do math transformations and return result in xvec */
{
    register int i, t;
    long int xint; 
    double xx = 0.0, yy = 0.0, *x = NULL;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;

    /* mvec contains vector of data to be transformed, result
       returned in xvec */

#ifdef GENR_DEBUG
    fprintf(stderr, "in domath with nt=%d, scalar=%d\n", nt, *scalar);
#endif

    if (*scalar) {
	for (t=0; t<t1; t++) xvec[t] = NADBL;
	for (t=t2+1; t<pdinfo->n; t++) xvec[t] = NADBL;
    }

    switch (nt) {

    case T_LOG:
    case T_LN:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    if (na(xx)) {
		xvec[t] = NADBL;
		continue;
	    }
	    else if (xx <= 0.0) return E_LOGS;
	    xvec[t] = log(xx);
#ifdef GENR_DEBUG
	    fprintf(stderr, "T_LN: set xvec[%d]=%g\n", t, xvec[t]);
#endif
	}
	break;

    case T_EXP:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    if (na(xx)) {
                xvec[t] = NADBL;
                continue;
	    }
	    if (exp(xx) == HUGE_VAL) {
		fprintf(stderr, "genr: T_EXP: exponent = %g\n", xx);
		return E_HIGH;
	    }
	    xvec[t] = exp(xx);
	}
	break;

    case T_SIN:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    xvec[t] = (na(xx))? NADBL: sin(xx);
	}
	break;

    case T_COS:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    xvec[t] = (na(xx))? NADBL: cos(xx);
	}
	break;

    case T_TAN:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    xvec[t] = (na(xx))? NADBL: tan(xx);
	}
	break;

    case T_ATAN:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    xvec[t] = (na(xx))? NADBL: atan(xx);
	}
	break;

    case T_DIFF:
	for (t=t1+1; t<=t2; t++) {
	    if (pdinfo->time_series == STACKED_TIME_SERIES &&
		panel_unit_first_obs(t, pdinfo)) {
		xvec[t] = NADBL;
		continue;
	    }
	    xx = mvec[t];
	    if (pdinfo->time_series == STACKED_CROSS_SECTION) 
		yy = (t - pdinfo->pd >= 0)? mvec[t-pdinfo->pd] : NADBL;
	    else 
		yy = mvec[t-1];
	    xvec[t] = (na(xx) || na(yy))? NADBL : xx - yy;
	}
	xvec[t1] = NADBL;
	break;

    case T_LDIFF:
	for (t=t1+1; t<=t2; t++) {
	    if (pdinfo->time_series == STACKED_TIME_SERIES &&
		panel_unit_first_obs(t, pdinfo)) {
		xvec[t] = NADBL;
		continue;
	    }
	    xx = mvec[t];
	    if (pdinfo->time_series == STACKED_CROSS_SECTION) 
		yy = (t - pdinfo->pd >= 0)? mvec[t-pdinfo->pd] : NADBL;
	    else 
		yy = mvec[t-1];
	    if (na(xx) || na(yy)) {
		xvec[t] = NADBL;
		continue;
	    }   
	    else if (xx <= 0.0 || yy <= 0.0) return E_LOGS;
	    xvec[t] = log(xx) - log(yy);
	}
	xvec[t1] = NADBL;
	break;

    case T_NOBS:
	i = 0;
	for (t=t1; t<=t2; t++) {
	    if (!na(mvec[t])) i++;
	}
	for (t=0; t<pdinfo->n; t++) xvec[t] = (double) i;
	break;

    case T_MEAN: 
    case T_SUM:
    case T_SD:
    case T_VAR:
    case T_SST:
    case T_MEDIAN:
    case T_MIN:
    case T_MAX:
    case T_SORT:
	x = malloc((t2 - t1 + 1) * sizeof *x);
	if (x == NULL) return E_ALLOC;

	i = -1;
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    if (na(xx)) continue;
	    x[++i] = xx;
	}

	if (nt == T_MEAN) {
	    xx = _esl_mean(0, i, x);
	}
	else if (nt == T_SUM) {
	    xx = _esl_mean(0, i, x);
	    xx *= (i + 1);
	}
	else if (nt == T_SD) {
	    xx = _esl_stddev(0, i, x);
	}
	else if (nt == T_VAR) {
	    xx = _esl_variance(0, i, x);
	}
	else if (nt == T_SST) {
	    xx = _esl_sst(0, i, x);
	}
	else if (nt == T_MEDIAN) {
	    xx = gretl_median(x, i+1);
	}
	else if (nt == T_MIN || nt == T_MAX) {
	    double min, max;

	    _minmax(0, i, x, &min, &max);
	    xx = (nt == T_MIN)? min : max;
	}

	if (nt == T_SORT) {
	    qsort(x, i+1, sizeof(double), _compare_doubles);
	    for (t=t1; t<=t2; t++) xvec[t] = x[t-t1];
	} else {
	    for (t=0; t<pdinfo->n; t++) xvec[t] = xx;
	}

	free(x);
	break;

    case T_INT:
	for (t=t1; t<=t2; t++) {
	    xint = (int) (mvec[t] + _VSMALL);
	    if (xint == -999) {
		xvec[t] = NADBL;
		continue;
	    }
	    xvec[t] = (double) xint;
	}
	break;

    case T_ABS:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    if (na(xx)) {
		xvec[t] = NADBL;
		continue;
	    }
	    xvec[t] = (xx < 0.0)? -xx : xx;
	}
	break;

    case T_SQRT:
	for (t=t1; t<=t2; t++) {
	    xx = mvec[t];
	    if (na(xx)) {
		xvec[t] = NADBL;
		continue;
	    }
	    else if (xx < 0.0) return E_SQRT;
	    xvec[t] = sqrt(xx);
	}
	break;

    case T_CUM:  /* cumulate, with "cum" function */
	xvec[t1] = (na(mvec[t1])) ? 0.0 : mvec[t1];
	for (t=t1+1; t<=t2; t++) {
	    if (na(mvec[t])) xvec[t] = xvec[t-1];
	    else xvec[t] = xvec[t-1] + mvec[t];
	}
	break;

    case T_MISSING:  /* check whether obs is missing or not */
	for (t=t1; t<=t2; t++) 
	    xvec[t] = (na(mvec[t])) ? 1.0 : 0.0;
	break;

    case T_MISSZERO:  /* change missing obs to zero */
	for (t=t1; t<=t2; t++) 
	    xvec[t] = (na(mvec[t])) ? 0.0 : mvec[t];
	break;

    case T_ZEROMISS:  /* change zero to missing obs */
	for (t=t1; t<=t2; t++) 
	    xvec[t] = (floateq(mvec[t], 0.0)) ? NADBL : mvec[t];
	break;

    }

    if (*scalar) {
#ifdef GENR_DEBUG
	fprintf(stderr, "domath: calling expand_vec() on xvec\n");
#endif
	expand_vec(xvec, pdinfo);
    }

    return 0;
}

/* .....................................................*/

static int evalexp (char *ss, int nt, double *mvec, double *xvec, 
		    double **Z, const DATAINFO *pdinfo, 
		    const MODEL *pmod, GENERATE *genr)
{
    char s3[MAXLEN], op2, op3;
    int v, expand;
    int *pscalar = &(genr->scalar);
    int err = 0;

    if (SCALAR_SCOPE(nt)) {
	pscalar = NULL;
	genr->scalar = 1;
    }

#ifdef GENR_DEBUG
    fprintf(stderr, "evalexp: ss='%s'\n", ss);
    fprintf(stderr, "         nt = %d, scalar_scope? %s\n", nt,
	   SCALAR_SCOPE(nt)? "Yes" : "No" );
#endif

    /* evaluate expression inside parentheses and value in xvec */

    op3 = '\0';
    do {
	/* don't expand real, vector variables */
	expand = genr->scalar;
	v = varindex(pdinfo, ss);
	if (v == UHATNUM || v == YHATNUM || v == TNUM || v == INDEXNUM ||
	    v == OBSBOOLNUM || (v < pdinfo->v && pdinfo->vector[v])) {
	    expand = 0;
	}
	getvar(ss, s3, &op2);
	if (op2 == '\0' || is_operator(op2)) {
	    err = getxvec(s3, mvec, Z, pdinfo, pmod, pscalar);
	    if (!err) {
#ifdef GENR_DEBUG
		fprintf(stderr, "evalexp, about to do cstack: ss='%s'\n", ss);
#endif	    
		cstack(xvec, mvec, op3, pdinfo, expand);
		op3 = op2;
	    }
        }
    } while (*ss != '\0' && !err);

    return err;
}

/* ........................................................ */

static void getvar (char *str, char *word, char *c)
     /*   
	  Scans str for the first occurrence of { } ( ) or math
	  operator. If found, copies the character into c; copies
	  string to the left into word; and deletes word from str. 
	  If no occurrence, sets word = str, str = '\0', and c = '\0'.
     */
{
    size_t i;

#ifdef GENR_DEBUG
    fprintf(stderr, "genr: getvar: working on '%s'\n", str);
#endif

    /* don't pick apart valid floating-point numbers */
    if (_isnumber(str)) {
#ifdef GENR_DEBUG
	fprintf(stderr, " is_number returned non-zero: set word='%s'\n", 
		word);
#endif
	strcpy(word, str);
	*str = '\0';
	*c = '\0';
	return;
    }	

    *word = '\0';
    for (i=0; i<strlen(str); i++)  { 
	if (str[i] == '{' || str[i] == '}' || str[i] == '(' ||
	    str[i] == ')' || is_operator(str[i])) {
	    *c = str[i];
	    copy(str, 0, i, word);
#ifdef GENR_DEBUG
	    fprintf(stderr, " word='%s', returning\n", word);
#endif
	    _delete(str, 0, i + 1);
	    return;
	}
    }

    strcpy(word, str);

#ifdef GENR_DEBUG
    fprintf(stderr, " word='%s', returning\n", word);
#endif

    *str = '\0';
    *c = '\0';
}

/* ...........................................................*/

static int check_modelstat (const MODEL *pmod, int type1)
{
    if (pmod == NULL || pmod->list == NULL) {
	switch (type1) {
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
	type1 == R_LNL) {
	strcpy(gretl_errmsg, 
	       _("$lnl (log-likelihood) is not available for the last model"));
	return 1;
    }

    return 0;
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

static int getxvec (char *s, double *xvec, 
		    double **Z, const DATAINFO *pdinfo, 
		    const MODEL *pmod, int *scalar)
     /* calculate and return the xvec vector of values */
{
    int type1 = strtype(s, pdinfo);
    int v, n = pdinfo->n;
    register int t;
    double value;

#ifdef GENR_DEBUG
    fprintf(stderr, "in getxvec() with s='%s', type1=%d\n", s, type1);
#endif

    if (check_modelstat(pmod, type1)) return 1;

    if (pmod && (pmod->ci == LOGIT || pmod->ci == PROBIT) &&
	(type1 == R_RSQ || type1 == R_ESS || type1 == R_SIGMA || 
	 type1 == R_TRSQ)) 
	return E_BADSTAT;

    if (pmod && pmod->ci == LAD && (type1 == R_RSQ || type1 == R_TRSQ))
	return E_BADSTAT;

    switch (type1) {  

    case R_ESS:
	for (t=0; t<n; t++) xvec[t] = pmod->ess;
	break;

    case R_NOBS:
	for (t=0; t<n; t++) 
	    xvec[t] = (double) (pdinfo->t2 - pdinfo->t1 + 1);
	break;

    case R_PD:
	for (t=0; t<n; t++) xvec[t] = (double) pdinfo->pd;
	break;

    case R_T:
	for (t=0; t<n; t++) xvec[t] = (double) pmod->nobs;
	break;	

    case R_RSQ:
	for (t=0; t<n; t++) xvec[t] = pmod->rsq;
	break;

    case R_LNL:
	for (t=0; t<n; t++) xvec[t] = pmod->lnL;
	break;

    case R_SIGMA:
	if (pmod->nwt) 
	    for (t=0; t<n; t++) xvec[t] = pmod->sigma_wt;
	else 
	    for (t=0; t<n; t++) xvec[t] = pmod->sigma;
	break;

    case R_TRSQ:
	for (t=0; t<n; t++) xvec[t] = pmod->nobs * pmod->rsq;
	break;

    case R_DF:
	for (t=0; t<n; t++) xvec[t] = (double) pmod->dfd;
	break;

    case R_NUMERIC:
	if (strcmp(s, "pi") == 0) {
	    value = M_PI;
	} else {
	    value = dot_atof(s);
	}
	for (t=0; t<n; t++) xvec[t] = value; 
	break;

    case R_VARNAME:
	v = varindex(pdinfo, s);
	if (v == UHATNUM) { /* model residual */
	    if (pmod->uhat == NULL) return 1;
	    if (pmod->t2 - pmod->t1 + 1 > n ||
		model_sample_issue(pmod, NULL, Z, pdinfo)) {
		strcpy(gretl_errmsg, 
		       _("Can't retrieve uhat: data set has changed"));
		return 1;
	    }	    
	    for (t=0; t<pmod->t1; t++) xvec[t] = NADBL;
	    if (pmod->data != NULL) {
		int t2 = pmod->t2 + get_misscount(pmod);

		for (t=pmod->t1; t<=t2; t++) xvec[t] = pmod->uhat[t]; 
		for (t=t2+1; t<n; t++) xvec[t] = NADBL;
	    } else {
		for (t=pmod->t1; t<=pmod->t2; t++) xvec[t] = pmod->uhat[t]; 
		for (t=pmod->t2 + 1; t<n; t++) xvec[t] = NADBL;
	    }
	    if (scalar != NULL) *scalar = 0;
	}
	else if (v == YHATNUM) { /* model fitted values */
	    if (pmod->yhat == NULL) return 1;
	    if (pmod->t2 - pmod->t1 + 1 > n ||
		model_sample_issue(pmod, NULL, Z, pdinfo)) {
		strcpy(gretl_errmsg, 
		       _("Can't retrieve yhat: data set has changed"));
		return 1;
	    }	    
	    for (t=0; t<pmod->t1; t++) xvec[t] = NADBL;
	    if (pmod->data != NULL) {
		int t2 = pmod->t2 + get_misscount(pmod);

		for (t=pmod->t1; t<=t2; t++) xvec[t] = pmod->yhat[t]; 
		for (t=t2+1; t<n; t++) xvec[t] = NADBL;
	    } else {
		for (t=pmod->t1; t<=pmod->t2; t++) xvec[t] = pmod->yhat[t]; 
		for (t=pmod->t2 + 1; t<n; t++) xvec[t] = NADBL;
	    }
	    if (scalar != NULL) *scalar = 0;
	}
	else if (v == INDEXNUM) { /* internal index variable */
	    int k = genr_scalar_index(0, 0);

	    for (t=0; t<n; t++) xvec[t] = (double) k;
	    if (scalar != NULL) *scalar = 0;
	}
	else if (v == TNUM) { /* auto trend/index variable */
	    if (pdinfo->time_series && pdinfo->pd == 1) {
		/* annual data: let 't' be the year */ 
		for (t=0; t<n; t++) xvec[t] = pdinfo->sd0 + t;
	    } else {
		/* let 't' be the 1-based observation number */
		for (t=0; t<n; t++) xvec[t] = (double) (t + 1);
	    }
	    *scalar = 0;

	}
	else if (v == OBSBOOLNUM) { /* auto-boolean based on obs label */
	    make_obs_dummy(xvec, pdinfo);
	    break;
	}
	else { /* a regular variable */
#ifdef GENR_DEBUG
	    fprintf(stderr, "get_xvec: R_VARNAME: v=%d, name=%s\n",
		    v, pdinfo->varname[v]);
#endif
	    for (t=0; t<n; t++) {
		xvec[t] = (pdinfo->vector[v])? Z[v][t] : Z[v][0];
	    }
	    if (pdinfo->vector[v]) {
		if (scalar != NULL) *scalar = 0;
	    }
	}
	break;

    case R_VAROBS:
	value = get_obs_value(s, Z, pdinfo);
	for (t=0; t<n; t++) {
	    xvec[t] = value;
	}
	break;

    case R_OBSNUM:
	value = obs_num(s, pdinfo);
	for (t=0; t<n; t++) {
	    xvec[t] = value;
	}
	break;

    case R_UNKNOWN:  return 1;

    default:
	if (*s != '\0') {
	    if (strncmp(s, "q#$", 3)) {
		sprintf(gretl_errmsg, _("Undefined variable name '%s' in genr"), s);
	    } else {
		sprintf(gretl_errmsg, _("Syntax error in genr formula"));
	    }
	    return 1;
	}
	break;
    } 

    return 0;
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

/* ......................................................  */

static int scanb (const char *s, char *word)
     /*  scan string right to left for + - * / ^ ( 
	 s is string, n is no. of chars in string, return word to
	 left of operator 
     */
{
    int n = strlen(s);
    int i = n - 1;

    *word = '\0';

    if (i < 0) return 0;

    if (s[i] == '(' || s[i] == '\0' || is_operator(s[i])) {
	word[0] = s[n-1];
	word[1] = '\0';
	return 0;
    }

    for (i=n-1; i>=0; i--) {
	if (s[i] == '(' || s[i] == '\0' || is_operator(s[i])) {
	    strcpy(word, s + i + 1);
	    return 1;
	}
    }

    if (i == -1) {
        strcpy(word, s);
        if (*s == '\0') return 0;
        else return 1;
    }

    return 0;
}

/* ......................................................   */

static int obs_num (const char *s, const DATAINFO *pdinfo)
{
    int t;

    if (pdinfo->markers && pdinfo->S != NULL) {
	for (t=0; t<pdinfo->n; t++) {
	    if (!strcmp(s, pdinfo->S[t])) return t + 1;
	}
    }

    if (pdinfo->time_series == TIME_SERIES) {
	t = dateton(s, pdinfo);
	if (t >= 0) return t + 1;
    }

    return 0;
}

/* ......................................................   */

static int varname_plus_obs (const char *ss, const DATAINFO *pdinfo)
{
    if (strchr(ss, '[') == NULL || strchr(ss, ']') == NULL) {
	return 0;
    } else {
	char vname[VNAMELEN], obs[OBSLEN];

	if (sscanf(ss, "%8[^[][%8[^]]]", vname, obs) != 2) {
	    return 0;
	} else {
	    int i = varindex(pdinfo, vname);

	    if (i == pdinfo->v || !pdinfo->vector[i]) {
		/* invalid varname or scalar */
		return 0; 
	    }
	    if (dateton(obs, pdinfo) < 0) {
		/* invalid observation string */
		return 0;
	    }
	    return 1;
	}
    }
    return 0;
}

/* ......................................................   */

static int strtype (char *ss, const DATAINFO *pdinfo)
     /*  checks whether ss is a number, variable name or transformation */
{
    int i;

#ifdef GENR_DEBUG
    fprintf(stderr, "genr: strtype: working on '%s'\n", ss);
#endif

    if (ss[0] == '$') {
        lower(ss);
        if (strcmp(ss, "$ess") == 0)  
	    return R_ESS;
        if (strcmp(ss, "$nobs") == 0) 
	    return R_NOBS;
        if (strcmp(ss, "$t") == 0) 
	    return R_T;
        if (strcmp(ss, "$rsq") == 0)  
	    return R_RSQ;
	if (strcmp(ss, "$sigma") == 0)  
	    return R_SIGMA;
        if (strcmp(ss, "$df") == 0)   
	    return R_DF;
        if (strcmp(ss, "$lnl") == 0)   
	    return R_LNL;
        if (strcmp(ss, "$nrsq") == 0 || strcmp(ss, "$trsq") == 0) 
	    return R_TRSQ;
	if (strcmp(ss, "$pd") == 0)
	    return R_PD;
    }

    if (_isnumber(ss)) {
#ifdef GENR_DEBUG
	fprintf(stderr, "genr: numeric string = '%s'\n", ss);
#endif
        i = strlen(ss) - 1;
        if (ss[i] == 'e') { 
	    sprintf(gretl_errmsg, 
		    _("Scientific notation not allowed for numbers"));
            return R_UNKNOWN;
        }
        else return R_NUMERIC;
    }

    if (strcmp(ss, "pi") == 0) return R_NUMERIC;

    if (get_function(ss)) return R_MATH;

    i = varindex(pdinfo, ss);
    if (i < pdinfo->v || i == UHATNUM || i == YHATNUM ||
	i == TNUM || i == INDEXNUM || i == OBSBOOLNUM) {
	return R_VARNAME; 
    }

    if (varname_plus_obs(ss, pdinfo)) {
	return R_VAROBS;
    }

    if (obs_num(ss, pdinfo)) {
	return R_OBSNUM;
    }    

    return 0;
}

/* ........................................................  */

static int whichtrans (const char *ss)
{
    register int i;

    for (i=0; ; i++) {
	if (math[i] == NULL) break;
        if (strcmp(ss, math[i]) == 0) return i+1;
    }
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

    if (pdinfo->markers && *varname == '"') {
	set_or_get_obs_marker(varname, OBS_SET);
	return OBSBOOLNUM;
    }

    return pdinfo->v;
}

/* ........................................................ */

static int createvar (double *xvec, char *snew, char *sleft, 
		      char *sright, int ssnum, double ***pZ, 
		      DATAINFO *pdinfo, int scalar)
{
    static char ss[10];
    int mv, t1 = pdinfo->t1, t2 = pdinfo->t2;
    register int t;

    sprintf(ss, "q#$%d", ssnum); 
    mv = varindex(pdinfo, ss);

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    strcpy(pdinfo->varname[mv], ss);
    if (scalar) {
#ifdef GENR_DEBUG
	fprintf(stderr, "createvar: added %s\n", ss);
#endif
	pdinfo->vector[mv] = 0;
	for (t=0; t<pdinfo->n; t++) { 
#ifdef GENR_DEBUG
	    fprintf(stderr, "   putting xvec[%d]=%g into pos %d\n", t, xvec[t], t);
#endif
	    (*pZ)[mv][t] = xvec[t];
	}
    } else {
	for (t=0; t<t1; t++) (*pZ)[mv][t] = NADBL;
	for (t=t1; t<=t2; t++) (*pZ)[mv][t] = xvec[t];
	for (t=t2+1; t<pdinfo->n; t++) (*pZ)[mv][t] = NADBL;
    }

    /* return a new string with the temporary variable name in
       place of the calculated expression */
    strcpy(snew, sleft);
    strcat(snew, ss);
    strcat(snew, sright);

    return 0;
}

/* ........................................................ */

static void genrfree (double ***pZ, DATAINFO *pdinfo, GENERATE *genr,
		      double *mstack, double *mvec, int nv)
{
    int s = pdinfo->v - nv;

    if (s > 0) dataset_drop_vars(s, pZ, pdinfo);
    if (mstack != NULL) free(mstack);
    if (mvec != NULL) free(mvec);
    if (genr != NULL) free(genr->xvec);
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

    err = mp_raise (Z[v], xvec, pdinfo->n, pwr);

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

static void genr_msg (GENERATE *genr, int nv)
{
    sprintf(gretl_msg, "%s %s %s (ID %d)", 
	    (genr->varnum < nv)? _("Replaced") : _("Generated"), 
	    (genr->scalar)? _("scalar") : _("vector"),
	    genr->varname, genr->varnum);
    if (genr->scalar && !na(genr->scalval)) {
	char numstr[24];

	sprintf(numstr, " = %g", genr->scalval);
	strcat(gretl_msg, numstr);
    }
}

/* ......................................................  */

static int ismatch (int lv, const int *list)
{
    int i;

    for (i=list[0]; i>=1; i--)
        if (lv == list[i]) return i;
    return 0;
}

/* .......................................................... */

static void varerror (const char *ss)
     /* print error message for variable not in name list */
{
    sprintf(gretl_errmsg, _("Undefined variable name '%s'"), ss);
    if (!strcmp(ss, "const")) 
        sprintf(gretl_errmsg, _("const cannot be used to store values"));
    else if (!strcmp(ss, "uhat")) 
        sprintf(gretl_errmsg,
		_("uhat can be used only in genr.  First use the command: "
		  "genr newname = uhat"));
    else if (*ss == '$') 
	sprintf(gretl_errmsg, _("Reserved var. names starting with "
				"$ can be used only in genr.\nFirst use the "
				"command:  genr newname = %s"), ss);
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

/* handle lags of the LHS variable in a genr formula */

#define COEFFLEN 16
#define CHUNKLEN 32

struct sim_term {
    int lag;
    char coeff[COEFFLEN];
};

static int has_lhs_lags (const char *str, const char *vname)
{
    char finder[VNAMELEN + 1];
    const char *p = str;
    int c, lag, nlags = 0;
    size_t n;

    sprintf(finder, "%s(", vname);
    n = strlen(finder);
    while ((p = strstr(p, finder))) {
	c = *(p - 1);
	p += n;
	if (c == '+' || c == '-' || c == '*' || c == ' ') {
	    if (sscanf(p, "-%d)", &lag)) nlags++;
	    else break;
	} 
    }

    return nlags;
}

static int is_lagvar_string (const char *s, const char *vname)
{
    char finder[VNAMELEN + 1];
    const char *p;
    size_t n;
    int lag = 0;

    sprintf(finder, "%s(", vname);
    n = strlen(finder);

    if (!strncmp(finder, s, n)) {
	sscanf(s + n, "-%d)", &lag);
    }

    /* check for unwanted remainder */
    p = strchr(s + n, ')');
    if (p == NULL || *(p + 1) != '\0') lag = 0;

    return lag;
}

static int is_numeric (const char *s)
{
    char *test;
    extern int errno;

    errno = 0;

    if (*s == '\0') return 0;

    strtod(s, &test);

    if (*test != '\0' || errno == ERANGE) {
	return 0;
    }

    return 1;
}

static int process_simterm (char *s, const char *vname, 
			    struct sim_term *term, 
			    const DATAINFO *pdinfo)
{
    char *p = s;
    int lag, mult = 0;
    int numcoeff = 0;
    int err = 0;

    term->lag = 0;
    *term->coeff = '\0';

#ifdef GENR_DEBUG
    fprintf(stderr, "process_simterm: looking at '%s'\n", s);
#endif

    while (*p && !mult) {  
	if (*p == '*') {
	    mult = 1;
	    *p = '\0';
	} else p++;
    }

    if ((lag = is_lagvar_string(s, vname))) {
	term->lag = lag;
    } else {
	numcoeff = is_numeric(s);
	strncat(term->coeff, s, COEFFLEN - 1);
    }   

    if (mult) {
	if ((lag = is_lagvar_string(p + 1, vname))) {
	    term->lag = lag;
	} else {
	    numcoeff = is_numeric(p + 1);
	    strncat(term->coeff, p + 1, COEFFLEN - 1);
	}
    }	

    if (*term->coeff == '+') {
	memmove(term->coeff, term->coeff + 1, strlen(term->coeff));
    }

    if (*term->coeff && !numcoeff) { /* is this a valid varname? */
	const char *test = term->coeff;

	if (*test == '-') test++;
	if (varindex(pdinfo, test) == pdinfo->v) {
#ifdef GENR_DEBUG
	    fprintf(stderr, "process_simterm: bad coeff '%s'\n", test);
#endif
	    err = 1;
	}
    }

    return err;
}

static int build_sim_string (const char *vname,
			     struct sim_term *terms, int nt,
			     char **psim)
{
    int i, j;
    int maxlag = 0, err = 0;
    char *sim = NULL;
    size_t simlen = 5;

    if (terms == NULL || nt == 0) return 1;

    /* check for max lag, and for duplicate lags */
    for (i=0; i<nt && !err; i++) {
	simlen += strlen(terms[i].coeff) + 1;
	if (terms[i].lag > maxlag) {
	    maxlag = terms[i].lag;
	}
	for (j=i+1; j<nt && !err; j++) {
	    if (terms[j].lag == terms[i].lag) {
		err = 1;
	    }
	}
    }

    if (err) {
	fprintf(stderr, "build_sim_string: lag %d is duplicated\n",
		terms[j-1].lag);
	return err;
    }

#ifdef GENR_DEBUG
    fprintf(stderr, "build_sim_string: nt = %d, maxlag = %d\n", nt, maxlag);
#endif

    simlen += strlen(vname);
    if (maxlag > nt) simlen += (maxlag - nt) * 2;

#ifdef GENR_DEBUG
    fprintf(stderr, "simlen = %d\n", (int) simlen);
#endif

    sim = malloc(simlen);
    if (sim == NULL) return 1;

    sprintf(sim, "sim %s", vname);
    for (i=0; i<=maxlag; i++) {
	int done = 0;

	for (j=0; j<nt; j++) {
	    if (terms[j].lag == i) {
		strcat(sim, " ");
		strcat(sim, terms[j].coeff);
		done = 1;
		break;
	    }
	}
	if (!done) strcat(sim, " 0");
    }

    *psim = sim;

    return 0;
}

static char *genr_to_sim (const char *s, int nlags, const char *vname,
			  const DATAINFO *pdinfo)
{
    int inparen = 0, err = 0, nt = 0;
    const char *q = s, *p = s;
    char chunk[CHUNKLEN];
    struct sim_term *terms = NULL;
    char *sim = NULL;

    terms = malloc((nlags + 1) * sizeof *terms);
    if (terms == NULL) return NULL;

    while (*p && !err) {
	if (*p == '(') inparen++;
	else if (*p == ')') inparen--;
	if (inparen < 0) err = 1;

	if (!err && !inparen && (p != s) && (*p == '+' || *p == '-')) {
	    if (p - q >= CHUNKLEN || nt > nlags) err = 1;
	    else {
		*chunk = '\0';
		strncat(chunk, q, p - q);
		err = process_simterm(chunk, vname, &terms[nt++], pdinfo);
		q = p;
	    }
	}
	p++;
    }

    if (!err && *q != '\0') {
	if (nt > nlags) err = 1;
	else {
	    *chunk = '\0';
	    strncat(chunk, q, CHUNKLEN - 1);
	    err = process_simterm(chunk, vname, &terms[nt++], pdinfo);
	}
    }

    if (nt != nlags && nt != nlags + 1) {
	fprintf(stderr, "genr_to_sim error: lags = %d but numbers of "
		"terms = %d\n", nlags, nt);
	err = 1;
    }

    if (!err) {
	build_sim_string(vname, terms, nt, &sim);
    } 

    free(terms);

    return sim;
}
