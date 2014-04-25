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

/* lexer module for 'genr' and related commands */

#include "genparse.h"
#include "gretl_func.h"
#include "uservar.h"
#include "gretl_string_table.h"
#include "gretl_bundle.h"

#include <glib.h>

#define NUMLEN 32
#define MAXQUOTE 64

#if GENDEBUG
# define LDEBUG 1
#else
# define LDEBUG 0
#endif

const char *wordchars = "abcdefghijklmnopqrstuvwxyz"
                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                        "0123456789_";

struct str_table {
    int id;
    const char *str;
};

struct str_table consts[] = {
    { CONST_PI,    "$pi" },
    { CONST_NA,    "NA" },
    { CONST_INF,   "$inf" },
    { CONST_WIN32, "WIN32" },
    { CONST_EPS,   "$macheps" },
    { CONST_HAVE_MPI, "$havempi" },
    { CONST_MPI_RANK, "$mpirank" },
    { CONST_MPI_SIZE, "$mpisize" },
    { CONST_N_PROC,   "$nproc" },
    { CONST_SYSINFO,  "$sysinfo" },
    { 0,        NULL }
};

struct str_table dummies[] = {
    { DUM_NULL,    "null" },
    { DUM_DIAG,    "diag" },
    { DUM_DATASET, "dataset" },
    { 0,        NULL }
};

struct str_table dvars[] = {
    { R_NOBS,      "$nobs" },
    { R_NVARS,     "$nvars" },
    { R_PD,        "$pd" },
    { R_DATATYPE,  "$datatype" },
    { R_TEST_STAT, "$test" },
    { R_TEST_PVAL, "$pvalue" },
    { R_TEST_LNL,  "$rlnl" },
    { R_KLNL,      "$kalman_lnl" },
    { R_KS2,       "$kalman_s2" },
    { R_KSTEP,     "$kalman_t" },
    { R_INDEX,     "obs" },
    { R_T1,        "$t1" },
    { R_T2,        "$t2" },
    { R_STOPWATCH, "$stopwatch" },
    { R_PUNIT,     "$unit" },
    { R_OBSMAJ,    "$obsmajor" },
    { R_OBSMIN,    "$obsminor" },
    { R_OBSMIC,    "$obsmicro" },
    { R_DATES,     "$obsdate" },
    { R_WINDOWS,   "$windows" },
    { R_VERSION,   "$version" },
    { R_ERRNO,     "$error" },
    { R_SEED,      "$seed" },
    { R_HUGE,      "$huge" },
    { 0,           NULL },
};

struct str_table mvars[] = {
    { M_ESS,     "$ess" },
    { M_T,       "$T" },
    { M_RSQ,     "$rsq" },
    { M_SIGMA,   "$sigma" },
    { M_DF,      "$df" },
    { M_NCOEFF,  "$ncoeff" },
    { M_LNL,     "$lnl" },
    { M_GMMCRIT, "$gmmcrit" },
    { M_AIC,     "$aic" },
    { M_BIC,     "$bic" },
    { M_HQC,     "$hqc" },
    { M_TRSQ,    "$trsq" },
    { M_DWPVAL,  "$dwpval" },
    { M_FSTT,    "$Fstat" },
    { M_CHISQ,   "$chisq" },
    { M_DIAGTEST, "$diagtest" },
    { M_DIAGPVAL, "$diagpval" },
    { M_UHAT,    "$uhat" },
    { M_YHAT,    "$yhat" },
    { M_LLT,     "$llt" },
    { M_AHAT,    "$ahat" },
    { M_SAMPLE,  "$sample" },
    { M_H,       "$h" },
    { M_COEFF,   "$coeff" },
    { M_SE,      "$stderr" },
    { M_VCV,     "$vcv" },
    { M_RHO,     "$rho" },
    { M_COMPAN,  "$compan" },
    { M_XTXINV,  "$xtxinv" },
    { M_VECG,    "$vecGamma" },
    { M_VMA,     "$vma" },
    { M_FEVD,    "$fevd" },
    { M_EVALS,   "$evals" },
    { M_JALPHA,  "$jalpha" }, 
    { M_JBETA,   "$jbeta" },
    { M_JVBETA,  "$jvbeta" },
    { M_JS00,    "$s00" },
    { M_JS11,    "$s11" },
    { M_JS01,    "$s01" },
    { M_EC,      "$ec" },
    { M_HAUSMAN, "$hausman" },
    { M_SARGAN,  "$sargan" },
    { M_SYSGAM,  "$sysGamma" },
    { M_SYSA,    "$sysA" },
    { M_SYSB,    "$sysB" },
    { M_FCAST,   "$fcast" },
    { M_FCERR,   "$fcerr" },
    { M_COEFF_CI,"$coeff_ci" },
    { M_KLLT,    "$kalman_llt" },
    { M_KUHAT,   "$kalman_uhat" },
    { M_EHAT,    "$ehat" },
    { M_MNLPROBS, "$mnlprobs" },
    { M_XLIST,   "$xlist" },
    { M_YLIST,   "$ylist" },
    { M_COMMAND, "$command" },
    { M_DEPVAR,  "$depvar" },
    { 0,         NULL }
};

struct str_table funcs[] = {
    { F_ABS,      "abs" },
    { F_SIN,      "sin" },
    { F_COS,      "cos" },
    { F_TAN,      "tan" },
    { F_ASIN,     "asin" },
    { F_ACOS,     "acos" },
    { F_ATAN,     "atan" },
    { F_SINH,     "sinh" },
    { F_COSH,     "cosh" },
    { F_TANH,     "tanh" },
    { F_ASINH,    "asinh" },
    { F_ACOSH,    "acosh" },
    { F_ATANH,    "atanh" },
    { F_LOG,      "log" },
    { F_LOG,      "ln" },
    { F_LOG10,    "log10" },
    { F_LOG2,     "log2" },
    { F_EXP,      "exp" },
    { F_SQRT,     "sqrt" },
    { F_DIFF,     "diff" },
    { F_LDIFF,    "ldiff" },
    { F_SDIFF,    "sdiff" },
    { F_LLAG,     "lags" },
    { F_TOINT,    "int" },
    { F_ROUND,    "round" },
    { F_CEIL,     "ceil" },
    { F_FLOOR,    "floor" },
    { F_SORT,     "sort" }, 
    { F_DSORT,    "dsort" }, 
    { F_SORTBY,   "sortby" }, 
    { F_RANKING,  "ranking" },
    { F_ODEV,     "orthdev" },
    { F_NOBS,     "nobs" },
    { F_T1,       "firstobs" },
    { F_T2,       "lastobs" },
    { F_RUNIFORM, "uniform" }, 
    { F_RNORMAL,  "normal" }, 
    { F_CUM,      "cum" }, 
    { F_MISSING,  "missing" },
    { F_DATAOK,   "ok" },        /* opposite of missing */
    { F_MISSZERO, "misszero" },
    { F_LRVAR,    "lrvar" },
    { F_QUANTILE, "quantile" },
    { F_MEDIAN,   "median" },
    { F_GINI,     "gini" },
    { F_ZEROMISS, "zeromiss" },
    { F_SUM,      "sum" },
    { F_SUMALL,   "sumall" },
    { F_MEAN,     "mean" },
    { F_MIN,      "min" },
    { F_MAX,      "max" },
    { F_SD,       "sd" },
    { F_VCE,      "var" },
    { F_SKEWNESS, "skewness" },
    { F_KURTOSIS, "kurtosis" },
    { F_SST,      "sst" },
    { F_CNORM,    "cnorm" },
    { F_DNORM,    "dnorm" },
    { F_QNORM,    "qnorm" },
    { F_GAMMA,    "gammafun" },
    { F_LNGAMMA,  "lngamma" },
    { F_DIGAMMA,  "digamma" },
    { F_RESAMPLE, "resample" },
    { F_PNOBS,    "pnobs" },     /* per-unit nobs in panels */
    { F_PMIN,     "pmin" },      /* panel min */
    { F_PMAX,     "pmax" },      /* panel max */
    { F_PSUM,     "psum" },      /* panel sum */
    { F_PMEAN,    "pmean" },     /* panel mean */
    { F_PXSUM,    "pxsum" },     /* panel x-sectional sum */
    { F_PSD,      "psd" },       /* panel std dev */
    { F_PSHRINK,  "pshrink" },
    { F_HPFILT,   "hpfilt" },    /* Hodrick-Prescott filter */
    { F_BKFILT,   "bkfilt" },    /* Baxter-King filter */
    { F_BWFILT,   "bwfilt" },    /* Butterworth filter */
    { F_FRACDIFF, "fracdiff" },  /* fractional difference */
    { F_BOXCOX,   "boxcox" },    /* Box-Cox transformation */
    { F_COV,      "cov" },
    { F_COR,      "corr" },
    { F_MOVAVG,   "movavg" },
    { F_IMAT,     "I" },
    { F_ZEROS,    "zeros" },
    { F_ONES,     "ones" },
    { F_SEQ,      "seq" },
    { F_REPLACE,  "replace" },
    { F_MUNIF,    "muniform" },
    { F_MNORM,    "mnormal" },
    { F_SUMC,     "sumc" },
    { F_SUMR,     "sumr" },
    { F_PRODC,    "prodc" },
    { F_PRODR,    "prodr" },
    { F_MEANC,    "meanc" },
    { F_MEANR,    "meanr" },
    { F_SDC,      "sdc" },
    { F_MINC,     "minc" },
    { F_MAXC,     "maxc" },
    { F_MINR,     "minr" },
    { F_MAXR,     "maxr" },
    { F_IMINC,    "iminc" },
    { F_IMAXC,    "imaxc" },
    { F_IMINR,    "iminr" },
    { F_IMAXR,    "imaxr" }, 
    { F_FFT,      "fft" },
    { F_FFTI,     "ffti" },
    { F_CMULT,    "cmult" },
    { F_HDPROD,   "hdprod" },
    { F_CDIV,     "cdiv" },
    { F_MCOV,     "mcov" },
    { F_MCORR,    "mcorr" },
    { F_MXTAB,    "mxtab" },
    { F_CDEMEAN,  "cdemean" },
    { F_CHOL,     "cholesky" },
    { F_PSDROOT,  "psdroot" },
    { F_INV,      "inv" },
    { F_INVPD,    "invpd" },
    { F_GINV,     "ginv" },
    { F_DIAG,     "diag" },
    { F_TRANSP,   "transp" },
    { F_VEC,      "vec" },
    { F_VECH,     "vech" },
    { F_UNVECH,   "unvech" },
    { F_UPPER,    "upper" },
    { F_LOWER,    "lower" },
    { F_ROWS,     "rows" },
    { F_COLS,     "cols" },
    { F_DET,      "det" },
    { F_LDET,     "ldet" },
    { F_TRACE,    "tr" },
    { F_NORM1,    "onenorm" },
    { F_INFNORM,  "infnorm" },
    { F_RCOND,    "rcond" },
    { F_RANK,     "rank" },
    { F_QFORM,    "qform" },
    { F_MLAG,     "mlag" },
    { F_QR,       "qrdecomp" },
    { F_EIGSYM,   "eigensym" },
    { F_EIGGEN,   "eigengen" },
    { F_EIGSOLVE, "eigsolve" },
    { F_NULLSPC,  "nullspace" },
    { F_PRINCOMP, "princomp" },
    { F_MEXP,     "mexp" },
    { F_FDJAC,    "fdjac" },
    { F_BFGSMAX,  "BFGSmax" },
    { F_NRMAX,    "NRmax" },
    { F_OBSNUM,   "obsnum" },
    { F_ISSERIES, "isseries" },
    { F_ISLIST,   "islist" },
    { F_ISSTRING, "isstring" },
    { F_ISNULL,   "isnull" },
    { F_LISTLEN,  "nelem" },
    { F_PDF,      "pdf" },
    { F_CDF,      "cdf" },
    { F_INVCDF,   "invcdf" },
    { F_PVAL,     "pvalue" },
    { F_CRIT,     "critical" },
    { F_RANDGEN,  "randgen" },
    { F_MRANDGEN, "mrandgen" },
    { F_RANDGEN1, "randgen1" },
    { F_URCPVAL,  "urcpval" },
    { F_VALUES,   "values" },
    { F_UNIQ,     "uniq" },
    { F_MSHAPE,   "mshape" },
    { F_SVD,      "svd" },
    { F_MOLS,     "mols" },
    { F_MPOLS,    "mpols" },
    { F_MRLS,     "mrls" },
    { F_MREAD,    "mread" },
    { F_MWRITE,   "mwrite" },
    { F_BREAD,    "bread" },
    { F_BWRITE,   "bwrite" },
    { F_MCSEL,    "selifc" },
    { F_MRSEL,    "selifr" },
    { F_POLROOTS, "polroots" },
    { F_DUMIFY,   "dummify" },
    { F_WMEAN,    "wmean" },
    { F_WVAR,     "wvar" },
    { F_WSD,      "wsd" },
    { F_XPX,      "xpx" },
    { F_FILTER,   "filter" },
    { F_KFILTER,  "kfilter" },
    { F_KSMOOTH,  "ksmooth" },
    { F_KSIMUL,   "ksimul" },
    { F_TRIMR,    "trimr" },
    { F_GETENV,   "getenv" },
    { F_NGETENV,  "ngetenv" },
    { F_ARGNAME,  "argname" },
    { F_OBSLABEL, "obslabel" },
    { F_READFILE, "readfile" },
    { F_BACKTICK, "grab" },
    { F_STRSTR,   "strstr" },
    { F_STRSTRIP, "strstrip" },
    { F_STRNCMP,  "strncmp" },
    { F_STRLEN,   "strlen" },
    { F_PRINTF,   "printf" },
    { F_SPRINTF,  "sprintf" },
    { F_SSCANF,   "sscanf" },
    { F_VARNAME,  "varname" },
    { F_VARNUM,   "varnum" },
    { F_TOLOWER,  "tolower" },
    { F_TOUPPER,  "toupper" },
    { F_COLNAMES, "colnames" },
    { F_ROWNAMES, "rownames" },
    { F_LJUNGBOX, "ljungbox" },
    { F_MSORTBY,  "msortby" },
    { F_LINCOMB,  "lincomb" },
    { F_IMHOF,    "imhof" },
    { F_TOEPSOLV, "toepsolv" },
    { F_DSUM,     "diagcat" },
    { F_XMIN,     "xmin" },
    { F_XMAX,     "xmax" },
    { F_CORRGM,   "corrgm" },
    { F_MCOVG,    "mcovg" },
    { F_FCSTATS,  "fcstats" },
    { F_BESSEL,   "bessel" },
    { F_FRACLAG,  "fraclag" },
    { F_MREVERSE, "mreverse" },
    { F_DESEAS,   "deseas" },
    { F_PERGM,    "pergm" },
    { F_IRR,      "irr" },
    { F_NPV,      "npv" },
    { F_LOGISTIC, "logistic" },
    { F_WEEKDAY,  "weekday" },
    { F_KDENSITY, "kdensity" },
    { F_MONTHLEN, "monthlen" },
    { F_EPOCHDAY, "epochday" },
    { F_SETNOTE,  "setnote" },
    { F_INVMILLS, "invmills" },
    { F_POLYFIT,  "polyfit" },
    { F_CHOWLIN,  "chowlin" },
    { F_VARSIMUL, "varsimul" },
    { F_STRSPLIT, "strsplit" },
    { F_INLIST,   "inlist" },
    { F_ERRMSG,   "errmsg" },
    { F_ISCONST,  "isconst" },
    { F_IRF,      "irf" },
    { F_INBUNDLE, "inbundle" },
    { F_STRSUB,   "strsub" },
    { F_REGSUB,   "regsub" },
    { F_COLNAME,  "colname" },
    { F_RANDINT,  "randint" },
    { F_NADARWAT, "nadarwat" },   /* Nadaraya-Watson */
    { F_SIMANN,   "simann" },     /* simulated annealing */
    { F_LOESS,    "loess" },
    { F_FREQ,     "freq" },
    { F_GHK,      "ghk" },
    { F_HALTON,   "halton" },
    { F_IWISHART, "iwishart" },
    { F_ISNAN,    "isnan" },
    { F_TYPESTR,  "typestr" },
    { F_QUADTAB,  "quadtable" },
    { F_AGGRBY,   "aggregate" },
    { F_REMOVE,   "remove" },
    { F_ISODATE,  "isodate" },
    { F_GETLINE,  "getline" },
    { F_TYPEOF,   "typeof" },
    { F_ATOF,     "atof" },
    { F_FIXNAME,  "fixname" },
    { F_ISOCONV,  "isoconv" },
    { F_SUBSTR,   "substr" },
    { F_MPI_SEND, "mpisend" },
    { F_MPI_RECV, "mpirecv" },
    { F_BCAST,    "mpibcast" },
    { F_REDUCE,   "mpireduce" },
    { F_ALLREDUCE, "mpiallred" },
    { F_SCATTER,   "mpiscatter" },
    { F_EASTER,    "easterday" },
    { F_GENSERIES, "genseries" },
    { 0,          NULL }
};

struct str_table func_alias[] = {
    { F_GAMMA,     "gammafunc" },
    { F_GAMMA,     "gamma" },
    { F_PVAL,      "pval" },
    { F_LOG,       "logs" },
    { F_OBSLABEL,  "date" },
    { F_BACKTICK,  "$" },
    { 0,          NULL }
};

int const_lookup (const char *s)
{
    int i;

    for (i=0; consts[i].id != 0; i++) {
	if (!strcmp(s, consts[i].str)) {
	    return consts[i].id;
	}
    }

    /* backward comatibility: to be removed eventually */
    if (!strcmp(s, "pi")) {
	return CONST_PI;
    } else if (!strcmp(s, "macheps")) {
	return CONST_EPS;
    }

    return 0;
}

const char *constname (int c)
{
    int i;

    for (i=0; consts[i].id != 0; i++) {
	if (c == consts[i].id) {
	    return consts[i].str;
	}
    }

    return "unknown";
}

static GHashTable *gretl_function_hash_init (void)
{
    GHashTable *ht;
    int i;

    ht = g_hash_table_new(g_str_hash, g_str_equal);

    for (i=0; funcs[i].str != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) funcs[i].str, 
			    GINT_TO_POINTER(funcs[i].id));
    }

    return ht;
}

enum {
    NO_ALIAS,
    ALLOW_ALIAS
};

static int real_function_lookup (const char *s, int a)
{
    static GHashTable *fht;
    gpointer p;
 
    if (s == NULL) {
	/* cleanup signal */
	if (fht != NULL) {
	    g_hash_table_destroy(fht);
	    fht = NULL;
	}
	return 0;
    }

    if (fht == NULL) {
	fht = gretl_function_hash_init();
    }
    
    p = g_hash_table_lookup(fht, s);
    if (p != NULL) {
	return GPOINTER_TO_INT(p);
    }

    if (a == ALLOW_ALIAS) {
	int i;

	for (i=0; func_alias[i].id != 0; i++) {
	    if (!strcmp(s, func_alias[i].str)) {
		return func_alias[i].id;
	    }
	} 
    }   

    return 0;
}

void gretl_function_hash_cleanup (void)
{
    real_function_lookup(NULL, 0);
}

int function_lookup (const char *s)
{
    return real_function_lookup(s, NO_ALIAS);
}

static int function_lookup_with_alias (const char *s)
{
    return real_function_lookup(s, ALLOW_ALIAS);
}

static const char *funname (int t)
{
    int i;

    for (i=0; funcs[i].id != 0; i++) {
	if (t == funcs[i].id) {
	    return funcs[i].str;
	}
    }

    return "unknown";
}

/* for external purposes (.lang file, manual) */

int gen_func_count (void)
{
    int i;

    for (i=0; funcs[i].id != 0; i++) ;
    return i;
}

const char *gen_func_name (int i)
{
    return funcs[i].str;
}

int model_var_count (void)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) ;
    return i;
}

const char *model_var_name (int i)
{
    return mvars[i].str;
}

int data_var_count (void)
{
    int i, n = 0;

    for (i=0; dvars[i].id != 0; i++) {
	if (dvars[i].str[0] == '$') {
	    n++;
	}
    }

    return n;
}

const char *data_var_name (int i)
{
    return dvars[i].str;
}

const char *gretl_function_complete (const char *s)
{
    size_t n = strlen(s);
    int i;

    for (i=0; funcs[i].str != NULL; i++) {
	if (!strncmp(s, funcs[i].str, n)) {
	    return funcs[i].str;
	}
    }

    return NULL;
}

/* end external stuff */

static int dummy_lookup (const char *s)
{
    int i;

    for (i=0; dummies[i].id != 0; i++) {
	if (!strcmp(s, dummies[i].str)) {
	    return dummies[i].id;
	}
    }

    return 0;
}

const char *dumname (int t)
{
    int i;

    for (i=0; dummies[i].id != 0; i++) {
	if (t == dummies[i].id) {
	    return dummies[i].str;
	}
    }

    return "unknown";
}

static int dvar_lookup (const char *s)
{
    int i;

    for (i=0; dvars[i].id != 0; i++) {
	if (!strcmp(s, dvars[i].str)) {
	    return dvars[i].id;
	}
    }

    return 0;
}

const char *dvarname (int t)
{
    int i;

    for (i=0; dvars[i].id != 0; i++) {
	if (t == dvars[i].id) {
	    return dvars[i].str;
	}
    }

    return "unknown";
}

static int mvar_lookup (const char *s)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {
	if (!strcmp(s, mvars[i].str)) {
	    return mvars[i].id;
	}
    }

    if (!strcmp(s, "$nrsq")) {
	/* alias */
	return M_TRSQ;
    }

    return 0;
}

const char *mvarname (int t)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {
	if (t == mvars[i].id) {
	    return mvars[i].str;
	}
    }

    return "unknown";
}

int genr_function_word (const char *s)
{
    int ret = 0;

    ret = real_function_lookup(s, NO_ALIAS);
    if (!ret) {
	ret = dvar_lookup(s);
    }
    if (!ret) {
	ret = mvar_lookup(s);
    }

    return ret;
}

void undefined_symbol_error (const char *s, parser *p)
{
    parser_print_input(p);

    if (p->ch == '.') {
	gretl_errmsg_sprintf(_("%s: no such object"), s);
    } else {
	gretl_errmsg_sprintf(_("The symbol '%s' is undefined"), s);
    }

    p->err = E_UNKVAR;
}

static void function_noargs_error (const char *s, parser *p)
{
    parser_print_input(p);

    pprintf(p->prn, _("'%s': no argument was given"), s);
    pputc(p->prn, '\n');
    gretl_errmsg_sprintf(_("'%s': no argument was given"), s);

    p->err = E_ARGS;
}

void context_error (int c, parser *p)
{
    if (c != 0) {
	parser_print_input(p);
	pprintf(p->prn, _("The symbol '%c' is not valid in this context\n"), c);
	if (c == '&') {
	    pputs(p->prn, _("(for logical AND, use '&&')\n"));
	} else if (c == '|') {
	    pputs(p->prn, _("(for logical OR, use '||')\n"));
	} else if (c == ',') {
	    p->err = E_PARSE;
	}
    } else {
	const char *s = getsymb(p->sym, p);

	if (s != NULL && *s != '\0') {
	    pprintf(p->prn, _("The symbol '%s' is not valid in this context\n"), 
		    getsymb(p->sym, p));
	}
    }

    if (!p->err) {
	p->err = E_PARSE;
    }
}

static char *get_quoted_string (parser *p)
{
    char *s = NULL;
    int n;

    /* look for a matching (non-escaped) double-quote */
    n = double_quote_position(p->point);

    if (n < 0) {
	/* backward compatibility */
	n = gretl_charpos('"', p->point);
    }

    if (n >= 0) {
	s = gretl_strndup(p->point, n);
	parser_advance(p, n + 1);
    } else {
	parser_print_input(p);
	pprintf(p->prn, _("Unmatched '%c'\n"), '"');
	p->err = E_PARSE;
    }

    if (!p->err) {
	if (p->ch == '.' && *p->point == '$') {
	    /* maybe quoted name of saved object followed by 
	       dollar variable? */
	    p->sym = OVAR;
	} else {
	    p->sym = STR;
	}
    }

    return s;
}

static int might_be_date_string (const char *s, int n)
{
    char test[12];
    int y, m, d;

#if LDEBUG
    fprintf(stderr, "might_be_date_string: s='%s', n=%d\n", s, n);
#endif
    
    if (n > 10) {
	return 0;
    }

    *test = 0;
    strncat(test, s, n);

    if (strspn(s, "1234567890") == n) {
	/* plain integer (FIXME?) */
	return 1;
    } else if (sscanf(s, "%d:%d", &y, &m) == 2) {
	/* quarterly, monthly date */
	return 1;
    } else if (sscanf(s, "%d-%d-%d", &y, &m, &d) == 3) {
	/* daily date? */
	return 1;
    } else if (sscanf(s, "%d/%d/%d", &y, &m, &d) == 3) {
	/* daily date? */
	return 1;
    }

    return 0;
}

NODE *obs_node (parser *p)
{
    NODE *ret = NULL;
    char word[OBSLEN + 2] = {0};
    const char *s = p->point - 1;
    int close;
    int special = 0;
    int t = -1;

    close = gretl_charpos(']', s);

#if LDEBUG
    fprintf(stderr, "obs_node: s='%s', ch='%c', close=%d\n", 
	    s, (char) p->ch, close);
#endif

    if (close == 0) {
	pprintf(p->prn, _("Empty observation []\n"));
	p->err = E_PARSE;
    } else if (close < 0) {
	pprintf(p->prn, _("Unmatched '%c'\n"), '[');
	p->err = E_PARSE;
    } else if (*s == '"' && close < OBSLEN + 2 &&
	       gretl_charpos('"', s+1) == close - 2) {
	/* quoted observation label? */
	strncat(word, s, close);
	special = 1;
    } else if (might_be_date_string(s, close)) {
	strncat(word, s, close);
	special = 1;
    } 

    if (special && !p->err) {
	t = get_t_from_obs_string(word, p->dset);
	if (t >= 0) {
	    /* convert to user-style 1-based index */
	    t++;
	}
    }

    if (t > 0) {
	parser_advance(p, close - 1);
	lex(p);
	ret = newdbl(t);
    } else if (!p->err) {
#if LDEBUG
	fprintf(stderr, "obs_node: first try failed, going for expr\n");
#endif
	lex(p);
	ret = expr(p);
    }

    return ret;
}

static void look_up_string_variable (const char *s, parser *p)
{
    char *val = get_string_by_name(s + 1);

    if (val != NULL) {
	p->idstr = gretl_strdup(s + 1);
	if (p->idstr == NULL) {
	    p->err = E_ALLOC;
	} else {
	    p->uval = val;
	    p->sym = USTR;
	}
    } else {
	undefined_symbol_error(s, p);
    }
}

int is_gretl_accessor (const char *s)
{
    int i, n;

    for (i=0; dvars[i].id != 0; i++) {
	n = strlen(dvars[i].str);
	if (!strncmp(s, dvars[i].str, n)) {
	    return !isalpha(s[n]);
	}
    }

    for (i=0; mvars[i].id != 0; i++) {
	n = strlen(mvars[i].str);
	if (!strncmp(s, mvars[i].str, n)) {
	    return !isalpha(s[n]);
	}
    }

    return 0;
}

static void look_up_dollar_word (const char *s, parser *p)
{
    if ((p->idnum = dvar_lookup(s)) > 0) {
	p->sym = DVAR;
    } else if ((p->idnum = const_lookup(s)) > 0) {
	if (p->idnum == CONST_SYSINFO) {
	    p->sym = BUNDLE;
	    p->idstr = gretl_strdup("$sysinfo");
	    p->uval = get_sysinfo_bundle(&p->err);
	} else {
	    p->sym = CON;
	}
    } else if ((p->idnum = mvar_lookup(s)) > 0) {
	p->sym = MVAR;
    } else {
	undefined_symbol_error(s, p);
    }

#if LDEBUG
    fprintf(stderr, "look_up_dollar_word: '%s' -> %d\n",
	    s, p->idnum);
#endif
}

#ifdef USE_RLIB
# include "gretl_foreign.h"

static int maybe_get_R_function (const char *s)
{
    if (strlen(s) >= 3 && !strncmp(s, "R.", 2)) {
	return get_R_function_by_name(s + 2);
    } else {
	return 0;
    }
}

#else /* !USE_RLIB */
# define maybe_get_R_function(s) (0)
#endif

/* @parsing_query: we want to keep track of the case
   where we're lexing/parsing the branches of a
   ternary "query" expression. When such an expression
   is evaluated, it's OK if the branch _not_ taken
   contains an undefined symbol; indeed, this can
   occur by design, as in

     scalar y = isnull(x) ? 0 : x

   when "x" is in fact undefined.

   We therefore use the "UNDEF" node type to defuse the
   error that would otherwise arise on parsing. An error
   is triggered only if the branch that references the
   UNDEF node is selected (attempting to evaluate an UNDEF
   node automatically throws an error.)
*/ 

static int parsing_query;

void set_parsing_query (int s)
{
    parsing_query = s;
}

static void look_up_word (const char *s, parser *p)
{
    int fsym, err = 0;

    fsym = p->sym = function_lookup_with_alias(s);

    if (p->sym == 0 || p->ch != '(') {
	p->idnum = const_lookup(s);
	if (p->idnum > 0) {
	    p->sym = CON;
	} else {
	    p->idnum = dummy_lookup(s);
	    if (p->idnum > 0) {
		p->sym = DUM;
	    } else {
		GretlType vtype = 0;
		char *bstr;

		if ((p->idnum = current_series_index(p->dset, s)) >= 0) {
		    p->sym = UVEC;
		    p->idstr = gretl_strdup(s);
		} else if (!strcmp(s, "time")) {
		    p->sym = DUM;
		    p->idnum = DUM_TREND;
		} else if ((p->uval = user_var_get_value_and_type(s, &vtype)) != NULL) {
		    if (vtype == GRETL_TYPE_DOUBLE) {
			p->sym = UNUM;
		    } else if (vtype == GRETL_TYPE_MATRIX) {
			p->sym = UMAT;
		    } else if (vtype == GRETL_TYPE_BUNDLE) {
			p->sym = BUNDLE;
		    } else if (vtype == GRETL_TYPE_STRING) {
			p->sym = USTR;
		    } else if (vtype == GRETL_TYPE_LIST) {
			p->sym = ULIST;
		    } else if (vtype == GRETL_TYPE_STRING) {
			p->sym = USTR;
		    }
		    p->idstr = gretl_strdup(s);
		} else if ((bstr = get_built_in_string_by_name(s))) {
		    /* FIXME should use $-accessors? */
		    p->sym = STR;
		    p->idstr = gretl_strdup(bstr);
		} else if (gretl_get_object_by_name(s)) {
		    p->sym = UOBJ;
		    p->idstr = gretl_strdup(s);
		} else if (get_user_function_by_name(s)) {
		    p->sym = UFUN;
		    p->idstr = gretl_strdup(s);
		} else if (p->targ == LIST && varname_match_any(p->dset, s)) {
		    p->sym = WLIST;
		    p->idstr = gretl_strdup(s);
		} else if (!strcmp(s, "t")) {
		    /* if "t" has not been otherwise defined, treat it
		       as an alias for "obs"
		    */
		    p->sym = DVAR;
		    p->idnum = R_INDEX;
		} else if (maybe_get_R_function(s)) {
		    /* note: all "native" types take precedence over this */
		    p->sym = RFUN;
		    p->idstr = gretl_strdup(s + 2);
		} else if (parsing_query) {
		    p->sym = UNDEF;
		    p->idstr = gretl_strdup(s);
		} else {
		    err = E_UNKVAR;
		}
	    }
	}
    }

    if (err) {
	if (fsym) {
	    function_noargs_error(s, p);
	} else {
	    undefined_symbol_error(s, p);
	}
    }
}

static void maybe_treat_as_postfix (parser *p)
{
    if (p->sym == UNUM) {
	const char *ok = ")]}+-*/%,:";
	int c = parser_next_nonspace_char(p, 1);

	/* Interpret as foo++ or foo-- ? Only if
	   the following character is suitable.
	*/
	if (c == 0 || strchr(ok, c)) {
	    p->sym = p->ch == '+'? UNUM_P : UNUM_M;
	    /* swallow the pluses or minuses */
	    parser_advance(p, 1);
	}
    }
}

#define dollar_series(t) (t > R_SCALAR_MAX && t < R_SERIES_MAX)

#define could_be_matrix(t) (model_data_matrix(t) || \
			    model_data_matrix_builder(t) || \
			    t == M_UHAT || t == M_YHAT || \
			    t == R_TEST_STAT || t == R_TEST_PVAL)

static void word_check_next_char (parser *p)
{
#if LDEBUG
    if (p->ch) fprintf(stderr, "word_check_next_char: ch = '%c'\n", p->ch);
    else fprintf(stderr, "word_check_next_char: ch = NUL\n");
#endif

    if (p->ch == '(') {
	/* series (lag) or function */
	if (p->sym == UVEC) {
	    if (p->idnum > 0 && p->idnum == p->lh.v) {
		p->flags |= P_AUTOREG;
	    }
	    p->sym = LAG;
	} else if (p->sym == MVAR && model_data_matrix(p->idnum)) {
	    /* old-style "$coeff(x1)" etc. */
	    p->sym = DMSTR;
	} else if (!func1_symb(p->sym) && 
		   !func2_symb(p->sym) &&
		   !func3_symb(p->sym) &&
		   !funcn_symb(p->sym) && 
		   p->sym != UFUN && 
		   p->sym != RFUN) {
	    p->err = E_PARSE;
	} 
    } else if (p->ch == '[') {
	if (p->sym == UMAT) {
	    /* slice of user matrix */
	    p->sym = MSL;
	} else if ((p->sym == MVAR || p->sym == DVAR) && 
		   could_be_matrix(p->idnum)) {
	    /* slice of $ matrix */
	    p->sym = DMSL;
	} else if (p->sym == UVEC) {
	    /* observation from series */
	    p->sym = OBS;
	} else if (p->sym == DVAR && dollar_series(p->idnum)) {
	    /* observation from "dollar" series */
	    p->sym = DOBS;
	} else if (p->sym == ULIST) {
	    /* element of list */
	    p->sym = LISTELEM;
	} else if (p->sym == MVAR && model_data_list(p->idnum)) {
	    /* element of accessor list */
	    p->sym = MLISTELEM;
	} else if (p->sym == BUNDLE) {
	    /* object from bundle */
	    p->sym = BOBJ;
	} else {
	    p->err = E_PARSE;
	} 
    } else if (p->ch == '.' && *p->point == '$') {
	if (p->sym == UOBJ) {
	    /* name of saved object followed by dollar variable? */
	    p->sym = OVAR;
	} else if (p->sym == STR) {
	    /* maybe quoted name of saved object followed by 
	       dollar variable? */
	    p->sym = OVAR;
	} else {
	    p->err = E_PARSE;
	}	    
    } else if (p->ch == '.' && isalpha(*p->point)) {
	if (p->sym == ULIST) {
	    p->sym = LISTVAR;
	} else if (p->sym == BUNDLE) {
	    p->sym = BMEMB;
	} else {
	    p->err = E_PARSE;
	}
    } else if (p->ch == '+' && *p->point == '+') {
	maybe_treat_as_postfix(p);
    } else if (p->ch == '-' && *p->point == '-') {
	maybe_treat_as_postfix(p);
    }	

    if (p->err) {
	context_error(p->ch, p);
    } 
}

static int is_word_char (parser *p)
{
    if (strchr(wordchars, p->ch) != NULL) {
	return 1;
    } else if (p->targ == LIST && p->ch == '*') {
	return 1;
    } 

    return 0;
}

static void getword (parser *p)
{  
    char word[32];
    int i = 0;

    /* we know the first char is acceptable (and might be '$' or '@') */
    word[i++] = p->ch;
    parser_getc(p);

#ifdef USE_RLIB
    /* allow for R.foo function namespace */
    if (*word == 'R' && p->ch == '.' && *p->point != '$') {
	word[i++] = p->ch;
	parser_getc(p);
    }
#endif

    while (p->ch != 0 && is_word_char(p) && i < 31) {
	word[i++] = p->ch;
	parser_getc(p);
    }

    word[i] = '\0';

#if LDEBUG
    fprintf(stderr, "getword: word = '%s'\n", word);
#endif

    while (p->ch != 0 && strchr(wordchars, p->ch) != NULL) {
	/* flush excess word characters */
	parser_getc(p);
    }

    if (p->flags & P_GETSTR) {
	/* uninterpreted string wanted */
	p->sym = STR;
	p->idstr = gretl_strdup(word);
	p->flags ^= P_GETSTR;
	return;
    }

    if ((*word == '$' && word[1]) || !strcmp(word, "obs")) {
	look_up_dollar_word(word, p);
    } else if (*word == '@') {
	/* do we actually want to do this? */
	look_up_string_variable(word, p);
    } else if (*word == '$' && word[1] == '\0' && p->ch == '[') {
	p->sym = BUNDLE;
	p->idstr = gretl_strdup("$");
    } else {
	look_up_word(word, p);
    }

    if (!p->err && *word != '@') {
	word_check_next_char(p);
    }

#if LDEBUG
    fprintf(stderr, "getword: p->err = %d\n", p->err);
#endif
}

static int colon_ok (parser *p, char *s, int n)
{
    int i;

    if (p->flags & P_SLICING) {
	/* calculating a matrix "slice": colon is a separator in 
	   this context, cannot be part of a time/panel observation 
	   string
	*/
#if LDEBUG
	fprintf(stderr, "colon_ok: doing matrix slice\n");
#endif
	return 0;
    }

    if (n != 1 && n != 3) {
	return 0;
    }

    for (i=0; i<=n; i++) {
	if (!isdigit(s[i])) {
	    return 0;
	}
    }

    return 1;
}

/* below: we're testing 'ch' for validity, given what we've already
   packed into string 's' up to element 'i'
*/

static int ok_dbl_char (parser *p, char *s, int i)
{
    if (i < 0) {
	return 1;
    }

    if (p->ch >= '0' && p->ch <= '9') {
	return 1;
    }

    switch (p->ch) {
    case '+':
    case '-':
	return s[i] == 'e' || s[i] == 'E';
    case '.':
	return !strchr(s, '.') && !strchr(s, ':') &&
	    !strchr(s, 'e') && !strchr(s, 'E');
    case 'e':
    case 'E':
	return !strchr(s, 'e') && !strchr(s, 'E') && 
	    !strchr(s, ':');
    case ':':
	/* allow for obs numbers in the form, e.g., "1995:10" */
	return colon_ok(p, s, i);
    default:
	break;
    }

    return 0;
}

static void parse_number (parser *p)
{
    char xstr[NUMLEN] = {0};
    int gotcol = 0;
    int i = 0;

    while (ok_dbl_char(p, xstr, i - 1) && i < NUMLEN - 1) {
	xstr[i++] = p->ch;
	if (p->ch == ':') {
	    gotcol = 1;
	}
	parser_getc(p);
    }  

    while (p->ch >= '0' && p->ch <= '9') {
	/* flush excess numeric characters */
	parser_getc(p);
    } 

#if LDEBUG
    fprintf(stderr, "parse_number: xstr = '%s'\n", xstr);
#endif
    
    if (gotcol) {
#if LDEBUG
	fprintf(stderr, " got colon: obs identifier?\n");
#endif
	if (p->dset == NULL || p->dset->n == 0) {
	    p->err = E_NODATA;
	} else if (p->dset->pd == 1) {
	    p->err = E_PDWRONG;
	} else if (dateton(xstr, p->dset) < 0) {
	    p->err = E_DATA;
	} else {
	    p->idstr = gretl_strdup(xstr);
	    p->sym = STR;
	}
    } else {
	p->xval = dot_atof(xstr);
	p->sym = NUM;
#if LDEBUG
	fprintf(stderr, " dot_atof gave %g\n", p->xval);
#endif
    }
}

#define word_start_special(c) (c == '$' || c == '@' || c == '_')

#define lag_range_sym(p) ((p->flags & P_LAGPRSE) && p->ch == 't' && \
                          *p->point == 'o' && \
			  *(p->point + 1) == ' ')

void lex (parser *p)
{
#if LDEBUG
    if (p->ch) fprintf(stderr, "lex: p->ch = '%c'\n", p->ch);
    else fprintf(stderr, "lex: p->ch = NUL\n");
#endif

    if (p->ch == 0) {
	p->sym = EOT;
	return;
    }

    while (p->ch != 0) {
	switch (p->ch) {
	case ' ':
	case '\t':
	case '\r':
        case '\n': 
	    parser_getc(p);
	    break;
        case '+': 
	    p->sym = B_ADD;
	    parser_getc(p);
	    return;
        case '-': 
	    p->sym = B_SUB;
	    parser_getc(p);
	    return;
        case '*': 
	    if (p->targ == LIST) {
		/* treat '*' as wildcard */
		getword(p);
		return;
	    }
	    parser_getc(p);
	    if (p->ch == '*') {
		p->sym = B_KRON;
		parser_getc(p);
	    } else {
		p->sym = B_MUL;
	    }
	    return;
	case '\'':
	    p->sym = B_TRMUL;
	    parser_getc(p);
	    return;
        case '/': 
	    p->sym = B_DIV;
	    parser_getc(p);
	    return;
        case '\\': 
	    p->sym = B_LDIV;
	    parser_getc(p);
	    return;
        case '%': 
	    p->sym = B_MOD;
	    parser_getc(p);
	    return;
        case '^': 
	    p->sym = B_POW;
	    parser_getc(p);
	    return;
        case '&': 
	    parser_getc(p);
	    if (p->ch == '&') {
		p->sym = B_AND;
		parser_getc(p);
	    } else {
		p->sym = U_ADDR;
	    }
	    return;
        case '|': 
	    parser_getc(p);
	    if (p->ch == '|') {
		p->sym = B_OR;
		parser_getc(p);
	    } else {
		p->sym = B_VCAT;
	    }
	    return;
        case '!': 
	    parser_getc(p);
	    if (p->ch == '=') {
		p->sym = B_NEQ;
		parser_getc(p);
	    } else {
		p->sym = U_NOT;
	    }
	    return;
        case '=': 
	    parser_getc(p);
	    if (p->ch == '=') {
		/* allow "==" as synonym for "=" */
		parser_getc(p);
	    }
	    p->sym = B_EQ;
	    return;
        case '>': 
	    parser_getc(p);
	    if (p->ch == '=') {
		p->sym = B_GTE;
		parser_getc(p);
	    } else {
		p->sym = B_GT;
	    }
	    return;
        case '<': 
	    parser_getc(p);
	    if (p->ch == '=') {
		p->sym = B_LTE;
		parser_getc(p);
	    } else if (p->ch == '>') {
		p->sym = B_NEQ;
		parser_getc(p);
	    } else {
		p->sym = B_LT;
	    }
	    return;
        case '(': 
	    p->sym = G_LPR;
	    parser_getc(p);
	    return;
        case ')': 
	    p->sym = G_RPR;
	    parser_getc(p);
	    return;
        case '[': 
	    p->sym = G_LBR;
	    parser_getc(p);
	    return;
        case '{': 
	    p->sym = G_LCB;
	    parser_getc(p);
	    return;
        case '}': 
	    p->sym = G_RCB;
	    parser_getc(p);
	    return;
        case ']': 
	    p->sym = G_RBR;
	    parser_getc(p);
	    return;
        case '~':
	    p->sym = B_HCAT;
	    parser_getc(p);
	    return;
        case ',': 
	    p->sym = P_COM;
	    parser_getc(p);
	    return;
        case ';': 
	    if (p->targ == LIST) {
		p->sym = B_JOIN;
	    } else {
		/* used in matrix definition */
		p->sym = P_SEMI;
	    }
	    parser_getc(p);
	    return;
        case ':': 
	    p->sym = P_COL;
	    parser_getc(p);
	    return;
        case '?': 
	    p->sym = QUERY;
	    parser_getc(p);
	    return;
	case '.':
	    if (*p->point == '$') {
		p->sym = P_DOT;
		parser_getc(p);
		return;
	    }
	    parser_getc(p);
	    if (p->ch == '*') {
		p->sym = B_DOTMULT;
		parser_getc(p);
		return;
	    } else if (p->ch == '/') {
		p->sym = B_DOTDIV;
		parser_getc(p);
		return;
	    } else if (p->ch == '^') {
		p->sym = B_DOTPOW;
		parser_getc(p);
		return;
	    } else if (p->ch == '+') {
		p->sym = B_DOTADD;
		parser_getc(p);
		return;
	    } else if (p->ch == '-') {
		p->sym = B_DOTSUB;
		parser_getc(p);
		return;
	    } else if (p->ch == '=') {
		p->sym = B_DOTEQ;
		parser_getc(p);
		return;
	    } else if (p->ch == '>') {
		p->sym = B_DOTGT;
		parser_getc(p);
		if (p->ch == '=') {
		    p->sym = B_DOTGTE;
		    parser_getc(p);
		}
		return;
	    } else if (p->ch == '<') {
		p->sym = B_DOTLT;
		parser_getc(p);
		if (p->ch == '=') {
		    p->sym = B_DOTLTE;
		    parser_getc(p);
		}
		return;
	    } else if (p->ch == '.') {
		p->sym = B_ELLIP;
		parser_getc(p);
		return;
	    } else {
		/* not a "dot operator", so back up */
		parser_ungetc(p);
	    }
        default: 
	    if (p->targ == LIST && lag_range_sym(p)) {
		p->sym = B_RANGE;
		parser_getc(p);
		parser_getc(p);
		return;
	    }
	    if (p->targ == LIST && *(p->point - 2) == ' ' && 
		(bare_data_type(p->sym) || closing_sym(p->sym) ||
		 (p->sym == LAG))) {
		/* may be forming a list, but only if there are 
		   spaces between the terms */
		p->sym = B_LCAT;
		return;
	    }
	    if (isdigit(p->ch) || (p->ch == '.' && isdigit(*p->point))) {
		parse_number(p);
		return;
	    } else if (islower(p->ch) || isupper(p->ch) || 
		       word_start_special(p->ch)) {
		getword(p);
		return;
	    } else if (p->ch == '"') {
		p->idstr = get_quoted_string(p);
		return;
	    } else {
		parser_print_input(p);
		pprintf(p->prn, _("Invalid character '%c'\n"), p->ch);
		p->err = E_PARSE;
		return;
	    }
	} /* end ch switch */
    } /* end while ch != 0 */
}

const char *getsymb (int t, const parser *p)
{  
    if ((t > F1_MIN && t < F1_MAX) ||
	(t > F1_MAX && t < F2_MAX) ||
	(t > F2_MAX && t < FN_MAX)) {
	return funname(t);
    }

    if (t == EOT) {
	return "";
    }

    /* yes, well */
    if (t == OBS || t == DOBS) {
	return "OBS";
    } else if (t == MSL) {
	return "MSL";
    } else if (t == DMSL) {
	return "DMSL";
    } else if (t == DMSTR) {
	return "DMSTR";
    } else if (t == MSL2) {
	return "MSL2";
    } else if (t == MSPEC) {
	return "MSPEC";
    } else if (t == SUBSL) {
	return "SUBSL";
    } else if (t == MDEF) {
	return "MDEF";
    } else if (t == FARGS) {
	return "FARGS";
    } else if (t == LIST || t == ULIST || t == WLIST) {
	return "LIST";
    } else if (t == OVAR) {
	return "OVAR";
    } else if (t == USTR) {
	return "USTR";
    } else if (t == EMPTY) {
	return "EMPTY";
    } else if (t == LISTVAR) {
	return "LISTVAR";
    } else if (t == BOBJ) {
	return "BOBJ";
    } else if (t == BMEMB) {
	return "BMEMB";
    } else if (t == LISTELEM) {
	return "LISTELEM";
    } else if (t == VEC) {
	return "VEC";
    } else if (t == MAT) {
	return "MAT";
    } else if (t == EROOT) {
	return "EROOT";
    } else if (t == UNDEF) {
	return "UNDEF";
    } else if (t == NUM) {
	return "NUM";
    }

    if (p != NULL) {
	if (t == UVEC) {
	    return p->dset->varname[p->idnum];
	} else if (t == UNUM) {
	    return p->idstr;
	} else if (t == BUNDLE) {
	    return p->idstr;
	} else if (t == UMAT || t == UOBJ) {
	    return p->idstr;
	} else if (t == CON) {
	    return constname(p->idnum);
	} else if (t == DUM) {
	    return dumname(p->idnum);
	} else if (t == DVAR) {
	    return dvarname(p->idnum);
	} else if (t == MVAR) {
	    return mvarname(p->idnum);
	} else if (t == UFUN || t == RFUN) {
	    return p->idstr;
	} else if (t == STR) {
	    return p->idstr;
	}
    } else {
	if (t == UVEC) {
	    return "UVEC";
	} else if (t == UNUM) {
	    return "UNUM";
	} else if (t == BUNDLE) {
	    return "BUNDLE";
	} else if (t == UMAT) {
	    return "UMAT";
	} else if (t == UOBJ) {
	    return "UOBJ";
	} else if (t == CON) {
	    return "CON";
	} else if (t == DUM) {
	    return "DUM";
	} else if (t == DVAR) {
	    return "DVAR";
	} else if (t == MVAR) {
	    return "MVAR";
	} else if (t == UFUN) {
	    return "UFUN";
	} else if (t == RFUN) {
	    return "RFUN";
	} else if (t == STR) {
	    return "STR";
	}
    }	

    switch (t) {
    case B_ASN:
	return "=";
    case B_ADD: 
    case U_POS:
	return "+";
    case B_SUB: 
    case U_NEG:
	return "-";
    case B_MUL: 
	return "*";
    case B_TRMUL: 
	return "'";
    case B_DIV: 
	return "/";
    case B_LDIV:
	return "\\";
    case B_MOD: 
	return "%";
    case B_POW: 
	return "^";
    case B_EQ: 
	return "=";
    case B_NEQ: 
	return "!=";
    case B_GT: 
	return ">";
    case B_LT: 
	return "<";
    case B_GTE: 
	return ">=";
    case B_LTE: 
	return "<=";
    case B_AND: 
	return "&&";
    case B_RANGE:
	return " to ";
    case U_ADDR:
	return "&";
    case B_OR: 
	return "||";	
    case U_NOT: 
	return "!";
    case G_LPR: 
	return "(";
    case G_RPR: 
	return ")";
    case G_LBR: 
	return "[";
    case G_RBR: 
	return "]";
    case G_LCB: 
	return "{";
    case G_RCB: 
	return "}";
    case B_DOTMULT: 
	return ".*";
    case B_DOTDIV: 
	return "./";
    case B_DOTPOW: 
	return ".^";
    case B_DOTADD: 
	return ".+";
    case B_DOTSUB: 
	return ".-";
    case B_DOTEQ: 
	return ".=";
    case B_DOTGT: 
	return ".>";
    case B_DOTLT: 
	return ".<";
    case B_DOTGTE: 
	return ".>=";
    case B_DOTLTE: 
	return ".<=";
    case B_KRON: 
	return "**";
    case B_HCAT: 
	return "~";
    case B_VCAT: 
	return "|";
    case B_LCAT:
	return "LCAT";
    case P_COM: 
	return ",";
    case P_DOT: 
	return ".";
    case P_SEMI: 
	return ";";
    case P_COL: 
	return ":";
    case QUERY: 
	return "query";
    case LAG:
	return "lag";
    default: 
	break;
    }

    return "unknown";
}
