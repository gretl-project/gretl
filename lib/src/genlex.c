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
#include "gretl_normal.h"
#include "gretl_bundle.h"
#include "uservar_priv.h"
#include "gretl_cmatrix.h"

#define NUMLEN 32
#define HEXLEN 11
#define MAXQUOTE 64

#if GENDEBUG
# define LDEBUG 1
#else
# define LDEBUG 0
#endif

static int parser_next_char (parser *p);

#define defining_list(p) (p->flags & P_LISTDEF)

#define bare_data_type(s) (s > PUNCT_MAX && s < DTYPE_MAX)

#define closing_sym(s) (s == G_RPR || s == G_RBR || s == G_RCB)

const char *wordchars = "abcdefghijklmnopqrstuvwxyz"
                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                        "0123456789_";

struct str_table {
    int id;
    const char *str;
};

struct str_table_ex {
    int id;
    const char *str;
    void *ptr;
};

struct str_table consts[] = {
    { CONST_PI,       "$pi" },
    { CONST_NA,       "NA" },
    { CONST_INF,      "$inf" },
    { CONST_NAN,      "$nan" },
    { CONST_WIN32,    "WIN32" },
    { CONST_EPS,      "$macheps" },
    { CONST_HAVE_MPI, "$havempi" },
    { CONST_MPI_RANK, "$mpirank" },
    { CONST_MPI_SIZE, "$mpisize" },
    { CONST_N_PROC,   "$nproc" },
    { CONST_TRUE,     "TRUE" },
    { CONST_FALSE,    "FALSE" },
    { 0,        NULL }
};

struct str_table dummies[] = {
    { DUM_NULL,    "null" },
    { DUM_EMPTY,   "empty" },
    { DUM_DIAG,    "diag" },
    { DUM_UPPER,   "upper" },
    { DUM_LOWER,   "lower" },
    { DUM_REAL,    "real" },
    { DUM_IMAG,    "imag" },
    { DUM_END,     "end" },
    { DUM_DATASET, "dataset" },
    { 0,        NULL }
};

/* Identify matrix-selection dummy constants:
   these can be valid only between '[' and ']'.
*/
#define MSEL_DUM(d) (d >= DUM_DIAG && d <= DUM_END)

/* dvars: dataset- and test-related accessors */

struct str_table dvars[] = {
    { R_NOBS,      "$nobs" },
    { R_NVARS,     "$nvars" },
    { R_PD,        "$pd" },
    { R_PANEL_PD,  "$panelpd" },
    { R_T1,        "$t1" },
    { R_T2,        "$t2" },
    { R_TMAX,      "$tmax" },
    { R_DATATYPE,  "$datatype" },
    { R_TEST_STAT, "$test" },
    { R_TEST_PVAL, "$pvalue" },
    { R_TEST_BRK,  "$qlrbreak" },
    { R_TEST_LNL,  "$rlnl" },
    { R_STOPWATCH, "$stopwatch" },
    { R_TIME,      "$time" },
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
    { R_NOW,       "$now" },
    { R_RESULT,    "$result" },
    { R_PNGFONT,   "$pngfont" },
    { R_MAPFILE,   "$mapfile" },
    { R_MAP,       "$map" },
    { R_INDEX,     "obs" },
    { R_LOGLEVEL,  "$loglevel" },
    { R_LOGSTAMP,  "$logstamp" },
    { 0,           NULL },
};

/* mvars: model-related accessors */

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
    { M_DW,      "$dw" },
    { M_DWPVAL,  "$dwpval" },
    { M_FSTT,    "$Fstat" },
    { M_CHISQ,   "$chisq" },
    { M_DIAGTEST, "$diagtest" },
    { M_DIAGPVAL, "$diagpval" },
    { M_PMANTEAU, "$pmanteau" },
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
    { M_FCSE,    "$fcse" },
    { M_COEFF_CI,"$coeff_ci" },
    { M_EHAT,    "$ehat" },
    { M_ODDSRATIOS, "$oddsratios" },
    { M_MNLPROBS, "$mnlprobs" }, /* legacy */
    { M_ALLPROBS, "$allprobs" },
    { M_XLIST,   "$xlist" },
    { M_YLIST,   "$ylist" },
    { M_COMMAND, "$command" },
    { M_DEPVAR,  "$depvar" },
    { M_PARNAMES, "$parnames" },
    { 0,         NULL }
};

/* bvars: bundle accessors */

struct str_table bvars[] = {
    { B_MODEL,   "$model" },
    { B_SYSTEM,  "$system" },
    { B_SYSINFO, "$sysinfo" },
    { 0,         NULL }
};

/* Below, @ptrfuncs: table of functions for which we wish to
   attach function-pointers to the relevant NODE. Nota bene:
   it's crucial that no function in @ptrfuncs is also listed
   in @funcs below!

   The order of function symbols in @ptrfuncs need not match
   the order in which they're listed in genparse.h, but it's
   crucial that every function with ID number <= FP_MAX has
   an entry in @ptrfuncs.

   "Crucial" -> certain crash on calling wrongly classified
   function!
*/

struct str_table_ex ptrfuncs[] = {
    { F_ABS,   "abs",   fabs },
    { F_CEIL,  "ceil",  ceil },
    { F_FLOOR, "floor", floor },
    { F_SIN,   "sin",   sin },
    { F_COS,   "cos",   cos },
    { F_TAN,   "tan",   tan },
    { F_ASIN,  "asin",  asin },
    { F_ACOS,  "acos",  acos },
    { F_ATAN,  "atan",  atan },
    { F_SINH,  "sinh",  sinh },
    { F_COSH,  "cosh",  cosh },
    { F_TANH,  "tanh",  tanh },
    { F_ASINH, "asinh", asinh },
    { F_ACOSH, "acosh", acosh },
    { F_ATANH, "atanh", atanh },
    { F_LOG,   "log",   log },
    { F_LOG10, "log10", log10 },
    { F_LOG2,  "log2",  log2 },
    { F_EXP,   "exp",   exp },
    { F_SQRT,  "sqrt",  sqrt },
    { F_GAMMA,    "gammafun", gammafun },
    { F_LNGAMMA,  "lngamma",  lngamma },
    { F_DIGAMMA,  "digamma",  digamma },
    { F_TRIGAMMA, "trigamma", trigamma },
    { F_INVMILLS, "invmills", invmills },
    { F_ROUND,    "round",    gretl_round },
    { F_SGN,      "sgn",      gretl_sgn },
    { F_CNORM, "cnorm", normal_cdf },
    { F_DNORM, "dnorm", normal_pdf },
    { F_QNORM, "qnorm", normal_cdf_inverse },
    { F_LOGISTIC, "logistic", logistic_cdf },
    { F_REAL,  "Re",    creal },
    { F_IMAG,  "Im",    cimag },
    { F_CARG,  "carg",  carg },
    { F_CMOD,  "cmod",  cabs },
    { F_CQUAD, "cquad", gretl_cquad },
    { 0, NULL, NULL }
};

struct str_table funcs[] = {
    { F_ATAN2,    "atan2" },
    { F_DIFF,     "diff" },
    { F_LDIFF,    "ldiff" },
    { F_SDIFF,    "sdiff" },
    { F_LLAG,     "lags" },
    { F_HFLAG,    "hflags" },
    { F_DROPCOLL, "dropcoll" },
    { F_TOINT,    "int" },
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
    { F_LRCOVAR,  "lrcovar" },
    { F_FEVD,     "fevd" },
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
    { F_RESAMPLE, "resample" },
    { F_PNOBS,    "pnobs" },     /* per-unit nobs in panels */
    { F_PMIN,     "pmin" },      /* panel min */
    { F_PMAX,     "pmax" },      /* panel max */
    { F_PSUM,     "psum" },      /* panel sum */
    { F_PMEAN,    "pmean" },     /* panel mean */
    { F_PXSUM,    "pxsum" },     /* panel x-sectional sum */
    { F_PXNOBS,   "pxnobs" },    /* panel x-sectional obs count */
    { F_PSD,      "psd" },       /* panel std dev */
    { F_PSHRINK,  "pshrink" },
    { F_PEXPAND,  "pexpand" },
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
    { F_MCNORM,   "mcnormal" },
    { F_SUMC,     "sumc" },
    { F_SUMR,     "sumr" },
    { F_PRODC,    "prodc" },
    { F_PRODR,    "prodr" },
    { F_MEANC,    "meanc" },
    { F_MEANR,    "meanr" },
    { F_ASORT,    "asort" },
    { F_CORRESP,  "corresp" },
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
    { F_HDPROD,   "hdprod" },
    { F_MCOV,     "mcov" },
    { F_MCORR,    "mcorr" },
    { F_MXTAB,    "mxtab" },
    { F_CDEMEAN,  "cdemean" },
    { F_CHOL,     "cholesky" },
    { F_PSDROOT,  "psdroot" },
    { F_INSTRINGS, "instrings" },
    { F_INV,      "inv" },
    { F_INVPD,    "invpd" },
    { F_GINV,     "ginv" },
    { F_DIAG,     "diag" },
    { F_TRANSP,   "transp" },
    { F_CTRANS,   "ctrans" },
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
    { F_EIGEN,    "eigen" },
    { F_EIGGEN,   "eigengen" }, /* legacy */
    { F_CMULT,    "cmult" },    /* legacy */
    { F_CDIV,     "cdiv" },     /* legacy */
    { F_SCHUR,    "schur" },
    { F_EIGSOLVE, "eigsolve" },
    { F_NULLSPC,  "nullspace" },
    { F_PRINCOMP, "princomp" },
    { F_MEXP,     "mexp" },
    { F_MLOG,     "mlog" },
    { F_FDJAC,    "fdjac" },
    { F_BFGSMAX,  "BFGSmax" },
    { F_BFGSCMAX, "BFGScmax" },
    { F_NRMAX,    "NRmax" },
    { F_NUMHESS,  "numhess" },
    { F_OBSNUM,   "obsnum" },
    { F_ISDISCR,  "isdiscrete" },
    { F_ISDUMMY,  "isdummy"},
    { F_TYPEOF,   "typeof" },
    { F_TYPENAME, "typename" },
    { F_EXISTS,   "exists" },
    { F_NELEM,    "nelem" },
    { F_PDF,      "pdf" },
    { F_CDF,      "cdf" },
    { F_INVCDF,   "invcdf" },
    { F_PVAL,     "pvalue" },
    { F_CRIT,     "critical" },
    { F_RANDGEN,  "randgen" },
    { F_MRANDGEN, "mrandgen" },
    { F_RANDGEN1, "randgen1" },
    { F_URCPVAL,  "urcpval" },
    { F_QLRPVAL,  "qlrpval" },
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
    { F_SQUARE,   "square" },
    { F_FILTER,   "filter" },
    { F_KFILTER,  "kfilter" },
    { F_KSMOOTH,  "ksmooth" },
    { F_KDSMOOTH, "kdsmooth" },
    { F_KSIMUL,   "ksimul" },
    { F_KSIMDATA, "ksimdata" },
    { F_TRIMR,    "trimr" },
    { F_GETENV,   "getenv" },
    { F_NGETENV,  "ngetenv" },
    { F_ARGNAME,  "argname" },
    { F_OBSLABEL, "obslabel" },
    { F_READFILE, "readfile" },
    { F_BACKTICK, "grab" },
    { F_STRSTR,   "strstr" },
    { F_INSTRING, "instring" },
    { F_STRSTRIP, "strstrip" },
    { F_STRNCMP,  "strncmp" },
    { F_STRLEN,   "strlen" },
    { F_PRINTF,   "printf" },
    { F_SPRINTF,  "sprintf" },
    { F_SSCANF,   "sscanf" },
    { F_VARNAME,  "varname" },
    { F_VARNAMES, "varnames" },
    { F_VARNUM,   "varnum" },
    { F_TOLOWER,  "tolower" },
    { F_TOUPPER,  "toupper" },
    { F_CNAMESET, "cnameset" },
    { F_RNAMESET, "rnameset" },
    { F_LJUNGBOX, "ljungbox" },
    { F_MSORTBY,  "msortby" },
    { F_LINCOMB,  "lincomb" },
    { F_IMHOF,    "imhof" },
    { F_TOEPSOLV, "toepsolv" },
    { F_RGBMIX,   "rgbmix" },
    { F_DSUM,     "diagcat" },
    { F_XMIN,     "xmin" },
    { F_XMAX,     "xmax" },
    { F_CORRGM,   "corrgm" },
    { F_MCOVG,    "mcovg" },
    { F_FCSTATS,  "fcstats" },
    { F_BESSEL,   "bessel" },
    { F_FRACLAG,  "fraclag" },
    { F_MREV,     "mreverse" },
    { F_DESEAS,   "deseas" },
    { F_TRAMOLIN, "linearize" },
    { F_PERGM,    "pergm" },
    { F_IRR,      "irr" },
    { F_NPV,      "npv" },
    { F_WEEKDAY,  "weekday" },
    { F_KDENSITY, "kdensity" },
    { F_MONTHLEN, "monthlen" },
    { F_EPOCHDAY, "epochday" },
    { F_SETNOTE,  "setnote" },
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
    { F_CNAMEGET, "cnameget" },
    { F_RNAMEGET, "rnameget" },
    { F_RANDINT,  "randint" },
    { F_RANDSTR,  "randstr" },
    { F_NADARWAT, "nadarwat" },   /* Nadaraya-Watson */
    { F_SIMANN,   "simann" },     /* simulated annealing */
    { F_LOESS,    "loess" },
    { F_GHK,      "ghk" },
    { F_HALTON,   "halton" },
    { F_IWISHART, "iwishart" },
    { F_ISNAN,    "isnan" },
    { F_TYPESTR,  "typestr" },
    { F_QUADTAB,  "quadtable" },
    { F_AGGRBY,   "aggregate" },
    { F_REMOVE,   "remove" },
    { F_ISODATE,  "isodate" },
    { F_ISOWEEK,  "isoweek" },
    { F_JULDATE,  "juldate" },
    { F_GETLINE,  "getline" },
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
    { F_BARRIER,   "mpibarrier" },
    { F_EASTER,    "easterday" },
    { F_GENSERIES, "genseries" },
    { F_CURL,      "curl" },
    { F_JSONGET,   "jsonget" },
    { F_JSONGETB,  "jsongetb" },
    { F_XMLGET,    "xmlget" },
    { F_NLINES,    "nlines" },
    { F_KPSSCRIT,  "kpsscrit" },
    { F_ARRAY,     "array" },
    { F_STRVALS,   "strvals" },
    { F_STRINGIFY, "stringify" },
    { F_STRVSORT,  "strvsort" },
    { F_BOOTCI,    "bootci" },
    { F_BOOTPVAL,  "bootpval" },
    { F_SEASONALS, "seasonals" },
    { F_DEFARRAY,  "defarray" },
    { F_DEFBUNDLE, "defbundle" },
    { F_DEFLIST,   "deflist" },
    { F_DEFARGS,   "_" },
    { F_KSETUP,    "ksetup" },
    { F_MWEIGHTS,  "mweights" },
    { F_MGRADIENT, "mgradient" },
    { F_MLINCOMB,  "mlincomb" },
    { F_MIDASMULT, "midasmult" },
    { F_HFDIFF,    "hfdiff" },
    { F_HFLDIFF,   "hfldiff" },
    { F_HFLIST,    "hflist" },
    { F_NMMAX,     "NMmax" },
    { F_GSSMAX,    "GSSmax" },
    { F_CNUMBER,   "cnumber" },
    { F_NAALEN,    "naalen" },
    { F_KMEIER,    "kmeier" },
    { F_NORMTEST,  "normtest" },
    { F_ECDF,      "ecdf" },
    { F_NPCORR,    "npcorr" },
    { F_DAYSPAN,   "dayspan" },
    { F_SMPLSPAN,  "smplspan" },
    { F_SLEEP,     "sleep" },
    { F_GETINFO,   "getinfo" },
    { F_CDUMIFY,   "cdummify" },
    { F_SVM,       "svm" },
    { F_GETKEYS,   "getkeys" },
    { F_FEVAL,     "feval" },
    { F_FEVALB,    "fevalb" },
    { F_BINPERMS,  "binperms" },
    { F_BRENAME,   "brename" },
    { F_CCODE,     "isocountry" },
    { F_LSOLVE,    "Lsolve" },
    { F_HYP2F1,    "hyp2f1" },
    { F_STRFTIME,  "strftime" },
    { F_STRFDAY,   "strfday" },
    { F_STRPTIME,  "strptime" },
    { F_STRPDAY,   "strpday" },
    { F_BKW,       "bkw" },
    { F_FZERO,     "fzero" },
    { F_CONV2D,    "conv2d" },
    { F_MSPLITBY,  "msplitby" },
    { F_FLATTEN,   "flatten" },
    { F_ERRORIF,   "errorif" },
    { F_ISCMPLX,   "iscomplex" },
    { F_COMPLEX,   "complex" },
    { F_CONJ,      "conj" },
    { F_CSWITCH,   "cswitch" },
    { F_RANDPERM,  "randperm" },
    { F_STDIZE,    "stdize" },
    { F_STACK,     "stack" },
    { F_GEOPLOT,   "geoplot" },
    { F_BINCOEFF,  "bincoeff" },
    { F_TDISAGG,   "tdisagg" },
    { F_ASSERT,    "assert" },
    { F_VMA,       "vma" },
    { F_BCHECK,    "bcheck" },
    { F_CONTAINS,  "contains" },
    { F_LPSOLVE,   "lpsolve" },
    { F_DISTANCE,  "distance" },
    { F_INTERPOL,  "interpol" },
    { F_MAT2LIST,  "mat2list" },
    { F_DEC2BIN,   "dec2bin" },
    { F_BIN2DEC,   "bin2dec" },
    { F_COMMUTE,   "commute" },
    { F_ACCESS,    "access" },
    { F_SPHCORR,   "sphericorr" },
    { 0,           NULL }
};

struct str_table func_alias[] = {
    { F_EIGEN,    "eiggen2" }, /* deprecated */
    { F_FFT,      "fft2" },    /* deprecated */
    { F_NMMAX,    "NMmin" },
    { F_NRMAX,    "NRmin" },
    { F_BFGSMAX,  "BFGSmin" },
    { F_BFGSCMAX, "BFGScmin" },
    { F_GSSMAX,   "GSSmin" },
    { F_GAMMA,    "gammafunc" },
    { F_GAMMA,    "gamma" },
    { F_LOG,      "logs" },
    { F_LOG,      "ln" },
    { F_SQUARE,   "xpx" },
    { F_BACKTICK, "$" },
    { F_CNAMESET, "colnames" },
    { F_RNAMESET, "rownames" },
    { F_CNAMEGET, "colname" },
    { F_RNAMEGET, "rowname" },
    { 0,          NULL }
};

struct str_table hidden_funcs[] = {
    { HF_JBTERMS,  "_jbterms" },
    { HF_LISTINFO, "_listinfo" },
    { HF_REGLS,    "_regls" },
    { HF_FELOGITR, "_felogit_rec" },
    { HF_FDEPTH,   "_fdepth" },
    { HF_GLASSO,   "_glasso" },
    { 0,           NULL }
};

int const_lookup (const char *s)
{
    int i;

    for (i=0; consts[i].id != 0; i++) {
	if (!strcmp(s, consts[i].str)) {
	    return consts[i].id;
	}
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

    for (i=0; ptrfuncs[i].str != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) ptrfuncs[i].str, &ptrfuncs[i]);
    }

    for (i=0; funcs[i].str != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) funcs[i].str, &funcs[i]);
    }

    for (i=0; hidden_funcs[i].str != NULL; i++) {
	g_hash_table_insert(ht, (gpointer) hidden_funcs[i].str,
			    &hidden_funcs[i]);
    }

    return ht;
}

static GHashTable *oht;

int install_function_override (const char *funname,
			       const char *pkgname,
			       gpointer data)
{
    if (funname == NULL) {
	/* cleanup signal */
	if (oht != NULL) {
	    g_hash_table_destroy(oht);
	    oht = NULL;
	}
	return 0;
    }

    if (oht == NULL) {
	oht = g_hash_table_new_full(g_str_hash, g_str_equal,
				    g_free, NULL);
    }

    if (oht != NULL) {
	gchar *key = g_strdup_printf("%s::%s", pkgname, funname);

	g_hash_table_insert(oht, (gpointer) key, data);
    }

    return 0;
}

int delete_function_override (const char *funname,
			      const char *pkgname)
{
    int ret = 0;

    if (oht != NULL) {
	gchar *key = g_strdup_printf("%s::%s", pkgname, funname);

	if (g_hash_table_remove(oht, key)) {
	    fprintf(stderr, "'%s': deleted override of built-in\n", key);
	    ret = 1;
	}
	g_free(key);
    }

    return ret;
}

static ufunc *get_function_override (const char *sf,
				     gpointer p)
{
    const char *sp = function_package_get_name(p);
    gchar *key = g_strdup_printf("%s::%s", sp, sf);
    ufunc *uf = g_hash_table_lookup(oht, key);

#if 0
    if (uf != NULL) {
	fprintf(stderr, "'%s': using package override\n", key);
    }
#endif
    g_free(key);

    return uf;
}

/* Attention: this function may be called from function_loopkup()
   below, and in that context @p will be NULL -- so don't
   dereference @p without checking it for nullity first!
*/

static int real_function_lookup (const char *s, int aliases,
				 parser *p)
{
    static GHashTable *fht;
    gpointer fnp;

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

    fnp = g_hash_table_lookup(fht, s);
    if (fnp != NULL) {
	struct str_table *st = (struct str_table *) fnp;

	if (p != NULL && st->id > 0 && st->id < FP_MAX) {
	    struct str_table_ex *sx = (struct str_table_ex *) fnp;

	    p->data = sx->ptr;
	}
#if 1
	if (p != NULL) {
	    /* note: point d'appui for deprecation of built-in function */
	    if (st->id == F_EIGGEN) {
		pprintf(p->prn, "*** Warning: %s() is obsolete, please use "
			"eigen() instead ***\n", st->str);
	    } else if (st->id == F_CHOWLIN) {
		pprintf(p->prn, "*** Warning: %s() is obsolete, please use "
			"tdisagg() instead ***\n", st->str);
	    }
	}
#endif
	return st->id;
    }

    if (aliases) {
	int i;

	for (i=0; func_alias[i].id != 0; i++) {
	    if (!strcmp(s, func_alias[i].str)) {
#if 1
		if (!strcmp(s, "fft2")) {
		    gretl_warnmsg_set(_("deprecated alias 'fft2': please call fft()"));
		}
#endif
		if (p != NULL) {
		    p->flags |= P_ALIASED;
		}
		return func_alias[i].id;
	    }
	}
    }

    return 0;
}

void gretl_function_hash_cleanup (void)
{
    real_function_lookup(NULL, 0, NULL);
    install_function_override(NULL, NULL, NULL);
}

int function_lookup (const char *s)
{
    return real_function_lookup(s, 0, NULL);
}

int is_function_alias (const char *s)
{
    int i;

    for (i=0; func_alias[i].id != 0; i++) {
	if (!strcmp(s, func_alias[i].str)) {
	    return 1;
	}
    }

    return 0;
}

void *get_genr_function_pointer (int f)
{
    int i;

    for (i=0; ptrfuncs[i].str != NULL; i++) {
	if (ptrfuncs[i].id == f) {
	    return ptrfuncs[i].ptr;
	}
    }

    return NULL;
}

static int function_lookup_with_alias (const char *s,
				       parser *p)
{
    if (oht != NULL) {
	/* we have a record of one or more package-private
	   functions whose names collide with built-ins
	*/
	gpointer pp = get_active_function_package(OPT_O);

	if (pp != NULL) {
	    ufunc *uf = get_function_override(s, pp);

	    if (uf != NULL) {
		p->idstr = gretl_strdup(s);
		p->data = uf;
		return UFUN;
	    }
	}
    }

    return real_function_lookup(s, 1, p);
}

static const char *funname (int t)
{
    int i;

    for (i=0; ptrfuncs[i].id != 0; i++) {
	if (t == ptrfuncs[i].id) {
	    return ptrfuncs[i].str;
	}
    }

    for (i=0; funcs[i].id != 0; i++) {
	if (t == funcs[i].id) {
	    return funcs[i].str;
	}
    }

    for (i=0; hidden_funcs[i].id != 0; i++) {
	if (t == hidden_funcs[i].id) {
	    return hidden_funcs[i].str;
	}
    }

    return "unknown";
}

static int show_alias (int i)
{
    if (strstr(func_alias[i].str, "min")) {
	return 1;
    } else {
	return 0;
    }
}

/* return the number of built-in functions */

int gen_func_count (void)
{
    int i, n = 0;

    for (i=0; ptrfuncs[i].id != 0; i++) {
	n++;
    }

    for (i=0; funcs[i].id != 0; i++) {
	n++;
    }

    for (i=0; func_alias[i].id != 0; i++) {
	if (show_alias(i)) {
	    n++;
	}
    }

    return n;
}

/* return the name of function @i, including aliases */

const char *gen_func_name (int i)
{
    int j, seq = -1;

    for (j=0; ptrfuncs[j].id != 0; j++) {
	seq++;
	if (seq == i) {
	    return ptrfuncs[i].str;
	}
    }

    for (j=0; funcs[j].id != 0; j++) {
	seq++;
	if (seq == i) {
	    return funcs[j].str;
	}
    }

    for (j=0; func_alias[j].id != 0; j++) {
	if (show_alias(j)) {
	    seq++;
	}
	if (seq == i) {
	    return func_alias[j].str;
	}
    }

    return NULL;
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

int bundle_var_count (void)
{
    int i;

    for (i=0; bvars[i].id != 0; i++) ;
    return i;
}

const char *bundle_var_name (int i)
{
    return bvars[i].str;
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

    for (i=0; ptrfuncs[i].str != NULL; i++) {
	if (!strncmp(s, ptrfuncs[i].str, n)) {
	    return ptrfuncs[i].str;
	}
    }

    for (i=0; funcs[i].str != NULL; i++) {
	if (!strncmp(s, funcs[i].str, n)) {
	    return funcs[i].str;
	}
    }

    return NULL;
}

int gretl_const_count (void)
{
    int i;

    for (i=0; consts[i].id != 0; i++) ;
    return i;
}

const char *gretl_const_name (int i)
{
    return consts[i].str;
}

/* end external stuff */

static int dummy_lookup (const char *s, parser *p)
{
    int i, d = 0;

    for (i=0; dummies[i].id != 0; i++) {
	if (!strcmp(s, dummies[i].str)) {
	    d = dummies[i].id;
	    break;
	}
    }

    if (d == DUM_END) {
	; /* let this pass: we'll check its validity later */
    } else if (MSEL_DUM(d) && parser_next_char(p) != ']') {
	/* most MSEL dummies are stand-alone; they can
	   be valid only if followed by ']'
	*/
	d = 0;
    } else if (d == DUM_EMPTY && strcmp(p->rhs, "empty")) {
        /* if "empty" is not the only term on the right-hand
           side, convert it to "null"
        */
	p->sym = EMPTY;
    }

    return d;
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

int mvar_lookup (const char *s)
{
    int i;

    for (i=0; mvars[i].id != 0; i++) {
	if (!strcmp(s, mvars[i].str)) {
	    return mvars[i].id;
	}
    }

    /* aliases */

    if (!strcmp(s, "$nrsq")) {
	return M_TRSQ;
    } else if (!strcmp(s, "$fcerr")) {
	return M_FCSE;
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

int bvar_lookup (const char *s)
{
    int i;

    for (i=0; bvars[i].id != 0; i++) {
	if (!strcmp(s, bvars[i].str)) {
	    return bvars[i].id;
	}
    }

    return 0;
}

const char *bvarname (int t)
{
    int i;

    for (i=0; bvars[i].id != 0; i++) {
	if (t == bvars[i].id) {
	    return bvars[i].str;
	}
    }

    return "unknown";
}

int genr_function_word (const char *s)
{
    int ret = 0;

    ret = real_function_lookup(s, 0, NULL);
    if (!ret) {
	ret = dvar_lookup(s);
    }
    if (!ret) {
	ret = mvar_lookup(s);
    }
    if (!ret) {
	ret = bvar_lookup(s);
    }
    if (!ret) {
	ret = const_lookup(s);
    }

    return ret;
}

int parser_ensure_error_buffer (parser *p)
{
    if (p->prn == NULL && p->errprn == NULL) {
	p->errprn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
	if (p->errprn != NULL) {
	    p->prn = p->errprn;
	    return 0;
	} else {
	    return E_ALLOC;
	}
    }

    return 0;
}

void undefined_symbol_error (const char *s, parser *p)
{
    parser_ensure_error_buffer(p);
    parser_print_input(p);

    if (p->ch == '.') {
	pprintf(p->prn, _("%s: no such object"), s);
    } else {
	pprintf(p->prn, _("The symbol '%s' is undefined"), s);
    }
    pputc(p->prn, '\n');
    p->err = E_DATA;
}

static void function_noargs_error (const char *s, parser *p)
{
    parser_ensure_error_buffer(p);
    parser_print_input(p);

    pprintf(p->prn, _("'%s': no argument was given"), s);
    pputc(p->prn, '\n');
    p->err = E_ARGS;
}

void context_error (int c, parser *p, const char *func)
{
#if LDEBUG
    if (func != NULL) {
	fprintf(stderr, "context error in %s()\n", func);
    }
#endif
    parser_ensure_error_buffer(p);
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
    } else if (p->sym == EOT) {
	parser_print_input(p);
	pputs(p->prn, _("Incomplete expression\n"));
    } else {
	const char *s = getsymb_full(p->sym, p);

	if (s != NULL && *s != '\0' && strcmp(s, "unknown")) {
	    pprintf(p->prn, _("The symbol '%s' is not valid in this context\n"),
		    getsymb_full(p->sym, p));
	} else {
	    parser_print_input(p);
	}
    }

    if (!p->err) {
	p->err = E_PARSE;
    }
}

/* @parsing_query: we want to keep track of the case where we're
   lexing/parsing the branches of a ternary "query" expression. When
   such an expression is evaluated, it's OK if the branch _not_ taken
   contains an undefined symbol; indeed, this can occur by design, as
   in

     scalar y = exists(x) ? x : 0

   when "x" is in fact undefined.

   We therefore use the "UNDEF" node type to defuse the error that
   would otherwise arise on parsing. An error is triggered only if the
   branch that references the UNDEF node is selected (attempting to
   evaluate an UNDEF node automatically throws an error.)
*/

static int parsing_query;

void set_parsing_query (int s)
{
    parsing_query = s;
}

/* 2023-09-15: apparatus for allowing \" as an embedded quote
   in string literals.
*/

static int alt_double_quote_pos (const char *s, int *esc)
{
    int i, ret = -1;

    for (i=0; s[i]; i++) {
	if (s[i] == '"') {
	    if (i == 0 || s[i-1] != '\\' || s[i+1] == '\0') {
		ret = i;
		break;
	    } else {
		*esc = 1;
	    }
	}
    }

    return ret;
}

static char *escape_quotes (char *s)
{
    int i;

    for (i=1; s[i]; i++) {
	if (s[i] == '"' && s[i-1] == '\\') {
	    shift_string_left(s + i - 1, 1);
	}
    }

    return s;
}

static char *get_quoted_string (parser *p, int prevsym)
{
    char *s = NULL;
    int esc = 0;
    int n;

#if LDEBUG
    fprintf(stderr, "get_quoted_string: sym = '%s', prevsym '%s'\n",
	    getsymb(p->sym), getsymb(prevsym));
    fprintf(stderr, " p->ch = '%c', p->point = '%s'\n", p->ch, p->point);
#endif

    if (prevsym == F_SPRINTF || prevsym == F_PRINTF) {
	/* look for a matching non-escaped double-quote,
	   allowance made for "\\" as itself an escape
	*/
	n = double_quote_position(p->point);
    } else {
	/* look for a matching non-escaped double-quote when
	   backslash is special only when preceding a double
	   quote
	*/
	n = alt_double_quote_pos(p->point, &esc);
    }

    if (n >= 0) {
	s = gretl_strndup(p->point, n);
	if (esc == 1) {
	    /* "\"" is accepted as an escape, but that's all */
	    escape_quotes(s);
	}
	parser_advance(p, n + 1);
    } else {
	parser_print_input(p);
	pprintf(p->prn, _("Unmatched '%c'\n"), '"');
	p->err = E_PARSE;
    }

    if (!p->err) {
	if (p->ch == '.' && *p->point == '$') {
	    /* maybe quoted name of saved model followed by
	       dollar variable? */
	    p->sym = MMEMB;
	} else {
	    p->sym = CSTR;
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

    for (i=0; bvars[i].id != 0; i++) {
	n = strlen(bvars[i].str);
	if (!strncmp(s, bvars[i].str, n)) {
	    return !isalpha(s[n]);
	}
    }

    return 0;
}

static void look_up_dollar_word (const char *s, parser *p)
{
    char *bstr;

    if ((p->idnum = dvar_lookup(s)) > 0) {
	p->sym = DVAR;
    } else if ((p->idnum = const_lookup(s)) > 0) {
	p->sym = CON;
    } else if ((p->idnum = mvar_lookup(s)) > 0) {
	p->sym = MVAR;
    } else if ((p->idnum = bvar_lookup(s)) > 0) {
	p->sym = DBUNDLE;
    } else if ((bstr = get_built_in_string_by_name(s+1))) {
	p->sym = CSTR;
	p->idstr = gretl_strdup(bstr);
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
# include "libset.h"

static int maybe_get_R_function (const char *s)
{
    if (libset_get_bool(R_FUNCTIONS) &&
	strlen(s) >= 3 && !strncmp(s, "R.", 2)) {
	return get_R_function_by_name(s + 2);
    } else {
	return 0;
    }
}

#else /* !USE_RLIB */
# define maybe_get_R_function(s) (0)
#endif

static int doing_genseries;

void set_doing_genseries (int s)
{
    doing_genseries = s;
}

/* Get the next non-space byte beyond what's already parsed:
   this will either be p->ch, or may be found at p->point
   or beyond.
*/

static int parser_next_char (parser *p)
{
    if (p->ch != ' ') {
	return p->ch;
    } else {
	const char *s = p->point;

	while (*s) {
	    if (!isspace(*s)) {
		return *s;
	    }
	    s++;
	}
	return 0;
    }
}

static int char_past_point (parser *p)
{
    if (*p->point != '\0') {
	int i;

	for (i=1; p->point[i] != '\0'; i++) {
	    if (!isspace(p->point[i])) {
		return p->point[i];
	    }
	}
    }
    return 0;
}

/* helpers specific to the stack() function -- dispensable if
   we decide not to support the old-style syntax any more
*/

static int get_oldstyle_stack_args (const char *s, char **arg,
				    char **opt1, char **opt2)
{
    const char *p1, *p2;
    int len;

    /* Let the marker for old-style uage be the presence of
       one or more legacy option flags */

    p1 = strstr(s, "--length=");
    if (p1 != NULL) {
	len = strcspn(p1 + 9, " \n");
	*opt1 = gretl_strndup(p1 + 9, len);
    }

    p2 = strstr(s, "--offset=");
    if (p2 != NULL) {
	len = strcspn(p2 + 9, " \n");
	*opt2 = gretl_strndup(p2 + 9, len);
    }

    if (p1 != NULL || p2 != NULL) {
	len = strcspn(s, ")");
	*arg = gretl_strndup(s, len);
	return 1;
    } else {
	return 0;
    }
}

static int stack_update_parser_input (parser *p)
{
    char *arg = NULL, *opt1 = NULL, *opt2 = NULL;
    char *s, *start;
    gchar *tmp = NULL;
    GString *gs;
    int offset;

    offset = p->point - p->input;
    start = gretl_strndup(p->input, offset);
    gs = g_string_new(start);
    free(start);

    s = strstr(p->input, "stack(") + 6;
    get_oldstyle_stack_args(s, &arg, &opt1, &opt2);
    if (arg != NULL) {
	gs = g_string_append(gs, arg);
	free(arg);
    }
    if (opt1 != NULL) {
	gs = g_string_append_c(gs, ',');
	gs = g_string_append(gs, opt1);
	free(opt1);
    }
    if (opt2 != NULL) {
	gs = g_string_append_c(gs, ',');
	gs = g_string_append(gs, opt2);
	free(opt2);
    }
    gs = g_string_append_c(gs, ')');
    tmp = g_string_free(gs, FALSE);

    p->input = tmp;
    p->point = p->input + offset;
    p->flags |= P_ALTINP;

    return 0;
}

static void handle_lpnext (const char *s, parser *p,
			   int have_dset)
{
    ufunc *u = get_user_function_by_name(s);
    int vnum = -1;

    if (have_dset) {
	vnum = current_series_index(p->dset, s);
    }

    if (u == NULL && vnum >= 0) {
	/* unambiguous: series */
	p->sym = SERIES;
    } else if (u != NULL && vnum < 0) {
	/* unambiguous: function */
	p->sym = UFUN;
    } else if (u != NULL) {
	/* ambiguous case! */
	if (gretl_function_depth() > 0) {
	    /* function writers should avoid collisions
	       when naming series
	    */
	    p->sym = SERIES;
	} else if (p->targ != UNK && p->targ != LIST && p->targ != SERIES) {
	    /* target not compatible with series lag? */
	    p->sym = UFUN;
	} else if (defining_list(p) || p->targ == SERIES) {
	    /* debatable */
	    p->sym = SERIES;
	} else {
	    /* debatable */
	    p->sym = UFUN;
	}
    }

    if (p->sym != 0) {
	p->idstr = gretl_strdup(s);
	if (p->sym == UFUN) {
	    p->data = u;
	} else {
	    p->idnum = vnum;
	    /* in case of any intervening space */
	    while (p->ch == ' ') {
		parser_getc(p);
	    }
	}
    }
}

static int is_function_word (const char *s)
{
    return function_lookup_with_alias(s, NULL) != 0 ||
	get_user_function_by_name(s) != NULL;
}

static void look_up_word (const char *s, parser *p)
{
    int have_dset = (p->dset != NULL && p->dset->v > 0);
    int prevsym = p->sym;
    int lpnext, err = 0;

#if LDEBUG
    fprintf(stderr, "look_up_word: s='%s', ch='%c', next='%c'\n",
	    s, p->ch, parser_next_char(p));
#endif

    /* is the next (or next non-space) character left paren? */
    lpnext = parser_next_char(p) == '(';

    /* initialize */
    p->sym = 0;
    p->data = NULL;

    if (lpnext) {
	/* identifier is immediately followed by left paren:
	   most likely a function call but could be the name
	   of a series followed by a lag specifier
	*/
	p->sym = function_lookup_with_alias(s, p);
	if (p->sym == 0) {
	    handle_lpnext(s, p, have_dset);
	} else if (p->sym == F_STACK) {
	    /* special! */
	    if (strstr(p->point, "--length") || strstr(p->point, "--offset")) {
		stack_update_parser_input(p);
	    }
	    p->flags |= P_STACK;
	}
    }

    if (p->sym == 0) {
	p->idnum = const_lookup(s);
	if (p->idnum > 0) {
	    p->sym = CON;
	} else {
	    p->idnum = dummy_lookup(s, p);
	    if (p->idnum > 0) {
		p->sym = DUM;
	    } else if (have_dset &&
		       (p->idnum = current_series_index(p->dset, s)) >= 0) {
		p->sym = SERIES;
		p->idstr = gretl_strdup(s);
	    } else if (have_dset && !strcmp(s, "time")) {
		p->sym = DUM;
		p->idnum = DUM_TREND;
	    } else if ((p->data = get_user_var_by_name(s)) != NULL) {
		user_var *u = p->data;

		if (u->type == GRETL_TYPE_DOUBLE) {
		    p->sym = NUM;
		} else if (u->type == GRETL_TYPE_MATRIX) {
		    p->sym = MAT;
		} else if (u->type == GRETL_TYPE_BUNDLE) {
		    p->sym = BUNDLE;
		} else if (u->type == GRETL_TYPE_ARRAY) {
		    p->sym = ARRAY;
		} else if (u->type == GRETL_TYPE_STRING) {
		    p->sym = STR;
		} else if (u->type == GRETL_TYPE_LIST) {
		    p->sym = LIST;
		}
		p->idstr = gretl_strdup(s);
	    } else if (defining_list(p) && varname_match_any(p->dset, s)) {
		p->sym = WLIST;
		p->idstr = gretl_strdup(s);
	    } else if (have_dset && !strcmp(s, "t")) {
		/* if "t" has not been otherwise defined, treat it
		   as an alias for "obs"
		*/
		p->sym = DVAR;
		p->idnum = R_INDEX;
	    } else if (maybe_get_R_function(s)) {
		/* note: all "native" types take precedence over this */
		p->sym = RFUN;
		p->idstr = gretl_strdup(s + 2);
	    } else if (parsing_query || prevsym == B_AND || prevsym == B_OR) {
		p->sym = UNDEF;
		p->idstr = gretl_strdup(s);
	    } else if (p->flags & (P_AND | P_OR)) {
		p->sym = UNDEF;
		p->idstr = gretl_strdup(s);
	    } else if (gretl_get_object_by_name(s)) {
		p->sym = UOBJ;
		p->idstr = gretl_strdup(s);
	    } else {
		err = E_UNKVAR;
	    }
	}
    }

    if (err) {
	if (is_function_word(s)) {
	    /* @s is a function identifier with no
	       following left paren */
	    function_noargs_error(s, p);
	} else if (object_is_function_arg(s)) {
	    /* @s is the name of a function parameter
	       with a null value */
	    p->sym = NULLARG;
	    p->idstr = gretl_strdup(s);
	} else if (p->sym != DUM_EMPTY) {
	    undefined_symbol_error(s, p);
	}
    }

#if LDEBUG
    fprintf(stderr, "look_up_word: at return err=%d, p->err=%d\n", err, p->err);
#endif
}

static void maybe_treat_as_postfix (parser *p)
{
    if (p->sym == NUM) {
	const char *ok = ")]}+-*/%,:";
	int c = char_past_point(p);

	/* Interpret as foo++ or foo-- ? Only if
	   the following character is suitable.
	*/
	if (c == 0 || strchr(ok, c)) {
	    p->sym = p->ch == '+'? NUM_P : NUM_M;
	    /* swallow the pluses or minuses */
	    parser_advance(p, 1);
	}
    }
}

#define dollar_series(t) (t > R_SCALAR_MAX && t < R_SERIES_MAX)

#define could_be_matrix(t) (model_data_matrix(t) || \
			    model_data_matrix_builder(t) || \
			    t == M_UHAT || t == M_YHAT || \
			    (t > R_SERIES_MAX && t < R_MAX))

#define could_be_bundle(t) (t == R_RESULT)

static void word_check_next_char (parser *p)
{
    char chk[2] = {p->ch, '\0'};

#if LDEBUG
    if (p->ch) {
	fprintf(stderr, "word_check_next_char: ch='%c', sym='%s'\n",
		p->ch, getsymb(p->sym));
    } else {
	fprintf(stderr, "word_check_next_char: ch = NUL\n");
    }
#endif
    p->upsym = 0;

    /* Here we're checking for validity and syntactical implications
       of one of "([.+-" immediately following a 'word' of some kind.
    */
    if (strspn(chk, "([.+-") == 0) {
	/* none of the above */
	return;
    }

    if (p->sym == UNDEF) {
	/* 2020-03-16: suspend disbelief, we might be in a branch
	   of "cond ? x : y" that ends up not being taken */
	return;
    }

    if (p->ch == '(') {
	/* series or list (lag), or function */
	if (p->sym == SERIES) {
	    if (p->idnum > 0 && p->idnum == p->lh.vnum) {
		p->flags |= P_AUTOREG;
	    }
	    p->upsym = p->sym;
	    p->sym = LAG;
	} else if (p->sym == LIST) {
	    p->upsym = p->sym;
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
	p->upsym = p->sym;
	if (p->sym == MAT || p->sym == ARRAY ||
	    p->sym == LIST || p->sym == STR) {
	    /* slice of sliceable object */
	    p->sym = OSL;
	} else if (*p->point != '"' &&
		   (p->sym == MVAR || p->sym == DVAR) &&
		   could_be_matrix(p->idnum)) {
	    /* slice of $-matrix? */
	    p->sym = OSL;
	} else if (p->sym == SERIES) {
	    /* observation from series */
	    p->sym = OBS;
	} else if (p->sym == DVAR && dollar_series(p->idnum)) {
	    /* observation from "dollar" series */
	    p->sym = OBS;
	} else if (p->sym == MVAR && model_data_list(p->idnum)) {
	    /* element/range of accessor list */
	    p->sym = OSL;
	} else if (p->sym == BUNDLE) {
	    /* member from bundle */
	    p->sym = BMEMB;
	} else if (p->sym == DBUNDLE) {
	    /* member from $ bundle */
	    p->sym = DBMEMB;
	} else if (p->sym == DVAR && could_be_bundle(p->idnum)) {
	    p->sym = DBMEMB;
	} else {
	    p->err = E_PARSE;
	}
    } else if (p->ch == '.' && *p->point == '$') {
	if (p->sym == UOBJ) {
	    /* name of saved model followed by dollar variable? */
	    p->sym = MMEMB;
	} else if (p->sym == CSTR) {
	    /* maybe quoted name of saved object followed by
	       dollar variable? */
	    p->sym = MMEMB;
	} else {
	    p->err = E_PARSE;
	}
    } else if (p->ch == '.' && isalpha(*p->point)) {
	if (p->sym == LIST) {
	    p->sym = LISTVAR;
	} else if (p->sym == BUNDLE) {
	    p->sym = BMEMB;
	} else if (p->sym == DBUNDLE) {
	    p->sym = DBMEMB;
	} else if (p->sym == DVAR && could_be_bundle(p->idnum)) {
	    p->sym = DBMEMB;
	} else {
	    p->err = E_PARSE;
	}
    } else if (p->ch == '+' && *p->point == '+') {
	maybe_treat_as_postfix(p);
    } else if (p->ch == '-' && *p->point == '-') {
	maybe_treat_as_postfix(p);
    }

    if (p->err) {
	context_error(p->ch, p, "word_check_next_char");
    }
}

static int is_word_char (parser *p)
{
    if (strchr(wordchars, p->ch) != NULL) {
	return 1;
    } else if (defining_list(p) && !doing_genseries &&
	       (p->ch == '*' || p->ch == '?')) {
	return 1;
    }

    return 0;
}

static void getword (parser *p, int greek)
{
    char word[32];
    int i = 0;

    if (greek) {
	/* we have a single (2-byte) UTF-8 greek letter in scope */
	for (i=0; i<2; i++) {
	    word[i] = p->ch;
	    parser_getc(p);
	}
    } else {
	/* we know the first char is acceptable (and might be '$' or '_') */
	word[i++] = p->ch;
	parser_getc(p);

#ifdef USE_RLIB
	/* allow for R.foo function namespace */
	if (*word == 'R' && p->ch == '.' && *p->point != '$') {
	    if (libset_get_bool(R_FUNCTIONS) && !gretl_is_bundle("R")) {
		word[i++] = p->ch;
		parser_getc(p);
	    }
	}
#endif
	while (p->ch != 0 && is_word_char(p) && i < 31) {
	    word[i++] = p->ch;
	    parser_getc(p);
	}
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
	p->sym = CSTR;
	p->idstr = gretl_strdup(word);
	p->flags ^= P_GETSTR;
	return; /* FIXME bundle-member name */
    } else if ((*word == '$' && word[1]) || !strcmp(word, "obs")) {
	look_up_dollar_word(word, p);
    } else if (*word == '$' && word[1] == '\0' && (p->ch == '[' || p->ch == '.')) {
	p->sym = DBUNDLE;
	p->idnum = B_MODEL;
    } else {
	look_up_word(word, p);
    }

    if (!p->err) {
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

/* Below: we're testing p->ch for validity, given what we've already
   packed into string @s up to element @i, and with some regard to
   the next character in line.
*/

static int ok_dbl_char (parser *p, char *s, int i)
{
    int ret = 0;

    if (i < 0 || (p->ch >= '0' && p->ch <= '9')) {
	return 1;
    }

    switch (p->ch) {
    case '+':
    case '-':
	ret = s[i] == 'e' || s[i] == 'E';
	break;
    case '.':
	ret = !strchr(s, '.') && !strchr(s, ':') &&
	    !strchr(s, 'e') && !strchr(s, 'E') &&
	    *p->point != '.';
	break;
    case 'e':
    case 'E':
	ret = !strchr(s, 'e') && !strchr(s, 'E') &&
	    !strchr(s, ':');
	break;
    case ':':
	/* allow for obs numbers in the form, e.g., "1995:10" */
	ret = colon_ok(p, s, i);
	break;
    default:
	break;
    }

    return ret;
}

static int ok_hex_char (parser *p, char *s, int i)
{
    if (i < 1) {
	return 1;
    } else if (p->ch >= '0' && p->ch <= '9') {
	return 1;
    } else if (p->ch >= 'a' && p->ch <= 'f') {
	return 1;
    } else if (p->ch >= 'A' && p->ch <= 'F') {
	return 1;
    }

    return 0;
}

static void parse_number (parser *p)
{
    char xstr[NUMLEN] = {0};
    int got_colon = 0;
    int got_hex = 0;
    int i = 0;

    if (p->ch == '0' && (p->point[0] == 'x' || p->point[0] == 'X')) {
	/* hexadecimal input */
	got_hex = 1;
	while (ok_hex_char(p, xstr, i - 1) && i < HEXLEN - 1) {
	    xstr[i++] = p->ch;
	    parser_getc(p);
	}
    } else {
	/* decimal input */
	while (ok_dbl_char(p, xstr, i - 1) && i < NUMLEN - 1) {
	    xstr[i++] = p->ch;
	    if (p->ch == ':') {
		got_colon = 1;
	    }
	    parser_getc(p);
	}
	while (p->ch >= '0' && p->ch <= '9') {
	    /* flush excess numeric characters */
	    parser_getc(p);
	}
    }

#if LDEBUG
    fprintf(stderr, "parse_number: xstr = '%s'\n", xstr);
#endif

    if (got_colon) {
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
	    p->sym = CSTR;
	}
    } else if (got_hex) {
	guint32 u = strtoul(xstr, NULL, 16);

	p->xval = (double) u;
	p->sym = CNUM;
    } else {
	p->xval = dot_atof(xstr);
	p->sym = CNUM;
#if LDEBUG
	fprintf(stderr, " dot_atof gave %g\n", p->xval);
#endif
    }
}

static int wildcard_special (parser *p)
{
    char cprev = *(p->point - 2);

    if (p->ch == '?') {
	char cnext = *p->point;

	if ((cprev == ' ' || cprev == ')') && cnext == ' ') {
	    /* '?' is presumably ternary operator */
	    return 0;
	}
    }

    if (cprev == ' ' &&
	(bare_data_type(p->sym) || closing_sym(p->sym) ||
	 (p->sym == LAG))) {
	p->sym = B_LCAT;
    } else {
	getword(p, 0);
    }

    return 1;
}

/* 0xE2 0x88 0x92 = UTF-8 minus */

static void lex_try_utf8 (parser *p)
{
    if ((unsigned char) *p->point == 0x88 &&
	(unsigned char) *(p->point + 1) == 0x92) {
	p->sym = B_SUB;
	parser_getc(p);
	parser_getc(p);
	parser_getc(p);
    } else {
	pprintf(p->prn, _("Unexpected byte 0x%x\n"),
		(unsigned char) p->ch);
	p->err = E_PARSE;
    }
}

#define word_start_special(c) (c == '$' || c == '_')

/* accept 'to', but only with spaces before and after */
#define lag_range_sym(p) ((p->flags & P_LAGPRSE) && p->ch == 't' && \
                          *p->point == 'o' && \
			  *(p->point - 2) == ' ' && \
			  *(p->point + 1) == ' ')

void lex (parser *p)
{
    static int prevsyms[2];

#if LDEBUG
    if (p->ch) {
	fprintf(stderr, "lex: p->ch='%c', point='%c'\n", p->ch, *p->point);
    } else {
	fprintf(stderr, "lex: p->ch is NUL\n");
    }
#endif
    prevsyms[0] = prevsyms[1];
    prevsyms[1] = p->sym;

    if (p->ch == 0) {
	p->sym = EOT;
	return;
    }

    while (p->ch != 0) {
	if ((unsigned char) p->ch == 0xE2) {
	    lex_try_utf8(p);
	    return;
	} else if (is_greek_letter(p->point - 1)) {
	    getword(p, 1);
	    return;
	}
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
	    if (defining_list(p) && !doing_genseries) {
		/* allow for '*' as wildcard */
		wildcard_special(p);
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
		parser_getc(p);
		p->sym = B_EQ;
	    } else {
		gretl_errmsg_set(_("If you meant to test for "
				   "equality, please use '=='"));
		p->err = E_PARSE;
	    }
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
	    if (defining_list(p)) {
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
	    if (defining_list(p) && !doing_genseries) {
		/* allow for '?' as wildcard */
		if (wildcard_special(p)) {
		    return;
		}
	    }
	    p->sym = QUERY;
	    parser_getc(p);
	    return;
	case '.':
	    if (*p->point == '$') {
		p->sym = P_DOT;
		parser_getc(p);
		return;
	    } else if (isalpha(*p->point)) {
		/* 2017-01-07 */
		p->sym = BMEMB;
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
	    } else if (p->ch == '!' && *p->point == '=') {
		p->sym = B_DOTNEQ;
		parser_getc(p);
		parser_getc(p);
		return;
	    } else if (p->ch == '.') {
		p->sym = B_ELLIP;
		parser_getc(p);
		return;
	    } else {
		/* not a "dot operator", so back up */
		parser_ungetc(p);
	    }
	    /* Falls through. */
        default:
	    if (defining_list(p) && lag_range_sym(p)) {
		p->sym = B_RANGE;
		parser_getc(p);
		parser_getc(p);
		return;
	    }
	    if (defining_list(p) && !doing_genseries &&
		(bare_data_type(p->sym) || closing_sym(p->sym) ||
		 p->sym == LAG) && *(p->point - 2) == ' ') {
		/* may be forming a list, but only if there are
		   spaces between the terms
		*/
		p->sym = B_LCAT;
		return;
	    }
	    if (isdigit(p->ch)) {
		parse_number(p);
		return;
	    } else if (p->ch == '.' && isdigit(*p->point)) {
		parse_number(p);
		return;
	    } else if (islower(p->ch) || isupper(p->ch) ||
		       word_start_special(p->ch)) {
		getword(p, 0);
		return;
	    } else if (p->ch == '"') {
		p->idstr = get_quoted_string(p, prevsyms[0]);
		return;
	    } else if (p->ch == '#' && prevsyms[0] == CSTR) {
		/* 2023-12-22: an uncaught inline comment? */
		p->ch = '\0';
		p->sym = EOT;
		return;
	    } else {
		parser_print_input(p);
		if (isprint(p->ch)) {
		    pprintf(p->prn, _("Invalid character '%c'\n"), p->ch);
		} else {
		    pprintf(p->prn, _("Unexpected byte 0x%x\n"),
			    (unsigned char) p->ch);
		}
		p->err = E_PARSE;
		return;
	    }
	} /* end ch switch */
    } /* end while ch != 0 */
}

const char *getsymb_full (int t, const parser *p)
{
    if (t == F_DEFMAT) {
	/* pseudo function */
	return "DEFMAT";
    } else if ((t > F1_MIN && t < F1_MAX) ||
	(t > F1_MAX && t < F2_MAX) ||
	(t > F2_MAX && t < FN_MAX)) {
	return funname(t);
    }

    if (t == EOT) {
	return "EOT";
    }

    /* yes, well */
    if (t == OBS) {
	return "OBS";
    } else if (t == OSL) {
	return "OSL";
    } else if (t == SUB_ADDR) {
	return "SUB_ADDR";
    } else if (t == DMSTR) {
	return "DMSTR";
    } else if (t == SLRAW) {
	return "SLRAW";
    } else if (t == MSPEC) {
	return "MSPEC";
    } else if (t == SUBSL) {
	return "SUBSL";
    } else if (t == FARGS) {
	return "FARGS";
    } else if (t == LIST || t == WLIST) {
	return "LIST";
    } else if (t == EMPTY) {
	return "EMPTY";
    } else if (t == LISTVAR) {
	return "LISTVAR";
    } else if (t == BMEMB) {
	return "BMEMB";
    } else if (t == SERIES) {
	return "SERIES";
    } else if (t == MAT) {
	return "MAT";
    } else if (t == UNDEF) {
	return "UNDEF";
    } else if (t == NUM) {
	return "NUM";
    } else if (t == CNUM) {
	return "CNUM";
    } else if (t == IVEC) {
	return "IVEC";
    } else if (t == NUM_P) {
	return "NUM_P";
    } else if (t == NUM_M) {
	return "NUM_M";
    } else if (t == DBUNDLE) {
	return "DBUNDLE";
    } else if (t == DBMEMB) {
	return "DBMEMB";
    } else if (t == MMEMB) {
	return "MMEMB";
    } else if (t == DUM_EMPTY) {
	return "empty";
    }

    if (p != NULL) {
	if (t == BUNDLE) {
	    return p->idstr;
	} else if (t == ARRAY) {
	    return p->idstr;
	} else if (t == UOBJ) {
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
	} else if (t == STR || t == CSTR) {
	    return p->idstr;
	}
    } else {
	if (t == BUNDLE) {
	    return "BUNDLE";
	} else if (t == ARRAY) {
	    return "ARRAY";
	} else if (t == UOBJ) {
	    return "UOBJ";
	} else if (t == CON) {
	    return "CON";
	} else if (t == DUM) {
	    return "dummy constant";
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
	} else if (t == CSTR) {
	    return "CSTR";
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
	return "==";
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
    case B_JOIN:
	return "JOIN";
    case B_RANGE:
	return " to ";
    case B_ELLIP:
	return "..";
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
    case B_DOTNEQ:
	return ".!=";
    case B_DOTASN:
	return "dot-assign";
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

const char *getsymb (int t)
{
    return getsymb_full(t, NULL);
}
