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

/* shared private header for all 'genr' related modules */

#ifndef GENPARSE_H
#define GENPARSE_H

#include "libgretl.h"
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_bundle.h"
#include "gretl_array.h"

#define GENDEBUG 0

/* operators, types, punctuation */

enum {
              U_NEG = 1,
              U_POS,
              U_NOT,
              U_ADDR,
              U_MAX,      /* SEPARATOR: end of unary operators */
              B_ASN,
              B_ADD,
              B_SUB,
              B_MUL,
  /* 10 */    B_DIV,
              B_MOD,
              B_POW,
              B_EQ,
              B_LT,
              B_GT,
              B_LTE,
              B_GTE,
              B_NEQ,
              B_AND,
  /* 20 */    B_OR,
              B_TRMUL,
	      B_RANGE,
              B_DOTMULT,
	      B_DOTDIV,
	      B_DOTPOW,
              B_DOTADD,
              B_DOTSUB,
              B_DOTEQ,
              B_DOTLT,
  /* 30 */    B_DOTGT,
	      B_DOTLTE,
	      B_DOTGTE,
	      B_DOTNEQ,
              B_DOTASN,
              B_KRON,     /* Kronecker product */
              B_HCAT,     /* horizontal concatenation */
              B_VCAT,     /* vertical concatenation */
	      B_LCAT,     /* list concatentation */
	      B_LDIV,     /* matrix left division */
  /* 40 */    B_ELLIP,    /* list-generating ellipsis */
	      B_JOIN,     /* list-joining with separator */
	      OP_MAX,     /* SEPARATOR: end of binary operators */
              G_LPR,      /* grouping: left paren */
	      G_RPR,      /* right paren */
              G_LBR,      /* left bracket */
              G_RBR,      /* right bracket */
              G_LCB,      /* left curly bracket */
	      G_RCB,      /* right curly bracket */
	      P_COM,	  /* punctuation: comma */
  /* 50 */    P_DOT,	  /* period */
	      P_SEMI,	  /* semi-colon */
	      P_COL,	  /* colon */
              PUNCT_MAX,  /* SEPARATOR: end of grouping and punctuation marks */
	      NUM,        /* scalar */
	      SERIES,	  /* series */
	      LIST,       /* list of series */
	      MAT,	  /* matrix */
	      BUNDLE,     /* gretl bundle (hash table) */
	      ARRAY,      /* generic array object */
  /* 60 */    STR,	  /* string */
	      CNUM,	  /* constant (literal) numeric value */
	      CSTR,       /* constant (literal) string */
	      CON,	  /* named numeric constant */
	      DUM,	  /* dummy variable */
	      UOBJ,	  /* user-defined object (e.g. model) */
	      NUM_P,      /* user scalar++ */
	      NUM_M,      /* user scalar-- */
	      OBS,	  /* observation from a series */
	      DMSTR,	  /* "dollar" matrix plus string subspec */
  /* 70 */    SLRAW,	  /* unevaluated "slice" specification */
	      MSPEC,	  /* evaluated matrix subspec */
	      SUBSL,	  /* row or column component of MSPEC */
              LAG,        /* variable plus lag length */
	      DVAR,	  /* $ "dataset" variable (mostly scalar or series) */
	      MVAR,	  /* $ model var (scalar, series, or matrix) */
	      LISTVAR,    /* variable in list, dot syntax */
	      DBUNDLE,    /* $ bundle accessor */
	      BMEMB,      /* member of bundle */
	      DBMEMB,     /* member of $ bundle */
  /* 80 */    MMEMB,      /* member of named model */
	      FARGS,	  /* set of n function arguments */
              WLIST,      /* wildcard list spec */
              EMPTY,      /* "null" or empty arg slot */
	      UNDEF,      /* undefined (allowed in "query" context only) */
	      NULLARG,    /* record of missing function argument */
	      DTYPE_MAX,  /* SEPARATOR: end of "bare" types */
	      UFUN,	  /* user-defined function */
	      RFUN,       /* GNU R function */
	      IVEC,       /* array of ints, not a varlist */
  /* 90 */    OSL,        /* "slice" of object (matrix, array, string) */
              USERIES,    /* named series (defined only for error reporting) */
	      SUB_ADDR,   /* "address" of (e.g.) array element */
              INC,        /* increment */
              DEC,        /* decrement */
	      QUERY,      /* ternary "?" expression */
	      EOT,	  /* end of transmission */
	      UNK
};

/* functions: don't collide with the enumeration above */

enum {
    F1_MIN = 1 << 8,
    F_ABS,
    F_SGN,
    F_CEIL,
    F_FLOOR,
    F_SIN,
    F_COS,
    F_TAN,
    F_ASIN,
    F_ACOS,
    F_ATAN,
    F_SINH,
    F_COSH,
    F_TANH,
    F_ASINH,
    F_ACOSH,
    F_ATANH,
    F_LOG,
    F_LOG10,
    F_LOG2,
    F_EXP,
    F_SQRT,
    F_GAMMA,
    F_LNGAMMA,
    F_DIGAMMA,
    F_TRIGAMMA,
    F_INVMILLS,
    F_ROUND,
    F_CNORM,
    F_DNORM,
    F_QNORM,
    F_LOGISTIC,
    F_REAL,
    F_IMAG,
    F_CARG,
    F_CMOD,
    F_CQUAD,
    FP_MAX,      /* separator: end of pointerized functions */
    F_CONJ,
    F_TOINT,
    F_DIFF,	  /* first difference */
    F_SDIFF,	  /* seasonal difference */
    F_SORT,	  /* ascending sort */
    F_DSORT,	  /* descending sort */
    F_RANKING,
    F_ODEV,	  /* orthogonal deviation */
    F_NOBS,
    F_CUM,
    F_MISSING,
    F_DATAOK,
    F_MISSZERO,
    F_ZEROMISS,
    F_MEDIAN,
    F_GINI,
    F_SUMALL,
    F_SKEWNESS,
    F_KURTOSIS,
    F_SST,
    F_CHOL,
    F_INV,
    F_DIAG,
    F_TRANSP,
    F_VEC,
    F_ROWS,
    F_COLS,
    F_DET,
    F_LDET,
    F_TRACE,
    F_NORM1,
    F_INFNORM,
    F_RCOND,
    F_OBSNUM,
    F_ISDISCR,
    F_ISDUMMY,
    F_TYPEOF,
    F_TYPENAME,
    F_EXISTS,
    F_NELEM,
    F_VALUES,
    F_UNIQ,
    F_NULLSPC,
    F_MEXP,
    F_FFT,
    F_FFTI,
    F_UPPER,
    F_LOWER,
    F_OBSLABEL,
    F_BACKTICK,
    F_STRLEN,
    F_VARNAME,
    F_VARNAMES,
    F_VARNUM,
    F_TOLOWER,
    F_TOUPPER,
    F_IRR,
    F_ERRMSG,
    F_GETENV,
    F_NGETENV,
    F_ISNAN,
    F_TYPESTR,
    F_STRSTRIP,
    F_REMOVE,
    F_ATOF,
    F_MPI_RECV,
    F_EASTER,
    F_CURL,
    F_NLINES,
    F_ARRAY,
    F_TRAMOLIN,
    F_CNUMBER,
    F_ECDF,
    F_SLEEP,
    F_GETINFO,
    F_CDUMIFY,
    F_GETKEYS,
    F_MCORR,
    F_ISCMPLX,
    F_CTRANS,
    F_MLOG,
    F_BARRIER,
    F_LPSOLVE,
    F_INTERPOL,
    F_DEC2BIN,
    F_BIN2DEC,
    F_ACCESS,
    F_POLROOTS,
    F_RANDSTR,
    HF_JBTERMS,
    HF_FDEPTH,
    F1_MAX,	  /* SEPARATOR: end of single-arg functions */
    HF_LISTINFO,
    F_MIN,
    F_MAX,
    F_RANK,
    F_GINV,
    F_SUM,
    F_MEAN,
    F_VCE,
    F_SD,
    F_ARGNAME,
    F_T1,
    F_T2,
    F_COV,
    F_MCOV,
    F_DUMIFY,
    F_SORTBY,
    F_RUNIFORM,
    F_RNORMAL,
    F_FRACDIFF,
    F_BOXCOX,
    F_ZEROS,
    F_ONES,
    F_MUNIF,
    F_QFORM,
    F_EIGSYM,
    F_QUANTILE,
    F_HDPROD,     /* horizontal direct product */
    F_MXTAB,
    F_MRSEL,
    F_MCSEL,
    F_CNAMESET,
    F_RNAMESET,
    F_LJUNGBOX,
    F_MSORTBY,
    F_LINCOMB,
    F_IMHOF,
    F_FRACLAG,
    F_MREV,
    F_PERGM,
    F_NPV,
    F_DSUM,
    F_POLYFIT,
    F_INLIST,
    F_ISCONST,
    F_INBUNDLE,
    F_CNAMEGET,
    F_RNAMEGET,
    F_PNOBS,
    F_PMIN,
    F_PMAX,
    F_PSUM,
    F_PMEAN,
    F_PXSUM,
    F_PXNOBS,
    F_PSD,
    F_PSHRINK,
    F_RANDINT,
    F_MREAD,
    F_BREAD,
    F_GETLINE,
    F_ISODATE,
    F_JULDATE,
    F_READFILE,
    F_PRINTF,
    F_SPRINTF,
    F_MPI_SEND,
    F_BCAST,
    F_ALLREDUCE,
    F_GENSERIES,
    F_KPSSCRIT,
    F_STRINGIFY,
    F_SQUARE,
    F_SEASONALS,
    F_DROPCOLL,
    F_KSIMDATA,
    F_HFDIFF,
    F_HFLDIFF,
    F_NAALEN,
    F_KMEIER,
    F_NORMTEST,
    F_COR,
    F_LRCOVAR,
    F_JSONGETB,
    F_FIXNAME,
    F_ATAN2,
    F_CCODE,
    F_LSOLVE,
    F_STRFDAY,
    F_STRPTIME,
    F_STRPDAY,
    F_CONV2D,
    F_FLATTEN,
    F_IMAT,
    F_COMPLEX,
    F_RANDPERM,
    F_CSWITCH,
    F_PSDROOT,
    F_STRVALS,
    F_ERRORIF,
    F_BINCOEFF,
    F_ASSERT,
    F_CONTAINS,
    F_VECH,
    F_UNVECH,
    F_DESEAS,
    F_MAT2LIST,
    F_CMULT,
    F_CDIV,
    F_PEXPAND,
    F_INVPD,
    F_MINC,
    F_MAXC,
    F_MINR,
    F_MAXR,
    F_IMINC,
    F_IMAXC,
    F_IMINR,
    F_IMAXR,
    F_SUMC,
    F_SUMR,
    F_PRODC,
    F_PRODR,
    F_MEANC,
    F_MEANR,
    F_ASORT,
    F_CORRESP,
    F_STRVSORT,
    F_FEVALB,
    F_BINPERMS,
    F_MNORM,
    HF_VCNORM,
    HF_GLASSO,
    F2_MAX,	  /* SEPARATOR: end of two-arg functions */
    F_WMEAN,
    F_WVAR,
    F_WSD,
    F_LLAG,
    F_HFLAG,
    F_PRINCOMP,
    F_BFGSMAX,
    F_MSHAPE,
    F_SVD,
    F_QR,
    F_TRIMR,
    F_CORRGM,
    F_SEQ,
    F_REPLACE,
    F_STRNCMP,
    F_BESSEL,
    F_WEEKDAY,
    F_MONTHLEN,
    F_EPOCHDAY,
    F_KDENSITY,
    F_SETNOTE,
    F_BWFILT,
    F_VARSIMUL,
    F_STRSUB,
    F_REGSUB,
    F_MLAG,
    F_EIGSOLVE,
    F_SIMANN,
    F_HALTON,
    F_MWRITE,
    F_BWRITE,
    F_AGGRBY,
    F_IWISHART,
    F_SSCANF,
    F_SUBSTR,
    F_REDUCE,
    F_SCATTER,
    F_MWEIGHTS,
    F_MGRADIENT,
    F_MLINCOMB,
    F_HFLIST,
    F_NMMAX,
    F_GSSMAX,
    F_NPCORR,
    F_DAYSPAN,
    F_SMPLSPAN,
    F_FDJAC,
    F_NUMHESS,
    F_STRSPLIT,
    F_HPFILT,
    F_XMLGET,
    F_JSONGET,
    F_FEVD,
    F_LRVAR,
    F_BRENAME,
    F_ISOWEEK,
    F_BKW,
    F_FZERO,
    F_EIGEN,
    F_EIGGEN, /* legacy */
    F_SCHUR,
    F_RESAMPLE,
    F_STACK,
    F_GEOPLOT,
    F_VMA,
    F_FCSTATS,
    F_BCHECK,
    F_MSPLITBY,
    F_DISTANCE,
    F_SPHCORR,
    F_STRSTR,
    F_INSTRING,
    F_STRFTIME,
    F_LDIFF,
    F_SDC,
    F_CDEMEAN,
    F_STDIZE,
    F_INSTRINGS,
    HF_REGLS,
    F3_MAX,       /* SEPARATOR: end of three-arg functions */
    F_URCPVAL,
    F_RANDGEN,
    F_MRANDGEN,
    F_RANDGEN1,
    F_PDF,
    F_PVAL,
    F_CDF,
    F_INVCDF,
    F_CRIT,
    F_BKFILT,
    F_MOLS,
    F_MPOLS,
    F_MRLS,
    F_FILTER,
    F_MCOVG,
    F_KFILTER,
    F_KSMOOTH,
    F_KDSMOOTH,
    F_KSIMUL,
    F_NRMAX,
    F_LOESS,
    F_GHK,
    F_QUADTAB,
    F_ISOCONV,
    F_QLRPVAL,
    F_BOOTCI,
    F_BOOTPVAL,
    F_MOVAVG,
    F_DEFARRAY,
    F_DEFBUNDLE,
    F_DEFLIST,
    F_DEFARGS,
    F_DEFMAT,
    F_KSETUP,
    F_BFGSCMAX,
    F_SVM,
    F_IRF,
    F_NADARWAT,
    F_FEVAL,
    F_CHOWLIN,
    F_TDISAGG,
    F_HYP2F1,
    F_MIDASMULT,
    F_COMMUTE,
    F_TOEPSOLV,
    F_RGBMIX,
    HF_FELOGITR,
    FN_MAX,	  /* SEPARATOR: end of n-arg functions */
};

enum {
    CONST_PI = 1,
    CONST_NA,
    CONST_INF,
    CONST_NAN,
    CONST_WIN32,
    CONST_EPS,
    CONST_HAVE_MPI,
    CONST_MPI_RANK,
    CONST_MPI_SIZE,
    CONST_N_PROC,
    CONST_TRUE,
    CONST_FALSE,
    CONST_SYSINFO
};

enum {
    DUM_NULL = 1,
    DUM_EMPTY,
    DUM_DIAG,
    DUM_UPPER,
    DUM_LOWER,
    DUM_REAL,
    DUM_IMAG,
    DUM_END,
    DUM_DATASET,
    DUM_TREND
};

#define GENSTRLEN 128
#define NO_VNUM -1

#define unary_op(s)  (s >= 1 && s < U_MAX)
#define binary_op(s) (s > U_MAX && s < OP_MAX)
#define bool_comp(s) (s >= B_EQ && s <= B_OR)
#define dot_op(s)    (s >= B_DOTMULT && s <= B_DOTNEQ)

#define func1_symb(s) (s > F1_MIN && s < F1_MAX)
#define func2_symb(s) (s > F1_MAX && s < F2_MAX)
#define func3_symb(s) (s > F2_MAX && s < F3_MAX)
#define funcn_symb(s) (s > F3_MAX && s < FN_MAX)

#define bnsym(s) (s == FARGS)

#define alias_reversed(n) (n->flags & ALS_NODE)

/* function with single string argument */
#define string_arg_func(s) (s == F_ISDISCR || s == F_OBSNUM || \
			    s == F_BACKTICK || s == F_VARNUM || \
			    s == F_REMOVE || s == F_EXISTS)

/* function with multiple args, string for first arg */
#define str0_func(s) (s == F_PVAL || s == F_CDF || s == F_INVCDF || \
		      s == F_CRIT || s == F_RANDGEN || s == F_PDF || \
		      s == F_BESSEL || s == F_MRANDGEN || s == F_RANDGEN1)

/* functions taking a string arg in last position */
#define string_last_func(s) (s == F_AGGRBY || s == F_GENSERIES || \
			     s == F_PRINTF || s == F_SPRINTF || \
			     s == F_ALLREDUCE || s == F_NORMTEST || \
			     s == F_SSCANF || s == F_NPCORR || \
			     s == F_INBUNDLE || s == F_ASORT)

/* functions taking string arg in middle position */
#define string_mid_func(s) (s == F_REDUCE || s == F_SCATTER)

/* functions taking one or more "fncall" (string) arguments */
#define fncall_func(s) (s == F_BFGSMAX || s == F_NRMAX || \
			s == F_FDJAC || s == F_SIMANN || \
			s == F_BFGSCMAX || s == F_NMMAX || \
			s == F_GSSMAX || s == F_NUMHESS || \
			s == F_FZERO)

/* functions with "reversing" aliases */
#define als_func(s) (s == F_BFGSMAX || s == F_NRMAX || \
		     s == F_SIMANN || s == F_BFGSCMAX || \
		     s == F_NMMAX || s == F_GSSMAX || \
		     s == F_EXISTS)

/* functions where the right-hand argument is actually a return
   location */
#define r_return(s) (s == F_QR || s == F_EIGSYM || s == F_EIGEN || \
                     s == F_MOLS || s == F_MPOLS || s == F_SVD || \
		     s == F_EIGGEN)

/* functions where the middle argument is actually a return
   location */
#define m_return(s) (s == F_SVD || s == F_EIGEN || s == F_QR)

#define undef_arg_ok(s) (s == F_TYPEOF || s == F_TYPENAME || \
			 s == F_ISCMPLX)

#define reusable(p) (p->flags & (P_COMPILE | P_EXEC))

typedef struct node NODE;

struct branchn {
    int n_nodes;
    NODE **n;
};

union val {
    struct branchn bn;
    int idnum;
    char *str;
    double xval;
    double *xvec;
    int *ivec;
    gretl_matrix *m;
    matrix_subspec *mspec;
    gretl_bundle *b;
    gretl_array *a;
    void *ptr;
};

enum node_flags {
    AUX_NODE = 1 << 0, /* auxiliary: free on exit */
    TMP_NODE = 1 << 1, /* temporary: free content on exit */
    PRX_NODE = 1 << 2, /* aux node is proxy (don't reuse!) */
    LHT_NODE = 1 << 3, /* node holds terminal of LHS */
    MUT_NODE = 1 << 4, /* node is inherently mutable in type */
    ALS_NODE = 1 << 5  /* function subject to "reversing" alias */
};

struct node {
    gint16 t;        /* type identifier */
    guint8 flags;    /* AUX_NODE etc., see above */
    int vnum;        /* associated series ID number */
    char *vname;     /* associated variable name */
    user_var *uv;    /* associated named variable */
    union val v;     /* value (of whatever type) */
    NODE *L, *M, *R; /* up to three child nodes */
    NODE *aux;       /* auxiliary (result) node */
    NODE *parent;    /* parent node (or NULL) */
    int refcount;    /* reference counter, used by aux nodes */
};

typedef enum {
    P_DISCARD = 1 <<  0, /* compute and print, don't save */
    P_START   = 1 <<  1, /* first round of evaluation */
    P_AUTOREG = 1 <<  2, /* expression is autoregressive */
    P_DECL    = 1 <<  3, /* statement is actually a declaration */
    P_PRIV    = 1 <<  4, /* generating a "private" or internal var */
    P_COMPILE = 1 <<  5, /* compiling the parse tree */
    P_EXEC    = 1 <<  6, /* evaluating a compiled tree */
    P_NATEST  = 1 <<  7, /* testing for NAs in expression */
    P_UFRET   = 1 <<  8, /* returning value generated by user function */
    P_QUIET   = 1 <<  9, /* don't print any messages or labels */
    P_GETSTR  = 1 << 10, /* state: flag acceptance of plain strings */
    P_MMASK   = 1 << 11, /* genr result is masked matrix */
    P_SLICING = 1 << 12, /* state: calculating object slice (temporary) */
    P_LAGPRSE = 1 << 13, /* state: parsing lag spec (temporary) */
    P_DELTAN  = 1 << 14, /* flag for change in series length */
    P_CATCH   = 1 << 15, /* "catch" is in force */
    P_NODECL  = 1 << 16, /* type of result was not specified */
    P_LISTDEF = 1 << 17, /* expression defines a list */
    P_ANON    = 1 << 18, /* generating an anonymous object */
    P_VOID    = 1 << 19, /* function call, no assignment */
    P_NOEXEC  = 1 << 20, /* just compile, don't evaluate */
    P_MSAVE   = 1 << 21, /* trying for reuse of an aux matrix */
    P_OBSVAL  = 1 << 22, /* generating value of observation in series */
    P_ALIASED = 1 << 23, /* state: handling aliased object (temporary) */
    P_AND     = 1 << 24, /* state: working on right-hand term of B_AND */
    P_OR      = 1 << 25, /* state: working on right-hand term of B_OR */
    P_STACK   = 1 << 26, /* executing stack() */
    P_ALTINP  = 1 << 27, /* the input string has been substituted */
    P_OBJQRY  = 1 << 28, /* querying the existence of an object */
    P_PRNLIST = 1 << 29  /* defining a list for "print" */
} genflags;

struct lhinfo {
    int t;                 /* type of pre-existing LHS variable, if any */
    char name[VNAMELEN];   /* name of LHS variable */
    char *label;           /* descriptive string for series */
    series_table *stab;    /* holds string values for series */
    int vnum;              /* ID number of pre-existing LHS series */
    user_var *uv;          /* address of pre-existing LHS variable */
    char *expr;            /* expression on left */
    GretlType gtype;       /* gretl type of LHS array, if any, or
			      of LHS bundle member */
    gretl_matrix *mret;    /* matrix output (possibly under bundle or array) */
};

typedef struct parser_ parser;

struct parser_ {
    const char *input; /* complete input string */
    const char *point; /* remaining unprocessed input */
    const char *rhs;   /* for use in labelling */
    DATASET *dset;     /* convenience pointer to dataset */
    PRN *prn;          /* for printing messages */
    PRN *errprn;       /* for storing error message in case @prn is NULL */
    genflags flags;    /* various attributes (see @genflags above) */
    int targ;          /* target type */
    int op;            /* assignment operator (possibly inflected) */
    struct lhinfo lh;  /* left-hand side info */
    NODE *lhtree;      /* LHS syntax tree, if needed */
    NODE *lhres;       /* result of eval() on @lhtree */
    NODE *tree;        /* RHS syntax tree */
    NODE *ret;         /* result of eval() on @tree */
    /* below: parser state variables */
    NODE *aux;         /* convenience pointer to current auxiliary node */
    int callcount;
    int dset_n;
    int obs;
    int sym;
    int upsym;
    int ch;
    double xval;
    int idnum;
    char *idstr;
    void *data;
    int err;
};

int parser_getc (parser *p);
void parser_ungetc (parser *p);
void parser_advance (parser *p, int n);
int parser_char_index (parser *p, int c);
int parser_print_input (parser *p);
void free_tree (NODE *t, parser *p, int code);
void lex (parser *s);
NODE *new_node (int t);
NODE *expr (parser *s);
NODE *newdbl (double x);
NODE *newempty (void);
NODE *newb2 (int t, NODE *l, NODE *r);
NODE *obs_node (parser *p);
NODE *bncopy (NODE *t, int *err);
const char *getsymb (int t);
const char *getsymb_full (int t, const parser *p);
void set_parsing_query (int s);
void set_doing_genseries (int s);

int parser_ensure_error_buffer (parser *p);
void context_error (int c, parser *p, const char *func);
void undefined_symbol_error (const char *s, parser *p);

int realgen (const char *s, parser *p, DATASET *dset,
	     PRN *prn, int flags, int targtype);
void gen_save_or_print (parser *p, PRN *prn);
void gen_cleanup (parser *p);

/* handling declarations of variables, wanted in genfuncs.c */
int check_declarations (char ***pS, parser *p);

/* in genfuncs.c, used only internally */
int cross_sectional_stat (double *x, const int *list,
			  const DATASET *dset,
			  int f, int partial_ok);
int x_sectional_weighted_stat (double *x, const int *list,
			       const int *wlist,
			       const DATASET *dset,
			       int f, int partial_ok);

/* in geneval.c, wanted in geneval.c */
double dvar_get_scalar (int i, const DATASET *dset);
int *node_get_list (NODE *n, parser *p);

/* in genlex.c, used only in geneval.c */
void *get_genr_function_pointer (int f);

#endif /* GENPARSE_H */
