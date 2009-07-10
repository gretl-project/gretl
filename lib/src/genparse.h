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

#include "libgretl.h"
#include "usermat.h"
#include "gretl_func.h"
#include "gretl_scalar.h"

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
              B_DOTGT,
  /* 30 */    B_DOTLT,
              B_KRON,     /* Kronecker product */
              B_HCAT,     /* horizontal concatenation */
              B_VCAT,     /* vertical concatenation */
	      B_LCAT,     /* list concatentation */
	      OP_MAX,     /* SEPARATOR: end of binary operators */
              G_LPR,      /* grouping: left paren */
              G_RPR,      /* right paren */
              G_LBR,      /* left bracket */
              G_RBR,      /* right bracket */
  /* 40 */    G_LCB,      /* left curly bracket */
	      G_RCB,      /* right curly bracket */  
	      P_COM,	  /* punctuation: comma */
	      P_DOT,	  /* period */
	      P_SEMI,	  /* semi-colon */
	      P_COL,	  /* colon */
	      PUNCT_MAX,  /* SEPARATOR: end of grouping and punctuation marks */
              CON,	  /* named constant */
              DUM,	  /* dummy variable */
	      USCALAR,	  /* user variable, scalar */
  /* 50 */    USERIES,    /* user variable, series */
	      UMAT,	  /* user-defined matrix */
	      UOBJ,	  /* user-defined object (e.g. model) */
	      NUM,	  /* scalar, evaluated */
	      VEC,	  /* series, evaluated */
	      IVEC,	  /* vector of integers, evaluated */
	      LVEC,       /* list, as array of ints */
	      MAT,	  /* matrix, evaluated */
	      OBS,	  /* observation from a series */
	      MSL,	  /* matrix plus subspec */
  /* 60 */    DMSL,	  /* "dollar" matrix plus subspec */
	      DMSTR,	  /* "dollar" matrix plus string subspec */
	      MSL2,	  /* unevaluated matrix subspec */
	      MSPEC,	  /* evaluated matrix subspec */
	      SUBSL,	  /* row or column component of MSPEC */
	      MDEF,	  /* explicit matrix definition {...} */
	      LAG,	  
	      DVAR,	  /* $ dataset variable (scalar or series) */
	      MVAR,	  /* $ model var (scalar, series, or matrix) */
	      OVAR,	  /* object variable: variable "under" an object */
  /* 70 */    LIST,	  /* reference to named list */
	      LISTVAR,    /* variable in list */
	      STR,	  /* string */
	      FARGS,	  /* set of n function arguments */
	      EMPTY,      /* "null" */
	      ABSENT,
	      DTYPE_MAX,  /* SEPARATOR: end of "bare" types */
	      EROOT,	  /* dummy root for (...) expression */
              UFUN,	  /* user-defined function */
	      RFUN,       /* GNU R function */
  /* 80 */    VSTR,       /* string variable */
              INC,   
	      DEC,
	      QUERY,
	      EOT,	  /* end of transmission */
	      UNK 
};

/* functions: don't collide with the enumeration above */

enum {
    F1_MIN = 1 << 8,
    F_ABS,
    F_TOINT,
    F_CEIL,
    F_FLOOR,
    F_ROUND,
    F_SIN,
    F_COS,
    F_TAN,
    F_ASIN,
    F_ACOS,
    F_ATAN,
    F_LOG,
    F_LOG10,
    F_LOG2,
    F_EXP,
    F_SQRT,
    F_DIFF,	  /* first difference */
    F_LDIFF,	  /* log difference */
    F_SDIFF,	  /* seasonal difference */
    F_SORT,	  /* ascending sort */
    F_DSORT,	  /* descending sort */
    F_RANKING,    
    F_ODEV,	  /* orthogonal deviation */
    F_NOBS,
    F_T1,
    F_T2,   
    F_CUM,
    F_MISSING,
    F_DATAOK,
    F_MISSZERO,
    F_ZEROMISS,
    F_MEDIAN,
    F_GINI,
    F_SUM,   
    F_MEAN,
    F_MIN,
    F_MAX,
    F_SD,
    F_VCE,	  /* variance */
    F_SST,
    F_CNORM,
    F_DNORM,
    F_QNORM,
    F_GAMMA,	
    F_LNGAMMA,
    F_HPFILT,
    F_BKFILT,
    F_PNOBS,
    F_PMIN,
    F_PMAX,
    F_PMEAN,
    F_PSD,
    F_IMAT,
    F_SUMR,
    F_SUMC,
    F_MEANR,	 
    F_MEANC,
    F_SDC,
    F_MCOV,
    F_MCORR,
    F_CDEMEAN,
    F_CHOL,
    F_PSDROOT,
    F_INV,
    F_INVPD,
    F_GINV,
    F_DIAG,
    F_TRANSP,
    F_VEC,
    F_VECH,	
    F_UNVECH,
    F_ROWS,
    F_COLS,
    F_DET,
    F_LDET,
    F_TRACE,
    F_NORM1,
    F_INFNORM,
    F_RCOND,
    F_RANK, 
    F_OBSNUM,
    F_ISSERIES,
    F_ISSCALAR,
    F_ISLIST,	 
    F_ISSTRING,
    F_ISNULL,
    F_LISTLEN,
    F_PDF,
    F_PVAL,
    F_CDF,
    F_INVCDF,
    F_CRIT, 
    F_URCPVAL,
    F_RANDGEN,
    F_RPOISSON,
    F_MAKEMASK,
    F_VALUES,
    F_NULLSPC,
    F_MEXP,
    F_MINC,
    F_MAXC,
    F_MINR,
    F_MAXR,
    F_IMINC, 
    F_IMAXC,
    F_IMINR,
    F_IMAXR,
    F_FFT,
    F_FFTI,
    F_UPPER,
    F_LOWER,
    F_MREAD,
    F_POLROOTS,
    F_XPX,
    F_GETENV,
    F_ARGNAME,
    F_OBSLABEL,
    F_READFILE,
    F_BACKTICK,
    F_STRLEN,
    F_VARNAME,
    F_VARNUM,
    F_TOLOWER,
    F1_MAX,	  /* SEPARATOR: end of single-arg functions */
    F_COR,
    F_COV,
    F_DUMIFY,
    F_SORTBY,
    F_RUNIFORM,
    F_RNORMAL,
    F_FRACDIFF,
    F_BOXCOX,
    F_ZEROS,
    F_ONES,
    F_MUNIF,
    F_MNORM,
    F_QFORM,
    F_MLAG,
    F_QR,
    F_EIGSYM,	 
    F_EIGGEN,
    F_FDJAC,
    F_BFGSMAX,
    F_LRVAR,
    F_PRINCOMP,
    F_QUANTILE,
    F_CMULT,	  /* complex multiplication */
    F_CDIV,	  /* complex division */
    F_MXTAB,
    F_MWRITE,
    F_MOVAVG,
    F_MRSEL,
    F_MCSEL,
    F_LLAG,
    F_WMEAN,
    F_WVAR,
    F_WSD,
    F_STRSTR,
    F_STRCMP,
    F_COLNAMES,
    F_LJUNGBOX,
    F_MSORTBY,
    F_LINCOMB,
    F_IMHOF,
    F_XMIN,
    F_XMAX,
    F_RESAMPLE,
    F_FCSTATS,
    F_FRACLAG,
    F_MREVERSE,
    F2_MAX,	  /* SEPARATOR: end of two-arg functions */
    F_MSHAPE,
    F_SVD,
    F_MOLS,
    F_MPOLS,
    F_TRIMR,
    F_TOEPSOLV,
    F_CORRGM,
    F_SEQ,
    F_REPLACE,
    F_STRNCMP,
    F_BESSEL,
    F3_MAX,       /* SEPARATOR: end of three-arg functions */
    F_FILTER,
    F_MCOVG,
    F_KFILTER,
    F_KSMOOTH,
    F_KSIMUL,
    FN_MAX,	  /* SEPARATOR: end of n-arg functions */
};

enum {
    CONST_PI = 1,
    CONST_NA,
    CONST_WIN32
};

enum {
    DUM_NULL = 1,
    DUM_DIAG,
    DUM_DATASET,
    DUM_TREND
};

#define GENSTRLEN 128
#define NO_VNUM -1

#define func1_symb(s) (s > F1_MIN && s < F1_MAX)
#define func2_symb(s) (s > F1_MAX && s < F2_MAX)
#define func3_symb(s) (s > F2_MAX && s < F3_MAX)
#define funcn_symb(s) (s > F3_MAX && s < FN_MAX)

#define string_arg_func(s) (s == F_ISSERIES || s == F_ISNULL || \
			    s == F_ISLIST || s == F_LISTLEN || \
			    s == F_ISSTRING || s == F_OBSNUM || \
			    s == F_BACKTICK || s == F_VARNUM)

#define string0_func(s) (s == F_PVAL || s == F_CDF || s == F_INVCDF || \
                         s == F_CRIT || s == F_RANDGEN || s == F_PDF || \
			 s == F_BESSEL)

#define unary_op(s)  (s >= 1 && s < U_MAX)
#define binary_op(s) (s > U_MAX && s < OP_MAX)
#define bool_comp(s) (s >= B_EQ && s <= B_OR)

#define evalb3(s) (func3_symb(s)) /* plus ternary op? */

#define evalb2(s) (binary_op(s) || func2_symb(s) || s == MSL || \
                   s == MSL2 || s == SUBSL || s == LAG || \
                   s == OBS)

#define b1sym(s) (unary_op(s) || func1_symb(s) || funcn_symb(s) || \
                  s == G_LPR || s == EROOT)

#define evalb1(s) (b1sym(s) && !(string0_func(s)) && s != U_ADDR && \
                   !func2_symb(s) && s != EROOT)

#define b2sym(s) (evalb2(s) || s == DMSTR || s == DMSL || \
                  s == OVAR || s == UFUN || s == RFUN || \
                  s == LISTVAR)

#define b3sym(s) (s == QUERY || func3_symb(s))

#define bnsym(s) (s == MDEF || s == FARGS)

#define freestr(s) (s == STR || s == USCALAR || s == UMAT || s == UOBJ || \
                    s == LIST || s == VSTR)

#define bare_data_type(s) (s > PUNCT_MAX && s < DTYPE_MAX)

#define closing_sym(s) (s == G_RPR || s == G_RBR || s == G_RCB)


/* functions where the right-hand argument is actually a return
   location */
#define r_return(s) (s == F_QR || s == F_EIGSYM || s == F_EIGGEN || \
                     s == F_MOLS || s == F_MPOLS || s == F_SVD)

/* functions where the middle argument is actually a return
   location */
#define m_return(s) (s == F_SVD)

#define dollar_node(n) (n->t == DVAR || n->t == MVAR || \
                        (n->t == MSL && n->v.b2.l->v.str[0] == '$'))

#define reusable(p) (p->flags & (P_COMPILE | P_EXEC))

typedef struct node NODE;

struct branch1 {
    NODE *b;  
};

struct branch2 {
    NODE *l, *r; 
};

struct branch3 {
    NODE *l, *m, *r; 
};

struct branchn {
    int n_nodes;
    NODE **n;
};

union val {
    struct branch1 b1; 
    struct branch2 b2; 
    struct branch3 b3; 
    struct branchn bn;
    int idnum;
    char *str;
    double xval; 
    double *xvec;
    int *ivec;
    gretl_matrix *m;
    matrix_subspec *mspec;
};

enum {
    AUX_NODE    = 1 << 0,
    TMP_NODE    = 1 << 1,
    TRANSP_NODE = 1 << 2
};

struct node {
    short t;       /* type indentifier */
    char flags;    /* AUX_NODE etc., see above */
    int vnum;      /* associated variable ID number */
    union val v;   /* value (of whatever type) */
};

enum {
    P_DISCARD = 1 <<  0, /* compute and print, don't save */
    P_START   = 1 <<  1, /* first round of evaluation */
    P_AUTOREG = 1 <<  2, /* expression is autoregressive */
    P_DECL    = 1 <<  3, /* command is actually a declaration */
    P_PRINT   = 1 <<  4, /* command just prints an existing var */
    P_SCALAR  = 1 <<  5, /* return a scalar result (only) */ 
    P_SERIES  = 1 <<  6, /* return a series result (only) */
    P_MATRIX  = 1 <<  7, /* return a matrix result (only) */
    P_STRING  = 1 <<  8, /* return a string result (only) */
    P_LIST    = 1 <<  9, /* return a list result (only) */
    P_PRIVATE = 1 << 10, /* generating a private or internal var */
    P_COMPILE = 1 << 11, /* just compiling tree, not evaluating */
    P_EXEC    = 1 << 12, /* evaluating pre-built tree */ 
    P_SLICE   = 1 << 13, /* compute matrix slice specification */
    P_UFUN    = 1 << 14, /* user function call with no assignment */
    P_SORT    = 1 << 15, /* doing a simple sort of a variable */
    P_NATEST  = 1 << 16, /* testing for NAs in expression */
    P_UFRET   = 1 << 17, /* returning value generated by user function */
    P_LHSCAL  = 1 << 18, /* there was a pre-existing LHS scalar */
    P_LHLIST  = 1 << 19, /* there was a pre-existing LHS list */
    P_LHSTR   = 1 << 20, /* there was a pre-existing LHS string */
    P_LHMAT   = 1 << 21, /* there was a pre-existing LHS matrix */
    P_QUIET   = 1 << 22, /* don't print any messages or labels */
    P_GETSTR  = 1 << 23, /* state: flag acceptance of plain strings */
    P_LOOP    = 1 << 24, /* executing pre-built genr in loop context */
    P_SLAVE   = 1 << 25  /* running as "slave" of NLS/MLE/GMM */
};

struct lhinfo {
    int t;                 /* type of result */
    char name[VNAMELEN];   /* name of LHS variable */   
    char label[MAXLABEL];  /* descriptive string for var */
    int v;                 /* ID number for variable */
    int obs;               /* specific obs number in series */
    gretl_matrix *m0;      /* original LHS matrix (or NULL) */
    gretl_matrix *m1;      /* computed LHS matrix */
    char *substr;          /* obs or matrix selection string */
    matrix_subspec *mspec; /* evaluated submatrix spec */
};

typedef struct parser_ parser;

struct parser_ {
    const char *input; /* complete input string */
    const char *point; /* remaining unprocessed input */
    const char *rhs;   /* for use in labelling */
    double ***Z;       /* convenience pointer to data array */
    DATAINFO *dinfo;   /* convenience pointer to data info */
    PRN *prn;          /* for printing messages */
    int flags;         /* various attributes (see above) */
    int targ;          /* target type */
    int op;            /* assignment operator (possibly inflected) */
    struct lhinfo lh;  /* left-hand side info */
    parser *subp;      /* left-hand side matrix subslice tree */
    NODE *tree;        /* parsed syntax tree */
    NODE *ret;         /* evaluated result node */
    NODE **aux;        /* auxiliary nodes used in evaluation */
    int n_aux;         /* the number of the above */
    int aux_i;         /* the current ID of the above */
    char warning[64];  /* to hold a warning, if needed */
    int loopline;      /* line number within loop, if applicable */
    /* below: parser state variables */
    int obs;
    int sym;
    int ch;
    double xval;
    int idnum;
    char *idstr;
    int err;
    int warn;
};

#define starting(p) (p->flags & P_START)
#define autoreg(p) (p->flags & P_AUTOREG)
#define simple_ufun_call(p) (p->flags & P_UFUN)

int parser_getc (parser *p);
void parser_ungetc (parser *p);
void parser_advance (parser *p, int n);
int parser_charpos (parser *p, int c);
int parser_next_char (parser *p);
void parser_print_input (parser *p);
void lex (parser *s);
NODE *expr (parser *s);
NODE *newdbl (double x);
NODE *obs_node (parser *p);
NODE *msl_node_direct (parser *p);
void context_error (int c, parser *p);
const char *getsymb (int t, const parser *p);

int realgen (const char *s, parser *p, double ***Z, 
	     DATAINFO *pdinfo, PRN *prn, int flags);
void gen_save_or_print (parser *p, PRN *prn);
void gen_cleanup (parser *p);
void parser_free_aux_nodes (parser *p);

void set_matrix_slice_on (void);
void set_matrix_slice_off (void);

/* name lookup functions */
const char *constname (int c);
const char *dvarname (int t);
const char *mvarname (int t);
const char *dumname (int t);

/* handling declarations of variables */
int check_declarations (char ***pS, parser *p);

/* handling sorted marker strings */
int maybe_pick_up_sorted_markers (parser *p);

/* in genfuncs.c, used only internally */
int cross_sectional_stat (double *x, const int *list, 
			  const double **Z, 
			  const DATAINFO *pdinfo,
			  int f);
int x_sectional_weighted_stat (double *x, const int *list, 
			       const int *wlist,
			       const double **Z, 
			       const DATAINFO *pdinfo,
			       int f);

/* helper functions for manual, gretl.lang file */
int gen_func_count (void);
const char *gen_func_name (int i);
int model_var_count (void);
const char *model_var_name (int i);
int data_var_count (void);
const char *data_var_name (int i);
