/*
 *   Copyright (c) by Allin Cottrell
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

/* shared private header for all 'genr' related modules */

#include "libgretl.h"
#include "usermat.h"
#include "gretl_func.h"

#define GENDEBUG 0

/* this may be a good idea, but perhaps risky? */
#define PRESERVE_AUX_NODES 0

enum {
    U_NEG = 1,
    U_POS,
    U_NOT,
    U_ADDR,
    U_MAX,   /* separator: end of unary operators */ 
    B_ASN,
    B_ADD,
    B_SUB,
    B_MUL,
    B_DIV,  /* 10 */
    B_MOD, 
    B_POW,
    B_EQ,
    B_LT,
    B_GT,
    B_LTE,
    B_GTE,
    B_NEQ,
    B_AND,
    B_OR,   /* 20 */
    B_TRMUL,
    LPR,     /* left paren */
    RPR,     /* right paren */
    LBR,     /* left bracket */
    RBR,     /* right bracket */
    LCB,     /* left curly bracket */
    RCB,     /* right curly bracket */
    DOTMULT,
    DOTDIV,
    DOTPOW, /* 30 */
    DOTADD,
    DOTSUB,
    DOTEQ,
    KRON,     /* Kronecker product */
    MCCAT,    /* matrix concatenation (columns) */
    MRCAT,    /* matrix concatenation (rows) */
    OP_MAX,   /* separator: end of operators */
    ABS,
    TOINT,
    SIN,    /* 40 */
    COS,
    TAN,
    ATAN,
    LOG,
    LOG10,
    LOG2,
    EXP,
    SQRT,
    DIF,      /* first difference */
    LDIF,   /* 50: log difference */
    SDIF,     /* seasonal difference */
    SORT,     /* ascending sort */
    DSORT,    /* descending sort */
    ODEV,     /* orthogonal deviation */
    NOBS,
    T1,
    T2,
    CHISQ,
    STUDENT,
    CUM,      /* 60 */  
    MISSING,
    OK,
    MISSZERO,
    ZEROMISS,
    MEDIAN,
    GINI,
    SUM,
    MEAN,
    MIN,
    MAX,    /* 70 */
    SD,
    VCE,      /* variance */
    SST,
    CNORM,
    DNORM,
    QNORM,
    GAMMA,
    LNGAMMA,
    HPFILT,  
    BKFILT,   /* 80 */
    RESAMPLE,
    PMEAN,
    PSD,
    IMAT,
    SUMR,
    SUMC,
    MEANR,
    MEANC,
    MCOV,    
    MCORR,     /* 90 */
    CDEMEAN,
    CHOL,
    INV,
    DIAG,
    TRANSP,
    TVEC,
    VECH,
    UNVECH,
    ROWS,   
    COLS,      /* 100 */
    DET,
    LDET,
    TRACE,
    NORM1,
    INFNORM,
    RCOND,
    VARNUM,
    OBSNUM,
    ISSERIES,
    ISLIST,    /* 110 */  
    ISSTRING,
    ISNULL,
    LISTLEN,
    PVAL,
    CDF,
    CRIT,
    GENPOIS,
    MAKEMASK,
    VALUES,
    NULLSPC,  /* 120 */
    MEXP,
    FUNC_MAX,   /*separator: end of single-arg functions */
    COR,
    COV,
    UNIFORM,
    NORMAL,   
    BINOMIAL,
    FRACDIF,
    ZEROS,
    ONES,    /* 130 */  
    MUNIF,
    MNORM,
    QFORM,
    MLAG,
    QR,
    EIGSYM,  
    EIGGEN,
    FDJAC,
    BFGSMAX,
    LRVAR,  /* 140: long-run variance */
    F2_MAX,   /* separator: end of two-arg functions */
    MSHAPE,
    SVD,
    FN_MAX,   /* separator: end of n-arg functions */
    COM,      /* comma */
    DOT,      /* period */
    SEMI,     /* semi-colon */
    COL,      /* colon */
    CON,      /* named constant */
    DUM,    /* 150: dummy variable */
    UVAR,     /* user variable (scalar or series) */
    UMAT,     /* user-defined matrix */
    UOBJ,     /* user-defined object (e.g. model) */
    NUM,      /* scalar, evaluated */
    VEC,      /* series, evaluated */
    IVEC,     /* vector of integers, evaluated */
    MAT,      /* matrix, evaluated */
    OBS,      /* observation from a series */
    MSL,      /* matrix plus subspec */
    DMSL,   /* 160: "dollar" matrix plus subspec */
    DMSTR,    /* "dollar" matrix plus old-style string subspec */
    MSL2,     /* unevaluated matrix subspec */
    MSPEC,    /* evaluated matrix subspec */
    SUBSL,    /* row or column component of MSPEC */
    MDEF,     /* explicit matrix definition {...} */
    LAG,
    DVAR,     /* $ dataset variable (scalar or series) */
    MVAR,     /* $ model var (scalar, series, or matrix) */
    OVAR,     /* object variable: variable "under" an object */
    LOOPIDX, /* 170: loop index variable */
    LIST,     /* reference to named list */
    STR,      /* string */
    EROOT,    /* dummy root for (...) expression */
    UFUN,     /* user-defined function */
    FARGS,    /* set of n function arguments */
    EMPTY,
    ABSENT, 
    INC,    
    DEC,
    QUERY,
    AST2,     /* prime facie ambiguous case of "**" */
    UNK
};

enum {
    CONST_PI = 1,
    CONST_NA
};

enum {
    DUM_NULL = 1,
    DUM_DIAG,
    DUM_DATASET
};

#define MAXSTR 128

#define func_symb(s) ((s > OP_MAX && s < FUNC_MAX) || \
                       s == LAG || s == OBS)
#define func2_symb(s) (s > FUNC_MAX && s < F2_MAX)
#define funcn_symb(s) (s > F2_MAX && s < FN_MAX)
#define string_arg_func(s) (s == VARNUM || s == ISSERIES || s == ISNULL || \
                            s == ISLIST || s == LISTLEN || s == ISSTRING || \
                            s == OBSNUM)

#define string0_func(s) (s == PVAL || s == CDF || s == CRIT)

#define unary_op(s) (s >= 1 && s < U_MAX)
#define binary_op(s) (s > U_MAX && s < OP_MAX)
#define bool_comp(s) (s >= B_EQ && s <= B_OR)

#define evalb2(s) (binary_op(s) || func2_symb(s) || s == MSL || \
                   s == MSL2 || s == SUBSL)

#define b1sym(s) (unary_op(s) || func_symb(s) || funcn_symb(s) || \
                  s == LPR || s == EROOT)

#define evalb1(s) (b1sym(s) && !(string0_func(s)) && s != U_ADDR && \
                   !func2_symb(s))

#define b2sym(s) (evalb2(s) || s == DMSTR || s == OVAR || s == UFUN)

#define b3sym(s) (s == QUERY)

#define bnsym(s) (s == MDEF || s == FARGS)

#define freestr(s) (s == STR || s == UMAT || s == UOBJ || \
                    s == LOOPIDX || s == LIST)

/* functions where the right-hand "argument" is actually a return
   location */
#define r_return(s) (s == QR || s == EIGSYM || s == EIGGEN)

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

struct node {
    int t;         /* type indentifier */
    int ext;       /* extra information for some cases */
    int tmp;       /* if non-zero, node holds temporary data which
		      should be freed on completion */
    int aux;       /* non-zero for an auxiliary node */
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
    P_PRIVATE = 1 <<  7, /* generating a private or internal var */
    P_COMPILE = 1 <<  8, /* just compiling tree, not evaluating */
    P_EXEC    = 1 <<  9, /* evaluating pre-built tree */ 
    P_SLICE   = 1 << 10, /* compute matrix slice specification */
    P_UFUN    = 1 << 11, /* user function call with no assignment */
    P_SORT    = 1 << 12, /* doing a simple sort of a variable */
    P_NATEST  = 1 << 13, /* testing for NAs in expression */
    P_UFRET   = 1 << 14  /* returning value generated by user function */
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
    NODE *tree;        /* parsed syntax tree */
    NODE *ret;         /* evaluated result node */
    NODE **aux;        /* auxiliary nodes used in evaluation */
    int n_aux;         /* the number of the above */
    int aux_i;         /* the current ID of the above */
    int ecount;        /* number of times this parser has been exec'd */
    char warning[64];  /* to hold a warning, if needed */
    /* below: parser state variables */
    int obs;
    int sym;
    int ch;
    double xval;
    int idnum;
    char *idstr;
    int getstr;
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
int function_lookup (const char *s);

int realgen (const char *s, parser *p, double ***Z, 
	     DATAINFO *pdinfo, PRN *prn, int flags);
void gen_save_or_print (parser *p, PRN *prn);
void gen_cleanup (parser *p);
void parser_free_aux_nodes (parser *p);

/* name lookup functions */
const char *constname (int c);
const char *dvarname (int t);
const char *mvarname (int t);
const char *dumname (int t);

/* for handling user-functions */
void fn_args_init (fnargs *args);
void fn_args_free (fnargs *args);
int push_fn_arg (fnargs *args, int type, void *p);

/* handling declarations of variables */
int check_declarations (char ***pS, parser *p);

/* handling sorted market strings */
int maybe_pick_up_sorted_markers (parser *p);

/* helper functions for manual, gretl.lang file */
int gen_func_count (void);
const char *gen_func_name (int i);
int model_var_count (void);
const char *model_var_name (int i);
int data_var_count (void);
const char *data_var_name (int i);
