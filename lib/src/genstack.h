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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GENSTACK_H
#define GENSTACK_H

#define GENR_DEBUG 0
#define GEN_MATRIX_DEBUG 0

#if (GENR_DEBUG || GEN_MATRIX_DEBUG)
void dprintf (const char *format, ...);
#endif

#if GENR_DEBUG
# define DPRINTF(x) dprintf x
#else 
# define DPRINTF(x)
#endif

#if GEN_MATRIX_DEBUG
# define MPRINTF(x) dprintf x
#else 
# define MPRINTF(x)
#endif

#define ATOMLEN 32         /* length of auxiliary string in genr atom */
#define MATRIX_ON_ATOM 333 /* int to indicate that matrix is temporarily
			      attached to a genatom */

#define matrix_is_on_atom(m) (MATRIX_ON_ATOM == gretl_matrix_get_int(m))
#define set_matrix_on_atom(m) (gretl_matrix_set_int(m, MATRIX_ON_ATOM))
#define unset_matrix_on_atom(m) (gretl_matrix_set_int(m, 0))

typedef enum {
    ATOM_SERIES   = 0,
    ATOM_SCALAR   = 1 << 0,
    ATOM_MATRIX   = 1 << 1,
    ATOM_MODELDAT = 1 << 2,
    ATOM_TMP      = 1 << 3
} GenAtomType;

typedef struct genatom_ genatom;
typedef struct atomset_ atomset;

struct genatom_ {
    char level;        /* evaluation priority level */
    char atype;        /* flags giving type of atom (see GenAtomType) */
    int varnum;        /* variable ID number */
    int varobs;        /* observation number, if applicable */
    int tmpvar;        /* ID number of temporary variable, if any */
    char lag;          /* lag, if lagged variable */
    double val;        /* numerical value, if applicable */
    char func;         /* function code, if atom is function */
    char op;           /* operator */
    char popped;       /* flag indicating whether or not popped */
    char str[ATOMLEN]; /* string info saved on atom */
    gretl_matrix *M;   /* matrix pointer, if matrix atom */
    genatom *parent;   /* pointer to parent, if child */
    atomset *aset;     /* pointer to set to which the atom belongs */
};

/* below: matches funcs[] in generate.c -- these are in addition to
   the (mostly) standard math functions in the GretlMathFunc
   enumeration
*/
enum transformations {
    T_DIFF = T_MATHMAX,
    T_LDIFF, 
    T_SDIFF,
    T_MEAN, 
    T_SD, 
    T_MIN,
    T_MAX,
    T_SORT, 
    T_DSORT,
    T_SUM, 
    T_NOBS,
    T_T1,
    T_T2,
    T_CUM, 
    T_MISSING,
    T_OK,
    T_MISSZERO,
    T_CORR,
    T_VAR,
    T_SST,
    T_COV,
    T_MEDIAN,
    T_GINI,
    T_ZEROMISS,
    T_PVALUE,
    T_CRIT,
    T_OBSNUM,
    T_MPOW,
#ifdef HAVE_MPFR
    T_MLOG,
#endif
    T_RESAMPLE,
    T_HPFILT,
    T_BKFILT,
    T_FRACDIFF,
    T_VARNUM,
    T_SERIES,
    T_ISLIST,
    T_NELEM,
    T_DET,
    T_INV,
    T_CHOL,
    T_CDMEAN,
    T_DIAG,
    T_QR,
    T_EIGSYM,
    T_EIGGEN,
    T_LDET,
    T_TRACE,
    T_1NORM,
    T_RCOND,
    T_ROWS,
    T_COLS,
    T_TRANSP,
    T_IMAT,
    T_ZEROS,
    T_ONES,
    T_IDENTITY
};

#define VALSTACK_SIZE 32
#define MATSTACK_SIZE 32

enum genr_flags {
    GENR_SAVE         = 1 << 0,
    GENR_SCALAR       = 1 << 1,
    GENR_FORCE_SERIES = 1 << 2,
    GENR_NEED_SCALAR  = 1 << 3,
    GENR_WARN         = 1 << 4,
    GENR_SIMPLE_SORT  = 1 << 5,
    GENR_PRIVATE      = 1 << 6,
    GENR_MATRIX       = 1 << 7,
    GENR_DEPOSIT      = 1 << 8
};

struct _GENERATOR {
    int err;                        /* error code */  
    int done;                       /* completion indicator */
    int flags;                      /* option flags */
    char orig_s[MAXLINE];           /* formula provided by user */
    char lhs[32];                   /* left-hand side of formula */
    double *xvec;                   /* temporary storage */
    int varnum;                     /* number of variable generated */
    int obs;                        /* observation number, if generating
				       a single observation */
    char varname[32];               /* name of variable generated */
    char label[MAXLABEL];           /* descriptive label for variable,
				       or submatrix specification */
    int tmpv;                       /* temporary variable ID */
    double **tmpZ;                  /* temporary dataset */
    DATAINFO *pdinfo;               /* pointer to "outer" data info */
    double ***pZ;                   /* pointer to "outer" dataset */
    atomset *aset;                  /* set of atomic compoenents */
    double valstack[VALSTACK_SIZE]; /* stack of temporary values */
    int nvals;                      /* number of temporary values */
    gretl_matrix **mstack;          /* stack of temporary matrices */
    int nmats;                      /* number of same */
    char **S;                       /* array of obs strings, used in
				       certain cases of sorting */
    PRN *prn;                       /* for printing messages */
};

int attach_atomset (GENERATOR *genr);
int push_atom (genatom *atom);
genatom *pop_atom (GENERATOR *genr);
genatom *pop_child_atom (genatom *atom);
genatom *peek_child_atom (genatom *atom);
void reset_atom_stack (GENERATOR *genr);
void destroy_atom_stack (GENERATOR *genr);
void atom_stack_set_parentage (GENERATOR *genr);
void atom_eat_children (genatom *atom);
void atom_stack_bookmark (GENERATOR *genr);
void atom_stack_resume (GENERATOR *genr);
int atom_stack_check_for_scalar (GENERATOR *genr);

gretl_matrix *atom_stack_get_matrix (GENERATOR *genr, const char *str);
void atom_stack_nullify_matrix (const gretl_matrix *M, GENERATOR *genr);

int arg_atom_available (genatom *atom);

int calc_push (double x, GENERATOR *genr);
double calc_pop (GENERATOR *genr);
void reset_calc_stack (GENERATOR *genr);

int matrix_calc_push (gretl_matrix *M, GENERATOR *genr);
gretl_matrix *matrix_calc_pop (GENERATOR *genr);
void reset_matrix_calc_stack (GENERATOR *genr);

const char *get_genr_func_word (int fnum);

#endif /* GENSTACK_H */
