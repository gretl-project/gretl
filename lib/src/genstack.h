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

#ifndef GENSTACK_H
#define GENSTACK_H

/* #define GENR_DEBUG */

#ifdef GENR_DEBUG
void dprintf (const char *format, ...);
# define DPRINTF(x) dprintf x
#else 
# define DPRINTF(x)
#endif /* GENR_DEBUG */

#define ATOMLEN 32  /* length of auxiliary string in genr atom */

typedef struct _genatom genatom;

struct _genatom {
    char level;
    char scalar;
    int varnum;
    int tmpvar;
    char lag;
    double val;
    char func;
    char op;
    char popped;
    char str[ATOMLEN];
    genatom *parent;
};

/* below: matches funcs[] in generate.c */
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
    T_SST,
    T_COV,
    T_MEDIAN,
    T_ZEROMISS,
    T_PVALUE,
    T_OBSNUM,
    T_MPOW,
    T_DNORM,
    T_CNORM,
    T_RESAMPLE,
#ifdef HAVE_MPFR
    T_MLOG,
#endif
    T_IDENTITY
};

int push_atom (genatom *atom);
genatom *pop_atom (void);
genatom *pop_child_atom (genatom *atom);
genatom *peek_child_atom (genatom *atom);
void reset_atom_stack (void);
void destroy_atom_stack (void);
void atom_stack_set_parentage (void);
void atom_eat_children (genatom *atom);
void atom_stack_bookmark (void);
void atom_stack_resume (void);
int atom_stack_check_for_scalar (void);

int calc_push (double x);
double calc_pop (void);
void reset_calc_stack (void);

#endif /* GENSTACK_H */
