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

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
# include "config.h" /* for HAVE_MPFR */
#endif

#ifdef _WIN32
# include "winconfig.h"
#endif

#include "libgretl.h"
#include "genstack.h"
#include "usermat.h"

enum {
    STACK_PUSH,
    STACK_POP,
    STACK_GET,
    STACK_RESET,
    STACK_DEPEND,
    STACK_POP_CHILDREN,
    STACK_EAT_CHILDREN,
    STACK_PEEK_CHILDREN,
    STACK_BOOKMARK,
    STACK_RESUME,
    STACK_SCALAR_CHECK,
    STACK_GET_MATRIX,
    STACK_DESTROY
};

struct atomset_ {
    genatom **atoms;
    int n_atoms;
    int n_popped;
    int bookmark;
};

#define must_know_children(f) (f == T_NOBS || f == T_MEAN || f == T_SUM || \
                               f == T_SD || f == T_VAR || f == T_SST || \
                               f == T_MEDIAN || f == T_MIN || f == T_MAX || \
                               f == T_DIFF || f == T_LDIFF || f == T_SDIFF || \
                               f == T_T1 || f == T_T2 || f == T_GINI || \
                               f == T_CUM || f == T_SORT || f == T_DET || \
                               f == T_INV || f == T_CDMEAN || f == T_RCOND || \
                               f == T_CHOL || f == T_QR || f == T_1NORM || \
                               f == T_LDET || f == T_TRACE || f == T_DIAG || \
                               f == T_ROWS || f == T_COLS || f == T_TRANSP || \
                               f == T_VARNUM || f == T_SERIES || \
                               f == T_ISLIST || f == T_NELEM || \
                               f == T_RESAMPLE || f == T_HPFILT || \
                               f == T_BKFILT || f == T_DSORT)

#define atom_is_scalar(a) ((a->atype & ATOM_SCALAR) || (a->atype & ATOM_TMP))

static int all_children_scalar (genatom **atoms, int n,
				int pos, int level)
{
    int j, sc = 1;

    if (pos == 0) return 0;

    for (j=pos-1; j>=0; j--) {
	/* below: unsure between "!= level + 1" and "<= level" */
	if (atoms[j]->level <= level) { 
	    if (j == pos - 1) sc = 0;
	    break;
	}
	if (!atom_is_scalar(atoms[j])) {
	    sc = 0; 
	    break;
	}
    }

    return sc;
}

static int real_check_for_scalar_result (genatom **atoms, int n)
{
    int i, k;
    int maxlevel = 0;
    int scalar = 1;

    /* find highest nesting depth */
    for (i=0; i<n; i++) {
	if (atoms[i]->level > maxlevel) {
	    maxlevel = atoms[i]->level;
	}
    }

    DPRINTF(("check_for_scalar: maxlevel = %d\n", maxlevel));

    /* recurse to identify scalar property of composites */
    for (k=maxlevel-1; k>=0; k--) {
	DPRINTF(("check_for_scalar: checking level %d\n", k));
	for (i=0; i<n; i++) {
	    if (atoms[i]->level == k && !atom_is_scalar(atoms[i])) {
		DPRINTF(("check_for_scalar: checking atom %d\n", i));
		if (all_children_scalar(atoms, n, i, k)) {
		    DPRINTF(("check_for_scalar: setting ATOM_TMP "
			     "for atom[%d]\n", i));
		    atoms[i]->atype |= ATOM_TMP;
		}
	    }
	}
    }

    for (i=0; i<n; i++) {
	if (atoms[i]->level == 0 && !atom_is_scalar(atoms[i])) {
	    scalar = 0;
	}
	if (atoms[i]->atype & ATOM_TMP) {
	    /* restore original atom type */
	    atoms[i]->atype &= ~ATOM_TMP;
	}
    }

    return scalar;
}

static void free_atom (genatom *atom, int i)
{
    if (atom->M != NULL) {
	if (!is_user_matrix(atom->M)) {
	    MPRINTF(("atom %d: freeing matrix at %p\n", i, (void *) atom->M));
	    gretl_matrix_free(atom->M);
	} else {
	    unset_matrix_on_atom(atom->M);
	}
    }

    free(atom);
}

static int real_stack_eat_children (genatom *parent, 
				    genatom **atoms, int n)
{
    int i, j, ndel = 0;

    for (i=0; i<n; i++) {
	if (atoms[i]->parent == parent) {
	    DPRINTF(("freeing child atom, pos %d\n", i));
	    free_atom(atoms[i], i);
	    for (j=i; j<n-1; j++) {
		atoms[j] = atoms[j+1];
	    }
	    atoms[n-1] = NULL;
	    ndel++;
	    i--;
	    n--;
	}
    }

    DPRINTF(("ate %d children\n", ndel));

    return ndel;
}

static void real_stack_set_parentage (genatom **atoms, int n)
{
    int i, j, level;

    for (i=n-1; i>=0; i--) {
	DPRINTF(("checking for children: looking at atom %d\n", i));
	if (must_know_children(atoms[i]->func)) {
	    level = atoms[i]->level;
	    DPRINTF((" got candidate (atom %d, level %d)...\n", i, level));
	    for (j=i-1; j>=0; j--) {
		DPRINTF(("  looking at atom %d\n", j));
		if (atoms[j]->level > level) {
		    DPRINTF(("    marking atom %d as parent of atom %d\n",
			    i, j));
		    atoms[j]->parent = atoms[i];
		} else {
		    DPRINTF(("    not a child, breaking\n"));
		    j++;
		    break;
		}
	    }
	    i = j;
	}
    }
}

static genatom *atom_stack (genatom *atom, atomset *aset, int op)
{
    genatom *ret = NULL;
    int i, j;

    if (op == STACK_PUSH && atom != NULL) {
	genatom **atoms = 
	    realloc(aset->atoms, (aset->n_atoms + 1) * sizeof *atoms);

	if (atoms != NULL) {
	    aset->atoms = atoms;
	    aset->atoms[aset->n_atoms] = atom;
	    aset->n_atoms += 1;
	    ret = atom;
	}
    } else if (op == STACK_POP) {
	if (aset->n_popped < aset->n_atoms) {
	    ret = aset->atoms[aset->n_popped];
	    aset->n_popped += 1;
	}
    } else if (op == STACK_RESET) {
	for (j=0; j<aset->n_atoms; j++) {
	    aset->atoms[j]->popped = 0;
	}
	aset->n_popped = 0;
    } else if (op == STACK_DESTROY) {
	for (j=0; j<aset->n_atoms; j++) {
	    if (aset->atoms[j]->M != NULL) {
		for (i=j+1; i<aset->n_atoms; i++) {
		    if (aset->atoms[i]->M == aset->atoms[j]->M) {
			aset->atoms[i]->M = NULL;
		    }
		}
	    }
	    free_atom(aset->atoms[j], j);
	}
	free(aset->atoms);
	aset->atoms = NULL;
	aset->n_atoms = 0;
	aset->n_popped = 0;
    } else if (op == STACK_POP_CHILDREN && atom != NULL) {
	for (j=0; j<aset->n_atoms; j++) {
	    if (aset->atoms[j]->parent == atom && !aset->atoms[j]->popped) {
		ret = aset->atoms[j];
		aset->atoms[j]->popped = 1;
		break;
	    }
	}
    } else if (op == STACK_PEEK_CHILDREN && atom != NULL) {
	for (j=0; j<aset->n_atoms; j++) {
	    if (aset->atoms[j]->parent == atom) {
		ret = aset->atoms[j];
		break;
	    }
	}	
    } else if (op == STACK_EAT_CHILDREN && atom != NULL) {
	int ndel = real_stack_eat_children(atom, aset->atoms, aset->n_atoms);

	aset->n_atoms -= ndel;
	aset->bookmark -= ndel; /* is this always right? */
    } else if (op == STACK_DEPEND) {
	real_stack_set_parentage(aset->atoms, aset->n_atoms);
    } else if (op == STACK_BOOKMARK) {
	aset->bookmark = aset->n_popped;
    } else if (op == STACK_RESUME) {
	aset->n_popped = aset->bookmark;
    } else if (op == STACK_SCALAR_CHECK) {
	if (real_check_for_scalar_result(aset->atoms, aset->n_atoms)) {
	    ret = aset->atoms[0];
	}
    } 

    return ret;
}

void atom_stack_nullify_matrix (const gretl_matrix *M, GENERATOR *genr)
{
    int i;

    if (M == NULL) return;

    for (i=0; i<genr->aset->n_atoms; i++) {
	if (genr->aset->atoms[i]->M == M) {
	    genr->aset->atoms[i]->M = NULL;
	    break;
	}
    }
}

gretl_matrix *atom_stack_get_matrix (GENERATOR *genr, const char *str)
{
    gretl_matrix *M = NULL;
    int j;

    if (genr->aset == NULL) {
	return NULL;
    }

    for (j=0; j<genr->aset->n_atoms; j++) {
	if (genr->aset->atoms[j]->M != NULL && 
	    !strcmp(genr->aset->atoms[j]->str, str)) {
	    M = genr->aset->atoms[j]->M;
	    break;
	}
    }

    return M;
}

int arg_atom_available (genatom *atom)
{
    int j, ret = 0;

    for (j=0; j<atom->aset->n_atoms; j++) {
	if (atom->aset->atoms[j] == atom) {
	    if (j > 0 && atom->aset->atoms[j-1]->level == atom->level + 1) {
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

genatom *atom_stack_get_current_func (GENERATOR *genr)
{
    genatom *a = NULL;
    int j;

    if (genr->aset == NULL) {
	return NULL;
    }

    fprintf(stderr, "n_atoms = %d\n", genr->aset->n_atoms);

    for (j = genr->aset->n_atoms - 1; j >= 0; j--) {
	if (genr->aset->atoms[j]->func != 0) {
	    a = genr->aset->atoms[j];
	    break;
	}
    }

    return a;    
}

int attach_atomset (GENERATOR *genr)
{
    atomset *aset = malloc(sizeof *aset);

    if (aset == NULL) return 1;
    
    aset->atoms = NULL;
    aset->n_atoms = 0;
    aset->n_popped = 0;
    aset->bookmark = 0;

    genr->aset = aset;

    return 0;
}

int push_atom (genatom *atom)
{
    return (atom_stack(atom, atom->aset, STACK_PUSH) == NULL);
}

genatom *pop_atom (GENERATOR *genr)
{
    return atom_stack(NULL, genr->aset, STACK_POP);
}

genatom *pop_child_atom (genatom *atom)
{
    return atom_stack(atom, atom->aset, STACK_POP_CHILDREN);
}

genatom *peek_child_atom (genatom *atom)
{
    return atom_stack(atom, atom->aset, STACK_PEEK_CHILDREN);
}

void reset_atom_stack (GENERATOR *genr)
{
    atom_stack(NULL, genr->aset, STACK_RESET);
}

void destroy_atom_stack (GENERATOR *genr)
{
    if (genr->aset == NULL) return;

    atom_stack(NULL, genr->aset, STACK_DESTROY);
    free(genr->aset);
    genr->aset = NULL;
}

void atom_stack_set_parentage (GENERATOR *genr)
{
    atom_stack(NULL, genr->aset, STACK_DEPEND);
}

void atom_eat_children (genatom *atom)
{
    atom_stack(atom, atom->aset, STACK_EAT_CHILDREN);
}

void atom_stack_bookmark (GENERATOR *genr)
{
    atom_stack(NULL, genr->aset, STACK_BOOKMARK);
}

void atom_stack_resume (GENERATOR *genr)
{
    atom_stack(NULL, genr->aset, STACK_RESUME);
}

int atom_stack_check_for_scalar (GENERATOR *genr)
{
    return (atom_stack(NULL, genr->aset, STACK_SCALAR_CHECK) != NULL);
}

static gretl_matrix *
matrix_calc_stack (gretl_matrix *M, int op, GENERATOR *genr)
{
    gretl_matrix **mstack = genr->mstack;
    gretl_matrix *R = NULL;
    int i, j;

    if (op == STACK_PUSH) {
	if (genr->nmats == MATSTACK_SIZE - 1) {
	    fprintf(stderr, "genr: matrix stack depth exceeded\n");
	    genr->err = 1;
	    return R;
	} else {
	    for (i=genr->nmats; i>0; i--) {
		mstack[i] = mstack[i-1];
	    }
	    mstack[0] = M;
	    genr->nmats += 1;
	    MPRINTF(("matrix_calc_stack: STACK_PUSH: added %p, nmats = %d\n",
		     M, genr->nmats));
	}
    } else if (op == STACK_POP && genr->nmats > 0) {
	R = mstack[0];
	for (i=0; i<MATSTACK_SIZE-1; i++) {
	    mstack[i] = mstack[i+1];
	}
	genr->nmats -= 1;
	MPRINTF(("matrix_calc_stack: STACK_POP: returning %p, nmats = %d\n", R,
		 genr->nmats));
    } else if (op == STACK_RESET) {
	MPRINTF(("matrix_calc_stack: STACK_RESET\n"));
	for (i=0; i<MATSTACK_SIZE; i++) {
	    if (mstack[i] != NULL) {
		if (!is_user_matrix(mstack[i]) && !matrix_is_on_atom(mstack[i])) {
		    MPRINTF(("freeing mstack[%d] at %p\n", i, 
			     (void *) mstack[i]));
		    for (j=i+1; j<MATSTACK_SIZE; j++) {
			/* insure against double-freeing */
			if (mstack[j] == mstack[i]) {
			    mstack[j] = NULL;
			}
		    }
		    gretl_matrix_free(mstack[i]);
		}
		mstack[i] = NULL;
	    }
	}
	genr->nmats = 0;
    }

    return R;
}

static double calc_stack (double val, int op, GENERATOR *genr)
{
    double *valstack = genr->valstack;
    int i;
    double x = 0.0;

    if (op == STACK_PUSH) {
	if (genr->nvals == VALSTACK_SIZE - 1) {
	    fprintf(stderr, "genr: stack depth exceeded\n");
	    genr->err = 1;
	    return x;
	} else {
	    for (i=genr->nvals; i>0; i--) {
		valstack[i] = valstack[i-1];
	    }
	    valstack[0] = val;
	    genr->nvals += 1;
	}
    } else if (op == STACK_POP && genr->nvals > 0) {
	x = valstack[0];
	for (i=0; i<VALSTACK_SIZE-1; i++) {
	    valstack[i] = valstack[i+1];
	}
	genr->nvals -= 1;
    } else if (op == STACK_RESET) {
	for (i=0; i<VALSTACK_SIZE; i++) {
	    valstack[i] = 0.0;
	}
	genr->nvals = 0;
    }

    return x;
}

int calc_push (double x, GENERATOR *genr)
{
    calc_stack(x, STACK_PUSH, genr);
    return genr->err;
}

double calc_pop (GENERATOR *genr)
{
    return calc_stack(0., STACK_POP, genr);
}

void reset_calc_stack (GENERATOR *genr)
{
    calc_stack(0., STACK_RESET, genr);
}

int matrix_calc_push (gretl_matrix *M, GENERATOR *genr)
{
    matrix_calc_stack(M, STACK_PUSH, genr);
    return genr->err;
}

gretl_matrix *matrix_calc_pop (GENERATOR *genr)
{
    return matrix_calc_stack(NULL, STACK_POP, genr);
}

void reset_matrix_calc_stack (GENERATOR *genr)
{
    matrix_calc_stack(NULL, STACK_RESET, genr);
}

#if (GENR_DEBUG || GEN_MATRIX_DEBUG)
# include <stdarg.h>
void dprintf (const char *format, ...)
{
   va_list args;

   va_start(args, format);
   vfprintf(stderr, format, args);
   va_end(args);

   return;
}
#endif 
