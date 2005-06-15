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
                               f == T_DIFF || f == T_LDIFF || f == T_T1 || f == T_T2 || \
                               f == T_CUM || f == T_SORT || \
                               f == T_VARNUM || f == T_VECTOR || \
                               f == T_ISLIST || f == T_NELEM || \
                               f == T_RESAMPLE || f == T_HPFILT || f == T_BKFILT)

static int all_children_scalar (genatom **atoms, int n,
				int pos, int level)
{
    int j, sc = 1;

    if (pos == 0) return 0;

    for (j=pos-1; j>=0; j--) {
	/* below: unsure between "!= level + 1" and "<= level" */
	if ((atoms[j])->level <= level) { 
	    if (j == pos - 1) sc = 0;
	    break;
	}
	if (!(atoms[j])->scalar) {
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
	if ((atoms[i])->level > maxlevel) {
	    maxlevel = (atoms[i])->level;
	}
    }

    DPRINTF(("check_for_scalar: maxlevel = %d\n", maxlevel));

    /* recurse to identify scalar property of composites */
    for (k=maxlevel-1; k>=0; k--) {
	DPRINTF(("check_for_scalar: checking level %d\n", k));
	for (i=0; i<n; i++) {
	    if ((atoms[i])->level == k && !(atoms[i])->scalar) {
		DPRINTF(("check_for_scalar: checking atom %d\n", i));
		if (all_children_scalar(atoms, n, i, k)) {
		    DPRINTF(("check_for_scalar: setting scalar=-1 "
			     "for atom[%d]\n", i));
		    (atoms[i])->scalar = -1; 
		}
	    }
	}
    }

    for (i=0; i<n; i++) {
	if ((atoms[i])->level == 0 && !(atoms[i])->scalar) {
	    scalar = 0;
	}
	if ((atoms[i])->scalar == -1) {
	    /* restore original property */
	    (atoms[i])->scalar = 0;
	}
    }

    return scalar;
}

static int real_stack_eat_children (genatom *parent, 
				    genatom **atoms, int n)
{
    int i, j, ndel = 0;

    for (i=0; i<n; i++) {
	if ((atoms[i])->parent == parent) {
	    free(atoms[i]);
	    DPRINTF(("freed child atom, pos %d\n", i));
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
	if (must_know_children((atoms[i])->func)) {
	    level = (atoms[i])->level;
	    DPRINTF((" got candidate (atom %d, level %d)...\n", i, level));
	    for (j=i-1; j>=0; j--) {
		DPRINTF(("  looking at atom %d\n", j));
		if ((atoms[j])->level > level) {
		    DPRINTF(("    marking atom %d as parent of atom %d\n",
			    i, j));
		    (atoms[j])->parent = atoms[i];
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

/* below: ths static stuff is a big problem if genr is
   called within aa genr. 
*/

static genatom *atom_stack (genatom *atom, atomset *aset, int op)
{
    genatom *ret = NULL;
    int j;

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
	    (aset->atoms[j])->popped = 0;
	}
	aset->n_popped = 0;
    } else if (op == STACK_DESTROY) {
	for (j=0; j<aset->n_atoms; j++) {
	    free(aset->atoms[j]);
	}
	free(aset->atoms);
	aset->atoms = NULL;
	aset->n_atoms = 0;
	aset->n_popped = 0;
    } else if (op == STACK_POP_CHILDREN && atom != NULL) {
	for (j=0; j<aset->n_atoms; j++) {
	    if ((aset->atoms[j])->parent == atom && !(aset->atoms[j])->popped) {
		ret = aset->atoms[j];
		(aset->atoms[j])->popped = 1;
		break;
	    }
	}
    } else if (op == STACK_PEEK_CHILDREN && atom != NULL) {
	for (j=0; j<aset->n_atoms; j++) {
	    if ((aset->atoms[j])->parent == atom) {
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

int attach_atomset (GENERATE *genr)
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

genatom *pop_atom (GENERATE *genr)
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

void reset_atom_stack (GENERATE *genr)
{
    atom_stack(NULL, genr->aset, STACK_RESET);
}

void destroy_atom_stack (GENERATE *genr)
{
    if (genr->aset == NULL) return;

    atom_stack(NULL, genr->aset, STACK_DESTROY);
    free(genr->aset);
    genr->aset = NULL;
}

void atom_stack_set_parentage (GENERATE *genr)
{
    atom_stack(NULL, genr->aset, STACK_DEPEND);
}

void atom_eat_children (genatom *atom)
{
    atom_stack(atom, atom->aset, STACK_EAT_CHILDREN);
}

void atom_stack_bookmark (GENERATE *genr)
{
    atom_stack(NULL, genr->aset, STACK_BOOKMARK);
}

void atom_stack_resume (GENERATE *genr)
{
    atom_stack(NULL, genr->aset, STACK_RESUME);
}

int atom_stack_check_for_scalar (GENERATE *genr)
{
    return (atom_stack(NULL, genr->aset, STACK_SCALAR_CHECK) != NULL);
}

static double calc_stack (double val, int op, GENERATE *genr)
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
    }
    else if (op == STACK_POP && genr->nvals > 0) {
	x = valstack[0];
	for (i=0; i<VALSTACK_SIZE-1; i++) {
	    valstack[i] = valstack[i+1];
	}
	genr->nvals -= 1;
    }
    else if (op == STACK_RESET) {
	for (i=0; i<VALSTACK_SIZE; i++) {
	    valstack[i] = 0.0;
	}
	genr->nvals = 0;
    }

    return x;
}

int calc_push (double x, GENERATE *genr)
{
    calc_stack(x, STACK_PUSH, genr);
    return genr->err;
}

double calc_pop (GENERATE *genr)
{
    return calc_stack(0., STACK_POP, genr);
}

void reset_calc_stack (GENERATE *genr)
{
    calc_stack(0., STACK_RESET, genr);
}

#if GENR_DEBUG
# include <stdarg.h>
void dprintf (const char *format, ...)
{
   va_list args;

   va_start(args, format);
   vfprintf(stderr, format, args);
   va_end(args);

   return;
}
#endif /* GENR_DEBUG */
