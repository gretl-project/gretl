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

#define must_know_children(f) (f == T_NOBS || f == T_MEAN || f == T_SUM || \
                               f == T_SD || f == T_VAR || f == T_SST || \
                               f == T_MEDIAN || f == T_MIN || f == T_MAX || \
                               f == T_DIFF || f == T_LDIFF || \
                               f == T_CUM || f == T_SORT || \
                               f == T_RESAMPLE || f == T_HPFILT)

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

static genatom *atom_stack (genatom *atom, int op)
{
    static genatom **atoms;
    static int n_atoms;
    static int n_popped;
    static int bookmark;
    genatom *ret = NULL;
    int j;

    if (op == STACK_PUSH && atom != NULL) {
	atoms = realloc(atoms, (n_atoms + 1) * sizeof *atoms);
	if (atoms != NULL) {
	    atoms[n_atoms++] = atom;
	    ret = atom;
	}
    } else if (op == STACK_POP) {
	if (n_popped < n_atoms) {
	    ret = atoms[n_popped++];
	}
    } else if (op == STACK_RESET) {
	for (j=0; j<n_atoms; j++) (atoms[j])->popped = 0;
	n_popped = 0;
    } else if (op == STACK_DESTROY) {
	for (j=0; j<n_atoms; j++) free(atoms[j]);
	free(atoms);
	atoms = NULL;
	n_atoms = 0;
	n_popped = 0;
    } else if (op == STACK_POP_CHILDREN && atom != NULL) {
	for (j=0; j<n_atoms; j++) {
	    if ((atoms[j])->parent == atom && !(atoms[j])->popped) {
		ret = atoms[j];
		(atoms[j])->popped = 1;
		break;
	    }
	}
    } else if (op == STACK_PEEK_CHILDREN && atom != NULL) {
	for (j=0; j<n_atoms; j++) {
	    if ((atoms[j])->parent == atom) {
		ret = atoms[j];
		break;
	    }
	}	
    } else if (op == STACK_EAT_CHILDREN && atom != NULL) {
	int ndel = real_stack_eat_children(atom, atoms, n_atoms);

	n_atoms -= ndel;
	bookmark -= ndel; /* is this always right? */
    } else if (op == STACK_DEPEND) {
	real_stack_set_parentage(atoms, n_atoms);
    } else if (op == STACK_BOOKMARK) {
	bookmark = n_popped;
    } else if (op == STACK_RESUME) {
	n_popped = bookmark;
    } else if (op == STACK_SCALAR_CHECK) {
	if (real_check_for_scalar_result(atoms, n_atoms)) {
	    ret = atoms[0];
	}
    }

    return ret;
}

int push_atom (genatom *atom)
{
    return (atom_stack(atom, STACK_PUSH) == NULL);
}

genatom *pop_atom (void)
{
    return atom_stack(NULL, STACK_POP);
}

genatom *pop_child_atom (genatom *atom)
{
    return atom_stack(atom, STACK_POP_CHILDREN);
}

genatom *peek_child_atom (genatom *atom)
{
    return atom_stack(atom, STACK_PEEK_CHILDREN);
}

void reset_atom_stack (void)
{
    atom_stack(NULL, STACK_RESET);
}

void destroy_atom_stack (void)
{
    atom_stack(NULL, STACK_DESTROY);
}

void atom_stack_set_parentage (void)
{
    atom_stack(NULL, STACK_DEPEND);
}

void atom_eat_children (genatom *atom)
{
    atom_stack(atom, STACK_EAT_CHILDREN);
}

void atom_stack_bookmark (void)
{
    atom_stack(NULL, STACK_BOOKMARK);
}

void atom_stack_resume (void)
{
    atom_stack(NULL, STACK_RESUME);
}

int atom_stack_check_for_scalar (void)
{
    return (atom_stack(NULL, STACK_SCALAR_CHECK) != NULL);
}

#define STACKSIZE 32

static double calc_stack (double val, int op, int *err)
{
    static double valstack[STACKSIZE]; /* how big should this be? */
    static int nvals;
    int i;
    double x = 0.0;

    if (op == STACK_PUSH) {
	if (nvals == STACKSIZE - 1) {
	    fprintf(stderr, "genr: stack depth exceeded\n");
	    *err = 1;
	    return x;
	} else {
	    for (i=nvals; i>0; i--) {
		valstack[i] = valstack[i-1];
	    }
	    valstack[0] = val;
	    nvals++;
	}
    }
    else if (op == STACK_POP && nvals > 0) {
	x = valstack[0];
	for (i=0; i<STACKSIZE-1; i++) {
	    valstack[i] = valstack[i+1];
	}
	nvals--;
    }
    else if (op == STACK_RESET) {
	for (i=0; i<STACKSIZE; i++) {
	    valstack[i] = 0.0;
	}
	nvals = 0;
    }

    return x;
}

int calc_push (double x)
{
    int err = 0;

    calc_stack(x, STACK_PUSH, &err);
    return err;
}

double calc_pop (void)
{
    return calc_stack(0., STACK_POP, NULL);
}

void reset_calc_stack (void)
{
    calc_stack(0., STACK_RESET, NULL);
}

#ifdef GENR_DEBUG
# include <stdarg.h>
void dprintf (const char *format, ...)
{
   va_list args;

   va_start(args, format);
   vfprintf(stderr, format, args);
   va_end(args);

   return;
}
# define DPRINTF(x) dprintf x
#else 
# define DPRINTF(x)
#endif /* GENR_DEBUG */
