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

#include "libgretl.h"
#include "gretl_func.h"
#include "gretl_scalar.h"

#define SDEBUG 0

typedef struct gretl_scalar_ gretl_scalar;

struct gretl_scalar_ {
    char name[VNAMELEN];
    double val;
    int level;
};

static gretl_scalar **scalars;
static int n_scalars;
static int scalar_imin;

#if SDEBUG

static void print_scalars (const char *s)
{
    if (n_scalars == 0) {
	fprintf(stderr, "%s: no scalars defined\n", s);
    } else {
	int i;

	fprintf(stderr, "%s: defined scalars:\n", s);
	for (i=0; i<n_scalars; i++) {
	    fprintf(stderr, "%15s = %g\n", scalars[i]->name,
		    scalars[i]->val);
	}
    }
}

#endif

static gretl_scalar *get_scalar_pointer (const char *s, int level)
{
    int i;

    for (i=scalar_imin; i<n_scalars; i++) {
	if (!strcmp(s, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    return scalars[i];
	}
    }

    return NULL;
}

static int gretl_scalar_push (gretl_scalar *s)
{
    gretl_scalar **tmp;
    int n = n_scalars + 1;

    tmp = realloc(scalars, n * sizeof *tmp);
    if (tmp == NULL) {
	free(s);
	return E_ALLOC;
    }

    scalars = tmp;
    scalars[n-1] = s;
    n_scalars++;

    return 0;
}

static int real_delete_scalar (int i)
{
    int n = n_scalars - 1;
    int err = 0;

    free(scalars[i]);

    if (n == 0) {
	free(scalars);
	scalars = NULL;
	n_scalars = 0;
    } else {
	gretl_scalar **tmp;
	int j;

	for (j=i; j<n; j++) {
	    scalars[j] = scalars[j+1];
	}
	tmp = realloc(scalars, n * sizeof *tmp);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    scalars = tmp;
	    n_scalars = n;
	}
    }

    return err;
}

#if 0

static double
real_get_scalar_by_name (const char *name, int slevel, int *err)
{
    int level, i;

    if (slevel == LEVEL_AUTO) {
	level = gretl_function_depth();
    } else {
	level = slevel;
    }

    for (i=0; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    return scalars[i]->val;
	}
    }

    *err = E_UNKVAR;

    return NADBL;
}

#endif

int gretl_is_scalar (const char *name)
{
#if SDEBUG
    print_scalars("gretl_is_scalar");
#endif

    return (get_scalar_pointer(name, gretl_function_depth()) != NULL);
}

double gretl_scalar_get_value (const char *name, int *err)
{
    gretl_scalar *s;

#if SDEBUG
    print_scalars("gretl_scalar_get_value");
#endif

    s = get_scalar_pointer(name, gretl_function_depth());
    if (s != NULL) {
	return s->val;
    } 

    if (err != NULL) {
	*err = E_UNKVAR;
    }

    return NADBL;
}

void gretl_scalar_set_value (const char *name, double val)
{
    gretl_scalar *s;

    s = get_scalar_pointer(name, gretl_function_depth());
    if (s != NULL) {
	s->val = val;
    }

#if SDEBUG
    print_scalars("gretl_scalar_set_value");
#endif
}

int gretl_scalar_add (const char *name, double val)
{
    gretl_scalar *s;
    int err;

    s = malloc(sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    strcpy(s->name, name);
    s->val = val;
    s->level = gretl_function_depth();

    err = gretl_scalar_push(s);

#if SDEBUG
    print_scalars("gretl_scalar_add");
#endif

    return err;
}

int gretl_scalar_delete (const char *name)
{
    int level = gretl_function_depth();
    int i, ret = E_UNKVAR;

    for (i=scalar_imin; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    ret = real_delete_scalar(i);
	    break;
	}
    }

#if SDEBUG
    print_scalars("gretl_scalar_delete");
#endif

    return ret;
}

/* "auxiliary scalars": this apparatus is used when we want to do
   "private" NLS estimation (e.g. in ARMA initialization).  It ensures
   that the scalar NLS parameters don't collide with the public scalar
   namespace.
*/

void set_auxiliary_scalars (void)
{
    scalar_imin = n_scalars;
}

void unset_auxiliary_scalars (void)
{
    if (scalar_imin == 0) {
	destroy_user_scalars();
    } else {
	gretl_scalar **tmp;
	int i;
    
	for (i=scalar_imin; i<n_scalars; i++) {
	    free(scalars[i]);
	}

	tmp = realloc(scalars, scalar_imin * sizeof *tmp);
	if (tmp != NULL) {
	    scalars = tmp;
	}
    }

    n_scalars = scalar_imin;
    scalar_imin = 0;
}

#define LEV_PRIVATE -1

static int levels_match (gretl_scalar *s, int lev)
{
    if (s->level == lev) {
	return 1;
    } else if (lev == LEV_PRIVATE && *s->name == '$') {
	return 1;
    }

    return 0;
}

int destroy_user_scalars_at_level (int level)
{
    gretl_scalar **tmp;
    int i, j, ns = 0;
    int err = 0;

    for (i=scalar_imin; i<n_scalars; i++) {
	if (scalars[i] == NULL) {
	    break;
	}
	if (levels_match(scalars[i], level)) {
	    free(scalars[i]);
	    for (j=i; j<n_scalars - 1; j++) {
		scalars[j] = scalars[j+1];
	    }
	    scalars[n_scalars - 1] = NULL;
	    i--;
	} else {
	    ns++;
	}
    }

    if (ns < n_scalars) {
	n_scalars = ns;
	if (ns == 0) {
	    free(scalars);
	    scalars = NULL;
	} else {
	    tmp = realloc(scalars, ns * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		scalars = tmp;
	    }
	}
    }

    return err;
}

/**
 * destroy_user_scalars:
 *
 * Frees all resources associated with the stack of user-
 * defined scalars.
 */

void destroy_user_scalars (void)
{
    int i;

    if (scalars == NULL) {
	return;
    }

    for (i=0; i<n_scalars; i++) {
	free(scalars[i]);
    }

    free(scalars);
    scalars = NULL;
    n_scalars = 0;
}





