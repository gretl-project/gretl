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

struct gretl_scalar_ {
    char name[VNAMELEN];
    double val;
    int level;
};

static gretl_scalar **scalars;
static int n_scalars;

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

static double
real_get_scalar_by_name (const char *name, int slevel,
			 int *err)
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

int gretl_is_scalar (const char *name)
{
    int level = gretl_function_depth();
    int i;

    for (i=0; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    return 1;
	}
    }

    return 0;
}

double gretl_scalar_get_value (const char *name, int *err)
{
    int level = gretl_function_depth();
    int i;

    for (i=0; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    return scalars[i]->val;
	}
    }

    *err = E_UNKVAR;

    return NADBL;
}

void gretl_scalar_set_value (const char *name, double val)
{
    int level = gretl_function_depth();
    int i;

    for (i=0; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    scalars[i]->val = val;
	    break;
	}
    }
}

int gretl_scalar_add (const char *name, double val)
{
    gretl_scalar *s;

    s = malloc(sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    strcpy(s->name, name);
    s->val = val;
    s->level = gretl_function_depth();

    return gretl_scalar_push(s);
}

int gretl_scalar_delete (const char *name)
{
    int level = gretl_function_depth();
    int i, ret = E_UNKVAR;

    for (i=0; i<n_scalars; i++) {
	if (!strcmp(name, scalars[i]->name) &&
	    scalars[i]->level == level) {
	    ret = real_delete_scalar(i);
	    break;
	}
    }

    return ret;
}




