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

#include "gretl.h"
#include "gfn_arglists.h"

/* apparatus for remembering the list of arguments given to
   packaged functions (gfn)
*/

struct arglist_ {
    char pkgname[32];
    const void *func;
    int argc;
    char **argv;
};

static arglist **arglists;
static int n_arglists;

static int arglist_push (arglist *a)
{
    int n = n_arglists + 1;
    arglist **A;

    A = realloc(arglists, n * sizeof *A);

    if (A == NULL) {
	return E_ALLOC;
    } else {
	arglists = A;
	arglists[n-1] = a;
	n_arglists = n;
    }

    return 0;
}

static void arglist_destroy (arglist *a)
{
    if (a != NULL) {
	strings_array_free(a->argv, a->argc);
	free(a);
    }
}

void arglist_cleanup (void)
{
    int i;
    
    for (i=0; i<n_arglists; i++) {
	arglist_destroy(arglists[i]);
    }

    free(arglists);
    arglists = NULL;
    n_arglists = 0;
}

arglist *arglist_new (const char *pkgname, const void *func, int argc)
{
    arglist *a = malloc(sizeof *a);

    if (a != NULL) {
	a->argv = strings_array_new(argc);
	if (a->argv == NULL) {
	    free(a);
	    a = NULL;
	} else {
	    int err;
	    
	    *a->pkgname = '\0';
	    strncat(a->pkgname, pkgname, 31);
	    a->func = func;
	    a->argc = argc;
	    err = arglist_push(a);
	    if (err) {
		arglist_destroy(a);
		a = NULL;
	    }
	}
    }

    return a;
}

arglist *arglist_lookup (const char *pkgname, const void *func)
{
    int i;

    for (i=0; i<n_arglists; i++) {
	if (!strcmp(arglists[i]->pkgname, pkgname) &&
	    arglists[i]->func == func) {
	    return arglists[i];
	}
    }

    return NULL;
}

int arglist_record_arg (arglist *a, int i, const char *val)
{
    if (a == NULL || i < 0 || i > a->argc) {
	return E_DATA;
    } else {
	free(a->argv[i]);
	if (strcmp(val, AUTOLIST)) {
	    a->argv[i] = gretl_strdup(val);
	} else {
	    a->argv[i] = NULL;
	}
	return 0;
    }
}

const char *arglist_lookup_val (arglist *a, int i)
{
    if (a == NULL || i < 0 || i > a->argc) {
	return NULL;
    } else {
	return a->argv[i];
    }
}

