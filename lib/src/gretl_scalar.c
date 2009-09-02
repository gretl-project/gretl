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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "gretl_func.h"
#include "libset.h"
#include "gretl_xml.h"
#include "gretl_string_table.h"
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

static void debug_print_scalars (const char *s)
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

static void set_n_scalars (int n)
{
#if SDEBUG
    fprintf(stderr, "n_scalars: %d -> %d\n", n_scalars, n);
#endif   

    n_scalars = n;

    if (n_scalars == 0) {
	free(scalars);
	scalars = NULL;
    }
}

static void free_scalar (int i)
{
#if SDEBUG
    fprintf(stderr, "n_scalars = %d, freeing scalars[%d]: ", n_scalars, i);
    fprintf(stderr, "'%s' = %g, level %d\n", scalars[i]->name, 
	    scalars[i]->val, scalars[i]->level);
#endif    

    free(scalars[i]);
    scalars[i] = NULL;
}

static int reallocate_scalars (int n)
{
    gretl_scalar **tmp;

    tmp = realloc(scalars, n * sizeof *tmp);

    if (tmp == NULL) {
	return E_ALLOC;
    } else {
	scalars = tmp;
	return 0;
    }
}

static gretl_scalar *get_scalar_pointer (const char *s, int level)
{
    int i;

    for (i=scalar_imin; i<n_scalars; i++) {
	if (scalars[i]->level == level &&
	    !strcmp(s, scalars[i]->name)) {
	    return scalars[i];
	}
    }

    return NULL;
}

static int gretl_scalar_push (gretl_scalar *s)
{
    int n = n_scalars + 1;

    if (reallocate_scalars(n)) {
	free(s);
	return E_ALLOC;
    }

    scalars[n-1] = s;
    set_n_scalars(n);

    return 0;
}

static int real_delete_scalar (int i)
{
    int n = n_scalars - 1;
    int err = 0;

    free_scalar(i);

    if (n == 0) {
	set_n_scalars(0);
    } else {
	int j;

	for (j=i; j<n; j++) {
	    scalars[j] = scalars[j+1];
	}
	if (reallocate_scalars(n)) {
	    err = E_ALLOC;
	} else {
	    set_n_scalars(n);
	}
    }

    return err;
}

int n_saved_scalars (void)
{
    return n_scalars;
}

int gretl_is_scalar (const char *name)
{
#if SDEBUG
    debug_print_scalars("gretl_is_scalar");
#endif

    if (name == NULL || *name == '\0') {
	return 0;
    }

    return (get_scalar_pointer(name, gretl_function_depth()) != NULL);
}

int gretl_scalar_get_index (const char *name, int *err)
{
    int level = gretl_function_depth();
    int i;

#if SDEBUG
    debug_print_scalars("gretl_scalar_get_index");
#endif

    for (i=0; i<n_scalars; i++) {
	if (level == scalars[i]->level &&
	    !strcmp(name, scalars[i]->name)) {
	    return i;
	}
    }

    *err = E_UNKVAR;

    return -1;
}

const char *gretl_scalar_get_name (int i)
{
    if (i >= 0 && i < n_scalars) {
	return scalars[i]->name;
    } else {
	return NULL;
    }
}

double gretl_scalar_get_value_by_index (int i)
{
    if (i >= 0 && i < n_scalars) {
	return scalars[i]->val;
    } else {
	return NADBL;
    }
}

double gretl_scalar_get_value (const char *name)
{
    gretl_scalar *s;

#if SDEBUG
    debug_print_scalars("gretl_scalar_get_value");
#endif

    s = get_scalar_pointer(name, gretl_function_depth());

    return (s != NULL)? s->val : NADBL;
}

void gretl_scalar_set_value (const char *name, double val)
{
    gretl_scalar *s;

    s = get_scalar_pointer(name, gretl_function_depth());
    if (s != NULL) {
	s->val = val;
    }

#if SDEBUG
    debug_print_scalars("gretl_scalar_set_value");
#endif
}

/* check that we're not colliding with a pre-existing object of
   different type */

static int gretl_scalar_check_name (const char *s, const DATAINFO *pdinfo)
{
    int v, err = 0;

    if ((v = series_index(pdinfo, s)) < pdinfo->v) {
	err = E_TYPES;
    } else if (get_matrix_by_name(s)) {
	err = E_TYPES;
    } else if (get_list_by_name(s)) {
	err = E_TYPES;
    } else if (get_string_by_name(s)) {
	err = E_TYPES;
    } else {
	err = check_varname(s);
    }

    return err;
}	

int gretl_scalar_add (const char *name, double val)
{
    int level = gretl_function_depth();
    gretl_scalar *s;
    int err;

    s = get_scalar_pointer(name, level);
    if (s != NULL) {
	fprintf(stderr, "*** gretl_scalar_add: there's already a '%s' at level %d (%.15g)\n", 
		name, s->level, s->val);
	s->val = val;
	return 0;
    }

    s = malloc(sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    strcpy(s->name, name);
    s->val = val;
    s->level = gretl_function_depth();

    err = gretl_scalar_push(s);

#if SDEBUG
    debug_print_scalars("gretl_scalar_add");
#endif

    return err;
}

int gretl_scalar_add_with_check (const char *name, double val,
				 const DATAINFO *pdinfo)
{
    int err = gretl_scalar_check_name(name, pdinfo);

    if (!err) {
	err = gretl_scalar_add(name, val);
    }

    return err;
}

int gretl_scalar_add_as_arg (const char *name, double val)
{
    gretl_scalar *s;
    int err;

    s = malloc(sizeof *s);
    if (s == NULL) {
	return E_ALLOC;
    }

    strcpy(s->name, name);
    s->val = val;
    s->level = gretl_function_depth() + 1;

    err = gretl_scalar_push(s);

#if SDEBUG
    debug_print_scalars("gretl_scalar_add_as");
#endif

    return err;
}

int gretl_scalar_set_local_name (int i, const char *name)
{
    if (i >= 0 && i < n_scalars) {
	scalars[i]->name[0] = '\0';
	strncat(scalars[i]->name, name, VNAMELEN - 1);
	scalars[i]->level += 1;
	return 0;
    } else {
	return E_DATA;
    }
}

int gretl_scalar_restore_name (int i, const char *name)
{
    if (i >= 0 && i < n_scalars) {
	scalars[i]->name[0] = '\0';
	strncat(scalars[i]->name, name, VNAMELEN - 1);
	scalars[i]->level -= 1;
	return 0;
    } else {
	return E_DATA;
    }
}

int gretl_scalar_delete (const char *name, PRN *prn)
{
    int level = gretl_function_depth();
    int i, err = E_UNKVAR;

    for (i=scalar_imin; i<n_scalars; i++) {
	if (scalars[i]->level == level && 
	    !strcmp(name, scalars[i]->name)) {
	    err = real_delete_scalar(i);
	    break;
	}
    }

    if (!err && prn != NULL && gretl_messages_on()) {
	pprintf(prn, _("Deleted scalar %s"), name);
	pputc(prn, '\n');
    }

#if SDEBUG
    debug_print_scalars("gretl_scalar_delete");
#endif

    return err;
}

static void print_scalar (gretl_scalar *s, gretlopt opt, PRN *prn, int n)
{
    if (n == 0) {
	pputc(prn, '\n');
    }

    pprintf(prn, "%15s = ", s->name);

    if (na(s->val)) {
	pputs(prn, "NA");
    } else {
	if (s->val >= 0.0) {
	    pputc(prn, ' ');
	}
	if (opt & OPT_L) {
	    pprintf(prn, "%#.*E", libset_get_int(LONGDIGITS), s->val);
	} else if (opt & OPT_T) {
	    pprintf(prn, "%#.10E", s->val);
	} else {
	    pprintf(prn, "%#.6g", s->val);
	}
    }

    pputc(prn, '\n');
}

void print_scalar_by_name (const char *name, gretlopt opt, PRN *prn)
{
    gretl_scalar *s = get_scalar_pointer(name, gretl_function_depth());

    if (s != NULL) {
	print_scalar(s, opt, prn, 0);
    }
}

void print_all_scalars (gretlopt opt, PRN *prn)
{
    int n = 0;

    if (n_scalars > 0) {
	int i, level = gretl_function_depth();

	for (i=0; i<n_scalars; i++) {
	    if (scalars[i]->level == level) {
		print_scalar(scalars[i], opt, prn, n++);
	    }
	}
    }

    if (n == 0) {
	pputs(prn, "No scalars are defined\n");
    }
}

/* "auxiliary scalars": this apparatus is used when we want to do
   "private" NLS estimation (e.g. in ARMA initialization).  It ensures
   that the scalar NLS parameters don't collide with the public scalar
   namespace.
*/

void set_auxiliary_scalars (void)
{
#if SDEBUG
    fprintf(stderr, "set_auxiliary_scalars: imin = %d\n", n_scalars);
#endif
    scalar_imin = n_scalars;
}

void unset_auxiliary_scalars (void)
{
#if SDEBUG
    fprintf(stderr, "unset_auxiliary_scalars: imin = %d, n_scalars = %d\n", 
	    scalar_imin, n_scalars);
#endif

    if (scalar_imin == 0) {
	destroy_user_scalars();
    } else {
	int i;
    
	for (i=scalar_imin; i<n_scalars; i++) {
	    free_scalar(i);
	}

	reallocate_scalars(scalar_imin);
    }

    set_n_scalars(scalar_imin);
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

/**
 * destroy_user_scalars_at_level.
 * @level: depth of function execution.
 *
 * Deletes any saved scalars at a particular level of function
 * execution; used on exiting a user-defined function.
 * 
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_scalars_at_level (int level)
{
    int ns = n_scalars;
    int smax = ns - 1;
    int i, j;
    int err = 0;

#if SDEBUG
    fprintf(stderr, "destroy_user_scalars_at_level: level = %d\n", level);
#endif

    for (i=smax; i>=scalar_imin; i--) {
	if (levels_match(scalars[i], level)) {
	    free_scalar(i);
	    for (j=i; j<smax; j++) {
		scalars[j] = scalars[j+1];
	    }
	    scalars[smax] = NULL;
	    ns--;
	} 
    }

    if (ns < n_scalars) {
	set_n_scalars(ns);
	if (ns > 0) {
	    err = reallocate_scalars(ns);
	}
    }

    return err;
}

/**
 * destroy_private_scalars:
 *
 * Gets rid of private or "internal" scalars whose
 * names begin with '$'.
 */

void destroy_private_scalars (void)
{
    destroy_user_scalars_at_level(LEV_PRIVATE);
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

#if SDEBUG
    fprintf(stderr, "destroy_user_scalars: scalars = %p, n_scalars = %d\n", 
	    (void *) scalars, n_scalars);
#endif

    for (i=0; i<n_scalars; i++) {
	free_scalar(i);
    }

    set_n_scalars(0);
}

/**
 * write_scalars_to_file:
 * @fp: stream to which to write.
 *
 * Prints information on any saved scalars as XML, for use
 * when saving a gretl session.
 */

void write_scalars_to_file (FILE *fp)
{
    int i;

    gretl_xml_header(fp);
    fputs("<gretl-scalars>\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n_scalars; i++) {
	fprintf(fp, " <gretl-scalar name=\"%s\" value=\"%.15g\"/>\n", 
		scalars[i]->name, scalars[i]->val);
    }

    gretl_pop_c_numeric_locale();

    fputs("</gretl-scalars>\n", fp);
}
