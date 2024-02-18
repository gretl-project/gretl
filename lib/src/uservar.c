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
  */

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "usermat.h"
#include "gretl_string_table.h"
#include "matrix_extra.h"
#include "libset.h"
#include "monte_carlo.h"
#include "gretl_typemap.h"
#include "uservar.h"
#include "uservar_priv.h"
#include "gretl_cmatrix.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#define UVDEBUG 0
#define HDEBUG 0

#if HDEBUG && defined(_OPENMP)
# include <omp.h>
#endif

#define LEVEL_AUTO -1
#define LEV_PRIVATE -1

static user_var **uvars;
static int n_vars;
static int n_alloc;
static int scalar_imin;

/* callback for the benefit of the edit scalars window
   in the gretl GUI */

static void (*scalar_edit_callback)(void);

/* callback for adding or deleting icons representing
   things in the GUI session window */

static USER_VAR_FUNC user_var_callback;

#define UV_CHUNK 32

#define var_is_private(u) ((u->flags & UV_PRIVATE) || *u->name == '$' || *u->name == '_')
#define var_is_shell(u)   (u->flags & UV_SHELL)

static double *na_ptr (void)
{
    double *px = malloc(sizeof *px);

    if (px != NULL) {
        *px = NADBL;
    }

    return px;
}

static user_var *user_var_new (const char *name, GretlType type,
                               void *value, int *err)
{
    user_var *u;

    if (type == GRETL_TYPE_NONE) {
        *err = E_DATA;
        fputs("user_var_new: type = GRETL_TYPE_NONE\n", stderr);
        return NULL;
    }

    u = malloc(sizeof *u);

    if (u == NULL) {
        *err = E_ALLOC;
    } else {
        u->type = type;
        u->level = gretl_function_depth();
        u->flags = (u->level == 0)? UV_MAIN : 0;
        *u->name = '\0';
        strncat(u->name, name, VNAMELEN - 1);
        u->ptr = NULL;

        if (type == GRETL_TYPE_MATRIX) {
            gretl_matrix *m = value;

            if (m == NULL) {
                u->ptr = gretl_null_matrix_new();
            } else if (get_user_var_by_data(m) != NULL) {
                /* this check should be redundant? */
                u->ptr = gretl_matrix_copy(m);
            } else {
                u->ptr = value;
            }
        } else if (type == GRETL_TYPE_BUNDLE) {
            if (value == NULL) {
                u->ptr = gretl_bundle_new();
            } else {
                u->ptr = value;
            }
        } else if (type == GRETL_TYPE_STRING) {
            if (value == NULL) {
                u->ptr = gretl_strdup("");
            } else {
                u->ptr = value;
            }
        } else if (type == GRETL_TYPE_LIST) {
            if (value == NULL) {
                u->ptr = gretl_null_list();
            } else {
                u->ptr = value;
            }
        } else if (type == GRETL_TYPE_DOUBLE) {
            if (value == NULL) {
                u->ptr = na_ptr();
            } else {
                u->ptr = value;
            }
        } else if (gretl_array_type(type) || type == GRETL_TYPE_ANY) {
            if (value == NULL) {
                u->ptr = gretl_array_new(type, 0, err);
            } else {
                u->ptr = value;
            }
            u->type = GRETL_TYPE_ARRAY;
        } else {
            fprintf(stderr, "user_var_new error, type=%d (%s)\n", type,
                    gretl_type_get_name(type));
            *err = E_DATA;
        }
    }

    if (u->ptr == NULL) {
        if (!*err) {
            *err = E_ALLOC;
        }
        free(u);
        u = NULL;
    }

    return u;
}

static void uvar_free_value (user_var *u)
{
    if (u->ptr == NULL) {
        return;
    } else if (u->type == GRETL_TYPE_MATRIX) {
        gretl_matrix_free(u->ptr);
    } else if (u->type == GRETL_TYPE_BUNDLE) {
        gretl_bundle_destroy(u->ptr);
    } else if (u->type == GRETL_TYPE_STRING) {
        bufgets_finalize(u->ptr);
        free(u->ptr);
    } else if (u->type == GRETL_TYPE_ARRAY) {
        gretl_array_destroy(u->ptr);
    } else {
        /* scalar, list */
        free(u->ptr);
    }
}

static GHashTable *uvh0;       /* for use at "main" exec level */
static GHashTable *uvh1;       /* for use within functions */
static GHashTable *uvars_hash; /* pointer to one or other of the above */
static int previous_d = -1;    /* record of previous "function depth" */

void set_previous_depth (int d)
{
    previous_d = d;
}

static int get_previous_depth (void)
{
    return previous_d;
}

void switch_uservar_hash (int level)
{
#if HDEBUG
    fprintf(stderr, "switch_uservar_hash: going to %s (level %d, recursing %d)\n",
            level == 0 ? "uvh0" : "uvh1", level, recursing);
    fprintf(stderr, " uvh0 = %p, uvh1 = %p\n", (void *) uvh0, (void *) uvh1);
#endif

    if (level == 0) {
        uvars_hash = uvh0;
        if (uvh1 != NULL) {
            g_hash_table_remove_all(uvh1);
        }
    } else {
        uvars_hash = uvh1;
    }
}

static void uvar_hash_destroy (void)
{
    if (uvh0 != NULL) {
#if HDEBUG
        fprintf(stderr, "uvar_hash_destroy: destroying uvh0\n");
#endif
        g_hash_table_destroy(uvh0);
        uvh0 = NULL;
    }

    if (uvh1 != NULL) {
#if HDEBUG
        fprintf(stderr, "uvar_hash_destroy: destroying uvh1\n");
#endif
        g_hash_table_destroy(uvh1);
        uvh1 = NULL;
    }

    /* also NULL the convenience pointer */
    uvars_hash = NULL;

    set_previous_depth(-1);
}

static void user_var_destroy (user_var *u)
{
#if HDEBUG > 1
    fprintf(stderr, "user_var_destroy: '%s' (level %d)\n", u->name, u->level);
#endif

    if (uvars_hash != NULL) {
# if HDEBUG > 1
        if (g_hash_table_remove(uvars_hash, u->name)) {
            fprintf(stderr, "removed '%s' from hash table at %p\n",
                    u->name, (void *) uvars_hash);
        }
# else
        g_hash_table_remove(uvars_hash, u->name);
# endif
    }

    if (!var_is_shell(u)) {
        uvar_free_value(u);
    }

    free(u);
}

static int resize_uvar_stack (int n)
{
    int err = 0;

    if (n > n_alloc) {
        int n_new = n_alloc + UV_CHUNK;
        user_var **tmp;

        tmp = realloc(uvars, n_new * sizeof *tmp);
        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            uvars = tmp;
            n_alloc = n_new;
        }
    }

    return err;
}

static void set_nvars (int n, const char *caller)
{
#if UVDEBUG
    fprintf(stderr, "%s: setting n_vars = %d (was %d)\n",
            caller, n, n_vars);
#endif
    n_vars = n;
}

static int bname_is_temp (const char *name)
{
    return !strncmp(name, "btmp___", 7) && isdigit(name[7]);
}

static int real_user_var_add (const char *name,
                              GretlType type,
                              void *value,
                              gretlopt opt,
			      user_var **pu)
{
    user_var *u;
    int err = 0;

    u = user_var_new(name, type, value, &err);

    if (u == NULL) {
        fprintf(stderr, "real_user_var_add: name='%s', value=%p, u=NULL\n",
                name, value);
	if (pu != NULL) {
	    *pu = NULL;
	}
        return err ? err : E_DATA;
    }

    /* We use OPT_P for a private variable, OPT_A
       when adding as a function argument, OPT_S
       when adding as a "shell" variable, OPT_C
       when we're auto-casting a 1 x 1 matrix result
       to a scalar.
    */

#if UVDEBUG
    fprintf(stderr, "real_user_var_add: '%s', level %d, err = %d\n",
            name, u->level, err);
#endif

    if (!err) {
        err = resize_uvar_stack(n_vars + 1);
        if (!err) {
            if (opt & OPT_P) {
                u->flags = UV_PRIVATE;
            } else if (opt & OPT_S) {
                u->flags = UV_SHELL;
            }
            if (opt & OPT_A) {
                u->flags &= ~UV_MAIN;
                u->level += 1;
            }
            if (opt & OPT_C) {
                u->flags |= UV_NODECL;
            }
            uvars[n_vars] = u;
            set_nvars(n_vars + 1, "user_var_add");
        }
    }

    if (!err && user_var_callback != NULL && u->level == 0 &&
        !(opt & (OPT_P | OPT_S)) && *name != '$' &&
        (type == GRETL_TYPE_MATRIX || type == GRETL_TYPE_BUNDLE) &&
        !(type == GRETL_TYPE_BUNDLE && bname_is_temp(name))) {
        return (*user_var_callback)(name, type, UVAR_ADD);
    }

    if (pu != NULL) {
	*pu = u;
    }

    return err;
}

/**
 * user_var_add:
 * @name: name to give the variable.
 * @type: the type of the variable.
 * @value: pointer to value for variable.
 *
 * Adds a new user-variable with the given characteristics.
 * Note that the user-variable takes ownership of the
 * supplied @value; this should be copied first if need be.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int user_var_add (const char *name, GretlType type, void *value)
{
    return real_user_var_add(name, type, value, OPT_NONE, NULL);
}

/**
 * alt_user_var_add:
 * @name: name to give the variable.
 * @type: the type of the variable.
 * @value: pointer to value for variable.
 *
 * Adds a new user_var with the given characteristics.
 * Note that the user_var takes ownership of the supplied
 * @value; this should be copied first if need be.
 *
 * Returns: pointer to the new user_var, or NULL on failure.
 */

user_var *alt_user_var_add (const char *name, GretlType type, void *value)
{
    user_var *uv = NULL;

    real_user_var_add(name, type, value, OPT_NONE, &uv);
    return uv;
}

int private_matrix_add (gretl_matrix *M, const char *name)
{
    return real_user_var_add(name, GRETL_TYPE_MATRIX,
			     M, OPT_P, NULL);
}

int private_scalar_add (double val, const char *name)
{
    double *px = malloc(sizeof *px);
    int err;

    if (px == NULL) {
        err = E_ALLOC;
    } else {
        *px = val;
        err = real_user_var_add(name, GRETL_TYPE_DOUBLE,
                                px, OPT_P, NULL);
    }

    return err;
}

/**
 * user_var_delete_by_name:
 * @name: name of the variable to delete.
 * @prn: pointer to gretl printer, or NULL.
 *
 * Deletes the specified user-variable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int user_var_delete_by_name (const char *name, PRN *prn)
{
    GretlType type = 0;
    int level = gretl_function_depth();
    user_var *targ = NULL;
    int i, j, k = 0;
    int err = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->level == level && !strcmp(uvars[i]->name, name)) {
            targ = uvars[i];
            k = i;
            break;
        }
    }

    if (targ == NULL) {
        return E_UNKVAR;
    }

    if (level > 0 && (targ->flags & UV_MAIN)) {
        gretl_errmsg_sprintf(_("%s: cannot be deleted here"), targ->name);
        return E_DATA;
    }

    if (user_var_callback != NULL && level == 0 &&
        !var_is_private(targ) &&
        (targ->type == GRETL_TYPE_MATRIX ||
         targ->type == GRETL_TYPE_BUNDLE)) {
        /* run this deletion through the GUI program to ensure
           that things stay in sync
        */
        return (*user_var_callback)(name, targ->type,
                                    UVAR_DELETE);
    }

    type = targ->type;
    user_var_destroy(targ);
    for (j=k; j<n_vars-1; j++) {
        uvars[j] = uvars[j+1];
    }
    resize_uvar_stack(n_vars - 1);
    set_nvars(n_vars - 1, "user_var_delete_by_name");

    if (prn != NULL && gretl_messages_on()) {
        pprintf(prn, _("Deleted %s"), name);
        pputc(prn, '\n');
    }
    if (level == 0 && type == GRETL_TYPE_DOUBLE &&
        scalar_edit_callback != NULL) {
        scalar_edit_callback();
    }

    return err;
}

int user_var_delete (user_var *uvar)
{
    int i, j, err = E_UNKVAR;

    for (i=0; i<n_vars; i++) {
        if (uvar == uvars[i]) {
            user_var_destroy(uvars[i]);
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            set_nvars(n_vars - 1, "user_var_delete");
            err = 0;
            break;
        }
    }

    return err;
}

#if HDEBUG > 1

static int uvar_index (user_var *u)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (u == uvars[i]) {
            return i;
        }
    }

    return -1;
}

#endif

/* Try to guess whether the currently-called function is big enough
   (number of lines of code) to make it worthwhile to construct a hash
   table for uservars at its level of execution, namely @uvh1, given
   that we'll have to empty the table on exit from the function.

   The number-of-lines threshold here is obviously kinda arbitrary;
   some systematic experimentation might be useful.
*/

static inline int use_uvh1 (void)
{
    /* maybe add: && !gretl_function_recursing() */
    return current_function_size() > 40;
}

user_var *get_user_var_of_type_by_name (const char *name,
                                        GretlType type)
{
    int prev_d = get_previous_depth();
    int d = gretl_function_depth();
    int i, imin = 0;
    user_var *u = NULL;

    if (name == NULL || *name == '\0') {
        return NULL;
    }

    if (type == GRETL_TYPE_DOUBLE) {
        /* support "auxiliary scalars" mechanism */
        imin = scalar_imin;
    }

#if HDEBUG > 1
    int hfound = 0;

    fprintf(stderr, "get user var: '%s', %s (n_vars=%d, level=%d, "
            "previous=%d, imin=%d)\n", name, gretl_type_get_name(type),
            n_vars, d, prev_d, imin);
# if HDEBUG > 2
    fputs("uvars list:\n", stderr);
    for (i=0; i<n_vars; i++) {
        fprintf(stderr, " %d: '%s', %s, level %d, ptr %p\n", i,
                uvars[i]->name, gretl_type_get_name(uvars[i]->type),
                uvars[i]->level, uvars[i]->ptr);
    }
# endif
#endif

    if (d != prev_d) {
        if (d == 0) {
            /* we're now at "main" level */
            if (uvh0 == NULL) {
                uvh0 = g_hash_table_new(g_str_hash, g_str_equal);
#if HDEBUG
                fprintf(stderr, "uvh0: d=0, allocated at %p\n", uvh0);
#endif
            }
            if (uvh1 != NULL) {
#if HDEBUG
                fprintf(stderr, "d=0, prev=%d: clear uvh1 at %p\n",
                        prev_d, uvh1);
#endif
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = uvh0;
        } else if (!use_uvh1()) {
            /* exec'ing a function, hash table not wanted */
            if (prev_d > 0 && uvh1 != NULL) {
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = NULL;
        } else {
            /* exec'ing a function, hash table wanted */
            if (uvh1 == NULL) {
                uvh1 = g_hash_table_new(g_str_hash, g_str_equal);
#if HDEBUG
                fprintf(stderr, "uvh1: d=%d, prev=%d, allocated at %p\n",
                        d, prev_d, uvh1);
#endif
            } else if (prev_d > 0 && uvh1 != NULL) {
#if HDEBUG
                fprintf(stderr, "d=%d, prev=%d: clear uvh1 at %p\n",
                        d, prev_d, uvh1);
#endif
                g_hash_table_remove_all(uvh1);
            }
            uvars_hash = uvh1;
        }
        set_previous_depth(d);
    }

    if (uvars_hash != NULL) {
        /* first resort: try a hash look-up */
        u = g_hash_table_lookup(uvars_hash, name);
        /* but verify type, if specified */
        if (u != NULL && type != GRETL_TYPE_ANY && u->type != type) {
            u = NULL;
        }
#if HDEBUG > 1
        if (u != NULL) hfound = 1;
#endif
    }

    if (u == NULL) {
        /* "On demand" hashing: if we're successful in looking
           up a variable in the traditional manner, then
           insert it into the uservars hash table.
        */
        for (i=imin; i<n_vars; i++) {
            if (uvars[i]->level == d &&
                (type == GRETL_TYPE_ANY || uvars[i]->type == type) &&
                !strcmp(uvars[i]->name, name)) {
                u = uvars[i];
                if (uvars_hash != NULL) {
                    g_hash_table_insert(uvars_hash, u->name, u);
                }
                break;
            }
        }
    }

#if HDEBUG > 1
    if (hfound)
        fprintf(stderr, "found at pos %d via hash (%s)\n\n", uvar_index(u),
                uvars_hash == uvh1 ? "uvh1" : "uvh0");
    else if (u != NULL)
        fprintf(stderr, "found at pos %d via regular search\n\n", uvar_index(u));
    else
        fprintf(stderr, "not found\n\n");
#endif

    return u;
}

user_var *get_user_var_by_name (const char *name)
{
    return get_user_var_of_type_by_name(name, GRETL_TYPE_ANY);
}

GretlType user_var_get_type_by_name (const char *name)
{
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_ANY);

    return u == NULL ? GRETL_TYPE_NONE : u->type;
}

GretlType user_var_get_specific_type (const char *name)
{
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_ANY);

    if (u != NULL && u->type == GRETL_TYPE_ARRAY) {
	return gretl_array_get_type(u->ptr);
    } else {
	return u == NULL ? GRETL_TYPE_NONE : u->type;
    }
}

void *user_var_get_value_and_type (const char *name,
                                   GretlType *type)
{
    void *ret = NULL;
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_ANY);

    if (u != NULL) {
        ret = u->ptr;
        *type = u->type;
    } else {
        *type = GRETL_TYPE_NONE;
    }

    return ret;
}

const char *uservar_name_complete (const char *s)
{
    const char *ret = NULL;

    if (uvars_hash != NULL) {
        GList *hk = g_hash_table_get_keys(uvars_hash);
        int n = strlen(s);

        while (hk != NULL) {
            if (!strncmp((const char *) hk->data, s, n)) {
                ret = (const char *) hk->data;
                break;
            }
            hk = hk->next;
        }
        g_list_free(hk);
    }

    return ret;
}

int gretl_is_user_var (const char *name)
{
    return get_user_var_by_name(name) != NULL;
}

user_var *get_user_var_by_data (const void *data)
{
    int i, d = gretl_function_depth();

    if (data == NULL) {
        return NULL;
    }

    for (i=0; i<n_vars; i++) {
        if (uvars[i] != NULL && uvars[i]->level == d &&
            uvars[i]->ptr == data) {
            return uvars[i];
        }
    }

    return NULL;
}

const char *user_var_get_name (user_var *uvar)
{
    return uvar == NULL ? NULL : uvar->name;
}

const char *user_var_get_name_by_data (const void *data)
{
    user_var *u = get_user_var_by_data(data);

    return u == NULL ? NULL : u->name;
}

int user_var_get_level (user_var *uvar)
{
    return (uvar == NULL)? -1 : uvar->level;
}

int user_var_get_flags (user_var *uvar)
{
    return (uvar == NULL)? 0 : (int) uvar->flags;
}

int user_var_set_flag (user_var *uvar, UVFlags flag)
{
    if (uvar != NULL) {
        uvar->flags |= flag;
        return 0;
    } else {
        return E_INVARG;
    }
}

int user_var_unset_flag (user_var *uvar, UVFlags flag)
{
    if (uvar != NULL) {
        uvar->flags &= ~flag;
        return 0;
    } else {
        return E_INVARG;
    }
}

void user_var_privatize_by_name (const char *name)
{
    user_var *u = get_user_var_by_name(name);

    if (u != NULL) {
        u->flags |= UV_PRIVATE;
    }
}

void *user_var_get_value (user_var *uvar)
{
    return (uvar == NULL)? NULL : uvar->ptr;
}

GretlType user_var_get_type (user_var *uvar)
{
    return (uvar == NULL)? 0 : uvar->type;
}

void *user_var_get_value_by_name (const char *name)
{
    user_var *u = get_user_var_by_name(name);

    return (u == NULL)? NULL : u->ptr;
}

/* special for scalars since user_var_get_value returns
   a pointer */

double user_var_get_scalar_value (user_var *uvar)
{
    if (uvar != NULL && uvar->type == GRETL_TYPE_DOUBLE) {
        return *(double *) uvar->ptr;
    } else {
        return NADBL;
    }
}

int user_var_set_scalar_value (user_var *uvar, double x)
{
    if (uvar != NULL && uvar->type == GRETL_TYPE_DOUBLE) {
        *(double *) uvar->ptr = x;
        return 0;
    } else {
        return E_DATA;
    }
}

int user_var_adjust_level (user_var *uvar, int adj)
{
    if (uvar == NULL) {
        return E_UNKVAR;
    } else {
        uvar->level += adj;
        return 0;
    }
}

/* Note: the following should be called only from internal
   contexts in which we know that the attempted renaming
   is not broken (e.g. trying to assign to @uvar a name
   that is already taken by some other object).
*/

int user_var_set_name (user_var *uvar, const char *name)
{
    int err = 0;

    if (uvar == NULL) {
        err = E_DATA;
    } else {
        *uvar->name = '\0';
        strncat(uvar->name, name, VNAMELEN - 1);
    }

    return err;
}

/**
 * user_var_localize:
 * @origname: name of variable at caller level.
 * @localname: name to be used within function.
 * @type: type of the variable.
 *
 * On entry to a function, renames the named variable (provided
 * as an argument) and sets its level so that is is accessible
 * within the function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int user_var_localize (const char *origname,
                       const char *localname,
                       GretlType type)
{
    user_var *u;
    int err = 0;

    if (gretl_is_array_ref_type(type)) {
        type = GRETL_TYPE_ARRAY;
    } else {
        type = gretl_type_get_plain_type(type);
    }

    u = get_user_var_of_type_by_name(origname, type);

    if (u == NULL) {
        err = E_DATA;
    } else {
        user_var_set_name(u, localname);
        u->level += 1;
    }

    return err;
}

static int user_var_count_for_type (GretlType type)
{
    int i, n = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == type) {
            n++;
        }
    }

    return n;
}

int n_user_matrices (void)
{
    return user_var_count_for_type(GRETL_TYPE_MATRIX);
}

int n_user_scalars (void)
{
    return user_var_count_for_type(GRETL_TYPE_DOUBLE);
}

int n_user_lists (void)
{
    return user_var_count_for_type(GRETL_TYPE_LIST);
}

int n_user_bundles (void)
{
    return user_var_count_for_type(GRETL_TYPE_BUNDLE);
}

/**
 * user_var_replace_value:
 * @uvar: user variable.
 * @value: the new value to place as the value or @uvar.
 * @type: the type of the replacement value.
 *
 * Replaces the value of @uvar; the existing value is
 * freed first.
 *
 * Returns: 0 on success, non-zero on error.
 */

int user_var_replace_value (user_var *uvar, void *value,
                            GretlType type)
{
    int err = 0;

    if (uvar == NULL) {
        err = E_UNKVAR;
    } else if (value != uvar->ptr && (uvar->flags & UV_NOREPL)) {
        gretl_errmsg_sprintf(_("The variable %s is read-only"), uvar->name);
        err = E_DATA;
    } else if (type != uvar->type) {
        err = E_TYPES; /* assume the worst */
        if (uvar->type == GRETL_TYPE_ARRAY && uvar->ptr != NULL) {
            /* but we might be OK */
            if (type == gretl_array_get_type(uvar->ptr)) {
                err = 0;
            }
        }
        if (err) {
            fputs("*** user_var_replace_value: type mismatch ***\n", stderr);
            fprintf(stderr, " (expected %s but got %s)\n",
                    gretl_type_get_name(uvar->type), gretl_type_get_name(type));
        }
    }

    if (!err && value != uvar->ptr) {
        if (uvar->ptr != NULL) {
            uvar_free_value(uvar);
        }
        uvar->ptr = value;
    }

    return err;
}

/**
 * user_var_set_pointer:
 * @uvar: user variable.
 * @ptr: the new pointer to set on @uvar.
 *
 * Unlike @user_var_replace_value, this function does not free
 * the existing data pointer on @uvar. It should therefore be
 * used only when another pointer to the prior data exists,
 * otherwise a memory leak would result. Plus there's no type
 * checking so use with care.
 *
 * Returns: 0 on success, error code if @uvar is NULL.
 */

int user_var_set_pointer (user_var *uvar, void *ptr)
{
    if (uvar != NULL) {
	uvar->ptr = ptr;
	return 0;
    } else {
	return E_DATA;
    }
}

char *user_string_resize (const char *name, size_t len, int *err)
{
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_STRING);

    if (u == NULL) {
        *err = E_INVARG;
        return NULL;
    } else {
        char *orig = u->ptr;

        if (orig == NULL || len > strlen(orig) + 1) {
            char *tmp = realloc(u->ptr, len);

            if (tmp == NULL) {
                *err = E_ALLOC;
            } else {
                u->ptr = tmp;
            }
        }
    }

    return (char *) u->ptr;
}

char *user_string_reset (const char *name, const char *repl, int *err)
{
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_STRING);

    if (u == NULL) {
        *err = E_INVARG;
        return NULL;
    } else {
        free(u->ptr);
        if (repl == NULL) {
            u->ptr = gretl_strdup("");
        } else {
            u->ptr = gretl_strdup(repl);
        }
        return (char *) u->ptr;
    }
}

static int check_array_type_compat (GretlType type,
                                    user_var *u)
{
    int err = 0;

    if (u->type != GRETL_TYPE_ARRAY) {
        err = E_TYPES;
    } else {
        /* we also need a more specific check here */
        if (type != gretl_array_get_type(u->ptr)) {
            err = E_TYPES;
        }
    }

    return err;
}

int user_var_add_or_replace (const char *name,
                             GretlType type,
                             void *value)
{
    user_var *u = get_user_var_by_name(name);
    int err = 0;

    if (u != NULL) {
        if (gretl_array_type(type)) {
            err = check_array_type_compat(type, u);
        } else if (u->type != type) {
            err = E_TYPES;
        }
        if (!err) {
            err = user_var_replace_value(u, value, type);
        }
    } else {
        err = real_user_var_add(name, type, value, OPT_NONE, NULL);
    }

    return err;
}

void *user_var_steal_value (user_var *uvar)
{
    void *ret = NULL;

    if (uvar != NULL) {
        ret = uvar->ptr;
        uvar->ptr = NULL;
    }

    return ret;
}

/* FIXME: are both the above and the below necessary? */

void *user_var_unstack_value (user_var *uvar)
{
    void *ret = NULL;
    int i, j;

    for (i=0; i<n_vars; i++) {
        if (uvar == uvars[i]) {
            ret = uvar->ptr;
            uvars[i]->ptr = NULL;
            user_var_destroy(uvars[i]);
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            set_nvars(n_vars - 1, "user_var_unstack_value");
            break;
        }
    }

    return ret;
}

int user_matrix_replace_matrix_by_name (const char *name,
                                        gretl_matrix *m)
{
    user_var *u = get_user_var_by_name(name);

    if (u != NULL) {
        return user_var_replace_value(u, m, GRETL_TYPE_MATRIX);
    } else {
        return E_DATA;
    }
}

GList *user_var_names_for_type (GretlType type)
{
    GList *list = NULL;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == type) {
            list = g_list_append(list, (gpointer) uvars[i]->name);
        }
    }

    return list;
}

GList *user_var_list_for_type (GretlType type)
{
    GList *list = NULL;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == type) {
            list = g_list_append(list, (gpointer) uvars[i]);
        }
    }

    return list;
}

/**
 * set_user_var_callback:
 * @callback: function function to put in place.
 *
 * Sets the callback function to be invoked when a user-defined
 * matrix is added to or removed from the stack of saved objects.
 * Intended for synchronizing the GUI program with the saved object
 * state.
 */

void set_user_var_callback (USER_VAR_FUNC callback)
{
    user_var_callback = callback;
}

void set_scalar_edit_callback (void (*callback))
{
    scalar_edit_callback = callback;
}

/**
 * arg_add_as_shell:
 * @name: the name to be given to the "shell" variable.
 * @type: the type of the variable.
 * @value: the value pointer
 *
 * The value in question is added to the stack of named
 * variables under the name @name with the shell flag
 * set. This is used (a) when an anonymous matrix is given as
 * a %const argument to a user-defined function and (b) when
 * an anonymous bundle is given as the argument corresponding
 * to a bundle-pointer parameter. The @value becomes
 * accessible by @name within the function, but is protected
 * from destruction on exit from the function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int arg_add_as_shell (const char *name, GretlType type,
                      void *value)
{
    return real_user_var_add(name, type, value, OPT_S | OPT_A, NULL);
}

/**
 * copy_matrix_as:
 * @m: the original matrix.
 * @newname: the name to be given to the copy.
 * @fnarg: 0 for regular use.
 *
 * A copy of matrix @m is added to the stack of saved matrices
 * under the name @newname.
 *
 * The @fnarg argument should be non-zero only if this function
 * is used to handle the case where a matrix is given as the argument
 * to a user-defined function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int copy_matrix_as (const gretl_matrix *m, const char *newname,
                    int fnarg)
{
    gretl_matrix *m2 = gretl_matrix_copy(m);
    int err = 0;

    if (m2 == NULL) {
        err = E_ALLOC;
    } else {
        gretlopt opt = fnarg ? OPT_A : OPT_NONE;

        err = real_user_var_add(newname, GRETL_TYPE_MATRIX,
				m2, opt, NULL);
    }

    return err;
}

int copy_as_arg (const char *param_name, GretlType type, void *value)
{
    void *copyval = NULL;
    GretlType cpytype = type;
    int err = 0;

    if (type == GRETL_TYPE_MATRIX) {
        gretl_matrix *mcpy = gretl_matrix_copy((gretl_matrix *) value);

        if (mcpy == NULL) {
            err = E_ALLOC;
        } else {
            copyval = mcpy;
        }
    } else if (type == GRETL_TYPE_LIST) {
        int *lcpy = gretl_list_copy((int *) value);

        if (lcpy == NULL) {
            err = E_ALLOC;
        } else {
            copyval = lcpy;
        }
    } else if (type == GRETL_TYPE_STRING) {
        char *scpy = gretl_strdup((char *) value);

        if (scpy == NULL) {
            err = E_ALLOC;
        } else {
            copyval = scpy;
        }
    } else if (type == GRETL_TYPE_DOUBLE) {
        double *px = malloc(sizeof *px);

        if (px == NULL) {
            err = E_ALLOC;
        } else {
            *px = *(double *) value;
            copyval = px;
        }
    } else if (type == GRETL_TYPE_BUNDLE) {
        gretl_bundle *bcpy = gretl_bundle_copy(value, &err);

        if (!err) {
            copyval = bcpy;
        }
    } else if (gretl_array_type(type)) {
        gretl_array *acpy = gretl_array_copy(value, &err);

        if (!err) {
            copyval = acpy;
            cpytype = gretl_array_get_type(acpy);
        }
    }

    if (!err) {
        err = real_user_var_add(param_name, cpytype,
				copyval, OPT_A, NULL);
    }

    return err;
}

int *copy_list_as_arg (const char *param_name, int *list,
                       int *err)
{
    int *ret = NULL;

    *err = copy_as_arg(param_name, GRETL_TYPE_LIST, list);
    if (!*err) {
        ret = uvars[n_vars-1]->ptr;
    }

    return ret;
}

void destroy_user_vars (void)
{
    int i, j;

#if HDEBUG
    fprintf(stderr, "destroy_user_vars, uvars_hash = %p (uvh0 %p, uvh1 %p)\n",
            (void *) uvars_hash, (void *) uvh0, (void *) uvh1);
#endif

    for (i=0; i<n_vars; i++) {
        if (uvars[i] == NULL) {
            break;
        }
        user_var_destroy(uvars[i]);
        for (j=i; j<n_vars-1; j++) {
            uvars[j] = uvars[j+1];
        }
        uvars[n_vars-1] = NULL;
        i--;
    }

    if (uvh0 != NULL || uvh1 != NULL) {
        uvar_hash_destroy();
    }

    set_nvars(0, "destroy_user_vars");

    free(uvars);
    uvars = NULL;
    n_alloc = 0;
}

static int uvar_levels_match (user_var *u, int level)
{
    int ret = 0;

    if (u->level == level) {
        ret = 1;
    } else if (level == LEV_PRIVATE && var_is_private(u)) {
        ret = 1;
    }

    return ret;
}

static int real_destroy_user_vars_at_level (int level, int type,
                                            int imin)
{
    int i, j, nv = imin;
    int err = 0;

#if HDEBUG
    fprintf(stderr, "real_destroy_user_vars_at_level: level %d, "
            "type %d (%s), imin=%d\n", level, type,
            gretl_type_get_name(type), imin);
#endif

    for (i=imin; i<n_vars; i++) {
        if (uvars[i] == NULL) {
            break;
        }
        if (type > 0 && uvars[i]->type != type) {
            /* preserve this variable */
            nv++;
            continue;
        }
        if (uvar_levels_match(uvars[i], level)) {
            user_var_destroy(uvars[i]);
            /* shuffle the remainder down one place */
            for (j=i; j<n_vars-1; j++) {
                uvars[j] = uvars[j+1];
            }
            uvars[n_vars-1] = NULL;
            i--;
        } else {
            /* preserving */
            nv++;
        }
    }

    set_nvars(nv, "real_destroy_user_vars_at_level");

    return err;
}

static int destroy_user_vars_via_callback (int type)
{
    user_var **delvars = NULL;
    int i, j, n = 0;
    int err = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->level == 0 && uvars[i]->type == type) {
            n++;
        }
    }

    if (n == 0) {
        return 0;
    }

    delvars = malloc(n * sizeof *delvars);
    if (delvars == NULL) {
        return E_ALLOC;
    }

    j = 0;
    for (i=0; i<n_vars; i++) {
        if (uvars[i]->level == 0 && uvars[i]->type == type) {
            delvars[j++] = uvars[i];
        }
    }

    for (j=0; j<n && !err; j++) {
        err = (*user_var_callback)(delvars[j]->name,
                                   delvars[j]->type,
                                   UVAR_DELETE);
    }

    free(delvars);

    return err;
}

/**
 * destroy_user_vars_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of user matrices all
 * matrices that were created at the given @level.  This is
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_vars_at_level (int level)
{
    return real_destroy_user_vars_at_level(level, 0, 0);
}

int destroy_private_uvars (void)
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE, 0, 0);
}

int destroy_private_matrices (void)
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE,
                                           GRETL_TYPE_MATRIX,
                                           0);
}

int destroy_private_lists (void)
{
    return real_destroy_user_vars_at_level(LEV_PRIVATE,
                                           GRETL_TYPE_LIST,
                                           0);
}

int delete_user_vars_of_type (GretlType type, PRN *prn)
{
    int err = 0;

    if (type == GRETL_TYPE_MATRIX ||
        type == GRETL_TYPE_BUNDLE ||
        type == GRETL_TYPE_ARRAY  ||
        type == GRETL_TYPE_STRING ||
        type == GRETL_TYPE_DOUBLE ||
        type == GRETL_TYPE_LIST) {
        int level = gretl_function_depth();

        if (level == 0 && user_var_callback != NULL &&
            (type == GRETL_TYPE_MATRIX || type == GRETL_TYPE_BUNDLE)) {
            err = destroy_user_vars_via_callback(type);
        } else {
            err = real_destroy_user_vars_at_level(level, type, 0);
        }

        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Deleted all variables of type %s\n"),
                    gretl_type_get_name(type));
        }
    } else {
        err = E_TYPES;
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
    real_destroy_user_vars_at_level(LEV_PRIVATE,
                                    GRETL_TYPE_DOUBLE,
                                    0);
}

char *temp_name_for_bundle (void)
{
    char tmpname[VNAMELEN];
    int i, nb = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_BUNDLE) {
            nb++;
        }
    }

    sprintf(tmpname, "btmp___%d", nb);
    return gretl_strdup(tmpname);
}

static void xml_put_user_matrix (user_var *u, PRN *prn)
{
    if (u != NULL && u->ptr != NULL) {
        gretl_matrix_serialize(u->ptr, u->name, prn);
    }
}

static void write_scalar_value (double x, const char *fmt, PRN *prn)
{
    if (na(x)) {
#ifdef WIN32
        win32_pprint_nonfinite(prn, x, '\n');
#else
        pprintf(prn, "%g\n", x);
#endif
    } else {
        pprintf(prn, fmt, x);
    }
}

static void serialize_scalar_value (double x, PRN *prn)
{
    if (na(x)) {
#ifdef WIN32
        win32_pprint_nonfinite(prn, x, 0);
#else
        pprintf(prn, "%g", x);
#endif
    } else {
        pprintf(prn, "%.16g", x);
    }
}

/**
 * print_scalars:
 * @prn: pointer to gretl printing struct.
 *
 * Prints names and values of any saved scalars.
 */

void print_scalars (PRN *prn)
{
    double x;
    int level = gretl_function_depth();
    int len, ns = 0, maxlen = 0;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_DOUBLE &&
            uvars[i]->level == level) {
            len = strlen(uvars[i]->name);
            if (len > maxlen) {
                maxlen = len;
            }
            ns++;
        }
    }

    if (ns == 0) {
        pprintf(prn, "%s\n", _("none"));
        return;
    }

    pputc(prn, '\n');

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_DOUBLE &&
            uvars[i]->level == level) {
            x = *(double *) uvars[i]->ptr;
            pprintf(prn, " %*s = ", maxlen, uvars[i]->name);
            write_scalar_value(x, "%.16g\n", prn);
        }
    }

    pputc(prn, '\n');
}

void print_scalar_by_name (const char *name, PRN *prn)
{
    user_var *u;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_DOUBLE);

    if (u != NULL) {
        double x = *(double *) u->ptr;

        pprintf(prn, "\n%15s = ", u->name);
        write_scalar_value(x, "% #.8g\n", prn);
    }
}

/* "auxiliary scalars": this apparatus is used when we want to do
   "private" NLS estimation (e.g. in ARMA initialization).  It ensures
   that the scalar NLS parameters don't collide with the public scalar
   namespace. FIXME.
*/

void set_auxiliary_scalars (void)
{
    scalar_imin = n_vars;
}

void unset_auxiliary_scalars (void)
{
    real_destroy_user_vars_at_level(0, GRETL_TYPE_DOUBLE, scalar_imin);
    scalar_imin = 0;
}

static int real_scalar_add (const char *name, double val,
                            gretlopt opt)
{
    user_var *u = get_user_var_by_name(name);
    int level = gretl_function_depth();
    int err = 0;

    if (u != NULL) {
        if (u->type == GRETL_TYPE_DOUBLE) {
            *(double *) u->ptr = val;
        } else {
            err = E_TYPES;
        }
        return err;
    } else {
        double *px = malloc(sizeof *px);

        if (px == NULL) {
            err = E_ALLOC;
        } else {
            *px = val;
            err = real_user_var_add(name, GRETL_TYPE_DOUBLE,
                                    px, opt, NULL);
        }

        if (!err && level == 0 && scalar_edit_callback != NULL) {
            scalar_edit_callback();
        }
    }

    return err;
}

int gretl_scalar_add (const char *name, double val)
{
    return real_scalar_add(name, val, OPT_NONE);
}

int gretl_scalar_add_mutable (const char *name, double val)
{
    return real_scalar_add(name, val, OPT_C);
}

int gretl_scalar_convert_to_matrix (user_var *u)
{
    gretl_matrix *m = NULL;

    if (u == NULL) {
        return E_UNKVAR;
    } else if (u->type != GRETL_TYPE_DOUBLE) {
        return E_TYPES;
    }

    m = gretl_matrix_alloc(1, 1);
    if (m == NULL) {
        return E_ALLOC;
    }

    m->val[0] = *(double *) u->ptr;
    free(u->ptr);
    u->ptr = m;
    u->type = GRETL_TYPE_MATRIX;

    if (gretl_function_depth() == 0) {
        if (scalar_edit_callback != NULL) {
            (*scalar_edit_callback)();
        }
        if (user_var_callback != NULL) {
            (*user_var_callback)(u->name, GRETL_TYPE_MATRIX, UVAR_ADD);
        }
    }

    return 0;
}

int add_auxiliary_scalar (const char *name, double val)
{
    double *px = malloc(sizeof *px);
    int err;

    /* Note that unlike gretl_scalar_add() above, this function
       adds a new scalar variable unconditionally; it never
       modifies the value of an existing scalar of the same
       name.
    */

    if (px == NULL) {
        err = E_ALLOC;
    } else {
        *px = val;
        err = real_user_var_add(name, GRETL_TYPE_DOUBLE,
                                px, OPT_NONE, NULL);
    }

    return err;
}

int gretl_scalar_set_value (const char *name, double val)
{
    user_var *u;
    int err = 0;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_DOUBLE);

    if (u == NULL) {
        gretl_errmsg_sprintf(_("%s: no such scalar"), name);
        err = E_DATA;
    } else if (scalar_is_read_only_index(name)) {
        err = E_DATA;
        gretl_errmsg_sprintf(_("The variable %s is currently read-only"), name);
    } else {
        *(double *) u->ptr = val;

        if (scalar_edit_callback != NULL) {
            scalar_edit_callback();
        }
    }

    return err;
}

/* get the value from a user variable of scalar type */

double gretl_scalar_get_value (const char *name, int *err)
{
    user_var *u;
    double ret = NADBL;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_DOUBLE);

    if (u != NULL) {
        ret = *(double *) u->ptr;
    } else {
        ret = get_const_by_name(name, err);
    }

    return ret;
}

static double maybe_get_bundled_scalar (const char *name, int *err)
{
    const char *p = strchr(name, '.');
    gretl_bundle *b = NULL;
    char bname[VNAMELEN];
    char key[VNAMELEN];
    double x = NADBL;

    *bname = '\0';
    strncat(bname, name, p - name);
    b = get_bundle_by_name(bname);

    if (b == NULL) {
        *err = E_INVARG;
    } else {
        *key = '\0';
        strncat(key, p + 1, VNAMELEN - 1);
        x = gretl_bundle_get_scalar(b, key, err);
    }

    return x;
}

/* more "permissive" than gretl_scalar_get_value(): allows
   for @name being the identifier for a 1 x 1 matrix, or
   bundle.member
*/

double get_scalar_value_by_name (const char *name, int *err)
{
    double ret = NADBL;
    user_var *u;

    if (strchr(name, '.')) {
        ret = maybe_get_bundled_scalar(name, err);
        goto bailout;
    }

    u = get_user_var_by_name(name);

    if (u != NULL) {
        if (u->type == GRETL_TYPE_DOUBLE) {
            ret = *(double *) u->ptr;
        } else if (u->type == GRETL_TYPE_MATRIX) {
            gretl_matrix *m = u->ptr;

            if (gretl_matrix_is_scalar(m)) {
                ret = m->val[0];
            } else {
                *err = E_TYPES;
            }
        } else {
            *err = E_TYPES;
        }
    } else {
        ret = get_const_by_name(name, err);
    }

 bailout:

    if (*err) {
        gretl_errmsg_sprintf(_("'%s': not a scalar"), name);
    }

    return ret;
}

int gretl_is_scalar (const char *name)
{
    int ret = 0;

    if (get_user_var_of_type_by_name(name, GRETL_TYPE_DOUBLE) != NULL) {
        ret = 1;
    }

    if (!ret) {
        ret = const_lookup(name);
    }

    return ret;
}

/**
 * get_string_by_name:
 * @name: the name of the string variable to access.
 *
 * Returns: the value of string variable @name, or %NULL
 * if there is no such variable. Note that this is the
 * actual string value, not a copy thereof, compare
 * copy_string_by_name().
 */

char *get_string_by_name (const char *name)
{
    user_var *u = NULL;

    if (name == NULL) {
        return NULL;
    }

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_STRING);
    if (u != NULL) {
        return (char *) u->ptr;
    } else {
        return get_built_in_string_by_name(name);
    }
}

/**
 * copy_string_by_name:
 * @name: the name of the string variable to access.
 * @err: location to receive error code.
 *
 * Returns: a copy of the value of string variable @name,
 * or %NULL on failure.
 */

char *copy_string_by_name (const char *name, int *err)
{
    user_var *u;
    const char *s;
    char *ret = NULL;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_STRING);

    if (u != NULL) {
        s = u->ptr;
    } else {
        s = get_built_in_string_by_name(name);
    }

    if (s == NULL) {
        *err = E_DATA;
    } else {
        ret = gretl_strdup(s);
        if (ret == NULL) {
            *err = E_ALLOC;
        }
    }

    return ret;
}

/**
 * gretl_is_string:
 * @name: name to test.
 *
 * Returns: 1 if @name is the name of a currently defined
 * string variable, otherwise 0.
 */

int gretl_is_string (const char *name)
{
    if (*name == '@' && *(name + 1) != '@') {
        name++;
    }

    if (get_user_var_of_type_by_name(name, GRETL_TYPE_STRING) != NULL) {
        return 1;
    } else if (get_built_in_string_by_name(name) != NULL) {
        return 1;
    } else {
        return 0;
    }
}

int is_user_string (const char *name)
{
    if (*name == '@' && *(name + 1) != '@') {
        name++;
    }

    if (get_user_var_of_type_by_name(name, GRETL_TYPE_STRING) != NULL) {
        return 1;
    } else {
        return 0;
    }
}

int max_varno_in_saved_lists (void)
{
    int *list;
    int i, j, vmax = 0;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_LIST) {
            list = uvars[i]->ptr;
            if (list != NULL) {
                for (j=1; j<=list[0]; j++) {
                    if (list[j] > vmax) {
                        vmax = list[j];
                    }
                }
            }
        }
    }

    return vmax;
}

static int var_is_deleted (const int *dlist, int dmin, int i)
{
    int v = dmin + i - 1;

    if (dlist != NULL) {
        return in_gretl_list(dlist, v);
    } else {
        return (v >= dmin);
    }
}

/**
 * gretl_lists_revise:
 * @dlist: list of variables to be deleted (or NULL).
 * @dmin: lowest ID number of deleted var (referenced only
 * if @dlist is NULL).
 *
 * Goes through any saved lists, adjusting the ID numbers
 * they contain to reflect the deletion from the dataset of
 * certain variables: those referenced in @dlist, if given,
 * or if @dlist is NULL, those variables with IDs greater
 * than or equal to @dmin.
 *
 * Returns: 0 on success, non-zero code on failure.
 */

int gretl_lists_revise (const int *dlist, int dmin)
{
    int *list, *maplist;
    int lmax = 0;
    int i, j, k;

    if (dlist != NULL) {
        /* determine lowest deleted ID */
        dmin = dlist[1];
        for (i=2; i<=dlist[0]; i++) {
            if (dlist[i] > 0 && dlist[i] < dmin) {
                dmin = dlist[i];
            }
        }
    }

    /* find highest ID ref'd in any saved list */
    for (j=0; j<n_vars; j++) {
        if (uvars[j]->type == GRETL_TYPE_LIST) {
            list = uvars[j]->ptr;
            if (list != NULL) {
                for (i=1; i<=list[0]; i++) {
                    if (list[i] > lmax) {
                        lmax = list[i];
                    }
                }
            }
        }
    }

    if (lmax < dmin) {
        /* nothing to be done */
        return 0;
    }

    /* make mapping from old to new IDs */

    maplist = gretl_list_new(lmax - dmin + 1);
    if (maplist == NULL) {
        return E_ALLOC;
    }

    j = dmin;

    for (i=1; i<=maplist[0]; i++) {
        if (var_is_deleted(dlist, dmin, i)) {
            maplist[i] = -1;
        } else {
            maplist[i] = j++;
        }
    }

    /* use mapping to revise saved lists */
    for (j=0; j<n_vars; j++) {
        if (uvars[j]->type == GRETL_TYPE_LIST) {
            list = uvars[j]->ptr;
            if (list != NULL) {
                for (i=list[0]; i>0; i--) {
                    k = list[i] - dmin + 1;
                    if (k >= 1) {
                        if (maplist[k] == -1) {
                            gretl_list_delete_at_pos(list, i);
                        } else {
                            list[i] = maplist[k];
                        }
                    }
                }
            }
        }
    }

    free(maplist);

    return 0;
}

/**
 * gretl_lists_cleanup:
 *
 * Frees all resources associated with the internal
 * apparatus for saving and retrieving named lists.
 */

void gretl_lists_cleanup (void)
{
    real_destroy_user_vars_at_level(0,
                                    GRETL_TYPE_LIST,
                                    0);
}

/* below: serialization of user vars to XML, plus de-serialization
   from XML -- for use in GUI session mechanism
*/

static void write_user_scalars (PRN *prn)
{
    double x;
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_DOUBLE) {
            x = *(double *) uvars[i]->ptr;
            pprintf(prn, " <gretl-scalar name=\"%s\" value=\"", uvars[i]->name);
            serialize_scalar_value(x, prn);
            pputs(prn, "\"/>\n");
        }
    }
}

static void write_user_matrices (PRN *prn)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_MATRIX) {
            xml_put_user_matrix(uvars[i], prn);
        }
    }
}

static void write_user_lists (PRN *prn)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_LIST) {
            gretl_list_serialize(uvars[i]->ptr,
                                 uvars[i]->name,
                                 prn);
        }
    }
}

static void write_user_bundles (PRN *prn)
{
    int i;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_BUNDLE) {
            gretl_bundle_serialize(uvars[i]->ptr,
                                   uvars[i]->name,
                                   prn);
        }
    }
}

static int read_user_scalars (xmlDocPtr doc, xmlNodePtr cur)
{
    char *name, *val;
    double x;
    int n, err = 0;

    cur = cur->xmlChildrenNode;

    gretl_push_c_numeric_locale();

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-scalar")) {
            name = (char *) xmlGetProp(cur, (XUC) "name");
            val = (char *) xmlGetProp(cur, (XUC) "value");
            if (name == NULL || val == NULL) {
                err = 1;
            } else {
                n = sscanf(val, "%lf", &x);
                if (n < 1) {
#ifdef WIN32
                    x = win32_sscan_nonfinite(val, &err);
#else
                    x = NADBL;
#endif
                }
                err = gretl_scalar_add(name, x);
            }
            free(name);
            free(val);
        }
        cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    return err;
}

static int read_user_matrices (xmlDocPtr doc, xmlNodePtr cur)
{
    gretl_matrix *m;
    char *name;
    int err = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-matrix")) {
            name = (char *) xmlGetProp(cur, (XUC) "name");
            if (name == NULL) {
                err = 1;
            } else {
                m = gretl_xml_get_matrix(cur, doc, &err);
                if (m != NULL) {
                    err = user_var_add(name, GRETL_TYPE_MATRIX, m);
                }
                free(name);
            }
        }
        cur = cur->next;
    }

    return err;
}

static int read_user_lists (xmlDocPtr doc, xmlNodePtr cur)
{
    int *list;
    char *name;
    int err = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "list")) {
            if (!gretl_xml_get_prop_as_string(cur, "name", &name)) {
                err = E_DATA;
            } else {
                list = gretl_xml_get_list(cur, doc, &err);
                if (!err) {
                    err = user_var_add(name, GRETL_TYPE_LIST, list);
                }
                free(name);
            }
        }
        cur = cur->next;
    }

    return err;
}

static int read_user_bundles (xmlDocPtr doc, xmlNodePtr cur)
{
    int err = 0;

    gretl_push_c_numeric_locale();

    cur = cur->xmlChildrenNode;

    while (cur != NULL && !err) {
        if (!xmlStrcmp(cur->name, (XUC) "gretl-bundle")) {
            char *name = (char *) xmlGetProp(cur, (XUC) "name");

            if (name == NULL) {
                err = 1;
            } else {
                char *creator = NULL;
                gretl_bundle *b;

                b = gretl_bundle_deserialize(cur, doc, &err);
                if (!err) {
                    creator = (char *) xmlGetProp(cur, (XUC) "creator");
                    gretl_bundle_set_creator(b, creator);
                    err = user_var_add(name, GRETL_TYPE_BUNDLE, b);
                }
                free(name);
                free(creator);
            }
        }
        cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    return err;
}

typedef void (*xml_write_func) (PRN *);
typedef int (*xml_read_func) (xmlDocPtr, xmlNodePtr);

struct uvar_file_ {
    GretlType type;
    const char *typestr;
    xml_write_func write_func;
    xml_read_func read_func;
};

typedef struct uvar_file_ uvar_file;

static uvar_file uvar_files[] = {
    { GRETL_TYPE_DOUBLE, "scalars",  write_user_scalars,  read_user_scalars },
    { GRETL_TYPE_MATRIX, "matrices", write_user_matrices, read_user_matrices },
    { GRETL_TYPE_LIST,   "lists",    write_user_lists,    read_user_lists },
    { GRETL_TYPE_BUNDLE, "bundles",  write_user_bundles,  read_user_bundles }
};

int serialize_user_vars (const char *dirname)
{
    GretlType type;
    const char *typestr;
    void (*write_func)(PRN *prn);
    char path[MAXLEN];
    PRN *prn;
    int i, n, ni;
    int err = 0;

    n = sizeof uvar_files / sizeof uvar_files[0];

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        type = uvar_files[i].type;
        ni = user_var_count_for_type(type);
        if (ni > 0) {
            int errp = 0;

            typestr = uvar_files[i].typestr;
            sprintf(path, "%s%c%s.xml", dirname, SLASH, typestr);
            write_func = uvar_files[i].write_func;
            prn = gretl_print_new_with_filename(path, &errp);
            if (prn == NULL) {
                err++;
                continue;
            }
            gretl_xml_header(prn);
            pprintf(prn, "<gretl-%s count=\"%d\">\n", typestr, ni);
            (*write_func)(prn);
            pprintf(prn, "</gretl-%s>\n", typestr);
            gretl_print_destroy(prn);
        }
    }

    gretl_pop_c_numeric_locale();

    if (err > 0) {
        fprintf(stderr, "Failed writing %d user_var files\n", err);
        err = E_FOPEN;
    }

    return err;
}

#define UDEBUG 0

int deserialize_user_vars (const char *dirname)
{
    xmlDocPtr doc = NULL;
    xmlNodePtr cur = NULL;
    const char *typestr;
    int (*read_func)(xmlDocPtr, xmlNodePtr);
    char root_name[16];
    char path[MAXLEN];
    FILE *fp;
    int i, n;
    int n_failed = 0;
    int err = 0;

    n = sizeof uvar_files / sizeof uvar_files[0];

#if UDEBUG
    fprintf(stderr, "deserialize_user_vars:\n");
#endif

    for (i=0; i<n; i++) {
        int err_i = 0;

        typestr = uvar_files[i].typestr;
        sprintf(path, "%s%c%s.xml", dirname, SLASH, typestr);

#if UDEBUG
        fprintf(stderr, " checking for '%s.xml'\n", typestr);
#endif
        fp = gretl_fopen(path, "r");
        if (fp == NULL) {
            /* OK, no user-vars of this type */
#if UDEBUG
            fprintf(stderr, "  not found\n");
#endif
            continue;
        }
        fclose(fp);
        sprintf(root_name, "gretl-%s", typestr);
        err_i = gretl_xml_open_doc_root(path, root_name, &doc, &cur);
        if (!err_i) {
            read_func = uvar_files[i].read_func;
#if UDEBUG
            fprintf(stderr, "  found, reading...\n");
#endif
            err_i = read_func(doc, cur);
#if UDEBUG
            fprintf(stderr, "  done.\n");
#endif
        }
        if (doc != NULL) {
            xmlFreeDoc(doc);
            doc = NULL;
        }
        if (err_i) {
            n_failed++;
            if (!err) {
                err = err_i;
            }
        }
    }

    if (n_failed > 0) {
        fprintf(stderr, "Failed reading %d user_var files\n", n_failed);
    }

    return err;
}

int print_user_var_by_name (const char *name,
                            const DATASET *dset,
                            gretlopt opt,
                            PRN *prn)
{
    user_var *u = get_user_var_by_name(name);
    int err = 0;

    if (u == NULL || u->ptr == NULL) {
        return E_DATA;
    }

    if (u->type == GRETL_TYPE_DOUBLE) {
        print_scalar_by_name(name, prn);
    } else if (u->type == GRETL_TYPE_MATRIX) {
        gretl_matrix *tmp = u->ptr;

        if (tmp->is_complex || opt & OPT_C) {
            err = gretl_cmatrix_print(tmp, name, prn);
	} else {
            gretl_matrix_print_to_prn(tmp, name, prn);
        }
    } else if (u->type == GRETL_TYPE_BUNDLE) {
        if (opt & OPT_T) {
            gretl_bundle_print_tree(u->ptr, prn);
        } else {
            gretl_bundle_print(u->ptr, prn);
        }
    } else if (u->type == GRETL_TYPE_ARRAY) {
        gretl_array_print(u->ptr, prn);
    } else if (u->type == GRETL_TYPE_LIST) {
        gretl_list_print(u->ptr, dset, prn);
    } else if (u->type == GRETL_TYPE_STRING) {
        pputs(prn, (char *) u->ptr);
        pputc(prn, '\n');
    }

    return err;
}

static int uvar_type_match (user_var *u, GretlType t)
{
    if (u->type == t) {
        return 1;
    } else if (u->type == GRETL_TYPE_ARRAY &&
               gretl_array_type(t)) {
        return t == gretl_array_get_type(u->ptr);
    } else {
        return 0;
    }
}

int list_user_vars_of_type (const DATASET *dset,
                            PRN *prn)
{
    const char *typename;
    GretlType t;

    typename = get_optval_string(VARLIST, OPT_T);
    if (typename == NULL) {
        return E_INVARG;
    }

    if (!strcmp(typename, "accessor")) {
        list_ok_dollar_vars((DATASET *) dset, prn);
        return 0;
    }

    t = gretl_type_from_string(typename);
    if (t == GRETL_TYPE_NONE) {
        return E_INVARG;
    }

    if (t == GRETL_TYPE_SERIES) {
        list_series(dset, OPT_NONE, prn);
    } else if (t == GRETL_TYPE_DOUBLE) {
        print_scalars(prn);
    } else if (t == GRETL_TYPE_LIST ||
               t == GRETL_TYPE_MATRIX ||
               t == GRETL_TYPE_BUNDLE ||
               t == GRETL_TYPE_ARRAY ||
               t == GRETL_TYPE_STRING ||
               gretl_array_type(t)) {
        int i, n = 0;

        pprintf(prn, _("variables of type %s:"), typename);
        for (i=0; i<n_vars; i++) {
            if (uvar_type_match(uvars[i], t)) {
                if (n == 0) {
                    pputc(prn, '\n');
                }
                if (uvars[i]->name[0] == '\0') {
                    pputs(prn, _("  (unnamed)\n"));
                } else if (t == GRETL_TYPE_ARRAY) {
                    GretlType at = gretl_array_get_type(uvars[i]->ptr);

                    pprintf(prn, "  %s (%s)\n", uvars[i]->name,
                            gretl_type_get_name(at));
                } else {
                    pprintf(prn, "  %s\n", uvars[i]->name);
                }
                n++;
            }
        }
        if (n == 0) {
            pprintf(prn, " %s\n", _("none"));
        }
        pputc(prn, '\n');
    } else {
        return E_INVARG;
    }

    return 0;
}

int leads_midas_list (int ID, const DATASET *dset,
                      char *listname)
{
    int level = gretl_function_depth();
    int *list;
    int i, ret = 0;

    for (i=0; i<n_vars && !ret; i++) {
        if (uvars[i]->type == GRETL_TYPE_LIST &&
            uvars[i]->level == level) {
            list = uvars[i]->ptr;
            if (list[0] > 2 && list[1] == ID) {
                ret = gretl_is_midas_list(list, dset);
                if (ret && listname != NULL) {
                    strcpy(listname, uvars[i]->name);
                }
            }
        }
    }

    return ret;
}

int in_midas_list (int ID, const DATASET *dset,
                   char *listname)
{
    int level = gretl_function_depth();
    int *list;
    int i, ret = 0;

    for (i=0; i<n_vars && !ret; i++) {
        if (uvars[i]->type == GRETL_TYPE_LIST &&
            uvars[i]->level == level) {
            list = uvars[i]->ptr;
            if (list[0] > 2 && in_gretl_list(list, ID)) {
                ret = gretl_is_midas_list(list, dset);
                if (ret && listname != NULL) {
                    strcpy(listname, uvars[i]->name);
                }
            }
        }
    }

    return ret;
}

const char *get_listname_by_consecutive_content (int l0, int l1)
{
    int level = gretl_function_depth();
    const char *ret = NULL;
    int i, j, *list;

    for (i=0; i<n_vars; i++) {
        if (uvars[i]->type == GRETL_TYPE_LIST &&
            uvars[i]->level == level) {
            list = uvars[i]->ptr;
            if (list[0] == l0 && list[1] == l1) {
                int found = 1;

                for (j=2; j<=l0; j++) {
                    if (list[j] != list[j-1] + 1) {
                        found = 0;
                        break;
                    }
                }
                if (found) {
                    return uvars[i]->name;
                }
            }
        }
    }

    return ret;
}

/* Dropping terms from the list @targ inside a function:
   this is tricky with regard to the auto-generated
   "time" variable, which may be a member of a list that
   was passed as an argument yet not "visible" by name
   within the function. Here we attempt to fix this by
   transcribing the ID number for "time" from the caller's
   namespace into the @drop list -- if the latter is
   trying to drop this variable.
*/

static void check_auto_time_var (const int *targ, int *drop,
                                 const DATASET *dset)
{
    int i, vi, tnum = 0;

    for (i=1; i<=targ[0]; i++) {
        vi = targ[i];
        if (!strcmp(dset->varname[vi], "time")) {
            tnum = vi;
            break;
        }
    }

    if (tnum > 0) {
        for (i=drop[0]; i>0; i--) {
            vi = drop[i];
            if (!strcmp(dset->varname[vi], "time")) {
                drop[i] = tnum;
            }
        }
    }
}

/* functions called from geneval.c when "editing" a named list */

int user_list_append (user_var *uvar, const int *add)
{
    int err = 0;

    if (uvar == NULL || user_var_get_type(uvar) != GRETL_TYPE_LIST) {
        err = E_DATA;
    } else {
        const int *list = user_var_get_value(uvar);
        int *tmp = gretl_list_copy(list);

        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            err = gretl_list_add_list(&tmp, add);
            if (!err) {
                user_var_replace_value(uvar, tmp, GRETL_TYPE_LIST);
            }
        }
    }

    return err;
}

int user_list_subtract (user_var *uvar, int *sub,
                        const DATASET *dset)
{
    int err = 0;

    if (uvar == NULL || user_var_get_type(uvar) != GRETL_TYPE_LIST) {
        err = E_DATA;
    } else {
        const int *list = user_var_get_value(uvar);
        int *tmp;

        if (gretl_function_depth() > 0) {
            check_auto_time_var(list, sub, dset);
        }
        tmp = gretl_list_drop(list, sub, &err);
        if (!err) {
            user_var_replace_value(uvar, tmp, GRETL_TYPE_LIST);
        }
    }

    return err;
}

int user_list_replace (user_var *uvar, const int *src)
{
    int err = 0;

    if (uvar == NULL || user_var_get_type(uvar) != GRETL_TYPE_LIST) {
        err = E_DATA;
    } else {
        int *tmp = gretl_list_copy(src);

        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            user_var_replace_value(uvar, tmp, GRETL_TYPE_LIST);
        }
    }

    return err;
}

/**
 * remember_list:
 * @list: array of integers, the first element being a count
 * of the following elements.
 * @name: name to be given to the list.
 * @prn: printing struct.
 *
 * Adds a copy of @list to the stack of saved lists and associates
 * it with @name, unless there is already a list with the given
 * name in which case the original list is replaced.  A status
 * message is printed to @prn.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int remember_list (const int *list, const char *name, PRN *prn)
{
    int *lcpy = gretl_list_copy(list);
    int err = 0;

    if (lcpy == NULL) {
        err = (list == NULL)? E_DATA : E_ALLOC;
    } else {
        user_var *orig;

        orig = get_user_var_of_type_by_name(name, GRETL_TYPE_LIST);

        if (orig != NULL) {
            /* replace existing list of same name */
            user_var_replace_value(orig, lcpy, GRETL_TYPE_LIST);
            if (prn != NULL && gretl_messages_on()) {
                pprintf(prn, _("Replaced list '%s'\n"), name);
            }
        } else {
            err = user_var_add(name, GRETL_TYPE_LIST, lcpy);
            if (!err && prn != NULL && gretl_messages_on()) {
                pprintf(prn, _("Added list '%s'\n"), name);
            }
        }
    }

    return err;
}

/**
 * get_list_by_name:
 * @name: the name of the list to be found.
 *
 * Looks up @name in the stack of saved variables, at the current level
 * of function execution, and retrieves the associated list.
 *
 * Returns: the list, or NULL if the lookup fails.
 */

int *get_list_by_name (const char *name)
{
    user_var *u;
    int *ret = NULL;

    u = get_user_var_of_type_by_name(name, GRETL_TYPE_LIST);

    if (u != NULL) {
        ret = user_var_get_value(u);
    }

    return ret;
}

/**
 * gretl_is_list:
 * @name: the name to test.
 *
 * Returns: 1 if @name is the name of a saved list, 0
 * otherwise.
 */

int gretl_is_list (const char *name)
{
    return get_user_var_of_type_by_name(name, GRETL_TYPE_LIST) != NULL;
}
