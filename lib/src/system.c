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

/* mechanisms for defining and handling systems of equations */

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "system.h"
#include "objstack.h"
#include "usermat.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "tsls.h"
#include "uservar.h"

#define SYSDEBUG 0

enum {
    OP_PLUS,
    OP_MINUS
} identity_ops;

enum {
    ENDOG_LIST,
    INSTR_LIST
} aux_list_types;

enum {
    SYS_TEST_NONE,
    SYS_TEST_LR,
    SYS_TEST_F,
    SYS_TEST_NOTIMP
} system_test_types;

enum {
    SYS_UHAT,
    SYS_YHAT
};

struct id_atom_ {
    int op;         /* operator (plus or minus) */
    int varnum;     /* ID number of variable to right of operator */
};

struct identity_ {
    int n_atoms;    /* number of "atomic" elements in identity */
    int depvar;     /* LHS variable in indentity */
    id_atom *atoms; /* pointer to RHS "atoms" */
};

struct predet_ {
    int id;         /* ID number of lag variable */
    int src;        /* ID number of "parent" */
    int lag;        /* lag order */
};

struct liml_data_ {
    double *lmin;   /* min. eigenvalues */
    double *ll;     /* per-equation log likelihood */
    int *idf;       /* overidentification df */
};

const char *gretl_system_method_strings[] = {
    "sur",
    "3sls",
    "fiml",
    "liml",
    "ols",
    "tsls",
    "wls",
    NULL
};

const char *gretl_system_short_strings[] = {
    N_("SUR"),
    N_("3SLS"),
    N_("FIML"),
    N_("LIML"),
    N_("OLS"),
    N_("TSLS"),
    N_("WLS"),
    NULL
};

const char *gretl_system_long_strings[] = {
    N_("Seemingly Unrelated Regressions"),
    N_("Three-Stage Least Squares"),
    N_("Full Information Maximum Likelihood"),
    N_("Limited Information Maximum Likelihood"),
    N_("Ordinary Least Squares"),
    N_("Two-Stage Least Squares"),
    N_("Weighted Least Squares"),
    NULL
};

const char *nosystem = N_("No system of equations has been defined");
const char *badsystem = N_("Unrecognized equation system type");
const char *toofew = N_("An equation system must have at least two equations");

static void destroy_ident (identity *ident);
static int
add_predet_to_sys (equation_system *sys, const DATASET *dset,
                   int id, int src, int lag);
static int sys_check_lists (equation_system *sys, const DATASET *dset);

#define sys_anonymous(s) (strcmp(s, "$system") == 0)

static GList *sysstack;

static equation_system *get_anon_system_at_depth (int fd)
{
    equation_system *sys;
    GList *tmp = sysstack;

    while (tmp != NULL) {
        sys = tmp->data;
        if (sys->fd == fd) {
            return sys;
        }
        tmp = tmp->next;
    }

    return NULL;
}

equation_system *get_anonymous_equation_system (void)
{
    int fd = gretl_function_depth();

    return get_anon_system_at_depth(fd);
}

/* The use-case for push_anon_system() is when a user wants to define
   a system and then estimate it by more than one method (or set it as
   a target for "restrict" before estimation) but does not care to
   give it a specific name via the "name <- system" mechanism -- or
   we're inside a user-defined function where this naming mechanism is
   not available.

   In this case we save the system "anonymously" instead (under a name
   of "$system").
*/

static void push_anon_system (equation_system *sys)
{
    equation_system *old = get_anon_system_at_depth(sys->fd);

    if (old != NULL) {
        sysstack = g_list_remove(sysstack, old);
        gretl_object_unref(old, GRETL_OBJ_SYS);
    }

    gretl_object_ref(sys, GRETL_OBJ_SYS);
    sysstack = g_list_append(sysstack, sys);
}

/* For use when terminating function execution: if an
   anonymous equation system was set up at the given
   level of function execution, trash it.
*/

void delete_anonymous_equation_system (int level)
{
    equation_system *sys = get_anon_system_at_depth(level);

    if (sys != NULL) {
        sysstack = g_list_remove(sysstack, sys);
        gretl_object_unref(sys, GRETL_OBJ_SYS);
    }
}

static void
print_system_equation (const int *list, const DATASET *dset,
                       PRN *prn)
{
    int i, v;

    pputs(prn, "equation");

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (v == LISTSEP) {
            pputs(prn, " ;");
        } else if (v > 0 && v < dset->v) {
            pprintf(prn, " %s", dset->varname[v]);
        } else {
            pprintf(prn, " %d", v);
        }
    }

    pputc(prn, '\n');
}

static void
print_system_identity (const identity *ident, const DATASET *dset,
                       gretlopt opt, PRN *prn)
{
    int i;

    if (opt & OPT_H) {
        pprintf(prn, "Identity: %s = %s ",
                dset->varname[ident->depvar],
                dset->varname[ident->atoms[0].varnum]);
    } else {
        pprintf(prn, "identity %s = %s ",
                dset->varname[ident->depvar],
                dset->varname[ident->atoms[0].varnum]);
    }

    for (i=1; i<ident->n_atoms; i++) {
        pprintf(prn, "%c %s ", (ident->atoms[i].op == OP_PLUS)? '+' : '-',
                dset->varname[ident->atoms[i].varnum]);
    }

    pputc(prn, '\n');
}

static int sys_max_predet_lag (const equation_system *sys)
{
    int i, m = 0;

    for (i=0; i<sys->plist[0]; i++) {
        if (sys->pre_vars[i].lag > m) {
            m = sys->pre_vars[i].lag;
        }
    }

    return m;
}

static int get_predet_parent (const equation_system *sys, int v, int *lag)
{
    int i;

    for (i=0; i<sys->plist[0]; i++) {
        if (v == sys->pre_vars[i].id) {
            *lag = sys->pre_vars[i].lag;
            return sys->pre_vars[i].src;
        }
    }

    *lag = 0;

    return -1;
}

/**
 * print_equation_system_info:
 * @sys: gretl equation system.
 * @dset: dataset information.
 * @opt: use OPT_H for printing in a form designed to appear
 * in the header when results of estimation are printed.
 * @prn: printing struct.
 *
 * Prints details of an equation system to @prn.
 */

void
print_equation_system_info (const equation_system *sys,
                            const DATASET *dset,
                            gretlopt opt, PRN *prn)
{
    int header = (opt & OPT_H);
    int i, vi, lag;

    if (header && sys->name != NULL && !sys_anonymous(sys->name)) {
        pprintf(prn, "%s %s\n", _("Equation system"), sys->name);
    }

    if (!header) {
        for (i=0; i<sys->neqns; i++) {
            print_system_equation(sys->lists[i], dset, prn);
        }
    }

    for (i=0; i<sys->nidents; i++) {
        print_system_identity(sys->idents[i], dset, opt, prn);
    }

    if (sys->ylist != NULL) {
        pputs(prn, (header)? _("Endogenous variables:") : "endog");
        for (i=1; i<=sys->ylist[0]; i++) {
            vi = sys->ylist[i];
            pprintf(prn, " %s", dset->varname[vi]);
        }
        pputc(prn, '\n');
    }

    if (header) {
        if (sys->pre_vars != NULL) {
            pputs(prn, _("Predetermined variables:"));
            for (i=0; i<sys->plist[0]; i++) {
                vi = sys->pre_vars[i].src;
                lag = sys->pre_vars[i].lag;
                pprintf(prn, " %s(-%d)", dset->varname[vi], lag);
            }
            pputc(prn, '\n');
        }
        if (sys->ilist != NULL && sys->ilist[0] > sys->plist[0]) {
            pputs(prn, _("Exogenous variables:"));
            for (i=1; i<=sys->ilist[0]; i++) {
                vi = sys->ilist[i];
                if (!in_gretl_list(sys->plist, vi)) {
                    pprintf(prn, " %s", dset->varname[vi]);
                }
            }
            pputc(prn, '\n');
        }
    } else if (sys->ilist != NULL) {
        pputs(prn, "instr");
        for (i=1; i<=sys->ilist[0]; i++) {
            vi = sys->ilist[i];
            pprintf(prn, " %s", dset->varname[vi]);
        }
        pputc(prn, '\n');
    }
}

int system_method_from_string (const char *s)
{
    int i = 0;

    while (gretl_system_method_strings[i] != NULL) {
        if (!strcmp(s, gretl_system_method_strings[i])) {
            return i;
        }
        i++;
    }

    return i;
}

const char *system_method_full_string (int method)
{
    if (method >= SYS_METHOD_SUR && method < SYS_METHOD_MAX) {
        return gretl_system_long_strings[method];
    } else {
        return NULL;
    }
}

const char *system_method_short_string (int method)
{
    if (method >= SYS_METHOD_SUR && method < SYS_METHOD_MAX) {
        return gretl_system_method_strings[method];
    } else {
        return NULL;
    }
}

static equation_system *
equation_system_new (int method, const char *name, int *err)
{
    equation_system *sys;

    if (method < 0 && name == NULL) {
        *err = E_DATA;
        return NULL;
    }

    sys = malloc(sizeof *sys);
    if (sys == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    sys->name = NULL;

    if (name != NULL) {
        equation_system_set_name(sys, name);
    }

    sys->refcount = 0;
    sys->fd = gretl_function_depth();
    sys->method = method;

    sys->t1 = sys->t2 = 0;
    sys->smpl_t1 = sys->smpl_t2 = 0;
    sys->T = sys->df = 0;

    sys->neqns = 0;
    sys->nidents = 0;

    sys->q = sys->R = NULL;

    sys->order = sys->flags = sys->iters = 0;

    sys->ll = sys->llu = sys->ldet = NADBL;
    sys->X2 = NADBL;
    sys->ess = NADBL;
    sys->diag_test = 0.0;
    sys->bdiff = 0.0;

    sys->b = sys->vcv = NULL;
    sys->E = sys->S = NULL;
    sys->yhat = NULL;

    sys->Gamma = NULL;
    sys->A = sys->B = NULL;
    sys->F = sys->Sr = NULL;

    sys->tslists = sys->lists = NULL;
    sys->ilist = sys->ylist = NULL;
    sys->plist = sys->xlist = NULL;
    sys->biglist = NULL;

    sys->pre_vars = NULL;
    sys->idents = NULL;

    sys->models = NULL;
    sys->ldata = NULL;

    return sys;
}

static void system_clear_restrictions (equation_system *sys)
{
    if (sys->R != NULL) {
        free(sys->R);
        sys->R = NULL;
    }

    if (sys->q != NULL) {
        free(sys->q);
        sys->q = NULL;
    }

    sys->flags &= ~SYSTEM_RESTRICT;
}

static void system_clear_results (equation_system *sys)
{
    sys->iters = 0;

    sys->t1 = sys->t2 = 0;
    sys->smpl_t1 = sys->smpl_t2 = 0;
    sys->T = sys->df = 0;

    sys->ll = NADBL;
    sys->llu = NADBL;
    sys->ldet = NADBL;
    sys->X2 = NADBL;
    sys->ess = NADBL;

    sys->diag_test = 0.0;
    sys->bdiff = 0.0;

    if (sys->b != NULL) {
        gretl_matrix_replace(&sys->b, NULL);
    }

    if (sys->vcv != NULL) {
        gretl_matrix_replace(&sys->vcv, NULL);
    }

    if (sys->S != NULL) {
        gretl_matrix_replace(&sys->S, NULL);
    }

    if (sys->E != NULL) {
        gretl_matrix_replace(&sys->E, NULL);
    }

    if (sys->yhat != NULL) {
        gretl_matrix_replace(&sys->yhat, NULL);
    }

    if (sys->Gamma != NULL) {
        gretl_matrix_replace(&sys->Gamma, NULL);
    }

    if (sys->B != NULL) {
        gretl_matrix_replace(&sys->B, NULL);
    }

    if (sys->A != NULL) {
        gretl_matrix_replace(&sys->A, NULL);
    }

    if (sys->Sr != NULL) {
        gretl_matrix_replace(&sys->Sr, NULL);
    }

    if (sys->F != NULL) {
        gretl_matrix_replace(&sys->F, NULL);
    }

    if (sys->ldata != NULL) {
        free(sys->ldata->lmin);
        free(sys->ldata->ll);
        free(sys->ldata->idf);
        free(sys->ldata);
        sys->ldata = NULL;
    }
}

void equation_system_destroy (equation_system *sys)
{
    int i;

#if SYSDEBUG
    fprintf(stderr, "equation_system_destroy: %p\n", (void *) sys);
#endif

    if (sys == NULL || sys->lists == NULL) {
        return;
    }

    sys->refcount -= 1;
    if (sys->refcount > 0) {
        return;
    }

    for (i=0; i<sys->neqns; i++) {
        free(sys->lists[i]);
    }

    free(sys->lists);
    sys->lists = NULL;

    for (i=0; i<sys->nidents; i++) {
        destroy_ident(sys->idents[i]);
    }

    free(sys->idents);

    free(sys->ylist);
    free(sys->ilist);
    free(sys->xlist);
    free(sys->plist);
    free(sys->biglist);

    free(sys->pre_vars);

    free(sys->name);

    gretl_matrix_free(sys->R);
    gretl_matrix_free(sys->q);

    system_clear_results(sys);

    free(sys);
}

#define SYS_PRUNE_LISTS 1

#if SYS_PRUNE_LISTS

/* 2020-12-09: When a system is defined using TSLS-style
   regression lists, with instruments specified per
   equation, there were problems estimating the system
   via estimators that don't take such a specification.
   The apparatus below is designed to handle this: we
   remove the instruments immediately prior to estimation,
   then restore them immediately afterwards.
*/

/* Make a copy of the incoming per-equation regression
   lists, with any instruments removed.
*/

static int lists_strip_insts (equation_system *sys)
{
    int i, *eqlist;

    sys->tslists = malloc(sys->neqns * sizeof *sys->tslists);

    for (i=0; i<sys->neqns; i++) {
	sys->tslists[i] = sys->lists[i];
	if (gretl_list_has_separator(sys->lists[i])) {
	    gretl_list_split_on_separator(sys->tslists[i], &eqlist, NULL);
	    sys->lists[i] = eqlist;
	} else {
	    sys->lists[i] = gretl_list_copy(sys->tslists[i]);
	}
    }

    return 0;
}

/* Restore the original per-equation regression lists,
   in case instruments were removed for estimation by
   OLS, WLS, SUR or FIML.
*/

static void lists_restore_insts (equation_system *sys)
{
    int i;

    for (i=0; i<sys->neqns; i++) {
	free(sys->lists[i]);
	sys->lists[i] = sys->tslists[i];
	sys->tslists[i] = NULL;
    }

    free(sys->tslists);
    sys->tslists = NULL;
}

#endif /* SYS_PRUNE_LISTS */

static int sys_rearrange_eqn_lists (equation_system *sys,
                                    const DATASET *dset)
{
    int i, err = 0;

    for (i=0; i<sys->neqns; i++) {
        reglist_check_for_const(sys->lists[i], dset);
    }

#if SYS_PRUNE_LISTS
    if (sys->method != SYS_METHOD_TSLS &&
        sys->method != SYS_METHOD_3SLS &&
	sys->method != SYS_METHOD_LIML) {
        /* we can't have ';' in equation lists */
	for (i=0; i<sys->neqns; i++) {
	    if (gretl_list_has_separator(sys->lists[i])) {
		lists_strip_insts(sys);
		break;
	    }
	}
    }
#else
    /* the status quo prior to 2020-12-09 */
    if (sys->method != SYS_METHOD_TSLS &&
        sys->method != SYS_METHOD_3SLS) {
        /* we can't have ';' in equation lists */
        int j;

        for (i=0; i<sys->neqns && !err; i++) {
            for (j=0; j<=sys->lists[i][0] && !err; j++) {
                if (sys->lists[i][j] == LISTSEP) {
                    gretl_errmsg_sprintf("%s: tsls-style lists not supported",
                                         gretl_system_short_strings[sys->method]);
                    err = E_DATA;
                }
            }
        }
    }
#endif

    return err;
}

/* Form a list from row @i of matrix @m, ignoring trailing
   zero elements. We check that series ID numbers are in
   bounds, and also that the list does not contain
   duplicated elements.
*/

static int *matrix_row_to_list (const gretl_matrix *m, int i,
                                const DATASET *dset,
                                int *err)
{
    int *list = NULL;
    int j, k, n = m->cols;

    /* figure how many elements to read */
    for (j=m->cols-1; j>=0; j--) {
        if (gretl_matrix_get(m, i, j) == 0) {
            n--;
        } else {
            break;
        }
    }

    if (n == 0) {
        *err = E_DATA;
        return NULL;
    }

    /* check for out-of-bounds values */
    for (j=0; j<n; j++) {
        k = gretl_matrix_get(m, i, j);
        if (k < 0 || k >= dset->v) {
            *err = E_UNKVAR;
            return NULL;
        }
    }

    /* construct the list */
    list = gretl_list_new(n);
    if (list == NULL) {
        *err = E_ALLOC;
    } else {
        k = 1;
        for (j=0; j<n; j++) {
            list[k++] = (int) gretl_matrix_get(m, i, j);
        }
    }

    if (!*err) {
        k = gretl_list_duplicates(list, EQUATION);
        if (k >= 0) {
            gretl_errmsg_sprintf(_("variable %d duplicated in the "
                                   "command list."), k);
            *err = E_DATA;
            free(list);
            list = NULL;
        }
    }

    return list;
}

static int sys_check_eqn_list (const int *list)
{
    int dupv = gretl_list_duplicates(list, EQUATION);

    if (dupv >= 0) {
        gretl_errmsg_sprintf(_("variable %d duplicated in the "
                               "command list."), dupv);
        return E_DATA;
    } else {
        return 0;
    }
}

/* @LY is a simple list of g regressands; either @LX is
   simple list of common regressors or @AX is an array of
   lists. We get either @LX or @AX as non-NULL, never both.
*/

static int add_equations_from_lists (equation_system *sys,
                                     const int *LY,
                                     const int *LX,
                                     gretl_array *AX,
                                     const DATASET *dset)
{
    int nx, n = sys->neqns;
    int g = LY[0];
    int i, j, err = 0;

    if (AX != NULL && gretl_array_get_length(AX) != g) {
        err = E_ARGS;
    } else {
        int **L = realloc(sys->lists, (n + g) * sizeof *sys->lists);

        if (L == NULL) {
            err = E_ALLOC;
        } else {
            sys->lists = L;
            for (i=0; i<g; i++) {
                sys->lists[n+i] = NULL;
            }
        }
    }

    if (err) {
        return err;
    }

    if (AX != NULL) {
        /* handle an array of regressor lists */
        int *rhs, *list;

        for (i=0; i<g && !err; i++) {
            list = NULL;
            rhs = gretl_array_get_element(AX, i, NULL, &err);
            if (!err) {
                nx = rhs[0];
                list = gretl_list_new(1 + nx);
                if (list == NULL) {
                    err = E_ALLOC;
                }
            }
            if (!err) {
                list[1] = LY[i+1];
                for (j=1; j<=nx; j++) {
                    list[j+1] = rhs[j];
                }
                err = sys_check_eqn_list(list);
            }
            if (err) {
                free(list);
            } else {
                sys->lists[n++] = list;
            }
        }
    } else {
        /* handle a list of common regressors */
        int nx = LX[0];
        int *list, *l0 = NULL;

        for (i=0; i<g && !err; i++) {
            if (i == 0) {
                l0 = list = gretl_list_new(1 + nx);
            } else {
                list = gretl_list_copy(l0);
            }
            if (list == NULL) {
                err = E_ALLOC;
            } else {
                list[1] = LY[i+1];
                if (i == 0) {
                    for (j=1; j<=nx; j++) {
                        list[j+1] = LX[j];
                    }
                }
                err = sys_check_eqn_list(list);
            }
            if (err) {
                free(list);
            } else {
                sys->lists[n++] = list;
            }
        }
    }

    if (!err) {
        sys->neqns += g;
    }

    return err;
}

static int add_equations_from_matrix (equation_system *sys,
                                      const gretl_matrix *m,
                                      const DATASET *dset)
{
    int *list, **L;
    int n = sys->neqns;
    int i, err = 0;

    L = realloc(sys->lists, (n + m->rows) * sizeof *sys->lists);
    if (L == NULL) {
        return E_ALLOC;
    }

    sys->lists = L;

    for (i=0; i<m->rows && !err; i++) {
        list = matrix_row_to_list(m, i, dset, &err);
        if (!err) {
            err = sys_check_eqn_list(list);
        }
        if (!err) {
            sys->lists[n++] = list;
        }
    }

    if (!err) {
        sys->neqns += m->rows;
    }

    return err;
}

/**
 * equation_system_append_multi:
 * @sys: initialized equation system.
 * @parm1: the name of a pre-defined matrix or list.
 * @parm2: the name of list or array of lists (or %NULL).
 * @dset: pointer to dataset struct.
 *
 * Adds one or more equations to @sys as follows.
 *
 * If @parm2 is %NULL then @parm1 is taken to be the
 * name of a matrix, and we interpret the rows of the
 * specified matrix as lists. Lists of differing length
 * can be accommodated by padding unused trailing elements
 * of short rows with zeros.
 *
 * If @parm2 is non-%NULL then @parm1 is taken to be the name
 * of a list of regressands, one per equation; @parm2,
 * pertaining to regressors, should be the name of a list
 * if the regressors are in common across the equations
 * or an array of lists (one per equation) otherwise.
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int equation_system_append_multi (equation_system *sys,
                                  const char *parm1,
                                  const char *parm2,
                                  const DATASET *dset)
{
    int err = 0;

    if (sys == NULL) {
        gretl_errmsg_set(_(nosystem));
        return E_DATA;
    } else if (parm1 == NULL) {
        return E_ARGS;
    }

    if (parm2 != NULL) {
        /* look for two lists, or list plus array */
        const int *LY = get_list_by_name(parm1);
        const int *LX = get_list_by_name(parm2);
        gretl_array *AX = NULL;

        if (LY == NULL) {
            gretl_errmsg_sprintf(_("'%s': no such list"), parm1);
            err = E_UNKVAR;
        } else if (LX == NULL) {
            /* try for an array on the right */
            AX = get_array_by_name(parm2);
            if (gretl_array_get_type(AX) != GRETL_TYPE_LISTS) {
                gretl_errmsg_sprintf(_("'%s': not a list or array of lists"),
                                     parm2);
                err = E_UNKVAR;
            }
        }
        if (!err) {
            err = add_equations_from_lists(sys, LY, LX, AX, dset);
        }
    } else {
        /* look for one matrix */
        const gretl_matrix *m = get_matrix_by_name(parm1);

        if (m == NULL) {
            gretl_errmsg_sprintf(_("'%s': no such matrix"), parm1);
            err = E_UNKVAR;
        } else if (m->rows == 0 || m->cols == 0) {
            err = E_DATA;
        } else {
            err = add_equations_from_matrix(sys, m, dset);
        }
    }

    if (err) {
        equation_system_destroy(sys);
    }

    return err;
}

/**
 * equation_system_append:
 * @sys: initialized equation system.
 * @list: list containing dependent variable and regressors.
 *
 * Adds an equation, as represented by @list, to @sys.
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int equation_system_append (equation_system *sys, const int *list)
{
    int err = 0;

    if (sys == NULL) {
        gretl_errmsg_set(_(nosystem));
        err = E_DATA;
    } else {
        int n = sys->neqns;
        int **lists;

        lists = realloc(sys->lists, (n + 1) * sizeof *sys->lists);
        if (lists == NULL) {
            err = E_ALLOC;
        } else {
            sys->lists = lists;
            sys->lists[n] = gretl_list_copy(list);
#if SYSDEBUG
            fprintf(stderr, "equation_system_append: added list %d\n", n);
            printlist(list, "newly added list");
#endif
            if (sys->lists[n] == NULL) {
                err = E_ALLOC;
            } else {
                sys->neqns += 1;
            }
        }

        if (err) {
            equation_system_destroy(sys);
        }
    }

    return err;
}

/* retrieve the name -- possibly quoted with embedded spaces -- for
   an equation system */

char *get_system_name_from_line (const char *s)
{
    const char *p = NULL;
    char *name = NULL;
    int pchars = 0;

    if (!strncmp(s, "method", 6)) {
        /* skip "method = whatever", with possible
           spaces around '='
        */
        char c = *(s+6);

        if (c == ' ' || c == '=') {
            p = s + 6;
            p += strspn(s, " ");
            if (*p == '=') p++;
            p += strspn(s, " ");
            p += strcspn(s, " ");
            p += strspn(s, " ");
            s = p;
        }
    }

#if SYSDEBUG
    fprintf(stderr, "get_system_name_from_line: s = '%s'\n", s);
#endif

    if (*s == '"') {
        if (*(s + 1) != '\0') s++;
        p = s;
        while (*p && *p != '"') {
            if (!isspace((unsigned char) *p)) pchars++;
            p++;
        }
        if (*p != '"') {
            /* no closing quote */
            pchars = 0;
        }
    } else {
        p = s;
        while (*p && !isspace((unsigned char) *p)) {
            pchars++;
            p++;
        }
    }

    if (pchars > 0) {
        name = gretl_strndup(s, p - s);
    }

    return name;
}

/* parse a system estimation method out of parameter @s */

static int sys_get_estimator (const char *s)
{
    const char *p = strstr(s, "method");
    int method = -1;

    if (p != NULL) {
        p += 6;
        p += strspn(s, " ");
        if (*p == '=') {
            char mstr[9];

            p++;
            p += strspn(s, " ");
            if (sscanf(p, "%8s", mstr) == 1) {
                gretl_lower(mstr);
                method = system_method_from_string(mstr);
            }
        }
    }

    return method;
}

static int parse_sys_start_param (const char *s,
                                  int *method,
                                  char **name)
{
    int len = strlen(s);
    int err = 0;

    if (len > 7 && !strncmp(s, "method=", 7)) {
        char mstr[9];

        *mstr = '\0';
        strncat(mstr, s + 7, 8);
        gretl_lower(mstr);
        *method = system_method_from_string(mstr);
        if (*method == SYS_METHOD_MAX) {
            /* invalid method was given */
            gretl_errmsg_set(_(badsystem));
            err = E_DATA;
        }
    } else if (len > 5 && !strncmp(s, "name=", 5)) {
        *name = gretl_strdup(s + 5);
        gretl_unquote(*name, &err);
    } else {
        err = E_PARSE;
    }

    return err;
}

/**
 * equation_system_start:
 * @param: may incude an estimation method spec.
 * @name: the name to be given to system, if any, otherwise NULL
 * or an empty string.
 * @opt: may include OPT_I for iterative estimation (will be
 * ignored if the the estimation method does not support it).
 * @err: location to receive error code.
 *
 * Start compiling an equation system. Either @param must contain
 * an estimation method or the system must be given a name.
 * If a method is given, the system will be estimated as soon as its
 * definition is complete. If a name is given, the system definition
 * is saved on a stack and it can subsequently be estimated via various
 * methods. If both a name and an estimation method are given, the system
 * is both estimated and saved.
 *
 * Returns: pointer to a new equation system, or %NULL on error.
 */

equation_system *equation_system_start (const char *param,
                                        const char *name,
                                        gretlopt opt,
                                        int *err)
{
    equation_system *sys = NULL;
    char *oldname = NULL;
    int method = -1;

#if SYSDEBUG > 1
    fprintf(stderr, "equation_system_start\n");
#endif

    if (param != NULL) {
        *err = parse_sys_start_param(param, &method, &oldname);
        if (*err) {
            return NULL;
        }
    }

    if (name != NULL && *name != '\0') {
        sys = equation_system_new(method, name, err);
    } else if (oldname != NULL) {
        sys = equation_system_new(method, oldname, err);
        free(oldname);
    } else {
        /* treat the system as "anonymous" */
        sys = equation_system_new(method, "$system", err);
        if (!*err) {
            push_anon_system(sys);
        }
    }

    if (sys != NULL) {
        if (opt & OPT_I) {
            sys->flags |= SYSTEM_ITERATE;
        }
        if (opt & OPT_Q) {
            sys->flags |= SYSTEM_QUIET;
        }
    }

#if SYSDEBUG > 1
    fprintf(stderr, "new system '%s' at %p, flags = %d\n", sys->name,
            (void *) sys, sys->flags);
#endif

    return sys;
}

/* determine the degrees of freedom for the unrestricted
   estimation of @sys
*/

static int system_get_dfu (const equation_system *sys)
{
    int dfu = sys->T * sys->neqns; /* total observations */
    int i, pos;

    for (i=0; i<sys->neqns; i++) {
        /* subtract the number of parameters */
        pos = gretl_list_separator_position(sys->lists[i]);
        if (pos > 0) {
            dfu -= pos - 2;
        } else {
            dfu -= sys->lists[i][0] - 1;
        }
    }

    return dfu;
}

static int get_eqn_ref (const equation_system *sys, int j)
{
    int i, pos, nparm = 0;

    for (i=0; i<sys->neqns; i++) {
        pos = gretl_list_separator_position(sys->lists[i]);
        if (pos > 0) {
            nparm += pos - 2;
        } else {
            nparm += sys->lists[i][0] - 1;
        }
        if (nparm > j) {
            return i;
        }
    }

    return -1;
}

/* Determine whether the restrictions on @sys embodied in @R
   reference just one equation, and if so return the degrees
   of freedom for that equation. Otherwise return 0.
*/

static int maybe_get_single_equation_dfu (const equation_system *sys,
                                          const gretl_matrix *R)
{
    int i, j, eq = -1, eq0 = -1;
    int df = 0;

    for (j=0; j<R->cols; j++) {
        for (i=0; i<R->rows; i++) {
            if (gretl_matrix_get(R, i, j) != 0) {
                eq = get_eqn_ref(sys, j);
                if (eq0 < 0) {
                    eq0 = eq;
                } else if (eq != eq0) {
                    /* multiple equations are referenced */
                    return 0;
                }
            }
        }
    }

    if (eq >= 0) {
        df = sys->T - (sys->lists[eq][0] - 1);
    }

    return df;
}

/* Asymptotic F-test, as in William Greene:

   (1/J) * (Rb-q)' * [R*Var(b)*R']^{-1} * (Rb-q)

   or Chi-square version if dfu = 0
*/

int multi_eqn_wald_test (const gretl_matrix *b,
                         const gretl_matrix *V,
                         const gretl_matrix *R,
                         const gretl_matrix *q,
                         int dfu, gretlopt opt,
                         PRN *prn)
{
    gretl_matrix *Rbq, *RvR;
    int Rrows, dfn;
    double test = NADBL;
    int err = 0;

    if (R == NULL || q == NULL || b == NULL || V == NULL) {
        pputs(prn, "Missing matrix in system Wald test!\n");
        return E_DATA;
    }

    Rrows = gretl_matrix_rows(R);
    dfn = gretl_matrix_rows(R);

    Rbq = gretl_matrix_alloc(Rrows, 1);
    RvR = gretl_matrix_alloc(Rrows, Rrows);

    if (Rbq == NULL || RvR == NULL) {
        err = E_ALLOC;
    } else {
        gretl_matrix_multiply(R, b, Rbq);
        gretl_matrix_subtract_from(Rbq, q);
        gretl_matrix_qform(R, GRETL_MOD_NONE, V,
                           RvR, GRETL_MOD_NONE);
        err = gretl_invert_symmetric_matrix(RvR);
    }

    if (!err) {
        test = gretl_scalar_qform(Rbq, RvR, &err);
    }

    if (!err) {
        double pval;

        if (dfu == 0) {
            pval = chisq_cdf_comp(dfn, test);
        } else {
            test /= dfn;
            pval = snedecor_cdf_comp(dfn, dfu, test);
        }

        record_test_result(test, pval);

        if (!(opt & OPT_Q)) {
            if (dfu == 0) {
                pprintf(prn, "%s:\n", _("Wald test for the specified restrictions"));
                pprintf(prn, "  %s(%d) = %g [%.4f]\n", _("Chi-square"),
                        dfn, test, pval);
            } else {
                pprintf(prn, "%s:\n", _("F test for the specified restrictions"));
                pprintf(prn, "  F(%d,%d) = %g [%.4f]\n", dfn, dfu, test, pval);
            }
            pputc(prn, '\n');
        }
    }

    gretl_matrix_free(Rbq);
    gretl_matrix_free(RvR);

    return err;
}

static int get_wald_dfu (const equation_system *sys,
                         const gretl_matrix *R)
{
    int dfu = 0;

    if (sys->method == SYS_METHOD_OLS || sys->method == SYS_METHOD_TSLS) {
        dfu = maybe_get_single_equation_dfu(sys, R);
    }
    if (dfu == 0) {
        dfu = system_get_dfu(sys);
    }

    return dfu;
}

int system_wald_test (const equation_system *sys,
                      const gretl_matrix *R,
                      const gretl_matrix *q,
                      gretlopt opt,
                      PRN *prn)
{
    const gretl_matrix *b = sys->b;
    const gretl_matrix *V = sys->vcv;
    int dfu = 0;

    if (b == NULL || V == NULL) {
        gretl_errmsg_set("restrict: no estimates are available");
        return E_DATA;
    }

    if (sys->method == SYS_METHOD_OLS || sys->method == SYS_METHOD_TSLS) {
        /* use the F-form? */
        dfu = maybe_get_single_equation_dfu(sys, R);
    }

    return multi_eqn_wald_test(b, V, R, q, dfu, opt, prn);
}

static int system_do_LR_test (const equation_system *sys,
                              double llu, gretlopt opt,
                              PRN *prn)
{
    double X2, pval, llr = sys->ll;
    int df = gretl_matrix_rows(sys->R);
    int err = 0;

    if (na(llr) || na(llu) || llr == 0.0 || llu == 0.0 || df <= 0) {
        fputs("bad or missing data in system LR test\n", stderr);
        return E_DATA;
    }

    X2 = 2.0 * (llu - llr);
    pval = chisq_cdf_comp(df, X2);

    if (!(opt & OPT_Q)) {
        pprintf(prn, "%s:\n", _("LR test for the specified restrictions"));
        pprintf(prn, "  %s = %g\n", _("Restricted log-likelihood"), llr);
        pprintf(prn, "  %s = %g\n", _("Unrestricted log-likelihood"), llu);
        pprintf(prn, "  %s(%d) = %g [%.4f]\n", _("Chi-square"), df, X2, pval);
        pputc(prn, '\n');
    }

    return err;
}

static int sys_test_type (equation_system *sys)
{
    int ret = SYS_TEST_NONE;

    if (system_n_restrictions(sys) > 0) {
        if (sys->method == SYS_METHOD_SUR ||
            sys->method == SYS_METHOD_WLS) {
            if (sys->flags & SYSTEM_ITERATE) {
                ret = SYS_TEST_LR;
            } else {
                ret = SYS_TEST_F;
            }
        } else if (sys->method == SYS_METHOD_OLS ||
                   sys->method == SYS_METHOD_TSLS ||
                   sys->method == SYS_METHOD_3SLS) {
            ret = SYS_TEST_F;
        } else if (sys->method == SYS_METHOD_LIML) {
            /* experimental */
            ret = SYS_TEST_F;
        } else if (sys->method == SYS_METHOD_FIML) {
            ret = SYS_TEST_LR;
        }
    }

    return ret;
}

/* Handle the case where sys->b and sys->vcv have been augmented in
   order to test some restrictions: we now have to trim these matrices
   down to their final size. The vector @b contains unrestricted
   coefficient estimates, so its length can be used to gauge the
   true number of coeffs.
*/

static int shrink_b_and_vcv (const gretl_matrix *b,
                             equation_system *sys)
{
    int nc = gretl_vector_get_length(b);
    gretl_matrix *V;
    double x;
    int i, j;

    if (sys->vcv == NULL) {
        return E_DATA;
    }

    if (sys->vcv->rows == nc) {
        /* no-op (shouldn't happen) */
        return 0;
    }

    V = gretl_matrix_alloc(nc, nc);
    if (V == NULL) {
        return E_ALLOC;
    }

    gretl_matrix_reuse(sys->b, nc, 1);

    for (i=0; i<nc; i++) {
        for (j=0; j<nc; j++) {
            x = gretl_matrix_get(sys->vcv, i, j);
            gretl_matrix_set(V, i, j, x);
        }
    }

    gretl_matrix_replace(&sys->vcv, V);

    return 0;
}

static int estimate_with_test (equation_system *sys, DATASET *dset,
                               int stest, int (*system_est)
                               (equation_system *, DATASET *,
                                gretlopt, PRN *),
                               gretlopt opt, PRN *prn)
{
    gretl_matrix *vcv = NULL;
    gretl_matrix *b = NULL;
    double llu = 0.0;
    int err = 0;

    /* estimate the unrestricted system first */
    sys->flags &= ~SYSTEM_RESTRICT;
    err = (*system_est) (sys, dset, opt | OPT_Q, prn);
    sys->flags ^= SYSTEM_RESTRICT;
    if (err) {
        goto bailout;
    }

    /* grab the data from unrestricted estimation */
    b = sys->b;
    sys->b = NULL;
    if (stest == SYS_TEST_LR) {
        llu = sys->ll;
    } else if (stest == SYS_TEST_F) {
        vcv = sys->vcv;
        sys->vcv = NULL;
    }

    /* now estimate the restricted system */
    system_clear_results(sys);
    err = (*system_est) (sys, dset, opt, prn);
    if (!err) {
        if (stest == SYS_TEST_LR) {
            err = system_do_LR_test(sys, llu, opt, prn);
        } else if (stest == SYS_TEST_F) {
            int dfu = get_wald_dfu(sys, sys->R);

            err = multi_eqn_wald_test(b, vcv, sys->R, sys->q,
                                      dfu, opt, prn);
        }
        shrink_b_and_vcv(b, sys);
    }

 bailout:

    gretl_matrix_free(b);
    gretl_matrix_free(vcv);

    return err;
}

static void
adjust_sys_flags_for_method (equation_system *sys, int method)
{
    char oldflags = sys->flags;

    sys->flags = 0;

    /* the robust option is for now OLS-only */
    if (oldflags & SYSTEM_ROBUST) {
        if (sys->method == SYS_METHOD_OLS) {
            sys->flags |= SYSTEM_ROBUST;
        }
    }

    /* the iterate option is only available for WLS, SUR or 3SLS */
    if (oldflags & SYSTEM_ITERATE) {
        if (sys->method == SYS_METHOD_WLS ||
            sys->method == SYS_METHOD_SUR ||
            sys->method == SYS_METHOD_3SLS) {
            sys->flags |= SYSTEM_ITERATE;
        }
    }

    /* by default, apply df correction for single-equation methods */
    if (sys->method == SYS_METHOD_OLS ||
        sys->method == SYS_METHOD_TSLS ||
        sys->method == SYS_METHOD_LIML) {
        if (oldflags & SYSTEM_DFCORR) {
            sys->flags |= SYSTEM_DFCORR;
        }
    }

    /* carry forward the GEOMEAN flag */
    if (oldflags & SYSTEM_VCV_GEOMEAN) {
        sys->flags |= SYSTEM_VCV_GEOMEAN;
    }
}

static int
set_sys_flags_from_opt (equation_system *sys, gretlopt opt)
{
    sys->flags = 0;

    /* the robust option is for now OLS-only */
    if (opt & OPT_R) {
        if (sys->method == SYS_METHOD_OLS) {
            sys->flags |= SYSTEM_ROBUST;
        } else {
            gretl_errmsg_sprintf(_("Invalid option '--%s'"), "robust");
            return E_BADOPT;
        }
    }

    /* the iterate option is available for WLS, SUR or 3SLS */
    if (opt & OPT_I) {
        if (sys->method == SYS_METHOD_WLS ||
            sys->method == SYS_METHOD_SUR ||
            sys->method == SYS_METHOD_3SLS) {
            sys->flags |= SYSTEM_ITERATE;
        }
    }

    /* by default, apply a df correction for single-equation methods */
    if (sys->method == SYS_METHOD_OLS ||
        sys->method == SYS_METHOD_TSLS ||
        sys->method == SYS_METHOD_LIML) {
        if (!(opt & OPT_N)) {
            sys->flags |= SYSTEM_DFCORR;
        }
    }

    if (opt & OPT_M) {
        sys->flags |= SYSTEM_VCV_GEOMEAN;
    }
    if (opt & OPT_Q) {
        sys->flags |= SYSTEM_QUIET;
    }
    if (opt & OPT_S) {
        /* estimating single equation via LIML */
        sys->flags |= SYSTEM_LIML1;
    }

    return 0;
}

/**
 * equation_system_estimate:
 * @sys: pre-defined equation system.
 * @dset: dataset struct.
 * @opt: may include OPT_V for more verbose operation.
 * @prn: printing struct.
 *
 * Estimate a pre-defined equation system and print the results
 * to @prn.
 *
 * Returns: 0 on success, non-zero on error.
 */

int equation_system_estimate (equation_system *sys, DATASET *dset,
                              gretlopt opt, PRN *prn)
{
    int (*system_est) (equation_system *, DATASET *,
                       gretlopt, PRN *);
    int err = 0;

#if SYSDEBUG
    fprintf(stderr, "*** equation_system_estimate\n");
#endif

    gretl_error_clear();

    if (sys->xlist == NULL || sys->biglist == NULL) {
        /* allow for the possibility that we're looking at a
           system restored from a session file
        */
        err = sys_check_lists(sys, dset);
        if (err) {
            return err;
        }
    }

    if (opt == OPT_UNSET) {
        opt = OPT_NONE;
        adjust_sys_flags_for_method(sys, sys->method);
    } else {
        err = set_sys_flags_from_opt(sys, opt);
        if (err) {
            return err;
        }
    }

    /* if restrictions are in place, reset the restricted flag */
    if (sys->R != NULL && sys->q != NULL) {
        sys->flags |= SYSTEM_RESTRICT;
    }

    /* in case we're re-estimating */
    system_clear_results(sys);
    err = sys_rearrange_eqn_lists(sys, dset);

    if (!err) {
        system_est = get_plugin_function("system_estimate");
        if (system_est == NULL) {
            err = 1;
        }
    }

    if (!err) {
        int stest = sys_test_type(sys);

        if (stest == SYS_TEST_NOTIMP) {
            pputs(prn, _("Sorry, command not available for this estimator"));
            pputc(prn, '\n');
            err = 1;
        } else if (stest != SYS_TEST_NONE) {
            err = estimate_with_test(sys, dset, stest,
                                     system_est, opt, prn);
        } else {
            err = (*system_est) (sys, dset, opt, prn);
        }
    }

    if (sys->tslists != NULL) {
	lists_restore_insts(sys);
    }

    if (!err) {
        sys->smpl_t1 = dset->t1;
        sys->smpl_t2 = dset->t2;
        if (!(sys->flags & SYSTEM_LIML1)) {
            set_as_last_model(sys, GRETL_OBJ_SYS);
        }
    }

    return err;
}

int system_adjust_t1t2 (equation_system *sys, const DATASET *dset)
{
    int err;

    if (sys->biglist == NULL) {
        fprintf(stderr, "system_adjust_t1t2: no 'biglist' present!\n");
        return E_DATA;
    }

    sys->t1 = dset->t1;
    sys->t2 = dset->t2;

    err = list_adjust_sample(sys->biglist, &sys->t1, &sys->t2, dset, NULL);
    if (!err) {
        sys->T = sys->t2 - sys->t1 + 1;
    }

    return err;
}

/* construct a list containing all the variables referenced in the
   system */

static int sys_make_biglist (equation_system *sys)
{
    int *biglist;
    const int *slist;
    const identity *ident;
    int bign = 0;
    int i, j, vj, k;

    for (i=0; i<sys->neqns; i++) {
        bign += sys->lists[i][0] - gretl_list_has_separator(sys->lists[i]);
    }

    for (i=0; i<sys->nidents; i++) {
        bign += 1 + sys->idents[i]->n_atoms;
    }

    if (sys->ilist != NULL) {
        bign += sys->ilist[0];
    }

    biglist = gretl_list_new(bign);
    if (biglist == NULL) {
        return E_ALLOC;
    }

    for (i=1; i<=bign; i++) {
        /* invalidate all elements */
        biglist[i] = -1;
    }

    k = 0;

    /* system equations */
    for (i=0; i<sys->neqns; i++) {
        slist = sys->lists[i];
        for (j=1; j<=slist[0]; j++) {
            vj = slist[j];
            if (vj != LISTSEP && !in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    /* system identities */
    for (i=0; i<sys->nidents; i++) {
        ident = sys->idents[i];
        for (j=0; j<=ident->n_atoms; j++) {
            if (j == 0) {
                vj = ident->depvar;
            } else {
                vj = ident->atoms[j-1].varnum;
            }
            if (!in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    if (sys->ilist != NULL) {
        /* system instruments */
        for (j=1; j<=sys->ilist[0]; j++) {
            vj = sys->ilist[j];
            if (!in_gretl_list(biglist, vj)) {
                biglist[++k] = vj;
            }
        }
    }

    biglist[0] = k;

#if SYSDEBUG
    printlist(biglist, "system biglist (all vars)");
#endif

    free(sys->biglist);
    sys->biglist = biglist;

    return 0;
}

static int sys_get_lag_src (const char *vname, const DATASET *dset)
{
    int fd, src = series_index(dset, vname);

    if (src == dset->v && vname != NULL &&
        (fd = gretl_function_depth()) > 0) {
        int i;

        for (i=1; i<dset->v; i++) {
            if (fd == series_get_stack_level(dset, i) &&
                series_is_listarg(dset, i, NULL) &&
                !strcmp(dset->varname[i], vname)) {
                src = i;
                break;
            }
        }
    }

    return src;
}

static int is_tsls_style_instrument (equation_system *sys, int v)
{
    int i, j, pos;

    for (i=0; i<sys->neqns; i++) {
        pos = gretl_list_separator_position(sys->lists[i]);
        if (pos > 0) {
            for (j=pos+1; j<=sys->lists[i][0]; j++) {
                if (sys->lists[i][j] == v) {
                    return 1;
                }
            }
        }
    }

    return 0;
}

/* Here we deal with the case where the user has specified the
   system equations "tsls-style": besides the left-hand side
   variables, we also need to identify as endogenous any
   regressors that do not appear as instruments.
*/

static int tsls_style_augment_ylist (equation_system *sys,
                                     int **pylist)
{
    const int *list;
    int pos, ny = (*pylist)[0];
    int i, j, vj;

    for (i=0; i<sys->neqns; i++) {
        list = sys->lists[i];
        pos = gretl_list_separator_position(list);
        if (pos == 0) {
            continue;
        }
        for (j=2; j<pos; j++) {
            vj = list[j];
            if (vj == 0 || in_gretl_list(*pylist, vj)) {
                continue;
            }
            if (!is_tsls_style_instrument(sys, vj)) {
                if (ny < sys->neqns) {
                    (*pylist)[++ny] = vj;
                    (*pylist)[0] = ny;
                } else {
                    gretl_list_append_term(pylist, vj);
                    if (*pylist == NULL) {
                        return E_ALLOC;
                    }
                    ny++;
                }
            }
        }
    }

    return ny >= sys->neqns ? 0 : E_DATA;
}

/* It's possible to specify instruments per equation, as in tsls,
   with a two-part regression list. But this cannot be mixed with
   any other means of specifying endogeneity/exogeneity of
   variables.
*/

static int check_for_tsls_style_lists (equation_system *sys,
                                       int *err)
{
    int i, tsls_style = 0;

    for (i=0; i<sys->neqns; i++) {
        if (gretl_list_has_separator(sys->lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (tsls_style) {
        if (sys->ylist != NULL || sys->ilist != NULL) {
            gretl_errmsg_set("You can't mix tsls-style lists with 'endog' or 'instr'");
            *err = E_DATA;
        } else if (sys->nidents > 0) {
            gretl_errmsg_set("You can't mix identities with tsls-style lists");
            *err = E_DATA;
        }
    }

    return tsls_style;
}

/* Given a "partial" system, with equations expressed in tsls
   style and with more endogenous variables than equations, shift
   the "excess" endogenous terms into @xplist. They will still
   be instrumented when the equations are estimated, but they
   will be treated "as if" exogenous when it comes to building
   the structural form of the system.
*/

static void tsls_style_shift_vars (equation_system *sys, int *xplist)
{
    int i, vi, nxp = xplist[0];

    for (i=sys->ylist[0]; i>sys->neqns; i--) {
        vi = sys->ylist[i];
        gretl_list_delete_at_pos(sys->ylist, i);
        xplist[++nxp] = vi;
    }

    xplist[0] = nxp;
}

/* For consistency with the way in which per-equation results
   are organized, we should ensure that if a constant is present
   among the exogenous terms in a system it appears as the
   first element of the system's "xlist".
*/

static void sys_xlist_reshuffle_const (int *list)
{
    int i, cpos = 0;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == 0) {
            cpos = i;
            break;
        }
    }

    if (cpos > 1) {
        for (i=cpos; i>1; i--) {
            list[i] = list[i-1];
        }
        list[1] = 0;
    }
}

/* prior to system estimation, get all the required lists of variables
   in order
*/

static int sys_check_lists (equation_system *sys,
                            const DATASET *dset)
{
    const int *slist;
    const char *vname;
    const identity *ident;
    int user_ylist = (sys->ylist != NULL);
    int user_ilist = (sys->ilist != NULL);
    int *ylist = NULL;
    int *xplist = NULL;
    int src, lag, nlhs;
    int tsls_style;
    int i, j, k, vj;
    int err = 0;

#if SYSDEBUG
    fprintf(stderr, "*** sys_check_lists\n");
#endif

    tsls_style = check_for_tsls_style_lists(sys, &err);
    if (err) {
        return err;
    }

    /* start an empty list for predetermined variables */

    sys->plist = gretl_null_list();
    if (sys->plist == NULL) {
        return E_ALLOC;
    }

    /* Compose a minimal (and possibly incomplete) list of endogenous
       variables, including all vars that appear on the left-hand side
       of structural equations or identities.
    */

    nlhs = sys->neqns + sys->nidents;
    ylist = gretl_list_new(nlhs);
    if (ylist == NULL) {
        return E_ALLOC;
    }

    k = 0;

    for (i=0; i<sys->neqns; i++) {
#if SYSDEBUG
        printlist(sys->lists[i], "incoming equation list");
#endif
        slist = sys->lists[i];
        vj = slist[1];
        if (!in_gretl_list(ylist, vj)) {
            ylist[++k] = vj;
        }
    }

    if (tsls_style) {
        ylist[0] = k;
        err = tsls_style_augment_ylist(sys, &ylist);
        if (err) {
            goto bailout;
        }
    } else {
        for (i=0; i<sys->nidents; i++) {
            ident = sys->idents[i];
            vj = ident->depvar;
            if (!in_gretl_list(ylist, vj)) {
                ylist[++k] = vj;
            }
        }
        ylist[0] = k;
    }

#if SYSDEBUG
    printlist(ylist, "system auto ylist");
#endif

    /* If the user gave a list of endogenous vars (in sys->ylist), check
       that it contains all the presumably endogenous variables we found
       above and recorded in ylist; otherwise use the ylist as computed
       above.
    */
    if (user_ylist) {
        for (j=1; j<=ylist[0]; j++) {
            vj = ylist[j];
            if (!in_gretl_list(sys->ylist, vj)) {
                gretl_errmsg_sprintf("%s appears on the left-hand side "
                                     "of an equation but is not marked as endogenous",
                                     dset->varname[vj]);
                err = E_DATA;
                goto bailout;
            }
        }
    } else {
        sys->ylist = ylist;
        ylist = NULL;
    }

    /* Now compose a big list of all the variables in the system
       (referenced in structural equations or identities, or specified
       as instruments).
    */
    err = sys_make_biglist(sys);
    if (err) {
        goto bailout;
    }

    /* If the user gave an endogenous list, check that it doesn't
       contain any variables that are not actually present in the
       system.
    */
    if (user_ylist) {
        for (j=1; j<=sys->ylist[0]; j++) {
            vj = sys->ylist[j];
            if (!in_gretl_list(sys->biglist, vj)) {
                gretl_errmsg_sprintf("%s is marked as endogenous but is "
                                     "not present in the system",
                                     dset->varname[vj]);
                err = E_DATA;
                goto bailout;
            }
        }
    }

    /* By elimination from biglist, using ylist, compose a maximal
       list of exogenous and predetermined variables, "xplist".
    */

    xplist = gretl_list_new(sys->biglist[0]);
    if (xplist == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    k = 0;
    for (j=1; j<=sys->biglist[0]; j++) {
        vj = sys->biglist[j];
        if (!in_gretl_list(sys->ylist, vj)) {
            xplist[++k] = vj;
        }
    }
    xplist[0] = k;

#if SYSDEBUG
    printlist(xplist, "system auto xplist (exog + predet)");
#endif

    /* If the user gave a list of instruments (sys->ilist), check
       that it does not contain anything that is in fact clearly
       endogenous; otherwise use the xplist we computed above.
    */
    if (user_ilist) {
        for (j=1; j<=sys->ilist[0]; j++) {
            vj = sys->ilist[j];
            if (in_gretl_list(sys->ylist, vj)) {
                gretl_errmsg_sprintf("%s is marked as an instrument "
                                     "but is endogenous",
                                     dset->varname[vj]);
                err = E_DATA;
                goto bailout;
            }
        }
    } else {
        sys->ilist = gretl_list_copy(xplist);
        if (sys->ilist == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && sys->ylist[0] != nlhs) {
        /* Note: check added 2009-08-10, modified 2012-04-05,
           moved in front of the following "xplist" block
           2015-10-29.
        */
        if (tsls_style && sys->ylist[0] > nlhs) {
            /* from the pov of the structural form, endogenous regressors
               without an equation should be treated "as if" exogenous?
            */
            tsls_style_shift_vars(sys, xplist);
        } else {
            gretl_errmsg_sprintf("Found %d endogenous variables but %d equations",
                                 sys->ylist[0], nlhs);
            err = E_DATA;
        }
    }

    /* Now narrow down "xplist", removing variables that are in fact
       lags of endogenous variables (and recording them as such): this
       should leave a list of truly exogenous variables.
    */
    for (j=1; j<=xplist[0] && !err; j++) {
        vj = xplist[j];
        lag = series_get_lag(dset, vj);
        if (lag > 0) {
            vname = series_get_parent_name(dset, vj);
            if (vname != NULL) {
                src = sys_get_lag_src(vname, dset);
                if (in_gretl_list(sys->ylist, src)) {
                    err = add_predet_to_sys(sys, dset, vj, src, lag);
                    gretl_list_delete_at_pos(xplist, j--);
                }
            }
        }
    }

    if (!err) {
        sys_xlist_reshuffle_const(xplist);
    }

#if SYSDEBUG
    printlist(xplist, "final system exog list");
#endif

    if (!err) {
        free(sys->xlist);
        sys->xlist = xplist;
        xplist = NULL;
    }

 bailout:

    free(ylist);
    free(xplist);

    return err;
}

static int sys_has_user_name (equation_system *sys)
{
    return (sys->name != NULL && *sys->name != '\0' &&
            !sys_anonymous(sys->name));
}

/**
 * equation_system_finalize:
 * @sys: pre-defined equation system.
 * @dset: dataset struct.
 * @opt: may include OPT_V for verbose operation, OPT_S
 * to permit estimation of a single equation, OPT_R for
 * a robust variance matrix (OLS only), OPT_N for no
 * df correction (OLS and TSLS only).
 * @prn: printing struct.
 *
 * Finalize an equation system, e.g. in response to "end system".
 * If the system has a specified name, we save it on a stack of
 * defined systems.  If it has a specified estimation method,
 * we go ahead and estimate it.  If it has both, we do both.
 *
 * Returns: 0 on success, non-zero on error.  If the system is
 * mal-formed, it is destroyed.
 */

int equation_system_finalize (equation_system *sys, DATASET *dset,
                              gretlopt opt, PRN *prn)
{
    int mineq = 1;
    int err = 0;

#if SYSDEBUG
    fprintf(stderr, "*** equation_system_finalize\n");
#endif

    gretl_error_clear();

    if (sys == NULL) {
        gretl_errmsg_set(_(nosystem));
        return 1;
    } else if (sys->neqns < mineq) {
        gretl_errmsg_set(_(toofew));
        equation_system_destroy(sys);
        return 1;
    } else if (sys->method >= SYS_METHOD_MAX) {
        gretl_errmsg_set(_(badsystem));
        equation_system_destroy(sys);
        return 1;
    }

    err = sys_check_lists(sys, dset);

    if (!err && !(opt & OPT_S) && sys_has_user_name(sys)) {
        /* save the system for subsequent estimation: but note that we
           should not do this if given OPT_S, for single-equation
           LIML
        */
        err = gretl_stack_object_as(sys, GRETL_OBJ_SYS, sys->name);
    }

    if (!err && sys->method >= 0) {
        if (sys->flags & SYSTEM_QUIET) {
            opt |= OPT_Q;
        }
        err = equation_system_estimate(sys, dset, opt, prn);
    }

    return err;
}

/* Implement the "estimate" command. The full form of the
   command line is

   estimate <sysname> method=<method>

   where <sysname> is the name of a previously defined
   equation system and <method> is a recognized estimation
   method.

   However, if the "last model" is an equation system, we
   can accept just

   estimate method=<method>

   where the system is not named but is implicitly the
   last model.

   The method=<method> field is optional, provided that the
   (named or last-model) system already has a specified
   estimation method; so the minimal form of the command
   line is simply "estimate".
*/

int estimate_named_system (const char *sysname, const char *param,
                           DATASET *dset, gretlopt opt, PRN *prn)
{
    equation_system *sys = NULL;
    int err = 0;

#if SYSDEBUG
    fprintf(stderr, "*** estimate_named_system: '%s'\n", sysname);
#endif

    if (sysname != NULL) {
        /* we got a name */
        if (sys_anonymous(sysname)) {
            sys = get_anonymous_equation_system();
        } else {
            sys = get_equation_system_by_name(sysname);
        }
        if (sys == NULL) {
            gretl_errmsg_sprintf(_("'%s': unrecognized name"), sysname);
            err = E_DATA;
        }
    } else {
        /* no name given: try "last model"? */
        GretlObjType type;
        void *ptr;

        ptr = get_last_model(&type);
        if (ptr == NULL || type != GRETL_OBJ_SYS) {
            gretl_errmsg_sprintf(_("%s: no system was specified"), "estimate");
            err = E_DATA;
        } else {
            sys = ptr;
        }
    }

#if 1 /* do we really want this? */
    if (err) {
        /* we haven't found a system to estimate yet */
        sys = get_anonymous_equation_system();
        if (sys != NULL) {
            gretl_error_clear();
            err = 0;
        }
    }
#endif

    if (!err) {
        int method;

        if (param != NULL) {
            method = sys_get_estimator(param);
        } else {
            method = sys->method;
        }

        if (method < 0 || method >= SYS_METHOD_MAX) {
            gretl_errmsg_set("estimate: no valid method was specified");
            err = E_DATA;
        } else {
            sys->method = method;
            err = equation_system_estimate(sys, dset, opt, prn);
        }
    }

#if SYSDEBUG
    fprintf(stderr, "*** estimate_named_system: returning %d\n", err);
#endif

    return err;
}

/* effective list length, allowance made for IVREG-style lists */

static int get_real_list_length (const int *list)
{
    int i, len = list[0];

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            len = i - 1;
            break;
        }
    }

    return len;
}

/* below: used in formulating restriction matrices for
   an equation system */

int system_get_list_length (const equation_system *sys, int i)
{
    if (i >= 0 && i < sys->neqns) {
        return get_real_list_length(sys->lists[i]);
    } else {
        return 0;
    }
}

int system_max_indep_vars (const equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->neqns; i++) {
        nvi = get_real_list_length(sys->lists[i]) - 1;
        if (nvi > nv) nv = nvi;
    }

    return nv;
}

int system_n_indep_vars (const equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->neqns; i++) {
        nvi = get_real_list_length(sys->lists[i]) - 1;
        nv += nvi;
    }

    return nv;
}

const char *system_short_string (const MODEL *pmod)
{
    int i = gretl_model_get_int(pmod, "method");

    return gretl_system_short_strings[i];
}

void equation_system_set_name (equation_system *sys, const char *name)
{
    if (name == sys->name) {
        return;
    }

    if (sys->name == NULL) {
        sys->name = malloc(MAXSAVENAME);
    }

    if (sys->name != NULL) {
        *sys->name = '\0';
        strncat(sys->name, name, MAXSAVENAME - 1);
    }
}

int *compose_ivreg_list (const equation_system *sys, int i)
{
    int *list;
    int j, k1, k2;

    if (i >= sys->neqns) {
        return NULL;
    }

    k1 = sys->lists[i][0];
    k2 = sys->ilist[0];

    list = gretl_list_new(k1 + k2 + 1);
    if (list == NULL) {
        return NULL;
    }

    for (j=1; j<=list[0]; j++) {
        if (j <= k1) {
            list[j] = sys->lists[i][j];
        } else if (j == k1 + 1) {
            list[j] = LISTSEP;
        } else {
            list[j] = sys->ilist[j - (k1 + 1)];
        }
    }

    return list;
}

int system_normality_test (const equation_system *sys,
                           gretlopt opt, PRN *prn)
{
    int err = 0;

    if (sys->E == NULL || sys->S == NULL) {
        err = 1;
    } else {
        err = multivariate_normality_test(sys->E,
                                          sys->S,
                                          opt,
                                          prn);
    }

    return err;
}

double *system_get_resid_series (equation_system *sys, int eqnum,
                                 DATASET *dset, int *err)
{
    double *u = NULL;
    int t;

    if (sys->E == NULL) {
        *err = E_DATA;
        return NULL;
    }

    u = malloc(dset->n * sizeof *u);
    if (u == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    for (t=0; t<dset->n; t++) {
        if (t < sys->t1 || t > sys->t2) {
            u[t] = NADBL;
        } else {
            u[t] = gretl_matrix_get(sys->E, t - sys->t1, eqnum);
        }
    }

    return u;
}

static gchar *system_get_full_string (const equation_system *sys)
{
    const char *lstr = gretl_system_long_strings[sys->method];
    gchar *ret;

    if (sys->flags & SYSTEM_ITERATE) {
	ret = g_strdup_printf(_("iterated %s"), _(lstr));
    } else {
	ret = g_strdup(_(lstr));
    }

    return ret;
}

/* simple accessor functions */

int system_want_df_corr (const equation_system *sys)
{
    return sys->flags & SYSTEM_DFCORR;
}

int system_n_restrictions (const equation_system *sys)
{
    int nr = 0;

    if (sys->R != NULL && (sys->flags & SYSTEM_RESTRICT)) {
        nr = gretl_matrix_rows(sys->R);
    }

    return nr;
}

int *system_get_list (const equation_system *sys, int i)
{
    if (i >= sys->neqns) return NULL;

    return sys->lists[i];
}

int system_get_depvar (const equation_system *sys, int i)
{
    if (i >= sys->neqns) return 0;

    return sys->lists[i][1];
}

int system_get_method (const equation_system *sys)
{
    return sys->method;
}

int *system_get_endog_vars (const equation_system *sys)
{
    return sys->ylist;
}

int *system_get_instr_vars (const equation_system *sys)
{
    return sys->ilist;
}

void system_attach_coeffs (equation_system *sys, gretl_matrix *b)
{
    gretl_matrix_replace(&sys->b, b);
}

void system_attach_vcv (equation_system *sys, gretl_matrix *vcv)
{
    gretl_matrix_replace(&sys->vcv, vcv);
}

void system_attach_sigma (equation_system *sys, gretl_matrix *S)
{
    gretl_matrix_replace(&sys->S, S);
}

void system_attach_uhat (equation_system *sys, gretl_matrix *E)
{
    gretl_matrix_replace(&sys->E, E);
}

MODEL *system_get_model (const equation_system *sys, int i)
{
    if (sys == NULL || sys->models == NULL || i >= sys->neqns) {
        return NULL;
    }

    return sys->models[i];
}

/* for case of applying df correction to cross-equation
   residual covariance matrix */

int system_vcv_geomean (const equation_system *sys)
{
    return sys->flags & SYSTEM_VCV_GEOMEAN;
}

static int sys_eqn_indep_coeffs (const equation_system *sys, int eq)
{
    int nc;

    if (eq >= sys->neqns) {
        nc = 0;
    } else if (sys->R == NULL) {
        nc = sys->lists[eq][0] - 1;
    } else {
        int nr = gretl_matrix_rows(sys->R);
        int cols = gretl_matrix_cols(sys->R);
        int cmax, cmin = 0;
        int single, i, j;

        nc = sys->lists[eq][0] - 1;

        for (i=0; i<eq; i++) {
            cmin += sys->lists[i][0] - 1;
        }

        cmax = cmin + nc;

        for (i=0; i<nr; i++) {
            single = 1;
            for (j=0; j<cols; j++) {
                if (j >= cmin && j < cmax) {
                    if (gretl_matrix_get(sys->R, i, j) != 0.0) {
                        single = 0;
                        break;
                    }
                }
            }
            if (single) {
                nc--;
            }
        }
    }

    return nc;
}

double system_vcv_denom (const equation_system *sys, int i, int j)
{
    double den = sys->T;

    if ((sys->flags & SYSTEM_VCV_GEOMEAN) &&
        i < sys->neqns && j < sys->neqns) {
        int ki = sys_eqn_indep_coeffs(sys, i);

        if (j == i) {
            den = sys->T - ki;
        } else {
            int kj = sys_eqn_indep_coeffs(sys, j);

            den = (sys->T - ki) * (sys->T - kj);
            den = sqrt(den);
        }
    }

    return den;
}

/* for system over-identification test */

int system_get_overid_df (const equation_system *sys)
{
    int gl, i, k = 0;

    gl = sys->neqns * sys->ilist[0];

    for (i=0; i<sys->neqns; i++) {
        k += sys->lists[i][0] - 1;
    }

    return gl - k;
}

/* dealing with identities (FIML, LIML)
   (really also needed for LIML? not clear why...) */

int rhs_var_in_identity (const equation_system *sys, int lhsvar,
                         int rhsvar)
{
    const identity *ident;
    int i, j;

    for (i=0; i<sys->nidents; i++) {
        ident = sys->idents[i];
        if (ident->depvar == lhsvar) {
            for (j=0; j<ident->n_atoms; j++) {
                if (ident->atoms[j].varnum == rhsvar) {
                    return (ident->atoms[j].op == OP_PLUS)? 1 : -1;
                }
            }
        }
    }

    return 0;
}

static void destroy_ident (identity *ident)
{
    free(ident->atoms);
    free(ident);
}

static identity *ident_new (int nv)
{
    identity *ident;

    ident = malloc(sizeof *ident);
    if (ident == NULL) return NULL;

    ident->n_atoms = nv;
    ident->atoms = malloc(nv * sizeof *ident->atoms);
    if (ident->atoms == NULL) {
        free(ident);
        ident = NULL;
    }

    return ident;
}

static identity *
parse_identity (const char *str, DATASET *dset, int *err)
{
    identity *ident;
    char *test;
    const char *p = str;
    char vname[VNAMELEN];
    int lag, gotop = 0;
    int i, v, nv, len;

#if SYSDEBUG
    fprintf(stderr, "parse_identity: input = '%s'\n", str);
#endif

    /* First pass: count the elements and check everything for validity.
       nv is the number of variables appearing on the right-hand side
       of the identity (nv >= 1).
    */

    i = 0;
    while (*p && !*err) {
        p += strspn(p, " ");
        if (i == 0) {
            /* left-hand side variable */
            *err = extract_varname(vname, p, &len);
            if (!*err) {
                v = current_series_index(dset, vname);
                if (v < 0) {
                    *err = E_UNKVAR;
                } else {
                    p += len;
                    p += strspn(p, " ");
                    if (*p != '=') {
                        *err = E_PARSE;
                    } else {
                        p++;
                        i++;
                        gotop = 1;
                    }
                }
            }
        } else if (*p == '+' || *p == '-') {
            gotop = 1;
            p++;
        } else if (i > 1 && !gotop) {
            /* rhs: no + or - given */
            *err = E_PARSE;
        } else {
            /* right-hand size variable (may be lag) */
            *err = extract_varname(vname, p, &len);
            if (!*err && gotop && len == 0) {
                /* dangling operator */
                *err = E_PARSE;
            }
            if (!*err) {
                v = current_series_index(dset, vname);
                if (v < 0) {
                    *err = E_UNKVAR;
                } else {
                    p += len;
                    if (*p == '(') {
                        lag = strtol(p + 1, &test, 10);
                        if (*test != ')') {
                            *err = E_PARSE;
                        } else if (lag >= 0) {
                            *err = E_DATA;
                        } else {
                            v = laggenr(v, -lag, dset);
                            if (v < 0) {
                                *err = E_DATA;
                            } else {
                                p = test + 1;
                            }
                        }
                    }
                    i++;
                    gotop = 0;
                }
            }
        }
    }

    if (gotop) {
        /* trailing operator */
        *err = E_PARSE;
    }

    if (!*err) {
        nv = i - 1;
        if (nv < 1) {
            *err = E_PARSE;
        } else {
            ident = ident_new(nv);
            if (ident == NULL) {
                *err = E_ALLOC;
                return NULL;
            }
        }
    }

    if (*err) {
        return NULL;
    }

    ident->atoms[0].op = OP_PLUS;

    /* Second pass: fill out the identity */

    p = str;
    i = 0;
    while (*p) {
        p += strspn(p, " ");
        if (i == 0) {
            extract_varname(vname, p, &len);
            ident->depvar = series_index(dset, vname);
            p = strchr(p, '=') + 1;
            i++;
        } else if (*p == '+' || *p == '-') {
            ident->atoms[i-1].op = (*p == '+')? OP_PLUS : OP_MINUS;
            p++;
        } else {
            extract_varname(vname, p, &len);
            v = series_index(dset, vname);
            p += len;
            if (*p == '(') {
                lag = strtol(p + 1, &test, 10);
                v = laggenr(v, -lag, dset);
                p = test + 1;
            }
            ident->atoms[i-1].varnum = v;
            i++;
        }
    }

    return ident;
}

/* add information regarding a predetermined regressor */

static int
add_predet_to_sys (equation_system *sys, const DATASET *dset,
                   int id, int src, int lag)
{
    int n = sys->plist[0];
    int *test;
    predet *pre;

    if (id < 0 || src < 0 || id >= dset->v || src >= dset->v) {
        /* something screwy */
        return E_DATA;
    }

    if (in_gretl_list(sys->plist, id)) {
        /* already present */
        return 0;
    }

    pre = realloc(sys->pre_vars, (n + 1) * sizeof *pre);
    if (pre == NULL) {
        return E_ALLOC;
    }

    sys->pre_vars = pre;
    sys->pre_vars[n].id = id;
    sys->pre_vars[n].src = src;
    sys->pre_vars[n].lag = lag;

    /* add to list of predetermined variable IDs */
    test = gretl_list_append_term(&sys->plist, id);
    if (test == NULL) {
        return E_ALLOC;
    }

    return 0;
}

static int
add_identity_to_sys (equation_system *sys, const char *line,
                     DATASET *dset)
{
    identity **pident;
    identity *ident;
    int ni = sys->nidents;
    int err = 0;

    ident = parse_identity(line, dset, &err);
    if (ident == NULL) return err;

    /* connect the identity to the equation system */
    pident = realloc(sys->idents, (ni + 1) * sizeof *sys->idents);
    if (pident == NULL) {
        destroy_ident(ident);
        return E_ALLOC;
    }

    sys->idents = pident;
    sys->idents[ni] = ident;
    sys->nidents += 1;

    return 0;
}

static int
add_aux_list_to_sys (equation_system *sys, const char *line,
                     DATASET *dset, int which)
{
    int *list;
    int err = 0;

    if (which == ENDOG_LIST) {
        if (sys->ylist != NULL) {
            gretl_errmsg_set("Only one list of endogenous variables may be given");
            return 1;
        }
    } else if (which == INSTR_LIST) {
        if (sys->ilist != NULL) {
            gretl_errmsg_set("Only one list of instruments may be given");
            return 1;
        }
    } else {
        return E_DATA;
    }

    line += strspn(line, " ");

    list = generate_list(line, dset, SYSTEM, &err);

    if (!err) {
        if (which == ENDOG_LIST) {
            sys->ylist = list;
        } else {
            sys->ilist = list;
        }
    }

    return err;
}

/**
 * system_parse_line:
 * @sys: initialized equation system.
 * @line: command line.
 * @dset: dataset struct.
 *
 * Modifies @sys according to the command supplied in @line,
 * which must start with "identity" (and supply an identity
 * to be added to the system, or "endog" (and supply a list
 * of endogenous variables), or "instr" (and supply a list of
 * instrumental variables).
 *
 * Returns: 0 on success, non-zero on failure, in which case
 * @sys is destroyed.
 */

int
system_parse_line (equation_system *sys, const char *line,
                   DATASET *dset)
{
    int err = 0;

    gretl_error_clear();

#if SYSDEBUG > 1
    fprintf(stderr, "*** system_parse_line: '%s'\n", line);
#endif

    line += strspn(line, " ");

    if (!strncmp(line, "identity ", 9)) {
        err = add_identity_to_sys(sys, line + 8, dset);
    } else if (!strncmp(line, "endog ", 6)) {
        err = add_aux_list_to_sys(sys, line + 5, dset, ENDOG_LIST);
    } else if (!strncmp(line, "instr ", 6)) {
        err = add_aux_list_to_sys(sys, line + 5, dset, INSTR_LIST);
    } else {
        err = E_PARSE;
    }

    if (err) {
        equation_system_destroy(sys);
    }

    return err;
}

void
system_set_restriction_matrices (equation_system *sys,
                                 gretl_matrix *R, gretl_matrix *q)
{
    system_clear_restrictions(sys);

    sys->R = R;
    sys->q = q;

    sys->flags |= SYSTEM_RESTRICT;
}

double *
equation_system_get_series (const equation_system *sys,
                            const DATASET *dset,
                            int idx, const char *key, int *err)
{
    const gretl_matrix *M = NULL;
    double *x = NULL;
    const char *msel;
    int t, col = 0;

    if (sys == NULL || (idx != M_UHAT && idx != M_YHAT)) {
        *err = E_BADSTAT;
        return NULL;
    }

    msel = strchr(key, '[');
    if (msel == NULL || sscanf(msel, "[,%d]", &col) != 1) {
        *err = E_PARSE;
    } else if (col <= 0 || col > sys->neqns) {
        *err = E_DATA;
    }

    if (!*err) {
        M = (idx == M_UHAT)? sys->E : sys->yhat;
        if (M == NULL) {
            *err = E_DATA;
        }
    }

    if (!*err) {
        x = malloc(dset->n * sizeof *x);
        if (x == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        int s = 0;

        col--; /* switch to 0-based */
        for (t=0; t<dset->n; t++) {
            if (t < sys->t1 || t > sys->t2) {
                x[t] = NADBL;
            } else {
                x[t] = gretl_matrix_get(M, s++, col);
            }
        }
    }

    return x;
}

static int system_add_yhat_matrix (equation_system *sys)
{
    double x, avc;
    int k = 0;
    int i, s, t;

    sys->yhat = gretl_matrix_alloc(sys->T, sys->neqns);
    if (sys->yhat == NULL) {
        return E_ALLOC;
    }

    gretl_matrix_set_t1(sys->yhat, sys->t1);
    gretl_matrix_set_t2(sys->yhat, sys->t2);

    for (i=0; i<sys->neqns; i++) {
        s = 0;
        for (t=sys->t1; t<=sys->t2; t++) {
            x = sys->models[i]->yhat[t];
            gretl_matrix_set(sys->yhat, s++, i, x);
        }
        k += sys->models[i]->ncoeff;
    }

    k -= system_n_restrictions(sys);

    avc = k / (double) sys->neqns;
    sys->df = sys->T - floor(avc);

    return 0;
}

static gretl_matrix *get_stderr_vec (const gretl_matrix *V)
{
    gretl_matrix *se;
    int k = V->rows;

    se = gretl_column_vector_alloc(k);

    if (se != NULL) {
        double x;
        int i;

        for (i=0; i<k; i++) {
            x = gretl_matrix_get(V, i, i);
            se->val[i] = sqrt(x);
        }
    }

    return se;
}

/* retrieve a copy of a specified matrix from an equation
   system */

gretl_matrix *
equation_system_get_matrix (const equation_system *sys, int idx,
                            int *err)
{
    gretl_matrix *M = NULL;

    if (sys == NULL) {
        *err = E_BADSTAT;
        return NULL;
    }

    switch (idx) {
    case M_COEFF:
        if (sys->b == NULL) {
            *err = E_BADSTAT;
        } else {
            M = gretl_matrix_copy(sys->b);
        }
        break;
    case M_UHAT:
        M = gretl_matrix_copy(sys->E);
        break;
    case M_YHAT:
        M = gretl_matrix_copy(sys->yhat);
        break;
    case M_VCV:
    case M_SE:
        if (sys->vcv == NULL) {
            *err = E_BADSTAT;
        } else if (idx == M_SE) {
            M = get_stderr_vec(sys->vcv);
        } else {
            M = gretl_matrix_copy(sys->vcv);
        }
        break;
    case M_SIGMA:
        M = gretl_matrix_copy(sys->S);
        break;
    case M_SYSGAM:
        if (sys->Gamma == NULL) {
            *err = E_BADSTAT;
        } else {
            M = gretl_matrix_copy(sys->Gamma);
        }
        break;
    case M_SYSA:
        if (sys->A == NULL) {
            *err = E_BADSTAT;
        } else {
            M = gretl_matrix_copy(sys->A);
        }
        break;
    case M_SYSB:
        if (sys->B == NULL) {
            *err = E_BADSTAT;
        } else {
            M = gretl_matrix_copy(sys->B);
        }
        break;
    default:
        *err = E_BADSTAT;
        break;
    }

    if (M == NULL && !*err) {
        *err = E_ALLOC;
    }

    return M;
}

int highest_numbered_var_in_system (const equation_system *sys,
                                    const DATASET *dset)
{
    int i, j, v, vmax = 0;

    if (sys->biglist != NULL) {
        for (j=1; j<=sys->biglist[0]; j++) {
            v = sys->biglist[j];
            if (v > vmax) {
                vmax = v;
            }
        }
    } else {
        /* should not happen */
        for (i=0; i<sys->neqns; i++) {
            for (j=1; j<=sys->lists[i][0]; j++) {
                v = sys->lists[i][j];
                if (v == LISTSEP || v >= dset->v) {
                    /* temporary variables, already gone? */
                    continue;
                }
                if (v > vmax) {
                    vmax = v;
                }
            }
        }
    }

    return vmax;
}

static identity *sys_retrieve_identity (xmlNodePtr node, int *err)
{
    identity *ident;
    xmlNodePtr cur;
    int n_atoms, depvar;
    int i, got = 0;

    got += gretl_xml_get_prop_as_int(node, "n_atoms", &n_atoms);
    got += gretl_xml_get_prop_as_int(node, "depvar", &depvar);
    if (got < 2) {
        *err = E_DATA;
        return NULL;
    }

    ident = ident_new(n_atoms);
    if (ident == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    ident->depvar = depvar;

    cur = node->xmlChildrenNode;

    i = 0;
    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "id_atom")) {
            got = gretl_xml_get_prop_as_int(cur, "op", &ident->atoms[i].op);
            got += gretl_xml_get_prop_as_int(cur, "varnum", &ident->atoms[i].varnum);
            if (got < 2) {
                *err = E_DATA;
            } else {
                i++;
            }
        }
        cur = cur->next;
    }

    if (!*err && i != n_atoms) {
        *err = E_DATA;
    }

    if (*err) {
        destroy_ident(ident);
        ident = NULL;
    }

    return ident;
}

equation_system *
equation_system_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    equation_system *sys;
    xmlNodePtr cur;
    char *sname;
    int method = 0;
    int i, j, got = 0;

    got += gretl_xml_get_prop_as_string(node, "name", &sname);
    got += gretl_xml_get_prop_as_int(node, "method", &method);

    if (got < 2) {
        *err = E_DATA;
        return NULL;
    }

    sys = equation_system_new(method, sname, err);
    if (*err) {
        return NULL;
    }

    got = 0;
    got += gretl_xml_get_prop_as_int(node, "n_equations", &sys->neqns);
    got += gretl_xml_get_prop_as_int(node, "nidents", &sys->nidents);
    got += gretl_xml_get_prop_as_int(node, "flags", &sys->flags);

    if (got < 3) {
        *err = E_DATA;
        goto bailout;
    }

    gretl_xml_get_prop_as_int(node, "order", &sys->order);

    sys->lists = malloc(sys->neqns * sizeof sys->lists);
    if (sys->lists == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    if (sys->nidents > 0) {
        sys->idents = malloc(sys->nidents * sizeof sys->idents);
        if (sys->idents == NULL) {
            *err = E_ALLOC;
            goto bailout;
        }
    }

    cur = node->xmlChildrenNode;

    i = j = 0;
    while (cur != NULL && !*err) {
        if (!xmlStrcmp(cur->name, (XUC) "eqnlist")) {
            sys->lists[i++] = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "endog_vars")) {
            sys->ylist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "instr_vars")) {
            sys->ilist = gretl_xml_get_list(cur, doc, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "identity")) {
            sys->idents[j++] = sys_retrieve_identity(cur, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "R")) {
            sys->R = gretl_xml_get_matrix(cur, doc, err);
        } else if (!xmlStrcmp(cur->name, (XUC) "q")) {
            sys->q = gretl_xml_get_matrix(cur, doc, err);
        }
        cur = cur->next;
    }

    if (!*err && (i != sys->neqns || j != sys->nidents)) {
        *err = E_DATA;
    }

    if (*err) {
        equation_system_destroy(sys);
        sys = NULL;
    }

 bailout:

    return sys;
}

static void xml_print_identity (identity *ident, PRN *prn)
{
    int i;

    pprintf(prn, "<identity n_atoms=\"%d\" depvar=\"%d\">\n",
            ident->n_atoms, ident->depvar);
    for (i=0; i<ident->n_atoms; i++) {
        pprintf(prn, " <id_atom op=\"%d\" varnum=\"%d\"/>\n",
                ident->atoms[i].op, ident->atoms[i].varnum);
    }
    pputs(prn, "</identity>\n");
}

int equation_system_serialize (equation_system *sys,
                               SavedObjectFlags flags,
                               PRN *prn)
{
    char *xmlname = NULL;
    int tsls_style = 0;
    int i, err = 0;

    if (sys->name == NULL || *sys->name == '\0') {
        xmlname = gretl_strdup("none");
    } else {
        xmlname = gretl_xml_encode(sys->name);
    }

    pprintf(prn, "<gretl-equation-system name=\"%s\" saveflags=\"%d\" method=\"%d\" ",
            xmlname, flags, sys->method);
    free(xmlname);

    pprintf(prn, "n_equations=\"%d\" nidents=\"%d\" flags=\"%d\" order=\"%d\">\n",
            sys->neqns, sys->nidents, sys->flags, sys->order);

    for (i=0; i<sys->neqns; i++) {
        gretl_xml_put_tagged_list("eqnlist", sys->lists[i], prn);
    }

    for (i=0; i<sys->neqns; i++) {
        if (gretl_list_has_separator(sys->lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (!tsls_style) {
        gretl_xml_put_tagged_list("endog_vars", sys->ylist, prn);
        gretl_xml_put_tagged_list("instr_vars", sys->ilist, prn);
    }

    for (i=0; i<sys->nidents; i++) {
        xml_print_identity(sys->idents[i], prn);
    }

    gretl_matrix_serialize(sys->R, "R", prn);
    gretl_matrix_serialize(sys->q, "q", prn);

    pputs(prn, "</gretl-equation-system>\n");

    return err;
}

static gretl_bundle *bundlize_diag_test (const equation_system *sys)
{
    gretl_bundle *b = NULL;
    double test = NADBL;
    double pval = NADBL;
    int k = sys->S->rows;
    int df = k * (k - 1) / 2;

    if (sys->method == SYS_METHOD_SUR && sys->iters > 0) {
        if (!na(sys->ldet) && sys->diag_test > 0.0) {
            test = sys->T * (sys->diag_test - sys->ldet);
            pval = chisq_cdf_comp(df, test);

        }
    } else {
        test = sys->diag_test;
        pval = chisq_cdf_comp(df, test);
    }

    if (!na(test) && !na(pval)) {
        b = gretl_bundle_new();
        if (b != NULL) {
            gretl_bundle_set_scalar(b, "test", test);
            gretl_bundle_set_scalar(b, "pvalue", pval);
            gretl_bundle_set_int(b, "df", df);
        }
    }

    return b;
}

int equation_system_bundlize (equation_system *sys,
                              gretl_bundle *b)
{
    const char *s;
    char lname[16];
    int tsls_style = 0;
    int i, err = 0;

    gretl_bundle_set_string(b, "command", "system");

    gretl_bundle_set_int(b, "neqns", sys->neqns);
    s = system_method_short_string(sys->method);
    if (s != NULL) {
        gretl_bundle_set_string(b, "method", s);
    }
    if (sys->flags & SYSTEM_ROBUST) {
        gretl_bundle_set_int(b, "robust", 1);
    }

    /* convert to 1-based for user-space */
    gretl_bundle_set_int(b, "t1", sys->t1 + 1);
    gretl_bundle_set_int(b, "t2", sys->t2 + 1);
    gretl_bundle_set_int(b, "sample_t1", sys->smpl_t1 + 1);
    gretl_bundle_set_int(b, "sample_t2", sys->smpl_t2 + 1);

    gretl_bundle_set_matrix(b, "coeff", sys->b);
    gretl_bundle_set_matrix(b, "vcv",   sys->vcv);
    gretl_bundle_set_matrix(b, "sigma", sys->S);
    gretl_bundle_set_matrix(b, "uhat",  sys->E);
    gretl_bundle_set_matrix(b, "yhat",  sys->yhat);

    if (sys->Gamma != NULL) {
        gretl_bundle_set_matrix(b, "Gamma", sys->Gamma);
    }
    if (sys->A != NULL) {
        gretl_bundle_set_matrix(b, "A", sys->A);
    }
    if (sys->B != NULL) {
        gretl_bundle_set_matrix(b, "B", sys->B);
    }

    /* per-equation model structs? */

    for (i=0; i<sys->neqns; i++) {
        sprintf(lname, "list%d", i+1);
        gretl_bundle_set_list(b, lname, sys->lists[i]);
    }

    for (i=0; i<sys->neqns; i++) {
        if (gretl_list_has_separator(sys->lists[i])) {
            tsls_style = 1;
            break;
        }
    }

    if (!tsls_style) {
        gretl_bundle_set_list(b, "endog_vars", sys->ylist);
        gretl_bundle_set_list(b, "instr_vars", sys->ilist);
    }

    if (sys->xlist != NULL) {
	gretl_bundle_set_list(b, "xlist", sys->xlist);
    }

    if (sys->plist != NULL) {
        gretl_bundle_set_list(b, "predet_vars", sys->plist);
    }

#if 0 /* FIXME */
    for (i=0; i<sys->nidents; i++) {
        xml_print_identity(sys->idents[i], fp);
    }
#endif

    /* restrictions? */
    if (sys->R != NULL) {
        gretl_bundle_set_matrix(b, "R", sys->R);
    }
    if (sys->q != NULL) {
        gretl_bundle_set_matrix(b, "q", sys->q);
    }

    if (sys->diag_test > 0) {
        gretl_bundle *dt = bundlize_diag_test(sys);

        if (dt != NULL) {
            gretl_bundle_donate_data(b, "diag_test", dt,
                                     GRETL_TYPE_BUNDLE, 0);
        }
    }

    return err;
}

static int sur_ols_diag (equation_system *sys)
{
    double s2, ls2sum = 0.0;
    int i, err = 0;

    for (i=0; i<sys->neqns; i++) {
        s2 = gretl_model_get_double(sys->models[i], "ols_sigma_squared");
        if (na(s2)) {
            err = 1;
            break;
        }
        ls2sum += log(s2);
    }

    if (!err) {
        sys->diag_test = ls2sum;
    }

    return err;
}

static void
print_system_overid_test (const equation_system *sys, PRN *prn)
{
    int tex = tex_format(prn);
    int df = system_get_overid_df(sys);
    double pv;

    if (sys->method == SYS_METHOD_FIML && df > 0) {
        double X2;

        if (na(sys->ll) || na(sys->llu) ||
            sys->ll == 0.0 || sys->llu == 0.0) {
            return;
        }

        X2 = 2.0 * (sys->llu - sys->ll);
        pv = chisq_cdf_comp(df, X2);

        if (tex) {
            pprintf(prn, "%s:\\\\\n", _("LR over-identification test"));
            if (sys->ll < 0) {
                pprintf(prn, "  %s = $-$%g", _("Restricted log-likelihood"), -sys->ll);
            } else {
                pprintf(prn, "  %s = %g", _("Restricted log-likelihood"), sys->ll);
            }
            gretl_prn_newline(prn);
            if (sys->llu < 0) {
                pprintf(prn, "  %s = $-$%g", _("Unrestricted log-likelihood"), -sys->llu);
            } else {
                pprintf(prn, "  %s = %g", _("Unrestricted log-likelihood"), sys->llu);
            }
            gretl_prn_newline(prn);
            pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]\n", df, X2, pv);
        } else {
            pprintf(prn, "%s:\n", _("LR over-identification test"));
            pprintf(prn, "  %s = %g\n", _("Restricted log-likelihood"), sys->ll);
            pprintf(prn, "  %s = %g\n", _("Unrestricted log-likelihood"), sys->llu);
            pprintf(prn, "  %s(%d) = %g [%.4f]\n\n", _("Chi-square"), df, X2, pv);
        }
    } else if ((sys->method == SYS_METHOD_3SLS ||
                sys->method == SYS_METHOD_SUR) && df > 0) {
        if (na(sys->X2) || sys->X2 <= 0.0) {
            if (!tex) {
                pputs(prn, _("Warning: the Hansen-Sargan over-identification test "
                             "failed.\nThis probably indicates that the estimation "
                             "problem is ill-conditioned.\n"));
                pputc(prn, '\n');
            }
            return;
        }

        pv = chisq_cdf_comp(df, sys->X2);

        if (tex) {
            pprintf(prn, "\\noindent %s:\\\\\n",
                    _("Hansen--Sargan over-identification test"));
            pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]\\\\\n", df, sys->X2, pv);
        } else {
            pprintf(prn, "%s:\n", _("Hansen-Sargan over-identification test"));
            pprintf(prn, "  %s(%d) = %g [%.4f]\n\n", _("Chi-square"),
                    df, sys->X2, pv);
        }
    }
}

int system_diag_test (const equation_system *sys, double *test,
                      double *pval)
{
    int k, df, err = 0;

    if (sys->S == NULL) {
        return E_BADSTAT;
    }

    k = sys->S->rows;
    df = k * (k - 1) / 2;

    if (sys->method == SYS_METHOD_SUR && sys->iters > 0) {
        /* iterated SUR */
        if (!na(sys->ldet) && sys->diag_test != 0.0) {
            double lr = sys->T * (sys->diag_test - sys->ldet);

            if (test != NULL) {
                *test = lr;
            }
            if (pval != NULL) {
                *pval = chisq_cdf_comp(df, lr);
            }
        } else {
            err = E_BADSTAT;
        }
    } else if (sys->diag_test > 0) {
        /* other estimators */
        if (test != NULL) {
            *test = sys->diag_test;
        }
        if (pval != NULL) {
            *pval = chisq_cdf_comp(df, sys->diag_test);
        }
    } else {
        err = E_BADSTAT;
    }

    return err;
}

int system_print_sigma (const equation_system *sys, PRN *prn)
{
    double test, pval;
    int tex = tex_format(prn);
    int k, df;

    if (sys->S == NULL) {
        return E_DATA;
    }

    k = sys->S->rows;
    df = k * (k - 1) / 2;

    print_contemp_covariance_matrix(sys->S, sys->ldet, prn);

    if (sys->method == SYS_METHOD_SUR && sys->iters > 0) {
        if (!na(sys->ldet) && sys->diag_test != 0.0) {
            test = sys->T * (sys->diag_test - sys->ldet);
            pval = chisq_cdf_comp(df, test);
            if (tex) {
                pprintf(prn, "%s:\\\\\n", _("LR test for diagonal covariance matrix"));
                pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]", df, test, pval);
                gretl_prn_newline(prn);
            } else {
                pprintf(prn, "%s:\n", _("LR test for diagonal covariance matrix"));
                pprintf(prn, "  %s(%d) = %g [%.4f]\n", _("Chi-square"),
                        df, test, pval);
            }
        }
    } else {
        const char *labels[] = {
            N_("Breusch--Pagan test for diagonal covariance matrix"),
            N_("Robust Breusch--Pagan test for diagonal covariance matrix"),
            N_("Breusch-Pagan test for diagonal covariance matrix"),
            N_("Robust Breusch-Pagan test for diagonal covariance matrix"),
        };
        const char *label;

        test = sys->diag_test;
        if (sys->flags & SYSTEM_ROBUST) {
            label = tex ? labels[1] : labels[3];
        } else {
            label = tex ? labels[0] : labels[2];
        }
        if (test > 0) {
            pval = chisq_cdf_comp(df, test);
            if (tex) {
                pprintf(prn, "%s:", _(label));
                gretl_prn_newline(prn);
                pprintf(prn, "  $\\chi^2(%d)$ = %g [%.4f]", df, test, pval);
                gretl_prn_newline(prn);
            } else {
                pprintf(prn, "%s:\n", _(label));
                pprintf(prn, "  %s(%d) = %g [%.4f]\n", _("Chi-square"),
                        df, test, pval);
            }
        }
    }

    pputc(prn, '\n');

    return 0;
}

enum {
    ENDOG,
    EXOG,
    PREDET
};

static int categorize_variable (int vnum, const equation_system *sys,
                                int *col, int *lag)
{
    int pos, ret = -1;

    *lag = 0;

    pos = in_gretl_list(sys->ylist, vnum);
    if (pos > 0) {
        ret = ENDOG;
    } else {
        pos = in_gretl_list(sys->xlist, vnum);
        if (pos > 0) {
            ret = EXOG;
        } else {
            int v = get_predet_parent(sys, vnum, lag);

            pos = in_gretl_list(sys->ylist, v);
            if (pos > 0) {
                ret = PREDET;
            }
        }
    }

    if (pos > 0) {
        *col = pos - 1;
    }

    return ret;
}

static int get_col_and_lag (int vnum, const equation_system *sys,
                            int *col, int *lag)
{
    int vp, pos, err = 0;

    vp = get_predet_parent(sys, vnum, lag);
    pos = in_gretl_list(sys->ylist, vp);
    if (pos > 0) {
        *col = pos - 1;
    } else {
        *col = 0;
        err = 1;
    }

    return err;
}

/* reduced-form error covariance matrix */

static int
sys_add_RF_covariance_matrix (equation_system *sys, int n)
{
    gretl_matrix *G = NULL, *S = NULL;
    int err = 0;

    if (!gretl_is_identity_matrix(sys->Gamma)) {
        G = gretl_matrix_copy(sys->Gamma);
        if (G == NULL) {
            return E_ALLOC;
        }
    }

    if (n > sys->S->rows) {
        S = gretl_zero_matrix_new(n, n);
        if (S == NULL) {
            gretl_matrix_free(G);
            return E_ALLOC;
        }
        gretl_matrix_inscribe_matrix(S, sys->S, 0, 0, GRETL_MOD_NONE);
    }

    sys->Sr = gretl_matrix_alloc(n, n);
    if (sys->Sr == NULL) {
        gretl_matrix_free(G);
        gretl_matrix_free(S);
        return E_ALLOC;
    }

    if (G != NULL) {
        err = gretl_SVD_invert_matrix(G);
        if (!err) {
            err = gretl_matrix_qform(G, GRETL_MOD_NONE,
                                     (S != NULL)? S : sys->S,
                                     sys->Sr, GRETL_MOD_NONE);
        }
        gretl_matrix_free(G);
    } else {
        gretl_matrix_copy_values(sys->Sr, (S != NULL)? S : sys->S);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->Sr, "G^{-1} S G^{-1}'");
#endif

    if (S != NULL) {
        gretl_matrix_free(S);
    }

    return err;
}

static int sys_add_structural_form (equation_system *sys,
                                    const DATASET *dset)
{
    const int *ylist = sys->ylist;
    const int *xlist = sys->xlist;
    int ne = sys->neqns;
    int ni = sys->nidents;
    int n = ne + ni;
    double x = 0.0;
    int type, col = 0;
    int i, j, vj, lag;
    int err = 0;

#if SYSDEBUG
    printlist(ylist, "endogenous vars");
    printlist(xlist, "exogenous vars");
#endif

    if (ylist[0] > n) {
        gretl_errmsg_set("system: can't add structural form");
        return E_DATA;
    }

    sys->order = sys_max_predet_lag(sys);

    /* allocate coefficient matrices */

    sys->Gamma = gretl_zero_matrix_new(n, n);
    if (sys->Gamma == NULL) {
        return E_ALLOC;
    }

    if (sys->order > 0) {
        sys->A = gretl_zero_matrix_new(n, sys->order * n);
        if (sys->A == NULL) {
            return E_ALLOC;
        }
    }

    if (xlist[0] > 0) {
        sys->B = gretl_zero_matrix_new(n, xlist[0]);
        if (sys->B == NULL) {
            return E_ALLOC;
        }
    }

    /* process stochastic equations */

    for (i=0; i<sys->neqns && !err; i++) {
        const MODEL *pmod = sys->models[i];
        const int *mlist = pmod->list;

        for (j=1; j<=mlist[0] && mlist[j]!=LISTSEP; j++) {
            vj = mlist[j];
            type = categorize_variable(vj, sys, &col, &lag);
            x = (j > 1)? pmod->coeff[j-2] : 1.0;
            if (type == ENDOG) {
                if (j == 1) {
                    /* left-hand side variable */
                    gretl_matrix_set(sys->Gamma, i, col, 1.0);
                } else {
                    gretl_matrix_set(sys->Gamma, i, col, -x);
                }
            } else if (type == EXOG) {
                gretl_matrix_set(sys->B, i, col, x);
            } else if (type == PREDET) {
                col += n * (lag - 1);
                gretl_matrix_set(sys->A, i, col, x);
            } else {
                gretl_errmsg_sprintf("structural form:\n couldn't categorize series "
                                     "%s in eqn %d\n", dset->varname[vj], i+1);
                err = E_DATA;
            }
        }
    }

    /* process identities */

    for (i=0; i<sys->nidents && !err; i++) {
        const identity *ident = sys->idents[i];

        for (j=0; j<=ident->n_atoms; j++) {
            if (j == 0) {
                vj = ident->depvar;
                x = 1.0;
            } else {
                vj = ident->atoms[j-1].varnum;
                x = (ident->atoms[j-1].op)? -1.0 : 1.0;
            }
            type = categorize_variable(vj, sys, &col, &lag);
            if (type == ENDOG) {
                if (j == 0) {
                    gretl_matrix_set(sys->Gamma, ne+i, col, 1.0);
                } else {
                    gretl_matrix_set(sys->Gamma, ne+i, col, -x);
                }
            } else if (type == EXOG) {
                gretl_matrix_set(sys->B, ne+i, col, x);
            } else if (type == PREDET) {
                col += n * (lag - 1);
                gretl_matrix_set(sys->A, ne+i, col, x);
            } else {
                gretl_errmsg_sprintf("structural form:\n couldn't categorize series "
                                     "%s in identity %d\n", dset->varname[vj], i+1);
                err = E_DATA;
            }
        }
    }

    if (!err) {
        err = sys_add_RF_covariance_matrix(sys, n);
        if (err) {
            fprintf(stderr, "error %d in sys_add_RF_covariance_matrix\n", err);
        }
    }

#if SYSDEBUG
    gretl_matrix_print(sys->Gamma, "sys->Gamma");

    if (sys->A != NULL) {
        gretl_matrix_print(sys->A, "sys->A");
    } else {
        fputs("sys->A: no lagged endog variables used as instruments\n\n", stderr);
    }

    if (sys->B != NULL) {
        gretl_matrix_print(sys->B, "sys->B");
    } else {
        fputs("sys->B: no truly exogenous variables present\n", stderr);
    }
#endif

    return err;
}

static gretl_matrix *sys_companion_matrix (equation_system *sys,
                                           int *err)
{
    int m = sys->A->rows;
    int n = sys->A->cols;
    gretl_matrix *C;

    if (m == n) {
        C = sys->A;
    } else {
        C = gretl_zero_matrix_new(n, n);
        if (C == NULL) {
            *err = E_ALLOC;
        } else {
            gretl_matrix_inscribe_matrix(C, sys->A, 0, 0,
                                         GRETL_MOD_NONE);
            gretl_matrix_inscribe_I(C, m, 0, n - m);
        }
    }

#if SYSDEBUG
    gretl_matrix_print(C, "system companion matrix");
#endif

    return C;
}

static gretl_matrix *sys_get_fcast_se (equation_system *sys,
                                       int periods, int *err)
{
    int n = sys->neqns + sys->nidents;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *V0 = NULL, *Vt = NULL;
    gretl_matrix *C = NULL, *se = NULL;
    double vti;
    int i, k, t;

    if (periods <= 0) {
        fprintf(stderr, "Invalid number of periods\n");
        *err = E_DATA;
        return NULL;
    }

    if (sys->A == NULL) {
        fprintf(stderr, "sys->A is NULL\n");
        *err = E_DATA;
        return NULL;
    }

    C = sys_companion_matrix(sys, err);
    if (*err) {
        return NULL;
    }

    se = gretl_zero_matrix_new(periods, n);
    if (se == NULL) {
        *err = E_ALLOC;
        if (C != sys->A) {
            gretl_matrix_free(C);
        }
        return NULL;
    }

    k = sys->A->cols;

    Vt = gretl_matrix_alloc(k, k);
    V0 = gretl_zero_matrix_new(k, k);
    Tmp = gretl_matrix_alloc(k, k);

    if (Vt == NULL || V0 == NULL || Tmp == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    for (t=0; t<periods; t++) {
        if (t == 0) {
            /* initial variance */
            gretl_matrix_inscribe_matrix(V0, sys->Sr, 0, 0, GRETL_MOD_NONE);
            gretl_matrix_copy_values(Vt, V0);
        } else {
            /* calculate further variances */
            gretl_matrix_copy_values(Tmp, Vt);
            gretl_matrix_qform(C, GRETL_MOD_NONE,
                               Tmp, Vt, GRETL_MOD_NONE);
            gretl_matrix_add_to(Vt, V0);
        }

        for (i=0; i<n; i++) {
            vti = gretl_matrix_get(Vt, i, i);
            gretl_matrix_set(se, t, i, sqrt(vti));
        }
    }

 bailout:

    gretl_matrix_free(V0);
    gretl_matrix_free(Vt);
    gretl_matrix_free(Tmp);

    if (C != sys->A) {
        gretl_matrix_free(C);
    }

    if (*err) {
        gretl_matrix_free(se);
        se = NULL;
    }

    return se;
}

static int sys_add_fcast_variance (equation_system *sys, gretl_matrix *F,
                                   int n_static)
{
    gretl_matrix *se = NULL;
    double ftj, vti;
    int n = sys->neqns + sys->nidents;
    int k = F->rows - n_static;
    int i, j, s;
    int err = 0;

    if (k > 0) {
        se = sys_get_fcast_se(sys, k, &err);
        if (err) {
            goto bailout;
        }
    }

    for (j=0; j<n; j++) {
        for (s=0; s<F->rows; s++) {
            ftj = gretl_matrix_get(F, s, j);
            if (na(ftj)) {
                gretl_matrix_set(F, s, n + j, NADBL);
            } else {
                i = s - n_static;
                if (i < 0) {
                    if (j < sys->neqns) {
                        vti = sqrt(gretl_matrix_get(sys->Sr, j, j));
                    } else {
                        /* LHS of identity */
                        vti = 0.0;
                    }
                } else {
                    vti = gretl_matrix_get(se, i, j);
                }
                gretl_matrix_set(F, s, n + j, vti);
            }
        }
    }

    if (se != NULL) {
        gretl_matrix_free(se);
    }

 bailout:

    if (err) {
        for (i=0; i<n; i++) {
            for (s=0; s<F->rows; s++) {
                gretl_matrix_set(F, s, n + i, NADBL);
            }
        }
    }

    return err;
}

/* Experimental: generate "fitted values " (in-sample, prior to the
   forecast period) for variables that do not appear on the left-hand
   side of any stochastic equation.  But I'm not quite sure how these
   should be defined.
*/

gretl_matrix *sys_get_fitted_values (equation_system *sys,
                                     int v, int t1, int t2,
                                     const DATASET *dset,
                                     int *err)
{
    gretl_matrix *F = NULL;
    gretl_matrix *G = NULL;
    gretl_matrix *yl = NULL, *x = NULL;
    gretl_matrix *y = NULL, *yh = NULL;
    const int *ylist = sys->ylist;
    const int *xlist = sys->xlist;
    const int *plist = sys->plist;
    int n = sys->neqns + sys->nidents;
    double xit;
    int col, T;
    int i, vi, s, t, lag;

    if (sys->Gamma == NULL) {
        *err = E_DATA;
        return NULL;
    }

#if SYSDEBUG
    printlist(ylist, "ylist");
    printlist(xlist, "xlist");
    printlist(plist, "plist");
    fprintf(stderr, "sys->order = %d\n", sys->order);
#endif

    if (!gretl_is_identity_matrix(sys->Gamma)) {
        G = gretl_matrix_copy(sys->Gamma);
        if (G == NULL) {
            *err = E_ALLOC;
        } else {
            *err = gretl_SVD_invert_matrix(G);
        }
        if (*err) {
            gretl_matrix_free(G);
            return NULL;
        }
    }

    T = t2 - t1 + 1;

    y = gretl_matrix_alloc(n, 1);
    yh = gretl_matrix_alloc(n, 1);
    F = gretl_zero_matrix_new(T, 1);

    if (y == NULL || yh == NULL || F == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    if (sys->order > 0) {
        yl = gretl_matrix_alloc(sys->order * ylist[0], 1);
        if (yl == NULL) {
            *err = E_ALLOC;
            goto bailout;
        }
    }

    if (xlist[0] > 0) {
        x = gretl_matrix_alloc(xlist[0], 1);
        if (x == NULL) {
            *err = E_ALLOC;
            goto bailout;
        }
    }

    for (t=t1, s=0; t<=t2; t++, s++) {
        int miss = 0;

        /* lags of endogenous vars */
        if (sys->order > 0) {
            gretl_matrix_zero(yl);
            for (i=1; i<=plist[0] && !miss; i++) {
                vi = plist[i];
                get_col_and_lag(vi, sys, &col, &lag);
                xit = dset->Z[vi][t];
                col += n * (lag - 1);
                if (na(xit)) {
                    miss = 1;
                } else {
                    yl->val[col] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply(sys->A, yl, y);
            }
        } else {
            gretl_matrix_zero(y);
        }

#if SYSDEBUG > 1
        gretl_matrix_print(yl, "yl");
#endif

        /* exogenous vars */
        if (xlist[0] > 0 && !miss) {
            for (i=1; i<=xlist[0] && !miss; i++) {
                xit = dset->Z[xlist[i]][t];
                if (na(xit)) {
                    miss = 1;
                } else {
                    x->val[i-1] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply_mod(sys->B, GRETL_MOD_NONE,
                                          x, GRETL_MOD_NONE,
                                          y, GRETL_MOD_CUMULATE);
            }
        }

        if (miss) {
            gretl_vector_set(sys->F, s, NADBL);
        } else {
            /* multiply by Gamma^{-1} */
            if (G != NULL) {
                gretl_matrix_multiply(G, y, yh);
            } else {
                gretl_matrix_multiply(sys->Gamma, y, yh);
            }
            gretl_vector_set(F, s, yh->val[v]);
        }
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(yh);
    gretl_matrix_free(yl);
    gretl_matrix_free(x);
    gretl_matrix_free(G);

    if (*err) {
        gretl_matrix_free(F);
        F = NULL;
    } else {
        gretl_matrix_set_t1(F, t1);
        gretl_matrix_set_t2(F, t2);
    }

#if SYSDEBUG
    gretl_matrix_print(F, "F: calculated fitted values");
#endif

    return F;
}

/* trim the first sys->order rows off sys->F */

static int sys_F_trim (equation_system *sys)
{
    int err, p = sys->order;
    int r = sys->F->rows - p;
    gretl_matrix *tmp = gretl_matrix_alloc(r, sys->F->cols);

    if (tmp == NULL) {
        err = E_ALLOC;
    } else {
        err = gretl_matrix_extract_matrix(tmp, sys->F, p, 0,
                                          GRETL_MOD_NONE);
    }

    if (!err) {
        gretl_matrix_free(sys->F);
        sys->F = tmp;
    }

    return err;
}

static int sys_add_forecast (equation_system *sys,
                             int t1, int t2,
                             const DATASET *dset,
                             gretlopt opt)
{
    gretl_matrix *G = NULL;
    gretl_matrix *yl = NULL, *x = NULL;
    gretl_matrix *y = NULL, *yh = NULL;
    const int *ylist = sys->ylist;
    const int *xlist = sys->xlist;
    const int *plist = sys->plist;
    int n = sys->neqns + sys->nidents;
    double xit, xitd;
    int tdyn, col, T, ncols;
    int i, vi, s, t, lag, s0 = 0;
    int err = 0;

    if (sys->Gamma == NULL || sys->Gamma->rows != n) {
        fprintf(stderr, "sys_add_forecast: Gamma is broken!\n");
        return E_DATA;
    }

#if SYSDEBUG
    fprintf(stderr, "*** sys_add_forecast\n");
    printlist(ylist, "ylist");
    printlist(xlist, "xlist");
    printlist(plist, "plist");
    fprintf(stderr, "sys->order = %d\n", sys->order);
    fprintf(stderr, "t1=%d, t2=%d, out-of-sample=%d\n",
            t1, t2, (opt & OPT_O)? 1 : 0);
#endif

    if (!gretl_is_identity_matrix(sys->Gamma)) {
        G = gretl_matrix_copy(sys->Gamma);
        if (G == NULL) {
            err = E_ALLOC;
        } else {
            err = gretl_SVD_invert_matrix(G);
        }
        if (err) {
            gretl_matrix_free(G);
            return err;
        }
    }

    /* At which observation (tdyn) does the forecast turn
       dynamic, if at all? */

    if ((opt & OPT_S) || sys->A == NULL) {
        /* got the --static option (or no dynamics): never */
        tdyn = t2 + 1;
    } else if (opt & OPT_D) {
        /* got the --dynamic option: from the start */
        tdyn = t1;
    } else {
        /* by default, for a model with dynamics: just after
           the estimation sample ends */
        tdyn = sys->t2 + 1;
    }

    T = t2 - t1 + 1;
    ncols = 2 * n;

    if (sys->order > 0 && (opt & OPT_O)) {
        /* out-of-sample option: we'll need some extra rows
           for pre-sample values
        */
        T += sys->order;
    }

    y = gretl_matrix_alloc(n, 1);
    yh = gretl_matrix_alloc(n, 1);
    sys->F = gretl_zero_matrix_new(T, ncols);

    if (y == NULL || yh == NULL || sys->F == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    if (sys->order > 0) {
        yl = gretl_matrix_alloc(sys->order * ylist[0], 1);
        if (yl == NULL) {
            err = E_ALLOC;
            goto bailout;
        }
        if (opt & OPT_O) {
            /* out-of-sample */
            int p = sys->order;
            double yti;

            for (t=t1-p, s=0; t<t1; t++, s++) {
                for (i=0; i<n; i++) {
                    yti = dset->Z[ylist[i+1]][t];
                    gretl_matrix_set(sys->F, s, i, yti);
                }
            }
            s0 = p;
        }
    }

    if (xlist[0] > 0) {
        x = gretl_matrix_alloc(xlist[0], 1);
        if (x == NULL) {
            err = E_ALLOC;
            goto bailout;
        }
    }

    for (t=t1, s=s0; t<=t2; t++, s++) {
        int miss = 0;

        /* lags of endogenous vars */
        if (sys->order > 0) {
            gretl_matrix_zero(yl);
            for (i=1; i<=plist[0] && !miss; i++) {
                vi = plist[i];
                get_col_and_lag(vi, sys, &col, &lag);
                xitd = NADBL;
                if (t < tdyn || s - lag < 0) {
                    /* pre-forecast value */
                    xit = dset->Z[vi][t];
                } else {
                    /* prior forecast value preferred */
                    if (s - lag >= 0) {
                        xitd = xit = dset->Z[vi][t];
                    }
                    xit = gretl_matrix_get(sys->F, s - lag, col);
                }
                col += n * (lag - 1);
                if (!na(xit)) {
                    yl->val[col] = xit;
                } else if (!na(xitd)) {
                    yl->val[col] = xitd;
                } else {
                    miss = 1;
                }
            }
            if (!miss) {
                gretl_matrix_multiply(sys->A, yl, y);
            }
        } else {
            gretl_matrix_zero(y);
        }

#if SYSDEBUG > 1
        gretl_matrix_print(yl, "yl");
#endif

        /* exogenous vars */
        if (xlist[0] > 0 && !miss) {
            for (i=1; i<=xlist[0] && !miss; i++) {
                xit = dset->Z[xlist[i]][t];
                if (na(xit)) {
                    miss = 1;
                } else {
                    x->val[i-1] = xit;
                }
            }
            if (!miss) {
                gretl_matrix_multiply_mod(sys->B, GRETL_MOD_NONE,
                                          x, GRETL_MOD_NONE,
                                          y, GRETL_MOD_CUMULATE);
            }
        }

        if (miss) {
            for (i=0; i<n; i++) {
                gretl_matrix_set(sys->F, s, i, NADBL);
            }
        } else {
            /* multiply by Gamma^{-1} */
            if (G != NULL) {
                gretl_matrix_multiply(G, y, yh);
            } else {
                gretl_matrix_multiply(sys->Gamma, y, yh);
            }
            for (i=0; i<n; i++) {
                gretl_matrix_set(sys->F, s, i, yh->val[i]);
            }
        }
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(yh);
    gretl_matrix_free(yl);
    gretl_matrix_free(x);
    gretl_matrix_free(G);

    if (!err && s0 > 0) {
        err = sys_F_trim(sys);
    }

    if (err) {
        gretl_matrix_free(sys->F);
        sys->F = NULL;
    } else {
        gretl_matrix_set_t1(sys->F, t1);
        gretl_matrix_set_t2(sys->F, t2);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->F, "sys->F, forecasts only");
#endif

    if (!err) {
        sys_add_fcast_variance(sys, sys->F, tdyn - t1);
    }

#if SYSDEBUG
    gretl_matrix_print(sys->F, "sys->F, with variance");
#endif

    return err;
}

const gretl_matrix *
system_get_forecast_matrix (equation_system *sys, int t1, int t2,
                            DATASET *dset, gretlopt opt,
                            int *err)
{
    if (sys->F != NULL) {
        gretl_matrix_replace(&sys->F, NULL);
    }

    *err = sys_add_forecast(sys, t1, t2, dset, opt);

    return sys->F;
}

/* attach to the system some extra statistics generated in the course
   of estimation via LIML */

static int sys_attach_ldata (equation_system *sys)
{
    int i, n = sys->neqns;
    int err = 0;

    sys->ldata = malloc(sizeof *sys->ldata);

    if (sys->ldata == NULL) {
        err = E_ALLOC;
    } else {
        sys->ldata->lmin = NULL;
        sys->ldata->ll = NULL;
        sys->ldata->idf = NULL;

        sys->ldata->lmin = malloc(n * sizeof *sys->ldata->lmin);
        sys->ldata->ll = malloc(n * sizeof *sys->ldata->ll);
        sys->ldata->idf = malloc(n * sizeof *sys->ldata->idf);

        if (sys->ldata->lmin == NULL ||
            sys->ldata->ll == NULL ||
            sys->ldata->idf == NULL) {
            free(sys->ldata->lmin);
            free(sys->ldata->ll);
            free(sys->ldata->idf);
            free(sys->ldata);
            sys->ldata = NULL;
            err = E_ALLOC;
        }
    }

    if (!err) {
        const MODEL *pmod;

        for (i=0; i<n; i++) {
            pmod = sys->models[i];
            sys->ldata->lmin[i] = gretl_model_get_double(pmod, "lmin");
            sys->ldata->ll[i] = pmod->lnL;
            sys->ldata->idf[i] = gretl_model_get_int(pmod, "idf");
        }
    }

    return err;
}

static int sys_print_reconstituted_models (const equation_system *sys,
                                           const DATASET *dset,
                                           PRN *prn)
{
    MODEL mod;
    const double *y;
    double x;
    int print_insts = 0;
    int ifc = 0, nc = 0, ncmax = 0;
    int i, j, t, k = 0;
    int err = 0;

    gretl_model_init(&mod, dset);

    mod.t1 = sys->t1;
    mod.t2 = sys->t2;
    mod.nobs = sys->T;
    mod.aux = AUX_SYS;
    gretl_model_set_int(&mod, "method", sys->method);

    mod.lnL = NADBL;

    if (sys->method == SYS_METHOD_TSLS ||
        sys->method == SYS_METHOD_3SLS ||
        sys->method == SYS_METHOD_FIML ||
        sys->method == SYS_METHOD_LIML) {
        mod.ci = IVREG;
    } else {
        mod.ci = OLS;
    }

    print_insts = sys->method == SYS_METHOD_TSLS ||
        sys->method == SYS_METHOD_3SLS;

    for (i=0; i<sys->neqns && !err; i++) {
        int *mlist = sys->lists[i];
        int freelist = 0;

        mod.ID = i;

        if (print_insts && !gretl_list_has_separator(mlist)) {
            mod.list = compose_ivreg_list(sys, i);
            if (mod.list == NULL) {
                err = E_ALLOC;
            } else {
                freelist = 1;
            }
        } else {
            mod.list = mlist;
        }

        if (!err) {
            ifc = nc = 0;
            for (j=2; j<=mod.list[0]; j++) {
                if (mod.list[j] == 0) {
                    ifc = 1;
                } else if (mod.list[j] == LISTSEP) {
                    break;
                }
                nc++;
            }
        }

        if (!err && nc > ncmax) {
            mod.coeff = realloc(mod.coeff, nc * sizeof *mod.coeff);
            mod.sderr = realloc(mod.sderr, nc * sizeof *mod.sderr);
            ncmax = nc;
        }

        if (mod.coeff == NULL || mod.sderr == NULL) {
            err = E_ALLOC;
            break;
        }

        mod.ncoeff = nc;
        mod.dfn = nc - ifc;

        if (sys->flags & SYSTEM_DFCORR) {
            gretl_model_set_int(&mod, "dfcorr", 1);
            mod.dfd = mod.nobs - nc;
        } else {
            mod.dfd = mod.nobs;
        }

        y = dset->Z[mod.list[1]];
        mod.ybar = gretl_mean(mod.t1, mod.t2, y);
        mod.sdy = gretl_stddev(mod.t1, mod.t2, y);

        mod.ess = 0.0;
        for (t=0; t<sys->T; t++) {
            x = gretl_matrix_get(sys->E, t, i);
            mod.ess += x * x;
        }

        if (sys->method == SYS_METHOD_OLS ||
            sys->method == SYS_METHOD_TSLS ||
            sys->method == SYS_METHOD_LIML) {
            /* single-equation methods */
            mod.sigma = sqrt(mod.ess / mod.dfd);
        } else {
            mod.sigma = sqrt(gretl_matrix_get(sys->S, i, i));
        }

        for (j=0; j<mod.ncoeff; j++) {
            mod.coeff[j] = sys->b->val[k];
            mod.sderr[j] = sqrt(gretl_matrix_get(sys->vcv, k, k));
            k++;
        }

        if (sys->method == SYS_METHOD_LIML && sys->ldata != NULL) {
            gretl_model_set_double(&mod, "lmin", sys->ldata->lmin[i]);
            mod.lnL = sys->ldata->ll[i];
            gretl_model_set_int(&mod, "idf", sys->ldata->idf[i]);
        }

        printmodel(&mod, dset, OPT_NONE, prn);

        if (freelist) {
            free(mod.list);
        }

        mod.list = NULL;
    }

    clear_model(&mod);

    return err;
}

int gretl_system_print (equation_system *sys, const DATASET *dset,
                        gretlopt opt, PRN *prn)
{
    const char *name = sys->name;
    gchar *sysstr = NULL;
    int tex = tex_format(prn);
    int nr = system_n_restrictions(sys);
    int i;

    if (sys->models != NULL &&
        sys->method == SYS_METHOD_LIML &&
        sys->ldata == NULL) {
        sys_attach_ldata(sys);
    }

    if (name != NULL && sys_anonymous(name)) {
        /* don't print internal reference name */
        name = NULL;
    }

    if (tex) {
        pputs(prn, "\\begin{center}\n");
    } else {
	pputc(prn, '\n');
    }

    sysstr = system_get_full_string(sys);
    if (name != NULL) {
	if (tex) {
	    pprintf(prn, "%s, %s\\\\\n", _("Equation system"), name);
	} else {
	    pprintf(prn, "%s, %s\n", _("Equation system"), name);
	}
	pprintf(prn, "%s: %s", _("Estimator"), sysstr);
    } else {
	pprintf(prn, "%s, %s", _("Equation system"), sysstr);
    }
    g_free(sysstr);

    if (sys->iters > 0) {
        gretl_prn_newline(prn);
	pprintf(prn, _("Convergence achieved after %d iterations\n"), sys->iters);
        if (sys->method == SYS_METHOD_SUR ||
            sys->method == SYS_METHOD_FIML) {
            if (tex) {
                gretl_prn_newline(prn);
                pprintf(prn, "%s = ", _("Log-likelihood"));
                if (sys->ll < 0) {
                    pprintf(prn, "$-$%g", -sys->ll);
                } else {
                    pprintf(prn, "%g", sys->ll);
                }
            } else {
                pprintf(prn, "%s = %g\n", _("Log-likelihood"), sys->ll);
            }
        }
    }

    if (tex) {
        pputs(prn, "\n\\end{center}\n\n");
    } else {
        pputc(prn, '\n');
    }

    if (sys->models != NULL) {
        for (i=0; i<sys->neqns; i++) {
            if (sys->flags & SYSTEM_DFCORR) {
                gretl_model_set_int(sys->models[i], "dfcorr", 1);
            }
            printmodel(sys->models[i], dset, OPT_NONE, prn);
        }
    } else {
        sys_print_reconstituted_models(sys, dset, prn);
    }

    system_print_sigma(sys, prn);

    if (nr == 0 && (sys->method == SYS_METHOD_FIML ||
                    sys->method == SYS_METHOD_3SLS ||
                    sys->method == SYS_METHOD_SUR)) {
        print_system_overid_test(sys, prn);
    }

    return 0;
}

static void ensure_asy_printout (equation_system *sys)
{
    int i;

    if (sys->models != NULL) {
        for (i=0; i<sys->neqns; i++) {
            sys->models[i]->ci = IVREG;
        }
    }
}

int
system_save_and_print_results (equation_system *sys, DATASET *dset,
                               gretlopt opt, PRN *prn)
{
    int nr = system_n_restrictions(sys);
    int err = 0;

    if (sys->E != NULL) {
        gretl_matrix_set_t1(sys->E, sys->t1);
        gretl_matrix_set_t2(sys->E, sys->t2);
    }

    if (sys->iters > 0 && sys->method == SYS_METHOD_SUR && nr == 0) {
        sur_ols_diag(sys);
    }

    sys->ldet = gretl_vcv_log_determinant(sys->S, &err);

    if (!err) {
        err = system_add_yhat_matrix(sys);
    }

    if (!err) {
        err = sys_add_structural_form(sys, dset);
    }

    if (!(opt & OPT_Q)) {
        if (sys->method == SYS_METHOD_FIML) {
            ensure_asy_printout(sys);
        }
        gretl_system_print(sys, dset, opt, prn);
    }

    return err;
}

int system_autocorrelation_test (equation_system *sys, int order,
                                 gretlopt opt, PRN *prn)
{
    double *u, lb;
    int i, err = 0;

    for (i=0; i<sys->neqns && !err; i++) {
        if (!(opt & OPT_Q)) {
            pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
        }
        u = sys->E->val + (i * sys->T);
        lb = ljung_box(order, 0, sys->T - 1, u, &err);
        if (!err && !((opt & OPT_Q))) {
            pprintf(prn, "%s: %s(%d) = %g [%.4f]\n\n",
                    _("Ljung-Box Q'"), _("Chi-square"), order,
                    lb, chisq_cdf_comp(order, lb));
        }
    }

    return err;
}

int system_arch_test (equation_system *sys, int order,
                      gretlopt opt, PRN *prn)
{
    const double *u;
    int i, err = 0;

    if (!(opt & OPT_I)) {
        pputc(prn, '\n');
    }

    for (i=0; i<sys->neqns && !err; i++) {
        if (!(opt & OPT_I)) {
            pprintf(prn, "%s %d:\n", _("Equation"), i + 1);
        }
        u = sys->E->val + (i * sys->T);
        err = array_arch_test(u, sys->T, order, opt | OPT_M, prn);
    }

    return err;
}

int system_supports_method (equation_system *sys, int method)
{
    if (sys != NULL) {
        int i;

        for (i=0; i<sys->neqns; i++) {
            if (gretl_list_has_separator(sys->lists[i])) {
                return method == SYS_METHOD_TSLS ||
                    method == SYS_METHOD_3SLS;
            }
        }
    }

    return 1;
}

static void finalize_liml_model (MODEL *pmod, equation_system *sys)
{
    *pmod = *sys->models[0];

#if SYSDEBUG
    display_model_data_items(pmod);
#endif

    gretl_model_destroy_data_item(pmod, "tslsX");
    gretl_model_destroy_data_item(pmod, "endog");
    gretl_model_destroy_data_item(pmod, "method");
    gretl_model_destroy_data_item(pmod, "liml_y");

#if SYSDEBUG
    display_model_data_items(pmod);
#endif

    free(sys->models[0]);
    free(sys->models);
    sys->models = NULL;

    pmod->aux = AUX_NONE;
    pmod->opt |= OPT_L;
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = NADBL;
    set_model_id(pmod, OPT_NONE);
}

/* implement the --liml option to "tsls" */

MODEL single_equation_liml (const int *list, DATASET *dset,
                            gretlopt opt)
{
    int *mlist = NULL, *ilist = NULL;
    equation_system *sys = NULL;
    MODEL model;
    int err = 0;

    gretl_model_init(&model, dset);

    err = ivreg_process_lists(list, &mlist, &ilist);

    if (!err) {
        sys = equation_system_new(SYS_METHOD_LIML, NULL, &err);
    }
    if (!err) {
        err = equation_system_append(sys, mlist);
    }
    if (!err) {
        sys->ilist = ilist;
        err = equation_system_finalize(sys, dset, OPT_S, NULL);
    }

    if (err) {
        model.errcode = err;
    } else {
        finalize_liml_model(&model, sys);
    }

    equation_system_destroy(sys);
    free(mlist);

    return model;
}

int gretl_system_get_sample (const equation_system *sys,
                             int *t1, int *t2)
{
    if (sys != NULL) {
        *t1 = sys->smpl_t1;
        *t2 = sys->smpl_t2;
        return 0;
    }

    return E_DATA;
}
