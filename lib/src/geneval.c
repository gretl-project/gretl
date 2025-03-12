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

/* syntax tree evaluator for "genr" */

#include "genparse.h"
#include "gen_public.h"
#include "monte_carlo.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "uservar.h"
#include "gretl_bfgs.h"
#include "gretl_normal.h"
#include "gretl_panel.h"
#include "kalman.h"
#include "libset.h"
#include "version.h"
#include "csvdata.h"
#include "uservar_priv.h"
#include "genr_optim.h"
#include "gretl_cmatrix.h"
#include "qr_estimate.h"
#include "gretl_foreign.h"
#include "gretl_midas.h"
#include "gretl_tdisagg.h"
#include "gretl_xml.h"
#include "gretl_mt.h"
#include "gretl_color.h"
#include "var.h"
#include "vartest.h"
#include "flow_control.h"
#include "mapinfo.h"

#include <time.h> /* for the $now accessor */

#ifdef USE_CURL
# include "gretl_www.h"
#endif

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#ifdef WIN32
# include "gretl_win32.h" /* for strptime() */
#endif

#include <errno.h>

#if GENDEBUG
# define EDEBUG GENDEBUG
# define LHDEBUG GENDEBUG
#else
# define EDEBUG 0
# define LHDEBUG 0
#endif

#if LHDEBUG || EDEBUG > 1
# define IN_GENEVAL
# include "mspec_debug.c"
#endif

#define AUX_NODES_DEBUG 0
#define AGGR_PREFER_MATRIX 1

#if AUX_NODES_DEBUG
# include <stdarg.h>
static void real_rndebug (const char *format, ...)
{
    va_list args;

    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}
# define rndebug(x) real_rndebug x
#else
# define rndebug(x)
#endif

#define ONE_BY_ONE_CAST 1

enum {
    FR_TREE = 1,
    FR_RET,
    FR_LHTREE,
    FR_LHRES,
    FR_ARET
};

#define is_aux_node(n) (n != NULL && (n->flags & AUX_NODE))
#define is_tmp_node(n) (n != NULL && (n->flags & TMP_NODE))
#define is_proxy_node(n) (n != NULL && (n->flags & PRX_NODE))

#define emptymat_ok(f) (f == F_GINV || f == F_DIAG || f == F_TRANSP || \
                        f == F_VEC || f == F_VECH || f == F_UNVECH ||   \
                        f == F_CHOL || f == F_UPPER || f == F_LOWER ||  \
                        f == F_SORT || f == F_DSORT || f == F_VALUES || \
                        f == F_MREV)

#define dataset_dum(n) (n->t == DUM && n->v.idnum == DUM_DATASET)

#define postfix_node(n) (n->t == NUM_P || n->t == NUM_M)

#define uscalar_node(n) ((n->t == NUM && n->vname != NULL) || postfix_node(n))

#define umatrix_node(n) (n->t == MAT && n->vname != NULL)
#define ubundle_node(n) (n->t == BUNDLE && n->vname != NULL)
#define uarray_node(n)  (n->t == ARRAY && n->vname != NULL)
#define ustring_node(n) (n->t == STR && n->vname != NULL)
#define ulist_node(n)   (n->t == LIST && n->vname != NULL)
#define useries_node(n) (n->t == SERIES && n->vnum >= 0)
#define uvar_node(n)    (n->vname != NULL)

#define scalar_matrix_node(n) (n->t == MAT && gretl_matrix_is_scalar(n->v.m))
#define scalar_node(n) (n->t == NUM || scalar_matrix_node(n))
#define ok_matrix_node(n) (n->t == MAT || n->t == NUM)
#define complex_node(n) (n->t == MAT && n->v.m->is_complex)
#define cscalar_node(n) (n->t == MAT && gretl_matrix_is_cscalar(n->v.m))

#define strvals_node(n,p) (n->vnum > 0 && is_string_valued(p->dset, n->vnum))
#define mutable_node(n) (n->flags & MUT_NODE)

#define null_node(n) (n == NULL || n->t == EMPTY)
#define null_or_scalar(n) (null_node(n) || scalar_node(n))
#define null_or_string(n) (null_node(n) || n->t == STR)

#define ok_bundled_type(t) (t == NUM || t == STR || t == MAT || t == LIST || \
                            t == SERIES || t == BUNDLE || t == ARRAY)

#define compiled(p) (p->flags & P_EXEC)
#define starting(p) (p->flags & P_START)
#define autoreg(p)  (p->flags & P_AUTOREG)
#define exestart(p) ((p->flags & P_EXEC) && (p->flags & P_START))

static void parser_init (parser *p, const char *str, DATASET *dset,
                         PRN *prn, int flags, int targtype, int *done);
static void parser_reinit (parser *p, DATASET *dset, PRN *prn);
static NODE *eval (NODE *t, parser *p);
static void node_type_error (int ntype, int argnum, int goodt,
                             NODE *bad, parser *p);
static int node_is_true (NODE *n, parser *p);
static gretl_matrix *series_to_matrix (const double *x, parser *p);
static void printnode (NODE *t, parser *p, int value);
static inline int attach_aux_node (NODE *t, NODE *ret, parser *p);
static char *get_opstr (int op);
static gretl_bundle *node_get_bundle (NODE *n, parser *p);
static int gen_type_is_arrayable (int gen_t);
static GretlType gretl_type_from_gen_type (int gen_t);

static int user_qsorting;

void set_user_qsorting (int s)
{
    if (s != 0) {
        user_qsorting = 1;
        gretl_iteration_push();
    } else {
        user_qsorting = 0;
        gretl_iteration_pop();
    }
}

/* ok_list_node: This is a first-pass assessment of whether a given
   node _may_ be interpretable as holding a LIST.  The follow-up is
   node_get_list(), which will determine whether the interpretation
   really works. The code below dates from 2023-05-11, when we
   reverted to an earlier version [a9ec40]. In the interim, version
   [b6ed06] had "if (null_node(n))" as the first test and "return
   scalar_node(n)" as the consequent of the last test, giving a
   somewhat different take on the case where @n is an EMPTY node (as
   opposed to NULL). It seems that the interim version broke some
   things.
*/

static int ok_list_node (NODE *n, parser *p)
{
    if (n == NULL) {
        return 0;
    } else if (n->t == LIST) {
        return 1;
    } else if (n->t == SERIES && n->vnum >= 0) {
        /* can interpret as singleton list */
        return 1;
    } else if (p->flags & P_LISTDEF) {
        /* when defining a list we can be a bit more accommodating */
        return null_or_scalar(n);
    }

    return 0;
}

/* a more "lenient" version of the above, to accommodate list
   expressions such as (L - 0), indicating the list that results from
   dropping the constant from L
*/

static int ok_list_node_plus (NODE *n)
{
    if (n->t == LIST) {
        return 1;
    } else if (n->t == SERIES && n->vnum >= 0) {
        return 1;
    } else if (n->t == NUM) {
        return 1;
    } else {
        return 0;
    }
}

static const char *typestr (int t)
{
    switch (t) {
    case NUM:
        return "scalar";
    case SERIES:
        return "series";
    case MAT:
        return "matrix";
    case STR:
        return "string";
    case U_ADDR:
        return "address";
    case LIST:
        return "list";
    case BUNDLE:
        return "bundle";
    case DBUNDLE:
        return "$-bundle";
    case ARRAY:
        return "array";
    case USERIES:
        return "named series";
    case EMPTY:
        return "empty";
    case UNDEF:
        return "undefined";
    default:
        return "?";
    }
}

static void free_mspec (matrix_subspec *spec, parser *p)
{
    if (spec != NULL) {
        free(spec->rslice);
        free(spec->cslice);
        free(spec);
    }
}

static void clear_mspec (matrix_subspec *spec, parser *p)
{
    free(spec->rslice);
    free(spec->cslice);

    memset(spec, 0, sizeof(*spec));
}

/* return 1 if @n equals @t or is found under @t, otherwise 0 */

int find_in_tree (NODE *t, NODE *n)
{
    if (t == NULL) {
        return 0;
    } else if (t == n) {
        return 1;
    } else if (find_in_tree(t->L, n)) {
        return 1;
    } else if (find_in_tree(t->R, n)) {
        return 1;
    } else if (find_in_tree(t->M, n)) {
        return 1;
    } else if (bnsym(t->t)) {
        int i;

        for (i=0; i<t->v.bn.n_nodes; i++) {
            if (find_in_tree(t->v.bn.n[i], n)) {
                return 1;
            }
        }
    }

    return 0;
}

#if EDEBUG || LHDEBUG

static char *flagstr (guint8 flags)
{
    static char ret[16];

    if (flags & AUX_NODE) {
        strcpy(ret, "aux");
        if (flags & TMP_NODE) {
            strcat(ret, ",tmp");
        }
    } else if (flags & TMP_NODE) {
        strcpy(ret, "tmp");
    } else if (flags & LHT_NODE) {
        strcpy(ret, "lht");
    } else {
        sprintf(ret, "%d", (int) flags);
    }

    return ret;
}

static void print_node_details (NODE *t)
{
    fprintf(stderr, "node at %p (type %03d, %s, flags %s",
            (void *) t, t->t, getsymb(t->t), flagstr(t->flags));
    if (t->t == STR) {
        fprintf(stderr, ", val '%s')\n", t->v.str);
    } else if (t->t == NUM) {
        fprintf(stderr, ", val %g)\n", t->v.xval);
    } else {
        fputs(")\n", stderr);
    }
}

static void print_tree (NODE *t, parser *p, int level, char pos)
{
    if (t == NULL) {
        fprintf(stderr, " %d: node is null\n", level);
        return;
    }

    if (bnsym(t->t)) {
        int i;

        for (i=0; i<t->v.bn.n_nodes; i++) {
            print_tree(t->v.bn.n[i], p, level+1, 0);
        }
    } else {
        if (t->L != NULL) {
            print_tree(t->L, p, level+1, 'L');
        }
        if (t->M != NULL) {
            print_tree(t->M, p, level+1, 'M');
        }
        if (t->R != NULL) {
            print_tree(t->R, p, level+1, 'R');
        }
    }

    if (pos != 0) {
        fprintf(stderr, " %d (%c): ", level, pos);
    } else {
        fprintf(stderr, " %d: ", level);
    }

    if (t->vname != NULL) {
        fprintf(stderr, "node at %p (type %03d, %s, flags %s), vname='%s'",
                (void *) t, t->t, getsymb(t->t), flagstr(t->flags), t->vname);
        if (t->t == NUM) {
            fprintf(stderr, ", val %g\n", t->v.xval);
        } else {
            fputc('\n', stderr);
        }
    } else {
        print_node_details(t);
    }

    if (t->aux != NULL) {
        fprintf(stderr, "  aux node at %p (type %03d, %s, flags %s)\n",
                (void *) t->aux, t->aux->t, getsymb(t->aux->t),
                flagstr(t->aux->flags));
    }
}

#endif /* EDEBUG || LHDEBUG */

#if EDEBUG

static const char *free_tree_tag (int t)
{
    if (t == FR_TREE) {
        return "free tree";
    } else if (t == FR_RET) {
        return "free ret";
    } else if (t == FR_LHTREE) {
        return "free lhtree";
    } else if (t == FR_LHRES) {
        return "free lhres";
    } else {
        return "free other";
    }
}

#endif /* EDEBUG */

/* used when we know that @t is a terminal node: skip
   the tests for attached tree */

static void free_node (NODE *t, parser *p)
{
    if (t->refcount > 1) {
        rndebug(("free node %p (%s): decrement refcount to %d\n",
                 (void *) t, getsymb(t->t), t->refcount - 1));
        t->refcount -= 1;
        return;
    }

    if (is_tmp_node(t)) {
#if EDEBUG
        fprintf(stderr, " tmp node: freeing attached data\n");
#endif
        if (t->t == SERIES) {
            free(t->v.xvec);
        } else if (t->t == LIST || t->t == IVEC) {
            free(t->v.ivec);
        } else if (t->t == MAT) {
            gretl_matrix_free(t->v.m);
        } else if (t->t == MSPEC) {
            free_mspec(t->v.mspec, p);
        } else if (t->t == BUNDLE) {
            gretl_bundle_destroy(t->v.b);
        } else if (t->t == ARRAY) {
            gretl_array_destroy(t->v.a);
        } else if (t->t == STR) {
            free(t->v.str);
        }
    }

    if (t->t == UOBJ || t->t == WLIST) {
        free(t->v.str);
    }

    if (t->vname != NULL) {
        free(t->vname);
    }

    if (p != NULL && t == p->ret) {
        p->ret = NULL;
    }

    free(t);
}

/* A word on "aux" nodes. These come in two sorts, which
   might be described as "robust" and "fragile" respectively.

   A robust node (identified by the TMP_NODE flag) is one
   whose data pointer is independently allocated. With such
   a node it's OK simply to "pass on" the pointer in
   assignment, and if it's not passed on it should be freed
   on completion of "genr". (So nota bene: if it's assigned
   elsewhere, the pointer on the aux node itself must then
   be set to NULL to avoid double-freeing.)

   A fragile node is one whose data pointer is not
   independently allocated; it actually "belongs to someone
   else". In assignment, then, it must be deeply copied,
   and it must _not_ be freed on completion of genr.

   Obviously, it's necessary to be careful in handling
   fragile nodes, but the advantage of allowing them
   is that they cut down on wasteful deep-copying of objects
   that may be used in calculation, without being modified,
   on the fly.
*/

void free_tree (NODE *t, parser *p, int code)
{
    if (t == NULL) {
        return;
    }

#if EDEBUG
    fprintf(stderr, "%-11s: starting at %p (type %03d, %s)\n",
            free_tree_tag(code), (void *) t, t->t,
            getsymb(t->t));
#endif

    /* free recursively */
    if (bnsym(t->t)) {
        int i;

        for (i=0; i<t->v.bn.n_nodes; i++) {
            free_tree(t->v.bn.n[i], p, code);
        }
        free(t->v.bn.n);
    } if (!(t->flags & LHT_NODE)) {
        free_tree(t->L, p, code);
        free_tree(t->M, p, code);
        free_tree(t->R, p, code);
    }

    if (t->aux != NULL && t->aux != p->ret && t->aux != p->lhres) {
        rndebug(("freeing aux node at %p (%s)\n", (void *) t->aux,
                 getsymb(t->aux->t)));
        free_node(t->aux, p);
    } else if (t->aux != NULL) {
        rndebug(("NOT freeing aux at %p (= p->ret)\n", (void *) t->aux));
        t->aux->refcount -= 1;
    }

#if EDEBUG
    fprintf(stderr, "%-11s: freeing node at %p (type %03d, %s, flags = %d)\n",
            free_tree_tag(code), (void *) t, t->t, getsymb(t->t),
            t->flags);
#endif

    free_node(t, p);
}

static void clear_uvnodes (NODE *t)
{
    if (t == NULL) {
        return;
    }

    if (bnsym(t->t)) {
        int i;

        for (i=0; i<t->v.bn.n_nodes; i++) {
            clear_uvnodes(t->v.bn.n[i]);
        }
    } else {
        clear_uvnodes(t->L);
        clear_uvnodes(t->M);
        clear_uvnodes(t->R);
    }

    if (t->t == SERIES) {
        if (t->vnum >= 0 || t->vname != NULL) {
#if EDEBUG
            fprintf(stderr, " clear_uvnode: series at %p\n", (void *) t);
#endif
            t->v.xvec = NULL;
        }
    } else if (t->uv != NULL) {
#if EDEBUG
        fprintf(stderr, " clear_uvnode: uvar '%s' at %p, uv %p\n", t->vname,
                (void *) t, (void *) t->uv);
#endif
        t->uv = NULL;
    }
}

#if AUX_NODES_DEBUG
static void reset_p_aux (parser *p, NODE *n)
{
    fprintf(stderr, "resetting p->aux = %p\n", (void *) n);
    p->aux = n;
}
#else
# define reset_p_aux(p, n) (p->aux = n)
#endif

static double *na_array (int n)
{
    double *x = malloc(n * sizeof *x);
    int i;

    if (x != NULL) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    return x;
}

/* new node to hold array of doubles */

static NODE *newseries (int n, int flags)
{
    NODE *b = new_node(SERIES);

    if (b != NULL) {
        b->flags = flags;
        if (n > 0) {
            b->v.xvec = na_array(n);
            if (b->v.xvec == NULL) {
                free(b);
                b = NULL;
            }
        } else {
            b->v.xvec = NULL;
        }
    }

    return b;
}

/* new node to hold array of @n ints */

static NODE *newivec (int n)
{
    NODE *b = new_node(IVEC);

    if (b != NULL) {
        b->flags = TMP_NODE;
        if (n > 0) {
            b->v.ivec = malloc(n * sizeof(int));
            if (b->v.ivec == NULL) {
                free(b);
                b = NULL;
            }
        } else {
            b->v.ivec = NULL;
        }
    }

    return b;
}

/* new node to hold a gretl_matrix */

static NODE *newmat (int flags)
{
    NODE *n = new_node(MAT);

    if (n != NULL) {
        n->flags = flags;
        n->v.m = NULL;
    }

    return n;
}

/* new node to hold a matrix specification */

static NODE *newmspec (void)
{
    NODE *n = new_node(MSPEC);

    if (n != NULL) {
        n->flags = TMP_NODE;
        n->v.mspec = NULL;
    }

    return n;
}

/* new node to hold a list */

static NODE *newlist (int flags)
{
    NODE *n = new_node(LIST);

    if (n != NULL) {
        n->flags = flags;
        n->v.ivec = NULL;
    }

    return n;
}

static NODE *newstring (int flags)
{
    NODE *n = new_node(STR);

    if (n != NULL) {
        n->flags = flags;
        n->v.str = NULL;
    }

    return n;
}

static NODE *newbundle (int flags)
{
    NODE *n = new_node(BUNDLE);

    if (n != NULL) {
        n->flags = flags;
        n->v.b = NULL;
    }

    return n;
}

static NODE *newarray (int flags)
{
    NODE *n = new_node(ARRAY);

    if (n != NULL) {
        n->flags = flags;
        n->v.a = NULL;
    }

    return n;
}

static void clear_tmp_node_data (NODE *n, parser *p)
{
    int nullify = 1;

#if 0
    /* 2022-09-19: The look-up below (added quite recently) slows
       things down significantly. But it now seems to be redundant;
       that is, the current code doesn't make the mistake of putting a
       named "uservar" into a tmp node (which is supposed to hold only
       disposable intermediate results). I'm leaving it here, def'd
       out, just in case that error somehow reappears.
    */
    user_var *uv = NULL;

    if ((uv = get_user_var_by_data(n->v.ptr)) != NULL) {
        /* should NOT happen */
        fprintf(stderr, "clear_tmp_node_data: n->v.ptr belongs to uvar %s (%p)\n",
                uv->name, (void *) uv);
        return;
    }
#endif

    if (n->t == LIST) {
        free(n->v.ivec);
    } else if (n->t == MAT) {
        /* (how) can we avoid doing this? */
        gretl_matrix_free(n->v.m);
    } else if (n->t == MSPEC) {
        if (n->v.mspec != NULL) {
            clear_mspec(n->v.mspec, p);
        }
        nullify = 0;
    } else if (n->t == BUNDLE) {
        gretl_bundle_destroy(n->v.b);
    } else if (n->t == ARRAY) {
        gretl_array_destroy(n->v.a);
    } else if (n->t == STR) {
        free(n->v.str);
    } else if (n->t == SERIES) {
        /* preserve any existing tmp series, unless the
           dataset series length has changed
        */
        if (p->flags & P_DELTAN) {
            free(n->v.xvec);
            n->v.xvec = NULL;
            if (p->dset_n > 0) {
                n->v.xvec = na_array(p->dset_n);
                if (n->v.xvec == NULL) {
                    p->err = E_ALLOC;
                }
            }
        } else {
            /* scrub any pre-existing values in the current
               sample range */
            int t;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                n->v.xvec[t] = NADBL;
            }
        }
        nullify = 0;
    } else {
        nullify = 0;
    }

    if (nullify) {
        n->v.ptr = NULL;
    }
}

static int mutate_bundle_member_node (NODE *n, int type,
                                      int flags, parser *p)
{
    int err = 0;

#if EDEBUG > 1
    fprintf(stderr, "mutate_bundle_member: %p, %s -> %s, tmp %d -> %d\n",
            (void *) n, getsymb(n->t), getsymb(type),
            (n->flags & TMP_NODE)? 1 : 0,
            (flags & TMP_NODE)? 1 : 0);
#endif

    if (is_tmp_node(n)) {
        /* some allocated storage should be freed */
        if (n->t == SERIES) {
            free(n->v.xvec);
        } else if (n->t == MAT) {
            gretl_matrix_free(n->v.m);
        }
        n->v.ptr = NULL;
    }

    if (type == SERIES) {
        /* switching to a series node : allocate xvec */
        n->v.xvec = malloc(p->dset->n * sizeof(double));
        if (n->v.xvec == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        n->t = type;
        n->flags = flags;
    }

    return err;
}

/* We allow here for some equivocation in type between
   1 x 1 matrices and scalars in the course of executing
   a compiled parser with saved aux nodes.
*/

static void maybe_switch_node_type (NODE *n, int type,
                                    int flags, parser *p)
{
    if (mutable_node(n)) {
        /* bundle members only */
        p->err = mutate_bundle_member_node(n, type, flags, p);
    } else if (n->t == MAT && type == NUM) {
        /* switch aux node @n from matrix to scalar */
        if (is_tmp_node(n)) {
            gretl_matrix_free(n->v.m);
        }
        n->t = NUM;
        n->v.xval = NADBL;
        n->flags = 0;
        n->vnum = NO_VNUM;
        n->vname = NULL;
    } else if (n->t == NUM && type == MAT) {
        /* switch @n from scalar to matrix */
        n->t = MAT;
        n->v.m = NULL;
        n->flags = flags;
    } else if (type == EMPTY) {
        ; /* LHS mechanism: OK */
    } else if (gen_type_is_arrayable(n->t) && type == ARRAY) {
        n->t = ARRAY;
        n->flags = flags;
    } else {
        /* any other discrepancy presumably means that
           things have gone badly wrong
        */
        fprintf(stderr, "aux node mismatch: expected %s but got %s, TMP_NODE=%d\n",
                getsymb(n->t), getsymb(type), (flags & TMP_NODE)? 1 : 0);
        gretl_errmsg_set("internal genr error: aux node mismatch");
        p->err = E_DATA;
    }
}

/* get an auxiliary node: if starting from scratch we allocate
   a new node, otherwise we look up an existing one */

static NODE *get_aux_node (parser *p, int t, int n, int flags)
{
    NODE *ret = p->aux;

#if EDEBUG
    fprintf(stderr, "get_aux_node: t=%s, tmp=%d, starting=%d, "
            "p->aux=%p\n", getsymb(t), (flags & TMP_NODE)? 1 : 0,
            starting(p) ? 1 : 0, (void *) p->aux);
#endif

    if (is_proxy_node(ret)) {
        /* this node will get freed later */
        ret = NULL;
    }

    if (ret != NULL) {
        /* got a pre-existing aux node */
        if (starting(p)) {
            if (ret->t != t) {
                maybe_switch_node_type(ret, t, flags, p);
            } else if (is_tmp_node(ret) && !(p->flags & P_MSAVE)) {
                clear_tmp_node_data(ret, p);
            }
        }
    } else {
        /* we need to create a new aux node */
        if (t == NUM) {
            ret = newdbl(NADBL);
        } else if (t == SERIES) {
            ret = newseries(n, flags);
        } else if (t == IVEC) {
            ret = newivec(n);
        } else if (t == LIST) {
            ret = newlist(flags);
        } else if (t == MAT) {
            ret = newmat(flags);
        } else if (t == MSPEC) {
            ret = newmspec();
        } else if (t == STR) {
            ret = newstring(flags);
        } else if (t == BUNDLE) {
            ret = newbundle(flags);
        } else if (t == ARRAY) {
            ret = newarray(flags);
        } else if (t == EMPTY) {
            ret = newempty();
        } else {
            /* invalid aux node spec */
            p->err = E_DATA;
        }

        if (!p->err && ret == NULL) {
            p->err = E_ALLOC;
        }

        if (!p->err) {
            ret->flags |= AUX_NODE;
        }
    }

    return ret;
}

static int e_types (NODE *n)
{
    if (n != NULL) {
        GretlType t;
        const char *s;

        if (n->t == ARRAY) {
            t = gretl_array_get_type(n->v.a);
            s = gretl_type_get_name(t);
        } else if (n->t == U_ADDR) {
            s = "pointer-to-object";
        } else if (n->t == EMPTY) {
            s = "null";
        } else {
            t = gretl_type_from_gen_type(n->t);
            s = gretl_type_get_name(t);
            if (!strcmp(s, "invalid")) {
                s = getsymb(n->t);
            }
        }
        gretl_errmsg_sprintf("%s is %s, not acceptable in context",
                             n->vname == NULL ? "[anonymous]" : n->vname, s);
    }

    return E_TYPES;
}

/* We come here by preference to the generic get_aux_node()
   (above) if we want an aux node holding an allocated matrix
   of known size (m x n). On the second or subsequent
   iterations of a loop, with any luck we may find that the
   aux node already holds a matrix of the required dimensions
   which can then be reused.
*/

static NODE *aux_sized_matrix_node (parser *p, int m, int n,
                                    int cmplx)
{
    NODE *ret = p->aux;

    if (is_proxy_node(ret)) {
        /* this node will get freed later */
        ret = NULL;
    }

    if (ret != NULL) {
        /* got a pre-existing node */
        if (ret->t == NUM) {
            /* switch @ret from scalar to matrix */
            ret->t = MAT;
            ret->v.m = NULL;
            ret->flags |= TMP_NODE;
        } else if (ret->t != MAT) {
            p->err = e_types(ret);
        } else {
            /* check for reusable matrix */
            gretl_matrix *a = ret->v.m;

            if (a != NULL) {
                if (a->is_complex + cmplx == 1) {
                    /* too difficult to reuse */
                    gretl_matrix_free(ret->v.m);
                    ret->v.m = NULL;
                } else if (a->rows != m || a->cols != n) {
                    p->err = gretl_matrix_realloc(ret->v.m, m, n);
                }
            }
        }
    } else {
        /* we need to create a new node */
        ret = newmat(TMP_NODE | AUX_NODE);
        if (ret == NULL) {
            p->err = E_ALLOC;
        }
    }

    if (!p->err && ret->v.m == NULL) {
        if (cmplx) {
            ret->v.m = gretl_cmatrix_new(m, n);
        } else {
            ret->v.m = gretl_matrix_alloc(m, n);
        }
        if (ret->v.m == NULL) {
            p->err = E_ALLOC;
        }
    }

    return ret;
}

static int no_data_error (parser *p)
{
    p->err = E_NODATA;
    return E_NODATA;
}

static NODE *aux_series_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
        no_data_error(p);
        return NULL;
    } else {
        return get_aux_node(p, SERIES, p->dset->n, TMP_NODE);
    }
}

static NODE *aux_empty_series_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
        no_data_error(p);
        return NULL;
    } else {
        return get_aux_node(p, SERIES, 0, TMP_NODE);
    }
}

static NODE *aux_list_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
        no_data_error(p);
        return NULL;
    } else {
        return get_aux_node(p, LIST, 0, TMP_NODE);
    }
}

#define aux_scalar_node(p) get_aux_node(p,NUM,0,0)
#define aux_ivec_node(p,n) get_aux_node(p,IVEC,n,TMP_NODE)
#define aux_matrix_node(p) get_aux_node(p,MAT,0,TMP_NODE)
#define matrix_pointer_node(p) get_aux_node(p,MAT,0,0)
#define aux_mspec_node(p) get_aux_node(p,MSPEC,0,TMP_NODE)
#define aux_string_node(p) get_aux_node(p,STR,0,TMP_NODE)
#define string_pointer_node(p) get_aux_node(p,STR,0,0)
#define aux_bundle_node(p) get_aux_node(p,BUNDLE,0,TMP_NODE)
#define bundle_pointer_node(p) get_aux_node(p,BUNDLE,0,0)
#define aux_array_node(p) get_aux_node(p,ARRAY,0,TMP_NODE)
#define array_pointer_node(p) get_aux_node(p,ARRAY,0,0)
#define aux_parent_node(p) get_aux_node(p,EMPTY,0,0)
#define aux_any_node(p) get_aux_node(p,0,0,0)
#define aux_empty_node(p) get_aux_node(p,EMPTY,0,0)

/* Start of functions that probably should not be needed in
   their present full form, but testing is required before
   they're slimmed down. The general idea is that we should
   already have the user_var pointer we're in need of without
   having to look up it by name, again. However, we'll fall
   back to name look-up (and squawk about it on stderr) if
   need be.
*/

static void *gen_get_lhs_var (parser *p, GretlType type)
{
    void *data = NULL;

    if (p->lh.uv != NULL && p->lh.uv->type == type) {
        data = p->lh.uv->ptr;
    } else {
        if (p->lh.uv == NULL) {
            fprintf(stderr, "*** get: LHS %s '%s' is NULL!\n",
                    gretl_type_get_name(type), p->lh.name);
        } else {
            fprintf(stderr, "*** get: LHS uv '%s' of wrong type!\n",
                    p->lh.name);
        }
        if (type == GRETL_TYPE_BUNDLE) {
            data = get_bundle_by_name(p->lh.name);
        } else if (type == GRETL_TYPE_ARRAY) {
            data = get_array_by_name(p->lh.name);
        } else if (type == GRETL_TYPE_MATRIX) {
            data = get_matrix_by_name(p->lh.name);
        } else if (type == GRETL_TYPE_STRING) {
            data = get_string_by_name(p->lh.name);
        } else if (type == GRETL_TYPE_LIST) {
            data = get_list_by_name(p->lh.name);
        }
    }

    return data;
}

struct typeconv {
    GretlType t;
    int gen_t;
};

struct typeconv conversions[] = {
    { GRETL_TYPE_NONE,   EMPTY },
    { GRETL_TYPE_DOUBLE, NUM },
    { GRETL_TYPE_SERIES, SERIES },
    { GRETL_TYPE_MATRIX, MAT },
    { GRETL_TYPE_LIST,   LIST },
    { GRETL_TYPE_STRING, STR },
    { GRETL_TYPE_BUNDLE, BUNDLE },
    { GRETL_TYPE_ARRAY,  ARRAY }
};

static int gen_type_from_gretl_type (GretlType t)
{
    int i, n = G_N_ELEMENTS(conversions);

    for (i=0; i<n; i++) {
        if (t == conversions[i].t) {
            return conversions[i].gen_t;
        }
    }

    if (t == GRETL_TYPE_INT) {
        return NUM;
    }

    return UNDEF;
}

static GretlType gretl_type_from_gen_type (int gen_t)
{
    int i, n = G_N_ELEMENTS(conversions);

    for (i=0; i<n; i++) {
        if (gen_t == conversions[i].gen_t) {
            return conversions[i].t;
        }
    }

    return GRETL_TYPE_NONE;
}

static int gen_type_is_arrayable (int gen_t)
{
    return gretl_is_arrayable_type(gretl_type_from_gen_type(gen_t));
}

static int type_func_in_parentage (NODE *n)
{
    while (n != NULL) {
        if (undef_arg_ok(n->t)) {
            return 1;
        }
        n = n->parent;
    }

    return 0;
}

static NODE *maybe_rescue_undef_node (NODE *n, parser *p)
{
    int v = current_series_index(p->dset, n->vname);
    user_var *uv = NULL;

    if (v >= 0) {
        n->t = SERIES;
        n->vnum = v;
        n->v.xvec = p->dset->Z[v];
    } else if ((uv = get_user_var_by_name(n->vname)) != NULL) {
        GretlType type = user_var_get_type(uv);

        n->t = gen_type_from_gretl_type(type);
        n->uv = uv;
        if (type == GRETL_TYPE_DOUBLE) {
            n->v.xval = *(double *) uv->ptr;
        } else {
            n->v.ptr = uv->ptr;
        }
    } else {
        undefined_symbol_error(n->vname, p);
    }

    return n;
}

static NODE *revive_undef_node (NODE *n, parser *p)
{
    p->err = 0;
    gretl_error_clear();

    if (n == NULL) {
        n = aux_empty_node(p);
    } else {
        if (n->vname != NULL) {
            free(n->vname);
            n->vname = NULL;
        }
        n->t = EMPTY;
    }

    return n;
}

static int gen_add_or_replace (parser *p, GretlType type, void *data)
{
    int err;

    if (p->lh.uv != NULL) {
        err = user_var_replace_value(p->lh.uv, data, type);
    } else {
        err = user_var_add_or_replace(p->lh.name, type, data);
    }

    return err;
}

static int gen_replace_lhs (parser *p, GretlType type, void *data)
{
    if (p->lh.uv == NULL) {
        fputs("*** gen_replace_lhs: lhs user_var is NULL ***\n", stderr);
        fprintf(stderr, " (type is specified as %s)\n",
                gretl_type_get_name(type));
        return E_DATA;
    } else {
        return user_var_replace_value(p->lh.uv, data, type);
    }
}

static int gen_add_uvar (parser *p, GretlType type, void *data)
{
    int err;

    err = user_var_add(p->lh.name, type, data);

    /* FIXME attach lh.uv pointer? */
    return err;
}

static int gen_edit_list (parser *p, int *list, int op)
{
    user_var *u;
    int err;

    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_LIST) {
        u = p->lh.uv;
    } else {
        if (p->lh.uv == NULL) {
            fprintf(stderr, "*** replace list: LHS uv is NULL!\n");
        } else {
            fprintf(stderr, "*** replace list: LHS uv of wrong type!\n");
        }
        u = get_user_var_of_type_by_name(p->lh.name, GRETL_TYPE_LIST);
    }

    if (op == B_ASN) {
        err = user_list_replace(u, list);
    } else if (op == B_ADD) {
        err = user_list_append(u, list);
    } else {
        /* must be B_SUB */
        err = user_list_subtract(u, list, p->dset);
    }

    return err;
}

static int node_replace_scalar (NODE *n, double x)
{
    int err = 0;

    if (n->uv != NULL && n->uv->type == GRETL_TYPE_DOUBLE) {
        uvar_set_scalar_fast(n->uv, x);
    } else {
        if (n->uv == NULL) {
            fprintf(stderr, "*** node_replace scalar: node uv is NULL!\n");
        } else {
            fprintf(stderr, "*** node_replace scalar: node uv of wrong type!\n");
        }
        err = gretl_scalar_set_value(n->vname, x);
    }

    return err;
}

static int gen_replace_scalar (parser *p, double x)
{
    int err = 0;

    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_DOUBLE) {
        uvar_set_scalar_fast(p->lh.uv, x);
    } else {
        if (p->lh.uv == NULL) {
            fprintf(stderr, "*** gen_replace scalar: LHS uv is NULL!\n");
        } else {
            fprintf(stderr, "***gen_ replace scalar: LHS uv of wrong type!\n");
        }
        err = gretl_scalar_set_value(p->lh.name, x);
    }

    return err;
}

/* end of functions that can probably be slimmed down */

static void eval_warning (parser *p, int op, int errnum)
{
    if (!check_gretl_warning()) {
        const char *w = (op == B_POW)? "pow" : getsymb(op);
        const char *s = (errnum)? gretl_strerror(errnum) : NULL;

        if (s != NULL) {
            gretl_warnmsg_sprintf("%s: %s", w, s);
        } else {
            gretl_warnmsg_set(w);
        }
    }
}

/* evaluation of binary operators (yielding x op y) for
   scalar operands (also increment/decrement operators)
*/

static double xy_calc (double x, double y, int op, int targ, parser *p)
{
    double z = NADBL;

#if EDEBUG > 1
    fprintf(stderr, "xy_calc: x = %g, y = %g, op = %d ('%s')\n",
            x, y, op, getsymb(op));
#endif

    /* assignment */
    if (op == B_ASN || op == B_DOTASN) {
        return y;
    }

    /* testing for presence of NAs? */
    if ((p->flags & P_NATEST) && (na(x) || na(y))) {
        return NADBL;
    }

    /* 0 times anything (even NA) = 0 ? But let's not do this
       for matrices */
    if (targ != MAT && op == B_MUL && (x == 0 || y == 0)) {
        return 0;
    }

    /* logical OR: if x or y is valid and non-zero, ignore NA for
       the other term */
    if (op == B_OR && ((!na(x) && x != 0) || (!na(y) && y != 0))) {
        return 1.0;
    }

    /* logical AND: if either x or y is false, the logical product
       should be false, even if the other term is NA */
    if (op == B_AND && (x == 0 || y == 0)) {
        return 0;
    }

    if (na(x) || na(y)) {
        /* NaN always propagates to the result */
        return NADBL;
    }

    errno = 0;

    switch (op) {
    case B_ADD:
    case INC:
        return x + y;
    case B_SUB:
    case DEC:
        return x - y;
    case B_MUL:
        return x * y;
    case B_DIV:
        return x / y;
    case B_MOD:
        return fmod(x, y);
    case B_AND:
        return x != 0 && y != 0;
    case B_OR:
        return x != 0 || y != 0;
    case B_EQ:
        return x == y;
    case B_NEQ:
        return x != y;
    case B_GT:
        return x > y;
    case B_LT:
        return x < y;
    case B_GTE:
        return x >= y;
    case B_LTE:
        return x <= y;
    case B_POW:
        z = pow(x, y);
        if (errno) {
            eval_warning(p, op, errno);
        }
        return z;
    default:
        return z;
    }
}

static int rmatrix_xy_calc (gretl_matrix *targ,
                            gretl_matrix *src,
                            double x, int xleft,
                            int op, parser *p)
{
    int i, n = targ->rows * targ->cols;

    if (xleft) {
        for (i=0; i<n; i++) {
            targ->val[i] = xy_calc(x, src->val[i], op, MAT, p);
        }
    } else {
        for (i=0; i<n; i++) {
            targ->val[i] = xy_calc(src->val[i], x, op, MAT, p);
        }
    }

    return p->err;
}

static int operator_real_only (int op)
{
    gretl_errmsg_sprintf("'%s': %s", getsymb(op),
                         _("complex operands are not supported"));
    return E_CMPLX;
}

static int function_real_only (int f)
{
    gretl_errmsg_sprintf("%s: %s", getsymb(f),
                         _("complex arguments are not supported"));
    return E_CMPLX;
}

static double complex c_xy_calc (double complex x,
                                 double complex y,
                                 int op, parser *p)
{
    if (op == B_ASN || op == B_DOTASN) {
        return y;
    }

    switch (op) {
    case B_ADD:
        return x + y;
    case B_SUB:
        return x - y;
    case B_MUL:
        return x * y;
    case B_DIV:
        return x / y;
    case B_EQ:
        return x == y;
    case B_NEQ:
        return x != y;
    default:
        p->err = operator_real_only(op);
        return NADBL;
    }
}

static int cmatrix_xy_calc (gretl_matrix *targ,
                            gretl_matrix *src,
                            double complex x,
                            int xleft, int op,
                            parser *p)
{
    int i, n = targ->rows * targ->cols;

    if (xleft) {
        for (i=0; i<n && !p->err; i++) {
            targ->z[i] = c_xy_calc(x, src->z[i], op, p);
        }
    } else {
        for (i=0; i<n && !p->err; i++) {
            targ->z[i] = c_xy_calc(src->z[i], x, op, p);
        }
    }

    return p->err;
}

static double cmatrix_xy_comp (gretl_matrix *m, double x,
                               int op, parser *p)
{
    int i, n = m->rows * m->cols;
    double complex zcond, z = x;
    double ret = 1;

    for (i=0; i<n && ret==1; i++) {
        zcond = c_xy_calc(m->z[i], z, op, p);
        if (p->err) {
            ret = NADBL;
        } else if (zcond == 0) {
            ret = 0;
        }
    }

    return ret;
}

static void n_args_error (int k, int nmin, int nmax, int f, parser *p)
{
    if (nmax == nmin) {
        gretl_errmsg_sprintf( _("Number of arguments (%d) does not "
                                "match the number of\nparameters for "
                                "function %s (%d)"), k, getsymb(f), nmin);
    } else {
        gretl_errmsg_sprintf( _("Number of arguments (%d) does not "
                                "match the number of\nparameters for "
                                "function %s (%d to %d)"), k, getsymb(f),
                              nmin, nmax);
    }
    p->err = 1;
}

#define randgen(f) (f == F_RANDGEN || f == F_MRANDGEN || f == F_RANDGEN1)

static int check_dist_count (int d, int f, int *np, int *argc)
{
    int err = 0;

    *np = *argc = 0;

    if (d == D_NC_T || d == D_NC_F || d == D_NC_CHISQ) {
        /* non-central t, chisq and F: only CDF, PDF and INVCDF supported */
        if (f == F_PDF || f == F_CDF || f == F_INVCDF) {
            *np = (d == D_NC_F) ? 3 : 2;
        } else {
            err = E_INVARG;
        }
    } else if (d == D_UNIFORM || d == D_UDISCRT) {
        /* only randgen is supported */
        if (randgen(f)) {
            *np = 2; /* min, max */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_NORMAL) {
        /* all functions supported */
        if (randgen(f)) {
            *np = 2; /* mu, sigma */
        } else {
            *np = 0; /* N(0,1) is assumed */
        }
    } else if (d == D_STUDENT) {
        /* Student t: all functions supported */
        *np = 1; /* df */
    } else if (d == D_CHISQ) {
        /* chi-square: all functions supported */
        *np = 1; /* df */
    } else if (d == D_SNEDECOR) {
        /* all functions supported */
        *np = 2; /* dfn, dfd */
    } else if (d == D_GAMMA) {
        /* partial support */
        if (f == F_CRIT) {
            err = 1;
        } else {
            *np = 2; /* shape, scale */
        }
    } else if (d == D_BINOMIAL) {
        *np = 2; /* prob, trials */
    } else if (d == D_BINORM) {
        /* bivariate normal: cdf only */
        if (f == F_CDF) {
            *np = 1; /* rho */
            *argc = 2; /* note: special */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_POISSON) {
        *np = 1;
    } else if (d == D_EXPON) {
        /* inverse cdf not supported */
        if (f == F_INVCDF) {
            err = E_INVARG;
        } else {
            *np = 1; /* scale */
        }
    } else if (d == D_WEIBULL) {
        /* inverse cdf not supported */
        if (f == F_INVCDF) {
            err = E_INVARG;
        } else {
            *np = 2; /* shape, scale */
        }
    } else if (d == D_GED) {
        *np = 1; /* shape */
    } else if (d == D_LAPLACE) {
        *np = 2; /* mean, scale */
    } else if (d == D_DW) {
        /* Durbin-Watson: only critical value */
        if (f == F_CRIT) {
            *np = 2; /* n, k */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_JOHANSEN) {
        /* Johansen trace test: only p-value */
        if (f == F_PVAL) {
            *np = 3;
        } else {
            err = E_INVARG;
        }
    } else if (d == D_BETA) {
        /* cdf, pdf, randgen only */
        if (f == F_CDF || f == F_PDF || randgen(f)) {
            *np = 2; /* shape1, shape2 */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_BETABIN) {
        /* randgen only */
        if (randgen(f)) {
            *np = 3; /* n, shape1, shape2 */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_LOGISTIC) {
        if (randgen(f)) {
            *np = 2; /* location, scale */
        } else if (f == F_CDF) {
            *np = 0; /* (0,1) assumed */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_DIRICHLET) {
        /* mrandgen only */
        if (randgen(f)) {
            *np = 1; /* alpha vector */
        } else {
            err = E_INVARG;
        }
    } else if (d == D_DISCRETE) {
        /* randgen and mrandgen only */
        if (randgen(f)) {
            *np = 1; /* prob vector */
        } else {
            err = E_INVARG;
        }
    } else {
        err = E_INVARG;
    }

    if (!err && !randgen(f) && *argc == 0) {
        *argc = 1;
    }

    return err;
}

static double scalar_pdist (int t, int d, const double *parm,
                            int np, double arg, parser *p)
{
    double x = NADBL;
    int i;

    for (i=0; i<np; i++) {
        if (na(parm[i])) {
            return NADBL;
        }
    }

    if (t == F_PVAL) {
        x = gretl_get_pvalue(d, parm, arg);
    } else if (t == F_PDF) {
        x = gretl_get_pdf(d, parm, arg);
    } else if (t == F_CDF) {
        x = gretl_get_cdf(d, parm, arg);
    } else if (t == F_INVCDF) {
        x = gretl_get_cdf_inverse(d, parm, arg);
    } else if (t == F_CRIT) {
        x = gretl_get_critval(d, parm, arg);
    } else {
        p->err = E_PARSE;
    }

    return x;
}

/* @parm contains an array of scalar parameters;
   @argvec contains a series of argument values.
*/

static int series_pdist (double *x, int f, int d,
                         double *parm, int np,
                         const double *argvec,
                         parser *p)
{
    int t;

    if (f == F_PDF) {
        /* fast treatment, for pdf only at this point */
        int n = sample_size(p->dset);

        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            x[t] = argvec[t];
        }
        gretl_fill_pdf_array(d, parm, x + p->dset->t1, n);
    } else {
        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            x[t] = scalar_pdist(f, d, parm, np, argvec[t], p);
        }
    }

    return 0;
}

/* @parm contains an array of zero to two scalar parameters;
   @argmat contains an array of argument values.
*/

static gretl_matrix *matrix_pdist (int f, int d,
                                   double *parm, int np,
                                   gretl_matrix *argmat,
                                   parser *p)
{
    gretl_matrix *m;
    double x;
    int i, n;

    if (gretl_is_null_matrix(argmat)) {
        return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(argmat->rows, argmat->cols);
    if (m == NULL) {
        p->err = E_ALLOC;
        return NULL;
    }

    n = m->rows * m->cols;

    for (i=0; i<n && !p->err; i++) {
        x = scalar_pdist(f, d, parm, np, argmat->val[i], p);
        if (na(x)) {
            p->err = E_MISSDATA;
        } else {
            m->val[i] = x;
        }
    }

    if (p->err) {
        gretl_matrix_free(m);
        m = NULL;
    }

    return m;
}

/* Gets a matrix from a node of type MAT or NUM. In the
   latter case it's a static matrix, good for use in
   calculation, but it should NOT be passed on or
   freed.
*/

#define N_STATIC 8

static gretl_matrix *node_get_matrix (NODE *n, parser *p,
                                      int i, int argnum)
{
    if (p->err) {
        /* don't compound prior error */
        return NULL;
    } else if (n->t == MAT) {
        return n->v.m;
    } else if (n->t != NUM) {
        if (argnum > 0) {
            gretl_errmsg_sprintf(_("arg %d is missing or of invalid type"),
                                 argnum);
        }
        p->err = E_INVARG;
        return NULL;
    } else if (i < 0 || i >= N_STATIC) {
        fprintf(stderr, "node_get_matrix: invalid static matrix index\n");
        p->err = E_DATA;
        return NULL;
    } else {
        /* "fake" a 1 x 1 matrix */
        static gretl_matrix mm[N_STATIC];
        static double mval[N_STATIC];
        static int started;
        gretl_matrix *ret;

        if (!started) {
            int j;

            for (j=0; j<N_STATIC; j++) {
                gretl_matrix_init(&mm[j]);
                mm[j].rows = mm[j].cols = 1;
            }
            started = 1;
        }
        mval[i] = n->v.xval;
        ret = &mm[i];
        ret->val = &mval[i];

        return ret;
    }
}

static gretl_matrix *node_get_real_matrix (NODE *n, parser *p,
                                           int i, int argnum)
{
    gretl_matrix *m = node_get_matrix(n, p, i, argnum);

    if (!p->err && m->is_complex) {
        p->err = E_CMPLX;
        return NULL;
    } else {
        return m;
    }
}

static double node_get_scalar (NODE *n, parser *p)
{
    if (n->t == NUM) {
        return n->v.xval;
    } else if (scalar_matrix_node(n)) {
        return n->v.m->val[0];
    } else if (n->t == EMPTY) {
        p->err = E_ARGS;
        return NADBL;
    } else {
        p->err = E_INVARG;
        return NADBL;
    }
}

static int node_get_int (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    if (p->err) {
        return -1;
    } else {
        return gretl_int_from_double(x, &p->err);
    }
}

static guint32 node_get_guint32 (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    if (p->err) {
        return 0;
    } else {
        return gretl_unsigned_from_double(x, &p->err);
    }
}

static int node_get_bool (NODE *n, parser *p, int deflt)
{
    int ret = -1;

    if (!null_node(n)) {
        int k = node_get_int(n, p);

        if (!p->err) {
            ret = (k != 0);
        }
    } else if (deflt == 0 || deflt == 1) {
        ret = deflt;
    } else {
        p->err = E_ARGS;
    }

    return ret;
}

static NODE *dec2bin_node (NODE *n, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (!p->err) {
        gretl_matrix *v = NULL;
        double x = NADBL;

        if (n->t == NUM) {
            x = n->v.xval;
        } else {
            v = n->v.m;
        }
        if (!p->err) {
            ret->v.m = dec2bin(x, v, &p->err);
        }
    }

    return ret;
}

static NODE *bin2dec_node (NODE *n, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (!p->err) {
        ret->v.m = bin2dec(n->v.m, &p->err);
    }

    return ret;
}

static NODE *access_series_node (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (!p->err) {
        int v0 = caller_series_index(p->dset, n->v.str);

        if (v0 > 0) {
            int t, v1 = p->dset->v;

            p->err = dataset_add_series(p->dset, 1);
            if (!p->err) {
                for (t=0; t<p->dset->n; t++) {
                    p->dset->Z[v1][t] = p->dset->Z[v0][t];
                }
                strcpy(p->dset->varname[v1], n->v.str);
                ret->v.ivec = gretl_list_new(1);
                ret->v.ivec[1] = v1;
            }
        } else {
            ret->v.ivec = gretl_list_new(0);
        }
    }

    return ret;
}

static NODE *DW_node (NODE *r, parser *p)
{
    NODE *ret = NULL;
    int i, parm[2] = {0};

    for (i=0; i<2 && !p->err; i++) {
        parm[i] = node_get_int(r->v.bn.n[i+1], p);
    }

    if (!p->err && (parm[0] < 6 || parm[1] < 0)) {
        p->err = E_INVARG;
    }

    if (!p->err) {
        ret = aux_matrix_node(p);
        if (ret != NULL) {
            ret->v.m = gretl_get_DW(parm[0], parm[1], &p->err);
        }
    }

    return ret;
}

static NODE *eval_urcpval (NODE *n, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        NODE *e;
        int i, m = r->v.bn.n_nodes;
        int iargs[3] = {0};
        double tau = NADBL;

        if (m != 4) {
            n_args_error(m, 4, 4, n->t, p);
        }

        /* need double, int, int, int */
        for (i=0; i<4 && !p->err; i++) {
            e = r->v.bn.n[i];
            if (scalar_node(e)) {
                if (i == 0) {
                    tau = node_get_scalar(e, p);
                } else {
                    iargs[i-1] = node_get_int(e, p);
                }
            } else {
                p->err = e_types(e);
            }
        }

        if (!p->err) {
            int nobs = iargs[0];
            int niv = iargs[1];
            int itv = iargs[2];

            ret = aux_scalar_node(p);
            if (ret != NULL) {
                ret->v.xval = get_urc_pvalue(tau, nobs, niv, itv);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *binperm_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        int n = node_get_int(l, p);
        int k = node_get_int(r, p);

        if (!p->err) {
            ret->v.m = bit_permutations(n, k, &p->err);
        }
    }

    return ret;
}

static int get_matrix_size (gretl_matrix *a, gretl_matrix *b,
                            int *r, int *c)
{
    int err = 0;

    /* if both matrices are present, they must be the
       same size */

    if (a != NULL) {
        *r = a->rows;
        *c = b->cols;
        if (b != NULL && (b->rows != *r || b->cols != *c)) {
            err = E_NONCONF;
        }
    } else if (b != NULL) {
        *r = b->rows;
        *c = b->cols;
    } else {
        *r = *c = 0;
    }

    return err;
}

static NODE *bvnorm_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        double *avec = NULL, *bvec = NULL;
        gretl_matrix *amat = NULL, *bmat = NULL;
        double a, b, args[2];
        double rho = NADBL;
        NODE *e;
        int i, mode = 0;

        for (i=0; i<3 && !p->err; i++) {
            e = n->v.bn.n[i+1];
            if (scalar_node(e)) {
                if (i == 0) {
                    rho = node_get_scalar(e, p);
                } else {
                    args[i-1] = node_get_scalar(e, p);
                }
            } else if (i == 1) {
                if (e->t == SERIES) {
                    avec = e->v.xvec;
                } else if (e->t == MAT) {
                    amat = e->v.m;
                }
            } else if (i == 2) {
                if (e->t == SERIES) {
                    bvec = e->v.xvec;
                } else if (e->t == MAT) {
                    bmat = e->v.m;
                }
            } else {
                node_type_error(F_CDF, i+1, NUM, e, p);
            }
        }

        if (!p->err) {
            if ((avec != NULL && bmat != NULL) ||
                (bvec != NULL && amat != NULL)) {
                p->err = E_INVARG;
            } else if (avec != NULL || bvec != NULL) {
                mode = 1;
                ret = aux_series_node(p);
            } else if (amat != NULL || bmat != NULL) {
                mode = 2;
                ret = aux_matrix_node(p);
            } else {
                mode = 0;
                ret = aux_scalar_node(p);
            }
        }

        if (p->err) {
            return ret;
        }

        if (mode == 0) {
            /* a, b are both scalars */
            ret->v.xval = bvnorm_cdf(rho, args[0], args[1]);
        } else if (mode == 1) {
            /* a and/or b are series */
            int t;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                a = (avec != NULL)? avec[t] : args[0];
                b = (bvec != NULL)? bvec[t] : args[1];
                if (na(a) || na(b)) {
                    ret->v.xvec[t] = NADBL;
                } else {
                    ret->v.xvec[t] = bvnorm_cdf(rho, a, b);
                }
            }
        } else if (mode == 2) {
            /* a and/or b are matrices */
            gretl_matrix *m = NULL;
            int r, c;

            p->err = get_matrix_size(amat, bmat, &r, &c);

            if (!p->err && r > 0 && c > 0) {
                m = gretl_matrix_alloc(r, c);
            }
            if (m != NULL) {
                int i, n = r * c;

                for (i=0; i<n && !p->err; i++) {
                    a = (amat != NULL)? amat->val[i] : args[0];
                    b = (bmat != NULL)? bmat->val[i] : args[1];
                    m->val[i] = bvnorm_cdf(rho, a, b);
                    if (na(m->val[i])) {
                        /* matrix: change NAs to NaNs */
                        m->val[i] = 0.0/0.0;
                    }
                }
            }
            if (ret->v.m != NULL) {
                gretl_matrix_free(ret->v.m);
            }
            ret->v.m = m;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *dirichlet_node (NODE *args, parser *p)
{
    NODE *ret = NULL;
    gretl_matrix *alpha = NULL;
    int n = 0;

    if (args->v.bn.n[1]->t == MAT) {
        alpha = args->v.bn.n[1]->v.m;
    }
    if (args->v.bn.n[2]->t == NUM) {
        n = node_get_int(args->v.bn.n[2], p);
    }
    if (alpha == NULL || n <= 0) {
        p->err = E_INVARG;
    } else {
        ret = aux_matrix_node(p);
        if (!p->err) {
            ret->v.m = gretl_rand_dirichlet(alpha, n, &p->err);
        }
    }

    return ret;
}

static NODE *discrete_rand_node (NODE *n, NODE *args, parser *p)
{
    NODE *ret = NULL;
    gretl_matrix *probs = NULL;
    int vsize = 0;

    if (args->v.bn.n[1]->t == MAT) {
        probs = args->v.bn.n[1]->v.m;
        if (gretl_vector_get_length(probs) == 0) {
            p->err = E_INVARG;
            return NULL;
        }
    }

    if (n->t == F_RANDGEN1) {
        ret = aux_scalar_node(p);
    } else if (n->t == F_RANDGEN) {
        ret = aux_series_node(p);
    } else {
        /* F_MRANDGEN */
        int r = node_get_int(args->v.bn.n[2], p);
        int c = node_get_int(args->v.bn.n[3], p);

        if (!p->err && (r <= 0 || c <= 0)) {
            p->err = E_INVARG;
        } else {
            ret = aux_sized_matrix_node(p, r, c, 0);
            vsize = r * c;
        }
    }

    if (!p->err && n->t == F_RANDGEN1) {
        p->err = gretl_rand_discrete(&ret->v.xval, 0, 0, probs);
    } else if (!p->err && n->t == F_RANDGEN) {
        p->err = gretl_rand_discrete(ret->v.xvec,
                                     p->dset->t1, p->dset->t2,
                                     probs);
    } else if (!p->err) {
        /* F_MRANDGEN */
        p->err = gretl_rand_discrete(ret->v.m->val, 0, vsize-1, probs);
    }

    return ret;
}

/* Return a node containing the evaluated result of a probability
   distribution function. Here, node @n is the function node and node
   @r holds the arguments, the number of which ir not known in
   advance.
*/

static NODE *eval_pdist (NODE *n, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        NODE *e, *s;
        int i, k, nr = r->v.bn.n_nodes;
        int rgen = (n->t == F_RANDGEN);
        int mrgen = (n->t == F_MRANDGEN);
        int rgen1 = (n->t == F_RANDGEN1);
        double parm[3] = {0};
        double argval = NADBL;
        double *parmvec[2] = { NULL };
        double *argvec = NULL;
        int pvlen[2] = {0};
        gretl_matrix *argmat = NULL;
        int rows = 0, cols = 0;
        int np, argc, bb, d = 0;

        if (nr < 2) {
            p->err = E_ARGS;
            goto disterr;
        }

        /* the first argument must hold a string identifying
           the distribution to be employed
        */
        s = r->v.bn.n[0];
        if (s->t == STR) {
            char *dstr = s->v.str;

            d = dist_code_from_string(dstr);
            if (d == 0) {
                dstr = get_string_by_name(dstr);
                if (dstr != NULL) {
                    d = dist_code_from_string(dstr);
                }
            }
            if (d == 0) {
                p->err = E_INVARG;
            }
        } else {
            node_type_error(n->t, 0, STR, s, p);
        }

        if (mrgen) {
            if (d == D_DIRICHLET) {
                if (nr != 3) {
                    n_args_error(nr, 3, 3, n->t, p);
                }
            } else if (d == D_DISCRETE) {
                if (nr != 4) {
                    n_args_error(nr, 4, 4, n->t, p);
                }
            } else if (nr < 4 || nr > 7) {
                n_args_error(nr, 4, 7, n->t, p);
            }
        } else if (nr > 5) {
            n_args_error(nr, 2, 5, n->t, p);
        } else {
            if (d == D_DIRICHLET) {
                /* mrandgen only */
                p->err = E_INVARG;
            }
        }

        if (p->err) {
            goto disterr;
        } else {
            p->err = check_dist_count(d, n->t, &np, &argc);
        }
        if (!p->err && d != D_DIRICHLET) {
            k = np + argc + 2 * mrgen;
            if (k != nr - 1) {
                p->err = E_INVARG;
            }
        }
        if (p->err) {
            goto disterr;
        }

        bb = (d == D_BETABIN);

        if (d == D_DW) {
            /* special: Durbin-Watson */
            return DW_node(r, p);
        } else if (d == D_BINORM) {
            /* special: bivariate normal */
            return bvnorm_node(r, p);
        } else if (d == D_DIRICHLET) {
            /* special: dirichlet */
            return dirichlet_node(r, p);
        } else if (d == D_DISCRETE) {
            return discrete_rand_node(n, r, p);
        }

        for (i=1; i<=k && !p->err; i++) {
            e = r->v.bn.n[i];
            if (scalar_node(e)) {
                /* scalars always acceptable */
                if (mrgen) {
                    if (i == k) {
                        cols = node_get_int(e, p);
                    } else if (i == k-1) {
                        rows = node_get_int(e, p);
                    } else {
                        parm[i-1] = node_get_scalar(e, p);
                    }
                } else if (i == k && argc > 0) {
                    argval = node_get_scalar(e, p);
                } else {
                    parm[i-1] = node_get_scalar(e, p);
                }
            } else if (i == k && e->t == SERIES) {
                /* a series in the last place? */
                if (bb) {
                    node_type_error(n->t, i, NUM, e, p);
                } else if (rgen) {
                    parmvec[i-1] = e->v.xvec;
                } else if (mrgen) {
                    node_type_error(n->t, i, NUM, e, p);
                } else {
                    argvec = e->v.xvec;
                }
            } else if (i == k && e->t == MAT) {
                /* a matrix in the last place? */
                if (mrgen) {
                    parmvec[i-1] = e->v.m->val;
                    pvlen[i-1] = e->v.m->rows * e->v.m->cols;
                } else if (rgen) {
                    node_type_error(n->t, i, NUM, e, p);
                } else {
                    argmat = e->v.m;
                }
            } else if (e->t == SERIES) {
                /* a series param for randgen? */
                if (rgen && !bb) {
                    parmvec[i-1] = e->v.xvec;
                } else {
                    node_type_error(n->t, i, NUM, e, p);
                }
            } else if (e->t == MAT) {
                /* a matrix param for mrandgen? */
                if (mrgen && !bb) {
                    parmvec[i-1] = e->v.m->val;
                    pvlen[i-1] = e->v.m->rows * e->v.m->cols;
                } else {
                    node_type_error(n->t, i, NUM, e, p);
                }
            } else {
                p->err = E_INVARG;
                fprintf(stderr, "eval_pdist: arg %d, bad type %d\n", i+1, e->t);
            }
        }

        if (mrgen) {
            int rlen = rows * cols;

            if ((parmvec[0] != NULL && pvlen[0] != rlen) ||
                (parmvec[1] != NULL && pvlen[1] != rlen)) {
                p->err = E_NONCONF;
            }
        }

        if (p->err) {
            goto disterr;
        }

        if (mrgen) {
            ret = aux_matrix_node(p);
        } else if (rgen || argvec != NULL) {
            ret = aux_series_node(p);
        } else if (argmat != NULL) {
            ret = aux_matrix_node(p);
        } else {
            ret = aux_scalar_node(p);
        }

        if (ret == NULL) {
            goto disterr;
        }

        if (rgen) {
            p->err = gretl_fill_random_series(ret->v.xvec, d, parm,
                                              parmvec[0], parmvec[1],
                                              p->dset);
        } else if (mrgen) {
            ret->v.m = gretl_get_random_matrix(d, parm,
                                               parmvec[0], parmvec[1],
                                               rows, cols,
                                               &p->err);
        } else if (rgen1) {
            ret->v.xval = gretl_get_random_scalar(d, parm, &p->err);
        } else if (argvec != NULL) {
            p->err = series_pdist(ret->v.xvec, n->t, d, parm, np,
                                  argvec, p);
        } else if (argmat != NULL) {
            ret->v.m = matrix_pdist(n->t, d, parm, np, argmat, p);
        } else {
            ret->v.xval = scalar_pdist(n->t, d, parm, np, argval, p);
        }
    } else {
        ret = aux_any_node(p);
    }

  disterr:

    return ret;
}

static int mpi_rank = -1;
static int mpi_size = 0;

void set_mpi_rank_and_size (int rank, int size)
{
    mpi_rank = rank;
    mpi_size = size;
}

static double get_const_by_id (int id)
{
    if (id == CONST_PI) {
        return M_PI;
    } else if (id == CONST_EPS) {
        /* see https://en.wikipedia.org/wiki/Machine_epsilon :
           we now use the (b) definition, as per Matlab, Gauss,
           R and others
        */
        return pow(2.0, -52);
    } else if (id == CONST_INF) {
#ifdef INFINITY
        return INFINITY;
#else
        return 1.0/0.0;
#endif
    } else if (id == CONST_NAN) {
#ifdef NAN
        return NAN;
#else
        return 0.0/0.0;
#endif
    } else if (id == CONST_WIN32) {
#ifdef WIN32
        return 1;
#else
        return 0;
#endif
    } else if (id == CONST_HAVE_MPI) {
#ifdef HAVE_MPI
        return check_for_mpiexec();
#else
        return 0;
#endif
    } else if (id == CONST_MPI_RANK) {
        return mpi_rank;
    } else if (id == CONST_MPI_SIZE) {
        return mpi_size;
    } else if (id == CONST_N_PROC) {
        return gretl_n_processors();
    } else if (id == CONST_TRUE) {
        return 1;
    } else if (id == CONST_FALSE) {
        return 0;
    } else {
        return NADBL;
    }
}

/* look up and return numerical values of symbolic constants */

static NODE *retrieve_const (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        ret->v.xval = get_const_by_id(n->v.idnum);
    }

    return ret;
}

double get_const_by_name (const char *name, int *err)
{
    int id = const_lookup(name);

    if (id > 0) {
        return get_const_by_id(id);
    } else {
        if (err != NULL) {
            *err = E_DATA;
        }
        return NADBL;
    }
}

#ifdef HAVE_MPI

#include "genmpi.c"

#else

static NODE *mpi_transfer_node (NODE *l, NODE *r, NODE *r2,
                                int f, parser *p)
{
    gretl_errmsg_set(_("MPI is not supported in this gretl build"));
    p->err = 1;
    return NULL;
}

static NODE *mpi_barrier_node (parser *p)
{
    gretl_errmsg_set(_("MPI is not supported in this gretl build"));
    p->err = 1;
    return NULL;
}

#endif /* !HAVE_MPI */

static NODE *scalar_calc (NODE *x, NODE *y, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        ret->v.xval = xy_calc(x->v.xval, y->v.xval, f, NUM, p);
    }

    return ret;
}

/* support addition of integer meaning subset of string */

static NODE *string_offset (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        int n = g_utf8_strlen(l->v.str, -1);
        int k = r->v.xval;

        if (k < 0) {
            p->err = E_DATA;
        } else if (k >= n) {
            ret->v.str = gretl_strdup("");
        } else {
            char *p = g_utf8_offset_to_pointer(l->v.str, k);

            ret->v.str = gretl_strdup(p);
        }
        if (!p->err && ret->v.str == NULL) {
            p->err = E_ALLOC;
        }
    }

    return ret;
}

/* support multiplication by integer meaning repetition of a string */

static NODE *string_mult (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        const char *s = NULL;
        int k = 0;

        if (l->t == STR && r->t == NUM) {
            s = l->v.str;
            k = node_get_int(r, p);
        } else if (l->t == NUM && r->t == STR) {
            s = r->v.str;
            k = node_get_int(l, p);
        } else {
            p->err = E_TYPES;
        }
        if (!p->err && k < 0) {
            p->err = E_INVARG;
        }
        if (!p->err && k == 0) {
            ret->v.str = gretl_strdup("");
        } else if (!p->err) {
            int i, n = strlen(s);

            ret->v.str = calloc(k * n + 1, 1);
            for (i=0; i<k; i++) {
                strcat(ret->v.str, s);
            }
        }
    }

    return ret;
}

static NODE *compare_strings (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        int s = strcmp(l->v.str, r->v.str);

        ret->v.xval = (f == B_EQ)? (s == 0) : (s != 0);
    }

    return ret;
}

/*
   We're looking at a comparison, with either a series on the left and
   a string on the right or vice versa.  This can work if the series
   in question is string-valued, as in

     Case 1: series foo = (x == "strval")

   It can also work if the string is an observation marker, as in

     Case 2: series foo = (obs >= "CA")
*/

static NODE *series_string_calc (NODE *l, NODE *r, int f, parser *p)
{
    double xt = NADBL, yt = NADBL;
    double *x = NULL, *y = NULL;
    double *alt;
    const char *strval;
    int vnum, t;
    NODE *ret;

    if (r->t == STR) {
        strval = r->v.str;
        vnum = l->vnum;
        x = l->v.xvec;
        alt = &yt;
    } else {
        strval = l->v.str;
        vnum = r->vnum;
        y = r->v.xvec;
        alt = &xt;
    }

    ret = aux_series_node(p);
    if (p->err) {
        return ret;
    }

    if (vnum > 0 && is_string_valued(p->dset, vnum)) {
        /* we must be in Case 1 */
        *alt = series_decode_string(p->dset, vnum, strval);
        if (na(*alt)) {
            /* @strval is not a string value of the given series */
            double xval = (f == B_EQ)? 0 : (f == B_NEQ)? 1 : NADBL;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                ret->v.xvec[t] = xval;
            }
            return ret; /* NA case handled */
        }
    } else {
        /* try interpreting @strval as an observation string */
        if (annual_data(p->dset)) {
            *alt = get_date_x(p->dset->pd, strval);
        } else {
            t = dateton(strval, p->dset);
            if (t >= 0) {
                *alt = t + 1;
            }
        }
    }

    if (na(*alt)) {
        gretl_errmsg_sprintf(_("got invalid field '%s'"), strval);
        p->err = E_TYPES;
        return NULL;
    }

    if (ret != NULL) {
        int t1 = autoreg(p) ? p->obs : p->dset->t1;
        int t2 = autoreg(p) ? p->obs : p->dset->t2;

        for (t=t1; t<=t2; t++) {
            if (x != NULL) {
                xt = x[t];
            } else if (y != NULL) {
                yt = y[t];
            }
            ret->v.xvec[t] = xy_calc(xt, yt, f, SERIES, p);
        }
    }

    return ret;
}

static double *list_node_get_series (NODE *n, parser *p)
{
    if (n->v.ivec[0] == 1) {
        int v = n->v.ivec[1];

        if (v >= 0 && v < p->dset->v) {
            return p->dset->Z[v];
        }
    }

    p->err = E_INVARG;
    return NULL;
}

/* At least one of the nodes is a series; the other may be a
   scalar or 1 x 1 matrix */

static NODE *series_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);
    const double *x = NULL, *y = NULL;
    double xt = 0, yt = 0;
    int tmax = p->dset->t2;

    if (ret == NULL) {
        return NULL;
    }

    if (p->dset->n > p->dset_n) {
        /* can arise when stack() is in the tree ->
           dataset gets extended on the fly
        */
        if (l->t == SERIES) {
            if (useries_node(l)) {
                l->v.xvec = p->dset->Z[l->vnum];
            } else {
                tmax = MIN(tmax, p->dset_n - 1);
            }
        }
        if (r->t == SERIES) {
            if (useries_node(r)) {
                r->v.xvec = p->dset->Z[r->vnum];
            } else {
                tmax = MIN(tmax, p->dset_n - 1);
            }
        }
    }

    if (l->t == SERIES) {
        x = l->v.xvec;
    } else if (l->t == LIST) {
        x = list_node_get_series(l, p);
    } else if (l->t == NUM) {
        xt = l->v.xval;
    } else if (l->t == MAT) {
        xt = l->v.m->val[0];
    }

    if (r->t == SERIES) {
        y = r->v.xvec;
    } else if (r->t == LIST) {
        y = list_node_get_series(r, p);
    } else if (r->t == NUM) {
        yt = r->v.xval;
    } else if (r->t == MAT) {
        yt = r->v.m->val[0];
    }

    if (!p->err) {
        int t1 = autoreg(p) ? p->obs : p->dset->t1;
        int t2 = autoreg(p) ? p->obs : tmax;
        int t;

        for (t=t1; t<=t2; t++) {
            if (x != NULL) {
                xt = x[t];
            }
            if (y != NULL) {
                yt = y[t];
            }
            ret->v.xvec[t] = xy_calc(xt, yt, f, SERIES, p);
        }
    }

    return ret;
}

/* Try to determine if it's OK for a function to return a
   string-valued series, given that such a series must be a member of
   a dataset. The node @f should hold the function in question. In
   case the node that's to hand in the caller is actually a top-level
   argument of the function, one can pass its "parent" pointer as @f.
   Returns an error code if this operation is not OK.
*/

static int strvals_return_check (NODE *f, parser *p)
{
    if (f != p->tree || (p->targ != SERIES && p->targ != UNK)) {
        /* we must be at the top of the tree, with the target of
           assignment either known to be a series or undetermined
           as yet
        */
        gretl_errmsg_sprintf(_("%s: series usage is valid only in direct assignment"),
                             getsymb(f->t));
        return E_DATA;
    } else if (dataset_is_subsampled(p->dset) && p->lh.t == SERIES &&
               !is_string_valued(p->dset, p->lh.vnum)) {
        /* can't assign string values to just part of a numeric series */
        gretl_errmsg_set(_("Sorry, can't do this with a sub-sampled dataset"));
        return E_DATA;
    } else {
        /* OK, we'll try it */
        return 0;
    }
}

/* Create and store a series table to be attached to a string-valued
   series return.
*/

static void prep_strvals_return (char **S, int ns, parser *p)
{
    series_table *st = series_table_new(S, ns, &p->err);

    if (!p->err) {
        p->lh.stab = st;
    } else {
        strings_array_free(S, ns);
    }
}

/* Both nodes are string-valued series. We support a limited
   set of operations.
*/

static NODE *strvals_calc (NODE *l, NODE *r, NODE *n, parser *p)
{
    NODE *ret = NULL;
    const char *sl, *sr;
    char **Sx = NULL;
    int nr = 0, nx = 0;
    int vl, vr, f = n->t;
    int i, t, eq;

    if (f == B_POW) {
        p->err = strvals_return_check(n, p);
    } else if (f != B_EQ && f != B_NEQ) {
        p->err = E_TYPES;
    }
    if (p->err) {
        return NULL;
    }

    ret = aux_series_node(p);
    if (ret == NULL) {
        return NULL;
    }

    vl = l->vnum;
    vr = r->vnum;

    if (f == B_POW) {
        /* "logical product" */
        char *slr, **Sl, **Sr;
        int nl, j, ll;

        Sl = series_get_string_vals(p->dset, vl, &nl, 1);
        Sr = series_get_string_vals(p->dset, vr, &nr, 1);
        nx = nl * nr;

        Sx = strings_array_new(nx);
        if (Sx == NULL) {
            p->err = E_ALLOC;
            return ret;
        }

        for (i=0; i<nl; i++) {
            ll = strlen(Sl[i]) + 2;
            for (j=0; j<nr; j++) {
                slr = calloc(ll + strlen(Sr[j]), 1);
                sprintf(slr, "%s.%s", Sl[i], Sr[j]);
                Sx[i*nr+j] = slr;
            }
        }
    }

    for (t=p->dset->t1; t<=p->dset->t2; t++) {
        sl = series_get_string_for_obs(p->dset, vl, t);
        sr = series_get_string_for_obs(p->dset, vr, t);
        if (sl == NULL || sr == NULL) {
            ret->v.xvec[t] = NADBL;
        } else if (f == B_POW) {
            int il = p->dset->Z[vl][t] - 1;
            int ir = p->dset->Z[vr][t] - 1;

            ret->v.xvec[t] = il*nr + ir + 1;
        } else {
            eq = strcmp(sl, sr) == 0;
            ret->v.xvec[t] = (f == B_EQ)? eq : !eq;
        }
    }

    if (f == B_POW && Sx != NULL) {
        prep_strvals_return(Sx, nx, p);
    }

    return ret;
}

static int op_symbol (int op)
{
    switch (op) {
    case B_DOTMULT: return '*';
    case B_DOTDIV:  return '/';
    case B_DOTPOW:  return '^';
    case B_DOTADD:  return '+';
    case B_DOTSUB:  return '-';
    case B_DOTEQ:   return '=';
    case B_DOTGT:   return '>';
    case B_DOTLT:   return '<';
    case B_DOTGTE:  return ']';
    case B_DOTLTE:  return '[';
    case B_DOTNEQ:  return '!';
    default: return 0;
    }
}

static gretl_matrix *nullmat_multiply (const gretl_matrix *A,
                                       const gretl_matrix *B,
                                       int op, int *err)
{
    gretl_matrix *C = NULL;

    if (A->rows == 0 && A->cols == 0 &&
        B->rows == 0 && B->cols == 0) {
        C = gretl_null_matrix_new();
    } else {
        int Lc = op == B_TRMUL ? A->rows : A->cols;
        int Cr = op == B_TRMUL ? A->cols : A->rows;
        int Cc = B->cols;

        if (Lc != B->rows) {
            *err = E_NONCONF;
        } else {
            if (Cr > 0 && Cc > 0) {
                C = gretl_zero_matrix_new(Cr, Cc);
            } else {
                C = gretl_matrix_alloc(Cr, Cc);
            }
            if (C == NULL) {
                *err = E_ALLOC;
            }
        }
    }

    return C;
}

static gretl_matrix *matrix_add_sub_scalar (const gretl_matrix *A,
                                            const gretl_matrix *B,
                                            int op)
{
    gretl_matrix *C;
    double xval, *xvec;
    int r, c;

    if (gretl_matrix_is_scalar(A)) {
        r = B->rows;
        c = B->cols;
        xval = A->val[0];
        xvec = B->val;
    } else {
        r = A->rows;
        c = A->cols;
        xval = B->val[0];
        xvec = A->val;
    }

    C = gretl_matrix_alloc(r, c);

    if (C != NULL) {
        int i, n = r * c;

        if (op == B_ADD) {
            for (i=0; i<n; i++) {
                C->val[i] = xvec[i] + xval;
            }
        } else {
            if (xvec == A->val) {
                for (i=0; i<n; i++) {
                    C->val[i] = xvec[i] - xval;
                }
            } else {
                for (i=0; i<n; i++) {
                    C->val[i] = xval - xvec[i];
                }
            }
        }
    }

    return C;
}

/* See if we can reuse an existing matrix on an
   auxiliary node. If so, return it; otherwise
   free it and return a newly allocated matrix.
*/

static gretl_matrix *calc_get_matrix (gretl_matrix **pM,
                                      int r, int c)
{
    if (*pM == NULL) {
        /* allocate from scratch */
        return gretl_matrix_alloc(r, c);
    } else if ((*pM)->rows == r && (*pM)->cols == c) {
        /* reusable as-is */
        return *pM;
    } else if ((*pM)->rows == c && (*pM)->cols == r) {
        /* reusable if reoriented */
        (*pM)->rows = r;
        (*pM)->cols = c;
        return *pM;
    } else {
        /* new matrix needed */
        gretl_matrix_free(*pM);
        *pM = NULL;
        return gretl_matrix_alloc(r, c);
    }
}

#define op_no_complex(o) ((o >= B_LT && o <= B_GTE) || \
                          (o >= B_DOTLT && o <= B_DOTGTE))

#define fn_no_complex(f) (f == F_QFORM || f == F_LSOLVE || \
                          f == F_CMULT || f == F_CDIV ||   \
                          f == F_CONV2D || f == F_SGN ||   \
                          f == F_MAX || f == F_MIN)

/* return allocated result of binary operation performed on
   two matrices */

static int real_matrix_calc (const gretl_matrix *A,
                             const gretl_matrix *B,
                             int op, gretl_matrix **pM)
{
    GretlMatrixMod mod;
    gretl_matrix *C = NULL;
    int ra, ca;
    int rb, cb;
    int r, c;
    int err = 0;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(B)) {
        if (op != B_HCAT && op != B_VCAT && op != F_DSUM &&
            op != B_MUL && op != B_TRMUL) {
            return E_NONCONF;
        }
        if (op == B_MUL || op == B_TRMUL) {
            C = nullmat_multiply(A, B, op, &err);
            goto finish;
        }
    }

    if (A->is_complex || B->is_complex) {
        /* gatekeeper for complex */
        if (op_no_complex(op)) {
            return operator_real_only(op);
        } else if (fn_no_complex(op)) {
            return function_real_only(op);
        }
    }

    switch (op) {
    case B_ADD:
    case B_SUB:
        if (A->is_complex || B->is_complex) {
            int sgn = (op == B_SUB)? -1 : 1;

            C = gretl_cmatrix_add_sub(A, B, sgn, &err);
        } else if (gretl_matrix_is_scalar(A) ||
                   gretl_matrix_is_scalar(B)) {
            C = matrix_add_sub_scalar(A, B, op);
            if (C == NULL) {
                err = E_ALLOC;
            }
        } else {
            C = calc_get_matrix(pM, A->rows, A->cols);
            if (C == NULL) {
                err = E_ALLOC;
            } else if (op == B_ADD) {
                err = gretl_matrix_add(A, B, C);
            } else {
                err = gretl_matrix_subtract(A, B, C);
            }
        }
        break;
    case B_HCAT:
    case B_VCAT:
        if (op == B_HCAT) {
            C = gretl_matrix_col_concat(A, B, &err);
        } else {
            C = gretl_matrix_row_concat(A, B, &err);
        }
        break;
    case F_DSUM:
        C = gretl_matrix_direct_sum(A, B, &err);
        break;
    case B_MUL:
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_multiply(A, B, &err);
        } else {
            ra = gretl_matrix_rows(A);
            ca = gretl_matrix_cols(A);
            rb = gretl_matrix_rows(B);
            cb = gretl_matrix_cols(B);
            r = (ra == 1 && ca == 1)? rb : ra;
            c = (rb == 1 && cb == 1)? ca : cb;

            C = calc_get_matrix(pM, r, c);
            if (C == NULL) {
                err = E_ALLOC;
            } else {
                err = gretl_matrix_multiply(A, B, C);
                if (!err) {
                    gretl_matrix_transcribe_obs_info(C, A);
                }
            }
        }
        break;
    case B_TRMUL:
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_AHB(A, B, &err);
        } else {
            ra = gretl_matrix_cols(A);
            ca = gretl_matrix_rows(A);
            rb = gretl_matrix_rows(B);
            cb = gretl_matrix_cols(B);

            r = (ra == 1 && ca == 1)? rb : ra;
            c = (rb == 1 && cb == 1)? ca : cb;

            C = calc_get_matrix(pM, r, c);
            if (C == NULL) {
                err = E_ALLOC;
            } else {
                err = gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
                                                B, GRETL_MOD_NONE,
                                                C, GRETL_MOD_NONE);
            }
        }
        break;
    case F_QFORM:
        /* quadratic form, A * B * A', for symmetric B */
        ra = gretl_matrix_rows(A);
        ca = gretl_matrix_cols(A);
        rb = gretl_matrix_rows(B);
        cb = gretl_matrix_cols(B);

        /* input checks */
        if (rb == 1 || cb == 1) {
            /* the diagonal B case */
            if (ca != gretl_vector_get_length(B)) {
                err = E_NONCONF;
            }
        } else {
            /* the "full B" case */
            if (ca != rb || cb != rb) {
                err = E_NONCONF;
            } else {
                gretl_matrix_set_equals_tolerance(1.0e-7);
                if (!gretl_matrix_is_symmetric(B)) {
                    gretl_errmsg_set(_("Matrix is not symmetric"));
                    err = E_NONCONF;
                }
                gretl_matrix_unset_equals_tolerance();
            }
        }
        if (!err) {
            mod = GRETL_MOD_NONE;
            C = calc_get_matrix(pM, ra, ra);
            if (C == NULL) {
                err = E_ALLOC;
            } else if (rb == 1 || cb == 1) {
                err = gretl_matrix_diag_qform(A, mod, B, C, mod);
            } else {
                err = gretl_matrix_qform(A, mod, B, C, mod);
            }
        }
        break;
    case B_LDIV:
    case B_DIV:
        /* Matrix left (A\B) or right (A/B) "division". Note that
           A/B = (B'\A')', which we handle by passing the transpose
           flag to gretl_{c}matrix_divide.
        */
        mod = (op == B_LDIV)? GRETL_MOD_NONE : GRETL_MOD_TRANSPOSE;
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_divide(A, B, mod, &err);
        } else {
            C = gretl_matrix_divide(A, B, mod, &err);
        }
        break;
    case F_LSOLVE:
        C = calc_get_matrix(pM, B->rows, B->cols);
        if (C == NULL) {
            err = E_ALLOC;
        } else {
            gretl_matrix_copy_values(C, B);
            err = gretl_cholesky_solve(A, C);
        }
        break;
    case B_DOTMULT:
    case B_DOTDIV:
    case B_DOTPOW:
    case B_DOTADD:
    case B_DOTSUB:
    case B_DOTEQ:
    case B_DOTGT:
    case B_DOTLT:
    case B_DOTGTE:
    case B_DOTLTE:
    case B_DOTNEQ:
        /* apply operator element-wise */
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_dot_op(A, B, op_symbol(op), &err);
        } else {
            C = gretl_matrix_dot_op(A, B, op_symbol(op), &err);
        }
        break;
    case B_KRON:
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_kronecker(A, B, &err);
        } else {
            C = gretl_matrix_kronecker_product_new(A, B, &err);
        }
        break;
    case F_HDPROD:
        if (A->is_complex || B->is_complex) {
            C = gretl_cmatrix_hdprod(A, B, &err);
        } else {
            C = gretl_matrix_hdproduct_new(A, B, &err);
        }
        break;
    case F_CMULT:
        C = gretl_matrix_complex_multiply(A, B, 0, &err);
        break;
    case F_CDIV:
        C = gretl_matrix_complex_divide(A, B, 0, &err);
        break;
    case F_MRSEL:
        C = gretl_matrix_bool_sel(A, B, 1, &err);
        break;
    case F_MCSEL:
        C = gretl_matrix_bool_sel(A, B, 0, &err);
        break;
    case F_CONV2D:
        C = gretl_matrix_2d_convolution(A, B, &err);
        break;
    default:
        err = E_TYPES;
        break;
    }

    if (err) {
        if (C != NULL) {
            if (pM != NULL && *pM == C) {
                *pM = NULL;
            }
            gretl_matrix_free(C);
            C = NULL;
        }
    } else {
        /* preserve data-row info? */
        int At1 = gretl_matrix_get_t1(A);
        int At2 = gretl_matrix_get_t2(A);
        int Bt1 = gretl_matrix_get_t1(B);
        int Bt2 = gretl_matrix_get_t2(B);

        if (C->rows == A->rows && At1 >= 0 && At2 > At1) {
            gretl_matrix_set_t1(C, At1);
            gretl_matrix_set_t2(C, At2);
        } else if (C->rows == B->rows && Bt1 >= 0 && Bt2 > Bt1) {
            gretl_matrix_set_t1(C, Bt1);
            gretl_matrix_set_t2(C, Bt2);
        }
    }

 finish:

    if (*pM != NULL && *pM != C) {
        /* we neither freed nor reused *pM */
        gretl_matrix_free(*pM);
    }

    *pM = C;

    return err;
}

static gretl_matrix *tmp_matrix_from_series (NODE *n, parser *p)
{
    int T = sample_size(p->dset);
    const double *x = n->v.xvec;
    gretl_matrix *m = NULL;

    m = gretl_column_vector_alloc(T);

    if (m == NULL) {
        p->err = E_ALLOC;
    } else {
        memcpy(m->val, x + p->dset->t1, T * sizeof *x);
    }

    return m;
}

/* "Fake" a series using a column vector: the vector must be
   of the same length as the current dataset.
*/

const double *get_colvec_as_series (NODE *n, int f, parser *p)
{
    double *ret = NULL;

    if (n->t != MAT) {
        node_type_error(f, 1, SERIES, n, p);
    } else {
        const gretl_matrix *m = n->v.m;

        if (m->rows == p->dset->n && m->cols == 1) {
            ret = m->val;
        } else {
            node_type_error(f, 1, SERIES, n, p);
        }
    }

    return ret;
}

/* As get_colvec_as_series(), above, except that we'll tolerate
   a column vector shorter than the length of the current dataset,
   provided it equals the length of the current sample range.
*/

const double *get_colvec_as_series_with_size (NODE *n,
                                              int *size,
                                              parser *p)
{
    double *ret = NULL;

    if (n->t != MAT) {
        p->err = e_types(n);
    } else {
        const gretl_matrix *m = n->v.m;

        if (m->cols != 1) {
            p->err = E_INVARG;
        } else if (m->rows == p->dset->n) {
            *size = p->dset->n;
            ret = m->val;
        } else if (m->rows == sample_size(p->dset)) {
            *size = sample_size(p->dset);
            ret = m->val;
        } else {
            p->err = E_TYPES;
        }
    }

    return ret;
}

/* One of the operands is a matrix (or scalar), the other
   a series: we "cast" the series to a matrix.
*/

static NODE *matrix_series_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *a, *b, *tmp;

        if (l->t == SERIES) {
            tmp = a = tmp_matrix_from_series(l, p);
            b = node_get_real_matrix(r, p, 0, 0);
        } else {
            a = node_get_real_matrix(l, p, 0, 0);
            tmp = b = tmp_matrix_from_series(r, p);
        }

        if (!p->err) {
            p->err = real_matrix_calc(a, b, op, &ret->v.m);
        }

        gretl_matrix_free(tmp);
    }

    return ret;
}

static NODE *array_str_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretl_array *a = l->v.a;

        if (gretl_array_get_type(a) != GRETL_TYPE_STRINGS) {
            p->err = e_types(l);
        } else {
            int i, n = gretl_array_get_length(a);
            const char *si;

            ret->v.m = gretl_zero_matrix_new(1, n);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
            } else {
                for (i=0; i<n; i++) {
                    si = gretl_array_get_data(a, i);
                    if (op == B_DOTEQ) {
                        if (si != NULL && !strcmp(si, r->v.str)) {
                            ret->v.m->val[i] = 1;
                        }
                    } else if (si == NULL || strcmp(si, r->v.str)) {
                        ret->v.m->val[i] = 1;
                    }
                }
            }
        }
    }

    return ret;
}

#define comparison_op(o) (o == B_EQ  || o == B_NEQ || \
                          o == B_LT  || o == B_GT ||  \
                          o == B_LTE || o == B_GTE)

/* Here we know have a scalar and a 1 x 1 matrix to work with,
   in either order */

static NODE *matrix_scalar_calc2 (NODE *l, NODE *r, int op,
                                  parser *p)
{
    NODE *ret;

    if (scalar_node(l) && scalar_node(r) && comparison_op(op)) {
        ret = aux_scalar_node(p);
    } else if (l->t == NUM && (op == B_MOD || op == B_POW)) {
        /* the matrix on the right is functioning as
           a scalar argument, so produce a scalar
        */
        ret = aux_scalar_node(p);
    } else {
        /* one of the operands is a matrix, albeit 1 x 1,
           so it's safer to produce a matrix result
        */
        ret = aux_sized_matrix_node(p, 1, 1, 0);
    }

    if (!p->err) {
        double x, y;

        if (l->t == NUM) {
            x = l->v.xval;
            y = r->v.m->val[0];
        } else {
            x = l->v.m->val[0];
            y = r->v.xval;
        }

        if (ret->t == NUM) {
            ret->v.xval = xy_calc(x, y, op, NUM, p);
        } else {
            ret->v.m->val[0] = xy_calc(x, y, op, MAT, p);
        }
    }

    return ret;
}

/* Mixed types: one of the operands is a matrix, the other a scalar,
   giving a matrix result unless we're looking at a comparison
   operator.
*/

static NODE *matrix_scalar_calc (NODE *l, NODE *r, int op, parser *p)
{
    gretl_matrix *m = NULL;
    int comp = comparison_op(op);
    double x;
    NODE *ret = NULL;

    /* Check for the simple case of scalar and
       1 x 1 matrix, either way round
    */
    if ((l->t == NUM && scalar_node(r)) ||
        (r->t == NUM && scalar_node(l))) {
        return matrix_scalar_calc2(l, r, op, p);
    }

    /* get a scalar @x and matrix @m */
    x = (l->t == NUM)? l->v.xval : r->v.xval;
    m = (l->t == MAT)? l->v.m : r->v.m;

    if (gretl_is_null_matrix(m)) {
        p->err = E_DATA;
        return NULL;
    }

    /* mod, pow: the right-hand term must be scalar */
    if ((op == B_MOD || op == B_POW) && !scalar_node(r)) {
        p->err = e_types(r);
        return NULL;
    }

    if (comp) {
        ret = aux_scalar_node(p);
    } else if (op == B_POW) {
        ret = aux_matrix_node(p);
    } else {
        ret = aux_sized_matrix_node(p, m->rows, m->cols, m->is_complex);
    }

    if (ret == NULL) {
        return NULL;
    }

    if (op == B_POW) {
        /* note: the (scalar, 1x1) and (1x1, scalar) cases are
           handled above
        */
        double s = node_get_scalar(r, p);

        if (!p->err) {
            ret->v.m = gretl_matrix_pow(m, s, &p->err);
        }
        return ret;
    } else {
        int i, n = m->rows * m->cols;

        if (comp && m->is_complex) {
            /* B_EQ and B_NOTEQ; needs special treatment */
            ret->v.xval = cmatrix_xy_comp(m, x, op, p);
        } else if (comp) {
            /* condition assumed true by until shown false */
            double cond = 1;

            for (i=0; i<n && cond==1; i++) {
                if (l->t == NUM) {
                    cond = xy_calc(x, m->val[i], op, MAT, p);
                } else {
                    cond = xy_calc(m->val[i], x, op, MAT, p);
                }
            }
            ret->v.xval = cond;
        } else {
            int xleft = (l->t == NUM);

            if (m->is_complex) {
                p->err = cmatrix_xy_calc(ret->v.m, m, x, xleft, op, p);
            } else {
                p->err = rmatrix_xy_calc(ret->v.m, m, x, xleft, op, p);
            }
            if (gretl_matrix_is_dated(m)) {
                gretl_matrix_set_t1(ret->v.m, gretl_matrix_get_t1(m));
                gretl_matrix_set_t2(ret->v.m, gretl_matrix_get_t2(m));
            }
        }
    }

    return ret;
}

static NODE *matrix_transpose_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        if (is_tmp_node(n)) {
            /* transpose temp matrix in place */
            if (n->v.m->is_complex) {
                p->err = gretl_ctrans_in_place(n->v.m);
            } else {
                p->err = gretl_matrix_transpose_in_place(n->v.m);
            }
            ret = n;
        } else {
            /* create transpose as new matrix */
            ret = aux_matrix_node(p);
            if (!p->err) {
                if (n->v.m->is_complex) {
                    ret->v.m = gretl_ctrans(n->v.m, 1, &p->err);
                } else {
                    ret->v.m = gretl_matrix_copy_transpose(n->v.m);
                }
                if (ret->v.m == NULL) {
                    p->err = E_ALLOC;
                }
            }
        }
    } else {
        ret = is_tmp_node(n) ? n : aux_matrix_node(p);
    }

    return ret;
}

/* We're looking at a string argument that is supposed to represent
   a function call: we'll do a rudimentary heuristic check here.
   FIXME this should be more rigorous?
*/

static int is_function_call (const char *s)
{
    if (!strchr(s, '(') || !strchr(s, ')')) {
        gretl_errmsg_sprintf("%s: expected a function call", s);
        return 0;
    } else {
        return 1;
    }
}

static NODE *numeric_jacobian_or_hessian (NODE *l, NODE *m, NODE *r,
                                          int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *s = m->v.str;
        double eps = 0.0;

        if (!is_function_call(s)) {
            p->err = E_TYPES;
            return NULL;
        }

        ret = aux_matrix_node(p);
        if (ret == NULL) {
            return NULL;
        }

        if (!null_node(r)) {
            eps = node_get_scalar(r, p);
        }

        if (!p->err) {
            if (f == F_FDJAC) {
                ret->v.m = user_fdjac(l->v.m, s, eps, p->dset, &p->err);
            } else {
                ret->v.m = user_numhess(l->v.m, s, eps, p->dset, &p->err);
            }
        }
    } else {
        ret = aux_matrix_node(p);
    }

    return ret;
}

/* try to get a (named) uservar of type @t from a U_ADDR node */

static user_var *ptr_node_get_uvar (NODE *n, int t, parser *p)
{
    user_var *uv = NULL;

    if (n->t == U_ADDR) {
        if (n->L->t != t) {
            p->err = e_types(n->L);
        } else {
            uv = n->L->uv;
            if (uv == NULL && n->L->vname != NULL) {
                uv = n->L->uv = get_user_var_by_name(n->L->vname);
            }
            if (uv == NULL) {
                p->err = E_DATA;
            }
        }
    } else {
        p->err = e_types(n);
    }

    return uv;
}

static gretl_matrix *ptr_node_get_matrix (NODE *n, parser *p)
{
    user_var *uv = ptr_node_get_uvar(n, MAT, p);

    return uv != NULL ? uv->ptr : NULL;
}

static gretl_bundle *ptr_node_get_bundle (NODE *n, parser *p)
{
    user_var *uv = ptr_node_get_uvar(n, BUNDLE, p);

    return uv != NULL ? uv->ptr : NULL;
}

/* To be called ONLY for a node @n whose verified type is U_ADDR */

static NODE *ptr_node_get_referent_node (NODE *n, parser *p)
{
    NODE *ret = n->L;

    if (ret == NULL) {
        p->err = E_DATA;
    } else if (ret->t == EMPTY) {
        return ret;
    } else if (ret->t == SERIES) {
        if (ret->v.xvec == NULL) {
            int id = current_series_index(p->dset, ret->vname);

            if (id > 0) {
                ret->v.xvec = p->dset->Z[id];
            }
        }
        if (ret->v.xvec == NULL) {
            p->err = E_DATA;
            ret = NULL;
        }
    } else {
        if (ret->uv == NULL) {
            ret->uv = get_user_var_by_name(ret->vname);
        }
        if (ret->uv == NULL) {
            p->err = E_DATA;
            ret = NULL;
        }
    }

    return ret;
}

/* note: allows @n to be either a regular matrix node or a
   matrix-pointer node
*/

static gretl_matrix *mat_node_get_real_matrix (NODE *n, parser *p)
{
    gretl_matrix *m = NULL;

    if (n == NULL) {
        p->err = E_DATA;
    } else if (n->t == MAT) {
        m = n->v.m;
    } else if (n->t == U_ADDR) {
        m = ptr_node_get_matrix(n, p);
    } else {
        p->err = e_types(n);
    }

    if (!p->err) {
        if (m == NULL) {
            p->err = E_DATA;
        } else if (m->is_complex) {
            p->err = E_CMPLX;
            m = NULL;
        }
    }

    return m;
}

static const char *node_get_fncall (NODE *n, parser *p)
{
    const char *ret = NULL;

    if (n->t != STR) {
        p->err = e_types(n);
    } else {
        ret = n->v.str;
        if (!is_function_call(ret)) {
            p->err = E_TYPES;
        }
    }

    return ret;
}

static NODE *BFGS_constrained_max (NODE *t, NODE *n, parser *p)
{
    NODE *ret = NULL;
    NODE *e = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *bounds = NULL;
    const char *sf = NULL;
    const char *sg = NULL;
    int i, k = n->v.bn.n_nodes;

    if (k < 3 || k > 4) {
        n_args_error(k, 3, 4, F_BFGSCMAX, p);
    }

    for (i=0; i<k && !p->err; i++) {
        e = n->v.bn.n[i];
        if (i == 0) {
            b = mat_node_get_real_matrix(e, p);
        } else if (i == 1) {
            bounds = mat_node_get_real_matrix(e, p);
        } else if (i == 2) {
            sf = node_get_fncall(e, p);
        } else if (i == 3 && !null_node(e)) {
            sg = node_get_fncall(e, p);
        }
    }

    if (!p->err) {
        ret = aux_scalar_node(p);
    }

    if (!p->err) {
        int minimize = alias_reversed(t) ? 1 : 0;

        ret->v.xval = user_BFGS(b, sf, sg, p->dset, bounds,
                                minimize, p->prn, &p->err);
    }

    return ret;
}

static NODE *BFGS_maximize (NODE *l, NODE *m, NODE *r,
                            parser *p, NODE *t)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *b;
        const char *sf = m->v.str;
        const char *sg = NULL;

        b = mat_node_get_real_matrix(l, p);

        if (!p->err) {
            if (r->t == STR) {
                sg = r->v.str;
            } else if (r->t != EMPTY) {
                p->err = e_types(r);
            }
        }
        if (!p->err && !is_function_call(sf)) {
            p->err = E_TYPES;
        }
        if (!p->err && sg != NULL && !is_function_call(sg)) {
            p->err = E_TYPES;
        }
        if (!p->err && gretl_is_null_matrix(b)) {
            p->err = E_DATA;
        }
        if (p->err) {
            return NULL;
        }

        ret = aux_scalar_node(p);
        if (ret != NULL) {
            int minimize = alias_reversed(t) ? 1 : 0;

            ret->v.xval = user_BFGS(b, sf, sg, p->dset, NULL,
                                    minimize, p->prn, &p->err);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *deriv_free_node (NODE *l, NODE *m, NODE *r,
                              parser *p, NODE *t)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *b = NULL;
        const char *fcall = m->v.str;
        double tol = NADBL;
        int maxit = 0;

        b = mat_node_get_real_matrix(l, p);
        if (!p->err) {
            if (gretl_is_null_matrix(b)) {
                p->err = E_DATA;
            } else if (!is_function_call(fcall)) {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            if (scalar_node(r)) {
                if (t->t == F_GSSMAX) {
                    tol = r->v.xval;
                } else {
                    maxit = node_get_int(r, p);
                }
            } else if (!null_node(r)) {
                p->err = e_types(r);
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (ret != NULL) {
            int minimize = alias_reversed(t) ? 1 : 0;
            MaxMethod method = SIMANN_MAX;

            if (t->t == F_NMMAX) {
                method = NM_MAX;
            } else if (t->t == F_GSSMAX) {
                method = GSS_MAX;
            }
            ret->v.xval = deriv_free_optimize(method, b, fcall, maxit, tol,
                                              minimize, p->dset, p->prn,
                                              &p->err);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *fzero_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *fcall = l->v.str;
        gretl_matrix *b = NULL;
        double tol = NADBL;
        int free_b = 0;

        if (null_or_scalar(m)) {
            b = gretl_matrix_alloc(1, 2);
            b->val[0] = null_node(m) ? NADBL : node_get_scalar(m, p);
            b->val[1] = NADBL;
            free_b = 1;
        } else if (m->t == MAT) {
            b = m->v.m;
        } else {
            p->err = e_types(m);
        }
        if (!p->err) {
            if (gretl_is_null_matrix(b)) {
                p->err = E_DATA;
            } else if (!is_function_call(fcall)) {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            if (scalar_node(r)) {
                tol = node_get_scalar(r, p);
            } else if (!null_node(r)) {
                p->err = e_types(r);
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (ret != NULL) {
            ret->v.xval = deriv_free_optimize(ROOT_FIND, b, fcall, 0,
                                              tol, 0, p->dset, p->prn,
                                              &p->err);
        }
        if (free_b) {
            gretl_matrix_free(b);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

static void lag_calc (double *y, const double *x,
                      int k, int t1, int t2,
                      int op, double mul,
                      parser *p)
{
    int s, t;

    for (t=t1; t<=t2; t++) {
        s = t - k;
        if (dataset_is_panel(p->dset)) {
            if (s / p->dset->pd != t / p->dset->pd) {
                /* s and t pertain to different units */
                s = -1;
            }
        }
        if (s >= 0 && s < p->dset->n) {
            if (op == B_ASN && mul == 1.0) {
                y[t] = x[s];
            } else if (op == B_ASN) {
                y[t] = mul * x[s];
            } else if (op == B_ADD) {
                y[t] += mul * x[s];
            } else {
                p->err = E_DATA;
            }
        }
    }
}

static NODE *matrix_file_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *a = node_get_matrix(l, p, 0, 1);
        const char *fname = m->v.str;
        int export = node_get_bool(r, p, 0);

        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            p->err = gretl_matrix_write_to_file(a, fname, export);
            ret->v.xval = p->err;
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *bundle_file_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *s = m->v.str;
        int control = 0;

        if (!null_node(r)) {
            control = node_get_int(r, p);
        }
        ret = aux_scalar_node(p);
        if (ret != NULL) {
            p->err = ret->v.xval = gretl_bundle_write_to_file(l->v.b, s, control);
        }
    } else {
        ret = aux_scalar_node(p);
    }

    return ret;
}

static int check_cswitch_param (gretl_matrix *m, int *k)
{
    int err = 0;

    if (*k < 0 || *k > 4) {
        /* out of bounds */
        err = E_INVARG;
    } else if (*k == 1) {
        /* real to complex, column-wise */
        if (m->cols % 2) {
            err = E_NONCONF;
        }
    } else if (*k == 3) {
        /* real to complex, row-wise */
        if (m->rows % 2) {
            err = E_NONCONF;
        }
    }

    return err;
}

/* matrix on left, scalar(s) on right: returns a matrix */

static NODE *matrix_scalar_func (NODE *l, NODE *r,
                                 int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *m = l->v.m;
        int k = 0;

        if (f == F_CSWITCH) {
            k = 1; /* default */
            if (!null_node(r)) {
                k = node_get_int(r, p);
            }
            if (!p->err) {
                p->err = check_cswitch_param(m, &k);
            }
        } else {
            k = node_get_int(r, p);
        }
        if (!p->err && gretl_is_null_matrix(m)) {
            p->err = E_INVARG;
        }
        if (!p->err && f == F_MSORTBY && m->is_complex) {
            p->err = E_CMPLX;
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (p->err) {
            return NULL;
        }

        if (f == F_MSORTBY) {
            ret->v.m = gretl_matrix_sort_by_column(m, k-1, &p->err);
        } else if (f == F_CSWITCH) {
            if (k > 2) {
                /* the old _setcmplx() */
                k = (k == 3)? 1 : 0;
                ret->v.m = gretl_matrix_copy(m);
                if (ret->v.m == NULL) {
                    p->err = E_ALLOC;
                } else {
                    p->err = gretl_matrix_set_complex_full(ret->v.m, k);
                }
            } else {
                k = (k == 1)? 1 : 0;
                ret->v.m = gretl_cmatrix_switch(m, k, &p->err);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *matrix_vector_func (NODE *l, NODE *m, NODE *r,
                                 int f, parser *p)
{
    NODE *ret = NULL;

    /* at present only F_MSPLITBY comes here */

    if (starting(p)) {
        gretl_matrix *a = node_get_matrix(l, p, 0, 1);
        gretl_matrix *v = node_get_matrix(m, p, 1, 2);
        int colwise = 0;

        if (!p->err) {
            colwise = node_get_bool(r, p, 0);
        }
        if (!p->err) {
            ret = aux_array_node(p);
        }
        if (ret != NULL) {
            int chunks = m->t == NUM; /* scalar */

            ret->v.a = gretl_matrix_split_by(a, v, colwise, chunks, &p->err);
        }
    } else {
        ret = aux_array_node(p);
    }

    return ret;
}

static NODE *matrix_string_func (NODE *l, NODE *m, NODE *r,
                                 int f, parser *p)
{
    NODE *ret = NULL;

    /* at present only F_DISTANCE comes here */

    if (starting(p)) {
        gretl_matrix *X = l->v.m;
        gretl_matrix *Y = NULL;
        const char *s = "e"; /* euclidean default */

        if (!null_node(m)) {
            if (m->t == STR) {
                s = m->v.str;
            } else {
                node_type_error(f, 2, STR, m, p);
            }
        }
        if (!null_node(r)) {
            if (r->t == MAT) {
                Y = r->v.m;
            } else {
                node_type_error(f, 3, MAT, r, p);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (ret != NULL) {
            ret->v.m = distance(X, Y, s, &p->err);
        }
    } else {
        ret = aux_matrix_node(p);
    }

    return ret;
}

/* both operands are known to be matrices or scalars */

static NODE *matrix_matrix_calc (NODE *l, NODE *r, int op, parser *p)
{
    gretl_matrix *ml = NULL, *mr = NULL;
    NODE *ret = NULL;

    if ((op == B_MUL || op == B_TRMUL || op == B_ADD || op == B_SUB) &&
        l->t == MAT && r->t == MAT) {
        ml = l->v.m;
        mr = r->v.m;
        if (ml->is_complex || mr->is_complex) {
            ret = aux_matrix_node(p);
            if (!p->err) {
                p->err = real_matrix_calc(ml, mr, op, &ret->v.m);
            }
            return ret;
        }
    }

    if (op == B_DOTPOW || op == B_POW) {
        if (op == B_POW) {
            if (scalar_node(l) && scalar_node(r)) {
                op = B_DOTPOW;
            } else if (!scalar_node(r)) {
                p->err = e_types(r);
                return NULL;
            }
        }
        ret = aux_matrix_node(p);
    } else {
        /* experiment: try reusing aux matrix */
        p->flags |= P_MSAVE;
        ret = get_aux_node(p, MAT, 0, TMP_NODE);
        p->flags ^= P_MSAVE;
    }

#if EDEBUG
    fprintf(stderr, "matrix_matrix_calc: l=%p, r=%p, ret=%p\n",
            (void *) l, (void *) r, (void *) ret);
#endif

    if (ml == NULL) {
        ml = node_get_matrix(l, p, 0, 1);
        if (op != B_POW) {
            mr = node_get_matrix(r, p, 1, 2);
        }
    }

    if (ret != NULL && starting(p)) {
        if (op == B_DOTPOW) {
            if (ml->is_complex) {
                ret->v.m = gretl_cmatrix_dot_op(ml, mr, '^', &p->err);
            } else {
                ret->v.m = gretl_matrix_dot_op(ml, mr, '^', &p->err);
            }
        } else if (op == B_POW) {
            int s = node_get_int(r, p);

            if (!p->err) {
                ret->v.m = gretl_matrix_pow(ml, s, &p->err);
            }
        } else {
            p->err = real_matrix_calc(ml, mr, op, &ret->v.m);
        }
    }

    return ret;
}

static NODE *matrix_and_or (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const gretl_matrix *a = l->v.m;
        const gretl_matrix *b = r->v.m;
        int i, n = a->rows * a->cols;

        if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
            p->err = E_NONCONF;
        } else if (a->rows != b->rows || a->cols != b->cols) {
            p->err = E_NONCONF;
        } else {
            ret->v.m = gretl_unit_matrix_new(a->rows, a->cols);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
                return NULL;
            }
            for (i=0; i<n; i++) {
                if (op == B_AND) {
                    if (a->val[i] == 0.0 || b->val[i] == 0.0) {
                        ret->v.m->val[i] = 0.0;
                    }
                } else if (op == B_OR) {
                    if (a->val[i] == 0.0 && b->val[i] == 0.0) {
                        ret->v.m->val[i] = 0.0;
                    }
                }
            }
        }
    }

    return ret;
}

/* both operands are matrices */

static NODE *matrix_bool (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret;

    if (op == B_OR || op == B_AND) {
        return matrix_and_or(l, r, op, p);
    }

    ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const gretl_matrix *a = l->v.m;
        const gretl_matrix *b = r->v.m;
        int i, n = a->rows * a->cols;

        if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
            ret->v.xval = NADBL;
        } else if (a->rows != b->rows || a->cols != b->cols) {
            ret->v.xval = NADBL;
        } else {
            ret->v.xval = op == B_NEQ ? 0 : 1;
            for (i=0; i<n; i++) {
                if (op == B_EQ && a->val[i] != b->val[i]) {
                    ret->v.xval = 0;
                    break;
                } else if (op == B_LT && a->val[i] >= b->val[i]) {
                    ret->v.xval = 0;
                    break;
                } else if (op == B_GT && a->val[i] <= b->val[i]) {
                    ret->v.xval = 0;
                    break;
                } else if (op == B_LTE && a->val[i] > b->val[i]) {
                    ret->v.xval = 0;
                    break;
                } else if (op == B_GTE && a->val[i] < b->val[i]) {
                    ret->v.xval = 0;
                    break;
                } else if (op == B_NEQ && a->val[i] != b->val[i]) {
                    ret->v.xval = 1;
                    break;
                }
            }
        }
    }

    return ret;
}

static void matrix_error (parser *p)
{
    if (p->err == 0) {
        p->err = 1;
    }

    if (gretl_errmsg_is_set()) {
        errmsg(p->err, p->prn);
    }
}

/* functions taking a matrix argument and returning a
   scalar result */

static NODE *matrix_to_scalar_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *m = node_get_matrix(l, p, 0, 0);

        if (m->is_complex && f != F_ROWS && f != F_COLS && f != F_RANK) {
            /* gatekeeper for complex */
            p->err = function_real_only(f);
            return ret;
        }

        switch (f) {
        case F_ROWS:
            ret->v.xval = m->rows;
            break;
        case F_COLS:
            ret->v.xval = m->cols;
            break;
        case F_NORM1:
            ret->v.xval = gretl_matrix_one_norm(m);
            break;
        case F_INFNORM:
            ret->v.xval = gretl_matrix_infinity_norm(m);
            break;
        case F_RCOND:
            ret->v.xval = gretl_matrix_rcond(m, &p->err);
            break;
        case F_CNUMBER:
            ret->v.xval = gretl_matrix_cond_index(m, &p->err);
            break;
        case F_RANK:
            if (m->is_complex) {
                ret->v.xval = gretl_cmatrix_rank(m, &p->err);
            } else {
                double eps = NADBL;
                if (!null_node(r)) {
                    eps = node_get_scalar(r, p);
                }
                if (!p->err) {
                    ret->v.xval = gretl_matrix_rank(m, eps, &p->err);
                }
            }
            break;
        default:
            p->err = E_PARSE;
            break;
        }

        if (p->err) {
            matrix_error(p);
        }
    }

    return ret;
}

/* Compute a value which will be a scalar for a real matrix
   but a complex scalar (2-vector) for a complex matrix:
   handles determinant and trace.
*/

static NODE *matrix_to_alt_node (NODE *n, int f, parser *p)
{
    gretl_matrix *m = node_get_matrix(n, p, 0, 0);
    NODE *ret = NULL;

    if (!p->err) {
        ret = m->is_complex ? aux_matrix_node(p) : aux_scalar_node(p);
    }

    if (!p->err) {
        if (m->is_complex) {
            if (f == F_TRACE) {
                ret->v.m = gretl_cmatrix_trace(m, &p->err);
            } else {
                ret->v.m = gretl_cmatrix_determinant(m, f==F_LDET, &p->err);
            }
        } else if (f == F_TRACE) {
            ret->v.xval = gretl_matrix_trace(m);
        } else {
            int tmpmat = (n->t == MAT && is_tmp_node(n));

            ret->v.xval =
                user_matrix_get_determinant(m, tmpmat, f==F_LDET, &p->err);
        }
    }

    return ret;
}

static NODE *matrix_add_names (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *m = l->v.m;
        int byrow = (f == F_RNAMESET);

        if (m->is_complex) {
            /* we could set column names for a complex matrix
               but they wouldn't show up on printing
            */
            p->err = E_CMPLX;
            return ret;
        }

        if (r->t == STR) {
            ret->v.xval = umatrix_set_names_from_string(m, r->v.str, byrow);
        } else if (r->t == ARRAY) {
            if (gretl_array_get_type(r->v.a) != GRETL_TYPE_STRINGS) {
                p->err = e_types(r);
            } else {
                p->err = ret->v.xval = umatrix_set_names_from_array(m, r->v.a,
                                                                    byrow);
            }
        } else {
            /* some sort of list-bearing node */
            int *list = node_get_list(r, p);

            if (p->err) {
                ret->v.xval = 1;
            } else {
                p->err = ret->v.xval = umatrix_set_names_from_list(m, list, p->dset,
                                                                   byrow);
            }
            free(list);
        }
    }

    return ret;
}

static NODE *matrix_get_col_or_row_name (int f, NODE *l, NODE *r,
                                         parser *p)
{
    int get_all = null_node(r);
    NODE *ret = get_all ? aux_array_node(p) : aux_string_node(p);

    if (ret != NULL && starting(p)) {
        if (get_all) {
            const char **S;
            int n = 0;

            if (f == F_CNAMEGET) {
                S = gretl_matrix_get_colnames(l->v.m);
                if (S != NULL) n = l->v.m->cols;
            } else {
                S = gretl_matrix_get_rownames(l->v.m);
                if (S != NULL) n = l->v.m->rows;
            }
            ret->v.a = gretl_array_from_strings((char **) S, n,
                                                1, &p->err);
        } else {
            int i = node_get_int(r, p);

            if (f == F_CNAMEGET) {
                ret->v.str = user_matrix_get_column_name(l->v.m, i, &p->err);
            } else {
                ret->v.str = user_matrix_get_row_name(l->v.m, i, &p->err);
            }
        }
    }

    return ret;
}

static NODE *matrix_imhof (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const gretl_matrix *m = l->v.m;
        double arg = node_get_scalar(r, p);

        ret->v.xval = imhof(m, arg, &p->err);
    }

    return ret;
}

static NODE *bkw_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const gretl_matrix *V = l->v.m;
        gretl_array *pnames = NULL;
        PRN *vprn = NULL;
        int ns = 0;

        if (!null_node(m)) {
            if (m->t == STR) {
                /* for compat with Lee's bkw() we expect comma-
                   separated parameter names here
                */
                char **S = gretl_string_split(m->v.str, &ns, ",");

                if (S == NULL) {
                    p->err = E_DATA;
                } else {
                    pnames = gretl_array_from_strings(S, ns, 0, &p->err);
                }
            } else if (m->t == ARRAY) {
                if (gretl_array_get_type(m->v.a) != GRETL_TYPE_STRINGS) {
                    p->err = e_types(m);
                } else {
                    pnames = gretl_array_copy(m->v.a, &p->err);
                }
            } else {
                p->err = e_types(m);
            }
        }

        if (node_get_bool(r, p, 0)) {
            /* optional verbose flag */
            vprn = p->prn;
        }

        if (!p->err) {
            gretl_matrix *(*bkwfunc) (const gretl_matrix *, gretl_array *,
                                      PRN *, int *);

            bkwfunc = get_plugin_function("bkw_matrix");
            if (bkwfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                ret->v.m = (*bkwfunc)(V, pnames, vprn, &p->err);
            }
        }

        gretl_array_destroy(pnames);
    }

    return ret;
}

/* Here we handle the case where the relevant libgretl C function
   overwrites its matrix argument. If @m is just an on-the-fly matrix
   it can be passed on as argument to the function in question, but if
   it's a named user-matrix we'll have to make a copy to pass.
*/

static gretl_matrix *do_matrix_ovwrite (gretl_matrix *m,
                                        int f, int parm,
                                        double xparm,
                                        int tmpmat,
                                        int *err)
{
    gretl_matrix *R = m;

    if (f == F_CHOL && !gretl_is_null_matrix(m) &&
        !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_set(_("Matrix is not symmetric"));
        *err = E_DATA;
        return NULL;
    }

    if (!tmpmat) {
        /* shouldn't overwrite @m */
        R = gretl_matrix_copy(m);
        if (R == NULL) {
            *err = E_ALLOC;
        }
    }

    if (R != NULL) {
        if (f == F_CHOL) {
            *err = gretl_matrix_cholesky_decomp(R);
        } else if (f == F_PSDROOT) {
            *err = gretl_matrix_psd_root(R, parm);
        } else if (f == F_GINV) {
            *err = gretl_matrix_moore_penrose(R, xparm);
        } else if (f == F_INV) {
            *err = gretl_invert_matrix(R);
        } else if (f == F_UPPER) {
            *err = gretl_matrix_zero_lower(R);
        } else if (f == F_LOWER) {
            *err = gretl_matrix_zero_upper(R);
        } else {
            *err = E_DATA;
        }
        if (*err && R != m) {
            gretl_matrix_free(R);
            R = NULL;
        }
    }

    return R;
}

static void matrix_minmax_indices (int f, int *mm, int *rc, int *idx)
{
    /* mm: 1 for maximum functions, 0 for minimum */
    *mm = (f == F_MAXR || f == F_MAXC || f == F_IMAXR || f == F_IMAXC);
    /* rc: 1 for column-wise functions, 0 for row-wise */
    *rc = (f == F_MINC || f == F_MAXC || f == F_IMINC || f == F_IMAXC);
    /* idx: 1 for index functions, 0 for values */
    *idx = (f == F_IMINR || f == F_IMINC || f == F_IMAXR || f == F_IMAXC);
}

#define mmf_does_complex(f) (f==F_INV || f==F_UPPER || f==F_LOWER || \
                             f==F_DIAG || f==F_TRANSP || f==F_CTRANS || \
                             f==F_VEC || f==F_VECH || f==F_UNVECH || \
                             f==F_MREV || f== F_FFT || f==F_FFTI || \
                             f==F_CUM || f==F_DIFF || f==F_SUMC || \
                             f==F_SUMR || f==F_PRODC || f==F_PRODR || \
                             f==F_MEANC || f==F_MEANR || f==F_GINV || \
                             f==F_MLOG || f==F_MEXP || f==F_CHOL)

static NODE *matrix_to_matrix_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *m = NULL;
        double xparm = NADBL;
        int tmpmat = 0;
        int parm = 0;
        int a = 0, b = 0, c = 0;

        /* Note: @parm is an integer parameter, required
           for some functions, optional for others. In the
           case of mcov() the default @parm value is 1,
           otherwise it's 0.
        */
        if (f == F_MCOV) {
            parm = 1;
        }

        m = node_get_matrix(n, p, 0, 0);
        tmpmat = n->t == MAT && is_tmp_node(n);

        if (!p->err && m != NULL && m->is_complex) {
            /* gatekeeper for complex */
            if (!mmf_does_complex(f)) {
                p->err = function_real_only(f);
            }
        }

        if (p->err) {
            goto finalize;
        }

        if (f == F_MREV || f == F_MCOV ||
            f == F_PSDROOT || f == F_VECH) {
            /* if present, the @r node should hold an integer */
            if (!null_or_scalar(r)) {
                node_type_error(f, 2, NUM, r, p);
            } else if (!null_node(r)) {
                parm = node_get_int(r, p);
            }
        } else if (f == F_UNVECH || f == F_GINV) {
            /* if present, the @r node should hold a scalar */
            if (!null_or_scalar(r)) {
                node_type_error(f, 2, NUM, r, p);
            } else if (!null_node(r)) {
                xparm = node_get_scalar(r, p);
            }
        } else if (f >= F_MINC && f <= F_MEANR) {
            /* set the skip_na boolean option */
            parm = node_get_bool(r, p, 0);
        } else if (f == F_RANKING) {
            if (gretl_vector_get_length(m) == 0) {
                /* m must be a vector */
                p->err = E_TYPES;
            }
        }

        if (!p->err && gretl_is_null_matrix(m) && !emptymat_ok(f)) {
            p->err = E_DATA;
        }

        if (p->err) {
            goto finalize;
        }

        gretl_error_clear();

        if (gretl_is_null_matrix(m)) {
            ret->v.m = gretl_null_matrix_new();
            goto finalize;
        }

        switch (f) {
        case F_SUMC:
            ret->v.m = gretl_matrix_vector_stat(m, V_SUM, 0, parm, &p->err);
            break;
        case F_SUMR:
            ret->v.m = gretl_matrix_vector_stat(m, V_SUM, 1, parm, &p->err);
            break;
        case F_PRODC:
            ret->v.m = gretl_matrix_vector_stat(m, V_PROD, 0, parm, &p->err);
            break;
        case F_PRODR:
            ret->v.m = gretl_matrix_vector_stat(m, V_PROD, 1, parm, &p->err);
            break;
        case F_MEANC:
            ret->v.m = gretl_matrix_vector_stat(m, V_MEAN, 0, parm, &p->err);
            break;
        case F_MEANR:
            ret->v.m = gretl_matrix_vector_stat(m, V_MEAN, 1, parm, &p->err);
            break;
        case F_MCOV:
            ret->v.m = gretl_covariance_matrix(m, 0, parm, &p->err);
            break;
        case F_MCORR:
            ret->v.m = gretl_covariance_matrix(m, 1, 1, &p->err);
            break;
        case F_CUM:
            ret->v.m = gretl_matrix_cumcol(m, &p->err);
            break;
        case F_DIFF:
            ret->v.m = gretl_matrix_diffcol(m, 0, &p->err);
            break;
        case F_DATAOK:
            ret->v.m = gretl_matrix_isfinite(m, &p->err);
            break;
        case F_INV:
            if (m->is_complex) {
                ret->v.m = gretl_cmatrix_inverse(m, &p->err);
            } else {
                ret->v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p->err);
            }
            break;
        case F_GINV:
            if (m->is_complex) {
                ret->v.m = gretl_cmatrix_ginv(m, &p->err);
            } else {
                ret->v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p->err);
            }
            break;
        case F_CHOL:
            if (m->is_complex) {
                ret->v.m = gretl_cmatrix_cholesky(m, &p->err);
            } else {
                ret->v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p->err);
            }
            break;
        case F_PSDROOT:
        case F_UPPER:
        case F_LOWER:
            ret->v.m = do_matrix_ovwrite(m, f, parm, xparm, tmpmat, &p->err);
            break;
        case F_DIAG:
            ret->v.m = gretl_matrix_get_diagonal(m, &p->err);
            break;
        case F_TRANSP:
            if (m->is_complex) {
                ret->v.m = gretl_ctrans(m, 0, &p->err);
            } else {
                ret->v.m = gretl_matrix_copy_transpose(m);
            }
            break;
        case F_VEC:
            ret->v.m = user_matrix_vec(m, &p->err);
            break;
        case F_VECH:
            ret->v.m = user_matrix_vech(m, parm, &p->err);
            break;
        case F_UNVECH:
            ret->v.m = user_matrix_unvech(m, xparm, &p->err);
            break;
        case F_MREV:
            if (parm != 0) {
                ret->v.m = gretl_matrix_reverse_cols(m, &p->err);
            } else {
                ret->v.m = gretl_matrix_reverse_rows(m, &p->err);
            }
            break;
        case F_NULLSPC:
            ret->v.m = gretl_matrix_right_nullspace(m, &p->err);
            break;
        case F_MEXP:
            if (m->is_complex) {
                ret->v.m = gretl_cmatrix_exp(m, &p->err);
            } else {
                ret->v.m = gretl_matrix_exp(m, &p->err);
            }
            break;
        case F_MLOG:
            ret->v.m = gretl_matrix_log(m, &p->err);
            break;
        case F_FFT:
            if (m->is_complex) {
                ret->v.m = gretl_cmatrix_fft(m, 0, &p->err);
            } else {
                ret->v.m = gretl_matrix_fft(m, &p->err);
            }
            break;
        case F_FFTI:
            ret->v.m = gretl_matrix_ffti(m, &p->err);
            break;
        case F_POLROOTS:
            ret->v.m = gretl_matrix_polroots(m, 1, 0, &p->err);
            break;
        case F_RANKING:
            ret->v.m = rank_vector(m, F_SORT, &p->err);
            break;
        case F_MINC:
        case F_MAXC:
        case F_MINR:
        case F_MAXR:
        case F_IMINC:
        case F_IMAXC:
        case F_IMINR:
        case F_IMAXR:
            matrix_minmax_indices(f, &a, &b, &c);
            ret->v.m = gretl_matrix_minmax(m, a, b, c, parm, &p->err);
            break;
        case F_CTRANS:
            ret->v.m = gretl_ctrans(m, 1, &p->err);
            break;
        default:
            break;
        }

        if (ret->v.m == m && n->t == MAT) {
            /* input matrix was recycled: avoid double-freeing */
            n->v.m = NULL;
        }

    finalize:

        if (ret->v.m == NULL) {
            matrix_error(p);
        }
    }

    return ret;
}

static NODE *glasso_node (NODE *l, NODE *r, parser *p)
{
    gretl_matrix *(*glasso) (const gretl_matrix *,
                             gretl_bundle *,
                             PRN *, int *);
    NODE *ret = NULL;

    glasso = get_plugin_function("gretl_glasso");
    if (glasso == NULL) {
        p->err = E_FOPEN;
    } else {
        ret = aux_matrix_node(p);
    }
    if (ret != NULL) {
        ret->v.m = (*glasso)(l->v.m, r->v.b, p->prn, &p->err);
    }

    return ret;
}

static NODE *list_reverse_node (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int i, nt = n->v.ivec[0];
        int *rev = gretl_list_new(nt);

        for (i=1; i<=nt; i++) {
            rev[i] = n->v.ivec[nt-i+1];
        }
        ret->v.ivec = rev;
    }

    return ret;
}

/* We come here if we got a ".csv" suffix for the argument
   to mread(). If we can open the file as-is, we check to
   make sure it's not a gretl-format matrix file with the
   wrong suffix.
*/

static int check_matrix_file (const char *fname, int *csv)
{
    char line[1024];
    FILE *fp;
    int r, c, n;

    *csv = 1;

    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
        /* just assume it's really CSV */
        return 0;
    }

    while (fgets(line, sizeof line, fp)) {
        if (*line != '#') {
            /* heuristic: if the non-comment portion of the file
               starts with two tab-separated integers, it's
               actually a native gretl .mat file regardless of
               the filename suffix?
            */
            n = sscanf(line, "%d\t%d", &r, &c);
            if (n == 2 && count_fields(line, "\t") == 2) {
                *csv = 0;
            }
            break;
        }
    }

    fclose(fp);

    return 0;
}

static NODE *read_object_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (f == F_MREAD) {
        ret = aux_matrix_node(p);
    } else {
        ret = aux_bundle_node(p);
    }

    if (ret != NULL && starting(p)) {
        const char *fname = n->v.str;
        const char *realpath = fname;
        gchar *tmp = NULL;
        int import = node_get_bool(r, p, 0);
        int csv = 0;
        int gdt = 0;

        gretl_error_clear();

        if (import) {
            tmp = gretl_make_dotpath(fname);
            realpath = tmp;
        }

        if (has_suffix(realpath, ".csv")) {
            p->err = check_matrix_file(realpath, &csv);
            if (p->err) {
                return ret;
            }
        } else if (has_suffix(realpath, ".gdt") ||
                   has_suffix(realpath, ".gdtb")) {
            gdt = 1;
        }

        switch (f) {
        case F_MREAD:
            if (csv) {
                ret->v.m = import_csv_as_matrix(realpath, &p->err);
            } else if (gdt) {
                set_dset_matrix_target(&ret->v.m);
                p->err = gretl_read_gdt(realpath, NULL, OPT_NONE, NULL);
                set_dset_matrix_target(NULL);
            } else {
                ret->v.m = gretl_matrix_read_from_file(realpath, 0, &p->err);
            }
            break;
        case F_BREAD:
            ret->v.b = gretl_bundle_read_from_file(realpath, 0, &p->err);
            break;
        default:
            break;
        }

        g_free(tmp);
    }

    return ret;
}

/* Build a node holding a complex matrix, given two scalars,
   two matrices, or matrix plus scalar.
*/

static NODE *complex_matrix_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL) {
        gretl_matrix *Re = l->t == MAT ? l->v.m : NULL;
        gretl_matrix *Im = r->t == MAT ? r->v.m : NULL;
        double x = l->t == NUM ? l->v.xval : 0;
        double y = r->t == NUM ? r->v.xval : 0;

        if (l->t == NUM && null_or_scalar(r)) {
            ret->v.m = gretl_cmatrix_from_scalar(x + y*I, &p->err);
        } else {
            ret->v.m = gretl_cmatrix_build(Re, Im, x, y, &p->err);
        }
    }

    return ret;
}

/* 2023-02-16: added optional scalar-pointer arg to invpd(),
   to retrieve the log determinant if wanted.
*/

static NODE *invpd_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *m = node_get_matrix(l, p, 0, 0);
        user_var *uv = NULL;

        if (!p->err && gretl_is_null_matrix(m)) {
            p->err = E_DATA;
        }
        if (!p->err && !null_node(r)) {
            uv = ptr_node_get_uvar(r, NUM, p);
        }
        if (!p->err) {
            if (l->t == MAT && is_tmp_node(l)) {
                /* OK to overwrite @m */
                ret->v.m = l->v.m;
                l->v.m = NULL;
            } else {
                /* don't destroy @m */
                ret->v.m = gretl_matrix_copy(m);
                if (ret->v.m == NULL) {
                    p->err = E_ALLOC;
                }
            }
        }
        if (!p->err) {
            if (uv != NULL) {
                double ldet;

                p->err = gretl_invert_symmetric_matrix2(ret->v.m, &ldet);
                user_var_set_scalar_value(uv, ldet);
            } else {
                p->err = gretl_invpd(ret->v.m);
            }
        }
        if (p->err && ret->v.m != NULL) {
            gretl_matrix_free(ret->v.m);
            ret->v.m = NULL;
        }
    }

    return ret;
}

static NODE *
matrix_to_matrix2_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *m1 = node_get_matrix(n, p, 0, 0);
        gretl_matrix *m2 = NULL;

        if (!p->err && gretl_is_null_matrix(m1)) {
            p->err = E_DATA;
        }
        if (!p->err && !null_node(r)) {
            m2 = ptr_node_get_matrix(r, p);
        }

        if (!p->err) {
            if (f == F_EIGSYM) {
                ret->v.m = user_matrix_eigensym(m1, m2, &p->err);
            } else if (f == F_HDPROD) {
                if (m1->is_complex) {
                    ret->v.m = gretl_cmatrix_hdprod(m1, NULL, &p->err);
                } else {
                    ret->v.m = gretl_matrix_hdproduct_new(m1, NULL, &p->err);
                }
            }
        }

        if (ret->v.m == NULL) {
            matrix_error(p);
        }
    }

    return ret;
}

static int ok_matrix_dim (int r, int c, int f)
{
    if (f == F_IMAT || f == F_ZEROS || f == F_ONES ||
        f == F_MUNIF || f == F_MNORM) {
        /* zero is OK for matrix creation functions, which then
           return an empty matrix
        */
        return (r >= 0 && c >= 0);
    } else {
        double xm = (double) r * (double) c;

        return (r > 0 && c > 0 && xm < INT_MAX);
    }
}

static NODE *matrix_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    int rows = node_get_int(l, p);
    int n = 0, cols = 0;
    NODE *ret = NULL;

    if (!p->err) {
        if (f == F_RANDPERM) {
            n = rows; /* switched interpretation of first arg */
            rows = 1; /* row vector, per Matlab */
            if (null_node(r)) {
                cols = n;
            } else {
                cols = node_get_int(r, p);
            }
        } else if (f == F_IMAT && null_node(r)) {
            /* default to square */
            cols = rows;
        } else if (null_node(r)) {
            /* default to a column vector */
            cols = 1;
        } else {
            cols = node_get_int(r, p);
        }
    }

    if (!p->err && !ok_matrix_dim(rows, cols, f)) {
        p->err = E_INVARG;
        matrix_error(p);
    }

    if (!p->err) {
        ret = aux_sized_matrix_node(p, rows, cols, 0);
    }

    if (p->err || rows * cols == 0) {
        return ret;
    }

    switch (f) {
    case F_IMAT:
        if (rows != cols) {
            gretl_matrix_zero(ret->v.m);
        }
        gretl_matrix_inscribe_I(ret->v.m, 0, 0, MIN(rows, cols));
        break;
    case F_ZEROS:
        gretl_matrix_fill(ret->v.m, 0.0);
        break;
    case F_ONES:
        gretl_matrix_fill(ret->v.m, 1.0);
        break;
    case F_MUNIF:
        gretl_matrix_random_fill(ret->v.m, D_UNIFORM);
        break;
    case F_MNORM:
        gretl_matrix_random_fill(ret->v.m, D_NORMAL);
        break;
    case F_RANDPERM:
        p->err = fill_permutation_vector(ret->v.m, n);
        break;
    default:
        break;
    }

    return ret;
}

static NODE *mvnormal_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    gretl_matrix *L = l->v.m;
    int prec;

    if (L == NULL || L->rows != L->cols) {
        p->err = E_INVARG;
    } else {
        prec = node_get_bool(r, p, -1);
    }
    if (!p->err) {
        ret = aux_sized_matrix_node(p, L->rows, 1, 0);
    }
    if (!p->err && L->rows > 0) {
        p->err = correlated_normal_vec(ret->v.m, L, prec);
    }

    return ret;
}

/* Putative row or column selection matrix: must be a vector;
   cannot contain zero; cannot have both positive and negative
   entries; and entries must be integer-valued.
*/

static int set_sel_vector (matrix_subspec *spec, int r,
                           gretl_matrix *m)
{
    int i, n = gretl_vector_get_length(m);
    int err = 0;

    if (n > 0) {
        double x;
        int nneg = 0;

        for (i=0; i<n && !err; i++) {
            x = m->val[i];
            if (x == 0 || na(x) || x != floor(x)) {
                err = E_DATA;
            }
            nneg += x < 0;
        }
        if (!err && nneg > 0 && nneg < n) {
            err = E_DATA;
        }
    } else {
        err = E_TYPES;
    }

    if (err) {
        gretl_errmsg_set("Invalid selection vector");
    } else if (r) {
        spec->rsel.m = m;
        spec->rtype = SEL_MATRIX;
    } else {
        spec->lsel.m = m;
        spec->ltype = SEL_MATRIX;
    }

    return err;
}

/* Compose a sub-matrix specification, from scalars and/or
   index matrices.
*/

static void build_mspec (NODE *targ, NODE *l, NODE *r, parser *p)
{
    matrix_subspec *spec = targ->v.mspec;
    int lscalar = 0;
    int rscalar = 0;
    int i = 0, j = 0;

    if (spec == NULL) {
        spec = matrix_subspec_new();
        if (spec == NULL) {
            p->err = E_ALLOC;
            return;
        }
    }

#if EDEBUG > 1
    fprintf(stderr, "build_mspec: l->t=%d (%s)\n", l->t, getsymb(l->t));
    if (r == NULL) {
        fprintf(stderr, " r = NULL\n");
    } else {
        fprintf(stderr, " r->t=%d (%s)\n", r->t, getsymb(r->t));
    }
#endif

    /* special case: bundle membership */
    if (l->t == STR) {
        if (r == NULL) {
            spec->ltype = SEL_STR;
            spec->rtype = SEL_NULL;
            spec->lsel.str = l->v.str;
        } else {
            p->err = e_types(r);
        }
        goto finished;
    }

    lscalar = scalar_node(l);
    rscalar = (r != NULL && scalar_node(r));

    if (lscalar) {
        i = node_get_int(l, p);
        if (!p->err && i == 0) {
            gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 0);
            p->err = E_INVARG;
        }
        if (!p->err && r == NULL && i > 0) {
            /* identify and flag the single index case */
            spec->ltype = SEL_SINGLE;
            spec->rtype = SEL_NULL;
            mspec_set_row_index(spec, i);
            goto finished;
        }
    }
    if (!p->err && rscalar) {
        j = node_get_int(r, p);
        if (!p->err && j == 0) {
            gretl_errmsg_sprintf(_("Index value %d is out of bounds"), 0);
            p->err = E_INVARG;
        }
    }

    if (l->t == DUM) {
        if (r != NULL) {
            p->err = E_INVARG;
        } else {
            spec->rtype = SEL_ALL;
            if (l->v.idnum == DUM_DIAG) {
                spec->ltype = SEL_DIAG;
            } else if (l->v.idnum == DUM_UPPER) {
                spec->ltype = SEL_UPPER;
            } else if (l->v.idnum == DUM_LOWER) {
                spec->ltype = SEL_LOWER;
            } else if (l->v.idnum == DUM_REAL) {
                spec->ltype = SEL_REAL;
            } else if (l->v.idnum == DUM_IMAG) {
                spec->ltype = SEL_IMAG;
            } else {
                p->err = E_TYPES;
            }
        }
        goto finished;
    } else if (i > 0 && j > 0) {
        spec->ltype = spec->rtype = SEL_ELEMENT;
        mspec_set_row_index(spec, i);
        mspec_set_col_index(spec, j);
        goto finished;
    } else if (lscalar) {
        spec->ltype = i > 0 ? SEL_RANGE : SEL_EXCL;
        mspec_set_row_index(spec, i);
    } else if (l->t == IVEC) {
        spec->ltype = SEL_RANGE;
        spec->lsel.range[0] = l->v.ivec[0];
        spec->lsel.range[1] = l->v.ivec[1];
    } else if (l->t == MAT) {
        p->err = set_sel_vector(spec, 0, l->v.m);
    } else if (null_node(l)) {
        spec->ltype = SEL_ALL;
    } else {
        p->err = e_types(l);
        goto finished;
    }

    if (r == NULL) {
        spec->rtype = SEL_NULL;
    } else if (rscalar) {
        spec->rtype = j > 0 ? SEL_RANGE : SEL_EXCL;
        mspec_set_col_index(spec, j);
    } else if (r->t == IVEC) {
        spec->rtype = SEL_RANGE;
        spec->rsel.range[0] = r->v.ivec[0];
        spec->rsel.range[1] = r->v.ivec[1];
    } else if (r->t == MAT) {
        p->err = set_sel_vector(spec, 1, r->v.m);
    } else if (null_node(r)) {
        spec->rtype = SEL_ALL;
    } else {
        p->err = e_types(r);
    }

 finished:

#if EDEBUG > 1
    print_mspec(spec);
#endif

    if (p->err && spec != NULL) {
        free(spec);
        spec = NULL;
    }

    targ->v.mspec = spec;
}

/* node holding evaluated result of matrix specification */

static NODE *mspec_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (l->t == UNDEF || (r != NULL && r->t == UNDEF)) {
        p->err = E_DATA;
    } else {
        ret = aux_mspec_node(p);
        if (ret != NULL && starting(p)) {
            build_mspec(ret, l, r, p);
        }
    }

    return ret;
}

static NODE *submatrix_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        matrix_subspec *spec = r->v.mspec;
        gretl_matrix *m = node_get_matrix(l, p, 0, 0);

        if (!p->err) {
            if (spec->ltype == SEL_STR) {
                p->err = E_TYPES;
            } else {
                p->err = check_matrix_subspec(spec, m);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }

        if (!p->err) {
            if (spec->ltype == SEL_CONTIG) {
                ret->v.m = matrix_get_chunk(m, spec, &p->err);
            } else if (spec->ltype == SEL_ELEMENT) {
                int i = mspec_get_element(spec);

                if (m->is_complex) {
                    ret->v.m = cmatrix_get_element(m, i, &p->err);
                } else {
                    /* 2020-12-29: don't collapse to scalar here */
                    ret->v.m = gretl_matrix_alloc(1,1);
                    ret->v.m->val[0] = m->val[i];
                }
            } else {
                ret->v.m = matrix_get_submatrix(m, spec, 1, &p->err);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Check a list that has been stored in a bundle or array to see
   if it can be interpreted as a list given the characteristics
   of the current dataset (or lack thereof). Return 1 if this
   is so, otherwise 0;
*/

static int stored_list_ok (const int *list, const DATASET *dset)
{
    int i;

    if (dset == NULL || dset->n == 0) {
        return 0;
    }

    for (i=1; i<=list[0]; i++) {
        if (list[i] >= dset->v ||
            (list[i] < 0 && list[i] != LISTSEP)) {
            return 0;
        }
    }

    return 1;
}

static NODE *array_element_node (gretl_array *a, int i,
                                 parser *p)
{
    NODE *ret = NULL;
    GretlType type = 0;
    int gen_t;
    void *data;

    data = gretl_array_get_element(a, i-1, &type, &p->err);

    if (p->err == E_INVARG) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds"), i);
    }

    gen_t = gen_type_from_gretl_type(type);
    ret = get_aux_node(p, gen_t, 0, 0);

    if (!p->err) {
        if (type == GRETL_TYPE_STRING) {
            ret->v.str = data;
         } else if (type == GRETL_TYPE_MATRIX) {
            ret->v.m = data;
        } else if (type == GRETL_TYPE_BUNDLE) {
            ret->v.b = data;
        } else if (type == GRETL_TYPE_ARRAY) {
            ret->v.a = data;
        } else if (type == GRETL_TYPE_DOUBLE) {
            /* FIXME not reached? */
            ret->v.xval = *(double *) data;
        } else if (type == GRETL_TYPE_LIST) {
            /* last revised 2022-05-06 */
            if (stored_list_ok((const int *) data, p->dset)) {
                ret->v.ivec = data;
            } else {
                /* fallback: extract list as row vector */
                ret->t = MAT;
                ret->flags |= TMP_NODE;
                ret->v.m = gretl_list_to_vector((const int *) data,
                                                &p->err);
            }
        }
    }

    return ret;
}

static NODE *array_subspec_node (gretl_array *a, int *list,
                                 parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL) {
        ret->v.a = gretl_array_copy_subspec(a, list, &p->err);
    }

    return ret;
}

static NODE *list_range_node (int *list, int r1, int r2, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        ret = aux_list_node(p);
        if (ret != NULL) {
            ret->v.ivec = gretl_list_sublist(list, r1, r2);
            if (ret->v.ivec == NULL) {
                p->err = E_ALLOC;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *string_range_node (const char *s, int r1, int r2, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        ret = aux_string_node(p);
        if (ret != NULL) {
            ret->v.str = gretl_substring(s, r1, r2, &p->err);
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Return a node holding a shallow copy of series @vnum in p->dset */

static NODE *real_list_series_node (int vnum, parser *p)
{
    NODE *ret = aux_empty_series_node(p);

    if (!p->err) {
        /* overwrite the TMP_NODE flag on @ret, since we're
           borrowing a dset->Z member!
        */
        ret->flags = AUX_NODE;
        ret->vnum = vnum;
        ret->v.xvec = p->dset->Z[vnum];
    }

    return ret;
}

/* Handling a context where we have a list under node @l,
   and position within the list in @pos.
*/

static NODE *list_member_by_pos (NODE *l, int pos, parser *p)
{
    NODE *ret = NULL;
    int *list = l->v.ivec;

#if EDEBUG
    fprintf(stderr, "list_member_by_pos: pos = %d\n", pos);
    printlist(list, "list");
#endif

    if (pos < 1 || pos > list[0]) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds"), pos);
        p->err = E_INVARG;
    } else {
        int v = list[pos];

        if (v < 0 || v >= p->dset->v) {
            gretl_errmsg_sprintf(_("Variable number %d is out of bounds"), v);
            p->err = E_INVARG;
        } else {
            ret = real_list_series_node(v, p);
        }
    }

    return ret;
}

/* Handling a context where we have a list under node @l,
   and the name of a (putative) list member under node @r.
*/

static NODE *list_member_by_name (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    int *list = l->v.ivec;
    const char *s = r->v.str;
    int i, v;

    if (p->dset == NULL) {
        p->err = E_DATA;
        return NULL;
    }

#if EDEBUG
    fprintf(stderr, "list_member_by_name: %s.%s\n", l->vname, s);
#endif

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        if (!strcmp(s, p->dset->varname[v])) {
            /* found the named member series */
            ret = real_list_series_node(v, p);
            break;
        }
    }

    if (ret == NULL) {
        p->err = E_UNKVAR;
    }

    return ret;
}

static int mspec_get_series_index (matrix_subspec *s,
                                   parser *p)
{
    int t = -1;

    if (s->ltype == SEL_SINGLE) {
        t = s->lsel.range[0];
    } else if (s->ltype == SEL_RANGE && s->rtype == SEL_NULL) {
        if (s->lsel.range[0] == s->lsel.range[1]) {
            t = s->lsel.range[0];
        } else {
            /* allow for dates such as "2008:4" */
            gchar *tmp;

            tmp = g_strdup_printf("%d:%d", s->lsel.range[0],
                                  s->lsel.range[1]);
            t = get_observation_number(tmp, p->dset);
            g_free(tmp);
        }
    }

    if (t < 1 || t > p->dset->n) {
        p->err = E_DATA;
    }

    return t;
}

static int mspec_get_array_index (matrix_subspec *spec,
                                  int *err)
{
    int idx = 0;

    if (spec->ltype == SEL_SINGLE) {
        idx = spec->lsel.range[0];
    } else if (spec->ltype == SEL_RANGE &&
        spec->rtype == SEL_NULL &&
        spec->lsel.range[0] == spec->lsel.range[1]) {
        idx = spec->lsel.range[0];
    } else {
        gretl_errmsg_set("Invalid left-hand side index value");
        *err = E_TYPES;
    }

    return idx;
}

/* stricter variant of test_for_single_range */

static int get_single_element (matrix_subspec *spec,
                               parser *p)
{
    int ret = 0;

    if (spec->ltype == SEL_SINGLE) {
        ret = spec->lsel.range[0];
    } else {
        if (p != NULL) {
            p->err = E_TYPES;
        }
        ret = -1;
    }

    return ret;
}

static void *sub_addr_get_data (NODE *t, GretlType *ptype,
                                user_var **puv)
{
    NODE *l = t->L, *r = t->R;
    gretl_array *a = l->v.ptr;
    GretlType type = 0;
    void *elem;
    int idx, err = 0;

    idx = get_single_element(r->v.mspec, NULL);
    elem = gretl_array_get_element(a, idx-1, &type, &err);
    *ptype = (elem == NULL)? GRETL_TYPE_NONE :
        gretl_type_get_ref_type(type);
    *puv = l->uv;
#if 0
    fprintf(stderr, "sub_addr_get_data: idx=%d, a=%p, type=%d, err=%d\n",
            idx, a, type, err);
#endif

    return elem;
}

static NODE *process_OSL_address (NODE *t, NODE *l, NODE *r, parser *p)
{
    int idx = get_single_element(r->v.mspec, p);
    NODE *lb = l->L;
    NODE *ret = NULL;

    if (lb->t != OSL || lb->uv == NULL || idx <= 0) {
        p->err = E_TYPES;
    } else {
        GretlType type = user_var_get_type(lb->uv);

        if (type == GRETL_TYPE_ARRAY) {
            ret = aux_parent_node(p);
            if (ret != NULL) {
                ret->t = SUB_ADDR;
                ret->L = lb; /* extracted left-hand */
                ret->R = r;  /* evaluated right-hand */
                /* prevent double-freeing of children @l and @r */
                ret->flags |= LHT_NODE;
            }
        } else {
            p->err = E_TYPES;
        }
    }

    if (p->err) {
        gretl_errmsg_set(_("Wrong type of operand for unary '&'"));
    }

    return ret;
}

static int want_singleton_array (NODE *n, parser *p)
{
    if (p->aux != NULL && p->aux->t == ARRAY) {
        GretlType t = gretl_array_get_type(n->v.a);

        /* We want to preserve the ARRAY type of the aux
           node associated with @n: this requires producing
           a singleton array unless we're looking at an
           array of arrays.
        */
        return t != GRETL_TYPE_ARRAYS;
    } else if (p->aux == NULL && p->targ == ARRAY &&
               p->tree == n->parent) {
        return 1;
    }

    return 0;
}

static int *array_subspec_list (NODE *l, NODE *r, parser *p)
{
    matrix_subspec *spec = r->v.mspec;
    int *list = NULL;

    if (spec->rtype != SEL_NULL) {
        /* array selection must be one-dimensional */
        p->err = E_INVARG;
    } else {
        int len;

        if (l->t == ARRAY) {
            len = gretl_array_get_length(l->v.a);
        } else if (l->t == STR) {
            len = g_utf8_strlen(l->v.str, -1);
        } else if (l->t == LIST) {
            len = l->v.ivec[0];
        } else {
            p->err = e_types(l);
            return NULL;
        }

        /* convert @spec to list of elements */
        list = mspec_make_list(spec->ltype, &spec->lsel,
                               len, &p->err);
    }

    return list;
}

static NODE *subobject_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        if (r == NULL || r->t != MSPEC) {
            p->err = e_types(r);
        } else if (l->t == MAT) {
            ret = submatrix_node(l, r, p);
        } else if (l->t == ARRAY) {
            int *vlist = array_subspec_list(l, r, p);

            if (!p->err && vlist[0] == 1) {
                if (want_singleton_array(l, p)) {
                    /* produce a 1-element array */
                    ret = array_subspec_node(l->v.a, vlist, p);
                } else {
                    /* extract an array element */
                    ret = array_element_node(l->v.a, vlist[1], p);
                }
            } else if (!p->err) {
                ret = array_subspec_node(l->v.a, vlist, p);
            }
            free(vlist);
        } else if (l->t == LIST || l->t == STR) {
            int *vlist = array_subspec_list(l, r, p);

            if (!p->err && vlist[0] == 1 && l->t == LIST) {
                /* singleton list selection */
                ret = list_member_by_pos(l, vlist[1], p);
            } else if (!p->err && gretl_list_is_consecutive(vlist)) {
                /* selected elements are consecutive */
                int r1 = vlist[1];
                int r2 = vlist[vlist[0]];

                if (l->t == LIST) {
                    ret = list_range_node(l->v.ivec, r1, r2, p);
                } else {
                    ret = string_range_node(l->v.str, r1, r2, p);
                }
            } else if (!p->err) {
                /* selection is not consecutive */
                if (l->t == STR) {
                    ret = aux_string_node(p);
                    ret->v.str = gretl_utf8_select(l->v.str, vlist);
                } else {
                    ret = aux_list_node(p);
                    ret->v.ivec = gretl_list_select(l->v.ivec, vlist, &p->err);
                }
            }
            free(vlist);
        } else if (l->t == SERIES) {
            int t = mspec_get_series_index(r->v.mspec, p);

            if (!p->err) {
                ret = aux_scalar_node(p);
                if (!p->err) {
                    ret->v.xval = l->v.xvec[t-1];
                }
            }
        } else if (l->t == BUNDLE) {
            /* the "mspec" must hold a single key string */
            const char *key = mspec_get_string(r->v.mspec, 0);
            GretlType type = GRETL_TYPE_NONE;
            void *val = NULL;
            int tmp = 0;
            int size = 0;

            if (key == NULL) {
                p->err = E_TYPES;
            } else {
                val = gretl_bundle_get_element(l->v.b, key, &type,
                                               &size, &tmp, &p->err);
            }
            if (!p->err) {
                int t = gen_type_from_gretl_type(type);

                if (t == NUM) {
                    ret = aux_scalar_node(p);
                    if (!p->err) {
                        ret->v.xval = *(double *) val;
                    }
                } else {
                    ret = get_aux_node(p, t, 0, tmp ? TMP_NODE : 0);
                    if (!p->err) {
                        ret->v.ptr = val;
                    }
                }
            }
        } else if (l->t == NUM) {
            /* allow "indexing into" a scalar, but only for a
               single index with value 1
            */
            int i = get_single_element(r->v.mspec, p);

            if (i == 1) {
                ret = aux_scalar_node(p);
                if (!p->err) {
                    ret->v.xval = l->v.xval;
                }
            } else {
                p->err = E_TYPES;
            }
        } else {
            p->err = e_types(l);
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *process_subslice (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        if (!scalar_node(l)) {
            p->err = e_types(l);
        } else if (!null_or_scalar(r)) {
            p->err = e_types(r);
        } else {
            ret = aux_ivec_node(p, 2);
            if (ret != NULL) {
                ret->v.ivec[0] = node_get_int(l, p);
                if (null_node(r)) {
                    ret->v.ivec[1] = MSEL_MAX; /* placeholder */
                } else {
                    ret->v.ivec[1] = node_get_int(r, p);
                }
            }
        }
    } else {
        ret = aux_ivec_node(p, 2);
    }

    return ret;
}

/* Note: many standard and a few non-standard math functions
   are not included in the switch below, because pointers to
   the functions are saved, obviating the need for repeated
   lookup. See the @ptrfuncs mechanism in genlex.c.
*/

static double real_apply_func (double x, int f, parser *p)
{
    double y;

    errno = 0;

    if (na(x)) {
        switch (f) {
        case F_MISSING:
            return 1.0;
        case F_DATAOK:
        case F_MISSZERO:
            return 0.0;
        default:
            return NADBL;
        }
    }

    switch (f) {
    case U_NEG:
        return -x;
    case U_POS:
        return x;
    case U_NOT:
        return x == 0;
    case F_TOINT:
        if (x > G_MAXINT32 || x < G_MININT32) {
            return NADBL;
        } else {
            return (double) (int) x;
        }
    case F_MISSING:
        return 0.0;
    case F_DATAOK:
        return 1.0;
    case F_MISSZERO:
        return x;
    case F_ZEROMISS:
        return (x == 0.0)? NADBL : x;
    case F_EASTER:
        y = easterdate(x);
        return y;
        /* below: functions that should already be mapped;
           it should be possible to delete them
        */
    case F_LOG: /* in case it's aliased */
        return log(x);
    default:
        return 0.0;
    }
}

/* @n must be of type NUM, MAT or SERIES, pre-checked */

static double node_get_double (NODE *n, int i, parser *p)
{
    if (n->t == NUM) {
        return n->v.xval;
    } else if (n->t == MAT) {
        return n->v.m->val[i];
    } else {
        return n->v.xvec[p->dset->t1 + i];
    }
}

/* @n must be of type NUM, MAT or SERIES, pre-checked */

static void node_set_double (NODE *n, int i, double x, parser *p)
{
    if (n->t == NUM) {
        n->v.xval = x;
    } else if (n->t == MAT) {
        n->v.m->val[i] = x;
    } else {
        n->v.xvec[p->dset->t1 + i] = x;
    }
}

static int is_integral (double d)
{
    return d == floor(d) || d == ceil(d);
}

static double bincoeff (double n, double k, int *err)
{
    double ret;

    if (n < k || k < 0) {
        *err = E_INVARG;
        return NADBL;
    }

    /* catch special cases first */
    if (n == k || k == 0) {
        ret = 1.0;
    } else if (n - k == 1 || k == 1) {
        ret = n;
    } else {
        ret = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
        if (is_integral(n) && is_integral(k)) {
            ret = nearbyint(exp(ret));
        } else {
            ret = exp(ret);
        }
    }

    return ret;
}

/* flexible_2arg_node() handles cases like atan2, where we have two
   possibly heterogeneous arguments (scalar, series, matrix) and the
   objective is to return a sensibly sized object.
*/

static NODE *flexible_2arg_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    int rettype = 0;
    int nmin, nmax;
    int nl = 0, nr = 0;

    if (l->t == NUM) {
        nl = 1;
    } else if (l->t == MAT) {
        nl = gretl_vector_get_length(l->v.m);
    } else {
        nl = sample_size(p->dset);
    }

    if (r->t == NUM) {
        nr = 1;
    } else if (r->t == MAT) {
        nr = gretl_vector_get_length(r->v.m);
    } else {
        nr = sample_size(p->dset);
    }

    nmin = nr < nl ? nr : nl;
    nmax = nr > nl ? nr : nl;

    if (nmin == 0 || (nmin > 1 && nmax != nmin)) {
        p->err = E_NONCONF;
        return NULL;
    }

    /* ordering is MAT > SERIES > NUM */
    rettype = r->t > l->t ? r->t : l->t;

    if (rettype == NUM) {
        ret = aux_scalar_node(p);
    } else if (rettype == MAT) {
        ret = aux_sized_matrix_node(p, nmax, 1, 0);
    } else {
        ret = aux_series_node(p);
    }

    if (ret != NULL && !p->err) {
        double x1, x2, y = NADBL;
        int i;

        for (i=0; i<nmax; i++) {
            x1 = node_get_double(l, i, p);
            x2 = node_get_double(r, i, p);
            if (f == F_ATAN2) {
                y = atan2(x1, x2);
            } else if (f == F_BINCOEFF) {
                y = bincoeff(x1, x2, &p->err);
                if (p->err) {
                    break;
                }
            }
            node_set_double(ret, i, y, p);
        }
    }

    return ret;
}

static NODE *apply_scalar_func (NODE *n, NODE *f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        double (*dfunc) (double) = f->v.ptr;

        if (dfunc != NULL) {
            ret->v.xval = dfunc(n->v.xval);
        } else {
            ret->v.xval = real_apply_func(n->v.xval, f->t, p);
        }
    }

    return ret;
}

static NODE *sleep_node (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        double secs = node_get_scalar(n, p);

        if (secs < 0 || secs * G_USEC_PER_SEC >= (double) G_MAXULONG) {
            p->err = E_INVARG;
        } else {
            g_usleep(G_USEC_PER_SEC * secs);
            ret->v.xval = 0;
        }
    }

    return ret;
}

static NODE *scalar_isnan_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        double x = node_get_scalar(n, p);

        if (!p->err) {
            ret->v.xval = isnan(x) != 0;
        }
    }

    return ret;
}

static NODE *matrix_isnan_node (NODE *n, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const gretl_matrix *m = n->v.m;

        if (m->rows == 0 || m->cols == 0) {
            p->err = E_DATA;
        } else {
            ret->v.m = gretl_matrix_alloc(m->rows, m->cols);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
            } else {
                int i, n = m->rows * m->cols;

                for (i=0; i<n; i++) {
                    ret->v.m->val[i] = isnan(m->val[i]) != 0;
                }
                gretl_matrix_set_complex(ret->v.m, m->is_complex);
            }
        }
    }

    return ret;
}

static NODE *apply_series_func (NODE *n, NODE *f, parser *p)
{
    NODE *ret = aux_series_node(p);
    int t;

    if (ret != NULL) {
        double (*dfunc) (double) = f->v.ptr;
        const double *x;

        if (n->t == SERIES) {
            if (n->v.xvec == NULL) {
                fprintf(stderr, "apply_series_func: got NULL xvec\n");
                p->err = E_DATA;
            } else {
                x = n->v.xvec;
            }
        } else {
            x = get_colvec_as_series(n, f->t, p);
        }

        if (!p->err) {
            if (autoreg(p)) {
                if (dfunc != NULL) {
                    ret->v.xvec[p->obs] = dfunc(x[p->obs]);
                } else {
                    ret->v.xvec[p->obs] = real_apply_func(x[p->obs], f->t, p);
                }
            } else if (dfunc != NULL) {
                for (t=p->dset->t1; t<=p->dset->t2; t++) {
                    ret->v.xvec[t] = dfunc(x[t]);
                }
            } else {
                for (t=p->dset->t1; t<=p->dset->t2; t++) {
                    ret->v.xvec[t] = real_apply_func(x[t], f->t, p);
                }
            }
        }
    }

    return ret;
}

/* argument is series or list; value returned is list */

static NODE *dummify_func (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list = NULL;
        double oddval = NADBL;

        if (!null_node(r)) {
            if (r->t != NUM) {
                p->err = e_types(r);
                return ret;
            } else {
                oddval = r->v.xval;
            }
        }

        if (l->t == LIST) {
            list = gretl_list_copy(l->v.ivec);
        } else if (useries_node(l)) {
            list = gretl_list_new(1);
            list[1] = l->vnum;
        } else {
            gretl_errmsg_set(_("The first argument must be a named series "
                               "in the current dataset"));
            p->err = E_INVARG;
        }

        if (p->err) {
            ; /* don't do anything more */
        } else if (list == NULL) {
            p->err = E_ALLOC;
        } else if (null_node(r)) {
            /* got just one argument */
            p->err = list_dumgenr(&list, p->dset, OPT_F);
            ret->v.ivec = list;
        } else if (list[0] > 1) {
            gretl_errmsg_set("dummify(x, y): first argument should be a single series");
            free(list);
            p->err = E_DATA;
        } else {
            p->err = dumgenr_with_oddval(&list, p->dset, oddval);
            ret->v.ivec = list;
        }
    }

    return ret;
}

/* argument is list; value returned is list */

static NODE *cdummify_func (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list = NULL;

        if (n->t == LIST) {
            list = gretl_list_copy(n->v.ivec);
        } else if (useries_node(n)) {
            list = gretl_list_new(1);
            list[1] = n->vnum;
        } else {
            p->err = e_types(n);
        }

        if (p->err) {
            ; /* don't do anything more */
        } else if (list == NULL) {
            p->err = E_ALLOC;
        } else {
            p->err = auto_dummify_list(&list, p->dset);
            ret->v.ivec = list;
        }
    }

    return ret;
}

static NODE *get_info_on_series (NODE *n, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL && starting(p)) {
        int v = 0;

        if (useries_node(n)) {
            v = n->vnum;
        } else {
            v = node_get_int(n, p);
        }

        if (!p->err) {
            ret->v.b = series_info_bundle(p->dset, v, &p->err);
        }
    }

    return ret;
}

static NODE *list_stdize (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list = NULL;
        int dfc = 1;

        if (!null_node(r)) {
            dfc = node_get_int(r, p);
        }
        if (!p->err) {
            list = node_get_list(l, p);
        }
        if (!p->err) {
            gretlopt opt;

            opt = dfc < 0 ? OPT_C : dfc == 0 ? OPT_N : OPT_NONE;
            if (list[0] > 0) {
                p->err = list_stdgenr(list, p->dset, opt);
            }
            ret->v.ivec = list;
        }
    }

    return ret;
}

static NODE *series_stdize (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int dfc = 1;

        if (!null_node(r)) {
            dfc = node_get_int(r, p);
        }
        if (!p->err) {
            p->err = standardize_series(l->v.xvec, ret->v.xvec,
                                        dfc, p->dset);
        }
    }

    return ret;
}

/* middle argument is series or list; value returned is list in
   either case */

static NODE *list_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        gretlopt opt = OPT_NONE;
        gretl_matrix *v = NULL;
        int *list = NULL;
        int k = 0;

        /* ordering of the results? */
        if (node_get_bool(r, p, 0) > 0 && !p->err) {
            opt = OPT_L; /* by lags */
        }

        if (!p->err) {
            /* scalar lag order or vector */
            if (l->t == NUM) {
                k = node_get_int(l, p);
            } else {
                v = l->v.m;
            }
        }

        if (!p->err) {
            list = node_get_list(m, p);
        }

        if (!p->err) {
            if (list[0] > 0) {
                p->err = list_laggenr(&list, 1, k, v,
                                      p->dset, 0, opt);
            }
            ret->v.ivec = list;
        }
    }

    return ret;
}

static NODE *matrix_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        gretlopt opt = OPT_NONE;
        gretl_matrix *kvec = NULL;
        gretl_matrix *src = NULL;

        if (node_get_bool(r, p, 0) > 0 && !p->err) {
            opt = OPT_L; /* by lags */
        }

        if (!p->err) {
            /* scalar max lag order or vector */
            if (l->t == NUM) {
                int i, k = node_get_int(l, p);

                if (!p->err && k <= 0) {
                    /* scalar k must be positive */
                    p->err = E_INVARG;
                } else if (!p->err) {
                    kvec = gretl_vector_alloc(k);
                    if (kvec == NULL) {
                        p->err = E_ALLOC;
                    } else {
                        for (i=0; i<k; i++) {
                            kvec->val[i] = i+1;
                        }
                    }
                }
            } else {
                kvec = l->v.m;
            }
        }

        if (!p->err) {
            src = m->v.m;
        }

        if (!p->err) {
            ret->v.m = gretl_matrix_lag(src, kvec, opt, 0.0);
        }

        if (kvec != l->v.m) {
            gretl_matrix_free(kvec);
        }
    }

    return ret;
}

/* args are minlag, maxlag, MIDAS-list-to-be-lagged */

static NODE *hf_list_make_lags (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list = NULL;
        int lmin, lmax = 0;
        int cfac = 0;

        lmin = node_get_int(l, p);
        if (!p->err) {
            lmax = node_get_int(m, p);
        }
        if (!p->err) {
            list = node_get_list(r, p);
        }

        if (!p->err) {
            /* compaction factor for high-frequency data */
            cfac = list[0];
            if (cfac < 2) {
                fprintf(stderr, "hflags: not a MIDAS list\n");
                p->err = E_INVARG;
                free(list);
            }
        }

        if (!p->err) {
            if (list[0] > 0) {
                p->err = list_laggenr(&list, lmin, lmax, NULL,
                                      p->dset, cfac, OPT_L);
            }
            ret->v.ivec = list;
        }
    }

    return ret;
}

/* In case we get either of the two optional arguments to
   ldiff(), implement them via post-processing.
*/

static void post_process_ldiff (NODE *ret, NODE *m, NODE *r,
                                parser *p)
{
    double d = NADBL;
    int fill, mul = 0;

    if (ret->t != SERIES && ret->t != LIST) {
        /* "can't happen" */
        p->err = E_DATA;
        return;
    }

    fill = node_get_bool(m, p, 0);
    if (!p->err && !null_node(r)) {
        d = node_get_scalar(r, p);
        mul = 1;
    }

    if (!p->err) {
        double *x = NULL;
        int *list = NULL;
        int imax = 1;
        int i, t;

        if (ret->t == SERIES) {
            x = ret->v.xvec;
        } else {
            list = ret->v.ivec;
            imax = list[0];
        }
        for (i=0; i<imax; i++) {
            if (list != NULL) {
                x = p->dset->Z[list[i+1]];
            }
            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                if (na(x[t])) {
                    if (t > 0 && fill) {
                        x[t] = 0.0;
                    }
                } else if (mul) {
                    x[t] *= d;
                }
            }
        }
    }
}

#define ok_list_func(f) (f == F_LOG || f == F_DIFF || \
                         f == F_LDIFF || f == F_SDIFF || \
                         f == F_SQUARE || f == F_ODEV || \
                         f == F_RESAMPLE || f == F_DROPCOLL || \
                         f == F_HFDIFF || f == F_HFLDIFF)

/* The following handles functions that are "basically" for series,
   but which can also be applied to lists -- except for F_DROPCOLL,
   F_HFDIFF and F_HDLDIFF, which require a list argument.
*/

static NODE *apply_list_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (!ok_list_func(f)) {
        p->err = E_TYPES;
        return ret;
    }

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(n, p);
        gretlopt opt = OPT_NONE;
        double parm = NADBL;
        int t = 0;

        if (f == F_SQUARE) {
            if (r != NULL && node_is_true(r, p)) {
                opt = OPT_O;
            }
        } else if (f == F_DROPCOLL || f == F_HFDIFF ||
                   f == F_HFLDIFF) {
            /* handle optional parameter */
            if (!null_node(r)) {
                parm = node_get_scalar(r, p);
                if (p->err) {
                    return ret;
                }
            }
        }

        /* note: @list is modified by the library functions
           called below */

        if (list != NULL) {
            /* note: an empty list argument produces an
               empty list return
            */
            if (list[0] > 0) {
                switch (f) {
                case F_LOG:
                    p->err = list_loggenr(list, p->dset);
                    break;
                case F_DIFF:
                case F_LDIFF:
                case F_SDIFF:
                    if (f == F_DIFF) t = DIFF;
                    else if (f == F_LDIFF) t = LDIFF;
                    else if (f == F_SDIFF) t = SDIFF;
                    p->err = list_diffgenr(list, t, p->dset);
                    break;
                case F_SQUARE:
                    p->err = list_xpxgenr(&list, p->dset, opt);
                    break;
                case F_ODEV:
                    p->err = list_orthdev(list, p->dset);
                    break;
                case F_RESAMPLE:
                    p->err = list_resample(list, p->dset);
                    break;
                case F_DROPCOLL:
                    p->err = list_dropcoll(list, parm, p->dset);
                    break;
                case F_HFDIFF:
                case F_HFLDIFF:
                    t = (f == F_HFDIFF)? DIFF : LDIFF;
                    p->err = hf_list_diffgenr(list, t, parm, p->dset);
                    break;
                default:
                    break;
                }
            }
            ret->v.ivec = list;
        }
    }

    return ret;
}

static NODE *hf_list_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    gretl_matrix *v = l->v.m;
    int f_ratio = node_get_int(m, p);
    char *pfx = r->v.str;
    NODE *ret = NULL;

    if (!p->err) {
        int n = gretl_vector_get_length(v);

        if (n == 0) {
            p->err = E_NONCONF;
        } else if (f_ratio < 3) {
            p->err = E_INVARG;
        } else if (*pfx == '\0' || !gretl_is_ascii(pfx) ||
                   strlen(pfx) > 24) {
            p->err = E_INVARG;
        } else {
            int T = sample_size(p->dset);

            if (n != f_ratio * T) {
                p->err = E_INVARG;
            }
        }
    }

    if (!p->err) {
        ret = aux_list_node(p);
        if (ret != NULL) {
            ret->v.ivec = vector_to_midas_list(v, f_ratio, pfx,
                                               p->dset, &p->err);
        }
    }

    return ret;
}

static NODE *dataset_list_node (parser *p)
{
    NODE *ret = NULL;

    if (gretl_function_depth() > 0) {
        gretl_errmsg_set("'dataset' is not recognized as a list within functions");
        p->err = E_DATA;
    } else {
        ret = aux_list_node(p);
    }

    if (ret != NULL && starting(p)) {
        int *list = full_var_list(p->dset, NULL);

        if (list == NULL) {
            list = gretl_null_list();
        }
        if (list == NULL) {
            p->err = E_DATA;
        }
        ret->v.ivec = list;
    }

    return ret;
}

static NODE *trend_node (parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        ret = aux_empty_series_node(p);
        if (!p->err) {
            p->err = gen_time(p->dset, 1, &ret->vnum);
            if (!p->err) {
                ret->v.xvec = p->dset->Z[ret->vnum];
                /* not TMP_NODE because we're borrowing a Z column */
                ret->flags &= ~TMP_NODE;
            }
        }
    }

    return ret;
}

static int object_get_size (NODE *n, int context, parser *p)
{
    int ret = 0;

    if (n->t == NUM) {
        ret = 1;
    } else if (n->t == MAT) {
        ret = n->v.m->rows * n->v.m->cols;
    } else if (n->t == ARRAY) {
        ret = gretl_array_get_length(n->v.a);
    } else if (n->t == BUNDLE) {
        ret = gretl_bundle_get_n_members(n->v.b);
    } else if (n->t == LIST) {
        ret = n->v.ivec[0];
    } else if (n->t == STR) {
        if (context == DUM_END) {
            ret = g_utf8_strlen(n->v.str, -1);
        } else {
            ret = strlen(n->v.str);
        }
    } else {
        p->err = e_types(n);
    }

    return ret;
}

#define END_DEBUG 0

/* Cash out 'end' where it is legit, namely indicating the last
   element of a matrix row or column, or an array, list or string. We
   need to find the object to which 'end' is supposed to refer (@obj),
   and we also need to determine whether 'end' occurs in the first or
   second of the two possible indexation slots under the SLRAW node in
   its parentage. The second slot (idx = 2) is valid only as the column
   index for a matrix.
*/

static int object_end_index (NODE *t, parser *p)
{
    NODE *pa = t->parent;
    NODE *obj = NULL;
    int idx = 0;
    int ret = -1;

#if END_DEBUG
    fprintf(stderr, "\nobject_end_index: t=%p\n", (void *) t);
#endif

    while (pa != NULL) {
        /* examine the parentage of @t */
        if (pa->t == SLRAW) {
            if (find_in_tree(pa->R, t)) {
                idx = 2;
            } else if (pa->R != NULL && find_in_tree(pa->L, t)) {
                idx = 1;
            }
        } else if (pa->t == OSL) {
            obj = (pa->L->aux != NULL)? pa->L->aux : pa->L;
            break;
        }
        pa = pa->parent;
    }

#if END_DEBUG
    fprintf(stderr, " objtype = %s, idx = %d\n",
            obj == NULL ? "none" : getsymb(obj->t), idx);
#endif

    if (obj == NULL || (idx != 0 && obj->t != MAT)) {
        /* either we didn't find a referent, or it's a
           case of invalid double indexation
        */
        p->err = E_INVARG;
    } else if (obj->t == MAT) {
        if (idx == 2) {
            /* index must refer to column */
            ret = obj->v.m->cols;
        } else if (idx == 1) {
            /* index must refer to row */
            ret = obj->v.m->rows;
        } else {
            /* singleton index: may refer to row or column
               if the MAT is a vector; otherwise the error
               will be caught later
            */
            ret = gretl_vector_get_length(obj->v.m);
            if (ret == 0) {
                ret = obj->v.m->rows;
            }
        }
    } else if (obj->t == ARRAY || obj->t == LIST || obj->t == STR) {
        ret = object_get_size(obj, DUM_END, p);
    } else {
        p->err = E_INVARG;
    }

    if (p->err) {
        gretl_errmsg_set("'end' is not defined in this context");
    }

    return ret;
}

static NODE *array_last_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        ret = aux_scalar_node(p);
        if (!p->err) {
            ret->v.xval = object_end_index(t, p);
            if (p->err) {
                ret->v.xval = NADBL;
            }
        }
    }

    return ret;
}

static NODE *seasonals_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (!dataset_is_seasonal(p->dset) &&
        !dataset_is_seasonal_panel(p->dset)) {
        p->err = E_PDWRONG;
    } else {
        int ref = 0, center = 0;

        if (!null_node(l)) {
            ref = node_get_int(l, p);
        }
        if (!null_node(r)) {
            center = node_is_true(r, p);
        }
        if (!p->err) {
            ret = aux_list_node(p);
        }
        if (ret != NULL) {
            ret->v.ivec = seasonals_list(p->dset, ref, center, &p->err);
        }
    }

    return ret;
}

static NODE *get_lag_list (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        int *list = NULL, *srclist = NULL;
        int i, imin = 1, imax = 1;
        int lv = 0;

        if (!useries_node(l) && l->t != LIST) {
            /* we need a named series or a list on the left */
            p->err = e_types(l);
        } else if (r->t != IVEC && r->t != NUM) {
            /* we need one or more integers on the right */
            p->err = e_types(r);
        }

        if (p->err) {
            return NULL;
        }

        if (l->t == LIST) {
            srclist = l->v.ivec;
            imax = srclist[0];
        } else {
            lv = l->vnum;
        }

        if (imax == 0) {
            /* empty list on input -> empty on output */
            list = gretl_list_new(0);
            if (list == NULL) {
                p->err = E_ALLOC;
            }
            goto loopdone;
        }

        for (i=imin; i<=imax && !p->err; i++) {
            if (srclist != NULL) {
                lv = srclist[i];
            }
            if (r->t == IVEC) {
                int fromlag = r->v.ivec[0];
                int tolag = r->v.ivec[1];

                if (list == NULL) {
                    list = laggenr_from_to(lv, fromlag, tolag,
                                           p->dset, &p->err);
                } else {
                    int *tmp;

                    tmp = laggenr_from_to(lv, fromlag, tolag,
                                           p->dset, &p->err);
                    if (!p->err) {
                        p->err = gretl_list_append_list(&list, tmp);
                        free(tmp);
                    }
                }
            } else {
                int lag = -r->v.xval;

                lv = laggenr(lv, lag, p->dset);
                if (lv > 0) {
                    list = gretl_list_append_term(&list, lv);
                    if (list == NULL) {
                        p->err = E_ALLOC;
                    }
                }
            }
        }

    loopdone:

        if (list != NULL) {
            ret = aux_list_node(p);
            if (ret != NULL) {
                ret->v.ivec = list;
            } else {
                free(list);
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

int *list_from_strings_array (gretl_array *a, parser *p)
{
    GretlType type = gretl_array_get_type(a);
    int *list = NULL;

    if (type != GRETL_TYPE_STRINGS) {
        p->err = E_TYPES;
    } else {
        int i, vi, n = 0;
        char **S = gretl_array_get_strings(a, &n);

        for (i=0; i<n && !p->err; i++) {
            vi = current_series_index(p->dset, S[i]);
            if (vi < 0) {
                gretl_errmsg_sprintf("'%s' is not a known series", S[i]);
                p->err = E_UNKVAR;
            }
        }

        if (!p->err) {
            list = gretl_list_new(n);
            if (list == NULL) {
                p->err = E_ALLOC;
            } else {
                for (i=0; i<n; i++) {
                    list[i+1] = current_series_index(p->dset, S[i]);
                }
            }
        }
    }

    return list;
}

/* get an *int LIST from node @n: note that the list is always
   newly allocated, and so should be freed by the caller if
   it's just for temporary use
*/

int *node_get_list (NODE *n, parser *p)
{
    int *list = NULL;
    int v = 0;

    if (null_node(n)) {
        if (p->flags & P_PRNLIST) {
            p->err = e_types(n);
        } else {
            list = gretl_null_list();
        }
    } else if (n->t == LIST) {
        list = gretl_list_copy(n->v.ivec);
    } else if (n->t == SERIES || n->t == NUM) {
        v = (n->t == SERIES)? n->vnum : node_get_int(n, p);
        if (!p->err) {
            if (v < 0 || v >= p->dset->v) {
                p->err = E_UNKVAR;
            } else {
                list = gretl_list_new(1);
                if (list == NULL) {
                    p->err = E_ALLOC;
                } else {
                    list[1] = v;
                }
            }
        }
    } else if (dataset_dum(n)) {
        list = full_var_list(p->dset, NULL);
    } else if (n->t == MAT) {
        list = gretl_list_from_vector(n->v.m, p->dset, &p->err);
    } else {
        p->err = e_types(n);
    }

    if (!p->err && list == NULL) {
        p->err = E_ALLOC;
    } else if (p->err == E_UNKVAR && v != 0) {
        gretl_errmsg_sprintf(_("Variable number %d is out of bounds"), v);
    }

    return list;
}

static NODE *eval_lcat (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list1, *list2 = NULL;

        list1 = node_get_list(l, p); /* note, copied */
        if (list1 != NULL) {
            list2 = node_get_list(r, p); /* copied */
        }
        if (list2 != NULL) {
            p->err = gretl_list_append_list(&list1, list2);
        }
        ret->v.ivec = list1;
        free(list2);
    }

    return ret;
}

static NODE *list_list_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *llist, *rlist = NULL;
        int *list = NULL;

        llist = node_get_list(l, p);
        if (llist != NULL) {
            rlist = node_get_list(r, p);
        }
        if (rlist != NULL) {
            if (f == B_AND) {
                list = gretl_list_intersection(llist, rlist, &p->err);
            } else if (f == B_OR) {
                list = gretl_list_union(llist, rlist, &p->err);
            } else if (f == B_SUB) {
                list = gretl_list_drop(llist, rlist, &p->err);
            } else if (f == B_POW) {
                list = gretl_list_product(llist, rlist, p->dset, &p->err);
            } else if (f == B_ADD) {
                list = gretl_list_plus(llist, rlist, &p->err);
            }
        }
        ret->v.ivec = list;
        free(llist);
        free(rlist);
    }

    return ret;
}

/* Binary operator applied to two bundles: at present only '+'
   (for union) is supported.
*/

static NODE *bundle_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_bundle *bl = l->v.b;
        gretl_bundle *br = r->v.b;

        if (f == B_ADD) {
            ret = aux_bundle_node(p);
        } else if (f == B_EQ || f == B_NEQ) {
            ret = aux_scalar_node(p);
        } else {
            p->err = E_TYPES;
        }
        if (!p->err && f == B_ADD) {
            ret->v.b = gretl_bundle_union(bl, br, &p->err);
        } else if (!p->err) {
            int eq = gretl_bundles_are_equal(bl, br);

            ret->v.xval = (f == B_EQ)? eq : !eq;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

/* Binary operator applied to two arrays: '+' (append),
   '||' (union) or '&&' (intersection). But the latter
   two are only for strings arrays.
*/

static NODE *array_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_array *al = l->v.a;
        gretl_array *ar = r->v.a;

        if (!p->err) {
            if (f == B_EQ || f == B_NEQ) {
                ret = aux_scalar_node(p);
            } else if (f == B_ADD || f == B_OR || f == B_AND) {
                ret = aux_array_node(p);
            } else {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            if (f == B_ADD) {
                ret->v.a = gretl_arrays_join(al, ar, &p->err);
            } else if (f == B_OR) {
                ret->v.a = gretl_arrays_union(al, ar, &p->err);
            } else if (f == B_AND) {
                ret->v.a = gretl_arrays_intersection(al, ar, &p->err);
            } else {
                int eq = gretl_arrays_are_equal(al, ar);

                ret->v.xval = (f == B_EQ)? eq : !eq;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *augment_array_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL && starting(p)) {
        GretlType lt = gretl_array_get_content_type(l->v.a);
        GretlType rt = gretl_type_from_gen_type(r->t);

        if (rt == lt) {
            ret->v.a = gretl_array_copy(l->v.a, &p->err);
            if (!p->err) {
                p->err = gretl_array_append_object(ret->v.a, r->v.ptr, 1);
            }
        } else {
            p->err = E_TYPES;
        }
    }

    return ret;
}

static NODE *subtract_from_array_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL && starting(p)) {
        if (gretl_array_get_type(l->v.a) == GRETL_TYPE_STRINGS &&
            r->t == STR) {
            ret->v.a = gretl_array_copy(l->v.a, &p->err);
            if (!p->err) {
                p->err = gretl_array_drop_string(ret->v.a, r->v.str);
            }
        }
    } else {
        p->err = E_TYPES;
    }

    return ret;
}

/* in case we switched the LHS and RHS in a boolean comparison */

static int reversed_comp (int f)
{
    if (f == B_GT) {
        return B_LT;
    } else if (f == B_LT) {
        return B_GT;
    } else if (f == B_GTE) {
        return B_LTE;
    } else if (f == B_LTE) {
        return B_GTE;
    } else {
        return f;
    }
}

/* Boolean test of all vars in list against a scalar or series, for
   each observation in the sample, hence generating a series.
   The list will always be on the left-hand node; the 'reversed'
   flag is set if the list was originally on the right.
*/

static NODE *list_bool_comp (NODE *l, NODE *r, int f, int reversed,
                             parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(l, p);
        double *x = ret->v.xvec;
        double xit, targ = NADBL;
        double *tvec = NULL;
        int i, t;

        if (r->t == NUM) {
            targ = r->v.xval;
        } else {
            tvec = r->v.xvec;
        }

        if (reversed) {
            f = reversed_comp(f);
        }

        if (list != NULL) {
            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                if (tvec != NULL) {
                    targ = tvec[t];
                }
                if (na(targ)) {
                    x[t] = NADBL;
                    continue;
                }
                x[t] = 1.0; /* assume 'true' */
                for (i=1; i<=list[0]; i++) {
                    xit = p->dset->Z[list[i]][t];
                    if (na(xit)) {
                        x[t] = NADBL;
                        break;
                    } else if (f == B_EQ && xit != targ) {
                        x[t] = 0.0;
                    } else if (f == B_NEQ && xit == targ) {
                        x[t] = 0.0;
                    } else if (f == B_LT && xit >= targ) {
                        x[t] = 0.0;
                    } else if (f == B_GT && xit <= targ) {
                        x[t] = 0.0;
                    } else if (f == B_LTE && xit > targ) {
                        x[t] = 0.0;
                    } else if (f == B_GTE && xit < targ) {
                        x[t] = 0.0;
                    }
                }
            }
            free(list);
        }
    }

    return ret;
}

/* Test for whether or not two lists are identical.  Note that
   using gretl_list_cmp() the order of the members matters.
   Perhaps the order shouldn't matter?
*/

static NODE *list_list_comp (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        int *llist = node_get_list(l, p);
        int *rlist = node_get_list(r, p);

        if (llist != NULL && rlist != NULL) {
            int d = gretl_list_cmp(llist, rlist);

            if (f == B_NEQ) {
                ret->v.xval = d;
            } else if (f == B_EQ) {
                ret->v.xval = !d;
            } else {
                p->err = E_TYPES;
            }
        }
        free(llist);
        free(rlist);
    }

    return ret;
}

/* argument is list; value returned is series */

static NODE *list_to_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int partial_ok = 0;
        int *list = NULL;

        if (f == F_MIN || f == F_MAX) {
            /* Accept partial observations? This is debatable */
            partial_ok = 1;
        }
        list = node_get_list(n, p);
        if (list != NULL) {
            p->err = cross_sectional_stat(ret->v.xvec, list,
                                          p->dset, f, partial_ok);
            free(list);
        }
    }

    return ret;
}

/* arguments are series on left, list on right: we add all members
   of list to series, or subtract all members */

static NODE *series_list_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(r, p);

        if (list != NULL) {
            double xt, xi;
            int i, t;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                xt = l->v.xvec[t];
                if (!na(xt)) {
                    for (i=1; i<=list[0]; i++) {
                        xi = p->dset->Z[list[i]][t];
                        if (na(xi)) {
                            xt = NADBL;
                            break;
                        } else if (f == B_ADD) {
                            xt += xi;
                        } else {
                            xt -= xi;
                        }
                    }
                }
                ret->v.xvec[t] = xt;
            }
            free(list);
        }
    }

    return ret;
}

static int node_get_midas_method (NODE *n, parser *p)
{
    int ret = -1;

    if (n->t == STR) {
        if (!strcmp(n->v.str, "umidas")) {
            ret = MIDAS_U;
        } else if (!strcmp(n->v.str, "nealmon")) {
            ret = MIDAS_NEALMON;
        } else if (!strcmp(n->v.str, "beta0")) {
            ret = MIDAS_BETA0;
        } else if (!strcmp(n->v.str, "betan")) {
            ret = MIDAS_BETAN;
        } else if (!strcmp(n->v.str, "almonp")) {
            ret = MIDAS_ALMONP;
        } else if (!strcmp(n->v.str, "beta1")) {
            ret = MIDAS_BETA1;
        }
    } else {
        ret = node_get_int(n, p);
    }

    if (ret < MIDAS_U || ret >= MIDAS_MAX) {
        p->err = E_INVARG;
    }

    return ret;
}

static NODE *lincomb_func (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(l, p);
        const gretl_matrix *b = node_get_real_matrix(m, p, 0, 2);
        int k = 0;

        if (!p->err && (list == NULL || gretl_is_null_matrix(b))) {
            p->err = E_DATA;
        }

        if (!p->err && f == F_MLINCOMB) {
            k = node_get_midas_method(r, p);
        }

        if (!p->err) {
            if (f == F_MLINCOMB) {
                p->err = midas_linear_combo(ret->v.xvec, list, b, k, p->dset);
            } else {
                p->err = list_linear_combo(ret->v.xvec, list, b, p->dset);
            }
        }

        free(list);
    }

    return ret;
}

static NODE *list_list_series_func (NODE *l1, NODE *l2, int f,
                                    NODE *o, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int partial_ok = node_get_bool(o, p, 0);
        int *llist = node_get_list(l1, p);
        int *rlist = node_get_list(l2, p);

        if (!p->err) {
            p->err = x_sectional_weighted_stat(ret->v.xvec, llist, rlist,
                                               p->dset, f, partial_ok);
        }
        free(llist);
        free(rlist);
    }

    return ret;
}

/* check for missing obs in a list of variables */

static NODE *list_ok_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        int *list = n->v.ivec;
        int i, vi, t;
        double x;

        if (list[0] == 0) {
            return ret;
        }

        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            x = (f == F_DATAOK)? 1 : 0;
            for (i=1; i<=list[0]; i++) {
                vi = list[i];
                if (na(p->dset->Z[vi][t])) {
                    x = (f == F_DATAOK)? 0 : 1;
                    break;
                }
            }
            ret->v.xvec[t] = x;
        }
    }

    return ret;
}

/* functions taking (up to) two scalars as arguments and
   returning a series result */

static NODE *
series_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        double x, y;

        x = null_node(l) ? NADBL : node_get_scalar(l, p);
        y = null_node(r) ? NADBL : node_get_scalar(r, p);

        switch (f) {
        case F_RUNIFORM:
            p->err = gretl_rand_uniform_minmax(ret->v.xvec,
                                               p->dset->t1,
                                               p->dset->t2,
                                               x, y);
            break;
        case F_RNORMAL:
            p->err = gretl_rand_normal_full(ret->v.xvec,
                                            p->dset->t1,
                                            p->dset->t2,
                                            x, y);
            break;
        default:
            break;
        }
    }

    return ret;
}

static gretl_matrix *fc_matrix_from_list (NODE *n, int n1,
                                          parser *p)
{
    gretl_matrix *ret = NULL;

    if (sample_size(p->dset) != n1) {
        p->err = E_NONCONF;
    } else {
        ret = gretl_matrix_data_subset(n->v.ivec,
                                       p->dset,
                                       p->dset->t1,
                                       p->dset->t2,
                                       M_MISSING_OK,
                                       &p->err);
    }

    return ret;
}

static NODE *fcstats_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (starting(p)) {
        gretl_matrix *Fmat = NULL;
        const double *x = NULL, *y = NULL;
        gretlopt Fopt = OPT_D;
        int U2, free_Fmat = 0;
        int n = 0, n2 = 0;

        if (l->t == SERIES || m->t == SERIES || m->t == LIST) {
            if (dataset_is_time_series(p->dset)) {
                Fopt |= OPT_T;
            }
        }
        U2 = node_get_bool(r, p, (Fopt & OPT_T) ? 1 : 0);
        if (U2) {
            Fopt |= OPT_T;
        } else {
            Fopt &= ~OPT_T;
        }

        if (!p->err) {
            if (l->t == SERIES) {
                n = sample_size(p->dset);
                x = l->v.xvec + p->dset->t1;
            } else {
                n = gretl_vector_get_length(l->v.m);
                if (n == 0) {
                    p->err = E_TYPES;
                } else {
                    x = l->v.m->val;
                }
            }
        }

        if (!p->err) {
            if (m->t == SERIES) {
                n2 = sample_size(p->dset);
                if (n2 != n) {
                    p->err = E_NONCONF;
                } else {
                    y = m->v.xvec + p->dset->t1;
                }
            } else if (m->t == MAT) {
                n2 = gretl_vector_get_length(m->v.m);
                if (n2 != n) {
                    if (m->v.m->rows == n) {
                        Fmat = m->v.m;
                    } else {
                        p->err = E_NONCONF;
                    }
                } else {
                    y = m->v.m->val;
                }
            } else {
                Fmat = fc_matrix_from_list(m, n, p);
                if (!p->err) {
                    free_Fmat = 1;
                }
            }
        }

        if (!p->err) {
            if (Fmat != NULL) {
                ret->v.m = matrix_fc_stats(x, Fmat, Fopt, &p->err);
            } else {
                ret->v.m = forecast_stats(x, y, 0, n-1, NULL, Fopt, &p->err);
            }
            if (free_Fmat) {
                gretl_matrix_free(Fmat);
            }
        }
    }

    return ret;
}

/* Functions taking two series as arguments and returning a scalar
   or matrix result. We also accept as arguments two matrices if
   they are vectors of the same length. In the case of F_NAALEN
   and F_KMEIER we can accept input with no @r node argument
   (meaning no censoring).
*/

static NODE *series_2_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const double *x = NULL, *y = NULL;
        int n = 0, n2 = 0;

        if (!p->err) {
            if (l->t == SERIES) {
                n = sample_size(p->dset);
                x = l->v.xvec + p->dset->t1;
            } else if (l->t == NUM) {
                n = 1;
                x = &l->v.xval;
            } else if (l->t == MAT) {
                n = gretl_vector_get_length(l->v.m);
                if (n == 0) {
                    p->err = E_TYPES;
                } else {
                    x = l->v.m->val;
                }
            } else {
                p->err = E_INVARG;
            }
        }

        if (!p->err) {
            if (null_node(r)) {
                ; /* OK for duration funcs only */
            } else if (r->t == SERIES) {
                /* series on right */
                n2 = sample_size(p->dset);
                if (n2 != n) {
                    p->err = E_NONCONF;
                } else {
                    y = r->v.xvec + p->dset->t1;
                }
            } else if (r->t == NUM) {
                n2 = 1;
                if (n2 != n) {
                    p->err = E_NONCONF;
                } else {
                    y = &r->v.xval;
                }
            } else if (r->t == MAT) {
                /* matrix on right */
                n2 = gretl_vector_get_length(r->v.m);
                if (n2 != n) {
                    p->err = E_NONCONF;
                } else {
                    y = r->v.m->val;
                }
            } else {
                p->err = e_types(r);
            }
        }

        if (p->err) {
            return ret;
        } else if (f == F_NAALEN || f == F_KMEIER) {
            ret = aux_matrix_node(p);
        } else {
            ret = aux_scalar_node(p);
        }
        if (ret == NULL) {
            return NULL;
        }

        switch (f) {
        case F_COR:
            ret->v.xval = gretl_corr(0, n-1, x, y, NULL);
            break;
        case F_COV:
            ret->v.xval = gretl_covar(0, n-1, x, y, NULL);
            break;
        case F_NAALEN:
            ret->v.m = duration_func(x, y, 0, n-1, OPT_NONE, &p->err);
            break;
        case F_KMEIER:
            ret->v.m = duration_func(x, y, 0, n-1, OPT_K, &p->err);
            break;
        case F_CORRESP:
            ret->v.xval = correspondence(x, y, n, &p->err);
            break;
        default:
            break;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static gretlopt get_npcorr_option (NODE *n, parser *p)
{
    gretlopt opt = OPT_NONE;

    if (null_node(n)) {
        ; /* OK */
    } else {
        /* screened already: must be string */
        const char *s = n->v.str;

        if (!strcmp(s, "kendall")) {
            opt = OPT_K;
        } else if (!strcmp(s, "spearman")) {
            opt = OPT_S;
        } else {
            p->err = E_INVARG;
        }
    }

    return opt;
}

static NODE *npcorr_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = NULL, *y = NULL;
        gretlopt opt = OPT_NONE;
        int n1 = 0, n2 = 0;

        if (l->t == SERIES) {
            x = l->v.xvec + p->dset->t1;
            n1 = sample_size(p->dset);
        } else {
            n1 = gretl_vector_get_length(l->v.m);
            if (n1 == 0) {
                p->err = E_INVARG;
            } else {
                x = l->v.m->val;
            }
        }

        if (!p->err && m->t == SERIES) {
            y = m->v.xvec + p->dset->t1;
            n2 = sample_size(p->dset);
        } else if (!p->err) {
            n2 = gretl_vector_get_length(m->v.m);
            if (n2 == 0) {
                p->err = E_INVARG;
            } else {
                y = m->v.m->val;
            }
        }

        if (!p->err && n1 != n2) {
            p->err = E_NONCONF;
        } else if (!p->err) {
            opt = get_npcorr_option(r, p);
        }

        if (!p->err) {
            if (opt & OPT_S) {
                ret->v.m = spearman_rho_func(x, y, n1, &p->err);
            } else {
                ret->v.m = kendall_tau_func(x, y, n1, &p->err);
            }
        }
    }

    return ret;
}

/* takes two series or two vectors as arguments */

static NODE *mxtab_func (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = NULL;
        const double *y = NULL;
        int nx = 0;
        int ny = 0;

        if (l->t == MAT && r->t == MAT) {
            nx = gretl_vector_get_length(l->v.m);
            ny = gretl_vector_get_length(r->v.m);
            x = l->v.m->val;
            y = r->v.m->val;
        } else if (l->t == SERIES && r->t == SERIES) {
            nx = ny = sample_size(p->dset);
            x = l->v.xvec + p->dset->t1;
            y = r->v.xvec + p->dset->t1;
        } else {
            p->err = E_TYPES;
        }
        if (!p->err && (nx == 0 || ny != nx)) {
            p->err = E_NONCONF;
        }
        if (!p->err) {
            ret->v.m = gretl_matrix_xtab(x, y, nx, &p->err);
        }
    }

    return ret;
}

/* remedial code to support <bundle>.<key> as passed to exists() */

static GretlType remedy_old_exists (const char *s)
{
    GretlType ret = GRETL_TYPE_NONE;
    char **S;
    int ns = 0;

    S = gretl_string_split(s, &ns, ".");

    if (ns == 2) {
        gretl_bundle *b = get_bundle_by_name(S[0]);

        if (b != NULL) {
            gretl_bundle_get_data(b, S[1], &ret, NULL, NULL);
        }
    }

    strings_array_free(S, ns);

    return ret;
}

static NODE *object_status (NODE *n, NODE *func, parser *p)
{
    NODE *ret = aux_scalar_node(p);
    int f = func->t;

    if (ret != NULL && starting(p)) {
        const char *s = n->v.str;

        ret->v.xval = NADBL;

        if (f == F_EXISTS) {
            GretlType type = user_var_get_type_by_name(s);

            if (type == 0 && gretl_is_series(s, p->dset)) {
                type = GRETL_TYPE_SERIES;
            } else if (type == 0 && strchr(s, '.')) {
                type = remedy_old_exists(s);
            }
            ret->v.xval = gretl_type_get_order(type);
        } else if (f == F_ISDISCR) {
            int v = current_series_index(p->dset, s);

            if (v >= 0) {
                ret->v.xval = series_is_discrete(p->dset, v);
            }
        } else if (f == F_OBSNUM) {
            int t = get_observation_number(s, p->dset);

            if (t > 0) {
                ret->v.xval = t;
            }
        } else if (f == F_STRLEN) {
            /* ret->v.xval = strlen(s); */
            ret->v.xval = g_utf8_strlen(s, -1);
        } else if (f == F_NLINES) {
            ret->v.xval = count_lines(s);
        } else if (f == F_REMOVE) {
            gretl_maybe_switch_dir(s);
            ret->v.xval = gretl_remove(s);
        }
    }

    return ret;
}

static NODE *multi_str_node (NODE *l, int f, parser *p)
{
    NODE *ret = NULL;

    if (l->t == SERIES) {
        if (!strvals_node(l, p)) {
            p->err = E_TYPES;
        } else {
            ret = aux_series_node(p);
        }
    } else if (l->t == ARRAY) {
        if (gretl_array_get_type(l->v.a) != GRETL_TYPE_STRINGS) {
            p->err = e_types(l);
        } else {
            ret = aux_matrix_node(p);
        }
    } else {
        p->err = e_types(l);
    }

    if (!p->err && l->t == SERIES) {
        series_table *st;
        const char *s;
        int t;

        st = series_get_string_table(p->dset, l->vnum);
        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            s = series_table_get_string(st, l->v.xvec[t]);
            ret->v.xvec[t] = (s == NULL)? NADBL : g_utf8_strlen(s, -1);
        }
    } else if (!p->err) {
        gretl_matrix *m = NULL;
        char **S;
        int i, ns = 0;

        S = gretl_array_get_strings(l->v.a, &ns);
        m = (ns == 0)? gretl_null_matrix_new() :
            gretl_matrix_alloc(ns, 1);
        if (m == NULL) {
            p->err = E_ALLOC;
        } else {
            for (i=0; i<ns; i++) {
                m->val[i] = g_utf8_strlen(S[i], -1);
            }
            ret->v.m = m;
        }
    }

    return ret;
}

static NODE *generic_typeof_node (NODE *n, NODE *func, parser *p)
{
    NODE *ret = NULL;

    if (func->t == F_TYPENAME) {
        ret = aux_string_node(p);
    } else {
        ret = aux_scalar_node(p);
    }

    if (ret != NULL) {
        GretlType t = gretl_type_from_gen_type(n->t);

        if (func->t == F_TYPENAME) {
            if (t == GRETL_TYPE_ARRAY) {
                t = gretl_array_get_type(n->v.a);
            }
            ret->v.str = gretl_strdup(gretl_type_get_name(t));
        } else if (func->t == F_ISCMPLX) {
            if (t == GRETL_TYPE_MATRIX) {
                ret->v.xval = matrix_is_complex(n->v.m);
            } else {
                ret->v.xval = NADBL;
            }
        } else if (alias_reversed(func)) {
            /* handle the "isnull" alias for exists() */
            ret->v.xval = (t == 0);
        } else {
            /* func->t == F_TYPEOF */
            ret->v.xval = gretl_type_get_order(t);
        }
    }

    return ret;
}

/* return scalar node holding the number of elements in
   the object associated with node @n
*/

static NODE *n_elements_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        ret->v.xval = object_get_size(n, 0, p);
    }

    return ret;
}

static int look_up_vname (const char *s, const DATASET *dset)
{
    int i;

    for (i=0; i<dset->v; i++) {
        if (!strcmp(s, dset->varname[i])) {
            return i;
        }
    }

    return -1;
}

/* return scalar node holding the position of the series
   associated with node @r in the list associated with node
   @l, or zero if the series is not present in the list
*/

static NODE *in_list_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (p->err == 0 && (p->dset == NULL || p->dset->v == 0)) {
        p->err = E_NODATA;
    }

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(l, p);

        if (list != NULL) {
            int k = -1;

            if (useries_node(r)) {
                k = r->vnum;
            } else if (r->t == NUM) {
                if (r->v.xval >= 0 && r->v.xval < p->dset->v) {
                    k = (int) r->v.xval;
                }
            } else if (r->t == STR) {
                k = look_up_vname(r->v.str, p->dset);
                if (k < 0) {
                    ret->v.xval = 0;
                }
            } else {
                node_type_error(F_INLIST, 2, SERIES, r, p);
            }
            if (k >= 0) {
                ret->v.xval = in_gretl_list(list, k);
            }
            free(list);
        }
    }

    return ret;
}

static NODE *list_info_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);
    int k = 0;

    if (p->err == 0 && (p->dset == NULL || p->dset->v == 0)) {
        p->err = E_NODATA;
    }

    if (!p->err && !null_node(r)) {
        k = node_get_int(r, p);
    }

    if (!p->err) {
        const int *list = l->v.ivec;
        gretlopt opt = OPT_NONE;

        if (k & 1) {
            opt |= OPT_C;
        }
        if (k & 2) {
            opt |= OPT_B;
        }
        ret->v.m = list_info_matrix(list, p->dset, opt, &p->err);
    }

    return ret;
}

static NODE *argname_from_uvar (NODE *n, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        const char *vname = NULL;

        if (!null_or_string(r)) {
            /* if @r is present it must hold a string */
            p->err = e_types(r);
            return ret;
        }

        if (n->t == SERIES) {
            vname = p->dset->varname[n->vnum];
        } else {
            vname = n->vname;
        }

        if (vname == NULL) {
            p->err = E_DATA;
        } else {
            ret->v.str = gretl_func_get_arg_name(vname, &p->err);
        }

        if (ret->v.str == NULL || ret->v.str[0] == '\0') {
            if (!null_node(r)) {
                ret->v.str = gretl_strdup(r->v.str);
            }
        }
    }

    return ret;
}

static NODE *varnum_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        if (n->t == STR) {
            int v = current_series_index(p->dset, n->v.str);

            ret->v.xval = (v >= 0)? v : NADBL;
        } else {
            p->err = E_DATA;
        }
    }

    return ret;
}

static NODE *int_to_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *v = NULL;
        int i = 0;

        if (scalar_node(n)) {
            i = node_get_int(n, p);
        } else if (n->t == SERIES && f == F_VARNAME) {
            i = n->vnum;
        } else if (n->t == MAT && f == F_OBSLABEL) {
            v = n->v.m;
        } else {
            node_type_error(f, 0, NUM, n, p);
        }

        if (!p->err && v != NULL) {
            ret = aux_array_node(p);
        } else if (!p->err) {
            ret = aux_string_node(p);
        }

        if (!p->err) {
            if (f == F_OBSLABEL && v != NULL) {
                ret->v.a = retrieve_date_strings(v, p->dset, &p->err);
            } else if (f == F_OBSLABEL) {
                ret->v.str = retrieve_date_string(i, p->dset, &p->err);
            } else if (f == F_VARNAME) {
                if (i >= 0 && i < p->dset->v) {
                    ret->v.str = gretl_strdup(p->dset->varname[i]);
                } else {
                    p->err = E_INVARG;
                }
            } else if (f == F_RANDSTR) {
                ret->v.str = gretl_rand_hex_string(i, &p->err);
            } else {
                p->err = E_DATA;
            }
        }

        if (!p->err && v == NULL && ret->v.str == NULL) {
            p->err = E_ALLOC;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *list_to_string_func (NODE *n, int f, parser *p)
{
    NODE *ret;

    if (f == F_VARNAMES) {
        ret = aux_array_node(p);
    } else {
        ret = aux_string_node(p);
    }

    if (ret != NULL && starting(p)) {
        int *list = node_get_list(n, p);

        if (p->err) {
            return ret;
        }

        if (f == F_VARNAME) {
            ret->v.str = gretl_list_get_names(list, p->dset,
                                              &p->err);
        } else if (f == F_VARNAMES) {
            char **S = gretl_list_get_names_array(list, p->dset,
                                                  &p->err);
            int ns = list[0];

            if (!p->err) {
                ret->v.a = gretl_array_from_strings(S, ns, 0, &p->err);
            }
        } else {
            p->err = E_DATA;
        }

        free(list);
    }

    return ret;
}

/* handles both getenv (string value of variable) and
   ngetenv (numerical value of variable)
*/

static NODE *do_getenv (NODE *l, int f, parser *p)
{
    NODE *ret = (f == F_GETENV)? aux_string_node(p) :
        aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        int defined = 0;
        char *estr;

        estr = gretl_getenv(l->v.str, &defined, &p->err);

        if (f == F_GETENV) {
            ret->v.str = estr;
        } else {
            /* ngetenv */
            if (defined) {
                char *test = NULL;
                double x;

                errno = 0;
                x = strtod(estr, &test);
                if (*test == '\0' && errno == 0) {
                    ret->v.xval = x;
                }
            }
            free(estr);
        }
    }

    return ret;
}

static void write_mpi_errmsg (const char *funcname, const char *s)
{
#ifdef HAVE_MPI
    gchar *tmp = gretl_make_dotpath("mpi.fail");
    FILE *fp = gretl_fopen(tmp, "wb");

    if (fp != NULL) {
        if (funcname != NULL) {
            fprintf(fp, _("Error message from %s():\n %s"),
                    funcname, s);
            fputc('\n', fp);
        } else {
            fprintf(fp, "Error message from gretlmpi: %s\n", s);
        }
        fclose(fp);
    }
    g_free(tmp);
#else
    return;
#endif
}

static NODE *do_errorif (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);
    int fd = gretl_function_depth();

    if (fd == 0 && !gretl_mpi_initialized()) {
        gretl_errmsg_sprintf("'%s': can only be used within a function",
                             "errorif");
        p->err = E_DATA;
    } else {
        int cond = node_get_bool(l, p, 1);

        if (cond) {
            const char *funcname = NULL;

            if (fd > 0) {
                current_function_info(&funcname, NULL);
            }
            if (gretl_mpi_initialized()) {
                if (mpi_rank == 0) {
                    write_mpi_errmsg(funcname, r->v.str);
                    p->err = E_FUNCERR;
                }
            } else {
                gretl_errmsg_sprintf(_("Error message from %s():\n %s"),
                                     funcname, r->v.str);
                p->err = E_FUNCERR;
            }
        }
    }

    if (ret != NULL) {
        ret->v.xval = 1;
    }

    return ret;
}

static NODE *do_assert (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);
    int assert_val = libset_get_int(GRETL_ASSERT);

    if (ret == NULL) {
        p->err = E_ALLOC;
        return NULL;
    }

    if (l->v.xval != 0 && !na(l->v.xval)) {
        /* non-zero, non-missing: success */
        ret->v.xval = 1;
    } else if (assert_val == 0) {
        /* flag but ignore failure */
        ret->v.xval = 0;
    } else if (assert_val == 1) {
        /* warn about failure */
        pprintf(p->prn, _("Warning: assertion '%s' failed"), r->v.str);
        pputc(p->prn, '\n');
        ret->v.xval = 0;
    } else if (assert_val == 2) {
        /* complain and throw an error on failure */
        p->err = 1;
        gretl_errmsg_sprintf(_("Assertion '%s' failed"), r->v.str);
        ret->v.xval = l->v.xval;
    } else {
        /* abort execution altogether */
        fprintf(stderr, "Assertion '%s' failed\n", r->v.str);
        abort();
    }

    return ret;
}

static NODE *contains_node (NODE *val, NODE *set, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const gretl_matrix *m = node_get_matrix(set, p, 0, 2);

        if (val->t == NUM) {
            ret = aux_scalar_node(p);
        } else if (val->t == SERIES) {
            ret = aux_series_node(p);
        } else if (val->t == MAT) {
            ret = aux_matrix_node(p);
            if (!p->err) {
                int r = val->v.m->rows;
                int c = val->v.m->cols;

                ret->v.m = gretl_zero_matrix_new(r, c);
            }
        }
        if (!p->err) {
            int i, n = m->rows * m->cols;
            double x;

            if (val->t == NUM) {
                x = val->v.xval;
                if (na(x)) {
                    ret->v.xval = NADBL;
                } else {
                    ret->v.xval = 0;
                    for (i=0; i<n; i++) {
                        if (m->val[i] == x) {
                            ret->v.xval = 1;
                            break;
                        }
                    }
                }
            } else if (val->t == SERIES) {
                int t;

                for (t=p->dset->t1; t<=p->dset->t2; t++) {
                    x = val->v.xvec[t];
                    if (na(x)) {
                        ret->v.xvec[t] = NADBL;
                    } else {
                        ret->v.xvec[t] = 0;
                        for (i=0; i<n; i++) {
                            if (m->val[i] == x) {
                                ret->v.xvec[t] = 1;
                                break;
                            }
                        }
                    }
                }
            } else if (val->t == MAT) {
                gretl_matrix *v = val->v.m;
                int k, nv = v->rows * v->cols;

                for (k=0; k<nv; k++) {
                    x = v->val[k];
                    if (na(x)) {
                        ret->v.m->val[k] = NADBL;
                    } else {
                        ret->v.m->val[k] = 0;
                        for (i=0; i<n; i++) {
                            if (m->val[i] == x) {
                                ret->v.m->val[k] = 1;
                                break;
                            }
                        }
                    }
                }
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *single_string_func (NODE *n, NODE *x, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        const char *s = n->v.str;

        if (f == F_ARGNAME) {
            char *deflt = NULL;

            if (!null_node(x)) {
                if (x->t == STR) {
                    deflt = x->v.str;
                } else {
                    p->err = e_types(x);
                }
            }
            s = (n->vname != NULL)? n->vname : s;
            ret->v.str = gretl_func_get_arg_name(s, &p->err);
            if (!p->err && ret->v.str[0] == '\0' && deflt != NULL) {
                ret->v.str = gretl_strdup(deflt);
            }
        } else if (f == F_BACKTICK) {
            ret->v.str = gretl_backtick(s, &p->err);
        } else if (f == F_STRSTRIP) {
            ret->v.str = gretl_strstrip_copy(s, &p->err);
        } else if (f == F_FIXNAME) {
            int uscore = 0;

            if (!null_node(x)) {
                uscore = node_get_bool(x, p, 0);
            }
            ret->v.str = calloc(VNAMELEN, 1);
            gretl_normalize_varname(ret->v.str, s, uscore, 0);
        } else {
            p->err = E_DATA;
        }
    }

    return ret;
}

static NODE *country_code_node (NODE *n, NODE *r, parser *p)
{
    NODE *ret = NULL;
    char *src = NULL;
    char *tmp = NULL;

    if (n->t == STR) {
        src = n->v.str;
        ret = aux_string_node(p);
    } else if (n->t == ARRAY || n->t == SERIES || n->t == MAT) {
        ret = aux_array_node(p);
    } else if (n->t == NUM) {
        int k = node_get_int(n, p);

        if (!p->err) {
            src = tmp = g_strdup_printf("%03d", k);
            ret = aux_string_node(p);
        }
    } else {
        p->err = E_INVARG;
    }

    if (!p->err) {
        int output = 0; /* default to automatic */

        if (!null_node(r)) {
            output = node_get_int(r, p);
        }
        if (!p->err && src != NULL) {
            char *(*cfunc) (const char *, int, PRN *, int *);

            cfunc = get_plugin_function("iso_country");
            if (cfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                ret->v.str = cfunc(src, output, p->prn, &p->err);
            }
        } else if (!p->err && n->t == ARRAY) {
            gretl_array *(*cfunc) (gretl_array *, int, PRN *, int *);

            cfunc = get_plugin_function("iso_country_array");
            if (cfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                ret->v.a = cfunc(n->v.a, output, p->prn, &p->err);
            }
        } else if (!p->err) {
            gretl_array *(*cfunc) (const double *, int, int, PRN *, int *);

            cfunc = get_plugin_function("iso_country_series");
            if (cfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                const double *x;
                int nx;

                if (n->t == SERIES) {
                    nx = sample_size(p->dset);
                    x = n->v.xvec + p->dset->t1;
                } else {
                    /* matrix */
                    nx = gretl_vector_get_length(n->v.m);
                    if (nx == 0) {
                        p->err = E_INVARG;
                    } else {
                        x = n->v.m->val;
                    }
                }
                if (!p->err) {
                    ret->v.a = cfunc(x, nx, output, p->prn, &p->err);
                }
            }
        }
    }

    g_free(tmp);

    return ret;
}

static NODE *readfile_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        const char *codeset = NULL;

        if (!null_node(r)) {
            if (r->t == STR) {
                codeset = r->v.str;
            } else {
                node_type_error(F_READFILE, 2, STR, r, p);
            }
        }
        if (!p->err) {
            ret->v.str = retrieve_file_content(l->v.str, codeset, &p->err);
        }
    }

    return ret;
}

#if 0
static void strstr_escape (char *s)
{
    int i, n = strlen(s);

    for (i=0; i<n; i++) {
        if (s[i] == '\\' && (i == 0 || s[i-1] != '\\')) {
            if (s[i+1] == 'n') {
                s[i] = '\n';
                shift_string_left(s + i + 1, 1);
                i++;
            } else if (s[i+1] == 't') {
                s[i] = '\t';
                shift_string_left(s + i + 1, 1);
                i++;
            }
        }
    }
}
#endif

static const char *strstr_insensitive (const gchar *sd1,
                                       const gchar *sd2,
                                       const char *sl)
{
    const char *ret = strstr(sd1, sd2);

    if (ret != NULL) {
        /* preserve case of @sl? */
        const char *chk = sl + (ret - sd1);

        if (g_utf8_validate(chk, -1, NULL)) {
            ret = chk;
        }
    }

    return ret;
}

static NODE *two_string_func (NODE *l, NODE *r, NODE *x,
                              int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *sl = l->v.str;
        const char *sr = NULL;

        if (f == F_JSONGETB) {
            ; /* checks done below */
        } else if (f == F_XMLGET && r->t == ARRAY) {
            if (gretl_array_get_type(r->v.a) != GRETL_TYPE_STRINGS) {
                p->err = E_TYPES;
            }
        } else if (r->t != STR) {
            p->err = e_types(r);
        }
        if (!p->err) {
            ret = (f == F_INSTRING)? aux_scalar_node(p) :
                (f == F_JSONGETB)? aux_bundle_node(p) :
                aux_string_node(p);
        }
        if (p->err) {
            return NULL;
        }

        if (r != NULL && r->t == STR) {
            sr = r->v.str;
        }

        if (f == F_STRSTR || f == F_INSTRING) {
            int ignore_case = node_get_bool(x, p, 0);
            const char *sret = NULL;

            if (ignore_case) {
                gchar *sd1 = g_utf8_strdown(sl, -1);
                gchar *sd2 = g_utf8_strdown(sr, -1);

                sret = strstr_insensitive(sd1, sd2, sl);
                g_free(sd1); g_free(sd2);
            } else {
                sret = strstr(sl, sr);
            }
            if (f == F_INSTRING) {
                ret->v.xval = sret != NULL;
            } else if (sret != NULL) {
                ret->v.str = gretl_strdup(sret);
            } else {
                ret->v.str = gretl_strdup("");
            }
        } else if (f == B_HCAT) {
            int n1 = strlen(l->v.str);
            int n2 = strlen(r->v.str);

            ret->v.str = malloc(n1 + n2 + 1);
            if (ret->v.str != NULL) {
                *ret->v.str = '\0';
                strcat(ret->v.str, l->v.str);
                strcat(ret->v.str, r->v.str);
            }
        } else if (f == F_JSONGET) {
            char *(*jfunc) (const char *, const char *,
                            int *, int *) = NULL;
            user_var *uv = NULL;

            if (!null_node(x)) {
                uv = ptr_node_get_uvar(x, NUM, p);
            }
            if (!p->err) {
                jfunc = get_plugin_function("json_get_string");
                if (jfunc == NULL) {
                    p->err = E_FOPEN;
                }
            }
            if (jfunc != NULL) {
                int nobj = 0;
                int *pnobj = uv == NULL ? NULL : &nobj;

                ret->v.str = jfunc(l->v.str, r->v.str, pnobj, &p->err);
                if (!p->err && uv != NULL) {
                    user_var_set_scalar_value(uv, (double) *pnobj);
                }
            }
        } else if (f == F_JSONGETB) {
            gretl_bundle *(*jfunc) (const char *, const char *, int *);
            const char *path = null_node(r) ? NULL: r->v.str;

            jfunc = get_plugin_function("json_get_bundle");
            if (jfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                ret->v.b = jfunc(l->v.str, path, &p->err);
            }
        } else if (f == F_XMLGET) {
            char *(*xfunc) (const char *, void *, GretlType,
                            int *, int *) = NULL;
            user_var *uv = NULL;

            if (!null_node(x)) {
                uv = ptr_node_get_uvar(x, NUM, p);
            }
            if (!p->err) {
                xfunc = get_plugin_function("xml_get");
                if (xfunc == NULL) {
                    p->err = E_FOPEN;
                }
            }
            if (xfunc != NULL) {
                int nobj = 0;
                int *pnobj = uv == NULL ? NULL : &nobj;

                if (r->t == ARRAY) {
                    ret->v.str = xfunc(l->v.str, r->v.a, GRETL_TYPE_ARRAY,
                                       pnobj, &p->err);
                } else {
                    ret->v.str = xfunc(l->v.str, r->v.str, GRETL_TYPE_STRING,
                                       pnobj, &p->err);
                }
                if (!p->err && uv != NULL) {
                    user_var_set_scalar_value(uv, (double) *pnobj);
                }
            }
        } else {
            p->err = E_DATA;
        }

        if (!p->err && f != F_INSTRING && ret->v.str == NULL) {
            p->err = E_ALLOC;
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *strseq_node (NODE *l, NODE *r, int rstr, parser *p)
{
    NODE *ret = NULL;
    NODE *mn = rstr ? l : r; /* the matrix or array node */
    NODE *sn = rstr ? r : l; /* the string node */
    const char *s = sn->v.str;
    int dim;

    if (mn->t == ARRAY) {
        dim = gretl_array_get_length(mn->v.a);
    } else {
        dim = gretl_vector_get_length(mn->v.m);
    }
    if (dim == 0) {
        p->err = E_INVARG;
        return NULL;
    }

    ret = aux_array_node(p);

    if (ret != NULL) {
        gretl_array *S = gretl_array_new(GRETL_TYPE_STRINGS, dim, &p->err);
        gchar *tmp, *ds;
        int i, d;

        for (i=0; i<dim && !p->err; i++) {
            ds = NULL;
            if (mn->t == ARRAY) {
                d = 0;
                ds = g_strdup((char *) gretl_array_get_data(mn->v.a, i));
            } else {
                d = floor(mn->v.m->val[i]);
                ds = g_strdup_printf("%d", d);
            }
            if (d < 0) {
                p->err = E_INVARG;
            } else {
                if (rstr) {
                    tmp = g_strdup_printf("%s%s", ds, s);
                } else {
                    tmp = g_strdup_printf("%s%s", s, ds);
                }
                gretl_array_set_string(S, i, tmp, 1);
                g_free(tmp);
            }
            g_free(ds);
        }
        ret->v.a = S;
    }

    return ret;
}

static int is_strings_array_node (NODE *n)
{
    return n->t == ARRAY &&
        gretl_array_get_type(n->v.a) == GRETL_TYPE_STRINGS;
}

static NODE *horizontal_concat_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (l->t == STR && r->t == STR) {
        ret = two_string_func(l, r, NULL, B_HCAT, p);
    } else if (ok_matrix_node(l) && ok_matrix_node(r)) {
        ret = matrix_matrix_calc(l, r, B_HCAT, p);
    } else if (l->t == STR && (ok_matrix_node(r) || is_strings_array_node(r))) {
        ret = strseq_node(l, r, 0, p);
    } else if ((ok_matrix_node(l) || is_strings_array_node(l)) && r->t == STR) {
        ret = strseq_node(l, r, 1, p);
    } else if (is_strings_array_node(l) && is_strings_array_node(r)) {
        ret = aux_array_node(p);
        if (ret != NULL) {
            ret->v.a = gretl_arrays_concat(l->v.a, r->v.a, &p->err);
        }
    } else {
        p->err = E_TYPES;
    }

    return ret;
}

#define is_matcher(f) (f == F_STRSUB || f == F_REGSUB)

/* helper function for variant_string_func(): multi-string case */

static char **get_transformed_strings (NODE *l, NODE *m, NODE *r,
                                       int f, parser *p, int *ns)
{
    char **S1 = NULL;
    char **S0;
    char *match = NULL;
    char *repl = NULL;
    int ini = 0;
    int fin = 0;

    if (l->t == ARRAY) {
        S0 = gretl_array_get_strings(l->v.a, ns);
    } else {
        S0 = series_get_string_vals(p->dset, l->vnum, ns, 0);
    }

    if (S0 == NULL) {
        p->err = E_DATA;
    } else {
        S1 = strings_array_new(*ns);
        if (S1 == NULL) {
            p->err = E_ALLOC;
        }
    }

    if (f == F_SUBSTR) {
        ini = node_get_int(m, p);
        fin = node_get_int(r, p);
    } else if (is_matcher(f)) {
        match = m->v.str;
        repl = r->v.str;
    }

    if (!p->err) {
        gchar *tmp = NULL;
        int i;

        /* FIXME compile and re-use regexps? */

        for (i=0; i<*ns && !p->err; i++) {
            if (f == F_SUBSTR) {
                S1[i] = gretl_substring(S0[i], ini, fin, &p->err);
            } else if (f == F_STRSUB) {
                S1[i] = gretl_literal_replace(S0[i], match, repl, &p->err);
            } else if (f == F_REGSUB) {
                tmp = gretl_regexp_replace(S0[i], match, repl, &p->err);
            } else if (f == F_TOLOWER) {
                tmp = g_utf8_strdown(S0[i], -1);
            } else {
                tmp = g_utf8_strup(S0[i], -1);
            }
            if (tmp != NULL) {
                S1[i] = gretl_strdup(tmp);
                g_free(tmp);
            }
        }
    }

    return S1;
}

/* helper function for variant_string_func(): single string case */

static char *get_transformed_string (const char *s,
                                     NODE *m, NODE *r,
                                     int f, parser *p)
{
    char *match = is_matcher(f) ? m->v.str : NULL;
    char *repl =  is_matcher(f) ? r->v.str : NULL;
    char *ret = NULL;
    gchar *tmp = NULL;

    if (f == F_SUBSTR) {
        int ini = node_get_int(m, p);
        int fin = node_get_int(r, p);

        if (!p->err) {
            ret = gretl_substring(s, ini, fin, &p->err);
        }
    } else if (f == F_STRSUB) {
        ret = gretl_literal_replace(s, match, repl, &p->err);
    } else if (f == F_REGSUB) {
        tmp = gretl_regexp_replace(s, match, repl, &p->err);
    } else if (f == F_TOLOWER) {
        tmp = g_utf8_strdown(s, -1);
    } else if (f == F_TOUPPER) {
        tmp = g_utf8_strup(s, -1);
    }
    if (tmp != NULL) {
        ret = gretl_strdup(tmp);
        g_free(tmp);
    }

    return ret;
}

/* A function that can be applied to a single string, an array of strings,
   or a string-valued series. If @m and @r are non-NULL they have already
   been verified as of type STR.
*/

static NODE *variant_string_func (NODE *l, NODE *m, NODE *r,
                                  int f, parser *p)
{
    NODE *ret = NULL;

    if (l->t == STR) {
        ret = aux_string_node(p);
    } else if (l->t == ARRAY) {
        if (gretl_array_get_type(l->v.a) != GRETL_TYPE_STRINGS) {
            p->err = E_TYPES;
        } else {
            ret = aux_array_node(p);
        }
    } else if (l->t == SERIES) {
        p->err = strvals_return_check(l->parent, p);
        if (!p->err) {
            ret = aux_series_node(p);
        }
    } else {
        p->err = E_DATA;
    }

    if (ret == NULL) {
        return NULL;
    }

    if (l->t == STR) {
        /* a single string target */
        ret->v.str = get_transformed_string(l->v.str, m, r, f, p);
    } else {
        /* array or series target */
        char **S = NULL;
        int ns = 0;

        S = get_transformed_strings(l, m, r, f, p, &ns);
        if (l->t == ARRAY) {
            ret->v.a = gretl_array_from_strings(S, ns, 0, &p->err);
        } else {
            double *y = ret->v.xvec + p->dset->t1;
            double *x = l->v.xvec + p->dset->t1;
            int n = sample_size(p->dset);

            p->err = series_from_string_transform(y, x, n, S, ns,
                                                  &p->lh.stab);
        }
    }

    return ret;
}

static char *escape_strsplit_sep (const char *s)
{
    char *ret = calloc(strlen(s) + 1, 1);
    int i = 0;

    while (*s) {
        if (*s == '\\') {
            if (*(s+1) == '\\') {
                ret[i++] = '\\';
                s++;
            } else if (*(s+1) == 'n') {
                ret[i++] = '\n';
                s++;
            } else if (*(s+1) == 'r') {
                ret[i++] = '\r';
                s++;
            } else if (*(s+1) == 't') {
                ret[i++] = '\t';
                s++;
            } else {
                ret[i++] = *s;
            }
        } else {
            ret[i++] = *s;
        }
        s++;
    }

    return ret;
}

static NODE *strsplit_node (int f, NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *s = l->v.str;
        const char *sep0 = NULL;
        char *sep = NULL;
        int k = 0;

        /* We'll accept the two trailing optional arguments,
           string separator and integer index, in either order.
        */
        if (m != NULL) {
            if (m->t == STR) {
                sep0 = m->v.str;
            } else if (m->t != EMPTY) {
                k = node_get_int(m, p);
                if (k < 1) {
                    p->err = E_INVARG;
                }
            }
        }
        if (!p->err && r != NULL) {
            if (r->t == STR) {
                if (sep0 == NULL) {
                    /* OK, didn't get @sep0 yet */
                    sep0 = r->v.str;
                } else {
                    p->err = E_INVARG;
                }
            } else if (r->t != EMPTY) {
                if (k == 0) {
                    /* OK, didn't get @k yet */
                    k = node_get_int(r, p);
                    if (k < 1) {
                        p->err = E_INVARG;
                    }
                } else {
                    p->err = E_INVARG;
                }
            }
        }

        if (!p->err) {
            ret = k > 0 ? aux_string_node(p) : aux_array_node(p);
        }

        if (!p->err) {
            if (sep0 == NULL) {
                /* default: split on white space */
                sep = gretl_strdup(" \t\r\n");
            } else if (strchr(sep0, '\\')) {
                sep = escape_strsplit_sep(sep0);
            } else {
                sep = gretl_strdup(sep0);
            }
        }

        if (!p->err) {
            char **S = NULL;
            int ns = 0;

            S = gretl_string_split(s, &ns, sep);
            if (!p->err) {
                if (k > 0) {
                    ret->v.str = gretl_strdup(k > ns ? "" : S[k-1]);
                    strings_array_free(S, ns);
                } else {
                    ret->v.a = gretl_array_from_strings(S, ns, 0, &p->err);
                }
            }
        }

        free(sep);
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *array_sort_node (NODE *n, int f, parser *p)
{
    NODE *ret = NULL;

    if (gretl_array_get_type(n->v.a) != GRETL_TYPE_STRINGS) {
        p->err = e_types(n);
    } else {
        ret = aux_array_node(p);
        if (!p->err) {
            ret->v.a = gretl_strings_sort(n->v.a, f == F_DSORT, &p->err);
        }
    }
    return ret;
}

static NODE *array_func_node (NODE *l, NODE *r, NODE *b, int f, parser *p)
{
    gretl_array *a = l->v.a;
    GretlType t = gretl_array_get_type(a);
    NODE *ret = NULL;

    if (f == F_INSTRINGS) {
        if (t != GRETL_TYPE_STRINGS || r->t != STR) {
            p->err = E_TYPES;
        } else {
            int simple = node_get_bool(b, p, 0);
            const char *s = r->v.str;

            if (!p->err && simple) {
                ret = aux_scalar_node(p);
                if (!p->err) {
                    ret->v.xval = gretl_strings_array_includes(a, s);
                }
            } else if (!p->err) {
                ret = aux_matrix_node(p);
                if (!p->err) {
                    ret->v.m = gretl_strings_array_pos(a, s, &p->err);
                }
            }
        }
    } else if (f == F_ASORT) {
        if (r->t != STR) {
            p->err = E_TYPES;
        } else {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            ret->v.xval = p->err = gretl_array_qsort(a, r->v.str,
                                                     p->dset, p->prn);
        }
    } else if (f == F_FLATTEN) {
        if (t == GRETL_TYPE_MATRICES) {
            int mode = 0;

            if (!null_node(r)) {
                mode = node_get_int(r, p);
                if (mode < 0 || mode > 2) {
                    p->err = E_INVARG;
                }
            }
            if (!p->err) {
                ret = aux_matrix_node(p);
            }
            if (!p->err) {
                ret->v.m = gretl_matrix_array_flatten(l->v.a, mode, &p->err);
            }
        } else if (t == GRETL_TYPE_STRINGS) {
            const char *sep = NULL;

            if (r->t == STR) {
                sep = r->v.str;
            } else {
                int space = node_get_bool(r, p, 0);

                if (!p->err) {
                    sep = space ? " " : "\n";
                }
            }
            if (!p->err) {
                ret = aux_string_node(p);
            }
            if (!p->err) {
                ret->v.str = gretl_strings_array_flatten(l->v.a, sep, &p->err);
            }
        }
    } else {
        p->err = E_TYPES;
    }

    return ret;
}

static NODE *errmsg_node (NODE *l, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
        const char *src = NULL;

        if (null_node(l)) {
            src = gretl_errmsg_get();
        } else {
            int errval = node_get_int(l, p);

            if (errval < 0 || errval >= E_MAX) {
                p->err = E_DATA;
            } else {
                src = errmsg_get_with_default(errval);
            }
        }

        if (src != NULL) {
            ret->v.str = gretl_strdup(src);
            if (ret->v.str == NULL) {
                p->err = E_ALLOC;
            }
        }
    }

    return ret;
}

static NODE *aux_sized_array_node (GretlType type, int n, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL) {
        ret->v.a = gretl_array_new(type, n, &p->err);
    }

    return ret;
}

/* strftime_node() implements both strftime() and strfday().  In both
   cases the @l node must hold a scalar, vector, or numeric series;
   for strftime() this is interpreted as a time_t, while for strfday
   it is interpreted as a GLib "julian day" or "epoch day".

   isodate() also piggy-backs here when a string-valued series is
   wanted on return.
*/

static NODE *strftime_node (NODE *l, NODE *r, NODE *o, int f,
                            int julian, parser *p)
{
    NODE *ret = NULL;
    const char *fmt = NULL;
    int *vmap = NULL;
    char **S = NULL;
    double offset = NADBL;
    int nv = 0;

    if (l->t == SERIES) {
        p->err = strvals_return_check(l->parent, p);
        if (!p->err) {
            ret = aux_series_node(p);
        }
    } else if (l->t == MAT) {
        nv = gretl_vector_get_length(l->v.m);
        if (nv == 0) {
            p->err = E_INVARG;
        } else if (nv == 1) {
            ret = aux_string_node(p);
        } else {
            ret = aux_sized_array_node(GRETL_TYPE_STRINGS, nv, p);
        }
    } else {
        ret = aux_string_node(p);
    }

    if (ret == NULL) {
        return NULL;
    }

    /* if @r isn't null or empty it should hold a format string */
    if (r != NULL) {
        if (r->t == STR) {
            fmt = r->v.str;
        } else if (r->t != EMPTY) {
            p->err = e_types(r);
        }
    }

    /* if @o isn't null or empty it should hold an offset in seconds */
    if (!null_node(o)) {
        offset = node_get_scalar(o, p);
    }

    if (!p->err && l->t == SERIES) {
        const double *x = l->v.xvec + p->dset->t1;
        int n = sample_size(p->dset);

        vmap = maybe_get_values_map(x, n, &nv, &p->err);
        if (!p->err) {
            S = strings_array_new(nv);
            if (S == NULL) {
                p->err = E_ALLOC;
            }
        }
    }

    if (!p->err) {
        char buf[64] = {0};
        char *dstr;
        int t, t1 = 0, t2 = 0;
        double tx;
        gint64 tt;
        guint32 ed;
        int bytes;
        int i = 0;

        if (fmt == NULL && (f == F_STRFDAY || f == F_ISODATE)) {
            /* default to ISO 8601 */
            fmt = "%Y-%m-%d";
        }
        if (l->t == SERIES) {
            t1 = p->dset->t1;
            t2 = p->dset->t2;
        } else if (l->t == MAT) {
            t2 = nv - 1;
        }

        for (t=t1; t<=t2 && !p->err; t++) {
            tx = (l->t == SERIES)? l->v.xvec[t] : (l->t == MAT)? l->v.m->val[t] :
                node_get_scalar(l, p);
            if (na(tx)) {
                if (l->t == SERIES) {
                    ret->v.xvec[t] = NADBL;
                } else if (l->t != MAT) {
                    ret->v.str = gretl_strdup("");
                }
            } else if (vmap != NULL && S[vmap[i]] != NULL) {
                /* a repeated date in series mode */
                ret->v.xvec[t] = vmap[i] + 1;
                i++;
            } else {
                bytes = 0;
                if (f == F_STRFDAY) {
                    ed = (guint32) floor(tx);
                    bytes = gretl_strfdate(buf, sizeof buf, fmt, ed);
                } else if (f == F_ISODATE) {
                    ed = (guint32) floor(tx);
                    bytes = gretl_alt_strfdate(buf, sizeof buf, julian, ed);
                } else {
                    tt = (gint64) floor(tx);
                    bytes = gretl_strftime(buf, sizeof buf, fmt, tt, offset);
                }
                if (bytes > 0) {
                    dstr = gretl_strdup(g_strchomp(buf));
                    if (l->t == SERIES) {
                        if (vmap != NULL) {
                            ret->v.xvec[t] = vmap[i] + 1;
                            S[vmap[i++]] = dstr;
                        } else {
                            ret->v.xvec[t] = i + 1;
                            S[i++] = dstr;
                        }
                    } else if (l->t == MAT && nv > 1) {
                        p->err = gretl_array_set_string(ret->v.a, t, dstr, 0);
                    } else {
                        ret->v.str = dstr;
                    }
                } else if (l->t == SERIES) {
                    ret->v.xvec[t] = NADBL;
                } else {
                    ret->v.str = gretl_strdup("");
                }
            }
        }
    }

    if (l->t == SERIES) {
        if (p->err && S != NULL) {
            strings_array_free(S, nv);
        } else if (!p->err) {
            prep_strvals_return(S, nv, p);
        }
        free(vmap);
    }

    return ret;
}

static NODE *isodate_node (NODE *l, NODE *r, int f, parser *p)
{
    int julian = (f == F_JULDATE);
    int as_string;
    int n = 0;
    NODE *ret = NULL;

    as_string = node_get_bool(r, p, 0);

    if (l->t == MAT) {
        n = gretl_vector_get_length(l->v.m);
        if (n == 0) {
            p->err = E_NONCONF;
        }
    } else if (!scalar_node(l) && l->t != SERIES) {
        node_type_error(f, 1, NUM, l, p);
    }

    if (p->err) {
        return NULL;
    } else if ((l->t == SERIES || l->t == MAT) && as_string) {
        return strftime_node(l, NULL, NULL, F_ISODATE, julian, p);
    }

    if (l->t == MAT) {
        ret = aux_sized_matrix_node(p, n, 1, 0);
    } else if (l->t == NUM) {
        ret = as_string ? aux_string_node(p) : aux_scalar_node(p);
    } else {
        ret = aux_series_node(p);
    }

    if (ret != NULL && scalar_node(l)) {
        double x = node_get_scalar(l, p);

        if (!as_string && na(x)) {
            ret->v.xval = NADBL;
        } else if (x < 1 || x > UINT_MAX) {
            p->err = E_INVARG;
        } else if (as_string) {
            ret->v.str = ymd_extended_from_epoch_day((guint32) x,
                                                     julian, &p->err);
        } else {
            ret->v.xval = ymd_basic_from_epoch_day((guint32) x,
                                                   julian, &p->err);
        }
    } else if (ret != NULL) {
        double xt, *targ;
        int t1 = l->t == MAT ? 0 : p->dset->t1;
        int t2 = l->t == MAT ? n-1 : p->dset->t2;
        int t;

        for (t=t1; t<=t2; t++) {
            targ = l->t == MAT ? &ret->v.m->val[t] : &ret->v.xvec[t];
            xt = l->t == MAT ? l->v.m->val[t] : l->v.xvec[t];
            if (na(xt)) {
                *targ = NADBL;
            } else if (xt >= 1 && xt <= UINT_MAX) {
                *targ = ymd_basic_from_epoch_day((guint32) xt,
                                                 julian, &p->err);
            } else {
                p->err = E_INVARG;
                break;
            }
        }
    }

    return ret;
}

static int strptime_error (const char *s, double x, int t, int ntype)
{
    if (ntype == STR) {
        gretl_errmsg_sprintf("'%s': argument doesn't match format", s);
    } else if (ntype == NUM) {
        gretl_errmsg_sprintf("Invalid date number %.12g", x);
    } else if (s != NULL) {
        gretl_errmsg_sprintf("'%s': argument doesn't match format at "
                             "obs %d\n", s, t);
    } else if (!na(x)) {
        gretl_errmsg_sprintf("Invalid date number %.12g at obs %d\n",
                             x, t);
    }

    return E_INVARG;
}

/* strptime_node() implements both strptime() and strpday(). In the
   primary use case the @l node contains a string, array of
   strings, or string-valued series. In the secondary case @l
   contains ISO 8601 "basic" dates (scalar, series or vector).
*/

static NODE *strptime_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    const char *fmt = NULL;
    const char *src = NULL;
    int str_series = 0;
    int v = -1;
    int n = 0;

    if (l->t == SERIES) {
        /* either string-valued or ISO basic */
        v = l->vnum;
        if (strvals_node(l, p)) {
            str_series = 1;
        }
        ret = aux_series_node(p);
    } else if (l->t == ARRAY) {
        /* must be strings array */
        if (gretl_array_get_type(l->v.a) != GRETL_TYPE_STRINGS) {
            p->err = E_INVARG;
        } else {
            n = gretl_array_get_length(l->v.a);
            ret = aux_sized_matrix_node(p, n, 1, 0);
        }
    } else if (l->t == MAT) {
        /* must be ISO basic */
        n = gretl_vector_get_length(l->v.m);
        if (n == 0) {
            p->err = E_INVARG;
        } else {
            ret = aux_sized_matrix_node(p, n, 1, 0);
        }
    } else {
        /* single string or ISO basic */
        if (l->t == STR) {
            src = l->v.str;
        }
        ret = aux_scalar_node(p);
    }

    if (ret == NULL) {
        return NULL;
    }

    if (src == NULL && !str_series && l->t != ARRAY) {
        /* without a string arg we won't accept a format string */
        if (r->t != EMPTY) {
            p->err = E_INVARG;
        }
    } else {
        /* if @r isn't empty it should hold a format string */
        if (r->t == STR) {
            fmt = r->v.str;
        } else if (r->t != EMPTY) {
            p->err = e_types(r);
        }
    }

    if (!p->err) {
        double *targ = NULL;
        int t, t1 = 0, t2 = 0;
        char buf[32];
        char *s;

        if (l->t == SERIES) {
            t1 = p->dset->t1;
            t2 = p->dset->t2;
        } else if (l->t == ARRAY) {
            t2 = gretl_array_get_length(l->v.a) - 1;
        } else if (l->t == MAT) {
            t2 = n - 1;
        } else {
            targ = &ret->v.xval;
        }

        for (t=t1; t<=t2; t++) {
            struct tm tm = {0,0,0,1,0,0,0,0,-1};
            int handled = 0;
            double x = NADBL;

            if (l->t == SERIES) {
                if (str_series) {
                    src = series_get_string_for_obs(p->dset, v, t);
                } else {
                    x = l->v.xvec[t];
                }
                targ = &ret->v.xvec[t];
            } else if (l->t == ARRAY) {
                if ((src = gretl_array_get_data(l->v.a, t)) == NULL) {
                    p->err = E_INVARG;
                    break;
                }
                targ = &ret->v.m->val[t];
            } else if (l->t == MAT) {
                x = l->v.m->val[t];
                targ = &ret->v.m->val[t];
            } else if (l->t == NUM) {
                x = l->v.xval;
            }

            if (src == NULL) {
                /* has to be ISO 8601 basic */
                if (na(x)) {
                    *targ = NADBL;
                    continue;
                } else if (x < 0 || x != floor(x) || x > 99991231) {
                    p->err = strptime_error(src, x, t, l->t);
                    break;
                } else {
                    sprintf(buf, "%08d", (int) x);
                    s = strptime(buf, "%Y%m%d", &tm);
                }
            } else {
                if (fmt == NULL) {
                    /* default to ISO 8601 extended */
                    fmt = "%Y-%m-%d";
                }
                if (f == F_STRPDAY) {
                    s = strptime(src, fmt, &tm);
                } else {
                    s = gretl_strptime(src, fmt, targ);
                    handled = 1;
                }
            }

            if (s == NULL) {
                /* strptime() failed */
                p->err = strptime_error(src, x, t, l->t);
                break;
            } else if (f == F_STRPDAY) {
                int y = tm.tm_year + 1900;
                int m = tm.tm_mon + 1;

                *targ = (double) epoch_day_from_ymd(y, m, tm.tm_mday);
                if (*targ == 0) {
                    *targ = NADBL;
                }
            } else if (!handled) {
                /* F_STRPTIME */
#ifdef WIN32
                *targ = win32_mktime(&tm);
#else
                *targ = (double) mktime(&tm);
#endif
            }
        }
    }

    return ret;
}

static NODE *atof_node (NODE *l, parser *p)
{
    NODE *ret = NULL;
    char *endptr = NULL;
    int warn = 0;

    errno = 0;

    gretl_push_c_numeric_locale();

    if (l->t == STR) {
        ret = aux_scalar_node(p);
        if (ret != NULL && starting(p)) {
            if (strchr(l->v.str, ',') != NULL) {
                warn = 1;
            }
            ret->v.xval = strtod(l->v.str, &endptr);
            if (errno || endptr == l->v.str) {
                errno = 0;
                ret->v.xval = NADBL;
            }
        }
    } else if (l->t == SERIES) {
        int v = l->vnum;

        if (!strvals_node(l, p)) {
            p->err = E_TYPES;
        } else {
            ret = aux_series_node(p);
        }
        if (ret != NULL && starting(p)) {
            const char *st;
            int t;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                st = series_get_string_for_obs(p->dset, v, t);
                if (st == NULL) {
                    /* happens if obs @t is missing */
                    ret->v.xvec[t] = NADBL;
                } else {
                    if (!warn && strchr(st, ',') != NULL) {
                        warn = 1;
                    }
                    ret->v.xvec[t] = strtod(st, &endptr);
                    if (errno || endptr == st) {
                        errno = 0;
                        ret->v.xvec[t] = NADBL;
                    }
                }
            }
        }
    }

    if (warn) {
        gretl_warnmsg_set(_("comma(s) found in atof() input"));
    }

    gretl_pop_c_numeric_locale();

    return ret;
}

static void strip_newline (char *s)
{
    if (s != NULL && *s != '\0') {
        int i, len = strlen(s);

        for (i=len-1; i>=0; i--) {
            if (s[i] == '\n' || s[i] == '\r') {
                s[i] = '\0';
            } else {
                break;
            }
        }
    }
}

static NODE *getline_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const char *buf = NULL;
        NODE *rs = NULL;

        if (l->vname == NULL) {
            gretl_errmsg_set("getline: the source must be a named string variable");
            p->err = E_INVARG;
        } else {
            buf = l->v.str;
            if (null_node(r)) {
                /* clean-up only */
                bufgets_finalize(buf);
                ret->v.xval = 0;
                return ret;
            }
            if (r->t == STR && r->vname != NULL) {
                rs = r;
            } else if (r->t == U_ADDR && r->L->t == STR) {
                rs = ptr_node_get_referent_node(r, p);
            } else {
                gretl_errmsg_set("getline: the target must be a named string variable");
                p->err = E_INVARG;
            }
        }

        if (!p->err) {
            p->err = query_bufgets_init(buf);
        }
        if (!p->err) {
            size_t len = bufgets_peek_line_length(buf);

            if (len == 0) {
                bufgets_finalize(buf);
                rs->v.str = user_string_reset(rs->vname, NULL, &p->err);
                ret->v.xval = 0;
            } else {
                rs->v.str = user_string_resize(rs->vname, len, &p->err);
                if (!p->err) {
                    bufgets(rs->v.str, len, buf);
                    strip_newline(rs->v.str);
                    ret->v.xval = 1;
                }
            }
        }
    }

    return ret;
}

static int series_get_start (int t1, int t2, const double *x)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            break;
        }
    }

    return t + 1;
}

static int series_get_end (int t1, int t2, const double *x)
{
    int t;

    for (t=t2; t>=t1; t--) {
        if (!na(x[t])) {
            break;
        }
    }

    return t + 1;
}

static void cast_to_series (NODE *n, int f, gretl_matrix **tmp,
                            int *t1, int *t2, parser *p)
{
    gretl_matrix *m = n->v.m;
    int len = gretl_vector_get_length(m);

    if (gretl_is_null_matrix(m)) {
        p->err = E_DATA;
    } else if (m->is_complex) {
        node_type_error(f, 1, SERIES, n, p);
    } else if (len > 0 && len == p->dset->n) {
        *tmp = m;
        n->v.xvec = m->val;
    } else if (len > 0 && t1 != NULL && t2 != NULL) {
        *tmp = m;
        n->v.xvec = m->val;
        *t1 = 0;
        *t2 = len - 1;
    } else {
        node_type_error(f, 1, SERIES, n, p);
    }
}

/* Functions taking a series (or matrix) as argument and returning a
   scalar; allowance is made for an additional boolean arg in some
   cases.
*/

static NODE *series_scalar_func (NODE *n, int f,
                                 NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        gretl_matrix *tmp = NULL;
        int t1 = p->dset->t1;
        int t2 = p->dset->t2;
        const double *x;

        if (n->t == MAT && gretl_is_complex(n->v.m)) {
            p->err = E_CMPLX;
        } else if (n->t == MAT) {
            if (f == F_SUM || f == F_MAX || f == F_MIN) {
                /* we'll sum, max, or min all elements of a matrix */
                if (f == F_SUM) {
                    ret->v.xval = gretl_matrix_global_sum(n->v.m, &p->err);
                } else {
                    int mm = (f == F_MAX);

                    ret->v.xval = gretl_matrix_global_minmax(n->v.m, mm,
                                                             &p->err);
                }
                return ret; /* handled */
            } else if (f == F_T1 || f == F_T2) {
                cast_to_series(n, f, &tmp, NULL, NULL, p);
            } else {
                cast_to_series(n, f, &tmp, &t1, &t2, p);
            }
            if (p->err) {
                return NULL;
            }
        }

        if (f == F_T1 || f == F_T2) {
            int insample = node_get_bool(r, p, 0);

            if (p->err) {
                return NULL;
            } else if (!insample) {
                t1 = 0;
                t2 = p->dset->n - 1;
            }
        }

        x = n->v.xvec;

        switch (f) {
        case F_SUM:
            ret->v.xval = gretl_sum(t1, t2, x);
            break;
        case F_SUMALL:
            ret->v.xval = series_sum_all(t1, t2, x);
            break;
        case F_MEAN:
            ret->v.xval = gretl_mean(t1, t2, x);
            break;
        case F_SD:
            ret->v.xval = gretl_stddev(t1, t2, x);
            break;
        case F_VCE:
            ret->v.xval = gretl_variance(t1, t2, x);
            break;
        case F_SST:
            ret->v.xval = gretl_sst(t1, t2, x);
            break;
        case F_SKEWNESS:
            ret->v.xval = gretl_skewness(t1, t2, x);
            break;
        case F_KURTOSIS:
            ret->v.xval = gretl_kurtosis(t1, t2, x);
            break;
        case F_MIN:
            ret->v.xval = gretl_min(t1, t2, x);
            break;
        case F_MAX:
            ret->v.xval = gretl_max(t1, t2, x);
            break;
        case F_MEDIAN:
            ret->v.xval = gretl_median(t1, t2, x);
            break;
        case F_GINI:
            ret->v.xval = gretl_gini(t1, t2, x);
            break;
        case F_NOBS:
            ret->v.xval = series_get_nobs(t1, t2, x);
            break;
        case F_ISCONST:
            ret->v.xval = gretl_isconst(t1, t2, x);
            break;
        case F_ISDUMMY:
            ret->v.xval = gretl_isdummy(t1, t2, x);
            break;
        case F_T1:
            ret->v.xval = series_get_start(t1, t2, x);
            break;
        case F_T2:
            ret->v.xval = series_get_end(t1, t2, x);
            break;
        default:
            break;
        }

        if (n->t == MAT) {
            n->v.m = tmp;
        }
    }

    return ret;
}

/* Functions normally taking a series or vector as argument and
   returning a scalar, but are evaluated on a scalar, so output is
   trivial.
*/

static NODE *pretend_matrix_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        if (f == F_SUM || f == F_SUMALL || f == F_MEAN || f == F_MAX ||
            f == F_MIN || f == F_MEDIAN) {
            ret->v.xval = n->v.xval;
        } else if (f == F_NOBS) {
            ret->v.xval = 1;
        } else if (f == F_SD || f == F_VCE || f == F_SST || f == F_GINI) {
            ret->v.xval = 0;
        } else if (f == F_SKEWNESS || f == F_KURTOSIS) {
            /* this is probably less intuitive than all the above:
               in "normal" cases we're returning NADBL when the variance
               is 0, so we're just being consistent here
            */
            ret->v.xval = NADBL;
        } else {
            /* any other cases not legit */
            node_type_error(f, 0, SERIES, n, p);
        }
    }

    return ret;
}

static gretlopt get_normtest_option (NODE *n, parser *p)
{
    gretlopt opt = OPT_NONE;

    if (null_node(n)) {
        ; /* OK */
    } else if (n->t == STR) {
        const char *s = n->v.str;

        if (!strcmp(s, "swilk")) {
            opt = OPT_W;
        } else if (!strcmp(s, "lillie")) {
            opt = OPT_L;
        } else if (!strcmp(s, "jbera")) {
            opt = OPT_J;
        } else if (!strcmp(s, "dhansen")) {
            opt = OPT_D;
        } else if (!strcmp(s, "all")) {
            opt = OPT_A;
        } else {
            p->err = E_INVARG;
        }
    } else {
        p->err = e_types(n);
    }

    return opt;
}

static NODE *series_matrix_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const double *x = NULL;
        int n = 0, t1 = 0, t2 = 0;
        gretlopt opt = OPT_NONE;

        if (f == F_NORMTEST) {
            opt = get_normtest_option(r, p);
            if (p->err) {
                return NULL;
            }
        }

        if (l->t == MAT) {
            n = gretl_vector_get_length(l->v.m);
            if (n == 0) {
                p->err = E_TYPES;
            } else {
                x = l->v.m->val;
                t2 = n - 1;
            }
        } else if (f == F_NORMTEST) {
            x = l->v.xvec;
            t1 = p->dset->t1;
            t2 = p->dset->t2;
        } else {
            x = l->v.xvec + p->dset->t1;
            n = sample_size(p->dset);
        }

        if (!p->err) {
            ret = aux_matrix_node(p);
        }

        if (!p->err) {
            if (f == F_NORMTEST) {
                ret->v.m = gretl_normtest_matrix(x, t1, t2, opt, &p->err);
            } else {
                /* F_ECDF */
                ret->v.m = empirical_cdf(x, n, &p->err);
            }
        }
    }

    return ret;
}

/* There must be a matrix in @l; @r may hold a vector or
   a scalar value */

static NODE *matrix_quantiles_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        gretl_matrix *pmat = node_get_real_matrix(r, p, 0, 0);

        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (ret != NULL) {
            ret->v.m = gretl_matrix_quantiles(l->v.m, pmat, &p->err);
        }
    }

    return ret;
}

/* functions taking a series and a scalar as arguments and returning
   a scalar
*/

static NODE *series_scalar_scalar_func (NODE *l, NODE *r,
                                        NODE *r2, int f,
                                        parser *p)
{
    NODE *ret = NULL;

    if (f == F_LRVAR && dataset_is_panel(p->dset)) {
        /* FIXME */
        gretl_errmsg_set("lrvar(): not yet supported for panel data");
        p->err = E_DATA;
        return ret;
    }

    if (starting(p)) {
        double rval = -1;
        double r2val = NADBL;
        const double *xvec;
        int t1 = p->dset->t1;
        int t2 = p->dset->t2;
        int pd = 1;

        if (f == F_LRVAR && null_node(r)) {
            ; /* OK, second arg is optional */
        } else {
            rval = node_get_scalar(r, p);
        }

        if (l->t == NUM) {
            t1 = 0;
            t2 = 0;
            xvec = &l->v.xval;
        } else if (l->t == MAT) {
            int n = gretl_vector_get_length(l->v.m);

            if (n == 0) {
                p->err = E_TYPES;
                return NULL;
            }
            t1 = 0;
            t2 = n - 1;
            xvec = l->v.m->val;
        } else {
            /* got a series on the left */
            pd = p->dset->pd;
            xvec = l->v.xvec;
        }

        if (f == F_LRVAR && !null_node(r2)) {
            /* optional third arg */
            r2val = node_get_scalar(r2, p);
        }

        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (p->err) {
            return ret;
        }

        switch (f) {
        case F_LRVAR:
            ret->v.xval = gretl_long_run_variance(t1, t2, xvec, (int) rval, r2val);
            break;
        case F_QUANTILE:
            ret->v.xval = gretl_quantile(t1, t2, xvec, rval, OPT_NONE, &p->err);
            break;
        case F_NPV:
            ret->v.xval = gretl_npv(t1, t2, xvec, rval, pd, &p->err);
            break;
        case F_ISCONST:
            ret->v.xval = panel_isconst(t1, t2, pd, xvec, (int) rval);
            break;
        default:
            break;
        }

    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *isconst_or_dum_node (NODE *l, NODE *r, parser *p, int f)
{
    if (f == F_ISDUMMY || null_node(r)) {
        return series_scalar_func(l, f, NULL, p);
    } else if (l->t == MAT) {
        node_type_error(f, 1, SERIES, l, p);
        return NULL;
    } else if (!dataset_is_panel(p->dset)) {
        p->err = E_PDWRONG;
        return NULL;
    } else {
        return series_scalar_scalar_func(l, r, NULL, f, p);
    }
}

/* Series on left, scalar or string on right, as in
   x[23] or somevar["CA"]. We return the selected
   scalar value from the series -- unless the series
   is string-valued, in which case we return the
   string value for the given observation.
*/

static NODE *series_obs (NODE *l, NODE *r, parser *p)
{
    int strval = strvals_node(l, p);
    NODE *ret;

    ret = strval ? aux_string_node(p) : aux_scalar_node(p);

    if (ret != NULL) {
        int t = -1; /* invalid */

        if (r->t == STR) {
            t = dateton(r->v.str, p->dset);
            if (t < 0) {
                if (dataset_has_markers(p->dset)) {
                    gretl_errmsg_sprintf_replace(_("Invalid observation specifier \"%s\""),
                                                 r->v.str);
                }
                p->err = E_DATA;
            }
        } else {
            /* plain integer */
            t = node_get_int(r, p);
            if (!p->err) {
                if (t > 0 && t <= p->dset->n) {
                    t--; /* OK, convert to zero based */
                } else {
                    gretl_errmsg_sprintf(_("Index value %d is out of bounds"), t);
                    p->err = E_DATA;
                }
            }
        }

        if (!p->err) {
            if (strval) {
                const char *s =
                    series_get_string_for_obs(p->dset, l->vnum, t);

                if (s == NULL) {
                    ret->v.str = gretl_strdup("");
                } else {
                    ret->v.str = gretl_strdup(s);
                }
            } else {
                ret->v.xval = l->v.xvec[t];
            }
        }
    }

    return ret;
}

static NODE *series_ljung_box (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = l->v.xvec;
        int k = node_get_int(r, p);
        int t1 = p->dset->t1;
        int t2 = p->dset->t2;

        if (!p->err && k <= 0) {
            gretl_errmsg_sprintf(_("Invalid lag order %d"), k);
            p->err = E_DATA;
        }

        if (!p->err) {
            p->err = series_adjust_sample(x, &t1, &t2);
        }

        if (!p->err) {
            ret->v.xval = ljung_box(k, t1, t2, x, &p->err);
        }
    }

    return ret;
}

static NODE *series_polyfit (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = l->v.xvec;
        int order = node_get_int(r, p);

        if (!p->err) {
            p->err = poly_trend(x, ret->v.xvec, p->dset, order);
        }
    }

    return ret;
}

static NODE *series_lag (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    const double *x = l->v.xvec;
    int k = -(node_get_int(r, p));
    int t1, t2;

    if (!p->err && l->vnum == 0) {
        gretl_errmsg_set(_("The constant cannot be lagged"));
        p->err = E_TYPES;
    }

    if (!p->err) {
        ret = aux_series_node(p);
    }

    if (ret == NULL) {
        return NULL;
    }

    t1 = autoreg(p) ? p->obs : p->dset->t1;
    t2 = autoreg(p) ? p->obs : p->dset->t2;

    lag_calc(ret->v.xvec, x, k, t1, t2, B_ASN, 1.0, p);

    return ret;
}

static NODE *series_sort_by (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        if (l->t == SERIES && r->t == SERIES) {
            p->err = gretl_sort_by(l->v.xvec, r->v.xvec, ret->v.xvec, p->dset);
        } else if (l->t != SERIES) {
            p->err = e_types(l);
        } else {
            p->err = e_types(r);
        }
    }

    return ret;
}

static NODE *vector_sort (NODE *l, int f, parser *p)
{
    NODE *ret = (l->t == SERIES)? aux_series_node(p) :
        aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        if (l->t == SERIES) {
            p->err = sort_series(l->v.xvec, ret->v.xvec, f, p->dset);
        } else if (l->t == NUM) {
            ret->v.m = gretl_matrix_from_scalar(l->v.xval);
        } else if (gretl_is_null_matrix(l->v.m)) {
            ret->v.m = gretl_null_matrix_new();
        } else if (l->v.m->is_complex) {
            p->err = E_CMPLX;
        } else {
            int descending = (f == F_DSORT);

            ret->v.m = gretl_vector_sort(l->v.m, descending, &p->err);
        }
    }

    return ret;
}

static NODE *vector_values (NODE *l, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = NULL;
        int n = 0;

        if (l->t == NUM) {
            n = 1;
            x = &l->v.xval;
        } else if (l->t == SERIES) {
            n = sample_size(p->dset);
            x = l->v.xvec + p->dset->t1;
        } else if (gretl_is_null_matrix(l->v.m)) {
            ret->v.m = gretl_null_matrix_new();
        } else {
            n = gretl_vector_get_length(l->v.m);
            x = l->v.m->val;
        }

        if (n > 0 && x != NULL) {
            gretlopt opt = (f == F_VALUES)? OPT_S : OPT_NONE;

            ret->v.m = gretl_matrix_values(x, n, opt, &p->err);
        } else if (ret->v.m == NULL) {
            p->err = E_DATA;
        }
    }

    return ret;
}

static NODE *do_irr (NODE *l, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const double *x = NULL;
        int pd = 1, n = 0;

        if (l->t == NUM) {
            n = 1;
            x = &l->v.xval;
        } else if (l->t == SERIES) {
            n = sample_size(p->dset);
            x = l->v.xvec + p->dset->t1;
            pd = p->dset->pd;
        } else if (!gretl_is_null_matrix(l->v.m)) {
            n = gretl_vector_get_length(l->v.m);
            x = l->v.m->val;
        }

        if (n > 0 && x != NULL) {
            ret->v.xval = gretl_irr(x, n, pd, &p->err);
        } else {
            p->err = E_DATA;
        }
    }

    return ret;
}

static int deseasonalize (const double *x, double *y,
                          const char *vname, int f,
                          NODE *r, parser *p)
{
    gretl_bundle *b = NULL;
    int tramo = 0;
    int err = 0;

    if (!null_node(r)) {
        if (r->t == STR) {
            /* legacy usage */
            if (!strcmp(r->v.str, "T")) {
                tramo = 1;
            } else if (strcmp(r->v.str, "X")) {
                err = E_INVARG;
            }
        } else {
            /* current style */
            b = node_get_bundle(r, p);
        }
    }
    if (!err) {
        err = seasonally_adjust_series(x, y, vname, p->dset,
                                       tramo, b, p->prn);
    }

    return err;
}

#define is_panel_stat(f) (f == F_PNOBS || \
                          f == F_PMIN ||  \
                          f == F_PMAX ||  \
                          f == F_PSUM || \
                          f == F_PMEAN || \
                          f == F_PXSUM ||  \
                          f == F_PXNOBS ||  \
                          f == F_PSD)

/* Functions taking a series as argument and returning a series.  Note
   that the @r node may be null or may contain an auxiliary parameter,
   as follows:

   @f == F_DESEAS: bundle or (backward compat) string
   @f == panel stats function: series
   otherwise: scalar

   In case @f == F_HPFILT, the @o node can be used to pass
   an additional optional argument.
*/

static NODE *series_series_func (NODE *l, NODE *r, NODE *o,
                                 int f, parser *p)
{
    NODE *ret = NULL;

    if ((f == F_SDIFF || f == F_DESEAS) &&
        !dataset_is_seasonal(p->dset)) {
        p->err = E_PDWRONG;
        return NULL;
    } else if ((f == F_FRACDIFF || f == F_FRACLAG ||
                f == F_HPFILT || f == F_TRAMOLIN) &&
               !dataset_is_time_series(p->dset)) {
        gretl_errmsg_set(_("This function requires time-series data"));
        p->err = E_DATA;
        return NULL;
    }

    /* check the @r argument */
    if ((f == F_FRACDIFF || f == F_FRACLAG) && null_node(r)) {
        /* fractional difference needed */
        p->err = E_ARGS;
    } else if (null_node(r)) {
        ; /* not present, OK */
    } else if (f == F_DESEAS) {
        if (r->t != STR && r->t != BUNDLE && r->t != U_ADDR) {
            p->err = e_types(r);
        }
    } else if (is_panel_stat(f)) {
        if (r->t != SERIES) {
            node_type_error(f, 2, SERIES, r, p);
        }
    } else if (!scalar_node(r)) {
        node_type_error(f, 2, NUM, r, p);
    }

    if (l->t == MAT && f == F_BOXCOX) {
        /* Do all columns of matrix input: this could
           be generalized to some other functions?
        */
        double d = node_get_scalar(r, p);

        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = boxcox_matrix(l->v.m, d, &p->err);
        }
        return ret;
    }

    ret = aux_series_node(p);

    if (ret != NULL) {
        gretl_matrix *tmp = NULL;
        const char *vname = NULL;
        double parm = NADBL;
        int oneside = 0;
        const double *z = NULL;
        const double *x;
        double *y;

        if (l->t == MAT) {
            cast_to_series(l, f, &tmp, NULL, NULL, p);
        } else if (l->vnum > 0) {
            vname = p->dset->varname[l->vnum];
        }
        if (!null_node(r)) {
            if (r->t == SERIES) {
                z = r->v.xvec;
            } else if (scalar_node(r)) {
                parm = node_get_scalar(r, p);
            }
        }
        if (p->err) {
            return NULL;
        }

        x = l->v.xvec;
        y = ret->v.xvec;

        switch (f) {
        case F_HPFILT:
            oneside = node_get_bool(o, p, 0);
            if (!p->err && oneside) {
                p->err = oshp_filter(x, y, p->dset, parm, OPT_NONE);
            } else if (!p->err) {
                p->err = hp_filter(x, y, p->dset, parm, OPT_NONE);
            }
            break;
        case F_FRACDIFF:
            p->err = fracdiff_series(x, y, parm, 1, autoreg(p) ? p->obs : -1, p->dset);
            break;
        case F_FRACLAG:
            p->err = fracdiff_series(x, y, parm, 0, autoreg(p) ? p->obs : -1, p->dset);
            break;
        case F_BOXCOX:
            p->err = boxcox_series(x, y, parm, p->dset);
            break;
        case F_DIFF:
        case F_LDIFF:
        case F_SDIFF:
            p->err = diff_series(x, y, f, p->dset);
            break;
        case F_ODEV:
            p->err = orthdev_series(x, y, p->dset);
            break;
        case F_CUM:
            p->err = cum_series(x, y, p->dset);
            break;
        case F_DESEAS:
            p->err = deseasonalize(x, y, vname, f, r, p);
            break;
        case F_TRAMOLIN:
            p->err = tramo_linearize_series(x, y, p->dset);
            break;
        case F_INTERPOL:
            p->err = interpolate_series(x, y, p->dset);
            break;
        case F_RESAMPLE:
            if (!na(parm)) {
                p->err = block_resample_series(x, y, parm, p->dset);
            } else {
                p->err = resample_series(x, y, p->dset);
            }
            break;
        case F_PNOBS:
        case F_PMIN:
        case F_PMAX:
        case F_PSUM:
        case F_PMEAN:
        case F_PXSUM:
        case F_PXNOBS:
        case F_PSD:
            p->err = panel_statistic(x, y, p->dset, f, z);
            break;
        case F_RANKING:
            p->err = rank_series(x, y, F_SORT, p->dset);
            break;
        default:
            break;
        }

        if (l->t == MAT) {
            l->v.m = tmp;
        }
    }

    return ret;
}

static NODE *do_panel_shrink (NODE *l, int noskip, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        ret->v.m = panel_shrink(l->v.xvec, noskip, p->dset, &p->err);
    }

    return ret;
}

static NODE *do_panel_expand (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_series_node(p);

    if (ret != NULL && starting(p)) {
        gretlopt opt = OPT_NONE;

        if (node_get_bool(r, p, 0)) {
            opt = OPT_X; /* expand by unit, not time */
        }
        if (!p->err) {
            p->err = panel_expand(l->v.m, ret->v.xvec, opt, p->dset);
        }
    }

    return ret;
}

/* pergm function takes series or column vector arg, returns matrix:
   if we come up with more functions on that pattern, the following
   could be extended
*/

static NODE *pergm_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (!null_or_scalar(r)) {
        /* optional 'r' node must be scalar */
        node_type_error(F_PERGM, 2, NUM, r, p);
    } else if (l->t == MAT && gretl_vector_get_length(l->v.m) == 0) {
        /* if 'l' node is not a series, must be a vector */
        node_type_error(F_PERGM, 1, SERIES, l, p);
    } else {
        ret = aux_matrix_node(p);
    }

    if (!p->err) {
        const double *x = NULL;
        int t1 = 0, t2 = 0;
        int width = -1;

        if (l->t == SERIES) {
            x = l->v.xvec;
            t1 = p->dset->t1;
            t2 = p->dset->t2;
        } else if (l->t == MAT) {
            x = l->v.m->val;
            t1 = 0;
            t2 = gretl_vector_get_length(l->v.m) - 1;
        }

        if (!null_node(r)) {
            width = node_get_int(r, p);
        }

        if (!p->err) {
            ret->v.m = periodogram_matrix(x, t1, t2, width, &p->err);
        }
    }

    return ret;
}

static void *get_complex_counterpart (void *func)
{
    if (func == acos) return cacos;
    if (func == asin) return casin;
    if (func == atan) return catan;
    if (func == cos)  return ccos;
    if (func == sin)  return csin;
    if (func == tan)  return ctan;
    if (func == acosh) return cacosh;
    if (func == asinh) return casinh;
    if (func == atanh) return catanh;
    if (func == cosh) return ccosh;
    if (func == sinh) return csinh;
    if (func == tanh) return ctanh;
    if (func == exp)  return cexp;
    if (func == log)  return clog;
    if (func == sqrt) return csqrt;

    return NULL;
}

#define cmplx_to_double(f) (f == F_CARG || f == F_CMOD || \
                            f == F_REAL || f == F_IMAG || \
                            f == F_CQUAD || f == F_ABS)

/* application of scalar function to each element of matrix */

static NODE *apply_matrix_func (NODE *t, NODE *f, parser *p)
{
    gretl_matrix *m = NULL;
    gretl_matrix tmp = {0};
    int ret_complex = 0;
    NODE *ret;

    if (t->t == NUM) {
        m = &tmp;
        gretl_matrix_init(m);
        m->rows = m->cols = 1;
        m->val = &t->v.xval;
    } else {
        m = t->v.m;
        ret_complex = m->is_complex;
    }

    if (m->is_complex && cmplx_to_double(f->t)) {
        ret_complex = 0;
    }

    ret = aux_sized_matrix_node(p, m->rows, m->cols, ret_complex);
    if (ret == NULL) {
        return ret;
    }

    if (m->is_complex) {
        if (f->t == F_ABS) {
            apply_cmatrix_dfunc(ret->v.m, m, cabs);
        } else if (!ret_complex) {
            apply_cmatrix_dfunc(ret->v.m, m, f->v.ptr);
        } else if (f->t == F_CONJ) {
            apply_cmatrix_cfunc(ret->v.m, m, conj);
        } else if (f->t == U_NEG || f->t == U_POS || f->t == U_NOT) {
            apply_cmatrix_unary_op(ret->v.m, m, f->t);
        } else {
            double complex (*cfunc) (double complex) = NULL;

            if (f->v.ptr != NULL) {
                cfunc = get_complex_counterpart(f->v.ptr);
            }
            if (cfunc == NULL) {
                /* gatekeeper for complex */
                p->err = function_real_only(f->t);
            } else {
                p->err = apply_cmatrix_cfunc(ret->v.m, m, cfunc);
            }
        }
    } else if (f->t == F_REAL || f->t == F_IMAG) {
        if (f->t == F_REAL) {
            size_t msize = m->rows * m->cols * sizeof(double);

            memcpy(ret->v.m->val, m->val, msize);
        } else {
            int i, n = m->rows * m->cols;

            for (i=0; i<n; i++) {
                ret->v.m->val[i] = 0;
            }
        }
    } else {
        double (*dfunc) (double) = f->v.ptr;
        int i, n = m->rows * m->cols;

        if (dfunc != NULL) {
            for (i=0; i<n && !p->err; i++) {
                ret->v.m->val[i] = dfunc(m->val[i]);
            }
        } else {
            for (i=0; i<n && !p->err; i++) {
                ret->v.m->val[i] = real_apply_func(m->val[i], f->t, p);
            }
        }
    }

#if 0
    if (f->t == F_SQRT) {
        fprintf(stderr, "apply_matrix_func, sqrt, ret = %p, ret->v.m = %p\n",
                (void *) ret, (void *) ret->v.m);
        fprintf(stderr, " input node %p\n", (void *) t);
        gretl_matrix_print(ret->v.m, "ret->v.m");
    }
#endif

    return ret;
}

static gretl_matrix *matrix_from_scalars (NODE **nn, int m,
                                          int nsep, int seppos,
                                          parser *p)
{
    gretl_matrix *M;
    NODE *n;
    int r = nsep + 1;
    int c = (seppos > 0)? seppos : m;
    int nelem = m - nsep;
    double x;
    int i, j, k;

    /* check that all rows are the same length */

    if (nelem != r * c) {
        p->err = E_PARSE;
    } else if (nsep > 0) {
        k = 0;
        for (i=0; i<m; i++) {
            n = nn[i];
            if (null_node(n)) {
                if (i - k != seppos) {
                    p->err = E_PARSE;
                    break;
                }
                k = i + 1;
            }
        }
    }

    if (p->err) {
        pprintf(p->prn, _("Matrix specification is not coherent"));
        pputc(p->prn, '\n');
        return NULL;
    }

#if EDEBUG
    fprintf(stderr, "matrix_from_scalars: m=%d, nsep=%d, seppos=%d, nelem=%d\n",
            m, nsep, seppos, nelem);
#endif

    M = gretl_matrix_alloc(r, c);
    if (M == NULL) {
        p->err = E_ALLOC;
    } else {
        k = 0;
        for (i=0; i<r && !p->err; i++) {
            for (j=0; j<c; j++) {
                n = nn[k++];
                if (null_node(n)) {
                    n = nn[k++];
                }
                x = node_get_scalar(n, p);
                gretl_matrix_set(M, i, j, x);
            }
        }
    }

    return M;
}

static int *full_series_list (const DATASET *dset, int *err)
{
    int *list = NULL;

    if (dset->v < 2) {
        *err = E_DATA;
        return NULL;
    }

    list = gretl_consecutive_list_new(1, dset->v - 1);
    if (list == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    return list;
}

static gretl_matrix *real_matrix_from_list (const int *list,
                                            const DATASET *dset,
                                            parser *p)
{
    gretl_matrix *M = NULL;

    if (list != NULL && list[0] == 0) {
        M = gretl_null_matrix_new();
    } else {
        const gretl_matrix *mmask = get_matrix_mask();

        if (mmask != NULL) {
            M = gretl_matrix_data_subset_special(list, dset,
                                                 mmask, &p->err);
        } else {
            int missop = M_MISSING_OK;

            if (libset_get_bool(SKIP_MISSING)) {
                missop = M_MISSING_SKIP;
            }
            M = gretl_matrix_data_subset(list, dset, dset->t1, dset->t2,
                                         missop, &p->err);
        }
    }

    return M;
}

static gretl_matrix *matrix_from_list (NODE *n, parser *p)
{
    gretl_matrix *M = NULL;
    int *list = NULL;
    int freelist = 0;

    if (n != NULL) {
        if (n->t == LIST) {
            list = n->v.ivec;
        } else {
            p->err = E_DATA;
        }
    } else {
        list = full_series_list(p->dset, &p->err);
        freelist = 1;
    }

    if (!p->err) {
        M = real_matrix_from_list(list, p->dset, p);
    }

    if (freelist) {
        free(list);
    }

    return M;
}

static void *arg_get_data (NODE *n, int ref, GretlType *type,
                           user_var **puv)
{
    void *data = NULL;

    *puv = n->uv;

    if (n->t == SERIES) {
        if (ref) {
            *type = GRETL_TYPE_SERIES_REF;
            data = &n->vnum;
        } else if (n->vname != NULL) {
            /* FIXME conditionality here? */
            *type = GRETL_TYPE_USERIES;
            data = &n->vnum;
        } else {
            *type = GRETL_TYPE_SERIES;
            data = n->v.xvec;
        }
    } else if (n->t == NUM) {
        *type = ref ? GRETL_TYPE_SCALAR_REF : GRETL_TYPE_DOUBLE;
        data = &n->v.xval;
    } else if (n->t == MAT) {
        *type = ref ? GRETL_TYPE_MATRIX_REF : GRETL_TYPE_MATRIX;
        data = n->v.m;
    } else if (n->t == BUNDLE) {
        *type = ref ? GRETL_TYPE_BUNDLE_REF : GRETL_TYPE_BUNDLE;
        data = n->v.b;
    } else if (n->t == ARRAY) {
        *type = ref ? GRETL_TYPE_ARRAY_REF : GRETL_TYPE_ARRAY;
        data = n->v.a;
    } else if (n->t == STR) {
        *type = ref ? GRETL_TYPE_STRING_REF : GRETL_TYPE_STRING;
        data = n->v.str;
    } else if (n->t == LIST) {
        *type = GRETL_TYPE_LIST;
        data = n->v.ivec;
    } else if (n->t == SUB_ADDR) {
        data = sub_addr_get_data(n, type, puv);
    } else {
        *type = GRETL_TYPE_NONE;
    }

    return data;
}

static NODE *suitable_ufunc_ret_node (parser *p,
                                      GretlType t)
{
    if (t == GRETL_TYPE_DOUBLE) {
        return aux_scalar_node(p);
    } else if (t == GRETL_TYPE_SERIES) {
        return aux_empty_series_node(p);
    } else if (t == GRETL_TYPE_MATRIX) {
        return aux_matrix_node(p);
    } else if (t == GRETL_TYPE_LIST) {
        return aux_list_node(p);
    } else if (t == GRETL_TYPE_STRING) {
        return aux_string_node(p);
    } else if (t == GRETL_TYPE_BUNDLE) {
        return aux_bundle_node(p);
    } else if (gretl_array_type(t)) {
        return aux_array_node(p);
    } else {
        p->err = E_TYPES;
        return NULL;
    }
}

#define ok_ufunc_sym(s) (s == NUM || s == SERIES || s == MAT || \
                         s == LIST || s == U_ADDR || s == DUM || \
                         s == STR || s == EMPTY || s == BUNDLE || \
                         s == ARRAY || s == SUB_ADDR)

/* evaluate a user-defined function */

static NODE *eval_ufunc (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    const char *funname = l->vname;
    ufunc *uf = l->v.ptr;
    fncall *fc = NULL;
    GretlType rtype = 0;
    int i, nparam, argc = 0;

    rtype = user_func_get_return_type(uf);

#if EDEBUG
    fprintf(stderr, "eval_ufunc (%s, depth %d, rtype %s):\n", funname,
            gretl_function_depth(), gretl_type_get_name(rtype));
    fprintf(stderr, " compiling %d, exec'ing %d\n", (p->flags & P_COMPILE)? 1 : 0,
            (p->flags & P_EXEC)? 1 : 0);
#endif

    if (!p->err && rtype == GRETL_TYPE_VOID) {
        if (p->targ == UNK && p->lh.name[0] == '\0' && p->lh.expr == NULL) {
            /* never reached? */
            p->targ = EMPTY;
        } else if (p->targ != EMPTY) {
            gretl_errmsg_sprintf(_("The function %s does not return any value"),
                                 funname);
            p->err = E_TYPES;
        }
    }

    if (!p->err) {
        /* get the argument and param counts */
        argc = r->v.bn.n_nodes;
        nparam = fn_n_params(uf);
        if (argc > nparam) {
            gretl_errmsg_sprintf(_("Number of arguments (%d) does not "
                                   "match the number of\nparameters for "
                                   "function %s (%d)"),
                                 argc, funname, nparam);
            p->err = E_DATA;
        }
    }

    if (p->err) {
        /* no sense in continuing */
        return NULL;
    }

    if (l == p->tree && (p->flags & P_CATCH)) {
        gretl_warnmsg_set(_("\"catch\" should not be used on calls to "
                            "user-defined functions"));
    }

    fc = user_func_get_fncall(uf);
    if (fc == NULL) {
        p->err = E_ALLOC;
        return NULL;
    }

    /* check the function argument nodes */

    for (i=0; i<argc && !p->err; i++) {
        NODE *arg = r->v.bn.n[i];
        int reftype = 0;

        if (!ok_ufunc_sym(arg->t)) {
            gretl_errmsg_sprintf("%s: invalid argument type %s", funname,
                                 typestr(arg->t));
            p->err = E_TYPES;
        }
#if EDEBUG
        fprintf(stderr, "%s: arg %d is of type %s (err=%d)\n",
                funname, i+1, arg == NULL ? "?" : getsymb(arg->t),
                p->err);
#endif
        if (p->err) {
            break;
        }
        if (arg->t == U_ADDR) {
            /* address node: switch to the 'content' sub-node */
            arg = ptr_node_get_referent_node(arg, p);
            if (!p->err) {
                reftype = 1;
            }
        } else if (arg->t == DUM && arg->v.idnum != DUM_NULL) {
            p->err = E_TYPES;
        }
        if (!p->err) {
            /* assemble info and push argument */
            user_var *uv = NULL;
            GretlType argt = 0;
            void *data;

            data = arg_get_data(arg, reftype, &argt, &uv);
            p->err = push_function_arg(fc, arg->vname, uv, argt, data);
        }
        if (p->err) {
            fprintf(stderr, "%s: error evaluating arg %d\n", funname, i+1);
        }
    }

    /* try sending args to function */

    if (!p->err) {
        char **pdescrip = NULL;
        series_table **ptab = NULL;
        void *altp = NULL;
        void *retp = &altp;

        /* special cases */
        if (rtype == GRETL_TYPE_VOID) {
            retp = NULL;
        } else if (rtype == GRETL_TYPE_MATRIX) {
            if (p->targ == UNK && p->tree == l) {
                /* target type not specified, and function returns
                   a matrix -> set target type to matrix
                */
                p->targ = MAT;
            }
        } else if (rtype == GRETL_TYPE_LIST) {
            if (p->targ == EMPTY && p->tree == l) {
                /* this function offers a list return, but the
                   caller hasn't assigned it and it's not being
                   used as an argument to a further function, so
                   ignore the return value
                */
                retp = NULL;
            }
        }

        if ((p->flags & P_UFRET) && rtype == GRETL_TYPE_SERIES) {
            /* arrange to pick up description of generated series, if any */
            pdescrip = &p->lh.label;
            /* and perhaps string values attached to series */
            ptab = &p->lh.stab;
        }

        p->err = gretl_function_exec_full(fc, rtype, p->dset, retp,
                                          pdescrip, ptab, p->prn);

        if (rtype == GRETL_TYPE_NUMERIC) {
            /* determine which numeric type we actually got */
            rtype = fncall_get_return_type(fc);
        }

        if (!p->err && retp != NULL) {
            ret = suitable_ufunc_ret_node(p, rtype);
        }

        if (!p->err && ret != NULL) {
            if (rtype == GRETL_TYPE_DOUBLE) {
                ret->v.xval = *(double *) altp;
            } else if (rtype == GRETL_TYPE_SERIES) {
                if (ret->v.xvec != NULL) {
                    free(ret->v.xvec);
                }
                ret->v.xvec = altp;
            } else if (rtype == GRETL_TYPE_MATRIX) {
                if (is_tmp_node(ret)) {
                    gretl_matrix_free(ret->v.m);
                }
                ret->v.m = altp;
            } else if (rtype == GRETL_TYPE_LIST) {
                if (is_tmp_node(ret)) {
                    free(ret->v.ivec);
                }
                if (altp != NULL) {
                    ret->v.ivec = altp;
                } else {
                    ret->v.ivec = gretl_list_new(0);
                }
            } else if (rtype == GRETL_TYPE_STRING) {
                if (is_tmp_node(ret)) {
                    free(ret->v.str);
                }
                ret->v.str = altp;
            } else if (rtype == GRETL_TYPE_BUNDLE) {
                if (is_tmp_node(ret)) {
                    gretl_bundle_destroy(ret->v.b);
                }
                ret->t = BUNDLE;
                ret->v.b = altp;
            } else if (gretl_array_type(rtype)) {
                if (is_tmp_node(ret)) {
                    gretl_array_destroy(ret->v.a);
                }
                ret->t = ARRAY;
                ret->v.a = altp;
            }
        }
    }

#if EDEBUG
    fprintf(stderr, "return from %s, depth %d, ret = %p\n",
            funname, gretl_function_depth(), (void *) ret);
#endif

    return ret;
}

#ifdef USE_RLIB

/* evaluate an R function */

static NODE *eval_Rfunc (NODE *l, NODE *r, parser *p)
{
    int i, argc = r->v.bn.n_nodes;
    const char *funname = l->vname;
    int rtype = GRETL_TYPE_NONE;
    NODE *ret = NULL;

    /* first find the function */
    p->err = gretl_R_get_call(funname, argc);
    if (p->err) {
        fprintf(stderr, "eval_Rfunc: can't find function %s\n", funname);
        return NULL;
    }

    /* evaluate the function arguments */
    for (i=0; i<argc && !p->err; i++) {
        NODE *arg = r->v.bn.n[i];
        GretlType type;

        if (arg == NULL) {
            fprintf(stderr, "%s: arg %d is NULL!\n", funname, i);
            p->err = E_DATA;
        } else if (!ok_ufunc_sym(arg->t)) {
            p->err = e_types(arg);
        }
        if (p->err) {
            break;
        }

#if EDEBUG
        fprintf(stderr, "%s: arg[%d] is of type %d\n", funname, i, arg->t);
#endif
        type = gretl_type_from_gen_type(arg->t);

        if (type == GRETL_TYPE_SERIES) {
            p->err = gretl_R_function_add_series(arg->v.xvec, p->dset, arg->vnum);
        } else if (type == GRETL_TYPE_DOUBLE) {
            p->err = gretl_R_function_add_arg(&arg->v.xval, type);
        } else if (type == GRETL_TYPE_MATRIX) {
            p->err = gretl_R_function_add_arg(arg->v.m, type);
        } else if (type == GRETL_TYPE_ARRAY) {
            p->err = gretl_R_function_add_arg(arg->v.a, type);
        } else if (type == GRETL_TYPE_STRING) {
            p->err = gretl_R_function_add_arg(arg->v.str, type);
        } else if (type == GRETL_TYPE_BUNDLE) {
            p->err = gretl_R_function_add_arg(arg->v.b, type);
        } else {
            fprintf(stderr, "eval_Rfunc: argument not supported\n");
            p->err = E_TYPES;
            return NULL;
        }
        if (p->err) {
            fprintf(stderr, "eval_Rfunc: error evaluating arg %d\n", i);
        }
    }

    /* try sending args to function */

    if (!p->err) {
        double xret = NADBL;
        void *retp = &xret;

        p->err = gretl_R_function_exec(funname, &rtype, &retp);

        if (!p->err) {
            if (gretl_scalar_type(rtype)) {
                ret = aux_scalar_node(p);
                if (ret != NULL) {
                    ret->v.xval = xret;
                }
            } else if (rtype == GRETL_TYPE_MATRIX) {
                ret = aux_matrix_node(p);
                if (ret != NULL) {
                    ret->v.m = (gretl_matrix *) retp;
                }
            } else if (rtype == GRETL_TYPE_STRING) {
                ret = aux_string_node(p);
                if (ret != NULL) {
                    if (is_tmp_node(ret)) {
                        free(ret->v.str);
                    }
                    ret->v.str = (char *) retp;
                }
            } else if (rtype == GRETL_TYPE_ARRAY) {
                ret = aux_array_node(p);
                if (ret != NULL) {
                    ret->v.a = (gretl_array *) retp;
                }
            } else if (rtype == GRETL_TYPE_BUNDLE) {
                ret = aux_bundle_node(p);
                if (ret != NULL) {
                    ret->v.b = (gretl_bundle *) retp;
                }
            } else if (rtype == GRETL_TYPE_NONE) {
                ; /* OK? */
            }
        }
    }

#if EDEBUG
    fprintf(stderr, "eval_Rfunc: p->err = %d, ret = %p\n",
            p->err, (void *) ret);
#endif

    return ret;
}

#else

static NODE *eval_Rfunc (NODE *t, NODE *r, parser *p)
{
    gretl_errmsg_set(_("R functions cannot be called via this gretl build"));
    p->err = 1;
    return NULL;
}

#endif /* USE_RLIB or not */

/* Extracting series data @x of length @n from a bundle: whether
   we can get it as a series depends on the current state of the
   dataset. Failing that, we'll get it as a vector.
*/

static NODE *bundled_series_node (parser *p, const double *x,
                                  int n, int *is_tmp)
{
    NODE *ret = NULL;
    int t;

    if (n <= p->dset->n) {
        ret = aux_series_node(p);
        if (!p->err) {
            for (t=p->dset->t1; t<=p->dset->t2 && t<n; t++) {
                ret->v.xvec[t] = x[t];
            }
            *is_tmp = 1;
        }
    } else if (n > 0) {
        ret = aux_matrix_node(p);
        if (!p->err) {
            ret->v.m = gretl_matrix_alloc(n, 1);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
            } else {
                for (t=0; t<n; t++) {
                    ret->v.m->val[t] = x[t];
                }
                *is_tmp = 1;
            }
        }
    } else {
        p->err = E_DATA;
    }

    return ret;
}

/* Getting an object from within a bundle: on the left is the
   bundle reference, on the right should be a string -- the
   key to look up to get content.
*/

static NODE *get_bundle_member (NODE *l, NODE *r, parser *p)
{
    char *key = r->v.str;
    GretlType type;
    int size = 0;
    int is_tmp = 0;
    int gen_t;
    void *val = NULL;
    NODE *ret = NULL;

#if EDEBUG
    fprintf(stderr, "get_bundle_member: %s[\"%s\"], query=%d\n",
            l->vname, key, (p->flags & P_OBJQRY)? 1 : 0);
#endif

    if ((p->flags & P_OBJQRY) && !gretl_bundle_has_key(l->v.b, key)) {
        ret = get_aux_node(p, EMPTY, 0, 0);
        ret->flags |= MUT_NODE;
        return ret;
    }

    val = gretl_bundle_get_element(l->v.b, key, &type, &size,
                                   &is_tmp, &p->err);
    if (p->err) {
        return ret;
    }

    gen_t = gen_type_from_gretl_type(type);
    if (type != GRETL_TYPE_SERIES) {
        ret = get_aux_node(p, gen_t, 0, 0);
    }

    if (type == GRETL_TYPE_INT) {
        ret->v.xval = *(int *) val;
    } else if (type == GRETL_TYPE_DOUBLE) {
        ret->v.xval = *(double *) val;
    } else if (type == GRETL_TYPE_STRING) {
        ret->v.str = (char *) val;
    } else if (type == GRETL_TYPE_MATRIX) {
        ret->v.m = (gretl_matrix *) val;
    } else if (type == GRETL_TYPE_BUNDLE) {
        ret->v.b = (gretl_bundle *) val;
    } else if (type == GRETL_TYPE_ARRAY) {
        ret->v.a = (gretl_array *) val;
    } else if (type == GRETL_TYPE_SERIES) {
        ret = bundled_series_node(p, (const double *) val, size, &is_tmp);
    } else if (type == GRETL_TYPE_LIST) {
        if (stored_list_ok((const int *) val, p->dset)) {
            /* extract as list */
            ret->v.ivec = (int *) val;
            is_tmp = 0;
        } else {
            /* fallback: extract list as row vector */
            ret->t = MAT;
            ret->v.m = gretl_list_to_vector((const int *) val, &p->err);
            is_tmp = 1;
        }
    } else {
        p->err = E_DATA;
    }

    if (ret != NULL) {
        if (is_tmp) {
            ret->flags |= TMP_NODE;
        } else {
            ret->flags &= ~TMP_NODE;
        }
        ret->flags |= MUT_NODE;
    }

    return ret;
}

/* implements the inbundle() function */

static NODE *test_bundle_key (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        gretl_bundle *bundle = l->v.b;
        const char *key = r->v.str;
        GretlType type;

        type = gretl_bundle_get_member_type(bundle, key, NULL);
        ret->v.xval = gretl_type_get_order(type);
    }

    return ret;
}

static NODE *get_bundle_array (NODE *n, int f, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (ret != NULL) {
        if (f == F_GETKEYS) {
            ret->v.a = gretl_bundle_get_keys(n->v.b, &p->err);
        } else {
            /* HF_JBTERMS */
            gretl_array *(*jfunc) (gretl_bundle *, int *);

            jfunc = get_plugin_function("json_bundle_get_terminals");
            if (jfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                ret->v.a = jfunc(n->v.b, &p->err);
            }
        }
    }

    return ret;
}

static const char *optional_bundle_get (gretl_bundle *b,
                                        const char *key,
                                        double *px,
                                        int *err)
{
    const char *s = NULL;

    if (!*err) {
        /* proceed only if we haven't already hit an error */
        if (px != NULL) {
            *px = gretl_bundle_get_scalar(b, key, err);
        } else {
            s = gretl_bundle_get_string(b, key, err);
        }
        if (*err == E_DATA) {
            /* non-existence of item (E_DATA) is OK, but
               wrong type (E_TYPES) is not
            */
            gretl_error_clear();
            *err = 0;
        }
    }

    return s;
}

static NODE *curl_bundle_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

#ifndef USE_CURL
    gretl_errmsg_set(_("Internet access not supported"));
    p->err = E_DATA;
#else
    if (ret != NULL) {
        gretl_bundle *b = NULL;
        int curl_err = 0;

        if (n->t != U_ADDR) {
            p->err = e_types(n);
        } else {
            /* switch to 'content' sub-node */
            n = ptr_node_get_referent_node(n, p);
            if (!p->err) {
                if (n->t != BUNDLE) {
                    p->err = e_types(n);
                } else {
                    b = n->v.b;
                }
            }
        }

        if (!p->err) {
            const char *url = NULL;
            const char *header = NULL;
            const char *postdata = NULL;
            char *output = NULL;
            char *errmsg = NULL;
            double xinclude = 0;

            url = gretl_bundle_get_string(b, "URL", &p->err);
            header = optional_bundle_get(b, "header", NULL, &p->err);
            postdata = optional_bundle_get(b, "postdata", NULL, &p->err);
            optional_bundle_get(b, "include", &xinclude, &p->err);

            if (!p->err) {
                int include = (xinclude == 1.0);

                curl_err = gretl_curl(url, header, postdata, include,
                                      &output, &errmsg);
            }
            if (output != NULL) {
                p->err = gretl_bundle_set_string(b, "output", output);
                free(output);
            } else if (errmsg != NULL) {
                p->err = gretl_bundle_set_string(b, "errmsg", errmsg);
                free(errmsg);
            }
        }

        if (!p->err) {
            ret->v.xval = curl_err;
        }
    }
#endif /* curl supported in libgretl */

    return ret;
}

static NODE *function_depth_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        ret->v.xval = gretl_function_depth();
    }

    return ret;
}

static NODE *lpsolve_bundle_node (NODE *n, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL) {
        gretl_bundle *(*lpfunc) (gretl_bundle *, PRN *, int *);

        lpfunc = get_plugin_function("gretl_lpsolve");
        if (lpfunc == NULL) {
            p->err = E_FOPEN;
        } else {
            ret->v.b = (*lpfunc)(n->v.b, p->prn, &p->err);
        }
    }

    return ret;
}

static gretl_bundle *node_get_bundle (NODE *n, parser *p)
{
    gretl_bundle *b = NULL;

    if (n->t == BUNDLE) {
        b = n->v.b;
    } else if (n->t == U_ADDR) {
        n = ptr_node_get_referent_node(n, p);
        if (!p->err) {
            if (n->t != BUNDLE) {
                p->err = e_types(n);
            } else {
                b = n->v.b;
            }
        }
    } else {
        p->err = e_types(n);
    }

    return b;
}

static NODE *svm_driver_node (NODE *t, NODE *n, parser *p)
{
    NODE *e, *ret = NULL;
    int *list = NULL;
    gretl_bundle *bparm = NULL;
    gretl_bundle *bmod = NULL;
    gretl_bundle *bprob = NULL;
    int i, k = n->v.bn.n_nodes;

    if (k < 2 || k > 4) {
        n_args_error(k, 2, 4, F_SVM, p);
    }

    for (i=0; i<k && !p->err; i++) {
        e = n->v.bn.n[i];
        if (i == 0) {
            list = node_get_list(e, p);
        } else if (i == 1) {
            bparm = node_get_bundle(e, p);
        } else if (i == 2) {
            if (!null_node(e)) {
                bmod = node_get_bundle(e, p);
            }
        } else {
            if (!null_node(e)) {
                bprob = node_get_bundle(e, p);
            }
        }
    }

    if (!p->err) {
        ret = aux_series_node(p);
    }

    if (ret != NULL) {
        int (*pfunc) (const int *, gretl_bundle *,
                      gretl_bundle *, gretl_bundle *,
                      double *, int *, DATASET *, PRN *);
        int got_yhat = 0;

        pfunc = get_plugin_function("gretl_svm_driver");
        if (pfunc == NULL) {
            p->err = E_FOPEN;
        } else {
            p->err = pfunc(list, bparm, bmod, bprob, ret->v.xvec,
                           &got_yhat, p->dset, p->prn);
            if (!p->err && !got_yhat) {
                /* change the return type to scalar NA */
                free(ret->v.xvec);
                ret->t = NUM;
                ret->v.xval = NADBL;
            }
        }
    }

    free(list);

    return ret;
}

static gretl_bundle *bvar_get_bundle (NODE *n, parser *p)
{
    gretl_bundle *b = NULL;

    if (n->v.idnum == B_MODEL) {
        b = bundle_from_model(NULL, p->dset, &p->err);
    } else if (n->v.idnum == B_SYSTEM) {
        b = bundle_from_system(NULL, 0, p->dset, &p->err);
    } else if (n->v.idnum == B_SYSINFO) {
        gretl_bundle *tmp = get_sysinfo_bundle(&p->err);

        if (!p->err) {
            b = gretl_bundle_copy(tmp, &p->err);
        }
    } else if (n->v.idnum == R_RESULT) {
        GretlType type = 0;
        void *ptr = get_last_result_data(&type, &p->err);

        if (type == GRETL_TYPE_BUNDLE) {
            b = ptr;
        } else if (!p->err) {
            p->err = E_TYPES;
        }
    } else {
        p->err = E_DATA;
    }

    return b;
}

static NODE *dollar_bundle_node (NODE *n, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL) {
        ret->v.b = bvar_get_bundle(n, p);
    }

    return ret;
}

static NODE *type_string_node (NODE *n, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL) {
        int t = node_get_int(n, p);

        if (!p->err) {
            const char *s = "null";

            if (t == 1) {
                s = "scalar";
            } else if (t == 2) {
                s = "series";
            } else if (t == 3) {
                s = "matrix";
            } else if (t == 4) {
                s = "string";
            } else if (t == 5) {
                s = "bundle";
            } else if (t == 6) {
                s = "array";
            } else if (t == 7) {
                s = "list";
            }

            ret->v.str = gretl_strdup(s);
            if (ret->v.str == NULL) {
                p->err = E_ALLOC;
            }
        }
    }

    return ret;
}

static double *scalar_to_series (NODE *n, parser *p)
{
    double *ret = NULL;
    int t;

    if (p->dset == NULL || p->dset->n == 0) {
        p->err = E_NODATA;
    } else {
        ret = malloc(p->dset->n * sizeof *ret);
        if (ret == NULL) {
            p->err = E_ALLOC;
        } else {
            for (t=0; t<p->dset->n; t++) {
                if (t >= p->dset->t1 && t <= p->dset->t2) {
                    ret[t] = n->v.xval;
                } else {
                    ret[t] = NADBL;
                }
            }
        }
    }

    return ret;
}

/* We come here only when setting a bundle-member or an element
   of an array -- otherwise we use get_check_return_type().
   Note 2021-08-12: in principle we could allow the case
   where @spec is an array type and @rhs is the singular of
   that type, and support auto-promotion of (e.g.) a string
   to an array of strings. But I'm not sure that's a good idea.
*/

static int lhs_type_check (GretlType spec, GretlType rhs, int t)
{
    int err = 0;

    if (spec != 0 && spec != rhs) {
        if (t == BUNDLE) {
            gretl_errmsg_sprintf(_("Expected %s but got %s"),
                                 gretl_type_get_name(spec),
                                 gretl_type_get_name(rhs));
        } else {
            gretl_errmsg_sprintf(_("Specified type %s does not match array type %s"),
                                 gretl_type_get_name(spec),
                                 gretl_type_get_name(rhs));
        }
        err = E_TYPES;
    }

#if LHDEBUG
    fprintf(stderr, "lhs_type_check: spec=%d, rhs=%d: err = %d\n", spec, rhs, err);
#endif

    return err;
}

static void *get_mod_assign_result (void *lp, GretlType ltype,
                                    NODE *r, int *size, parser *p)
{
    void *ret = NULL;
    NODE *l, *op;

    if (p->op == INC || p->op == DEC) {
        /* handle increment/decrement postfix operator */
        if (ltype == GRETL_TYPE_DOUBLE) {
            double x = *(double *) lp;

            if (!na(x)) {
                x += (p->op == INC)? 1 : -1;
                *(double *) lp = x;
            }
            ret = lp;
        } else {
            p->err = E_TYPES;
        }
        return ret; /* handled */
    }

    /* create binary tree to hold left and right */
    l = newempty();
    op = newb2(p->op, l, r);

    /* if that went OK, put the relevant type specifier
       and pointer onto the LHS sub-node
    */
    if (op == NULL || l == NULL) {
        p->err = E_ALLOC;
    } else if (ltype == GRETL_TYPE_MATRIX) {
        l->t = MAT;
        l->v.m = lp;
    } else if (ltype == GRETL_TYPE_DOUBLE) {
        l->t = NUM;
        l->v.xval = *(double *) lp;
    } else if (ltype == GRETL_TYPE_STRING) {
        l->t = STR;
        l->v.str = lp;
    } else if (ltype == GRETL_TYPE_BUNDLE) {
        l->t = BUNDLE;
        l->v.b = lp;
    } else if (ltype == GRETL_TYPE_ARRAY) {
        l->t = ARRAY;
        l->v.a = lp;
    } else if (ltype == GRETL_TYPE_LIST) {
        l->t = LIST;
        l->v.ivec = lp;
        if (p->op == B_ADD) {
            /* reinterpret '+' when appending to a list */
            op->t = B_LCAT;
        }
    } else if (ltype == GRETL_TYPE_SERIES) {
        l->t = SERIES;
        l->v.xvec = lp;
    } else {
        p->err = E_TYPES;
    }

    if (!p->err) {
        /* FIXME parser state variables? */
        int saveflags = p->flags;
        int savetarg = p->targ;
        NODE *ev;

#if LHDEBUG
        fputs("*** op tree, before ***\n", stderr);
        print_tree(op, p, 0, 0);
#endif
        p->targ = l->t;
        p->flags = P_START;
        if (saveflags & P_LISTDEF) {
            p->flags |= P_LISTDEF;
        }
        ev = eval(op, p);
#if LHDEBUG
        fputs("*** ev tree, after ***\n", stderr);
        print_tree(ev, p, 0, 0);
#endif

        if (!p->err) {
            /* get @ret off node @ev and clean up */
            if (ev->t == MAT) {
                ret = ev->v.m;
            } else if (ev->t == BUNDLE) {
                ret = ev->v.b;
            } else if (ev->t == STR) {
                ret = ev->v.str;
            } else if (ev->t == ARRAY) {
                ret = ev->v.a;
            } else if (ev->t == LIST) {
                ret = ev->v.ivec;
            } else if (ev->t == SERIES) {
                /* FIXME sample range? */
                if (size != NULL) {
                    *size = p->dset->n;
                }
                ret = ev->v.xvec;
            } else if (ev->t == NUM) {
                ret = lp;
                *(double *) lp = ev->v.xval;
            } else {
                p->err = e_types(ev);
            }
        }

        p->targ = savetarg;
        p->flags = saveflags;
        free(ev);
    }

    /* thought: if @p is reusable, should we try preserving the
       nodes allocated here?
    */

    /* trash temporary nodes */
    free(op);
    free(l);

    if (ret == NULL && !p->err) {
        p->err = E_DATA;
    }

    return ret;
}

/* ".=" : we need a scalar (possibly complex) on the RHS */

static int dot_assign_to_matrix (gretl_matrix *m, parser *p)
{
    NODE *n = p->ret;
    int err = 0;

    if (scalar_node(n)) {
        gretl_matrix_fill(m, node_get_scalar(n, p));
    } else if (cscalar_node(n)) {
        err = gretl_cmatrix_fill(m, n->v.m->z[0]);
    } else {
        err = e_types(n);
    }

    return err;
}

#define empty_rhs_ok(t) (t==GRETL_TYPE_BUNDLE || gretl_is_array_type(t))

/* Setting an object in a bundle under a given key string. We get here
   only if p->lh.expr is non-NULL.
*/

static int set_bundle_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    GretlType targ = 0;
    GretlType type = 0;
    gretl_bundle *bundle;
    void *ptr = NULL;
    char *key = NULL;
    int size = 0;
    int donate = 0;
    int err = 0;

    if (lh1->t != BUNDLE) {
        return E_DATA;
    } else if (lh2->t != STR && lh2->t != MSPEC) {
        return E_DATA;
    }

    bundle = lh1->v.b;
    key = lh2->t == STR ? lh2->v.str : lh2->v.mspec->lsel.str;

    if (bundle == NULL || key == NULL) {
        return E_DATA;
    }

#if LHDEBUG
    fprintf(stderr, "set_bundle_value: bundle = %p, key = '%s'\n",
            (void *) bundle, key);
#endif

    if ((p->flags & P_PRIV) && p->op == B_ASN && null_node(rhs)) {
        /* this is an internal "special" that implements removal
           of a bundle member via the "delete" command
        */
        return gretl_bundle_delete_data(bundle, key);
    }

    if (p->op != B_ASN) {
        /* We must have an existing bundle member under @key, and
           its type will determine the type of the result of
           inflected assignment.
        */
        GretlType ltype = 0;
        void *lp;

        lp = gretl_bundle_get_target(bundle, key, &ltype, &size, &err);
        if (!err) {
            targ = gretl_type_from_gen_type(p->targ);
            err = lhs_type_check(targ, ltype, BUNDLE);
        }
        if (p->op == B_DOTASN) {
            /* accepted only for matrices */
            if (!err) {
                if (ltype == GRETL_TYPE_MATRIX) {
                    err = dot_assign_to_matrix(lp, p);
                } else {
                    err = E_TYPES;
                }
            }
            return err; /* handled */
        }
        if (!err) {
            ptr = get_mod_assign_result(lp, ltype, rhs, &size, p);
            err = p->err;
            if (p->op == INC || p->op == DEC) {
                return err; /* handled */
            }
        }
        if (!err) {
            type = ltype;
            if (ptr != lp) {
                donate = 1; /* donate: is this always right? */
            }
        }
        goto push_data;
    }

    /* Note: @targ is the gretl type specified by the caller for
       the bundle member (if any, this need not be supplied), and
       @type is the gretl type of the object arising on the RHS.
       It's an error if @targ is non-zero and @type does not
       agree with it -- except for the case where @targ is given
       as "series" and we get a suitable matrix on the right. As
       of 2015-10-03, when we get a request to put a series into
       a bundle we actually put in a matrix, which in fact makes it
       easier to get a series back out again.
    */

    if (p->targ == ARRAY) {
        targ = p->lh.gtype;
    } else {
        targ = gretl_type_from_gen_type(p->targ);
    }

    if (targ == GRETL_TYPE_NONE && null_node(rhs)) {
        /* at this point @targ is indeterminate, but maybe there's
           an existing member to fix its value?
        */
        void *lp = gretl_bundle_get_target(bundle, key, &targ, NULL, NULL);

        if (targ == GRETL_TYPE_ARRAY) {
            targ = gretl_array_get_type(lp);
        }
    }

    if (!err && targ == GRETL_TYPE_LIST) {
        ptr = node_get_list(rhs, p);
        err = p->err;
        if (!err) {
            type = GRETL_TYPE_LIST;
            donate = 1;
        }
    } else if (!err) {
        switch (rhs->t) {
        case NUM:
            if (targ == GRETL_TYPE_SERIES) {
                ptr = scalar_to_series(rhs, p);
                if (p->err) {
                    err = p->err;
                } else {
                    type = GRETL_TYPE_SERIES;
                    size = p->dset->n;
                    donate = 1;
                }
            } else if (targ == GRETL_TYPE_MATRIX) {
                ptr = gretl_matrix_from_scalar(rhs->v.xval);
                type = GRETL_TYPE_MATRIX;
                donate = 1;
            } else {
                ptr = &rhs->v.xval;
                type = GRETL_TYPE_DOUBLE;
            }
            break;
        case STR:
            ptr = rhs->v.str;
            type = GRETL_TYPE_STRING;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case MAT:
            if (targ == GRETL_TYPE_DOUBLE && scalar_matrix_node(rhs)) {
                ptr = &rhs->v.m->val[0];
                type = GRETL_TYPE_DOUBLE;
            } else if (targ == GRETL_TYPE_SERIES) {
                ptr = (double *) get_colvec_as_series_with_size(rhs, &size, p);
                if (!p->err) {
                    type = GRETL_TYPE_SERIES;
                } else {
                    err = E_TYPES;
                }
            } else {
                ptr = rhs->v.m;
                type = GRETL_TYPE_MATRIX;
                donate = is_tmp_node(rhs);
            }
            break;
        case SERIES:
            ptr = rhs->v.xvec;
            type = GRETL_TYPE_SERIES;
            size = p->dset->n;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case BUNDLE:
            ptr = rhs->v.b;
            type = GRETL_TYPE_BUNDLE;
            break;
        case ARRAY:
            ptr = rhs->v.a;
            /* get more specific type for comparison with
               what the user specified (if anything)
            */
            type = gretl_array_get_type(rhs->v.a);
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case LIST:
            if (targ == GRETL_TYPE_MATRIX) {
                ptr = gretl_list_to_vector(rhs->v.ivec, &p->err);
                if (!p->err) {
                    type = GRETL_TYPE_MATRIX;
                    donate = 1;
                }
            } else {
                ptr = rhs->v.ivec;
                type = GRETL_TYPE_LIST;
                donate = is_tmp_node(rhs);
            }
            break;
        case EMPTY:
            /* "null" is OK as (re-)initializer for bundle or array */
            if (empty_rhs_ok(targ)) {
                if (targ == GRETL_TYPE_BUNDLE) {
                    ptr = gretl_bundle_new();
                    if (ptr == NULL) {
                        err = E_ALLOC;
                    }
                } else {
                    ptr = gretl_array_new(targ, 0, &err);
                }
                if (!err) {
                    type = targ;
                    donate = 1;
                }
            } else {
                err = E_TYPES;
            }
            break;
        default:
            err = E_TYPES;
            break;
        }
    }

    if (!err) {
        /* check for result type-incompatible with user's spec */
        err = lhs_type_check(targ, type, BUNDLE);
    }

 push_data:

    if (!err) {
        if (gretl_is_array_type(type)) {
            /* revert to generic array type for the functions below */
            type = GRETL_TYPE_ARRAY;
        }
        if (donate) {
            /* it's OK to hand over the data pointer */
            err = gretl_bundle_donate_data(bundle, key, ptr, type, size);
            if (ptr == rhs->v.ptr) {
                rhs->v.ptr = NULL; /* avoid freeing! */
            }
        } else {
            /* the data must be copied into the bundle */
            err = gretl_bundle_set_data(bundle, key, ptr, type, size);
        }
        if (!err && type == GRETL_TYPE_MATRIX) {
            /* for use by genr_get_output_matrix() */
            p->lh.mret = ptr;
        }
    }

    return err;
}

static int set_array_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    GretlType atype = 0;
    GretlType type = 0;
    GretlType targ = 0;
    gretl_array *array = NULL;
    void *ptr = NULL;
    int idx = 0;
    int donate = 0;
    int err = 0;

    if (lh1->t != ARRAY) {
        return E_TYPES;
    }

    if (lh2->t == MSPEC) {
        /* FIXME allow a range here? */
        idx = mspec_get_array_index(lh2->v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = lh2->v.xval;
    }

    array = lh1->v.a;
    if (array == NULL) {
        return E_DATA;
    } else if (idx <= 0 || idx > gretl_array_get_length(array)) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
        return E_DATA;
    }

    atype = gretl_array_get_content_type(array);
    targ = gretl_type_from_gen_type(p->targ);
    err = lhs_type_check(targ, atype, ARRAY);

#if LHDEBUG
    fprintf(stderr, "set_array_value: array %p, idx=%d, atype=%s, targ=%s, err=%d\n",
            (void *) array, idx, gretl_type_get_name(atype),
            gretl_type_get_name(targ), err);
#endif

    idx--; /* convert index to 0-based */

    if (!err && p->op != B_ASN) {
        GretlType ltype = 0;
        void *lp;

        lp = gretl_array_get_element(array, idx, &ltype, &err);
        if (p->op == B_DOTASN) {
            if (!err) {
                if (ltype == GRETL_TYPE_MATRIX) {
                    err = dot_assign_to_matrix(lp, p);
                } else {
                    err = E_TYPES;
                }
            }
            return err; /* handled */
        }
        if (!err) {
            ptr = get_mod_assign_result(lp, ltype, rhs, NULL, p);
            err = p->err;
            if (p->op == INC || p->op == DEC) {
                return err; /* handled */
            }
        }
        if (!err) {
            type = ltype;
            if (ptr != lp) {
                donate = 1; /* donate: always right? */
            }
        }
        goto push_data;
    }

    if (!err && atype == GRETL_TYPE_LIST) {
        ptr = node_get_list(rhs, p);
        err = p->err;
        if (!err) {
            type = GRETL_TYPE_LIST;
            donate = 1;
        }
    } else if (!err) {
        switch (rhs->t) {
        case NUM:
            if (atype == GRETL_TYPE_MATRIX) {
                ptr = gretl_matrix_from_scalar(rhs->v.xval);
                type = GRETL_TYPE_MATRIX;
                donate = 1;
            }
            break;
        case STR:
            ptr = rhs->v.str;
            type = GRETL_TYPE_STRING;
            donate = !reusable(p) && is_tmp_node(rhs);
            break;
        case MAT:
            ptr = rhs->v.m;
            type = GRETL_TYPE_MATRIX;
            donate = is_tmp_node(rhs);
            break;
        case BUNDLE:
            ptr = rhs->v.b;
            type = GRETL_TYPE_BUNDLE;
            break;
        case LIST:
            ptr = rhs->v.ivec;
            type = GRETL_TYPE_LIST;
            donate = is_tmp_node(rhs);
            break;
        case ARRAY:
            ptr = rhs->v.a;
            type = GRETL_TYPE_ARRAY;
            donate = is_tmp_node(rhs);
            break;
        default:
            err = E_TYPES;
            break;
        }
    }

    if (!err && type != atype) {
        err = E_TYPES;
    }

 push_data:

    if (!err) {
        if (donate) {
            /* it's OK to hand over the data pointer */
            err = gretl_array_set_element(array, idx, ptr, type, 0);
            if (ptr == rhs->v.ptr) {
                rhs->v.ptr = NULL; /* gone! */
            }
        } else {
            /* the data must be copied into the array */
            err = gretl_array_set_element(array, idx, ptr, type, 1);
        }
        if (!err && type == GRETL_TYPE_MATRIX) {
            /* for use by genr_get_output_matrix() */
            p->lh.mret = ptr;
        }
    }

    return err;
}

/* setting member of list: only straight assignment is accepted */

static int set_list_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    int *list = NULL;
    int idx = 0, v = -1;
    int err = 0;

    if (p->op != B_ASN) {
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type"),
                             get_opstr(p->op));
        return E_TYPES;
    }

    if (lh2->t == MSPEC) {
        idx = mspec_get_array_index(lh2->v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = lh2->v.xval;
    }

    list = lh1->v.ivec;
    if (list == NULL) {
        return E_DATA;
    } else if (idx < 1 || idx > list[0]) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
        return E_DATA;
    }

#if LHDEBUG
    fprintf(stderr, "set_list_value: list = %p, idx = %d\n",
            (void *) list, idx);
#endif

    if (rhs->t == NUM) {
        v = node_get_int(rhs, p);
    } else if (rhs->t == SERIES) {
        v = rhs->vnum;
    } else {
        p->err = e_types(rhs);
    }

    if (!p->err) {
        if (v < 0 || v >= p->dset->v) {
            gretl_errmsg_set(_("Invalid list element"));
            p->err = E_DATA;
        }
    }

#if 0 /* we're not applying this check (yet?) */
    if (!p->err && gretl_function_depth() > 0) {
        if (!series_is_accessible_in_function(v, p->dset)) {
            p->err = E_DATA;
        }
    }
#endif

    if (!p->err) {
        list[idx] = v;
    }

    return p->err;
}

/* setting element of string: only straight assignment accepted */

static int set_string_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    char *s1 = NULL;
    const char *s2 = NULL;
    int idx = 0, err = 0;

    if (p->op != B_ASN) {
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type"),
                             get_opstr(p->op));
        return E_TYPES;
    } else if (rhs->t != STR) {
        return e_types(rhs);
    }

    if (lh2->t == MSPEC) {
        idx = mspec_get_array_index(lh2->v.mspec, &err);
        if (err) {
            return err;
        }
    } else {
        idx = lh2->v.xval;
    }

    /* @s1 may be modified below */
    s1 = lh1->v.str;

    if (s1 == NULL) {
        return E_DATA;
    } else if (idx < 1 || idx > g_utf8_strlen(s1, -1)) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds"), idx);
        return E_DATA;
    }

    s2 = rhs->v.str;
    if (g_utf8_strlen(s2, -1) != 1) {
        return E_INVARG;
    } else if (g_utf8_strlen(s1, -1) == strlen(s1) && strlen(s2) == 1) {
        /* simple: no multibyte characters */
        s1[idx-1] = s2[0];
    } else {
        /* handle the multibyte case */
        char *tmp = gretl_utf8_replace_char(s1, s2, idx - 1);

        if (strlen(tmp) <= strlen(s1)) {
            strcpy(s1, tmp);
            free(tmp);
        } else {
            user_var *uv = get_user_var_by_data(s1);

            if (uv != NULL) {
                p->err = user_var_replace_value(uv, tmp,
                                                GRETL_TYPE_STRING);
            } else {
                p->err = E_DATA;
                free(tmp);
            }
        }
    }

    return p->err;
}

static int set_series_obs_value (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    double **Z = p->dset->Z;
    char *label = NULL;
    double x = NADBL;
    int op = p->op;
    int v, t;

    if (lh1->t == SERIES && (lh2->t == NUM || lh2->t == MSPEC)) {
        ; /* OK */
    } else {
        return E_TYPES;
    }

    v = lh1->vnum;
    if (v <= 0 || v >= p->dset->v) {
        return E_DATA;
    } else if (object_is_const(NULL, v)) {
        return overwrite_err(p->dset->varname[v]);
    }

    if (lh2->t == MSPEC) {
        t = mspec_get_series_index(lh2->v.mspec, p);
    } else {
        t = node_get_int(lh2, p);
    }
    if (t < 1 || t > p->dset->n) {
        return E_DATA;
    }
    /* convert to 0-based */
    t--;

    if (rhs == NULL) {
        if (p->op == INC) {
            x = 1;
            op = B_ADD;
        } else if (p->op == DEC) {
            x = 1;
            op = B_SUB;
        } else {
            return E_TYPES; /* ? */
        }
    } else if (rhs->t == STR) {
        if (is_string_valued(p->dset, v)) {
            label = rhs->v.str;
        } else {
            return E_TYPES;
        }
    } else {
        x = node_get_scalar(rhs, p);
        if (p->err) {
            return p->err;
        }
    }

    if (is_string_valued(p->dset, v)) {
        if (label != NULL) {
            if (op != B_ASN) {
                p->err = E_TYPES;
            } else {
                p->err = series_set_string_val(p->dset, v, t, label);
            }
        } else {
            x = xy_calc(Z[v][t], x, op, NUM, p);
            if (!p->err) {
                p->err = string_series_assign_value(p->dset, v, t, x);
            }
        }
    } else {
        Z[v][t] = xy_calc(Z[v][t], x, op, NUM, p);
    }

    if (p->err == 0) {
        /* made a change to an element of a series */
        p->flags |= P_OBSVAL;
        set_dataset_is_changed(p->dset, 1);
    }

    return p->err;
}

/* Here we're replacing a submatrix, by either straight or
   inflected assignment.

   @lhs must be a binary node holding the target matrix
   on its L branch and a matrix subspec on its R branch.
   @rhs must hold the replacement value: either a matrix
   or a scalar (or a series standing in for a matrix).
*/

static int set_matrix_chunk (NODE *lhs, NODE *rhs, parser *p)
{
    NODE *lh1 = lhs->L;
    NODE *lh2 = lhs->R;
    gretl_matrix *m1, *m2 = NULL;
    matrix_subspec *spec;
    double rhs_x = NADBL;
    double complex rhs_z = NADBL;
    int rhs_scalar = 0;
    int rhs_cscalar = 0;
    int inflected = 0;
    int free_m2 = 0;

    if (p->op == B_HCAT || p->op == B_VCAT) {
        /* can't do these things on a submatrix */
        gretl_errmsg_sprintf(_("The operator '%s' is not valid in this context"),
                             get_opstr(p->op));
        return E_TYPES;
    } else if (lh1->t != MAT) {
        /* is this ever possible? */
        fprintf(stderr, "set_matrix_chunk: got %s, not matrix!\n",
                getsymb(lh1->t));
        return E_DATA;
    }

    /* set up the target */
    m1 = lh1->v.m;
    spec = lh2->v.mspec;
    if (m1 == NULL || spec == NULL) {
        return E_DATA;
    }

    /* check the validity of the subspec we got, and
       adjust it if need be in the light of the
       dimensions of @m.
    */
    if (!spec->checked) {
        p->err = check_matrix_subspec(spec, m1);
        if (p->err) {
            fprintf(stderr, "set_matrix_chunk: check_matrix_subspec failed\n");
            return p->err;
        }
    }

#if EDEBUG > 1
    gretl_matrix_print(m1, "m1, in set_matrix_chunk");
    fprintf(stderr, "op = '%s'\n", getsymb(p->op));
    print_mspec(spec);
    if (rhs != NULL) {
        fprintf(stderr, "rhs type %s\n", getsymb(rhs->t));
        if (rhs->t == NUM) fprintf(stderr, " value %g\n", rhs->v.xval);
    } else {
        fprintf(stderr, "rhs NULL\n");
    }
#endif

    /* Is the assignment straight or inflected?  Note that in
       this context there's no distinction between '=' and '.='
       and the latter doesn't count as inflected.
    */
    if (p->op != B_ASN && p->op != B_DOTASN) {
        inflected = 1;
    }

    if (p->op == INC || p->op == DEC) {
        /* treat as add or subtract */
        rhs_x = 1;
        rhs_z = rhs_x;
        rhs_scalar = 1;
    } else if (scalar_node(rhs)) {
        /* single value (could be 1 x 1 matrix) on RHS */
        rhs_x = (rhs->t == NUM)? rhs->v.xval: rhs->v.m->val[0];
        rhs_z = rhs_x;
        rhs_scalar = 1;
    } else if (cscalar_node(rhs)) {
        if (!m1->is_complex) {
            gretl_errmsg_set("Cannot assign complex values to a real matrix");
            p->err = E_TYPES;
        } else {
            m2 = rhs->v.m;
            rhs_z = rhs->v.m->z[0];
            rhs_cscalar = 1;
        }
    } else if (rhs->t == MAT) {
        /* not a scalar: get the RHS matrix */
        m2 = rhs->v.m;
        if (m2->is_complex && !m1->is_complex) {
            gretl_errmsg_set("Cannot assign complex values to a real matrix");
            p->err = E_TYPES;
        }
    } else if (rhs->t == SERIES) {
        /* legacy: this has long been accepted */
        m2 = series_to_matrix(rhs->v.xvec, p);
        free_m2 = 1; /* flag temporary status of @m2 */
    } else {
        p->err = E_TYPES;
    }

    if (p->err) {
        return p->err;
    }

    if (spec->ltype == SEL_ELEMENT) {
        /* assignment, plain or inflected, to a single
           element of target matrix.
        */
        int i = mspec_get_element(spec);

        if (rhs_cscalar || (rhs_scalar && m1->is_complex)) {
            if (!inflected) {
                m1->z[i] = rhs_z;
            } else {
                m1->z[i] = c_xy_calc(m1->z[i], rhs_z, p->op, p);
            }
        } else if (rhs_scalar) {
            if (!inflected) {
                m1->val[i] = rhs_x;
            } else {
                m1->val[i] = xy_calc(m1->val[i], rhs_x, p->op, MAT, p);
            }
        } else {
            /* here the RHS must be 1 x 1 */
            p->err = E_NONCONF;
        }
        return p->err; /* we're done */
    }

    if (!inflected) {
        if (rhs_cscalar) {
            return assign_scalar_to_submatrix(m1, m2, 0, spec);
        } else if (rhs_scalar) {
            return assign_scalar_to_submatrix(m1, NULL, rhs_x, spec);
        } else if (is_sel_dummy(spec->ltype)) {
            return gretl_matrix_set_part(m1, m2, 0, spec->ltype);
        }
    }

    if (inflected) {
        /* Here we're doing '+=' or some such, in which case a new
           submatrix must be calculated using the original
           submatrix @a and the newly generated matrix (or
           scalar value).
        */
        gretl_matrix *a = matrix_get_submatrix(m1, spec, 1, &p->err);

        if (!p->err) {
            if (rhs_scalar || rhs_cscalar) {
                if (a->is_complex) {
                    cmatrix_xy_calc(a, a, rhs_z, 0, p->op, p);
                } else {
                    rmatrix_xy_calc(a, a, rhs_x, 0, p->op, p);
                }
                m2 = a; /* assign computed matrix to m2 */
                free_m2 = 1;
            } else {
                gretl_matrix *b = NULL;

                p->err = real_matrix_calc(a, m2, p->op, &b);
                gretl_matrix_free(a);
                /* replace RHS m2 with computed result */
                if (free_m2) {
                    /* m2 was temp result of series conversion */
                    gretl_matrix_free(m2);
                }
                m2 = b;
                free_m2 = 1;
            }
            /* we now proceed to matrix_replace_submatrix() */
        }
    }

    if (!p->err) {
        /* Write new submatrix @m2 into place: note that we come here
           directly if none of the special conditions above are
           satisfied -- for example, if the newly generated value
           is a matrix and the task is straight assignment. Also
           check for numerical "breakage" in the replacement
           submatrix.
        */
        p->err = matrix_replace_submatrix(m1, m2, spec);
    }

    if (free_m2) {
        gretl_matrix_free(m2);
    }

    return p->err;
}

static gretl_matrix *get_corrgm_matrix (NODE *l,
                                        NODE *m,
                                        NODE *r,
                                        parser *p)
{
    int xcf = (r->t != EMPTY);
    int *list = NULL;
    gretl_matrix *A = NULL;
    int k;

    /* ensure we've got an order */
    k = node_get_int(m, p);
    if (p->err) {
        return NULL;
    }

    /* hook up list if arg1 is list */
    if (l->t == LIST) {
        list = l->v.ivec;
    }

    /* if third node is matrix, must be real col vector */
    if (r->t == MAT) {
        if (r->v.m->cols != 1) {
            p->err = E_NONCONF;
            return NULL;
        } else if (r->v.m->is_complex) {
            p->err = E_NONCONF;
            return NULL;
        }
    }

    if (!xcf) {
        /* acf/pacf */
        if (l->t == SERIES) {
            A = acf_matrix(l->v.xvec, k, p->dset, 0, &p->err);
        } else if (l->t == MAT) {
            A = multi_acf(l->v.m, NULL, NULL, k, &p->err);
        } else {
            /* it must be a list */
            A = multi_acf(NULL, list, p->dset, k, &p->err);
        }
    } else {
        /* cross-correlogram */
        void *px = NULL, *py = NULL;
        int xtype = SERIES;

        if (list != NULL) {
            px = list;
            xtype = LIST;
        } else if (l->t == MAT) {
            px = l->v.m;
            xtype = MAT;
        } else {
            px = l->v.xvec;
        }

        py = (r->t == MAT)? (void *) r->v.m : (void *) r->v.xvec;

        A = multi_xcf(px, xtype, py, r->t, p->dset, k, &p->err);
    }

    return A;
}

static gretl_matrix *get_density_matrix (NODE *t, double bws,
                                         int ctrl, parser *p)
{
    gretl_matrix *(*kdfunc1) (const double *, int, double,
                              gretlopt, int *);
    gretl_matrix *(*kdfunc2) (const gretl_matrix *, double,
                              gretlopt, int *);
    gretlopt opt = ctrl ? OPT_O : OPT_NONE;
    gretl_matrix *m = NULL;
    gretl_matrix *X = NULL;
    const double *x = NULL;
    int free_X = 0;
    int n = 0;

    kdfunc1 = NULL;
    kdfunc2 = NULL;

    if (t->t == SERIES) {
        n = sample_size(p->dset);
        x = t->v.xvec + p->dset->t1;
    } else if (t->t == LIST) {
        X = gretl_matrix_data_subset(t->v.ivec, p->dset,
                                     p->dset->t1, p->dset->t2,
                                     M_MISSING_SKIP, &p->err);
        free_X = 1;
    } else {
        /* matrix */
        if (t->v.m->is_complex) {
            p->err = E_CMPLX;
        } else {
            n = gretl_vector_get_length(t->v.m);
            if (n > 0) {
                /* vector */
                x = t->v.m->val;
            } else {
                X = t->v.m;
            }
        }
    }

    if (!p->err) {
        if (X != NULL) {
            kdfunc2 = get_plugin_function("multiple_kd_matrix");
        } else if (!p->err) {
            kdfunc1 = get_plugin_function("kernel_density_matrix");
        }
        if (kdfunc2 == NULL && kdfunc1 == NULL) {
            p->err = E_FOPEN;
        } else if (X != NULL) {
            m = (*kdfunc2)(X, bws, opt, &p->err);
        } else {
            m = (*kdfunc1)(x, n, bws, opt, &p->err);
        }
    }

    if (free_X) {
        gretl_matrix_free(X);
    }

    return m;
}

static int aggregate_discrete_check (const int *list, const DATASET *dset)
{
    int i, vi;

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (!accept_as_discrete(dset, vi, 0)) {
            gretl_errmsg_sprintf(_("The variable '%s' is not discrete"),
                                 dset->varname[vi]);
            return E_DATA;
        }
    }

    return 0;
}

static gretl_matrix *mshape_scalar (double x, int r, int c, int *err)
{
    gretl_matrix *m = gretl_matrix_alloc(r, c);

    if (m == NULL) {
        *err = E_ALLOC;
    } else {
        int i, n = r * c;

        for (i=0; i<n; i++) {
            m->val[i] = x;
        }
    }

    return m;
}

static void node_get_int_or_series (int *ip, double **vecp,
                                    NODE *n, parser *p)
{
    if (scalar_node(n)) {
        *ip = node_get_int(n, p);
    } else if (n->t == SERIES) {
        *vecp = n->v.xvec;
    } else {
        p->err = e_types(n);
    }
}

/* eval_3args_func: evaluate a built-in function that has three
   arguments. The @post_process flag is a convenience for the
   case where the function in question returns a matrix: it
   centralizes the creation or retrieval of an "aux node" of the
   right type, and attaches the computed matrix @A to it. This
   flag must be set to zero for all cases where the function does
   NOT return a matrix. Conversely, when a function DOES return
   a matrix it should assign this to @A, and leave the aux node
   business to the post-processor.
*/

static NODE *eval_3args_func (NODE *l, NODE *m, NODE *r,
                              int f, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *ret = NULL;
    gretl_matrix *A = NULL;
    int post_process = 1;

    if (f == F_MSHAPE) {
        if (l->t != MAT && l->t != NUM) {
            node_type_error(f, 1, MAT, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, NUM, m, p);
        } else if (!null_node(r) && !scalar_node(r)) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            int n, k2, k1 = node_get_int(m, p);

            if (scalar_node(r)) {
                k2 = node_get_int(r, p);
            } else if (l->t == NUM) {
                k2 = 1;
            } else {
                n = l->v.m->rows * l->v.m->cols;
                if (n % k1 == 0) {
                    k2 = n / k1;
                } else {
                    p->err = E_INVARG;
                }
            }
            if (!p->err) {
                if (l->t == NUM) {
                    A = mshape_scalar(l->v.xval, k1, k2, &p->err);
                } else {
                    A = gretl_matrix_shape(l->v.m, k1, k2, &p->err);
                }
            }
        }
    } else if (f == F_TRIMR) {
        if (l->t != MAT) {
            node_type_error(f, 1, MAT, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, NUM, m, p);
        } else if (!scalar_node(r)) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            int k1 = node_get_int(m, p);
            int k2 = node_get_int(r, p);

            if (!p->err) {
                A = gretl_matrix_trim_rows(l->v.m, k1, k2, &p->err);
            }
        }
    } else if (f == F_SVD) {
        /* note: the complex case is supported */
        gretl_matrix *lm = node_get_matrix(l, p, 0, 1);
        gretl_matrix *U = NULL;
        gretl_matrix *V = NULL;

        if (!p->err) {
            if (m->t == U_ADDR) {
                U = ptr_node_get_matrix(m, p);
            } else if (m->t != EMPTY) {
                node_type_error(f, 2, U_ADDR, m, p);
            }
        }
        if (!p->err) {
            if (r->t == U_ADDR) {
                V = ptr_node_get_matrix(r, p);
            } else if (r->t != EMPTY) {
                node_type_error(f, 3, U_ADDR, r, p);
            }
        }
        if (!p->err) {
            A = user_matrix_SVD(lm, U, V, &p->err);
        }
    } else if (f == F_QR) {
        gretl_matrix *X = node_get_matrix(l, p, 0, 1);
        gretl_matrix *R = NULL;
        gretl_matrix *P = NULL;

        if (!p->err) {
            if (m->t == U_ADDR) {
                R = ptr_node_get_matrix(m, p);
            } else if (m->t != EMPTY) {
                node_type_error(f, 2, U_ADDR, m, p);
            }
        }
        if (!p->err) {
            if (r->t == U_ADDR) {
                P = ptr_node_get_matrix(r, p);
            } else if (r->t != EMPTY) {
                node_type_error(f, 3, U_ADDR, r, p);
            }
        }
        if (!p->err) {
            if (X->is_complex) {
                if (P != NULL) {
                    A = gretl_cmatrix_QR_pivot_decomp(X, R, P, &p->err);
                } else {
                    A = gretl_cmatrix_QR_decomp(X, R, &p->err);
                }
            } else {
                A = user_matrix_QR_decomp(X, R, P, &p->err);
            }
        }
    } else if (f == F_VARSIMUL) {
        gretl_matrix *m1 = node_get_real_matrix(l, p, 0, 1);
        gretl_matrix *m2 = node_get_real_matrix(m, p, 1, 2);
        gretl_matrix *m3 = node_get_real_matrix(r, p, 2, 3);

        if (!p->err) {
            A = gretl_matrix_varsimul(m1, m2, m3, &p->err);
        }
    } else if (f == F_EIGEN) {
        gretl_matrix *lm = node_get_matrix(l, p, 0, 1);
        gretl_matrix *v1 = NULL, *v2 = NULL;

        if (l->t != MAT) {
            node_type_error(f, 1, MAT, l, p);
        } else {
            if (!null_node(m)) {
                v1 = ptr_node_get_matrix(m, p);
            }
            if (!null_node(r)) {
                v2 = ptr_node_get_matrix(r, p);
            }
        }
        if (!p->err) {
            if (lm->is_complex) {
                A = gretl_zgeev(lm, v1, v2, &p->err);
            } else {
                A = gretl_dgeev(lm, v1, v2, &p->err);
            }
        }
    } else if (f == F_EIGGEN) {
        /* legacy: get rid of this asap */
        gretl_matrix *lm = node_get_matrix(l, p, 0, 1);
        gretl_matrix *v1 = NULL, *v2 = NULL;

        if (l->t != MAT) {
            node_type_error(f, 1, MAT, l, p);
        } else {
            if (!null_node(m)) {
                v1 = ptr_node_get_matrix(m, p);
            }
            if (!null_node(r)) {
                v2 = ptr_node_get_matrix(r, p);
            }
        }
        if (!p->err) {
            if (lm->is_complex) {
                p->err = E_CMPLX;
            } else {
                A = old_eigengen(lm, v1, v2, &p->err);
            }
        }
    } else if (f == F_SCHUR) {
        gretl_matrix *Z = NULL;
        gretl_matrix *W = NULL;

        if (l->t != MAT || !l->v.m->is_complex) {
            node_type_error(f, 1, MAT, l, p);
        } else {
            if (!null_node(m)) {
                Z = ptr_node_get_matrix(m, p);
            }
            if (!null_node(r)) {
                W = ptr_node_get_matrix(r, p);
            }
        }
        if (!p->err) {
            A = gretl_zgees(l->v.m, Z, W, &p->err);
        }
    } else if (f == F_CORRGM) {
        if (l->t != SERIES && l->t != MAT && l->t != LIST) {
            node_type_error(f, 1, SERIES, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != EMPTY && r->t != SERIES && r->t != MAT) {
            node_type_error(f, 3, SERIES, r, p);
        } else {
            A = get_corrgm_matrix(l, m, r, p);
        }
    } else if (f == F_SEQ) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, NUM, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, NUM, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            double start = node_get_scalar(l, p);
            double end = node_get_scalar(m, p);
            double step = (!null_node(r))? node_get_scalar(r, p) : 1.0;

            if (!p->err) {
                A = gretl_matrix_seq(start, end, step, &p->err);
            }
        }
    } else if (f == F_STRNCMP) {
        post_process = 0;
        if (l->t != STR) {
            node_type_error(f, 1, STR, l, p);
        } else if (m->t != STR) {
            node_type_error(f, 2, STR, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                if (!null_node(r)) {
                    int len = node_get_int(r, p);

                    if (!p->err) {
                        ret->v.xval = strncmp(l->v.str, m->v.str, len);
                    }
                } else {
                    ret->v.xval = strcmp(l->v.str, m->v.str);
                }
            }
        }
    } else if (f == F_WEEKDAY || f == F_ISOWEEK) {
        int series_mode = 0;
        int basic_mode = 0;

        post_process = 0;

        if (null_node(m) && null_node(r)) {
            /* got a single argument */
            if (l->t == SERIES) {
                basic_mode = series_mode = 1;
            } else if (scalar_node(l)) {
                basic_mode = 1;
            } else {
                p->err = e_types(l);
            }
        } else if (l->t == SERIES && m->t == SERIES && r->t == SERIES) {
            series_mode = 1;
        } else if (!(scalar_node(l) && scalar_node(m) && scalar_node(r))) {
            p->err = E_TYPES;
        }

        if (!p->err && series_mode) {
            const double *mx = basic_mode ? NULL : m->v.xvec;
            const double *rx = basic_mode ? NULL : r->v.xvec;

            reset_p_aux(p, save_aux);
            ret = aux_series_node(p);
            if (ret != NULL && f == F_WEEKDAY) {
                p->err = fill_day_of_week_array(ret->v.xvec, l->v.xvec,
                                                mx, rx, p->dset);
            } else if (ret != NULL) {
                p->err = fill_isoweek_array(ret->v.xvec, l->v.xvec,
                                            mx, rx, p->dset);
            }
        } else if (!p->err) {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                int mo, day, yr = node_get_int(l, p);

                if (!p->err && basic_mode) {
                    char tstr[16];
                    int n;

                    sprintf(tstr, "%d", yr);
                    n = sscanf(tstr, "%4d%2d%2d", &yr, &mo, &day);
                    if (n != 3) {
                        p->err = E_INVARG;
                    }
                } else if (!p->err) {
                    mo = node_get_int(m, p);
                    day = node_get_int(r, p);
                }
                if (!p->err && f == F_WEEKDAY) {
                    int julian = 0;

                    if (yr < 0) {
                        yr = -yr;
                        julian = 1;
                    }
                    ret->v.xval = legacy_day_of_week(yr, mo, day, julian, &p->err);
                } else if (!p->err) {
                    ret->v.xval = iso_week_number(yr, mo, day, &p->err);
                }
            }
        }
    } else if (f == F_DAYSPAN) {
        guint32 ed1 = node_get_guint32(l, p);
        guint32 ed2 = node_get_guint32(m, p);
        int wkdays = node_get_int(r, p);

        post_process = 0;
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            ret->v.xval = day_span(ed1, ed2, wkdays, &p->err);
        }
    } else if (f == F_SMPLSPAN) {
        post_process = 0;
        if (l->t == STR && m->t == STR && r->t == NUM) {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                int pd = node_get_int(r, p);

                if (!p->err) {
                    ret->v.xval = sample_span(l->v.str, m->v.str,
                                              pd, &p->err);
                }
            }
        } else {
            p->err = E_TYPES;
        }
    } else if (f == F_KDENSITY) {
        if (l->t != SERIES && l->t != LIST && l->t != MAT) {
            node_type_error(f, 1, SERIES, l, p);
        } else if (m->t != NUM && m->t != EMPTY) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != NUM && r->t != EMPTY) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            double bws = (m->t != EMPTY)? m->v.xval : 1.0;
            int ctrl = (r->t != EMPTY)? (int) r->v.xval : 0;

            A = get_density_matrix(l, bws, ctrl, p);
        }
    } else if (f == F_MONTHLEN) {
        double *movec = NULL, *yrvec = NULL;
        int wkdays, julian = 0;
        int mo = 0, yr = 0;
        int rettype = NUM;

        post_process = 0;
        wkdays = node_get_int(r, p);
        if (!p->err && (wkdays < 5 || wkdays > 7)) {
            p->err = E_INVARG;
        }
        if (!p->err) {
            node_get_int_or_series(&mo, &movec, l, p);
            if (!p->err) {
                if (movec != NULL) {
                    rettype = SERIES;
                } else if (mo < 1 || mo > 12) {
                    p->err = E_INVARG;
                }
            }
        }
        if (!p->err) {
            node_get_int_or_series(&yr, &yrvec, m, p);
            if (!p->err) {
                if (yrvec != NULL) {
                    rettype = SERIES;
                } else if (yr < 0) {
                    yr = -yr;
                    julian = 1;
                }
            }
        }
        reset_p_aux(p, save_aux);
        if (!p->err && rettype == NUM) {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                ret->v.xval = get_days_in_month(mo, yr, wkdays, julian);
            }
        } else if (!p->err) {
            ret = aux_series_node(p);
            if (ret != NULL) {
                p->err = fill_monthlen_array(ret->v.xvec,
                                             p->dset->t1, p->dset->t2,
                                             wkdays, mo, yr,
                                             movec, yrvec,
                                             julian);
            }
         }
    } else if (f == F_SETNOTE || f == F_BRENAME) {
        post_process = 0;
        if (l->t != BUNDLE) {
            node_type_error(f, 1, BUNDLE, l, p);
        } else if (m->t != STR) {
            node_type_error(f, 2, STR, m, p);
        } else if (r->t != STR) {
            node_type_error(f, 3, STR, r, p);
        } else {
            reset_p_aux(p, save_aux);
            ret = aux_scalar_node(p);
            if (!p->err && f == F_SETNOTE) {
                ret->v.xval = gretl_bundle_set_note(l->v.b, m->v.str, r->v.str);
            } else if (!p->err) {
                p->err = gretl_bundle_rekey_data(l->v.b, m->v.str, r->v.str);
                if (!p->err) {
                    ret->v.xval = 0;
                }
            }
        }
    } else if (f == F_BWFILT) {
        gretl_matrix *tmp = NULL;

        post_process = 0;
        if (l->t != SERIES) {
            if (l->t == MAT) {
                cast_to_series(l, f, &tmp, NULL, NULL, p);
            } else {
                node_type_error(f, 1, SERIES, l, p);
            }
        } else if (m->t != NUM) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != NUM) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            ret = aux_series_node(p);
            if (!p->err) {
                p->err = butterworth_filter(l->v.xvec, ret->v.xvec, p->dset,
                                            m->v.xval, r->v.xval);
            }
        }
        if (tmp != NULL) {
            l->v.m = tmp;
        }
    } else if (f == F_MLAG) {
        gretl_matrix *m1 = node_get_real_matrix(l, p, 0, 1);
        gretl_matrix *m2 = node_get_real_matrix(m, p, 1, 2);

        if (p->err) {
            ; /* skip the rest */
        } else if (r->t != NUM && r->t != EMPTY) {
            /* optional scalar */
            node_type_error(f, 3, NUM, r, p);
        } else {
            double missval = (r->t == NUM)? r->v.xval : 0.0;

            A = gretl_matrix_lag(m1, m2, OPT_L, missval);
        }
    } else if (f == F_LRCOVAR) {
        gretl_matrix *mc = node_get_real_matrix(l, p, 0, 1);
        int d = 1; /* demean the matrix arg? */

        if (!p->err) {
            d = node_get_bool(r, p, d);
        }
        if (!p->err) {
            A = long_run_covariance(mc, d, &p->err);
        }
    } else if (f == F_EIGSOLVE) {
        gretl_matrix *m1 = node_get_real_matrix(l, p, 0, 1);
        gretl_matrix *m2 = node_get_real_matrix(m, p, 1, 2);

        if (p->err) {
            ; /* skip the rest */
        } else if (r->t != EMPTY && r->t != U_ADDR) {
            /* optional matrix-pointer */
            node_type_error(f, 3, U_ADDR, r, p);
        } else {
            gretl_matrix *V = NULL;

            if (r->t == U_ADDR) {
                V = ptr_node_get_matrix(r, p);
            }
            if (!p->err) {
                A = user_gensymm_eigenvals(m1, m2, V, &p->err);
            }
        }
    } else if (f == F_PRINCOMP) {
        if (l->t != MAT) {
            node_type_error(f, 1, MAT, l, p);
        } else if (m->t != NUM) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != EMPTY && r->t != NUM) {
            /* optional boolean */
            node_type_error(f, 3, NUM, r, p);
        } else {
            int cov = null_node(r) ? 0 : node_get_int(r, p);
            gretlopt opt = cov ? OPT_V : OPT_NONE;
            int k = node_get_int(m, p);

            if (!p->err) {
                A = gretl_matrix_pca(l->v.m, k, opt, &p->err);
            }
        }
    } else if (f == F_HALTON) {
        if (l->t != NUM) {
            node_type_error(f, 1, NUM, l, p);
        } else if (m->t != NUM) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != EMPTY && r->t != NUM) {
            /* optional offset */
            node_type_error(f, 3, NUM, r, p);
        } else {
            int offset = null_node(r) ? 10 : node_get_int(r, p);
            int rows = node_get_int(l, p);
            int cols = node_get_int(m, p);

            if (!p->err) {
                A = halton_matrix(rows, cols, offset, &p->err);
            }
        }
    } else if (f == F_IWISHART) {
        if (l->t != MAT && l->t != NUM) {
            node_type_error(f, 1, MAT, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(f, 2, NUM, m, p);
        } else if (r->t != EMPTY && r->t != NUM) {
            /* optional number of replications */
            node_type_error(f, 3, NUM, r, p);
        } else {
            gretl_matrix *S = node_get_real_matrix(l, p, 0, 0);
            int v = node_get_int(m, p);
            int N = null_node(r) ? 0 : node_get_int(r, p);

            if (!p->err) {
                if (N == 0) {
                    A = inverse_wishart_matrix(S, v, &p->err);
                } else {
                    A = inverse_wishart_sequence(S, v, N, &p->err);
                }
            }
        }
    } else if (f == F_AGGRBY) {
        DATASET *dset = p->dset;
        DATASET *mdset = NULL;
        const char *fnname = NULL;
        const double *x = NULL;
        const double *y = NULL;
        int *xlist = NULL;
        int *ylist = NULL;
        gretl_matrix *xm = NULL;
        gretl_matrix *ym = NULL;
        int matrix_version = 0;
        int data_args = 0;
        int mat_args = 0;

        /* if present, the optional r node must hold a string */
        if (!null_or_string(r)) {
            node_type_error(f, 3, STR, r, p);
        } else if (!null_node(r)) {
            fnname = r->v.str;
        }
        if (!p->err && !null_node(l)) {
            /* the [optional] l node holds @x */
            if (l->t == SERIES) {
                x = l->v.xvec; data_args++;
            } else if (l->t == LIST) {
                xlist = l->v.ivec; data_args++;
            } else if (l->t == MAT) {
                xm = l->v.m; data_args++;
                mat_args++;
            } else {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            /* the [required] m node holds @byvar */
            if (m->t == SERIES) {
                y = m->v.xvec; data_args++;
            } else if (m->t == LIST) {
                ylist = m->v.ivec; data_args++;
            } else if (m->t == MAT) {
                ym = m->v.m; data_args++;
                mat_args++;
            } else {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            if (mat_args > 0 && mat_args < data_args) {
                /* matrix and series/list arguments can't be mixed */
                p->err = E_TYPES;
            } else if (mat_args > 0) {
#if AGGR_PREFER_MATRIX
                matrix_version = 1;
#else
                /* convert to dataset (with borrowed Z) */
                dset = mdset = matrix_dset_plus_lists(xm, ym, &xlist,
                                                      &ylist, &p->err);
#endif
            }
        }
        if (!p->err && ylist != NULL) {
            p->err = aggregate_discrete_check(ylist, dset);
        }
        if (!p->err && matrix_version) {
            A = matrix_aggregate(xm, ym, fnname, &p->err);
        } else if (!p->err) {
            A = aggregate_by(x, y, xlist, ylist, fnname, dset, &p->err);
        }
        if (mdset != NULL) {
            destroy_dataset(mdset);
            free(xlist);
            free(ylist);
        }
    } else if (f == F_MWEIGHTS) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, NUM, l, p);
        } else if (m->t != MAT) {
            node_type_error(f, 2, MAT, m, p);
        } else if (!scalar_node(r) && r->t != STR) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            int length = node_get_int(l, p);
            int method = node_get_midas_method(r, p);
            gretl_matrix *wm = node_get_real_matrix(m, p, 1, 2);

            if (!p->err) {
                A = midas_weights(length, wm, method, &p->err);
            }
        }
    } else if (f == F_MGRADIENT) {
        if (!scalar_node(l)) {
            node_type_error(f, 1, NUM, l, p);
        } else if (m->t != MAT) {
            node_type_error(f, 2, MAT, m, p);
        } else if (!scalar_node(r) && r->t != STR) {
            node_type_error(f, 3, NUM, r, p);
        } else {
            int length = node_get_int(l, p);
            int method = node_get_midas_method(r, p);
            gretl_matrix *gm = node_get_real_matrix(m, p, 1, 2);

            if (!p->err) {
                A = midas_gradient(length, gm, method, &p->err);
            }
        }
    } else if (f == F_RESAMPLE) {
        int blocklen = 0, draws = 0;

        if (l->t != MAT) {
            node_type_error(f, 1, MAT, l, p);
        }
        if (!p->err && !null_node(m)) {
            blocklen = node_get_int(m, p);
        }
        if (!p->err && !null_node(r)) {
            draws = node_get_int(r, p);
        }
        if (!p->err) {
            if (blocklen != 0) {
                A = gretl_matrix_block_resample(l->v.m, blocklen, draws, &p->err);
            } else {
                A = gretl_matrix_resample(l->v.m, draws, &p->err);
            }
        }
    } else if (f == HF_REGLS) {
        post_process = 0;
        if (null_node(l) && null_node(m) && null_node(r)) {
            /* doing automatic MPI: no args needed */
            int (*regfunc) (PRN *);

            regfunc = get_plugin_function("regls_xv_mpi");
            if (regfunc == NULL) {
                p->err = E_FOPEN;
            } else {
                p->err = regfunc(p->prn);
            }
        } else {
            int (*regfunc) (const gretl_matrix *, const gretl_matrix *,
                            gretl_bundle *, PRN *);

            regfunc = get_plugin_function("gretl_regls");
            if (regfunc == NULL) {
                p->err = E_FOPEN;
            } else if (l->t != MAT || m->t != MAT) {
                p->err = E_TYPES;
            } else {
                if (r->t == U_ADDR) {
                    r = ptr_node_get_referent_node(r, p);
                }
                if (r->t != BUNDLE) {
                    p->err = E_TYPES;
                }
            }
            if (!p->err) {
                p->err = regfunc(l->v.m, m->v.m, r->v.b, p->prn);
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
            ret->v.xval = 0;
        }
    } else if (f == F_STACK) {
        int length = 0;
        int offset = 0;
        int *list = NULL;

        post_process = 0;
        ret = aux_empty_series_node(p);
        list = node_get_list(l, p);
        if (null_node(m)) {
            p->err = E_ARGS;
        } else {
            length = node_get_int(m, p);
        }
        if (!p->err && !null_node(r)) {
            offset = node_get_int(r, p);
        }
        if (!p->err) {
            p->err = build_stacked_series(&ret->v.xvec, list, length, offset,
                                          p->dset);
        }
        free(list);
    } else if (f == F_VMA) {
        if (l->t != MAT) {
            /* matrix A, required */
            node_type_error(f, 1, MAT, l, p);
        } else if (m->t != MAT && m->t != EMPTY) {
            /* matrix C, optional */
            node_type_error(f, 2, MAT, m, p);
        } else if (r->t != NUM && r->t != EMPTY) {
            /* horizon, optional */
            node_type_error(f, 3, NUM, r, p);
        } else {
            gretl_matrix *compan_top = node_get_real_matrix(l, p, 0, 1);
            int horizon = null_node(r) ? 24: node_get_int(r, p);
            int n = compan_top->rows;
            gretl_matrix *C = NULL;

            if (m->t != EMPTY) {
                C = node_get_real_matrix(m, p, 1, 2);
                if (C->rows != n || C->cols !=n) {
                    p->err = E_NONCONF;
                }
            }
            if (!p->err) {
                A = vma_rep(compan_top, C, horizon, &p->err);
            }
        }
    } else if (f == F_BCHECK) {
        gretl_bundle *lb = NULL;
        gretl_bundle *mb = NULL;
        gretl_array *reqd = NULL;

        post_process = 0;
        if ((lb = ptr_node_get_bundle(l, p)) == NULL) {
            node_type_error(f, 1, BUNDLE, l, p);
        } else if (!null_node(m) && m->t != BUNDLE) {
            node_type_error(f, 2, BUNDLE, m, p);
        } else if (!null_node(r) && r->t != ARRAY) {
            node_type_error(f, 3, ARRAY, r, p);
        } else {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            if (!null_node(m)) {
                mb = m->v.b;
            }
            if (!null_node(r)) {
                reqd = r->v.a;
            }
        }
        if (!p->err) {
            ret->v.xval = gretl_bundle_extract_args(lb, mb, reqd,
                                                    NULL, p->prn,
                                                    &p->err);
        }
    } else if (f == F_SPHCORR) {
        int mode = null_node(m) ? 0 : node_get_int(m, p);
        gretl_matrix *J = NULL;

        if (!p->err) {
            if (mode < 0 || mode > 2) {
                p->err = E_INVARG;
            } else if (l->t != MAT) {
                node_type_error(f, 1, MAT, l, p);
            } else if (!null_node(r) && r->t != U_ADDR) {
                node_type_error(f, 3, U_ADDR, r, p);
            } else if (!null_node(r)) {
                if (null_node(m)) {
                    mode = 1;
                }
                J = ptr_node_get_matrix(r, p);
            }
        }
        if (!p->err) {
            if (mode == 0) {
                A = omega_from_R(l->v.m, &p->err);
            } else {
                A = R_from_omega(l->v.m, mode == 2, J, &p->err);
            }
        }
    } else if (f == F_SDC || f == F_CDEMEAN || f == F_STDIZE) {
        gretl_matrix *X = NULL;
        int param = f == F_STDIZE ? 1 : 0;
        int skip_na = 0;

        if (l->t != MAT) {
            p->err = E_TYPES;
        } else if (f == F_SDC) {
            X = l->v.m;
        } else {
            A = gretl_matrix_copy(l->v.m);
            if (A == NULL) {
                p->err = E_ALLOC;
            }
        }
        if (!p->err && !null_node(m)) {
            param = node_get_int(m, p);
        }
        if (!p->err && !null_node(r)) {
            skip_na = node_get_bool(r, p, 0);
        }
        if (!p->err) {
            if (f == F_SDC) {
                A = gretl_matrix_column_sd(X, param, skip_na, &p->err);
            } else if (f == F_CDEMEAN && param == 0) {
                p->err = gretl_matrix_center(A, skip_na);
            } else if (f == F_CDEMEAN) {
                p->err = gretl_matrix_standardize(A, 1, skip_na);
            } else if (f == F_STDIZE && param == -1) {
                p->err = gretl_matrix_center(A, skip_na);
            } else {
                p->err = gretl_matrix_standardize(A, param, skip_na);
            }
        }
    }

    if (post_process) {
        if (!p->err) {
            reset_p_aux(p, save_aux);
            ret = aux_matrix_node(p);
            if (!p->err) {
                ret->v.m = A;
            }
        }
        if (p->err) {
            /* don't leak memory on error */
            gretl_matrix_free(A);
        }
    }

    return ret;
}

static NODE *geoplot_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (!p->err) {
        const char *mapfile = NULL;
        gretl_bundle *mapbun = NULL;
        int plv = -1;       /* payload series ID */
        double *plx = NULL; /* payload data */
        gretl_bundle *opts = NULL;

        if (null_node(l) || l->t == STR || l->t == BUNDLE) {
            /* [map-fname-or-bundle] [,series] [,options] */
            if (!null_node(l)) {
                if (l->t == STR) {
                    mapfile = l->v.str;
                } else {
                    mapbun = l->v.b;
                }
            }
            if (!null_node(m)) {
                if (m->t == SERIES) {
                    plv = m->vnum;
                    plx = m->v.xvec;
                } else if (m->t == BUNDLE) {
                    opts = m->v.b;
                } else {
                    p->err = E_INVARG;
                }
            }
            if (!p->err && !null_node(r)) {
                if (opts == NULL) {
                    opts = r->v.b;
                } else {
                    p->err = E_INVARG;
                }
            }
        } else if (l->t == SERIES) {
            /* series [,options] (map is implicit) */
            plv = l->vnum;
            plx = l->v.xvec;
            if (!null_node(m)) {
                if (m->t == BUNDLE) {
                    opts = m->v.b;
                } else {
                    p->err = E_INVARG;
                }
            }
            if (!p->err && !null_node(r)) {
                p->err = E_INVARG;
            }
        } else if (null_node(l) && null_node(m) && null_node(r)) {
            ; /* implicit map, no payload, no options */
        } else {
            p->err = E_INVARG;
        }

        if (!p->err) {
            p->err = ret->v.xval = geoplot_driver(mapfile, mapbun, plv, plx,
                                                  p->dset, opts);
        }
    }

    return ret;
}

static int scan_to_vector (NODE *n, const char *fmt,
                           const char *arg, int *err)
{
    gretl_matrix *m = NULL;
    user_var *uvar = NULL;
    char **S = NULL;
    int ns = 0;
    int nmax = 0;

    uvar = get_user_var_of_type_by_name(arg, GRETL_TYPE_MATRIX);
    if (uvar == NULL) {
        *err = E_INVARG;
        return 0;
    }

    if (gretl_array_get_type(n->v.a) != GRETL_TYPE_STRINGS) {
        *err = e_types(n);
    } else {
        S = gretl_array_get_strings(n->v.a, &ns);
        m = vector_from_strings(S, ns, fmt, &nmax, err);
    }

    if (!*err) {
        *err = user_var_replace_value(uvar, m, GRETL_TYPE_MATRIX);
    }

    return nmax;
}

/* Code to handle printf() and sprintf(). The @l node is a
   format string; the @r node holds varargs or may be NULL.
*/

static NODE *print_func_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (f == F_SPRINTF) {
        ret = aux_string_node(p);
    } else {
        ret = aux_scalar_node(p);
    }

    if (ret != NULL) {
        const char *fmt = l->v.str;
        const char *args = NULL;
        int n = 0;

        if (!null_node(r)) {
            args = r->v.str;
        }
        if (f == F_SPRINTF) {
            ret->v.str = do_sprintf_function(fmt, args, p->dset, &p->err);
        } else {
            p->err = do_printf(fmt, args, p->dset, p->prn, &n);
        }
        if (!p->err && f == F_PRINTF) {
            ret->v.xval = n;
        }
    }

    return ret;
}

/* Code to handle sscanf(). The @l node must be a string, or array of
   strings; @m holds a format string; and @r holds varargs.
*/

static NODE *sscanf_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        const char *fmt = m->v.str;
        const char *args = r->v.str;
        const char *src = NULL;
        int n = 0;

        if (l->t == ARRAY) {
            /* scanning array of strings to vector */
            n = scan_to_vector(l, fmt, r->v.str, &p->err);
        } else {
            src = l->v.str;
            p->err = do_sscanf(src, fmt, args, p->dset, &n);
        }
        if (!p->err) {
            ret->v.xval = n;
        }
    }

    return ret;
}

static int x_to_period (double x, char c, int *julian, int *err)
{
    if (julian != NULL && c == 'y') {
        if (x < 0) {
            *julian = 1;
            x = -x;
        } else {
            *julian = 0;
        }
    }

    if (na(x)) {
        /* note: error not flagged here */
        return -1;
    } else if (x < 0 || fabs(x) > INT_MAX) {
        *err = E_INVARG;
        return -1;
    } else {
        int k = x;
        int ret = x;

        if (c == 'y' && k <= 0) {
            ret = -1;
        } else if (c == 'm' && (k < 1 || k > 12)) {
            ret = -1;
        } else if (c == 'd' && (k < 1 || k > 31)) {
            ret = -1;
        }

        if (ret <= 0) {
            fprintf(stderr, "epochday: got %c = %d!\n", c, k);
            *err = E_INVARG;
        }

        return ret;
    }
}

static void fill_xymd (double *targ, double x)
{
    int rem;

    targ[0] = floor(x / 10000);
    rem = x - 10000 * targ[0];
    targ[1] = floor(rem / 100);
    targ[2] = rem - 100 * targ[1];
}

static void bad_date_message (int y, int m, int d)
{
    gretl_warnmsg_sprintf("%04d-%02d-%02d: %s", y, m, d,
                          _("non-existent date"));
}

/* epochday policy: NAs for year, month or day give an NA result;
   non-NA but inherently out-of-bounds values for y, m or d (for
   example, negative y, m > 12, d > 31) will produce an error,
   and otherwise non-existent dates produce NA.
*/

static NODE *eval_epochday (NODE *ny, NODE *nm, NODE *nd, parser *p)
{
    NODE *ret = NULL;
    double *x[3] = {NULL, NULL, NULL};
    double xymd[3];
    int ymd[3] = {-1, -1, -1};
    const char *code = "ymd";
    int basic_input = 0;
    int n_series = 0;
    int julian = 0;
    double sval;
    int i;

    if (null_node(nm) && null_node(nd)) {
        /* try for ISO 8601 basic input */
        basic_input = 1;
        if (scalar_node(ny)) {
            sval = node_get_scalar(ny, p);
            if (!p->err) {
                fill_xymd(xymd, sval);
                for (i=0; i<3 && !p->err; i++) {
                    ymd[i] = x_to_period(xymd[i], code[i], NULL, &p->err);
                }
            }
        } else if (ny->t == SERIES) {
            x[0] = ny->v.xvec;
            n_series = 1;
        } else {
            node_type_error(F_EPOCHDAY, 1, NUM, ny, p);
        }
    } else {
        NODE *nodes[3] = {ny, nm, nd};

        for (i=0; i<3 && !p->err; i++) {
            if (scalar_node(nodes[i])) {
                sval = node_get_scalar(nodes[i], p);
                if (!p->err) {
                    ymd[i] = x_to_period(sval, code[i], &julian, &p->err);
                }
            } else if (nodes[i]->t == SERIES) {
                x[i] = nodes[i]->v.xvec;
                n_series++;
            } else {
                node_type_error(F_EPOCHDAY, i+1, NUM, nodes[i], p);
            }
        }
    }

    if (!p->err) {
        double edt;
        int t, t1, t2;
        int y = ymd[0];
        int m = ymd[1];
        int d = ymd[2];

        if (n_series > 0) {
            t1 = p->dset->t1;
            t2 = p->dset->t2;
            ret = aux_series_node(p);
        } else {
            t1 = t2 = 0;
            ret = aux_scalar_node(p);
        }

        for (t=t1; t<=t2 && !p->err; t++) {
            if (basic_input) {
                if (x[0] != NULL) {
                    fill_xymd(xymd, x[0][t]);
                    y = x_to_period(xymd[0], 'y', NULL, &p->err);
                    m = x_to_period(xymd[1], 'm', NULL, &p->err);
                    d = x_to_period(xymd[2], 'd', NULL, &p->err);
                }
            } else {
                y = (x[0] == NULL)? y : x_to_period(x[0][t], 'y', &julian, &p->err);
                m = (x[1] == NULL)? m : x_to_period(x[1][t], 'm', NULL, &p->err);
                d = (x[2] == NULL)? d : x_to_period(x[2][t], 'd', NULL, &p->err);
            }
            if (p->err) {
                break;
            } else if (y < 0 || m < 0 || d < 0) {
                /* got an NA somewhere */
                edt = NADBL;
            } else {
                if (julian) {
                    edt = epoch_day_from_julian_ymd(y, m, d);
                } else {
                    edt = epoch_day_from_ymd(y, m, d);
                }
                if (edt <= 0) {
                    bad_date_message(y, m, d);
                    edt = NADBL;
                }
            }
            if (n_series > 0) {
                ret->v.xvec[t] = edt;
            } else {
                ret->v.xval = edt;
            }
        }
    }

    return ret;
}

/* Bessel function handler: the @r node can be of scalar, series or
   matrix type. But the @m node (order) must hold a scalar.
*/

static NODE *eval_bessel_func (NODE *l, NODE *m, NODE *r, parser *p)
{
    char ftype;
    double v;
    NODE *ret = NULL;

    if (!starting(p) && r->t != SERIES) {
        return aux_any_node(p);
    }

    ftype = l->v.str[0];
    v = node_get_scalar(m, p);

    if (r->t == NUM) {
        double x = r->v.xval;

        ret = aux_scalar_node(p);
        if (ret != NULL) {
            ret->v.xval = gretl_bessel(ftype, v, x, &p->err);
        }
    } else if (r->t == MAT) {
        const gretl_matrix *x = r->v.m;
        int i, n = x->rows * x->cols;

        ret = aux_sized_matrix_node(p, x->rows, x->cols, 0);
        if (ret != NULL) {
            for (i=0; i<n && !p->err; i++) {
                ret->v.m->val[i] = gretl_bessel(ftype, v, x->val[i], &p->err);
            }
        }
    } else if (r->t == SERIES) {
        const double *x = r->v.xvec;
        int t1 = autoreg(p) ? p->obs : p->dset->t1;
        int t2 = autoreg(p) ? p->obs : p->dset->t2;
        int t;

        ret = aux_series_node(p);
        if (ret != NULL) {
            for (t=t1; t<=t2 && !p->err; t++) {
                ret->v.xvec[t] = gretl_bessel(ftype, v, x[t], &p->err);
            }
        }
    }

    return ret;
}

/* replace_value: non-interactive search-and-replace for series and
   matrices.  @src holds the series or matrix of which we want a
   modified copy; @n0 holds the value (or vector of values) to be
   replaced; and @n1 holds the replacement value(s). It would be nice
   to extend this to lists.
*/

static NODE *replace_value (NODE *src, NODE *n0, NODE *n1, parser *p)
{
    gretl_vector *vx0 = NULL;
    gretl_vector *vx1 = NULL;
    double x0 = 0, x1 = 0;
    int k0 = -1, k1 = -1;
    NODE *ret = NULL;

    if (!starting(p)) {
        return aux_any_node(p);
    }

    /* n0: the original value, to be replaced */
    if (n0->t == NUM) {
        x0 = n0->v.xval;
    } else if (n0->t == MAT) {
        vx0 = n0->v.m;
        if (gretl_is_null_matrix(vx0)) {
            p->err = E_DATA;
        } else if ((k0 = gretl_vector_get_length(vx0)) == 0) {
            p->err = E_NONCONF;
        }
    } else {
        node_type_error(F_REPLACE, 1, NUM, n0, p);
    }

    if (p->err) {
        return NULL;
    }

    /* n1: the replacement value */
    if (n1->t == NUM) {
        x1 = n1->v.xval;
    } else if (n1->t == MAT) {
        vx1 = n1->v.m;
        if (gretl_is_null_matrix(vx1)) {
            p->err = E_DATA;
        } else if ((k1 = gretl_vector_get_length(vx1)) == 0) {
            p->err = E_NONCONF;
        }
    } else {
        node_type_error(F_REPLACE, 2, NUM, n1, p);
    }

    if (!p->err) {
        if (n0->t == NUM && n1->t == MAT) {
            /* can't replace scalar with vector */
            p->err = E_TYPES;
        } else if (k0 > 0 && k1 > 0 && k0 != k1) {
            /* if they're both vectors, they must be
               the same length */
            p->err = E_NONCONF;
        }
    }

    if (!p->err) {
        if (src->t == SERIES) {
            ret = aux_series_node(p);
        } else if (src->t == MAT) {
            ret = aux_matrix_node(p);
        } else {
            node_type_error(F_REPLACE, 3, SERIES, src, p);
        }
    }

    if (!p->err) {
        double *px0 = (vx0 != NULL)? vx0->val : &x0;
        double *px1 = (vx1 != NULL)? vx1->val : &x1;
        gretl_matrix *m = NULL;
        const double *x = NULL;
        double *targ = NULL;
        int n = 0;

        if (k0 < 0) k0 = 1;
        if (k1 < 0) k1 = 1;

        if (src->t == SERIES) {
            n = sample_size(p->dset);
            x = src->v.xvec + p->dset->t1;    /* source array */
            targ = ret->v.xvec + p->dset->t1; /* target array */
        } else if (src->t == MAT) {
            m = src->v.m;
            ret->v.m = gretl_matrix_copy(m);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
            } else {
                n = m->rows * m->cols;
                x = m->val;           /* source array */
                targ = ret->v.m->val; /* target array */
            }
        }

        if (!p->err) {
            substitute_values(targ, x, n, px0, k0, px1, k1);
        }
    }

    return ret;
}

static NODE *isoconv_node (NODE *t, NODE *n, parser *p)
{
    NODE *e, *ret = NULL;
    const double *x = NULL;
    double *ymd[3] = {NULL, NULL, NULL};
    int i, k = n->v.bn.n_nodes;

    if (p->dset == NULL) {
        p->err = E_NODATA;
    } else if (k < 3 || k > 4) {
        n_args_error(k, 3, 4, t->t, p);
    }

    if (p->err) {
        return NULL;
    }

    for (i=0; i<k && !p->err; i++) {
        e = n->v.bn.n[i];
        if (i == 0) {
            if (e->t == SERIES) {
                x = e->v.xvec + p->dset->t1;
            } else {
                node_type_error(t->t, 1, SERIES, e, p);
            }
        } else if (i == 3 && null_node(e)) {
            ; /* OK for the last arg to be omitted */
        } else if (e->t != U_ADDR) {
            node_type_error(t->t, i+1, U_ADDR, e, p);
        } else {
            e = ptr_node_get_referent_node(e, p);
            if (e == NULL || e->t != SERIES) {
                node_type_error(t->t, i+1, SERIES, e, p);
            } else {
                ymd[i-1] = p->dset->Z[e->vnum] + p->dset->t1;
            }
        }
    }

    if (!p->err) {
        ret = aux_scalar_node(p);
    }

    if (!p->err) {
        int n = sample_size(p->dset);

        ret->v.xval = iso_basic_to_extended(x, ymd[0], ymd[1], ymd[2], n);
    }

    return ret;
}

/* The arguments here are:

   @A: array to which an element is to be added
   @n: node holding candidate array element
*/

static int check_array_element_type (gretl_array *A, NODE *n)
{
    GretlType t = gretl_array_get_type(A);
    int ok = 0;

    if (t == GRETL_TYPE_ANY) {
        /* The array type is not yet determinate; this will be
           the case when when we're looking at the first element.
           If the type n->t is acceptable we use it to set the
           type of @A.
        */
        t = 0;
        if (n->t == MAT || n->t == NUM) {
            t = GRETL_TYPE_MATRICES;
        } else if (n->t == STR) {
            t = GRETL_TYPE_STRINGS;
        } else if (n->t == BUNDLE) {
            t = GRETL_TYPE_BUNDLES;
        } else if (n->t == LIST) {
            t = GRETL_TYPE_LISTS;
        } else if (n->t == ARRAY) {
            t = GRETL_TYPE_ARRAYS;
        }
        if (t > 0) {
            gretl_array_set_type(A, t);
            ok = 1;
        }
    } else {
        /* We're looking for a match between the array type
           and the type of the candidate element.
        */
        if (t == GRETL_TYPE_MATRICES) {
            ok = (n->t == MAT || n->t == NUM);
        } else if (t == GRETL_TYPE_STRINGS) {
            ok = n->t == STR;
        } else if (t == GRETL_TYPE_BUNDLES) {
            ok = n->t == BUNDLE;
        } else if (t == GRETL_TYPE_LISTS) {
            ok = n->t == LIST;
        } else if (t == GRETL_TYPE_ARRAYS) {
            ok = n->t == ARRAY;
        }
    }

    return ok ? 0 : E_TYPES;
}

static void node_nullify_ptr (NODE *n)
{
    if      (n->t == MAT)    n->v.m = NULL;
    else if (n->t == STR)    n->v.str = NULL;
    else if (n->t == BUNDLE) n->v.b = NULL;
    else if (n->t == LIST)   n->v.ivec = NULL;
    else if (n->t == ARRAY)  n->v.a = NULL;
    else if (n->t == SERIES) n->v.xvec = NULL;
}

/* 2024-11-05: in defbundle() and _() we've been turning a
   1 x 1 matrix argument into a scalar bundle member. It
   seems this is actually not a good idea so now we don't
   do it -- but to re-establish that behavior, define
   BUNDLED_MATRIX_SPECIAL to non-zero.
*/
#define BUNDLED_MATRIX_SPECIAL 0

/* supports retrieval of data for candidate array elements
   or bundle members
*/

static void *node_get_ptr (NODE *n, int f, parser *p, int *donate)
{
    void *ptr = NULL;
    int t = n->t;

    /* default to copying the node's data */
    *donate = 0;

    if (f == F_DEFBUNDLE || f == F_DEFARGS) {
        /* specific to bundles */
        if (t == ARRAY) {
            ptr = n->v.a;
        } else if (t == SERIES) {
            ptr = n->v.xvec;
        } else if (t == NUM) {
            ptr = &n->v.xval;
        }
#if BUNDLED_MATRIX_SPECIAL
        else if (scalar_matrix_node(n)) {
            ptr = n->v.m->val;
            t = NUM;
        }
#endif
    }

    if (ptr == NULL) {
        /* common to arrays and bundles */
        if (t == MAT) {
            ptr = n->v.m;
        } else if (t == STR) {
            ptr = n->v.str;
        } else if (t == BUNDLE) {
            ptr = n->v.b;
        } else if (t == LIST) {
            ptr = n->v.ivec;
        } else if (t == ARRAY) {
            ptr = n->v.a;
        }
    }

    if (t == NUM) {
        *donate = 1;
    } else if (!reusable(p) && is_tmp_node(n)) {
        *donate = 1;
        node_nullify_ptr(n);
    }

    return ptr;
}

/* given an FARGS node, detect if the first argument
   is a pointer to bundle */

static int bundle_pointer_arg0 (NODE *t)
{
    NODE *n = t->L;

    if (n->v.bn.n_nodes > 0) {
        NODE *n0 = n->v.bn.n[0];

        if (n0->t == U_ADDR && ubundle_node(n0->L)) {
            return 1;
        }
    }

    return 0;
}

struct fn_argspec {
    int f;
    int argc_min;
    int argc_max;
};

/* Check the number of arguments, @k, supplied for the multi-arg
   function @f.
*/

static int check_argc (int f, int k, parser *p)
{
    struct fn_argspec fn_argspecs[] = {
        { F_BKFILT,    1, 4 },
        { F_FILTER,    1, 5 },
        { F_MCOVG,     4, 4 },
        { F_MRLS,      4, 6 },
        { F_LOESS,     2, 6 },
        { F_GHK,       4, 5 },
        { F_QUADTAB,   1, 4 },
        { F_QLRPVAL,   4, 4 },
        { F_BOOTCI,    1, 5 },
        { F_BOOTPVAL,  1, 3 },
        { F_NADARWAT,  2, 5 },
        { F_HYP2F1,    4, 4 },
        { F_CHOWLIN,   2, 3 },
        { F_MIDASMULT, 1, 3 },
        { F_TDISAGG,   3, 5 },
        { F_COMMUTE,   2, 5 },
        { F_TOEPSOLV,  3, 4 },
        { F_RGBMIX,    3, 4 }
    };

    int argc_min = 2;
    int argc_max = 4;
    int i;

    for (i=0; i<G_N_ELEMENTS(fn_argspecs); i++) {
        if (fn_argspecs[i].f == f) {
            argc_min = fn_argspecs[i].argc_min;
            argc_max = fn_argspecs[i].argc_max;
            break;
        }
    }

    if (k < argc_min || k > argc_max) {
        n_args_error(k, argc_min, argc_max, f, p);
    }

    return p->err;
}

static int no_ts_data (parser *p)
{
    if (!dataset_is_time_series(p->dset)) {
        gretl_errmsg_set(_("This function requires time-series data"));
        p->err = E_DATA;
    }

    return p->err;
}

#define nargs_needs_ts(f) (f == F_BKFILT)

/* evaluate a built-in function that has more than three arguments */

static NODE *eval_nargs_func (NODE *t, NODE *n, parser *p)
{
    NODE *e, *ret = NULL;
    int i, k = n->v.bn.n_nodes;

    check_argc(t->t, k, p);
    if (p->err) {
        return NULL;
    }

    if (nargs_needs_ts(t->t) && no_ts_data(p)) {
        return NULL;
    }

    if (t->t == F_BKFILT) {
        const double *x = NULL;
        int bk[3] = {0};

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                if (e->t != SERIES) {
                    node_type_error(t->t, 1, SERIES, e, p);
                } else {
                    x = e->v.xvec;
                }
            } else if (null_node(e)) {
                ; /* NULL arguments are OK */
            } else {
                bk[i-1] = node_get_int(e, p);
            }
        }
        if (!p->err) {
            ret = aux_series_node(p);
        }
        if (!p->err) {
            p->err = bkbp_filter(x, ret->v.xvec, p->dset, bk[0], bk[1], bk[2]);
        }
    } else if (t->t == F_FILTER) {
        const double *x = NULL; /* series */
        gretl_matrix *X = NULL;
        gretl_matrix *C = NULL;
        gretl_matrix *A = NULL;
        gretl_matrix *x0 = NULL;
        double y0 = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                /* the series or matrix to filter */
                if (e->t != SERIES && e->t != MAT) {
                   node_type_error(t->t, i+1, 0, e, p);
                } else if (e->t == SERIES && no_ts_data(p)) {
                    ; /* p->err is set */
                } else if (e->t == SERIES) {
                    x = e->v.xvec;
                } else {
                    X = e->v.m;
                }
            } else if (i == 1) {
                /* matrix for MA polynomial (but we'll take a scalar) */
                if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
                    node_type_error(t->t, i+1, MAT, e, p);
                } else if (e->t != EMPTY) {
                    C = node_get_real_matrix(e, p, 0, 2);
                }
            } else if (i == 2) {
                /* matrix for AR polynomial (but we'll take a scalar) */
                if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
                    node_type_error(t->t, i+1, MAT, e, p);
                } else if (e->t != EMPTY) {
                    A = node_get_real_matrix(e, p, 1, 3);
                }
            } else if (i == 3) {
                /* initial (optional scalar) value for output series */
                if (e->t != EMPTY && !scalar_node(e)) {
                    node_type_error(t->t, i+1, NUM, e, p);
                } else if (e->t != EMPTY) {
                    y0 = node_get_scalar(e, p);
                    if (!p->err && na(y0)) {
                        p->err = E_MISSDATA;
                    }
                }
            } else if (i == 4) {
                /* optional pre-sample x vector */
                if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
                    node_type_error(t->t, i+1, MAT, e, p);
                } else if (e->t != EMPTY) {
                    x0 = node_get_real_matrix(e, p, 1, 5);
                }
            }
        }
        if (!p->err) {
            if (X != NULL) {
                /* matrix output wanted */
                ret = aux_matrix_node(p);
                if (!p->err) {
                    ret->v.m = filter_matrix(X, A, C, y0, x0, &p->err);
                }
            } else if (x != NULL) {
                /* series output */
                ret = aux_series_node(p);
                if (!p->err) {
                    p->err = filter_series(x, ret->v.xvec, p->dset,
                                           A, C, y0, x0);
                }
            }
        }
    } else if (t->t == F_MCOVG) {
        gretl_matrix *X = NULL;
        gretl_vector *u = NULL;
        gretl_vector *w = NULL;
        int targ, maxlag = 0;

        for (i=0; i<k && !p->err; i++) {
            targ = (i == 3)? NUM : MAT;
            e = n->v.bn.n[i];
            if ((i == 1 || i == 2) && null_node(e)) {
                ; /* for u or w, NULL is acceptable */
            } else if (e->t != targ) {
                node_type_error(t->t, i+1, targ, e, p);
            } else if (i == 0) {
                X = mat_node_get_real_matrix(e, p);
            } else if (i == 1) {
                u = mat_node_get_real_matrix(e, p);
            } else if (i == 2) {
                w = mat_node_get_real_matrix(e, p);
            } else if (i == 3) {
                maxlag = e->v.xval;
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = gretl_matrix_covariogram(X, u, w, maxlag, &p->err);
        }
    } else if (t->t == F_MOLS || t->t == F_MPOLS) {
        gretlopt opt = (t->t == F_MPOLS)? OPT_M : OPT_NONE;
        gretl_matrix *M[2] = {NULL};
        gretl_matrix *U = NULL;
        gretl_matrix *V = NULL;
        char freemat[2] = {0};

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 2) {
                if (e->t == SERIES) {
                    M[i] = tmp_matrix_from_series(e, p);
                    freemat[i] = 1;
                } else {
                    M[i] = node_get_real_matrix(e, p, i, i+1);
                }
            } else {
                if (null_node(e)) {
                    ; /* OK */
                } else if (e->t != U_ADDR) {
                    node_type_error(t->t, i+1, U_ADDR, e, p);
                } else if (i == 2) {
                    U = ptr_node_get_matrix(e, p);
                } else {
                    V = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = user_matrix_ols(M[0], M[1], U, V, opt, &p->err);
        }
        if (freemat[0]) gretl_matrix_free(M[0]);
        if (freemat[1]) gretl_matrix_free(M[1]);
    } else if (t->t == F_MRLS) {
        gretl_matrix *M[4] = {NULL};
        gretl_matrix *U = NULL;
        gretl_matrix *V = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 4) {
                M[i] = node_get_real_matrix(e, p, i, i+1);
            } else {
                if (null_node(e)) {
                    ; /* OK */
                } else if (e->t != U_ADDR) {
                    node_type_error(t->t, i+1, U_ADDR, e, p);
                } else if (i == 4) {
                    U = ptr_node_get_matrix(e, p);
                } else {
                    V = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = user_matrix_rls(M[0], M[1], M[2], M[3], U, V, &p->err);
        }
    } else if (t->t == F_NRMAX) {
        gretl_matrix *b = NULL;
        const char *sf = NULL;
        const char *sg = NULL;
        const char *sh = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                b = mat_node_get_real_matrix(e, p);
            } else if (i == 1) {
                if (e->t != STR) {
                    node_type_error(t->t, i+1, STR, e, p);
                } else {
                    sf = e->v.str;
                }
            } else if (null_node(e)) {
                ; /* OK */
            } else if (e->t != STR) {
                node_type_error(t->t, i+1, STR, e, p);
            } else if (i == 2) {
                sg = e->v.str;
            } else {
                sh = e->v.str;
            }
        }
        if (!p->err) {
            if (!gretl_vector_get_length(b)) {
                p->err = E_TYPES;
            } else if (!is_function_call(sf) ||
                       (sg != NULL && !is_function_call(sg)) ||
                       (sh != NULL && !is_function_call(sh))) {
                p->err = E_TYPES;
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            int minimize = alias_reversed(t) ? 1 : 0;

            ret->v.xval = user_NR(b, sf, sg, sh, p->dset,
                                  minimize, p->prn, &p->err);
        }
    } else if (t->t == F_LOESS) {
        const double *y = NULL, *x = NULL;
        double bandwidth = 0.5;
        int poly_order = 1;
        gretlopt opt = OPT_NONE;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 2) {
                if (e->t != SERIES) {
                    node_type_error(t->t, i+1, SERIES, e, p);
                } else if (i == 0) {
                    y = e->v.xvec;
                } else {
                    x = e->v.xvec;
                }
            } else if (i == 2 || i == 3) {
                if (e->t != NUM && e->t != EMPTY) {
                    node_type_error(t->t, i+1, NUM, e, p);
                } else if (i == 2 && e->t == NUM) {
                    poly_order = node_get_int(e, p);
                } else if (e->t == NUM) {
                    bandwidth = e->v.xval;
                }
            } else {
                if (e->t != EMPTY && e->t != NUM) {
                    node_type_error(t->t, i+1, NUM, e, p);
                } else {
                    int ival = node_get_int(e, p);

                    if (!p->err && ival != 0) {
                        if (i == 4) {
                            opt |= OPT_R;
                        } else {
                            opt |= OPT_O;
                        }
                    }
                }
            }
        }
        if (!p->err) {
            ret = aux_series_node(p);
            if (ret != NULL) {
                p->err = gretl_loess(y, x, poly_order, bandwidth,
                                     opt, p->dset, ret->v.xvec);
            }
        }
    } else if (t->t == F_GHK) {
        gretl_matrix *M[4] = {NULL};
        gretl_matrix *dP = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 4) {
                M[i] = node_get_real_matrix(e, p, i, i+1);
            } else {
                /* the optional last argument */
                if (null_node(e)) {
                    ; /* OK */
                } else if (e->t != U_ADDR) {
                    node_type_error(t->t, i+1, U_ADDR, e, p);
                } else {
                    dP = ptr_node_get_matrix(e, p);
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            if (dP == NULL) {
                ret->v.m = gretl_GHK(M[0], M[1], M[2], M[3], &p->err);
            } else {
                ret->v.m = user_matrix_GHK(M[0], M[1], M[2], M[3],
                                           dP, &p->err);
            }
        }
    } else if (t->t == F_QUADTAB) {
        int order = -1, method = 1;
        double a = NADBL;
        double b = NADBL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                order = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t->t, i+1, NUM, e, p);
            } else if (i == 1) {
                method = node_get_int(e, p);
            } else if (i == 2) {
                a = node_get_scalar(e, p);
            } else {
                b = node_get_scalar(e, p);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = gretl_quadrule_matrix_new(order, method,
                                                 a, b, &p->err);
        }
    } else if (t->t == F_IRF) {
        int targ = 0, shock = 0;
        double alpha = 0.0;
        gretl_bundle *vb = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                targ = node_get_int(e, p);
            } else if (i == 1) {
                shock = node_get_int(e, p);
            } else if (i == 2) {
                /* optional bootstrap alpha */
                if (e->t != EMPTY) {
                    alpha = node_get_scalar(e, p);
                }
            } else {
                /* final optional arg must be bundle */
                if (e->t != EMPTY && e->t != BUNDLE) {
                    node_type_error(t->t, 4, BUNDLE, e, p);
                } else if (e->t == BUNDLE) {
                    vb = e->v.b;
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            /* convert indices to zero-based */
            targ--;
            shock--;
            if (vb != NULL) {
                ret->v.m = gretl_IRF_from_bundle(vb, targ, shock, alpha,
                                                 p->dset, &p->err);
            } else {
                ret->v.m = last_model_get_irf_matrix(targ, shock, alpha,
                                                     p->dset, &p->err);
            }
        }
    } else if (t->t == F_QLRPVAL) {
        double X2 = NADBL;
        double p1 = 0, p2 = 0;
        int df = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                X2 = node_get_scalar(e, p);
            } else if (i == 1) {
                df = node_get_int(e, p);
            } else if (i == 2) {
                p1 = node_get_scalar(e, p);
            } else {
                p2 = node_get_scalar(e, p);
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            ret->v.xval = QLR_pval(X2, df, p1, p2);
        }
    } else if (t->t == F_BOOTCI) {
        int cnum = -1, method = 1, B = 0;
        int studentize = 1;
        double alpha = NADBL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                cnum = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t->t, i+1, NUM, e, p);
            } else if (i == 1) {
                B = node_get_int(e, p);
            } else if (i == 2) {
                alpha = node_get_scalar(e, p);
            } else if (i == 3) {
                method = node_get_int(e, p);
            } else {
                studentize = node_get_int(e, p);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = last_model_get_boot_ci(cnum, p->dset, B, alpha, method,
                                              studentize, &p->err);
        }
    } else if (t->t == F_BOOTPVAL) {
        int cnum = -1, method = 1, B = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                cnum = node_get_int(e, p);
            } else if (!null_or_scalar(e)) {
                node_type_error(t->t, i+1, NUM, e, p);
            } else if (i == 1) {
                B = node_get_int(e, p);
            } else {
                method = node_get_int(e, p);
            }
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            ret->v.xval = last_model_get_boot_pval(cnum, p->dset, B,
                                                   method, &p->err);
        }
    } else if (t->t == F_MOVAVG) {
        const double *x = NULL;
        double d = 0, y0 = NADBL;
        int len = 0, ctrl = -9999;
        int EMA = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i > 1 && null_node(e)) {
                continue; /* OK */
             } else if (i == 0) {
                if (e->t == SERIES) {
                    x = e->v.xvec;
                } else {
                    node_type_error(t->t, i+1, SERIES, e, p);
                }
            } else if (i == 1) {
                d = node_get_scalar(e, p);
                if (d < 1.0 && d > 0.0) {
                    EMA = 1;
                } else if (d < 1.0) {
                    p->err = E_INVARG;
                } else {
                    len = node_get_int(e, p);
                }
            } else if (i == 2) {
                ctrl = node_get_int(e, p);
            } else {
                y0 = node_get_scalar(e, p);
            }
        }
        if (!p->err) {
            ret = aux_series_node(p);
        }
        if (!p->err) {
            if (ctrl == -9999) {
                /* set the respective defaults */
                ctrl = EMA ? 1 : 0;
            }
            if (EMA) {
                p->err = exponential_movavg_series(x, ret->v.xvec, p->dset,
                                                   d, ctrl, y0);
            } else {
                p->err = movavg_series(x, ret->v.xvec, p->dset, len, ctrl);
            }
        }
    } else if (t->t == F_NADARWAT) {
        const double *x = NULL;
        const double *y = NULL;
        int LOO = 0;
        double h = 0;
        double trim = NADBL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 2 && e->t != SERIES) {
                node_type_error(t->t, i+1, SERIES, e, p);
            }
            if (p->err) {
                break;
            }
            if (i == 0) {
                x = e->v.xvec;
            } else if (i == 1) {
                y = e->v.xvec;
            } else if (i == 2) {
                /* set bandwidth? */
                if (null_node(e)) {
                    ; /* OK: will use the default */
                } else {
                    h = node_get_scalar(e, p);
                    if (h < 0 && k > 3) {
                        gretl_errmsg_sprintf(_("Bandwidth must be non-negative"));
                        p->err = E_INVARG;
                        break;
                    } else if (h < 0) {
                        /* it's a legacy thing */
                        LOO = 1;
                        h = -h;
                    }
                }
            } else if (i == 3) {
                /* Leave-one-out */
                LOO = node_get_bool(e, p, 0);
            } else if (i == 4) {
                /* trim spec? (overrides legacy "set" value) */
                trim = node_get_scalar(e, p);
            }
        }
        if (!p->err) {
            ret = aux_series_node(p);
        }
        if (!p->err) {
            if (na(trim)) {
                trim = libset_get_double(NADARWAT_TRIM);
            }
            p->err = nadaraya_watson(x, y, h, p->dset, LOO,
                                     trim, ret->v.xvec);
        }
    } else if (t->t == F_HYP2F1) {
        gretl_matrix *x = NULL;
        double a[3];

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 3) {
                a[i] = node_get_scalar(e, p);
            } else {
                x = node_get_real_matrix(e, p, 0, 3);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = gretl_matrix_alloc(x->rows, x->cols);
            if (ret->v.m == NULL) {
                p->err = E_ALLOC;
            } else {
                int n = x->rows * x->cols;
                double xi, y;

                for (i=0; i<n; i++) {
                    xi = x->val[i];
                    if (xi < -1.0 || xi > 1.0) {
                        y = NADBL;
                    } else {
                        y = hypergeo(a[0], a[1], a[2], xi);
                    }
                    ret->v.m->val[i] = y;
                }
            }
        }
    } else if (t->t == F_CHOWLIN) {
        gretl_matrix *Y = NULL;
        gretl_matrix *X = NULL;
        int fac = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                /* Y matrix */
                if (e->t == MAT) {
                    Y = e->v.m;
                } else {
                    p->err = e_types(e);
                }
            } else if (i == 1) {
                /* expansion factor */
                fac = node_get_int(e, p);
            } else if (i == 2) {
                /* optional X matrix  */
                if (e->t == MAT) {
                    X = e->v.m;
                } else if (!null_node(e)) {
                    p->err = e_types(e);
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = matrix_chowlin(Y, X, fac, &p->err);
        }
    } else if (t->t == F_MIDASMULT) {
        gretl_bundle *mb = NULL;
        int cum = 0;
        int idx = 1;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                if (e->t == BUNDLE) {
                    mb = e->v.b;
                } else {
                    p->err = e_types(e);
                }
            } else if (i == 1) {
                /* cumulate? */
                cum = node_get_bool(e, p, 0);
            } else if (!null_node(e)) {
                /* index of MIDAS term? */
                idx = node_get_int(e, p);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = midas_multipliers(mb, cum, idx, &p->err);
        }
    } else if (t->t == F_TDISAGG) {
        struct tdisagg_info tdi = {0};
        gretl_bundle *b = NULL;
        gretl_bundle *r = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                /* Y: matrix, series or list */
                if (e->t == MAT) {
                    tdi.Y = e->v.m;
                } else if (e->t == SERIES) {
                    tdi.ynum = e->vnum;
                    tdi.yval = e->v.xvec;
                } else if (e->t == LIST) {
                    tdi.ylist = e->v.ivec;
                } else {
                    p->err = e_types(e);
                }
            } else if (i == 1) {
                /* X: matrix, series, list or null */
                if (e->t == MAT) {
                    tdi.X = e->v.m;
                } else if (e->t == SERIES) {
                    tdi.xnum = e->vnum;
                    tdi.xval = e->v.xvec;
                } else if (e->t == LIST) {
                    if (gretl_is_midas_list(e->v.ivec, p->dset)) {
                        tdi.xlist = e->v.ivec;
                        tdi.xmidas = 1;
                    } else {
                        tdi.xlist = e->v.ivec;
                    }
                } else if (!null_node(e)) {
                    p->err = e_types(e);
                }
            } else if (i == 2) {
                /* integer expansion factor */
                tdi.efac = node_get_int(e, p);
            } else if (i == 3) {
                /* optional options bundle */
                if (e->t == BUNDLE) {
                    b = e->v.b;
                } else if (!null_node(e)) {
                    p->err = e_types(e);
                }
            } else if (e->t == BUNDLE) {
                /* optional retrieval bundle */
                r = e->v.b;
            } else if (!null_node(e)) {
                p->err = e_types(e);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            tdi.targ = gretl_type_from_gen_type(p->targ);
            ret->v.m = get_tdisagg_matrix(&tdi, p->dset, b, r,
                                          p->prn, &p->err);
        }
    } else if (t->t == F_COMMUTE) {
        int rowdim = 0, coldim = 0;
        int post = 0, add_id = 0;
        gretl_matrix *A = NULL;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                if (e->t == MAT) {
                    A = e->v.m;
                } else {
                    p->err = e_types(e);
                }
            } else if (i == 1) {
                /* row dimension */
                rowdim = node_get_int(e, p);
                coldim = rowdim;
            } else if (i == 2 && !null_node(e)) {
                /* column dimension */
                coldim = node_get_int(e, p);
            } else if (i == 3) {
                /* postmultiply instead of premultiply? */
                post = node_get_bool(e, p, 0);
            } else if (i == 4) {
                /* add identity matrix? */
                add_id = node_get_bool(e, p, 0);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = gretl_matrix_commute(A, rowdim, coldim, !post, add_id, &p->err);
        }
    } else if (t->t == F_TOEPSOLV) {
        gretl_matrix *m1 = node_get_real_matrix(n->v.bn.n[0], p, 0, 1);
        gretl_matrix *m2 = node_get_real_matrix(n->v.bn.n[1], p, 1, 2);
        gretl_matrix *m3 = node_get_real_matrix(n->v.bn.n[2], p, 2, 3);
        user_var *uv = NULL;

        if (k == 4) {
            uv = ptr_node_get_uvar(n->v.bn.n[3], NUM, p);
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            double d, *pdet = (uv != NULL)? &d : NULL;

            ret->v.m = gretl_toeplitz_solve(m1, m2, m3, pdet, &p->err);
            if (!p->err && uv != NULL) {
                user_var_set_scalar_value(uv, d);
            }
        }
    } else if (t->t == F_RGBMIX) {
        guint32 c[2] = {0};
        double *f = NULL;
        int do_plot = 0;
        int nf = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i < 2) {
                /* color specs */
                if (e->t == NUM) {
                    c[i] = node_get_guint32(e, p);
                } else if (e->t == STR) {
                    c[i] = numeric_color_from_string(e->v.str, &p->err);
                } else {
                    p->err = E_TYPES;
                }
            } else if (i == 2) {
                /* mixing spec(s) */
                if (e->t == MAT) {
                    nf = gretl_vector_get_length(e->v.m);
                    if (nf == 0) {
                        p->err = E_INVARG;
                    } else {
                        f = e->v.m->val;
                    }
                } else {
                    p->err = E_TYPES;
                }
            } else {
                /* plot option */
                do_plot = node_get_bool(e, p, 0);
            }
        }
        if (!p->err) {
            ret = aux_array_node(p);
        }
        if (!p->err) {
            ret->v.a = colormix_array(c[0], c[1], f, nf, do_plot, &p->err);
        }
    } else if (t->t == HF_FELOGITR) {
        gretl_matrix *U = NULL;
        gretl_matrix *X = NULL;
        int t = 0, T = 0;

        for (i=0; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 0) {
                /* number of ones for unit i */
                t = node_get_int(e, p);
            } else if (i == 1) {
                /* number of observations for unit i */
                T = node_get_int(e, p);
            } else if (i == 2) {
                /* T x 1 index function */
                if (e->t == MAT) {
                    U = e->v.m;
                } else {
                    p->err = e_types(e);
                }
            } else if (i == 3) {
                /* T x k matrix of covariates */
                if (e->t == MAT) {
                    X = e->v.m;
                } else {
                    p->err = e_types(e);
                }
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
            ret->v.m = felogit_rec_loglik(t, T, U, X);
        }
    }

    return ret;
}

/* evaluate an object definition function (bundle, array, list) */

static NODE *object_def_node (NODE *t, NODE *n, parser *p)
{
    NODE *e, *ret = NULL;
    int i, f = t->t;
    int k = n->v.bn.n_nodes;

    if ((f == F_DEFBUNDLE || f == F_DEFARGS) && k % 2 != 0) {
        gretl_errmsg_sprintf(_("%s: requires an even number of arguments"),
                             getsymb(f));
        p->err = E_INVARG;
        return NULL;
    }

    if (f == F_DEFARRAY) {
        gretl_array *A = gretl_array_new(GRETL_TYPE_ANY, 0, &p->err);
        int donate;
        void *ptr;

        for (i=0; i<k && !p->err; i++) {
            donate = 0;
            e = n->v.bn.n[i];
            p->err = check_array_element_type(A, e);
            if (p->err) {
                break;
            }
            if (e->t == NUM) {
                ptr = gretl_matrix_from_scalar(e->v.xval);
                donate = 1;
            } else {
                ptr = node_get_ptr(e, t->t, p, &donate);
            }
            if (donate) {
                /* copy not required */
                p->err = gretl_array_append_object(A, ptr, 0);
            } else {
                p->err = gretl_array_append_object(A, ptr, 1);
            }
        }
        if (p->err) {
            gretl_array_destroy(A);
        } else {
            ret = aux_array_node(p);
            if (ret != NULL) {
                ret->v.a = A;
            }
        }
    } else if (f == F_DEFBUNDLE || f == F_DEFARGS) {
        gretl_bundle *b = gretl_bundle_new();
        GretlType gtype;
        char *key = NULL;
        void *ptr;
        int donate;

        if (b == NULL) {
            p->err = E_ALLOC;
        }
        for (i=0; i<k && !p->err; i++) {
            ptr = NULL;
            e = n->v.bn.n[i];
            if (i == 0 || i % 2 == 0) {
                /* we need a key string */
                if (e->t == STR) {
                    key = e->v.str;
                } else {
                    p->err = e_types(e);
                }
            } else {
                /* we need some valid content */
                gtype = gretl_type_from_gen_type(e->t);
                if (type_can_be_bundled(gtype)) {
                    int size = 0;

                    if (e->t == SERIES) {
                        size = p->dset->n;
                    }
#if BUNDLED_MATRIX_SPECIAL
                    else if (scalar_matrix_node(e)) {
                        gtype = GRETL_TYPE_DOUBLE;
                    }
#endif
                    ptr = node_get_ptr(e, t->t, p, &donate);
                    if (donate) {
                        gretl_bundle_donate_data(b, key, ptr, gtype, size);
                    } else {
                        gretl_bundle_set_data(b, key, ptr, gtype, size);
                    }
                } else {
                    p->err = e_types(e);
                }
            }
        }
        if (p->err) {
            gretl_bundle_destroy(b);
        } else {
            ret = aux_bundle_node(p);
            if (ret != NULL) {
                ret->v.b = b;
            }
        }
    } else if (f == F_DEFLIST) {
        int *li, *full_list = gretl_list_new(0);

        for (i=0; i<k && !p->err; i++) {
            li = NULL;
            e = n->v.bn.n[i];
            if (k == 1 && e->t == ARRAY) {
                li = list_from_strings_array(e->v.a, p);
            } else if (ok_list_node_plus(e)) {
                li = node_get_list(e, p);
            } else if (e->t == MAT) {
                li = node_get_list(e, p);
            } else {
                p->err = e_types(e);
            }
            if (!p->err && li[0] > 0) {
                p->err = gretl_list_append_list(&full_list, li);
            }
            free(li);
        }
        if (p->err) {
            free(full_list);
        } else {
            ret = aux_list_node(p);
            if (ret != NULL) {
                ret->v.ivec = full_list;
            }
        }
    }

    return ret;
}

/* helper for fevalb() in case a bundled series is found */

static NODE *fevalb_get_bundled_series (gretl_bundle *b,
                                        const char *key,
                                        parser *p)
{
    NODE *ret = NULL;
    double *x = NULL;
    int n, s;

    s = sample_size(p->dset);
    if (s == 0) {
        no_data_error(p);
        return NULL;
    }

    x = gretl_bundle_get_series(b, key, &n, NULL);
    if (x == NULL) {
        p->err = E_DATA;
    } else if (n != p->dset->n && n != s) {
        p->err = E_NONCONF;
    }

    if (!p->err) {
        int v = current_series_index(p->dset, key);
        size_t sz = s * sizeof *x;

        if (v < 0) {
            /* we need a new series */
            p->err = dataset_add_NA_series(p->dset, 1);
            if (!p->err) {
                p->err = rename_series(p->dset, p->dset->v - 1,
                                       key, OPT_NONE);
            }
        }
        if (!p->err) {
            memcpy(p->dset->Z[v] + p->dset->t1, x, sz);
            ret = new_node(SERIES);
            ret->vnum = v;
            ret->v.xvec = p->dset->Z[v];
            ret->vname = gretl_strdup(key);
        }
    }

    return ret;
}

/* helper function for fevalb() */

static NODE **multi_node_from_bundle (gretl_bundle *b, int *pargc,
                                      parser *p)
{
    gretl_array *keys = NULL;
    gretl_array *args = NULL;
    NODE **nn = NULL;
    int argc;
    int i;

    keys = gretl_bundle_get_keys(b, &p->err);
    if (keys == NULL) {
        return NULL;
    }
    argc = gretl_array_get_length(keys);

    /* check for ordering array */
    if (gretl_bundle_has_key(b, "arglist")) {
        args = gretl_bundle_get_array(b, "arglist", &p->err);
        if (!p->err && !arglist_validate(keys, args)) {
            p->err = E_INVARG;
        }
        if (p->err) {
            gretl_array_destroy(keys);
            return NULL;
        }
        argc -= 1;
    }

    nn = calloc(argc, sizeof *nn);

    if (nn == NULL) {
        p->err = E_ALLOC;
    } else {
        NODE *save_aux = p->aux;
        char **S = NULL;
        GretlType t;
        NODE bn = {0};
        NODE sn = {0};
        int ns = 0;

        bn.t = BUNDLE;
        bn.v.b = b;
        sn.t = STR;
        p->aux = NULL;

        if (args != NULL) {
            S = gretl_array_get_strings(args, &ns);
        } else {
            S = gretl_array_get_strings(keys, &ns);
            strings_array_sort(&S, &ns, OPT_NONE);
        }

        p->aux = NULL;
        for (i=0; i<ns; i++) {
            t = gretl_bundle_get_member_type(b, S[i], NULL);
            if (t == GRETL_TYPE_SERIES) {
                nn[i] = fevalb_get_bundled_series(b, S[i], p);
            } else {
                sn.v.str = S[i];
                nn[i] = get_bundle_member(&bn, &sn, p);
            }
            if (p->err) {
                break;
            }
        }
        p->aux = save_aux;
    }

    if (!p->err) {
        *pargc = argc;
    } else if (nn != NULL) {
        for (i=0; i<argc; i++) {
            if (nn[i] != NULL) {
                free_node(nn[i], p);
            }
        }
        free(nn);
        nn = NULL;
    }

    gretl_array_destroy(keys);

    return nn;
}

/* In the "traditional" feval() case @l is a multi-node holding the
   function name plus arguments and @r is not referenced.

   In the bundlized fevalb() case, @l holds the name of the function
   to be called and @r holds a bundle of arguments,
*/

static NODE *eval_feval (int f, NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    NODE **nn = NULL;
    char *fname = NULL;
    ufunc *u = NULL;
    int argc = 0;
    int fid = 0;

    if (f == F_FEVALB) {
        gretl_bundle *argb = r->v.b;

        fname = l->v.str;
        nn = multi_node_from_bundle(argb, &argc, p);
    } else {
        int k = l->v.bn.n_nodes;
        NODE *e = NULL;

        if (k < 1) {
            p->err = E_ARGS;
        } else {
            e = l->v.bn.n[0];
            if (e->t != STR) {
                node_type_error(f, 1, STR, e, p);
            }
        }
        if (!p->err) {
            argc = k - 1;
            fname = e->v.str;
            nn = l->v.bn.n + 1; /* array of arg nodes */
        }
    }

    if (p->err) {
        return NULL;
    }

    /* first try for a built-in function */
    fid = function_lookup(fname);
    if (fid != 0) {
        NODE tmp = {0};
        NODE aux = {0};
        NODE mn = {0};
        int np;

        /* temporary function node */
        tmp.t = fid;
        if (fid < FP_MAX) {
            tmp.v.ptr = get_genr_function_pointer(fid);
        }

        np = func1_symb(fid) ? 1 : func2_symb(fid) ? 2 :
            func3_symb(fid) ? 3 : -1;
        aux.t = EMPTY;

        if (np > 0) {
            /* known max number of arguments */
            if (argc > np) {
                gretl_errmsg_sprintf("%s: too many arguments", fname);
                p->err = E_DATA;
            } else if (np == 1) {
                tmp.L = argc > 0 ? nn[0] : &aux;
            } else if (np == 2) {
                tmp.L = argc > 0 ? nn[0] : &aux;
                tmp.R = argc > 1 ? nn[1] : &aux;
            } else if (np == 3) {
                tmp.L = argc > 0 ? nn[0] : &aux;
                tmp.M = argc > 1 ? nn[1] : &aux;
                tmp.R = argc > 2 ? nn[2] : &aux;
            }
        } else {
            /* multi-arg function */
            mn.t = FARGS;
            mn.v.bn.n_nodes = argc;
            mn.v.bn.n = nn;
            tmp.L = &mn;
        }
        if (!p->err) {
            ret = eval(&tmp, p);
        }
        if (tmp.aux != NULL) {
            /* transfer guardianship to the parent of @l */
            if (l->parent->aux != NULL) {
                free_node(l->parent->aux, p);
            }
            l->parent->aux = tmp.aux;
            tmp.aux = NULL;
        }
    }

    if (!p->err && fid == 0) {
        /* try for a user function */
        u = get_user_function_by_name(fname);
        if (u != NULL) {
            NODE tmp = {0};
            NODE mn = {0};

            tmp.t = UFUN;
            tmp.vname = fname;
            tmp.v.ptr = u;
            mn.v.bn.n_nodes = argc;
            mn.v.bn.n = nn;
            tmp.R = &mn;
            ret = eval_ufunc(&tmp, &mn, p);
        }
    }

    if (!p->err && fid == 0 && u == NULL) {
        gretl_errmsg_sprintf(_("%s: function not found"), fname);
        p->err = E_DATA;
    }

    if (f == F_FEVALB) {
        int i;

        for (i=0; i<argc; i++) {
            free_node(nn[i], p);
        }
        free(nn);
    }

    return ret;
}

/* try to get a matrix from @n, even if it's not in fact a
   matrix node as such, provided we can make a matrix out
   of its content
*/

static gretl_matrix *node_get_matrix_lenient (NODE *n,
                                              int ok,
                                              parser *p)
{
    gretl_matrix *m = NULL;

    if (n->t == NUM && (ok == NUM)) {
        m = gretl_matrix_from_scalar(n->v.xval);
    } else if (n->t == SERIES && (ok == SERIES)) {
        m = gretl_vector_from_series(n->v.xvec, p->dset->t1,
                                     p->dset->t2);
    } else if (n->t == LIST && (ok == SERIES)) {
        m = gretl_matrix_data_subset(n->v.ivec, p->dset,
                                     p->dset->t1, p->dset->t2,
                                     M_MISSING_OK, &p->err);
    } else {
        p->err = e_types(n);
    }

    if (!p->err && m == NULL) {
        p->err = E_ALLOC;
    }

    return m;
}

static gretl_bundle *get_kalman_bundle_arg (NODE *n, parser *p)
{
    gretl_bundle *b = NULL;
    NODE *e = NULL;

    if (n->t == FARGS) {
        /* multi-arguments node */
        e = n->v.bn.n[0];
        e = e->L;
    } else if (n->t == U_ADDR) {
        e = ptr_node_get_referent_node(n, p);
    } else if (n->t == BUNDLE) {
        e = n;
    }

    if (e == NULL || e->t != BUNDLE) {
        p->err = e_types(e);
    } else {
        b = e->v.b; /* get the actual bundle */
        if (gretl_bundle_get_type(b) != BUNDLE_KALMAN ||
            gretl_bundle_get_private_data(b) == NULL) {
            p->err = E_TYPES;
            b = NULL;
        }
    }

    if (p->err) {
        gretl_errmsg_set("Argument 1 must point to a state-space bundle");
    }

    return b;
}

static NODE *eval_kalman_bundle_func (NODE *t, NODE *n, parser *p)
{
    NODE *e, *ret = NULL;
    int i, k = n->v.bn.n_nodes;

    if (t->t == F_KSETUP) {
        gretl_matrix *M[5] = {NULL};
        int copy[5] = {0};
        int dkvar = 0;
        int nmats = 0;

        if (k < 4) {
            n_args_error(k, 4, 4, t->t, p);
        } else if (k >= 5) {
            nmats = 5;
            if (k > 5) {
                dkvar = node_get_bool(n->v.bn.n[5], p, 0);
            }
        } else {
            nmats = 4;
        }
        for (i=0; i<nmats && !p->err; i++) {
            e = n->v.bn.n[i];
            if (e->t == MAT) {
                M[i] = mat_node_get_real_matrix(e, p);
                if (!p->err) {
                    if (is_tmp_node(e)) {
                        e->v.m = NULL;
                    } else {
                        copy[i] = 1;
                    }
                }
            } else if (i == 0) {
                /* obsy: accept series or list */
                M[i] = node_get_matrix_lenient(e, SERIES, p);
            } else {
                /* system matrices, state variance */
                M[i] = node_get_matrix_lenient(e, NUM, p);
            }
        }
        if (!p->err) {
            gretl_bundle *b = kalman_bundle_new(M, copy, nmats, dkvar, &p->err);

            if (!p->err) {
                ret = aux_bundle_node(p);
                if (ret != NULL) {
                    ret->v.b = b;
                }
            }
        }
    } else if (t->t == F_KFILTER) {
        gretl_bundle *b = get_kalman_bundle_arg(n, p);

        if (!p->err && k != 1) {
            n_args_error(k, 1, 1, t->t, p);
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
        }
        if (!p->err) {
            ret->v.xval = kalman_bundle_filter(b, p->prn, &p->err);
        }
    } else if (t->t == F_KDSMOOTH) {
        gretl_bundle *b = get_kalman_bundle_arg(n, p);
        int param = 1;
        int dkstyle = 0;

        if (!p->err) {
            if (k == 2) {
                e = n->v.bn.n[1];
                dkstyle = node_get_int(e, p);
            } else if (k < 1 || k > 2) {
                n_args_error(k, 1, 2, t->t, p);
            }
        }
        if (!p->err) {
            param += dkstyle != 0;
            ret = aux_scalar_node(p);
            if (!p->err) {
                ret->v.xval = kalman_bundle_smooth(b, param, p->prn);
            }
        }
    } else if (t->t == F_KSMOOTH) {
        gretl_bundle *b = get_kalman_bundle_arg(n, p);

        if (!p->err && k != 1) {
            n_args_error(k, 1, 1, t->t, p);
        }
        if (!p->err) {
            ret = aux_scalar_node(p);
            if (!p->err) {
                ret->v.xval = kalman_bundle_smooth(b, 0, p->prn);
            }
        }
    } else if (t->t == F_KSIMUL) {
        /* we need a bundle pointer, a matrix,
           and perhaps an optional boolean
        */
        gretl_bundle *b = get_kalman_bundle_arg(n, p);
        gretl_matrix *U = NULL;
        int freeU = 0;
        int get_state = 0;

        if (!p->err && k != 2 && k != 3) {
            n_args_error(k, 2, 3, t->t, p);
        }
        for (i=1; i<k && !p->err; i++) {
            e = n->v.bn.n[i];
            if (i == 1) {
                if (e->t == MAT) {
                    U = mat_node_get_real_matrix(e, p);
                } else {
                    U = node_get_matrix_lenient(e, SERIES, p);
                    if (U != NULL) {
                        freeU = 1;
                    }
                }
            } else {
                get_state = node_get_int(e, p);
            }
        }
        if (!p->err) {
            ret = aux_matrix_node(p);
        }
        if (!p->err) {
            ret->v.m = kalman_bundle_simulate(b, U, get_state,
                                              p->prn, &p->err);
        }
        if (freeU) {
            gretl_matrix_free(U);
        }
    }

    return ret;
}

static NODE *kalman_data_node (NODE *l, NODE *r, parser *p)
{
    gretl_bundle *b = get_kalman_bundle_arg(l, p);
    NODE *ret = NULL;

    if (!p->err) {
        ret = aux_matrix_node(p);
    }
    if (!p->err) {
        ret->v.m = kalman_bundle_simdata(b, r->v.m, p->prn, &p->err);
    }

    return ret;
}

/* Create a matrix using selected series, or a mixture of series and
   lists, or more than one list.  We proceed by setting up a "dummy"
   dataset and constructing a list that indexes into it.  (We can't
   use a regular list, in the general case, since some of the series
   may be temporary variables that are not part of the "real"
   dataset.)
*/

static gretl_matrix *assemble_matrix (NODE **nn, int nnodes, parser *p)
{
    NODE *n;
    gretl_matrix *m = NULL;
    const int *list;
    double **Z = NULL;
    int *dumlist;
    int i, j, k = 0;

#if EDEBUG
    fprintf(stderr, "assemble_matrix...\n");
#endif

    if (nnodes == 1 && get_matrix_mask() == NULL) {
        /* take a shortcut if we just got a single series
           and there's no "matrix mask" in place
        */
        n = nn[0];
        if (n->t == SERIES) {
            m = series_to_matrix(n->v.xvec, p);
            return m;
        }
    }

    /* how many columns will we need? */
    for (i=0; i<nnodes; i++) {
        n = nn[i];
        if (n->t == LIST) {
            k += n->v.ivec[0];
        } else if (n->t == SERIES) {
            k++;
        }
    }

    /* create dummy data array */
    Z = malloc(k * sizeof *Z);
    if (Z == NULL) {
        p->err = E_ALLOC;
        return NULL;
    }

#if EDEBUG
    fprintf(stderr, " got %d columns, Z at %p\n", k, (void *) Z);
#endif

    /* and a list associated with Z */
    dumlist = gretl_consecutive_list_new(0, k-1);
    if (dumlist == NULL) {
        p->err = E_ALLOC;
        free(Z);
        return NULL;
    }

    /* attach series pointers to Z */
    k = 0;
    for (i=0; i<nnodes; i++) {
        n = nn[i];
        if (n->t == LIST) {
            list = n->v.ivec;
            for (j=1; j<=list[0]; j++) {
                Z[k++] = p->dset->Z[list[j]];
            }
        } else if (n->t == SERIES) {
            Z[k++] = n->v.xvec;
        }
    }

    if (!p->err) {
        DATASET dumset = {0};

        dumset.Z = Z;
        dumset.v = k;
        dumset.n = p->dset->n;
        dumset.t1 = p->dset->t1;
        dumset.t2 = p->dset->t2;

        m = real_matrix_from_list(dumlist, &dumset, p);
    }

    free(dumlist);
    free(Z);

    return m;
}

#define ok_matdef_sym(s) (s == NUM || s == SERIES || s == EMPTY || \
                          s == DUM || s == LIST || s == ARRAY)

/* composing a matrix from scalars, series or lists */

static NODE *matrix_def_node (NODE *t, NODE *n, parser *p)
{
    gretl_matrix *M = NULL;
    NODE *e, *ret = NULL;
    int k = n->v.bn.n_nodes;
    int nnum = 0, nvec = 0;
    int dum = 0, nsep = 0;
    int nlist = 0;
    int seppos = -1;
    int i;

    if (autoreg(p)) {
        fprintf(stderr, "You can't define a matrix in this context\n");
        p->err = E_TYPES;
        return NULL;
    }

#if EDEBUG
    fprintf(stderr, "Processing DEFMAT...\n");
#endif

    for (i=0; i<k && !p->err; i++) {
        e = n->v.bn.n[i];
        if (!ok_matdef_sym(e->t) && !scalar_matrix_node(e)) {
            p->err = e_types(e);
            break;
        }
        if (scalar_node(e)) {
            nnum++;
        } else if (e->t == SERIES) {
            nvec++;
        } else if (e->t == DUM) {
            dum++;
        } else if (e->t == LIST) {
            nlist++;
        } else if (e->t == EMPTY) {
            if (nsep == 0) {
                seppos = i;
            }
            nsep++;
        }
        if (dum && k != 1) {
            /* dummy, array must be singleton nodes */
            p->err = E_TYPES;
        } else if ((nvec || nlist) && nnum) {
            /* can't mix series/lists with scalars */
            p->err = E_TYPES;
        } else if ((nvec || nlist) && nsep) {
            /* can't have row separators in a matrix
               composed of series or lists */
            p->err = E_TYPES;
        }
    }

    if (!p->err) {
        NODE **nn = n->v.bn.n;

        if (nvec > 0 || nlist > 1) {
            M = assemble_matrix(nn, k, p);
        } else if (nnum > 0) {
            M = matrix_from_scalars(nn, k, nsep, seppos, p);
        } else if (nlist) {
            M = matrix_from_list(nn[0], p);
        } else if (dum) {
            n = nn[0];
            if (n->v.idnum == DUM_DATASET) {
                if (gretl_function_depth() > 0) {
                    gretl_errmsg_set("'dataset' is not recognized as a list within functions");
                    p->err = E_DATA;
                } else {
                    M = matrix_from_list(NULL, p);
                }
            } else {
                pprintf(p->prn, "Wrong sort of dummy var\n");
                p->err = E_TYPES;
            }
        } else {
            /* empty matrix def */
            M = gretl_null_matrix_new();
        }
    }

    if (p->err) {
        if (M != NULL) {
            gretl_matrix_free(M);
        }
    } else {
        ret = aux_matrix_node(p);
        if (ret != NULL) {
            ret->v.m = M;
        } else {
            gretl_matrix_free(M);
        }
    }

    return ret;
}

static NODE *gen_series_from_string (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    gchar *line;
    int vnum = -1;
    int err = 0;

    line = g_strdup_printf("%s=%s", l->v.str, r->v.str);
    err = generate(line, p->dset, GRETL_TYPE_SERIES,
                   OPT_U, p->prn);

    if (!err) {
        vnum = current_series_index(p->dset, l->v.str);
    }

    ret = aux_scalar_node(p);
    if (ret != NULL) {
        ret->v.xval = vnum;
    }

    g_free(line);

    return ret;
}

static const double *xvec_from_matrix (gretl_matrix *m,
                                       parser *p,
                                       int *subset,
                                       int *err)
{
    const double *ret = NULL;

    if (gretl_is_null_matrix(m) || m->cols != 1) {
        *err = E_TYPES;
    } else if (m->rows == p->dset->n) {
        ret = m->val;
    } else if (m->rows == sample_size(p->dset)) {
        *subset = 1;
        ret = m->val;
    } else {
        *err = E_TYPES;
    }

    return ret;
}

static NODE *gen_series_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (p->dset == NULL || p->dset->n == 0) {
        no_data_error(p);
    } else if (l->t == STR && r->t == STR) {
        return gen_series_from_string(l, r, p);
    } else if (l->t != STR || (r->t != SERIES && r->t != MAT && r->t != NUM)) {
        p->err = E_TYPES;
    } else {
        char *vname = l->v.str;
        int v = current_series_index(p->dset, vname);
        const double *xvec = NULL;
        double xt = NADBL;
        int subset = 0;

        if (v == 0) {
            /* can't overwite the constant */
            p->err = E_INVARG;
        } else if (v == -1) {
            /* no series called @vname: can we create one? */
            if (user_var_get_type_by_name(vname) != GRETL_TYPE_NONE) {
                p->err = E_TYPES;
            } else {
                p->err = check_varname(vname);
            }
        }

        if (!p->err) {
            /* get the data source */
            if (r->t == SERIES) {
                xvec = r->v.xvec;
            } else if (r->t == MAT) {
                xvec = xvec_from_matrix(r->v.m, p, &subset, &p->err);
            } else {
                xt = r->v.xval;
            }
        }

        if (!p->err && v == -1) {
            /* add a new series if required */
            p->err = dataset_add_NA_series(p->dset, 1);
            if (!p->err) {
                v = p->dset->v - 1;
                strcpy(p->dset->varname[v], vname);
            }
        }

        if (!p->err) {
            /* copy data across */
            int t, i = 0;

            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                if (xvec != NULL) {
                    xt = subset ? xvec[i++] : xvec[t];
                }
                p->dset->Z[v][t] = xt;
            }
        }

        if (!p->err) {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                ret->v.xval = v;
            }
        }
    }

    return ret;
}

static NODE *gen_array_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (!null_or_scalar(n)) {
        p->err = e_types(n);
    } else if (p->lh.gtype == 0) {
        gretl_errmsg_set(_("array: no type was specified"));
        p->err = E_DATA;
    } else {
        int len = 0;

        if (!null_node(n)) {
            len = node_get_int(n, p);
        }

        if (!p->err) {
            ret = aux_array_node(p);
            if (!p->err) {
                ret->v.a = gretl_array_new(p->lh.gtype, len, &p->err);
            }
        }
    }

    return ret;
}

static NODE *get_series_stringvals (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_array_node(p);

    if (!p->err) {
        int v = l->vnum;

        if (strvals_node(l, p)) {
            int sub = node_get_bool(r, p, 0);
            int n_strs = 0;
            char **S;

            S = series_get_string_vals(p->dset, v, &n_strs, sub);
            ret->v.a = gretl_array_from_strings(S, n_strs, 1,
                                                &p->err);
        } else {
            ret->v.a = gretl_array_new(GRETL_TYPE_STRINGS, 0, &p->err);
        }
    }

    return ret;
}

static NODE *stringify_series (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
        if (f == F_STRVSORT && !null_node(r)) {
            ret->v.xval = series_reorder_strings(p->dset, l->vnum, r->v.a);
        } else if (f == F_STRVSORT) {
            ret->v.xval = series_alphabetize_strings(p->dset, l->vnum);
        } else if (null_node(r)) {
            ret->v.xval = series_destroy_string_table(p->dset, l->vnum);
        } else {
            ret->v.xval = series_set_string_vals(p->dset, l->vnum, r->v.a);
        }
        if (ret->v.xval != 0.0) {
            p->err = E_DATA;
        }
    }

    return ret;
}

enum {
    FORK_L,
    FORK_R,
    FORK_BOTH,
    FORK_NONE
};

/* Determine whether or not a series is constant in boolean terms,
   i.e. all elements zero, or all non-zero, over the relevant range.
   If so, return FORK_L (all 1) or FORK_R (all 0), otherwise
   return FORK_UNK.
*/

static int vec_branch (const double *c, parser *p)
{
    int c1, t, t1, t2;
    int ret;

    t1 = autoreg(p) ? p->obs : p->dset->t1;
    t2 = autoreg(p) ? p->obs : p->dset->t2;

    c1 = (c[t1] != 0.0);
    ret = (c1)? FORK_L : FORK_R;

    for (t=t1; t<=t2; t++) {
        if (!na(c[t])) {
            if ((c1 && c[t] == 0) || (!c1 && c[t] != 0)) {
                ret = FORK_BOTH;
                break;
            }
        }
    }

    return ret;
}

/* Given a series condition in a ternary "?" expression, return the
   evaluated counterpart.  We evaluate both forks and select based on
   the value of the condition at each observation.  We accept only
   scalar (NUM or 1x1 matrix) and series (SERIES) types on input, and
   always produce a SERIES type on output.
*/

static NODE *query_eval_series (const double *c, NODE *n, parser *p)
{
    NODE *l = NULL, *r = NULL, *ret = NULL;
    NODE *save_aux = p->aux;
    double *xvec = NULL, *yvec = NULL;
    double x = NADBL, y = NADBL;
    double xt, yt;
    int t, t1, t2;
    int branch;

    branch = vec_branch(c, p);

    if (autoreg(p) || branch != FORK_R) {
        l = eval(n->M, p);
        if (p->err) {
            return NULL;
        }
        if (l->t == SERIES) {
            xvec = l->v.xvec;
        } else if (scalar_node(l)) {
            x = node_get_scalar(l, p);
        } else {
            p->err = e_types(l);
            return NULL;
        }
    }

    if (autoreg(p) || branch != FORK_L) {
        r = eval(n->R, p);
        if (p->err) {
            return NULL;
        }
        if (r->t == SERIES) {
            yvec = r->v.xvec;
        } else if (scalar_node(r)) {
            y = node_get_scalar(r, p);
        } else {
            p->err = e_types(r);
            return NULL;
        }
    }

    reset_p_aux(p, save_aux);
    ret = aux_series_node(p);

    t1 = autoreg(p) ? p->obs : p->dset->t1;
    t2 = autoreg(p) ? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
        if (na(c[t])) {
            ret->v.xvec[t] = NADBL;
        } else {
            xt = (xvec != NULL)? xvec[t] : x;
            yt = (yvec != NULL)? yvec[t] : y;
            ret->v.xvec[t] = (c[t] != 0.0)? xt : yt;
        }
    }

    return ret;
}

/* The condition in the ternary query operator is a scalar,
   which has been evaluated to @x, and which must now be
   interpreted as a boolean.
*/

static NODE *query_eval_scalar (double x, NODE *n, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *l = NULL, *r = NULL, *ret = NULL;
    int branch, gen_t = 0;

    branch = na(x) ? FORK_NONE : (x != 0 ? FORK_L : FORK_R);

    if (autoreg(p) || branch != FORK_R) {
        l = eval(n->M, p);
        if (p->err) {
            return NULL;
        } else {
            gen_t = l->t;
        }
    }

    if (autoreg(p) || branch != FORK_L) {
        r = eval(n->R, p);
        if (p->err) {
            return NULL;
        } else {
            gen_t = r->t;
        }
    }

    if (branch == FORK_NONE) {
        if (gen_t == NUM) {
            reset_p_aux(p, save_aux);
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                ret->v.xval = NADBL;
            }
        } else {
            gretl_errmsg_set(_("indeterminate condition for '?' operator"));
            p->err = E_DATA;
        }
    } else if (branch == FORK_L) {
        ret = l;
    } else if (branch == FORK_R) {
        ret = r;
    }

    return ret;
}

/* The following allows for @n to hold a scalar, real matrix
   or complex matrix, and also allows for the case where
   @n's payload is real-valued but a complex result is required,
   signalled by @need_z.
*/

static int query_term_get_value (NODE *n, int i, int j,
                                 double *py, double complex *pz,
                                 int need_z)
{
    if (n->t == MAT) {
        if (n->v.m == NULL) {
            fprintf(stderr, "missing matrix in query_term_get_value: n->v.m = %p\n",
                    (void *) n->v.m);
            return E_DATA;
        } else if (n->v.m->is_complex) {
            *pz = gretl_cmatrix_get(n->v.m, i, j);
        } else {
            *py = gretl_matrix_get(n->v.m, i, j);
            if (need_z) {
                *pz = *py + 0 * I;
            }
        }
    } else {
        *py = n->v.xval;
        if (need_z) {
            *pz = *py + 0 * I;
        }
    }

    return 0;
}

/* the condition in the ternary query operator is a matrix */

static NODE *query_eval_matrix (gretl_matrix *m, NODE *n, parser *p)
{
    NODE *save_aux = p->aux;
    NODE *l, *r, *ret = NULL;
    gretl_matrix *mret;
    int lcomplex = 0;
    int rcomplex = 0;

    if (gretl_is_null_matrix(m)) {
        p->err = E_TYPES;
        return NULL;
    }

    l = eval(n->M, p);
    if (!p->err) {
        r = eval(n->R, p);
    }
    if (p->err) {
        return NULL;
    }

    if ((l->t != NUM && l->t != MAT) ||
        (r->t != NUM && r->t != MAT)) {
        if (gretl_matrix_is_scalar(m)) {
            double x = m->val[0];

            if (na(x)) {
                gretl_errmsg_set(_("indeterminate condition for '?' operator"));
                p->err = E_DATA;
            } else if (x != 0) {
                ret = l;
            } else {
                ret = r;
            }
        } else {
            p->err = E_TYPES;
        }
        return ret;
    }

    if (l->t == MAT) {
        if (l->v.m->cols != m->cols || l->v.m->rows != m->rows) {
            p->err = E_NONCONF;
            return NULL;
        } else if (l->v.m->is_complex) {
            lcomplex = 1;
        }
    } else if (r->t == MAT) {
        if (r->v.m->cols != m->cols || r->v.m->rows != m->rows) {
            p->err = E_NONCONF;
            return NULL;
        } else if (r->v.m->is_complex) {
            rcomplex = 1;
        }
    }

    if (lcomplex || rcomplex) {
        mret = gretl_cmatrix_build(m, NULL, 0, 0, &p->err);
    } else {
        mret = gretl_matrix_copy(m);
        if (mret == NULL) {
            p->err = E_ALLOC;
        }
    }
    if (p->err) {
        return NULL;
    }

    reset_p_aux(p, save_aux);
    ret = aux_matrix_node(p);

    if (!p->err) {
        int need_z = mret->is_complex;
        double complex z = {0};
        double x, y = 0;
        int j, i;

        for (j=0; j<m->cols; j++) {
            for (i=0; i<m->rows; i++) {
                x = gretl_matrix_get(m, i, j);
                if (isnan(x)) {
                    if (mret->is_complex) {
                        gretl_cmatrix_set(mret, i, j, x + x * I);
                    } else {
                        gretl_matrix_set(mret, i, j, x);
                    }
                } else if (x != 0.0) {
                    p->err = query_term_get_value(l, i, j, &y, &z, need_z);
                    if (p->err) {
                        break;
                    }
                    if (mret->is_complex) {
                        gretl_cmatrix_set(mret, i, j, z);
                    } else {
                        gretl_matrix_set(mret, i, j, y);
                    }
                } else {
                    p->err = query_term_get_value(r, i, j, &y, &z, need_z);
                    if (p->err) {
                        break;
                    }
                    if (mret->is_complex) {
                        gretl_cmatrix_set(mret, i, j, z);
                    } else {
                        gretl_matrix_set(mret, i, j, y);
                    }
                }
            }
        }
        ret->v.m = mret;
    }

    return ret;
}

/* Evaluate a ternary "query" expression: C ? X : Y.  The condition C
   may be a scalar, series or matrix. On entry to this function @c will
   hold the evaluated C; the alternates are represented by the M and R
   sub-nodes of @t, not yet evaluated.
*/


static NODE *eval_query (NODE *t, NODE *c, parser *p)
{
    NODE *ret = NULL;

#if EDEBUG
    if (c->t == NUM) {
        fprintf(stderr, " condition type=NUM, value=%g\n", c->v.xval);
    } else {
        fprintf(stderr, " condition type=%s\n", getsymb(c->t));
    }
#endif

    if (c->t == NUM) {
        ret = query_eval_scalar(c->v.xval, t, p);
    } else if (c->t == SERIES) {
        ret = query_eval_series(c->v.xvec, t, p);
    } else if (c->t == MAT) {
        if (0 /* gretl_matrix_is_scalar(c->v.m) */) {
            /* Watch out! This breaks compilation. */
            ret = query_eval_scalar(c->v.m->val[0], t, p);
        } else {
            ret = query_eval_matrix(c->v.m, t, p);
        }
    } else {
        /* invalid type for boolean condition */
        p->err = e_types(c);
    }

#if EDEBUG
    fprintf(stderr, "eval_query return: ret = %p, err %d\n", (void *) ret, p->err);
#endif

    return ret;
}

#define dvar_scalar(i) (i > 0 && i < R_SCALAR_MAX)
#define dvar_series(i) (i > R_SCALAR_MAX && i < R_SERIES_MAX)
#define dvar_matrix(i) (i == R_NOW)
#define dvar_variant1(i) (i == R_TEST_STAT || i == R_TEST_PVAL)
#define dvar_variant2(i) (i == R_RESULT)

#define no_data(p) (p == NULL || p->n == 0)

double dvar_get_scalar (int i, const DATASET *dset)
{
    switch (i) {
    case R_NOBS:
        return (dset == NULL) ? NADBL :
        (dset->n == 0 ? 0 : sample_size(dset));
    case R_NVARS:
        return (dset == NULL)? NADBL : dset->v;
    case R_PD:
        return (no_data(dset))? NADBL : dset->pd;
    case R_PANEL_PD:
        return (dataset_is_panel(dset))? dset->panel_pd : NADBL;
    case R_T1:
        return (no_data(dset))? NADBL : dset->t1 + 1;
    case R_T2:
        return (no_data(dset))? NADBL : dset->t2 + 1;
    case R_TMAX:
        if (no_data(dset)) {
            return NADBL;
        } else {
            int tmax;

            sample_range_get_extrema(dset, NULL, &tmax);
            return tmax + 1;
        }
    case R_DATATYPE:
        return dataset_get_structure(dset);
    case R_TEST_PVAL:
        return get_last_pvalue();
    case R_TEST_STAT:
        return get_last_test_statistic();
    case R_TEST_LNL:
        return get_last_lnl();
    case R_TEST_BRK:
        return get_last_break();
    case R_STOPWATCH:
        return gretl_stopwatch();
    case R_WINDOWS:
#ifdef WIN32
        return 1;
#else
        return 0;
#endif
    case R_VERSION:
        return gretl_version_number(GRETL_VERSION);
    case R_ERRNO:
        return get_gretl_errno();
    case R_SEED:
        return gretl_rand_get_seed();
    case R_HUGE:
        return libset_get_double(CONV_HUGE);
    case R_LOGLEVEL:
        return libset_get_int(LOGLEVEL);
    case R_LOGSTAMP:
        return libset_get_bool(LOGSTAMP);
    default:
        return NADBL;
    }
}

static int date_series_ok (const DATASET *dset)
{
    if (calendar_data(dset)) {
        return 1;
    } else if (quarterly_or_monthly(dset)) {
        return 1;
    } else if (annual_data(dset) || decennial_data(dset)) {
        return 1;
    } else if (dataset_has_panel_time(dset)) {
        return 1;
    } else {
        return 0;
    }
}

static int dvar_get_series (double *x, int i, const DATASET *dset)
{
    int tseries = dataset_is_time_series(dset);
    int panel = dataset_is_panel(dset);
    int YMD = calendar_data(dset);
    int t, err = 0;

    if (dset == NULL || dset->n == 0) {
        return E_NODATA;
    }

    if (i == R_OBSMIN && dset->pd < 2) {
        return E_PDWRONG;
    }

    if (i == R_OBSMIC && !YMD) {
        return E_PDWRONG;
    }

    if (i == R_PUNIT && !panel) {
        return E_PDWRONG;
    }

    if (i == R_TIME && !panel && !tseries) {
        return E_PDWRONG;
    }

    if (i == R_DATES && !date_series_ok(dset)) {
        return E_PDWRONG;
    }

    if (i == R_OBSMAJ) {
        if (dset->pd == 1 && !tseries) {
            i = R_INDEX;
        } else if (panel) {
            i = R_PUNIT;
        }
    }

    if (i == R_INDEX) {
        for (t=0; t<dset->n; t++) {
            x[t] = t + 1;
        }
    } else if (YMD && i != R_INDEX && i != R_DATES) {
        /* Watch out: we're handling most calendar-data cases
           here, so we have to explicitly exclude cases that
           require different treatment.
        */
        char obs[12];
        int y, m, d;

        for (t=0; t<dset->n && !err; t++) {
            ntolabel(obs, t, dset);
            if (sscanf(obs, YMD_READ_FMT, &y, &m, &d) != 3) {
                err = E_DATA;
            } else if (i == R_OBSMAJ) {
                x[t] = y;
            } else if (i == R_OBSMIN) {
                x[t] = m;
            } else {
                x[t] = d;
            }
        }
    } else if (i == R_TIME) {
        if (panel) {
            make_panel_time_var(x, dset);
        } else {
            for (t=0; t<dset->n; t++) {
                x[t] = (double) (t + 1);
            }
        }
    } else if (i == R_PUNIT) {
        for (t=0; t<dset->n; t++) {
            x[t] = t / dset->pd + 1;
        }
    } else if (i == R_OBSMAJ) {
        int maj;

        for (t=0; t<dset->n; t++) {
            date_maj_min(t, dset, &maj, NULL);
            x[t] = maj;
        }
    } else if (i == R_OBSMIN) {
        int min;

        for (t=0; t<dset->n; t++) {
            date_maj_min(t, dset, NULL, &min);
            x[t] = min;
        }
    } else if (i == R_DATES) {
        err = fill_dataset_dates_series(dset, x);
    } else {
        err = E_DATA;
    }

    return err;
}

static gretl_matrix *dvar_get_matrix (int i, int *err)
{
    gretl_matrix *m = NULL;

    switch (i) {
    case R_TEST_STAT:
        m = get_last_test_matrix(err);
        break;
    case R_TEST_PVAL:
        m = get_last_pvals_matrix(err);
        break;
    case R_NOW:
        m = gretl_matrix_alloc(1, 2);
        if (m == NULL) {
            *err = E_ALLOC;
        } else {
            /* package epoch second and ISO 8601 basic date */
            time_t t = time(NULL);
            struct tm tm;
            int y, mon, d;

#ifdef WIN32
            tm = *localtime(&t);
#else
            localtime_r(&t, &tm);
#endif
            y = tm.tm_year + 1900;
            mon = tm.tm_mon + 1;
            d = tm.tm_mday;
            m->val[0] = (double) t;
            m->val[1] = y * 10000 + mon * 100 + d;
        }
        break;
    default:
        *err = E_DATA;
        break;
    }

    return m;
}

static NODE *dollar_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        int idx = t->v.idnum;

        if (dvar_scalar(idx)) {
            ret = aux_scalar_node(p);
            if (ret != NULL) {
                ret->v.xval = dvar_get_scalar(idx, p->dset);
            }
        } else if (dvar_series(idx)) {
            ret = aux_series_node(p);
            if (ret != NULL) {
                p->err = dvar_get_series(ret->v.xvec, idx, p->dset);
            }
        } else if (dvar_matrix(idx)) {
            ret = aux_matrix_node(p);
            if (ret != NULL) {
                ret->v.m = dvar_get_matrix(idx, &p->err);
            }
        } else if (idx == R_PNGFONT) {
            ret = aux_string_node(p);
            if (ret != NULL) {
                ret->v.str = gretl_png_font_string();
            }
        } else if (idx == R_MAPFILE) {
            const char *fname = dataset_get_mapfile(p->dset);

            ret = aux_string_node(p);
            if (!p->err) {
                if (fname == NULL) {
                    ret->v.str = gretl_strdup("");
                } else {
                    ret->v.str = gretl_strdup(fname);
                }
            }
        } else if (idx == R_MAP) {
            ret = aux_bundle_node(p);
            if (!p->err) {
                ret->v.b = get_current_map(p->dset, NULL, &p->err);
            }
        } else if (dvar_variant1(idx)) {
            GretlType type = get_last_test_type();

            if (type == GRETL_TYPE_MATRIX) {
                ret = aux_matrix_node(p);
                if (ret != NULL) {
                    ret->v.m = dvar_get_matrix(idx, &p->err);
                }
            } else {
                /* scalar or none */
                ret = aux_scalar_node(p);
                if (ret != NULL) {
                    ret->v.xval = dvar_get_scalar(idx, p->dset);
                }
            }
        } else if (dvar_variant2(idx)) {
            GretlType type = 0;
            void *ptr = get_last_result_data(&type, &p->err);

            if (type == GRETL_TYPE_MATRIX) {
                ret = aux_matrix_node(p);
                if (ret != NULL) {
                    ret->v.m = ptr;
                }
            } else if (type == GRETL_TYPE_BUNDLE) {
                ret = aux_bundle_node(p);
                if (ret != NULL) {
                    ret->v.b = ptr;
                }
            } else if (!p->err) {
                p->err = E_TYPES;
            }
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *fevd_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL) {
        int targ = node_get_int(l, p);
        int shock = 0;

        if (!p->err && !null_node(m)) {
            shock = node_get_int(m, p);
        }
        if (!p->err && !null_node(r) && r->t != BUNDLE) {
            p->err = E_INVARG;
        }

        if (!p->err) {
            /* convert @targ, @shock to zero-based */
            targ -= 1;
            shock -= 1;
            if (r->t == BUNDLE) {
                ret->v.m = gretl_FEVD_from_bundle(r->v.b, targ, shock,
                                                  p->dset, &p->err);
            } else {
                GretlObjType otype;
                void *ptr;

                ptr = get_last_model(&otype);
                if (ptr == NULL || otype != GRETL_OBJ_VAR) {
                    p->err = E_BADSTAT;
                } else {
                    ret->v.m = gretl_VAR_get_FEVD_matrix(ptr, targ, shock, 0,
                                                         p->dset, &p->err);
                }
            }
        }
    }

    return ret;
}

static GretlType object_var_type (int idx, const char *oname,
                                  int *needs_data)
{
    GretlType vtype = GRETL_TYPE_NONE;

    *needs_data = 0;

    if (model_data_scalar(idx)) {
        vtype = GRETL_TYPE_DOUBLE;
    } else if (model_data_series(idx)) {
        vtype = GRETL_TYPE_SERIES;
    } else if (model_data_matrix(idx)) {
        vtype = GRETL_TYPE_MATRIX;
    } else if (model_data_matrix_builder(idx)) {
        vtype = GRETL_TYPE_MATRIX;
        *needs_data = 1;
    } else if (model_data_list(idx)) {
        vtype = GRETL_TYPE_LIST;
    } else if (model_data_string(idx)) {
        vtype = GRETL_TYPE_STRING;
    } else if (model_data_array(idx)) {
        vtype = GRETL_TYPE_ARRAY;
    } else if (idx == B_MODEL) {
        vtype = GRETL_TYPE_BUNDLE;
    } else if (idx == B_SYSTEM) {
        vtype = GRETL_TYPE_BUNDLE;
    }

    if (idx == M_UHAT || idx == M_YHAT || idx == M_SIGMA) {
        /* could be a matrix */
        vtype = saved_object_get_data_type(oname, idx);
    }

    return vtype;
}

/* For example, $coeff(sqft); or model1.$coeff(const) */

static NODE *dollar_str_node (NODE *t, MODEL *pmod, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        const char *str = NULL;
        NODE *l = t->L;
        NODE *r = t->R;

        if (r->t == STR) {
            str = r->v.str;
        } else {
            /* could be element of strings array? */
            NODE *e = eval(r, p);

            if (e->t == STR) {
                str = e->v.str;
            } else {
                p->err = e_types(e);
            }
        }
        if (!p->err) {
            ret->v.xval = gretl_model_get_data_element(pmod, l->v.idnum, str,
                                                       p->dset, &p->err);
        }
        if (p->err && r->t == STR) {
            pprintf(p->prn, _("'%s': invalid argument for %s()\n"),
                    r->v.str, mvarname(l->v.idnum));
        }
    }

    return ret;
}

/* Retrieve a data item from a model (single equation or system of
   some kind. We're not sure in advance here of the type of the data
   item (scalar, matrix, etc.): we look that up with object_var_type().

   If retrieval is from a named model, the model name will be given
   by the @vname member of @t and the key will be on node @k.
   Otherwise @k will be NULL and @t will already hold the index of
   the required data item.
*/

static NODE *model_var_node (NODE *t, NODE *k, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
        const char *oname = NULL;
        int idx, needs_data = 0;
        GretlType vtype;

        if (t->t == DBUNDLE) {
            /* pseudo-bundle: holds name of model */
            oname = t->vname;
            idx = mvar_lookup(k->v.str);
            if (idx == 0) {
                p->err = E_DATA;
                return NULL;
            }
        } else {
            /* @t already holds the data-item index */
            idx = t->v.idnum;
        }

        /* determine the type of the requested data */
        vtype = object_var_type(idx, oname, &needs_data);

#if EDEBUG
        fprintf(stderr, "model_var_node: idx = %d, vtype = %d\n", idx, vtype);
#endif

        if (vtype == GRETL_TYPE_DOUBLE) {
            ret = aux_scalar_node(p);
        } else if (vtype == GRETL_TYPE_SERIES) {
            ret = aux_series_node(p);
        } else if (vtype == GRETL_TYPE_LIST) {
            ret = aux_list_node(p);
        } else if (vtype == GRETL_TYPE_STRING) {
            ret = aux_string_node(p);
        } else if (vtype == GRETL_TYPE_BUNDLE) {
            ret = aux_bundle_node(p);
        } else if (vtype == GRETL_TYPE_ARRAY) {
            ret = aux_array_node(p);
        } else {
            ret = aux_matrix_node(p);
        }

        if (ret == NULL) {
            return ret;
        } else if (vtype == GRETL_TYPE_DOUBLE) {
            ret->v.xval = saved_object_get_scalar(oname, idx, p->dset,
                                                  &p->err);
        } else if (vtype == GRETL_TYPE_SERIES) {
            p->err = saved_object_get_series(ret->v.xvec, oname, idx, p->dset);
        } else if (vtype == GRETL_TYPE_LIST) {
            ret->v.ivec = saved_object_get_list(oname, idx, &p->err);
        } else if (vtype == GRETL_TYPE_STRING) {
            ret->v.str = saved_object_get_string(oname, idx, p->dset,
                                                 &p->err);
        } else if (vtype == GRETL_TYPE_BUNDLE) {
            ret->v.b = bundle_from_model(NULL, p->dset, &p->err);
        } else if (vtype == GRETL_TYPE_MATRIX) {
            if (needs_data) {
                ret->v.m = saved_object_build_matrix(oname, idx,
                                                     p->dset, &p->err);
            } else {
                ret->v.m = saved_object_get_matrix(oname, idx, &p->err);
            }
        } else if (vtype == GRETL_TYPE_ARRAY) {
            ret->v.a = saved_object_get_array(oname, idx, p->dset, &p->err);
        }
    } else {
        ret = aux_any_node(p);
    }

    return ret;
}

static NODE *wildlist_node (NODE *n, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *list = varname_match_list(p->dset, n->v.str,
                                       &p->err);

        ret->v.ivec = list;
    }

    return ret;
}

static NODE *ellipsis_list_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int v1 = useries_node(l) ? l->vnum : l->v.xval;
        int v2 = useries_node(r) ? r->vnum : r->v.xval;

        ret->v.ivec = ellipsis_list(p->dset, v1, v2, &p->err);
    }

    return ret;
}

static NODE *mat2list_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL) {
        const char *pfx = null_node(r) ? NULL : r->v.str;

        ret->v.ivec = list_from_matrix(l->v.m, pfx, p->dset, &p->err);
    }

    return ret;
}

/* see if a plain NUM node can be interpreted as holding a
   series ID, in the context of creating a list */

static int could_be_series_id (NODE *n, parser *p)
{
    if (n->t == NUM && p->dset != NULL && p->dset->n > 0) {
        int k = node_get_int(n, p);

        if (!p->err && k >= 0 && k < p->dset->v) {
            return 1;
        }
    }

    return 0;
}

static NODE *list_join_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
        int *L1 = node_get_list(l, p);
        int *L2 = node_get_list(r, p);

        if (!p->err) {
            ret->v.ivec = gretl_lists_join_with_separator(L1, L2);
            if (ret->v.ivec == NULL) {
                p->err = E_ALLOC;
            }
        }

        free(L1);
        free(L2);
    }

    return ret;
}

static NODE *two_scalars_func (NODE *l, NODE *r, int t, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        double xl = node_get_scalar(l, p);
        double xr = node_get_scalar(r, p);

        if (!na(xl) && !na(xr)) {
            if (t == F_MIN) {
                ret->v.xval = (xl < xr)? xl : xr;
            } else if (t == F_MAX) {
                ret->v.xval = (xl > xr)? xl : xr;
            } else if (t == F_RANDINT) {
                int k;

                p->err = gretl_rand_int_minmax(&k, 1, xl, xr);
                if (!p->err) {
                    ret->v.xval = k;
                }
            }
        }
    }

    return ret;
}

/* Return the main or max value "within" the argument attached
   to the node @n.
*/

static NODE *within_minmax_node (NODE *n, int f, parser *p)
{
    NODE *ret = NULL;

    if (n->t == NUM) {
        ret = pretend_matrix_scalar_func(n, f, p);
    } else if (n->t == SERIES) {
        ret = series_scalar_func(n, f, NULL, p);
    } else if (n->t == MAT) {
        if (n->v.m->is_complex) {
            p->err = E_CMPLX;
        } else {
            ret = series_scalar_func(n, f, NULL, p);
        }
    } else if (ok_list_node(n, p)) {
        ret = list_to_series_func(n, f, p);
    } else {
        p->err = E_INVARG;
    }

    return ret;
}

/* Return the main or max value "between" the arguments attached
   to the nodes @l and @r.
*/

static NODE *between_minmax_node (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (scalar_node(l) && scalar_node(r)) {
        ret = two_scalars_func(l, r, f, p);
    } else if (l->t != r->t) {
        p->err = E_TYPES;
    } else if (l->t == SERIES) {
        ret = aux_series_node(p);
        if (ret != NULL) {
            p->err = two_series_minmax(l->v.xvec, r->v.xvec,
                                       ret->v.xvec, f, p->dset);
        }
    } else if (l->t == MAT) {
        if (l->v.m->is_complex || r->v.m->is_complex) {
            p->err = E_CMPLX;
        } else {
            ret = aux_matrix_node(p);
        }
        if (ret != NULL) {
            ret->v.m = two_matrices_minmax(l->v.m, r->v.m,
                                           f, &p->err);
        }
    } else {
        p->err = E_INVARG;
    }

    return ret;
}

static NODE *kpss_crit_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
        int T = node_get_int(l, p);
        int trend = node_get_int(r, p);

        if (!p->err) {
            ret->v.m = kpss_critvals(T, trend, &p->err);
        }
    }

    return ret;
}

static NODE *scalar_postfix_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
        if (n->vname != NULL && p->lh.name[0] != '\0' &&
            strcmp(n->vname, p->lh.name) == 0) {
            /* undefined behavior */
            fprintf(stderr, "BAD NUM_P/NUM_M\n");
            gretl_errmsg_sprintf(_("The result for %s is not well defined"),
                                 n->vname);
            p->err = E_DATA;
            ret->v.xval = NADBL;
        } else {
            double x = n->v.xval;

            ret->v.xval = x;
            if (n->t == NUM_P) {
                p->err = node_replace_scalar(n, x + 1.0);
            } else {
                p->err = node_replace_scalar(n, x - 1.0);
            }
        }
    }

    return ret;
}

static int series_calc_nodes (NODE *l, NODE *r)
{
    int ret = 0;

    if (l->t == SERIES) {
        ret = (r->t == SERIES || r->t == NUM || scalar_matrix_node(r));
    } else if (r->t == SERIES) {
        ret = scalar_node(l);
    }

    return ret;
}

static int cast_series_to_list (parser *p, NODE *n, short f)
{
    if (p->tree->t == F_GENSERIES || p->tree->t == UFUN) {
        /* any other cases when we shouldn't do this "cast"? */
        return 0;
    } else if (p->targ == LIST && useries_node(n)) {
        return (f == F_LOG || f == F_DIFF ||
                f == F_LDIFF || f == F_SDIFF ||
                f == F_ODEV);
    } else {
        return 0;
    }
}

static NODE *lhs_terminal_node (NODE *t, NODE *l, NODE *r,
                                parser *p)
{
    /* Pass through eval'd @l and @r subnodes, but don't eval
       the parent @t itself */
    NODE *ret = aux_parent_node(p);

    ret->t = t->t; /* transcribe type */
    ret->L = l;    /* evaluated left-hand */
    ret->R = r;    /* evaluated right-hand */

    /* prevent double-freeing of children @l and @r */
    ret->flags |= LHT_NODE;

    return ret;
}

/* reattach_series: on successive executions of a given
   compiled "genr", the "xvec" pointer recorded on a
   SERIES node will have become invalid if:

   (1) the series has been renamed,
   (2) the dataset has been differently sub-sampled, or
   (3) the series has been deleted (should be impossible).

   In case (2) the ID number of the series should still
   be valid, so it ought to be sufficient to reconnect the
   xvec pointer as in the second branch below. In case (1),
   however, it's necessary to look up the ID number by name
   again (as things stand).

   Note that n->v.xvec will be NULL when this function is
   reached only in case a genr is attached to a loop that
   is saved across function calls -- then the pointer is
   reset to NULL on each call to the function (but not
   on each iteration of the loop itself).
*/

static void reattach_series (NODE *n, parser *p)
{
    /* FIXME relax the "1" below? */
    if (1 || n->v.xvec == NULL || get_loop_renaming()) {
        /* do a full reset */
        n->vnum = current_series_index(p->dset, n->vname);
        if (n->vnum < 0) {
            gretl_errmsg_sprintf("'%s': not a series", n->vname);
            p->err = E_DATA;
        } else {
            n->v.xvec = p->dset->Z[n->vnum];
        }
    } else if (n->vnum >= 0 && n->vnum < p->dset->v) {
        /* handles the case of the dataset moving due to
           sub-sampling, provided that no series are deleted
           (the name -> ID mapping remains unchanged)
        */
        n->v.xvec = p->dset->Z[n->vnum];
    }
}

static void reattach_data_error (NODE *n, parser *p)
{
    GString *msg = g_string_new(NULL);

    g_string_printf(msg, "'%s': expected %s", n->vname, getsymb(n->t));

    if (n->uv == NULL) {
        g_string_append(msg, " but name look-up failed");
        p->err = E_DATA;
    } else {
        const char *s = gretl_type_get_name(n->uv->type);

        if (n->uv->ptr == NULL) {
            g_string_append_printf(msg, ", found %s with NULL data", s);
            p->err = E_DATA;
        } else {
            g_string_append_printf(msg, " but found %s", s);
            p->err = E_TYPES;
        }
    }

    gretl_errmsg_set(msg->str);
    g_string_free(msg, TRUE);
}

/* Deep debugging to detect any cases where we've failed to NULL-ify
   uservars (n->uv) that may have been relocated in memory between
   visits to geneval. Turn this on in case of "weird" genr errors.
*/
#define UVDEBUG 0

static void node_reattach_data (NODE *n, parser *p)
{
    if (n->t == SERIES) {
        reattach_series(n, p);
    } else {
#if UVDEBUG
        void *p0 = n->uv;
#endif
        GretlType type = 0;
        void *data = NULL;

#if UVDEBUG
        if (n->uv == NULL || gretl_iterating() || gretl_function_recursing()) {
            n->uv = get_user_var_by_name(n->vname);
            if (p0 != NULL && n->uv != p0) {
                fprintf(stderr, "node_reattach_data ('%s'): MOVED %p -> %p (ptr %p)\n",
                        n->vname, p0, (void *) n->uv, (void *) n->uv->ptr);
            }
        }
#else
        if (n->uv == NULL || user_qsorting || (n->t == LIST && gretl_iterating())) {
            n->uv = get_user_var_by_name(n->vname);
        }
#endif

        if (n->uv != NULL) {
            data = n->uv->ptr;
            type = n->uv->type;
        }

        if (data == NULL) {
            reattach_data_error(n, p);
        } else if (uscalar_node(n)) {
            if (type == GRETL_TYPE_DOUBLE) {
                n->v.xval = *(double *) data;
            } else if (type == GRETL_TYPE_MATRIX) {
                /* allow type-mutation */
                n->t = MAT;
                n->v.m = (gretl_matrix *) data;
            } else {
                reattach_data_error(n, p);
            }
        } else if (n->t == MAT) {
            if (type == GRETL_TYPE_MATRIX) {
                n->v.m = data;
            } else if (type == GRETL_TYPE_DOUBLE) {
                /* allow type-mutation */
                n->t = NUM;
                n->v.xval = *(double *) data;
            } else {
                reattach_data_error(n, p);
            }
        } else if (n->t == LIST && type == GRETL_TYPE_LIST) {
            n->v.ivec = data;
        } else if (n->t == BUNDLE && type == GRETL_TYPE_BUNDLE) {
            n->v.b = data;
        } else if (n->t == STR && type == GRETL_TYPE_STRING) {
            n->v.str = data;
        } else if (n->t == ARRAY && type == GRETL_TYPE_ARRAY) {
            n->v.a = data;
        } else if (n->t == OSL) {
            n->v.ptr = data;
        } else {
            reattach_data_error(n, p);
        }
    }
}

static void node_type_error (int ntype, int argnum, int goodt,
                             NODE *bad, parser *p)
{
    const char *nstr;

    if (ntype == 0) {
        p->err = E_TYPES;
        return;
    } else if (bad != NULL && bad->t == UNDEF) {
        p->err = E_UNKVAR;
        return;
    }

    parser_ensure_error_buffer(p);

    if (ntype == LAG) {
        nstr = (goodt == NUM)? "lag order" : "lag variable";
    } else {
        nstr = getsymb(ntype);
    }

    if (goodt == 0 && bad != NULL) {
        if (null_node(bad)) {
            pprintf(p->prn, _("%s: insufficient arguments"), nstr);
        } else if (argnum <= 0) {
            pprintf(p->prn, _("%s: invalid argument type %s"),
                    nstr, typestr(bad->t));
        } else {
            pprintf(p->prn, _("%s: argument %d: invalid type %s"),
                    nstr, argnum, typestr(bad->t));
        }
        pputc(p->prn, '\n');
        p->err = E_TYPES;
        return;
    }

    if (ntype < OP_MAX) {
        if (argnum <= 0) {
            pprintf(p->prn, _("%s: operand should be %s"),
                    nstr, typestr(goodt));
        } else {
            pprintf(p->prn, _("%s: operand %d should be %s"),
                    nstr, argnum, typestr(goodt));
        }
    } else {
        if (argnum <= 0) {
            pprintf(p->prn, _("%s: argument should be %s"),
                    nstr, typestr(goodt));
        } else {
            pprintf(p->prn, _("%s: argument %d should be %s"),
                    nstr, argnum, typestr(goodt));
        }
    }

    if (bad != NULL) {
        pprintf(p->prn, _(", is %s"), typestr(bad->t));
        if (bad->t == STR) {
            pprintf(p->prn, " (%s)", bad->v.str);
        }
    }
    pputc(p->prn, '\n');

    if (!strcmp(nstr, "&")) {
        pputs(p->prn, "(for logical AND, please use \"&&\")\n");
    } else if (!strcmp(nstr, "|")) {
        pputs(p->prn, "(for logical OR, please use \"||\")\n");
    }

    p->err = E_TYPES;
}

static int node_is_true (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    return !na(x) && x != 0.0;
}

static int node_is_false (NODE *n, parser *p)
{
    return (node_get_scalar(n, p) == 0.0);
}

/* clean up a node holding evaluated sub-nodes of FARGS */

static void destroy_multi (NODE *m)
{
    free(m->v.bn.n);
    free(m);
}

static NODE *evaluate_multi_args (NODE *nn, parser *p)
{
    NODE *ret = bncopy(nn, &p->err);
    int i;

    for (i=0; i<nn->v.bn.n_nodes && !p->err; i++) {
        ret->v.bn.n[i] = eval(nn->v.bn.n[i], p);
    }

    return ret;
}

#define ok_nullarg_sym(s) (s == FARGS || s > FP_MAX || s == U_ADDR)

/* determine if an EMPTY node is acceptable in context */

static int empty_ok (NODE *n, parser *p)
{
    if (n->vname == NULL) {
        /* regular EMPTY node, OK */
        return 1;
    } else if (n->parent != NULL && ok_nullarg_sym(n->parent->t)) {
        /* null argument pass-through, OK */
        return 1;
    } else {
        p->err = E_TYPES;
        gretl_errmsg_sprintf("'%s': no such object", n->vname);
        return 0;
    }
}

/* core function: evaluate the parsed syntax tree */

static NODE *eval (NODE *t, parser *p)
{
    NODE *l = NULL, *m = NULL, *r = NULL;
    NODE *multi = NULL;
    NODE *ret = NULL;
    int ldone = 0;
    int rdone = 0;

    if (t == NULL) {
        /* catch NULL node right away */
        return NULL;
    }

#if EDEBUG
    if (t->vname != NULL) {
        fprintf(stderr, "eval: t = %p ('%s', vname=%s)",
                (void *) t, getsymb(t->t), t->vname);
    } else {
        fprintf(stderr, "eval: t = %p ('%s')", (void *) t, getsymb(t->t));
    }
    if (t->parent != NULL) {
        fprintf(stderr, ", parent %s", getsymb(t->parent->t));
    }
    if (p->flags & P_COMPILE) {
        fputs(", compile", stderr);
    }
    if (p->flags & P_EXEC) {
        fputs(", exec", stderr);
    }
    fputc('\n', stderr);
#endif

    /* handle terminals first */
    if (t->t == MSPEC) {
        return t;
    } else if (t->t == EMPTY) {
        if (empty_ok(t, p)) {
            return t;
        } else {
            return NULL;
        }
    } else if (t->t >= NUM && t->t <= STR) {
        if (exestart(p) && uvar_node(t)) {
            node_reattach_data(t, p);
        }
        return t;
    }

    if (t->t == QUERY) {
        /* needs special treatment, see eval_query() */
        l = eval(t->L, p);
        if (p->err) {
            goto bailout;
        } else {
            goto do_switch;
        }
    }

    /* handle multi-argument L or R subnodes */
    if (t->L != NULL && bnsym(t->L->t)) {
        t->L->parent = t;
        multi = evaluate_multi_args(t->L, p);
        ldone = 1;
    } else if (t->R != NULL && bnsym(t->R->t)) {
        t->R->parent = t;
        multi = evaluate_multi_args(t->R, p);
        rdone = 1;
    }

    if (!p->err && t->L != NULL && !ldone) {
        t->L->parent = t;
        if (type_func_in_parentage(t)) {
            p->flags |= P_OBJQRY;
            l = eval(t->L, p);
            p->flags ^= P_OBJQRY;
            if (p->err) {
                l = revive_undef_node(l, p);
                t->L->aux = l;
            }
        } else {
            l = eval(t->L, p);
        }
        if (l == NULL && !p->err) {
            p->err = 1;
        }
    }

    if (!p->err && t->M != NULL) {
        t->M->parent = t;
        if (m_return(t->t)) {
            m = t->M;
        } else {
            m = eval(t->M, p);
            if (m == NULL && !p->err) {
                p->err = 1;
            }
        }
    }

    if (!p->err && t->R != NULL && !rdone) {
        t->R->parent = t;
        if (r_return(t->t)) {
            r = t->R;
        } else {
            if (t->t == B_AND || t->t == B_OR) {
                /* logical operators: avoid redundant evaluation */
                if (l != NULL && l->t == NUM) {
                    if ((t->t == B_AND && node_is_false(l, p)) ||
                        (t->t == B_OR && node_is_true(l, p))) {
                        /* no need to evaluate @r */
                        ret = l;
                        if (is_aux_node(ret)) {
                            /* mark return node as proxy */
                            ret->flags |= PRX_NODE;
                        }
                        goto finish;
                    }
                }
            } else if (t->t == F_GENSERIES && t->R->t == UFUN) {
                p->flags |= P_UFRET;
            }
            if (r == NULL && !p->err) {
                r = eval(t->R, p);
                if (r == NULL && !p->err) {
                    p->err = 1;
                }
            }
        }
    }

    if (p->err) {
        goto bailout;
    }

 do_switch:

    /* establish convenience pointer */
    p->aux = t->aux;

    switch (t->t) {
    case DBUNDLE:
        if (t->vname != NULL) {
            /* pseudo-bundle: name of model */
            ret = t;
        } else {
            /* built-in bundle identifed by idnum */
            ret = dollar_bundle_node(t, p);
        }
        break;
    case UNDEF:
        if (type_func_in_parentage(t)) {
            ret = t;
        } else {
            ret = maybe_rescue_undef_node(t, p);
        }
        break;
    case U_ADDR:
        if (!uvar_node(t->L) && t->L->t != OSL) {
            p->err = E_DATA;
        } else if (exestart(p)) {
            node_reattach_data(t->L, p);
        }
        ret = t;
        break;
    case WLIST:
        ret = wildlist_node(t, p);
        break;
    case DUM:
        if (t->v.idnum == DUM_DATASET) {
            ret = dataset_list_node(p);
        } else if (t->v.idnum == DUM_TREND) {
            ret = trend_node(p);
        } else if (t->v.idnum == DUM_END) {
            ret = array_last_node(t, p);
        } else {
            /* otherwise treat as terminal */
            ret = t;
        }
        break;
    case NUM_P:
    case NUM_M:
        if (exestart(p)) {
            node_reattach_data(t, p);
        }
        ret = scalar_postfix_node(t, p);
        break;
    case B_ADD:
    case B_SUB:
    case B_MUL:
    case B_DIV:
    case B_MOD:
    case B_POW:
    case B_AND:
    case B_OR:
    case B_EQ:
    case B_NEQ:
    case B_GT:
    case B_LT:
    case B_GTE:
    case B_LTE:
        /* arithmetic and logical binary operators: be as
           flexible as possible with regard to argument types
        */
        if (t->t == B_ADD && l->t == STR && r->t == NUM) {
            ret = string_offset(l, r, p);
        } else if ((t->t == B_EQ || t->t == B_NEQ) &&
                   l->t == STR && r->t == STR) {
            ret = compare_strings(l, r, t->t, p);
        } else if (l->t == NUM && r->t == NUM) {
            ret = scalar_calc(l, r, t->t, p);
        } else if (l->t == BUNDLE && r->t == BUNDLE) {
            ret = bundle_op(l, r, t->t, p);
        } else if (l->t == ARRAY && r->t == ARRAY) {
            ret = array_op(l, r, t->t, p);
        } else if (strvals_node(l, p) && strvals_node(r, p)) {
            ret = strvals_calc(l, r, t, p);
        } else if (series_calc_nodes(l, r)) {
            ret = series_calc(l, r, t->t, p);
        } else if (l->t == MAT && r->t == MAT) {
            if (bool_comp(t->t)) {
                ret = matrix_bool(l, r, t->t, p);
            } else {
                ret = matrix_matrix_calc(l, r, t->t, p);
            }
        } else if ((l->t == MAT && r->t == NUM) ||
                   (l->t == NUM && r->t == MAT)) {
            ret = matrix_scalar_calc(l, r, t->t, p);
        } else if ((l->t == MAT && r->t == SERIES) ||
                   (l->t == SERIES && r->t == MAT)) {
            ret = matrix_series_calc(l, r, t->t, p);
        } else if (t->t >= B_EQ && t->t <= B_NEQ &&
                   ((l->t == SERIES && r->t == STR) ||
                    (l->t == STR && r->t == SERIES))) {
            ret = series_string_calc(l, r, t->t, p);
        } else if ((t->t == B_AND || t->t == B_OR || t->t == B_ADD || t->t == B_SUB) &&
                   ok_list_node_plus(l) && ok_list_node_plus(r)) {
            ret = list_list_op(l, r, t->t, p);
        } else if (t->t == B_POW && ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = list_list_op(l, r, t->t, p);
        } else if (bool_comp(t->t)) {
            if (ok_list_node(l, p) && (r->t == NUM || r->t == SERIES)) {
                ret = list_bool_comp(l, r, t->t, 0, p);
            } else if (ok_list_node(r, p) && (l->t == NUM || l->t == SERIES)) {
                ret = list_bool_comp(r, l, t->t, 1, p);
            } else if (ok_list_node(l, p) && ok_list_node(r, p)) {
                ret = list_list_comp(r, l, t->t, p);
            } else {
                p->err = E_TYPES;
            }
        } else if ((t->t == B_ADD || t->t == B_SUB) &&
                   l->t == SERIES && ok_list_node(r, p)) {
            ret = series_list_calc(l, r, t->t, p);
        } else if (t->t == B_ADD && l->t == ARRAY) {
            ret = augment_array_node(l, r, p);
        } else if (t->t == B_SUB && l->t == ARRAY) {
            ret = subtract_from_array_node(l, r, p);
        } else if (t->t == B_MUL && (l->t == STR || r->t == STR)) {
            ret = string_mult(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case B_TRMUL:
        /* matrix on left, otherwise be flexible */
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            ret = matrix_matrix_calc(l, r, t->t, p);
        } else if (l->t == MAT && r->t == SERIES) {
            ret = matrix_series_calc(l, r, t->t, p);
        } else if (l->t == MAT && r->t == EMPTY) {
            ret = matrix_transpose_node(l, p);
        } else if (l->t == NUM && r->t == EMPTY) {
            ret = l;
        } else {
            p->err = E_TYPES;
        }
        break;
    case B_DOTMULT:
    case B_DOTDIV:
    case B_DOTPOW:
    case B_DOTADD:
    case B_DOTSUB:
    case B_DOTEQ:
    case B_DOTGT:
    case B_DOTLT:
    case B_DOTGTE:
    case B_DOTLTE:
    case B_DOTNEQ:
        /* matrix-matrix or matrix-scalar binary operators:
           in addition we permit scalar-scalar to allow for
           the possibility that results that could be taken
           to be 1 x 1 matrix results have been registered
           internally as scalars.
        */
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            ret = matrix_matrix_calc(l, r, t->t, p);
        } else if ((ok_matrix_node(l) && r->t == SERIES) ||
                   (l->t == SERIES && ok_matrix_node(r))) {
            ret = matrix_series_calc(l, r, t->t, p);
        } else if ((t->t == B_DOTEQ || t->t == B_DOTNEQ) &&
                   l->t == ARRAY && r->t == STR) {
            ret = array_str_calc(l, r, t->t, p);
        } else {
            node_type_error(t->t, (l->t == MAT)? 2 : 1,
                            MAT, (l->t == MAT)? r : l, p);
        }
        break;
    case SUBSL:
    case B_RANGE:
        /* matrix sub-slice, x:y, or lag range, 'p to q' */
        ret = process_subslice(l, r, p);
        break;
    case B_VCAT:
    case F_QFORM:
    case F_HDPROD:
    case F_CMULT:
    case F_CDIV:
    case F_LSOLVE:
    case F_MRSEL:
    case F_MCSEL:
    case F_DSUM:
    case B_LDIV:
    case B_KRON:
    case F_CONV2D:
        if (ok_matrix_node(l) && ok_matrix_node(r)) {
            /* matrix-only binary operators (but promote scalars) */
            ret = matrix_matrix_calc(l, r, t->t, p);
        } else if (ok_matrix_node(l) && null_node(r) && t->t == F_HDPROD) {
            ret = matrix_to_matrix2_func(l, r, t->t, p);
        } else {
            node_type_error(t->t, (l->t == MAT)? 2 : 1,
                            MAT, (l->t == MAT)? r : l, p);
        }
        break;
    case B_HCAT:
        ret = horizontal_concat_node(l, r, p);
        break;
    case B_ELLIP:
        /* list-making ellipsis */
        if ((useries_node(l) || could_be_series_id(l, p)) &&
            (useries_node(r) || could_be_series_id(r, p))) {
            ret = ellipsis_list_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_MAT2LIST:
        /* parlay matrix into list of series */
        if (l->t != MAT) {
            node_type_error(t->t, 1, MAT, l, p);
        } else if (!null_node(r) && r->t != STR) {
            node_type_error(t->t, 2, STR, r, p);
        } else {
            ret = mat2list_node(l, r, p);
        }
        break;
    case B_JOIN:
        /* list join with separator */
        if (ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = list_join_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_MSORTBY:
    case F_CSWITCH:
        /* matrix on left, scalar on right */
        if (l->t == MAT && null_or_scalar(r)) {
            ret = matrix_scalar_func(l, r, t->t, p);
        } else if (l->t == MAT) {
            node_type_error(t->t, 2, NUM, r, p);
        } else {
            node_type_error(t->t, 1, MAT, l, p);
        }
        break;
    case F_DISTANCE:
        /* matrix on left, plus optional terms */
        if (l->t == MAT) {
            ret = matrix_string_func(l, m, r, t->t, p);
        } else {
            node_type_error(t->t, 1, MAT, l, p);
        }
        break;
    case F_MSPLITBY:
        /* matrix on left, vector, optional boolean */
        if (ok_matrix_node(l) && ok_matrix_node(m)) {
            ret = matrix_vector_func(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_LLAG:
        if (null_node(m)) {
            p->err = E_ARGS;
        } else if (ok_matrix_node(l) && m->t != MAT && ok_list_node(m, p)) {
            ret = list_make_lags(l, m, r, p);
        } else if (ok_matrix_node(l) && m->t == MAT) {
            ret = matrix_make_lags(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_STDIZE:
        if (!null_or_scalar(m)) {
            p->err = e_types(m);
        } else if (l->t == SERIES) {
            ret = series_stdize(l, m, p);
        } else if (l->t == LIST) {
            ret = list_stdize(l, m, p);
        } else if (l->t == MAT) {
            ret = eval_3args_func(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_HFLAG:
        if (scalar_node(l) && scalar_node(m) && ok_list_node(r, p)) {
            ret = hf_list_make_lags(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_HFLIST:
        if (l->t == MAT && scalar_node(m) && r->t == STR) {
            ret = hf_list_node(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_HFDIFF:
    case F_HFLDIFF:
        if (ok_list_node(l, p) && null_or_scalar(r)) {
            ret = apply_list_func(l, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case U_NEG:
    case U_POS:
    case U_NOT:
    case F_ABS:
    case F_SGN:
    case F_TOINT:
    case F_CEIL:
    case F_FLOOR:
    case F_ROUND:
    case F_SIN:
    case F_COS:
    case F_TAN:
    case F_ASIN:
    case F_ACOS:
    case F_ATAN:
    case F_SINH:
    case F_COSH:
    case F_TANH:
    case F_ASINH:
    case F_ACOSH:
    case F_ATANH:
    case F_LOG:
    case F_LOG10:
    case F_LOG2:
    case F_EXP:
    case F_SQRT:
    case F_CNORM:
    case F_DNORM:
    case F_QNORM:
    case F_LOGISTIC:
    case F_GAMMA:
    case F_LNGAMMA:
    case F_DIGAMMA:
    case F_TRIGAMMA:
    case F_INVMILLS:
    case F_EASTER:
        /* functions taking one numeric argument, any type */
        if (l->t == NUM) {
            ret = apply_scalar_func(l, t, p);
        } else if (l->t == SERIES) {
            if (cast_series_to_list(p, l, t->t)) {
                ret = apply_list_func(l, NULL, t->t, p);
            } else {
                ret = apply_series_func(l, t, p);
            }
        } else if (l->t == MAT) {
            ret = apply_matrix_func(l, t, p);
        } else if (ok_list_node(l, p) && t->t == F_LOG) {
            ret = apply_list_func(l, NULL, t->t, p);
        } else {
            p->err = e_types(l);
        }
        break;
    case F_REAL:
    case F_IMAG:
        if (l->t == MAT || l->t == NUM) {
            ret = apply_matrix_func(l, t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_CARG:
    case F_CONJ:
    case F_CMOD:
    case F_CQUAD:
        if (complex_node(l)) {
            ret = apply_matrix_func(l, t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_ATAN2:
    case F_BINCOEFF:
        if ((l->t == NUM || l->t == MAT || l->t == SERIES) &&
            (r->t == NUM || r->t == MAT || r->t == SERIES)) {
            ret = flexible_2arg_node(l, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_DUMIFY:
    case F_CDUMIFY:
        /* series or list argument wanted */
        if (l->t == SERIES || l->t == LIST) {
            if (t->t == F_DUMIFY) {
                ret = dummify_func(l, r, p);
            } else {
                ret = cdummify_func(l, p);
            }
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_GETINFO:
        /* named series (or ID) argument wanted */
        if (useries_node(l) || l->t == NUM) {
            ret = get_info_on_series(l, p);
        } else {
            p->err = e_types(l);
        }
        break;
    case F_SEASONALS:
        /* two optional args: int, bool */
        if (!null_or_scalar(l)) {
            node_type_error(t->t, 1, NUM, l, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t->t, 2, NUM, r, p);
        } else {
            ret = seasonals_node(l, r, p);
        }
        break;
    case F_MISSZERO:
    case F_ZEROMISS:
        /* one series, scalar or matrix argument needed */
        if (l->t == SERIES) {
            ret = apply_series_func(l, t, p);
        } else if (l->t == NUM) {
            ret = apply_scalar_func(l, t, p);
        } else if (l->t == MAT) {
            ret = apply_matrix_func(l, t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_MISSING:
    case F_DATAOK:
        /* series, scalar or list argument needed */
        if (l->t == MAT) {
            if (t->t == F_DATAOK) {
                ret = matrix_to_matrix_func(l, NULL, t->t, p);
            } else {
                ret = matrix_isnan_node(l, p);
            }
        } else if (l->t == SERIES) {
            ret = apply_series_func(l, t, p);
        } else if (l->t == NUM) {
            ret = apply_scalar_func(l, t, p);
        } else if (l->t == LIST) {
            ret = list_ok_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_ISNAN:
        /* scalar or matrix */
        if (scalar_node(l)) {
            ret = scalar_isnan_node(l, p);
        } else if (l->t == MAT) {
            ret = matrix_isnan_node(l, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_SLEEP:
        if (scalar_node(l)) {
            ret = sleep_node(l, t->t, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_BARRIER:
        if (l->t == EMPTY) {
            ret = mpi_barrier_node(p);
        } else {
            node_type_error(t->t, 0, EMPTY, l, p);
        }
        break;
    case LAG:
        if (p->targ == LIST) {
            ret = get_lag_list(l, r, p);
        } else if (l->t == SERIES && scalar_node(r)) {
            ret = series_lag(l, r, p);
        } else if (l->t != SERIES) {
            node_type_error(t->t, 1, SERIES, l, p);
        } else {
            node_type_error(t->t, 2, NUM, r, p);
        }
        break;
    case F_LJUNGBOX:
    case F_POLYFIT:
        /* series on left, scalar on right */
        if (l->t != SERIES) {
            node_type_error(t->t, 1, SERIES, l, p);
        } else if (!scalar_node(r)) {
            node_type_error(t->t, 2, NUM, r, p);
        } else if (t->t == F_LJUNGBOX) {
            ret = series_ljung_box(l, r, p);
        } else if (t->t == F_POLYFIT) {
            ret = series_polyfit(l, r, p);
        }
        break;
    case OBS:
        if (l->t != SERIES) {
            node_type_error(t->t, 1, SERIES, l, p);
        } else if (!scalar_node(r) && r->t != STR) {
            node_type_error(t->t, 2, NUM, r, p);
        } else if (t->flags & LHT_NODE) {
            ret = lhs_terminal_node(t, l, r, p);
        } else {
            ret = series_obs(l, r, p);
        }
        break;
    case OSL:
        /* object plus subspec */
        if (t->flags & LHT_NODE) {
            ret = lhs_terminal_node(t, l, r, p);
        } else if (l->t == U_ADDR) {
            ret = process_OSL_address(t, l, r, p);
        } else {
            ret = subobject_node(l, r, p);
        }
        break;
    case SLRAW:
        /* unevaluated object slice spec */
        ret = mspec_node(l, r, p);
        break;
    case LISTVAR:
        ret = list_member_by_name(l, r, p);
        break;
    case BMEMB:
    case F_INBUNDLE:
        /* name of bundle plus string */
        if (l->t == BUNDLE && r->t == STR) {
            if (t->t == BMEMB) {
                if (t->flags & LHT_NODE) {
                    ret = lhs_terminal_node(t, l, r, p);
                } else {
                    ret = get_bundle_member(l, r, p);
                }
            } else {
                ret = test_bundle_key(l, r, p);
            }
        } else if (l->t == BUNDLE) {
            node_type_error(t->t, 2, STR, r, p);
        } else {
            node_type_error(t->t, 1, BUNDLE, l, p);
        }
        break;
    case F_GETKEYS:
    case HF_JBTERMS:
        if (l->t == BUNDLE) {
            ret = get_bundle_array(l, t->t, p);
        } else {
            node_type_error(t->t, 0, BUNDLE, l, p);
        }
        break;
    case DBMEMB:
        /* name of $-bundle plus string */
        if (l->t == BUNDLE && r->t == STR) {
            ret = get_bundle_member(l, r, p);
        } else if (l->t == DBUNDLE && r->t == STR) {
            ret = model_var_node(l, r, p);
        } else if (r->t != STR) {
            node_type_error(t->t, 2, STR, r, p);
        } else {
            node_type_error(t->t, 1, BUNDLE, l, p);
        }
        break;
    case HF_FDEPTH:
        ret = function_depth_node(l, p);
        break;
    case F_CURL:
        ret = curl_bundle_node(l, p);
        break;
    case F_LPSOLVE:
        if (l->t == BUNDLE) {
            ret = lpsolve_bundle_node(l, p);
        } else {
            node_type_error(t->t, 0, BUNDLE, l, p);
        }
        break;
    case F_SVM:
        ret = svm_driver_node(t, multi, p);
        break;
    case F_TYPESTR:
        /* numerical type code to string */
        if (scalar_node(l)) {
            ret = type_string_node(l, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_LDIFF:
    case F_SDIFF:
    case F_ODEV:
    case F_INTERPOL:
        if (l->t == SERIES && cast_series_to_list(p, l, t->t)) {
            ret = apply_list_func(l, NULL, t->t, p);
        } else if (l->t == SERIES ||
                   (t->t != F_ODEV && t->t != F_INTERPOL && l->t == MAT)) {
            ret = series_series_func(l, r, NULL, t->t, p);
        } else if (ok_list_node(l, p)) {
            ret = apply_list_func(l, NULL, t->t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_DROPCOLL:
        /* list argument is required on left, optional scalar
           on the right */
        if (l->t == LIST) {
            ret = apply_list_func(l, r, t->t, p);
        } else {
            node_type_error(t->t, 0, LIST, l, p);
        }
        break;
    case F_HPFILT:
    case F_FRACDIFF:
    case F_FRACLAG:
    case F_BOXCOX:
    case F_PNOBS:
    case F_PMIN:
    case F_PMAX:
    case F_PSUM:
    case F_PMEAN:
    case F_PXSUM:
    case F_PXNOBS:
    case F_PSD:
    case F_DESEAS:
    case F_TRAMOLIN:
        /* series argument needed */
        if (l->t == SERIES || l->t == MAT) {
            if (t->t == F_HPFILT) {
                ret = series_series_func(l, m, r, t->t, p);
            } else {
                ret = series_series_func(l, r, NULL, t->t, p);
            }
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_PSHRINK:
        if (l->t == SERIES) {
            int noskip = node_get_bool(r, p, 0);

            if (!p->err) {
                ret = do_panel_shrink(l, noskip, p);
            }
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_PEXPAND:
        if (l->t == MAT) {
            ret = do_panel_expand(l, r, p);
        } else {
            node_type_error(t->t, 0, MAT, l, p);
        }
        break;
    case F_CUM:
    case F_DIFF:
    case F_RESAMPLE:
    case F_RANKING:
        /* series or matrix argument */
        if (l->t == SERIES) {
            if (cast_series_to_list(p, l, t->t)) {
                ret = apply_list_func(l, NULL, t->t, p);
            } else {
                ret = series_series_func(l, r, NULL, t->t, p);
            }
        } else if (l->t == MAT) {
            if (t->t == F_RESAMPLE) {
                ret = eval_3args_func(l, m, r, t->t, p);
            } else {
                ret = matrix_to_matrix_func(l, r, t->t, p);
            }
        } else if (t->t == F_DIFF && ok_list_node(l, p)) {
            ret = apply_list_func(l, NULL, t->t, p);
        } else if (t->t == F_RESAMPLE && ok_list_node(l, p) &&
                   null_node(r)) {
            ret = apply_list_func(l, NULL, t->t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_SORT:
    case F_DSORT:
        /* series or vector or string array argument needed */
        if (l->t == SERIES || l->t == MAT || l->t == NUM) {
            ret = vector_sort(l, t->t, p);
        } else if (l->t == ARRAY) {
            ret = array_sort_node(l, t->t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_FLATTEN:
    case F_INSTRINGS:
    case F_ASORT:
        if (l->t == ARRAY) {
            if (t->t == F_INSTRINGS) {
                ret = array_func_node(l, m, r, t->t, p);;
            } else {
                ret = array_func_node(l, r, NULL, t->t, p);
            }
        } else {
            node_type_error(t->t, 0, ARRAY, l, p);
        }
        break;
    case F_VALUES:
    case F_UNIQ:
    case F_PERGM:
    case F_IRR:
        /* series or vector argument needed */
        if (l->t == SERIES || l->t == MAT || l->t == NUM) {
            if (t->t == F_PERGM) {
                ret = pergm_node(l, r, p);
            } else if (t->t == F_VALUES || t->t == F_UNIQ) {
                ret = vector_values(l, t->t, p);
            } else if (t->t == F_IRR) {
                ret = do_irr(l, p);
            } else {
                ret = vector_sort(l, t->t, p);
            }
        } else {
            p->err = e_types(l);
        }
        break;
    case F_SUM:
    case F_SUMALL:
    case F_MEAN:
    case F_SD:
    case F_VCE:
    case F_SST:
    case F_SKEWNESS:
    case F_KURTOSIS:
    case F_MEDIAN:
    case F_GINI:
    case F_NOBS:
    case F_T1:
    case F_T2:
        /* functions taking series arg (mostly), returning scalar */
        if (l->t == SERIES || l->t == MAT) {
            ret = series_scalar_func(l, t->t, r, p);
        } else if ((t->t == F_MEAN || t->t == F_SD ||
                    t->t == F_VCE || t->t == F_SUM ||
                    t->t == F_MEDIAN || t->t == F_NOBS)
                   && ok_list_node(l, p)) {
            /* list -> series also acceptable for these cases */
            ret = list_to_series_func(l, t->t, p);
        } else if (l->t == NUM) {
            ret = pretend_matrix_scalar_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, SERIES, l, p);
        }
        break;
    case F_MIN:
    case F_MAX:
        if (null_node(r)) {
            /* "within arg" min or max */
            ret = within_minmax_node(l, t->t, p);
        } else {
            /* "between args" min or max */
            ret = between_minmax_node(l, r, t->t, p);
        }
        break;
    case F_ECDF:
    case F_NORMTEST:
        /* series or vector (plus optional string arg for
           normtest); returns matrix */
        if (l->t != SERIES && l->t != MAT) {
            p->err = e_types(l);
        } else if (null_or_string(r)) {
            ret = series_matrix_node(l, r, t->t, p);
        } else {
            p->err = e_types(r);
        }
        break;
    case F_LRVAR:
    case F_ISCONST:
    case F_ISDUMMY:
        /* takes series and scalar arg, returns scalar */
        if (l->t == SERIES || l->t == MAT) {
            if (t->t == F_ISCONST || t->t == F_ISDUMMY ) {
                ret = isconst_or_dum_node(l, r, p, t->t);
            } else if (t->t == F_LRVAR) {
                ret = series_scalar_scalar_func(l, m, r, t->t, p);
            }
        } else {
            node_type_error(t->t, 1, SERIES, l, p);
        }
        break;
    case F_NPV:
        if (l->t != SERIES && l->t != MAT && l->t != NUM) {
            node_type_error(t->t, 1, SERIES, l, p);
        } else if (!scalar_node(r)) {
            node_type_error(t->t, 2, NUM, r, p);
        } else {
            ret = series_scalar_scalar_func(l, r, NULL, t->t, p);
        }
        break;
    case F_QUANTILE:
        if (l->t == SERIES) {
            if (scalar_node(r)) {
                ret = series_scalar_scalar_func(l, r, NULL, t->t, p);
            } else {
                node_type_error(t->t, 2, NUM, r, p);
            }
        } else if (l->t == MAT) {
            if (r->t == MAT || scalar_node(r)) {
                ret = matrix_quantiles_node(l, r, p);
            } else {
                node_type_error(t->t, 2, MAT, r, p);
            }
        } else {
            node_type_error(t->t, 1, (r->t == MAT)? MAT : SERIES,
                            l, p);
        }
        break;
    case F_RUNIFORM:
    case F_RNORMAL:
        /* functions taking zero or two scalars as args */
        if (scalar_node(l) && scalar_node(r)) {
            ret = series_fill_func(l, r, t->t, p);
        } else if (l->t == EMPTY && r->t == EMPTY) {
            ret = series_fill_func(l, r, t->t, p);
        } else {
            node_type_error(t->t, (l->t == NUM)? 2 : 1,
                            NUM, (l->t == NUM)? r : l, p);
        }
        break;
    case F_COV:
    case F_COR:
    case F_NAALEN:
    case F_KMEIER:
    case F_CORRESP:
        /* functions taking two series/vectors as args, mostly */
        if ((l->t == SERIES || l->t == MAT || l->t == NUM) &&
            (r->t == SERIES || r->t == MAT || r->t == NUM)) {
            ret = series_2_func(l, r, t->t, p);
        } else if ((l->t == SERIES || l->t == MAT) &&
                   null_node(r) &&
                   (t->t == F_NAALEN || t->t == F_KMEIER)) {
            ret = series_2_func(l, NULL, t->t, p);
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_FCSTATS:
        /* two series or vectors, plus optional boolean */
        if ((l->t == SERIES || l->t == MAT) &&
            (m->t == SERIES || m->t == MAT || m->t == LIST)) {
            ret = fcstats_node(l, m, r, p);
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_NPCORR:
        /* two series or vectors, plus optional control string */
        if ((l->t == SERIES || l->t == MAT) &&
            (m->t == SERIES || m->t == MAT) &&
            null_or_string(r)) {
            ret = npcorr_node(l, m, r, p);
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_MXTAB:
        /* functions taking two series or matrices as args and returning
           a matrix */
        if ((l->t == SERIES && r->t == SERIES) || (l->t == MAT && r->t == MAT)) {
            ret = mxtab_func(l, r, p);
        } else {
            node_type_error(t->t, (l->t == SERIES)? 2 : 1,
                            SERIES, (l->t == SERIES)? r : l, p);
        }
        break;
    case F_SORTBY:
        /* takes two series as args, returns series */
        if (l->t == SERIES && r->t == SERIES) {
            ret = series_sort_by(l, r, p);
        } else {
            node_type_error(t->t, (l->t == SERIES)? 2 : 1,
                            SERIES, (l->t == SERIES)? r : l, p);
        }
        break;
    case F_IMAT:
    case F_ZEROS:
    case F_ONES:
    case F_MUNIF:
    case F_MNORM:
    case F_RANDPERM:
        /* matrix-creation functions */
        if (scalar_node(l) && null_or_scalar(r)) {
            ret = matrix_fill_func(l, r, t->t, p);
        } else if (!scalar_node(l)) {
            node_type_error(t->t, 1, NUM, l, p);
        } else {
            node_type_error(t->t, 2, NUM, r, p);
        }
        break;
    case HF_VCNORM:
        if (l->t == MAT && scalar_node(r)) {
            ret = mvnormal_node(l, r, p);
        } else if (l->t != MAT) {
            node_type_error(t->t, 1, MAT, l, p);
        } else {
            node_type_error(t->t, 2, NUM, r, p);
        }
        break;
    case F_SUMC:
    case F_SUMR:
    case F_PRODC:
    case F_PRODR:
    case F_MEANC:
    case F_MEANR:
    case F_MCOV:
    case F_MCORR:
    case F_CHOL:
    case F_PSDROOT:
    case F_INV:
    case F_GINV:
    case F_DIAG:
    case F_TRANSP:
    case F_MREV:
    case F_VEC:
    case F_VECH:
    case F_UNVECH:
    case F_UPPER:
    case F_LOWER:
    case F_NULLSPC:
    case F_MEXP:
    case F_MLOG:
    case F_MINC:
    case F_MAXC:
    case F_MINR:
    case F_MAXR:
    case F_IMINC:
    case F_IMAXC:
    case F_IMINR:
    case F_IMAXR:
    case F_FFT:
    case F_FFTI:
    case F_POLROOTS:
    case F_CTRANS:
        /* matrix -> matrix functions */
        if (l->t == MAT || l->t == NUM) {
            ret = matrix_to_matrix_func(l, r, t->t, p);
        } else if (t->t == F_MREV && l->t == LIST) {
            ret = list_reverse_node(l, p);
        } else {
            node_type_error(t->t, 1, MAT, l, p);
        }
        break;
    case HF_GLASSO:
        if (l->t == MAT && r->t == BUNDLE) {
            ret = glasso_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_INVPD:
        /* matrix -> matrix with optional scalar pointer */
        if (l->t == MAT || l->t == NUM) {
            ret = invpd_node(l, r, p);
        } else {
            node_type_error(t->t, 1, MAT, l, p);
        }
        break;
    case F_ROWS:
    case F_COLS:
    case F_NORM1:
    case F_INFNORM:
    case F_RCOND:
    case F_CNUMBER:
    case F_RANK:
        /* matrix -> scalar functions */
        if (l->t == MAT || l->t == NUM) {
            ret = matrix_to_scalar_func(l, r, t->t, p);
        } else {
            node_type_error(t->t, 0, MAT, l, p);
        }
        break;
    case F_DET:
    case F_LDET:
    case F_TRACE:
        if (l->t == MAT || l->t == NUM) {
            ret = matrix_to_alt_node(l, t->t, p);
        }
        break;
    case F_MREAD:
    case F_BREAD:
        if (l->t != STR) {
            node_type_error(t->t, 1, STR, l, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t->t, 2, NUM, r, p);
        } else {
            ret = read_object_func(l, r, t->t, p);
        }
        break;
    case F_EIGSYM:
        /* matrix -> matrix functions, with single indirect return */
        if (l->t != MAT && l->t != NUM) {
            node_type_error(t->t, 1, MAT, l, p);
        } else if (r->t != U_ADDR && r->t != EMPTY) {
            node_type_error(t->t, 2, U_ADDR, r, p);
        } else {
            ret = matrix_to_matrix2_func(l, r, t->t, p);
        }
        break;
    case F_COMPLEX:
        if ((l->t == MAT || l->t == NUM) &&
            (r->t == MAT || null_or_scalar(r))) {
            ret = complex_matrix_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_FDJAC:
    case F_NUMHESS:
        /* matrix, fncall, optional scalar */
        if (l->t == MAT && m->t == STR) {
            ret = numeric_jacobian_or_hessian(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_MWRITE:
        /* matrix, with string as second arg */
        if ((l->t == MAT || l->t == NUM) &&
            m->t == STR && null_or_scalar(r)) {
            ret = matrix_file_write(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_BWRITE:
        /* bundle, with string as second arg */
        if (l->t == BUNDLE && m->t == STR && null_or_scalar(r)) {
            ret = bundle_file_write(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_BFGSMAX:
        /* matrix-pointer, plus one or two string args */
        if ((l->t == U_ADDR || l->t == MAT) && m->t == STR) {
            ret = BFGS_maximize(l, m, r, p, t);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_BFGSCMAX:
        ret = BFGS_constrained_max(t, multi, p);
        break;
    case F_SIMANN:
    case F_NMMAX:
    case F_GSSMAX:
        /* matrix(-pointer), plus string and scalar args */
        if ((l->t == U_ADDR || l->t == MAT) && m->t == STR) {
            ret = deriv_free_node(l, m, r, p, t);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_FZERO:
        if (l->t == STR) {
            ret = fzero_node(l, m, r, p);
        } else {
            p->err = e_types(l);
        }
        break;
    case F_IMHOF:
        /* matrix, scalar as second arg */
        if (l->t == MAT && scalar_node(r)) {
            ret = matrix_imhof(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_BKW:
        /* matrix, string(s) as second optional arg,
           quiet flag as optional third arg */
        if (l->t == MAT) {
            ret = bkw_node(l, m, r, p);
        } else {
            p->err = e_types(l);
        }
        break;
    case F_FEVD:
        /* integer target, source plus optional bundle */
        if (scalar_node(l) && null_or_scalar(m)) {
            ret = fevd_node(l, m, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_CNAMESET:
    case F_RNAMESET:
        /* matrix, with (list, string or strings array) as second arg */
        if (l->t == MAT && (ok_list_node(r, p) || r->t == STR || r->t == ARRAY)) {
            ret = matrix_add_names(l, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_CNAMEGET:
    case F_RNAMEGET:
        /* matrix, scalar as second arg */
        if (l->t == MAT && null_or_scalar(r)) {
            ret = matrix_get_col_or_row_name(t->t, l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_RANDINT:
    case F_KPSSCRIT:
        /* two scalars */
        if (scalar_node(l) && scalar_node(r)) {
            if (t->t == F_KPSSCRIT) {
                ret = kpss_crit_node(l, r, p);
            } else {
                ret = two_scalars_func(l, r, t->t, p);
            }
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_MSHAPE:
    case F_SVD:
    case F_QR:
    case F_EIGEN:
    case F_EIGGEN:
    case F_SCHUR:
    case F_TRIMR:
    case F_CORRGM:
    case F_SEQ:
    case F_REPLACE:
    case F_STRNCMP:
    case F_WEEKDAY:
    case F_DAYSPAN:
    case F_SMPLSPAN:
    case F_MONTHLEN:
    case F_EPOCHDAY:
    case F_KDENSITY:
    case F_SETNOTE:
    case F_BWFILT:
    case F_VARSIMUL:
    case F_MLAG:
    case F_EIGSOLVE:
    case F_PRINCOMP:
    case F_HALTON:
    case F_AGGRBY:
    case F_IWISHART:
    case F_MWEIGHTS:
    case F_MGRADIENT:
    case F_LRCOVAR:
    case F_BRENAME:
    case F_ISOWEEK:
    case F_STACK:
    case F_VMA:
    case F_BCHECK:
    case F_SPHCORR:
    case F_SDC:
    case F_CDEMEAN:
    case HF_REGLS:
        /* built-in functions taking three args */
        if (t->t == F_REPLACE) {
            ret = replace_value(l, m, r, p);
        } else if (t->t == F_EPOCHDAY) {
            ret = eval_epochday(l, m, r, p);
        } else {
            ret = eval_3args_func(l, m, r, t->t, p);
        }
        break;
    case F_GEOPLOT:
        ret = geoplot_node(l, m, r, p);
        break;
    case F_PRINTF:
    case F_SPRINTF:
        if (l->t == STR && null_or_string(r)) {
            ret = print_func_node(l, r, t->t, p);
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_SSCANF:
        if ((l->t == STR || l->t == ARRAY) &&
            (m->t == STR && r->t == STR)) {
            ret = sscanf_node(l, m, r, p);
        } else {
            p->err = E_INVARG;
        }
        break;
    case F_BESSEL:
        /* functions taking one char, one scalar and one
           scalar/series/matrix as args */
        if (l->t != STR) {
            node_type_error(t->t, 1, STR, l, p);
        } else if (!scalar_node(m)) {
            node_type_error(t->t, 2, NUM, m, p);
        } else if (r->t != NUM && r->t != SERIES && r->t != MAT) {
            node_type_error(t->t, 3, NUM, r, p);
        } else {
            ret = eval_bessel_func(l, m, r, p);
        }
        break;
    case F_BKFILT:
    case F_MOLS:
    case F_MPOLS:
    case F_MRLS:
    case F_FILTER:
    case F_MCOVG:
    case F_NRMAX:
    case F_LOESS:
    case F_GHK:
    case F_QUADTAB:
    case F_QLRPVAL:
    case F_BOOTCI:
    case F_BOOTPVAL:
    case F_MOVAVG:
    case F_IRF:
    case F_NADARWAT:
    case F_FEVAL:
    case F_CHOWLIN:
    case F_HYP2F1:
    case F_TDISAGG:
    case F_MIDASMULT:
    case F_COMMUTE:
    case F_TOEPSOLV:
    case F_RGBMIX:
    case HF_FELOGITR:
        /* built-in functions taking more than three args */
        if (multi == NULL) {
            fprintf(stderr, "INTERNAL ERROR: @multi is NULL\n");
            p->err = E_DATA;
        } else if (t->t == F_FEVAL) {
            multi->parent = t;
            ret = eval_feval(t->t, multi, NULL, p);
        } else {
            ret = eval_nargs_func(t, multi, p);
        }
        break;
    case F_FEVALB:
        if (l->t == STR && r->t == BUNDLE) {
            l->parent = t;
            ret = eval_feval(t->t, l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_DEFARRAY:
    case F_DEFBUNDLE:
    case F_DEFARGS:
    case F_DEFLIST:
        ret = object_def_node(t, multi, p);
        break;
    case F_DEFMAT:
        ret = matrix_def_node(t, multi, p);
        break;
    case F_KSETUP:
    case F_KFILTER:
    case F_KSMOOTH:
    case F_KSIMUL:
    case F_KDSMOOTH:
        if (t->t == F_KSETUP || bundle_pointer_arg0(t)) {
            ret = eval_kalman_bundle_func(t, multi, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_KSIMDATA:
        if (r->t == MAT) {
            ret = kalman_data_node(l, r, p);
        } else {
            node_type_error(t->t, 2, MAT, r, p);
        }
        break;
    case F_ISOCONV:
        ret = isoconv_node(t, multi, p);
        break;
    case MVAR:
        /* variable "under" model */
        ret = model_var_node(t, NULL, p);
        break;
    case DMSTR:
        ret = dollar_str_node(t, NULL, p);
        break;
    case DVAR:
        /* dataset "dollar" variable */
        ret = dollar_var_node(t, p);
        break;
    case F_OBSNUM:
    case F_ISDISCR:
    case F_NLINES:
    case F_REMOVE:
    case F_EXISTS:
        if (l->t == STR) {
            ret = object_status(l, t, p);
        } else {
            node_type_error(t->t, 1, STR, l, p);
        }
        break;
    case F_STRLEN:
        if (l->t == STR) {
            ret = object_status(l, t, p);
        } else if (useries_node(l) || l->t == ARRAY) {
            ret = multi_str_node(l, t->t, p);
        } else {
            node_type_error(t->t, 1, STR, l, p);
        }
        break;
    case F_TYPEOF:
    case F_TYPENAME:
    case F_ISCMPLX:
        ret = generic_typeof_node(l, t, p);
        break;
    case F_NELEM:
        ret = n_elements_node(l, p);
        break;
    case F_INLIST:
    case HF_LISTINFO:
        if (l->t == LIST && t->t == HF_LISTINFO) {
            ret = list_info_node(l, r, p);
        } else if (ok_list_node(l, p)) {
            ret = in_list_node(l, r, p);
        } else {
            node_type_error(t->t, 1, LIST, l, p);
        }
        break;
    case F_PDF:
    case F_CDF:
    case F_INVCDF:
    case F_CRIT:
    case F_PVAL:
    case F_RANDGEN:
    case F_MRANDGEN:
    case F_RANDGEN1:
    case F_URCPVAL:
        if (t->t == F_URCPVAL) {
            ret = eval_urcpval(t, multi, p);
        } else {
            ret = eval_pdist(t, multi, p);
        }
        break;
    case CON:
        /* built-in constant */
        ret = retrieve_const(t, p);
        break;
    case UFUN:
        ret = eval_ufunc(t, multi, p);
        break;
    case RFUN:
        ret = eval_Rfunc(t, multi, p);
        break;
    case QUERY:
        ret = eval_query(t, l, p);
        break;
    case B_LCAT:
        /* list concatenation */
        if (ok_list_node(l, p) && ok_list_node(r, p)) {
            ret = eval_lcat(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_SQUARE:
        if (ok_list_node(l, p) && null_or_scalar(r)) {
            ret = apply_list_func(l, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_WMEAN:
    case F_WVAR:
    case F_WSD:
        /* two lists -> series, with optional boolean */
        if (ok_list_node(l, p) && ok_list_node(m, p)) {
            ret = list_list_series_func(l, m, t->t, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_LINCOMB:
        /* list + matrix -> series */
        if (ok_list_node(l, p) && ok_matrix_node(r)) {
            ret = lincomb_func(l, r, NULL, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_MLINCOMB:
        /* list + matrix + int -> series */
        if (ok_list_node(l, p) && m->t == MAT &&
            (scalar_node(r) || r->t == STR)) {
            ret = lincomb_func(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_ARGNAME:
    case F_BACKTICK:
    case F_STRSTRIP:
    case F_FIXNAME:
        if (l->t == STR) {
            ret = single_string_func(l, r, t->t, p);
        } else if (t->t == F_ARGNAME && uvar_node(l)) {
            ret = argname_from_uvar(l, r, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    case F_CCODE:
        ret = country_code_node(l, r, p);
        break;
    case F_READFILE:
        if (l->t == STR) {
            ret = readfile_node(l, r, p);
        } else {
            node_type_error(t->t, 1, STR, l, p);
        }
        break;
    case F_GETENV:
    case F_NGETENV:
        if (l->t == STR) {
            ret = do_getenv(l, t->t, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    case F_ERRORIF:
        if (r->t != STR) {
            p->err = e_types(r);
        } else {
            ret = do_errorif(l, r, p);
        }
        break;
    case F_ASSERT:
        if (l->t == NUM && r->t == STR) {
            ret = do_assert(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_CONTAINS:
        if ((r->t == MAT || r->t == NUM) &&
            (l->t == NUM || l->t == SERIES || l->t == MAT)) {
            ret = contains_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_OBSLABEL:
        if (l->t == NUM || l->t == MAT) {
            ret = int_to_string_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_VARNAME:
        if (l->t == NUM || l->t == MAT || l->t == SERIES) {
            ret = int_to_string_func(l, t->t, p);
        } else if (l->t == LIST) {
            ret = list_to_string_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_RANDSTR:
        if (scalar_node(l)) {
            ret = int_to_string_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_VARNAMES:
        if (l->t == LIST) {
            ret = list_to_string_func(l, t->t, p);
        } else {
            node_type_error(t->t, 0, LIST, l, p);
        }
        break;
    case F_VARNUM:
        if (l->t == STR) {
            ret = varnum_node(l, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    case F_TOLOWER:
    case F_TOUPPER:
        if (l->t == STR || l->t == ARRAY || useries_node(l)) {
            ret = variant_string_func(l, NULL, NULL, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_STRSUB:
    case F_REGSUB:
        if ((l->t == STR || l->t == ARRAY || strvals_node(l, p)) &&
            m->t == STR && r->t == STR) {
            ret = variant_string_func(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_SUBSTR:
        if ((l->t == STR || l->t == ARRAY || strvals_node(l, p)) &&
            scalar_node(m) && scalar_node(r)) {
            ret = variant_string_func(l, m, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_JSONGET:
    case F_XMLGET:
        if (l->t == STR && m->t == STR) {
            ret = two_string_func(l, m, r, t->t, p);
        } else if (l->t == STR && m->t == ARRAY && t->t == F_XMLGET) {
            ret = two_string_func(l, m, r, t->t, p);
        } else {
            node_type_error(t->t, (l->t == STR)? 2 : 1,
                            STR, (l->t == STR)? m : l, p);
        }
        break;
    case F_JSONGETB:
        if (l->t == STR && null_or_string(r)) {
            ret = two_string_func(l, r, NULL, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_STRSTR:
    case F_INSTRING:
        if (l->t == STR && m->t == STR) {
            ret = two_string_func(l, m, r, t->t, p);
        } else {
            node_type_error(t->t, (l->t == STR)? 2 : 1,
                            STR, (l->t == STR)? m : l, p);
        }
        break;
    case F_STRSPLIT:
        if (l->t == STR) {
            ret = strsplit_node(t->t, l, m, r, p);
        } else {
            node_type_error(t->t, 1, STR, l, p);
        }
        break;
    case F_GETLINE:
        if (l->t == STR && (null_or_string(r) || r->t == U_ADDR)) {
            ret = getline_node(l, r, p);
        } else {
            node_type_error(t->t, (l->t == STR)? 2 : 1,
                            STR, (l->t == STR)? r : l, p);
        }
        break;
    case F_ERRMSG:
        if (null_or_scalar(l)) {
            ret = errmsg_node(l, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_ISODATE:
    case F_JULDATE:
        ret = isodate_node(l, r, t->t, p);
        break;
    case F_STRFTIME:
    case F_STRFDAY:
        if (l->t == NUM || l->t == SERIES || l->t == MAT) {
            if (t->t == F_STRFTIME) {
                ret = strftime_node(l, m, r, t->t, 0, p);
            } else {
                ret = strftime_node(l, r, NULL, t->t, 0, p);
            }
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_STRPTIME:
    case F_STRPDAY:
        if (l->t == STR || l->t == SERIES || l->t == ARRAY ||
            l->t == MAT || l->t == NUM) {
            ret = strptime_node(l, r, t->t, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    case F_ATOF:
        if (l->t == STR || l->t == SERIES) {
            ret = atof_node(l, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    case F_MPI_RECV:
        ret = mpi_transfer_node(l, NULL, NULL, t->t, p);
        break;
    case F_MPI_SEND:
    case F_BCAST:
    case F_ALLREDUCE:
        if (t->t == F_ALLREDUCE && r->t != STR) {
            node_type_error(t->t, 2, STR, r, p);
        } else {
            ret = mpi_transfer_node(l, r, NULL, t->t, p);
        }
        break;
    case F_REDUCE:
    case F_SCATTER:
        if (m->t != STR) {
            node_type_error(t->t, 2, STR, m, p);
        } else if (!null_or_scalar(r)) {
            node_type_error(t->t, 3, NUM, r, p);
        } else {
            ret = mpi_transfer_node(l, m, r, t->t, p);
        }
        break;
    case F_GENSERIES:
        ret = gen_series_node(l, r, p);
        break;
    case F_ARRAY:
        ret = gen_array_node(l, p);
        break;
    case F_STRVALS:
        if (!useries_node(l)) {
            node_type_error(t->t, 0, USERIES, l, p);
        } else {
            ret = get_series_stringvals(l, r, p);
        }
        break;
    case F_STRINGIFY:
    case F_STRVSORT:
        if (useries_node(l) && (r->t == ARRAY || null_node(r))) {
            ret = stringify_series(l, r, t->t, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_DEC2BIN:
        if (l->t == NUM || l->t == MAT) {
            ret = dec2bin_node(l, p);
        } else {
            node_type_error(t->t, 0, NUM, l, p);
        }
        break;
    case F_BIN2DEC:
        if (l->t == MAT) {
            ret = bin2dec_node(l, p);
        } else {
            node_type_error(t->t, 0, MAT, l, p);
        }
        break;
    case F_BINPERMS:
        if (scalar_node(l) && scalar_node(r)) {
            ret = binperm_node(l, r, p);
        } else {
            p->err = E_TYPES;
        }
        break;
    case F_ACCESS:
        if (l->t == STR) {
            ret = access_series_node(l, p);
        } else {
            node_type_error(t->t, 0, STR, l, p);
        }
        break;
    default:
        fprintf(stderr, "eval: weird node %s (t->t = %d)\n",
                getsymb(t->t), t->t);
        p->err = E_TYPES;
        break;
    }

 finish:

    if ((p->flags & P_PRNLIST) && ret != NULL &&
        scalar_node(ret) && is_aux_node(ret)) {
        /* don't try to interpret this scalar as a series ID */
        p->err = E_DATA;
    }

    if (!p->err && ret != NULL && ret != t && is_aux_node(ret)) {
        if ((t->t == F_FEVAL || t->t == F_FEVALB) && ret->refcount > 0) {
            ; /* special case: don't attach here */
        } else {
            p->err = attach_aux_node(t, ret, p);
        }
    }

    if (multi != NULL) {
        destroy_multi(multi);
    }

    if (!p->err && t->t == F_LDIFF) {
        if (!null_node(m) || !null_node(r)) {
            post_process_ldiff(ret, m, r, p);
        }
    }

 bailout:

#if EDEBUG
    fprintf(stderr, "eval (t->t = %03d, %s): returning NODE %s at %p, err %d\n",
            t->t, getsymb(t->t), ret == NULL ? "nil" : getsymb(ret->t),
            (void *) ret, p->err);
    if (t->t == SERIES) {
        fprintf(stderr, " (SERIES node, xvec at %p, vnum = %d)\n",
                (void *) t->v.xvec, t->vnum);
    } else if (ret != NULL && ret->t == NUM) {
        fprintf(stderr, " (NUM node, return value %g)\n", ret->v.xval);
    }
#endif

    return ret;
}

#if !AUX_NODES_DEBUG

/* non-debugging variant: easier to see what's going on */

static inline int attach_aux_node (NODE *t, NODE *ret, parser *p)
{
    if (t->aux == NULL) {
        t->aux = ret;
        ret->refcount += 1;
    } else if (t->aux != ret) {
        if (t->t == QUERY || t->t == B_AND || t->t == B_OR) {
            /* OK to switch aux node in these cases */
            free_node(t->aux, p);
            t->aux = ret;
            ret->refcount += 1;
        } else if (is_proxy_node(t->aux)) {
            /* an extension to the above */
            free_node(t->aux, p);
            t->aux = ret;
            ret->refcount += 1;
        } else {
            /* otherwise if we're trying to switch aux node,
               something must have gone wrong
            */
            fprintf(stderr, "! %s node already has aux node attached\n",
                    getsymb(t->t));
            fprintf(stderr, " orig aux %p (%s), incoming %p (%s)\n",
                    (void *) t->aux, getsymb(t->aux->t),
                    (void *) ret, getsymb(ret->t));
            return E_DATA;
        }
    }

    return 0;
}

#else

/* variant with lots of debugging spew */

static inline int attach_aux_node (NODE *t, NODE *ret, parser *p)
{
    if (t->aux == NULL) {
        fprintf(stderr, "++ attach aux node %p (%s) to node %p (%s)\n",
                ret, getsymb(ret->t), t, getsymb(t->t));
        if (ret->refcount > 0) {
            fprintf(stderr, "   note: refcount on %p = %d\n",
                    (void *) ret, ret->refcount);
        }
        t->aux = ret;
        ret->refcount += 1;
    } else if (t->aux != ret) {
        if (t->t == QUERY || t->t == B_AND || t->t == B_OR) {
            /* the result node may switch in these cases (only?) */
            fprintf(stderr, "boolean: freeing %p (%s)\n",
                    (void *) t->aux, getsymb(t->aux->t));
            free_node(t->aux, p);
            fprintf(stderr, "boolean: attaching %p (%s)\n",
                    (void *) ret, getsymb(ret->t));
            t->aux = ret;
            ret->refcount += 1;
        } else if (is_proxy_node(t->aux)) {
            /* an extension to the above */
            fprintf(stderr, "proxy: freeing %p\n", (void *) t->aux);
            free_node(t->aux, p);
            fprintf(stderr, "proxy: attaching %p\n", (void *) ret);
            t->aux = ret;
            ret->refcount += 1;
        } else {
            fprintf(stderr, "!! node %p (%s) has aux node %p (%s) attached\n"
                    "   not attaching ret at %p (%s)\n", t, getsymb(t->t),
                    t->aux, getsymb(t->aux->t), ret, getsymb(ret->t));
            return E_DATA;
        }
    }

    return 0;
}

#endif

static int more_input (const char *s)
{
    while (*s) {
        if (!isspace((unsigned char) *s)) {
            return 1;
        }
        s++;
    }

    return 0;
}

/* get the next input character for the lexer */

int parser_getc (parser *p)
{
#if EDEBUG > 1
    fprintf(stderr, "parser_getc: src='%s'\n", p->point);
#endif

    p->ch = 0;

    if (more_input(p->point)) {
        p->ch = *p->point;
        p->point += 1;
    }

#if EDEBUG > 1
    if (p->ch) {
        fprintf(stderr, "parser_getc: returning '%c'\n", p->ch);
    }
#endif

    return p->ch;
}

/* advance the read position by n characters */

void parser_advance (parser *p, int n)
{
    p->point += n;
    p->ch = *p->point;
    p->point += 1;
}

/* throw back the last-read character */

void parser_ungetc (parser *p)
{
    p->point -= 1;
    p->ch = *(p->point - 1);
}

/* Look ahead for the first occurrence of a given character in
   the remaining input stream; return its 0-based index or
   -1 if not found.
*/

int parser_char_index (parser *p, int c)
{
    int i;

    for (i=0; p->point[i] != '\0'; i++) {
        if (p->point[i] == c) {
            return i;
        }
    }

    return -1;
}

/* For error reporting: print the input up to the current
   parse point, unless it's not valid UTF-8. Return 0
   if the input is printed OK, otherwise non-zero.
*/

int parser_print_input (parser *p)
{
    int len = p->point - p->input;
    char *s = gretl_strndup(p->input, len);
    int err = 0;

    if (s != NULL) {
        if (g_utf8_validate(s, -1, NULL)) {
            pprintf(p->prn, "> %s\n", s);
        } else {
            err = 1;
        }
        free(s);
    } else {
        err = 1;
    }

    return err;
}

/* "pretty print" syntactic nodes and symbols */

static void printsymb (int symb, const parser *p)
{
    pputs(p->prn, getsymb(symb));
}

static void printnode (NODE *t, parser *p, int value)
{
    if (t == NULL) {
        pputs(p->prn, "NULL");
    } else if (!value && useries_node(t)) {
        pprintf(p->prn, "%s", p->dset->varname[t->vnum]);
    } else if (!value && uscalar_node(t)) {
        pprintf(p->prn, "%s", t->vname);
    } else if (t->t == NUM) {
        if (na(t->v.xval)) {
            pputs(p->prn, "NA");
        } else {
            pprintf(p->prn, "%.8g", t->v.xval);
        }
    } else if (t->t == SERIES) {
        const double *x = t->v.xvec;
        int i, j = 1;

        if (p->lh.vnum > 0 && p->lh.vnum < p->dset->v) {
            pprintf(p->prn, "%s\n", p->dset->varname[p->lh.vnum]);
        }

        for (i=p->dset->t1; i<=p->dset->t2; i++, j++) {
            if (na(x[i])) {
                pputs(p->prn, "NA");
            } else {
                pprintf(p->prn, "%g", x[i]);
            }
            if (j % 8 == 0) {
                pputc(p->prn, '\n');
            } else if (i < p->dset->t2) {
                pputc(p->prn, ' ');
            }
        }
    } else if (t->t == MAT) {
        if (t->vname != NULL) {
            pputs(p->prn, t->vname);
        } else {
            gretl_matrix_print_to_prn(t->v.m, NULL, p->prn);
        }
    } else if (t->t == BUNDLE) {
        gretl_bundle_print(t->v.b, p->prn);
    } else if (t->t == DBUNDLE) {
        pputs(p->prn, bvarname(t->v.idnum));
    } else if (t->t == ARRAY) {
        gretl_array_print(t->v.a, p->prn);
    } else if (t->t == UOBJ) {
        pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == DVAR) {
        pputs(p->prn, dvarname(t->v.idnum));
    } else if (t->t == MVAR) {
        pputs(p->prn, mvarname(t->v.idnum));
    } else if (t->t == CON) {
        pputs(p->prn, constname(t->v.idnum));
    } else if (t->t == DUM) {
        pputs(p->prn, dumname(t->v.idnum));
    } else if (binary_op(t->t)) {
        pputc(p->prn, '(');
        printnode(t->L, p, 0);
        printsymb(t->t, p);
        printnode(t->R, p, 0);
        pputc(p->prn, ')');
    } else if (t->t == OSL) {
        printnode(t->L, p, 0);
        pputc(p->prn, '[');
        printnode(t->R, p, 0);
        pputc(p->prn, ']');
    } else if (t->t == SLRAW) {
        pputs(p->prn, "SLRAW");
    } else if (t->t == SUBSL) {
        pputs(p->prn, "SUBSL");
    } else if (func1_symb(t->t)) {
        printsymb(t->t, p);
        pputc(p->prn, '(');
        printnode(t->L, p, 0);
        pputc(p->prn, ')');
    } else if (unary_op(t->t)) {
        printsymb(t->t, p);
        printnode(t->L, p, 0);
    } else if (func2_symb(t->t)) {
        printsymb(t->t, p);
        pputc(p->prn, '(');
        printnode(t->L, p, 0);
        if (t->R->t != EMPTY) {
            pputc(p->prn, ',');
        }
        printnode(t->R, p, 0);
        pputc(p->prn, ')');
    } else if (t->t == STR) {
        pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == UFUN) {
        pprintf(p->prn, "%s", t->vname);
    } else if (t->t == DMSTR) {
        printnode(t->L, p, 0);
        pputc(p->prn, '(');
        printnode(t->R, p, 0);
        pputc(p->prn, ')');
    } else if (t->t == LISTVAR) {
        pprintf(p->prn, "%s.%s", t->L->v.str, t->R->v.str);
    } else if (t->t == LIST) {
        pputs(p->prn, "LIST");
    } else if (t->t == LAG) {
        pputs(p->prn, "LAG");
    } else if (t->t != EMPTY) {
        pputs(p->prn, "weird tree - ");
        printsymb(t->t, p);
    }
}

/* which modified assignment operators of the type '+='
   will we accept, when generating various types of
   result? */
#define ok_matrix_op(o) (o == B_ASN  || o == B_DOTASN || \
                         o == B_ADD  || o == B_SUB || \
                         o == B_MUL  || o == B_DIV || \
                         o == B_HCAT || o == B_VCAT)
#define ok_list_op(o) (o == B_ASN || o == B_ADD || o == B_SUB)
#define ok_string_op(o) (o == B_ASN || o == B_ADD || \
                         o == B_HCAT || o == INC)
#define ok_array_op(o) (o == B_ASN || o == B_ADD || o == B_SUB)
#define ok_bundle_op(o) (o == B_ASN || o == B_ADD)

struct mod_assign {
    int c;
    int op;
};

/* supported "inflections" of assignment */

struct mod_assign m_assign[] = {
    { '+', B_ADD },
    { '-', B_SUB },
    { '*', B_MUL },
    { '/', B_DIV },
    { '%', B_MOD},
    { '^', B_POW },
    { '~', B_HCAT },
    { '|', B_VCAT },
    { '.', B_DOTASN },
    { 0, 0}
};

/* read operator from formula: this is either
   simple assignment or something like '+=' */

static int get_op (char *s)
{
    if (s[0] == '=') {
        s[1] = '\0';
        return B_ASN;
    } else if (!strcmp(s, "++")) {
        return INC;
    } else if (!strcmp(s, "--")) {
        return DEC;
    } else if (s[1] == '=') {
        int i;

        for (i=0; m_assign[i].c; i++) {
            if (s[0] == m_assign[i].c) {
                return m_assign[i].op;
            }
        }
    }

    return 0;
}

static char *get_opstr (int op)
{
    static char opstr[4] = {0};

    if (op == B_ASN) {
        return "=";
    } else if (op == INC) {
        return "++";
    } else if (op == DEC) {
        return "--";
    } else {
        int i;

        for (i=0; m_assign[i].c; i++) {
            if (op == m_assign[i].op) {
                opstr[0] = m_assign[i].c;
                opstr[1] = '=';
                return opstr;
            }
        }
        return "??";
    }
}

/* implement the declaration of new variables */

static void do_declaration (parser *p)
{
    char **S = NULL;
    int i, v, n;

    n = check_declarations(&S, p);

    if (n == 0) {
        return;
    }

    for (i=0; i<n && !p->err; i++) {
        if (S[i] != NULL) {
            if (p->targ == SERIES) {
                p->err = dataset_add_NA_series(p->dset, 1);
                if (!p->err) {
                    v = p->dset->v - 1;
                    strcpy(p->dset->varname[v], S[i]);
                }
            } else {
                GretlType type = 0;

                if (p->targ == MAT) {
                    type = GRETL_TYPE_MATRIX;
                } else if (p->targ == NUM) {
                    type = GRETL_TYPE_DOUBLE;
                } else if (p->targ == STR) {
                    type = GRETL_TYPE_STRING;
                } else if (p->targ == BUNDLE) {
                    type = GRETL_TYPE_BUNDLE;
                } else if (p->targ == LIST) {
                    type = GRETL_TYPE_LIST;
                } else if (p->targ == ARRAY) {
                    type = p->lh.gtype;
                } else {
                    p->err = E_DATA;
                }
                if (!p->err) {
                    p->err = user_var_add(S[i], type, NULL);
                }
            }
        }
    }

    strings_array_free(S, n);
}

/* The expression supplied for evaluation does not contain an '=':
   can we interpret it as an implicit request to print the value
   of an existing variable?
*/

static void parser_try_print (parser *p, const char *s, int *done)
{
    if (p->lh.t != 0 && p->lh.expr == NULL) {
        p->flags |= P_DISCARD;
        p->point = s;
    } else {
        p->err = E_EQN;
    }
}

/* Here we try to parse out the LHS of the statement and also the
   operator. If we find a unitary LHS (simply an identifier) we write
   it into p->lh.name, but if we find a compound LHS (such as a
   sub-matrix specification) we save it as p->lh.expr. The content of
   @ps is advanced to the first position beyond the operator.
*/

static int extract_lhs_and_op (const char **ps, parser *p,
                               char *opstr)
{
    const char *s = *ps;
    int quoted = 0;
    int i, n = 0;
    int err = 0;

#if LHDEBUG
    fprintf(stderr, "extract_lhs_and_op: input='%s', targ=%s\n",
            s, getsymb(p->targ));
#endif

    if (p->targ != UNK && strchr(s, '=') == NULL) {
        /* we got a type specification but no assignment,
           so should be variable declaration(s) ?
        */
        p->flags |= P_DECL;
        p->lh.expr = gretl_strdup(s);
        goto done;
    }

    /* Count bytes preceding first unquoted '='. Note that
       the "unquoted" condition is required only because
       a string-literal bundle key might contain an equals
       sign, as in b["foo=bar"] = ...
    */
    for (i=0; s[i] != '\0'; i++) {
        if (s[i] == '"') {
            quoted = !quoted;
        } else if (!quoted && s[i] == '=') {
            break;
        }
        n++;
    }

    if (n > 0) {
        char *lhs = NULL;
        int lhlen = n;

        if (s[n] == '=') {
            /* we actually reached an '=' */
            if (s[n+1] == '=') {
                /* 2022-04-05: should this be an error? */
                goto done;
            } else if (strspn(s + n - 1, "+-*/%^~|.") == 1) {
                /* preceded by a modifier: inflected assignment */
                lhlen--;
                opstr[0] = s[n-1];
                opstr[1] = '=';
            } else {
                /* straight assignment */
                opstr[0] = '=';
            }
            n++; /* plus 1 for '=' */
        }

        if (lhlen > 0) {
            lhs = gretl_strndup(s, lhlen);
            tailstrip(lhs);
            lhlen = strlen(lhs);
        }

        if (opstr[0] == '\0' && lhlen > 2) {
            /* check for postfix operator */
            char *test = lhs + lhlen - 2;

            if (!strcmp(test, "++") || !strcmp(test, "--")) {
                strcpy(opstr, test);
                *test = '\0';
                lhlen -= 2;
            }
        }

        if (lhlen > 0) {
            if (lhlen == gretl_namechar_spn(lhs)) {
                /* a straight identifier? */
                if (lhlen >= VNAMELEN) {
                    pprintf(p->prn, _("'%s': name is too long (max %d characters)\n"),
                            lhs, VNAMELEN - 1);
                    err = E_PARSE;
                } else {
                    strcpy(p->lh.name, lhs);
                }
            } else if ((p->flags & P_PRIV) && (*lhs == '$' || *lhs == '_') &&
                       gretl_namechar_spn(lhs + 1) == lhlen - 1) {
                /* "private" genr of the form $foo=expr or _foo=expr */
                strcpy(p->lh.name, lhs);
            } else {
                /* treat as an expression to be evaluated */
                p->lh.expr = lhs;
                lhs = NULL; /* protect against freeing */
            }
        } else {
            /* nothing relevant found */
            err = E_PARSE;
        }

        if (!err && opstr[0] != '\0') {
            p->op = get_op(opstr);
        }

        free(lhs);
        *ps = s + n;
    }

 done:

#if LHDEBUG
    fprintf(stderr, "extract: lh.name='%s', lh.expr='%s', op='%s', err=%d, s='%s'\n",
            p->lh.name, p->lh.expr ? p->lh.expr : "NULL", opstr, err, *ps);
#endif

    if (!p->err && !(p->flags & P_DECL) && p->lh.name[0] == '\0' && p->op == 0) {
        /* added 2021-05-29 */
        p->flags |= P_DISCARD;
    }

    return err;
}

static void maybe_do_type_errmsg (const char *name, int t)
{
    const char *tstr = typestr(t);

    if (tstr != NULL && strcmp(tstr, "?")) {
        if (name != NULL && *name != '\0') {
            gretl_errmsg_sprintf(_("The variable %s is of type %s, "
                                   "not acceptable in context"),
                                 name, tstr);
        } else {
            gretl_errmsg_sprintf(_("A variable of type %s is not "
                                   "acceptable in context"), tstr);
        }
    }
}

static void assignment_type_errmsg (int targ, int rhs, int op)
{
    const char *rhstr = typestr(rhs);

    if (*rhstr == '?') {
        rhstr = getsymb(rhs);
    } else if (targ == LIST && rhs == SERIES) {
        rhstr = "anonymous series";
    }
    gretl_errmsg_sprintf(_("Incompatible types in assignment: "
                           "%s %s %s"), typestr(targ), get_opstr(op),
                         rhstr);
}

static int overwrite_type_check (parser *p)
{
    int err = 0;

    /* FIXME check for series/function collision here */

    if (p->targ != p->lh.t) {
        /* don't overwrite one type with another */
        maybe_do_type_errmsg(p->lh.name, p->lh.t);
        err = E_TYPES;
    }

    return err;
}

static int overwrite_const_check (const char *name, int vnum)
{
    if (object_is_const(name, vnum)) {
        return overwrite_err(name);
    } else {
        return 0;
    }
}

/* Check that we're not trying to modify a const object
   via a compound LHS expression; and while we're at
   it, check whether we should be generating a list
   (a list member of a bundle or an element of an array
   of lists).
*/

static int compound_const_check (NODE *lhs, parser *p)
{
    NODE *n = lhs;
    int i = 0, err = 0;

    if (n->t == BMEMB && n->R != NULL && n->R->t == STR) {
        GretlType t;

        t = gretl_bundle_get_member_type(n->L->v.b, n->R->v.str, NULL);
        if (t == GRETL_TYPE_LIST) {
            p->flags |= P_LISTDEF;
        }
    }

    while (n->t == OSL || n->t == OBS || n->t == BMEMB) {
        n = n->L;
        if (i == 0 && lhs->t == OSL && n->t == ARRAY) {
            if (gretl_array_get_type(n->v.a) == GRETL_TYPE_LISTS) {
                p->flags |= P_LISTDEF;
            }
        }
        i++;
    }

    /* do we have a const object at the tip of the tree? */
    if (n->vname != NULL) {
        err = overwrite_const_check(n->vname, n->vnum);
    }

    return err;
}

static int ok_array_decl (parser *p, const char *s)
{
    p->lh.gtype = 0;

    if (!strncmp(s, "strings ", 8)) {
        p->lh.gtype = GRETL_TYPE_STRINGS;
    } else if (!strncmp(s, "matrices ", 9)) {
        p->lh.gtype = GRETL_TYPE_MATRICES;
    } else if (!strncmp(s, "bundles ", 8)) {
        p->lh.gtype = GRETL_TYPE_BUNDLES;
    } else if (!strncmp(s, "lists ", 6)) {
        p->lh.gtype = GRETL_TYPE_LISTS;
    }

    return p->lh.gtype != 0;
}

/* Given an existing LHS variable, whose type is recorded in
   p->lh.t, check that the specified operator is supported
   for the type. Return error code if not.
*/

static int check_operator_validity (parser *p, const char *opstr)
{
    if (p->lh.t == MAT && !ok_matrix_op(p->op)) {
        /* matrices: we accept only a limited range of
           modified assignment operators */
        gretl_errmsg_sprintf(_("'%s' : not implemented for matrices"), opstr);
        return E_PARSE;
    } else if (p->lh.t == LIST && !ok_list_op(p->op)) {
        /* lists: same story as matrices */
        gretl_errmsg_sprintf(_("'%s' : not implemented for lists"), opstr);
        return E_PARSE;
    } else if (p->lh.t == STR && !ok_string_op(p->op)) {
        /* strings: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for strings"), opstr);
        return E_PARSE;
    } else if (p->lh.t == ARRAY && !ok_array_op(p->op)) {
        /* arrays: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for arrays"), opstr);
        return E_PARSE;
    } else if (p->lh.t == BUNDLE && !ok_bundle_op(p->op)) {
        /* bundles: ditto */
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
        return E_PARSE;
    } else if (p->lh.t != MAT && (p->op == B_VCAT || p->op == B_DOTASN)) {
        /* vertical concat: only OK for matrices */
        gretl_errmsg_sprintf(_("'%s' : only defined for matrices"), opstr);
        return E_PARSE;
    } else if (p->lh.t != MAT && p->lh.t != STR && p->op == B_HCAT) {
        /* horizontal concat: only OK for matrices, strings */
        gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
        return E_PARSE;
    }

    /* otherwise OK? */
    return 0;
}

/* Do we have an inline type specification preceding the
   statement proper? In most cases we shouldn't, since it
   will already have been handled by the tokenizer (and
   the type will now be recorded in p->targ). But we allow
   for finding a typespec here in case of "genrs" within
   nls/mle/gmm blocks, where the statement bypasses the
   regular tokenizer. (FIXME?)
*/

static void check_for_inline_typespec (const char **ps, parser *p)
{
    const char *s = *ps;

    if (!strncmp(s, "scalar ", 7)) {
        p->targ = NUM;
        s += 7;
    } else if (!strncmp(s, "series ", 7)) {
        p->targ = SERIES;
        s += 7;
    } else if (!strncmp(s, "matrix ", 7)) {
        p->targ = MAT;
        s += 7;
    } else if (!strncmp(s, "list ", 5)) {
        p->targ = LIST;
        s += 5;
    } else if (!strncmp(s, "string ", 7)) {
        p->targ = STR;
        s += 7;
    } else if (!strncmp(s, "bundle ", 7)) {
        p->targ = BUNDLE;
        s += 7;
    } else if (ok_array_decl(p, s)) {
        p->targ = ARRAY;
        s += strcspn(s, " ") + 1;
    }

    /* advance pointer */
    *ps = s;
}

/* Check @p->lh.name for the name of an existing series or
   user_var of some kind. If found, record the relevant
   info in p->lh.t, and also either p->lh.vnum (series) or
   p->lh.uv (other types).
*/

static int check_existing_lhs_type (parser *p, int *newvar)
{
    user_var *uvar;
    int v, err = 0;

    if ((err = gretl_reserved_word(p->lh.name))) {
        return err;
    }

    v = current_series_index(p->dset, p->lh.name);
    if (v >= 0) {
        p->lh.vnum = v;
        p->lh.t = SERIES;
        *newvar = 0;
        return 0;
    }

    uvar = get_user_var_by_name(p->lh.name);

    if (uvar != NULL) {
        GretlType vtype = uvar->type;

        p->lh.uv = uvar;
        *newvar = 0;

        if (vtype == GRETL_TYPE_MATRIX) {
            p->lh.t = MAT;
        } else if (vtype == GRETL_TYPE_DOUBLE) {
            if (uvar->flags & UV_NODECL) {
                if (p->targ == UNK) {
                    p->flags |= P_NODECL;
                } else {
                    uvar->flags &= ~UV_NODECL;
                }
            }
            p->lh.t = NUM;
        } else if (vtype == GRETL_TYPE_LIST) {
            p->lh.t = LIST;
        } else if (vtype == GRETL_TYPE_STRING) {
            p->lh.t = STR;
        } else if (vtype == GRETL_TYPE_BUNDLE) {
            p->lh.t = BUNDLE;
        } else if (vtype == GRETL_TYPE_ARRAY) {
            p->lh.gtype = gretl_array_get_type(uvar->ptr);
            p->lh.t = ARRAY;
        }
    }

    return err;
}

/* pre-process a "genr" statement */

static void gen_preprocess (parser *p, int flags, int *done)
{
    const char *s = p->input;
    char opstr[3] = {0};
    int newvar = 1;

    while (isspace(*s)) s++;

    /* skip leading command word, if any */
    if (!strncmp(s, "genr ", 5)) {
        s += 5;
    } else if (!strncmp(s, "print ", 6)) {
        /* allow this within (e.g.) mle block */
        p->flags |= P_DISCARD;
        s += 6;
    }

    while (isspace(*s)) s++;

    if (p->targ == UNK) {
        check_for_inline_typespec(&s, p);
    } else if (gretl_array_type(p->targ)) {
        /* record a plural type spec such as "matrices" under
           the "lh" member of @p.
        */
        p->lh.gtype = p->targ;
        p->targ = ARRAY;
    }

    /* check for types that cannot be generated in the
       absence of a dataset */
    if ((p->targ == SERIES || p->targ == LIST) &&
        (p->dset == NULL || p->dset_n == 0)) {
        no_data_error(p);
        return;
    }

    if (p->flags & P_DISCARD) {
        /* doing a simple "eval" */
        p->point = s;
        return;
    }

    /* extract LHS expression and operator, and test for a declaration */
    p->err = extract_lhs_and_op(&s, p, opstr);
    if (p->err || (p->flags & (P_DECL | P_DISCARD))) {
        return;
    }

    /* record next read position */
    p->point = s;

    if (p->lh.expr != NULL) {
        /* create syntax tree for the LHS expression */
        const char *savepoint = p->point;

        p->point = p->lh.expr;
        p->ch = parser_getc(p);
        lex(p);
        p->lhtree = expr(p);
        p->point = savepoint;
        p->ch = 0;
        if (!p->err) {
            p->err = compound_const_check(p->lhtree, p);
        }
        if (p->err) {
            return;
        } else {
            goto get_rhs;
        }
    }

    /* find out if the LHS var already exists, and if
       so, what type it is */
    if (!p->err) {
        p->err = check_existing_lhs_type(p, &newvar);
    }

#if LHDEBUG
    fprintf(stderr, "newvar=%d, err=%d\n", newvar, p->err);
#endif

    if (p->err) {
        return;
    }

    if (newvar) {
        /* new variable: check name for legality */
        if (!(flags & P_PRIV)) {
            p->err = check_identifier(p->lh.name);
        }
    } else {
        /* pre-existing var: check for const-ness */
        p->err = overwrite_const_check(p->lh.name, p->lh.vnum);
    }

    if (p->err) {
        return;
    }

    if (p->lh.t != 0) {
        if (p->targ == UNK) {
            /* when a result type is not specified, set this
               from existing LHS variable, if present
            */
            p->targ = p->lh.t;
        } else if (overwrite_type_check(p)) {
            /* don't overwrite one type with another */
            p->err = E_TYPES;
            return;
        }
    }

 get_rhs:

    /* advance past white space */
    while (isspace(*s)) s++;
    p->point = p->rhs = s;

    if (p->lh.expr != NULL) {
        goto alt_set_targ;
    }

    /* expression ends here with no operator: a call to print? */
    if (*s == '\0' && p->op == 0) {
        parser_try_print(p, p->lh.name, done);
        return;
    }

    /* if the LHS variable does not already exist, then
       we can't do '+=' or anything of that sort, only
       simple assignment, B_ASN
    */
    if (newvar && p->op != B_ASN) {
        undefined_symbol_error(p->lh.name, p);
        return;
    }

    if (p->op) {
        p->err = check_operator_validity(p, opstr);
        if (p->err) {
            return;
        }
    }

 alt_set_targ:

    if (p->targ == UNK && *p->rhs == '{') {
        /* if the target type is still unknown and the RHS
           expression is wrapped in '{' and '}', make the target
           a matrix
        */
        p->targ = MAT;
    } else if (p->targ == LIST) {
        /* flag presence of list target to parser */
        p->flags |= P_LISTDEF;
    }
}

/* tests for saving variable */

static int matrix_may_be_masked (const gretl_matrix *m, int n,
                                 parser *p)
{
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int fullrows = mt2 - mt1 + 1;
    int nobs = get_matrix_mask_nobs();

    if (n == nobs && fullrows > n) {
        p->flags |= P_MMASK;
        return 1;
    } else {
        return 0;
    }
}

/* check whether a matrix result can be assigned to a series
   on return */

static int series_compatible_matrix (const gretl_matrix *m, parser *p)
{
    int n = gretl_vector_get_length(m);
    int mt2 = gretl_matrix_get_t2(m);
    int T = sample_size(p->dset);
    int ok = 0;

    if (mt2 > 0) {
        int mt1 = gretl_matrix_get_t1(m);

        if (n == mt2 - mt1 + 1) {
            /* sample is recorded on matrix */
            ok = 1;
        } else if (matrix_may_be_masked(m, n, p)) {
            ok = 1;
        }
    } else if (n == T) {
        /* length matches current sample */
        ok = 1;
    } else if (n == p->dset->n) {
        /* length matches full series length */
        ok = 1;
    } else if (n == 1) {
        /* scalar: can be expanded */
        ok = 1;
    }

    return ok;
}

/* Check for assignment of an array of strings to a series,
   which thereby becomes string-valued. The size of the array
   must equal the full length of the dataset, which must not
   be sub-sampled, and if there's a pre-existing series it
   must not be string-valued already.
*/

static int series_compatible_array (gretl_array *a, parser *p)
{
    if (gretl_array_get_type(a) != GRETL_TYPE_STRINGS) {
        return 0;
    } else if (gretl_array_get_length(a) != p->dset->n) {
        return 0;
    } else if (dataset_is_subsampled(p->dset)) {
        return 0;
    } else if (is_string_valued(p->dset, p->lh.vnum)) {
        return 0;
    } else {
        /* all tests passed */
        return 1;
    }
}

/* This function converts a series to a column vector,
   respecting the value of the set-variable "skip_missing".
   In that respect it differs from tmp_matrix_from_series(),
   which always just grabs the entire sample range.
*/

static gretl_matrix *series_to_matrix (const double *x,
                                       parser *p)
{
    int t, n = sample_size(p->dset);
    int t1 = p->dset->t1;
    int t2 = p->dset->t2;
    int skip_na, contiguous = 1;
    gretl_matrix *v;

    skip_na = libset_get_bool(SKIP_MISSING);

    if (skip_na) {
        int err = series_adjust_sample(x, &t1, &t2);

        if (!err) {
            /* no interior NAs, just use (possibly) revised sample */
            n = t2 - t1 + 1;
        } else {
            /* we have to count the non-missing values */
            n = 0;
            for (t=t1; t<=t2; t++) {
                if (!na(x[t])) n++;
            }
            /* the values we want are not contiguous */
            contiguous = 0;
        }
    }

    if (n == 0) {
        v = gretl_null_matrix_new();
    } else {
        v = gretl_column_vector_alloc(n);
    }

    if (v == NULL) {
        p->err = E_ALLOC;
    } else if (n > 0) {
        if (contiguous) {
            memcpy(v->val, x + t1, n * sizeof *x);
        } else {
            int i = 0;

            for (t=t1; t<=t2; t++) {
                if (na(x[t])) {
                    if (!skip_na) {
                        v->val[i++] = x[t];
                    }
                } else {
                    v->val[i++] = x[t];
                }
            }
        }
        if (contiguous) {
            gretl_matrix_set_t1(v, t1);
            gretl_matrix_set_t2(v, t2);
        }
    }

    return v;
}

static gretl_matrix *retrieve_matrix_result (parser *p)
{
    NODE *r = p->ret;
    gretl_matrix *m = NULL;

#if EDEBUG
    fprintf(stderr, "retrieve_matrix_result: r->t = %d\n", r->t);
#endif

    if (r->t == NUM) {
        m = gretl_matrix_from_scalar(r->v.xval);
        if (m == NULL) {
            p->err = E_ALLOC;
        } else if (na(r->v.xval)) {
            set_gretl_warning(W_GENNAN);
        }
    } else if (r->t == SERIES) {
        m = series_to_matrix(r->v.xvec, p);
    } else if (r->t == LIST) {
        m = gretl_list_to_vector(r->v.ivec, &p->err);
    } else if (r->t == MAT && is_tmp_node(r)) {
        /* result matrix is newly allocated, steal it */
#if EDEBUG
        fprintf(stderr, "matrix result (%p) is tmp, stealing it\n",
                (void *) r->v.m);
#endif
        m = r->v.m;
        r->v.m = NULL; /* avoid double-freeing */
    } else if (r->t == MAT) {
        /* r->v.m is an existing user matrix (or bundled matrix):
           must make a copy to keep pointers distinct
        */
        m = gretl_matrix_copy(r->v.m);
#if EDEBUG
        fprintf(stderr, "matrix result (%p) is pre-existing, copied to %p\n",
                (void *) r->v.m, (void *) m);
#endif
        if (m == NULL) {
            p->err = E_ALLOC;
        }
    } else if (null_node(r)) {
        m = gretl_null_matrix_new();
    } else {
        fprintf(stderr, "Looking for matrix, but r->t = %s\n", getsymb(r->t));
        p->err = E_TYPES;
    }

    return m;
}

/* Check to see if the existing LHS matrix is of the
   same dimensions as the RHS result */

static int LHS_matrix_reusable (parser *p, gretl_matrix **pm,
                                gretl_matrix *tmp)
{
    gretl_matrix *m = gen_get_lhs_var(p, GRETL_TYPE_MATRIX);
    int ok = 0;

    if (m == NULL) {
        return 0;
    } else if (p->ret->t == NUM) {
        ok = (m->rows == 1 && m->cols == 1);
    } else if (p->ret->t == SERIES) {
        ok = (m->rows == tmp->rows && m->cols == 1);
    } else if (p->ret->t == MAT) {
        gretl_matrix *retm = p->ret->v.m;

        ok = (retm != NULL &&
              m->rows == retm->rows &&
              m->cols == retm->cols &&
              m->is_complex == retm->is_complex);
    }

    *pm = m;

    return ok;
}

/* Generating a matrix, and there's a pre-existing LHS matrix:
   we re-use the left-hand side matrix if possible.
*/

static gretl_matrix *assign_to_matrix (parser *p)
{
    gretl_matrix *m = NULL;
    gretl_matrix *tmp = NULL;
    int free_tmp = 1;
    double x;

    if (p->ret->t == SERIES) {
        /* a legacy thing */
        tmp = series_to_matrix(p->ret->v.xvec, p);
        if (p->err) {
            return NULL;
        }
    }

    if (LHS_matrix_reusable(p, &m, tmp)) {
        /* The result is of the same dimensions as the LHS matrix:
           this means that we don't need to construct an RHS
           matrix if it doesn't already exist as such, nor do we
           need to copy it if it does already exist.
        */
#if EDEBUG
        fprintf(stderr, "assign_to_matrix: reusing LHS\n");
#endif
        if (p->ret->t == NUM) {
            /* using RHS scalar */
            m->val[0] = x = p->ret->v.xval;
        } else if (p->ret->t == SERIES) {
            /* using RHS series, converted to @tmp */
            p->err = gretl_matrix_copy_data(m, tmp);
        } else {
            /* using RHS matrix: just copy data across */
            p->err = gretl_matrix_copy_data(m, p->ret->v.m);
        }
    } else {
        /* Dimensions differ: replace the LHS matrix */
#if EDEBUG
        fprintf(stderr, "assign_to_matrix: replacing\n");
#endif
        if (tmp != NULL) {
            p->err = gen_replace_lhs(p, GRETL_TYPE_MATRIX, tmp);
            free_tmp = 0; /* @tmp is the return value */
        } else {
            m = retrieve_matrix_result(p);
            if (!p->err) {
                p->err = gen_replace_lhs(p, GRETL_TYPE_MATRIX, m);
            }
        }
    }

    if (tmp != NULL && free_tmp) {
        gretl_matrix_free(tmp);
    }

    return m;
}

/* Assigning to an existing (whole) LHS matrix, but using '+='
   or some such modified/inflected assignment. Note that
   save_generated_var() is the only caller.
*/

static gretl_matrix *assign_to_matrix_mod (gretl_matrix *m1,
                                           parser *p)
{
    gretl_matrix *m2 = NULL;
    int mcat;

    if (m1 == NULL) {
        p->err = E_DATA;
    }

    /* In most cases we can take a shortcut when the RHS
       value is scalar, but we can't do that when the
       inflection is one of the matrix concatenation
       operators: here we record that fact.
    */
    mcat = (p->op == B_HCAT || p->op == B_VCAT);

    if (!p->err) {
        if (p->op == B_DOTASN) {
            p->err = dot_assign_to_matrix(m1, p);
            m2 = m1; /* no change in matrix pointer */
        } else if (!mcat && scalar_node(p->ret)) {
            double x = node_get_scalar(p->ret, p);

            if (m1->is_complex) {
                cmatrix_xy_calc(m1, m1, x, 0, p->op, p);
            } else {
                rmatrix_xy_calc(m1, m1, x, 0, p->op, p);
            }
            m2 = m1; /* no change in matrix pointer */
        } else if (!mcat && cscalar_node(p->ret)) {
            if (m1->is_complex) {
                double complex z = p->ret->v.m->z[0];

                cmatrix_xy_calc(m1, m1, z, 0, p->op, p);
                m2 = m1; /* no change in matrix pointer */
            } else {
                p->err = E_TYPES;
            }
        } else {
            gretl_matrix *tmp = retrieve_matrix_result(p);

            if (tmp != NULL) {
                p->err = real_matrix_calc(m1, tmp, p->op, &m2);
                gretl_matrix_free(tmp);
            }
        }
    }

    return m2;
}

static int should_copy_into_array (NODE *n, parser *p)
{
    int donate = is_tmp_node(n);

    /* In general if @n is a tmp node we can 'donate' its content
       to an array rather than copying it in, but there's a tricky
       case: we're compiling a "genr", or running a compiled one,
       and @n holds a string literal. If we don't copy this, it
       will have disappeared by the next call to the genr. There
       may be a simpler and more elegant way of identifying this
       case, but for now the following works (2022-08-10).
    */
    if (donate && !is_aux_node(n) && n->t == STR &&
        (p->flags & (P_COMPILE | P_EXEC))) {
        donate = 0; /* we _should_ copy! */
    }

    return !donate;
}

static void do_array_append (parser *p)
{
    gretl_array *A = NULL;
    GretlType atype;
    NODE *rhs = p->ret;
    void *ptr = NULL;

    A = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);
    if (A == NULL) {
        p->err = E_DATA;
        return;
    }

    atype = gretl_array_get_content_type(A);

    if (atype == GRETL_TYPE_STRING && rhs->t == STR) {
        ptr = rhs->v.str;
    } else if (atype == GRETL_TYPE_MATRIX && rhs->t == MAT) {
        ptr = rhs->v.m;
    } else if (atype == GRETL_TYPE_BUNDLE && rhs->t == BUNDLE) {
        ptr = rhs->v.b;
    } else if (atype == GRETL_TYPE_LIST && rhs->t == LIST) {
        ptr = rhs->v.ivec;
    } else if (atype == GRETL_TYPE_ARRAY && rhs->t == ARRAY) {
        ptr = rhs->v.a;
    } else if (rhs->t == ARRAY) {
        /* special: not actually appending an _element_;
           stick rhs array onto end of lhs array
        */
        p->err = gretl_array_copy_into(A, rhs->v.a);
    } else {
        p->err = E_TYPES;
    }

#if EDEBUG
    fprintf(stderr, "do_array_append: node to append\n");
    print_node_details(rhs);
#endif

    if (!p->err && ptr != NULL) {
        int copy = should_copy_into_array(rhs, p);

        p->err = gretl_array_append_object(A, ptr, copy);
        if (!copy && !p->err) {
            rhs->v.ptr = NULL;
        }
    }
}

static void do_array_subtract (parser *p)
{
    gretl_array *A;
    NODE *rhs = p->ret;

    A = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);
    if (A == NULL) {
        p->err = E_DATA;
    } else if (gretl_array_get_type(A) == GRETL_TYPE_STRINGS && rhs->t == STR) {
        p->err = gretl_array_drop_string(A, rhs->v.str);
    } else if (rhs->t == EMPTY) {
        p->err = gretl_array_drop_null(A);
    } else {
        p->err = E_TYPES;
    }
}

static void do_bundle_append (parser *p)
{
    gretl_bundle *bl = NULL;
    gretl_bundle *br = NULL;
    NODE *rhs = p->ret;

    bl = gen_get_lhs_var(p, GRETL_TYPE_BUNDLE);
    if (rhs->t == BUNDLE) {
        br = rhs->v.b;
    }
    if (bl == NULL || br == NULL) {
        p->err = E_TYPES;
    } else {
        p->err = gretl_bundle_append(bl, br);
    }
}

static int create_or_edit_string (parser *p)
{
    const char *src = NULL;
    const char *orig = NULL;
    char *newstr = NULL;
    user_var *uvar;

    if (p->ret->t == NUM) {
        /* OK only in case of "+=" */
        if (p->op != B_ADD) {
            p->err = E_TYPES;
            return p->err;
        }
    } else if (p->ret->t == LIST) {
        /* string representation of list: simple assign only */
        if (p->op != B_ASN) {
            p->err = E_TYPES;
        } else {
            newstr = gretl_list_to_compact_string(p->ret->v.ivec,
                                                  p->dset, 0, &p->err);
        }
        if (p->err) {
            return p->err;
        }
        src = g_strchug(newstr);
    } else if (null_node(p->ret) || p->ret->v.str == NULL) {
        src = "";
    } else {
        src = p->ret->v.str;
    }

#if EDEBUG
    fprintf(stderr, "edit_string: src='%s'\n", src);
#endif

    uvar = p->lh.uv;

    if (uvar != NULL) {
        orig = uvar->ptr;
    } else if (p->op != B_ASN) {
        /* without an existing LHS string we can only assign */
        p->err = E_DATA;
        return p->err;
    }

    if (p->ret->t == NUM) {
        /* taking an offset into an existing string */
        int len = g_utf8_strlen(orig, -1);
        int adj = p->ret->v.xval;

        if (adj < 0) {
            p->err = E_DATA;
        } else if (adj == 0) {
            ; /* no-op */
        } else {
            if (adj < len) {
                src = g_utf8_offset_to_pointer(orig, adj);
            } else {
                src = "";
            }
            newstr = gretl_strdup(src);
            if (newstr == NULL) {
                p->err = E_ALLOC;
            } else {
                gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
            }
        }
    } else if (src == NULL) {
        ; /* no-op -- e.g. argname() didn't get anything */
    } else if (p->op == B_ASN) {
        /* simple assignment */
        if (newstr == NULL) {
            newstr = gretl_strdup(src);
            if (newstr == NULL) {
                p->err = E_ALLOC;
            }
        }
        if (!p->err) {
            if (uvar == NULL) {
                gen_add_uvar(p, GRETL_TYPE_STRING, newstr);
            } else {
                gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
            }
        }
    } else if (p->op == B_HCAT || p->op == B_ADD) {
        /* string concatenation */
        if (*src == '\0') {
            ; /* no-op */
        } else {
            newstr = malloc(strlen(orig) + strlen(src) + 1);
            if (newstr == NULL) {
                p->err = E_ALLOC;
            } else {
                strcpy(newstr, orig);
                strcat(newstr, src);
                gen_replace_lhs(p, GRETL_TYPE_STRING, newstr);
            }
        }
    }

    return p->err;
}

static int create_or_edit_list (parser *p)
{
    int *list = NULL;

    if (p->ret->t == MAT && gretl_vector_get_length(p->ret->v.m) == 0) {
        /* special case, list from matrix */
        const char *prefix;

        prefix = p->ret->vname != NULL ? p->ret->vname : p->lh.name;
        list = gretl_list_from_matrix(p->ret->v.m, prefix,
                                      p->dset, &p->err);
    } else {
        list = node_get_list(p->ret, p); /* note: copied */
    }

#if EDEBUG
    printlist(list, "RHS list in edit_list()");
#endif

#if 0 /* we're not applying the following check (yet?) */
    if (gretl_function_depth() > 0) {
        int i, vi;

        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi < 0 || vi >= p->dset->v) {
                /* this error will be caught below */
                break;
            }
            if (!series_is_accessible_in_function(vi, p->dset)) {
                p->err = E_DATA;
                break;
            }
        }
    }
#endif

    if (!p->err) {
        if (p->lh.t != LIST) {
            /* no pre-existing LHS list: must be simple assignment */
            p->err = remember_list(list, p->lh.name, NULL);
        } else if (p->op == B_ASN || p->op == B_ADD || p->op == B_SUB) {
            /* replace, append or subtract list members */
            p->err = gen_edit_list(p, list, p->op);
        } else {
            p->err = E_TYPES;
        }
    }

#if 0
    if (!p->err) {
        /* 2020-05-29: is this right, for list? */
        set_dataset_is_changed(p->dset, 1);
    }
#endif

    free(list);

    return p->err;
}

#define ok_return_type(t) (t == NUM || t == SERIES || t == MAT ||       \
                           t == LIST || t == DUM || t == EMPTY ||       \
                           t == STR || t == BUNDLE || t == ARRAY ||     \
                           t == U_ADDR || t == DBUNDLE)

/* Note: we're doing this only in relation to "primary" types
   (excluding bundle members, array elements, matrix sub-
   specs).
*/

static int gen_check_return_type (parser *p)
{
    NODE *r = p->ret;
    int err = 0;

    if (r == NULL) {
        fprintf(stderr, "gen_check_return_type: p->ret = NULL!\n");
        return E_DATA;
    }

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: targ=%s; ret at %p, type %s\n",
            getsymb(p->targ), (void *) r, getsymb(r->t));
#endif

    if (!ok_return_type(r->t)) {
        return E_TYPES;
    }

    if (r->t == SERIES && r->v.xvec == NULL) {
        fprintf(stderr, "got SERIES return with xvec = NULL!\n");
        return E_DATA;
    }

    if (p->targ == NUM) {
        if (r->t == NUM || scalar_matrix_node(r)) {
            ; /* scalar or 1 x 1 matrix: OK */
        } else if (r->t == MAT && (p->flags & P_NODECL)) {
            ; /* morphing to matrix may be OK */
        } else {
            err = E_TYPES;
        }
    } else if (p->targ == SERIES) {
        /* result must be scalar, series, or conformable matrix
           or array of strings
        */
        if (r->t == NUM || r->t == SERIES) {
            ; /* OK */
        } else if (r->t == MAT) {
            if (!series_compatible_matrix(r->v.m, p)) {
                err = E_TYPES;
            }
        } else if (r->t == ARRAY) {
            if (p->op != B_ASN || !series_compatible_array(r->v.a, p)) {
                err = E_TYPES;
            }
        } else {
            err = E_TYPES;
        }
    } else if (p->targ == MAT) {
        ; /* no-op: handled later */
    } else if (p->targ == LIST) {
        if (r->t != EMPTY && r->t != MAT && !ok_list_node(r, p)) {
            err = E_TYPES;
        }
    } else if (p->targ == STR) {
        if (r->t != EMPTY && r->t != STR && r->t != NUM && r->t != LIST) {
            err = E_TYPES;
        }
    } else if (p->targ == BUNDLE) {
        if (p->op == B_ASN) {
            /* plain assignment: bundle or null */
            if (r->t != BUNDLE && r->t != DBUNDLE && r->t != EMPTY) {
                err = E_TYPES;
            }
        } else {
            /* the only other assignment possibility is "+=",
               in which case we'll only accept a bundle
            */
            if (r->t != BUNDLE) {
                err = E_TYPES;
            }
        }
    } else if (p->targ == ARRAY) {
        if (p->op == B_ASN) {
            /* plain assignment: array, null, or suitable series */
            if (gen_type_is_arrayable(r->t) || r->t == EMPTY) {
                ; /* OK */
            } else if (strvals_node(r, p)) {
                ; /* possibly OK */
            } else {
                err = E_TYPES;
            }
        } else {
            /* arrays: the only other assignment possibility is "+=",
               in which case we'll only accept an array or an
               object which matches the content type on the left
               (but the matching check is deferred)
            */
            if (!gen_type_is_arrayable(r->t) && r->t != EMPTY) {
                err = E_TYPES;
            }
        }
    }

    if (err == E_TYPES) {
        assignment_type_errmsg(p->targ, r->t, p->op);
    }

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: returning with p->err = %d\n",
            err);
#endif

    return err;
}

/* Allocate storage if saving a new series to the dataset. We also
   check for collision with the name of a function and issue a warning
   if need be.
*/

static int gen_allocate_storage (parser *p)
{
    if (p->dset == NULL || p->dset->Z == NULL) {
        p->err = E_DATA;
    } else {
        p->err = dataset_add_NA_series(p->dset, 1);
        if (!p->err) {
            p->lh.vnum = p->dset->v - 1;
        }
    }
    if (!p->err && gretl_function_depth() == 0 &&
        get_user_function_by_name(p->lh.name) != NULL) {
        gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name"),
                              p->lh.name);
    } else if (!p->err && function_lookup(p->lh.name)) {
        gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name"),
                              p->lh.name);
    }

    return p->err;
}

static void align_matrix_to_series (double *y, const gretl_matrix *m,
                                    parser *p)
{
    const gretl_matrix *mask = get_matrix_mask();
    int t, s = 0;

    if (mask == NULL || mask->rows != p->dset->n) {
        p->err = E_DATA;
        return;
    }

    for (t=0; t<p->dset->n; t++) {
        if (mask->val[t] != 0.0) {
            if (t >= p->dset->t1 && t <= p->dset->t2) {
                y[t] = xy_calc(y[t], m->val[s], p->op, SERIES, p);
            }
            s++;
        }
    }
}

static int assign_null_to_bundle (parser *p)
{
    gretl_bundle *b;
    int err = 0;

    if (p->lh.t == BUNDLE) {
        b = gen_get_lhs_var(p, GRETL_TYPE_BUNDLE);
        gretl_bundle_void_content(b);
    } else {
        b = gretl_bundle_new();
        if (b == NULL) {
            err = E_ALLOC;
        } else {
            err = gen_add_uvar(p, GRETL_TYPE_BUNDLE, b);
        }
    }

    return err;
}

static int assign_null_to_array (parser *p)
{
    gretl_array *a;
    int err = 0;

    if (p->lh.t == ARRAY) {
        a = gen_get_lhs_var(p, GRETL_TYPE_ARRAY);
        gretl_array_void_content(a);
    } else {
        a = gretl_array_new(p->lh.gtype, 0, &err);
        if (!err) {
            err = gen_add_uvar(p, p->lh.gtype, a);
        }
    }

    return err;
}

static double get_lhs_scalar (parser *p)
{
    if (p->lh.uv == NULL) {
        p->lh.uv = get_user_var_by_name(p->lh.name);
    }
    if (p->lh.uv != NULL) {
        return uvar_get_scalar_value(p->lh.uv);
    } else {
        p->err = E_DATA;
        return NADBL;
    }
}

/* apply postfix '++' or '--' to LHS scalar, or '++' to
   LHS string (only) */

static int do_incr_decr (parser *p)
{
    if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_DOUBLE) {
        double x = get_lhs_scalar(p);

        if (!na(x)) {
            x += (p->op == INC)? 1 : -1;
            uvar_set_scalar_fast(p->lh.uv, x);
        }
    } else if (p->lh.uv != NULL && p->lh.uv->type == GRETL_TYPE_STRING) {
        if (p->op == DEC) {
            p->err = E_TYPES;
        } else {
            char *s = p->lh.uv->ptr;

            if (*s != '\0') {
                char *smod = gretl_strdup(s + 1);

                gen_replace_lhs(p, GRETL_TYPE_STRING, smod);
            }
        }
    } else {
        p->err = E_TYPES;
    }

    return p->err;
}

#define has_aux_mat(n) (n->aux != NULL && n->aux->t == MAT)

static int explore_node (NODE *t, int lev, NODE *prev,
                         parser *p)
{
    NODE pms = {0};
    int save_op;
    int err = 0;

#if LHDEBUG
    fprintf(stderr, "%d: %s %p, prev %p", lev, getsymb(t->t), (void *) t,
            (void *) prev);
    fprintf(stderr, " (aux %s)", (t->aux != NULL)? getsymb(t->aux->t) : "null");
    if (t->R != NULL) {
        fprintf(stderr, ", R %s", getsymb(t->R->t));
        if (t->R->aux != NULL) {
            fprintf(stderr, " (aux %s)", getsymb(t->R->aux->t));
            if (t->R->aux->t == MSPEC) {
                fputc('\n', stderr);
                print_mspec(t->R->aux->v.mspec);
            }
        }
    }
    if (t->t == MAT) {
        fputc('\n', stderr);
        gretl_matrix_print(t->v.m, "t->v.m");
    } else if (t->aux != NULL && t->aux->t == MAT) {
        gretl_matrix_print(t->aux->v.m, "t->aux->v.m");
    } else {
        fputc('\n', stderr);
    }
#endif
    if (prev != NULL && (t->t == MAT || has_aux_mat(t))) {
#if LHDEBUG
        fprintf(stderr, "doing ASSIGN to %s\n\n",
                t->t == MAT ? "MAT" : "aux MAT");
#endif
        pms.t = OSL;
        /* pms.L: node holding target matrix */
        pms.L = t->t == MAT ? t : t->aux;
        /* pms.R: node holding mspec */
        pms.R = prev->R->aux;
        save_op = p->op;
        p->op = B_ASN;
        /* prev->aux holds the replacement matrix */
        err = set_matrix_chunk(&pms, prev->aux, p);
        p->op = save_op;
    }

    return err;
}

static int traverse_left (parser *p)
{
    NODE *t = p->lhtree;
    NODE *prev = NULL;
    int level = 0;
    int err = 0;

    while (t && !err) {
        err = explore_node(t, level, prev, p);
        if (t->aux != NULL && t->aux->t == MAT) {
            prev = t;
        } else {
            prev = NULL;
        }
        t = t->L;
        level++;
    }

    return err;
}

/* set_nested_matrix_value(): this and its helper above,
   traverse_left(), require a little comment.

   We come here when a hansl statement modifies a matrix that
   is "under" something else. That something could be a
   bundle or array, as in

   # Case 0
   b.m[diag] = x     # under bundle b
   a[3][1:2,1:2] = y # under array a

   In such cases the first invocation of set_matrix_chunk
   below is sufficient. However, we also come here when the
   matrix is "under" another matrix -- that is, we have a
   double index or subspec, as in these examples for a complex
   matrix, C:

   # Case 1
   C[real][1:2,1:2] = x
   C[3,3][real] = y

   and also in these more extended examples where the
   complex matrix is itself "under" something else:

   # Case 2
   b.C[real][1:2,1:2] = x
   a[3][i,j][imag] = y

   To handle such cases we have to crawl the parser's
   "lhtree" (left-hand side tree) to find the matrix that
   ultimately has to be modified, executing further calls
   to set_matrix_chunk(). In Case 1 above, the matrix
   we're looking for will be at depth 1 in the lhtree,
   while in Case 2 it will be at depth 2.
*/

static int set_nested_matrix_value (NODE *lhs,
                                    NODE *rhs,
                                    parser *p)
{
#if LHDEBUG
    int err;
    gretl_matrix_print(lhs->L->v.m, "L->v.m, before set matrix chunk");
    err = set_matrix_chunk(lhs, rhs, p);
    gretl_matrix_print(lhs->L->v.m, "L->v.m, after set matrix chunk");
#else
    int err = set_matrix_chunk(lhs, rhs, p);
#endif

    if (!err) {
        err = traverse_left(p);
    }

    return err;
}

/* Support direct assignment from an array of strings to a series
   (which will become string-valued). The source array must hold
   a number of strings equal to the length of the dataset.
*/

static int series_from_strings_array (DATASET *dset, int v,
                                      gretl_array *a)
{
    char **S;
    int ns;

    S = gretl_array_get_strings(a, &ns);
    return series_from_strings(dset, v, S, ns);
}

static void series_from_matrix (double *y, const gretl_matrix *m,
                                parser *p)
{
    int k = gretl_vector_get_length(m);
    int mt1 = gretl_matrix_get_t1(m);
    int mt2 = gretl_matrix_get_t2(m);
    int s, t;

    if (p->flags & P_MMASK) {
        /* result needs special alignment */
        align_matrix_to_series(y, m, p);
    } else if (k == 1) {
        /* result is effectively a scalar */
        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            y[t] = xy_calc(y[t], m->val[0], p->op, SERIES, p);
        }
    } else if (k == p->dset->n) {
        /* treat result as full-length series */
        for (t=p->dset->t1; t<=p->dset->t2; t++) {
            y[t] = xy_calc(y[t], m->val[t], p->op, SERIES, p);
        }
    } else if (k == sample_size(p->dset)) {
        /* treat as series of current sample length */
        for (t=p->dset->t1, s=0; t<=p->dset->t2; t++, s++) {
            y[t] = xy_calc(y[t], m->val[s], p->op, SERIES, p);
        }
    } else if (mt2 > 0 && mt2 >= mt1) {
        /* align using matrix "t1" value */
        for (t=mt1; t<mt1 + k && t<=p->dset->t2; t++) {
            if (t >= p->dset->t1) {
                y[t] = xy_calc(y[t], m->val[t-mt1], p->op,
                                  SERIES, p);
            }
        }
    }
}

enum {
    OVW_STRINGS = 1,
    OVW_NUMERIC
};

/* Determine if it seems OK to overwrite a string-valued series.
   This action was banned altogether up to gretl 2023b, but
   we're now allowing it subject to two conditions:

   * It must be a case of simple (not inflected) assignment
   * The RHS must be a series

   If the RHS series is not string-valued we'll treat that as
   OK provisionally. The suitability of the RHS values will be
   checked later, when we call assign_numeric_to_strvar():
   they must be compatible with the encoding on the LHS.

   If the case is not rejected as unsuitable, we return either
   OVW_STRINGS to indicate a string-valued RHS or OVW_NUMERIC
   to indicate an overwrite using plain numerical input.
*/

static int strv_overwrite_ok (parser *p)
{
    int lv = p->lh.vnum;   /* LHS series ID */
    int rv = p->ret->vnum; /* RHS series ID */
    int ok = 1;
    int ret = 0;

    if (p->op != B_ASN) {
        /* must be a case of straight assignment */
        ok = 0;
    } else if (p->ret->t != SERIES) {
        /* RHS must be a series */
        ok = 0;
    } else if (lv == rv) {
        /* and not identical to LHS */
        ok = 0;
    } else if (is_string_valued(p->dset, rv)) {
        ret = OVW_STRINGS;
    } else if (p->lh.stab != NULL) {
        ret = OVW_STRINGS;
    } else {
        ret = OVW_NUMERIC;
    }

    if (!ok) {
        gretl_errmsg_set("Cannot overwrite string-valued series in this context");
        p->err = E_TYPES;
    }

    return ret;
}

static inline int savegen_retval (int err)
{
#if EDEBUG
    fprintf(stderr, "save_generated_var: returning %d\n", err);
#endif
    return err;
}

#if ONE_BY_ONE_CAST

static void handle_one_by_one_cast (parser *p, NODE *r,
                                    int *no_decl)
{
    if (p->targ == UNK) {
        if (scalar_matrix_node(r)) {
            /* "cast" 1 x 1 matrix to scalar */
            *no_decl = 1;
            p->targ = NUM;
            p->flags |= P_NODECL;
        } else {
            p->targ = r->t;
        }
    } else if (p->targ == NUM && r->t == MAT && (p->flags & P_NODECL)) {
        /* We're looking at a @targ that was previously
           set to NUM by the "auto-cast" mechanism: allow
           it to morph to matrix if need be.
        */
        if (scalar_matrix_node(r)) {
            ; /* not a problem */
        } else if (p->lh.t == 0) {
            /* no pre-existing scalar var */
            p->targ = MAT;
        } else if (p->lh.t == NUM) {
            /* type-convert existing scalar */
            p->err = gretl_scalar_convert_to_matrix(p->lh.uv);
            if (!p->err) {
                p->targ = MAT;
            }
        }
    }
}

#endif

static int handle_compound_target (parser *p)
{
    NODE *r = p->ret;
    int compound_t;

    p->lhtree->flags |= LHT_NODE;
    p->flags |= P_START;
#if LHDEBUG > 1
    fprintf(stderr, "\n*** lhtree before eval ***\n");
    print_tree(p->lhtree, p, 0, 0);
#endif
    p->lhres = eval(p->lhtree, p);
#if LHDEBUG > 1
    if (p->lhres != NULL) {
        fprintf(stderr, "\n*** lhres post-eval ***\n");
        print_tree(p->lhres, p, 0, 0);
        fprintf(stderr, "\n*** lhtree post-eval ***\n");
        print_tree(p->lhtree, p, 0, 0);
        fputc('\n', stderr);
    }
#endif
    if (p->err) {
        return savegen_retval(p->err);
    }

    compound_t = p->lhres->t;
#if LHDEBUG
    fprintf(stderr, "save_generated_var: type = %s\n",
            getsymb(compound_t));
#endif
    if (compound_t == BMEMB) {
        p->err = set_bundle_value(p->lhres, r, p);
    } else if (compound_t == OBS) {
        p->err = set_series_obs_value(p->lhres, r, p);
    } else if (compound_t == OSL) {
        NODE *lh1 = p->lhres->L;

#if LHDEBUG
        fprintf(stderr, "OSL save: lh1 type = %s\n", getsymb(lh1->t));
#endif
        if (p->targ == SERIES && lh1->t != LIST && lh1->t != BUNDLE) {
            /* a sub-object can be of series type only if the "parent"
               object is a list or bundle
            */
            p->err = E_TYPES;
        } else if (lh1->t == ARRAY) {
            p->err = set_array_value(p->lhres, r, p);
        } else if (lh1->t == LIST) {
            p->err = set_list_value(p->lhres, r, p);
        } else if (lh1->t == STR) {
            p->err = set_string_value(p->lhres, r, p);
        } else if (lh1->t == BUNDLE) {
            p->err = set_bundle_value(p->lhres, r, p);
        } else if (lh1->t == SERIES) {
            p->err = set_series_obs_value(p->lhres, r, p);
        } else if (lh1->t == MAT) {
            p->err = set_nested_matrix_value(p->lhres, r, p);
        } else {
            gretl_errmsg_set(_("Invalid left-hand side expression"));
            p->err = E_TYPES;
        }
    } else {
        gretl_errmsg_set(_("Invalid left-hand side expression"));
        p->err = E_TYPES;
    }

    return savegen_retval(p->err);
}

static int handle_overwrite_strvals (parser *p, int targ, NODE *r)
{
    series_table *st;
    int copy = 0;
    int err;

    if (r->vnum == NO_VNUM) {
        /* we're assigning from an on-the-fly series, whose
           temporary series_table is attached to @p->lh
        */
        st = p->lh.stab;
    } else {
        /* assigning from a string-valued member of the dataset,
           whose series_table must be preserved as is
        */
        st = series_get_string_table(p->dset, r->vnum);
        copy = 1;
    }

    err = assign_strings_to_strvar(p->dset, targ, r->v.xvec, st, copy);

    if (st == p->lh.stab) {
        /* @st has done its work, and must now be set to NULL */
        if (dataset_is_subsampled(p->dset)) {
            /* in the subsampled case @st is not "donated"
               to the LHS series, so we have to free it here
               to avoid leaking memory
            */
            series_table_destroy(p->lh.stab);
        }
        p->lh.stab = NULL;
    }

    return err;
}

static int save_generated_var (parser *p, PRN *prn)
{
    NODE *r = p->ret;
    double **Z = NULL;
    double x;
    int no_decl = 0;
    int strv_ovwrite = 0;
    int t, v = 0;

#if EDEBUG
    fprintf(stderr, "\nsave (%s): lhname='%s', callcount %d\n"
            "  lh.t=%s, targ=%s, no_decl=%d, r->t=%s, r=%p\n",
            p->lhtree != NULL ? "compound" : "unitary",
            p->lh.name, p->callcount, getsymb(p->lh.t),
            getsymb(p->targ), (p->flags & P_NODECL)? 1 : 0,
            (r == NULL)? "none" : getsymb(r->t), (void *) r);
#endif

    if (p->lh.t == SERIES && is_string_valued(p->dset, p->lh.vnum) &&
        p->lhtree == NULL) {
        strv_ovwrite = strv_overwrite_ok(p);
        if (!strv_ovwrite) {
            return savegen_retval(p->err);
        }
    }

    /* some special cases */
    if (p->lhtree != NULL) {
        return handle_compound_target(p);
    } else if (p->op == INC || p->op == DEC) {
        return do_incr_decr(p);
    }

    if (p->callcount < 2) {
        /* first exec: test for type mismatch errors */
        p->err = gen_check_return_type(p);
        if (p->err) {
            return savegen_retval(p->err);
        }
    }

#if ONE_BY_ONE_CAST
    handle_one_by_one_cast(p, r, &no_decl);
#else
    if (p->targ == UNK) {
        p->targ = r->t;
    }
#endif

#if EDEBUG
    fprintf(stderr, "after preliminaries: targ=%s, op='%s'\n",
            getsymb(p->targ), getsymb(p->op));
#endif

    if (p->targ == UNK || p->targ == EMPTY) {
        p->err = E_UNKVAR;
        return savegen_retval(p->err);
    }

    if (p->targ == SERIES && (unsigned char) p->lh.name[0] > 126) {
        /* can't allow Greek letters for series names */
        gretl_errmsg_sprintf("Invalid series name '%s'", p->lh.name);
        p->err = E_DATA;
        return savegen_retval(p->err);
    }

    /* allocate dataset storage, if needed */
    if (p->targ == SERIES && p->lh.vnum <= 0) {
        gen_allocate_storage(p);
        if (p->err) {
            return savegen_retval(p->err);
        }
    }

    if (p->dset != NULL && p->dset->Z != NULL) {
        /* convenience notation */
        Z = p->dset->Z;
        v = p->lh.vnum;
    }

    /* put the generated data into place */

    if (p->targ == NUM) {
        if (p->lh.t == NUM) {
            /* modifying an existing scalar */
            if (r->t == NUM) {
                x = r->v.xval;
            } else if (scalar_matrix_node(r)) {
                x = r->v.m->val[0];
            } else {
                p->err = E_TYPES;
            }
            if (!p->err && p->op != B_ASN) {
                double x0 = get_lhs_scalar(p);

                x = xy_calc(x0, x, p->op, NUM, p);
            }
            if (!p->err) {
                p->err = gen_replace_scalar(p, x);
            }
        } else {
            /* a new scalar */
            if (r->t == NUM) {
                x = r->v.xval;
            } else if (scalar_matrix_node(r)) {
                x = r->v.m->val[0];
            } else {
                p->err = E_TYPES;
            }
            if (!p->err) {
                if (no_decl) {
                    p->err = gretl_scalar_add_mutable(p->lh.name, x);
                } else {
                    p->err = gretl_scalar_add(p->lh.name, x);
                }
            }
        }
    } else if (p->targ == SERIES) {
        /* writing a series */
        if (strv_ovwrite == OVW_STRINGS) {
            p->err = handle_overwrite_strvals(p, v, r);
        } else if (strv_ovwrite == OVW_NUMERIC) {
            p->err = assign_numeric_to_strvar(p->dset, v, r->v.xvec);
        } else if (r->t == SERIES) {
            const double *x = r->v.xvec;

            if (p->op == B_ASN) {
                /* avoid multiple calls to xy_calc */
                if (Z[v] != x) {
                    size_t sz = sample_size(p->dset) * sizeof *x;

                    memcpy(Z[v] + p->dset->t1, x + p->dset->t1, sz);
                }
            } else {
                for (t=p->dset->t1; t<=p->dset->t2; t++) {
                    Z[v][t] = xy_calc(Z[v][t], x[t], p->op, SERIES, p);
                }
            }
        } else if (r->t == NUM) {
            for (t=p->dset->t1; t<=p->dset->t2; t++) {
                Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, SERIES, p);
            }
        } else if (r->t == MAT) {
            series_from_matrix(Z[v], r->v.m, p);
        } else if (r->t == ARRAY) {
            p->err = series_from_strings_array(p->dset, v, r->v.a);
            if (p->err) {
                return savegen_retval(p->err);
            }
        }
        strcpy(p->dset->varname[v], p->lh.name);
        series_unset_orig_pd(p->dset, v); /* 2020-09-27 */
#if EDEBUG
        fprintf(stderr, "var %d: gave generated series the name '%s'\n",
                v, p->lh.name);
        fprintf(stderr, " value[1] = %g\n", p->dset->Z[v][1]);
#endif
        if (!p->err) {
            /* (probably) changed or added a series */
            set_dataset_is_changed(p->dset, 1);
        }
    } else if (p->targ == MAT) {
        /* we're writing a matrix */
        gretl_matrix *m = NULL;

        if (p->lh.uv == NULL) {
            /* there's no pre-existing left-hand side matrix */
            m = retrieve_matrix_result(p);
            if (!p->err) {
                p->err = gen_add_uvar(p, GRETL_TYPE_MATRIX, m);
            }
        } else if (p->op == B_ASN) {
            /* uninflected assignment to an existing matrix */
            m = assign_to_matrix(p);
        } else {
            /* inflected assignment to entire existing matrix */
            gretl_matrix *m1 = gen_get_lhs_var(p, GRETL_TYPE_MATRIX);

            m = assign_to_matrix_mod(m1, p);
            if (!p->err) {
                p->err = gen_replace_lhs(p, GRETL_TYPE_MATRIX, m);
            }
        }
        /* note: for use by genr_get_output_matrix() */
        p->lh.mret = m;
    } else if (p->targ == LIST) {
        create_or_edit_list(p);
    } else if (p->targ == STR) {
        create_or_edit_string(p);
    } else if (p->targ == BUNDLE) {
        if (null_node(r)) {
            /* as in "bundle b = null" */
            p->err = assign_null_to_bundle(p);
        } else if (p->op != B_ASN) {
            do_bundle_append(p);
        } else {
            /* full assignment of RHS bundle */
            gretl_bundle *b;

            if (r->t == DBUNDLE) {
                b = bvar_get_bundle(r, p);
            } else if (is_tmp_node(r) || (p->flags & P_UFRET)) {
                /* grabbing r->v.b is OK */
                b = r->v.b;
            } else {
                /* we need to make a copy */
                b = gretl_bundle_copy(r->v.b, &p->err);
            }

            if (!p->err) {
                p->err = gen_add_or_replace(p, GRETL_TYPE_BUNDLE, b);
                if (!p->err && r->t != DBUNDLE && b == r->v.b) {
                    /* avoid destroying the assigned bundle */
                    r->v.b = NULL;
                }
            }
        }
    } else if (p->targ == ARRAY) {
        if (p->op == B_ADD) {
            do_array_append(p);
        } else if (p->op == B_SUB) {
            do_array_subtract(p);
        } else if (null_node(r)) {
            /* as in, e.g., "strings A = null" */
            p->err = assign_null_to_array(p);
        } else if (r->t == ARRAY) {
            /* full assignment of RHS array */
            GretlType atype = gretl_array_get_type(r->v.a);
            gretl_array *a = NULL;

            if (p->lh.gtype > 0 && atype != p->lh.gtype) {
                if (atype == GRETL_TYPE_ANY) {
                    gretl_array_set_type(r->v.a, p->lh.gtype);
                } else {
                    p->err = E_TYPES;
                }
            } else if (is_tmp_node(r) || (p->flags & P_UFRET)) {
                /* grabbing r->v.a is OK */
                a = r->v.a;
            } else {
                /* we need to make a copy */
                a = gretl_array_copy(r->v.a, &p->err);
            }
            if (!p->err) {
                p->err = gen_add_or_replace(p, atype, a);
                if (!p->err && a == r->v.a) {
                    /* avoid destroying the assigned array */
                    r->v.a = NULL;
                }
            }
        } else if (r->t == SERIES) {
            /* verified in advance as string-valued */
            GretlType atype = GRETL_TYPE_STRINGS;

            if (p->lh.gtype > 0 && p->lh.gtype != atype) {
                p->err = E_TYPES;
            } else {
                gretl_array *a =
                    get_strings_array_from_series(p->dset, r->vnum, &p->err);

                if (!p->err) {
                    p->err = gen_add_or_replace(p, atype, a);
                }
            }
        } else {
            /* Allow promotion of a single object to an array of
               size 1? Note 2021-08-12: not sure this is actually
               a good idea.
            */
            GretlType rtype = gretl_type_from_gen_type(r->t);
            GretlType atype = p->lh.gtype;

            if (rtype == gretl_type_get_singular(atype)) {
                gretl_array *a = gretl_singleton_array(r->v.ptr, atype,
                                                       1, &p->err);

                if (!p->err) {
                    p->err = gen_add_or_replace(p, atype, a);
                }
            } else {
                p->err = e_types(r);
            }
        }
    }

    return savegen_retval(p->err);
}

static void maybe_update_lhs_uvar (parser *p, GretlType *type)
{
    if (p->targ == SERIES) {
        /* targetting a series */
        int v = p->lh.vnum;

        if (get_loop_renaming() || v <= 0 || v >= p->dset->v) {
            p->lh.vnum = current_series_index(p->dset, p->lh.name);
        }
        if (p->lh.vnum < 0) {
            p->lh.vnum = 0;
        }
        return;
    }

    if (p->lh.uv == NULL) {
        p->lh.uv = get_user_var_by_name(p->lh.name);
    }

    if (p->lh.uv != NULL) {
        *type = p->lh.uv->type;
    }

    switch (*type) {
    case GRETL_TYPE_DOUBLE:
        p->lh.t = NUM;
        break;
    case GRETL_TYPE_MATRIX:
        p->lh.t = MAT;
        if (p->targ == NUM) {
            p->targ = MAT;
        }
        break;
    case GRETL_TYPE_LIST:
        p->lh.t = LIST;
        break;
    case GRETL_TYPE_STRING:
        p->lh.t = STR;
        break;
    case GRETL_TYPE_BUNDLE:
        p->lh.t = BUNDLE;
        break;
    case GRETL_TYPE_ARRAY:
        p->lh.t = ARRAY;
        break;
    default:
        p->lh.t = 0;
        break;
    }
}

static void parser_reinit (parser *p, DATASET *dset, PRN *prn)
{
    /* flags that should be reinstated if they were
       set at compile time, or in previous execution
    */
    int saveflags[] = {
        P_NATEST, P_AUTOREG, P_DISCARD, P_NODECL,
        P_LISTDEF, P_PRIV, P_UFRET, 0
    };
    int i, prevflags = p->flags;
    GretlType lhtype = 0;
    int dset_n;

    if (p->callcount > 1) {
        /* the flags should basically have stabilized by now */
        p->flags |= P_START;
        p->flags &= ~P_DELTAN;
    } else {
        p->flags = (P_START | P_EXEC);
        for (i=0; saveflags[i] > 0; i++) {
            if (prevflags & saveflags[i]) {
                p->flags |= saveflags[i];
            }
        }
    }

    p->dset = dset;
    p->prn = prn;

    p->obs = 0;
    p->sym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->data = NULL;
    p->errprn = NULL;

    p->ret = NULL;
    p->lhres = NULL;
    p->err = 0;

#if EDEBUG
    fprintf(stderr, "parser_reinit: targ=%s, lhname='%s', op='%s', "
            "callcount=%d, compiled=%d, priv=%d\n",
            getsymb(p->targ), p->lh.name, getsymb(p->op),
            p->callcount, compiled(p) ? 1 : 0,
            (p->flags & P_PRIV)? 1 : 0);
#endif

    if (*p->lh.name != '\0') {
        maybe_update_lhs_uvar(p, &lhtype);
    }

    /* allow for change in length of dataset */
    dset_n = dset != NULL ? dset->n : 0;
    if (dset_n != p->dset_n) {
        p->dset_n = dset_n;
        p->flags |= P_DELTAN;
    }
}

static void parser_init (parser *p, const char *str,
                         DATASET *dset, PRN *prn,
                         int flags, int targtype,
                         int *done)
{
    p->point = p->rhs = p->input = str;
    p->dset = dset;
    p->dset_n = dset != NULL ? dset->n : 0;
    p->prn = prn;
    p->errprn = NULL;
    p->flags = flags | P_START;
    p->targ = targtype;
    p->op = 0;

    p->lhtree = NULL;
    p->lhres = NULL;
    p->tree = NULL;
    p->ret = NULL;

    /* left-hand side info */
    p->lh.t = 0;
    p->lh.name[0] = '\0';
    p->lh.label = NULL;
    p->lh.stab = NULL;
    p->lh.vnum = 0;
    p->lh.uv = NULL;
    p->lh.expr = NULL;
    p->lh.gtype = 0;
    p->lh.mret = NULL;

    /* auxiliary apparatus */
    p->aux = NULL;

    p->callcount = 0;
    p->obs = 0;
    p->sym = 0;
    p->upsym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->err = 0;

    if (p->input == NULL) {
        p->err = E_DATA;
        return;
    }

    if (p->flags & P_VOID) {
        p->flags |= P_DISCARD;
    } else if (p->targ == UNK || !(p->flags & P_ANON)) {
        gen_preprocess(p, flags, done);
    } else if (p->targ == LIST) {
        p->flags |= P_LISTDEF;
    }

    if (!p->err) {
        p->ch = parser_getc(p);
    }
}

/* called from genmain.c (only!) */

void gen_save_or_print (parser *p, PRN *prn)
{
    if (autoreg(p)) {
        /* no transcription required */
        return;
    }
    if (p->flags & P_DISCARD) {
        /* doing "eval" */
        if (p->ret == NULL) {
            return;
        } else if (p->ret->t == EMPTY) {
            pputs(prn, "[null]\n");
        } else if (p->ret->t == MAT) {
            if (p->ret->v.m->is_complex) {
                gretl_cmatrix_print(p->ret->v.m, p->lh.name, p->prn);
            } else {
                gretl_matrix_print_to_prn(p->ret->v.m, p->lh.name, p->prn);
            }
        } else if (p->ret->t == LIST) {
            if (p->lh.name[0] != '\0') {
                gretl_list_print(get_list_by_name(p->lh.name),
                                 p->dset, p->prn);
            } else {
                gretl_list_print(p->ret->v.ivec, p->dset, p->prn);
            }
        } else if (p->ret->t == STR) {
            if (p->lh.name[0] != '\0') {
                pprintf(p->prn, "%s\n", gen_get_lhs_var(p, GRETL_TYPE_STRING));
            } else {
                pprintf(p->prn, "%s\n", p->ret->v.str);
            }
        } else if (p->ret->t == BUNDLE) {
            gretl_bundle_print(p->ret->v.b, prn);
        } else if (p->ret->t == ARRAY) {
            gretl_array_print(p->ret->v.a, prn);
        } else {
            /* scalar, series */
            printnode(p->ret, p, 1);
            pputc(p->prn, '\n');
        }
    } else if (p->flags & P_DECL) {
        do_declaration(p);
    } else {
        save_generated_var(p, prn);
    }
}

void gen_cleanup (parser *p)
{
    int save = (p->flags & (P_COMPILE | P_EXEC));

#if EDEBUG
    fprintf(stderr, "gen cleanup on %p: save = %d\n",
            p, save ? 1 : 0);
#endif

    if (p->lh.label != NULL) {
        free(p->lh.label);
        p->lh.label = NULL;
    }

    if (p->lh.stab != NULL) {
        series_table_destroy(p->lh.stab);
        p->lh.stab = NULL;
    }

    if (p->flags & P_ALTINP) {
        free((char *) p->input);
        p->input = NULL;
    }

    if (p->err && (p->flags & P_COMPILE)) {
        save = 0;
    }

    if (!save) {
        if (p->lhtree != NULL) {
            if (p->lhtree != p->lhres) {
                /* we have to scrub the LHT_NODE flag on p->lhtree,
                   or else its children will not get freed and we'll
                   leak memory
                */
                p->lhtree->flags &= ~LHT_NODE;
                rndebug(("freeing p->lhtree %p\n", (void *) p->lhtree));
                free_tree(p->lhtree, p, FR_LHTREE);
            }
            if (p->lhres != NULL) {
                rndebug(("freeing p->lhres %p\n", (void *) p->lhres));
                free_tree(p->lhres, p, FR_LHRES);
            }
        }

        if (p->tree != p->ret) {
            rndebug(("freeing p->tree %p\n", (void *) p->tree));
            free_tree(p->tree, p, FR_TREE);
        }

        rndebug(("freeing p->ret %p\n", (void *) p->ret));
        free_tree(p->ret, p, FR_RET);

        free(p->lh.expr);
    }

#if EDEBUG
    fprintf(stderr, "gen cleanup finished\n");
#endif
}

#define LS_DEBUG 0

static void real_reset_uvars (parser *p)
{
    if (0 && p->err) {
        /* 2022-09-13: not sure about this */
        return;
    }

#if LS_DEBUG
    fprintf(stderr, "\nreal_reset_uvars (%s '%s') *\n",
            getsymb(p->targ), p->lh.name);
#endif

    clear_uvnodes(p->tree);

    if (p->lhtree != NULL) {
        clear_uvnodes(p->lhtree);
    }

    p->lh.uv = NULL;
    p->lh.vnum = 0;
}

void genr_reset_uvars (parser *p)
{
    real_reset_uvars(p);
}

static void maybe_set_return_flags (parser *p)
{
    NODE *t = p->tree;

    if (t != NULL && t->t == UFUN) {
        p->flags |= P_UFRET;
    }
}

static int decl_check (parser *p, int flags)
{
    if (flags & P_COMPILE) {
        p->err = E_EQN;
        gretl_errmsg_sprintf("%s:\n> '%s'",
                             _("Bare declarations are not allowed here"),
                             p->input);
    }

    return p->err;
}

static void autoreg_error (parser *p, int t)
{
    fprintf(stderr, "*** autoreg error at obs t = %d (t1 = %d):\n",
            t, p->dset->t1);

    if (p->ret != NULL && p->ret->t != SERIES) {
        fprintf(stderr, " ret type != SERIES (=%d), p->err = %d\n",
                p->ret->t, p->err);
    } else if (p->ret == NULL) {
        fprintf(stderr, " ret = NULL, p->err = %d\n", p->err);
    }

    fprintf(stderr, " input = '%s'\n", p->input);

    if (!p->err) {
        p->err = E_DATA;
    }
}

#if EDEBUG

void autoreg_genr_report (const double *x, const double *y,
                          int initted, parser *p)
{
    int t = p->obs;

    fprintf(stderr, "\n*** autoreg: p->obs = %d\n", t);
    if (!initted && na(x[t])) {
        fprintf(stderr, "skipping xvec[%d], leaving y[%d] = %g\n",
                t, t, y[t]);
    } else if (p->op == B_ASN) {
        fprintf(stderr, "writing xvec[%d] = %g into y[%d] (was %g)\n",
                t, x[t], t, y[t]);
    } else {
        fprintf(stderr, "using xvec[%d] = %g to modify y[%d] (was %g)\n",
                t, x[t], t, y[t]);
    }
}

#endif

int realgen (const char *s, parser *p, DATASET *dset, PRN *prn,
             int flags, int targtype)
{
#if LHDEBUG || EDEBUG || AUX_NODES_DEBUG
    fprintf(stderr, "\n*** realgen: task = %s, depth %d, priv %d\n",
            (flags & P_COMPILE)? "compile" : (flags & P_EXEC)? "exec" : "normal",
            gretl_function_depth(), (flags & P_PRIV)? 1 : 0);
    if (s != NULL) {
        fprintf(stderr, "targ=%d (%s), input='%s'\n", targtype,
                (targtype < PUNCT_MAX)? gretl_type_get_name(targtype) :
                getsymb(targtype), s);
    }
#endif

    if (flags & P_EXEC) {
#if EDEBUG
        fprintf(stderr, "*** printing p->tree (before reinit)\n");
        print_tree(p->tree, p, 0, 0);
#endif
        parser_reinit(p, dset, prn);
        if (p->err) {
            fprintf(stderr, "error in parser_reinit\n");
            goto gen_finish;
        } else if (p->op == INC || p->op == DEC) {
            /* more or less a no-op: the work is done by
               save_generated_var()
            */
            goto gen_finish;
        } else {
            goto starteval;
        }
    } else {
        int done = 0;

        parser_init(p, s, dset, prn, flags, targtype, &done);
        if (p->err) {
            if (gretl_function_depth() == 0) {
                errmsg(p->err, prn);
            }
            goto gen_finish;
        } else if (done) {
            goto gen_finish;
        }
    }

#if EDEBUG
    fprintf(stderr, "after parser %s, p->err = %d (decl? %s)\n",
            (flags & P_EXEC)? "reinit" : "init", p->err,
            (p->flags & P_DECL)? "yes" : "no");
#endif

    if (p->flags & P_DECL) {
        /* check validity of declaration(s) */
        decl_check(p, flags);
        goto gen_finish;
    }

    if (p->op == INC || p->op == DEC) {
        /* implemented via save_generated_var() */
        goto gen_finish;
    }

    /* fire up the lexer */
    lex(p);
    if (p->err) {
#if EDEBUG
        fprintf(stderr, "realgen %p ('%s'): got on lex() error %d\n",
                (void *) p, s, p->err);
#endif
        goto gen_finish;
    }

    /* build the syntax tree */
    p->tree = expr(p);
    if (p->err) {
        goto gen_finish;
    }

#if EDEBUG
    if (p->tree != NULL) {
        fprintf(stderr, "realgen: p->tree at %p, type %d (%s)\n", (void *) p->tree,
                p->tree->t, getsymb(p->tree->t));
    }
    if (p->ch == '\0') {
        fprintf(stderr, " p->ch = NUL, p->sym = %d\n", p->sym);
    } else {
        fprintf(stderr, " p->ch = '%c', p->sym = %d\n", p->ch, p->sym);
    }
#endif

    if (p->sym != EOT || p->ch != 0) {
        int c = p->ch;

        if (c == ' ') {
            c = 0;
        } else if (c != 0) {
            parser_ungetc(p);
            c = p->ch;
        }
        context_error(c, p, "realgen");
        goto gen_finish;
    }

    if (flags & P_NOEXEC) {
        /* we're done at this point */
        goto gen_finish;
    }

    if (!p->err) {
        /* set P_UFRET here if relevant */
        maybe_set_return_flags(p);
    }

 starteval:

#if EDEBUG
    if (flags & P_EXEC) {
        fprintf(stderr, "*** printing p->tree (about to start eval)\n");
        print_tree(p->tree, p, 0, 0);
    }
#endif

    if (autoreg(p)) {
        /* e.g. y = b*y(-1) : evaluate dynamically */
        double *y = p->dset->Z[p->lh.vnum];
        const double *x;
        int t, initted = 0;

        for (t=p->dset->t1; t<=p->dset->t2 && !p->err; t++) {
            /* initialize for this observation */
            p->obs = t;
            if (dataset_is_panel(p->dset) && t % p->dset->pd == 0) {
                initted = 0;
            }
            p->ret = eval(p->tree, p);
            if (p->ret != NULL && p->ret->t == SERIES) {
                x = p->ret->v.xvec;
#if EDEBUG
                autoreg_genr_report(x, y, initted, p);
#endif
                if (!initted && na(x[t])) {
                    ; /* don't overwrite initializer */
                } else {
                    if (p->op == B_ASN) {
                        y[t] = x[t];
                    } else {
                        y[t] = xy_calc(y[t], x[t], p->op, SERIES, p);
                    }
                    initted = 1;
                }
            } else {
                autoreg_error(p, t);
            }
            if (t == p->dset->t1) {
                p->flags &= ~P_START;
            }
        }
    } else {
        /* standard non-dynamic evaluation */
        p->ret = eval(p->tree, p);
    }

    if (p->flags & P_EXEC) {
        p->callcount += 1;
    }

 gen_finish:

    if (p->errprn != NULL) {
        /* Pick and forward any error message that may not be
           seen if realgen was invoked with a NULL value for
           the printer @prn.
        */
        const char *buf = gretl_print_get_buffer(p->errprn);

        if (buf != NULL && *buf != '\0') {
            gretl_errmsg_set(buf);
        }
        gretl_print_destroy(p->errprn);
        p->errprn = NULL;
    }

#if EDEBUG
    fprintf(stderr, "realgen: at finish, err = %d\n", p->err);
# if EDEBUG > 1
    printnode(p->ret, p, 0);
    pputc(prn, '\n');
# endif
#endif

    return p->err;
}
