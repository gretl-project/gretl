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

/* monte_carlo.c - loop procedures */

#include "libgretl.h"
#include "monte_carlo.h"
#include "libset.h"
#include "cmd_private.h"
#include "var.h"
#include "objstack.h"
#include "gretl_func.h"
#include "uservar.h"
#include "uservar_priv.h"
#include "flow_control.h"
#include "system.h"
#include "gen_public.h"
#include "gretl_string_table.h"
#include "gretl_typemap.h"
#include "genr_optim.h"

#include <time.h>
#include <unistd.h>

#define LOOP_DEBUG 0
#define SUBST_DEBUG 0
#define COMP_DEBUG 0

enum loop_types {
    COUNT_LOOP,
    WHILE_LOOP,
    INDEX_LOOP,
    DATED_LOOP,
    FOR_LOOP,
    EACH_LOOP
};

#define DEFAULT_NOBS 512

#define indexed_loop(l) (l->type == INDEX_LOOP || \
                         l->type == DATED_LOOP || \
                         l->type == EACH_LOOP)

#if HAVE_GMP
# include <gmp.h>
typedef mpf_t bigval;
# include "prog_loop.h"
#endif

typedef enum {
    LOOP_PROGRESSIVE = 1 << 0,
    LOOP_VERBOSE     = 1 << 1,
    LOOP_DELVAR      = 1 << 2,
    LOOP_ATTACHED    = 1 << 3,
    LOOP_RENAMING    = 1 << 4,
    LOOP_ERR_CAUGHT  = 1 << 5,
    LOOP_CONDITIONAL = 1 << 6,
    LOOP_DECREMENT   = 1 << 7
} LoopFlags;

struct controller_ {
    double val;            /* evaluated value */
    char vname[VNAMELEN];  /* name of (scalar) variable, if used */
    user_var *uv;          /* pointer to scalar variable */
    int vsign;             /* 1 or -1, if vname is used */
    char *expr;            /* expression to pass to genr, if used */
    GENERATOR *genr;       /* compiled generator */
    int subst;             /* expression uses string substitution? */
};

typedef struct controller_ controller;

typedef enum {
    LOOP_LINE_LIT     = 1 << 0, /* literal printing */
    LOOP_LINE_DOLLAR  = 1 << 1, /* line contains '$' */
    LOOP_LINE_AT      = 1 << 2, /* line contains '@' */
    LOOP_LINE_CATCH   = 1 << 3, /* "catch" flag present */
    LOOP_LINE_PDONE   = 1 << 4, /* progressive loop command started */
    LOOP_LINE_NOCOMP  = 1 << 5, /* not compilable */
    LOOP_LINE_QUIET   = 1 << 6, /* non-printing model in progressive loop */
} LoopLineFlags;

/* local convenience typedef */
typedef struct stmt_ loop_line;

/* information on model estimated within loop */
struct model_record_ {
    int lno;
    MODEL *pmod;
    gretlopt opt;
};

typedef struct model_record_ model_record;

struct LOOPSET_ {
    /* basic characteristics */
    char type;
    LoopFlags flags;
    int level;
    int err;

    /* iterations */
    int itermax;
    int iter;
    int index;

    /* index/foreach control variables */
    char idxname[VNAMELEN];
    user_var *idxvar;
    int idxval;
    char eachname[VNAMELEN];
    GretlType eachtype;

    /* break, continue signals */
    char brk;
    char cont;

    /* control structures */
    controller init;
    controller test;
    controller delta;
    controller final;

    /* numbers of various subsidiary objects */
    int n_lines;
    int n_models;
    int n_children;

    /* subsidiary objects */
    loop_line *lines;      /* saved command info */
    char **eachstrs;       /* for use with "foreach" loop */
    model_record *mrecords;
    LOOPSET *parent;
    LOOPSET **children;

#if HAVE_GMP
    /* "progressive" objects and counts thereof */
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    LOOP_STORE store;
    int n_loop_models;
    int n_prints;
#endif
};

#define loop_is_progressive(l)  (l->flags & LOOP_PROGRESSIVE)
#define loop_set_progressive(l) (l->flags |= LOOP_PROGRESSIVE)
#define loop_is_verbose(l)      (l->flags & LOOP_VERBOSE)
#define loop_set_verbose(l)     (l->flags |= LOOP_VERBOSE)
#define loop_is_attached(l)     (l->flags & LOOP_ATTACHED)
#define loop_set_attached(l)    (l->flags |= LOOP_ATTACHED)
#define loop_is_renaming(l)     (l->flags & LOOP_RENAMING)
#define loop_set_renaming(l)    (l->flags |= LOOP_RENAMING)
#define loop_has_cond(l)        (l->flags & LOOP_CONDITIONAL)
#define loop_set_has_cond(l)    (l->flags |= LOOP_CONDITIONAL)
#define loop_decrement(l)       (l->flags & LOOP_DECREMENT)
#define loop_set_decrement(l)   (l->flags |= LOOP_DECREMENT)

#define loop_line_catch(ll)  (ll->flags & LOOP_LINE_CATCH)
#define loop_line_nosub(ll)  (!(ll->flags & (LOOP_LINE_AT | LOOP_LINE_DOLLAR)))
#define loop_line_quiet(ll)  (ll->flags & LOOP_LINE_QUIET)

#define line_is_compiled(ll) (ll->ptr != NULL || ll->ci == ELSE || \
                              ll->ci == ENDIF || ll->ci == BREAK || \
                              ll->ci == CONTINUE)
#define line_has_dollar(ll)  (ll->flags & LOOP_LINE_DOLLAR)
#define line_has_at_sign(ll) (ll->flags & LOOP_LINE_AT)
#define line_has_genr(l,j) (l->lines[j].ptr != NULL && l->lines[j].ci != LOOP)

#define may_be_compilable(ll)  (!(ll->flags & LOOP_LINE_NOCOMP))
#define set_non_compilable(ll) (ll->flags |= LOOP_LINE_NOCOMP)

#define model_print_deferred(o) (o & OPT_F)

static void controller_init (controller *clr);
static int gretl_loop_prepare (LOOPSET *loop);
static void controller_free (controller *clr);
static void destroy_loop_stack (LOOPSET *loop);

static int
make_dollar_substitutions (char *str, int maxlen,
                           const LOOPSET *loop,
                           const DATASET *dset,
                           int *subst,
                           gretlopt opt);

#define LOOP_BLOCK 32

/* record of state, and communication of state with outside world */

static LOOPSET *currloop;

static int compile_level;
static int loop_execute;
static int loop_renaming;

int gretl_compiling_loop (void)
{
    return compile_level;
}

void gretl_abort_compiling_loop (void)
{
    if (currloop != NULL) {
        destroy_loop_stack(currloop);
    }
}

int gretl_execute_loop (void)
{
    return loop_execute;
}

int get_loop_renaming (void)
{
    return loop_renaming;
}

/* Test for a "while" or "for" expression in a loop controller:
   if it involves string substitution we can't compile it.
*/

static int controller_does_string_sub (const char *s,
				       LOOPSET *loop,
				       DATASET *dset)
{
    int subst = 0;

    if (strchr(s, '@')) {
        subst = 1;
    } else if (strchr(s, '$')) {
        char test[64];

        *test = '\0';
        strncat(test, s, 63);
        make_dollar_substitutions(test, 63, loop, dset,
                                  &subst, OPT_T);
    }

    return subst;
}

static double get_loop_index_value (double x)
{
    double xi = nearbyint(x);

    if (fabs(x - xi) < 1.0e-15) {
	return xi;
    } else {
	return x;
    }
}

/* For indexed loops: get a value from a loop "limit" element (lower
   or upper).  If we got the name of a scalar variable at setup time,
   look up its current value (and modify the sign if wanted).  Or if
   we got a "genr" expression, evaluate it.  Otherwise we should have
   got a numerical constant at setup, in which case we just return
   that value.
*/

static double controller_get_val (controller *clr,
                                  LOOPSET *loop,
                                  DATASET *dset,
                                  int *err)
{
    /* note: check for "compiled" variants first */
    if (clr->uv != NULL) {
        clr->val = uvar_get_scalar_value(clr->uv) * clr->vsign;
    } else if (clr->genr != NULL) {
        clr->val = evaluate_scalar_genr(clr->genr, dset, NULL, err);
    } else if (clr->vname[0] != '\0') {
        if (clr->vname[0] == '$') {
            /* built-in scalar constant */
            clr->val = get_const_by_name(clr->vname, err) * clr->vsign;
        } else {
            /* should be scalar uservar */
            if (clr->uv == NULL) {
                clr->uv = get_user_var_of_type_by_name(clr->vname, GRETL_TYPE_DOUBLE);
            }
            if (clr->uv == NULL) {
                gretl_errmsg_sprintf(_("'%s': not a scalar"), clr->vname);
                *err = E_TYPES;
            } else {
		double d = uvar_get_scalar_value(clr->uv) * clr->vsign;

                clr->val = get_loop_index_value(d);
            }
        }
    } else if (clr->expr != NULL && clr->subst) {
        int done = 0;

        if (strchr(clr->expr, '@')) {
            /* the expression needs string substitution? */
            char expr[64];

            *expr = '\0';
            strncat(expr, clr->expr, 63);
            *err = substitute_named_strings(expr, &clr->subst);
            if (!*err && clr->subst) {
                clr->val = generate_scalar(expr, dset, err);
                done = 1;
            }
        }
        if (!done && !*err && strchr(clr->expr, '$')) {
            /* the expression needs dollar substitution? */
            char expr[64];

            *expr = '\0';
            strncat(expr, clr->expr, 63);
            *err = make_dollar_substitutions(expr, 63, loop, dset,
                                             &clr->subst, OPT_T);
            if (!*err && clr->subst) {
                clr->val = generate_scalar(expr, dset, err);
                done = 1;
            }
        }
        if (!*err && !done) {
            clr->subst = 0;
            clr->val = generate_scalar(clr->expr, dset, err);
        }
    } else if (clr->expr != NULL) {
        /* expression with no string substitution */
        if (clr->genr == NULL) {
            clr->genr = genr_compile(clr->expr, dset, GRETL_TYPE_DOUBLE,
                                     OPT_P | OPT_N | OPT_A, NULL, err);
        }
        if (clr->genr != NULL) {
            clr->val = evaluate_scalar_genr(clr->genr, dset, NULL, err);
        } else {
            /* fallback: or should we just flag an error? */
            *err = 0;
            clr->val = generate_scalar(clr->expr, dset, err);
        }
    }

    if (*err && clr->expr != NULL) {
        gchar *msg;

        msg = g_strdup_printf("Bad loop-control expression '%s'", clr->expr);
        gretl_errmsg_append(msg, *err);
        g_free(msg);
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "clr->val: loop %p, vname='%s', uv=%p, expr='%s' val=%g, err=%d\n",
            (void *) loop, clr->vname, (void *) clr->uv, clr->expr, clr->val, *err);
#endif

    return clr->val;
}

/* apply initialization in case of for-loop */

static void
forloop_init (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->init.expr;

    if (expr != NULL) {
        *err = generate(expr, dset, GRETL_TYPE_ANY, OPT_Q, NULL);
        if (*err) {
            gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
                                 expr);
        }
    }
}

/* evaluate boolean condition in for-loop or while-loop */

static int
loop_testval (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->test.expr;
    int ret = 1;

    if (expr != NULL) {
        double x = NADBL;

        if (loop->test.subst < 0) {
            /* not checked yet */
            loop->test.subst = controller_does_string_sub(expr, loop, dset);
        }

        if (!loop->test.subst && loop->test.genr == NULL) {
            loop->test.genr = genr_compile(expr, dset,
                                           GRETL_TYPE_BOOL,
                                           OPT_P | OPT_N,
                                           NULL, err);
        }

        if (loop->test.genr != NULL) {
            x = evaluate_if_cond(loop->test.genr, dset, NULL, err);
        } else if (!*err) {
            x = generate_scalar(expr, dset, err);
        }

        if (!*err && na(x)) {
            *err = E_DATA;
            ret = 0;
        } else {
            ret = x;
        }
        if (*err) {
            gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
                                 expr);
        }
    }

    return ret;
}

/* evaluate third expression in for-loop, if any */

static void
loop_delta (LOOPSET *loop, DATASET *dset, int *err)
{
    const char *expr = loop->delta.expr;

    if (expr != NULL) {
        if (loop->delta.subst < 0) {
            /* not checked yet */
            loop->delta.subst = controller_does_string_sub(expr, loop, dset);
        }

        if (!loop->delta.subst && loop->delta.genr == NULL) {
            loop->delta.genr = genr_compile(expr, dset,
                                            GRETL_TYPE_ANY,
                                            OPT_N,
                                            NULL, err);
        }

        if (loop->delta.genr != NULL) {
            *err = execute_genr(loop->delta.genr, dset, NULL);
        } else if (!*err) {
            *err = generate(expr, dset, GRETL_TYPE_ANY, OPT_Q, NULL);
        }
        if (*err) {
            gretl_errmsg_sprintf("%s: '%s'", _("error evaluating loop condition"),
                                 expr);
        }
    }
}

static void set_loop_opts (LOOPSET *loop, gretlopt opt)
{
    if (opt & OPT_P) {
        loop_set_progressive(loop);
    }
    if (opt & OPT_V) {
        loop_set_verbose(loop);
    }
    if (opt & OPT_D) {
        loop_set_decrement(loop);
    }
}

#define plain_model_ci(c) (MODEL_COMMAND(c) && \
                           c != NLS && \
                           c != MLE && \
                           c != GMM)

/**
 * ok_in_loop:
 * @ci: command index.
 *
 * Returns: 1 if the given command is acceptable inside the loop construct,
 * 0 otherwise.
 */

int ok_in_loop (int c)
{
    /* here are the commands we _don't_ currently allow */

    if (c == FUNC ||
        c == INCLUDE ||
        c == NULLDATA ||
        c == RUN ||
        c == SETMISS ||
        c == QUIT) {
        return 0;
    }

#if 0
    if (c == RENAME) return 0;
#endif

    return 1;
}

static int loop_attach_child (LOOPSET *loop, LOOPSET *child)
{
    LOOPSET **children;
    int nc = loop->n_children;

    children = realloc(loop->children, (nc + 1) * sizeof *children);
    if (children == NULL) {
        return E_ALLOC;
    }

    loop->children = children;
    loop->children[nc] = child;
    child->parent = loop;
    child->level = loop->level + 1;

#if LOOP_DEBUG
    fprintf(stderr, "child loop %p has parent %p\n",
            (void *) child, (void *) child->parent);
#endif

    loop->n_children += 1;

    return 0;
}

static void gretl_loop_init (LOOPSET *loop)
{
#if LOOP_DEBUG > 1
    fprintf(stderr, "gretl_loop_init: initing loop at %p\n", (void *) loop);
#endif

    loop->flags = 0;
    loop->level = 0;

    loop->itermax = 0;
    loop->iter = 0;
    loop->err = 0;
    *loop->idxname = '\0';
    loop->idxvar = NULL;
    loop->idxval = 0;
    loop->brk = loop->cont = 0;
    *loop->eachname = '\0';
    loop->eachtype = 0;
    loop->eachstrs = NULL;

    controller_init(&loop->init);
    controller_init(&loop->test);
    controller_init(&loop->delta);
    controller_init(&loop->final);

    loop->n_lines = 0;
    loop->lines = NULL;
    loop->n_models = 0;
    loop->mrecords = NULL;

    loop->parent = NULL;
    loop->children = NULL;
    loop->n_children = 0;

#if HAVE_GMP
    /* "progressive" apparatus */
    loop->n_loop_models = 0;
    loop->lmodels = NULL;
    loop->n_prints = 0;
    loop->prns = NULL;
    loop_store_init(&loop->store);
#endif
}

static LOOPSET *gretl_loop_new (LOOPSET *parent)
{
    LOOPSET *loop = malloc(sizeof *loop);

    if (loop == NULL) {
        return NULL;
    }

    gretl_loop_init(loop);

    if (parent != NULL) {
        int err = loop_attach_child(parent, loop);

        if (err) {
            free(loop);
            loop = NULL;
        }
    }

    return loop;
}

void gretl_loop_destroy (LOOPSET *loop)
{
    int i;

    if (loop == NULL) {
        return;
    }

#if COMP_DEBUG || LOOP_DEBUG
    fprintf(stderr, "destroying LOOPSET at %p\n", (void *) loop);
#endif

    if (loop->children != NULL) {
        for (i=0; i<loop->n_children; i++) {
            gretl_loop_destroy(loop->children[i]);
            loop->children[i] = NULL;
        }
        free(loop->children);
        loop->children = NULL;
    }

    controller_free(&loop->init);
    controller_free(&loop->test);
    controller_free(&loop->delta);
    controller_free(&loop->final);

    if (loop->lines != NULL) {
        for (i=0; i<loop->n_lines; i++) {
            free(loop->lines[i].s);
            if (line_has_genr(loop, i)) {
                destroy_genr(loop->lines[i].ptr);
            }
        }
        free(loop->lines);
    }

    if (loop->mrecords != NULL) {
        free(loop->mrecords);
    }

    if (loop->eachstrs != NULL && loop->eachtype != GRETL_TYPE_STRINGS) {
        strings_array_free(loop->eachstrs, loop->itermax);
    }

#if HAVE_GMP
    if (loop->lmodels != NULL) {
        for (i=0; i<loop->n_loop_models; i++) {
            loop_model_free(&loop->lmodels[i]);
        }
        free(loop->lmodels);
    }
    if (loop->prns != NULL) {
        for (i=0; i<loop->n_prints; i++) {
            loop_print_free(&loop->prns[i]);
        }
        free(loop->prns);
    }
    loop_store_free(&loop->store);
#endif

    if (loop->flags & LOOP_DELVAR) {
        user_var_delete_by_name(loop->idxname, NULL);
    }

    free(loop);
}

static void destroy_loop_stack (LOOPSET *loop)
{
    if (loop == NULL) {
        return;
    }

    /* find the origin of the stack */
    while (loop->parent != NULL) {
        loop = loop->parent;
    }

    /* and destroy recursively */
    gretl_loop_destroy(loop);

    compile_level = 0;
    loop_renaming = 0;
    set_loop_off();
    currloop = NULL;
}

static int parse_as_while_loop (LOOPSET *loop, const char *s)
{
    int err = 0;

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_as_while_loop: cond = '%s'\n", s);
#endif

    if (s == NULL || *s == '\0') {
        err = E_PARSE;
    } else {
        loop->type = WHILE_LOOP;
        loop->test.expr = gretl_strdup(s);
        if (loop->test.expr == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static int check_index_in_parentage (LOOPSET *loop, const char *vname)
{
    int thistype = loop->type;

    while ((loop = loop->parent) != NULL) {
        if ((loop->type != FOR_LOOP || loop->type != thistype) &&
            strcmp(vname, loop->idxname) == 0) {
            gretl_errmsg_sprintf(_("Using the same index variable (%s) for nested loops:\n"
                                   "this is acceptable only with \"for\" loops."), vname);
            return E_DATA;
        }
    }

    return 0;
}

static user_var *get_local_scalar_by_name (const char *s, int *err)
{
    user_var *u = get_user_var_by_name(s);

    if (u == NULL) {
        /* no pre-existing var, OK */
        return NULL;
    } else if (u->type != GRETL_TYPE_DOUBLE) {
        gretl_errmsg_set("loop index must be a scalar");
        *err = E_TYPES;
        return NULL;
    } else {
        return u;
    }
}

/* The following is called only once, at the point of initial
   "compilation" of a loop.
*/

static int loop_attach_index_var (LOOPSET *loop,
                                  const char *vname,
                                  DATASET *dset)
{
    int err = 0;

    if (loop->parent != NULL) {
        err = check_index_in_parentage(loop, vname);
        if (err) {
            return err;
        }
    }

    loop->idxvar = get_local_scalar_by_name(vname, &err);

    if (loop->idxvar != NULL) {
        strcpy(loop->idxname, vname);
        uvar_set_scalar_fast(loop->idxvar, loop->init.val);
    } else if (!err) {
        /* create index var from scratch */
        char genline[64];

        if (na(loop->init.val)) {
            sprintf(genline, "%s=NA", vname);
        } else {
            gretl_push_c_numeric_locale();
            sprintf(genline, "%s=%g", vname, loop->init.val);
            gretl_pop_c_numeric_locale();
        }

        err = generate(genline, dset, GRETL_TYPE_DOUBLE, OPT_Q, NULL);

        if (!err) {
            /* automatic index variable */
            strcpy(loop->idxname, vname);
            loop->idxvar = get_user_var_by_name(vname);
            loop->flags |= LOOP_DELVAR;
        }
    }

    return err;
}

/* for a loop control expression such as "j=start..end", get the
   initial or final value from the string @s (we also use this to get
   the count for a simple count loop).
*/

static int index_get_limit (LOOPSET *loop, controller *clr,
                            const char *s, DATASET *dset)
{
    int v, err = 0;

    if (integer_string(s)) {
        /* plain numerical value */
        clr->val = atoi(s);
    } else {
        if (*s == '-') {
            /* negative of variable? */
            clr->vsign = -1;
            s++;
        }
        if (gretl_is_scalar(s)) {
            *clr->vname = '\0';
            strncat(clr->vname, s, VNAMELEN - 1);
	    clr->val = get_loop_index_value(gretl_scalar_get_value(s, NULL));
        } else if ((v = current_series_index(dset, s)) >= 0) {
            /* found a series by the name of @s */
            gretl_errmsg_sprintf(_("'%s': not a scalar"), s);
        } else if (loop->parent != NULL && strlen(s) == gretl_namechar_spn(s)) {
            /* potentially valid varname, but unknown at present */
            *clr->vname = '\0';
            strncat(clr->vname, s, VNAMELEN - 1);
        } else {
            /* expression to be evaluated to scalar? */
            clr->expr = gretl_strdup(s);
            if (clr->expr == NULL) {
                err = E_ALLOC;
            }
        }
    }

    return err;
}

#define maybe_date(s) (strchr(s, ':') || strchr(s, '/'))

static int parse_as_indexed_loop (LOOPSET *loop,
                                  DATASET *dset,
                                  const char *lvar,
                                  const char *start,
                                  const char *end)
{
    int err = 0;

    /* starting and ending values: the order in which we try
       for valid values is: dates, numeric constants,
       named scalars, scalar expressions.
    */

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_as_indexed_loop: start='%s', end='%s'\n", start, end);
#endif

    if (maybe_date(start)) {
        loop->init.val = dateton(start, dset);
        if (loop->init.val < 0) {
            err = E_DATA;
        } else {
            loop->init.val += 1;
            loop->final.val = dateton(end, dset);
            if (loop->final.val < 0) {
                err = E_DATA;
            } else {
                loop->final.val += 1;
                loop->type = DATED_LOOP;
            }
        }
    } else {
        err = index_get_limit(loop, &loop->init, start, dset);
        if (!err) {
            err = index_get_limit(loop, &loop->final, end, dset);
        }
        if (!err) {
            loop->type = INDEX_LOOP;
        }
    }

    if (!err) {
        err = loop_attach_index_var(loop, lvar, dset);
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "indexed_loop: init.val=%.16f, final.val=%.16f, err=%d\n",
            loop->init.val, loop->final.val, err);
#endif

    return err;
}

/* for example, "loop 100" or "loop K" */

static int parse_as_count_loop (LOOPSET *loop,
                                DATASET *dset,
                                const char *s)
{
    int err;

    err = index_get_limit(loop, &loop->final, s, dset);

    if (!err) {
        loop->init.val = 1;
        loop->type = COUNT_LOOP;
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_as_count_loop: init.val=%g, final.val=%g\n",
            loop->init.val, loop->final.val);
#endif

    return err;
}

static int set_forloop_element (char *s, LOOPSET *loop, int i)
{
    controller *clr = (i == 0)? &loop->init :
        (i == 1)? &loop->test : &loop->delta;
    int len, err = 0;

#if LOOP_DEBUG > 1
    fprintf(stderr, "set_forloop_element: i=%d: '%s'\n", i, s);
#endif

    if (s == NULL || *s == '\0') {
        /* an empty "for" field */
        if (i == 1) {
            /* test is implicitly always true */
            clr->val = 1;
        } else {
            /* no-op */
            clr->val = 0;
        }
        return 0;
    }

    clr->expr = gretl_strdup(s);
    if (clr->expr == NULL) {
        err = E_ALLOC;
    }

    if (!err && i == 0) {
        /* initialization: look for varname for possible substitution */
        err = extract_varname(clr->vname, s, &len);
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, " expr='%s', vname='%s'\n", clr->expr, clr->vname);
#endif

    return err;
}

static int allocate_each_strings (LOOPSET *loop, int n)
{
    loop->eachstrs = strings_array_new(n);

    return (loop->eachstrs == NULL)? E_ALLOC : 0;
}

static int list_vars_to_strings (LOOPSET *loop, const int *list,
                                 const DATASET *dset)
{
    int i, vi;
    int err;

#if LOOP_DEBUG > 1
    fprintf(stderr, "list_vars_to_strings: adding %d strings\n", list[0]);
#endif

    err = allocate_each_strings(loop, list[0]);

    for (i=0; i<list[0] && !err; i++) {
        vi = list[i+1];
        if (vi < 0 || vi >= dset->v) {
            err = E_DATA;
        } else {
            loop->eachstrs[i] = gretl_strdup(dset->varname[vi]);
            if (loop->eachstrs[i] == NULL) {
                err = E_ALLOC;
            }
        }
    }

    return err;
}

static gretl_bundle *look_up_bundle (const char *s)
{
    gretl_bundle *b = get_bundle_by_name(s);

    if (b == NULL) {
	int err = 0;
	b = genr_get_pointer(s, GRETL_TYPE_BUNDLE, &err);
    }

    return b;
}

static gretl_array *look_up_array (const char *s)
{
    gretl_array *a = get_array_by_name(s);

    if (a == NULL) {
	int err = 0;
	a = genr_get_pointer(s, GRETL_TYPE_ARRAY, &err);
    }

    return a;
}

static gretl_array *look_up_strings_array (const char *s)
{
    gretl_array *a = get_strings_array_by_name(s);

    if (a == NULL) {
	a = look_up_array(s);
	if (a != NULL) {
	    if (gretl_array_get_type(a) != GRETL_TYPE_STRINGS) {
		a = NULL;
	    }
	}
    }

    return a;
}

static void *get_eachvar_by_name (const char *s, GretlType *t)
{
    void *ptr = NULL;

    if (*t == GRETL_TYPE_LIST) {
        ptr = get_list_by_name(s);
    } else if (*t == GRETL_TYPE_STRINGS) {
        ptr = look_up_strings_array(s);
    } else if (*t == GRETL_TYPE_BUNDLE) {
        ptr = look_up_bundle(s);
    } else {
        /* type not yet determined */
        if ((ptr = get_list_by_name(s)) != NULL) {
            *t = GRETL_TYPE_LIST;
        } else if ((ptr = look_up_strings_array(s)) != NULL) {
            *t = GRETL_TYPE_STRINGS;
        } else if ((ptr = look_up_bundle(s)) != NULL) {
            *t = GRETL_TYPE_BUNDLE;
        }
    }

    return ptr;
}

/* At loop runtime, check the named list or array and insert the names
   (or numbers) of the variables as "eachstrs"; flag an error if the
   source has disappeared. We also have to handle the case where the
   name of the loop-controlling list is subject to $-substitution.
*/

static int loop_list_refresh (LOOPSET *loop, DATASET *dset)
{
    void *eachvar = NULL;
    const char *strval = NULL;
    int err = 0;

    if (strchr(loop->eachname, '$') != NULL) {
        /* $-string substitution required */
        char vname[VNAMELEN];

        strcpy(vname, loop->eachname);
        err = make_dollar_substitutions(vname, VNAMELEN, loop,
                                        dset, NULL, OPT_T);
        if (!err) {
            eachvar = get_eachvar_by_name(vname, &loop->eachtype);
        }
    } else if (*loop->eachname == '@') {
        /* @-string substitution required */
        strval = get_string_by_name(loop->eachname + 1);
        if (strval != NULL && strlen(strval) < VNAMELEN) {
            eachvar = get_eachvar_by_name(strval, &loop->eachtype);
        }
    } else {
        /* no string substitution needed */
        eachvar = get_eachvar_by_name(loop->eachname, &loop->eachtype);
    }

    /* note: if @eachvar is an array of strings then loop->eachstrs
       will be borrowed data and should not be freed!
    */
    if (loop->eachstrs != NULL) {
        if (loop->eachtype != GRETL_TYPE_STRINGS) {
            strings_array_free(loop->eachstrs, loop->itermax);
        }
        loop->eachstrs = NULL;
    }

    loop->itermax = loop->final.val = 0;

    if (loop->eachtype != GRETL_TYPE_NONE && eachvar == NULL) {
        /* foreach variable has disappeared? */
        err = E_DATA;
    } else if (loop->eachtype == GRETL_TYPE_LIST) {
        int *list = eachvar;

        if (list[0] > 0) {
            err = list_vars_to_strings(loop, list, dset);
            if (!err) {
                loop->final.val = list[0];
            }
        }
    } else if (loop->eachtype == GRETL_TYPE_STRINGS) {
        gretl_array *a = eachvar;
        int n = gretl_array_get_length(a);

        if (n > 0) {
            loop->eachstrs = gretl_array_get_strings(a, &n);
            loop->final.val = n;
        }
    } else if (loop->eachtype == GRETL_TYPE_BUNDLE) {
        gretl_bundle *b = eachvar;
        int n = gretl_bundle_get_n_keys(b);

        if (n > 0) {
            loop->eachstrs = gretl_bundle_get_keys_raw(b, &n);
            loop->final.val = n;
        }
    } else if (!err) {
        /* FIXME do/should we ever come here? */
        if (strval != NULL) {
            /* maybe space separated strings? */
            int nf = 0;

            loop->eachstrs = gretl_string_split_quoted(strval, &nf, NULL, &err);
            if (!err) {
                loop->final.val = nf;
            }
        } else {
            err = E_UNKVAR;
        }
    }

    return err;
}

static GretlType find_target_in_parentage (LOOPSET *loop,
                                           const char *s)
{
    char lfmt[16], afmt[18], vname[VNAMELEN];
    int i;

    sprintf(lfmt, "list %%%d[^ =]", VNAMELEN-1);
    sprintf(afmt, "strings %%%d[^ =]", VNAMELEN-1);

    while ((loop = loop->parent) != NULL) {
        for (i=0; i<loop->n_lines; i++) {
            if (sscanf(loop->lines[i].s, lfmt, vname)) {
                if (!strcmp(vname, s)) {
                    return GRETL_TYPE_LIST;
                }
            } else if (sscanf(loop->lines[i].s, afmt, vname)) {
                if (!strcmp(vname, s)) {
                    return GRETL_TYPE_STRINGS;
                }
            }
        }
    }

    return GRETL_TYPE_NONE;
}

/* We're looking at a "foreach" loop with a single field after the
   index variable, so it's most likely a loop over a list or array.

   We begin by looking for a currently existing named list, but if
   this fails we don't give up immediately.  If we're working on an
   embedded loop, the list may be created within a parent loop whose
   commands have not yet been executed, so we search upward among the
   ancestors of this loop (if any) for a relevant list-creation
   command.

   Even if we find an already-existing list, we do not yet fill out
   the variable-name (or variable-number) strings: these will be set
   when the loop is actually run, since the list may have changed in
   the meantime.

   Besides the possibilities mentioned above, the single field
   may be an @-string that cashes out into one or more "words".
*/

static int list_loop_setup (LOOPSET *loop, char *s, int *nf,
                            int *idxmax)
{
    GretlType t = 0;
    gretl_array *a = NULL;
    gretl_bundle *b = NULL;
    int *list = NULL;
    int len = 0;
    int err = 0;

    while (isspace(*s)) s++;
    tailstrip(s);

    if (*s == '@') {
        /* tricksy: got a list-name that needs string subst? */
        *loop->eachname = '\0';
        strncat(loop->eachname, s, VNAMELEN - 1);
        *nf = 0;
        return 0;
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "list_loop_setup: s = '%s'\n", s);
#endif

    if ((list = get_list_by_name(s)) != NULL) {
        t = GRETL_TYPE_LIST;
        len = list[0];
    } else if ((a = look_up_array(s)) != NULL) {
        t = gretl_array_get_type(a);
        len = gretl_array_get_length(a);
        if (t != GRETL_TYPE_STRINGS) {
            *idxmax = len;
            return 0;
        }
    } else if ((b = look_up_bundle(s)) != NULL) {
        t = GRETL_TYPE_BUNDLE;
        len = gretl_bundle_get_n_keys(b);
    } else {
        t = find_target_in_parentage(loop, s);
    }

    if (t == GRETL_TYPE_NONE) {
	err = E_UNKVAR;
    } else {
        loop->eachtype = t;
        *loop->eachname = '\0';
        strncat(loop->eachname, s, VNAMELEN - 1);
        *nf = len;
    }

    return err;
}

enum {
    DOTTED_LIST,
    WILDCARD_LIST
};

static int
each_strings_from_list_of_vars (LOOPSET *loop, const DATASET *dset,
                                char *s, int *pnf, int type)
{
    int *list = NULL;
    int err = 0;

    if (type == WILDCARD_LIST) {
        s += strspn(s, " \t");
        list = varname_match_list(dset, s, &err);
    } else {
        char vn1[VNAMELEN], vn2[VNAMELEN];
        char fmt[16];

        gretl_delchar(' ', s);
        sprintf(fmt, "%%%d[^.]..%%%ds", VNAMELEN-1, VNAMELEN-1);

        if (sscanf(s, fmt, vn1, vn2) != 2) {
            err = E_PARSE;
        } else {
            int v1 = current_series_index(dset, vn1);
            int v2 = current_series_index(dset, vn2);

            if (v1 < 0 || v2 < 0) {
                err = E_UNKVAR;
            } else if (v2 - v1 + 1 <= 0) {
                err = E_DATA;
            } else {
                list = gretl_consecutive_list_new(v1, v2);
                if (list == NULL) {
                    err = E_ALLOC;
                }
            }
        }
        if (err) {
            *pnf = 0;
        }
    }

    if (list != NULL) {
        int i, vi;

        err = allocate_each_strings(loop, list[0]);
        if (!err) {
            for (i=1; i<=list[0] && !err; i++) {
                vi = list[i];
                loop->eachstrs[i-1] = gretl_strdup(dset->varname[vi]);
                if (loop->eachstrs[i-1] == NULL) {
                    strings_array_free(loop->eachstrs, list[0]);
                    loop->eachstrs = NULL;
                    err = E_ALLOC;
                }
            }
        }
        if (!err) {
            *pnf = list[0];
        }
        free(list);
    }

    return err;
}

/* in context of "foreach" loop, split a string variable by
   both spaces and newlines */

static int count_each_fields (const char *s)
{
    int nf = 0;

    if (s != NULL && *s != '\0') {
        const char *p;

        s += strspn(s, " ");

        if (*s != '\0' && *s != '\n') {
            s++;
            nf++;
        }

        while (*s) {
            p = strpbrk(s, " \n");
            if (p != NULL) {
                s = p + strspn(p, " \n");
                if (*s) {
                    nf++;
                }
            } else {
                break;
            }
        }
    }

    return nf;
}

/* Implement "foreach" for arrays other than strings:
   convert to index loop with automatic max value set
   to the length of the array.
*/

static int set_alt_each_loop (LOOPSET *loop, DATASET *dset,
                              const char *ivar, int len)
{
    loop->type = INDEX_LOOP;
    loop->init.val = 1;
    loop->final.val = len;
    return loop_attach_index_var(loop, ivar, dset);
}

static int
parse_as_each_loop (LOOPSET *loop, DATASET *dset, char *s)
{
    char ivar[VNAMELEN] = {0};
    int done = 0;
    int nf, err = 0;

    /* we're looking at the string that follows "loop foreach" */
    if (*s == '\0') {
        return E_PARSE;
    }

    s += strspn(s, " "); /* skip any spaces */

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_as_each_loop: s = '%s'\n", s);
#endif

    /* get the index variable name (as in "foreach i") */
    if (gretl_scan_varname(s, ivar) != 1) {
        return E_PARSE;
    }

    s += strlen(ivar);
    nf = count_each_fields(s);

#if LOOP_DEBUG > 1
    fprintf(stderr, " number of fields = %d\n", nf);
#endif

    if (nf == 0) {
        return E_PARSE;
    }

    if (nf <= 3 && strstr(s, "..") != NULL) {
        /* range of values, foo..quux */
        err = each_strings_from_list_of_vars(loop, dset, s, &nf,
                                             DOTTED_LIST);
        done = 1;
    } else if (nf == 1 && strchr(s, '*')) {
        err = each_strings_from_list_of_vars(loop, dset, s, &nf,
                                             WILDCARD_LIST);
        done = (err == 0);
    }

    if (!done && nf == 1) {
        /* try for a named list or array? */
        int nelem = -1;

        err = list_loop_setup(loop, s, &nf, &nelem);
        if (!err && nelem >= 0) {
            /* got an array, but not of strings */
            return set_alt_each_loop(loop, dset, ivar, nelem);
        }
        done = (err == 0);
    }

    if (!done) {
        /* simple array of strings: allow for quoted substrings */
        loop->eachstrs = gretl_string_split_quoted(s, &nf, NULL, &err);
    }

    if (!err) {
        loop->type = EACH_LOOP;
        loop->init.val = 1;
        loop->final.val = nf;
        loop->itermax = nf;
        err = loop_attach_index_var(loop, ivar, dset);
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_as_each_loop: final.val=%g\n", loop->final.val);
#endif

    return err;
}

/* try to parse out (expr1; expr2; expr3) */

static int parse_as_for_loop (LOOPSET *loop, char *s)
{
    char *tmp, *q;
    int i, j, len;
    int sc = 0;
    int err = 0;

    s += strcspn(s, "(");
    if (*s != '(') {
        return E_PARSE;
    }

    s++;
    q = strrchr(s, ')');
    if (q == NULL) {
        return E_PARSE;
    }

    len = q - s;
    if (len < 2) { /* minimal OK string is ";;" */
        return E_PARSE;
    }

    tmp = malloc(len + 1);
    if (tmp == NULL) {
        return E_ALLOC;
    }

    for (j=0; j<3 && s!=q && !err; j++) {
        /* make a compressed copy of field j */
        i = 0;
        while (s != q) {
            if (*s == ';') {
                sc++;
                s++;
                break; /* onto next field */
            }
            if (*s != ' ') {
                tmp[i++] = *s;
            }
            s++;
        }
        tmp[i] = '\0';
        err = set_forloop_element(tmp, loop, j);
    }

    if (!err && (sc != 2 || s != q)) {
        /* we've reached the reached rightmost ')' but have not
           found two semi-colons */
        err = E_PARSE;
    }

    free(tmp);

    if (!err) {
        loop->type = FOR_LOOP;
    }

    return err;
}

static int is_indexed_loop (const char *s,
                            char *lvar,
                            char **start,
                            char **stop)
{
    int n;

    /* must conform to the pattern

      lvar = start..stop

      where @start and/or @stop may be compound terms,
      with or without whitespace between terms
    */

    s += strspn(s, " ");
    n = gretl_namechar_spn(s);

    if (n > 0 && n < VNAMELEN) {
        const char *s0 = s;
        const char *p = strstr(s, "..");

        if (p != NULL) {
            s += n;
            s += strspn(s, " ");
            if (*s == '=' && *(s+1) != '=') {
                *lvar = '\0';
                strncat(lvar, s0, n);
                /* skip any space after '=' */
                s++;
                s += strspn(s, " ");
                *start = gretl_strndup(s, p - s);
                g_strchomp(*start);
                p += 2; /* skip ".." */
                p += strspn(p, " ");
                *stop = gretl_strdup(p);
                g_strchomp(*stop);
                return 1;
            }
        }
    }

    return 0;
}

static int parse_first_loopline (char *s, LOOPSET *loop,
                                 DATASET *dset)
{
    char vname[VNAMELEN];
    char *start = NULL;
    char *stop = NULL;
    int err = 0;

    /* skip preliminary string */
    while (isspace(*s)) s++;
    if (!strncmp(s, "loop", 4)) {
        s += 4;
        while (isspace(*s)) s++;
    }

    /* syntactic slop: accept "for i=lo..hi" -> "i=lo..hi" */
    if (!strncmp(s, "for ", 4) && !strchr(s, ';')) {
        s += 4;
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_first_loopline: '%s'\n", s);
#endif

    if (!strncmp(s, "foreach ", 8)) {
        err = parse_as_each_loop(loop, dset, s + 8);
    } else if (!strncmp(s, "for ", 4)) {
        err = parse_as_for_loop(loop, s + 4);
    } else if (!strncmp(s, "while ", 6)) {
        err = parse_as_while_loop(loop, s + 6);
    } else if (is_indexed_loop(s, vname, &start, &stop)) {
        err = parse_as_indexed_loop(loop, dset, vname, start, stop);
        free(start);
        free(stop);
    } else {
        /* must be a count loop, or erroneous */
        err = parse_as_count_loop(loop, dset, s);
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "parse_first_loopline: returning %d\n", err);
#endif

    return err;
}

/**
 * start_new_loop:
 * @s: loop specification line.
 * @inloop: current loop struct pointer, or %NULL.
 * @dset: dataset struct.
 * @opt: options associated with new loop.
 * @nested: location to receive info on whether a new
 * loop was created, nested within the input loop.
 * @err: location to receive error code.
 *
 * Create a new LOOPSET based on the input line; this may or
 * may not be a child of @inloop.
 *
 * Returns: loop pointer on successful completion, %NULL on error.
 */

static LOOPSET *start_new_loop (char *s, LOOPSET *inloop,
                                DATASET *dset,
                                gretlopt opt,
                                int *nested,
                                int *err)
{
    LOOPSET *loop = NULL;

    gretl_error_clear();

#if COMP_DEBUG || LOOP_DEBUG
    fprintf(stderr, "start_new_loop: inloop=%p, line='%s'\n",
            (void *) inloop, s);
#endif

    if (inloop == NULL || compile_level <= inloop->level) {
        loop = gretl_loop_new(NULL);
    } else {
        loop = gretl_loop_new(inloop);
        *nested = 1;
    }

    if (loop == NULL) {
        gretl_errmsg_set(_("Out of memory!"));
        *err = E_ALLOC;
        return NULL;
    }

#if COMP_DEBUG || LOOP_DEBUG
    fprintf(stderr, " added loop at %p (%s)\n", (void *) loop,
            (*nested)? "nested" : "independent");
#endif

    *err = parse_first_loopline(s, loop, dset);

    if (!*err) {
        *err = gretl_loop_prepare(loop);
    }

    if (*err) {
#if LOOP_DEBUG
        fprintf(stderr, "start_new_loop: aborting on error\n");
#endif
        destroy_loop_stack(loop);
        loop = NULL;
    }

    return loop;
}

#if LOOP_DEBUG
# define MAX_FOR_TIMES  10
#else
# define MAX_FOR_TIMES  50000000
#endif

static int loop_count_too_high (LOOPSET *loop)
{
    int nt = loop->iter + 1;

    if (loop->type == FOR_LOOP) {
        if (nt > MAX_FOR_TIMES) {
            gretl_errmsg_sprintf(_("Reached maximum iterations, %d"),
                                 MAX_FOR_TIMES);
            loop->err = 1;
        }
    } else {
        int maxit = libset_get_int(LOOP_MAXITER);

        if (maxit > 0 && nt > maxit) {
            gretl_errmsg_sprintf(_("Reached maximum iterations, %d"),
                                 maxit);
            gretl_errmsg_append(_("You can use \"set loop_maxiter\" "
                                  "to increase the limit"), 0);
            loop->err = 1;
        }
    }

    return loop->err;
}

/**
 * loop_condition:
 * @loop: pointer to loop commands struct.
 * @dset: data information struct.
 * @err: location to receive error code.
 *
 * Check whether a loop continuation condition is still satisfied.
 *
 * Returns: 1 to indicate looping should continue, 0 to terminate.
 */

static int loop_condition (LOOPSET *loop, DATASET *dset, int *err)
{
    int ok = 0;

    if (loop->brk) {
        /* got "break" comand */
        loop->brk = 0;
        ok = 0;
    } else if (loop->type == COUNT_LOOP || indexed_loop(loop)) {
        if (loop->iter < loop->itermax) {
            ok = 1;
            if (indexed_loop(loop) && loop->iter > 0) {
                if (loop->flags & LOOP_DECREMENT) {
                    loop->idxval -= 1;
                } else {
                    loop->idxval += 1;
                }
                if (loop->idxvar == NULL) {
                    loop->idxvar = get_user_var_by_name(loop->idxname);
                }
                if (loop->idxvar != NULL) {
                    uvar_set_scalar_fast(loop->idxvar, loop->idxval);
                } else {
                    fprintf(stderr, "loop_condition: couldn't find index variable!\n");
                    ok = 0;
                }
            }
        }
    } else if (!loop_count_too_high(loop)) {
        /* more complex forms of control (for, while) */
        if (loop->type == FOR_LOOP) {
            if (loop->iter > 0) {
                loop_delta(loop, dset, err);
            }
            ok = loop_testval(loop, dset, err);
        } else if (loop->type == WHILE_LOOP) {
            ok = loop_testval(loop, dset, err);
        }
    }

    return ok;
}

static void controller_init (controller *clr)
{
    clr->val = NADBL;
    clr->vname[0] = '\0';
    clr->uv = NULL;
    clr->vsign = 1;
    clr->expr = NULL;
    clr->genr = NULL;
    clr->subst = -1;
}

static void controller_free (controller *clr)
{
    if (clr->expr != NULL) {
        free(clr->expr);
        clr->expr = NULL;
    }
    if (clr->genr != NULL) {
        destroy_genr(clr->genr);
        clr->genr = NULL;
    }
}

static void loop_lines_init (LOOPSET *loop, int i1, int i2)
{
    int i;

    for (i=i1; i<i2; i++) {
        loop->lines[i].s = NULL;
        loop->lines[i].ci = 0;
        loop->lines[i].next = 0;
	loop->lines[i].idx = 0;
        loop->lines[i].ptr = NULL;
        loop->lines[i].flags = 0;
    }
}

static int gretl_loop_prepare (LOOPSET *loop)
{
#if HAVE_GMP
    mpf_set_default_prec(256);
#endif

    /* allocate some initial lines/commands for loop */
    loop->lines = malloc(LOOP_BLOCK * sizeof *loop->lines);

    if (loop->lines == NULL) {
        return E_ALLOC;
    } else {
        loop_lines_init(loop, 0, LOOP_BLOCK);
    }

    return 0;
}

/* See if we already have a model recorder in place for the command on
   line @lno of the loop.  If so, fetch it, otherwise create a new one
   and return it.
*/

static MODEL *get_model_record_by_line (LOOPSET *loop, int lno,
                                        gretlopt opt, int *err)
{
    model_record *mrec;
    MODEL *pmod = NULL;
    int n = loop->n_models;
    int i;

    for (i=0; i<n; i++) {
        if (lno == loop->mrecords[i].lno) {
            return loop->mrecords[i].pmod;
        }
    }

    mrec = realloc(loop->mrecords, (n + 1) * sizeof *mrec);
    if (mrec == NULL) {
        *err = E_ALLOC;
    } else {
        loop->mrecords = mrec;
        pmod = gretl_model_new();
        if (pmod == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        /* 2016-10-24: I think this is right, AC. Note
           that there's a matching "unprotect" when a loop
           is destroyed.
        */
        gretl_model_protect(pmod);
        mrec = &loop->mrecords[n];
        mrec->lno = lno;
        mrec->opt = opt;
        mrec->pmod = pmod;
        pmod->ID = n + 1;
        loop->n_models += 1;
    }

    return pmod;
}

#if HAVE_GMP
/* functions specific to "progressive" loops, which require GMP */
# include "prog_loop.c"
#endif

int model_is_in_loop (const MODEL *pmod)
{
    LOOPSET *loop = currloop;
    int i;

    while (loop != NULL) {
        for (i=0; i<loop->n_models; i++) {
            if (pmod == loop->mrecords[i].pmod) {
                return 1;
            }
        }
        loop = loop->parent;
    }

    return 0;
}

static int loop_get_structure (LOOPSET *loop)
{
    return statements_get_structure(loop->lines,
				    loop->n_lines,
				    "loop", NULL);
}

static int add_more_loop_lines (LOOPSET *loop)
{
    int nb = 1 + (loop->n_lines + 1) / LOOP_BLOCK;
    int totlines = nb * LOOP_BLOCK;
    loop_line *lines;

    /* in case we ran out of space */
    lines = realloc(loop->lines, totlines * sizeof *lines);

    if (lines == NULL) {
        return E_ALLOC;
    }

    loop->lines = lines;
    loop_lines_init(loop, loop->n_lines, totlines);

    return 0;
}

static int real_append_line (ExecState *s, LOOPSET *loop,
                             LOOPSET *child)
{
    int n = loop->n_lines;
    int err = 0;

#if LOOP_DEBUG > 1
    fprintf(stderr, "real_append_line: s->line = '%s'\n", s->line);
#endif

    if ((n + 1) % LOOP_BLOCK == 0) {
        if (add_more_loop_lines(loop)) {
            return E_ALLOC;
        }
    }

    loop->lines[n].s = gretl_strdup(s->line);

    if (loop->lines[n].s == NULL) {
        err = E_ALLOC;
    } else {
	if (strchr(s->line, '$')) {
	    loop->lines[n].flags |= LOOP_LINE_DOLLAR;
	}
	if (strchr(s->line, '@')) {
	    loop->lines[n].flags |= LOOP_LINE_AT;
	}
        if (strncmp(s->line, "catch ", 6) == 0) {
            loop->lines[n].flags |= LOOP_LINE_CATCH;
        }
        if (s->cmd->ci == PRINT) {
            if (!loop_is_progressive(loop) || strchr(s->line, '"')) {
                /* printing a literal string, not a variable's value */
                loop->lines[n].flags |= LOOP_LINE_LIT;
            }
        } else if (s->cmd->ci == RENAME || s->cmd->ci == OPEN) {
            loop_set_renaming(loop);
        } else if (s->cmd->ci == IF) {
            loop_set_has_cond(loop);
        }
        loop->lines[n].ci = s->cmd->ci;
        loop->lines[n].ptr = child;
        loop->n_lines += 1;
    }

#if LOOP_DEBUG > 1
    fprintf(stderr, "loop %p: n_lines=%d, line[%d]='%s', ci=%s\n",
            (void *) loop, loop->n_lines, n, loop->lines[n].s,
            gretl_command_word(loop->lines[n].ci));
#endif

    return err;
}

/**
 * gretl_loop_append_line:
 * @s: program execution state.
 * @dset: dataset struct.
 * @ploop: optional location to receive address of completed loop.
 *
 * Add the command line @s->line to accumulated loop buffer.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_loop_append_line_full (ExecState *s, DATASET *dset,
				 LOOPSET **ploop)
{
    LOOPSET *loop = currloop;
    LOOPSET *newloop = currloop;
    LOOPSET *child = NULL;
    int err = 0;

    warnmsg(s->prn); /* catch "end loop" if present */
    gretl_error_clear();

#if COMP_DEBUG > 1 || LOOP_DEBUG > 1
    fprintf(stderr, "gretl_loop_append_line: currloop = %p, line = '%s'\n",
            (void *) loop, s->line);
#endif

    if (!ok_in_loop(s->cmd->ci)) {
        gretl_errmsg_sprintf(_("The '%s' command is not available in loop mode"),
                             gretl_command_word(s->cmd->ci));
        destroy_loop_stack(loop);
        return E_NOTIMP;
    }

    if (s->cmd->ci == LOOP) {
        /* starting from scratch */
        char *spec = s->cmd->param;
        gretlopt opt = s->cmd->opt;
        int nested = 0;

        if (spec == NULL) {
            fprintf(stderr, "GRETL_ERROR: loop line is unparsed\n");
            err = E_DATA;
        }

#if !HAVE_GMP
        if (opt & OPT_P) {
            gretl_errmsg_set("The progressive option is not available "
                             "in this build");
            err = E_BADOPT;
        }
#endif

        if (!err) {
            newloop = start_new_loop(spec, loop, dset, opt,
                                     &nested, &err);
#if GLOBAL_TRACE || LOOP_DEBUG
            fprintf(stderr, "got LOOP: newloop at %p (err = %d)\n",
                    (void *) newloop, err);
#endif
            if (newloop == NULL) {
                return err;
            } else {
                set_loop_opts(newloop, opt);
                compile_level++;
                if (!nested) {
                    currloop = newloop;
                    return 0; /* done */
                }
                child = newloop;
            }
        }
    } else if (s->cmd->ci == ENDLOOP) {
        /* got to the end */
        compile_level--;
#if GLOBAL_TRACE || LOOP_DEBUG
        fprintf(stderr, "got ENDLOOP, compile_level now %d\n",
                compile_level);
#endif
        if (loop_has_cond(loop)) {
            err = loop_get_structure(loop);
        }
        if (!err && compile_level == 0) {
            /* set flag to run the loop */
	    if (ploop != NULL) {
		/* signal to caller that a loop is ready */
		*ploop = loop;
		loop_set_attached(loop);
	    }
            loop_execute = 1;
        } else {
            /* back up a level */
            newloop = loop->parent;
        }
    }

    if (!err && loop != NULL && s->cmd->ci != ENDLOOP) {
        err = real_append_line(s, loop, child);
    }

    if (err) {
        if (loop != NULL) {
            gretl_loop_destroy(loop);
            compile_level = 0;
        }
    } else {
        currloop = newloop;
    }

    return err;
}

int gretl_loop_append_line (ExecState *s, DATASET *dset)
{
    return gretl_loop_append_line_full(s, dset, NULL);
}

/**
 * print_loop_results:
 * @loop: pointer to loop struct.
 * @dset: data information struct.
 * @prn: gretl printing struct.
 *
 * Print out the results after completion of the loop @loop.
 * This function is not called for progressive loops.
 */

static void print_loop_results (LOOPSET *loop, const DATASET *dset,
                                PRN *prn)
{
    model_record *mrec;
    gretlopt popt;
    int i;

    for (i=0; i<loop->n_models; i++) {
        mrec = &loop->mrecords[i];
        if (model_print_deferred(mrec->opt)) {
            set_model_id(mrec->pmod, OPT_NONE);
            popt = get_printmodel_opt(mrec->pmod, mrec->opt);
            printmodel(mrec->pmod, dset, popt, prn);
        }
    }
}

static int substitute_dollar_targ (char *str, int maxlen,
                                   const LOOPSET *loop,
                                   const DATASET *dset,
                                   int *subst)
{
    char insert[32], targ[VNAMELEN + 3] = {0};
    char *p, *ins, *q, *s;
    int targlen, inslen, idx = 0;
    int incr, cumlen = 0;
    int err = 0;

#if SUBST_DEBUG
    fprintf(stderr, "subst_dollar_targ:\n original: '%s'\n", str);
#endif

    /* construct the target for substitution */

    if (loop->type == FOR_LOOP) {
        if (!gretl_is_scalar(loop->init.vname)) {
            /* nothing to substitute */
            return 0;
        }
        sprintf(targ, "$%s", loop->init.vname);
        targlen = strlen(targ);
    } else if (indexed_loop(loop)) {
        sprintf(targ, "$%s", loop->idxname);
        targlen = strlen(targ);
        idx = loop->init.val + loop->iter;
    } else {
        /* shouldn't be here! */
        return 1;
    }

#if SUBST_DEBUG
    fprintf(stderr, " target = '%s', idx = %d\n", targ, idx);
#endif

    if (strstr(str, targ) == NULL) {
        /* nothing to be done */
        return 0;
    }

    ins = insert;

    /* prepare the substitute string */

    if (loop->type == FOR_LOOP) {
        double x = gretl_scalar_get_value(loop->init.vname, NULL);

        if (na(x)) {
            strcpy(insert, "NA");
        } else {
            sprintf(insert, "%g", x);
        }
    } else if (loop->type == INDEX_LOOP) {
        sprintf(insert, "%d", idx);
    } else if (loop->type == DATED_LOOP) {
        /* note: ntolabel is 0-based */
        ntolabel(insert, idx - 1, dset);
    } else if (loop->type == EACH_LOOP) {
        ins = loop->eachstrs[idx - 1];
    }

    inslen = strlen(ins);
    incr = inslen - targlen;
    if (incr > 0) {
        /* substitution will lengthen the string */
        cumlen = strlen(str);
    }

    q = malloc(strlen(strstr(str, targ)));
    if (q == NULL) {
        err = E_ALLOC;
    }

    /* crawl along str, replacing targ with ins */

    s = str;
    while ((p = strstr(s, targ)) != NULL && !err) {
        if (is_gretl_accessor(p)) {
            s++;
            continue;
        }
        if (incr > 0) {
            cumlen += incr;
            if (cumlen >= maxlen) {
                /* substitution would cause overflow */
                err = (maxlen == VNAMELEN)? E_UNKVAR : E_TOOLONG;
                break;
            }
        }
        strcpy(q, p + targlen);
        strcpy(p, ins);
        strcpy(p + inslen, q);
        if (subst != NULL) {
            *subst = 1;
        }
        s++; /* += strlen(ins)? */
    }

    free(q);

#if SUBST_DEBUG
    fprintf(stderr, " after: '%s'\n", str);
#endif

    return err;
}

/* When re-executing a loop that has been saved onto its
   calling function, the loop index variable may have been
   destroyed, in which case it has to be recreated.
*/

static int loop_reattach_index_var (LOOPSET *loop, DATASET *dset)
{
    char genline[64];
    int err = 0;

    if (na(loop->init.val)) {
        sprintf(genline, "%s=NA", loop->idxname);
    } else {
        gretl_push_c_numeric_locale();
        sprintf(genline, "%s=%g", loop->idxname, loop->init.val);
        gretl_pop_c_numeric_locale();
    }

    err = generate(genline, dset, GRETL_TYPE_DOUBLE, OPT_Q, NULL);

    if (!err) {
        loop->idxvar = get_user_var_by_name(loop->idxname);
    }

    return err;
}

/* Called at the start of iteration for a given loop */

static int top_of_loop (LOOPSET *loop, DATASET *dset)
{
    int err = 0;

    loop->iter = 0;

    if (loop->eachname[0] != '\0') {
        err = loop_list_refresh(loop, dset);
    } else if (loop->type == INDEX_LOOP) {
        loop->init.val = controller_get_val(&loop->init, loop, dset, &err);
    } else if (loop->type == FOR_LOOP) {
        forloop_init(loop, dset, &err);
    }

    if (!err && loop->idxname[0] != '\0' && loop->idxvar == NULL) {
        err = loop_reattach_index_var(loop, dset);
    }

    if (!err && (loop->type == COUNT_LOOP || indexed_loop(loop))) {
        loop->final.val = controller_get_val(&loop->final, loop, dset, &err);
        if (na(loop->init.val) || na(loop->final.val)) {
            gretl_errmsg_set(_("error evaluating loop condition"));
            fprintf(stderr, "loop: got NA for init and/or final value\n");
            err = E_DATA;
        } else if ((loop->type == INDEX_LOOP || loop->type == DATED_LOOP) &&
                   loop_decrement(loop)) {
            loop->itermax = loop->init.val - loop->final.val + 1;
        } else {
            loop->itermax = loop->final.val - loop->init.val + 1;
            loop->flags &= ~LOOP_DECREMENT;
        }
    }

    if (!err) {
        if (indexed_loop(loop)) {
            loop->idxval = loop->init.val;
            uvar_set_scalar_fast(loop->idxvar, loop->idxval);
        }
        /* initialization, in case this loop is being run more than
           once (i.e. it's embedded in an outer loop)
        */
#if HAVE_GMP
        if (loop_is_progressive(loop)) {
            progressive_loop_zero(loop);
        } else {
            free(loop->mrecords);
            loop->mrecords = NULL;
            loop->n_models = 0;
        }
#else
        free(loop->mrecords);
        loop->mrecords = NULL;
        loop->n_models = 0;
#endif /* HAVE_GMP or not */
    }

    return err;
}

static const LOOPSET *subst_loop_in_parentage (const LOOPSET *loop)
{
    while ((loop = loop->parent) != NULL) {
        if (indexed_loop(loop) || loop->type == FOR_LOOP) break;
    }

    return loop;
}

static int
make_dollar_substitutions (char *str, int maxlen,
                           const LOOPSET *loop,
                           const DATASET *dset,
                           int *subst,
                           gretlopt opt)
{
    int err = 0;

    if (subst != NULL) {
        *subst = 0;
    }

    /* if (opt & OPT_T) we're just processing a variable name, at the top
       of a loop, so we can skip to the "parentage" bit
    */

    if (!(opt & OPT_T) && (indexed_loop(loop) || loop->type == FOR_LOOP)) {
        err = substitute_dollar_targ(str, maxlen, loop, dset, subst);
    }

    while (!err && (loop = subst_loop_in_parentage(loop)) != NULL) {
        err = substitute_dollar_targ(str, maxlen, loop, dset, subst);
    }

    return err;
}

int scalar_is_read_only_index (const char *name)
{
    const LOOPSET *loop = currloop;

    while (loop != NULL) {
        if (indexed_loop(loop) && !strcmp(name, loop->idxname)) {
            return 1;
        }
        loop = loop->parent;
    }

    return 0;
}

static int try_add_loop_genr (LOOPSET *loop,
                              int lno,
                              CMD *cmd,
                              DATASET *dset,
                              PRN *prn)
{
    GretlType gtype = cmd->gtype;
    const char *line = cmd->vstart;
    gretlopt gopt = OPT_NONE;
    int err = 0;

#if !COMPILE_FEVAL
    /* see comment on COMPILE_FEVAL in cmd_private.h */
    if (strstr(line, "feval")) {
	return 0;
    }
#endif

    if (cmd->opt & OPT_O) {
        gopt |= OPT_O;
    }

    loop->lines[lno].ptr = genr_compile(line, dset, gtype,
                                        gopt, prn, &err);
    if (err == E_EQN) {
        /* may be a non-compilable special such as "genr time",
	   a non-fatal thing */
        err = 0;
    }

    return err;
}

static int loop_print_save_model (MODEL *pmod, DATASET *dset,
                                  PRN *prn, ExecState *s)
{
    int err = pmod->errcode;

    if (!err) {
        int havename = *s->cmd->savename != '\0';
        int window = (s->cmd->opt & OPT_W) != 0;

        set_gretl_errno(0);
        if (!(s->cmd->opt & OPT_Q)) {
            gretlopt popt = get_printmodel_opt(pmod, s->cmd->opt);

            printmodel(pmod, dset, popt, prn);
        }
        attach_subsample_to_model(pmod, dset);
        s->pmod = maybe_stack_model(pmod, s->cmd, prn, &err);
        if (!err && gretl_in_gui_mode() && s->callback != NULL &&
            (havename || window)) {
            s->callback(s, s->pmod, GRETL_OBJ_EQN);
        }
    }

    return err;
}

static int loop_process_error (LOOPSET *loop, loop_line *ll,
                               int err, PRN *prn)
{
    if (loop_line_catch(ll)) {
        set_gretl_errno(err);
	if (err != E_STOP) {
	    loop->flags |= LOOP_ERR_CAUGHT;
	    err = 0;
	}
    }

    return err;
}

/* Based on the stored flags in the loop-line record, set
   or unset some flags for the command parser: this can
   reduce the amount of work the parser has to do on each
   iteration of a loop (maybe some of this obsolete?).
*/

static inline void loop_info_to_cmd (LOOPSET *loop,
                                     loop_line *ll,
                                     CMD *cmd)
{
    if (loop_is_progressive(loop)) {
        cmd->flags |= CMD_PROG;
    } else {
        cmd->flags &= ~CMD_PROG;
    }

    /* don't carry this over */
    *cmd->savename = '\0';

    if (line_has_at_sign(ll)) {
	cmd->flags &= ~CMD_NOSUB;
    } else {
	/* tell parser not to bother trying for @-substitution */
	cmd->flags |= CMD_NOSUB;
    }

    /* readjust "catch" for commands that are not being
       sent through the parser again */
    if (loop_line_catch(ll)) {
        cmd->flags |= CMD_CATCH;
    } else if (!cmd->context) {
        cmd->flags &= ~CMD_CATCH;
    }
}

/* Based on the parsed info in @cmd, maybe modify some flags in
   the current loop-line record.
*/

static inline void cmd_info_to_loop (LOOPSET *loop, int j, CMD *cmd)
{
    loop_line *ll = &loop->lines[j];

    if (line_has_at_sign(ll) && !cmd_subst(cmd)) {
        /* This loop line has not already been marked as
           free of @-substitution, but has been found to
	   be such.
        */
	ll->flags &= ~LOOP_LINE_AT;
    }

    if (cmd->flags & CMD_CATCH) {
        ll->flags |= LOOP_LINE_CATCH;
    }
}

/* We come here when the --force option has been applied to
   the "delete" command, trying to prevent deletion of the
   index variable for the loop: even --force can't allow that,
   on penalty of crashing.
*/

static int loop_check_deletion (LOOPSET *loop, const char *param,
                                PRN *prn)
{
    user_var *uv = get_user_var_by_name(param);

    if (uv != NULL) {
        while (loop != NULL) {
            if (loop->idxvar == uv) {
                pprintf(prn, _("delete %s: not allowed\n"), param);
                return 1;
            }
            loop = loop->parent;
        }
    }

    return 0;
}

/* We come here if the --force option has not been applied to
   the "delete" command, and we'll be conservative.
*/

static int loop_delete_object (LOOPSET *loop, CMD *cmd, PRN *prn)
{
    int err = 0;

    if (cmd->list != NULL && cmd->list[0] > 0) {
        /* too dangerous! */
        pputs(prn, _("You cannot delete series in this context\n"));
        err = 1;
    } else if (gretl_is_scalar(cmd->param)) {
        /* could delete loop index */
        pputs(prn, _("You cannot delete scalars in this context\n"));
        err = 1;
    } else if (loop->parent != NULL || loop->n_children > 0) {
        /* not a "singleton" loop */
        pprintf(prn, _("delete %s: not allowed\n"), cmd->param);
        err = 1;
    } else {
        /* check for compiled genrs on board: don't let these
           get screwed up by deletion of variables of any kind
        */
        int i, ok = 1;

        for (i=0; i<loop->n_lines; i++) {
            if (line_has_genr(loop, i)) {
                ok = 0;
                break;
            }
        }
        if (ok) {
            err = gretl_delete_var_by_name(cmd->param, prn);
        } else {
            pprintf(prn, _("delete %s: not allowed\n"), cmd->param);
            err = 1;
        }
    }

    return err;
}

static char *inner_errline;

static int loop_report_error (LOOPSET *loop, int err,
                              char *errline,
                              ExecState *state,
                              PRN *prn)
{
    int fd = gretl_function_depth();

    if (fd > 0 && inner_errline != NULL) {
        errline = inner_errline;
    }

    if (err) {
        if (fd == 0) {
            errmsg(err, prn);
            if (errline != NULL && *errline != '\0') {
                pprintf(prn, ">> %s\n", errline);
            }
        }
    } else if (loop->err) {
        if (fd == 0) {
            errmsg(loop->err, prn);
        }
        err = loop->err;
    }

    if (fd > 0 && err && errline != NULL && *errline != '\0') {
        strcpy(state->line, errline);
    }

    return err;
}

static void loop_reset_error (void)
{
    if (inner_errline != NULL) {
        free(inner_errline);
        inner_errline = NULL;
    }
}

#if !HAVE_GMP

static int model_command_post_process (ExecState *s,
                                       DATASET *dset,
                                       LOOPSET *loop,
                                       int j)
{
    int moderr = check_gretl_errno();
    gretlopt opt = s->cmd->opt;
    int err = 0;

    if (moderr) {
        if (model_print_deferred(opt)) {
            err = moderr;
        } else {
            errmsg(moderr, s->prn);
        }
    } else if (model_print_deferred(opt)) {
        MODEL *pmod = get_model_record_by_line(loop, j, opt, &err);

        if (!err) {
            swap_models(s->model, pmod);
            pmod->ID = j + 1;
            set_as_last_model(pmod, GRETL_OBJ_EQN);
            model_count_minus(NULL);
        }
    } else {
        loop_print_save_model(s->model, dset, s->prn, s);
    }

    return err;
}

#endif /* !HAVE_GMP */

/* loop_reset_uvars(): called on exit from a function onto
   which one or more "compiled" loops have been attached.
   The point is to reset to NULL the stored addresses of
   any "uservars" that have been recorded in the context
   of the loop, since in general on a subsequent invocation
   of the function a variable of a given name will occupy a
   different memory address. A reset to NULL will force a
   new lookup of these variables by name, both within "genr"
   and within the loop machinery.
*/

void loop_reset_uvars (LOOPSET *loop)
{
    controller *clrs[4] = {
        &loop->init, &loop->test,
        &loop->delta, &loop->final
    };
    int i;

    for (i=0; i<loop->n_children; i++) {
        loop_reset_uvars(loop->children[i]);
    }

    /* stored references within "genrs" */
    if (loop->lines != NULL) {
        for (i=0; i<loop->n_lines; i++) {
            if (line_has_genr(loop, i)) {
                genr_reset_uvars(loop->lines[i].ptr);
            }
        }
    }

    /* stored refs in controllers? */
    for (i=0; i<4; i++) {
        if (clrs[i]->genr != NULL) {
            genr_reset_uvars(clrs[i]->genr);
        }
        clrs[i]->uv = NULL;
    }

    /* another (possibly) stored reference */
    loop->idxvar = NULL;
}

/* We come here only if @ll has the LOOP_LINE_DOLLAR flag
   set. We scrub this flag if '$' turns out not to call
   string substitution.
*/

static int check_for_dollar_subst (loop_line *ll,
				   char *line,
				   LOOPSET *loop,
				   DATASET *dset)
{
    int err = 0;
    int subst = 0;

    /* handle loop-specific $-string substitution */
    err = make_dollar_substitutions(line, MAXLINE, loop,
				    dset, &subst, OPT_NONE);
    if (!err && !subst) {
	ll->flags &= ~LOOP_LINE_DOLLAR;
    } else if (subst) {
	ll->flags |= LOOP_LINE_NOCOMP;
    }

    return err;
}

static int abort_loop_execution (ExecState *s)
{
    *s->cmd->savename = '\0';
    gretl_cmd_destroy_context(s->cmd);
    return E_STOP;
}

static int block_model (CMD *cmd)
{
    return cmd->ci == END &&
        (!strcmp(cmd->param, "mle") ||
         !strcmp(cmd->param, "nls") ||
         !strcmp(cmd->param, "gmm"));
}

#define LTRACE 0

int gretl_loop_exec (ExecState *s, DATASET *dset)
{
    LOOPSET *loop = NULL;
    char *line = s->line;
    CMD *cmd = s->cmd;
    PRN *prn = s->prn;
    char *currline = NULL;
    char *showline = NULL;
    int indent0;
    int gui_mode, echo;
    int show_activity = 0;
    int prev_messages;
#if HAVE_GMP
    int progressive;
#endif
    int err = 0;

    if (s->loop) {
        loop = s->loop;
	s->loop = NULL; /* avoid stale data */
    } else {
        loop = currloop;
    }

    /* for the benefit of the caller: register the fact that execution
       of this loop is already under way */
    loop_execute = 0;

    if (loop == NULL) {
        pputs(prn, "Got a NULL loop\n");
        set_loop_off();
        return E_DATA;
    }

    gui_mode = gretl_in_gui_mode();
    echo = gretl_echo_on();
    prev_messages = gretl_messages_on();
    if (!prev_messages && loop_is_verbose(loop)) {
        set_gretl_messages(1);
    }
    indent0 = gretl_if_state_record();
    set_loop_on();
#if HAVE_GMP
    progressive = loop_is_progressive(loop);
#endif

#if LOOP_DEBUG
    fprintf(stderr, "loop_exec: loop = %p\n", (void *) loop);
#endif

    err = top_of_loop(loop, dset);

    if (!err) {
        if (loop_is_renaming(loop)) {
            loop_renaming = 1;
        }
        if (gui_mode) {
            show_activity = show_activity_func_installed();
        }
    }

    while (!err && loop_condition(loop, dset, &err)) {
        /* respective iterations of a given loop */
        int j;

#if LOOP_DEBUG > 1
        fprintf(stderr, "*** top of loop: iter = %d\n", loop->iter);
#endif
        if (gui_mode && loop->iter % 10 == 0 && get_user_stop()) {
            /* the user called a halt */
            err = abort_loop_execution(s);
            break;
        }

        for (j=0; j<loop->n_lines && !err; j++) {
            /* execute commands on this iteration */
            loop_line *ll = &loop->lines[j];
            int compiled = line_is_compiled(ll);
            int ci = ll->ci;
            int parse = 1;

            loop_renaming = loop_is_renaming(loop);
            if (loop_renaming) {
                set_non_compilable(ll);
            }

            currline = ll->s;
            if (compiled) {
                /* just for "echo" purposes */
                showline = currline;
            } else {
                /* line may be modified below */
                showline = strcpy(line, currline);
            }

#if LTRACE || (LOOP_DEBUG > 1)
            fprintf(stderr, "iter=%d, j=%d, line='%s', ci=%s, compiled=%d\n",
                    loop->iter, j, showline, gretl_command_word(ci),
                    compiled);
#endif

	    if (gretl_if_state_true() && (ci == ELSE || ci == ELIF)) {
		if (ll->next > 0) {
		    j = ll->next - 1;
		}
		continue;
	    }

	    /* deal with "compiled" cases */
	    if (ci == GENR && compiled) {
                if (echo && loop_is_verbose(loop)) {
                    pprintf(prn, "? %s\n", showline);
                }
                err = execute_genr(ll->ptr, dset, prn);
                if (ll->next > 0) {
                    j = ll->next - 1;
                }
	    } else if ((ci == IF || ci == ELIF) && compiled) {
                cmd->ci = ci;
                if (s->cmd->vstart != NULL) {
                    free(s->cmd->vstart);
                    s->cmd->vstart = NULL;
                }
                flow_control(s, dset, &ll->ptr);
                if (cmd->err) {
                    /* we hit an error evaluating the if state */
                    err = cmd->err;
                } else if (gretl_if_state_false() && ll->next > 0) {
		    j = ll->next - 1;
                }
	    } else if (ci == ELSE || ci == ENDIF) {
		/* these don't need compilation */
		cmd->ci = ci;
		flow_control(s, dset, &ll->ptr);
		if (gretl_if_state_false() && ll->next > 0) {
		    j = ll->next - 1;
		}
            } else if (ci == BREAK || ci == CONTINUE) {
                if (ci == BREAK) {
                    loop->brk = 1;
                } else {
                    loop->cont = 1;
                }
                break;
	    } else if (ci == LOOP) {
                cmd->ci = ci;
                goto child_loop;
            }

	    if (err) {
		goto handle_err;
	    } else if (compiled) {
		continue;
	    }

	    /* check for $-substitution? */
	    if (line_has_dollar(ll)) {
		err = check_for_dollar_subst(ll, line, loop, dset);
		if (err) {
		    goto handle_err;
		}
	    }

            /* transcribe saved loop info -> cmd */
            loop_info_to_cmd(loop, ll, cmd);

	    if (ci == GENR) {
		/* Looking at a "genr" line that's not yet compiled:
		   it might not really be a case of genr, for example
		   part of a "restrict" block.
		*/
		err = parse_command_line(s, dset, NULL);
		if (err) {
		    goto handle_err;
		} else if (cmd->ci != GENR) {
		    ci = ll->ci = cmd->ci;
		} else if (s->cmd->flags & CMD_SUBST) {
		    set_non_compilable(ll);
		}
		parse = 0;
	    }

	    if (ci == GENR) {
		if (may_be_compilable(ll)) {
		    err = try_add_loop_genr(loop, j, cmd, dset, prn);
		    if (ll->ptr == NULL && !err) {
                        /* fallback */
			set_non_compilable(ll);
                        err = generate(cmd->vstart, dset, cmd->gtype,
                                       cmd->opt, prn);
                    }
		} else {
		    /* string substitution or a "genr special" */
		    if (!err) {
			if (!loop_is_verbose(loop)) {
			    cmd->opt |= OPT_Q;
			}
			err = generate(cmd->vstart, dset, cmd->gtype,
				       cmd->opt, prn);
		    }
		}
		if (err) {
		    goto handle_err;
		} else {
                    if (ll->next > 0) {
                        j = ll->next - 1;
                    }
		    continue;
		}
	    } else if (ci == IF || ci == ELIF) {
		/* looking at a conditional that's not yet compiled */
		if (may_be_compilable(ll)) {
		    err = parse_command_line(s, dset, &ll->ptr);
		    if (s->cmd->flags & CMD_SUBST) {
			set_non_compilable(ll);
		    }
		} else {
		    err = parse_command_line(s, dset, NULL);
		}
		if (gretl_if_state_false() && ll->next > 0) {
		    j = ll->next - 1;
		}
		continue;
	    }

	    if (prog_cmd_started(loop, j)) {
                cmd->ci = ci;
                if (loop_line_nosub(ll)) {
                    parse = 0;
                }
            }

            if (parse && !err) {
                err = parse_command_line(s, dset, NULL);
#if LOOP_DEBUG > 1
                fprintf(stderr, "    after parse: '%s', ci=%s\n", line,
			gretl_command_word(cmd->ci));
                fprintf(stderr, "    cmd->savename = '%s'\n", cmd->savename);
                fprintf(stderr, "    err from parse_command_line: %d\n", err);
#endif
            }

        handle_err:

            if (err) {
                cmd_info_to_loop(loop, j, cmd);
                cmd->err = err = loop_process_error(loop, ll, err, prn);
                if (err) {
                    break;
                } else {
                    continue;
                }
            } else {
                gretl_exec_state_transcribe_flags(s, cmd);
                cmd_info_to_loop(loop, j, cmd);
            }

            if (echo && loop_is_verbose(loop)) {
		gretl_echo_command(cmd, showline, prn);
            }

            /* now branch based on the command index: some commands
               require special treatment in loop context
            */

        child_loop:

            if (cmd->ci == LOOP) {
                currloop = ll->ptr;
                if (loop_is_attached(loop) && currloop != NULL) {
                    loop_set_attached(currloop);
                }
                err = gretl_loop_exec(s, dset);
            } else if (cmd->ci == FUNCRET) {
                /* The following clause added 2016-11-20: just in case
                   the return value is, or references, an automatic
                   loop index scalar.
                */
                loop->flags &= ~LOOP_DELVAR;
                err = set_function_should_return(line);
                loop->brk = 1;
                break;
            } else if (cmd->ci == DELEET && !(cmd->opt & (OPT_F | OPT_T))) {
                err = loop_delete_object(loop, cmd, prn);
#if HAVE_GMP
            } else if (progressive && handle_prog_command(loop, j, cmd, &err)) {
                ; /* OK, or not */
#endif
            } else {
                /* send command to the regular processor */
                int catch = cmd->flags & CMD_CATCH;

                if (cmd->ci == DELEET && cmd->param != NULL) {
                    /* don't delete loop indices! */
                    err = loop_check_deletion(loop, cmd->param, prn);
                }
                if (!err) {
                    err = gretl_cmd_exec(s, dset);
                }
                if (catch) {
                    /* ensure "catch" hasn't been scrubbed */
                    cmd->flags |= CMD_CATCH;
                }
                if (!err && plain_model_ci(cmd->ci)) {
                    err = model_command_post_process(s, dset, loop, j);
                } else if (!err && !check_gretl_errno() && block_model(cmd)) {
                    /* NLS, etc. */
                    loop_print_save_model(s->model, dset, prn, s);
                }
                if (ll->next > 0) {
                    j = ll->next - 1;
                }
            }
            if (err && (cmd->flags & CMD_CATCH)) {
                set_gretl_errno(err);
                cmd->flags ^= CMD_CATCH;
                err = 0;
            }
        } /* end execution of commands within loop */

        if (err || loop->brk || loop->cont) {
            /* 2022-02-06: gretl_if_state_clear() was called if @err */
            gretl_if_state_reset(indent0);
            if (loop->cont) {
                loop->cont = 0;
            }
        } else {
            err = gretl_if_state_check(indent0);
        }

        if (!err && !loop->brk) {
            loop->iter += 1;
            if (show_activity && (loop->iter % 10 == 0)) {
                show_activity_callback();
            }
        }

        if (err && inner_errline == NULL) {
            inner_errline = gretl_strdup(currline);
        }
    } /* end iterations of loop */

    cmd->flags &= ~CMD_NOSUB;

    if (loop->brk) {
        /* turn off break flag */
        loop->brk = 0;
    }

    if (err || loop->err) {
        err = loop_report_error(loop, err, currline, s, prn);
    }

    if (!err && loop->iter > 0) {
#if HAVE_GMP
        if (progressive) {
            progressive_loop_finalize(loop, dset, prn);
        } else if (loop->mrecords != NULL) {
            print_loop_results(loop, dset, prn);
        }
#else
        if (loop->mrecords != NULL) {
            print_loop_results(loop, dset, prn);
        }
#endif
    }

    if (loop->n_models > 0) {
        /* we need to update models[0] */
        GretlObjType type;
        void *ptr = get_last_model(&type);
        MODEL *pmod;
        int i;

        if (type == GRETL_OBJ_EQN && s->model != ptr) {
            pmod = loop->mrecords[loop->n_models - 1].pmod;
            swap_models(s->model, pmod);
            set_as_last_model(s->model, GRETL_OBJ_EQN);
        }
        for (i=0; i<loop->n_models; i++) {
            pmod = loop->mrecords[i].pmod;
            gretl_model_unprotect(pmod);
            gretl_model_free(pmod);
        }
    }

    if (err && gretl_function_depth() > 0) {
        ; /* leave 'line' alone */
    } else if (line != NULL) {
        *line = '\0';
    }

    /* be sure to clear some loop-special parser flags */
    cmd->flags &= ~CMD_PROG;

    if (err) {
        err = process_command_error(s, err);
    }

    set_gretl_messages(prev_messages);

    if (loop->parent == NULL) {
        /* reached top of stack: clean up */
        currloop = NULL;
        loop_renaming = 0;
        set_loop_off();
        loop_reset_error();
        if (loop_is_attached(loop)) {
            /* prevent destruction of saved loop */
            loop = NULL;
        }
        gretl_loop_destroy(loop);
    }

    return err;
}
