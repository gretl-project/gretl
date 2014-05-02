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

/* syntax tree evaluator for 'genr' and related commands */

#include "genparse.h"
#include "monte_carlo.h"
#include "gretl_string_table.h"
#include "matrix_extra.h"
#include "usermat.h"
#include "uservar.h"
#include "gretl_bfgs.h"
#include "gretl_fft.h"
#include "gretl_panel.h"
#include "kalman.h"
#include "libset.h"
#include "version.h"
#include "csvdata.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
#endif

#ifdef USE_RLIB
# include "gretl_foreign.h"
#endif

#include <errno.h>

#if GENDEBUG
# define EDEBUG GENDEBUG
# define LHDEBUG GENDEBUG
#else
# define EDEBUG 0
# define LHDEBUG 0
#endif

#define SCALARS_ENSURE_FINITE 1 /* debatable, but watch out for read/write */
#define SERIES_ENSURE_FINITE 1  /* debatable */

#define is_aux_node(n) (n != NULL && (n->flags & AUX_NODE))
#define is_tmp_node(n) (n != NULL && (n->flags & TMP_NODE))

#define nullmat_ok(f) (f == F_ROWS || f == F_COLS || f == F_DET || \
		       f == F_LDET || f == F_DIAG || f == F_TRANSP || \
		       f == F_VEC || f == F_VECH || f == F_UNVECH || \
		       f == F_DATAOK)

#define dataset_dum(n) (n->t == DUM && n->v.idnum == DUM_DATASET)

#define postfix_node(n) (n->t == UNUM_P || n->t == UNUM_M)

#define ok_list_node(n) (n->t == LIST || n->t == WLIST || n->t == NUM || \
			 n->t == MAT || n->t == EMPTY || \
			 (n->t == VEC && n->vnum >= 0))

#define uscalar_node(n) ((n->t == NUM && n->vname != NULL) || postfix_node(n))

#define umatrix_node(n) (n->t == MAT && n->vname != NULL)
#define ubundle_node(n) (n->t == BUNDLE && n->vname != NULL)
#define ulist_node(n)   (n->t == LIST && n->vname != NULL)
#define ustring_node(n) (n->t == STR && n->vname != NULL)
#define useries_node(n) (n->t == VEC && n->vnum >= 0)
#define uvar_node(n) (n->vname != NULL)

#define scalar_matrix_node(n) (n->t == MAT && gretl_matrix_is_scalar(n->v.m))
#define scalar_node(n) (n->t == NUM || scalar_matrix_node(n))

#define stringvec_node(n) (n->flags & SVL_NODE)

#define empty_or_num(n) (n == NULL || n->t == EMPTY || n->t == NUM)
#define empty_or_string(n) (n == NULL || n->t == EMPTY || n->t == STR)
#define null_or_empty(n) (n == NULL || n->t == EMPTY)

#define ok_bundled_type(t) (t == NUM || t == STR || t == MAT || \
			    t == VEC || t == BUNDLE || t == U_ADDR) 

#define lhscalar(p) (p->flags & P_LHSCAL)
#define lhlist(p) (p->flags & P_LHLIST)
#define lhstr(p) (p->flags & P_LHSTR)
#define lhbundle(p) (p->flags & P_LHBUN)

#define compiled(p) (p->flags & P_EXEC)

static void parser_init (parser *p, const char *str, DATASET *dset, 
			 PRN *prn, int flags);
static void parser_reinit (parser *p, DATASET *dset, PRN *prn);
static void printnode (NODE *t, parser *p);
static NODE *eval (NODE *t, parser *p);
static void node_type_error (int ntype, int argnum, int goodt, 
			     NODE *bad, parser *p);
static int *node_get_list (NODE *n, parser *p);
static void reattach_series (NODE *n, parser *p);

static const char *typestr (int t)
{
    switch (t) {
    case NUM:
	return "scalar";
    case VEC:
	return "series";
    case MAT:
	return "matrix";
    case UMAT:
	return "user matrix";
    case STR:
	return "string";
    case U_ADDR:
	return "address";
    case LIST:
	return "list";
    case BUNDLE:
	return "bundle";
    case EMPTY:
	return "empty";
    default:
	return "?";
    }
}

static void free_mspec (matrix_subspec *spec)
{
    if (spec != NULL) {
	free(spec->rslice);
	free(spec->cslice);
	free(spec);
    }
}

#if EDEBUG

static void print_tree (NODE *t, parser *p, int level)
{
    if (t == NULL) {
	return;
    }

    if (bnsym(t->t)) {
	int i;

	level++;
	for (i=0; i<t->v.bn.n_nodes; i++) {
	    print_tree(t->v.bn.n[i], p, level);
	}
    } else if (b3sym(t->t)) {
	level++;
	print_tree(t->v.b3.l, p, level);
	print_tree(t->v.b3.m, p, level);
	print_tree(t->v.b3.r, p, level);
    } else if (b2sym(t->t)) {
	level++;
	print_tree(t->v.b2.l, p, level);
	print_tree(t->v.b2.r, p, level);
    } else if (b1sym(t->t)) {
	level++;
	print_tree(t->v.b1.b, p, level);
    }

    fprintf(stderr, "%d: node at %p (type %03d, %s, flags %d), vname='%s'\n",
	    level, (void *) t, t->t, getsymb(t->t, NULL), t->flags, t->vname);
}

#endif

/* determine if node @n is located anywhere in tree @t */

static int in_tree (NODE *t, NODE *n)
{
    if (t == NULL || n == NULL) {
	return 0;
    }

    if (n == t) {
	return 1;
    }

    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    if (in_tree(t->v.bn.n[i], n)) {
		return 1;
	    }
	}
    } else if (b3sym(t->t)) {
	if (in_tree(t->v.b3.l, n)) return 1;
	if (in_tree(t->v.b3.m, n)) return 1;
	if (in_tree(t->v.b3.r, n)) return 1;
    } else if (b2sym(t->t)) {
	if (in_tree(t->v.b2.l, n)) return 1;
	if (in_tree(t->v.b2.r, n)) return 1;
    } else if (b1sym(t->t)) {
	if (in_tree(t->v.b1.b, n)) return 1;
    }

    return 0;
}

static void free_tree (NODE *t, parser *p, const char *msg)
{
    if (t == NULL) {
	return;
    }

#if EDEBUG
    fprintf(stderr, "%-8s: starting with t at %p (type %03d, %s)\n", msg, 
	    (void *) t, t->t, getsymb(t->t, NULL));
#endif

    /* free recursively */
    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    free_tree(t->v.bn.n[i], p, msg);
	}
	free(t->v.bn.n);
    } else if (b3sym(t->t)) {
	free_tree(t->v.b3.l, p, msg);
	free_tree(t->v.b3.m, p, msg);
	free_tree(t->v.b3.r, p, msg);
    } else if (b2sym(t->t)) {
	free_tree(t->v.b2.l, p, msg);
	free_tree(t->v.b2.r, p, msg);
    } else if (b1sym(t->t)) {
	free_tree(t->v.b1.b, p, msg);
    } 

#if EDEBUG
    fprintf(stderr, "%-8s: freeing node at %p (type %03d, %s, flags = %d)\n", msg, 
	    (void *) t, t->t, getsymb(t->t, NULL), t->flags);
#endif

    if (is_tmp_node(t)) {
#if EDEBUG
	fprintf(stderr, " tmp node: freeing attached data\n");
#endif
	if (t->t == VEC) {
	    free(t->v.xvec);
	} else if (t->t == LIST || t->t == IVEC) {
	    free(t->v.ivec);
	} else if (t->t == MAT) {
	    gretl_matrix_free(t->v.m);
	} else if (t->t == MSPEC) {
	    free_mspec(t->v.mspec);
	} else if (t->t == BUNDLE) {
	    gretl_bundle_destroy(t->v.b);
	}
    }

    if (t->t == UOBJ || t->t == WLIST) {
	free(t->v.str);
    } else if (t->t == STR && t->vname == NULL) {
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

static void parser_aux_init (parser *p)
{
    p->aux = NULL;
    p->n_aux = 0;
    p->aux_i = 0;
}

void parser_free_aux_nodes (parser *p)
{
    int i;

    if (p->aux != NULL) {
	for (i=0; i<p->n_aux; i++) {
	    if (p->aux[i] != p->ret) {
		free_tree(p->aux[i], p, "Aux");
	    }
	}
	free(p->aux);
    }
}

static NODE *newmdef (int k)
{  
    NODE *n = new_node(MDEF);

    if (n != NULL) {
	int i;

	if (k > 0) {
	    n->v.bn.n = malloc(k * sizeof n);
	    if (n->v.bn.n != NULL) {
		for (i=0; i<k; i++) {
		    n->v.bn.n[i] = NULL;
		}
	    } else {
		free(n);
		n = NULL;
	    }
	} else {
	    n->v.bn.n = NULL;
	}
	if (n != NULL) {
	    n->v.bn.n_nodes = k;
	}
    }

    return n;
}

/* new node to hold array of doubles */

static NODE *newvec (int n, int tmp)
{  
    NODE *b = new_node(VEC);

    if (b != NULL) {
	int i;

	b->flags = (tmp)? TMP_NODE : 0;
	b->v.xvec = NULL;
	if (n > 0) {
	    b->v.xvec = malloc(n * sizeof *b->v.xvec);
	    if (b->v.xvec == NULL) {
		free(b);
		b = NULL;
	    } else {
		for (i=0; i<n; i++) {
		    b->v.xvec[i] = NADBL;
		}
	    }		
	}
    }

    return b;
}

/* new node to hold array of ints */

static NODE *newivec (int n, int type)
{  
    NODE *b = new_node(type);

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

static NODE *newmat (int tmp)
{  
    NODE *n = new_node(MAT);

    if (n != NULL) {
	n->flags = (tmp)? TMP_NODE : 0;
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

static NODE *newlist (void)
{  
    NODE *n = new_node(LIST);

    if (n != NULL) {
	n->flags = TMP_NODE;
	n->v.str = NULL;
    }    

    return n;
}

static int node_allocate_matrix (NODE *t, int m, int n, parser *p)
{
    if (m == 0 || n == 0) {
	t->v.m = gretl_null_matrix_new();
    } else {
	t->v.m = gretl_matrix_alloc(m, n);
    }
    if (t->v.m == NULL) {
	p->err = E_ALLOC;
    }

    return p->err;
}

static NODE *newstring (void)
{  
    NODE *n = new_node(STR);

    if (n != NULL) {
	n->flags = TMP_NODE;
	n->v.str = NULL;
    }

    return n;
}

static NODE *newbundle (void)
{  
    NODE *n = new_node(BUNDLE);

    if (n != NULL) {
	n->flags = TMP_NODE;
	n->v.b = NULL;
    }

    return n;
}

/* push an auxiliary evaluation node onto the stack of
   such nodes */

static int add_aux_node (parser *p, NODE *t)
{
    NODE **aux;

    aux = realloc(p->aux, (p->n_aux + 1) * sizeof *aux);
    
    if (aux == NULL) {
	p->err = E_ALLOC;
    } else {
	t->flags |= AUX_NODE;
	aux[p->n_aux] = t;
	p->aux = aux;
	p->aux_i = p->n_aux;
	p->n_aux += 1;
    }

    return p->err;
}

/* get an auxiliary node: if starting from scratch we allocate
   a new node, otherwise we look up an existing one */

static NODE *get_aux_node (parser *p, int t, int n, int tmp)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (t == NUM) {
	    ret = newdbl(NADBL);
	} else if (t == VEC) {
	    ret = newvec(n, tmp);
	} else if (t == IVEC) {
	    ret = newivec(n, IVEC);
	} else if (t == LIST) {
	    ret = newivec(n, LIST);
	} else if (t == MAT) {
	    ret = newmat(tmp);
	} else if (t == MSPEC) {
	    ret = newmspec();
	} else if (t == MDEF) {
	    ret = newmdef(n);
	} else if (t == LIST) {
	    ret = newlist();
	} else if (t == STR) {
	    ret = newstring();
	} else if (t == BUNDLE) {
	    ret = newbundle();
	} else if (t == EMPTY) {
	    ret = newempty();
	}

	if (ret == NULL) {
	    p->err = (t == 0)? E_DATA : E_ALLOC;
	} else if (add_aux_node(p, ret)) {
	    free_tree(ret, p, "On error");
	    ret = NULL;
	} 
    } else if (p->aux == NULL) {
	p->err = E_DATA;
    } else {
	while (p->aux[p->aux_i] == NULL) {
 	    p->aux_i += 1;
 	}
 	ret = p->aux[p->aux_i];
	p->aux_i += 1;
    }

    return ret;
}

static NODE *aux_scalar_node (parser *p)
{
    return get_aux_node(p, NUM, 0, 0);
}

static void no_data_error (parser *p)
{
    p->err = E_NODATA;
}

static NODE *aux_vec_node (parser *p, int n)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, VEC, n, 1);
    }
}

static NODE *aux_series_node (parser *p, int n)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, VEC, n, 0);
    }
}

static NODE *aux_ivec_node (parser *p, int n)
{
    return get_aux_node(p, IVEC, n, 1);
}

static NODE *aux_list_node (parser *p)
{
    if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
	return NULL;
    } else {
	return get_aux_node(p, LIST, 0, 1);
    }
}

static NODE *aux_matrix_node (parser *p)
{
    return get_aux_node(p, MAT, 0, 1);
}

static NODE *matrix_pointer_node (parser *p)
{
    return get_aux_node(p, MAT, 0, 0);
}

static NODE *aux_mspec_node (parser *p)
{
    return get_aux_node(p, MSPEC, 0, 0);
}

static NODE *aux_string_node (parser *p)
{
    return get_aux_node(p, STR, 0, 0);
}

static NODE *aux_bundle_node (parser *p)
{
    return get_aux_node(p, BUNDLE, 0, 0);
}

static NODE *aux_empty_node (parser *p)
{
    return get_aux_node(p, EMPTY, 0, 0);
}

static NODE *aux_any_node (parser *p)
{
    return get_aux_node(p, 0, 0, 0);
}

static void eval_warning (parser *p, int op, int errnum)
{
    if (!check_gretl_warning()) {
	const char *s = NULL;
	const char *w = "";

	if (op == B_POW) {
	    w = "pow";
	} else if (op == F_LOG) {
	    w = "log";
	} else if (op == F_SQRT) {
	    w = "sqrt";
	} else if (op == F_EXP) {
	    w = "exp";
	} else if (op == F_GAMMA) {
	    w = "gammafun";
	} else if (op == F_LNGAMMA) {
	    w = "lgamma";
	} else if (op == F_DIGAMMA) {
	    w = "digamma";
	}

	if (errnum) {
	    s = strerror(errnum);
	}

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
	    x, y, op, getsymb(op, NULL));
#endif

#if SERIES_ENSURE_FINITE
    if (targ == VEC) {
	/* this may be questionable */
	if (isnan(x)) x = NADBL;
	if (isnan(y)) y = NADBL;
    }
#endif

    /* assignment */
    if (op == B_ASN) {
	return y;
    }    

    /* testing for presence of NAs? */
    if ((p->flags & P_NATEST) && (na(x) || na(y))) {
	return NADBL;
    }

    /* 0 * NA = NA * 0 = 0 */
    if (op == B_MUL && ((na(x) && y == 0) || (na(y) && x == 0))) {
	return 0;
    }

    /* logical OR: if x is valid and non-zero, ignore NA for y */
    if (op == B_OR && !na(x) && x != 0) {
	return 1.0;
    }  

    /* otherwise NA propagates to the result */
    if (na(x) || na(y)) {
	return NADBL;
    }

    errno = 0;

    switch (op) {
    case B_ADD: 
	return x + y;
    case B_SUB: 
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
    case INC:
	return x + 1.0;
    case DEC:
	return x - 1.0;
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

#define randgen(f) (f == F_RANDGEN || f == F_MRANDGEN || f == F_RANDGEN1)

static int check_dist_count (int d, int f, int *np, int *argc)
{
    int err = 0;

    *np = *argc = 0;

    if (d == D_UNIFORM || d == D_UDISCRT) {
	/* only RANDGEN is supported */
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
    } else if (d == D_WEIBULL) {
	/* inverse cdf not supported */
	if (f == F_INVCDF) {
	    err = E_INVARG;
	} else {
	    *np = 2; /* shape, scale */
	}
    } else if (d == D_GED) {
	/* GED: critical values not supported */
	if (f == F_CRIT) {
	    err = E_INVARG;
	} else {
	    *np = 1; /* shape */
	}	
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
	/* randgen only */
	if (randgen(f)) {
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

static double *full_length_NA_vec (parser *p)
{
    double *x = malloc(p->dset->n * sizeof *x);
    int t;

    if (x == NULL) {
	p->err = E_ALLOC;
    } else {
	for (t=0; t<p->dset->n; t++) {
	    x[t] = NADBL;
	}
    }

    return x;
}

/* @parm contains an array of scalar parameters;
   @argvec contains a series of argument values.
*/

static double *series_pdist (int f, int d, 
			     double *parm, int np,
			     const double *argvec,
			     parser *p)
{
    double *xvec;
    int t;

    xvec = full_length_NA_vec(p);
    if (xvec == NULL) {
	return NULL;
    }

    if (f == F_PDF) {
	/* fast treatment, for pdf only at this point */
	int n = sample_size(p->dset);

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    xvec[t] = argvec[t];
	}
	gretl_fill_pdf_array(d, parm, xvec + p->dset->t1, n);
    } else {
	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    xvec[t] = scalar_pdist(f, d, parm, np, argvec[t], p);
	}
    }

    return xvec;
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

static double node_get_scalar (NODE *n, parser *p)
{
    if (n->t == NUM) {
	return n->v.xval;
    } else if (scalar_matrix_node(n)) {
	return n->v.m->val[0];
    } else {
	p->err = E_INVARG;
	return NADBL;
    }
}

static int node_get_int (NODE *n, parser *p)
{
    double x = node_get_scalar(n, p);

    if (p->err == 0 && (na(x) || fabs(x) > INT_MAX)) {
	p->err = E_INVARG;
	return -1;
    } else {
	return (int) x;
    }
}

static NODE *DW_node (NODE *r, parser *p)
{
    NODE *s, *e, *ret = NULL;
    int i, parm[2] = {0};

    for (i=0; i<2 && !p->err; i++) {
	s = r->v.bn.n[i+1];
	if (scalar_node(s)) {
	    parm[i] = node_get_int(s, p);
	} else {
	    e = eval(s, p);
	    if (!p->err) {
		if (scalar_node(e)) {
		    parm[i] = node_get_int(e, p);
		} else {
		    p->err = E_INVARG;
		}
	    }
	}
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

static NODE *eval_urcpval (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *s, *e, *r = n->v.b1.b;
	int i, m = r->v.bn.n_nodes;
	double x[4];

	if (m != 4) {
	    p->err = E_INVARG;
	}

	/* need double, int, int, int */

	for (i=0; i<4 && !p->err; i++) {
	    s = r->v.bn.n[i];
	    if (scalar_node(s)) {
		if (i == 0) {
		    x[i] = node_get_scalar(s, p);
		} else {
		    x[i] = node_get_int(s, p);
		}
	    } else {
		e = eval(s, p);
		if (!p->err) {
		    if (scalar_node(e)) {
			if (i == 0) {
			    x[i] = node_get_scalar(e, p);
			} else {
			    x[i] = node_get_int(e, p);
			}
		    } else {
			p->err = E_TYPES;
		    }
		    if (!reusable(p)) {
			free_tree(s, p, "Pdist");
			r->v.bn.n[i] = NULL;
		    }		    
		}
	    }
	}

	if (!p->err) {
	    double tau = x[0];
	    int nobs = (int) x[1];
	    int niv = (int) x[2];
	    int itv = (int) x[3];

	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		ret->v.xval = get_urc_pvalue(tau, nobs, niv, 
					     itv, OPT_NONE);
	    }
	}
    } else {
	ret = aux_any_node(p);
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
	    e = eval(n->v.bn.n[i+1], p);
	    if (p->err) {
		break;
	    } 
	    if (scalar_node(e)) {
		if (i == 0) {
		    rho = node_get_scalar(e, p);
		} else {
		    args[i-1] = node_get_scalar(e, p);
		}
	    } else if (i == 1) {
		if (e->t == VEC) {
		    avec = e->v.xvec;
		} else if (e->t == MAT) {
		    amat = e->v.m;
		}
	    } else if (i == 2) {
		if (e->t == VEC) {
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
		ret = aux_vec_node(p, p->dset->n);
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

/* return a node containing the evaluated result of a
   probability distribution function */

static NODE *eval_pdist (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *e, *s, *r = n->v.b1.b;
	int i, k, m = r->v.bn.n_nodes;
	int rgen = (n->t == F_RANDGEN);
	int mrgen = (n->t == F_MRANDGEN);
	int rgen1 = (n->t == F_RANDGEN1);
	double parm[3] = {0};
	double argval = NADBL;
	double *parmvec[2] = { NULL };
	double *argvec = NULL;
	gretl_matrix *argmat = NULL;
	int rows = 0, cols = 0;
	int d, np, argc, bb;

	if (mrgen) {
	    if (m < 4 || m > 7) {
		p->err = E_INVARG;
		goto disterr;
	    }
	} else if (m < 2 || m > 5) {
	    p->err = E_INVARG;
	    goto disterr;
	}

	s = r->v.bn.n[0];
	if (s->t == STR) {
	    d = dist_code_from_string(s->v.str);
	    if (d == 0) {
		p->err = E_INVARG;
		goto disterr;
	    }
	} else {
	    node_type_error(n->t, 0, STR, s, p);
	    goto disterr;
	}

	p->err = check_dist_count(d, n->t, &np, &argc);
	k = np + argc + 2 * mrgen;
	if (!p->err && k != m - 1) {
	    p->err = E_INVARG;
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
	} 

	for (i=1; i<=k && !p->err; i++) {
	    s = r->v.bn.n[i];
	    e = eval(s, p);
	    if (p->err) {
		break;
	    }	    
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
	    } else if (i == k && e->t == VEC) {
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
		if (rgen || mrgen) {
		    node_type_error(n->t, i, NUM, e, p);
		} else {
		    argmat = e->v.m;
		}
	    } else if (e->t == VEC) {
		/* a series param for randgen? */
		if (rgen && !bb) {
		    parmvec[i-1] = e->v.xvec;
		} else {
		    node_type_error(n->t, i, NUM, e, p);
		}
	    } else {
		p->err = E_INVARG;
		fprintf(stderr, "eval_pdist: arg %d, bad type %d\n", i+1, e->t);
	    }

	    if (!reusable(p)) { 
		free_tree(s, p, "Pdist");
		r->v.bn.n[i] = NULL;
	    }		    
	}

	if (p->err) {
	    goto disterr;
	}	
	
	if (mrgen) {
	    ret = aux_matrix_node(p);
	} else if (rgen || argvec != NULL) {
	    ret = aux_vec_node(p, 0);
	} else if (argmat != NULL) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = aux_scalar_node(p);
	}

	if (ret == NULL) {
	    goto disterr;
	}

	if (rgen) {
	    ret->v.xvec = gretl_get_random_series(d, parm, 
						  parmvec[0], parmvec[1], 
						  p->dset, &p->err);
	} else if (mrgen) {
	    ret->v.m = gretl_get_random_matrix(d, parm, rows, cols, 
					       &p->err);
	} else if (rgen1) {
	    ret->v.xval = gretl_get_random_scalar(d, parm, &p->err);
	} else if (argvec != NULL) {
	    ret->v.xvec = series_pdist(n->t, d, parm, np, argvec, p);
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
	/* IEEE 754 - 2008, double precision */
	return pow(2.0, -53);
    } else if (id == CONST_INF) {
#ifdef INFINITY
	return INFINITY;
#else
	return 1.0/0.0;
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
	    *err = E_UNKVAR;
	}
	return NADBL;
    }
}

#ifdef HAVE_MPI

static Gretl_MPI_Op reduce_op_from_string (const char *s)
{
    if (!strcmp(s, "sum")) {
	return GRETL_MPI_SUM;
    } else if (!strcmp(s, "prod")) {
	return GRETL_MPI_PROD;
    } else if (!strcmp(s, "max")) {
	return GRETL_MPI_MAX;
    } else if (!strcmp(s, "min")) {
	return GRETL_MPI_MIN;
    } else if (!strcmp(s, "hcat")) {
	return GRETL_MPI_HCAT;
    } else if (!strcmp(s, "vcat")) {
	return GRETL_MPI_VCAT;
    } else {
	return 0;
    }
}

static Gretl_MPI_Op scatter_op_from_string (const char *s)
{
    if (!strcmp(s, "bycols")) {
	return GRETL_MPI_HSPLIT;
    } else if (!strcmp(s, "byrows")) {
	return GRETL_MPI_VSPLIT;
    } else {
	return 0;
    }
}

#endif

static NODE *mpi_transfer_node (NODE *l, NODE *r, NODE *r2,
				int f, parser *p)
{
#ifndef HAVE_MPI
    gretl_errmsg_set(_("MPI is not supported in this gretl build"));
    p->err = 1;
    return NULL;
#else
    NODE *ret = NULL;
    GretlType type = 0;
    int root = 0;
    int id = 0;

    if (!gretl_mpi_initialized()) {
	gretl_errmsg_set(_("The MPI library is not loaded"));
	p->err = 1;
	return NULL;
    }

    if (f == F_MPI_SEND) {
	/* we support sending a matrix or scalar; we need
	   the destination id as second argument 
	*/
	if (l->t == MAT) {
	    type = GRETL_TYPE_MATRIX;
	} else if (l->t == NUM) {
	    type = GRETL_TYPE_DOUBLE;
	} else {
	    p->err = E_TYPES;
	}
	if (!p->err) {
	    /* destination id */
	    id = node_get_int(r, p);
	}
    } else if (f == F_MPI_RECV) {
	/* the single argument is the source id */
	id = node_get_int(l, p);
    } else if (f == F_BCAST || f == F_REDUCE || 
	       f == F_ALLREDUCE || f == F_SCATTER) {
	/* we need a variable's address on the left */
	if (l->t != U_ADDR) {
	    p->err = E_TYPES;
	} else {
	    /* switch to 'content' sub-node */
	    l = l->v.b1.b;
	    if (umatrix_node(l)) {
		type = GRETL_TYPE_MATRIX;
	    } else if (f != F_SCATTER && uscalar_node(l)) {
		type = GRETL_TYPE_DOUBLE;
	    } else {
		p->err = E_TYPES;
	    }
	}
	if (!p->err && f != F_ALLREDUCE) {
	    /* optional root specification */
	    NODE *rootspec = (f == F_BCAST)? r : r2;

	    if (!null_or_empty(rootspec)) {
		root = node_get_int(rootspec, p);
	    }
	}
	if (!p->err) {
	    /* "self" id */
	    id = gretl_mpi_rank();
	}	    
    }

    if (p->err) {
	return NULL;
    } else if (f == F_MPI_SEND) {
	void *sendp;

	if (type == GRETL_TYPE_MATRIX) {
	    sendp = l->v.m;
	} else {
	    sendp = &l->v.xval;
	}
	ret = aux_scalar_node(p);
	if (!p->err) {
	    p->err = ret->v.xval = gretl_mpi_send(sendp, type, id);
	}
    } else if (f == F_MPI_RECV) {
	gretl_matrix *m = NULL;
	double x = NADBL;

	p->err = gretl_mpi_receive(id, &type, &m, &x);

	if (!p->err) {
	    if (type == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (!p->err) {
		    ret->v.m = m;
		}
	    } else {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = x;
		}
	    }
	}
    } else if (f == F_BCAST) {
	void *bcastp;
	gretl_matrix *m = NULL;
	double x = NADBL;

	if (type == GRETL_TYPE_MATRIX) {
	    if (id == root) {
		m = l->v.m;
	    }
	    bcastp = &m;
	} else {
	    x = l->v.xval;
	    bcastp = &x;
	}	    

	ret = aux_scalar_node(p);
	if (!p->err) {
	    p->err = gretl_mpi_bcast(bcastp, type, root);
	    if (!p->err && id != root) {
		if (type == GRETL_TYPE_MATRIX) {
		    p->err = user_matrix_replace_matrix_by_name(l->vname, m);
		} else {
		    p->err = gretl_scalar_set_value(l->vname, x);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_REDUCE || f == F_ALLREDUCE) {
	ret = aux_scalar_node(p);
	if (!p->err) {
	    Gretl_MPI_Op op = reduce_op_from_string(r->v.str);
	    gretlopt opt = (f == F_REDUCE)? OPT_NONE : OPT_A;
	    gretl_matrix *m = NULL;
	    double x = NADBL;

	    if (type == GRETL_TYPE_MATRIX) { 
		p->err = gretl_matrix_mpi_reduce(l->v.m, &m, op, root, opt);
	    } else {
		p->err = gretl_scalar_mpi_reduce(l->v.xval, &x, op, root, opt);
	    }
	    if (!p->err && (id == root || f == F_ALLREDUCE)) {
		if (type == GRETL_TYPE_MATRIX) {
		    p->err = user_matrix_replace_matrix_by_name(l->vname, m);
		} else {
		    p->err = gretl_scalar_set_value(l->vname, x);
		}
	    }
	    ret->v.xval = p->err;
	}
    } else if (f == F_SCATTER) {
	ret = aux_scalar_node(p);
	if (!p->err) {
	    Gretl_MPI_Op op = scatter_op_from_string(r->v.str);
	    gretl_matrix *m = NULL;

	    p->err = ret->v.xval = gretl_matrix_mpi_scatter(l->v.m, &m, 
							    op, root);
	    if (!p->err) {
		p->err = user_matrix_replace_matrix_by_name(l->vname, m);
	    }
	}	
    } else {
	gretl_errmsg_set("MPI function not yet supported");
	p->err = 1;
    }

    return ret;
#endif /* HAVE_MPI */
}

static NODE *scalar_calc (NODE *x, NODE *y, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.xval = xy_calc(x->v.xval, y->v.xval, f, NUM, p);
    }

    return ret;
}

static NODE *string_offset (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int n = strlen(l->v.str);
	int k = r->v.xval;

	if (k < 0) {
	    p->err = E_DATA;
	} else if (k >= n) {
	    ret->v.str = gretl_strdup("");
	} else {
	    ret->v.str = gretl_strdup(l->v.str + k);
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
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

     series foo = x=="strval"

   It can also work if the string is an observation marker, as in

     series foo = obs>="CA"
*/

static NODE *series_string_calc (NODE *l, NODE *r, int f, parser *p)
{
    double xt = NADBL, yt = NADBL;
    double *x = NULL, *y = NULL;
    double *alt;
    const char *strval;
    int vnum, t, t1, t2;
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

    if (vnum > 0 && is_string_valued(p->dset, vnum)) {
	*alt = series_decode_string(p->dset, vnum, strval);
    }

    if (na(*alt)) {
	/* try for an observation string */
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

    ret = aux_vec_node(p, p->dset->n);
    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
	if (x != NULL) {
	    xt = x[t];
	} else if (y != NULL) {
	    yt = y[t];
	}
	ret->v.xvec[t] = xy_calc(xt, yt, f, VEC, p);
    }

    return ret;
}

/* At least one of the nodes is a series; the other may be a
   scalar or 1 x 1 matrix */

static NODE *series_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;
    const double *x = NULL, *y = NULL;
    double xt = 0, yt = 0;
    int t, t1, t2;

    ret = aux_vec_node(p, p->dset->n);
    if (ret == NULL) {
	return NULL;
    }

    if (l->t == VEC) {
	x = l->v.xvec;
    } else if (l->t == NUM) {
	xt = l->v.xval;
    } else if (l->t == MAT) {
	xt = l->v.m->val[0];
    }

    if (r->t == VEC) {
	y = r->v.xvec;
    } else if (r->t == NUM) {
	yt = r->v.xval;
    } else if (r->t == MAT) {
	yt = r->v.m->val[0];
    } 

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
	if (x != NULL) {
	    xt = x[t];
	}
	if (y != NULL) {
	    yt = y[t];
	}
	ret->v.xvec[t] = xy_calc(xt, yt, f, VEC, p);
    }

    return ret;
}

/* Both nodes are string-valued series */

static NODE *stringvec_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    const char *sl, *sr;
    int t, eq;

    if (f != B_EQ && f != B_NEQ) {
	p->err = E_TYPES;
	return NULL;
    }

    ret = aux_vec_node(p, p->dset->n);
    if (ret == NULL) {
	return NULL;
    }

    for (t=p->dset->t1; t<=p->dset->t2; t++) {
	sl = series_get_string_for_obs(p->dset, l->vnum, t);
	sr = series_get_string_for_obs(p->dset, r->vnum, t);
	if (sl == NULL || sr == NULL) {
	    ret->v.xvec[t] = NADBL;
	} else {
	    eq = strcmp(sl, sr) == 0;
	    ret->v.xvec[t] = (f == B_EQ)? eq : !eq;
	}
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

static gretl_matrix *
matrix_add_sub_scalar (const gretl_matrix *A,
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

/* return allocated result of binary operation performed on
   two matrices */

static gretl_matrix *real_matrix_calc (const gretl_matrix *A, 
				       const gretl_matrix *B, 
				       int op, int *err) 
{
    gretl_matrix *C = NULL;
    int ra, ca;
    int rb, cb;
    int r, c;

    if (gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(B)) {
	if (op != B_HCAT && op != B_VCAT && op != F_DSUM &&
	    op != B_MUL && op != B_TRMUL) {
	    *err = E_NONCONF;
	    return NULL;
	}
	if (op == B_MUL || op == B_TRMUL) {
	    return nullmat_multiply(A, B, op, err);
	}
    }

    switch (op) {
    case B_ADD:
    case B_SUB:
	if (gretl_matrix_is_scalar(A) ||
	    gretl_matrix_is_scalar(B)) {
	    C = matrix_add_sub_scalar(A, B, op);
	    if (C == NULL) {
		*err = E_ALLOC;
	    }
	} else {
	    C = gretl_matrix_alloc(A->rows, A->cols);
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else if (op == B_ADD) {
		*err = gretl_matrix_add(A, B, C);
	    } else {
		*err = gretl_matrix_subtract(A, B, C);
	    }
	}
	break;
    case B_HCAT:
	C = gretl_matrix_col_concat(A, B, err);
	break;
    case B_VCAT:
	C = gretl_matrix_row_concat(A, B, err);
	break;
    case F_DSUM:
	C = gretl_matrix_direct_sum(A, B, err);
	break;
    case B_MUL:
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	r = (ra == 1 && ca == 1)? rb : ra;
	c = (rb == 1 && cb == 1)? ca : cb;

	C = gretl_matrix_alloc(r, c);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = gretl_matrix_multiply(A, B, C);
	    if (!*err) {
		gretl_matrix_transcribe_obs_info(C, A);
	    }
	}	
	break;
    case B_TRMUL:
	ra = gretl_matrix_cols(A);
	ca = gretl_matrix_rows(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	r = (ra == 1 && ca == 1)? rb : ra;
	c = (rb == 1 && cb == 1)? ca : cb;

	C = gretl_matrix_alloc(r, c);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else {
	    *err = gretl_matrix_multiply_mod(A, GRETL_MOD_TRANSPOSE,
					     B, GRETL_MOD_NONE,
					     C, GRETL_MOD_NONE);
	}	
	break;
    case F_QFORM:
	/* quadratic form, A * B * A', for symmetric B */
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	if (ca != rb || cb != rb) {
	    *err = E_NONCONF;
	} else if (!gretl_matrix_is_symmetric(B)) {
	    *err = E_NONCONF;
	} else {
	    C = gretl_matrix_alloc(ra, ra);
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else {
		*err = gretl_matrix_qform(A, GRETL_MOD_NONE, B,
					  C, GRETL_MOD_NONE);
	    }
	}
	break;
    case B_DIV:
    case B_LDIV:
	/* matrix right or left "division" */
	if (op == B_LDIV) {
	    C = gretl_matrix_divide(A, B, GRETL_MOD_NONE, err);
	} else {
	    /* A/B = (B'\A')' */
	    C = gretl_matrix_divide(A, B, GRETL_MOD_TRANSPOSE, err);
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
	/* apply operator element-wise */
	C = gretl_matrix_dot_op(A, B, op_symbol(op), err);
	break;
    case B_KRON:
	/* Kronecker product */
	C = gretl_matrix_kronecker_product_new(A, B, err);
	break;
    case F_HDPROD:
	C = gretl_matrix_hdproduct_new(A, B, err);
	break;    
    case F_CMULT:
	C = gretl_matrix_complex_multiply(A, B, err);
	break;
    case F_CDIV:
	C = gretl_matrix_complex_divide(A, B, err);
	break;
    case F_MRSEL:
	C = gretl_matrix_bool_sel(A, B, 1, err);
	break;
    case F_MCSEL:
	C = gretl_matrix_bool_sel(A, B, 0, err);
	break;
    default:
	*err = E_TYPES;
	break;
    } 

    if (!*err) {
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

    if (*err && C != NULL) {
	gretl_matrix_free(C);
	C = NULL;
    }

    return C;
}

static gretl_matrix *tmp_matrix_from_series (NODE *n, parser *p)
{
    int t, T = sample_size(p->dset);
    const double *x = n->v.xvec;
    gretl_matrix *m = NULL;

    m = gretl_column_vector_alloc(T);

    if (m == NULL) {
	p->err = E_ALLOC;
    } else {
	int i = 0;

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    if (na(x[t])) {
		m->val[i++] = M_NA;
	    } else {
		m->val[i++] = x[t];
	    }
	}
    }

    return m;
}

const double *get_colvec_as_series (NODE *n, int f, parser *p)
{
    if (n->t != MAT) {
	node_type_error(f, 1, VEC, n, p);
	return NULL;
    } else {	
	const gretl_matrix *m = n->v.m;

	if (m->rows == p->dset->n && m->cols == 1) {
	    return m->val;
	} else {
	    node_type_error(f, 1, VEC, n, p);
	    return NULL;
	}
    } 
}

/* One of the operands is a matrix, the other a series: we
   try "casting" the series to a matrix.
*/

static NODE *matrix_series_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *a = NULL;
	gretl_matrix *b = NULL;
	gretl_matrix *c = NULL;

	if (l->t == VEC) {
	    a = tmp_matrix_from_series(l, p);
	    c = a;
	    b = r->v.m;
	} else {
	    a = l->v.m;
	    b = tmp_matrix_from_series(r, p);
	    c = b;
	}

	if (!p->err) {
	    ret->v.m = real_matrix_calc(a, b, op, &p->err);
	}

	gretl_matrix_free(c);
    }

    return ret;
}

static int 
matrix_pow_check (int t, double x, const gretl_matrix *m, parser *p)
{
    if (t != MAT) {
	p->err = E_TYPES;
    } else if (gretl_is_null_matrix(m)) {
	p->err = E_DATA;
    } else if (m->rows != m->cols) {
	p->err = E_NONCONF;
    } else if (x < 0 || x > (double) INT_MAX || floor(x) != x) {
	p->err = E_DATA;
    } 

    return p->err;
}

#define comparison_op(o) (o == B_EQ  || o == B_NEQ || \
			  o == B_LT  || o == B_GT ||  \
			  o == B_LTE || o == B_GTE)

/* one of the operands is a matrix, the other a scalar, giving a
   matrix result unless we're looking at a comparison operator.
*/

static NODE *matrix_scalar_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = NULL;

    if (op == B_KRON) {
	p->err = E_TYPES;
	return NULL;
    }

    if (starting(p)) {
	const gretl_matrix *m = NULL;
	int comp = comparison_op(op);
	double y, x = 0.0;
	int i, n = 0;

	x = (l->t == NUM)? l->v.xval : r->v.xval;
	m = (l->t == MAT)? l->v.m : r->v.m;

	if (gretl_is_null_matrix(m)) {
	    p->err = E_DATA;
	    return NULL;
	}

	n = m->rows * m->cols;

	/* special: raising a matrix to an integer power */
	if (n > 1 && op == B_POW && matrix_pow_check(l->t, x, m, p)) {
	    fprintf(stderr, "matrix_pow_check failed\n");
	    return NULL;
	}

	/* mod: scalar must be on the right */
	if (op == B_MOD && l->t == NUM) {
	    p->err = E_TYPES;
	    return NULL;
	}

	if (comp || (op == B_POW && n == 1)) {
	    ret = aux_scalar_node(p);
	} else {
	    ret = aux_matrix_node(p);
	}

	if (ret == NULL) { 
	    return NULL;
	}

	if (op == B_POW) {
	    if (n > 1) {
		ret->v.m = gretl_matrix_pow(m, (int) x, &p->err);
	    } else if (l->t == NUM) {
		ret->v.xval = xy_calc(x, m->val[0], op, MAT, p);
	    } else {
		ret->v.xval = xy_calc(m->val[0], x, op, MAT, p);
	    }
	    return ret;
	} else if (op == B_TRMUL) {
	    gretl_matrix *tmp;

	    if (l->t == NUM) {
		tmp = gretl_matrix_copy(m);
	    } else {
		tmp = gretl_matrix_copy_transpose(m);
	    }
	    if (tmp == NULL) {
		p->err = E_ALLOC;
	    } else {
		gretl_matrix_multiply_by_scalar(tmp, x);
		ret->v.m = tmp;
	    }
	    return ret;
	}

	if (comp) {
	    ret->v.xval = 1;
	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    if (xy_calc(x, m->val[i], op, MAT, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}
	    } else {
		for (i=0; i<n; i++) {
		    if (xy_calc(m->val[i], x, op, MAT, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}		
	    }
	} else {
	    if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
		return NULL;
	    }

	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    y = xy_calc(x, m->val[i], op, MAT, p);
		    ret->v.m->val[i] = y;
		}
	    } else {
		for (i=0; i<n; i++) {
		    y = xy_calc(m->val[i], x, op, MAT, p);
		    ret->v.m->val[i] = y;
		}	
	    }
	} 
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *matrix_transpose_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (is_tmp_node(n)) {
	    /* transpose temp matrix in place */
	    p->err = gretl_matrix_transpose_in_place(n->v.m);
	    ret = n;
	} else {
	    /* create transpose as new matrix */
	    ret = aux_matrix_node(p);
	    if (!p->err) {
		ret->v.m = gretl_matrix_copy_transpose(n->v.m);
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

static NODE *scalar_transpose_node (NODE *n, parser *p)
{
    NODE *ret = n;

    if (reusable(p)) {
	/* we must duplicate the node to avoid double-freeing */
	ret = aux_scalar_node(p);
	if (starting(p) && !p->err) {
	    ret->v.xval = n->v.xval;
	}
    }

    return ret;
}

/* We're looking at a string argument that is supposed to represent
   a function call: we'll do a rudimentary heuristic check here.
   FIXME this should be more rigorous.
*/

static int is_function_call (const char *s)
{
    if (!strchr(s, '(')) {
	return 0;
    } else {
	return 1;
    }
}

static NODE *numeric_jacobian (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = r->v.str;

	if (!is_function_call(s)) {
	    p->err = E_TYPES;
	    return NULL;
	}

	ret = aux_matrix_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.m = fdjac(l->v.m, s, p->dset, &p->err);
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static NODE *BFGS_maximize (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *b = l->v.m;
	const char *sf = m->v.str;
	const char *sg = NULL;

	if (r->t == STR) {
	    sg = r->v.str;
	} else if (r->t != EMPTY) {
	    p->err = E_TYPES;
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
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = user_BFGS(b, sf, sg, p->dset, 
				p->prn, &p->err);
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *simann_node (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *b = l->v.m;
	const char *sf = m->v.str;
	int maxit = 0;

	if (gretl_is_null_matrix(b)) {
	    p->err = E_DATA;
	} else if (!is_function_call(sf)) {
	    p->err = E_TYPES;
	} 

	if (!p->err) {
	    if (scalar_node(r)) {
		maxit = node_get_int(r, p);
	    } else if (!null_or_empty(r)) {
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    return NULL;
	}

	ret = aux_scalar_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = user_simann(b, sf, maxit, p->dset, 
				  p->prn, &p->err);
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
	if (dated_daily_data(p->dset)) {
	    if (s >= 0 && s < p->dset->n) {
		while (s >= 0 && xna(x[s])) {
		    s--;
		}
	    }
	} else if (p->dset->structure == STACKED_TIME_SERIES) {
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

static NODE *matrix_text_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = m->v.str;
	int export = 0;
	
	if (!null_or_empty(r)) {
	    export = (r->v.xval != 0);
	}

	ret = aux_scalar_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = gretl_matrix_write_as_text(l->v.m, s, export);
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *bundle_text_write (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = m->v.str;
	int export = 0;
	
	if (!null_or_empty(r)) {
	    export = (r->v.xval != 0);
	}

	ret = aux_scalar_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = gretl_bundle_write_as_xml(l->v.b, s, export);
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

/* matrix on left, scalar on right */

static NODE *matrix_scalar_func (NODE *l, NODE *r, 
				 int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *m = l->v.m;
	int k = node_get_int(r, p);

	if (gretl_is_null_matrix(m)) {
	    p->err = E_INVARG;
	}

	if (p->err) {
	    return NULL;
	}

	ret = aux_matrix_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	if (f == F_MSORTBY) {
	    ret->v.m = gretl_matrix_sort_by_column(m, k-1, &p->err);
	} 
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static gretl_matrix *make_scalar_matrix (double x)
{
    gretl_matrix *m = gretl_matrix_alloc(1, 1);

    if (m != NULL) {
	m->val[0] = x;
    }

    return m;
}

/* both operands are known to be matrices or scalars */

static NODE *matrix_matrix_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);
    gretl_matrix *ml, *mr;

    if (op == B_POW) {
	p->err = E_TYPES;
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "matrix_matrix_calc: l=%p, r=%p, ret=%p\n",
	    (void *) l, (void *) r, (void *) ret);
#endif

    if (l->t == NUM) {
	ml = make_scalar_matrix(l->v.xval);
    } else {
	ml = l->v.m;
    }

    if (r->t == NUM) {
	mr = make_scalar_matrix(r->v.xval);
    } else {
	mr = r->v.m;
    }

    if (ret != NULL && starting(p)) {
	ret->v.m = real_matrix_calc(ml, mr, op, &p->err);
    }

    if (l->t == NUM) gretl_matrix_free(ml);
    if (r->t == NUM) gretl_matrix_free(mr);

    return ret;
}

static NODE *matrix_bundle_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);
    gretl_matrix *ml, *mr;
    gretl_matrix *mtmp = NULL;

    if (ret == NULL) {
	return NULL;
    }

    if (l->t == BUNDLE) {
	mtmp = ml = gretl_bundle_get_payload_matrix(l->v.b, &p->err);
	mr = r->v.m;
    } else {
	ml = l->v.m;
	mtmp = mr = gretl_bundle_get_payload_matrix(r->v.b, &p->err);
    }

    if (!p->err) {
	ret->v.m = real_matrix_calc(ml, mr, op, &p->err);
    }

    gretl_matrix_free(mtmp);

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
    NODE *ret = aux_scalar_node(p);

    if (op == B_OR || op == B_AND) {
	return matrix_and_or(l, r, op, p);
    }

    if (ret != NULL && starting(p)) {
	const gretl_matrix *a = l->v.m;
	const gretl_matrix *b = r->v.m;
	int i, n = a->rows * a->cols;

	if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
	    ret->v.xval = NADBL;
	} else if (a->rows != b->rows || a->cols != b->cols) {
	    ret->v.xval = NADBL;
	} else {
	    ret->v.xval = 1;
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
		} else if (op == B_NEQ && a->val[i] == b->val[i]) {
		    ret->v.xval = 0;
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

static NODE *matrix_to_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m;
	int tmpmat = 0;

	if (n->t == NUM) {
	    m = gretl_matrix_from_scalar(node_get_scalar(n, p));
	    tmpmat = 1;
	} else {
	    tmpmat = is_tmp_node(n);
	    m = n->v.m;
	}

	switch (f) {
	case F_ROWS:
	    ret->v.xval = m->rows;
	    break;
	case F_COLS:
	    ret->v.xval = m->cols;
	    break;
	case F_DET:
	case F_LDET:
	    ret->v.xval = user_matrix_get_determinant(m, tmpmat, f, &p->err);
	    break;
	case F_TRACE:
	    ret->v.xval = gretl_matrix_trace(m);
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
	case F_RANK:
	    ret->v.xval = gretl_matrix_rank(m, &p->err);
	    break;
	default:
	    p->err = E_PARSE;
	    break;
	}

	if (n->t == NUM) {
	    gretl_matrix_free(m);
	}

	if (p->err) {
	    matrix_error(p);
	}    
    }

    return ret;
}

static NODE *matrix_add_names (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m = l->v.m;
	int byrow = (f == F_ROWNAMES);

	if (r->t == STR) {
	    ret->v.xval = umatrix_set_names_from_string(m, r->v.str, byrow);
	} else {
	    int *list = node_get_list(r, p);

	    if (p->err) {
		ret->v.xval = 1;
	    } else {
		ret->v.xval = umatrix_set_names_from_list(m, list, p->dset,
							  byrow);
	    }
	    free(list);
	}
    }

    return ret;
}

static NODE *matrix_get_colname (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int c = node_get_int(r, p);

	ret->v.str = user_matrix_get_column_name(l->v.m, c, &p->err);
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

static void matrix_minmax_indices (int f, int *mm, int *rc, int *idx)
{
    *mm = (f == F_MAXR || f == F_MAXC || f == F_IMAXR || f == F_IMAXC);
    *rc = (f == F_MINC || f == F_MAXC || f == F_IMINC || f == F_IMAXC);
    *idx = (f == F_IMINR || f == F_IMINC || f == F_IMAXR || f == F_IMAXC);
}

static NODE *matrix_to_matrix_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m = NULL;
	int tmpmat = 0;
	int a, b, c;

	if (n->t == NUM) {
	    m = gretl_matrix_from_scalar(node_get_scalar(n, p));
	    tmpmat = 1;
	} else {
	    tmpmat = is_tmp_node(n);
	    m = n->v.m;
	}

	if (gretl_is_null_matrix(m) && !nullmat_ok(f)) {
	    p->err = E_DATA;
	} else if (f == F_RESAMPLE || f == F_MREVERSE || f == F_SDC) {
	    /* the r node may be absent, but if present it should
	       hold a scalar */
	    if (!empty_or_num(r)) {
		node_type_error(f, 2, NUM, r, p);
	    }
	} else if (f == F_RANKING) {
	    if (gretl_vector_get_length(m) == 0) {
		/* m must be a vector */
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    goto finalize;
	}

	gretl_error_clear();

	switch (f) {
	case F_SUMC:
	    ret->v.m = gretl_matrix_column_sum(m, &p->err);
	    break;
	case F_SUMR:
	    ret->v.m = gretl_matrix_row_sum(m, &p->err);
	    break;
	case F_PRODC:
	    ret->v.m = gretl_matrix_column_prod(m, &p->err);
	    break;
	case F_PRODR:
	    ret->v.m = gretl_matrix_row_prod(m, &p->err);
	    break;
	case F_MEANC:
	    ret->v.m = gretl_matrix_column_mean(m, &p->err);
	    break;
	case F_MEANR:
	    ret->v.m = gretl_matrix_row_mean(m, &p->err);
	    break;
	case F_SD:
	    ret->v.m = gretl_matrix_column_sd(m, &p->err);
	    break;
	case F_SDC:
	    if (r != NULL && r->t == NUM) {
		ret->v.m = gretl_matrix_column_sd2(m, r->v.xval, &p->err);
	    } else {
		ret->v.m = gretl_matrix_column_sd(m, &p->err);
	    }
	    break;
	case F_MCOV:
	    ret->v.m = gretl_covariance_matrix(m, 0, &p->err);
	    break;
	case F_MCORR:
	    ret->v.m = gretl_covariance_matrix(m, 1, &p->err);
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
	case F_RESAMPLE:
	    if (r != NULL && r->t == NUM) {
		ret->v.m = gretl_matrix_block_resample(m, r->v.xval, &p->err);
	    } else {
		ret->v.m = gretl_matrix_resample(m, &p->err);
	    }
	    break;
	case F_CDEMEAN:
	case F_CHOL:
	case F_PSDROOT:
	case F_INV:
	case F_INVPD:
	case F_GINV:
	case F_UPPER:
	case F_LOWER:
	    ret->v.m = user_matrix_matrix_func(m, tmpmat, f, &p->err);
	    break;
	case F_DIAG:
	    ret->v.m = gretl_matrix_get_diagonal(m, &p->err);
	    break;
	case F_TRANSP:
	    ret->v.m = gretl_matrix_copy_transpose(m);
	    break;
	case F_VEC:
	    ret->v.m = user_matrix_vec(m, &p->err);
	    break;
	case F_VECH:
	    ret->v.m = user_matrix_vech(m, &p->err);
	    break;
	case F_UNVECH:
	    ret->v.m = user_matrix_unvech(m, &p->err);
	    break;
	case F_MREVERSE:
	    if (r != NULL && r->t == NUM && r->v.xval != 0) {
		ret->v.m = gretl_matrix_reverse_cols(m);
	    } else {
		ret->v.m = gretl_matrix_reverse_rows(m);
	    }
	    break;
	case F_NULLSPC:
	    ret->v.m = gretl_matrix_right_nullspace(m, &p->err);
	    break;
	case F_MEXP:
	    ret->v.m = gretl_matrix_exp(m, &p->err);
	    break;
	case F_FFT:
	    ret->v.m = gretl_matrix_fft(m, &p->err);
	    break;
	case F_FFTI:
	    ret->v.m = gretl_matrix_ffti(m, &p->err);
	    break;
	case F_POLROOTS:
	    ret->v.m = gretl_matrix_polroots(m, &p->err);
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
	    ret->v.m = gretl_matrix_minmax(m, a, b, c, &p->err);
	    break;
	default:
	    break;
	}

	if (ret->v.m == m) {
	    /* input matrix was recycled: avoid double-freeing */
	    if (n->t == MAT) {
		n->v.m = NULL;
	    } else {
		m = NULL;
	    }
	}	

    finalize:

	if (n->t == NUM) {
	    /* trash the on-the-fly scalar matrix */
	    gretl_matrix_free(m);
	}

	if (ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
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
	int import = 0;

	if (!null_or_empty(r)) {
	    import = (r->v.xval != 0);
	}

	gretl_error_clear();

	switch (f) {
	case F_MREAD:
	    ret->v.m = gretl_matrix_read_from_text(n->v.str, import, 
						   &p->err);
	    break;
	case F_BREAD:
	    ret->v.b = gretl_bundle_read_from_xml(n->v.str, import, 
						  &p->err);
	    break;
	default:
	    break;
	}

	if (f == F_MREAD && p->err == E_FOPEN) {
	    /* should we be doing this?? */
	    ret->v.m = gretl_null_matrix_new();
	    p->err = 0;
	}

	if (f == F_MREAD && ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
}

static NODE *
matrix_to_matrix2_func (NODE *n, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;
	const char *rname;

	if (gretl_is_null_matrix(m)) {
	    p->err = E_DATA;
	    goto finalize;
	}

	gretl_error_clear();

	/* on the right: address of matrix or null */
	if (r->t == EMPTY) {
	    rname = "null";
	} else {
	    /* note: switch to the 'content' sub-node */
	    r = r->v.b1.b;
	    if (umatrix_node(r)) {
		rname = r->vname;
	    } else {
		p->err = E_PARSE;
		gretl_errmsg_set("Expected the address of a matrix");
		return ret;
	    }
	}

	switch (f) {
	case F_QR:
	    ret->v.m = user_matrix_QR_decomp(m, rname, &p->err);
	    break;
	case F_EIGSYM:
	    ret->v.m = user_matrix_eigen_analysis(m, rname, 1, &p->err);
	    break;
	case F_EIGGEN:
	    ret->v.m = user_matrix_eigen_analysis(m, rname, 0, &p->err);
	    break;
	}

    finalize:

	if (ret->v.m == NULL) {
	    matrix_error(p);
	}
    }

    return ret;
}

static int ok_matrix_dim (int r, int c, int f)
{
    if (f == F_IMAT || f == F_ZEROS || f == F_ONES || f == F_MUNIF || \
	f == F_MNORM) {
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
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	int cols = 0, rows = node_get_int(l, p);

	if (!p->err) {
	    if (f == F_IMAT) {
		cols = rows;
	    } else {
		cols = node_get_int(r, p);
	    }
	}

	if (!p->err && !ok_matrix_dim(rows, cols, f)) {
	    p->err = E_INVARG;
	    matrix_error(p);
	}

	if (p->err) {
	    return ret;
	}

	switch (f) {
	case F_IMAT:
	    ret->v.m = gretl_identity_matrix_new(rows);
	    break;
	case F_ZEROS:
	    ret->v.m = gretl_zero_matrix_new(rows, cols);
	    break;
	case F_ONES:
	    ret->v.m = gretl_unit_matrix_new(rows, cols);
	    break;
	case F_MUNIF:
	    ret->v.m = gretl_random_matrix_new(rows, cols, 
					       D_UNIFORM);
	    break;
	case F_MNORM:
	    ret->v.m = gretl_random_matrix_new(rows, cols,
					       D_NORMAL);
	    break;
	default:
	    break;
	}
    }

    return ret;
}

#if LHDEBUG > 1

static void print_mspec (matrix_subspec *mspec)
{
    int i;

    fprintf(stderr, "mspec at %p:\n", (void *) mspec);

    if (mspec == NULL) {
	return;
    }

    for (i=0; i<2; i++) {
	fprintf(stderr, "type[%d] = %d\n", i, mspec->type[i]);
	if (mspec->type[i] == SEL_RANGE) {
	    fprintf(stderr, "sel[%d].range[0] = %d\n", i, mspec->sel[i].range[0]);
	    fprintf(stderr, "sel[%d].range[1] = %d\n", i, mspec->sel[i].range[1]);
	} else if (mspec->type[i] == SEL_MATRIX) {
	    gretl_matrix_print(mspec->sel[i].m, "sel matrix");
	}
    }
}

#endif

/* compose a sub-matrix specification, from scalars and/or
   index matrices */

static matrix_subspec *build_mspec (NODE *l, NODE *r, int *err)
{
    matrix_subspec *spec;

    spec = malloc(sizeof *spec);
    if (spec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    spec->rslice = spec->cslice = NULL;

    if (l->t == DUM) {
	if (l->v.idnum == DUM_DIAG) {
	    spec->type[0] = SEL_DIAG;
	    spec->type[1] = SEL_ALL;
	    return spec;
	} else {
	    *err = E_TYPES;
	    goto bailout;
	}
    } else if (l->t == NUM && r != NULL && r->t == NUM) {
	spec->type[0] = spec->type[1] = SEL_ELEMENT;
	mspec_set_row_index(spec, l->v.xval);
	mspec_set_col_index(spec, r->v.xval);
	return spec;
    }

    if (l->t == NUM) {
	spec->type[0] = SEL_RANGE;
	mspec_set_row_index(spec, l->v.xval);
    } else if (l->t == IVEC) {
	spec->type[0] = SEL_RANGE;
	spec->sel[0].range[0] = l->v.ivec[0];
	spec->sel[0].range[1] = l->v.ivec[1];
    } else if (l->t == MAT) {
	spec->type[0] = SEL_MATRIX;
	spec->sel[0].m = l->v.m;
    } else if (l->t == EMPTY) {
	spec->type[0] = SEL_ALL;
    } else {
	*err = E_TYPES;
	goto bailout;
    }

    if (r->t == ABSENT) {
	spec->type[1] = SEL_NULL;
    } else if (r->t == NUM) {
	spec->type[1] = SEL_RANGE;
	mspec_set_col_index(spec, r->v.xval);
    } else if (r->t == IVEC) {
	spec->type[1] = SEL_RANGE;
	spec->sel[1].range[0] = r->v.ivec[0];
	spec->sel[1].range[1] = r->v.ivec[1];
    } else if (r->t == MAT) {
	spec->type[1] = SEL_MATRIX;
	spec->sel[1].m = r->v.m;
    } else if (r->t == EMPTY) {
	spec->type[1] = SEL_ALL;
    } else {
	*err = E_TYPES;
	goto bailout;
    }

 bailout:
    
    if (*err) {
	free(spec);
	spec = NULL;
    }

    return spec;
}

/* node holding evaluated result of matrix specification */

static NODE *mspec_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_mspec_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.mspec = build_mspec(l, r, &p->err);
    }

    return ret;
}

static NODE *submatrix_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *a = NULL;

	if (r->t != MSPEC) {
	    fprintf(stderr, "submatrix_node: couldn't find mspec\n");
	    p->err = E_TYPES;
	    return NULL;
	}

	if (l->t == MAT) {
	    a = matrix_get_submatrix(l->v.m, r->v.mspec, 0, &p->err);
	} else if (l->t == STR) {
	    a = user_matrix_get_submatrix(l->v.str, r->v.mspec, &p->err);
	} else {
	    p->err = E_TYPES;
	}

	if (a != NULL) {
	    ret = aux_matrix_node(p);
	    if (ret == NULL) {
		gretl_matrix_free(a);
	    } else {
		ret->v.m = a;
	    }
	}
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static NODE *process_subslice (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (scalar_node(l) && (scalar_node(r) || r->t == EMPTY)) {
	    ret = aux_ivec_node(p, 2);
	    if (ret != NULL) {
		ret->v.ivec[0] = node_get_int(l, p);
		ret->v.ivec[1] = (r->t == EMPTY)? 
		    MSEL_MAX : node_get_int(r, p);
	    }
	} else {
	    p->err = E_TYPES;
	}
    } else {
	ret = aux_ivec_node(p, 2);
    }

    return ret;
}

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
	    if (na(x)) {
		return NADBL;
	    }
	}
    } 

    switch (f) {
    case U_NEG: 
	return -x;
    case U_POS: 
	return x;
    case U_NOT:
	return x == 0;
    case F_ABS:
	return fabs(x);
    case F_TOINT:
	return (double) (int) x;
    case F_CEIL:
	return ceil(x);
    case F_FLOOR:
	return floor(x);
    case F_ROUND:
	return gretl_round(x);
    case F_SIN:
	return sin(x);
    case F_COS:
	return cos(x);
    case F_TAN:
	return tan(x);
    case F_ASIN:
	return asin(x);
    case F_ACOS:
	return acos(x);
    case F_ATAN:
	return atan(x);
    case F_SINH:
	return sinh(x);
    case F_COSH:
	return cosh(x);
    case F_TANH:
	return tanh(x);
    case F_ASINH:
	return asinh(x);
    case F_ACOSH:
	return acosh(x);
    case F_ATANH:
	return atanh(x);
    case F_CNORM:
	return normal_cdf(x);
    case F_DNORM:
	return normal_pdf(x);
    case F_QNORM:
	return normal_cdf_inverse(x);
    case F_LOGISTIC:
	return logistic_cdf(x);
    case F_GAMMA:
	y = gamma_function(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}	
	return y;
    case F_LNGAMMA:
	y = ln_gamma(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_DIGAMMA:
	y = digamma(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_MISSING:
	return 0.0;
    case F_DATAOK:
	return 1.0;
    case F_MISSZERO:
	return x;
    case F_ZEROMISS:
	return (x == 0.0)? NADBL : x;
    case F_SQRT:
	y = sqrt(x);
	if (errno) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_LOG:
    case F_LOG10:
    case F_LOG2:
	y = log(x);
	if (errno) {
	    eval_warning(p, F_LOG, errno);
	} else {	    
	    if (f == F_LOG10) {
		y /= log(10.0);
	    } else if (f == F_LOG2) {
		y /= log(2.0);
	    }
	}
	return y;
    case F_EXP:
	y = exp(x);
	if (errno) {
	    eval_warning(p, F_EXP, errno);
	}
	return y;
    case F_INVMILLS:
	y = invmills(x);
	if (na(y)) {
	    eval_warning(p, f, errno);
	}
	return y;
    case F_EASTER:
	y = easterdate(x);
	return y;
    default:
	return 0.0;
    }
}

static NODE *apply_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	ret->v.xval = real_apply_func(n->v.xval, f, p);
    }

    return ret;
}

static NODE *scalar_isnan_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	double x = node_get_scalar(n, p);

	if (!p->err) {
	    ret->v.xval = isnan(x);
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
		    ret->v.m->val[i] = isnan(m->val[i]);
		}
	    }
	}
    }

    return ret;
}

static NODE *apply_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);
    int t, t1, t2;

    if (ret != NULL) {
	/* AC: changed for autoreg case, 2007/7/1 */
	const double *x;

	if (n->t == VEC) {
	    x = n->v.xvec;
	} else {
	    x = get_colvec_as_series(n, f, p);
	}

	if (!p->err) {
	    t1 = (autoreg(p))? p->obs : p->dset->t1;
	    t2 = (autoreg(p))? p->obs : p->dset->t2;
	    for (t=t1; t<=t2; t++) {
		ret->v.xvec[t] = real_apply_func(x[t], f, p);
	    }
	}
    }

    return ret;
}

/* argument is series; value returned is list */

static NODE *dummify_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	int *list = NULL;
	double oddval = NADBL;

	if (r->t != EMPTY) {
	    if (r->t != NUM) {
		p->err = E_TYPES;
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
	    p->err = E_TYPES;
	}

	if (p->err) {
	    ; /* don't do anything more */
	} else if (list == NULL) {
	    p->err = E_ALLOC;
	} else if (r->t == EMPTY) {
	    /* got just one argument */
	    p->err = list_dumgenr(&list, p->dset, OPT_F);
	    ret->v.ivec = list;
	} else if (list[0] > 1) {
	    gretl_errmsg_set("dummify(x, y): first argument should be a single variable");
	    free(list);
	    p->err = E_DATA;
	} else {
	    p->err = dumgenr_with_oddval(&list, p->dset, oddval);
	    ret->v.ivec = list;
	}
    }

    return ret;
}

/* argument is series or list; value returned is list in either
   case */

static NODE *list_make_lags (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (ret != NULL && starting(p)) {
	gretlopt opt = OPT_NONE;
	int k, *list = NULL;

	if (!null_or_empty(r)) {
	    k = node_get_int(r, p);
	    if (!p->err && k != 0) {
		opt = OPT_L;
	    }
	}

	if (!p->err) {
	    if (m->t == LIST) {
		list = gretl_list_copy(m->v.ivec);
	    } else if (useries_node(m)) {
		list = gretl_list_new(1);
		list[1] = m->vnum;
	    } else {
		p->err = E_TYPES;
	    }
	}

	if (!p->err) {
	    if (list == NULL) {
		p->err = E_ALLOC;
	    } else {
		k = node_get_int(l, p);
		if (!p->err) {
		    p->err = list_laggenr(&list, k, p->dset, opt);
		}
		ret->v.ivec = list;
	    }
	}
    }

    return ret;
}

#define ok_list_func(f) (f == F_LOG || f == F_DIFF || \
			 f == F_LDIFF || f == F_SDIFF || \
			 f == F_XPX || f == F_ODEV)

/* functions that are "basically" for series, but which
   can also be applied to lists */

static NODE *apply_list_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_list_node(p);

    if (!ok_list_func(f)) {
	p->err = E_TYPES;
	return ret;
    }

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);
	int t = 0;

	/* note: list is modified below */

	if (list != NULL) {
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
	    case F_XPX:
		p->err = list_xpxgenr(&list, p->dset, OPT_O);
		break;
	    case F_ODEV:
		p->err = list_orthdev(list, p->dset);
		break;
	    default:
		break;
	    }
	    ret->v.ivec = list;
	}
    }

    return ret;
}

static NODE *dataset_list_node (parser *p)
{
    NODE *ret = aux_list_node(p);

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
	if (p->dset == NULL || p->dset->n == 0) {
	    no_data_error(p);
	} else {
	    ret = aux_series_node(p, 0);
	}

	if (ret != NULL) {
	    p->err = gen_time(p->dset, 1);
	    if (!p->err) {
		ret->vname = gretl_strdup("time");
		ret->vnum = series_index(p->dset, "time");
		ret->v.xvec = p->dset->Z[ret->vnum];
		ret->flags = 0;
	    }
	}
    }

    return ret;
}

static NODE *get_lag_list (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	int *list = NULL;
	int lv;

	if (l->t != VEC || l->vnum < 0 || 
	    (r->t != IVEC && r->t != NUM)) {
	    p->err = E_TYPES;
	    return NULL;
	}

	lv = l->vnum;

	if (r->t == IVEC) {
	    int minlag = -r->v.ivec[0];
	    int maxlag = -r->v.ivec[1];

	    list = laggenr_from_to(lv, minlag, maxlag,
				   p->dset, &p->err);
	} else {
	    int lag = -r->v.xval;

	    lv = laggenr(lv, lag, p->dset);
	    if (lv > 0) {
		list = gretl_list_new(1);
		if (list != NULL) {
		    list[1] = lv;
		}
	    }
	}

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

/* get an *int list from node @n: note that the list is always
   newly allocated, and so should be freed by the caller if
   it's just for temporary use
*/

static int *node_get_list (NODE *n, parser *p)
{
    int *list = NULL;
    int v = 0;

    if (n->t == LIST) {
	list = gretl_list_copy(n->v.ivec);
    } else if (n->t == WLIST) {
	/* handle wildcard */
	list = varname_match_list(p->dset, n->v.str, &p->err);
    } else if (n->t == VEC || n->t == NUM) {
	v = (n->t == VEC)? n->vnum : n->v.xval;
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
    } else if (n->t == EMPTY) {
	list = gretl_null_list();
    } else if (dataset_dum(n)) {
	list = full_var_list(p->dset, NULL);
    } else if (n->t == MAT) {
	list = list_from_matrix(n->v.m, p->dset, &p->err);
    }

    if (p->err == E_UNKVAR && v != 0) {
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
	    p->err = gretl_list_add_list(&list1, list2);
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
	    }
	}
	ret->v.ivec = list;
	free(llist);
	free(rlist);
    }

    return ret;
}

static NODE *get_named_list_element (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(l, p);
	int i = node_get_int(r, p);

	if (!p->err) {
	    if (i < 1 || i > list[0]) {
		gretl_errmsg_sprintf("Index value %d is out of bounds", i);
		p->err = E_DATA;
	    } else {
		ret->v.xval = list[i];
	    }
	}
	free(list);
    }

    return ret;
}    

/* Binary operator applied to two bundles: at present only '+'
   (for union) is supported.
*/

static NODE *bundle_op (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_bundle_node(p);

    if (ret != NULL && starting(p)) {
	gretl_bundle *bl = l->v.b;
	gretl_bundle *br = r->v.b;

	if (!p->err) {
	    if (f == B_ADD) {
		ret->v.b = gretl_bundle_union(bl, br, &p->err);
	    } else {
		p->err = E_TYPES;
	    }
	}
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
    NODE *ret = aux_vec_node(p, p->dset->n);

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

static NODE *num_string_comp (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	NODE *xnode = l->t == NUM ? l : r;
	NODE *snode = l->t == STR ? l : r;
	const char *s = snode->v.str;
	int v = xnode->vnum;

	if (v <= 0 || !is_string_valued(p->dset, v)) {
	    p->err = E_TYPES;
	} else {
	    double sx = series_decode_string(p->dset, v, s);

	    if (na(sx)) {
		ret->v.xval = NADBL;
	    } else if (f == B_EQ) {
		ret->v.xval = (sx == xnode->v.xval);
	    } else {
		ret->v.xval = (sx != xnode->v.xval);
	    }
	}
    }

    return ret;
}

/* argument is list; value returned is series */

static NODE *list_to_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);

	if (list != NULL) {
	    p->err = cross_sectional_stat(ret->v.xvec, list,
					  p->dset, f);
	    free(list);
	}
    }

    return ret;
}

/* arguments are series on left, list on right: we add all members 
   of list to series, or subtract all members */

static NODE *series_list_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

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

/* arguments are list, matrix; return is series */

static NODE *list_matrix_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(l, p);
	const gretl_matrix *b = r->v.m;

	if (list != NULL && !gretl_is_null_matrix(b)) {
	    p->err = list_linear_combo(ret->v.xvec, list, b, p->dset);
	}
	free(list);
    }

    return ret;
}

static NODE *list_list_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	int *llist = node_get_list(l, p);
	int *rlist = node_get_list(r, p);

	if (llist != NULL && rlist != NULL) {
	    p->err = x_sectional_weighted_stat(ret->v.xvec,
					       llist, rlist,
					       p->dset, f);
	}
	free(llist);
	free(rlist);
    }

    return ret;
}

/* check for missing obs in a list of variables */

static NODE *list_ok_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);
	int i, v, t;
	double x;

	if (list == NULL || list[0] == 0) {
	    free(list);
	    return ret;
	}

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    x = (f == F_DATAOK)? 1 : 0;
	    for (i=1; i<=list[0]; i++) {
		v = list[i];
		if (na(p->dset->Z[v][t])) {
		    x = (f == F_DATAOK)? 0 : 1;
		    break;
		}
	    }
	    ret->v.xvec[t] = x;
	}

	free(list);
    }

    return ret;
}

/* functions taking (up to) two scalars as arguments and 
   returning a series result */

static NODE *
series_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	double x, y;

	x = (l->t == EMPTY)? NADBL : node_get_scalar(l, p);
	y = (r->t == EMPTY)? NADBL : node_get_scalar(r, p);

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

/* Functions taking two series as arguments and returning a scalar
   or matrix result. We also accept as arguments two matrices if 
   they are vectors of the same length.
*/

static NODE *series_2_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const double *x = NULL, *y = NULL;
	int t1 = 0, t2 = 0;

	if (l->t == VEC && r->t == VEC) {
	    /* two series */
	    x = l->v.xvec;
	    y = r->v.xvec;
	    t1 = p->dset->t1;
	    t2 = p->dset->t2;
	} else {
	    /* two matrices */
	    int n1 = gretl_vector_get_length(l->v.m);
	    int n2 = gretl_vector_get_length(r->v.m);

	    if (n1 == 0 || n1 != n2) {
		p->err = E_TYPES;
	    } else {
		x = l->v.m->val;
		y = r->v.m->val;
		t1 = 0;
		t2 = n1 - 1;
	    }
	}

	if (p->err) {
	    return ret;
	}

	if (f == F_FCSTATS) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = aux_scalar_node(p);
	}

	if (ret == NULL) {
	    return NULL;
	}

	switch (f) {
	case F_COR:
	    ret->v.xval = gretl_corr(t1, t2, x, y, NULL);
	    break;
	case F_COV:
	    ret->v.xval = gretl_covar(t1, t2, x, y, NULL);
	    break;
	case F_FCSTATS:
	    ret->v.m = forecast_stats(x, y, t1, t2, OPT_D, &p->err);
	    break;
	default:
	    break;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

/* takes two series or two matrices as arguments */

static NODE *mxtab_func (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == MAT && r->t == MAT) {
	    ret->v.m = matrix_matrix_xtab(l->v.m, r->v.m, &p->err);
	} else if (l->t == VEC && r->t == VEC) {
	    const double *x = l->v.xvec;
	    const double *y = r->v.xvec;
	    
	    ret->v.m = gretl_matrix_xtab(p->dset->t1, p->dset->t2, 
					 x, y, &p->err);
	} else {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static int type_translate_to_int (GretlType type)
{
    if (type == GRETL_TYPE_DOUBLE) {
	return 1;
    } else if (type == GRETL_TYPE_SERIES) {
	return 2;
    } else if (type == GRETL_TYPE_MATRIX) {
	return 3;
    } else if (type == GRETL_TYPE_STRING) {
	return 4;
    } else if (type == GRETL_TYPE_BUNDLE) {
	return 5;
    } else if (type == GRETL_TYPE_MATRIX_REF) {
	return 6;
    } else {
	return 0;
    }
}

static NODE *object_status (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = n->v.str;

	ret->v.xval = NADBL;

	if (f == F_ISSERIES) {
	    ret->v.xval = gretl_is_series(s, p->dset);
	} else if (f == F_ISLIST) {
	    ret->v.xval = gretl_is_list(s);
	} else if (f == F_ISSTRING) {
	    ret->v.xval = gretl_is_string(s);
	} else if (f == F_ISNULL) {
	    ret->v.xval = 1;
	    if (gretl_is_series(s, p->dset)) {
		ret->v.xval = 0.0;
	    } else if (gretl_is_user_var(s)) {
		ret->v.xval = 0.0;
	    }
	} else if (f == F_OBSNUM) {
	    int t = get_observation_number(s, p->dset);

	    if (t > 0) {
		ret->v.xval = t;
	    }
	} else if (f == F_STRLEN) {
	    /* ret->v.xval = strlen(s); */
	    ret->v.xval = g_utf8_strlen(s, -1);
	} else if (f == F_REMOVE) {
	    ret->v.xval = gretl_remove(s);
	} else if (f == F_TYPEOF) {
	    GretlType type = 0;

	    if (current_series_index(p->dset, s) >= 0) {
		type = GRETL_TYPE_SERIES;
	    } else {
		user_var_get_value_and_type(s, &type);
	    }
	    ret->v.xval = type_translate_to_int(type);
	}
    }

    return ret;
}

/* return scalar node holding the length of the list associated
   with node @n
*/

static NODE *list_length_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	if (n->t == STR) {
	    int *list = get_list_by_name(n->v.str);

	    if (list != NULL) {
		ret->v.xval = list[0];
	    } else {
		node_type_error(F_LISTLEN, 1, LIST, n, p);
	    }
	} else {
	    int *list = node_get_list(n, p);

	    if (list != NULL) {
		ret->v.xval = list[0];
		free(list);
	    }
	}
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
		node_type_error(F_INLIST, 2, VEC, r, p);
	    }
	    if (k >= 0) {
		ret->v.xval = in_gretl_list(list, k);
	    }
	    free(list);
	}
    }	

    return ret;
}

static NODE *argname_from_uvar (NODE *n, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *vname = NULL;

	if (n->t == VEC) {
	    vname = p->dset->varname[n->vnum];
	} else if (n->t == NUM) {
	    vname = n->vname;
	}

	if (vname == NULL) {
	    p->err = E_DATA;
	} else {
	    ret->v.str = gretl_func_get_arg_name(vname, &p->err);
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
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int i;

	if (scalar_node(n)) {
	    i = node_get_int(n, p);
	} else {
	    node_type_error(f, 0, NUM, n, p);
	    return NULL;
	}

	if (f == F_OBSLABEL) {
	    ret->v.str = retrieve_date_string(i, p->dset, &p->err);
	} else if (f == F_VARNAME) {
	    if (i >= 0 && i < p->dset->v) {
		ret->v.str = gretl_strdup(p->dset->varname[i]);
	    } else {
		p->err = E_DATA;
	    }
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	} 	
    }

    return ret;
}

static NODE *list_to_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	int *list = node_get_list(n, p);

	if (p->err) {
	    return ret;
	}

	if (f == F_VARNAME) {
	    ret->v.str = gretl_list_get_names(list, p->dset,
					      &p->err);
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

static NODE *single_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = n->v.str;

	if (f == F_ARGNAME) {
	    ret->v.str = gretl_func_get_arg_name(s, &p->err);
	} else if (f == F_BACKTICK) {
	    ret->v.str = gretl_backtick(s, &p->err);
	} else if (f == F_STRSTRIP) {
	    ret->v.str = gretl_strstrip_copy(s, &p->err);
	} else if (f == F_FIXNAME) {
	    ret->v.str = calloc(VNAMELEN, 1);
	    normalize_join_colname(ret->v.str, s, 0);
	} else {
	    p->err = E_DATA;
	}
    }

    return ret;
}

static NODE *readfile_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *codeset = NULL;

	if (!null_or_empty(r)) {
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

static NODE *two_string_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;

    if (starting(p)) {
	const char *sl = l->v.str;
	const char *sr = r->v.str;

	ret = aux_string_node(p);

	if (ret == NULL) {
	    return NULL;
	}

	if (f == F_STRSTR) {
	    char *sret, *tmp = gretl_strdup(sr);

	    if (tmp != NULL) {
		strstr_escape(tmp);
		sret = strstr(sl, tmp);
		if (sret != NULL) {
		    ret->v.str = gretl_strdup(sret);
		} else {
		    ret->v.str = gretl_strdup("");
		}
		free(tmp);
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
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *one_string_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	char *s;

	if (f == F_TOLOWER) {
	    s = ret->v.str = gretl_strdup(n->v.str);
	    while (s && *s) {
		*s = tolower(*s);
		s++;
	    }
	} else if (f == F_TOUPPER) {
	    s = ret->v.str = gretl_strdup(n->v.str);
	    while (s && *s) {
		*s = toupper(*s);
		s++;
	    }	    
	} else {
	    p->err = E_DATA;
	}

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

static const char *advance_to_split (const char *s)
{
    int n = strcspn(s, " \t\r\n");
    int m = strlen(s);

    return (n == m)? NULL : s + n;
}

static NODE *gretl_strsplit (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = l->v.str;
	int k = node_get_int(r, p);

	if (k < 1) {
	    p->err = E_DATA;
	} else {
	    const char *q;
	    int i;

	    for (i=1; i<k; i++) {
		q = advance_to_split(s);
		if (q != NULL) {
		    q += strspn(q, " \t\r\n");
		    s = q;
		} else {
		    s = "";
		    break;
		}
	    }
	    ret->v.str = gretl_strndup(s, strcspn(s, " \t\r\n"));
	} 

	if (!p->err && ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

static NODE *errmsg_node (NODE *l, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *src = NULL;

	if (l->t == NUM) {
	    int errval = l->v.xval;

	    if (errval < 0 || errval >= E_MAX) {
		p->err = E_DATA;
	    } else {
		src = errmsg_get_with_default(errval);
	    }
	} else {
	    /* empty input node */
	    src = gretl_errmsg_get();
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

static NODE *isodate_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (l->t != NUM && l->t != VEC) {
	node_type_error(F_ISODATE, 1, NUM, l, p);
    } else if (r->t != NUM && r->t != EMPTY) {
	node_type_error(F_ISODATE, 2, NUM, r, p);
    }

    if (!p->err) {
	if (l->t == NUM) {
	    /* epoch day node is scalar */
	    int as_string = (r->t == NUM)? node_get_int(r, p) : 0;

	    if (!p->err) {
		ret = as_string ? aux_string_node(p) : aux_scalar_node(p);
	    }
	    if (ret != NULL) {
		double x = l->v.xval;

		if (!as_string && na(x)) {
		    ret->v.xval = NADBL;
		} else if (x >= 1 && x <= LONG_MAX) {
		    if (as_string) {
			ret->v.str = ymd_extended_from_epoch_day((long) x, &p->err);
		    } else {
			ret->v.xval = ymd_basic_from_epoch_day((long) x, &p->err);
		    }
		} else {
		    p->err = E_INVARG;
		}
	    }
	} else {
	    /* epoch day node is series */
	    ret = aux_vec_node(p, p->dset->n);
	    if (ret != NULL) {
		double xt;
		int t;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    xt = l->v.xvec[t];
		    if (na(xt)) {
			ret->v.xvec[t] = NADBL;
		    } else if (xt >= 1 && xt <= LONG_MAX) {
			ret->v.xvec[t] = ymd_basic_from_epoch_day((long) xt, &p->err);
		    } else {
			p->err = E_INVARG;
			break;
		    }
		}
	    }
	}
    }

    return ret;
}

static NODE *atof_node (NODE *l, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	char *endptr = NULL;

	errno = 0;
	gretl_push_c_numeric_locale();
	ret->v.xval = strtod(l->v.str, &endptr);
	if (errno || endptr == l->v.str) {
	    errno = 0;
	    ret->v.xval = NADBL;
	}
	gretl_pop_c_numeric_locale();
    }

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
	const char *buf = l->v.str;

	if (l->vname == NULL) {
	    gretl_errmsg_set("getline: the source must be a string variable");
	    p->err = E_TYPES;
	} else if (null_or_empty(r)) {
	    bufgets_finalize(buf);
	} else if (r->vname == NULL) {
	    gretl_errmsg_set("getline: the target must be a string variable");
	} else {
	    p->err = query_bufgets_init(buf);
	    if (!p->err) {
		size_t len = bufgets_peek_line_length(buf);

		if (len == 0) {
		    bufgets_finalize(buf);
		    strcpy(r->v.str, "");
		    ret->v.xval = 0;
		} else {
		    r->v.str = user_string_resize(r->vname, len, &p->err);
		    if (!p->err) {
			bufgets(r->v.str, len, buf);
			strip_newline(r->v.str);
			ret->v.xval = 1;
		    }
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
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

static int series_get_end (int t1, int t2, const double *x)
{
    int t;

    for (t=t2; t>=t1; t--) {
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

#define series_cast_optional(f) (f == F_SD) 

static void cast_to_series (NODE *n, int f, gretl_matrix **tmp, 
			    int *t1, int *t2, parser *p)
{
    gretl_matrix *m = n->v.m;
    int len = gretl_vector_get_length(m);

    if (gretl_is_null_matrix(m)) {
	p->err = E_DATA;
    } else if (len == p->dset->n) {
	*tmp = m;
	n->v.xvec = m->val;
    } else {
	if (series_cast_optional(f)) {
	    p->err = E_TYPES; /* temporary error flag */
	} else if (len > 0 && t1 != NULL && t2 != NULL) {
	    *tmp = m;
	    n->v.xvec = m->val;
	    *t1 = 0;
	    *t2 = len - 1;
	} else {
	    node_type_error(f, 1, VEC, n, p);
	}
    } 
}

/* functions taking a series as argument and returning a scalar */

static NODE *
series_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *tmp = NULL;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;
	const double *x;

	if (n->t == MAT) {
	    if (f == F_SUM) {
		/* we'll sum all elements of a matrix */
		tmp = n->v.m;
		t1 = 0;
		t2 = tmp->rows * tmp->cols - 1;
		n->v.xvec = tmp->val;
	    } else if (f == F_T1 || f == F_T2) {
		cast_to_series(n, f, &tmp, NULL, NULL, p);
	    } else {
		cast_to_series(n, f, &tmp, &t1, &t2, p);
	    }
	    if (p->err) {
		if (f == F_SD) {
		    /* offer column s.d. instead */
		    p->err = 0;
		    return matrix_to_matrix_func(n, NULL, f, p);
		} else {
		    return NULL;
		}
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
	case F_T1:
	    ret->v.xval = series_get_start(0, p->dset->n - 1, x);
	    break;
	case F_T2:
	    ret->v.xval = series_get_end(0, p->dset->n - 1, x);
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

/* There must be a matrix in @l; @r may hold a vector or
   a scalar value */

static NODE *matrix_quantiles_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *pmat = NULL;
	int free_pmat = 0;

	if (r->t == MAT) {
	    pmat = r->v.m;
	} else {
	    double x = node_get_scalar(r, p);

	    pmat = gretl_matrix_from_scalar(x);
	    free_pmat = 1;
	} 

	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    ret->v.m = gretl_matrix_quantiles(l->v.m, pmat, &p->err);
	}

	if (free_pmat) {
	    gretl_matrix_free(pmat);
	}
    }

    return ret;
}

/* functions taking a series and a scalar as arguments and returning 
   a scalar
*/

static NODE *
series_scalar_scalar_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	double rval = node_get_scalar(r, p);
	const double *xvec;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;
	int pd = 1;

	if (l->t == MAT) {
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

	ret = aux_scalar_node(p);
	if (p->err) {
	    return ret;
	}

	switch (f) {
	case F_LRVAR:
	    ret->v.xval = gretl_long_run_variance(t1, t2, xvec, (int) rval);
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

static NODE *isconst_node (NODE *l, NODE *r, parser *p)
{
    if (r->t == EMPTY) {
	return series_scalar_func(l, F_ISCONST, p);
    } else if (l->t == MAT) {
	node_type_error(F_ISCONST, 1, VEC, l, p);
	return NULL;
    } else if (!dataset_is_panel(p->dset)) {
	p->err = E_PDWRONG;
	return NULL;
    } else {
	return series_scalar_scalar_func(l, r, F_ISCONST, p);
    }
}

/* Series on left, scalar or string on right, as in
   x[23] or somevar["CA"]. We return the selected
   scalar value from the series.
*/

static NODE *series_obs (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	int t = -1; /* invalid */

	if (r->t == STR) {
	    t = dateton(r->v.str, p->dset);
	} else {
	    /* plain integer */
	    t = node_get_int(r, p);
	    if (!p->err && t > 0 && t <= p->dset->n) {
		t--; /* convert to zero based */
	    }
	}

	if (p->err) {
	    return ret;
	}

	if (t >= 0 && t < p->dset->n) {
	    ret->v.xval = l->v.xvec[t];
	    ret->vnum = l->vnum; /* added 2013-09-14 */
	} else {
	    ret->v.xval = NADBL;
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
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	const double *x = l->v.xvec;
	int order = node_get_int(r, p);

	if (!p->err) {
	    p->err = poly_trend(x, ret->v.xvec, p->dset, order);
	}
    }

    return ret;
}

static NODE *series_movavg (NODE *l, NODE *m, NODE *r, parser *p)
{
    NODE *ret;
    const double *x = l->v.xvec;
    double d = node_get_scalar(m, p);
    int ctrl = 0;
    int k = 0;

    if (p->err) {
	/* from node_get_scalar */
	return NULL;
    }

    if (d < 1.0) {
	/* exponential MA */
	if (d < 0.0) {
	    p->err = E_INVARG;
	    return NULL;
	} else if (dataset_is_panel(p->dset)) {
	    p->err = E_PDWRONG;
	    return NULL;
	}
	k = -1;
    } else {
	/* regular MA */  
	if (d <= 0 || d > INT_MAX) {
	    p->err = E_INVARG;
	    return NULL;
	}
	k = (int) d;
	d = -1.0;
    } 

    if (r->t != EMPTY) {
	/* optional control argument */
	ctrl = node_get_int(r, p);
	if (p->err) {
	    return NULL;
	}
    } else if (d > 0) {
	/* EMA default: initialize with one obs */
	ctrl = 1;
    }

    ret = aux_vec_node(p, p->dset->n);
    if (ret == NULL) {
	return NULL;
    }

    if (d > 0) {
	p->err = exponential_movavg_series(x, ret->v.xvec, 
					   p->dset, d, ctrl);
    } else {
	p->err = movavg_series(x, ret->v.xvec, p->dset, k, ctrl);
    }

    return ret;
}

static NODE *series_lag (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;
    const double *x = l->v.xvec;
    int k = -(node_get_int(r, p));
    int t1, t2;

    if (!p->err) {
	ret = aux_vec_node(p, p->dset->n);
    }

    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    lag_calc(ret->v.xvec, x, k, t1, t2, B_ASN, 1.0, p);

    return ret;
}

static NODE *series_sort_by (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL && starting(p)) {
	if (l->t == VEC && r->t == VEC) {
	    p->err = gretl_sort_by(l->v.xvec, r->v.xvec, ret->v.xvec, p->dset); 
	} else {
	    p->err = E_TYPES;
	}
    } 

    return ret;
}

static NODE *vector_sort (NODE *l, int f, parser *p)
{
    NODE *ret = (l->t == VEC)? aux_vec_node(p, p->dset->n) :
	aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == VEC) {
	    p->err = sort_series(l->v.xvec, ret->v.xvec, f, p->dset); 
	} else if (gretl_is_null_matrix(l->v.m)) {
	    p->err = E_DATA;
	} else {
	    int n = gretl_vector_get_length(l->v.m);

	    if (n > 0) {
		ret->v.m = gretl_matrix_copy(l->v.m);
		if (ret->v.m == NULL) {
		    p->err = E_ALLOC;
		} else {
		    double *x = ret->v.m->val;

		    qsort(x, n, sizeof *x, (f == F_SORT)? gretl_compare_doubles :
			  gretl_inverse_compare_doubles);
		}
	    } else {
		p->err = E_TYPES;
	    }
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

	if (l->t == VEC) {
	    n = sample_size(p->dset);
	    x = l->v.xvec + p->dset->t1;
	} else if (!gretl_is_null_matrix(l->v.m)) {
	    n = gretl_vector_get_length(l->v.m);
	    x = l->v.m->val;
	}

	if (n > 0 && x != NULL) {
	    gretlopt opt = (f == F_VALUES)? OPT_S : OPT_NONE;

	    ret->v.m = gretl_matrix_values(x, n, opt, &p->err);
	} else {
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

	if (l->t == VEC) {
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

/* Takes a series as argument and returns a matrix: 
   right now only F_FREQ does this 
*/

static NODE *series_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL) {
	gretl_matrix *tmp = NULL;
	int t1 = p->dset->t1;
	int t2 = p->dset->t2;

	if (n->t == MAT) {
	    cast_to_series(n, f, &tmp, &t1, &t2, p);
	} 

	if (!p->err) {
	    ret->v.m = freqdist_matrix(n->v.xvec, t1, t2, &p->err);
	    if (n->t == MAT) {
		/* restore matrix on @n after "cast" above */
		n->v.m = tmp;
	    }
	}
    }

    return ret;
}

static int get_logtrans (const char *s)
{
    if (s != NULL) {
	if (*s != 'T' && strchr(s, 'l')) {
	    return 1;
	}
    }

    return 0;
}

#define use_tramo(s) (s != NULL && (s[0] == 't' || s[0] == 'T'))

#define is_panel_stat(f) (f == F_PNOBS || \
			  f == F_PMIN ||  \
			  f == F_PMAX ||  \
			  f == F_PSUM || \
			  f == F_PMEAN || \
			  f == F_PXSUM ||  \
			  f == F_PSD)

/* Functions taking a series as argument and returning a series.
   Note that the 'r' node may contain an auxiliary parameter;
   in that case the aux value should be a scalar, unless 
   we're doing F_DESEAS, in which case it should be a string,
   or one of the panel stats functions, in which case it should
   be a series.
*/

static NODE *series_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    int rtype = NUM; /* the optional right-node type */

    if (f == F_SDIFF && !dataset_is_seasonal(p->dset)) {
	p->err = E_PDWRONG;
	return NULL;
    }

    if (f == F_DESEAS) {
	rtype = STR;
    } else if (is_panel_stat(f)) {
	rtype = VEC;
    } 

    if (null_or_empty(r)) {
	rtype = 0; /* OK */
    } else if (r->t != rtype) {
	node_type_error(f, 2, rtype, r, p);
	return NULL;
    }

    ret = aux_vec_node(p, p->dset->n);

    if (ret != NULL) {
	gretl_matrix *tmp = NULL;
	double parm = NADBL;
	const double *z = NULL;
	const double *x;
	double *y;

	if (l->t == MAT) {
	    cast_to_series(l, f, &tmp, NULL, NULL, p);
	}

	if (rtype == VEC) {
	    z = r->v.xvec;
	} else if (rtype == NUM) {
	    parm = node_get_scalar(r, p);
	}

	if (p->err) {
	    return NULL;
	}	

	x = l->v.xvec;
	y = ret->v.xvec;

	switch (f) {
	case F_HPFILT:
	    p->err = hp_filter(x, y, p->dset, parm, OPT_NONE);
	    break;
	case F_FRACDIFF:
	    p->err = fracdiff_series(x, y, parm, 1, (autoreg(p))? p->obs : -1, p->dset);
	    break;
	case F_FRACLAG:
	    p->err = fracdiff_series(x, y, parm, 0, (autoreg(p))? p->obs : -1, p->dset);
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
	    if (rtype == STR) {
		int tramo = use_tramo(r->v.str);
		int logt = get_logtrans(r->v.str);

		p->err = seasonally_adjust_series(x, y, p->dset, tramo, logt); 
	    } else {
		p->err = seasonally_adjust_series(x, y, p->dset, 0, 0);
	    }
	    break;
	case F_RESAMPLE:
	    if (rtype == NUM) {
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

static NODE *do_panel_shrink (NODE *l, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.m = panel_shrink(l->v.xvec, p->dset, &p->err);
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

    if (!empty_or_num(r)) {
	/* optional 'r' node must be scalar */
	node_type_error(F_PERGM, 2, NUM, r, p);
    } else if (l->t == MAT && gretl_vector_get_length(l->v.m) == 0) {
	/* if 'l' node is not a series, must be a vector */
	node_type_error(F_PERGM, 1, VEC, l, p);
    } else {
	ret = aux_matrix_node(p);
    }

    if (!p->err) {
	const double *x = NULL;
	int t1 = 0, t2 = 0;
	int width = -1;

	if (l->t == VEC) {
	    x = l->v.xvec;
	    t1 = p->dset->t1;
	    t2 = p->dset->t2;
	} else if (l->t == MAT) {
	    x = l->v.m->val;
	    t1 = 0;
	    t2 = gretl_vector_get_length(l->v.m) - 1;
	} 

	if (r != NULL && r->t == NUM) {
	    width = r->v.xval;
	}

	ret->v.m = periodogram_matrix(x, t1, t2, width, &p->err);
    }

    return ret;
}

/* application of scalar function to each element of matrix */

static NODE *apply_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;
	int i, n = m->rows * m->cols;
	double x;

	if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
	    return NULL;
	}

	for (i=0; i<n && !p->err; i++) {
	    /* FIXME error handling? */ 
	    x = real_apply_func(m->val[i], f, p);
	    ret->v.m->val[i] = x;
	}
    }

    return ret;
}

/* node holding a cashed-out string variable */

static NODE *string_var_node (NODE *t, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL && starting(p)) {
	const char *sval, *sname = t->v.str;

	if (*sname == '@') sname++;
	sval = get_string_by_name(sname);
	if (sval == NULL) {
	    p->err = E_UNKVAR;
	} else {
	    ret->v.str = gretl_strdup(sval);
	    if (ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }
	}
    }

    return ret;
}

static gretl_matrix *matrix_from_scalars (NODE *t, int m,
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
	    n = t->v.bn.n[i];
	    if (n->t == EMPTY) {
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
		n = t->v.bn.n[k++];
		if (n->t == EMPTY) {
		    n = t->v.bn.n[k++];
		}
		x = node_get_scalar(n, p);
		if (na(x)) {
		    x = M_NA;
		}
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
    gretl_matrix *M;

    if (list != NULL && list[0] == 0) {
	M = gretl_null_matrix_new();
    } else {
	const gretl_matrix *mmask = get_matrix_mask();

	if (mmask != NULL) {
	    M = gretl_matrix_data_subset_special(list, dset, 
						 mmask, &p->err);
	} else {
	    int missop = (libset_get_bool(SKIP_MISSING))? M_MISSING_SKIP :
		M_MISSING_OK;

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

static void *arg_get_data (NODE *n, int ref, GretlType *type)
{
    void *data = NULL;

    if (n->t == VEC) {
	if (ref) {
	    *type = GRETL_TYPE_SERIES_REF;
	    data = &n->vnum;
	} else if (n->vname != NULL) {
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
    } else if (n->t == STR) {
	*type = GRETL_TYPE_STRING;
	data = n->v.str;
    } else if (n->t == LIST) {
	*type = GRETL_TYPE_LIST;
	data = n->v.ivec;
    } else {
	*type = GRETL_TYPE_NONE;
    }

    return data;
}

static int check_uaddr_type (NODE *u, parser *p)
{
    NODE *n = u->v.b1.b;
    int err = 0;

    if (uscalar_node(n)) {
	; /* OK */
    } else if (useries_node(n)) {
	; /* OK */
    } else if (umatrix_node(n)) {
	; /* OK */
    } else if (ubundle_node(n)) {
	; /* OK */	
    } else {
	pputs(p->prn, _("Wrong type of operand for unary '&'"));
	pputc(p->prn, '\n');
	err = E_TYPES;
    }

    return err;
}

#define ok_ufunc_sym(s) (s == NUM || s == VEC || s == MAT || \
                         s == LIST || s == U_ADDR || s == DUM || \
                         s == STR || s == EMPTY || s == BUNDLE)

/* evaluate a user-defined function */

static NODE *eval_ufunc (NODE *t, parser *p)
{
    NODE *l = t->v.b2.l;
    NODE *r = t->v.b2.r;
    NODE *ret = NULL;
    const char *funname = l->v.str;
    fnargs *args = NULL;
    ufunc *uf = NULL;
    int i, nparam, argc = 0;
    int rtype = GRETL_TYPE_NONE;

    /* first find the function */
    uf = get_user_function_by_name(funname);
    if (uf == NULL) {
	fprintf(stderr, "%s: couldn't find a function of this name\n", funname);
	p->err = E_DATA;
    }

    if (!p->err && !(p->flags & P_VOID)) {
	/* check that the function returns something suitable */
	rtype = user_func_get_return_type(uf);
	if (!ok_function_return_type(rtype) || rtype == GRETL_TYPE_VOID) {
	    fprintf(stderr, "eval_ufunc: %s: invalid return type %d\n", 
		    funname, rtype);
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

    if (!p->err) {
	/* allocate an arguments array */
	args = fn_args_new(argc);
	if (args == NULL) {
	    p->err = E_ALLOC;
	}
    }

    if (p->err) {
	/* no sense in continuing */
	return NULL;
    }

    /* evaluate the function argument nodes */

    for (i=0; i<argc && !p->err; i++) {
	NODE *arg, *ni = r->v.bn.n[i];
	GretlType argt = 0;
	int reftype = 0;
	void *data;

	if (starting(p)) {
	    /* evaluate all nodes */
	    arg = eval(ni, p);
	} else if (ni->vname != NULL) {
	    /* otherwise let named variables through "as is" */
	    arg = ni;
	} else {
	    arg = eval(ni, p);
	}

	if (p->err || arg == NULL) {
	    fprintf(stderr, "%s: failed to evaluate arg %d\n", funname, i);
	    p->err = (p->err == 0)? E_DATA : p->err;
	} else if (!ok_ufunc_sym(arg->t)) {
	    gretl_errmsg_sprintf("%s: invalid argument type %s", funname, 
				 typestr(arg->t));
	    p->err = E_TYPES;
	}

#if EDEBUG
	fprintf(stderr, "%s: arg %d is of type %d\n", funname, i, 
		arg == NULL? -1 : arg->t);
#endif

	if (!p->err && arg->t == U_ADDR) {
	    /* address node: switch to the 'content' sub-node */
	    reftype = 1;
	    arg = arg->v.b1.b;
	}

	if (!p->err && arg->t == DUM && arg->v.idnum != DUM_NULL) {
	    p->err = E_TYPES;
	}

	if (!p->err) {
	    /* assemble info and push argument */
	    data = arg_get_data(arg, reftype, &argt);
	    p->err = push_fn_arg(args, arg->vname, argt, data);
	}

	if (p->err) {
	    fprintf(stderr, "%s: error evaluating arg %d\n", funname, i);
	}
    }

    /* try sending args to function */

    if (!p->err) {
	char *descrip = NULL;
	char **pdescrip = NULL;
	double xret = NADBL;
	double *Xret = NULL;
	gretl_matrix *mret = NULL;
	gretl_bundle *bret = NULL;
	char *sret = NULL;
	int *iret = NULL;
	void *retp = NULL;

	if (rtype == GRETL_TYPE_DOUBLE) {
	    retp = &xret;
	} else if (rtype == GRETL_TYPE_SERIES) {
	    retp = &Xret;
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    retp = &mret;
	} else if (rtype == GRETL_TYPE_LIST) {
	    retp = &iret;
	} else if (rtype == GRETL_TYPE_STRING) {
	    retp = &sret;
	} else if (rtype == GRETL_TYPE_BUNDLE) {
	    retp = &bret;
	}

	if ((p->flags & P_UFRET) && 
	    (rtype == GRETL_TYPE_DOUBLE || rtype == GRETL_TYPE_SERIES)) {
	    /* pick up description of generated var, if any */
	    pdescrip = &descrip;
	}

	p->err = gretl_function_exec(uf, args, rtype, p->dset, 
				     retp, pdescrip, p->prn);

	if (!p->err) {
	    if (rtype == GRETL_TYPE_DOUBLE) {
		ret = aux_scalar_node(p);
		if (ret != NULL) {
		    ret->v.xval = xret;
		}
	    } else if (rtype == GRETL_TYPE_SERIES) {
		ret = aux_vec_node(p, 0);
		if (ret != NULL) {
		    if (ret->v.xvec != NULL) {
			free(ret->v.xvec);
		    }
		    ret->v.xvec = Xret;
		}
	    } else if (rtype == GRETL_TYPE_MATRIX) {
		ret = aux_matrix_node(p);
		if (ret != NULL) {
		    if (is_tmp_node(ret)) {
			gretl_matrix_free(ret->v.m);
		    }
		    ret->v.m = mret;
		}
	    } else if (rtype == GRETL_TYPE_LIST) {
		ret = aux_list_node(p);
		if (ret != NULL) {
		    if (is_tmp_node(ret)) {
			free(ret->v.ivec);
		    }
		    ret->v.ivec = iret;
		}
	    } else if (rtype == GRETL_TYPE_STRING) {
		ret = aux_string_node(p);
		if (ret != NULL) {
		    if (is_tmp_node(ret)) {
			free(ret->v.str);
		    }
		    ret->v.str = sret;
		}
	    } else if (rtype == GRETL_TYPE_BUNDLE) {
		ret = aux_bundle_node(p);
		if (ret != NULL) {
		    if (is_tmp_node(ret)) {
			gretl_bundle_destroy(ret->v.b);
		    }
		    ret->t = BUNDLE;
		    ret->v.b = bret;
		}
	    }		
	}

	if (descrip != NULL) {
	    strcpy(p->lh.label, descrip);
	    free(descrip);
	}
    }

    fn_args_free(args);

#if EDEBUG
    fprintf(stderr, "eval_ufunc: p->err = %d, ret = %p\n", 
	    p->err, (void *) ret);
#endif

    return ret;
}

#ifdef USE_RLIB

/* evaluate an R function */

static NODE *eval_Rfunc (NODE *t, parser *p)
{
    NODE *l = t->v.b2.l;
    NODE *r = t->v.b2.r;
    int i, argc = r->v.bn.n_nodes;
    const char *funname = l->v.str;
    int rtype = GRETL_TYPE_NONE;
    NODE *ret = NULL;

    /* first find the function */
    p->err = gretl_R_get_call(funname, argc);
    if (p->err) {
	return NULL;
    }

    /* evaluate the function arguments */
    for (i=0; i<argc && !p->err; i++) {
	NODE *arg = eval(r->v.bn.n[i], p);

	if (arg == NULL) {
	    fprintf(stderr, "%s: failed to evaluate arg %d\n", funname, i); 
	} else if (!ok_ufunc_sym(arg->t)) {
	    fprintf(stderr, "%s: node type %d: not OK\n", funname, arg->t);
	    p->err = E_TYPES;
	}

	if (p->err) {
	    break;
	}

#if EDEBUG
	fprintf(stderr, "%s: arg[%d] is of type %d\n", funname, i, arg->t);
#endif

	if (arg->t == NUM) {
	    p->err = gretl_R_function_add_scalar(arg->v.xval);
	} else if (arg->t == VEC) {
	    gretl_matrix *m = tmp_matrix_from_series(arg, p);

	    if (m != NULL) {
		p->err = gretl_R_function_add_matrix(m);
		gretl_matrix_free(m);
	    }
	} else if (arg->t == MAT) {
	    p->err = gretl_R_function_add_matrix(arg->v.m);
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
		    if (is_tmp_node(ret)) {
			gretl_matrix_free(ret->v.m);
		    }
		    ret->v.m = (gretl_matrix *) retp;
		}
	    } 
	}
    }

#if EDEBUG
    fprintf(stderr, "eval_Rfunc: p->err = %d, ret = %p\n", 
	    p->err, (void *) ret);
#endif

    return ret;
}

#endif

static gretl_matrix *complex_array_to_matrix (cmplx *c, int sz,
					      parser *p)
{
    gretl_matrix *m = NULL;
    int i, n = sz / sizeof *c;

    if (n <= 0) {
	p->err = E_DATA;
    } else {
	m = gretl_matrix_alloc(n, 2);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		gretl_matrix_set(m, i, 0, c[i].r);
		gretl_matrix_set(m, i, 1, c[i].i);
	    }
	}
    }

    return m;
}

static void process_bundle_member_subspec (NODE *n, parser *p,
					   char *spec)
{
    if (n->t == MAT) {
	gretl_matrix *m = NULL;
	gchar *tmp;

	tmp = g_strdup_printf("gentmpmat___%s", spec);
	p->err = private_matrix_add(n->v.m, "gentmpmat___");
	if (!p->err) {
	    m = generate_matrix(tmp, p->dset, &p->err);
	}
	user_var_delete_by_name("gentmpmat___", NULL);
	g_free(tmp);
	if (!p->err) {
	    /* replace the outgoing matrix */
	    n->v.m = m;
	    n->flags &= ~PTR_NODE;
	    n->flags |= TMP_NODE;
	}
    } else if (n->t == VEC) {
	char *s = strrchr(spec, ']');
	double xt;
	int t;

	*s = '\0';
	t = get_observation_number(spec + 1, p->dset);
	if (t > 0) {
	    /* revise outgoing node to hold scalar */
	    xt = n->v.xvec[t-1];
	    free(n->v.xvec);
	    n->t = NUM;
	    n->v.xval = xt;
	    n->flags &= ~TMP_NODE;
	} else {
	    p->err = E_DATA;
	}
    } else {
	context_error('[', p);
    }
}

/* Getting an object from within a bundle: on the left is the
   bundle reference, on the right should be a string -- the
   key to look up to get content. 
*/

static NODE *get_named_bundle_value (NODE *l, NODE *r, parser *p)
{
    char *key = r->v.str;
    char *s, *extra = NULL;
    GretlType type;
    int size = 0;
    void *val = NULL;
    NODE *ret = NULL;

#if EDEBUG
    fprintf(stderr, "get_named_bundle_value: %s[\"%s\"]\n", l->vname, key);
#endif

    if ((s = strchr(key, '[')) != NULL) {
	/* a matrix slice spec is appended to the key? */
	extra = gretl_strdup(s);
	*s = '\0';
    }

    if (!strcmp(l->vname, "$")) {
	/* special: treat the 'last model' as a bundle */
	val = last_model_get_data(key, &type, &size, &p->err);
    } else {
	val = gretl_bundle_get_data(l->v.b, key, &type, &size, &p->err);
    }

    if (p->err) {
	return ret;
    }

    if (type == GRETL_TYPE_INT) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    int *ip = val;
		
	    ret->v.xval = *ip;
	}    
    } else if (type == GRETL_TYPE_DOUBLE) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    double *dp = val;
		
	    ret->v.xval = *dp;
	}
    } else if (type == GRETL_TYPE_STRING) {
	ret = aux_string_node(p);
	if (ret != NULL) {
	    ret->v.str = gretl_strdup((char *) val);
	    if (ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }
	}
    } else if (type == GRETL_TYPE_MATRIX) {
	ret = matrix_pointer_node(p);
	if (ret != NULL) {
	    ret->v.m = (gretl_matrix *) val;
	    ret->flags &= ~TMP_NODE; /* don't free content! */
	}
    } else if (type == GRETL_TYPE_MATRIX_REF) {
	ret = matrix_pointer_node(p);
	if (ret != NULL) {
	    ret->v.m = (gretl_matrix *) val;
	    ret->flags |= PTR_NODE;
	}
    } else if (type == GRETL_TYPE_BUNDLE) {
	ret = aux_bundle_node(p);
	if (ret != NULL) {
	    ret->v.b = gretl_bundle_copy((gretl_bundle *) val,
					 &p->err);
	} 
    } else if (type == GRETL_TYPE_CMPLX_ARRAY) {
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    ret->v.m = complex_array_to_matrix((cmplx *) val, size, p);
	}
    } else if (type == GRETL_TYPE_SERIES) {
	const double *x = val;

	if (size == p->dset->n) {
	    ret = aux_vec_node(p, p->dset->n);
	    if (ret != NULL) {
		int t;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    ret->v.xvec[t] = x[t];
		}
	    }
	} else if (size > 0) {
	    ret = aux_matrix_node(p);
	    if (ret != NULL) {
		ret->v.m = gretl_vector_from_array(x, size, 
						   GRETL_MOD_NONE);
		if (ret->v.m == NULL) {
		    p->err = E_ALLOC;
		}
	    }
	} else {
	    p->err = E_DATA;
	}
    } else {
	p->err = E_DATA;
    }

    if (extra != NULL) {
	process_bundle_member_subspec(ret, p, extra);
	free(extra);
	*s = '['; /* reinstate full incoming string */
    }

    return ret;
}

static NODE *test_bundle_key (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	gretl_bundle *bundle = l->v.b;
	const char *key = r->v.str;
	GretlType type = 0;

	gretl_bundle_get_data(bundle, key, &type, NULL, NULL);
	ret->v.xval = type_translate_to_int(type);
    }

    return ret;
}

static NODE *type_string_node (NODE *n, parser *p)
{
    NODE *ret = aux_string_node(p);

    if (ret != NULL) {
	int t = n->v.xval;
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
	    s = "matrixref";
	}

	ret->v.str = gretl_strdup(s);
	if (ret->v.str == NULL) {
	    p->err = E_ALLOC;
	}
    }

    return ret;
}

/* Setting an object in a bundle under a given key string. We get here
   only if p->lh.substr is non-NULL. That "substr" may be a string
   literal, or it may be the name of a string variable. In the latter
   case we wait till this point to cash out the string, since we may
   be in a context (e.g. a loop) where the value of the string
   variable changes from one invocation of the generator to the
   next.
*/

static int set_named_bundle_value (const char *name, NODE *n, parser *p)
{
    gretl_bundle *bundle;
    GretlType type;
    void *ptr = NULL;
    char *key = NULL;
    int free_key = 0;
    int size = 0;
    int err = 0;

    bundle = get_bundle_by_name(name);

    if (bundle == NULL) {
	err = E_UNKVAR;
    } else if (p->targ == BMEMB) {
	key = p->lh.substr;
    } else if (*p->lh.substr == '"') {
	key = gretl_strdup(p->lh.substr);
	if (key == NULL) {
	    err = E_ALLOC;
	} else {
	    gretl_unquote(key, &err);
	    free_key = 1;
	}
    } else if (gretl_is_string(p->lh.substr)) {
	key = (char *) get_string_by_name(p->lh.substr);
    } else {
	err = E_DATA;
    }

#if EDEBUG
    fprintf(stderr, "set_named_bundle_value: %s[\"%s\"]\n", name, key);
#endif

    if (!err && p->op != B_ASN) {
	/* e.g. bundle.member += foo: note that right now we don't
	   reach here: for bundles and their members, we choke off 
	   modified assignment at the initial parse stage.
	   AC, 2013-05-12
	*/
	type = gretl_bundle_get_type(bundle, key, &err);
    }

    if (!err) {
	switch (n->t) {
	case NUM:
	    ptr = &n->v.xval;
	    type = GRETL_TYPE_DOUBLE;
	    break;
	case STR:
	    ptr = n->v.str;
	    type = GRETL_TYPE_STRING;
	    break;
	case MAT:
	    ptr = n->v.m;
	    type = GRETL_TYPE_MATRIX;
	    break;
	case U_ADDR:
	    n = n->v.b1.b;
	    if (umatrix_node(n)) {
		ptr = n->v.m;
		type = GRETL_TYPE_MATRIX_REF;
	    } else {
		err = E_TYPES;
	    }	 
	    break;
	case VEC:
	    ptr = n->v.xvec;
	    type = GRETL_TYPE_SERIES;
	    size = p->dset->n;
	    break;
	case BUNDLE:
	    ptr = n->v.b;
	    type = GRETL_TYPE_BUNDLE;
	    err = p->err;
	    break;
	default:
	    err = E_DATA;
	    break;
	}
    }

    if (!err) {
	err = gretl_bundle_set_data(bundle, key, ptr, type, size);
    }

    if (free_key) {
	free(key);
    }

    return err;
}

static int set_named_bundle_note (const char *name, const char *key,
				  const char *note, parser *p)
{
    gretl_bundle *bundle;
    int err = 0;

    bundle = get_bundle_by_name(name);
    if (bundle == NULL) {
	p->err = E_UNKVAR;
    } else {
	err = gretl_bundle_set_note(bundle, key, note);
    }

    return err;
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

    /* if we're supposed to have a list, check that we
       actually have one */
    if (l->t != VEC && l->t != MAT) {
	list = node_get_list(l, p);
	if (p->err) {
	    return NULL;
	}
    }

    /* if third node is matrix, must be col vector */
    if (r->t == MAT) {
	if (r->v.m->cols != 1) {
	    p->err = E_NONCONF;
	    return NULL;
	}
    }

    if (!xcf) {
	/* acf/pacf */
	if (l->t == VEC) {
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
	int xtype = VEC;
	
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

    free(list);

    return A;
}

static const char *ptr_node_get_matrix_name (NODE *t, parser *p)
{
    const char *name = NULL;

    if (t->t == U_ADDR) {
	NODE *n = t->v.b1.b;

	if (umatrix_node(n)) {
	    name = n->vname;
	} else {
	    p->err = E_TYPES;
	}
    } 

    return name;
}

static gretl_matrix *get_density_matrix (const double *x, 
					 const DATASET *dset, 
					 double bws, int ctrl,
					 int *err)
{
    gretl_matrix *(*kdfunc) (const double *, const DATASET *,
			     double, gretlopt, int *);
    void *handle;
    gretl_matrix *m = NULL;
    gretlopt opt;

    kdfunc = get_plugin_function("kernel_density_matrix", &handle);
    if (kdfunc == NULL) {
	*err = E_FOPEN;
	return NULL;
    }

    opt = ctrl ? OPT_O : OPT_NONE; 
    m = (*kdfunc)(x, dset, bws, opt, err);
    close_plugin(handle);

    return m;
}

static int aggregate_discrete_check (const int *list, const DATASET *dset)
{
    int i, vi;

    for (i=1; i<=list[0]; i++) {
	vi = list[i];
	if (!series_is_discrete(dset, vi) && 
	    !gretl_isdiscrete(dset->t1, dset->t2, dset->Z[vi])) {
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

/* evaluate a built-in function that has three arguments */

static NODE *eval_3args_func (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;
    gretl_matrix *A = NULL;
    int post_process = 1;

    if (f == F_MSHAPE) {
	if (l->t != MAT && l->t != NUM) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int k1 = node_get_int(m, p);
	    int k2 = node_get_int(r, p);

	    if (!p->err) {
		if (l->t == NUM) {
		    A = mshape_scalar(l->v.xval, k1, k2, &p->err);
		} else {
		    A = gretl_matrix_shape(l->v.m, k1, k2);
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
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != U_ADDR && m->t != EMPTY) {
	    node_type_error(f, 2, U_ADDR, m, p);
	} else if (r->t != U_ADDR && r->t != EMPTY) {
	    node_type_error(f, 3, U_ADDR, r, p);
	} else {
	    const char *uname, *vname;

	    uname = ptr_node_get_matrix_name(m, p);
	    vname = ptr_node_get_matrix_name(r, p);
	    A = user_matrix_SVD(l->v.m, uname, vname, &p->err);
	}
    } else if (f == F_TOEPSOLV || f == F_VARSIMUL) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != MAT) {
	    node_type_error(f, 3, MAT, r, p);
	} else {
	    if (f == F_TOEPSOLV) {
		A = gretl_toeplitz_solve(l->v.m, m->v.m, r->v.m, &p->err);
	    } else { 
		A = gretl_matrix_varsimul(l->v.m, m->v.m, r->v.m, &p->err);
	    }
	} 
    } else if (f == F_CORRGM) {
	if (l->t != VEC && l->t != MAT && !ok_list_node(l)) {
	    node_type_error(f, 1, VEC, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != VEC && r->t != MAT) {
	    node_type_error(f, 3, VEC, r, p);
	} else {
	    A = get_corrgm_matrix(l, m, r, p);
	}
    } else if (f == F_SEQ) {
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r) && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int start = node_get_int(l, p);
	    int end = node_get_int(m, p);
	    int step = (r->t != EMPTY)? node_get_int(r, p) : 1;

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
	} else if (!empty_or_num(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		if (r != NULL && r->t == NUM) {
		    ret->v.xval = strncmp(l->v.str, m->v.str, 
					  (int) r->v.xval);
		} else {
		    ret->v.xval = strcmp(l->v.str, m->v.str);
		}
	    }
	}
    } else if (f == F_WEEKDAY) {
	post_process = 0;
	if (l->t != NUM) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		int yr = l->v.xval;
		int mo = m->v.xval;
		int day = r->v.xval;

		ret->v.xval = day_of_week(yr, mo, day, &p->err);
	    }
	}
    } else if (f == F_KDENSITY) {
	if (l->t != VEC) {
	    node_type_error(f, 1, VEC, l, p);
	} else if (m->t != NUM && m->t != EMPTY) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    const double *x = l->v.xvec;
	    double bws = (m->t != EMPTY)? m->v.xval : 1.0;
	    int ctrl = (r->t != EMPTY)? (int) r->v.xval : 0;

	    A = get_density_matrix(x, p->dset, bws,
				   ctrl, &p->err);
	}
    } else if (f == F_MONTHLEN) {
	post_process = 0;
	if (!scalar_node(l)) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (!scalar_node(r)) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    int mo = node_get_int(l, p);
	    int yr = node_get_int(m, p);
	    int wk = node_get_int(r, p);

	    if (p->err) {
		; /* from node_get_int() */
	    } else if (yr < 0 || mo < 1 || mo > 12 ||
		(wk != 5 && wk != 6 && wk != 7)) {
		p->err = E_INVARG;
	    } else {
		ret = aux_scalar_node(p);
		if (!p->err) {
		    ret->v.xval = get_days_in_month(mo, yr, wk);
		}
	    }
	}
    } else if (f == F_SETNOTE) {
	post_process = 0;
	if (l->t != BUNDLE) {
	    node_type_error(f, 1, BUNDLE, l, p);
	} else if (m->t != STR) {
	    node_type_error(f, 2, STR, m, p);
	} else if (r->t != STR) {
	    node_type_error(f, 3, STR, r, p);
	} else {
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		ret->v.xval = set_named_bundle_note(l->vname, m->v.str, 
						    r->v.str, p);
	    }
	}
    } else if (f == F_BWFILT) {
	gretl_matrix *tmp = NULL;

	post_process = 0;
	if (l->t != VEC) {
	    if (l->t == MAT) {
		cast_to_series(l, f, &tmp, NULL, NULL, p);
	    } else {
		node_type_error(f, 1, VEC, l, p);
	    }
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_vec_node(p, p->dset->n);
	    if (!p->err) {
		p->err = butterworth_filter(l->v.xvec, ret->v.xvec, p->dset,
					    m->v.xval, r->v.xval);
	    }
	}
	if (tmp != NULL) {
	    l->v.m = tmp;
	}
    } else if (f == F_CHOWLIN) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != MAT && r->t != EMPTY) {
	    node_type_error(f, 3, MAT, r, p);
	} else {
	    const gretl_matrix *X = (r->t == MAT)? r->v.m : NULL;

	    A = matrix_chowlin(l->v.m, X, m->v.xval, &p->err);
	}
    } else if (f == F_IRF) {
	if (l->t != NUM) {
	    node_type_error(f, 1, NUM, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double alpha = (r->t == NUM)? r->v.xval : 0.0;
	    int targ = (int) l->v.xval - 1;
	    int shock = (int) m->v.xval - 1;

	    A = last_model_get_irf_matrix(targ, shock, alpha,
					  p->dset, &p->err);
	}
    } else if (f == F_MLAG) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT && m->t != NUM) {
	    /* scalar or vector */
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != NUM && r->t != EMPTY) {
	    /* optional scalar */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    double missval = (r->t == NUM)? r->v.xval : 0.0;
	    gretl_matrix *mm;

	    if (m->t == NUM) {
		/* promote arg2 if scalar */
		mm = make_scalar_matrix(m->v.xval);
	    } else {
		mm = m->v.m;
	    }
	    A = gretl_matrix_lag(l->v.m, mm, missval);
	    if (m->t == NUM && mm != NULL) {
		gretl_matrix_free(mm);
	    }
	}
    } else if (f == F_EIGSOLVE) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != MAT) {
	    node_type_error(f, 2, MAT, m, p);
	} else if (r->t != EMPTY && r->t != U_ADDR) {
	    /* optional matrix-pointer */
	    node_type_error(f, 3, U_ADDR, r, p);
	} else {
	    const char *rname;

	    rname = (r->t == U_ADDR)? ptr_node_get_matrix_name(r, p) : "null";
	    A = user_gensymm_eigenvals(l->v.m, m->v.m, rname, &p->err);
	}
    } else if (f == F_NADARWAT) {
	post_process = 0;
	if (l->t != VEC) {
	    node_type_error(f, 1, VEC, l, p);
	} else if (m->t != VEC) {
	    node_type_error(f, 2, VEC, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    ret = aux_vec_node(p, p->dset->n);
	    if (!p->err) {
		p->err = nadaraya_watson(l->v.xvec, m->v.xvec,
					 r->v.xval, p->dset, 
					 ret->v.xvec);
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
	    int cov = (r->t == EMPTY)? 0 : node_get_int(r, p);
	    int k = node_get_int(m, p);

	    if (!p->err) {
		A = gretl_matrix_pca(l->v.m, k, 
				     cov ? OPT_C : OPT_NONE,
				     &p->err);
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
	    int offset = (r->t == EMPTY)? 10 : node_get_int(r, p);
	    int rows = node_get_int(l, p);
	    int cols = node_get_int(m, p);

	    if (!p->err) {
		A = halton_matrix(rows, cols, offset, &p->err);
	    }
	}
    } else if (f == F_IWISHART) {
	if (l->t != MAT) {
	    node_type_error(f, 1, MAT, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != EMPTY && r->t != NUM) {
	    /* optional number of replications */
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    const gretl_matrix *S = l->v.m;
	    int v = node_get_int(m, p);
	    int N = (r->t == EMPTY)? 0 : node_get_int(r, p);

	    if (!p->err) {
		if (N == 0) {
		    A = inverse_wishart_matrix(S, v, &p->err);
		} else {
		    A = inverse_wishart_sequence(S, v, N, &p->err);
		}
	    }
	}
    } else if (f == F_AGGRBY) {
	if (l->t != VEC && l->t != LIST && !null_or_empty(l)) {
	    node_type_error(f, 1, VEC, l, p);
	} else if (m->t != VEC && m->t != LIST) {
	    node_type_error(f, 2, VEC, m, p);
	} else if (r->t != STR && !null_or_empty(r)) {
	    node_type_error(f, 3, STR, r, p);
	} else {
	    const char *fncall = NULL;
	    const double *x = NULL;
	    const double *y = NULL;
	    const int *xlist = NULL;
	    const int *ylist = NULL;

	    if (r->t == STR) {
		fncall = r->v.str;
	    }
	    if (l->t == VEC) {
		x = l->v.xvec;
	    } else if (l->t == LIST) {
		xlist = l->v.ivec;
	    }	    
	    if (m->t == VEC) {
		y = m->v.xvec;
	    } else {
		ylist = m->v.ivec;
		p->err = aggregate_discrete_check(ylist, p->dset);
	    }
	    
	    if (!p->err) {
		A = aggregate_by(x, y, xlist, ylist, fncall, 
				 p->dset, &p->err);
	    }
	}
    } else if (f == F_SUBSTR) {
	if (l->t != STR) {
	    node_type_error(f, 1, STR, l, p);
	} else if (m->t != NUM) {
	    node_type_error(f, 2, NUM, m, p);
	} else if (r->t != NUM) {
	    node_type_error(f, 3, NUM, r, p);
	} else {
	    post_process = 0;
	    ret = aux_string_node(p);

	    if (ret != NULL) {
		int ini = node_get_int(m, p);
		int fin = node_get_int(r, p);

		if (!p->err) {
		    ret->v.str = gretl_substring(l->v.str, ini, fin, &p->err);
		}
	    }
	}
    }

    if (!p->err && post_process) {
	ret = aux_matrix_node(p);
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = A;
	}
    }

    return ret;
}

static NODE *eval_print_scan (NODE *l, NODE *m, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL) {
	const char *fmt = m->v.str;
	const char *args = NULL;
	const char *lstr = NULL;
	int n = 0;

	if (l != NULL) {
	    /* note: this doesn't apply for printf */
	    if (f == F_SPRINTF) {
		if (ustring_node(l)) {
		    /* for the first argument to sprintf we need to
		       pass the name of a string variable, not its
		       value
		    */
		    lstr = l->vname;
		} else {
		    p->err = E_TYPES;
		}
	    } else {
		/* but for sscanf it's the value we want */
		lstr = l->v.str;
	    }
	}

	if (!p->err) {
	    if (!null_or_empty(r)) {
		args = r->v.str;
	    }
	    if (f == F_SSCANF) {
		p->err = do_sscanf(lstr, fmt, args, p->dset, &n);
	    } else if (f == F_SPRINTF) {
		p->err = do_printf(lstr, fmt, args, p->dset, p->prn, &n);
	    } else {
		p->err = do_printf(NULL, fmt, args, p->dset, p->prn, &n);
	    }
	    if (!p->err) {
		ret->v.xval = n;
	    }
	}
    }

    return ret;
}

static int x_to_period (double x, char c, int *err)
{
    if (na(x)) {
	return -1;
    } else if (x < 0 || fabs(x) > INT_MAX) {
	*err = E_INVARG;
	return -1;
    } else {
	int k = x;
	int ret = x;

	if (c == 'y' && k < 0) {
	    ret = -1;
	} else if (c == 'm' && (k < 1 || k > 12)) {
	    ret = -1;
	} else if (c == 'd' && (k < 1 || k > 31)) {
	    ret = -1;
	}

	if (ret < 0) {
	    fprintf(stderr, "epochday: got %c = %d!\n", c, k);
	    *err = E_INVARG;
	}

	return ret;
    }
}

static NODE *eval_epochday (NODE *ny, NODE *nm, NODE *nd, parser *p)
{
    NODE *ret = NULL;
    NODE *nodes[3] = {ny, nm, nd};
    double *x[3] = {NULL, NULL, NULL};
    int ymd[3] = {-1, -1, -1};
    const char *code = "ymd";
    double sval;
    int i;

    /* Policy: NA for y, m, or d will give an NA result;
       non-NA but out-of-bounds values for y, m or d produce
       an error.  
    */

    for (i=0; i<3 && !p->err; i++) {
	if (scalar_node(nodes[i])) {
	    sval = node_get_scalar(nodes[i], p);
	    if (!p->err) {
		ymd[i] = x_to_period(sval, code[i], &p->err);
	    }
	} else if (nodes[i]->t == VEC) {
	    x[i] = nodes[i]->v.xvec;
	} else {
	    node_type_error(F_EPOCHDAY, i+1, NUM, nodes[i], p);
	}
    }

    if (!p->err) {
	int y = ymd[0];
	int m = ymd[1];
	int d = ymd[2];

	/* From this point, -1 for y, m or d indicates either an NA
	   or the fact that a series and not a scalar was given for
	   that date-element.
	*/

	if (x[0] == NULL && x[1] == NULL && x[2] == NULL) {
	    /* y, m and d all given as scalars */
	    ret = aux_scalar_node(p);
	    if (!p->err) {
		if (y < 0 || m < 0 || d < 0) {
		    ret->v.xval = NADBL;
		} else {
		    ret->v.xval = epoch_day_from_ymd(y, m, d);
		    if (ret->v.xval < 0) {
			p->err = E_INVARG;
		    }
		}
	    }
	} else {
	    ret = aux_vec_node(p, p->dset->n);
	    if (!p->err) {
		int t;

		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    y = (x[0] == NULL)? y : x_to_period(x[0][t], 'y', &p->err);
		    m = (x[1] == NULL)? m : x_to_period(x[1][t], 'm', &p->err);
		    d = (x[2] == NULL)? d : x_to_period(x[2][t], 'd', &p->err);
		    if (p->err) {
			break;
		    } else if (y < 0 || m < 0 || d < 0) {
			/* got an NA somewhere */
			ret->v.xvec[t] = NADBL;
		    } else {
			ret->v.xvec[t] = epoch_day_from_ymd(y, m, d);
			if (ret->v.xvec[t] < 0) {
			    p->err = E_INVARG;
			}
		    }		    
		}
	    }
	}
    }

    return ret;
}

/* Bessel function handler: the 'r' node can be of scalar, series or
   matrix type.  Right now, this only supports scalar order ('m'
   node).
*/

static NODE *eval_bessel_func (NODE *l, NODE *m, NODE *r, parser *p) 
{
    char ftype;
    double v;
    NODE *ret = NULL;

    if (!starting(p) && r->t != VEC) {
	return aux_any_node(p);
    }

    ftype = l->v.str[0];
    v = node_get_scalar(m, p);

    if (r->t == NUM) {
	ret = aux_scalar_node(p);
    	if (ret != NULL) {
	    double x = r->v.xval;

	    ret->v.xval = gretl_bessel(ftype, v, x, &p->err);
    	}
    } else if (r->t == MAT) {
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    const gretl_matrix *x = r->v.m;
	    int i, n = x->rows * x->cols;

	    if (node_allocate_matrix(ret, x->rows, x->cols, p)) {
		return NULL;
	    }

	    for (i=0; i<n && !p->err; i++) {
		ret->v.m->val[i] = gretl_bessel(ftype, v, x->val[i], &p->err);
	    }
	}
    } else if (r->t == VEC) {
	ret = aux_vec_node(p, p->dset->n);
	if (ret != NULL) {
	    const double *x = r->v.xvec;
	    int t1 = (autoreg(p))? p->obs : p->dset->t1;
	    int t2 = (autoreg(p))? p->obs : p->dset->t2;
	    int t;

	    for (t=t1; t<=t2 && !p->err; t++) {
		ret->v.xvec[t] = gretl_bessel(ftype, v, x[t], &p->err);
	    }
	}
    }

    return ret;
}

/* Given an original value @x, see if it matches any of the @n0 values
   in @x0.  If so, return the substitute value from @x1, otherwise
   return the original.
*/

static double subst_val (double x, const double *x0, int n0,
			 const double *x1, int n1)
{
    int i;

    for (i=0; i<n0; i++) {
	if (x == x0[i]) {
	    return (n1 == 1)? *x1 : x1[i];
	}
    }

    return x;
}

/* String search and replace: return a node containing a copy
   of the string on node @src in which all occurrences of
   the string on @n0 are replaced by the string on @n1.
   This is literal string replacement if @f is F_STRSUB,
   regular expression replacement if @f is F_REGSUB.
*/

static NODE *string_replace (NODE *src, NODE *n0, NODE *n1, int f,
			     parser *p)
{
    if (!starting(p)) {
	return aux_string_node(p);
    } else {
	NODE *ret = NULL;
	NODE *n[3] = {src, n0, n1};
	char const *S[3];
	int i;
	
	for (i=0; i<3; i++) {
	    /* all nodes must be of string type */
	    if (n[i]->t != STR) {
		node_type_error(f, i, STR, n[i], p);
		return NULL;
	    } else {
		S[i] = n[i]->v.str;
	    }
	}

	ret = aux_string_node(p);
	if (p->err) {
	    return NULL;
	}

	if (f == F_REGSUB) {
	    ret->v.str = gretl_regexp_replace(S[0], S[1], S[2], &p->err);
	} else {
	    ret->v.str = gretl_literal_replace(S[0], S[1], S[2], &p->err);
	}
	
	return ret;
    }
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
	if (src->t == VEC) {
	    ret = aux_vec_node(p, p->dset->n);
	} else if (src->t == MAT) {
	    ret = aux_matrix_node(p);
	} else {
	    node_type_error(F_REPLACE, 3, VEC, src, p);
	}
    }

    if (!p->err) {
	double *px0 = (vx0 != NULL)? vx0->val : &x0;
	double *px1 = (vx1 != NULL)? vx1->val : &x1;
	double xt;
	int t;

	if (k0 < 0) k0 = 1;
	if (k1 < 0) k1 = 1;

	if (src->t == VEC) {
	    for (t=p->dset->t1; t<=p->dset->t2; t++) {
		xt = src->v.xvec[t];
		ret->v.xvec[t] = subst_val(xt, px0, k0, px1, k1);
	    }
	} else if (src->t == MAT) {
	    gretl_matrix *m = src->v.m;
	    int n = gretl_matrix_rows(m) * gretl_matrix_cols(m);

	    ret->v.m = gretl_matrix_copy(m);
	    if (ret->v.m == NULL) {
		p->err = E_ALLOC;
	    } else {
		for (t=0; t<n; t++) {
		    xt = m->val[t];
		    ret->v.m->val[t] = subst_val(xt, px0, k0, px1, k1);
		}
	    }
	}
    }

    return ret;
}

static void n_args_error (int k, int n, int f, parser *p)
{
    gretl_errmsg_sprintf( _("Number of arguments (%d) does not "
			    "match the number of\nparameters for "
			    "function %s (%d)"), k, getsymb(f, p), n);
    p->err = 1;
}

static NODE *isoconv_node (NODE *t, parser *p)
{
    NODE *e, *n = t->v.b1.b;
    NODE *ret = NULL;
    const double *x = NULL;
    double *ymd[3] = {NULL, NULL, NULL};
    int i, k = n->v.bn.n_nodes;

    if (p->dset == NULL) {
	p->err = E_NODATA;
	return NULL;
    }

    if (k < 3 || k > 4) {
	n_args_error(k, 4, t->t, p);
    } else {
	/* evaluate the first (series) argument */
	e = eval(n->v.bn.n[0], p);
	if (!p->err && e->t != VEC) {
	    node_type_error(t->t, 1, VEC, e, p);
	} else {
	    x = e->v.xvec + p->dset->t1;
	}
    }

    for (i=1; i<k && !p->err; i++) {
	/* the remaining args must be addresses of series */
	e = n->v.bn.n[i];
	if (i == 3 && null_or_empty(e)) {
	    ; /* OK for the last one to be omitted */
	} else if (e->t != U_ADDR) {
	    node_type_error(t->t, i+1, U_ADDR, e, p);
	} else {
	    e = e->v.b1.b;
	    if (e->t != VEC) {
		node_type_error(t->t, i+1, VEC, e, p);
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

/* evaluate a built-in function that has more than three arguments */

static NODE *eval_nargs_func (NODE *t, parser *p)
{
    NODE *e, *n = t->v.b1.b;
    NODE *ret = NULL;
    int i, k = n->v.bn.n_nodes;

    if (t->t == F_BKFILT) {
	const double *x = NULL;
	int bk[3] = {0};

	if (k < 1 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	} 

	/* evaluate the first (series) argument */
	e = eval(n->v.bn.n[0], p);
	if (!p->err && e->t != VEC) {
	    node_type_error(t->t, 1, VEC, e, p);
	}

	if (!p->err) {
	    x = e->v.xvec;
	}

	for (i=1; i<k && !p->err; i++) {
	    e = n->v.bn.n[i];
	    if (e->t == EMPTY) {
		; /* NULL arguments are OK */
	    } else {
		e = eval(n->v.bn.n[i], p);
		if (e == NULL) {
		    fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
		} else {
		    bk[i] = node_get_int(e, p);
		}
	    }
	}

	if (!p->err) {
	    ret = aux_vec_node(p, p->dset->n);
	}

	if (!p->err) {
	    p->err = bkbp_filter(x, ret->v.xvec, p->dset, bk[0], bk[1], bk[2]);
	} 
    } else if (t->t == F_FILTER) {
	const double *x = NULL;
	gretl_matrix *X = NULL;
	gretl_matrix *C = NULL;
	gretl_matrix *A = NULL;
	int freeA = 0, freeC = 0;
	double y0 = 0;

	if (k < 2 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		/* the series or matrix to filter */
		if (e->t != VEC && e->t != MAT) {
		   node_type_error(t->t, i+1, 0, e, p);
		} else if (e->t == VEC) {
		   x = e->v.xvec;
		} else {
		   X = e->v.m;
		} 
	    } else if (i == 1) {
		/* matrix for MA polynomial (but we'll take a scalar) */
		if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (e->t == MAT) {
		    C = e->v.m;
		} else if (e->t == NUM) {
		    C = gretl_matrix_from_scalar(e->v.xval);
		    if (C == NULL) {
			p->err = E_ALLOC;
		    } else {
			freeC = 1;
		    }		    
		}
	    } else if (i == 2) {
		/* matrix for AR polynomial (but we'll take a scalar) */
		if (e->t != MAT && e->t != NUM && e->t != EMPTY) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (e->t == MAT) {
		    A = e->v.m;
		} else if (e->t == NUM) {
		    A = gretl_matrix_from_scalar(e->v.xval);
		    if (A == NULL) {
			p->err = E_ALLOC;
		    } else {
			freeA = 1;
		    }
		}
	    } else if (i == 3) {
		/* initial (scalar) value for output series */
		if (!scalar_node(e)) {
		    node_type_error(t->t, i+1, NUM, e, p);
		} else {
		    y0 = node_get_scalar(e, p);
		    if (!p->err && na(y0)) {
			p->err = E_MISSDATA;
		    }
		} 
	    }
	} 
	
	if (!p->err) {
	    if (X != NULL) {
		/* matrix output wanted */
		ret = aux_matrix_node(p);
		if (!p->err) {
		    ret->v.m = filter_matrix(X, A, C, y0, &(p->err));
		}
	    } else if (x != NULL) {
		/* series output */
		ret = aux_vec_node(p, p->dset->n);
		if (!p->err) {
		    p->err = filter_series(x, ret->v.xvec, p->dset, A, C, y0);
		}
	    }
	}

	if (freeA) gretl_matrix_free(A);
	if (freeC) gretl_matrix_free(C);

    } else if (t->t == F_MCOVG) {
	gretl_matrix *X = NULL;
	gretl_vector *u = NULL;
	gretl_vector *w = NULL;
	int targ, maxlag = 0;

	if (k != 4) {
	    n_args_error(k, 4, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    targ = (i == 3)? NUM : MAT;
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if ((i == 1 || i == 2) && e->t == EMPTY) {
		; /* for u or w, NULL is acceptable */
	    } else if (e->t != targ) {
		node_type_error(t->t, i+1, targ, e, p);
	    } else if (i == 0) {
		X = e->v.m;
	    } else if (i == 1) {
		u = e->v.m;
	    } else if (i == 2) {
		w = e->v.m;
	    } else if (i == 3) {
		maxlag = e->v.xval;
	    }
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	    
	    ret->v.m = gretl_matrix_covariogram(X, u, w, maxlag, &p->err);
	} 
    } else if (t->t == F_KFILTER) {
	const char *E = NULL;
	const char *V = NULL;
	const char *S = NULL;
	const char *P = NULL;
	const char *G = NULL;

	if (k > 5) {
	    n_args_error(k, 5, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e->t == EMPTY) {
		; /* NULL arguments are OK */
	    } else if (e->t != U_ADDR) {
		node_type_error(t->t, i+1, U_ADDR, e, p);
	    } else if (i == 0) {
		E = ptr_node_get_matrix_name(e, p);
	    } else if (i == 1) {
		V = ptr_node_get_matrix_name(e, p);
	    } else if (i == 2) {
		S = ptr_node_get_matrix_name(e, p);
	    } else if (i == 3) {
		P = ptr_node_get_matrix_name(e, p);
	    } else if (i == 4) {
		G = ptr_node_get_matrix_name(e, p);
	    }
	}

	if (!p->err) {
	    ret = aux_scalar_node(p);
	}

	if (!p->err) {
	    ret->v.xval = user_kalman_run(E, V, S, P, G, 
					  p->dset, p->prn, 
					  &p->err);
	} 
    } else if (t->t == F_KSMOOTH) {
	const char *P = NULL;
	const char *U = NULL;

	if (k > 2) {
	    n_args_error(k, 2, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e->t == EMPTY) {
		; /* NULL arguments are acceptable */
	    } else if (e->t != U_ADDR) {
		node_type_error(t->t, i+1, U_ADDR, e, p);
	    } else if (i == 0) {
		P = ptr_node_get_matrix_name(e, p);
	    } else if (i == 1) {
		U = ptr_node_get_matrix_name(e, p);
	    }
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	    
	    ret->v.m = user_kalman_smooth(P, U, &p->err);
	} 
    } else if (t->t == F_KSIMUL) {
	gretl_matrix *V = NULL;
	gretl_matrix *W = NULL;
	const char *S = NULL;
	int freeV = 0, freeW = 0;
	
	if (k < 1 || k > 3) {
	    n_args_error(k, 1, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i == 0) {
		if (e->t == VEC) {
		    V = tmp_matrix_from_series(e, p);
		    freeV = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    V = e->v.m;
		}
	    } else if (i == 1) {
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t == VEC) {
		    W = tmp_matrix_from_series(e, p);
		    freeW = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    W = e->v.m;
		}
	    } else {		
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else {
		    S = ptr_node_get_matrix_name(e, p);
		}
	    } 
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	    
	    ret->v.m = user_kalman_simulate(V, W, S, p->prn, &p->err);
	}

	if (freeV) gretl_matrix_free(V);
	if (freeW) gretl_matrix_free(W);
    } else if (t->t == F_MOLS || t->t == F_MPOLS) {
	gretlopt opt = (t->t == F_MPOLS)? OPT_M : OPT_NONE;
	gretl_matrix *Y = NULL;
	gretl_matrix *X = NULL;
	const char *SU = NULL;
	const char *SV = NULL;
	int freeY = 0, freeX = 0;

	if (k < 2 || k > 4) {
	    n_args_error(k, 1, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i == 0) {
		if (e->t == VEC) {
		    Y = tmp_matrix_from_series(e, p);
		    freeY = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    Y = e->v.m;
		}
	    } else if (i == 1) {
		if (e->t == VEC) {
		    X = tmp_matrix_from_series(e, p);
		    freeX = 1;
		} else if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    X = e->v.m;
		}
	    } else {		
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else if (i == 2) {
		    SU = ptr_node_get_matrix_name(e, p);
		} else {
		    SV = ptr_node_get_matrix_name(e, p);
		}
	    } 
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	
	    ret->v.m = user_matrix_ols(Y, X, SU, SV, opt, &p->err);
	}

	if (freeY) gretl_matrix_free(Y);
	if (freeX) gretl_matrix_free(X);
    } else if (t->t == F_MRLS) {
	gretl_matrix *Y = NULL;
	gretl_matrix *X = NULL;
	gretl_matrix *R = NULL;
	gretl_matrix *Q = NULL;
	const char *SU = NULL;
	const char *SV = NULL;

	if (k < 4 || k > 6) {
	    n_args_error(k, 1, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i < 4) {
		if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else if (i == 0) {
		    Y = e->v.m;
		} else if (i == 1) {
		    X = e->v.m;
		} else if (i == 2) {
		    R = e->v.m;
		} else if (i == 3) {
		    Q = e->v.m;
		}
	    } else {
		if (e->t == EMPTY) {
		    ; /* OK */
		} else if (e->t != U_ADDR) {
		    node_type_error(t->t, i+1, U_ADDR, e, p);
		} else if (i == 4) {
		    SU = ptr_node_get_matrix_name(e, p);
		} else {
		    SV = ptr_node_get_matrix_name(e, p);
		}
	    } 
	}

	if (!p->err) {
	    ret = aux_matrix_node(p);
	}

	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	
	    ret->v.m = user_matrix_rls(Y, X, R, Q, SU, SV, &p->err);
	}
    } else if (t->t == F_NRMAX) {
	gretl_matrix *b = NULL;
	const char *sf = NULL;
	const char *sg = NULL;
	const char *sh = NULL;

	if (k < 2 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i == 0) {
		if (e->t != MAT) {
		    node_type_error(t->t, i+1, MAT, e, p);
		} else {
		    b = e->v.m;
		}
	    } else if (i == 1) {
		if (e->t != STR) {
		    node_type_error(t->t, i+1, STR, e, p);
		} else {
		    sf = e->v.str;
		}
	    } else if (e->t == EMPTY) {
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
	    ret->v.xval = user_NR(b, sf, sg, sh, p->dset, 
				  p->prn, &p->err);
	}
    } else if (t->t == F_LOESS) {
	const double *y = NULL, *x = NULL;
	double bandwidth = 0.5;
	int poly_order = 1;
	gretlopt opt = OPT_NONE;

	if (k < 2 || k > 6) {
	    n_args_error(k, 5, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (p->err) {
		break;
	    }
	    if (i < 2) {
		if (e->t != VEC) {
		    node_type_error(t->t, i+1, VEC, e, p);
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
	    ret = aux_vec_node(p, p->dset->n);
	    if (ret != NULL) {
		p->err = gretl_loess(y, x, poly_order, bandwidth,
				     opt, p->dset, ret->v.xvec);
	    }
	}
    } else if (t->t == F_GHK) {
	gretl_matrix *C = NULL;
	gretl_matrix *A = NULL;
	gretl_matrix *B = NULL;
	gretl_matrix *U = NULL;

	if (k != 4) {
	    n_args_error(k, 4, t->t, p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (e->t != MAT) {
		node_type_error(t->t, i+1, MAT, e, p);
	    } else if (i == 0) {
		C = e->v.m;
	    } else if (i == 1) {
		A = e->v.m;
	    } else if (i == 2) {
		B = e->v.m;
	    } else if (i == 3) {
		U = e->v.m;
	    }
	}
	if (!p->err) {
	    ret = aux_matrix_node(p);
	}
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	    
	    ret->v.m = gretl_GHK(C, A, B, U, &p->err);
	} 
    } else if (t->t == F_QUADTAB) {
	int order = -1, method = 0;
	double a = NADBL;
	double b = NADBL;

	if (k < 1 || k > 4) {
	    n_args_error(k, 4, t->t, p);
	} 
	
	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
	    } else if (i == 0) {
		order = node_get_int(e, p);
	    } else if (!empty_or_num(e)) {
		node_type_error(t->t, i+1, NUM, e, p);
	    } else if (i == 1) {
		method = node_get_int(e, p);
	    } else if (i == 2) {
		a = e->v.xval;
	    } else {
		b = e->v.xval;
	    }
	}
	if (!p->err) {
	    ret = aux_matrix_node(p);
	}
	if (!p->err) {
	    if (ret->v.m != NULL) {
		gretl_matrix_free(ret->v.m);
	    }	    
	    ret->v.m = gretl_quadrule_matrix_new(order, method, 
						 a, b, &p->err);
	} 
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

static gretl_matrix *assemble_matrix (NODE *nn, int nnodes, parser *p)
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

    /* how many columns will we need? */
    for (i=0; i<nnodes; i++) {
	n = nn->v.bn.n[i];
	if (n->t == LIST) {
	    k += n->v.ivec[0];
	} else if (n->t == VEC) {
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
	n = nn->v.bn.n[i];
	if (n->t == LIST) {
	    list = n->v.ivec;
	    for (j=1; j<=list[0]; j++) {
		Z[k++] = p->dset->Z[list[j]];
	    }	    
	} else if (n->t == VEC) {
	    if (useries_node(n) && compiled(p)) {
		reattach_series(n, p);
	    } 
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

static NODE **tmp_node_holder (NODE *n, parser *p)
{
    int i, m = n->v.bn.n_nodes;
    NODE **t = malloc(m * sizeof *t);

    if (t == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<m; i++) {
	t[i] = n->v.bn.n[i];
    }

    return t;
}

#define ok_matdef_sym(s) (s == NUM || s == VEC || s == EMPTY || \
                          s == DUM || s == LIST)

/* composing a matrix from scalars, series or lists */

static NODE *matrix_def_node (NODE *nn, parser *p)
{
    gretl_matrix *M = NULL;
    NODE *n, *ret = NULL;
    NODE **nntmp = NULL;
    int m = nn->v.bn.n_nodes;
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

    if (reusable(p)) {
	nntmp = tmp_node_holder(nn, p);
	if (p->err) {
	    return NULL;
	}
    }

#if EDEBUG
    fprintf(stderr, "Processing MDEF...\n");
#endif

    for (i=0; i<m && !p->err; i++) {
	NODE *ni = nn->v.bn.n[i];

	if (!compiled(p) && (ok_matdef_sym(ni->t) || scalar_matrix_node(ni))) {
	    /* pass through */
	    n = ni;
	} else {
	    n = eval(ni, p);
	    if (n == NULL && !p->err) {
		p->err = E_UNSPEC; /* "can't happen" */
	    }
	    if (p->err) {
		break;
	    }
	    if (!ok_matdef_sym(n->t) && !scalar_matrix_node(n)) {
		fprintf(stderr, "matrix_def_node: node type %d: not OK\n", n->t);
		p->err = E_TYPES;
		break;
	    } else if (n != ni) {
		if (nntmp == NULL) {
		    free_tree(nn->v.bn.n[i], p, "MatDef");
		}
		nn->v.bn.n[i] = n;
	    }
	}
	if (scalar_node(n)) {
	    nnum++;
	} else if (n->t == VEC) {
	    nvec++;
	} else if (n->t == DUM) {
	    dum++;
	} else if (n->t == LIST) {
	    nlist++;
	} else if (n->t == EMPTY) {
	    if (nsep == 0) {
		seppos = i;
	    }
	    nsep++;
	}

	if (dum && m != 1) {
	    /* dummy must be singleton node */
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
	if (nvec > 0 || nlist > 1) {
	    M = assemble_matrix(nn, m, p);
	} else if (nnum > 0) {
	    M = matrix_from_scalars(nn, m, nsep, seppos, p);
	} else if (nlist) {
	    M = matrix_from_list(nn->v.bn.n[0], p);
	} else if (dum) {
	    n = nn->v.bn.n[0];
	    if (n->v.idnum == DUM_DATASET) {
		M = matrix_from_list(NULL, p);
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

    if (nntmp != NULL) {
	/* restore the original subnodes */
	for (i=0; i<m; i++) {
	    nn->v.bn.n[i] = nntmp[i];
	}
	free(nntmp);
    } else {
	/* forestall double-freeing: null out any aux nodes */
	for (i=0; i<m; i++) {
	    if (is_aux_node(nn->v.bn.n[i])) {
		nn->v.bn.n[i] = NULL;
	    }
	}
    }
	
    return ret;
}

static NODE *gen_series_node (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (l->t != STR || r->t != STR) {
	p->err = E_TYPES;
    } else if (p->dset == NULL || p->dset->n == 0) {
	no_data_error(p);
    } else {
	const char *vname = l->v.str;
	const char *formula = r->v.str;
	char *line;

	line = malloc(9 + strlen(vname) + strlen(formula));
	if (line == NULL) {
	    p->err = E_ALLOC;
	} else {
	    int vnum = -1;
	    int err = 0;

	    sprintf(line, "series %s=%s", vname, formula);
	    err = generate(line, p->dset, OPT_NONE, p->prn);
	    if (!err) {
		vnum = current_series_index(p->dset, vname);
	    }
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		ret->v.xval = vnum;
	    }
	    free(line);
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
   If so, return FORK_L (all 1) or FORK_R (all 0), othewise
   return FORK_UNK.
*/

static int vec_branch (const double *c, parser *p)
{
    int c1, t, t1, t2;
    int ret;

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    c1 = (c[t1] != 0.0); 
    ret = (c1)? FORK_L : FORK_R;

    for (t=t1; t<=t2; t++) {
	if (!xna(c[t])) {
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
   scalar (NUM) and series (VEC) types on input, and always produce
   a VEC type on output.
*/

static NODE *query_eval_vec (const double *c, NODE *n, parser *p)
{
    NODE *l = NULL, *r = NULL, *ret = NULL;
    double *xvec = NULL, *yvec = NULL;
    double x = NADBL, y = NADBL;
    double xt, yt;
    int t, t1, t2;
    int branch;

    branch = vec_branch(c, p);

    if (autoreg(p) || branch != FORK_R) {
	l = eval(n->v.b3.m, p);
	if (p->err) {
	    return NULL;
	}
	if (l->t == VEC) {
	    xvec = l->v.xvec;
	} else if (l->t == NUM) {
	    x = l->v.xval;
	} else {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    if (autoreg(p) || branch != FORK_L) {
	r = eval(n->v.b3.r, p);
	if (p->err) {
	    return NULL;
	}
	if (r->t == VEC) {
	    yvec = r->v.xvec;
	} else if (r->t == NUM) {
	    y = r->v.xval;
	} else {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    ret = aux_vec_node(p, p->dset->n);

    t1 = (autoreg(p))? p->obs : p->dset->t1;
    t2 = (autoreg(p))? p->obs : p->dset->t2;

    for (t=t1; t<=t2; t++) {
	if (xna(c[t])) {
	    ret->v.xvec[t] = NADBL;
	} else {
	    xt = (xvec != NULL)? xvec[t] : x;
	    yt = (yvec != NULL)? yvec[t] : y;
	    ret->v.xvec[t] = (c[t] != 0.0)? xt : yt;
	}
    }

    return ret;
}

static NODE *query_eval_scalar (double x, NODE *n, parser *p)
{
    NODE *l = NULL, *r = NULL, *ret = NULL;
    int indef = na(x) || isnan(x);
    int branch;

    branch = indef ? FORK_NONE : (x != 0 ? FORK_L : FORK_R);

    if (autoreg(p) || branch != FORK_R) {
	l = eval(n->v.b3.m, p);
	if (p->err) {
	    return NULL;
	}
    }

    if (autoreg(p) || branch != FORK_L) {
	r = eval(n->v.b3.r, p);
	if (p->err) {
	    return NULL;
	}
    }

    if (branch == FORK_NONE) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = NADBL;
	}
    } else if (branch == FORK_L) {
	ret = l;
    } else if (branch == FORK_R) {
	ret = r;
    }

    return ret;
}

static NODE *query_eval_matrix (gretl_matrix *m, NODE *n, parser *p)
{
    NODE *ret, *l, *r;
    gretl_matrix *mret;

    if (gretl_is_null_matrix(m)) {
	p->err = E_TYPES;
	return NULL;
    }

    l = eval(n->v.b3.m, p);
    if (p->err) {
	return NULL;
    }

    r = eval(n->v.b3.r, p);
    if (p->err) {
	return NULL;
    }

    if ((l->t != NUM && l->t != MAT) ||
	(r->t != NUM && r->t != MAT)) {
	p->err = E_TYPES;
	return NULL;
    }

    if (l->t == MAT && (l->v.m->cols != m->cols ||
			l->v.m->rows != m->rows)) {
	p->err = E_NONCONF;
	return NULL;
    } else if (r->t == MAT && (r->v.m->cols != m->cols ||
			       r->v.m->rows != m->rows)) {
	p->err = E_NONCONF;
	return NULL;
    }	

    mret = gretl_matrix_copy(m);
    if (mret == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    ret = aux_matrix_node(p);

    if (!p->err) {
	double x, y;
	int j, i;

	for (j=0; j<m->cols; j++) {
	    for (i=0; i<m->rows; i++) {
		x = gretl_matrix_get(m, i, j);
		if (isnan(x)) {
		    gretl_matrix_set(mret, i, j, x);
		} else if (x != 0.0) {
		    y = l->t == NUM ? l->v.xval :
			gretl_matrix_get(l->v.m, i, j);
		    gretl_matrix_set(mret, i, j, y);
		} else {
		    y = r->t == NUM ? r->v.xval :
			gretl_matrix_get(r->v.m, i, j);
		    gretl_matrix_set(mret, i, j, y);
		}
	    }
	}
	ret->v.m = mret;
    }

    return ret;
}

/* Handle the case where a ternary "query" expression has produced one
   of its own child nodes as output: we duplicate the information in an
   auxiliary node so as to avoid double-freeing of the result.
*/

static NODE *ternary_return_node (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (n->t == NUM) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = n->v.xval;
	}
    } else if (n->t == VEC) {
	int t, T = p->dset->n;

	ret = aux_vec_node(p, T);
	if (ret != NULL) {
	    for (t=0; t<T; t++) {
		ret->v.xvec[t] = n->v.xvec[t];
	    }
	}
    } else if (n->t == MAT) {
	ret = aux_matrix_node(p);
	if (ret != NULL) {
	    if (is_tmp_node(ret)) {
		gretl_matrix_free(ret->v.m);
	    }
	    ret->v.m = gretl_matrix_copy(n->v.m);
	    if (ret->v.m == NULL) {
		p->err = E_ALLOC;
	    }
	}
    } else if (n->t == STR) {
	ret = aux_string_node(p);
	if (ret != NULL) {
	    if (is_tmp_node(ret)) {
		free(ret->v.str);
	    }
	    ret->v.str = gretl_strdup(n->v.str);
	    if (ret->v.str == NULL) {
		p->err = E_ALLOC;
	    }	    
	}
    } else if (n->t == LIST) {
	ret = aux_list_node(p);
	if (ret != NULL) {
	    if (is_tmp_node(ret)) {
		free(ret->v.ivec);
	    }
	    ret->v.ivec = gretl_list_copy(n->v.ivec);
	    if (ret->v.ivec == NULL) {
		p->err = E_ALLOC;
	    }	    
	}
    } else if (n->t == EMPTY) {
	ret = aux_empty_node(p);
	if (ret == NULL) {
	    p->err = E_ALLOC;
	}	
    } else {
	fprintf(stderr, "ternary_return_node: unhandled type '%s'\n", 
		getsymb(n->t, p));
	p->err = E_TYPES;
    }

    return ret;
}

/* Evaluate a ternary "query" expression: (C)? X : Y.  The condition C
   must be a scalar, series or matrix.  The relevant sub-nodes of @t
   are named "l" (left, the condition), "m" and "r" (middle and right
   respectively, the two alternates).
*/

static NODE *eval_query (NODE *t, parser *p)
{
    NODE *e, *ret = NULL;
    double *vec = NULL;
    gretl_matrix *m = NULL;
    double x = NADBL;

#if EDEBUG
    fprintf(stderr, "eval_query: t=%p, l=%p, m=%p, r=%p\n", 
	    (void *) t, (void *) t->v.b3.l, (void *) t->v.b3.m,
	    (void *) t->v.b3.r);
#endif

    /* evaluate and check the condition */

    e = eval(t->v.b3.l, p);

    if (!p->err) {
	if (e->t == NUM) {
	    x = e->v.xval;
	} else if (e->t == VEC) {
	    vec = e->v.xvec;
	} else if (e->t == MAT) {
	    m = e->v.m;
	} else {
	    p->err = E_TYPES;
	}
    }

    if (p->err) {
	return NULL;
    }

    if (vec != NULL) {
	ret = query_eval_vec(vec, t, p);
    } else if (m != NULL) {
	ret = query_eval_matrix(m, t, p);
    } else {
	ret = query_eval_scalar(x, t, p);
    }

    if (ret != NULL && (ret == t->v.b3.m || ret == t->v.b3.r)) {
	/* forestall double-freeing */
	ret = ternary_return_node(ret, p);
    }

    return ret;
}

static int get_version_as_scalar (void)
{
    int x, y, z;

    sscanf(GRETL_VERSION, "%d.%d.%d", &x, &y, &z);
    return 10000 * x + 100 * y + z;
}

#define dvar_scalar(i) (i > 0 && i < R_SCALAR_MAX)
#define dvar_series(i) (i > R_SCALAR_MAX && i < R_SERIES_MAX)
#define dvar_variant(i) (i > R_SERIES_MAX && i < R_MAX)

#define no_data(p) (p == NULL || p->n == 0)

double dvar_get_scalar (int i, const DATASET *dset,
			char *label)
{
    int ival;

    switch (i) {
    case R_NOBS:
	return (dset == NULL) ? NADBL : 
	(dset->n == 0 ? 0 : sample_size(dset));
    case R_NVARS:
	return (dset == NULL)? NADBL : dset->v;
    case R_PD:
	return (no_data(dset))? NADBL : dset->pd;
    case R_T1:
	return (no_data(dset))? NADBL : dset->t1 + 1;
    case R_T2:
	return (no_data(dset))? NADBL : dset->t2 + 1;
    case R_DATATYPE:
	return (no_data(dset))? NADBL : dataset_get_structure(dset);
    case R_TEST_PVAL:
	return get_last_pvalue(label);
    case R_TEST_STAT:
	return get_last_test_statistic(label);
    case R_TEST_LNL:
	return get_last_lnl(label);
    case R_KLNL:
	return user_kalman_get_loglik();
    case R_KS2:
	return user_kalman_get_s2();
    case R_KSTEP:
	ival = user_kalman_get_time_step();

	return (double) ival;
    case R_STOPWATCH:
	return gretl_stopwatch();
    case R_WINDOWS:
#ifdef WIN32
	return 1;
#else
	return 0;
#endif
    case R_VERSION:
	return get_version_as_scalar();
    case R_ERRNO:
	return get_gretl_errno();
    case R_SEED:
	return gretl_rand_get_seed();
    case R_HUGE:
	return libset_get_double(CONV_HUGE);
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

static double *dvar_get_series (int i, const DATASET *dset, 
				int *err)
{
    double *x = NULL;
    int YMD = calendar_data(dset);
    int t;

    if (dset == NULL || dset->n == 0) {
	*err = E_NODATA;
	return NULL;
    }

    if (i == R_OBSMIN && dset->pd < 2) {
	*err = E_PDWRONG;
	return NULL;
    }    

    if (i == R_OBSMIC && !YMD) {
	*err = E_PDWRONG;
	return NULL;
    }

    if (i == R_PUNIT && !dataset_is_panel(dset)) {
	*err = E_PDWRONG;
	return NULL;
    }

    if (i == R_DATES && !date_series_ok(dset)) {
	*err = E_PDWRONG;
	return NULL;
    }

    if (i == R_OBSMAJ) {
	if (dset->pd == 1 && !dataset_is_time_series(dset)) {
	    i = R_INDEX;
	} else if (dataset_is_panel(dset)) {
	    i = R_PUNIT;
	}
    }

    x = malloc(dset->n * sizeof *x);
    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
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

	for (t=0; t<dset->n && !*err; t++) {
	    ntodate(obs, t, dset);
	    if (sscanf(obs, YMD_READ_FMT, &y, &m, &d) != 3) {
		*err = E_DATA;
	    } else if (i ==  R_OBSMAJ) {
		x[t] = y;
	    } else if (i == R_OBSMIN) {
		x[t] = m;
	    } else {
		x[t] = d;
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
	*err = fill_dataset_dates_series(dset, x);
	if (*err) {
	    free(x);
	    x = NULL;
	}
    } else {
	*err = E_DATA;
    }

    return x;
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
    default:
	*err = E_DATA;
	break;
    }

    return m;
}

static gretl_matrix *submat_postprocess (gretl_matrix *M,
					 NODE *n, parser *p)
{
    gretl_matrix *S = matrix_get_submatrix(M, n->v.mspec, 0, &p->err);

    gretl_matrix_free(M);

    return S;
}

static gretl_matrix *dvar_get_submatrix (int idx, NODE *t, parser *p)
{
    NODE *r = eval(t->v.b2.r, p);
    gretl_matrix *M, *S = NULL;

    if (r == NULL || r->t != MSPEC) {
	if (!p->err) {
	    p->err = E_TYPES;
	}
	return NULL;
    }

    M = dvar_get_matrix(idx, &p->err);

    if (M != NULL) {
	S = submat_postprocess(M, r, p);
    }

    return S;
}

static NODE *dollar_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	int idx, mslice = (t->t == DMSL);

	if (mslice) {
	    /* "slice" spec on right subnode, index on left */
	    idx = t->v.b2.l->v.idnum;
	} else {
	    idx = t->v.idnum;
	}

	if (mslice) {
	    ret = aux_matrix_node(p);
	    if (ret != NULL) {
		ret->v.m = dvar_get_submatrix(idx, t, p);
	    }
	} else if (dvar_scalar(idx)) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL) {
		ret->v.xval = dvar_get_scalar(idx, p->dset, 
					      p->lh.label);
	    }
	} else if (dvar_series(idx)) {
	    ret = aux_vec_node(p, 0);
	    if (ret != NULL) {
		ret->v.xvec = dvar_get_series(idx, p->dset,
					      &p->err);
	    }
	} else if (dvar_variant(idx)) {
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
		    ret->v.xval = dvar_get_scalar(idx, p->dset, 
						  p->lh.label);
		}
	    } 
	}	    
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

/* The incoming node here is binary: matrix ID on the left, and
   specification of a matrix sub-slice on the right.
*/

static gretl_matrix *
object_var_get_submatrix (const char *oname, int idx, NODE *t, parser *p,
			  int needs_data)
{
    NODE *r = eval(t->v.b2.r, p);
    gretl_matrix *M, *S = NULL;

    if (r == NULL || r->t != MSPEC) {
	if (!p->err) {
	    p->err = E_TYPES;
	}
	return NULL;
    }

    if (needs_data) {
	M = saved_object_build_matrix(oname, idx, p->dset, &p->err);
    } else {
	M = saved_object_get_matrix(oname, idx, &p->err);
    }

    if (M != NULL) {
	S = submat_postprocess(M, r, p);
    }

    return S;
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
    }

    if (idx == M_UHAT || idx == M_YHAT || idx == M_SIGMA) {
	/* could be a matrix */
	int ci = 0;
	GretlObjType otype = gretl_model_get_type_and_ci(oname, &ci);

	if (otype != GRETL_OBJ_EQN) {
	    vtype = GRETL_TYPE_MATRIX;
	} else if ((idx == M_UHAT || idx == M_YHAT) && ci == BIPROBIT) {
	    vtype = GRETL_TYPE_MATRIX;
	}
    } 

    return vtype;
}

/* For example, $coeff(sqft); or model1.$coeff(const) */

static NODE *dollar_str_node (NODE *t, MODEL *pmod, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	NODE *l = t->v.b2.l;
	NODE *r = t->v.b2.r;

	ret->v.xval = gretl_model_get_data_element(pmod, l->v.idnum, r->v.str, 
						   p->dset, &p->err);

	if (na(ret->v.xval)) {
	    const char *s = get_string_by_name(r->v.str);

	    if (s != NULL) {
		p->err = 0;
		ret->v.xval = gretl_model_get_data_element(pmod, l->v.idnum, s, 
							   p->dset, &p->err);
	    }
	}

	if (na(ret->v.xval)) {
	    p->err = E_INVARG;
	    pprintf(p->prn, _("'%s': invalid argument for %s()\n"), 
		    r->v.str, mvarname(l->v.idnum));
	}
    }

    return ret;
}

/* Retrieve a data item from an object (typically, a model). We're not
   sure in advance here of the type of the data item (scalar, matrix,
   etc.): we look that up with object_var_type().

   This function handles two cases: the input @t may be a binary node
   with the name of an object on its left-hand subnode (as in
   "mymodel.$vcv"), or the object may be implicit (as when accessing
   data from the last model). In the latter case @t is itself the data
   item specification, while in the former the data item spec will be
   found on the right-hand subnode of @t. 

   And there's another thing: in the case where the data item to be
   retrieved is a matrix, we handle the possibility that the user
   actually wants a sub-slice of that matrix. Handling that in
   this function may be the wrong thing to do -- perhaps it should be
   subsumed under the more general handling of <matrix>[<subspec>].
*/

static NODE *object_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

#if EDEBUG
    fprintf(stderr, "object_var_node: t->t = %d\n", t->t);
#endif

    if (starting(p)) {
	NODE *r = NULL; /* the data spec node */
	const char *oname = NULL;
	GretlType vtype;
	int idx, mslice;
	int needs_data = 0;

	if (t->t == OVAR) {
	    /* objectname.<stuff> */
	    oname = t->v.b2.l->v.str;
	    r = t->v.b2.r;
	} else {
	    /* plain <stuff> */
	    r = t;
	}

	if (oname != NULL && gretl_get_object_by_name(oname) == NULL) {
	    gretl_errmsg_sprintf(_("%s: no such object\n"), oname);
	    p->err = E_UNKVAR;
	    return NULL;
	}

	if (oname != NULL && r->t == DMSTR) {
	    MODEL *pmod = get_model_by_name(oname);
	    
	    if (pmod == NULL) {
		p->err = E_INVARG;
		return NULL;
	    } else {
		return dollar_str_node(r, pmod, p);
	    }
	}

	mslice = (r->t == DMSL);

	/* find the index which identifies the data item
	   the user wants */

	if (mslice) {
	    /* slice spec is on right subnode, index on left */
	    idx = r->v.b2.l->v.idnum;
	} else {
	    idx = r->v.idnum;
	}

	if (mslice && dvar_variant(idx)) {
	    /* we're in the wrong place */
	    return dollar_var_node(r, p);
	}

	vtype = object_var_type(idx, oname, &needs_data);

#if EDEBUG
	fprintf(stderr, "object_var_node: t->t = %d (%s), r->t = %d (%s)\n", 
		t->t, getsymb(t->t, NULL), r->t, getsymb(r->t, NULL));
	fprintf(stderr, "idx = %d, vtype = %d, mslice = %d\n", 
		idx, vtype, mslice);
#endif

	if (vtype == GRETL_TYPE_DOUBLE) {
	    ret = aux_scalar_node(p);
	} else if (vtype == GRETL_TYPE_SERIES) {
	    ret = aux_vec_node(p, 0);
	} else if (vtype == GRETL_TYPE_LIST) {
	    ret = aux_list_node(p);
	} else if (vtype == GRETL_TYPE_STRING) {
	    ret = aux_string_node(p);
	} else {
	    ret = aux_matrix_node(p);
	}

	if (ret == NULL) {
	    return ret;
	} else if (vtype == GRETL_TYPE_DOUBLE) {
	    ret->v.xval = saved_object_get_scalar(oname, idx, p->dset, 
						  &p->err);
	} else if (vtype == GRETL_TYPE_SERIES) {
	    ret->v.xvec = saved_object_get_series(oname, idx, p->dset,
						  &p->err);
	} else if (vtype == GRETL_TYPE_LIST) {
	    ret->v.ivec = saved_object_get_list(oname, idx, &p->err);
	} else if (vtype == GRETL_TYPE_STRING) {
	    ret->v.str = saved_object_get_string(oname, idx, p->dset,
						 &p->err);
	} else if (mslice) {
	    /* the right-hand subnode needs more work */
	    ret->v.m = object_var_get_submatrix(oname, idx, r, p, needs_data);
	} else if (vtype == GRETL_TYPE_MATRIX) {
	    if (needs_data) {
		ret->v.m = saved_object_build_matrix(oname, idx,
						     p->dset, &p->err);
	    } else {
		ret->v.m = saved_object_get_matrix(oname, idx, &p->err);
	    }
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
	int *list = ellipsis_list(p->dset, l->vnum, r->vnum,
				  &p->err);

	ret->v.ivec = list;
    }

    return ret;
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
	double xl = l->v.xval;
	double xr = r->v.xval;

	if (!na(xl) && !na(xr)) {
	    if (t == F_XMIN) {
		ret->v.xval = (xl < xr)? xl : xr;
	    } else if (t == F_XMAX) {
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

static NODE *scalar_postfix_node (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	double x = n->v.xval;

	ret->v.xval = x;
	if (n->t == UNUM_P) {
	    p->err = gretl_scalar_set_value(n->vname, x + 1.0);
	} else {
	    p->err = gretl_scalar_set_value(n->vname, x - 1.0);
	}
    }

    return ret;
}

static int series_calc_nodes (NODE *l, NODE *r)
{
    int ret = 0;

    if (l->t == VEC) {
	ret = (r->t == VEC || r->t == NUM || scalar_matrix_node(r));
    } else if (r->t == VEC) {
	ret = scalar_node(l);
    }

    return ret;
}

static void reattach_series (NODE *n, parser *p)
{
    int v = n->vnum;

    if (v >= p->dset->v) {
	fprintf(stderr, "VEC node, ID = %d but p->dset->v = %d\n", v, p->dset->v);
	p->err = E_DATA;
    } else {
	n->v.xvec = p->dset->Z[v];
    }
}

static void node_reattach_data (NODE *n, parser *p)
{
    void *data = n;

    if (uscalar_node(n)) {
	n->v.xval = gretl_scalar_get_value(n->vname, &p->err);
    } else if (umatrix_node(n)) {
	data = n->v.m = get_matrix_by_name(n->vname);
    } else if (useries_node(n)) {
	reattach_series(n, p);
    } else if (ulist_node(n)) {
	data = n->v.ivec = get_list_by_name(n->vname);
    } else if (ubundle_node(n)) {
	if (n->vname[0] == '$') {
	    /* built-in bundle */
	    return;
	}
	data = n->v.b = get_bundle_by_name(n->vname);
    } else if (ustring_node(n)) {
	data = n->v.str = get_string_by_name(n->vname);
    }

    if (!p->err && data == NULL) {
	p->err = E_UNKVAR;
    }
}

static void node_type_error (int ntype, int argnum, int goodt, 
			     NODE *bad, parser *p)
{
    const char *nstr;

    if (ntype == LAG) {
	nstr = (goodt == NUM)? "lag order" : "lag variable";
    } else {
	nstr = getsymb(ntype, NULL);
    }

    if (goodt == 0 && bad != NULL) {
	if (bad->t == EMPTY) {
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
    } 
    pputc(p->prn, '\n');

    if (!strcmp(nstr, "&")) {
	pputs(p->prn, "(for logical AND, please use \"&&\")\n");
    } else if (!strcmp(nstr, "|")) {
	pputs(p->prn, "(for logical OR, please use \"||\")\n");
    }

    p->err = E_TYPES;
}

static NODE *bool_node (int s, parser *p)
{
    NODE *n = aux_scalar_node(p);

    if (n != NULL) {
	n->v.xval = s;
    }

    return n;
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

#define eval_left(t)   (evalb1(t) || evalb2(t) || evalb3(t))
#define eval_middle(t) (evalb3(t))
#define eval_right(t)  (evalb2(t) || evalb3(t))

static NODE *raw_node (NODE *t, int i)
{
    if (evalb1(t->t)) {
	return t->v.b1.b;
    } else if (evalb2(t->t)) {
	return (i==0)? t->v.b2.l : t->v.b2.r;
    } else {
	/* 3-place node */
	return (i==0)? t->v.b3.l : (i==1)? t->v.b3.m : t->v.b3.r;
    }
}

/* core function: evaluate the parsed syntax tree */

static NODE *eval (NODE *t, parser *p)
{  
    NODE *l = NULL, *m = NULL, *r = NULL;
    NODE *ret = NULL;

    if (t == NULL) {
	fprintf(stderr, "eval: got NULL input node\n");
	p->err = E_ALLOC;
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "eval: incoming node %p ('%s', vname=%s)\n", 
	    (void *) t, getsymb(t->t, p), t->vname);
#endif

    if (!p->err && eval_left(t->t)) {
	l = eval(raw_node(t, 0), p);
	if (l == NULL && !p->err) {
	    p->err = 1;
	}
    }

    if (!p->err && eval_middle(t->t)) {
	if (m_return(t->t)) {
	    m = t->v.b3.m;
	} else {
	    m = eval(t->v.b3.m, p);
	    if (m == NULL && !p->err) {
		p->err = 1;
	    }
	}
    }	    

    if (!p->err && eval_right(t->t)) {
	if (r_return(t->t)) {
	    r = raw_node(t, 2);
	} else {
	    if ((t->t == B_AND || t->t == B_OR) && l != NULL && l->t == NUM) {
		/* logical operators: avoid redundant evaluations */
		if (t->t == B_AND && node_is_false(l, p)) {
		    r = bool_node(0, p);
		} else if (t->t == B_OR && node_is_true(l, p)) {
		    r = bool_node(1, p);
		}
	    }
	    if (r == NULL && !p->err) {
		r = eval(raw_node(t, 2), p);
	    }
	    if (r == NULL && !p->err) {
		p->err = 1;
	    }
	}
    }

    if (p->err) {
	goto bailout;
    }

    switch (t->t) {
    case MSPEC:
    case EMPTY:
    case ABSENT:
	ret = t;
	break;
    case NUM:
    case VEC:
    case MAT:
    case STR:
    case LIST:
    case BUNDLE:
	if (compiled(p) && starting(p) && uvar_node(t)) {
	    node_reattach_data(t, p);
	}
	ret = t;
	break;
    case UNDEF:
	undefined_symbol_error(t->vname, p);
	break;
    case U_ADDR:
	p->err = check_uaddr_type(t, p);
	if (!p->err && compiled(p) && starting(p)) {
	    node_reattach_data(t->v.b1.b, p);
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
	} else {
	    /* otherwise treat as terminal */
	    ret = t;
	}
	break;
    case UNUM_P:
    case UNUM_M:
	if (compiled(p) && starting(p)) {
	    node_reattach_data(t, p);
	}
	ret = scalar_postfix_node(t, p);
	break;
    case FARGS:
	/* will be evaluated in context */
	ret = t;
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
	} else if ((t->t == B_EQ || t->t == B_NEQ) && l->t == STR && r->t == STR) {
	    ret = compare_strings(l, r, t->t, p);
	} else if (l->t == NUM && r->t == NUM) {
	    ret = scalar_calc(l, r, t->t, p);
	} else if (l->t == BUNDLE && r->t == BUNDLE) {
	    ret = bundle_op(l, r, t->t, p);
	} else if (stringvec_node(l) && stringvec_node(r)) {
	    ret = stringvec_calc(l, r, t->t, p);
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
	} else if ((l->t == MAT && r->t == VEC) ||
		   (l->t == VEC && r->t == MAT)) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else if (t->t >= B_EQ && t->t <= B_NEQ &&
		   ((l->t == VEC && r->t == STR) ||
		    (l->t == STR && r->t == VEC))) {
	    ret = series_string_calc(l, r, t->t, p);
	} else if ((t->t == B_AND || t->t == B_OR || t->t == B_SUB) &&
		   ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_op(l, r, t->t, p);
	} else if (t->t == B_POW && ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_op(l, r, t->t, p);
	} else if ((t->t == B_EQ || t->t == B_NEQ) && 
		   ((l->t == NUM && r->t == STR) ||
		    (l->t == STR && r->t == NUM))) {
	    ret = num_string_comp(l, r, t->t, p);
	} else if (bool_comp(t->t)) {
	    if (ok_list_node(l) && (r->t == NUM || r->t == VEC)) {
		ret = list_bool_comp(l, r, t->t, 0, p);
	    } else if (ok_list_node(r) && (l->t == NUM || l->t == VEC)) {
		ret = list_bool_comp(r, l, t->t, 1, p);
	    } else if (ok_list_node(l) && ok_list_node(r)) {
		ret = list_list_comp(r, l, t->t, p);
	    } else {
		p->err = E_TYPES;
	    }
	} else if ((t->t == B_ADD || t->t == B_SUB) && 
		   l->t == VEC && ok_list_node(r)) {
	    ret = series_list_calc(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case B_TRMUL:
	/* matrix on left, otherwise be flexible */
	if (l->t == MAT && r->t == MAT) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if (l->t == MAT && r->t == VEC) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == NUM) ||
		   (l->t == NUM && r->t == MAT)) {
	    ret = matrix_scalar_calc(l, r, t->t, p);
	} else if (l->t == MAT && r->t == EMPTY) {
	    ret = matrix_transpose_node(l, p);
	} else if (l->t == NUM && r->t == EMPTY) {
	    ret = scalar_transpose_node(l, p);
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
	/* matrix-matrix or matrix-scalar binary operators */
	if ((l->t == MAT && r->t == MAT) ||
	    (l->t == MAT && r->t == NUM) ||
	    (l->t == NUM && r->t == MAT)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == VEC) ||
		   (l->t == VEC && r->t == MAT)) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == MAT)? 2 : 1,
			    MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case B_HCAT:
    case B_VCAT:
    case F_QFORM:
    case F_HDPROD:	
    case F_CMULT:
    case F_CDIV:
    case F_MRSEL:
    case F_MCSEL:
    case F_DSUM:
    case B_LDIV:	
	/* matrix-only binary operators (but promote scalars) */
	if ((l->t == MAT || l->t == NUM) && 
	    (r->t == MAT || r->t == NUM)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if (t->t == B_HCAT && l->t == STR && r->t == STR) {
	    /* exception: string concatenation */
	    ret = two_string_func(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == BUNDLE) ||
		   (l->t == BUNDLE && r->t == MAT)) {
	    ret = matrix_bundle_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == MAT)? 2 : 1,
			    MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case B_KRON:
	/* Kronecker product ("**"): insist on matrices only */
	if (l->t == MAT && r->t == MAT) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == MAT)? 2 : 1,
			    MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case B_ELLIP:
	/* list-making ellipsis */
	if (useries_node(l) && useries_node(r)) {
	    ret = ellipsis_list_node(l, r, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;
    case B_JOIN:
	/* list join with separator */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = list_join_node(l, r, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;	
    case F_MSORTBY:
	/* matrix on left, scalar on right */
	if (l->t == MAT && scalar_node(r)) {
	    ret = matrix_scalar_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;
    case F_LLAG:
	if (scalar_node(l) && ok_list_node(m)) {
	    ret = list_make_lags(l, m, r, t->t, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;
    case U_NEG: 
    case U_POS:
    case U_NOT:
    case F_ABS:
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
    case F_INVMILLS:
    case F_EASTER:
	/* functions taking one argument, any type */
	if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (l->t == VEC) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == MAT) {
	    ret = apply_matrix_func(l, t->t, p);
	} else if (ok_list_node(l) && t->t == F_LOG) {
	    ret = apply_list_func(l, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_DUMIFY:
	/* series argument wanted */ 
	if (ok_list_node(l)) {
	    ret = dummify_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_MISSZERO:
    case F_ZEROMISS:
	/* one series or scalar argument needed */
	if (l->t == VEC || l->t == MAT) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	}
	break;
    case F_MISSING:	
    case F_DATAOK:
	/* series, scalar or list argument needed */
	if (l->t == MAT) {
	    if (t->t == F_DATAOK) {
		ret = matrix_to_matrix_func(l, NULL, t->t, p);
	    } else {
		node_type_error(t->t, 0, VEC, l, p);
	    }
	} else if (l->t == VEC) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (ok_list_node(l)) {
	    ret = list_ok_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
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
    case LAG:
	if (p->targ == LIST) {
	    ret = get_lag_list(l, r, p);
	    break;
	}
	/* note: otherwise fall through */
    case F_LJUNGBOX:
    case F_POLYFIT:
	/* series on left, scalar on right */
	if (l->t != VEC) {
	    node_type_error(t->t, 1, VEC, l, p);
	} else if (!scalar_node(r)) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else if (t->t == LAG) {
	    ret = series_lag(l, r, p); 
	} else if (t->t == F_LJUNGBOX) {
	    ret = series_ljung_box(l, r, p); 
	} else if (t->t == F_POLYFIT) {
	    ret = series_polyfit(l, r, p);
	} 
	break;
    case OBS:
    case DOBS:
	if (l->t != VEC) {
	    node_type_error(t->t, 1, VEC, l, p);
	} else if (!scalar_node(r) && r->t != STR) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    ret = series_obs(l, r, p); 
	}
	break;
    case F_MOVAVG:
	/* series on left, plus one or two scalars */
	if (l->t != VEC) {
	    node_type_error(t->t, 1, VEC, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(t->t, 2, NUM, m, p);
	} else if (!empty_or_num(r)) {
	    node_type_error(t->t, 3, NUM, r, p);
	} else {
	    ret = series_movavg(l, m, r, p);
	}
	break;
    case MSL:
	/* user matrix plus subspec */
	ret = submatrix_node(l, r, p);
	break;
    case MSL2:
	/* unevaluated matrix subspec */
	ret = mspec_node(l, r, p);
	break;
    case SUBSL:
    case B_RANGE:
	/* matrix sub-slice, x:y, or lag range, 'p to q' */
	ret = process_subslice(l, r, p);
	break;
    case LISTELEM:
	/* list plus scalar */
	if (r->t != NUM) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    ret = get_named_list_element(l, r, p);
	}
	break;
    case BOBJ:
    case BMEMB:
    case F_INBUNDLE:
	/* name of bundle plus string */
	if (l->t == BUNDLE && r->t == STR) {
	    if (t->t == BOBJ || t->t == BMEMB) {
		ret = get_named_bundle_value(l, r, p);
	    } else {
		ret = test_bundle_key(l, r, p);
	    }
	} else if (l->t == BUNDLE) {
	    node_type_error(t->t, 1, STR, r, p);
	} else {
	    node_type_error(t->t, 0, BUNDLE, l, p);
	}	    
	break;
    case F_TYPESTR:
	/* numerical type code to string */
	if (l->t == NUM) {
	    ret = type_string_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_LDIFF:
    case F_SDIFF:
    case F_ODEV:	
	if (l->t == VEC || (t->t != F_ODEV && l->t == MAT)) {
	    ret = series_series_func(l, r, t->t, p);
	} else if (ok_list_node(l)) {
	    ret = apply_list_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
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
    case F_PSD:
    case F_DESEAS:
	/* series argument needed */
	if (l->t == VEC || l->t == MAT) {
	    ret = series_series_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	} 
	break;
    case F_FREQ:
	/* series -> matrix */
	if (l->t == VEC || l->t == MAT) {
	    ret = series_matrix_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	} 
	break;	
    case F_PSHRINK:
	if (l->t == VEC) {
	    ret = do_panel_shrink(l, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	}
	break;
    case F_CUM:
    case F_DIFF:
    case F_RESAMPLE:
    case F_RANKING:
	/* series or matrix argument */
	if (l->t == VEC) {
	    ret = series_series_func(l, r, t->t, p);
	} else if (l->t == MAT) {
	    ret = matrix_to_matrix_func(l, r, t->t, p);
	} else if (t->t == F_DIFF && ok_list_node(l)) {
	    ret = apply_list_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	}
	break;
    case F_SORT:
    case F_DSORT:
    case F_VALUES:
    case F_UNIQ:
    case F_PERGM:
    case F_IRR:
	/* series or vector argument needed */
	if (l->t == VEC || l->t == MAT) {
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
	    node_type_error(t->t, 0, VEC, l, p);
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
    case F_MIN:
    case F_MAX:
    case F_MEDIAN:
    case F_GINI:
    case F_NOBS:
    case F_T1:
    case F_T2:
	/* functions taking series arg, returning scalar */
	if (l->t == VEC || l->t == MAT) {
	    ret = series_scalar_func(l, t->t, p);
	} else if ((t->t == F_MEAN || t->t == F_SD || 
		    t->t == F_VCE || t->t == F_MIN ||
		    t->t == F_MAX || t->t == F_SUM) 
		   && ok_list_node(l)) {
	    /* list -> series also acceptable for these cases */
	    ret = list_to_series_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, VEC, l, p);
	} 
	break;	
    case F_LRVAR:
    case F_NPV:
    case F_ISCONST:
	/* takes series and scalar arg, returns scalar */
	if (l->t == VEC || l->t == MAT) {
	    if (t->t == F_ISCONST) {
		ret = isconst_node(l, r, p);
	    } else if (scalar_node(r)) {
		if (t->t == F_QUANTILE && l->t == MAT) {
		    ret = matrix_quantiles_node(l, r, p);
		} else {
		    ret = series_scalar_scalar_func(l, r, t->t, p);
		}
	    } else {
		node_type_error(t->t, 2, NUM, r, p);
	    } 
	} else {
	    node_type_error(t->t, 1, VEC, l, p);
	}
	break;
    case F_QUANTILE:
	if (l->t == VEC) {
	    if (scalar_node(r)) {
		ret = series_scalar_scalar_func(l, r, t->t, p);
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
	    node_type_error(t->t, 1, (r->t == MAT)? MAT : VEC,
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
    case F_COR:
    case F_COV:
    case F_FCSTATS:
	/* functions taking two series/vectors as args */
	if (l->t == VEC && r->t == VEC) {
	    ret = series_2_func(l, r, t->t, p);
	} else if (l->t == MAT && r->t == MAT) {
	    ret = series_2_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == VEC)? 2 : 1,
			    VEC, (l->t == VEC)? r : l, p);
	} 
	break;
    case F_MXTAB:
	/* functions taking two series or matrices as args and returning 
	   a matrix */
	if ((l->t == VEC && r->t == VEC) || (l->t == MAT && r->t == MAT)) {
	    ret = mxtab_func(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == VEC)? 2 : 1,
			    VEC, (l->t == VEC)? r : l, p);
	} 
	break;
    case F_SORTBY:
	/* takes two series as args, returns series */
	if (l->t == VEC && r->t == VEC) {
	    ret = series_sort_by(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == VEC)? 2 : 1,
			    VEC, (l->t == VEC)? r : l, p);
	} 
	break;	
    case F_IMAT:
    case F_ZEROS:
    case F_ONES:
    case F_MUNIF:
    case F_MNORM:
	/* matrix-creation functions */
	if (scalar_node(l) && (r == NULL || scalar_node(r))) {
	    ret = matrix_fill_func(l, r, t->t, p);
	} else if (!scalar_node(l)) {
	    node_type_error(t->t, 1, NUM, l, p);
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
    case F_SDC:
    case F_MCOV:
    case F_MCORR:
    case F_CDEMEAN:
    case F_CHOL:
    case F_PSDROOT:
    case F_INV:
    case F_INVPD:
    case F_GINV:
    case F_DIAG:
    case F_TRANSP:
    case F_MREVERSE:
    case F_VEC:
    case F_VECH:
    case F_UNVECH:
    case F_UPPER:
    case F_LOWER:
    case F_NULLSPC:
    case F_MEXP:
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
	/* matrix -> matrix functions */
	if (l->t == MAT || l->t == NUM) {
	    ret = matrix_to_matrix_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, 1, MAT, l, p);
	}
	break;
    case F_ROWS:
    case F_COLS:
    case F_DET:
    case F_LDET:
    case F_TRACE:
    case F_NORM1:
    case F_INFNORM:
    case F_RCOND:
    case F_RANK:
	/* matrix -> scalar functions */
	if (l->t == MAT || l->t == NUM) {
	    ret = matrix_to_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, MAT, l, p);
	}
	break;
    case F_MREAD:
    case F_BREAD:
	if (l->t != STR) {
	    node_type_error(t->t, 1, STR, l, p);	
	} else if (!empty_or_num(r)) {
	    node_type_error(t->t, 2, NUM, r, p);
	} else {
	    ret = read_object_func(l, r, t->t, p);
	}
	break;
    case F_QR:
    case F_EIGSYM:
    case F_EIGGEN:
	/* matrix -> matrix functions, with indirect return */
	if (l->t != MAT) {
	    node_type_error(t->t, 1, MAT, l, p);
	} else if (r->t != U_ADDR && r->t != EMPTY) {
	    node_type_error(t->t, 2, U_ADDR, r, p);
	} else {
	    ret = matrix_to_matrix2_func(l, r, t->t, p);
	}
	break;
    case F_FDJAC:
    case F_MWRITE:
	/* matrix, with string as second arg */
	if (t->t == F_FDJAC && l->t == MAT && r->t == STR) {
	    ret = numeric_jacobian(l, r, p);
	} else if (t->t == F_MWRITE && l->t == MAT && 
		   m->t == STR && empty_or_num(r)) {
	    ret = matrix_text_write(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;
    case F_BWRITE:
	/* bundle, with string as second arg */
	if (l->t == BUNDLE && m->t == STR && empty_or_num(r)) {
	    ret = bundle_text_write(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_BFGSMAX:
	/* matrix, plus one or two string args */
	if (l->t == MAT && m->t == STR) {
	    ret = BFGS_maximize(l, m, r, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;
    case F_SIMANN:
	/* matrix, plus string and int args */
	if (l->t == MAT && m->t == STR) {
	    ret = simann_node(l, m, r, p);
	} else {
	    p->err = E_TYPES;
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
    case F_COLNAMES:
    case F_ROWNAMES:
	/* matrix, list or string as second arg */
	if (l->t == MAT && (ok_list_node(r) || r->t == STR)) {
	    ret = matrix_add_names(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;
    case F_COLNAME:
	/* matrix, scalar as second arg */
	if (l->t == MAT && r->t == NUM) {
	    ret = matrix_get_colname(l, r, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;
    case F_XMIN:
    case F_XMAX:
    case F_RANDINT:
	/* two scalars */
	if (l->t == NUM && r->t == NUM) {
	    ret = two_scalars_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;	
    case F_MSHAPE:
    case F_SVD:
    case F_TRIMR:
    case F_TOEPSOLV:
    case F_CORRGM:
    case F_SEQ:
    case F_REPLACE:
    case F_STRNCMP:
    case F_WEEKDAY:
    case F_MONTHLEN:
    case F_EPOCHDAY:
    case F_KDENSITY:
    case F_SETNOTE:
    case F_BWFILT:
    case F_CHOWLIN:
    case F_VARSIMUL:
    case F_IRF:
    case F_STRSUB:
    case F_REGSUB:
    case F_MLAG:
    case F_EIGSOLVE:
    case F_NADARWAT:
    case F_PRINCOMP:
    case F_HALTON:
    case F_AGGRBY:
    case F_IWISHART:
    case F_SUBSTR:
	/* built-in functions taking three args */
	if (t->t == F_REPLACE) {
	    ret = replace_value(l, m, r, p);
	} else if (t->t == F_STRSUB || t->t == F_REGSUB) {
	    ret = string_replace(l, m, r, t->t, p);
	} else if (t->t == F_EPOCHDAY) {
	    ret = eval_epochday(l, m, r, p);
	} else {
	    ret = eval_3args_func(l, m, r, t->t, p);
	}
	break;
    case F_PRINTF:
	if (l->t == STR && empty_or_string(r)) {
	    ret = eval_print_scan(NULL, l, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, NULL, p);
	}
	break;	
    case F_SPRINTF:
	if (l->t == STR && m->t == STR && empty_or_string(r)) {
	    ret = eval_print_scan(l, m, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, NULL, p);
	}
	break;		
    case F_SSCANF:
	if (l->t == STR && m->t == STR && r->t == STR) {
	    ret = eval_print_scan(l, m, r, t->t, p);
	} else {
	    node_type_error(t->t, 0, STR, NULL, p);
	}
	break;
    case F_BESSEL:
	/* functions taking one char, one scalar/series and one 
	   matrix/series/scalar as args */
	if (l->t != STR) {
	    node_type_error(t->t, 1, STR, l, p);
	} else if (!scalar_node(m)) {
	    node_type_error(t->t, 2, NUM, m, p);
	} else if (r->t != NUM && r->t != VEC && r->t != MAT) {
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
    case F_KFILTER:
    case F_KSMOOTH:
    case F_KSIMUL:
    case F_NRMAX:
    case F_LOESS:
    case F_GHK:
    case F_QUADTAB:
	/* built-in functions taking more than three args */
	ret = eval_nargs_func(t, p);
	break;
    case F_ISOCONV:
	ret = isoconv_node(t, p);
	break;
    case USTR:
	/* string variable */
	ret = string_var_node(t, p);
	break;
    case OVAR:
    case MVAR:
    case DMSL:
	/* variable "under" user-defined object, or last object */
	ret = object_var_node(t, p);
	break;
    case DMSTR:
	ret = dollar_str_node(t, NULL, p);
	break;
    case DVAR: 
	/* dataset "dollar" variable */
	ret = dollar_var_node(t, p);
	break;
    case MDEF:
	/* matrix definition */
	ret = matrix_def_node(t, p);
	break;
    case F_OBSNUM:
    case F_ISSERIES:
    case F_ISLIST:
    case F_ISSTRING:
    case F_ISNULL:
    case F_STRLEN:
    case F_REMOVE:
    case F_TYPEOF:
	if (l->t == STR) {
	    ret = object_status(l, t->t, p);
	} else {
	    node_type_error(t->t, 1, STR, l, p);
	}
	break;
    case F_LISTLEN:
	if (ok_list_node(l) || l->t == STR) {
	    ret = list_length_node(l, p);
	} else {
	    node_type_error(t->t, 1, LIST, l, p);
	}
	break;
    case F_INLIST:
	if (ok_list_node(l)) {
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
	if (t->v.b1.b->t == FARGS) {
	    if (t->t == F_URCPVAL) {
		ret = eval_urcpval(t, p);
	    } else {
		ret = eval_pdist(t, p);
	    }
	} else {
	    node_type_error(t->t, 0, FARGS, t->v.b1.b, p);
	}
	break;	
    case CON:
	/* built-in constant */
	ret = retrieve_const(t, p);
	break;
    case EROOT:
	ret = eval(t->v.b1.b, p);
	break;
    case UFUN:
	ret = eval_ufunc(t, p);
	break;
#ifdef USE_RLIB
    case RFUN:
	ret = eval_Rfunc(t, p);
	break;
#endif
    case QUERY:
	ret = eval_query(t, p);
	break;
    case B_LCAT:
	/* list concatenation */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = eval_lcat(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case F_XPX:
	if (ok_list_node(l)) {
	    ret = apply_list_func(l, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;	
    case F_WMEAN:
    case F_WVAR:
    case F_WSD:
	/* two lists -> series */
	if (ok_list_node(l) && ok_list_node(r)) {
	    ret = list_list_series_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;	
    case F_LINCOMB:
	/* list + matrix -> series */
	if (ok_list_node(l) && r->t == MAT) {
	    ret = list_matrix_series_func(l, r, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;	
    case F_ARGNAME:
    case F_BACKTICK:
    case F_STRSTRIP:
    case F_FIXNAME:
	if (l->t == STR) {
	    ret = single_string_func(l, t->t, p);
	} else if (t->t == F_ARGNAME && (uscalar_node(l) || useries_node(l))) {
	    ret = argname_from_uvar(l, p);
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
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
    case F_OBSLABEL:
	if (l->t == NUM || l->t == MAT) {
	    ret = int_to_string_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_VARNAME:
	if (l->t == NUM || l->t == MAT) {
	    ret = int_to_string_func(l, t->t, p);
	} else if (l->t == LIST) {
	    ret = list_to_string_func(l, t->t, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_VARNUM:
    case F_TOLOWER:
    case F_TOUPPER:
	if (l->t == STR) {
	    if (t->t == F_TOLOWER || t->t == F_TOUPPER) {
		ret = one_string_func(l, t->t, p);
	    } else {
		ret = varnum_node(l, p);
	    }
	} else {
	    node_type_error(t->t, 0, STR, l, p);
	}
	break;
    case F_STRSTR:
	if (l->t == STR && r->t == STR) {
	    ret = two_string_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    STR, (l->t == STR)? r : l, p);
	}
	break;
    case F_STRSPLIT:
	if (l->t == STR && r->t == NUM) {
	    ret = gretl_strsplit(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    (l->t == STR)? NUM : STR,
			    (l->t == STR)? r : l, p);
	}
	break;
    case F_GETLINE:
	if (l->t == STR && (r->t == STR || null_or_empty(r))) {
	    ret = getline_node(l, r, p);
	} else {
	    node_type_error(t->t, (l->t == STR)? 2 : 1,
			    STR, (l->t == STR)? r : l, p);
	}
	break;
    case F_ERRMSG:
	if (l->t == NUM || l->t == EMPTY) {
	    ret = errmsg_node(l, p);
	} else {
	    node_type_error(t->t, 0, NUM, l, p);
	}
	break;
    case F_ISODATE:
	ret = isodate_node(l, r, p);
	break;
    case F_ATOF:
	if (l->t == STR) {
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
	} else if (!empty_or_num(r)) {
	    node_type_error(t->t, 3, NUM, r, p);
	} else {
	    ret = mpi_transfer_node(l, m, r, t->t, p);
	}
	break;
    case F_GENSERIES:
	ret = gen_series_node(l, r, p);
	break;
    default: 
	printf("eval: weird node %s (t->t = %d)\n", getsymb(t->t, NULL),
	       t->t);
	p->err = E_PARSE;
	break;
    }

 bailout:

#if EDEBUG
    fprintf(stderr, "eval (t->t = %03d, %s): returning NODE at %p\n", 
	    t->t, getsymb(t->t, NULL), (void *) ret);
    if (t->t == VEC) 
	fprintf(stderr, " (VEC node, xvec at %p, vnum = %d)\n", 
		(void *) t->v.xvec, t->vnum);
#endif

    return ret;
}

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

/* Look ahead to the next non-space character in the
   parser stream and return it; if @skip then start at
   offset 1 beyond the current p->point.
*/

int parser_next_nonspace_char (parser *p, int skip)
{
    int i, offset = skip ? 1 : 0;

    if (p->point == '\0') {
	return 0;
    }

    for (i=offset; p->point[i] != '\0'; i++) {
	if (!isspace(p->point[i])) {
	    return p->point[i];
	}
    }

    return 0;
}

/* gets the next character in the input stream */

int parser_next_char (parser *p)
{
    return p->point[0];
}

/* for error reporting: print the input up to the current
   parse point */

void parser_print_input (parser *p)
{
    int pos = p->point - p->input;
    char *s;

    s = gretl_strndup(p->input, pos);
    if (s != NULL) {
	pprintf(p->prn, "> %s\n", s);
	free(s);
    }
}

/* "pretty print" syntatic nodes and symbols */

static void printsymb (int symb, const parser *p)
{
    pputs(p->prn, getsymb(symb, NULL));
}

static void printnode (NODE *t, parser *p)
{  
    if (t == NULL) {
	pputs(p->prn, "NULL"); 
    } else if (useries_node(t)) {
	pprintf(p->prn, "%s", p->dset->varname[t->vnum]);
    } else if (uscalar_node(t)) {
	pprintf(p->prn, "%s", t->vname);
    } else if (t->t == NUM) {
	if (na(t->v.xval)) {
	    pputs(p->prn, "NA");
	} else {
	    pprintf(p->prn, "%.8g", t->v.xval);
	}
    } else if (t->t == VEC) {
	const double *x = t->v.xvec;
	int i, j = 1;

	if (p->lh.v > 0 && p->lh.v < p->dset->v) {
	    pprintf(p->prn, "%s\n", p->dset->varname[p->lh.v]);
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
	printnode(t->v.b2.l, p);
	printsymb(t->t, p);
	printnode(t->v.b2.r, p);
	pputc(p->prn, ')');
    } else if (t->t == MSL || t->t == DMSL) {
	printnode(t->v.b2.l, p);
	pputc(p->prn, '[');
	printnode(t->v.b2.r, p);
	pputc(p->prn, ']');
    } else if (t->t == MSL2) {
	pputs(p->prn, "MSL2");
    } else if (t->t == SUBSL) {
	pputs(p->prn, "SUBSL");
    } else if (t->t == OVAR) {
	printnode(t->v.b2.l, p);
	pputc(p->prn, '.');
	printnode(t->v.b2.r, p);
    } else if (func1_symb(t->t)) {
	printsymb(t->t, p);
	pputc(p->prn, '(');
	printnode(t->v.b1.b, p);
	pputc(p->prn, ')');
    } else if (unary_op(t->t)) {
	printsymb(t->t, p);
	printnode(t->v.b1.b, p);
    } else if (t->t == EROOT) {
	printnode(t->v.b1.b, p);
    } else if (func2_symb(t->t)) {
	printsymb(t->t, p);
	pputc(p->prn, '(');
	printnode(t->v.b2.l, p);
	if (t->v.b2.r->t != EMPTY) {
	    pputc(p->prn, ',');
	}
	printnode(t->v.b2.r, p);
	pputc(p->prn, ')');
    } else if (t->t == STR || t->t == USTR) {
	pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == MDEF) {
	pprintf(p->prn, "{ MDEF }");
    } else if (t->t == DMSTR || t->t == UFUN) {
	printnode(t->v.b2.l, p);
	pputc(p->prn, '(');
	printnode(t->v.b2.r, p);
	pputc(p->prn, ')');
    } else if (t->t == LISTVAR) {
	pprintf(p->prn, "%s.%s", t->v.b2.l->v.str, t->v.b2.r->v.str);
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
   will we accept, when generating a matrix? */

#define ok_matrix_op(o) (o == B_ASN || o == B_ADD || \
			 o == B_SUB || o == B_MUL || \
			 o == B_DIV || o == INC || \
			 o == DEC || o == B_HCAT || \
                         o == B_VCAT || o == B_DOTASN)
#define ok_list_op(o) (o == B_ASN || o == B_ADD || o == B_SUB)
#define ok_string_op(o) (o == B_ASN || o == B_ADD || \
			 o == B_HCAT || o == INC) 

struct mod_assign {
    int c;
    int op;
};

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

/* read operator from "genr" formula: this is either
   simple assignment or something like '+=' */

static int get_op (char *s)
{
    if (s[0] == '=') {
	s[1] = '\0';
	return B_ASN;
    }

    if (!strcmp(s, "++")) {
	return INC;
    }

    if (!strcmp(s, "--")) {
	return DEC;
    }

    if (s[1] == '=') {
	int i;

	for (i=0; m_assign[i].c; i++) {
	    if (s[0] == m_assign[i].c) {
		return m_assign[i].op;
	    }
	}
    }

    return 0;
}

/* extract a substring [...] from the left-hand side
   of a genr expression */

static void get_lhs_substr (char *str, parser *p)
{
    char *q, *s = strchr(str, '[');

#if EDEBUG
    fprintf(stderr, "get_lhs_substr: str = '%s'\n", str);
#endif

    q = gretl_strdup(s + 1);
    if (q == NULL) {
	p->err = E_ALLOC;
    } else {
	int n = strlen(q);

	if (q[n-1] != ']') {
	    /* error message */
	    p->err = E_PARSE;
	} else {
	    q[n-1] = '\0';
	}
	p->lh.substr = q;
    }

    *s = '\0';
}

static void get_lhs_dotstr (char *str, parser *p)
{
    char *q, *s = strchr(str, '.');

    q = gretl_strdup(s + 1);
    if (q == NULL) {
	p->err = E_ALLOC;
    } else {
	p->lh.substr = q;
    }

#if EDEBUG
    fprintf(stderr, "get_lhs_dotstr: '%s'\n", p->lh.substr);
#endif

    *s = '\0';
}

static void nullify_aux_return (parser *p)
{
    int i;

    if (p->aux != NULL) {
	for (i=0; i<p->n_aux; i++) {
	    if (p->aux[i] == p->ret) {
		p->aux[i] = NULL;
		return;
	    }
	}
    }
}

/* Given a string [...], parse and evaluate it as a
   sub-matrix specification.  This is for the case where
   assignment is to a submatrix, as in m[spec] = foo.
*/

static void get_lh_mspec (parser *p)
{
    parser *subp = p->subp;

#if LHDEBUG
    fprintf(stderr, "get_lh_mspec: %s\n", (p->flags & P_COMPILE)?
	    "compiling" : "running");
#endif

    if (subp != NULL) {
	/* executing a previously compiled parser */
	parser_free_aux_nodes(subp);
	parser_aux_init(subp);
	parser_reinit(subp, p->dset, p->prn);
    } else {
	/* starting from scratch */
	int subflags = P_SLICE;
	char *s;

	p->subp = subp = malloc(sizeof *p->subp);
	if (subp == NULL) {
	    p->err = E_ALLOC;
	    return;
	}
	s = malloc(strlen(p->lh.substr) + 3);
	if (s == NULL) {
	    p->err = E_ALLOC;
	    return;
	}
	sprintf(s, "[%s]", p->lh.substr);
	if (p->flags & P_COMPILE) {
	    subflags |= P_COMPILE;
	}
	parser_init(subp, s, p->dset, p->prn, subflags);
	parser_aux_init(subp);
	subp->targ = MSPEC;
# if LHDEBUG
	fprintf(stderr, "subp->input='%s'\n", subp->input);
# endif
	subp->tree = msl_node_direct(subp);
	free(s);
    }

    p->err = subp->err;

    if (subp->tree != NULL && !(p->flags & P_COMPILE)) {
	/* evaluate subp to get a matrix subspec */
	parser_aux_init(subp);
	subp->ret = eval(subp->tree, subp);

	if (subp->err) {
	    printf("Error in subp eval = %d\n", subp->err);
	    p->err = subp->err;
	} else {
	    /* free previous result, if any, and appropriate the
	       mspec result from current evaluation */
	    if (p->lh.mspec != NULL) {
		free_mspec(p->lh.mspec);
	    }
	    p->lh.mspec = subp->ret->v.mspec;
#if LHDEBUG > 1
	    print_mspec(p->lh.mspec);
#endif
	    nullify_aux_return(subp); /* don't double-free */
	    free(subp->ret);
	    subp->ret = NULL;
	}
    } 
}

/* Given a string [...], parse and evaluate it as a series observation
   index.  This is for the case where assignment is to a specific
   observation, as in y[obs] = foo.
*/

static void get_lh_obsnum (parser *p)
{
    int done = 0;

    if (p->lh.substr[0] != '"') {
	int err = 0;

	p->lh.obs = generate_scalar(p->lh.substr, p->dset, &err);
	if (!err) {
	    p->lh.obs -= 1; /* convert to 0-based for internal use */
	    done = 1;
	}
    }

    if (!done) {
	p->lh.obs = get_t_from_obs_string(p->lh.substr, p->dset);
    }
    
    if (p->lh.obs < 0 || p->lh.obs >= p->dset->n) {
	gretl_errmsg_sprintf("'[%s]': bad observation specifier", p->lh.substr);
	p->err = E_PARSE;
    } else {
	gretl_error_clear();
	p->targ = NUM;
    }
}

/* check validity of "[...]" on the LHS, and evaluate
   the expression if needed */

static void process_lhs_substr (const char *lname, char c, parser *p)
{
#if LHDEBUG || EDEBUG
    fprintf(stderr, "process_lhs_substr: p->lh.t = %d, substr = '%s'\n", 
	    p->lh.t, p->lh.substr);
#endif

    /* FIXME support stringvar[] ? */

    if (p->lh.t == VEC) {
	get_lh_obsnum(p);
    } else if (p->lh.t == MAT) {
	get_lh_mspec(p);
    } else if (p->lh.t == BUNDLE) {
	/* substr should be key string; handled later */
	if (c != 0) {
	    p->targ = (c == '.')? BMEMB : BOBJ;
	}
    } else if (p->lh.t == UNK) {
	if (lname != NULL) {
	    undefined_symbol_error(lname, p);
	}
	p->err = E_UNKVAR;
    } else {
	if (c != 0) {
	    gretl_errmsg_sprintf(_("The symbol '%c' is not valid in this context\n"), c);
	}
	p->err = E_DATA;
    }	
}

#if EDEBUG
static void parser_print_result (parser *p, PRN *prn)
{
    if (p->targ == VEC) {
	int list[2] = { 1, p->lh.v };

	printdata(list, NULL, p->dset, OPT_NONE, prn);
    } else if (p->targ == NUM) {
	printdata(NULL, p->lh.name, p->dset, OPT_NONE, prn);
    } else if (p->targ == MAT) {
	gretl_matrix_print_to_prn(p->lh.m1, p->lh.name, prn);
    }
}
#endif

/* implement the declaration of new variables */

static void do_decl (parser *p)
{
    char **S = NULL;
    int i, v, n;

    n = check_declarations(&S, p);

    if (n == 0) {
	return;
    }

    for (i=0; i<n && !p->err; i++) {
	if (S[i] != NULL) {
	    if (p->targ == VEC) {
		p->err = dataset_add_NA_series(p->dset);
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
		} else {
		    p->err = E_DATA;
		}
		if (!p->err) {
		    p->err = create_user_var(S[i], type);
		}
	    }
	}
    }

    strings_array_free(S, n);
}

/* create a dummy node to facilitate (a) printing an
   existing variable, or (b) incrementing or decrementing
   that variable
*/

static NODE *lhs_copy_node (parser *p)
{
    NODE *n = new_node(p->targ);

    if (n != NULL) {
	if (p->targ == NUM) {
	    n->v.xval = gretl_scalar_get_value(p->lh.name, NULL);
	} else if (p->targ == VEC) {
	    n->v.xvec = p->dset->Z[p->lh.v];
	} else if (p->targ == STR) {
	    n->v.str = gretl_strdup(get_string_by_name(p->lh.name));
	} else {
	    n->v.m = p->lh.m0;
	}
    }

    return n;
}

/* The expression supplied for evaluation does not contain an '=': can
   we parse it as a declaration of a new variable, or as an implicit
   request to print the value of an existing variable?
*/

static void parser_try_print (parser *p)
{
    if (p->lh.v == 0 && p->lh.m0 == NULL && !lhlist(p) && 
	!lhstr(p) && !lhscalar(p) && !lhbundle(p)) {
	/* varname on left is not the name of a current variable */
	p->err = E_EQN;
    } else if (p->lh.substr != NULL) {
	/* could perhaps be construed as a valid print request? */
	p->err = E_EQN;
    } else if (p->targ != p->lh.t) {
	/* attempt to re-declare a variable with a different type */
	p->err = E_TYPES;
    } else {
	/* e.g. "series x", for an existing series named 'x' */
	p->flags |= (P_PRINT | P_DISCARD);
#if EDEBUG
	fprintf(stderr, "parser_try_print: set print/discard flags\n");
#endif
    }
}

static int extract_LHS_string (const char *s, char *lhs, parser *p)
{
    int n, b = 0;

    *lhs = '\0';

    if (p->targ != UNK && strchr(s, '=') == NULL) {
	/* variable declaration(s) ? */
	p->flags |= P_DECL;
	p->lh.substr = gretl_strdup(s);
	return 0;
    }

    n = strcspn(s, "+-*/%^~|([= ");

    if (n > 0) {
	if (*(s+n) == '[') {
	    const char *q = s + n;

	    while (*q) {
		if (*q == '[') {
		    b++;
		} else if (*q == ']') {
		    b--;
		}
		n++;
		if (b == 0) {
		    break;
		}
		q++;
	    }
	    if (b != 0) {
		pprintf(p->prn, "> %s\n", s);
		pprintf(p->prn, _("Unmatched '%c'\n"), '[');
	    }
	}
    }

    if (n > 0 && n < GENSTRLEN && b == 0) {
	strncat(lhs, s, n);
    }

    return (*lhs == '\0')? E_PARSE : 0;
}

/* in the case of a "private" genr we allow ourselves some
   more latitude in variable names, so as not to collide
   with userspace names: specifically, we can use '$'
*/

static int check_private_varname (const char *s)
{
    const char *ok = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"0123456789_$";
    int n = 0, err = 0;

    if (isalpha(*s) || *s == '$') {
	n = strspn(s, ok);
    }

    if (n != strlen(s)) {
	err = E_PARSE;
    }

    return err;
}

static void maybe_do_type_errmsg (const char *name, int t)
{
    if (name != NULL && *name != '\0') {
	const char *tstr = NULL;

	if (t == NUM) {
	    tstr = "scalar";
	} else if (t == VEC || t == UVEC) {
	    tstr = "series";
	} else if (t == MAT) {
	    tstr = "matrix";
	} else if (t == STR) {
	    tstr = "string";
	} else if (t == BUNDLE) {
	    tstr = "bundle";
	} else if (t == LIST) {
	    tstr = "list";
	}

	if (tstr != NULL) {
	    gretl_errmsg_sprintf(_("The variable %s is of type %s"), 
				 name, tstr);
	}
    }
}

static int overwrite_type_check (parser *p)
{
    int err = 0;

    if (p->targ == NUM && p->lh.t == VEC && p->lh.obs >= 0) {
	; /* OK */
    } else if ((p->targ == BOBJ || p->targ == BMEMB) && p->lh.t == BUNDLE) {
	; /* OK */
    } else if (p->targ != p->lh.t) {
	/* don't overwrite one type with another */
	maybe_do_type_errmsg(p->lh.name, p->lh.t);
	err = E_TYPES;
    }

    return err;
}

static int overwrite_const_check (const char *s, parser *p)
{
    if (object_is_const(s)) {
	p->err = overwrite_err(s);
	return p->err;
    } else {
	return 0;
    }
}

static void maybe_set_matrix_target (parser *p)
{
    int n = strlen(p->rhs);

    if (n > 1 && p->rhs[n-1] == '}') {
	p->targ = MAT;
    }
}

/* process the left-hand side of a genr formula */

static void pre_process (parser *p, int flags)
{
    const char *s = p->input;
    char test[GENSTRLEN];
    char opstr[3] = {0};
    char subchar = 0;
    int v, newvar = 1;

    while (isspace(*s)) s++;

    /* skip leading command word, if any */
    if (!strncmp(s, "genr ", 5)) {
	s += 5;
    } else if (!strncmp(s, "eval ", 5)) {
	p->flags |= P_DISCARD;
	s += 5;
    } else if (!strncmp(s, "print ", 6)) {
	p->flags |= P_PRINT;
	s += 6;
    }

    while (isspace(*s)) s++;

    /* do we have a type specification? */
    if (flags & P_SCALAR) {
	p->targ = NUM;
    } else if (flags & P_SERIES) {
	p->targ = VEC;
    } else if (flags & P_MATRIX) {
	p->targ = MAT;
    } else if (flags & P_STRING) {
	p->targ = STR;
    } else if (flags & P_LIST) {
	p->targ = LIST;
    } else if (!strncmp(s, "scalar ", 7)) {
	p->targ = NUM;
	s += 7;
    } else if (!strncmp(s, "series ", 7)) {
	p->targ = VEC;
	s += 7;
    } else if (!strncmp(s, "matrix *", 8)) {
	p->targ = MAT;
	p->flags |= P_LHPTR;
	s += 8;
	s += strspn(s, " ");
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
    }

    if (p->targ == VEC && p->dset->n == 0) {
	no_data_error(p);
	return;
    }

    if (p->flags & P_DISCARD) {
	/* doing a simple "eval" */
	p->point = s;
	return;
    }

    /* LHS varname (possibly with substring) */
    p->err = extract_LHS_string(s, test, p);
    if (p->err) {
	return;
    } 

    if (p->flags & P_DECL) {
	return;
    }

    /* record next read position */
    p->point = s + strlen(test);

    /* grab LHS obs string, matrix slice, or bundle element, 
       if present */
    if (strchr(test, '[') != NULL) {
	subchar = '[';
	get_lhs_substr(test, p);
    } else if (strchr(test, '.') != NULL) {
	subchar = '.';
	get_lhs_dotstr(test, p);
    }

    if (p->err) {
	return;
    }    

#if LHDEBUG || EDEBUG
    fprintf(stderr, "LHS: %s", test);
    if (p->lh.substr != NULL) {
	fprintf(stderr, " substr: '%s'\n", p->lh.substr);
    } else {
	fputc('\n', stderr);
    }
#endif

    if (strlen(test) > VNAMELEN - 1) {
	pprintf(p->prn, _("'%s': name is too long (max %d characters)\n"), 
		test, VNAMELEN - 1);
	p->err = E_DATA;
	return;
    }

    /* find out if the LHS var already exists, and if
       so, what type it is */
    if ((v = current_series_index(p->dset, test)) >= 0) {
	p->lh.v = v;
	p->lh.t = VEC;
	newvar = 0;
    } else {
	GretlType vtype = user_var_get_type_by_name(test);

	if (vtype == GRETL_TYPE_MATRIX) {
	    p->lh.m0 = get_matrix_by_name(test);
	    p->flags |= P_LHMAT;
	    p->lh.t = MAT;
	    newvar = 0;
	} else if (vtype == GRETL_TYPE_DOUBLE) {
	    p->flags |= P_LHSCAL;
	    p->lh.t = NUM;
	    newvar = 0;
	} else if (vtype == GRETL_TYPE_LIST) {
	    p->flags |= P_LHLIST;
	    p->lh.t = LIST;
	    newvar = 0;
	} else if (vtype == GRETL_TYPE_STRING) {
	    p->flags |= P_LHSTR;
	    p->lh.t = STR;
	    newvar = 0;
	} else if (vtype == GRETL_TYPE_BUNDLE) {
	    p->flags |= P_LHBUN;
	    p->lh.t = BUNDLE;
	    newvar = 0;
	} 
    }

    /* if pre-existing var, check for const-ness */
    if (!newvar && overwrite_const_check(test, p)) {
	return;
    }

    /* if new public variable, check name for legality */
    if (newvar) {
	if (flags & P_PRIVATE) {
	    p->err = check_private_varname(test);
	} else {
	    p->err = check_varname(test);
	}
	if (p->err) {
	    return;
	}
    }

    if (p->lh.substr != NULL) {
	process_lhs_substr(test, subchar, p);
	if (p->err) {
	    return;
	}
    }

    strcpy(p->lh.name, test);

    if (p->lh.t != UNK) {
	if (p->targ == UNK) {
	    /* when a type is not specified, set from existing
	       variable, if present */
	    p->targ = p->lh.t;
	} else if (overwrite_type_check(p)) {
	    /* don't overwrite one type with another */
	    p->err = E_TYPES;
	    return;
	}
    }

    /* advance past varname */
    s = p->point;
    while (isspace(*s)) s++;

    /* expression ends here: a call to print? */
    if (*s == '\0' || !strcmp(s, "print")) {
	parser_try_print(p);
	return;
    }

    /* operator: '=' or '+=' etc. */
    strncat(opstr, s, 2);
    if ((p->op = get_op(opstr)) == 0) {
	p->err = E_EQN;
	return;
    }

    /* if the LHS variable does not already exist, then
       we can't do '+=' or anything of that sort, only 
       simple assignment, B_ASN
    */
    if (newvar && p->op != B_ASN) {
	/* error message */
	pprintf(p->prn, "%s: unknown variable\n", test);
	p->err = E_UNKVAR;
	return;
    }

    /* matrices: we accept only a limited range of
       modified assignment operators */
    if (p->lh.t == MAT && !ok_matrix_op(p->op)) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for matrices"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* lists: same story as matrices */
    if (p->lh.t == LIST && !ok_list_op(p->op)) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for lists"), opstr);
	p->err = E_PARSE;
	return;
    }	

    /* strings: ditto */
    if (p->lh.t == STR && !ok_string_op(p->op)) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for strings"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* bundles: we can't do any sort of modified assignment (yet) */
    if (p->lh.t == BUNDLE && p->op != B_ASN) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* vertical concat: only OK for matrices */
    if (p->lh.t != MAT && (p->op == B_VCAT || p->op == B_DOTASN)) {
	gretl_errmsg_sprintf(_("'%s' : only defined for matrices"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* horizontal concat: only OK for matrices, strings */
    if (p->lh.t != MAT && p->lh.t != STR && p->op == B_HCAT) {
	gretl_errmsg_sprintf(_("'%s' : not implemented for this type"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* string-valued series: do not overwrite wholesale */
    if (p->lh.t == VEC && p->lh.obs < 0 && 
	is_string_valued(p->dset, p->lh.v)) {
	gretl_errmsg_set("Cannot overwrite entire string-valued series");
	p->err = E_TYPES;
	return;
    }

    /* advance past operator */
    s += strlen(opstr);
    while (isspace(*s)) s++;

    /* set starting point for RHS parser, and also
       for a possible label */
    p->point = p->rhs = s;

    /* if the target type is still unknown, and the RHS expression
       is wrapped in '{' and '}', make the target a matrix */
    if (p->targ == UNK && *p->rhs == '{') {
	maybe_set_matrix_target(p);
    }

    /* unary increment/decrement operators */
    if ((p->op == INC || p->op == DEC) && *s != '\0') {
	p->err = E_PARSE;
    }

    /* pointer for use in label */
    if (p->op == B_ASN) {
	p->rhs = s;
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

static int series_compatible (const gretl_matrix *m, parser *p)
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

/* below: note that node we're checking is p->ret, which may differ in
   type from the target to which it will be assigned/converted
*/

static void gen_check_errvals (parser *p)
{
    NODE *n = p->ret;

    if (n == NULL || (n->t == VEC && n->v.xvec == NULL)) {
	return;
    }

    if (p->targ == MAT) {
	/* the matrix target is handled separately */
	return;
    }

    if (n->t == NUM) {
	if (!isfinite(n->v.xval)) {
#if SCALARS_ENSURE_FINITE
	    n->v.xval = NADBL;
	    set_gretl_warning(W_GENMISS);
#else
	    set_gretl_warning(W_GENNAN);
#endif
	}
    } else if (n->t == VEC) {
	int t;

	for (t=p->dset->t1; t<=p->dset->t2; t++) {
	    if (!isfinite(n->v.xvec[t])) {
#if SERIES_ENSURE_FINITE
		n->v.xvec[t] = NADBL;
		set_gretl_warning(W_GENMISS);
#else
		set_gretl_warning(W_GENNAN);
#endif
		break;
	    }
	}
    } else if (n->t == MAT) {
	const gretl_matrix *m = n->v.m;
	int i, k = gretl_matrix_rows(m) * gretl_matrix_cols(m);
	
	if (p->targ == NUM && k == 1) {
	    if (!isfinite(m->val[0])) {
#if SCALARS_ENSURE_FINITE
		m->val[0] = NADBL;
		set_gretl_warning(W_GENMISS);
#else
		set_gretl_warning(W_GENNAN);
#endif
	    }		
	} else {
	    /* convert any NAs to NaNs */
	    for (i=0; i<k; i++) {
		if (na(m->val[i])) {
		    m->val[i] = M_NA;
		    set_gretl_warning(W_GENNAN);
		} else if (!isfinite(m->val[i])) {
		    set_gretl_warning(W_GENNAN);
		}
	    }
	}
    }
}

static gretl_matrix *list_to_matrix (const int *list, int *err)
{
    gretl_matrix *v = NULL;

    if (list == NULL) {
	*err = E_DATA;
    } else {
	int i, n = list[0];

	if (n == 0) {
	    v = gretl_null_matrix_new();
	    if (v == NULL) {
		*err = E_ALLOC;
	    }	    
	} else if (n > 0) {
	    v = gretl_vector_alloc(n);
	    if (v == NULL) {
		*err = E_ALLOC;
	    } else {
		for (i=0; i<n; i++) {
		    v->val[i] = list[i+1];
		}
	    }
	} else {
	    *err = E_DATA;
	}
    }

    return v;
}

static gretl_matrix *series_to_matrix (const double *x, 
				       parser *p)
{
    gretl_matrix *v = NULL;
    int i, n = sample_size(p->dset);

    v = gretl_column_vector_alloc(n);
    if (v == NULL) {
	p->err = E_ALLOC;
    } else {
	for (i=0; i<n; i++) {
	    v->val[i] = x[i + p->dset->t1];
	}
    }

    return v;
}

static gretl_matrix *grab_or_copy_matrix_result (parser *p,
						 int *prechecked)
{
    NODE *r = p->ret;
    gretl_matrix *m = NULL;

#if EDEBUG
    fprintf(stderr, "grab_or_copy_matrix_result: r->t = %d\n", r->t);
#endif

    if (r->t == NUM) {
	m = gretl_matrix_from_scalar(r->v.xval);
    } else if (r->t == VEC) {
	m = series_to_matrix(r->v.xvec, p);
    } else if (r->t == LIST) {
	m = list_to_matrix(r->v.ivec, &p->err);
    } else if (r->t == MAT && is_tmp_node(r)) {
	/* result matrix is newly allocated, steal it */
#if EDEBUG
	fprintf(stderr, "matrix result (%p) is tmp, stealing it\n", 
		(void *) r->v.m);
#endif
	m = r->v.m;
	r->v.m = NULL; /* avoid double-freeing */
    } else if (r->t == MAT) {
	/* r->v.m is an existing user matrix (or bundled matrix) */
	if (p->flags & P_LHPTR) {
	    if (r->flags & PTR_NODE) {
		/* OK, we'll share the matrix pointer */
		m = r->v.m;
	    } else {
		p->err = E_TYPES;
	    }
	} else {
	    /* must make a copy to keep pointers distinct */
	    m = gretl_matrix_copy(r->v.m);
#if EDEBUG
	    fprintf(stderr, "matrix result (%p) is pre-existing, copied to %p\n",
		    (void *) r->v.m, (void *) m);
#endif
	    if (m == NULL) {
		p->err = E_ALLOC;
	    }
	}
	if (p->err == 0 && prechecked != NULL) {
	    *prechecked = 1;
	}
    } else if (r->t == BUNDLE) {
	/* can the bundle be "cast" to a matrix? */
	m = gretl_bundle_get_payload_matrix(r->v.b, &p->err);
	if (p->err == 0 && prechecked != NULL) {
	    *prechecked = 1;
	}	
    } else {
	fprintf(stderr, "Looking for matrix, but r->t = %d\n", r->t);
	p->err = E_TYPES;
    }

    return m;
}

/* generating a matrix, no pre-existing matrix of that name */

static gretl_matrix *matrix_from_scratch (parser *p, int tmp,
					  int *prechecked)
{
    gretl_matrix *m = NULL;

#if EDEBUG
    fprintf(stderr, "matrix_from_scratch: reusable = %d\n",
	    reusable(p));
#endif

    if (p->ret->t == NUM) {
	m = gretl_matrix_from_scalar(p->ret->v.xval);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} 
    } else {
	m = grab_or_copy_matrix_result(p, prechecked);
    }

    if (!tmp && !p->err) {
	int adj = p->lh.m0 != NULL;

	p->err = user_var_add(p->lh.name, GRETL_TYPE_MATRIX, m);
	if (adj) {
	    p->lh.m0 = m; /* make the lh matrix pointer valid */
	}
    }

    p->lh.m1 = m;
    
    return m;
}

static int LHS_matrix_reusable (parser *p)
{
    gretl_matrix *m = p->lh.m0;
    int ok = 0;

    if (p->ret->t == NUM) {
	ok = (m->rows == 1 && m->cols == 1);
    } else if (p->ret->t == VEC) {
	int T = sample_size(p->dset);

	ok = (m->rows == T && m->cols == 1);
    } else if (p->ret->t == MAT) {
	gretl_matrix *rm = p->ret->v.m;

	ok = (rm != NULL &&
	      m->rows == rm->rows && 
	      m->cols == rm->cols);
    }

    return ok;
}

/* Generating a matrix, and there's a pre-existing LHS matrix:
   we re-use the left-hand side matrix if possible.
*/

static void assign_to_matrix (parser *p, int *prechecked)
{
    gretl_matrix *m;
    double x;

    if (LHS_matrix_reusable(p)) {
	/* result is conformable with original matrix */
#if EDEBUG
	fprintf(stderr, "assign_to_matrix: reusing LHS\n");
#endif
	m = p->lh.m0;
	if (p->ret->t == NUM) {
	    x = p->ret->v.xval;
	    m->val[0] = na(x)? M_NA : x;
	} else if (p->ret->t == VEC) {
	    int i, s = p->dset->t1;

	    for (i=0; i<m->rows; i++) {
		x = p->ret->v.xvec[s++];
		m->val[i] = na(x)? M_NA : x;
	    }
	} else {
	    p->err = gretl_matrix_copy_data(m, p->ret->v.m);
	}
    } else {
	/* replace the old matrix with result */
#if EDEBUG
	fprintf(stderr, "assign_to_matrix: replacing\n");
#endif
	m = grab_or_copy_matrix_result(p, prechecked);
	if (!p->err) {
	    p->err = user_matrix_replace_matrix_by_name(p->lh.name, m);
	    p->lh.m0 = NULL; /* invalidate pointer */
	}
    }

    p->lh.m1 = m;
}

/* assigning to an existing (whole) LHS matrix, but using '+=' or
   some such modified/inflected assignment */

static void assign_to_matrix_mod (parser *p)
{
    user_var *uvar;
    gretl_matrix *m = NULL;
    gretl_matrix *a, *b;

    uvar = get_user_var_of_type_by_name(p->lh.name, GRETL_TYPE_MATRIX);

    if (uvar == NULL) {
	p->err = E_DATA;
    } else {
	a = user_var_get_value(uvar);
	if (a == NULL) {
	    p->err = E_DATA;
	}
    }

    if (!p->err) {
	if (p->op == B_DOTASN) {
	    if (p->ret->t == NUM) {
		gretl_matrix_fill(a, p->ret->v.xval);
		return;
	    } else {
		p->err = E_TYPES;
	    }
	} else {
	    b = matrix_from_scratch(p, 1, NULL);
	    if (b != NULL) {
		m = real_matrix_calc(a, b, p->op, &p->err);
		gretl_matrix_free(b);
	    }
	}
    }
 
    p->lh.m1 = m;

    if (!p->err) {
	p->err = user_var_replace_value(uvar, m);
    }
}

/* Here we're replacing a sub-matrix of the original LHS matrix, by
   either straight or inflected assignment. The value that we're
   using for replacement will be either a matrix or a scalar.
*/

static void edit_matrix (parser *p)
{
    matrix_subspec *spec;
    gretl_matrix *m = NULL;

    if (p->ret->t != NUM) {
	/* not a scalar: get the replacement matrix */
	m = grab_or_copy_matrix_result(p, NULL);
	if (m == NULL) {
	    return;
	}
    }

#if EDEBUG
    fprintf(stderr, "edit_matrix: m = %p\n", (void *) m);
#endif

    spec = p->lh.mspec;

    /* check the validity of the subspec we got */
    p->err = check_matrix_subspec(spec, p->lh.m0);
    if (p->err) {
	return;
    }

    if (p->ret->t == NUM && spec->type[0] == SEL_ELEMENT) {
	/* Assignment (possibly "inflected") of a scalar value 
	   to a single element of an existing matrix
	*/
	int i = mspec_get_row_index(spec);
	int j = mspec_get_col_index(spec);
	double x = matrix_get_element(p->lh.m0, i, j, &p->err);

	if (!p->err) {
	    x = xy_calc(x, p->ret->v.xval, p->op, MAT, p);
	    if (xna(x)) {
		if (na(x)) {
		    x = M_NA;
		}
		set_gretl_warning(W_GENNAN);
	    }
	    gretl_matrix_set(p->lh.m0, i-1, j-1, x);
	    /* Flag the fact that we produced a matrix, even
	       though it's the one that was present on input.
	    */
	    p->lh.m1 = p->lh.m0;
	}
	return; /* note, we're done */
    } 

    if (p->ret->t == NUM && p->op == B_ASN) {
	/* Straight assignment of a scalar value to non-scalar
	   submatrix */
	double x = p->ret->v.xval;

	if (xna(x)) {
	    if (na(x)) {
		x = M_NA;
	    }
	    set_gretl_warning(W_GENNAN);
	}
	p->err = assign_scalar_to_submatrix(p->lh.m0, x, spec);
	p->lh.m1 = p->lh.m0;
	return; /* note, we're done */
    }

    if (p->op != B_ASN) {
	/* Here we're doing '+=' or some such, in which case a new
	   submatrix must be calculated using the original
	   submatrix 'a' and the newly generated matrix (or
	   scalar value).
	*/
	gretl_matrix *a = matrix_get_submatrix(p->lh.m0, spec, 
					       1, &p->err);

	if (!p->err) {
	    if (p->ret->t == NUM) {
		int i, n = a->rows * a->cols;
		double x = p->ret->v.xval;

		for (i=0; i<n; i++) {
		    a->val[i] = xy_calc(a->val[i], x, p->op, MAT, p);
		}
		/* assign computed matrix to m */
		m = a;
	    } else {
		gretl_matrix *b = real_matrix_calc(a, m, p->op, &p->err);

		gretl_matrix_free(a);
		/* replace existing m with computed result */
		gretl_matrix_free(m);
		m = b;
	    }
	}
    } 

    if (!p->err) {
	/* Write new submatrix m into place: note that we come here
	   directly if none of the special conditions above are
	   satisfied -- for example, if the newly generated value
	   is a matrix and the task is straight assignment. Also
	   check for numerical "breakage" in the replacement
	   submatrix.
	*/
	p->err = user_matrix_replace_submatrix(p->lh.name, m, spec);
	if (!p->err && gretl_matrix_xna_check(m)) {
	    set_gretl_warning(W_GENNAN);
	}
	gretl_matrix_free(m);
	if (p->ret->t == MAT) {
	    p->ret->v.m = NULL; /* ?? */
	}
	p->lh.m1 = get_matrix_by_name(p->lh.name);
    }
}

static int edit_string (parser *p)
{
    const char *src = NULL;
    const char *orig = NULL;
    char *newstr = NULL;
    user_var *uvar;

    if (p->ret->t == NUM && p->op != B_ADD) {
	p->err = E_TYPES;
	return p->err;
    }

    if (p->ret->t == EMPTY) {
	src = "";
    } else {
	src = p->ret->v.str;
    }

#if EDEBUG
    fprintf(stderr, "edit_string: src='%s'\n", src);
#endif

    uvar = get_user_var_of_type_by_name(p->lh.name, GRETL_TYPE_STRING);

    if (uvar != NULL) {
	orig = user_var_get_value(uvar);
    } else if (p->op != B_ASN) {
	/* without a left-hand string we can only assign */
	p->err = E_DATA;
	return p->err;
    }

    if (p->ret->t == NUM) {
	/* taking an offset into an existing string ("+=") */
	int len = strlen(orig);
	int adj = p->ret->v.xval;
		
	if (adj < 0) {
	    p->err = E_DATA;
	} else if (adj == 0) {
	    ; /* no-op */
	} else {
	    src = (adj < len)? (orig + adj) : "";
	    newstr = gretl_strdup(src);
	    if (newstr == NULL) {
		p->err = E_ALLOC;
	    } else {
		user_var_replace_value(uvar, newstr);
	    }
	}
    } else if (src == NULL) {
	; /* no-op -- e.g. argname() didn't get anything */
    } else if (p->op == B_ASN) {
	/* simple assignment */
	newstr = gretl_strdup(src);
	if (newstr == NULL) {
	    p->err = E_ALLOC;
	} else if (uvar == NULL) {
	    user_var_add(p->lh.name, GRETL_TYPE_STRING, newstr);
	} else {
	    user_var_replace_value(uvar, newstr);
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
		user_var_replace_value(uvar, newstr);
	    }
	}
    } else if (p->op == INC) {
	/* string++ */
	size_t len = strlen(orig);

	if (len < 2) {
	    newstr = gretl_strdup("");
	} else {
	    newstr = malloc(len);
	}
	if (newstr == NULL) {
	    p->err = E_ALLOC;
	} else {
	    *newstr = '\0';
	    strncat(newstr, orig + 1, len - 1);
	    user_var_replace_value(uvar, newstr); 
	}
    }

    return p->err;
}

static int edit_list (parser *p)
{
    int *list = node_get_list(p->ret, p); /* copied */

#if EDEBUG
    printlist(list, "incoming list in edit_list()");
#endif

    if (!p->err) {
	if (!lhlist(p)) {
	    /* no pre-existing LHS list: must be simple assignment */
	    p->err = remember_list(list, p->lh.name, NULL);
	} else if (p->op == B_ASN) {
	    /* assign to (i.e. replace) existing LHS list */
	    p->err = replace_list_by_name(p->lh.name, list);
	} else if (p->op == B_ADD) {
	    /* add to existing LHS list */
	    p->err = append_to_list_by_name(p->lh.name, list);
	} else if (p->op == B_SUB) {
	    /* remove elements from existing LHS list */
	    p->err = subtract_from_list_by_name(p->lh.name, list);
	} else {
	    p->err = E_TYPES;
	}
    }

    free(list);

    return p->err;
}

#define ok_return_type(t) (t == NUM || t == VEC || t == MAT || \
			   t == LIST || t == DUM || t == EMPTY || \
                           t == STR || t == BUNDLE || t == U_ADDR)

static int gen_check_return_type (parser *p)
{
    NODE *r = p->ret;

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: targ = %d; ret at %p, type %d\n", 
	    p->targ, (void *) r, (r == NULL)? -999 : r->t);
#endif

    if (r == NULL) {
	fprintf(stderr, "gen_check_return_type: p->ret = NULL!\n");
	return (p->err = E_DATA);
    }

    if ((p->dset == NULL || p->dset->n == 0) && 
	r->t != MAT && r->t != NUM && 
	r->t != STR && r->t != BUNDLE &&
	r->t != U_ADDR && r->t != EMPTY) {
	no_data_error(p);
	return p->err;
    }

    if (!ok_return_type(r->t)) {
	return (p->err = E_TYPES);
    }

    if (r->t == VEC && r->v.xvec == NULL) {
	fprintf(stderr, "got VEC return with xvec == NULL!\n");
	return (p->err = E_DATA);
    }

    if (p->targ == NUM) {
	if (r->t == NUM || scalar_matrix_node(r)) {
	    ; /* scalar or 1 x 1 matrix: OK */
	} else if (r->t == STR && p->lh.obs >= 0) {
	    ; /* a string value might be acceptable */
	} else {
	    p->err = E_TYPES;
	}
    } else if (p->targ == VEC) {
	/* result must be scalar, series, or conformable matrix */
	if (r->t == NUM || r->t == VEC) {
	    ; /* OK */
	} else if (r->t == MAT) {
	    if (!series_compatible(r->v.m, p)) {
		p->err = E_TYPES;
	    }
	} else {
	    p->err = E_TYPES;
	}
    } else if (p->targ == MAT) {
	; /* no-op: leave targ alone */
    } else if (p->targ == LIST) {
	if (r->t != EMPTY && !ok_list_node(r)) {
	    p->err = E_TYPES;
	} 
    } else if (p->targ == STR) {
	if (r->t != EMPTY && r->t != STR && r->t != NUM) {
	    p->err = E_TYPES;
	}
    } else if (p->targ == BUNDLE) {
	if (r->t != BUNDLE && r->t != EMPTY) {
	    p->err = E_TYPES;
	}
    } else if (p->targ == BOBJ || p->targ == BMEMB) {
	if (!ok_bundled_type(r->t)) {
	    p->err = E_TYPES;
	}
    } else {
	/* target type was not specified: set it now, based
	   on the type of the object we computed */
	if (r->t == LIST) {
	    p->targ = LIST;
	} else if (scalar_matrix_node(r)) {
	    /* cast a 1 x 1 matrix to a scalar */
	    p->targ = NUM;
	} else {
	    p->targ = r->t;
	}
    }

    if (p->err == E_TYPES) {
	maybe_do_type_errmsg(p->lh.name, p->lh.t);
    }	

#if EDEBUG
    fprintf(stderr, "gen_check_return_type: returning with p->err = %d\n", 
	    p->err);
#endif

    return p->err;
}

/* allocate storage if saving a series to the dataset: 
   lh.v == 0 means that the LHS variable does not already 
   exist
*/

static int gen_allocate_storage (parser *p)
{
    if (p->lh.v == 0) {
	if (p->dset == NULL || p->dset->Z == NULL) {
	    p->err = E_DATA;
	} else {
	    p->err = dataset_add_series(p->dset, 1);
	}
	if (!p->err) {
	    int t;

	    p->lh.v = p->dset->v - 1;
	    for (t=0; t<p->dset->n; t++) {
		p->dset->Z[p->lh.v][t] = NADBL;
	    }
#if EDEBUG
	    fprintf(stderr, "gen_allocate_storage: added series #%d (%s)\n",
		    p->lh.v, p->lh.name);
#endif
	}
    }

    return p->err;
}

#if SERIES_ENSURE_FINITE

static void series_make_finite (double *x, int n)
{
    int i;

    /* This may be questionable (converting NaNs and infinities to
       NA), but is necessary for backward compatibility (e.g. when
       taking logs of a series that contains zeros).
    */

    for (i=0; i<n; i++) {
	if (!isfinite(x[i])) {
	    x[i] = NADBL;
	    set_gretl_warning(W_GENMISS);
	}
    }
}

#endif

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
		y[t] = xy_calc(y[t], m->val[s], p->op, VEC, p);
	    }
	    s++;
	}
    }
}

static int assign_null_to_bundle (parser *p)
{
    gretl_bundle *b;
    int err = 0;

    if (p->flags & P_LHBUN) {
	b = get_bundle_by_name(p->lh.name);
	gretl_bundle_void_content(b);
    } else {
	b = gretl_bundle_new();
	if (b == NULL) {
	    err = E_ALLOC;
	} else {
	    err = user_var_add_or_replace(p->lh.name,
					  GRETL_TYPE_BUNDLE,
					  b);
	}
    }

    return err;
}

static int save_generated_var (parser *p, PRN *prn)
{
    NODE *r = p->ret;
    double **Z = NULL;
    double x;
    int t, v = 0;

    /* test for type mismatch errors */
    gen_check_return_type(p);
    if (p->err) {
	return p->err;
    }

#if EDEBUG
    fputs("save_generated_var:\n ", stderr);
    fprintf(stderr, "targ = %d (%s), ret = %d (%s), op = %d (%s)\n",
	    p->targ, getsymb(p->targ, NULL), 
	    p->ret->t, getsymb(p->ret->t, NULL), 
	    p->op, getsymb(p->op, NULL));
#endif

    /* allocate dataset storage, if needed */
    if (p->targ == VEC) {
	gen_allocate_storage(p);
	if (p->err) {
	    return p->err;
	}
    }

    if (p->dset != NULL && p->dset->Z != NULL) {
	Z = p->dset->Z;
	v = p->lh.v;
    }

    /* put the generated data into place */
    
    if (p->targ == NUM) {
	/* writing a scalar */
	if (p->lh.obs >= 0) {
	    /* target is actually a specific observation in a series */
	    t = p->lh.obs;
	    if (is_string_valued(p->dset, v)) {
		if (r->t == STR) {
		    p->err = series_set_string_val(p->dset, v, t, r->v.str);
		} else {
		    gretl_errmsg_sprintf(_("type error: %s is string-valued"), 
					 p->lh.name);
		    p->err = E_TYPES;
		}
	    } else if (r->t == NUM) {
		Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, NUM, p);
	    } else if (r->t == MAT) {
		Z[v][t] = xy_calc(Z[v][t], r->v.m->val[0], p->op, NUM, p);
	    }
	    if (p->err == 0) {
		strcpy(p->dset->varname[v], p->lh.name);
		set_dataset_is_changed();
	    }
	} else if (p->flags & P_LHSCAL) {
	    /* modifying existing scalar */
	    x = gretl_scalar_get_value(p->lh.name, NULL);
	    if (r->t == NUM) {
		x = xy_calc(x, r->v.xval, p->op, NUM, p);
	    } else if (scalar_matrix_node(r)) {
		x = xy_calc(x, r->v.m->val[0], p->op, NUM, p);
	    } else {
		p->err = E_TYPES;
	    }
	    if (!p->err) {
		p->err = gretl_scalar_set_value(p->lh.name, x);
	    }
	} else {
	    /* a new scalar */
	    x = (r->t == MAT)? r->v.m->val[0] : r->v.xval;
	    p->err = gretl_scalar_add(p->lh.name, x);
	}
    } else if (p->targ == VEC) {
	/* writing a series */
	if (r->t == NUM) {
	    for (t=p->dset->t1; t<=p->dset->t2; t++) { 
		Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, VEC, p);
	    }
	} else if (r->t == VEC) {
	    const double *x = r->v.xvec;
	    int t1 = p->dset->t1;

	    if (autoreg(p) && p->op == B_ASN) {
		while (xna(x[t1]) && t1 <= p->dset->t2) {
		    t1++;
		}
	    }
	    if (p->op == B_ASN) {
		/* avoid multiple calls to xy_calc */
		if (Z[v] != x) {
		    size_t sz = (p->dset->t2 - t1 + 1) * sizeof *x;

		    memcpy(Z[v] + t1, x + t1, sz);
		}
	    } else {
		for (t=t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], x[t], p->op, VEC, p);
		}
	    }
	} else if (r->t == MAT) {
	    const gretl_matrix *m = r->v.m;
	    int k = gretl_vector_get_length(m);
	    int mt1 = gretl_matrix_get_t1(m);
	    int s;

	    if (p->flags & P_MMASK) {
		/* result needs special alignment */
		align_matrix_to_series(Z[v], m, p);
	    } else if (k == 1) {
		/* result is effectively a scalar */
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[0], p->op, VEC, p);
		}
	    } else if (k == p->dset->n) {
		/* treat result as full-length series */
		for (t=p->dset->t1; t<=p->dset->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[t], p->op, VEC, p);
		}
	    } else if (k == sample_size(p->dset) && mt1 == 0) {
		/* treat as series of current sample length */
		for (t=p->dset->t1, s=0; t<=p->dset->t2; t++, s++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[s], p->op, VEC, p);
		}
	    } else {
		/* align using matrix "t1" value */
		for (t=mt1; t<mt1 + k && t<=p->dset->t2; t++) {
		    if (t >= p->dset->t1) {
			Z[v][t] = xy_calc(Z[v][t], m->val[t - mt1], p->op, VEC, p);
		    }
		}
	    }
	}

#if SERIES_ENSURE_FINITE
	if (!p->err) {
	    series_make_finite(Z[v], p->dset->n);
	}
#endif
	strcpy(p->dset->varname[v], p->lh.name);
#if EDEBUG
	fprintf(stderr, "var %d: gave generated series the name '%s'\n", 
		v, p->lh.name);
#endif
	if (!p->err) {
	    set_dataset_is_changed();
	}
    } else if (p->targ == MAT) {
	/* we're writing a matrix */
	int prechecked = 0;

	if (p->lh.m0 == NULL) {
	    /* there's no pre-existing left-hand side matrix */
	    matrix_from_scratch(p, 0, &prechecked);
	} else if (p->lh.substr == NULL && p->op == B_ASN) {
	    /* uninflected assignment to an existing matrix */
	    assign_to_matrix(p, &prechecked);
	} else if (p->lh.substr == NULL) {
	    /* inflected assignment to entire existing matrix */
	    assign_to_matrix_mod(p);
	} else {
	    /* assignment to submatrix of original */
	    edit_matrix(p);
	    prechecked = 1;
	}
	if (!prechecked && gretl_matrix_xna_check(p->lh.m1)) {
	    set_gretl_warning(W_GENNAN);
	}
    } else if (p->targ == LIST) {
	edit_list(p);
    } else if (p->targ == STR) {
	edit_string(p);
    } else if (p->targ == BUNDLE) {
	if ((r->flags & TMP_NODE) || (p->flags & P_UFRET)) {
	    /* bundle created on the fly */
	    p->err = user_var_add_or_replace(p->lh.name,
					     GRETL_TYPE_BUNDLE,
					     r->v.b);
	    if (!p->err) {
		/* avoid destroying the returned bundle */
		r->v.b = NULL;
	    }
	} else if (r->t == EMPTY) {
	    /* as in "bundle b = null" */
	    p->err = assign_null_to_bundle(p);
	} else {
	    /* assignment from pre-existing named bundle */
	    p->err = gretl_bundle_copy_as(p->rhs, p->lh.name);
	}
    } else if (p->targ == BOBJ || p->targ == BMEMB) {
	/* saving an object into a bundle */
	p->err = set_named_bundle_value(p->lh.name, r, p);
    }

#if EDEBUG
    if (!p->err) {
	parser_print_result(p, prn);
    } else {
	fprintf(stderr, "save_generated_var: returning p->err = %d\n",
		p->err);
    }
#endif

    return p->err;
}

static void parser_reinit (parser *p, DATASET *dset, PRN *prn) 
{
    /* flags that should be reinstated if they were
       present at compile time */
    int repflags[] = { P_PRINT, P_NATEST, P_AUTOREG,
		       P_SLAVE, P_SLICE, P_VOID,
		       P_LHPTR, P_DISCARD, P_SCALAR, 
		       0 };
    int i, saveflags = p->flags;

    /* P_LHSCAL, P_LHLIST, P_LHSTR ? */

    p->flags = (P_START | P_PRIVATE | P_EXEC);

    for (i=0; repflags[i] > 0; i++) {
	if (saveflags & repflags[i]) {
	    p->flags |= repflags[i];
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
    p->uval = NULL;

    p->ret = NULL;
    p->err = 0;

#if EDEBUG
    fprintf(stderr, "parser_reinit: p->subp=%p, targ=%d, lhname='%s', op=%d\n", 
	    (void *) p->subp, p->targ, p->lh.name, p->op);
#endif

    if (p->targ == MAT) {
	/* matrix target: check the LH name again */
	if (*p->lh.name != '\0') {
	    p->lh.m0 = get_matrix_by_name(p->lh.name);
	}
    } else if (p->targ == NUM) {
	/* scalar target, also check LH name */
	if (*p->lh.name != '\0' && gretl_is_scalar(p->lh.name)) {
	    p->flags |= P_LHSCAL;
	}
	if (p->lh.substr != NULL) {
	    /* reevaluate obs string */
	    process_lhs_substr(NULL, 0, p);
	}
    } else if (p->targ == LIST) {
	/* list target, check LH name */
	if (*p->lh.name != '\0' && get_list_by_name(p->lh.name)) {
	    p->flags |= P_LHLIST;
	}
    } else if (p->targ == VEC) {
	if (p->lh.v >= p->dset->v) {
	    /* recorded series ID is no longer valid */
	    p->lh.v = 0;
	}
    } else if (p->targ == STR) {
	/* FIXME is this wanted */
	if (*p->lh.name != '\0' && get_string_by_name(p->lh.name)) {
	    p->flags |= P_LHSTR;
	}
    }	

    /* LHS matrix subspec: re-evaluate */
    if (p->subp != NULL) {
	get_lh_mspec(p);
    }    
}

static void parser_init (parser *p, const char *str, 
			 DATASET *dset, PRN *prn, 
			 int flags)
{
    p->point = p->rhs = p->input = str;
    p->dset = dset;
    p->prn = prn;
    p->flags = flags | P_START;
    p->targ = UNK;
    p->op = 0;

    p->tree = NULL;
    p->ret = NULL;

    p->lh.t = UNK;
    p->lh.name[0] = '\0';
    p->lh.label[0] = '\0';
    p->lh.v = 0;
    p->lh.obs = -1;
    p->lh.m0 = NULL;
    p->lh.m1 = NULL;
    p->lh.substr = NULL;
    p->lh.mspec = NULL;

    p->subp = NULL;

    p->obs = 0;
    p->sym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->err = 0;

    if (p->input == NULL) {
	p->err = E_DATA;
	return;
    }

    if (p->flags & P_SLICE) {
	p->lh.t = MAT;
    } else if (p->flags & P_SCALAR) {
	p->targ = NUM;
    } else if (p->flags & P_SERIES) {
	p->targ = VEC;
    } else if (p->flags & P_MATRIX) {
	p->targ = MAT;
    } else if (p->flags & P_STRING) {
	p->targ = STR;
    } else if (p->flags & P_LIST) {
	p->targ = LIST;
    } else if (p->flags & P_VOID) {
	p->targ = EMPTY;
    } else {
	pre_process(p, flags);
    }

    if (!p->err) {
	p->ch = parser_getc(p);
    }
}

/* called from genmain.c (only!) */

void gen_save_or_print (parser *p, PRN *prn)
{
    if (p->err == 0) {
	if (p->flags & (P_DISCARD | P_PRINT)) {
	    if (p->ret->t == MAT) {
		gretl_matrix_print_to_prn(p->ret->v.m, p->lh.name, p->prn);
	    } else if (p->ret->t == LIST) {
		gretl_list_print(p->lh.name, p->dset, p->prn);
	    } else if (p->ret->t == STR) {
		if (p->lh.name[0] != '\0') {
		    pprintf(p->prn, "%s\n", get_string_by_name(p->lh.name));
		} else {
		    pprintf(p->prn, "%s\n", p->ret->v.str);
		}
	    } else if (p->ret->t == BUNDLE) {
		gretl_bundle_print(get_bundle_by_name(p->lh.name), prn);
	    } else {
		printnode(p->ret, p);
		pputc(p->prn, '\n');
	    }
	} else if (p->flags & P_DECL) {
	    do_decl(p);
	} else {
	    save_generated_var(p, prn);
	} 
    }

#if 0
    if (p->targ == MAT) {
	fprintf(stderr, "genr exec (%s): "
		" m0 = %p, m1 = %p\n", p->lh.name, (void *) p->lh.m0, 
		(void *) p->lh.m1);
    } else if (p->lh.v == 0) {
	fprintf(stderr, "genr exec (%s): p->lh.v = %d\n", p->lh.name,
		p->lh.v);
    }
#endif    
}

void gen_cleanup (parser *p)
{
    int protect = reusable(p);

#if EDEBUG
    fprintf(stderr, "gen cleanup on %p: reusable = %d, err = %d\n", 
	    p, reusable(p) ? 1 : 0, p->err);
#endif

    if (p->err && (p->flags & P_COMPILE)) {
	protect = 0;
    }

    if (protect) {
	/* just do limited cleanup */
	if (!in_tree(p->tree, p->ret)) {
	    free_tree(p->ret, p, "p->ret");
	    p->ret = NULL;
	}
    } else {
	if (p->ret != p->tree) {
	    free_tree(p->tree, p, "p->tree");
	}

	free_tree(p->ret, p, "p->ret");
	free(p->lh.substr);
	free_mspec(p->lh.mspec);

	if (p->subp != NULL) {
	    /* since the parent genr will not be run again,
	       wipe the "reusable" flags from its subp */
	    p->subp->flags &= ~P_COMPILE;
	    p->subp->flags &= ~P_EXEC;
	    parser_free_aux_nodes(p->subp);
	    gen_cleanup(p->subp);
	    free(p->subp);
	    p->subp = NULL;
	}  
    }
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
	p->err = E_PARSE;
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

    if (p->ret != NULL && p->ret->t != VEC) {
	fprintf(stderr, " ret type != VEC (=%d), p->err = %d\n", p->ret->t, p->err);
    } else if (p->ret == NULL) {
	fprintf(stderr, " ret = NULL, p->err = %d\n", p->err);
    }

    fprintf(stderr, " input = '%s'\n", p->input);
    
    if (!p->err) {
	p->err = E_DATA;
    }
}

int realgen (const char *s, parser *p, DATASET *dset, PRN *prn, 
	     int flags)
{
    int t;

#if LHDEBUG || EDEBUG
    fprintf(stderr, "*** realgen: task = %s\n", (flags & P_COMPILE)?
	    "compile" : (flags & P_EXEC)? "exec" : "normal");
#endif

    if (flags & P_EXEC) {
#if EDEBUG
	fprintf(stderr, "*** printing p->tree (before reinit)\n");
	print_tree(p->tree, p, 0);
#endif
	parser_reinit(p, dset, prn);
	if (p->err) {
	    fprintf(stderr, "error in parser_reinit\n");
	    return p->err;
	} else if (p->op == INC || p->op == DEC || (p->flags & P_PRINT)) {
	    p->ret = lhs_copy_node(p);
	    return p->err;
	} else {
	    goto starteval;
	}
    } else {
	parser_init(p, s, dset, prn, flags);
	if (p->err) {
	    if (gretl_function_depth() == 0) {
		errmsg(p->err, prn);
	    }
	    return p->err;
	}
    }

#if EDEBUG
    fprintf(stderr, "after parser (re-)init, p->err = %d (decl? %s)\n", 
	    p->err, (p->flags & P_DECL)? "yes" : "no");
#endif

    if (p->flags & P_DECL) {
	decl_check(p, flags);
	return p->err;
    }

    if (p->op == INC || p->op == DEC || (p->flags & P_PRINT)) {
	if (!(p->flags & P_COMPILE)) {
	    p->ret = lhs_copy_node(p);
	}
	return p->err;
    }

    lex(p);
    if (p->err) {
	fprintf(stderr, "realgen: exiting on lex() error %d\n", p->err);
	return p->err;
    }

    p->tree = expr(p);
    if (p->err) {
	fprintf(stderr, "realgen: exiting on expr() error %d\n", p->err);
	return p->err;
    }

#if EDEBUG
    fprintf(stderr, "realgen: p->tree at %p, type %d\n", (void *) p->tree, 
	    p->tree->t);
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
	context_error(c, p);
	return p->err;
    }    

    if (flags & P_COMPILE) {
	return p->err;
    }

    /* set "simple sort" or other flags here if relevant */
    if (!p->err) {
	maybe_set_return_flags(p);
    }

 starteval:

    parser_aux_init(p);

    if (autoreg(p)) {
	/* e.g. y = b*y(-1) : evaluate dynamically */
	for (t=p->dset->t1; t<p->dset->t2 && !p->err; t++) {
	    const double *x;

	    p->aux_i = 0;
	    p->obs = t;
#if EDEBUG
	    fprintf(stderr, "\n*** autoreg: p->obs = %d\n", p->obs);
#endif
	    p->ret = eval(p->tree, p);
	    if (p->ret != NULL && p->ret->t == VEC) {
		x = p->ret->v.xvec;
		if (!na(x[t])) { 
#if EDEBUG
		    fprintf(stderr, "writing xvec[%d] = %g into Z[%d][%d]\n",
			    t, x[t], p->lh.v, t);
#endif
		    p->dset->Z[p->lh.v][t] = x[t];
		} 
	    } else {
		autoreg_error(p, t);
	    }
	    if (t == p->dset->t1) {
		p->flags &= ~P_START;
	    } 
	}
	p->obs = t;
    }

    p->aux_i = 0;

#if EDEBUG
    if (flags & P_EXEC) {
	fprintf(stderr, "*** printing p->tree (about to start eval)\n");
	print_tree(p->tree, p, 0);
    }
#endif

    if (!p->err) {
	p->ret = eval(p->tree, p);
    }

#if EDEBUG
    fprintf(stderr, "realgen: post-eval, err = %d\n", p->err);
#endif

#if EDEBUG > 1
    printnode(p->ret, p);
    pputc(prn, '\n');
#endif

#if EDEBUG
    fprintf(stderr, "calling parser_free_aux_nodes\n");
#endif
    parser_free_aux_nodes(p);

#if 1
    gen_check_errvals(p);
#endif

    return p->err;
}
