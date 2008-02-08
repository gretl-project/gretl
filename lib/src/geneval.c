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
#include "gretl_fft.h"

#include <errno.h>

#if GENDEBUG
# define EDEBUG GENDEBUG
#else
# define EDEBUG 0
#endif

#define is_aux_node(n) (n != NULL && (n->flags & AUX_NODE))
#define is_tmp_node(n) (n != NULL && (n->flags & TMP_NODE))

#define nullmat_ok(f) (f == ROWS || f == COLS || f == DET || \
		       f == LDET || f == DIAG || f == TRANSP || \
		       f == TVEC || f == VECH || f == UNVECH)

#define ok_lcat_node(n) (n->t == LIST || n->t == LVEC || n->t == NUM || \
			 (n->t == VEC && (n->flags & UVEC_NODE)))

#define getvec(n,p) ((n->flags & UVEC_NODE)? \
		     (*p->Z)[n->v.idnum] : n->v.xvec)

static void parser_init (parser *p, const char *str, 
			 double ***pZ, DATAINFO *dinfo,
			 PRN *prn, int flags);

static void printnode (const NODE *t, const parser *p);

static NODE *eval (NODE *t, parser *p);

static void node_type_error (int ntype, int goodt, NODE *bad, parser *p);

static const char *typestr (int t)
{
    switch (t) {
    case NUM:
	return "scalar";
    case VEC:
	return "series";
    case MAT:
    case UMAT:
	return "matrix";
    case STR:
	return "string";
    case U_ADDR:
	return "address";
    default:
	return "?";
    }
}

static void free_tree (NODE *t, const char *msg)
{
    if (t == NULL) {
	return;
    }

#if EDEBUG
    fprintf(stderr, "%-8s: starting with t at %p (type %03d)\n", msg, 
	    (void *) t, t->t);
#endif

    /* free recursively */
    if (bnsym(t->t)) {
	int i;

	for (i=0; i<t->v.bn.n_nodes; i++) {
	    free_tree(t->v.bn.n[i], msg);
	}
	free(t->v.bn.n);
    } else if (b3sym(t->t)) {
	free_tree(t->v.b3.l, msg);
	free_tree(t->v.b3.m, msg);
	free_tree(t->v.b3.r, msg);
    } else if (b2sym(t->t)) {
	free_tree(t->v.b2.l, msg);
	free_tree(t->v.b2.r, msg);
    } else if (b1sym(t->t)) {
	free_tree(t->v.b1.b, msg);
    } 

#if EDEBUG
    fprintf(stderr, "%-8s: freeing node at %p (type %03d, flags = %d)\n", msg, 
	    (void *) t, t->t, t->flags);
#endif

    if (is_tmp_node(t)) {
	if (t->t == VEC) {
	    free(t->v.xvec);
	} else if (t->t == IVEC || t->t == LVEC) {
	    free(t->v.ivec);
	} else if (t->t == MAT) {
	    gretl_matrix_free(t->v.m);
	} else if (t->t == MSPEC) {
	    free(t->v.mspec);
	}
    }

    if (freestr(t->t)) {
	free(t->v.str);
    }

    free(t);
}

static void parser_aux_init (parser *p)
{
    if (p->ecount == 0) {
	p->aux = NULL;
	p->n_aux = 0;
    }
    p->aux_i = 0;
}

void parser_free_aux_nodes (parser *p)
{
    int i;

    if (p->aux != NULL) {
	for (i=0; i<p->n_aux; i++) {
	    if (p->aux[i] != p->ret) {
		free_tree(p->aux[i], "Aux");
	    }
	}
	free(p->aux);
    }
}

static NODE *newmdef (int k)
{  
    NODE *n = malloc(sizeof *n);
    int i;

#if MDEBUG
    fprintf(stderr, "newmdef: allocated node at %p\n", (void *) n);
#endif

    if (n == NULL) {
	return NULL;
    }

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
	n->t = MDEF;
	n->v.bn.n_nodes = k;
	n->flags = 0;
    }

    return n;
}

/* new node to hold array of doubles */

static NODE *newvec (int n, int tmp)
{  
    NODE *b = malloc(sizeof *b);
    int i;

#if EDEBUG
    fprintf(stderr, "newvec: allocated node at %p\n", (void *) b);
#endif

    if (b != NULL) {
	b->t = VEC;
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
    NODE *b = malloc(sizeof *b);

#if EDEBUG
    fprintf(stderr, "newivec: allocated node at %p\n", (void *) b);
#endif

    if (b != NULL) {
	b->t = type;
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
    NODE *b = malloc(sizeof *b);

#if EDEBUG
    fprintf(stderr, "newmat: allocated node at %p\n", (void *) b);
#endif

    if (b != NULL) {
	b->t = MAT;
	b->flags = (tmp)? TMP_NODE : 0;
	b->v.m = NULL;
    }

    return b;
}

/* new node to hold a matrix specification */

static NODE *newmspec (void)
{
    NODE *b = malloc(sizeof *b);

#if EDEBUG
    fprintf(stderr, "newmspec: allocated node at %p\n", (void *) b);
#endif

    if (b != NULL) {
	b->t = MSPEC;
	b->flags = TMP_NODE;
	b->v.mspec = NULL;
    }

    return b;
}

/* new node to hold a list reference */

static NODE *newlist (void)
{  
    NODE *b = malloc(sizeof *b);

#if EDEBUG
    fprintf(stderr, "newlist: allocated node at %p\n", (void *) b); 
#endif

    if (b != NULL) {
	b->t = LIST;
	b->flags = TMP_NODE;
	b->v.str = NULL;
    }    

    return b;
}

static int node_allocate_matrix (NODE *t, int m, int n, parser *p)
{
    t->v.m = gretl_matrix_alloc(m, n);
    if (t->v.m == NULL) {
	p->err = E_ALLOC;
    }

    return p->err;
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

#define repeat_exec(p) (p->ecount > 0)

/* get an auxiliary node: if starting from scratch we allocate
   a new node, otherwise we look up an existing one */

static NODE *get_aux_node (parser *p, int t, int n, int tmp)
{
    NODE *ret = NULL;

    if (starting(p) && !repeat_exec(p)) {
	if (t == NUM) {
	    ret = newdbl(NADBL);
	} else if (t == VEC) {
	    ret = newvec(n, tmp);
	} else if (t == IVEC) {
	    ret = newivec(n, IVEC);
	} else if (t == LVEC) {
	    ret = newivec(n, LVEC);
	} else if (t == MAT) {
	    ret = newmat(tmp);
	} else if (t == MSPEC) {
	    ret = newmspec();
	} else if (t == MDEF) {
	    ret = newmdef(n);
	} else if (t == LIST) {
	    ret = newlist();
	}

	if (ret == NULL) {
	    p->err = (t == 0)? E_DATA : E_ALLOC;
	} else if (add_aux_node(p, ret)) {
	    free_tree(ret, "On error");
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

static NODE *aux_vec_node (parser *p, int n)
{
    if (p->dinfo->n == 0) {
	gretl_errmsg_set(_("No dataset is in place"));
	return NULL;
    } else {
	return get_aux_node(p, VEC, n, 1);
    }
}

static NODE *aux_ivec_node (parser *p, int n)
{
    if (p->dinfo->n == 0) {
	gretl_errmsg_set(_("No dataset is in place"));
	return NULL;
    } else {
	return get_aux_node(p, IVEC, n, 1);
    }
}

static NODE *aux_lvec_node (parser *p, int n)
{
    if (p->dinfo->n == 0) {
	gretl_errmsg_set(_("No dataset is in place"));
	return NULL;
    } else {
	return get_aux_node(p, LVEC, n, 1);
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

static NODE *aux_mdef_node (parser *p, int n)
{
    return get_aux_node(p, MDEF, n, 0);
}

static NODE *aux_list_node (parser *p)
{
    return get_aux_node(p, LIST, 0, 0);
}

static NODE *aux_any_node (parser *p)
{
    return get_aux_node(p, 0, 0, 0);
}

static void eval_warning (parser *p, int op)
{
    if (*p->warning == '\0') {
	if (op == B_POW) {
	    strcpy(p->warning, _("invalid operands for '^'"));
	} else if (op == LOG) {
	    strcpy(p->warning, _("invalid argument for log()"));
	} else if (op == SQRT) {
	    strcpy(p->warning, _("invalid argument for sqrt()"));
	} else if (op == EXP) {
	    strcpy(p->warning, _("invalid argument for exp()"));
	}
    }
}

/* implementation of binary operators for scalar operands
   (also increment/decrement operators) */

static double xy_calc (double x, double y, int op, parser *p)
{
    double z = NADBL;

#if EDEBUG > 1
    fprintf(stderr, "xy_calc: x = %g, y = %g, op = '%s'\n",
	    x, y, getsymb(op, NULL));
#endif

    /* assignment */
    if (op == B_ASN) {
	return y;
    }    

    /* testing for presence of NAs? */
    if ((p->flags & P_NATEST) && (na(x) || na(y))) {
	return NADBL;
    }

    /* special case: 0 * anything (including even NA) = 0 */
    if (op == B_MUL && (x == 0.0 || y == 0.0)) {
	return 0.0;
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
	return (int) x % (int) y;
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
	    eval_warning(p, op);
	}
	return z;
    default: 
	return z;
    }
}

static int dist_argc (char *s, int f)
{
    if (strlen(s) > 1) {
	return 0;
    }

    switch (s[0]) {
    case 'u':
    case 'U':
	s[0] = 'u';
	return (f == RANDGEN)? 2 : 0;
    case '1':
    case 'z':
    case 'n':
    case 'N':
	s[0] = 'z';
	return (f == RANDGEN)? 2 : 1;
    case '2':
    case 't':
	s[0] = 't';
	return (f == RANDGEN)? 1 : 2;
    case '3':
    case 'c':
    case 'x':
    case 'X':
	s[0] = 'X';
	return (f == RANDGEN)? 1 : 2;
    case '4':
    case 'f':
    case 'F':
	s[0] = 'F';
	return (f == RANDGEN)? 2 : 3;
    case '5':
    case 'g':
    case 'G':
	s[0] = 'G';
	return (f == CRIT || f == INVCDF)? 0 : 
	    (f == RANDGEN)? 2 : 3;
    case '6':
    case 'b':
    case 'B':
	s[0] = 'B';
	return (f == RANDGEN)? 2 : 3;
    case '7':
    case 'D':
	s[0] = 'D';
	return (f == CDF)? 3 : 0;
    case '8':
    case 'p':
    case 'P':
	s[0] = 'P';
	return (f == RANDGEN)? 1 : 2;
    }

    return 0;
}

/* make a column vector containing the 1-based observation numbers
   corresponding to the non-zero entries in the series under node n
*/

static NODE *make_series_mask (NODE *n, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = getvec(n, p);
	gretl_matrix *v;
	int t, s, T = 0;

	for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	    if (x[t] != 0) {
		T++;
	    }
	}

	if (T == 0) {
	    p->err = E_DATA;
	    return NULL;
	}

	v = gretl_column_vector_alloc(T);
	if (v == NULL) {
	    p->err = E_ALLOC;
	    return NULL;
	}

	s = 0;
	for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	    if (x[t] != 0) {
		gretl_vector_set(v, s++, t + 1);
	    }
	}
	
	ret->v.m = v;
    }

    return ret;
}

static double scalar_pdist (int t, char d, double *parm,
			    int np, parser *p)
{
    double x = NADBL;
    int i;

    for (i=0; i<np; i++) {
	if (na(parm[i])) {
	    return NADBL;
	}
    }

    if (t == PVAL) {
	x = gretl_get_pvalue(d, parm);
    } else if (t == CDF) {
	x = gretl_get_cdf(d, parm);
    } else if (t == INVCDF) {
	x = gretl_get_cdf_inverse(d, parm);
    } else if (t == CRIT) {
	x = gretl_get_critval(d, parm);
    } else {
	p->err = 1;
    }

    return x;
}

static double *series_pdist (int t, char d, double *parm,
			     double *pvec, double *bvec,
			     int np, parser *p)
{
    double *xvec;
    int s;

    xvec = malloc(p->dinfo->n * sizeof *xvec);
    if (xvec == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    for (s=0; s<p->dinfo->n; s++) {
	if (s < p->dinfo->t1 || s > p->dinfo->t2 || pvec[s] == NADBL) {
	    xvec[s] = NADBL;
	} else {
	    if (pvec != NULL) {
		parm[np-1] = pvec[s];
	    }
	    if (bvec != NULL) {
		parm[np-2] = bvec[s];
	    }
	    xvec[s] = scalar_pdist(t, d, parm, np, p);
	}
    }

    return xvec;
}

static gretl_matrix *matrix_pdist (int t, char d, double *parm,
				   gretl_matrix *pmat, 
				   gretl_matrix *bmat,
				   int np, parser *p)
{
    gretl_matrix *m;
    gretl_matrix *a;
    double x;
    int i, n;

    if (pmat != NULL && bmat != NULL) {
	if (pmat->rows != bmat->rows ||
	    pmat->cols != bmat->cols) {
	    p->err = E_NONCONF;
	    return NULL;
	}
    }

    a = (pmat != NULL)? pmat : bmat;

    m = gretl_matrix_alloc(a->rows, a->cols);
    if (m == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    n = a->rows * a->cols;

    for (i=0; i<n; i++) {
	if (pmat != NULL) {
	    parm[np-1] = pmat->val[i];
	}
	if (bmat != NULL) {
	    parm[np-2] = bmat->val[i];
	}
	x = scalar_pdist(t, d, parm, np, p);
	if (na(x)) {
	    p->err = E_MISSDATA;
	    gretl_matrix_free(m);
	    m = NULL;
	    break;
	} else {
	    m->val[i] = x;
	}
    }

    return m;
}

/* return a node containing the evaluated result of a
   probability distriution function */

static NODE *eval_pdist (NODE *n, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *e, *s, *r = n->v.b1.b;
	int i, k, argc, m = r->v.bn.n_nodes;
	int rgen = (n->t == RANDGEN);
	double *pvec = NULL, *bvec = NULL;
	gretl_matrix *pmat = NULL, *bmat = NULL;
	double parm[3];
	char d, *dist, code[2];

	if (m < 2 || m > 4) {
	    p->err = E_INVARG;
	    goto disterr;
	}

	s = r->v.bn.n[0];
	if (s->t == STR) {
	    dist = s->v.str;
	} else if (s->t == NUM && s->v.xval > 0 && s->v.xval < 10) {
	    sprintf(code, "%d", (int) s->v.xval);
	    dist = code;
	} else {
	    p->err = E_INVARG;
	    goto disterr;
	}

	argc = dist_argc(dist, n->t);

	if (argc != m - 1) {
	    p->err = E_INVARG;
	    goto disterr;
	}

	d = dist[0];
	k = argc - 1;

	for (i=0; i<argc && !p->err; i++) {
	    s = r->v.bn.n[i+1];
	    if (s->t == NUM) {
		parm[i] = s->v.xval;
	    } else if (i == k && !rgen && s->t == VEC && bmat == NULL) {
		pvec = getvec(s, p);
	    } else if (i == k && !rgen && s->t == MAT && bvec == NULL) {
		pmat = s->v.m;
	    } else if (i == k-1 && d == 'D' && s->t == VEC) {
		bvec = getvec(s, p);
	    } else if (i == k-1 && d == 'D' && s->t == MAT) {
		bmat = s->v.m;
	    } else {
		e = eval(s, p);
		if (p->err) {
		    goto disterr;
		}
		if (e->t == NUM) {
		    parm[i] = e->v.xval;
		    free_tree(s, "Pdist");
		    r->v.bn.n[i+1] = NULL;
		} else if (i == k && !rgen && e->t == VEC && bmat == NULL) {
		    pvec = getvec(e, p);
		    free_tree(s, "Pdist");
		    r->v.bn.n[i+1] = NULL;
		} else if (i == k && !rgen && e->t == MAT && bvec == NULL) {
		    pmat = e->v.m;
		    free_tree(s, "Pdist");
		    r->v.bn.n[i+1] = NULL;
		} else if (i == k-1 && d == 'D' && e->t == VEC) {
		    bvec = getvec(e, p);
		    free_tree(s, "Pdist");
		    r->v.bn.n[i+1] = NULL;
		} else if (i == k-1 && d == 'D' && e->t == MAT) {
		    bmat = e->v.m;
		    free_tree(s, "Pdist");
		    r->v.bn.n[i+1] = NULL;
		} else {
		    p->err = E_INVARG;
		    goto disterr;
		}
	    }
	}

	if (rgen || pvec != NULL || bvec != NULL) {
	    ret = aux_vec_node(p, 0);
	} else if (pmat != NULL || bmat != NULL) {
	    ret = aux_matrix_node(p);
	} else {
	    ret = aux_scalar_node(p);
	}

	if (ret == NULL) {
	    goto disterr;
	}

	if (rgen) {
	    ret->v.xvec = gretl_get_random_series(d, parm, p->dinfo,
						  &p->err);
	} else if (pvec != NULL || bvec != NULL) {
	    ret->v.xvec = series_pdist(n->t, d, parm, pvec, bvec, argc, p);
	} else if (pmat != NULL || bmat != NULL) {
	    ret->v.m = matrix_pdist(n->t, d, parm, pmat, bmat, argc, p);
	} else {
	    ret->v.xval = scalar_pdist(n->t, d, parm, argc, p);
	}

    } else {
	ret = aux_any_node(p);
    }

  disterr:  

    return ret;
}

/* look up and return numerical values of symbolic constants */

static NODE *retrieve_const (NODE *n, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	switch (n->v.idnum) {
	case CONST_PI:
	    ret->v.xval = M_PI;
	    break;
	case CONST_NA:
	    ret->v.xval = NADBL;
	    break;
	case CONST_WIN32:
#ifdef WIN32
	    ret->v.xval = 1;
#else
	    ret->v.xval = 0;
#endif
	    break;
	}
    }

    return ret;
}

static NODE *scalar_calc (NODE *x, NODE *y, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.xval = xy_calc(x->v.xval, y->v.xval, f, p);
    }

    return ret;
}

/* try interpreting the string on the right as identifying
   an observation number */

static NODE *number_string_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;
    double yt, xt = (l->t == NUM)? l->v.xval : 0.0;
    double *x = NULL;
    int t, t1, t2;

    t = dateton(r->v.str, p->dinfo);
    if (t >= 0) {
	yt = t + 1;
    } else {
	p->err = 1;
	return NULL;
    }

    ret = aux_vec_node(p, p->dinfo->n);
    if (ret == NULL) {
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "number_string_calc: l=%p, r=%p, ret=%p\n", 
	    (void *) l, (void *) r, (void *) ret);
#endif

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

    if (l->t == VEC) x = getvec(l, p);

    for (t=t1; t<=t2; t++) {
	if (x != NULL) {
	    xt = x[t];
	}
	ret->v.xvec[t] = xy_calc(xt, yt, f, p);
    }

    return ret;
}

static NODE *series_calc (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret;
    double xt = (l->t == NUM)? l->v.xval : 0.0;
    double yt = (r->t == NUM)? r->v.xval : 0.0;
    const double *x = NULL, *y = NULL;
    int t, t1, t2;

    ret = aux_vec_node(p, p->dinfo->n);
    if (ret == NULL) {
	return NULL;
    }

#if EDEBUG
    fprintf(stderr, "series_calc: l=%p, r=%p, ret=%p\n", 
	    (void *) l, (void *) r, (void *) ret);
#endif

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

    if (l->t == VEC) x = getvec(l, p);
    if (r->t == VEC) y = getvec(r, p);

    for (t=t1; t<=t2; t++) {
	if (x != NULL) {
	    xt = x[t];
	}
	if (y != NULL) {
	    yt = y[t];
	}
	ret->v.xvec[t] = xy_calc(xt, yt, f, p);
    }

    return ret;
}

static int op_symbol (int op)
{
    switch (op) {
    case DOTMULT: return '*';
    case DOTDIV:  return '/';
    case DOTPOW:  return '^';
    case DOTADD:  return '+';
    case DOTSUB:  return '-';
    case DOTEQ:   return '=';
    case DOTGT:   return '>';
    case DOTLT:   return '<';
    default: return 0;
    }
}

/* return allocated result of binary operation performed on
   two matrices */

static gretl_matrix *real_matrix_calc (const gretl_matrix *A, 
				       const gretl_matrix *B, 
				       int op, int *err) 
{
    gretl_matrix *C = NULL;
    gretl_matrix *D = NULL;
    int ra, ca;
    int rb, cb;
    int r, c;

    if (gretl_is_null_matrix(A) ||
	gretl_is_null_matrix(B)) {
	if (op != MCCAT && op != MRCAT) {
	    *err = E_NONCONF;
	    return NULL;
	}
    }

    switch (op) {
    case B_ADD:
    case B_SUB:
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);

	if (ra == 1 && ca == 1) {
	    C = gretl_matrix_copy(B);
	} else {
	    C = gretl_matrix_copy(A);
	}
	if (C == NULL) {
	    *err = E_ALLOC;
	} else if (ra == 1 && ca == 1) {
	    if (op == B_ADD) {
		*err = gretl_matrix_add_to(C, A);
	    } else {
		gretl_matrix_switch_sign(C);
		*err = gretl_matrix_add_to(C, A);
	    }
	} else {
	    if (op == B_ADD) {
		*err = gretl_matrix_add_to(C, B);
	    } else {
		*err = gretl_matrix_subtract_from(C, B);
	    }
	}
	break;
    case MCCAT:
	C = gretl_matrix_col_concat(A, B, err);
	break;
    case MRCAT:
	C = gretl_matrix_row_concat(A, B, err);
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
    case QFORM:
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
	/* matrix "division" */
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	C = gretl_matrix_copy(A);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else {
	    if (rb == 1 && cb == 1) {
		*err = gretl_matrix_divide_by_scalar(C, B->val[0]);
	    } else {
		D = gretl_matrix_copy(B);
		if (D == NULL) {
		    gretl_matrix_free(C);
		    C = NULL;
		    *err = E_ALLOC;
		} else {	
		    *err = gretl_LU_solve(D, C);
		    gretl_matrix_free(D);
		}
	    }
	}
	break;
    case DOTMULT:
    case DOTDIV:
    case DOTPOW:
    case DOTADD:
    case DOTSUB:
    case DOTEQ:
    case DOTGT:
    case DOTLT:
	/* apply operator element-wise */
	C = gretl_matrix_dot_op(A, B, op_symbol(op), err);
	break;
    case KRON:
	/* Kronecker product */
	C = gretl_matrix_kronecker_product_new(A, B, err);
	break;
    case CMULT:
	C = gretl_matrix_complex_multiply(A, B, err);
	break;
    case CDIV:
	C = gretl_matrix_complex_divide(A, B, err);
	break;
    case MRSEL:
	C = gretl_matrix_bool_sel(A, B, 1, err);
	break;
    case MCSEL:
	C = gretl_matrix_bool_sel(A, B, 0, err);
	break;
    default:
	*err = E_TYPES;
	break;
    } 

    if (*err && C != NULL) {
	gretl_matrix_free(C);
	C = NULL;
    }

    return C;
}

static gretl_matrix *
tmp_matrix_from_series (const double *x, const DATAINFO *pdinfo,
			int *err)
{
    gretl_matrix *m = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int i, t;

    /* FIXME autoregressive case? */

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (xna(x[t])) {
	    *err = E_MISSDATA;
	    return NULL;
	}
    }

    m = gretl_column_vector_alloc(T);
    if (m == NULL) {
	*err = E_ALLOC;
    } else {
	i = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    m->val[i++] = x[t];
	}
    }

    return m;
}

/* one of the operands is a matrix, the other a series,
   which gets "promoted" to a (temporary) matrix if possible
*/

static NODE *matrix_series_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *a = NULL;
	gretl_matrix *b = NULL;
	gretl_matrix *c = NULL;
	const double *xvec;

	if (l->t == VEC) {
	    xvec = getvec(l, p);
	    a = tmp_matrix_from_series(xvec, p->dinfo, &p->err);
	    c = a;
	    b = r->v.m;
	} else {
	    xvec = getvec(r, p);
	    a = l->v.m;
	    b = tmp_matrix_from_series(xvec, p->dinfo, &p->err);
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

#define comp_op(o) (o == B_EQ  || o == B_NEQ || \
                    o == B_LT  || o == B_GT || \
                    o == B_LTE || o == B_GTE)

/* one of the operands is a matrix, the other a scalar, giving a
   matrix result unless we're looking at a comparison operator.
*/

static NODE *matrix_scalar_calc (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = NULL;

    if (op == KRON) {
	p->err = E_TYPES;
	return NULL;
    }

    if (starting(p)) {
	const gretl_matrix *m = NULL;
	int comp = comp_op(op);
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
	if (op == B_POW && matrix_pow_check(l->t, x, m, p)) {
	    return NULL;
	}

	/* mod: scalar must be on the right */
	if (op == B_MOD && l->t == NUM) {
	    p->err = E_TYPES;
	    return NULL;
	}

	if (comp) {
	    ret = aux_scalar_node(p);
	} else {
	    ret = aux_matrix_node(p);
	}

	if (ret == NULL) { 
	    return NULL;
	}

	if (op == B_POW) {
	    ret->v.m = gretl_matrix_pow(m, (int) x, &p->err);
	    return ret;
	}

	if (comp) {
	    ret->v.xval = 1;
	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    if (xy_calc(x, m->val[i], op, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}
	    } else {
		for (i=0; i<n; i++) {
		    if (xy_calc(m->val[i], x, op, p) == 0) {
			ret->v.xval = 0;
			break;
		    }
		}		
	    }
	} else {
	    if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
		free_tree(ret, "On error");
		return NULL;
	    }

	    if (l->t == NUM) {
		for (i=0; i<n; i++) {
		    y = xy_calc(x, m->val[i], op, p);
		    ret->v.m->val[i] = y;
		}
	    } else {
		for (i=0; i<n; i++) {
		    y = xy_calc(m->val[i], x, op, p);
		    ret->v.m->val[i] = y;
		}	
	    }
	} 
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *numeric_jacobian (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = r->v.str;

	ret = aux_matrix_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.m = fdjac(l->v.m, s, p->Z, p->dinfo, &p->err);
    } else {
	ret = aux_matrix_node(p);
    }

    return ret;
}

static NODE *BFGS_maximize (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *m = l->v.m;
	const char *s = r->v.str;

	if (gretl_is_null_matrix(m)) {
	    p->err = E_DATA;
	    return NULL;
	}

	ret = aux_scalar_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = user_BFGS(m, s, p->Z, p->dinfo, 
				p->prn, &p->err);
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *matrix_csv_write (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	const char *s = r->v.str;

	ret = aux_scalar_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.xval = gretl_matrix_write_as_text(l->v.m, s);
    } else {
	ret = aux_scalar_node(p);
    }

    return ret;
}

static NODE *matrix_lag (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *m = l->v.m;
	int k = r->v.xval;

	if (gretl_is_null_matrix(m)) {
	    p->err = E_DATA;
	    return NULL;
	}	

	ret = aux_matrix_node(p);
	if (ret == NULL) { 
	    return NULL;
	}

	ret->v.m = gretl_matrix_lag(m, k, 0.0);
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

/* both operands are matrices */

static NODE *matrix_bool (NODE *l, NODE *r, int op, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (op == B_OR || op == B_AND) {
	/* not meaningful? */
	p->err = E_TYPES;
	return NULL;
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

    if (*gretl_errmsg != '\0') {
	errmsg(p->err, p->prn);
    }
}

/* functions taking a matrix argument and returning a
   scalar result */

static NODE *matrix_to_scalar_func (NODE *n, int f, parser *p)
{
    gretl_matrix *m = n->v.m;
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {

	switch (f) {
	case ROWS:
	    ret->v.xval = m->rows;
	    break;
	case COLS:
	    ret->v.xval = m->cols;
	    break;
	case DET:
	case LDET:
	    ret->v.xval = user_matrix_get_determinant(m, f, &p->err);
	    break;
	case TRACE:
	    ret->v.xval = gretl_matrix_trace(m, &p->err);
	    break;
	case NORM1:
	    ret->v.xval = gretl_matrix_one_norm(m);
	    break;
	case INFNORM:
	    ret->v.xval = gretl_matrix_infinity_norm(m);
	    break;
	case RCOND:
	    ret->v.xval = gretl_symmetric_matrix_rcond(m, &p->err);
	    break;
	case RANK:
	    ret->v.xval = gretl_matrix_rank(m, &p->err);
	    break;
	default:
	    p->err = 1;
	    break;
	}

	if (p->err) {
	    matrix_error(p);
	}    
    }

    return ret;
}

static NODE *matrix_princomp (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = l->v.m;
	int k = r->v.xval;

	ret->v.m = gretl_matrix_pca(m, k, &p->err);
    }

    return ret;
}

static void matrix_minmax_indices (int f, int *mm, int *rc, int *idx)
{
    *mm = (f == MAXR || f == MAXC || f == IMAXR || f == IMAXC);
    *rc = (f == MINC || f == MAXC || f == IMINC || f == IMAXC);
    *idx = (f == IMINR || f == IMINC || f == IMAXR || f == IMAXC);
}

static NODE *matrix_to_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *m = n->v.m;
	int a, b, c;

	if (gretl_is_null_matrix(m) && !nullmat_ok(f)) {
	    p->err = E_DATA;
	    goto finalize;
	}

	gretl_error_clear();

	switch (f) {
	case SUMC:
	    ret->v.m = gretl_matrix_column_sum(m, &p->err);
	    break;
	case SUMR:
	    ret->v.m = gretl_matrix_row_sum(m, &p->err);
	    break;
	case MEANC:
	    ret->v.m = gretl_matrix_column_mean(m, &p->err);
	    break;
	case MEANR:
	    ret->v.m = gretl_matrix_row_mean(m, &p->err);
	    break;
	case SD:
	    ret->v.m = gretl_matrix_column_sd(m, &p->err);
	    break;
	case MCOV:
	    ret->v.m = gretl_covariance_matrix(m, 0, &p->err);
	    break;
	case MCORR:
	    ret->v.m = gretl_covariance_matrix(m, 1, &p->err);
	    break;
	case CUM:
	    ret->v.m = gretl_matrix_cumcol(m, &p->err);
	    break;
	case DIF:
	    ret->v.m = gretl_matrix_diffcol(m, 0, &p->err);
	    break;
	case RESAMPLE:
	    ret->v.m = gretl_matrix_resample(m, 0, &p->err);
	    break;
	case CDEMEAN:
	case CHOL:
	case INV:
	case INVPD:
	case GINV:
	case UPPER:
	case LOWER:
	    ret->v.m = user_matrix_matrix_func(m, f, &p->err);
	    break;
	case DIAG:
	    ret->v.m = gretl_matrix_get_diagonal(m, &p->err);
	    break;
	case TRANSP:
	    ret->v.m = gretl_matrix_copy_transpose(m);
	    break;
	case TVEC:
	    ret->v.m = user_matrix_vec(m, &p->err);
	    break;
	case VECH:
	    ret->v.m = user_matrix_vech(m, &p->err);
	    break;
	case UNVECH:
	    ret->v.m = user_matrix_unvech(m, &p->err);
	    break;
	case NULLSPC:
	    ret->v.m = gretl_matrix_right_nullspace(m, &p->err);
	    break;
	case MEXP:
	    ret->v.m = gretl_matrix_exp(m, &p->err);
	    break;
	case FFT:
	    ret->v.m = gretl_matrix_fft(m, &p->err);
	    break;
	case FFTI:
	    ret->v.m = gretl_matrix_ffti(m, &p->err);
	    break;
	case POLROOTS:
	    ret->v.m = gretl_matrix_polroots(m, &p->err);
	    break;
	case MINC:
	case MAXC:
	case MINR:
	case MAXR:
	case IMINC:
	case IMAXC:
	case IMINR:
	case IMAXR:  
	    matrix_minmax_indices(f, &a, &b, &c);
	    ret->v.m = gretl_matrix_minmax(m, a, b, c, &p->err);
	    break;
	default:
	    break;
	}

	if (ret->v.m == n->v.m) {
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

static NODE *string_to_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	gretl_error_clear();

	switch (f) {
	case MREAD:
	    ret->v.m = gretl_matrix_read_from_text(n->v.str, &p->err);
	    break;
	default:
	    break;
	}

	if (ret->v.m == NULL) {
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
	    r = r->v.b1.b;
	    if (r->t == UMAT) {
		rname = r->v.str;
	    } else {
		p->err = 1;
		strcpy(gretl_errmsg, "Expected the address of a matrix");
		return ret;
	    }
	}

	switch (f) {
	case QR:
	    ret->v.m = user_matrix_QR_decomp(m, rname, &p->err);
	    break;
	case EIGSYM:
	    ret->v.m = user_matrix_eigen_analysis(m, rname, 1, &p->err);
	    break;
	case EIGGEN:
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

static int ok_matrix_dim (double xr, double xc, int f)
{
    double xm, imax = (double) INT_MAX;

    if (f == SEQ) {
	/* negative parameters are OK */
	return (fabs(xr) < imax && 
		fabs(xc) < imax);
    }

    xm = xr * xc;

    if (f == IMAT || f == ZEROS || f == ONES || f == MUNIF || \
	f == MNORM) {
	/* zero is OK for matrix creation functions, which then 
	   return an empty matrix 
	*/
	return (xr >= 0 && xr <= imax && 
		xc >= 0 && xc <= imax &&
		xm >= 0 && xm <= imax);
    } else {
	return (xr > 0 && xr <= imax && 
		xc > 0 && xc <= imax &&
		xm > 0 && xm <= imax);
    }
}

static NODE *matrix_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	double xr = l->v.xval;
	double xc = (f == IMAT)? l->v.xval : r->v.xval;
	int rows, cols;

	gretl_error_clear();

	if (!ok_matrix_dim(xr, xc, f)) {
	    p->err = E_DATA;
	    matrix_error(p);
	    return ret;
	}

	rows = xr;
	cols = xc;

	switch (f) {
	case IMAT:
	    ret->v.m = gretl_identity_matrix_new(rows);
	    break;
	case ZEROS:
	    ret->v.m = gretl_zero_matrix_new(rows, cols);
	    break;
	case ONES:
	    ret->v.m = gretl_unit_matrix_new(rows, cols);
	    break;
	case SEQ:
	    ret->v.m = gretl_matrix_seq(rows, cols);
	    break;
	case MUNIF:
	    ret->v.m = gretl_random_matrix_new(rows, cols, 
					       D_UNIFORM);
	    break;
	case MNORM:
	    ret->v.m = gretl_random_matrix_new(rows, cols,
					       D_NORMAL);
	    break;
	default:
	    break;
	}
    }

    return ret;
}

/* compose a sub-matrix specification, from scalars and/or
   index matrices */

static matrix_subspec *build_mspec (NODE *l, NODE *r, int *err)
{
    matrix_subspec *mspec;

    mspec = malloc(sizeof *mspec);
    if (mspec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    if (l->t == DUM) {
	if (l->v.idnum == DUM_DIAG) {
	    mspec->type[0] = SEL_DIAG;
	    mspec->type[1] = SEL_ALL;
	    return mspec;
	} else {
	    *err = E_TYPES;
	    goto bailout;
	}
    }

    if (l->t == NUM) {
	mspec->type[0] = SEL_RANGE;
	mspec->sel[0].range[0] = l->v.xval;
	mspec->sel[0].range[1] = l->v.xval;
    } else if (l->t == IVEC) {
	mspec->type[0] = SEL_RANGE;
	mspec->sel[0].range[0] = l->v.ivec[0];
	mspec->sel[0].range[1] = l->v.ivec[1];
    } else if (l->t == MAT) {
	mspec->type[0] = SEL_MATRIX;
	mspec->sel[0].m = l->v.m;
    } else if (l->t == EMPTY) {
	mspec->type[0] = SEL_ALL;
    } else {
	*err = E_TYPES;
	goto bailout;
    }

    if (r->t == ABSENT) {
	mspec->type[1] = SEL_NULL;
    } else if (r->t == NUM) {
	mspec->type[1] = SEL_RANGE;
	mspec->sel[1].range[0] = r->v.xval;
	mspec->sel[1].range[1] = r->v.xval;
    } else if (r->t == IVEC) {
	mspec->type[1] = SEL_RANGE;
	mspec->sel[1].range[0] = r->v.ivec[0];
	mspec->sel[1].range[1] = r->v.ivec[1];
    } else if (r->t == MAT) {
	mspec->type[1] = SEL_MATRIX;
	mspec->sel[1].m = r->v.m;
    } else if (r->t == EMPTY) {
	mspec->type[1] = SEL_ALL;
    } else {
	*err = E_TYPES;
	goto bailout;
    }

 bailout:
    
    if (*err) {
	free(mspec);
	mspec = NULL;
    }

    return mspec;
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

static NODE *get_submatrix (NODE *l, NODE *r, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *a = NULL;

	if (r->t != MSPEC) {
	    fprintf(stderr, "get_submatrix: couldn't find mspec\n");
	    p->err = E_TYPES;
	    return NULL;
	}

	if (l->t == MAT) {
	    a = matrix_get_submatrix(l->v.m, r->v.mspec, &p->err);
	} else if (l->t == STR) {
	    a = user_matrix_get_submatrix(l->v.str, r->v.mspec, &p->err);
	} else {
	    p->err = E_TYPES;
	}

	if (a != NULL) {
	    if (gretl_matrix_is_scalar(a)) {
		ret = aux_scalar_node(p);
		if (ret != NULL) {
		    ret->v.xval = a->val[0];
		} 
		gretl_matrix_free(a);
	    } else {
		ret = aux_matrix_node(p);
		if (ret == NULL) {
		    gretl_matrix_free(a);
		} else {
		    ret->v.m = a;
		}
	    }
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
	if (l->t == NUM && r->t == NUM) {
	    ret = aux_ivec_node(p, 2);
	    if (ret != NULL) {
		ret->v.ivec[0] = (int) l->v.xval;
		ret->v.ivec[1] = (int) r->v.xval;
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

    if (xna(x)) {
	switch (f) {
	case MISSING:
	    return 1.0;
	case OK:
	case MISSZERO:
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
    case ABS:
	return fabs(x);
    case TOINT:
	return (double) (int) x;
    case CEIL:
	return ceil(x);
    case FLOOR:
	return floor(x);
    case ROUND:
	if (x < 0) {
	    return (x - floor(x) <= 0.5)? floor(x) : ceil(x);
	} else {
	    return (x - floor(x) < 0.5)? floor(x) : ceil(x);
	}
    case SIN:
	return sin(x);
    case COS:
	return cos(x);
    case TAN:
	return tan(x);
    case ASIN:
	return asin(x);
    case ACOS:
	return acos(x);
    case ATAN:
	return atan(x);
    case CNORM:
	return normal_cdf(x);
    case DNORM:
	return normal_pdf(x);
    case QNORM:
	return normal_cdf_inverse(x);
    case GAMMA:
	return cephes_gamma(x);
    case LNGAMMA:
	return cephes_lgamma(x);
    case MISSING:
	return 0.0;
    case OK:
	return 1.0;
    case MISSZERO:
	return x;
    case ZEROMISS:
	return (x == 0.0)? NADBL : x;
    case SQRT:
	y = sqrt(x);
	if (errno) {
	    eval_warning(p, SQRT);
	}
	return y;
    case LOG:
    case LOG10:
    case LOG2:
	y = log(x);
	if (f == LOG10 && !errno) {
	    y /= log(10.0);
	} else if (f == LOG2 && !errno) {
	    y /= log(2.0);
	}
	if (errno) {
	    eval_warning(p, LOG);
	}
	return y;
    case EXP:
	y = exp(x);
	if (errno) {
	    eval_warning(p, EXP);
	}
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

static NODE *apply_series_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dinfo->n);
    int t, t1, t2;

    if (ret != NULL) {
	/* AC: changed for autoreg case, 2007/7/1 */
	t1 = (autoreg(p))? p->obs : p->dinfo->t1;
	t2 = (autoreg(p))? p->obs : p->dinfo->t2;
	const double *xvec = getvec(n, p);

	for (t=t1; t<=t2; t++) {
	    ret->v.xvec[t] = real_apply_func(xvec[t], f, p);
	}
    }

    return ret;
}

static int *node_get_list (NODE *n, parser *p)
{
    int *list = NULL;

    if (n->t == LVEC || n->t == LIST) {
	int *src = (n->t == LVEC)? n->v.ivec : get_list_by_name(n->v.str);

	if (src == NULL) {
	    p->err = E_DATA;
	} else {
	    list = gretl_list_copy(src);
	}
    } else if (n->t == VEC || n->t == NUM) {
	int v = (n->t == USERIES)? n->v.idnum : n->v.xval;

	if (v < 0 || v >= p->dinfo->v) {
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

    return list;
}

static NODE *eval_lcat (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_lvec_node(p, 0);

    if (ret != NULL && starting(p)) {
	int *llist, *rlist = NULL;

	llist = node_get_list(l, p);
	if (llist != NULL) {
	    rlist = node_get_list(r, p);
	}
	if (rlist != NULL) {
	    p->err = gretl_list_add_list(&llist, rlist);
	}
	ret->v.ivec = llist;
	free(rlist);
    }

    return ret;
}

/* check for missing obs in a list of variables */

static NODE *list_ok_func (NODE *n, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dinfo->n);

    if (ret != NULL && starting(p)) {
	const int *list = get_list_by_name(n->v.str);
	int allbad = 0;
	int nseries = 0;
	int i, v, t;
	double x;

	if (list == NULL || list[0] == 0) {
	    return ret;
	}

	for (i=1; i<=list[0]; i++) {
	    v = list[i];
	    if (var_is_scalar(p->dinfo, v)) {
		if (na((*p->Z)[v][0])) {
		    allbad = 1;
		    break;
		}
	    } else {
		nseries++;
	    }
	}

	for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	    if (allbad) {
		ret->v.xvec[t] = 0;
	    } else if (nseries == 0) {
		ret->v.xvec[t] = 1;
	    } else {
		x = 1;
		for (i=1; i<=list[0]; i++) {
		    v = list[i];
		    if (var_is_series(p->dinfo, v) && na((*p->Z)[v][t])) {
			x = 0;
			break;
		    }
		}
		ret->v.xvec[t] = x;
	    }
	}
    }

    return ret;
}

/* functions taking (up to) two scalars as arguments and 
   returning a series result */

static NODE *
series_fill_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dinfo->n);

    if (ret != NULL && starting(p)) {
	const double *vx = NULL;
	double x = 0.0;
	double y = 0.0;
	int v = 0;

	if (f == RPOISSON) {
	    if (l->t == VEC) {
		vx = getvec(l, p);
		v = 1;
	    } else {
		x = l->v.xval;
	    }
	} else if (f == RUNIFORM || f == RNORMAL) {
	    x = (l->t == EMPTY)? NADBL : l->v.xval;
	    y = (r->t == EMPTY)? NADBL : r->v.xval;
	} else {
	    v = l->v.xval;
	}

	switch (f) {
	case RUNIFORM:
	    p->err = gretl_rand_uniform_minmax(ret->v.xvec, 
					       p->dinfo->t1, 
					       p->dinfo->t2,
					       x, y);
	    break;
	case RNORMAL:
	    p->err = gretl_rand_normal_full(ret->v.xvec, 
					    p->dinfo->t1, 
					    p->dinfo->t2,
					    x, y);
	    break;
	case RPOISSON:
	    gretl_rand_poisson(ret->v.xvec, p->dinfo->t1, p->dinfo->t2,
			       (v)? vx : &x, v);
	    break;
	default:
	    break;
	}
    }

    return ret;
}

/* functions taking two series as arguments and returning a scalar
   result */

static NODE *series_2_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = getvec(l, p);
	const double *y = getvec(r, p);

	switch (f) {
	case COR:
	    ret->v.xval = gretl_corr(p->dinfo->t1, p->dinfo->t2, x, y, NULL);
	    break;
	case COV:
	    ret->v.xval = gretl_covar(p->dinfo->t1, p->dinfo->t2, x, y);
	    break;
	default:
	    break;
	}
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
	    const double *x = getvec(l, p);
	    const double *y = getvec(r, p);
	    
	    ret->v.m = gretl_matrix_xtab(p->dinfo->t1, p->dinfo->t2, 
					 x, y, &p->err);
	} else {
	    p->err = E_TYPES;
	}
    }

    return ret;
}

static NODE *object_status (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	const char *s = n->v.str;

	ret->v.xval = NADBL;
	if (f == ISSERIES) {
	    int v = varindex(p->dinfo, s);

	    if (v < p->dinfo->v) {
		ret->v.xval = var_is_series(p->dinfo, v);
	    }
	} else if (f == ISLIST || f == LISTLEN) {
	    int *list = get_list_by_name(s);

	    if (list != NULL) {
		ret->v.xval = (f == ISLIST)? 1.0 : list[0];
	    } else if (f == ISLIST) {
		ret->v.xval = 0;
	    }
	} else if (f == ISSTRING) {
	    ret->v.xval = string_is_defined(s);
	} else if (f == ISNULL) {
	    ret->v.xval = 1;
	    if (varindex(p->dinfo, s) < p->dinfo->v) {
		ret->v.xval = 0.0;
	    } else if (get_matrix_by_name(s)) {
		ret->v.xval = 0.0;
	    } else if (get_list_by_name(s)) {
		ret->v.xval = 0.0;
	    } 
	} else if (f == OBSNUM) {
	    int t = get_observation_number(s, p->dinfo);

	    if (t > 0) {
		ret->v.xval = t;
	    }
	}
    }

    return ret;
}

static int series_get_nobs (int t1, int t2, const double *x)
{
    int t, n = 0;

    for (t=t1; t<=t2; t++) {
	if (!xna(x[t])) n++;
    }

    return n;
}

static int series_get_start (int n, const double *x)
{
    int t;

    for (t=0; t<n; t++) {
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

static int series_get_end (int n, const double *x)
{
    int t;

    for (t=n-1; t>=0; t--) {
	if (!xna(x[t])) {
	    break;
	}
    }

    return t + 1;
}

#define series_cast_optional(f) (f == SD) 

static void cast_to_series (NODE *n, int f, gretl_matrix **tmp, parser *p)
{
    gretl_matrix *m = n->v.m;

    if (gretl_is_null_matrix(m)) {
	p->err = E_DATA;
    } else if (gretl_vector_get_length(m) != p->dinfo->n) {
	if (series_cast_optional(f)) {
	    p->err = 1;
	} else {
	    node_type_error(f, VEC, n, p);
	}
    } else {
	*tmp = m;
	n->v.xvec = m->val;
    }
}

/* functions taking a series as argument and returning a scalar */

static NODE *
series_scalar_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	gretl_matrix *tmp = NULL;
	const double *x;

	if (n->t == MAT) {
	    cast_to_series(n, f, &tmp, p);
	    if (p->err) {
		if (f == SD) {
		    /* offer column s.d. instead */
		    p->err = 0;
		    return matrix_to_matrix_func(n, f, p);
		} else {
		    return NULL;
		}
	    }
	}

	x = getvec(n, p);

	switch (f) {
	case SUM:
	    ret->v.xval = gretl_sum(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case MEAN:
	    ret->v.xval = gretl_mean(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case SD:
	    ret->v.xval = gretl_stddev(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case VCE:
	    ret->v.xval = gretl_variance(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case SST:
	    ret->v.xval = gretl_sst(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case MIN:
	    ret->v.xval = gretl_min(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case MAX: 
	    ret->v.xval = gretl_max(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case MEDIAN:
	    ret->v.xval = gretl_median(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case GINI:
	    ret->v.xval = gretl_gini(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case NOBS:
	    ret->v.xval = series_get_nobs(p->dinfo->t1, p->dinfo->t2, x);
	    break;
	case T1:
	    ret->v.xval = series_get_start(p->dinfo->n, x);
	    break;
	case T2:
	    ret->v.xval = series_get_end(p->dinfo->n, x);
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

/* 
   functions taking a series and a scalar as arguments and returning 
   a scalar
 */

static NODE *
series_scalar_scalar_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	gretl_matrix *tmp = NULL;
	const double *xvec;

	if (l->t == MAT) {
	    if (f == QUANTILE) {
		ret = aux_matrix_node(p);
		if (ret != NULL) {
		    ret->v.m = gretl_matrix_quantiles(l->v.m, r->v.xval,
						      &p->err);
		}
		return ret;
	    } else {
		cast_to_series(l, f, &tmp, p);
		if (p->err) {
		    return NULL;
		}
	    }
	}

	ret = aux_scalar_node(p);
	if (p->err) {
	    return ret;
	}

	xvec = getvec(l, p);

	switch (f) {
	case LRVAR:
	    ret->v.xval = gretl_long_run_variance(p->dinfo->t1, p->dinfo->t2, 
						  xvec, (int) (r->v.xval));
	    break;
	case QUANTILE:
	    ret->v.xval = gretl_quantile(p->dinfo->t1, p->dinfo->t2, xvec, 
					 r->v.xval);
	    break;
	default:
	    break;
	}

	if (l->t == MAT) {
	    l->v.m = tmp;
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static NODE *series_obs (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_scalar_node(p);
    

    if (ret != NULL) {
	const double *x = getvec(l, p);
	char word[16];
	int t;

	if (r->v.xval < 0 || r->v.xval > (double) INT_MAX) {
	    ret->v.xval = NADBL;
	    return ret;
	}

	/* convert to 0-based, and allow for dates */
	t = (int) r->v.xval;
	sprintf(word, "%d", t);
	t = get_t_from_obs_string(word, (const double **) *p->Z, 
				  p->dinfo);

	if (t >= 0 && t < p->dinfo->n) {
	    ret->v.xval = x[t];
	} else {
	    ret->v.xval = NADBL;
	}
    }

    return ret;
}

static NODE *series_movavg (NODE *l, NODE *r, parser *p)
{
    NODE *ret;
    const double *x = getvec(l, p);
    int k = (int) r->v.xval;
    int t, t1, t2;

    ret = aux_vec_node(p, p->dinfo->n);
    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

    for (t=t1; t<=t2; t++) {
	double xs, msum = 0.0;
	int i, s;

	for (i=0; i<k; i++) {
	    s = t - i;
	    if (p->dinfo->structure == STACKED_TIME_SERIES) {
		if (s >= 0 && s < p->dinfo->n && 
		    p->dinfo->paninfo->unit[s] != 
		    p->dinfo->paninfo->unit[t]) {
		    s = -1;
		}
	    }

	    if (s >= 0) {
		xs = x[s];
	    } else {
		xs = NADBL;
	    }

	    if (na(xs)) {
		msum = NADBL;
		break;
	    } else {
		msum += x[s];
	    }
	}

	if (!na(msum)) {
	    ret->v.xvec[t] = msum / k;
	} 
    }

    return ret;
}

static NODE *series_lag (NODE *l, NODE *r, parser *p)
{
    NODE *ret;
    const double *x = getvec(l, p);
    int k = (int) -r->v.xval;
    int t, s, t1, t2;

    ret = aux_vec_node(p, p->dinfo->n);
    if (ret == NULL) {
	return NULL;
    }

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

    for (t=t1; t<=t2; t++) {
	s = t - k;
	if (dated_daily_data(p->dinfo)) {
	    if (s >= 0 && s < p->dinfo->n) {
		while (s >= 0 && xna(x[s])) {
		    s--;
		}
	    }
	} else if (p->dinfo->structure == STACKED_TIME_SERIES) {
	    if (s >= 0 && s < p->dinfo->n && 
		p->dinfo->paninfo->unit[s] != 
		p->dinfo->paninfo->unit[t]) {
		s = -1;
	    }
	}
	if (s >= 0 && s < p->dinfo->n) {
	    ret->v.xvec[t] = x[s];
	} 
    }

    return ret;
}

static NODE *series_sort_by (NODE *l, NODE *r, parser *p)
{
    NODE *ret = aux_vec_node(p, p->dinfo->n);

    if (ret != NULL && starting(p)) {
	if (l->t == VEC && r->t == VEC) {
	    const double *x = getvec(l, p);
	    const double *y = getvec(r, p);

	    p->err = gretl_sort_by(x, y, ret->v.xvec, p->dinfo); 
	} else {
	    p->err = E_TYPES;
	}

	if (p->err) {
	    free(ret);
	    ret = NULL;
	}
    } 

    return ret;
}

static NODE *vector_sort (NODE *l, int f, parser *p)
{
    NODE *ret = (l->t == VEC)? aux_vec_node(p, p->dinfo->n) :
	aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	if (l->t == VEC) {
	    const double *x = getvec(l, p);

	    p->err = sort_series(x, ret->v.xvec, f, p->dinfo); 
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

		    qsort(x, n, sizeof *x, (f == SORT)? gretl_compare_doubles :
			  gretl_inverse_compare_doubles);
		}
	    } else {
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    free(ret);
	    ret = NULL;
	}
    } 

    return ret;
}

static NODE *vector_values (NODE *l, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const double *x = NULL;
	int n = 0;

	if (l->t == VEC) {
	    n = p->dinfo->t2 - p->dinfo->t1 + 1;
	    x = getvec(l, p) + p->dinfo->t1;
	} else if (!gretl_is_null_matrix(l->v.m)) {
	    n = gretl_vector_get_length(l->v.m);
	    x = l->v.m->val;
	}

	if (n > 0 && x != NULL) {
	    ret->v.m = gretl_matrix_values(x, n, &p->err);
	} else {
	    p->err = E_DATA;
	}

	if (p->err) {
	    free(ret);
	    ret = NULL;
	}
    } 

    return ret;
}

/* functions taking a series as argument and returning a series */

static NODE *series_series_func (NODE *l, NODE *r, int f, parser *p)
{
    NODE *ret = NULL;

    if (f == SDIF && !dataset_is_seasonal(p->dinfo)) {
	p->err = E_PDWRONG;
	return NULL;
    } else {
	ret = aux_vec_node(p, p->dinfo->n);
    }

    if (ret != NULL && starting(p)) {
	gretl_matrix *tmp = NULL;
	const double *x;
	double *y;

	if (l->t == MAT) {
	    cast_to_series(l, f, &tmp, p);
	    if (p->err) {
		return NULL;
	    }
	}

	x = getvec(l, p);
	y = ret->v.xvec;

	switch (f) {
	case HPFILT:
	    p->err = hp_filter(x, y, p->dinfo, OPT_NONE);
	    break;
	case BKFILT:
	    p->err = bkbp_filter(x, y, p->dinfo);
	    break;
	case FRACDIF:
	    p->err = fracdiff_series(x, y, r->v.xval, p->dinfo);
	    break;
	case DIF:
	case LDIF:
	case SDIF:
	    p->err = diff_series(x, y, f, p->dinfo); 
	    break;
	case ODEV:
	    p->err = orthdev_series(x, y, p->dinfo); 
	    break;
	case CUM:
	    p->err = cum_series(x, y, p->dinfo); 
	    break;
	case RESAMPLE:
	    p->err = resample_series(x, y, p->dinfo); 
	    break;
	case PMEAN:
	    p->err = panel_mean_series(x, y, p->dinfo); 
	    break;
	case PSD:
	    p->err = panel_sd_series(x, y, p->dinfo); 
	    break;
	case RANKING:
	    p->err = rank_series(x, y, SORT, p->dinfo); 
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

/* application of scalar function to each element of matrix */

static NODE *apply_matrix_func (NODE *n, int f, parser *p)
{
    NODE *ret = aux_matrix_node(p);

    if (ret != NULL && starting(p)) {
	const gretl_matrix *m = n->v.m;
	int i, n = m->rows * m->cols;
	double x;

	if (node_allocate_matrix(ret, m->rows, m->cols, p)) {
	    free_tree(ret, "On error");
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

/* node holding a user-defined variable, either a scalar
   or a series */

static NODE *uvar_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (t->t == USCLR) {
	ret = aux_scalar_node(p);
	if (ret != NULL) {
	    ret->v.xval = (*p->Z)[t->v.idnum][0];
	}
    } else if (t->t == USERIES) {
	ret = t;
	t->t = VEC;
	t->flags |= UVEC_NODE;
    }

    return ret;
}

/* node holding a user-defined matrix */

static NODE *umatrix_node (NODE *t, parser *p)
{
    NODE *ret = matrix_pointer_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.m = get_matrix_by_name(t->v.str);
    }

    return ret;
}

static NODE *loop_index_node (NODE *t, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	ret->v.xval = loop_scalar_read(*t->v.str);
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
    int i, j, k;

    /* check that all rows are the same length */

    if (nelem != r * c) {
	p->err = 1;
    } else if (nsep > 0) {
	k = 0;
	for (i=0; i<m; i++) {
	    n = t->v.bn.n[i];
	    if (n->t == EMPTY) {
		if (i - k != seppos) {
		    p->err = 1;
		    break;
		}
		k = i + 1;
	    }
	}
    }

    if (p->err) {
	pprintf(p->prn, "matrix specification is not coherent\n");
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
		gretl_matrix_set(M, i, j, n->v.xval);
	    }
	}
    }

    return M;
}

static int *full_series_list (const DATAINFO *pdinfo, int *err)
{
    int *list = NULL;
    int i, j, n = 0;

    for (i=1; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    n++;
	}
    }

    if (n == 0) {
	*err = E_DATA;
	return NULL;
    }

    list = gretl_list_new(n);
    if (list == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=1, j=1; i<pdinfo->v; i++) {
	if (var_is_series(pdinfo, i)) {
	    list[j++] = i;
	}
    }

    return list;
}

#define MATRIX_SKIP_MISSING 1

static gretl_matrix *matrix_from_list (NODE *t, parser *p)
{
    gretl_matrix *M;
    int *list = NULL;
    int freelist = 0;

    if (t != NULL) {
	list = get_list_by_name(t->v.str);
	if (list == NULL) {
	    p->err = E_DATA;
	}
    } else {
	list = full_series_list(p->dinfo, &p->err);
	freelist = 1;
    }

    if (p->err) {
	return NULL;
    }

#if MATRIX_SKIP_MISSING
    M = gretl_matrix_data_subset_skip_missing(list, (const double **) *p->Z, 
					      p->dinfo->t1, p->dinfo->t2, 
					      &p->err);
#else
    M = gretl_matrix_data_subset_no_missing(list, (const double **) *p->Z, 
					    p->dinfo->t1, p->dinfo->t2, 
					    &p->err);
#endif

    if (freelist) {
	free(list);
    }

    return M;
}

#define ok_ufunc_sym(s) (s == NUM || s == VEC || s == MAT || \
                         s == LIST || s == U_ADDR || s == DUM || \
                         s == STR || s == EMPTY)

/* evaluate a user-defined function */

static NODE *eval_ufunc (NODE *t, parser *p)
{
    fnargs *args = NULL;
    ufunc *uf = NULL;
    int argc = 0;
    NODE *l = t->v.b2.l;
    NODE *r = t->v.b2.r;
    int i, m = r->v.bn.n_nodes;
    const char *funname = l->v.str;
    int rtype = GRETL_TYPE_NONE;
    NODE *n, *ret = NULL;

    /* find the function */
    uf = get_user_function_by_name(funname);
    if (uf == NULL) {
	fprintf(stderr, "%s: couldn't find a function of this name\n", funname);
	p->err = 1;
	return NULL;
    }

    /* check that the function returns something suitable, if required */
    if (!simple_ufun_call(p)) {
	rtype = user_func_get_return_type(uf);
	if (rtype != GRETL_TYPE_DOUBLE && rtype != GRETL_TYPE_SERIES &&
	    rtype != GRETL_TYPE_MATRIX && rtype != GRETL_TYPE_LIST) {
	    fprintf(stderr, "%s: invalid return type %d\n", funname, rtype);
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    /* check the argument count */
    argc = fn_n_params(uf);
    if (m > argc) {
	pprintf(p->prn, _("Number of arguments (%d) does not "
			  "match the number of\nparameters for "
			  "function %s (%d)"),
		m, funname, argc);
	p->err = 1;
	return NULL;
    }

    /* make an arguments array */
    args = fn_args_new();
    if (args == NULL) {
	fprintf(stderr, "%s: invalid return type %d\n", funname, rtype);
	p->err = E_ALLOC;
	return NULL;
    }

    /* evaluate the function arguments */

    for (i=0; i<m && !p->err; i++) {
	int type = r->v.bn.n[i]->t;

	if (type == USCLR || type == USERIES) {
	    /* let dataset variables through "as is" */
	    n = r->v.bn.n[i];
	} else {
	    n = eval(r->v.bn.n[i], p);
	    if (n == NULL) {
		fprintf(stderr, "%s: failed to evaluate arg %d\n", funname, i); 
	    } else if (!ok_ufunc_sym(n->t)) {
		fprintf(stderr, "%s: node type %d: not OK\n", funname, n->t);
		p->err = E_TYPES;
	    }
	}

	if (p->err) {
	    break;
	}

#if EDEBUG
	fprintf(stderr, "%s: arg[%d] is of type %d\n", funname, i, n->t);
#endif

	if (n->t == U_ADDR) {
	    NODE *u = n->v.b1.b;

	    if (u->t == USCLR) {
		p->err = push_fn_arg(args, GRETL_TYPE_SCALAR_REF, &u->v.idnum);
	    } else if (u->t == USERIES) {
		p->err = push_fn_arg(args, GRETL_TYPE_SERIES_REF, &u->v.idnum);
	    } else if (u->t == UMAT) {
		user_matrix *m = get_user_matrix_by_name(u->v.str);

		p->err = push_fn_arg(args, GRETL_TYPE_MATRIX_REF, m);
	    } else {
		pputs(p->prn, "Wrong type of operand for unary '&'\n");
		p->err = 1;
	    }
	} else if (n->t == DUM) {
	    if (n->v.idnum == DUM_NULL) {
		p->err = push_fn_arg(args, GRETL_TYPE_NONE, NULL);
	    } else {
		p->err = E_TYPES;
	    }
	} else if (n->t == EMPTY) {
	    p->err = push_fn_arg(args, GRETL_TYPE_NONE, NULL);
	} else if (n->t == NUM) {
	    p->err = push_fn_arg(args, GRETL_TYPE_DOUBLE, &n->v.xval);
	} else if (n->t == VEC) {
	    p->err = push_fn_arg(args, GRETL_TYPE_SERIES, getvec(n, p));
	} else if (n->t == MAT) {
	    p->err = push_fn_arg(args, GRETL_TYPE_MATRIX, n->v.m);
	} else if (n->t == LIST) {
	    p->err = push_fn_arg(args, GRETL_TYPE_LIST, n->v.str);
	} else if (n->t == STR) {
	    p->err = push_fn_arg(args, GRETL_TYPE_STRING, n->v.str);
	} else if (n->t == USCLR || n->t == USERIES) {
	    p->err = push_fn_arg(args, GRETL_TYPE_UVAR, &n->v.idnum);
	} 

	if (p->err) {
	    fprintf(stderr, "eval_ufunc: error evaluating arg %d\n", i);
	}
    }

    /* try sending args to function */

    if (!p->err) {
	char *descrip = NULL;
	char **pdescrip = NULL;
	double xret = NADBL;
	double *Xret = NULL;
	gretl_matrix *mret = NULL;
	char *sret = NULL;
	void *retp = NULL;

	if (rtype == GRETL_TYPE_DOUBLE) {
	    retp = &xret;
	} else if (rtype == GRETL_TYPE_SERIES) {
	    retp = &Xret;
	} else if (rtype == GRETL_TYPE_MATRIX) {
	    retp = &mret;
	} else if (rtype == GRETL_TYPE_LIST) {
	    retp = &sret;
	}

	if ((p->flags & P_UFRET) && 
	    (rtype == GRETL_TYPE_DOUBLE || rtype == GRETL_TYPE_SERIES)) {
	    /* pick up description of generated var, if any */
	    pdescrip = &descrip;
	}

	p->err = gretl_function_exec(uf, args, rtype, p->Z, p->dinfo, 
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
			free(ret->v.str);
		    }
		    ret->v.str = sret;
		}
	    }
	}

	if (descrip != NULL) {
	    strcpy(p->lh.label, descrip);
	    free(descrip);
	}
    }

    fn_args_free(args);

    return ret;
}

static void n_args_error (int k, int n, const char *s, parser *p)
{
    pprintf(p->prn, _("Number of arguments (%d) does not "
		      "match the number of\nparameters for "
		      "function %s (%d)"),
	    k, "mshape", n);
    p->err = 1;
}

/* evaluate a built-in function that has more than two arguments */

static NODE *eval_nargs_func (NODE *t, parser *p)
{
    NODE *e, *n = t->v.b1.b;
    NODE *ret = NULL;
    gretl_matrix *m = NULL;
    int i, k = n->v.bn.n_nodes;

    if (t->t == MSHAPE) {
	gretl_matrix *A = NULL;
	int r = 0, c = 0;

	if (k != 3) {
	    n_args_error(k, 3, "mshape", p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    e = eval(n->v.bn.n[i], p);
	    if (e == NULL) {
		fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i); 
	    } else if (i == 0) {
		if (e->t != MAT) {
		    p->err = E_TYPES;
		} else {
		    A = e->v.m;
		}
	    } else if (e->t != NUM) {
		p->err = E_TYPES;
	    } else if (i == 1) {
		r = e->v.xval;
	    } else {
		c = e->v.xval;
	    }
	}

	if (!p->err) {
	    m = gretl_matrix_shape(A, r, c);
	}

    } else if (t->t == SVD) {
	gretl_matrix *A = NULL;
	const char *lname = NULL;
	const char *rname = NULL;
	const char **name;

	if (k != 3) {
	    n_args_error(k, 3, "sdv", p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    if (i == 0) {
		e = eval(n->v.bn.n[i], p);
		if (e == NULL) {
		    fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
		} else {
		    A = e->v.m;
		}
	    } else {
		name = (i == 1)? &lname : &rname;
		e = n->v.bn.n[i];
		if (e->t == EMPTY) {
		    *name = "null";
		} else if (e->t == U_ADDR) {
		    e = e->v.b1.b;
		    if (e->t == UMAT) {
			*name = e->v.str;
		    } else {
			p->err = E_TYPES;
			strcpy(gretl_errmsg, "Expected the address of a matrix");
		    }
		} else {
		    p->err = E_TYPES;
		}
	    }
	}

	if (!p->err) {
	    m = user_matrix_SVD(A, lname, rname, &p->err);
	}
    } else if (t->t == MOLS) {
	gretl_matrix *Y = NULL;
	gretl_matrix *X = NULL;
	const char *Uname = NULL;

	if (k != 3) {
	    n_args_error(k, 3, "sdv", p);
	} 

	for (i=0; i<k && !p->err; i++) {
	    if (i < 2) {
		e = eval(n->v.bn.n[i], p);
		if (e == NULL) {
		    fprintf(stderr, "eval_nargs_func: failed to evaluate arg %d\n", i);
		} else if (i == 0) {
		    Y = e->v.m;
		} else {
		    X = e->v.m;
		}
	    } else {
		e = n->v.bn.n[i];
		if (e->t == EMPTY) {
		    Uname = "null";
		} else if (e->t == U_ADDR) {
		    e = e->v.b1.b;
		    if (e->t == UMAT) {
			Uname = e->v.str;
		    } else {
			p->err = E_TYPES;
			strcpy(gretl_errmsg, "Expected the address of a matrix");
		    }
		} else {
		    p->err = E_TYPES;
		}
	    }
	}

	if (!p->err) {
	    m = user_matrix_ols(Y, X, Uname, &p->err);
	}
    }	

    if (!p->err) {
	ret = aux_matrix_node(p);
    }

    if (!p->err) {
	if (ret->v.m != NULL) {
	    gretl_matrix_free(ret->v.m);
	}
	ret->v.m = m;
    }

    return ret;
}

/* Create a matrix using selected series, or a mixture of series and
   lists, or more than one list.  Note that we can't use an augmented
   list here, because the series are not necessarily members of the
   dataset: they could be auxiliary series.
*/

static gretl_matrix *assemble_matrix (NODE *nn, int nnodes, parser *p)
{
    NODE *n;
    gretl_matrix *m = NULL;
    int *list;
    double **X = NULL;
    int t, T, k = 0;
    int i, j, s;

    for (i=0; i<nnodes; i++) {
	n = nn->v.bn.n[i];
	if (n->t == LIST) {
	    list = get_list_by_name(n->v.str);
	    if (list == NULL) {
		p->err = E_DATA;
		return NULL;
	    } 
	    k += list[0];
	} else if (n->t == VEC) {
	    k++;
	}
    }

    X = malloc(k * sizeof *X);
    if (X == NULL) {
	return NULL;
    }

    s = 0;
    for (i=0; i<nnodes; i++) {
	n = nn->v.bn.n[i];
	if (n->t == LIST) {
	    list = get_list_by_name(n->v.str);
	    for (j=1; j<=list[0]; j++) {
		X[s++] = (*p->Z)[list[j]];
	    }
	} else if (n->t == VEC) {
	    X[s++] = getvec(n, p);
	}
    }

    T = p->dinfo->t2 - p->dinfo->t1 + 1;

    for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	for (i=0; i<k; i++) {
	    if (na(X[i][t])) {
#if MATRIX_SKIP_MISSING
		T--;
		break;
#else
		free(X);
		p->err = E_MISSDATA;
		return NULL;
#endif
	    }
	}
    }

    if (T == 0) {
	free(X);
	p->err = E_DATA;
	return NULL;
    }

    m = gretl_matrix_alloc(T, k);
    if (m == NULL) {
	p->err = E_ALLOC;
    } else {
	int skip;

	i = 0;
	for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	    skip = 0;
	    for (j=0; j<k; j++) {
		if (na(X[j][t])) {
		    skip = 1;
		    break;
		}
	    }
	    if (!skip) {
		for (j=0; j<k; j++) {
		    gretl_matrix_set(m, i, j, X[j][t]);
		}
		if (i == 0) {
		    gretl_matrix_set_t1(m, t);
		} else if (i == T - 1) {
		    gretl_matrix_set_t2(m, t);
		}
		i++;
	    }
	}
    }

    free(X);

    return m;
}

#define ok_matdef_sym(s) (s == NUM || s == VEC || s == EMPTY || \
                          s == DUM || s == LIST)

/* composing a matrix from scalars, series or lists */

static NODE *matrix_def_node (NODE *t, parser *p)
{
    gretl_matrix *M = NULL;
    NODE *nn, *n, *ret = NULL;
    int m = t->v.bn.n_nodes;
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
	nn = aux_mdef_node(p, m);
	if (nn == NULL) {
	    return NULL;
	}
    } else {
	nn = t;
    }

#if EDEBUG
    fprintf(stderr, "Processing MDEF...\n");
#endif

    for (i=0; i<m && !p->err; i++) {
	n = t->v.bn.n[i];
	if (ok_matdef_sym(n->t)) {
	    nn->v.bn.n[i] = n;
	} else {
	    n = eval(n, p);
	    if (p->err) {
		break;
	    }
	    if (ok_matdef_sym(n->t)) {
		if (nn == t) {
		    free_tree(t->v.bn.n[i], "MatDef");
		}
		nn->v.bn.n[i] = n;
	    } else {
		fprintf(stderr, "matrix_def_node: node type %d: not OK\n", n->t);
		p->err = E_TYPES;
		break;
	    }
	}
	if (n->t == NUM) {
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
	}
    }

    for (i=0; i<m; i++) {
	/* forestall double-freeing: null out any aux nodes */
	if (is_aux_node(nn->v.bn.n[i])) {
	    nn->v.bn.n[i] = NULL;
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

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

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
	    xvec = getvec(l, p);
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
	    yvec = getvec(r, p);
	} else if (r->t == NUM) {
	    y = r->v.xval;
	} else {
	    p->err = E_TYPES;
	    return NULL;
	}
    }

    ret = aux_vec_node(p, p->dinfo->n);

    t1 = (autoreg(p))? p->obs : p->dinfo->t1;
    t2 = (autoreg(p))? p->obs : p->dinfo->t2;

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
    int branch;

    branch = (xna(x))? FORK_NONE : (x != 0)? FORK_L : FORK_R;

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
	const double *x = getvec(n, p);
	int t, T = p->dinfo->n;

	ret = aux_vec_node(p, T);
	if (ret != NULL) {
	    for (t=0; t<T; t++) {
		ret->v.xvec[t] = x[t];
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
    } else {
	p->err = E_TYPES;
    }

    return ret;
}

/* Evaluate a ternary "query" expression: (C)? X : Y.  The condition C
   must be a scalar or a series.  The relevant sub-nodes of t are
   named "l" (left, the condition), "m" and "r" (middle and right
   respectively, the two alternates).
*/

static NODE *eval_query (NODE *t, parser *p)
{
    NODE *e, *ret = NULL;
    double *vec = NULL;
    double x = NADBL;

#if EDEBUG
    fprintf(stderr, "eval_query: t=%p, l=%p, m=%p, r=%p\n", 
	    (void *) t, (void *) t->v.b3.l, (void *) t->v.b3.m,
	    (void *) t->v.b3.r);
#endif

    /* evaluate and check the condition */

    e = eval(t->v.b3.l, p);
    if (!p->err) {
	if (e->t != NUM && e->t != VEC) {
	    p->err = E_TYPES;
	} else if (e->t == NUM) {
	    x = e->v.xval;
	} else {
	    vec = getvec(e, p);
	}
    }

    if (p->err) {
	return NULL;
    }

    if (vec != NULL) {
	ret = query_eval_vec(vec, t, p);
    } else {
	ret = query_eval_scalar(x, t, p);
    }

    if (ret != NULL && (ret == t->v.b3.m || ret == t->v.b3.r)) {
	/* forestall double-freeing */
	ret = ternary_return_node(ret, p);
    }

    return ret;
}

#define dvar_scalar(i) (i < R_SCALAR_MAX)
#define dvar_series(i) (i == R_INDEX)

static double dvar_get_value (int i, parser *p)
{
    switch (i) {
    case R_NOBS:
	return p->dinfo->t2 - p->dinfo->t1 + 1;
    case R_NVARS:
	return p->dinfo->v;
    case R_PD:
	return p->dinfo->pd;
    case R_T1:
	return p->dinfo->t1 + 1;
    case R_T2:
	return p->dinfo->t2 + 1;
    case R_TEST_PVAL:
	return get_last_pvalue(p->lh.label);
    case R_TEST_STAT:
	return get_last_test_statistic(p->lh.label);
    case R_TEST_LNL:
	return get_last_lnl(p->lh.label);
    case R_STOPWATCH:
	return gretl_stopwatch();
    case R_NSCAN:
	return n_scanned_items();
    default:
	return NADBL;
    }
}

static double *dvar_get_series (int i, parser *p)
{
    double *x = NULL;

    switch (i) {
    case R_INDEX:
	x = malloc(p->dinfo->n * sizeof *x);
	if (x != NULL) {
	    int t, yr = p->dinfo->structure == TIME_SERIES 
		&& p->dinfo->pd == 1;

	    for (t=0; t<p->dinfo->n; t++) {
		x[t] = (yr)? p->dinfo->sd0 + t : t + 1;
	    }
	}
	break;
    default:
	break;
    }

    return x;
}

static NODE *dollar_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	if (dvar_scalar(t->v.idnum)) {
	    ret = aux_scalar_node(p);
	    if (ret != NULL && starting(p)) {
		ret->v.xval = dvar_get_value(t->v.idnum, p);
	    }
	} else if (dvar_series(t->v.idnum)) {
	    ret = aux_vec_node(p, 0);
	    if (ret != NULL && starting(p)) {
		ret->v.xvec = dvar_get_series(t->v.idnum, p);
	    }
	}
    } else {
	ret = aux_any_node(p);
    }

    return ret;
}

static gretl_matrix *
object_var_get_submatrix (const char *oname, NODE *t, parser *p)
{
    NODE *r = eval(t->v.b2.r, p);
    gretl_matrix *M, *S = NULL;
    int idx;

    if (r == NULL || r->t != MSPEC) {
	if (!p->err) {
	    p->err = E_TYPES;
	}
	return NULL;
    }

    /* the sort of matrix we want (e.g. $coeff) */
    idx = t->v.b2.l->v.idnum;
    M = saved_object_get_matrix(oname, idx, &p->err);

    if (M != NULL) {
	S = matrix_get_submatrix(M, r->v.mspec, &p->err);
	gretl_matrix_free(M);
    }

    return S;
}

static GretlType object_var_type (int idx, const char *oname)
{
    GretlType vtype = GRETL_TYPE_NONE;
    
    if (model_data_scalar(idx)) {
	vtype = GRETL_TYPE_DOUBLE;
    } else if (model_data_series(idx)) {
	vtype = GRETL_TYPE_SERIES;
    } else if (model_data_matrix(idx)) {
	vtype = GRETL_TYPE_MATRIX;
    }
    
    if (idx == M_UHAT || idx == M_YHAT || idx == M_SIGMA) {
	/* could be a matrix */
	GretlObjType otype = gretl_model_get_type(oname);

	if (otype != GRETL_OBJ_EQN) {
	    vtype = GRETL_TYPE_MATRIX;
	}
    }

    return vtype;
}

/* the left-hand subnode holds the name of the object in
   question; on the right is a specification of what we
   want from that object 
*/

static NODE *object_var_node (NODE *t, parser *p)
{
    NODE *ret = NULL;

    if (starting(p)) {
	NODE *r = (t->t == MVAR || t->t == DMSL)? t : t->v.b2.r;
	const char *oname = (t->t == MVAR || t->t == DMSL)?
	    NULL : t->v.b2.l->v.str;
	int mslice = r->t == DMSL;
	GretlType vtype;

	vtype = object_var_type(r->v.idnum, oname);

#if EDEBUG
	fprintf(stderr, "object_var_node: t->t = %d (%s), r->t = %d (%s)\n", 
		t->t, getsymb(t->t, NULL), r->t, getsymb(r->t, NULL));
	fprintf(stderr, "vtype = %d, mslice = %d\n", vtype, mslice);
#endif

	if (vtype == GRETL_TYPE_DOUBLE) {
	    ret = aux_scalar_node(p);
	} else if (vtype == GRETL_TYPE_SERIES) {
	    ret = aux_vec_node(p, 0);
	} else {
	    ret = aux_matrix_node(p);
	}
	if (ret != NULL) {
	    if (vtype == GRETL_TYPE_DOUBLE) {
		ret->v.xval = saved_object_get_scalar(oname, r->v.idnum, &p->err);
	    } else if (vtype == GRETL_TYPE_SERIES) {
		ret->v.xvec = saved_object_get_series(oname, r->v.idnum, p->dinfo,
						      &p->err);
	    } else if (mslice) {
		/* the right-hand subnode needs more work */
		ret->v.m = object_var_get_submatrix(oname, r, p);
	    } else if (vtype == GRETL_TYPE_MATRIX) {
		ret->v.m = saved_object_get_matrix(oname, r->v.idnum, &p->err);
	    } 
	}
    } else {
	ret = aux_any_node(p);
    }
    
    return ret;
}

/* FIXME below: hook this into OVAR, to allow for, e.g.
   "m1.coeff(sqft)" ? */

static NODE *dollar_str_node (NODE *t, parser *p)
{
    NODE *ret = aux_scalar_node(p);

    if (ret != NULL && starting(p)) {
	NODE *l = t->v.b2.l;
	NODE *r = t->v.b2.r;

	ret->v.xval = gretl_model_get_data_element(NULL, l->v.idnum, r->v.str, 
						   p->dinfo, &p->err);

	if (na(ret->v.xval)) {
	    p->err = 1;
	    pprintf(p->prn, _("'%s': invalid argument for %s()\n"), 
		    r->v.str, mvarname(l->v.idnum));
	}
    }

    return ret;
}

static void transpose_matrix_result (NODE *n, parser *p)
{
    if (n == NULL || p->err) {
	return;
    }

    if (n->t == MAT) {
	gretl_matrix *m = n->v.m;

	n->v.m = gretl_matrix_copy_transpose(m);
	if (is_tmp_node(n)) {
	    gretl_matrix_free(m);
	}
	n->flags |= TMP_NODE;
    } else {
	p->err = E_TYPES;
    }
}

static void node_type_error (int ntype, int goodt, NODE *bad, parser *p)
{
    const char *nstr;

    if (ntype == LAG) {
	nstr = (goodt == NUM)? "lag order" : "lag variable";
    } else {
	nstr = getsymb(ntype, NULL);
    }

    pprintf(p->prn, _("Wrong type argument for %s: should be %s"),
	    nstr, typestr(goodt));

    if (bad != NULL) {
	pprintf(p->prn, ", is %s\n", typestr(bad->t));
    } else {
	pputc(p->prn, '\n');
    }

    p->err = E_TYPES;
}

/* core function: evaluate the parsed syntax tree */

static NODE *eval (NODE *t, parser *p)
{  
    NODE *l = NULL, *r = NULL;
    NODE *ret = NULL;

    if (t == NULL) {
	p->err = E_ALLOC;
	return NULL;
    }

    if (evalb2(t->t)) {
	l = eval(t->v.b2.l, p);
	if (l == NULL && p->err == 0) {
	    p->err = 1;
	} else {
	    if (r_return(t->t)) {
		r = t->v.b2.r;
	    } else {
		r = eval(t->v.b2.r, p);
		if (r == NULL && p->err == 0) {
		    p->err = 1;
		}
	    }
	}
    } else if (evalb1(t->t)) {
	l = eval(t->v.b1.b, p);
	if (l == NULL && p->err == 0) {
	    p->err = 1;
	}
    }

    if (p->err) {
	goto bailout;
    }

    switch (t->t) {
    case NUM: 
    case VEC:
    case MAT:
    case STR:
    case DUM:
    case MSPEC:
    case EMPTY:
    case ABSENT:
    case U_ADDR:
    case LIST:
    case LVEC:
	/* terminal symbol: pass on through */
	ret = t;
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
	if (l->t == NUM && r->t == NUM) {
	    ret = scalar_calc(l, r, t->t, p);
	} else if ((l->t == VEC && r->t == VEC) ||
		   (l->t == VEC && r->t == NUM) ||
		   (l->t == NUM && r->t == VEC)) {
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
		   (l->t == VEC || l->t == NUM) && 
		   r->t == STR) {
	    ret = number_string_calc(l, r, t->t, p);
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
	} else if (l->t == MAT && r->t == NUM) {
	    ret = matrix_scalar_calc(l, r, t->t, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;
    case DOTMULT:
    case DOTDIV:
    case DOTPOW:
    case DOTADD:
    case DOTSUB:
    case DOTEQ:
    case DOTGT:
    case DOTLT:
	/* matrix-matrix or matrix-scalar binary operators */
	if ((l->t == MAT && r->t == MAT) ||
	    (l->t == MAT && r->t == NUM) ||
	    (l->t == NUM && r->t == MAT)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else if ((l->t == MAT && r->t == VEC) ||
		   (l->t == VEC && r->t == MAT)) {
	    ret = matrix_series_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case KRON:
    case MCCAT:
    case MRCAT:
    case QFORM:
    case CMULT:
    case CDIV:
    case MRSEL:
    case MCSEL:
	/* matrix-only binary operators (but promote scalars) */
	if ((l->t == MAT || l->t == NUM) && 
	    (r->t == MAT || r->t == NUM)) {
	    ret = matrix_matrix_calc(l, r, t->t, p);
	} else {
	    node_type_error(t->t, MAT, (l->t == MAT)? r : l, p);
	}
	break;
    case MLAG:
	if (l->t == MAT && r->t == NUM) {
	    ret = matrix_lag(l, r, p);
	} else {
	    p->err = E_TYPES; 
	}
	break;
    case U_NEG: 
    case U_POS:
    case U_NOT:
    case ABS:
    case TOINT:
    case CEIL:
    case FLOOR:
    case ROUND:
    case SIN:
    case COS:
    case TAN:
    case ASIN:
    case ACOS:
    case ATAN:
    case LOG:
    case LOG10:
    case LOG2:
    case EXP:
    case SQRT:
    case CNORM:
    case DNORM:
    case QNORM:
    case GAMMA:
    case LNGAMMA:
	/* functions taking one argument, any type */
	if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (l->t == VEC) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == MAT) {
	    ret = apply_matrix_func(l, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case MISSING:
    case MISSZERO:
    case ZEROMISS:
	/* one series or scalar argument needed */
	if (l->t == VEC) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	}
	break;
    case OK:
	/* series, scalar or list argument needed */
	if (l->t == VEC) {
	    ret = apply_series_func(l, t->t, p);
	} else if (l->t == NUM) {
	    ret = apply_scalar_func(l, t->t, p);
	} else if (l->t == LIST) {
	    ret = list_ok_func(l, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	}
	break;
    case MAKEMASK:
	/* one series argument needed: vector output */
	if (l->t == VEC) {
	    ret = make_series_mask(l, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	}
	break;
    case LAG:
    case OBS:
    case MOVAVG:
	/* series on left, scalar on right */
	if (l->t != VEC) {
	    node_type_error(t->t, VEC, l, p);
	} else if (r->t != NUM) {
	    node_type_error(t->t, NUM, r, p);
	} else if (t->t == LAG) {
	    ret = series_lag(l, r, p); 
	} else if (t->t == OBS) {
	    ret = series_obs(l, r, p); 
	} else if (t->t == MOVAVG) {
	    ret = series_movavg(l, r, p); 
	}
	break;
    case MSL:
	/* user matrix plus subspec */
	ret = get_submatrix(l, r, p);
	break;
    case MSL2:
	/* unevaluated matrix subspec */
	ret = mspec_node(l, r, p);
	break;
    case SUBSL:
	/* matrix sub-slice, x:y */
	ret = process_subslice(l, r, p);
	break;
    case LDIF:
    case SDIF:
    case ODEV:
    case HPFILT:
    case BKFILT:
    case FRACDIF:
    case PMEAN:
    case PSD:
    case RANKING:
	/* series argument needed */
	if (l->t == VEC || l->t == MAT) {
	    ret = series_series_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	} 
	break;
    case CUM:
    case DIF:
    case RESAMPLE:
	/* series or matrix argument */
	if (l->t == VEC) {
	    ret = series_series_func(l, r, t->t, p);
	} else if (l->t == MAT) {
	    ret = matrix_to_matrix_func(l, t->t, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	}
	break;
    case SORT:
    case DSORT:
    case VALUES:
	/* series or vector argument needed */
	if (l->t == VEC || l->t == MAT) {
	    if (t->t == VALUES) {
		ret = vector_values(l, p);
	    } else {
		ret = vector_sort(l, t->t, p);
	    }
	} else {
	    node_type_error(t->t, VEC, l, p);
	} 
	break;
    case SUM:
    case MEAN:
    case SD:
    case VCE:
    case SST:
    case MIN:
    case MAX:
    case MEDIAN:
    case GINI:
    case NOBS:
    case T1:
    case T2:
	/* functions taking series arg, returning scalar */
	if (l->t == VEC || l->t == MAT) {
	    ret = series_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	} 
	break;	
    case LRVAR:
    case QUANTILE:	
	/* takes series and scalar arg, returns scalar */
	if (l->t == VEC || l->t == MAT) {
	    if (r->t == NUM) {
		ret = series_scalar_scalar_func(l, r, t->t, p);
	    } else {
		node_type_error(t->t, NUM, r, p);
	    } 
	} else {
	    node_type_error(t->t, VEC, l, p);
	}
	break;
    case RUNIFORM:
    case RNORMAL:
	/* functions taking zero or two scalars as args */
	if (l->t == NUM && r->t == NUM) {
	    ret = series_fill_func(l, r, t->t, p);
	} else if (l->t == EMPTY && r->t == EMPTY) {
	    ret = series_fill_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, NUM, (l->t == NUM)? r : l, p);
	} 
	break;
    case RPOISSON:
	/* one arg: scalar or series */
	if (l->t == NUM || l->t == VEC) {
	    ret = series_fill_func(l, NULL, t->t, p);
	} else {
	    node_type_error(t->t, VEC, l, p);
	} 
	break;
    case COR:
    case COV:
	/* functions taking two series as args */
	if (l->t == VEC && r->t == VEC) {
	    ret = series_2_func(l, r, t->t, p);
	} else {
	    node_type_error(t->t, VEC, (l->t == VEC)? r : l, p);
	} 
	break;
    case MXTAB:
	/* functions taking two series or matrices as args and returning 
	   a matrix */
	if ((l->t == VEC && r->t == VEC) || (l->t == MAT && r->t == MAT)) {
	    ret = mxtab_func(l, r, p);
	} else {
	    node_type_error(t->t, VEC, (l->t == VEC)? r : l, p);
	} 
	break;
    case SORTBY:
	/* takes two series as args, returns series */
	if (l->t == VEC && r->t == VEC) {
	    ret = series_sort_by(l, r, p);
	} else {
	    node_type_error(t->t, VEC, (l->t == VEC)? r : l, p);
	} 
	break;	
    case IMAT:
    case ZEROS:
    case ONES:
    case SEQ:
    case MUNIF:
    case MNORM:
	/* matrix-creation functions */
	if (l->t != NUM || (r != NULL && r->t != NUM)) {
	    node_type_error(t->t, NUM, NULL, p);
	} else {
	    ret = matrix_fill_func(l, r, t->t, p);
	}
	break;
    case SUMC:
    case SUMR:
    case MEANC:
    case MEANR:
    case MCOV:
    case MCORR:
    case CDEMEAN:
    case CHOL:
    case INV:
    case INVPD:
    case GINV:
    case DIAG:
    case TRANSP:
    case TVEC:
    case VECH:
    case UNVECH:
    case UPPER:
    case LOWER:
    case NULLSPC:
    case MEXP:
    case MINC:
    case MAXC:
    case MINR:
    case MAXR:
    case IMINC:
    case IMAXC:
    case IMINR:
    case IMAXR: 
    case FFT:
    case FFTI:
    case POLROOTS:
	/* matrix -> matrix functions */
	if (l->t == MAT) {
	    ret = matrix_to_matrix_func(l, t->t, p);
	} else {
	    node_type_error(t->t, MAT, l, p);
	}
	break;
    case ROWS:
    case COLS:
    case DET:
    case LDET:
    case TRACE:
    case NORM1:
    case INFNORM:
    case RCOND:
    case RANK:
	/* matrix -> scalar functions */
	if (l->t == MAT) {
	    ret = matrix_to_scalar_func(l, t->t, p);
	} else {
	    node_type_error(t->t, MAT, l, p);
	}
	break;
    case MREAD:
	/* string -> matrix functions */
	if (l->t == STR) {
	    ret = string_to_matrix_func(l, t->t, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    case QR:
    case EIGSYM:
    case EIGGEN:
	/* matrix -> matrix functions, with indirect return */
	if (l->t != MAT) {
	    node_type_error(t->t, MAT, l, p);
	} else if (r->t != U_ADDR && r->t != EMPTY) {
	    node_type_error(t->t, U_ADDR, r, p);
	} else {
	    ret = matrix_to_matrix2_func(l, r, t->t, p);
	}
	break;
    case BFGSMAX:
    case FDJAC:
    case MWRITE:
	/* matrix, with string as second arg */
	if (l->t == MAT && r->t == STR) {
	    if (t->t == BFGSMAX) {
		ret = BFGS_maximize(l, r, p);
	    } else if (t->t == FDJAC) {
		ret = numeric_jacobian(l, r, p);
	    } else {
		ret = matrix_csv_write(l, r, p);
	    }
	} else {
	    p->err = E_TYPES;
	} 
	break;
    case PRINCOMP:
	/* matrix, scalar as second arg */
	if (l->t == MAT && r->t == NUM) {
	    ret = matrix_princomp(l, r, p);
	} else {
	    p->err = E_TYPES;
	} 
	break;	
    case MSHAPE:
    case SVD:
    case MOLS:
	/* built-in functions taking more than two args */
	ret = eval_nargs_func(t, p);
	break;
    case USCLR: 
    case USERIES:
	/* user-defined variable */
	ret = uvar_node(t, p);
	break;
    case UMAT:
	/* user-defined matrix */
	ret = umatrix_node(t, p);
	break;
    case OVAR:
    case MVAR:
    case DMSL:
	/* variable "under" user-defined object, or last object */
	ret = object_var_node(t, p);
	break;
    case DMSTR:
	ret = dollar_str_node(t, p);
	break;
    case DVAR: 
	/* dataset "dollar" variable */
	ret = dollar_var_node(t, p);
	break;
    case MDEF:
	/* matrix definition */
	ret = matrix_def_node(t, p);
	break;
    case LOOPIDX:
	ret = loop_index_node(t, p);
	break;
    case OBSNUM:
    case ISSERIES:
    case ISLIST:
    case ISSTRING:
    case ISNULL:
    case LISTLEN:
	if (l->t == STR) {
	    ret = object_status(l, t->t, p);
	} else {
	    node_type_error(t->t, STR, l, p);
	}
	break;
    case CDF:
    case INVCDF:
    case CRIT:
    case PVAL:
    case RANDGEN:
	if (t->v.b1.b->t == FARGS) {
	    ret = eval_pdist(t, p);
	} else {
	    node_type_error(t->t, FARGS, t->v.b1.b, p);
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
    case QUERY:
	ret = eval_query(t, p);
	break;
    case LCAT:
	/* list concatenation */
	if (p->lh.t != LIST) {
	    p->err = E_PARSE;
	} else if (ok_lcat_node(l) && ok_lcat_node(r)) {
	    ret = eval_lcat(l, r, p);
	} else {
	    p->err = E_TYPES;
	}
	break;
    default: 
	printf("EVAL: weird node %s (t->t = %d)\n", getsymb(t->t, NULL),
	       t->t);
	p->err = E_PARSE;
	break;
    }

 bailout:

    if (t->flags & TRANSP_NODE) { /* "starting"? */
	transpose_matrix_result(ret, p);
    }

#if EDEBUG
    fprintf(stderr, "eval (t->t = %03d): returning NODE at %p\n", 
	    t->t, (void *) ret);
#endif

    return ret;
}

/* get the next input character for the lexer */

int parser_getc (parser *p)
{
#if EDEBUG > 1
    fprintf(stderr, "parser_getc: src='%s'\n", p->point);
#endif

    p->ch = 0;

    if (*p->point) {
	p->ch = *p->point;
	p->point += 1;
    }

#if EDEBUG > 1
    fprintf(stderr, "parser_getc: returning '%c'\n", p->ch);
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

/* look ahead to the position of a given character in
   the remaining input stream */

int parser_charpos (parser *p, int c)
{
    int i;

    for (i=0; p->point[i] != '\0'; i++) {
	if (p->point[i] == c) {
	    return i;
	}
    }

    return -1;
}

/* look ahead to the next non-space character and return it */

int parser_next_char (parser *p)
{
    int i;

    for (i=0; p->point[i] != '\0'; i++) {
	if (!isspace(p->point[i])) {
	    return p->point[i];
	}
    }

    return 0;
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

static void printnode (const NODE *t, const parser *p)
{  
    if (t == NULL) {
	pputs(p->prn, "NULL"); 
    } else if (t->t == NUM) {
	if (na(t->v.xval)) {
	    pputs(p->prn, "NA");
	} else {
	    pprintf(p->prn, "%.8g", t->v.xval);
	}
    } else if (t->t == VEC) {
	const double *x = getvec(t, p);
	int i, j = 1;

	if (p->lh.v > 0 && p->lh.v < p->dinfo->v) {
	    pprintf(p->prn, "%s\n", p->dinfo->varname[p->lh.v]);
	}

	for (i=p->dinfo->t1; i<=p->dinfo->t2; i++, j++) {
	    if (na(x[i])) {
		pputs(p->prn, "NA");
	    } else {
		pprintf(p->prn, "%g", x[i]);
	    }
	    if (j % 8 == 0) {
		pputc(p->prn, '\n');
	    } else if (i < p->dinfo->t2) {
		pputc(p->prn, ' ');
	    }
	}
    } else if (t->t == MAT) {
	gretl_matrix_print_to_prn(t->v.m, NULL, p->prn);
    } else if (t->t == USCLR || t->t == USERIES) {
	pprintf(p->prn, "%s", p->dinfo->varname[t->v.idnum]);
    } else if (t->t == UMAT || t->t == UOBJ) {
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
    } else if (func_symb(t->t)) {
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
    } else if (t->t == STR) {
	pprintf(p->prn, "%s", t->v.str);
    } else if (t->t == MDEF) {
	pprintf(p->prn, "{ MDEF }");
    } else if (t->t == DMSTR || t->t == UFUN) {
	printnode(t->v.b2.l, p);
	pputc(p->prn, '(');
	printnode(t->v.b2.r, p);
	pputc(p->prn, ')');
    } else if (t->t == LIST) {
	pputs(p->prn, "LIST");
    } else if (t->t == LVEC) {
	pputs(p->prn, "LVEC");
    } else if (t->t != EMPTY) {
	pputs(p->prn, "weird tree - ");
	printsymb(t->t, p);
    }
}

/* which modified assignment operators of the type '+=' 
   will we accept, when generating a matrix? */

static int ok_matrix_op (int op)
{
    if (op == B_ASN || op == B_ADD || op == B_SUB ||
	op == B_MUL || op == B_DIV || 
	op == INC || op == DEC) {
	return 1;
    } else {
	return 0;
    }
}

static int ok_list_op (int op)
{
    if (op == B_ASN || op == B_ADD || op == B_SUB) {
	return 1;
    } else {
	return 0;
    }
}

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
	switch (s[0]) {
	case '+': return B_ADD;
	case '-': return B_SUB;
	case '*': return B_MUL;
	case '/': return B_DIV;
	case '^': return B_POW;
	case '&': return B_AND;
	case '|': return B_OR;
	}
    }

    return 0;
}

/* extract a substring [...] from the left-hand side
   of a genr expression */

static void get_lhs_substr (char *str, parser *p)
{
    char *s = strchr(str, '[');
    char *q;

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

/* given a substring [...], parse and evaluate it as a
   sub-matrix specification */

static void get_lh_mspec (parser *p)
{
    parser subp;
    char *s;

    s = malloc(strlen(p->lh.substr) + 3);
    if (s == NULL) {
	p->err = E_ALLOC;
	return;
    }

    sprintf(s, "[%s]", p->lh.substr);
    parser_init(&subp, s, p->Z, p->dinfo, p->prn, P_SLICE);

#if EDEBUG
    fprintf(stderr, "subp.input='%s'\n", subp.input);
#endif

    subp.tree = msl_node_direct(&subp);
    p->err = subp.err;

    if (subp.tree != NULL) {
	parser_aux_init(&subp);
	subp.ret = eval(subp.tree, &subp);

	if (subp.err) {
	    printf("Error in subp eval = %d\n", subp.err);
	    p->err = subp.err;
	} else {
	    p->lh.mspec = subp.ret->v.mspec;
	    subp.ret->v.mspec = NULL;
	}

	parser_free_aux_nodes(&subp);
	gen_cleanup(&subp);
    } 

    free(s);
}

/* check validity of "[...]" on the LHS, and evaluate
   the expression if needed */

static void process_lhs_substr (parser *p)
{
    if (p->lh.t == NUM) {
	p->err = E_TYPES;
    } else if (p->lh.t == VEC) {
	p->lh.obs = get_t_from_obs_string(p->lh.substr, (const double **) *p->Z, 
					  p->dinfo); 
	if (p->lh.obs < 0) {
	    p->err = E_PARSE; /* FIXME message */
	} else {
	    p->lh.t = NUM;
	}
    } else if (p->lh.t == MAT) {
	get_lh_mspec(p);
    }
}

#if EDEBUG
static void parser_print_result (parser *p, PRN *prn)
{
    if (p->targ == NUM || p->targ == VEC) {
	int list[2] = { 1, p->lh.v };

	printdata(list, NULL, (const double **) *p->Z, p->dinfo, OPT_NONE, prn);
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

    if (n > 0) {
	for (i=0; i<n && !p->err; i++) {
	    if (S[i] != NULL) {
		if (p->targ == MAT) {
		    gretl_matrix *m = gretl_null_matrix_new();

		    if (m == NULL) {
			p->err = E_ALLOC;
		    } else {
			p->err = user_matrix_add(m, S[i]);
		    }
		} else {
		    if (p->targ == NUM) {
			p->err = dataset_add_scalar(p->Z, p->dinfo);
		    } else if (p->targ == VEC) {
			p->err = dataset_add_series(1, p->Z, p->dinfo);
		    }
		    if (!p->err) {
			v = p->dinfo->v - 1;
			strcpy(p->dinfo->varname[v], S[i]);
		    }
		} 
	    }
	}
    }

    free_strings_array(S, n);
}

/* create a dummy node to facilitate (a) printing an
   existing variable, or (b) incrementing or decrementing
   that variable
*/

static NODE *lhs_copy_node (parser *p)
{
    NODE *n = malloc(sizeof *n);

    if (n == NULL) {
	return NULL;
    }

    n->t = p->targ;
    n->flags = 0;

    if (p->targ == NUM) {
	n->v.xval = (*p->Z)[p->lh.v][0];
    } else if (p->targ == VEC) {
	n->v.xvec = (*p->Z)[p->lh.v];
    } else {
	n->v.m = p->lh.m0;
    }

    return n;
}

/* The expression supplied for evaluation does not contain an '=': can
   we parse it as a declaration of a new variable, or as an implicit
   request to print the value of an existing variable?
*/

static void parser_try_print (parser *p)
{
    if (p->lh.v == 0 && p->lh.m0 == NULL) {
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

    n = strcspn(s, "+-([= ");

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

/* process the left-hand side of a genr formula */

static void pre_process (parser *p, int flags)
{
    const char *s = p->input;
    char test[GENSTRLEN];
    char opstr[3] = {0};
    int newvar = 1;

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
    } else if (!strncmp(s, "scalar ", 7)) {
	p->targ = NUM;
	s += 7;
    } else if (!strncmp(s, "series ", 7)) {
	p->targ = VEC;
	s += 7;
    } else if (!strncmp(s, "matrix ", 7)) {
	p->targ = MAT;
	s += 7;
    } else if (!strncmp(s, "list ", 5)) {
	p->targ = LIST;
	s += 5;
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

    /* grab LHS obs string or matrix slice if present */
    if (strchr(test, '[') != NULL) {
	get_lhs_substr(test, p);
	if (p->err) {
	    return;
	}
    }

#if EDEBUG
    fprintf(stderr, "LHS: %s", test);
    if (p->lh.substr != NULL) {
	fprintf(stderr, "[%s]\n", p->lh.substr);
    } else {
	fputc('\n', stderr);
    }
#endif

    if (strlen(test) > VNAMELEN - 1) {
	pprintf(p->prn, "'%s': name is too long and will be truncated\n", test);
	test[VNAMELEN - 1] = '\0';
    }

    /* find out if the LHS var already exists, and if
       so, what type it is */
    p->lh.v = varindex(p->dinfo, test);
    if (p->lh.v >= p->dinfo->v) {
	/* not a variable: try a matrix? */
	p->lh.v = 0;
	p->lh.m0 = get_matrix_by_name(test);
	if (p->lh.m0 != NULL) {
	    p->lh.t = MAT;
	    newvar = 0;
	} else if (get_list_by_name(test)) {
	    p->lh.islist = 1;
	    p->lh.t = LIST;
	    newvar = 0;
	}
    } else if (var_is_scalar(p->dinfo, p->lh.v)) {
	p->lh.t = NUM;
	newvar = 0;
    } else if (var_is_series(p->dinfo, p->lh.v)) {
	p->lh.t = VEC;
	newvar = 0;
    }

    /* if pre-existing var, check for const-ness */
    if (!newvar && (p->lh.t == NUM || p->lh.t == VEC)) {
	if (var_is_const(p->dinfo, p->lh.v)) {
	    p->err = overwrite_err(p->dinfo, p->lh.v);
	    return;
	}
    }

    /* if new public variable, check name for legality */
    if (newvar && !(flags & P_PRIVATE)) {
	p->err = check_varname(test);
	if (p->err) {
	    return;
	}
    }

    if (p->lh.substr != NULL) {
	process_lhs_substr(p);
	if (p->err) {
	    return;
	}
    }

    if (p->lh.t != UNK) {
	if (p->targ == UNK) {
	    /* when a type is not specified, set from existing
	       variable, if present */
	    p->targ = p->lh.t;
	} else if (p->targ != p->lh.t) {
	    /* don't overwrite one type with another */
	    p->err = E_TYPES;
	    return;
	}
    }

    strcpy(p->lh.name, test);

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
	/* error message */
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
	sprintf(gretl_errmsg, _("'%s' : not implemented for matrices"), opstr);
	p->err = E_PARSE;
	return;
    }

    /* lists: same story as matrices */
    if (p->lh.t == LIST && !ok_list_op(p->op)) {
	sprintf(gretl_errmsg, _("'%s' : not implemented for lists"), opstr);
	p->err = E_PARSE;
	return;
    }	

    /* advance past operator */
    s += strlen(opstr);
    while (isspace(*s)) s++;

    /* set starting point for RHS parser, and also
       for a possible label */
    p->point = p->rhs = s;

    /* increment/decrement operators */
    if ((p->op == INC || p->op == DEC) && *s != '\0') {
	p->err = E_PARSE;
    }

    /* pointer for use in label */
    if (p->op == B_ASN) {
	p->rhs = s;
    } 
}

/* tests for saving variable */

static int non_scalar_matrix (const NODE *r)
{
    return r->t == MAT && (r->v.m->rows != 1 || r->v.m->cols != 1);
}

static int series_compatible (const gretl_matrix *m,
			      const DATAINFO *pdinfo)
{
    int n = gretl_vector_get_length(m);
    int T = pdinfo->t2 - pdinfo->t1 + 1;

    return n == T || n == pdinfo->n || n == 1 || 
	(m->t1 > 0 && m->t1 + n <= pdinfo->n);
}

static int has_missvals (const double *x, const DATAINFO *pdinfo)
{
    int t;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (xna(x[t])) {
	    return 1;
	}
    }

    return 0;
}

static void gen_check_errvals (parser *p)
{
    if (p->ret == NULL || 
	(p->ret->t == VEC && getvec(p->ret, p) == NULL)) {
	return;
    }

    if (p->ret->t == NUM) {
	if (!isfinite(p->ret->v.xval)) {
	    p->ret->v.xval = NADBL;
	    p->warn = 1;
	}
    } else if (p->ret->t == VEC) {
	double *x = getvec(p->ret, p);
	int t;

	for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
	    if (!isfinite(x[t])) {
		x[t] = NADBL;
		p->warn = 1;
	    }
	}
    }
}

static gretl_matrix *list_to_matrix (const char *name, int *err)
{
    const int *list = get_list_by_name(name);
    gretl_matrix *v = NULL;
    int i;

    if (list == NULL) {
	*err = E_DATA;
	return NULL;
    }

    if (list[0] > 0) {
	v = gretl_vector_alloc(list[0]);
	if (v == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	} 
	for (i=0; i<list[0]; i++) {
	    v->val[i] = list[i+1];
	}
    } else {
	*err = E_DATA;
    }

    return v;
}

static gretl_matrix *grab_or_copy_matrix_result (parser *p)
{
    NODE *r = p->ret;
    gretl_matrix *m = NULL;

    if (r->t == NUM) {
	/* result was a scalar, not a matrix */
	m = gretl_matrix_alloc(1, 1);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} else {
	    m->val[0] = r->v.xval;
	}	
    } else if (r->t == VEC) {
	/* result was a series, not a matrix */
	int i, n = p->dinfo->t2 - p->dinfo->t1 + 1;
	const double *x = getvec(r, p);

	m = gretl_column_vector_alloc(n);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} else {
	    for (i=0; i<n; i++) {
		m->val[i] = x[i + p->dinfo->t1];
	    }
	}
    } else if (r->t == MAT && is_tmp_node(r)) {
	/* result r->v.m is newly allocated, steal it */
#if EDEBUG
	fprintf(stderr, "matrix result (%p) is tmp, stealing it\n", 
		(void *) r->v.m);
#endif
	m = r->v.m;
	r->v.m = NULL; /* avoid double-freeing */
    } else if (r->t == MAT) {
	/* r->v.m is an existing user matrix, copy it */
#if EDEBUG
	fprintf(stderr, "matrix result (%p) is pre-existing, copying it\n",
		(void *) r->v.m);
#endif
	m = gretl_matrix_copy(r->v.m);
	if (m == NULL) {
	    p->err = E_ALLOC;
	}
    } else if (r->t == LIST) {
	m = list_to_matrix(r->v.str, &p->err);
    } else {
	fprintf(stderr, "Looking for matrix, but r->t = %d\n", r->t);
	p->err = E_TYPES;
    }

    return m;
}

/* generating a matrix, no pre-existing matrix of that name */

static gretl_matrix *matrix_from_scratch (parser *p, int tmp)
{
    gretl_matrix *m = NULL;

    if (p->ret->t == NUM) {
	m = gretl_matrix_alloc(1, 1);
	if (m == NULL) {
	    p->err = E_ALLOC;
	} else {
	    m->val[0] = p->ret->v.xval;
	}
    } else {
	m = grab_or_copy_matrix_result(p);
    }

    if (!tmp && !p->err) {
	int adj = p->lh.m0 != NULL;

	p->err = user_matrix_add(m, p->lh.name);
	if (adj) {
	    p->lh.m0 = m; /* make the lh matrix pointer valid */
	}
    }

    p->lh.m1 = m;
    
    return m;
}

static gretl_matrix *copy_old_matrix (parser *p)
{
    gretl_matrix *old = get_matrix_by_name(p->lh.name);
    gretl_matrix *m;

    if (old == NULL) {
	p->err = E_UNKVAR;
	return NULL;
    }

    m = gretl_matrix_copy(old);
    if (m == NULL) {
	p->err = E_ALLOC;
    }

    return m;
}

static int LHS_matrix_reusable (parser *p)
{
    gretl_matrix *m = p->lh.m0;
    int ok = 0;

    if (p->ret->t == NUM) {
	ok = (m->rows == 1 && m->cols == 1);
    } else if (p->ret->t == VEC) {
	int T = p->dinfo->t2 - p->dinfo->t1 + 1;

	ok = (m->rows == T && m->cols == 1);
    } else if (p->ret->t == MAT) {
	gretl_matrix *rm = p->ret->v.m;

	ok = (m->rows == rm->rows && m->cols == rm->cols);
    }

    return ok;
}

/* generating a matrix: there's a pre-existing LHS matrix */

static void assign_to_matrix (parser *p)
{
    gretl_matrix *m;

    if (LHS_matrix_reusable(p)) {
	/* result is conformable with original matrix */
	m = p->lh.m0;
	if (p->ret->t == NUM) {
	    m->val[0] = p->ret->v.xval;
	} else if (p->ret->t == VEC) {
	    const double *x = getvec(p->ret, p);
	    int i, s = p->dinfo->t1;

	    for (i=0; i<m->rows; i++) {
		m->val[i] = x[s++];
	    }
	} else {
	    gretl_matrix_copy_values(m, p->ret->v.m);
	}
    } else {
	/* replace the old matrix with result */
	m = grab_or_copy_matrix_result(p);
	p->err = user_matrix_replace_matrix_by_name(p->lh.name, m);
    }

    p->lh.m1 = m;
}

/* assigning to an existing (whole) LHS matrix, but using '+=' or
   some such modified/inflected assignment */

static void assign_to_matrix_mod (parser *p)
{
    gretl_matrix *m = NULL;

    if (p->ret->t == NUM) {
	/* copy original matrix and calculate */
	m = copy_old_matrix(p);
	if (m != NULL) {
	    int i, n = m->rows * m->cols;

	    for (i=0; i<n; i++) {
		m->val[i] = xy_calc(m->val[i], p->ret->v.xval, p->op, p);
	    }
	}
    } else {
	/* get a pointer to original matrix and calculate */
	gretl_matrix *a = get_matrix_by_name(p->lh.name);
	gretl_matrix *b = NULL;

	if (a == NULL) {
	    p->err = E_UNKVAR;
	} else {
	    b = matrix_from_scratch(p, 1);
	}

	if (a != NULL && b != NULL) {
	    m = real_matrix_calc(a, b, p->op, &p->err);
	}
	gretl_matrix_free(b);
    }

    if (!p->err) {
	p->err = user_matrix_replace_matrix_by_name(p->lh.name, m);
	p->lh.m1 = m;
    }
}

/* replacing a sub-matrix of the original LHS matrix, by
   either straight or modified assignment */

static void matrix_edit (parser *p)
{
    gretl_matrix *m = NULL;

    if (p->ret->t != NUM) {
	m = grab_or_copy_matrix_result(p);
	if (m == NULL) {
	    return;
	}
    }

#if EDEBUG
    fprintf(stderr, "matrix_edit: m = %p\n", (void *) m);
#endif

    if (p->op != B_ASN || p->ret->t == NUM) {
	/* doing '+=' or some such: new submatrix 'b' must
	   be calculated using original submatrix 'a' and
	   generated matrix 'm'.
	*/
	gretl_matrix *a = NULL;
	gretl_matrix *b = NULL;

	a = user_matrix_get_submatrix(p->lh.name, p->lh.mspec, &p->err);
	if (!p->err) {
	    if (p->ret->t == NUM) {
		int i, n = a->rows * a->cols;

		for (i=0; i<n; i++) {
		    a->val[i] = xy_calc(a->val[i], p->ret->v.xval, p->op, p);
		}
		m = a; /* preserve modified submatrix */
	    } else {
		b = real_matrix_calc(a, m, p->op, &p->err);
		gretl_matrix_free(a);
		gretl_matrix_free(m);
		m = b; /* replace 'm' with fully computed result */
	    }
	}
    } 

    if (!p->err) {
	/* write new submatrix into place */
	p->err = user_matrix_replace_submatrix(p->lh.name, m,
					       p->lh.mspec);
	gretl_matrix_free(m);
	p->ret->v.m = NULL; /* ?? */
	p->lh.m1 = get_matrix_by_name(p->lh.name);
    }
}

static int edit_list (parser *p)
{
    int tmplist[2] = {1,0};
    const int *list = NULL;
    int v = -1;

    if (p->lh.islist == 0) {
	/* creating a list from scratch */
	if (p->ret->t == LIST) {
	    if (strcmp(p->lh.name, p->ret->v.str)) {
		rename_saved_list(p->ret->v.str, p->lh.name);
	    }
	    return 0;
	} 
    }

    if (p->ret->t == LIST) {
	list = get_list_by_name(p->ret->v.str);
    } else if (p->ret->t == LVEC) {
	list = p->ret->v.ivec;
    } else if (p->ret->t == VEC) {
	v = p->ret->v.idnum;
    } else if (p->ret->t == NUM) {
	v = p->ret->v.xval;
    }

    if (list == NULL) {
	/* v should be ID number of series */
	if (v < 0 || v >= p->dinfo->v) {
	    p->err = E_UNKVAR;
	} else if (var_is_scalar(p->dinfo, v)) {
	    p->err = E_TYPES;
	} else {
	    list = tmplist;
	    tmplist[1] = v;
	}
    }

    if (!p->err) {
	if (p->lh.islist == 0) {
	    remember_list(list, p->lh.name, NULL);
	} else if (p->op == B_ASN) {
	    p->err = replace_list_by_name(p->lh.name, list);
	} else if (p->op == B_ADD) {
	    p->err = append_to_list_by_name(p->lh.name, list);
	} else if (p->op == B_SUB) {
	    p->err = subtract_from_list_by_name(p->lh.name, list);
	} else {
	    p->err = E_TYPES;
	}
    }

    return p->err;
}

static int matrix_missvals (const gretl_matrix *m)
{
    if (m == NULL) {
	return 1;
    } else {
	int i, n = m->rows * m->cols;

	for (i=0; i<n; i++) {
	    if (na(m->val[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

#define scalar_matrix(n) (n->t == MAT && n->v.m->rows == 1 && \
			  n->v.m->cols == 1)

#define ok_return_type(t) (t == NUM || t == VEC || t == MAT || \
			   t == LIST || t == USERIES || t == LVEC)

static int gen_check_return_type (parser *p)
{
    NODE *r = p->ret;

    if (r == NULL) {
	fprintf(stderr, "gen_check_return_type: p->ret = NULL!\n");
	return (p->err = E_DATA);
    }

    if (p->dinfo->n == 0 && r->t != MAT && r->t != NUM) {
	gretl_errmsg_set(_("No dataset is in place"));
	return (p->err = E_DATA);
    }

    if (!ok_return_type(r->t)) {
	return (p->err = E_TYPES);
    }

    if (p->targ == NUM) {
	/* result must be scalar or 1 x 1 matrix */
	if (r->t == VEC || r->t == LIST || non_scalar_matrix(r)) {
	    p->err = E_TYPES;
	} 
    } else if (p->targ == VEC) {
	/* error if result is matrix of wrong dim */
	if ((r->t == MAT && !series_compatible(r->v.m, p->dinfo)) ||
	    r->t == LIST) {
	    p->err = E_TYPES;
	}
    } else if (p->targ == MAT) {
	/* error if result contains NAs */
	if (r->t == VEC && has_missvals(getvec(r, p), p->dinfo)) {
	    p->err = E_MISSDATA;
	} else if (r->t == NUM && xna(r->v.xval)) {
	    p->err = E_MISSDATA;
	} else if (r->t == MAT && matrix_missvals(r->v.m)) {
	    p->err = E_MISSDATA;
	}
    } else if (p->targ == LIST) {
	if (!ok_lcat_node(r)) {
	    p->err = E_TYPES;
	} 
    } else {
	/* target type was not specified: set it now, based
	   on the type of the object we computed */
	p->targ = r->t;
    }

    return p->err;
}

/* allocate storage if saving scalar or series to dataset: 
   lh.v == 0 means that the LHS variable does not already 
   exist
*/

static int gen_allocate_storage (parser *p)
{
    if (p->targ == NUM && p->lh.v == 0) {
	p->err = dataset_add_scalar(p->Z, p->dinfo);
	if (!p->err) {
	    p->lh.v = p->dinfo->v - 1;
#if EDEBUG
	    fprintf(stderr, "gen_allocate_storage: added scalar #%d (%s)\n",
		    p->lh.v, p->lh.name);
#endif
	}
    } else if (p->targ == VEC && p->lh.v == 0) {
	p->err = dataset_add_series(1, p->Z, p->dinfo);
	if (!p->err) {
	    int t;

	    p->lh.v = p->dinfo->v - 1;
	    for (t=0; t<p->dinfo->n; t++) {
		(*p->Z)[p->lh.v][t] = NADBL;
	    }
#if EDEBUG
	    fprintf(stderr, "gen_allocate_storage: added series #%d (%s)\n",
		    p->lh.v, p->lh.name);
#endif
	}
    }

    return p->err;
}

static int save_generated_var (parser *p, PRN *prn)
{
    NODE *r = p->ret;
    double **Z = NULL;
    int t, v;

    /* test for type mismatch errors */
    gen_check_return_type(p);
    if (p->err) {
	return p->err;
    }

#if EDEBUG
    fprintf(stderr, "save_generated_var: targ = %d, ret = %d, op = %d\n",
	    p->targ, p->ret->t, p->op);
#endif

    /* allocate dataset storage, if needed */
    gen_allocate_storage(p);
    if (p->err) {
	return p->err;
    }

    /* put the generated data into place */
    Z = *p->Z;
    v = p->lh.v;
    
    if (p->targ == NUM) {
	/* writing a scalar */
	t = p->lh.obs;
	if (r->t == NUM) {
	    Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, p);
	} else if (r->t == MAT) {
	    Z[v][t] = xy_calc(Z[v][t], r->v.m->val[0], p->op, p);
	}
	strcpy(p->dinfo->varname[v], p->lh.name);
    } else if (p->targ == VEC) {
	/* writing a series */
	if (r->t == NUM) {
	    for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) { 
		Z[v][t] = xy_calc(Z[v][t], r->v.xval, p->op, p);
	    }
	} else if (r->t == VEC) {
	    const double *x = getvec(r, p);
	    int t1 = p->dinfo->t1;

	    if (autoreg(p) && p->op == B_ASN) {
		while (xna(x[t1]) && t1 <= p->dinfo->t2) {
		    t1++;
		}
	    }
	    for (t=t1; t<=p->dinfo->t2; t++) {
		Z[v][t] = xy_calc(Z[v][t], x[t], p->op, p);
	    }
	} else if (r->t == MAT) {
	    const gretl_matrix *m = r->v.m;
	    int s, k = gretl_vector_get_length(m);

	    if (k == 1) {
		/* result is effectively a scalar */
		for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[0], p->op, p);
		}
	    } else if (k == p->dinfo->n) {
		/* treat result as full-length series */
		for (t=p->dinfo->t1; t<=p->dinfo->t2; t++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[t], p->op, p);
		}
	    } else if (k == p->dinfo->t2 - p->dinfo->t1 + 1 && m->t1 == 0) {
		/* treat as series of current sample length */
		for (t=p->dinfo->t1, s=0; t<=p->dinfo->t2; t++, s++) {
		    Z[v][t] = xy_calc(Z[v][t], m->val[s], p->op, p);
		}
	    } else {
		/* align using m->t1 */
		for (t=m->t1; t<m->t1 + k && t<=p->dinfo->t2; t++) {
		    if (t >= p->dinfo->t1) {
			Z[v][t] = xy_calc(Z[v][t], m->val[t - m->t1], p->op, p);
		    }
		}
	    }
	}
	strcpy(p->dinfo->varname[v], p->lh.name);
#if EDEBUG
	fprintf(stderr, "var %d: gave name '%s'\n", v, p->lh.name);
#endif
    } else if (p->targ == MAT) {
	/* writing a matrix */
	if (p->lh.m0 == NULL) {
	    /* no pre-existing LHS: substr must be NULL */
	    matrix_from_scratch(p, 0);
	} else if (p->lh.substr == NULL && p->op == B_ASN) {
	    /* unmodified assignment to existing matrix */
	    assign_to_matrix(p);
	} else if (p->lh.substr == NULL) {
	    /* modified assignment to whole existing matrix */
	    assign_to_matrix_mod(p);
	} else {
	    /* assignment to submatrix of original */
	    matrix_edit(p);
	}
    } else if (p->targ == LIST) {
	edit_list(p);
    }

#if EDEBUG /* print results */
    if (!p->err) {
	parser_print_result(p, prn);
    }
#endif

    return p->err;
}

static void parser_reinit (parser *p, double ***pZ, 
			   DATAINFO *dinfo, PRN *prn) 
{
    int saveflags = p->flags;

    p->flags = (P_START | P_PRIVATE | P_EXEC);

    if (saveflags & P_PRINT) {
	p->flags |= P_PRINT;
    }

    if (saveflags & P_NATEST) {
	p->flags |= P_NATEST;
    }

    if (saveflags & P_AUTOREG) {
	p->flags |= P_AUTOREG;
    }

    p->Z = pZ;
    p->dinfo = dinfo;
    p->prn = prn;

    p->obs = 0;
    p->sym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->getstr = 0;

    p->ret = NULL;
    p->err = 0;
    p->warn = 0;

    *p->warning = '\0';

    /* matrix: check the LH name again */
    if (p->targ == MAT && p->lh.m0 == NULL) {
	p->lh.m0 = get_matrix_by_name(p->lh.name);
    }
}

static void parser_init (parser *p, const char *str, 
			 double ***pZ, DATAINFO *dinfo,
			 PRN *prn, int flags)
{
    p->input = str;
    p->point = p->rhs = p->input;
    p->Z = pZ;
    p->dinfo = dinfo;
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
    p->lh.obs = 0;
    p->lh.m0 = NULL;
    p->lh.m1 = NULL;
    p->lh.substr = NULL;
    p->lh.mspec = NULL;
    p->lh.islist = 0;

    p->obs = 0;
    p->sym = 0;
    p->ch = 0;
    p->xval = 0.0;
    p->idnum = 0;
    p->idstr = NULL;
    p->getstr = 0;
    p->err = 0;
    p->warn = 0;
    p->ecount = 0;

    *p->warning = '\0';

    if (p->flags & P_SLICE) {
	p->lh.t = MAT;
    } else if (p->flags & P_SCALAR) {
	p->targ = NUM;
    } else if (p->flags & P_SERIES) {
	p->targ = VEC;
    } else if (p->flags & P_UFUN) {
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
	    } else {
		printnode(p->ret, p);
		pputc(p->prn, '\n');
	    }
	} else if (p->flags & (P_SCALAR | P_SERIES)) {
	    gen_check_return_type(p);
	} else if (p->flags & P_DECL) {
	    do_decl(p);
	} else if (p->Z != NULL) {
	    save_generated_var(p, prn);
	} 
    }

#if 0
    if (p->targ == MAT) {
	fprintf(stderr, "genr exec (%s): "
		" m0 = %p, m1 = %p\n", p->lh.name, (void *) p->lh.m0, 
		(void *) p->lh.m1);
    } else {
	if (p->lh.v == 0) {
	    fprintf(stderr, "genr exec (%s): p->lh.v = %d\n", p->lh.name,
		    p->lh.v);
	}
    }
#endif    
}

void gen_cleanup (parser *p)
{
    if (reusable(p)) {
#if PRESERVE_AUX_NODES
	if (p->ret != p->tree && !is_aux_node(p->ret)) {
	    free_tree(p->ret, "p->ret");
	    p->ret = NULL;
	}
#else
	if (p->ret != p->tree) {
	    free_tree(p->ret, "p->ret");
	    p->ret = NULL;
	}
#endif
    } else {
	if (p->ret != p->tree) {
	    free_tree(p->tree, "p->tree");
	}
	free_tree(p->ret, "p->ret");
	free(p->lh.substr);
	free(p->lh.mspec);
    }
}

static void maybe_set_return_flags (parser *p)
{
    NODE *t = p->tree;

    if (t != NULL && (t->t == SORT || t->t == DSORT)) {
	NODE *l = t->v.b1.b;

	if (l->t == USERIES) {
	    p->flags |= P_SORT;
	}
    } else if (t != NULL && t->t == UFUN) {
	p->flags |= P_UFRET;
    }
}

static int decl_check (parser *p, int flags)
{
    if (flags & P_COMPILE) {
	p->err = E_PARSE;
	sprintf(gretl_errmsg, "Bare declarations are not allowed here:\n> '%s'",
		p->input);
    }

    return p->err;
}

static void autoreg_error (parser *p, int t)
{
    fprintf(stderr, "*** autoreg error at obs t = %d (t1 = %d):\n", 
	    t, p->dinfo->t1);

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

int realgen (const char *s, parser *p, double ***pZ, 
	     DATAINFO *pdinfo, PRN *prn, int flags)
{
    int t;

    if (flags & P_EXEC) {
	parser_reinit(p, pZ, pdinfo, prn);
	if (p->flags & P_PRINT) {
	    p->ret = lhs_copy_node(p);
	    return p->err;
	} else {
	    goto starteval;
	}
    } else {
	parser_init(p, s, pZ, pdinfo, prn, flags);
	if (p->err) {
	    errmsg(p->err, prn);
	    return 1;
	}
    }

    if (p->flags & P_DECL) {
	decl_check(p, flags);
	return p->err;
    }

    if (p->op == INC || p->op == DEC || (p->flags & P_PRINT)) {
	p->ret = lhs_copy_node(p);
	return p->err;
    }

    lex(p);
    if (p->err) {
	fprintf(stderr, "exiting on lex() error (p->err = %d)\n",
		p->err);
	return p->err;
    }

    p->tree = expr(p);
    if (p->err) {
	fprintf(stderr, "exiting on expr() error\n");
	return p->err;
    }

#if EDEBUG
    fprintf(stderr, "after expr, p->tree->type = %d\n", p->tree->t);
#endif

    if (p->ch != 0) {
	parser_ungetc(p);
	context_error(p->ch, p);
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

    if (p->flags & P_AUTOREG) {
	/* e.g. y = b*y(-1) : evaluate dynamically */
	for (t=p->dinfo->t1; t<p->dinfo->t2 && !p->err; t++) {
	    const double *x;

	    p->aux_i = 0;
	    p->obs = t;
#if EDEBUG
	    fprintf(stderr, "\n*** autoreg: p->obs = %d\n", p->obs);
#endif
	    p->ret = eval(p->tree, p);
	    if (p->ret != NULL && p->ret->t == VEC) {
		x = getvec(p->ret, p);
		if (!na(x[t])) { 
#if EDEBUG
		    fprintf(stderr, "writing xvec[%d] = %g into Z[%d][%d]\n",
			    t, x[t], p->lh.v, t);
#endif
		    (*p->Z)[p->lh.v][t] = x[t];
		} 
	    } else {
		autoreg_error(p, t);
	    }
	    if (t == p->dinfo->t1) {
		p->flags &= ~P_START;
	    } 
	}
	p->obs = t;
    } 

    p->aux_i = 0;
    if (!p->err) {
	p->ret = eval(p->tree, p);
    }

#if EDEBUG > 1
    printnode(p->ret, p);
    pputc(prn, '\n');
#endif

#if PRESERVE_AUX_NODES
    if (!reusable(p)) {
	parser_free_aux_nodes(p);
    }
#else
    parser_free_aux_nodes(p);
#endif

    gen_check_errvals(p);

    if (flags & P_EXEC) {
	/* context is NLS/MLE or similar */
	if (p->warn != 0 && p->err == 0) {
	    /* warnings re. NAs become errors */
	    p->err = p->warn;
	}
    }

    return p->err;
}
