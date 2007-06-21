/*
 *  Copyright (c) 2004 by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include "system.h"
#include "var.h"
#include "objstack.h"
#include "matrix_extra.h"
#include "gretl_restrict.h"
#include "bootstrap.h"

#define RDEBUG 0

typedef struct restriction_ restriction;

struct restriction_ {
    int nterms;      /* number of terms in restriction */
    double *mult;    /* array of numerical multipliers on coeffs */
    int *eq;         /* array of equation numbers (for multi-equation case) */
    int *bnum;       /* array of coeff numbers */
    double rhs;      /* numerical value on right-hand side */
};

/* FIXME: "cross" below is a misnomer.  At least for VECM purposes we
   need to distinguish between multi-equation restrictions (which
   will be OK before too long) and true cross-equation (which we're
   not yet ready to handle).  What's designated as "cross" below is
   mostly just the multi-equation case.
*/

struct restriction_set_ {
    int k;                        /* number of restrictions (rows) */
    int cross;                    /* includes cross-equation restrictions? */
    gretl_matrix *R;              /* LHS restriction matrix */
    gretl_matrix *q;              /* RHS restriction matrix */
    char *mask;                   /* selection mask for coeffs */
    restriction **restrictions;
    void *obj;
    GretlObjType type;
    gretlopt opt;
    double test;
    double pval;
    double bsum;
    double bsd;
    int code;
};

static int check_R_matrix (const gretl_matrix *R)
{
    gretl_matrix *m;
    int k = gretl_matrix_rows(R);
    int err = 0;

    m = gretl_matrix_alloc(k, k);
    if (m == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      m, GRETL_MOD_NONE);

    err = gretl_invert_general_matrix(m);
    
    gretl_matrix_free(m);

    return err;
}

static int 
get_R_vecm_column (const gretl_restriction_set *rset, int i, int j)
{
    const restriction *r = rset->restrictions[i];
    GRETL_VAR *var = rset->obj;
    int nb = gretl_VECM_n_beta(var);
    int col = r->bnum[j];
    int k;

    if (r->eq != NULL) {
	for (k=0; k<r->eq[j]; k++) {
	    col += nb;
	}
    }

    return col;
}

static int 
get_R_sys_column (const gretl_restriction_set *rset, int i, int j)
{
    const restriction *r = rset->restrictions[i];
    gretl_equation_system *sys = rset->obj;
    const int *list;
    int col = r->bnum[j];
    int k;

    for (k=0; k<r->eq[j]; k++) {
	list = system_get_list(sys, k);
	col += list[0] - 1;
    }

    return col;
}

static int get_R_column (const gretl_restriction_set *rset, int i, int j)
{
    if (rset->type == GRETL_OBJ_SYS) {
	return get_R_sys_column(rset, i, j);
    } else if (rset->type == GRETL_OBJ_VAR) {
	return get_R_vecm_column(rset, i, j);
    } else {
	return 0;
    }
}

static double get_restriction_param (const restriction *r, int k)
{
    double x = 0.0;
    int i;    

    for (i=0; i<r->nterms; i++) {
	if (r->bnum[i] == k) {
	    x = r->mult[i];
	    break;
	}
    }

    return x;
}

/* FIXME this can be done more elegantly */

static int multiple_eq (gretl_restriction_set *rset)
{
    int i;

    for (i=0; i<rset->k; i++) { 
	if (rset->restrictions[i]->eq != NULL) {
	    return 1;
	}
    }

    return 0;
}

static int R_n_columns (gretl_restriction_set *rset)
{
    int i, cols = 0;

    if (rset->type == GRETL_OBJ_VAR) {
	cols = gretl_VECM_n_beta(rset->obj);
	if (multiple_eq(rset)) {
	    cols *= gretl_VECM_rank(rset->obj);
	}
    } else if (rset->type == GRETL_OBJ_SYS) {
	cols = system_n_indep_vars(rset->obj);
    } else {
	MODEL *pmod = rset->obj;
	
	if (rset->mask == NULL) {
	    cols = pmod->ncoeff;
	} else {
	    for (i=0; i<pmod->ncoeff; i++) {
		if (rset->mask[i]) {
		    cols++;
		}
	    }
	}
    }

    return cols;
}

/* used when generating restricted estimates (single-equation OLS
   only) */

static int 
restriction_set_form_full_matrices (gretl_restriction_set *rset)
{
    MODEL *pmod;
    gretl_matrix *R = NULL;
    gretl_vector *q = NULL;
    restriction *r;
    double x;
    int i, j, k;

    if (rset->type != GRETL_OBJ_EQN || rset->obj == NULL) {
	return 1;
    }

    pmod = rset->obj;
    k = pmod->ncoeff;

    R = gretl_matrix_alloc(rset->k, k);
    if (R == NULL) {
	return E_ALLOC;
    }

    q = gretl_column_vector_alloc(rset->k);
    if (q == NULL) {
	gretl_matrix_free(R);
	return E_ALLOC;
    }

    gretl_matrix_zero(R);
    gretl_matrix_zero(q);

    for (i=0; i<rset->k; i++) { 
	r = rset->restrictions[i];
	for (j=0; j<k; j++) {
	    if (rset->mask[j]) {
		x = get_restriction_param(r, j);
		gretl_matrix_set(R, i, j, x);
	    }
	}
	gretl_vector_set(q, i, r->rhs);
    }

    rset->R = R;
    rset->q = q;

    return 0;
}

/* create the matrices needed for testing a set of restrictions */

static int 
restriction_set_form_matrices (gretl_restriction_set *rset)
{
    MODEL *pmod = NULL;
    gretl_matrix *R = NULL;
    gretl_vector *q = NULL;
    restriction *r;
    double x;
    int col, i, j;
    int err = 0;

    R = gretl_zero_matrix_new(rset->k, R_n_columns(rset));
    if (R == NULL) {
	return E_ALLOC;
    }

    q = gretl_zero_matrix_new(rset->k, 1);
    if (q == NULL) {
	gretl_matrix_free(R);
	return E_ALLOC;
    }

    if (rset->type == GRETL_OBJ_EQN) {
	pmod = rset->obj;
    }

    for (i=0; i<rset->k; i++) { 
	r = rset->restrictions[i];
	if (pmod != NULL) {
	    col = 0;
	    for (j=0; j<pmod->ncoeff; j++) {
		if (rset->mask[j]) {
		    x = get_restriction_param(r, j);
		    gretl_matrix_set(R, i, col++, x);
		}
	    }
	} else {	    
	    for (j=0; j<r->nterms; j++) {
		col = get_R_column(rset, i, j);
		x = r->mult[j];
		gretl_matrix_set(R, i, col, x);
	    }
	} 
	gretl_vector_set(q, i, r->rhs);
    }

#if RDEBUG
    gretl_matrix_print(R, "R");
    gretl_matrix_print(q, "q");
#endif

    rset->R = R;
    rset->q = q;

    return err;
}

/* Make a mask with 1s in positions in the array of coeffs where
   a coeff is referenced in one or more restrictions, 0s otherwise.
   We do this only for single-equation restriction sets.
*/

static int restriction_set_make_mask (gretl_restriction_set *rset)
{
    MODEL *pmod;
    restriction *r;
    int i, j;

    if (rset->type != GRETL_OBJ_EQN || rset->obj == NULL) {
	return 1;
    }

    pmod = rset->obj;
    rset->mask = calloc(pmod->ncoeff, 1);

    if (rset->mask == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    for (i=0; i<rset->k; i++) {
	r = rset->restrictions[i];
	for (j=0; j<r->nterms; j++) {
	    rset->mask[r->bnum[j]] = 1;
	}	
    }

    return 0;
}

static int count_ops (const char *p)
{
    int n = 0 ;

    while (*p) {
	if (*p == '+' || *p == '-') n++;
	if (*p == '=') break;
	p++;
    }

    return n;
}

/* Given the dataset position of a variable as the identifier for a
   parameter, try to retrieve the corresponding 0-based coefficient
   number in the model to be restricted.  We don't yet attempt this
   for anything other than single-equation models.
*/

static int 
bnum_from_vnum (gretl_restriction_set *r, int v, const DATAINFO *pdinfo)
{
    const MODEL *pmod;
    int k;

    if (r->type != GRETL_OBJ_EQN || r->obj == NULL) {
	return -1;
    }

    pmod = r->obj;

    k = gretl_model_get_param_number(pmod, pdinfo, pdinfo->varname[v]);

#if RDEBUG
    fprintf(stderr, "bnum_from_vnum: vnum = %d (%s) -> bnum = %d (coeff = %g)\n", 
	    v, pdinfo->varname[v], k, pmod->coeff[k]);
#endif

    /* convert to 1-based for compatibility with numbers read directly:
       the index will be converted to 0-base below */

    return k + 1; 
}

/* Pick apart strings of the form "b[X]" or "b[X,Y]".  If the ",Y" is
   present the "X" element must be an equation number, and the "Y" may
   be a coefficient number or the name of a variable.  If the ",Y" is
   not present, "X" may be a coefficient number or the name of a
   variable.  This function is actually fed the string in question at
   an offset of 1 beyond the "[".  In parsing, we skip any white space
   in the string.
*/

static int pick_apart (gretl_restriction_set *r, const char *s, 
		       int *eq, int *bnum,
		       const DATAINFO *pdinfo)
{
    char s1[16] = {0};
    char s2[16] = {0};
    char *targ = s1;
    int i, j, k;
    int vnum = -1;

#if RDEBUG
    fprintf(stderr, "pick_apart: looking at '%s'\n", s);
#endif

    *eq = *bnum = -1;

    k = haschar(']', s);
    if (k <= 0 || k > 30) {
	return E_PARSE;
    }

    j = 0;
    for (i=0; i<k; i++) {
	if (s[i] == ',') {
	    targ = s2;
	    j = 0;
	} else if (!isspace(s[i])) {
	    if (j == 15) {
		return E_PARSE;
	    }
	    targ[j++] = s[i];
	}
    }

#if RDEBUG
    fprintf(stderr, " s1 = '%s', s2 = '%s'\n", s1, s2);
#endif

    if (targ == s2) {
	/* got a comma separator: [eqn,bnum] */
	*eq = positive_int_from_string(s1);
	if (*eq <= 0) {
	    return E_PARSE;
	}
	if (isdigit(*s2)) {
	    *bnum = positive_int_from_string(s2);
	} else if (pdinfo != NULL) {
	    vnum = varindex(pdinfo, s2);
	}
    } else {
	/* only one field: [bnum] */
	*eq = 1;
	if (isdigit(*s1)) {
	    *bnum = positive_int_from_string(s1);
	} else if (pdinfo != NULL) {
	    vnum = varindex(pdinfo, s1);
	}	    
    }

    if (pdinfo != NULL && vnum >= 0 && vnum < pdinfo->v) {
	*bnum = bnum_from_vnum(r, vnum, pdinfo);
    }

    return 0;
}

static int parse_b_bit (gretl_restriction_set *r, const char *s, 
			int *eq, int *bnum,
			const DATAINFO *pdinfo)
{
    int err = E_PARSE;

    if (isdigit((unsigned char) *s)) {
	sscanf(s, "%d", bnum);
	err = 0;
    } else if (*s == '[') {
	err = pick_apart(r, s + 1, eq, bnum, pdinfo);
    }

    if (*bnum < 1) {
	sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		*bnum);
	err = 1;
    } else {
	*bnum -= 1; /* convert to zero base */
    }

    if (*eq < 1) {
	sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		*eq);
	err = 1;
    } else {
	*eq -= 1; /* convert to zero base */
    }

    return err;
}

static int 
parse_coeff_chunk (gretl_restriction_set *r, const char *s, double *x, 
		   int *eq, int *bnum, const DATAINFO *pdinfo)
{
    const char *s0 = s;
    int err = E_PARSE;

    *eq = 1;

    while (isspace((unsigned char) *s)) s++;

    if (*s == 'b') {
	s++;
	err = parse_b_bit(r, s, eq, bnum, pdinfo);
	*x = 1.0;
    } else if (sscanf(s, "%lf", x)) {
	s += strspn(s, " ");
	s += strcspn(s, " *");
	s += strspn(s, " *");
	if (*s == 'b') {
	    s++;
	    err = parse_b_bit(r, s, eq, bnum, pdinfo);
	}
    }

    if (err && *gretl_errmsg == '\0') {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), s0);
    } 

#if RDEBUG
    fprintf(stderr, "parse_coeff_chunk: x=%g, eq=%d, bnum=%d\n", 
	    *x, *eq, *bnum);
#endif

    return err;
}

static void destroy_restriction (restriction *r)
{
    if (r == NULL) return;

    free(r->mult);
    free(r->eq);
    free(r->bnum);
    free(r);
}

void destroy_restriction_set (gretl_restriction_set *rset)
{
    int i;

    for (i=0; i<rset->k; i++) {
	destroy_restriction(rset->restrictions[i]);
    }

    free(rset->restrictions);
    free(rset->mask);
    
    gretl_matrix_free(rset->R);
    gretl_matrix_free(rset->q);

    free(rset);
}

static restriction *restriction_new (int n, int cross)
{
    restriction *r;
    int i;

    r = malloc(sizeof *r);
    if (r == NULL) {
	return NULL;
    }

    r->mult = NULL;
    r->eq = NULL;
    r->bnum = NULL;
	
    r->mult = malloc(n * sizeof *r->mult);
    r->bnum = malloc(n * sizeof *r->bnum);
    if (r->mult == NULL || r->bnum == NULL) {
	destroy_restriction(r);
	return NULL;
    }

    for (i=0; i<n; i++) {
	r->mult[i] = 0.0;
	r->bnum[i] = 0;
    }

    if (cross) {
	r->eq = malloc(n * sizeof *r->eq);
	if (r->eq == NULL) {
	    destroy_restriction(r);
	    return NULL;
	}
    }

    r->nterms = n;
    r->rhs = 0.0;

    return r;
}

static restriction *
augment_restriction_set (gretl_restriction_set *rset, int n_terms)
{
    restriction **rlist = NULL;
    int n = rset->k;

    rlist = realloc(rset->restrictions, (n + 1) * sizeof *rlist);
    if (rlist == NULL) {
	return NULL;
    }

    rset->restrictions = rlist;

    rset->restrictions[n] = restriction_new(n_terms, rset->cross);
    if (rset->restrictions[n] == NULL) {
	return NULL;
    }

    rset->k += 1;

    return rset->restrictions[n];
}

static void print_mult (double mult, int first, PRN *prn)
{
    if (mult == 1.0) {
	if (!first) pputs(prn, " + ");
    } else if (mult == -1.0) {
	if (first) pputs(prn, "-");
	else pputs(prn, " - ");
    } else if (mult > 0.0) {
	if (first) pprintf(prn, "%g*", mult);
	else pprintf(prn, " + %g*", mult);
    } else if (mult < 0.0) {
	if (first) pprintf(prn, "%g*", mult);
	else pprintf(prn, " - %g*", fabs(mult));
    }	
}

static void print_restriction (const gretl_restriction_set *rset,
			       int j, const DATAINFO *pdinfo, 
			       PRN *prn)
{
    const restriction *r = rset->restrictions[j];
    char vname[24];
    int i, k;

    for (i=0; i<r->nterms; i++) {
	k = r->bnum[i];
	print_mult(r->mult[i], i == 0, prn);
	if (rset->cross) {
	    pprintf(prn, "b[%d,%d]", r->eq[i] + 1, k + 1);
	} else if (rset->type == GRETL_OBJ_VAR) {
	    pprintf(prn, "b[%d]", k + 1);
	} else {
	    MODEL *pmod = rset->obj;

	    gretl_model_get_param_name(pmod, pdinfo, k, vname);
	    pprintf(prn, "b[%s]", vname);
	}
    }

    pprintf(prn, " = %g\n", r->rhs);
}

static void 
print_restriction_set (const gretl_restriction_set *rset, 
		       const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    if (rset->k > 1) {
	pputs(prn, _("Restriction set"));
    } else {
	pprintf(prn, "%s:", _("Restriction"));
    }
    pputc(prn, '\n');

    for (i=0; i<rset->k; i++) {
	if (rset->k > 1) {
	    pprintf(prn, " %d: ", i + 1);
	} else {
	    pputc(prn, ' ');
	}
	print_restriction(rset, i, pdinfo, prn);
    }
}

static int 
add_term_to_restriction (restriction *r, double mult, int eq, int bnum, int i)
{
    int j;

    for (j=0; j<i; j++) {
	if (bnum == r->bnum[j] && (r->eq == NULL || eq == r->eq[j])) {
	    /* additional reference to a previously referenced coeff */
	    r->mult[j] += mult;
	    r->nterms -= 1;
	    return 0;
	}
    }

    r->mult[i] = mult;
    r->bnum[i] = bnum;

    if (r->eq != NULL) {
	r->eq[i] = eq;
    }

    return 0;
}

static gretl_restriction_set *
real_restriction_set_start (void *ptr, GretlObjType type,
			    gretlopt opt)
{
    gretl_restriction_set *rset;

    rset = malloc(sizeof *rset);
    if (rset == NULL) return NULL;

    rset->obj = ptr;
    rset->type = type;
    rset->opt = opt;

    rset->test = NADBL;
    rset->pval = NADBL;
    rset->bsum = NADBL;
    rset->bsd = NADBL;

    rset->k = 0;
    rset->R = NULL;
    rset->q = NULL;
    rset->mask = NULL;
    rset->restrictions = NULL;
    rset->cross = 0;
    rset->code = GRETL_STAT_NONE;

    if (rset->type == GRETL_OBJ_SYS) {
	rset->cross = 1;
    } else if (rset->type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = ptr;

	if (var != NULL && gretl_VECM_rank(var) > 1) {
	    rset->cross = 1;
	}
    } 

    return rset;
}

/* check that the coefficients referenced in a restriction are
   within bounds, relative to the equation or system that is
   to be restricted */

static int bnum_out_of_bounds (const gretl_restriction_set *rset,
			       int bnum, int eq)
{
    int ret = 1;

    if (rset->type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = rset->obj;

	if (eq >= gretl_VECM_rank(var)) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= gretl_VECM_n_beta(var)) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		    bnum + 1);
	} else {
	    ret = 0;
	}
    } else if (rset->type == GRETL_OBJ_SYS) {
	gretl_equation_system *sys = rset->obj;
	const int *list = system_get_list(sys, eq);

	if (list == NULL) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= list[0] - 1) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) out of range "
				    "for equation %d"), bnum + 1, eq + 1);
	} else {
	    ret = 0;
	}
    } else {
	MODEL *pmod = rset->obj;

	if (eq > 0) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= pmod->ncoeff || bnum < 0) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		    bnum + 1);
	} else {
	    ret = 0;
	}
    }

    return ret;
}

static int 
real_restriction_set_parse_line (gretl_restriction_set *rset, 
				 const char *line,
				 const DATAINFO *pdinfo,
				 int first)
{
    const char *p = line;
    restriction *r;
    int sgn = 1;
    int i, nt, err = 0;

#if RDEBUG
    fprintf(stderr, "parse restriction line: got '%s'\n", line);
#endif

    if (!strncmp(p, "restrict", 8)) {
	if (strlen(line) == 8) {
	    if (first) {
		return 0;
	    } else {
		return E_PARSE;
	    }
	}
	p += 8;
	while (isspace((unsigned char) *p)) p++;
    }

    if (*p == '+' || *p == '-') {
	sgn = (*p == '+')? 1 : -1;
	p++;
    }

    nt = 1 + count_ops(p);

#if RDEBUG
    fprintf(stderr, "restriction line: assuming %d terms\n", nt);
#endif

    r = augment_restriction_set(rset, nt);

    if (r == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    for (i=0; i<nt; i++) {
	char chunk[32];
	int len, bnum = 1, eq = 1;
	double mult;

	len = strcspn(p, "+-=");
	if (len > 31) {
	    err = 1;
	    break;
	}

	*chunk = 0;
	strncat(chunk, p, len);
	p += len;

#if RDEBUG
	fprintf(stderr, " working on chunk %d, '%s'\n", i, chunk);
#endif

	err = parse_coeff_chunk(rset, chunk, &mult, &eq, &bnum, pdinfo);
	if (err) {
	    break;
	} else if (bnum_out_of_bounds(rset, bnum, eq)) {
	    err = E_DATA;
	    break;
	}

	mult *= sgn;
	add_term_to_restriction(r, mult, eq, bnum, i);

	if (*p == '+') {
	    sgn = 1.0;
	    p++;
	} else if (*p == '-') {
	    sgn = -1.0;
	    p++;
	}
    }

    if (!err) {
	if (!sscanf(p, " = %lf", &r->rhs)) {
	    err = E_PARSE;
	} else if (rset->type == GRETL_OBJ_VAR && r->rhs != 0.0) {
	    strcpy(gretl_errmsg, "VECM restrictions: the equations must "
		   "be homogeneous");
	    err = 1;
	}
    }

    if (err) {
	destroy_restriction_set(rset);
    } 
    
    return err;
}

int 
restriction_set_parse_line (gretl_restriction_set *rset, const char *line,
			    const DATAINFO *pdinfo)
{
    int nx = 0;

    if (rset->type == GRETL_OBJ_EQN) {
	MODEL *pmod = rset->obj;

	nx = pmod->ncoeff;
    } else if (rset->type == GRETL_OBJ_SYS) {
	nx = system_n_indep_vars(rset->obj);
    } else if (rset->type == GRETL_OBJ_VAR) {
	nx = gretl_VECM_n_beta(rset->obj);
    }

    if (rset->k >= nx) {
	sprintf(gretl_errmsg, _("Too many restrictions (maximum is %d)"), 
		nx - 1);
	destroy_restriction_set(rset);
	return 1;
    }

    return real_restriction_set_parse_line(rset, line, pdinfo, 0);
}

/* set-up for a set of restrictions for a VAR (vecm, actually) */

gretl_restriction_set *
var_restriction_set_start (const char *line, GRETL_VAR *var)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(var, GRETL_OBJ_VAR, OPT_NONE);
    if (rset == NULL) {
	strcpy(gretl_errmsg, _("Out of memory!"));
	return NULL;
    }

    gretl_error_clear();

    if (real_restriction_set_parse_line(rset, line, NULL, 1)) {
	if (*gretl_errmsg == '\0') {
	    sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	}
	return NULL;
    }

    return rset;
}

/* set-up for a set of cross-equation restrictions, for a system of
   simultaneous equations */

gretl_restriction_set *
cross_restriction_set_start (const char *line, gretl_equation_system *sys)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(sys, GRETL_OBJ_SYS, OPT_NONE);
    if (rset == NULL) {
	strcpy(gretl_errmsg, _("Out of memory!"));
	return NULL;
    }

    if (real_restriction_set_parse_line(rset, line, NULL, 1)) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	return NULL;
    }

    return rset;
}

/* set-up for a set of restrictions on a single equation */

gretl_restriction_set *
eqn_restriction_set_start (const char *line, MODEL *pmod, gretlopt opt)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(pmod, GRETL_OBJ_EQN, opt);
    if (rset == NULL) {
	strcpy(gretl_errmsg, _("Out of memory!"));
	return NULL;
    }

    if (real_restriction_set_parse_line(rset, line, NULL, 1)) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	return NULL;
    }

    return rset;
}

gretl_restriction_set *
restriction_set_start (const char *line, gretlopt opt, int *err)
{
    gretl_restriction_set *rset = NULL;
    char *name = NULL;
    GretlObjType type;
    void *ptr = NULL;

#if RDEBUG
    fprintf(stderr, "restriction_set_start: line='%s'\n", line);
#endif

    if (!strncmp(line, "restrict", 8)) {
	name = get_system_name_from_line(line);
    }

    if (name != NULL) {
	/* get pointer to named object */
	*err = gretl_get_object_and_type(name, &ptr, &type);
	if (ptr == NULL) {
	    sprintf(gretl_errmsg, "'%s': unrecognized name", name);
	}
    } else {
	/* get pointer to last-created object */
	ptr = get_last_model(&type);  
    }

    if (ptr == NULL) {
	*err = E_DATA;
	goto bailout;
    }

    if (type != GRETL_OBJ_EQN && type != GRETL_OBJ_SYS &&
	type != GRETL_OBJ_VAR) {
	*err = E_DATA;
	goto bailout;
    }

    if (type == GRETL_OBJ_VAR && gretl_is_restricted_VECM(ptr)) {
	strcpy(gretl_errmsg, _("The model is already restricted"));
	*err = E_DATA;
	goto bailout;
    }

    rset = real_restriction_set_start(ptr, type, opt);
    if (rset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err && name == NULL) {
	*err = real_restriction_set_parse_line(rset, line, NULL, 1);
	if (*err) {
	    rset = NULL;
	    if (*err == E_PARSE) {
		sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	    }
	}
    }

 bailout:

    free(name);

    return rset;
}

static int 
print_restricted_coeff (const MODEL *pmod, int i,
			double coeff, double sderr, int k,
			const DATAINFO *pdinfo, 
			PRN *prn)
{
    model_coeff mc;
    int gotnan = 0;

    model_coeff_init(&mc);

    if (xna(coeff)) {
	gotnan = 1;
    }

    mc.b = coeff;
    mc.se = sderr;

    if (!xna(coeff) && !xna(sderr) && sderr > 0) {
	mc.tval = coeff / sderr;
	mc.pval = coeff_pval(pmod->ci, mc.tval, pmod->dfd + k);
    }

    gretl_model_get_param_name(pmod, pdinfo, i, mc.name);

    print_coeff(&mc, prn);

    return gotnan;
}

static void coeff_header (const MODEL *pmod, PRN *prn)
{
    int use_param = pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM;

    print_coeff_heading(use_param, prn);
}

/* generate full restricted estimates: this function is used
   only for single-equation models, estimated via OLS */

static int 
do_restricted_estimates (gretl_restriction_set *rset,
			 const double **Z, const DATAINFO *pdinfo,
			 PRN *prn)
{
    MODEL *pmod = rset->obj;
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *S = NULL;
    int *xlist = NULL;
    double s2 = 0.0;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t;
    int yno, err = 0;

    X = gretl_matrix_alloc(T, k);
    y = gretl_matrix_alloc(T, 1);
    b = gretl_matrix_alloc(k, 1);
    S = gretl_matrix_alloc(k, k);

    if (X == NULL || y == NULL || b == NULL || S == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    err = restriction_set_form_full_matrices(rset);
    if (err) {
	goto bailout;
    }

    yno = gretl_model_get_depvar(pmod);
    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	gretl_vector_set(y, s, Z[yno][t]);
	for (i=0; i<k; i++) {
	    gretl_matrix_set(X, s, i, Z[xlist[i+1]][t]);
	}
	s++;
    }

#if RDEBUG
    gretl_matrix_print(rset->R, "R");
    gretl_matrix_print(rset->q, "q");
#endif

    err = gretl_matrix_restricted_ols(y, X, rset->R, rset->q, 
				      b, S, NULL, &s2);

    if (!err) {
	double v, coeff, se;

	pprintf(prn, "%s:\n\n", _("Restricted estimates"));
	coeff_header(pmod, prn);
	for (i=0; i<k; i++) {
	    coeff = gretl_vector_get(b, i);
	    v = gretl_matrix_get(S, i, i);
	    se = (v > 1.0e-16)? sqrt(v) : 0.0;
	    print_restricted_coeff(pmod, i, coeff, se, rset->k, pdinfo, prn);
	}
	pputc(prn, '\n');
	pprintf(prn, "  %s = %.*g\n", _("Standard error of residuals"), 
		GRETL_DIGITS, sqrt(s2));
    } 

 bailout:
    
    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(b);
    gretl_matrix_free(S);

    free(xlist);

    return err;
}

/* print result, single equation */

static void 
restriction_set_print_result (gretl_restriction_set *rset, 
			      const double **Z, const DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL *pmod = rset->obj;
    int robust, asym;

    robust = gretl_model_get_int(pmod, "robust");
    asym = ASYMPTOTIC_MODEL(pmod->ci);

    if (asym) {
	rset->code = GRETL_STAT_WALD_CHISQ;
	rset->pval = chisq_cdf_comp(rset->test, rset->k);
	pprintf(prn, "\n%s: %s(%d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust chi^2"): "chi^2",
		rset->k, rset->test);
    } else {
	rset->code = GRETL_STAT_F;
	rset->test /= rset->k;
	rset->pval = f_cdf_comp(rset->test, rset->k, pmod->dfd);
	pprintf(prn, "\n%s: %s(%d, %d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust F"): "F",
		rset->k, pmod->dfd, rset->test);
    }

    pprintf(prn, _("with p-value = %g\n"), rset->pval);
    pputc(prn, '\n');

    if (!(rset->opt & OPT_C)) {
	record_test_result(rset->test, rset->pval, _("restriction"));
    }

    if (pmod != NULL && Z != NULL && !(rset->opt & OPT_Q) 
	&& pmod->ci == OLS) {
	do_restricted_estimates(rset, Z, pdinfo, prn);
    }
}

/* execute the test, for a single equation */

static int 
test_restriction_set (gretl_restriction_set *rset, PRN *prn)
{
    MODEL *pmod = rset->obj;
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    gretl_vector *br = NULL;
    gretl_matrix *RvR = NULL;
    int err, freeRvR = 1;

    gretl_error_clear();

    err = restriction_set_form_matrices(rset);
    if (err) {
	return err;
    }

#if RDEBUG
    gretl_matrix_print(rset->R, "R matrix");
    gretl_matrix_print(rset->q, "q vector");
#endif

    if ((err = check_R_matrix(rset->R))) {
	if (err == E_SINGULAR) {
	    pputs(prn, _("Matrix inversion failed:\n"
			 " restrictions may be inconsistent or redundant\n"));
	} else {
	    err = E_ALLOC;
	}
	goto bailout;
    }

    b = gretl_coeff_vector_from_model(pmod, rset->mask);
    vcv = gretl_vcv_matrix_from_model(pmod, rset->mask);
    if (b == NULL || vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    br = gretl_column_vector_alloc(rset->k);
    if (br == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

#if RDEBUG
    gretl_matrix_print(vcv, "VCV matrix");
    gretl_matrix_print(b, "coeff vector");
#endif  

    err = gretl_matrix_multiply(rset->R, b, br);
    if (err) {
	fprintf(stderr, "Failed: gretl_matrix_multiply(R, b, br)\n");
	goto bailout;
    }

#if RDEBUG
    gretl_matrix_print(br, "br");
#endif  

    if (rset->opt & OPT_C) {
	rset->bsum = br->val[0];
    }

    if (!gretl_is_zero_matrix(rset->q)) {
	err = gretl_matrix_subtract_from(br, rset->q);
	if (err) {
	    fprintf(stderr, "Failed: gretl_matrix_subtract_from(br, q)\n");
	    goto bailout;
	}
    }

    if (gretl_is_identity_matrix(rset->R)) {
#if RDEBUG
	fprintf(stderr, "R is identity matrix: taking shortcut\n");
#endif  
	RvR = vcv;
	freeRvR = 0;
    } else {
	RvR = gretl_matrix_alloc(rset->R->rows, rset->R->rows);
	if (RvR == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
	gretl_matrix_qform(rset->R, GRETL_MOD_NONE, vcv,
			   RvR, GRETL_MOD_NONE);
#if RDEBUG
	gretl_matrix_print(RvR, "RvR");
#endif  
	if (rset->opt & OPT_C) {
	    rset->bsd = sqrt(RvR->val[0]);
	}
    }

    err = gretl_invert_symmetric_matrix(RvR);
    if (err) {
	pputs(prn, _("Matrix inversion failed:\n"
		     " restrictions may be inconsistent or redundant\n"));
	goto bailout;
    }
    
    rset->test = gretl_scalar_qform(br, RvR, &err);
    if (err) {
	pputs(prn, _("Failed to compute test statistic\n"));
	goto bailout;
    }

 bailout:

    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    gretl_vector_free(br);
    
    if (freeRvR) {
	gretl_matrix_free(RvR);
    }

    return err;
}

static int rset_expand_R (gretl_restriction_set *rset, int k)
{
    gretl_matrix *R = NULL;
    double rij;
    int i, j, jj;

    R = gretl_zero_matrix_new(rset->k, k);
    if (R == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<rset->k; i++) {
	jj = 0;
	for (j=0; j<k; j++) {
	    if (rset->mask[j]) {
		rij = gretl_matrix_get(rset->R, i, jj);
		gretl_matrix_set(R, i, j, rij);
		jj++;
	    }
	}
    }

    gretl_matrix_free(rset->R);
    rset->R = R;

    return 0;
}

static int do_single_equation_test (gretl_restriction_set *rset,
				    const double **Z,
				    const DATAINFO *pdinfo,
				    PRN *prn)
{
    int err, done = 0;

    if (rset->opt & OPT_B) {
	MODEL *pmod = rset->obj;

	if (!bootstrap_ok(pmod->ci)) {
	    pputs(prn, "Sorry, the bootstrap option is not supported for this test");
	} else {
	    restriction *r = rset->restrictions[0];

	    if (rset->k == 1 && r->nterms == 1 && r->rhs == 0) {
		/* a simple zero restriction */
		gretlopt bopt = OPT_P | OPT_R;
		int B = 0;

		gretl_restriction_get_boot_params(&B, &bopt);
		err = bootstrap_analysis(pmod, r->bnum[0], B, Z, pdinfo, 
					 bopt, prn);
		done = 1;
	    } else {
		/* a more complex restriction */
		err = test_restriction_set(rset, prn);
		if (!err) {
		    rset->test /= rset->k;
		    rset_expand_R(rset, pmod->ncoeff);
		    err = bootstrap_test_restriction(pmod, rset->R, rset->q,
						     rset->test, rset->k, Z, 
						     pdinfo, prn);
		}
	    }
	    done = 1;
	}
    }

    if (!done) {
	err = test_restriction_set(rset, prn);
	if (!err) {
	    restriction_set_print_result(rset, Z, pdinfo, prn);
	}
    }

    return err;
}

GRETL_VAR *
gretl_restricted_vecm (gretl_restriction_set *rset, 
		       double ***pZ,
		       DATAINFO *pdinfo,
		       PRN *prn,
		       int *err)
{
    gretl_matrix *D = NULL;
    GRETL_VAR *jvar = NULL;

    if (rset == NULL || rset->type != GRETL_OBJ_VAR) {
	*err = E_DATA;
	return NULL;
    }

    *err = restriction_set_form_matrices(rset);

    if (!*err) {
	D = gretl_matrix_right_nullspace(rset->R, err);
    }

    if (!*err) {
	print_restriction_set(rset, pdinfo, prn);
	jvar = real_gretl_restricted_vecm(rset->obj, rset->R, D, pZ, pdinfo, 
					  prn, err);
	rset->R = NULL; /* transfer R to vecm */
    }

    destroy_restriction_set(rset);
    gretl_matrix_free(D);

    return jvar;
}

/* Respond to "end restrict": in the case of a single equation, go
   ahead and do the test; in the case of a system of equations,
   form the restriction matrices R and q and attach these to the
   equation system.
*/

int
gretl_restriction_set_finalize (gretl_restriction_set *rset, 
				const double **Z,
				const DATAINFO *pdinfo,
				PRN *prn)
{
    int err = 0;

    if (rset == NULL) {
	return 1;
    }

    if (rset->type == GRETL_OBJ_VAR) {
	gretl_matrix *D = NULL;

	err = restriction_set_form_matrices(rset);
	if (!err) {
	    D = gretl_matrix_right_nullspace(rset->R, &err);
	}
	if (!err) {
	    print_restriction_set(rset, pdinfo, prn);
	    gretl_VECM_test_beta(rset->obj, D, pdinfo, prn);
	    gretl_matrix_free(D);
	}
	destroy_restriction_set(rset);
    } else if (rset->type == GRETL_OBJ_SYS) {
	/* simultaneous equations system */
#if RDEBUG
	print_restriction_set(rset, pdinfo, prn);
#endif
	err = restriction_set_form_matrices(rset);
	if (!err) {
	    err = check_R_matrix(rset->R);
	    if (err == E_SINGULAR) {
		pputs(prn, _("Matrix inversion failed:\n"
			     " restrictions may be inconsistent or redundant\n"));
	    } else if (err) {
		err = E_ALLOC;
	    }
	}	
	if (!err) {
	    system_set_restriction_matrices(rset->obj, rset->R, rset->q);
	    rset->R = NULL;
	    rset->q = NULL;
	    destroy_restriction_set(rset);
	}
    } else {
	/* single-equation model */
	err = restriction_set_make_mask(rset);
	if (!err) {
	    print_restriction_set(rset, pdinfo, prn);
	    err = do_single_equation_test(rset, Z, pdinfo, prn);
	    if (!(rset->opt & OPT_C)) {
		destroy_restriction_set(rset);
	    }
	}
    }	

    return err;
}

/**
 * gretl_sum_test:
 * @list: list of variables to use.
 * @pmod: pointer to model.
 * @pdinfo: information on the data set.
 * @prn: gretl printing struct.
 * 
 * Calculates the sum of the coefficients, relative to the given model, 
 * for the variables given in @list.  Prints this estimate along 
 * with its standard error.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int 
gretl_sum_test (const int *list, MODEL *pmod, DATAINFO *pdinfo,
		PRN *prn)
{
    gretl_restriction_set *r;
    char line[MAXLEN];
    char bstr[24];
    int i, len, err = 0;

    if (list[0] < 2) {
	pprintf(prn, _("Invalid input\n"));
	return E_DATA;
    }

    if (!command_ok_for_model(COEFFSUM, 0, pmod->ci)) {
	return E_NOTIMP;
    }

    r = real_restriction_set_start(pmod, GRETL_OBJ_EQN, OPT_Q | OPT_C);
    if (r == NULL) {
	return 1;
    }

    *line = '\0';
    len = 0;

    for (i=1; i<=list[0]; i++) {
	sprintf(bstr, "b[%s]", pdinfo->varname[list[i]]);
	len += strlen(bstr) + 4;
	if (len >= MAXLEN - 1) {
	    err = E_PARSE;
	    break;
	}
	strcat(line, bstr);
	if (i < list[0]) {
	    strcat(line, " + ");
	} else {
	    strcat(line, " = 0");
	}
    }

    if (!err) {
	err = real_restriction_set_parse_line(r, line, pdinfo, 1); 
    }

    if (!err) {
	err = gretl_restriction_set_finalize(r, NULL, pdinfo, NULL);
    }

    if (!err) {
	double test;

	pprintf(prn, "\n%s: ", _("Variables"));

	for (i=1; i<=list[0]; i++) {
	    pprintf(prn, "%s ", pdinfo->varname[list[i]]);
	}

	pprintf(prn, "\n   %s = %g\n", _("Sum of coefficients"), r->bsum);

	if (r->code == GRETL_STAT_F) {
	    pprintf(prn, "   %s = %g\n", _("Standard error"), r->bsd);
	    test = sqrt(r->test);
	    if (r->bsum < 0) {
		test = -test;
	    }
	    pprintf(prn, "   t(%d) = %g ", pmod->dfd, test);
	    pprintf(prn, _("with p-value = %g\n"), r->pval);
	    record_test_result(test, r->pval, _("sum")); 
	} else if (r->code == GRETL_STAT_WALD_CHISQ) {
	    pprintf(prn, "   %s = %g\n", _("Standard error"), r->bsd);
	    test = sqrt(r->test);
	    if (r->bsum < 0) {
		test = -test;
	    }
	    r->pval = normal_pvalue_2(test);
	    pprintf(prn, "   z = %g ", test);
	    pprintf(prn, _("with p-value = %g\n"), r->pval);
	    record_test_result(test, r->pval, _("sum")); 
	}	    

	destroy_restriction_set(r);
    }

    return err;
}

static int restrict_B;
static gretlopt rboot_opt;

int gretl_restriction_set_boot_params (int B, gretlopt opt)
{
    int err = 0;

    rboot_opt = opt;

    if (B > 0) {
	restrict_B = B;
    } else {
	err = E_DATA;
    }

    return err;
}

void gretl_restriction_get_boot_params (int *pB, gretlopt *popt)
{
    *pB = restrict_B;
    *popt |= rboot_opt;

    /* these are ad hoc values */
    restrict_B = 0;
    rboot_opt = OPT_NONE;
}

gretlopt gretl_restriction_get_options (const gretl_restriction_set *rset)
{
    return (rset != NULL)? rset->opt : OPT_NONE;
}
