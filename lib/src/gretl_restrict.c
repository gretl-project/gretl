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
#include "gretl_restrict.h"

#define RDEBUG 0

typedef struct restriction_ restriction;

struct restriction_ {
    int nterms;      /* number of terms in restriction */
    double *mult;    /* array of numerical multipliers on coeffs */
    int *eq;         /* array of equation numbers (for cross-equation case) */
    int *coeff;      /* array of coeff numbers */
    double rhs;      /* numerical value on right-hand side */
};

struct restriction_set_ {
    int k;                        /* number of restrictions (rows) */
    int cross;                    /* includes cross-equation restrictions? */
    char *mask;                   /* selection mask for coeffs */
    restriction **restrictions;
    void *obj;
    GretlObjType type;
    gretlopt opt;
};

static void destroy_restriction_set (gretl_restriction_set *rset);

static int check_R_matrix (const gretl_matrix *R)
{
    gretl_matrix *m;
    int k = gretl_matrix_rows(R);
    int err = 0;

    m = gretl_matrix_alloc(k, k);
    if (m == NULL) return E_ALLOC;

    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      m);

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
    int col = r->coeff[j];
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
    int col = r->coeff[j];
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
	if (r->coeff[i] == k) {
	    x = r->mult[i];
	    break;
	}
    }

    return x;
}

static int R_n_columns (gretl_restriction_set *rset)
{
    int i, cols = 0;

    if (rset->type == GRETL_OBJ_VAR) {
	cols = gretl_VECM_n_beta(rset->obj);
    } else if (rset->type == GRETL_OBJ_SYS) {
	cols = system_n_indep_vars(rset->obj);
    } else {
	MODEL *pmod = rset->obj;

	for (i=0; i<pmod->ncoeff; i++) {
	    if (rset->mask[i]) {
		cols++;
	    }
	}
    }

    return cols;
}

static int 
restriction_set_form_full_matrices (gretl_restriction_set *rset,
				    gretl_matrix **Rin,
				    gretl_vector **qin)
{
    MODEL *pmod;
    gretl_matrix *R;
    gretl_vector *q;
    restriction *r;
    double x;
    int i, j, k;

    if (rset->type != GRETL_OBJ_EQN ||
	rset->obj == NULL) {
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

    *Rin = R;
    *qin = q;

    return 0;
}

static int 
restriction_set_form_matrices (gretl_restriction_set *rset,
			       gretl_matrix **Rin,
			       gretl_vector **qin)
{
    MODEL *pmod = NULL;
    gretl_matrix *R;
    gretl_vector *q;
    restriction *r;
    double x;
    int col, i, j;

    R = gretl_matrix_alloc(rset->k, R_n_columns(rset));

    if (R == NULL) return E_ALLOC;

    q = gretl_column_vector_alloc(rset->k);

    if (q == NULL) {
	gretl_matrix_free(R);
	return E_ALLOC;
    }

    gretl_matrix_zero(R);
    gretl_matrix_zero(q);

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

    if (rset->type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = rset->obj;

	gretl_matrix *D = gretl_matrix_right_nullspace(R);
#if RDEBUG
	gretl_matrix_print(D, "D");
#endif
	gretl_VAR_attach_restrictions(var, D);
	gretl_matrix_free(R);
	gretl_matrix_free(q);
    } else {
	*Rin = R;
	*qin = q;
    }

    return 0;
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

    if (rset->type != GRETL_OBJ_EQN ||
	rset->obj == NULL) {
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
	    rset->mask[r->coeff[j]] = 1;
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

static int parse_b_bit (const char *s, int *eq, int *bnum)
{
    int err = E_PARSE;

    if (isdigit((unsigned char) *s)) {
	sscanf(s, "%d", bnum);
	err = 0;
    } else if (*s == '[') {
	if (sscanf(s, "[%d,%d]", eq, bnum) == 2) {
	    err = 0;
	} else if (sscanf(s, "[%d]", bnum)) {
	    *eq = 0;
	    err = 0;
	}	
    }

    if (*bnum > 0) {
	*bnum -= 1;
    }

    if (*eq > 0) {
	*eq -= 1;
    }

    return err;
}

static int 
parse_coeff_chunk (const char *s, double *x, int *eq, int *bnum)
{
    const char *s0 = s;
    int err = E_PARSE;

    *eq = 0;

    while (isspace((unsigned char) *s)) s++;

    if (*s == 'b') {
	s++;
	err = parse_b_bit(s, eq, bnum);
	*x = 1.0;
    } else if (sscanf(s, "%lf", x)) {
	s += strspn(s, " ");
	s += strcspn(s, " *");
	s += strspn(s, " *");
	if (*s == 'b') {
	    s++;
	    err = parse_b_bit(s, eq, bnum);
	}
    }

    if (err) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), s0);
    } 

#if RDEBUG
    printf("parse_coeff_chunk: x=%g, eq=%d, bnum=%d\n", 
	   *x, *eq, *bnum);
#endif

    return err;
}

static void destroy_restriction (restriction *r)
{
    if (r == NULL) return;

    free(r->mult);
    free(r->eq);
    free(r->coeff);
    free(r);
}

static void destroy_restriction_set (gretl_restriction_set *rset)
{
    int i;

    for (i=0; i<rset->k; i++) {
	destroy_restriction(rset->restrictions[i]);
    }

    free(rset->restrictions);
    free(rset->mask);
    free(rset);
}

static restriction *restriction_new (int n, int cross)
{
    restriction *r;
    int i;

    r = malloc(sizeof *r);
    if (r == NULL) return NULL;

    r->mult = NULL;
    r->eq = NULL;
    r->coeff = NULL;
	
    r->mult = malloc(n * sizeof *r->mult);
    r->coeff = malloc(n * sizeof *r->coeff);
    if (r->mult == NULL || r->coeff == NULL) {
	destroy_restriction(r);
	return NULL;
    }

    for (i=0; i<n; i++) {
	r->mult[i] = 0.0;
	r->coeff[i] = 0;
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
    int i;

    for (i=0; i<r->nterms; i++) {
	print_mult(r->mult[i], i == 0, prn);
	if (rset->cross) {
	    pprintf(prn, "b[%d,%d]", r->eq[i] + 1, r->coeff[i] + 1);
	} else if (rset->type == GRETL_OBJ_VAR) {
	    GRETL_VAR *var = rset->obj;
	    const int *list = gretl_VECM_list(var);
	    int li = (list == NULL)? 0 : r->coeff[i];

	    if (li > 0 && li <= list[0]) {
		pprintf(prn, "b[%s]", pdinfo->varname[list[li]]);
	    } else {
		pprintf(prn, "b[%d]", r->coeff[i] + 1);
	    }
	} else {
	    MODEL *pmod = rset->obj;

	    gretl_model_get_param_name(pmod, pdinfo, 
				       r->coeff[i], vname);
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

#if RDEBUG
    if (rset->pmod != NULL) {
	pprintf(prn, "Selection mask for coefficients:");
	for (i=0; i<rset->pmod->ncoeff; i++) {
	    pprintf(prn, " %d", (int) rset->mask[i]);
	}
	pputc(prn, '\n');
    }
#endif
}

static int 
add_term_to_restriction (restriction *r, double mult, int eq, int bnum, int i)
{
    int j;

    for (j=0; j<i; j++) {
	if (bnum == r->coeff[j] && (r->eq == NULL || eq == r->eq[j])) {
	    /* additional reference to a previously referenced coeff */
	    r->mult[j] += mult;
	    r->nterms -= 1;
	    return 0;
	}
    }

    r->mult[i] = mult;
    r->coeff[i] = bnum;

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

    rset->k = 0;
    rset->mask = NULL;
    rset->restrictions = NULL;
    rset->cross = 0;

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
	} else if (bnum >= pmod->ncoeff) {
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
				 int first)
{
    const char *p = line;
    restriction *r;
    int sgn = 1;
    int i, nt, err = 0;

#if RDEBUG
    printf("parse restriction line: got '%s'\n", line);
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
    printf("restriction line: assuming %d terms\n", nt);
#endif

    r = augment_restriction_set(rset, nt);

    if (r == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    for (i=0; i<nt; i++) {
	char chunk[32];
	int len, bnum, eq;
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
	printf(" working on chunk %d, '%s'\n", i, chunk);
#endif

	err = parse_coeff_chunk(chunk, &mult, &eq, &bnum);
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
restriction_set_parse_line (gretl_restriction_set *rset, const char *line)
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

    return real_restriction_set_parse_line(rset, line, 0);
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

    *gretl_errmsg = '\0';

    if (real_restriction_set_parse_line(rset, line, 1)) {
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

    if (real_restriction_set_parse_line(rset, line, 1)) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	return NULL;
    }

    return rset;
}

/* set-up for a set of restrictions on a single equation */

gretl_restriction_set *
eqn_restriction_set_start (const char *line, MODEL *pmod)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(pmod, GRETL_OBJ_EQN, OPT_NONE);
    if (rset == NULL) {
	strcpy(gretl_errmsg, _("Out of memory!"));
	return NULL;
    }

    if (real_restriction_set_parse_line(rset, line, 1)) {
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
	return NULL;
    }

    if (type != GRETL_OBJ_EQN && type != GRETL_OBJ_SYS &&
	type != GRETL_OBJ_VAR) {
	*err = E_DATA;
	return NULL;
    }

    rset = real_restriction_set_start(ptr, type, opt);
    if (rset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err && name == NULL) {
	*err = real_restriction_set_parse_line(rset, line, 1);
	if (*err) {
	    rset = NULL;
	    if (*err == E_PARSE) {
		sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	    }
	}
    }	

    return rset;
}

static void print_pval_str (double pval, char *str)
{
    if (pval < .00001) {
	sprintf(str, "<%.5f", 0.00001);
    } else {
	sprintf(str, "%.5f", pval);
    }
}

static int print_coeff (const MODEL *pmod, int i,
			double coeff, double sderr, int k,
			const DATAINFO *pdinfo, 
			PRN *prn)
{
    int do_pval = 1;
    double t, pvalue = 999.0;
    int gotnan = 0;
    char varname[24];

    gretl_model_get_param_name(pmod, pdinfo, i, varname);
    pprintf(prn, "  %-15s ", varname);
    
    if (isnan(coeff) || na(coeff)) {
	pprintf(prn, "%*s", UTF_WIDTH(_("undefined"), 17), _("undefined"));
	gotnan = 1;
    } else {
	gretl_print_value(coeff, prn);
    }

    if (isnan(sderr) || na(sderr)) {
	pprintf(prn, "%*s\n", UTF_WIDTH(_("undefined"), 16), _("undefined"));
	return 1;
    }

    gretl_print_value(sderr, prn); 

    if (sderr > 0.0) {
	t = coeff / sderr;
	if (fabs(t) >= 1000.0) {
	    char numstr[9];

	    sprintf(numstr, "%#8.2G", t);
	    pprintf(prn, " %8s", numstr);
	} else {
	    pprintf(prn, " %7.3f", t);
	}

	if (do_pval) {
	    char pvalstr[16];
	    int dfd = pmod->dfd + k;

	    pvalue = coeff_pval(pmod, t, dfd);
	    print_pval_str(pvalue, pvalstr);
	    pprintf(prn, "%*s", UTF_WIDTH(pvalstr, 10), pvalstr);
	}
    } else if (do_pval) { 
	do_pval = 0;
	pprintf(prn, "     %*s", UTF_WIDTH(_("undefined"), 10), _("undefined"));
    }

    if (do_pval) {
	if (pvalue < 0.01) {
	    pputs(prn, " ***");
	} else if (pvalue < 0.05) {
	    pputs(prn, " **");
	} else if (pvalue < 0.10) {
	    pputs(prn, " *");
	}
    } 

    pputc(prn, '\n');

    return gotnan;
}

static void coeff_header (const MODEL *pmod, PRN *prn)
{
    int use_param = pmod->ci == NLS || pmod->ci == MLE;

    if (use_param) {
	pputs(prn, _("      PARAMETER       ESTIMATE          STDERROR"
		     "      T STAT   P-VALUE\n\n"));
    } else {
	pputs(prn, _("      VARIABLE       COEFFICIENT        STDERROR"
		     "      T STAT   P-VALUE\n\n"));
    }
}

static int 
do_restricted_estimates (gretl_restriction_set *rset,
			 const double **Z, const DATAINFO *pdinfo,
			 PRN *prn)
{
    MODEL *pmod = rset->obj;
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *q = NULL;
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

    err = restriction_set_form_full_matrices(rset, &R, &q);
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
    gretl_matrix_print(R, "R");
    gretl_matrix_print(q, "q");
#endif

    err = gretl_matrix_restricted_ols(y, X, R, q, b, S, &s2);

    if (!err) {
	double v, coeff, se;

	pprintf(prn, "%s:\n\n", _("Restricted estimates"));
	coeff_header(pmod, prn);
	for (i=0; i<k; i++) {
	    coeff = gretl_vector_get(b, i);
	    v = gretl_matrix_get(S, i, i);
	    se = (v > 1.0e-16)? sqrt(v) : 0.0;
	    print_coeff(pmod, i, coeff, se, rset->k, pdinfo, prn);
	}
	pputc(prn, '\n');
	pprintf(prn, "  %s = %.*g\n", _("Standard error of residuals"), 
		GRETL_DIGITS, sqrt(s2));
    } 

 bailout:
    
    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(R);
    gretl_matrix_free(q);
    gretl_matrix_free(b);
    gretl_matrix_free(S);

    free(xlist);

    return err;
}

/* execute the test, for a single equation */

static int test_restriction_set (gretl_restriction_set *rset, 
				 const double **Z,
				 const DATAINFO *pdinfo,
				 PRN *prn)
{
    MODEL *pmod = rset->obj;
    gretl_matrix *R;
    gretl_vector *q;
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    gretl_vector *br = NULL;
    gretl_matrix *Rv = NULL;
    double test_stat, pval;
    int err, robust, freeRv = 1;

    int asym = ASYMPTOTIC_MODEL(pmod->ci);

    *gretl_errmsg = '\0';

    err = restriction_set_form_matrices(rset, &R, &q);
    if (err) {
	return err;
    }

#if RDEBUG
    gretl_matrix_print(R, "R matrix");
    gretl_matrix_print(q, "q vector");
#endif

    if ((err = check_R_matrix(R))) {
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

    err = gretl_matrix_multiply(R, b, br);
    if (err) {
	fprintf(stderr, "Failed: gretl_matrix_multiply(R, b, br)\n");
	goto bailout;
    }

#if RDEBUG
    gretl_matrix_print(br, "br");
#endif  

    if (!gretl_is_zero_matrix(q)) {
	err = gretl_matrix_subtract_from(br, q);
	if (err) {
	    fprintf(stderr, "Failed: gretl_matrix_subtract_from(br, q)\n");
	    goto bailout;
	}
    }

    if (gretl_is_identity_matrix(R)) {
#if RDEBUG
	fprintf(stderr, "R is identity matrix: taking shortcut\n");
#endif  
	Rv = vcv;
	freeRv = 0;
    } else {
	Rv = gretl_matrix_A_X_A(R, GRETL_MOD_NONE, vcv, &err);
	if (err) goto bailout;
#if RDEBUG
	gretl_matrix_print(Rv, "Rv");
#endif  
    }

    err = gretl_invert_symmetric_matrix(Rv);
    if (err) {
	pputs(prn, _("Matrix inversion failed:\n"
		     " restrictions may be inconsistent or redundant\n"));
	goto bailout;
    }
    
    test_stat = gretl_scalar_b_X_b(br, GRETL_MOD_TRANSPOSE, Rv, &err);
    if (err) {
	pputs(prn, _("Failed to compute test statistic\n"));
	goto bailout;
    }

    robust = gretl_model_get_int(pmod, "robust");

    if (asym) {
	pval = chisq_cdf_comp(test_stat, rset->k);
	pprintf(prn, "\n%s: %s(%d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust chi^2"): "chi^2",
		rset->k, test_stat);
    } else {
	test_stat /= rset->k;
	pval = f_cdf_comp(test_stat, rset->k, pmod->dfd);
	pprintf(prn, "\n%s: %s(%d, %d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust F"): "F",
		rset->k, pmod->dfd, test_stat);
    }

    pprintf(prn, _("with p-value = %g\n"), pval);
    pputc(prn, '\n');

    record_test_result(test_stat, pval, _("restriction"));

    if (pmod != NULL && !(rset->opt & OPT_Q)
	&& pmod->ci != PANEL) {
	do_restricted_estimates(rset, Z, pdinfo, prn);
    }

 bailout:

    gretl_matrix_free(R);
    gretl_vector_free(q);
    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    gretl_vector_free(br);
    
    if (freeRv) {
	gretl_matrix_free(Rv);
    }

    return err;
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
	/* vecm */
	err = restriction_set_form_matrices(rset, NULL, NULL);
	if (!err) {
	    print_restriction_set(rset, pdinfo, prn);
	    gretl_VECM_test_beta(rset->obj, prn);
	}
    } else if (rset->type == GRETL_OBJ_SYS) {
	/* simultaneous equations system */
	gretl_matrix *R;
	gretl_matrix *q;

#if RDEBUG
	print_restriction_set(rset, pdinfo, prn);
#endif
	err = restriction_set_form_matrices(rset, &R, &q);
	if (!err) {
	    err = check_R_matrix(R);
	    if (err == E_SINGULAR) {
		pputs(prn, _("Matrix inversion failed:\n"
			     " restrictions may be inconsistent or redundant\n"));
	    } else if (err) {
		err = E_ALLOC;
	    }
	}	
	if (!err) {
	    system_set_restriction_matrices(rset->obj, R, q);
	    destroy_restriction_set(rset);
	}
    } else {
	/* single model */
	err = restriction_set_make_mask(rset);
	if (!err) {
	    print_restriction_set(rset, pdinfo, prn);
	    test_restriction_set(rset, Z, pdinfo, prn);
	    destroy_restriction_set(rset);
	}
    }	

    return err;
}

