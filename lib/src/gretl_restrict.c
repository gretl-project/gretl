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
    MODEL *pmod;
    gretl_equation_system *sys;
    GRETL_VAR *var;
    const DATAINFO *pdinfo;
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
    int nb = gretl_VECM_n_beta(rset->var);
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
    const int *list;
    int col = r->coeff[j];
    int k;

    for (k=0; k<r->eq[j]; k++) {
	list = system_get_list(rset->sys, k);
	col += list[0] - 1;
    }

    return col;
}

static int get_R_column (const gretl_restriction_set *rset, int i, int j)
{
    if (rset->sys != NULL) {
	return get_R_sys_column(rset, i, j);
    } else if (rset->var != NULL) {
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

    if (rset->var) {
	cols = gretl_VECM_n_beta(rset->var);
    } else if (rset->sys) {
	cols = system_n_indep_vars(rset->sys);
    } else {
	for (i=0; i<rset->pmod->ncoeff; i++) {
	    if (rset->mask[i]) {
		cols++;
	    }
	}
    }

    return cols;
}

static int 
restriction_set_form_matrices (gretl_restriction_set *rset,
			       gretl_matrix **Rin,
			       gretl_vector **qin)
{
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

    for (i=0; i<rset->k; i++) { 
	r = rset->restrictions[i];
	if (rset->pmod != NULL) {
	    col = 0;
	    for (j=0; j<rset->pmod->ncoeff; j++) {
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

    if (rset->var != NULL) {
	gretl_matrix *D = gretl_matrix_right_nullspace(R);

#if RDEBUG
	gretl_matrix_print(D, "D");
#endif
	gretl_VAR_attach_restrictions(rset->var, D);
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
    restriction *r;
    int i, j;

    if (rset->pmod == NULL) {
	return 1;
    }

    rset->mask = calloc(rset->pmod->ncoeff, 1);

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
	    err = 0;
	}	
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
    if (rlist == NULL) return NULL;

    rset->restrictions = rlist;

    rset->restrictions[n] = restriction_new(n_terms, rset->cross);
    if (rset->restrictions[n] == NULL) return NULL;

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
	    pprintf(prn, "b[%d,%d]", r->eq[i] + 1, r->coeff[i]);
	} else if (rset->var) {
	    const int *list = gretl_VECM_list(rset->var);
	    int li = (list == NULL)? 0 : r->coeff[i] + 1;

	    if (li > 0 && li <= list[0]) {
		pprintf(prn, "b[%s]", pdinfo->varname[list[li]]);
	    } else {
		pprintf(prn, "b[%d]", r->coeff[i]);
	    }
	} else {
	    gretl_model_get_param_name(rset->pmod, pdinfo, 
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
    r->mult[i] = mult;
    r->coeff[i] = bnum;

    if (r->eq != NULL) {
	r->eq[i] = eq;
    }

    return 0;
}

static gretl_restriction_set *
real_restriction_set_start (MODEL *pmod, const DATAINFO *pdinfo,
			    gretl_equation_system *sys,
			    GRETL_VAR *var)
{
    gretl_restriction_set *rset;

    rset = malloc(sizeof *rset);
    if (rset == NULL) return NULL;

    rset->pmod = pmod;
    rset->pdinfo = pdinfo;
    rset->sys = sys;
    rset->var = var;

    rset->k = 0;
    rset->mask = NULL;
    rset->restrictions = NULL;

    if (rset->sys != NULL) {
	rset->cross = 1;
    } else if (rset->var != NULL && gretl_VECM_rank(rset->var) > 1) {
	rset->cross = 1;
    } else {
	rset->cross = 0;
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

    if (rset->var) {
	if (eq >= gretl_VECM_rank(rset->var)) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= gretl_VECM_n_beta(rset->var)) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		    bnum);
	} else {
	    ret = 0;
	}
    } else if (rset->sys) {
	const int *list = system_get_list(rset->sys, eq);

	if (list == NULL) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= list[0] - 1) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) out of range "
				    "for equation %d"), bnum, eq + 1);
	} else {
	    ret = 0;
	}
    } else {
	if (eq > 0) {
	    sprintf(gretl_errmsg, _("Equation number (%d) is out of range"), 
		    eq + 1);
	} else if (bnum >= rset->pmod->ncoeff) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		    bnum);
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
	} else if (rset->var != NULL && r->rhs != 0.0) {
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

    if (rset->pmod != NULL) {
	nx = rset->pmod->ncoeff;
    } else if (rset->sys != NULL) {
	nx = system_n_indep_vars(rset->sys);
    } else if (rset->var != NULL) {
	nx = gretl_VECM_n_beta(rset->var);
    }

    if (rset->k >= nx) {
	sprintf(gretl_errmsg, _("Too many restrictions (maximum is %d)"), 
		nx - 1);
	destroy_restriction_set(rset);
	return 1;
    }

    return real_restriction_set_parse_line(rset, line, 0);
}

/* special set-up for a set of restrictions for a VAR (vecm,
   actually) */

gretl_restriction_set *
var_restriction_set_start (const char *line, GRETL_VAR *var)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(NULL, NULL, NULL, var);
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

/* special set-up for a set of cross-equation restrictions, for
   a system of simultaneous equations */

gretl_restriction_set *
cross_restriction_set_start (const char *line, gretl_equation_system *sys)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(NULL, NULL, sys, NULL);
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
restriction_set_start (const char *line, MODEL *pmod, const DATAINFO *pdinfo)
{
    gretl_restriction_set *rset = NULL;
    char *sname = NULL;

#if RDEBUG
    fprintf(stderr, "restriction_set_start: line='%s'\n", line);
#endif

    if (!strncmp(line, "restrict", 8)) {
	sname = get_system_name_from_line(line);
    }

    /* are we applying a restriction to a named system of equations? */
    if (sname != NULL) {
	gretl_equation_system *sys = NULL;
	GRETL_VAR *var = NULL;

	sys = get_equation_system_by_name(sname);
	if (sys != NULL) {
	    rset = real_restriction_set_start(NULL, NULL, sys, NULL);
	    if (rset == NULL) {
		strcpy(gretl_errmsg, _("Out of memory!"));
	    }	    
	} else {
	    var = get_VECM_by_name(sname);
	    if (var != NULL) {
		rset = real_restriction_set_start(NULL, NULL, NULL, var);
		if (rset == NULL) {
		    strcpy(gretl_errmsg, _("Out of memory!"));
		}
	    }
	}

	if (sys == NULL && var == NULL) {
	    sprintf(gretl_errmsg, "'%s': unrecognized name", sname);
	}

	free(sname);
    } else {
	rset = real_restriction_set_start(pmod, pdinfo, NULL, NULL);

	if (rset == NULL) {
	    strcpy(gretl_errmsg, _("Out of memory!"));
	} else if (real_restriction_set_parse_line(rset, line, 1)) {
	    sprintf(gretl_errmsg, _("parse error in '%s'\n"), line);
	    rset = NULL;
	}
    }

    return rset;
}

#if 0

static int 
estimate_restricted_model (const MODEL *pmod, const char *mask,
			   gretl_matrix *Rin, gretl_matrix *qin,
			   const double **Z, const DATAINFO *pdinfo,
			   PRN *prn)
{
    gretl_matrix *X = NULL;
    gretl_matrix *y = NULL;
    gretl_matrix *R = NULL;
    gretl_matrix *q = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *S = NULL;
    const int *list = pmod->list;
    double s2 = 0.0;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t;
    int err = 0;

    X = gretl_matrix_alloc(T, k);
    y = gretl_matrix_alloc(T, 1);
    b = gretl_matrix_alloc(k, 1);
    S = gretl_matrix_alloc(k, k);

    if (X == NULL || y == NULL || b == NULL || S == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	gretl_vector_set(y, s, Z[list[1]][t]);
	for (i=0; i<k; i++) {
	    gretl_matrix_set(X, s, i, Z[list[i+2]][t]);
	}
	s++;
    }

    gretl_matrix_print(R, "R");
    gretl_matrix_print(q, "q");
    

    err = gretl_matrix_restricted_ols(y, X, R, q, b, S, &s2);

    if (!err) {
	pputs(prn, "Restricted estimates:\n\n");
	for (i=0; i<k; i++) {
	    pprintf(prn, "b[%d] = %10.5g (%10.5g)\n", 
		    gretl_vector_get(b, i),
		    sqrt(gretl_matrix_get(S, i, i)));
	}
    } else {
	fprintf(stderr, "gretl_matrix_restricted_ols: err = %d\n", err);
    }

 bailout:
    
    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(b);
    gretl_matrix_free(S);

    return err;
}

#endif

/* execute the test, for a single equation */

static int test_restriction_set (gretl_restriction_set *rset, 
				 const double **Z,
				 const DATAINFO *pdinfo,
				 PRN *prn)
{
    gretl_matrix *R;
    gretl_vector *q;
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    gretl_vector *br = NULL;
    gretl_matrix *Rv = NULL;
    double test_stat, pval;
    int err, robust, freeRv = 1;

    int asym = ASYMPTOTIC_MODEL(rset->pmod->ci);

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

    b = gretl_coeff_vector_from_model(rset->pmod, rset->mask);
    vcv = gretl_vcv_matrix_from_model(rset->pmod, rset->mask);
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

    if (!gretl_is_zero_vector(q)) {
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

    robust = gretl_model_get_int(rset->pmod, "robust");

    if (asym) {
	pval = chisq(test_stat, rset->k);
	pprintf(prn, "\n%s: %s(%d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust chi^2"): "chi^2",
		rset->k, test_stat);
    } else {
	test_stat /= rset->k;
	pval = fdist(test_stat, rset->k, rset->pmod->dfd);
	pprintf(prn, "\n%s: %s(%d, %d) = %g, ", _("Test statistic"), 
		(robust)? _("Robust F"): "F",
		rset->k, rset->pmod->dfd, test_stat);
    }

    pprintf(prn, _("with p-value = %g\n"), pval);
    pputc(prn, '\n');

    record_test_result(test_stat, pval, _("restriction"));

#if 0
    estimate_restricted_model(rset->pmod, R, q, Z, pdinfo, prn);
#endif

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

    if (rset->var != NULL) {
	/* vecm */
	err = restriction_set_form_matrices(rset, NULL, NULL);
	if (!err) {
	    print_restriction_set(rset, pdinfo, prn);
	    gretl_VECM_test_beta(rset->var, prn);
	}
    } else if (rset->sys != NULL) {
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
	    system_set_restriction_matrices(rset->sys, R, q);
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

