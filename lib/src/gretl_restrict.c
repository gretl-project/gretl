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
#include "gretl_restrict.h"

#undef RDEBUG

typedef struct restriction_ restriction;

struct restriction_ {
    int nterms;
    double *mult;
    int *coeff;
    double rhs;
};

struct restriction_set_ {
    int n_equations;
    int ncoeff;
    char *select;
    restriction **restrictions;
    MODEL *pmod;
    const DATAINFO *pdinfo;
};


static void destroy_restriction_set (gretl_restriction_set *rset);

static void 
augment_selection_from_restriction (char *select, const restriction *r)
{
    int i;

    for (i=0; i<r->nterms; i++) {
	select[r->coeff[i]] = 1;
    }
}

static double get_restriction_param (const restriction *r, int k)
{
    int i;
    double x = 0.0;

    for (i=0; i<r->nterms; i++) {
	if (r->coeff[i] == k) {
	    x = r->mult[i];
	    break;
	}
    }

    return x;
}

static int check_R_matrix (const gretl_matrix *R)
{
    gretl_matrix *m;
    int k = gretl_matrix_rows(R);
    int err = 0;

    m = gretl_matrix_alloc(k, k);
    if (m == NULL) return GRETL_MATRIX_NOMEM;

    gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
			      R, GRETL_MOD_TRANSPOSE,
			      m);

    err = gretl_invert_general_matrix(m);
    
    gretl_matrix_free(m);

    return err;
}

static int 
restriction_set_form_matrices (gretl_restriction_set *rset,
			       gretl_matrix **Rin,
			       gretl_vector **rin)
{
    gretl_matrix *R;
    gretl_vector *r;
    restriction *rst;
    double x;
    int i, j, k;

    R = gretl_matrix_alloc(rset->n_equations, rset->ncoeff);
    if (R == NULL) return E_ALLOC;

    r = gretl_column_vector_alloc(rset->n_equations);
    if (r == NULL) {
	gretl_matrix_free(R);
	return E_ALLOC;
    }

    for (i=0; i<rset->n_equations; i++) { 
	rst = rset->restrictions[i];
	j = 0;
	for (k=0; k<rset->pmod->ncoeff; k++) {
	    if (rset->select[k]) {
		x = get_restriction_param(rst, k);
		gretl_matrix_set(R, i, j++, x);
	    }
	}
	gretl_vector_set(r, i, rst->rhs);
    }

    *Rin = R;
    *rin = r;

    return 0;
}

static int restriction_set_form_selection (gretl_restriction_set *rset)
{
    int i;

    rset->select = malloc(rset->pmod->ncoeff);
    if (rset->select == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    memset(rset->select, 0, rset->pmod->ncoeff);

    for (i=0; i<rset->n_equations; i++) {
	augment_selection_from_restriction(rset->select, 
					   rset->restrictions[i]);
    }

    rset->ncoeff = 0;
    for (i=0; i<rset->pmod->ncoeff; i++) {
	if (rset->select[i]) rset->ncoeff += 1;
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

static int parse_chunk (const char *s, double *x)
{
    int bnum = -1;

    while (isspace((unsigned char) *s)) s++;

    if (*s == 'b') {
	sscanf(s, " b%d", &bnum);
	*x = 1.0;
    } else {
	int got;

	got = sscanf(s, "%lf b%d", x, &bnum);
	if (got != 2) {
	    got = sscanf(s, "%lf * b%d", x, &bnum);
	}
    }

    if (bnum < 0) {
	sprintf(gretl_errmsg, _("parse error in '%s'\n"), s);
    }

    return bnum;
}

static void destroy_restriction (restriction *r)
{
    if (r == NULL) return;

    free(r->mult);
    free(r->coeff);
    free(r);
}

static void destroy_restriction_set (gretl_restriction_set *rset)
{
    int i;

    for (i=0; i<rset->n_equations; i++) {
	destroy_restriction(rset->restrictions[i]);
    }

    free(rset->restrictions);
    free(rset->select);
    free(rset);
}

static restriction *restriction_new (int n)
{
    restriction *r;

    r = malloc(sizeof *r);
    if (r == NULL) return NULL;

    r->mult = NULL;
    r->coeff = NULL;
	
    r->mult = malloc(n * sizeof *r->mult);
    r->coeff = malloc(n * sizeof *r->coeff);
    if (r->mult == NULL || r->coeff == NULL) {
	destroy_restriction(r);
	return NULL;
    }

    r->nterms = n;
    r->rhs = 0.0;

    return r;
}

static restriction *
augment_restriction_set (gretl_restriction_set *rset, int n_terms)
{
    int n = rset->n_equations;
    restriction **rlist = NULL;

    rlist = realloc(rset->restrictions, (n + 1) * sizeof *rlist);
    if (rlist == NULL) return NULL;

    rset->restrictions = rlist;

    rset->restrictions[n] = restriction_new(n_terms);
    if (rset->restrictions[n] == NULL) return NULL;

    rset->n_equations += 1;

    return rset->restrictions[n];
}

static const char *
get_varname (const gretl_restriction_set *rset, int cnum)
{
    int vnum;

    /* FIXME: garch -- vnum is wrong (probably ARMA too) 
       and TSLS?? */
    if (rset->pmod->ci == GARCH) {
	vnum = rset->pmod->list[cnum + 2];
    } else {
	vnum = rset->pmod->list[cnum + 2];
    }

    return rset->pdinfo->varname[vnum];
}

static void print_mult (double mult, int first, PRN *prn)
{
    if (mult == 1.0) {
	if (!first) pputs(prn, " + ");
    }
    else if (mult == -1.0) {
	if (first) pputs(prn, "-");
	else pputs(prn, " - ");
    }
    else if (mult > 0.0) {
	if (first) pprintf(prn, "%g*", mult);
	else pprintf(prn, " + %g*", mult);
    }
    else if (mult < 0.0) {
	if (first) pprintf(prn, "%g*", mult);
	else pprintf(prn, " - %g*", fabs(mult));
    }	
}

static void print_restriction (const gretl_restriction_set *rset,
			       int j, PRN *prn)
{
    const restriction *r = rset->restrictions[j];
    int i;

    /* FIXME: reinstate get_varname when it's corrected */

    for (i=0; i<r->nterms; i++) {
	print_mult(r->mult[i], i == 0, prn);
	if (rset->pmod->ci == ARMA ||
	    rset->pmod->ci == GARCH) {
	    pprintf(prn, "b%d", r->coeff[i]);
	} else {
	    pprintf(prn, "b(%s)", get_varname(rset, r->coeff[i]));
	}
    }
    pprintf(prn, " = %g\n", r->rhs);
}

static void 
print_restriction_set (const gretl_restriction_set *rset, PRN *prn)
{
    int i;

    if (rset->n_equations > 1) {
	pputs(prn, _("Restriction set"));
    } else {
	pprintf(prn, "%s:", _("Restriction"));
    }
    pputc(prn, '\n');

    for (i=0; i<rset->n_equations; i++) {
	if (rset->n_equations > 1) {
	    pprintf(prn, " %d: ", i + 1);
	} else {
	    pputs(prn, " ");
	}
	print_restriction(rset, i, prn);
    }

#ifdef RDEBUG
    pprintf(prn, "Selection mask for coefficients:");
    for (i=0; i<rset->pmod->ncoeff; i++) {
	pprintf(prn, " %d", rset->select[i]);
    }
    pputc(prn, '\n');

    pprintf(prn, "Number of coefficients referenced = %d\n",
	    rset->ncoeff);
#endif
}

static int 
add_term_to_restriction (restriction *r, double mult, int bnum, int i)
{
    r->mult[i] = mult;
    r->coeff[i] = bnum;

    return 0;
}

static gretl_restriction_set *
real_restriction_set_start (MODEL *pmod, const DATAINFO *pdinfo)
{
    gretl_restriction_set *rset;

    rset = malloc(sizeof *rset);
    if (rset == NULL) return NULL;

    rset->pmod = pmod;
    rset->pdinfo = pdinfo;

    rset->n_equations = 0;
    rset->ncoeff = 0;
    rset->select = NULL;
    rset->restrictions = NULL;
    
    return rset;
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

    if (!strncmp(p, "restrict", 8)) {
	if (strlen(line) == 8) {
	    if (first) return 0;
	    else return E_PARSE;
	}
	p += 8;
	while (isspace((unsigned char) *p)) p++;
    }

    if (*p == '+') p++;
    else if (*p == '-') {
	sgn = -1;
	p++;
    }

    nt = 1 + count_ops(p);

    r = augment_restriction_set(rset, nt);
    if (r == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    for (i=0; i<nt; i++) {
	char chunk[32];
	int len, bnum;
	double mult;

	len = strcspn(p, "+-=");
	if (len > 31) {
	    err = 1;
	    break;
	}
	*chunk = 0;
	strncat(chunk, p, len);
	p += len;

	bnum = parse_chunk(chunk, &mult);
	if (bnum < 0) {
	    err = E_PARSE;
	    break;
	} else if (bnum >= rset->pmod->ncoeff) {
	    sprintf(gretl_errmsg, _("Coefficient number (%d) is out of range"), 
		    bnum);
	    err = E_DATA;
	    break;
	}

	mult *= sgn;
	add_term_to_restriction(r, mult, bnum, i);

	if (*p == '+') {
	    sgn = 1.0;
	    p++;
	} else if (*p == '-') {
	    sgn = -1.0;
	    p++;
	}
    }

    if (!err) {
	if (!sscanf(p, " = %lf", &r->rhs)) err = E_PARSE;
    }

    if (err) {
	destroy_restriction_set(rset);
    } 
    
    return err;
}

int 
restriction_set_parse_line (gretl_restriction_set *rset, const char *line)
{
    if (rset->n_equations == rset->pmod->ncoeff) {
	sprintf(gretl_errmsg, _("Too many restrictions (maximum is %d)"), 
		rset->pmod->ncoeff);
	destroy_restriction_set(rset);
	return 1;
    }

    return real_restriction_set_parse_line(rset, line, 0);
}

gretl_restriction_set *
restriction_set_start (const char *line, MODEL *pmod, const DATAINFO *pdinfo)
{
    gretl_restriction_set *rset;

    rset = real_restriction_set_start(pmod, pdinfo);
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

static int test_restriction_set (gretl_restriction_set *rset, PRN *prn)
{
    gretl_matrix *R;
    gretl_vector *r;
    gretl_matrix *vcv = NULL;
    gretl_vector *b = NULL;
    gretl_vector *br = NULL;
    gretl_matrix *Rv = NULL;
    double F;
    int err, robust, freeRv = 1;

    *gretl_errmsg = '\0';

    err = restriction_set_form_matrices(rset, &R, &r);
    if (err) return err;

#ifdef RDEBUG
    gretl_matrix_print(R, "R matrix", prn);
    gretl_matrix_print(r, "r vector", prn);
#endif

    if ((err = check_R_matrix(R))) {
	if (err == GRETL_MATRIX_SINGULAR) {
	    pputs(prn, _("Matrix inversion failed:\n"
			 " restrictions may be inconsistent or redundant\n"));
	} else {
	    err = E_ALLOC;
	}
	goto bailout;
    }	

    b = gretl_coeff_vector_from_model(rset->pmod, rset->select);
    vcv = gretl_vcv_matrix_from_model(rset->pmod, rset->select);
    if (b == NULL || vcv == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    br = gretl_column_vector_alloc(rset->n_equations);
    if (br == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

#ifdef RDEBUG
    gretl_matrix_print(vcv, "VCV matrix", prn);
    gretl_matrix_print(b, "coeff vector", prn);
#endif  

    err = gretl_matrix_multiply(R, b, br);
    if (err) {
	fprintf(stderr, "Failed: gretl_matrix_multiply(R, b, br)\n");
	goto bailout;
    }

#ifdef RDEBUG
    gretl_matrix_print(br, "br", prn);
#endif  

    if (!gretl_is_zero_vector(r)) {
	err = gretl_matrix_subtract_from(br, r);
	if (err) {
	    fprintf(stderr, "Failed: gretl_matrix_subtract_from(R, b, br)\n");
	    goto bailout;
	}
    }

    if (gretl_is_identity_matrix(R)) {
#ifdef RDEBUG
	fprintf(stderr, "R is identity matrix: taking shortcut\n");
#endif  
	Rv = vcv;
	freeRv = 0;
    } else {
	Rv = gretl_matrix_A_X_A_prime(R, vcv, &err);
	if (err) goto bailout;
#ifdef RDEBUG
	gretl_matrix_print(Rv, "Rv", prn);
#endif  
    }

    err = gretl_invert_symmetric_matrix(Rv);
    if (err) {
	pputs(prn, _("Matrix inversion failed:\n"
		     " restrictions may be inconsistent or redundant\n"));
	goto bailout;
    }
    
    F = gretl_scalar_b_prime_X_b(br, Rv, &err);
    if (err) {
	pputs(prn, _("Failed to compute F statistic for test\n"));
	goto bailout;
    }

    F /= rset->n_equations;
    robust = gretl_model_get_int(rset->pmod, "robust");
    pprintf(prn, "\n%s: %s(%d, %d) = %g, ", _("Test statistic"), 
	    (robust)? _("Robust F"): "F",
	    rset->n_equations, rset->pmod->dfd, F);
    pprintf(prn, _("with p-value = %g\n"), 
	    fdist(F, rset->n_equations, rset->pmod->dfd));
    pputc(prn, '\n');

 bailout:

    gretl_matrix_free(R);
    gretl_vector_free(r);
    gretl_matrix_free(vcv);
    gretl_vector_free(b);
    gretl_vector_free(br);
    
    if (freeRv) gretl_matrix_free(Rv);

    return err;
}

int
gretl_restriction_set_finalize (gretl_restriction_set *rset, PRN *prn)
{
    int err;

    err = restriction_set_form_selection(rset);

    if (!err) {
	print_restriction_set(rset, prn);
	test_restriction_set(rset, prn);
	destroy_restriction_set(rset);
    }	

    return err;
}

