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

#include "libgretl.h"
#include "system.h"
#include "var.h"
#include "objstack.h"
#include "uservar.h"
#include "matrix_extra.h"
#include "gretl_restrict.h"
#include "bootstrap.h"
#include "gretl_func.h"
#include "gretl_bfgs.h"
#include "cmd_private.h"

#define RDEBUG 0

#define R_UNSPEC -9

enum {
    VECM_NONE = 0,
    VECM_B = 1 << 0,
    VECM_A = 1 << 1
};

typedef struct rrow_ rrow;

struct rrow_ {
    int nterms;      /* number of terms in restriction */
    double *mult;    /* array of numerical multipliers on coeffs */
    int *eq;         /* array of equation numbers (for multi-equation case) */
    int *bnum;       /* array of coeff numbers */
    char *letter;    /* array of coeff letters */
    double rhs;      /* numerical value on right-hand side */
};

struct gretl_restriction_ {
    int g;                        /* number of restrictions (rows of R) */
    int gmax;                     /* max. possible restrictions */
    int bmulti;                   /* flag for coeffs need two indices */
    int amulti;                   /* VECM only */
    int gb, ga;                   /* restrictions on beta, alpha (VECM only) */
    int bcols, acols;             /* VECM only */
    int vecm;                     /* pertains to VECM beta or alpha? */
    gretl_matrix *R;              /* LHS restriction matrix */
    gretl_matrix *q;              /* RHS restriction matrix */
    gretl_matrix *Ra;             /* second LHS restriction matrix */
    gretl_matrix *qa;             /* second RHS restriction matrix */
    char *mask;                   /* selection mask for coeffs */
    rrow **rows;                  /* "atomic" restrictions */
    void *obj;                    /* pointer to model, system */
    GretlObjType otype;           /* type of model, system */
    gretlopt opt;                 /* OPT_C is used for "coeffsum" */
    char *rfunc;                  /* name of nonlinear restriction function */
    double test;                  /* test statistic */
    double pval;                  /* p-value of test statistic */
    double lnl;
    double bsum;
    double bse;
    int code;
};

#define targ_specified(r,i,j,c) (((c=='b' && r->bmulti) || (c=='a' && r->amulti)) \
				 && r->rows[i]->eq != NULL		\
				 && r->rows[i]->eq[j] != R_UNSPEC)

/* @R is putatively the matrix in R \beta -q = 0.  Check that it makes
   sense in that role, i.e., that RR' is non-singular. */

static int check_R_matrix (const gretl_matrix *R)
{
    gretl_matrix *m;
    int err = 0;

    m = gretl_matrix_alloc(R->rows, R->rows);

    if (m == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_multiply_mod(R, GRETL_MOD_NONE,
				  R, GRETL_MOD_TRANSPOSE,
				  m, GRETL_MOD_NONE);

	err = gretl_invert_general_matrix(m);

	if (err == E_SINGULAR) {
	    gretl_errmsg_set(_("Matrix inversion failed: restrictions may be "
			       "inconsistent or redundant"));
	}

	gretl_matrix_free(m);
    }

    return err;
}

/* Here we're adding to an existing restriction on a VECM, either in
   respect of beta (@letter = 'b') or in respect of alpha 
   (@letter = 'a').
*/

static int augment_vecm_restriction (gretl_restriction *rset,
				     char letter)
{
    GRETL_VAR *vecm = rset->obj;
    const gretl_matrix *R0, *q0;
    gretl_matrix *R2 = NULL, *q2 = NULL;
    int err = 0;

    if (letter == 'b') {
	R0 = gretl_VECM_R_matrix(vecm);
	q0 = gretl_VECM_q_matrix(vecm);
    } else {
	R0 = gretl_VECM_Ra_matrix(vecm);
	q0 = gretl_VECM_qa_matrix(vecm);
    }	

    R2 = gretl_matrix_row_concat(R0, rset->R, &err);
    if (err) {
	return err;
    }

    err = check_R_matrix(R2);
    if (err) {
	gretl_matrix_free(R2);
	return err;
    }

    if (q0 == NULL) {
	q2 = gretl_column_vector_alloc(R2->rows);
	if (q2 == NULL) {
	    err = E_ALLOC;
	} else {
	    int i, n = R0->rows;

	    for (i=0; i<R2->rows; i++) {
		q2->val[i] = (i < n)? 0.0 : rset->q->val[i-n];
	    }
	}
    } else {
	q2 = gretl_matrix_row_concat(q0, rset->q, &err);
    }

    if (err) {
	gretl_matrix_free(R2);
    } else if (letter == 'b') {
	gretl_matrix_replace(&rset->R, R2);
	gretl_matrix_replace(&rset->q, q2);
    } else {
	gretl_matrix_replace(&rset->Ra, R2);
	gretl_matrix_replace(&rset->qa, q2);
    }	

    return err;
}

/* Adding an existing restriction on a VECM, in respect of beta or
   alpha, to a new restriction in respect of the other matrix.
*/

static int add_old_vecm_restriction (gretl_restriction *rset,
				     char letter)
{
    GRETL_VAR *vecm = rset->obj; 
    int err = 0;

    if (letter == 'b') {
	const gretl_matrix *R = gretl_VECM_R_matrix(vecm);
	const gretl_matrix *q = gretl_VECM_q_matrix(vecm);

	rset->R = gretl_matrix_copy(R);
	if (rset->R == NULL) {
	    err = E_ALLOC;
	} else if (q != NULL) {	
	    rset->q = gretl_matrix_copy(q);
	    if (rset->q == NULL) {
		err = E_ALLOC;
	    }
	}
    } else {
	const gretl_matrix *Ra = gretl_VECM_Ra_matrix(vecm);
	const gretl_matrix *qa = gretl_VECM_qa_matrix(vecm);

	rset->Ra = gretl_matrix_copy(Ra);
	if (rset->Ra == NULL) {
	    err = E_ALLOC;
	} else if (qa != NULL) {	
	    rset->qa = gretl_matrix_copy(qa);
	    if (rset->qa == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    return err;
}

/* We set things up here such that the restrictions are
   R_b * vec(beta) = q; R_a * vec(alpha') = 0.
*/

static int get_R_vecm_column (const gretl_restriction *rset, 
			      int i, int j, char letter)
{
    const rrow *r = rset->rows[i];
    GRETL_VAR *var = rset->obj;
    int col = 0;

    if (letter == 'b') {
	col = r->bnum[j];
	if (rset->bmulti) {
	    col += r->eq[j] * gretl_VECM_n_beta(var);
	}
    } else if (letter == 'a') {
	col = r->bnum[j];
	if (rset->amulti) {
	    col += r->eq[j] * gretl_VECM_n_alpha(var);
	} 
    }

    return col;
}

static int get_R_sys_column (const gretl_restriction *rset, 
			     int i, int j)
{
    const rrow *r = rset->rows[i];
    equation_system *sys = rset->obj;
    int col = r->bnum[j];
    int k, n;

    for (k=0; k<r->eq[j]; k++) {
	n = system_get_list_length(sys, k);
	if (n > 0) {
	    col += n - 1;
	}
    }

    return col;
}

static double get_restriction_param (const rrow *r, int k)
{
    int i;    

    for (i=0; i<r->nterms; i++) {
	if (r->bnum[i] == k) {
	    return r->mult[i];
	}
    }

    return 0.0;
}

static void vecm_cross_error (void)
{
    gretl_errmsg_set("VECM: beta/alpha cross restrictions are "
		     "not handled");
}

/* See if we have both beta and alpha terms; check that there
   are no cross beta/alpha restrictions. */

static int vecm_ab_check (gretl_restriction *rset)
{
    int atotal = 0, btotal = 0;
    int i, j, err = 0;

    for (i=0; i<rset->g && !err; i++) {
	rrow *r = rset->rows[i];
	int a = 0, b = 0;

	for (j=0; j<r->nterms && !err; j++) {
	    if (r->letter[j] == 'a') {
		a++;
		atotal++;
	    } else {
		b++;
		btotal++;
	    }
	    if (a > 0 && b > 0) {
		vecm_cross_error();
		err = E_NOTIMP;
	    }
	}
    }

    if (!err) {
	/* incoming default is VECM_B */
	if (atotal > 0) {
	    rset->vecm |= VECM_A;
	}
	if (btotal == 0) {
	    rset->vecm &= ~VECM_B;
	}
    }

    return err;
}

/* check integrity of vecm restrictions, either beta or alpha */

static int vecm_x_check (gretl_restriction *rset, char letter)
{
    int *multi = (letter == 'a')? &rset->amulti : &rset->bmulti;
    int *ki = (letter == 'a')? &rset->ga : &rset->gb;
    int unspec = 0, anyspec = 0;
    int i, j, err = 0;

    for (i=0; i<rset->g && !err; i++) {
	rrow *r = rset->rows[i];

	if (r->letter[0] != letter) {
	    continue;
	}

	*ki += 1;

	for (j=0; j<r->nterms && !err; j++) {
	    if (r->eq != NULL) {
		if (r->eq[j] == R_UNSPEC) {
		    unspec = 1;
		} else if (r->eq[j] >= 0) {
		    anyspec = 1;
		} 
	    }
	    if (anyspec && unspec) {
		err = E_PARSE;
	    } 
	}
    }

    if (!err && unspec && *multi) {
	for (i=0; i<rset->g; i++) {
	    if (rset->rows[i]->letter[0] == letter) {
		free(rset->rows[i]->eq);
		rset->rows[i]->eq = NULL;
	    }
	}
	*multi = 0;
    }

    return err;
}

/* Check the validity of a set of VECM restrictions */

static int vecm_restriction_check (gretl_restriction *rset)
{
    int r = gretl_VECM_rank(rset->obj);
    int err;

    err = vecm_ab_check(rset);

    if (!err && (rset->vecm & VECM_B)) {
	err = vecm_x_check(rset, 'b');
    }

    if (!err && (rset->vecm & VECM_A)) {
	err = vecm_x_check(rset, 'a');
    }   

    if (!err && (rset->vecm & VECM_B)) {
	rset->bcols = gretl_VECM_n_beta(rset->obj);
	if (rset->bmulti) {
	    rset->bcols *= r;
	}
    }

    if (!err && (rset->vecm & VECM_A)) {
	rset->acols = gretl_VECM_n_alpha(rset->obj);
	if (rset->amulti) {
	    rset->acols *= r;
	}
    }    

    if (!err) {
	if (rset->gb > rset->bcols) {
	    gretl_errmsg_set("Too many restrictions");
	    err = E_NONCONF;
	} else if (rset->ga > rset->acols) {
	    fprintf(stderr, "rset->ga = %d, rset->acols = %d\n",
		    rset->ga, rset->acols);
	    gretl_errmsg_set("Too many restrictions");
	    err = E_NONCONF;
	}
    }

    return err;
}

/* If we have been given restrictions in row-wise fashion,
   allocate the required R and q matrices */

static int rset_allocate_R_q (gretl_restriction *rset, int nc)
{
    rset->R = gretl_zero_matrix_new(rset->g, nc);
    rset->q = gretl_zero_matrix_new(rset->g, 1);

#if RDEBUG
    fprintf(stderr, "rset_allocate_R_q: on output R=%p (%dx%d), q=%p (%dx%d)\n",
	    (void *) rset->R, rset->g, nc, (void *) rset->q, rset->g, 1);
#endif

    if (rset->R == NULL || rset->q == NULL) {
	gretl_matrix_free(rset->R);
	gretl_matrix_free(rset->q);
	rset->R = rset->q = NULL;
	return E_ALLOC;
    }

    return 0;
}

/* For the single-equation case, allocate and fill out the
   R and q matrices. */

static int equation_form_matrices (gretl_restriction *rset)
{
    MODEL *pmod = rset->obj;
    int nc = 0;
    int col, i, j;

    if (rset->mask == NULL) {
	nc = pmod->ncoeff;
    } else {
	for (i=0; i<pmod->ncoeff; i++) {
	    if (rset->mask[i]) {
		nc++;
	    }
	}
    }

    if (rset_allocate_R_q(rset, nc)) {
	return E_ALLOC;
    }

    for (i=0; i<rset->g; i++) { 
	rrow *r = rset->rows[i];
	double x;

	col = 0;
	for (j=0; j<pmod->ncoeff; j++) {
	    if (rset->mask[j]) {
		x = get_restriction_param(r, j);
		gretl_matrix_set(rset->R, i, col++, x);
	    }
	}
	gretl_vector_set(rset->q, i, r->rhs);
    }

    return 0;
}

/* For the equation-system case, allocate and fill out the
   R and q matrices. */

static int sys_form_matrices (gretl_restriction *rset)
{
    int nc = system_n_indep_vars(rset->obj);
    int col, i, j;
    int err = 0;

    if (rset_allocate_R_q(rset, nc)) {
	return E_ALLOC;
    }

    for (i=0; i<rset->g; i++) { 
	rrow *r = rset->rows[i];

	for (j=0; j<r->nterms; j++) {
	    col = get_R_sys_column(rset, i, j);
	    gretl_matrix_set(rset->R, i, col, r->mult[j]);
	}
	gretl_vector_set(rset->q, i, r->rhs);
    } 

    err = check_R_matrix(rset->R);

    return err;
}

/* For the VECM case, allocate and fill out the restriction
   matrices. */

static int vecm_form_matrices (gretl_restriction *rset)
{
    int i, j, m, err;

    err = vecm_restriction_check(rset);
    if (err) {
	return err;
    }

    if (rset->vecm & VECM_B) {
	rset->R = gretl_zero_matrix_new(rset->gb, rset->bcols);
	rset->q = gretl_zero_matrix_new(rset->gb, 1);
	if (rset->R == NULL || rset->q == NULL) {
	    return E_ALLOC;
	}
    }

    if (rset->vecm & VECM_A) {
	rset->Ra = gretl_zero_matrix_new(rset->ga, rset->acols);
	rset->qa = gretl_zero_matrix_new(rset->ga, 1);
	if (rset->Ra == NULL || rset->qa == NULL) {
	    return E_ALLOC;
	}
    }

    /* construct the restriction matrices, beta first then alpha 
       (if both are given) */

    for (m=0; m<2; m++) {
	gretl_matrix *R = (m == 0)? rset->R : rset->Ra;
	gretl_matrix *q = (m == 0)? rset->q : rset->qa;
	char letter = (m == 0)? 'b' : 'a';
	int col, row = 0;

	if ((letter == 'b' && rset->vecm == VECM_A) ||
	    (letter == 'a' && rset->vecm == VECM_B)) {
	    continue;
	}

	for (i=0; i<rset->g; i++) { 
	    rrow *r = rset->rows[i];

	    if (r->letter[0] == letter) {
		for (j=0; j<r->nterms; j++) {
		    col = get_R_vecm_column(rset, i, j, letter);
		    gretl_matrix_set(R, row, col, r->mult[j]);
		}
		gretl_vector_set(q, row, r->rhs);
		row++;
	    }
	}

#if RDEBUG
	if (letter == 'a') {
	    gretl_matrix_print(R, "Ra, constructed");
	    gretl_matrix_print(q, "qa, constructed");
	}
#endif
    }

    /* cumulate prior restrictions if needed */

    if (!err && beta_restricted_VECM(rset->obj)) {
	if (rset->vecm & VECM_B) {
	    err = augment_vecm_restriction(rset, 'b');
	} else {
	    err = add_old_vecm_restriction(rset, 'b');
	}
    }

    if (!err && alpha_restricted_VECM(rset->obj)) {
	if (rset->vecm & VECM_A) {
	    err = augment_vecm_restriction(rset, 'a');
	} else {
	    err = add_old_vecm_restriction(rset, 'a');
	}
    }    

    if (!err && rset->R != NULL) {
	err = check_R_matrix(rset->R);
    }

    if (!err && rset->Ra != NULL) {
	err = check_R_matrix(rset->Ra);
    }

#if RDEBUG
    fprintf(stderr, "vecm_form_matrices: returning %d\n", err);
#endif

    return err;
}

/* We were given R and/or q directly by the user: check them
   for sanity. */

static int test_user_matrices (gretl_restriction *rset)
{
    gretl_matrix *R = rset->R;
    gretl_matrix *q = rset->q;

    if (R == NULL || q == NULL) {
	/* we didn't get both parts */
	return E_DATA;
    }

    if (R->rows != q->rows || q->cols != 1) {
	/* R and q don't work */
	return E_NONCONF;
    }

    if (R->rows > rset->gmax) {
	/* too many restrictions */
	return E_NONCONF;
    }    

    if (rset->vecm) {
	rset->bcols = R->cols;
    }

    if (R->cols != rset->gmax) {
	if (rset->vecm) {
	    int nb = gretl_VECM_n_beta(rset->obj);

	    if (R->cols == nb && R->rows <= nb) {
		rset->bmulti = 0;
	    } else {
		return E_NONCONF;
	    }
	} else {
	    return E_NONCONF;
	}
    }

    return 0;
}

/* Set up the matrices needed for testing a set of restrictions:
   general driver function that calls more specific functions
   depending on the case. */

static int restriction_set_form_matrices (gretl_restriction *rset)
{
    int err = 0;

    if (rset->R != NULL || rset->q != NULL) {
	err = test_user_matrices(rset);
    } else if (rset->otype == GRETL_OBJ_EQN) {
	err = equation_form_matrices(rset);
    } else if (rset->otype == GRETL_OBJ_VAR) {
	err = vecm_form_matrices(rset);
    } else if (rset->g > 0) {
	/* mutivariate system */
	err = sys_form_matrices(rset);
    }

#if RDEBUG
    gretl_matrix_print(rset->R, "R");
    gretl_matrix_print(rset->q, "q");
#endif

    return err;
}

/* Make a mask with 1s in positions in the array of coeffs where
   a coeff is referenced in one or more restrictions, 0s otherwise.
   We do this only for single-equation restriction sets, and
   only when the restrictions are given in row-wise form.
*/

static int restriction_set_make_mask (gretl_restriction *rset)
{
    MODEL *pmod;
    rrow *r;
    int i, j;

    if (rset->otype != GRETL_OBJ_EQN || rset->obj == NULL) {
	return E_DATA;
    }

    pmod = rset->obj;
    rset->mask = calloc(pmod->ncoeff, 1);

    if (rset->mask == NULL) {
	destroy_restriction_set(rset);
	return E_ALLOC;
    }

    for (i=0; i<rset->g; i++) {
	r = rset->rows[i];
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

/* Try to retrieve the coefficient number in a model to be restricted,
   based on the parameter name.  We attempt this only for
   single-equation models.
*/

static int 
bnum_from_name (gretl_restriction *r, const DATASET *dset,
		const char *s)
{
    int k = -1;

    if (dset == NULL || r->obj == NULL) {
	gretl_errmsg_set(_("Please give a coefficient number"));
    } else if (r->otype == GRETL_OBJ_VAR) {
	GRETL_VAR *var = r->obj;
	int i, v;

	for (i=1; i<=var->ylist[0]; i++) {
	    v = var->ylist[i];
	    if (v > 0 && v < dset->v && !strcmp(dset->varname[v], s)) {
		k = i;
		break;
	    }
	}
	if (k < 0) {
	    gretl_errmsg_set(_("Please give a coefficient number"));
	}
    } else if (r->otype == GRETL_OBJ_EQN) {
	const MODEL *pmod = r->obj;

	k = gretl_model_get_param_number(pmod, dset, s);

	if (k < 0) {
	    gretl_errmsg_sprintf(_("%s: not a valid parameter name"), s);
	} else {
	    /* convert to 1-based for compatibility with numbers read
	       directly: the index will be converted to 0-base below
	    */
	    k++;
	}
    }

#if RDEBUG
    fprintf(stderr, "bnum_from_name: %s -> bnum = %d\n", s, k);
#endif

    return k; 
}

/* Pick apart strings of the form "b[X]" or "b[X,Y]".  If the ",Y" is
   present the "X" element must be an equation number, and the "Y" may
   be a coefficient number or the name of a variable.  If the ",Y" is
   not present, "X" may be a coefficient number or the name of a
   variable.  This function is actually fed the string in question at
   an offset of 1 beyond the "[".  We skip any white space in the 
   string.

   Such terms may be preceded by a multiplier, but they should not
   be followed by anything, so we flag an error if the closing "]"
   is not the end of the string.
*/

static int pick_apart (gretl_restriction *r, const char *s, 
		       int *eq, int *bnum,
		       const DATASET *dset)
{
    const char *junk;
    char s1[VNAMELEN] = {0};
    char s2[VNAMELEN] = {0};
    char *targ = s1;
    int i, j, k;
    int err = 0;

#if RDEBUG
    fprintf(stderr, "pick_apart: looking at '%s'\n", s);
#endif

    *eq = *bnum = -1;

    k = gretl_charpos(']', s);
    if (k <= 0 || k > 2*(VNAMELEN-1)) {
	return E_PARSE;
    }
    
    junk = s + k + 1;
    while (isspace(*junk)) junk++;
    if (*junk != '\0') {
	fprintf(stderr, "Got trailing junk '%s'\n", junk);
	return E_PARSE;
    }

    j = 0;
    for (i=0; i<k; i++) {
	if (s[i] == ',') {
	    targ = s2;
	    j = 0;
	} else if (!isspace(s[i])) {
	    if (j == VNAMELEN-1) {
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
	    err = E_PARSE;
	}
	if (isdigit(*s2)) {
	    *bnum = positive_int_from_string(s2);
	} else {
	    *bnum = bnum_from_name(r, dset, s2);
	}
    } else {
	/* only one field: [bnum] */
	*eq = R_UNSPEC;
	if (isdigit(*s1)) {
	    *bnum = positive_int_from_string(s1);
	} else {
	    *bnum = bnum_from_name(r, dset, s1);
	}	    
    }

    return err;
}

static int bbit_trailing_garbage (const char *s)
{
    if (*s) {
	while (isspace(*s)) s++;
	if (*s) {
	    if (*s == '/') {
		gretl_errmsg_sprintf(_("%s: division not allowed here"), s-1);
	    }
	    return 1;
	}
    }

    return 0;
}

/* We got a leading letter ('a' or 'b') and now we parse
   what follows.  The simplest case is just a number, giving
   e.g., "b1".  Then we assess the validity of what we've
   found.
*/

static int parse_b_bit (gretl_restriction *r, const char *s, 
			int *eq, int *bnum,
			const DATASET *dset)
{
    int err = E_PARSE;

    if (isdigit((unsigned char) *s)) {
	char *test;

	*bnum = strtol(s, &test, 10);
	if (bbit_trailing_garbage(test)) {
	    return err;
	}
	if (r->otype == GRETL_OBJ_VAR) {
	    *eq = R_UNSPEC;
	}
	err = 0;
    } else if (*s == '[') {
	err = pick_apart(r, s + 1, eq, bnum, dset);
    }

    if (*bnum < 1) {
	if (!gretl_errmsg_is_set() && *bnum >= 0) {
	    gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range"), 
				 *bnum);
	}
	err = 1;
    } else {
	*bnum -= 1; /* convert to zero base */
    }

    if (*eq == R_UNSPEC) {
	/* didn't get a specific equation number */
	if (r->otype == GRETL_OBJ_EQN) {
	    *eq = 0;
	} else if (r->otype != GRETL_OBJ_VAR) {
	    err = E_PARSE;
	}
    } else if (*eq == 0) {
	gretl_errmsg_sprintf(_("Equation number (%d) is out of range"), 
			     *eq);
	err = 1;
    } else {
	*eq -= 1; /* convert to zero base */
    }

    return err;
}

static int r_get_scalar (const char *s, double *x)
{
    if (sscanf(s, "%lf", x)) {
	return 1;
    } else {
	int n;

	s += strspn(s, " ");
	n = gretl_namechar_spn(s);

	if (n > 0) {
	    char tmp[VNAMELEN];

	    *tmp = '\0';
	    strncat(tmp, s, n);
	    if (gretl_is_scalar(tmp)) {
		*x = gretl_scalar_get_value(tmp, NULL);
		if (!na(*x)) {
		    return 1;
		}
	    } 
	}
    }

    return 0;
}

static int get_named_param (const char *s, MODEL *pmod,
			    int *bnum) 
{
    if (NONLIST_MODEL(pmod->ci) && pmod->params != NULL) {
	char test[VNAMELEN];
	int i;

	*test = '\0';
	strncat(test, s, VNAMELEN - 1);
	tailstrip(test);

	for (i=0; i<pmod->nparams; i++) {
	    if (!strcmp(test, pmod->params[i])) {
		*bnum = i;
		return 1;
	    }
	}
    }

    return 0;
}

#define b_start(r,s) ((*s == 'b' || (r->vecm && *s == 'a')) && \
		      (isdigit(*(s+1)) || *(s+1) == '['))

#define could_be_param(r,s) (r->otype == GRETL_OBJ_EQN && \
	                     !isdigit(*s))

static int 
parse_coeff_chunk (gretl_restriction *r, const char *s, double *x, 
		   int *eq, int *bnum, char *letter, 
		   const DATASET *dset)
{
    const char *s0 = s;
    int done = 0;
    int err = E_PARSE;

    *eq = 1;

    while (isspace((unsigned char) *s)) s++;

#if RDEBUG
    fprintf(stderr, "parse_coeff_chunk: s='%s'\n", s);
#endif

    if (b_start(r, s)) {
	/* something like "b[" or "b2" */
	*letter = *s;
	s++;
#if RDEBUG
	fprintf(stderr, " first branch: s='%s'\n", s);
#endif
	err = parse_b_bit(r, s, eq, bnum, dset);
	*x = 1.0;
	done = 1;
    } else if (could_be_param(r, s)) {
	/* e.g. "beta" for MLE or GMM */
	done = get_named_param(s, r->obj, bnum);
	if (done) {
	    *letter = 'b';
	    *x = 1.0;
	    *eq = 0;
	    err = 0;
	}
    }

    if (!done && r_get_scalar(s, x)) {
	s += strspn(s, " ");
	s += strcspn(s, " *");
	s += strspn(s, " *");
	if (b_start(r, s)) {
	    *letter = *s;
	    s++;
	    err = parse_b_bit(r, s, eq, bnum, dset);
	} else if (*s == '\0') {
	    /* plain numeric value, saved in *x */
	    *bnum = 0;
	    *letter = '\0';
	    err = 0;
	}
    }

    if (err) {
	gretl_errmsg_sprintf(_("parse error in '%s'\n"), s0);
    } 

#if RDEBUG
    fprintf(stderr, "parse_coeff_chunk: x=%g, eq=%d, bnum=%d, letter=%c\n", 
	    *x, *eq, *bnum, *letter);
#endif

    return err;
}

static void destroy_restriction (rrow *r)
{
    if (r != NULL) {
	free(r->mult);
	free(r->eq);
	free(r->bnum);
	free(r->letter);
	free(r);
    }
}

void destroy_restriction_set (gretl_restriction *rset)
{
    if (rset->rows != NULL) {
	int i;

	for (i=0; i<rset->g; i++) {
	    destroy_restriction(rset->rows[i]);
	}
	free(rset->rows);
    }

    free(rset->mask);

#if RDEBUG
    fprintf(stderr, "destroy_restriction_set: R at %p, q at %p\n",
	    (void *) rset->R, (void *) rset->q);
    fprintf(stderr, " Ra at %p, qa at %p\n",
	    (void *) rset->Ra, (void *) rset->qa);
#endif
    
    gretl_matrix_free(rset->R);
    gretl_matrix_free(rset->q);
    gretl_matrix_free(rset->Ra);
    gretl_matrix_free(rset->qa);

    free(rset->rfunc);

    free(rset);
}

static rrow *restriction_new (int n, int multi, int vecm)
{
    rrow *r;
    int i;

    r = malloc(sizeof *r);
    if (r == NULL) {
	return NULL;
    }

    r->mult = NULL;
    r->eq = NULL;
    r->bnum = NULL;
    r->letter = NULL;
	
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

    if (multi) {
	r->eq = malloc(n * sizeof *r->eq);
	if (r->eq == NULL) {
	    destroy_restriction(r);
	    return NULL;
	}
    }

    if (vecm) {
	r->letter = malloc(n * sizeof *r->letter);
	if (r->letter == NULL) {
	    destroy_restriction(r);
	    return NULL;
	}	
    }

    r->nterms = n;
    r->rhs = 0.0;

    return r;
}

static rrow *restriction_set_add_row (gretl_restriction *rset, 
				      int n_terms)
{
    rrow **rows = NULL;
    int n = rset->g;

    rows = realloc(rset->rows, (n + 1) * sizeof *rows);
    if (rows == NULL) {
	return NULL;
    }

    rset->rows = rows;

    rset->rows[n] = restriction_new(n_terms, rset->bmulti, rset->vecm);
    if (rset->rows[n] == NULL) {
	return NULL;
    }

    rset->g += 1;

    return rset->rows[n];
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

static void print_restriction (const gretl_restriction *rset,
			       int i, const DATASET *dset, 
			       PRN *prn)
{
    const rrow *r = rset->rows[i];
    char letter;
    char vname[24];
    int j, k;

    for (j=0; j<r->nterms; j++) {
	letter = (r->letter != NULL)? r->letter[j] : 'b';
	k = r->bnum[j];
	print_mult(r->mult[j], j == 0, prn);
	if (targ_specified(rset, i, j, letter)) {
	    pprintf(prn, "%c[%d,%d]", letter, r->eq[j] + 1, k + 1);
	} else if (rset->otype == GRETL_OBJ_VAR) {
	    pprintf(prn, "%c[%d]", letter, k + 1);
	} else {
	    MODEL *pmod = rset->obj;

	    gretl_model_get_param_name(pmod, dset, k, vname);
	    if (NONLIST_MODEL(pmod->ci)) {
		pputs(prn, vname);
	    } else {
		pprintf(prn, "b[%s]", vname);
	    }
	}
    }

    pprintf(prn, " = %g\n", r->rhs);
}

/* Pretty-print the set of restrictions.  We do this only if
   the restrictions were given row-wise by the user.
*/

static void print_restriction_set (const gretl_restriction *rset, 
				   const DATASET *dset, 
				   gretlopt opt, PRN *prn)
{
    int i;

    if (rset->g > 1) {
	pputs(prn, _("Restriction set"));
    } else {
	pprintf(prn, "%s:", _("Restriction"));
    }
    pputc(prn, '\n');

    for (i=0; i<rset->g; i++) {
	if (rset->g > 1) {
	    pprintf(prn, " %d: ", i + 1);
	} else {
	    pputc(prn, ' ');
	}
	print_restriction(rset, i, dset, prn);
    }

    if (opt & OPT_V) {
	if (rset->R != NULL || rset->q != NULL ||
	    rset->Ra != NULL || rset->qa != NULL) {
	    pputc(prn, '\n');
	}
	if (rset->R != NULL) {
	    gretl_matrix_print_to_prn(rset->R, "R", prn);
	}
	if (rset->q != NULL) {
	    gretl_matrix_print_to_prn(rset->q, "q", prn);
	}
	if (rset->Ra != NULL) {
	    gretl_matrix_print_to_prn(rset->Ra, "R(alpha)", prn);
	}
	if (rset->qa != NULL) {
	    gretl_matrix_print_to_prn(rset->qa, "q(alpha)", prn);
	}
    }	
}

void print_restriction_from_matrices (const gretl_matrix *R,
				      const gretl_matrix *q,
				      char letter, int npar, 
				      PRN *prn)
{
    double x;
    int eqn, coeff, started;
    int i, j;

    for (i=0; i<R->rows; i++) {
	started = 0;
	coeff = 1;
	eqn = (R->cols > npar)? 1 : 0;
	for (j=0; j<R->cols; j++) {
	    x = gretl_matrix_get(R, i, j);
	    if (x != 0.0) {
		if (!started) {
		    pputs(prn, "  ");
		}
		if (x == 1.0) {
		    if (started) {
			pputs(prn, " + ");
		    }
		} else if (x == -1.0) {
		    if (started) {
			pputs(prn, " - ");
		    } else {
			pputc(prn, '-');
		    }
		} else if (x > 0.0) {
		    if (started) {
			pprintf(prn, " + %g*", x);
		    } else {
			pprintf(prn, "%g*", x);
		    }
		} else if (x < 0.0) {
		    if (started) {
			pprintf(prn, " - %g*", -x);
		    } else {
			pprintf(prn, "%g*", x);
		    }
		}
		if (eqn > 0) {
		    pprintf(prn, "%c[%d,%d]", letter, eqn, coeff);
		} else {
		    pprintf(prn, "%c%d", letter, coeff);
		}
		started = 1;
	    }
	    if ((j + 1) % npar == 0) {
		eqn++;
		coeff = 1;
	    } else {
		coeff++;
	    }
	}
	pprintf(prn, " = %g\n", (q == NULL)? 0.0 : q->val[i]);
    }
}

static int 
add_term_to_restriction (rrow *r, double mult, int eq, int bnum, 
			 char letter, int i)
{
    int j;

    for (j=0; j<i; j++) {
	if (r->letter != NULL && r->letter[j] != letter) {
	    continue;
	}
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

    if (r->letter != NULL) {
	r->letter[i] = letter;
    }

    return 0;
}

static gretl_restriction *restriction_set_new (void *ptr, 
					       GretlObjType type,
					       gretlopt opt)
{
    gretl_restriction *rset;

    rset = malloc(sizeof *rset);
    if (rset == NULL) return NULL;

    rset->obj = ptr;
    rset->otype = type;
    rset->opt = opt;

    rset->test = NADBL;
    rset->pval = NADBL;
    rset->lnl = NADBL;
    rset->bsum = NADBL;
    rset->bse = NADBL;

    rset->g = rset->gmax = 0;
    rset->gb = rset->ga = 0;
    rset->R = NULL;
    rset->q = NULL;
    rset->Ra = NULL;
    rset->qa = NULL;
    rset->mask = NULL;
    rset->rows = NULL;

    rset->rfunc = NULL;

    rset->bmulti = 0;
    rset->amulti = 0;
    rset->bcols = 0;
    rset->acols = 0;
    rset->vecm = 0;
    rset->code = GRETL_STAT_NONE;

    if (rset->otype == GRETL_OBJ_EQN) {
	MODEL *pmod = ptr;

	rset->gmax = pmod->ncoeff;
    } else if (rset->otype == GRETL_OBJ_SYS) {
	rset->gmax = system_n_indep_vars(ptr);
	rset->bmulti = 1;
    } else if (rset->otype == GRETL_OBJ_VAR) {
	GRETL_VAR *var = ptr;

	if (var != NULL && gretl_VECM_rank(var) > 1) {
	    rset->bmulti = 1;
	    rset->amulti = 1;
	}
	rset->vecm = VECM_B;
	rset->gmax = gretl_VECM_n_beta(var) *
	    gretl_VECM_rank(var);
	rset->gmax += gretl_VECM_n_alpha(var) *
	    gretl_VECM_rank(var);
    } 

    return rset;
}

gretl_restriction *rset_from_VECM (GRETL_VAR *var, int *err)
{
    const gretl_matrix *R = gretl_VECM_R_matrix(var);
    const gretl_matrix *q = gretl_VECM_q_matrix(var);
    const gretl_matrix *Ra = gretl_VECM_Ra_matrix(var);
    const gretl_matrix *qa = gretl_VECM_qa_matrix(var);
    gretl_restriction *rset;

    if (R == NULL && q == NULL && Ra == NULL && qa == NULL) {
	/* no restriction in place */
	return NULL;
    }

    rset = calloc(1, sizeof *rset);
    if (rset == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    rset->obj = var;
    rset->otype = GRETL_OBJ_VAR;
    rset->opt = OPT_B;

    rset->R = (gretl_matrix *) R;
    rset->q = (gretl_matrix *) q;
    rset->Ra = (gretl_matrix *) Ra;
    rset->qa = (gretl_matrix *) qa;

    return rset;
}

/* check that the coefficients referenced in a restriction are
   within bounds, relative to the equation or system that is
   to be restricted */

static int bnum_out_of_bounds (const gretl_restriction *rset,
			       int i, int j, char letter)
{
    int ret = 1;

    if (rset->otype == GRETL_OBJ_VAR) {
	GRETL_VAR *var = rset->obj;

	if (i >= gretl_VECM_rank(var)) {
	    gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range"), 
				 i + 1);
	} else if ((letter == 'b' && j >= gretl_VECM_n_beta(var)) ||
		   (letter == 'a' && j >= gretl_VECM_n_alpha(var))) {
	    gretl_errmsg_sprintf(_("Equation number (%d) is out of range"), 
				 j + 1);
	} else {
	    ret = 0;
	}
    } else if (rset->otype == GRETL_OBJ_SYS) {
	equation_system *sys = rset->obj;
	const int *list = system_get_list(sys, i);

	if (list == NULL) {
	    gretl_errmsg_sprintf(_("Equation number (%d) is out of range"), 
				 i + 1);
	} else if (j >= list[0] - 1) {
	    gretl_errmsg_sprintf(_("Coefficient number (%d) out of range "
				   "for equation %d"), j + 1, i + 1);
	} else {
	    ret = 0;
	}
    } else {
	MODEL *pmod = rset->obj;

	if (i > 0) {
	    gretl_errmsg_sprintf(_("Equation number (%d) is out of range"), 
		    i + 1);
	} else if (j >= pmod->ncoeff || j < 0) {
	    if (!gretl_errmsg_is_set()) {
		gretl_errmsg_sprintf(_("Coefficient number (%d) is out of range"), 
				     j + 1);
	    }
	} else {
	    ret = 0;
	}
    }

    return ret;
}

/* Parse a matrix specification for a restriction set.  This should
   take the form "R = <matrix>" or "q = <matrix>", where <matrix> is
   the name of a pre-defined matrix.  On this approach we want exactly
   one R definition and exactly one q definition, and no other sort
   of specification can be given.
*/

static int read_matrix_line (const char *s, gretl_restriction *rset)
{
    const gretl_matrix *m;
    char mname[VNAMELEN];
    char job = *s;
    int err = 0;

    if (rset->rows != NULL) {
	/* we already have "row-wise" restrictions */
	return E_PARSE;
    } else if (job == 'R' && rset->R != NULL) {
	/* duplicate entry for R matrix */
	return E_PARSE;
    } else if (job == 'q' && rset->q != NULL) {
	/* duplicate entry for q matrix */
	return E_PARSE;
    }

    s++;
    while (isspace((unsigned char) *s)) s++;
    if (*s != '=') {
	return E_PARSE;
    }

    s++;
    while (isspace((unsigned char) *s)) s++;
    if (gretl_scan_varname(s, mname) != 1) {
	return E_PARSE;
    }

    m = get_matrix_by_name(mname);
    if (m == NULL) {
	return E_UNKVAR;
    } 

    if (job == 'R') {
	rset->R = gretl_matrix_copy(m);
	if (rset->R == NULL) {
	    err = E_ALLOC;
	}
    } else if (job == 'q') {
	rset->q = gretl_matrix_copy(m);
	if (rset->q == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

/* Try to read the name of a function to be called for testing
   a nonlinear restriction, as in "rfunc = <fname>".  This is
   admissable only if it's not mixed with any other sort
   of restriction.
*/

static int read_fncall_line (const char *s, gretl_restriction *rset)
{
    char fname[FN_NAMELEN];

    if (rset->rows != NULL || rset->R != NULL || rset->q != NULL) {
	/* other stuff is in the way! */
	return E_PARSE;
    }

    s += strspn(s, " ");

    if (*s != '=') {
	return E_PARSE;
    }

    s++;
    s += strspn(s, " ");

    if (sscanf(s, "%31s", fname) != 1) {
	return E_PARSE;
    }

    s += strlen(fname);

    if (!string_is_blank(s)) {
	/* there should be no trailing junk */
	return E_PARSE;
    }

    rset->rfunc = gretl_strdup(fname);

    return 0;
}

static int parse_restriction_row (gretl_restriction *rset,
				  const char *s,
				  const DATASET *dset)
{ 
    rrow *row;
    int nt, sgn = 1;
    int i, j;
    int err = 0;

    if (rset->R != NULL || rset->q != NULL) {
	/* can't mix row-wise spec with matrix spec */
	return E_PARSE;
    }

    if (*s == '+' || *s == '-') {
	sgn = (*s == '+')? 1 : -1;
	s++;
    }

    nt = 1 + count_ops(s);

#if RDEBUG
    fprintf(stderr, "restriction line: assuming %d terms\n", nt);
#endif

    row = restriction_set_add_row(rset, nt);

    if (row == NULL) {
	return E_ALLOC;
    }

    j = 0;

    for (i=0; i<nt; i++) {
	char chunk[32];
	int len, bnum = 1, eq = 1;
	int numeric = 0;
	char letter = 'b';
	double mult;

	len = strcspn(s, "+-=");
	if (len > 31) {
	    err = 1;
	    break;
	}

	*chunk = 0;
	strncat(chunk, s, len);
	s += len;

#if RDEBUG
	fprintf(stderr, " working on chunk %d, '%s'\n", i, chunk);
#endif

	err = parse_coeff_chunk(rset, chunk, &mult, &eq, &bnum, 
				&letter, dset);
	if (err) {
	    break;
	} else if (letter == '\0') {
	    numeric = 1;
	} else if (bnum_out_of_bounds(rset, eq, bnum, letter)) {
	    err = E_DATA;
	    break;
	}

	mult *= sgn;

	if (numeric) {
	    /* got a numeric constant */
	    row->rhs -= mult;
	    row->nterms -= 1;
	} else {
	    add_term_to_restriction(row, mult, eq, bnum, letter, j++);
	}

	if (*s == '+') {
	    sgn = 1.0;
	    s++;
	} else if (*s == '-') {
	    sgn = -1.0;
	    s++;
	}
    }

    if (!err) {
	char rhstr[32];
	double rhs = 0.0;

	s += strspn(s, " ");

	if (*s != '=') {
	    err = E_PARSE;
	} else {
	    s++;
	    if (!sscanf(s, "%31s", rhstr)) {
		err = E_PARSE;
	    } else if (sscanf(s, "%lf", &rhs) == 1) {
		row->rhs += rhs;
	    } else if (gretl_is_scalar(rhstr)) {
		rhs = gretl_scalar_get_value(rhstr, NULL);
		if (na(rhs)) {
		    err = E_DATA;
		} else {
		    row->rhs += rhs;
		}
	    }
	}
    }

#if RDEBUG
    fprintf(stderr, "restriction row: nterms = %d, rhs = %g\n", 
	    row->nterms, row->rhs);
#endif

    if (!err && row->nterms == 0) {
	gretl_errmsg_set(_("Invalid restriction"));
	err = E_PARSE;
    }

    return err;
}

static int 
real_restriction_set_parse_line (gretl_restriction *rset, 
				 const char *line,
				 const DATASET *dset,
				 int first)
{
    const char *s = line;
    int err = 0;

    gretl_error_clear();

#if RDEBUG
    fprintf(stderr, "parse restriction line: got '%s'\n", line);
#endif

    /* if we have a function name specified already, nothing else
       should be supplied */
    if (rset->rfunc != NULL) {
	destroy_restriction_set(rset);
	return E_PARSE;
    }

    if (!strncmp(s, "restrict", 8)) {
	if (strlen(line) == 8) {
	    if (first) {
		return 0;
	    } else {
		return E_PARSE;
	    }
	}
	s += 8;
	while (isspace((unsigned char) *s)) s++;
    }

    if (!strncmp(s, "rfunc", 5)) {
	/* nonlinear function given */
	err = read_fncall_line(s + 5, rset);
    } else if (*s == 'R' || *s == 'q') {
	/* restrictions given in matrix form */
	err = read_matrix_line(s, rset);
    } else {
	/* a regular linear restriction */
	err = parse_restriction_row(rset, s, dset);
    }

    if (err) {
	destroy_restriction_set(rset);
    }    
    
    return err;
}

int 
restriction_set_parse_line (gretl_restriction *rset, const char *line,
			    const DATASET *dset)
{
    if (rset->g > rset->gmax) {
	gretl_errmsg_sprintf(_("Too many restrictions (maximum is %d)"), 
			     rset->gmax);
	destroy_restriction_set(rset);
	return E_DATA;
    }

    return real_restriction_set_parse_line(rset, line, dset, 0);
}

/* set-up for a set of restrictions for a VAR (vecm, actually) */

gretl_restriction *
var_restriction_set_start (const char *line, GRETL_VAR *var)
{
    gretl_restriction *rset;

    rset = restriction_set_new(var, GRETL_OBJ_VAR, OPT_NONE);
    if (rset == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    gretl_error_clear();

    if (real_restriction_set_parse_line(rset, line, NULL, 1)) {
	if (!gretl_errmsg_is_set()) {
	    gretl_errmsg_sprintf(_("parse error in '%s'\n"), line);
	}
	return NULL;
    }

    return rset;
}

/* set-up for a set of (possibly) cross-equation restrictions, for a
   system of simultaneous equations */

gretl_restriction *
cross_restriction_set_start (const char *line, equation_system *sys)
{
    gretl_restriction *rset;

    rset = restriction_set_new(sys, GRETL_OBJ_SYS, OPT_NONE);
    if (rset == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    if (real_restriction_set_parse_line(rset, line, NULL, 1)) {
	gretl_errmsg_sprintf(_("parse error in '%s'\n"), line);
	return NULL;
    }

    return rset;
}

/* set-up for a set of restrictions on a single equation */

gretl_restriction *
eqn_restriction_set_start (const char *line, MODEL *pmod, 
			   const DATASET *dset,
			   gretlopt opt)
{
    gretl_restriction *rset;

    rset = restriction_set_new(pmod, GRETL_OBJ_EQN, opt);
    if (rset == NULL) {
	gretl_errmsg_set(_("Out of memory!"));
	return NULL;
    }

    if (real_restriction_set_parse_line(rset, line, dset, 1)) {
	gretl_errmsg_sprintf(_("parse error in '%s'\n"), line);
	return NULL;
    }

    return rset;
}

gretl_restriction *
restriction_set_start (const char *line, gretlopt opt, int *err)
{
    gretl_restriction *rset = NULL;
    char *name = NULL;
    GretlObjType type;
    void *ptr = NULL;

#if RDEBUG
    fprintf(stderr, "restriction_set_start: line='%s'\n", line);
#endif

    if (!strncmp(line, "restrict ", 9)) {
	name = get_system_name_from_line(line + 9);
    }

    if (name != NULL) {
	/* get pointer to named object */
	*err = gretl_get_object_and_type(name, &ptr, &type);
	if (ptr == NULL) {
	    gretl_errmsg_sprintf("'%s': unrecognized name", name);
	}
    } else {
	/* get pointer to last-estimated model */
	ptr = get_last_model(&type);
    }

    if (ptr == NULL) {
	*err = E_DATA;
	goto bailout;
    }

#if RDEBUG
    fprintf(stderr, " restriction: ptr = %p, type = %d\n", ptr, type);
#endif

    if (type != GRETL_OBJ_EQN && type != GRETL_OBJ_SYS &&
	type != GRETL_OBJ_VAR) {
	*err = E_DATA;
	goto bailout;
    }

    rset = restriction_set_new(ptr, type, opt);
    if (rset == NULL) {
	*err = E_ALLOC;
    }

    if (!*err && name == NULL) {
	*err = real_restriction_set_parse_line(rset, line, NULL, 1);
	if (*err) {
	    rset = NULL;
	    if (*err == E_PARSE) {
		gretl_errmsg_sprintf(_("parse error in '%s'\n"), line);
	    }
	}
    }

 bailout:

    free(name);

    return rset;
}

static int print_restricted_estimates (MODEL *pmod,
				       const DATASET *dset,
				       gretl_matrix *vb,
				       gretl_matrix *S,
				       double s2,
				       int nc, int rk, 
				       PRN *prn)
{
    const double *b = vb->val;
    char **names;
    double v, *se;
    int i;

    names = strings_array_new_with_length(nc, VNAMELEN);
    if (names == NULL) {
	return E_ALLOC;
    }

    se = malloc(nc * sizeof *se);
    if (se == NULL) {
	strings_array_free(names, nc);
	return E_ALLOC;
    }

    pprintf(prn, "%s:\n\n", _("Restricted estimates"));

    for (i=0; i<nc; i++) {
	v = gretl_matrix_get(S, i, i);
	if (fabs(v) < 1.0e-16) {
	    se[i] = 0.0;
	} else {
	    se[i] = sqrt(v);
	}
	gretl_model_get_param_name(pmod, dset, i, names[i]);
    }

    print_coeffs(b, se, (const char **) names, nc, pmod->dfd + rk, 
		 pmod->ci, prn);

    pputc(prn, '\n');
    pprintf(prn, "  %s = %.*g\n", _("Standard error of the regression"), 
	    GRETL_DIGITS, sqrt(s2));
    
    strings_array_free(names, nc);
    free(se);

    return 0;
}

/* Used when generating restricted estimates or bootstrapping, which
   we do for single-equation OLS/WLS only. */

static int rset_expand_R (gretl_restriction *rset, int k)
{
    gretl_matrix *R = NULL;
    double rij;
    int i, j, jj;

    R = gretl_zero_matrix_new(rset->g, k);
    if (R == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<rset->g; i++) {
	jj = 0;
	for (j=0; j<k; j++) {
	    if (rset->mask[j]) {
		rij = gretl_matrix_get(rset->R, i, jj);
		gretl_matrix_set(R, i, j, rij);
		jj++;
	    }
	}
    }

    gretl_matrix_replace(&rset->R, R);

    return 0;
}

static void rmod_check_ifc (MODEL *rmod, const gretl_matrix *y)
{
    double x1 = 0.0, x2 = 0.0;
    double reldiff;
    int t, s = 0;

    for (t=rmod->t1; t<=rmod->t2; t++) {
	if (!na(rmod->yhat[t])) {
	    x1 += rmod->yhat[t];
	    x2 += y->val[s++];
	}
    }

    reldiff = fabs((x1 - x2) / x2);
    rmod->ifc = reldiff < 9.0e-16;
}

/* respond to the --full option in case of a restriction on
   a model estimated via OLS: arrange to replace the "last model"
   with the restricted OLS estimates (and to print the model
   if wanted)
*/

static int save_restricted_model (ExecState *state,
				  MODEL *pmod,
				  const DATASET *dset,
				  const gretl_matrix *y,
				  const gretl_matrix *b,
				  const gretl_matrix *S,
				  const gretl_matrix *u,
				  double s2,
				  gretl_restriction *rset,
				  PRN *prn)
{
    MODEL *rmod;
    int i, err = 0;

    if (state == NULL) {
	return E_DATA;
    }

    rmod = gretl_model_new();
    if (rmod == NULL) {
	return E_ALLOC;
    }    

    rmod->ci = OLS;
    rmod->ncoeff = gretl_vector_get_length(b);
    rmod->nobs = gretl_vector_get_length(y);
    rmod->full_n = dset->n;
    rmod->t1 = pmod->t1;
    rmod->t2 = pmod->t2;

    rmod->dfn = pmod->dfn - rset->g;
    rmod->dfd = pmod->dfd + rset->g;
    rmod->ybar = pmod->ybar;
    rmod->sdy = pmod->sdy;
    rmod->tss = pmod->tss;

    gretl_model_smpl_init(rmod, dset);

    err = gretl_model_allocate_storage(rmod);
    
    if (!err) {
	for (i=0; i<rmod->ncoeff; i++) {
	    rmod->coeff[i] = b->val[i];
	}
	gretl_model_set_int(rmod, "restricted", 1);
	err = gretl_model_write_vcv(rmod, S);
    }

    if (!err) {
	rmod->list = gretl_list_copy(pmod->list);
	if (rmod->list == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	double yt, ut;
	int t, s = 0;

	rmod->sigma = sqrt(s2);
	rmod->ess = 0.0;

	for (t=0; t<dset->n; t++) {
	    if (!na(pmod->uhat[t])) {
		yt = gretl_vector_get(y, s);
		ut = gretl_vector_get(u, s);
		rmod->uhat[t] = ut;
		rmod->yhat[t] = yt - ut;
		rmod->ess += ut * ut;
		s++;
	    }
	}

	rmod_check_ifc(rmod, y); /* FIXME? */

	if (rmod->ifc) {
	    double tmp;

	    rmod->rsq = 1.0 - (rmod->ess / rmod->tss);
	    if (rmod->rsq < 0.0) {
		rmod->rsq = 0.0;
	    }
	    if (rmod->dfn == 0) {
		rmod->adjrsq = 0.0;
		rmod->fstt = NADBL;
	    } else {
		tmp = rmod->tss * rmod->dfd;
		rmod->adjrsq = 1 - (rmod->ess * (rmod->nobs - 1) / tmp);
		tmp = rmod->dfd / (double) rmod->dfn;
		rmod->fstt = tmp * (rmod->tss - rmod->ess) / rmod->ess;
		if (rmod->fstt < 0) {
		    rmod->fstt = NADBL;
		}
	    }
	}

	rmod->ncoeff -= rset->g;
	ls_criteria(rmod);
	rmod->ncoeff += rset->g;
    }

    if (err) {
	gretl_model_free(rmod);
    } else {
	set_model_id(rmod);
	gretl_exec_state_set_model(state, rmod);
    } 

    return err;
}

/* generate full restricted estimates: this function is used
   only for single-equation models, estimated via OLS */

static int do_restricted_estimates (ExecState *state,
				    gretl_restriction *rset,
				    const DATASET *dset,
				    PRN *prn)
{
    MODEL *pmod = rset->obj;
    gretl_matrix_block *B;
    gretl_matrix *X, *y, *b, *S;
    gretl_matrix *u = NULL;
    int *xlist = NULL;
    double s2 = 0.0;
    int T = pmod->nobs;
    int k = pmod->ncoeff;
    int i, s, t;
    int yno, err = 0;

    B = gretl_matrix_block_new(&X, T, k,
			       &y, T, 1,
			       &b, k, 1,
			       &S, k, k,
			       NULL);

    if (B == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(rset->R) != k) {
	err = rset_expand_R(rset, k);
	if (err) {
	    goto bailout;
	}
    }

    yno = gretl_model_get_depvar(pmod);
    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    if (rset->opt & OPT_F) {
	/* full save wanted */
	u = gretl_matrix_alloc(T, 1);
	if (u == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}
    }

    s = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    continue;
	}
	gretl_vector_set(y, s, dset->Z[yno][t]);
	for (i=0; i<k; i++) {
	    gretl_matrix_set(X, s, i, dset->Z[xlist[i+1]][t]);
	}
	s++;
    }

#if RDEBUG
    gretl_matrix_print(rset->R, "R (for restricted estimates)");
    gretl_matrix_print(rset->q, "q (for restricted estimates)");
#endif

    err = gretl_matrix_restricted_ols(y, X, rset->R, rset->q, 
				      b, S, u, &s2);

    if (!err) {
	if (rset->opt & OPT_F) {
	    err = save_restricted_model(state, pmod, dset, y, b, S,
					u, s2, rset, prn);
	} else {
	    print_restricted_estimates(pmod, dset, b, S, s2,
				       k, rset->g, prn);
	}
    }

 bailout:
    
    gretl_matrix_block_destroy(B);
    gretl_matrix_free(u);
    free(xlist);

    return err;
}

/* print or save restricted estimates, single equation */

static void print_or_save_result (ExecState *state,
				  gretl_restriction *rset, 
				  const DATASET *dset,
				  PRN *prn)
{
    MODEL *pmod = rset->obj;
    int robust, asym = 0;
    int quiet  = (rset->opt & OPT_Q);
    int silent = (rset->opt & OPT_S);

    robust = (pmod->opt & OPT_R);

    if (ASYMPTOTIC_MODEL(pmod->ci) || rset->rfunc != NULL) {
	asym = 1;
    }

    if (asym) {
	rset->code = GRETL_STAT_WALD_CHISQ;
	rset->pval = chisq_cdf_comp(rset->g, rset->test);
	if (!silent) {
	    pprintf(prn, "\n%s: %s(%d) = %g, ", _("Test statistic"), 
		    (robust)? _("Robust chi^2"): "chi^2",
		    rset->g, rset->test);
	}
    } else {
	rset->code = GRETL_STAT_F;
	rset->test /= rset->g;
	rset->pval = snedecor_cdf_comp(rset->g, pmod->dfd, rset->test);
	if (!silent) {
	    pprintf(prn, "\n%s: %s(%d, %d) = %g, ", _("Test statistic"), 
		    (robust)? _("Robust F"): "F",
		    rset->g, pmod->dfd, rset->test);
	}
    }

    if (!na(rset->pval) && !silent) {
	pprintf(prn, _("with p-value = %g\n"), rset->pval);
    }

    if (!silent) {
	pputc(prn, '\n');
    }

    if (!(rset->opt & OPT_C)) {
	record_test_result(rset->test, rset->pval, _("restriction"));
    }

    if (pmod != NULL && pmod->ci == OLS && 
	dset != NULL && dset->Z != NULL) {
	if ((rset->opt & OPT_F) || !(quiet || silent)) {
	    /* print and/or save restricted estimates */
	    do_restricted_estimates(state, rset, dset, prn);
	}
    }
}

/* execute the test, for a single equation */

static int test_restriction_set (gretl_restriction *rset, PRN *prn)
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

    err = check_R_matrix(rset->R);

    if (!err) {
	b = gretl_coeff_vector_from_model(pmod, rset->mask, &err);
    }

    if (!err) {
	vcv = gretl_vcv_matrix_from_model(pmod, rset->mask, &err);
    }

    if (err) {
	goto bailout;
    }

    if (rset->g == 0) {
	rset->g = gretl_matrix_rows(rset->q);
    }

    br = gretl_column_vector_alloc(rset->g);
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
	    rset->bse = sqrt(RvR->val[0]);
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

static int bootstrap_zero_restriction (gretl_restriction *rset,
				       MODEL *pmod,
				       const DATASET *dset,
				       PRN *prn)
{
    rrow *r = rset->rows[0];
    gretlopt bopt = OPT_P | OPT_R;
    int B = 0;
    int err = 0;

    if (rset->opt & OPT_S) {
	/* silent */
	bopt |= OPT_S;
    }

    gretl_restriction_get_boot_params(&B, &bopt);
    err = bootstrap_analysis(pmod, r->bnum[0], B, dset, 
			     bopt, prn);
    return err;
}

/* Are we able to do a bootstrap test of the restrictions on
   the given model?  Return non-zero if not. */

static int check_bootstrap (MODEL *pmod, PRN *prn)
{
    int ok = bootstrap_ok(pmod->ci);

    if (!ok) {
	pputs(prn, "Sorry, the bootstrap option is not supported for this test");
	pputc(prn, '\n');
    }

    return (ok)? 0 : E_NOTIMP;
}

/* Do we have a simple b[i] = 0? */

static int is_simple_zero_restriction (gretl_restriction *rset)
{
    if (rset->rows != NULL && rset->g == 1) {
	return (rset->rows[0]->nterms == 1 && rset->rows[0]->rhs == 0);
    } else {
	return 0;
    }
}

static int do_single_equation_test (ExecState *state,
				    gretl_restriction *rset,
				    const DATASET *dset,
				    PRN *prn)
{
    int err = 0;

    if (rset->rows != NULL) {
	err = restriction_set_make_mask(rset);
	if (err) {
	    return err;
	}
    }

    if (rset->opt & OPT_B) {
	/* bootstrapping, if possible */
	MODEL *pmod = rset->obj;
	
	err = check_bootstrap(pmod, prn);
	if (err) {
	    return err;
	}

	if (is_simple_zero_restriction(rset)) {
	    err = bootstrap_zero_restriction(rset, pmod, dset, prn);
	} else {
	    err = test_restriction_set(rset, prn);
	    if (!err && rset->R->cols != pmod->ncoeff) {
		err = rset_expand_R(rset, pmod->ncoeff);
	    }
	    if (!err) {
		rset->test /= rset->g;
		err = bootstrap_test_restriction(pmod, rset->R, rset->q,
						 rset->test, rset->g, 
						 dset, rset->opt, prn);
	    }
	}
    } else {
	err = test_restriction_set(rset, prn);
	if (!err) {
	    print_or_save_result(state, rset, dset, prn);
	}
    }

    return err;
}

GRETL_VAR *gretl_restricted_vecm (gretl_restriction *rset, 
				  const DATASET *dset,
				  gretlopt opt,
				  PRN *prn,
				  int *err)
{
    GRETL_VAR *jvar = NULL;

#if RDEBUG
    fprintf(stderr, "gretl_restricted_vecm: starting\n");
#endif

    if (rset == NULL || rset->otype != GRETL_OBJ_VAR) {
	*err = E_DATA;
	return NULL;
    }

    rset->opt |= opt;

    *err = restriction_set_form_matrices(rset);

    if (rset->rows != NULL && !(opt & OPT_S)) {
	print_restriction_set(rset, dset, opt, prn);
    } 

    if (!*err) {
	jvar = real_gretl_restricted_vecm(rset->obj, rset, dset, 
					  prn, err);
    }

    destroy_restriction_set(rset);

    if (jvar != NULL && *err != 0) {
	gretl_VAR_free(jvar);
	jvar = NULL;
    }

    return jvar;
}

#define RCOEFFNAME "restr__b"

static int nonlinear_wald_test (gretl_restriction *rset, gretlopt opt,
				PRN *prn)
{
    gretl_matrix *coeff = NULL;
    gretl_matrix *vcv = NULL;
    gretl_matrix *J = NULL;
    gretl_matrix *bread = NULL;
    gretl_matrix *ham = NULL;
    char fncall[64];
    int published = 0;
    int err = 0;

    if (get_user_function_by_name(rset->rfunc) == NULL) {
	gretl_errmsg_sprintf(_("The symbol '%s' is undefined\n"), rset->rfunc);
	return E_UNKVAR;
    }

    if (rset->otype == GRETL_OBJ_EQN) {
	MODEL *pmod = rset->obj;

	coeff = gretl_coeff_vector_from_model(pmod, NULL, &err);
	if (!err) {
	    vcv = gretl_vcv_matrix_from_model(pmod, NULL, &err);
	}
    } else if (rset->otype == GRETL_OBJ_SYS) {
	equation_system *sys = rset->obj;

	coeff = equation_system_get_matrix(sys, M_COEFF, &err);
	if (!err) {
	    vcv = equation_system_get_matrix(sys, M_VCV, &err);
	}
    } else {
	return E_NOTIMP; /* relax this later? */
    }

    if (!err) {
	/* "publish" the coeff matrix temporarily */
	err = private_matrix_add(coeff, RCOEFFNAME);
	if (!err) {
	    published = 1;
	}
    }

    if (!err) {
	/* formulate function call string and make the call:
	   should get a vector */
	sprintf(fncall, "%s(%s)", rset->rfunc, RCOEFFNAME);
	bread = generate_matrix(fncall, NULL, &err);
    }

#if RDEBUG
    fprintf(stderr, "nonlinear_wald_test: fncall = '%s'\n", fncall);
    gretl_matrix_print(bread, "'bread'");
#endif

    if (!err) {
	rset->g = gretl_vector_get_length(bread);
	if (rset->g == 0) {
	    err = E_NONCONF;
	}
    }

    if (!err) {
	J = fdjac(coeff, fncall, NULL, &err);
    }

#if RDEBUG
    gretl_matrix_print(J, "J");
#endif

    if (!err) {
	ham = gretl_matrix_alloc(J->rows, J->rows);
	if (ham == NULL) {
	    err = E_ALLOC;
	} else {
	    err = gretl_matrix_qform(J, GRETL_MOD_NONE, vcv, 
				     ham, GRETL_MOD_NONE);
	}
    }

    if (!err) {
	err = gretl_invpd(ham);
    }

    if (!err) {
	rset->test = gretl_scalar_qform(bread, ham, &err);
    }

    if (!err) {
	print_or_save_result(NULL, rset, NULL, prn);
    }

    gretl_matrix_free(vcv);
    gretl_matrix_free(J);
    gretl_matrix_free(bread);
    gretl_matrix_free(ham);

    if (published) {
	destroy_private_matrices();
    } else {
	gretl_matrix_free(coeff);
    }

    return err;
}

/* Respond to "end restrict": in the case of a single equation, go
   ahead and do the test; for a VECM, hand off to the driver in
   var.c; for non-VAR systems, form the restriction matrices R 
   and q and attach these to the equation system.
*/

int
gretl_restriction_finalize_full (ExecState *state,
				 gretl_restriction *rset, 
				 const DATASET *dset,
				 gretlopt opt,
				 PRN *prn)
{
    int silent, err = 0;

#if RDEBUG
    fprintf(stderr, "gretl_restriction_finalize: starting\n");
#endif

    if (rset == NULL) {
	return 1;
    }

    rset->opt |= opt;
    silent = (rset->opt & OPT_S);

    if (rset->otype == GRETL_OBJ_EQN && (opt & OPT_F)) {
	/* --full option is for OLS only */
	MODEL *pmod = rset->obj;

	if (pmod->ci != OLS) {
	    return E_BADOPT;
	}
    }

    if (rset->rfunc != NULL) {
	err = nonlinear_wald_test(rset, opt, prn);
	destroy_restriction_set(rset);
	return err;
    }

    if (rset->otype != GRETL_OBJ_EQN) {
	err = restriction_set_form_matrices(rset);
	if (err) {
	    destroy_restriction_set(rset);
	    return err;
	}
    }

    if (!silent && rset->rows != NULL) {
	print_restriction_set(rset, dset, opt, prn);
    }

    if (rset->otype == GRETL_OBJ_VAR) {
	err = gretl_VECM_test(rset->obj, rset, dset, rset->opt, prn);
    } else if (rset->otype == GRETL_OBJ_SYS) {
	if (rset->opt & OPT_W) {
	    err = system_wald_test(rset->obj, rset->R, rset->q, opt, prn);
	} else {
	    system_set_restriction_matrices(rset->obj, rset->R, rset->q);
	    rset->R = NULL;
	    rset->q = NULL;
	}
    } else {
	/* single-equation model */
	err = do_single_equation_test(state, rset, dset, prn);
    }

    if (!(rset->opt & OPT_C)) {
	destroy_restriction_set(rset);
    }

    return err;
}

int gretl_restriction_finalize (gretl_restriction *rset, 
				const DATASET *dset,
				gretlopt opt,
				PRN *prn)
{
    return gretl_restriction_finalize_full(NULL, rset, dset,
					   opt, prn);
}

/**
 * gretl_sum_test:
 * @list: list of variables to use.
 * @pmod: pointer to model.
 * @dset: information on the data set.
 * @prn: gretl printing struct.
 * 
 * Calculates the sum of the coefficients, relative to the given model, 
 * for the variables given in @list.  Prints this estimate along 
 * with its standard error.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int 
gretl_sum_test (const int *list, MODEL *pmod, DATASET *dset,
		PRN *prn)
{
    gretl_restriction *r;
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

    if (exact_fit_check(pmod, prn)) {
	return 0;
    }

    r = restriction_set_new(pmod, GRETL_OBJ_EQN, OPT_Q | OPT_C);
    if (r == NULL) {
	return 1;
    }

    *line = '\0';
    len = 0;

    for (i=1; i<=list[0]; i++) {
	sprintf(bstr, "b[%s]", dset->varname[list[i]]);
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
	err = real_restriction_set_parse_line(r, line, dset, 1); 
    }

    if (!err) {
	err = gretl_restriction_finalize(r, dset, OPT_NONE, NULL);
    }

    if (!err) {
	double test;

	pprintf(prn, "\n%s: ", _("Variables"));

	for (i=1; i<=list[0]; i++) {
	    pprintf(prn, "%s ", dset->varname[list[i]]);
	}

	pprintf(prn, "\n   %s = %g\n", _("Sum of coefficients"), r->bsum);

	if (r->code == GRETL_STAT_F) {
	    pprintf(prn, "   %s = %g\n", _("Standard error"), r->bse);
	    test = sqrt(r->test);
	    if (r->bsum < 0) {
		test = -test;
	    }
	    pprintf(prn, "   t(%d) = %g ", pmod->dfd, test);
	    pprintf(prn, _("with p-value = %g\n"), r->pval);
	    record_test_result(test, r->pval, _("sum")); 
	} else if (r->code == GRETL_STAT_WALD_CHISQ) {
	    pprintf(prn, "   %s = %g\n", _("Standard error"), r->bse);
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

gretlopt gretl_restriction_get_options (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->opt : OPT_NONE;
}

GretlObjType gretl_restriction_get_type (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->otype : GRETL_OBJ_NULL;
}

const gretl_matrix *rset_get_R_matrix (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->R : NULL;
}

const gretl_matrix *rset_get_q_matrix (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->q : NULL;
}

const gretl_matrix *rset_get_Ra_matrix (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->Ra : NULL;
}

const gretl_matrix *rset_get_qa_matrix (const gretl_restriction *rset)
{
    return (rset != NULL)? rset->qa : NULL;
}

int rset_VECM_bcols (const gretl_restriction *rset)
{
    return (rset == NULL || rset->R == NULL)? 0 : rset->R->cols;
}

int rset_VECM_acols (const gretl_restriction *rset)
{
    return (rset == NULL || rset->Ra == NULL)? 0 : rset->Ra->cols;
}

void rset_add_results (gretl_restriction *rset,
		       double test, double pval,
		       double lnl)
{
    rset->test = test;
    rset->pval = pval;
    rset->lnl = lnl;
}

void rset_record_LR_result (gretl_restriction *rset)
{
    record_LR_test_result(rset->test, rset->pval, rset->lnl,
			  "LR");
}

