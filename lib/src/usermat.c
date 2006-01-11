/*
 *  Copyright (c) by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"
#include "gretl_func.h"
#include "usermat.h"

#include <errno.h>

#define MDEBUG 0

#define MNAMELEN 32

typedef struct user_matrix_ user_matrix;

struct user_matrix_ {
    gretl_matrix *M;
    int level;
    char name[MNAMELEN];
};

static user_matrix **matrices;
static int n_matrices;

static const char *varchars = "acdefghijklmnopqrstuvwxyz"
                              "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                              "0123456789_";

static user_matrix *user_matrix_new (gretl_matrix *M, const char *name)
{
    user_matrix *u;

    u = malloc(sizeof *u);
    if (u == NULL) {
	return NULL;
    }

    u->M = M;

    u->level = gretl_function_stack_depth();

    *u->name = '\0';
    strncat(u->name, name, MNAMELEN - 1);

    return u;
}

static int add_user_matrix (gretl_matrix *M, const char *name)
{
    user_matrix **tmp;

    if (!strcmp(name, "I") ||
	!strcmp(name, "zeros") ||
	!strcmp(name, "ones")) {
	sprintf(gretl_errmsg, "The matrix name '%s' is reserved", name);
	return 1;
    }

    if (check_varname(name)) {
	return E_DATA;
    }

    tmp = realloc(matrices, (n_matrices + 1) * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    matrices = tmp;
    matrices[n_matrices] = user_matrix_new(M, name);

    if (matrices[n_matrices] == NULL) {
	return E_ALLOC;
    }

    n_matrices++;

#if MDEBUG
    fprintf(stderr, "add_user_matrix: allocated '%s' at %p (M at %p, %dx%d)\n",
	    name, (void *) matrices[n_matrices-1], (void *) M,
	    gretl_matrix_rows(M), gretl_matrix_cols(M));
#endif

    return 0;
}

static int replace_user_matrix (user_matrix *u, gretl_matrix *M,
				gretl_matrix **R)
{
    if (R != NULL) {
	*R = u->M;
    }

    gretl_matrix_free(u->M);
    u->M = M;

    return 0;
}

int is_user_matrix (gretl_matrix *m)
{
    int i;

    for (i=0; i<n_matrices; i++) {
	if (m == matrices[i]->M) {
	    return 1;
	}
    }

    return 0;
}

static gretl_matrix *
real_get_matrix_by_name (const char *name, int transpose_ok)
{
    char test[MNAMELEN];
    int level = gretl_function_stack_depth();
    int transp = 0;
    int i;

    *test = '\0';
    strncat(test, name, MNAMELEN - 1);

    if (transpose_ok) {
	if (test[strlen(test) - 1] == '\'') {
	    test[strlen(test) - 1] = '\0';
	    transp = 1;
	}
    } 

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(test, matrices[i]->name) &&
	    matrices[i]->level == level) {
	    if (transp) {
		return gretl_matrix_copy_transpose(matrices[i]->M);
	    } else {
		return matrices[i]->M;
	    }
	}
    }

    return NULL;
}

/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

gretl_matrix *get_matrix_by_name (const char *name)
{
    return real_get_matrix_by_name(name, 1);
}

int named_matrix_get_scalar (const char *name, double *x)
{
    gretl_matrix *m = real_get_matrix_by_name(name, 0);
    int err = 0;

    if (m == NULL || gretl_matrix_rows(m) != 1 ||
	gretl_matrix_cols(m) != 1) {
	err = E_UNKVAR;
    } else {
	*x = gretl_matrix_get(m, 0, 0);
    }

    return err;
}

static gretl_matrix *original_matrix_by_name (const char *name)
{
    return real_get_matrix_by_name(name, 0);
}

static user_matrix *get_user_matrix_by_name (const char *name)
{
    int level = gretl_function_stack_depth();
    int i;

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(name, matrices[i]->name) &&
	    matrices[i]->level == level) {
	    return matrices[i];
	}
    }

    return NULL;
}

static int 
row_or_column_selection (const char *s, int *idx, int *err)
{
    int n = 0;

    if (strchr(s, ':')) {
	if (sscanf(s, "%d:%d", &idx[0], &idx[1]) != 2) {
	    *err = 1;
	} else if (idx[0] <= 0 || idx[1] <= 0) {
	    *err = 1;
	} else {
	    n = idx[1] - idx[0] + 1;
	}
    } else if (sscanf(s, "%d", &idx[0]) != 1) {
	*err = 1;
    } else if (idx[0] <= 0) {
	*err = 1;
    } else {
	idx[1] = idx[0];
	n = 1;
    }

    return n;
}

/* This supports the extraction of scalars, vectors or sub-matrices.
   E.g. B[1,2] extracts a scalar, B[1,] extracts a row vector, and
   B[,2] extracts a column vector.  B[1:2,2:3] extracts a sub-matrix
   composed of the intersection of rows 1 and 2 with columns 2 and 3.
*/

static gretl_matrix *
real_matrix_get_slice (gretl_matrix *M, const char *s, int *err)
{
    gretl_matrix *S = NULL;
    int m = gretl_matrix_rows(M);
    int n = gretl_matrix_cols(M);
    char tmp[32] = {0};
    int r[2] = {-1, -1};
    int c[2] = {-1, -1};
    int nr = 0, nc = 0;
    const char *p;
    char *q = NULL;
    int len;

    p = strrchr(s, ']');
    if (p == NULL) {
	*err = E_SYNTAX;
	return NULL;
    }

    if (strchr(s, ',') == NULL) {
	*err = E_SYNTAX;
	return NULL;
    }    

    s = strchr(s, '[') + 1;
    len = p - s;

    if (len > 31) {
	*err = 1;
	return NULL;
    } 

    strncat(tmp, s, len);
    p = tmp;
    while (*p == ' ') p++;

    if (*p != ',') {
	/* left-hand (row selection) block is present */
	q = strchr(p, ',');
	*q = '\0';
	nr = row_or_column_selection(p, r, err);
    } else {
	/* default: whole range */
	r[0] = 1;
	r[1] = m;
	nr = m;
    }

    if (q != NULL) {
	/* undo masking of right-hand expression */
	*q = ',';
    }

    if (!*err) {
	/* skip to right-hand block, if any */
	p += strspn(p, ":0123456789");
	while (*p == ' ' || *p == ',') p++;
	if (*p) {
	    /* right-hand (column selection) block is present */
	    nc = row_or_column_selection(p, c, err);
	} else {
	    /* default: whole range */
	    c[0] = 1;
	    c[1] = n;
	    nc = n;
	}
    }

    if (nr < 1 || nc < 1) {
	*err = 1;
    } else if (r[0] > m || c[0] > n) {
	*err = 1;
    } else if (r[1] > m || c[1] > n) {
	*err = 1;
    }

    if (!*err) {
	S = gretl_matrix_alloc(nr, nc);
	if (S == NULL) {
	    *err = E_ALLOC;	
	}
    }

    if (S != NULL) {
	double x;
	int i, j, l, k = 0;

	for (i=r[0]-1; i<r[1]; i++) {
	    l = 0;
	    for (j=c[0]-1; j<c[1]; j++) {
		x = gretl_matrix_get(M, i, j);
		gretl_matrix_set(S, k, l++, x);
	    }
	    k++;
	}
    }
	
    return S;
}

gretl_matrix *user_matrix_get_slice (const char *s, int *err)
{
    gretl_matrix *M = NULL;
    gretl_matrix *S = NULL;
    char test[MNAMELEN];
    int len = strcspn(s, "[");

    if (len < MNAMELEN) {
	*test = '\0';
	strncat(test, s, len);
	M = real_get_matrix_by_name(test, 0);
	if (M != NULL) {
	    S = real_matrix_get_slice(M, s, err);
	}
    }

    return S;
}

/**
 * add_or_replace_user_matrix:
 * @M: gretl matrix.
 * @name: name for the matrix.
 * @R: location to receive address of matrix that was
 * replaced, if any (or %NULL).
 *
 * Checks whether a matrix of the given @name already exists.
 * If so, the original matrix is replaced by @M; if not, the
 * the matrix @M is added to the stack of user-defined
 * matrices.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int add_or_replace_user_matrix (gretl_matrix *M, const char *name,
				gretl_matrix **R)
{
    user_matrix *u;
    int err = 0;

    u = get_user_matrix_by_name(name);
    if (u != NULL) {
	err = replace_user_matrix(u, M, R);
    } else {
	err = add_user_matrix(M, name);
    }

    return err;
}

static void destroy_user_matrix (user_matrix *u)
{
#if 0
    fprintf(stderr, "destroy_user_matrix: dummy call\n");
#else
#if MDEBUG
    fprintf(stderr, "destroy_user_matrix: freeing matrix at %p\n", 
	    (void *) u->M);
#endif
    gretl_matrix_free(u->M);
    free(u);
#endif
}

/**
 * destroy_user_matrices_at_level:
 * @level: stack level of function execution.
 *
 * Destroys and removes from the stack of user matrices all
 * matrices that were created at the given @level.  This is 
 * part of the cleanup that is performed when a user-defined
 * function terminates.
 *
 * Returns: 0 on success, non-zero on error.
 */

int destroy_user_matrices_at_level (int level)
{
    user_matrix **tmp;
    int i, j, nm = 0;
    int err = 0;

    for (i=0; i<n_matrices; i++) {
	if (matrices[i]->level == level) {
	    destroy_user_matrix(matrices[i]);
	    for (j=i; j<n_matrices - 1; j++) {
		matrices[j] = matrices[j+1];
	    }
	    matrices[n_matrices - 1] = NULL;
	} else {
	    nm++;
	}
    }

    if (nm < n_matrices) {
	n_matrices = nm;
	if (nm == 0) {
	    free(matrices);
	    matrices = NULL;
	} else {
	    tmp = realloc(matrices, nm * sizeof *tmp);
	    if (tmp == NULL) {
		err = E_ALLOC;
	    } else {
		matrices = tmp;
	    }
	}
    }

    return err;
}

/**
 * destroy_user_matrices:
 *
 * Frees all resources associated with the stack of user-
 * defined matrices.
 */

void destroy_user_matrices (void)
{
    int i;

#if MDEBUG
    fprintf(stderr, "destroy_user_matrices called, n_matrices = %d\n",
	    n_matrices);
#endif

    if (matrices == NULL) {
	return;
    }

    for (i=0; i<n_matrices; i++) {
#if MDEBUG
	fprintf(stderr, "destroying user_matrix %d (%s) at %p\n", i,
		matrices[i]->name, (void *) matrices[i]);
#endif
	destroy_user_matrix(matrices[i]);
    }

    free(matrices);
    matrices = NULL;
    n_matrices = 0;
}

static int scalar_varnum (const char *s, const DATAINFO *pdinfo)
{
    char vname[VNAMELEN];
    int len = strspn(s, varchars);
    int v = 0;

    if (len > 0 && len < VNAMELEN) {
	*vname = '\0';
	strncat(vname, s, len);
	v = varindex(pdinfo, vname);
	if (v == pdinfo->v || pdinfo->vector[v]) {
	    v = 0;
	} 
    }

    return v;
}

#define numeric_start(s) (strcspn(s, "+-0123456789") == 0)

static int 
get_rows_cols (const char *str, const DATAINFO *pdinfo,
	       int *r, int *c, int *scalars)
{
    const char *p = str;
    char *s;
    int nf0 = 0, nf1;
    int i, len;
    int err = 0;

    *r = 1;
    while (*p) {
	if (*p == ';') {
	    *r += 1;
	}
	p++;
    }

    for (i=0; i<*r && !err; i++) {
	len = strcspn(str, ";");

	s = gretl_strndup(str, len);
	if (s == NULL) {
	    err = E_ALLOC;
	    break;
	}

	/* clean up */
	charsub(s, ',', ' ');
	charsub(s, '}', ' ');
	charsub(s, '\'', ' ');

	nf1 = count_fields(s);
	if (i > 0 && nf1 != nf0) {
	    strcpy(gretl_errmsg, "Inconsistent matrix specification");
	    err = 1;
	    break;
	}

	nf0 = nf1;

	/* peek at first field: is it a scalar? */
	if (i == 0) {
	    p = s;
	    while (isspace(*p)) p++;
	    if (numeric_start(p)) {
		*scalars = 1;
	    } else if (isalpha(*p)) {
		*scalars = scalar_varnum(p, pdinfo);
	    }
	}

	str += len + 1;

	free(s);
    }

    if (!err) {
	*c = nf0;
    }

    return err;
}

static int 
get_varnum (const char **s, const DATAINFO *pdinfo, int *err)
{
    char vname[12];
    int v = 0;

    if (sscanf(*s, "%11s", vname)) {
	int len = strlen(vname);
	char *p = strrchr(vname, '}');
	
	if (p != NULL) {
	    *p = '\0';
	}

	v = varindex(pdinfo, vname);
	if (v < pdinfo->v) {
	    *s += len;
	} else {
	    sprintf(gretl_errmsg, _("Unknown variable '%s'"), vname);
	    *err = 1;
	} 
    } else {
	*err = 1;
    }

    return v;
}

static double get_var_double (const char **s, const double **Z,
			      const DATAINFO *pdinfo, int *err)
{
    char vname[VNAMELEN];
    double x = NADBL;
    int v, len = strspn(*s, varchars);

    if (len > VNAMELEN - 1) {
	*err = 1;
    } else {
	*vname = '\0';
	strncat(vname, *s, len);
	v = varindex(pdinfo, vname);
	if (v == pdinfo->v) {
	    *err = E_UNKVAR;
	} else if (pdinfo->vector[v]) {
	    *err = E_DATA;
	} else {
	    x = Z[v][0];
	    if (na(x)) {
		*err = E_MISSDATA;
	    } else {
		*s += len;
	    }
	}
    }

    return x;
}

static double get_numeric_double (const char **s, int *err)
{
    double x = NADBL;
    char *p;

    x = strtod(*s, &p);

    if (!strcmp(*s, p)) {
	sprintf(gretl_errmsg, _("'%s' -- no numeric conversion performed!"), *s);
	*err = 1;
    } else if (*p != '\0' && *p != ',' && *p != ';' && *p != ' ' && *p != '}') {
	if (isprint(*p)) {
	    sprintf(gretl_errmsg, _("Extraneous character '%c' in data"), *p);
	} else {
	    sprintf(gretl_errmsg, _("Extraneous character (0x%x) in data"), *p);
	}
	*err = 1;
    } else if (errno == ERANGE) {
	sprintf(gretl_errmsg, _("'%s' -- number out of range!"), *s);
	*err = 1;
    }

    *s = p;

    return x;
}

static int 
fill_matrix_from_scalars (gretl_matrix *M, const char *s, 
			  int r, int c, int transp,
			  const double **Z, const DATAINFO *pdinfo)
{
    double x;
    int i, j;
    int err = 0;

    if (transp) {
	for (j=0; j<c && !err; j++) {
	    for (i=0; i<r && !err; i++) {
		if (isalpha(*s)) {
		    x = get_var_double(&s, Z, pdinfo, &err);
		} else {
		    x = get_numeric_double(&s, &err);
		}
		if (!err) {
		    gretl_matrix_set(M, i, j, x);
		    s += strspn(s, " ,;");
		}
	    }
	}
    } else {
	for (i=0; i<r && !err; i++) {
	    for (j=0; j<c && !err; j++) {
		if (isalpha(*s)) {
		    x = get_var_double(&s, Z, pdinfo, &err);
		} else {
		    x = get_numeric_double(&s, &err);
		}
		if (!err) {
		    gretl_matrix_set(M, i, j, x);	
		    s += strspn(s, " ,;");
		}
	    }
	}
    }

    return err;
}

/* Fill a matrix with values from one or more series: since
   gretl's matrix methods cannot handle missing values, it is
   an error if any missing values are encountered in the
   given series.  Each series occupies a column by default
   (unless transp is non-zero).
*/

static int fill_matrix_from_series (gretl_matrix *M, const char *s,
				    int r, int c, int transp,
				    const double **Z, const DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int nvr = r, nvc = c;
    int i, j, k, t, v;
    double x;
    int err = 0;

    s += strspn(s, " ");

    if (transp) {
	nvr /= T;
	for (j=0; j<c && !err; j++) {
	    for (i=0; i<nvr && !err; i++) {
		v = get_varnum(&s, pdinfo, &err);
		if (!err) {
		    for (t=0; t<T; t++) {
			if (pdinfo->vector[v]) {
			    x = Z[v][t + pdinfo->t1];
			} else {
			    x = Z[v][0];
			}
			if (na(x)) {
			    err = E_MISSDATA;
			} else {
			    k = i * T + t;
			    gretl_matrix_set(M, k, j, x);
			}
		    }
		    s += strspn(s, " ,;");
		}
	    }
	}
    } else {
	nvc /= T;
	for (i=0; i<r && !err; i++) {
	    for (j=0; j<nvc && !err; j++) {
		v = get_varnum(&s, pdinfo, &err);
		if (!err) {
		    for (t=0; t<T; t++) {
			if (pdinfo->vector[v]) {
			    x = Z[v][t + pdinfo->t1]; 
			} else {
			    x = Z[v][0];
			}
			if (na(x)) {
			    err = E_MISSDATA;
			} else {			
			    k = j * T + t;
			    gretl_matrix_set(M, i, k, x);
			}
		    }
		    s += strspn(s, " ,;");
		}
	    }
	}
    } 

    return err;
}

static int matrix_genr (const char *name, const char *s, 
			double ***pZ, DATAINFO *pdinfo)
{
    char genline[MAXLINE];
    int err;

    if (strlen(name) + strlen(s) + 7 > MAXLINE - 1) {
	err = 1;
    } else {
	sprintf(genline, "genr %s %s", name, s);
	err = generate(genline, pZ, pdinfo, OPT_M);
    }

    return err;
}

gretl_matrix *fill_matrix_from_list (const char *s, const double **Z,
				     const DATAINFO *pdinfo, int transp,
				     int *err)
{
    gretl_matrix *M = NULL;
    char word[32];
    char *mask = NULL;
    const int *list;
    int len;

    while (isspace(*s)) s++;

    len = strspn(s, varchars);
    if (len == 0 || len > 31) {
	return NULL;
    }

    *word = '\0';
    strncat(word, s, len);
    list = get_list_by_name(word);
    if (list == NULL) {
	return NULL;
    }

    M = gretl_matrix_data_subset(list, Z, pdinfo->t1, pdinfo->t2, &mask);

    if (M == NULL) {
	*err = 1;
    }

    if (mask != NULL) {
	*err = E_MISSDATA;
	free(mask);
	gretl_matrix_free(M);
	M = NULL;
    }

    if (M != NULL && transp) {
	gretl_matrix *R = gretl_matrix_copy_transpose(M);

	if (R == NULL) {
	    *err = E_ALLOC;
	    gretl_matrix_free(M);
	    M = NULL;
	} else {
	    gretl_matrix_free(M);
	    M = R;
	}
    }

    return M;
}

/* Currently we can create a user matrix in any one of four ways (but
   we can't mix these in a single matrix specification).

   1. Full specification of scalar elements, either numerical values or
      by reference to scalar variables.

   2. Specification of individual data series to place in the matrix.

   3. Specification of one named list of variables to place in matrix.

   4. Use of a "genr"-type expression referring to existing matrices
      and/or the I(n) function which gives an n x n identity matrix.
*/

static int create_matrix (const char *name, const char *s, 
			  double ***pZ, DATAINFO *pdinfo,
			  PRN *prn)
{
    gretl_matrix *M = NULL;
    char *p;
    int nm = n_matrices;
    int scalars = 0;
    int transp = 0;
    int r = 0, c = 0;
    int err = 0;

    p = strchr(s, '{');
    if (p == NULL) {
	err = matrix_genr(name, s, pZ, pdinfo);
	goto finalize;
    }
    s = p + 1;

    if (!err) {
	p = strchr(s, '}');
	if (p == NULL) {
	    err = 1;
	}
	if (*(p+1) == '\'') {
	    transp = 1;
	}
    }

    if (!err) {
	M = fill_matrix_from_list(s, (const double **) *pZ, pdinfo, 
				  transp, &err);
	if (err) {
	    goto finalize;
	} else if (M != NULL) {
	    err = add_or_replace_user_matrix(M, name, NULL);
	    goto finalize;
	}
    }

    if (!err) {
	err = get_rows_cols(s, pdinfo, &r, &c, &scalars);
	if (!err && c == 0) {
	    err = 1;
	}
    }

    if (!err && !scalars) {
	c *= pdinfo->t2 - pdinfo->t1 + 1;
    }

    if (!err && transp) {
	int tmp = r;

	r = c;
	c = tmp;
    }

#if MDEBUG
    fprintf(stderr, "r=%d, c=%d, transp=%d, scalars=%d\n",
	    r, c, transp, scalars);
#endif

    if (!err) {
	M = gretl_matrix_alloc(r, c);
	if (M == NULL) {
	   err = E_ALLOC;
	} 
    }

    if (!err) {
	if (scalars) {
	    err = fill_matrix_from_scalars(M, s, r, c, transp,
					   (const double **) *pZ,
					   pdinfo);
	} else {
	    err = fill_matrix_from_series(M, s, r, c, transp, 
					  (const double **) *pZ, 
					  pdinfo);
	}
    }
    
    if (!err) {
	err = add_or_replace_user_matrix(M, name, NULL);
    }

 finalize:

    if (!err) {
	if (n_matrices > nm) {
	    pprintf(prn, "Added matrix '%s'\n", name);
	} else {
	    pprintf(prn, "Replaced matrix '%s'\n", name);
	}
    } else {
	pprintf(prn, "Error adding matrix '%s'\n", name);
    }

    return err;
}

static int print_matrix_by_name (const char *name, PRN *prn)
{
    gretl_matrix *M;
    int err = 0;

    M = get_matrix_by_name(name);
    if (M == NULL) {
	pprintf(prn, _("'%s': no such matrix\n"), name);
	err = 1;
    } else {
	gretl_matrix_print_to_prn(M, name, prn);
	if (!original_matrix_by_name(name)) {
	    /* we got a transpose, created on the fly */
	    gretl_matrix_free(M);
	}
    }

    return err;
}

/**
 * matrix_command:
 * @s: string that specifies matrix command.
 * @prn: pointer to printing struct.
 *
 * To be written.
 * 
 * Returns: 0 on success, non-zero code on failure.
 */

int matrix_command (const char *line, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    char name[MNAMELEN];
    char word[9];
    int err = 0;

    if (!strncmp(line, "matrix ", 7)) line += 7;
    while (isspace(*line)) line++;

    if (!sscanf(line, "%31s", name)) {
	return 1;
    } 

    line += strlen(name);
    while (isspace(*line)) line++;

    if (*line == '=') {
	/* defining a matrix */
	err = create_matrix(name, line, pZ, pdinfo, prn);
    } else {
	*word = '\0';
	sscanf(line, "%8s", word);
	if (*word == '\0' || !strcmp(word, "print")) {
	    err = print_matrix_by_name(name, prn);
	} else {
	    /* no other commands available yet */
	    err = 1;
	}
    } 

    return err;
}

/* for use in genr, for matrices */

gretl_matrix *matrix_calc_AB (gretl_matrix *A, gretl_matrix *B, 
			      char op, int *err) 
{
    gretl_matrix *C = NULL;
    double x;
    int ra, ca;
    int rb, cb;

    *err = 0;

#if MDEBUG
    fprintf(stderr, "\n*** matrix_calc_AB: A = %p, B = %p, ", 
	    (void *) A, (void *) B);
    if (isprint(op)) fprintf(stderr, "op='%c'\n", op);
    else fprintf(stderr, "op=%d\n", op);
    debug_print_matrix(A, "input A");
    debug_print_matrix(B, "input B");
#endif

    switch (op) {
    case '\0':
	C = B;
	break;
    case '+':
    case '-':
	C = gretl_matrix_copy(A);
	if (C == NULL) {
	    *err = E_ALLOC;
	} else if (op == '+') {
	    *err = gretl_matrix_add_to(C, B);
	} else {
	    *err = gretl_matrix_subtract_from(C, B);
	}
	break;
    case '~':
	/* column-wise concatenation */
	C = gretl_matrix_col_concat(A, B, err);
	break;
    case '*':
	ra = gretl_matrix_rows(A);
	ca = gretl_matrix_cols(A);
	rb = gretl_matrix_rows(B);
	cb = gretl_matrix_cols(B);

	if (ra == 1 && ca == 1) {
	    C = gretl_matrix_copy(B);
	    if (C == NULL) {
		*err = E_ALLOC;	
	    } else {
		x = gretl_matrix_get(A, 0, 0);
		gretl_matrix_multiply_by_scalar(C, x);
	    }
	} else if (rb == 1 && cb == 1) {
	    C = gretl_matrix_copy(A);
	    if (C == NULL) {
		*err = E_ALLOC;	
	    } else {
		x = gretl_matrix_get(B, 0, 0);
		gretl_matrix_multiply_by_scalar(C, x);
	    }
	} else {
	    C = gretl_matrix_alloc(ra, cb);
#if MDEBUG
	    fprintf(stderr, "matrix_calc_AB: allocated 'C' (%dx%d) at %p\n", 
		    ra, cb, (void *) C);
#endif
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else {
		*err = gretl_matrix_multiply(A, B, C);
	    }
	}
	break;
    case OP_DOTMULT:
	/* element-wise multiplication */
	C = gretl_matrix_dot_multiply(A, B, err);
	break;
    case OP_DOTDIV:
	/* element-wise division */
	C = gretl_matrix_dot_divide(A, B, err);
	break;
    case OP_DOTPOW:
	/* element-wise exponentiation */
	if (gretl_matrix_rows(B) != 1 ||
	    gretl_matrix_cols(B) != 1) {
	    *err = 1;
	} else {
	    C = gretl_matrix_copy(A);
	    if (C == NULL) {
		*err = E_ALLOC;
	    } else {
		x = gretl_matrix_get(B, 0, 0);
		gretl_matrix_dot_pow(C, x);
	    }
	}
	break;
    case OP_KRON:
	/* Kronecker product */
	C = gretl_matrix_kronecker_product(A, B);
	if (C == NULL) {
	    *err = E_ALLOC;
	}
	break;
    default:
	*err = 1;
	break;
    } 

    return C;
}

static double 
real_user_matrix_get_determinant (gretl_matrix *m, int log)
{
    double d = NADBL;

    if (m != NULL) {
	gretl_matrix *tmp = gretl_matrix_copy(m);

	if (tmp != NULL) {
	    if (log) {
		d = gretl_matrix_log_determinant(tmp);
	    } else {
		d = gretl_matrix_determinant(tmp);
	    }
	    gretl_matrix_free(tmp);
	}
    }

    return d;
}

double user_matrix_get_determinant (gretl_matrix *m)
{
    return real_user_matrix_get_determinant(m, 0);
}

double user_matrix_get_log_determinant (gretl_matrix *m)
{
    return real_user_matrix_get_determinant(m, 1);
}

static gretl_matrix *
real_user_matrix_get_determinant_as_matrix (gretl_matrix *m, int log)
{
    gretl_matrix *dm = NULL;
    double d;

    if (log) {
	d = user_matrix_get_log_determinant(m);
    } else {
	d = user_matrix_get_determinant(m);
    }

    if (!na(d)) {
	dm = gretl_matrix_alloc(1, 1);
	if (dm != NULL) {
	    gretl_matrix_set(dm, 0, 0, d);
	}
    }

    return dm;
}

gretl_matrix *user_matrix_get_determinant_as_matrix (gretl_matrix *m)
{
    return real_user_matrix_get_determinant_as_matrix(m, 0);
}

gretl_matrix *user_matrix_get_log_determinant_as_matrix (gretl_matrix *m)
{
    return real_user_matrix_get_determinant_as_matrix(m, 1);
}

gretl_matrix *user_matrix_get_inverse (gretl_matrix *m)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    if (gretl_invert_general_matrix(R)) {
		gretl_matrix_free(R);
		R = NULL;
	    }
	}
    }

    return R;
}

gretl_matrix *
user_matrix_get_transformation (gretl_matrix *m, GretlMathFunc fn)
{
    gretl_matrix *R = NULL;

    if (m != NULL) {
	R = gretl_matrix_copy(m);
	if (R != NULL) {
	    if (gretl_matrix_transform_elements(R, fn)) {
		gretl_matrix_free(R);
		R = NULL;
	    }
	}
    }

    return R;
}  

/* move tranpose symbol ' in front of parenthesized
   matrix expression so genr can handle it as a function
*/

int reposition_transpose_symbol (char *s)
{
    int pc, len = strlen(s);
    int offset;
    int i, j, sz;
    int err = 0;

    for (i=3; i<len; i++) {
	if (s[i] == '\'' && s[i-1] == ')') {
	    pc = sz = 1;
	    /* back up to matching left paren */
	    for (j=i-2; j>=0; j--) {
		if (s[j] == ')') {
		    pc++;
		} else if (s[j] == '(') {
		    pc--;
		}
		sz++;
		if (pc == 0) {
		    offset = i - sz;
		    memmove(s + offset + 1, s + offset, sz);
		    s[offset] = '\'';
		    i++;
		    break;
		}
	    }
	    if (j <= 0 && pc != 0) {
		err = E_UNBAL;
		break;
	    }
	}
    }

    return err;
}

		    


