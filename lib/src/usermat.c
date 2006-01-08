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
#include "usermat.h"

#include <errno.h>

#define MNAMELEN 32

typedef struct user_matrix_ user_matrix;

struct user_matrix_ {
    gretl_matrix *M;
    char name[MNAMELEN];
};

user_matrix **matrices;
int n_matrices;

static user_matrix *user_matrix_new (gretl_matrix *M, const char *name)
{
    user_matrix *u;

    u = malloc(sizeof *u);
    if (u == NULL) {
	return NULL;
    }

    u->M = M;

    *u->name = '\0';
    strncat(u->name, name, MNAMELEN - 1);

    return u;
}

static int add_user_matrix (gretl_matrix *M, const char *name)
{
    user_matrix **tmp;

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

    return 0;
}

static int replace_user_matrix (user_matrix *u, gretl_matrix *M)
{
    gretl_matrix_free(u->M);
    u->M = M;

    return 0;
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
    int i;

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(name, matrices[i]->name)) {
	    return matrices[i]->M;
	}
    }

    return NULL;
}

static user_matrix *get_user_matrix_by_name (const char *name)
{
    int i;

    for (i=0; i<n_matrices; i++) {
	if (!strcmp(name, matrices[i]->name)) {
	    return matrices[i];
	}
    }

    return NULL;
}

/**
 * add_or_replace_user_matrix:
 * @M: gretl matrix.
 * @name: name for the matrix.
 *
 * Checks whether a matrix of the given @name already exists.
 * If so, the original matrix is replaced by @M; if not, the
 * the matrix @M is added to the stack of user-defined
 * matrices.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int add_or_replace_user_matrix (gretl_matrix *M, const char *name)
{
    user_matrix *u;
    int err = 0;

    u = get_user_matrix_by_name(name);
    if (u != NULL) {
	err = replace_user_matrix(u, M);
    } else {
	err = add_user_matrix(M, name);
    }

    return err;
}

static void destroy_user_matrix (user_matrix *u)
{
    gretl_matrix_free(u->M);
    free(u);
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

    for (i=0; i<n_matrices; i++) {
	destroy_user_matrix(matrices[i]);
    }

    free(matrices);
    matrices = NULL;
    n_matrices = 0;
}

static int get_rows_cols (char *str, int *r, int *c)
{
    char *s = str;
    int len;

    *r = 1;
    while (*s) {
	if (*s == ';') {
	    *r += 1;
	}
	s++;
    }

    len = strcspn(str, ";");

    s = gretl_strndup(str, len);
    if (s == NULL) {
	return E_ALLOC;
    }

    /* impose space-separation */
    charsub(s, ',', ' ');

    *c = count_fields(s);

    free(s);

    return 0;
}

static double get_double (char **s, int *err)
{
    double x = NADBL;
    char *p;

    x = strtod(*s, &p);

    if (!strcmp(*s, p)) {
	sprintf(gretl_errmsg, _("'%s' -- no numeric conversion performed!"), *s);
	*err = 1;
    } else if (*p != '\0' && *p != ',' && *p != ';' && *p != ' ') {
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

static int create_matrix (const char *name, char *s, PRN *prn)
{
    gretl_matrix *M = NULL;
    char *p;
    double x;
    int nm = n_matrices;
    int transp = 0;
    int r = 0, c = 0;
    int i, j;
    int err = 0;

    p = strchr(s, '{');
    if (p == NULL) {
	err = 1;
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
	*p = '\0';
    }

    if (!err) {
	err = get_rows_cols(s, &r, &c);
	if (!err && c == 0) {
	    err = 1;
	}
    }

    if (!err && transp) {
	int tmp = r;

	r = c;
	c = tmp;
    }

    if (!err) {
	M = gretl_matrix_alloc(r, c);
	if (M == NULL) {
	   err = E_ALLOC;
	} 
    }
    
    if (!err) {
	if (transp) {
	    for (j=0; j<c && !err; j++) {
		for (i=0; i<r && !err; i++) {
		    x = get_double(&s, &err);
		    if (!err) {
			gretl_matrix_set(M, i, j, x);
			s += strspn(s, " ,;");
		    }
		}
	    }
	} else {
	    for (i=0; i<r && !err; i++) {
		for (j=0; j<c && !err; j++) {
		    x = get_double(&s, &err);
		    if (!err) {
			gretl_matrix_set(M, i, j, x);	
			s += strspn(s, " ,;");
		    }
		}
	    }
	}	    
    }

    if (!err) {
	err = add_or_replace_user_matrix(M, name);
    }

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

int matrix_command (char *s, PRN *prn)
{
    char name[MNAMELEN];
    char word[9];
    int err = 0;

    if (!strncmp(s, "matrix ", 7)) s += 7;

    while (isspace(*s)) s++;

    if (!sscanf(s, "%31s", name)) {
	err = 1;
    } else {
	s += strlen(name);
	while (isspace(*s)) s++;
    }

    if (err) {
	return err;
    }

    if (*s == '=') {
	/* defining a matrix */
	err = create_matrix(name, s, prn);
    } else {
	*word = '\0';
	sscanf(s, "%8s", word);
	if (*word == '\0' || !strcmp(word, "print")) {
	    err = print_matrix_by_name(name, prn);
	} else {
	    err = 1;
	}
    } 

    return err;
}


