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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_string_table.h"
#include "gretl_private.h"

typedef struct _col_table col_table;

struct _col_table {
    int idx;
    int n_strs;
    char **strs;
};

struct _gretl_string_table {
    int n_cols;
    col_table **cols;
};

static col_table *col_table_new (int colnum)
{
    col_table *ct;

    ct = malloc(sizeof *ct);
    if (ct == NULL) return NULL;

    ct->strs = NULL;
    ct->n_strs = 0;
    ct->idx = colnum;

    return ct;
}

gretl_string_table *gretl_string_table_new (void)
{
    gretl_string_table *st;

    st = malloc(sizeof *st);
    if (st == NULL) return NULL;

    st->cols = NULL;
    st->n_cols = 0;

    return st;
}

static int col_table_get_index (const col_table *ct, const char *s)
{
    int i;

    for (i=0; i<ct->n_strs; i++) {
	if (!strcmp(s, ct->strs[i])) return i;
    }

    return -1;
}

static int 
col_table_add_string (col_table *ct, const char *s)
{
    char **strs;
    int n = ct->n_strs + 1;

    strs = realloc(ct->strs, n * sizeof *strs);
    if (strs == NULL) return -1;

    ct->strs = strs;
    strs[n-1] = malloc(strlen(s) + 1);
    if (strs[n-1] == NULL) return -1;

    strcpy(strs[n-1], s);
    ct->n_strs += 1;

    return n - 1;
}

static col_table *
gretl_string_table_add_column (gretl_string_table *st, int colnum)
{
    col_table **cols;
    int n = st->n_cols + 1;

    cols = realloc(st->cols, n * sizeof *cols);
    if (cols == NULL) return NULL;

    st->cols = cols;
    cols[n-1] = col_table_new(colnum);
    if (cols[n-1] == NULL) return NULL;

    st->n_cols += 1;

    return cols[n-1];
}

int 
gretl_string_table_index (gretl_string_table *st, const char *s, int col,
			  int addcol, PRN *prn)
{
    col_table *ct = NULL;
    int i, idx = -1;

    if (st == NULL) return idx;

    for (i=0; i<st->n_cols; i++) {
	if ((st->cols[i])->idx == col) {
	    ct = st->cols[i];
	    break;
	}
    }

    if (ct != NULL) {
	/* there's a table for this column already */
	idx = col_table_get_index(ct, s);
    } else if (addcol) {
	/* no table for this column yet: start one now */
	ct = gretl_string_table_add_column(st, col);
	if (ct != NULL) {
	    pprintf(prn, M_("variable %d: translating from strings to code numbers\n"), 
		    col);
	}
    }

    if (idx < 0 && ct != NULL) {
	idx = col_table_add_string(ct, s);
    }
	    
    return idx;
}

static void col_table_destroy (col_table *ct)
{
    int i;

    if (ct == NULL) return;

    for (i=0; i<ct->n_strs; i++) {
	free(ct->strs[i]);
    }
    free(ct->strs);
    free(ct);
}

void gretl_string_table_destroy (gretl_string_table *st)
{
    int i;

    if (st == NULL) return;

    for (i=0; i<st->n_cols; i++) {
	col_table_destroy(st->cols[i]);
    }
    free(st->cols);
    free(st);
}

int gretl_string_table_print (gretl_string_table *st, PATHS *ppaths, 
			       PRN *prn)
{
    int i, j;
    const col_table *ct;
    char stname[MAXLEN];
    FILE *fp;
    int err = 0;

    if (st == NULL) return 1;

    strcpy(stname, "string_table.txt");
    path_append(stname, ppaths->userdir);
    fp = fopen(stname, "w");
    if (fp == NULL) err = 1;

    for (i=0; i<st->n_cols; i++) {
	ct = st->cols[i];
	if (!err) {
	    fprintf(fp, M_("String code table for variable %d:\n"), ct->idx);
	}
	pprintf(prn, M_("String code table for variable %d:\n"), ct->idx);
	for (j=0; j<ct->n_strs; j++) {
	    if (!err) {
		fprintf(fp, "%3d = '%s'\n", j, ct->strs[j]);
	    }
	    pprintf(prn, "%3d = '%s'\n", j, ct->strs[j]);
	}
    }

    if (fp != NULL) fclose(fp);

    gretl_string_table_destroy(st);

    return err;
}


