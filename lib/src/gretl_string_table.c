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
    col_table *ct = malloc(sizeof *ct);

    if (ct != NULL) {
	ct->strs = NULL;
	ct->n_strs = 0;
	ct->idx = colnum;
    }

    return ct;
}

gretl_string_table *gretl_string_table_new (void)
{
    gretl_string_table *st = malloc(sizeof *st);

    if (st != NULL) {
	st->cols = NULL;
	st->n_cols = 0;
    }

    return st;
}

gretl_string_table *string_table_new_from_cols_list (int *list)
{
    gretl_string_table *st;
    int ncols = list[0];
    int i, j;

    st = malloc(sizeof *st);
    if (st == NULL) return NULL;

    st->cols = malloc(ncols * sizeof *st->cols);
    if (st->cols == NULL) {
	free(st);
	st = NULL;
    } else {
	st->n_cols = ncols;
	for (i=0; i<ncols; i++) {
	    st->cols[i] = col_table_new(list[i+1]);
	    if (st->cols[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(st->cols[j]);
		}
		free(st->cols);
		free(st);
		st = NULL;
	    } 
	}
    }

    return st;
}

static int col_table_get_index (const col_table *ct, const char *s)
{
    int ret = -1;
    int i;

    for (i=0; i<ct->n_strs; i++) {
	if (!strcmp(s, ct->strs[i])) {
	    ret = i + 1;
	    break;
	}
    }

    return ret;
}

static int 
col_table_add_string (col_table *ct, const char *s)
{
    char **strs;
    int n = ct->n_strs + 1;
    int ret = n;

    strs = realloc(ct->strs, n * sizeof *strs);
    if (strs == NULL) {
	ret = -1;
    } else {
	ct->strs = strs;
	strs[n-1] = gretl_strdup(s);

	if (strs[n-1] == NULL) {
	    ret = -1;
	} else {
	    ct->n_strs += 1;
	}
    }

    return ret;
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

int gretl_string_table_print (gretl_string_table *st, DATAINFO *pdinfo,
			      const char *fname, PRN *prn)
{
    int i, j;
    const col_table *ct;
    const char *fshort;
    char stname[MAXLEN];
    FILE *fp;
    int err = 0;

    if (st == NULL) return 1;

    strcpy(stname, "string_table.txt");
    gretl_path_prepend(stname, gretl_user_dir());

    fp = gretl_fopen(stname, "w");
    if (fp == NULL) {
	err = E_FOPEN;
	goto bailout;
    }

    fshort = strrchr(fname, SLASH);
    if (fshort != NULL) {
	fprintf(fp, "%s\n\n", fshort + 1);
    } else {
	fprintf(fp, "%s\n\n", fname);
    }

    fputs(M_("One or more non-numeric variables were found.\n"
	     "Gretl cannot handle such variables directly, so they\n"
	     "have been given numeric codes as follows.\n\n"), fp);

    for (i=0; i<st->n_cols; i++) {
	ct = st->cols[i];
	if (!err) {
	    fprintf(fp, M_("String code table for variable %d (%s):\n"), 
		    ct->idx, pdinfo->varname[ct->idx]);
	} else {
	    pprintf(prn, M_("String code table for variable %d (%s):\n"), 
		    ct->idx, pdinfo->varname[ct->idx]);
	}
	for (j=0; j<ct->n_strs; j++) {
	    if (!err) {
		fprintf(fp, "%3d = '%s'\n", j+1, ct->strs[j]);
	    } else {
		pprintf(prn, "%3d = '%s'\n", j+1, ct->strs[j]);
	    }
	}
    }

    if (fp != NULL) {
	pprintf(prn, M_("String code table written to\n %s\n"), stname);
	fclose(fp);
	set_string_table_written();
    }

 bailout:

    gretl_string_table_destroy(st);

    return err;
}


