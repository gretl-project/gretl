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

typedef struct ufunc_ ufunc;

struct ufunc_ {
    char *name;
    int n_args;
    char **args;
    int n_lines;
    char **lines;
};

static int n_ufuns;
static ufunc **ufuns;

static ufunc *ufunc_new (void)
{
    ufunc *func = malloc(sizeof *func);

    if (func == NULL) return NULL;

    func->name = NULL;
    func->n_args = func_nlines = 0;
    func->args = func->lines = NULL;

    return func;
}

static ufunc *add_ufunc (void)
{
    int nf = n_ufuns;
    ufunc **myfuns;

    myfuns = realloc(ufuns, (nf + 1) * sizeof *myfuns);
    if (myfuns == NULL) {
	return NULL;
    }
    ufuns = myfuns;

    ufuns[nf] = ufunc_new();
    if (ufuns[nf] == NULL) {
	return NULL;
    }

    n_ufuns++;

    return ufuns[nf];
}

static void free_ufunc (ufunc *fun)
{
    int i;

    for (i=0; i<fun->n_args; i++) {
	free(fun->args[i]);
    }
    free(fun->args);

    for (i=0; i<fun->n_lines; i++) {
	free(fun->lines[i]);
    }
    free(fun->lines);

    free(fun->name);
    free(fun);
}

static int delete_ufunc_by_name (const char *fname)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    free_ufunc(ufuns[i]);
	    ufuns[i] = NULL;
	    return 0;
	}
    }

    return 1;
}

static int check_func_name (const char *fname)
{
    int i;

    if (!isalpha((unsigned char) *fname)) {
	strcpy(gretl_errmsg, "function names must start with a letter");
	return 1;
    }

    if (gretl_command_number(fname)) {
	sprintf(gretl_errmsg, "'%s' is the name of a gretl command",
		fname);
	return 1;
    }

    for (i=0; i<n_ufuns; i++) {
	if (!strcmp(fname, (ufuns[i])->name)) {
	    sprintf(gretl_errmsg, "'%s': function is already defined",
		    fname);
	    return 1;
	}
    }

    return 0;
}

static int count_args (char *line)
{
    int na;

    return na;
}

int start_function_parse (const char *line)
{
    char *linecpy, *p, *fname;
    ufunc *fun = NULL;
    int n_args, i;

    if (strncmp(line, "function ", 9)) return 1;

    linecpy = gretl_strdup(line + 9);
    if (linecpy == NULL) {
	return E_ALLOC;
    }

    p = linecpy;

    namelen = strcspn(p, " \t");
    fname = malloc(namelen + 1);
    if (fname == NULL) {
	free(linecpy);
	return E_ALLOC;
    }
    *fname = '\0';
    strncat(fname, namelen);
    p += namelen + 1;

    if (check_func_name(fname)) {
	free(fname);
	free(linecpy);
	return E_ALLOC;
    }
    
    p += namelen + 1;
    n_args = count_args(p);

    /* do stuff */

    fun = add_ufunc();
    if (fun == NULL) {
	; /*error */
    }

    fun->name = fname;
    fun->n_args = n_args;
    
    /* etc */

    free(linecpy);

    return 0;
}
