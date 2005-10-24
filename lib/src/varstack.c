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
#include "johansen.h"
#include "varstack.h"

static GRETL_VAR **var_stack;
static int n_vars;

static GRETL_VAR *
real_get_VAR_by_name (const char *vname, int *vnum)
{
    int i;

    if (vname == NULL) {
	return NULL;
    }

    for (i=0; i<n_vars; i++) {
	if (!strcmp(vname, var_stack[i]->name)) {
	    if (vnum != NULL) {
		*vnum = i;
	    }
	    return var_stack[i];
	}
    }

    return NULL;
}

GRETL_VAR *get_VAR_by_name (const char *vname)
{
    return real_get_VAR_by_name(vname, NULL);
}

GRETL_VAR *get_VECM_by_name (const char *vname)
{
    GRETL_VAR *var = real_get_VAR_by_name(vname, NULL);

    if (var != NULL && var->ci == VECM) {
	return var;
    } else {
	return NULL;
    }
}

static int stack_VAR (GRETL_VAR *var, const char *vname, PRN *prn)
{
    GRETL_VAR *orig;
    int vnum;

    if (var == NULL) {
	return 1;
    }

    if (vname != NULL) {
	var->name = gretl_strdup(vname);
    }

    if (var->name == NULL) {
	return 1;
    }

    orig = real_get_VAR_by_name(var->name, &vnum);

    if (orig != NULL) {
	/* replace existing VAR of same name */
	gretl_VAR_free(orig);
	var_stack[vnum] = var;
	pprintf(prn, "Replaced VAR '%s'\n", var->name);
    } else {
	GRETL_VAR **vstack;

	vstack = realloc(var_stack, (n_vars + 1) * sizeof *vstack);
	if (vstack == NULL) {
	    return E_ALLOC;
	}
	var_stack = vstack;
	var_stack[n_vars++] = var;
	pprintf(prn, "Added VAR '%s'\n", var->name);
    }

    return 0;
}

int stack_VAR_as (GRETL_VAR *var, const char *vname)
{
    return stack_VAR(var, vname, NULL);
}

int maybe_stack_var (GRETL_VAR *var, const CMD *cmd)
{
    const char *vname = gretl_cmd_get_name(cmd);
    int ret = 0;

    if (*vname) {
	ret = stack_VAR(var, vname, NULL);
    } else {
	gretl_VAR_free(var);
    }

    return ret;
}

void gretl_VARs_cleanup (void)
{
    int i;

    for (i=0; i<n_vars; i++) {
	gretl_VAR_free(var_stack[i]);
    }

    free(var_stack);
    var_stack = NULL;
    n_vars = 0;
}
