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
#include "var.h"
#include "johansen.h"
#include "system.h"
#include "objstack.h"

typedef struct stacker_ stacker;

struct stacker_ {
    int type;
    void *ptr;
};

static stacker *obj_stack;
static int n_obj;

static void gretl_saved_object_free (stacker *s)
{
    if (s->type == EQUATION) {
	free_model(s->ptr);
    } else if (s->type == VAR) {
	gretl_VAR_free(s->ptr);
    } else if (s->type == SYSTEM) {
	gretl_equation_system_destroy(s->ptr);
    }
}

static const char *saved_object_get_name (stacker *s)
{
    if (s->type == EQUATION) {
	MODEL *pmod = s->ptr;
	return pmod->name;
    } else if (s->type == VAR) {
	GRETL_VAR *var = s->ptr;
	return var->name;
    } else if (s->type == SYSTEM) {
	gretl_equation_system *sys = s->ptr;
	return gretl_system_get_name(sys);
    }

    return NULL;
}

static void *get_object_by_name (const char *oname, int type, int *onum)
{
    const char *test;
    int i;

    if (oname == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	if (type != obj_stack[i].type) {
	    continue;
	}
	test = saved_object_get_name(&obj_stack[i]);
	if (!strcmp(oname, test)) {
	    if (onum != NULL) {
		*onum = i;
	    }
	    return obj_stack[i].ptr;
	}
    }

    return NULL;
}

GRETL_VAR *get_VAR_by_name (const char *vname)
{
    return get_object_by_name(vname, VAR, NULL);
}

GRETL_VAR *get_VECM_by_name (const char *vname)
{
    GRETL_VAR *var = get_object_by_name(vname, VAR, NULL);

    if (var != NULL && var->ci == VECM) {
	return var;
    } else {
	return NULL;
    }
}

gretl_equation_system *get_equation_system_by_name (const char *sname)
{
    return get_object_by_name(sname, SYSTEM, NULL);
}

static int gretl_object_set_name (void *p, int type, const char *oname)
{
    int err = 0;

    if (type == EQUATION) {
	gretl_model_set_name((MODEL *) p, oname);
    } else if (type == VAR) {
	gretl_VAR_set_name((GRETL_VAR *) p, oname);
    } else if (type == SYSTEM) {
	gretl_system_set_name((gretl_equation_system *) p, oname);
    } else {
	err = 1;
    }

    return err;
}

static int stack_object (void *p, int type, const char *oname, PRN *prn)
{
    stacker *orig;
    int onum, err = 0;

    if (p == NULL) {
	return 1;
    }

    if (oname != NULL) {
	err = gretl_object_set_name(p, type, oname);
	if (err) {
	    return err;
	}
    }

    orig = get_object_by_name(oname, type, &onum);

    if (orig != NULL) {
	/* replace existing object of same name */
	gretl_saved_object_free(orig);
	obj_stack[onum].ptr = p;
	obj_stack[onum].type = type;
	pprintf(prn, "Replaced object '%s'\n", oname);
    } else {
	stacker *ostack;

	ostack = realloc(obj_stack, (n_obj + 1) * sizeof *ostack);
	if (ostack == NULL) {
	    return E_ALLOC;
	}
	obj_stack = ostack;
	obj_stack[n_obj].ptr = p;
	obj_stack[n_obj].type = type;
	n_obj++;
	pprintf(prn, "Added object '%s'\n", oname);
    }

    return 0;
}

int stack_VAR_as (GRETL_VAR *var, const char *vname)
{
    return stack_object(var, VAR, vname, NULL);
}

int stack_system_as (gretl_equation_system *sys, const char *sname)
{
    return stack_object(sys, SYSTEM, sname, NULL);
}

int stack_system (gretl_equation_system *sys, PRN *prn)
{
    return stack_object(sys, SYSTEM, NULL, prn);
}

int maybe_stack_var (GRETL_VAR *var, const CMD *cmd)
{
    const char *vname = gretl_cmd_get_name(cmd);
    int ret = 0;

    if (*vname) {
	ret = stack_object(var, VAR, vname, NULL);
    } else {
	gretl_VAR_free(var);
    }

    return ret;
}

double saved_object_get_value (const char *oname, const char *valname)
{
    stacker *smatch = NULL;
    const char *test;
    double ret = NADBL;
    int i;

    for (i=0; i<n_obj; i++) {
	test = saved_object_get_name(&obj_stack[i]);
	if (!strcmp(oname, test)) {
	    smatch = &obj_stack[i];
	    break;
	}
    }

    if (smatch != NULL) {
	fprintf(stderr, "'%s': matched object %p, type %d\n", 
		oname, smatch->ptr, smatch->type);
    }

    return ret;
}

void gretl_saved_objects_cleanup (void)
{
    int i;

    for (i=0; i<n_obj; i++) {
	gretl_saved_object_free(&obj_stack[i]);
    }

    free(obj_stack);
    obj_stack = NULL;
    n_obj = 0;
}
