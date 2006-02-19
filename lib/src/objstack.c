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
#include "system_private.h"
#include "objstack.h"
#include "genrfuncs.h"
#include "usermat.h"

#define ODEBUG 0

typedef struct stacker_ stacker;

struct stacker_ {
    int type;
    void *ptr;
};

static stacker *obj_stack;
static int n_obj;
static int n_sys;
static int n_vars;

static stacker last_model;

typedef enum {
    OBJ_FREE_NORMAL,
    OBJ_FREE_FINAL
} FreeCode;

typedef enum {
    REFCOUNT_ONE,
    REFCOUNT_PLUS
} RefCountCode;

static void gretl_object_set_refcount (stacker *s, RefCountCode code)
{
    if (s->type == EQUATION) {
	MODEL *pmod = (MODEL *) s->ptr;

	if (pmod != NULL) {
	    if (code == REFCOUNT_ONE) {
		pmod->refcount = 1;
	    } else if (code == REFCOUNT_PLUS) {
		pmod->refcount += 1;
	    }
	}
    } else if (s->type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) s->ptr;

	if (var != NULL) {
	    if (code == REFCOUNT_ONE) {
		var->refcount = 1;
	    } else if (code == REFCOUNT_PLUS) {
		var->refcount += 1;
	    }
	}
    } else if (s->type == SYSTEM) {
	gretl_equation_system *sys = (gretl_equation_system *) s->ptr;

	if (sys != NULL) {
	    if (code == REFCOUNT_ONE) {
		sys->refcount = 1;
	    } else if (code == REFCOUNT_PLUS) {
		sys->refcount += 1;
	    }
	}
    }
}

void gretl_object_ref (void *ptr, int type)
{
    stacker s;

    s.ptr = ptr;
    s.type = type;

    gretl_object_set_refcount(&s, REFCOUNT_PLUS);
}

static int gretl_object_get_refcount (stacker *s)
{
    int rc = -1;

#if ODEBUG
    fprintf(stderr, "gretl_object_get_refcount: s = %p, s->ptr = %p\n", 
	    (void *) s, (void *) s->ptr);
#endif   

    if (s->ptr == NULL) {
	return 0;
    }

    if (s->type == EQUATION) {
	MODEL *pmod = (MODEL *) s->ptr;
	rc = pmod->refcount;
    } else if (s->type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) s->ptr;
	rc = var->refcount;
    } else if (s->type == SYSTEM) {
	gretl_equation_system *sys = (gretl_equation_system *) s->ptr;
        fprintf(stderr, "getting refcount on %p\n", (void *) s->ptr);
	fflush(stderr);
	rc = sys->refcount;
    }

    return rc;
}

static MODEL *protected_models[4] = {
    NULL, NULL, NULL, NULL
};

void gretl_model_protect (MODEL *pmod)
{
    int i;

    for (i=0; i<4; i++) {
	if (protected_models[i] == NULL) {
	    protected_models[i] = pmod;
#if ODEBUG
	    fprintf(stderr, "model at %p: setting as protected\n", (void *) pmod);
#endif
	    break;
	}
    }
}

static int model_is_protected (MODEL *pmod)
{
    int i;

    for (i=0; i<4; i++) {
	if (pmod == protected_models[i]) {
#if ODEBUG
	    fprintf(stderr, "model_is_protected: returning 1 for pmod = %p\n", 
		    (void *) pmod);
#endif
	    return 1;
	}
    }

    return 0;
}

/* returns 1 if actual freeing will take place (the other possibility
   is just a reduction in the object's refcount)
*/

static int saved_object_free (stacker *s, FreeCode code)
{
    int rc;

    if (code == OBJ_FREE_FINAL) {
	gretl_object_set_refcount(s, REFCOUNT_ONE);
	rc = 1;
    } else {
	rc = gretl_object_get_refcount(s);
    }

#if ODEBUG
    fprintf(stderr, "saved_object_free: ptr=%p, refcount=%d\n", 
	    (void *) s->ptr, rc);    
#endif

    if (s->type == EQUATION) {
	if (!model_is_protected(s->ptr)) {
	    gretl_model_free(s->ptr);
	}
    } else if (s->type == VAR) {
	gretl_VAR_free(s->ptr);
    } else if (s->type == SYSTEM) {
	gretl_equation_system_destroy(s->ptr);
    }

    if (rc <= 1) {
	/* we actually freed the object */
	s->ptr = NULL;
	s->type = 0;
    }

    return rc <= 1;
}

void *get_last_model (int *type)
{
    if (last_model.type == 0) {
	*type = 0;
	return NULL;
    } else {
	*type = last_model.type;
	return last_model.ptr;
    }
}

void set_as_last_model (void *ptr, int type)
{
#if ODEBUG
    fprintf(stderr, "set_as_last_model: ptr=%p, type=%d; old ptr=%p\n",
	    ptr, type, last_model.ptr);
#endif

    if (ptr != last_model.ptr) {
#if ODEBUG
	if (last_model.ptr) {
	    fprintf(stderr, " kicking old object at %p (type %d) off the stack\n",
		last_model.ptr, last_model.type);
	}
	fprintf(stderr, " calling saved_object_free\n");
#endif
	saved_object_free(&last_model, OBJ_FREE_NORMAL);
	last_model.ptr = ptr;
	last_model.type = type;
	if (type != EQUATION) {
	    gretl_object_ref(ptr, type);
	}
    }

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) ptr;
	pmod->refcount = 2;
    } 

#if ODEBUG
    fprintf(stderr, " refcount on \"last_model\" = %d\n",
	    gretl_object_get_refcount(&last_model));
#endif
}

void set_as_last_model_if_unnamed (void *ptr, int type)
{
    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) ptr;

#if ODEBUG
	fprintf(stderr, "set_as_last_model_if_unnamed: pmod->name = %s\n",
		(pmod->name == NULL)? "NULL" : pmod->name);
#endif
	if (pmod->name == NULL) {
	    set_as_last_model(ptr, type);
	}
    } else if (type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	if (var->name == NULL) {
	    set_as_last_model(ptr, type);
	}
    }
}

/* Response to swap_models() in gretl_model.c.  This is done (only)
   when we're swapping the content of the pointer models[0] with
   something else, so that models[0] will continue to point to
   the model most recently estimated.
*/

void maybe_swap_into_last_model (MODEL *new, MODEL *old)
{
#if ODEBUG
    fprintf(stderr, "\nmaybe_swap_into_last_model: new=%p, old=%p\n",
	    (void *) new, (void *) old);
#endif
    if (last_model.ptr == old) {
#if ODEBUG
	fprintf(stderr, " doing swap 1\n");
#endif
	last_model.ptr = new;
	if (new->refcount < 2) {
	    new->refcount = 2;
	}
    } else if (last_model.ptr == new) {
#if ODEBUG
	fprintf(stderr, " doing swap 2\n");
#endif
	last_model.ptr = old;
    } else {
	fprintf(stderr, " No swap done\n");
    }

#if ODEBUG
    fprintf(stderr, " refcount on last_model = %d\n",
	    gretl_object_get_refcount(&last_model));
#endif
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
	return sys->name;
    }

    return NULL;
}

static void *get_object_by_name (const char *oname, int type, int *onum)
{
    void *ptr = NULL;
    const char *test;
    int i;

    if (oname == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	if (type >= 0 && type != obj_stack[i].type) {
	    continue;
	}
	test = saved_object_get_name(&obj_stack[i]);
	if (!strcmp(oname, test)) {
	    if (onum != NULL) {
		*onum = i;
	    }
	    ptr = obj_stack[i].ptr;
	    break;
	}
    }

#if ODEBUG
    fprintf(stderr, "get_object_by_name: name='%s', type=%d, ptr=%p\n",
	    oname, type, (void *) ptr);
#endif    

    return ptr;
}

void *gretl_get_object_by_name (const char *name)
{
    return get_object_by_name(name, -1, NULL);
}

int gretl_get_object_and_type (const char *name, void **pp, int *type)
{
    const char *test;
    int i, err = E_DATA;

    *pp = NULL;
    *type = 0;

    if (name == NULL) {
	return err;
    }

    for (i=0; i<n_obj; i++) {
	test = saved_object_get_name(&obj_stack[i]);
	if (!strcmp(name, test)) {
	    *pp = obj_stack[i].ptr;
	    *type = obj_stack[i].type;
	    err = 0;
	    break;
	}
    }

    return err;
}

/* Called from session.c when the user chooses to delete a model that
   is stored in the session as an icon; also called in cli program in
   response to "Obj.free".  If the object in question is also
   referenced as last_model, it won't actually be destroyed though it
   will be removed from the list of named objects.
*/

void gretl_delete_saved_object (void *p)
{
    int i, pos = -1;

#if ODEBUG
    fprintf(stderr, "gretl_delete_saved_object: got %p\n", p);
#endif

    for (i=0; i<n_obj; i++) {
	if (p == obj_stack[i].ptr) {
	    pos = i;
	    break;
	}
    }
    
    if (pos >= 0) {
	stacker *new_stack;

	saved_object_free(&obj_stack[pos], OBJ_FREE_NORMAL);

	for (i=pos; i<n_obj-1; i++) {
	    obj_stack[i] = obj_stack[i+1];
	}

	n_obj--;

	new_stack = realloc(obj_stack, (n_obj - 1) * sizeof *new_stack);
	if (new_stack != NULL) {
	    obj_stack = new_stack;
	}
    }
}

MODEL *get_model_by_name (const char *mname)
{
    return get_object_by_name(mname, EQUATION, NULL);
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

static int gretl_object_auto_assign_name (void *p, int type)
{
    char name[32];
    int err = 0;

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) p;

	sprintf(name, "%s %d", _("Model"), pmod->ID);
	gretl_model_set_name(pmod, name);
    } else if (type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (var->ci == VAR) {
	    sprintf(name, "%s %d", _("VAR"), ++n_vars);
	} else {
	    sprintf(name, "%s %d", _("VECM"), gretl_VECM_id(var));
	}
	gretl_VAR_set_name(var, name);
    } else if (type == SYSTEM) {
	gretl_equation_system *sys = (gretl_equation_system *) p;

	sprintf(name, "%s %d", _("System"), ++n_sys);
	gretl_system_set_name(sys, name);
    } else {
	err = 1;
    }

    return err;
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

void gretl_rename_saved_object (void *p, const char *name)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (p == obj_stack[i].ptr) {
	    gretl_object_set_name(p, obj_stack[i].type, name);
	    break;
	}
    }
}

void remove_model_from_stack (MODEL *pmod)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (pmod == obj_stack[i].ptr) {
	    obj_stack[i].ptr = NULL;
	    obj_stack[i].type = 0;
	    break;
	}
    }

    if (last_model.ptr == pmod) {
	last_model.ptr = NULL;
	last_model.type = 0;
    }
}

static int object_on_stack (void *p)
{
    int i, ret = -1;

    for (i=0; i<n_obj; i++) {
	if (p == obj_stack[i].ptr) {
	    ret = i;
	    break;
	}
    }

#if ODEBUG
    if (ret >= 0) {
	fprintf(stderr, "object_on_stack: object at %p already stacked at pos %d\n",
		p, ret);
    }
#endif

    return ret;
}

static int stack_object (void *p, int type, const char *oname, PRN *prn)
{
    stacker *orig;
    int onum, err = 0;

    if (p == NULL) {
	return 1;
    }

    if (object_on_stack(p) >= 0) {
	return 0;
    }

    if (oname == NULL) {
	err = gretl_object_auto_assign_name(p, type);
    } else {
	err = gretl_object_set_name(p, type, oname);
    }

    if (err) return err;

    orig = get_object_by_name(oname, type, &onum);

    if (orig != NULL) {
	/* replace existing object of same name */
#if ODEBUG
	fprintf(stderr, "stack_object: calling saved_object_free\n");
#endif
	saved_object_free(orig, OBJ_FREE_NORMAL);
	obj_stack[onum].ptr = p;
	obj_stack[onum].type = type;
	pprintf(prn, "Replaced object '%s'\n", oname);
	gretl_object_set_refcount(&obj_stack[onum], REFCOUNT_PLUS);
    } else {
	stacker *ostack;

	ostack = realloc(obj_stack, (n_obj + 1) * sizeof *ostack);
	if (ostack == NULL) {
	    return E_ALLOC;
	}
	obj_stack = ostack;
	obj_stack[n_obj].ptr = p;
	obj_stack[n_obj].type = type;
	gretl_object_set_refcount(&obj_stack[n_obj], REFCOUNT_PLUS);
	n_obj++;
	pprintf(prn, "Added object '%s'\n", oname);

    }

#if ODEBUG
    fprintf(stderr, "stack_object: name='%s', type=%d, ptr=%p (n_obj=%d)\n",
	    oname, type, (void *) p, n_obj);
#endif  

    return 0;
}

int stack_model (MODEL *pmod)
{
    return stack_object(pmod, EQUATION, NULL, NULL);
}

int stack_model_as (MODEL *pmod, const char *mname)
{
    return stack_object(pmod, EQUATION, mname, NULL);
}

int stack_VAR (GRETL_VAR *var)
{
    return stack_object(var, VAR, NULL, NULL);
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
    char vname[MAXSAVENAME];
    int ret = 0;

    gretl_cmd_get_savename(vname);

#if ODEBUG
    fprintf(stderr, "\nmaybe_stack_var: initial refcount = %d\n", var->refcount);
#endif

    set_as_last_model(var, VAR);

#if ODEBUG
    fprintf(stderr, "maybe_stack_var: set %p as last model, refcount = %d\n", 
	    (void *) var, var->refcount);
#endif

    if (*vname) {
	ret = stack_object(var, VAR, vname, NULL);
    } 

    return ret;
}

/* Called in gretlcli.c, after sucessful estimation of a
   (single-equation) model.  We automatically put the model in place
   as the "last model" for reference purposes (e.g. in genr).  In
   addition, if the model has been assigned a "savename" (via
   something like "mymod <- ols 1 0 2"), we make a copy and add it to
   the stack of saved objects.  We need to do the copying so that
   if/when the original model pointer is reassigned to, via a further
   estimation command, we do not lose the saved (named) model content.

   Reference counting: the model that is passed in should have its
   refcount raised to 2 on account of set_as_last_model.  The refcount of
   the copy (if applicable) should be 1, since the only reference to
   this model pointer is the one on the stack of named objects.
*/

int maybe_stack_model (MODEL *pmod, const CMD *cmd, PRN *prn)
{
    char mname[MAXSAVENAME];
    MODEL *cpy = NULL;
    int err = 0;

    gretl_cmd_get_savename(mname);

#if ODEBUG
    fprintf(stderr, "\nmaybe_stack_model: initial refcount = %d\n", 
	    pmod->refcount);
#endif

    set_as_last_model(pmod, EQUATION);

#if ODEBUG
    fprintf(stderr, "maybe_stack_model: set %p as last model, refcount %d\n",
	    (void *) pmod, pmod->refcount);
#endif

    if (*mname == 0) {
	return 0;
    }

    cpy = gretl_model_copy(pmod);
    if (cpy == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = stack_model_as(cpy, mname);
    }

    if (!err) {
	pprintf(prn, _("%s saved\n"), mname);
    }

    return err;
}

#define INVALID_STAT -999.999

/* retrieve from an object some value that is stored on the object in
   the form of a simple scalar */

static double real_get_obj_scalar (void *p, int type, int idx)
{
    double x = INVALID_STAT;
    int err = 0;
    
    if (idx <= 0) {
	return x;
    }

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) p;

	x = gretl_model_get_scalar(pmod, idx, &err);
	if (err) {
	    x = INVALID_STAT;
	}
    } else if (type == SYSTEM) {
	gretl_equation_system *sys = (gretl_equation_system *) p;

	if (idx == M_LNL) {
	    x = sys->ll;
	} else if (idx == M_ESS) {
	    x = sys->ess;
	}
    } else if (type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (idx == M_LNL) {
	    x = var->ll;
	} else if (idx == M_AIC) {
	    x = var->AIC;
	} else if (idx == M_BIC) {
	    x = var->BIC;	
	} else if (idx == M_HQC) {
	    x = var->HQC;
	}
    }

    return x;
}

/* retrieve from an object some value that is stored on the object in
   the form of a scalar element of an array */

static double 
real_get_obj_scalar_element (void *p, int type, int idx, const char *key, 
			     const DATAINFO *pdinfo, int *err)
{
    double x = INVALID_STAT;
    
    if (idx <= 0) {
	*err = 1;
	return x;
    }

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) p;

	x = get_model_data_element(pmod, idx, key, pdinfo, err);
	if (*err) {
	    x = INVALID_STAT;
	}
    } 

    /* FIXME add accessors for systems, VARs */

    return x;
}

static double *real_get_obj_series (void *p, int type, int idx, 
				    const DATAINFO *pdinfo,
				    int *err)
{
    double *x = NULL;

    if (idx <= 0) {
	*err = 1;
	return x;
    }

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) p;

	x = gretl_model_get_series(pmod, pdinfo, idx, err);
    }

    return x;
}

static gretl_matrix *real_get_obj_matrix (void *p, int type, int idx, int *err)
{
    gretl_matrix *M = NULL;
    
    if (idx <= 0) {
	*err = 1;
	return M;
    }

    if (type == EQUATION) {
	MODEL *pmod = (MODEL *) p;

	M = gretl_model_get_matrix(pmod, idx, err);
    } else if (type == SYSTEM) {
	gretl_equation_system *sys = (gretl_equation_system *) p;

	M = gretl_equation_system_get_matrix(sys, idx, err);
    } else if (type == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	M = gretl_VAR_get_matrix(var, idx, err);
    }

    return M;
}

static stacker *find_smatch (const char *oname)
{
    stacker *smatch = NULL;

    if (oname == NULL || *oname == '\0') {
	smatch = &last_model;
    } else {
	const char *test;
	int i;
	
	for (i=0; i<n_obj; i++) {
	    test = saved_object_get_name(&obj_stack[i]);
	    if (!strcmp(oname, test)) {
		smatch = &obj_stack[i];
		break;
	    }
	}
    }

    return smatch;
}

double saved_object_get_scalar (const char *oname, const char *key,
				int *err)
{
    double ret = INVALID_STAT;
    stacker *smatch;
    int idx;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	idx = gretl_model_data_index(key);
	ret = real_get_obj_scalar(smatch->ptr, smatch->type, idx);
    }

    if (ret == INVALID_STAT) {
	*err = 1;
    }

    return ret;
}

double saved_object_get_scalar_element (const char *oname, const char *key,
					const DATAINFO *pdinfo, int *err)
{
    double ret = INVALID_STAT;
    stacker *smatch;
    int idx;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	idx = gretl_model_data_index(key);
	ret = real_get_obj_scalar_element(smatch->ptr, smatch->type, idx, 
					  key, pdinfo, err);
    }

    if (ret == INVALID_STAT && !*err) {
	*err = 1;
    }

    return ret;
}

double *saved_object_get_series (const char *oname, const char *key, 
				 const DATAINFO *pdinfo, int *err)
{
    double *x = NULL;
    stacker *smatch;
    int idx;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	idx = gretl_model_data_index(key);
	x = real_get_obj_series(smatch->ptr, smatch->type, idx, 
				pdinfo, err);
    }

    if (x == NULL && !*err) {
	*err = 1;
    }

    return x;
}

gretl_matrix *
saved_object_get_matrix (const char *oname, const char *key,
			 const double **Z, const DATAINFO *pdinfo,
			 int *err)
{
    gretl_matrix *M = NULL;
    const char *mspec = NULL;
    stacker *smatch;
    int idx;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	idx = gretl_model_data_index(key);
	M = real_get_obj_matrix(smatch->ptr, smatch->type, idx, err);
	mspec = strchr(key, '[');
    }

    if (M != NULL && mspec != NULL) {
	gretl_matrix *S = matrix_get_submatrix(M, mspec, Z, pdinfo, err);
	
	gretl_matrix_free(M);
	M = S;
    }

    if (M == NULL && !*err) {
	*err = 1;
    }    

    return M;
}

int 
saved_object_print_scalar (const char *oname, const char *key, PRN *prn)
{
    int err = 0;
    double val = saved_object_get_scalar(oname, key, &err);

    if (err) {
	pprintf(prn, _("%s: no data for '%s'\n"), oname, key);
    } else {
	pprintf(prn, "%s: %s = %.8g\n", oname, key + 1, val);
    }

    return err;
}

#define OPDEBUG 0

/* try to parse an "object-priented" command, such as
   MyModel.free or "System 1".show */

int parse_object_command (const char *s, char *name, char **cmd)
{
    int len, start = 0;
    int quoted = 0;
    const char *p;
    int err = 0;

    *name = 0;
    *cmd = 0;

    /* skip any leading whitespace */
    while (*s && isspace(*s)) {
	s++; start++;
    }

    /* skip an opening quote */
    if (*s == '"') {
	s++;
	quoted = 1;
    }

    p = s;

    /* crawl to end of (possibly quoted) "name" */
    len = 0;
    while (*s) {
	if (*s == '"') {
	    quoted = 0;
	    len--;
	}
	if (!quoted && (isspace(*s) || *s == '.')) break;
	s++; len++;
    }

    if (len == 0) {
	return 0;
    }    

    if (len > MAXSAVENAME - 1) {
	len = MAXSAVENAME - 1;
    }

    strncat(name, p, len);

    /* is an object command embedded? */
    if (*s == '.') {
	s++;
	if (*s && !isspace(*s)) {
	    *cmd = gretl_strdup(s);
	    if (*cmd == NULL) {
		err = 1;
	    }
	}
    }

#if OPDEBUG
    if (*cmd != NULL) {
	fprintf(stderr, "name='%s', cmd='%s'\n", name, *cmd);
    } else {
	fprintf(stderr, "name='%s'\n", name);
    }
#endif
    
    return err;
}

void gretl_saved_objects_cleanup (void)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (obj_stack[i].ptr == last_model.ptr) {
	    /* don't double-free! */
	    last_model.ptr = NULL;
	}
#if ODEBUG
	fprintf(stderr, "gretl_saved_objects_cleanup:\n"
		" calling saved_object_free on obj_stack[%d]\n", i);
#endif
	saved_object_free(&obj_stack[i], OBJ_FREE_FINAL);
    }

    free(obj_stack);
    obj_stack = NULL;

    n_obj = 0;
    n_sys = 0;
    n_vars = 0;

    if (last_model.ptr != NULL) {
#if ODEBUG
	fprintf(stderr, "gretl_saved_objects_cleanup:\n"
		" calling saved_object_free on last_model\n");
#endif
	saved_object_free(&last_model, OBJ_FREE_FINAL);
	last_model.ptr = NULL;
	last_model.type = 0;
    }
}
