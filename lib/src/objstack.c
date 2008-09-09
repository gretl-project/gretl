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
#include "monte_carlo.h"
#include "var.h"
#include "johansen.h"
#include "system.h"
#include "objstack.h"
#include "usermat.h"
#include "modelspec.h"
#include "forecast.h"

#define ODEBUG 0

typedef struct stacker_ stacker;

struct stacker_ {
    int type;
    void *ptr;
};

static stacker *ostack;
static int n_obj;
static int n_sys;
static int n_vars;

static stacker last_model;

static GretlObjType get_stacked_type_by_data (void *ptr)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (ostack[i].ptr == ptr) {
	    return ostack[i].type;
	}
    }

    return GRETL_OBJ_NULL;
}

/**
 * gretl_object_ref:
 * @ptr: pointer to gretl obejct (e.g. #MODEL).
 * @type: type of object.
 *
 * Augments the reference count for the object represented by
 * @ptr, of type @type.
 */

void gretl_object_ref (void *ptr, GretlObjType type)
{
    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	if (pmod != NULL) {
	    pmod->refcount += 1;
#if ODEBUG
	    fprintf(stderr, "gretl_object_ref: refcount on %p is now %d\n",
		    (void *) pmod, pmod->refcount);
#endif
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	if (var != NULL) {
	    var->refcount += 1;
	}
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) ptr;

	if (sys != NULL) {
	    sys->refcount += 1;
	}
    }
}

static void gretl_object_unstack (void *ptr)
{
    int i, pos = -1;

    for (i=0; i<n_obj; i++) {
	if (ptr == ostack[i].ptr) {
	    pos = i;
	    break;
	}
    }
    
    if (pos >= 0) {
	n_obj--;
	if (n_obj == 0) {
	    free(ostack);
	    ostack = NULL;
	} else {
	    stacker *new_stack;

	    for (i=pos; i<n_obj; i++) {
		ostack[i] = ostack[i+1];
	    }

	    new_stack = realloc(ostack, n_obj * sizeof *new_stack);
	    if (new_stack != NULL) {
		ostack = new_stack;
	    }
	}
    }
}

static void gretl_object_destroy (void *ptr, GretlObjType type)
{
#if ODEBUG
    fprintf(stderr, "gretl_object_destroy: ptr %p, type %d\n",
	    ptr, type);
#endif

    gretl_object_unstack(ptr);
    
    if (type == GRETL_OBJ_EQN) {
	gretl_model_free(ptr);
    } else if (type == GRETL_OBJ_VAR) {
	gretl_VAR_free(ptr);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system_destroy(ptr);
    }
}

static int gretl_object_get_refcount (void *ptr, GretlObjType type)
{
    int rc = -999;

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	if (pmod != NULL) {
	    rc = pmod->refcount;
	}
	rc = pmod->refcount;
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	if (var != NULL) {
	    rc = var->refcount;
	}
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) ptr;

	if (sys != NULL) {
	    rc = sys->refcount;
	}
    }

    return rc;
}    

/* The stuff below: Note that we can't "protect" a model (against
   deletion) simply by setting its refcount to some special value,
   when models are being reused, since the refcount will go to 0 every
   time the model is assigned to!  Hence we need to set up this
   "protected species" list.
*/

static MODEL **protected_models;
static int n_prot;

int gretl_model_protect (MODEL *pmod)
{
    MODEL **prmod;
    int err = 0;

    prmod = realloc(protected_models, (n_prot + 1) * sizeof *prmod);

    if (prmod == NULL) {
	fprintf(stderr, "gretl_model_protect: out of memory!\n");
	err = E_ALLOC;
    } else {
	protected_models = prmod;
	protected_models[n_prot++] = pmod;
    }

    return err;
}

static int gretl_model_unprotect (MODEL *pmod)
{
    MODEL **prmod;
    int match = 0;
    int i, j, err = 0;

    for (i=0; i<n_prot; i++) {
	if (protected_models[i] == pmod) {
	    match = 1;
	    for (j=i; j<n_prot-1; j++) {
		protected_models[j] = protected_models[j+1];
	    }
	    break;
	}
    }

    if (match) {
	if (n_prot == 1) {
	    free(protected_models);
	    protected_models = NULL;
	    n_prot = 0;
	} else {
	    prmod = realloc(protected_models, (n_prot - 1) * sizeof *prmod);
	    if (prmod == NULL) {
		fprintf(stderr, "gretl_model_unprotect: out of memory!\n");
		err = E_ALLOC;
	    } else {
		protected_models = prmod;
		n_prot--;
	    }
	}
    }

    return err;
}

static int model_is_protected (MODEL *pmod)
{
    int i, prot = 0;

    for (i=0; i<n_prot; i++) {
	if (pmod == protected_models[i]) {
	    prot = 1;
	    break;
	}
    }

    if (!prot) {
	prot = model_is_in_loop(pmod);
    }

#if ODEBUG
    if (prot) {
	fprintf(stderr, "model at %p is protected\n", (void *) pmod);
    }
#endif

    return prot;
}

/**
 * gretl_object_unref:
 * @ptr: pointer to gretl object (e.g. #MODEL).
 * @type: type of object.
 *
 * Decrements the reference count for the object represented by
 * @ptr, of type @type. When the count reaches zero the object
 * is destroyed.
 */

void gretl_object_unref (void *ptr, GretlObjType type)
{
    int rc = 999;

    if (ptr == NULL) {
	/* no-op */
	return;
    }

    if (type == GRETL_OBJ_ANY) {
	type = get_stacked_type_by_data(ptr);
    } 

#if ODEBUG
    fprintf(stderr, "gretl_object_unref: %p (incoming count = %d)\n", 
	    ptr, gretl_object_get_refcount(ptr, type));
#endif

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	if (pmod != NULL) {
	    if (model_is_protected(pmod)) {
		return; /* note */
	    }
	    pmod->refcount -= 1;
	    rc = pmod->refcount;
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	if (var != NULL) {
	    var->refcount -= 1;
	    rc = var->refcount;
	}
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) ptr;

	if (sys != NULL) {
	    sys->refcount -= 1;
	    rc = sys->refcount;
	}
    }

    if (rc <= 0) {
	gretl_object_destroy(ptr, type);
    }
}

/**
 * set_as_last_model:
 * @ptr: pointer to gretl object (e.g. #MODEL).
 * @type: type of object.
 *
 * Puts @ptr in place as the "last model" (which will be
 * accessed by default via accessors such as "$uhat").
 */

void set_as_last_model (void *ptr, GretlObjType type)
{
#if ODEBUG
    fprintf(stderr, "set_as_last_model: ptr=%p, type=%d; 'old' ptr=%p\n",
	    ptr, type, last_model.ptr);
#endif

    if (last_model.ptr != ptr && last_model.ptr != NULL) {
#if ODEBUG
	fprintf(stderr, " kicking old object at %p (type %d) off the stack\n",
		last_model.ptr, last_model.type);
#endif
	gretl_object_unref(last_model.ptr, last_model.type);
    }

    if (last_model.ptr != ptr || last_model.type != type) {
	last_model.ptr = ptr;
	last_model.type = type;
	if (ptr != NULL) {
	    gretl_object_ref(ptr, type);
	}
    }

#if ODEBUG
    if (last_model.ptr != NULL) {
	fprintf(stderr, " refcount on \"last_model\" = %d\n",
		gretl_object_get_refcount(last_model.ptr, last_model.type));
    }
#endif
}

/**
 * get_last_model:
 * @type: location to receive type of last model, or %NULL.
 *
 * Returns: pointer to the last model estimated.  Note that
 * this may be %NULL if no model has been estimated.
 */

void *get_last_model (GretlObjType *type)
{
    if (type != NULL) {
	*type = last_model.type;
    }

    return last_model.ptr;
}

/**
 * gretl_object_get_name:
 * @p: pointer to gretl object (e.g. #MODEL).
 * @type: type of object.
 *
 * Returns: the name of the object of type @type with
 * location @p, or %NULL if the object is not found.
 * The return value may be ovewritten (up to
 * MAXSAVENAME-1 characters), but must not be freed.
 */

char *gretl_object_get_name (void *p, GretlObjType type)
{
    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = p;

	if (pmod->name == NULL) {
	    pmod->name = calloc(MAXSAVENAME, 1);
	}
	return pmod->name;
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = p;

	if (var->name == NULL) {
	    var->name = calloc(MAXSAVENAME, 1);
	}
	return var->name;
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = p;

	if (sys->name == NULL) {
	    sys->name = calloc(MAXSAVENAME, 1);
	}
	return sys->name;
    }

    return NULL;
}

static stacker *
get_stacker_by_name (const char *oname, GretlObjType type, int *onum)
{
    stacker *s = NULL;
    const char *test;
    int i;

    if (oname == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	if (type == GRETL_OBJ_ANY || type == ostack[i].type) {
	    test = gretl_object_get_name(ostack[i].ptr, ostack[i].type);
	    if (!strcmp(oname, test)) {
		if (onum != NULL) {
		    *onum = i;
		}
		s = &ostack[i];
		break;
	    }
	}
    }

#if ODEBUG
    fprintf(stderr, "get_stacker_by_name: name='%s', type=%d: got s=%p\n",
	    oname, type, (void *) s);
#endif
    
    return s;
}

static void *
get_object_by_name (const char *oname, GretlObjType type, int *onum)
{
    void *ptr = NULL;
    const char *test;
    int i;

    if (oname == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	if (type == GRETL_OBJ_ANY || type == ostack[i].type) {
	    test = gretl_object_get_name(ostack[i].ptr, ostack[i].type);
	    if (!strcmp(oname, test)) {
		if (onum != NULL) {
		    *onum = i;
		}
		ptr = ostack[i].ptr;
		break;
	    }
	}
    }

#if ODEBUG
    fprintf(stderr, "get_object_by_name: name='%s', type=%d, ptr=%p\n",
	    oname, type, (void *) ptr);
    if (ptr == NULL) {
	if (n_obj == 0) {
	    fprintf(stderr, "(no currently saved objects)\n");
	} else {
	    fprintf(stderr, "failed: current objects:\n");
	    for (i=0; i<n_obj; i++) {
		fprintf(stderr, " %02d: '%s' (%p, type %d)\n", i,
			gretl_object_get_name(ostack[i].ptr, ostack[i].type),
			ostack[i].ptr, ostack[i].type);
	    }
	}
    }
#endif    

    return ptr;
}

void *gretl_get_object_by_name (const char *name)
{
    return get_object_by_name(name, GRETL_OBJ_ANY, NULL);
}

int 
gretl_get_object_and_type (const char *name, void **pp, GretlObjType *type)
{
    const char *test;
    int i, err = E_DATA;

    *pp = NULL;
    *type = 0;

    if (name == NULL) {
	return err;
    }

    for (i=0; i<n_obj; i++) {
	test = gretl_object_get_name(ostack[i].ptr, ostack[i].type);
	if (!strcmp(name, test)) {
	    *pp = ostack[i].ptr;
	    *type = ostack[i].type;
	    err = 0;
	    break;
	}
    }

    return err;
}

MODEL *get_model_by_ID (int ID) 
{
    MODEL *pmod;
    int i;

    for (i=0; i<n_obj; i++) {
	if (ostack[i].type == GRETL_OBJ_EQN) {
	    pmod = (MODEL *) ostack[i].ptr;
	    if (pmod->ID == ID) {
		return pmod;
	    }
	}
    }

    return NULL;
}

MODEL *get_model_by_name (const char *mname)
{
    return get_object_by_name(mname, GRETL_OBJ_EQN, NULL);
}

GRETL_VAR *get_VAR_by_name (const char *vname)
{
    return get_object_by_name(vname, GRETL_OBJ_VAR, NULL);
}

GRETL_VAR *get_VECM_by_name (const char *vname)
{
    GRETL_VAR *var = get_object_by_name(vname, GRETL_OBJ_VAR, NULL);

    if (var != NULL && var->ci == VECM) {
	return var;
    } else {
	return NULL;
    }
}

equation_system *get_equation_system_by_name (const char *sname)
{
    return get_object_by_name(sname, GRETL_OBJ_SYS, NULL);
}

int gretl_object_compose_name (void *p, GretlObjType type)
{
    char name[32];
    int err = 0;

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	sprintf(name, "%s %d", _("Model"), pmod->ID);
	gretl_model_set_name(pmod, name);
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (var->ci == VAR) {
	    sprintf(name, "%s %d", _("VAR"), ++n_vars);
	} else {
	    sprintf(name, "%s %d", _("VECM"), gretl_VECM_id(var));
	}
	gretl_VAR_set_name(var, name);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) p;

	sprintf(name, "%s %d", _("System"), ++n_sys);
	equation_system_set_name(sys, name);
    } else {
	err = 1;
    }

    return err;
}

int gretl_object_rename (void *p, GretlObjType type, const char *oname)
{
    int err = 0;

    if (type == GRETL_OBJ_EQN) {
	gretl_model_set_name((MODEL *) p, oname);
    } else if (type == GRETL_OBJ_VAR) {
	gretl_VAR_set_name((GRETL_VAR *) p, oname);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system_set_name((equation_system *) p, oname);
    } else {
	err = 1;
    }

    return err;
}

/* safety measure prior to freeing models on exit */

void remove_model_from_stack_on_exit (MODEL *pmod)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (pmod == ostack[i].ptr) {
	    ostack[i].ptr = NULL;
	    ostack[i].type = 0;
	    break;
	}
    }

    if (last_model.ptr == pmod) {
	last_model.ptr = NULL;
	last_model.type = 0;
    }

    gretl_model_unprotect(pmod);
}

void gretl_object_remove_from_stack (void *ptr, GretlObjType type)
{
    gretl_object_unstack(ptr);
    gretl_object_unref(ptr, type);
}

static int object_stack_index (const void *p)
{
    int i, ret = -1;

    for (i=0; i<n_obj; i++) {
	if (p == ostack[i].ptr) {
	    ret = i;
	    break;
	}
    }

#if ODEBUG
    if (ret >= 0) {
	fprintf(stderr, "object_on_stack: object at %p already stacked "
		"at pos %d\n", p, ret);
    }
#endif

    return ret;
}

int object_is_on_stack (const void *ptr)
{
    return (object_stack_index(ptr) >= 0);
}

static int 
real_stack_object (void *p, GretlObjType type, const char *name, PRN *prn)
{
    stacker *orig;
    int onum, err = 0;

#if ODEBUG
    fprintf(stderr, "real_stack_object: on entry, p=%p, type=%d, name='%s'\n", 
	    (void *) p, type, name);
#endif

    if (p == NULL) {
	return 1;
    }

    if (object_stack_index(p) >= 0) {
	return 0;
    }

    if (name == NULL || *name == '\0') {
	err = gretl_object_compose_name(p, type);
    } else {
	err = gretl_object_rename(p, type, name);
    }

    if (err) {
	return err;
    }

    orig = get_stacker_by_name(name, type, &onum);

    if (orig != NULL) {
	/* replace existing object of same name */
#if ODEBUG
	fprintf(stderr, "stack_object: replacing at %p\n", orig->ptr);
#endif
	gretl_object_unref(orig->ptr, orig->type);
	ostack[onum].ptr = p;
	ostack[onum].type = type;
	pprintf(prn, "Replaced object '%s'\n", name);
	gretl_object_ref(p, type);
    } else {
	stacker *tmp;

	tmp = realloc(ostack, (n_obj + 1) * sizeof *ostack);
	if (tmp == NULL) {
	    return E_ALLOC;
	}
	ostack = tmp;
	ostack[n_obj].ptr = p;
	ostack[n_obj].type = type;
	gretl_object_ref(p, type);
	n_obj++;
	pprintf(prn, "Added object '%s'\n", name);

    }

#if ODEBUG
    fprintf(stderr, "stack_object: name='%s', type=%d, ptr=%p (n_obj=%d)\n",
	    name, type, (void *) p, n_obj);
#endif  

    return 0;
}

int gretl_stack_object (void *ptr, GretlObjType type)
{
    char *name = gretl_object_get_name(ptr, type);

    return real_stack_object(ptr, type, name, NULL);
}

int gretl_stack_object_as (void *ptr, GretlObjType type, const char *name)
{
    return real_stack_object(ptr, type, name, NULL);
}

int maybe_stack_var (GRETL_VAR *var, const CMD *cmd)
{
    char vname[MAXSAVENAME];
    int ret = 0;

    if (var == NULL) {
	return 0;
    }

    gretl_cmd_get_savename(vname);

#if ODEBUG
    fprintf(stderr, "\nmaybe_stack_var: initial refcount = %d\n", var->refcount);
#endif

    set_as_last_model(var, GRETL_OBJ_VAR);

#if ODEBUG
    fprintf(stderr, "maybe_stack_var: set %p as last model, refcount = %d\n", 
	    (void *) var, var->refcount);
#endif

    if (*vname) {
	ret = real_stack_object(var, GRETL_OBJ_VAR, vname, NULL);
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

    set_as_last_model(pmod, GRETL_OBJ_EQN);

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
	err = real_stack_object(cpy, GRETL_OBJ_EQN, mname, NULL);
    }

    if (!err) {
	pprintf(prn, _("%s saved\n"), mname);
    }

    return err;
}

#define INVALID_STAT -999.999

/* retrieve from an object some value that is stored on the object in
   the form of a simple scalar */

static double real_get_obj_scalar (void *p, GretlObjType type, int idx)
{
    double x = INVALID_STAT;
    int err = 0;
    
    if (idx <= 0) {
	return x;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	x = gretl_model_get_scalar(pmod, idx, &err);
	if (err) {
	    x = INVALID_STAT;
	}
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) p;

	if (idx == M_T) {
	    x = sys->T;
	} else if (idx == M_LNL) {
	    x = sys->ll;
	} else if (idx == M_ESS) {
	    x = sys->ess;
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (idx == M_T) {
	    x = var->T;
	} else if (idx == M_LNL) {
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

static double *
real_get_obj_series (void *p, GretlObjType type, int idx,
		     const DATAINFO *pdinfo, int *err)
{
    double *x = NULL;

    if (idx <= 0) {
	*err = 1;
	return x;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	x = gretl_model_get_series(pmod, pdinfo, idx, err);
    } 

    return x;
}

/* find out what sort of object we're dealing with, and call
   the appropriate function to get the requested matrix 
*/

static gretl_matrix *
real_get_obj_matrix (void *p, GretlObjType type, int idx, int *err)
{
    gretl_matrix *M = NULL;
    
    if (idx <= 0) {
	*err = 1;
	return M;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	M = gretl_model_get_matrix(pmod, idx, err);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) p;

	M = equation_system_get_matrix(sys, idx, err);
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	M = gretl_VAR_get_matrix(var, idx, err);
    }

    return M;
}

static int *
real_get_obj_list (void *p, GretlObjType type, int idx, int *err)
{
    int *list = NULL;
    
    if (idx <= 0) {
	*err = 1;
	return list;
    }

    if (type == GRETL_OBJ_EQN && idx == M_XLIST) {
	MODEL *pmod = (MODEL *) p;

	list = gretl_model_get_x_list(pmod);
    } else {
	*err = E_BADSTAT;
    }

    return list;
}

static stacker genr_model;

void set_genr_model (MODEL *pmod)
{
    genr_model.type = GRETL_OBJ_EQN;
    genr_model.ptr = pmod;
}

void unset_genr_model (void)
{
    genr_model.type = GRETL_OBJ_NULL;
    genr_model.ptr = NULL;
}

void *get_genr_model (GretlObjType *type)
{
    if (type != NULL) {
	*type = genr_model.type;
    }

    return genr_model.ptr;
}

static stacker *find_smatch (const char *oname)
{
    stacker *smatch = NULL;

    if (oname == NULL || *oname == '\0') {
	if (genr_model.ptr != NULL) {
	    smatch = &genr_model;
	} else {
	    smatch = &last_model;
	}
    } else {
	const char *test;
	int i;
	
	for (i=0; i<n_obj; i++) {
	    test = gretl_object_get_name(ostack[i].ptr, ostack[i].type);
	    if (!strcmp(oname, test)) {
		smatch = &ostack[i];
		break;
	    }
	}
    }

    return smatch;
}

/* retrieve the type and command index of a named model, or of 
   the last model if @name is NULL */

GretlObjType gretl_model_get_type_and_ci (const char *name,
					  int *ci)
{
    stacker *smatch = find_smatch(name);
    GretlObjType ret = GRETL_OBJ_NULL;

    *ci = 0;

    if (smatch != NULL) {
	ret = smatch->type;
	if (ret == GRETL_OBJ_EQN) {
	    MODEL *pmod = smatch->ptr;

	    *ci = pmod->ci;
	}
    }

    return ret;
}

int *saved_object_get_list (const char *oname, int idx, int *err)
{
    int *ret = NULL;
    stacker *smatch;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	ret = real_get_obj_list(smatch->ptr, smatch->type, idx, err);
    }

    return ret;
}

double saved_object_get_scalar (const char *oname, int idx, int *err)
{
    double ret = INVALID_STAT;
    stacker *smatch;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	ret = real_get_obj_scalar(smatch->ptr, smatch->type, idx);
    }

    if (ret == INVALID_STAT) {
	*err = E_BADSTAT;
    }

    return ret;
}

double *saved_object_get_series (const char *oname, int idx,
				 const DATAINFO *pdinfo, 
				 int *err)
{
    double *x = NULL;
    stacker *smatch;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	x = real_get_obj_series(smatch->ptr, smatch->type, idx, 
				pdinfo, err);
    }

    if (x == NULL && !*err) {
	*err = E_BADSTAT;
    }

    return x;
}

/* starting point for getting a matrix from a saved model
   of one kind or another.  We start by looking for the right
   model: this is either a match for @oname, or if @oname
   is NULL or blank, the last model estimated.
*/

gretl_matrix *
saved_object_get_matrix (const char *oname, int idx, int *err)
{
    gretl_matrix *M = NULL;

    if (idx == M_FCAST || idx == M_FCERR) {
	M = get_forecast_matrix(idx, err);
    } else {
	stacker *smatch = find_smatch(oname);

	if (smatch != NULL) {
	    M = real_get_obj_matrix(smatch->ptr, smatch->type, idx, err);
	}
    }

    if (M == NULL && !*err) {
	*err = E_BADSTAT;
    }    

    return M;
}

static int namechar_spn_with_space (const char *s)
{
    const char *ok = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"0123456789_ ";
    int ret = 0;

    if (isalpha(*s)) {
	ret = strspn(s, ok);
    }

    return ret;
}

#define OPDEBUG 0

/* try to parse an "object-oriented" command, such as
   MyModel.free or "System 1".show */

int parse_object_command (const char *s, char *name, char **cmd)
{
    int len, start = 0;
    int quoted = 0;
    int err = 0;

    *name = 0;
    *cmd = 0;

    /* skip any leading whitespace */
    while (*s && isspace(*s)) {
	s++; start++;
    }

    /* skip an opening quote */
    if (*s == '"') {
	quoted = 1;
	s++;
    }

    if (quoted) {
	len = namechar_spn_with_space(s);
    } else {
	len = gretl_namechar_spn(s);
    }

    if (len == 0) {
	return 0;
    } 

    if (len > MAXSAVENAME - 1) {
	len = MAXSAVENAME - 1;
    }

    strncat(name, s, len);
    s += len;

    if (quoted && *s == '"') {
	s++;
    }

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

int match_object_command (const char *s, GretlObjType type)
{
    if (type == GRETL_OBJ_EQN) {
	if (*s == 0) return OBJ_ACTION_SHOW; /* default */
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SHOW;
	if (strncmp(s, "add", 3) == 0) return OBJ_ACTION_ADD;
	if (strncmp(s, "omit", 4) == 0) return OBJ_ACTION_OMIT;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE;
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    } 

    if (type == GRETL_OBJ_VAR) {
	if (*s == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "irf") == 0)  return OBJ_ACTION_IRF;
	if (strncmp(s, "omit", 4) == 0) return OBJ_ACTION_OMIT;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE; 
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    }

    if (type == GRETL_OBJ_SYS) {
	if (*s == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "free") == 0) return OBJ_ACTION_FREE; 
	if (*s == '$') return OBJ_ACTION_SHOW_STAT;
    } 

    if (type == GRETL_OBJ_GRAPH) {
	if (*s == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SHOW;
    } 

    if (type == GRETL_OBJ_TEXT) {
	if (*s == 0) return OBJ_ACTION_SHOW;
	if (strcmp(s, "show") == 0) return OBJ_ACTION_SHOW;
    }  

    return OBJ_ACTION_INVALID;
}

static void saved_object_free (stacker *s)
{
    if (s->type == GRETL_OBJ_EQN) {
	if (!model_is_protected(s->ptr)) {
	    gretl_model_free(s->ptr);
	}
    } else if (s->type == GRETL_OBJ_VAR) {
	gretl_VAR_free(s->ptr);
    } else if (s->type == GRETL_OBJ_SYS) {
	equation_system_destroy(s->ptr);
    }
}

int last_model_test_ok (int ci, gretlopt opt, const DATAINFO *pdinfo, 
			PRN *prn)
{
    GretlObjType type;
    void *ptr;
    int err = 0;

    ptr = get_last_model(&type);  
    if (ptr == NULL) {
	pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	return 1;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;
  
	if (!model_test_ok(ci, opt, pmod, pdinfo)) {
	    err = E_NOTIMP;
	}
	if (model_sample_problem(pmod, pdinfo)) {
	    pputs(prn, _("Can't do: the current data set is different from "
			 "the one on which\nthe reference model was estimated\n"));
	    err = 1;
	}
    } else if (type == GRETL_OBJ_SYS) {
	err = E_NOTIMP;
	if (ci == RESTRICT || ci == TESTUHAT || ci == FCAST) {
	    err = 0;
	} else if (ci == LMTEST && ((opt & OPT_A) || (opt & OPT_H))) {
	    err = 0;
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;
	int r = gretl_VECM_rank(var);

	err = E_NOTIMP;

	if (ci == RESTRICT && r > 0) {
	    err = 0;
	} else if (ci == TESTUHAT || ci == FCAST) {
	    err = 0;
	} else if (ci == LMTEST && ((opt & OPT_A) || (opt & OPT_H))) {
	    err = 0;
	} 
    }

    if (err == E_NOTIMP) {
	pputs(prn, _("Sorry, command not available for this estimator"));
	pputc(prn, '\n');
    }

    return err;
}

int last_model_test_uhat (double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    GretlObjType type;
    void *ptr;
    int err = 0;

    ptr = get_last_model(&type);  
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = ptr;

	if ((pmod->ci == LOGIT || pmod->ci == PROBIT) &&
	    gretl_model_get_int(pmod, "ordered")) {
	    err = E_NOTIMP;
	} else if (!exact_fit_check(pmod, prn)) {
	    err = model_error_dist(ptr, pZ, pdinfo, prn);
	}
    } else if (type == GRETL_OBJ_SYS) {
	err = system_normality_test(ptr, prn);
    } else if (type == GRETL_OBJ_VAR) {
	err = gretl_VAR_normality_test(ptr, prn);
    } else {
	err = E_DATA;
    }

    return err;
}

int highest_numbered_var_in_saved_object (const DATAINFO *pdinfo)
{
    GretlObjType type;
    void *ptr;
    int i, mvm, vmax = 0;

    for (i=-1; i<n_obj; i++) {
	if (i < 0) {
	    ptr = get_last_model(&type);
	} else {
	    ptr = ostack[i].ptr;
	    type = ostack[i].type;
	}
	if (ptr == NULL) {
	    continue;
	}
	if (type == GRETL_OBJ_EQN) {
	    mvm = highest_numbered_var_in_model((MODEL *) ptr, pdinfo);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_VAR) {
	    mvm = gretl_VAR_get_highest_variable((GRETL_VAR *) ptr);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_SYS) {
	    mvm = highest_numbered_var_in_system((equation_system *) ptr, pdinfo);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	}
    }

    return vmax;
}

int check_variable_deletion_list (int *list, const DATAINFO *pdinfo)
{
    int pruned = 0;
    int i, vsave;

    vsave = highest_numbered_var_in_saved_object(pdinfo);

    for (i=list[0]; i>0; i--) {
	if (list[i] <= vsave) {
	    gretl_list_delete_at_pos(list, i);
	    pruned = 1;
	}
    }

    return pruned;
}

void gretl_saved_objects_cleanup (void)
{
    int i;

    for (i=0; i<n_obj; i++) {
	if (ostack[i].ptr == last_model.ptr) {
	    if (gretl_object_get_refcount(ostack[i].ptr, ostack[i].type) == 1) {
		/* don't double-free! */
		last_model.ptr = NULL;
	    }
	}
#if ODEBUG
	fprintf(stderr, "gretl_saved_objects_cleanup:\n"
		" calling saved_object_free on ostack[%d]\n", i);
#endif
	saved_object_free(&ostack[i]);
    }

    free(ostack);
    ostack = NULL;

    n_obj = 0;
    n_sys = 0;
    n_vars = 0;

    if (last_model.ptr != NULL) {
	if (last_model.type != GRETL_OBJ_EQN || 
	    !model_is_protected(last_model.ptr)) {
#if ODEBUG
	    fprintf(stderr, "gretl_saved_objects_cleanup:\n"
		    " calling gretl_object_destroy on last_model\n");
#endif
	    gretl_object_destroy(last_model.ptr, last_model.type);
	}
	last_model.ptr = NULL;
	last_model.type = 0;
    }
}

