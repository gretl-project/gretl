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
#include "forecast.h"
#include "bootstrap.h"
#include "libset.h"

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
static stacker genr_model;

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

#if ODEBUG

static int object_get_refcount (void *ptr, GretlObjType type)
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

#endif

#if ODEBUG > 1

static void print_ostack (void)
{
    fprintf(stderr, "\n*** object stack: n_obj = %d\n", n_obj);

    if (n_obj > 0) {
	void *p;
	int i, t;

	for (i=0; i<n_obj; i++) {
	    p = ostack[i].ptr;
	    t = ostack[i].type;
	    fprintf(stderr, " %d: %p (type=%d, refcount=%d)\n", i, p,
		    t, object_get_refcount(p, t));
	}
	fputc('\n', stderr);
    }
}

#endif

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

enum {
    UNSTACK_REMOVE,
    UNSTACK_DESTROY
};

static int unstack_replace;

static void gretl_object_unstack (void *ptr, int action)
{
    int i, pos = -1;

#if ODEBUG
    fprintf(stderr, "gretl_object_unstack: ptr=%p, action=%d\n", ptr, action);
#endif

    if (action == UNSTACK_DESTROY && ptr == last_model.ptr) {
#if ODEBUG
	fprintf(stderr, " %p is 'last_model'\n", ptr);
#endif
	/* avoid double-freeing */
	last_model.ptr = NULL;
	last_model.type = GRETL_OBJ_NULL;
    }

    if (unstack_replace) {
	return;
    }

    for (i=0; i<n_obj; i++) {
	if (ptr == ostack[i].ptr) {
	    pos = i;
	    break;
	}
    }

#if ODEBUG
    fprintf(stderr, " stack pos for %p = %d, n_obj = %d\n",
	    ptr, pos, n_obj);
#endif

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

    gretl_object_unstack(ptr, UNSTACK_DESTROY);

    if (type == GRETL_OBJ_EQN) {
	gretl_model_free(ptr);
    } else if (type == GRETL_OBJ_VAR) {
	gretl_VAR_free(ptr);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system_destroy(ptr);
    }
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

int gretl_model_unprotect (MODEL *pmod)
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
    int *rc = NULL;

    if (ptr == NULL) {
	/* no-op */
	return;
    }

    if (type == GRETL_OBJ_ANY) {
	type = get_stacked_type_by_data(ptr);
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	if (pmod != NULL) {
	    if (model_is_protected(pmod)) {
		return; /* note! */
	    }
	    rc = &pmod->refcount;
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	if (var != NULL) {
	    rc = &var->refcount;
	}
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) ptr;

	if (sys != NULL) {
	    rc = &sys->refcount;
	}
    } else {
	fprintf(stderr, "gretl_object_unref: %p: bad object type\n",
		ptr);
	return;
    }

#if ODEBUG
    fprintf(stderr, "gretl_object_unref: %p (incoming count = %d)\n",
	    ptr, *rc);
#endif

    if (rc != NULL) {
	*rc -= 1;
	if (*rc <= 0) {
	    gretl_object_destroy(ptr, type);
	}
    }

#if ODEBUG > 1
    print_ostack();
#endif
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
	fprintf(stderr, " unrefing old object at %p (type %d)\n",
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

    /* Note: ensure that a newly estimated model supercedes
       the "genr_model" set by the GUI (fncall.c), if present.
    */
    if (ptr != NULL && genr_model.ptr != NULL && genr_model.ptr != ptr) {
	unset_genr_model();
    }

#if ODEBUG
    if (last_model.ptr != NULL) {
	int rc = object_get_refcount(last_model.ptr, last_model.type);

	fprintf(stderr, " refcount on \"last_model\" = %d\n", rc);
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
 * get_last_model_type:
 *
 * Returns: the type indentifier for the last model estimated.
 */

GretlObjType get_last_model_type (void)
{
    return last_model.type;
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

void *get_model_object_and_type (const char *name,
				 GretlObjType *type)
{
    void *ret = NULL;
    GretlObjType ti;
    const char *test;
    int i;

    *type = GRETL_OBJ_NULL;

    if (name == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	ti = ostack[i].type;
	if (ti == GRETL_OBJ_EQN || ti == GRETL_OBJ_VAR ||
	    ti == GRETL_OBJ_SYS) {
	    test = gretl_object_get_name(ostack[i].ptr, ti);
	    if (!strcmp(name, test)) {
		ret = ostack[i].ptr;
		*type = ti;
		break;
	    }
	}
    }

    return ret;
}

void *gretl_get_object_and_type (const char *name,
				 GretlObjType *type)
{
    void *ret = NULL;
    const char *test;
    int i;

    *type = GRETL_OBJ_NULL;

    if (name == NULL) {
	return NULL;
    }

    for (i=0; i<n_obj; i++) {
	test = gretl_object_get_name(ostack[i].ptr, ostack[i].type);
	if (!strcmp(name, test)) {
	    ret = ostack[i].ptr;
	    *type = ostack[i].type;
	    break;
	}
    }

    return ret;
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
    char name[48];
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

int gretl_object_compose_unique_name (void *p, GretlObjType type)
{
    char name[48];
    int id, err = 0;

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	id = pmod->ID;
	sprintf(name, "%s %d", _("Model"), id);
	while (get_model_by_name(name) != NULL) {
	    sprintf(name, "%s %d", _("Model"), ++id);
	}
	gretl_model_set_name(pmod, name);
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (var->ci == VAR) {
	    char *vstr = _("VAR");

	    if (strlen(vstr) > 3) {
		vstr = "VAR";
	    }
	    id = ++n_vars;
	    sprintf(name, "%s %d", vstr, id);
	    while (get_VAR_by_name(name) != NULL) {
		sprintf(name, "%s %d", vstr, ++id);
	    }
	} else {
	    char *vstr = _("VECM");

	    if (strlen(vstr) > 4) {
		vstr = "VECM";
	    }
	    id = gretl_VECM_id(var);
	    sprintf(name, "%s %d", vstr, id);
	    while (get_VECM_by_name(name) != NULL) {
		sprintf(name, "%s %d", vstr, ++id);
	    }
	}
	gretl_VAR_set_name(var, name);
    } else if (type == GRETL_OBJ_SYS) {
	equation_system *sys = (equation_system *) p;

	id = ++n_sys;
	sprintf(name, "%s %d", _("System"), id);
	while (get_equation_system_by_name(name) != NULL) {
	    sprintf(name, "%s %d", _("System"), ++id);
	}
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
#if ODEBUG
    fprintf(stderr, "gretl_object_remove_from_stack\n");
#endif
    gretl_object_unstack(ptr, UNSTACK_REMOVE);
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
#if ODEBUG
	fprintf(stderr, " real_stack_object: done, no-op\n");
#endif
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
	fprintf(stderr, "  replacing at %p (onum = %d, ptr = %p)\n",
		orig, onum, orig->ptr);
#endif
	unstack_replace = 1;
	gretl_object_unref(orig->ptr, orig->type);
	unstack_replace = 0;
	ostack[onum].ptr = p;
	ostack[onum].type = type;
	pprintf(prn, _("Replaced object '%s'\n"), name);
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
	pprintf(prn, _("Added object '%s'\n"), name);
    }

#if ODEBUG
    fprintf(stderr, " real_stack_object, on exit: '%s', type=%d, ptr=%p (n_obj=%d)\n",
	    name, type, (void *) p, n_obj);
#endif

#if ODEBUG > 1
    print_ostack();
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

int maybe_stack_var (GRETL_VAR *var, CMD *cmd)
{
    int err = 0;

    if (var != NULL) {
	const char *name = gretl_cmd_get_savename(cmd);

	set_as_last_model(var, GRETL_OBJ_VAR);
	if (*name != '\0') {
	    err = real_stack_object(var, GRETL_OBJ_VAR, name, NULL);
	}
    }

    return err;
}

/* Called from interact.c, after sucessful estimation of a
   (single-equation) model.  We automatically put the model in place
   as the "last model" for reference purposes (e.g. in genr).  In
   addition, if the model has been assigned a "savename" (via
   something like "mymod <- ols 1 0 2"), we make a copy and add it to
   the stack of saved objects.  We need to do the copying so that
   if/when the original model pointer is reassigned to, via a further
   estimation command, we do not lose the saved (named) model content.

   Reference counting: the model that is passed in should have its
   refcount raised by 1 on account of set_as_last_model.  The refcount of
   the copy (if applicable) should be 1, since the only reference to
   this model pointer is the one on the stack of named objects.

   We return a pointer to the model as stacked, in case any further
   action has to be taken (e.g. in the GUI).
*/

MODEL *maybe_stack_model (MODEL *pmod, CMD *cmd, PRN *prn, int *err)
{
    const char *name = gretl_cmd_get_savename(cmd);
    gretlopt opt = gretl_cmd_get_opt(cmd);
    MODEL *smod = NULL;

    if (*name != '\0' || (opt & OPT_W)) {
	MODEL *cpy = gretl_model_copy(pmod);

	if (cpy == NULL) {
	    *err = E_ALLOC;
	} else if (*name != '\0') {
	    *err = real_stack_object(cpy, GRETL_OBJ_EQN, name, NULL);
	}
	if (!*err) {
	    set_as_last_model(cpy, GRETL_OBJ_EQN);
	    if (*name != '\0') {
		pprintf(prn, _("%s saved\n"), name);
	    }
	} else {
	    errmsg(*err, prn);
	}
	smod = cpy;
    } else {
	set_as_last_model(pmod, GRETL_OBJ_EQN);
	smod = pmod;
    }

    return smod;
}

#define INVALID_STAT -999.999

/* retrieve from an object some value that is stored on the object in
   the form of a simple scalar */

static double real_get_obj_scalar (void *p, GretlObjType type,
				   DATASET *dset, int idx)
{
    double x = INVALID_STAT;
    int err = 0;

    if (idx <= 0) {
	return x;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) p;

	x = gretl_model_get_scalar(pmod, idx, dset, &err);
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
	} else if (idx == M_DF) {
	    x = sys->df;
	} else if (idx == M_DIAGTEST) {
	    err = system_diag_test(sys, &x, NULL);
	} else if (idx == M_DIAGPVAL) {
	    err = system_diag_test(sys, NULL, &x);
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) p;

	if (idx == M_T) {
	    x = var->T;
	} else if (idx == M_DF) {
	    x = var->df;
	} else if (idx == M_NCOEFF) {
	    x = var->ncoeff;
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

#define list_carrying_type(t) (t == GRETL_OBJ_EQN || \
			       t == GRETL_OBJ_VAR || \
			       t == GRETL_OBJ_SYS)

static int *
real_get_obj_list (void *p, GretlObjType type, int idx, int *err)
{
    const MODEL *pmod = NULL;
    const GRETL_VAR *var = NULL;
    const equation_system *sys = NULL;
    const int *list = NULL;
    int *ret = NULL;

    if (idx <= 0 || p == NULL || !list_carrying_type(type)) {
	*err = E_DATA;
	return NULL;
    }

    if (type == GRETL_OBJ_EQN) {
	pmod = p;
    } else if (type == GRETL_OBJ_VAR) {
	var = p;
    } else {
	sys = p;
    }

    if (idx == M_XLIST) {
	if (type == GRETL_OBJ_EQN) {
	    /* the list is already a copy */
	    ret = gretl_model_get_x_list(pmod);
	} else {
	    if (type == GRETL_OBJ_VAR) {
		list = var->xlist;
	    } else {
		list = sys->xlist;
	    }
	    if (list == NULL) {
		*err = E_BADSTAT;
	    } else {
		ret = gretl_list_copy(list);
	    }
	}
    } else if (idx == M_YLIST) {
	if (type == GRETL_OBJ_EQN) {
	    ret = gretl_model_get_y_list(pmod);
	} else {
	    if (type == GRETL_OBJ_VAR) {
		list = gretl_VAR_get_endo_list(var);
	    } else {
		list = system_get_endog_vars(sys);
	    }
	    if (list == NULL) {
		*err = E_BADSTAT;
	    } else {
		ret = gretl_list_copy(list);
	    }
	}
    } else {
	*err = E_BADSTAT;
    }

    if (ret == NULL && !*err) {
	*err = E_ALLOC;
    }

    return ret;
}

static char *
real_get_obj_string (void *p, GretlObjType type, int idx,
		     const DATASET *dset, int *err)
{
    char *str = NULL;

    if (idx <= 0) {
	*err = 1;
	return str;
    }

    if (idx == M_COMMAND) {
	if (type == GRETL_OBJ_EQN) {
	    MODEL *pmod = (MODEL *) p;

	    str = gretl_strdup(gretl_command_word(pmod->ci));
	} else if (type == GRETL_OBJ_SYS) {
	    str = gretl_strdup(gretl_command_word(SYSTEM));
	} else if (type == GRETL_OBJ_VAR) {
	    GRETL_VAR *var = (GRETL_VAR *) p;

	    str = gretl_strdup(gretl_command_word(var->ci));
	}
    } else if (idx == M_DEPVAR && type == GRETL_OBJ_EQN) {
	const char *s = gretl_model_get_depvar_name((MODEL *) p,
						    dset);
	str = gretl_strdup(s);
    }

    if (str == NULL) {
	*err = E_BADSTAT;
    }

    return str;
}

void set_genr_model (void *ptr, GretlObjType type)
{
    genr_model.ptr = ptr;
    genr_model.type = type;
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

int get_genr_model_ID (void)
{
    if (genr_model.ptr != NULL && genr_model.type == GRETL_OBJ_EQN) {
	MODEL *pmod = genr_model.ptr;

	return pmod->ID;
    } else {
	return 0;
    }
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

    if (smatch != NULL && smatch->type == GRETL_OBJ_EQN) {
	MODEL *pmod = smatch->ptr;

	if (pmod == NULL || pmod->errcode) {
	    fprintf(stderr, "find_smatch: duff model data!\n");
	    smatch = NULL;
	}
    }

    return smatch;
}

GretlType saved_object_get_data_type (const char *name,
				      ModelDataIndex idx)
{
    stacker *smatch = find_smatch(name);
    GretlType ret = GRETL_TYPE_NONE;

    /* note: handles M_UHAT, M_YHAT, M_SIGMA */

    if (smatch != NULL) {
	if (smatch->type == GRETL_OBJ_EQN) {
	    if (idx == M_SIGMA) {
		ret = GRETL_TYPE_DOUBLE;
	    } else {
		MODEL *pmod = smatch->ptr;

		if (pmod->ci == BIPROBIT ||
		    gretl_is_between_model(pmod)) {
		    ret = GRETL_TYPE_MATRIX;
		} else {
		    ret = GRETL_TYPE_SERIES;
		}
	    }
	} else {
	    /* VAR, VECM, system */
	    ret = GRETL_TYPE_MATRIX;
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

char *saved_object_get_string (const char *oname, int idx,
			       const DATASET *dset, int *err)
{
    char *ret = NULL;
    stacker *smatch;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	ret = real_get_obj_string(smatch->ptr, smatch->type, idx,
				  dset, err);
    }

    return ret;
}

double saved_object_get_scalar (const char *oname, int idx,
				DATASET *dset, int *err)
{
    double ret = INVALID_STAT;
    stacker *smatch;

    smatch = find_smatch(oname);

    if (smatch != NULL) {
	ret = real_get_obj_scalar(smatch->ptr, smatch->type,
				  dset, idx);
    }

    if (ret == INVALID_STAT) {
	*err = E_BADSTAT;
    }

    return ret;
}

int saved_object_get_series (double *x, const char *oname,
			     int idx, const DATASET *dset)
{
    int err = 0;

    if (idx <= 0) {
	err = E_DATA;
    } else {
	stacker *smatch = find_smatch(oname);

	if (smatch == NULL || smatch->type != GRETL_OBJ_EQN) {
	    err = E_BADSTAT;
	} else {
	    MODEL *pmod = (MODEL *) smatch->ptr;

	    err = gretl_model_get_series(x, pmod, dset, idx);
	}
    }

    return err;
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

    if (idx == M_FCAST || idx == M_FCSE) {
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

static gretl_matrix *get_all_probs (MODEL *pmod, int idx,
				    const DATASET *dset,
				    int *err)
{
    if (pmod->ci == LOGIT &&
	gretl_model_get_int(pmod, "multinom")) {
	/* M_MNLPROBS or M_ALLPROBS */
	return mn_logit_probabilities(pmod, pmod->t1, pmod->t2,
				      dset, err);
    }
    if (idx == M_ALLPROBS &&
	(pmod->ci == LOGIT || pmod->ci == PROBIT) &&
	gretl_model_get_int(pmod, "ordered")) {
	return ordered_probabilities(pmod, pmod->yhat,
				     pmod->t1, pmod->t2,
				     dset, err);
    }

    *err = E_BADSTAT;
    return NULL;
}

/* similar to the above, but in this case the matrix to be
   retrieved is not pre-built, and requires (or may require)
   access to the current dataset for its building.
*/

gretl_matrix *
saved_object_build_matrix (const char *oname, int idx,
			   const DATASET *dset,
			   int *err)
{
    stacker *smatch = find_smatch(oname);
    gretl_matrix *M = NULL;

    if (smatch == NULL) {
	*err = E_DATA;
    } else if (idx == M_EC && smatch->type == GRETL_OBJ_VAR) {
	M = VECM_get_EC_matrix(smatch->ptr, dset, err);
    } else if (idx == M_VMA && smatch->type == GRETL_OBJ_VAR) {
	M = gretl_VAR_get_vma_matrix(smatch->ptr, dset, err);
    } else if (idx == M_FEVD && smatch->type == GRETL_OBJ_VAR) {
	M = gretl_VAR_get_full_FEVD_matrix(smatch->ptr, dset, err);
    } else if ((idx == M_MNLPROBS || idx == M_ALLPROBS) &&
	       (smatch->type == GRETL_OBJ_EQN)) {
	M = get_all_probs(smatch->ptr, idx, dset, err);
    } else {
	*err = E_BADSTAT;
    }

    return M;
}

void *saved_object_get_array (const char *oname, int idx,
			      const DATASET *dset,
			      int *err)
{
    stacker *smatch = find_smatch(oname);
    gretl_array *A = NULL;

    if (smatch == NULL) {
	*err = E_DATA;
    } else if (idx == M_PARNAMES && smatch->type == GRETL_OBJ_EQN) {
	A = gretl_model_get_param_names(smatch->ptr, dset, err);
    } else {
	*err = E_BADSTAT;
    }

    return A;
}

gretl_matrix *
last_model_get_irf_matrix (int targ, int shock, double alpha,
			   const DATASET *dset, int *err)
{
    stacker *smatch = find_smatch(NULL);
    gretl_matrix *M = NULL;

    if (smatch == NULL || smatch->type != GRETL_OBJ_VAR) {
	*err = E_BADSTAT;
    } else {
	M = gretl_VAR_get_impulse_response(smatch->ptr, targ, shock, 0,
					   alpha, dset, err);
    }

    return M;
}

gretl_matrix *last_model_get_boot_ci (int cnum,
				      const DATASET *dset,
				      int B,
				      double alpha,
				      int method,
				      int studentize,
				      int *err)
{
    stacker *smatch = find_smatch(NULL);
    gretl_matrix *ret = NULL;

    if (smatch == NULL || smatch->type != GRETL_OBJ_EQN) {
	*err = E_DATA;
    } else {
	MODEL *pmod = smatch->ptr;

	ret = bootstrap_ci_matrix(pmod, dset, cnum, B, alpha,
				  method, studentize, err);
    }

    return ret;
}

double last_model_get_boot_pval (int cnum,
				 const DATASET *dset,
				 int B,
				 int method,
				 int *err)
{
    stacker *smatch = find_smatch(NULL);
    double ret = NADBL;

    if (smatch == NULL || smatch->type != GRETL_OBJ_EQN) {
	*err = E_DATA;
    } else {
	MODEL *pmod = smatch->ptr;

	ret = bootstrap_pvalue(pmod, dset, cnum, B, method, err);
    }

    return ret;
}

void *last_model_get_data (const char *key, GretlType *type,
			   int *size, int *copied, int *err)
{
    stacker *smatch = find_smatch(NULL);
    void *ret = NULL;

    if (smatch == NULL || smatch->type != GRETL_OBJ_EQN) {
	*err = E_DATA;
    } else {
	const MODEL *pmod = smatch->ptr;
	size_t sz = 0;

	ret = gretl_model_get_data_full(pmod, key, type, copied, &sz);
	if (ret == NULL) {
	    *err = E_DATA;
	} else if (size != NULL) {
	    *size = sz;
	}
    }

    if (*err) {
	gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
    }

    return ret;
}

char *last_model_get_vcv_type (void)
{
    stacker *smatch = find_smatch(NULL);
    static char ret[16];

    *ret = '\0';

    if (smatch != NULL && smatch->type == GRETL_OBJ_EQN) {
	const MODEL *pmod = smatch->ptr;
	VCVInfo *vi;

	vi = gretl_model_get_data(pmod, "vcv_info");

	if (vi != NULL && vi->vmaj == VCV_ML) {
	    if (vi->vmin == ML_HESSIAN) {
		strcpy(ret, "Hessian");
	    } else if (vi->vmin == ML_OP) {
		strcpy(ret, "OPG");
	    } else if (vi->vmin == ML_QML) {
		strcpy(ret, "Sandwich");
	    }
	}
    }

    return (*ret != '\0')? ret : NULL;
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
    int len, quoted = 0;
    int err = 0;

    *name = 0;
    *cmd = 0;

    /* skip any leading whitespace */
    while (*s && isspace(*s)) {
	s++;
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

int match_object_command (const char *s)
{
    if (strcmp(s, "show") == 0) {
	return OBJ_ACTION_SHOW;
    } else if (strcmp(s, "free") == 0) {
	return OBJ_ACTION_FREE;
    } else {
	return OBJ_ACTION_INVALID;
    }
}

static void saved_object_free (stacker *s)
{
    if (s == NULL) {
	return;
    }

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

#define sys_modtest_opt_ok(o) (o & (OPT_A | OPT_H | OPT_N))

int last_model_test_ok (int ci, gretlopt opt, const DATASET *dset,
			PRN *prn)
{
    GretlObjType type;
    void *ptr;
    int err = 0;

    ptr = get_last_model(&type);
    if (ptr == NULL) {
	pputs(prn, _("Can't do this: no model has been estimated yet\n"));
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN) {
	MODEL *pmod = (MODEL *) ptr;

	if (pmod->errcode) {
	    err = E_DATA;
	} else if (!model_test_ok(ci, opt, pmod, dset)) {
	    err = E_NOTIMP;
	} else if (ci == FCAST) {
	    err = fcast_not_feasible(pmod, dset);
	} else if (model_sample_problem(pmod, dset)) {
	    pputs(prn, _("Can't do: the current data set is different from "
			 "the one on which\nthe reference model was estimated\n"));
	    err = E_DATA;
	}
    } else if (type == GRETL_OBJ_SYS) {
	err = E_NOTIMP;
	if (ci == RESTRICT || ci == FCAST) {
	    err = 0;
	} else if (ci == MODTEST && sys_modtest_opt_ok(opt)) {
	    err = 0;
	}
    } else if (type == GRETL_OBJ_VAR) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;
	int r = gretl_VECM_rank(var);

	err = E_NOTIMP;

	if (ci == RESTRICT) {
	    err = 0;
	} else if (ci == FCAST) {
	    err = 0;
	} else if (ci == MODTEST && sys_modtest_opt_ok(opt)) {
	    err = 0;
	} else if (ci == OMIT && r == 0 && !(opt & OPT_A)) {
	    /* we can handle a test on exogenous terms in a VAR */
	    err = 0;
	}
    }

    return err;
}

int last_model_test_uhat (DATASET *dset, gretlopt opt, PRN *prn)
{
    GretlObjType type;
    void *ptr;
    int err = 0;

    ptr = get_last_model(&type);
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN) {
	err = model_error_dist(ptr, dset, opt, prn);
    } else if (type == GRETL_OBJ_SYS) {
	err = system_normality_test(ptr, opt, prn);
    } else if (type == GRETL_OBJ_VAR) {
	err = gretl_VAR_normality_test(ptr, opt, prn);
    } else {
	err = E_DATA;
    }

    return err;
}

int highest_numbered_var_in_saved_object (const DATASET *dset)
{
    GretlObjType type;
    void *ptr, *lmp = NULL;
    int i, mvm, vmax = 0;

    for (i=-1; i<n_obj; i++) {
	if (i < 0) {
	    lmp = ptr = get_last_model(&type);
	} else {
	    ptr = ostack[i].ptr;
	    type = ostack[i].type;
	}
	if (ptr == NULL || (i >= 0 && ptr == lmp)) {
	    continue;
	}
	if (type == GRETL_OBJ_EQN) {
	    mvm = highest_numbered_var_in_model((MODEL *) ptr, dset);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_VAR) {
	    mvm = gretl_VAR_get_highest_variable((GRETL_VAR *) ptr);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	} else if (type == GRETL_OBJ_SYS) {
	    mvm = highest_numbered_var_in_system((equation_system *) ptr, dset);
	    if (mvm > vmax) {
		vmax = mvm;
	    }
	}
    }

    return vmax;
}

int check_variable_deletion_list (int *list, const DATASET *dset)
{
    int pruned = 0;
    int i, vsave;

    vsave = highest_numbered_var_in_saved_object(dset);

    for (i=list[0]; i>0; i--) {
	if (list[i] <= vsave) {
	    gretl_list_delete_at_pos(list, i);
	    pruned = 1;
	}
    }

    return pruned;
}

static GList *(*get_or_send_gui_models)(GList *);

void set_gui_model_list_callback (GList *(*callback)())
{
    get_or_send_gui_models = callback;
}

/* The following function is called from subsample.c when
   the user has called for a permanent subsampling (represented
   by @newmask) of the original dataset.

   If @ndropped is non-NULL that means we're just checking
   whether the proposed subsampling is inconsistent with
   preservation of any saved models -- and this argument gets
   filled out with the number of such models. This variant
   call is issued only in gui mode.

   Otherwise, for each saved model, we either mark it for
   deletion (via the GUI) or revise its sample information
   appropriately in light of the shrinkage of the dataset.
*/

int check_models_for_subsample (char *newmask, int *ndropped)
{
    GList *fromgui = NULL;
    GretlObjType type;
    MODEL *pmod;
    void *ptr;
    int moderr;
    int err = 0;

    if (get_or_send_gui_models != NULL) {
	fromgui = (*get_or_send_gui_models)(NULL);
    }

    if (ndropped != NULL) {
	/* gui-only precheck: how many models are problematic? */
	*ndropped = 0;
	while (fromgui != NULL) {
	    pmod = fromgui->data;
	    moderr = subsample_check_model(pmod, newmask);
	    if (moderr) {
		*ndropped += 1;
		err = E_CANCEL;
	    }
	    if (fromgui->next != NULL) {
		fromgui = fromgui->next;
	    } else {
		break;
	    }
	}
	fprintf(stderr, "gui-precheck: ndropped = %d\n", *ndropped);
    } else {
	/* delete or fix-up affected models */
	GList *togui = NULL;

	ptr = get_last_model(&type);
	if (ptr != NULL && type == GRETL_OBJ_EQN) {
	    pmod = ptr;
	    moderr = subsample_check_model(pmod, newmask);
	    if (moderr) {
		set_as_last_model(NULL, GRETL_OBJ_NULL);
	    } else if (fromgui == NULL || g_list_find(fromgui, ptr) == NULL) {
		/* don't do this twice on a given model */
		remove_model_subsample_info(pmod);
	    }
	}

	while (fromgui != NULL) {
	    pmod = fromgui->data;
	    fprintf(stderr, "finalizing model %p\n", fromgui->data);
	    moderr = subsample_check_model(pmod, newmask);
	    if (moderr) {
		fprintf(stderr, " error: adding to deletion list\n");
		togui = g_list_append(togui, fromgui->data);
	    } else {
		fprintf(stderr, " OK: fixing sample info\n");
		remove_model_subsample_info(pmod);
	    }
	    if (fromgui->next != NULL) {
		fromgui = fromgui->next;
	    } else {
		break;
	    }
	}

	if (togui != NULL && get_or_send_gui_models != NULL) {
	    (*get_or_send_gui_models)(togui);
	}
    }

    if (fromgui != NULL) {
	g_list_free(fromgui);
    }

    return err;
}

int n_stacked_models (void)
{
    GList *list = NULL;
    int n = 0;

    if (get_or_send_gui_models != NULL) {
	list = (*get_or_send_gui_models)(NULL);
	n = g_list_length(list);
	g_list_free(list);
    } else {
	GretlObjType type;
	void *ptr;
	int i;

	for (i=0; i<n_obj; i++) {
	    ptr = ostack[i].ptr;
	    type = ostack[i].type;
	    if (ptr != NULL && type == GRETL_OBJ_EQN) {
		n++;
	    }
	}
    }

    return n;
}

void gretl_saved_objects_cleanup (void)
{
    void *lmp = last_model.ptr;
    int lmt = last_model.type;
    int i;

#if ODEBUG
    fprintf(stderr, "gretl_saved_objects_cleanup, n_obj = %d\n", n_obj);
#endif

    for (i=0; i<n_obj; i++) {
	if (ostack[i].ptr == lmp) {
	    /* stacked object i also occupies the "last model" place;
	       so we drop the associated refcount and nullify the
	       last model pointer to guard against double-freeing
	    */
#if ODEBUG
	    fprintf(stderr, " ostack[%d] = last model, unrefing\n", i);
#endif
	    gretl_object_unref(lmp, lmt);
	    last_model.ptr = NULL;
	    last_model.type = 0;
	    lmp = NULL;
	}
#if ODEBUG
	fprintf(stderr, " calling saved_object_free on ostack[%d] (%p)\n",
		i, ostack[i].ptr);
#endif
	saved_object_free(&ostack[i]);
    }

    free(ostack);
    ostack = NULL;

    n_obj = 0;
    n_sys = 0;
    n_vars = 0;

    if (lmp != NULL) {
	if (lmt != GRETL_OBJ_EQN || !model_is_protected(lmp)) {
#if ODEBUG
	    fprintf(stderr, "gretl_saved_objects_cleanup:\n"
		    " calling gretl_object_destroy on last_model (%p)\n", lmp);
#endif
	    gretl_object_destroy(lmp, lmt);
	}
	last_model.ptr = NULL;
	last_model.type = 0;
    }
}
