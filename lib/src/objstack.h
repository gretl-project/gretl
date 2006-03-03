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

#ifndef OBJSTACK_H
#define OBJSTACK_H

#include "system.h"

typedef enum {
    GRETL_OBJ_NONE,
    GRETL_OBJ_EQN,
    GRETL_OBJ_SYS,
    GRETL_OBJ_VAR,
    GRETL_OBJ_MAX
} GretlObjType;

enum {
    OBJ_ACTION_NONE,
    OBJ_ACTION_INVALID,
    OBJ_ACTION_NULL,
    OBJ_ACTION_FREE,
    OBJ_ACTION_MODEL_SHOW,
    OBJ_ACTION_MODEL_FREE,
    OBJ_ACTION_VAR_SHOW,
    OBJ_ACTION_VAR_IRF,
    OBJ_ACTION_VAR_FREE,
    OBJ_ACTION_GRAPH_SHOW,
    OBJ_ACTION_GRAPH_FREE,
    OBJ_ACTION_TEXT_SHOW,
    OBJ_ACTION_TEXT_FREE,
    OBJ_ACTION_SYS_SHOW,
    OBJ_ACTION_SYS_FREE,
    OBJ_ACTION_SHOW_STAT,
    OBJ_ACTION_MODEL_ADD,
    OBJ_ACTION_MODEL_OMIT,
    OBJ_ACTION_VAR_OMIT
};

#define obj_action_free(a) (a == OBJ_ACTION_MODEL_FREE || \
                            a == OBJ_ACTION_VAR_FREE || \
                            a == OBJ_ACTION_SYS_FREE || \
                            a == OBJ_ACTION_FREE)

#define GRETL_OBJ_PROTECTED 666

void set_as_last_model (void *ptr, GretlObjType type);

void maybe_swap_into_last_model (MODEL *new, MODEL *old);

void gretl_model_protect (MODEL *pmod);

MODEL *get_model_by_name (const char *mname);

GRETL_VAR *get_VAR_by_name (const char *vname);

GRETL_VAR *get_VECM_by_name (const char *vname);

gretl_equation_system *get_equation_system_by_name (const char *sname);

void *gretl_get_object_by_name (const char *name);

int gretl_get_object_and_type (const char *name, void **pp, 
			       GretlObjType *type);

int stack_model (MODEL *pmod);

int stack_model_as (MODEL *pmod, const char *mname);

void remove_model_from_stack (MODEL *pmod);

int stack_system (gretl_equation_system *sys, PRN *prn);

int stack_system_as (gretl_equation_system *sys, const char *sname);

int stack_VAR (GRETL_VAR *var);

int stack_VAR_as (GRETL_VAR *var, const char *vname);

int maybe_stack_model (MODEL *pmod, const CMD *cmd, PRN *prn);

int maybe_stack_var (GRETL_VAR *var, const CMD *cmd);

void gretl_object_ref (void *ptr, GretlObjType type);

void gretl_object_unref (void *ptr, GretlObjType type);

double saved_object_get_scalar (const char *oname, const char *key, int *err);

int saved_object_print_scalar (const char *oname, const char *key, PRN *prn);

double saved_object_get_scalar_element (const char *oname, const char *key,
					const DATAINFO *pdinfo, int *err);

double *saved_object_get_series (const char *oname, const char *key, 
				 const DATAINFO *pdinfo, int *err);

gretl_matrix *
saved_object_get_matrix (const char *oname, const char *key,
			 const double **Z, const DATAINFO *pdinfo,
			 int *err);

void gretl_rename_saved_object (void *p, const char *name);

int parse_object_command (const char *s, char *name, char **cmd);

void gretl_saved_objects_cleanup (void);

#endif

