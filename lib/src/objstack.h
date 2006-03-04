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
    GRETL_OBJ_ANY,
    GRETL_OBJ_EQN,
    GRETL_OBJ_SYS,
    GRETL_OBJ_VAR,
    GRETL_OBJ_DSET,
    GRETL_OBJ_INFO,
    GRETL_OBJ_STATS,
    GRETL_OBJ_CORR,
    GRETL_OBJ_SCRIPT,
    GRETL_OBJ_NOTES,
    GRETL_OBJ_MODTAB,
    GRETL_OBJ_GPAGE,
    GRETL_OBJ_GRAPH,
    GRETL_OBJ_PLOT,
    GRETL_OBJ_TEXT,
    GRETL_OBJ_UNKNOWN,
    GRETL_OBJ_MAX
} GretlObjType;

enum {
    OBJ_ACTION_NONE,
    OBJ_ACTION_INVALID,
    OBJ_ACTION_NULL,
    OBJ_ACTION_SHOW,
    OBJ_ACTION_FREE,
    OBJ_ACTION_SHOW_STAT,
    OBJ_ACTION_ADD,
    OBJ_ACTION_OMIT,
    OBJ_ACTION_IRF
};

void set_as_last_model (void *ptr, GretlObjType type);

void gretl_model_protect (MODEL *pmod);

void maybe_swap_into_last_model (MODEL *new, MODEL *old);

void *get_last_model (GretlObjType *type);

MODEL *get_model_by_name (const char *mname);

GRETL_VAR *get_VAR_by_name (const char *vname);

GRETL_VAR *get_VECM_by_name (const char *vname);

gretl_equation_system *get_equation_system_by_name (const char *sname);

void *gretl_get_object_by_name (const char *name);

int gretl_get_object_and_type (const char *name, void **pp, 
			       GretlObjType *type);

int gretl_stack_object (void *ptr, GretlObjType type);

int gretl_stack_object_as (void *ptr, GretlObjType type, const char *name);

void remove_model_from_stack (MODEL *pmod);

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

int gretl_object_rename (void *p, GretlObjType type, const char *oname);

int gretl_object_compose_name (void *p, GretlObjType type);

const char *gretl_object_get_name (void *p, GretlObjType type);

int parse_object_command (const char *s, char *name, char **cmd);

int match_object_command (const char *s, GretlObjType type);

void gretl_saved_objects_cleanup (void);

#endif

