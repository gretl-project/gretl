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

#ifndef OBJSTACK_H
#define OBJSTACK_H

/* Note the following enumerarion may be added to, but the
   existing entries up to GRETL_OBJ_SCALARS must not be
   moved.
*/

typedef enum {
    GRETL_OBJ_NULL,
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
    GRETL_OBJ_MATRIX,
    GRETL_OBJ_SCALARS,
    GRETL_OBJ_ANY,
    GRETL_OBJ_MAX
} GretlObjType;

typedef enum {
    IN_GUI_SESSION = 1 << 0,
    IN_NAMED_STACK = 1 << 1,
    IN_MODEL_TABLE = 1 << 2,
    IS_LAST_MODEL  = 1 << 3,
    IN_GRAPH_PAGE  = 1 << 4
} SavedObjectFlags;

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

int gretl_model_protect (MODEL *pmod);

void *get_last_model (GretlObjType *type);

void *get_genr_model (GretlObjType *type);

MODEL *get_model_by_name (const char *mname);

MODEL *get_model_by_ID (int ID);

GRETL_VAR *get_VAR_by_name (const char *vname);

GRETL_VAR *get_VECM_by_name (const char *vname);

equation_system *get_equation_system_by_name (const char *sname);

void *gretl_get_object_by_name (const char *name);

int gretl_get_object_and_type (const char *name, void **pp, 
			       GretlObjType *type);

GretlObjType gretl_model_get_type_and_ci (const char *name,
					  int *ci);

int object_is_on_stack (const void *ptr);

int gretl_stack_object (void *ptr, GretlObjType type);

int gretl_stack_object_as (void *ptr, GretlObjType type, const char *name);

void gretl_object_remove_from_stack (void *ptr, GretlObjType type);

void remove_model_from_stack_on_exit (MODEL *pmod);

int maybe_stack_model (MODEL *pmod, const CMD *cmd, PRN *prn);

int maybe_stack_var (GRETL_VAR *var, const CMD *cmd);

void gretl_object_ref (void *ptr, GretlObjType type);

void gretl_object_unref (void *ptr, GretlObjType type);

double saved_object_get_scalar (const char *oname, int idx,
				double ***pZ, DATAINFO *pdinfo,
				int *err);

int saved_object_print_scalar (const char *oname, const char *key, PRN *prn);

double saved_object_get_scalar_element (const char *oname, const char *key,
					const DATAINFO *pdinfo, int *err);

double *saved_object_get_series (const char *oname, int idx, 
				 const DATAINFO *pdinfo, 
				 int *err);

gretl_matrix *
saved_object_get_matrix (const char *oname, int idx, int *err);

int *saved_object_get_list (const char *oname, int idx, int *err);

int gretl_object_rename (void *p, GretlObjType type, const char *oname);

int gretl_object_compose_name (void *p, GretlObjType type);

char *gretl_object_get_name (void *p, GretlObjType type);

int parse_object_command (const char *s, char *name, char **cmd);

int match_object_command (const char *s, GretlObjType type);

int last_model_test_ok (int ci, gretlopt opt, const DATAINFO *pdinfo, 
			PRN *prn);

int last_model_test_uhat (double ***pZ, DATAINFO *pdinfo, PRN *prn);

void set_genr_model (MODEL *pmod);

void unset_genr_model (void);

int highest_numbered_var_in_saved_object (const DATAINFO *pdinfo);

int check_variable_deletion_list (int *list, const DATAINFO *pdinfo);

void gretl_saved_objects_cleanup (void);

#endif

