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

#ifdef  __cplusplus
extern "C" {
#endif

/* Note: the following enumeration may be added to, but the
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
    GRETL_OBJ_BUNDLE,
    GRETL_OBJ_SESSION,
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
    OBJ_ACTION_FREE
};

void set_as_last_model (void *ptr, GretlObjType type);

int gretl_model_protect (MODEL *pmod);

int gretl_model_unprotect (MODEL *pmod);

void *get_last_model (GretlObjType *type);

GretlObjType get_last_model_type (void);

void *get_genr_model (GretlObjType *type);

int get_genr_model_ID (void);

MODEL *get_model_by_name (const char *mname);

MODEL *get_model_by_ID (int ID);

GRETL_VAR *get_VAR_by_name (const char *vname);

GRETL_VAR *get_VECM_by_name (const char *vname);

equation_system *get_equation_system_by_name (const char *sname);

void *gretl_get_object_by_name (const char *name);

void *get_model_object_and_type (const char *name,
				 GretlObjType *type);

void *gretl_get_object_and_type (const char *name,
				 GretlObjType *type);

GretlType saved_object_get_data_type (const char *name,
				      ModelDataIndex idx);

int object_is_on_stack (const void *ptr);

int gretl_stack_object (void *ptr, GretlObjType type);

int gretl_stack_object_as (void *ptr, GretlObjType type, const char *name);

void gretl_object_remove_from_stack (void *ptr, GretlObjType type);

void remove_model_from_stack_on_exit (MODEL *pmod);

MODEL *maybe_stack_model (MODEL *pmod, CMD *cmd, PRN *prn, int *err);

int maybe_stack_var (GRETL_VAR *var, CMD *cmd);

void gretl_object_ref (void *ptr, GretlObjType type);

void gretl_object_unref (void *ptr, GretlObjType type);

double saved_object_get_scalar (const char *oname, int idx,
				DATASET *dset, int *err);

int saved_object_print_scalar (const char *oname, const char *key, PRN *prn);

int saved_object_get_series (double *x, const char *oname,
			     int idx, const DATASET *dset);

gretl_matrix *
saved_object_get_matrix (const char *oname, int idx, int *err);

gretl_matrix *
saved_object_build_matrix (const char *oname, int idx,
			   const DATASET *dset, int *err);

gretl_matrix *
last_model_get_irf_matrix (int targ, int shock, double alpha,
			   const DATASET *dset, int *err);

void *saved_object_get_array (const char *oname, int idx,
			      const DATASET *dset,
			      int *err);

gretl_matrix *last_model_get_boot_ci (int cnum, const DATASET *dset,
				      int B, double alpha,
				      int method, int studentize,
				      int *err);

double last_model_get_boot_pval (int cnum,
				 const DATASET *dset,
				 int B, int method,
				 int *err);

void *last_model_get_data (const char *key, GretlType *type,
			   int *size, int *copied, int *err);

char *last_model_get_vcv_type (void);

int *saved_object_get_list (const char *oname, int idx, int *err);

char *saved_object_get_string (const char *oname, int idx,
			       const DATASET *dset, int *err);

int gretl_object_rename (void *p, GretlObjType type, const char *oname);

int gretl_object_compose_name (void *p, GretlObjType type);

int gretl_object_compose_unique_name (void *p, GretlObjType type);

char *gretl_object_get_name (void *p, GretlObjType type);

int parse_object_command (const char *s, char *name, char **cmd);

int match_object_command (const char *s);

int last_model_test_ok (int ci, gretlopt opt, const DATASET *dset,
			PRN *prn);

int last_model_test_uhat (DATASET *dset, gretlopt opt, PRN *prn);

void set_genr_model (void *ptr, GretlObjType type);

void unset_genr_model (void);

int highest_numbered_var_in_saved_object (const DATASET *dset);

int check_variable_deletion_list (int *list, const DATASET *dset);

int check_models_for_subsample (char *newmask, int *ndropped);

int n_stacked_models (void);

void set_gui_model_list_callback (GList *(*callback)(GList *));

void gretl_saved_objects_cleanup (void);

#ifdef  __cplusplus
}
#endif

#endif
