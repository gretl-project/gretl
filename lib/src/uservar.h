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

#ifndef USERVAR_H_
#define USERVAR_H_

typedef enum {
    UVAR_ADD = 1,
    UVAR_DELETE
} UvarAction;

typedef enum {
    UV_PRIVATE = 1 << 0,
    UV_SHELL   = 1 << 1,
    UV_MAIN    = 1 << 2,
    UV_NODECL  = 1 << 3,
    UV_NOREPL  = 1 << 4
} UVFlags;

typedef int (*USER_VAR_FUNC) (const char *, GretlType, int);

typedef struct user_var_ user_var;

int user_var_add (const char *name, GretlType type, void *value);

user_var *alt_user_var_add (const char *name, GretlType type,
			    void *value);

int user_var_add_or_replace (const char *name,
			     GretlType type,
			     void *value);

int user_var_delete (user_var *uvar);

int user_var_delete_by_name (const char *name, PRN *prn);

int gretl_is_user_var (const char *name);

const char *uservar_name_complete (const char *s);

GretlType user_var_get_type_by_name (const char *name);

GretlType user_var_get_specific_type (const char *name);

user_var *get_user_var_by_name (const char *name);

user_var *get_user_var_of_type_by_name (const char *name,
					GretlType type);

user_var *get_user_var_by_data (const void *data);

const char *user_var_get_name (user_var *uvar);

const char *user_var_get_name_by_data (const void *data);

void *user_var_get_value (user_var *uvar);

GretlType user_var_get_type (user_var *uvar);

double user_var_get_scalar_value (user_var *uvar);

int user_var_set_scalar_value (user_var *uvar, double x);

void *user_var_get_value_by_name (const char *name);

void *user_var_get_value_and_type (const char *name,
				   GretlType *type);

char *user_string_resize (const char *name, size_t len, int *err);

char *user_string_reset (const char *name, const char *repl,
			 int *err);

void *user_var_steal_value (user_var *uvar);

void *user_var_unstack_value (user_var *uvar);

int user_var_get_level (user_var *uvar);

int user_var_get_flags (user_var *uvar);

int user_var_set_flag (user_var *uvar, UVFlags flag);

int user_var_unset_flag (user_var *uvar, UVFlags flag);

void user_var_privatize_by_name (const char *name);

int user_var_set_name (user_var *uvar, const char *name);

int user_var_adjust_level (user_var *uvar, int adj);

int user_var_localize (const char *origname,
		       const char *localname,
		       GretlType type);

void set_previous_depth (int d);

void switch_uservar_hash (int level);

int copy_as_arg (const char *param_name, GretlType type, 
		 void *value);

int arg_add_as_shell (const char *name,
		      GretlType type,
		      void *value);

int *copy_list_as_arg (const char *param_name, int *list,
		       int *err);

int user_var_replace_value (user_var *uvar, void *value,
			    GretlType type);

int user_var_set_pointer (user_var *uvar, void *ptr);

int user_matrix_replace_matrix_by_name (const char *name, 
					gretl_matrix *m);

void destroy_user_vars (void);

int destroy_user_vars_at_level (int level);

int destroy_private_matrices (void);

int destroy_private_lists (void);

int destroy_private_uvars (void);

int delete_user_vars_of_type (GretlType type, PRN *prn);

int n_user_matrices (void);

int n_user_scalars (void);

int n_user_lists (void);

int n_user_bundles (void);

GList *user_var_names_for_type (GretlType type);

GList *user_var_list_for_type (GretlType type);

void set_user_var_callback (USER_VAR_FUNC);

int copy_matrix_as (const gretl_matrix *m, const char *newname,
		    int fnarg);

int private_matrix_add (gretl_matrix *M, const char *name);

double get_scalar_value_by_name (const char *name, int *err);

int gretl_is_scalar (const char *name);

double gretl_scalar_get_value (const char *name, int *err);

int gretl_scalar_set_value (const char *name, double val);

int gretl_scalar_add (const char *name, double val);

int gretl_scalar_add_mutable (const char *name, double val);

int gretl_scalar_convert_to_matrix (user_var *uvar);

int private_scalar_add (double val, const char *name);

int add_auxiliary_scalar (const char *name, double val);

void set_auxiliary_scalars (void);

void unset_auxiliary_scalars (void);

void set_scalar_edit_callback (void (*callback));

void destroy_private_scalars (void);

void print_scalars (PRN *prn);

void print_scalar_by_name (const char *name, PRN *prn);

char *get_string_by_name (const char *name);

char *copy_string_by_name (const char *name, int *err);

int gretl_is_string (const char *name);

int is_user_string (const char *name);

char *temp_name_for_bundle (void);

int max_varno_in_saved_lists (void);

int gretl_lists_revise (const int *dlist, int dmin);

void gretl_lists_cleanup (void);

int serialize_user_vars (const char *dirname);

int deserialize_user_vars (const char *dirname);

int print_user_var_by_name (const char *name,
			    const DATASET *dset,
			    gretlopt opt,
			    PRN *prn);

int list_user_vars_of_type (const DATASET *dset,
			    PRN *prn);

int leads_midas_list (int ID, const DATASET *dset,
		      char *listname);

int in_midas_list (int ID, const DATASET *dset,
		   char *listname);

const char *
get_listname_by_consecutive_content (int l0, int l1);

int gretl_is_list (const char *name);

int *get_list_by_name (const char *name);

int user_list_append (user_var *uvar, const int *add);

int user_list_subtract (user_var *uvar, int *sub,
			const DATASET *dset);

int user_list_replace (user_var *uvar, const int *src);

int remember_list (const int *list, const char *name,
		   PRN *prn);

#endif /* USERVAR_H_ */
