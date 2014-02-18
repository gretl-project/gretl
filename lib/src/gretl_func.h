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

#ifndef GRETL_FUNC_H
#define GRETL_FUNC_H

#include "usermat.h"

typedef enum {
    FN_NEEDS_TS = 1, /* function requires time-series data */
    FN_NEEDS_QM,     /* function requires quarterly or monthly data */
    FN_NEEDS_PANEL,  /* function requires panel data */
    FN_NODATA_OK     /* function does not require a dataset */
} FuncDataReq;

typedef enum {
    UFUN_ROLE_NONE,
    UFUN_BUNDLE_PRINT,
    UFUN_BUNDLE_PLOT,
    UFUN_BUNDLE_TEST,
    UFUN_BUNDLE_FCAST,
    UFUN_BUNDLE_EXTRA,
    UFUN_GUI_MAIN,
    UFUN_GUI_PRECHECK,
    UFUN_ROLE_MAX
} UfunRole; 

#define NEEDS_TS    "needs-time-series-data"
#define NEEDS_QM    "needs-qm-data"
#define NEEDS_PANEL "needs-panel-data"
#define NO_DATA_OK  "no-data-ok"

#define FN_NAMELEN 32

#define ok_function_return_type(r) (r == GRETL_TYPE_DOUBLE || \
				    r == GRETL_TYPE_SERIES || \
				    r == GRETL_TYPE_MATRIX || \
				    r == GRETL_TYPE_LIST ||   \
				    r == GRETL_TYPE_STRING || \
				    r == GRETL_TYPE_BUNDLE || \
				    r == GRETL_TYPE_VOID)

typedef struct ufunc_ ufunc;
typedef struct fnpkg_ fnpkg;
typedef struct fnargs_ fnargs;

int n_user_functions (void);

int n_free_functions (void);

ufunc *get_user_function_by_name (const char *name);

const ufunc *get_user_function_by_index (int idx);

int fn_n_params (const ufunc *fun);

int fn_param_type (const ufunc *fun, int i);

const char *fn_param_name (const ufunc *fun, int i);

const char *fn_param_descrip (const ufunc *fun, int i);

const char **fn_param_value_labels (const ufunc *fun, int i, 
				    int *n);

int fn_param_has_default (const ufunc *fun, int i);

double fn_param_default (const ufunc *fun, int i);

double fn_param_minval (const ufunc *fun, int i);

double fn_param_maxval (const ufunc *fun, int i);

double fn_param_step (const ufunc *fun, int i);

int fn_param_optional (const ufunc *fun, int i);

int fn_param_uses_xlist (const ufunc *fun, int i);

int user_func_get_return_type (const ufunc *fun);

int user_func_is_noprint (const ufunc *fun);

const char *user_function_name_by_index (int i);

int user_function_index_by_name (const char *name, 
				 fnpkg *pkg);

int user_function_set_debug (const char *name, int debug);

void function_names_init (void);

const char *next_available_function_name (int *idxp);

int gretl_compiling_function (void);

int gretl_compiling_python (const char *line);

int gretl_function_depth (void);

int gretl_start_compiling_function (const char *line, PRN *prn);

int gretl_function_append_line (const char *line);

int gretl_is_public_user_function (const char *name);

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 DATASET *dset, void *ret, 
			 char **descrip, PRN *prn);

char *gretl_func_get_arg_name (const char *argvar, int *err);

int object_is_const (const char *name);

int object_is_function_arg (const char *name);

void sample_range_get_extrema (const DATASET *dset, int *t1, int *t2);

int function_return_type_from_string (const char *s);

int gretl_function_print_code (ufunc *u, PRN *prn);

void set_current_function_package (fnpkg *pkg);

fnpkg *function_package_new (const char *fname, 
			     char **pubnames, int n_pub,
			     char **privnames, int n_priv, 
			     int *err);

int function_package_connect_funcs (fnpkg *pkg, 
				    char **pubnames, int n_pub,
				    char **privnames, int n_priv);

int function_package_set_properties (fnpkg *pkg, ...);

int function_package_get_properties (fnpkg *pkg, ...);

const char *function_package_get_name (fnpkg *pkg);

int function_package_write_file (fnpkg *pkg);

int create_and_write_function_package (const char *fname, 
				       gretlopt opt,
				       PRN *prn);

int function_set_package_role (const char *name, fnpkg *pkg,
			       const char *attr,
			       gretlopt opt, PRN *prn);

const char *package_role_get_key (int flag);

int check_function_needs (const DATASET *dset, FuncDataReq dreq,
			  int minver);

int package_version_ok (int minver, char *reqstr);

int write_session_functions_file (const char *fname);

int read_session_functions_file (const char *fname);

fnpkg *get_function_package_by_name (const char *pkgname);

fnpkg *get_function_package_by_filename (const char *fname, int *err);

int load_function_package_by_filename (const char *fname, PRN *prn);

int function_package_is_loaded (const char *fname);

void function_package_unload_by_filename (const char *fname);

void function_package_unload_full_by_filename (const char *fname);

int print_function_package_info (const char *fname, PRN *prn);

int print_function_package_code (const char *fname, PRN *prn);

ufunc *get_function_from_package (const char *funname, fnpkg *pkg);

int get_function_file_header (const char *fname, char **pdesc, char **pver);

int user_function_help (const char *fnname, gretlopt opt, PRN *prn);

int function_package_has_PDF_doc (fnpkg *pkg, char **pdfname);

void gretl_functions_cleanup (void);

const char *gretl_arg_type_name (GretlType type);

GretlType gretl_type_from_string (const char *s);

fnargs *fn_args_new (int argc);

void fn_args_free (fnargs *args);

int push_fn_arg (fnargs *args, const char *name, GretlType type,
		 void *value);

void adjust_indent (const char *line, int *this_indent,
		    int *next_indent);

#endif /* GRETL_FUNC_H */
