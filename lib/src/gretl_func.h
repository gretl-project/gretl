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
    FN_NEEDS_DATA,  /* needs some (any) sort of dataset */
    FN_NEEDS_TS,    /* function requires time-series data */
    FN_NEEDS_QM,    /* function requires quarterly or monthly data */
    FN_NEEDS_PANEL, /* function requires panel data */
    FN_NODATA_OK    /* function does not require a dataset */
} DataReq;

typedef enum {
    UFUN_ROLE_NONE,
    UFUN_BUNDLE_PRINT,
    UFUN_BUNDLE_PLOT,
    UFUN_BUNDLE_TEST,
    UFUN_BUNDLE_FCAST,
    UFUN_BUNDLE_EXTRA,
    UFUN_GUI_MAIN,
    UFUN_GUI_PRECHECK,
    UFUN_PLOT_PRECHECK,
    UFUN_LIST_MAKER,
    UFUN_R_SETUP,
    UFUN_UI_MAKER,
    UFUN_ROLE_MAX
} UfunRole;

typedef enum {
    UFUN_PRIVATE   = 1 << 0, /* is private to a package */
    UFUN_NOPRINT   = 1 << 1, /* offers no printed output */
    UFUN_MENU_ONLY = 1 << 2, /* is GUI-only */
    UFUN_USES_SET  = 1 << 3, /* includes the "set" command */
    UFUN_HAS_FLOW  = 1 << 4, /* includes flow-control (ifs, loops) */
    UFUN_RECURSES  = 1 << 5  /* includes a call to itself */
} UfunAttrs;

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
				    r == GRETL_TYPE_STRINGS ||  \
				    r == GRETL_TYPE_MATRICES || \
				    r == GRETL_TYPE_BUNDLES ||	\
				    r == GRETL_TYPE_LISTS ||	\
				    r == GRETL_TYPE_ARRAYS ||	\
				    r == GRETL_TYPE_VOID || \
				    r == GRETL_TYPE_NUMERIC)

typedef struct ufunc_ ufunc;
typedef struct fncall_ fncall;
typedef struct fnpkg_ fnpkg;

int n_user_functions (void);

int n_free_functions (void);

ufunc *get_user_function_by_name (const char *name);

int is_user_function (const char *name);

const ufunc *get_user_function_by_index (int idx);

fncall *fncall_new (ufunc *fun, int preserve);

void fncall_destroy (void *ptr);

GretlType fncall_get_return_type (fncall *fc);

fncall *get_pkg_function_call (const char *funcname,
			       const char *pkgname,
			       const char *pkgpath);

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

int fn_param_automatic (const ufunc *fun, int i);

int fn_param_optional (const ufunc *fun, int i);

int fn_param_uses_xlist (const ufunc *fun, int i);

int fn_param_uses_mylist (const ufunc *fun, int i);

int user_func_get_return_type (const ufunc *fun);

fncall *user_func_get_fncall (ufunc *fun);

int user_func_is_noprint (const ufunc *fun);

int user_func_is_menu_only (const ufunc *fun);

const char *user_function_name_by_index (int i);

int user_function_index_by_name (const char *name, 
				 fnpkg *pkg);

void function_names_init (void);

const char *next_available_function_name (fnpkg *pkg,
					  int *idxp);

int gretl_compiling_function (void);

int gretl_compiling_python (const char *line);

int gretl_function_depth (void);

int gretl_function_recursing (void);

int function_is_executing (const char *funcname);

void current_function_info (char const **funcname,
			    char const **pkgname);

fnpkg *get_active_function_package (gretlopt opt);

fnpkg *gretl_function_get_package (const ufunc *fun);

int gretl_start_compiling_function (const char *line,
				    const DATASET *dset,
				    PRN *prn);

int gretl_function_append_line (ExecState *s);

int gretl_is_public_user_function (const char *name);

int gretl_function_exec_full (fncall *call, int rtype,
			      DATASET *dset, void *ret,
			      char **descrip, series_table **stab,
			      PRN *prn);

int gretl_function_exec (fncall *call, int rtype, DATASET *dset,
			 void *ret, PRN *prn);

int set_function_should_return (const char *line);

int current_function_size (void);

char *gretl_func_get_arg_name (const char *argvar, int *err);

int object_is_const (const char *name, int vnum);

int object_is_function_arg (const char *name);

void allow_full_data_access (int s);

int series_is_accessible_in_function (int ID, const DATASET *dset);

const char *series_get_list_parent (int ID);

void sample_range_get_extrema (const DATASET *dset, int *t1, int *t2);

void extend_function_sample_range (int addobs);

int function_return_type_from_string (const char *s);

int gretl_function_print_code (ufunc *u, int tabwidth, PRN *prn);

char **gretl_function_retrieve_code (ufunc *u, int *nlines);

int print_function_package_sample (const char *fname, int tabwidth,
				   PRN *prn);

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

const char *function_package_get_string (fnpkg *pkg,
					 const char *id);

int function_package_set_data_files (fnpkg *pkg, char **S, int n);

char **function_package_get_data_files (fnpkg *pkg, int *n);

int function_package_set_depends (fnpkg *pkg, char **S, int n);

char **function_package_get_depends (fnpkg *pkg, int *n);

const char *function_package_get_name (fnpkg *pkg);

double function_package_get_version (fnpkg *pkg);

int function_package_write_file (fnpkg *pkg);

int create_and_write_function_package (const char *fname, 
				       gretlopt opt,
				       PRN *prn);

int function_set_package_role (const char *name, fnpkg *pkg,
			       const char *attr, PRN *prn);

int function_ok_for_package_role (const char *name,
				  int role);

const char *package_role_get_key (int flag);

int check_function_needs (const DATASET *dset, DataReq dreq,
			  int minver);

int package_version_ok (int minver, char *reqstr);

int write_loaded_functions_file (const char *fname, int mpicall);

int read_session_functions_file (const char *fname);

fnpkg *get_function_package_by_name (const char *pkgname);

fnpkg *get_function_package_by_filename (const char *fname, int *err);

const char *get_function_package_path_by_name (const char *pkgname);

char **package_peek_dependencies (const char *fname, int *ndeps);

int grab_package_sample (const char *pkgname, char **pscript);

int include_gfn (const char *fname, gretlopt opt, PRN *prn);

int function_package_is_loaded (const char *fname,
				const char **version);

int gfn_is_loaded (const char *gfnname);

void function_package_unload_by_filename (const char *fname);

void function_package_unload_full_by_filename (const char *fname);

int function_package_unload_full (const char *pkgname);

int delete_function_package (const char *gfnname);

int uninstall_function_package (const char *package, gretlopt opt,
				PRN *prn);

int print_function_package_info (const char *fname, int gui_mode,
				 PRN *prn);

int bundle_function_package_info (const char *fname, gretl_bundle *b);

int print_function_package_code (const char *fname, int tabwidth,
				 PRN *prn);

int print_function_package_help (const char *fname, char **pbuf,
                                 PRN *prn);

ufunc *get_function_from_package (const char *funname, fnpkg *pkg);

int get_function_file_header (const char *fname, char **pdesc,
			      char **pver, char **pdate,
			      char **pauthor, int *pdfdoc);

double gfn_file_get_version (const char *fname);

int user_function_help (const char *fnname, gretlopt opt, PRN *prn);

int function_package_has_PDF_doc (fnpkg *pkg, char **pdfname);

int function_package_has_gui_help (fnpkg *pkg);

void function_package_set_editor (fnpkg *pkg, void *editor);

void *function_package_get_editor (fnpkg *pkg);

int package_has_menu_attachment (const char *fname,
				 char **pkgname,
				 char **attach,
				 char **label);

int package_needs_zipping (const char *fname,
			   int *pdfdoc,
			   char ***datafiles,
			   int *n_files);

void gretl_functions_cleanup (void);

int push_function_arg (fncall *fc, const char *name,
		       void *uvar, GretlType type,
		       void *value);

int push_anon_function_arg (fncall *fc, GretlType type,
			    void *value);

int set_anon_function_arg (fncall *fc, int i, GretlType type,
			   void *value);

int push_function_args (fncall *fc, ...);

void adjust_indent (const char *line, int *this_indent,
		    int *next_indent);

void normalize_hansl (const char *buf, int tabwidth,
		      PRN *prn);

#endif /* GRETL_FUNC_H */
