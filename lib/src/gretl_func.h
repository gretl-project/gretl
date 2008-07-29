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

#define NEEDS_TS    "needs-time-series-data"
#define NEEDS_QM    "needs-qm-data"
#define NEEDS_PANEL "needs-panel-data"
#define NO_DATA_OK  "no-data-ok"

#define FN_NAMELEN 32

typedef struct ufunc_ ufunc;
typedef struct fnpkg_ fnpkg;
typedef struct fnargs_ fnargs;

int n_free_functions (void);

ufunc *get_user_function_by_name (const char *name);

const ufunc *get_user_function_by_index (int idx);

int fn_n_params (const ufunc *fun);

int fn_param_type (const ufunc *fun, int i);

const char *fn_param_name (const ufunc *fun, int i);

const char *fn_param_descrip (const ufunc *fun, int i);

double fn_param_default (const ufunc *fun, int i);

double fn_param_minval (const ufunc *fun, int i);

double fn_param_maxval (const ufunc *fun, int i);

int fn_param_optional (const ufunc *fun, int i);

int user_func_get_return_type (const ufunc *fun);

const char *user_function_name_by_index (int i);

int user_function_index_by_name (const char *name);

void function_names_init (void);

const char *next_free_function_name (void);

int gretl_compiling_function (void);

int gretl_function_depth (void);

int gretl_start_compiling_function (const char *line, PRN *prn);

int gretl_function_append_line (const char *line);

int gretl_is_public_user_function (const char *name);

int gretl_function_exec (ufunc *u, fnargs *args, int rtype,
			 double ***pZ, DATAINFO *pdinfo,
			 void *ret, char **descrip,
			 PRN *prn);

char *gretl_func_get_arg_name (const char *argvar, int *err);

int object_is_const (const char *name);

const char *get_funcerr_message (void);

int gretl_function_set_info (int i, const char *help);

int gretl_function_get_info (int i, const char *key, char const **value);

int gretl_function_print_code (int i, PRN *prn);

void gretl_function_set_private (int i, int priv);

int write_function_package (fnpkg *pkg,
			    const char *fname,
			    int pub, 
			    const int *privlist, 
			    const char *author,
			    const char *version,
			    const char *date,
			    const char *descrip,
			    char *sample,
			    FuncDataReq dreq,
			    int minver);

int function_package_get_info (const char *fname,
			       fnpkg **ppkg,
			       int *pub, 
			       int **privlist,
			       char **author,
			       char **version,
			       char **date,
			       char **descrip,
			       char **sample,
			       FuncDataReq *dreq,
			       int *minver);

int check_function_needs (const DATAINFO *pdinfo, FuncDataReq dreq,
			  int minver);

int write_user_function_file (const char *fname);

int function_package_is_loaded (const char *fname);

const char *function_package_description (const char *fname);

int read_session_functions_file (const char *fname);

int load_user_function_file (const char *fname);

int get_function_file_info (const char *fname, PRN *prn, char **pname);

int get_function_file_code (const char *fname, PRN *prn, char **pname);

char *get_function_file_header (const char *fname, char **pver,
				int *err);

int update_function_from_script (const char *fname, int idx);

int user_function_help (const char *fnname, PRN *prn);

void gretl_functions_cleanup (void);

fnargs *fn_args_new (void);

void fn_args_free (fnargs *args);

int push_fn_arg (fnargs *args, int type, void *p);

void adjust_indent (const char *line, int *this_indent,
		    int *next_indent);


#endif /* GRETL_FUNC_H */
