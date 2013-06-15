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

#ifndef FNCALL_H
#define FNCALL_H

enum {
    FN_NO_LOAD,
    FN_NO_DATA
} fncall_errors;

void call_function_package (const char *fname, windata_t *vwin,
			    int *loaderr);

void function_call_cleanup (void);

void gui_define_list (void);

void fncall_register_genr (int addv, gpointer p);

gchar *get_bundle_plot_function (gretl_bundle *b);

int exec_bundle_plot_function (gretl_bundle *b, const char *funname);

int try_exec_bundle_print_function (gretl_bundle *b, PRN *prn);

void maybe_add_packages_to_model_menus (windata_t *vwin);

void maybe_add_packages_to_menus (windata_t *vwin);

int gui_add_package_to_menu (const char *path, gboolean prechecked,
			     GtkWidget *parent);

int maybe_handle_pkg_menu_option (const char *path, 
				  GtkWidget *parent);

int package_is_available_for_menu (const gchar *pkgname,
				   const char *path);

void unregister_function_package (const gchar *pkgname);

int query_addons (void);

int download_addon (const char *pkgname, char **local_path);

char *installed_addon_status_string (const char *path,
				     const char *svstr);

int revise_package_status (const gchar *pkgname,
			   const gchar *label,
			   const gchar *mpath,
			   gboolean uses_subdir,
			   gboolean maybe_edit,
			   gboolean installing);

#endif /* FNCALL_H */
