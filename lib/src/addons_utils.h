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

int is_gretl_addon (const char *pkgname);

const char **get_addon_names (int *n);

int get_addon_basic_info (const char *addon,
			  char **date,
			  char **descrip,
			  char **path);

int update_addons_index (PRN *prn);

int gretl_is_updated (const char *prev_build);

char *gretl_addon_get_path (const char *addon);

char *get_addon_examples_dir (const char *addon);

char *get_addon_pdf_path (const char *addon);
