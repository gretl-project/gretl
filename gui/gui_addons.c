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

#include "gretl.h"
#include "gui_addons.h"
#include "gretl_www.h"
#include "addons_utils.h"
#include "gretl_untar.h"

/* @basename and @filepath can be given if the caller wants the path
   to a specific newly installed file on success, otherwise NULLs
   are OK.

   Note: the return value is non-zero iff download and installation
   are successful.
*/

DLCode maybe_download_addons (GtkWidget *parent,
			      const char *basename,
			      char **filepath)
{
    const char *msg = N_("You have selected an action that requires access to\n"
			 "the gretl addons. But these packages are missing,\n"
			 "incomplete, or not up to date.\n\n"
			 "Do you want to download and install the current\n"
			 "addons now?");
    int resp = yes_no_dialog(NULL, _(msg), parent);
    DLCode ret = DL_CANCEL;
    int err = 0;

    if (resp == GRETL_YES) {
	gchar *dlpath = get_download_path("addons.tar.gz", &err);

	if (!err) {
	    err = retrieve_addons_package(dlpath);
	}
	if (!err) {
	    err = unpack_files_collection(dlpath);
	}
	if (!err) {
	    update_addons_index(NULL);
	}

	if (err) {
	    ret = DL_FAIL;
	} else {
	    ret = DL_SUCCESS;
	    if (basename != NULL && filepath != NULL) {
		if (has_suffix(basename, ".pdf")) {
		    *filepath = get_addon_pdf_path(basename);
		} else {
		    *filepath = gretl_addon_get_path(basename);
		}
	    }
	}

	 gretl_remove(dlpath);
	 g_free(dlpath);
    }

    return ret;
}
