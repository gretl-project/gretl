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

#include <stdio.h>
#include <stdlib.h>
#include <glib.h>

#include "zipunzip.h"

/* We'll try to take charge at this level of ensuring that filenames
   are in the right encoding for the OS, so we can use the
   regular stdio functions in the rest of the zipunzip apparatus.
*/

static char *get_sys_fname (const char *fname, GError **gerr)
{
    const gchar *cset;
    gchar *sfname = NULL;
    gsize bytes;

    if (g_get_charset(&cset)) {
	sfname = g_strdup(fname);
    } else {
#if defined(G_OS_WIN32) && GLIB_MAJOR_VERSION == 2 && GLIB_MINOR_VERSION >= 10
	sfname = g_locale_from_utf8(fname, -1, NULL, &bytes, gerr);
#else
	sfname = g_filename_from_utf8(fname, -1, NULL, &bytes, gerr);
#endif
    }

    return sfname;
}

int gretl_unzip_file (const char *fname, GError **gerr)
{
    char *sysfname;
    int err = 0;

    sysfname = get_sys_fname(fname, gerr);

    if (sysfname == NULL) {
	err = 1;
    } else {
	err = zipfile_extract_files(sysfname, NULL, 0, gerr);
	g_free(sysfname);
    }
    
    return err;
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    char *sysfname = NULL;
    char *sysdname = NULL;
    const char *array[2];
    int err = 0;

    sysfname = get_sys_fname(fname, gerr);
    if (sysfname == NULL) {
        err = 1;
    } else {
        sysdname = get_sys_fname(path, gerr);
        if (sysdname == NULL) {
            err = 1;
        }
    }

    if (!err) {
	array[0] = sysdname;
	array[1] = NULL;

	err = zipfile_archive_files(sysfname, array, 9, 
				    ZIP_RECURSE_DIRS,
				    gerr);
    }

    g_free(sysfname);
    g_free(sysdname);

    return err;
}

int gretl_is_zipfile (const char *fname)
{
    char *sysfname;
    zipinfo *zinfo;
    int ret = 0;

    sysfname = get_sys_fname(fname, NULL);
    
    if (sysfname != NULL) {
	zinfo = zipfile_get_info(sysfname, 0, NULL);
	if (zinfo != NULL) {
	    zipinfo_destroy(zinfo);
	    ret = 1;
	}
	g_free(sysfname);
    }

    return ret;
}

