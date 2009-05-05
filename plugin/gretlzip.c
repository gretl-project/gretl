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
#include <string.h>
#include <glib.h>

#include "zipunzip.h"

int gretl_unzip_file (const char *fname, GError **gerr)
{
    /* for verbose operation, make 3rd arg ZIP_VERBOSE or ZIP_TRACE */
    return zipfile_extract_files(fname, NULL, 0, gerr);
}

/*
 * @fname: full path to zipfile to be created.
 * @path: path relative to userdir for files to be picked up
 * and zipped.
 */

int gretl_make_zipfile (const char *fname, const char *path,
			GError **gerr)
{
    const char *array[2] = { path, NULL };
    int err;

    err = zipfile_archive_files(fname, array, 9, 
				ZIP_RECURSE_DIRS,
				gerr);
    return err;
}

int gretl_is_zipfile (const char *fname)
{
    zipinfo *zinfo;
    int ret = 0;

    zinfo = zipfile_get_info(fname, 0, NULL);
    if (zinfo != NULL) {
	zipinfo_destroy(zinfo);
	ret = 1;
    }

    return ret;
}

gchar *gretl_zipfile_get_topdir (const char *fname)
{
    zipinfo *zinfo;
    gchar *topdir = NULL;

    zinfo = zipfile_get_info(fname, 0, NULL);

    if (zinfo != NULL) {
	int i, n, gotit = 0;
	const gchar *s;

	for (i=0; i<zinfo->nfiles && !gotit; i++) {
	    s = zinfo->fnames[i];
	    if (s != NULL) {
		n = strlen(s);
		if (n > 13 && !strcmp(s + n - 11, "session.xml")) {
		    topdir = g_strndup(s, n - 11);
		    if (topdir != NULL) {
			n = strlen(topdir);
			if (topdir[n-1] == '/' || topdir[n-1] == '\\') {
			    topdir[n-1] = '\0';
			}
		    }
		}
	    }
	}
	zipinfo_destroy(zinfo);
    }

    return topdir;
}


