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
#include "gui_recode.h"

#include <glib/gstdio.h>

/* Used when, e.g. loading a script into a GTK window: the
   script might be in a locale encoding other than UTF-8
   if it was created in a third-party editor.
*/

static gchar *real_my_locale_to_utf8 (const gchar *src,
				      int starting)
{
    static int errcount;
    const gchar *charset = NULL;
    gsize bytes;
    GError *err = NULL;
    GError **errp = NULL;
    gchar *ret;

    if (starting) {
	errcount = 0;
    }

    /* avoid multiple repetition of error message */
    errp = errcount == 0 ? &err : NULL;

    if (g_get_charset(&charset)) {
	/* we're in a UTF-8 locale */
	ret = g_convert(src, -1, "UTF-8", "ISO-8859-15",
			NULL, &bytes, errp);
    } else {
	ret = g_locale_to_utf8(src, -1, NULL, &bytes, errp);
    }

    if (err) {
	errbox(err->message);
	g_error_free(err);
	errcount++;
    }

    return ret;
}

/* wrappers to avoid repeated error messages where one
   will make the point
*/

gchar *my_locale_to_utf8 (const gchar *src)
{
    return real_my_locale_to_utf8(src, 1);
}

gchar *my_locale_to_utf8_next (const gchar *src)
{
    return real_my_locale_to_utf8(src, 0);
}
