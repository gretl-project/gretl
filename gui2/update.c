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
#include "gretl_www.h"

static int real_update_query (int verbose)
{
    char *getbuf = NULL;
    int err;

    err = get_update_info(&getbuf, verbose);

    if (err || getbuf == NULL) {
	return 1;
    }

    if (strncmp(getbuf, "message:", 8) == 0) {
	infobox(getbuf + 9);
    } else if (verbose) {
	infobox(_("No new files"));
    }

    free(getbuf);

    return err;
}

#if 0
int silent_update_query (void)
{
    /* not used right now, but we may want to reactivate it
       in some form */
    return real_update_query(0);
}
#endif

int update_query (void)
{
    return real_update_query(1);
}
