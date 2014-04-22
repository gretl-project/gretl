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

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

/* E.g. Sun Mar 16 13:50:52 EST 2003 */

static time_t get_time_from_stamp_file (const char *fname)
{
    FILE *fp;
    struct tm stime;
    char wday[4], mon[4];
    int i;
    const char *months[] = {
        "Jan", "Feb", "Mar",
        "Apr", "May", "Jun",
        "Jul", "Aug", "Sep",
        "Oct", "Nov", "Dec"
    };

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	fprintf(stderr, "get_time_from_stamp_file:\n"
		" couldn't open '%s'\n", fname);
	return (time_t) 0;
    }

    if (fscanf(fp, "%3s %3s %d %d:%d:%d %*s %d", 
               wday, mon, &stime.tm_mday, &stime.tm_hour,
               &stime.tm_min, &stime.tm_sec, &stime.tm_year) != 7) {
	fclose(fp);
        return (time_t) 0;
    }

    fclose(fp);

    stime.tm_mon = 20;
    for (i=0; i<12; i++) {
        if (!strcmp(mon, months[i])) {
            stime.tm_mon = i;
            break;
        }
    }

    if (stime.tm_mon == 20) {
	return (time_t) 0;
    }

    stime.tm_year -= 1900;

    return mktime(&stime);
}

static void new_files_response (int admin, const char *testfile,
				const char *hometest)
{
    char infotxt[512];
    FILE *fp;

    if (admin) {
	strcpy(infotxt, _("New files are available from the gretl web site\n"
			  "http://gretl.sourceforge.net/"));
	fp = fopen(testfile, "w");
    } else {
	strcpy(infotxt, _("You may want to let the system administrator know\n"
			  "that new files are available from the gretl web site\n"
			  "http://gretl.sourceforge.net/"));
	fp = fopen(hometest, "w");
    }

    if (fp != NULL) {
	fprintf(fp, _("This file is part of the gretl update notification "
		      "system\n"));
	fclose(fp);
    }

    infobox(infotxt);
}

static int real_update_query (int queryopt)
{
    char *getbuf = NULL;
    char testfile[MAXLEN];
    char hometest[MAXLEN];
    struct stat fbuf;
    time_t filedate = (time_t) 0;
    int admin = 0;
    int err = 0;

    build_path(testfile, gretl_home(), "gretl.stamp", NULL);

    if (gretl_stat(testfile, &fbuf)) {
	fprintf(stderr, "update_query: couldn't stat testfile '%s'\n", 
		testfile);
	return 1;
    } else {
	filedate = get_time_from_stamp_file(testfile);
	*hometest = '\0';
	if (getuid() != fbuf.st_uid) { 
	    /* the user is not the owner of gretl.stamp */
	    build_path(hometest, gretl_dotdir(), "gretl.stamp", NULL);
	    if (gretl_stat(hometest, &fbuf) == 0) {
		filedate = get_time_from_stamp_file(hometest);
	    }
	} else {
	    admin = 1;
	}
    }

    if (filedate == (time_t) 0) {
	fprintf(stderr, "update_query: couldn't get time from stamp file\n"); 
	return 1;
    }

    err = get_update_info(&getbuf, filedate, queryopt);

    if (err || getbuf == NULL) {
	return 1;
    }

    if (strncmp(getbuf, "message:", 8) == 0) {
	infobox(getbuf + 9);
    } else if (strncmp(getbuf, "No new", 6)) {
	new_files_response(admin, testfile, hometest);
    } else if (queryopt == QUERY_VERBOSE) {
	infobox(_("No new files"));
    }

    free(getbuf);

    return err;
}

int silent_update_query (void)
{
    return real_update_query(QUERY_SILENT);
}

int update_query (void)
{
    return real_update_query(QUERY_VERBOSE);
}
