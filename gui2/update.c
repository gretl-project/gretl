/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "gretl.h"
#include "gretl_www.h"

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

# ifdef WIN32

static size_t get_bufsize (const char *buf)
{
    size_t i, newsize = 0L;
    int pos;
    char line[60];

    while ((pos = haschar('\n', buf)) > 0) {
	strncpy(line, buf, pos);
	line[pos] = 0;
	sscanf(line, "%*s %u", &i);
	newsize += i;
	buf += pos + 1;
    }

    return newsize;
}

# endif /* WIN32 */

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

#ifdef WIN32

static void maybe_fork_updater (char *msg)
{
    int resp;

    resp = yes_no_dialog("gretl", msg, 0);

    if (resp == GRETL_YES) {
	gchar *ud;
	size_t n = strlen(paths.gretldir);

	if (paths.gretldir[n-1] != SLASH) {
	    ud = g_strdup_printf("%s\\gretl_updater.exe -g", paths.gretldir);
	} else {
	    ud = g_strdup_printf("%sgretl_updater.exe -g", paths.gretldir);
	}
	WinExec(ud, SW_SHOWNORMAL);
	exit(EXIT_SUCCESS);
    }
}

static void win_new_files_response (const char *buf)
{
    char infotxt[512];

    sprintf(infotxt, _("New files are available from the gretl web site.\n"
		       "These files have a combined size of %u bytes.\n\nWould "
		       "you like to exit from gretl and update your installation now?\n"
		       "(You can run gretl_updater.exe later if you prefer.)"),
	    get_bufsize(buf));
    maybe_fork_updater(infotxt);
}

#else

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

#endif

static int real_update_query (int queryopt)
{
    int err = 0;
    char *getbuf = NULL;
    char testfile[MAXLEN];
#ifndef WIN32
    int admin = 0;
    char hometest[MAXLEN];
#endif
    struct stat fbuf;
    time_t filedate = (time_t) 0;

    build_path(testfile, paths.gretldir, "gretl.stamp", NULL);

    if (stat(testfile, &fbuf)) {
	fprintf(stderr, "update_query: couldn't stat testfile '%s'\n", 
		testfile);
	return 1;
    } else {
	filedate = get_time_from_stamp_file(testfile);
#ifndef WIN32
	*hometest = '\0';
	if (getuid() != fbuf.st_uid) { 
	    /* user is not owner of gretl.stamp */
	    build_path(hometest, paths.userdir, "gretl.stamp", NULL);
	    if (!stat(hometest, &fbuf)) {
		filedate = get_time_from_stamp_file(hometest);
	    }
	} else {
	    admin = 1;
	}
#endif /* !WIN32 */
    }

    if (filedate == (time_t) 0) {
	fprintf(stderr, "update_query: couldn't get time from stamp file\n"); 
	return 1;
    }

    err = get_update_info(&getbuf, filedate, queryopt);

    if (err || getbuf == NULL) return 1;

    if (strncmp(getbuf, "message:", 8) == 0) {
	infobox(getbuf + 9);
    } else if (strncmp(getbuf, "No new", 6)) {
#ifdef WIN32
	win_new_files_response(getbuf);
#else
	new_files_response(admin, testfile, hometest);
#endif 
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
