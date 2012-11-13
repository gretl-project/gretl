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
#include "cmdstack.h"
#include "lib_private.h"
#include "session.h"

/* This module contains apparatus for handling a log of commands
   generated via point-and-click, or via the GUI console.  If
   the user is working in the context of a saved session file, the
   log goes to "session.inp" inside the session directory,
   otherwise it goes to a temporary file in the user's personal
   "dotdir".

   The command log may well be unusable as a script without
   modification, due to possible "nonlinearities" in the GUI (e.g.
   having several model windows open at once and performing
   tests on the models out of the order of estimation), but
   nonetheless it may be useful as a record and for pedagogical
   purposes.

   We flag errors via the GUI only in response to specific requests
   to view or update the log.
*/

#define CMD_DEBUG 0

static char logname[FILENAME_MAX]; /* filename for log */
static PRN *logprn;                /* log printer */
static int n_cmds;                 /* number of commands logged */
static int prev_ID;                /* keep track of model ID */
static int session_open;           /* are we doing the session file
				      thing? (0/1) */

/* Close down the logfile writing apparatus.  If we were writing to a
   temporary file, this will be deleted by gretl_print_destroy.
   Designed so that it doesn't hurt to call this function more
   times than is strictly necessary (to ensure cleanup).
*/

void free_command_stack (void)
{
#if CMD_DEBUG
    fprintf(stderr, "free_command_stack\n");
#endif

    if (logprn != NULL) {
	gretl_print_destroy(logprn);
	logprn = NULL;
    }

    n_cmds = prev_ID = 0;

    *logname = '\0';
}

/* For use with an existing session log file, whose name is registered
   in 'logname' */

static int session_logfile_init (void)
{
    char timestr[64];
    FILE *fp;

    fp = gretl_fopen(logname, "a");
    if (fp == NULL) {
	return E_FOPEN;
    }
    
    logprn = gretl_print_new_with_stream(fp); 
    if (logprn == NULL) {
	fclose(fp);
	return E_FOPEN;
    }

    prev_ID = 0;
    
#if CMD_DEBUG
    fprintf(stderr, "session_logfile_init: open prn for '%s'\n", logname);
#endif

    print_time(timestr);
    pprintf(logprn, _("# Log re-started %s\n"), timestr);

    return 0;
}

/* For use when the current gretl session is not associated
   with a session file (starting from scratch): we use a
   PRN with a temporary file */

static int scratch_logfile_init (void)
{
    const char *fname;
    char timestr[64];
    int err = 0;

    n_cmds = prev_ID = 0;
    *logname = '\0';

    logprn = gretl_print_new_with_tempfile(&err); 
    if (err) {
	return err;
    }

    fname = gretl_print_get_tempfile_name(logprn);

    if (fname == NULL) {
	gretl_print_destroy(logprn);
	logprn = NULL;
	return E_FOPEN; /* a bit vague */
    }

    strcpy(logname, fname);

#if CMD_DEBUG
    fprintf(stderr, "logfile_init: open prn for '%s'\n", logname);
#endif

    print_time(timestr);
    pprintf(logprn, _("# Log started %s\n"), timestr);

    pputs(logprn, _("# Record of session commands.  Please note that this will\n"
		    "# likely require editing if it is to be run as a script.\n"));

    return 0;
}

static int logfile_init (void)
{
    if (session_open && *logname != '\0') {
	return session_logfile_init();
    } else {
	return scratch_logfile_init();
    }
}

static int real_write_log_entry (const char *s)
{
    int err = 0;

#if CMD_DEBUG
    fprintf(stderr, "real_write_log_entry: logname='%s', logprn=%p\n",
	    logname, (void *) logprn);
#endif

    if (logprn == NULL) {
	err = logfile_init();
    }

    if (!err) {
	int n = strlen(s);

	pputs(logprn, s);
	if (s[n-1] != '\n') {
	    pputc(logprn, '\n');
	}
	gretl_print_flush_stream(logprn);
    }

    return err;
}

/* for a given GUI command (not associated with a model): place it in
   the log buffer */

int add_command_to_stack (const char *s)
{
    int err;

    if (s == NULL || *s == '\0') {
	return 1;
    }

    /* not a model command, so zero out record of
       previous model ID number */
    prev_ID = 0;

    err = real_write_log_entry(s);

    if (!err) {
#if CMD_DEBUG
	fprintf(stderr, "Written to log: '%s'\n", s);
#endif 
	n_cmds++;

	if (strlen(s) > 2 && *s != '#' && 
	    strncmp(s, "help", 4) &&
	    strncmp(s, "info", 4) &&
	    strncmp(s, "list", 4) &&
	    strncmp(s, "open", 4) &&
	    strncmp(s, "quit", 4)) {
	    set_commands_recorded();
	}
    }

    return err;
}

int add_model_command_to_stack (const char *s, int model_ID)
{
    int err = 0;

    if (s == NULL || *s == '\0') {
	return 1;
    }

    if (logprn == NULL) {
	err = logfile_init();
    } 

    if (!err) {
	if (model_ID != prev_ID) {
	    pprintf(logprn, "# %s %d\n", _("model"), model_ID);
	    prev_ID = model_ID;
	}
	err = real_write_log_entry(s);
	if (!err) {
	    set_commands_recorded();
	}
    }

    return err;
}

/* Called in response to the refresh/reload button in the viewer
   window for the command log: retrieve the updated log content.
   Display of any error messages is handled by the caller, which
   is also responsible for freeing the value returned by this
   function.
*/

gchar *get_logfile_content (int *err)
{
    gchar *s = NULL;

    if (n_cmds > 0) {
	*err = gretl_file_get_contents(logname, &s, NULL);
    }

    return s;
}

/* called from the main window /Tools menu */

void view_command_log (void)
{
    if (!session_open && n_cmds == 0) {
	warnbox(_("The command log is empty"));
    } else {
	view_file(logname, 0, 0, 78, 370, VIEW_LOG);
    }
}

/* 
   This function gets called from session.c when when saving a
   session, opening a session file, or closing a session.

   On first saving a session, we shift the existing logfile (if any)
   from the user's dotdir into the session directory so it'll get
   saved along with the other materials.  If the save is of a session
   that has already been saved to file, however, the logfile location
   should already be correct.

   On opening a session file, we set the logfile name so the re-opened
   session log can be displayed and added to.

   On closing a session, we close the current session logfile and
   redirect the log to a tempfile in the user's "dotdir".
*/

void set_session_log (const char *dirname, int code)
{
    char tmp[FILENAME_MAX];
    int err;

#if CMD_DEBUG
    fprintf(stderr, "set_session_log: dirname = '%s'\n", dirname);
    fprintf(stderr, "session_open = %d\n", session_open);
#endif

    if (code == LOG_SAVE_AS && logprn == NULL) {
	/* previous log file has been closed; we'll
	   let it reopen as and when needed
	*/
	strcpy(logname, dirname);
	strcat(logname, "session.inp");
	session_open = 1;
    } else if (code == LOG_SAVE || code == LOG_SAVE_AS) {
	if (gretl_print_has_tempfile(logprn)) {
	    /* if logprn is NOT a tempfile-prn, it will already
	       have been renamed, via the renaming of the
	       session directory on save
	    */
	    strcpy(tmp, dirname);
	    strcat(tmp, "session.inp");
	    if (strcmp(logname, tmp)) {
		if (logprn != NULL) {
		    err = gretl_print_rename_file(logprn, logname, tmp);
		    if (err) {
			free_command_stack();
		    }
		}
		strcpy(logname, tmp);
		session_open = 1;
	    }
	}
    } else if (code == LOG_OPEN) {
	free_command_stack();
	strcpy(logname, dirname);
	strcat(logname, "session.inp");
	session_open = 1;
    } else if (code == LOG_CLOSE) {
	free_command_stack();
	session_open = 0;
    } else if (code == LOG_NULL) {
	session_open = 0;
    }
}

void maybe_suspend_session_log (void)
{
    /* If we're saving an open session under another name, we should
       close the old session.inp file, if present (this will be
       reopened later). On Windows we get "permission denied" on
       trying to rename the session directory if it contains an open
       stream.

       The test for !gretl_print_has_tempfile catches the case
       where the current logprn is a tempfile: that means it's
       in the user's dotdir, not inside an old session directory.
       In that case it'll be moved by set_session_log().
    */
    if (logprn != NULL && !gretl_print_has_tempfile(logprn)) {
	gretl_print_close_stream(logprn);
	gretl_print_destroy(logprn);
	logprn = NULL;
    }
}

