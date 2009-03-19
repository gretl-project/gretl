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
static char *logline;              /* buffered command line */
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
    if (logline != NULL) {
	free(logline);
	logline = NULL;
    }

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
    time_t logtime;
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

    logtime = time(NULL);
    pprintf(logprn, _("# Log re-started %s\n"), print_time(&logtime));

    return 0;
}

/* For use when the current gretl session is not associated
   with a session file (starting from scratch): we use a
   PRN with a temporary file */

static int scratch_logfile_init (void)
{
    const char *fname;
    time_t logtime;
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

    logtime = time(NULL);
    pprintf(logprn, _("# Log started %s\n"), print_time(&logtime));

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

/* write out the last stacked command line, if any, and
   flush the log stream */

static int flush_logfile (void)
{
    int n_cmds_save;
    int err;

    if (n_cmds == 0) {
	/* nothing to be done */
	return 0;
    }

    gretl_error_clear();
    n_cmds_save = n_cmds;

#if CMD_DEBUG
    fprintf(stderr, "flush_logfile: logname='%s', logprn=%p\n",
	    logname, (void *) logprn);
#endif

    if (logprn == NULL) {
	err = logfile_init();
	if (err) {
	    return err;
	}
    }

    n_cmds = n_cmds_save;

#if CMD_DEBUG
    fprintf(stderr, "n_cmds = %d, logline='%s'\n", n_cmds, logline);
#endif

    n_cmds = n_cmds_save;

    if (logline != NULL) {
	int n = strlen(logline);

	pputs(logprn, logline);
	if (logline[n-1] != '\n') {
	    pputc(logprn, '\n');
	}
	free(logline);
	logline = NULL;
    }

    gretl_print_flush_stream(logprn);

    return 0;
}

/* in case an error was discovered after buffering a given command,
   strike it from the record */

void delete_last_command (void)
{
    if (logline != NULL) {
	free(logline);
	logline = NULL;
	n_cmds--;
    }
}

/* for a given GUI command (not associated with a model): place it in
   the log buffer */

int add_command_to_stack (const char *s)
{
    int err;

    if (s == NULL || *s == '\0') {
	return 1;
    }

    /* not a model command, so delete record of
       previous model ID number */
    prev_ID = 0;

    /* is there a previous buffered command? if so, send it to
       the logfile first */
    if (logline != NULL) {
	err = flush_logfile();
	if (err) {
	    return err;
	}
    }

    /* buffer the current line in case we need to delete it */
    logline = gretl_strdup(s);
    if (logline == NULL) {
	return E_ALLOC;
    }

#if CMD_DEBUG
    fprintf(stderr, "Added to stack: '%s'\n", s);
#endif 

    n_cmds++;

    if (strlen(s) > 2 && 
	strncmp(s, "help", 4) &&
	strncmp(s, "info", 4) &&
	strncmp(s, "list", 4) &&
	strncmp(s, "open", 4) &&
	strncmp(s, "quit", 4)) {
	mark_session_changed();
    }

    return 0;
}

/* save to the command log a command associated with a particular
   model */

int model_command_init (int model_ID)
{
    char *line = get_lib_cmdline();
    CMD *libcmd = get_lib_cmd();
    int err = 0;

#if CMD_DEBUG
    fprintf(stderr, "model_command_init: line='%s'\n", line);
#endif

    /* pre-process the line */
    if (check_specific_command(line)) {
	return 1;
    }

    n_cmds++;

    err = flush_logfile();

    if (err) {
	n_cmds--;
    } else {
	if (model_ID != prev_ID) {
	    pprintf(logprn, "# %s %d\n", _("model"), model_ID);
	    prev_ID = model_ID;
	}
	echo_cmd(libcmd, datainfo, line, CMD_RECORDING, logprn);
	mark_session_changed();
    } 

    return err;
}

/* Called in response to the refresh/reload button in the viewer
   window for the command log: retrieve the updated log content.
   Display of any error messages is handled by the caller.
*/

gchar *get_logfile_content (int *err)
{
    gchar *s = NULL;

    if (n_cmds == 0) {
	return NULL;
    }

    *err = flush_logfile();

    if (!*err) {
	*err = gretl_file_get_contents(logname, &s);
    }

    return s;
}

/* called from main menu: /Tools/View command log */

void view_command_log (void)
{
    if (!session_open && n_cmds == 0) {
	warnbox(_("The command log is empty"));
    } else {
	int err = flush_logfile();

	if (err) {
	    gui_errmsg(err);
	} else {
	    view_file(logname, 0, 0, 78, 370, VIEW_LOG);
	}
    }
}

/* 
   This function gets called from session.c when when saving a
   session, opening a session file, or closing a session.

   On saving, we shift the existing logfile (if any) into the session
   directory so it'll get saved along with the other materials -- if
   required.  (But if the save is of a session that has already been
   saved to file, the logfile location will already be correct.)

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

    flush_logfile();

    if (code == LOG_SAVE) {
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
