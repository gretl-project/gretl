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

#define CMD_DEBUG 0

static char logname[FILENAME_MAX];
static time_t logtime;
static char *logline;
static PRN *logprn;
static int n_cmds;
static int prev_ID;
static int session_open;

/* For use with existing session log file */

static int logfile_reinit (void)
{
    FILE *fp;

    fp = gretl_fopen(logname, "a");
    if (fp == NULL) {
	file_write_errbox(logname);
	return 1;
    }
    
    logprn = gretl_print_new_with_stream(fp); 
    if (logprn == NULL) {
	file_write_errbox(logname);
	return 1;
    }

    prev_ID = 0;
    logtime = time(NULL);
    
#if CMD_DEBUG
    fprintf(stderr, "logfile_reinit: open prn for '%s'\n", logname);
#endif

    pprintf(logprn, "Log re-started %s\n", print_time(&logtime));

    return 0;
}

/* make the name of the logfile, open a stream for writing, 
   and write header text */

static int logfile_init (void)
{
    int err = 0;

    if (session_open && *logname != '\0') {
	return logfile_reinit();
    }
    
    n_cmds = prev_ID = 0;

    if (*logname == '\0') {
	strcpy(logname, paths.dotdir);
	strcat(logname, "session.inp");
    }

    logtime = time(NULL);

    logprn = gretl_print_new_with_filename(logname, &err); 
    
    if (err) {
	file_write_errbox(logname);
	return err;
    }

#if CMD_DEBUG
    fprintf(stderr, "logfile_init: open prn for '%s'\n", logname);
#endif

    pprintf(logprn, "Log started %s\n", print_time(&logtime));

    pputs(logprn, "# Record of commands issued in the current session.  Please note\n"
	  "# that this may require editing if it is to be run as a script.\n");

    return 0;
}

/* close down the logfile writing apparatus */

void free_command_stack (void)
{
    if (logline != NULL) {
	free(logline);
	logline = NULL;
    }

    if (logprn != NULL) {
	/* close but do not remove file */
	gretl_print_destroy(logprn);
	logprn = NULL;
    }

    n_cmds = prev_ID = 0;
}

/* write out the last stacked command line, if any, and
   flush the log stream */

static int flush_logfile (void)
{
    int err;

    if (n_cmds == 0) {
	/* nothing to be done */
	return 0;
    }

    if (logprn == NULL) {
	err = logfile_init();
	if (err) {
	    return err;
	}
    }

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

int add_command_to_stack (const char *s)
{
    int err;

    if (s == NULL || *s == '\0') {
	return 1;
    }

#if 0
    fprintf(stderr, "session dir = '%s'\n", get_session_dirname());
#endif

    /* not a model command, so delete record of
       previous model ID number */
    prev_ID = 0;

    /* is there a previous stacked command? if so, send it to
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

void delete_last_command (void)
{
    if (logline != NULL) {
	free(logline);
	logline = NULL;
    }
}

/* save to the command log a command associated with a 
   particular model */

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

    err = flush_logfile();

    if (!err) {
	if (model_ID != prev_ID) {
	    pprintf(logprn, "# model %d\n", model_ID);
	    prev_ID = model_ID;
	}
	echo_cmd(libcmd, datainfo, line, CMD_RECORDING, logprn);
	mark_session_changed();
	n_cmds++;
    }

    return err;
}

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

void view_command_log (void)
{
    if (!session_open && n_cmds == 0) {
	warnbox(_("The command log is empty"));
    } else {
	int err = flush_logfile();

	if (!err) {
	    view_file(logname, 0, 0, 78, 370, VIEW_LOG);
	}
    }
}

/* 
   This function gets called when when saving a session, opening a
   session file, or closing a session.

   On saving, we shift the logfile into the session directory so it'll
   get saved along with the other materials.  On opening a session, we
   set the logfile name so the re-opened log can be displayed.  On
   closing a session, we redirect the log to the default location in
   the user's "dotdir"; this is signalled by a NULL value for the
   dirname argument.
*/

void set_session_log (const char *dirname, int code)
{
    char tmp[FILENAME_MAX];
    int err;

    flush_logfile();

    if (code == SAVE_SESSION) {
	strcpy(tmp, dirname);
	strcat(tmp, "session.inp");
	if (strcmp(logname, tmp)) {
	    err = gretl_print_rename_file(logprn, logname, tmp);
	    if (err) {
		free_command_stack();
	    }
	    strcpy(logname, tmp);
	    session_open = 1;
	}
    } else {
	/* closing or opening session */
	free_command_stack();
	if (dirname != NULL) {
	    /* opening */
	    strcpy(logname, dirname);
	    strcat(logname, "session.inp");
	    session_open = 1;
	} else {
	    /* closing session */
	    *logname = '\0';
	    session_open = 0;
	}
    }
}
