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

#define CMD_DEBUG 1

static char **cmd_stack;
static int n_cmds;

void free_command_stack (void)
{
    free_strings_array(cmd_stack, n_cmds);
    cmd_stack = NULL;
    n_cmds = 0;
}

int add_command_to_stack (const char *s)
{
    int err = 0;

    if (s == NULL || *s == '\0') {
	return 1;
    }

    err = strings_array_add(&cmd_stack, &n_cmds, s);

#if CMD_DEBUG
    fprintf(stderr, "Added to stack: '%s'\n", s);
    fprintf(stderr, "err = %d, n_cmds = %d\n", err, n_cmds);    
#endif    

    if (strlen(s) > 2 && 
	strncmp(s, "help", 4) &&
	strncmp(s, "info", 4) &&
	strncmp(s, "list", 4) &&
	strncmp(s, "open", 4) &&
	strncmp(s, "quit", 4)) {
	mark_session_changed();
    }

    return err;
}

void delete_last_command (void)
{
    if (n_cmds > 0) {
	free(cmd_stack[n_cmds - 1]);
	n_cmds--;
    }
}

/* print to the command log a command associated with a 
   particular model */

int model_command_init (int model_ID)
{
    char *line;
    CMD *libcmd;
    PRN *prn;
    int err = 0;

    line = get_lib_cmdline();
    libcmd = get_lib_cmd();

    /* pre-process the line */
    if (check_specific_command(line)) {
	return 1;
    }

    if (bufopen(&prn)) {
	return 1;
    }

    pprintf(prn, "# model %d\n", model_ID);
    echo_cmd(libcmd, datainfo, line, CMD_RECORDING, prn);

    err = add_command_to_stack(gretl_print_get_buffer(prn));
    
    gretl_print_destroy(prn);

    return err;
}

static char logname[FILENAME_MAX];
static time_t logtime;

/* ship out the stack of commands entered in the current session */

static int update_logfile (void)
{
    const char *s;
    FILE *fp;
    int i, n;

    fp = gretl_fopen(logname, "w"); 
    
    if (fp == NULL) {
	file_write_errbox(logname);
	return 1;
    }

    fprintf(fp, "Record started %s\n", print_time(&logtime));

    fputs("# Record of commands issued in the current session.  Please note\n"
	  "# that this may require editing if it is to be run as a script.\n", 
	  fp);

    for (i=0; i<n_cmds; i++) {
	s = cmd_stack[i];
	n = strlen(s);
	fputs(s, fp);
	if (s[n-1] != '\n') {
	    fputc('\n', fp);
	}
    }

    fclose(fp);

    return 0;
}

gchar *get_logfile_content (int *err)
{
    gchar *s = NULL;

    *err = update_logfile();

    if (!*err) {
	*err = gretl_file_get_contents(logname, &s);
    }

    return s;
}

void view_command_log (void)
{
    if (n_cmds == 0) {
	warnbox(_("The command log is empty"));
	return;
    }

    if (*logname == '\0') {
	strcpy(logname, paths.dotdir);
	strcat(logname, "session.inp");
	logtime = time(NULL);
    }

    if (update_logfile()) {
	return;
    }

    view_file(logname, 0, 0, 78, 370, VIEW_LOG);
}
