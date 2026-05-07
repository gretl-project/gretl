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

/* load_functions.c for gretl_edit */

#include "gretl.h"
#include "libset.h"
#include "cmd_private.h"
#include "load_functions.h"

static CMD edcmd;
static int cmd_init_done;

static char *get_input_line (char *line,
			     const char *buf,
			     int *err)
{
    char *got;
    int n;

    *line = '\0';
    got = bufgets(line, MAXLINE, buf);

    if (got != NULL) {
        n = strlen(line);
        if (n > MAXLINE - 2  && line[n-1] != '\n') {
            *err = E_TOOLONG;
        }
    }

    return got;
}

static int exec_line (ExecState *s)
{
    int err = 0;

    if (gretl_compiling_function()) {
        return gretl_function_append_line(s);
    } else if (!strncmp(s->line, "function ", 9)) {
	err = parse_command_line(s, NULL, NULL);
    }

    if (!err && s->cmd->ci == FUNC) {
	err = gretl_cmd_exec(s, NULL);
    }

    return err;
}

int load_functions (const char *buf)
{
    ExecState state;
    char line[MAXLINE] = {0};
    char tmp[MAXLINE] = {0};
    int exec_err = 0;

    if (buf == NULL || *buf == '\0') {
	errbox(_("No functions to load"));
	return -1;
    }    

    if (!cmd_init_done) {
	gretl_cmd_init(&edcmd);
	cmd_init_done = 1;
    }

    gretl_set_batch_mode(1);
    bufgets_init(buf);
    gretl_exec_state_init(&state, 0, line, &edcmd, NULL, NULL);

    while (edcmd.ci != QUIT) {
	char *gotline = NULL;

	gotline = get_input_line(line, buf, &exec_err);
	if (gotline == NULL || exec_err) {
	    break;
	}
	if (state.in_comment) {
	    tailstrip(line);
	} else {
	    int contd = top_n_tail(line, sizeof line, &exec_err);
		
	    while (contd && !state.in_comment && !exec_err) {
		/* handle continued lines */
		gotline = get_input_line(tmp, buf, &exec_err);
		if (gotline == NULL) {
		    break;
		}
		if (!exec_err && *tmp != '\0') {
		    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
			exec_err = E_TOOLONG;
			break;
		    } else {
			strcat(line, tmp);
			compress_spaces(line);
		    }
		}
		contd = top_n_tail(line, sizeof line, &exec_err);
	    }
	}

	if (!exec_err) {
	    state.flags = SCRIPT_EXEC;
	    exec_err = exec_line(&state);
	}

	if (exec_err) {
	    gui_errmsg(exec_err);
	    break;
	}
    }

    bufgets_finalize(buf);

    return exec_err;
}
