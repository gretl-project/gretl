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
#include "textbuf.h"
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
        err = gretl_function_append_line(s);
    } else if (!strncmp(s->line, "function ", 9)) {
	err = parse_command_line(s, NULL, NULL);
    }

    if (!err && s->cmd->ci == FUNC) {
	err = gretl_cmd_exec(s, NULL);
    }
    if (err) {
	gretl_errmsg_prepend(s->line, err);
    }

    return err;
}

static int real_load_functions (const gchar *buf,
				ExecState *state,
				char *line,
				char *tmp)
{
    int err = 0;

    bufgets_init(buf);

    while (state->cmd->ci != QUIT) {
	char *gotline = NULL;

	gotline = get_input_line(line, buf, &err);
	if (gotline == NULL || err) {
	    break;
	}
	if (state->in_comment) {
	    tailstrip(line);
	} else {
	    int contd = top_n_tail(line, sizeof line, &err);

	    while (contd && !state->in_comment && !err) {
		/* handle continued lines */
		gotline = get_input_line(tmp, buf, &err);
		if (gotline == NULL) {
		    break;
		}
		if (!err && *tmp != '\0') {
		    if (strlen(line) + strlen(tmp) > MAXLINE - 1) {
			err = E_TOOLONG;
			break;
		    } else {
			strcat(line, tmp);
			compress_spaces(line);
		    }
		}
		contd = top_n_tail(line, sizeof line, &err);
	    }
	}

	if (!err) {
	    state->flags = SCRIPT_EXEC;
	    err = exec_line(state);
	}

	if (err) {
	    gui_errmsg(err);
	    break;
	}
    }

    bufgets_finalize(buf);

    return err;
}

int load_functions (windata_t *vwin, gboolean all_tabs)
{
    ExecState state;
    char line[MAXLINE] = {0};
    char tmp[MAXLINE] = {0};
    gchar *buf = NULL;
    int err = 0;

    if (!cmd_init_done) {
	gretl_cmd_init(&edcmd);
	cmd_init_done = 1;
    }

    gretl_set_batch_mode(1);
    gretl_exec_state_init(&state, 0, line, &edcmd, NULL, NULL);

    if (all_tabs) {
	/* load functions from all hansl tabs */
	GtkNotebook *book = GTK_NOTEBOOK(editor_get_tabs(vwin));
	int np = gtk_notebook_get_n_pages(book);
	GtkWidget *tab;
	windata_t *viewer;
	GtkTextView *tview;
	int i;

	for (i=0; i<np; i++) {
	    tab = gtk_notebook_get_nth_page(book, i);
	    viewer = g_object_get_data(G_OBJECT(tab), "vwin");
	    tview = GTK_TEXT_VIEW(viewer->text);
	    if (viewer->role == EDIT_HANSL && textview_has_functions(tview)) {
		buf = textview_get_hansl(tview, 1);
		err = real_load_functions(buf, &state, line, tmp);
		g_free(buf);
	    }
	}
    } else {
	/* load functions from "this" tab only */
	buf = textview_get_hansl(GTK_TEXT_VIEW(vwin->text), 1);
	err = real_load_functions(buf, &state, line, tmp);
	g_free(buf);
    }

    return err;
}

/* If the user starts an action associated with look-up of user
   functions, we check if any are currently loaded. If not, and if
   @tview contains an instance of a line-start that indicates definition
   of a function, we'll automatically load any function definitions in
   the current buffer.
*/

void maybe_load_functions (windata_t *vwin)
{
    if (n_user_functions() == 0) {
	load_functions(vwin, TRUE);
    }
}

gboolean textview_has_functions (GtkTextView *tview)
{
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tview);
    GtkTextIter start, match_start, match_end;
    gboolean ret = FALSE;

    gtk_text_buffer_get_start_iter(tbuf, &start);

    while (gtk_text_iter_forward_search(&start, "function ",
					GTK_TEXT_SEARCH_TEXT_ONLY,
					&match_start, &match_end,
					NULL)) {
	if (gtk_text_iter_starts_line(&match_start)) {
	    ret = TRUE;
	    break;
	}
	start = match_end; /* continue the search */
    }

    return ret;
}
