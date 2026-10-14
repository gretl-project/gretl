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
   functions, we check if any are currently loaded. If not, we
   automatically load any function definitions in the current editor
   tabs.
*/

void maybe_load_functions (windata_t *vwin)
{
    if (n_user_functions() == 0) {
	load_functions(vwin, TRUE);
    }
}

/* Determine whether @tview contains any function definitions, using a
   basic heuristic.
*/

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

/* Apparatus for assembling a list of functions whose definitions appear
   in the buffer, along with the point in the buffer at which they
   start.
*/

typedef struct {
    gchar *name;  /* name of function */
    gint offset;  /* character offset of the name in the buffer */
} fndef;

static GRegex *fndef_re;

static void init_fndef_regex (void)
{
    GError *err = NULL;

    fndef_re = g_regex_new("^function[ \\t]+"
                           "(?:void|scalar|series|matrix|list|string|bundle"
                           "|strings|matrices|bundles|lists|arrays|numeric)"
                           "[ \\t]+(?<name>[A-Za-z][A-Za-z0-9_]{0,30})(?![A-Za-z0-9_])",
                           G_REGEX_MULTILINE | G_REGEX_OPTIMIZE, 0, &err);
    if (err != NULL) {
        g_warning("%s", err->message);
        g_error_free(err);
    }
}

/* Collect all function definitions in @buf, in order of appearance */

static GArray *get_fn_definitions (GtkTextBuffer *buf)
{
    GArray *defs = g_array_new(FALSE, FALSE, sizeof(fndef));
    GtkTextIter start, end;
    GMatchInfo *mi;
    gchar *text, *prev;
    gint offset = 0;

    if (fndef_re == NULL) {
	init_fndef_regex();
    }

    gtk_text_buffer_get_bounds(buf, &start, &end);
    text = gtk_text_buffer_get_text(buf, &start, &end, TRUE);
    prev = text;

    g_regex_match(fndef_re, text, 0, &mi);

    while (g_match_info_matches(mi)) {
        fndef d;
        gint bstart;

        g_match_info_fetch_named_pos(mi, "name", &bstart, NULL);
        /* convert byte offset to char offset incrementally */
        offset += g_utf8_pointer_to_offset(prev, text + bstart);
        prev = text + bstart;
        d.name = g_match_info_fetch_named(mi, "name");
        d.offset = offset;
        g_array_append_val(defs, d);
        g_match_info_next(mi, NULL);
    }

    g_match_info_free(mi);
    g_free(text);

    return defs;
}

static void free_fn_definitions (GArray *defs)
{
    guint i;

    for (i=0; i<defs->len; i++) {
        g_free(g_array_index(defs, fndef, i).name);
    }
    g_array_free(defs, TRUE);
}

static void jump_to_def (GtkMenuItem *item, GtkTextView *view)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    gint offset = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item), "offset"));
    GtkTextIter it;

    gtk_text_buffer_get_iter_at_offset(buf, &it, offset);
    gtk_text_buffer_place_cursor(buf, &it);
    gtk_text_view_scroll_to_mark(view, gtk_text_buffer_get_insert(buf),
                                 0.0, TRUE, 0.0, 0.3);
    gtk_widget_grab_focus(GTK_WIDGET(view));
}

void show_fn_menu (GtkWidget *button, GtkTextView *view)
{
    GArray *defs = get_fn_definitions(gtk_text_view_get_buffer(view));
    GtkWidget *menu = gtk_menu_new();
    GtkWidget *item;
    guint i;

    if (defs->len == 0) {
        item = gtk_menu_item_new_with_label(_("No functions defined"));
        gtk_widget_set_sensitive(item, FALSE);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    for (i=0; i<defs->len; i++) {
        fndef *d = &g_array_index(defs, fndef, i);

        item = gtk_menu_item_new_with_label(d->name);
        g_object_set_data(G_OBJECT(item), "offset", GINT_TO_POINTER(d->offset));
        g_signal_connect(item, "activate", G_CALLBACK(jump_to_def), view);
        gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }
    free_fn_definitions(defs);

    /* storing the menu on the button destroys the previous one */
    g_object_set_data_full(G_OBJECT(button), "fn-menu", menu,
                           (GDestroyNotify) gtk_widget_destroy);
    gtk_widget_show_all(menu);
    gtk_menu_popup_at_widget(GTK_MENU(menu), button, GDK_GRAVITY_SOUTH_WEST,
                             GDK_GRAVITY_NORTH_WEST, NULL);
}
