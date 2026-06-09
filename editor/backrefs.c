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
#include "textbuf.h"
#include "gretl_func.h"
#include "gretl_typemap.h"
#include "backrefs.h"

/* Apparatus to support multi-level "Go back" references in a
   GtkTextView object.
*/

typedef struct backrefs_ {
    GtkTextMark **marks;
    int n_marks;
    int n_slots;
} backrefs;

/* Allocate a new backrefs struct. */

static backrefs *backrefs_new (void)
{
    backrefs *refs = malloc(sizeof *refs);

    refs->marks = malloc(sizeof *refs->marks);
    refs->n_marks = 0;
    refs->n_slots = 1;

    return refs;
}

/* Free a backrefs struct when its parent widget is destroyed. */

static void backrefs_destroy (gpointer data)
{
    if (data != NULL) {
	backrefs *refs = (backrefs *) data;

	free(refs->marks);
	free(refs);
    }
}

/* Push @mark onto a stack of backward references for @buf. */

static void push_backref (GtkTextBuffer *buf, GtkTextMark *mark)
{
    backrefs *refs = g_object_get_data(G_OBJECT(buf), "backrefs");

    if (refs == NULL) {
	refs = backrefs_new();
	g_object_set_data_full(G_OBJECT(buf), "backrefs", refs,
			       backrefs_destroy);
    }

    if (refs->n_marks == refs->n_slots) {
	/* @refs is full already */
	int ns = refs->n_slots + 1;

	refs->marks = realloc(refs->marks, ns * sizeof *refs->marks);
	refs->n_slots = ns;
    }

    refs->marks[refs->n_marks] = mark;
    refs->n_marks += 1;
}

/* Pop a GtkTextMark off the stack of backward references for @buf, if
   such a stack exists and is not empty.
*/

static GtkTextMark *pop_backref (GtkTextBuffer *buf)
{
    backrefs *refs = g_object_get_data(G_OBJECT(buf), "backrefs");
    GtkTextMark *ret = NULL;

    if (refs != NULL && refs->n_marks > 0) {
	int n = refs->n_marks - 1;

	ret = refs->marks[n];
	refs->marks[n] = NULL;
	refs->n_marks = n;
    }

    return ret;
}

/* Add a mark at the current insertion point and push it
   onto the stack.
*/

static void textbuf_set_backref (GtkTextBuffer *buf)
{
    GtkTextIter point;
    GtkTextMark *mark;

    gtk_text_buffer_get_iter_at_mark(buf, &point,
				     gtk_text_buffer_get_insert(buf));
    mark = gtk_text_mark_new(NULL, FALSE);
    gtk_text_buffer_add_mark(buf, mark, &point);
    push_backref(buf, mark);
}

/* Respond to Alt-, in gretl_edit */

void textview_go_back (GtkTextView *tview)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(tview);
    GtkTextMark *mark = pop_backref(buf);

    if (mark != NULL) {
	GtkTextIter iter;

	gtk_text_buffer_get_iter_at_mark(buf, &iter, mark);
	gtk_text_buffer_place_cursor(buf, &iter);
	gtk_text_view_scroll_to_mark(tview, mark, 0.0, TRUE, 0, 0.1);
	gtk_text_buffer_delete_mark(buf, mark);
    }
}

/* Determine whether the GtkTextBuffer in @vwin has mark to
   which "Go back"
*/

int textview_has_backref (GtkTextView *tview)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(tview);
    backrefs *refs = g_object_get_data(G_OBJECT(buf), "backrefs");

    return (refs != NULL && refs->n_marks > 0);
}

static void find_function_def (GtkTextView *tview,
			       gchar *sigstart)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, match;
    gboolean found;

    tbuf = gtk_text_view_get_buffer(tview);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    found = gtk_text_iter_forward_search(&start, sigstart,
					 GTK_TEXT_SEARCH_TEXT_ONLY,
					 &match, NULL, NULL);
    if (found) {
	GtkTextMark *targ;

	/* first set a mark for going back */
	textbuf_set_backref(tbuf);

	/* then move to the function definition */
	gtk_text_buffer_place_cursor(tbuf, &match);
	targ = gtk_text_buffer_create_mark(tbuf, "targ", &match, FALSE);
	gtk_text_view_scroll_to_mark(tview, targ, 0.05, FALSE, 0, 0);
    }
}

void alt_dot_find (GtkTextView *tview)
{
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tview);
    int role = FUNC_HELP;
    gchar *id = NULL;

    id = get_identifier_at_cursor(tbuf, &role);

    if (id != NULL && *id != '\0') {
	ufunc *uf = get_user_function_by_name(id);

	if (uf == NULL) {
	    warnbox(_("Function was not found"));
	} else {
	    GretlType t = user_func_get_return_type(uf);
	    const char *tstr = gretl_type_get_name(t);
	    gchar *needle;

	    needle = g_strdup_printf("function %s %s", tstr, id);
	    find_function_def(tview, needle);
	    g_free(needle);
	}
    }

    g_free(id);
}

void find_funcdef_callback (GtkWidget *w, gpointer data)
{
    gchar *needle = g_object_get_data(G_OBJECT(w), "needle");
    windata_t *vwin = g_object_get_data(G_OBJECT(w), "searchwin");

    find_function_def(GTK_TEXT_VIEW(vwin->text), needle);
    gtk_widget_destroy(gtk_widget_get_toplevel(w));
}
