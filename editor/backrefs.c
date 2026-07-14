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
#include "gretl_edit.h"
#include "textbuf.h"
#include "load_functions.h"
#include "gretl_func.h"
#include "gretl_typemap.h"
#include "backrefs.h"

/* Apparatus to support multi-level "Go back" references in a
   GtkTextView object.
*/

typedef struct backrefs_ {
    GtkNotebook *book;
    GtkTextMark **marks;
    int *pages;
    int n_marks;
    int n_slots;
} backrefs;

static backrefs *editor_refs;

/* Allocate a backrefs struct with space for one reference. */

static backrefs *backrefs_new (void)
{
    backrefs *refs = malloc(sizeof *refs);

    refs->marks = malloc(sizeof *refs->marks);
    refs->pages = malloc(sizeof *refs->pages);
    refs->n_marks = 0;
    refs->n_slots = 1;

    return refs;
}

#if 0 /* In case we decide a cleanup function is needed */

/* Free the backrefs. */

static void backrefs_destroy (void)
{
    if (editor_refs != NULL) {
	free(editor_refs->marks);
	free(editor_refs->pages);
	free(editor_refs);
	editor_refs = NULL;
    }
}

#endif

/* Push @mark onto the stack of backward references. */

static void push_backref (GtkTextMark *mark, int page)
{
    backrefs *refs = editor_refs;

    if (refs == NULL) {
	refs = editor_refs = backrefs_new();
    }

    if (refs->n_marks == refs->n_slots) {
	/* @refs is full already */
	int ns = refs->n_slots + 1;

	refs->marks = realloc(refs->marks, ns * sizeof *refs->marks);
	refs->pages = realloc(refs->pages, ns * sizeof *refs->pages);
	refs->n_slots = ns;
    }

    refs->marks[refs->n_marks] = mark;
    refs->pages[refs->n_marks] = page;
    refs->n_marks += 1;
}

/* Pop a GtkTextMark off the stack of backward references, if such a
   stack exists and is not empty.
*/

static GtkTextMark *pop_backref (int *page)
{
    backrefs *refs = editor_refs;
    GtkTextMark *mark = NULL;

    if (refs != NULL && refs->n_marks > 0) {
	int n = refs->n_marks - 1;

	*page = refs->pages[n];
	mark = refs->marks[n];
	refs->marks[n] = NULL;
	refs->n_marks = n;
    }

    return mark;
}

/* Add a GtkTextMark at the current insertion point and push it. */

static void textbuf_set_backref (GtkTextBuffer *tbuf, int page)
{
    GtkTextIter point;
    GtkTextMark *mark;

    gtk_text_buffer_get_iter_at_mark(tbuf, &point,
				     gtk_text_buffer_get_insert(tbuf));
    mark = gtk_text_mark_new(NULL, FALSE);
    gtk_text_buffer_add_mark(tbuf, mark, &point);

    push_backref(mark, page);
}

static GtkTextMark *back_mark;

static gboolean scroll_mark_onscreen_idle (gpointer data)
{
    GtkTextView *tview = GTK_TEXT_VIEW(data);
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tview);

    if (back_mark != NULL) {
        gtk_text_view_scroll_mark_onscreen(tview, back_mark);
	gtk_text_buffer_delete_mark(tbuf, back_mark);
	back_mark = NULL;
    }

    return G_SOURCE_REMOVE;
}

/* Respond to Alt-comma : return to the point at which the last search
   for a function definition was initiated.
*/

void editor_go_back (void)
{
    int page = 0;

    back_mark = pop_backref(&page);

    if (back_mark != NULL) {
	GtkNotebook *book;
	windata_t *vwin;
	GtkTextView *tview;
	GtkTextBuffer *tbuf;
	GtkTextIter iter;
	GtkWidget *w;
	int cp;

	book = GTK_NOTEBOOK(get_notebook());
	cp = gtk_notebook_get_current_page(book);
	if (cp != page) {
	    /* switch back to the page containing @mark */
	    gtk_notebook_set_current_page(book, page);
	}
	w = gtk_notebook_get_nth_page(book, page);
	vwin = g_object_get_data(G_OBJECT(w), "vwin");
	tview = GTK_TEXT_VIEW(vwin->text);
	tbuf = gtk_text_view_get_buffer(tview);

	gtk_text_buffer_get_iter_at_mark(tbuf, &iter, back_mark);
	gtk_text_buffer_place_cursor(tbuf, &iter);
	g_idle_add(scroll_mark_onscreen_idle, tview);
    }
}

/* Check whether there's a previously set mark to which we can return */

int editor_has_backref (void)
{
    backrefs *refs = editor_refs;

    return (refs != NULL && refs->n_marks > 0);
}

static gboolean find_function_def (windata_t *target,
				   const gchar *sigstart,
				   GtkNotebook *book,
				   int page,
				   windata_t *src)
{
    GtkTextView *tview;
    GtkTextBuffer *tbuf;
    GtkTextIter start, match;
    gboolean found = FALSE;

    maybe_load_functions(target);
    tview = GTK_TEXT_VIEW(target->text);
    tbuf = gtk_text_view_get_buffer(tview);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    found = gtk_text_iter_forward_search(&start, sigstart,
					 GTK_TEXT_SEARCH_TEXT_ONLY,
					 &match, NULL, NULL);
    if (found) {
	int cp = gtk_notebook_get_current_page(book);
	GtkTextView *orig_view = GTK_TEXT_VIEW(src->text);
	GtkTextBuffer *orig_buf = gtk_text_view_get_buffer(orig_view);
	GtkTextMark *targ;

	/* first set a mark for going back */
	textbuf_set_backref(orig_buf, cp);

	/* switch the page, if necessary */
	if (page != cp) {
	    gtk_notebook_set_current_page(book, page);
	}

	/* then move to the function definition */
	gtk_text_buffer_place_cursor(tbuf, &match);
	targ = gtk_text_buffer_create_mark(tbuf, "targ", &match, FALSE);
	gtk_text_view_scroll_mark_onscreen(tview, targ);
    }

    return found;
}

static gboolean find_in_notebook (windata_t *vwin,
				  const gchar *needle)
{
    GtkNotebook *book = GTK_NOTEBOOK(get_notebook());
    int np = gtk_notebook_get_n_pages(book);
    GtkWidget *tab;
    windata_t *target;
    gboolean found = FALSE;
    int i;

    for (i=0; i<np && !found; i++) {
	tab = gtk_notebook_get_nth_page(book, i);
	target = g_object_get_data(G_OBJECT(tab), "vwin");
	if (target->role == EDIT_HANSL) {
	    found = find_function_def(target, needle, book, i, vwin);
	    if (found) {
		fprintf(stderr, "found '%s' on nb page %d\n", needle, i);
	    }
	}
    }

    return found;
}

static void find_funcdef_callback (GtkWidget *w, gpointer data)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(w), "vwin");
    gchar *needle = g_object_get_data(G_OBJECT(w), "needle");

    find_in_notebook(vwin, needle);
    gtk_widget_destroy(gtk_widget_get_toplevel(w));
}

static gchar *get_alt_dot_needle (const gchar *id)
{
    ufunc *uf = get_user_function_by_name(id);
    gchar *needle = NULL;

    if (uf == NULL) {
	warnbox(_("Function was not found"));
    } else {
	GretlType t = user_func_get_return_type(uf);
	const char *tstr = gretl_type_get_name(t);

	needle = g_strdup_printf("function %s %s", tstr, id);
    }

    return needle;
}

void alt_dot_find (windata_t *vwin)
{
    GtkTextView *tview = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tview);
    int role = FUNC_HELP;
    gchar *id = NULL;

    maybe_load_functions(vwin);
    id = get_identifier_at_cursor(tbuf, &role);

    if (id != NULL && *id != '\0') {
	gchar *needle = get_alt_dot_needle(id);

	if (needle != NULL) {
	    find_in_notebook(vwin, needle);
	    g_free(needle);
	}
    }

    g_free(id);
}

/* Called from viewers.c, appending to a pop-up window showing the
   signature of a function.
*/

void add_funcdef_finder (const char *sig,
			 windata_t *popwin,
			 windata_t *vwin)
{
    const gchar *p = strchr(sig, '(');
    gchar *needle = g_strndup(sig, p - sig);
    GtkWidget *ebox = gtk_event_box_new();
    GtkWidget *label;
    gchar *fmt = NULL;
    gchar *buf = NULL;

    /* create an event box with a label inside */
    gtk_event_box_set_visible_window(GTK_EVENT_BOX(ebox), FALSE);
    gtk_widget_set_can_default(ebox, FALSE);
    fmt = g_strdup_printf("<span color=\"%s\">%%s</span>", blue_for_text());
    buf = g_markup_printf_escaped(fmt, _("Find definition"));
    label = gtk_label_new(NULL);
    gtk_label_set_markup(GTK_LABEL(label), buf);
    gtk_label_set_justify(GTK_LABEL(label), GTK_JUSTIFY_CENTER);
    // set_popup_bg(label);
    gtk_container_add(GTK_CONTAINER(ebox), label);
    g_free(buf);
    g_free(fmt);

    /* pack the event box into @popwin and connect signals */
    gtk_box_pack_start(GTK_BOX(popwin->vbox), ebox, FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(ebox), "vwin", vwin);
    g_object_set_data_full(G_OBJECT(ebox), "needle", needle, g_free);
    g_signal_connect(ebox, "button-release-event",
		     G_CALLBACK(find_funcdef_callback), NULL);
    g_signal_connect(ebox, "enter-notify-event",
		     G_CALLBACK(show_link_cursor), NULL);
    g_signal_connect(ebox, "leave-notify-event",
		     G_CALLBACK(revert_cursor), NULL);
}
