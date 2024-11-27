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
#include "toolbar.h"
#include "dlgutils.h"
#include "winstack.h"
#include "tabwin.h"
#include "gui_recode.h"
#include "gretl_func.h"
#include "addons_utils.h"

#ifndef GRETL_EDIT
#include "library.h"
#include "guiprint.h"
#include "datafiles.h"
#include "database.h"
#include "fncall.h"
#endif

#ifdef G_OS_WIN32
# include "gretlwin32.h" /* for browser_open() */
#endif

#if GTKSOURCEVIEW_VERSION > 2
# define GTK_IS_SOURCE_VIEW GTK_SOURCE_IS_VIEW
#else /* using GtkSourceView 2 */
# include <gtksourceview/gtksourcelanguagemanager.h>
# include <gtksourceview/gtksourceprintcompositor.h>
# include <gtksourceview/gtksourcestyleschememanager.h>
#endif

#ifdef HAVE_GTKSV_COMPLETION
# include "completions.h"
#endif

#define TABDEBUG 0
#define KDEBUG 0

/* Dummy "page" numbers for use in hyperlinks: these
   must be greater than the number of gretl commands
   and built-in functions to avoid collisions.
*/

#define GUIDE_PAGE  999
#define SCRIPT_PAGE 998
#define GFR_PAGE    997
#define BIB_PAGE    996
#define EXT_PAGE    995
#define PDF_PAGE    994
#define MNU_PAGE    993
#define DBN_PAGE    992
#define DBS_PAGE    991
#define NEXT_PAGE   990

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

#define gui_help(r) (r == GUI_HELP || r == GUI_HELP_EN)
#define function_help(r) (r == FUNC_HELP || r == FUNC_HELP_EN)
#define foreign_script_role(r) (r == EDIT_GP || \
				r == EDIT_R || \
				r == EDIT_OX || \
				r == EDIT_OCTAVE || \
				r == EDIT_PYTHON || \
				r == EDIT_STATA ||  \
				r == EDIT_JULIA || \
				r == EDIT_DYNARE || \
				r == EDIT_LPSOLVE)
#define medium_markup(r) (r == VIEW_PKG_INFO || \
			  r == EDIT_PKG_HELP || \
			  r == EDIT_PKG_GHLP)

/* globals accessed in settings.c */
int tabwidth = 4;
int smarttab = 1;
int script_line_numbers = 0;
int script_auto_bracket = 0;

/* file-scope constant */
static const char *hidden_marker = "hidden region\n";

static gboolean script_electric_enter (windata_t *vwin, int alt);
static gboolean script_tab_handler (windata_t *vwin, GdkEvent *event);
static gboolean
script_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p);
static gboolean
insert_text_with_markup (GtkTextBuffer *tbuf, GtkTextIter *iter,
			 const char *s, int role);
static void connect_link_signals (windata_t *vwin);
static void auto_indent_script (GtkWidget *w, windata_t *vwin);
static int maybe_insert_smart_tab (windata_t *vwin, int *comp_ok);
#ifndef GRETL_EDIT
static gchar *textview_get_current_line_with_newline (GtkWidget *view);
#endif

static int object_get_int (gpointer p, const char *key)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), key));
}

static void object_set_int (gpointer p, const char *key, int k)
{
    g_object_set_data(G_OBJECT(p), key, GINT_TO_POINTER(k));
}

#define DELETE_DEBUG 0 /* activate if needed */

/* for use when pasting text obtained via the clipboard */

void text_delete_invisibles (gchar *s)
{
    gchar *p = s;
    gunichar gu;
    int n, m;

#if DELETE_DEBUG
    fprintf(stderr, "text_delete_invisibles 1 len %ld\n", strlen(s));
#endif
    while (*p != '\0') {
        p = g_utf8_find_next_char(p, NULL);
        gu = g_utf8_get_char(p);
        if (g_unichar_iszerowidth(gu)) {
            /* delete bytes representing a zero-width character */
            n = g_utf8_find_next_char(p, NULL) - p;
            m = strlen(p+n) + 1;
            memmove(p, p+n, m);
        }
    }
#if DELETE_DEBUG
    fprintf(stderr, "text_delete_invisibles 2 len %ld\n", strlen(s));
#endif
}

void text_paste (GtkWidget *w, windata_t *vwin)
{
    GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);
    gchar *src = gtk_clipboard_wait_for_text(cb);

    if (src != NULL) {
        text_delete_invisibles(src);
        textview_insert_text(vwin->text, src);
        g_free(src);
    }
}

/* replaced by the above, 2024-11-23

void text_paste (GtkWidget *w, windata_t *vwin)
{
    gtk_text_buffer_paste_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)),
				    gtk_clipboard_get(GDK_NONE),
				    NULL, TRUE);
}

*/

void text_set_cursor (GtkWidget *w, GdkCursorType cspec)
{
    static GdkCursor *question_cursor;
    GdkWindow *win = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
                                              GTK_TEXT_WINDOW_TEXT);

    if (cspec == 0) {
	gdk_window_set_cursor(win, NULL);
    } else if (cspec == GDK_QUESTION_ARROW) {
	if (question_cursor == NULL) {
	    question_cursor = gdk_cursor_new(GDK_QUESTION_ARROW);
	}
	gdk_window_set_cursor(win, question_cursor);
    } else {
	GdkCursor *cursor = gdk_cursor_new(cspec);

	if (cursor != NULL) {
	    gdk_window_set_cursor(win, cursor);
	    gdk_cursor_unref(cursor);
	}
    }
}

void cursor_to_top (windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter start;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_to_mark(view, mark, 0.0, FALSE, 0, 0);
    gtk_text_buffer_delete_mark(buf, mark);
}

void cursor_to_end (windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter end;

    gtk_text_buffer_get_end_iter(buf, &end);
    gtk_text_buffer_place_cursor(buf, &end);
}

void scroll_to_foot (windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter end;
    GtkTextMark *mark;

    gtk_text_buffer_get_end_iter(buf, &end);
    mark = gtk_text_buffer_create_mark(buf, NULL, &end, FALSE);
    gtk_text_view_scroll_to_mark(view, mark, 0.0, FALSE, 0, 0);
    gtk_text_buffer_delete_mark(buf, mark);
}

void cursor_to_mark (windata_t *vwin, GtkTextMark *mark)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter iter;

    gtk_text_buffer_get_iter_at_mark(buf, &iter, mark);
    gtk_text_buffer_place_cursor(buf, &iter);
    gtk_text_view_scroll_to_mark(view, mark, 0.0, TRUE, 0, 0.1);
}

void scroll_to_line (windata_t *vwin, int line)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter iter;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &iter);
    gtk_text_iter_forward_lines(&iter, line - 1);
    gtk_text_buffer_place_cursor(buf, &iter);
    mark = gtk_text_buffer_create_mark(buf, NULL, &iter, FALSE);
    gtk_text_view_scroll_to_mark(view, mark, 0.0, FALSE, 0, 0);
    gtk_text_buffer_delete_mark(buf, mark);
}

int textbuf_get_n_lines (windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkTextBuffer *buf = gtk_text_view_get_buffer(view);
    GtkTextIter end;

    gtk_text_buffer_get_end_iter(buf, &end);
    return gtk_text_iter_get_line(&end) + 1;
}

static void get_char_width_and_height (GtkWidget *widget,
				       int *width,
				       int *height)
{
    PangoContext *pc;
    PangoLayout *pl;
    int w = 0, h = 0;

    pc = gtk_widget_get_pango_context(widget);
    pango_context_set_font_description(pc, fixed_font);
    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, "X", -1);
    pango_layout_get_pixel_size(pl, &w, &h);
    g_object_unref(pl);

    if (width != NULL) {
	*width = w;
    }
    if (height != NULL) {
	*height = h;
    }
}

gint get_char_width (GtkWidget *widget)
{
    int width;

    get_char_width_and_height(widget, &width, NULL);
    return width;
}

gint get_char_height (GtkWidget *widget)
{
    int height;

    get_char_width_and_height(widget, NULL, &height);
    return height;
}

gchar *textview_get_text (GtkWidget *view)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);

    return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
}

gchar *textview_get_trimmed_text (GtkWidget *view)
{
    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    return g_strchug(g_strchomp(textview_get_text(view)));
}

static gchar *strip_hidden_placeholders (gchar *content)
{
    int len = strlen(hidden_marker);
    gchar *s = content;

    while ((s = strstr(s, hidden_marker)) != NULL) {
        shift_string_left(s, len);
        s += len;
    }

    return content;
}

/* Note: a non-zero value for @save means that we're grabbing the
   GtkTextBuffer content to save it to file; otherwise we want it for
   the purpose of execution. In the former case we include any
   currently invisible regions (by passing TRUE as the last argument
   to gtk_text_buffer_get_text); in the latter we'll exclude them.  In
   both cases we need to exclude any hidden region placeholders.
*/

gchar *textview_get_hansl (GtkWidget *view, int save)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    gboolean include_invisible;
    int n_hidden;
    gchar *content;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);

    n_hidden = object_get_int(tbuf, "n_hidden");
    include_invisible = (save && n_hidden > 0);
    content = gtk_text_buffer_get_text(tbuf, &start, &end, include_invisible);

    if (n_hidden > 0) {
        strip_hidden_placeholders(content);
    }

    return content;
}

static gchar *normalize_line (const char *s)
{
    int i, j = 0, n = strlen(s);
    gchar *ret = g_malloc0(n);

     for (i=0; s[i]; i++) {
	if (isspace(s[i])) {
	    ret[j++] = ' ';
	    while (isspace(s[i])) i++;
	}
        ret[j++] = s[i];
    }

    return ret;
}

gchar *textview_get_normalized_line (GtkWidget *view)
{
    gchar *ret = NULL;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    ret = g_strchug(g_strchomp(textview_get_text(view)));
    if (strchr(ret, '\n') != NULL || strstr(ret, "  ") != NULL) {
	gchar *alt = normalize_line(ret);

	g_free(ret);
	ret = alt;
    }

    return ret;
}

/* Special: handle the case where text has been line-wrapped
   in an editor window and we want to save the text as
   wrapped -- that is, forcibly to truncate excessively
   long lines. We use this for function-package help text.
*/

gchar *textview_get_wrapped_text (GtkWidget *view)
{
    GtkTextView *tview;
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    GString *str;
    gchar *line;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tview = GTK_TEXT_VIEW(view);
    tbuf = gtk_text_view_get_buffer(tview);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    end = start;

    /* first detect and handle the case where no wrapping has
       occurred */
    if (!gtk_text_view_forward_display_line(tview, &end)) {
	gtk_text_buffer_get_end_iter(tbuf, &end);
	return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    }

    str = g_string_new(NULL);

    end = start;
    while (gtk_text_view_forward_display_line(tview, &end)) {
	line = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	g_strchomp(line);
	g_string_append(str, line);
	g_string_append_c(str, '\n');
	g_free(line);
	start = end;
    }

    if (!gtk_text_iter_is_end(&start)) {
	/* there's some residual text */
	gtk_text_buffer_get_end_iter(tbuf, &end);
	line = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	g_strchomp(line);
	g_string_append(str, line);
	g_string_append_c(str, '\n');
	g_free(line);
    }

    return str == NULL ? NULL : g_string_free(str, FALSE);
}

gchar *textview_get_selection_or_all (GtkWidget *view,
				      gboolean *selection)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    gchar *content;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    if (tbuf == NULL) {
	return NULL;
    }

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	*selection = TRUE;
    } else {
	*selection = FALSE;
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
    }

    content = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

    if (*selection == FALSE && object_get_int(tbuf, "n_hidden") > 0) {
        strip_hidden_placeholders(content);
    }

    return content;
}

static int real_textview_set_text (GtkWidget *view,
				   const gchar *text,
				   gboolean select)
{
    GtkTextBuffer *tbuf;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), 1);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    g_return_val_if_fail(tbuf != NULL, 1);

    if (text != NULL && select) {
	GtkTextIter start, end;

	gtk_text_buffer_set_text(tbuf, text, -1);
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
	gtk_text_buffer_select_range(tbuf, &start, &end);
    } else if (text != NULL) {
	gtk_text_buffer_set_text(tbuf, text, -1);
    } else {
	gtk_text_buffer_set_text(tbuf, "", -1);
    }

    return 0;
}

int textview_set_text (GtkWidget *view, const gchar *text)
{
    return real_textview_set_text(view, text, FALSE);
}

int textview_set_text_selected (GtkWidget *view, const gchar *text)
{
    return real_textview_set_text(view, text, TRUE);
}

int textview_set_cursor_at_line (GtkWidget *view, int line)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), 1);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    g_return_val_if_fail(tbuf != NULL, 1);

    gtk_text_buffer_get_iter_at_line(tbuf, &iter, line);
    gtk_text_buffer_place_cursor(tbuf, &iter);

    return 0;
}

int viewer_char_count (windata_t *vwin)
{
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    return gtk_text_buffer_get_char_count(tbuf);
}

void text_redo (GtkWidget *w, windata_t *vwin)
{
    if (vwin->sbuf != NULL && gtk_source_buffer_can_redo(vwin->sbuf)) {
	gtk_source_buffer_redo(vwin->sbuf);
    } else {
	warnbox(_("No redo information available"));
    }
}

void text_undo (GtkWidget *w, windata_t *vwin)
{
    if (vwin->sbuf != NULL && gtk_source_buffer_can_undo(vwin->sbuf)) {
	gtk_source_buffer_undo(vwin->sbuf);
    } else {
	warnbox(_("No undo information available"));
    }
}

int text_can_undo (windata_t *vwin)
{
    if (vwin->sbuf != NULL) {
	return gtk_source_buffer_can_undo(vwin->sbuf);
    } else {
	return 0;
    }
}

static int source_buffer_load_file (GtkSourceBuffer *sbuf,
				    int role,
				    const char *fname)
{
    GtkTextBuffer *tbuf = GTK_TEXT_BUFFER(sbuf);
    GtkTextIter iter;
    gchar *buf = NULL;
    gsize sz = 0;

    gtk_source_buffer_begin_not_undoable_action(sbuf);
    gtk_text_buffer_set_text(tbuf, "", -1);
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    gretl_file_get_contents(fname, &buf, &sz);

    if (buf != NULL) {
	gchar *trbuf = NULL;

	if (!g_utf8_validate(buf, -1, NULL)) {
	    trbuf = my_locale_to_utf8(buf);
	    if (trbuf != NULL) {
		gtk_text_buffer_insert(tbuf, &iter, trbuf, -1);
		g_free(trbuf);
	    }
	} else {
	    gtk_text_buffer_insert(tbuf, &iter, buf, -1);
	}
	g_free(buf);
    }

    gtk_source_buffer_end_not_undoable_action(sbuf);
    gtk_text_buffer_set_modified(tbuf, role == EDIT_PKG_SAMPLE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(tbuf, &iter);
    gtk_text_buffer_place_cursor(tbuf, &iter);

    return 0;
}

static int source_buffer_load_buf (GtkSourceBuffer *sbuf, const char *buf)
{
    GtkTextBuffer *tbuf = GTK_TEXT_BUFFER(sbuf);
    GtkTextIter iter;

    gtk_source_buffer_begin_not_undoable_action(sbuf);
    gtk_text_buffer_set_text(tbuf, buf, -1);
    gtk_source_buffer_end_not_undoable_action(sbuf);
    gtk_text_buffer_set_modified(tbuf, FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(tbuf, &iter);
    gtk_text_buffer_place_cursor(tbuf, &iter);

    return 0;
}

struct sv_lang {
    int role;
    const char *str;
};

static void sourceview_apply_language (windata_t *vwin)
{
    struct sv_lang sv_langs[] = {
	{ EDIT_GP,      "gnuplot" },
	{ EDIT_R,       "r" },
	{ EDIT_OX,      "cpp" },
	{ EDIT_OCTAVE,  "octave" },
	{ EDIT_PYTHON,  "python" },
	{ EDIT_JULIA,   "julia" },
	{ EDIT_STATA,   "stata" },
	{ EDIT_DYNARE,  "cpp" },
	{ EDIT_LPSOLVE, "lpsolve" },
	{ EDIT_SPEC,    "gfnspec" }
    };
    GtkSourceLanguageManager *lm;
    GtkSourceLanguage *lang = NULL;
    const char *lname = "gretl";
    int i;

    lm = g_object_get_data(G_OBJECT(vwin->sbuf), "languages-manager");
    if (lm == NULL) {
	return;
    }

    for (i=0; i<G_N_ELEMENTS(sv_langs); i++) {
	if (vwin->role == sv_langs[i].role) {
	    lname = sv_langs[i].str;
	}
    }

    lang = gtk_source_language_manager_get_language(lm, lname);

    if (lang == NULL) {
	const gchar * const *S = gtk_source_language_manager_get_search_path(lm);

	fprintf(stderr, "*** gtksourceview: lang is NULL for '%s'\n", lname);
	if (S != NULL) {
	    fprintf(stderr, "the gtksourceview search path:\n");
	    while (*S != NULL) {
		fprintf(stderr, " %s\n", *S);
		S++;
	    }
	} else {
	    fprintf(stderr, "the gtksourceview search path is NULL!\n");
	}
    } else {
	gtk_source_buffer_set_language(vwin->sbuf, lang);
    }
}

#define SV_PRINT_DEBUG 0

#if SV_PRINT_DEBUG

static void show_print_context (GtkPrintContext *context)
{
    gdouble x = gtk_print_context_get_width(context);
    gdouble y = gtk_print_context_get_height(context);
    GtkPageSetup *psu;

    fprintf(stderr, "begin_print: context pixel size: %g x %g\n", x, y);
    x = gtk_print_context_get_dpi_x(context);
    y = gtk_print_context_get_dpi_y(context);
    fprintf(stderr, "  context dpi: %g, %g\n", x, y);

    psu = gtk_print_context_get_page_setup(context);
    if (psu != NULL) {
	x = gtk_page_setup_get_paper_width(psu, GTK_UNIT_INCH);
	y = gtk_page_setup_get_paper_width(psu, GTK_UNIT_POINTS);
	fprintf(stderr, "  paper width: %g in, %g points\n", x, y);
	x = gtk_page_setup_get_paper_height(psu, GTK_UNIT_INCH);
	y = gtk_page_setup_get_paper_height(psu, GTK_UNIT_POINTS);
	fprintf(stderr, "  paper height: %g in, %g points\n", x, y);
    }
}

#endif

static void begin_print (GtkPrintOperation *operation,
                         GtkPrintContext *context,
                         gpointer data)
{
    GtkSourcePrintCompositor *comp;
    int n_pages;

    comp = GTK_SOURCE_PRINT_COMPOSITOR(data);

#if SV_PRINT_DEBUG
    GtkPrintSettings *st = gtk_print_operation_get_print_settings(operation);
    if (st != NULL) {
	int rx = gtk_print_settings_get_resolution_x(st);
	int ry = gtk_print_settings_get_resolution_y(st);

        fprintf(stderr, "settings: resolution %d,%d\n", rx, ry);
    }
    show_print_context(context);
#endif

    while (!gtk_source_print_compositor_paginate(comp, context));

    n_pages = gtk_source_print_compositor_get_n_pages(comp);
    gtk_print_operation_set_n_pages(operation, n_pages);
}

static void draw_page (GtkPrintOperation *operation,
		       GtkPrintContext *context,
		       gint page_nr,
		       gpointer data)
{
    GtkSourcePrintCompositor *comp =
	GTK_SOURCE_PRINT_COMPOSITOR(data);

    gtk_source_print_compositor_draw_page(comp, context,
					  page_nr);
}

void sourceview_print (windata_t *vwin)
{
    GtkSourceView *view = GTK_SOURCE_VIEW(vwin->text);
    GtkSourcePrintCompositor *comp;
    GtkPrintOperation *print;
    GtkPrintOperationResult res;
    GtkWidget *mainwin;
    GError *error = NULL;

    comp = gtk_source_print_compositor_new_from_view(view);
    print = gtk_print_operation_new();

#ifdef G_OS_WIN32
    /* the units are wacky if we don't set this */
    gtk_print_operation_set_unit(print, GTK_UNIT_POINTS);
#endif

    gtk_source_print_compositor_set_right_margin(comp, 60, GTK_UNIT_POINTS);
    gtk_source_print_compositor_set_left_margin(comp, 60, GTK_UNIT_POINTS);
    gtk_source_print_compositor_set_top_margin(comp, 54, GTK_UNIT_POINTS);
    gtk_source_print_compositor_set_bottom_margin(comp, 72, GTK_UNIT_POINTS);
    gtk_source_print_compositor_set_wrap_mode(comp, GTK_WRAP_WORD);
    gtk_source_print_compositor_set_body_font_name(comp, "Monospace 9");

    g_signal_connect(G_OBJECT(print), "begin_print", G_CALLBACK(begin_print), comp);
    g_signal_connect(G_OBJECT(print), "draw_page", G_CALLBACK(draw_page), comp);

    mainwin = vwin_toplevel(vwin);
    res = gtk_print_operation_run(print,
				  GTK_PRINT_OPERATION_ACTION_PRINT_DIALOG,
				  GTK_WINDOW(mainwin),
				  &error);

    if (res == GTK_PRINT_OPERATION_RESULT_ERROR) {
	GtkWidget *dlg;

	dlg = gtk_message_dialog_new(GTK_WINDOW(mainwin),
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     GTK_MESSAGE_ERROR,
				     GTK_BUTTONS_CLOSE,
				     "Error printing file:\n%s",
				     error->message);
	g_signal_connect(G_OBJECT(dlg), "response",
			 G_CALLBACK(gtk_widget_destroy), NULL);
	gtk_widget_show(dlg);
	g_error_free(error);
    } else if (res == GTK_PRINT_OPERATION_RESULT_APPLY) {
	; /* OK: maybe save the settings? */
    }

    if (print != NULL) {
	g_object_unref(print);
    }
}

void sourceview_insert_file (windata_t *vwin, const char *fname)
{
    sourceview_apply_language(vwin);
    if (fname != NULL) {
	source_buffer_load_file(vwin->sbuf, vwin->role, fname);
    }
}

void sourceview_insert_buffer (windata_t *vwin, const char *buf)
{
    sourceview_apply_language(vwin);
    source_buffer_load_buf(vwin->sbuf, buf);
}

static void set_source_tabs (GtkWidget *w, int cw)
{
    static PangoTabArray *ta;

    if (ta == NULL) {
	int tabw = tabwidth * cw;
	gint i, loc = tabw;

	ta = pango_tab_array_new(10, TRUE);
	for (i=0; i<10; i++) {
	    pango_tab_array_set_tab(ta, i, PANGO_TAB_LEFT, loc);
	    loc += tabw;
	}
    }

    gtk_text_view_set_tabs(GTK_TEXT_VIEW(w), ta);
}

#define tabkey(k) (k == GDK_Tab || \
		   k == GDK_ISO_Left_Tab || \
		   k == GDK_KP_Tab)

/* Special keystrokes in native script window: Ctrl-Return sends the
   current line for execution; Ctrl-R sends the whole script for
   execution (keyboard equivalent of the "execute" button);
   Ctrl-I does auto-indentation.
*/

static gint script_key_handler (GtkWidget *w,
				GdkEvent *event,
				windata_t *vwin)
{
    guint keyval = ((GdkEventKey *) event)->keyval;
    guint state = ((GdkEventKey *) event)->state;
    gboolean ret = FALSE;

#if KDEBUG
    fprintf(stderr, "HERE script_key_handler (keyval %u, %s)\n",
	    keyval, gdk_keyval_name(keyval));
#endif

    if (state & GDK_CONTROL_MASK) {
	if (keyval == GDK_R) {
	    /* Ctrl-Shift-r */
#ifdef GRETL_EDIT
	    do_run_script(w, vwin);
#else
	    run_script_silent(w, vwin);
#endif
	    ret = TRUE;
	} else if (keyval == GDK_r) {
	    /* plain Ctrl-r */
	    do_run_script(w, vwin);
	    ret = TRUE;
        } else if (keyval == GDK_v) {
            text_paste(w, vwin);
            ret = TRUE;
        }
#ifndef GRETL_EDIT
	else if (keyval == GDK_Return) {
	    gchar *str = textview_get_current_line_with_newline(w);

	    if (str != NULL) {
		if (!string_is_blank(str)) {
		    run_script_fragment(vwin, str);
		}
		g_free(str);
	    }
	    ret = TRUE;
	}
#endif
	else if (keyval == GDK_i) {
	    auto_indent_script(w, vwin);
	    ret = TRUE;
	}
    } else if (keyval == GDK_F1) {
	set_window_help_active(vwin);
	interactive_script_help(NULL, NULL, vwin);
	ret = TRUE;
    } else if (editing_hansl(vwin->role)) {
	if (keyval == GDK_Return) {
	    ret = script_electric_enter(vwin, state & GDK_MOD1_MASK);
	} else if (tabkey(keyval)) {
#if TABDEBUG
	    fprintf(stderr, "*** calling script_tab_handler ***\n");
#endif
	    ret = script_tab_handler(vwin, event);
	} else if (script_auto_bracket && lbracket(keyval)) {
	    ret = script_bracket_handler(vwin, keyval);
	}
    }

    return ret;
}

static gint
foreign_script_key_handler (GtkWidget *w, GdkEvent *event, windata_t *vwin)
{
    guint keyval = ((GdkEventKey *) event)->keyval;
    gboolean ret = FALSE;

    if (((GdkEventKey *) event)->state & GDK_CONTROL_MASK) {
	if (keyval == GDK_r)  {
	    do_run_script(w, vwin);
	    ret = TRUE;
	}
    }

    return ret;
}

#ifdef PKGBUILD

# ifdef G_OS_WIN32

static gchar *ensure_utf8_path (gchar *path)
{
    if (!g_utf8_validate(path, -1, NULL)) {
	gchar *tmp;
	gsize bytes;

	tmp = g_locale_to_utf8(path, -1, NULL, &bytes, NULL);
	if (tmp != NULL) {
	    g_free(path);
	    path = tmp;
	}
    }

    return path;
}

# endif

/* Packages for Windows and macOS: gtksourceview needs to
   be told where to find its language-specs and style
   files: these live under gtksourceview inside the package.

   On Windows we need to ensure that the "set_search_path"
   functions are fed a UTF-8 path, since gtksourceview uses
   g_open() internally and the Glib filename encoding is
   always UTF-8 on Windows.
*/

static void ensure_sourceview_path (GtkSourceLanguageManager *lm)
{
    static int done;

    if (!done && lm == NULL) {
	lm = gtk_source_language_manager_get_default();
    }

    if (!done && lm != NULL) {
	GtkSourceStyleSchemeManager *mgr;
	gchar *dirs[2] = {NULL, NULL};

	dirs[0] = g_strdup_printf("%sgtksourceview", gretl_home());
# ifdef G_OS_WIN32
	dirs[0] = ensure_utf8_path(dirs[0]);
# endif
	gtk_source_language_manager_set_search_path(lm, dirs);

	mgr = gtk_source_style_scheme_manager_get_default();
	gtk_source_style_scheme_manager_set_search_path(mgr, dirs);
	gtk_source_style_scheme_manager_force_rescan(mgr);

	g_free(dirs[0]);
	done = 1;
    }
}

#else /* not PKGBUILD */

/* gtksourceview needs to be told to search both its own "native"
   paths and @prefix/share/gretl/gtksourceview for both language
   file and style files
*/

static void ensure_sourceview_path (GtkSourceLanguageManager *lm)
{
    static int done;

    if (!done && lm == NULL) {
	lm = gtk_source_language_manager_get_default();
    }

    if (!done && lm != NULL) {
	GtkSourceStyleSchemeManager *mgr;
	gchar *dirs[3] = {NULL, NULL, NULL};

	/* languages: we need to set path, can't just append */
	dirs[0] = g_strdup_printf("%sgtksourceview", gretl_home());
#if GTKSOURCEVIEW_VERSION > 3
	dirs[1] = g_strdup_printf("%s/share/gtksourceview-%d/language-specs",
				  SVPREFIX, GTKSOURCEVIEW_VERSION);
#else
	dirs[1] = g_strdup_printf("%s/share/gtksourceview-%d.0/language-specs",
				  SVPREFIX, GTKSOURCEVIEW_VERSION);
#endif
	gtk_source_language_manager_set_search_path(lm, dirs);

	/* styles: we can just append to the default path */
	mgr = gtk_source_style_scheme_manager_get_default();
	gtk_source_style_scheme_manager_append_search_path(mgr, dirs[0]);
	gtk_source_style_scheme_manager_force_rescan(mgr);

	g_free(dirs[0]);
	g_free(dirs[1]);

	done = 1;
    }
}

#endif

static void set_console_output_style (GtkSourceBuffer *sbuf,
				      GtkSourceStyleScheme *scheme)
{
    GtkSourceStyle *style = NULL;
    GtkTextTag *tag = NULL;
    GtkTextTagTable *tt;
    int done = 0;

    tt = gtk_text_buffer_get_tag_table(GTK_TEXT_BUFFER(sbuf));
    if (tt != NULL) {
	tag = gtk_text_tag_table_lookup(tt, "output");
    }
    if (tag != NULL) {
	style = gtk_source_style_scheme_get_style(scheme, "text");
    }
    if (style != NULL) {
	gchar *fg = NULL, *bg = NULL;

	g_object_get(style, "foreground", &fg, "background", &bg, NULL);
	if (fg != NULL && bg != NULL) {
	    g_object_set(tag, "foreground", fg, "background", bg, NULL);
	    done = 1;
	}
	g_free(fg);
	g_free(bg);
    }
    if (tag != NULL && !done) {
	/* fallback */
	if (dark_theme_active()) {
	    g_object_set(tag, "foreground", "#d3d7cf", "background", "#242424", NULL);
	} else {
	    g_object_set(tag, "foreground", "black", "background", "white", NULL);
	}
    }
}

static void set_style_for_buffer (GtkSourceBuffer *sbuf,
				  const char *id,
				  int role)
{
    GtkSourceStyleSchemeManager *mgr;
    GtkSourceStyleScheme *scheme;

    if (id == NULL || *id == '\0') {
	return;
    }

    mgr = gtk_source_style_scheme_manager_get_default();
    scheme = gtk_source_style_scheme_manager_get_scheme(mgr, id);
    if (scheme != NULL) {
	gtk_source_buffer_set_style_scheme(sbuf, scheme);
	if (role == CONSOLE) {
	    set_console_output_style(sbuf, scheme);
	}
    }
}

void set_style_for_textview (GtkWidget *text, const char *id)
{
    GtkSourceStyleSchemeManager *mgr;
    GtkSourceStyleScheme *scheme;
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text));
    mgr = gtk_source_style_scheme_manager_get_default();
    scheme = gtk_source_style_scheme_manager_get_scheme(mgr, id);
    if (scheme != NULL) {
	gtk_source_buffer_set_style_scheme(GTK_SOURCE_BUFFER(tbuf), scheme);
    }
}

#define gretl_script_role(r) (r == EDIT_HANSL || \
			      r == VIEW_SCRIPT || \
			      r == EDIT_PKG_CODE || \
			      r == EDIT_PKG_SAMPLE || \
			      r == VIEW_PKG_SAMPLE || \
                              r == VIEW_PKG_CODE)

void create_source (windata_t *vwin, int hsize, int vsize,
		    gboolean editable)
{
    GtkSourceLanguageManager *lm = NULL;
    GtkSourceBuffer *sbuf;
    GtkTextView *view;
    int cw;

    if (textview_use_highlighting(vwin->role)) {
	lm = gtk_source_language_manager_get_default();
	ensure_sourceview_path(lm);
    }

    if (editing_hansl(vwin->role)) {
	GtkTextTagTable *table = gtk_text_tag_table_new();
	GtkTextTag *htag = gtk_text_tag_new("hidden");

	g_object_set(htag, "invisible", 1, NULL);
	gtk_text_tag_table_add(table, htag);
	sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(table));
    } else {
	sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(NULL));
    }

    if (lm != NULL) {
	g_object_set_data(G_OBJECT(sbuf), "languages-manager", lm);
    }
    gtk_source_buffer_set_highlight_matching_brackets(sbuf, TRUE);

    vwin->text = gtk_source_view_new_with_buffer(sbuf);
    vwin->sbuf = sbuf;

#ifdef HAVE_GTKSV_COMPLETION
    if (editing_hansl(vwin->role) || vwin->role == CONSOLE) {
	set_sv_completion(vwin);
    }
#endif

    view = GTK_TEXT_VIEW(vwin->text);
    gtk_text_view_set_wrap_mode(view, GTK_WRAP_NONE);
    gtk_text_view_set_left_margin(view, 4);
    gtk_text_view_set_right_margin(view, 4);
#if GTK_MAJOR_VERSION > 2
    gtk_text_view_set_bottom_margin(view, 10);
#endif

    gtk_widget_modify_font(GTK_WIDGET(vwin->text), fixed_font);

    cw = get_char_width(vwin->text);
    set_source_tabs(vwin->text, cw);

    if (hsize > 0) {
	hsize *= cw;
	hsize += 48;
    }

    if (!(vwin->flags & VWIN_SWALLOW) && hsize > 0 && vsize > 0) {
	GtkWidget *vmain = vwin_toplevel(vwin);

	if (window_is_tab(vwin)) {
	    vsize += 15;
	}
	if (vsize < 0.62 * hsize) {
	    /* approx golden ratio */
	    vsize = 0.62 * hsize;
	}
	gtk_window_set_default_size(GTK_WINDOW(vmain), hsize, vsize);
    }

    gtk_text_view_set_editable(view, editable);
    gtk_text_view_set_cursor_visible(view, editable);

    if (vwin->role != EDIT_HEADER) {
	gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(vwin->text),
					      script_line_numbers);
    }

    if (lm != NULL) {
	set_style_for_buffer(sbuf, get_sourceview_style(), vwin->role);
    }

    if (!(vwin->flags & WVIN_KEY_SIGNAL_SET)) {
	g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
			 G_CALLBACK(catch_viewer_key), vwin);
    }

    if (gretl_script_role(vwin->role)) {
	g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
			 G_CALLBACK(script_key_handler), vwin);
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(script_popup_handler),
			 vwin);
	g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
			 G_CALLBACK(interactive_script_help), vwin);
    } else if (foreign_script_role(vwin->role)) {
	g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
			 G_CALLBACK(foreign_script_key_handler), vwin);
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(script_popup_handler),
			 vwin);
    } else if (vwin->role == VIEW_LOG) {
	g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
			 G_CALLBACK(interactive_script_help), vwin);
    }

    if (editing_hansl(vwin->role)) {
	connect_link_signals(vwin);
    }
}

/* Manufacture a little sampler sourceview for use in the
   Editor tab of the gretl preferences dialog
*/

GtkWidget *create_sample_source (const char *style)
{
    GtkSourceLanguageManager *lm;
    GtkSourceLanguage *lang;
    GtkSourceBuffer *sbuf;
    GtkTextView *view;
    GtkWidget *text;

    const gchar *sample =
	"# sample of highlighting style\n"
	"open wages.gdt\n"
	"series l_wage = log(wage)\n"
	"ols l_wage 0 male school exper --robust\n";

    lm = gtk_source_language_manager_get_default();
    if (lm == NULL) {
	return NULL;
    }

    ensure_sourceview_path(lm);

    lang = gtk_source_language_manager_get_language(lm, "gretl");
    if (lang == NULL) {
	return NULL;
    }

    sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new_with_language(lang));
    text = gtk_source_view_new_with_buffer(sbuf);
    view = GTK_TEXT_VIEW(text);

    gtk_text_view_set_left_margin(view, 4);
    gtk_text_view_set_right_margin(view, 4);
    gtk_widget_modify_font(text, fixed_font);
    gtk_text_view_set_editable(view, FALSE);
    gtk_text_view_set_cursor_visible(view, FALSE);

    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), sample, -1);
    gtk_source_buffer_set_highlight_syntax(sbuf, TRUE);
    set_style_for_buffer(sbuf, style, 0);

    return text;
}

/* callback after changes made in preferences dialog */

void update_script_editor_options (windata_t *vwin)
{
    ensure_sourceview_path(NULL);

    if (vwin->role != CONSOLE) {
	gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(vwin->text),
					      script_line_numbers);
    }
    set_style_for_buffer(vwin->sbuf, get_sourceview_style(), vwin->role);

#ifdef HAVE_GTKSV_COMPLETION
    if (vwin->role == CONSOLE || editing_hansl(vwin->role)) {
	set_sv_completion(vwin);
    }
#endif
}

/* Modify the tag-table members that specify a monospaced font in
   response to the user's changing the font size for a given window.
   Otherwise elements tagged in this way will not change size with
   the rest of the text.
*/

static void revise_mono_tag (GtkTextTag *tag, gpointer data)
{
    PangoFontDescription *pfd = data;

    if (widget_get_int(tag, "mono")) {
	g_object_set(tag, "font-desc", pfd, NULL);
    }
}

static void revise_mono_tags (GtkTextTagTable *TT,
			      PangoFontDescription *pfd)
{
    gtk_text_tag_table_foreach(TT, revise_mono_tag, pfd);
}

static int text_change_size (windata_t *vwin, int delta)
{
    static PangoFontDescription *hpf;
    static int fsize0;
    GtkTextTagTable *TT;
    GtkTextBuffer *tbuf;
    int fsize;

    if (hpf == NULL) {
	hpf = pango_font_description_copy(fixed_font);
	fsize0 = pango_font_description_get_size(hpf) / PANGO_SCALE;
    }

    if (vwin == NULL) {
	return fsize0;
    }

    fsize = widget_get_int(vwin->text, "fsize");
    if (fsize == 0) {
	fsize = fsize0;
    }

    fsize += delta;

    pango_font_description_set_size(hpf, fsize * PANGO_SCALE);
    gtk_widget_modify_font(vwin->text, hpf);
    widget_set_int(vwin->text, "fsize", fsize);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    TT = gtk_text_buffer_get_tag_table(tbuf);
    if (TT != NULL) {
	revise_mono_tags(TT, hpf);
    }

    return fsize;
}

void text_larger (GtkWidget *w, gpointer data)
{
    text_change_size((windata_t *) data, 1);
}

void text_smaller (GtkWidget *w, gpointer data)
{
    text_change_size((windata_t *) data, -1);
}

#ifdef OS_OSX
# define helpfont "Geneva"
#else
# define helpfont "sans"
#endif

static GtkTextTagTable *gretl_tags_new (int role)
{
    const char *code_bg;
    GtkTextTagTable *table;
    GtkTextTag *tag;
    int bigsize;
    int smallsize;
    int tags_level = 1;

    if (help_role(role)) {
	tags_level = 3;
    } else if (medium_markup(role)) {
	tags_level = 2;
    }

#if 0
    fprintf(stderr, "gretl_tags_new: role %d, tags_level %d (cf. VIEW_FILE=%d)\n",
	    role, tags_level, VIEW_FILE);
#endif

    code_bg = dark_theme_active() ? "#216fb6" : "#e6f3ff";

    bigsize = 15 * PANGO_SCALE;
    smallsize = 8 * PANGO_SCALE;

    table = gtk_text_tag_table_new();

    tag = gtk_text_tag_new("bluetext");
    g_object_set(tag, "foreground", blue_for_text(), NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("redtext");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("greentext");
    g_object_set(tag, "foreground", "green", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("bold");
    g_object_set(tag, "family", helpfont,
		 "weight", PANGO_WEIGHT_BOLD,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("italic");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_ITALIC,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("replaceable");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_ITALIC,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("text");
    g_object_set(tag, "family", helpfont, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("mono");
    g_object_set(tag, "font-desc", fixed_font, NULL);
    gtk_text_tag_table_add(table, tag);
    widget_set_int(tag, "mono", 1);

    if (tags_level == 1) {
	/* we shouldn't need the rest of these tags */
	return table;
    }

    tag = gtk_text_tag_new("literal");
    g_object_set(tag, "font-desc", fixed_font, NULL);
    gtk_text_tag_table_add(table, tag);
    widget_set_int(tag, "mono", 1);

    tag = gtk_text_tag_new("indented");
    g_object_set(tag, "left_margin", 16, "indent", -12, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("code");
    g_object_set(tag, "font-desc", fixed_font,
		 "paragraph-background", code_bg, NULL);
    gtk_text_tag_table_add(table, tag);
    widget_set_int(tag, "mono", 1);

    if (tags_level == 2) {
	/* the tags below are used only in the "online" help */
	return table;
    }

    tag = gtk_text_tag_new("title");
    g_object_set(tag, "justification", GTK_JUSTIFY_CENTER,
		 "pixels_above_lines", 15,
		 "family", helpfont,
		 "size", bigsize, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("subhead");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_ITALIC,
		 "pixels_below_lines", 8,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("grayhead");
    g_object_set(tag, "family", helpfont,
		 "foreground", "gray",
		 "pixels_below_lines", 5,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("superscript");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_NORMAL,
		 "rise", 4 * PANGO_SCALE,
		 "size", smallsize,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("subscript");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_ITALIC,
		 "rise", -3 * PANGO_SCALE,
		 "size", smallsize,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("subscript-numeral");
    g_object_set(tag, "family", helpfont,
		 "style", PANGO_STYLE_NORMAL,
		 "rise", -3 * PANGO_SCALE,
		 "size", smallsize,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("optflag");
    g_object_set(tag, "font-desc", fixed_font,
		 "foreground", "#396d60", NULL);
    gtk_text_tag_table_add(table, tag);
    widget_set_int(tag, "mono", 1);

    return table;
}

/* be conservative for now */
#define no_tags(r) ((r >= VIEW_SERIES &&	\
		     r <= VIEW_PKG_SAMPLE &&	\
		     r != VIEW_PKG_INFO) ||	\
		    r < NC)

static GtkTextBuffer *gretl_text_buf_new (int role)
{
    GtkTextTagTable *TT = NULL;
    GtkTextBuffer *tbuf;

    if (!no_tags(role)) {
	TT = gretl_tags_new(role);
    }

    tbuf = gtk_text_buffer_new(TT);
    if (TT != NULL) {
	g_object_unref(TT);
    }

    return tbuf;
}

/* net out the effect of one or more instances of '\r' in @buf */

gchar *unctrlr (const char *buf)
{
    gchar *ret = g_malloc0(strlen(buf) + 1);
    int i, j = 0, k = 0;

    for (i=0; buf[i]; i++) {
        if (buf[i] == '\r') {
            j = k;
        } else if (buf[i] == '\n') {
            strcat(ret, "\n");
            k = j = strlen(ret);
        } else {
            ret[j++] = buf[i];
        }
    }

    return ret;
}

static void
real_textview_add_colorized (GtkWidget *view, const char *buf,
			     int append, int trim)
{
    GtkTextBuffer *tb;
    GtkTextIter iter;
    int nextcolor, thiscolor = PLAIN_TEXT;
    int in_comment = 0;
    int blanks = 0;
    char readbuf[4096];
    int console;
    int i = 0;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tb = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));

    if (append) {
	gtk_text_buffer_get_end_iter(tb, &iter);
    } else {
	gtk_text_buffer_get_iter_at_offset(tb, &iter, 0);
    }

    if (strchr(buf, '\r')) {
        gchar *ubuf = unctrlr(buf);

        gtk_text_buffer_insert(tb, &iter, ubuf, -1);
        g_free(ubuf);
        return;
    }

    console = widget_get_int(view, "console");
    bufgets_init(buf);

    while (bufgets(readbuf, sizeof readbuf, buf)) {
	if (trim && i++ < 2) {
	    continue;
	}

	/* try to avoid successive blank lines */
	if (*readbuf == '\n' && readbuf[1] == '\0') {
	    if (++blanks > 1) {
		continue;
	    }
	} else {
	    blanks = 0;
	}

	if (ends_with_backslash(readbuf)) {
	    nextcolor = thiscolor;
	} else {
	    nextcolor = PLAIN_TEXT;
	}

	if (*readbuf == '#' || *readbuf == '?' ||
	    *readbuf == '>' || in_comment) {
	    thiscolor = BLUE_TEXT;
	} else if (!strncmp(readbuf, "/*", 2)) {
	    in_comment = 1;
	    thiscolor = nextcolor = BLUE_TEXT;
	}

	if (strstr(readbuf, "*/")) {
	    in_comment = 0;
	    nextcolor = PLAIN_TEXT;
	}

	if (thiscolor == BLUE_TEXT) {
	    gtk_text_buffer_insert_with_tags_by_name(tb, &iter,
						     readbuf, -1,
						     "bluetext", NULL);
	} else if (console) {
	    gtk_text_buffer_insert_with_tags_by_name(tb, &iter,
						     readbuf, -1,
						     "output", NULL);
	} else {
	    gtk_text_buffer_insert(tb, &iter, readbuf, -1);
	}

	thiscolor = nextcolor;
    }

    bufgets_finalize(buf);
}

void textview_set_text_colorized (GtkWidget *view, const char *buf)
{
    real_textview_add_colorized(view, buf, 0, 0);
}

void textview_append_text_colorized (GtkWidget *view, const char *buf, int trim)
{
    real_textview_add_colorized(view, buf, 1, trim);
}

void textview_set_text_report (GtkWidget *view, const char *buf)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    const char *p;

    /* plain text except for "<@ok>" represented as "OK" in
       green and "<@fail>" as "failed" in red
    */

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    while ((p = strstr(buf, "<@"))) {
	gtk_text_buffer_insert(tbuf, &iter, buf, p - buf);
	if (!strncmp(p, "<@ok>", 5)) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, _("OK"),
						     -1, "greentext", NULL);
	    buf = p + 5;
	} else if (!strncmp(p, "<@fail>", 7)) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter, _("failed"),
						     -1, "redtext", NULL);
	    buf = p + 7;
	}
    }

    gtk_text_buffer_insert(tbuf, &iter, buf, -1);
}

static const char *semicolon_pos (const char *s)
{
    while (*s != '\0' && *s != '"') {
	if (*s == ';') {
	    return s;
	}
	s++;
    }
    return NULL;
}

void textview_set_text_dbsearch (windata_t *vwin, const char *buf)
{
    GtkTextBuffer *tbuf;
    GtkTextTagTable *tab;
    GtkTextIter iter;
    GtkTextTag *tag;
    gchar *dsname = NULL;
    gchar *show = NULL;
    int page;
    const char *p, *q;

    /* plain text except for "<@dbn>" tags for links */

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    tab = gtk_text_buffer_get_tag_table(tbuf);
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    while ((p = strstr(buf, "<@dbn"))) {
	gtk_text_buffer_insert(tbuf, &iter, buf, p - buf);
	p += 7;
	q = semicolon_pos(p);
	if (q != NULL) {
	    /* should be show;tagname */
	    page = DBS_PAGE;
	    show = g_strndup(p, q-p);
	    p = q + 1;
	    q = strchr(p, '"');
	    dsname = g_strndup(p, q-p);
	} else {
	    q = strchr(p, '"');
	    dsname = g_strndup(p, q-p);
	    if (!strcmp(dsname, "_NEXT_")) {
		page = NEXT_PAGE;
		show = g_strdup(_("Next results"));
	    } else {
		page = DBN_PAGE;
		show = NULL;
	    }
	}
	tag = gtk_text_tag_table_lookup(tab, dsname);
	if (tag == NULL) {
	    tag = gtk_text_buffer_create_tag(tbuf, dsname, "foreground",
                                             blue_for_text(),
					     NULL);
	    object_set_int(tag, "page", page);
	}
	gtk_text_buffer_insert_with_tags(tbuf, &iter,
					 show == NULL ? dsname : show,
					 -1, tag, NULL);
	g_free(dsname);
	if (show != NULL) {
	    g_free(show);
	}
	buf = q + 2;
    }

    gtk_text_buffer_insert(tbuf, &iter, buf, -1);

    connect_link_signals(vwin);
}

void textview_delete_processing_message (GtkWidget *view)
{
    GtkTextBuffer *tbuf;
    GtkTextMark *m0, *m1;
    GtkTextIter i0, i1;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    m0 = gtk_text_buffer_get_mark(tbuf, "pstart");
    m1 = gtk_text_buffer_get_mark(tbuf, "pstop");
    if (m0 != NULL && m1 != NULL) {
	gtk_text_buffer_get_iter_at_mark(tbuf, &i0, m0);
	gtk_text_buffer_get_iter_at_mark(tbuf, &i1, m1);
	gtk_text_buffer_delete(tbuf, &i0, &i1);
	gtk_text_buffer_delete_mark(tbuf, m0);
	gtk_text_buffer_delete_mark(tbuf, m1);
    }
}

void textview_add_processing_message (GtkWidget *view)
{
    const char *msg = N_("processing...\n");
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    GtkTextMark *mstop;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_end_iter(tbuf, &iter);
    gtk_text_buffer_create_mark(tbuf, "pstart", &iter, TRUE);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     _(msg), -1,
					     "redtext", NULL);
    mstop = gtk_text_buffer_create_mark(tbuf, "pstop", &iter, TRUE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(view), mstop, 0.0,
				 FALSE, 0, 0);
}

void textview_append_text (GtkWidget *view, const char *text)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_end_iter(tbuf, &iter);
    gtk_text_buffer_insert(tbuf, &iter, text, -1);
}

void textview_insert_text (GtkWidget *view, const char *text)
{
    GtkTextBuffer *tbuf;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_insert_at_cursor(tbuf, text, -1);
}

void textview_clear_text (GtkWidget *view)
{
    GtkTextBuffer *tbuf;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_set_text (tbuf, "", -1);
}

void textview_insert_file (windata_t *vwin, const char *fname)
{
    FILE *fp;
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    int thiscolor, nextcolor;
    char fline[MAXSTR], *chunk;
    int links = 0;
    int i = 0;

    g_return_if_fail(GTK_IS_TEXT_VIEW(vwin->text));

    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	file_read_errbox(fname);
	return;
    }

    thiscolor = nextcolor = PLAIN_TEXT;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    memset(fline, 0, sizeof fline);

    while (fgets(fline, sizeof fline, fp)) {
	if (!g_utf8_validate(fline, -1, NULL)) {
	    if (i == 0) {
		chunk = my_locale_to_utf8(fline);
		i++;
	    } else {
		chunk = my_locale_to_utf8_next(fline);
	    }
	    if (chunk == NULL) {
		continue;
	    }
	} else {
	    chunk = fline;
	}

	nextcolor = PLAIN_TEXT;

	if (vwin->role == VIEW_DOC && strchr(chunk, '<')) {
	    if (!links) {
		links = insert_text_with_markup(tbuf, &iter, chunk, vwin->role);
	    } else {
		insert_text_with_markup(tbuf, &iter, chunk, vwin->role);
	    }
	} else {
	    if (vwin->role == SCRIPT_OUT && ends_with_backslash(chunk)) {
		nextcolor = BLUE_TEXT;
	    }

	    if (*chunk == '?') {
		thiscolor = (vwin->role == CONSOLE)? RED_TEXT : BLUE_TEXT;
	    } else if (*chunk == '#') {
		thiscolor = BLUE_TEXT;
	    }

	    switch (thiscolor) {
	    case PLAIN_TEXT:
		gtk_text_buffer_insert(tbuf, &iter, chunk, -1);
		break;
	    case BLUE_TEXT:
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
							 chunk, -1,
							 "bluetext", NULL);
		break;
	    case RED_TEXT:
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
							 chunk, -1,
							 "redtext", NULL);
		break;
	    }
	}

	if (chunk != fline) {
	    g_free(chunk);
	}

	thiscolor = nextcolor;
	memset(fline, 0, sizeof fline);
    }

    fclose(fp);

    if (links) {
	connect_link_signals(vwin);
    }
}

static gchar *get_mnu_string (const char *key)
{
    const char *s;
    gchar *ret;

    if (!strcmp(key, "LocalGfn")) {
	s = _("On _local machine...");
    } else if (!strcmp(key, "RemoteGfn")) {
	s = _("On _server...");
    } else if (!strcmp(key, "Pkgbook")) {
	s = _("_Function package guide");
    } else if (!strcmp(key, "Addons")) {
	s = _("Check for _addons");
    } else if (!strcmp(key, "Registry")) {
	s = _("Package registry");
    } else if (!strcmp(key, "GeoplotDoc")) {
	s = _("Creating maps");
    } else if (!strcmp(key, "gretlMPI")) {
	s = _("gretl + MPI");
    } else if (!strcmp(key, "gretlSVM")) {
	s = _("gretl + SVM");
    } else if (!strcmp(key, "SetSeed")) {
	s = _("_Seed for random numbers");
    } else if (!strcmp(key, "gretlDBN")) {
	s = _("dbnomics for gretl");
    } else if (!strcmp(key, "gretlLpsolve")) {
	s = _("gretl + lpsolve");
    } else {
	s = key;
    }

    ret = g_strdup(s);

    gretl_delchar('_', ret);
    gretl_delchar('.', ret);

    return ret;
}

#define TAGLEN 128

static gboolean insert_link (GtkTextBuffer *tbuf, GtkTextIter *iter,
			     const char *text, gint page,
			     const char *indent)
{
    GtkTextTagTable *tab = gtk_text_buffer_get_tag_table(tbuf);
    GtkTextTag *tag;
    gchar *show = NULL;
    gchar tagname[TAGLEN];

    if (page == GUIDE_PAGE) {
	char *p = strrchr(text, '#');

	if (p != NULL) {
	    show = g_strndup(text, p - text);
	    strcpy(tagname, p + 1);
	} else {
	    strcpy(tagname, "tag:guide");
	}
    } else if (page == MNU_PAGE) {
	show = get_mnu_string(text);
	strcpy(tagname, text);
    } else if (page == SCRIPT_PAGE || page == EXT_PAGE) {
	*tagname = '\0';
	strncat(tagname, text, TAGLEN-1);
    } else if (page == PDF_PAGE) {
	const char *p = path_last_slash_const(text);

	*tagname = '\0';
	strncat(tagname, text, TAGLEN-1);
	if (p != NULL) {
	    show = g_strdup(p + 1);
	} else {
	    show = g_strdup(text); /* OK? */
	}
    } else if (page == BIB_PAGE) {
	char *p = strrchr(text, ';');

	if (p != NULL) {
	    strcpy(tagname, p + 1);
	    show = g_strndup(text, p - text);
	} else {
	    strcpy(tagname, text);
	}
    } else {
	sprintf(tagname, "tag:p%d", page);
    }

    tag = gtk_text_tag_table_lookup(tab, tagname);

    if (tag == NULL) {
        const char *myblue = blue_for_text();

	if (page == GUIDE_PAGE || page == BIB_PAGE || page == MNU_PAGE) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue,
					     "family", helpfont, NULL);
	} else if (page == SCRIPT_PAGE || page == EXT_PAGE ||
		   page == PDF_PAGE) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue,
					     "font-desc", fixed_font, NULL);
	    widget_set_int(tag, "mono", 1);
	} else if (indent != NULL) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue,
					     "left_margin", 30, NULL);
	} else {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue, NULL);
	}
	object_set_int(tag, "page", page);
    }

    if (show != NULL) {
	gtk_text_buffer_insert_with_tags(tbuf, iter, show, -1, tag, NULL);
	g_free(show);
    } else {
	gtk_text_buffer_insert_with_tags(tbuf, iter, text, -1, tag, NULL);
    }

    return TRUE;
}

static gboolean insert_xlink (GtkTextBuffer *tbuf, GtkTextIter *iter,
			      const char *text, gint page,
			      const char *indent)
{
    GtkTextTagTable *tab = gtk_text_buffer_get_tag_table(tbuf);
    GtkTextTag *tag;
    int gfr = 0;
    gchar tagname[32];

    if (page == GFR_PAGE) {
	strcpy(tagname, "tag:gfr");
	gfr = 1;
	page = 0;
    } else {
	sprintf(tagname, "xtag:p%d", page);
    }

    tag = gtk_text_tag_table_lookup(tab, tagname);

    if (tag == NULL) {
	/* the required tag is not already in the table */
        const char *myblue = blue_for_text();

	if (gfr) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue,
					     "family", helpfont, NULL);
	} else if (indent != NULL) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue,
					     "left_margin", 30, NULL);
	} else {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", myblue, NULL);
	}
	object_set_int(tag, "page", page);
	object_set_int(tag, "xref", 1);
    }

    gtk_text_buffer_insert_with_tags(tbuf, iter, text, -1, tag, NULL);

    return TRUE;
}

static void open_script_link (GtkTextTag *tag)
{
    char fullname[MAXLEN] = {0};
    gchar *fname = NULL;
    int err;

    g_object_get(G_OBJECT(tag), "name", &fname, NULL);
    err = get_full_filename(fname, fullname, OPT_S);
    if (err) {
	errbox_printf(_("Couldn't find %s"), fname);
    } else {
	err = gretl_test_fopen(fullname, "r");
	if (err) {
	    errbox_printf(_("Couldn't read %s"), fullname);
	}
    }
    g_free(fname);

    if (!err) {
	view_script(fullname, 0, VIEW_SCRIPT);
    }
}

static void make_bibitem_window (const char *buf,
				 GtkWidget *tview)
{
    windata_t *vwin;
    GtkWidget *top, *vmain;

    vwin = view_formatted_text_buffer(NULL, buf, 64, 100, VIEW_BIBITEM);
    vmain = vwin_toplevel(vwin);
    top = gtk_widget_get_toplevel(tview);
    gtk_window_set_transient_for(GTK_WINDOW(vmain), GTK_WINDOW(top));
    gtk_window_set_destroy_with_parent(GTK_WINDOW(vmain), TRUE);
    gtk_window_set_position(GTK_WINDOW(vmain),
			    GTK_WIN_POS_CENTER_ON_PARENT);
    gtk_widget_show(vmain);
}

static void open_bibitem_link (GtkTextTag *tag, GtkWidget *tview)
{
    const char *gretldir = gretl_home();
    gchar *key = NULL;
    char fullname[MAXLEN];
    FILE *fp;

    g_object_get(G_OBJECT(tag), "name", &key, NULL);
    sprintf(fullname, "%sgretlhelp.refs", gretldir);
    fp = gretl_fopen(fullname, "r");

    if (fp != NULL) {
	char *buf, line[4096];
	gchar *p, *modbuf;
	int n = strlen(key);

	while (fgets(line, sizeof line, fp)) {
	    if (!strncmp(line, "<@key=\"", 7)) {
		if (!strncmp(line + 7, key, n)) {
		    p = strchr(line + 7, '>');
		    if (p != NULL) {
			buf = p + 1;
			if ((p = strstr(buf, "<@url")) != NULL) {
			    /* put bibitem URL on new line */
			    n = p - buf;
			    modbuf = g_strdup_printf("%.*s\n%s", n, buf, p);
			    make_bibitem_window(modbuf, tview);
			    g_free(modbuf);
			} else {
			    make_bibitem_window(buf, tview);
			}
		    }
		    break;
		}
	    }
	}

	fclose(fp);
    }

    g_free(key);
}

static void open_external_link (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	if (strncmp(name, "http://", 7) &&
	    strncmp(name, "https://", 8)) {
	    gchar *url = g_strdup_printf("http://%s", name);

	    browser_open(url);
	    g_free(url);
	} else {
	    browser_open(name);
	}
	g_free(name);
    }
}

#ifdef GRETL_EDIT

static void open_menu_item (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	/* should be a PDF help file */
	static GtkAction *action;

	if (action == NULL) {
	    action = gtk_action_new(name, NULL, NULL, NULL);
	}
	display_pdf_help(action);
	g_free(name);
    }
}

#else

static void open_menu_item (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	if (!strcmp(name, "RemoteGfn")) {
	    display_files(REMOTE_FUNC_FILES, NULL);
	} else if (!strcmp(name, "LocalGfn")) {
	    display_files(FUNC_FILES, NULL);
	} else if (!strcmp(name, "Addons")) {
	    display_files(ADDONS_FILES, NULL);
	} else if (!strcmp(name, "Registry")) {
	    display_files(PKG_REGISTRY, NULL);
	} else if (!strcmp(name, "SetSeed")) {
	    rand_seed_dialog();
	} else {
	    /* should be a PDF help file */
	    static GtkAction *action;

	    if (action == NULL) {
		action = gtk_action_new(name, NULL, NULL, NULL);
	    }
	    display_pdf_help(action);
	}
	g_free(name);
    }
}

/* opening a series-listing window, coming from a dbnomics
   dataset window */

static void open_dbn_link (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	display_files(DBNOMICS_SERIES, name);
	g_free(name);
    }
}

/* opening a specific series info window, coming from a
   dbnomics dataset-search window */

static void open_dbs_link (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	dbnomics_get_series_call(name);
	g_free(name);
    }
}

/* opening next "page" of dbnomics search results */

static void open_next_link (GtkTextTag *tag, GtkWidget *tview)
{
    windata_t *vwin;

    vwin = g_object_get_data(G_OBJECT(tview), "vwin");
    if (vwin != NULL) {
	dbnomics_search(NULL, vwin);
    }
}

#endif /* GRETL_EDIT or not */

static void open_pdf_file (GtkTextTag *tag)
{
    gchar *name = NULL;

    g_object_get(G_OBJECT(tag), "name", &name, NULL);

    if (name != NULL) {
	int warn = 0;

	if (strchr(name, '/') == NULL && strchr(name, '\\') == NULL) {
	    char *path = get_addon_pdf_path(name);

	    if (path != NULL) {
		gretl_show_pdf(path, NULL);
		free(path);
	    } else {
		/* not an addon file */
		char fname[FILENAME_MAX] = {0};

		warn = get_pdf_path(name, fname);
		if (!warn) {
		    gretl_show_pdf(fname, NULL);
		}
	    }
	} else if (gretl_stat(name, NULL) == 0) {
	    gretl_show_pdf(name, NULL);
	} else {
	    warn = 1;
	}
	if (warn) {
	    warnbox_printf(_("Couldn't open %s"), name);
	}
	g_free(name);
    }
}

static void follow_if_link (GtkWidget *tview, GtkTextIter *iter,
			    windata_t *hwin)
{
    GSList *tags = NULL, *tagp = NULL;
    int en = english_help_role(hwin->role);

    tags = gtk_text_iter_get_tags(iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	GtkTextTag *tag = tagp->data;
	gint page = object_get_int(tag, "page");
	gint xref = object_get_int(tag, "xref");

	if (page != 0 || xref != 0) {
	    if (page == GUIDE_PAGE) {
		gchar *name = NULL;

		g_object_get(tag, "name", &name, NULL);
		if (name != NULL && strstr(name, "chap:")) {
		    display_guide_chapter(name);
		} else {
		    display_pdf_help(NULL);
		}
		g_free(name);
	    } else if (page == SCRIPT_PAGE) {
		open_script_link(tag);
	    } else if (page == BIB_PAGE) {
		open_bibitem_link(tag, tview);
	    } else if (page == EXT_PAGE) {
		open_external_link(tag);
	    } else if (page == PDF_PAGE) {
		open_pdf_file(tag);
	    } else if (page == MNU_PAGE) {
		open_menu_item(tag);
	    }
#ifndef GRETL_EDIT
	    else if (page == DBN_PAGE) {
		open_dbn_link(tag);
	    } else if (page == DBS_PAGE) {
		open_dbs_link(tag);
	    } else if (page == NEXT_PAGE) {
		open_next_link(tag, tview);
	    }
#endif
	    else {
		int role = object_get_int(tview, "role");

		if (function_help(role)) {
		    if (xref) {
			command_help_callback(page, en);
		    } else {
			function_help_callback(page, en);
		    }
		} else {
		    /* commands help */
		    if (xref) {
			function_help_callback(page, en);
		    } else {
			command_help_callback(page, en);
		    }
		}
	    }
	    break;
	}
    }

    if (tags) {
	g_slist_free(tags);
    }
}

/* Help links can be activated by pressing Enter */

static gboolean cmdref_key_press (GtkWidget *tview, GdkEventKey *ev,
				  windata_t *vwin)
{
    GtkTextIter iter;
    GtkTextBuffer *tbuf;

    switch (ev->keyval) {
    case GDK_Return:
    case GDK_KP_Enter:
	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));
	gtk_text_buffer_get_iter_at_mark(tbuf, &iter,
					 gtk_text_buffer_get_insert(tbuf));
	follow_if_link(tview, &iter, vwin);
	break;
    default:
	break;
    }

    return FALSE;
}

static GtkTextTag *active_tag;

static void show_hidden_region (GtkTextBuffer *buffer,
				GtkTextTag *tag,
				GtkTextIter *iter)
{
    GtkTextMark *mark1 = g_object_get_data(G_OBJECT(tag), "mark1");
    GtkTextMark *mark2 = g_object_get_data(G_OBJECT(tag), "mark2");
    GtkTextTagTable *tt;
    GtkTextIter start = *iter;
    GtkTextIter end = *iter;
    GSList *tags, *list = NULL;
    int got = 0, inv = 0;
    int n_hidden;

    g_return_if_fail(mark1 != NULL && mark2 != NULL);

    /* get the bounds of the "show" placeholder */
    got += gtk_text_iter_begins_tag(&start, tag) ||
	gtk_text_iter_backward_to_tag_toggle(&start, tag);
    got += gtk_text_iter_ends_tag(&end, tag) ||
	gtk_text_iter_forward_to_tag_toggle(&end, tag);
    g_return_if_fail(got == 2);

    /* remove the placeholder and its single-use tag */
    tt = gtk_text_buffer_get_tag_table(buffer);
    gtk_text_buffer_remove_tag(buffer, tag, &start, &end);
    gtk_text_buffer_delete(buffer, &start, &end);
    gtk_text_tag_table_remove(tt, tag);

    /* get the bounds of the region to restore */
    gtk_text_buffer_get_iter_at_mark(buffer, &start, mark1);
    gtk_text_buffer_get_iter_at_mark(buffer, &end, mark2);
    /* get the list of applied tags */
    tags = list = gtk_text_iter_get_tags(&start);
    /* find the "invisible" tag and remove it */
    while (tags != NULL) {
	tag = tags->data;
	g_object_get(tag, "invisible", &inv, NULL);
	if (inv) {
	    gtk_text_buffer_remove_tag(buffer, tag, &start, &end);
	    break;
	}
	tags = tags->next;
    }
    g_slist_free(list);
    gtk_text_buffer_delete_mark(buffer, mark1);
    gtk_text_buffer_delete_mark(buffer, mark2);
    n_hidden = object_get_int(buffer, "n_hidden");
    object_set_int(buffer, "n_hidden", --n_hidden);
}

/* Clicking: help links can be activated; hidden regions can be shown */

static gboolean textview_event_after (GtkWidget *w, GdkEvent *ev,
				      windata_t *vwin)
{
    GtkTextIter start, end, iter;
    GtkTextView *view;
    GtkTextBuffer *buffer;
    GdkEventButton *event;
    gint x, y;

    if (ev->type != GDK_BUTTON_RELEASE) {
	return FALSE;
    }

    event = (GdkEventButton *) ev;

    if (event->button != 1) {
	return FALSE;
    }

    view = GTK_TEXT_VIEW(w);
    buffer = gtk_text_view_get_buffer(view);

    /* don't follow a link if the user has selected something */
    gtk_text_buffer_get_selection_bounds(buffer, &start, &end);
    if (gtk_text_iter_get_offset(&start) != gtk_text_iter_get_offset(&end)) {
	return FALSE;
    }

    gtk_text_view_window_to_buffer_coords(view, GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);
    gtk_text_view_get_iter_at_location(view, &iter, x, y);

    if (active_tag != NULL && object_get_int(active_tag, "show")) {
	show_hidden_region(buffer, active_tag, &iter);
    } else {
	follow_if_link(w, &iter, vwin);
    }

    return FALSE;
}

static GdkCursor *hand_cursor = NULL;
static GdkCursor *regular_cursor = NULL;

static void ensure_text_cursors (void)
{
    if (hand_cursor == NULL) {
	hand_cursor = gdk_cursor_new(GDK_HAND2);
    }
    if (regular_cursor == NULL) {
	regular_cursor = gdk_cursor_new(GDK_XTERM);
    }
}

static void
set_cursor_if_appropriate (GtkTextView *view, gint x, gint y,
			   windata_t *vwin)
{
    static gboolean hovering_over_link = FALSE;
    GSList *tags = NULL, *tagp = NULL;
    GtkTextTag *tag = NULL;
    GtkTextIter iter;
    gboolean hovering = FALSE;

    gtk_text_view_get_iter_at_location(view, &iter, x, y);
    tags = gtk_text_iter_get_tags(&iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	tag = tagp->data;
	if (object_get_int(tag, "page") ||
	    object_get_int(tag, "xref") ||
	    object_get_int(tag, "show")) {
	    hovering = TRUE;
	    break;
	}
    }

    if (hovering != hovering_over_link) {
	hovering_over_link = hovering;
	if (hovering_over_link) {
	    active_tag = tag;
	    gdk_window_set_cursor(gtk_text_view_get_window(view, GTK_TEXT_WINDOW_TEXT),
				  hand_cursor);
	} else {
	    active_tag = NULL;
	    gdk_window_set_cursor(gtk_text_view_get_window(view, GTK_TEXT_WINDOW_TEXT),
				  regular_cursor);
	}
    }

    if (tags) {
	g_slist_free(tags);
    }
}

static gboolean
textview_motion_notify (GtkWidget *w, GdkEventMotion *event,
			windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(w);
    gint x, y;

    gtk_text_view_window_to_buffer_coords(view, GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);
    set_cursor_if_appropriate(view, x, y, vwin);

    return FALSE;
}

static gboolean
textview_visibility_notify (GtkWidget *w,  GdkEventVisibility *e,
			    windata_t *vwin)
{
    GtkTextView *view = GTK_TEXT_VIEW(w);
    gint wx, wy, bx, by;

    widget_get_pointer_info(w, &wx, &wy, NULL);
    gtk_text_view_window_to_buffer_coords(view, GTK_TEXT_WINDOW_WIDGET,
					  wx, wy, &bx, &by);
    set_cursor_if_appropriate(view, bx, by, vwin);

    return FALSE;
}

static void connect_link_signals (windata_t *vwin)
{
    ensure_text_cursors();
    if (!vwin_is_editing(vwin)) {
	g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
			 G_CALLBACK(cmdref_key_press), vwin);
    }
    g_signal_connect(G_OBJECT(vwin->text), "event-after",
		     G_CALLBACK(textview_event_after), vwin);
    g_signal_connect(G_OBJECT(vwin->text), "motion-notify-event",
		     G_CALLBACK(textview_motion_notify), vwin);
    g_signal_connect(G_OBJECT(vwin->text), "visibility-notify-event",
		     G_CALLBACK(textview_visibility_notify), vwin);
}

static void maybe_connect_help_signals (windata_t *hwin)
{
    int done = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->text),
						 "sigs_connected"));

    ensure_text_cursors();

    if (!done) {
	g_signal_connect(G_OBJECT(hwin->text), "key-press-event",
			 G_CALLBACK(cmdref_key_press), hwin);
	g_signal_connect(G_OBJECT(hwin->text), "event-after",
			 G_CALLBACK(textview_event_after), hwin);
	g_signal_connect(G_OBJECT(hwin->text), "motion-notify-event",
			 G_CALLBACK(textview_motion_notify), NULL);
	g_signal_connect(G_OBJECT(hwin->text), "visibility-notify-event",
			 G_CALLBACK(textview_visibility_notify), NULL);
        object_set_int(G_OBJECT(hwin->text), "sigs_connected", 1);
    }
}

static void maybe_set_help_tabs (windata_t *hwin)
{
    int done = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->text),
						 "tabs_set"));

    if (!done) {
	PangoTabArray *tabs;

	tabs = pango_tab_array_new(1, TRUE);
	pango_tab_array_set_tab(tabs, 0, PANGO_TAB_LEFT, 50);
	gtk_text_view_set_tabs(GTK_TEXT_VIEW(hwin->text), tabs);
	pango_tab_array_free(tabs);
	object_set_int(hwin->text, "tabs_set", 1);
    }
}

/* Construct the index page for the gretl command reference.
   Note: we assume here that the maximum length of a gretl
   command word is 8 characters.
*/

static void cmdref_index_page (windata_t *hwin, GtkTextBuffer *tbuf)
{
    const char *header = N_("Gretl Command Reference");
    const gchar *s = (const gchar *) hwin->data;
    GtkTextIter iter;
    char word[10];
    int llen, llen_max = 6;
    int idx, j, n, en;

    en = english_help_role(hwin->role);
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     (en)? header : _(header), -1,
					     "title", NULL);
    gtk_text_buffer_insert(tbuf, &iter, "\n\n", -1);

    llen = 0;

    while (*s) {
	if (*s == '\n' && *(s+1) == '#' && *(s+2) != '\0') {
	    if (sscanf(s + 2, "%8s", word)) {
		idx = gretl_help_index(word);
		insert_link(tbuf, &iter, word, idx, NULL);
		if (++llen == llen_max) {
		    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
		    llen = 0;
		} else {
		    n = 10 - strlen(word);
		    for (j=0; j<n; j++) {
			gtk_text_buffer_insert(tbuf, &iter, " ", -1);
		    }
		}
	    }
	}
	s++;
    }

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->text), tbuf);

    maybe_connect_help_signals(hwin);
    maybe_set_help_tabs(hwin);
}

/* construct the index page for the gretl function reference */

static void funcref_index_page (windata_t *hwin, GtkTextBuffer *tbuf)
{
    const char *header = N_("Gretl Function Reference");
    const char *heads[] = {
	N_("Accessors"),
	N_("Built-in strings"),
	N_("Functions proper")
    };
    const gchar *s = (const gchar *) hwin->data;
    gchar *hstr;
    GtkTextIter iter;
    char word[12];
    int llen, llen_max = 5;
    int i, j, k, n, en;

    en = english_help_role(hwin->role);
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     (en)? header : _(header), -1,
					     "title", NULL);
    gtk_text_buffer_insert(tbuf, &iter, "\n\n", -1);

    i = 1;
    k = 0;
    llen = 0;

    while (*s) {
	if (*s == '\n' && s[1] == '#' && s[2] != '\0') {
	    if (s[2] == '#') {
		/* insert section heading */
		if (i > 1) {
		    gtk_text_buffer_insert(tbuf, &iter, "\n\n", -1);
		}
		hstr = g_strdup_printf("%s\n", en ? heads[k] : _(heads[k]));
		gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
							 hstr, -1,
							 "grayhead", NULL);
		g_free(hstr);
		llen = 0;
		s += 2;
		k++;
	    } else if (sscanf(s + 2, "%10s", word)) {
		/* got a function name */
		insert_link(tbuf, &iter, word, i, NULL);
		if (++llen == llen_max) {
		    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
		    llen = 0;
		} else {
		    n = 12 - strlen(word);
		    for (j=0; j<n; j++) {
			gtk_text_buffer_insert(tbuf, &iter, " ", -1);
		    }
		}
		i++;
	    }
	}
	s++;
    }

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->text), tbuf);

    maybe_connect_help_signals(hwin);
    maybe_set_help_tabs(hwin);
}

/* apparatus to support the 'Back' popup menu item */

static void push_backpage (GtkWidget *w, int pg)
{
    object_set_int(w, "backpage", pg);
}

static int pop_backpage (GtkWidget *w)
{
    return object_get_int(w, "backpage");
}

static gint help_popup_click (GtkWidget *w, gpointer p)
{
    windata_t *hwin = (windata_t *) p;
    int action = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    int en = (hwin->role == CMD_HELP_EN || hwin->role == FUNC_HELP_EN);
    int page = 0;

    if (action == 2) {
	page = pop_backpage(hwin->text);
    }

    if (function_help(hwin->role)) {
	function_help_callback(page, en);
    } else {
	command_help_callback(page, en);
    }

    return FALSE;
}

static GtkWidget *build_help_popup (windata_t *hwin)
{
    const char *items[] = {
	N_("Index"),
	N_("Back")
    };
    GtkWidget *pmenu = gtk_menu_new();
    GtkWidget *item;
    int i, imin = 0, imax = 2;

    if (hwin->active_var == 0) {
	/* don't offer "Index" if we're in the index */
	imin = 1;
    }

    if (pop_backpage(hwin->text) == 0) {
	/* don't offer "Back" if we haven't been anywhere */
	imax = 1;
    }

    for (i=imin; i<imax; i++) {
	item = gtk_menu_item_new_with_label(_(items[i]));
	object_set_int(item, "action", i+1);
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(help_popup_click),
			 hwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

    return pmenu;
}

gboolean
help_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p)
{
    if (right_click(event)) {
	windata_t *hwin = (windata_t *) p;

	if (hwin->active_var == 0 && pop_backpage(w) == 0) {
	    return TRUE;
	}

	if (hwin->popup) {
	    gtk_widget_destroy(hwin->popup);
	    hwin->popup = NULL;
	}

	hwin->popup = build_help_popup(hwin);

	if (hwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(hwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    g_signal_connect(G_OBJECT(hwin->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed),
			     &hwin->popup);
	}

	return TRUE;
    }

    return FALSE;
}

static void reformat_para (char *buf, int maxlen)
{
    char *p = buf;
    char *line;
    int i, n;

    g_strchomp(g_strchug(buf));

    /* normalize spaces, removing any existing line breaks */

    while (*p) {
	if (*p == ' ' || *p == '\n' || *p == '\t') {
	    *p = ' ';
	    n = strspn(p + 1, " \t\n");
	    if (n > 0) {
		g_strchug(p + 1);
		p += n;
	    }
	}
	p++;
    }

    line = p = buf;
    n = 0;

    /* insert line breaks to give lines of length up
       to @maxlen */

    while (*p++) {
	n++;
	if (n > maxlen) {
	    /* back up to first available break-point */
	    for (i=n-1; i>0; i--) {
		if (line[i] == ' ') {
		    line[i] = '\n';
		    p = line = &line[i+1];
		    n = 0;
		}
		if (n == 0) {
		    break;
		}
	    }
	}
    }
}

static gboolean prev_double_nl (GtkTextIter *pos,
				GtkTextIter *start)
{
    GtkTextIter cpos = *pos;
    int nlcount = 0;
    gunichar c;
    gboolean ret = 0;

    if (gtk_text_iter_get_char(pos) == '\n') {
	nlcount = 1;
    }

    while (gtk_text_iter_backward_char(&cpos)) {
	c = gtk_text_iter_get_char(&cpos);
	if (c == '\n') {
	    if (++nlcount == 2) {
		*start = cpos;
		gtk_text_iter_forward_chars(start, 2);
		ret = 1;
		break;
	    }
	} else if (!isspace(c)) {
	    nlcount = 0;
	}
    }

    return ret;
}

static gboolean next_double_nl (GtkTextIter *pos,
				GtkTextIter *end)
{
    GtkTextIter cpos = *pos;
    int nlcount = 0;
    gunichar c;
    gboolean ret = 0;

    if (gtk_text_iter_get_char(pos) == '\n') {
	nlcount = 1;
    }

    while (gtk_text_iter_forward_char(&cpos)) {
	c = gtk_text_iter_get_char(&cpos);
	if (c == '\n') {
	    if (++nlcount == 2) {
		*end = cpos;
		gtk_text_iter_backward_char(end);
		ret = 1;
		break;
	    }
	} else if (!isspace(c)) {
	    nlcount = 0;
	}
    }

    return ret;
}

static int not_in_para (GtkTextBuffer *buf,
			GtkTextIter *pos)
{
    GtkTextIter cpos = *pos;
    int got_text = 0;
    gunichar c;

    /* We're "not in a paragraph" if the current
       cursor position is bracketed by newlines,
       with no non-space character intervening.
    */

    c = gtk_text_iter_get_char(&cpos);

    if (!isspace(c)) {
	got_text = 1;
    } else if (c == '\n') {
	/* crawl backwards to newline or text */
	while (gtk_text_iter_backward_char(&cpos)) {
	    c = gtk_text_iter_get_char(&cpos);
	    if (c == '\n') {
		break;
	    } else if (!isspace(c)) {
		got_text = 1;
		break;
	    }
	}
    } else if (isspace(c)) {
	/* crawl forward to newline or text */
	while (gtk_text_iter_forward_char(&cpos)) {
	    c = gtk_text_iter_get_char(&cpos);
	    if (c == '\n') {
		break;
	    } else if (!isspace(c)) {
		got_text = 1;
		break;
	    }
	}
	if (!got_text) {
	    /* OK, try backwards */
	    cpos = *pos;
	    while (gtk_text_iter_backward_char(&cpos)) {
		c = gtk_text_iter_get_char(&cpos);
		if (c == '\n') {
		    break;
		} else if (!isspace(c)) {
		    got_text = 1;
		    break;
		}
	    }
	}
    }

    return !got_text;
}

static gboolean textbuf_get_para_limits (GtkTextBuffer *buf,
					 GtkTextIter *pos,
					 GtkTextIter *start,
					 GtkTextIter *end)
{
    if (not_in_para(buf, pos)) {
	return FALSE;
    }

    if (!prev_double_nl(pos, start)) {
	gtk_text_buffer_get_start_iter(buf, start);
    }

    if (!next_double_nl(pos, end)) {
	gtk_text_buffer_get_end_iter(buf, end);
    }

    return TRUE;
}

void textview_format_paragraph (GtkWidget *view)
{
    GtkTextBuffer *buf;
    GtkTextIter pos, start, end;
    gchar *para = NULL;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));

    /* find where the cursor is */
    gtk_text_buffer_get_iter_at_mark(buf, &pos,
				     gtk_text_buffer_get_insert(buf));

    /* find start and end of paragraph, if we're in one */
    if (!textbuf_get_para_limits(buf, &pos, &start, &end)) {
	return;
    }

    /* grab the para text */
    para = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

    if (para != NULL && !string_is_blank(para)) {
	reformat_para(para, 72);
	gtk_text_buffer_begin_user_action(buf);
	gtk_text_buffer_delete(buf, &start, &end);
	gtk_text_buffer_insert(buf, &start, para, -1);
	gtk_text_buffer_end_user_action(buf);
	g_free(para);

    }
}

static gchar *textview_get_current_line (GtkWidget *view, int allow_blank)
{
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    gchar *ret = NULL;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_iter_at_mark(buf, &start,
				     gtk_text_buffer_get_insert(buf));
    gtk_text_iter_set_line_offset(&start, 0);
    gtk_text_buffer_get_iter_at_mark(buf, &end,
				     gtk_text_buffer_get_insert(buf));
    if (!gtk_text_iter_ends_line(&end)) {
	/* N.B. don't skip on to the end of the _next_ line */
	gtk_text_iter_forward_to_line_end(&end);
    }

    ret = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

    if (!allow_blank && string_is_blank(ret)) {
	g_free(ret);
	ret = NULL;
    }

    return ret;
}

#ifndef GRETL_EDIT

static gchar *textview_get_current_line_with_newline (GtkWidget *view)
{
    gchar *s = textview_get_current_line(view, 0);

    if (s != NULL && *s != '\0' && s[strlen(s)-1] != '\n') {
	gchar *tmp = g_strdup_printf("%s\n", s);

	g_free(s);
	s = tmp;
    }

    return s;
}

#endif

/* Determine whether or not any of the lines in a chunk of text
   are indented, via spaces or tabs.
*/

static int text_is_indented (const gchar *s)
{
    int leading = 1;

    if (s == NULL) {
	return 0;
    }

    while (*s) {
	if (*s == '\n') {
	    leading = 1;
	} else if (*s != ' ' && *s != '\t') {
	    leading = 0;
	}
	if (leading && (*s == ' ' || *s == '\t')) {
	    return 1;
	}
	s++;
    }

    return 0;
}

/* Determine whether or not a chunk of text is commented, in the form
   of each line beginning with '#' (with possible leading white
   space).  If some lines are commented and others are not, return -1,
   which blocks the comment/uncomment menu items.
*/

static int text_is_commented (const gchar *s)
{
    int gotc = 0, comm = 0;
    int lines = 1;

    if (s == NULL) {
	return -1;
    }

    while (*s) {
	if (!gotc) {
	    if (*s == '#') {
		comm++;
		gotc = 1;
	    } else if (!isspace(*s)) {
		gotc = 1;
	    }
	} else if (*s == '\n') {
	    gotc = 0;
	    if (*(s+1)) {
		lines++;
	    }
	}
	s++;
    }

    if (comm > 0 && comm < lines) {
	/* mixed */
	comm = -1;
    }

    return comm;
}

struct textbit {
    windata_t *vwin;
    GtkTextBuffer *buf;
    GtkTextIter start;
    GtkTextIter end;
    gchar *chunk;
    int commented;
    int selected;
};

/* either insert or remove '#' comment markers at the start of the
   line(s) of a chunk of text
*/

static void comment_or_uncomment_text (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;
    gchar *s;

    gtk_text_buffer_delete(tb->buf, &tb->start, &tb->end);

    if (tb->selected) {
	char line[1024];

	bufgets_init(tb->chunk);
	while (bufgets(line, sizeof line, tb->chunk)) {
	    if (tb->commented) {
		s = strchr(line, '#');
		if (s != NULL) {
		    s++;
		    if (*s == ' ') s++;
		    gtk_text_buffer_insert(tb->buf, &tb->start, s, -1);
		}
	    } else {
		gtk_text_buffer_insert(tb->buf, &tb->start, "# ", -1);
		gtk_text_buffer_insert(tb->buf, &tb->start, line, -1);
	    }
	}
	bufgets_finalize(tb->chunk);
    } else {
	if (tb->commented) {
	    s = strchr(tb->chunk, '#');
	    if (s != NULL) {
		s++;
		if (*s == ' ') s++;
		gtk_text_buffer_insert(tb->buf, &tb->start, s, -1);
	    }
	} else {
	    gtk_text_buffer_insert(tb->buf, &tb->start, "# ", -1);
	    gtk_text_buffer_insert(tb->buf, &tb->start, tb->chunk, -1);
	}
    }
}

enum {
    TAB_NEXT,
    TAB_PREV
};

static int spaces_to_tab_stop (const char *s, int targ)
{
    int ret, n = 0;

    while (*s) {
	if (*s == ' ') {
	    n++;
	} else if (*s == '\t') {
	    n += tabwidth;
	} else {
	    break;
	}
	s++;
    }

    if (targ == TAB_NEXT) {
	ret = tabwidth - (n % tabwidth);
    } else {
	if (n % tabwidth == 0) {
	    ret = n - tabwidth;
	    if (ret < 0) ret = 0;
	} else {
	    ret = (n / tabwidth) * tabwidth;
	}
    }

    return ret;
}

/* This function is called for both indenting an entire script and
   indenting a selected region. The @start and @end iters will differ
   in the two cases, as will the size of &buf.

   The selected-region case is marked by the last argument. Note that
   in the whole-script case we automatically append a newline if the
   script is not so terminated, but we shouldn't do that in the region
   case unless the region runs right to the end of the buffer.
*/

static void normalize_indent (GtkTextBuffer *tbuf,
			      const gchar *buf,
			      GtkTextIter *start,
			      GtkTextIter *end,
			      int region)
{
    int strip_nl = 0;
    char *ins;
    PRN *prn;

    if (buf == NULL) {
	return;
    }

    if (region) {
	/* identify the case where a region ended prior to
	   the end of the script, and did not include a
	   trailing newline
	*/
	GtkTextIter final;

	gtk_text_buffer_get_end_iter(tbuf, &final);
	if (!gtk_text_iter_equal(end, &final) &&
	    buf[strlen(buf)-1] != '\n') {
	    strip_nl = 1;
	}
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    normalize_hansl(buf, tabwidth, prn);
    ins = gretl_print_steal_buffer(prn);
    if (strip_nl) {
	/* strip an appended newline if need be */
	int len = strlen(ins);

	if (ins[len-1] == '\n') {
	    ins[len-1] = '\0';
	}
    }
    gtk_text_buffer_delete(tbuf, start, end);
    gtk_text_buffer_insert(tbuf, start, ins, -1);
    free(ins);
    gretl_print_destroy(prn);
}

static int in_foreign_land (GtkWidget *text_widget)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    gchar *buf;
    char *s, line[1024];
    int inforeign = 0;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text_widget));
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_iter_at_mark(tbuf, &end,
				     gtk_text_buffer_get_insert(tbuf));
    buf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	s = line + strspn(line, " \t");
	if (!strncmp(s, "foreign ", 8)) {
	    inforeign = 1;
	} else if (!strncmp(s, "end foreign", 11)) {
	    inforeign = 0;
	}
    }

    bufgets_finalize(buf);
    g_free(buf);

    return inforeign;
}

static void auto_indent_script (GtkWidget *w, windata_t *vwin)
{
    GtkAdjustment *adj;
    GtkTextBuffer *tbuf;
    GtkTextMark *mark;
    GtkTextIter here, start, end;
    gchar *buf;
    gint line;
    gboolean ends;
    gdouble pos;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    /* record scrolling position */
    adj = gtk_text_view_get_vadjustment(GTK_TEXT_VIEW(vwin->text));
    pos = gtk_adjustment_get_value(adj);

    /* record cursor position (line, offset) */
    mark = gtk_text_buffer_get_insert(tbuf);
    gtk_text_buffer_get_iter_at_mark(tbuf, &here, mark);
    line = gtk_text_iter_get_line(&here);
    ends = gtk_text_iter_ends_line(&here);

    /* grab and revise the text */
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);
    buf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    normalize_indent(tbuf, buf, &start, &end, 0);
    g_free(buf);

    /* restore cursor position */
    gtk_text_buffer_get_iter_at_line(tbuf, &here, line);
    if (ends && !gtk_text_iter_ends_line(&here)) {
	gtk_text_iter_forward_to_line_end(&here);
    }
    gtk_text_buffer_place_cursor(tbuf, &here);

    /* restore scrolling position */
    gtk_adjustment_set_value(adj, pos);
    gtk_adjustment_value_changed(adj);
}

static void indent_region (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;

    if (smarttab) {
	normalize_indent(tb->buf, tb->chunk, &tb->start, &tb->end, 1);
    } else {
	char line[1024];
	int i, n;

	gtk_text_buffer_delete(tb->buf, &tb->start, &tb->end);

	bufgets_init(tb->chunk);

	while (bufgets(line, sizeof line, tb->chunk)) {
	    n = spaces_to_tab_stop(line, TAB_NEXT);
	    for (i=0; i<n; i++) {
		gtk_text_buffer_insert(tb->buf, &tb->start, " ", -1);
	    }
	    gtk_text_buffer_insert(tb->buf, &tb->start, line, -1);
	}

	bufgets_finalize(tb->chunk);
    }
}

static void hide_region (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;
    GtkTextMark *mark1, *mark2;
    GtkTextTag *htag = NULL;
    GtkTextTag *gtag = NULL;
    GtkTextTagTable *tt;

    tt = gtk_text_buffer_get_tag_table(tb->buf);
    if (tt != NULL) {
	htag = gtk_text_tag_table_lookup(tt, "hidden");
	/* make single-use placeholder tag */
	gtag = gtk_text_tag_new(NULL);
	g_object_set(gtag, "foreground", "gray", "editable", 0, NULL);
	gtk_text_tag_table_add(tt, gtag);
    }

    if (htag != NULL && gtag != NULL) {
        int n_hidden = object_get_int(tb->buf, "n_hidden");

	gtk_text_buffer_apply_tag(tb->buf, htag, &tb->start, &tb->end);
	mark1 = gtk_text_buffer_create_mark(tb->buf, NULL, &tb->start, FALSE);
	mark2 = gtk_text_buffer_create_mark(tb->buf, NULL, &tb->end, FALSE);
	g_object_set_data(G_OBJECT(gtag), "mark1", mark1);
	g_object_set_data(G_OBJECT(gtag), "mark2", mark2);
        object_set_int(gtag, "show", 1);
	gtk_text_buffer_insert_with_tags(tb->buf, &tb->start, hidden_marker,
					 -1, gtag, NULL);
        object_set_int(tb->buf, "n_hidden", ++n_hidden);
	gtk_text_buffer_place_cursor(tb->buf, &tb->start);
    }
}

void indent_hansl (GtkWidget *w, windata_t *vwin)
{
    auto_indent_script(w, vwin);
}

static void unindent_region (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;
    char line[1024];
    char *ins;
    int i, n;

    gtk_text_buffer_delete(tb->buf, &tb->start, &tb->end);

    bufgets_init(tb->chunk);

    while (bufgets(line, sizeof line, tb->chunk)) {
	n = spaces_to_tab_stop(line, TAB_PREV);
	ins = line + strspn(line, " \t");
	for (i=0; i<n; i++) {
	    gtk_text_buffer_insert(tb->buf, &tb->start, " ", -1);
	}
	gtk_text_buffer_insert(tb->buf, &tb->start, ins, -1);
    }

    bufgets_finalize(tb->chunk);
}

#ifndef GRETL_EDIT

static void exec_script_text (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;

    if (editing_hansl(tb->vwin->role)) {
	run_native_script(tb->vwin, tb->chunk, NULL, 0);
    } else {
	run_script_fragment(tb->vwin, tb->chunk);
    }
    g_free(tb->chunk);
    tb->chunk = NULL;
}

#endif

enum {
    AUTO_SELECT_NONE,
    AUTO_SELECT_LINE
};

static struct textbit *vwin_get_textbit (windata_t *vwin, int mode)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    int selected = 0;
    struct textbit *tb;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	selected = 1;
    }

    if (!selected && mode != AUTO_SELECT_LINE) {
	return NULL;
    }

    tb = malloc(sizeof *tb);
    if (tb == NULL) {
	return NULL;
    }

    tb->vwin = vwin;
    tb->buf = tbuf;
    tb->start = start;
    tb->end = end;
    tb->selected = selected;
    tb->commented = 0;
    tb->chunk = NULL;

    if (selected) {
	int endpos;

	gtk_text_iter_set_line_offset(&tb->start, 0);
	endpos = gtk_text_iter_get_line_offset(&tb->end);
	if (endpos > 0 && !gtk_text_iter_ends_line(&tb->end)) {
	    gtk_text_iter_forward_to_line_end(&tb->end);
	}
	tb->chunk = gtk_text_buffer_get_text(tb->buf, &tb->start, &tb->end, FALSE);
    } else {
	gtk_text_buffer_get_iter_at_mark(tb->buf, &tb->start,
					 gtk_text_buffer_get_insert(tb->buf));
	gtk_text_iter_set_line_offset(&tb->start, 0);
	gtk_text_buffer_get_iter_at_mark(tb->buf, &tb->end,
					 gtk_text_buffer_get_insert(tb->buf));
	gtk_text_iter_forward_to_line_end(&tb->end);
	tb->chunk = gtk_text_buffer_get_text(tb->buf, &tb->start, &tb->end, FALSE);
    }

    return tb;
}

static int count_leading_spaces (const char *s)
{
    int n = 0;

    while (*s) {
	if (*s == ' ') {
	    n++;
	} else if (*s == '\t') {
	    n += tabwidth;
	} else {
	    break;
	}
	s++;
    }

    return n;
}

static int lparen_offset (const char *s)
{
    const char *p = strchr(s, '(');

    if (p != NULL && strchr(p, ')') == NULL) {
	return p - s;
    } else {
	return 0;
    }
}

/* Given what's presumed to be a start-of-line iter, find how many
   leading spaces are on the line, counting tabs as multiple spaces.
*/

static int leading_spaces_at_iter (GtkTextBuffer *tbuf,
				   GtkTextIter *start,
				   const char *word)
{
    GtkTextIter end = *start;
    gchar *s;
    int n = 0;

    gtk_text_iter_forward_to_line_end(&end);
    s = gtk_text_buffer_get_text(tbuf, start, &end, FALSE);
    if (s != NULL) {
	if (!strcmp(word, "function")) {
	    n = lparen_offset(s);
	    if (n > 0) {
		n--;
	    } else {
		n = count_leading_spaces(s);
	    }
	} else {
	    n = count_leading_spaces(s);
	}
	g_free(s);
    }

    return n;
}

static int incremental_leading_spaces (const char *prevword,
				       const char *thisword)
{
    int this_indent = 0;
    int next_indent = 0;

    if (*prevword != '\0') {
	int prev_indent = 0;

	adjust_indent(prevword, &this_indent, &next_indent);
#if TABDEBUG > 1
	fprintf(stderr, "adjust_indent 1: prevword='%s', this=%d, next=%d\n",
		prevword, this_indent, next_indent);
#endif
	prev_indent = this_indent;
	if (*thisword != '\0') {
	    adjust_indent(thisword, &this_indent, &next_indent);
#if TABDEBUG > 1
	    fprintf(stderr, "adjust_indent 2: thisword='%s', this=%d\n",
		    thisword, this_indent);
#endif
	    this_indent -= prev_indent;
#if TABDEBUG > 1
	    fprintf(stderr, "adjust_indent 3: this=%d\n", this_indent);
#endif
	} else {
	    this_indent = next_indent - this_indent;
	}
    }

#if TABDEBUG > 1
    fprintf(stderr, "incremental_leading_spaces: returning %d*%d = %d\n",
	    this_indent, tabwidth, this_indent * tabwidth);
#endif

    return this_indent * tabwidth;
}

static int line_continues (const gchar *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>=0; i--) {
	if (s[i] != ' ' && s[i] != '\t') {
	    return (s[i] == '\\' || s[i] == ',');
	}
    }

    return 0;
}

static int get_word_and_cont (const char *s, char *word, int *contd)
{
    /* don't move onto next line */
    if (*s != '\n' && *s != '\r') {
	s += strspn(s, " \t");
	if (sscanf(s, "%*s <- %8s", word) != 1) {
	    sscanf(s, "%8s", word);
	}
	if (*word != '#' && contd != NULL) {
	    *contd = line_continues(s);
	}
    }

    return *word != '\0';
}

static int line_continues_previous (GtkTextBuffer *tbuf,
				    GtkTextIter iter)
{
    GtkTextIter end, prev = iter;
    gchar *s;
    int ret = 0;

    if (gtk_text_iter_backward_line(&prev)) {
	end = prev;
	if (gtk_text_iter_forward_to_line_end(&end)) {
	    s = gtk_text_buffer_get_text(tbuf, &prev, &end, FALSE);
	    if (s != NULL) {
		if (*s != '\n' && *s != '\r') {
		    ret = line_continues(s);
		}
		g_free(s);
	    }
	}
    }

    return ret;
}

/* get "command word", max 8 characters: work backwards up script
   to find this */

static char *get_previous_line_start_word (char *word,
					   GtkTextBuffer *tbuf,
					   GtkTextIter iter,
					   int *leadspace,
					   int *contd)
{
    GtkTextIter end, prev = iter;
    int *pcont = contd;
    int rparen = 0;
    int i = 0;
    gchar *s;

    *word = '\0';

    while (*word == '\0' && gtk_text_iter_backward_line(&prev)) {
	end = prev;
	if (gtk_text_iter_forward_to_line_end(&end)) {
	    s = gtk_text_buffer_get_text(tbuf, &prev, &end, FALSE);
	    if (s != NULL) {
		if (i == 0 && s[strlen(s)-1] == ')') {
		    rparen = 1;
		}
		if (get_word_and_cont(s, word, pcont)) {
		    pcont = NULL;
		}
		g_free(s);
	    }
	}

	if (line_continues_previous(tbuf, prev)) {
	    /* back up one line further */
	    *word = '\0';
	}

	if (*word != '\0' && leadspace != NULL) {
	    *leadspace = leading_spaces_at_iter(tbuf, &prev, word);
	}
	i++;
    }

    if (rparen && !strcmp(word, "function")) {
	/* Revise our judgement: looks like we must be on the
	   first line following a function signature, so we do
	   not want the deep indent suitable for a continuation
	   signature line.
	*/
	*leadspace = 0;
    }

    return word;
}

/* Is the insertion point at the start of a line, or in a white-space
   field to the left of any non-space characters?  If so, we'll trying
   inserting a "smart" soft tab in response to the Tab key. If not,
   and @comp_ok is non-NULL, we'll check whether gtksourceview
   completion seems possible at the current insertion point.

   Plus a special case (a la emacs): if we're at the end of line which
   should end preceding indentation (such as "endif" or "endloop"),
   let Tab at the end of the line adjust its indentation.
*/

static int maybe_insert_smart_tab (windata_t *vwin, int *comp_ok)
{
    GtkTextBuffer *tbuf;
    GtkTextMark *mark;
    GtkTextIter start, end;
    int curr_nsp = 0;
    int pos, ret = 0;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	/* don't do this if there's a selection in place */
	return 0;
    }

    /* find @pos, the line-index of the insertion point */
    mark = gtk_text_buffer_get_insert(tbuf);
    gtk_text_buffer_get_iter_at_mark(tbuf, &end, mark);
    pos = gtk_text_iter_get_line_offset(&end);

    if (pos == 0) {
	/* we're at the left margin, OK */
	ret = 1;
    } else {
	gchar *chunk;

	start = end;
	gtk_text_iter_set_line_offset(&start, 0);
	/* grab the text between line start and insertion point */
	chunk = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	if (strspn(chunk, " \t") == pos) {
	    /* set @ret if this chunk is just white space */
	    ret = 1;
	} else if (comp_ok != NULL && pos > 1) {
	    /* follow-up: is the context OK for completion? */
	    *comp_ok = !isspace(chunk[pos-1]) && !isspace(chunk[pos-2]);
	}
	g_free(chunk);
    }

    if (ret) {
	/* OK, let's insert a smart tab */
	GtkTextIter prev = start;
	char *s, thisword[9] = {0};
	char prevword[9];
	int contd = 0;
	int i, nsp = 0;

	s = textview_get_current_line(vwin->text, 1);
	if (s != NULL) {
#if TABDEBUG > 1
	    fprintf(stderr, "*** maybe_insert_smart_tab: "
		    "current line = '%s'\n", s);
#endif
	    sscanf(s, "%8s", thisword);
	    curr_nsp = strspn(s, " \t");
	    g_free(s);
	}

	get_previous_line_start_word(prevword, tbuf, prev, &nsp, &contd);

	if (contd) {
	    nsp += 2;
	} else {
	    nsp += incremental_leading_spaces(prevword, thisword);
#if TABDEBUG > 1
	    fprintf(stderr, "    leading spaces: nsp + incr = %d, curr_nsp %d\n",
		    nsp, curr_nsp);
#endif
	}

	if (curr_nsp > 0) {
	    end = start;
	    gtk_text_iter_set_line_offset(&end, curr_nsp);
	    gtk_text_buffer_delete(tbuf, &start, &end);
	}
	gtk_text_iter_set_line_offset(&start, 0);
	for (i=0; i<nsp; i++) {
	    gtk_text_buffer_insert(tbuf, &start, " ", -1);
	}
    }

    return ret;
}

#ifdef HAVE_GTKSV_COMPLETION

/* Is the insertion point directly preceded by at least two
   non-space characters? If so we have a possible candidate for
   completion via gtksourceview. We come here only if "smart
   tab" is not in force; otherwise this check is rolled into
   the check for inserting a smart tab.

   This function is public because it's called from console.c
   as well as internally.
*/

int maybe_try_completion (windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    GtkTextMark *mark;
    GtkTextIter start, end;
    int pos, ret = 0;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	/* don't do this if there's a selection in place? */
	return 0;
    }

    /* find @pos, the line-index of the insertion point */
    mark = gtk_text_buffer_get_insert(tbuf);
    gtk_text_buffer_get_iter_at_mark(tbuf, &end, mark);
    pos = gtk_text_iter_get_line_offset(&end);

    if (pos > 1) {
	const char *test;
	gchar *chunk;

	start = end;
	gtk_text_iter_set_line_offset(&start, 0);
	test = chunk = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	if (vwin->role == CONSOLE && !strncmp(chunk, "? ", 2)) {
	    test += 2;
	    pos -= 2;
	}
	ret = !isspace(test[pos-1]) && !isspace(test[pos-2]);
	g_free(chunk);
    }

    return ret;
}

#endif /* HAVE_GTKSV_COMPLETION */

static char leftchar (guint k)
{
    return k == GDK_parenleft ? '(' :
	k == GDK_bracketleft ? '[' : '{';
}

static char rightchar (guint k)
{
    return k == GDK_parenleft ? ')' :
	k == GDK_bracketleft ? ']' : '}';
}

/* Is the insertion point at the end of a line? If so, we'll
   auto-insert a matching right bracket and move the cursor
   back before it.
*/

static int maybe_insert_auto_bracket (windata_t *vwin,
				      guint keyval)
{
    GtkTextBuffer *tbuf;
    GtkTextMark *mark;
    GtkTextIter start;
    gchar *chunk = NULL;
    int ret = 0;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_has_selection(tbuf)) {
	return 0;
    }

    mark = gtk_text_buffer_get_insert(tbuf);
    gtk_text_buffer_get_iter_at_mark(tbuf, &start, mark);

    if (gtk_text_iter_ends_line(&start)) {
	ret = 1;
    } else {
	GtkTextIter end = start;

	gtk_text_iter_forward_to_line_end(&end);
	chunk = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	if (chunk != NULL) {
	    ret = strspn(chunk, " \t\n") == strlen(chunk);
	    g_free(chunk);
	}
    }

    if (ret) {
	char s[2] = {0};

	s[0] = leftchar(keyval);
	gtk_text_buffer_insert(tbuf, &start, s, -1);
	s[0] = rightchar(keyval);
	gtk_text_buffer_insert(tbuf, &start, s, -1);
	gtk_text_iter_backward_char(&start);
	gtk_text_buffer_place_cursor(tbuf, &start);
    }

    return ret;
}

/* On "Enter" in script editing, try to compute the correct indent
   level for the current line, and make an adjustment if it's not
   already right. We also attempt to place the cursor at the
   appropriate indent on the next, new line -- unless @alt is
   non-zero, in which case we don't insert a newline.
*/

static gboolean script_electric_enter (windata_t *vwin, int alt)
{
    char *s = NULL;
    int targsp = 0;
    gboolean ret = FALSE;

    if (!smarttab || in_foreign_land(vwin->text)) {
	return FALSE;
    }

#if TABDEBUG
    fprintf(stderr, "*** script_electric_enter: alt = %d\n", alt);
#endif

    s = textview_get_current_line(vwin->text, 0);

    if (s != NULL && *s != '\0') {
	GtkTextBuffer *tbuf;
	GtkTextMark *mark;
	GtkTextIter start, end;
	char thisword[9];
	char prevword[9];
	int i, diff, nsp, incr;
	int ipos, k, contd = 0;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_get_iter_at_mark(tbuf, &start, mark);
	ipos = gtk_text_iter_get_line_offset(&start);
	gtk_text_iter_set_line_offset(&start, 0);

	*thisword = '\0';
	sscanf(s, "%8s", thisword);
	nsp = count_leading_spaces(s);
	get_previous_line_start_word(prevword, tbuf, start, &targsp, &contd);

#if TABDEBUG
	if (contd) {
	    fprintf(stderr, "prevword='%s', leading space = %d\n",
		    prevword, targsp);
	    fprintf(stderr, "got line continuation\n");
	} else {
	    fprintf(stderr, "thisword='%s', leading space = %d\n",
		    thisword, nsp);
	    fprintf(stderr, "prevword='%s', leading space = %d\n",
		    prevword, targsp);
	}
#endif

	if (contd) {
	    incr = 2;
	} else {
#if TABDEBUG > 1
	    fprintf(stderr, "getting leading spaces ('%s', '%s')\n",
		    prevword, thisword);
#endif
	    incr = incremental_leading_spaces(prevword, thisword);
	}

	targsp += incr;
	if (targsp < 0) {
	    /* indentation messed up? */
	    targsp = 0;
	}

	diff = nsp - targsp;

#if TABDEBUG
	fprintf(stderr, "incr = %d: after increment targsp = %d, diff = %d\n",
		incr, targsp, diff);
#endif
	if (diff > 0) {
	    end = start;
	    gtk_text_iter_forward_chars(&end, diff);
	    gtk_text_buffer_delete(tbuf, &start, &end);
	} else if (diff < 0) {
	    diff = -diff;
	    for (i=0; i<diff; i++) {
		gtk_text_buffer_insert(tbuf, &start, " ", -1);
	    }
	}

	if (!alt) {
	    /* try to arrange correct indent on the new line? */
	    if (ipos == 0) {
		k = targsp + incremental_leading_spaces(prevword, "");
	    } else {
		k = targsp + incremental_leading_spaces(thisword, "");
	    }
#if TABDEBUG
	    fprintf(stderr, "new line indent: k = %d\n", k);
#endif
	    gtk_text_buffer_begin_user_action(tbuf);
	    gtk_text_buffer_get_iter_at_mark(tbuf, &start, mark);
	    gtk_text_buffer_insert(tbuf, &start, "\n", -1);
	    for (i=0; i<k; i++) {
		gtk_text_buffer_insert(tbuf, &start, " ", -1);
	    }
	    gtk_text_buffer_place_cursor(tbuf, &start);
	    gtk_text_buffer_end_user_action(tbuf);
	    mark = gtk_text_buffer_get_insert(tbuf);
	    gtk_text_view_scroll_mark_onscreen(GTK_TEXT_VIEW(vwin->text), mark);
	}
	ret = TRUE;
    }

    g_free(s);

    return ret;
}

/* Handler for the Tab key when editing a gretl script: we
   may want to insert a "smart tab", or take Tab as a
   request for gtksourceview completion.
*/

#ifdef HAVE_GTKSV_COMPLETION

static gboolean script_tab_handler (windata_t *vwin, GdkEvent *event)
{
    int ucomp;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(vwin->text), FALSE);

    if (in_foreign_land(vwin->text)) {
	return FALSE;
    } else if (((GdkEventKey *) event)->state & GDK_SHIFT_MASK) {
	return FALSE;
    }

    ucomp = (hansl_completion == COMPLETE_USER);

    if (smarttab) {
	int comp_ok = 0;
	int *ptr = ucomp ? &comp_ok : NULL;

	if (maybe_insert_smart_tab(vwin, ptr)) {
	    return TRUE;
	} else if (comp_ok) {
	    call_user_completion(vwin->text);
	    return TRUE;
	}
    } else if (ucomp && maybe_try_completion(vwin)) {
	call_user_completion(vwin->text);
	return TRUE;
    }

    return FALSE;
}

#else

static gboolean script_tab_handler (windata_t *vwin, GdkEvent *event)
{
    g_return_val_if_fail(GTK_IS_TEXT_VIEW(vwin->text), FALSE);

    if (in_foreign_land(vwin->text)) {
	return FALSE;
    } else if (((GdkEventKey *) event)->state & GDK_SHIFT_MASK) {
	return FALSE;
    }

    if (smarttab && maybe_insert_smart_tab(vwin, NULL)) {
	return TRUE;
    }

    return FALSE;
}

#endif /* HAVE_GTKSV_COMPLETION or not */

gboolean script_bracket_handler (windata_t *vwin, guint keyval)
{
    if (maybe_insert_auto_bracket(vwin, keyval)) {
	return TRUE;
    } else {
	return FALSE;
    }
}

/* Return a listing of the available gtksourceview style
   ids for use in the gretl preferences dialog.
*/

const char **get_sourceview_style_ids (int *n)
{
    GtkSourceStyleSchemeManager *mgr;
    const gchar * const *ids = NULL;

    ensure_sourceview_path(NULL);

    *n = 0;

    mgr = gtk_source_style_scheme_manager_get_default();
    if (mgr != NULL) {
	int i = 0;

	ids = gtk_source_style_scheme_manager_get_scheme_ids(mgr);
	if (ids != NULL) {
	    while (ids[i] != NULL) i++;
	    *n = i;
	}
    }

    return (const char **) ids;
}

const char **get_graph_theme_ids (int *n)
{
    static char **S = NULL;
    static int n_found;

    if (S != NULL) {
	*n = n_found;
    } else {
	gchar *path;
	GDir *dir;

	*n = 0;
	path = g_build_filename(gretl_home(), "data", "gnuplot", NULL);
	dir = gretl_opendir(path);

	S = strings_array_new(1);
	S[0] = gretl_strdup("classic");
	*n = 1;

	if (dir != NULL) {
	    const gchar *fname;
	    gchar *tmp, *p;
	    int err = 0;

	    while (!err && (fname = g_dir_read_name(dir))) {
		if (!strncmp(fname, "default.", 8) ||
		    !strncmp(fname, "classic.", 8)) {
		    continue;
		}
		if (has_suffix(fname, ".gpsty")) {
		    tmp = g_strdup(fname);
		    p = strstr(tmp, ".gpsty");
		    *p = '\0';
		    err = strings_array_add(&S, n, tmp);
		    g_free(tmp);
		}
	    }
	    g_dir_close(dir);
	}
	n_found = *n;
	g_free(path);
    }

    return (const char **) S;
}

static void call_prefs_dialog (GtkWidget *w, windata_t *vwin)
{
    preferences_dialog(TAB_EDITOR, NULL, vwin_toplevel(vwin));
}

static GtkWidget *
build_script_popup (windata_t *vwin, struct textbit **ptb)
{
    const char *items[] = {
	N_("Comment line"),
	N_("Uncomment line"),
	N_("Comment region"),
	N_("Uncomment region"),
    };
    GtkWidget *pmenu = NULL;
    struct textbit *tb = NULL;
    GtkWidget *item;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(vwin->text), NULL);

    /* "generic" text window menu -- we may add to this */
    pmenu = build_text_popup(vwin);

    if (foreign_script_role(vwin->role)) {
	*ptb = NULL;
	goto dock_undock;
    }

    tb = vwin_get_textbit(vwin, AUTO_SELECT_LINE);
    if (tb == NULL) {
	*ptb = NULL;
	return pmenu;
    }

    tb->commented = text_is_commented(tb->chunk);

    if (tb->commented > 0 && !editing_hansl(vwin->role)) {
	g_free(tb->chunk);
	free(tb);
	*ptb = NULL;
	goto dock_undock;
    }

    *ptb = tb;

#ifndef GRETL_EDIT
    if (tb->commented <= 0 && vwin->role != EDIT_PKG_CODE) {
	/* we have some uncommented material: allow exec option */
	if (tb->selected) {
	    item = gtk_menu_item_new_with_label(_("Execute region"));
	} else {
	    item = gtk_menu_item_new_with_label(_("Execute line"));
	}
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(exec_script_text),
			 *ptb);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }
#endif

    if (editing_hansl(vwin->role) && tb->commented >= 0) {
	/* material is either all commented or all uncommented:
	   allow comment/uncomment option
	*/
	int i = (tb->selected && !tb->commented)? 2 :
	    (tb->selected && tb->commented)? 3 :
	    (!tb->selected && !tb->commented)? 0 : 1;

	item = gtk_menu_item_new_with_label(_(items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(comment_or_uncomment_text),
			 *ptb);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

    if (editing_hansl(vwin->role)) {
	if (tb->selected) {
	    /* indentation of selection */
	    item = gtk_menu_item_new_with_label(smarttab?
						_("Auto-indent region") :
						_("Indent region"));
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(indent_region),
			     *ptb);
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
	    /* visibility of selection */
	    item = gtk_menu_item_new_with_label(_("Hide region"));
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(hide_region),
			     *ptb);
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
	}
	if (!smarttab && tb->selected && text_is_indented(tb->chunk)) {
	    item = gtk_menu_item_new_with_label(_("Unindent region"));
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(unindent_region),
			     *ptb);
	    gtk_widget_show(item);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
	}

	item = gtk_menu_item_new_with_label(_("Auto-indent script (Ctrl+I)"));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(auto_indent_script),
			 vwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

 dock_undock:

    if (GTK_IS_SOURCE_VIEW(vwin->text)) {
	item = gtk_menu_item_new_with_label(_("Preferences..."));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(call_prefs_dialog),
			 vwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

    if (window_is_undockable(vwin)) {
	add_undock_popup_item(pmenu, vwin);
    } else if (window_is_dockable(vwin)) {
	add_dock_popup_item(pmenu, vwin);
    }

    return pmenu;
}

static gboolean destroy_textbit (GtkWidget **pw, struct textbit *tc)
{
    if (tc != NULL) {
	tc->vwin->popup = NULL;
	g_free(tc->chunk);
	free(tc);
    }

    return FALSE;
}

static gboolean
script_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p)
{
    if (right_click(event)) {
	windata_t *vwin = (windata_t *) p;
	struct textbit *tc = NULL;

	if (vwin->popup) {
	    gtk_widget_destroy(vwin->popup);
	    vwin->popup = NULL;
	}

	vwin->popup = build_script_popup(vwin, &tc);

	if (vwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(vwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    g_signal_connect(G_OBJECT(vwin->popup), "destroy",
			     G_CALLBACK(destroy_textbit),
			     tc);
	}

	return TRUE;
    }

    return FALSE;
}

enum {
    INSERT_NONE,
    INSERT_REF,
    INSERT_XREF,
    INSERT_FIG,
    INSERT_REPL,
    INSERT_LIT,
    INSERT_URL,
    INSERT_OPT,
    INSERT_ITAL,
    INSERT_SUP,
    INSERT_SUB,
    INSERT_TEXT,
    INSERT_MATH,
    INSERT_GUGLINK,
    INSERT_INPLINK,
    INSERT_GFRLINK,
    INSERT_BIBLINK,
    INSERT_ADBLINK,
    INSERT_MNULINK,
    INSERT_BOLD
};

static void insert_help_figure (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *fig)
{
    char figfile[FILENAME_MAX];
    GdkPixbuf *pixbuf;

    sprintf(figfile, "%shelpfigs%c%s.png", gretl_home(),
	    SLASH, fig);

    pixbuf = gdk_pixbuf_new_from_file(figfile, NULL);

    if (pixbuf != NULL) {
	gtk_text_buffer_insert_pixbuf(tbuf, iter, pixbuf);
	g_object_unref(pixbuf);
    }
}

static void insert_math_content (GtkTextBuffer *tbuf, GtkTextIter *iter,
				 const char *s, const char *indent)
{
    static char minus[4];
    gchar ubuf[6];
    gunichar c;
    int i, n;

    if (*minus == '\0') {
	/* find the best representation of minus */
	PangoFontDescription *pfd;

	pfd = pango_font_description_from_string(helpfont);

	if (font_has_symbol(pfd, 0x2212)) {
	    /* preferred: unicode minus */
	    minus[0] = 0xE2;
	    minus[1] = 0x88;
	    minus[2] = 0x92;
	} else if (font_has_symbol(pfd, 0x2013)) {
	    /* fallback: unicode endash */
	    minus[0] = 0xE2;
	    minus[1] = 0x80;
	    minus[2] = 0x93;
	} else {
	    /* otherwise: plain old dash */
	    minus[0] = '-';
	}
	pango_font_description_free(pfd);
    }

    n = g_utf8_strlen(s, -1);

    for (i=0; i<n; i++) {
	c = g_utf8_get_char(s);
	if (*s == '-') {
	    gtk_text_buffer_insert(tbuf, iter, minus, -1);
	} else {
	    memset(ubuf, 0, sizeof ubuf);
	    g_unichar_to_utf8(c, ubuf);
	    if (g_unichar_isalpha(c)) {
		gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, ubuf, -1,
							 "italic", indent, NULL);
	    } else {
		gtk_text_buffer_insert(tbuf, iter, ubuf, -1);
	    }
	}
	if (i < n-1) {
	    s = g_utf8_find_next_char(s, NULL);
	}
    }
}

static void insert_tagged_text (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *s, int ins, const char *indent)
{
    const char *ftag = NULL;

    switch (ins) {
    case INSERT_ITAL:
    case INSERT_MATH: /* FIXME */
	ftag = "italic";
	break;
    case INSERT_REPL:
	ftag = "replaceable";
	break;
    case INSERT_LIT:
	ftag = "literal";
	break;
    case INSERT_OPT:
	ftag = "optflag";
	break;
    case INSERT_SUP:
	ftag = "superscript";
	break;
    case INSERT_SUB:
	if (integer_string(s)) {
	    ftag = "subscript-numeral";
	} else {
	    ftag = "subscript";
	}
	break;
    case INSERT_BOLD:
	ftag = "bold";
	break;
    default:
	break;
    }

#ifdef G_OS_WIN32
    if (ins == INSERT_OPT) {
	/* Unicode word joiner not supported? Try zero width
	   non breaking space instead */
	char tmp[32];

	strcpy(tmp, s);
	tmp[2] = 0xEF;
	tmp[3] = 0xBB;
	tmp[4] = 0xBF;

	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, tmp, -1,
						 ftag, indent, NULL);
	return;
    }
#endif

    if (ftag != NULL) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 ftag, indent, NULL);
    }
}

static gchar *get_string_and_instruction (const char *p, int *ins)
{
    gchar *str = NULL;

    *ins = INSERT_NONE;

    if (!strncmp(p, "ref", 3)) {
	*ins = INSERT_REF;
    } else if (!strncmp(p, "xrf", 3)) {
	*ins = INSERT_XREF;
    } else if (!strncmp(p, "fig", 3)) {
	*ins = INSERT_FIG;
    } else if (!strncmp(p, "itl", 3)) {
	*ins = INSERT_ITAL;
    } else if (!strncmp(p, "bld", 3)) {
	*ins = INSERT_BOLD;
    } else if (!strncmp(p, "var", 3)) {
	*ins = INSERT_REPL;
    } else if (!strncmp(p, "lit", 3)) {
	*ins = INSERT_LIT;
    } else if (!strncmp(p, "url", 3)) {
	*ins = INSERT_URL;
    } else if (!strncmp(p, "opt", 3)) {
	*ins = INSERT_OPT;
    } else if (!strncmp(p, "sup", 3)) {
	*ins = INSERT_SUP;
    } else if (!strncmp(p, "sub", 3)) {
	*ins = INSERT_SUB;
    } else if (!strncmp(p, "mth", 3)) {
	*ins = INSERT_MATH;
    } else if (!strncmp(p, "pdf", 3)) {
	*ins = INSERT_GUGLINK;
    } else if (!strncmp(p, "inp", 3)) {
	*ins = INSERT_INPLINK;
    } else if (!strncmp(p, "gfr", 3)) {
	*ins = INSERT_GFRLINK;
    } else if (!strncmp(p, "bib", 3)) {
	*ins = INSERT_BIBLINK;
    } else if (!strncmp(p, "adb", 3)) {
	*ins = INSERT_ADBLINK;
    } else if (!strncmp(p, "mnu", 3)) {
	*ins = INSERT_MNULINK;
    } else if (!strncmp(p, "hd1", 3)) {
	*ins = INSERT_BOLD;
    }

    if (*ins != INSERT_NONE) {
	int i;

	p += 5; /* skip 'tag="' */
	for (i=0; p[i]; i++) {
	    if (p[i] == '"' && p[i+1] == '>') {
		break;
	    }
	}
	str = g_strndup(p, i);
    }

    return str;
}

static int get_code_skip (const char *s)
{
    int skip = 5;

    while (*s) {
	if (*s == '\n') {
	    skip++;
	    break;
	} else if (isspace(*s)) {
	    skip++;
	}
	s++;
    }

    return skip;
}

static int command_word_index (char *s)
{
    int i = gretl_help_index(s);

    if (i == 0) {
	i = extra_command_number(s);
	if (i < 0) {
	    i = 0;
	} else {
	    /* e.g. "MIDAS_list" -> "MIDAS list", for
	       display purposes */
	    gretl_charsub(s, '_', ' ');
	}
    }

    return i;
}

/* note: subhead must be a one-liner */

static int handle_subhead (GtkTextBuffer *tbuf,
			   GtkTextIter *iter,
			   const char *p)
{
    const char *q;

    q = strstr(p, "</sub");
    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, p, q-p,
					     "subhead", NULL);
    q += 10; /* skip "</subhead>" */
    while (*q == ' ') {
	/* skip any trailing spaces */
	q++;
    }

    /* add 1 to swallow a newline */
    return 1 + q - p;
}

/* return non-zero if we inserted any hyperlinks */

static gboolean
insert_text_with_markup (GtkTextBuffer *tbuf, GtkTextIter *iter,
			 const char *s, int role)
{
    gboolean ret = FALSE;
    gchar *targ = NULL;
    const char *indent = NULL;
    const char *p;
    int code = 0;
    int mono = 0;
    int itarg, ins;

    while ((p = strstr(s, "<"))) {
	int skip = 0;

	if (code) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "code", indent, NULL);
	} else if (mono) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "mono", indent, NULL);
	} else {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "text", indent, NULL);
	}

	p++;

	if (*p == '@') {
	    /* "atomic" markup */
	    targ = get_string_and_instruction(p + 1, &ins);
	    if (ins == INSERT_REF) {
		if (function_help(role)) {
		    /* FIXME? */
		    itarg = function_help_index_from_word(targ, role);
		} else {
		    itarg = command_word_index(targ);
		}
		ret = insert_link(tbuf, iter, targ, itarg, indent);
	    } else if (ins == INSERT_XREF) {
		if (function_help(role)) {
		    itarg = command_word_index(targ);
		} else {
		    itarg = function_help_index_from_word(targ, FUNC_HELP);
		}
		ret = insert_xlink(tbuf, iter, targ, itarg, indent);
	    } else if (ins == INSERT_URL) {
		ret = insert_link(tbuf, iter, targ, EXT_PAGE, indent);
	    } else if (ins == INSERT_GUGLINK) {
		ret = insert_link(tbuf, iter, targ, GUIDE_PAGE, indent);
	    } else if (ins == INSERT_INPLINK) {
		ret = insert_link(tbuf, iter, targ, SCRIPT_PAGE, indent);
	    } else if (ins == INSERT_BIBLINK) {
		ret = insert_link(tbuf, iter, targ, BIB_PAGE, indent);
	    } else if (ins == INSERT_GFRLINK) {
		ret = insert_xlink(tbuf, iter, targ, GFR_PAGE, indent);
	    } else if (ins == INSERT_ADBLINK) {
		ret = insert_link(tbuf, iter, targ, PDF_PAGE, indent);
	    } else if (ins == INSERT_MNULINK) {
		ret = insert_link(tbuf, iter, targ, MNU_PAGE, indent);
	    } else if (ins == INSERT_FIG) {
		insert_help_figure(tbuf, iter, targ);
	    } else if (ins == INSERT_MATH) {
		insert_math_content(tbuf, iter, targ, indent);
	    } else if (ins != INSERT_NONE) {
		insert_tagged_text(tbuf, iter, targ, ins, indent);
	    }
	    skip = 8 + strlen(targ);
	} else if (!strncmp(p, "indent", 6)) {
	    indent = "indented";
	    skip = 7 + (*(p+7) == '\n');
	} else if (!strncmp(p, "/indent", 7)) {
	    indent = NULL;
	    skip = 8 + (*(p+8) == '\n');
	} else if (!strncmp(p, "code", 4)) {
	    code = 1;
	    skip = get_code_skip(p + 5);
	} else if (!strncmp(p, "/code", 5)) {
	    code = 0;
	    skip = 6 + (*(p+6) == '\n');
	} else if (!strncmp(p, "mono", 4)) {
	    mono = 1;
	    skip = get_code_skip(p + 5);
	} else if (!strncmp(p, "/mono", 5)) {
	    mono = 0;
	    skip = 6 + (*(p+6) == '\n');
	} else if (!strncmp(p, "subhead", 7)) {
	    p += 8;
	    skip = handle_subhead(tbuf, iter, p);
	} else {
	    /* literal "<" */
	    gtk_text_buffer_insert(tbuf, iter, "<", 1);
	}

	if (targ != NULL) {
	    g_free(targ);
	    targ = NULL;
	}
	s = p + skip;
    }

    if (code) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "code", indent, NULL);
    } else if (mono) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "mono", indent, NULL);
    } else {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "text", indent, NULL);
    }

    return ret;
}

static char *grab_topic_buffer (const char *s)
{
    const char *p = strstr(s, "\n#");
    char *buf;

    if (p != NULL) {
	buf = g_strndup(s, p - s);
    } else {
	buf = g_strdup(s);
    }

    return buf;
}

/* Pull the appropriate chunk of help text out of the buffer attached
   to the help viewer and display it, if possible.  Return >= 0
   if we did OK, < 0 on failure.
*/

int set_help_topic_buffer (windata_t *hwin, int pos)
{
    GtkTextBuffer *textb;
    GtkTextIter iter;
    char line[256];
    gchar *hbuf;
    char *buf;

    push_backpage(hwin->text, hwin->active_var);

    textb = gretl_text_buf_new(hwin->role);

    if (pos == 0) {
	/* no topic selected */
	if (function_help(hwin->role)) {
	    funcref_index_page(hwin, textb);
	} else {
	    cmdref_index_page(hwin, textb);
	}
	cursor_to_top(hwin);
	return 0;
    }

    /* OK, so pos is non-zero */
    maybe_connect_help_signals(hwin);
    maybe_set_help_tabs(hwin);

    gtk_text_buffer_get_iter_at_offset(textb, &iter, 0);

    hbuf = (gchar *) hwin->data + pos;

    bufgets_init(hbuf);
    buf = bufgets(line, sizeof line, hbuf);
    bufgets_finalize(hbuf);

    if (buf == NULL) {
	return -1;
    }

    tailstrip(line);

    if (gui_help(hwin->role)) {
	/* topic heading: descriptive string */
	gchar *p = quoted_help_string(line);

	gtk_text_buffer_insert_with_tags_by_name(textb, &iter,
						 p, -1,
						 "bold", NULL);
	free(p);
    } else {
	/* topic heading: plain command word */
	char hword[12];

	sscanf(line + 2, "%11s", hword);
	gtk_text_buffer_insert_with_tags_by_name(textb, &iter,
						 hword, -1,
						 "redtext", NULL);
    }

    if (function_help(hwin->role)) {
	gtk_text_buffer_insert(textb, &iter, "\n\n", 2);
    } else {
	gtk_text_buffer_insert(textb, &iter, "\n", 1);
    }

    buf = grab_topic_buffer(hbuf + strlen(line) + 1);
    if (buf == NULL) {
	return -1;
    }

    insert_text_with_markup(textb, &iter, buf, hwin->role);
    free(buf);

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->text), textb);
    maybe_connect_help_signals(hwin);
    cursor_to_top(hwin);

    return 1;
}

void gretl_viewer_set_formatted_buffer (windata_t *vwin,
					const char *buf,
					int role)
{
    GtkTextBuffer *textb;
    GtkTextIter iter;
    gboolean links;

    textb = gretl_text_buf_new(vwin->role);
    gtk_text_buffer_get_iter_at_offset(textb, &iter, 0);
    links = insert_text_with_markup(textb, &iter, buf, FUNC_HELP);

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(vwin->text), textb);
    cursor_to_top(vwin);

    if (links) {
	connect_link_signals(vwin);
    }
}

void gretl_viewer_insert_formatted_buffer (windata_t *vwin, const char *buf)
{
    GtkTextBuffer *textb;
    GtkTextIter start, end;

    textb = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_bounds(textb, &start, &end);
    gtk_text_buffer_delete(textb, &start, &end);
    insert_text_with_markup(textb, &start, buf, FUNC_HELP);
}

int get_screen_height (void)
{
    static int screen_height;

    if (screen_height == 0) {
	GdkScreen *s = gdk_screen_get_default();

	if (s != NULL) {
	    screen_height = gdk_screen_get_height(s);
	}
    }

    return screen_height;
}

static void set_max_text_width (windata_t *vwin,
				int width,
				int height)
{
    GdkGeometry hints = {0};

    hints.max_width = width;
    hints.max_height = height * 3;

    gtk_window_set_geometry_hints(GTK_WINDOW(vwin->main),
				  GTK_WIDGET(vwin->main),
				  &hints,
				  GDK_HINT_MAX_SIZE);
}

#define HDEBUG 0
#define HELP_WRAP 1

void create_text (windata_t *vwin, int hsize, int vsize,
		  int nlines, gboolean editable)
{
    GtkTextBuffer *tbuf = gretl_text_buf_new(vwin->role);
    GtkWidget *w = gtk_text_view_new_with_buffer(tbuf);
    int role = vwin->role;

    vwin->text = w;

    /* in which cases should we do text wrapping? */

    if (help_role(role) || role == VIEW_PKG_INFO ||
	role == VIEW_BIBITEM || role == VIEW_CODEBOOK ||
#if HELP_WRAP
	role == EDIT_PKG_HELP || role == EDIT_PKG_GHLP ||
#endif
	role == VIEW_DBNOMICS || role == IMPORT ||
	role == VIEW_DBSEARCH || role == INFO) {
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    } else {
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_NONE);
    }

    if (role == VIEW_DBSEARCH) {
	/* make @vwin discoverable via @w */
	g_object_set_data(G_OBJECT(w), "vwin", vwin);
    }

    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(w), 4);

    gtk_widget_modify_font(GTK_WIDGET(w), fixed_font);

#if HDEBUG
    /* the incoming @hsize is expressed in characters */
    fprintf(stderr, "create_text: initial hsize = %d\n", hsize);
#endif

    if (hsize > 0 || nlines > 0) {
        int sv = get_screen_height();
	int px, py;

	get_char_width_and_height(w, &px, &py);
	if (hsize > 0) {
	    hsize *= px;
	    hsize += 3 * px;
	}
#if HDEBUG
	fprintf(stderr, " px = %d, py = %d; hsize now = %d, nlines = %d\n",
		px, py, hsize, nlines);
#endif
	if (nlines > 0) {
	    /* Perhaps adjust how tall the window is? */
	    int v1 = (nlines + 2) * py;

	    if (v1 > 0.85 * vsize && v1 <= 0.7 * sv) {
		vsize = v1;
	    } else if (v1 > 0.7 * sv) {
		vsize = 0.7 * sv;
	    }
        } else if (role == SCRIPT_OUT && vsize < 0.65 * sv) {
            vsize = 0.65 * sv;
	} else if (role != VIEW_BIBITEM && vsize < 0.62 * hsize) {
	    vsize = 0.62 * hsize;
	}
    }

    if (hsize > 0 && vsize > 0) {
	GtkWidget *vmain = vwin_toplevel(vwin);

	if (window_is_tab(vwin)) {
	    vsize += 15;
	}
#if HDEBUG
	fprintf(stderr, " setting default size (%d, %d)\n", hsize, vsize);
#endif
	gtk_window_set_default_size(GTK_WINDOW(vmain), hsize, vsize);
	if (role == EDIT_PKG_HELP || role == EDIT_PKG_GHLP) {
	    /* for editing help files: limit the width of the window
	       to discourage use of excessively long lines
	    */
	    set_max_text_width(vwin, hsize, vsize);
	}
    }

    gtk_text_view_set_editable(GTK_TEXT_VIEW(w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), editable);
}

static GtkTextTagTable *gretl_console_tags_new (void)
{
    GtkTextTagTable *table;
    GtkTextTag *tag;

    table = gtk_text_tag_table_new();

    tag = gtk_text_tag_new("prompt");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("output");
    g_object_set(tag, "foreground", "black",
		 "weight", PANGO_WEIGHT_NORMAL, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("bluetext");
    g_object_set(tag, "foreground", blue_for_text(),
		 "weight", PANGO_WEIGHT_NORMAL, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("redtext");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    return table;
}

void create_console (windata_t *vwin, int hsize, int vsize)
{
    static GtkTextTagTable *console_tags = NULL;
    GtkSourceLanguageManager *lm = NULL;
    GtkSourceBuffer *sbuf;
    GtkTextView *view;
    int cw;

    if (console_tags == NULL) {
	console_tags = gretl_console_tags_new();
    }

    /* since 2018-06-09: use syntax highlighting */
    lm = gtk_source_language_manager_get_default();
    ensure_sourceview_path(lm);

    sbuf = gtk_source_buffer_new(console_tags);
    gtk_source_buffer_set_highlight_matching_brackets(sbuf, TRUE);
    if (lm != NULL) {
	g_object_set_data(G_OBJECT(sbuf), "languages-manager", lm);
	set_style_for_buffer(sbuf, get_sourceview_style(), CONSOLE);
    }

    vwin->text = gtk_source_view_new_with_buffer(sbuf);
    widget_set_int(vwin->text, "console", 1);
    vwin->sbuf = sbuf;

    view = GTK_TEXT_VIEW(vwin->text);
    gtk_text_view_set_wrap_mode(view, GTK_WRAP_NONE);
    gtk_text_view_set_left_margin(view, 4);
    gtk_text_view_set_right_margin(view, 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->text), fixed_font);

    cw = get_char_width(vwin->text);
    set_source_tabs(vwin->text, cw);

#ifdef HAVE_GTKSV_COMPLETION
    if (console_completion) {
	/* since 2021-06-11 */
	set_sv_completion(vwin);
    }
#endif

    if (hsize > 0) {
	hsize *= cw;
	hsize += 48; /* ?? */
    }
    if (vsize < 0.62 * hsize) {
	vsize = 0.62 * hsize;
    }

    if (vwin->main != vwin->vbox && hsize > 0 && vsize > 0) {
	GtkWidget *vmain = vwin_toplevel(vwin);

	gtk_window_set_default_size(GTK_WINDOW(vmain), hsize, vsize);
    }

    sourceview_apply_language(vwin);
    gtk_text_view_set_editable(view, TRUE);
    gtk_text_view_set_cursor_visible(view, TRUE);
}

void text_set_word_wrap (GtkWidget *w, gboolean wrap)
{
    if (wrap) {
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    } else {
	gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_NONE);
    }
}

void text_table_setup (GtkWidget *vbox, GtkWidget *w)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vbox), sw, TRUE, TRUE, FALSE);
    g_object_set_data(G_OBJECT(vbox), "sw", sw);

    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), w);
    gtk_widget_show(w);
    gtk_widget_show(sw);
}

static void set_pane_text_properties (GtkWidget *w2,
				      GtkWidget *w1)
{
    GtkTextView *tv2 = GTK_TEXT_VIEW(w2);
    GtkTextView *tv1 = GTK_TEXT_VIEW(w1);
    gboolean s;

    gtk_text_view_set_wrap_mode(tv2, 0);
    gtk_text_view_set_left_margin(tv2, 4);
    gtk_text_view_set_right_margin(tv2, 4);
    gtk_widget_modify_font(w2, fixed_font);

    s = gtk_text_view_get_editable(tv1);
    gtk_text_view_set_editable(tv2, s);
    gtk_text_view_set_cursor_visible(tv2, s);

    /* sourceview line numbering? */
    if (GTK_IS_SOURCE_VIEW(w2)) {
	s = gtk_source_view_get_show_line_numbers(GTK_SOURCE_VIEW(w1));
	gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(w2), s);
    }
}

/* divide a text window into two panes */

void viewer_split_pane (windata_t *vwin, int vertical)
{
    GtkWidget *vbox = vwin->vbox;
    GtkWidget *view1 = vwin->text;
    GtkWidget *sw, *paned, *view2;
    GtkWidget *vmain;
    GtkTextBuffer *tbuf;
    gint width, height;

    vmain = vwin_toplevel(vwin);
    gtk_window_get_size(GTK_WINDOW(vmain), &width, &height);

    sw = g_object_get_data(G_OBJECT(vbox), "sw");
    g_object_ref(sw);
    gtk_container_remove(GTK_CONTAINER(vwin->vbox), sw);

    if (vertical) {
	paned = gtk_hpaned_new();
    } else {
	paned = gtk_vpaned_new();
    }

    gtk_container_set_border_width(GTK_CONTAINER(paned), 0);
    gtk_container_add(GTK_CONTAINER(vbox), paned);

    g_object_set_data(G_OBJECT(vwin->vbox), "paned", paned);
    g_object_set_data(G_OBJECT(vwin->vbox), "sw", NULL);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view1));

    if (GTK_IS_SOURCE_VIEW(view1)) {
	view2 = gtk_source_view_new_with_buffer(GTK_SOURCE_BUFFER(tbuf));
    } else {
	view2 = gtk_text_view_new_with_buffer(tbuf);
    }

    set_pane_text_properties(view2, view1);

    g_signal_connect(G_OBJECT(view2), "button-press-event",
		     G_CALLBACK(text_popup_handler), vwin);

    gtk_paned_add1(GTK_PANED(paned), sw);
    g_object_unref(sw);

    sw = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_paned_add2(GTK_PANED(paned), sw);
    gtk_container_add(GTK_CONTAINER(sw), view2);

    if (vertical) {
	gtk_paned_set_position(GTK_PANED(paned), width / 2 - 5);
    } else {
	gtk_paned_set_position(GTK_PANED(paned), height / 2 - 24);
    }

    gtk_widget_show_all(paned);
}

/* script output window: revert to a single pane */

void viewer_close_pane (windata_t *vwin)
{
    GtkWidget *sw, *paned;

    paned = g_object_get_data(G_OBJECT(vwin->vbox), "paned");

    /* grab the first child and reference it */
    sw = gtk_paned_get_child1(GTK_PANED(paned));
    g_object_ref(sw);
    gtk_container_remove(GTK_CONTAINER(paned), sw);

    /* remove the "paned" widget */
    gtk_widget_destroy(paned);
    g_object_set_data(G_OBJECT(vwin->vbox), "paned", NULL);

    /* and repack the scrolled window */
    gtk_container_add(GTK_CONTAINER(vwin->vbox), sw);
    g_object_set_data(G_OBJECT(vwin->vbox), "sw", sw);
    gtk_widget_show(sw);
    g_object_unref(sw);
}
