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
#include "toolbar.h"
#include "dlgutils.h"

#include <gtksourceview/gtksourceview.h>
#include <gtksourceview/gtksourcelanguage.h>
#ifdef USE_GTKSOURCEVIEW_2
# include <gtksourceview/gtksourcelanguagemanager.h>
#else
# include <gtksourceview/gtksourcelanguagesmanager.h>
#endif

#define GUIDE_PAGE  999
#define SCRIPT_PAGE 998
#define GFR_PAGE    997

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

#define gui_help(r) (r == GUI_HELP || r == GUI_HELP_EN)
#define foreign_script_role(r) (r == EDIT_GP || r == EDIT_R || r == EDIT_OX)

/* globals accessed in settings.c */
int tabwidth = 4;
int smarttab = 1;

static gboolean script_electric_enter (windata_t *vwin);
static gboolean script_tab_handler (windata_t *vwin, GdkModifierType mods);
static gboolean 
script_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p);
static gchar *textview_get_current_line_with_newline (GtkWidget *view);

void text_set_cursor (GtkWidget *w, GdkCursorType cspec)
{
    GdkWindow *win = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
                                              GTK_TEXT_WINDOW_TEXT);

    if (cspec == 0) {
	gdk_window_set_cursor(win, NULL);
    } else {
	GdkCursor *cursor = gdk_cursor_new(cspec);

	gdk_window_set_cursor(win, cursor);
	gdk_cursor_unref(cursor);
    } 
}

void cursor_to_top (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)); 
    GtkTextIter start;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->text), 
				 mark, 0.0, FALSE, 0, 0);
}

void cursor_to_mark (windata_t *vwin, GtkTextMark *mark)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)); 
    GtkTextIter iter;

    gtk_text_buffer_get_iter_at_mark(buf, &iter, mark);
    gtk_text_buffer_place_cursor(buf, &iter);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->text), 
				 mark, 0.0, TRUE, 0, 0.1);
}

static void get_char_width_and_height (GtkWidget *widget,
				       int *width,
				       int *height)
{
    PangoLayout *pl;
    PangoContext *pc;
    GtkRcStyle *style;
    int w = 0, h = 0;

    pc = gtk_widget_get_pango_context(widget);
    style = gtk_widget_get_modifier_style(widget);
    pango_context_set_font_description(pc, style->font_desc);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, "X", -1);
    pango_layout_get_pixel_size(pl, &w, &h);
    g_object_unref(G_OBJECT(pl));

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

gchar *textview_get_selection_or_all (GtkWidget *view,
				      int *sel)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), NULL);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    if (tbuf == NULL) {
	return NULL;
    }

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	*sel = 1;
    } else {
	*sel = 0;
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
    }

    return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
}

int textview_set_text (GtkWidget *view, const gchar *text)
{
    GtkTextBuffer *tbuf;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(view), 1);

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    g_return_val_if_fail(tbuf != NULL, 1);

    if (text != NULL) {
	gtk_text_buffer_set_text(tbuf, text, -1);
    } else {
	gtk_text_buffer_set_text(tbuf, "", -1);
    }

    return 0;
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

void text_paste (GtkWidget *w, windata_t *vwin)
{
    gchar *undo_buf = textview_get_text(vwin->text);
    gchar *old;

    old = g_object_get_data(G_OBJECT(vwin->text), "undo");
    g_free(old);

    g_object_set_data(G_OBJECT(vwin->text), "undo", undo_buf);

    gtk_text_buffer_paste_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)),
				    gtk_clipboard_get(GDK_NONE),
				    NULL, TRUE);
}

void text_undo (GtkWidget *w, windata_t *vwin)
{
    gchar *old = NULL;

    if (vwin->sbuf != NULL) {
	if (gtk_source_buffer_can_undo(vwin->sbuf)) {
	    gtk_source_buffer_undo(vwin->sbuf);
	} else {
	    errbox(_("No undo information available"));
	}
	return;
    }
    
    old = g_object_steal_data(G_OBJECT(vwin->text), "undo");

    if (old == NULL) {
	errbox(_("No undo information available"));
    } else {
	GtkTextBuffer *buf;
	GtkTextIter start, end;
	GtkTextMark *ins;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	ins = gtk_text_buffer_get_insert(buf);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_delete(buf, &start, &end);

	gtk_text_buffer_insert(buf, &start, old, -1);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->text), 
				     ins, 0.0, TRUE, 0.1, 0.0);
	g_free(old);
    }
}

int text_can_undo (windata_t *vwin)
{
    if (vwin->sbuf != NULL) {
	return gtk_source_buffer_can_undo(vwin->sbuf);
    } else {
	gchar *old = g_object_get_data(G_OBJECT(vwin->text), "undo");

	return old != NULL;
    }
}

static void strip_CRLF (char *s)
{
    int n = strlen(s);

    if (n >= 2 && s[n-2] == '\r') {
	s[n-2] = '\n';
	s[n-1] = '\0';
    }
}

static int source_buffer_load_file (GtkSourceBuffer *sbuf, 
				    int role, FILE *fp)
{
    char fline[MAXSTR];
    gchar *chunk = NULL;
    GtkTextIter iter;
    int i = 0;

    gtk_source_buffer_begin_not_undoable_action(sbuf);

    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), "", -1);
    gtk_text_buffer_get_iter_at_offset(GTK_TEXT_BUFFER(sbuf), &iter, 0);

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

	strip_CRLF(chunk);

	gtk_text_buffer_insert(GTK_TEXT_BUFFER(sbuf), &iter, chunk, -1);
	memset(fline, 0, sizeof fline);

	if (chunk != fline) {
	    g_free(chunk);
	    chunk = NULL;
	}
    }

    gtk_source_buffer_end_not_undoable_action(sbuf);
    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sbuf), FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(sbuf), &iter);
    gtk_text_buffer_place_cursor(GTK_TEXT_BUFFER(sbuf), &iter);

    return 0;
}

static int source_buffer_load_buf (GtkSourceBuffer *sbuf, const char *buf)
{
    char line[MAXLINE];
    GtkTextIter iter;   

    gtk_source_buffer_begin_not_undoable_action(sbuf);
    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), "", -1);
    gtk_text_buffer_get_iter_at_offset(GTK_TEXT_BUFFER(sbuf), &iter, 0);

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	gtk_text_buffer_insert(GTK_TEXT_BUFFER(sbuf), &iter, line, -1);
    }

    bufgets_finalize(buf);

    gtk_source_buffer_end_not_undoable_action(sbuf);
    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sbuf), FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(sbuf), &iter);
    gtk_text_buffer_place_cursor(GTK_TEXT_BUFFER(sbuf), &iter);

    return 0;
}

#ifdef USE_GTKSOURCEVIEW_2

static void sourceview_apply_language (windata_t *vwin)
{
    GtkSourceLanguageManager *lm; 
    GtkSourceLanguage *lang = NULL;
    const char *id = NULL;

    lm = g_object_get_data(G_OBJECT(vwin->sbuf), "languages-manager");

    if (vwin->role == EDIT_GP) {
	id = "gnuplot";
    } else if (vwin->role == EDIT_R) {
	id = "r";
    } else if (vwin->role == EDIT_OX) {
	id = "cpp";
    } else {
	id = "gretl";
    }

    lang = gtk_source_language_manager_get_language(lm, id);
    if (lang == NULL) {
	fprintf(stderr, "*** gtksourceview: lang is NULL for id='%s'\n", id);
    } else {
	gtk_source_buffer_set_language(vwin->sbuf, lang);
    }
}

#else /* use gtksourceview-1.0 API */

static void sourceview_apply_language (windata_t *vwin)
{
    GtkSourceLanguagesManager *lm; 
    GtkSourceLanguage *lang = NULL;
    const char *mtype = NULL;

    lm = g_object_get_data(G_OBJECT(vwin->sbuf), "languages-manager");

    if (vwin->role == EDIT_GP) {
	mtype = "application/x-gnuplot";
    } else if (vwin->role == EDIT_R) {
	mtype = "text/x-R";
    } else if (vwin->role == EDIT_OX) {
	mtype = "text/x-c++src";
    } else {
	mtype = "application/x-gretlscript";
    }

    lang = gtk_source_languages_manager_get_language_from_mime_type(lm, mtype);
    if (lang == NULL) {
	g_object_set(G_OBJECT(vwin->sbuf), "highlight", FALSE, NULL);
    } else {
	g_object_set(G_OBJECT(vwin->sbuf), "highlight", TRUE, NULL);
	gtk_source_buffer_set_language(vwin->sbuf, lang);
    }
}

#endif

static void 
real_sourceview_insert (windata_t *vwin, const char *fname, const char *buf)
{
    FILE *fp = NULL;

    if (fname != NULL) {
	fp = gretl_fopen(fname, "rb");
	if (fp == NULL) {
	    file_read_errbox(fname);
	    return;
	}
    }

    sourceview_apply_language(vwin);

    if (fp != NULL) {
	source_buffer_load_file(vwin->sbuf, vwin->role, fp);
	fclose(fp);
    } else {
	source_buffer_load_buf(vwin->sbuf, buf);
    }
}

void sourceview_insert_file (windata_t *vwin, const char *fname)
{
    real_sourceview_insert(vwin, fname, NULL);
}

void sourceview_insert_buffer (windata_t *vwin, const char *buf)
{
    real_sourceview_insert(vwin, NULL, buf);
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

#define script_editing(r) (r == EDIT_SCRIPT || r == EDIT_FUNC_CODE)

/* Special keystrokes in script window: Ctrl-Return sends the current
   line for execution; Ctrl-R sends the whole script for execution
   (i.e. is the keyboard equivalent of the "execute" icon).
*/

static gint 
script_key_handler (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods = widget_get_pointer_mask(w);
    gboolean ret = FALSE;

    if (mods & GDK_CONTROL_MASK) {
	if (key->keyval == GDK_r)  {
	    do_run_script(w, vwin);
	    ret = TRUE;
	} else if (key->keyval == GDK_Return) {
	    gchar *str = textview_get_current_line_with_newline(w);

	    if (str != NULL) {
		if (!string_is_blank(str)) {
		    run_script_fragment(vwin, str);
		}
		g_free(str);
	    }
	    ret = TRUE;
	}
    } else {
	if (key->keyval == GDK_F1) {
	    set_window_help_active(vwin);
	    interactive_script_help(NULL, NULL, vwin);
	} else if (script_editing(vwin->role)) {    
	    if (key->keyval == GDK_Return) {
		ret = script_electric_enter(vwin);
	    } else if (tabkey(key->keyval)) {
		ret = script_tab_handler(vwin, mods);
	    }
	}
    } 

    return ret;
}

#define gretl_script_role(r) (r == EDIT_SCRIPT || \
			      r == VIEW_SCRIPT || \
			      r == EDIT_FUNC_CODE)

void create_source (windata_t *vwin, int hsize, int vsize, 
		    gboolean editable)
{
#ifdef USE_GTKSOURCEVIEW_2
    GtkSourceLanguageManager *lm = gtk_source_language_manager_new();
#else
    GtkSourceLanguagesManager *lm = gtk_source_languages_manager_new();
#endif
    GtkSourceBuffer *sbuf;
    int cw;
    
    sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(NULL));
    g_object_ref(lm);
    g_object_set_data_full(G_OBJECT(sbuf), "languages-manager",
			   lm, (GDestroyNotify) g_object_unref); 
    g_object_unref(lm); 

#ifdef USE_GTKSOURCEVIEW_2
    gtk_source_buffer_set_highlight_matching_brackets(sbuf, TRUE);
#else
    gtk_source_buffer_set_check_brackets(sbuf, TRUE);
#endif

    vwin->text = gtk_source_view_new_with_buffer(sbuf);
    vwin->sbuf = sbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->text), GTK_WRAP_NONE);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->text), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->text), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->text), fixed_font);

    cw = get_char_width(vwin->text);
    hsize *= cw;
    hsize += 48;

    set_source_tabs(vwin->text, cw);

    gtk_window_set_default_size(GTK_WINDOW(vwin->main), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->text), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->text), editable);

    if (gretl_script_role(vwin->role)) {
	g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
			 G_CALLBACK(script_key_handler), vwin);
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(script_popup_handler), 
			 vwin);
	g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
			 G_CALLBACK(interactive_script_help), vwin);
    } else if (foreign_script_role(vwin->role)) {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(script_popup_handler), 
			 vwin);
    } else if (vwin->role == VIEW_LOG) {
	g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
			 G_CALLBACK(interactive_script_help), vwin);
    }	
}

void text_zoom (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    windata_t *vwin = (windata_t *) data;
    GtkTextBuffer *tbuf;
    GtkTextTagTable *table;
    static PangoFontDescription *hpf;
    static gint fsize;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    table = gtk_text_buffer_get_tag_table(tbuf);

    if (hpf == NULL) {
	hpf = pango_font_description_copy(fixed_font);
	fsize = pango_font_description_get_size(hpf) / PANGO_SCALE;
    }

    if (!strcmp(s, "ZoomIn")) {
	fsize++;
    } else if (!strcmp(s, "ZoomOut")) {
	fsize--;
    } 

    pango_font_description_set_size(hpf, fsize * PANGO_SCALE);
    gtk_widget_modify_font(vwin->text, hpf);
}

static GtkTextTagTable *gretl_tags_new (void)
{
    GtkTextTagTable *table;
    GtkTextTag *tag;

    table = gtk_text_tag_table_new(); 

    tag = gtk_text_tag_new("bluetext");
    g_object_set(tag, "foreground", "blue", NULL);
    gtk_text_tag_table_add(table, tag);
    
    tag = gtk_text_tag_new("redtext");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("title");
    g_object_set(tag, "justification", GTK_JUSTIFY_CENTER,
		 "pixels_above_lines", 15,
		 "family", "sans",
		 "size", 15 * PANGO_SCALE, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("sansbold");
    g_object_set(tag, "family", "sans", 
		 "weight", PANGO_WEIGHT_BOLD, 
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("italic");
    g_object_set(tag, "family", "sans",
		 "style", PANGO_STYLE_ITALIC, 
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("replaceable");
    g_object_set(tag, "family", "sans",
		 "style", PANGO_STYLE_ITALIC, 
		 NULL);

    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("superscript");
    g_object_set(tag, "style", PANGO_STYLE_NORMAL,
		 "rise", 4 * PANGO_SCALE,
		 "size", 8 * PANGO_SCALE,
		 NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("subscript");
    g_object_set(tag, "family", "sans",
		 "style", PANGO_STYLE_ITALIC,
		 "rise", -3 * PANGO_SCALE,
		 "size", 8 * PANGO_SCALE,
		 NULL);
    gtk_text_tag_table_add(table, tag);
		 
    tag = gtk_text_tag_new("literal");
    g_object_set(tag, "family", "monospace", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("text");
    g_object_set(tag, "family", "sans", NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("indented");
    g_object_set(tag, "left_margin", 16, "indent", -12, NULL);
    gtk_text_tag_table_add(table, tag);

    tag = gtk_text_tag_new("code");
    g_object_set(tag, "family", "monospace", NULL);
    gtk_text_tag_table_add(table, tag);

    return table;
}

GtkTextBuffer *gretl_text_buf_new (void)
{
    static GtkTextTagTable *tags = NULL;
    GtkTextBuffer *tbuf; 

    if (tags == NULL) {
	tags = gretl_tags_new();
    }

    tbuf = gtk_text_buffer_new(tags);

    return tbuf;
}

static void 
real_textview_add_colorized (GtkWidget *view, const char *buf,
			     int append, int trim)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter; 
    int nextcolor, thiscolor = PLAIN_TEXT;
    int in_comment = 0;
    char readbuf[MAXSTR];
    int i = 0;

    g_return_if_fail(GTK_IS_TEXT_VIEW(view));

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));

    if (append) {
	gtk_text_buffer_get_end_iter(tbuf, &iter);
    } else {
	gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    }

    bufgets_init(buf);

    while (bufgets(readbuf, sizeof readbuf, buf)) {
	if (trim && i++ < 2) {
	    continue;
	}

	if (ends_with_backslash(readbuf)) {
	    nextcolor = BLUE_TEXT;
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
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
						     readbuf, -1,
						     "bluetext", NULL);
	} else {
	    gtk_text_buffer_insert(tbuf, &iter, readbuf, -1);
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

void textview_insert_file (windata_t *vwin, const char *fname)
{
    FILE *fp;
    GtkTextBuffer *tbuf;
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    char fline[MAXSTR], *chunk;
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

	if (chunk != fline) {
	    g_free(chunk);
	}

	thiscolor = nextcolor;
	memset(fline, 0, sizeof fline);
    }

    fclose(fp);
}

void textview_insert_from_tempfile (windata_t *vwin, PRN *prn)
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;    
    char readbuf[MAXSTR];
    FILE *fp;

    fp = gretl_print_read_tempfile(prn);
    if (fp == NULL) return;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, -1);
    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf, fp)) {
	gtk_text_buffer_insert(tbuf, &iter, readbuf, -1);
	memset(readbuf, 0, sizeof readbuf);
    }

    gretl_print_stop_tempfile_read(prn, fp);

    while (gtk_events_pending()) {
        gtk_main_iteration();
    }
}

static void insert_link (GtkTextBuffer *tbuf, GtkTextIter *iter, 
			 const char *text, gint page, 
			 const char *indent)
{
    GtkTextTagTable *tab = gtk_text_buffer_get_tag_table(tbuf);
    GtkTextTag *tag;
    gchar tagname[32];

    if (page == GUIDE_PAGE) {
	strcpy(tagname, "tag:guide");
    } else if (page == SCRIPT_PAGE) {
	strcpy(tagname, text);
    } else {
	sprintf(tagname, "tag:p%d", page);
    }

    tag = gtk_text_tag_table_lookup(tab, tagname);

    if (tag == NULL) {
	if (page == GUIDE_PAGE) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", 
					     "family", "sans", NULL);
	} else if (page == SCRIPT_PAGE) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", 
					     "family", "monospace", NULL);
	    g_object_set_data_full(G_OBJECT(tag), "fname", g_strdup(text), 
				   g_free);
	} else if (indent != NULL) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", 
					     "left_margin", 30, NULL);
	} else {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", NULL);
	}
	g_object_set_data(G_OBJECT(tag), "page", GINT_TO_POINTER(page));
    } 

    gtk_text_buffer_insert_with_tags(tbuf, iter, text, -1, tag, NULL);
}

static void insert_xlink (GtkTextBuffer *tbuf, GtkTextIter *iter, 
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
	if (gfr) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", 
					     "family", "sans", NULL);
	} else if (indent != NULL) {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", 
					     "left_margin", 30, NULL);
	} else {
	    tag = gtk_text_buffer_create_tag(tbuf, tagname, "foreground", "blue", NULL);
	}
	g_object_set_data(G_OBJECT(tag), "page", GINT_TO_POINTER(page));
	g_object_set_data(G_OBJECT(tag), "xref", GINT_TO_POINTER(1));
    } 

    gtk_text_buffer_insert_with_tags(tbuf, iter, text, -1, tag, NULL);
}

static void link_open_script (GtkTextTag *tag)
{
    const char *fname = g_object_get_data(G_OBJECT(tag), "fname");
    char fullname[MAXLEN];

    sprintf(fullname, "%sscripts%cmisc%c%s", paths.gretldir, 
	    SLASH, SLASH, fname);
    view_file(fullname, 0, 0, 78, 370, VIEW_SCRIPT);
}

static int object_get_int (gpointer p, const char *key)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(p), key));
}

static void follow_if_link (GtkWidget *tview, GtkTextIter *iter, gpointer p)
{
    GSList *tags = NULL, *tagp = NULL;

    tags = gtk_text_iter_get_tags(iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	GtkTextTag *tag = tagp->data;
	gint page = object_get_int(tag, "page");
	gint xref = object_get_int(tag, "xref");
	gint en = GPOINTER_TO_INT(p);

	if (page != 0 || xref != 0) {
	    if (page == GUIDE_PAGE) {
		display_pdf_help(NULL);
	    } else if (page == SCRIPT_PAGE) {
		link_open_script(tag);
	    } else {
		int role = object_get_int(tview, "role");

		if (role == FUNCS_HELP) {
		    if (xref) {
			command_help_callback(page, en); 
		    } else {
			function_help_callback(page);
		    }
		} else {
		    if (xref) {
			function_help_callback(page);
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
				  gpointer p)
{
    GtkTextIter iter;
    GtkTextBuffer *tbuf;

    switch (ev->keyval) {
    case GDK_Return: 
    case GDK_KP_Enter:
	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));
	gtk_text_buffer_get_iter_at_mark(tbuf, &iter, 
					 gtk_text_buffer_get_insert(tbuf));
	follow_if_link(tview, &iter, p);
	break;
    default:
	break;
    }

    return FALSE;
}

/* Help links can be activated by clicking */

static gboolean cmdref_event_after (GtkWidget *tview, GdkEvent *ev,
				    gpointer p)
{
    GtkTextIter start, end, iter;
    GtkTextBuffer *buffer;
    GdkEventButton *event;
    gint x, y;

    if (ev->type != GDK_BUTTON_RELEASE) {
	return FALSE;
    }

    event = (GdkEventButton *) ev;

    if (event->button != 1)
	return FALSE;

    buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tview));

    /* don't follow a link if the user has selected something */
    gtk_text_buffer_get_selection_bounds(buffer, &start, &end);
    if (gtk_text_iter_get_offset(&start) != gtk_text_iter_get_offset(&end))
	return FALSE;

    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);

    gtk_text_view_get_iter_at_location(GTK_TEXT_VIEW(tview), &iter, x, y);

    follow_if_link(tview, &iter, p);

    return FALSE;
}

static GdkCursor *hand_cursor = NULL;
static GdkCursor *regular_cursor = NULL;

static void
set_cursor_if_appropriate (GtkTextView *tview, gint x, gint y)
{
    static gboolean hovering_over_link = FALSE;
    GSList *tags = NULL, *tagp = NULL;
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    gboolean hovering = FALSE;

    tbuf = gtk_text_view_get_buffer(tview);
    gtk_text_view_get_iter_at_location(tview, &iter, x, y);
  
    tags = gtk_text_iter_get_tags(&iter);

    for (tagp = tags; tagp != NULL; tagp = tagp->next) {
	GtkTextTag *tag = tagp->data;
	gint page = object_get_int(tag, "page");
	gint xref = object_get_int(tag, "xref");

	if (page != 0 || xref != 0) { 
	    hovering = TRUE;
	    break;
        }
    }

    if (hovering != hovering_over_link) {
	hovering_over_link = hovering;

	if (hovering_over_link) {
	    gdk_window_set_cursor(gtk_text_view_get_window(tview, GTK_TEXT_WINDOW_TEXT), 
				  hand_cursor);
	} else {
	    gdk_window_set_cursor(gtk_text_view_get_window(tview, GTK_TEXT_WINDOW_TEXT), 
				  regular_cursor);
	}
    }

    if (tags) {
	g_slist_free(tags);
    }
}

static gboolean 
cmdref_motion_notify (GtkWidget *tview, GdkEventMotion *event)
{
    gint x, y;

    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  event->x, event->y, &x, &y);
    set_cursor_if_appropriate(GTK_TEXT_VIEW(tview), x, y);
    widget_get_pointer_mask(tview);

    return FALSE;
}

static gboolean
cmdref_visibility_notify (GtkWidget *tview,  GdkEventVisibility *e)
{
    gint wx, wy, bx, by;

    widget_get_pointer_info(tview, &wx, &wy, NULL);
    gtk_text_view_window_to_buffer_coords(GTK_TEXT_VIEW(tview), 
					  GTK_TEXT_WINDOW_WIDGET,
					  wx, wy, &bx, &by);
    set_cursor_if_appropriate(GTK_TEXT_VIEW(tview), bx, by);

    return FALSE;
}

static void maybe_connect_help_signals (windata_t *hwin, int en)
{
    int done = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(hwin->text), 
						 "sigs_connected"));

    if (hand_cursor == NULL) {
	hand_cursor = gdk_cursor_new(GDK_HAND2);
    }

    if (regular_cursor == NULL) {
	regular_cursor = gdk_cursor_new(GDK_XTERM);
    }    

    if (!done) {
	gpointer en_ptr = GINT_TO_POINTER(en);

	g_signal_connect(hwin->text, "key-press-event", 
			 G_CALLBACK(cmdref_key_press), en_ptr);
	g_signal_connect(hwin->text, "event-after", 
			 G_CALLBACK(cmdref_event_after), en_ptr);
	g_signal_connect(hwin->text, "motion-notify-event", 
			 G_CALLBACK(cmdref_motion_notify), NULL);
	g_signal_connect(hwin->text, "visibility-notify-event", 
			 G_CALLBACK(cmdref_visibility_notify), NULL);
	g_object_set_data(G_OBJECT(hwin->text), "sigs_connected", 
			  GINT_TO_POINTER(1));
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
	g_object_set_data(G_OBJECT(hwin->text), "tabs_set", GINT_TO_POINTER(1));
    }
}

static void cmdref_index_page (windata_t *hwin, GtkTextBuffer *tbuf, int en)
{
    const char *header = N_("Gretl Command Reference");
    GtkTextIter iter;
    int jmax = 6;
    int i, j, k;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     (en)? header : _(header), -1,
					     "title", NULL);
    gtk_text_buffer_insert(tbuf, &iter, "\n\n", -1);

    j = 1;

    for (i=1; i<NC; i++) {
	const char *word;

	if (HIDDEN_COMMAND(i)) {
	    continue;
	}

	word = gretl_command_word(i);
	insert_link(tbuf, &iter, gretl_command_word(i), i, NULL);
	if (j++ % jmax == 0) {
	    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
	} else {
	    int n = 10 - strlen(word);

	    for (k=0; k<n; k++) {
		gtk_text_buffer_insert(tbuf, &iter, " ", -1);
	    }
	}
    }

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->text), tbuf);

    maybe_connect_help_signals(hwin, en);
    maybe_set_help_tabs(hwin);
}

static void funcref_index_page (windata_t *hwin, GtkTextBuffer *tbuf, int en)
{
    const char *header = N_("Gretl Function Reference");
    const gchar *s;
    GtkTextIter iter;
    char funword[12];
    int llen_max = 5;
    int llen;
    int i, j, n;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    gtk_text_buffer_insert_with_tags_by_name(tbuf, &iter,
					     (en)? header : _(header), -1,
					     "title", NULL);
    gtk_text_buffer_insert(tbuf, &iter, "\n\n", -1);

    s = (const gchar *) hwin->data;
    i = 1;
    llen = 0;

    while (*s) {
	if (*s == '\n' && *(s+1) == '#' && *(s+2) != '\0') {
	    if (*(s+2) == '#') {
		/* category divider */
		if (i > 1) {
		    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
		    if (llen < llen_max) {
			gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
			llen = 0;
		    }
		}
		s += 2;
	    } else if (sscanf(s + 2, "%10s", funword)) {
		/* function name */
		insert_link(tbuf, &iter, funword, i, NULL);
		if (++llen == llen_max) {
		    gtk_text_buffer_insert(tbuf, &iter, "\n", -1);
		    llen = 0;
		} else {
		    n = 12 - strlen(funword);
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

    maybe_connect_help_signals(hwin, en);
    maybe_set_help_tabs(hwin);
}

static void push_backpage (GtkWidget *w, int pg)
{
    gpointer p = GINT_TO_POINTER(pg);

    g_object_set_data(G_OBJECT(w), "backpage", p);
}

static int pop_backpage (GtkWidget *w)
{
    gpointer p = g_object_get_data(G_OBJECT(w), "backpage");

    return GPOINTER_TO_INT(p);
}

static gint help_popup_click (GtkWidget *w, gpointer p)
{
    windata_t *hwin = (windata_t *) p;
    int action = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));
    int en = (hwin->role == CLI_HELP_EN);
    int page = 0;

    if (action == 2) {
	page = pop_backpage(hwin->text); 
    }

    if (hwin->role == FUNCS_HELP) {
	function_help_callback(page);
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
    int i, imin = 0;

    if (hwin->active_var == 0) {
	imin = 1;
    }

    for (i=imin; i<2; i++) {
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_object_set_data(G_OBJECT(item), "action", GINT_TO_POINTER(i+1));
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
    GdkModifierType mods = widget_get_pointer_mask(w);

    if (mods & GDK_BUTTON3_MASK) {
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

gchar *textview_get_current_line (GtkWidget *view)
{
    GtkTextBuffer *buf;
    GtkTextIter start, end;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
    gtk_text_buffer_get_iter_at_mark(buf, &start, 
				     gtk_text_buffer_get_insert(buf));
    gtk_text_iter_set_line_offset(&start, 0);
    gtk_text_buffer_get_iter_at_mark(buf, &end, 
				     gtk_text_buffer_get_insert(buf));
    gtk_text_iter_forward_to_line_end(&end);

    return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}

static gchar *textview_get_current_line_with_newline (GtkWidget *view)
{
    gchar *s = textview_get_current_line(view);

    if (s != NULL && *s != '\0' && s[strlen(s)-1] != '\n') {
	gchar *tmp = g_strdup_printf("%s\n", s);

	g_free(s);
	s = tmp;
    }

    return s;
}

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

static void get_cmdword (const char *s, char *word)
{
    if (sscanf(s, "%*s <- %8s", word) != 1) {
	sscanf(s, "%8s", word);
    }
#if 0
    if (*word == '\0') {
	int i;

	fprintf(stderr, "get_cmdword: s = '%s'\n", s);
	for (i=0; i<strlen(s); i++) {
	    fprintf(stderr, "s[%d] = %d (%c)\n",
		    i, (int) s[i], s[i]);
	}
    }
#endif
}

#define bare_quote(p,s)   (*p == '"' && (p-s==0 || *(p-1) != '\\'))
#define starts_comment(p) (*p == '/' && *(p+1) == '*')
#define ends_comment(p)   (*p == '*' && *(p+1) == '/')

static void check_for_comment (const char *s, int *incomm)
{
    int commbak = *incomm;
    const char *p = s;
    int quoted = 0;

    while (*p) {
	if (!quoted && !*incomm && *p == '#') {
	    break;
	}
	if (!*incomm && bare_quote(p, s)) {
	    quoted = !quoted;
	}
	if (!quoted) {
	    if (starts_comment(p)) {
		*incomm = 1;
		p += 2;
	    } else if (ends_comment(p)) {
		*incomm = 0;
		p += 2;
		p += strspn(p, " ");
	    }
	}
	if (*p) {
	    p++;
	}
    }

    if (*incomm && commbak) {
	/* on the second or subsequent line of a multiline
	   comment */
	*incomm = 2;
    }
}

static void normalize_indent (GtkTextBuffer *tbuf, 
			      const gchar *buf,
			      GtkTextIter *start,
			      GtkTextIter *end)
{
    int this_indent = 0;
    int next_indent = 0;
    char word[9], line[1024];
    const char *ins;
    int incomment = 0;
    int i, nsp;

    if (buf == NULL) {
	return;
    }    

    gtk_text_buffer_delete(tbuf, start, end);

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (string_is_blank(line)) {
	    gtk_text_buffer_insert(tbuf, start, line, -1);
	    continue;
	}
	check_for_comment(line, &incomment);
#if 0
	if (incomment) {
	    /* in multiline comment */
	    gtk_text_buffer_insert(tbuf, start, line, -1);
	    continue;
	}
#endif
	ins = line + strspn(line, " \t");
	if (!incomment) {
	    *word = '\0';
	    get_cmdword(ins, word);
	    adjust_indent(word, &this_indent, &next_indent);
	}
	nsp = this_indent * tabwidth;
	if (incomment == 2) {
	    nsp += 3;
	}
	for (i=0; i<nsp; i++) {
	    gtk_text_buffer_insert(tbuf, start, " ", -1);
	}
	gtk_text_buffer_insert(tbuf, start, ins, -1);
    }

    bufgets_finalize(buf); 
}

static void auto_indent_script (GtkWidget *w, windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;
    gchar *buf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);
    buf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    normalize_indent(tbuf, buf, &start, &end);
    g_free(buf);
}

static void indent_region (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;

    if (smarttab) {
	normalize_indent(tb->buf, tb->chunk, &tb->start, &tb->end);
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

void script_tabs_dialog (GtkWidget *w, windata_t *vwin)
{
    const char *title = _("gretl: configure tabs");
    const char *spintxt = _("Spaces per tab");
    const char *opt = _("Use \"smart\" tabs");
    int tsp = tabwidth;
    int smt = smarttab;
    int resp;

    resp = checks_dialog(title, NULL, &opt, 1, &smt,
			 0, NULL, /* no radio buttons */
			 &tsp, spintxt, 2, 8, 0);

    if (resp != GRETL_CANCEL) {
	tabwidth = tsp;
	smarttab = smt;
    }
}

static void exec_script_text (GtkWidget *w, gpointer p)
{
    struct textbit *tb = (struct textbit *) p;

    run_script_fragment(tb->vwin, tb->chunk);
    tb->chunk = NULL; /* will be freed already */
}

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

/* Given what's presumed to be a start-of-line iter, find how many
   leading spaces are on the line, counting tabs as multiple spaces.
*/

static int leading_spaces_at_iter (GtkTextBuffer *tbuf, GtkTextIter *start)
{
    GtkTextIter end = *start;
    gchar *s;
    int n = 0;

    gtk_text_iter_forward_to_line_end(&end);
    s = gtk_text_buffer_get_text(tbuf, start, &end, FALSE);
    if (s != NULL) {
	n = count_leading_spaces(s);
	g_free(s);
    }

    return n;
}

#define TABDEBUG 0

static int incremental_leading_spaces (const char *prevword,
				       const char *thisword)
{
    int this_indent = 0;
    int next_indent = 0;

    if (*prevword != '\0') {
	int prev_indent = 0;

	adjust_indent(prevword, &this_indent, &next_indent);
#if TABDEBUG > 1
	fprintf(stderr, "adjust_indent 1: this=%d, next=%d\n",
		this_indent, next_indent);
#endif
	prev_indent = this_indent;
	if (*thisword != '\0') {
	    adjust_indent(thisword, &this_indent, &next_indent);
#if TABDEBUG > 1
	    fprintf(stderr, "adjust_indent 2: this=%d\n", this_indent);
#endif
	    this_indent -= prev_indent;
#if TABDEBUG > 1
	    fprintf(stderr, "adjust_indent 2: this=%d\n", this_indent);
#endif
	} else {
	    this_indent = next_indent - this_indent;
	}
    }

    return this_indent * tabwidth;
}

static int line_continues (const gchar *s)
{
    int i, n = strlen(s);

    for (i=n-1; i>=0; i--) {
	if (s[i] != ' ' && s[i] != '\t') {
	    return (s[i] == '\\');
	}
    } 

    return 0;
}

static int get_word_and_cont (const char *s, char *word, int *contd)
{
    /* don't move onto next line */
    if (*s != '\n' && *s != '\r') {
	if (contd != NULL) {
	    *contd = line_continues(s);
	}
	s += strspn(s, " \t");
	if (sscanf(s, "%*s <- %8s", word) != 1) {
	    sscanf(s, "%8s", word);
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
   to find this
 */

static char *get_previous_line_start_word (char *word, 
					   GtkTextBuffer *tbuf,
					   GtkTextIter iter,
					   int *leadspace,
					   int *contd)
{
    GtkTextIter end, prev = iter;
    int *pcont = contd;
    gchar *s;

    *word = '\0';

    while (*word == '\0' && gtk_text_iter_backward_line(&prev)) {
	end = prev;
	if (gtk_text_iter_forward_to_line_end(&end)) {
	    s = gtk_text_buffer_get_text(tbuf, &prev, &end, FALSE);
	    if (s != NULL) {
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
	    *leadspace = leading_spaces_at_iter(tbuf, &prev);
	}
    }

    return word;
}

/* Is the insertion point at the start of a line, or in a white-space
   field to the left of any non-space characters?  If so, we'll trying
   inserting a "smart" soft tab in response to the Tab key.
*/

static int maybe_insert_smart_tab (windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    GtkTextMark *mark;
    GtkTextIter start, end;
    gchar *chunk = NULL;
    int pos = 0, ret = 0;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	return 0;
    }

    mark = gtk_text_buffer_get_insert(tbuf);
    gtk_text_buffer_get_iter_at_mark(tbuf, &end, mark);
    pos = gtk_text_iter_get_line_offset(&end);

    if (pos == 0) {
	ret = 1;
    } else {
	start = end;
	gtk_text_iter_set_line_offset(&start, 0);
	chunk = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
	ret = strspn(chunk, " \t") == strlen(chunk);
    }

    if (ret) {
	GtkTextIter prev = start;
	char *s, thisword[9];
	char prevword[9];
	int i, nsp = 0, contd = 0;

	*prevword = *thisword = '\0';

	s = textview_get_current_line(vwin->text);
	if (s != NULL) {
	    sscanf(s, "%8s", thisword);
	    g_free(s);
	} 

	get_previous_line_start_word(prevword, tbuf, prev, &nsp, &contd);

	if (contd) {
	    nsp += 2;
	} else {
	    nsp += incremental_leading_spaces(prevword, thisword);
	}

	if (pos > 0) {
	    gtk_text_buffer_delete(tbuf, &start, &end);
	}
	for (i=0; i<nsp; i++) {
	    gtk_text_buffer_insert(tbuf, &start, " ", -1);
	}	
	if (pos > 0) {
	    s = chunk + strspn(chunk, " \t");
	    gtk_text_buffer_insert(tbuf, &start, s, -1);
	}
    }

    if (chunk != NULL) {
	g_free(chunk);
    }

    return ret;
}

/* On "Enter" in script editing, try to compute the correct indent
   level for the current line, and make an adjustment if it's not
   already in place
*/

static gboolean script_electric_enter (windata_t *vwin)
{
    char *s;

    if (!smarttab) {
	return FALSE;
    }

    s = textview_get_current_line(vwin->text);

    if (s == NULL) {
	return FALSE;
    } else if (*s == '\0') {
	g_free(s);
	return FALSE;
    } else {
	GtkTextBuffer *tbuf;
	GtkTextMark *mark;
	GtkTextIter start, end;
	char thisword[9];
	char prevword[9];
	int diff, nsp, incr;
	int targsp = 0, contd = 0;

	*thisword = *prevword = '\0';

	sscanf(s, "%8s", thisword);
	nsp = count_leading_spaces(s);
	g_free(s);

	if (*thisword == '\0') {
	    return FALSE;
	}

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	mark = gtk_text_buffer_get_insert(tbuf);
	gtk_text_buffer_get_iter_at_mark(tbuf, &start, mark);
	gtk_text_iter_set_line_offset(&start, 0);

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
	    int i;

	    diff = -diff;
	    for (i=0; i<diff; i++) {
		gtk_text_buffer_insert(tbuf, &start, " ", -1);
	    }
	}
    }

    return FALSE;
}

/* handler for the user pressing the Tab key when editing a script */

static gboolean script_tab_handler (windata_t *vwin, GdkModifierType mods)
{
    struct textbit *tb;
    gboolean ret = FALSE;

    g_return_val_if_fail(GTK_IS_TEXT_VIEW(vwin->text), FALSE);

    if (smarttab && !(mods & GDK_SHIFT_MASK)) {
	if (maybe_insert_smart_tab(vwin)) {
	    return TRUE;
	}
    }

    /* do we really want the rest of this? */

    tb = vwin_get_textbit(vwin, AUTO_SELECT_NONE);
    if (tb == NULL) {
	return FALSE;
    }

    if (tb->selected) {
	if (mods & GDK_SHIFT_MASK) {
	    unindent_region(NULL, tb);
	} else {
	    indent_region(NULL, tb);
	}
	ret = TRUE;
    }

    g_free(tb->chunk);
    free(tb);

    return ret;
}

static void line_numbers_cb (GtkWidget *w, windata_t *vwin)
{
    int s = gtk_source_view_get_show_line_numbers(GTK_SOURCE_VIEW(vwin->text));

    gtk_source_view_set_show_line_numbers(GTK_SOURCE_VIEW(vwin->text), !s);
}

#define editing_code(r) (r == EDIT_SCRIPT || r == EDIT_FUNC_CODE)

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
	goto line_nums;
    }

    tb = vwin_get_textbit(vwin, AUTO_SELECT_LINE);
    if (tb == NULL) {
	*ptb = NULL;
	return pmenu;
    }

    tb->commented = text_is_commented(tb->chunk);

    if (tb->commented > 0 && !editing_code(vwin->role)) {
	g_free(tb->chunk);
	free(tb);
	*ptb = NULL;
	goto line_nums;
    }

    *ptb = tb;

    if (tb->commented <= 0 && vwin->role != EDIT_FUNC_CODE) {
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

    if (editing_code(vwin->role) && tb->commented >= 0) {
	/* material is either all commented or all uncommented:
	   allow comment/uncomment option */
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

    if (editing_code(vwin->role)) {
	if (tb->selected) {
	    item = gtk_menu_item_new_with_label(smarttab? 
						_("Auto-indent region") :
						_("Indent region"));
	    g_signal_connect(G_OBJECT(item), "activate",
			     G_CALLBACK(indent_region),
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

	item = gtk_menu_item_new_with_label(_("Auto-indent script"));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(auto_indent_script),
			 vwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
    }

 line_nums:	

    if (GTK_IS_SOURCE_VIEW(vwin->text)) {
	item = gtk_menu_item_new_with_label(_("Toggle line numbers"));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(line_numbers_cb),
			 vwin);
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), item);
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
    GdkModifierType mods = widget_get_pointer_mask(w);

    if (mods & GDK_BUTTON3_MASK) {
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
    INSERT_ITAL,
    INSERT_SUP,
    INSERT_SUB,
    INSERT_TEXT,
    INSERT_PDFLINK,
    INSERT_INPLINK,
    INSERT_GFRLINK
};

static void insert_help_figure (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *fig)
{
    char figfile[FILENAME_MAX];
    GdkPixbuf *pixbuf;

    sprintf(figfile, "%shelpfigs%c%s.png", paths.gretldir,
	    SLASH, fig);

    pixbuf = gdk_pixbuf_new_from_file(figfile, NULL);

    if (pixbuf != NULL) {
	gtk_text_buffer_insert_pixbuf(tbuf, iter, pixbuf);
	g_object_unref(G_OBJECT(pixbuf));
    }
}

static void insert_tagged_text (GtkTextBuffer *tbuf, GtkTextIter *iter,
				const char *s, int ins, const char *indent)
{
    const char *ftag = NULL;

    switch (ins) {
    case INSERT_ITAL:
	ftag = "italic";
	break;
    case INSERT_REPL:
	ftag = "replaceable";
	break;
    case INSERT_LIT:
	ftag = "literal";
	break;
    case INSERT_SUP:
	ftag = "superscript";
	break;
    case INSERT_SUB:
	ftag = "subscript";
	break;
    default:
	break;
    }

    if (ftag != NULL) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 ftag, indent, NULL);
    }
}

static int get_instruction_and_string (const char *p, char *str)
{
    int ins = INSERT_NONE;
    *str = '\0';

    if (!strncmp(p, "ref", 3)) {
	ins = INSERT_REF;
    } else if (!strncmp(p, "xrf", 3)) {
	ins = INSERT_XREF;
    } else if (!strncmp(p, "fig", 3)) {
	ins = INSERT_FIG;
    } else if (!strncmp(p, "itl", 3)) {
	ins = INSERT_ITAL;
    } else if (!strncmp(p, "var", 3)) {
	ins = INSERT_REPL;
    } else if (!strncmp(p, "lit", 3)) {
	ins = INSERT_LIT;
    } else if (!strncmp(p, "sup", 3)) {
	ins = INSERT_SUP;
    } else if (!strncmp(p, "sub", 3)) {
	ins = INSERT_SUB;
    } else if (!strncmp(p, "pdf", 3)) {
	ins = INSERT_PDFLINK;
    } else if (!strncmp(p, "inp", 3)) {
	ins = INSERT_INPLINK;
    } else if (!strncmp(p, "gfr", 3)) {
	ins = INSERT_GFRLINK;
    }

    if (ins != INSERT_NONE) {
	int i = 0;

	p += 5;
	while (*p) {
	    if (*p == '"' && *(p+1) == '>') {
		str[i] = '\0';
		break;
	    } else {
		str[i++] = *p++;
	    }
	}
    }

    return ins;
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

static void
insert_text_with_markup (GtkTextBuffer *tbuf, GtkTextIter *iter,
			 const char *s, int role)
{
    static char targ[128];
    const char *indent = NULL;
    const char *code = NULL;
    const char *p;
    int itarg, ins;

    while ((p = strstr(s, "<"))) {
	int skip = 0;

	if (code != NULL) {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "code", indent, NULL);
	} else {
	    gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, p - s,
						     "text", indent, NULL);
	}

	p++;

	if (*p == '@') {
	    /* "atomic" markup */
	    ins = get_instruction_and_string(p + 1, targ);
	    if (ins == INSERT_REF) {
		if (role == FUNCS_HELP) {
		    itarg = function_help_index_from_word(targ);
		} else {
		    itarg = gretl_command_number(targ);
		}
		insert_link(tbuf, iter, targ, itarg, indent);
	    } else if (ins == INSERT_XREF) {
		if (role == FUNCS_HELP) {
		    itarg = gretl_command_number(targ);
		} else {
		    itarg = function_help_index_from_word(targ);
		}
		insert_xlink(tbuf, iter, targ, itarg, indent);
	    } else if (ins == INSERT_PDFLINK) {
		insert_link(tbuf, iter, targ, GUIDE_PAGE, indent);
	    } else if (ins == INSERT_INPLINK) {
		insert_link(tbuf, iter, targ, SCRIPT_PAGE, indent);
	    } else if (ins == INSERT_GFRLINK) {
		insert_xlink(tbuf, iter, targ, GFR_PAGE, indent);
	    } else if (ins == INSERT_FIG) {
		insert_help_figure(tbuf, iter, targ);
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
	    code = "code";
	    skip = get_code_skip(p + 5);
	} else if (!strncmp(p, "/code", 5)) {
	    code = NULL;
	    skip = 6 + (*(p+6) == '\n');
	} else {
	    /* literal "<" */
	    gtk_text_buffer_insert(tbuf, iter, "<", 1);
	}

	s = p + skip;
    }

    if (code != NULL) {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "code", indent, NULL);
    } else {
	gtk_text_buffer_insert_with_tags_by_name(tbuf, iter, s, -1,
						 "text", indent, NULL);
    }
}

static char *grab_topic_buffer (const char *s)
{
    const char *p = strstr(s, "\n# ");
    char *buf;

    if (p != NULL) {
	buf = g_strndup(s, p - s);
    } else {
	buf = g_strdup(s);
    }

    return buf;
}

/* Pull the appropriate chunk of help text out of the buffer attached
   to the help viewer and display it.  Also set the active_var member
   of hwin to represent the topic displayed.
*/

void set_help_topic_buffer (windata_t *hwin, int hcode, int pos, int en)
{
    GtkTextBuffer *textb;
    GtkTextIter iter;
    char line[256];
    gchar *hbuf;
    char *buf;

    textb = gretl_text_buf_new();

    if (pos == 0) {
	/* no topic selected */
	if (hwin->role == FUNCS_HELP) {
	    funcref_index_page(hwin, textb, en);
	} else {
	    cmdref_index_page(hwin, textb, en);
	}
	cursor_to_top(hwin);
	hwin->active_var = 0;
	return;
    }

    /* OK, pos is non-zero */

    maybe_connect_help_signals(hwin, en);
    maybe_set_help_tabs(hwin);
    
    gtk_text_buffer_get_iter_at_offset(textb, &iter, 0);

    hbuf = (gchar *) hwin->data + pos;

    bufgets_init(hbuf);
    buf = bufgets(line, sizeof line, hbuf);
    bufgets_finalize(hbuf);

    if (buf == NULL) {
	return;
    }

    tailstrip(line);

    if (gui_help(hwin->role)) {
	/* topic heading: descriptive string */
	gchar *p = quoted_help_string(line);

	gtk_text_buffer_insert_with_tags_by_name(textb, &iter,
						 p, -1,
						 "sansbold", NULL);
	free(p);
    } else {
	/* topic heading: plain command word */
	char hword[12];

	sscanf(line + 2, "%11s", hword);
	gtk_text_buffer_insert_with_tags_by_name(textb, &iter,
						 hword, -1,
						 "redtext", NULL);
    }

    if (hwin->role == FUNCS_HELP) {
	gtk_text_buffer_insert(textb, &iter, "\n\n", 2);
    } else {
	gtk_text_buffer_insert(textb, &iter, "\n", 1);
    }

    buf = grab_topic_buffer(hbuf + strlen(line) + 1);
    if (buf == NULL) {
	return;
    }

    insert_text_with_markup(textb, &iter, buf, hwin->role);
    free(buf);

    gtk_text_view_set_buffer(GTK_TEXT_VIEW(hwin->text), textb);
    push_backpage(hwin->text, hwin->active_var);
    maybe_connect_help_signals(hwin, en);
    cursor_to_top(hwin);
    hwin->active_var = hcode;
}

static int get_screen_height (void)
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

void create_text (windata_t *vwin, int hsize, int vsize, 
		  int nlines, gboolean editable)
{
    GtkTextBuffer *tbuf = gretl_text_buf_new();
    GtkWidget *w = gtk_text_view_new_with_buffer(tbuf);

    vwin->text = w;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(w), 4);

    gtk_widget_modify_font(GTK_WIDGET(w), fixed_font);

    if (hsize > 0 || nlines > 0) {
	int px, py;

	get_char_width_and_height(w, &px, &py);

	if (hsize > 0) {
	    hsize *= px;
	    hsize += 48;
	}

	if (nlines > 0) {
	    double v1 = (nlines + 2) * py;
	    int sv = get_screen_height();

	    if (v1 > vsize / 1.2 && v1 < vsize * 1.2 && v1 <= .9 * sv) {
		vsize = v1;
	    }
	}
    }

    gtk_window_set_default_size(GTK_WINDOW(vwin->main), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), editable);
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
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(sw),
				   GTK_POLICY_AUTOMATIC,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(sw),
					GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(sw), w); 
    gtk_widget_show(w);
    gtk_widget_show(sw);
}

