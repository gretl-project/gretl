/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "gretl.h"
#include "textbuf.h"

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
# include <gtksourceview/gtksourcelanguage.h>
# include <gtksourceview/gtksourcelanguagesmanager.h>
#endif

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

/* .................................................................. */

void text_set_cursor (GtkWidget *w, GdkCursorType cspec)
{
    GdkWindow *win = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
                                              GTK_TEXT_WINDOW_TEXT);

    if (cspec == 0) {
	gdk_window_set_cursor(win, NULL);
    } else {
	GdkCursor *cursor = gdk_cursor_new(cspec);

	gdk_window_set_cursor(win, cursor);
	gdk_cursor_destroy(cursor);
    } 
}

/* ........................................................... */

void cursor_to_top (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w)); 
    GtkTextIter start;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				 mark, 0.0, FALSE, 0, 0);
}

/* .................................................................. */

gint get_char_width (GtkWidget *widget)
{
    PangoLayout *pl;
    PangoContext *pc;
    GtkRcStyle *style;
    int width;

    pc = gtk_widget_get_pango_context(widget);
    style = gtk_widget_get_modifier_style(widget);
    pango_context_set_font_description(pc, style->font_desc);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, "X", 1);
    pango_layout_get_pixel_size(pl, &width, NULL);

    g_object_unref(G_OBJECT(pl));

    return width;
}

/* .................................................................. */

gchar *textview_get_text (GtkTextView *view)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    tbuf = gtk_text_view_get_buffer(view);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);

    return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
}

/* .................................................................. */

void text_paste (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old;
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gchar *undo_buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));

    old = g_object_get_data(G_OBJECT(vwin->w), "undo");
    g_free(old);

    g_object_set_data(G_OBJECT(vwin->w), "undo", undo_buf);

    gtk_text_buffer_paste_clipboard(buf, gtk_clipboard_get(GDK_NONE),
				    NULL, TRUE);
}

/* .................................................................. */

void text_undo (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old = NULL;

#ifdef USE_GTKSOURCEVIEW
    if (vwin->sbuf != NULL) {
	if (gtk_source_buffer_can_undo(vwin->sbuf)) {
	    gtk_source_buffer_undo(vwin->sbuf);
	} else {
	    errbox(_("No undo information available"));
	}
	return;
    }
#endif
    
    old = g_object_steal_data(G_OBJECT(vwin->w), "undo");

    if (old == NULL) {
	errbox(_("No undo information available"));
    } else {
	GtkTextBuffer *buf;
	GtkTextIter start, end;
	GtkTextMark *ins;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	ins = gtk_text_buffer_get_insert(buf);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_delete(buf, &start, &end);

	gtk_text_buffer_insert(buf, &start, old, strlen(old));
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				     ins, 0.0, TRUE, 0.1, 0.0);
	g_free(old);
    }
}

#ifdef USE_GTKSOURCEVIEW

static int 
gtk_source_buffer_load_file (GtkSourceBuffer *sbuf, 
			     const char *fname)
{
    FILE *fp;
    GtkTextIter iter;    
    char readbuf[MAXSTR], *chunk = NULL;

    fp = fopen(fname, "rb");
    if (fp == NULL) return 1;

    gtk_source_buffer_begin_not_undoable_action (sbuf);

    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), "", 0);
    gtk_text_buffer_get_iter_at_offset(GTK_TEXT_BUFFER(sbuf), &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf, fp)) {
	int len;

# ifdef ENABLE_NLS
	if (!g_utf8_validate(readbuf, -1, NULL)) {
	    chunk = my_locale_to_utf8(readbuf);
	    if (chunk == NULL) {
		continue;
	    }
	} else {
	    chunk = readbuf;
	}
# else
	chunk = readbuf;
# endif /* ENABLE_NLS */

	/* check that this works */
	len = strlen(chunk);
	if (chunk[len - 2] == '\r') {
	    chunk[len - 2] = '\n';
	    chunk[len - 1] = '\0';
	}

	gtk_text_buffer_insert(GTK_TEXT_BUFFER(sbuf), &iter, chunk, -1);
	memset(readbuf, 0, sizeof readbuf);
	if (chunk != readbuf) {
	    g_free(chunk);
	    chunk = NULL;
	}
    }

    fclose(fp);
	
    gtk_source_buffer_end_not_undoable_action(sbuf);

    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sbuf), FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(sbuf), &iter);
    gtk_text_buffer_place_cursor(GTK_TEXT_BUFFER(sbuf), &iter);

    return 0;
}

void source_buffer_insert_file (GtkSourceBuffer *sbuf, const char *filename,
				int role)
{
    GtkSourceLanguagesManager *manager;    
    GtkSourceLanguage *language = NULL;
		
    manager = g_object_get_data(G_OBJECT (sbuf), "languages-manager");

    if (role == GR_PLOT) {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gnuplot");
    } else {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gretlsession");
    }

    if (language == NULL) {
	g_object_set(G_OBJECT(sbuf), "highlight", FALSE, NULL);
    } else {
	g_object_set(G_OBJECT(sbuf), "highlight", TRUE, NULL);
	gtk_source_buffer_set_language(sbuf, language);
    }

    gtk_source_buffer_load_file(sbuf, filename);
}

void create_source (windata_t *vwin, GtkSourceBuffer **buf, 
		    int hsize, int vsize, gboolean editable)
{
    GtkSourceLanguagesManager *lm;
    GtkSourceBuffer *sbuf;
    GtkSourceTagStyle *tagstyle;
    GdkColor blue;

    blue.green = blue.red = 0;
    blue.blue = 65535.0;

    lm = gtk_source_languages_manager_new ();
    tagstyle = gtk_source_tag_style_new ();
    
    sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(NULL));
    g_object_ref (lm);
    g_object_set_data_full (G_OBJECT (sbuf), "languages-manager",
			    lm, (GDestroyNotify) g_object_unref); 
    g_object_unref (lm); 

    tagstyle->mask = GTK_SOURCE_TAG_STYLE_USE_FOREGROUND;
    tagstyle->foreground = blue;
    g_object_set_data_full (G_OBJECT (sbuf), "tag-style",
			    tagstyle, 
			    (GDestroyNotify) gtk_source_tag_style_free); 
    gtk_source_buffer_set_bracket_match_style(sbuf, tagstyle);
    gtk_source_buffer_set_check_brackets(sbuf, TRUE);

    vwin->w = gtk_source_view_new_with_buffer(sbuf);
    *buf = sbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->w), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->w), fixed_font);

    hsize *= get_char_width(vwin->w);
    hsize += 48;

    gtk_window_set_default_size (GTK_WINDOW(vwin->dialog), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);
}

#endif /* USE_GTKSOURCEVIEW */

static gchar *my_utf_string (char *t)
{
    static gchar *s = NULL;
    GError *error = NULL;
    gsize r_bytes, w_bytes;
    unsigned char *c;
    const char *fc;
    const char *smb;
    gchar *from_codeset = NULL;
    
    if (t == NULL || *t == '\0') return t;

    if (g_utf8_validate(t, -1, NULL)) return t;   
    
    /* so we got a non-UTF-8 */

    smb = getenv("SMB_CODESET");
    if (smb != NULL && *smb != '\0') {
	from_codeset = g_strdup(smb);
    } else {
    	g_get_charset(&fc);
    	if (fc) {
	    from_codeset = g_strdup(fc);
    	} else {
	    from_codeset = g_strdup("ISO-8859-1");
	}
    }
    
    if (!strcmp(from_codeset, "ISO-")) {
	g_free(from_codeset);
	from_codeset = g_strdup("ISO-8859-1");
    }  
  
    if (s) g_free(s);

    for (c = (unsigned char *)t; *c != 0; c++) {
	if (*c < 32 && *c != '\n') {
	    *c = ' ';
	}
    }

    s = g_convert(t, strlen(t), "UTF-8", from_codeset, &r_bytes, &w_bytes,
		  &error);

    g_free(from_codeset);

    if (s == NULL) {
	s = g_strdup(t);
	for (c = s; *c != 0; c++) if (*c > 128) *c = '?';
    }

    if (error) {
        printf("DBG: %s. Codeset for system is: %s\n",
	       error->message,from_codeset);
        printf("DBG: You should set the environment variable "
	       "SMB_CODESET to ISO-8859-1\n");
	g_error_free(error);
    }

    return s;
}

void text_buffer_insert_colorized_buffer (GtkTextBuffer *tbuf, PRN *prn)
{
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    char readbuf[MAXSTR];

    thiscolor = PLAIN_TEXT;
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    bufgets(NULL, 0, prn->buf);

    while (bufgets(readbuf, sizeof readbuf, prn->buf)) {

	if (ends_with_backslash(readbuf)) {
	    nextcolor = BLUE_TEXT;
	} else {
	    nextcolor = PLAIN_TEXT;
	}

	if (*readbuf == '#' || *readbuf == '?') {
	    thiscolor = BLUE_TEXT;
	} 

	if (thiscolor == BLUE_TEXT) {
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      readbuf, -1,
						      "bluetext", NULL);
	} else {
	    gtk_text_buffer_insert(tbuf, &iter, readbuf, -1);
	}

	/* bufgets strips newlines */
	gtk_text_buffer_insert(tbuf, &iter, "\n", 1);

	thiscolor = nextcolor;
    }
}

void text_buffer_insert_file (GtkTextBuffer *tbuf, const char *fname, 
			      int role)
{
    FILE *fp;
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    char readbuf[MAXSTR], *chunk;
    int started = 0;
    int newhelp = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) return;

    thiscolor = nextcolor = PLAIN_TEXT;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    if (role == GUI_HELP || role == GUI_HELP_ENGLISH) {
	newhelp = new_style_gui_help(fp);
    }

    while (fgets(readbuf, sizeof readbuf, fp)) {
#ifdef ENABLE_NLS
	chunk = my_utf_string(readbuf);
#else
	chunk = readbuf;
#endif

	if (help_role(role) && *chunk == '@') continue;

	nextcolor = PLAIN_TEXT;
	
	if (newhelp && *chunk == '#') { /* help topic line */
	    char *p;

	    if (started) {
		gtk_text_buffer_insert(tbuf, &iter, "\n", 1);
	    } else {
		started = 1;
	    }
	    p = quoted_help_string(chunk);
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      p, -1,
						      "redtext", NULL);
	    gtk_text_buffer_insert(tbuf, &iter, "\n", 1);
	    free(p);
	    continue;
	}

	if (role == SCRIPT_OUT && ends_with_backslash(chunk)) {
	    nextcolor = BLUE_TEXT;
	}

	if (*chunk == '?') {
	    thiscolor = (role == CONSOLE)? RED_TEXT : BLUE_TEXT;
	} 
	else if (*chunk == '#') {
	    if (help_role(role)) {
		*chunk = ' ';
		nextcolor = RED_TEXT;
	    } else {
		thiscolor = BLUE_TEXT;
	    }
	} 

	switch (thiscolor) {
	case PLAIN_TEXT:
	    gtk_text_buffer_insert(tbuf, &iter, chunk, -1);
	    break;
	case BLUE_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      chunk, -1,
						      "bluetext", NULL);
	    break;
	case RED_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      chunk, -1,
						      "redtext", NULL);
	    break;
	}

	thiscolor = nextcolor;
	memset(readbuf, 0, sizeof readbuf);
    }

    fclose(fp);
}

/* ........................................................... */

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

    return table;
}

#ifndef USE_GTKSOURCEVIEW

void correct_line_color (windata_t *vwin)
{
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    gint linelen;
    gchar *txt;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gtk_text_buffer_get_iter_at_mark(buf, &end, 
				     gtk_text_buffer_get_insert(buf));
    linelen = gtk_text_iter_get_chars_in_line(&end);
    start = end;
    gtk_text_iter_backward_chars(&start, linelen);

    txt = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

    if (*txt == '#') {
	gtk_text_buffer_apply_tag_by_name (buf, "bluetext",
					   &start, &end);
    }
    g_free(txt);
}

#endif /* not USE_GTKSOURCEVIEW */

/* ........................................................... */

void create_text (windata_t *vwin, GtkTextBuffer **buf, 
		  int hsize, int vsize, gboolean editable)
{
    static GtkTextTagTable *tags = NULL;
    GtkTextBuffer *tbuf; 

    if (tags == NULL) tags = gretl_tags_new();

    tbuf = gtk_text_buffer_new(tags);
    vwin->w = gtk_text_view_new_with_buffer(tbuf);
    *buf = tbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->w), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->w), fixed_font);

    hsize *= get_char_width(vwin->w);
    hsize += 48;

    gtk_window_set_default_size (GTK_WINDOW(vwin->dialog), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);
}

/* ........................................................... */

void text_table_setup (windata_t *vwin)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
				    GTK_POLICY_AUTOMATIC,
				    GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (sw),
					 GTK_SHADOW_IN);
    gtk_container_add (GTK_CONTAINER(sw), vwin->w); 
    gtk_widget_show(vwin->w);
    gtk_widget_show(sw);
}
