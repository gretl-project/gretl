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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* console.c for gretl */

#include "gretl.h"

/* from library.c */
extern int gui_exec_line (char *line, 
			  LOOPSET *plp, int *plstack, int *plrun, 
			  SESSION *psession, SESSIONBUILD *rebuild,
			  PRN *prn, int exec_code, 
			  const char *myname);
extern const char *grab_last_cmd (void);

GtkItemFactoryEntry console_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/File/Save _As..."), NULL, file_save, SAVE_CONSOLE, NULL },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

static GtkWidget *console_view;
static char errline[MAXLEN] = ""; /* for use with the console */


static void console_exec (void)
{
    PRN *prn;
    int loopstack = 0, looprun = 0;
    gchar *c_line; 
    char execline[MAXLEN];
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    GtkTextMark *mark;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_end_iter(buf, &end);
    start = end;
    gtk_text_iter_set_line_offset(&start, 2);
    c_line = gtk_text_buffer_get_text(buf, &start, &end, FALSE);
    top_n_tail(c_line);

    if (strcmp(c_line, "quit") == 0 || strcmp(c_line, "q") == 0) {
	gtk_widget_destroy(console_view->parent->parent->parent);
	g_free(c_line);
	return;
    }

    if (bufopen(&prn)) {
	g_free(c_line);
	return;
    }

    strncpy(execline, c_line, MAXLEN - 1);
    g_free(c_line);
    gui_exec_line(execline, NULL, &loopstack, &looprun, NULL, NULL, 
		  prn, CONSOLE_EXEC, NULL);

    /* put results into console window */
    gtk_text_buffer_get_end_iter(buf, &start);
    gtk_text_buffer_insert(buf, &start, "\n", 1);
    gtk_text_buffer_insert(buf, &start, prn->buf, strlen(prn->buf));

    gretl_print_destroy(prn);

    gtk_text_buffer_insert_with_tags_by_name(buf, &start, 
					     "\n? ", 3,
					     "redtext", NULL);

    /* go to end */
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW(console_view),
					mark);
}

void console (void)
{
    PRN *prn;
    char fname[MAXLEN];
    windata_t *vwin;
    GtkTextBuffer *buf;
    GtkTextIter end;

    if (console_view != NULL) {
	gdk_window_show(console_view->parent->window);
	gdk_window_raise(console_view->parent->window);
	return;
    }

    if (!user_fopen("console_tmp", fname, &prn)) return;

    pprintf(prn, _("gretl console: type 'help' for a list of commands\n? "));
    gretl_print_destroy(prn);
    vwin = view_file(fname, 1, 0, 78, 400, CONSOLE, console_items);
    console_view = vwin->w;
    g_signal_connect(G_OBJECT(console_view), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &console_view);

    /* go to end of last line of text */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_end_iter(buf, &end);
    gtk_text_buffer_place_cursor(buf, &end);

    gtk_widget_grab_focus (console_view);
}

#define IS_BACKKEY(k) (k == GDK_BackSpace || k == GDK_Left)

/* ........................................................... */

gboolean console_handler (GtkWidget *w, GdkEventKey *key, gpointer d)
{
    static int lastkey;
    int savekey = key->keyval;
    gint last_line, curr_line, line_pos;
    GdkModifierType mods;
    GtkTextBuffer *buf;
    GtkTextIter iter;

    /* where are we? */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));
    curr_line = gtk_text_iter_get_line(&iter);    
    line_pos = gtk_text_iter_get_line_index(&iter);

    /* null action if not at prompt */
    last_line = gtk_text_buffer_get_line_count(buf) - 1;
    if (curr_line < last_line || (IS_BACKKEY(key->keyval) && line_pos < 3)) {
	return TRUE;
    }

    /* make return key execute the command */
    if (key->keyval == GDK_Return) {
	console_exec();
	key->keyval = GDK_End;
	lastkey = savekey;
	return FALSE;
    }

    /* make up-arrow recall last command entered */
    if (key->keyval == GDK_Up) {
#ifdef CMD_DEBUG
	fprintf(stderr, "errline: '%s'\n", errline);
#endif
	key->keyval = GDK_VoidSymbol;
	if (lastkey != GDK_Up) {
	    if (*errline != 0) {
		gtk_text_buffer_insert(buf, &iter, errline, strlen(errline));
	    } else {
		const char *lastcmd = grab_last_cmd();

		if (lastcmd != NULL) {
		    gtk_text_buffer_insert(buf, &iter, lastcmd, 
					   strlen(lastcmd) - 1);
		}
	    }
	}
	lastkey = savekey;
	return TRUE;
    }

    /* down-arrow clears line */
    if (lastkey == GDK_Up && key->keyval == GDK_Down) {
	GtkTextIter start, end;

	start = end = iter;
	gtk_text_iter_set_line_index(&start, 2);
	gtk_text_iter_forward_to_line_end(&end);
	gtk_text_buffer_delete(buf, &start, &end);
	lastkey = savekey;
	return TRUE;
    } 

    /* Ctrl-A: go to start of typing area */
    gdk_window_get_pointer(console_view->window, NULL, NULL, &mods);
    if (mods & GDK_CONTROL_MASK && 
	gdk_keyval_to_upper(key->keyval) == GDK_A) {
	GtkTextIter where = iter;

	gtk_text_iter_set_line_index(&where, 2);	
	gtk_text_buffer_place_cursor(buf, &where);
	lastkey = savekey;
	return TRUE;
    }

    lastkey = savekey;
    return FALSE;
}

void set_errline (int err, const char *line)
{
    if (err) strcpy(errline, line);
    else *errline = 0;
}
