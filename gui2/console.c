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
#include "console.h"

/* from library.c */
extern int gui_exec_line (char *line, 
			  LOOPSET *plp, int *plstack, int *plrun, 
			  PRN *prn, int exec_code, 
			  const char *myname);

GtkItemFactoryEntry console_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL }, 
    { N_("/File/Save _As..."), NULL, file_save, SAVE_CONSOLE, NULL, GNULL },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL, GNULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

static GtkWidget *console_view;

#define DEFAULT_HLINES 32

static char **cmd_history;
static int hl, hlmax, hlines;

static int gretl_console_init (void)
{
    int i;
    char *hstr;

    hlines = 0;

    hstr = getenv("GRETL_HISTORY_LINES");
    if (hstr != NULL) hlines = atoi(hstr);

    if (hlines <= 2 || hlines > 128) 
	hlines = DEFAULT_HLINES;

    cmd_history = mymalloc(hlines * sizeof *cmd_history);
    if (cmd_history == NULL) return E_ALLOC;

    for (i=0; i<hlines; i++) 
	cmd_history[i] = NULL;

    hlmax = hl = 0;

    return 0;
}

static void gretl_console_free (GtkWidget *w, gpointer p)
{
    int i;

    if (cmd_history != NULL) {
	for (i=0; i<hlines; i++) free(cmd_history[i]);
	free(cmd_history);
	cmd_history = NULL;
    }

    hlines = hlmax = hl = 0;
}

static int push_history_line (const char *line)
{
    int i;

    /* drop last entry */
    free(cmd_history[hlines-1]);

    /* push the other lines down one place */
    for (i=hlines-1; i>0; i--) {
	cmd_history[i] = cmd_history[i-1];
    }

    /* add the new line */
    cmd_history[0] = g_strdup(line);
    
    if (hlmax < hlines) hlmax++;

    hl = 0;

    return 0;
}

static void beep (void)
{
    putchar('\a');
    fflush(stdout);
}

static char *pop_history_line (int keyval)
{
    static int blank;

    if (keyval == GDK_Up) {
	if (hl < 0) hl = (hlmax > 1 && !blank)? 1 : 0;
	if (hl == hlmax) {
	    beep();
	    return NULL;
	}
	blank = 0;
	return cmd_history[hl++];
    }

    if (keyval == GDK_Down) {
	if (hl == hlmax) hl = hlmax - 2;
	if (hl < 0) {
	    blank = 1;
	    beep();
	    return NULL;
	}
	blank = 0;
	return cmd_history[hl--];
    }

    return NULL;
}

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

    *execline = 0;
    strncat(execline, c_line, MAXLEN - 1);
    g_free(c_line);
    push_history_line(execline);    
    gui_exec_line(execline, NULL, &loopstack, &looprun,  prn, 
		  CONSOLE_EXEC, NULL);

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

void show_gretl_console (void)
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
    if (gretl_console_init()) return;

    pprintf(prn, _("gretl console: type 'help' for a list of commands\n? "));
    gretl_print_destroy(prn);
    vwin = view_file(fname, 1, 0, 78, 400, CONSOLE, console_items);
    console_view = vwin->w;

    g_signal_connect(G_OBJECT(console_view), "destroy",
		     G_CALLBACK(gretl_console_free),
		     NULL);
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

gint console_key_handler (GtkWidget *w, GdkEventKey *key, gpointer d)
{
    gint last_line, curr_line, line_pos;
    GdkModifierType mods;
    GtkTextBuffer *buf;
    GtkTextIter iter;

    /* where are we? */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));
    curr_line = gtk_text_iter_get_line(&iter);    
    line_pos = gtk_text_iter_get_line_index(&iter);
    last_line = gtk_text_buffer_get_line_count(buf) - 1;

    /* if at start of command line, backspacing does nothing */
    if (IS_BACKKEY(key->keyval) && line_pos < 3) return TRUE;

    /* if not on prompt line, return to (the end of) it */
    if (curr_line < last_line) {
	GtkTextIter end;

	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_place_cursor(buf, &end);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(console_view), TRUE);	
	return FALSE;
    }

    /* make return key execute the command */
    if (key->keyval == GDK_Return) {
	console_exec();
	key->keyval = GDK_End;
	return FALSE;
    }

    /* up and down arrows: command history */
    if (key->keyval == GDK_Up || key->keyval == GDK_Down) {
	char *histline;

	histline = pop_history_line(key->keyval);

	if (histline != NULL || key->keyval == GDK_Down) {
	    GtkTextIter start, end;

	    start = end = iter;
	    gtk_text_iter_set_line_index(&start, 2);
	    gtk_text_iter_forward_to_line_end(&end);
	    gtk_text_buffer_delete(buf, &start, &end);
	    if (histline != NULL) {
		gtk_text_buffer_insert(buf, &start, histline, strlen(histline));
	    }
	}
	
	return TRUE;
    }

    /* Ctrl-A: go to start of typing area */
    gdk_window_get_pointer(console_view->window, NULL, NULL, &mods);
    if (mods & GDK_CONTROL_MASK && 
	gdk_keyval_to_upper(key->keyval) == GDK_A) {
	GtkTextIter where = iter;

	gtk_text_iter_set_line_index(&where, 2);	
	gtk_text_buffer_place_cursor(buf, &where);
	return TRUE;
    }

    /* tab completion for gretl commands */
    if (key->keyval == GDK_Tab) {
	GtkTextIter start, end;	
	gchar *bit;

	start = end = iter;
	gtk_text_iter_set_line_index(&start, 2);
	gtk_text_iter_forward_to_line_end(&end);
	bit = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

	if (bit != NULL && *bit != 0) {
	    const char *complete = NULL;
	    int i, len = strlen(bit);

	    for (i=0; i<NC; i++) {
		if (strncmp(bit, gretl_commands[i], len) == 0)
		    complete = gretl_commands[i];
	    }

	    g_free(bit);

	    if (complete != NULL) {
		gtk_text_buffer_delete(buf, &start, &end);
		gtk_text_buffer_insert(buf, &start, complete, strlen(complete));
	    }
	}
	return TRUE;
    }

    return FALSE;
}

static gint on_last_line (void)
{
    GtkTextBuffer *buf;
    GtkTextIter iter;
    
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));

    return (gtk_text_iter_get_line(&iter) == 
	    gtk_text_buffer_get_line_count(buf) - 1);
}

gint console_mouse_handler (GtkWidget *w, GdkEventButton *event,
			    gpointer p)
{
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), on_last_line());

    return FALSE;
}

