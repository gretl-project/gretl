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
    { N_("/_File"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/File/Save _As..."), NULL, file_save, SAVE_CONSOLE, NULL },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { NULL, NULL, NULL, 0, NULL }
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

static int last_console_line_len (int len)
{
    int i, c;

    for (i=len; i>0; i--) {
	c = GTK_TEXT_INDEX(GTK_TEXT(console_view), i - 1);
	if (c == '\n') break; 
    }
    return len - i - 2;
}

static void console_exec (void)
{
    PRN *prn;
    int len, loopstack = 0, looprun = 0;
    gchar *c_line; 
    char execline[MAXLEN];
    extern GdkColor red;

    len = gtk_text_get_length (GTK_TEXT(console_view));
    c_line = gtk_editable_get_chars(GTK_EDITABLE(console_view), 
				    len - last_console_line_len(len), len);
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
    push_history_line(execline);
    gui_exec_line(execline, NULL, &loopstack, &looprun, 
		  prn, CONSOLE_EXEC, NULL);

    /* put results into console window */
    gtk_text_freeze(GTK_TEXT(console_view));
    gtk_editable_insert_text(GTK_EDITABLE(console_view), 
			     "\n", 1, &len);
    gtk_text_insert(GTK_TEXT(console_view), fixed_font, 
		    NULL, NULL, prn->buf, strlen(prn->buf));
    gretl_print_destroy(prn);

    gtk_text_insert(GTK_TEXT(console_view), fixed_font,
		    &red, NULL, "\n? ", 3);
    gtk_text_thaw(GTK_TEXT(console_view));

    /* go to end */
    len = gtk_text_get_length(GTK_TEXT(console_view));
    gtk_editable_set_position(GTK_EDITABLE(console_view), len);
}

void show_gretl_console (void)
{
    PRN *prn;
    char fname[MAXLEN];
    windata_t *vwin;
    const char *intro = 
	N_("gretl console: type 'help' for a list of commands\n? ");
    size_t introlen;

    if (console_view != NULL) {
	gdk_window_show(console_view->parent->window);
	gdk_window_raise(console_view->parent->window);
	return;
    }

    if (!user_fopen("console_tmp", fname, &prn)) return;
    if (gretl_console_init()) return;

    pputs(prn, _(intro));
    gretl_print_destroy(prn);
    introlen = strlen(_(intro));

    vwin = view_file(fname, 1, 0, 78, 400, CONSOLE, console_items);
    console_view = vwin->w;

    gtk_signal_connect(GTK_OBJECT(console_view), "destroy",
		       GTK_SIGNAL_FUNC(gretl_console_free),
		       NULL);
    gtk_signal_connect(GTK_OBJECT(console_view), "destroy",
		       GTK_SIGNAL_FUNC(gtk_widget_destroyed),
		       &console_view);
    
    gtk_text_set_point(GTK_TEXT(console_view), introlen); 

    GTK_TEXT(console_view)->cursor_pos_x = 2 * gdk_char_width(fixed_font, 'X');
    gtk_editable_set_position(GTK_EDITABLE(console_view), introlen);
    gtk_editable_changed(GTK_EDITABLE(console_view));

    gtk_widget_grab_focus (console_view);
}

#define IS_BACKKEY(k) (k == GDK_BackSpace || k == GDK_Left)

static void replace_command_line (int cw, const char *repl)
{
    int len = gtk_text_get_length(GTK_TEXT(console_view));
    int adjust;

    adjust = GTK_TEXT(console_view)->cursor_pos_x / cw - 2;
    gtk_editable_delete_text(GTK_EDITABLE(console_view), 
			     len - adjust, len);
    len = gtk_text_get_length(GTK_TEXT(console_view));
    if (repl != NULL) {
	gtk_editable_insert_text(GTK_EDITABLE(console_view), 
				 repl, strlen(repl),
				 &len);
    }
}

gint console_key_handler (GtkWidget *w, GdkEventKey *key, gpointer d)
{
    int cw = gdk_char_width(fixed_font, 'X'); 
    int len = gtk_text_get_length(GTK_TEXT(console_view));
    int currpos, linelen, xpos;
    GdkModifierType mods;

    /* where are we? */
    currpos = GTK_EDITABLE(console_view)->current_pos;
    xpos = GTK_TEXT(console_view)->cursor_pos_x / cw; 
    linelen = last_console_line_len(len);

    /* if at start of command line, backspacing does nothing */
    if (IS_BACKKEY(key->keyval) && xpos < 3) {
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }

    /* if not on prompt line, return to (the end of) it */
    if (currpos < len - linelen) {
	gtk_text_set_point(GTK_TEXT(console_view), len);
	gtk_editable_set_position(GTK_EDITABLE(console_view), len);
	gtk_text_set_editable(GTK_TEXT(console_view), TRUE);
	return TRUE;
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
	    replace_command_line(cw, histline);
	}
	key->keyval = GDK_VoidSymbol;
	return TRUE;
    }

    /* response to Ctrl-A: go to start of typing area */
    gdk_window_get_pointer(console_view->window, NULL, NULL, &mods);
    if (mods & GDK_CONTROL_MASK && 
	gdk_keyval_to_upper(key->keyval) == GDK_A) {
	gtk_editable_set_position(GTK_EDITABLE(console_view), 
				  len - last_console_line_len(len));	
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }

    /* tab completion for gretl commands */
    if (key->keyval == GDK_Tab) {
	gchar *bit;
	int start;

	start = len - last_console_line_len(len);
	bit = gtk_editable_get_chars(GTK_EDITABLE(console_view), 
				     start, -1);

	if (bit != NULL && *bit != 0) {
	    const char *complete = NULL;
	    int i, len = strlen(bit);

	    for (i=0; i<NC; i++) {
		if (strncmp(bit, gretl_commands[i], len) == 0)
		    complete = gretl_commands[i];
	    }

	    g_free(bit);

	    if (complete != NULL) replace_command_line(cw, complete);
	}
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }    

    return FALSE;
}

static int on_last_line (void)
{
    int len = gtk_text_get_length(GTK_TEXT(console_view));
    int currpos = GTK_EDITABLE(console_view)->current_pos;

    return (currpos >= len - last_console_line_len(len));
}

gint console_mouse_handler (GtkWidget *w, GdkEventButton *event,
			    gpointer p)
{
    gtk_text_set_editable(GTK_TEXT(w), on_last_line());
    return FALSE;
}

