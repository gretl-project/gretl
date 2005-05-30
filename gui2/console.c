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
#include "menustate.h"

#include "libset.h"
#include "gretl_func.h"

static GtkWidget *console_view;
static PRN *console_prn;
static gchar *cbuf;

#define DEFAULT_HLINES 32

static char **cmd_history;
static int hl, hlmax, hlines;

static LOOPSET *loop;
static int loopstack, looprun;

static int gretl_console_init (void)
{
    char *hstr;
    int i;

    hlines = 0;

    hstr = getenv("GRETL_HISTORY_LINES");
    if (hstr != NULL) {
	hlines = atoi(hstr);
    }

    if (hlines <= 2 || hlines > 128) {
	hlines = DEFAULT_HLINES;
    }

    cmd_history = mymalloc(hlines * sizeof *cmd_history);
    if (cmd_history == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<hlines; i++) {
	cmd_history[i] = NULL;
    }

    hlmax = 0;
    hl = -1;

    cbuf = NULL;

    loop = NULL;
    loopstack = looprun = 0;

    set_gretl_echo(1);

    return 0;
}

static void gretl_console_free (GtkWidget *w, gpointer p)
{
    int i;

    if (cmd_history != NULL) {
	for (i=0; i<hlines; i++) {
	    free(cmd_history[i]);
	}
	free(cmd_history);
	cmd_history = NULL;
    }

    hlines = hlmax = 0;
    hl = -1;

    if (console_prn != NULL) {
	infobox(_("Closing redirected output file"));
	gretl_print_destroy(console_prn);
	console_prn = NULL;
    }

    free(cbuf);
    cbuf = NULL;
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
    
    if (hlmax < hlines) {
	hlmax++;
    }

    hl = -1;

    return 0;
}

static void beep (void)
{
#ifdef G_OS_WIN32
    MessageBeep(MB_ICONEXCLAMATION);
#else
    putchar('\a');
    fflush(stdout);
#endif
}

static char *pop_history_line (int keyval)
{
    static int beeptime;
    char *ret = NULL;

    if (keyval == GDK_Up) {
	if (hl < hlmax) hl++;
	if (hl == hlmax) {
	    hl--;
	    beep();
	} else {
	    ret = cmd_history[hl];
	}
    }

    else if (keyval == GDK_Down) {
	if (hl >= 0) hl--;
	if (hl < 0) {
	    if (beeptime) {
		beep();
	    }
	    beeptime = 1;
	} else {
	    beeptime = 0;
	    ret = cmd_history[hl];
	}
    }

    return ret;
}

#ifdef OLD_GTK

static int last_console_line_len (int len)
{
    int i, c;

    for (i=len; i>0; i--) {
	c = GTK_TEXT_INDEX(GTK_TEXT(console_view), i - 1);
	if (c == '\n') {
	    break; 
	}
    }
    return len - i - 2;
}

static void console_scroll_to_end (void)
{
    int len = gtk_text_get_length(GTK_TEXT(console_view));

    gtk_editable_set_position(GTK_EDITABLE(console_view), len);
}

#else

static void console_scroll_to_end (GtkTextBuffer *buf, 
				   GtkTextIter *start)
{
    GtkTextMark *mark;

    gtk_text_buffer_place_cursor(buf, start);
    mark = gtk_text_buffer_create_mark(buf, NULL, start, FALSE);
    gtk_text_view_scroll_mark_onscreen (GTK_TEXT_VIEW(console_view),
					mark);
}

#endif

enum {
    SAMPLE_RECORD,
    SAMPLE_CHECK
};

static int console_sample_handler (const DATAINFO *pdinfo, int code)
{
    static int pd, t1, t2, ts;
    static double sd0;

    int ret = 0;

    if (code == SAMPLE_RECORD) {
	pd = pdinfo->pd;
	t1 = pdinfo->t1;
	t2 = pdinfo->t2;
	ts = pdinfo->structure;
	sd0 = pdinfo->sd0;
    } else if (code == SAMPLE_CHECK) {
	if (pdinfo->pd != pd ||
	    pdinfo->t1 != t1 ||
	    pdinfo->t2 != t2 ||
	    pdinfo->structure != ts ||
	    pdinfo->sd0 != sd0) {
	    ret = 1;
	}
    }

    return ret;
}

void console_record_sample (const DATAINFO *pdinfo)
{
    console_sample_handler(pdinfo, SAMPLE_RECORD);
}

static int sample_changed (const DATAINFO *pdinfo)
{
    return console_sample_handler(pdinfo, SAMPLE_CHECK);
}

static int console_function_exec (char *execline)
{
    char *gotline;
    int err = 0;

    while (!looprun) {
	gotline = gretl_function_get_line(execline, MAXLINE, &Z, datainfo);
	if (gotline == NULL || *gotline == '\0') {
	    break;
	}
	err = gui_exec_line(execline, &loop, &loopstack, &looprun, console_prn, 
			    SCRIPT_EXEC, NULL);
    }

    return err;
}

static void console_exec (void)
{
    static int redirected;
    int oldv = datainfo->v;
    char execline[MAXLINE];
#ifndef OLD_GTK
    GtkTextBuffer *buf;
    GtkTextIter start, end;
#else
    int len;
    extern GdkColor red;
#endif

#ifndef OLD_GTK
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_end_iter(buf, &end);
    start = end;
    gtk_text_iter_set_line_offset(&start, 2);
#else
    len = gtk_text_get_length(GTK_TEXT(console_view));
#endif

    top_n_tail(cbuf);

    if (strcmp(cbuf, "quit") == 0 || 
	strcmp(cbuf, "q") == 0 ||
	strcmp(cbuf, "exit") == 0 ||
	strcmp(cbuf, "x") == 0) {
	gtk_widget_destroy(console_view->parent->parent->parent);
	g_free(cbuf);
	cbuf = NULL;
	return;
    }

    if (console_prn == NULL && bufopen(&console_prn)) {
	g_free(cbuf);
	cbuf = NULL;
	return;
    }

    *execline = 0;
    strncat(execline, cbuf, MAXLINE - 1);
    g_free(cbuf);
    cbuf = NULL;

    console_record_sample(datainfo);

    push_history_line(execline); 

    /* actually execute the command line */
    gui_exec_line(execline, &loop, &loopstack, &looprun, console_prn, 
		  CONSOLE_EXEC, NULL);

    /* the control structure below is rather weird and needs more
       testing/thought.  The issue is the nesting of loops inside
       functions or vice versa */

    if (looprun) { 
    fn_run_loop:
	loop_exec(loop, execline, &Z, &datainfo, models, console_prn);
	gretl_loop_destroy(loop);
	loop = NULL;
	looprun = 0;
	if (gretl_executing_function_or_macro()) {
	    console_function_exec(execline);
	}
    } else if (gretl_executing_function_or_macro()) {
	console_function_exec(execline);
	if (looprun) {
	    /* the function we are exec'ing includes a loop */
	    goto fn_run_loop;
	}
    } 

    redirected = printing_is_redirected(console_prn);

#ifndef OLD_GTK
    gtk_text_buffer_get_end_iter(buf, &start);
    gtk_text_buffer_insert(buf, &start, "\n", 1);
#else
    gtk_text_freeze(GTK_TEXT(console_view));
    gtk_editable_insert_text(GTK_EDITABLE(console_view), 
			     "\n", 1, &len);
#endif

    if (!redirected) {
	const char *cbuf = gretl_print_get_buffer(console_prn);

        /* put results into console window */
#ifndef OLD_GTK
	if (!g_utf8_validate(cbuf, -1, NULL)) {
	    fprintf(stderr, "text did not validate as utf8:\n'%s'\n", 
		    cbuf);
	} else {
	    gtk_text_buffer_insert(buf, &start, cbuf, strlen(cbuf));
	}
#else
	gtk_text_insert(GTK_TEXT(console_view), fixed_font, 
			NULL, NULL, cbuf, strlen(cbuf));
#endif
	gretl_print_destroy(console_prn);
	console_prn = NULL;
    } else {
	gretl_print_reset_buffer(console_prn);
    }

#ifndef OLD_GTK
    gtk_text_buffer_insert_with_tags_by_name(buf, &start, 
					     (loopstack)? "\n> " : "\n? ", 3,
					     "redtext", NULL);
#else
    gtk_text_insert(GTK_TEXT(console_view), fixed_font,
		    &red, NULL, (loopstack)? "\n> " : "\n? ", 3);
    gtk_text_thaw(GTK_TEXT(console_view));
#endif

    /* scroll to end of buffer */
#ifndef OLD_GTK
    console_scroll_to_end(buf, &start);
#else
    console_scroll_to_end();
#endif

    /* update variable listing in main window if needed */
    if (datainfo->v != oldv || !strncmp(execline, "rename", 6)) {
	populate_varlist();
    }

    /* update sample info and options if needed */
    if (sample_changed(datainfo)) {
	set_sample_label(datainfo);
    }
}

void show_gretl_console (void)
{
    PRN *prn;
    char fname[MAXLEN];
    windata_t *vwin;
#ifndef OLD_GTK
    GtkTextBuffer *buf;
    GtkTextIter end;
#else
    size_t introlen;
#endif
    const char *intro = 
	N_("gretl console: type 'help' for a list of commands\n? ");

    if (console_view != NULL) {
	gdk_window_show(console_view->parent->window);
	gdk_window_raise(console_view->parent->window);
	return;
    }

    if (!user_fopen("console_tmp", fname, &prn)) {
	return;
    }

    if (gretl_console_init()) {
	return;
    }

    pputs(prn, _(intro));
    gretl_print_destroy(prn);
#ifdef OLD_GTK
    introlen = strlen(_(intro));
#endif

    vwin = view_file(fname, 1, 0, 78, 400, CONSOLE);
    console_view = vwin->w;

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(console_view), "destroy",
		     G_CALLBACK(gretl_console_free),
		     NULL);
    g_signal_connect(G_OBJECT(console_view), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &console_view);
#else
    gtk_signal_connect(GTK_OBJECT(console_view), "destroy",
		       GTK_SIGNAL_FUNC(gretl_console_free),
		       NULL);
    gtk_signal_connect(GTK_OBJECT(console_view), "destroy",
		       GTK_SIGNAL_FUNC(gtk_widget_destroyed),
		       &console_view);
#endif

    /* go to end of last line of text */
#ifndef OLD_GTK
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_end_iter(buf, &end);
    gtk_text_buffer_place_cursor(buf, &end);
#else
    gtk_text_set_point(GTK_TEXT(console_view), introlen); 
    GTK_TEXT(console_view)->cursor_pos_x = 2 * gdk_char_width(fixed_font, 'X');
    gtk_editable_set_position(GTK_EDITABLE(console_view), introlen);
    gtk_editable_changed(GTK_EDITABLE(console_view));
#endif

    gtk_widget_grab_focus(console_view);
}

#define IS_BACKKEY(k) (k == GDK_BackSpace || k == GDK_Left)

static int bslash_cont (const gchar *line)
{
    int bslash = ends_with_backslash(line);

    if (cbuf == NULL) {
	cbuf = g_strdup(line);
    } else {
	cbuf = g_realloc(cbuf, strlen(cbuf) + 1 + strlen(line));
	if (cbuf != NULL) {
	    strcat(cbuf, line);
	}
    }

    if (cbuf == NULL) {
	return 0;
    }

    if (bslash) {
	char *p = strrchr(cbuf, '\\');

	if (p - cbuf > 0 && !isspace(*(p - 1))) {
	    *p = ' ';
	} else {
	    *p = 0;
	}
    }

    return bslash;
}

#ifndef OLD_GTK

gint console_key_handler (GtkWidget *w, GdkEventKey *key, gpointer d)
{
    gint last_line, curr_line, line_pos;
    GdkModifierType mods;
    GtkTextBuffer *buf;
    GtkTextIter iter;

 start_again:

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
	goto start_again;
    }

    /* make return key execute the command, unless backslash-
       continuation is happening */
    if (key->keyval == GDK_Return) {
	GtkTextIter start, end;	
	gchar *bit;
	int cont = 0;

	start = end = iter;
	gtk_text_iter_set_line_index(&start, 2);
	gtk_text_iter_forward_to_line_end(&end);
	bit = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

	if (bit != NULL) {
	    cont = bslash_cont(bit);
	    g_free(bit);
	}

	if (cont) {
	    gtk_text_buffer_insert_with_tags_by_name(buf, &end, 
						     "\n> ", 3,
						     "redtext", NULL);
	    console_scroll_to_end(buf, &end);
	} else {
	    console_exec();
	}

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

	if (bit != NULL) {
	    if (*bit != 0) {
		const char *complete = gretl_command_complete(bit);

		if (complete != NULL) {
		    gtk_text_buffer_delete(buf, &start, &end);
		    gtk_text_buffer_insert(buf, &start, complete, 
				       strlen(complete));
		}
	    }
	    g_free(bit);
	}
	return TRUE;
    }

    return FALSE;
}

#else /* alternate gtk versions */

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

    /* make return key execute the command, unless backslash-
       continuation is happening */
    if (key->keyval == GDK_Return) {
	gchar *bit;
	int start, cont = 0;
	extern GdkColor red;

	start = len - last_console_line_len(len);
	bit = gtk_editable_get_chars(GTK_EDITABLE(console_view), 
				     start, -1);
	if (bit != NULL) {
	    cont = bslash_cont(bit);
	    g_free(bit);
	}

	if (cont) {
	    gtk_text_insert(GTK_TEXT(console_view), fixed_font,
			    &red, NULL, "\n> ", 3);
	    console_scroll_to_end();
	} else {
	    console_exec();
	}

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
	    const char *complete = gretl_command_complete(bit);

	    g_free(bit);
	    if (complete != NULL) {
		replace_command_line(cw, complete);
	    }
	}
	gtk_signal_emit_stop_by_name(GTK_OBJECT(w), "key-press-event");
	return TRUE;
    }    

    return FALSE;
}

#endif

static gint on_last_line (void)
{
#ifndef OLD_GTK
    GtkTextBuffer *buf;
    GtkTextIter iter;
    
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(console_view));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));

    return (gtk_text_iter_get_line(&iter) == 
	    gtk_text_buffer_get_line_count(buf) - 1);
#else
    int len = gtk_text_get_length(GTK_TEXT(console_view));
    int currpos = GTK_EDITABLE(console_view)->current_pos;

    return (currpos >= len - last_console_line_len(len));
#endif
}

gint console_mouse_handler (GtkWidget *w, GdkEventButton *event,
			    gpointer p)
{
    int lline = on_last_line();

#ifndef OLD_GTK
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(w), lline);
#else
    gtk_text_set_editable(GTK_TEXT(w), lline);
#endif

    return FALSE;
}


