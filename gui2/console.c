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

/* console.c for gretl */

#include "gretl.h"
#include "console.h"
#include "menustate.h"
#include "dlgutils.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include "libset.h"
#include "monte_carlo.h"
#include "gretl_func.h"
#include "cmd_private.h"

static char **cmd_history;
static int hpos, hlines;

static gint console_key_handler (GtkWidget *cview, GdkEventKey *key, 
				 char *cbuf);

static void command_history_init (void)
{
    hlines = hpos = 0;
}

static void command_history_destroy (void)
{
    if (cmd_history != NULL) {
	free_strings_array(cmd_history, hlines);
	cmd_history = NULL;
    }

    hlines = hpos = 0;
}

static ExecState *gretl_console_init (void)
{
    ExecState *s;
    PRN *prn;

    s = mymalloc(sizeof *s);
    if (s == NULL) {
	return NULL;
    }   

    if (bufopen(&prn)) {
	free(s);
	return NULL;
    }

    set_gretl_echo(1);

    gretl_exec_state_init(s, CONSOLE_EXEC, NULL, 
			  get_lib_cmd(), models, prn);

    command_history_init();

    return s;
}

static void gretl_console_free (GtkWidget *cview, ExecState *s)
{
    command_history_destroy();

    gretl_print_destroy(s->prn);
    free(s);

    gtk_main_quit();
}

static void push_history_line (const char *line)
{
    strings_array_add(&cmd_history, &hlines, line);
    hpos = hlines;
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

static const char *fetch_history_line (int keyval)
{
    static int dobeep;
    char *ret = NULL;

    if (keyval == GDK_Up) {
	if (hpos > 0) {
	    /* up is OK */
	    ret = cmd_history[--hpos];
	    dobeep = 0;
	} else {
	    /* can't go up */
	    beep();
	}
    } else if (keyval == GDK_Down) {
	if (hlines == 0) {
	    /* no history yet */
	    hpos = 0;
	    beep();
	} else if (hpos < hlines - 1) {
	    /* down is OK */
	    ret = cmd_history[++hpos];
	    dobeep = 0;
	} else {
	    /* can't go down */
	    hpos = hlines;
	    if (dobeep) {
		beep();
	    } else {
		dobeep = 1;
	    }
	}
    }

    return ret;
}

static void console_scroll_to_end (GtkWidget *cview,
				   GtkTextBuffer *buf, 
				   GtkTextIter *end)
{
    GtkTextMark *mark;

    mark = gtk_text_buffer_create_mark(buf, NULL, end, FALSE);
    gtk_text_view_scroll_mark_onscreen(GTK_TEXT_VIEW(cview), mark);
}

static gint on_last_line (GtkWidget *cview)
{
    GtkTextBuffer *buf;
    GtkTextIter iter;
    
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));

    return (gtk_text_iter_get_line(&iter) == 
	    gtk_text_buffer_get_line_count(buf) - 1);
}

static gint console_mouse_handler (GtkWidget *cview, GdkEventButton *event,
				   gpointer p)
{
    int lline = on_last_line(cview);

    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(cview), lline);

    return FALSE;
}

static gint console_paste_text (GtkWidget *cview, GdkAtom atom)
{
    GtkClipboard *cb = gtk_clipboard_get(atom);
    gchar *src = gtk_clipboard_wait_for_text(cb);

    if (src != NULL) {
	GtkTextBuffer *buf;
	GtkTextIter iter;
	char *p;

	p = strchr(src, '\n');
	if (p != NULL) {
	    /* no newlines allowed! */
	    *p = '\0';
	}

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));
	gtk_text_buffer_get_end_iter(buf, &iter);
	gtk_text_buffer_insert(buf, &iter, src, -1);

	g_free(src);
    }

    return TRUE;
}

static void console_paste_handler (GtkWidget *w, gpointer p)
{
    /* we don't accept pasted material, other than via 
       the X selection */
    return;
}

/* paste from X selection onto the command line */

static gint console_click_handler (GtkWidget *w, 
				   GdkEventButton *event,
				   gpointer p)
{
    GdkModifierType mods = widget_get_pointer_mask(w);

    if (mods & GDK_BUTTON2_MASK) {
	return console_paste_text(w, GDK_SELECTION_PRIMARY);
    }

    return FALSE;
}

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
    } else if (code == SAMPLE_CHECK && pdinfo->v > 0) {
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

int console_sample_changed (const DATAINFO *pdinfo)
{
    return console_sample_handler(pdinfo, SAMPLE_CHECK);
}

static void print_result_to_console (GtkTextBuffer *buf,
				     GtkTextIter *start,
				     ExecState *s)
{
    const char *getbuf = gretl_print_get_buffer(s->prn);

    if (!g_utf8_validate(getbuf, -1, NULL)) {
	gchar *trbuf = my_locale_to_utf8(getbuf);

	fprintf(stderr, "console text did not validate as utf8\n");
	if (trbuf != NULL) {
	    gtk_text_buffer_insert(buf, start, trbuf, -1);
	    g_free(trbuf);
	}
    } else {
	gtk_text_buffer_insert(buf, start, getbuf, -1);
    }
}

static void console_insert_prompt (GtkTextBuffer *buf,
				   GtkTextIter *iter,
				   const char *prompt)
{
    gtk_text_buffer_insert_with_tags_by_name(buf, iter, prompt, -1,
					     "redtext", NULL);  
    gtk_text_buffer_place_cursor(buf, iter);
}

/* callback from Return/Enter key in gretl console */

static void console_exec (GtkWidget *cview, char *cbuf)
{
    ExecState *state;
    char execline[MAXLINE];
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    int coding = 0;
    int err = 0;

    state = g_object_get_data(G_OBJECT(cview), "ExecState");

    /* get into printing position */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));
    gtk_text_buffer_get_end_iter(buf, &end);
    start = end;
    gtk_text_iter_set_line_offset(&start, 2);

    /* 'cbuf' contains the current command line: transcribe
       it and zero it out */
    strcpy(execline, cbuf);
    *cbuf = '\0';

    console_record_sample(datainfo);
    push_history_line(execline);
    state->line = execline;
    state->flags = CONSOLE_EXEC;

    /* actually execute the command line */
    err = gui_exec_line(state, &Z, datainfo);

    while (!err && gretl_execute_loop()) {
	err = gretl_loop_exec(state, &Z, datainfo);
    }

    gtk_text_buffer_get_end_iter(buf, &start);
    gtk_text_buffer_insert(buf, &start, "\n", 1);

    if (!printing_is_redirected(state->prn)) {
	print_result_to_console(buf, &start, state);
    }

    gretl_print_reset_buffer(state->prn);

    if (state->cmd->ci == QUIT) {
	gtk_widget_destroy(gtk_widget_get_toplevel(cview));
	return;
    }

    coding = gretl_compiling_loop() || gretl_compiling_function();

    /* set up prompt for next command and scroll to it */
    console_insert_prompt(buf, &start, (coding)? "> " : "? ");
    console_scroll_to_end(cview, buf, &start);

    /* update variable listing in main window if needed */
    if (check_dataset_is_changed()) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    /* update sample info and options if needed */
    if (console_sample_changed(datainfo)) {
	set_sample_label(datainfo);
    }
}

static GtkWidget *console_main;

int maybe_raise_console (void)
{
    if (console_main != NULL) {
	gtk_window_present(GTK_WINDOW(console_main));
	return 1;
    } else {
	return 0;
    }
}

/* callback from menu/button: launch the console */

void show_gretl_console (void)
{
    char cbuf[MAXLINE];
    windata_t *vwin;
    GtkTextBuffer *buf;
    GtkTextIter iter;
    ExecState *cstate;
    const gchar *intro = 
	N_("gretl console: type 'help' for a list of commands");

    if (maybe_raise_console()) {
	return;
    }

    cstate = gretl_console_init();
    if (cstate == NULL) {
	return;
    }

    *cbuf = '\0';

    vwin = console_window(78, 400);

    console_main = vwin->main;
    g_object_set_data(G_OBJECT(vwin->text), "ExecState", cstate);

    g_signal_connect(G_OBJECT(vwin->text), "paste-clipboard",
		     G_CALLBACK(console_paste_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(console_click_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
		     G_CALLBACK(console_mouse_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
		     G_CALLBACK(console_key_handler), cbuf);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(gretl_console_free), cstate);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &console_main);

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_start_iter(buf, &iter);  

    /* insert intro string and first prompt */
    gtk_text_buffer_insert(buf, &iter, _(intro), -1);
    console_insert_prompt(buf, &iter, "\n? ");

    gtk_widget_grab_focus(vwin->text);

    /* Establish a main loop that exits when the console is destroyed.
       We need the blocking to ensure that the state variables in this
       function remain on the stack while the console is being used.
    */
    gtk_main();
}

/* handle backslash continuation of console command line */

static int 
command_continues (char *targ, const gchar *src, int *err)
{
    int contd = 0;

    if (strlen(targ) + strlen(src) + 1 > MAXLINE) {
	*err = E_TOOLONG;
	return 0;
    } else {
	contd = ends_with_backslash(src);
	strcat(targ, src);
	if (contd) {
	    char *p = strrchr(targ, '\\');

	    if (p - targ > 0 && !isspace(*(p - 1))) {
		*p = ' ';
	    } else {
		*p = 0;
	    }
	}
    }	    

    return contd;
}

const char *console_varname_complete (const char *s)
{
    size_t n = strlen(s);
    int i;

    for (i=0; i<datainfo->v; i++) {
	if (!strncmp(s, datainfo->varname[i], n)) {
	    return datainfo->varname[i];
	}
    }

    return NULL;
}  

#define IS_BACKKEY(k) (k == GDK_BackSpace || k == GDK_Left)

static gint console_key_handler (GtkWidget *cview, GdkEventKey *key, 
				 char *cbuf)
{
    gint last_line, curr_line, line_pos;
    GtkTextIter iter, start, end;
    GtkTextBuffer *buf;
    gint ret = FALSE;

    /* retreat point in case we're not on the prompt line */
 start_again:

    /* where are we? */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));
    gtk_text_buffer_get_iter_at_mark(buf, &iter, gtk_text_buffer_get_insert(buf));
    curr_line = gtk_text_iter_get_line(&iter);    
    line_pos = gtk_text_iter_get_line_index(&iter);
    last_line = gtk_text_buffer_get_line_count(buf) - 1;

    if (IS_BACKKEY(key->keyval) && line_pos < 3) {
	/* if at start of command line, block backspacing */
	return TRUE;
    }

    if (curr_line < last_line) {
	/* if not on prompt line, return to (the end of) it */
	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_place_cursor(buf, &end);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(cview), TRUE);	
	goto start_again;
    }

    if (key->keyval == GDK_Return) {
	/* executes the command, unless backslash-continuation 
	   is happening */
	int contd = 0, err = 0;
	gchar *line;

	start = end = iter;
	gtk_text_iter_set_line_index(&start, 2);
	gtk_text_iter_forward_to_line_end(&end);
	line = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

	if (line != NULL) {
	    g_strstrip(line);
	    contd = command_continues(cbuf, line, &err);
	    g_free(line);
	}

	if (err) {
	    gui_errmsg(err);
	} else if (contd) {
	    console_insert_prompt(buf, &end, "\n> ");
	    console_scroll_to_end(cview, buf, &end);
	} else {
	    /* execute the completed command */
	    console_exec(cview, cbuf);
#ifdef G_OS_WIN32
	    gtk_window_present(GTK_WINDOW(gtk_widget_get_toplevel(cview)));
	    gtk_widget_grab_focus(cview);
#endif
	}

	key->keyval = GDK_End;
    } else if (key->keyval == GDK_Up || key->keyval == GDK_Down) {
	/* up/down arrows: navigate the command history */
	const char *histline;

	histline = fetch_history_line(key->keyval);

	if (histline != NULL || key->keyval == GDK_Down) {
	    start = end = iter;
	    gtk_text_iter_set_line_index(&start, 2);
	    gtk_text_iter_forward_to_line_end(&end);
	    gtk_text_buffer_delete(buf, &start, &end);
	    if (histline != NULL) {
		gtk_text_buffer_insert(buf, &start, histline, -1);
	    }
	}
	
	ret = TRUE;
    } else if (key->keyval == GDK_Tab) {
	/* tab completion for gretl commands, variable names */
	const char *targ = NULL;
	gchar *src = NULL;

	start = end = iter;

	if (!gtk_text_iter_starts_word(&start)) {
	    gtk_text_iter_backward_word_start(&start);
	}

	if (!gtk_text_iter_ends_word(&end)) {
	    gtk_text_iter_forward_word_end(&end);
	}

	src = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

	if (src != NULL && *src != '\0') {
	    if (gtk_text_iter_get_line_offset(&start) == 2) {
		/* first word on line */
		targ = gretl_command_complete(src);
	    } else {
		targ = console_varname_complete(src);
	    }
	    if (targ != NULL) {
		gtk_text_buffer_delete(buf, &start, &end);
		gtk_text_buffer_insert(buf, &start, targ, -1);
	    }
	}
	g_free(src);
	ret = TRUE;
    } else {
	GdkModifierType mods = widget_get_pointer_mask(cview);

	if (mods & GDK_CONTROL_MASK && 
	    gdk_keyval_to_upper(key->keyval) == GDK_A) {
	    /* Ctrl-A: go to start of typing area */
	    GtkTextIter where = iter;

	    gtk_text_iter_set_line_index(&where, 2);	
	    gtk_text_buffer_place_cursor(buf, &where);
	    ret = TRUE;
	}
    }

    return ret;
}

