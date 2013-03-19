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

#define CDEBUG 0

/* file-scope globals */
static char **cmd_history;
static int hpos, hlines;
static gchar *hist0;
static ExecState *console_state;
static int command_entered;
static GtkWidget *console_main;
static GtkWidget *console_text;

static gint console_key_handler (GtkWidget *cview, GdkEventKey *key, 
				 gpointer p);

static void reset_console_globals (void)
{
    console_state = NULL;
    command_entered = 0;
}

static void command_history_init (void)
{
    hlines = hpos = 0;
    hist0 = NULL;
}

static void command_history_destroy (void)
{
    if (cmd_history != NULL) {
	strings_array_free(cmd_history, hlines);
	cmd_history = NULL;
    }
    
    g_free(hist0);
    hist0 = NULL;

    hlines = hpos = 0;
}

static ExecState *gretl_console_init (char *cbuf)
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

    gretl_exec_state_init(s, CONSOLE_EXEC, cbuf, 
			  get_lib_cmd(), model, prn);

    command_history_init();

    return s;
}

static void push_history_line (const char *line)
{
    strings_array_add(&cmd_history, &hlines, line);
    hpos = hlines;
}

static void console_beep (void)
{
#ifdef G_OS_WIN32
    MessageBeep(MB_ICONEXCLAMATION);
#else
    gdk_beep();
#endif
}

static const char *fetch_history_line (int keyval)
{
    static int beep;
    char *ret = NULL;

    if (keyval == GDK_Up) {
	if (hpos > 0) {
	    /* up is OK */
	    ret = cmd_history[--hpos];
	    beep = 0;
	} else {
	    /* can't go up */
	    console_beep();
	}
    } else if (keyval == GDK_Down) {
	if (hlines == 0) {
	    /* no history yet */
	    hpos = 0;
	    console_beep();
	} else if (hpos < hlines - 1) {
	    /* down is OK */
	    ret = cmd_history[++hpos];
	    beep = 0;
	} else {
	    /* can't go down */
	    hpos = hlines;
	    if (beep) {
		console_beep();
	    } else {
		beep = 1;
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
    gtk_text_buffer_delete_mark(buf, mark);
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
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(cview), 
				     on_last_line(cview));
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
    if (event->button == 2) {
	return console_paste_text(w, GDK_SELECTION_PRIMARY);
    }

    return FALSE;
}

enum {
    SAMPLE_RECORD,
    SAMPLE_CHECK
};

/* mechanism to check if a console action has altered the
   current sample information */

static int console_sample_handler (const DATASET *pdinfo, int code)
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

/* the two functions below are public because they are
   also used with the command 'minibuffer' */

void console_record_sample (const DATASET *pdinfo)
{
    console_sample_handler(pdinfo, SAMPLE_RECORD);
}

int console_sample_changed (const DATASET *pdinfo)
{
    return console_sample_handler(pdinfo, SAMPLE_CHECK);
}

static void print_result_to_console (GtkTextBuffer *buf,
				     GtkTextIter *iter,
				     ExecState *state)
{
    const char *prnbuf = gretl_print_get_buffer(state->prn);

    if (g_utf8_validate(prnbuf, -1, NULL)) {
	gtk_text_buffer_insert(buf, iter, prnbuf, -1);
    } else {
	gchar *trbuf = my_locale_to_utf8(prnbuf);

	fprintf(stderr, "console text did not validate as utf8\n");
	if (trbuf != NULL) {
	    gtk_text_buffer_insert(buf, iter, trbuf, -1);
	    g_free(trbuf);
	}
    } 

    gretl_print_reset_buffer(state->prn);   
}

static void console_insert_prompt (GtkTextBuffer *buf,
				   GtkTextIter *iter,
				   const char *prompt)
{
    gtk_text_buffer_insert_with_tags_by_name(buf, iter, prompt, -1,
					     "redtext", NULL);  
    gtk_text_buffer_place_cursor(buf, iter);
}

static int real_console_exec (ExecState *state)
{
    int err = 0;

#if CDEBUG
    fprintf(stderr, "*** real_console_exec: '%s'\n", state->line);
#endif

    push_history_line(state->line);

    state->flags = CONSOLE_EXEC;
    err = gui_exec_line(state, dataset);

    while (!err && gretl_execute_loop()) {
	err = gretl_loop_exec(state, dataset);
    }

#if CDEBUG
    fprintf(stderr, "*** real_console_exec returning %d\n", err);
#endif

    return err;
}

static const char *console_prompt (ExecState *s)
{
    if (s != console_state) {
	return "$ ";
    } else if (gretl_compiling_function() ||
	       gretl_compiling_loop()) {
	return "> ";
    } else {
	return "? ";
    }
}

/* called on receipt of a completed command line */

static void update_console (ExecState *state, GtkWidget *cview)
{
    GtkTextBuffer *buf;
    GtkTextIter iter;

    console_record_sample(dataset);

    if (state == console_state) {
	real_console_exec(state);
	if (state->cmd->ci == QUIT) {
	    return;
	}
    }

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));

    gtk_text_buffer_get_end_iter(buf, &iter);
    gtk_text_buffer_insert(buf, &iter, "\n", 1);

    if (!printing_is_redirected(state->prn)) {
	print_result_to_console(buf, &iter, state);
    }
    
    /* set up prompt for next command and scroll to it */
    console_insert_prompt(buf, &iter, console_prompt(state));
    console_scroll_to_end(cview, buf, &iter);

    /* update variable listing in main window if needed */
    if (check_dataset_is_changed()) {
	mark_dataset_as_modified();
	populate_varlist();
    }

    /* update sample info if needed */
    if (console_sample_changed(dataset)) {
	set_sample_label(dataset);
    }

#ifdef G_OS_WIN32
    gtk_window_present(GTK_WINDOW(gtk_widget_get_toplevel(cview)));
    gtk_widget_grab_focus(cview);
#endif
}

/* callback for function debugger to flush output */

static int console_update_callback (void *p)
{
    ExecState *state = (ExecState *) p;

#if CDEBUG
    fprintf(stderr, "*** console_update_callback\n");
#endif

    if (console_text == NULL) {
	return 1;
    } else {
	update_console(state, console_text); 
	return 0;
    }
}

/* command-getter, used here but also passed as callback to 
   the function debugger */

static int console_get_line (void *p)
{
    ExecState *state = (ExecState *) p;

#if CDEBUG
    fprintf(stderr, "*** console_get_line entered\n");
#endif

    if (console_text == NULL) {
	if (state != console_state) {
	    /* send 'continue' to debugger */
	    strcpy(state->line, "c");
	}
	return 1;
    }

    if (state != g_object_get_data(G_OBJECT(console_text), "ExecState")) {
#if CDEBUG
	if (state != console_state) {    
	    fprintf(stderr, "*** entered debugger\n");
	} else {
	    fprintf(stderr, "*** exited debugger\n");
	}
#endif
	g_object_set_data(G_OBJECT(console_text), "ExecState", state);
    }

    *state->line = '\0';

    /* wait for a command to be entered via the console */
    while (!command_entered) {
	if (gtk_events_pending()) {
	    gtk_main_iteration();
	}
	g_usleep(1000);
    }

    command_entered = 0;

#if CDEBUG
    fprintf(stderr, "*** console_get_line: '%s'\n", state->line);
#endif

    return 0;
}

int console_is_busy (void)
{
    if (console_main != NULL) {
	if (hlines > 0) {
	    gtk_window_present(GTK_WINDOW(console_main));
	    return 1;
	} else {
	    gtk_widget_destroy(console_main);
	    return 0;
	}
    }

    return 0;
}

static gboolean console_quit (GtkWidget *w, ExecState *state)
{
    if (state->cmd->ci != QUIT) {
	/* we're still in the command loop: defer the destroy */
	ExecState *curr;

	curr = g_object_get_data(G_OBJECT(console_text), "ExecState");
	if (curr != state) {
	    /* we're paused in the debugger */
	    strcpy(curr->line, "c");
	}	
	state->cmd->ci = QUIT;
	command_entered = 1;
	return TRUE;
    } else {
	/* OK, really destroy the console window */
	return FALSE;
    }
}

/* callback from menu/button: launches the console and remains
   in a command loop until done */

void gretl_console (void)
{
    char cbuf[MAXLINE];
    windata_t *vwin;
    GtkTextBuffer *buf;
    GtkTextIter iter;
    ExecState *state;
    const gchar *intro = 
	N_("gretl console: type 'help' for a list of commands");

    if (console_main != NULL) {
	gtk_window_present(GTK_WINDOW(console_main));
	return;
    }

    state = gretl_console_init(cbuf);
    if (state == NULL) {
	return;
    }

    console_state = state;

    vwin = console_window(78, 400);

    console_text = vwin->text;
    console_main = vwin->main;

    g_signal_connect(G_OBJECT(vwin->text), "paste-clipboard",
		     G_CALLBACK(console_paste_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
		     G_CALLBACK(console_click_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "button-release-event",
		     G_CALLBACK(console_mouse_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->text), "key-press-event",
		     G_CALLBACK(console_key_handler), NULL);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(console_quit), state);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), 
		     &console_main);
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(gtk_widget_destroyed), 
		     &console_text);

    g_object_set_data(G_OBJECT(vwin->text), "ExecState", state);

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_get_start_iter(buf, &iter);  

    /* insert intro string and first prompt */
    gtk_text_buffer_insert(buf, &iter, _(intro), -1);
    console_insert_prompt(buf, &iter, "\n? ");

    gtk_widget_grab_focus(vwin->text);

    set_debug_read_func(console_get_line);
    set_debug_output_func(console_update_callback);

    /* console command loop */
    while (state->cmd->ci != QUIT) {
	console_get_line(state);
	if (state->cmd->ci != QUIT) {
	    update_console(state, vwin->text);
	}
    }

    if (console_main != NULL) {
	/* the user actually typed quit/exit */
	gtk_widget_destroy(vwin->main);
    }

    command_history_destroy();
    gretl_print_destroy(state->prn);
    free(state);  

    set_debug_read_func(NULL);
    set_debug_output_func(NULL);
    reset_console_globals();

#if CDEBUG
    fprintf(stderr, "gretl_console: returning\n");
#endif
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
		*p = '\0';
	    }
	}
    }	    

    return contd;
}

const char *console_varname_complete (const char *s)
{
    size_t n = strlen(s);
    int i;

    for (i=0; i<dataset->v; i++) {
	if (!strncmp(s, dataset->varname[i], n)) {
	    return dataset->varname[i];
	}
    }

    return NULL;
}  

static gint console_complete_word (GtkTextBuffer *buf,
				   GtkTextIter *iter)
{
    GtkTextIter start, end;
    const char *targ = NULL;
    gchar *src;

    start = end = *iter;

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
	} else {
	    console_beep();
	}
    } 

    g_free(src);

    return TRUE;
}

static gchar *console_get_current_line (GtkTextBuffer *buf,
					GtkTextIter *iter)
{
    GtkTextIter start, end;

    start = end = *iter;
    gtk_text_iter_set_line_index(&start, 2);
    gtk_text_iter_forward_to_line_end(&end);

    return gtk_text_buffer_get_text(buf, &start, &end, FALSE);
}

#define IS_BACKKEY(k) (k == GDK_BackSpace || k == GDK_Left)

static gint console_key_handler (GtkWidget *cview, GdkEventKey *event, 
				 gpointer p)
{
    guint keyval = event->keyval;
    guint upkey = gdk_keyval_to_upper(keyval);
    GtkTextIter ins, end;
    GtkTextBuffer *buf;
    GtkTextMark *mark;
    gint ctrl = 0;

#ifdef MAC_NATIVE
    if (cmd_key(event)) {
	if (upkey == GDK_C || upkey == GDK_X) {
	    /* allow regular copy/cut behavior */
	    return FALSE;
	}
    }	
#endif

    if (event->state & GDK_CONTROL_MASK) {
	if (keyval == GDK_Control_L || keyval == GDK_Control_R) {
	    return FALSE;
	} else if (upkey == GDK_C || upkey == GDK_X) {
	    /* allow regular copy/cut behavior */
	    return FALSE;
	} else {
	    ctrl = 1;
	}
    }

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(cview));

    /* first find out where the insertion point and end are */
    mark = gtk_text_buffer_get_insert(buf);
    gtk_text_buffer_get_iter_at_mark(buf, &ins, mark);
    gtk_text_buffer_get_end_iter(buf, &end);

    /* if the insertion point is not on the last line, move it */
    if (gtk_text_iter_get_line(&ins) != gtk_text_iter_get_line(&end)) {
	gtk_text_buffer_place_cursor(buf, &end);
	gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(cview), TRUE);
	gtk_text_buffer_get_end_iter(buf, &ins);
    }

    if (keyval == GDK_Home && (event->state & GDK_SHIFT_MASK)) {
	/* "select to start of line" */
	GtkTextIter start = ins;

	gtk_text_iter_set_line_index(&start, 2);	
	gtk_text_buffer_select_range(buf, &start, &ins);
	return TRUE;
    }

    if (IS_BACKKEY(keyval)) {
	/* if we're at the start of the input line, block backspacing */
	if (gtk_text_iter_get_line_index(&ins) < 3) {
	    return TRUE;
	}
    } else if (keyval == GDK_Home || (ctrl && upkey == GDK_A)) {
	/* go to start of typing area */
	gtk_text_iter_set_line_index(&ins, 2);	
	gtk_text_buffer_place_cursor(buf, &ins);
	return TRUE;
    } 

    /* At this point 'ins' indicates the insertion point and
       'end' points to the end of the current line of input,
       These may or may not be the same thing.
    */

    if (keyval == GDK_Return) {
	/* execute the command, unless backslash-continuation 
	   is happening */
	ExecState *state;
	gchar *thisline;
	int contd = 0, err = 0;

	state = g_object_get_data(G_OBJECT(cview), "ExecState");
	thisline = console_get_current_line(buf, &ins);

	if (thisline != NULL) {
	    g_strstrip(thisline);
	    contd = command_continues(state->line, thisline, &err);
	    g_free(thisline);
	}

	if (err) {
	    gui_errmsg(err);
	} else if (contd) {
	    console_insert_prompt(buf, &end, "\n> ");
	    console_scroll_to_end(cview, buf, &end);
	} else {
	    /* request execution of the completed command */
	    command_entered = 1;
	}

	event->keyval = GDK_End;
	return FALSE;
    }

    if (keyval == GDK_Up || keyval == GDK_Down) {
	/* up/down arrows: navigate the command history */
	GtkTextIter start = ins;
	const char *histline;

	if (hpos == hlines && keyval == GDK_Up) {
	    g_free(hist0);
	    hist0 = console_get_current_line(buf, &ins);
	}

	histline = fetch_history_line(keyval);

	if (histline != NULL || keyval == GDK_Down) {
	    gtk_text_iter_set_line_index(&start, 2);
	    gtk_text_buffer_delete(buf, &start, &end);
	    if (histline != NULL) {
		gtk_text_buffer_insert(buf, &start, histline, -1);
	    } else if (hpos == hlines && hist0 != NULL) {
		gtk_text_buffer_insert(buf, &start, hist0, -1);
	    }
	}

	return TRUE;
    }

    if (keyval == GDK_Tab) {
	/* tab completion for gretl commands, variable names */
	return console_complete_word(buf, &ins);
    } 

    return FALSE;
}


