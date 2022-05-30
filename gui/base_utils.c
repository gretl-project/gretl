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
#include "gretl_enums.h"
#include "gretl_foreign.h"
#include "textbuf.h"
#include "base_utils.h"

#ifndef GRETL_EDIT
# include "library.h"
# include "filelists.h"
# include "fnsave.h"
# include "gpt_control.h"
#endif

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#ifdef OS_OSX
# include "osx_open.h"
#endif

gchar *gretl_window_title (const char *s)
{
    if (s != NULL) {
	return g_strdup_printf("gretl: %s", s);
    } else {
	return g_strdup("gretl: untitled");
    }
}

int copyfile (const char *src, const char *dest)
{
    FILE *srcfd, *destfd;
    char buf[GRETL_BUFSIZE];
    size_t n;

    if (!strcmp(src, dest)) {
	return 0;
    }

    if ((srcfd = gretl_fopen(src, "rb")) == NULL) {
	file_read_errbox(src);
	return E_FOPEN;
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
	file_write_errbox(dest);
	fclose(srcfd);
	return E_FOPEN;
    }

    while ((n = fread(buf, 1, sizeof buf, srcfd)) > 0) {
	fwrite(buf, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}

FILE *gretl_tempfile_open (char *fname)
{
    FILE *fp;

    strcat(fname, ".XXXXXX");
    fp = gretl_mktemp(fname, "w+");

    if (fp == NULL) {
	errbox(_("Couldn't open temp file"));
    }

    return fp;
}

int bufopen (PRN **pprn)
{
    static int has_minus = -1;
    int err = 0;

    *pprn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
        gui_errmsg(err);
    } else {
        if (has_minus < 0) {
            /* check for Unicode minus sign, U+2212 */
            has_minus = font_has_symbol(fixed_font, 0x2212);
        }
        if (has_minus > 0) {
            gretl_print_set_has_minus(*pprn);
        }
    }

    return err;
}

void set_wait_cursor (GdkWindow **pcwin)
{
    GdkDisplay *disp = gdk_display_get_default();
    GdkWindow *w;

    if (*pcwin == NULL) {
	gint x, y;

	*pcwin = w = gdk_display_get_window_at_pointer(disp, &x, &y);
    } else {
	w = *pcwin;
    }

    if (w != NULL) {
	GdkCursor *cursor = gdk_cursor_new(GDK_WATCH);

	if (cursor != NULL) {
	    gdk_window_set_cursor(w, cursor);
	    gdk_display_sync(disp);
	    gdk_cursor_unref(cursor);
	}
    }
}

void unset_wait_cursor (GdkWindow *cwin)
{
    if (cwin != NULL) {
	gdk_window_set_cursor(cwin, NULL);
    }
}

/* by using gretl_set_window_modal() we make the main
   window visibly insensitive */

static int modcount;

static void increment_modal_count (GtkWidget *w)
{
    if (modcount == 0) {
        gtk_widget_set_sensitive(mdata->main, FALSE);
    }

    modcount++;
}

static void decrement_modal_count (GtkWidget *w, gpointer p)
{
    if (modcount > 0) {
        modcount--;
    }

    if (modcount == 0) {
        gtk_widget_set_sensitive(mdata->main, TRUE);
    }
}

void gretl_set_window_modal (GtkWidget *w)
{
    gtk_window_set_modal(GTK_WINDOW(w), TRUE);
    increment_modal_count(w);
    g_signal_connect(G_OBJECT(w), "destroy",
                     G_CALLBACK(decrement_modal_count),
                     NULL);
}

void gretl_set_window_quasi_modal (GtkWidget *w)
{
    increment_modal_count(w);
    g_signal_connect(G_OBJECT(w), "destroy",
                     G_CALLBACK(decrement_modal_count),
                     NULL);
}

void dummy_call (void)
{
    errbox(_("Sorry, this item not yet implemented!"));
}

void nomem (void)
{
    errbox(_("Out of memory!"));
}

void *mymalloc (size_t size)
{
    void *mem;

    if ((mem = malloc(size)) == NULL) {
	nomem();
    }

    return mem;
}

void *myrealloc (void *ptr, size_t size)
{
    void *mem;

    if ((mem = realloc(ptr, size)) == NULL) {
	nomem();
    }

    return mem;
}

const char *path_last_slash_const (const char *path)
{
    return (const char *) strrslash(path);
}

char *gretl_basename (char *dest, const char *src, int addscore)
{
    const char *p = path_last_slash_const(src);
    const char *s = (p != NULL)? p + 1 : src;

    if (dest == NULL) {
	/* allocate return value */
	size_t len = strlen(s) + 1;

	dest = calloc(len, 1);
	addscore = 0;
    }

    if (dest != NULL) {
	strcpy(dest, s);
    }

    if (addscore) {
	/* double any underscores in @dest */
	char mod[MAXSTR];
	int n = strlen(dest);
	int i, j = 0;

	for (i=0; i<=n; i++) {
	    if (dest[i] == '_') {
		mod[j++] = '_';
	    }
	    mod[j++] = dest[i];
	}
	strcpy(dest, mod);
    }

    return dest;
}

char *double_underscores (char *targ, const char *src)
{
    char *p = targ;

    while (*src) {
	if (*src == '_') {
	    *p++ = '_';
	    *p++ = '_';
	} else {
	    *p++ = *src;
	}
	src++;
    }
    *p = '\0';

    return targ;
}

gchar *double_underscores_new (const char *src)
{
    const char *s = src;
    gchar *ret;
    int n = 0;

    while (*s && (s = strchr(s, '_')) != NULL) {
	n++;
	s++;
    }

    ret = g_malloc(strlen(src) + n + 1);

    return double_underscores(ret, src);
}

char *adjust_fontspec_string (char *targ, const char *src,
			      int mod)
{
    char *p, c0, c1;

    strcpy(targ, src);

    if (mod == ADD_COMMA) {
	c0 = ' ';
	c1 = ',';
    } else {
	c0 = ',';
	c1 = ' ';
    }

    p = strrchr(targ, c0);
    if (p != NULL && isdigit(p[1])) {
	*p = c1;
    }

    return targ;
}

int get_string_width (const gchar *str)
{
    GtkWidget *w;
    PangoLayout *pl;
    PangoContext *pc;
    gint width;

    w = gtk_label_new(NULL);
    pc = gtk_widget_get_pango_context(w);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, str, -1);
    pango_layout_get_pixel_size(pl, &width, NULL);

    gtk_widget_destroy(w);
    g_object_unref(G_OBJECT(pl));

    return width;
}

void *gui_get_plugin_function (const char *funcname)
{
    void *func;

    func = get_plugin_function(funcname);
    if (func == NULL) {
	errbox(gretl_errmsg_get());
    }

    return func;
}

static int got_printable_output (PRN *prn)
{
    int ret = 0;

    if (prn == NULL) {
	warnbox(_("No output was produced"));
    } else {
	const char *buf = gretl_print_get_buffer(prn);

	if (string_is_blank(buf)) {
	    warnbox(_("No output was produced"));
	    gretl_print_destroy(prn);
	} else {
	    ret = 1;
	}
    }

    return ret;
}

#ifdef G_OS_WIN32

/* MS Windows variants of functions to exec some third-party
   programs */

# ifndef GRETL_EDIT

static void win32_run_R_sync (const char *buf, gretlopt opt)
{
    PRN *prn = NULL;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    err = execute_R_buffer(buf, dataset, opt, prn);

    if (err) {
	gui_errmsg(err);
    } else {
	view_buffer(prn, 78, 350, _("gretl: script output"),
		    PRINT, NULL);
    }
}

# endif

void win32_execute_script (gchar *cmd, int lang, windata_t *scriptwin)
{
    PRN *prn = NULL;
    int err = 0;

    if (bufopen(&prn)) {
	return;
    }

    if (lang == LANG_STATA) {
	gchar *buf = NULL;

	gretl_chdir(gretl_workdir());
	remove("gretltmp.log");
	err = gretl_spawn(cmd);

	if (g_file_get_contents("gretltmp.log", &buf, NULL, NULL)) {
	    pputs(prn, buf);
	    g_free(buf);
	    pputc(prn, '\n');
	}
    } else {
	err = gretl_win32_pipe_output(cmd, gretl_dotdir(), prn);
    }

    if (got_printable_output(prn)) {
	/* note: this check destroys @prn on failure */
	view_buffer(prn, 78, 350, _("gretl: script output"), PRINT, NULL);
    } else if (err) {
	gui_errmsg(err);
    }

    g_free(cmd);
}

void run_foreign_script (gchar *buf, int lang, gretlopt opt)
{
#ifdef GRETL_EDIT
    DATASET *dataset = NULL;
#endif
    const char *fname = NULL;
    int err;

    opt |= OPT_G;

    /* note: as things stand, the @fname we obtain here
       (composed in gretl_foreign.c) will be in the locale
       encoding, ready to pass on the Windows command line
       as in "foreign.exe fname"; this composite string is
       given to gretl_spawn() or gretl_win32_pipe_output()
       below. 2022-04-21: this comment is outdated, no?
    */

    err = write_gretl_foreign_script(buf, lang, opt, dataset, &fname);

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *cmd = NULL;

	if (lang == LANG_OX) {
	    cmd = g_strdup_printf("\"%s\" \"%s\"", gretl_oxl_path(), fname);
	} else if (lang == LANG_PYTHON) {
	    cmd = g_strdup_printf("\"%s\" \"%s\"", gretl_python_path(), fname);
	} else if (lang == LANG_JULIA) {
	    cmd = g_strdup_printf("\"%s\" \"%s\"", gretl_julia_path(), fname);
	} else if (lang == LANG_STATA) {
	    cmd = g_strdup_printf("\"%s\" /e do \"%s\"", gretl_stata_path(), fname);
	} else if (lang == LANG_OCTAVE) {
	    cmd = g_strdup_printf("\"%s\" -q \"%s\"", gretl_octave_path(), fname);
	}

	win32_execute_script(cmd, lang, NULL);
	g_free(cmd);
    }
}

#else /* some non-Windows functions follow */

int browser_open (const char *url)
{
# if defined(OS_OSX)
    return osx_open_url(url);
# else
    return gretl_fork("Browser", url, NULL);
# endif
}

#ifndef GRETL_EDIT

/* Start an R session in asynchronous (interactive) mode.
   Note that there's a separate win32 function for this
   in gretlwin32.c. We don't do interactive when in editor
   mode.
*/

static void start_R_async (void)
{
    char *s0 = NULL, *s1 = NULL, *s2 = NULL;
    int n = -1;

    s0 = mymalloc(64);
    s1 = mymalloc(32);
    s2 = mymalloc(32);

    if (s0 != NULL && s1 != NULL && s2 != NULL) {
	*s0 = *s1 = *s2 = '\0';
	/* probably "xterm -e R" or similar */
	n = sscanf(Rcommand, "%63s %31s %31s", s0, s1, s2);
    }

    if (n == 0) {
	errbox(_("No command was supplied to start R"));
    } else if (n > 0) {
	char *supp1 = "--no-init-file";
	char *supp2 = "--no-restore-data";
	gchar *argv[6];
	GError *error = NULL;
	gboolean ok;
	int i = 0;

	argv[i++] = s0;
	if (n > 1) {
	    argv[i++] = s1;
	}
	if (n > 2) {
	    argv[i++] = s2;
	}
	argv[i++] = supp1;
	argv[i++] = supp2;
	argv[i++] = NULL;

	ok = g_spawn_async(NULL,
			   argv,
			   NULL,
			   G_SPAWN_SEARCH_PATH,
			   NULL,
			   NULL,
			   NULL,
			   &error);

	if (error != NULL) {
	    errbox(error->message);
	    g_error_free(error);
	} else if (!ok) {
	    gui_errmsg(E_EXTERNAL);
	    g_error_free(error);
	}
    }

    free(s0);
    free(s1);
    free(s2);
}

#endif /* not GRETL_EDIT */

/* Run R, Ox, etc., in synchronous (batch) mode and display the
   results in a gretl window: non-Windows variant.
*/

static void run_prog_sync (char **argv, int lang)
{
    gchar *sout = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    PRN *prn = NULL;

    if (lang == LANG_STATA) {
	/* control location of Stata log file */
	gretl_chdir(gretl_workdir());
	remove("gretltmp.log");
    }

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		 NULL, NULL, &sout, &errout,
		 &status, &gerr);

    if (gerr != NULL) {
	errbox(gerr->message);
	g_error_free(gerr);
    } else if (status != 0) {
	if (errout != NULL && *errout != '\0') {
	    if (strlen(errout) < MAXLEN) {
		errbox(errout);
	    } else {
		bufopen(&prn);
		pputs(prn, errout);
	    }
	} else if (sout != NULL && *sout != '\0') {
	    if (strlen(sout) < MAXLEN) {
		errbox(sout);
	    } else {
		bufopen(&prn);
		pputs(prn, sout);
	    }
	} else {
	    errbox_printf("%s exited with status %d", argv[0], status);
	}
    } else if (lang == LANG_STATA) {
	/* read log file */
	gchar *buf = NULL;

	if (g_file_get_contents("gretltmp.log", &buf, NULL, NULL)) {
	    bufopen(&prn);
	    pputs(prn, buf);
	    g_free(buf);
	    pputc(prn, '\n');
	}
    } else {
	/* read good old stdout */
	if (!string_is_blank(sout)) {
	    bufopen(&prn);
	    pputs(prn, sout);
	}
    }

    if (got_printable_output(prn)) {
	view_buffer(prn, 78, 350, _("gretl: script output"), PRINT, NULL);
    }

    g_free(sout);
    g_free(errout);
}

#ifndef GRETL_EDIT

static void run_R_sync (void)
{
    gchar *argv[] = {
	"R",
	"--no-save",
	"--no-init-file",
	"--no-restore-data",
	"--slave",
	NULL
    };

    run_prog_sync(argv, LANG_R);
}

#endif

void run_foreign_script (gchar *buf, int lang, gretlopt opt)
{
#ifdef GRETL_EDIT
    DATASET *dataset = NULL;
#endif
    const char *fname = NULL;
    int err;

    opt |= OPT_G;

    err = write_gretl_foreign_script(buf, lang, opt, dataset, &fname);

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *argv[6];

	if (lang == LANG_OCTAVE && (opt & OPT_Y)) {
	    gretl_chdir(gretl_workdir());
	}

	if (lang == LANG_OX) {
	    argv[0] = (gchar *) gretl_oxl_path();
	    argv[1] = (gchar *) fname;
	    argv[2] = NULL;
	} else if (lang == LANG_PYTHON) {
	    argv[0] = (gchar *) gretl_python_path();
	    argv[1] = (gchar *) fname;
	    argv[2] = NULL;
	} else if (lang == LANG_JULIA) {
	    argv[0] = (gchar *) gretl_julia_path();
	    argv[1] = (gchar *) fname;
	    argv[2] = NULL;
	} else if (lang == LANG_STATA) {
	    argv[0] = (gchar *) gretl_stata_path();
	    argv[1] = (gchar *) "-q";
	    argv[2] = (gchar *) "-b";
	    argv[3] = (gchar *) "do";
	    argv[4] = (gchar *) fname;
	    argv[5] = NULL;
	} else if (lang == LANG_OCTAVE) {
	    argv[0] = (gchar *) gretl_octave_path();
	    argv[1] = (gchar *) "-q";
	    argv[2] = (gchar *) fname;
	    argv[3] = NULL;
	}

	run_prog_sync(argv, lang);
    }
}

#endif /* !G_OS_WIN32 */

static void modify_exec_button (windata_t *vwin, int to_spinner)
{
    GtkToolItem *eb = g_object_get_data(G_OBJECT(vwin->mbar), "exec_button");
    GtkToolItem *si = g_object_get_data(G_OBJECT(vwin->mbar), "spin_item");
    GtkWidget *sp = gtk_tool_button_get_icon_widget(GTK_TOOL_BUTTON(si));
    int idx;

    if (to_spinner) {
	/* replace exit button with "wait" spinner */
	idx = gtk_toolbar_get_item_index(GTK_TOOLBAR(vwin->mbar), eb);
	gtk_container_remove(GTK_CONTAINER(vwin->mbar), GTK_WIDGET(eb));
	gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), si, idx);
	gtk_widget_show_all(GTK_WIDGET(si));
	gtk_spinner_start(GTK_SPINNER(sp));
    } else {
	/* reinstate the exec button */
	idx = gtk_toolbar_get_item_index(GTK_TOOLBAR(vwin->mbar), si);
	gtk_spinner_stop(GTK_SPINNER(sp));
	gtk_container_remove(GTK_CONTAINER(vwin->mbar), GTK_WIDGET(si));
	gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), eb, idx);
    }
}

#ifdef G_OS_WIN32 /* back to Windows-specific */

typedef struct {
    gchar *cmd;
    gchar *fname;
    windata_t *vwin;
    PRN *prn;
    int err;
} exec_info;

static void gretlcli_done (GObject *obj,
			   GAsyncResult *res,
			   gpointer data)
{
    exec_info *ei = data;

    modify_exec_button(ei->scriptwin, 0);
    if (ei->err == 0 /* got_printable_output(ei->prn) */) {
	view_buffer(ei->prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, ei->vwin);
	ei->prn = NULL;
    } else if (ei->err) {
	gui_errmsg(ei->err);
    }

    gretl_remove(ei->fname);
    g_free(ei->fname);
    g_free(ei->cmd);
    gretl_print_destroy(ei->prn);
    free(ei);
}

static void exec_script_thread (GTask *task,
				gpointer obj,
				gpointer task_data,
				GCancellable *c)
{
    exec_info *ei = g_task_get_task_data(task);

    ei->err = gretl_win32_pipe_output(ei->cmd, gretl_dotdir(), ei->prn);
}

void win32_run_gretlcli_async (gchar *cmd,
			       gchar *fname,
			       windata_t *vwin)
{
    exec_info *ei = malloc(sizeof *ei);
    GTask *task;
    PRN *prn;

    bufopen(&prn);

    ei->cmd = cmd;
    ei->fname = g_strdup(fname);
    ei->vwin = vwin;
    ei->prn = prn;
    ei->err = 0;

    task = g_task_new(NULL, NULL, gretlcli_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}

#else /* non-Windows variant */

typedef struct {
    windata_t *scriptwin;
    gint fout;
    gint ferr;
    gchar *fname;
} exec_info;

static void gretlcli_done (GPid pid, gint status, gpointer p)
{
    exec_info *ei = p;

    modify_exec_button(ei->scriptwin, 0);

    if (ei != NULL) {
	char buf[4096];
	size_t bs = sizeof buf - 1;
	ssize_t got;
	PRN *prn;

	bufopen(&prn);

	if (ei->fout > 0) {
	    while ((got = read(ei->fout, buf, bs)) > 0) {
		buf[got] = '\0';
		pputs(prn, buf);
	    }
	    close(ei->fout);
	}
	if (ei->ferr > 0) {
	    if (status != 0 && (got = read(ei->ferr, buf, bs)) > 0) {
		buf[got] = '\0';
		pputs(prn, buf);
	    }
	    close(ei->ferr);
	}
	view_buffer(prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, ei->scriptwin);
	gretl_remove(ei->fname);
	g_free(ei->fname);
	free(ei);
    }

    g_spawn_close_pid(pid);
}

void run_gretlcli_async (char **argv, windata_t *scriptwin)
{
    gboolean run;
    gint ferr = -1;
    gint fout = -1;
    GPid pid = 0;
    GError *gerr = NULL;

    run = g_spawn_async_with_pipes(NULL, argv, NULL,
				   G_SPAWN_SEARCH_PATH |
				   G_SPAWN_DO_NOT_REAP_CHILD,
				   NULL, NULL, &pid, NULL,
				   &fout, &ferr,
				   &gerr);

    if (gerr != NULL) {
	errbox(gerr->message);
	g_error_free(gerr);
    } else if (!run) {
	errbox(_("gretlcli command failed"));
    } else if (pid > 0) {
	exec_info *ei = calloc(1, sizeof *ei);

	ei->scriptwin = scriptwin;
	ei->fout = fout;
	ei->ferr = ferr;
	ei->fname = g_strdup(argv[2]);
	g_child_watch_add(pid, gretlcli_done, ei);
    }
}

#endif /* run_gretlcli_async variants */

#ifndef GRETL_EDIT

/* driver for starting R, either interactive or in batch mode,
   not used when in editor mode
*/

void start_R (const char *buf, int send_data, int interactive)
{
    gretlopt Ropt = OPT_G;
    int err = 0;

    if (send_data && !data_status) {
	warnbox(_("Please open a data file first"));
	return;
    }

    if (interactive) {
	Ropt |= OPT_I;
    }

    if (send_data) {
	Ropt |= OPT_D;
	if (annual_data(dataset) || quarterly_or_monthly(dataset)) {
	    const char *opts[] = {
		N_("multiple time series object"),
		N_("data frame")
	    };
	    int resp;

	    resp = radio_dialog(NULL, _("Send data as"), opts, 2, 0, 0, NULL);
	    if (resp < 0) {
		return;
	    } else if (resp == 1) {
		Ropt |= OPT_F;
	    }
	}
    }

    /* On Windows in non-interactive mode, don't write
       these files here; that will be handled later
    */
#ifdef G_OS_WIN32
    if (interactive) {
	err = write_gretl_R_files(buf, dataset, Ropt);
    }
#else
    err = write_gretl_R_files(buf, dataset, Ropt);
#endif

    if (err) {
	gui_errmsg(err);
	delete_gretl_R_files();
    } else if (interactive) {
#ifdef G_OS_WIN32
	win32_start_R_async();
#else
	start_R_async();
#endif
    } else {
	/* non-interactive */
#ifdef G_OS_WIN32
	win32_run_R_sync(buf, Ropt);
#else
	run_R_sync();
#endif
    }
}

#endif /* not GRETL_EDIT */

static void verbose_gerror_report (GError *gerr, const char *src)
{
    fprintf(stderr, "GError details from %s\n"
	    " message: '%s'\n domain = %d, code = %d\n",
	    src, gerr->message, gerr->domain, gerr->code);
}

int gretl_file_get_contents (const gchar *fname, gchar **contents,
			     gsize *size)
{
    GError *gerr = NULL;
    gboolean ok = 0;

    ok = g_file_get_contents(fname, contents, size, &gerr);

    if (gerr != NULL) {
	verbose_gerror_report(gerr, "g_file_get_contents");
	errbox(gerr->message);
	g_error_free(gerr);
    }

    return ok ? 0 : E_FOPEN;
}

const char *print_today (void)
{
    static char timestr[16];
    struct tm *local;
    time_t t;

    t = time(NULL);
    local = localtime(&t);
    strftime(timestr, 15, "%Y-%m-%d", local);

    return timestr;
}

/* We mostly use this for checking whether the font described by @desc
   has the Unicode minus sign (0x2212), which looks better than a
   simple dash if it's available.
*/

int font_has_symbol (PangoFontDescription *desc, int symbol)
{
    GtkWidget *widget;
    PangoContext *context = NULL;
    PangoLayout *layout = NULL;
    PangoLanguage *lang = NULL;
    PangoCoverage *coverage = NULL;
    int ret = 0;

    if (desc == NULL) {
	return 0;
    }

    widget = gtk_label_new("");
    if (g_object_is_floating(widget)) {
	g_object_ref_sink(widget);
    }

    context = gtk_widget_get_pango_context(widget);
    if (context == NULL) {
	gtk_widget_destroy(widget);
	return 0;
    }

    layout = pango_layout_new(context);
    lang = pango_language_from_string("eng");

    if (layout != NULL && lang != NULL) {
	PangoFont *font = pango_context_load_font(context, desc);

	if (font != NULL) {
	    coverage = pango_font_get_coverage(font, lang);
	    if (coverage != NULL) {
		ret = (pango_coverage_get(coverage, symbol) == PANGO_COVERAGE_EXACT);
		pango_coverage_unref(coverage);
	    }
	    g_object_unref(font);
	}
    }

    g_object_unref(G_OBJECT(layout));
    g_object_unref(G_OBJECT(context));
    gtk_widget_destroy(widget);

    return ret;
}

#ifdef GRETL_EDIT

/* code specific to editor mode but in common between Windows and others */

static void editor_run_R_script (const char *buf, gretlopt opt)
{
    PRN *prn = NULL;

    if (bufopen(&prn)) {
	return;
    }

    execute_R_buffer(buf, NULL, OPT_G | OPT_T, prn);
    if (got_printable_output(prn)) {
	view_buffer(prn, 78, 350, _("gretl: script output"), PRINT, NULL);
    }
}

#endif

#ifndef G_OS_WIN32 /* non-Windows variants */

int gretl_fork (const char *progvar, const char *arg,
		const char *opt)
{
    const char *prog = NULL;
    gchar *argv[4] = {NULL, NULL, NULL, NULL};
    GError *err = NULL;
    gboolean run;

#ifndef GRETL_EDIT
    if (!strcmp(progvar, "calculator")) {
	prog = calculator;
	goto do_argv;
    }
#endif

#ifndef OS_OSX
    if (!strcmp(progvar, "Browser")) {
	prog = Browser;
    } else if (!strcmp(progvar, "viewpdf")) {
	prog = viewpdf;
    } else if (!strcmp(progvar, "viewps")) {
	prog = viewps;
    } else {
	prog = progvar;
    }
#endif

    if (prog == NULL) {
	errbox_printf("Internal error: variable %s is undefined", progvar);
	return 1;
    }

#ifndef GRETL_EDIT
 do_argv:
#endif

    argv[0] = g_strdup(prog);
    if (opt != NULL) {
	argv[1] = g_strdup(arg);
	argv[2] = g_strdup(opt);
    } else if (arg != NULL) {
	argv[1] = g_strdup(arg);
    }

    run = g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
			NULL, NULL, NULL, &err);

    if (err != NULL) {
	errbox(err->message);
	if (err->domain == G_SPAWN_ERROR &&
	    err->code == G_SPAWN_ERROR_NOENT) {
	    preferences_dialog(TAB_PROGS, progvar, mdata->main);
	}
	g_error_free(err);
    }

    g_free(argv[0]);
    g_free(argv[1]);
    g_free(argv[2]);

    return !run;
}

#endif /* end Windows vs non-Windows */

/* do_new_script(): passing a non-NULL @scriptname is a means
   of creating a new script with a name pre-given by the user;
   this applies only when we get the name of a non-existent
   script on the command line. Otherwise the new script gets
   a temporary name in the user's dotdir.
*/

void do_new_script (int code, const char *buf,
		    const char *scriptname)
{
    int action = (code == FUNC)? EDIT_HANSL : code;
    int vsize = SCRIPT_HEIGHT;
    fmode mode = 0;
    windata_t *vwin;
    gchar *fname;

    if (scriptname == NULL) {
	/* the usual case */
	FILE *fp;

	fname = gretl_make_dotpath("script_tmp");
	fp = gretl_tempfile_open(fname);
	if (fp == NULL) {
	    g_free(fname);
	    return;
	} else {
	    if (buf != NULL) {
		fputs(buf, fp);
	    } else if (code == FUNC) {
		fputs("function \n\nend function\n", fp);
	    } else if (code == EDIT_OX) {
		fputs("#include <oxstd.h>\n\n", fp);
		fputs("main()\n{\n\n}\n", fp);
	    }
	    fclose(fp);
	}
	mode = TMP_FILE;
    } else {
	/* special startup case */
	fname = g_strdup(scriptname);
	mode = NULL_FILE;
    }

    if (action == EDIT_HANSL) {
        strcpy(scriptfile, fname);
    }

#ifdef GRETL_EDIT
    vsize *= 1.5;
#endif

    vwin = view_file(fname, 1, mode, SCRIPT_WIDTH, vsize, action);
    g_free(fname);

    if (buf != NULL && *buf != '\0') {
        mark_vwin_content_changed(vwin);
    }
}

void new_script_callback (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    int etype = EDIT_HANSL;

    if (!strcmp(s, "GnuplotScript")) {
        etype = EDIT_GP;
    } else if (!strcmp(s, "RScript")) {
        etype = EDIT_R;
    } else if (!strcmp(s, "OxScript")) {
        etype = EDIT_OX;
    } else if (!strcmp(s, "OctaveScript")) {
        etype = EDIT_OCTAVE;
    } else if (!strcmp(s, "PyScript")) {
        etype = EDIT_PYTHON;
    } else if (!strcmp(s, "StataScript")) {
        etype = EDIT_STATA;
    } else if (!strcmp(s, "JuliaScript")) {
        etype = EDIT_JULIA;
    } else if (!strcmp(s, "DynareScript")) {
        etype = EDIT_DYNARE; /* FIXME not reached */
    } else if (!strcmp(s, "lpsolveScript")) {
	etype = EDIT_LPSOLVE;
    }

    do_new_script(etype, NULL, NULL);
}

gboolean do_open_script (int action)
{
    char *fname = get_tryfile();
    int err = gretl_test_fopen(fname, "r");

    if (err) {
	file_read_errbox(fname);
#ifndef GRETL_EDIT
        if (action == EDIT_HANSL) {
            delete_from_filelist(FILE_LIST_SESSION, fname);
            delete_from_filelist(FILE_LIST_SCRIPT, fname);
        }
#endif
        return FALSE;
    }

    if (action == EDIT_HANSL) {
        strcpy(scriptfile, fname);
#ifndef GRETL_EDIT
	mkfilelist(FILE_LIST_SCRIPT, scriptfile, 1);
	gretl_set_script_dir(scriptfile);
#endif
        if (has_system_prefix(scriptfile, SCRIPT_SEARCH)) {
            view_script(scriptfile, 0, VIEW_SCRIPT);
        } else {
            view_script(scriptfile, 1, EDIT_HANSL);
        }
    } else {
        view_script(fname, 1, action);
    }

    return TRUE;
}

#ifdef GRETL_EDIT

/* edit mode special: execute a script (which may be entirely native
   or may contain "foreign" code) via gretlcli
*/

static int gretlcli_exec_script (windata_t *vwin, gchar *buf)
{
    gchar *clipath = g_strdup_printf("%sgretlcli", gretl_bindir());
    gchar *inpname = gretl_make_dotpath("cli_tmp.inp");
    FILE *fp = gretl_fopen(inpname, "wb");
    int err = 0;

    if (fp == NULL) {
	file_read_errbox(inpname);
	err = E_FOPEN;
    } else {
	fputs(buf, fp);
	fclose(fp);
    }

    modify_exec_button(vwin, 1);

    if (!err) {
#ifdef G_OS_WIN32
	gchar *cmd;

	cmd = g_strdup_printf("\"%s\" -x \"%s\"", clipath, inpname);
	win32_run_gretlcli_async(cmd, inpname, vwin);
#else
	gchar *argv[4];

	argv[0] = (gchar *) clipath;
	argv[1] = (gchar *) "-x";
	argv[2] = (gchar *) inpname;
	argv[3] = NULL;
	run_gretlcli_async(argv, vwin);
#endif
    }

    g_free(clipath);
    g_free(inpname);

    return err;
}

#endif

static void ensure_newline_termination (gchar **ps)
{
    gchar *s = *ps;

    if (s[strlen(s)-1] != '\n') {
        gchar *tmp = g_strdup_printf("%s\n", s);

        g_free(s);
        *ps = tmp;
    }
}

static void real_run_script (GtkWidget *w, windata_t *vwin,
			     int silent)
{
#ifndef GRETL_EDIT
    gboolean selection = FALSE;
#endif
    gretlopt opt = OPT_NONE;
    gchar *prev_workdir = NULL;
    gchar *currdir = NULL;
    gchar *buf = NULL;

#ifdef GRETL_EDIT
    buf = textview_get_text(vwin->text);
#else
    if (vwin->role == EDIT_GP ||
        vwin->role == EDIT_R ||
        vwin->role == EDIT_OX ||
        vwin->role == EDIT_OCTAVE ||
        vwin->role == EDIT_PYTHON ||
        vwin->role == EDIT_JULIA ||
        vwin->role == EDIT_DYNARE ||
	vwin->role == EDIT_LPSOLVE ||
        vwin->role == EDIT_STATA ||
        vwin->role == EDIT_X12A) {
        buf = textview_get_text(vwin->text);
    } else if (vwin->role == EDIT_PKG_SAMPLE) {
        buf = package_sample_get_script(vwin);
    } else {
        buf = textview_get_selection_or_all(vwin->text, &selection);
    }
#endif

    if (buf == NULL || *buf == '\0') {
        warnbox("No commands to execute");
	g_free(buf);
        return;
    }

    if (vwin->fname[0] != '\0' &&
        strstr(vwin->fname, "script_tmp") == NULL &&
        g_path_is_absolute(vwin->fname)) {
        /* There's a "real" full filename in place */
        if (editing_alt_script(vwin->role)) {
            /* For an "alt" script we'll temporarily reset
               workdir to its location so we're able to pick up
               any data files it may reference. We'll also arrange
               to revert the working directory once we're done.
            */
            gchar *dname = g_path_get_dirname(vwin->fname);

            currdir = g_get_current_dir();
            prev_workdir = g_strdup(gretl_workdir());
            gretl_set_path_by_name("workdir", dname);
            gretl_chdir(dname);
            g_free(dname);
        } else if (vwin->role != EDIT_GP && vwin->role != EDIT_PKG_SAMPLE) {
            /* native script */
            gretl_set_script_dir(vwin->fname);
        }
    }

    if (vwin->role != EDIT_PKG_SAMPLE) {
        ensure_newline_termination(&buf);
    }

    if (vwin->role == EDIT_DYNARE) {
        opt = OPT_Y;
    }

#ifdef GRETL_EDIT
    if (vwin->role == EDIT_R) {
	editor_run_R_script(buf, opt);
    } else if (vwin->role == EDIT_OX) {
        run_foreign_script(buf, LANG_OX, opt);
    } else if (vwin->role == EDIT_OCTAVE) {
        run_foreign_script(buf, LANG_OCTAVE, opt);
    } else if (vwin->role == EDIT_PYTHON) {
        run_foreign_script(buf, LANG_PYTHON, opt);
    } else if (vwin->role == EDIT_JULIA) {
        run_foreign_script(buf, LANG_JULIA, opt);
    } else if (vwin->role == EDIT_DYNARE) {
        run_foreign_script(buf, LANG_OCTAVE, opt);
    } else if (vwin->role == EDIT_STATA) {
        run_foreign_script(buf, LANG_STATA, opt);
    } else {
	gretlcli_exec_script(vwin, buf);
    }
#else
    if (vwin->role == EDIT_GP) {
        run_gnuplot_script(buf, vwin);
    } else if (vwin->role == EDIT_R) {
        run_R_script(buf, vwin);
    } else if (vwin->role == EDIT_OX) {
        run_foreign_script(buf, LANG_OX, opt);
    } else if (vwin->role == EDIT_OCTAVE) {
        run_foreign_script(buf, LANG_OCTAVE, opt);
    } else if (vwin->role == EDIT_PYTHON) {
        run_foreign_script(buf, LANG_PYTHON, opt);
    } else if (vwin->role == EDIT_JULIA) {
        run_foreign_script(buf, LANG_JULIA, opt);
    } else if (vwin->role == EDIT_DYNARE) {
        run_foreign_script(buf, LANG_OCTAVE, opt);
    } else if (vwin->role == EDIT_LPSOLVE) {
	call_lpsolve_function(buf, vwin->fname, opt);
    } else if (vwin->role == EDIT_STATA) {
        run_foreign_script(buf, LANG_STATA, opt);
    } else if (vwin->role == EDIT_X12A) {
        run_x12a_script(buf);
    } else if (selection) {
	run_script_fragment(vwin, buf);
    } else {
        run_native_script(vwin, buf, silent);
    }
#endif

    g_free(buf);

    if (prev_workdir != NULL) {
        gretl_set_path_by_name("workdir", prev_workdir);
        g_free(prev_workdir);
    }
    if (currdir != NULL) {
        gretl_chdir(currdir);
        g_free(currdir);
    }
}

void do_run_script (GtkWidget *w, windata_t *vwin)
{
    real_run_script(w, vwin, 0);
}

void run_script_silent (GtkWidget *w, windata_t *vwin)
{
    real_run_script(w, vwin, 1);
}
