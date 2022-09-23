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

/* code specific to gretl_edit: execute a hansl script via gretlcli */

#include "dlgutils.h"

#ifndef G_OS_WIN32

static void argv_free (gchar **argv)
{
    if (argv != NULL) {
	int i;

	for (i=0; argv[i] != NULL; i++) {
	    g_free(argv[i]);
	}
	g_free(argv);
    }
}

static gchar **argv_copy (gchar **argv)
{
    gchar **cpy = NULL;

    if (argv != NULL) {
	int i, n = 0;

	for (i=0; argv[i] != NULL; i++) {
	    n++;
	}
	cpy = g_malloc0((n + 1) * sizeof *cpy);
	for (i=0; i<n; i++) {
	    cpy[i] = g_strdup(argv[i]);
	}
    }

    return cpy;
}

#endif

/* switch between EXEC and STOP icons */

static void modify_exec_button (windata_t *vwin, int to_killer)
{
    GtkWidget *eb = g_object_get_data(G_OBJECT(vwin->mbar), "exec_item");
    gchar *tip = NULL;

    if (to_killer) {
	widget_set_int(vwin->mbar, "exec_is_kill", 1);
	gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(eb), GTK_STOCK_STOP);
        tip = g_strdup_printf("%s (%s)", _("Stop script"), "Ctrl+Alt+K");
	gtk_tool_item_set_tooltip_text(GTK_TOOL_ITEM(eb), tip);
    } else {
	widget_set_int(vwin->mbar, "exec_is_kill", 0);
	gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(eb), GTK_STOCK_EXECUTE);
        tip = g_strdup_printf("%s (%s)", _("Run"), "Ctrl+R");
        gtk_tool_item_set_tooltip_text(GTK_TOOL_ITEM(eb), tip);
    }

    g_free(tip);
}

typedef struct {
    gchar *cmd;           /* command-line (for Windows) */
    gchar **argv;         /* argument vector (non-Windows) */
    gchar *exepath;       /* path to executable */
    gchar *fname;         /* input filename */
    gchar *buf;           /* script buffer (for R usage) */
    windata_t *scriptwin; /* source script-editor */
    PRN *prn;             /* for grabbing output */
    int lang;             /* language of script */
    int err;              /* error flag */
} exec_info;

static void exec_info_init (exec_info *ei,
			    gchar *cmd,
			    gchar **argv,
			    gchar *exepath,
			    gchar *fname,
			    gchar *buf,
			    windata_t *vwin)
{
    ei->cmd = cmd;
    ei->argv = argv;
    ei->exepath = exepath;
    ei->fname = fname;
    ei->buf = buf;
    ei->scriptwin = vwin;
    bufopen(&ei->prn);
}

static int grab_stata_log (exec_info *ei)
{
    gchar *buf = NULL;

    if (g_file_get_contents("gretltmp.log", &buf, NULL, NULL)) {
	pputs(ei->prn, buf);
	g_free(buf);
	pputc(ei->prn, '\n');
	return 1;
    } else {
	file_read_errbox("gretltmp.log");
	gretl_print_destroy(ei->prn);
	return 0;
    }
}

static void trash_stata_log (void)
{
    const char *wdir = gretl_workdir();
    char *logname = gretl_build_path((char *) wdir, "gretltmp.log");

    gretl_remove(logname);
    free(logname);
}

/* callback on completion of script execution */

static void exec_script_done (GObject *obj,
			      GAsyncResult *res,
			      gpointer data)
{
    exec_info *ei = data;

    modify_exec_button(ei->scriptwin, 0);

    if (ei->lang == LANG_STATA && !grab_stata_log(ei)) {
	goto no_output;
    }

    if (got_printable_output(ei->prn)) {
	view_buffer(ei->prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, ei->scriptwin);
	ei->prn = NULL;
    } else if (ei->err) {
	gui_errmsg(ei->err);
    }

 no_output:

    if (ei->fname != NULL) {
	gretl_remove(ei->fname);
    }
#ifdef G_OS_WIN32
    g_free(ei->cmd);
    g_free(ei->fname);
    g_free(ei->exepath);
#else
    argv_free(ei->argv); /* note: ei->clipath and ei->fname are included */
#endif
    g_free(ei->buf);
    free(ei);
}

static void exec_script_thread (GTask *task,
				gpointer obj,
				gpointer task_data,
				GCancellable *c)
{
    exec_info *ei = g_task_get_task_data(task);

    if (ei->lang == LANG_R) {
	ei->err = execute_R_buffer(ei->buf, NULL, OPT_G | OPT_T, ei->prn);
    } else {
	if (ei->exepath != NULL) {
	    pprintf(ei->prn, "using %s\n", ei->exepath);
	}
#ifdef G_OS_WIN32
	ei->err = gretl_win32_pipe_output(ei->cmd, gretl_workdir(), ei->prn);
#else
	ei->err = gretl_pipe_output(ei->argv, gretl_workdir(), ei->prn);
#endif
    }
}

static void task_cancel_callback (void)
{
    gchar *fname = gretl_make_dotpath("exec.pid");
    FILE *fp = fopen(fname, "r");

    if (fp != NULL) {
        long pid;

        if (fscanf(fp, "%ld", &pid) == 1) {
#ifdef G_OS_WIN32
            DWORD dw = (DWORD) pid;
            HANDLE h;

            h = OpenProcess(PROCESS_TERMINATE, FALSE, dw);
            if (h != NULL) {
                TerminateProcess(h, ERROR_SUCCESS);
            }
#else
            kill(pid, SIGKILL);
#endif
        }
        fclose(fp);
        gretl_remove(fname);
    }

    g_free(fname);
}

void cancel_run_script (void)
{
    GCancellable *stopper = g_cancellable_get_current();

    if (stopper != NULL) {
        g_cancellable_cancel(stopper);
        g_cancellable_pop_current(stopper);
    }
}

static void run_script_async (gchar *cmd,
			      gchar **argv,
			      gchar *exepath,
			      gchar *fname,
			      windata_t *vwin)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GCancellable *stopper;
    GTask *task;

    stopper = g_cancellable_new();
    g_cancellable_connect(stopper, G_CALLBACK(task_cancel_callback),
                          NULL, NULL);
    g_cancellable_push_current(stopper);

    exec_info_init(ei, cmd, argv, exepath, fname, NULL, vwin);
    ei->lang = 0; /* hansl */
    ei->err = 0;

    task = g_task_new(NULL, stopper, exec_script_done, ei);
    g_task_set_check_cancellable(task, TRUE);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}

/* list to accommodate builds of gretl other than the installed one */
static GList *gretlcli_paths;

/* For the benefit of settings.c: populate a combo box with
   the installed gretlcli path plus any others specified in
   the current session, with the most recently used path
   in the first position.
*/

void populate_gretlcli_path_combo (GtkWidget *box)
{
    GList *L;

    if (gretlcli_paths == NULL) {
        gchar *s = g_strdup_printf("%sgretlcli", gretl_bindir());

        gretlcli_paths = g_list_prepend(gretlcli_paths, s);
    }

    L = gretlcli_paths;
    while (L != NULL) {
        combo_box_append_text(box, L->data);
        L = L->next;
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(box), 0);
}

/* Called from settings.c: put whatever gretlcli path was
   selected, via the Editor tab under Preferences/General,
   into first position in @gretlcli_paths.
*/

void set_gretlcli_path (GtkWidget *box)
{
    gchar *path = combo_box_get_active_text(box);
    GList *L = g_list_first(gretlcli_paths);
    int i = 0;

    while (L != NULL) {
        if (!strcmp(path, (gchar *) L->data)) {
            if (i > 0) {
                /* move @L to first position */
                gretlcli_paths = g_list_remove_link(gretlcli_paths, L);
                gretlcli_paths = g_list_concat(L, gretlcli_paths);
            }
            g_free(path);
            path = NULL;
            break;
        }
        i++;
        L = L->next;
    }

    if (path != NULL) {
        /* @path is not yet recorded in @gretlcli_paths */
        gretlcli_paths = g_list_prepend(gretlcli_paths, path);
    }
}

static int gretlcli_exec_script (windata_t *vwin, gchar *buf)
{
    gchar *inpname = gretl_make_dotpath("cli_XXXXXX.inp");
    FILE *fp = gretl_mktemp(inpname, "wb");
    gchar *clipath;
    int err = 0;

    if (gretlcli_paths != NULL) {
        GList *L = g_list_first(gretlcli_paths);

        clipath = g_strdup(L->data);
    } else {
        clipath = g_strdup_printf("%sgretlcli", gretl_bindir());
    }

    if (fp == NULL) {
	file_read_errbox(inpname);
	err = E_FOPEN;
    } else {
	fputs(buf, fp);
	fclose(fp);
    }

    if (!err) {
#ifdef G_OS_WIN32
	gchar *cmd;

	cmd = g_strdup_printf("\"%s\" -x \"%s\"", clipath, inpname);
	run_script_async(cmd, NULL, clipath, inpname, vwin);
#else
	gchar **argv = g_malloc(4 * sizeof *argv);

	argv[0] = clipath;
	argv[1] = g_strdup("-x");
	argv[2] = inpname;
	argv[3] = NULL;
	run_script_async(NULL, argv, clipath, inpname, vwin);
#endif
    }

    return err;
}

static void editor_run_R_script (windata_t *vwin, gchar *buf)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GTask *task;

    exec_info_init(ei, NULL, NULL, NULL, NULL, buf, vwin);
    ei->lang = LANG_R;
    ei->err = 0;

    task = g_task_new(NULL, NULL, exec_script_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}

static void editor_run_other_script (windata_t *vwin,
				     gchar *cmd,
				     gchar **argv,
				     int lang)
{
    exec_info *ei = calloc(1, sizeof *ei);
    gchar **my_argv = NULL;
    GTask *task;

#ifndef G_OS_WIN32
    my_argv = argv_copy(argv);
#endif

    if (lang == LANG_STATA) {
	/* ensure we don't get stale output */
	trash_stata_log();
    }

    exec_info_init(ei, cmd, my_argv, NULL, NULL, NULL, vwin);
    ei->lang = lang;
    ei->err = 0;

    task = g_task_new(NULL, NULL, exec_script_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}
