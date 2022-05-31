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

#ifdef G_OS_WIN32 /* Windows-specific variant */

typedef struct {
    gchar *cmd;
    gchar *fname;
    windata_t *scriptwin;
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
	view_buffer(ei->prn, SCRIPT_WIDTH, 450, NULL, SCRIPT_OUT, ei->scriptwin);
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

static void win32_run_gretlcli_async (gchar *cmd,
				      gchar *fname,
				      windata_t *vwin)
{
    exec_info *ei = malloc(sizeof *ei);
    GTask *task;
    PRN *prn;

    bufopen(&prn);

    ei->cmd = cmd;
    ei->fname = g_strdup(fname);
    ei->scriptwin = vwin;
    ei->prn = prn;
    ei->err = 0;

    task = g_task_new(NULL, NULL, gretlcli_done, ei);
    g_task_set_task_data(task, ei, NULL);
    g_task_run_in_thread(task, exec_script_thread);
}

#else /* non-Windows variant */

typedef struct {
    gchar *fname;
    windata_t *scriptwin;
    gint fout;
    gint ferr;
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

static void run_gretlcli_async (char **argv, windata_t *scriptwin)
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

/* OS-invariant driver code */

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

    if (!err) {
#ifdef G_OS_WIN32
	gchar *cmd;

	cmd = g_strdup_printf("\"%s\" -x \"%s\"", clipath, inpname);
	modify_exec_button(vwin, 1);
	win32_run_gretlcli_async(cmd, inpname, vwin);
#else
	gchar *argv[4];

	argv[0] = (gchar *) clipath;
	argv[1] = (gchar *) "-x";
	argv[2] = (gchar *) inpname;
	argv[3] = NULL;
	modify_exec_button(vwin, 1);
	run_gretlcli_async(argv, vwin);
#endif
    }

    g_free(clipath);
    g_free(inpname);

    return err;
}

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
