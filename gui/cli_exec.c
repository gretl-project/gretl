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

#if GTK_MAJOR_VERSION > 2

static void modify_exec_button (windata_t *vwin, int to_spinner)
{
    GtkToolItem *eb = g_object_get_data(G_OBJECT(vwin->mbar), "exec_item");
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
	/* remove spinner and reinstate the exec button */
	idx = gtk_toolbar_get_item_index(GTK_TOOLBAR(vwin->mbar), si);
	gtk_spinner_stop(GTK_SPINNER(sp));
	gtk_container_remove(GTK_CONTAINER(vwin->mbar), GTK_WIDGET(si));
	gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), eb, idx);
    }
}

#else

/* GTK2: just toggle the sensitivity of the exec button */

static void modify_exec_button (windata_t *vwin, int to_spinner)
{
    GtkWidget *eb = g_object_get_data(G_OBJECT(vwin->mbar), "exec_item");

    if (to_spinner) {
	gtk_widget_set_sensitive(eb, FALSE);
    } else {
	gtk_widget_set_sensitive(eb, TRUE);
    }
}

#endif

typedef struct {
    gchar *cmd;           /* command-line (for Windows) */
    gchar **argv;         /* argument vector (non-Windows */
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
			    gchar *fname,
			    gchar *buf,
			    windata_t *vwin)
{
    ei->cmd = cmd;
    ei->argv = argv;
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
#else
    argv_free(ei->argv); /* note: ei->fname is part of this! */
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
#ifdef G_OS_WIN32
	ei->err = gretl_win32_pipe_output(ei->cmd, gretl_workdir(), ei->prn);
#else
	ei->err = gretl_pipe_output(ei->argv, gretl_workdir(), ei->prn);
#endif
    }
}

static void run_script_async (gchar *cmd,
			      gchar **argv,
			      gchar *fname,
			      windata_t *vwin)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GTask *task;

    exec_info_init(ei, cmd, argv, fname, NULL, vwin);
    ei->lang = 0; /* hansl */
    ei->err = 0;

    task = g_task_new(NULL, NULL, exec_script_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}

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
	run_script_async(cmd, NULL, inpname, vwin);
#else
	gchar **argv = g_malloc(4 * sizeof *argv);

	argv[0] = clipath;
	argv[1] = g_strdup("-x");
	argv[2] = inpname;
	argv[3] = NULL;
	run_script_async(NULL, argv, inpname, vwin);
#endif
    }

    return err;
}

static void editor_run_R_script (windata_t *vwin, gchar *buf)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GTask *task;

    exec_info_init(ei, NULL, NULL, NULL, buf, vwin);
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

    exec_info_init(ei, cmd, my_argv, NULL, NULL, vwin);
    ei->lang = lang;
    ei->err = 0;

    task = g_task_new(NULL, NULL, exec_script_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}