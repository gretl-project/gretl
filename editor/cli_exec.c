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
#include "winstack.h"

#ifndef G_OS_WIN32

static gchar **argv_copy (gchar **argv)
{
    gchar **cpy = NULL;

    if (argv != NULL) {
	int i, n = g_strv_length(argv);

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
    gchar **env;          /* environment for executable, or NULL */
    gchar *fname;         /* input filename */
    gchar *buf;           /* script buffer (for R usage) */
    windata_t *scriptwin; /* source script-editor */
    PRN *prn;             /* for grabbing output */
    gchar *errout;        /* for grabbing stderr */
    int lang;             /* language of script */
    int err;              /* error code */
    int cancelled;        /* cancelled flag */
} exec_info;

static void exec_info_init (exec_info *ei,
			    gchar *cmd,
			    gchar **argv,
			    gchar *exepath,
			    gchar **env,
			    gchar *fname,
			    gchar *buf,
			    windata_t *vwin)
{
    ei->cmd = cmd;
    ei->argv = argv;
    ei->exepath = exepath;
    ei->env = env;
    ei->fname = fname;
    ei->buf = buf;
    ei->scriptwin = vwin;
    bufopen(&ei->prn);
    ei->cancelled = 0;
    ei->errout = NULL;
}

void exec_info_destroy (exec_info *ei)
{
    if (ei == NULL) {
	return;
    }

#ifdef G_OS_WIN32
    /* @argv is NULL in this case */
    g_free(ei->cmd);
    g_free(ei->fname);
    g_free(ei->exepath);
#else
    /* ei->clipath and ei->fname are included in @argv */
    g_strfreev(ei->argv);
#endif
    if (ei->env != NULL) {
	g_strfreev(ei->env);
    }
    g_free(ei->buf);
    free(ei);
}

gboolean viewer_has_stderr (windata_t *vwin)
{
    gboolean ret = FALSE;

    if (vwin != NULL && vwin->data != NULL) {
	exec_info *ei = vwin->data;

	ret = (ei->errout != NULL);
    }

    return ret;
}

void viewer_show_stderr (windata_t *vwin)
{
    const gchar *title = "gretl: stderr";
    exec_info *ei = vwin->data;
    PRN *prn;

    prn = gretl_print_new_with_gchar_buffer(ei->errout);
    ei->errout = NULL;
    view_buffer(prn, 84, 480, title, VIEW_STDERR, NULL);
}

static void reuse_editor_output_viewer (windata_t *vwin, PRN *prn)
{
    const char *newtext = gretl_print_get_buffer(prn);
    GtkTextBuffer *buf;

    /* replace previous content */
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_set_text(buf, "", -1);
    textview_set_text_colorized(vwin->text, newtext);
    cursor_to_top(vwin);

    gretl_print_destroy(prn);
    gtk_window_present(GTK_WINDOW(vwin->main));
}

static int strvs_differ (gchar **S0, gchar **S1)
{
    int i, n = g_strv_length(S0);

    if (g_strv_length(S1) != n) {
	return 1;
    } else {
	for (i=0; i<n; i++) {
	    if (strcmp(S0[i], S1[i])) {
		return 1;
	    }
	}
    }

    return 0;
}

static int exec_info_match (exec_info *ei0, exec_info *ei1)
{
    int match = 1;

    if (strcmp(ei0->exepath, ei1->exepath)) {
	/* using different builds of gretlcli */
	match = 0;
    } else {
	/* check for differing environments */
	int n = (ei0->env != NULL) + (ei1->env != NULL);

	if (n == 1) {
	    /* one env is special, the other is not */
	    match = 0;
	} else if (n == 2 && strvs_differ(ei0->env, ei1->env)) {
	    /* we have two differing special envs */
	    match = 0;
	}
    }

    return match;
}

/* Try for an existing output-viewer child that matches
   the current exec_info specification (gretlcli variant
   and environment setting, if any).
*/

static windata_t *get_matching_viewer (exec_info *ei)
{
    windata_t *parent = ei->scriptwin;
    windata_t *child = NULL;
    int i;

    for (i=0; i<parent->n_gretl_children; i++) {
	child = parent->gretl_children[i];
	if (child != NULL && child->data != NULL) {
	    if (exec_info_match(ei, child->data)) {
		return child;
	    }
	}
    }

    return NULL;
}

static void view_script_output (exec_info *ei)
{
    windata_t *vwin = get_matching_viewer(ei);
    GtkWidget *b;

    if (vwin != NULL) {
	/* found a specification-matching viewer */
	exec_info_destroy((exec_info *) vwin->data);
	vwin->data = ei;
	reuse_editor_output_viewer(vwin, ei->prn);
    } else {
	/* we need a new viewer */
	vwin = view_buffer(ei->prn, SCRIPT_WIDTH, 450, NULL,
			   SCRIPT_OUT, ei->scriptwin);
	if (vwin != NULL) {
	    vwin->data = ei;
	}
    }

    b = g_object_get_data(G_OBJECT(vwin->mbar), "stderr_item");
    gtk_widget_set_sensitive(b, ei->errout != NULL);
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

/* function called on completion of script execution */

static void exec_script_done (GObject *obj,
			      GAsyncResult *res,
			      gpointer data)
{
    exec_info *ei = data;
    int preserve_ei = 0;

    /* restore the "Run" icon in the editor window */
    modify_exec_button(ei->scriptwin, 0);

    if (ei->lang == LANG_STATA && !grab_stata_log(ei)) {
	; /* nothing to show, and error already reported */
    } else if (ei->cancelled) {
	infobox(_("Execution was cancelled"));
    } else if (got_printable_output(ei->prn)) {
	view_script_output(ei);
	preserve_ei = 1;
	ei->prn = NULL;
    } else if (ei->err) {
	gui_errmsg(ei->err);
    }

    /* clean up */
    if (ei->fname != NULL) {
	gretl_remove(ei->fname);
    }
    if (!preserve_ei) {
	exec_info_destroy(ei);
    }
}

static void print_env (exec_info *ei)
{
    int i;

    for (i=0; ei->env[i] != NULL; i++) {
        if (strncmp(ei->env[i], "HOME=", 5) &&
            strncmp(ei->env[i], "R_HOME=", 7)) {
            pprintf(ei->prn, "%s \\\n", ei->env[i]);
        }
    }
    pprintf(ei->prn, "%s\n", ei->exepath);
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
        if (ei->env != NULL) {
            print_env(ei);
        } else if (ei->exepath != NULL) {
	    pprintf(ei->prn, "using %s\n", ei->exepath);
	}
#ifdef G_OS_WIN32
	/* note: ei->env is _not_ handled here */
	ei->err = gretl_win32_pipe_output(ei->cmd, gretl_workdir(), ei->prn);
#else
	ei->err = gretl_pipe_output(ei->argv, ei->env, gretl_workdir(),
				    ei->prn, &ei->errout);
#endif
    }
}

static void task_cancel_callback (GCancellable *stopper,
				  gpointer data)
{
    gchar *fname = gretl_make_dotpath("exec.pid");
    FILE *fp = fopen(fname, "r");

    if (fp != NULL) {
        long pid;

        if (fscanf(fp, "%ld", &pid) == 1) {
	    exec_info *ei = data;
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
	    ei->cancelled = 1;
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
			      gchar **envp,
			      gchar *fname,
			      windata_t *vwin)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GCancellable *stopper;
    GTask *task;

    stopper = g_cancellable_new();
    g_cancellable_connect(stopper, G_CALLBACK(task_cancel_callback),
                          ei, NULL);
    g_cancellable_push_current(stopper);

    exec_info_init(ei, cmd, argv, exepath, envp, fname, NULL, vwin);
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

/* and to accommodate different environment settings */
static GList *gretlcli_envs;

static void populate_gretlcli_combo (GtkWidget *box,
				     GList **pL,
				     int paths)
{
    GList *L = *pL;

    if (L == NULL) {
        gchar *s;

	if (paths) {
	    s = g_strdup_printf("%sgretlcli", gretl_bindir());
	} else {
	    s = g_strdup("default");
	}
        *pL = g_list_prepend(*pL, s);

    }

    L = *pL;
    while (L != NULL) {
        combo_box_append_text(box, L->data);
        L = L->next;
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(box), 0);
}

/* For the benefit of settings.c: populate a combo box with
   the installed gretlcli path plus any others specified in
   the current session, with the most recently used path
   in the first position.
*/

void populate_gretlcli_path_combo (GtkWidget *box)
{
    populate_gretlcli_combo(box, &gretlcli_paths, 1);
}

/* Also for settings.c: populate a combo box with "default"
   plus any other set of environment variable settings
   specified in the current session, with the most recently
   used value in first position.
*/

void populate_gretlcli_env_combo (GtkWidget *box)
{
    populate_gretlcli_combo(box, &gretlcli_envs, 0);
}

static void set_gretlcli_param (GtkWidget *box, GList **pL)
{
    gchar *s = g_strstrip(combo_box_get_active_text(box));
    GList *L = g_list_first(*pL);
    int i = 0;

    while (L != NULL) {
        if (!strcmp(s, (gchar *) L->data)) {
            if (i > 0) {
                /* move @L to first position */
                *pL = g_list_remove_link(*pL, L);
                *pL = g_list_concat(L, *pL);
            }
            g_free(s);
            s = NULL;
            break;
        }
        i++;
        L = L->next;
    }

    if (s != NULL) {
        /* @s is not yet recorded in list */
        *pL = g_list_prepend(*pL, s);
    }
}

/* Called from settings.c: put whatever gretlcli path was
   selected, via the Editor tab under Preferences/General,
   into first position in @gretlcli_paths.
*/

void set_gretlcli_path (GtkWidget *box)
{
    set_gretlcli_param(box, &gretlcli_paths);
}

/* Called from settings.c: put whatever gretlcli env was
   selected, via the Editor tab under Preferences/General,
   into first position in @gretlcli_envs.
*/

void set_gretlcli_env (GtkWidget *box)
{
    set_gretlcli_param(box, &gretlcli_envs);
}

static gchar *get_cli_path (void)
{
    if (gretlcli_paths != NULL) {
        GList *L = g_list_first(gretlcli_paths);

        return g_strdup(L->data);
    } else {
        return g_strdup_printf("%sgretlcli", gretl_bindir());
    }
}

static char *should_add_basic (const char *s, gchar **envp)
{
    int i, n = g_strv_length(envp);
    int k = strlen(s);

    for (i=0; i<n; i++) {
	if (!strncmp(envp[i], s, k) && envp[i][k] == '=') {
	    return NULL; /* present already */
	}
    }

    return getenv(s);
}

static gchar **maybe_append_env_basics (gchar **envp)
{
    gchar **ret = envp;
    gchar *sh = NULL;
    gchar *sr = NULL;

    /* We'll assume it's safer to carry these two definitions across
       from the parent environment if they're not specified in the
       specific-to-gretlcli env.
    */
    if ((sh = should_add_basic("HOME", envp)) != NULL) {
	sh = g_strdup_printf("HOME=%s", sh);
    }
    if ((sr = should_add_basic("R_HOME", envp)) != NULL) {
	sr = g_strdup_printf("R_HOME=%s", sr);
    }

    if (sh != NULL || sr != NULL) {
	int n = g_strv_length(envp);
	int m = n + 1 + (sh != NULL) + (sr != NULL);

	ret = g_realloc(envp, m * sizeof *ret);
	if (sh != NULL) {
	    ret[n++] = sh;
	}
	if (sr != NULL) {
	    ret[n++] = sr;
	}
	ret[n] = NULL;
    }

    return ret;
}

static gchar **get_cli_env (void)
{
    gchar **ret = NULL;

    if (gretlcli_envs != NULL) {
        GList *L = g_list_first(gretlcli_envs);
	gchar *s = L->data;

	if (s != NULL && *s != '\0' && strcmp(s, "default")) {
	    ret = g_strsplit(s, " ", -1);
	    if (ret != NULL) {
		ret = maybe_append_env_basics(ret);
	    }
	}
    }

    return ret;
}

static int gretlcli_exec_script (windata_t *vwin, gchar *buf)
{
    gchar *inpname = gretl_make_dotpath("cli_XXXXXX.inp");
    FILE *fp = gretl_mktemp(inpname, "wb");
    gchar *clipath = get_cli_path();
    gchar **envp = get_cli_env();
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
	run_script_async(cmd, NULL, clipath, NULL, inpname, vwin);
#else
	gchar **argv = g_malloc(4 * sizeof *argv);

	argv[0] = clipath;
	argv[1] = g_strdup("-x");
	argv[2] = inpname;
	argv[3] = NULL;
	run_script_async(NULL, argv, clipath, envp, inpname, vwin);
#endif
    }

    return err;
}

static void editor_run_R_script (windata_t *vwin, gchar *buf)
{
    exec_info *ei = calloc(1, sizeof *ei);
    GTask *task;

    exec_info_init(ei, NULL, NULL, NULL, NULL, NULL, buf, vwin);
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

    exec_info_init(ei, cmd, my_argv, NULL, NULL, NULL, NULL, vwin);
    ei->lang = lang;
    ei->err = 0;

    task = g_task_new(NULL, NULL, exec_script_done, ei);
    g_task_set_task_data(task, ei, NULL);
    modify_exec_button(vwin, 1);
    g_task_run_in_thread(task, exec_script_thread);
}
