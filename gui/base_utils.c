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

    if (!has_suffix(fname, ".XXXXXX")) {
	strcat(fname, ".XXXXXX");
    }

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

PRN *gui_prn_new (void)
{
    int err = 0;
    PRN *prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);

    if (err) {
        gui_errmsg(err);
    }

    return prn;
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

int gretl_file_get_contents (const gchar *fname, gchar **contents,
			     gsize *size)
{
    GError *gerr = NULL;
    gboolean ok = 0;

    ok = g_file_get_contents(fname, contents, size, &gerr);

    if (gerr != NULL) {
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

#ifdef GRETL_EDIT

/* 2022-06-07: at present this is specific to gretl_edit,
   but it may be good to use it in the main gretl GUI.
*/

static void editor_run_other_script (windata_t *vwin,
				     gchar *cmd,
				     gchar **argv,
				     int lang);

#endif

#ifdef G_OS_WIN32

/* MS Windows variants of functions to exec some third-party
   programs */

static void win32_execute_script (gchar *cmd, int lang)
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
}

/* Windows version of run_foreign_script() */

static void run_foreign_script (windata_t *vwin, gchar *buf,
				int lang, gretlopt opt)
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

#ifdef GRETL_EDIT /* experimental */
	editor_run_other_script(vwin, cmd, NULL, lang);
#else
	win32_execute_script(cmd, lang);
	g_free(cmd);
#endif
	g_free(cmd);
    }
}

#else /* some non-Windows functions follow */

# ifndef GRETL_EDIT

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

# endif /* not GRETL_EDIT */

static void run_foreign_script (windata_t *vwin, gchar *buf,
				int lang, gretlopt opt)
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
	gchar *argv[6] = {0};

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

#ifdef GRETL_EDIT /* experimental */
	editor_run_other_script(vwin, NULL, argv, lang);
#else
	run_prog_sync(argv, lang);
#endif
    }
}

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

int browser_open (const char *url)
{
# if defined(OS_OSX)
    return osx_open_url(url);
# else
    return gretl_fork("Browser", url, NULL);
# endif
}

#endif /* !G_OS_WIN32 */

#ifdef GRETL_EDIT
# include "cli_exec.c"
#else
# include "R_exec.c"
#endif

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
    char fname[MAXLEN];
    fmode mode = 0;
    windata_t *vwin;

    if (scriptname == NULL) {
	/* the usual case */
	FILE *fp;

	sprintf(fname, "%sscript_tmp.XXXXXX", gretl_dotdir());
	fp = gretl_tempfile_open(fname);
	if (fp == NULL) {
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
	strcpy(fname, scriptname);
	mode = NULL_FILE;
    }

    if (action == EDIT_HANSL) {
        strcpy(scriptfile, fname);
    }

#ifdef GRETL_EDIT
    vsize *= 1.5;
#endif

    vwin = view_file(fname, 1, mode, SCRIPT_WIDTH, vsize, action);

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

static void ensure_newline_termination (gchar **ps)
{
    gchar *s = *ps;

    if (s[strlen(s)-1] != '\n') {
        gchar *tmp = g_strdup_printf("%s\n", s);

        g_free(s);
        *ps = tmp;
    }
}

static int lang_from_role (int role)
{
    struct {
	int role;
	int lang;
    } langmap[] = {
	{EDIT_R, LANG_R},
	{EDIT_OX, LANG_OX},
	{EDIT_OCTAVE, LANG_OCTAVE},
	{EDIT_PYTHON, LANG_PYTHON},
	{EDIT_JULIA, LANG_JULIA},
	{EDIT_STATA, LANG_STATA},
	{EDIT_DYNARE, LANG_OCTAVE}
    };
    int i;

    for (i=0; i<G_N_ELEMENTS(langmap); i++) {
	if (role == langmap[i].role) {
	    return langmap[i].lang;
	}
    }

    return 0;
}

static void real_run_script (windata_t *vwin, int silent)
{
#ifndef GRETL_EDIT
    gboolean selection = FALSE;
#endif
    gretlopt opt = OPT_NONE;
    gchar *prev_workdir = NULL;
    gchar *currdir = NULL;
    gchar *buf = NULL;
    int lang;

#ifdef GRETL_EDIT
    if (vwin->role == EDIT_HANSL) {
	buf = textview_get_hansl(vwin->text, 0);
    } else {
	buf = textview_get_text(vwin->text);
    }
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
            /* native script (possibly a gfn sample script) */
            gretl_set_script_dir(vwin->fname);
        }
    }

    if (vwin->role != EDIT_PKG_SAMPLE) {
        ensure_newline_termination(&buf);
    }

    if (vwin->role == EDIT_DYNARE) {
        opt = OPT_Y;
    }

    lang = lang_from_role(vwin->role);

#ifdef GRETL_EDIT
    if (vwin->role == EDIT_HANSL) {
	gretlcli_exec_script(vwin, buf);
    } else if (vwin->role == EDIT_R) {
	editor_run_R_script(vwin, buf);
	/* editor_run... takes ownership of @buf */
	buf = NULL;
    } else {
	run_foreign_script(vwin, buf, lang, opt);
    }
#else
    if (editing_hansl(vwin->role)) {
	run_native_script(vwin, buf, NULL, silent);
    } else if (vwin->role == EDIT_GP) {
        run_gnuplot_script(buf, vwin);
    } else if (vwin->role == EDIT_R) {
        run_R_script(buf, vwin);
    } else if (lang > 0) {
	run_foreign_script(vwin, buf, lang, opt);
    } else if (vwin->role == EDIT_LPSOLVE) {
	call_lpsolve_function(buf, vwin->fname, opt);
    } else if (vwin->role == EDIT_X12A) {
        run_x12a_script(buf);
    } else if (selection) {
	run_script_fragment(vwin, buf);
    } else {
        run_native_script(vwin, buf, NULL, silent);
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
    real_run_script(vwin, 0);
}

void run_script_silent (GtkWidget *w, windata_t *vwin)
{
    real_run_script(vwin, 1);
}
