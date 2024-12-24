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

/* gretl_edit.c : main for gretl script editor */

#include "gretl.h"
#include "libset.h"

#include "toolbar.h"
#include "winstack.h"
#include "dlgutils.h"
#include "tabwin.h"

#ifndef G_OS_WIN32
# include <unistd.h>
# include <sys/types.h>
# include "../pixmaps/gretl_edit.xpm"  /* program icon for X */
#else
# include <windows.h>
# include "gretlwin32.h"
#endif

#ifdef MAC_INTEGRATION
# include <gtkosxapplication.h>
#endif

#define GUI_DEBUG 0
#define WIN32_DEBUG 0

#if GUI_DEBUG
# include "version.h"
# include "build.h"
#endif

#if defined(__APPLE__) && defined(PKGBUILD)
# define ALT_MAC_STARTUP
#endif

/* update.c */
extern int update_query (void);

GtkTargetEntry gretl_drag_targets[] = {
    { "text/uri-list",  0, GRETL_FILENAME }
};

static int optver;
#ifdef G_OS_WIN32
static int optdebug;
#endif

static gchar *param_msg =
    N_("\nYou may supply the name of a data or script file on the command line");

static GOptionEntry options[] = {
    { "version", 'v', 0, G_OPTION_ARG_NONE, &optver,
      N_("print version information"), NULL },
    { NULL, '\0', 0, 0, NULL, NULL, NULL }
};

windata_t *mdata;
GtkWidget *editor;
char scriptfile[MAXLEN];
float gui_scale;

/* defaults for some options */
int winsize = FALSE;
int main_x = -1;
int main_y = -1;
int mainwin_width = 520;
int mainwin_height = 420;

#if defined(G_OS_WIN32)
char calculator[MAXSTR] = "calc.exe";
char latex[MAXSTR] = "pdflatex.exe";
char Rcommand[MAXSTR] = "RGui.exe";
#elif defined(__APPLE__)
char calculator[MAXSTR] = "/Applications/Calculator.app/Contents/MacOS/Calculator";
char latex[MAXSTR] = "pdflatex";
char Rcommand[MAXSTR] = "/Applications/R.app/Contents/MacOS/R";
#else
char Browser[MAXSTR] = "mozilla";
char calculator[MAXSTR] = "xcalc";
char latex[MAXSTR] = "pdflatex";
char viewpdf[MAXSTR] = "acroread";
char viewps[MAXSTR] = "gv";
char Rcommand[MAXSTR] = "xterm -e R";
#endif

static char tryfile[MAXLEN];

void set_tryfile (const char *fname)
{
    tryfile[0] = '\0';
    strncat(tryfile, fname, MAXLEN - 1);
}

char *get_tryfile (void)
{
    return tryfile;
}

void clear_tryfile (void)
{
    tryfile[0] = '\0';
}

int tryfile_is_set (void)
{
    return tryfile[0] != '\0';
}

static int script_type (const char *fname)
{
    if (has_suffix(fname, ".inp")) {
	return EDIT_HANSL;
    } else if (has_suffix(fname, ".R")) {
	return EDIT_R;
    } else if (has_suffix(fname, ".plt") ||
	       has_suffix(fname, ".gp")) {
	return EDIT_GP;
    } else if (has_suffix(fname, ".ox")) {
	return EDIT_OX;
    } else if (has_suffix(fname, ".m")) {
	return EDIT_OCTAVE;
    } else if (has_suffix(fname, ".py")) {
	return EDIT_PYTHON;
    } else if (has_suffix(fname, ".jl")) {
	return EDIT_JULIA;
    } else if (has_suffix(fname, ".do")) {
	return EDIT_STATA;
    } else if (has_suffix(fname, ".mod")) {
	return EDIT_DYNARE;
    } else if (has_suffix(fname, ".lp")) {
	return EDIT_LPSOLVE;
    } else {
	return 0;
    }
}

#if !defined(ENABLE_NLS)

static void real_nls_init (void)
{
    return;
}

#elif defined(G_OS_WIN32) && defined(PKGBUILD)

static void real_nls_init (void)
{
    char localedir[MAXSTR];

    gretl_build_path(localedir, gretl_home(), "locale", NULL);
    record_win32_locale(setlocale(LC_ALL, ""));
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#elif defined(__APPLE__)

#define LOCALE_CHECK 1

#if LOCALE_CHECK

#include <CoreFoundation/CoreFoundation.h>

/* Use this to check what we get from setlocale() ? */

static void macos_check_locale (void)
{
    CFLocaleRef cfloc = CFLocaleCopyCurrent();
    CFStringRef cfprop;
    const char *s;

    cfprop = (CFStringRef) CFLocaleGetValue(cfloc, kCFLocaleIdentifier);
    s = CFStringGetCStringPtr(cfprop, kCFStringEncodingASCII);
    if (s != NULL) {
	fprintf(stderr, "macos_check_locale: CF gave ID '%s'\n", s);
    }

    CFRelease(cfloc);
}

#endif /* LOCALE_CHECK */

static void real_nls_init (void)
{
    char *gretlhome = getenv("GRETL_HOME");
    char localedir[MAXSTR];
    char *p;

    if (gretlhome == NULL) {
	return;
    }

    strcpy(localedir, gretlhome);
    p = strstr(localedir, "share/gretl");
    if (p != NULL) {
	strcpy(p, "share/locale");
    }

    p = setlocale(LC_ALL, "");
    fprintf(stderr, "NLS init: setlocale() gave '%s'\n", p);
#if LOCALE_CHECK
    macos_check_locale();
#endif
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#else /* end OS X specific code */

/* regular *nix treatment of NLS -- also applies
   for non-package MSYS2 build on Windows */

static void real_nls_init (void)
{
    setlocale(LC_ALL, "");
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#endif /* NLS init variants */

void gui_nls_init (void)
{
    char *mylang = getenv("GRETL_LANG");

    if (mylang != NULL) {
	if (!g_ascii_strcasecmp(mylang, "english") ||
	    !g_ascii_strcasecmp(mylang, "C")) {
	    /* don't set up translation */
	    return;
	}
    }

    real_nls_init();
}

static void record_filearg (char *targ, const char *src)
{
    if (*src == '.') {
	gchar *cdir = g_get_current_dir();
	gchar *tmp = g_build_filename(cdir, src, NULL);

	strcpy(targ, tmp);
	g_free(cdir);
	g_free(tmp);
    } else {
	strcpy(targ, src);
    }
}

#if !defined(G_OS_WIN32) && GTK_MAJOR_VERSION == 3

/* cut out annoying runtime spew from GLib-GObject
   warning about deprecated stuff
*/

static void logtrap (const gchar *domain,
		     GLogLevelFlags level,
		     const gchar *msg,
		     gpointer p)
{
    if (strstr(msg, "deprecat") == NULL) {
	g_log_default_handler(domain, level, msg, p);
    }
}

static void quell_glib_spew (void)
{
    g_log_set_handler("GLib-GObject", G_LOG_LEVEL_WARNING,
		      (GLogFunc) logtrap, NULL);
}

#endif

#ifdef MAC_INTEGRATION

static GtkosxApplication *MacApp;

static gboolean app_should_quit_cb (GtkosxApplication *App, gpointer p)
{
    gboolean ret;

    fprintf(stderr, "app_should_quit_cb (editor %p)...\n", (void *) editor);
    ret = tabwin_exit_check(editor);
    fprintf(stderr, "  ret = %d\n", ret);
    return ret;
}

static void app_will_quit_cb (GtkosxApplication *App, gpointer p)
{
    gtk_main_quit();
}

static gboolean app_open_file_cb (GtkosxApplication *app,
				  gchar *path, gpointer p)
{
    if (path != NULL) {
	if (!strcmp(tryfile, path)) {
	    /* we're already on it? */
	    return TRUE;
	}
	set_tryfile(path);
	return open_tryfile(FALSE);
    } else {
	clear_tryfile();
	return TRUE;
    }
}

static void install_mac_signals (GtkosxApplication *App)
{
    g_signal_connect(App, "NSApplicationBlockTermination",
		     G_CALLBACK(app_should_quit_cb), NULL);
    g_signal_connect(App, "NSApplicationWillTerminate",
		     G_CALLBACK(app_will_quit_cb), NULL);
    g_signal_connect(App, "NSApplicationOpenFile",
		     G_CALLBACK(app_open_file_cb), NULL);
}

#endif /* MAC_INTEGRATION */

#if !defined(G_OS_WIN32) && !defined(__APPLE__)

static void protect_against_ubuntu (void)
{
    FILE *fp = fopen("/etc/os-release", "r");

    if (fp != NULL) {
	char line[80];

	while (fgets(line, sizeof line, fp)) {
	    if (strstr(line, "buntu")) {
		setenv("UBUNTU_MENUPROXY", "0", 1);
		break;
	    }
	}
	fclose(fp);
    }
}

#endif /* end Linux-specific */

#ifdef G_OS_WIN32

/* The point of the following special code: when gretl is
   invoked by the OS (via double-click on a file associated
   with gretl in the registry) the command-line may contain
   a "mixed language" filename that is not representable in
   the locale Code Page. Such a filename will appear in
   mangled form in the argv array, and we need to call on
   Windows APIs to get a UTF-16 version of this array.
*/

static void alt_gtk_init (int *pargc,
			  char ***pargv,
			  char *filearg,
			  GError **popterr)
{
    int argc_w = 0;
    int initted = 0;
    LPWSTR *argv_w;

    /* get args as UTF-16 */
    argv_w = CommandLineToArgvW(GetCommandLineW(), &argc_w);

    if (argv_w != NULL) {
	gchar **argv_u8 = calloc(argc_w, sizeof *argv_u8);
	gchar **origp = argv_u8; /* for use with g_free */
	int n_u8 = argc_w;
	int i, uerr = 0;

	/* for GTK, convert args to UTF-8 */
	for (i=0; i<argc_w && !uerr; i++) {
	    argv_u8[i] = g_utf16_to_utf8(argv_w[i], -1, NULL, NULL, NULL);
	    if (argv_u8[i] == NULL) {
		uerr = 1;
	    }
	}
	if (!uerr) {
	    gtk_init_with_args(&argc_w, &argv_u8, _(param_msg),
			       options, "gretl", popterr);
	    if (argc_w > 1 && *filearg == '\0') {
		strncat(filearg, argv_u8[1], MAXLEN - 1);
	    }
	    *pargc = argc_w; /* update (residual) arg count */
	    initted = 1;
	}
	/* clean up */
	for (i=0; i<n_u8; i++) {
	    g_free(origp[i]);
	}
	g_free(origp);
	LocalFree(argv_w);
    }

    if (!initted) {
	/* try fallback? */
	gtk_init_with_args(pargc, pargv, _(param_msg), options,
			   "gretl", popterr);
    }
}

#endif /* G_OS_WIN32 */

#ifdef ALT_MAC_STARTUP

#include "osx_env.c"

#endif /* specific to Mac package */

int main (int argc, char **argv)
{
    char auxname[MAXLEN];
    char filearg[MAXLEN];
    GError *opterr = NULL;

#if defined(G_OS_WIN32)
    /* this must come before NLS initialization */
    win32_set_gretldir();
#elif defined(ALT_MAC_STARTUP)
    osx_setup_paths();
#elif !defined(__APPLE__)
    /* Linux-specific */
    protect_against_ubuntu();
#endif

    gui_nls_init();

    *tryfile = '\0';
    *scriptfile = '\0';
    *auxname = '\0';
    *filearg = '\0';

#ifdef G_OS_WIN32
    alt_gtk_init(&argc, &argv, filearg, &opterr);
#else
    gtk_init_with_args(&argc, &argv, _(param_msg), options, "gretl", &opterr);
#endif
    if (opterr != NULL) {
	g_print("%s\n", opterr->message);
	exit(EXIT_FAILURE);
    }

#ifdef MAC_INTEGRATION
    MacApp = g_object_new(GTKOSX_TYPE_APPLICATION, NULL);
    install_mac_signals(MacApp);
    gtkosx_application_set_use_quartz_accelerators(MacApp, FALSE);
#endif

#ifdef G_OS_WIN32
    /* let's call this before doing libgretl_init */
# if WIN32_DEBUG
    gretl_win32_debug_init(1);
# else
    gretl_win32_debug_init(optdebug);
# endif
#elif GTK_MAJOR_VERSION == 3
    quell_glib_spew();
#endif

    libgretl_init();
    gretl_set_gui_mode();

#ifdef G_OS_WIN32
    gretl_win32_init(optdebug, 0);
#else
    gretl_config_init(0);
#endif

    if (optver) {
	gui_logo(NULL);
	exit(EXIT_SUCCESS);
    }

    helpfile_init();

    if (argc > 1 && *filearg == '\0') {
	/* If we have a residual unhandled command-line argument,
	   it should be the name of a file to be opened.
	*/
	strncat(filearg, argv[1], MAXLEN - 1);
    }

    if (*filearg != '\0') {
	/* Record what is presumably a filename argument
	   given on the command line.
	*/
	record_filearg(tryfile, filearg);
    }

#if defined(G_OS_WIN32)
    set_up_windows_look();
#elif defined(__APPLE__) && defined(HAVE_MAC_THEMES)
    set_up_mac_look();
#endif

    /* create the GUI */
    set_fixed_font(NULL, 1);
    set_app_font(NULL, 1);
    gretl_stock_icons_init();

    if (tryfile_is_set()) {
	open_tryfile(TRUE);
    } else {
	do_new_script(EDIT_HANSL, NULL, NULL);
    }

    /* Enter the event loop */
    gtk_main();

    libgretl_cleanup();

    return EXIT_SUCCESS;
}

#ifndef G_OS_WIN32

void set_wm_icon (GtkWidget *w)
{
    GdkPixbuf *icon = gdk_pixbuf_new_from_xpm_data(gretl_edit_xpm);

# ifdef MAC_INTEGRATION
    if (icon != NULL) {
	gtkosx_application_set_dock_icon_pixbuf(MacApp, icon);
	gtk_window_set_icon(GTK_WINDOW(w), icon);
	g_object_unref(icon);
    }
# else
    if (icon != NULL) {
	gtk_window_set_icon(GTK_WINDOW(w), icon);
	g_object_unref(icon);
    }
# endif
}

#endif

/* At start-up only: if we can't open a script file specified on the
   command line, give the user a choice between creating a new file
   (of the given name, if writable) and quitting the program.
*/

static gboolean maybe_open_script (int stype)
{
    if (gretl_test_fopen(tryfile, "r") == 0) {
	/* file can be opened OK */
	return do_open_script(stype);
    } else {
	int noname = 0;
	gchar *msg;
	gint resp;

	if (gretl_test_fopen(tryfile, "w") == 0) {
	    /* OK, @tryfile is at least writable */
	    msg = g_strdup_printf(_("Couldn't open %s"), tryfile);
	} else {
	    msg = g_strdup_printf(_("Couldn't write to %s"), tryfile);
	    noname = 1;
	}
	resp = script_start_dialog(msg);
	g_free(msg);
	if (resp == GTK_RESPONSE_OK) {
	    /* "New script" */
	    do_new_script(stype, NULL, noname ? NULL : tryfile);
	    return TRUE;
	} else {
	    /* "Quit" */
	    exit(EXIT_FAILURE);
	}
    }

    return FALSE;
}

gboolean open_tryfile (gboolean startup)
{
    int stype = script_type(tryfile);
    gboolean ret = FALSE;

    if (stype > 0) {
	/* We're looking at some sort of script */
	if (startup) {
	    return maybe_open_script(stype);
	} else {
	    return do_open_script(stype);
	}
    } else {
	/* only script files are acceptable */
	do_new_script(EDIT_HANSL, NULL, NULL);
	warnbox_printf(_("%s: not a recognized script file"), tryfile);
	return TRUE;
    }

    return ret;
}

void set_editor (GtkWidget *w)
{
    editor = w;
}
