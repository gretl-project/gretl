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

/* settings.c for gretl */

#include "gretl.h"
#include "filelists.h"
#include "gretl_www.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "menustate.h"
#include "session.h"
#include "textbuf.h"
#include "ssheet.h"

#include "libset.h"
#include "version.h"
#include "texprint.h"

#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#else
# include <sys/stat.h>
# include <fcntl.h>
# include <errno.h>
# include "gtkfontselhack.h"
#endif

#if !defined(G_OS_WIN32)
char rcfile[MAXLEN];
#endif

char dbproxy[21];
int use_proxy;

static void make_prefs_tab (GtkWidget *notebook, int tab);
static void apply_changes (GtkWidget *widget, gpointer data);
#ifndef G_OS_WIN32
static void read_rc (void);
#endif

/* font handling */
#ifdef G_OS_WIN32
static char fixedfontname[MAXLEN] = "Courier New 10";
#else
# ifdef OSX_BUILD
static char fixedfontname[MAXLEN] = "Luxi Mono 12";
# else
static char fixedfontname[MAXLEN] = "Monospace 10";
# endif
#endif

#if defined(G_OS_WIN32)
static char appfontname[MAXLEN] = "tahoma 8";
#elif !defined(USE_GNOME)
# ifdef OSX_BUILD
static char appfontname[MAXLEN] = "Luxi Sans 12";
# else
static char appfontname[MAXLEN] = "Sans 10";
# endif
#endif

PangoFontDescription *fixed_font;

static int usecwd;
static int shellok;
static int manpref;
char gpcolors[64];
static char datapage[24];
static char scriptpage[24];

static int hc_by_default;
static int langpref;
static char hc_xsect[5] = "HC1";
static char hc_tseri[5] = "HAC";
static char hc_panel[9] = "Arellano";
static char hc_garch[5] = "QML";

#ifdef ENABLE_NLS
static int lcnumeric = 1;
#endif

#ifdef G_OS_WIN32
extern int use_wimp;
#endif

#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
char midiplayer[MAXSTR] = "timidity -ig";
#endif

typedef enum {
    ROOTSET  = 1 << 0,
    USERSET  = 1 << 1, 
    BOOLSET  = 1 << 2,
    INTSET   = 1 << 3,
    LISTSET  = 1 << 4,
    RADIOSET = 1 << 5,
    INVISET  = 1 << 6,
    FIXSET   = 1 << 7,  /* setting fixed by admin (Windows network use) */
    MACHSET  = 1 << 8,  /* "local machine" setting */
    BROWSER  = 1 << 9   /* wants "Browse" button */
} rcflags;

typedef struct {
    char *key;         /* config file variable name */
    char *description; /* How the field will show up in the options dialog */
    char *link;        /* in case of radio button pair, alternate string */
    void *var;         /* pointer to variable */
    rcflags flags;     /* ROOTSET user string
			  USERSET root string
			  MACHSET local machine string
			  BOOLSET boolean (user)
                          INTSET integer (user)
			  LISTSET user string, from fixed menu
                          RADIOSET user int, from fixed menu
			  INVISET not visible in Prefs dialog
		       */
    int len;           /* storage size for string variable (also see Note) */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
} RCVAR;

/* Note: actually "len" above is overloaded: if an rc_var is of
   type BOOLSET and not part of a radio group, then a non-zero value
   for len will link the var's toggle button with the sensitivity of
   the preceding rc_var's entry field.  For example, the "use_proxy"
   button controls the sensitivity of the "dbproxy" entry widget.
*/

RCVAR rc_vars[] = {
    { "gretldir", N_("Main gretl directory"), NULL, paths.gretldir, 
      MACHSET | BROWSER, MAXLEN, TAB_MAIN, NULL },
    { "userdir", N_("User's gretl directory"), NULL, paths.workdir, 
      INVISET, MAXLEN, TAB_MAIN, NULL },
    { "updater", N_("Tell me about gretl updates"), NULL, &updater, 
      BOOLSET, 0, TAB_MAIN, NULL },
#ifndef G_OS_WIN32
    { "winsize", N_("Remember main window size"), NULL, &winsize, 
      BOOLSET, 0, TAB_MAIN, NULL },
#endif
#ifdef ENABLE_NLS
    { "lcnumeric", N_("Use locale setting for decimal point"), NULL, &lcnumeric, 
      BOOLSET, 0, TAB_MAIN, NULL },
#endif
#ifdef G_OS_WIN32 	 
    { "wimp", N_("Emulate Windows look"), NULL, &use_wimp, 	 
      BOOLSET, 0, TAB_MAIN, NULL }, 	 
#endif
#if !defined(G_OS_WIN32) && !defined(OSX_BUILD)
    { "browser", N_("Web browser"), NULL, Browser, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
    { "shellok", N_("Allow shell commands"), NULL, &shellok, 
      BOOLSET, 0, TAB_MAIN, NULL },
    { "usecwd", N_("Set working directory from shell"), NULL, &usecwd, 
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
#ifdef ENABLE_NLS
    { "langpref", N_("Language preference"), NULL, &langpref, 
      LISTSET | INTSET, 0, TAB_MAIN, NULL },
#endif
#ifndef G_OS_WIN32 
    { "gnuplot", N_("Command to launch gnuplot"), NULL, paths.gnuplot, 
      MACHSET | BROWSER, MAXLEN, TAB_PROGS, NULL },
#endif
    { "Rcommand", N_("Command to launch GNU R"), NULL, Rcommand, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
    { "latex", N_("Command to compile TeX files"), NULL, latex, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
    { "viewdvi", N_("Command to view DVI files"), NULL, viewdvi, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#if !defined(G_OS_WIN32) && !defined(OSX_BUILD)
    { "viewps", N_("Command to view postscript files"), NULL, viewps, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
    { "viewpdf", N_("Command to view PDF files"), NULL, viewpdf, 
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
#if defined(HAVE_AUDIO) && !defined(G_OS_WIN32)
    { "midiplayer", N_("Program to play MIDI files"), NULL, midiplayer, 
      USERSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
    { "calculator", N_("Calculator"), NULL, calculator, 
      USERSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#ifdef HAVE_X12A
    { "x12a", N_("path to x12arima"), NULL, paths.x12a, 
      ROOTSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
#ifdef HAVE_TRAMO
    { "tramo", N_("path to tramo"), NULL, paths.tramo, 
      ROOTSET | BROWSER, MAXSTR, TAB_PROGS, NULL},
#endif
    { "binbase", N_("gretl database directory"), NULL, paths.binbase, 
      USERSET | BROWSER, MAXLEN, TAB_DBS, NULL },
    { "ratsbase", N_("RATS data directory"), NULL, paths.ratsbase, 
      USERSET | BROWSER, MAXLEN, TAB_DBS, NULL },
    { "dbhost", N_("Database server name"), NULL, paths.dbhost, 
      USERSET, 32, TAB_DBS, NULL },
    { "dbproxy", N_("HTTP proxy (ipnumber:port)"), NULL, dbproxy, 
      USERSET, 21, TAB_DBS, NULL },
    { "useproxy", N_("Use HTTP proxy"), NULL, &use_proxy, 
      BOOLSET, 1, TAB_DBS, NULL },
    { "Fixed_font", N_("Fixed font"), NULL, fixedfontname, 
      USERSET, MAXLEN, TAB_NONE, NULL },
#if !defined(USE_GNOME)
    { "App_font", N_("Menu font"), NULL, appfontname, 
      USERSET, MAXLEN, TAB_NONE, NULL },
#endif
    { "DataPage", "Default data page", NULL, datapage, 
      INVISET, sizeof datapage, TAB_NONE, NULL },
    { "ScriptPage", "Default script page", NULL, scriptpage, 
      INVISET, sizeof scriptpage, TAB_NONE, NULL },    
    { "Png_font", N_("PNG graph font"), NULL, paths.pngfont, 
      INVISET, 32, TAB_NONE, NULL },
    { "Gp_colors", N_("Gnuplot colors"), NULL, gpcolors, 
      INVISET, sizeof gpcolors, TAB_NONE, NULL },
    { "tabwidth", "spaces per tab", NULL, &tabwidth, 
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "smarttab", "\"Smart\" Tab and Enter", NULL, &smarttab, 
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
    { "main_width", "main window width", NULL, &mainwin_width, 
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_height", "main window height", NULL, &mainwin_height, 
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_x", "main window x position", NULL, &main_x, 
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_y", "main window y position", NULL, &main_y, 
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "HC_by_default", N_("Use robust covariance matrix by default"), NULL,
      &hc_by_default, BOOLSET, 0, TAB_VCV, NULL },
    { "HC_xsect", N_("For cross-sectional data"), NULL, hc_xsect, 
      LISTSET, 5, TAB_VCV, NULL },
    { "HC_tseri", N_("For time-series data"), NULL, hc_tseri, 
      LISTSET, 5, TAB_VCV, NULL },
    { "HC_panel", N_("For panel data"), NULL, hc_panel, 
      LISTSET, 9, TAB_VCV, NULL },
    { "HC_garch", N_("For GARCH estimation"), NULL, hc_garch, 
      LISTSET, 5, TAB_VCV, NULL },
    { "manpref", N_("PDF manual preference"), NULL, &manpref, 
      RADIOSET | INTSET, 0, TAB_MAN, NULL },
    { NULL, NULL, NULL, NULL, 0, 0, TAB_NONE, NULL }
};

/* accessor functions */

int using_hc_by_default (void)
{
    return hc_by_default;
}

int get_manpref (void)
{
    return manpref;
}

void set_datapage (const char *str)
{
    strcpy(datapage, str);
}

void set_scriptpage (const char *str)
{
    strcpy(scriptpage, str);
}

const char *get_datapage (void)
{
    return datapage;
}

const char *get_scriptpage (void)
{
    return scriptpage;
}

static gretlopt set_paths_opt = OPT_X;

void force_english_help (void)
{
    set_paths_opt |= OPT_N;
    gretl_set_paths(&paths, set_paths_opt);
}

void set_fixed_font (void)
{
    if (fixed_font != NULL) {
	pango_font_description_free(fixed_font);
    }

    fixed_font = pango_font_description_from_string(fixedfontname);
}

#ifndef G_OS_WIN32

/* remedial setting of font path for libgd, in case it's not set right
   (which, I'm sorry to say, seems often to be the case on Linux)
*/

static void set_gd_fontpath (void)
{
    char *gdpath = NULL;
    char *newpath = NULL;

    if (gnuplot_has_ttf(0)) {
	/* we're OK, don't mess */
	return;
    }

    gdpath = getenv("GDFONTPATH");

    if (gdpath != NULL) {
	if (strstr(gdpath, "gretl") == NULL &&
	    strstr(gdpath, "GRETL") == NULL) {
	    newpath = g_strdup_printf("%s:%sfonts", gdpath, 
				      paths.gretldir);
	}
    } else {
	newpath = g_strdup_printf("%sfonts", paths.gretldir);
    }

    if (newpath != NULL) {
	setenv("GDFONTPATH", newpath, 1);
	g_free(newpath);
	gnuplot_has_ttf(1);
    }
}

static void record_shell_opt (void)
{
    char shellstamp[FILENAME_MAX];

    sprintf(shellstamp, "%s.gretl_shell_stamp", paths.dotdir);

    if (shellok) {
	FILE *fp = gretl_fopen(shellstamp, "w");

	if (fp != NULL) {
	    fputs("ok\n", fp);
	    fclose(fp);
	}
    } else {
	gretl_remove(shellstamp);
    }
}

#endif

#ifndef USE_GNOME

const char *get_app_fontname (void)
{
    return appfontname;
}

void set_app_font (const char *fontname)
{
    GtkSettings *settings;

    if (fontname != NULL && *fontname == 0) return;

    settings = gtk_settings_get_default();

    if (fontname == NULL) {
	g_object_set(G_OBJECT(settings), "gtk-font-name", appfontname, NULL);
    } else {
	GtkWidget *w;
	PangoFontDescription *pfd;
	PangoContext *pc;
	PangoFont *pfont;

	w = gtk_label_new("");
	pfd = pango_font_description_from_string(fontname);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);

	if (pfont != NULL) {
	    strcpy(appfontname, fontname);
	    g_object_set(G_OBJECT(settings), "gtk-font-name", appfontname, NULL);
	}

	gtk_widget_destroy(w);
	pango_font_description_free(pfd);
    }
}

#endif

static void slash_terminate (char *path)
{
    if (path == NULL || *path == '\0') return;

    if (path[strlen(path) - 1] != SLASH) {
	strcat(path, SLASHSTR);
    }
}

static void get_functions_dir (char *dirname)
{
    char *target = NULL;
    FILE *fp;
    int err, ok = 0;

    sprintf(dirname, "%sfunctions", paths.gretldir);
    err = gretl_mkdir(dirname);
    if (!err) {
	target = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest");
	if (target != NULL) {
	    fp = gretl_fopen(target, "w");
	    if (fp != NULL) {
		ok = 1;
		fclose(fp);
		gretl_remove(target);
	    }
	    free(target);
	}
    } 
    
    if (ok) return;

    sprintf(dirname, "%sfunctions", paths.dotdir); /* ?? */
    err = gretl_mkdir(dirname);
    if (!err) {
	target = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest");
	if (target != NULL) {
	    fp = gretl_fopen(target, "w");
	    if (fp != NULL) {
		ok = 1;
		fclose(fp);
		gretl_remove(target);
	    }
	    free(target);
	}
    }

    if (ok) return;

    *dirname = '\0';
}

static char startdir[MAXLEN];

void set_gretl_startdir (void)
{
    char *test = getenv("GRETL_STARTDIR");

    if (test != NULL) {
	*startdir = '\0';
	strncat(startdir, test, MAXLEN - 1);
    } else {
	test = getcwd(startdir, MAXLEN);
	if (test == NULL) {
	    *startdir = '\0';
	}
    }

#if 0
    fprintf(stderr, "startdir = '%s'\n", startdir);
#endif

    if (usecwd && *startdir != '\0') {
	int err = set_gretl_work_dir(startdir, &paths);

	if (err) {
	    fprintf(stderr, "%s\n", gretl_errmsg_get());
	} else {
	    fprintf(stderr, "working dir = '%s'\n", paths.workdir);
	}
    }
}

const char *get_gretl_startdir (void)
{
    return startdir;
}

void get_default_dir (char *s, int action)
{
    *s = '\0';

    if (action == SAVE_FUNCTIONS) {
	get_functions_dir(s);
    } else if (action == OPEN_RATS_DB) {
	strcpy(s, paths.ratsbase);
    } else {
	strcpy(s, paths.workdir);
    }

    slash_terminate(s);
}

#ifdef G_OS_WIN32
static const char *get_reg_base (const char *key)
{
    if (!strncmp(key, "x12a", 4)) {
        return "x12arima";
    } else if (!strncmp(key, "tramo", 5)) {
        return "tramo";
    } else {
	return "gretl";
    }
}
#endif

#ifdef OSX_BUILD
static int alt_ok (const char *prog)
{ 
    char *p, test[MAXSTR];
    int tr = strstr(prog, "tramo") != NULL;
    int ok;
    
    strcpy(test, paths.gretldir);
    p = strstr(test, "share/gretl");
    if (p != NULL) {
    	*p = 0;
    }
    
    if (tr) {
         strcat(test, "tramo/tramo");
    } else {
         strcat(test, "x12arima/x12a");
    }

    ok = check_for_prog(test);
    
    if (ok) {
        if (tr) {
            strcpy(paths.tramo, test);
	} else {
	    strcpy(paths.x12a, test);
	}
    }
    
    return ok;
}
#endif

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)

static void set_tramo_x12a_status (void)
{
    int ok;

#ifdef HAVE_TRAMO
    ok = 0;
    if (*paths.tramodir != '\0') {
	ok = check_for_prog(paths.tramo);
# ifdef OSX_BUILD
	if (!ok) {
	    ok = alt_ok(paths.tramo);
	}
# endif    
    }

    if (mdata != NULL) {
	flip(mdata->ui, "/Variable/TRAMO analysis", ok);
    }
#endif /* TRAMO */

#ifdef HAVE_X12A
    ok = 0;
    if (*paths.x12adir != '\0') {
	ok = check_for_prog(paths.x12a);
# ifdef OSX_BUILD    
	if (!ok) {
	    ok = alt_ok(paths.x12a);
	}
# endif  
    }

    if (mdata != NULL) {
	flip(mdata->ui, "/Variable/X-12-ARIMA analysis", ok);
    }    
#endif /* X12A */
}

#endif /* tramo || x12a */

#ifdef G_OS_WIN32

int check_for_prog (const char *prog)
{
    char tmp[MAXLEN];
    WIN32_FIND_DATA find_data;
    HANDLE hfind;
    int ret = 1;

    if (prog == NULL || *prog == '\0') {
	return 0;
    }

    hfind = FindFirstFile(prog, &find_data);
    if (hfind == INVALID_HANDLE_VALUE) {
	ret = 0;
    }
    FindClose(hfind);

    if (ret == 0) {
	char *p;

	ret = SearchPath(NULL, prog, NULL, MAXLEN, tmp, &p);
    }

    return ret;
}

#else

int is_executable (const char *s, uid_t myid, gid_t mygrp)
{
    struct stat buf;
    int ok = 0;

    if (stat(s, &buf) == 0 && (buf.st_mode & (S_IFREG|S_IFLNK))) {
	if (buf.st_uid == myid && (buf.st_mode & S_IXUSR)) {
	    ok = 1;
	} else if (buf.st_gid == mygrp && (buf.st_mode & S_IXGRP)) {
	    ok = 1;
	} else if (buf.st_uid != myid && buf.st_gid != mygrp &&
		   (buf.st_mode & S_IXOTH)) {
	    ok = 1;
	}
    }

    return ok;
}

int check_for_prog (const char *prog)
{
    uid_t myid = getuid();
    gid_t mygrp = getgid();

    char *path;
    char *pathcpy;
    char **dirs;
    char *fullpath;
    char *p;

    int max_dlen = 0;
    int found = 0;
    int i, ndirs;

    if (prog != NULL && *prog == '/') {
	return is_executable(prog, myid, mygrp);
    }

    path = getenv("PATH");
    if (path == NULL || *path == '\0') {
	return 0;
    }

    pathcpy = gretl_strdup(path);
    if (pathcpy == NULL) {
	return 0;
    }

    ndirs = 1;
    p = pathcpy;
    while (*p) {
	if (*p == ':') ndirs++;
	p++;
    }

    dirs = malloc(ndirs * sizeof *dirs);
    if (dirs == NULL) {
	free(pathcpy);
	return 0;
    }

    if (ndirs == 1) {
	dirs[0] = pathcpy;
	max_dlen = strlen(pathcpy);
    } else {
	for (i=0; i<ndirs; i++) {
	    int dlen;

	    dirs[i] = strtok((i == 0)? pathcpy : NULL, ":");
	    if (dirs[i] == NULL) {
		ndirs = i;
		break;
	    }
	    dlen = strlen(dirs[i]);
	    if (dlen > max_dlen) {
		max_dlen = dlen;
	    }
	}
    }

    if (ndirs == 0 || 
	(fullpath = malloc(max_dlen + strlen(prog) + 2)) == NULL) {
	free(dirs);
	free(pathcpy);
	return 0;
    }

    for (i=0; i<ndirs && !found; i++) { 
	sprintf(fullpath, "%s/%s", dirs[i], prog);
	found = is_executable(fullpath, myid, mygrp);
    }

    free(dirs);
    free(pathcpy);
    free(fullpath);

    return found;
}

static void root_check (void)
{
    if (getuid() == 0) {
	int resp;

	resp = yes_no_dialog ("gretl", _("You seem to be running gretl " 
			      "as root.  Do you really want to do this?"), 
			      0);
	if (resp == GRETL_NO) {
	    exit(EXIT_FAILURE);
	}
    }
}

void gretl_config_init (void)
{
    char *custprof = getenv("GRETL_PROFILE");

    if (custprof == NULL) {
	sprintf(rcfile, "%s/.gretl2rc", getenv("HOME"));
    } else {
	strcpy(rcfile, custprof);
    }

    read_rc();
    set_gretl_startdir();
    set_gd_fontpath();
    root_check();
}

#endif /* *nix versus Windows */

static void highlight_options_entry (const char *vname)
{
    GtkWidget *w;
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (!strcmp(vname, rc_vars[i].key)) {
	    w = rc_vars[i].widget;
	    if (w != NULL && GTK_IS_ENTRY(w)) {
		gtk_editable_select_region(GTK_EDITABLE(w), 0, -1);
		gtk_widget_grab_focus(w);
	    }
	    break;
	}
    }
}

static void option_dialog_canceled (GtkWidget *w, int *c)
{
    *c = 1;
}

int options_dialog (int page, const char *varname) 
{
    static GtkWidget *dialog;
    GtkWidget *notebook;
    GtkWidget *button;
    GtkWidget *hbox;
    int canceled = 0;

    if (dialog != NULL) {
	gtk_window_present(GTK_WINDOW(dialog));
	return 0;
    }

    dialog = gretl_dialog_new(_("gretl: options"), mdata->main, GRETL_DLG_BLOCK);
    gtk_dialog_set_has_separator(GTK_DIALOG(dialog), FALSE);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 2);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(gtk_widget_destroyed), 
		     &dialog);

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), notebook, 
		       TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    make_prefs_tab(notebook, TAB_MAIN);
    make_prefs_tab(notebook, TAB_DBS);
    make_prefs_tab(notebook, TAB_PROGS);
    make_prefs_tab(notebook, TAB_VCV);
    make_prefs_tab(notebook, TAB_MAN);

    hbox = GTK_DIALOG(dialog)->action_area;

    /* Apply button */
    button = apply_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_changes), NULL);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* Cancel button */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(option_dialog_canceled), 
		     &canceled);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_show(button);

    /* OK button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_changes), NULL);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);
    gtk_widget_show(button);

    if (page > 1 && page < TAB_MAX) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), page - 1);
    }

    if (varname != NULL) {
	highlight_options_entry(varname);
    }

    gtk_widget_show(dialog);

    return canceled;
}

static void flip_sensitive (GtkWidget *w, gpointer data)
{
    GtkWidget *entry = GTK_WIDGET(data);
    
    gtk_widget_set_sensitive(entry, GTK_TOGGLE_BUTTON(w)->active);
}

void set_path_callback (char *setvar, char *setting)
{
    int i = 0;

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].var == (void *) setvar) {
	    /* FIXME: utf-8 issues here? */
	    if (rc_vars[i].widget != NULL) {
		gtk_entry_set_text(GTK_ENTRY(rc_vars[i].widget), 
				   setting);
	    }
	    break;
	}
	i++;
    }
}

static void browse_button_callback (GtkWidget *w, RCVAR *rc)
{
    GtkWidget *top = g_object_get_data(G_OBJECT(w), "parent");
    int code = SET_PROG;

    if (strstr(rc->description, "directory") != NULL) {
	code = SET_DIR;
    }

    file_selector_with_parent(_(rc->description), code, 
			      FSEL_DATA_MISC, rc->var, top);
}

static GtkWidget *make_path_browse_button (RCVAR *rc, GtkWidget *w)
{
    GtkWidget *top = gtk_widget_get_toplevel(w);
    GtkWidget *b;

    b = gtk_button_new_with_label(_("Browse..."));
    g_object_set_data(G_OBJECT(b), "parent", top);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(browse_button_callback), 
		     rc);
    return b;
}

static gboolean takes_effect_on_restart (void)
{
    infobox(_("This change will take effect when you restart gretl"));
    return FALSE;
}

#define HIDE_SPANISH_MANUAL 1

static const char **get_radio_setting_strings (void *var, int *n)
{
    static const char *man_strs[] = {
        N_("English (US letter paper)"),
        N_("English (A4 paper)"),
        N_("Italian"),
	N_("Spanish")
    };
    const char **strs = NULL;

    *n = 0;

    if (var == &manpref) {
	strs = man_strs;
	*n = sizeof man_strs / sizeof man_strs[0];
#if HIDE_SPANISH_MANUAL
	*n -= 1;
#endif
    } 

    return strs;
}

static const char **get_list_setting_strings (void *var, int *n)
{
    static const char *hc_strs[] = {
	"HC0", "HC1", "HC2", "HC3", "HC3a", "HAC"
    };
    static const char *hc_panel_strs[] = {
	"Arellano", "PCSE"
    };
    static const char *garch_strs[] = {
	"QML", "BW"
    };
    const char **strs = NULL;

    *n = 0;

    if (var == hc_xsect || var == hc_tseri) {
	strs = hc_strs;
	*n = sizeof hc_strs / sizeof hc_strs[0];
	if (var == hc_xsect) *n -= 1;
    } else if (var == hc_panel) {
	strs = hc_panel_strs;
	*n = sizeof hc_panel_strs / sizeof hc_panel_strs[0];
    } else if (var == hc_garch) {
	strs = garch_strs;
	*n = sizeof garch_strs / sizeof garch_strs[0];
    } 

    return strs;
}

static void 
get_table_sizes (int page, int *n_str, int *n_bool, int *n_browse,
		 int *n_list)
{
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].tab == page) {
	    if (rc_vars[i].flags & BROWSER) {
		*n_browse += 1;
	    } else if (rc_vars[i].flags & LISTSET) {
		*n_list += 1;
	    }
	    if (rc_vars[i].flags & BOOLSET) {
		*n_bool += 1;
	    } else if (rc_vars[i].flags & RADIOSET) {
		*n_bool += 1;
	    } else if (!(rc_vars[i].flags & INVISET)) {
		*n_str += 1;
	    } 
	}
    }
}

static void radio_change_value (GtkWidget *w, int *v)
{
    if (GTK_TOGGLE_BUTTON(w)->active) {
	gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));

	*v = i;
    }
}

static void make_prefs_tab (GtkWidget *notebook, int tab) 
{
    GtkWidget *b_table = NULL, *s_table = NULL;
    GtkWidget *l_table = NULL;
    GtkWidget *box, *w = NULL;
    int s_len = 1, b_len = 0, l_len = 1;
    int s_cols, b_cols = 0, l_cols = 0;
    int b_col = 0;
    int n_str = 0;
    int n_bool = 0;
    int n_browse = 0;
    int n_list = 0;
    int langs = 0;
    RCVAR *rc;
    int i;
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    if (tab == TAB_MAIN) {
	w = gtk_label_new(_("General"));
    } else if (tab == TAB_DBS) {
	w = gtk_label_new(_("Databases"));
    } else if (tab == TAB_PROGS) {
	w = gtk_label_new(_("Programs"));
    } else if (tab == TAB_VCV) {
	w = gtk_label_new(_("HCCME"));
    } else if (tab == TAB_MAN) {
	w = gtk_label_new(_("Manuals"));
    }
    
    gtk_widget_show(w);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, w);   

    get_table_sizes(tab, &n_str, &n_bool, &n_browse, &n_list);

    s_cols = (n_browse > 0)? 3 : 2;

    if (tab == TAB_VCV && n_list > 0) {
	/* VCV tab -- put the list entries first, right aligned */
	l_cols = 2;
	l_table = gtk_table_new(l_len, l_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(l_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(l_table), 5);
	gtk_box_pack_start(GTK_BOX(box), l_table, FALSE, FALSE, 0);
	gtk_widget_show(l_table);
    }	

    if (n_str > 0) {
	s_table = gtk_table_new(s_len, s_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(s_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(s_table), 5);
	gtk_box_pack_start(GTK_BOX(box), s_table, FALSE, FALSE, 0);
	gtk_widget_show(s_table);
    }

    if (n_bool > 0) {
	b_cols = 2;
	b_table = gtk_table_new(1, b_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(b_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(b_table), 5);
	gtk_box_pack_start(GTK_BOX(box), b_table, FALSE, FALSE, 10);
	gtk_widget_show(b_table);
    }

    if (tab != TAB_VCV && n_list > 0) {
	/* non-VCV tab: language list comes last */
	l_cols = 2;
	langs = 1;
	l_table = gtk_table_new(l_len, l_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(l_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(l_table), 5);
	gtk_box_pack_start(GTK_BOX(box), l_table, FALSE, FALSE, 0);
	gtk_widget_show(l_table);
    }

    for (i=0; rc_vars[i].key != NULL; i++) {
	rc = &rc_vars[i];

	if (rc->tab != tab) {
	    /* the item is not on this page */
	    continue;
	}

	if ((rc->flags & BOOLSET) && rc->link == NULL) { 
	    /* simple boolean variable (check box) */
	    int rcval = *(int *) (rc->var);

	    rc->widget = gtk_check_button_new_with_label(_(rc->description));
	    gtk_table_attach_defaults(GTK_TABLE (b_table), rc->widget, 
				      b_col, b_col + 1, b_len, b_len + 1);

	    if (rcval) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), TRUE);
	    } else {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), FALSE);
	    }

	    /* special case: deferred processing */
	    if (!strcmp(rc->key, "lcnumeric") ||
		!strcmp(rc->key, "usecwd") ||
		!strcmp(rc->key, "wimp")) {
		g_signal_connect(G_OBJECT(rc->widget), "toggled",
				 G_CALLBACK(takes_effect_on_restart), 
				 NULL);
	    }

	    /* special case: link between toggle and preceding entry */
	    if (rc->len && !(rc->flags & FIXSET)) {
		gtk_widget_set_sensitive(rc_vars[i-1].widget,
					 GTK_TOGGLE_BUTTON(rc->widget)->active);
		g_signal_connect(G_OBJECT(rc->widget), "clicked",
				 G_CALLBACK(flip_sensitive),
				 rc_vars[i-1].widget);
	    } 

	    gtk_widget_show(rc->widget);
	    b_col++;

	    if (b_col == 2) {
		b_col = 0;
		b_len++;
		gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);
	    }

	    if (rc->flags & FIXSET) {
		gtk_widget_set_sensitive(rc->widget, FALSE);
		if (rc->len) {
		    gtk_widget_set_sensitive(rc_vars[i-1].widget, FALSE);
		}
	    }

	} else if (rc->flags & BOOLSET) { 
	    /* radio-button dichotomy */
	    int rcval = *(int *) (rc->var);
	    GtkWidget *button;
	    GSList *group = NULL;

	    /* do we have some padding to do? */
	    if (b_col == 1) {
		w = gtk_label_new("   ");
		gtk_table_attach_defaults(GTK_TABLE(b_table), w, 
					  b_col, b_col + 1, 
					  b_len, b_len + 1);
		b_col = 0;
		b_len++;
		gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);
	    }

	    b_len += 3;
	    gtk_table_resize(GTK_TABLE(b_table), b_len + 1, 2);

	    /* separator for the group? */
	    w = gtk_hseparator_new();
	    gtk_table_attach_defaults(GTK_TABLE(b_table), w, 
				      b_col, b_col + 1, 
				      b_len - 3, b_len - 2);  
	    gtk_widget_show(w);

	    /* then a first button */
	    button = gtk_radio_button_new_with_label(group, _(rc->link));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), button, 
				      b_col, b_col + 1, 
				      b_len - 2, b_len - 1);    
	    if (!rcval) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    }
	    gtk_widget_show(button);

	    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(button));

	    /* and a second button */
	    rc->widget = gtk_radio_button_new_with_label(group, 
							 _(rc->description));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), rc->widget, 
				      b_col, b_col + 1, 
				      b_len - 1, b_len);  
	    if (rcval) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rc->widget), TRUE);
	    }
	    gtk_widget_show(rc->widget);

	    if (rc->flags & FIXSET) {
		gtk_widget_set_sensitive(button, FALSE);
		gtk_widget_set_sensitive(rc->widget, FALSE);
	    }
	} else if (rc->flags & RADIOSET) {
	    int nopt, j, rcval = *(int *) (rc->var);
	    int rcol;
	    GtkWidget *b;
	    GSList *group = NULL;
	    const char **strs;

	    if (b_len > 0) {
		/* there are buttons above: add a separator and
		   make this section full-width */
		rcol = b_cols;
		b_len++;
		w = gtk_hseparator_new();
		gtk_table_attach_defaults(GTK_TABLE(b_table), w, 
					  0, rcol, b_len - 1, b_len);  
		gtk_widget_show(w);
	    } else {
		rcol = b_col + 1;
	    }

	    b_len++;
	    b = gtk_label_new(_(rc->description));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), b, 
				      b_col, rcol, 
				      b_len - 1, b_len);
	    gtk_widget_show(b);

	    strs = get_radio_setting_strings(rc->var, &nopt);

	    for (j=0; j<nopt; j++) {
		b_len++;
		gtk_table_resize(GTK_TABLE(b_table), b_len, 2);
		b = gtk_radio_button_new_with_label(group, _(strs[j]));
		gtk_table_attach_defaults(GTK_TABLE(b_table), b, 
					  b_col, rcol, 
					  b_len - 1, b_len);
		g_object_set_data(G_OBJECT(b), "action", 
				  GINT_TO_POINTER(j));
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), j == rcval);
		gtk_widget_show(b);
		group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
		g_signal_connect(G_OBJECT(b), "clicked",
				 G_CALLBACK(radio_change_value),
				 rc->var);
	    }
	} else if (rc->flags & LISTSET) {
	    int j, active = 0;

	    l_len++;

	    gtk_table_resize(GTK_TABLE(l_table), l_len, l_cols);
	    w = gtk_label_new(_(rc->description));
	    if (langs) {
		gtk_table_attach(GTK_TABLE(l_table), w, 
				 0, 1, l_len - 1, l_len,
				 0, 0, 0, 0);
	    } else {		
		gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
		gtk_table_attach_defaults(GTK_TABLE(l_table), 
					  w, 0, 1, l_len - 1, l_len);
	    } 
	    gtk_widget_show(w);

	    rc->widget = gtk_combo_box_new_text();
	    gtk_table_attach(GTK_TABLE(l_table), rc->widget, 
			     1, 2, l_len - 1, l_len,
			     0, 0, 0, 0);
	    if (langs) {
		active = *(int *) rc->var;
#ifdef ENABLE_NLS
		for (j=LANG_AUTO; j<LANG_MAX; j++) {
		    gtk_combo_box_append_text(GTK_COMBO_BOX(rc->widget), 
					      lang_string_from_id(j));
		}
#endif
	    } else {
		char *strvar = (char *) rc->var;
		const char **strs;
		int nopt;
		
		strs = get_list_setting_strings(rc->var, &nopt);
		for (j=0; j<nopt; j++) {
		    gtk_combo_box_append_text(GTK_COMBO_BOX(rc->widget), 
					      strs[j]);
		    if (!strcmp(strs[j], strvar)) {
			active = j;
		    }
		}
	    }
	    if (tab == TAB_VCV) {
		int ww = get_string_width("XXArellanoXXXXX");

		gtk_widget_set_size_request(rc->widget, ww, -1);
	    }
	    gtk_combo_box_set_active(GTK_COMBO_BOX(rc->widget), active);
	    gtk_widget_show(rc->widget);
	    if (langs) {
		g_signal_connect(G_OBJECT(rc->widget), "changed",
				 G_CALLBACK(takes_effect_on_restart), 
				 NULL);
	    }
	} else if (!(rc->flags & INVISET)) { 
	    /* visible string variable */
	    char *strvar = (char *) rc->var;

	    s_len++;

	    gtk_table_resize(GTK_TABLE(s_table), s_len, s_cols);
	    w = gtk_label_new(_(rc->description));
	    gtk_misc_set_alignment(GTK_MISC(w), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(s_table), 
				      w, 0, 1, s_len - 1, s_len);
	    gtk_widget_show(w);

	    rc->widget = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(s_table), 
				      rc->widget, 1, 2, s_len - 1, s_len);
	    gtk_entry_set_text(GTK_ENTRY(rc->widget), strvar);
	    gtk_widget_show(rc->widget);

	    if (rc->flags & BROWSER) {
		/* add path browse button */
		w = make_path_browse_button(rc, notebook);
		gtk_table_attach_defaults(GTK_TABLE(s_table), 
					  w, 2, 3, s_len - 1, s_len);
		gtk_widget_show(w);
	    }

	    if (rc->flags & FIXSET) {
		gtk_widget_set_sensitive(rc->widget, FALSE);
		gtk_widget_set_sensitive(w, FALSE);
	    }
	} 
    }

    if (tab == TAB_VCV) {
	/* we need a help button */
	GtkWidget *hb = gtk_hbox_new(FALSE, 0);

	w = gtk_label_new("");
	gtk_box_pack_start(GTK_BOX(hb), w, TRUE, TRUE, 0);
	gtk_widget_show(w);

	w = gtk_button_new_from_stock(GTK_STOCK_HELP);
	gtk_box_pack_start(GTK_BOX(hb), w, FALSE, FALSE, 0);
	g_signal_connect(G_OBJECT(w), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(HCCME));
	gtk_widget_show(w);

	gtk_box_pack_start(GTK_BOX(box), hb, FALSE, FALSE, 0);
	gtk_widget_show(hb);
    }
}

static void set_gp_colors (void)
{
    const char *s = gpcolors;
    char cstr[N_GP_COLORS][8];
    int i, nc = 0;

    for (i=0; i<N_GP_COLORS; i++) {
	if (sscanf(s, "%7s", cstr[i]) == 1) {
	    nc++;
	    s += 7;
	    if (*s == ' ') {
		s++;
	    } else {
		break;
	    }
	} else {
	    *cstr[i] = '\0';
	    break;
	}
    }

    if (nc == 4) {
	/* old-style */
	for (i=0; i<3; i++) {
	    set_graph_palette_from_string(i, cstr[i]);
	}
	set_graph_palette_from_string(BOXCOLOR, cstr[3]);
    } else {
	for (i=0; i<nc; i++) {
	    set_graph_palette_from_string(i, cstr[i]);
	}
    }
}

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)

static void maybe_revise_tramo_x12a_status (void)
{
    int doit = 0;

# ifdef HAVE_TRAMO
    if (strcmp(paths.tramo, gretl_tramo())) {
	doit = 1;
    } 
# endif

# ifdef HAVE_X12A
    if (strcmp(paths.x12a, gretl_x12_arima())) {
	doit = 1;
    }
# endif

    if (doit) {
	set_tramo_x12a_status();
    }
}

#endif

/* register and react to changes from Preferences dialog */

static void apply_changes (GtkWidget *widget, gpointer data) 
{
    const gchar *str;
    char *strvar;
    int *intvar;
    int i = 0;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].widget != NULL) {
	    if (rc_vars[i].flags & BOOLSET) {
		intvar = (int *) rc_vars[i].var;
		if (GTK_TOGGLE_BUTTON(rc_vars[i].widget)->active) {
		    *intvar = TRUE;
		} else {
		    *intvar = FALSE;
		}
	    } else if ((rc_vars[i].flags & USERSET) || 
		       (rc_vars[i].flags & ROOTSET) ||
		       (rc_vars[i].flags & MACHSET)) {
		str = gtk_entry_get_text(GTK_ENTRY(rc_vars[i].widget));
		if (str != NULL && *str != '\0') { 
		    strvar = (char *) rc_vars[i].var;
		    *strvar = '\0';
		    strncat(strvar, str, rc_vars[i].len - 1);
		}
	    } else if (rc_vars[i].flags & LISTSET) {
		GtkComboBox *box = GTK_COMBO_BOX(rc_vars[i].widget);

		if (rc_vars[i].flags & INTSET) {
		    int boxval = gtk_combo_box_get_active(box);
		    
		    intvar = (int *) rc_vars[i].var;
		    *intvar = boxval;
		} else {
		    gchar *boxval;

		    boxval = gtk_combo_box_get_active_text(box);
		    strvar = (char *) rc_vars[i].var;
		    *strvar = '\0';
		    strcat(strvar, boxval);
		    g_free(boxval);
		}
	    }
	}
    }

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    maybe_revise_tramo_x12a_status();
#endif

    write_rc(); /* note: calls gretl_set_paths */

    /* register these for session using libset apparatus */
    libset_set_bool(SHELL_OK, shellok);
    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_panel_hccme(hc_panel);
    set_garch_robust_vcv(hc_garch);

    gretl_www_init(paths.dbhost, dbproxy, use_proxy);
}

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *) b) {
	strcpy(s, "true");
    } else {
	strcpy(s, "false");
    }
}

#ifndef G_OS_WIN32

static int write_plain_text_rc (void) 
{
    FILE *rc;
    char val[6];
    char *strvar;
    int i;

    rc = gretl_fopen(rcfile, "w");
    if (rc == NULL) {
	file_write_errbox(rcfile);
	return E_FOPEN;
    }

    fprintf(rc, "# gretl config file\n");

    for (i=0; rc_vars[i].var != NULL; i++) {
	fprintf(rc, "# %s\n", rc_vars[i].description);
	if (rc_vars[i].flags & BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, val);
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, val);
	} else if (rc_vars[i].flags & INTSET) {
	    fprintf(rc, "%s = %d\n", rc_vars[i].key, *(int *) rc_vars[i].var);
	} else {
	    strvar = (char *) rc_vars[i].var;
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, strvar);
	}
    }

    rc_save_file_lists(rc);

    fclose(rc);

    return 0;
}

#endif /* !G_OS_WIN32 */

static void str_to_boolvar (const char *s, void *b)
{
    int *bvar = (int *) b;

    if (s == NULL) return;

    if (strcmp(s, "true") == 0 || strcmp(s, "1") == 0) {
	*bvar = TRUE;
    } else {
	*bvar = FALSE;
    }	
}

static void str_to_int (const char *s, void *b)
{
    int *ivar = (int *) b;

    if (s == NULL) return;

    if (sscanf(s, "%d", ivar) != 1) {
	*ivar = 0;
    }
}

static int dir_exists (const char *dname, FILE *fp)
{
    DIR *test;
    int ok = 0;

#ifdef G_OS_WIN32
    test = win32_opendir(dname);
#else
    test = opendir(dname);
#endif

    if (test != NULL) {
	ok = 1;
	if (fp != NULL) {
	    fprintf(fp, "Directory '%s' exists, OK\n", dname);
	}
	closedir(test);
    } else if (fp != NULL) {
	fprintf(fp, "Directory '%s' does not exist\n", dname);
    }

    return ok;
}

#if !defined(G_OS_WIN32) && !defined(OSX_BUILD)

static void maybe_fix_viewpdf (void)
{
    gchar *prog = NULL;
    const gchar *viewers[] = {
	"xpdf",
	"acroread",
	"evince",
	"kpdf",
	"gpdf",
	NULL
    };
    int i;

    prog = g_find_program_in_path(viewpdf);

    if (prog != NULL) {
	g_free(prog);
    } else {
	for (i=0; viewers[i] != NULL; i++) {
	    if (strcmp(viewers[i], viewpdf)) {
		prog = g_find_program_in_path(viewers[i]);
		if (prog != NULL) {
		    strcpy(viewpdf, viewers[i]);
		    g_free(prog);
		    break;
		}
	    }
	}
    }
}

#endif

static int common_read_rc_setup (void)
{
    int err = 0;

    libset_set_bool(SHELL_OK, shellok);
    libset_set_bool(USE_CWD, usecwd);
    set_gp_colors();

    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_panel_hccme(hc_panel);
    set_garch_robust_vcv(hc_garch);

    err = gretl_set_paths(&paths, set_paths_opt);

    gretl_www_init(paths.dbhost, dbproxy, use_proxy);
    set_tex_use_pdf(latex);

#if !defined(G_OS_WIN32) && !defined(OSX_BUILD)
    maybe_fix_viewpdf();
#endif

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_status();
# endif

# ifdef ENABLE_NLS
    set_lcnumeric(langpref, lcnumeric);
    if (langpref > 0) {
	force_language(langpref);
	if (langpref == LANG_C) {
	    force_english_help();
	}
    } 
# endif

    return err;
}

#ifdef G_OS_WIN32

void write_rc (void) 
{
    char bval[6];
    char ival[16];
    const char *strval;
    int i = 0, err = 0;

    for (i=0; rc_vars[i].key != NULL; i++) {

	if (rc_vars[i].flags & (FIXSET | ROOTSET)) {
	    /* read-only variables */
	    continue;
	}

	if (rc_vars[i].flags & BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, bval);
	    err += write_reg_val(HKEY_CURRENT_USER, 
				 "gretl", 
				 rc_vars[i].key, 
				 bval);
	} else if (rc_vars[i].flags & INTSET) {
	    sprintf(ival, "%d", *(int *) rc_vars[i].var);
	    err += write_reg_val(HKEY_CURRENT_USER, 
				 "gretl", 
				 rc_vars[i].key, 
				 ival);	    
	} else if (rc_vars[i].flags & MACHSET) {
	    strval = (char *) rc_vars[i].var;
	    err += write_reg_val(HKEY_LOCAL_MACHINE, 
				 get_reg_base(rc_vars[i].key),
				 rc_vars[i].key, 
				 strval);
	} else {
	    strval = (char *) rc_vars[i].var;
	    err += write_reg_val(HKEY_CURRENT_USER, 
				 get_reg_base(rc_vars[i].key),
				 rc_vars[i].key, 
				 strval);
	}
    }

    save_file_lists();
    gretl_set_paths(&paths, set_paths_opt);
}

static int get_network_settings (void)
{
    const char *inifile;
    FILE *fp;
    char *strvar;
    int gotini = 0;

    inifile = get_network_cfg_filename();

    if (*inifile && (fp = gretl_fopen(inifile, "r"))) {
	int j, calldrive = tolower(inifile[0]);
	char line[MAXLEN], key[32], linevar[MAXLEN];

	while (fgets(line, MAXLEN, fp)) {
	    int gotvar = 0;
	    char *p = line;

	    while (isspace(*p)) p++;
	    if (*p == '#') continue;

	    if (sscanf(p, "%31s", key) == 1) {
		strcpy(linevar, p + strlen(key) + 3); 
		chopstr(linevar); 
		gotvar = 0;
		for (j=0; rc_vars[j].key != NULL; j++) {
		    if (!strcmp(key, rc_vars[j].key)) {
			if (rc_vars[j].flags & BOOLSET) {
			    str_to_boolvar(linevar, rc_vars[j].var);
			} else if (rc_vars[j].flags & INTSET) {
			    str_to_int(linevar, rc_vars[j].var);
			} else {
			    if (!strcmp(key, "gretldir") && 
				tolower(linevar[0]) != calldrive) {
				gotini = 0;
				goto network_quit;
			    } else {
				strvar = (char *) rc_vars[j].var;
				*strvar = '\0';
				strncat(strvar, linevar, rc_vars[j].len - 1);
			    }
			}
			rc_vars[j].flags |= FIXSET;
			gotvar = gotini = 1;
		    }
		    if (gotvar) break;
		}
	    }
	}
    network_quit:
	fclose(fp);
    }

    return gotini;
}

void read_rc (void) 
{
    char value[MAXSTR];
    char *strvar;
    int i;

    if (chinese_locale()) {
	strcpy(fixedfontname, "MS Gothic 10");
    }

    if (get_network_settings() && *paths.workdir != '\0') {
	int err = set_gretl_work_dir(paths.workdir, &paths);

	if (err) {
	    gui_errmsg(err);
	}
	for (i=0; rc_vars[i].key != NULL; i++) {
	    if (rc_vars[i].var == paths.tramodir ||
		rc_vars[i].var == paths.x12adir) {
		rc_vars[i].flags |= FIXSET;
	    }
	}
    } 

    for (i=0; rc_vars[i].key != NULL; i++) {
	int err = 0;

	if (rc_vars[i].flags & FIXSET) {
	    continue;
	}

	*value = '\0';

	if (rc_vars[i].flags & ROOTSET) {
	    err = read_reg_val_with_fallback(HKEY_LOCAL_MACHINE, 
					     HKEY_CLASSES_ROOT,
					     get_reg_base(rc_vars[i].key),
					     rc_vars[i].key, 
					     value);
	} else if (rc_vars[i].flags & MACHSET) {
	    err = read_reg_val(HKEY_LOCAL_MACHINE, 
			       get_reg_base(rc_vars[i].key),
			       rc_vars[i].key, 
			       value);
	} else {
	    err = read_reg_val(HKEY_CURRENT_USER, 
			       get_reg_base(rc_vars[i].key),
			       rc_vars[i].key, 
			       value);
	}
	    
	if (!err && *value != '\0') {
	    if (rc_vars[i].flags & BOOLSET) {
		str_to_boolvar(value, rc_vars[i].var);
	    } else if (rc_vars[i].flags & INTSET) {
		str_to_int(value, rc_vars[i].var);
	    } else {
		strvar = (char *) rc_vars[i].var;
		*strvar = '\0';
		strncat(strvar, value, rc_vars[i].len - 1);
	    }
	}
    }

    read_file_lists();
    common_read_rc_setup();
    set_fixed_font();
    set_app_font(NULL);
}

#else /* end of win32 version, now plain GTK */

void write_rc (void)
{
    int err;

    err = write_plain_text_rc();
    if (!err) {
	gretl_set_paths(&paths, set_paths_opt);
	record_shell_opt();
    }
}

static void find_and_write_var (const char *key, const char *val)
{
    char *strvar;
    int j;

    for (j=0; rc_vars[j].key != NULL; j++) {
	if (!strcmp(key, rc_vars[j].key)) {
	    if (rc_vars[j].flags & BOOLSET) {
		str_to_boolvar(val, rc_vars[j].var);
	    } else if (rc_vars[j].flags & INTSET) {
		str_to_int(val, rc_vars[j].var);
	    } else {
		strvar = (char *) rc_vars[j].var;
		*strvar = '\0';
		strncat(strvar, val, rc_vars[j].len - 1);
	    }
	    break;
	}
    }
}

static void read_rc (void) 
{
    FILE *fp;
    char line[MAXLEN], key[32], linevar[MAXLEN];
    int i;

    fp = fopen(rcfile, "r");
    if (fp == NULL) {
	fprintf(stderr, "Couldn't read %s\n", rcfile);
	return;
    }

    i = 0;
    while (rc_vars[i].var != NULL) {
	if (fgets(line, MAXLEN, fp) == NULL) {
	    break;
	}
	if (line[0] == '#') {
	    continue;
	}
	if (!strncmp(line, "recent", 6)) {
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    strcpy(linevar, line + strlen(key) + 3); 
	    chopstr(linevar); 
	    if (*linevar) {
		find_and_write_var(key, linevar);
	    }
	}
	i++;
    }

    read_file_lists(fp, line);
    fclose(fp);
    common_read_rc_setup();
}

#endif /* end of non-Windows versions of read_rc, write_rc */

static int fontsel_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    
    if (!strcmp(s, "MenuFont")) {
	return APP_FONT_SELECTION;
    } else {
	return FIXED_FONT_SELECTION;
    }
}

/* font selection: non-Windows, gtk-2.0 version first */

#ifndef G_OS_WIN32

static void font_selection_ok (GtkWidget *w, GtkFontselHackDialog *fs)
{
    guint which = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(fs), "which"));
    gchar *fontname;

    fontname = gtk_fontsel_hack_dialog_get_font_name(fs);

    if (*fontname == 0) {
	g_free(fontname);
	gtk_widget_destroy(GTK_WIDGET(fs));
	return;
    }

    if (which == FIXED_FONT_SELECTION) {
	strcpy(fixedfontname, fontname);
	set_fixed_font();
	write_rc();
    } 

# ifndef USE_GNOME /* gnome handles the app font */
    else if (which == APP_FONT_SELECTION) {
	set_app_font(fontname);
	write_rc();
    }
# endif

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fs));
}

void font_selector (GtkAction *action)
{
    static GtkWidget *fontsel = NULL;
    int which = fontsel_code(action);
    int filter = GTK_FONT_HACK_LATIN;
    char *title = NULL;
    const char *fontname = NULL;

# ifdef USE_GNOME
    if (which == APP_FONT_SELECTION) return; /* shouldn't happen */
# endif

    if (fontsel != NULL) {
	gtk_window_present(GTK_WINDOW(fontsel));
        return;
    }

    if (which == FIXED_FONT_SELECTION) {
	title = _("Font for gretl output windows");
	filter = GTK_FONT_HACK_LATIN_MONO;
	fontname = fixedfontname;
    } 

# ifndef USE_GNOME
    else if (which == APP_FONT_SELECTION) {
	title = _("Font for menus and labels");
	fontname = appfontname;
    }
# endif

    fontsel = gtk_fontsel_hack_dialog_new(title);
    gtk_fontsel_hack_dialog_set_filter
	(GTK_FONTSEL_HACK_DIALOG(fontsel), filter);
    gtk_fontsel_hack_dialog_set_font_name 
	(GTK_FONTSEL_HACK_DIALOG(fontsel), fontname); 
    g_object_set_data(G_OBJECT(fontsel), "which", GINT_TO_POINTER(which));

    gtk_window_set_position(GTK_WINDOW(fontsel), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &fontsel);
    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_main_quit),
		     NULL);
    g_signal_connect(G_OBJECT(gtk_fontsel_hack_dialog_ok_button(fontsel)),
		     "clicked", 
		     G_CALLBACK(font_selection_ok),
		     fontsel);
    g_signal_connect(G_OBJECT(gtk_fontsel_hack_dialog_cancel_button(fontsel)),
		     "clicked", 
		     G_CALLBACK(delete_widget),
		     fontsel);

    gtk_widget_show(fontsel);

    gtk_main();
}

#else /* end non-win32 font selection, start win32 */

void font_selector (GtkAction *action)
{
    int flag = fontsel_code(action);
    char fontname[128];

    if (flag == FIXED_FONT_SELECTION) {
	strcpy(fontname, fixedfontname);
    } else {
	strcpy(fontname, appfontname);
    }

    win32_font_selector(fontname, flag);

    if (*fontname != '\0') {
	if (flag == FIXED_FONT_SELECTION) {
	    strcpy(fixedfontname, fontname);
	    set_fixed_font();
	    write_rc();
	} else {
	    set_app_font(fontname);
	    write_rc();
	} 	    
    }
}

#endif /* end win32 font selection */

void update_persistent_graph_colors (void)
{
    print_palette_string(gpcolors);
}

void dump_rc (void) 
{
    char dumper[MAXLEN];
    FILE *fp;
    char *tmp;
    char val[6];
    int i;

    sprintf(dumper, "%sconfig-dump.txt", paths.workdir);

    fp = gretl_fopen(dumper, "w");
    if (fp == NULL) {
	file_write_errbox(dumper);
	return;
    }

    fprintf(fp, "gretl version %s\n", GRETL_VERSION);
    fprintf(fp, "built with GTK %d.%d.%d\n", GTK_MAJOR_VERSION, GTK_MINOR_VERSION,
	    GTK_MICRO_VERSION);

    tmp = getenv("HOME");
    if (tmp != NULL) {
	fprintf(fp, "HOME='%s'\n", tmp);
    } else {
	fputs("HOME not set\n", fp);
    }

    for (i=0; rc_vars[i].var != NULL; i++) {
	fprintf(fp, "# %s\n", rc_vars[i].description);
	if (rc_vars[i].flags & BOOLSET) {
	    boolvar_to_str(rc_vars[i].var, val);
	    fprintf(fp, "%s = %s\n", rc_vars[i].key, val);
	} else if (rc_vars[i].flags & INTSET) {
	    fprintf(fp, "%s = %d\n", rc_vars[i].key, *(int *) rc_vars[i].var);
	} else {
	    fprintf(fp, "%s = %s\n", rc_vars[i].key, (char *) rc_vars[i].var);
	}
    }

    dir_exists(paths.gretldir, fp);

    printf("Config info written to %s\n", dumper);

    fclose(fp);
}

int gui_set_working_dir (char *dirname)
{
    int err = set_gretl_work_dir(dirname, &paths);

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_WDIR, dirname);
    } else {
	mkfilelist(FILE_LIST_WDIR, dirname);
    }

    return err;
}

struct wdir_setter {
    GtkWidget *dialog;
    GtkWidget *wdir_combo;
    GtkWidget *r1, *r2;
};

/* callback from the file selector */

void set_working_dir_callback (GtkWidget *w, char *path)
{
    set_combo_box_default_text(GTK_COMBO_BOX(w), path);
}

static void wdir_browse_callback (GtkWidget *w, struct wdir_setter *wset)
{
    GtkWidget *combo = wset->wdir_combo;
    const char *s = N_("gretl working directory");

    file_selector_with_parent(_(s), SET_WDIR, FSEL_DATA_MISC, combo, 
			      wset->dialog);
}

static void free_fname (gchar *s, gpointer p)
{
    g_free(s);
}

static void 
add_wdir_content (GtkWidget *dialog, struct wdir_setter *wset) 
{
    GtkWidget *hbox, *vbox, *w = NULL;
    GSList *group = NULL;
    GList *list = NULL;
    char tmp[MAXLEN];

    vbox = GTK_DIALOG(dialog)->vbox;
    list = get_working_dir_list();

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Working directory:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5); 

    strcpy(tmp, paths.workdir); 
    trim_slash(tmp);
 
    /* combo + browse button for current working dir */
    w = gtk_combo_box_entry_new_text();
    gtk_container_add(GTK_CONTAINER(hbox), w);
    set_combo_box_strings_from_list(GTK_COMBO_BOX(w), list);
    set_combo_box_default_text(GTK_COMBO_BOX(w), tmp);
    gtk_entry_set_width_chars(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(w))), 32);
    wset->wdir_combo = w;
    w = gtk_button_new_with_label(_("Browse..."));
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(wdir_browse_callback), wset);
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, 0, 0, 5); 

    vbox_add_hsep(vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("On start-up, gretl should use:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, 0, 0, 5);

    /* radio 1 for "next time" */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_radio_button_new_with_label(group, 
					_("the directory selected above"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    wset->r1 = w;

    /* radio 2 for "next time" */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_radio_button_new_with_label(group, _("the current directory "
						 "as determined via the shell"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), usecwd);
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    wset->r2 = w;

    g_list_foreach(list, (GFunc) free_fname, NULL);
    g_list_free(list);
}

static void 
apply_wdir_changes (GtkWidget *w, struct wdir_setter *wset)
{ 
    char tmp[MAXLEN];
    gchar *str;
    int err;

    str = gtk_combo_box_get_active_text(GTK_COMBO_BOX(wset->wdir_combo));
    *tmp = '\0';
    if (str != NULL) {
	strncat(tmp, str, MAXLEN - 2);
	g_free(str);
    }

    err = set_gretl_work_dir(tmp, &paths);
    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_WDIR, tmp);
    } else {
	mkfilelist(FILE_LIST_WDIR, tmp);
    }

    usecwd = GTK_TOGGLE_BUTTON(wset->r2)->active;
}

void working_dir_dialog (void) 
{
    static GtkWidget *dialog;
    struct wdir_setter wset;
    GtkWidget *button;
    GtkWidget *hbox;

    if (dialog != NULL) {
	gtk_window_present(GTK_WINDOW(dialog));
	return;
    }

    dialog = gretl_dialog_new(_("gretl: working directory"), 
			      mdata->main, GRETL_DLG_BLOCK);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(gtk_widget_destroyed), 
		     &dialog);

    wset.dialog = dialog;
    add_wdir_content(dialog, &wset);

    hbox = GTK_DIALOG(dialog)->action_area;

    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);

    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_wdir_changes), &wset);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), 
		     dialog);

    context_help_button(hbox, WORKDIR);

    gtk_widget_show_all(dialog);
}
