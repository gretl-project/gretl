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
#include "toolbar.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "menustate.h"
#include "session.h"

#include "libset.h"
#include "version.h"
#include "texprint.h"

#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>

#ifdef USE_GNOME
# define USE_GCONF
# include <gconf/gconf-client.h>
#endif

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#else
# include <sys/stat.h>
# include <fcntl.h>
# include <errno.h>
# include "gtkfontselhack.h"
#endif

#if !defined(G_OS_WIN32) && !defined(USE_GNOME)
char rcfile[MAXLEN];
#endif

extern int want_toolbar;
extern char Rcommand[MAXSTR];

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
static int useqr;
static int shellok;
static int manpref;
char gpcolors[32];
static char datapage[24];
static char scriptpage[24];

static int hc_by_default;
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
char midiplayer[MAXSTR];
#endif

enum {
    TAB_NONE = 0,
    TAB_MAIN,
    TAB_DBS,
    TAB_PROGS,
    TAB_SAVE,
    TAB_VCV,
    TAB_MAN
};

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
			  INVISET "invisible" (user) string 
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
      USERSET | BROWSER, MAXLEN, TAB_MAIN, NULL },
    { "expert", N_("Expert mode (no warnings)"), NULL, &expert, 
      BOOLSET, 0, TAB_MAIN, NULL },
    { "updater", N_("Tell me about gretl updates"), NULL, &updater, 
      BOOLSET, 0, TAB_MAIN, NULL },
    { "toolbar", N_("Show gretl toolbar"), NULL, &want_toolbar, 
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
    { "gnuplot", N_("Command to launch gnuplot"), NULL, paths.gnuplot, 
      MACHSET | BROWSER, MAXLEN, TAB_PROGS, NULL },
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
    { "usecwd", N_("Use current working directory as default"), 
      N_("Use gretl user directory as default"), &usecwd, 
      BOOLSET, 0, TAB_SAVE, NULL },
    { "useqr", N_("Use QR decomposition"), N_("Use Cholesky decomposition"), &useqr, 
      BOOLSET, 0, TAB_MAIN, NULL },
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
	FILE *fp = fopen(shellstamp, "w");

	if (fp != NULL) {
	    fputs("ok\n", fp);
	    fclose(fp);
	}
    } else {
	remove(shellstamp);
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
		remove(target);
	    }
	    free(target);
	}
    } 
    
    if (ok) return;

    sprintf(dirname, "%sfunctions", paths.workdir);
    err = gretl_mkdir(dirname);
    if (!err) {
	target = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest");
	if (target != NULL) {
	    fp = gretl_fopen(target, "w");
	    if (fp != NULL) {
		ok = 1;
		fclose(fp);
		remove(target);
	    }
	    free(target);
	}
    }

    if (ok) return;

    *dirname = '\0';
}

static char startdir[MAXLEN];

void set_program_startdir (const char *callname)
{
    char *test = getcwd(startdir, MAXLEN);

    if (test == NULL) {
	*startdir = '\0';
    }

#if 0
    fprintf(stderr, "Starting in '%s'\n", startdir);
    fprintf(stderr, "Called as '%s'\n", callname);
#endif
}

void get_default_dir (char *s, int action)
{
    *s = '\0';

    if (action == SAVE_FUNCTIONS) {
	get_functions_dir(s);
	if (*s != '\0') {
	    slash_terminate(s);
	    return;
	}
    } else if (action == OPEN_RATS_DB) {
	strcpy(s, paths.ratsbase);
	return;
    }

    if (usecwd && action != SAVE_DBDATA) {
	if (*startdir != '\0') {
	    strcpy(s, startdir);
	} else {
	    const char *sdir = get_session_dirname();
	    char *test = getcwd(s, MAXLEN);

	    if (test == NULL || (*sdir != '\0' && strstr(s, sdir))) {
		strcpy(s, paths.workdir);
	    }
	} 
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
	flip(mdata->ifac, "/Variable/TRAMO analysis", ok);
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
	flip(mdata->ifac, "/Variable/X-12-ARIMA analysis", ok);
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
#ifndef USE_GCONF
    sprintf(rcfile, "%s/.gretl2rc", getenv("HOME"));
#endif

    read_rc();
    set_gd_fontpath();

    if (!expert) {
	root_check();
    }
}

#endif /* *nix versus Windows */

static void option_dialog_canceled (GtkWidget *w, int *c)
{
    *c = 1;
}

int options_dialog (int page) 
{
    static GtkWidget *dialog;

    GtkWidget *notebook;
    GtkWidget *button;
    GtkWidget *hbox;
    int canceled = 0;

    if (dialog != NULL) {
	gdk_window_raise(dialog->window);
	return 0;
    }

    dialog = gretl_dialog_new(_("gretl: options"), mdata->w, GRETL_DLG_BLOCK);
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
    make_prefs_tab(notebook, TAB_SAVE);
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

    if (page > 0) {
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), page);
    }

    gtk_widget_show(dialog);

    return canceled;
}

void options_dialog_callback (gpointer p, guint u, GtkWidget *w)
{
    options_dialog(u);
}

static void flip_sensitive (GtkWidget *w, gpointer data)
{
    GtkWidget *entry = GTK_WIDGET(data);
    
    gtk_widget_set_sensitive(entry, GTK_TOGGLE_BUTTON(w)->active);
}

void filesel_set_path_callback (const char *setting, char *strvar)
{
    int i = 0;

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].var == (void *) strvar) {
	    /* FIXME: utf-8 issues here?? */
	    gtk_entry_set_text(GTK_ENTRY(rc_vars[i].widget), 
			       setting);
	    break;
	}
	i++;
    }
}

static void browse_button_callback (GtkWidget *w, RCVAR *rc)
{
    int code = SET_PROG;

    if (strstr(rc->description, "directory") != NULL) {
	code = SET_DIR;
    }

    file_selector(_(rc->description), code, FSEL_DATA_MISC, rc->var);
}

static GtkWidget *make_path_browse_button (RCVAR *rc)
{
    GtkWidget *b;

    b = gtk_button_new_with_label(_("Browse..."));
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

static GList *get_settings_list (void *var, int *nopt)
{
    char *hc_strs[] = {
	"HC0", "HC1", "HC2", "HC3", "HC3a", "HAC"
    };
    char *hc_panel_strs[] = {
	"Arellano", "PCSE"
    };
    char *garch_strs[] = {
	"QML", "BW"
    };
    const char *man_strs[] = {
        N_("English (US letter paper)"),
        N_("English (A4 paper)"),
        N_("Italian"),
	N_("Spanish")
    };
    GList *list = NULL;
    int i, n = 0;

    if (var == hc_xsect || var == hc_tseri) {
	n = sizeof hc_strs / sizeof hc_strs[0];
	if (var == hc_xsect) n--;
	for (i=0; i<n; i++) {
	    list = g_list_append(list, hc_strs[i]);
	}
    } else if (var == hc_panel) {
	n = sizeof hc_panel_strs / sizeof hc_panel_strs[0];
	for (i=0; i<n; i++) {
	    list = g_list_append(list, hc_panel_strs[i]);
	}
    } else if (var == hc_garch) {
	n = sizeof garch_strs / sizeof garch_strs[0];
	for (i=0; i<n; i++) {
	    list = g_list_append(list, garch_strs[i]);
	}
    } else if (var == &manpref) {
	n = sizeof man_strs / sizeof man_strs[0];
#if HIDE_SPANISH_MANUAL
	n--;
#endif
	for (i=0; i<n; i++) {
	    list = g_list_append(list, _(man_strs[i]));
	}
    }

    if (nopt != NULL) {
	*nopt = n;
    }

    return list;
}

static void 
get_table_sizes (int page, int *n_str, int *n_bool, int *n_browse)
{
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].tab != page) {
	    continue;
	}

	if (rc_vars[i].flags & BROWSER) {
	    *n_browse += 1;
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
    GtkWidget *box, *w = NULL;
    int s_len = 1, b_len = 0, b_col = 0;
    int n_str = 0;
    int n_bool = 0;
    int n_browse = 0;
    int s_cols;
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
    } else if (tab == TAB_SAVE) {
	w = gtk_label_new(_("File Open/Save"));
    } else if (tab == TAB_VCV) {
	w = gtk_label_new(_("HCCME"));
    } else if (tab == TAB_MAN) {
	w = gtk_label_new(_("Manuals"));
    }
    
    gtk_widget_show(w);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, w);   

    get_table_sizes(tab, &n_str, &n_bool, &n_browse);

    s_cols = (n_browse > 0)? 3 : 2;

    if (n_str > 0) {
	s_table = gtk_table_new(s_len, s_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(s_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(s_table), 5);
	gtk_box_pack_start(GTK_BOX(box), s_table, FALSE, FALSE, 0);
	gtk_widget_show(s_table);
    }
    
    if (n_bool > 0) {
	b_table = gtk_table_new(1, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(b_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(b_table), 5);
	gtk_box_pack_start(GTK_BOX(box), b_table, FALSE, FALSE, 10);
	gtk_widget_show(b_table);
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

	    /* special case: warning */
	    if (!strcmp(rc->key, "lcnumeric") || !strcmp(rc->key, "wimp")) {
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
	} else if (rc->flags & LISTSET) {
	    char *strvar = (char *) rc->var;
	    GList *list;

	    s_len++;

	    gtk_table_resize(GTK_TABLE(s_table), s_len, s_cols);
	    w = gtk_label_new(_(rc->description));
	    gtk_misc_set_alignment(GTK_MISC(w), 0.75, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(s_table), 
				      w, 0, 1, s_len - 1, s_len);
	    gtk_widget_show(w);

	    rc->widget = gtk_combo_new();
	    gtk_table_attach(GTK_TABLE(s_table), rc->widget, 
			     1, 2, s_len-1, s_len,
			     0, 0, 0, 0);

	    list = get_settings_list(rc->var, NULL);
	    gtk_combo_set_popdown_strings(GTK_COMBO(rc->widget), list);
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(rc->widget)->entry), strvar);
	    gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(rc->widget)->entry), 8);
	    gtk_editable_set_editable(GTK_EDITABLE(GTK_COMBO(rc->widget)->entry), 
				      FALSE);
	    gtk_widget_show(rc->widget);
	} else if (rc->flags & RADIOSET) {
	    int i, rcval = *(int *) (rc->var);
	    int nopt = 0;
	    GtkWidget *b;
	    GSList *group = NULL;
	    GList *list, *mylist;

	    b_len++;
	    b = gtk_label_new(_(rc->description));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), b, 
				      b_col, b_col + 1, 
				      b_len - 1, b_len);
	    gtk_widget_show(b);

	    mylist = list = get_settings_list(rc->var, &nopt);

	    for (i=0; i<nopt; i++) {
		b_len++;
		gtk_table_resize(GTK_TABLE(b_table), b_len, 2);
		b = gtk_radio_button_new_with_label(group, mylist->data);
		gtk_table_attach_defaults(GTK_TABLE(b_table), b, 
					  b_col, b_col + 1, 
					  b_len - 1, b_len);
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), i == rcval);
		g_object_set_data(G_OBJECT(b), "action", 
				  GINT_TO_POINTER(i));
		g_signal_connect(G_OBJECT(b), "clicked",
				 G_CALLBACK(radio_change_value),
				 rc->var);
		gtk_widget_show(b);
		group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
		mylist = g_list_next(mylist);
	    }
	    g_list_free(list);
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
				      rc->widget, 1, 2, s_len-1, s_len);
	    gtk_entry_set_text(GTK_ENTRY(rc->widget), strvar);
	    gtk_widget_show(rc->widget);

	    /* path browse button */
	    if (rc->flags & BROWSER) {
		w = make_path_browse_button(rc);
		gtk_table_attach_defaults(GTK_TABLE(s_table), 
					  w, 2, 3, s_len-1, s_len);
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

#ifdef ENABLE_NLS

# ifdef G_OS_WIN32

struct langname {
    const char *abbr;
    const char *full;
};

static void set_lcnumeric (void)
{
    if (lcnumeric) {
	struct langname names[] = {
	    { "es", "Spanish" },
	    { "eu", "Basque" },
	    { "fr", "French" },
	    { "it", "Italian" },
	    { "pl", "Polish" },
	    { "de", "German" },
	    { "pt", "Portuguese" },
	    { NULL, NULL }
	};
	char *lang = getenv("LANG");
	char *set = NULL;
	int i;

	if (lang != NULL) {
	    for (i=0; names[i].abbr != NULL; i++) {
		if (!strncmp(lang, names[i].abbr, 2)) {
		    set = setlocale(LC_NUMERIC, names[i].full);
		    if (set == NULL) {
			set = setlocale(LC_NUMERIC, names[i].abbr);
		    }
		    if (set != NULL) {
			break;
		    }
		}
	    }
	}

	if (set == NULL) {
	    setlocale(LC_NUMERIC, "");
	    putenv("LC_NUMERIC=");
	}
    } else {
	setlocale(LC_NUMERIC, "C");
	putenv("LC_NUMERIC=C");
    }

    reset_local_decpoint();
}

# else /* ! G_OS_WIN32 */

static void set_lcnumeric (void)
{
    if (lcnumeric) {
	char *lang = getenv("LANG");
	char *set = NULL;

	if (lang != NULL) {
	    set = setlocale(LC_NUMERIC, lang);
	    fprintf(stderr, "setlocale(LC_NUMERIC, \"%s\") returned %s\n", 
		    lang, set);
	} else {
	    fprintf(stderr, "set_lcnumeric: getenv(\"LANG\") gave NULL\n");
	}
	if (set == NULL) {
	    setlocale(LC_NUMERIC, "");
	    putenv("LC_NUMERIC=");
	}
    } else {
	setlocale(LC_NUMERIC, "C");
	putenv("LC_NUMERIC=C");
    }

    reset_local_decpoint();
}

# endif

#endif /* ENABLE_NLS */

static void set_gp_colors (void)
{
    char cstr[N_GP_COLORS][8];
    int i, nc;

    *cstr[0] = *cstr[1] = *cstr[2] = *cstr[3] = '\0';

    nc = sscanf(gpcolors, "%7s %7s %7s %7s", 
		cstr[0], cstr[1], cstr[2], cstr[3]);

    for (i=0; i<nc; i++) {
	set_graph_palette_from_string(i, cstr[i]);
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
    int i = 0;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].widget != NULL) {
	    if (rc_vars[i].flags & BOOLSET) {
		if (GTK_TOGGLE_BUTTON(rc_vars[i].widget)->active) {
		    *(int *) (rc_vars[i].var) = TRUE;
		} else {
		    *(int *) (rc_vars[i].var) = FALSE;
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
		GtkWidget *entry = GTK_COMBO(rc_vars[i].widget)->entry;

		str = gtk_entry_get_text(GTK_ENTRY(entry));
		if (str != NULL && *str != '\0') { 
		    strvar = (char *) rc_vars[i].var;
		    *strvar = '\0';
		    strncat(strvar, str, rc_vars[i].len - 1);
		}
	    }
	}
    }

#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    maybe_revise_tramo_x12a_status();
#endif

    write_rc(); /* note: calls gretl_set_paths */

    show_or_hide_toolbar(want_toolbar);

    /* register these for session using libset apparatus */
    libset_set_bool(USE_QR, useqr);
    libset_set_bool(USE_CWD, usecwd);
    libset_set_bool(SHELL_OK, shellok);
    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_panel_hccme(hc_panel);
    set_garch_robust_vcv(hc_garch);

    gretl_www_init(paths.dbhost, dbproxy, use_proxy);
}

#ifndef USE_GCONF

static void str_to_boolvar (char *s, void *b)
{
    int *bvar = (int *) b;

    if (s == NULL) return;

    if (strcmp(s, "true") == 0 || strcmp(s, "1") == 0) {
	*bvar = TRUE;
    } else {
	*bvar = FALSE;
    }	
}

static void str_to_int (char *s, void *b)
{
    int *ivar = (int *) b;

    if (s == NULL) return;

    if (sscanf(s, "%d", ivar) != 1) {
	*ivar = 0;
    }
}

#endif

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *) b) {
	strcpy(s, "true");
    } else {
	strcpy(s, "false");
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

static int common_read_rc_setup (void)
{
    int err = 0;

    libset_set_bool(USE_QR, useqr);
    libset_set_bool(USE_CWD, usecwd);
    libset_set_bool(SHELL_OK, shellok);
    set_gp_colors();
    
    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_garch_robust_vcv(hc_garch);

    err = gretl_set_paths(&paths, set_paths_opt);

    gretl_www_init(paths.dbhost, dbproxy, use_proxy);
    set_tex_use_pdf(latex);

# if defined(HAVE_TRAMO) || defined(HAVE_X12A)
    set_tramo_x12a_status();
# endif

# ifdef ENABLE_NLS
    set_lcnumeric();
# endif

    return err;
}

/* next section: variant versions of write_rc and read_rc, depending
   on both gconf presence and platform
*/

#ifdef USE_GCONF 

void write_rc (void) 
{
    GConfClient *client;   
    char key[MAXSTR];
    gboolean bval;
    int ival;
    char *strvar;
    int i;

    client = gconf_client_get_default();

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(key, "/apps/gretl/%s", rc_vars[i].key);
	if (rc_vars[i].flags & BOOLSET) {
	    bval = *(gboolean *) rc_vars[i].var;
	    gconf_client_set_bool(client, key, bval, NULL);
	} else if (rc_vars[i].flags & INTSET) {
	    ival = *(int *) rc_vars[i].var;
	    gconf_client_set_int(client, key, ival, NULL);
	} else {
	    strvar = (char *) rc_vars[i].var;
	    gconf_client_set_string(client, key, strvar, NULL);
	}
    }

    save_file_lists(client);
    g_object_unref(G_OBJECT(client));
    gretl_set_paths(&paths, set_paths_opt);
    record_shell_opt();
}

static void read_rc (void) 
{
    GConfClient *client;
    GError *error = NULL;
    gboolean bval;
    int ival;
    gchar *strval;
    char key[MAXSTR];
    int i;

    client = gconf_client_get_default();

    for (i=0; rc_vars[i].key != NULL; i++) {
	sprintf(key, "/apps/gretl/%s", rc_vars[i].key);
	if (rc_vars[i].flags & BOOLSET) {
	    bval = gconf_client_get_bool(client, key, &error);
	    if (error) {
		fprintf(stderr, "Error reading %s\n", rc_vars[i].key);
		g_clear_error(&error);
	    } else {
		*(int *) rc_vars[i].var = (int) bval;
	    }
	} else if (rc_vars[i].flags & INTSET) {
	    ival = gconf_client_get_int(client, key, &error);
	    if (error) {
		fprintf(stderr, "Error reading %s\n", rc_vars[i].key);
		g_clear_error(&error);
	    } else {
		*(int *) rc_vars[i].var = ival;
	    }	    
	} else {
	    strval = gconf_client_get_string(client, key, &error);
	    if (error) {
		fprintf(stderr, "Error reading %s\n", rc_vars[i].key);
		g_clear_error(&error);
	    } else if (strval != NULL) {
		if (*strval != '\0') {
		    char *strvar = (char *) rc_vars[i].var;

		    *strvar = '\0';
		    strncat(strvar, strval, rc_vars[i].len - 1);
		}
		g_free(strval);
	    }
	}
    }

    read_file_lists(client);
    g_object_unref(G_OBJECT(client));
    common_read_rc_setup();
}

/* end of gconf version, now win32 */

#elif defined(G_OS_WIN32)

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

#if 0
    if (err) {
	win_show_last_error();
    }
#endif

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

    /* FIXME was paths.userdir below */

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

#else /* end of gconf and win32 versions, now plain GTK */

void write_rc (void) 
{
    FILE *rc;
    char val[6];
    char *strvar;
    int i;

    rc = fopen(rcfile, "w");
    if (rc == NULL) {
	file_write_errbox(rcfile);
	return;
    }

    fprintf(rc, "# gretl config file (note: not used by gnome version)\n");

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

    save_file_lists(rc);
    fclose(rc);
    gretl_set_paths(&paths, set_paths_opt);
    record_shell_opt();
}

static void read_rc (void) 
{
    FILE *fp;
    char line[MAXLEN], key[32], linevar[MAXLEN];
    char *strvar;
    int i, j;

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
	    for (j=0; rc_vars[j].key != NULL; j++) {
		if (!strcmp(key, rc_vars[j].key)) {
		    if (rc_vars[j].flags & BOOLSET) {
			str_to_boolvar(linevar, rc_vars[j].var);
		    } else if (rc_vars[j].flags & INTSET) {
			str_to_int(linevar, rc_vars[j].var);
		    } else {
			strvar = (char *) rc_vars[j].var;
			*strvar = '\0';
			strncat(strvar, linevar, rc_vars[j].len - 1);
		    }
		    break;
		}
	    }
	}
	i++;
    }

    read_file_lists(fp, line);
    fclose(fp);
    common_read_rc_setup();
}

#endif /* end of "plain gtk" versions of read_rc, write_rc */

/* font selection: non-Windows, gtk-2.0 version first */

#ifndef G_OS_WIN32

static void font_selection_ok (GtkWidget *w, GtkFontSelectionHackDialog *fs)
{
    guint which = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(fs), "which"));
    gchar *fontname;

    fontname = gtk_font_selection_hack_dialog_get_font_name(fs);

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

# ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	GtkWidget *fentry = g_object_get_data(G_OBJECT(fs), "font_entry");

	gtk_entry_set_text(GTK_ENTRY(fentry), fontname);
    }
# endif

# ifndef USE_GNOME /* gnome handles the app font */
    else if (which == APP_FONT_SELECTION) {
	set_app_font(fontname);
	write_rc();
    }
# endif

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fs));
}

static void fontsel_quit (GtkWidget *w, gpointer p)
{
    gtk_main_quit();
}

void font_selector (gpointer data, guint which, GtkWidget *widget)
{
    static GtkWidget *fontsel = NULL;
    int filter = GTK_FONT_HACK_LATIN;
    char *title = NULL;
    const char *fontname = NULL;

# ifdef USE_GNOME
    if (which == APP_FONT_SELECTION) return; /* shouldn't happen */
# endif

    if (fontsel != NULL) {
	if (!GTK_WIDGET_VISIBLE(fontsel)) gtk_widget_show (fontsel);
        gdk_window_raise(fontsel->window);
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
# ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	fontname = paths.pngfont;
    }
# endif

    fontsel = gtk_font_selection_hack_dialog_new(title);
    gtk_font_selection_hack_dialog_set_filter
	(GTK_FONT_SELECTION_HACK_DIALOG (fontsel), filter);
    gtk_font_selection_hack_dialog_set_font_name 
	(GTK_FONT_SELECTION_HACK_DIALOG (fontsel), fontname); 
    g_object_set_data(G_OBJECT(fontsel), "which", GINT_TO_POINTER(which));

# ifdef PLOT_FONT_SELECTOR
    if (which == GRAPH_FONT_SELECTION) {
	 g_object_set_data(G_OBJECT(fontsel), "font_entry", widget);
    }
# endif

    gtk_window_set_position (GTK_WINDOW (fontsel), GTK_WIN_POS_MOUSE);

    g_signal_connect (G_OBJECT(fontsel), "destroy",
		      G_CALLBACK(gtk_widget_destroyed),
		      &fontsel);
    g_signal_connect (G_OBJECT(fontsel), "destroy",
		      G_CALLBACK(fontsel_quit),
		      NULL);
    g_signal_connect (G_OBJECT(GTK_FONT_SELECTION_HACK_DIALOG(fontsel)->ok_button),
		      "clicked", 
		      G_CALLBACK(font_selection_ok),
		      fontsel);
    g_signal_connect(G_OBJECT(GTK_FONT_SELECTION_HACK_DIALOG(fontsel)->cancel_button),
		     "clicked", 
		     G_CALLBACK(delete_widget),
		     fontsel);

    gtk_widget_show (fontsel);

    gtk_main();
}

#else /* end non-win32 font selection, start win32 */

static const char *font_weight_string (int weight)
{
    if (weight >= FW_THIN && weight <= FW_LIGHT) {
	return " Light";
    }
    if (weight >= FW_NORMAL && weight <= FW_DEMIBOLD) {
	return "";
    }
    if (weight >= FW_BOLD) {
	return " Bold";
    }
    return "";
}

static void fontname_to_win32 (const char *src, int fixed,
			       char *name, int *pts)
{
    int sz;

    if (sscanf(src, "%31[^0123456789]%d", name, &sz) == 2) {
	size_t i, n = strlen(name);

	for (i=n-1; i>0; i--) {
	    if (name[i] == ' ') name[i] = 0;
	    else break;
	}
	*pts = sz * 10; /* measured in tenths of a point */
    } else {
	*name = 0;
	strncat(name, src, 31);
	if (fixed) *pts = 100;
	else *pts = 80;
    }
}

void font_selector (gpointer data, guint which, GtkWidget *widget)
{
    CHOOSEFONT cf;            /* common dialog box structure */
    LOGFONT lf;               /* logical font structure */
    char fontname[48];

    ZeroMemory(&cf, sizeof cf);
    cf.lStructSize = sizeof cf;
    cf.Flags = CF_SCREENFONTS | CF_TTONLY | CF_LIMITSIZE | CF_INITTOLOGFONTSTRUCT;
    cf.nSizeMin = 6;
    cf.nSizeMax = 24;

    ZeroMemory(&lf, sizeof lf);

    cf.Flags |= CF_NOSCRIPTSEL;
    
    lf.lfWeight = FW_REGULAR;
    lf.lfCharSet = DEFAULT_CHARSET;

    if (which == FIXED_FONT_SELECTION) {
	cf.Flags |= CF_FIXEDPITCHONLY;
	fontname_to_win32(fixedfontname, 1, lf.lfFaceName, &(cf.iPointSize));
    } 
    else if (which == APP_FONT_SELECTION) {
	fontname_to_win32(appfontname, 0, lf.lfFaceName, &(cf.iPointSize));
    } 
# ifdef PLOT_FONT_SELECTOR
    else if (which == GRAPH_FONT_SELECTION) {
	fontname_to_win32(paths.pngfont, 0, lf.lfFaceName, &(cf.iPointSize));
    }
# endif

    cf.lpLogFont = &lf;

    if (ChooseFont(&cf) == TRUE && *(cf.lpLogFont->lfFaceName)) {
	sprintf(fontname, "%s%s %d", cf.lpLogFont->lfFaceName, 
		font_weight_string(cf.lpLogFont->lfWeight), 
		cf.iPointSize / 10);
	if (which == FIXED_FONT_SELECTION) {
	    strcpy(fixedfontname, fontname);
	    set_fixed_font();
	    write_rc();
	} 
	else if (which == APP_FONT_SELECTION) {
	    set_app_font(fontname);
	    write_rc();
	}
# ifdef PLOT_FONT_SELECTOR
	else if (which == GRAPH_FONT_SELECTION) {
	    gtk_entry_set_text(GTK_ENTRY(widget), fontname);
	}
# endif
    }
}

#endif /* end win32 font selection */

/* graph color selection apparatus */

#define XPMROWS 19
#define XPMCOLS 17

#define scale_round(v) ((v) * 255.0 / 65535.0)

static GtkWidget *get_image_for_color (const gretlRGB *color)
{
    static char **xpm = NULL;
    GdkPixbuf *icon;
    GtkWidget *image;
    char colstr[8] = {0};
    int i;

    if (color == NULL) {
	return NULL;
    }

    if (xpm == NULL) {
	xpm = strings_array_new_with_length(XPMROWS, XPMCOLS);
	if (xpm != NULL) {
	    for (i=0; i<XPMROWS; i++) {
		if (i == 0) {
		    strcpy(xpm[i], "16 16 2 1");
		} else if (i == 1) {
		    strcpy(xpm[i], "X      c #000000");
		} else if (i == 2) {
		    strcpy(xpm[i], ".      c #000000");
		} else if (i == 3 || i == XPMROWS - 1) {
		    strcpy(xpm[i], "................");
		} else {
		    strcpy(xpm[i], ".XXXXXXXXXXXXXX.");
		}
	    }
	}
    }

    if (xpm == NULL) {
	return NULL;
    }

    print_rgb_hash(colstr, color);

    for (i=0; i<6; i++) {
	xpm[1][10+i] = colstr[i+1];
    }    

    icon = gdk_pixbuf_new_from_xpm_data((const char **) xpm);
    image = gtk_image_new_from_pixbuf(icon);
    
    return image;
}

static void color_select_callback (GtkWidget *button, GtkWidget *w)
{
    GtkWidget *csel;
    GtkWidget *color_button, *image;
    GdkColor gcolor;
    gretlRGB rgb;
    gint i;

    color_button = g_object_get_data(G_OBJECT(w), "color_button");
    csel = GTK_COLOR_SELECTION_DIALOG(w)->colorsel;

    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(csel), &gcolor);

    rgb.r = (unsigned char) (scale_round(gcolor.red));
    rgb.g = (unsigned char) (scale_round(gcolor.red));
    rgb.b = (unsigned char) (scale_round(gcolor.red));

    i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "colnum"));

    set_graph_palette(i, rgb);
    print_palette_string(gpcolors);

    /* update the "image" widget */
    image = g_object_get_data(G_OBJECT(color_button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(&rgb);
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(color_button), image);
    g_object_set_data(G_OBJECT(color_button), "image", image);
  
    gtk_widget_destroy(w);
}

static void color_cancel (GtkWidget *button, GtkWidget *w)
{
    gtk_widget_destroy(w);
}

GtkWidget *color_patch_button (int cnum)
{
    GtkWidget *image, *button;

    image = get_image_for_color(get_graph_color(cnum));

    if (image == NULL) {
	button = gtk_button_new_with_label(_("Select color"));
    } else {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
	g_object_set_data(G_OBJECT(button), "image", image);
    }	

    return button;
}

/* reset color patch button after selecting the option to
   restore the default plot colors */

void color_patch_button_reset (GtkWidget *button, int cnum)
{
    GtkWidget *image;

    image = g_object_get_data(G_OBJECT(button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(get_graph_color(cnum));
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(button), image);
    g_object_set_data(G_OBJECT(button), "image", image);

    if (cnum == BOXCOLOR || cnum == BOXCOLOR - 1) {
	print_palette_string(gpcolors);
    }
}

void graph_color_selector (GtkWidget *w, gpointer p)
{
    GtkWidget *cdlg;
    GtkWidget *button;
    gint i = GPOINTER_TO_INT(p);
    char colstr[8];
    const gretlRGB *rgb;
    GdkColor gcolor;

    rgb = get_graph_color(i);
    if (rgb == NULL) {
	fprintf(stderr, "graph_get_color(%d) gave NULL\n", i);
	return;
    }

    print_rgb_hash(colstr, rgb);
    gdk_color_parse(colstr, &gcolor);

    cdlg = gtk_color_selection_dialog_new(_("gretl: graph color selection"));

    g_object_set_data(G_OBJECT(cdlg), "colnum", GINT_TO_POINTER(i));
    g_object_set_data(G_OBJECT(cdlg), "color_button", w);

    gtk_color_selection_set_current_color(GTK_COLOR_SELECTION
					  (GTK_COLOR_SELECTION_DIALOG(cdlg)->colorsel),
					  &gcolor);					  

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->ok_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_select_callback), cdlg);

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->cancel_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_cancel), cdlg);

    gtk_widget_show(cdlg);
    gtk_window_set_modal(GTK_WINDOW(cdlg), TRUE);
}

/* end graph color selection apparatus */

#ifndef G_OS_WIN32

/* FIXME! */

static void real_set_workdir (GtkWidget *widget, dialog_t *dlg)
{
    const gchar *dirname;
    int err;

    dirname = edit_dialog_get_text(dlg);
    err = set_gretl_work_dir(dirname, &paths);

    if (err) {
	gui_errmsg(err);
	return;
    } else {
#if defined(HAVE_TRAMO) || defined(HAVE_X12A)
	set_tramo_x12a_status();
#endif
	close_dialog(dlg);
    }
}

void first_time_set_user_dir (void)
{
    DIR *test;

    /* see if an already-specified workdir exists */
    if (*paths.workdir != '\0') {
	test = opendir(paths.workdir);
	if (test != NULL) {
	    closedir(test);
	    return;
	}
    }
	
    /* user dir is not specified, or doesn't exist */
    edit_dialog (_("gretl: working directory"), 
                 _("You seem to be using gretl for the first time.\n"
		   "Please enter a directory for gretl user files."),
                 paths.workdir, 
                 real_set_workdir, NULL, 
                 CREATE_USERDIR, VARCLICK_NONE,
		 NULL);
}

#endif /* G_OS_WIN32 */

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

