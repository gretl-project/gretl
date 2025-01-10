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
#include "version.h"
#include "fileselect.h"
#include "gretl_www.h"
#include "dlgutils.h"
#include "textbuf.h"
#include "tabwin.h"
#include "build.h"
#include "addons_utils.h"

#ifndef GRETL_EDIT
#include "filelists.h"
#include "menustate.h"
#include "session.h"
#include "ssheet.h"
#include "selector.h"
#include "gpt_control.h"
#else
# if __linux__
#  define DEVEL_OPTS 1
#  include "gretl_edit.h"
# endif
#endif

#ifdef HAVE_GTKSV_COMPLETION
# include "completions.h"
#endif

#include "libset.h"
#include "texprint.h"
#include "uservar.h"
#include "gretl_foreign.h"

#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>

#if GTK_MAJOR_VERSION == 3
# define HAVE_GTK_FONT_CHOOSER 1
#else
# define HAVE_GTK_FONT_CHOOSER 0
#endif

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

#ifdef __APPLE__
# include "osx_open.h"
#endif

#if HAVE_GTK_FONT_CHOOSER
# include "fontfilter.h"
#else
# include "gtkfontselhack.h"
#endif

#if defined(__APPLE__) && defined(HAVE_MAC_THEMES)
# define MAC_THEMING
#endif

static char rcfile[FILENAME_MAX];
static char http_proxy[128];
static int use_proxy;

static ConfigPaths paths;

static void make_prefs_tab (GtkWidget *notebook, int tab, int console);
static void apply_prefs_changes (GtkWidget *widget, GtkWidget *parent);
static int common_read_rc_setup (int updated);

#ifndef G_OS_WIN32
static int read_gretlrc (void);
#endif

/* font handling */

static int tmpfontscale;
static char system_appfont[64];

#if defined(G_OS_WIN32)
# if defined(_WIN64)
static char fixedfontname[MAXLEN] = "Consolas 10";
static char default_fixedfont[64] = "Consolas 10";
# else
static char fixedfontname[MAXLEN] = "Courier New 10";
static char default_fixedfont[64] = "Courier New 10";
# endif
#elif defined(__APPLE__)
static char fixedfontname[MAXLEN] = "Menlo 13";
static char default_fixedfont[64] = "Menlo 13";
#else
static char fixedfontname[MAXLEN] = "monospace 10";
static char default_fixedfont[64] = "monospace 10";
#endif

#if defined(G_OS_WIN32)
static char appfontname[MAXLEN] = "";
#elif defined(__APPLE__)
static char appfontname[MAXLEN] = "Lucida Grande 13";
#else
static char appfontname[MAXLEN] = "sans 10";
#endif

PangoFontDescription *fixed_font;

/* end font handling */

/* status flag for console actually swallowed */
int swallow = 0;
/* status flag for console-swallowing selected */
static int swallow_pref = 0;

static int usecwd;
static int shellok;
static int manpref;
static int robust_z;
static int autoicon = 1;
static int session_prompt = 1;
static int keep_folder = 1;
static int tabbed_editor = 1;
static int tabbed_models = 0;
static int auto_collect = 0;

static int script_output_policy;
static char datapage[24] = "Gretl";
static char scriptpage[24] = "Gretl";
static char author_mail[32];
static char sview_style[32] = "classic";
static char graph_theme[24] = "dark2";
static char gpcolors[24];

static int hc_by_default;
static char langpref[32];
static char hc_xsect[5] = "HC1";
static char hc_tseri[5] = "HAC";
static char hc_panel[9] = "Arellano";
static char hc_garch[5] = "QML";

#ifdef HAVE_MPI
# ifdef G_OS_WIN32
static char mpi_pref[8] = "MS-MPI";
# else
static char mpi_pref[8] = "OpenMPI";
# endif
#endif

static int lcnumeric = 1;
static int icon_sizing = ICON_SIZE_AUTO;

static double graph_scale = 1.0;
static double graph_scales[] = {
    0.8, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0
};
static int n_graph_scales = G_N_ELEMENTS(graph_scales);

#if defined(MAC_THEMING)
static char themepref[16] = "Adwaita";
#elif defined(G_OS_WIN32)
static char themepref[16] = "Windows-10";
#endif

/* model table display variables */
static int modtab_colheads;
static gboolean modtab_tstats;
static gboolean modtab_pvalues;
static gboolean modtab_asterisks = TRUE;
static int modtab_digits = 4;
static gboolean modtab_decimals;

typedef enum {
    USERSET  = 1 << 0,  /* user-level variable */
    BOOLSET  = 1 << 1,  /* boolean value (user) */
    INTSET   = 1 << 2,  /* integer value (user) */
    FLOATSET = 1 << 3,  /* floating point value (user) */
    LISTSET  = 1 << 4,  /* user selection from fixed menu */
    RADIOSET = 1 << 5,  /* user int, from fixed menu */
    INVISET  = 1 << 6,  /* not visible in preferences dialog */
    FIXSET   = 1 << 7,  /* setting fixed by admin (Windows network use) */
    MACHSET  = 1 << 8,  /* "local machine" setting */
    BROWSER  = 1 << 9,  /* wants "Browse" button */
    RESTART  = 1 << 10, /* needs program restart to take effect */
    GOTSET   = 1 << 11, /* dynamic: already found a setting */
    SKIPSET  = 1 << 12, /* for string value, skip empty value */
    SPINSET  = 1 << 13  /* integer value: represented by spin-button */
} rcflags;

typedef struct {
    char *key;         /* config file variable name */
    char *description; /* string shown in the preferences dialog */
    char *link;        /* in case of radio button pair, alternate string */
    void *var;         /* pointer to variable */
    rcflags flags;     /* see above */
    int len;           /* storage size for string variable (also see Note) */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
} RCVAR;

/* Note: actually "len" above is overloaded: if an rc_var is of
   type BOOLSET and not part of a radio group, then a non-zero value
   for len will link the var's toggle button with the sensitivity of
   the preceding rc_var's entry field.  For example, the "use_proxy"
   button controls the sensitivity of the "http_proxy" entry widget.
*/

RCVAR rc_vars[] = {
    { "gretldir", N_("Main gretl directory"), NULL, paths.gretldir,
      MACHSET | BROWSER, sizeof paths.gretldir, TAB_MAIN, NULL },
    { "workdir", N_("User's gretl directory"), NULL, paths.workdir,
      INVISET, sizeof paths.workdir, TAB_MAIN, NULL },
#ifndef G_OS_WIN32
    { "winsize", N_("Remember main window size"), NULL, &winsize,
      BOOLSET, 0, TAB_MAIN, NULL },
#endif
#ifdef ENABLE_NLS
    { "lcnumeric", N_("Use locale setting for decimal point"), NULL, &lcnumeric,
      BOOLSET, 0, TAB_MAIN, NULL },
#endif
#if defined(MAC_THEMING) || defined(G_OS_WIN32)
    { "themepref", N_("Theme preference"), NULL, themepref,
      LISTSET | RESTART, 16, TAB_MAIN, NULL },
#endif
#if !defined(G_OS_WIN32) && !defined(__APPLE__)
    { "browser", N_("Web browser"), NULL, Browser,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
    { "shellok", N_("Allow shell commands"), NULL, &shellok,
      BOOLSET, 0, TAB_MAIN, NULL },
    { "autoicon", N_("Show icon view automatically"), NULL, &autoicon,
      BOOLSET, 0, TAB_MAIN, NULL },
    { "tabmodels", N_("Model viewer uses tabs"), NULL, &tabbed_models,
      BOOLSET, 0, TAB_MAIN, NULL },
    { "session_prompt", N_("Prompt to save session"), NULL, &session_prompt,
      BOOLSET, 0, TAB_MAIN, NULL },
    { "collect_plots", N_("Enable collecting plots"), NULL, &auto_collect,
      BOOLSET, 0, TAB_MAIN, NULL },
    { "swallow_console", N_("Main window includes console"), NULL, &swallow_pref,
      BOOLSET | RESTART, 0, TAB_MAIN, NULL },
    { "icon_sizing", N_("Toolbar icon size"), NULL, &icon_sizing,
      LISTSET | INTSET | RESTART, 0, TAB_MAIN, NULL },
    { "usecwd", N_("Set working directory from shell"), NULL, &usecwd,
      INVISET | BOOLSET | RESTART, 0, TAB_NONE, NULL },
    { "keepfolder", N_("File selector remembers folder"), NULL, &keep_folder,
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
#ifdef ENABLE_NLS
    { "langpref", N_("Language preference"), NULL, langpref,
      LISTSET | RESTART, 32, TAB_MAIN, NULL },
#endif
    { "graph_scale", N_("Default graph scale"), NULL, &graph_scale,
      LISTSET | FLOATSET, 0, TAB_MAIN, NULL },
    { "graph_theme", N_("Graph theme"), NULL, graph_theme,
      LISTSET, sizeof graph_theme, TAB_MAIN, NULL },
#if !defined(G_OS_WIN32) || !defined(PKGBUILD)
    { "gnuplot", N_("Command to launch gnuplot"), NULL, paths.gnuplot,
      MACHSET | BROWSER, MAXLEN, TAB_PROGS, NULL },
#endif
    { "Rcommand", N_("Command to launch GNU R"), NULL, Rcommand,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#ifdef G_OS_WIN32
    { "Rbin", N_("Path to R.exe"), NULL, paths.rbinpath,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
    { "latex", N_("Command to compile TeX files"), NULL, latex,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#if !defined(G_OS_WIN32) && !defined(__APPLE__)
    { "viewps", N_("Command to view postscript files"), NULL, viewps,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
    { "viewpdf", N_("Command to view PDF files"), NULL, viewpdf,
      MACHSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#endif
    { "calculator", N_("Calculator"), NULL, calculator,
      USERSET | BROWSER, MAXSTR, TAB_PROGS, NULL },
#ifdef HAVE_X12A
    { "x12a", N_("Path to x12arima (or x13)"), NULL, paths.x12a,
      MACHSET | BROWSER, sizeof paths.x12a, TAB_PROGS, NULL },
#endif
#ifdef HAVE_TRAMO
    { "tramo", N_("Path to tramo"), NULL, paths.tramo,
      MACHSET | BROWSER, sizeof paths.tramo, TAB_PROGS, NULL},
#endif
#ifdef USE_RLIB
    { "Rlib", N_("Path to R library"), NULL, paths.rlibpath,
      MACHSET | BROWSER, sizeof paths.rlibpath, TAB_PROGS, NULL},
#endif
    { "ox", N_("Path to oxl executable"), NULL, paths.oxlpath,
      MACHSET | BROWSER, sizeof paths.oxlpath, TAB_PROGS, NULL},
    { "octave", N_("Path to octave executable"), NULL, paths.octpath,
      MACHSET | BROWSER, sizeof paths.octpath, TAB_PROGS, NULL},
    { "stata", N_("Path to Stata executable"), NULL, paths.statapath,
      MACHSET | BROWSER, sizeof paths.statapath, TAB_PROGS, NULL},
    { "python", N_("Path to Python executable"), NULL, paths.pypath,
      MACHSET | BROWSER, sizeof paths.pypath, TAB_PROGS, NULL},
    { "julia", N_("Path to Julia executable"), NULL, paths.jlpath,
      MACHSET | BROWSER, sizeof paths.jlpath, TAB_PROGS, NULL},
#ifndef PKGBUILD
    { "lpsolve", N_("Path to lpsolve library"), NULL, paths.lppath,
      MACHSET | BROWSER, sizeof paths.lppath, TAB_PROGS, NULL},
#endif
#ifdef HAVE_MPI
    { "mpiexec", N_("Path to mpiexec"), NULL, paths.mpiexec,
      MACHSET | BROWSER, sizeof paths.mpiexec, TAB_MPI, NULL},
    { "mpi_hosts", N_("Path to MPI hosts file"), NULL, paths.mpi_hosts,
      MACHSET | BROWSER, sizeof paths.mpi_hosts, TAB_MPI, NULL},
    { "mpi_pref", N_("Installed MPI variant"), NULL, mpi_pref,
      LISTSET, 8, TAB_MPI, NULL},
#endif
    { "dbproxy", N_("HTTP proxy"), NULL, http_proxy,
      USERSET, sizeof http_proxy, TAB_NET, NULL },
    { "useproxy", N_("Use HTTP proxy"), NULL, &use_proxy,
      BOOLSET, 1, TAB_NET, NULL },
    { "Fixed_font", N_("Monospaced font"), NULL, fixedfontname,
      USERSET, sizeof fixedfontname, TAB_NONE, NULL },
    { "App_font", N_("Menu font"), NULL, appfontname,
      USERSET, sizeof appfontname, TAB_NONE, NULL },
    { "DataPage", "Default data page", NULL, datapage,
      INVISET, sizeof datapage, TAB_NONE, NULL },
    { "ScriptPage", "Default script page", NULL, scriptpage,
      INVISET, sizeof scriptpage, TAB_NONE, NULL },
    { "Png_font", N_("PNG graph font"), NULL, paths.pngfont,
      INVISET, sizeof paths.pngfont, TAB_NONE, NULL },
    { "Gp_extra_colors", N_("Gnuplot extra colors"), NULL, gpcolors,
      INVISET, sizeof gpcolors, TAB_NONE, NULL },
    { "tabwidth", N_("Number of spaces per tab"), NULL, &tabwidth,
      INTSET | SPINSET, 0, TAB_EDITOR, NULL },
    { "smarttab", N_("\"Smart\" Tab and Enter"), NULL, &smarttab,
      BOOLSET, 0, TAB_EDITOR, NULL },
    { "script_line_numbers", N_("Show line numbers"), NULL, &script_line_numbers,
      BOOLSET, 0, TAB_EDITOR, NULL },
    { "tabedit", N_("Script editor uses tabs"), NULL, &tabbed_editor,
      BOOLSET, 0, TAB_EDITOR, NULL },
#ifdef HAVE_GTKSV_COMPLETION
    { "hansl_completion", N_("Auto-completion"), NULL, &hansl_completion,
      LISTSET | INTSET, 0, TAB_EDITOR, NULL },
    { "console_completion", N_("Auto-completion"), NULL, &console_completion,
      LISTSET | INTSET | INVISET, 0, TAB_EDITOR, NULL },
#endif
    { "script_auto_bracket", N_("Enable auto-brackets"), NULL, &script_auto_bracket,
      BOOLSET, 0, TAB_EDITOR, NULL },
    { "sview_style", N_("Highlighting style"), NULL, &sview_style,
      LISTSET, sizeof sview_style, TAB_EDITOR, NULL },
    { "main_width", "main window width", NULL, &mainwin_width,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_height", "main window height", NULL, &mainwin_height,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_x", "main window x position", NULL, &main_x,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "main_y", "main window y position", NULL, &main_y,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "script_output_policy", "stickiness of output", NULL, &script_output_policy,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "HC_by_default", N_("Use robust covariance matrix by default"), NULL,
      &hc_by_default, BOOLSET, 0, TAB_VCV, NULL },
    { "robust_z", N_("Use the normal distribution for robust p-values"), NULL,
      &robust_z, BOOLSET, 0, TAB_VCV, NULL },
    { "HC_xsect", N_("For cross-sectional data"), NULL, hc_xsect,
      LISTSET, 5, TAB_VCV, NULL },
    { "HC_tseri", N_("For time-series data"), NULL, hc_tseri,
      LISTSET, 5, TAB_VCV, NULL },
    { "HC_garch", N_("For GARCH estimation"), NULL, hc_garch,
      LISTSET, 5, TAB_VCV, NULL },
    { "HC_panel", N_("For panel data"), NULL, hc_panel,
      LISTSET, 9, TAB_VCV, NULL },
    { "manpref", N_("PDF manual preference"), NULL, &manpref,
      LISTSET | INTSET, 0, TAB_MAIN, NULL },
    { "modtab_colheads", "Model table column heads", NULL, &modtab_colheads,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "modtab_tstats", "Model table t-ratios", NULL, &modtab_tstats,
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
    { "modtab_pvalues", "Model table p-values", NULL, &modtab_pvalues,
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
    { "modtab_asterisks", "Model table asterisks", NULL, &modtab_asterisks,
      INVISET | BOOLSET, 0, TAB_NONE, NULL },
    { "modtab_digits", "Model table digits", NULL, &modtab_digits,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "modtab_decimals", "Model table decimal places", NULL, &modtab_decimals,
      INVISET | INTSET, 0, TAB_NONE, NULL },
    { "author_mail", "Package author email", NULL, &author_mail,
      INVISET | SKIPSET, sizeof author_mail, TAB_NONE, NULL },
    { NULL, NULL, NULL, NULL, 0, 0, TAB_NONE, NULL }
};

#ifdef DEVEL_OPTS /* gretl_edit on Linux */

GtkWidget *cli_selector;
GtkWidget *env_selector;

#endif

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

void set_author_mail (const char *s)
{
    if (s != NULL && strlen(s) < sizeof author_mail) {
	strcpy(author_mail, s);
    } else {
	author_mail[0] = '\0';
    }
}

const char *get_author_mail (void)
{
    return author_mail;
}

const char *get_sourceview_style (void)
{
    return sview_style;
}

#ifndef GRETL_EDIT

int autoicon_on (void)
{
    if (dataset != NULL && dataset->v > 0) {
	return autoicon;
    } else if (n_user_matrices() > 0) {
	return autoicon;
    } else if (n_user_bundles() > 0) {
	return autoicon;
    } else {
	return 0;
    }
}

#endif

int dark_theme_active (void)
{
    int ret = 0;

#if defined(G_OS_WIN32) || defined(MAC_THEMING)
    if (strstr(themepref, "Dark") ||
        strstr(themepref, "dark")) {
	ret = 1;
    }
#else
    /* Linux, etc. */
    GSettings *settings;

    settings = g_settings_new("org.gnome.desktop.interface");
    if (settings != NULL) {
	gchar *theme = g_settings_get_string(settings, "gtk-theme");

	if (theme != NULL) {
	    if (strstr(theme, "Dark") || strstr(theme, "dark")) {
		ret = 1;
	    }
	    g_free(theme);
	}
	g_object_unref(settings);
    }
#endif

    return ret;
}

const char *blue_for_text (void)
{
    static char blue[8];

    if (*blue == '\0') {
	if (dark_theme_active()) {
	    strcpy(blue, "#B0EFEE");
	} else {
	    strcpy(blue, "blue");
	}
    }

    return blue;
}

int get_icon_sizing (void)
{
    return icon_sizing;
}

int use_tabbed_editor (void)
{
    return tabbed_editor;
}

int use_tabbed_model_viewer (void)
{
    return tabbed_models;
}

int session_prompt_on (void)
{
    return session_prompt;
}

void set_session_prompt (int val)
{
    session_prompt = val;
}

int get_keep_folder (void)
{
    return keep_folder;
}

void set_script_output_policy (int p, windata_t *vwin)
{
    script_output_policy = p;

    if (script_output_policy < 0 ||
	script_output_policy > OUTPUT_POLICY_NEW_WINDOW) {
	/* invalid setting */
	script_output_policy = OUTPUT_POLICY_REPLACE;
    }

    set_reuseable_output_window(p, vwin);
}

int get_script_output_policy (void)
{
    return script_output_policy;
}

void get_model_table_prefs (int *colheads,
			    int *use_tstats,
			    int *do_pvals,
			    int *do_asts,
			    int *figs,
			    char *fmt)
{
    *colheads   = modtab_colheads;
    *use_tstats = modtab_tstats;
    *do_pvals   = modtab_pvalues;
    *do_asts    = modtab_asterisks;
    *figs       = modtab_digits;
    *fmt        = modtab_decimals ? 'f' : 'g';
}

void set_model_table_prefs (int colheads,
			    int use_tstats,
			    int do_pvals,
			    int do_asts,
			    int figs,
			    char fmt)
{
    modtab_colheads  = colheads;
    modtab_tstats    = use_tstats;
    modtab_pvalues   = do_pvals;
    modtab_asterisks = do_asts;
    modtab_digits    = figs;
    modtab_decimals  = (fmt == 'f')? 1 : 0;
}

static gretlopt update_paths_opt = OPT_NONE;

void force_english_help (void)
{
    update_paths_opt |= OPT_N;
    gretl_update_paths(&paths, update_paths_opt);
}

static int fontname_get_size (const char *fontname)
{
    char *p = strrchr(fontname, ' ');

    return (p != NULL)? atoi(p+1) : 10;
}

static void fontname_set_size (char *fontname, int size)
{
    char *p = strrchr(fontname, ' ');

    if (p != NULL) {
	sprintf(p, " %d", size);
    }
}

static int font_is_changed (const char *f_new,
			    const char *f_old)
{
    if (strcmp(f_new, f_old)) {
	return 1;
    } else if (tmpfontscale > 0 &&
	       fontname_get_size(f_new) != tmpfontscale) {
	return 1;
    } else {
	return 0;
    }
}

void set_fixed_font (const char *fontname, int remember)
{
    if (fontname == NULL) {
	/* initial set-up */
	fixed_font = pango_font_description_from_string(fixedfontname);
    } else if (font_is_changed(fontname, fixedfontname)) {
	/* changed via the GUI */
	if (fixed_font != NULL) {
	    pango_font_description_free(fixed_font);
	}
	fixed_font = pango_font_description_from_string(fontname);
	if (remember) {
	    strcpy(fixedfontname, fontname);
	}
	infobox(_("This change will apply to newly opened windows"));
    }
}

void update_persistent_graph_font (void)
{
    strcpy(paths.pngfont, gretl_png_font());
}

const char *get_app_fontname (void)
{
    return appfontname;
}

const char *get_fixed_fontname (void)
{
    return fixedfontname;
}

static void record_system_appfont (GtkSettings *settings,
				   gchar **pfont)
{
    g_object_get(G_OBJECT(settings), "gtk-font-name", pfont, NULL);
#if defined(G_OS_WIN32)
    get_default_windows_app_font(system_appfont);
#else
    if (*pfont != NULL) {
	strcpy(system_appfont, *pfont);
    } else {
# if defined(__APPLE__)
	strcpy(system_appfont, "Lucida Grande 13");
# else
	strcpy(system_appfont, "sans 10");
# endif
    }
#endif
}

#ifdef G_OS_WIN32

static void win32_set_font (const char *fontname,
			    GtkSettings *settings)
{
    gchar *rc;

    rc = g_strdup_printf("style \"myfont\" {\n"
			 "  font_name = \"%s\"\n}\n"
			 "widget_class \"*\" style \"myfont\"\n"
			 "gtk-font-name = \"%s\"\n",
			 fontname, fontname);
    gtk_rc_parse_string(rc);
    g_object_set(G_OBJECT(settings), "gtk-font-name", fontname, NULL);
    g_free(rc);
}

#endif

void set_app_font (const char *fontname, int remember)
{
    static int default_recorded;
    GtkSettings *settings;
    gchar *deffont = NULL;

#if 0
    fprintf(stderr, "set_app_font: fontname='%s', remember=%d, "
	    "appfontname='%s'\n", fontname, remember, appfontname);
#endif

    if (fontname != NULL && *fontname == '\0') {
	return;
    }

    settings = gtk_settings_get_default();

    /* not font-related but, dammit, we want these! */
    g_object_set(G_OBJECT(settings), "gtk-menu-images", TRUE, NULL);

    if (!default_recorded) {
	record_system_appfont(settings, &deffont);
#if 0
	fprintf(stderr, "record app font default: system '%s', gtk '%s'\n",
		system_appfont, deffont);
#endif
	default_recorded = 1;
    }

#ifdef G_OS_WIN32
    if (fontname == NULL && *appfontname == '\0') {
	/* as we're called at initial startup from gretl.c */
	strcpy(appfontname, system_appfont);
    }
#endif

    /* check for nothing else to be done */
    if (deffont != NULL) {
	const char *test = fontname ? fontname : appfontname;
	int noop = 0;

	if (!font_is_changed(test, deffont)) {
	    noop = 1;
	}
	g_free(deffont);
	if (noop) {
	    return;
	}
    }

    if (fontname == NULL) {
	/* just loading @appfontname (pre-checked) */
#ifdef G_OS_WIN32
	win32_set_font(appfontname, settings);
#else
	g_object_set(G_OBJECT(settings), "gtk-font-name", appfontname, NULL);
#endif
    } else {
	/* loading a user-specified font: check that it works */
	GtkWidget *w;
	PangoFontDescription *pfd;
	PangoContext *pc;
	PangoFont *font;

	w = gtk_label_new("text");
	pfd = pango_font_description_from_string(fontname);
	pc = gtk_widget_get_pango_context(w);
	font = pango_context_load_font(pc, pfd);

	if (font != NULL) {
	    /* OK, found it */
	    if (remember) {
		strcpy(appfontname, fontname);
	    }
#ifdef G_OS_WIN32
	    win32_set_font(fontname, settings);
#else
	    g_object_set(G_OBJECT(settings), "gtk-font-name", fontname, NULL);
#endif
	    g_object_unref(font);
	}

	gtk_widget_destroy(w);
	pango_font_description_free(pfd);
    }
}

static int write_OK (gchar *dirname)
{
    int ok = 0;

    if (gretl_mkdir(dirname) == 0) {
	gchar *test = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest");

	if (test != NULL) {
	    ok = (gretl_test_fopen(test, "w") == 0);
	    g_free(test);
	}
    }

    return ok;
}

void set_gretl_startdir (void)
{
    if (usecwd) {
	char *test = getenv("GRETL_STARTDIR");
	gchar *startdir = NULL;

	/* the environment variable check is mostly for the macOS
	   package */
	if (test != NULL) {
	    startdir = g_strdup(test);
	} else {
	    startdir = g_get_current_dir();
	}

	if (startdir != NULL) {
	    int err = gretl_set_path_by_name("workdir", startdir);

	    if (err) {
		fprintf(stderr, "%s\n", gretl_errmsg_get());
	    } else {
		fprintf(stderr, "working dir = '%s'\n", startdir);
	    }
	    g_free(startdir);
	}
    }
}

static void get_pkg_save_dir (char *dirname, int action)
{

    const char *subdir = NULL;
    int try_sysdir = 1;
    int ok = 0;

    if (action == SAVE_FUNCTIONS) {
	subdir = "functions";
    } else if (action == SAVE_DATA_PKG) {
	subdir = "data";
    } else if (action == SAVE_REMOTE_DB) {
	subdir = "db";
    } else {
	return;
    }

#ifdef G_OS_WIN32
    try_sysdir = 0;
#endif

    if (try_sysdir) {
	/* try 'system' location first */
	sprintf(dirname, "%s%s", gretl_home(), subdir);
	ok = write_OK(dirname);
    }

    if (!ok) {
	/* go to user's dotdir */
	sprintf(dirname, "%s%s", gretl_dotdir(), subdir);
	ok = write_OK(dirname);
    }

    if (!ok) {
	*dirname = '\0';
    }
}

void get_default_dir_for_action (char *s, int action)
{
    *s = '\0';

    if (action == SAVE_FUNCTIONS ||
	action == SAVE_DATA_PKG ||
	action == SAVE_REMOTE_DB) {
	get_pkg_save_dir(s, action);
    } else {
	strcpy(s, gretl_workdir());
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

#ifdef __APPLE__

static int alt_ok (const char *prog)
{
    char *p, test[MAXSTR];
    int tr = strstr(prog, "tramo") != NULL;
    int ok;

    strcpy(test, gretl_home());
    p = strstr(test, "share/gretl");
    if (p != NULL) {
    	*p = 0;
    }

    if (tr) {
         strcat(test, "tramo/tramo");
    } else {
         strcat(test, "x12arima/x12a");
    }

    ok = check_for_program(test);

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

#if defined(HAVE_TRAMO) && !defined(GRETL_EDIT)

#define tramo_ts(d) ((d)->structure == TIME_SERIES && \
                     (d->pd == 1 || d->pd == 4 || d->pd == 12))

static int tramo_ok = 0;

int get_tramo_ok (void)
{
    return tramo_ok && tramo_ts(dataset);
}

static void set_tramo_status (void)
{
    const char *tramodir = gretl_tramo_dir();
    int gui_up = (mdata != NULL);
    int ok = 0;

    if (*tramodir != '\0') {
	const char *tramo = gretl_tramo();

	ok = check_for_program(tramo);
# ifdef __APPLE__
	if (!ok) {
	    ok = alt_ok(tramo);
	}
# endif
    }

    if (tramo_ok && !ok && gui_up) {
	warnbox_printf("Invalid path for %s", "TRAMO");
    }

    tramo_ok = ok;

#ifndef GRETL_EDIT
    if (gui_up) {
	flip(mdata->ui, "/menubar/Variable/Tramo", get_tramo_ok());
    }
#endif
}

#endif /* HAVE_TRAMO */

#if defined(HAVE_X12A) && !defined(GRETL_EDIT)

#define x12_ts(d) ((d)->structure == TIME_SERIES && \
                   (d->pd == 4 || d->pd == 12))

static int x12a_ok = 0;

int get_x12a_ok (void)
{
    return x12a_ok && x12_ts(dataset);
}

static void set_x12a_status (void)
{
    const char *x12adir = gretl_x12_arima_dir();
    int gui_up = (mdata != NULL);
    int ok = 0;

    if (*x12adir != '\0') {
	const char *x12a = gretl_x12_arima();

	ok = check_for_program(x12a);
# ifdef __APPLE__
	if (!ok) {
	    ok = alt_ok(x12a);
	}
# endif
    }

    if (x12a_ok && !ok && gui_up) {
	warnbox_printf("Invalid path for %s", "X-12-ARIMA");
    }

    x12a_ok = ok;

#ifndef GRETL_EDIT
    if (gui_up) {
	flip(mdata->ui, "/menubar/Variable/X12A",
	     get_x12a_ok());
    }
#endif
}

#endif /* HAVE_X12A */

#ifndef G_OS_WIN32

static void root_check (void)
{
    if (getuid() == 0) {
	int resp;

	resp = yes_no_dialog("gretl", _("You seem to be running gretl "
					"as root.  Do you really want to do this?"),
			     NULL);
	if (resp == GRETL_NO) {
	    exit(EXIT_FAILURE);
	}
    }
}

int gretl_config_init (int ignore_rc)
{
    int err = 0;

    get_gretl_rc_path(rcfile);
    if (ignore_rc) {
	err = common_read_rc_setup(0);
    } else {
	err = read_gretlrc();
    }
    set_gretl_startdir();
    root_check();

    return err;
}

#endif /* !G_OS_WIN32 */

static void highlight_preferences_entry (const char *vname)
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

static void preferences_dialog_canceled (GtkWidget *w, int *c)
{
    *c = 1;
}

static int page_has_help (int p)
{
    return p == TAB_EDITOR || p == TAB_VCV;
}

static void sensitize_prefs_help (GtkNotebook *book,
				  gpointer arg1,
				  guint newpg,
				  GtkWidget *button)
{
    gtk_widget_set_sensitive(button, page_has_help(newpg + 1));
}

static void show_prefs_help (GtkWidget *w, GtkWidget *notebook)
{
    gint page;

    page = gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));

    if (page + 1 == TAB_EDITOR) {
	show_gui_help(EDITOR);
    } else if (page + 1 == TAB_VCV) {
	show_gui_help(HCCME);
    }
}

/* To record state of preferences dialogs, and avoid opening
   both the 'global' one and the one specific to the console
   simultaneously, which would lead to bad effects.
*/
static GtkWidget *all_prefs;
static GtkWidget *console_prefs;

static void preferences_dialog_destroyed (GtkWidget *w,
					  GtkWidget **pw)
{
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	rc_vars[i].widget = NULL;
    }
    if (pw != NULL) {
	*pw = NULL;
    }
}

int preferences_dialog (int page, const char *varname, GtkWidget *parent)
{
    GtkWidget *dialog = all_prefs;
    GtkWidget *notebook;
    GtkWidget *button;
    GtkWidget *hbox;
    GtkWidget *vbox;
    int canceled = 0;

    if (dialog != NULL) {
	gtk_window_present(GTK_WINDOW(dialog));
	return 0;
    } else if (console_prefs != NULL) {
	gtk_widget_destroy(console_prefs);
    }

    dialog = gretl_dialog_new(_("gretl: preferences"), parent,
			      GRETL_DLG_RESIZE | GRETL_DLG_BLOCK);
    all_prefs = dialog;
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dialog), FALSE);
#endif
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_set_spacing(GTK_BOX(vbox), 2);

    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(preferences_dialog_destroyed),
		     &all_prefs);

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);

    make_prefs_tab(notebook, TAB_MAIN, 0);
    make_prefs_tab(notebook, TAB_PROGS, 0);
    make_prefs_tab(notebook, TAB_EDITOR, 0);
    make_prefs_tab(notebook, TAB_NET, 0);
#ifndef GRETL_EDIT
    make_prefs_tab(notebook, TAB_VCV, 0);
#endif
#ifdef HAVE_MPI
    make_prefs_tab(notebook, TAB_MPI, 0);
#endif

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Apply button */
    button = apply_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(apply_prefs_changes), parent);
    gtk_widget_grab_default(button);

    /* Cancel button */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(preferences_dialog_canceled),
		     &canceled);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    /* OK button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(apply_prefs_changes), parent);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    /* Help button */
    button = context_help_button(hbox, -1);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(show_prefs_help),
		     notebook);
    gtk_widget_set_sensitive(button, page_has_help(page));
    g_signal_connect(G_OBJECT(notebook), "switch-page",
		     G_CALLBACK(sensitize_prefs_help),
		     button);

    if (page > 1 && page < TAB_MAX) {
	page--;
	/* "show" the target page first (see GtkNoteBook API doc) */
	gtk_widget_show(gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook), page));
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), page);
    }

    if (varname != NULL) {
	highlight_preferences_entry(varname);
    }

    gtk_widget_show_all(dialog);

    return canceled;
}

static void refocus_console (GtkWidget *widget, GtkWidget *caller)
{
    gtk_widget_grab_focus(caller);
}

int console_prefs_dialog (GtkWidget *caller)
{
    GtkWidget *dialog = console_prefs;
    GtkWidget *parent;
    GtkWidget *button;
    GtkWidget *hbox;
    GtkWidget *vbox;
    int canceled = 0;

    if (dialog != NULL) {
	gtk_window_present(GTK_WINDOW(dialog));
    } else if (all_prefs != NULL) {
	gtk_widget_destroy(all_prefs);
    }

    parent = swallow ? mdata->main : caller;
    dialog = gretl_dialog_new(_("gretl: preferences"), parent,
			      GRETL_DLG_RESIZE | GRETL_DLG_BLOCK);
    console_prefs = dialog;
#if GTK_MAJOR_VERSION < 3
    gtk_dialog_set_has_separator(GTK_DIALOG(dialog), FALSE);
#endif
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_set_spacing(GTK_BOX(vbox), 2);

    if (swallow) {
	g_signal_connect(G_OBJECT(dialog), "destroy",
			 G_CALLBACK(refocus_console),
			 caller);
    }
    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(preferences_dialog_destroyed),
		     &console_prefs);

    make_prefs_tab(vbox, TAB_EDITOR, 1);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(preferences_dialog_canceled),
		     &canceled);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    /* OK button */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(apply_prefs_changes), caller);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    gtk_widget_show_all(dialog);

    return canceled;
}

static void flip_sensitive (GtkWidget *w, gpointer data)
{
    GtkWidget *entry = GTK_WIDGET(data);

    gtk_widget_set_sensitive(entry, button_is_active(w));
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
    GtkWidget *parent = g_object_get_data(G_OBJECT(w), "parent");
    int code = SET_PROG;

#ifdef HAVE_MPI
    if (!strcmp(rc->key, "mpi_hosts")) {
	file_selector_with_parent(SET_OTHER, FSEL_DATA_MISC, rc->var, parent);
	return;
    }
#endif

    if (strstr(rc->description, "directory") != NULL) {
	code = SET_DIR;
    }

    file_selector_with_parent(code, FSEL_DATA_MISC, rc->var, parent);
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

static gboolean try_switch_locale (GtkComboBox *box, gpointer p)
{
    static int lasterr;
    gchar *langstr;
    int err;

    if (lasterr) {
	lasterr = 0;
	return FALSE;
    }

    langstr = combo_box_get_active_text(box);
    err = test_locale(langstr);
    g_free(langstr);

    if (err) {
	lasterr = err;
	gui_errmsg(err);
	gretl_error_clear();
	gtk_combo_box_set_active(box, 0);
    }

    return FALSE;
}

static gboolean try_switch_style (GtkComboBox *box, GtkWidget *text)
{
    gchar *style = combo_box_get_active_text(box);

    set_style_for_textview(text, style);
    g_free(style);

    return FALSE;
}

static GtkWidget *embed_style_sampler (GtkWidget *vbox)
{
    GtkWidget *hbox, *text;

    hbox = gtk_hbox_new(TRUE, 5);
    text = create_sample_source(sview_style);
    gtk_box_pack_start(GTK_BOX(hbox), text, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 15);

    return text;
}

static const char *hc_strs[] = {
    "HC0", "HC1", "HC2", "HC3", "HC3a", "HAC"
};

static const char *hc_panel_strs[] = {
    "Arellano", "PCSE", "SCC"
};

static const char **get_list_setting_strings (void *var, int *n)
{
    static const char *garch_strs[] = {
	"QML", "BW"
    };
    static const char *manpref_strs[] = {
        N_("English (US letter paper)"),
        N_("English (A4 paper)"),
        N_("Translation, if available")
    };
#ifdef HAVE_GTKSV_COMPLETION
    static const char *completion_strs[] = {
	N_("none"),
	N_("automatic, as you type"),
	N_("on demand, via Tab")
    };
#endif
    static const char *icon_sizing_strs[] = {
	N_("Automatic"),
	N_("small"),
	N_("medium")
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
    } else if (var == &manpref) {
	strs = manpref_strs;
	*n = sizeof manpref_strs / sizeof manpref_strs[0];
#ifdef HAVE_GTKSV_COMPLETION
    } else if (var == &hansl_completion || var == &console_completion) {
	strs = completion_strs;
	*n = sizeof completion_strs / sizeof completion_strs[0];
#endif
    } else if (var == &icon_sizing) {
	strs = icon_sizing_strs;
	*n = sizeof icon_sizing_strs / sizeof icon_sizing_strs[0];
    } else if (var == sview_style) {
	strs = get_sourceview_style_ids(n);
    } else if (var == graph_theme) {
	strs = get_graph_theme_ids(n);
    }

#ifdef HAVE_MPI
    else if (var == mpi_pref) {
# ifdef G_OS_WIN32
	static const char *mpi_strs[] = {
	    "MS-MPI" /* FIXME? */
	};
# else
	static const char *mpi_strs[] = {
	    "OpenMPI", "MPICH"
	};
#endif

	strs = mpi_strs;
	*n = sizeof mpi_strs / sizeof mpi_strs[0];
    }
#endif

#if defined(MAC_THEMING) && GTK_MAJOR_VERSION < 3
    else if (var == themepref) {
	static const char *theme_strs[] = {
	    "Adwaita", "Adwaita-dark", "Clearlooks", "Raleigh"
	};

	strs = theme_strs;
	*n = sizeof theme_strs / sizeof theme_strs[0];
    }
#elif defined(MAC_THEMING) && GTK_MAJOR_VERSION == 3
    else if (var == themepref) {
	static const char *theme_strs[] = {
	    "Adwaita", "Adwaita-dark"
	};

	strs = theme_strs;
	*n = sizeof theme_strs / sizeof theme_strs[0];
    }
#elif defined(G_OS_WIN32) && GTK_MAJOR_VERSION < 3
    else if (var == themepref) {
	static const char *theme_strs[] = {
            "Windows-10", "Windows-10-Dark", "MS-Windows",
            "Clearlooks", "Raleigh"
	};

	strs = theme_strs;
	*n = sizeof theme_strs / sizeof theme_strs[0];
    }
#elif defined(G_OS_WIN32) && GTK_MAJOR_VERSION == 3
    else if (var == themepref) {
	static const char *theme_strs[] = {
            "Windows-10", "Windows-10-Dark", "Adwaita",
            "Windows 7"
	};

	strs = theme_strs;
	*n = sizeof theme_strs / sizeof theme_strs[0];
    }
#endif

    return strs;
}

static const char **get_radio_setting_strings (void *var, int *n)
{
    /* unused at present, but may be activated again at
       some point? */
    *n = 0;
    return NULL;
}

#ifndef GRETL_EDIT

const char *get_default_hc_string (int ci)
{
    if (ci == GARCH) {
	int k = libset_get_int(GARCH_ALT_VCV);

	return (k == ML_BW)? "BW" : "QML";
    } else if (!robust_conf(ci)) {
	return "QML";
    } else {
	int xsect = dataset_is_cross_section(dataset);
	int tseries = dataset_is_time_series(dataset);

	if (tseries && libset_get_bool(FORCE_HC)) {
	    xsect = 1;
	} else if (ci == VAR) {
	    xsect = 1;
	}

	if (xsect) {
	    return hc_strs[libset_get_int(HC_VERSION)];
	} else if (tseries) {
	    /* (and not forced to an HC variant) */
	    return "HAC";
	} else {
	    /* panel */
	    return hc_panel_strs[libset_get_int(PANEL_ROBUST)];
	}
    }
}

#endif /* not GRETL_EDIT */

static int non_console_var (void *ptr)
{
#ifdef HAVE_GTKSV_COMPLETION
    return (ptr == &smarttab || ptr == &script_line_numbers ||
	    ptr == &tabbed_editor || ptr == &tabwidth ||
	    ptr == &hansl_completion);
#else
    return (ptr == &smarttab || ptr == &script_line_numbers ||
	    ptr == &tabbed_editor || ptr == &tabwidth);
#endif
}

static int console_only (void *ptr)
{
#ifdef HAVE_GTKSV_COMPLETION
    return (ptr == &console_completion);
#else
    return 0;
#endif
}

static void
get_table_sizes (int page, int *n_str, int *n_bool, int *n_browse,
		 int *n_list, int console)
{
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	if (rc_vars[i].tab == page) {
	    if (console && non_console_var(rc_vars[i].var)) {
		continue;
	    } else if (!console && console_only(rc_vars[i].var)) {
		continue;
	    }
	    if (rc_vars[i].flags & BROWSER) {
		*n_browse += 1;
	    } else if (rc_vars[i].flags & LISTSET) {
		*n_list += 1;
	    } else if (rc_vars[i].flags & SPINSET) {
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
    if (button_is_active(w)) {
	gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "action"));

	*v = i;
    }
}

static void table_attach_fixed (GtkTable *table,
				GtkWidget *child,
				guint left, guint right,
				guint top, guint bottom)
{
    gtk_table_attach(table, child, left, right, top, bottom,
		     0, 0, 0, 0);
}

static void themes_page (GtkButton *button, gpointer p)
{
    if (browser_open("http://gretl.sourceforge.net/pub/gretl/plots/")) {
	errbox("Failed to open URL");
    }
}

static void add_themes_examples_button (GtkWidget *hbox)
{
    GtkWidget *b;

    b = gtk_button_new_with_label(_("Examples"));
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 5);
    g_signal_connect(G_OBJECT(b), "clicked",
		     G_CALLBACK(themes_page), NULL);
}

static GtkWidget *scroller_page (GtkWidget *vbox)
{
    GtkWidget *scroller;

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				   GTK_POLICY_NEVER,
				   GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller),
					  vbox);
    return scroller;
}

static void make_prefs_tab (GtkWidget *notebook, int tab,
			    int console)
{
    GtkWidget *b_table = NULL; /* widget table for boolean switches */
    GtkWidget *s_table = NULL; /* widget table for string variables */
    GtkWidget *l_table = NULL; /* widget table for list selections */
    GtkWidget *vbox, *w = NULL;
    GtkWidget *page;
    int s_len = 1, b_len = 0, l_len = 1;
    int s_cols, b_cols = 0, l_cols = 0;
    int b_col = 0;
    int n_str = 0;
    int n_bool = 0;
    int n_browse = 0;
    int n_list = 0;
    RCVAR *rc;
    int i;

    if (console) {
	vbox = notebook; /* not really a notebook! */
    } else {
	vbox = gtk_vbox_new(FALSE, 0);
	gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);

	if (tab == TAB_MAIN) {
	    w = gtk_label_new(_("General"));
	} else if (tab == TAB_PROGS) {
	    w = gtk_label_new(_("Programs"));
	} else if (tab == TAB_EDITOR) {
	    w = gtk_label_new(_("Editor"));
	} else if (tab == TAB_NET) {
	    w = gtk_label_new(_("Network"));
	} else if (tab == TAB_VCV) {
	    w = gtk_label_new(_("HCCME"));
#ifdef HAVE_MPI
	} else if (tab == TAB_MPI) {
	    w = gtk_label_new(_("MPI"));
#endif
	}
	if (tab == TAB_PROGS) {
	    page = scroller_page(vbox);
	} else {
	    page = vbox;
	}

	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, w);
    }

    get_table_sizes(tab, &n_str, &n_bool, &n_browse, &n_list, console);

    s_cols = (n_browse > 0)? 3 : 2;

    if (tab == TAB_VCV && n_list > 0) {
	/* VCV tab -- put the list entries first, right aligned */
	l_cols = 2;
	l_table = gtk_table_new(l_len, l_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(l_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(l_table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), l_table, FALSE, FALSE, 0);
    }

    if (n_str > 0) {
	s_table = gtk_table_new(s_len, s_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(s_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(s_table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), s_table, FALSE, FALSE, 0);
    }

    if (n_bool > 0) {
	b_cols = 2;
	b_table = gtk_table_new(1, b_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(b_table), 5);
	gtk_table_set_col_spacings(GTK_TABLE(b_table), 5);
	gtk_box_pack_start(GTK_BOX(vbox), b_table, FALSE, FALSE, 10);
    }

    if (tab != TAB_VCV && n_list > 0) {
	/* non-VCV tab -- list entries come last, and we
	   use an hbox to pack them left */
	GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

	l_cols = 2;
	l_table = gtk_table_new(l_len, l_cols, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(l_table), 10);
	gtk_table_set_col_spacings(GTK_TABLE(l_table), 10);
	gtk_box_pack_start(GTK_BOX(hbox), l_table, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);
    }

    for (i=0; rc_vars[i].key != NULL; i++) {
	rc = &rc_vars[i];

	if (rc->tab != tab) {
	    /* the item is not on this page */
	    continue;
	} else if (console && non_console_var(rc->var)) {
	    rc->widget = NULL;
	    continue;
	} else if (!console && console_only(rc->var)) {
	    rc->widget = NULL;
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

	    /* special case: link between toggle and preceding entry */
	    if (rc->len && !(rc->flags & FIXSET)) {
		gtk_widget_set_sensitive(rc_vars[i-1].widget,
					 button_is_active(rc->widget));
		g_signal_connect(G_OBJECT(rc->widget), "clicked",
				 G_CALLBACK(flip_sensitive),
				 rc_vars[i-1].widget);
	    }

	    b_col++;

	    if (tab == TAB_VCV || b_col == 2) {
		/* boolean strings under TAB_VCV are too long to
		   appear column-wise
		*/
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

	    /* then a first button */
	    button = gtk_radio_button_new_with_label(group, _(rc->link));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), button,
				      b_col, b_col + 1,
				      b_len - 2, b_len - 1);
	    if (!rcval) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
	    }

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
	    } else {
		rcol = b_col + 1;
	    }

	    b_len++;
	    b = gtk_label_new(_(rc->description));
	    gtk_table_attach_defaults(GTK_TABLE(b_table), b,
				      b_col, rcol,
				      b_len - 1, b_len);

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
		group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
		g_signal_connect(G_OBJECT(b), "clicked",
				 G_CALLBACK(radio_change_value),
				 rc->var);
	    }
	} else if (rc->flags & LISTSET) {
	    int langs = rc->var == langpref;
	    int gtheme = rc->var == graph_theme;
	    int j, active = 0;

	    l_len++;
	    gtk_table_resize(GTK_TABLE(l_table), l_len, l_cols);
	    w = gtk_label_new(_(rc->description));
	    if (tab == TAB_MAIN) {
		gtk_misc_set_alignment(GTK_MISC(w), 0.0, 0.5);
	    } else {
		gtk_misc_set_alignment(GTK_MISC(w), 1.0, 0.5);
	    }
	    gtk_table_attach_defaults(GTK_TABLE(l_table),
				      w, 0, 1, l_len - 1, l_len);

	    rc->widget = gtk_combo_box_text_new();

	    if (tab == TAB_MAIN) {
		GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

		gtk_box_pack_start(GTK_BOX(hbox), rc->widget, FALSE, FALSE, 5);
		if (gtheme) {
		    add_themes_examples_button(hbox);
		}
		gtk_table_attach(GTK_TABLE(l_table), hbox,
				 1, 2, l_len - 1, l_len,
				 GTK_EXPAND | GTK_FILL, 0, 0, 0);
	    } else {
		gtk_table_attach(GTK_TABLE(l_table), rc->widget,
				 1, 2, l_len - 1, l_len,
				 0, 0, 0, 0);
	    }

	    if (rc->flags & FLOATSET) {
		/* special: graph scale */
		double *xvar = (double *) rc->var;
		char numstr[8];

		for (j=0; j<n_graph_scales; j++) {
		    sprintf(numstr, "%.1f", graph_scales[j]);
		    combo_box_append_text(rc->widget, numstr);
		    if (graph_scales[j] == *xvar) {
			active = j;
		    }
		}
	    } else if (langs) {
		char *strvar = (char *) rc->var;
		const char *str;
		int jj = 0;

		for (j=LANG_AUTO; j<LANG_MAX; j++) {
		    str = gretl_lang_string_from_id(j);
		    if (str != NULL) {
			combo_box_append_text(rc->widget, str);
			if (!strcmp(str, strvar)) {
			    active = jj;
			}
			jj++;
		    }
		}
	    } else {
		char *strvar = NULL;
		int *intvar = NULL;
		const char **strs;
		int nopt;

		if (rc->flags & INTSET) {
		    intvar = (int *) rc->var;
		} else {
		    strvar = (char *) rc->var;
		}

		strs = get_list_setting_strings(rc->var, &nopt);
		for (j=0; j<nopt; j++) {
		    combo_box_append_text(rc->widget, _(strs[j]));
		    if (strvar != NULL && !strcmp(strs[j], strvar)) {
			active = j;
		    } else if (intvar != NULL && j == *intvar) {
			active = j;
		    }
		}
	    }
	    if (tab == TAB_VCV) {
		int ww = get_string_width("XXArellanoXXXXX");

		gtk_widget_set_size_request(rc->widget, ww, -1);
	    }
	    gtk_combo_box_set_active(GTK_COMBO_BOX(rc->widget), active);
	    if (langs) {
		g_signal_connect(G_OBJECT(rc->widget), "changed",
				 G_CALLBACK(try_switch_locale),
				 NULL);
	    } else if (rc->var == &sview_style) {
		GtkWidget *sampler;

		sampler = embed_style_sampler(vbox);
		g_signal_connect(G_OBJECT(rc->widget), "changed",
				 G_CALLBACK(try_switch_style),
				 sampler);
	    }
	} else if (rc->flags & SPINSET) {
	    int *intvar = (int *) rc->var;

	    l_len++;
	    gtk_table_resize(GTK_TABLE(l_table), l_len, l_cols);
	    w = gtk_label_new(_(rc->description));
	    gtk_misc_set_alignment(GTK_MISC(w), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(l_table),
				      w, 0, 1, l_len - 1, l_len);

	    /* for now, this is specific to tab-spaces */
	    rc->widget = gtk_spin_button_new_with_range(2, 8, 1);
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(rc->widget), *intvar);
	    gtk_table_attach_defaults(GTK_TABLE(l_table),
				      rc->widget, 1, 2, l_len - 1, l_len);
	} else if (!(rc->flags & INVISET)) {
	    /* visible string variable */
	    char *strvar = (char *) rc->var;

	    s_len++;
	    gtk_table_resize(GTK_TABLE(s_table), s_len, s_cols);
	    w = gtk_label_new(_(rc->description));
	    gtk_misc_set_alignment(GTK_MISC(w), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(s_table), w,
				      0, 1, s_len - 1, s_len);

	    rc->widget = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(s_table),
				      rc->widget, 1, 2, s_len - 1, s_len);
	    gtk_entry_set_text(GTK_ENTRY(rc->widget), strvar);

	    if (rc->flags & BROWSER) {
		/* add path browse button */
		w = make_path_browse_button(rc, notebook);
		table_attach_fixed(GTK_TABLE(s_table), w,
				   2, 3, s_len - 1, s_len);
	    }

	    if (rc->flags & FIXSET) {
		gtk_widget_set_sensitive(rc->widget, FALSE);
		gtk_widget_set_sensitive(w, FALSE);
	    }
	}
    }

#ifdef DEVEL_OPTS
    if (tab == TAB_EDITOR) {
	/* additional material for the Editor tab if we're making gretl_edit */
        GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
        GtkWidget *label, *e_table;

	label = gtk_label_new("Developer options:");
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

	hbox = gtk_hbox_new(FALSE, 5);
	e_table = gtk_table_new(2, 2, FALSE);
	gtk_table_set_row_spacings(GTK_TABLE(e_table), 10);
	gtk_table_set_col_spacings(GTK_TABLE(e_table), 10);
	gtk_box_pack_start(GTK_BOX(hbox), e_table, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

        label = gtk_label_new("gretlcli");
        gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
        gtk_table_attach_defaults(GTK_TABLE(e_table), label,
                                  0, 1, 0, 1);
        cli_selector = combo_box_text_new_with_entry();
        populate_gretlcli_path_combo(cli_selector);
        gtk_table_attach_defaults(GTK_TABLE(e_table), cli_selector,
                                  1, 2, 0, 1);

        label = gtk_label_new("environment");
        gtk_misc_set_alignment(GTK_MISC(label), 1.0, 0.5);
        gtk_table_attach_defaults(GTK_TABLE(e_table), label,
                                  0, 1, 1, 2);
        env_selector = combo_box_text_new_with_entry();
        populate_gretlcli_env_combo(env_selector);
        gtk_table_attach_defaults(GTK_TABLE(e_table), env_selector,
                                  1, 2, 1, 2);
        gtk_widget_show_all(e_table);
    }
#endif
}

static void set_gp_colors (void)
{
    const char *s = gpcolors;
    char cstr[2][12];

    *cstr[0] = *cstr[1] = '\0';

    if (sscanf(s, "%10s %10s", cstr[0], cstr[1]) == 2) {
	set_graph_color_from_string(0, cstr[0]);
	set_graph_color_from_string(1, cstr[1]);
    }
}

static void set_gp_scale (void)
{
    set_default_png_scale(graph_scale);
}

static void set_gp_theme (void)
{
    set_plotstyle(graph_theme);
}

double next_graph_scale (double s, int mod)
{
    int i;

    for (i=0; i<n_graph_scales; i++) {
	if (s == graph_scales[i]) {
	    if (mod > 0) {
		return i < n_graph_scales - 1 ? graph_scales[i+1] : NADBL;
	    } else {
		return i > 0 ? graph_scales[i-1] : NADBL;
	    }
	}
    }

    return NADBL;
}

double min_graph_scale (void)
{
    return graph_scales[0];
}

double max_graph_scale (void)
{
    return graph_scales[n_graph_scales-1];
}

#if (defined(HAVE_TRAMO) || defined(HAVE_X12A)) && !defined(GRETL_EDIT)

static void maybe_revise_tramo_x12a_status (void)
{
# ifdef HAVE_TRAMO
    if (strcmp(paths.tramo, gretl_tramo())) {
	set_tramo_status();
    }
# endif

# ifdef HAVE_X12A
    if (strcmp(paths.x12a, gretl_x12_arima())) {
	set_x12a_status();
    }
# endif
}

#endif

static void flag_changed (RCVAR *rcvar, int *changed)
{
    if (rcvar->flags & RESTART) {
	*changed = 2;
    } else if (*changed == 0) {
	*changed = 1;
    }
}

static void rcvar_set_int (RCVAR *rcvar, int ival, int *changed)
{
    int *intvar = (int *) rcvar->var;

    if (ival != *intvar) {
	flag_changed(rcvar, changed);
	*intvar = ival;
	if (intvar == &lcnumeric) {
	    int langid = gretl_lang_id_from_name(langpref);

#ifdef G_OS_WIN32
            if (langid == LANG_AUTO && !getenv("LANG")) {
                langid = win32_lang_id_from_locale();
            }
#endif
	    set_lcnumeric(langid, lcnumeric);
	}
    }
}

static int blank_ok (RCVAR *rcvar)
{
    if (!strcmp(rcvar->key, "mpi_hosts")) {
	return 1;
    } else if (!strcmp(rcvar->key, "dbproxy")) {
	return 1;
    } else {
	return 0;
    }
}

static void rcvar_set_string (RCVAR *rcvar, const char *sval, int *changed)
{
    char *strvar = (char *) rcvar->var;

    if (sval == NULL && blank_ok(rcvar)) {
	/* allow "blanking out" of the item */
	if (*strvar != '\0') {
	    flag_changed(rcvar, changed);
	    *strvar = '\0';
	}
	return;
    }

    if (sval != NULL && *sval != '\0' && strcmp(sval, strvar)) {
	flag_changed(rcvar, changed);
	*strvar = '\0';
	strncat(strvar, sval, rcvar->len - 1);
    }
}

static void rcvar_set_double (RCVAR *rcvar, const char *sval, int *changed)
{
    double *xvar = (double *) rcvar->var;

    if (sval != NULL && *sval != '\0') {
	double xval = atof(sval);

	if (xval != *xvar) {
	    flag_changed(rcvar, changed);
	    *xvar = xval;
	}
    }
}

static void restart_message (void)
{
    infobox(_("This change will take effect when you restart gretl"));
}

/* Register and react to changes from the main Preferences dialog
   or the console-specific preferences.
*/

static void apply_prefs_changes (GtkWidget *widget, GtkWidget *parent)
{
    RCVAR *rcvar;
    GtkWidget *w;
    int changed = 0;
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	rcvar = &rc_vars[i];
	w = rcvar->widget;
	if (w == NULL) {
	    continue;
	}
	if (rcvar->flags & BOOLSET) {
	    int bval = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

	    rcvar_set_int(rcvar, bval, &changed);
	} else if (rcvar->flags & SPINSET) {
	    int ival = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(w));

	    rcvar_set_int(rcvar, ival, &changed);
	} else if (rcvar->flags & (USERSET | MACHSET)) {
	    gchar *sval = entry_box_get_trimmed_text(w);

	    rcvar_set_string(rcvar, sval, &changed);
	    g_free(sval);
	} else if (rcvar->flags & LISTSET) {
	    GtkComboBox *box = GTK_COMBO_BOX(w);

	    if (rcvar->flags & INTSET) {
		int ival = gtk_combo_box_get_active(box);

		rcvar_set_int(rcvar, ival, &changed);
	    } else {
		gchar *sval = combo_box_get_active_text(box);

		if (rcvar->flags & FLOATSET) {
		    rcvar_set_double(rcvar, sval, &changed);
		} else {
		    rcvar_set_string(rcvar, sval, &changed);
		}
		g_free(sval);
	    }
	}
    }

    if (changed > 1) {
	restart_message();
    }

#if (defined(HAVE_TRAMO) || defined(HAVE_X12A)) && !defined(GRETL_EDIT)
    maybe_revise_tramo_x12a_status();
#endif

#ifdef HAVE_MPI
    set_mpi_variant(mpi_pref);
#endif

#ifdef DEVEL_OPTS
    set_gretlcli_path(cli_selector);
    set_gretlcli_env(env_selector);
#endif

    if (use_proxy && *http_proxy == '\0') {
	/* fix inconsistency */
	use_proxy = 0;
    }

    /* update graphing info */
    set_default_png_scale(graph_scale);
    set_plotstyle(graph_theme);

    write_rc(OPT_NONE); /* note: calls gretl_update_paths */

    /* register these settings for the current session using
       the "libset" apparatus
    */
    libset_set_bool(SHELL_OK, shellok);
    libset_set_bool(ROBUST_Z, robust_z);

#ifndef GRETL_EDIT
    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_garch_alt_vcv(hc_garch);
    set_panel_hccme(hc_panel);
    selector_register_hc_choice();
#endif

    if (parent != NULL) {
	windata_t *vwin = window_get_active_vwin(parent);

	if (vwin != NULL && vwin->sbuf != NULL) {
	    /* called from sourceview window or tab, or console */
	    update_script_editor_options(vwin);
	}
    }

    gretl_www_init(http_proxy, use_proxy);
}

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *) b) {
	strcpy(s, "true");
    } else {
	strcpy(s, "false");
    }
}

/* non-static because also called from dialogs.c on exit */

int write_rc (gretlopt opt)
{
    RCVAR *rcvar;
    FILE *fp;
    char val[6];
    char *strvar;
    int i;

    fp = gretl_fopen(rcfile, "w");
    if (fp == NULL) {
	file_write_errbox(rcfile);
	return E_FOPEN;
    }

    fprintf(fp, "# gretl config file\n");
    fputs("# build date\n", fp);
    fprintf(fp, "build_date = %s\n", BUILD_DATE);

    for (i=0; rc_vars[i].var != NULL; i++) {
	rcvar = &rc_vars[i];
#ifdef G_OS_WIN32
	/* let the registry or gretlnet.txt handle this */
	if (!strcmp(rcvar->key, "gretldir")) {
	    continue;
	}
#endif
	if (rcvar->flags & SKIPSET) {
	    /* don't bother writing out an empty entry */
	    strvar = (char *) rcvar->var;
	    if (*strvar == '\0') {
		continue;
	    }
	}
	fprintf(fp, "# %s\n", rcvar->description);
	if (rcvar->flags & BOOLSET) {
	    boolvar_to_str(rcvar->var, val);
	    fprintf(fp, "%s = %s\n", rcvar->key, val);
	} else if (rcvar->flags & INTSET) {
	    fprintf(fp, "%s = %d\n", rcvar->key, *(int *) rcvar->var);
	} else if (rcvar->flags & FLOATSET) {
	    gretl_push_c_numeric_locale();
	    fprintf(fp, "%s = %g\n", rcvar->key, *(double *) rcvar->var);
	    gretl_pop_c_numeric_locale();
	} else {
	    strvar = (char *) rcvar->var;
	    fprintf(fp, "%s = %s\n", rcvar->key, strvar);
	}
    }

#ifndef GRETL_EDIT
    rc_save_file_lists(fp);
#endif
    fclose(fp);

    if (!(opt & OPT_N)) {
	gretl_update_paths(&paths, update_paths_opt);
    }

    return 0;
}

void sync_path_from_lib (const char *path_id)
{
    if (!strcmp(path_id, "tramo")) {
	strcpy(paths.tramo, gretl_tramo());
	write_rc(OPT_N);
	restart_message();
    } else if (!strcmp(path_id, "x12a")) {
	strcpy(paths.x12a, gretl_x12_arima());
	write_rc(OPT_N);
	restart_message();
    } else {
	fprintf(stderr, "sync_path_from_lib: '%s' ??\n", path_id);
    }
}

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
	if (!strcmp(s, "true")) {
	    /* remedy for ex-boolean */
	    *ivar = 1;
	} else {
	    *ivar = 0;
	}
    }
}

static void str_to_double (const char *s, void *b)
{
    double *xvar = (double *) b;

    if (s != NULL && *s != '\0') {
	*xvar = dot_atof(s);
    }
}

#if !defined(G_OS_WIN32) && !defined(__APPLE__)

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

/* The various things we need to do after reading gretl's
   configuration info -- or perhaps after failing to do so.
   We do this only at start-up time. This function is
   called by the config-file reading function, either
   read_gretlrc() or on MS Windows, read_win32_config().
*/

static int common_read_rc_setup (int updated)
{
    int langid = 0;
    int err = 0;

    libset_set_bool(SHELL_OK, shellok);
    libset_set_bool(USE_CWD, usecwd);
    libset_set_bool(ROBUST_Z, robust_z);
    set_gp_colors();
    set_gp_scale();

    set_xsect_hccme(hc_xsect);
    set_tseries_hccme(hc_tseri);
    set_panel_hccme(hc_panel);
    set_garch_alt_vcv(hc_garch);

    err = gretl_set_paths(&paths);
    if (err) {
	/* tell the user, then turn off the special alarm */
	gui_errmsg(err);
	set_gretl_alarm(0);
    }

    gretl_www_init(http_proxy, use_proxy);
    set_tex_use_pdf(latex);
    set_gp_theme();

#if !defined(G_OS_WIN32) && !defined(__APPLE__)
    maybe_fix_viewpdf();
#endif

#ifndef GRETL_EDIT
# ifdef HAVE_TRAMO
    set_tramo_status();
# endif
# ifdef HAVE_X12A
    set_x12a_status();
# endif
#endif

    langid = gretl_lang_id_from_name(langpref);
#if 0
    fprintf(stderr, "rc_setup: langpref='%s', langid=%d, lcnumeric=%d\n",
	    langpref, langid, lcnumeric);
#endif
    force_language(langid);
    if (langid == LANG_C) {
	force_english_help();
    }
    set_lcnumeric(langid, lcnumeric);

    if (updated) {
	update_addons_index(NULL);
    }

    return err;
}

/* On reading from text rc file, look up the key in the "key = value"
   line we just read, and set the value of the associated variable.
*/

static void find_and_set_rc_var (const char *key, const char *val)
{
    RCVAR *rcvar;
    char *strvar;
    int i;

    for (i=0; rc_vars[i].key != NULL; i++) {
	rcvar = &rc_vars[i];
	if (!strcmp(key, rcvar->key)) {
	    if (!(rcvar->flags & FIXSET)) {
		if (rcvar->flags & BOOLSET) {
		    str_to_boolvar(val, rcvar->var);
		    if (!strcmp(key, "collect_plots")) {
			/* special: set to "off" or "auto" */
			int *bvar = (int *) rcvar->var;

			libset_set_int(PLOT_COLLECT, *bvar);
		    }
		} else if (rcvar->flags & INTSET) {
		    str_to_int(val, rcvar->var);
		} else if (rcvar->flags & FLOATSET) {
		    str_to_double(val, rcvar->var);
		} else {
		    strvar = (char *) rcvar->var;
		    *strvar = '\0';
		    strncat(strvar, val, rcvar->len - 1);
		}
		rcvar->flags |= GOTSET;
	    }
	    if (rcvar->var == &swallow_pref) {
		swallow = swallow_pref;
	    }
	    break;
	}
    }
}

#ifdef G_OS_WIN32

static int maybe_get_network_settings (void)
{
    const char *netfile;
    FILE *fp;
    char *strvar;
    int gotnet = 0;

    netfile = get_gretlnet_filename();
    if (netfile == NULL) {
	return 0;
    }

    fp = gretl_fopen(netfile, "r");

    if (fp != NULL) {
	char line[MAXLEN], key[32], linevar[MAXLEN];
	int j, calldrive = tolower(netfile[0]);

	while (fgets(line, MAXLEN, fp)) {
	    int gotvar = 0;
	    char *p = line;

	    while (isspace(*p)) p++;
	    if (*p == '#') continue;

	    if (sscanf(p, "%31s", key) == 1) {
		strcpy(linevar, p + strlen(key) + 3);
		gretl_strstrip(linevar);
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
				gotnet = 0;
				goto network_quit;
			    } else {
				strvar = (char *) rc_vars[j].var;
				*strvar = '\0';
				strncat(strvar, linevar, rc_vars[j].len - 1);
			    }
			}
			rc_vars[j].flags |= FIXSET;
			gotvar = gotnet = 1;
		    }
		    if (gotvar) break;
		}
	    }
	}
    network_quit:
	fclose(fp);
    }

    return gotnet;
}

/* Try reading user settings from .gretl2rc in appdata directory: if
   this succeeds it will pre-empt reading from the registry --
   except that we'll respect the registry entry (or perhaps
   argv[0] at startup) for gretldir.
*/

static void win32_read_gretlrc (int *updated)
{
    char line[MAXLEN], key[32], linevar[MAXLEN];
    FILE *fp;

    if (*rcfile == '\0') {
	/* shouldn't happen */
	return;
    }

    fp = gretl_fopen(rcfile, "r");

#if 0
    fprintf(stderr, "rcfile: '%s' (%s)\n", rcfile,
	    fp == NULL ? "not found" : "found");
#endif

    if (fp == NULL) {
	/* not necessarily an error: may be starting from scratch */
	return;
    }

    while (fgets(line, sizeof line, fp) != NULL) {
	if (line[0] == '#') {
	    continue;
	}
	if (!strncmp(line, "recent", 6)) {
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    /* note: don't take gretldir from here */
	    if (strcmp(key, "gretldir")) {
		strcpy(linevar, line + strlen(key) + 3);
		gretl_strstrip(linevar);
		if (*linevar != '\0') {
		    if (!strcmp(key, "build_date")) {
			*updated = gretl_is_updated(linevar);
		    } else if (!strcmp(key, "userdir")) {
			/* legacy */
			find_and_set_rc_var("workdir", linevar);
		    } else {
			find_and_set_rc_var(key, linevar);
		    }
		}
	    }
	}
    }

#ifndef GRETL_EDIT
    if (!strncmp(line, "recent", 6)) {
	rc_read_file_lists(fp, line);
    }
#endif

    fclose(fp);
}

/* This function is not static since it is called from
   gretl_win32_init() in gretlwin32.c
*/

int read_win32_config (int debug, int ignore_rc)
{
    RCVAR *rcvar;
    char value[MAXSTR];
    char *appdata;
    char *strvar;
    int updated = 0;
    int i, err = 0;

    if (chinese_locale()) {
	strcpy(fixedfontname, "NSimSun 10");
	strcpy(default_fixedfont, "NSimSun 10");
    } else if (japanese_locale()) {
	strcpy(fixedfontname, "MS Gothic 10");
	strcpy(default_fixedfont, "MS Gothic 10");
    }

    rcfile[0] = '\0';

#ifndef PKGBUILD
    /* try "HOME" first */
    if (rcfile[0] == '\0') {
	char *home = getenv("HOME");

	if (home != NULL) {
	    strcpy(rcfile, home);
	    slash_terminate(rcfile);
	    strcat(rcfile, ".gretl2rc");
	}
    }
#endif

    if (rcfile[0] == '\0') {
	appdata = appdata_path();
	if (appdata != NULL) {
	    sprintf(rcfile, "%s\\gretl\\.gretl2rc", appdata);
	    free(appdata);
	}
    }

    /* see if we have a gretlnet.txt in place, and if so,
       read config from it */
    maybe_get_network_settings();

    if (!ignore_rc) {
	/* read from user config file */
	win32_read_gretlrc(&updated);
    }

    /* now read from registry for a few items, if they're
       not already set */

    for (i=0; rc_vars[i].key != NULL; i++) {
	int regerr = 0;

	rcvar = &rc_vars[i];

	if (rcvar->flags & (FIXSET | GOTSET)) {
	    /* already set via gretlnet.txt or user rcfile */
	    continue;
	}

	*value = '\0';

	/* note: by now we have already determined gretldir as
	   best we can; don't overwrite its setting
	*/
	if ((rcvar->flags & MACHSET) && strcmp(rcvar->key, "gretldir")) {
	    regerr = read_reg_val(HKEY_LOCAL_MACHINE,
				  get_reg_base(rcvar->key),
				  rcvar->key,
				  value);
	}

	if (debug && *value != '\0') {
	    fprintf(stderr, "reg: err = %d, '%s' -> '%s'\n", regerr, rcvar->key,
		    value);
	}

	if (!regerr && *value != '\0') {
	    /* replace defaults only if we actually got something */
	    if (rcvar->flags & BOOLSET) {
		str_to_boolvar(value, rcvar->var);
	    } else if (rcvar->flags & INTSET) {
		str_to_int(value, rcvar->var);
	    } else {
		strvar = (char *) rcvar->var;
		*strvar = '\0';
		strncat(strvar, value, rcvar->len - 1);
	    }
	}
    }

    err = common_read_rc_setup(updated);

    if (debug) {
	fprintf(stderr, "read_win32_config: returning %d\n", err);
    }

    return err;
}

#else /* end of win32 version, now plain GTK */

/* Note: even if we fail to read the rc file, we should still do
   common_read_rc_setup(), since that will establish defaults for
   basic paths, etc., and give gretl a chance of running OK.
*/

static int read_gretlrc (void)
{
    FILE *fp = gretl_fopen(rcfile, "r");
    int updated = 0;

    if (fp == NULL) {
	fprintf(stderr, "Couldn't read %s\n", rcfile);
    } else {
	char line[MAXLEN], key[32], linevar[MAXLEN];

	while (fgets(line, sizeof line, fp)) {
	    if (*line == '#') {
		continue;
	    }
	    if (!strncmp(line, "recent", 6)) {
		break;
	    }
	    if (sscanf(line, "%31s", key) == 1) {
		*linevar = '\0';
		strncat(linevar, line + strlen(key) + 3, MAXLEN-1);
		gretl_strstrip(linevar);
		if (*linevar != '\0') {
		    if (!strcmp(key, "userdir")) {
			find_and_set_rc_var("workdir", linevar);
		    } else if (!strcmp(key, "build_date")) {
			updated = gretl_is_updated(linevar);
		    } else {
			find_and_set_rc_var(key, linevar);
		    }
		}
	    }
	}

#ifndef GRETL_EDIT
	if (!strncmp(line, "recent", 6)) {
	    rc_read_file_lists(fp, line);
	}
#endif

	fclose(fp);
    }

    return common_read_rc_setup(updated);
}

#endif /* end of non-Windows versions */

static int fontsel_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "MenuFont")) {
	return APP_FONT_SELECTION;
    } else {
	return FIXED_FONT_SELECTION;
    }
}

/* font selection via GtkFontChooser */

#if HAVE_GTK_FONT_CHOOSER

gboolean latin_font_filter (PangoFontFamily *family,
			    PangoFontFace *face,
			    gpointer data)
{
    const char *facename = pango_font_face_get_face_name(face);

    if (strstr(facename, "Italic") || strstr(facename, "Bold")) {
	return FALSE;
    } else {
	return validate_single_font(family, FONT_FILTER_LATIN);
    }
}

gboolean mono_font_filter (PangoFontFamily *family,
			   PangoFontFace *face,
			   gpointer data)
{
    const char *facename = pango_font_face_get_face_name(face);

    if (strstr(facename, "Italic") || strstr(facename, "Bold") ||
	strstr(facename, "Oblique")) {
	return FALSE;
    } else {
	return validate_single_font(family, FONT_FILTER_LATIN_MONO);
    }
}

static void close_font_chooser (GtkWidget *w, GtkFontChooser *fc)
{
    gretl_font_filter_cleanup();
    gtk_widget_destroy(GTK_WIDGET(fc));
}

static void font_selection_ok (GtkWidget *w, GtkFontChooser *fc)
{
    gchar *fontname = gtk_font_chooser_get_font(fc);

    gtk_widget_hide(GTK_WIDGET(fc));

    if (fontname != NULL && *fontname != '\0') {
	int mono = widget_get_int(fc, "mono");

	if (mono) {
	    set_fixed_font(fontname, 1);
	} else {
	    set_app_font(fontname, 1);
	}
	write_rc(OPT_NONE);
    }

    g_free(fontname);
    gretl_font_filter_cleanup();
    gtk_widget_destroy(GTK_WIDGET(fc));
}

static void font_selection_reset (GtkWidget *w, GtkFontChooser *fc)
{
    int mono = widget_get_int(fc, "mono");

    gtk_widget_hide(GTK_WIDGET(fc));

    if (mono) {
	set_fixed_font(default_fixedfont, 1);
    } else {
	set_app_font(system_appfont, 1);
    }
    write_rc(OPT_NONE);

    gretl_font_filter_cleanup();
    gtk_widget_destroy(GTK_WIDGET(fc));
}

static void chooser_font_selector (GtkAction *action)
{
    static GtkWidget *fc = NULL;
    GtkFontFilterFunc filter;
    GtkWidget *hbox, *button;
    int which = fontsel_code(action);
    char *title = NULL;
    const char *fontname = NULL;
    int err = 0;

    if (fc != NULL) {
	gtk_window_present(GTK_WINDOW(fc));
        return;
    }

    if (which == FIXED_FONT_SELECTION) {
	title = _("Font for gretl output windows");
	filter = (GtkFontFilterFunc) &mono_font_filter;
	fontname = fixedfontname;
    } else if (which == APP_FONT_SELECTION) {
	title = _("Font for menus and labels");
	filter = (GtkFontFilterFunc) &latin_font_filter;
	fontname = appfontname;
    }

    err = gretl_font_filter_init();
    if (err) {
	errbox("Failed to initialize font filter");
	return;
    }

    fc = gtk_font_chooser_dialog_new(title, GTK_WINDOW(mdata->main));
    gtk_font_chooser_set_font(GTK_FONT_CHOOSER(fc), fontname);
    gtk_font_chooser_set_filter_func(GTK_FONT_CHOOSER(fc),
				     filter, NULL, NULL);
    gtk_window_set_position(GTK_WINDOW(fc), GTK_WIN_POS_MOUSE);

    if (which == FIXED_FONT_SELECTION) {
	widget_set_int(fc, "mono", 1);
    }

    g_signal_connect(G_OBJECT(fc), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &fc);

    button = gtk_dialog_get_widget_for_response(GTK_DIALOG(fc),
						GTK_RESPONSE_OK);
    if (button != NULL) {
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(font_selection_ok),
			 fc);
    }

    button = gtk_dialog_get_widget_for_response(GTK_DIALOG(fc),
						GTK_RESPONSE_CANCEL);
    if (button != NULL) {
	g_signal_connect(G_OBJECT(button), "clicked",
			 G_CALLBACK(close_font_chooser),
			 fc);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(fc));
    button = gtk_button_new_with_label(_("Reset to default"));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(font_selection_reset), fc);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_reorder_child(GTK_BOX(hbox), button, 0);

    gtk_widget_show_all(fc);
}

#else /* GTK, not using GtkFontChooser */

static void font_selection_ok (GtkWidget *w, GtkFontselHackDialog *fs)
{
    gchar *fontname;

    fontname = gtk_fontsel_hack_dialog_get_font_name(fs);
    gtk_widget_hide(GTK_WIDGET(fs));

    if (fontname != NULL && *fontname != '\0') {
	int filter = gtk_fontsel_hack_dialog_get_filter(fs);

	if (filter == FONT_HACK_LATIN_MONO) {
	    set_fixed_font(fontname, 1);
	} else if (filter == FONT_HACK_LATIN) {
	    set_app_font(fontname, 1);
	}
	write_rc(OPT_NONE);
    }

    g_free(fontname);
    gtk_widget_destroy(GTK_WIDGET(fs));
}

static void font_selection_reset (GtkWidget *w, GtkFontselHackDialog *fs)
{
    int filter = gtk_fontsel_hack_dialog_get_filter(fs);

    gtk_widget_hide(GTK_WIDGET(fs));

    if (filter == FONT_HACK_LATIN_MONO) {
	set_fixed_font(default_fixedfont, 1);
    } else if (filter == FONT_HACK_LATIN) {
	set_app_font(system_appfont, 1);
    }
    write_rc(OPT_NONE);

    gtk_widget_destroy(GTK_WIDGET(fs));
}

static void gtk2_font_selector (GtkAction *action)
{
    static GtkWidget *fontsel = NULL;
    GtkWidget *hbox, *button;
    int filter, which = fontsel_code(action);
    char *title = NULL;
    const char *fontname = NULL;

    if (fontsel != NULL) {
	gtk_window_present(GTK_WINDOW(fontsel));
        return;
    }

    if (which == FIXED_FONT_SELECTION) {
	title = _("Font for gretl output windows");
	filter = FONT_HACK_LATIN_MONO;
	fontname = fixedfontname;
    } else if (which == APP_FONT_SELECTION) {
	title = _("Font for menus and labels");
	filter = FONT_HACK_LATIN;
	fontname = appfontname;
    }

    fontsel = gtk_fontsel_hack_dialog_new(title);

    gtk_fontsel_hack_dialog_set_filter(GTK_FONTSEL_HACK_DIALOG(fontsel),
				       filter);
    gtk_fontsel_hack_dialog_set_font_name(GTK_FONTSEL_HACK_DIALOG(fontsel),
					  fontname);

    gtk_window_set_position(GTK_WINDOW(fontsel), GTK_WIN_POS_MOUSE);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(fontsel));
    button = gtk_button_new_with_label(_("Reset to default"));
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(font_selection_reset),
		     fontsel);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_reorder_child(GTK_BOX(hbox), button, 0);
    gtk_widget_show(button);
    gtk_widget_show(hbox);

    g_signal_connect(G_OBJECT(fontsel), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &fontsel);
    g_signal_connect(G_OBJECT(gtk_fontsel_hack_dialog_ok_button(fontsel)),
		     "clicked", G_CALLBACK(font_selection_ok),
		     fontsel);
    g_signal_connect_swapped(G_OBJECT(gtk_fontsel_hack_dialog_cancel_button(fontsel)),
			     "clicked", G_CALLBACK(gtk_widget_destroy),
			     fontsel);

    gtk_widget_show(fontsel);
}

#endif /* end font-selection dialog variants */

void font_selector (GtkAction *action)
{
#if HAVE_GTK_FONT_CHOOSER
    chooser_font_selector(action);
#else
    gtk2_font_selector(action);
#endif
}

static void impose_font_scale (int scale, int remember)
{
    char fontname[64];

    strcpy(fontname, fixedfontname);
    fontname_set_size(fontname, scale);
    set_fixed_font(fontname, remember);

#ifdef G_OS_WIN32
    if (*appfontname == '\0') {
	get_default_windows_app_font(appfontname);
    }
#endif

    strcpy(fontname, appfontname);
    fontname_set_size(fontname, scale);
    set_app_font(fontname, remember);
}

#define FSCALE_DEFAULT 99

static void set_fscale_default (GtkWidget *w, int *resp)
{
    GtkWidget *dlg = g_object_get_data(G_OBJECT(w), "dlg");

    *resp = FSCALE_DEFAULT;
    gtk_widget_destroy(dlg);
}

void font_scale_selector (GtkAction *action)
{
    const char *opt = N_("Remember this setting");
    GtkWidget *dlg, *hbox, *button;
    int fscale = tmpfontscale;
    int remember = 0;
    int resp = GRETL_CANCEL;

    if (fscale == 0) {
	fscale = fontname_get_size(fixedfontname);
    }

    dlg = build_checks_dialog(_("gretl: font scale"), NULL,
			      &opt,
			      1, &remember,
			      0, 0,
			      0, NULL,
			      &fscale, _("Scale for monospaced and menu fonts"),
			      8, 24,
			      0, mdata->main,
			      &resp);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));
    button = gtk_button_new_with_label(_("Reset to default"));
    g_object_set_data(G_OBJECT(button), "dlg", dlg);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_fscale_default), &resp);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_box_reorder_child(GTK_BOX(hbox), button, 0);

    gtk_widget_show_all(dlg);

    if (resp == GRETL_CANCEL) {
	return;
    } else if (resp == FSCALE_DEFAULT) {
	set_fixed_font(default_fixedfont, 1);
	set_app_font(system_appfont, 1);
	write_rc(OPT_NONE);
    } else if (fscale > 0) {
	impose_font_scale(fscale, remember);
	tmpfontscale = fscale;
    }
}

void update_persistent_graph_colors (void)
{
    print_palette_string(gpcolors);
}

void dump_rc (void)
{
    char dumper[MAXLEN];
    const char *hname;
    GDir *test;
    FILE *fp;
    char *tmp;
    char val[6];
    int i;

    sprintf(dumper, "%sconfig-dump.txt", gretl_workdir());

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

    hname = gretl_home();
    test = gretl_opendir(hname);

    if (test != NULL) {
	fprintf(fp, "Directory '%s' exists, OK\n", hname);
	g_dir_close(test);
    } else {
	fprintf(fp, "Directory '%s' does not exist\n", hname);
    }

    printf("Config info written to %s\n", dumper);

    fclose(fp);
}

#ifndef GRETL_EDIT

int gui_set_working_dir (char *dirname)
{
    int err = gretl_set_path_by_name("workdir", dirname);

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_WDIR, dirname);
    } else {
	mkfilelist(FILE_LIST_WDIR, dirname, 0);
	set_workdir_label();
    }

    return err;
}

struct wdir_setter {
    GtkWidget *dialog;
    GtkWidget *wdir_combo;
    GtkWidget *cwd_radio;
    GtkWidget *keep_radio;
    GtkWidget *ok_button;
};

/* callback from the file selector */

void set_working_dir_callback (GtkWidget *w, char *path)
{
    set_combo_box_default_text(GTK_COMBO_BOX(w), path);
}

static void wdir_browse_callback (GtkWidget *w, struct wdir_setter *wset)
{
    GtkWidget *combo = wset->wdir_combo;

    file_selector_with_parent(SET_WDIR, FSEL_DATA_MISC, combo,
			      wset->dialog);
}

static void open_wdir (GtkButton *b, gpointer p)
{
#if defined(G_OS_WIN32)
    win32_open_file(gretl_workdir());
#elif defined(__APPLE__)
    osx_open_file(gretl_workdir());
#else
    gretl_fork("xdg-open", gretl_workdir(), NULL);
#endif
}

static void
add_wdir_content (GtkWidget *dialog, struct wdir_setter *wset)
{
    GtkWidget *hbox, *vbox, *w = NULL;
    GtkWidget *entry;
    GSList *group = NULL;
    GList *list = NULL;
    gchar *deflt;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    list = get_working_dir_list();

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Working directory:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);

    deflt = g_strdup(gretl_workdir());
    trim_slash(deflt);

    /* combo + browse button for current working dir */
    w = combo_box_text_new_with_entry();
    gtk_container_add(GTK_CONTAINER(hbox), w);
    set_combo_box_strings_from_list(w, list);
    if (deflt != NULL) {
	set_combo_box_default_text(GTK_COMBO_BOX(w), deflt);
    }
    entry = gtk_bin_get_child(GTK_BIN(w));
    gtk_entry_set_width_chars(GTK_ENTRY(entry), 32);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
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

    /* radio 2 for "next time" */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_radio_button_new_with_label(group, _("the current directory "
						 "as determined via the shell"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), usecwd);
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    wset->cwd_radio = w;

    vbox_add_hsep(vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("The file selection dialog should:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, 0, 0, 5);

    /* radio 1 for "remember folder" */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_radio_button_new_with_label(NULL,
					_("remember the last-opened folder"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(w));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    wset->keep_radio = w;

    /* radio 2 for "remember folder" */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_radio_button_new_with_label(group, _("always start in the "
						 "working directory"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), !keep_folder);
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);

    /* button to open working directory via OS */
    vbox_add_hsep(vbox);
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_button_new_with_label(_("Open working directory"));
    gtk_box_pack_start(GTK_BOX(hbox), w, 0, 0, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    g_signal_connect(w, "clicked", G_CALLBACK(open_wdir), NULL);

    g_list_free(list);
    g_free(deflt);
}

static void
apply_wdir_changes (GtkWidget *w, struct wdir_setter *wset)
{
    char tmp[MAXLEN];
    gchar *str;
    int err;

    str = combo_box_get_active_text(GTK_COMBO_BOX(wset->wdir_combo));
    *tmp = '\0';
    if (str != NULL) {
	strncat(tmp, str, MAXLEN - 2);
	g_free(str);
    }

    err = gretl_set_path_by_name("workdir", tmp);

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_WDIR, tmp);
    } else {
	/* sync with "local copy" */
	strcpy(paths.workdir, gretl_workdir());
	mkfilelist(FILE_LIST_WDIR, tmp, 0);
    }

    usecwd = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(wset->cwd_radio));
    keep_folder = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(wset->keep_radio));

    if (!err) {
	if (w == wset->ok_button) {
	    gtk_widget_destroy(wset->dialog);
	}
	set_workdir_label();
    }
}

static void workdir_dialog (int from_wlabel)
{
    static GtkWidget *dialog;
    struct wdir_setter wset;
    GtkWidget *button;
    GtkWidget *hbox, *vbox;

    if (dialog != NULL) {
	gtk_window_present(GTK_WINDOW(dialog));
	return;
    }

    dialog = gretl_dialog_new(_("gretl: working directory"),
			      mdata->main, GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_set_spacing(GTK_BOX(vbox), 5);
    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(gtk_widget_destroyed),
		     &dialog);

    wset.dialog = dialog;
    add_wdir_content(dialog, &wset);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    button = apply_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(apply_wdir_changes), &wset);

    button = cancel_button(hbox);
    g_signal_connect_swapped(G_OBJECT(button), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);

    wset.ok_button = button = ok_button(hbox);
    gtk_widget_grab_default(button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(apply_wdir_changes), &wset);

    context_help_button(hbox, WORKDIR);

    gtk_widget_show_all(dialog);
}

void workdir_dialog0 (void)
{
    workdir_dialog(0);
}

void workdir_dialog1 (void)
{
    workdir_dialog(1);
}

#endif /* not GRETL_EDIT */

#if defined(MAC_THEMING) && GTK_MAJOR_VERSION < 3

void set_up_mac_look (void)
{
    if (!strcmp(themepref, "Lion-like")) {
        /* 2024-12-19: fallback for broken GTK2 theme */
        strcpy(themepref, "Adwaita");
    }

    if (!strncmp(themepref, "Adwaita", 7) ||
	!strcmp(themepref, "Clearlooks")) {
	char *topdir = getenv("GTK_DATA_PREFIX");
	gchar *gtkrc;

	if (topdir != NULL) {
	    gtkrc = g_strdup_printf("%s/share/themes/%s/gtk-2.0/gtkrc",
				    topdir, themepref);
	    gtk_rc_parse(gtkrc);
	    g_free(gtkrc);
	} else {
#if defined(GTK_PREFIX)
	    /* go with the build-time prefix */
	    gtkrc = g_strdup_printf("%s/share/themes/%s/gtk-2.0/gtkrc",
				    GTK_PREFIX, themepref);
#else
	    /* hard-wired? */
	    const char *path = "/Library/Frameworks/gretl-dev.framework/Resources";

	    gtkrc = g_strdup_printf("%s/share/themes/%s/gtk-2.0/gtkrc",
				    path, themepref);
#endif
	    gtk_rc_parse(gtkrc);
	    g_free(gtkrc);
	}
    }
}

#elif defined(MAC_THEMING) && GTK_MAJOR_VERSION == 3

void set_up_mac_look (void)
{
    GtkSettings *settings = gtk_settings_get_default();

    if (strstr(themepref, "dark") != NULL) {
        g_object_set(G_OBJECT(settings), "gtk-application-prefer-dark-theme",
                     TRUE, NULL);
    }
}

#elif defined (G_OS_WIN32) && GTK_MAJOR_VERSION < 3

void set_up_windows_look (void)
{
    if (!strcmp(themepref, "Windows-10") ||
        !strcmp(themepref, "Windows-10-Dark") ||
	!strcmp(themepref, "MS-Windows") ||
	!strcmp(themepref, "Clearlooks")) {
	const char *prefix;
	char sl[2] = {0};
	gchar *gtkrc;
	size_t n;
	int needslash;

# ifdef PKGBUILD
	prefix = gretl_home();
# else
	prefix = GTK_PREFIX; /* defined at build-time */
# endif
	n = strlen(prefix);
	needslash = (prefix[n-1] != '\\' && prefix[n-1] != '/');
	sl[0] = strchr(prefix, '/') ? '/' : '\\';
	gtkrc = g_strdup_printf("%s%sshare%sthemes%s%s%sgtk-2.0%sgtkrc",
				prefix, (needslash)? sl : "", sl, sl,
				themepref, sl, sl);
	fprintf(stderr, "gtkrc = '%s'\n", gtkrc);
	gtk_rc_parse(gtkrc);
	g_free(gtkrc);
    } else {
	GtkSettings *settings = gtk_settings_get_default();

	g_object_set(G_OBJECT(settings), "gtk-theme-name", "Raleigh", NULL);
    }
}

void set_wimp_preferred (int s)
{
    if (s) {
	strcpy(themepref, "Windows-10");
    } else {
	strcpy(themepref, "Clearlooks");
    }
}

#elif defined (G_OS_WIN32) && GTK_MAJOR_VERSION == 3

void set_up_windows_look (void)
{
    GtkSettings *settings = gtk_settings_get_default();
    const char *theme_name;

    if (!strncmp(themepref, "Windows-10", 10)) {
        g_setenv("GTK_CSD", "0", 1);
    } else {
        g_setenv("GTK_CSD", "1", 1);
    }

    if (!strcmp(themepref, "Windows 7")) {
        theme_name = "win32";
    } else {
        theme_name = themepref;
    }

    g_object_set(G_OBJECT(settings), "gtk-theme-name", theme_name, NULL);
}

#endif
