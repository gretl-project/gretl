/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* gretl.c : main for gretl */

#include "gretl.h"
#include "treeutils.h"
#include "ssheet.h"
#include "gpt_control.h"

#include <dirent.h>

#ifndef G_OS_WIN32
# include <unistd.h>
# include <signal.h>
# include "../pixmaps/gretl.xpm"  /* program icon for X */
#else
# include <windows.h>
#endif

/* pixmaps for gretl toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.edit.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.netscape.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.ofolder.xpm"
#ifndef GNUPLOT_PNG
# include "../pixmaps/mini.camera.xpm"
#endif

/* functions from other gretl GUI files */

extern void free_modelspec (void);    /* lib.c */
extern void stats_calculator (gpointer data, guint ttest, GtkWidget *widget);
extern void bool_subsample (gpointer data, guint dropmiss, GtkWidget *w);
extern void free_command_stack (void);
extern void open_named_db_list (char *dbname);
extern void open_named_remote_db_list (char *dbname);
extern void gui_set_panel_structure (gpointer data, guint u, GtkWidget *w);
extern void time_series_dialog (gpointer data, guint u, GtkWidget *w);

/* functions private to gretl.c */

static void populate_list (GtkWidget *widget, DATAINFO *datainfo);
static void sort_varlist (gpointer p, guint col, GtkWidget *w);
static void make_toolbar (GtkWidget *w, GtkWidget *box);
static void clip_init (GtkWidget *w);
static GtkWidget *make_main_window (int gui_get_data);

static void build_var_popup (windata_t *win);
static void build_selection_popup (void);
static gint var_popup_click (GtkWidget *widget, gpointer data);
static gboolean main_popup_handler (GtkWidget *widget, GdkEvent *event,
				    gpointer data);
static void mdata_edit (gpointer data, guint colnum, GtkWidget *w);

static void check_for_extra_data (void);
static void set_up_main_menu (void);
static void startR (gpointer p, guint opt, GtkWidget *w);
static void Rcleanup (void);
static void auto_store (void);

GtkWidget *toolbar_box = NULL; /* shared with settings.c */

static GtkWidget *datalabel;
static GtkWidget *main_vbox;
static GtkWidget *gretl_toolbar;
static GtkWidget *selection_popup;
GtkTooltips *gretl_tips;

int *default_list = NULL;

static GtkTargetEntry target_table[] = {
    {"text/uri-list", 0, 1},
};

static void  
drag_data_received  (GtkWidget          *widget,
		     GdkDragContext     *dc,
		     gint                x,
		     gint                y,
		     GtkSelectionData   *data,
		     guint               info,
		     guint               time,
		     gpointer            p);

#ifdef USE_GNOME
static char *optrun, *optdb;

static const struct poptOption options[] = {
    {"run", 'r', POPT_ARG_STRING, &optrun, 0, 
     N_("open a script file on startup"), "SCRIPT"},
    {"db", 'd', POPT_ARG_STRING, &optdb, 0, 
     N_("open a database on startup"), "DATABASE"},
    {"webdb", 'w', POPT_ARG_STRING, &optdb, 0, 
     N_("open a remote (web) database on startup"), "REMOTE_DB"},
    {NULL, '\0', 0, NULL, 0}
};
#endif /* USE_GNOME */

windata_t *mdata;
DATAINFO *datainfo;
DATAINFO *subinfo;
DATAINFO *fullinfo;
char *errtext;
char cmdfile[MAXLEN], scriptfile[MAXLEN];
char trydatfile[MAXLEN], tryscript[MAXLEN];
char line[1024];
PATHS paths;                /* useful paths */
CMD command;                /* gretl command struct */
double **Z;                 /* data set */
double **subZ;              /* sub-sampled data set */
double **fullZ;             /* convenience pointer */
MODEL **models;             /* gretl models structs */
SESSION session;            /* hold models, graphs */
SESSIONBUILD rebuild;       /* rebuild session later */

int plot_count, data_status, orig_vars;
PRN *cmds;
gchar *clipboard_buf; /* for copying models as LaTeX */

/* Is NLS translation in effect? */
int nls_on;

/* defaults for some options */
int expert = FALSE; 
int updater = FALSE;
int want_toolbar = TRUE;
char dbproxy[21];

#ifdef G_OS_WIN32
char Rcommand[MAXSTR] = "RGui.exe";
char editor[MAXSTR] = "winword.exe";
char calculator[MAXSTR] = "calc.exe";
char viewdvi[MAXSTR] = "windvi.exe";
#else
char editor[MAXSTR] = "emacs";
char calculator[MAXSTR] = "xcalc";
char viewdvi[MAXSTR] = "xdvi";
# ifdef USE_GNOME
char Rcommand[MAXSTR] = "R --gui=gnome";
extern const char *version_string;
# else
char Rcommand[MAXSTR] = "xterm -e R";
# endif
#endif

static void spreadsheet_edit (gpointer p, guint u, GtkWidget *w) 
{
    show_spreadsheet(NULL);
}

#if defined(USE_GNOME)

static void gnome_help (void)
{
    GError *error = NULL;

    gnome_help_display ("gretl.xml", NULL, &error);
        
    if (error != NULL) {
	g_warning (error->message);
	g_error_free (error);
    }
}

#elif defined(G_OS_WIN32)
static void win_help (void)
{
    char hlpfile[MAXLEN];

    sprintf(hlpfile, "hh.exe \"%s\\gretl.chm\"", paths.gretldir);
    if (WinExec(hlpfile, SW_SHOWNORMAL) < 32)
        errbox(_("Couldn't access help file"));
}

static int unmangle (const char *dosname, char *longname);
#endif

extern void find_var (gpointer p, guint u, GtkWidget *w); /* gui_utils.c */

GtkItemFactoryEntry data_items[] = {

    /* File menu */
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },

    /* File, Open data */
    { N_("/File/_Open data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open data/user file..."), NULL, open_data, OPEN_DATA, NULL },
    { N_("/File/Open data/sample file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open data/sample file/Ramanathan..."), NULL, 
      display_files, RAMU_DATA, NULL },
    { N_("/File/Open data/sample file/Greene..."), NULL, 
      display_files, GREENE_DATA, NULL },
    { N_("/File/Open data/sample file/Wooldridge..."), NULL,
      open_data, OPEN_DES, NULL },
    { N_("/File/Open data/sample file/Penn World Table..."), NULL, 
      display_files, PWT_DATA, NULL },
    { N_("/File/Open data/sep1"), NULL, NULL, 0, "<Separator>" },    
    { N_("/File/Open data/import CSV..."), NULL, open_data, OPEN_CSV, NULL },
    { N_("/File/Open data/import BOX..."), NULL, open_data, OPEN_BOX, NULL },
    { N_("/File/Open data/import Gnumeric..."), NULL, open_data, 
      OPEN_GNUMERIC, NULL },
    { N_("/File/Open data/import Excel..."), NULL, open_data, 
      OPEN_EXCEL, NULL },

    /* File, Append data */
    { N_("/File/_Append data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Append data/from CSV..."), NULL, open_data, APPEND_CSV, NULL },
    { N_("/File/Append data/from Gnumeric..."), NULL, open_data, 
      APPEND_GNUMERIC, NULL },
    { N_("/File/Append data/from Excel..."), NULL, open_data, 
      APPEND_EXCEL, NULL },

    /* File, Save data */
    { N_("/File/_Save data"), NULL, auto_store, 0, NULL },
    { N_("/File/Save data _as"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Save data as/_standard format..."), NULL, file_save, 
      SAVE_DATA, NULL },
    { N_("/File/Save data as/_gzipped..."), NULL, 
      file_save, SAVE_GZDATA, NULL },
#ifdef notdef
    { N_("/File/Save data as/_alternative formats/_single precision binary..."), 
      NULL, file_save, SAVE_BIN1, NULL },
    { N_("/File/Save data as/_alternative formats/_double precision binary..."),
      NULL, file_save, SAVE_BIN2, NULL },
#endif 
    /* File, Export data */
    { N_("/File/_Export data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Export data/_CSV..."), NULL, file_save, EXPORT_CSV, NULL },
    { N_("/File/Export data/GNU _R..."), NULL, file_save, EXPORT_R, NULL },
    { N_("/File/Export data/GNU _octave..."), NULL, file_save, 
      EXPORT_OCTAVE, NULL },
    { N_("/File/C_lear data set"), NULL, verify_clear_data, 0, NULL },
    { N_("/File/sep0"), NULL, NULL, 0, "<Separator>" },

    /* File, Browse databases */
    { N_("/File/_Browse databases"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Browse databases/_gretl native"), NULL, display_files, 
      NATIVE_DB, NULL },
    { N_("/File/Browse databases/_RATS 4"), NULL, display_files, 
      RATS_DB, NULL },
    { N_("/File/Browse databases/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/File/Browse databases/on database _server"), NULL, display_files, 
      REMOTE_DB, NULL },

    /* File, Create dataset */
    { N_("/File/_Create data set"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Create data set/time-series"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Create data set/time-series/annual"), 
      NULL, newdata_dialog, 1, NULL },    
    { N_("/File/Create data set/time-series/quarterly"), 
      NULL, newdata_dialog, 4, NULL },    
    { N_("/File/Create data set/time-series/monthly"), 
      NULL, newdata_dialog, 12, NULL },
    { N_("/File/Create data set/time-series/high frequency"), 
      NULL, NULL, 0, "<Branch>" },
    { N_("/File/Create data set/time-series/high frequency/weekly"), 
      NULL, newdata_dialog, 52, NULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (5-day week)"), 
      NULL, newdata_dialog, 5, NULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (7-day week)"), 
      NULL, newdata_dialog, 7, NULL }, 
    { N_("/File/Create data set/time-series/high frequency/hourly"), 
      NULL, newdata_dialog, 24, NULL }, 
    { N_("/File/Create data set/cross-sectional"), 
      NULL, newdata_dialog, 0, NULL }, 
#ifdef notdef  
    { N_("/File/Create data set/panel"), 
      NULL, start_panel_dialog, 0, NULL }, 
#endif 
    { N_("/File/Create data set/simulation"), NULL, gretl_callback, 
      NULLDATA, NULL },
    { N_("/File/sep1"), NULL, NULL, 0, "<Separator>" },
#ifndef GNUPLOT_PNG
    { N_("/File/Save last graph"), NULL, gpt_save_dialog, 0, NULL }, 
    { N_("/File/sep2"), NULL, NULL, 0, "<Separator>" },
#endif
    { N_("/File/_View command log"), NULL, view_log, 0, NULL },
    { N_("/File/sep2a"), NULL, NULL, 0, "<Separator>" },

    /* File, command files */
    { N_("/File/Open command file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open command file/user file..."), NULL, open_script, 
      OPEN_SCRIPT, NULL },
    { N_("/File/Open command file/practice file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open command file/practice file/Ramanathan..."), NULL, 
      display_files, RAMU_PS, NULL },
    { N_("/File/Open command file/practice file/Greene..."), NULL, 
      display_files, GREENE_PS, NULL },
    { N_("/File/Open command file/practice file/Penn World Table..."), NULL, 
      display_files, PWT_PS, NULL },
    { N_("/File/New command file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/New command file/regular script"), NULL, do_new_script, 0, NULL },
    { N_("/File/New command file/Monte Carlo loop"), NULL, 
      do_new_script, 1, NULL },
    { N_("/File/sep3"), NULL, NULL, 0, "<Separator>" },

    /* File, preferences */
    { N_("/File/_Preferences"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/_Preferences/_General..."), NULL, options_dialog, 0, NULL },
    { N_("/File/Preferences/_Fixed font..."), NULL, font_selector, 1, NULL },
#ifndef USE_GNOME
    { N_("/File/Preferences/_Menu font..."), NULL, font_selector, 0, NULL },
#endif
    { N_("/File/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/File/E_xit"), NULL, menu_exit_check, 0, NULL },

    /* Utilities menu */
    { N_("/_Utilities"), NULL, NULL, 0, "<Branch>" },
    { N_("/Utilities/Statistical tables"), NULL, stats_calculator, 1, NULL },
    { N_("/Utilities/p-value finder"), NULL, stats_calculator, 0, NULL },
    { N_("/Utilities/Test statistic calculator"), NULL, stats_calculator, 2, NULL },
    { N_("/Utilities/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Gretl console"), NULL, console, 0, NULL },
    { N_("/Utilities/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Start GNU R"), NULL, startR, 0, NULL },

    /* Session menu */
    { N_("/_Session"), NULL, NULL, 0, "<Branch>" },
    { N_("/Session/_Icon view"), NULL, view_session, 0, NULL },
#ifndef GNUPLOT_PNG
    { N_("/Session/_Add last graph"), NULL, add_last_graph, 0, NULL },
#endif
    { N_("/Session/sep0"), NULL, NULL, 0, "<Separator>" },
    { N_("/Session/_Open..."), NULL, open_script, OPEN_SESSION, NULL },
    { N_("/Session/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Session/_Save"), NULL, save_session_callback, SAVE_AS_IS, NULL },
    { N_("/Session/Save _as..."), NULL, save_session_callback, SAVE_RENAME, NULL },
#ifdef not_yet
    { N_("/Session/_Delete"), NULL, delete_session_callback, 0, NULL },
#endif

    /* Data menu */
    { N_("/_Data"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/_Display values"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Display values/_all variables"), NULL, display_data, 0, NULL },
    { N_("/Data/Display values/_selected variables..."), 
      NULL, display_selected, 0, NULL },
    { N_("/Data/_Edit values"), NULL, spreadsheet_edit, 0, NULL },
    { "/Data/sepsort", NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Sort variables"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Sort variables/by ID number"), NULL, sort_varlist, 0, NULL },
    { N_("/Data/Sort variables/by name"), NULL, sort_varlist, 1, NULL },
    { N_("/Data/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/_Graph specified vars"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Graph specified vars/Time-series plot..."), 
      NULL, selector_callback, GR_PLOT, NULL },
    { N_("/Data/Graph specified vars/X-Y scatter..."), 
      NULL, selector_callback, GR_XY, NULL },
    { N_("/Data/Graph specified vars/X-Y with impulses..."), 
      NULL, selector_callback, GR_IMP, NULL },
    { N_("/Data/Graph specified vars/X-Y with factor separation..."), 
      NULL, selector_callback, GR_DUMMY, NULL },
    { N_("/Data/Graph specified vars/Boxplots..."), 
      NULL, gretl_callback, GR_BOX, NULL },
    { N_("/Data/Graph specified vars/Notched boxplots..."), 
      NULL, gretl_callback, GR_NBOX, NULL },
    { N_("/Data/_Multiple scatterplots..."), 
      NULL, selector_callback, SCATTERS, NULL},
    { N_("/Data/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/_Read info"), NULL, open_info, 0, NULL },
    { N_("/Data/Edit _info"), NULL, edit_header, 0, NULL },
    { N_("/Data/Print description"), NULL, print_report, 0, NULL },
    { N_("/Data/sep3"), NULL, NULL, 0, "<Separator>" },

    { N_("/Data/_Summary statistics"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/_Summary statistics/_all variables"), NULL, 
      do_menu_op, SUMMARY, NULL },
    { N_("/Data/_Summary statistics/_selected variables"), NULL, 
      do_menu_op, SUMMARY_SELECTED, NULL },

    { N_("/Data/_Correlation matrix"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/_Correlation matrix/_all variables"), 
      NULL, do_menu_op, CORR, NULL },
    { N_("/Data/_Correlation matrix/_selected variables"), NULL, do_menu_op, 
      CORR_SELECTED, NULL },

    { N_("/Data/sep4"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Difference of means"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Difference of means/assuming equal variances..."), NULL, 
      gretl_callback, MEANTEST, NULL },
    { N_("/Data/Difference of means/assuming unequal variances..."), NULL, 
      gretl_callback, MEANTEST2, NULL },
    { N_("/Data/Difference of variances..."), NULL, gretl_callback, VARTEST, NULL },
    { N_("/Data/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Add variables/time trend"), NULL, add_time, 0, NULL },
    { N_("/Data/Add variables/index variable"), NULL, add_time, 1, NULL },
    { N_("/Data/Add variables/logs of selected variables"), NULL, 
      add_logs_etc, LOGS, NULL },
    { N_("/Data/Add variables/lags of selected variables"), NULL, 
      add_logs_etc, LAGS, NULL },
    { N_("/Data/Add variables/squares of selected variables"), NULL, 
      add_logs_etc, SQUARE, NULL },
    { N_("/Data/Add variables/first differences of selected variables"), NULL, 
      add_logs_etc, DIFF, NULL },
    { N_("/Data/Add variables/log differences of selected variables"), NULL, 
      add_logs_etc, LDIFF, NULL },
    { N_("/Data/Add variables/periodic dummies"), NULL, add_dummies, 0, NULL },
    { N_("/Data/Add variables/panel dummies"), NULL, add_dummies, 1, NULL },
    { N_("/Data/Add variables/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables/random normal..."), NULL, 
      random_dialog, GENR_NORMAL, NULL },
    { N_("/Data/Add variables/random uniform..."), NULL, 
      random_dialog, GENR_UNIFORM, NULL },
    { N_("/Data/Add variables/seed generator..."), NULL, gretl_callback, 
      SEED, NULL },
    { N_("/Data/Refresh window"), NULL, refresh_data, 0, NULL },

    /* Sample menu */
    { N_("/_Sample"), NULL, NULL, 0, "<Branch>" },
    { N_("/Sample/_Set range..."), NULL, gretl_callback, SMPL, NULL },
    { N_("/Sample/_Restore full range"), NULL, restore_sample, 1, NULL },
    { N_("/Sample/sep1"), NULL, NULL, 0, "<Separator>" },    
    { N_("/Sample/Set _frequency, startobs..."), NULL, gretl_callback, 
      SETOBS, NULL },
    { N_("/Sample/sep2"), NULL, NULL, 0, "<Separator>" },   
    { N_("/Sample/_Define, based on dummy..."), NULL, gretl_callback, 
      SMPLDUM, NULL },
    { N_("/Sample/_Restrict, based on criterion..."), NULL, gretl_callback, 
      SMPLBOOL, NULL },
    { N_("/Sample/sep3"), NULL, NULL, 0, "<Separator>" },  
    { N_("/Sample/Drop all obs with _missing values"), NULL, bool_subsample, 
      0, NULL },
    { N_("/Sample/_Count missing values"), NULL, count_missing, 0, NULL },
    { N_("/Sample/Set missing _value code..."), NULL, gretl_callback, 
      GSETMISS, NULL },
    { N_("/Sample/sep4"), NULL, NULL, 0, "<Separator>" },  
    { N_("/Sample/_Add case markers..."), NULL, gretl_callback, MARKERS, NULL },
    { N_("/Sample/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/Sample/_Interpret as time series..."), NULL, time_series_dialog, 0, NULL },
    { N_("/Sample/Interpret as _panel..."), NULL, gui_set_panel_structure, 0, NULL },

    /* Variable menu */
    { N_("/_Variable"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Find..."), NULL, find_var, 0, NULL },
    { N_("/Variable/_Display values"), NULL, display_var, 0, NULL },
    { N_("/Variable/_Summary statistics"), NULL, do_menu_op, 
      VAR_SUMMARY, NULL },
    { N_("/Variable/_Time series plot"), NULL, do_graph_var, 0, NULL },
    { N_("/Variable/_Frequency distribution"), NULL, do_menu_op, 
      FREQ, NULL },
    { N_("/Variable/Frequency plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Frequency plot/simple"), NULL, do_freqplot, 0, NULL },
    { N_("/Variable/Frequency plot/against Normal"), NULL, do_freqplot, 
      NORMAL, NULL },
    { N_("/Variable/Frequency plot/against Gamma"), NULL, do_freqplot, 
      GAMMA, NULL },
    { N_("/Variable/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Correlogram"), NULL, gretl_callback, CORRGM, NULL },
    { N_("/Variable/Spectrum"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Spectrum/sample periodogram"), NULL, do_pergm, 0, NULL }, 
    { N_("/Variable/Spectrum/Bartlett lag window"), NULL, do_pergm, 1, NULL }, 
    { N_("/Variable/_Augmented Dickey-Fuller test"), NULL, gretl_callback, 
      ADF, NULL },
    { N_("/Variable/TRAMO analysis"), NULL, do_tramo, 0, NULL }, 
    { N_("/Variable/Range-mean graph"), NULL, do_range_mean, 0, NULL }, 
    { N_("/Variable/Runs test"), NULL, do_menu_op, RUNS, NULL }, 
    { N_("/Variable/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/_Rename"), NULL, mdata_edit, RENAME, NULL },
    { N_("/Variable/_Edit label"), NULL, mdata_edit, RELABEL, NULL },
#ifdef notdef
    { N_("/Variable/_Rename"), NULL, gretl_callback, RENAME, NULL },
    { N_("/Variable/_Edit label"), NULL, gretl_callback, RELABEL, NULL },
#endif
    { N_("/Variable/Set missing value code..."), NULL, gretl_callback, 
      VSETMISS, NULL },
    { N_("/Variable/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Simulate..."), NULL, gretl_callback, SIM, NULL },
    { N_("/Variable/Define _new variable..."), NULL, gretl_callback, GENR, NULL },
    { N_("/Variable/Delete last variable"), NULL, delete_var, 0, NULL },

    /* Model menu */
    { N_("/_Model"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model/_Ordinary Least Squares..."), NULL, model_callback, OLS, NULL },
    { N_("/Model/_Weighted Least Squares..."), NULL, model_callback, WLS, NULL },
    { N_("/Model/sep1"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/HCC_M..."), NULL, model_callback, HCCM, NULL },
    { N_("/Model/H_eteroskedasticity corrected..."), NULL, model_callback, 
      HSK, NULL },
    { N_("/Model/sep2"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Cochrane-Orcutt..."), NULL, model_callback, CORC, NULL },
    { N_("/Model/_Hildreth-Lu..."), NULL, model_callback, HILU, NULL },
    { N_("/Model/_Autoregressive estimation..."), NULL, model_callback, AR, NULL },
    { N_("/Model/sep3"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Vector Autoregression..."), NULL, model_callback, VAR, NULL },
    { N_("/Model/Cointe_gration test..."), NULL, selector_callback, COINT, NULL },
    { N_("/Model/_Two-Stage Least Squares..."), NULL, model_callback, TSLS, NULL },
    { N_("/Model/sep4"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Logit..."), NULL, model_callback, LOGIT, NULL },
    { N_("/Model/_Probit..."), NULL, model_callback, PROBIT, NULL },
    { N_("/Model/_Rank correlation..."), NULL, gretl_callback, SPEARMAN, NULL },
    { N_("/Model/_Pooled OLS (panel)..."), NULL, model_callback, POOLED, NULL },
#ifdef ENABLE_GMP
    { N_("/Model/High precision OLS..."), NULL, mp_ols_callback, MPOLS, NULL },
#endif

    /* Help menu */
    { N_("/_Help"), NULL, NULL, 0, "<LastBranch>" },
    { N_("/Help/_GUI commands"), NULL, do_gui_help, 0, NULL },
    { N_("/Help/_Script commands syntax"), NULL, do_script_help, 1, NULL },
    { N_("/Help/sep1"), NULL, NULL, 0, "<Separator>" },
#if defined(USE_GNOME)
    { N_("/Help/Manual in HTML"), NULL, gnome_help, 0, NULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>" },
#elif defined(G_OS_WIN32)
    { N_("/Help/Manual in HTML"), NULL, win_help, 0, NULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>" },
#endif
    { N_("/Help/_About gretl"), NULL, about_dialog, 0, NULL }
};

#ifndef G_OS_WIN32
static void make_userdir (PATHS *ppaths) 
{
    DIR *test;
    char buf[MAXLEN];
    
    if ((test = opendir(ppaths->userdir)) == NULL) {
	sprintf(buf, "mkdir -p %s", ppaths->userdir);
	system(buf);
	fprintf(stderr, _("Created user directory %s\n"
			  "If you prefer to use a different directory for "
			  "gretl user files, please make changes under\n"
			  "File, Preferences, General...\n"), ppaths->userdir);
    } else 
	closedir(test);
}
#endif

static void gui_usage (void)
{
    gui_logo(stdout);
    printf(_("You may supply the name of a data file on the command line.\n"));
    printf(_("Or you may do \"gretl -r script_file\" to open a script.\n"));
    printf(_("Or you may do \"gretl -d database\" to open a gretl database.\n"));
    exit(0);
}

static void noalloc (const char *str)
{
    fprintf(stderr, _("Couldn't allocate memory for %s\n"), str);
    exit(EXIT_FAILURE);
}

static void get_runfile (char *fname)
{
    int i;

    *tryscript = 0;
#ifdef G_OS_WIN32
    if (unmangle(fname, tryscript)) return;
#else
    strncat(tryscript, fname, MAXLEN-1);
#endif
    if (addpath(tryscript, &paths, 1) == NULL) {
	fprintf(stderr, _("Couldn't find script '%s'\n"), tryscript);
	exit(EXIT_FAILURE);
    } else {
	fprintf(stderr, _("%s found\n"), tryscript);
	i = slashpos(tryscript);
	if (i) {
	    paths.currdir[0] = 0;
	    strncat(paths.currdir, tryscript, i);
	}
	strcat(paths.currdir, SLASHSTR);
    }
}

static void fix_dbname (char *db)
{
    FILE *fp;

    if (strstr(db, ".bin") == NULL) strcat(db, ".bin");

    fp = fopen(db, "r");
    if (fp == NULL && strstr(db, paths.binbase) == NULL) {
	char tmp[MAXLEN];

	strcpy(tmp, db);
	build_path(paths.binbase, tmp, db, NULL);
    }
    if (fp != NULL) fclose(fp);
}

static void destroy (GtkWidget *widget, gpointer data)
{
    gtk_main_quit();
}

#ifdef G_OS_WIN32
extern int ws_startup (void);

void dummy_output_handler (const gchar *log_domain,
                           GLogLevelFlags log_level,
                           const gchar *message,
                           gpointer user_data)
{
    return;
}
#endif /* G_OS_WIN32 */

#ifdef ENABLE_NLS
void nls_init (void)
{
# ifdef G_OS_WIN32
    char gretldir[MAXSTR], LOCALEDIR[MAXSTR];

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretldir", gretldir))
	return;
    build_path(gretldir, "locale", LOCALEDIR, NULL);
#endif /* G_OS_WIN32 */

    setlocale (LC_ALL, "");
    bindtextdomain (PACKAGE, LOCALEDIR);
    textdomain (PACKAGE);
    bind_textdomain_codeset (PACKAGE, "UTF-8");
    nls_on = doing_nls();
}
#endif /* ENABLE_NLS */


int main (int argc, char *argv[])
{
    int err = 0, gui_get_data = 0;
    char dbname[MAXLEN];
#ifdef USE_GNOME
    GnomeProgram *program;
#endif

#ifdef ENABLE_NLS
    nls_init();
#endif       

    if ((errtext = malloc(MAXLEN)) == NULL) 
	noalloc(_("startup"));

    tryscript[0] = '\0';
    scriptfile[0] = '\0';
    paths.datfile[0] = '\0';

    /* Initialize gnome or GTK */
#ifdef USE_GNOME
    program = gnome_program_init ("gretl", version_string,
				  LIBGNOMEUI_MODULE, argc, argv,
				  GNOME_PARAM_POPT_TABLE, options,
				  GNOME_PARAM_HUMAN_READABLE_NAME,
				  _("The GNOME 2.0 econometrics package"),
				  GNOME_PARAM_APP_DATADIR, DATADIR,
				  GNOME_PARAM_NONE);
#else
    gtk_init(&argc, &argv);
#endif

    set_paths(&paths, 1, 1); /* 1 = defaults, 1 = gui */

#ifdef G_OS_WIN32
    read_rc(); /* get config info from registry */
    g_log_set_handler ("Gtk",
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    g_log_set_handler ("Gdk",
		       G_LOG_LEVEL_WARNING,
		       (GLogFunc) dummy_output_handler,
		       NULL);
    ws_startup(); 
    atexit(write_rc);
#else 
    set_rcfile();
    make_userdir(&paths);
#endif/* G_OS_WIN32 */

    if (argc > 1) {
	int opt = parseopt(argv[1]);

	switch (opt) {
	case OPT_HELP:
	    gui_usage();
	    break;
	case OPT_VERSION:
	    gui_logo(stdout);
	    exit(EXIT_SUCCESS);
	    break;
	case OPT_RUNIT:
#ifdef USE_GNOME
	    get_runfile(optrun);
#else
	    if (argc != 3) gui_usage();
	    get_runfile(argv[2]);
#endif
	    gui_get_data = 1;
	    break;
	case OPT_DBOPEN:
	case OPT_WEBDB:
#ifdef USE_GNOME
	    strncpy(dbname, optdb, MAXLEN-1);
#else
	    if (argc != 3) gui_usage();
	    strncpy(dbname, argv[2], MAXLEN-1);
#endif
	    if (opt == OPT_DBOPEN) fix_dbname(dbname);
	    gui_get_data = opt;
	    break;
	default:
	    /* let's suppose the argument is a data file */
	    break;
	}
    } else 
	gui_get_data = 1;

    strcpy(cmdfile, paths.userdir);
    strcat(cmdfile, "session.inp");
    cmds = gretl_print_new(GRETL_PRINT_FILE, cmdfile);
    if (cmds == NULL) {
	fprintf(stderr, _("Can't open file to save commands\n"));
	return EXIT_FAILURE;
    }
    fclose(cmds->fp);

    /* allocate data information struct */
    datainfo = datainfo_new();
    if (datainfo == NULL)
	noalloc(_("data information"));

    /* allocate memory for models */
    models = malloc(3 * sizeof *models);
    if (models == NULL) noalloc(_("models")); 
    models[0] = gretl_model_new(datainfo);
    models[1] = gretl_model_new(datainfo);
    models[2] = gretl_model_new(datainfo);
    if (models[0] == NULL || models[1] == NULL || models[2] == NULL) 
	noalloc(_("models")); 

    command.list = malloc(sizeof(int));
    command.param = malloc(1);
    if (command.list == NULL || command.param == NULL)  
	noalloc(_("command list")); 

    /* initialize random number generator */
    srand((unsigned) time(NULL));

    helpfile_init();
    session_init();

    /* get the data file, if specified on the command line */
    if (!(gui_get_data)) {
	int ftype;
	PRN *prn; 

	prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
	if (prn == NULL) 
	    exit(EXIT_FAILURE);
	paths.datfile[0] = '\0';
#ifdef G_OS_WIN32
	if (unmangle(argv[1], paths.datfile)) 
	    exit(EXIT_FAILURE);
#else
	strncat(paths.datfile, argv[1], MAXLEN-1);
#endif
	ftype = detect_filetype(paths.datfile, &paths, prn);

	switch (ftype) {
	case GRETL_UNRECOGNIZED:
	    exit(EXIT_FAILURE);
	case GRETL_NATIVE_DATA:
	    err = get_data(&Z, datainfo, paths.datfile, &paths, data_status, 
			   prn);
	    break;
	case GRETL_XML_DATA:
	    err = get_xmldata(&Z, datainfo, paths.datfile, &paths, data_status, 
			      prn, 0);
	    break;
	case GRETL_CSV_DATA:
	    err = import_csv(&Z, datainfo, paths.datfile, prn);
	    break;
	case GRETL_BOX_DATA:
	    err = import_box(&Z, datainfo, paths.datfile, prn);
	    break;
	case GRETL_SCRIPT:
	    gui_get_data = 1;
	    get_runfile(paths.datfile);
	    paths.datfile[0] = '\0';
	    break;
	}
	if (ftype != GRETL_SCRIPT && err) {
	    errmsg(err, prn);
	    exit(EXIT_FAILURE);
	}
	gretl_print_destroy(prn);
    }

    /* create the GUI */
    gretl_tips = gtk_tooltips_new();

    /* create main window */
    if ((mdata = mymalloc(sizeof(windata_t))) == NULL)
	noalloc(_("GUI"));
    if ((datalabel = make_main_window(gui_get_data)) == NULL) 
	noalloc(_("main window"));
    if (!gui_get_data) set_sample_label(datainfo);

    /* enable special copying to clipboard */
    clip_init(mdata->w);

    init_fileptrs();
    add_files_to_menu(1);
    add_files_to_menu(2);
    add_files_to_menu(3);
#ifndef GNUPLOT_PNG
    graphmenu_state(FALSE);
#endif
    session_menu_state(FALSE);
    restore_sample_state(FALSE);
    main_menubar_state(FALSE);
			  
    check_for_extra_data();

    if (!gui_get_data)
	register_data(paths.datfile, 1);

    /* opening a script from the command line? */
    if (tryscript[0] != '\0') { 
	do_open_script();
    }

    /* check for program updates? */
    proxy_init(dbproxy);
    if (updater) update_query(); 

    /* try opening specified database */
    if (gui_get_data == OPT_DBOPEN)
	open_named_db_list(dbname);
    else if (gui_get_data == OPT_WEBDB)
	open_named_remote_db_list(dbname);

    /* Enter the event loop */
    gtk_main();

    /* clean up before exiting */
    /* if (mdata) free_windata(NULL, mdata); */
    free_session();
    if (Z) free_Z(Z, datainfo);
    if (fullZ) free_Z(fullZ, fullinfo);
    free_model(models[0]);
    free_model(models[1]);
    free_model(models[2]);
    free(models);
    if (command.list != NULL) free(command.list);
    if (command.param != NULL) free(command.param);
    if (data_status) free_datainfo(datainfo);
    if (fullinfo) {
	clear_datainfo(fullinfo, CLEAR_SUBSAMPLE);
	free(fullinfo);
    }

    free_command_stack();
    free_modelspec();

    remove(paths.plotfile);
    Rcleanup();

    return EXIT_SUCCESS;
}

/* ........................................................... */

void refresh_data (void)
{
    if (data_status)
	populate_list(mdata->listbox, datainfo);
}

/* ........................................................... */

void main_menubar_state (gboolean s)
{
    if (mdata->ifac == NULL) return;

    flip(mdata->ifac, "/File/Append data", s);
    flip(mdata->ifac, "/File/Clear data set", s);
    flip(mdata->ifac, "/File/Save data", s);
    flip(mdata->ifac, "/File/Save data as", s);
    flip(mdata->ifac, "/File/Export data", s);
    flip(mdata->ifac, "/File/Create data set", !s);
    flip(mdata->ifac, "/Data", s);
    flip(mdata->ifac, "/Sample", s);
    flip(mdata->ifac, "/Variable", s);
    flip(mdata->ifac, "/Model", s);

    if (s && (data_status & BOOK_DATA))
	flip(mdata->ifac, "/Data/Edit info", 0);
}

/* ........................................................... */

void graphmenu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/File/Save last graph", s);
	flip(mdata->ifac, "/Session/Add last graph", s);
    }
}

/* ........................................................... */

void panel_menu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Model/Pooled OLS (panel)...", s);
	flip(mdata->ifac, "/Data/Add variables/panel dummies", s);
    }
}

/* ........................................................... */

void session_menu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Session/Icon view", s);
	flip(mdata->ifac, "/Session/Save", s);
	flip(mdata->ifac, "/Session/Save as...", s);
    }	
}

/* ........................................................... */

void restore_sample_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Sample/Restore full range", s);
	flip(mdata->ifac, "/Variable/Delete last variable", !s);
    }
}

/* ........................................................... */

static void populate_list (GtkWidget *widget, DATAINFO *datainfo)
{
    GtkListStore *store;
    GtkTreeSelection *select;
    GtkTreeIter iter;    
    char id[4];
    gint i;

    /* find and clear the existing list */
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(widget)));
    gtk_list_store_clear (store);

    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);

    for (i=0; i<datainfo->v; i++) {
	if (hidden_var(i, datainfo)) continue;
	gtk_list_store_append(store, &iter);
	sprintf(id, "%d", i);
	gtk_list_store_set (store, &iter, 
			    0, id, 
			    1, datainfo->varname[i],
			    2, datainfo->label[i],
			    3, FALSE, /* (i > 0)? TRUE : FALSE */
			    -1);
    }    

    mdata->active_var = 1;
    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);
    gtk_tree_model_iter_next (GTK_TREE_MODEL(store), &iter);
    select = gtk_tree_view_get_selection (GTK_TREE_VIEW(widget));
    gtk_tree_selection_select_iter (select, &iter);

    if (mdata->popup == NULL) {
	build_var_popup(mdata);
	build_selection_popup();
	g_signal_connect(G_OBJECT(mdata->listbox), "button_press_event",
			 G_CALLBACK(main_popup_handler), 
			 mdata->popup);
	g_signal_connect (G_OBJECT(mdata->listbox), "button_press_event",
			  G_CALLBACK(main_varclick),
			  mdata);
    }
}

/* ......................................................... */

static void sort_varlist (gpointer p, guint col, GtkWidget *w)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_sortable_set_sort_column_id (GTK_TREE_SORTABLE(model), 
					  col, GTK_SORT_ASCENDING);
}

/* ......................................................... */

void populate_varlist (void)
{
    populate_list (mdata->listbox, datainfo);
}

/* ......................................................... */

void clear_varlist (GtkWidget *widget)
{
    GtkListStore *store;

    store = GTK_LIST_STORE(gtk_tree_view_get_model (GTK_TREE_VIEW(widget)));
    gtk_list_store_clear (store);
}

/* ......................................................... */

void clear_sample_label (void)
{
    gtk_label_set_text(GTK_LABEL(mdata->status), "");
    gtk_label_set_text(GTK_LABEL(datalabel), _(" No datafile loaded "));
}

/* ......................................................... */

void set_sample_label (DATAINFO *pdinfo)
{
    char startdate[9], enddate[9], pdstr[10];
    char labeltxt[80];

    ntodate(startdate, pdinfo->t1, pdinfo);
    ntodate(enddate, pdinfo->t2, pdinfo);

    if (dataset_is_time_series(pdinfo)) {
	switch (pdinfo->pd) {
	case 1:
	    strcpy(pdstr, _("Annual")); break;
	case 4:
	    strcpy(pdstr, _("Quarterly")); break;
	case 12:
	    strcpy(pdstr, _("Monthly")); break;
	case 24:
	    strcpy(pdstr, _("Hourly")); break;
	case 52:
	    strcpy(pdstr, _("Weekly")); break;
	case 5:
	    strcpy(pdstr, _("Daily")); break;
	case 7:
	    strcpy(pdstr, _("Daily")); break;
	default:
	    strcpy(pdstr, _("Unknown")); break;
	}
    } 
    else if (dataset_is_panel(pdinfo)) 
	strcpy(pdstr, _("Panel"));
    else 
	strcpy(pdstr, _("Undated"));

    panel_menu_state(dataset_is_panel(pdinfo));

    flip(mdata->ifac, "/Sample/Interpret as time series...", 
	 !(dataset_is_time_series(pdinfo)));

    flip(mdata->ifac, "/Sample/Interpret as panel...", 
	 !(pdinfo->pd == 1));

    sprintf(labeltxt, _("%s: Full range %s - %s; current sample"
			" %s - %s"), pdstr, pdinfo->stobs, pdinfo->endobs,
	    startdate, enddate);
    gtk_label_set_text(GTK_LABEL(mdata->status), labeltxt);

    if (strlen(paths.datfile) > 2) {
	if (strrchr(paths.datfile, SLASH) == NULL)
	    sprintf(labeltxt, " %s ", paths.datfile);
	else
	    sprintf(labeltxt, " %s ", 
		    strrchr(paths.datfile, SLASH) + 1);
	if (data_status & MODIFIED_DATA) 
	    strcat(labeltxt, "* ");
	if (datalabel != NULL)
	    gtk_label_set_text(GTK_LABEL(datalabel), labeltxt);
    } 
    else if (data_status & MODIFIED_DATA) {
	strcpy(labeltxt, _(" Unsaved data "));
	gtk_label_set_text(GTK_LABEL(datalabel), labeltxt);
    }
}

/* ......................................................... */

#ifdef USE_WINFONT

#define NAME_BUFFER_LEN 32

static void try_to_get_windows_font (void)
{
    HDC h_dc;
    HGDIOBJ h_font;
    TEXTMETRIC tm;
    char name[NAME_BUFFER_LEN], regfont[MAXLEN];

    /* don't override user's choice of font */
    if (read_reg_val(HKEY_CURRENT_USER, "App_font", regfont) == 0)
	return;

    h_dc = CreateDC("DISPLAY", NULL, NULL, NULL);
    if (h_dc == NULL) return;

    h_font = GetStockObject(DEFAULT_GUI_FONT); 
    if (h_font == NULL || !SelectObject(h_dc, h_font)) {
	DeleteDC(h_dc);
	return;
    }

    if (GetTextFace(h_dc, NAME_BUFFER_LEN, name) <= 0) {
	DeleteDC(h_dc);
	return;
    }

    if (!GetTextMetrics(h_dc, &tm)) {
	DeleteDC(h_dc);
	return;
    } else {
	HDC screen = GetDC(0);
	double scaleY = GetDeviceCaps(screen, LOGPIXELSY) / 96.0;
	int pix_height = (int) (tm.tmHeight * scaleY);
	int match = 0;
	PangoFontDescription *pfd;
	PangoFont *pfont;
	PangoContext *pc;
	GtkWidget *w;
	gchar *fontname;

	ReleaseDC(0, screen);
	DeleteDC(h_dc);

	fontname = g_strdup_printf("%s %d", name, pix_height);
	pfd = pango_font_description_from_string(fontname);

	w = gtk_label_new(NULL);
	pc = gtk_widget_get_pango_context(w);
	pfont = pango_context_load_font(pc, pfd);
	match = (pfont != NULL);

	pango_font_description_free(pfd);
	g_object_unref(G_OBJECT(pc));
	gtk_widget_destroy(w);

	if (match) set_app_font(fontname);
	g_free(fontname);
    }
}
#endif /* USE_WINFONT */

/* ......................................................... */

static GtkWidget *make_main_window (int gui_get_data) 
{
    GtkWidget *box, *dlabel, *align;
    const char *titles[] = {
	_("ID #"), 
	_("Variable name"), 
	_("Descriptive label")
    };
    int listbox_data_width = 500;
    int listbox_file_height = 300;

    mdata->data = NULL;  
    mdata->listbox = NULL;
    mdata->popup = NULL;
    mdata->role = MAINWIN;

#ifdef USE_GNOME
    mdata->w = gnome_app_new("gretl", _("Econometrics program"));
#else
    mdata->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
#endif

#ifdef USE_WINFONT
    try_to_get_windows_font();
#endif

    g_signal_connect (G_OBJECT (mdata->w), "delete_event",
		      G_CALLBACK (exit_check), NULL);
    g_signal_connect (G_OBJECT (mdata->w), "destroy",
		      G_CALLBACK (destroy), NULL);

    gtk_window_set_title(GTK_WINDOW (mdata->w), "gretl");
    gtk_window_set_resizable(GTK_WINDOW (mdata->w), TRUE);
#ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(mdata->w), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif

    main_vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER (main_vbox), 10);

#ifdef USE_GNOME
    gnome_app_set_contents (GNOME_APP (mdata->w), main_vbox);
#else
    gtk_container_add (GTK_CONTAINER (mdata->w), main_vbox);
#endif

    set_up_main_menu();
    gtk_box_pack_start(GTK_BOX (main_vbox), mdata->mbar, FALSE, TRUE, 0);
    gtk_widget_show(mdata->mbar);

    dlabel = gtk_label_new(_(" No datafile loaded ")); 
    gtk_widget_show(dlabel);

    box = gtk_vbox_new (FALSE, 0);
    gtk_widget_set_size_request(box, listbox_data_width, listbox_file_height);
    /* gtk_container_set_border_width (GTK_CONTAINER (box), 5); */
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(box), align, FALSE, FALSE, 0);
    gtk_widget_show(align);
    gtk_container_add(GTK_CONTAINER(align), dlabel);
   
    mdata->listbox = list_box_create (mdata, GTK_BOX(box), 3, titles);
    gtk_widget_show(box);

    gtk_drag_dest_set (mdata->listbox,
		       GTK_DEST_DEFAULT_ALL,
		       target_table, 1,
		       GDK_ACTION_COPY);

    g_signal_connect (G_OBJECT(mdata->listbox), "drag_data_received",
		      G_CALLBACK(drag_data_received),
		      NULL);

    gtk_box_pack_start(GTK_BOX (main_vbox), box, TRUE, TRUE, 0);

    mdata->status = gtk_label_new("");
    
    gtk_box_pack_start (GTK_BOX (main_vbox), mdata->status, FALSE, TRUE, 0);

    /* put stuff into list box, activate menus */
    if (!gui_get_data) populate_list(mdata->listbox, datainfo);

    /* create gretl toolbar */
    if (want_toolbar) make_toolbar(mdata->w, main_vbox);

    /* get a monospaced font for various windows */
    set_fixed_font();

    /* and a proportional font for menus, etc */
#ifndef USE_GNOME
    set_app_font(NULL);
#endif

    gtk_widget_show_all(mdata->w); 

    return dlabel;
}

/* ........................................................... */

static void set_up_main_menu (void)
{
    gint n_items = sizeof data_items / sizeof data_items[0];

    mdata->ifac = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>", 
					NULL);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(mdata->ifac, menu_translate, NULL, NULL);
#endif    
    gtk_item_factory_create_items (mdata->ifac, n_items, data_items, NULL);
    mdata->mbar = gtk_item_factory_get_widget (mdata->ifac, "<main>");
}

/* ........................................................... */

static void mdata_edit (gpointer data, guint action, GtkWidget *w)
{
    gint row, colnum;
    GtkTreeViewColumn *col;

    if (action == RENAME) colnum = 1;
    else if (action == RELABEL) colnum = 2;
    else return;

    row = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(mdata->listbox),
					    "active_row"));
    col = gtk_tree_view_get_column(GTK_TREE_VIEW(mdata->listbox), colnum);

    if (col != NULL) {
	GtkTreeIter iter;
	GtkTreeModel *model;
	GtkListStore *store;
	GtkTreePath *path;
	gchar rowstr[8];

	model = gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));
	store = GTK_LIST_STORE(model);

	sprintf(rowstr, "%d", row);
	path = gtk_tree_path_new_from_string(rowstr);
	gtk_tree_model_get_iter(model, &iter, path);
	gtk_list_store_set(store, &iter, 3, TRUE, -1);

	gtk_tree_view_row_activated(GTK_TREE_VIEW(mdata->listbox), path, col);
	gtk_tree_view_set_cursor(GTK_TREE_VIEW(mdata->listbox),
				 path, col, TRUE);
	gtk_tree_path_free(path);
    }
}

extern void delete_var_by_id (int id); /* callbacks.c */

static gint var_popup_click (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Display values"))) display_var();
    if (!strcmp(item, _("Descriptive statistics"))) 
	do_menu_op(NULL, VAR_SUMMARY, NULL);
    else if (!strcmp(item, _("Time series plot"))) do_graph_var();
    else if (!strcmp(item, _("Frequency distribution"))) 
	do_menu_op(NULL, FREQ, NULL);
    else if (!strcmp(item, _("Frequency plot"))) do_freqplot(NULL, 0, NULL);
    else if (!strcmp(item, _("Boxplot")))
	do_boxplot_var();
    else if (!strcmp(item, _("Correlogram"))) 
	gretl_callback(NULL, CORRGM, NULL);
    else if (!strcmp(item, _("Spectrum"))) 
	do_pergm(NULL, 0, NULL);
    else if (!strcmp(item, _("Dickey-Fuller test"))) 
	gretl_callback(NULL, ADF, NULL);
    else if (!strcmp(item, _("Runs test"))) 
	do_menu_op(NULL, RUNS, NULL);
    else if (!strcmp(item, _("Rename"))) 
	mdata_edit(NULL, RENAME, NULL);
    else if (!strcmp(item, _("Edit label")))  
	mdata_edit(NULL, RELABEL, NULL);
    else if (!strcmp(item, _("Delete"))) 
	delete_var_by_id(mdata->active_var);
    else if (!strcmp(item, _("Simulate..."))) 
	gretl_callback(NULL, SIM, NULL);
    else if (!strcmp(item, _("Define new variable..."))) 
	gretl_callback(NULL, GENR, NULL);

    return TRUE;
}

static void build_var_popup (windata_t *win)
{
    const char *var_items[]={
	N_("Display values"),
	N_("Descriptive statistics"),
	N_("Time series plot"),
	N_("Frequency distribution"),
	N_("Frequency plot"),
	N_("Boxplot"),
	N_("Correlogram"),
	N_("Spectrum"),
	N_("Dickey-Fuller test"),
	N_("Runs test"),
	N_("Rename"),
	N_("Edit label"),
	N_("Delete"),
	N_("Simulate..."),
	N_("Define new variable...")
    };

    GtkWidget *var_item;
    int i;

    win->popup = gtk_menu_new();

    for (i=0; i<(sizeof var_items / sizeof var_items[0]); i++) {
	var_item = gtk_menu_item_new_with_label(_(var_items[i]));
	g_signal_connect(G_OBJECT(var_item), "activate",
			 G_CALLBACK(var_popup_click),
			 _(var_items[i]));
	gtk_widget_show(var_item);
	gtk_menu_shell_append(GTK_MENU_SHELL(win->popup), var_item);
    }
}

/* ........................................................... */

static gint selection_popup_click (GtkWidget *widget, gpointer data)
{
    gchar *item = (gchar *) data;

    if (!strcmp(item, _("Display values"))) 
	display_selected(NULL, 0, NULL); 
    else if (!strcmp(item, _("Descriptive statistics"))) 
	do_menu_op(NULL, SUMMARY_SELECTED, NULL);
    else if (!strcmp(item, _("Correlation matrix"))) 
	do_menu_op(NULL, CORR_SELECTED, NULL);
    else if (!strcmp(item, _("Time series plot"))) 
	plot_from_selection(NULL, GR_PLOT, NULL);
    return TRUE;
}

static void build_selection_popup (void)
{
    const char *items[] = {
	N_("Display values"),
	N_("Descriptive statistics"),
	N_("Correlation matrix"),
	N_("Time series plot")
    };

    GtkWidget *item;
    int i;

    selection_popup = gtk_menu_new();

    for (i=0; i<(sizeof items / sizeof items[0]); i++) {
	item = gtk_menu_item_new_with_label(_(items[i]));
	g_signal_connect(G_OBJECT(item), "activate",
			 G_CALLBACK(selection_popup_click),
			 _(items[i]));
	gtk_widget_show(item);
	gtk_menu_shell_append(GTK_MENU_SHELL(selection_popup), item);
    }
}

/* ........................................................... */

static void check_for_extra_data (void)
{
    DIR *dir;
    extern char pwtpath[MAXLEN]; /* datafiles.c */
    extern char woolpath[MAXLEN]; /* fileselect.c */
    const char *pwt_menu_paths[] = {
	"/File/Open data/sample file/Penn World Table...",
	"/File/Open command file/practice file/Penn World Table..."
    };
    const char *wool_menu_path =
	"/File/Open data/sample file/Wooldridge...";
    int gotpwt = 0, gotwool = 0;

    /* first check for Penn World Table */
    build_path(paths.datadir, "pwt56", pwtpath, NULL); 
    /* try at system level */
    if ((dir = opendir(pwtpath)) != NULL) {
        closedir(dir);
        gotpwt = 1;
    } else {
        build_path(paths.userdir, "pwt56", pwtpath, NULL); 
	/* and at user level */
        if ((dir = opendir(pwtpath)) != NULL) {
            closedir(dir);
            gotpwt = 1; 
        }
    }

    if (!gotpwt) {
        flip (mdata->ifac, pwt_menu_paths[0], FALSE);
        flip (mdata->ifac, pwt_menu_paths[1], FALSE);
    }

    /* then check for Wooldridge data */
    build_path(paths.datadir, "wooldridge", woolpath, NULL); 
    /* try at system level */
    if ((dir = opendir(woolpath)) != NULL) {
        closedir(dir);
        gotwool = 1;
    } else {
        build_path(paths.userdir, "wooldridge", woolpath, NULL); 
	/* and at user level */
        if ((dir = opendir(woolpath)) != NULL) {
            closedir(dir);
            gotwool = 1;
        }
    }

    if (!gotwool) flip (mdata->ifac, wool_menu_path, FALSE);
}

/* ........................................................... */

void restore_sample (gpointer data, int verbose, GtkWidget *w)
{
    int err = 0;

    err = restore_full_sample(&subZ, &fullZ, &Z,
			      &subinfo, &fullinfo, &datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }
    if (verbose) {
	infobox(_("Full sample range restored"));
	set_sample_label(datainfo);    
	restore_sample_state(FALSE);
	strcpy(line, "smpl full");
	check_cmd(line);
	cmd_init(line);
    }
}

/* ........................................................... */

#ifdef G_OS_WIN32

BOOL CreateChildProcess (char *prog) 
{ 
    PROCESS_INFORMATION piProcInfo; 
    STARTUPINFO siStartInfo; 
 
    ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));
    ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
    siStartInfo.cb = sizeof(STARTUPINFO); 
 
    return CreateProcess(
			 NULL, 
			 prog,          /* command line */
			 NULL,          /* process security attributes  */
			 NULL,          /* primary thread security attributes */ 
			 TRUE,          /* handles are inherited  */
			 0,             /* creation flags  */
			 NULL,          /* use parent's environment  */
			 NULL,          /* use parent's current directory  */
			 &siStartInfo,  /* STARTUPINFO pointer */ 
			 &piProcInfo);  /* receives PROCESS_INFORMATION  */
}

#else

void gretl_fork (const char *prog, const char *arg)
{
    pid_t pid;

    signal(SIGCLD, SIG_IGN);

    pid = fork();
    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return;
    } else if (pid == 0) {
	if (arg != NULL) 
	    execlp(prog, prog, arg, NULL);
	else
	    execlp(prog, prog, NULL);
	perror("execlp");
	_exit(EXIT_FAILURE);
    }
}

#endif	

/* ........................................................... */

static void startR (gpointer p, guint opt, GtkWidget *w)
{
    char Rdata[MAXLEN], line[MAXLEN];
    FILE *fp;
#ifndef G_OS_WIN32
    int i;
    char *s0, *s1, *s2;
    pid_t pid;
#endif

    if (!data_status) {
	errbox(_("Please open a data file first"));
	return;
    }

    fp = fopen(".Rprofile", "r");
    if (fp != NULL) {
	fclose(fp);
	if (copyfile(".Rprofile", ".Rprofile.gretltmp")) {
	    errbox(_("Couldn't move existing .Rprofile out of the way"));
	    return;
	}
    }
    fp = fopen(".Rprofile", "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }
    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);
    sprintf(line, "store -r %s", Rdata); 
    if (check_cmd(line) || cmd_init(line) ||
	write_data(Rdata, command.list, Z, datainfo, GRETL_DATA_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }
    if (dataset_is_time_series(datainfo)) {
	fprintf(fp, "source(\"%s\")\n", Rdata);
	fprintf(fp, "ls()\n");
    } else {
	fprintf(fp, "gretldata <- read.table(\"%s\")\n", Rdata);
	fprintf(fp, "attach(gretldata)\n");
	fprintf(fp, "ls(2)\n");
    }
    fclose(fp);

#ifdef G_OS_WIN32
    CreateChildProcess(Rcommand);
#else
    s0 = mymalloc(64);
    s1 = mymalloc(32);
    s2 = mymalloc(32);
    if (s0 == NULL || s1 == NULL || s2 == NULL)
	return;
    *s0 = *s1 = *s2 = '\0';
    i = sscanf(Rcommand, "%63s %31s %31s", s0, s1, s2);
    if (i == 0) {
	errbox(_("No command was supplied to start R"));
	free(s0); free(s1); free(s2);
	return;
    }

    signal(SIGCLD, SIG_IGN); 
    pid = fork();

    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return;
    } else if (pid == 0) {  
	if (i == 1)
	    execlp(s0, s0, NULL);
	else if (i == 2)
	    execlp(s0, s0, s1, NULL);
	else if (i == 3)
	    execlp(s0, s0, s1, s2, NULL);
	perror("execlp");
	_exit(EXIT_FAILURE);
    }
    free(s0); 
    free(s1); 
    free(s2);
#endif 
}

/* ........................................................... */

static void Rcleanup (void)
{
    FILE *fp;
    char Rdata[MAXLEN];

    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);
    remove(Rdata);

    fp = fopen(".Rprofile.gretltmp", "r");
    if (fp != NULL) {
	fclose(fp);
	if (copyfile(".Rprofile.gretltmp", ".Rprofile")) 
	    errbox(_("Error restoring .Rprofile from\n"
		     "the temporary copy, .Rprofile.gretltmp"));
	else 
	    remove(".Rprofile.gretltmp");
    }
}

/* ........................................................... */

static void show_calc (void)
{
#ifdef G_OS_WIN32
    CreateChildProcess(calculator);
#else
    gretl_fork(calculator, NULL);
#endif 
}

/* ........................................................... */

static void show_edit (void)
{
#ifdef G_OS_WIN32
    CreateChildProcess(editor);
#else
    gretl_fork(editor, NULL);
#endif 
}

/* ........................................................... */

static void open_ramudata (void)
{
    display_files(NULL, RAMU_DATA, NULL);
}

/* ........................................................... */

void show_toolbar (void)
{
    make_toolbar(mdata->w, main_vbox);
}

/* ........................................................... */

#ifdef G_OS_WIN32
extern int goto_url (const char *url);
#else
static void netscape_open (const char *url)
{
#ifdef USE_GNOME
    gnome_url_show(url, NULL);   
#else
    int err;
    char ns_cmd[128];

    sprintf(ns_cmd, "netscape -remote \"openURLNewWindow(%s)\"", url);
    err = system(ns_cmd);
    if (err) gretl_fork("netscape", url);
#endif /* USE_GNOME */
}
#endif /* G_OS_WIN32 */

static void gretl_website (void)
{
#ifdef G_OS_WIN32
    if (goto_url("http://gretl.sourceforge.net/"))
	errbox("Failed to open URL");
#else
    netscape_open("http://gretl.sourceforge.net/");
#endif
}

static void gretl_pdf (void)
{
#ifdef G_OS_WIN32
    if (goto_url("http://gretl.sourceforge.net/manual.pdf"))
	errbox(_("Failed to open URL"));
#else
    netscape_open("http://gretl.sourceforge.net/manual.pdf");
#endif
}

static void xy_graph (void)
{
    if (data_status)
	selector_callback(NULL, GR_XY, NULL);
    else
	errbox(_("Please open a data file first"));
}

static void go_session (void)
{
    if (data_status)
	view_session();
    else
	errbox(_("Please open a data file first"));
}

/* ........................................................... */

static GtkWidget *image_button_new (GdkPixmap *pix, GdkBitmap *mask,
				    void (*toolfunc)())
{
    GtkWidget *image = gtk_image_new_from_pixmap(pix, mask);
    GtkWidget *button = gtk_button_new();

    gtk_widget_set_size_request(button, 24, 24);

    gtk_container_add (GTK_CONTAINER(button), image);
    g_signal_connect (G_OBJECT(button), "clicked",
                      G_CALLBACK(toolfunc), NULL);

    return button;
}

/* ........................................................... */

static void make_toolbar (GtkWidget *w, GtkWidget *box)
{
    GtkWidget *button, *hbox;
    GdkPixmap *icon;
    GdkBitmap *mask;
    GdkColormap *cmap;
    int i;
    const char *toolstrings[] = {
	N_("launch calculator"), 
	N_("launch editor"), 
	N_("open gretl console"),
	N_("session icon view"),
	N_("gretl website"), 
	N_("gretl manual (PDF)"),
	N_("show help"), 
	N_("X-Y graph"), 
	N_("Capture last graph for editing"),
	N_("open dataset"),
	NULL
    };
    gchar **toolxpm = NULL;
    void (*toolfunc)() = NULL;
    const char *toolstr;

    cmap = gdk_colormap_get_system();

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);

    toolbar_box = gtk_handle_box_new();
    gtk_handle_box_set_shadow_type(GTK_HANDLE_BOX(toolbar_box), NONE);
    gtk_box_pack_start(GTK_BOX(hbox), toolbar_box, FALSE, FALSE, 0);

    gretl_toolbar = gtk_toolbar_new();

    gtk_container_set_border_width(GTK_CONTAINER(gretl_toolbar), 0);
    gtk_container_add(GTK_CONTAINER(toolbar_box), gretl_toolbar);

    for (i=0; toolstrings[i] != NULL; i++) {
	switch (i) {
	case 0:
	    toolxpm = mini_calc_xpm;
	    toolfunc = show_calc;
	    break;
	case 1:
	    toolxpm = mini_edit_xpm;
	    toolfunc = show_edit;
	    break;
	case 2:
	    toolxpm = mini_sh_xpm;
	    toolfunc = console;
	    break;
	case 3:
	    toolxpm = mini_session_xpm;
	    toolfunc = go_session;
	    break;
	case 4:
	    toolxpm = mini_netscape_xpm;
	    toolfunc = gretl_website;
	    break;  
	case 5:
	    toolxpm = mini_pdf_xpm;
	    toolfunc = gretl_pdf;
	    break;    
	case 6:
	    toolxpm = mini_manual_xpm;
	    toolfunc = do_gui_help;
	    break;
	case 7:
	    toolxpm = mini_plot_xpm;
	    toolfunc = xy_graph;
	    break;
	case 8:
#ifndef GNUPLOT_PNG
	    toolxpm = mini_camera_xpm;
	    toolfunc = add_last_graph;
#endif
	    break;
	case 9:
	    toolxpm = mini_ofolder_xpm;
	    toolfunc = open_ramudata;
	    break;
	default:
	    break;
	}

#ifdef GNUPLOT_PNG
	if (i == 8) continue;
#endif

	toolstr = _(toolstrings[i]);
	icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, 
						     NULL, toolxpm);
	button = image_button_new(icon, mask, toolfunc);
	gtk_toolbar_append_widget(GTK_TOOLBAR(gretl_toolbar), button,
				  toolstr, NULL);
    }

    gtk_widget_show_all (hbox);
}

/* Icon handling for X */
#ifndef G_OS_WIN32
void set_wm_icon (GtkWidget *w, gpointer data)
{
    GdkPixmap *icon;

    icon = gdk_pixmap_create_from_xpm_d(w->window, NULL, NULL, gretl_xpm);
    gdk_window_set_icon(w->window, NULL, icon, NULL);
}
#endif

/* Drag 'n' drop */
static void  
drag_data_received  (GtkWidget *widget,
		     GdkDragContext *context,
		     gint x,
		     gint y,
		     GtkSelectionData *data,
		     guint info,
		     guint time,
		     gpointer p)
{
    gchar *dfname;
    char *suff = NULL, tmp[MAXLEN];
    int pos, skip = 5;

    /* ignore the wrong sort of data */
    if (data == NULL || (dfname = data->data) == NULL || 
	strlen(dfname) <= 5 || strncmp(dfname, "file:", 5))
	return;

    if (strncmp(dfname, "file://", 7)) skip = 7;

    /* there may be multiple files: we ignore all but the first */
    tmp[0] ='\0';
    if ((pos = haschar('\r', dfname)) > 0 || 
	(pos = haschar('\n', dfname) > 0)) {
	strncat(tmp, dfname + skip, pos - skip);
    } else
	strcat(tmp, dfname + skip);

#ifdef G_OS_WIN32
    if (unmangle(tmp, tryscript)) return;
    strcpy(tmp, tryscript);
    tryscript[0] = '\0';
#endif

    suff = strrchr(tmp, '.');
    if (suff && (!strncmp(suff, ".gretl", 6) || 
		 !strncmp(suff, ".inp", 4) ||
		 !strncmp(suff, ".GRE", 4) ||
		 !strncmp(suff, ".INP", 4))) {
	strcpy(tryscript, tmp);
	verify_open_session(NULL);
    } else {
	strcpy(trydatfile, tmp);
	verify_open_data(NULL, 0);
    }	
}

/* ........................................................... */

static gint 
special_selection_get (GtkWidget *widget,
		       GtkSelectionData *selection_data,
		       guint info,
		       guint time)
{
    gchar *str;
    gint length;

    str = clipboard_buf;
    if (str == NULL) return TRUE;
    length = strlen(str);
  
    if (info == TARGET_STRING) {
	gtk_selection_data_set (selection_data,
				GDK_SELECTION_TYPE_STRING,
				8 * sizeof(gchar), 
				(guchar *) str, 
				length);
    } else if (info == TARGET_TEXT || info == TARGET_COMPOUND_TEXT) {
	guchar *text;
	gchar c;
	GdkAtom encoding;
	gint format;
	gint new_length;

	c = str[length];
	str[length] = '\0';
	gdk_string_to_compound_text(str, &encoding, &format, 
				    &text, &new_length);
	gtk_selection_data_set(selection_data, encoding, format, 
			       text, new_length);
	gdk_free_compound_text(text);
	str[length] = c;
    }
    g_free(str);
    clipboard_buf = NULL;
    return TRUE;
}

/* ........................................................... */

static void clip_init (GtkWidget *w)
{
    GdkAtom clipboard_atom = GDK_NONE;
    GtkTargetEntry targets[] = {
	{ "STRING", 0, TARGET_STRING },
	{ "TEXT",   0, TARGET_TEXT }, 
	{ "COMPOUND_TEXT", 0, TARGET_COMPOUND_TEXT }
    };

    gint n_targets = sizeof(targets) / sizeof(targets[0]);
  
    clipboard_atom = gdk_atom_intern("CLIPBOARD", FALSE);
    gtk_selection_add_targets(w, GDK_SELECTION_PRIMARY,
			      targets, n_targets);
    gtk_selection_add_targets(w, clipboard_atom,
			      targets, n_targets);
    g_signal_connect (G_OBJECT(mdata->w), "selection_get",
		      G_CALLBACK(special_selection_get), NULL);    
}

/* ........................................................... */

static void auto_store (void)
{
    int opt = 0;

    if (make_default_storelist()) return;

    /* if there's already a datafile, and it's gzipped, then
       arrange for the new store to be gzipped too */
    if (strlen(paths.datfile) && is_gzipped(paths.datfile))
	opt = OPT_Z;

    if (data_status & USER_DATA)
	do_store(paths.datfile, opt, 1);
    else
	file_selector(_("Save data file"), SAVE_DATA, NULL);	
}

/* ........................................................... */

#ifdef G_OS_WIN32

static int old_windows (void) {
    OSVERSIONINFO *winver;
    static int old = 1;

    if (!old) return 0; /* do only one look up */

    winver = mymalloc(sizeof *winver);
    if (winver == NULL) return old;
    winver->dwOSVersionInfoSize = sizeof *winver;
    GetVersionEx(winver);
    switch (winver->dwPlatformId) {
    case VER_PLATFORM_WIN32_WINDOWS:
        if (winver->dwMinorVersion >= 10) /* win98 or higher */
	    old = 0;
        break;
    case VER_PLATFORM_WIN32_NT:
        if (winver->dwMajorVersion > 4) /* win2000 or higher */
	    old = 0;
        break;
    }
    free(winver);
    return old;
}

static int unmangle (const char *dosname, char *longname)
{
    if (old_windows()) {
	/* sorry but I really can't be bothered */
	strcpy(longname, dosname);
	return 0;
    } else {
	int err;
	void *handle;
	void (*real_unmangle)(const char *, char *, int, int *); 
	
	if (gui_open_plugin("longname", &handle)) return 1;

	real_unmangle = get_plugin_function("real_unmangle", handle);
	if (real_unmangle == NULL) return 1;
	(*real_unmangle)(dosname, longname, MAXLEN, &err);
	close_plugin(handle);
	return err;
    }
}

#endif

static void count_selections (GtkTreeModel *model, GtkTreePath *path,
			      GtkTreeIter *iter, int *selcount)
{
    *selcount += 1;
}

static gboolean main_popup_handler (GtkWidget *widget, GdkEvent *event,
				    gpointer data)
{
    GdkModifierType mods;
    GtkTreeSelection *select; 
    GtkMenu *menu = GTK_MENU(data);
    int selcount = 0;

    /* if no selection, don't do anything special */
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(widget));
    if (select == NULL) return FALSE;
    gtk_tree_selection_selected_foreach (select, 
                                         (GtkTreeSelectionForeachFunc) 
                                         count_selections,
                                         &selcount);
    if (selcount == 0) return FALSE;

    /* ignore all but right-clicks */
    gdk_window_get_pointer(widget->window, NULL, NULL, &mods);
    if (mods & GDK_BUTTON3_MASK && event->type == GDK_BUTTON_PRESS) {
        GdkEventButton *bevent = (GdkEventButton *) event; 

	if (selcount > 1) menu = GTK_MENU(selection_popup);
	gtk_menu_popup (menu, NULL, NULL, NULL, NULL,
			bevent->button, bevent->time);
	return (selcount > 1);
    }

    return FALSE;
}
