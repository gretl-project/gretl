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
#include "console.h"
#include "session.h"
#include "webget.h"
#include "database.h"
#include "datafiles.h"
#include "cmdstack.h"
#include "filelists.h"
#include "toolbar.h"
#include "menustate.h"

#include <dirent.h>

#ifndef G_OS_WIN32
# include <unistd.h>
# include "../pixmaps/gretl.xpm"  /* program icon for X */
#else
# include <windows.h>
# include "gretlwin32.h"
#endif

/* #define WINDEBUG 1 */

#ifdef WINDEBUG
FILE *dbg;
#endif

/* functions private to gretl.c */
static void sort_varlist (gpointer p, guint col, GtkWidget *w);
static GtkWidget *make_main_window (int gui_get_data);
#ifndef G_OS_WIN32
static void clip_init (GtkWidget *w);
#endif

static gboolean main_popup_handler (GtkWidget *w, GdkEventButton *event,
				    gpointer data);
static int selection_count (GtkTreeSelection *select, int *row);
static void set_up_main_menu (void);
static void startRcallback (gpointer p, guint opt, GtkWidget *w);
static void auto_store (void);
static void restore_sample_callback (gpointer p, int verbose, GtkWidget *w);

GtkTargetEntry gretl_drag_targets[] = {
    { "text/uri-list", 0, GRETL_FILENAME },
    { "db_pointer", GTK_TARGET_SAME_APP, GRETL_POINTER }    
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
static int opteng, optdump;

static const struct poptOption options[] = {
    { "run", 'r', POPT_ARG_STRING, &optrun, 0, 
      N_("open a script file on startup"), "SCRIPT" },
    { "db", 'd', POPT_ARG_STRING, &optdb, 0, 
      N_("open a database on startup"), "DATABASE" },
    { "webdb", 'w', POPT_ARG_STRING, &optdb, 0, 
      N_("open a remote (web) database on startup"), "REMOTE_DB" },
    { "english", 'e', POPT_ARG_NONE, &opteng, 0, 
      N_("force use of English"), NULL },
    { "dump", 'c', POPT_ARG_NONE, &optdump, 0, 
      N_("dump gretl configuration to file"), NULL },
    { NULL, '\0', 0, NULL, 0, NULL, NULL },
};
#endif /* USE_GNOME */

windata_t *mdata;
DATAINFO *datainfo;
char *errtext;
char cmdfile[MAXLEN], scriptfile[MAXLEN];
char trydatfile[MAXLEN], tryscript[MAXLEN];
char line[1024];
PATHS paths;                /* useful paths */
double **Z;                 /* data set */
MODEL **models;             /* gretl models structs */

int plot_count, data_status, orig_vars;
gchar *clipboard_buf; /* for copying models as LaTeX */
float gui_scale;

/* defaults for some options */
int expert = FALSE; 
int updater = FALSE;
int want_toolbar = TRUE;

char dbproxy[21];

#ifdef G_OS_WIN32
char Rcommand[MAXSTR] = "RGui.exe";
char calculator[MAXSTR] = "calc.exe";
char viewdvi[MAXSTR] = "windvi.exe";
char viewps[MAXSTR] = "gsview.exe";
#else
char calculator[MAXSTR] = "xcalc";
char viewdvi[MAXSTR] = "xdvi";
char viewps[MAXSTR] = "gv";
# ifdef USE_GNOME
char Rcommand[MAXSTR] = "R --gui=gnome";
extern const char *version_string;
# else
char Rcommand[MAXSTR] = "xterm -e R";
# endif
#endif

#ifdef HAVE_TRAMO
char tramo[MAXSTR] = "tramo";
char tramodir[MAXSTR] = "";
#endif

static void spreadsheet_edit (gpointer p, guint u, GtkWidget *w) 
{
    show_spreadsheet(NULL);
}

static void manual_update_query (gpointer p, guint u, GtkWidget *w)
{
    update_query();
}

#ifdef USE_GNOME
static void gnome_help (void)
{
    GError *error = NULL;

    gnome_help_display ("gretl.xml", NULL, &error);
        
    if (error != NULL) {
	g_warning (error->message);
	g_error_free (error);
    }
}
#endif 

#ifdef OSX_PKG
static void osx_help (void)
{
    char *prefix, *syscmd;
   
    prefix = getenv("GTK_EXE_PREFIX");
    if (prefix == NULL) {
        errbox("Couldn't find the manual");
	return;
    }
    
    syscmd = g_strdup_printf("%s/bin/manual.sh", prefix);
    system(syscmd);
    g_free(syscmd);
}
#endif 

extern void find_var (gpointer p, guint u, GtkWidget *w); /* gui_utils.c */

static void varinfo_callback (gpointer p, guint u, GtkWidget *w)
{
    varinfo_dialog(mdata->active_var, 1);
}

GtkItemFactoryEntry data_items[] = {

    /* File menu */
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },

    /* File, Open data */
    { N_("/File/_Open data"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Open data/user file..."), NULL, open_data, OPEN_DATA, 
      "<StockItem>", GTK_STOCK_OPEN },
    { N_("/File/Open data/sample file..."), "", display_files, TEXTBOOK_DATA, 
      "<StockItem>", GTK_STOCK_OPEN },
    { N_("/File/Open data/sep1"), NULL, NULL, 0, "<Separator>", GNULL },    
    { N_("/File/Open data/import CSV..."), NULL, open_data, OPEN_CSV, NULL, GNULL },
    { N_("/File/Open data/import ASCII..."), NULL, open_data, OPEN_ASCII, NULL, GNULL },
    { N_("/File/Open data/import BOX..."), NULL, open_data, OPEN_BOX, NULL, GNULL },
#ifdef G_OS_WIN32
    { N_("/File/Open data/import Excel..."), NULL, open_data, 
      OPEN_EXCEL, NULL, GNULL },
    { N_("/File/Open data/import Gnumeric..."), NULL, open_data, 
      OPEN_GNUMERIC, NULL, GNULL },
#else
    { N_("/File/Open data/import Gnumeric..."), NULL, open_data, 
      OPEN_GNUMERIC, NULL, GNULL },
    { N_("/File/Open data/import Excel..."), NULL, open_data, 
      OPEN_EXCEL, NULL, GNULL },
#endif

    /* File, Append data */
    { N_("/File/_Append data"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Append data/standard format..."), NULL, open_data, APPEND_DATA, NULL, GNULL },
    { N_("/File/Append data/from CSV..."), NULL, open_data, APPEND_CSV, NULL, GNULL },
    { N_("/File/Append data/from ASCII..."), NULL, open_data, APPEND_ASCII, NULL, GNULL },
    { N_("/File/Append data/from Gnumeric..."), NULL, open_data, 
      APPEND_GNUMERIC, NULL, GNULL },
    { N_("/File/Append data/from Excel..."), NULL, open_data, 
      APPEND_EXCEL, NULL, GNULL },

    /* File, Save data */
    { N_("/File/_Save data"), "<control>S", auto_store, 0, "<StockItem>", GTK_STOCK_SAVE },
    { N_("/File/Save data _as"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Save data as/_standard format..."), NULL, file_save, SAVE_DATA_AS, 
      "<StockItem>", GTK_STOCK_SAVE_AS },
    { N_("/File/Save data as/_gzipped..."), NULL, file_save, SAVE_GZDATA, 
      "<StockItem>", GTK_STOCK_SAVE_AS },

    /* File, Export data */
    { N_("/File/_Export data"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Export data/_CSV..."), NULL, file_save, EXPORT_CSV, NULL, GNULL },
    { N_("/File/Export data/GNU _R..."), NULL, file_save, EXPORT_R, NULL, GNULL },
    { N_("/File/Export data/GNU _octave..."), NULL, file_save, 
      EXPORT_OCTAVE, NULL, GNULL },
    { N_("/File/Export data/_PcGive..."), NULL, file_save, EXPORT_DAT, NULL, GNULL },
    { N_("/File/C_lear data set"), NULL, verify_clear_data, 0, 
      "<StockItem>", GTK_STOCK_CLEAR },
    { N_("/File/sep0"), NULL, NULL, 0, "<Separator>", GNULL },

    /* File, Browse databases */
    { N_("/File/_Browse databases"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Browse databases/_gretl native"), NULL, display_files, 
      NATIVE_DB, NULL, GNULL },
    { N_("/File/Browse databases/_RATS 4"), NULL, display_files, 
      RATS_DB, NULL, GNULL },
    { N_("/File/Browse databases/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/File/Browse databases/on database _server"), NULL, display_files, 
      REMOTE_DB, NULL, GNULL },

    /* File, Create dataset */
    { N_("/File/_Create data set"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Create data set/time-series"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Create data set/time-series/annual"), 
      NULL, newdata_callback, 1, NULL, GNULL },    
    { N_("/File/Create data set/time-series/quarterly"), 
      NULL, newdata_callback, 4, NULL, GNULL },    
    { N_("/File/Create data set/time-series/monthly"), 
      NULL, newdata_callback, 12, NULL, GNULL },
    { N_("/File/Create data set/time-series/high frequency"), 
      NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Create data set/time-series/high frequency/weekly"), 
      NULL, newdata_callback, 52, NULL, GNULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (5-day week)"), 
      NULL, newdata_callback, 5, NULL, GNULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (6-day week)"), 
      NULL, newdata_callback, 6, NULL, GNULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (7-day week)"), 
      NULL, newdata_callback, 7, NULL, GNULL }, 
    { N_("/File/Create data set/time-series/high frequency/hourly"), 
      NULL, newdata_callback, 24, NULL, GNULL }, 
    { N_("/File/Create data set/cross-sectional"), 
      NULL, newdata_callback, 0, NULL, GNULL }, 
    { N_("/File/Create data set/simulation"), NULL, gretl_callback, 
      NULLDATA, NULL, GNULL },
    { N_("/File/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/File/_View command log"), NULL, view_command_log, 0, NULL, GNULL },
    { N_("/File/sep2a"), NULL, NULL, 0, "<Separator>", GNULL },

    /* File, command files */
    { N_("/File/Open command file"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Open command file/user file..."), "", open_script, OPEN_SCRIPT, 
      "<StockItem>", GTK_STOCK_OPEN },
    { N_("/File/Open command file/practice file..."), "", display_files, PS_FILES, 
      "<StockItem>", GTK_STOCK_OPEN },
    { N_("/File/New command file"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/New command file/regular script"), NULL, do_new_script, 0, 
      "<StockItem>", GTK_STOCK_NEW },
    { N_("/File/New command file/Monte Carlo loop"), "", do_new_script, 1, 
      "<StockItem>", GTK_STOCK_NEW },
    { N_("/File/sep3"), NULL, NULL, 0, "<Separator>", GNULL },

    /* File, preferences */
    { N_("/File/_Preferences"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Preferences/_General..."), NULL, options_dialog, 0, 
      "<StockItem>", GTK_STOCK_PREFERENCES },
    { N_("/File/Preferences/_Fixed font..."), NULL, font_selector, 
      FIXED_FONT_SELECTION, "<StockItem>", GTK_STOCK_SELECT_FONT },
#ifndef USE_GNOME
    { N_("/File/Preferences/_Menu font..."), NULL, font_selector, 
      APP_FONT_SELECTION, "<StockItem>", GTK_STOCK_SELECT_FONT },
#endif
    { N_("/File/sep5"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/File/E_xit"), "<control>X", menu_exit_check, 0, "<StockItem>", GTK_STOCK_QUIT },

    /* Utilities menu */
    { N_("/_Utilities"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Utilities/Statistical tables"), NULL, stats_calculator, CALC_DIST, NULL, GNULL },
    { N_("/Utilities/p-value finder"), NULL, stats_calculator, CALC_PVAL, NULL, GNULL },
    { N_("/Utilities/Test statistic calculator"), NULL, stats_calculator, CALC_TEST, 
      NULL, GNULL },
    { N_("/Utilities/sep"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Utilities/Gretl console"), NULL, show_gretl_console, 0, NULL, GNULL },
    { N_("/Utilities/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Utilities/Start GNU R"), NULL, startRcallback, 0, NULL, GNULL },
    { "/Utilities/sep3", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Utilities/NIST test suite"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Utilities/NIST test suite/basic"), NULL, do_nistcheck, 0, NULL, GNULL },
    { N_("/Utilities/NIST test suite/verbose"), NULL, do_nistcheck, 1, NULL, GNULL },
    { N_("/Utilities/NIST test suite/very verbose"), NULL, do_nistcheck, 2, NULL, GNULL },

    /* Session menu */
    { N_("/_Session"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Session/_Icon view"), NULL, view_session, 0, NULL, GNULL },
    { N_("/Session/sep0"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Session/_Open..."), "", open_script, OPEN_SESSION, "<StockItem>", GTK_STOCK_OPEN },
    { N_("/Session/sep1"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Session/_Save"), "", save_session_callback, SAVE_AS_IS, 
      "<StockItem>", GTK_STOCK_SAVE },
    { N_("/Session/Save _as..."), "", save_session_callback, SAVE_RENAME, 
      "<StockItem>", GTK_STOCK_SAVE_AS },

    /* Data menu */
    { N_("/_Data"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/_Display values"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/Display values/_all variables"), NULL, display_data, 0, NULL, GNULL },
    { N_("/Data/Display values/_selected variables..."), 
      NULL, display_selected, 0, NULL, GNULL },
    { N_("/Data/_Edit values"), NULL, spreadsheet_edit, 0, NULL, GNULL },
    { "/Data/sepsort", NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/Sort variables"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/Sort variables/by ID number"), NULL, sort_varlist, 0, NULL, GNULL },
    { N_("/Data/Sort variables/by name"), NULL, sort_varlist, 1, NULL, GNULL },
    { N_("/Data/sep1"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/_Graph specified vars"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/Graph specified vars/Time series plot..."), 
      NULL, selector_callback, GR_PLOT, NULL, GNULL },
    { N_("/Data/Graph specified vars/X-Y scatter..."), 
      NULL, selector_callback, GR_XY, NULL, GNULL },
    { N_("/Data/Graph specified vars/X-Y with impulses..."), 
      NULL, selector_callback, GR_IMP, NULL, GNULL },
    { N_("/Data/Graph specified vars/X-Y with factor separation..."), 
      NULL, selector_callback, GR_DUMMY, NULL, GNULL },
    { N_("/Data/Graph specified vars/Boxplots..."), 
      NULL, gretl_callback, GR_BOX, NULL, GNULL },
    { N_("/Data/Graph specified vars/Notched boxplots..."), 
      NULL, gretl_callback, GR_NBOX, NULL, GNULL },
    { N_("/Data/Graph specified vars/3D plot..."), 
      NULL, selector_callback, GR_3D, NULL, GNULL },
    { N_("/Data/_Multiple scatterplots..."), 
      NULL, selector_callback, SCATTERS, NULL, GNULL },
    { N_("/Data/sep2"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/_Read info"), NULL, open_info, 0, NULL, GNULL },
    { N_("/Data/Edit _info"), NULL, edit_header, 0, NULL, GNULL },
    { N_("/Data/Print description"), NULL, print_report, 0, NULL, GNULL },
    { N_("/Data/sep3"), NULL, NULL, 0, "<Separator>", NULL },

    { N_("/Data/_Summary statistics"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/_Summary statistics/_all variables"), NULL, 
      do_menu_op, SUMMARY, NULL, GNULL },
    { N_("/Data/_Summary statistics/_selected variables"), NULL, 
      do_menu_op, SUMMARY_SELECTED, NULL, GNULL },

    { N_("/Data/_Correlation matrix"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/_Correlation matrix/_all variables"), 
      NULL, do_menu_op, CORR, NULL, GNULL },
    { N_("/Data/_Correlation matrix/_selected variables"), NULL, do_menu_op, 
      CORR_SELECTED, NULL, GNULL },

    { N_("/Data/_Principal components"), NULL, do_menu_op, PCA, NULL, GNULL },

    { N_("/Data/sep4"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/Difference of means"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/Difference of means/assuming equal variances..."), NULL, 
      gretl_callback, MEANTEST, NULL, GNULL },
    { N_("/Data/Difference of means/assuming unequal variances..."), NULL, 
      gretl_callback, MEANTEST2, NULL, GNULL },
    { N_("/Data/Difference of variances..."), NULL, gretl_callback, VARTEST, NULL, GNULL },
    { N_("/Data/sep5"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/Add variables"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Data/Add variables/time trend"), NULL, add_index, 1, NULL, GNULL },
    { N_("/Data/Add variables/index variable"), NULL, add_index, 0, NULL, GNULL },
    { N_("/Data/Add variables/logs of selected variables"), NULL, 
      add_logs_etc, LOGS, NULL, GNULL },
    { N_("/Data/Add variables/lags of selected variables"), NULL, 
      add_logs_etc, LAGS, NULL, GNULL },
    { N_("/Data/Add variables/squares of selected variables"), NULL, 
      add_logs_etc, SQUARE, NULL, GNULL },
    { N_("/Data/Add variables/first differences of selected variables"), NULL, 
      add_logs_etc, DIFF, NULL, GNULL },
    { N_("/Data/Add variables/log differences of selected variables"), NULL, 
      add_logs_etc, LDIFF, NULL, GNULL },
    { N_("/Data/Add variables/periodic dummies"), NULL, add_dummies, 0, NULL, GNULL },
    { N_("/Data/Add variables/unit dummies"), NULL, add_dummies, 1, NULL, GNULL },
    { N_("/Data/Add variables/panel dummies"), NULL, add_dummies, 2, NULL, GNULL },
    { N_("/Data/Add variables/sep"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/Add variables/random normal..."), NULL, 
      add_random_callback, GENR_NORMAL, NULL, GNULL },
    { N_("/Data/Add variables/random uniform..."), NULL, 
      add_random_callback, GENR_UNIFORM, NULL, GNULL },
    { N_("/Data/Add variables/seed generator..."), NULL, gretl_callback, 
      SETSEED, NULL, GNULL },
    { N_("/Data/Add variables/sep2"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Data/Add variables/Define _new variable..."), NULL, gretl_callback, 
      GENR, NULL, GNULL },
    { N_("/Data/Refresh window"), NULL, refresh_data, 0, NULL, GNULL },

    /* Sample menu */
    { N_("/_Sample"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Sample/_Set range..."), NULL, sample_range_dialog, SMPL, NULL, GNULL },
    { N_("/Sample/_Restore full range"), NULL, restore_sample_callback, 1, NULL, GNULL },
    { N_("/Sample/sep1"), NULL, NULL, 0, "<Separator>", NULL }, 
    { N_("/Sample/Dataset structure..."), NULL, data_structure_wizard, 0, NULL, GNULL },
    { N_("/Sample/Compact data..."), NULL, do_compact_data_set, 0, NULL, GNULL },
    { N_("/Sample/sep2"), NULL, NULL, 0, "<Separator>", NULL },   
    { N_("/Sample/_Define, based on dummy..."), NULL, sample_range_dialog, 
      SMPLDUM, NULL, GNULL },
    { N_("/Sample/_Restrict, based on criterion..."), NULL, gretl_callback, 
      SMPLBOOL, NULL, GNULL },
    { N_("/Sample/R_andom sub-sample..."), NULL, sample_range_dialog, SMPLRAND, NULL, GNULL },
    { N_("/Sample/sep3"), NULL, NULL, 0, "<Separator>", NULL },  
    { N_("/Sample/Drop all obs with _missing values"), NULL, drop_all_missing, 
      0, NULL, GNULL },
    { N_("/Sample/_Count missing values"), NULL, count_missing, 0, NULL, GNULL },
    { N_("/Sample/Set missing _value code..."), NULL, gretl_callback, 
      GSETMISS, NULL, GNULL },
    { N_("/Sample/sep4"), NULL, NULL, 0, "<Separator>", NULL },  
    { N_("/Sample/_Add case markers..."), NULL, gretl_callback, MARKERS, NULL, GNULL },
    { N_("/Sample/Remove case _markers"), NULL, do_remove_markers, 0, NULL, GNULL },
    { N_("/Sample/sep5"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Sample/Restructure panel..."), NULL, panel_restructure_dialog, 0, NULL, GNULL },
    { N_("/Sample/Transpose data..."), NULL, gui_transpose_data, 0, NULL, GNULL },

    /* Variable menu */
    { N_("/_Variable"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Variable/Find..."), NULL, find_var, 0, "<StockItem>", GTK_STOCK_FIND },
    { N_("/Variable/_Display values"), NULL, display_var, 0, NULL, GNULL },
    { N_("/Variable/_Summary statistics"), NULL, do_menu_op, 
      VAR_SUMMARY, NULL, GNULL },
    { N_("/Variable/_Frequency distribution"), NULL, do_menu_op, 
      FREQ, NULL, GNULL },
    { N_("/Variable/Frequency plot"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Variable/Frequency plot/simple"), NULL, do_freqplot, 0, NULL, GNULL },
    { N_("/Variable/Frequency plot/against Normal"), NULL, do_freqplot, 
      NORMAL, NULL, GNULL },
    { N_("/Variable/Frequency plot/against Gamma"), NULL, do_freqplot, 
      GAMMA, NULL, GNULL },
    { N_("/Variable/Estimated density plot..."), NULL, do_kernel, 0, NULL, GNULL },
    { N_("/Variable/Range-mean graph"), NULL, do_range_mean, 0, NULL, GNULL }, 
    { N_("/Variable/sep1"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Variable/_Time series plot"), NULL, ts_plot_var, 0, NULL, GNULL },
    { N_("/Variable/Correlogram"), NULL, gretl_callback, CORRGM, NULL, GNULL },
    { N_("/Variable/Spectrum"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Variable/Spectrum/sample periodogram"), NULL, do_pergm, 0, NULL, GNULL }, 
    { N_("/Variable/Spectrum/Bartlett lag window"), NULL, do_pergm, 1, NULL, GNULL }, 
    { N_("/Variable/_Augmented Dickey-Fuller test"), NULL, unit_root_test, ADF, NULL, GNULL },
    { N_("/Variable/_KPSS test"), NULL, unit_root_test, KPSS, NULL, GNULL },
    { N_("/Variable/ARMA model"), NULL, arma_options_dialog, 0, NULL, GNULL },
#ifdef HAVE_X12A
    { N_("/Variable/X-12-ARIMA analysis"), NULL, do_tramo_x12a, X12A, NULL, GNULL },
#endif
#ifdef HAVE_TRAMO
    { N_("/Variable/TRAMO analysis"), NULL, do_tramo_x12a, TRAMO, NULL, GNULL },
#endif
    { N_("/Variable/Runs test"), NULL, do_menu_op, RUNS, NULL, GNULL }, 
    { N_("/Variable/sep2"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Variable/_Edit attributes"), NULL, varinfo_callback, 0, NULL, GNULL },
    { N_("/Variable/Set missing value code..."), NULL, gretl_callback, 
      VSETMISS, NULL, GNULL },
    { N_("/Variable/sep3"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Variable/Define _new variable..."), NULL, gretl_callback, GENR, NULL, GNULL },

    /* Model menu */
    { N_("/_Model"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Model/_Ordinary Least Squares..."), NULL, model_callback, OLS, NULL, GNULL },
    { N_("/Model/_Weighted Least Squares..."), NULL, model_callback, WLS, NULL, GNULL },
    { N_("/Model/sep1"),  NULL, NULL, 0, "<Separator>", GNULL },
#if 0
    { N_("/Model/HCC_M..."), NULL, model_callback, HCCM, NULL, GNULL },
#endif
    { N_("/Model/H_eteroskedasticity corrected..."), NULL, model_callback, 
      HSK, NULL, GNULL },
    { N_("/Model/sep2"),  NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Model/Time series"), NULL, NULL, 0, "<Branch>", NULL },
    { N_("/Model/sep3"),  NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Model/_Two-Stage Least Squares..."), NULL, model_callback, TSLS, NULL, GNULL },
    { "/Model/sep4",  NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Model/_Logit..."), NULL, model_callback, LOGIT, NULL, GNULL },
    { N_("/Model/_Probit..."), NULL, model_callback, PROBIT, NULL, GNULL },
    { N_("/Model/To_bit..."), NULL, model_callback, TOBIT, NULL, GNULL },
    { "/Model/sep5",  NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Model/Lo_gistic..."), NULL, model_callback, LOGISTIC, NULL, GNULL },
    { N_("/Model/Least _Absolute Deviation..."), NULL, model_callback, LAD, NULL, GNULL },
    { N_("/Model/_Rank correlation..."), NULL, gretl_callback, SPEARMAN, NULL, GNULL },
    { N_("/Model/_Pooled OLS (panel)..."), NULL, model_callback, POOLED, NULL, GNULL },
    { N_("/Model/Nonlinear Least Squares..."), NULL, gretl_callback, NLS, NULL, GNULL },
#ifdef ENABLE_GMP
    { N_("/Model/High precision OLS..."), NULL, mp_ols_callback, MPOLS, NULL, GNULL },
#endif

    /* Help menu */
    { N_("/_Help"), NULL, NULL, 0, "<LastBranch>", NULL },
    { N_("/Help/_GUI commands"), NULL, do_gui_help, 0, NULL, GNULL },
    { N_("/Help/_Script commands syntax"), NULL, do_script_help, 1, NULL, GNULL },
    { N_("/Help/sep1"), NULL, NULL, 0, "<Separator>", NULL },
#if defined(USE_GNOME)
    { N_("/Help/Manual in HTML"), NULL, gnome_help, 0, NULL, GNULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
#elif defined(G_OS_WIN32)
    { N_("/Help/Manual in HTML"), NULL, win_help, 0, NULL, GNULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
#elif defined(OSX_PKG)
    { N_("/Help/Manual in PDF"), NULL, osx_help, 0, NULL, GNULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
#endif
    { N_("/Help/_Check for updates"), NULL, manual_update_query, 0, NULL, GNULL },
    { N_("/Help/sep3"), NULL, NULL, 0, "<Separator>", NULL },
    { N_("/Help/_About gretl"), NULL, about_dialog, 0, NULL, GNULL }
};

static void gui_usage (void)
{
    gui_logo(stdout);
    printf(I_("You may supply the name of a data file on the command line.\n"));
    printf(I_("Or you may do \"gretl -r script_file\" to open a script.\n"));
    printf(I_("Or you may do \"gretl -d database\" to open a gretl database.\n"));
    printf(I_("You may do \"gretl -e\" to force use of English.\n"));
    exit(0);
}

static void noalloc (const char *str)
{
    fprintf(stderr, I_("Couldn't allocate memory for %s\n"), str);
    exit(EXIT_FAILURE);
}

static void get_runfile (char *fname)
{
    int i;

    *tryscript = '\0';
#ifdef G_OS_WIN32
    if (unmangle(fname, tryscript)) return;
#else
    strncat(tryscript, fname, MAXLEN-1);
#endif
    if (addpath(tryscript, &paths, 1) == NULL) {
	fprintf(stderr, I_("Couldn't find script '%s'\n"), tryscript);
	exit(EXIT_FAILURE);
    } else {
	fprintf(stderr, I_("%s found\n"), tryscript);
	i = slashpos(tryscript);
	if (i) {
	    paths.currdir[0] = '\0';
	    strncat(paths.currdir, tryscript, i);
	}
	strcat(paths.currdir, SLASHSTR);
    }
}

static void fix_dbname (char *db)
{
    FILE *fp = NULL;

    if (strstr(db, ".bin") == NULL &&
	strstr(db, ".rat") == NULL) {
	strcat(db, ".bin");
    }

    fp = fopen(db, "rb");

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

#ifdef ENABLE_NLS

# if defined(G_OS_WIN32)

static void real_nls_init (void)
{
    char gretldir[MAXSTR], localedir[MAXSTR];
    char *loc;

    if (read_reg_val(HKEY_CLASSES_ROOT, "gretl", "gretldir", gretldir)) {
	return;
    }

    build_path(gretldir, "locale", localedir, NULL);
    loc = setlocale (LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain (PACKAGE, localedir);
    textdomain (PACKAGE);
    bind_textdomain_codeset (PACKAGE, "UTF-8");
}

# elif defined(OSX_PKG)

static void real_nls_init (void)
{
    char *prefix = getenv("GTK_EXE_PREFIX");
    char *localedir;
    char *loc;

    if (prefix == NULL) return;
    
    localedir = g_strdup_printf("%s/share/locale", prefix);
    loc = setlocale (LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain (PACKAGE, localedir);
    textdomain (PACKAGE);
    bind_textdomain_codeset (PACKAGE, "UTF-8");
    g_free(localedir);
}

# else /* regular *nix treatment */

static void real_nls_init (void)
{
    char *loc;

    loc = setlocale (LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain (PACKAGE, LOCALEDIR);
    textdomain (PACKAGE);
    bind_textdomain_codeset (PACKAGE, "UTF-8");
}

# endif

void nls_init (void)
{
    char *mylang = getenv("GRETL_LANG");

    if (mylang != NULL) {
	if (!g_ascii_strcasecmp(mylang, "english") ||
	    !g_ascii_strcasecmp(mylang, "C")) return;
    } 

    real_nls_init();
}

static void force_english (void)
{
    setlocale (LC_ALL, "C");

# ifdef G_OS_WIN32
    SetEnvironmentVariable("LC_ALL", "C");
    putenv("LC_ALL=C");
    textdomain("none");
# endif

    force_english_help();
}

#endif /* ENABLE_NLS */


int main (int argc, char *argv[])
{
    int err = 0, gui_get_data = 0;
    char dbname[MAXLEN];
    char filearg[MAXLEN];
#ifdef USE_GNOME
    GnomeProgram *program;
#endif

#ifdef WINDEBUG
    dbg = fopen("debug.txt", "w");
#endif    

#ifdef ENABLE_NLS
    nls_init();
#endif  

#ifdef G_OS_WIN32
    putenv("PANGO_WIN32_NO_UNISCRIBE=a");
    putenv("G_FILENAME_ENCODING=@locale");
#endif     

    if ((errtext = malloc(MAXLEN)) == NULL) {
	noalloc(_("startup"));
    }

    *tryscript = '\0';
    *scriptfile = '\0';
    *paths.datfile = '\0';
    *dbname = '\0';

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
#endif /* USE_GNOME */

    set_paths(&paths, OPT_D | OPT_X); /* defaults, gui */

#ifdef G_OS_WIN32
    gretl_win32_init(argv[0]);
#else 
    set_rcfile(); /* also calls read_rc() */
#endif/* G_OS_WIN32 */

    if (argc > 1) {
	int english = 0;
	int opt = parseopt((const char **) argv, argc, filearg, &english);

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
	    if (*filearg == '\0') gui_usage();
	    get_runfile(filearg);
#endif
	    gui_get_data = 1;
	    break;
	case OPT_DBOPEN:
	case OPT_WEBDB:
#ifdef USE_GNOME
	    strncpy(dbname, optdb, MAXLEN-1);
#else
	    if (*filearg == '\0') gui_usage();
	    strcpy(dbname, filearg);
#endif
	    if (opt == OPT_DBOPEN) fix_dbname(dbname);
	    gui_get_data = opt;
	    break;
	case OPT_DUMP:
	    dump_rc();
	    exit(EXIT_SUCCESS);
	    break;
	default:
	    /* let's suppose the argument is a data file */
	    break;
	}

#ifdef ENABLE_NLS
	if (english) {
	    force_english();
	    if (argc == 2) {
		gui_get_data = 1;
	    }	
	}
#endif
    } else {
	gui_get_data = 1;
    }

    strcpy(cmdfile, paths.userdir);
    strcat(cmdfile, "session.inp");

    /* allocate data information struct */
    datainfo = datainfo_new();
    if (datainfo == NULL) {
	noalloc(_("data information"));
    }

    /* allocate memory for models */
    models = malloc(3 * sizeof *models);
    if (models == NULL) noalloc(_("models")); 

    models[0] = gretl_model_new();
    models[1] = gretl_model_new();
    models[2] = gretl_model_new();

    if (models[0] == NULL || models[1] == NULL || models[2] == NULL) 
	noalloc(_("models")); 

    library_command_init();

    helpfile_init();
    session_init();
    init_fileptrs();

    /* get the data file, if specified on the command line */
    if (!gui_get_data) {
	int ftype;
	PRN *prn; 

	prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
	if (prn == NULL) exit(EXIT_FAILURE);

	*paths.datfile = '\0';

#ifdef G_OS_WIN32
	if (unmangle(filearg, paths.datfile)) 
	    exit(EXIT_FAILURE);
#else
	strcpy(paths.datfile, filearg);
#endif
	/* record the name the user supplied */
	strcpy(trydatfile, paths.datfile);

	ftype = detect_filetype(paths.datfile, &paths, prn);

	switch (ftype) {
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(&Z, &datainfo, paths.datfile, &paths, DATA_NONE, 
				 prn);
	    break;
	case GRETL_XML_DATA:
	    err = get_xmldata(&Z, &datainfo, paths.datfile, &paths, DATA_NONE, 
			      prn, 0);
	    break;
	case GRETL_CSV_DATA:
	    err = import_csv(&Z, &datainfo, paths.datfile, &paths, prn);
	    break;
	case GRETL_BOX_DATA:
	    err = import_box(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_EXCEL:
	case GRETL_GNUMERIC:
	    err = get_worksheet_data(paths.datfile, ftype, 0, &gui_get_data);
	    break;
	case GRETL_SCRIPT:
	    gui_get_data = 1;
	    get_runfile(paths.datfile);
	    *paths.datfile = '\0';
	    break;
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:    
	    strcpy(dbname, paths.datfile);
	    *paths.datfile = '\0';
	    fix_dbname(dbname);
	    gui_get_data = OPT_DBOPEN;
	    break;
	case GRETL_UNRECOGNIZED:
	    gui_usage();
	default:
	    exit(EXIT_FAILURE);
	    break;
	}

#ifdef notdef
	if (paths.datfile[0] != 0) {
	    my_filename_to_utf8(paths.datfile);
	}
#endif

	if (ftype != GRETL_SCRIPT && err) {
	    errmsg(err, prn);
	    exit(EXIT_FAILURE);
	}
	gretl_print_destroy(prn);
    }

#ifdef WINDEBUG
    fprintf(dbg, "starting on GUI building\n");
    fclose(dbg);
#endif

    /* create the GUI */
    gretl_tooltips_init();

    /* create main window */
    if ((mdata = mymalloc(sizeof *mdata)) == NULL)
	noalloc(_("GUI"));
    if (make_main_window(gui_get_data) == NULL) 
	noalloc(_("main window"));
    if (!gui_get_data) {
	set_sample_label(datainfo);
    }

#ifndef G_OS_WIN32
    /* Let a first-time user set the working dir */
    first_time_set_user_dir(); 

    /* enable special copying to clipboard */
    clip_init(mdata->w);
#endif

    add_files_to_menus();

    session_menu_state(FALSE);
    restore_sample_state(FALSE);
    main_menubar_state(FALSE);

#ifdef HAVE_TRAMO
    set_tramo_ok(-1);
#endif
#ifdef HAVE_X12A
    set_x12a_ok(-1);
#endif

    if (!gui_get_data) {
	register_data(paths.datfile, trydatfile, 1);
	*trydatfile = 0;
    }

    /* opening a script from the command line? */
    if (tryscript[0] != '\0') { 
	do_open_script();
    }

    /* check for program updates? */
    proxy_init(dbproxy);
    if (updater) {
	silent_update_query(); 
    }

    /* try opening specified database */
    if (gui_get_data == OPT_DBOPEN) {
	open_named_db_list(dbname);
    }
    else if (gui_get_data == OPT_WEBDB) {
	open_named_remote_db_list(dbname);
    }

    /* Enter the event loop */
    gtk_main();

    /* clean up before exiting */
    free_session();

    if (Z != NULL) {
	free_Z(Z, datainfo);
    }

    free_model(models[0]);
    free_model(models[1]);
    free_model(models[2]);
    free(models);

    library_command_free();

    if (data_status) {
	free_datainfo(datainfo);
    }

    free_command_stack();
    exit_free_modelspec();

    destroy_file_collections();

    return EXIT_SUCCESS;
}

/* ........................................................... */

static void check_varmenu_state (GtkTreeSelection *select, gpointer p)
{
    if (mdata->ifac != NULL) {
	int selcount = selection_count(select, NULL);

	flip(mdata->ifac, "/Variable", (selcount == 1));
	flip(mdata->ifac, "/Data/Correlation matrix/selected variables", 
	     (selcount > 1));
	flip(mdata->ifac, "/Data/Principal components", 
	     (selcount > 1));
    }
}

/* ........................................................... */

static gint catch_mdata_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    if (key->keyval == GDK_a) {
	audio_render_window(vwin, AUDIO_LISTBOX);
    } else if (key->keyval == GDK_x) {
	audio_render_window(NULL, AUDIO_LISTBOX);
    }
#endif

    if (key->keyval == GDK_e || key->keyval == GDK_t) {
	int selcount, row;

	selcount = 
	    selection_count(gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox)),
			    &row);
	if (selcount == 1 && row != 0) {
	    gchar *varnum;

	    tree_view_get_string(GTK_TREE_VIEW(mdata->listbox), row, 0, &varnum);
	    mdata->active_var = atoi(varnum);
	    g_free(varnum);

	    if (key->keyval == GDK_e) {
		varinfo_dialog(mdata->active_var, 1);
	    } 
	    else if (key->keyval == GDK_t) {
		do_graph_var(mdata->active_var);
	    } 
	}
    } 

    return FALSE;
}

/* ........................................................... */

void populate_varlist (void)
{
    GtkListStore *store;
    GtkTreeSelection *select;
    GtkTreeIter iter;    
    char id[4];
    gint i;
    static gint check_connected;
    static gint click_connected;

    /* find and clear the existing list */
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox)));
    gtk_list_store_clear(store);

    gtk_tree_model_get_iter_first (GTK_TREE_MODEL(store), &iter);

    for (i=0; i<datainfo->v; i++) {
	if (hidden_var(i, datainfo)) continue;
	gtk_list_store_append(store, &iter);
	sprintf(id, "%d", i);
	gtk_list_store_set (store, &iter, 
			    0, id, 
			    1, datainfo->varname[i],
			    2, VARLABEL(datainfo, i),
			    -1);
    } 

    mdata->active_var = 1;

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_select_iter(select, &iter);

    if (datainfo->v > 1) {
	GtkTreePath *path = gtk_tree_path_new_from_string("1");

	gtk_tree_view_set_cursor(GTK_TREE_VIEW(mdata->listbox), path, NULL, FALSE);
	gtk_tree_path_free(path);
    }

    if (!check_connected) {
	g_signal_connect (G_OBJECT(select), "changed",
			  G_CALLBACK(check_varmenu_state),
			  mdata);
	check_connected = 1;
    }

    if (!click_connected) {
	g_signal_connect (G_OBJECT(mdata->listbox), "button_press_event",
			  G_CALLBACK(main_varclick),
			  mdata);
	g_signal_connect(G_OBJECT(mdata->listbox), "button_press_event",
			 G_CALLBACK(main_popup_handler), 
			 mdata);
	g_signal_connect (G_OBJECT(mdata->listbox), "key_press_event",
			  G_CALLBACK(catch_mdata_key),
			  mdata);
	click_connected = 1;
    }
}

/* ......................................................... */

static gint 
compare_var_ids (GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		 gpointer p)
{
    gchar *t1, *t2;
    int i1, i2;

    gtk_tree_model_get(model, a, 0, &t1, -1);
    gtk_tree_model_get(model, b, 0, &t2, -1);

    i1 = atoi(t1);
    i2 = atoi(t2);

    g_free(t1);
    g_free(t2);

    return ((i1 < i2) ? -1 : (i1 > i2) ? 1 : 0);    
}

static gint 
compare_varnames (GtkTreeModel *model, GtkTreeIter *a, GtkTreeIter *b,
		  gpointer p)
{
    gchar *t1, *t2;
    gint ret;

    gtk_tree_model_get(model, a, 1, &t1, -1);
    gtk_tree_model_get(model, b, 1, &t2, -1);

    if (!strcmp(t1, "const")) return 0;
    if (!strcmp(t2, "const")) return 1;

    ret = strcmp(t1, t2);
    g_free(t1);
    g_free(t2);

    return ret;    
}

static void sort_varlist (gpointer p, guint col, GtkWidget *w)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));

    gtk_tree_sortable_set_sort_func (GTK_TREE_SORTABLE(model), 0,
				     compare_var_ids, NULL, NULL);
    gtk_tree_sortable_set_sort_func (GTK_TREE_SORTABLE(model), 1,
				     compare_varnames, NULL, NULL);
    gtk_tree_sortable_set_sort_column_id (GTK_TREE_SORTABLE(model), 
					  col, GTK_SORT_ASCENDING);
}

/* ......................................................... */

void clear_varlist (GtkWidget *widget)
{
    GtkListStore *store;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(widget)));
    gtk_list_store_clear(store);
}

/* ......................................................... */

static float get_gui_scale (void)
{
    GtkSettings *settings;
    gchar *fontname = NULL;
    int fsize;
    float scale = 1.0;

    settings = gtk_settings_get_default();

    g_object_get(G_OBJECT(settings), "gtk-font-name", &fontname, NULL);

    if (fontname != NULL) {
	if (sscanf(fontname, "%*s %d", &fsize) == 1) {
	    if (fsize > 10.0) {
		scale = fsize / 10.0;
	    }
	}
	g_free(fontname);
    }

    return scale;
}

/* ......................................................... */

static GtkWidget *make_main_window (int gui_get_data) 
{
    GtkWidget *main_vbox;
    GtkWidget *box, *dlabel, *align;
    const char *titles[] = {
	_("ID #"), 
	_("Variable name"), 
	_("Descriptive label")
    };
    int mainwin_width = 550;  
    int mainwin_height = 420; 

    gui_scale = get_gui_scale();

    mainwin_width *= gui_scale;
    mainwin_height *= gui_scale;

    mdata->data = NULL;  
    mdata->listbox = NULL;
    mdata->popup = NULL;
    mdata->dialog = NULL;
    mdata->role = MAINWIN;

#ifdef USE_GNOME
    mdata->w = gnome_app_new("gretl", _("Econometrics program"));
#else
    mdata->w = gtk_window_new(GTK_WINDOW_TOPLEVEL);
#endif

#ifdef G_OS_WIN32
    set_up_windows_look();
#endif

    g_signal_connect (G_OBJECT (mdata->w), "delete_event",
		      G_CALLBACK (exit_check), NULL);
    g_signal_connect (G_OBJECT (mdata->w), "destroy",
		      G_CALLBACK (destroy), NULL);

    gtk_window_set_title(GTK_WINDOW (mdata->w), "gretl");
    gtk_window_set_default_size(GTK_WINDOW (mdata->w), 
				mainwin_width, mainwin_height);
#ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(mdata->w), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif

    main_vbox = gtk_vbox_new(FALSE, 4);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 8);

#ifdef USE_GNOME
    gnome_app_set_contents(GNOME_APP(mdata->w), main_vbox);
#else
    gtk_container_add(GTK_CONTAINER(mdata->w), main_vbox);
#endif

    g_object_set_data(G_OBJECT(mdata->w), "vbox", main_vbox);

    set_up_main_menu();
    gtk_box_pack_start(GTK_BOX(main_vbox), mdata->mbar, FALSE, TRUE, 0);
    gtk_widget_show(mdata->mbar);

    dlabel = gtk_label_new(_(" No datafile loaded ")); 
    gtk_widget_show(dlabel);

    g_object_set_data(G_OBJECT(mdata->w), "dlabel", dlabel);

    box = gtk_vbox_new (FALSE, 0);
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(box), align, FALSE, FALSE, 0);
    gtk_widget_show(align);
    gtk_container_add(GTK_CONTAINER(align), dlabel);
   
    mdata->listbox = list_box_create (mdata, GTK_BOX(box), 3, 0, titles);
    gtk_widget_show(box);

    gtk_drag_dest_set (mdata->listbox,
		       GTK_DEST_DEFAULT_ALL,
		       gretl_drag_targets, 2,
		       GDK_ACTION_COPY);

    g_signal_connect (G_OBJECT(mdata->listbox), "drag_data_received",
		      G_CALLBACK(drag_data_received),
		      NULL);

    gtk_box_pack_start(GTK_BOX(main_vbox), box, TRUE, TRUE, 0);

    mdata->status = gtk_label_new("");
    
    gtk_box_pack_start (GTK_BOX(main_vbox), mdata->status, FALSE, TRUE, 0);

    /* put stuff into list box, activate menus */
    if (!gui_get_data) {
	populate_varlist();
    }

    /* get a monospaced font for various windows */
    set_fixed_font();

    /* and a proportional font for menus, etc */
#ifndef USE_GNOME
    set_app_font(NULL);
#endif

    gtk_widget_show_all(mdata->w); 

    /* create gretl toolbar? */
    show_or_hide_toolbar(want_toolbar);

    return main_vbox;
}

/* ........................................................... */

static void set_up_main_menu (void)
{
    GtkAccelGroup *accel_group;
    gint n_items = sizeof data_items / sizeof data_items[0];

    accel_group = gtk_accel_group_new ();
    mdata->ifac = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>", 
					accel_group);
    gtk_window_add_accel_group (GTK_WINDOW (mdata->w), accel_group);
    g_object_unref (accel_group);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(mdata->ifac, menu_translate, NULL, NULL);
#endif    
    gtk_item_factory_create_items (mdata->ifac, n_items, data_items, NULL);
    mdata->mbar = gtk_item_factory_get_widget (mdata->ifac, "<main>");
}

/* ........................................................... */

int restore_sample (gretlopt opt)
{
    int err;

    err = restore_full_sample(&Z, &datainfo, opt);
    if (err) {
	gui_errmsg(err);
    } else {
	restore_sample_state(FALSE);
    }

    return err;
}

static void restore_sample_callback (gpointer p, int verbose, GtkWidget *w)
{
    restore_sample(OPT_C); 

    if (verbose) {
	infobox(_("Full sample range restored"));
	set_sample_label(datainfo);    
	restore_sample_state(FALSE);
	strcpy(line, "smpl full");
	verify_and_record_command(line);
    }
}

static void startRcallback (gpointer p, guint opt, GtkWidget *w)
{
    startR(Rcommand);
}

/* ........................................................... */

#ifndef G_OS_WIN32

int gretl_fork (const char *prog, const char *arg)
{
    gchar *argv[3];
    gboolean run;
    
    argv[0] = g_strdup(prog);
    if (arg != NULL) {
	argv[1] = g_strdup(arg);
	argv[2] = NULL;
    } else {
	argv[1] = NULL;
    }
    
    run = g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH, 
			NULL, NULL, NULL, NULL);

    if (!run) {
	sprintf(errtext, "%s: %s", _("Command failed"), prog);
	errbox(errtext);
    }

    g_free(argv[0]);
    g_free(argv[1]);

    return !run;
}

#endif	

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

    /* handle drag of pointer from database window */
    if (info == GRETL_POINTER && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	import_db_series(*(void **) data->data);
	return;
    }

    /* ignore the wrong sort of data */
    if (data == NULL || (dfname = data->data) == NULL || 
	strlen(dfname) <= 5 || strncmp(dfname, "file:", 5))
	return;

    if (strncmp(dfname, "file://", 7) == 0) skip = 7;
#ifdef G_OS_WIN32
    if (strncmp(dfname, "file:///", 8) == 0) skip = 8;
#endif

    /* there may be multiple files: we ignore all but the first */
    *tmp = 0;
    if ((pos = haschar('\r', dfname)) > 0 || 
	(pos = haschar('\n', dfname) > 0)) {
	strncat(tmp, dfname + skip, pos - skip);
    } else {
	strcat(tmp, dfname + skip);
    }

    /* handle spaces in filenames */
    unescape_url(tmp);

#ifdef G_OS_WIN32
    slash_convert(tmp, TO_BACKSLASH);
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

#ifndef G_OS_WIN32

static void gretl_clipboard_get (GtkClipboard *clip,
				 GtkSelectionData *selection_data,
				 guint info,
				 gpointer p)
{
    gchar *str;
    gint length;

    str = clipboard_buf; /* global */
    if (str == NULL) return;
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
	GdkAtom seltype;
	gint format;
	gint new_length;

	c = str[length];
	str[length] = '\0';
	gdk_string_to_compound_text(str, &seltype, &format, 
				    &text, &new_length);
	gtk_selection_data_set(selection_data, seltype, format, 
			       text, new_length);
	gdk_free_compound_text(text);
	str[length] = c;
    }
}

/* ........................................................... */

static void gretl_clipboard_clear (GtkClipboard *clip, gpointer p)
{
    free(clipboard_buf);
    clipboard_buf = NULL;
}

/* ........................................................... */

static void clip_init (GtkWidget *w)
{
    GtkClipboard *clip;
    GtkTargetEntry targets[] = {
	{ "STRING", 0, TARGET_STRING },
	{ "TEXT",   0, TARGET_TEXT }, 
	{ "COMPOUND_TEXT", 0, TARGET_COMPOUND_TEXT }
    };

    gint n_targets = sizeof targets / sizeof targets[0];

    clip = gtk_clipboard_get(GDK_SELECTION_CLIPBOARD);

    if (!gtk_clipboard_set_with_owner(clip,
				      targets, n_targets,
				      gretl_clipboard_get,
				      gretl_clipboard_clear,
				      G_OBJECT(w))) {
	fprintf(stderr, "Failed to initialize clipboard\n");
    }
}

#endif /* G_OS_WIN32 */

/* ........................................................... */

static int native_datafile (void)
{
    int n = strlen(paths.datfile);
    
    if (n > 4) {
	if (!strcmp(paths.datfile + n - 4, ".gdt")) {
	    return 1;
	}
	if (using_olddat() && !strcmp(paths.datfile + n - 4, ".dat")) {
	    return 1;
	}
    }

    return 0;
}

static void auto_store (void)
{
    gretlopt oflag = OPT_NONE;

    /* if there's already a datafile, and it's gzipped, then
       arrange for the new store to be gzipped too */
    if (*paths.datfile && is_gzipped(paths.datfile)) {
	oflag = OPT_Z;
    }

    if ((data_status & USER_DATA) && native_datafile()) {
	do_store(paths.datfile, oflag, 1);
    } else {
	file_selector(_("Save data file"), SAVE_DATA, NULL);
    }	
}

/* ........................................................... */

static void get_selection_row (GtkTreeModel *model, GtkTreePath *path,
			      GtkTreeIter *iter, int *row)
{
    gint *indices = gtk_tree_path_get_indices(path);
    *row = indices[0];
}

static void count_selections (GtkTreeModel *model, GtkTreePath *path,
			      GtkTreeIter *iter, int *selcount)
{
    *selcount += 1;
}

static int selection_count (GtkTreeSelection *select, int *row)
{
    int selcount = 0;

    if (select != NULL) {
	gtk_tree_selection_selected_foreach (select, 
					     (GtkTreeSelectionForeachFunc) 
					     count_selections,
					     &selcount);
    }
    
    if (row != NULL && selcount == 1) {
	gtk_tree_selection_selected_foreach (select, 
					     (GtkTreeSelectionForeachFunc) 
					     get_selection_row,
					     row);	
    }

    return selcount;
}

int mdata_selection_count (void)
{
    return selection_count(gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox)),
			   NULL);
}

static gboolean 
main_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer data)
{
    GdkModifierType mods;

    /* ignore all but right-clicks */
    gdk_window_get_pointer(w->window, NULL, NULL, &mods);
    if (mods & GDK_BUTTON3_MASK) {
	int selcount = mdata_selection_count();

	if (mdata->popup) {
	    gtk_widget_destroy(mdata->popup);
	    mdata->popup = NULL;
	}

	if (selcount == 1) {
	    mdata->popup = build_var_popup();
	} else if (selcount > 1) {
	    mdata->popup = build_selection_popup();
	}

	if (mdata->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(mdata->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    gtk_signal_connect(GTK_OBJECT(mdata->popup), "destroy",
			       GTK_SIGNAL_FUNC(gtk_widget_destroyed), 
			       &mdata->popup);
	}

	return (selcount > 1);
    }

    return FALSE;
}
