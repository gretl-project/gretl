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
#include "version.h"
#include "ssheet.h"
#include "console.h"
#include "session.h"
#include "webget.h"
#include "database.h"
#include "datafiles.h"
#include "filelists.h"
#include "settings.h"
#include "cmdstack.h"
#include "toolbar.h"
#include "menustate.h"
#include "fileselect.h"

#include <dirent.h>
#include <unistd.h>
#include <wait.h>
#include <errno.h>

#include "../pixmaps/gretl.xpm"  /* program icon for X */

/* functions private to gretl.c */
static GtkWidget *make_main_window (int gui_get_data);

static void set_up_main_menu (void);
static void startRcallback (gpointer p, guint opt, GtkWidget *w);
static void auto_store (void);
static void sort_varlist (gpointer p, guint col, GtkWidget *w);
static void restore_sample_callback (gpointer p, int verbose, GtkWidget *w);

GdkColor red, blue, gray;

static int click_connected;

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
static int opteng, optbasque, optdump;

static const struct poptOption options[] = {
    { "run", 'r', POPT_ARG_STRING, &optrun, 0, 
      N_("open a script file on startup"), "SCRIPT" },
    { "db", 'd', POPT_ARG_STRING, &optdb, 0, 
      N_("open a database on startup"), "DATABASE" },
    { "webdb", 'w', POPT_ARG_STRING, &optdb, 0, 
      N_("open a remote (web) database on startup"), "REMOTE_DB" },
    { "english", 'e', POPT_ARG_NONE, &opteng, 0, 
      N_("force use of English"), NULL },
    { "basque", 'q', POPT_ARG_NONE, &optbasque, 0, 
      N_("force use of Basque"), NULL },
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
PATHS paths;                /* useful paths */
double **Z;                 /* data set */
MODEL **models;             /* gretl models structs */

int plot_count, data_status, orig_vars;
float gui_scale;

/* defaults for some options */
int expert = FALSE; 
int updater = FALSE;
int want_toolbar = TRUE;
int mainwin_width = 520;
int mainwin_height = 420;

char dbproxy[21];

char editor[MAXSTR] = "emacs";
char calculator[MAXSTR] = "xcalc";
char latex[MAXSTR] = "latex";
char viewdvi[MAXSTR] = "xdvi";
char viewps[MAXSTR] = "gv";
char viewpdf[MAXSTR] = "acroread";
char Rcommand[MAXSTR] = "xterm -e R";
char Browser[MAXSTR] = "mozilla";

#ifdef HAVE_TRAMO
char tramo[MAXSTR] = "tramo";
char tramodir[MAXSTR] = "";
#endif

static void spreadsheet_edit (gpointer p, guint u, GtkWidget *w) 
{
    show_spreadsheet(SHEET_EDIT_DATASET);
}

static void manual_update_query (gpointer p, guint u, GtkWidget *w)
{
    update_query();
}

#ifdef USE_GNOME
static void gnome_help (void)
{
    static GnomeHelpMenuEntry help_entry = { "gretl", "index.html" };

    gnome_help_display(NULL, &help_entry);
}
#endif /* USE_GNOME */

extern void find_var (gpointer p, guint u, GtkWidget *w); /* gui_utils.c */

static void varinfo_callback (gpointer p, guint u, GtkWidget *w)
{
    varinfo_dialog(mdata->active_var, 1);
}

#ifdef ENABLE_MAILER
static void email_data (gpointer p, guint u, GtkWidget *w)
{
    send_file(paths.datfile);
}
#endif

GtkItemFactoryEntry data_items[] = {
    /* File menu */
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },

    /* File, Open data */
    { N_("/File/_Open data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open data/user file..."), NULL, open_data, OPEN_DATA, NULL },
    { N_("/File/Open data/sample file..."), NULL, display_files, TEXTBOOK_DATA, NULL },
    { N_("/File/Open data/sep1"), NULL, NULL, 0, "<Separator>" },    
    { N_("/File/Open data/Import"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open data/Import/CSV..."), NULL, open_data, OPEN_CSV, NULL },
    { N_("/File/Open data/Import/ASCII..."), NULL, open_data, OPEN_ASCII, NULL },
    { N_("/File/Open data/Import/Octave..."), NULL, open_data, OPEN_OCTAVE, NULL },
    { N_("/File/Open data/Import/Gnumeric..."), NULL, open_data, OPEN_GNUMERIC, NULL },
    { N_("/File/Open data/Import/Excel..."), NULL, open_data, OPEN_EXCEL, NULL },
    { N_("/File/Open data/Import/Eviews..."), NULL, open_data, OPEN_WF1, NULL },
    { N_("/File/Open data/Import/Stata..."), NULL, open_data, OPEN_DTA, NULL },
    { N_("/File/Open data/Import/BOX..."), NULL, open_data, OPEN_BOX, NULL },

    /* File, Append data */
    { N_("/File/_Append data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Append data/standard format..."), NULL, open_data, APPEND_DATA, NULL },
    { N_("/File/Append data/from CSV..."), NULL, open_data, APPEND_CSV, NULL },
    { N_("/File/Append data/from ASCII..."), NULL, open_data, APPEND_ASCII, NULL },
    { N_("/File/Append data/from Octave..."), NULL, open_data, APPEND_OCTAVE, NULL },
    { N_("/File/Append data/from Gnumeric..."), NULL, open_data, APPEND_GNUMERIC, NULL },
    { N_("/File/Append data/from Excel..."), NULL, open_data, APPEND_EXCEL, NULL },
    { N_("/File/Append data/from Eviews..."), NULL, open_data, APPEND_WF1, NULL },
    { N_("/File/Append data/from Stata..."), NULL, open_data, APPEND_DTA, NULL },

    /* File, Save data */
    { N_("/File/_Save data"), NULL, auto_store, 0, NULL },
    { N_("/File/Save data _as"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Save data as/_standard format..."), NULL, file_save, SAVE_DATA_AS, NULL },
    { N_("/File/Save data as/_gzipped..."), NULL, file_save, SAVE_GZDATA, NULL },
    { N_("/File/Save data as/_database..."), NULL, file_save, SAVE_DBDATA, NULL },

    /* File, Export data */
    { N_("/File/_Export data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Export data/_CSV..."), NULL, file_save, EXPORT_CSV, NULL },
    { N_("/File/Export data/GNU _R..."), NULL, file_save, EXPORT_R, NULL },
    { N_("/File/Export data/GNU _octave..."), NULL, file_save, EXPORT_OCTAVE, NULL },
    { N_("/File/Export data/_PcGive..."), NULL, file_save, EXPORT_DAT, NULL },
#ifdef ENABLE_MAILER
    { N_("/File/Send To..."), NULL, email_data, OPEN_DATA, NULL },
#endif
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
      NULL, newdata_callback, 1, NULL },    
    { N_("/File/Create data set/time-series/quarterly"), 
      NULL, newdata_callback, 4, NULL },    
    { N_("/File/Create data set/time-series/monthly"), 
      NULL, newdata_callback, 12, NULL },
    { N_("/File/Create data set/time-series/high frequency"), 
      NULL, NULL, 0, "<Branch>" },
    { N_("/File/Create data set/time-series/high frequency/weekly"), 
      NULL, newdata_callback, 52, NULL }, 
    { N_("/File/Create data set/time-series/high frequency/daily (5-day week)"), 
      NULL, newdata_callback, 5, NULL },
    { N_("/File/Create data set/time-series/high frequency/daily (6-day week)"), 
      NULL, newdata_callback, 6, NULL },
    { N_("/File/Create data set/time-series/high frequency/daily (7-day week)"), 
      NULL, newdata_callback, 7, NULL }, 
    { N_("/File/Create data set/time-series/high frequency/hourly"), 
      NULL, newdata_callback, 24, NULL }, 
    { N_("/File/Create data set/cross-sectional"), 
      NULL, newdata_callback, 0, NULL }, 
    { N_("/File/Create data set/simulation"), NULL, gretl_callback, 
      NULLDATA, NULL },
    { N_("/File/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/File/_View command log"), NULL, view_command_log, 0, NULL },
    { N_("/File/sep2a"), NULL, NULL, 0, "<Separator>" },

    /* File, command files */
    { N_("/File/Open command file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open command file/user file..."), NULL, open_script, 
      OPEN_SCRIPT, NULL },
    { N_("/File/Open command file/practice file..."), NULL, 
      display_files, PS_FILES, NULL },
    { N_("/File/New command file"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/New command file/regular script"), NULL, do_new_script, 0, NULL },
    { N_("/File/New command file/Monte Carlo loop"), NULL, 
      do_new_script, 1, NULL },
    { N_("/File/sep3"), NULL, NULL, 0, "<Separator>" },

    /* File, preferences */
    { N_("/File/_Preferences"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/_Preferences/_General..."), NULL, options_dialog, 0, NULL },
    { N_("/File/Preferences/_Fixed font..."), NULL, font_selector, 0, NULL },
    { N_("/File/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/File/E_xit"), NULL, menu_exit_check, 0, NULL },

    /* Utilities menu */
    { N_("/_Utilities"), NULL, NULL, 0, "<Branch>" },
    { N_("/Utilities/Statistical tables"), NULL, stats_calculator, CALC_DIST, NULL },
    { N_("/Utilities/p-value finder"), NULL, stats_calculator, CALC_PVAL, NULL },
    { N_("/Utilities/Test statistic calculator"), NULL, stats_calculator, CALC_TEST, NULL },
    { N_("/Utilities/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Gretl console"), NULL, show_gretl_console, 0, NULL },
    { N_("/Utilities/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Start GNU R"), NULL, startRcallback, 0, NULL },
    { N_("/Utilities/NIST test suite"), NULL, NULL, 0, "<Branch>" },
    { N_("/Utilities/NIST test suite/basic"), NULL, do_nistcheck, 0, NULL },
    { N_("/Utilities/NIST test suite/verbose"), NULL, do_nistcheck, 1, NULL },
    { N_("/Utilities/NIST test suite/very verbose"), NULL, do_nistcheck, 2, NULL },

    /* Session menu */
    { N_("/_Session"), NULL, NULL, 0, "<Branch>" },
    { N_("/Session/_Icon view"), NULL, view_session, 0, NULL },
    { N_("/Session/sep0"), NULL, NULL, 0, "<Separator>" },
    { N_("/Session/_Open..."), NULL, open_script, OPEN_SESSION, NULL },
    { N_("/Session/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Session/_Save"), NULL, save_session_callback, 0, NULL },
    { N_("/Session/Save _as..."), NULL, save_session_callback, 1, NULL },

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
    { N_("/Data/Graph specified vars/Time series plot..."), 
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
    { N_("/Data/Graph specified vars/3D plot..."), 
      NULL, selector_callback, GR_3D, NULL },
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

    { N_("/Data/_Principal components"), NULL, do_menu_op, PCA, NULL },
    { N_("/Data/_Mahalanobis distances"), NULL, do_menu_op, MAHAL, NULL },

    { N_("/Data/sep4"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Difference of means"), NULL, NULL, 0, "<Branch>" },
    { N_("/Data/Difference of means/assuming equal variances..."), NULL, 
      selector_callback, MEANTEST, NULL },
    { N_("/Data/Difference of means/assuming unequal variances..."), NULL, 
      selector_callback, MEANTEST2, NULL },
    { N_("/Data/Difference of variances..."), NULL, selector_callback, VARTEST, NULL },
    { N_("/Data/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables"), NULL, NULL, 0, "<Branch>", },
    { N_("/Data/Add variables/time trend"), NULL, add_index, 1, NULL },
    { N_("/Data/Add variables/index variable"), NULL, add_index, 0, NULL },
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
    { N_("/Data/Add variables/unit dummies"), NULL, add_dummies, 1, NULL },
    { N_("/Data/Add variables/panel dummies"), NULL, add_dummies, 2, NULL },
    { N_("/Data/Add variables/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables/random normal..."), NULL, 
      add_random_callback, GENR_NORMAL, NULL },
    { N_("/Data/Add variables/random uniform..."), NULL, 
      add_random_callback, GENR_UNIFORM, NULL },
    { N_("/Data/Add variables/seed generator..."), NULL, gretl_callback, 
      SETSEED, NULL },
    { N_("/Data/Add variables/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables/Define _new variable..."), NULL, gretl_callback, GENR, NULL },
    { N_("/Data/Add observations..."), NULL, do_add_obs, 0, NULL },
    { N_("/Data/Remove extra observations"), NULL, do_remove_obs, 0, NULL },
    { N_("/Data/Refresh window"), NULL, refresh_data, 0, NULL },

    /* Sample menu */
    { N_("/_Sample"), NULL, NULL, 0, "<Branch>" },
    { N_("/Sample/_Set range..."), NULL, sample_range_dialog, SMPL, NULL },
    { N_("/Sample/_Restore full range"), NULL, restore_sample_callback, 1, NULL },
    { N_("/Sample/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Sample/Dataset structure..."), NULL, data_structure_wizard, 0, NULL },
    { N_("/Sample/Compact data..."), NULL, do_compact_data_set, 0, NULL },
    { N_("/Sample/Expand data..."), NULL, do_expand_data_set, 0, NULL },
    { N_("/Sample/sep2"), NULL, NULL, 0, "<Separator>" },   
    { N_("/Sample/_Define, based on dummy..."), NULL, sample_range_dialog, 
      SMPLDUM, NULL },
    { N_("/Sample/_Restrict, based on criterion..."), NULL, gretl_callback, 
      SMPLBOOL, NULL },
    { N_("/Sample/R_andom sub-sample..."), NULL, sample_range_dialog, SMPLRAND, NULL },
    { N_("/Sample/sep3"), NULL, NULL, 0, "<Separator>" },  
    { N_("/Sample/Drop all obs with _missing values"), NULL, drop_all_missing, 
      0, NULL },
    { N_("/Sample/_Count missing values"), NULL, count_missing, 0, NULL },
    { N_("/Sample/Set missing _value code..."), NULL, gretl_callback, 
      GSETMISS, NULL },
    { N_("/Sample/sep4"), NULL, NULL, 0, "<Separator>" },  
    { N_("/Sample/_Add case markers..."), NULL, open_data, OPEN_MARKERS, NULL },
    { N_("/Sample/Remove case _markers"), NULL, do_remove_markers, 0, NULL },
    { N_("/Sample/sep5"), NULL, NULL, 0, "<Separator>" },
    { N_("/Sample/Restructure panel..."), NULL, panel_restructure_dialog, 0, NULL },
    { N_("/Sample/Transpose data..."), NULL, gui_transpose_data, 0, NULL },

    /* Variable menu */
    { N_("/_Variable"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Find..."), NULL, find_var, 0, NULL },
    { N_("/Variable/_Display values"), NULL, display_var, 0, NULL },
    { N_("/Variable/_Summary statistics"), NULL, do_menu_op, 
      VAR_SUMMARY, NULL },
    { N_("/Variable/_Frequency distribution"), NULL, do_menu_op, 
      FREQ, NULL },
    { N_("/Variable/Frequency plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Frequency plot/simple"), NULL, do_freqplot, 0, NULL },
    { N_("/Variable/Frequency plot/against Normal"), NULL, do_freqplot, 
      DIST_NORMAL, NULL },
    { N_("/Variable/Frequency plot/against Gamma"), NULL, do_freqplot, 
      DIST_GAMMA, NULL },
    { N_("/Variable/Estimated density plot..."), NULL, do_kernel, 0, NULL },
    { N_("/Variable/Range-mean graph"), NULL, do_range_mean, 0, NULL }, 
    { N_("/Variable/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/_Time series plot"), NULL, ts_plot_var, 0, NULL },
    { N_("/Variable/Correlogram"), NULL, do_corrgm, CORRGM, NULL },
    { N_("/Variable/Spectrum"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Spectrum/sample periodogram"), NULL, do_pergm, 0, NULL }, 
    { N_("/Variable/Spectrum/Bartlett lag window"), NULL, do_pergm, 1, NULL }, 
    { N_("/Variable/_Augmented Dickey-Fuller test"), NULL, unit_root_test, ADF, NULL },
    { N_("/Variable/_KPSS test"), NULL, unit_root_test, KPSS, NULL },
#ifdef HAVE_X12A
    { N_("/Variable/X-12-ARIMA analysis"), NULL, do_tramo_x12a, X12A, NULL },
#endif
#ifdef HAVE_TRAMO  
    { N_("/Variable/TRAMO analysis"), NULL, do_tramo_x12a, TRAMO, NULL }, 
#endif
    { N_("/Variable/Hurst exponent"), NULL, do_hurst, 0, NULL }, 
    { N_("/Variable/Runs test"), NULL, do_menu_op, RUNS, NULL }, 
    { N_("/Variable/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Set missing value code..."), NULL, gretl_callback, 
      VSETMISS, NULL },
    { N_("/Variable/_Edit attributes"), NULL, varinfo_callback, 0, NULL },
    { N_("/Variable/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Define _new variable..."), NULL, gretl_callback, GENR, NULL },
    { N_("/_Model"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model/_Ordinary Least Squares..."), NULL, model_callback, OLS, NULL },
    { N_("/Model/_Weighted Least Squares..."), NULL, model_callback, WLS, NULL },
    { N_("/Model/sep1"),  NULL, NULL, 0, "<Separator>" },
#if 0
    { N_("/Model/HCC_M..."), NULL, model_callback, HCCM, NULL },
#endif
    { N_("/Model/H_eteroskedasticity corrected..."), NULL, model_callback, 
      HSK, NULL },
    { N_("/Model/sep2"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/Time series"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model/sep3"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Two-Stage Least Squares..."), NULL, model_callback, TSLS, NULL },
    { "/Model/sep4",  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Logit..."), NULL, model_callback, LOGIT, NULL },
    { N_("/Model/_Probit..."), NULL, model_callback, PROBIT, NULL },
    { N_("/Model/To_bit..."), NULL, model_callback, TOBIT, NULL },
    { "/Model/sep5",  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/Poi_sson..."), NULL, model_callback, POISSON, NULL },
    { N_("/Model/Lo_gistic..."), NULL, model_callback, LOGISTIC, NULL },
    { N_("/Model/Least _Absolute Deviation..."), NULL, model_callback, LAD, NULL },
    { N_("/Model/_Rank correlation..."), NULL, selector_callback, SPEARMAN, NULL },
    { N_("/Model/_Pooled OLS (panel)..."), NULL, model_callback, POOLED, NULL },
    { N_("/Model/Nonlinear Least Squares..."), NULL, gretl_callback, NLS, NULL },
    { N_("/Model/Maximum likelihood..."), NULL, gretl_callback, MLE, NULL },
#ifdef ENABLE_GMP
    { N_("/Model/High precision OLS..."), NULL, mp_ols_callback, MPOLS, NULL },
#endif
    { N_("/_Help"), NULL, NULL, 0, "<LastBranch>" },
    { N_("/Help/_GUI commands"), NULL, do_gui_help, 0, NULL },
    { N_("/Help/_Script commands syntax"), NULL, do_script_help, 1, NULL },
    { N_("/Help/sep1"), NULL, NULL, 0, "<Separator>" },
#if defined(USE_GNOME)
    { N_("/Help/Manual in HTML"), NULL, gnome_help, 0, NULL },
    { N_("/Help/sep2"), NULL, NULL, 0, "<Separator>" },
#endif
    { N_("/Help/_Check for updates"), NULL, manual_update_query, 0, NULL },
    { N_("/Help/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Help/_About gretl"), NULL, about_dialog, 0, NULL }
};

static void gui_usage (void)
{
    gui_logo(NULL);
    printf(I_("You may supply the name of a data file on the command line.\n"));
    printf(I_("Or you may do \"gretl -r script_file\" to open a script.\n"));
    printf(I_("Or you may do \"gretl -d database\" to open a gretl database.\n"));
    printf(I_("You may do \"gretl -e\" to force use of English.\n"));
    printf(I_("You may do \"gretl -q\" to force use of Basque.\n"));
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
    strncat(tryscript, fname, MAXLEN-1);

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

static void setup_column (GtkWidget *listbox, int column, int width) 
{
    if (width == 0) {
	gtk_clist_set_column_auto_resize (GTK_CLIST (listbox), column, TRUE);
    } else if (width == -1) {
	gtk_clist_set_column_visibility (GTK_CLIST (listbox), column, FALSE);
    } else {
	gtk_clist_set_column_width (GTK_CLIST (listbox), column, width);
    }
}

static void destroy (GtkWidget *widget, gpointer data)
{
    gtk_main_quit();
}

#ifdef ENABLE_NLS

void nls_init (void)
{
    char *mylang = getenv("GRETL_LANG");
    char *loc;

    if (mylang != NULL) {
	if (!g_strcasecmp(mylang, "english") ||
	    !g_strcasecmp(mylang, "C")) return;
    }

    loc = setlocale(LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
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

static void force_language (int f)
{
    if (f == ENGLISH) {
	setlocale(LC_ALL, "C");
	force_english_help();
    } else if (f == BASQUE) {
	setlocale(LC_ALL, "eu_ES");
    } 
}

#endif /* ENABLE_NLS */

int main (int argc, char *argv[])
{
    int err = 0, gui_get_data = 0;
    char filearg[MAXLEN];
    char dbname[MAXLEN];

#ifdef ENABLE_NLS
    nls_init();
#endif       

    if ((errtext = malloc(MAXLEN)) == NULL) 
	noalloc(_("startup"));

    *tryscript = '\0';
    *scriptfile = '\0';
    *paths.datfile = '\0';

    /* Initialize gnome or GTK */
#ifdef USE_GNOME
    gnome_init_with_popt_table("gretl", GRETL_VERSION, argc, argv,
			       options, 0, NULL);
#else
    gtk_init(&argc, &argv);
#endif

    libgretl_init();
    set_paths(&paths, OPT_D | OPT_X); /* defaults, gui */

    set_rcfile();
    if (!expert) {
	root_check();
    }
    init_fileptrs();

    if (argc > 1) {
	int force_lang = 0;
	int opt = parseopt((const char **) argv, argc, filearg, &force_lang);

	switch (opt) {
	case OPT_HELP:
	    gui_usage();
	    break;
	case OPT_VERSION:
	    gui_logo(NULL);
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
	if (force_lang) {
	    force_language(force_lang);
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
    if (datainfo == NULL)
	noalloc(_("data information"));

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

    /* get the data file, if specified on the command line */
    if (!(gui_get_data)) {
	int ftype;
	PRN *prn; 

	prn = gretl_print_new(GRETL_PRINT_STDERR);
	if (prn == NULL) {
	    exit(EXIT_FAILURE);
	}

	strcpy(paths.datfile, filearg);

	/* record the filename the user gave */
	strcpy(trydatfile, paths.datfile);

	ftype = detect_filetype(paths.datfile, &paths, prn);

	switch (ftype) {
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(&Z, &datainfo, paths.datfile, &paths, data_status, 
				 prn);
	    break;
	case GRETL_XML_DATA:
	    err = get_xmldata(&Z, &datainfo, paths.datfile, &paths, data_status, 
			      prn, 0);
	    break;
	case GRETL_CSV_DATA:
	    err = import_csv(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_OCTAVE:
	    err = import_octave(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_BOX_DATA:
	    err = import_box(&Z, &datainfo, paths.datfile, prn);
	    break;
	case GRETL_EXCEL:
	case GRETL_GNUMERIC:
	case GRETL_WF1:
	case GRETL_DTA:
	    err = get_worksheet_data(paths.datfile, ftype, 0, &gui_get_data);
	    break;
	case GRETL_SCRIPT:
	    gui_get_data = 1;
	    get_runfile(paths.datfile);
	    paths.datfile[0] = '\0';
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

	if (ftype != GRETL_SCRIPT && err) {
	    errmsg(err, prn);
	    exit(EXIT_FAILURE);
	}
	gretl_print_destroy(prn);
    }

    /* create the GUI */
    gretl_tooltips_init();

    /* make red, blue, gray available globally for colorizing text */
    gdk_color_parse("#ff0000", &red);
    gdk_color_parse("#0000ff", &blue);
    gdk_color_parse("#eeeeee", &gray);
    if (!gdk_colormap_alloc_color(gdk_colormap_get_system(), &red, FALSE, TRUE) ||
	!gdk_colormap_alloc_color(gdk_colormap_get_system(), &blue, FALSE, TRUE) ||
	!gdk_colormap_alloc_color(gdk_colormap_get_system(), &gray, FALSE, TRUE)) 
	noalloc(_("colors"));

    /* create main window */
    if ((mdata = mymalloc(sizeof *mdata)) == NULL)
	noalloc(_("GUI"));
    if (make_main_window(gui_get_data) == NULL) 
	noalloc(_("main window"));
    if (!gui_get_data) set_sample_label(datainfo);

    /* Let a first-time user set the working dir */
    first_time_set_user_dir(); 

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
	maybe_display_string_table();
	*trydatfile = 0;
    }

    /* opening a script from the command line? */
    if (tryscript[0] != '\0') { 
	do_open_script();
    }

    /* check for program updates? */
    proxy_init(dbproxy);
    if (updater) silent_update_query(); 

    /* try opening specified database */
    if (gui_get_data == OPT_DBOPEN)
	open_named_db_list(dbname);
    else if (gui_get_data == OPT_WEBDB)
	open_named_remote_db_list(dbname);

    /* Enter the event loop */
    gtk_main();

    /* clean up before exiting */
    free_session();

    if (Z) free_Z(Z, datainfo);

    free_model(models[0]);
    free_model(models[1]);
    free_model(models[2]);
    free(models);

    library_command_free();
    libgretl_cleanup();

    if (data_status) {
	free_datainfo(datainfo);
    }

    free_command_stack();
    exit_free_modelspec();

    return EXIT_SUCCESS;
}

/* ........................................................... */

static void augment_selection_count (gint i, gint *nsel)
{
    *nsel += 1;
}

int mdata_selection_count (void)
{
    GList *mylist = GTK_CLIST(mdata->listbox)->selection;
    gint selcount = 0;

    if (mylist != NULL) {
	g_list_foreach(mylist, (GFunc) augment_selection_count, 
		       &selcount);
    }

    return selcount;
}

/* ........................................................... */

gint main_varclick (GtkWidget *widget, GdkEventButton *event, 
		    gpointer data)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    if (mdata->active_var == 0) {
	return FALSE;
    }

    topwin = gtk_widget_get_parent_window(mdata->listbox);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON3_MASK) { 
	gint selcount = mdata_selection_count();

	if (mdata->popup) {
	    gtk_widget_destroy(mdata->popup);
	    mdata->popup = NULL;
	}	

	if (selcount == 1) {
	    mdata->popup = build_var_popup();
	} else if (selcount > 1) {
	    mdata->popup = build_selection_popup();
	}

	gtk_menu_popup(GTK_MENU(mdata->popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
	gtk_signal_connect(GTK_OBJECT(mdata->popup), "destroy",
			   GTK_SIGNAL_FUNC(gtk_widget_destroyed), 
			   &mdata->popup);
    }

    return TRUE;
}

/* ........................................................... */

static void check_varmenu_state (GtkCList *list, gint i, gint j,
				 GdkEventButton *event, gpointer p)
{
    if (mdata->ifac != NULL) {
	gint selcount = mdata_selection_count();

	flip(mdata->ifac, "/Variable", (selcount == 1));
	flip(mdata->ifac, "/Data/Correlation matrix/selected variables", 
	     (selcount > 1));
	flip(mdata->ifac, "/Data/Principal components", 
	     (selcount > 1));
	flip(mdata->ifac, "/Data/Mahalanobis distances", 
	     (selcount > 1));
    }
}

/* ........................................................... */

static gint catch_mdata_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
#ifdef HAVE_FLITE
    if (key->keyval == GDK_a) {
	audio_render_window(vwin, AUDIO_LISTBOX);
    } else if (key->keyval == GDK_x) {
	audio_render_window(NULL, AUDIO_LISTBOX);
    }
#endif

    if (key->keyval == GDK_e) {
	varinfo_dialog(mdata->active_var, 1);
    } 
    else if (key->keyval == GDK_t) {
	do_graph_var(mdata->active_var);
    } 

    return 1;
}

/* ........................................................... */

void populate_varlist (void)
{
    char id[4];
    char *row[3];
    gint i, j;
    static gint check_connected;

    gtk_clist_freeze(GTK_CLIST(mdata->listbox));

    gtk_clist_clear(GTK_CLIST(mdata->listbox));

    j = 0;
    for (i=0; i<datainfo->v; i++) {
	if (is_hidden_variable(i, datainfo) ||
	    is_standard_lag(i, datainfo)) {
	    continue;
	}
	sprintf(id, "%d", i);
	row[0] = id;
	row[1] = datainfo->varname[i];
	row[2] = VARLABEL(datainfo, i);
	gtk_clist_append(GTK_CLIST(mdata->listbox), row);
	if (j % 2) {
	    gtk_clist_set_background(GTK_CLIST(mdata->listbox), j, &gray);
	}
	j++;
    }

    if (j >= 2) {
	mdata->active_var = 1;	
	/* both select and focus the first real variable */
	gtk_clist_select_row(GTK_CLIST(mdata->listbox), mdata->active_var, 1); 
	(GTK_CLIST(mdata->listbox))->focus_row = mdata->active_var;
    } else {
	mdata->active_var = 0;
    }

    gtk_clist_set_selectable(GTK_CLIST(mdata->listbox), 
			     0, FALSE);

    gtk_clist_thaw(GTK_CLIST(mdata->listbox));

    if (!check_connected) {
	gtk_signal_connect(GTK_OBJECT(GTK_CLIST(mdata->listbox)),
			   "select-row",
			   check_varmenu_state, NULL);
	gtk_signal_connect(GTK_OBJECT(GTK_CLIST(mdata->listbox)),
			   "unselect-row",
			   check_varmenu_state, NULL);
	check_connected = 1;
    }


    if (!click_connected) {
	gtk_signal_connect(GTK_OBJECT(mdata->listbox),
			   "button_press_event",
			   GTK_SIGNAL_FUNC(main_varclick), NULL);
	gtk_signal_connect (GTK_OBJECT(mdata->listbox), "key_press_event",
			    GTK_SIGNAL_FUNC(catch_mdata_key),
			    mdata);
	click_connected = 1;
    }
}

/* ........................................................... */

void clear_varlist (GtkWidget *widget)
{
    gtk_clist_clear(GTK_CLIST(widget));
    if (click_connected) {
	gtk_signal_disconnect_by_func(GTK_OBJECT(mdata->listbox),
				      GTK_SIGNAL_FUNC(main_varclick), 
				      NULL);
	click_connected = 0;
    }
}

/* ......................................................... */

static float get_gui_scale (void)
{
    GtkStyle *style = NULL;
    int fsize = 0;
    float scale = 1.0;

    style = gtk_rc_get_style(mdata->w);

    if (style != NULL && style->font != NULL) {
	fsize = gdk_char_width(style->font, 'X');
    }

    if (fsize > 0) scale = fsize / 9.0;
    
    return scale;
}

/* .................................................................. */

static GtkWidget *list_box_create (GtkBox *box, char *titles[]) 
{
    GtkWidget *view, *scroller;
    int listbox_id_width = 30;
    int listbox_varname_width = 90;

    if (strcmp(titles[1], "Variable name")) {
	listbox_varname_width = 110;
    }

    view = gtk_clist_new_with_titles (3, titles);
    gtk_clist_column_titles_passive(GTK_CLIST(view));
    gtk_clist_set_selection_mode (GTK_CLIST (view), 
				  GTK_SELECTION_EXTENDED); /* was BROWSE */
    gtk_clist_set_column_width (GTK_CLIST (view), 
				0, listbox_id_width * gui_scale);
    gtk_clist_set_column_justification (GTK_CLIST (view), 0, 
					GTK_JUSTIFY_LEFT);

    setup_column(view, 1, listbox_varname_width * gui_scale);
    setup_column(view, 2, 0); 

    gtk_signal_connect_after(GTK_OBJECT(view), "select_row",
			     GTK_SIGNAL_FUNC(selectrow), (gpointer) mdata);
    gtk_signal_connect_after(GTK_OBJECT(view), "unselect_row",
			     GTK_SIGNAL_FUNC(unselectrow), (gpointer) mdata);

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(box, scroller, TRUE, TRUE, TRUE);

    gtk_widget_show(view);
    gtk_widget_show(scroller);

    return view;
}

/* ......................................................... */

static GtkWidget *make_main_window (int gui_get_data) 
{
    GtkWidget *main_vbox;
    GtkWidget *box, *dlabel, *align;
    char *titles[] = {
	_("ID #"), 
	_("Variable name"), 
	_("Descriptive label")
    };

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

    gui_scale = get_gui_scale();

    if (mainwin_width <= 200 || mainwin_height <= 200) {
	mainwin_width = 580 * gui_scale;
	mainwin_height = 420 * gui_scale;
    }

    gtk_signal_connect (GTK_OBJECT (mdata->w), "delete_event",
			GTK_SIGNAL_FUNC (exit_check), NULL);
    gtk_signal_connect (GTK_OBJECT (mdata->w), "destroy",
			GTK_SIGNAL_FUNC (destroy), NULL);

    gtk_window_set_title(GTK_WINDOW (mdata->w), "gretl");
    gtk_window_set_default_size(GTK_WINDOW (mdata->w), 
				mainwin_width, mainwin_height);
    gtk_signal_connect_after(GTK_OBJECT(mdata->w), "realize", 
                             GTK_SIGNAL_FUNC(set_wm_icon), 
                             NULL);

    main_vbox = gtk_vbox_new(FALSE, 5);
    gtk_container_set_border_width(GTK_CONTAINER (main_vbox), 10);

#ifdef USE_GNOME
    gnome_app_set_contents(GNOME_APP(mdata->w), main_vbox);
#else
    gtk_container_add(GTK_CONTAINER(mdata->w), main_vbox);
#endif

    gtk_object_set_data(GTK_OBJECT(mdata->w), "vbox", main_vbox);

    set_up_main_menu();
    gtk_box_pack_start(GTK_BOX(main_vbox), mdata->mbar, FALSE, TRUE, 0);
    gtk_widget_show(mdata->mbar);

    dlabel = gtk_label_new(_(" No datafile loaded ")); 
    gtk_widget_show(dlabel);

    gtk_object_set_data(GTK_OBJECT(mdata->w), "dlabel", dlabel);

    box = gtk_vbox_new (FALSE, 0);
    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(box), align, FALSE, FALSE, 0);
    gtk_widget_show(align);
    gtk_container_add(GTK_CONTAINER(align), dlabel);
   
    mdata->listbox = list_box_create (GTK_BOX(box), titles);
    gtk_widget_show(box);

    gtk_drag_dest_set (mdata->listbox,
		       GTK_DEST_DEFAULT_ALL,
		       gretl_drag_targets, 2,
		       GDK_ACTION_COPY);

    gtk_signal_connect (GTK_OBJECT(mdata->listbox), "drag_data_received",
			GTK_SIGNAL_FUNC(drag_data_received),
			NULL);

    gtk_box_pack_start(GTK_BOX (main_vbox), box, TRUE, TRUE, 0);

    mdata->status = gtk_label_new("");
    
    gtk_box_pack_start (GTK_BOX (main_vbox), mdata->status, FALSE, TRUE, 0);

    /* put stuff into list box, activate menus */
    if (!gui_get_data) populate_varlist();

    /* create gretl toolbar? */
    show_or_hide_toolbar(want_toolbar);

    /* get a monospaced font for various windows */
    set_fixed_font();

    gtk_widget_show_all(mdata->w); 

    return main_vbox;
}

/* ........................................................... */

static void set_up_main_menu (void)
{
    GtkAccelGroup *main_accel;
   
    gint n_items = sizeof data_items / sizeof data_items[0];

    main_accel = gtk_accel_group_new();
    mdata->ifac = gtk_item_factory_new (GTK_TYPE_MENU_BAR, "<main>", 
					main_accel);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(mdata->ifac, menu_translate, NULL, NULL);
#endif    
    gtk_item_factory_create_items (mdata->ifac, n_items, data_items, NULL);
    mdata->mbar = gtk_item_factory_get_widget (mdata->ifac, "<main>");
    gtk_accel_group_attach(main_accel, GTK_OBJECT (mdata->w));
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
    int err = restore_sample(OPT_NONE); /* FIXME OPT_C? */

    if (verbose && !err) {
	set_sample_label(datainfo);    
	gretl_command_strcpy("smpl full");
	check_and_record_command();
    }
}

static void startRcallback (gpointer p, guint opt, GtkWidget *w)
{
    startR(Rcommand);
}

static volatile int fork_err;

static void fork_err_set (int signum)
{
    fork_err = 1;
}

int gretl_fork (const char *prog, const char *arg)
{
    pid_t pid;

    fork_err = 0;

    signal(SIGUSR1, fork_err_set);

    pid = fork();
    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return 1;
    } else if (pid == 0) {
	if (arg != NULL) { 
	    execlp(prog, prog, arg, NULL);
	} else {
	    execlp(prog, prog, NULL);
	}
	perror("execlp");
	kill(getppid(), SIGUSR1);
	_exit(EXIT_FAILURE);
    }
    
    sleep(1);

    if (fork_err) {
	sprintf(errtext, "%s: %s", _("Command failed"), prog);
	errbox(errtext);
    }

    return fork_err;
}

/* Icon handling for X */
void set_wm_icon (GtkWidget *w, gpointer data)
{
    GdkPixmap *icon;

    icon = gdk_pixmap_create_from_xpm_d(w->window, NULL, NULL, gretl_xpm);
    gdk_window_set_icon(w->window, NULL, icon, NULL);
}

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
    gchar *dfname = NULL;
    char *suff = NULL, tmp[MAXLEN];
    int pos, skip = 5;

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

    /* there may be multiple files: we ignore all but the first */
    *tmp = 0;
    if ((pos = haschar('\r', dfname)) > 0 || 
	(pos = haschar('\n', dfname) > 0)) {
	strncat(tmp, dfname + skip, pos - skip);
    } else {
	strcat(tmp, dfname + skip);
    }

    unescape_url(tmp);

    suff = strrchr(tmp, '.');
    if (suff && (!strncmp(suff, ".gretl", 6) || 
		 !strncmp(suff, ".inp", 4) ||
		 !strncmp(suff, ".GRE", 4) ||
		 !strncmp(suff, ".INP", 4))) {
	strcpy(tryscript, tmp);
	fprintf(stderr, "tryscript='%s'\n", tryscript);
	verify_open_session(NULL);
    } else {
	strcpy(trydatfile, tmp);
	verify_open_data(NULL, 0);
    }	
}

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
    if (strlen(paths.datfile) && is_gzipped(paths.datfile)) {
	oflag = OPT_Z;
    }

    if ((data_status & USER_DATA) && native_datafile()) {
	do_store(paths.datfile, oflag, 1);
    } else {
	file_selector(_("Save data file"), SAVE_DATA, FSEL_DATA_NONE, NULL);
    }	
}

/* ........................................................... */

static gint 
compare_var_ids (GtkCList *cl, gconstpointer p1, gconstpointer p2)
{
    const GtkCListRow *row1 = (const GtkCListRow *) p1;
    const GtkCListRow *row2 = (const GtkCListRow *) p2;
    gchar *t1, *t2;
    int i1, i2;

    t1 = GTK_CELL_TEXT(row1->cell[cl->sort_column])->text;
    t2 = GTK_CELL_TEXT(row2->cell[cl->sort_column])->text;    

    i1 = atoi(t1);
    i2 = atoi(t2);

    return ((i1 < i2) ? -1 : (i1 > i2) ? 1 : 0);
}

static gint 
compare_varnames (GtkCList *cl, gconstpointer p1, gconstpointer p2)
{
    const GtkCListRow *row1 = (const GtkCListRow *) p1;
    const GtkCListRow *row2 = (const GtkCListRow *) p2;
    gchar *t1, *t2;

    t1 = GTK_CELL_TEXT(row1->cell[cl->sort_column])->text;
    t2 = GTK_CELL_TEXT(row2->cell[cl->sort_column])->text;    

    if (!strcmp(t1, "const")) return 0;
    if (!strcmp(t2, "const")) return 1;

    return strcmp(t1, t2);
}

static void sort_varlist (gpointer p, guint col, GtkWidget *w)
{
    if (col == 0) {
	gtk_clist_set_compare_func(GTK_CLIST(mdata->listbox), compare_var_ids);
    } else {
	gtk_clist_set_compare_func(GTK_CLIST(mdata->listbox), compare_varnames);
    }
    gtk_clist_set_sort_column(GTK_CLIST(mdata->listbox), col);
    gtk_clist_sort(GTK_CLIST(mdata->listbox));
}
