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
#include "guiprint.h"
#include "console.h"
#include "session.h"

#include <dirent.h>
#include <unistd.h>
#include <signal.h>

#include "../pixmaps/gretl.xpm"  /* program icon for X */

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
extern void allocate_fileptrs (void); /* settings.c */
extern void stats_calculator (gpointer data, guint ttest, GtkWidget *widget);
extern void bool_subsample (gpointer data, guint dropmiss, GtkWidget *w);
extern void free_command_stack (void);
extern void open_named_db_clist (char *dbname);
extern void open_named_remote_clist (char *dbname);
extern void gui_set_panel_structure (gpointer data, guint u, GtkWidget *w);
extern void time_series_dialog (gpointer data, guint u, GtkWidget *w);
extern void panel_restructure_dialog (gpointer data, guint u, GtkWidget *w);

/* functions private to gretl.c */
static void make_toolbar (GtkWidget *w, GtkWidget *box);
static void clip_init (GtkWidget *w);
static GtkWidget *make_main_window (int gui_get_data);
static GtkWidget *build_var_popup (void);
static void build_selection_popup (void);
static gint popup_activated (GtkWidget *widget, gpointer data);
static void check_for_extra_data (void);
static void set_up_main_menu (void);
static void startR (gpointer p, guint opt, GtkWidget *w);
static void auto_store (void);
static void sort_varlist (gpointer p, guint col, GtkWidget *w);
static void restore_sample_callback (gpointer p, int verbose, GtkWidget *w);

GtkWidget *toolbar_box = NULL; /* shared with settings.c */

static GtkWidget *datalabel;
static GtkWidget *main_vbox;
static GtkWidget *gretl_toolbar;
static GtkWidget *selection_popup;

GdkColor red, blue, gray;

static int popup_connected;
int *default_list = NULL;

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
char *optrun = NULL, *optdb = NULL;

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

int plot_count, data_status, orig_vars;
gchar *clipboard_buf; /* for copying models as LaTeX */
float gui_scale;

/* defaults for some options */
int expert = FALSE; 
int updater = FALSE;
int want_toolbar = TRUE;
char dbproxy[21];

char editor[MAXSTR] = "emacs";
char calculator[MAXSTR] = "xcalc";
char viewdvi[MAXSTR] = "xdvi";

#ifdef USE_GNOME
char Rcommand[MAXSTR] = "R --gui=gnome";
extern const char *version_string;
#else
char Rcommand[MAXSTR] = "xterm -e R";
#endif

#ifdef HAVE_TRAMO
char tramo[MAXSTR] = "tramo";
char tramodir[MAXSTR] = "";
#endif

#ifdef HAVE_X12A
char x12a[MAXSTR] = "x12a";
char x12adir[MAXSTR] = "";
#endif

static void spreadsheet_edit (gpointer p, guint u, GtkWidget *w) 
{
    extern void show_spreadsheet (DATAINFO *pdinfo);  
  
    show_spreadsheet(NULL);
}

static void manual_update_query (gpointer p, guint u, GtkWidget *w)
{
    update_query(1);
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
    varinfo_dialog(mdata->active_var);
}

GtkItemFactoryEntry data_items[] = {
    /* File menu */
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },
    /* File, Open data */
    { N_("/File/_Open data"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/Open data/user file..."), NULL, open_data, OPEN_DATA, NULL },
    { N_("/File/Open data/sample file..."), NULL, display_files, TEXTBOOK_DATA, NULL },
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
      SAVE_DATA_AS, NULL },
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
    { N_("/Utilities/Statistical tables"), NULL, stats_calculator, 1, NULL },
    { N_("/Utilities/p-value finder"), NULL, stats_calculator, 0, NULL },
    { N_("/Utilities/Test statistic calculator"), NULL, stats_calculator, 2, NULL },
    { N_("/Utilities/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Gretl console"), NULL, show_gretl_console, 0, NULL },
    { N_("/Utilities/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Utilities/Start GNU R"), NULL, startR, 0, NULL },

    /* Session menu */
    { N_("/_Session"), NULL, NULL, 0, "<Branch>" },
    { N_("/Session/_Icon view"), NULL, view_session, 0, NULL },
#ifndef GNUPLOT_PNG
    { N_("/Session/_Add last graph"), NULL, add_graph_to_session, 
      GRETL_GNUPLOT_GRAPH, NULL },
#endif
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
    { N_("/Data/_Multiple scatterplots..."), 
      NULL, selector_callback, SCATTERS, NULL},
    { N_("/Data/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/_Read info"), NULL, open_info, 0, NULL },
    { N_("/Data/Edit _info"), NULL, edit_header, 0, NULL },
    { N_("/Data/Print description"), NULL, print_report, 0, NULL },
    { N_("/Data/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/_Summary statistics"), NULL, do_menu_op, SUMMARY, NULL },
    { N_("/Data/_Correlation matrix"), NULL, do_menu_op, CORR, NULL },
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
    { N_("/Data/Add variables/logs of variables..."), NULL, 
      addvars_dialog, LOGS, NULL },
    { N_("/Data/Add variables/lags of variables..."), NULL, 
      addvars_dialog, LAGS, NULL },
    { N_("/Data/Add variables/squares of variables..."), NULL, 
      addvars_dialog, SQUARE, NULL },
    { N_("/Data/Add variables/periodic dummies"), NULL, add_dummies, 0, NULL },
    { N_("/Data/Add variables/panel dummies"), NULL, add_dummies, 1, NULL },
    { N_("/Data/Add variables/first differences..."), NULL, 
      addvars_dialog, DIFF, NULL },
    { N_("/Data/Add variables/log differences..."), NULL, 
      addvars_dialog, LDIFF, NULL },
    { N_("/Data/Add variables/sep"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables/random normal..."), NULL, 
      random_dialog, GENR_NORMAL, NULL },
    { N_("/Data/Add variables/random uniform..."), NULL, 
      random_dialog, GENR_UNIFORM, NULL },
    { N_("/Data/Add variables/seed generator..."), NULL, gretl_callback, 
      SEED, NULL },
    { N_("/Data/Add variables/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Data/Add variables/Define _new variable..."), NULL, gretl_callback, GENR, NULL },
    { N_("/Data/Refresh window"), NULL, refresh_data, 0, NULL },

    /* Sample menu */
    { N_("/_Sample"), NULL, NULL, 0, "<Branch>" },
    { N_("/Sample/_Set range..."), NULL, gretl_callback, SMPL, NULL },
    { N_("/Sample/_Restore full range"), NULL, restore_sample_callback, 1, NULL },
    { N_("/Sample/sep1"), NULL, NULL, 0, "<Separator>" },    
    { N_("/Sample/Set _frequency, startobs..."), NULL, gretl_callback, 
      SETOBS, NULL },
    { N_("/Sample/Compact data..."), NULL, compact_data_set, 0, NULL },
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
    { N_("/Sample/Restructure panel..."), NULL, panel_restructure_dialog, 0, NULL },

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
      NORMAL, NULL },
    { N_("/Variable/Frequency plot/against Gamma"), NULL, do_freqplot, 
      GAMMA, NULL },
    { N_("/Variable/Range-mean graph"), NULL, do_range_mean, 0, NULL }, 
    { N_("/Variable/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/_Time series plot"), NULL, ts_plot_var, 0, NULL },
    { N_("/Variable/Correlogram"), NULL, gretl_callback, CORRGM, NULL },
    { N_("/Variable/Spectrum"), NULL, NULL, 0, "<Branch>" },
    { N_("/Variable/Spectrum/sample periodogram"), NULL, do_pergm, 0, NULL }, 
    { N_("/Variable/Spectrum/Bartlett lag window"), NULL, do_pergm, 1, NULL }, 
    { N_("/Variable/_Augmented Dickey-Fuller test"), NULL, gretl_callback, 
      ADF, NULL },
#ifdef HAVE_X12A
    { N_("/Variable/X-12-ARIMA analysis"), NULL, do_tramo_x12a, X12A, NULL },
#endif
#ifdef HAVE_TRAMO  
    { N_("/Variable/TRAMO analysis"), NULL, do_tramo_x12a, TRAMO, NULL }, 
#endif
    { N_("/Variable/Runs test"), NULL, do_menu_op, RUNS, NULL }, 
    { N_("/Variable/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Set missing value code..."), NULL, gretl_callback, 
      VSETMISS, NULL },
    { N_("/Variable/_Edit attributes"), NULL, varinfo_callback, 0, NULL },
    { N_("/Variable/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Variable/Simulate..."), NULL, gretl_callback, SIM, NULL },
    { N_("/Variable/Define _new variable..."), NULL, gretl_callback, GENR, NULL },
    { N_("/Variable/Delete last variable"), NULL, delete_var, 0, NULL },
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
    { N_("/Model/_Vector Autoregression..."), NULL, model_callback, VAR, NULL },
    { N_("/Model/Cointegration test"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model/Cointegration test/Engle-Granger..."), NULL, 
      selector_callback, COINT, NULL },
    { N_("/Model/Cointegration test/Johansen..."), NULL, 
      selector_callback, COINT2, NULL },
    { N_("/Model/sep3"),  NULL, NULL, 0, "<Separator>" },
    { N_("/Model/_Two-Stage Least Squares..."), NULL, model_callback, TSLS, NULL },
    { N_("/Model/_Logit..."), NULL, model_callback, LOGIT, NULL },
    { N_("/Model/_Probit..."), NULL, model_callback, PROBIT, NULL },
    { N_("/Model/Least _Absolute Deviation..."), NULL, model_callback, LAD, NULL },
    { N_("/Model/_Rank correlation..."), NULL, gretl_callback, SPEARMAN, NULL },
    { N_("/Model/_Pooled OLS (panel)..."), NULL, model_callback, POOLED, NULL },
    { N_("/Model/Nonlinear Least Squares..."), NULL, gretl_callback, NLS, NULL },
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

    tryscript[0] = '\0';
    strncat(tryscript, fname, MAXLEN-1);

    if (addpath(tryscript, &paths, 1) == NULL) {
	fprintf(stderr, _("Couldn't find script '%s'\n"), tryscript);
	exit(EXIT_FAILURE);
    } else {
	fprintf(stderr, _("%s found\n"), tryscript);
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

#ifdef ENABLE_NLS
void nls_init (void)
{
    setlocale (LC_ALL, "");
    bindtextdomain (PACKAGE, LOCALEDIR);
    textdomain (PACKAGE);
}
#endif /* ENABLE_NLS */

int main (int argc, char *argv[])
{
    int err = 0, gui_get_data = 0;
    char dbname[MAXLEN];

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
    gnome_init_with_popt_table("gretl", version_string, argc, argv,
			       options, 0, NULL);
#else
    gtk_init(&argc, &argv);
#endif

    set_paths(&paths, 1, 1); /* 1 = defaults, 1 = gui */
    set_rcfile();
    make_userdir(&paths);

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

    gretl_rand_init();
    helpfile_init();
    session_init();

    /* get the data file, if specified on the command line */
    if (!(gui_get_data)) {
	int ftype;
	PRN *prn; 

	prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
	if (prn == NULL) 
	    exit(EXIT_FAILURE);

	*paths.datfile = '\0';
	strncat(paths.datfile, argv[1], MAXLEN-1);

	/* record the filename the user gave */
	strcpy(trydatfile, paths.datfile);

	ftype = detect_filetype(paths.datfile, &paths, prn);

	switch (ftype) {
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
	case GRETL_EXCEL:
	case GRETL_GNUMERIC:
	    err = get_worksheet_data(paths.datfile, ftype, 0);
	    break;
	case GRETL_SCRIPT:
	    gui_get_data = 1;
	    get_runfile(paths.datfile);
	    paths.datfile[0] = '\0';
	    break;
	case GRETL_UNRECOGNIZED:
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:    
	    exit(EXIT_FAILURE);
	    break;
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
    if ((mdata = mymalloc(sizeof(windata_t))) == NULL)
	noalloc(_("GUI"));
    if ((datalabel = make_main_window(gui_get_data)) == NULL) 
	noalloc(_("main window"));
    if (!gui_get_data) set_sample_label(datainfo);

    /* enable special copying to clipboard */
    clip_init(mdata->w);

    allocate_fileptrs();
    add_files_to_menu(FILE_LIST_DATA);
    add_files_to_menu(FILE_LIST_SESSION);
    add_files_to_menu(FILE_LIST_SCRIPT);
#ifndef GNUPLOT_PNG
    graphmenu_state(FALSE);
#endif
    session_menu_state(FALSE);
    restore_sample_state(FALSE);
    menubar_state(FALSE);
			  
    check_for_extra_data();
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
	do_open_script(NULL, NULL);
    }

    /* check for program updates? */
    proxy_init(dbproxy);
    if (updater) update_query(0); 

    /* try opening specified database */
    if (gui_get_data == OPT_DBOPEN)
	open_named_db_clist(dbname);
    else if (gui_get_data == OPT_WEBDB)
	open_named_remote_clist(dbname);

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

    return EXIT_SUCCESS;
}

/* ........................................................... */

void refresh_data (void)
{
    if (data_status) {
	populate_main_varlist();
    }
}

/* ........................................................... */

void menubar_state (gboolean s)
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

#ifndef GNUPLOT_PNG
void graphmenu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/File/Save last graph", s);
	flip(mdata->ifac, "/Session/Add last graph", s);
    }
}
#endif

/* ........................................................... */

static void time_series_menu_state (gboolean s)
{
    if (mdata->ifac != NULL) {
	flip(mdata->ifac, "/Data/Graph specified vars/Time series plot...", s);
	flip(mdata->ifac, "/Variable/Time series plot", s);
	flip(mdata->ifac, "/Variable/Correlogram", s);
	flip(mdata->ifac, "/Variable/Spectrum", s);
	flip(mdata->ifac, "/Variable/Augmented Dickey-Fuller test", s);
#ifdef HAVE_X12A
	flip(mdata->ifac, "/Variable/X-12-ARIMA analysis", s);
#endif
#ifdef HAVE_TRAMO
	flip(mdata->ifac, "/Variable/TRAMO analysis", s);
#endif
	flip(mdata->ifac, "/Variable/Runs test", s);
	flip(mdata->ifac, "/Model/Cochrane-Orcutt...", s);
	flip(mdata->ifac, "/Model/Hildreth-Lu...", s);
	flip(mdata->ifac, "/Model/Autoregressive estimation...", s);
	flip(mdata->ifac, "/Model/Vector Autoregression...", s);
	flip(mdata->ifac, "/Model/Cointegration test", s);
	flip(mdata->ifac, "/Sample/Compact data...", 
	     s && (datainfo->pd == 4 || datainfo->pd == 12));
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

static void augment_selection_count (gint i, gint *nsel)
{
    *nsel += 1;
}

static gint get_mdata_selection (void)
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

gint main_popup (GtkWidget *widget, GdkEventButton *event, 
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
	gint selcount = get_mdata_selection();

	if (selcount == 1) {
	    if (mdata->popup) {
		gtk_widget_destroy(mdata->popup);
	    }
	    mdata->popup = build_var_popup();
	    gtk_menu_popup(GTK_MENU(mdata->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	} else if (selcount > 1) {
	    if (selection_popup) {
		gtk_widget_destroy(selection_popup);
	    }
	    build_selection_popup();
	    gtk_menu_popup(GTK_MENU(selection_popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	}
    }
    return TRUE;
}

/* ........................................................... */

void check_varmenu_state (GtkCList *list, gint i, gint j,
			  GdkEventButton *event, gpointer p)
{
    if (mdata->ifac != NULL) {
	gint selcount = get_mdata_selection();

	flip(mdata->ifac, "/Variable", (selcount == 1));
    }
}

/* ........................................................... */

gint populate_main_varlist (void)
{
    char id[4];
    char *row[3];
    gint i;
    static gint check_connected;

    gtk_clist_clear(GTK_CLIST(mdata->listbox));

    for (i=0; i<datainfo->v; i++) {
	if (hidden_var(i, datainfo)) continue;
	sprintf(id, "%d", i);
	row[0] = id;
	row[1] = datainfo->varname[i];
	row[2] = VARLABEL(datainfo, i);
	gtk_clist_append(GTK_CLIST(mdata->listbox), row);
	if (i % 2) {
	    gtk_clist_set_background(GTK_CLIST(mdata->listbox), i, &gray);
	}
    }

    mdata->active_var = 1;
    if (mdata->active_var > datainfo->v - 1) {
	mdata->active_var -= 1;
    }

    gtk_clist_set_selectable(GTK_CLIST(mdata->listbox), 
			     0, FALSE);

    gtk_clist_select_row 
	(GTK_CLIST(mdata->listbox), mdata->active_var, 1);  

    if (!check_connected) {
	gtk_signal_connect(GTK_OBJECT(GTK_CLIST(mdata->listbox)),
			   "select-row",
			   check_varmenu_state, NULL);
	check_connected = 1;
    }

    if (!popup_connected) {
	gtk_signal_connect(GTK_OBJECT(mdata->listbox),
			   "button_press_event",
			   GTK_SIGNAL_FUNC(main_popup), NULL);
	popup_connected = 1;
    }

    return 0;
}

/* ........................................................... */

void clear_clist (GtkWidget *widget)
{
    gtk_clist_clear(GTK_CLIST(widget));
    if (popup_connected) {
	gtk_signal_disconnect_by_func(GTK_OBJECT(mdata->listbox),
				      (GtkSignalFunc) main_popup, 
				      NULL);
	popup_connected = 0;
    }
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
    time_series_menu_state(dataset_is_time_series(pdinfo));

    flip(mdata->ifac, "/Sample/Interpret as time series...", 
	 !(dataset_is_time_series(pdinfo)));

    flip(mdata->ifac, "/Sample/Interpret as panel...", 
	 !(pdinfo->pd == 1));

    flip(mdata->ifac, "/Sample/Restructure panel...", 
	 pdinfo->time_series == STACKED_CROSS_SECTION);

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
    int listbox_label_width = 400;

    if (strcmp(titles[1], "Variable name")) listbox_varname_width = 110;

    view = gtk_clist_new_with_titles (3, titles);
    gtk_clist_column_titles_passive(GTK_CLIST(view));
    gtk_clist_set_selection_mode (GTK_CLIST (view), 
				  GTK_SELECTION_EXTENDED); /* was BROWSE */
    gtk_clist_set_column_width (GTK_CLIST (view), 
				0, listbox_id_width * gui_scale);
    gtk_clist_set_column_justification (GTK_CLIST (view), 0, 
					GTK_JUSTIFY_LEFT);
    setup_column(view, 1, listbox_varname_width * gui_scale);
    setup_column(view, 2, listbox_label_width * gui_scale);

    gtk_signal_connect_after (GTK_OBJECT (view), "select_row",
			      GTK_SIGNAL_FUNC (selectrow), (gpointer) mdata);

    scroller = gtk_scrolled_window_new (NULL, NULL);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (scroller),
				    GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_container_add (GTK_CONTAINER(scroller), view);

    gtk_box_pack_start (box, scroller, TRUE, TRUE, TRUE);

    gtk_widget_show(view);
    gtk_widget_show(scroller);

    return view;
}

/* ......................................................... */

static GtkWidget *make_main_window (int gui_get_data) 
{
    GtkWidget *box, *dlabel, *align;
    char *titles[] = {
	_("ID #"), 
	_("Variable name"), 
	_("Descriptive label")
    };
    int mainwin_width = 520;
    int mainwin_height = 400;

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

    gtk_signal_connect (GTK_OBJECT (mdata->w), "delete_event",
			GTK_SIGNAL_FUNC (exit_check), NULL);
    gtk_signal_connect (GTK_OBJECT (mdata->w), "destroy",
			GTK_SIGNAL_FUNC (destroy), NULL);

    gui_scale = get_gui_scale();

    mainwin_width *= gui_scale;
    mainwin_height *= gui_scale;

    gtk_window_set_title(GTK_WINDOW (mdata->w), "gretl");
    gtk_window_set_default_size(GTK_WINDOW (mdata->w), 
				mainwin_width, mainwin_height);
    gtk_signal_connect_after(GTK_OBJECT(mdata->w), "realize", 
                             GTK_SIGNAL_FUNC(set_wm_icon), 
                             NULL);

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
    if (!gui_get_data) populate_main_varlist();

    /* create gretl toolbar */
    if (want_toolbar) make_toolbar(mdata->w, main_vbox);

    /* get a monospaced font for various windows */
    load_fixed_font();

    /* build popup menu for multiple selections */
    build_selection_popup();

    gtk_widget_show_all(mdata->w); 

    return dlabel;
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

extern void delete_var_by_id (int id); /* callbacks.c */

static gint popup_activated (GtkWidget *widget, gpointer data)
{
    gchar *item;

    item = (gchar *) data;
    if (!strcmp(item, _("Display values"))) 
	display_var();
    if (!strcmp(item, _("Descriptive statistics"))) 
	do_menu_op(NULL, VAR_SUMMARY, NULL);
    else if (!strcmp(item, _("Time series plot"))) 
	do_graph_var(mdata->active_var);
    else if (!strcmp(item, _("Frequency distribution"))) 
	do_menu_op(NULL, FREQ, NULL);
    else if (!strcmp(item, _("Frequency plot"))) 
	do_freqplot(NULL, 0, NULL);
    else if (!strcmp(item, _("Boxplot")))
	do_boxplot_var(mdata->active_var);
    else if (!strcmp(item, _("Correlogram"))) 
	gretl_callback(NULL, CORRGM, NULL);
    else if (!strcmp(item, _("Spectrum"))) 
	do_pergm(NULL, 0, NULL);
    else if (!strcmp(item, _("Dickey-Fuller test"))) 
	gretl_callback(NULL, ADF, NULL);
    else if (!strcmp(item, _("Runs test"))) 
	do_menu_op(NULL, RUNS, NULL);
    else if (!strcmp(item, _("Edit attributes"))) 
	varinfo_dialog(mdata->active_var);
    else if (!strcmp(item, _("Delete"))) 
	delete_var_by_id(mdata->active_var);
    else if (!strcmp(item, _("Simulate..."))) 
	gretl_callback(NULL, SIM, NULL);
    else if (!strcmp(item, _("Define new variable..."))) 
	gretl_callback(NULL, GENR, NULL);
    gtk_widget_destroy(mdata->popup);
    return TRUE;
}

/* ........................................................... */

static GtkWidget *build_var_popup (void)
{
    static char *var_items[]={
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
	N_("Edit attributes"),
	N_("Delete"),
	N_("Simulate..."),
	N_("Define new variable...")
    };

    GtkWidget *var_menu;
    GtkWidget *var_item;
    int i, n_items = sizeof var_items / sizeof var_items[0];

    var_menu = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (i == 2 && !dataset_is_time_series(datainfo) &&
	    datainfo->time_series != STACKED_TIME_SERIES) {
	    continue;
	}
	if (!dataset_is_time_series(datainfo) && 
	    (i == 6 || i == 7 || i == 8 || i == 9)) {
	    continue;
	}
	var_item = gtk_menu_item_new_with_label(_(var_items[i]));
	gtk_signal_connect(GTK_OBJECT(var_item), "activate",
			   (GtkSignalFunc) popup_activated,
			   (gpointer) _(var_items[i]));
	GTK_WIDGET_SET_FLAGS (var_item, GTK_SENSITIVE | GTK_CAN_FOCUS);
	gtk_widget_show(var_item);
	gtk_menu_append(GTK_MENU(var_menu), var_item);
    }
    return var_menu;
}

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
    else if (!strcmp(item, _("Copy to clipboard"))) 
	csv_selected_to_clipboard();
    return TRUE;
}

static void build_selection_popup (void)
{
    const char *items[] = {
	N_("Display values"),
	N_("Descriptive statistics"),
	N_("Correlation matrix"),
	N_("Time series plot"),
	N_("Copy to clipboard"),
    };

    GtkWidget *item;
    int i, n_items = sizeof items / sizeof items[0];

    selection_popup = gtk_menu_new();

    for (i=0; i<n_items; i++) {
	if (!dataset_is_time_series(datainfo) && i == 3) {
	    continue;
	}
	item = gtk_menu_item_new_with_label(_(items[i]));
	gtk_signal_connect(GTK_OBJECT(item), "activate",
			   GTK_SIGNAL_FUNC(selection_popup_click),
			   _(items[i]));
	gtk_widget_show(item);
	gtk_menu_append(GTK_MENU(selection_popup), item);
    }
}

/* ........................................................... */

static void check_for_extra_data (void)
{
    DIR *dir;
    extern char pwtpath[MAXLEN]; /* datafiles.c */
    extern char jwpath[MAXLEN];  /* datafiles.c */
    extern char dgpath[MAXLEN];  /* datafiles.c */
    int gotpwt = 0, gotwool = 0, gotguj = 0;

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

    if (!gotpwt) *pwtpath = 0;

    /* then check for Wooldridge data */
    build_path(paths.datadir, "wooldridge", jwpath, NULL); 
    /* try at system level */
    if ((dir = opendir(jwpath)) != NULL) {
	closedir(dir);
	gotwool = 1;
    } else {
	build_path(paths.userdir, "wooldridge", jwpath, NULL); 
	/* and at user level */
	if ((dir = opendir(jwpath)) != NULL) {
	    closedir(dir);
	    gotwool = 1;
	}
    }

    if (!gotwool) *jwpath = 0;

    /* and Gujarati data */
    build_path(paths.datadir, "gujarati", dgpath, NULL); 
    /* try at system level */
    if ((dir = opendir(dgpath)) != NULL) {
	closedir(dir);
	gotguj = 1;
    } else {
	build_path(paths.userdir, "gujarati", dgpath, NULL); 
	/* and at user level */
	if ((dir = opendir(dgpath)) != NULL) {
	    closedir(dir);
	    gotguj = 1;
	}
    }

    if (!gotguj) *dgpath = 0;
}

/* ........................................................... */

void restore_sample (void)
{
    int err;

    err = restore_full_sample(&subZ, &fullZ, &Z,
			      &subinfo, &fullinfo, &datainfo);
    if (err) {
	gui_errmsg(err);
	return;
    }
}

static void restore_sample_callback (gpointer p, int verbose, GtkWidget *w)
{
    restore_sample();

    if (verbose) {
	infobox(_("Full sample range restored"));
	set_sample_label(datainfo);    
	restore_sample_state(FALSE);
	strcpy(line, "smpl full");
	verify_and_record_command(line);
    }
}

void gretl_fork (const char *prog, const char *arg)
{
    pid_t pid;

    signal(SIGCHLD, SIG_IGN);

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

static void startR (gpointer p, guint opt, GtkWidget *w)
{
    char Rprofile[MAXLEN], Rdata[MAXLEN], line[MAXLEN];
    const char *suppress = "--no-init-file";
    FILE *fp;
    int i;
    char *s0, *s1, *s2;
    pid_t pid;

    if (!data_status) {
	errbox(_("Please open a data file first"));
	return;
    }

    build_path(paths.userdir, "gretl.Rprofile", Rprofile, NULL);
    fp = fopen(Rprofile, "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }

    if (setenv("R_PROFILE", Rprofile, 1)) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 

    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);
    sprintf(line, "store \"%s\" -r", Rdata); 
    if (verify_and_record_command(line) ||
	write_data(Rdata, command.list, Z, datainfo, GRETL_DATA_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }

    if (dataset_is_time_series(datainfo)) {
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rdata);
    } else {
	char Rtmp[MAXLEN];
	FILE *fq;

	build_path(paths.userdir, "Rtmp", Rtmp, NULL);
	fq = fopen(Rtmp, "w");
	fprintf(fq, "gretldata <- read.table(\"%s\")\n", Rdata);
	fprintf(fq, "attach(gretldata)\n");
	fclose(fq);

	fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rtmp);
    }

    fclose(fp);

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

    signal(SIGCHLD, SIG_IGN); 
    pid = fork();

    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return;
    } else if (pid == 0) {  
#if 1
	if (i == 1)
	    execlp(s0, s0, suppress, NULL);
	else if (i == 2)
	    execlp(s0, s0, s1, suppress, NULL);
	else if (i == 3)
	    execlp(s0, s0, s1, s2, suppress, NULL);
#else
	if (i == 1)
	    execlp(s0, s0, NULL);
	else if (i == 2)
	    execlp(s0, s0, s1, NULL);
	else if (i == 3)
	    execlp(s0, s0, s1, s2, NULL);
#endif
	perror("execlp");
	_exit(EXIT_FAILURE);
    }
    free(s0); 
    free(s1); 
    free(s2);
}

/* ........................................................... */

static void show_calc (void)
{
    gretl_fork(calculator, NULL);
}

/* ........................................................... */

static void show_edit (void)
{
    gretl_fork(editor, NULL);
}

/* ........................................................... */

static void open_textbook_data (void)
{
    display_files(NULL, TEXTBOOK_DATA, NULL);
}

/* ........................................................... */

void show_toolbar (void)
{
    make_toolbar(mdata->w, main_vbox);
}

/* ........................................................... */

static void netscape_open (const char *url)
{
#ifdef USE_GNOME
    gnome_url_show(url);   
#else
    int err;
    char ns_cmd[128];

    sprintf(ns_cmd, "netscape -remote \"openURLNewWindow(%s)\"", url);
    err = system(ns_cmd);
    if (err) gretl_fork("netscape", url);
#endif /* USE_GNOME */
}

static void gretl_website (void)
{
    netscape_open("http://gretl.sourceforge.net/");
}

static void gretl_pdf (void)
{
    netscape_open("http://gretl.sourceforge.net/manual.pdf");
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

static void make_toolbar (GtkWidget *w, GtkWidget *box)
{
    GtkWidget *iconw, *button;
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
    toolbar_box = gtk_handle_box_new();
    gtk_handle_box_set_shadow_type(GTK_HANDLE_BOX(toolbar_box), 
				   GTK_SHADOW_NONE);
    gtk_box_pack_start(GTK_BOX(box), toolbar_box, FALSE, FALSE, 0);

    gretl_toolbar = gtk_toolbar_new(GTK_ORIENTATION_HORIZONTAL,
				    GTK_TOOLBAR_ICONS);
    gtk_container_set_border_width(GTK_CONTAINER(gretl_toolbar), 0);
    gtk_toolbar_set_space_size(GTK_TOOLBAR(gretl_toolbar), 0);
    gtk_container_add(GTK_CONTAINER(toolbar_box), gretl_toolbar);

    colorize_tooltips(GTK_TOOLBAR(gretl_toolbar)->tooltips);

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
	    toolfunc = show_gretl_console;
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
	    toolfunc = add_graph_to_session;
#endif
	    break;
	case 9:
	    toolxpm = mini_ofolder_xpm;
	    toolfunc = open_textbook_data;
	    break;
	default:
	    break;
	}

#ifdef GNUPLOT_PNG
	if (i == 8) continue;
#endif

	icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, 
						     NULL, toolxpm);
	iconw = gtk_pixmap_new(icon, mask);
	toolstr = _(toolstrings[i]);
	button = gtk_toolbar_append_item(GTK_TOOLBAR(gretl_toolbar),
					 NULL, toolstr, NULL,
					 iconw, toolfunc, NULL);
    }
    gtk_widget_show(gretl_toolbar);
    gtk_widget_show(toolbar_box);
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

#if 0
    if (info == GRETL_POINTER) {
	fprintf(stderr, "drag data: info = GRETL_POINTER\n");
    }
    if (info == GRETL_FILENAME) {
	fprintf(stderr, "drag data: info = GRETL_FILENAME\n");
    }
#endif

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
#if 0
    g_free(str);
    clipboard_buf = NULL;
#endif
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
    gtk_signal_connect (GTK_OBJECT(mdata->w), "selection_get",
			GTK_SIGNAL_FUNC(special_selection_get), NULL);    
}

/* ........................................................... */

static void auto_store (void)
{
    int opt = 0;

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

static void sort_varlist (gpointer p, guint col, GtkWidget *w)
{
    gtk_clist_set_compare_func(GTK_CLIST(mdata->listbox), NULL);
    gtk_clist_set_sort_column(GTK_CLIST(mdata->listbox), col);
    gtk_clist_sort(GTK_CLIST(mdata->listbox));
}
