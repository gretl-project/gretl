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

/* gretl.c : main for gretl */

#include "gretl.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "libset.h"
#include "gretl_bundle.h"

#include "treeutils.h"
#include "ssheet.h"
#include "console.h"
#include "session.h"
#include "database.h"
#include "datafiles.h"
#include "cmdstack.h"
#include "filelists.h"
#include "toolbar.h"
#include "menustate.h"
#include "fileselect.h"
#include "filters.h"
#include "calculator.h"
#include "fnsave.h"
#include "winstack.h"
#include "datawiz.h"
#include "varinfo.h"
#include "dlgutils.h"
#include "fncall.h"
#include "selector.h"

#include <dirent.h>

#ifndef G_OS_WIN32
# include <unistd.h>
# include <sys/types.h>
# include <signal.h>
# include "../pixmaps/gretl.xpm"  /* program icon for X */
#else
# include <windows.h>
# include "gretlwin32.h"
#endif

#define GUI_DEBUG 0

#if GUI_DEBUG
# include "version.h"
# include "build.h"
#endif

/* update.c */
extern int silent_update_query (void);
extern int update_query (void); 

/* functions private to gretl.c */
static void make_main_window (void);

static gboolean main_popup_handler (GtkWidget *w, GdkEventButton *event,
				    gpointer data);
static int set_up_main_menu (void);
static void start_R_callback (void);
static void auto_store (void);
static void restore_sample_callback (void);
static void show_sample_callback (void);
static void mdata_select_all (void);
static void mdata_select_list (void);
static int gui_query_stop (void);

GtkTargetEntry gretl_drag_targets[] = {
    { "text/uri-list",  0, GRETL_FILENAME },
    { "db_series_ptr",  GTK_TARGET_SAME_APP, GRETL_DBSERIES_PTR },
    { "model_ptr",      GTK_TARGET_SAME_APP, GRETL_MODEL_PTR },
    { "remote_db_ptr",  GTK_TARGET_SAME_APP, GRETL_REMOTE_DB_PTR },
    { "remote_pkg_ptr", GTK_TARGET_SAME_APP, GRETL_REMOTE_FNPKG_PTR },
    { "graph_file",     GTK_TARGET_SAME_APP, GRETL_GRAPH_FILE }
};

static void  
mdata_handle_drag  (GtkWidget          *widget,
		    GdkDragContext     *dc,
		    gint                x,
		    gint                y,
		    GtkSelectionData   *data,
		    guint               info,
		    guint               time,
		    gpointer            p);

static char *optrun, *optdb, *optwebdb, *optpkg;
static int opteng, optbasque, optdump, optver;
#ifdef G_OS_WIN32
static int optdebug;
#endif

static gchar *param_msg = 
    N_("\nYou may supply the name of a data file on the command line");

static GOptionEntry options[] = {
    { "run", 'r', 0, G_OPTION_ARG_FILENAME, &optrun, 
      N_("open a script file on startup"), "SCRIPT" },
    { "db", 'd', 0, G_OPTION_ARG_STRING, &optdb, 
      N_("open a database on startup"), "DATABASE" },
    { "webdb", 'w', 0, G_OPTION_ARG_STRING, &optwebdb, 
      N_("open a remote (web) database on startup"), "REMOTE_DB" },
    { "pkg", 'p', 0, G_OPTION_ARG_STRING, &optpkg, 
      N_("open (edit) a function package on startup"), "FUNCPKG" },
    { "english", 'e', 0, G_OPTION_ARG_NONE, &opteng, 
      N_("force use of English"), NULL },
    { "basque", 'q', 0, G_OPTION_ARG_NONE, &optbasque, 
      N_("force use of Basque"), NULL },
    { "dump", 'c', 0, G_OPTION_ARG_NONE, &optdump, 
      N_("dump gretl configuration to file"), NULL },
#ifdef G_OS_WIN32
    { "debug", 'b', 0, G_OPTION_ARG_NONE, &optdebug, 
      N_("send debugging info to console"), NULL }, 
#endif
    { "version", 'v', 0, G_OPTION_ARG_NONE, &optver, 
      N_("print version information"), NULL }, 
    { NULL, '\0', 0, 0, NULL, NULL, NULL },
};

windata_t *mdata;
DATASET *dataset;

char datafile[MAXLEN];
char scriptfile[MAXLEN];
char tryfile[MAXLEN];

MODEL *model;             /* gretl models struct */

int data_status, orig_vars;
float gui_scale;

/* defaults for some options */
int updater = FALSE;
int winsize = FALSE;
int main_x = -1;
int main_y = -1;
int mainwin_width = 520;
int mainwin_height = 420;
int ox_support = FALSE;

#if defined(G_OS_WIN32)
char calculator[MAXSTR] = "calc.exe";
char latex[MAXSTR] = "pdflatex.exe";
char viewdvi[MAXSTR] = "windvi.exe";
char Rcommand[MAXSTR] = "RGui.exe";
#elif defined(OSX_BUILD)
char calculator[MAXSTR] = "/Applications/Calculator.app/Contents/MacOS/Calculator";
char latex[MAXSTR] = "pdflatex";
char viewdvi[MAXSTR] = "xdvi";
char Rcommand[MAXSTR] = "/usr/X11R6/bin/xterm -e R";
#else
char Browser[MAXSTR] = "mozilla";
char calculator[MAXSTR] = "xcalc";
char latex[MAXSTR] = "latex";
char viewdvi[MAXSTR] = "xdvi";
char viewpdf[MAXSTR] = "acroread";
char viewps[MAXSTR] = "gv";
char Rcommand[MAXSTR] = "xterm -e R";
#endif

static void spreadsheet_edit (void) 
{
    show_spreadsheet(SHEET_EDIT_VARLIST);
}

static void varinfo_callback (void)
{
    varinfo_dialog(mdata->active_var);
}

static void options_dialog_callback (void)
{
    options_dialog(0, NULL, mdata->main);
}

static void open_script_callback (void)
{
    file_selector(OPEN_SCRIPT, FSEL_DATA_NONE, NULL);
}

static void open_session_callback (void)
{
    file_selector(OPEN_SESSION, FSEL_DATA_NONE, NULL);
}

static void edit_package_callback (GtkAction *action, gpointer p)
{
    file_selector(OPEN_GFN, FSEL_DATA_NONE, NULL);
}

static void new_package_callback (GtkAction *action, gpointer p)
{
    if (no_user_functions_check()) {
	return;
    } else {
	fsave_callback(action, p);
    }
}

#ifdef ENABLE_MAILER
static void email_data (gpointer p, guint u, GtkWidget *w)
{
    send_file(datafile);
}
#endif

static void noalloc (void)
{
    fputs(I_("Out of memory!\n"), stderr);
    exit(EXIT_FAILURE);
}

static void get_runfile (char *fname)
{
    *tryfile = '\0';

#ifdef G_OS_WIN32
    if (filename_to_win32(tryfile, fname)) {
	return;
    }
#else
    strncat(tryfile, fname, MAXLEN - 1);
#endif

    if (gretl_addpath(tryfile, 1) == NULL) {
	fprintf(stderr, I_("Couldn't find script '%s'\n"), tryfile);
	exit(EXIT_FAILURE);
    } else {
	fprintf(stderr, I_("%s found\n"), tryfile);
	gretl_set_current_dir(tryfile);
    }
}

static int script_type (const char *fname)
{
    if (has_suffix(fname, ".inp")) {
	return EDIT_SCRIPT;
    } else if (has_suffix(fname, ".R")) {
	return EDIT_R;
    } else if (has_suffix(fname, ".plt") ||
	       has_suffix(fname, ".gp")) {
	return EDIT_GP;
    } else if (ox_support && has_suffix(fname, ".ox")) {
	return EDIT_OX;
    } else if (has_suffix(fname, ".m")) {
	return EDIT_OCTAVE;
    } else {
	return 0;
    }
}

static void maybe_fix_dbname (char *dbname)
{
    FILE *fp = NULL;

    if (strstr(dbname, ".bin") == NULL &&
	strstr(dbname, ".rat") == NULL &&
	strstr(dbname, ".RAT") == NULL &&
	strstr(dbname, ".bn7") == NULL) {
	strcat(dbname, ".bin");
    }

    fp = gretl_fopen(dbname, "rb");

    if (fp != NULL) {
	fclose(fp);
    } else if (!g_path_is_absolute(dbname)) {
	gchar *tmp = g_build_path(SLASHSTR, gretl_home(), "db",
				  dbname, NULL);

	fp = gretl_fopen(tmp, "rb");
	if (fp != NULL) {
	    fclose(fp);
	    strcpy(dbname, tmp);
	}
	g_free(tmp);
    }
}

#if !defined(ENABLE_NLS)

static void real_nls_init (void)
{
    return;
}

#elif defined(G_OS_WIN32) 

static void real_nls_init (void)
{
    char localedir[MAXSTR];

    build_path(localedir, gretl_home(), "locale", NULL);
    setlocale(LC_ALL, "");
    set_gretl_charset();
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#elif defined(OSX_BUILD) 

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

    setlocale(LC_ALL, "");
    set_gretl_charset();
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#else /* regular *nix treatment of NLS */

static void real_nls_init (void)
{
    setlocale(LC_ALL, "");
    set_gretl_charset();
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

#ifndef G_OS_WIN32

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

static void usr1_handler (int signal) 
{
    fprintf(stderr, "Hello from gretl's signal handler: got SIGUSR1\n");
}

#endif

/* callback from within potentially lengthy libgretl
   operations: try to avoid having the GUI become
   totally unresponsive
*/

static void gui_show_activity (void)
{
    while (gtk_events_pending()) {
	gtk_main_iteration();
    }
}

static int have_data (void)
{
    return dataset != NULL && dataset->v > 0;
}

int main (int argc, char **argv)
{
#ifdef G_OS_WIN32
    char *callname = argv[0];
#endif
    int ftype = 0;
    char auxname[MAXLEN];
    char filearg[MAXLEN];
    GError *opterr = NULL;

#ifdef G_OS_WIN32
    win32_set_gretldir(callname);
#else
    signal(SIGUSR1, usr1_handler);
#endif

    gui_nls_init();

    *tryfile = '\0';
    *scriptfile = '\0';
    *datafile = '\0';
    *auxname = '\0';
    *filearg = '\0';

#if GUI_DEBUG
    fprintf(stderr, "starting gretl %s, %s\n", GRETL_VERSION, BUILD_DATE);
#endif

    gtk_init_with_args(&argc, &argv, _(param_msg), options, "gretl", &opterr);
    if (opterr != NULL) {
	g_print("%s\n", opterr->message);
	exit(EXIT_FAILURE);
    }

#ifdef G_OS_WIN32
    /* let's call this before doing libgretl_init */
    gretl_win32_debug_init(optdebug);
#endif

    libgretl_init();

#ifdef G_OS_WIN32
    gretl_win32_init(callname, optdebug);
#else 
    gretl_config_init();
#endif

    if (optver) {
	gui_logo(NULL);
	exit(EXIT_SUCCESS);
    } else if (optdump) {
	dump_rc();
	exit(EXIT_SUCCESS);
    } else if (optrun) {
	get_runfile(optrun);
    } else if (optdb != NULL) {
	strncat(auxname, optdb, MAXLEN - 1);
	maybe_fix_dbname(auxname);
    } else if (optwebdb != NULL) {
	strncat(auxname, optdb, MAXLEN - 1);
    } else if (optpkg != NULL) {
	strncat(auxname, optpkg, MAXLEN - 1);
    }

    if (opteng) {
	force_language(LANG_C);
	force_english_help();
    } else if (optbasque) {
	force_language(LANG_EU);
    }

#if GUI_DEBUG
    fprintf(stderr, "finished option processing\n");
#endif

    /* set libgretl callbacks */
    set_workdir_callback(gui_set_working_dir);
    set_show_activity_func(gui_show_activity);
    set_query_stop_func(gui_query_stop);

    /* allocate data information struct */
    dataset = datainfo_new();
    if (dataset == NULL) {
	noalloc();
    }

    /* allocate memory for models */
    model = allocate_working_model();
    if (model == NULL) {
	noalloc();
    } 

    library_command_init();

    helpfile_init();
    session_init();
    init_fileptrs();

#if GUI_DEBUG
    fprintf(stderr, "finished miscellaneous init functions\n");
#endif

    if (argc > 1) {
	/* Process what is presumably a filename argument
	   given on the command line (by now any options will
	   have been extracted from the argv array).
	*/
	PRN *prn; 
	int err = 0;

	prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
	if (err) {
	    exit(EXIT_FAILURE);
	}

	if (*filearg == '\0') {
	    /* not already registered */
	    strncat(filearg, argv[1], MAXLEN - 1);
	}

	*datafile = '\0';

#ifdef G_OS_WIN32
	if (filename_to_win32(datafile, filearg)) {
	    exit(EXIT_FAILURE);
	}
#else
	record_filearg(datafile, filearg);
#endif

	/* keep a copy of input filename */
	strcpy(tryfile, datafile);

	ftype = detect_filetype(datafile, OPT_P);

	switch (ftype) {
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(datafile, dataset, OPT_NONE, prn);
	    break;
	case GRETL_XML_DATA:
	    err = gretl_read_gdt(datafile, dataset, OPT_NONE, prn);
	    break;
	case GRETL_CSV:
	    err = import_csv(datafile, dataset, OPT_NONE, prn);
	    break;
	case GRETL_XLS:
	case GRETL_XLSX:    
	case GRETL_GNUMERIC:
	case GRETL_ODS:
	case GRETL_DTA:
	case GRETL_SAV:
	case GRETL_SAS:
	case GRETL_JMULTI:
	case GRETL_OCTAVE:
	case GRETL_WF1:
	    err = get_imported_data(datafile, ftype, 0);
	    break;
	case GRETL_SCRIPT:
	case GRETL_SESSION:
	    get_runfile(datafile);
	    *datafile = '\0';
	    break;
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:  
	case GRETL_PCGIVE_DB:
	    strcpy(auxname, datafile);
	    *tryfile = '\0';
	    *datafile = '\0';
	    maybe_fix_dbname(auxname);
	    optdb = auxname;
	    break;
	case GRETL_UNRECOGNIZED:
	default:
	    fprintf(stderr, "%s: unrecognized file type", tryfile);
	    exit(EXIT_FAILURE);
	    break;
	}

	if (err == E_CANCEL) {
	    err = 0;
	    ftype = 0;
	    *tryfile = '\0';
	    *datafile = '\0';
	}

	if (ftype != GRETL_SCRIPT && err) {
	    errmsg(err, prn);
	    exit(EXIT_FAILURE);
	}

	gretl_print_destroy(prn);
    }

#if GUI_DEBUG
    fprintf(stderr, "about to build GUI...\n");
#endif

    /* create the GUI */
    gretl_stock_icons_init();
#if GUI_DEBUG
    fprintf(stderr, " done gretl_stock_icons_init\n");
#endif
    make_main_window();

#if GUI_DEBUG
    fprintf(stderr, "done make_main_window\n");
#endif

    if (have_data()) {
	/* redundant? */
	set_sample_label(dataset);
    }

    add_files_to_menus();

#if GUI_DEBUG
    fprintf(stderr, "done add_files_to_menus\n");
#endif

    session_menu_state(FALSE);
    restore_sample_state(FALSE);
    dataset_menubar_state(FALSE);

#if GUI_DEBUG
    fprintf(stderr, "done setting GUI state\n");
#endif

    /* FIXME run init script, if found? */

    if (have_data()) {
	register_startup_data(tryfile);
	maybe_display_string_table();
	*tryfile = '\0';
    }

    /* opening a script or session from the command line? */
    if (*tryfile != '\0') { 
	if (gretl_is_pkzip_file(tryfile)) {
	    ftype = GRETL_SESSION;
	}
	if (ftype == GRETL_SESSION) {
	    do_open_session();
	} else if ((ftype = script_type(tryfile))) {
	    do_open_script(ftype);
	} else {
	    do_open_script(EDIT_SCRIPT);
	}
    }

    /* check for program updates? */
    if (updater) {
	silent_update_query(); 
    }

#if 0
    fprintf(stderr, " at last step");
    if (optdb != NULL) {
	fprintf(stderr, ": optdb != NULL\n");
    } else if (optwebdb != NULL) {
	fprintf(stderr, ": optwebdb != NULL\n");
    } else if (optpkg != NULL) {
	fprintf(stderr, ": optpkg != NULL\n");
    } else {
	fputc('\n', stderr);
    }
#endif

    /* try opening specified database or package */
    if (optdb != NULL) {
	open_named_db_index(auxname);
    } else if (optwebdb != NULL) {
	open_named_remote_db_index(auxname);
    } else if (optpkg != NULL) {
	edit_package_at_startup(auxname);
    }

#if GUI_DEBUG
    fprintf(stderr, "calling gtk_main()\n");
#endif

    /* Enter the event loop */
    gtk_main();

    /* clean up before exiting */
    free_session();

    destroy_working_model(model);

    library_command_free();
    libgretl_cleanup();

    if (data_status) {
	destroy_dataset(dataset);
    }

    destroy_file_collections();
    free_command_stack();

    return EXIT_SUCCESS;
}

static void check_varmenu_state (GtkTreeSelection *select, gpointer p)
{
    if (mdata->ui != NULL) {
	int vnum = 0;
	int sc = tree_selection_count(select, &vnum);

	if (sc == 1 && vnum > 0) {
	    mdata->active_var = vnum;
	    maybe_reset_varinfo_dialog();
	}

	variable_menu_state(sc == 1);
    }
}

/* if a keystroke (e.g. page up) would take us to row 0, countermand
   this and go to row 1 instead
*/

static void mdata_avoid_zero (GtkTreeView *view, gpointer p)
{
    GtkTreePath *path = NULL;
    int i;

    gtk_tree_view_get_cursor(view, &path, NULL);

    if (path != NULL) {
	i = gtk_tree_path_get_indices(path)[0];
	if (i == 0) {
	    GtkTreePath *newp;

	    newp = gtk_tree_path_new_from_indices(1, -1);
	    gtk_tree_view_set_cursor(view, newp, NULL, FALSE);
	    gtk_tree_path_free(newp);
	}
	gtk_tree_path_free(path);
    }
}

/* keystrokes recognized in the main gretl window */

static gint catch_mdata_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods = widget_get_pointer_mask(w);
    int k = key->keyval;

    if (k == GDK_h || k == GDK_F1) {
	/* invoke help */
	plain_text_cmdref(NULL);
	return TRUE;
    } else if (k == GDK_g) {
	/* invoke genr */
	genr_callback();
	return TRUE;
    } else if (k == GDK_c) {
	/* launch the console */
	gretl_console();
	return TRUE;
    } else if (mods & GDK_MOD1_MASK) {
	if (k == GDK_x) {
	    /* Alt-x: invoke command minibuffer */
	    minibuf_callback();
	    return TRUE;
	} else if (k == GDK_w) {
	    /* Alt-w: window list */
	    window_list_popup(vwin->main);
	    return TRUE;
	}
    }

    if (dataset->v == 0) {
	goto suppress;
    }

#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    if (!(mods & GDK_MOD1_MASK)) {
	if (k == GDK_a) {
	    audio_render_window(vwin, AUDIO_LISTBOX);
	    return TRUE;
	} else if (k == GDK_x) {
	    stop_talking();
	    return TRUE;
	}
    }
#endif

    if (k == GDK_r) {
	refresh_data();
	return TRUE;
    }

    if (k == GDK_Return              /* display variable(s) */
	|| k == GDK_Delete           /* delete variable(s) */
	|| k == GDK_e || k == GDK_F2 /* edit variable's info */
	|| k == GDK_t                /* graph variable */
	) {
	int selcount, vnum = 0;

	selcount = vwin_selection_count(mdata, &vnum);

	if (selcount == 1 && vnum != 0) {
	    mdata->active_var = vnum;
	    if (k == GDK_e || k == GDK_F2) {
		varinfo_dialog(mdata->active_var);
	    } else if (k == GDK_t) {
		do_graph_var(mdata->active_var);
	    } else if (k == GDK_Return) {
		display_var();
	    } else if (k == GDK_Delete) {
		delete_single_var(mdata->active_var);
	    }
	} else if (selcount > 1) {
	    if (k == GDK_Delete) {
		delete_selected_vars();
	    } else if (k == GDK_Return) {
		display_selected();
	    }
	}

	return TRUE;
    } 

 suppress:

    /* suppress echo of useless keystrokes */
    if (k != GDK_Up && k != GDK_Down &&
	k != GDK_Page_Up && k != GDK_Page_Down &&
	k != GDK_Home && k != GDK_End) {
	return TRUE;
    }

    return FALSE;
}

static int lagvar_get_parent_iter (int pv, GtkTreeIter *parent)
{
    GtkTreeModel *model = 
	gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox));
    GtkTreeIter iter;
    gchar *idstr;
    int ret = 0;

    if (!gtk_tree_model_get_iter_first(model, &iter)) {
	return 0;
    }
    
    while (1) {
	gtk_tree_model_get(model, &iter, 0, &idstr, -1);
	if (atoi(idstr) == pv) {
	    *parent = iter;
	    ret = 1;
	}
	g_free(idstr);
	if (ret || !gtk_tree_model_iter_next(model, &iter)) {
	    break;
	}
    }

    return ret;
}

static void mdata_select_all (void)
{
    GtkTreeSelection *select;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_select_all(select);
}

static int get_line_pos (GtkTreeModel *mod)
{
    GtkTreeIter iter;
    gchar *idstr;
    int i = 1, pos = 0;

    if (gtk_tree_model_get_iter_first(mod, &iter)) {
	while (gtk_tree_model_iter_next(mod, &iter)) {
	    gtk_tree_model_get(mod, &iter, 0, &idstr, -1);
	    if (idstr != NULL && atoi(idstr) == mdata->active_var) {
		pos = i;
	    }
	    g_free(idstr);
	    if (pos) {
		break;
	    }
	    i++;
	}
    }

    return pos;
}

/* populate the list of series in the main gretl window */

void populate_varlist (void)
{
    static gint check_connected;
    static gint click_connected;
    GtkTreeView *view = GTK_TREE_VIEW(mdata->listbox);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeStore *store = GTK_TREE_STORE(model);
    GtkTreeSelection *select;
    GtkTreeIter iter;    
    char id[8];
    int i, pos = 0;

    if (store != NULL) {
	/* record line position? */
	pos = get_line_pos(model);
    }
    
    gtk_tree_store_clear(store);
    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<dataset->v; i++) {
	int pv = 0;

	if (var_is_hidden(dataset, i)) {
	    continue;
	}
	if (i > 0 && (is_standard_lag(i, dataset, &pv) ||
		      is_dummy_child(i, dataset, &pv))) {
	    if (pv > 0) {
		GtkTreeIter child_iter, parent_iter;

		if (lagvar_get_parent_iter(pv, &parent_iter)) {
		    gtk_tree_store_insert_before(store, &child_iter, 
						 &parent_iter, NULL);
		    sprintf(id, "%d", i);
		    gtk_tree_store_set(store, &child_iter, 
				       0, id, 
				       1, dataset->varname[i],
				       2, VARLABEL(dataset, i),
				       -1);	
		}	
	    }
	}

	if (pv == 0) {
	    gtk_tree_store_append(store, &iter, NULL);
	    sprintf(id, "%d", i);
	    gtk_tree_store_set(store, &iter, 
			       0, id, 
			       1, dataset->varname[i],
			       2, VARLABEL(dataset, i),
			       -1);
	}
    } 

    gtk_tree_model_get_iter_first(model, &iter);

    if (pos == 0) {
	/* no saved position */
	pos = 1;
	gtk_tree_model_iter_next(model, &iter);
    } else {	
	/* try to return to previous position */
	GtkTreeIter last;

	i = 1;
	while (1) {
	    last = iter;
	    if (!gtk_tree_model_iter_next(model, &iter)) {
		/* reached the end! */
		iter = last;
		break;
	    } else if (i == pos) {
		/* found previous position */
		break;
	    }
	    i++;
	}
	pos = i;
    } 

    mdata->active_var = pos;
    select = gtk_tree_view_get_selection(view);
    gtk_tree_selection_select_iter(select, &iter);

    if (dataset->v > 1) {
	GtkTreePath *path;

	sprintf(id, "%d", pos);
	path = gtk_tree_path_new_from_string(id);
	gtk_tree_view_set_cursor(view, path, NULL, FALSE);
	gtk_tree_path_free(path);
    }

    if (!check_connected) {
	g_signal_connect(G_OBJECT(select), "changed",
			 G_CALLBACK(check_varmenu_state),
			 mdata);
	g_signal_connect(G_OBJECT(mdata->listbox), "cursor-changed",
			 G_CALLBACK(mdata_avoid_zero),
			 NULL);
	check_connected = 1;
    }

    if (!click_connected) {
	g_signal_connect(G_OBJECT(mdata->listbox), "button-press-event",
			 G_CALLBACK(main_popup_handler), 
			 mdata);
	g_signal_connect(G_OBJECT(mdata->listbox), "button-press-event",
			 G_CALLBACK(main_varclick),
			 mdata);
	click_connected = 1;
    }

    variable_menu_state(TRUE);
}

void mdata_select_last_var (void)
{
    GtkTreeIter iter, last;
    GtkTreeView *view;
    GtkTreeModel *model;
    GtkTreeSelection *select;

    view = GTK_TREE_VIEW(mdata->listbox);
    model = gtk_tree_view_get_model(view);
    gtk_tree_model_get_iter_first(model, &iter);

    while (1) {
	last = iter;
	if (!gtk_tree_model_iter_next(model, &iter)) {
	    iter = last;
	    break;
	}
    }

    select = gtk_tree_view_get_selection(view);
    gtk_tree_selection_unselect_all(select);
    gtk_tree_selection_select_iter(select, &iter);
}

static void real_select_list (const int *list)
{
    GtkTreeIter iter;
    GtkTreeView *view;
    GtkTreeModel *model;
    GtkTreeSelection *select;
    gchar *idstr;
    int nsel = 0;

    view = GTK_TREE_VIEW(mdata->listbox);
    model = gtk_tree_view_get_model(view);
    select = gtk_tree_view_get_selection(view);
    gtk_tree_selection_unselect_all(select);

    gtk_tree_model_get_iter_first(model, &iter);

    while (nsel < list[0]) {
	if (!gtk_tree_model_iter_next(model, &iter)) {
	    break;
	}	
	gtk_tree_model_get(model, &iter, 0, &idstr, -1);
	if (in_gretl_list(list, atoi(idstr))) {
	    gtk_tree_selection_select_iter(select, &iter);
	    nsel++;
	} 
	g_free(idstr);
    }
}

static void mdata_select_list (void)
{
    int nl = n_saved_lists();

    if (nl == 0) {
	warnbox(_("No lists are currently defined"));
	return;
    } else {
	char lname[32];
	int resp;

	resp = select_list_dialog(lname);

	if (!canceled(resp)) {
	    int *list = get_list_by_name(lname);

	    if (list != NULL) {
		real_select_list(list);
	    }
	}
    }
}

void clear_varlist (GtkWidget *widget)
{
    GtkTreeModel *model;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(widget));

    if (GTK_IS_TREE_STORE(model)) {
	gtk_tree_store_clear(GTK_TREE_STORE(model));
    } else if (GTK_IS_LIST_STORE(model)) {
	gtk_list_store_clear(GTK_LIST_STORE(model));
    }
}

static float get_gui_scale (void)
{
    GtkSettings *settings = gtk_settings_get_default();
    float scale = 1.0;

    if (settings != NULL) {
	gchar *fontname = NULL;
	int fsize;

	g_object_get(G_OBJECT(settings), "gtk-font-name", &fontname, NULL);

	if (fontname != NULL) {
	    if (sscanf(fontname, "%*s %d", &fsize) == 1) {
		if (fsize > 10 && fsize < 100) {
		    scale = fsize / 10.0;
		}
	    }
	    g_free(fontname);
	}
    }

    return scale;
}

static gboolean 
mainwin_config (GtkWidget *w, GdkEventConfigure *event, gpointer p)
{
    mainwin_width = event->width;
    mainwin_height = event->height;

    gdk_window_get_root_origin(gtk_widget_get_window(mdata->main), 
			       &main_x, &main_y);

    return FALSE;
}

/* scale up the main window if it seems to be too tiny in relation to
   the screen dimensions
*/

static void set_main_window_scale (void)
{
    GdkScreen *s = gdk_screen_get_default();

    if (s != NULL) {
	int w = gdk_screen_get_width(s);
	int h = gdk_screen_get_height(s);
	double aspect = 1.25;
	double hfac = 2.10;
 
	if (mainwin_height < h / hfac) {
	    mainwin_height = h / hfac;
	    if ((double) w / h > 1.35) {
		/* widescreen */
		aspect = 1.4;
	    }
	    w = aspect * mainwin_height;
	    if (mainwin_width < w) {
		mainwin_width = w;
	    }
	}
    }
}

static void make_main_window (void) 
{
    GtkWidget *main_vbox;
    GtkWidget *box, *dlabel;
    GtkWidget *align;
    const char *titles[] = {
	N_("ID #"), 
	N_("Variable name"), 
	N_("Descriptive label")
    };
    GType types[] = {
	G_TYPE_STRING,
	G_TYPE_STRING,
	G_TYPE_STRING
    };

    mdata = gretl_viewer_new(MAINWIN, "gretl", NULL, 0);
    if (mdata == NULL) {
	noalloc();
    }
    
    gui_scale = get_gui_scale();

    if (!winsize || mainwin_width <= 200 || mainwin_height <= 200) {
	/* set default window size */
	mainwin_width = 580 * gui_scale;
	mainwin_height = 420 * gui_scale;
	set_main_window_scale();
    }

    mdata->main = gtk_window_new(GTK_WINDOW_TOPLEVEL);

#ifdef G_OS_WIN32
    set_up_windows_look();
#endif

    g_signal_connect(G_OBJECT(mdata->main), "configure-event",
		     G_CALLBACK(mainwin_config), NULL);
    g_signal_connect(G_OBJECT(mdata->main), "delete-event",
		     G_CALLBACK(exit_check), NULL);
    g_signal_connect(G_OBJECT(mdata->main), "destroy",
		     G_CALLBACK(gtk_main_quit), NULL);

    gtk_window_set_title(GTK_WINDOW(mdata->main), "gretl");
    gtk_window_set_default_size(GTK_WINDOW(mdata->main), 
				mainwin_width, mainwin_height);

#ifndef G_OS_WIN32
    set_wm_icon(mdata->main);
#endif

    main_vbox = gtk_vbox_new(FALSE, 4);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 4);
    gtk_container_add(GTK_CONTAINER(mdata->main), main_vbox);
    g_object_set_data(G_OBJECT(mdata->main), "vbox", main_vbox);

    if (set_up_main_menu()) {
	exit(EXIT_FAILURE);
    }

    /* put the main menu bar in place */
    gtk_box_pack_start(GTK_BOX(main_vbox), mdata->mbar, FALSE, TRUE, 0);

    /* label for name of datafile */
    dlabel = gtk_label_new(_(" No datafile loaded ")); 
    g_object_set_data(G_OBJECT(mdata->main), "dlabel", dlabel);

    /* this will hold the list of variables */
    box = gtk_vbox_new(FALSE, 0);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(box), align, FALSE, FALSE, 0);
    gtk_container_add(GTK_CONTAINER(align), dlabel);

#if GUI_DEBUG
    fprintf(stderr, " adding main-window listbox...\n");
#endif
   
    vwin_add_list_box(mdata, GTK_BOX(box), 3, FALSE, types, titles, 1);

    gtk_drag_dest_set(mdata->listbox,
		      GTK_DEST_DEFAULT_ALL,
		      gretl_drag_targets, 2,
		      GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(mdata->listbox), "drag-data-received",
		     G_CALLBACK(mdata_handle_drag),
		     NULL);

    g_signal_connect(G_OBJECT(mdata->listbox), "key-press-event",
		     G_CALLBACK(catch_mdata_key),
		     mdata);

    gtk_box_pack_start(GTK_BOX(main_vbox), box, TRUE, TRUE, 0);
    mdata->status = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(main_vbox), mdata->status, FALSE, TRUE, 0);

#if GUI_DEBUG
    fprintf(stderr, " finalizing main window...\n");
#endif

    /* put stuff into list box, activate menus */
    if (have_data()) {
	populate_varlist();
    }

#if GUI_DEBUG
    fprintf(stderr, "  step 1 done\n");
#endif

    /* get a monospaced font for various windows */
    set_fixed_font();

#if GUI_DEBUG
    fprintf(stderr, "  set_fixed_font done\n");
#endif

    /* and a proportional font for menus, etc. */
    set_app_font(NULL);

#if GUI_DEBUG
    fprintf(stderr, "  set_app_font done\n");
#endif

    add_mainwin_toolbar(main_vbox);

#if GUI_DEBUG
    fprintf(stderr, "  add_mainwin_toolbar done\n");
#endif

    gtk_widget_show_all(mdata->main); 

#if GUI_DEBUG
    fprintf(stderr, "  gtk_widget_show_all done\n");
#endif

    if (winsize && main_x >= 0 && main_y >= 0) {
	gtk_window_move(GTK_WINDOW(mdata->main), main_x, main_y);
    }

#if GUI_DEBUG
    fprintf(stderr, "  possible window_move done\n");
#endif
}

GtkActionEntry main_entries[] = {
    /* File */
    { "File",         NULL, N_("_File"), NULL, NULL, NULL }, 
    { "OpenData",     NULL, N_("_Open data"), NULL, NULL, NULL }, 
    { "OpenGdt",      GTK_STOCK_OPEN, N_("_User file..."), NULL, NULL, G_CALLBACK(open_data) },
    { "DisplayDataFiles", GTK_STOCK_OPEN, N_("_Sample file..."), "", NULL, G_CALLBACK(show_files) },
    { "ImportData",   NULL, N_("_Import"), NULL, NULL, NULL }, 
    { "OpenCSV",      NULL, N_("_text/CSV..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenOctave",   NULL, N_("_Octave..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenGnumeric", NULL, N_("_Gnumeric..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenXLS",      NULL, N_("_Excel..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenODS",      NULL, N_("_Open Document..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenWF1",      NULL, N_("_Eviews..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenDTA",      NULL, N_("_Stata..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenSAV",      NULL, N_("_SPSS..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenSAS",      NULL, N_("_SAS (xport)..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenJMulTi",   NULL, N_("_JMulTi..."), NULL, NULL, G_CALLBACK(open_data) },

    { "AppendData",     NULL, N_("_Append data"), NULL, NULL, NULL }, 
    { "AppendGdt",      NULL, N_("_Standard format..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendCSV",      NULL, N_("_text/CSV..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendOctave",   NULL, N_("_Octave..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendGnumeric", NULL, N_("_Gnumeric..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendXLS",      NULL, N_("_Excel..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendODS",      NULL, N_("_Open Document..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendWF1",      NULL, N_("_Eviews..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendDTA",      NULL, N_("_Stata..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendSAV",      NULL, N_("_SPSS..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendSAS",      NULL, N_("_SAS (xport)..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendJMulTi",   NULL, N_("_JMulTi..."), NULL, NULL, G_CALLBACK(open_data) },

    { "SaveData",  GTK_STOCK_SAVE, N_("_Save data"), "<control>S", NULL, G_CALLBACK(auto_store) },
    { "SaveDataAs", GTK_STOCK_SAVE_AS, N_("Save data _as..."), NULL, NULL, G_CALLBACK(fsave_callback) }, 
    { "ExportData", NULL, N_("_Export data..."), NULL, NULL, G_CALLBACK(fsave_callback) },

    { "MailData", GRETL_STOCK_MAIL, N_("Send To..."), NULL, NULL, G_CALLBACK(email_data) },
    { "NewData", GTK_STOCK_NEW, N_("_New data set"), NULL, NULL, G_CALLBACK(newdata_callback) },
    { "ClearData", GTK_STOCK_CLEAR, N_("C_lear data set"), NULL, NULL, G_CALLBACK(verify_clear_data) },
    { "WorkingDir", NULL, N_("_Working directory..."), NULL, NULL, G_CALLBACK(working_dir_dialog) },
    { "ScriptFiles", NULL, N_("_Script files"), NULL, NULL, NULL },
    { "OpenScript", GTK_STOCK_OPEN, N_("_User file..."), "", NULL, G_CALLBACK(open_script_callback) },
    { "DisplayScripts", GTK_STOCK_OPEN, N_("_Practice file..."), "", NULL, G_CALLBACK(show_files) },
    { "NewScript", GTK_STOCK_NEW, N_("_New script"), "", NULL, NULL },
    { "GretlScript", NULL, N_("gretl script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "GnuplotScript", NULL, N_("gnuplot script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "RScript", NULL, N_("R script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "OctaveScript", NULL, N_("Octave script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "OxScript", NULL, N_("Ox program"), NULL, NULL, G_CALLBACK(new_script_callback) },

    { "SessionFiles", NULL, N_("_Session files"), NULL, NULL, NULL },
    { "OpenSession", GTK_STOCK_OPEN, N_("_Open session..."), "", NULL, G_CALLBACK(open_session_callback) },
    { "SaveSession", GTK_STOCK_SAVE, N_("_Save session"), "", NULL, 
      G_CALLBACK(save_session_callback) },
    { "SaveSessionAs", GTK_STOCK_SAVE_AS, N_("Save session _as..."), NULL, NULL, 
      G_CALLBACK(save_session_callback) },

    { "Databases", NULL, N_("_Databases"), NULL, NULL, NULL },
    { "NativeDB", GTK_STOCK_OPEN, N_("_Gretl native..."), "", NULL, G_CALLBACK(show_files) },
    { "RATSDB", GTK_STOCK_OPEN, N_("_RATS 4..."), "", NULL, G_CALLBACK(open_data) },
    { "PcGiveDB", GTK_STOCK_OPEN, N_("_PcGive..."), "", NULL, G_CALLBACK(open_data) },
    { "RemoteDB", GTK_STOCK_NETWORK, N_("On database _server..."), NULL, NULL, G_CALLBACK(show_files) },

    { "FunctionFiles", NULL, N_("_Function files"), NULL, NULL, NULL },
    { "LocalGfn", GTK_STOCK_OPEN, N_("On _local machine..."), "", NULL, G_CALLBACK(show_files) },
    { "RemoteGfn", GTK_STOCK_NETWORK, N_("On _server..."), NULL, NULL, G_CALLBACK(show_files) },
    { "EditGfn", GTK_STOCK_EDIT, N_("Edit package..."), NULL, NULL, G_CALLBACK(edit_package_callback) },
    { "NewGfn", GTK_STOCK_NEW, N_("_New package"), "", NULL, G_CALLBACK(new_package_callback) },

    { "Quit", GTK_STOCK_QUIT, N_("E_xit"), "<control>X", NULL,  G_CALLBACK(menu_exit_check)},

    /* Tools */
    { "Tools", NULL, N_("_Tools"), NULL, NULL, NULL },
    { "StatsTables", NULL, N_("_Statistical tables"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "PValues", NULL, N_("_P-value finder"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "DistGraphs", NULL, N_("_Distribution graphs"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "PlotCurve", NULL, N_("_Plot a curve"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "TestStats", NULL, N_("_Test statistic calculator"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "NonparamTests", NULL, N_("_Nonparametric tests"), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "SetSeed", NULL, N_("_Seed for random numbers"), NULL, NULL, G_CALLBACK(rand_seed_dialog) },
    { "CommandLog", NULL, N_("_Command log"), NULL, NULL, G_CALLBACK(view_command_log) },
    { "ShowConsole", NULL, N_("_Gretl console"), NULL, NULL, G_CALLBACK(gretl_console) },
    { "StartR", NULL, N_("Start GNU _R"), NULL, NULL, G_CALLBACK(start_R_callback) },
    { "NistTest", NULL, N_("_NIST test suite"), NULL, NULL, NULL },
    { "NistBasic", NULL, N_("_Basic"), NULL, NULL, G_CALLBACK(do_nistcheck) },
    { "NistVerbose", NULL, N_("_Verbose"), NULL, NULL, G_CALLBACK(do_nistcheck) },
    { "NistVVerbose", NULL, N_("V_ery verbose"), NULL, NULL, G_CALLBACK(do_nistcheck) },
    { "Preferences", NULL, N_("_Preferences"), NULL, NULL, NULL },    
    { "PrefsGeneral", GTK_STOCK_PREFERENCES, N_("_General..."), NULL, NULL, 
      G_CALLBACK(options_dialog_callback) },
    { "FixedFont", GTK_STOCK_SELECT_FONT, N_("_Fixed font..."), NULL, NULL, 
      G_CALLBACK(font_selector) },
    { "MenuFont", GTK_STOCK_SELECT_FONT, N_("_Menu font..."), NULL, NULL, 
      G_CALLBACK(font_selector) },

    /* Data */
    { "Data", NULL, N_("_Data"), NULL, NULL, NULL },
    { "DataSelectAll", NULL, N_("Select _all"), "<control>A", NULL, G_CALLBACK(mdata_select_all) },
    { "DefineList", NULL, N_("Define or edit _list..."), NULL, NULL, G_CALLBACK(gui_define_list) },
    { "SelectList", NULL, N_("_Set selection from list..."), NULL, NULL, G_CALLBACK(mdata_select_list) },
    { "DisplayValues", NULL, N_("_Display values"), NULL, NULL, G_CALLBACK(display_selected) },
    { "EditValues", NULL, N_("_Edit values"), NULL, NULL, G_CALLBACK(spreadsheet_edit) },
    { "AddObs", NULL, N_("_Add observations..."), NULL, NULL, G_CALLBACK(do_add_obs) },
    { "RemoveObs", NULL, N_("_Remove extra observations"), NULL, NULL, G_CALLBACK(do_remove_obs) },
    { "DataInfo", NULL, N_("_Dataset info"), NULL, NULL, G_CALLBACK(dataset_info) },
    { "DataMarkers", NULL, N_("_Observation markers..."), NULL, NULL, G_CALLBACK(markers_callback) },
    { "VarLabels", NULL, N_("_Variable labels..."), NULL, NULL, G_CALLBACK(labels_callback) },
    { "DataStructure", NULL, N_("Dataset _structure..."), NULL, NULL, G_CALLBACK(data_structure_dialog) },
    { "DataCompact", NULL, N_("_Compact data..."), NULL, NULL, G_CALLBACK(do_compact_data_set) },
    { "DataExpand", NULL, N_("_Expand data..."), NULL, NULL, G_CALLBACK(do_expand_data_set) },
    { "DataTranspose", NULL, N_("_Transpose data..."), NULL, NULL, G_CALLBACK(gui_transpose_data) },
    { "DataSort", NULL, N_("_Sort data..."), NULL, NULL, G_CALLBACK(gui_sort_data) },
    { "GSETMISS", NULL, N_("Set missing _value code..."), NULL, NULL, G_CALLBACK(gretl_callback) },
    { "VarFind", GTK_STOCK_FIND, N_("_Find variable..."), NULL, NULL, G_CALLBACK(listbox_find) },
    { "DataRefresh", NULL, N_("_Refresh window"), NULL, NULL, G_CALLBACK(refresh_data) },

    /* View */
    { "View", NULL, N_("_View"), NULL, NULL, NULL },
    { "IconView", NULL, N_("_Icon view"), NULL, NULL, G_CALLBACK(view_session) },
    { "Windows", NULL, N_("_Windows"), NULL, NULL, NULL },
    { "GraphVars", NULL, N_("_Graph specified vars"), NULL, NULL, NULL },
    { "TSPlot", NULL, N_("_Time series plot..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "ScatterPlot", NULL, N_("X-Y _scatter..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "ImpulsePlot", NULL, N_("X-Y with _impulses..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "FactorPlot", NULL, N_("X-Y with _factor separation..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "FrischPlot", NULL, N_("X-Y with _control..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "GR_BOX", NULL, N_("_Boxplots..."), NULL, NULL, G_CALLBACK(gretl_callback) },
    { "GR_FBOX", NULL, N_("_Factorized boxplot..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "GR_QQ", NULL, N_("_Q-Q plot..."), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "ThreeDPlot", NULL, N_("_3D plot..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "MultiPlots", NULL, N_("_Multiple graphs"), NULL, NULL, NULL },
    { "MultiXY", NULL, N_("X-Y _scatters..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "MultiTS", NULL, N_("_Time series..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "summary", NULL, N_("_Summary statistics"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "corr", NULL, N_("_Correlation matrix"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "xtab", NULL, N_("Cross _Tabulation"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "pca", NULL, N_("_Principal components"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "mahal", NULL, N_("_Mahalanobis distances"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "xcorrgm", NULL, N_("C_ross-correlogram"), NULL, NULL, G_CALLBACK(xcorrgm_callback) },

    /* Add */
    { "Add", NULL, N_("_Add"), NULL, NULL, NULL },
    { "logs", NULL, N_("_Logs of selected variables"), NULL, NULL, G_CALLBACK(logs_etc_callback) },
    { "square", NULL, N_("_Squares of selected variables"), NULL, NULL, G_CALLBACK(logs_etc_callback) },
    { "lags", NULL, N_("_Lags of selected variables"), NULL, NULL, G_CALLBACK(logs_etc_callback) },
    { "diff", NULL, N_("_First differences of selected variables"), NULL, NULL, 
      G_CALLBACK(logs_etc_callback) },
    { "ldiff", NULL, N_("_Log differences of selected variables"), NULL, NULL, 
      G_CALLBACK(logs_etc_callback) },
    { "sdiff", NULL, N_("_Seasonal differences of selected variables"), NULL, NULL, 
      G_CALLBACK(logs_etc_callback) },
    { "AddIndex", NULL, N_("_Index variable"), NULL, NULL, G_CALLBACK(add_index) },
    { "AddTime", NULL, N_("_Time trend"), NULL, NULL, G_CALLBACK(add_index) },
    { "AddRandom", NULL, N_("_Random variable..."), NULL, NULL, G_CALLBACK(stats_calculator) },
    { "PeriodDums", NULL, N_("_Periodic dummies"), NULL, NULL, G_CALLBACK(add_dummies) },
    { "UnitDums", NULL, N_("_Unit dummies"), NULL, NULL, G_CALLBACK(add_dummies) },
    { "TimeDums", NULL, N_("_Time dummies"), NULL, NULL, G_CALLBACK(add_dummies) },
    { "dummify", NULL, N_("Dummies for selected _discrete variables"), NULL, NULL, 
      G_CALLBACK(logs_etc_callback) },
    { "NewMatrix", NULL, N_("_Define matrix..."), NULL, NULL, G_CALLBACK(gui_new_matrix) },

    /* Sample */
    { "Sample", NULL, N_("_Sample"), NULL, NULL, NULL },
    { "SMPL", NULL, N_("_Set range..."), NULL, NULL, G_CALLBACK(sample_range_dialog) },
    { "FullRange", NULL, N_("_Restore full range"), NULL, NULL, G_CALLBACK(restore_sample_callback) },
    { "ShowSample", NULL, N_("_Show status"), NULL, NULL, G_CALLBACK(show_sample_callback) },
    { "SMPLBOOL", NULL, N_("_Restrict, based on criterion..."), NULL, NULL, 
      G_CALLBACK(sample_restrict_dialog) },
    { "SMPLRAND", NULL, N_("R_andom sub-sample..."), NULL, NULL, G_CALLBACK(sample_range_dialog) },
    { "SampleWReplace", NULL, N_("_Resample with replacement..."), NULL, NULL, 
      G_CALLBACK(gui_resample_data) },
    { "DropMissing", NULL, N_("Drop all obs with _missing values"), NULL, NULL, 
      G_CALLBACK(drop_all_missing) },
    { "CountMissing", NULL, N_("_Count missing values"), NULL, NULL, G_CALLBACK(count_missing) },

    /* Variable */
    { "Variable", NULL, N_("_Variable"), NULL, NULL, NULL },
    { "VarDisplay", NULL, N_("_Display values"), NULL, NULL, G_CALLBACK(display_var) },
    { "VarSummary", NULL, N_("_Summary statistics"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "normtest", NULL, N_("_Normality test"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "FreqDist", NULL, N_("_Frequency distribution..."), NULL, NULL, G_CALLBACK(do_freq_dist) },
    { "Density", NULL, N_("Estimated _density plot..."), NULL, NULL, G_CALLBACK(do_kernel) },
    { "qqplot", NULL, N_("Normal _Q-Q plot..."), NULL, NULL, G_CALLBACK(do_qqplot) },
    { "Gini", NULL, N_("_Gini coefficient"), NULL, NULL, G_CALLBACK(do_gini) },
    { "rmplot", NULL, N_("_Range-mean graph"), NULL, NULL, G_CALLBACK(do_range_mean) },
    { "VarTSPlot", NULL, N_("_Time series plot"), NULL, NULL, G_CALLBACK(ts_plot_callback) },
    { "PanPlot", NULL, N_("_Panel plot..."), NULL, NULL, G_CALLBACK(ts_plot_callback) },
    { "URTests", NULL, N_("_Unit root tests"), NULL, NULL, NULL },
    { "adf", NULL, N_("_Augmented Dickey-Fuller test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "dfgls", NULL, N_("ADF-GLS test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "kpss", NULL, N_("_KPSS test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "levinlin", NULL, N_("_Levin-Lin-Chu test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "fractint", NULL, N_("_Fractional integration"), NULL, NULL, G_CALLBACK(do_fractint) },
    { "corrgm", NULL, N_("_Correlogram"), NULL, NULL, G_CALLBACK(do_corrgm) },
    { "pergm", NULL, N_("_Periodogram"), NULL, NULL, G_CALLBACK(do_pergm) },
    { "Filter", NULL, N_("_Filter"), NULL, NULL, NULL },
    { "FilterSMA", NULL, N_("_Simple moving average"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterEMA", NULL, N_("_Exponential moving average"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterHP", NULL, N_("_Hodrick-Prescott"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterBK", NULL, N_("_Baxter-King"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterBW", NULL, N_("_Butterworth"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterPoly", NULL, N_("_Polynomial trend"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterFD", NULL, N_("_Fractional difference"), NULL, NULL, G_CALLBACK(filter_callback) },
    
#ifdef HAVE_X12A
    { "X12A", NULL, N_("_X-12-ARIMA analysis"), NULL, NULL, G_CALLBACK(do_tramo_x12a) },    
#endif
#ifdef HAVE_TRAMO
    { "Tramo", NULL, N_("_TRAMO analysis"), NULL, NULL, G_CALLBACK(do_tramo_x12a) },  
#endif
    { "Hurst", NULL, N_("_Hurst exponent"), NULL, NULL, G_CALLBACK(do_hurst) },     
    { "EditAttrs", NULL, N_("_Edit attributes"), NULL, NULL, G_CALLBACK(varinfo_callback) },     
    { "VSETMISS", NULL, N_("Set missing _value code..."), NULL, NULL, G_CALLBACK(gretl_callback) },     
    { "GENR", NULL, N_("Define _new variable..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 

    /* Model */
    { "Model", NULL, N_("_Model"), NULL, NULL, NULL },
    { "ols", NULL, N_("_Ordinary Least Squares..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "ivreg", NULL, N_("_Instrumental variables"), NULL, NULL, NULL },
    { "tsls", NULL, N_("_Two-Stage Least Squares..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "iv-liml", NULL, N_("_LIML..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "iv-gmm", NULL, N_("_GMM..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "LinearModels", NULL, N_("Other _linear models"), NULL, NULL, NULL },
    { "wls", NULL, N_("_Weighted Least Squares..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "hsk", NULL, N_("H_eteroskedasticity corrected..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "mpols", NULL, N_("High _precision OLS..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "anova", NULL, N_("ANOVA..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "TSModels", NULL, N_("_Time series"), NULL, NULL, NULL },
    { "CORC", NULL, N_("_Cochrane-Orcutt..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "HILU", NULL, N_("_Hildreth-Lu..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "PWE", NULL, N_("_Prais-Winsten..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "ar", NULL, N_("_Autoregressive estimation..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "arima", NULL, N_("ARI_MA..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "garch", NULL, N_("_GARCH..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "var", NULL, N_("_Vector Autoregression..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "VLAGSEL", NULL, N_("VAR _lag selection..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "vecm", NULL, N_("V_ECM..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "CointMenu", NULL, N_("_Cointegration test"), NULL, NULL, NULL },
    { "coint", NULL, N_("_Engle-Granger..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "coint2", NULL, N_("_Johansen..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "PanelModels", NULL, N_("_Panel"), NULL, NULL, NULL },
    { "panel", NULL, N_("_Fixed or random effects..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "PANEL_WLS", NULL, N_("_Weighted least squares..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "PANEL_B", NULL, N_("_Between model..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "dpanel", NULL, N_("_Dynamic panel model..."), NULL, NULL, G_CALLBACK(model_callback) },

    { "NonlinearModels", NULL, N_("_Nonlinear models"), NULL, NULL, NULL },
    { "logit", NULL, N_("_Logit"), NULL, NULL, NULL }, 
    { "blogit", NULL, N_("_Binary..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "ologit", NULL, N_("_Ordered..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "mlogit", NULL, N_("_Multinomial..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "probit", NULL, N_("_Probit"), NULL, NULL, NULL }, 
    { "bprobit", NULL, N_("_Binary..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "oprobit", NULL, N_("_Ordered..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "biprobit", NULL, N_("Bi_variate..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "tobit", NULL, N_("To_bit..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "heckit", NULL, N_("_Heckit..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "countmod", NULL, N_("_Count data..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "duration", NULL, N_("_Duration data..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "logistic", NULL, N_("Lo_gistic..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "nls", NULL, N_("_Nonlinear Least Squares..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "intreg", NULL, N_("_Interval regression..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "RobustModels", NULL, N_("_Robust estimation"), NULL, NULL, NULL },
    { "lad", NULL, N_("Least _Absolute Deviation..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "quantreg", NULL, N_("_Quantile regression..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "spearman", NULL, N_("_Rank correlation..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "loess", NULL, N_("_Loess..."), NULL, NULL, G_CALLBACK(selector_callback) }, 
    { "nadarwat", NULL, N_("_Nadaraya-Watson..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "mle", NULL, N_("_Maximum likelihood..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "gmm", NULL, N_("_GMM..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "system", NULL, N_("_Simultaneous equations..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 

    /* Help */
    { "Help", NULL, N_("_Help"), NULL, NULL, NULL },
    { "TextCmdRef", GRETL_STOCK_BOOK, N_("_Command reference"), NULL, NULL, G_CALLBACK(plain_text_cmdref) },
    { "FuncRef",   NULL, N_("_Function reference"), NULL, NULL, G_CALLBACK(genr_funcs_ref) },
    { "UserGuide", GRETL_STOCK_PDF, N_("_User's guide"), NULL, NULL, G_CALLBACK(display_pdf_help) },
    { "PDFCmdRef", GRETL_STOCK_PDF, N_("_Command reference"), NULL, NULL, G_CALLBACK(display_pdf_help) },
    { "UpdateCheck", GTK_STOCK_NETWORK, N_("Check for _updates"), NULL, NULL, G_CALLBACK(update_query) },
    { "SFAddons", NULL, N_("Check for _addons"), NULL, NULL, G_CALLBACK(show_files) },
    { "About", GTK_STOCK_ABOUT, N_("_About gretl"), NULL, NULL, G_CALLBACK(about_dialog) }
};

static void add_conditional_items (windata_t *vwin)
{
    GtkUIManager *ui = vwin->ui;
    int add_appfont = 1;

#ifdef G_OS_WIN32
    if (use_wimp) {
	add_appfont = 0;
    }
#endif

    if (add_appfont) {
	gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			      "/menubar/Tools/Preferences",
			      N_("_Menu font..."),
			      "MenuFont",
			      GTK_UI_MANAGER_MENUITEM, 
			      FALSE);
    }

#ifdef HAVE_X12A
    gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			  "/menubar/Variable/X12A",
			  N_("_X-12-ARIMA analysis"),
			  "X12A",
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);
#endif

#ifdef HAVE_TRAMO
    gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			  "/menubar/Variable/Tramo",
			  N_("_TRAMO analysis"),
			  "Tramo",
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);
#endif

    if (ox_support) {
	gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			      "/menubar/File/ScriptFiles/NewScript/OxScript",
			      N_("Ox program"),
			      "OxScript",
			      GTK_UI_MANAGER_MENUITEM, 
			      FALSE);
    }

    maybe_add_packages_to_menus(vwin);
}

/* retrieve the XML description of the main window menus */

static gchar *get_main_ui (void)
{
    gchar *main_ui = NULL;
    gchar *fname;
    int err;

    fname = g_strdup_printf("%sui%cgretlmain.xml", gretl_home(), SLASH);
    err = gretl_file_get_contents(fname, &main_ui, NULL);
    g_free(fname);

    return err ? NULL : main_ui;
}

static int set_up_main_menu (void)
{
    GtkActionGroup *actions;
    gchar *main_ui = NULL;
    GError *error = NULL;
    int err = 0;

    main_ui = get_main_ui();
    if (main_ui == NULL) {
	return 1;
    }

    mdata->ui = gtk_ui_manager_new();
    actions = gtk_action_group_new("Actions");
    gtk_action_group_set_translation_domain(actions, "gretl");
    gtk_action_group_add_actions(actions, main_entries, 
				 G_N_ELEMENTS(main_entries), mdata);

    gtk_ui_manager_insert_action_group(mdata->ui, actions, 0);
    g_object_unref(actions);

    gtk_window_add_accel_group(GTK_WINDOW(mdata->main), 
			       gtk_ui_manager_get_accel_group(mdata->ui));

    if (!gtk_ui_manager_add_ui_from_string(mdata->ui, main_ui, -1, &error)) {
	g_message("building menus failed: %s", error->message);
	g_error_free(error);
	err = 1;
    } else {
	add_conditional_items(mdata);
	mdata->mbar = gtk_ui_manager_get_widget(mdata->ui, "/menubar");
	if (mdata->mbar == NULL) {
	    fprintf(stderr, "/menubar widget is NULL!\n");
	    err = 1;
	} else {
	    GtkWidget *dataitem;

	    dataitem = gtk_ui_manager_get_widget(mdata->ui, "/menubar/Data");
	    if (dataitem == NULL) {
		fprintf(stderr, "/menubar/Data widget is NULL!\n");
		err = 1;
	    } else {
		g_signal_connect(G_OBJECT(dataitem), "activate",
				 G_CALLBACK(check_var_labels_state), mdata);
		g_free(main_ui);
	    }
	}
    }

    return err;
}

void raise_main_window (void)
{
    gtk_window_present(GTK_WINDOW(mdata->main));
}

static void gretl_window_raise (GtkAction *action, gpointer data)
{
    const gchar *s = gtk_action_get_name(action);
    unsigned long lptr = strtoul(s, NULL, 16);
    GtkWidget *wptr = (GtkWidget *) lptr;

    if (GTK_IS_WINDOW(wptr)) {
	gtk_window_present(GTK_WINDOW(wptr));
    }
}

static const gchar *window_list_icon (int role)
{
    const gchar *id = NULL;

    if (role == VIEW_MODEL) {
	id = GRETL_STOCK_MODEL;
    } else if (role == CONSOLE) {
	id = GRETL_STOCK_CONSOLE;
    } else if (role >= EDIT_HEADER && 
	       role < EDIT_MAX) {
	id = GTK_STOCK_EDIT;
    } else if (role == GNUPLOT) {
	id = GRETL_STOCK_SCATTER;
    } else if (BROWSER_ROLE(role)) {
	id = GTK_STOCK_INDEX;
    } else if (HELP_ROLE(role)) {
	id = GRETL_STOCK_BOOK;
    } else if (role == STAT_TABLE) {
	id = GRETL_STOCK_CALC;
    } else if (role == VIEW_SCRIPT) {
	id = GTK_STOCK_EXECUTE;
    } else if (role == OPEN_SESSION) {
	id = GRETL_STOCK_ICONS;
    }

    return id;
}

static int n_listed_windows;
static GtkActionGroup *window_list;

int get_n_listed_windows (void)
{
    return n_listed_windows;
}

static const gchar *get_window_title (GtkWidget *w)
{
    const gchar *s = NULL;

    if (GTK_IS_WINDOW(w)) {
	s = gtk_window_get_title(GTK_WINDOW(w));

	if (s != NULL && !strncmp(s, "gretl", 5)) {
	    s += 5;
	    s += strspn(s, " ");
	    if (*s == ':') {
		s++;
		s += strspn(s, " ");
	    }
	}
    }

    return s;
}

static void destroy_window_action (GtkWidget *w, GtkActionGroup *group)
{
    guint merge_id = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "merge_id"));

    if (merge_id > 0) {
	GtkAction *action;
	gchar *aname;

	/* remove the ui element */
	gtk_ui_manager_remove_ui(mdata->ui, merge_id);
	n_listed_windows--;
	if (n_listed_windows == 0) {
	    window_list_state(FALSE);
	}

	/* and also the underlying action */
	aname = g_strdup_printf("%p", (void *) w);
	action = gtk_action_group_get_action(group, aname);
	if (action != NULL) {
	    gtk_action_group_remove_action(group, action);
	}
	g_free(aname);
    }    
}

static gchar *get_parent_path (GtkWidget *w)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(w), "vwin");
    gchar *path = NULL;

    if (vwin != NULL && vwin->gretl_parent != NULL) {
	/* got a viewer with a 'parent' */
	path = g_strdup_printf("/menubar/View/Windows/%p", 
			       vwin->gretl_parent->main);
    } 
    
    return path;
}

void add_window_list_item (GtkWidget *w, int role)
{
    GtkActionEntry entry = { 
	NULL, NULL, NULL, NULL, NULL, G_CALLBACK(gretl_window_raise) 
    };
    const gchar *s = get_window_title(w);
    gchar *aname;
    guint merge_id;

    if (window_list == NULL) {
	/* create the window_list group */
	window_list = gtk_action_group_new("WindowList");
	gtk_ui_manager_insert_action_group(mdata->ui, window_list, 0);
    }

    /* set up an action entry */
    aname = g_strdup_printf("%p", (void *) w);
    entry.name = aname;
    entry.stock_id = window_list_icon(role);
    entry.label = s;

    /* add new action entry to group, and to ui */
    gtk_action_group_add_actions(window_list, &entry, 1, mdata);
    merge_id = gtk_ui_manager_new_merge_id(mdata->ui);

    if (merge_id > 0) {
	gchar *apath = get_parent_path(w);

	gtk_ui_manager_add_ui(mdata->ui, merge_id, 
			      (apath != NULL)? apath : "/menubar/View/Windows",
			      aname, aname,
			      GTK_UI_MANAGER_MENUITEM, FALSE);
	g_object_set_data(G_OBJECT(w), "merge_id", GINT_TO_POINTER(merge_id));
	g_signal_connect(G_OBJECT(w), "destroy", 
			 G_CALLBACK(destroy_window_action), 
			 window_list);

	window_list_state(TRUE);
	n_listed_windows++;
	g_free(apath);
    }

    g_free(aname);
}

static gint sort_window_items (gconstpointer a, gconstpointer b)
{
    const gchar *aname = gtk_action_get_name((GtkAction *) a);
    const gchar *bname = gtk_action_get_name((GtkAction *) b);
    GtkWidget *wa = (GtkWidget *) strtol(aname, NULL, 16);
    GtkWidget *wb = (GtkWidget *) strtol(bname, NULL, 16);
    windata_t *va = g_object_get_data(G_OBJECT(wa), "vwin");
    windata_t *vb = g_object_get_data(G_OBJECT(wb), "vwin");
    gint ret = 0;

    if (va != NULL && vb != NULL) {
	if (va->gretl_parent == vb) {
	    ret = 1;
	} else if (vb->gretl_parent == va) {
	    ret = -1;
	}
    }
    
    return ret;
}

void window_list_popup (GtkWidget *src)
{
    static GtkWidget *menu;    

    if (menu != NULL) {
	gtk_widget_destroy(menu);
	menu = NULL;
    }    

    if (n_listed_windows > 0) {
	GList *wlist = gtk_action_group_list_actions(window_list);
	GList *list = g_list_sort(wlist, sort_window_items);
	GtkAction *a;
	GtkWidget *w;
	int n_items = 0;

	menu = gtk_menu_new();

	while (list) {
	    if (n_items == 0) {
		/* at top: item for main window */
		w = gtk_menu_item_new_with_label(_("Main window"));
		g_signal_connect(G_OBJECT(w), "activate",
				 G_CALLBACK(raise_main_window), NULL);
		gtk_widget_show(w);
		gtk_menu_shell_append(GTK_MENU_SHELL(menu), w);
		n_items++;
	    }
	    a = (GtkAction *) list->data;
	    w = gtk_action_create_menu_item(a);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), w);
	    n_items++;
	    list = list->next;
	}

	g_list_free(wlist);

	if (n_items > 0) {
	    gtk_menu_popup(GTK_MENU(menu), NULL, NULL, NULL, NULL,
			   0, gtk_get_current_event_time());
	} else {
	    gtk_widget_destroy(menu);
	    menu = NULL;
	}
    }
}

/* on exiting, check for any editing windows with unsaved
   changes, and if we find any give the user a chance to
   save the changes, or to cancel the exit
*/

gboolean window_list_exit_check (void)
{
    gboolean ret = FALSE;

    if (n_listed_windows > 0) {
	GList *list = gtk_action_group_list_actions(window_list);
	const gchar *name;
	GtkAction *action;
	windata_t *vwin;
	GtkWidget *w;

	while (list) {
	    action = (GtkAction *) list->data;
	    name = gtk_action_get_name(action);
	    w = (GtkWidget *) strtol(name, NULL, 16);
	    vwin = g_object_get_data(G_OBJECT(w), "vwin");

	    if (vwin != NULL) {
		if (vwin_is_editing(vwin) && vwin_content_changed(vwin)) {
		    gtk_window_present(GTK_WINDOW(vwin->main));
		    ret = query_save_text(NULL, NULL, vwin);
		}
	    }
	    list = list->next;
	}

	g_list_free(list);
    }

    return ret;
}

int gui_restore_sample (DATASET *dset)
{
    int err = 0;

    if (dset != NULL && dset->Z != NULL) {
	err = restore_full_sample(dset, NULL);
	if (err) {
	    gui_errmsg(err);
	} else {
	    restore_sample_state(FALSE);
	}
    }

    return err;
}

static void restore_sample_callback (void)
{
    int err = gui_restore_sample(dataset); 

    if (!err) {
	set_sample_label(dataset);    
	lib_command_strcpy("smpl --full");
	record_command_verbatim();
    }
}

static void show_sample_callback (void)
{
    gchar *title;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    print_sample_status(dataset, prn);
    title = g_strdup_printf("gretl: %s", _("sample status"));
    view_buffer(prn, 78, 360, title, PRINT, NULL);
    g_free(title);
}

static void start_R_callback (void)
{
    start_R(NULL, 1, 1);
}

#ifndef G_OS_WIN32

int gretl_fork (const char *progvar, const char *arg)
{
    const char *prog = NULL;
    gchar *argv[3];
    GError *err = NULL;
    gboolean run;

#ifdef OSX_BUILD
    if (!strcmp(progvar, "calculator")) {
	prog = calculator;
    } else if (!strcmp(progvar, "viewdvi")) {
	prog = viewdvi;
    } 
#else
    if (!strcmp(progvar, "Browser")) {
	prog = Browser;
    } else if (!strcmp(progvar, "calculator")) {
	prog = calculator;
    } else if (!strcmp(progvar, "viewdvi")) {
	prog = viewdvi;
    } else if (!strcmp(progvar, "viewpdf")) {
	prog = viewpdf;
    } else if (!strcmp(progvar, "viewps")) {
	prog = viewps;
    }    
#endif

    if (prog == NULL) {
	errbox("Internal error: variable %s is undefined", progvar);
	return 1;
    }
    
    argv[0] = g_strdup(prog);
    if (arg != NULL) {
	argv[1] = g_strdup(arg);
	argv[2] = NULL;
    } else {
	argv[1] = NULL;
    }
    
    run = g_spawn_async(NULL, argv, NULL, G_SPAWN_SEARCH_PATH, 
			NULL, NULL, NULL, &err);

    if (err != NULL) {
	errbox(err->message);
	if (err->domain == G_SPAWN_ERROR &&
	    err->code == G_SPAWN_ERROR_NOENT) {
	    options_dialog(TAB_PROGS, progvar, mdata->main);
	}
	g_error_free(err);
    }

    g_free(argv[0]);
    g_free(argv[1]);

    return !run;
}

#endif	

/* Icon handling for X11 */

#ifndef G_OS_WIN32

void set_wm_icon (GtkWidget *w)
{
    GdkPixbuf *icon = gdk_pixbuf_new_from_xpm_data(gretl_xpm);

    if (icon != NULL) {
	gtk_window_set_icon(GTK_WINDOW(w), icon);
    } else {
	fprintf(stderr, "Couldn't create icon for gretl_xpm\n");
    }
}

#endif

static int has_db_suffix (const char *fname)
{
    return has_suffix(fname, ".bin") ||
	has_suffix(fname, ".rat") ||
	has_suffix(fname, ".bn7");
}

/* Drag 'n' drop: respond to data dropped into the main window */

static void  
mdata_handle_drag  (GtkWidget *widget,
		    GdkDragContext *context,
		    gint x,
		    gint y,
		    GtkSelectionData *data,
		    guint info,
		    guint time,
		    gpointer p)
{
    GdkAtom type = gtk_selection_data_get_data_type(data);
    const guchar *seldata = NULL;
    gchar *dfname;
    char tmp[MAXLEN];
    int pos, skip = 5;
    int ftype = 0;

    if (data != NULL) {
	seldata = gtk_selection_data_get_data(data);
    }

    /* handle drag of pointer from database window */
    if (info == GRETL_DBSERIES_PTR && data != NULL && 
	type == GDK_SELECTION_TYPE_INTEGER) {
	import_db_series(*(void **) seldata);
	return;
    }

    if (info != GRETL_FILENAME) {
	return;
    }

    /* ignore the wrong sort of data */
    if (data == NULL || (dfname = (gchar *) seldata) == NULL || 
	strlen(dfname) <= 5 || strncmp(dfname, "file:", 5)) {
	return;
    }

    if (strncmp(dfname, "file://", 7) == 0) skip = 7;
#ifdef G_OS_WIN32
    if (strncmp(dfname, "file:///", 8) == 0) skip = 8;
#endif

    /* there may be multiple files: we ignore all but the first */
    *tmp = 0;
    if ((pos = charpos('\r', dfname)) > 0 || 
	(pos = charpos('\n', dfname) > 0)) {
	strncat(tmp, dfname + skip, pos - skip);
    } else {
	strcat(tmp, dfname + skip);
    }

    /* handle spaces and such */
    unescape_url(tmp);

#ifdef G_OS_WIN32
    filename_to_win32(tryfile, tmp);
#else
    strcpy(tryfile, tmp);
#endif

    if (has_db_suffix(tryfile)) {
	open_named_db_index(tryfile);
    } else if (has_suffix(tryfile, ".gretl") && gretl_is_pkzip_file(tryfile)) {
	verify_open_session();
    } else if ((ftype = script_type(tryfile))) {
	do_open_script(ftype);
    } else {
	verify_open_data(NULL, 0);
    }
}

/* the callback for Save Data (Ctrl-S) in main window */

static void auto_store (void)
{
    if (data_status & SESSION_DATA) {
	/* the data file is embedded in a session file */
	save_session_dataset();
    } else {
	/* ensure there's no stale selection around */
	set_selector_storelist(NULL);

	if ((data_status & USER_DATA) && has_suffix(datafile, ".gdt")) {
	    /* bypass filename selection */
	    do_store(datafile, SAVE_DATA, NULL);
	} else {
	    file_selector(SAVE_DATA, FSEL_DATA_NONE, NULL);
	}
    }	
}

int mdata_selection_count (void)
{
    return vwin_selection_count(mdata, NULL);
}

int mdata_active_var (void)
{
    int selcount, v = 0;

    selcount = vwin_selection_count(mdata, &v);

    if (selcount == 1 && v != 0) {
	mdata->active_var = v;
    }

    return mdata->active_var;
}

static gboolean 
main_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer data)
{
    GdkModifierType mods = widget_get_pointer_mask(w);

    if (RIGHT_CLICK(mods)) {
	/* ignore all but right-clicks */
	int selcount = vwin_selection_count(mdata, NULL);

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
	    g_signal_connect(G_OBJECT(mdata->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed), 
			     &mdata->popup);
	}

	return TRUE;
    }

    return FALSE;
}

static int script_stopper (int set)
{
    static int stop;
    int ret = 0;

    if (set) {
	/* set the stop signal */
	stop = 1;
    } else if (stop) {
	warnbox(_("Execution aborted by request"));
	ret = 1;
	stop = 0;
    }

    return ret;
}

static int gui_query_stop (void)
{
    return script_stopper(0);
}

void do_stop_script (GtkWidget *w, windata_t *vwin)
{
    script_stopper(1);
}
