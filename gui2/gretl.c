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
#include "version.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "libset.h"

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

#include <dirent.h>

#ifndef G_OS_WIN32
# include <unistd.h>
# include <sys/types.h>
# include "../pixmaps/gretl.xpm"  /* program icon for X */
#else
# include <windows.h>
# include "gretlwin32.h"
#endif

#if GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6
# define OLD_GTK
#endif

/* update.c */
extern int silent_update_query (void);
extern int update_query (void); 

/* functions private to gretl.c */
static GtkWidget *make_main_window (void);

static gboolean main_popup_handler (GtkWidget *w, GdkEventButton *event,
				    gpointer data);
static int set_up_main_menu (void);
static void start_R_callback (void);
static void auto_store (void);
static void restore_sample_callback (void);
static void mdata_select_all (void);

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

static char *optrun, *optdb, *optwebdb;
static int opteng, optdump, optver;
#ifndef OLD_GTK
static int optswitch;
#endif
#ifdef G_OS_WIN32
static int optdebug;
#endif

static gchar *param_msg = 
    N_("\nYou may supply the name of a data file on the command line");

#ifdef OLD_GTK

struct start_opts {
    const char *longopt;
    char shortopt;
    const char *optstr;
    const char *filetype;
};

static struct start_opts options[] = {
    { "run",     'r', N_("open a script file on startup"), "SCRIPT" },
    { "db",      'd', N_("open a database on startup"),    "DATABASE"},
    { "webdb",   'w', N_("open a remote (web) database on startup"), "REMOTE_DB" },
    { "english", 'e', N_("force use of English"), NULL },
    { "dump",    'c', N_("dump gretl configuration to file"), NULL },
    { "version", 'v', N_("print version information"), NULL }, 
    { NULL, 0, NULL, NULL },
};

static void gui_help (void) 
{
    int i;

    putchar('\n');

    for (i=0; options[i].longopt != NULL; i++) {
	if (options[i].filetype != NULL) {
	    printf("gretl --%s %s : %s\n", options[i].longopt, options[i].filetype, 
		   _(options[i].optstr));
	} else {
	    printf("gretl --%s : %s\n", options[i].longopt, _(options[i].optstr));
	}
    }

    printf("%s\n", _(param_msg));
}

static void old_gretl_init (int *pargc, char ***pargv, char *filearg)
{
    static char fname[FILENAME_MAX];
    gretlopt opt = 0;
    int err;

    err = parseopt(pargc, pargv, &opt, fname);

    if (opt & OPT_ENGLISH) {
	opteng = 1;
    } 

    if (err || (opt & (OPT_HELP | OPT_VERSION))) {
	int ecode = (err)? EXIT_FAILURE : EXIT_SUCCESS;

	gui_logo(NULL);
	if (opt != OPT_VERSION) {
	    gui_help();
	}
	exit(ecode);
    } else if (opt & OPT_DUMP) {
	optdump = 1;
    } else if (opt & OPT_RUNIT) {
	optrun = fname;
    } else if (opt & OPT_DBOPEN) {
	optdb = fname;
    } else if (opt & OPT_WEBDB) {
	optwebdb = fname;
    } else if (*fname != '\0') {
	/* got a data file argument? */
	*pargc += 1;
	strncat(filearg, fname, MAXLEN - 1);
    }
}

#else

static GOptionEntry options[] = {
    { "run", 'r', 0, G_OPTION_ARG_FILENAME, &optrun, 
      N_("open a script file on startup"), "SCRIPT" },
    { "db", 'd', 0, G_OPTION_ARG_STRING, &optdb, 
      N_("open a database on startup"), "DATABASE" },
    { "webdb", 'w', 0, G_OPTION_ARG_STRING, &optwebdb, 
      N_("open a remote (web) database on startup"), "REMOTE_DB" },
    { "english", 'e', 0, G_OPTION_ARG_NONE, &opteng, 
      N_("force use of English"), NULL },
    { "dump", 'c', 0, G_OPTION_ARG_NONE, &optdump, 
      N_("dump gretl configuration to file"), NULL },
#ifdef G_OS_WIN32
    { "debug", 'b', 0, G_OPTION_ARG_NONE, &optdebug, 
      N_("send debugging info to console"), NULL }, 
#endif
    { "version", 'v', 0, G_OPTION_ARG_NONE, &optver, 
      N_("print version information"), NULL }, 
    { "switch", 's', 0, G_OPTION_ARG_INT, &optswitch,
      N_("pass integer value to script"), "value" },
    { NULL, '\0', 0, 0, NULL, NULL, NULL },
};

#endif

windata_t *mdata;
DATAINFO *datainfo;

char scriptfile[MAXLEN];
char tryfile[MAXLEN];

PATHS paths;                /* useful paths */
double **Z;                 /* data set */
MODEL **models;             /* gretl models structs */

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
char calculator[MAXSTR] = "xcalc";
char latex[MAXSTR] = "pdflatex";
char viewdvi[MAXSTR] = "xdvi";
char Rcommand[MAXSTR] = "xterm -e R";
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
    varinfo_dialog(mdata->active_var, 1);
}

static void options_dialog_callback (void)
{
    options_dialog(0, NULL, mdata->main);
}

static void new_function_pkg_callback (GtkAction *action, gpointer p)
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
    send_file(paths.datfile);
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

    if (addpath(tryfile, &paths, 1) == NULL) {
	fprintf(stderr, I_("Couldn't find script '%s'\n"), tryfile);
	exit(EXIT_FAILURE);
    } else {
	fprintf(stderr, I_("%s found\n"), tryfile);
	set_currdir_from_filename(tryfile);
    }
}

static int foreign_script_type (const char *fname)
{
    if (has_suffix(fname, ".R")) {
	return EDIT_R;
    } else if (has_suffix(fname, ".plt") ||
	       has_suffix(fname, ".gp")) {
	return EDIT_GP;
    } else if (ox_support && has_suffix(fname, ".ox")) {
	return EDIT_OX;
    } else {
	return 0;
    }
}

static void fix_dbname (char *db)
{
    FILE *fp = NULL;

    if (strstr(db, ".bin") == NULL &&
	strstr(db, ".rat") == NULL &&
	strstr(db, ".RAT") == NULL &&
	strstr(db, ".bn7") == NULL) {
	strcat(db, ".bin");
    }

    fp = gretl_fopen(db, "rb");

    if (fp == NULL && strstr(db, paths.binbase) == NULL) {
	char tmp[MAXLEN];

	strcpy(tmp, db);
	build_path(db, paths.binbase, tmp, NULL);
    }

    if (fp != NULL) {
	fclose(fp);
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
    char gretldir[MAXSTR], localedir[MAXSTR];
    char *loc;

    if (read_reg_val(HKEY_LOCAL_MACHINE, "gretl", "gretldir", gretldir)) {
	return;
    }

    build_path(localedir, gretldir, "locale", NULL);
    loc = setlocale(LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#elif defined(OSX_BUILD) 

static void real_nls_init (void)
{
    char *gretlhome = getenv("GRETL_HOME");
    char localedir[MAXSTR];
    char *p, *loc;

    if (gretlhome == NULL) {
	return;
    }

    strcpy(localedir, gretlhome);
    p = strstr(localedir, "share/gretl");
    if (p != NULL) {
	strcpy(p, "share/locale");
    }

    /* FIXME GUI language choice? */
    loc = setlocale(LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain(PACKAGE, localedir);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#else /* regular *nix treatment */

static void real_nls_init (void)
{
    char *loc;

    loc = setlocale(LC_ALL, "");
    set_gretl_charset(loc);
    bindtextdomain(PACKAGE, LOCALEDIR);
    textdomain(PACKAGE);
    bind_textdomain_codeset(PACKAGE, "UTF-8");
}

#endif /* NLS init variants */

void nls_init (void)
{
    char *mylang = getenv("GRETL_LANG");

    if (mylang != NULL) {
	if (!g_ascii_strcasecmp(mylang, "english") ||
	    !g_ascii_strcasecmp(mylang, "C")) return;
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

#endif

static int have_data (void)
{
    return datainfo != NULL && datainfo->v > 0;
}

int main (int argc, char **argv)
{
#ifdef G_OS_WIN32
    char *callname = argv[0];
#endif
    int ftype = 0;
    char dbname[MAXLEN];
    char filearg[MAXLEN];
#ifndef OLD_GTK
    GError *opterr = NULL;
#endif
    int err = 0;

    nls_init();

    *tryfile = '\0';
    *scriptfile = '\0';
    *paths.datfile = '\0';
    *dbname = '\0';
    *filearg = '\0';

#ifdef OLD_GTK 
    gtk_init(&argc, &argv);
    old_gretl_init(&argc, &argv, filearg);
#else
    gtk_init_with_args(&argc, &argv, param_msg, options, "gretl", &opterr);
    if (opterr != NULL) {
	g_print("%s\n", opterr->message);
	exit(EXIT_FAILURE);
    }
    if (optswitch) {
	set_script_switch(optswitch);
    }
#endif

    libgretl_init();
    gretl_set_paths(&paths, OPT_D | OPT_X); /* defaults, gui */

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
	strncat(dbname, optdb, MAXLEN - 1);
	fix_dbname(dbname);
    } else if (optwebdb != NULL) {
	strncat(dbname, optdb, MAXLEN - 1);
    } 

    if (opteng) {
	force_language(LANG_C);
	force_english_help();
    } 

    set_workdir_callback(gui_set_working_dir);

    /* allocate data information struct */
    datainfo = datainfo_new();
    if (datainfo == NULL) {
	noalloc();
    }

    /* allocate memory for models */
    models = allocate_working_models(3);
    if (models == NULL) {
	noalloc();
    } 

    library_command_init();

    helpfile_init();
    session_init();
    init_fileptrs();

    if (argc > 1) {
	/* Process what is presumably a filename argument
	   given on the command line (by now any options will
	   have been extracted from the argv array).
	*/
	PRN *prn; 

	prn = gretl_print_new(GRETL_PRINT_STDERR, &err);
	if (err) exit(EXIT_FAILURE);

	if (*filearg == '\0') {
	    /* not already registered */
	    strncat(filearg, argv[1], MAXLEN - 1);
	}

	*paths.datfile = '\0';

#ifdef G_OS_WIN32
	if (filename_to_win32(paths.datfile, filearg)) {
	    exit(EXIT_FAILURE);
	}
#else
	record_filearg(paths.datfile, filearg);
#endif

	/* keep a copy of input filename */
	strcpy(tryfile, paths.datfile);

	ftype = detect_filetype(paths.datfile, &paths);

	switch (ftype) {
	case GRETL_NATIVE_DATA:
	    err = gretl_get_data(paths.datfile, &paths, &Z, datainfo, 
				 OPT_NONE, prn);
	    break;
	case GRETL_XML_DATA:
	    err = gretl_read_gdt(paths.datfile, &paths, &Z, datainfo, 
				 OPT_NONE, prn);
	    break;
	case GRETL_CSV:
	    err = import_csv(paths.datfile, &Z, datainfo, 
			     OPT_NONE, prn);
	    break;
	case GRETL_XLS:
	case GRETL_GNUMERIC:
	case GRETL_ODS:
	case GRETL_DTA:
	case GRETL_SAV:
	case GRETL_JMULTI:
	case GRETL_OCTAVE:
	case GRETL_WF1:
	    err = get_imported_data(paths.datfile, ftype, 0);
	    break;
	case GRETL_SCRIPT:
	case GRETL_SESSION:
	    get_runfile(paths.datfile);
	    *paths.datfile = '\0';
	    break;
	case GRETL_NATIVE_DB:
	case GRETL_RATS_DB:  
	case GRETL_PCGIVE_DB:
	    strcpy(dbname, paths.datfile);
	    *tryfile = '\0';
	    *paths.datfile = '\0';
	    fix_dbname(dbname);
	    optdb = dbname;
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
	    *paths.datfile = '\0';
	}

	if (ftype != GRETL_SCRIPT && err) {
	    errmsg(err, prn);
	    exit(EXIT_FAILURE);
	}

	gretl_print_destroy(prn);
    }

    /* create the GUI */
    gretl_tooltips_init();
    gretl_stock_icons_init();
    make_main_window();

    if (have_data()) {
	/* redundant? */
	set_sample_label(datainfo);
    }

    add_files_to_menus();

    session_menu_state(FALSE);
    restore_sample_state(FALSE);
    main_menubar_state(FALSE);

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
	} else if ((ftype = foreign_script_type(tryfile))) {
	    do_open_script(ftype);
	} else {
	    do_open_script(EDIT_SCRIPT);
	}
    }

    /* check for program updates? */
    if (updater) {
	silent_update_query(); 
    }

    /* try opening specified database */
    if (optdb != NULL) {
	open_named_db_index(dbname);
    } else if (optwebdb != NULL) {
	open_named_remote_db_index(dbname);
    }

    /* Enter the event loop */
    gtk_main();

    /* clean up before exiting */
    free_session();

    if (Z != NULL) {
	free_Z(Z, datainfo);
    }

    destroy_working_models(models, 3);

    library_command_free();
    libgretl_cleanup();

    if (data_status) {
	free_datainfo(datainfo);
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

/* if a keystroke would take us to row 0 (e.g. page up), countermand
   and go to row 1 instead
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

static gint catch_mdata_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods = widget_get_pointer_mask(w);
    int k = key->keyval;

#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    if (!(mods & GDK_MOD1_MASK)) {
	if (k == GDK_a) {
	    audio_render_window(vwin, AUDIO_LISTBOX);
	} else if (k == GDK_x) {
	    stop_talking();
	}
    }
#endif

    if (k == GDK_h || k == GDK_F1) {
	/* invoke help */
	plain_text_cmdref(NULL);
	return FALSE;
    }

    if (k == GDK_g) {
	/* invoke genr */
	genr_callback();
	return FALSE;
    }

    if (k == GDK_r) {
	refresh_data();
	return FALSE;
    }

    if (k == GDK_x && (mods & GDK_MOD1_MASK)) {
	/* invoke command minibuffer */
	minibuf_callback();
	return FALSE;
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
		varinfo_dialog(mdata->active_var, 1);
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
    } 

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

void populate_varlist (void)
{
    static gint check_connected;
    static gint click_connected;
    GtkTreeStore *store;
    GtkTreeSelection *select;
    GtkTreeIter iter;    
    char id[8];
    int pos = 0;
    int i;

    store = GTK_TREE_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox)));

    if (store != NULL) {
	/* record line position? */
	pos = get_line_pos(GTK_TREE_MODEL(store));
    }
    
    gtk_tree_store_clear(store);

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<datainfo->v; i++) {
	int pv = 0;

	if (var_is_hidden(datainfo, i)) {
	    continue;
	}
	if (i > 0 && (is_standard_lag(i, datainfo, &pv) ||
		      is_dummy_child(i, datainfo, &pv))) {
	    if (pv > 0) {
		GtkTreeIter child_iter, parent_iter;

		if (lagvar_get_parent_iter(pv, &parent_iter)) {
		    gtk_tree_store_insert_before(store, &child_iter, 
						 &parent_iter, NULL);
		    sprintf(id, "%d", i);
		    gtk_tree_store_set(store, &child_iter, 
				       0, id, 
				       1, datainfo->varname[i],
				       2, VARLABEL(datainfo, i),
				       -1);	
		}	
	    }
	}

	if (pv == 0) {
	    gtk_tree_store_append(store, &iter, NULL);
	    sprintf(id, "%d", i);
	    gtk_tree_store_set(store, &iter, 
			       0, id, 
			       1, datainfo->varname[i],
			       2, VARLABEL(datainfo, i),
			       -1);
	}
    } 

    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (pos == 0) {
	/* no saved position */
	pos = 1;
	gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
    } else {	
	/* try to return to previous position */
	GtkTreeIter last;

	i = 1;
	while (1) {
	    last = iter;
	    if (!gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter)) {
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
    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_select_iter(select, &iter);

    if (datainfo->v > 1) {
	GtkTreePath *path;

	sprintf(id, "%d", pos);
	path = gtk_tree_path_new_from_string(id);
	gtk_tree_view_set_cursor(GTK_TREE_VIEW(mdata->listbox), path, NULL, FALSE);
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
	g_signal_connect(G_OBJECT(mdata->listbox), "key-press-event",
			 G_CALLBACK(catch_mdata_key),
			 mdata);
	click_connected = 1;
    }

    variable_menu_state(TRUE);
}

void mdata_select_last_var (void)
{
    GtkTreeIter iter, last;
    GtkTreeStore *store;
    GtkTreeSelection *select;

    store = GTK_TREE_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(mdata->listbox)));
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    while (1) {
	last = iter;
	if (!gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter)) {
	    iter = last;
	    break;
	}
    }

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_unselect_all(select);
    gtk_tree_selection_select_iter(select, &iter);
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
    GtkSettings *settings;
    gchar *fontname = NULL;
    int fsize;
    float scale = 1.0;

    settings = gtk_settings_get_default();

    g_object_get(G_OBJECT(settings), "gtk-font-name", &fontname, NULL);

    if (fontname != NULL) {
	if (sscanf(fontname, "%*s %d", &fsize) == 1) {
	    if (fsize > 10 && fsize < 100) {
		scale = fsize / 10.0;
	    }
	}
	g_free(fontname);
    }

    return scale;
}

static gboolean 
mainwin_config (GtkWidget *w, GdkEventConfigure *event, gpointer p)
{
    mainwin_width = event->width;
    mainwin_height = event->height;

    gdk_window_get_root_origin(mdata->main->window, &main_x, &main_y);

    return FALSE;
}

/* scale up the main window if it seems to be too tiny in relation to
   the screen dimensions
*/

static void scale_main_window (void)
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

static GtkWidget *make_main_window (void) 
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
	scale_main_window();
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
    g_signal_connect_after(G_OBJECT(mdata->main), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif

    main_vbox = gtk_vbox_new(FALSE, 4);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 8);
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

    /* will hold the list of variables */
    box = gtk_vbox_new(FALSE, 0);

    align = gtk_alignment_new(0, 0, 0, 0);
    gtk_box_pack_start(GTK_BOX(box), align, FALSE, FALSE, 0);
    gtk_widget_show(align);
    gtk_container_add(GTK_CONTAINER(align), dlabel);
   
    vwin_add_list_box(mdata, GTK_BOX(box), 3, FALSE, types, titles, 1);

    gtk_drag_dest_set(mdata->listbox,
		      GTK_DEST_DEFAULT_ALL,
		      gretl_drag_targets, 2,
		      GDK_ACTION_COPY);

    g_signal_connect(G_OBJECT(mdata->listbox), "drag-data-received",
		     G_CALLBACK(mdata_handle_drag),
		     NULL);

    gtk_box_pack_start(GTK_BOX(main_vbox), box, TRUE, TRUE, 0);
    mdata->status = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(main_vbox), mdata->status, FALSE, TRUE, 0);

    /* put stuff into list box, activate menus */
    if (have_data()) {
	populate_varlist();
    }

    /* get a monospaced font for various windows */
    set_fixed_font();

    /* and a proportional font for menus, etc */
    set_app_font(NULL);

    gtk_widget_show_all(mdata->main); 

    /* create gretl toolbar */
    show_toolbar();

    if (winsize && main_x >= 0 && main_y >= 0) {
	gtk_window_move(GTK_WINDOW(mdata->main), main_x, main_y);
    }

    return main_vbox;
}

static void iconview_callback (void)
{
    view_session(NULL);
}

GtkActionEntry main_entries[] = {
    /* File */
    { "File",         NULL, N_("_File"), NULL, NULL, NULL }, 
    { "OpenData",       NULL, N_("_Open data"), NULL, NULL, NULL }, 
    { "OpenGdt",        GTK_STOCK_OPEN, N_("_User file..."), NULL, NULL, G_CALLBACK(open_data) },
    { "DisplayDataFiles", GTK_STOCK_OPEN, N_("_Sample file..."), "", NULL, G_CALLBACK(show_files) },
    { "ImportData",   NULL, N_("_Import"), NULL, NULL, NULL }, 
    { "OpenCSV",      NULL, N_("_CSV..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenASCII",    NULL, N_("_ASCII..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenOctave",   NULL, N_("_Octave..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenGnumeric", NULL, N_("_Gnumeric..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenXLS",      NULL, N_("_Excel..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenODS",      NULL, N_("_Open Document..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenWF1",      NULL, N_("_Eviews..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenDTA",      NULL, N_("_Stata..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenSAV",      NULL, N_("_SPSS..."), NULL, NULL, G_CALLBACK(open_data) },
    { "OpenJMulTi",   NULL, N_("_JMulTi..."), NULL, NULL, G_CALLBACK(open_data) },

    { "AppendData",     NULL, N_("_Append data"), NULL, NULL, NULL }, 
    { "AppendGdt",      NULL, N_("_Standard format..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendCSV",      NULL, N_("_CSV..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendASCII",    NULL, N_("_ASCII..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendOctave",   NULL, N_("_Octave..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendGnumeric", NULL, N_("_Gnumeric..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendXLS",      NULL, N_("_Excel..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendODS",      NULL, N_("_Open Document..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendWF1",      NULL, N_("_Eviews..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendDTA",      NULL, N_("_Stata..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendSAV",      NULL, N_("_SPSS..."), NULL, NULL, G_CALLBACK(open_data) },
    { "AppendJMulTi",   NULL, N_("_JMulTi..."), NULL, NULL, G_CALLBACK(open_data) },

    { "SaveData",  GTK_STOCK_SAVE, N_("_Save data"), "<control>S", NULL, G_CALLBACK(auto_store) },
    { "SaveDataAs", GTK_STOCK_SAVE_AS, N_("Save data _as"), NULL, NULL, NULL }, 
    { "SaveAsGdt", NULL, N_("_Standard format..."), NULL, NULL, G_CALLBACK(fsave_callback) },
    { "SaveAsDb",  NULL, N_("_Database..."), NULL, NULL, G_CALLBACK(fsave_callback) },

    { "ExportData", NULL, N_("_Export data"), NULL, NULL, NULL },
    { "ExportCSV",    NULL, N_("_CSV..."), NULL, NULL, G_CALLBACK(fsave_callback) },
    { "ExportR",      NULL, N_("GNU _R..."), NULL, NULL, G_CALLBACK(fsave_callback) },
    { "ExportOctave", NULL, N_("_Octave..."), NULL, NULL, G_CALLBACK(fsave_callback) },
    { "ExportJMulTi", NULL, N_("_JMulTi..."), NULL, NULL, G_CALLBACK(fsave_callback) },
    { "ExportPcGive", NULL, N_("_PcGive..."), NULL, NULL, G_CALLBACK(fsave_callback) },

    { "MailData", GRETL_STOCK_MAIL, N_("Send To..."), NULL, NULL, G_CALLBACK(email_data) },
    { "NewData", GTK_STOCK_NEW, N_("_New data set"), NULL, NULL, G_CALLBACK(newdata_callback) },
    { "ClearData", GTK_STOCK_CLEAR, N_("C_lear data set"), NULL, NULL, G_CALLBACK(verify_clear_data) },
    { "WorkingDir", NULL, N_("_Working directory..."), NULL, NULL, G_CALLBACK(working_dir_dialog) },

    { "ScriptFiles", NULL, N_("_Script files"), NULL, NULL, NULL },
    { "OpenScript", GTK_STOCK_OPEN, N_("_User file..."), "", NULL, G_CALLBACK(open_script) },
    { "DisplayScripts", GTK_STOCK_OPEN, N_("_Practice file..."), "", NULL, G_CALLBACK(show_files) },
    { "NewScript", GTK_STOCK_NEW, N_("_New script"), "", NULL, NULL },
    { "GretlScript", NULL, N_("gretl script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "GnuplotScript", NULL, N_("gnuplot script"), NULL, NULL, G_CALLBACK(new_script_callback) },
    { "RScript", NULL, N_("R script"), NULL, NULL, G_CALLBACK(new_script_callback) },
#ifdef USE_OX
    { "OxScript", NULL, N_("Ox program"), NULL, NULL, G_CALLBACK(new_script_callback) },
#endif

    { "SessionFiles", NULL, N_("_Session files"), NULL, NULL, NULL },
    { "OpenSession", GTK_STOCK_OPEN, N_("_Open session..."), "", NULL, G_CALLBACK(open_script) },
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
    { "NewGfn", GTK_STOCK_NEW, N_("_New package"), "", NULL, G_CALLBACK(new_function_pkg_callback) },

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
    { "ShowConsole", NULL, N_("_Gretl console"), NULL, NULL, G_CALLBACK(show_gretl_console) },
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
    { "DisplayValues", NULL, N_("_Display values"), NULL, NULL, G_CALLBACK(display_selected) },
    { "EditValues", NULL, N_("_Edit values"), NULL, NULL, G_CALLBACK(spreadsheet_edit) },
    { "AddObs", NULL, N_("_Add observations..."), NULL, NULL, G_CALLBACK(do_add_obs) },
    { "RemoveObs", NULL, N_("_Remove extra observations"), NULL, NULL, G_CALLBACK(do_remove_obs) },
    { "ReadInfo", NULL, N_("_Read info"), NULL, NULL, G_CALLBACK(open_info) },
    { "EditInfo", NULL, N_("Edit _info"), NULL, NULL, G_CALLBACK(edit_header) },
    { "PrintReport", NULL, N_("_Print description"), NULL, NULL, G_CALLBACK(print_report) },
    { "AddMarkers", NULL, N_("Add _case markers..."), NULL, NULL, G_CALLBACK(open_data) },
    { "RemoveMarkers", NULL, N_("Remove case _markers"), NULL, NULL, G_CALLBACK(do_remove_markers) },
    { "DataStructure", NULL, N_("Dataset _structure..."), NULL, NULL, G_CALLBACK(data_structure_dialog) },
    { "DataCompact", NULL, N_("_Compact data..."), NULL, NULL, G_CALLBACK(do_compact_data_set) },
    { "DataExpand", NULL, N_("_Expand data..."), NULL, NULL, G_CALLBACK(do_expand_data_set) },
    { "DataTranspose", NULL, N_("_Transpose data..."), NULL, NULL, G_CALLBACK(gui_transpose_data) },
    { "DataSort", NULL, N_("_Sort data..."), NULL, NULL, G_CALLBACK(gui_sort_data) },
    { "DataRefresh", NULL, N_("_Refresh window"), NULL, NULL, G_CALLBACK(refresh_data) },

    /* View */
    { "View", NULL, N_("_View"), NULL, NULL, NULL },
    { "IconView", NULL, N_("_Icon view"), NULL, NULL, G_CALLBACK(iconview_callback) },
    { "EditScalars", NULL, N_("_Scalars"), NULL, NULL, G_CALLBACK(edit_scalars) },
    { "GraphVars", NULL, N_("_Graph specified vars"), NULL, NULL, NULL },
    { "TSPlot", NULL, N_("_Time series plot..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "ScatterPlot", NULL, N_("X-Y _scatter..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "ImpulsePlot", NULL, N_("X-Y with _impulses..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "FactorPlot", NULL, N_("X-Y with _factor separation..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "FrischPlot", NULL, N_("X-Y with _control..."), NULL, NULL, G_CALLBACK(selector_callback) },
    { "GR_BOX", NULL, N_("_Boxplots..."), NULL, NULL, G_CALLBACK(gretl_callback) },
    { "GR_NBOX", NULL, N_("_Notched boxplots..."), NULL, NULL, G_CALLBACK(gretl_callback) },
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
    { "SMPLDUM", NULL, N_("_Define, based on dummy..."), NULL, NULL, G_CALLBACK(sample_range_dialog) },
    { "SampleRestrict", NULL, N_("_Restrict, based on criterion..."), NULL, NULL, 
      G_CALLBACK(gretl_callback) },
    { "SMPLRAND", NULL, N_("R_andom sub-sample..."), NULL, NULL, G_CALLBACK(sample_range_dialog) },
    { "SampleWReplace", NULL, N_("_Resample with replacement..."), NULL, NULL, 
      G_CALLBACK(gui_resample_data) },
    { "DropMissing", NULL, N_("Drop all obs with _missing values"), NULL, NULL, 
      G_CALLBACK(drop_all_missing) },
    { "CountMissing", NULL, N_("_Count missing values"), NULL, NULL, G_CALLBACK(count_missing) },
    { "GSETMISS", NULL, N_("Set missing _value code..."), NULL, NULL, G_CALLBACK(gretl_callback) },

    /* Variable */
    { "Variable", NULL, N_("_Variable"), NULL, NULL, NULL },
    { "VarFind", GTK_STOCK_FIND, N_("_Find..."), NULL, NULL, G_CALLBACK(listbox_find) },
    { "VarDisplay", NULL, N_("_Display values"), NULL, NULL, G_CALLBACK(display_var) },
    { "VarSummary", NULL, N_("_Summary statistics"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "normtest", NULL, N_("_Normality test"), NULL, NULL, G_CALLBACK(menu_op_action) },
    { "FreqDist", NULL, N_("_Frequency distribution"), NULL, NULL, G_CALLBACK(freq_callback) },
    { "FreqPlot", NULL, N_("Frequency _plot"), NULL, NULL, G_CALLBACK(freq_callback) },
    { "Density", NULL, N_("Estimated _density plot..."), NULL, NULL, G_CALLBACK(do_kernel) },
    { "Gini", NULL, N_("_Gini coefficient"), NULL, NULL, G_CALLBACK(do_gini) },
    { "rmplot", NULL, N_("_Range-mean graph"), NULL, NULL, G_CALLBACK(do_range_mean) },
    { "VarTSPlot", NULL, N_("_Time series plot"), NULL, NULL, G_CALLBACK(ts_plot_callback) },
    { "corrgm", NULL, N_("_Correlogram"), NULL, NULL, G_CALLBACK(do_corrgm) },
    { "Spectrum", NULL, N_("_Spectrum"), NULL, NULL, NULL },
    { "Pergm", NULL, N_("Sample _periodogram"), NULL, NULL, G_CALLBACK(do_pergm) },
    { "Bartlett", NULL, N_("_Bartlett lag window"), NULL, NULL, G_CALLBACK(do_pergm) },
    { "ADF", NULL, N_("_Augmented Dickey-Fuller test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "DFGLS", NULL, N_("ADF-GLS test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "KPSS", NULL, N_("_KPSS test"), NULL, NULL, G_CALLBACK(ur_callback) },
    { "Filter", NULL, N_("_Filter"), NULL, NULL, NULL },
    { "FilterSMA", NULL, N_("_Simple moving average"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterEMA", NULL, N_("_Exponential moving average"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterHP", NULL, N_("_Hodrick-Prescott"), NULL, NULL, G_CALLBACK(filter_callback) },
    { "FilterBK", NULL, N_("_Baxter-King"), NULL, NULL, G_CALLBACK(filter_callback) },
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
#ifdef ENABLE_GMP
    { "mpols", NULL, N_("High _precision OLS..."), NULL, NULL, G_CALLBACK(model_callback) }, 
#endif
    { "anova", NULL, N_("ANOVA..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "TSModels", NULL, N_("_Time series"), NULL, NULL, NULL },
    { "CORC", NULL, N_("_Cochrane-Orcutt..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "HILU", NULL, N_("_Hildreth-Lu..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "PWE", NULL, N_("_Prais-Winsten..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "ar", NULL, N_("_Autoregressive estimation..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "arima", NULL, N_("ARI_MA..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "arch", NULL, N_("_ARCH..."), NULL, NULL, G_CALLBACK(model_callback) },
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
    { "arbond", NULL, N_("_Arellano-Bond..."), NULL, NULL, G_CALLBACK(model_callback) },
    { "NonlinearModels", NULL, N_("_Nonlinear models"), NULL, NULL, NULL },
    { "logit", NULL, N_("_Logit"), NULL, NULL, NULL }, 
    { "blogit", NULL, N_("_Binary..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "ologit", NULL, N_("_Ordered..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "mlogit", NULL, N_("_Multinomial..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "probit", NULL, N_("_Probit"), NULL, NULL, NULL }, 
    { "bprobit", NULL, N_("_Binary..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "oprobit", NULL, N_("_Ordered..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "tobit", NULL, N_("To_bit..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "heckit", NULL, N_("_Heckit..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "poisson", NULL, N_("Poi_sson..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "logistic", NULL, N_("Lo_gistic..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "nls", NULL, N_("_Nonlinear Least Squares..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "intreg", NULL, N_("_Interval regression..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "RobustModels", NULL, N_("_Robust estimation"), NULL, NULL, NULL },
    { "lad", NULL, N_("Least _Absolute Deviation..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "quantreg", NULL, N_("_Quantile regression..."), NULL, NULL, G_CALLBACK(model_callback) }, 
    { "spearman", NULL, N_("_Rank correlation..."), NULL, NULL, G_CALLBACK(selector_callback) }, 
    { "mle", NULL, N_("_Maximum likelihood..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "gmm", NULL, N_("_GMM..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 
    { "system", NULL, N_("_Simultaneous equations..."), NULL, NULL, G_CALLBACK(gretl_callback) }, 

    /* Help */
    { "Help", NULL, N_("_Help"), NULL, NULL, NULL },
    { "CommandRef", NULL, N_("_Command reference"), NULL, NULL, NULL },
    { "TextCmdRef", GRETL_STOCK_BOOK, N_("Plain _text"), NULL, NULL, G_CALLBACK(plain_text_cmdref) },
    { "PDFCmdRef", GRETL_STOCK_PDF, N_("_PDF"), NULL, NULL, G_CALLBACK(display_pdf_help) },
    { "FuncRef",   NULL, N_("Function reference"), NULL, NULL, G_CALLBACK(genr_funcs_ref) },
    { "UserGuide", GRETL_STOCK_PDF, N_("_User's guide"), NULL, NULL, G_CALLBACK(display_pdf_help) },
    { "UpdateCheck", GTK_STOCK_NETWORK, N_("Check for _updates"), NULL, NULL, G_CALLBACK(update_query) },
    { "About", GTK_STOCK_ABOUT, N_("_About gretl"), NULL, NULL, G_CALLBACK(about_dialog) }
};

static void add_conditional_items (GtkUIManager *ui)
{
    gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			  "/MenuBar/Tools/Preferences",
			  N_("_Menu font..."),
			  "MenuFont",
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);

#ifdef HAVE_X12A
    gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			  "/MenuBar/Variable/X12A",
			  N_("_X-12-ARIMA analysis"),
			  "X12A",
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);
#endif

#ifdef HAVE_TRAMO
    gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			  "/MenuBar/Variable/Tramo",
			  N_("_TRAMO analysis"),
			  "Tramo",
			  GTK_UI_MANAGER_MENUITEM, 
			  FALSE);
#endif

#ifdef USE_OX
    if (ox_support) {
	gtk_ui_manager_add_ui(ui, gtk_ui_manager_new_merge_id(ui),
			      "/MenuBar/File/ScriptFiles/NewScript/OxScript",
			      N_("Ox program"),
			      "OxScript",
			      GTK_UI_MANAGER_MENUITEM, 
			      FALSE);
    }
#endif
}

static gchar *get_main_ui (void)
{
    char fname[FILENAME_MAX];
    gchar *main_ui = NULL;

    sprintf(fname, "%sui%cgretlmain.xml", paths.gretldir, SLASH);
    gretl_file_get_contents(fname, &main_ui);

    return main_ui;
}

static int set_up_main_menu (void)
{
    GtkActionGroup *actions;
    gchar *main_ui = NULL;
    GError *error = NULL;

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
    } else {
	add_conditional_items(mdata->ui);
    }

    g_free(main_ui);
    mdata->mbar = gtk_ui_manager_get_widget(mdata->ui, "/MenuBar");

    return 0;
}

int gui_restore_sample (double ***pZ, DATAINFO *pdinfo)
{
    int err = 0;

    if (pZ != NULL && *pZ != NULL) {
	err = restore_full_sample(pZ, pdinfo, NULL);
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
    int err = gui_restore_sample(&Z, datainfo); 

    if (!err) {
	set_sample_label(datainfo);    
	gretl_command_strcpy("smpl --full");
	check_and_record_command();
    }
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

void set_wm_icon (GtkWidget *w, gpointer data)
{
    GdkPixmap *icon;

    icon = gdk_pixmap_create_from_xpm_d(w->window, NULL, NULL, gretl_xpm);
    gdk_window_set_icon(w->window, NULL, icon, NULL);
}

#endif

/* Drag 'n' drop */

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
    gchar *dfname;
    char tmp[MAXLEN];
    int pos, skip = 5;
    int ftype = 0;

    /* handle drag of pointer from database window */
    if (info == GRETL_DBSERIES_PTR && data != NULL && 
	data->type == GDK_SELECTION_TYPE_INTEGER) {
	import_db_series(*(void **) data->data);
	return;
    }

    if (info != GRETL_FILENAME) {
	return;
    }

    /* ignore the wrong sort of data */
    if (data == NULL || (dfname = (gchar *) data->data) == NULL || 
	strlen(dfname) <= 5 || strncmp(dfname, "file:", 5)) {
	return;
    }

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

    /* handle spaces and such */
    unescape_url(tmp);

#ifdef G_OS_WIN32
    filename_to_win32(tryfile, tmp);
#else
    strcpy(tryfile, tmp);
#endif

    if (has_suffix(tryfile, ".gretl") && gretl_is_pkzip_file(tryfile)) {
	verify_open_session();
    } else if ((ftype = foreign_script_type(tryfile))) {
	do_open_script(ftype);
    } else {
	verify_open_data(NULL, 0);
    }
}

static void auto_store (void)
{
    /* by default, use gzip compression */
    gretlopt oflag = OPT_Z;

    /* but if there's already a datafile, and it's not gzipped, then
       arrange for the new file to be uncompressed too
    */
    if (*paths.datfile && !is_gzipped(paths.datfile)) {
	oflag = OPT_NONE;
    }

    if ((data_status & USER_DATA) && has_suffix(paths.datfile, ".gdt")) {
	do_store(paths.datfile, oflag);
    } else {
	file_selector(_("Save data file"), SAVE_DATA, FSEL_DATA_NONE, NULL);
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

    if (mods & GDK_BUTTON3_MASK) {
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
