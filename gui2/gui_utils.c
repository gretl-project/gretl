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

/* gui_utils.c for gretl */

#include "gretl.h"

#include <sys/stat.h>
#include <unistd.h>

#include "guiprint.h"
#include "series_view.h"
#include "console.h"
#include "session.h"
#include "model_table.h"

#include "../pixmaps/mini.tex.xpm"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
# include <gtksourceview/gtksourcelanguage.h>
# include <gtksourceview/gtksourcelanguagesmanager.h>
#endif

char *storelist = NULL;
GtkWidget *active_edit_id = NULL;
GtkWidget *active_edit_name = NULL;
GtkWidget *active_edit_text = NULL;

extern int session_saved;
extern GtkWidget *mysheet;

#define MARK_SCRIPT_CHANGED(w) (w->active_var = 1)
#define MARK_SCRIPT_SAVED(w) (w->active_var = 0)
#define SCRIPT_IS_CHANGED(w) (w->active_var == 1)

#define MULTI_COPY_ENABLED(c) (c == SUMMARY || c == VAR_SUMMARY \
	                      || c == CORR || c == FCASTERR \
	                      || c == FCAST || c == COEFFINT \
	                      || c == COVAR || c == VIEW_MODEL \
                              || c == VIEW_MODELTABLE)

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static void file_viewer_save (GtkWidget *widget, windata_t *vwin);
static gint query_save_script (GtkWidget *w, GdkEvent *event, windata_t *vwin);
static void add_vars_to_plot_menu (windata_t *vwin);
static void add_dummies_to_plot_menu (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);
static void buf_edit_save (GtkWidget *widget, gpointer data);
static void model_copy_callback (gpointer p, guint u, GtkWidget *w);

#ifndef USE_GTKSOURCEVIEW
static void correct_line_color (windata_t *vwin);
#endif

extern void do_coeff_intervals (gpointer data, guint i, GtkWidget *w);
extern void save_plot (char *fname, GPT_SPEC *plot);
extern void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w);
extern void do_leverage (gpointer data, guint u, GtkWidget *w);


GtkItemFactoryEntry model_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, 
      "<StockItem>", GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as icon"), NULL, remember_model, 0, NULL, GNULL },
    { N_("/File/Save as icon and close"), NULL, remember_model, 1, NULL, GNULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL, GNULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, "<StockItem>", GTK_STOCK_COPY },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/Tests/omit variables"), NULL, selector_callback, OMIT, NULL, GNULL },
    { N_("/Tests/add variables"), NULL, selector_callback, ADD, NULL, GNULL },
    { N_("/Tests/sum of coefficients"), NULL, selector_callback, COEFFSUM, NULL, GNULL },
    { N_("/Tests/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/non-linearity (squares)"), NULL, do_lmtest, AUX_SQ, NULL, GNULL },
    { N_("/Tests/non-linearity (logs)"), NULL, do_lmtest, AUX_LOG, NULL, GNULL },
    { N_("/Tests/Ramsey's RESET"), NULL, do_reset, 0, NULL, GNULL },
    { N_("/Tests/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/autocorrelation"), NULL, model_test_callback, LMTEST, NULL, GNULL },
    { N_("/Tests/heteroskedasticity"), NULL, do_lmtest, AUX_WHITE, NULL, GNULL },
    { N_("/Tests/influential observations"), NULL, do_leverage, 0, NULL, GNULL },
    { N_("/Tests/Chow test"), NULL, model_test_callback, CHOW, NULL, GNULL },
    { N_("/Tests/CUSUM test"), NULL, do_cusum, 0, NULL, GNULL },
    { N_("/Tests/ARCH"), NULL, model_test_callback, ARCH, NULL, GNULL },
    { N_("/Tests/normality of residual"), NULL, do_resid_freq, 0, NULL, GNULL },
    { N_("/Tests/panel diagnostics"), NULL, do_panel_diagnostics, 0, NULL, GNULL },
    { N_("/_Graphs"), NULL, NULL, 0, "<Branch>", GNULL }, 
    { N_("/Graphs/residual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Graphs/fitted, actual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/_Model data"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Model data/Display actual, fitted, residual"), NULL, 
      display_fit_resid, 0, NULL, GNULL },
    { N_("/Model data/Forecasts with standard errors"), NULL, 
      model_test_callback, FCASTERR, NULL, GNULL },
    { N_("/Model data/Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL, GNULL },
    { N_("/Model data/coefficient covariance matrix"), NULL, 
      do_outcovmx, 0, NULL, GNULL },
    { N_("/Model data/Add to data set"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Model data/Add to data set/fitted values"), NULL, 
      fit_resid_callback, 1, NULL, GNULL },
    { N_("/Model data/Add to data set/residuals"), NULL, 
      fit_resid_callback, 0, NULL, GNULL },
    { N_("/Model data/Add to data set/squared residuals"), NULL, 
      fit_resid_callback, 2, NULL, GNULL },
    { N_("/Model data/Add to data set/error sum of squares"), NULL, 
      model_stat_callback, ESS, NULL, GNULL },
    { N_("/Model data/Add to data set/standard error of residuals"), NULL, 
      model_stat_callback, SIGMA, NULL, GNULL },
    { N_("/Model data/Add to data set/R-squared"), NULL, 
      model_stat_callback, R2, NULL, GNULL },
    { N_("/Model data/Add to data set/T*R-squared"), NULL, 
      model_stat_callback, TR2, NULL, GNULL },
    { N_("/Model data/Add to data set/log likelihood"), NULL, 
      model_stat_callback, LNL, NULL, GNULL },
    { N_("/Model data/Add to data set/degrees of freedom"), NULL, 
      model_stat_callback, DF, NULL, GNULL },
    { N_("/Model data/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Model data/Define new variable..."), NULL, model_test_callback, 
      MODEL_GENR, NULL, GNULL },
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/_View"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/View/_Tabular"), NULL, view_latex, 
      LATEX_VIEW_TABULAR, NULL, GNULL },
    { N_("/LaTeX/View/_Equation"), NULL, view_latex, 
      LATEX_VIEW_EQUATION, NULL, GNULL },
    { N_("/LaTeX/_Save"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Save/_Tabular"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Save/Tabular/as _document"), NULL, file_save, 
      SAVE_TEX_TAB, NULL, GNULL },
    { N_("/LaTeX/Save/Tabular/as _fragment"), NULL, file_save, 
      SAVE_TEX_TAB_FRAG, NULL, GNULL },
    { N_("/LaTeX/Save/_Equation"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Save/Equation/as _document"), NULL, file_save, 
      SAVE_TEX_EQ, NULL, GNULL },
    { N_("/LaTeX/Save/Equation/as _fragment"), NULL, file_save, 
      SAVE_TEX_EQ_FRAG, NULL, GNULL },
    { N_("/LaTeX/_Copy"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Copy/_Tabular"), NULL, text_copy, COPY_LATEX, NULL, GNULL },
    { N_("/LaTeX/Copy/_Equation"), NULL, text_copy, COPY_LATEX_EQUATION, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

static void model_copy_callback (gpointer p, guint u, GtkWidget *w)
{
    copy_format_dialog((windata_t *) p);
}

#ifdef ENABLE_NLS
gchar *menu_translate (const gchar *path, gpointer p)
{
    return (_(path));
}
#endif

/* ........................................................... */

int copyfile (const char *src, const char *dest) 
{
    FILE *srcfd, *destfd;
    char buf[GRETL_BUFSIZE];
    size_t n;

    if (!strcmp(src, dest)) return 1;
   
    if ((srcfd = fopen(src, "rb")) == NULL) {
	sprintf(errtext, _("Couldn't open %s"), src);
	errbox(errtext);
	return 1; 
    }
    if ((destfd = fopen(dest, "wb")) == NULL) {
	sprintf(errtext, _("Couldn't write to %s"), dest);
	errbox(errtext);
	fclose(srcfd);
	return 1;
    }
    while ((n = fread(buf, 1, sizeof buf, srcfd)) > 0) {
	fwrite(buf, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}

/* ........................................................... */

int isdir (const char *path)
{
    struct stat buf;

    if (stat(path, &buf) == 0 && S_ISDIR(buf.st_mode)) 
	return 1;
    else 
	return 0;
}

/* ........................................................... */

int getbufline (char *buf, char *line, int init)
{
    static int pos;
    int i = 0;

    if (init) pos = 0;
    else {
	while (buf[i+pos] != '\n') {
	    line[i] = buf[i+pos];
	    if (buf[i+pos] == 0)
		return 0;
	    i++;
	}
	pos += i + 1;
	line[i] = 0;
    }
    return i;
}

/* ........................................................... */

void append_dir (char *fname, const char *dir)
{
    size_t len;

    if (dir == NULL) return;

    len = strlen(fname);
    if (fname[len - 1] == '/' || fname[len - 1] == '\\')
        strcat(fname, dir);
    else {
        strcat(fname, SLASHSTR);
        strcat(fname, dir);
    }
    strcat(fname, SLASHSTR);
}

/* ........................................................... */

/* Below: Keep a record of (most) windows that are open, so they 
   can be destroyed en masse when a new data file is opened, to
   prevent weirdness that could arise if (e.g.) a model window
   that pertains to a previously opened data file remains open
   after the data set has been changed.  Script windows are
   exempt, otherwise they are likely to disappear when their
   "run" control is activated, which we don't want.
*/

enum winstack_codes {
    STACK_INIT,
    STACK_ADD,
    STACK_REMOVE,
    STACK_DESTROY,
    STACK_QUERY
};

static int winstack (int code, GtkWidget *w, gpointer ptest)
{
    static int n_windows;
    static GtkWidget **wstack;
    int found = 0;
    int i;

    switch (code) {

    case STACK_DESTROY:	
	for (i=0; i<n_windows; i++) 
	    if (wstack[i] != NULL) 
		gtk_widget_destroy(wstack[i]);
	free(wstack);
	/* fall-through intended */

    case STACK_INIT:
	wstack = NULL;
	n_windows = 0;
	break;

    case STACK_ADD:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] == NULL) {
		wstack[i] = w;
		break;
	    }
	}
	if (i == n_windows) {
	    n_windows++;
	    wstack = myrealloc(wstack, n_windows * sizeof *wstack);
	    if (wstack != NULL) 
		wstack[n_windows-1] = w;
	}
	break;

    case STACK_REMOVE:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] == w) {
		wstack[i] = NULL;
		break;
	    }
	}
	break;

    case STACK_QUERY:
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] != NULL) {
		gpointer p = g_object_get_data(G_OBJECT(wstack[i]), "object");

		if (p == ptest) {
		    found = 1;
		    break;
		}
	    }
	}
	break;

    default:
	break;
    }

    return found;
}

void winstack_init (void)
{
    winstack(STACK_INIT, NULL, NULL);
}
    
void winstack_destroy (void)
{
    winstack(STACK_DESTROY, NULL, NULL);
}

int winstack_match_data (gpointer p)
{
    return winstack(STACK_QUERY, NULL, p);
}

static void winstack_add (GtkWidget *w)
{
    winstack(STACK_ADD, w, NULL);
}

static void winstack_remove (GtkWidget *w)
{
    winstack(STACK_REMOVE, w, NULL);
}

/* ........................................................... */

static void delete_file (GtkWidget *widget, char *fname) 
{
    remove(fname);
    g_free(fname);
}

/* ........................................................... */

static void delete_file_viewer (GtkWidget *widget, gpointer data) 
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->role == EDIT_SCRIPT && SCRIPT_IS_CHANGED(vwin)) {
	gint resp;

	resp = query_save_script(NULL, NULL, vwin);
	if (!resp) gtk_widget_destroy(vwin->dialog);
    } else {
	gtk_widget_destroy(vwin->dialog);
    }
}

/* ........................................................... */

static void delete_unnamed_model (GtkWidget *widget, gpointer data) 
{
    MODEL *pmod = (MODEL *) data;

    if (pmod->dataset != NULL) {
	free_model_dataset(pmod);
    }

    if (pmod->name == NULL) {
	free_model(pmod);
	pmod = NULL;
    }
}

/* ........................................................... */

void delete_widget (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
}

/* ........................................................... */

static gint catch_button_3 (GtkWidget *w, GdkEventButton *event)
{
    GdkModifierType mods;

    gdk_window_get_pointer (w->window, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK) {
	return TRUE;
    }
    return FALSE;
}

/* ........................................................... */

#ifdef G_OS_WIN32
static void win_ctrl_c (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    if (gtk_text_buffer_get_selection_bounds(buf, NULL, NULL)) {
	text_copy(vwin, COPY_SELECTION, NULL);
    } else if (MULTI_COPY_ENABLED(vwin->role)) {
	text_copy(vwin, COPY_RTF, NULL);
    } else {
	text_copy(vwin, COPY_TEXT, NULL);
    }
}
#endif

/* ........................................................... */

static gint catch_view_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    }
    else if (key->keyval == GDK_s && Z != NULL && vwin->role == VIEW_MODEL) {
	remember_model(vwin, 1, NULL);
    }
#ifdef G_OS_WIN32
    else if (key->keyval == GDK_c) {
	GdkModifierType mods;

	gdk_window_get_pointer(w->window, NULL, NULL, &mods); 
	if (mods & GDK_CONTROL_MASK) {
	    win_ctrl_c(vwin);
	    return TRUE;
	}	
    }
#endif
    return FALSE;
}

/* ........................................................... */

static gint catch_edit_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);

    if (key->keyval == GDK_F1 && vwin->role == EDIT_SCRIPT) { 
	vwin->help_active = 1;
	edit_script_help(NULL, NULL, vwin);
    }

#ifndef USE_GTKSOURCEVIEW
    else if (key->keyval == GDK_Return) {
	/* newline: correct line color */
	correct_line_color(vwin);
    }
#endif

    else if (mods & GDK_CONTROL_MASK) {
	if (gdk_keyval_to_upper(key->keyval) == GDK_S) { 
	    if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		buf_edit_save(NULL, vwin);
	    } else {
		file_viewer_save(NULL, vwin);
	    }
	} 
	else if (gdk_keyval_to_upper(key->keyval) == GDK_Q) {
	    if (vwin->role == EDIT_SCRIPT && SCRIPT_IS_CHANGED(vwin)) {
		gint resp;

		resp = query_save_script(NULL, NULL, vwin);
		if (!resp) gtk_widget_destroy(vwin->dialog);
	    } else { 
		gtk_widget_destroy(w);
	    }
	}
#ifdef G_OS_WIN32 
	else if (key->keyval == GDK_c) {
	    win_ctrl_c(vwin);
	    return TRUE;
	}
#endif
    }

    return FALSE;
}

/* ........................................................... */

void *mymalloc (size_t size) 
{
    void *mem;
   
    if((mem = malloc(size)) == NULL) 
	errbox(_("Out of memory!"));
    return mem;
}

/* ........................................................... */

void *myrealloc (void *ptr, size_t size) 
{
    void *mem;
   
    if ((mem = realloc(ptr, size)) == NULL) 
	errbox(_("Out of memory!"));
    return mem;
}

/* ........................................................... */

void mark_dataset_as_modified (void)
{
    data_status |= MODIFIED_DATA;
    set_sample_label(datainfo);
}

/* ........................................................... */

void register_data (const char *fname, const char *user_fname,
		    int record)
{    
    char datacmd[MAXLEN];

    /* basic accounting */
    data_status |= HAVE_DATA;
    orig_vars = datainfo->v;

    /* set appropriate data_status bits */
    if (fname == NULL) {
	data_status |= GUI_DATA;
	mark_dataset_as_modified();
    } else if (!(data_status & IMPORT_DATA)) {
	if (strstr(paths.datfile, paths.datadir) != NULL) {
	    data_status |= BOOK_DATA;
	} else {
	    data_status |= USER_DATA; 
	}
    }

    /* sync main window with datafile */
    populate_varlist();
    set_sample_label(datainfo);
    main_menubar_state(TRUE);
    session_menu_state(TRUE);

    /* record opening of data file in command log */
    if (record && fname != NULL) {
	mkfilelist(FILE_LIST_DATA, fname);
	sprintf(datacmd, "open %s", user_fname ? user_fname : fname);
	check_cmd(datacmd);
	cmd_init(datacmd); 
    } 
}

#define APPENDING(action) (action == APPEND_CSV || \
                           action == APPEND_GNUMERIC || \
                           action == APPEND_EXCEL || \
                           action == APPEND_ASCII)

/* ........................................................... */

int get_worksheet_data (const char *fname, int datatype, int append)
{
    int err;
    void *handle;
    PRN *errprn;
    int (*sheet_get_data)(const char*, double ***, DATAINFO *, PRN *);

    if (datatype == GRETL_GNUMERIC) {
	if (gui_open_plugin("gnumeric_import", &handle)) return 1;
	sheet_get_data = get_plugin_function("wbook_get_data", handle);
    }
    else if (datatype == GRETL_EXCEL) {
	if (gui_open_plugin("excel_import", &handle)) return 1;
	sheet_get_data = get_plugin_function("excel_get_data", handle);
    }
    else {
	errbox(_("Unrecognized data type"));
	return 1;
    }

    if (sheet_get_data == NULL) {
        errbox(_("Couldn't load plugin function"));
        close_plugin(handle);
        return 1;
    }

    if (bufopen(&errprn)) {
	close_plugin(handle);
	return 1;
    }

    err = (*sheet_get_data)(fname, &Z, datainfo, errprn);
    close_plugin(handle);

    if (err == -1) /* the user canceled the import */
	return 0;

    if (err) {
	if (*errprn->buf != '\0') {
	    errbox(errprn->buf);
	}
	else {
	    errbox(_("Failed to import spreadsheet data"));
	}
	return 1;
    } else {
	if (*errprn->buf != '\0') infobox(errprn->buf);
    }

    gretl_print_destroy(errprn);

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	data_status |= IMPORT_DATA;
	strcpy(paths.datfile, fname);
	register_data(fname, NULL, 1);
    }

    return 0;
}

/* ........................................................... */

void do_open_data (GtkWidget *w, gpointer data, int code)
     /* cases: 
	- called from dialog: user has said Yes to opening data file,
	although a data file is already open (or user wants to append
	data)
	- reached without dialog, in expert mode or when no datafile
	is open yet
     */
{
    gint datatype, err = 0;
    dialog_t *d = NULL;
    windata_t *fwin = NULL;
    int append = APPENDING(code);

    if (data != NULL) {    
	if (w == NULL) { /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	} else {
	    d = (dialog_t *) data;
	    fwin = (windata_t *) d->data;
	}
    }

    if (code == OPEN_CSV || code == APPEND_CSV || code == OPEN_ASCII ||
	code == APPEND_ASCII) {
	datatype = GRETL_CSV_DATA;
    }
    else if (code == OPEN_GNUMERIC || code == APPEND_GNUMERIC) {
	datatype = GRETL_GNUMERIC;
    }
    else if (code == OPEN_EXCEL || code == APPEND_EXCEL) {
	datatype = GRETL_EXCEL;
    }
    else if (code == OPEN_BOX) {
	datatype = GRETL_BOX_DATA;
    }
    else if (code == OPEN_DES) {
        datatype = GRETL_DES_DATA;
    } else {
	/* no filetype specified: have to guess */
	PRN *prn;	

	if (bufopen(&prn)) return;
	datatype = detect_filetype(trydatfile, &paths, prn);
	gretl_print_destroy(prn);
    }

    /* destroy the current data set, etc., unless we're explicitly appending */
    if (!append) close_session();

    if (datatype == GRETL_GNUMERIC || datatype == GRETL_EXCEL) {
	get_worksheet_data(trydatfile, datatype, append);
	return;
    }
    else if (datatype == GRETL_CSV_DATA) {
	do_open_csv_box(trydatfile, OPEN_CSV, append);
	return;
    }
    else if (datatype == GRETL_BOX_DATA) {
	do_open_csv_box(trydatfile, OPEN_BOX, 0);
	return;
    }
    else if (datatype == GRETL_DES_DATA) {
        get_worksheet_data(trydatfile, datatype, 0);
        return;
    }    
    else { /* native data */
	PRN prn;

	gretl_print_attach_file(&prn, stderr);
	if (datatype == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, datainfo, trydatfile, &paths, 
			      data_status, &prn, 1);
	} else {
	    err = get_data(&Z, datainfo, trydatfile, &paths, data_status, &prn);
	}
    }

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_DATA, trydatfile);
	return;
    }	

    /* trash the practice files window that launched the query? */
    if (fwin != NULL) gtk_widget_destroy(fwin->w);

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	strcpy(paths.datfile, trydatfile);
	register_data(paths.datfile, NULL, 1);
    } 
}

/* ........................................................... */

void verify_open_data (gpointer userdata, int code)
     /* give user choice of not opening selected datafile,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert) {
	int resp = 
	    yes_no_dialog (_("gretl: open data"), 
			   _("Opening a new data file will automatically\n"
			     "close the current one.  Any unsaved work\n"
			     "will be lost.  Proceed to open data file?"), 0);

	if (resp != GRETL_YES) return;
    } 

    do_open_data(NULL, userdata, code);
}

/* ........................................................... */

void verify_open_session (gpointer userdata)
     /* give user choice of not opening session file,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert) {
	int resp = 
	    yes_no_dialog (_("gretl: open session"), 
			   _("Opening a new session file will automatically\n"
			     "close the current session.  Any unsaved work\n"
			     "will be lost.  Proceed to open session file?"), 0);

	if (resp != GRETL_YES) return;
    }

    do_open_session(NULL, userdata);
}

/* ........................................................... */

void save_session (char *fname) 
{
    int spos;
    char msg[MAXLEN], savedir[MAXLEN], fname2[MAXLEN];
    char session_base[MAXLEN];
    FILE *fp;
    PRN *prn;

    spos = slashpos(fname);
    if (spos) 
	safecpy(savedir, fname, spos);
    else *savedir = 0;

#ifdef CMD_DEBUG
    dump_cmd_stack("stderr", 0);
#endif

    /* save commands, by dumping the command stack */
    if (haschar('.', fname) < 0) strcat(fname, ".gretl");
    if (dump_cmd_stack(fname, 1)) return;

    get_base(session_base, fname, '.');

    /* get ready to save "session" */
    fp = fopen(fname, "a");
    if (fp == NULL) {
	sprintf(errtext, _("Couldn't open session file %s"), fname);
	errbox(errtext);
	return;
    }

    print_saved_object_specs(session_base, fp);

    fclose(fp);

    /* delete any extraneous graph files */
    session_file_manager(REALLY_DELETE_ALL, NULL);

    switch_ext(fname2, fname, "Notes");
    if (print_session_notes(fname2)) {
	errbox(_("Couldn't write session notes file"));
    }

    /* save output */
    switch_ext(fname2, fname, "txt");
    prn = gretl_print_new(GRETL_PRINT_FILE, fname2);
    if (prn == NULL) {
	errbox(_("Couldn't open output file for writing"));
	return;
    }

    gui_logo(prn->fp);
    session_time(prn->fp);
    pprintf(prn, _("Output from %s\n"), fname);
    execute_script(fname, NULL, prn, SAVE_SESSION_EXEC); 
    gretl_print_destroy(prn);

    sprintf(msg, _("session saved to %s -\n"), savedir);
    strcat(msg, _("commands: "));
    strcat(msg, (spos)? fname + spos + 1 : fname);
    strcat(msg, _("\noutput: "));
    spos = slashpos(fname2);
    strcat(msg, (spos)? fname2 + spos + 1 : fname2);
    infobox(msg);

    mkfilelist(FILE_LIST_SESSION, fname);
    session_saved = 1;
    session_changed(0);

    return;
}

/* ........................................................... */

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor (vwin->w, GDK_QUESTION_ARROW);
    vwin->help_active = 1;
}

/* ........................................................... */

static void buf_edit_save (GtkWidget *widget, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *text;
    char **pbuf = (char **) vwin->data;

    text = textview_get_text(GTK_TEXT_VIEW(vwin->w));

    if (text == NULL || !strlen(text)) {
	errbox(_("Buffer is empty"));
	g_free(text);
	return;
    }

    /* swap the edited text into the buffer */
    free(*pbuf); 
    *pbuf = text;

    if (vwin->role == EDIT_HEADER) {
	infobox(_("Data info saved"));
	mark_dataset_as_modified();
    } 
    else if (vwin->role == EDIT_NOTES) {
	infobox(_("Notes saved"));
	session_changed(1);
    }
}

/* ........................................................... */

static void file_viewer_save (GtkWidget *widget, windata_t *vwin)
{
    /* special case: a newly created script */
    if (strstr(vwin->fname, "script_tmp") || !strlen(vwin->fname)) {
	file_save(vwin, SAVE_SCRIPT, NULL);
	strcpy(vwin->fname, scriptfile);
    } else {
	char buf[MAXLEN];
	FILE *fp;
	gchar *text;

	if ((fp = fopen(vwin->fname, "w")) == NULL) {
	    errbox(_("Can't open file for writing"));
	    return;
	} else {
	    text = textview_get_text(GTK_TEXT_VIEW(vwin->w));
	    fprintf(fp, "%s", text);
	    fclose(fp);
	    g_free(text);
	    sprintf(buf, _("Saved %s\n"), vwin->fname);
	    infobox(buf);
	    if (vwin->role == EDIT_SCRIPT) { 
		MARK_SCRIPT_SAVED(vwin);
	    }
	}
    }
}

/* .................................................................. */

void windata_init (windata_t *vwin)
{
    vwin->dialog = NULL;
    vwin->listbox = NULL;
    vwin->vbox = NULL;
    vwin->mbar = NULL;
    vwin->w = NULL;
    vwin->status = NULL;
    vwin->popup = NULL;
    vwin->ifac = NULL;
    vwin->data = NULL;
    vwin->fname[0] = '\0';
    vwin->role = 0;
    vwin->active_var = 0;
    vwin->help_active = 0;
#ifdef USE_GTKSOURCEVIEW
    vwin->sbuf = NULL;
#endif
}

/* .................................................................. */

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin) {
	if (vwin->w) { /* gtktextview? */
	    gchar *undo = g_object_get_data(G_OBJECT(vwin->w), "undo");
	    
	    if (undo) g_free(undo);
	}

	/* menu stuff */
	if (vwin->popup) 
	    gtk_widget_destroy(GTK_WIDGET(vwin->popup));
	if (vwin->ifac) 
	    g_object_unref(G_OBJECT(vwin->ifac));

	/* data specific to certain windows */
	if (vwin->role == SUMMARY || vwin->role == VAR_SUMMARY)
	    free_summary(vwin->data); 
	else if (vwin->role == CORR)
	    free_corrmat(vwin->data);
	else if (vwin->role == FCASTERR || vwin->role == FCAST)
	    free_fit_resid(vwin->data);
	else if (vwin->role == COEFFINT)
	    free_confint(vwin->data);
	else if (vwin->role == COVAR)
	    free_vcv(vwin->data);
	else if (vwin->role == MPOLS)
	    free_gretl_mp_results(vwin->data);
	else if (vwin->role == VIEW_SERIES)
	    free_series_view(vwin->data);

	if (vwin->dialog)
	    winstack_remove(vwin->dialog);
	free(vwin);
    }
}

static void modeltable_tex_view (void)
{
    tex_print_model_table (NULL, 1, NULL);
}

static int tex_icon_init (void)
{
    static GtkIconFactory *ifac;

    if (ifac == NULL) {
	GtkIconSet *iset;
	GdkPixbuf *pbuf;

	pbuf = gdk_pixbuf_new_from_xpm_data((const char **) mini_tex_xpm);
	iset = gtk_icon_set_new_from_pixbuf(pbuf);
	ifac = gtk_icon_factory_new();
	gtk_icon_factory_add(ifac, "STOCK_TEX", iset);
	gtk_icon_factory_add_default(ifac);
    }

    return 0;
}

#if defined(G_OS_WIN32) || defined(USE_GNOME) 
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(vwin, 0, w);
}
#endif

static void choose_copy_format_callback (GtkWidget *w, windata_t *vwin)
{
    copy_format_dialog(vwin);
}

/* ........................................................... */

static void make_viewbar (windata_t *vwin, int text_out)
{
    GtkWidget *button, *viewbar, *hbox;
    int i;
    static char *viewstrings[] = {
	N_("Save"),
	N_("Save as..."),
	N_("Send to gnuplot"),
	N_("Print..."),
	N_("Run"),
	N_("Copy"), 
	N_("Paste"),
	N_("Find..."),
	N_("Replace..."),
	N_("Undo"),
	N_("Help on command"),
	N_("LaTeX"),
	N_("Close"),
	NULL
    };
    const gchar *stockicon = NULL;
    void (*toolfunc)() = NULL;
    gchar *toolstr;

    int run_ok = (vwin->role == EDIT_SCRIPT ||
		  vwin->role == VIEW_SCRIPT ||
		  vwin->role == VIEW_LOG);

    int edit_ok = (vwin->role == EDIT_SCRIPT ||
		   vwin->role == EDIT_HEADER ||
		   vwin->role == EDIT_NOTES ||
		   vwin->role == GR_PLOT || 
		   vwin->role == GR_BOX ||
		   vwin->role == SCRIPT_OUT);

    int save_as_ok = (vwin->role != EDIT_HEADER && 
		      vwin->role != EDIT_NOTES);

#if defined(G_OS_WIN32) || defined(USE_GNOME)
    int print_ok = 1;
#else
    int print_ok = 0;
#endif

    if (vwin->role == VIEW_MODELTABLE) tex_icon_init();

    if (text_out || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->dialog), "text_out", GINT_TO_POINTER(1));
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    viewbar = gtk_toolbar_new();
    gtk_box_pack_start(GTK_BOX(hbox), viewbar, FALSE, FALSE, 0);

    for (i=0; viewstrings[i] != NULL; i++) {
	switch (i) {
	case 0:
	    if (edit_ok && vwin->role != SCRIPT_OUT) {
		stockicon = GTK_STOCK_SAVE;
		if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		    toolfunc = buf_edit_save;
		} else if (vwin->role == GR_PLOT) {
		    toolfunc = save_plot_commands_callback;
		} else {
		    toolfunc = file_viewer_save;
		}
	    } else
		toolfunc = NULL;
	    break;
	case 1:
	    if (save_as_ok) {
		stockicon = GTK_STOCK_SAVE_AS;
		toolfunc = file_save_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 2:
	    if (vwin->role == GR_PLOT) {
		stockicon = GTK_STOCK_EXECUTE;
		toolfunc = gp_send_callback;
	    } else
		toolfunc = NULL;
	    break;	    
	case 3:
	    if (print_ok) {
#if defined(G_OS_WIN32) || defined(USE_GNOME)
		stockicon = GTK_STOCK_PRINT;
		toolfunc = window_print_callback;
#endif
	    } else
		toolfunc = NULL;
	    break;
	case 4:
	    if (run_ok) {
		stockicon = GTK_STOCK_EXECUTE;
		toolfunc = run_script_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 5:
	    stockicon = GTK_STOCK_COPY;
	    if (MULTI_COPY_ENABLED(vwin->role)) {
		toolfunc = choose_copy_format_callback;
	    } else {
		toolfunc = text_copy_callback;
	    }
	    break;
	case 6:
	    if (edit_ok) {
		stockicon = GTK_STOCK_PASTE;
		toolfunc = text_paste_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 7:
	    if (vwin->role != VIEW_MODELTABLE) {
		stockicon = GTK_STOCK_FIND;
		toolfunc = text_find_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 8:
	    if (edit_ok) {
		stockicon = GTK_STOCK_FIND_AND_REPLACE;
		toolfunc = text_replace_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 9:
	    if (edit_ok) {
		stockicon = GTK_STOCK_UNDO;
		toolfunc = text_undo_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 10:
	    if (run_ok) {
		stockicon = GTK_STOCK_HELP;
		toolfunc = activate_script_help;
	    } else
		toolfunc = NULL;
	    break;
	case 11:
	    if (vwin->role == VIEW_MODELTABLE) {
		stockicon = "STOCK_TEX";
		toolfunc = modeltable_tex_view;
	    } else
		toolfunc = NULL;
	    break;
	case 12:
	    stockicon = GTK_STOCK_CLOSE;
	    toolfunc = delete_file_viewer;
	    break;
	default:
	    break;
	}

	if (toolfunc == NULL) continue;

	toolstr = _(viewstrings[i]);

	button = gtk_image_new();
	gtk_image_set_from_stock(GTK_IMAGE(button), stockicon, 
				 GTK_ICON_SIZE_MENU);
        gtk_toolbar_append_item(GTK_TOOLBAR(viewbar),
				NULL, toolstr, NULL,
				button, toolfunc, vwin);
    }
    gtk_widget_show(viewbar);
    gtk_widget_show(hbox);
}


/* ........................................................... */

static gchar *make_viewer_title (int role, const char *fname)
{
    gchar *title = NULL;

    switch (role) {
    case GUI_HELP: 
	title = g_strdup(_("gretl: help")); break;
    case CLI_HELP:
	title = g_strdup(_("gretl: command syntax")); break;
    case GUI_HELP_ENGLISH: 
	title = g_strdup("gretl: help"); break;
    case CLI_HELP_ENGLISH:
	title = g_strdup("gretl: command syntax"); break;
    case VIEW_LOG:
	title = g_strdup(_("gretl: command log")); break;
    case CONSOLE:
	title = g_strdup(_("gretl console")); break;
    case EDIT_SCRIPT:
    case VIEW_SCRIPT:	
    case VIEW_FILE:
	if (strstr(fname, "script_tmp") || strstr(fname, "session.inp"))
	    title = g_strdup(_("gretl: command script"));
	else {
	    gchar *p = strrchr(fname, SLASH);
	    title = g_strdup_printf("gretl: %s", p? p + 1 : fname);
	} 
	break;
    case EDIT_NOTES:
	title = g_strdup(_("gretl: session notes")); break;
    case GR_PLOT:
    case GR_BOX:
	title = g_strdup(_("gretl: edit plot commands")); break;
    case SCRIPT_OUT:
	title = g_strdup(_("gretl: script output")); break;
    case VIEW_DATA:
	title = g_strdup(_("gretl: display data")); break;
    default:
	break;
    }
    return title;
}

/* ........................................................... */

static gint script_changed (GtkWidget *w, windata_t *vwin)
{
    MARK_SCRIPT_CHANGED(vwin);
    return FALSE;
}

/* ........................................................... */

static void auto_save_script (windata_t *vwin)
{
    FILE *fp;
    char msg[MAXLEN];
    gchar *savestuff;

    if (strstr(vwin->fname, "script_tmp") || !strlen(vwin->fname)) {
	file_save(vwin, SAVE_SCRIPT, NULL);
	strcpy(vwin->fname, scriptfile);
    }

    if ((fp = fopen(vwin->fname, "w")) == NULL) {
	sprintf(msg, _("Couldn't write to %s"), vwin->fname);
	errbox(msg); 
	return;
    }
    savestuff = textview_get_text(GTK_TEXT_VIEW(vwin->w));
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);
    infobox(_("script saved"));
    MARK_SCRIPT_SAVED(vwin);
}

/* ........................................................... */

static gint query_save_script (GtkWidget *w, GdkEvent *event, windata_t *vwin)
{
    if (SCRIPT_IS_CHANGED(vwin)) {
	int resp = 
	    yes_no_dialog(_("gretl: script"), 
			  _("Save changes?"), 1);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	}
	if (resp == GRETL_YES) {
	    auto_save_script(vwin);
	}
    }
    return FALSE;
}

/* ........................................................... */

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

static GtkTextTagTable *gretl_tags_new (void)
{
    GtkTextTagTable *table;
    GtkTextTag *tag;

    table = gtk_text_tag_table_new(); 

    tag = gtk_text_tag_new("bluetext");
    g_object_set(tag, "foreground", "blue", NULL);
    gtk_text_tag_table_add(table, tag);
    
    tag = gtk_text_tag_new("redtext");
    g_object_set(tag, "foreground", "red", NULL);
    gtk_text_tag_table_add(table, tag);

    return table;
}

#ifndef USE_GTKSOURCEVIEW

static void correct_line_color (windata_t *vwin)
{
    GtkTextBuffer *buf;
    GtkTextIter start, end;
    gint linelen;
    gchar *txt;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gtk_text_buffer_get_iter_at_mark(buf, &end, 
				     gtk_text_buffer_get_insert(buf));
    linelen = gtk_text_iter_get_chars_in_line(&end);
    start = end;
    gtk_text_iter_backward_chars(&start, linelen);

    txt = gtk_text_buffer_get_text(buf, &start, &end, FALSE);

    if (*txt == '#') {
	gtk_text_buffer_apply_tag_by_name (buf, "bluetext",
					   &start, &end);
    }
    g_free(txt);
}

#endif /* not USE_GTKSOURCEVIEW */

/* ........................................................... */

static windata_t *common_viewer_new (int role, const char *title, 
				     gpointer data, int record)
{
    windata_t *vwin;

    if ((vwin = mymalloc(sizeof *vwin)) == NULL) return NULL;

    windata_init(vwin);
    vwin->role = role;
    vwin->data = data;
    vwin->dialog = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(vwin->dialog), title);

    if (record) {
	g_object_set_data(G_OBJECT(vwin->dialog), "object", data);
	winstack_add(vwin->dialog);
    }

    return vwin;
}

/* ........................................................... */

#ifdef USE_GTKSOURCEVIEW

static void create_source (windata_t *vwin, GtkSourceBuffer **buf, 
			   int hsize, int vsize, gboolean editable)
{
    GtkSourceLanguagesManager *lm;
    GtkSourceBuffer *sbuf;
    GtkSourceTagStyle *tagstyle;
    GdkColor blue;

    blue.green = blue.red = 0;
    blue.blue = 65535.0;

    lm = gtk_source_languages_manager_new ();
    tagstyle = gtk_source_tag_style_new ();
    
    sbuf = GTK_SOURCE_BUFFER(gtk_source_buffer_new(NULL));
    g_object_ref (lm);
    g_object_set_data_full (G_OBJECT (sbuf), "languages-manager",
			    lm, (GDestroyNotify) g_object_unref); 
    g_object_unref (lm); 

    tagstyle->mask = GTK_SOURCE_TAG_STYLE_USE_FOREGROUND;
    tagstyle->foreground = blue;
    g_object_set_data_full (G_OBJECT (sbuf), "tag-style",
			    tagstyle, 
			    (GDestroyNotify) gtk_source_tag_style_free); 
    gtk_source_buffer_set_bracket_match_style(sbuf, tagstyle);
    gtk_source_buffer_set_check_brackets(sbuf, TRUE);

    vwin->w = gtk_source_view_new_with_buffer(sbuf);
    *buf = sbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->w), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->w), fixed_font);
    hsize *= get_char_width(vwin->w);
    hsize += 48;
    gtk_window_set_default_size (GTK_WINDOW(vwin->dialog), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);
}

#endif /* USE_GTKSOURCEVIEW */

/* ........................................................... */

static void create_text (windata_t *vwin, GtkTextBuffer **buf, 
			 int hsize, int vsize, gboolean editable)
{
    static GtkTextTagTable *tags = NULL;
    GtkTextBuffer *tbuf; 

    if (tags == NULL) tags = gretl_tags_new();

    tbuf = gtk_text_buffer_new(tags);
    vwin->w = gtk_text_view_new_with_buffer(tbuf);
    *buf = tbuf;

    gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(vwin->w), GTK_WRAP_WORD);
    gtk_text_view_set_left_margin(GTK_TEXT_VIEW(vwin->w), 4);
    gtk_text_view_set_right_margin(GTK_TEXT_VIEW(vwin->w), 4);

    gtk_widget_modify_font(GTK_WIDGET(vwin->w), fixed_font);
    hsize *= get_char_width(vwin->w);
    hsize += 48;
    gtk_window_set_default_size (GTK_WINDOW(vwin->dialog), hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);
}

/* ........................................................... */

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width (GTK_CONTAINER(vwin->vbox), 4);
    gtk_box_set_spacing (GTK_BOX(vwin->vbox), 4);
#ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(vwin->dialog), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif
    gtk_container_add(GTK_CONTAINER(vwin->dialog), vwin->vbox);
}

/* ........................................................... */

static void dialog_table_setup (windata_t *vwin)
{
    GtkWidget *sw;

    sw = gtk_scrolled_window_new (NULL, NULL);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), 
		       sw, TRUE, TRUE, FALSE);
    gtk_scrolled_window_set_policy (GTK_SCROLLED_WINDOW (sw),
				    GTK_POLICY_AUTOMATIC,
				    GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type (GTK_SCROLLED_WINDOW (sw),
					 GTK_SHADOW_IN);
    gtk_container_add (GTK_CONTAINER(sw), vwin->w); 
    gtk_widget_show(vwin->w);
    gtk_widget_show(sw);
}

/* ........................................................... */

static void cursor_to_top (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w)); 
    GtkTextIter start;
    GtkTextMark *mark;

    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_place_cursor(buf, &start);
    mark = gtk_text_buffer_create_mark(buf, NULL, &start, FALSE);
    gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				 mark, 0.0, FALSE, 0, 0);
}

/* ........................................................... */

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			char *title, int role, 
			GtkItemFactoryEntry menu_items[]) 
{
    GtkWidget *close;
    GtkTextBuffer *tbuf;
    windata_t *vwin;

    vwin = common_viewer_new(role, title, NULL, 1);
    if (vwin == NULL) return NULL;

    viewer_box_config(vwin);

    if (menu_items != NULL) {
	set_up_viewer_menu(vwin->dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    } else if (role != IMPORT) {
	make_viewbar(vwin, 1);
    }

    create_text (vwin, &tbuf, hsize, vsize, FALSE);
    
    dialog_table_setup (vwin);

    /* arrange for clean-up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    /* close button */
    close = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), 
		       close, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(close), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(close);

    /* insert and then free the text buffer */
    gtk_text_buffer_set_text(tbuf, prn->buf, -1);
    gretl_print_destroy(prn);
    
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_view_key), vwin);

    g_signal_connect (G_OBJECT(vwin->w), "button_press_event", 
		      G_CALLBACK(catch_button_3), vwin->w);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    cursor_to_top(vwin);

    return vwin;
}

#define help_role(r) (r == CLI_HELP || \
                      r == GUI_HELP || \
                      r == CLI_HELP_ENGLISH || \
                      r == GUI_HELP_ENGLISH)

#ifdef USE_GTKSOURCEVIEW

static int 
gtk_source_buffer_load_file (GtkSourceBuffer *sbuf, 
			     const char *fname,
			     int role)
{
    FILE *fp;
    GtkTextIter iter;    
    char readbuf[MAXSTR], *chunk = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL) return 1;

    gtk_source_buffer_begin_not_undoable_action (sbuf);

    gtk_text_buffer_set_text(GTK_TEXT_BUFFER(sbuf), "", 0);
    gtk_text_buffer_get_iter_at_offset(GTK_TEXT_BUFFER(sbuf), &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf - 1, fp)) {
#ifdef ENABLE_NLS
	if (!g_utf8_validate(readbuf, sizeof readbuf, NULL)) {
	    gsize bytes;

	    chunk = g_locale_to_utf8(readbuf, -1, NULL, &bytes, NULL);
	} else {
	    chunk = readbuf;
	}
#else
	chunk = readbuf;
#endif
	gtk_text_buffer_insert(GTK_TEXT_BUFFER(sbuf), &iter, chunk, -1);
	memset(readbuf, 0, sizeof readbuf);
	if (chunk != NULL && chunk != readbuf) {
	    g_free(chunk);
	    chunk = NULL;
	}
    }
    fclose(fp);
	
    gtk_source_buffer_end_not_undoable_action(sbuf);

    gtk_text_buffer_set_modified(GTK_TEXT_BUFFER(sbuf), FALSE);

    /* move cursor to the beginning */
    gtk_text_buffer_get_start_iter(GTK_TEXT_BUFFER(sbuf), &iter);
    gtk_text_buffer_place_cursor(GTK_TEXT_BUFFER(sbuf), &iter);

    return 0;
}

static void source_buffer_insert_file (GtkSourceBuffer *sbuf, 
				       const char *filename,
				       int role)
{
    GtkSourceLanguagesManager *manager;    
    GtkSourceLanguage *language = NULL;
		
    manager = g_object_get_data(G_OBJECT (sbuf), "languages-manager");

    if (role == GR_PLOT) {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gnuplot");
    } else {
	language = 
	    gtk_source_languages_manager_get_language_from_mime_type 
	    (manager, "application/x-gretlsession");
    }

    if (language == NULL) {
	g_object_set(G_OBJECT(sbuf), "highlight", FALSE, NULL);
    } else {
	g_object_set(G_OBJECT(sbuf), "highlight", TRUE, NULL);
	gtk_source_buffer_set_language(sbuf, language);
    }

    gtk_source_buffer_load_file(sbuf, filename, role);
}

#endif

static void text_buffer_insert_file (GtkTextBuffer *tbuf, const char *fname, 
				     int role)
{
    FILE *fp;
    GtkTextIter iter;    
    int thiscolor, nextcolor;
    char readbuf[MAXSTR], *chunk = NULL;

    fp = fopen(fname, "r");
    if (fp == NULL) return;

    thiscolor = nextcolor = PLAIN_TEXT;

    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);

    memset(readbuf, 0, sizeof readbuf);

    while (fgets(readbuf, sizeof readbuf - 1, fp)) {
#ifdef ENABLE_NLS
	if (!g_utf8_validate(readbuf, sizeof readbuf, NULL)) {
	    gsize bytes;

	    chunk = g_locale_to_utf8(readbuf, -1, NULL, &bytes, NULL);
	} else chunk = readbuf;
#else
	chunk = readbuf;
#endif
	if (*chunk == '@') continue;
	if (*chunk == '?') 
	    thiscolor = (role == CONSOLE)? RED_TEXT : BLUE_TEXT;
	if (*chunk == '#') {
	    if (help_role(role)) {
		*chunk = ' ';
		nextcolor = RED_TEXT;
	    } else
		thiscolor = BLUE_TEXT;
	} else
	    nextcolor = PLAIN_TEXT;
	switch (thiscolor) {
	case PLAIN_TEXT:
	    gtk_text_buffer_insert(tbuf, &iter, chunk, -1);
	    break;
	case BLUE_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      chunk, -1,
						      "bluetext", NULL);
	    break;
	case RED_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      chunk, -1,
						      "redtext", NULL);
	    break;
	}
	thiscolor = nextcolor;
	memset(readbuf, 0, sizeof readbuf);
	if (chunk != NULL && chunk != readbuf) {
	    free(chunk);
	    chunk = NULL;
	}
    }
    fclose(fp);
}

/* ........................................................... */

windata_t *view_file (char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role, 
		      GtkItemFactoryEntry menu_items[]) 
{
    GtkTextBuffer *tbuf = NULL;
#ifdef USE_GTKSOURCEVIEW
    GtkSourceBuffer *sbuf = NULL;
#endif
    char *fname = NULL;
    FILE *fp = NULL;
    windata_t *vwin;
    gchar *title;
    int show_viewbar = (role != CONSOLE &&
			role != VIEW_DATA &&
			role != VIEW_FILE &&
			!help_role(role));
    int doing_script = (role == EDIT_SCRIPT ||
			role == VIEW_SCRIPT ||
			role == VIEW_LOG);

    fp = fopen(filename, "r");
    if (fp == NULL) {
	sprintf(errtext, _("Can't open %s for reading"), filename);
	errbox(errtext);
	return NULL;
    } else {
	fclose(fp);
    }

    title = make_viewer_title(role, filename);
    vwin = common_viewer_new(role, (title != NULL)? title : filename, 
			     NULL, !doing_script && role != CONSOLE);

    if (title != NULL) g_free(title);
    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);

    if (menu_items != NULL) {
	set_up_viewer_menu(vwin->dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    } else if (show_viewbar) { 
	make_viewbar(vwin, 0);
	show_viewbar = 0;
    }

#ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	create_source(vwin, &sbuf, hsize, vsize, editable);
	tbuf = GTK_TEXT_BUFFER(sbuf);
	vwin->sbuf = sbuf;
    } else {
	create_text(vwin, &tbuf, hsize, vsize, editable);
    }
#else
    create_text(vwin, &tbuf, hsize, vsize, editable);
#endif

    dialog_table_setup (vwin);

    /* special case: the gretl console */
    if (role == CONSOLE) {
	g_signal_connect(G_OBJECT(vwin->w), "button_release_event",
			 G_CALLBACK(console_mouse_handler), NULL);
	g_signal_connect(G_OBJECT(vwin->w), "key_press_event",
			 G_CALLBACK(console_key_handler), NULL);
    } 

    if (doing_script) {
	g_signal_connect(G_OBJECT(vwin->w), "button_release_event",
			 G_CALLBACK(edit_script_help), vwin);
    } 

    /* is the file to be deleted after viewing? */
    if (del_file) {
	fname = g_strdup(filename);
    }

    /* should we show a toolbar at the foot of the window? */
    if (show_viewbar) { 
	make_viewbar(vwin, 0);
    } else { 
	/* make a simple Close button instead */
	GtkWidget *close = gtk_button_new_with_label(_("Close"));

	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   close, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(close), "clicked", 
			 G_CALLBACK(delete_file_viewer), vwin);
	gtk_widget_show(close);
    }

#ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	source_buffer_insert_file(sbuf, filename, role);
    } else {
	text_buffer_insert_file(tbuf, filename, role);
    }
#else
    text_buffer_insert_file(tbuf, filename, role);
#endif

    /* grab the "changed" signal when editing a script */
    if (role == EDIT_SCRIPT) {
	g_signal_connect(G_OBJECT(tbuf), "changed", 
			 G_CALLBACK(script_changed), vwin);
    }

    /* catch some keystrokes */
    if (!editable) {
	g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
			 G_CALLBACK(catch_view_key), vwin);
    } else {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
	g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
			 G_CALLBACK(catch_edit_key), vwin);	
    } 

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);

    /* offer chance to save script on exit */
    if (role == EDIT_SCRIPT)
	g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
			 G_CALLBACK(query_save_script), vwin);

    /* clean up when dialog is destroyed */
    if (del_file) {
	g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
			 G_CALLBACK(delete_file), (gpointer) fname);
    }
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);
    cursor_to_top(vwin);

    return vwin;
}

/* ........................................................... */

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role) 
{
    GtkWidget *close;
    GtkTextBuffer *tbuf;
    windata_t *vwin;

    vwin = common_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) return NULL;

    viewer_box_config(vwin); 

    /* add a menu bar */
    make_viewbar(vwin, 0);

    create_text(vwin, &tbuf, hsize, vsize, TRUE);

    dialog_table_setup (vwin);
    
    /* insert the buffer text */
    if (*pbuf) gtk_text_buffer_set_text(tbuf, *pbuf, -1);

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);

    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_edit_key), vwin);	

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    /* close button */
    close = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), 
		       close, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(close), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(close);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);
    cursor_to_top(vwin);

    return vwin;
}

/* ........................................................... */

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title) 
{
    windata_t *vwin;
    GtkWidget *close;
    GtkTextBuffer *tbuf;

    vwin = common_viewer_new(VIEW_MODEL, title, pmod, 1);
    if (vwin == NULL) return 1;

    viewer_box_config(vwin);

    set_up_viewer_menu(vwin->dialog, vwin, model_items);
    add_vars_to_plot_menu(vwin);
    add_dummies_to_plot_menu(vwin);
    g_signal_connect(G_OBJECT(vwin->mbar), "button_press_event", 
		       G_CALLBACK(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    create_text(vwin, &tbuf, hsize, vsize, FALSE);

    dialog_table_setup (vwin);

    /* close button */
    close = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), close, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(close), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(close);

    /* insert and then free the model buffer */
    gtk_text_buffer_set_text(tbuf, prn->buf, strlen(prn->buf));
    gretl_print_destroy(prn);

    if (pmod->ci != NLS) copylist(&default_list, pmod->list);

    /* attach shortcuts */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_view_key), vwin);

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(delete_unnamed_model), 
		     vwin->data);
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), 
		     vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show_all(vwin->dialog);
    cursor_to_top(vwin);

    return 0;
}

/* ........................................................... */

void flip (GtkItemFactory *ifac, const char *path, gboolean s)
{
    if (ifac != NULL) {
	GtkWidget *w = gtk_item_factory_get_item(ifac, path);

	if (w != NULL) {
	    gtk_widget_set_sensitive(w, s);
	} else {
	    fprintf(stderr, I_("Failed to flip state of \"%s\"\n"), path);
	}
    }
}

/* ........................................................... */

static void model_equation_copy_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/LaTeX/View/Equation", s);
    flip(ifac, "/LaTeX/Save/Equation", s);
    flip(ifac, "/LaTeX/Copy/Equation", s);
}

/* ........................................................... */

static void model_arch_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/ARCH", s);
}

/* ........................................................... */

static void model_panel_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/panel diagnostics", s);
}

/* ........................................................... */

static void model_ml_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Model data/Add to data set/log likelihood", s);
}

/* ........................................................... */

static void model_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/non-linearity (squares)", s);
    flip(ifac, "/Tests/non-linearity (logs)", s);
    flip(ifac, "/Tests/autocorrelation", s);
    flip(ifac, "/Tests/heteroskedasticity", s);
    flip(ifac, "/Tests/Chow test", s);
    flip(ifac, "/Tests/CUSUM test", s);
    flip(ifac, "/Tests/ARCH", s);
    flip(ifac, "/Tests/normality of residual", s);
    flip(ifac, "/Graphs", s);
    flip(ifac, "/Model data/Display actual, fitted, residual", s);
    flip(ifac, "/Model data/Forecasts with standard errors", s);
    flip(ifac, "/Model data/Add to data set/residuals", s);
    flip(ifac, "/Model data/Add to data set/squared residuals", s);
    flip(ifac, "/Model data/Add to data set/error sum of squares", s);
    flip(ifac, "/Model data/Add to data set/standard error of residuals", s);
    flip(ifac, "/Model data/Add to data set/R-squared", s);
    flip(ifac, "/Model data/Add to data set/T*R-squared", s);    
}

/* ........................................................... */

static void ols_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/non-linearity (squares)", s);
    flip(ifac, "/Tests/non-linearity (logs)", s);
    flip(ifac, "/Tests/autocorrelation", s);
    flip(ifac, "/Tests/heteroskedasticity", s);
    flip(ifac, "/Tests/Chow test", s);
    flip(ifac, "/Tests/CUSUM test", s);
    flip(ifac, "/Tests/ARCH", s);
    flip(ifac, "/Tests/influential observations", s);
}

/* ........................................................... */

static void lad_menu_mod (GtkItemFactory *ifac)
{
    flip(ifac, "/Tests", FALSE);
    flip(ifac, "/Model data/Forecasts with standard errors", FALSE);
    flip(ifac, "/Model data/coefficient covariance matrix", FALSE);
    flip(ifac, "/Model data/Add to data set/R-squared", FALSE);
    flip(ifac, "/Model data/Add to data set/T*R-squared", FALSE);
}

/* ........................................................... */

static void nls_menu_mod (GtkItemFactory *ifac)
{
    flip(ifac, "/Tests/omit variables", FALSE);
    flip(ifac, "/Tests/add variables", FALSE);
    flip(ifac, "/Tests/sum of coefficients", FALSE);
    flip(ifac, "/Tests/Ramsey's RESET", FALSE);
    flip(ifac, "/Model data/Forecasts with standard errors", FALSE);
    flip(ifac, "/Model data/Confidence intervals for coefficients", FALSE);
}

/* ........................................................... */

static void latex_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/LaTeX", s);
}

/* ........................................................... */

static void model_save_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/File/Save to session as icon", s);
    flip(ifac, "/File/Save as icon and close", s);
}

/* ........................................................... */

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[])
{
    gint n_items = 0;

    while (items[n_items].path != NULL) n_items++;

    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", 
				      NULL);
#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(vwin->ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items(vwin->ifac, n_items, items, vwin);
    vwin->mbar = gtk_item_factory_get_widget(vwin->ifac, "<main>");

    if (vwin->role == VIEW_MODEL && vwin->data != NULL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	latex_menu_state(vwin->ifac, !pmod->errcode);

	model_equation_copy_state(vwin->ifac, 
				  !pmod->errcode && pmod->ci != NLS);

	model_panel_menu_state(vwin->ifac, pmod->ci == POOLED);

	ols_menu_state(vwin->ifac, pmod->ci == OLS || pmod->ci == POOLED);

	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    model_menu_state(vwin->ifac, FALSE);
	    model_ml_menu_state(vwin->ifac, TRUE);
	} else {
	    model_ml_menu_state(vwin->ifac, FALSE);
	}

	if (pmod->name) model_save_state(vwin->ifac, FALSE);

	if (pmod->ci == LAD) lad_menu_mod(vwin->ifac);
	else if (pmod->ci == NLS) nls_menu_mod(vwin->ifac);

	if (dataset_is_panel(datainfo)) {
	    model_arch_menu_state(vwin->ifac, FALSE);
	}
    }
}

/* .................................................................. */

static void add_vars_to_plot_menu (windata_t *vwin)
{
    int i, j, varstart;
    GtkItemFactoryEntry varitem;
    const gchar *mpath[] = {
	N_("/Graphs/residual plot"), 
	N_("/Graphs/fitted, actual plot")
    };
    MODEL *pmod = vwin->data;

    for (i=0; i<2; i++) {
	varitem.accelerator = NULL;
	varitem.callback_action = 0; 
	varitem.item_type = NULL;
	if (dataset_is_time_series(datainfo)) {
	    varitem.path = 
		g_strdup_printf(_("%s/against time"), mpath[i]);
	} else {
	    varitem.path = 
		g_strdup_printf(_("%s/by observation number"), mpath[i]);	
	}
	varitem.callback = (i==0)? resid_plot : fit_actual_plot;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	varstart = (i == 0)? 1 : 2;

	/* put the indep vars on the menu list */
	for (j=varstart; pmod->ci != NLS && j<=pmod->list[0]; j++) {
	    if (pmod->list[j] == 0) continue;
	    if (!strcmp(datainfo->varname[pmod->list[j]], "time")) 
		continue;
	    varitem.accelerator = NULL;
	    varitem.callback_action = pmod->list[j]; 
	    varitem.item_type = NULL;
	    varitem.path = 
		g_strdup_printf(_("%s/against %s"), mpath[i],
				datainfo->varname[pmod->list[j]]);
	    varitem.callback = (i==0)? resid_plot : fit_actual_plot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}

	/* if the model has two indepdent vars, offer a 3-D fitted
	   versus actual plot */
	if (i == 1 && pmod->ifc && pmod->ncoeff == 3) {
	    varitem.accelerator = NULL;
	    varitem.callback_action = 0;
	    varitem.item_type = NULL;
	    varitem.path =
		g_strdup_printf(_("%s/against %s and %s"),
				mpath[i], 
				datainfo->varname[pmod->list[3]],
				datainfo->varname[pmod->list[4]]);
	    varitem.callback = fit_actual_splot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}	
    }
}

/* .................................................................. */

static void plot_dummy_call (gpointer data, guint v, GtkWidget *widget)
{
    GtkCheckMenuItem *item = GTK_CHECK_MENU_ITEM(widget);
    windata_t *vwin = (windata_t *) data;

    if (item->active) vwin->active_var = v; 
}

/* .................................................................. */

static void add_dummies_to_plot_menu (windata_t *vwin)
{
    int i, dums = 0;
    GtkItemFactoryEntry dumitem;
    MODEL *pmod = vwin->data;
    const gchar *mpath[] = {
	N_("/Graphs/dumsep"), 
	N_("/Graphs/Separation")
    };
    gchar *radiopath = NULL;

    dumitem.path = NULL;

    /* put the dummy independent vars on the menu list */
    for (i=2; i<pmod->list[0]; i++) {

	if (pmod->list[i] == 0) continue;
	if (pmod->list[i] == LISTSEP) break;

	if (!isdummy(Z[pmod->list[i]], datainfo->t1, datainfo->t2)) {
	    continue;
	}

	if (!dums) { /* add separator, branch and "none" */
	    dumitem.path = mymalloc(64);
	    sprintf(dumitem.path, _("%s"), mpath[0]);
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Separator>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* menu branch */
	    sprintf(dumitem.path, _("%s"), mpath[1]);
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Branch>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* "none" option */
	    sprintf(dumitem.path, _("%s/none"), mpath[1]);
	    radiopath = g_strdup(dumitem.path);
	    dumitem.callback = plot_dummy_call;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<RadioItem>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    dums = 1;
	} 

	dumitem.callback_action = pmod->list[i]; 
	sprintf(dumitem.path, _("%s/by %s"), mpath[1],  
		datainfo->varname[pmod->list[i]]);
	dumitem.callback = plot_dummy_call;	    
	dumitem.accelerator = NULL;
	dumitem.item_type = radiopath;
	gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);

    }

    free(dumitem.path);
    free(radiopath);
}

/* ........................................................... */

#define ALLOW_MODEL_DATASETS

static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data)
{
    windata_t *mwin = (windata_t *) data;
    MODEL *pmod = mwin->data;
    extern int quiet_sample_check (MODEL *pmod); /* library.c */
    int s, ok = 1;

    if (Z == NULL) {
	flip(mwin->ifac, "/File/Save to session as icon", FALSE);
	flip(mwin->ifac, "/File/Save as icon and close", FALSE);
	flip(mwin->ifac, "/Edit/Copy all", FALSE);
	flip(mwin->ifac, "/Model data", FALSE);
	flip(mwin->ifac, "/Tests", FALSE);
	flip(mwin->ifac, "/Graphs", FALSE);
	flip(mwin->ifac, "/Model data", FALSE);
	flip(mwin->ifac, "/LaTeX", FALSE);
	return FALSE;
    }

    if (quiet_sample_check(pmod)) ok = 0;

    s = GTK_WIDGET_IS_SENSITIVE
	(gtk_item_factory_get_item(mwin->ifac, "/Tests"));
    if ((s && ok) || (!s && !ok)) return FALSE;
    s = !s;

#ifdef ALLOW_MODEL_DATASETS
    if (!ok && dataset_added_to_model(pmod)) {
	flip(mwin->ifac, "/Tests", s);
	flip(mwin->ifac, "/Model data/Display actual, fitted, residual", s);
	if (pmod->ci != LAD) {
	    flip(mwin->ifac, "/Model data/Forecasts with standard errors", s);
	}
	flip(mwin->ifac, "/Model data/Confidence intervals for coefficients", s);
	flip(mwin->ifac, "/Model data/Add to data set/fitted values", s);
	flip(mwin->ifac, "/Model data/Add to data set/residuals", s);
	flip(mwin->ifac, "/Model data/Add to data set/squared residuals", s);
	flip(mwin->ifac, "/Model data/Define new variable...", s);
	infobox(get_gretl_errmsg());
	return FALSE;
    }
#endif

    flip(mwin->ifac, "/Tests", s);
    flip(mwin->ifac, "/Graphs", s);
    flip(mwin->ifac, "/Model data/Display actual, fitted, residual", s);
    if (pmod->ci != LAD) {
	flip(mwin->ifac, "/Model data/Forecasts with standard errors", s);
    }
    flip(mwin->ifac, "/Model data/Confidence intervals for coefficients", s);
    flip(mwin->ifac, "/Model data/Add to data set/fitted values", s);
    flip(mwin->ifac, "/Model data/Add to data set/residuals", s);
    flip(mwin->ifac, "/Model data/Add to data set/squared residuals", s);
    flip(mwin->ifac, "/Model data/Define new variable...", s);

    if (!ok) {
	infobox(get_gretl_errmsg());
    }    

    return FALSE;
}

/* ........................................................... */

#if defined(G_OS_WIN32)

static void msgbox (const char *msg, int err)
{
    gchar *trmsg = NULL;

    if (nls_on) {
	gint wrote;

	trmsg = g_locale_from_utf8 (msg, -1, NULL, &wrote, NULL);
    } 

    if (err) 
	MessageBox(NULL, (nls_on)? trmsg : msg, "gretl", 
		   MB_OK | MB_ICONERROR);
    else
	MessageBox(NULL, (nls_on)? trmsg : msg, "gretl", 
		   MB_OK | MB_ICONINFORMATION);

    if (nls_on) g_free(trmsg);
}

#else /* GTK native */

static void msgbox (const char *msg, int err)
{
    GtkWidget *dialog;

    dialog = gtk_message_dialog_new (NULL, /* GTK_WINDOW(mdata->w), */
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     (err)? GTK_MESSAGE_ERROR : GTK_MESSAGE_INFO,
				     GTK_BUTTONS_CLOSE,
				     msg);
    gtk_dialog_run (GTK_DIALOG (dialog));
    gtk_widget_destroy (dialog);
}

#endif

/* ........................................................... */

void errbox (const char *msg) 
{
    msgbox(msg, 1);
}

/* ........................................................... */

void infobox (const char *msg) 
{
    msgbox(msg, 0);
}

/* ........................................................... */

int validate_varname (const char *varname)
{
    int i, n = strlen(varname);
    char namebit[9];
    unsigned char c;

    *namebit = 0;
    
    if (n > 8) {
	strncat(namebit, varname, 8);
	sprintf(errtext, _("Variable name %s... is too long\n"
	       "(the max is 8 characters)"), namebit);
	errbox(errtext);
	return 1;
    }
    if (!(isalpha(*varname))) {
	sprintf(errtext, _("First char of name ('%c') is bad\n"
	       "(first must be alphabetical)"), *varname);
	errbox(errtext);
	return 1;
    }
    for (i=1; i<n; i++) {
	c = (unsigned char) varname[i];
	
	if ((!(isalpha(c)) && !(isdigit(c)) && c != '_') || c > 127) {
	    sprintf(errtext, _("Name contains an illegal char (in place %d)\n"
		    "Use only unaccented letters, digits and underscore"), i + 1);
	    errbox(errtext);
	    return 1;
	}
    }
    return 0;
}	

/* .................................................................. */

#if defined(G_OS_WIN32)

int prn_to_clipboard (PRN *prn, int copycode)
{
    return win_copy_buf(prn->buf, copycode, 0);
}

#elif defined(ENABLE_NLS)

int prn_to_clipboard (PRN *prn, int copycode)
{
    if (prn->buf == NULL) return 0;

    if (clipboard_buf) g_free(clipboard_buf);
    clipboard_buf = NULL;

    if (copycode == COPY_TEXT) { /* need to convert from utf8 */
	gchar *trbuf;
	gsize bytes;
	
	trbuf = g_locale_from_utf8(prn->buf, -1, NULL, &bytes, NULL);
	if (bytes > 0) {
	    clipboard_buf = mymalloc(bytes + 1);
	    if (clipboard_buf == NULL) {
		g_free(trbuf);
		return 1;
	    }
	    memcpy(clipboard_buf, trbuf, bytes + 1);
	    g_free(trbuf);
	}
    } else { /* copying TeX, RTF or CSV */
	size_t len;

	len = strlen(prn->buf);
	fprintf(stderr, "Copying to clipboard, %d bytes\n", (int) len);
	clipboard_buf = mymalloc(len + 1);
	if (clipboard_buf == NULL) return 1;
	memcpy(clipboard_buf, prn->buf, len + 1);
    }

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY, 
			    GDK_CURRENT_TIME);
    return 0;
}

#else /* plain GTK, no NLS */

int prn_to_clipboard (PRN *prn, int copycode)
{
    size_t len;
    
    if (prn->buf == NULL) return 0;
    len = strlen(prn->buf);
    if (len == 0) return 0;

    if (clipboard_buf) g_free(clipboard_buf);
    clipboard_buf = mymalloc(len + 1);
    if (clipboard_buf == NULL) return 1;

    memcpy(clipboard_buf, prn->buf, len + 1);

    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY,
			    GDK_CURRENT_TIME);
    return 0;
}

#endif /* switch for prn_to_clipboard */

/* .................................................................. */

#define SPECIAL_COPY(h) (h == COPY_LATEX || h == COPY_RTF)

void text_copy (gpointer data, guint how, GtkWidget *widget) 
{
    windata_t *vwin = (windata_t *) data;
    gchar *msg = NULL;
    PRN *prn;

    /* descriptive statistics */
    if ((vwin->role == SUMMARY || vwin->role == VAR_SUMMARY)
	&& SPECIAL_COPY(how)) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) vwin->data;
	
	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) {
	    texprint_summary(summ, datainfo, prn);
	} else if (how == COPY_RTF) { 
	    rtfprint_summary(summ, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* correlation matrix */
    else if (vwin->role == CORR && SPECIAL_COPY(how)) {
	CORRMAT *corr = (CORRMAT *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_corrmat(corr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_corrmat(corr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* display for fitted, actual, resid */
    else if (vwin->role == FCAST && SPECIAL_COPY(how)) {
	FITRESID *fr = (FITRESID *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_fit_resid(fr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_fit_resid(fr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }   

    /* forecasts with standard errors */
    else if (vwin->role == FCASTERR && SPECIAL_COPY(how)) {
	FITRESID *fr = (FITRESID *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_fcast_with_errs(fr, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_fcast_with_errs(fr, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }  

    /* coefficient confidence intervals */
    else if (vwin->role == COEFFINT && SPECIAL_COPY(how)) {
	CONFINT *cf = (CONFINT *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_confints(cf, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_confints(cf, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }  

    /* coefficient covariance matrix */
    else if (vwin->role == COVAR && SPECIAL_COPY(how)) {
	VCV *vcv = (VCV *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    texprint_vcv(vcv, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    rtfprint_vcv(vcv, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }      

    /* multiple-precision OLS */
    else if (vwin->role == MPOLS && SPECIAL_COPY(how)) {
	mp_results *mpvals = (mp_results *) vwin->data;

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) { 
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    print_mpols_results (mpvals, datainfo, prn);
	} 
	else if (how == COPY_RTF) { 
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    print_mpols_results (mpvals, datainfo, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* or it's a model window we're copying from? */
    else if (vwin->role == VIEW_MODEL &&
	(how == COPY_RTF || how == COPY_LATEX ||
	 how == COPY_LATEX_EQUATION)) {
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->errcode) { 
	    errbox("Couldn't format model");
	    return;
	}

	if (bufopen(&prn)) return;

	if (how == COPY_RTF) {
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    printmodel(pmod, datainfo, prn);
	}
	else if (how == COPY_LATEX) {
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    printmodel (pmod, datainfo, prn);
	}
	else if (how == COPY_LATEX_EQUATION) {
	    tex_print_equation(pmod, datainfo, 0, prn);
	}

	prn_to_clipboard(prn, how);
	gretl_print_destroy(prn);
    }

    /* or from the  model table? */
    else if (vwin->role == VIEW_MODELTABLE && SPECIAL_COPY(how)) {
	if (how == COPY_LATEX) {
	    tex_print_model_table(NULL, 0, NULL);
	}
	else if (how == COPY_RTF) {
	    rtf_print_model_table();
	} 
    }

    /* copying plain text from window */
    else if (how == COPY_TEXT || how == COPY_SELECTION) {
	GtkTextBuffer *textbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	PRN textprn;

	if (gtk_text_buffer_get_selection_bounds(textbuf, NULL, NULL)) {
	    /* there is a selection in place */
	    GtkTextIter selstart, selend;
	    gchar *selbuf;

	    gtk_text_buffer_get_selection_bounds(textbuf, &selstart, &selend);
	    selbuf = gtk_text_buffer_get_text(textbuf, &selstart, &selend, FALSE);
	    gretl_print_attach_buffer(&textprn, selbuf);
	    prn_to_clipboard(&textprn, COPY_TEXT);
	    g_free(selbuf);
	    infobox(_("Copied selection to clipboard"));
	    return;
	} else {
	    /* no selection: copy everything */
	    gretl_print_attach_buffer(&textprn,
				      textview_get_text(GTK_TEXT_VIEW(vwin->w))); 
	    prn_to_clipboard(&textprn, COPY_TEXT);
	    g_free(textprn.buf);
	}
    }

    msg = g_strdup_printf(_("Copied contents of window as %s"),
			  (how == COPY_LATEX)? "LaTeX" :
			  (how == COPY_RTF)? "RTF" : _("plain text"));
    infobox(msg);
    g_free(msg);
}

/* .................................................................. */

#if defined(G_OS_WIN32) || defined (USE_GNOME)

void window_print (windata_t *vwin, guint u, GtkWidget *widget) 
{
    char *buf, *selbuf = NULL;
    GtkTextView *tedit = GTK_TEXT_VIEW(vwin->w);
    GtkTextBuffer *tbuf = gtk_text_view_get_buffer(tedit);
    GtkTextIter start, end;

    buf = textview_get_text(tedit);

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end)) {
	selbuf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
    }

    winprint(buf, selbuf);
}

#endif

/* .................................................................. */

void text_undo (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old = NULL;

#ifdef USE_GTKSOURCEVIEW
    if (vwin->sbuf != NULL) {
	if (gtk_source_buffer_can_undo(vwin->sbuf)) {
	    gtk_source_buffer_undo(vwin->sbuf);
	} else {
	    errbox(_("No undo information available"));
	}
	return;
    }
#endif
    
    old = g_object_steal_data(G_OBJECT(vwin->w), "undo");

    if (old == NULL) {
	errbox(_("No undo information available"));
    } else {
	GtkTextBuffer *buf;
	GtkTextIter start, end;
	GtkTextMark *ins;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	ins = gtk_text_buffer_get_insert(buf);

	gtk_text_buffer_get_start_iter(buf, &start);
	gtk_text_buffer_get_end_iter(buf, &end);
	gtk_text_buffer_delete(buf, &start, &end);

	gtk_text_buffer_insert(buf, &start, old, strlen(old));
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				     ins, 0.0, TRUE, 0.1, 0.0);
	g_free(old);
    }
}

/* .................................................................. */

void text_paste (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old;
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    gchar *undo_buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));

    old = g_object_get_data(G_OBJECT(vwin->w), "undo");
    g_free(old);

    g_object_set_data(G_OBJECT(vwin->w), "undo", undo_buf);

    gtk_text_buffer_paste_clipboard(buf, gtk_clipboard_get(GDK_NONE),
				    NULL, TRUE);
}

/* ......................................................... */

gint popup_menu_handler (GtkWidget *widget, GdkEvent *event,
			 gpointer data)
{
    GdkModifierType mods;

    gdk_window_get_pointer(widget->window, NULL, NULL, &mods);
    
    if (mods & GDK_BUTTON3_MASK && event->type == GDK_BUTTON_PRESS) {
	GdkEventButton *bevent = (GdkEventButton *) event; 
	gtk_menu_popup (GTK_MENU(data), NULL, NULL, NULL, NULL,
			bevent->button, bevent->time);
	return TRUE;
    }
    return FALSE;
}

/* .................................................................. */

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback, gpointer data)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(label);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    gtk_widget_show(item);
    g_signal_connect (G_OBJECT(item), "activate",
		      G_CALLBACK(callback), data);
}

/* .................................................................. */

#ifndef G_OS_WIN32
# include <dlfcn.h>
#endif

int gui_open_plugin (const char *plugin, void **handle)
{
    char pluginpath[MAXLEN];

    strcpy(pluginpath, fetch_gretl_lib_path());

#ifdef G_OS_WIN32
    append_dir(pluginpath, "plugins");
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".dll");
    *handle = LoadLibrary(pluginpath);
    if (*handle == NULL) {
	sprintf(errtext, _("Couldn't load plugin %s"), pluginpath);
	errbox(errtext);
	return 1;
    }
#else
    strcat(pluginpath, plugin);
    strcat(pluginpath, ".so");
    *handle = dlopen(pluginpath, RTLD_LAZY);
    if (*handle == NULL) {
	sprintf(errtext, _("Failed to load plugin: %s"), pluginpath);
	errbox(errtext);
	return 1;
    } 
#endif 
    return 0;
}

/* .................................................................. */

void text_set_cursor (GtkWidget *w, GdkCursorType cspec)
{
    GdkWindow *win = gtk_text_view_get_window(GTK_TEXT_VIEW(w),
                                              GTK_TEXT_WINDOW_TEXT);

    if (cspec == 0) {
	gdk_window_set_cursor(win, NULL);
    } else {
	GdkCursor *cursor = gdk_cursor_new(cspec);

	gdk_window_set_cursor(win, cursor);
	gdk_cursor_destroy(cursor);
    } 
}

/* .................................................................. */

gint get_char_width (GtkWidget *widget)
{
    PangoLayout *pl;
    PangoContext *pc;
    GtkRcStyle *style;
    int width;

    pc = gtk_widget_get_pango_context(widget);
    style = gtk_widget_get_modifier_style(widget);
    pango_context_set_font_description(pc, style->font_desc);

    pl = pango_layout_new(pc);
    pango_layout_set_text(pl, "X", 1);
    pango_layout_get_pixel_size(pl, &width, NULL);

    g_object_unref(G_OBJECT(pl));

    return width;
}

/* .................................................................. */

gchar *textview_get_text (GtkTextView *view)
{
    GtkTextBuffer *tbuf;
    GtkTextIter start, end;

    tbuf = gtk_text_view_get_buffer(view);
    gtk_text_buffer_get_start_iter(tbuf, &start);
    gtk_text_buffer_get_end_iter(tbuf, &end);

    return gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);
}

/* .................................................................. */

int build_path (const char *dir, const char *fname, char *path, 
		const char *ext)
{
    size_t len;
    *path = 0;

    if (dir == NULL || fname == NULL || path == NULL) return 1;

    strcat(path, dir);
    len = strlen(path);
    if (len == 0) return 1;

    if (path[len - 1] == '/' || path[len - 1] == '\\') {
        /* dir is already properly terminated */
        strcat(path, fname);
    } else {
        /* otherwise put a separator in */
        strcat(path, SLASHSTR);
        strcat(path, fname);
    }

    if (ext != NULL) strcat(path, ext);

    return 0;
}


