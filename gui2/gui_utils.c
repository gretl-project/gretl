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

#ifdef G_OS_WIN32
# include <windows.h>
#endif

char *storelist = NULL;
GtkWidget *active_edit_id = NULL;
GtkWidget *active_edit_name = NULL;

extern int session_saved;
extern GtkWidget *mysheet;
extern char *space_to_score (char *str);

#define SCRIPT_CHANGED(w) w->active_var = 1
#define SCRIPT_SAVED(w) w->active_var = 0
#define SCRIPT_IS_CHANGED(w) w->active_var == 1

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static void file_viewer_save (GtkWidget *widget, windata_t *vwin);
static gint query_save_script (GtkWidget *w, GdkEvent *event, windata_t *vwin);
static void add_vars_to_plot_menu (windata_t *vwin);
static void add_dummies_to_plot_menu (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);

extern void do_coeff_intervals (gpointer data, guint i, GtkWidget *w);
extern void save_plot (char *fname, GPT_SPEC *plot);
extern gboolean console_handler (GtkWidget *w, GdkEventKey *key, 
				 gpointer user_data);
extern void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w);


GtkItemFactoryEntry model_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, NULL },
    { N_("/File/Save to session as icon"), NULL, remember_model, 0, NULL },
    { N_("/File/Save as icon and close"), NULL, remember_model, 1, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/Copy all/as plain _text"), NULL, text_copy, COPY_TEXT, NULL },
    { N_("/Edit/Copy all/as _LaTeX"), NULL, text_copy, COPY_LATEX, NULL },
    { N_("/Edit/Copy all/as _RTF"), NULL, text_copy, COPY_RTF, NULL },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>" },    
    { N_("/Tests/omit variables"), NULL, selector_callback, OMIT, NULL },
    { N_("/Tests/add variables"), NULL, selector_callback, ADD, NULL },
    { N_("/Tests/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/non-linearity (squares)"), NULL, do_lmtest, AUX_SQ, NULL },
    { N_("/Tests/non-linearity (logs)"), NULL, do_lmtest, AUX_LOG, NULL },
    { N_("/Tests/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/autocorrelation"), NULL, model_test_callback, LMTEST, NULL },
    { N_("/Tests/heteroskedasticity"), NULL, do_lmtest, AUX_WHITE, NULL },
    { N_("/Tests/Chow test"), NULL, model_test_callback, CHOW, NULL },
    { N_("/Tests/CUSUM test"), NULL, do_cusum, 0, NULL },
    { N_("/Tests/ARCH"), NULL, model_test_callback, ARCH, NULL },
    { N_("/Tests/normality of residual"), NULL, do_resid_freq, 0, NULL },
    { N_("/Tests/panel diagnostics"), NULL, do_panel_diagnostics, 0, NULL },
    { N_("/_Graphs"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/Graphs/residual plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/Graphs/fitted, actual plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/_Model data"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model data/Display actual, fitted, residual"), NULL, 
      display_fit_resid, 0, NULL },
    { N_("/Model data/Forecasts with standard errors"), NULL, 
      model_test_callback, FCAST, NULL },
    { N_("/Model data/Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL },
    { N_("/Model data/Add to data set"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model data/Add to data set/fitted values"), NULL, 
      fit_resid_callback, 1, NULL },
    { N_("/Model data/Add to data set/residuals"), NULL, 
      fit_resid_callback, 0, NULL },
    { N_("/Model data/Add to data set/squared residuals"), NULL, 
      fit_resid_callback, 2, NULL },
    { N_("/Model data/Add to data set/error sum of squares"), NULL, 
      model_stat_callback, ESS, NULL },
    { N_("/Model data/Add to data set/standard error of residuals"), NULL, 
      model_stat_callback, SIGMA, NULL },
    { N_("/Model data/Add to data set/R-squared"), NULL, 
      model_stat_callback, R2, NULL },
    { N_("/Model data/Add to data set/T*R-squared"), NULL, 
      model_stat_callback, TR2, NULL },
    { N_("/Model data/Add to data set/log likelihood"), NULL, 
      model_stat_callback, LNL, NULL },
    { N_("/Model data/Add to data set/degrees of freedom"), NULL, 
      model_stat_callback, DF, NULL },
    { N_("/Model data/coefficient covariance matrix"), NULL, 
      do_outcovmx, 0, NULL },
    { N_("/Model data/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Model data/Define new variable..."), NULL, model_test_callback, 
      MODEL_GENR, NULL },
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/_View"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/View/_Tabular"), NULL, view_latex, 0, NULL },
    { N_("/LaTeX/View/_Equation"), NULL, view_latex, 1, NULL },
    { N_("/LaTeX/_Save"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Save/_Tabular"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Save/Tabular/as _document"), NULL, file_save, 
      SAVE_TEX_TAB, NULL },
    { N_("/LaTeX/Save/Tabular/as _fragment"), NULL, file_save, 
      SAVE_TEX_TAB_FRAG, NULL },
    { N_("/LaTeX/Save/_Equation"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Save/Equation/as _document"), NULL, file_save, 
      SAVE_TEX_EQ, NULL },
    { N_("/LaTeX/Save/Equation/as _fragment"), NULL, file_save, 
      SAVE_TEX_EQ_FRAG, NULL },
    { N_("/LaTeX/_Copy"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Copy/_Tabular"), NULL, text_copy, COPY_LATEX, NULL },
    { N_("/LaTeX/Copy/_Equation"), NULL, text_copy, COPY_LATEX_EQUATION, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

GtkItemFactoryEntry edit_items[] = {
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif    
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, text_copy, COPY_TEXT, NULL },
    { N_("/Edit/_Paste"), NULL, text_paste, 0, NULL },
    { N_("/Edit/_Replace..."), NULL, text_replace, 0, NULL },
    { N_("/Edit/_Undo"), NULL, text_undo, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

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
    char buf[8192];
    size_t n;
   
    if ((srcfd = fopen(src, "rb")) == NULL) {
	return 1; 
    }
    if ((destfd = fopen(dest, "wb")) == NULL) {
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
    STACK_DESTROY
};

static void winstack (int code, GtkWidget *w)
{
    static int n_windows;
    static GtkWidget **wstack;
    int i;

    switch (code) {
    case STACK_DESTROY:	
	for (i=0; i<n_windows; i++) 
	    if (wstack[i] != NULL) 
		gtk_widget_destroy(wstack[i]);
	free(wstack);
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
    default:
	break;
    }
}

void winstack_init (void)
{
    winstack(STACK_INIT, NULL);
}
    
void winstack_destroy (void)
{
    winstack(STACK_DESTROY, NULL);
}

static void winstack_add (GtkWidget *w)
{
    winstack(STACK_ADD, w);
}

static void winstack_remove (GtkWidget *w)
{
    winstack(STACK_REMOVE, w);
}

/* ........................................................... */

static void delete_file (GtkWidget *widget, char *fle) 
{
    remove(fle);
    g_free(fle);
}

/* ........................................................... */

static void delete_file_viewer (GtkWidget *widget, gpointer data) 
{
    windata_t *vwin = (windata_t *) data;

    if (vwin->role == EDIT_SCRIPT && SCRIPT_IS_CHANGED(vwin)) {
	gint resp;

	resp = query_save_script(NULL, NULL, vwin);
	if (!resp) gtk_widget_destroy(vwin->dialog);
    } else 
	gtk_widget_destroy(vwin->dialog);
}

/* ........................................................... */

void delete_model (GtkWidget *widget, gpointer data) 
{
    MODEL *pmod = (MODEL *) data;
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

gint catch_key (GtkWidget *w, GdkEventKey *key)
{
    
    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    }
    else if (key->keyval == GDK_s) {
	windata_t *vwin = g_object_get_data(G_OBJECT(w), "ddata");

	if (Z != NULL && vwin != NULL && vwin->role == VIEW_MODEL)
	    remember_model(vwin, 1, NULL);
    }
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

    else if (mods & GDK_CONTROL_MASK) {
	if (gdk_keyval_to_upper(key->keyval) == GDK_S) 
	    file_viewer_save(NULL, vwin);
	else if (gdk_keyval_to_upper(key->keyval) == GDK_Q) {
	    if (vwin->role == EDIT_SCRIPT && SCRIPT_IS_CHANGED(vwin)) {
		gint resp;

		resp = query_save_script(NULL, NULL, vwin);
		if (!resp) gtk_widget_destroy(vwin->dialog);
	    } else 
		gtk_widget_destroy(w);
	}
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

void register_data (const char *fname, int record)
{    
    char datacmd[MAXLEN];

    /* basic accounting */
    data_status |= HAVE_DATA;
    orig_vars = datainfo->v;

    /* set appropriate data_status bits */
    if (fname == NULL) {
	data_status |= (GUI_DATA|MODIFIED_DATA);
    } else if (!(data_status & IMPORT_DATA)) {
	if (strstr(paths.datfile, paths.datadir) != NULL) 
	    data_status |= BOOK_DATA;
	else
	    data_status |= USER_DATA; 
    }

    /* sync main window with datafile */
    populate_varlist();
    set_sample_label(datainfo);
    main_menubar_state(TRUE);
    session_menu_state(TRUE);

    /* record opening of data file in command log */
    if (record && fname != NULL) {
	mkfilelist(1, fname);
	sprintf(datacmd, "open %s", fname);
	check_cmd(datacmd);
	cmd_init(datacmd); 
    } 
}

#define APPENDING(action) (action == APPEND_CSV || \
                           action == APPEND_GNUMERIC || \
                           action == APPEND_EXCEL)

/* ........................................................... */

static void get_worksheet_data (const char *fname, int datatype,
				int append)
{
    int err;
    void *handle;
    int (*sheet_get_data)(const char*, double ***, DATAINFO *, char *);

    if (datatype == GRETL_GNUMERIC) {
	if (gui_open_plugin("gnumeric_import-2", &handle)) return;
	sheet_get_data = get_plugin_function("wbook_get_data", handle);
    }
    else if (datatype == GRETL_EXCEL) {
	if (gui_open_plugin("excel_import-2", &handle)) return;
	sheet_get_data = get_plugin_function("excel_get_data", handle);
    }
    else if (datatype == GRETL_DES_DATA) {
        if (gui_open_plugin("des_import", &handle)) return;
        sheet_get_data = get_plugin_function("des_get_data", handle);
    }
    else {
	errbox(_("Unrecognized data type"));
	return;
    }

    if (sheet_get_data == NULL) {
        errbox(_("Couldn't load plugin function"));
        close_plugin(handle);
        return;
    }

    err = (*sheet_get_data)(fname, &Z, datainfo, errtext);
    close_plugin(handle);

    if (err == -1) /* the user canceled the import */
	return;

    if (err) {
	if (strlen(errtext)) errbox(errtext);
	else errbox(_("Failed to import spreadsheet data"));
	return;
    }

    if (append) {
	infobox(_("Data appended OK"));
	data_status |= MODIFIED_DATA;
	register_data(fname, 0);
    } else {
	data_status |= IMPORT_DATA;
	strcpy(paths.datfile, fname);
	register_data(fname, 1);
    }
}

/* ........................................................... */

void do_open_data (GtkWidget *w, gpointer data, int code)
     /* cases: 
	- called from dialog: user has said Yes to opening data file,
	although a data file is already open
	- reached without dialog, in expert mode or when no datafile
	is open yet
     */
{
    gint datatype, err;
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

    if (code == OPEN_CSV || code == APPEND_CSV) {
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
	prn.buf = NULL; prn.fp = stderr;
	if (datatype == GRETL_XML_DATA)
	    err = get_xmldata(&Z, datainfo, trydatfile, &paths, 
			      data_status, &prn, 1);
	else
	    err = get_data(&Z, datainfo, trydatfile, &paths, data_status, &prn);
    }

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(1, trydatfile);
	return;
    }	

    /* trash the practice files window that launched the query? */
    if (fwin != NULL) {
	gtk_widget_destroy(fwin->w); 
    }

    strcpy(paths.datfile, trydatfile);

    register_data(paths.datfile, 1);
}

/* ........................................................... */

void verify_open_data (gpointer userdata, int code)
     /* give user choice of not opening selected datafile,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert && 
	yes_no_dialog (_("gretl: open data"), 
		       _("Opening a new data file will automatically\n"
		       "close the current one.  Any unsaved work\n"
			 "will be lost.  Proceed to open data file?"), 0)) {
	return;
    } else {
	do_open_data(NULL, userdata, code);
    }
}

/* ........................................................... */

void verify_open_session (gpointer userdata)
     /* give user choice of not opening session file,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert &&
	yes_no_dialog (_("gretl: open session"), 
		       _("Opening a new session file will automatically\n"
		       "close the current session.  Any unsaved work\n"
		       "will be lost.  Proceed to open session file?"), 0))
	return;
    else 
	do_open_session(NULL, userdata);
}

/* ........................................................... */

void save_session (char *fname) 
{
    int i, spos;
    char msg[MAXLEN], savedir[MAXLEN], fname2[MAXLEN];
    char session_base[MAXLEN], tmp[MAXLEN], grftmp[64];
    FILE *fp;
    PRN *prn;

    spos = slashpos(fname);
    if (spos) 
	safecpy(savedir, fname, spos);
    else *savedir = 0;

#ifdef CMD_DEBUG
    dump_cmd_stack("stderr");
#endif

    /* save commands, by dumping the command stack */
    if (haschar('.', fname) < 0)
	strcat(fname, ".gretl");
    if (dump_cmd_stack(fname)) return;

    get_base(session_base, fname, '.');

    /* get ready to save "session" */
    fp = fopen(fname, "a");
    if (fp == NULL) {
	sprintf(errtext, _("Couldn't open session file %s"), fname);
	errbox(errtext);
	return;
    }
    fprintf(fp, "(* saved objects:\n");

    /* save session models */
    for (i=0; i<session.nmodels; i++) {
	fprintf(fp, "model %d \"%s\"\n", 
		(session.models[i])->ID, 
		(session.models[i])->name);
    }

    /* save session graphs */
    for (i=0; i<session.ngraphs; i++) {
	/* formulate save name for graph */
	strcpy(grftmp, (session.graphs[i])->name);
	sprintf(tmp, "%s%s", session_base, space_to_score(grftmp));
	/* does the constructed filename differ from the
	   current one? */
	if (strcmp((session.graphs[i])->fname, tmp)) {
	    if (copyfile((session.graphs[i])->fname, tmp)) {
		errbox(_("Couldn't copy graph file"));
		continue;
	    } else {
		remove((session.graphs[i])->fname);
		strcpy((session.graphs[i])->fname, tmp);
	    }
	}
	fprintf(fp, "%s %d \"%s\" %s\n", 
		((session.graphs[i])->name[0] == 'G')? "graph" : "plot",
		(session.graphs[i])->ID, 
		(session.graphs[i])->name, 
		(session.graphs[i])->fname);
    }

    fprintf(fp, "*)\n");
    fclose(fp);

    /* delete any extraneous graph files */
    session_file_manager(REALLY_DELETE_ALL, NULL);

    /* save session notes, if any */
    if (session.notes != NULL && strlen(session.notes)) {
	switch_ext(fname2, fname, "Notes");
	fp = fopen(fname2, "w");
	if (fp != NULL) {
	    fprintf(fp, "%s", session.notes);
	    fclose(fp);
	} else
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
    execute_script(fname, NULL, NULL, NULL, prn, SAVE_SESSION_EXEC); 
    gretl_print_destroy(prn);

    sprintf(msg, _("session saved to %s -\n"), savedir);
    strcat(msg, _("commands: "));
    strcat(msg, (spos)? fname + spos + 1 : fname);
    strcat(msg, _("\noutput: "));
    spos = slashpos(fname2);
    strcat(msg, (spos)? fname2 + spos + 1 : fname2);
    infobox(msg);

    mkfilelist(2, fname);
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
	data_status |= MODIFIED_DATA;
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
	    if (vwin->role == EDIT_SCRIPT) 
		SCRIPT_SAVED(vwin);
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
}

/* .................................................................. */

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin) {
	if (vwin->w) {
	    gchar *undo = g_object_get_data(G_OBJECT(vwin->w), "undo");
	    
	    if (undo) g_free(undo);
	}
	if (vwin->listbox) 
	    gtk_widget_destroy(GTK_WIDGET(vwin->listbox));
	if (vwin->mbar) 
	    gtk_widget_destroy(GTK_WIDGET(vwin->mbar));
	if (vwin->status) 
	    gtk_widget_destroy(GTK_WIDGET(vwin->status));
	if (vwin->popup) 
	    gtk_widget_destroy(GTK_WIDGET(vwin->popup));
	if (vwin->ifac) 
	    g_object_unref(G_OBJECT(vwin->ifac));  
	if (vwin->role == SUMMARY || vwin->role == VAR_SUMMARY)
	    free_summary(vwin->data); 
	else if (vwin->role == CORR)
	    free_corrmat(vwin->data);
	else if (vwin->role == MPOLS)
	    free_gretl_mp_results(vwin->data);
	if (vwin->dialog)
	    winstack_remove(vwin->dialog);
	free(vwin);
	vwin = NULL;
    }
}

#if defined(G_OS_WIN32) || defined(USE_GNOME) 
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(vwin, 0, w);
}
#endif

/* ........................................................... */

static void make_viewbar (windata_t *vwin)
{
    GtkWidget *button, *viewbar, *hbox;
    int i;
    static char *viewstrings[] = {
	N_("Save"),
	N_("Save as..."),
	N_("Print..."),
	N_("Run"),
	N_("Copy selection"), 
	N_("Paste"),
	N_("Find..."),
	N_("Replace..."),
	N_("Undo"),
	N_("Help on command"),
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
		   vwin->role == SCRIPT_OUT);

    int save_as_ok = (vwin->role != EDIT_HEADER && 
		      vwin->role != EDIT_NOTES);

#if defined(G_OS_WIN32) || defined(USE_GNOME)
    int print_ok = 1;
#else
    int print_ok = 0;
#endif

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    viewbar = gtk_toolbar_new();
    gtk_box_pack_start(GTK_BOX(hbox), viewbar, FALSE, FALSE, 0);

    for (i=0; viewstrings[i] != NULL; i++) {
	switch (i) {
	case 0:
	    if (edit_ok && vwin->role != SCRIPT_OUT) {
		stockicon = GTK_STOCK_SAVE;
		if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) 
		    toolfunc = buf_edit_save;
		else
		    toolfunc = file_viewer_save;
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
	    if (print_ok) {
#if defined(G_OS_WIN32) || defined(USE_GNOME)
		stockicon = GTK_STOCK_PRINT;
		toolfunc = window_print_callback;
#endif
	    } else
		toolfunc = NULL;
	    break;
	case 3:
	    if (run_ok) {
		stockicon = GTK_STOCK_EXECUTE;
		toolfunc = run_script_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 4:
	    stockicon = GTK_STOCK_COPY;
	    toolfunc = text_copy_callback;
	    break;
	case 5:
	    if (edit_ok) {
		stockicon = GTK_STOCK_PASTE;
		toolfunc = text_paste_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 6:
	    stockicon = GTK_STOCK_FIND;
	    toolfunc = text_find_callback;
	    break;
	case 7:
	    if (edit_ok) {
		stockicon = GTK_STOCK_FIND_AND_REPLACE;
		toolfunc = text_replace_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 8:
	    if (edit_ok) {
		stockicon = GTK_STOCK_UNDO;
		toolfunc = text_undo_callback;
	    } else
		toolfunc = NULL;
	    break;
	case 9:
	    if (run_ok) {
		stockicon = GTK_STOCK_HELP;
		toolfunc = activate_script_help;
	    } else
		toolfunc = NULL;
	    break;
	case 10:
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
				 GTK_ICON_SIZE_LARGE_TOOLBAR);
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
    case HELP: 
	title = g_strdup(_("gretl: help")); break;
    case CLI_HELP:
	title = g_strdup(_("gretl: command syntax")); break;
    case VIEW_LOG:
	title = g_strdup(_("gretl: command log")); break;
    case CONSOLE:
	title = g_strdup(_("gretl console")); break;
    case EDIT_SCRIPT:
    case VIEW_SCRIPT:	
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
	title = g_strdup(_("gretl: edit plot commands")); break;
    case SCRIPT_OUT:
	title = g_strdup(_("gretl: script output")); break;
    default:
	break;
    }
    return title;
}

/* ........................................................... */

static gint script_changed (GtkWidget *w, windata_t *vwin)
{
    SCRIPT_CHANGED(vwin);
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
    SCRIPT_SAVED(vwin);
}

/* ........................................................... */

static gint query_save_script (GtkWidget *w, GdkEvent *event, windata_t *vwin)
{
    if (SCRIPT_IS_CHANGED(vwin)) {
	int button;

	button = yes_no_dialog(_("gretl: script"), 
			       _("Save changes?"), 1);

	if (button == CANCEL_BUTTON)
	    return TRUE;
	if (button == YES_BUTTON)
	    auto_save_script(vwin);
    }
    return FALSE;
}

/* ........................................................... */

enum {
    PLAIN_TEXT,
    BLUE_TEXT,
    RED_TEXT
};

/* ........................................................... */

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

    if (record) winstack_add(vwin->dialog);

    return vwin;
}

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
    gtk_widget_set_size_request (vwin->dialog, hsize, vsize); 
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), editable);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), editable);
}

/* ........................................................... */

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_container_set_border_width (GTK_CONTAINER(vwin->vbox), 5);
    gtk_box_set_spacing (GTK_BOX(vwin->vbox), 5);
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
    }

    create_text (vwin, &tbuf, hsize, vsize, FALSE);
    
    dialog_table_setup (vwin);

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
		     G_CALLBACK(catch_key), vwin->dialog);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);
    cursor_to_top(vwin);
    return vwin;
}

/* ........................................................... */

windata_t *view_file (char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role, 
		      GtkItemFactoryEntry menu_items[]) 
{
    GtkTextBuffer *tbuf;
    GtkTextIter iter;
    int thiscolor, nextcolor;
    char tempstr[MAXSTR], *fle = NULL;
    FILE *fd = NULL;
    windata_t *vwin;
    gchar *title;
    int show_viewbar = (role != CONSOLE &&
			role != HELP &&
			role != CLI_HELP);
    int doing_script = (role == EDIT_SCRIPT ||
			role == VIEW_SCRIPT ||
			role == VIEW_LOG);

    fd = fopen(filename, "r");
    if (fd == NULL) {
	sprintf(errtext, _("Can't open %s for reading"), filename);
	errbox(errtext);
	return NULL;
    }

    title = make_viewer_title(role, filename);
    vwin = common_viewer_new(role, title, NULL, !doing_script);
    g_free(title);
    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);

    if (menu_items != NULL) {
	set_up_viewer_menu(vwin->dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    }

    create_text(vwin, &tbuf, hsize, vsize, editable);

    dialog_table_setup (vwin);

    /* special case: the gretl console */
    if (role == CONSOLE) {
	g_signal_connect(G_OBJECT(vwin->w), "key_press_event",
			 G_CALLBACK(console_handler), NULL);
    } 

    if (doing_script) {
	g_signal_connect(G_OBJECT(vwin->w), "button_release_event",
			 G_CALLBACK(edit_script_help), vwin);
    } 

    /* is the file to be deleted after viewing? */
    if (del_file) {
	if ((fle = mymalloc(strlen(filename) + 1)) == NULL) {
	    return NULL;
	}
	strcpy(fle, filename);
    }

    /* should we show a toolbar? */
    if (show_viewbar) { 
	make_viewbar(vwin);
    } else { /* make a simple Close button instead */
	GtkWidget *close = gtk_button_new_with_label(_("Close"));

	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   close, FALSE, TRUE, 0);
	g_signal_connect(G_OBJECT(close), "clicked", 
			 G_CALLBACK(delete_file_viewer), vwin);
	gtk_widget_show(close);
    }

    /* insert the file text */
    thiscolor = nextcolor = PLAIN_TEXT;
    gtk_text_buffer_get_iter_at_offset(tbuf, &iter, 0);
    memset(tempstr, 0, sizeof tempstr);
    while (fgets(tempstr, sizeof tempstr - 1, fd)) {
	if (tempstr[0] == '@') continue;
	if (tempstr[0] == '?') 
	    thiscolor = (role == CONSOLE)? RED_TEXT : BLUE_TEXT;
	if (tempstr[0] == '#') {
	    if (role == HELP || role == CLI_HELP) {
		tempstr[0] = ' ';
		nextcolor = RED_TEXT;
	    } else
		thiscolor = BLUE_TEXT;
	} else
	    nextcolor = PLAIN_TEXT;
	switch (thiscolor) {
	case PLAIN_TEXT:
	    gtk_text_buffer_insert(tbuf, &iter, tempstr, -1);
	    break;
	case BLUE_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      tempstr, -1,
						      "bluetext", NULL);
	    break;
	case RED_TEXT:
	    gtk_text_buffer_insert_with_tags_by_name (tbuf, &iter,
						      tempstr, -1,
						      "redtext", NULL);
	    break;
	}
	thiscolor = nextcolor;
	memset(tempstr, 0, sizeof tempstr);
    }
    fclose(fd);

    /* grab the "changed" signal when editing a script */
    if (role == EDIT_SCRIPT) {
	g_signal_connect(G_OBJECT(tbuf), "changed", 
			 G_CALLBACK(script_changed), vwin);
    }

    /* catch some keystrokes */
    if (!editable) {
	g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
			 G_CALLBACK(catch_key), vwin->dialog);
    } else {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
	g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
			 G_CALLBACK(catch_edit_key), vwin);	
    }  

    /* offer chance to save script on exit */
    if (role == EDIT_SCRIPT)
	g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
			 G_CALLBACK(query_save_script), vwin);

    /* clean up when dialog is destroyed */
    if (del_file) {
	g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
			 G_CALLBACK(delete_file), (gpointer) fle);
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
    GtkTextBuffer *tbuf;
    windata_t *vwin;

    vwin = common_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) return NULL;

    viewer_box_config(vwin); 

    /* add a menu bar */
    set_up_viewer_menu(vwin->dialog, vwin, edit_items);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    create_text(vwin, &tbuf, hsize, vsize, TRUE);

    dialog_table_setup (vwin);
    
    /* add an editing bar */
    make_viewbar(vwin);    

    /* insert the buffer text */
    if (*pbuf) gtk_text_buffer_set_text(tbuf, *pbuf, -1);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

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

    copylist(&default_list, pmod->list);

    /* attach shortcuts */
    g_object_set_data(G_OBJECT(vwin->dialog), "ddata", vwin);
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_key), vwin->dialog);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(delete_model), 
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

	if (w != NULL) 
	    gtk_widget_set_sensitive(w, s);
	else
	    fprintf(stderr, _("Failed to flip state of \"%s\"\n"), path);
    }
}

/* ........................................................... */

static void model_rtf_copy_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Edit/Copy all/as RTF", s);
}

/* ........................................................... */

static void model_latex_copy_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Edit/Copy all/as LaTeX", s);
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

static void lmmenu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/non-linearity (squares)", s);
    flip(ifac, "/Tests/non-linearity (logs)", s);
    flip(ifac, "/Tests/autocorrelation", s);
    flip(ifac, "/Tests/heteroskedasticity", s);
    flip(ifac, "/Tests/Chow test", s);
    flip(ifac, "/Tests/CUSUM test", s);
    flip(ifac, "/Tests/ARCH", s);
}

/* ........................................................... */

static void lad_menu (GtkItemFactory *ifac)
{
    /* FIXME: this is actually unnecessarily restrictive */
    flip(ifac, "/Tests", FALSE);
    flip(ifac, "/Model data/Add to data set", FALSE);
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

    if (vwin->role == SUMMARY || vwin->role == VAR_SUMMARY
	|| vwin->role == CORR || vwin->role == MPOLS) {
	augment_copy_menu(vwin);
	return;
    }

    if (vwin->role == VIEW_MODEL && vwin->data != NULL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	model_rtf_copy_state(vwin->ifac, !pmod->errcode);
	model_latex_copy_state(vwin->ifac, !pmod->errcode);
	latex_menu_state(vwin->ifac, !pmod->errcode);

	model_panel_menu_state(vwin->ifac, pmod->ci == POOLED);

	lmmenu_state(vwin->ifac, pmod->ci == OLS || pmod->ci == POOLED);

	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    model_menu_state(vwin->ifac, FALSE);
	    model_ml_menu_state(vwin->ifac, TRUE);
	} else {
	    model_ml_menu_state(vwin->ifac, FALSE);
	}
	if (pmod->name) {
	    model_save_state(vwin->ifac, FALSE);
	}

	if (pmod->ci == LAD) lad_menu(vwin->ifac);

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

    varitem.path = NULL;

    for (i=0; i<2; i++) {
	varitem.path = mymalloc(64);
	varitem.accelerator = NULL;
	varitem.callback_action = 0; 
	varitem.item_type = NULL;
	if (dataset_is_time_series(datainfo))
	    sprintf(varitem.path, _("%s/against time"), mpath[i]);
	else
	    sprintf(varitem.path, _("%s/by observation number"), mpath[i]);
	if (i == 0)
	    varitem.callback = resid_plot; 
	else
	    varitem.callback = fit_actual_plot;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);

	varstart = (i == 0)? 1 : 2;

	/* put the indep vars on the menu list */
	for (j=varstart; j<pmod->list[0]; j++) {
	    if (pmod->list[j] == 0) continue;
	    if (varitem.path == NULL)
		varitem.path = mymalloc(64);
	    varitem.accelerator = NULL;
	    varitem.callback_action = pmod->list[j]; 
	    varitem.item_type = NULL;
	    sprintf(varitem.path, _("%s/against %s"), mpath[i], 
		    datainfo->varname[pmod->list[j]]);
	    if (i == 0)
		varitem.callback = resid_plot; 
	    else
		varitem.callback = fit_actual_plot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	}
    }
    free(varitem.path);
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

    dumitem.path = NULL;

    /* put the dummy independent vars on the menu list */
    for (i=2; i<pmod->list[0]; i++) {
	if (pmod->list[i] == 0) continue;
	if (!isdummy(pmod->list[i], datainfo->t1, datainfo->t2, Z))
	    continue;
	if (!dums) { /* add separator, branch and "none" */
	    dumitem.path = mymalloc(64);
	    sprintf(dumitem.path, _("/Graphs/dumsep"));
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Separator>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* menu branch */
	    sprintf(dumitem.path, _("/Graphs/Separation"));
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Branch>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* "none" option */
	    sprintf(dumitem.path, _("/Graphs/Separation/none"));
	    dumitem.callback = plot_dummy_call;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<RadioItem>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    dums = 1;
	} 
	dumitem.callback_action = pmod->list[i]; 
	sprintf(dumitem.path, _("/Graphs/Separation/by %s"),  
		datainfo->varname[pmod->list[i]]);
	dumitem.callback = plot_dummy_call;	    
	dumitem.accelerator = NULL;
	dumitem.item_type = _("/Graphs/Separation/none");
	gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
    }
    free(dumitem.path);
}

/* ........................................................... */

static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data)
{
    windata_t *mwin = (windata_t *) data;
    MODEL *pmod = mwin->data;
    extern int quiet_sample_check (MODEL *pmod);
    int s, ok = 1;

    if (Z == NULL) {
	flip(mwin->ifac, "/File/Save to sesssion as icon", FALSE);
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
	(gtk_item_factory_get_item(mwin->ifac, "/Tests/omit variables"));
    if ((s && ok) || (!s && !ok)) return FALSE;
    s = !s;

    flip(mwin->ifac, "/Tests/omit variables", s);
    flip(mwin->ifac, "/Tests/add variables", s);
    flip(mwin->ifac, "/Tests/non-linearity (squares)", s);
    flip(mwin->ifac, "/Tests/non-linearity (logs)", s);
    flip(mwin->ifac, "/Tests/autocorrelation", s);
    flip(mwin->ifac, "/Tests/heteroskedasticity", s);
    flip(mwin->ifac, "/Tests/Chow test", s);
    flip(mwin->ifac, "/Tests/CUSUM test", s);
    flip(mwin->ifac, "/Tests/ARCH", s);
    flip(mwin->ifac, "/Graphs", s);
    flip(mwin->ifac, "/Model data/Display actual, fitted, residual", s);
    flip(mwin->ifac, "/Model data/Forecasts with standard errors", s);
    flip(mwin->ifac, "/Model data/Confidence intervals for coefficients", s);
    flip(mwin->ifac, "/Model data/Add to data set/fitted values", s);
    flip(mwin->ifac, "/Model data/Add to data set/residuals", s);
    flip(mwin->ifac, "/Model data/Add to data set/squared residuals", s);
    flip(mwin->ifac, "/Model data/Define new variable...", s);

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
    if (!(isalpha(varname[0]))) {
	sprintf(errtext, _("First char of name ('%c') is bad\n"
	       "(first must be alphabetical)"), varname[0]);
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

void prn_to_clipboard (PRN *prn)
{
    size_t len;

    if (prn->buf == NULL) return;
    len = strlen(prn->buf);

    if (clipboard_buf) g_free(clipboard_buf);
    clipboard_buf = mymalloc(len + 1);

    memcpy(clipboard_buf, prn->buf, len + 1);
    gtk_selection_owner_set(mdata->w,
			    GDK_SELECTION_PRIMARY,
			    GDK_CURRENT_TIME);
}

/* .................................................................. */

#define SPECIAL_COPY(h) (h == COPY_LATEX || h == COPY_RTF)

void text_copy (gpointer data, guint how, GtkWidget *widget) 
{
    windata_t *vwin = (windata_t *) data;
    PRN *prn;

    /* descriptive statistics */
    if ((vwin->role == SUMMARY || vwin->role == VAR_SUMMARY)
	&& SPECIAL_COPY(how)) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) vwin->data;
	
	if (bufopen(&prn)) return;
	if (how == COPY_LATEX) {
	    texprint_summary(summ, datainfo, prn);
	    prn_to_clipboard(prn);
	} else { /* RTF */
	    rtfprint_summary(summ, datainfo, prn);
#ifdef G_OS_WIN32
	    win_copy_text(prn, COPY_RTF);
#else
	    prn_to_clipboard(prn);
#endif
	}
	gretl_print_destroy(prn);
	return;
    }

    /* correlation matrix */
    if (vwin->role == CORR && SPECIAL_COPY(how)) {
	CORRMAT *corr = (CORRMAT *) vwin->data;

	if (bufopen(&prn)) return;
	if (how == COPY_LATEX) { 
	    texprint_corrmat(corr, datainfo, prn);
	    prn_to_clipboard(prn);
	} 
	else { /* RTF */
	    rtfprint_corrmat(corr, datainfo, prn);
#ifdef G_OS_WIN32
	    win_copy_text(prn, COPY_RTF);
#else
	    prn_to_clipboard(prn);
#endif
	}
	gretl_print_destroy(prn);
	return;
    }

    /* multiple-precision OLS */
    if (vwin->role == MPOLS && SPECIAL_COPY(how)) {
	mp_results *mpvals = (mp_results *) vwin->data;

	if (bufopen(&prn)) return;
	if (how == COPY_LATEX) { 
	    prn->format = GRETL_PRINT_FORMAT_TEX;
	    print_mpols_results (mpvals, datainfo, prn);
	    prn_to_clipboard(prn);
	} 
	else { /* RTF */
	    prn->format = GRETL_PRINT_FORMAT_RTF;
	    print_mpols_results (mpvals, datainfo, prn);
#ifdef G_OS_WIN32
	    win_copy_text(prn, COPY_RTF);
#else
	    prn_to_clipboard(prn);
#endif
	}
	gretl_print_destroy(prn);
	return;
    }

    /* or it's a model window we're copying from? */
    if (vwin->role == VIEW_MODEL &&
	(how == COPY_RTF || how == COPY_LATEX ||
	 how == COPY_LATEX_EQUATION)) {
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->errcode) { 
	    errbox("Couldn't format model");
	    return;
	}

	if (how == COPY_RTF) {
	    model_to_rtf(pmod);
	    return;
	}

	if (bufopen(&prn)) return;

	if (how == COPY_LATEX) {
	    tex_print_model(pmod, datainfo, 0, prn);
	}
	else if (how == COPY_LATEX_EQUATION) {
	    tex_print_equation(pmod, datainfo, 0, prn);
	}

	prn_to_clipboard(prn);
	gretl_print_destroy(prn);
	return;
    }

    /* otherwise copying plain text from window */
    if (how == COPY_TEXT) {
	PRN textprn;

	textprn.fp = NULL;
	textprn.buf = textview_get_text(GTK_TEXT_VIEW(vwin->w));
# ifdef G_OS_WIN32
	win_copy_text(&textprn, COPY_TEXT);	
# else
	prn_to_clipboard(&textprn);
# endif /* G_OS_WIN32 */
	g_free(textprn.buf);
    } else { /* COPY_SELECTION */
	gtk_text_buffer_copy_clipboard (gtk_text_view_get_buffer
					(GTK_TEXT_VIEW(vwin->w)),
					 gtk_clipboard_get(GDK_NONE));
    }
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

    if (gtk_text_buffer_get_selection_bounds(tbuf, &start, &end))
	selbuf = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

    winprint(buf, selbuf);
}

#endif

/* .................................................................. */

void text_undo (windata_t *vwin, guint u, GtkWidget *widget)
{
    gchar *old = g_object_steal_data(G_OBJECT(vwin->w), "undo");

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

/* .................................................................. */

/* #ifndef G_OS_WIN32 */
#define USE_GTK_STOCK_BUTTONS
/* #endif */

#ifndef USE_GTK_STOCK_BUTTONS
struct button_type {
    const gchar *flag;
    const gchar *text;
};
#endif

GtkWidget *standard_button (const gchar *flag)
{
#ifdef USE_GTK_STOCK_BUTTONS
    return gtk_button_new_from_stock(flag);
#else
    struct button_type mybuttons[] = {
	{ GTK_STOCK_APPLY, N_("Apply") },
	{ GTK_STOCK_CANCEL, N_("Cancel") },
	{ GTK_STOCK_CLEAR, N_("Clear") },
	{ GTK_STOCK_CLOSE, N_("Close") },
	{ GTK_STOCK_HELP, N_("Help") }, 
	{ GTK_STOCK_OK, N_("OK") },
	{ GTK_STOCK_SAVE, N_("Save") },
	{ NULL, NULL }
    };

    int i;

    for (i=0; mybuttons[i].flag != NULL; i++) {
	if (!strcmp(flag, mybuttons[i].flag)) {
	    return gtk_button_new_with_label(_(mybuttons[i].text));
	}
    }

    return NULL;
#endif
}
