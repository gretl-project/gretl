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
#include "htmlprint.h"
#include "guiprint.h"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#if !defined(G_OS_WIN32) && !defined(USE_GNOME)
char rcfile[MAXLEN];
#endif

char *storelist = NULL;

static GtkWidget *help_view = NULL;

extern GtkTooltips *gretl_tips;
extern int session_saved;
extern GtkWidget *mysheet;
extern GtkWidget *toolbar_box;
extern char *space_to_score (char *str);

extern int want_toolbar;
extern char calculator[MAXSTR];
extern char editor[MAXSTR];
extern char Rcommand[MAXSTR];
extern char dbproxy[21];
int use_proxy;

/* filelist stuff */
#define MAXRECENT 4

#ifdef USE_GNOME
static void gnome_printfilelist (int filetype);
#else
# ifdef G_OS_WIN32
static void win_printfilelist (int filetype);
# endif
static void printfilelist (int filetype, FILE *fp);
#endif

static char datalist[MAXRECENT][MAXSTR], *datap[MAXRECENT];
static char sessionlist[MAXRECENT][MAXSTR], *sessionp[MAXRECENT];
static char scriptlist[MAXRECENT][MAXSTR], *scriptp[MAXRECENT];

/* helpfile vars */
static int help_length, gui_help_length, script_help_length;

/* searching stuff */
static int look_for_string (char *haystack, char *needle, int nStart);
static void close_find_dialog (GtkWidget *widget, gpointer data);
static void find_in_help (GtkWidget *widget, gpointer data);
static void find_in_clist (GtkWidget *widget, gpointer data);
static void cancel_find (GtkWidget *widget, gpointer data);
static void find_string_dialog (void (*YesFunc)(), void (*NoFunc)(),
				gpointer data);
static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static GtkWidget *find_window = NULL;
static GtkWidget *find_entry;
static char *needle;

static void make_prefs_tab (GtkWidget *notebook, int tab);
static void apply_changes (GtkWidget *widget, gpointer data);
#ifndef G_OS_WIN32
static void read_rc (void);
#endif

extern void do_coeff_intervals (gpointer data, guint i, GtkWidget *w);
extern void save_plot (char *fname, GPT_SPEC *plot);
extern gboolean console_handler (GtkWidget *w, GdkEventKey *key, 
				 gpointer user_data);
extern void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w);

/* font handling */
static char fontspec[MAXLEN] = 
"-b&h-lucidatypewriter-medium-r-normal-sans-12-*-*-*-*-*-*-*";
GdkFont *fixed_font;

static int usecwd;
int olddat;

typedef struct {
    char *key;         /* config file variable name */
    char *description; /* How the field will show up in the options dialog */
    char *link;        /* in case of radio button pair, alternate string */
    void *var;         /* pointer to variable */
    char type;         /* 'U' (user) or 'R' (root) for string, 'B' for boolean */
    int len;           /* storage size for string variable (also see Note) */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
} RCVARS;

/* Note: actually "len" above is overloaded: if an rc_var is of type 'B'
   (boolean) and not part of a radio group, then a non-zero value for
   len will link the var's toggle button with the sensitivity of the
   preceding rc_var's entry field.  For example, the "use_proxy" button
   controls the sensitivity of the "dbproxy" entry widget. */

RCVARS rc_vars[] = {
    {"gretldir", "Main gretl directory", NULL, paths.gretldir, 
     'R', MAXLEN, 1, NULL},
    {"userdir", "User's gretl directory", NULL, paths.userdir, 
     'U', MAXLEN, 1, NULL},
    {"gnuplot", "Command to launch gnuplot", NULL, paths.gnuplot, 
     'R', MAXLEN, 1, NULL},
    {"Rcommand", "Command to launch GNU R", NULL, Rcommand, 
     'R', MAXSTR, 1, NULL},
    {"expert", "Expert mode (no warnings)", NULL, &expert, 
     'B', 0, 1, NULL},
    {"updater", "Tell me about gretl updates", NULL, &updater, 
     'B', 0, 1, NULL},
    {"binbase", "gretl database directory", NULL, paths.binbase, 
     'U', MAXLEN, 2, NULL},
    {"ratsbase", "RATS data directory", NULL, paths.ratsbase, 
     'U', MAXLEN, 2, NULL},
    {"dbhost_ip", "Database server IP", NULL, paths.dbhost_ip, 
     'U', 16, 2, NULL},
    {"dbproxy", "HTTP proxy (ipnumber:port)", NULL, dbproxy, 
     'U', 21, 2, NULL},
    {"useproxy", "Use HTTP proxy", NULL, &use_proxy, 
     'B', 1, 2, NULL},
    {"calculator", "Calculator", NULL, calculator, 
     'U', MAXSTR, 3, NULL},
    {"editor", "Editor", NULL, editor, 
     'U', MAXSTR, 3, NULL},
    {"toolbar", "Show gretl toolbar", NULL, &want_toolbar, 
     'B', 0, 3, NULL},
    {"usecwd", "Use current working directory as default", 
     "Use gretl user directory as default", &usecwd, 'B', 0, 4, NULL},
    {"olddat", "Use \".dat\" as default datafile suffix", 
     "Use \".gdt\" as default suffix", &olddat, 'B', 0, 5, NULL},
    {"fontspec", "Fixed font", NULL, fontspec, 'U', MAXLEN, 0, NULL},
    {NULL, NULL, NULL, NULL, 0, 0, 0, NULL}   
};

GtkItemFactoryEntry model_items[] = {
    { "/_File", NULL, NULL, 0, "<Branch>" },
    { "/File/_Save as text...", NULL, file_save, SAVE_MODEL, NULL },
    { "/File/Save to session as icon", NULL, remember_model, 0, NULL },
    { "/File/Save as icon and close", NULL, remember_model, 1, NULL },
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { "/File/_Print...", NULL, window_print, 0, NULL },
#endif
    { "/_Edit", NULL, NULL, 0, "<Branch>" },
    { "/Edit/_Copy selection", NULL, text_copy, COPY_SELECTION, NULL },
    { "/Edit/Copy _all", NULL, NULL, 0, "<Branch>" },
    { "/Edit/Copy _all/as plain _text", NULL, text_copy, COPY_TEXT, NULL },
    { "/Edit/Copy _all/as _HTML", NULL, text_copy, COPY_HTML, NULL },
    { "/Edit/Copy _all/as _LaTeX", NULL, text_copy, COPY_LATEX, NULL },
    { "/Edit/Copy _all/as _RTF", NULL, text_copy, COPY_RTF, NULL },
    { "/_Tests", NULL, NULL, 0, "<Branch>" },    
    { "/Tests/omit variables", NULL, model_test_callback, OMIT, NULL },
    { "/Tests/add variables", NULL, model_test_callback, ADD, NULL },
    { "/Tests/sep1", NULL, NULL, 0, "<Separator>" },
    { "/Tests/non-linearity (squares)", NULL, do_lmtest, AUX_SQ, NULL },
    { "/Tests/non-linearity (logs)", NULL, do_lmtest, AUX_LOG, NULL },
    { "/Tests/sep2", NULL, NULL, 0, "<Separator>" },
    { "/Tests/autocorrelation", NULL, do_lmtest, AUX_AR, NULL },
    { "/Tests/heteroskedasticity", NULL, do_lmtest, AUX_WHITE, NULL },
    { "/Tests/Chow test", NULL, model_test_callback, CHOW, NULL },
    { "/Tests/CUSUM test", NULL, do_cusum, 0, NULL },
    { "/Tests/ARCH", NULL, model_test_callback, ARCH, NULL },
    { "/Tests/normality of residual", NULL, do_resid_freq, 0, NULL },
    { "/Tests/panel diagnostics", NULL, do_panel_diagnostics, 0, NULL },
    { "/_Graphs", NULL, NULL, 0, "<Branch>" }, 
    { "/Graphs/residual plot", NULL, NULL, 0, "<Branch>" },
    { "/Graphs/fitted, actual plot", NULL, NULL, 0, "<Branch>" },
    { "/_Model data", NULL, NULL, 0, "<Branch>" },
    { "/_Model data/Display actual, fitted, residual", NULL, 
      display_fit_resid, 0, NULL },
    { "/_Model data/Forecasts with standard errors", NULL, 
      model_test_callback, FCAST, NULL },
    { "/_Model data/Confidence intervals for coefficients", NULL, 
      do_coeff_intervals, 0, NULL },
    { "/_Model data/Add to data set/fitted values", NULL, 
      fit_resid_callback, 1, NULL },
    { "/_Model data/Add to data set/residuals", NULL, 
      fit_resid_callback, 0, NULL },
    { "/_Model data/Add to data set/squared residuals", NULL, 
      fit_resid_callback, 2, NULL },
    { "/_Model data/Add to data set/error sum of squares", NULL, 
      model_stat_callback, ESS, NULL },
    { "/_Model data/Add to data set/standard error of residuals", NULL, 
      model_stat_callback, SIGMA, NULL },
    { "/_Model data/Add to data set/R-squared", NULL, 
      model_stat_callback, R2, NULL },
    { "/_Model data/Add to data set/T*R-squared", NULL, 
      model_stat_callback, TR2, NULL },
    { "/_Model data/Add to data set/degrees of freedom", NULL, 
      model_stat_callback, DF, NULL },
    { "/_Model data/coefficient covariance matrix", NULL, 
      do_outcovmx, 0, NULL },
    { "/_Model data/sep1", NULL, NULL, 0, "<Separator>" },
    { "/_Model data/Define new variable...", NULL, model_test_callback, 
      MODEL_GENR, NULL },
    { "/_LaTeX", NULL, NULL, 0, "<Branch>" },
    { "/LaTeX/_View", NULL, NULL, 0, "<Branch>" },
    { "/LaTeX/View/_Tabular", NULL, view_latex, 0, NULL },
    { "/LaTeX/View/_Equation", NULL, view_latex, 1, NULL },
    { "/LaTeX/_Save", NULL, NULL, 0, "<Branch>" },
    { "/LaTeX/Save/_Tabular", NULL, file_save, SAVE_TEX_TAB, NULL },
    { "/LaTeX/Save/_Equation", NULL, file_save, SAVE_TEX_EQ, NULL },
    { "/LaTeX/_Copy", NULL, NULL, 0, "<Branch>" },
    { "/LaTeX/Copy/_Tabular", NULL, text_copy, COPY_LATEX, NULL },
    { "/LaTeX/Copy/_Equation", NULL, text_copy, COPY_LATEX_EQUATION, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

GtkItemFactoryEntry help_items[] = {
    { "/_Topics", NULL, NULL, 0, "<Branch>" },    
    { "/Topics/Generate variable syntax", NULL, do_help, GENR, NULL },
    { "/Topics/sep1", NULL, NULL, 0, "<Separator>" },
    { "/Topics/Graphing", NULL, do_help, GNUPLOT, NULL },
    { "/Topics/sep2", NULL, NULL, 0, "<Separator>" },
    { "/Topics/Estimation/Ordinary Least Squares", 
      NULL, do_help, OLS, NULL },
    { "/Topics/Estimation/Weighted Least Squares", 
      NULL, do_help, WLS, NULL },
    { "/Topics/Estimation/Cochrane-Orcutt", NULL, do_help, CORC, NULL },
    { "/Topics/Estimation/HCCM", NULL, do_help, HCCM, NULL },
    { "/Topics/Estimation/Hildreth-Lu", NULL, do_help, HILU, NULL },
    { "/Topics/Estimation/Heteroskedasticity", NULL, do_help, HSK, NULL },
    { "/Topics/Estimation/Autoregressive estimation", NULL, do_help, AR, NULL },
    { "/Topics/Estimation/Vector Autoregression", NULL, do_help, VAR, NULL },
    { "/Topics/Estimation/Two-Stage Least Squares", NULL, do_help, TSLS, NULL }, 
    { "/Topics/Estimation/_Logit", NULL, do_help, LOGIT, NULL }, 
    { "/Topics/Estimation/_Probit", NULL, do_help, PROBIT, NULL }, 
    { "/Topics/Estimation/_Rank Correlation", NULL, do_help, SPEARMAN, NULL },
    { "/Topics/Estimation/Pooled OLS (panel)", NULL, do_help, POOLED, NULL }, 
    { "/Topics/sep3", NULL, NULL, 0, "<Separator>" },
    { "/Topics/Hypothesis tests/omit variables", NULL, do_help, OMIT, NULL },
    { "/Topics/Hypothesis tests/add variables", NULL, do_help, ADD, NULL },
    { "/Topics/Hypothesis tests/LM test", NULL, do_help, LMTEST, NULL },
    { "/Topics/Hypothesis tests/Dickey-Fuller test", 
      NULL, do_help, ADF, NULL },
    { "/Topics/Hypothesis tests/Chow test", NULL, do_help, CHOW, NULL },
    { "/Topics/Hypothesis tests/Cointegration test", NULL, do_help, 
      COINT, NULL },
    { "/Topics/Hypothesis tests/_Panel diagnostics", NULL, do_help, POOLED, NULL },
    { "/_Find", NULL, menu_find, 0, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

GtkItemFactoryEntry script_help_items[] = {
    { "/_Find", NULL, menu_find, 0, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

GtkItemFactoryEntry edit_items[] = {
#if defined(G_OS_WIN32) || defined(USE_GNOME)
    { "/File/_Print...", NULL, window_print, 0, NULL },
#endif    
    { "/_Edit", NULL, NULL, 0, "<Branch>" },
    { "/Edit/_Copy selection", NULL, text_copy, COPY_SELECTION, NULL },
    { "/Edit/Copy _all", NULL, text_copy, COPY_TEXT, NULL },
    { "/Edit/_Paste", NULL, text_paste, 0, NULL },
    { "/Edit/_Replace...", NULL, text_replace, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

/* ........................................................... */

void load_fixed_font (void)
{
    /* get a monospaced font for various windows */
    fixed_font = gdk_font_load(fontspec);
}

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
    if (dir != NULL) strcat(fname, dir);
    strcat(fname, SLASHSTR);
}

/* ........................................................... */

#if !defined(G_OS_WIN32) && !defined(USE_GNOME)
void set_rcfile (void) 
{
    char *tmp;

    tmp = getenv("HOME");
    strcpy(rcfile, tmp);
    strcat(rcfile, "/.gretlrc");
    read_rc(); 
}
#endif

#ifdef USE_GNOME
void set_rcfile (void)
{
    read_rc();
}
#endif

/* ........................................................... */

static void delete_file (GtkWidget *widget, char *fle) 
{
    remove(fle);
    g_free(fle);
}

/* ........................................................... */

static void delete_file_viewer (GtkWidget *widget, gpointer data) 
{
    gtk_widget_destroy((GtkWidget *) data);
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
    gtk_widget_destroy(data);
}

/* ........................................................... */

void catch_key (GtkWidget *w, GdkEventKey *key)
{
    if (key->keyval == GDK_q) 
        gtk_widget_destroy(w);
    else if (key->keyval == GDK_s) 
	remember_model 
	    (gtk_object_get_data(GTK_OBJECT(w), "ddata"), 1, NULL);
}

/* ........................................................... */

void *mymalloc (size_t size) 
{
    void *mem;
   
    if((mem = malloc(size)) == NULL) 
	errbox("Out of memory!");
    return mem;
}

/* ........................................................... */

void *myrealloc (void *ptr, size_t size) 
{
    void *mem;
   
    if ((mem = realloc(ptr, size)) == NULL) 
	errbox("Out of memory!");
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
    if (fname == NULL)
	data_status |= (GUI_DATA|MODIFIED_DATA);
    else if (!(data_status & IMPORT_DATA)) {
	if (strstr(paths.datfile, paths.datadir) != NULL) 
	    data_status |= BOOK_DATA;
	else
	    data_status |= USER_DATA; 
    }

    /* sync main window with datafile */
    populate_clist(mdata->listbox, datainfo);
    set_sample_label(datainfo);
    menubar_state(TRUE);
    session_state(TRUE);

    /* record opening of data file in command log */
    if (record && fname != NULL) {
	mkfilelist(1, fname);
	sprintf(datacmd, "open %s", fname);
	check_cmd(datacmd);
	cmd_init(datacmd); 
    } 
}

/* ........................................................... */

void do_open_data (GtkWidget *w, gpointer data)
     /* cases: 
	- called from dialog: user has said Yes to opening data file,
	although a data file is already open
	- reached without dialog, in expert mode or when no datafile
	is open yet
     */
{
    gint datatype, err;
    PRN *prn;
    dialog_t *d = NULL;
    windata_t *fwin = NULL;

    if (data != NULL) {    
	if (w == NULL) { /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	} else {
	    d = (dialog_t *) data;
	    fwin = (windata_t *) d->data;
	}
    }

    /* check file type first */
    if (bufopen(&prn)) return;
    datatype = detect_filetype(trydatfile, &paths, prn);
    gretl_print_destroy(prn);

    /* will this work right? */
    close_session();

    if (datatype == GRETL_CSV_DATA) {
	do_open_csv_box(trydatfile, OPEN_CSV);
	return;
    }
    else if (datatype == GRETL_BOX_DATA) {
	do_open_csv_box(trydatfile, OPEN_BOX);
	return;
    }
    else { /* native data */
	PRN prn;
	prn.buf = NULL; prn.fp = stderr;
	if (datatype == GRETL_XML_DATA)
	    err = get_xmldata(&Z, datainfo, trydatfile, &paths, data_status, &prn);
	else
	    err = get_data(&Z, datainfo, trydatfile, &paths, data_status, &prn);
    }

    if (err) {
	gui_errmsg(err);
	return;
    }	

    /* trash the practice files window that launched the query? */
    if (fwin) gtk_widget_destroy(fwin->w); 

    strcpy(paths.datfile, trydatfile);

    register_data(paths.datfile, 1);
}

/* ........................................................... */

void verify_open_data (gpointer userdata)
     /* give user choice of not opening selected datafile,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert && 
	yes_no_dialog ("gretl: open data", 
		       "Opening a new data file will automatically\n"
		       "close the current one.  Any unsaved work\n"
		       "will be lost.  Proceed to open data file?", 0))
	return;
    else 
	do_open_data(NULL, userdata);
}

/* ........................................................... */

void verify_open_session (gpointer userdata)
     /* give user choice of not opening session file,
	if there's already a datafile open and we're not
	in "expert" mode */
{
    if (data_status && !expert &&
	yes_no_dialog ("gretl: open session", 
		       "Opening a new session file will automatically\n"
		       "close the current session.  Any unsaved work\n"
		       "will be lost.  Proceed to open session file?", 0))
	return;
    else 
	do_open_session(NULL, userdata);
}

/* ........................................................... */

static void set_data_from_filelist (gpointer data, guint i, 
				    GtkWidget *widget)
{
    strcpy(trydatfile, datap[i]); 
    verify_open_data(NULL);
}

/* ........................................................... */

static void set_session_from_filelist (gpointer data, guint i, 
				       GtkWidget *widget)
{
    strcpy(tryscript, sessionp[i]);
    verify_open_session(NULL);
}

/* ........................................................... */

static void set_script_from_filelist (gpointer data, guint i, 
				      GtkWidget *widget)
{
    strcpy(tryscript, scriptp[i]);
    do_open_script(NULL, NULL);
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
	sprintf(errtext, "Couldn't open session file %s", fname);
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
		errbox("Couldn't copy graph file");
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
    session_saved = 1;
    mkfilelist(2, fname);

    /* save session notes, if any */
    sprintf(tmp, "%ssession.Notes", savedir);
    fp = fopen(tmp, "r"); 
    if (fp != NULL) {
	char test[5];

	if (fgets(test, 4, fp)) { /* don't save empty notes */
	    fclose(fp);
	    switch_ext(fname2, fname, "Notes");
	    copyfile(tmp, fname2);
	} else 
	    fclose(fp);
	remove(tmp);
    }

    /* save output */
    switch_ext(fname2, fname, "txt");
    prn = gretl_print_new(GRETL_PRINT_FILE, fname2);
    if (prn == NULL) {
	errbox("Couldn't open output file for writing");
	return;
    }

    gui_logo(prn->fp);
    session_time(prn->fp);
    pprintf(prn, "Output from %s\n", fname);
    execute_script(fname, NULL, NULL, prn, SAVE_SESSION_EXEC); 
    gretl_print_destroy(prn);

    sprintf(msg, "session saved to %s -\n", savedir);
    strcat(msg, "commands: ");
    strcat(msg, (spos)? fname + spos + 1 : fname);
    strcat(msg, "\noutput: ");
    spos = slashpos(fname2);
    strcat(msg, (spos)? fname2 + spos + 1 : fname2);
    infobox(msg);

    return;
}

/* ......................................................... */

void helpfile_init (void)
{
    FILE *fp;
    char testline[MAXLEN];

    fp = fopen(paths.helpfile, "r");
    if (fp != NULL) { 
	while (fgets(testline, MAXLEN-1, fp)) gui_help_length++;
	fclose(fp);
    } else fprintf(stderr, "help file %s is not accessible\n", 
		   paths.helpfile);

    fp = fopen(paths.cmd_helpfile, "r");
    if (fp != NULL) { 
	while (fgets(testline, MAXLEN-1, fp)) script_help_length++;
	fclose(fp);
    } else fprintf(stderr, "help file %s is not accessible\n", 
		   paths.cmd_helpfile);

    help_length = gui_help_length;
}

/* ........................................................... */

static windata_t *helpwin (int script) 
{
    windata_t *vwin = NULL;

    if (script) {
	help_length = script_help_length;
	vwin = view_file(paths.cmd_helpfile, 0, 0, 77, 400, 
			 "gretl: command syntax help", script_help_items);
    } else {
	help_length = gui_help_length;
	vwin = view_file(paths.helpfile, 0, 0, 77, 400, "gretl: help", 
			 help_items);
    }
    return vwin;
}

/* ........................................................... */

static int help_index (const char *str)
{
    gchar helpline[84];
    gchar word[32];
    int found = 0, pos = 0;
    FILE *fp;

    fp = fopen(paths.helpfile, "r");
    while (!found) {
	if (fgets(helpline, 83, fp) == NULL) break;
	if (helpline[0] == '#') {
	    fgets(helpline, 83, fp);
	    sscanf(helpline, "%s", word);
	    if (strcmp(word, str) == 0) {
		found = 1;
	    } else pos++;
	}
	pos++;
    }
    fclose(fp);

    if (!found) return -1;
    return pos - 1;
}

/* ........................................................... */

void menu_find (gpointer data, guint db, GtkWidget *widget)
{
    if (db) 
	find_string_dialog(find_in_clist, cancel_find, data);
    else 
	find_string_dialog(find_in_help, cancel_find, data);
}

/* ........................................................... */

void datafile_find (GtkWidget *widget, gpointer data)
{
    find_string_dialog(find_in_clist, cancel_find, data);
}

/* ........................................................... */

void context_help (GtkWidget *widget, gpointer data)
{
    int help_code = GPOINTER_TO_INT(data);

    do_help(NULL, (guint) help_code, NULL);
}

/* ........................................................... */

void help_show (gpointer data, guint cli, GtkWidget *widget)
{
    if (help_view == NULL) {
	windata_t *vwin = helpwin(cli);

	if (vwin != NULL) help_view = vwin->w;
	gtk_signal_connect(GTK_OBJECT(help_view), "destroy",
			   GTK_SIGNAL_FUNC(gtk_widget_destroyed),
			   &help_view);	
    } else {
	gdk_window_show(help_view->parent->window);
	gdk_window_raise(help_view->parent->window);
	gtk_adjustment_set_value(GTK_TEXT(help_view)->vadj, 0.0);
    }
} 

/* ........................................................... */

void do_help (gpointer data, guint code, GtkWidget *widget) 
{
    int pos = 0;
    char cmdstr[12];

    if (code == MEANTEST2) code = MEANTEST;
    else if (code == MODEL_GENR) code = GENR;

    if (code == GR_PLOT || code == GR_XY || code == GNUPLOT)
	strcpy(cmdstr, "graphing");
    else if (code == GR_DUMMY)
	strcpy(cmdstr, "factorized");
    else if (code == GR_BOX || code == GR_NBOX)
	strcpy(cmdstr, "boxplots");
    else if (code == ONLINE)
	strcpy(cmdstr, "online");
    else if (code == MARKERS)
	strcpy(cmdstr, "markers");
    else if (code == EXPORT)
	strcpy(cmdstr, "export");
    else if (code == SMPLBOOL || code == SMPLDUM)
	strcpy(cmdstr, "sampling");
    else if (code == PANEL) 
	strcpy(cmdstr, "panel");
    else if (code == COMPACT) 
	strcpy(cmdstr, "compact");
    else if (code < NC)
	strcpy(cmdstr, commands[code]);
    else
	pos = -1;

    if (!pos) pos = help_index(cmdstr);
    if (pos == -1) {
	errbox("Sorry, no help is available on this topic");
	return;
    }
    if (help_view == NULL) {
	windata_t *vwin = helpwin(0);

	if (vwin != NULL) help_view = vwin->w;
	gtk_signal_connect(GTK_OBJECT(help_view), "destroy",
			   GTK_SIGNAL_FUNC(gtk_widget_destroyed),
			   &help_view);	
    } else {
	gdk_window_show(help_view->parent->window);
	gdk_window_raise(help_view->parent->window);
    }

    gtk_adjustment_set_value(GTK_TEXT(help_view)->vadj, 
			     (gfloat) pos * 
			     GTK_TEXT(help_view)->vadj->upper / help_length);
}

/* ........................................................... */

static void buf_edit_save (GtkWidget *widget, gpointer data)
{
    windata_t *mydata = (windata_t *) data;
    gchar *text;
    char **pbuf = (char **) mydata->data;

    text = gtk_editable_get_chars(GTK_EDITABLE(mydata->w), 0, -1);
    if (text != NULL && strlen(text) > 0) {
	free(*pbuf); 
	*pbuf = text;
	infobox("Data info saved");
	data_status |= MODIFIED_DATA;
    } else if (strlen(text))
	g_free(text);
}

/* ........................................................... */

static void file_viewer_save (GtkWidget *widget, windata_t *mydata)
{
    /* special case: a newly created script */
    if (strstr(mydata->fname, "script_tmp") || !strlen(mydata->fname)) {
	file_save(mydata, SAVE_SCRIPT, NULL);
	strcpy(mydata->fname, scriptfile);
    } else {
	char buf[MAXLEN];
	FILE *fp;
	gchar *text;

	if ((fp = fopen(mydata->fname, "w")) == NULL) {
	    errbox("Can't open file for writing");
	    return;
	} else {
	    text = gtk_editable_get_chars (GTK_EDITABLE(mydata->w), 0, -1);
	    fprintf(fp, "%s", text);
	    fclose(fp);
	    g_free(text);
	    sprintf(buf, "Saved %s\n", mydata->fname);
	    infobox(buf);
	}
    }
} 

/* .................................................................. */

void windata_init (windata_t *mydata)
{
    mydata->listbox = NULL;
    mydata->mbar = NULL;
    mydata->w = NULL;
    mydata->status = NULL;
    mydata->popup = NULL;
    mydata->ifac = NULL;
    mydata->data = NULL;
    mydata->fname[0] = '\0';
    mydata->action = -1;
    mydata->active_var = 0;
}

/* .................................................................. */

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *mydata = (windata_t *) data;

    if (mydata) {
	if (mydata->listbox) 
	    gtk_widget_destroy(GTK_WIDGET(mydata->listbox));
	if (mydata->mbar) 
	    gtk_widget_destroy(GTK_WIDGET(mydata->mbar));
	if (mydata->status) 
	    gtk_widget_destroy(GTK_WIDGET(mydata->status));
	if (mydata->ifac) 
	    gtk_object_unref(GTK_OBJECT(mydata->ifac));  
	if (mydata->popup) 
	    gtk_object_unref(GTK_OBJECT(mydata->popup));
	if (mydata->action == SUMMARY || mydata->action == VAR_SUMMARY)
	    free_summary(mydata->data); 
	if (mydata->action == CORR)
	    free_corrmat(mydata->data); 
	free(mydata);
	mydata = NULL;
    }
}

/* ........................................................... */

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			char *title, int action,
			GtkItemFactoryEntry menu_items[]) 
{
    GtkWidget *dialog, *close, *table;
    GtkWidget *vscrollbar; 
    windata_t *vwin;

    if ((vwin = mymalloc(sizeof *vwin)) == NULL) return NULL;
    windata_init(vwin);
    vwin->action = action;

    hsize *= gdk_char_width(fixed_font, 'W');
    hsize += 48;

    dialog = gtk_dialog_new();
    gtk_widget_set_usize (dialog, hsize, vsize);
    gtk_window_set_title(GTK_WINDOW(dialog), title);
    gtk_container_border_width (GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), 5);
    gtk_container_border_width 
	(GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dialog)->action_area), TRUE);
#ifndef G_OS_WIN32
    gtk_signal_connect_after(GTK_OBJECT(dialog), "realize", 
			     GTK_SIGNAL_FUNC(set_wm_icon), 
			     NULL);
#endif

    if (menu_items != NULL) {
	set_up_viewer_menu(dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    }

    table = gtk_table_new(1, 2, FALSE);
    gtk_widget_set_usize(table, 500, 400);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       table, TRUE, TRUE, FALSE);

    vwin->w = gtk_text_new(NULL, NULL);

    gtk_text_set_editable(GTK_TEXT(vwin->w), FALSE);

    gtk_text_set_word_wrap(GTK_TEXT(vwin->w), TRUE);
    gtk_table_attach(GTK_TABLE(table), vwin->w, 0, 1, 0, 1,
		     GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND | 
		     GTK_SHRINK, 0, 0);
    gtk_widget_show(vwin->w);

    vscrollbar = gtk_vscrollbar_new(GTK_TEXT (vwin->w)->vadj);
    gtk_table_attach (GTK_TABLE (table), 
		      vscrollbar, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);
    gtk_widget_show (vscrollbar);

    gtk_widget_show(table);

    /* close button */
    close = gtk_button_new_with_label("Close");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       close, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(close), "clicked", 
		       GTK_SIGNAL_FUNC(delete_file_viewer), 
		       (gpointer) dialog);
    gtk_widget_show(close);

    /* insert and then free the text buffer */
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, prn->buf, 
		    strlen(prn->buf));
    gretl_print_destroy(prn);
    
    /* clean up when dialog is destroyed */
    gtk_signal_connect(GTK_OBJECT(dialog), "key_press_event", 
			   GTK_SIGNAL_FUNC(catch_key), dialog);

    gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
		       GTK_SIGNAL_FUNC(free_windata), vwin);

    gtk_widget_show(dialog);
    return vwin;
}

/* ........................................................... */

windata_t *view_file (char *filename, int editable, int del_file, 
		      int hsize, int vsize, char *title, 
		      GtkItemFactoryEntry menu_items[]) 
{
    GtkWidget *dialog, *close, *save = NULL, *table;
    GtkWidget *vscrollbar; 
    extern GdkColor red, blue;
    void *colptr = NULL, *nextcolor = NULL;
    char tempstr[MAXSTR], *fle = NULL;
    FILE *fd = NULL;
    windata_t *vwin;
    int console = 0;
    static GtkStyle *style;

    fd = fopen(filename, "r");
    if (fd == NULL) {
	sprintf(tempstr, "Can't open %s for reading", filename);
	errbox(tempstr);
	return NULL;
    }

    if ((vwin = mymalloc(sizeof *vwin)) == NULL)
	return NULL;
    windata_init(vwin);
    strcpy(vwin->fname, filename);

    hsize *= gdk_char_width(fixed_font, 'W');
    hsize += 48;

    dialog = gtk_dialog_new();
    gtk_widget_set_usize (dialog, hsize, vsize);
    gtk_window_set_title(GTK_WINDOW(dialog), title);
    gtk_container_border_width (GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), 5);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dialog)->action_area), TRUE);
#ifndef G_OS_WIN32
    gtk_signal_connect_after(GTK_OBJECT(dialog), "realize", 
			     GTK_SIGNAL_FUNC(set_wm_icon), 
			     NULL);
#endif

    if (menu_items != NULL) {
	set_up_viewer_menu(dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    }

    table = gtk_table_new(1, 2, FALSE);
    gtk_widget_set_usize(table, 500, 400);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       table, TRUE, TRUE, FALSE);

    vwin->w = gtk_text_new(NULL, NULL);

    if (style == NULL) {
	style = gtk_style_new();
	gdk_font_unref(style->font);
	style->font = fixed_font;
    }
    gtk_widget_set_style(GTK_WIDGET(vwin->w), style);

    if (editable) 
	gtk_text_set_editable(GTK_TEXT(vwin->w), TRUE);
    else 
	gtk_text_set_editable(GTK_TEXT(vwin->w), FALSE);

    /* special case: the gretl console */
    if (strcmp(title, "gretl console") == 0) console = 1;
    if (console) {
	gtk_signal_connect(GTK_OBJECT(vwin->w), "key_press_event",
			   (GtkSignalFunc) console_handler, NULL);
    }

    gtk_text_set_word_wrap(GTK_TEXT(vwin->w), TRUE);
    gtk_table_attach(GTK_TABLE(table), vwin->w, 0, 1, 0, 1,
		     GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND | 
		     GTK_SHRINK, 0, 0);
    gtk_widget_show(vwin->w);


    vscrollbar = gtk_vscrollbar_new (GTK_TEXT (vwin->w)->vadj);
    gtk_table_attach (GTK_TABLE (table), 
		      vscrollbar, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);
    gtk_widget_show (vscrollbar);

    gtk_widget_show(table);

    /* is the file to be deleted after viewing? */
    if (del_file) {
	if ((fle = mymalloc(strlen(filename) + 1)) == NULL)
	    return NULL;
	strcpy(fle, filename);
    }

    /* add a "save" button for editable files */
    if (editable && !console)  {
	save = gtk_button_new_with_label("Save");
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
			   save, FALSE, TRUE, 0);
	gtk_signal_connect(GTK_OBJECT(save), "clicked", 
			   GTK_SIGNAL_FUNC(file_viewer_save), 
			   vwin);
	gtk_widget_show(save);
    }

    /* close button for all uses */
    close = gtk_button_new_with_label("Close");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       close, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(close), "clicked", 
		       GTK_SIGNAL_FUNC(delete_file_viewer), 
		       (gpointer) dialog);
    gtk_widget_show(close);

    /* insert the file text */
    memset(tempstr, 0, sizeof tempstr);
    while (fgets(tempstr, sizeof tempstr - 1, fd)) {
	if (tempstr[0] == '?') 
	    colptr = (console)? &red : &blue;
	if (tempstr[0] == '#') {
	    tempstr[0] = ' ';
	    nextcolor = &red;
	} else
	    nextcolor = NULL;
	gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
			colptr, NULL, tempstr, 
			strlen(tempstr));
	colptr = nextcolor;
	memset(tempstr, 0, sizeof tempstr);
    }
    fclose(fd);

    /* clean up when dialog is destroyed */
    if (del_file) {
	gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
			   GTK_SIGNAL_FUNC(delete_file), (gpointer) fle);
    }
    if (!editable) {
	gtk_signal_connect(GTK_OBJECT(dialog), "key_press_event", 
			   GTK_SIGNAL_FUNC(catch_key), dialog);
    }
    gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
		       GTK_SIGNAL_FUNC(free_windata), vwin);

    gtk_widget_show(dialog);

    return vwin;
}

/* ........................................................... */

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, char *title) 
{
    GtkWidget *dialog, *close, *save, *table;
    GtkWidget *vscrollbar; 
    windata_t *vwin;

    if ((vwin = mymalloc(sizeof *vwin)) == NULL)
	return NULL;
    windata_init(vwin);
    vwin->data = pbuf;
    vwin->action = EDIT_BUFFER;

    hsize *= gdk_char_width(fixed_font, 'W');
    hsize += 48;

    dialog = gtk_dialog_new();
    gtk_widget_set_usize (dialog, hsize, vsize);
    gtk_window_set_title(GTK_WINDOW(dialog), title);
    gtk_container_border_width (GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), 5);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dialog)->action_area), TRUE);
#ifndef G_OS_WIN32
    gtk_signal_connect_after(GTK_OBJECT(dialog), "realize", 
			     GTK_SIGNAL_FUNC(set_wm_icon), 
			     NULL);
#endif

    /* add a menu bar */
    set_up_viewer_menu(dialog, vwin, edit_items);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    table = gtk_table_new(1, 2, FALSE);
    gtk_widget_set_usize(table, 500, 400);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       table, TRUE, TRUE, FALSE);

    vwin->w = gtk_text_new(NULL, NULL);

    gtk_text_set_editable(GTK_TEXT(vwin->w), TRUE);
    gtk_text_set_word_wrap(GTK_TEXT(vwin->w), TRUE);

    gtk_table_attach(GTK_TABLE(table), vwin->w, 0, 1, 0, 1,
		     GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND | 
		     GTK_SHRINK, 0, 0);
    gtk_widget_show(vwin->w);

    vscrollbar = gtk_vscrollbar_new (GTK_TEXT (vwin->w)->vadj);
    gtk_table_attach (GTK_TABLE (table), 
		      vscrollbar, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);
    gtk_widget_show (vscrollbar);

    gtk_widget_show(table);

    /* add a "save" button */
    save = gtk_button_new_with_label("Save");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       save, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(save), "clicked", 
		       GTK_SIGNAL_FUNC(buf_edit_save), 
		       vwin);
    gtk_widget_show(save);

    /* and a close button */
    close = gtk_button_new_with_label("Close");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       close, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(close), "clicked", 
		       GTK_SIGNAL_FUNC(delete_file_viewer), 
		       (gpointer) dialog);
    gtk_widget_show(close);

    /* insert the buffer text */
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, *pbuf, strlen(*pbuf));

    /* clean up when dialog is destroyed */
    gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
		       GTK_SIGNAL_FUNC(free_windata), vwin);

    gtk_widget_show(dialog);

    return vwin;
}

/* ........................................................... */

void flip (GtkItemFactory *ifac, char *path, gboolean s)
{
    if (ifac != NULL) {
	GtkWidget *w = gtk_item_factory_get_item(ifac, path);

	if (w != NULL) 
	    gtk_widget_set_sensitive(w, s);
	else
	    fprintf(stderr, "Failed to flip state of \"%s\"\n", path);
    }
}

/* ........................................................... */

static void model_panel_menu_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/Tests/panel diagnostics", s);
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
    flip(ifac, "/Model data/Add to data set/error sum of squares", s);
    flip(ifac, "/Model data/Add to data set/standard error of residuals", s);
    flip(ifac, "/Model data/Add to data set/R-squared", s);
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
    GtkAccelGroup *accel;
    gint n_items = 0;

    while (items[n_items].path != NULL) n_items++;

    accel = gtk_accel_group_new();
    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", 
				      accel);
    gtk_item_factory_create_items(vwin->ifac, n_items, items, vwin);
    vwin->mbar = gtk_item_factory_get_widget(vwin->ifac, "<main>");
    gtk_accel_group_attach(accel, GTK_OBJECT (window));

    if (vwin->action == SUMMARY || vwin->action == VAR_SUMMARY
	|| vwin->action == CORR) {
	augment_copy_menu(vwin);
	return;
    }

    if (vwin->action == VIEW_MODEL && vwin->data != NULL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->ci == POOLED) 
	    model_panel_menu_state(vwin->ifac, TRUE);
	else model_panel_menu_state(vwin->ifac, FALSE);
	if (pmod->ci != OLS && pmod->ci != POOLED) { 
	    lmmenu_state(vwin->ifac, FALSE);
	    latex_menu_state(vwin->ifac, FALSE);
	}
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    model_menu_state(vwin->ifac, FALSE);
	}
	if (pmod->name)
	    model_save_state(vwin->ifac, FALSE);
    }
}

/* .................................................................. */

static void add_vars_to_plot_menu (windata_t *vwin)
{
    int i, j;
    GtkItemFactoryEntry varitem;
    gchar *mpath[] = {"/Graphs/residual plot", 
		      "/Graphs/fitted, actual plot"};
    MODEL *pmod = vwin->data;

    varitem.path = NULL;

    for (i=0; i<2; i++) {
	varitem.path = mymalloc(64);
	varitem.accelerator = NULL;
	varitem.callback_action = 0; 
	varitem.item_type = NULL;
	if (dataset_is_time_series(datainfo))
	    sprintf(varitem.path, "%s/against time", mpath[i]);
	else
	    sprintf(varitem.path, "%s/by observation number", mpath[i]);
	if (i == 0)
	    varitem.callback = resid_plot; 
	else
	    varitem.callback = fit_actual_plot;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);

	/* put the indep vars on the menu list */
	for (j=2; j<pmod->list[0]; j++) {
	    if (pmod->list[j] == 0) continue;
	    if (varitem.path == NULL)
		varitem.path = mymalloc(64);
	    varitem.accelerator = NULL;
	    varitem.callback_action = pmod->list[j]; 
	    varitem.item_type = NULL;
	    sprintf(varitem.path, "%s/against %s", mpath[i], 
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
    windata_t *mydata = (windata_t *) data;

    if (item->active) mydata->active_var = v; 
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
	if (!isdummy(pmod->list[i], datainfo->t1, datainfo->t2, 
		     Z, datainfo->n))
	    continue;
	if (!dums) { /* add separator, branch and "none" */
	    dumitem.path = mymalloc(64);
	    sprintf(dumitem.path, "/Graphs/dumsep");
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Separator>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* menu branch */
	    sprintf(dumitem.path, "/Graphs/Separation");
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Branch>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    /* "none" option */
	    sprintf(dumitem.path, "/Graphs/Separation/none");
	    dumitem.callback = plot_dummy_call;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<RadioItem>";
	    dumitem.accelerator = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    dums = 1;
	} 
	dumitem.callback_action = pmod->list[i]; 
	sprintf(dumitem.path, "/Graphs/Separation/by %s",  
		datainfo->varname[pmod->list[i]]);
	dumitem.callback = plot_dummy_call;	    
	dumitem.accelerator = NULL;
	dumitem.item_type = "/Graphs/Separation/none";
	gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
    }
    free(dumitem.path);
}

/* ........................................................... */

static void check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data)
{
    windata_t *mwin = (windata_t *) data;
    MODEL *pmod = mwin->data;
    extern int quiet_sample_check (MODEL *pmod);
    int s, ok = 1;

    if (quiet_sample_check(pmod)) ok = 0;
    s = GTK_WIDGET_IS_SENSITIVE
	(gtk_item_factory_get_item(mwin->ifac, "/Tests/omit variables"));
    if ((s && ok) || (!s && !ok)) return;
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
}

/* ........................................................... */

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title) 
{
    windata_t *vwin;
    GtkWidget *dialog, *close, *table, *scroller;

    if ((vwin = mymalloc(sizeof *vwin)) == NULL) return 1;
    windata_init(vwin);

    hsize *= gdk_char_width (fixed_font, 'W');
    hsize += 48;

    vwin->data = pmod;
    vwin->action = VIEW_MODEL;
    dialog = gtk_dialog_new();
    gtk_widget_set_usize (dialog, hsize, vsize);
    gtk_window_set_title(GTK_WINDOW(dialog), title);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG(dialog)->vbox), 5);
    gtk_container_border_width 
	(GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    gtk_box_set_homogeneous(GTK_BOX 
			    (GTK_DIALOG(dialog)->action_area), TRUE);
#ifndef G_OS_WIN32
    gtk_signal_connect_after(GTK_OBJECT(dialog), "realize", 
			     GTK_SIGNAL_FUNC(set_wm_icon), 
			     NULL);
#endif

    set_up_viewer_menu(dialog, vwin, model_items);

    /* add menu of indep vars, against which to plot resid */
    add_vars_to_plot_menu(vwin);
    add_dummies_to_plot_menu(vwin);
    gtk_signal_connect(GTK_OBJECT(vwin->mbar), "button_press_event", 
		       GTK_SIGNAL_FUNC(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    table = gtk_table_new(1, 2, FALSE);
    gtk_widget_set_usize(table, hsize, vsize); 
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), 
		       table, TRUE, TRUE, FALSE);

    vwin->w = gtk_text_new(NULL, NULL);
    gtk_text_set_editable(GTK_TEXT(vwin->w), FALSE);
    gtk_text_set_word_wrap(GTK_TEXT(vwin->w), TRUE);
    gtk_table_attach(GTK_TABLE(table), vwin->w, 0, 1, 0, 1,
		     GTK_FILL | GTK_EXPAND, GTK_FILL | GTK_EXPAND | 
		     GTK_SHRINK, 0, 0);
    gtk_widget_show(vwin->w);
    scroller = gtk_vscrollbar_new(GTK_TEXT(vwin->w)->vadj);
    gtk_table_attach (GTK_TABLE(table), 
		      scroller, 1, 2, 0, 1,
		      GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);
    gtk_widget_show(scroller);

    gtk_widget_show(table);

    /* close button */
    close = gtk_button_new_with_label("Close");
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
		       close, FALSE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(close), "clicked", 
		       GTK_SIGNAL_FUNC(delete_file_viewer), 
		       (gpointer) dialog);
    gtk_widget_show(close);

    /* insert and then free the model buffer */
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, prn->buf, 
		    strlen(prn->buf));
    gretl_print_destroy(prn);

    copylist(&default_list, pmod->list);

    /* attach shortcuts */
    gtk_object_set_data(GTK_OBJECT(dialog), "ddata", vwin);
    gtk_signal_connect(GTK_OBJECT(dialog), "key_press_event", 
		       GTK_SIGNAL_FUNC(catch_key), 
		       dialog);
    gtk_signal_connect(GTK_OBJECT(dialog), "key_press_event", 
		       GTK_SIGNAL_FUNC(catch_key), 
		       dialog);

    /* clean up when dialog is destroyed */
    gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
		       GTK_SIGNAL_FUNC(delete_model), 
		       vwin->data);
    gtk_signal_connect(GTK_OBJECT(dialog), "destroy", 
		       GTK_SIGNAL_FUNC(free_windata), 
		       vwin);

    gtk_widget_show_all(dialog);
    return 0;
}

/* ......................................................... */

void setup_column (GtkWidget *listbox, int column, int width) 
{
    if (width == 0) 
	gtk_clist_set_column_auto_resize (GTK_CLIST (listbox), column, TRUE);
    else if (width == -1) 
	gtk_clist_set_column_visibility (GTK_CLIST (listbox), column, FALSE);
    else 
	gtk_clist_set_column_width (GTK_CLIST (listbox), column, width);
}

/* ........................................................... */

#ifdef USE_GNOME

static void msgbox (const char *msg, int err)
{
    if (err) gnome_app_error(GNOME_APP(mdata->w), msg);
    else gnome_app_message(GNOME_APP(mdata->w), msg);
}

#else
#ifdef G_OS_WIN32

static void msgbox (const char *msg, int err)
{
    if (err) 
	MessageBox(NULL, msg, "gretl", MB_OK | MB_ICONERROR);
    else
	MessageBox(NULL, msg, "gretl", MB_OK | MB_ICONINFORMATION);
}

#else /* win32 */

static void msgbox (const char *msg, int err) 
{
    GtkWidget *w, *label, *button, *table;
    char labeltext[MAXLEN];

    if (err)
	sprintf(labeltext, "Error:\n%s\n", msg);
    else
	sprintf(labeltext, "Info:\n%s\n", msg);
    w = gtk_window_new(GTK_WINDOW_DIALOG);
    gtk_container_border_width(GTK_CONTAINER(w), 5);
    gtk_window_position (GTK_WINDOW(w), GTK_WIN_POS_MOUSE);
    gtk_window_set_title (GTK_WINDOW (w), (err)? "gretl error" : 
			  "gretl info");  
  
    table = gtk_table_new(2, 3, FALSE);
    gtk_container_add(GTK_CONTAINER(w), table);
  
    label = gtk_label_new(labeltext);
    gtk_table_attach_defaults(GTK_TABLE(table), label, 0, 3, 0, 1);

    if (err)
	button = gtk_button_new_with_label("Close");
    else
	button = gtk_button_new_with_label("OK");
    gtk_table_attach_defaults(GTK_TABLE(table), button, 1, 2, 1, 2);
  
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(delete_widget), w);
    gtk_widget_show(button);
    gtk_widget_show(label);
    gtk_widget_show(table);
    gtk_widget_show(w);  
}

#endif
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
    
    if (n > 8) {
	safecpy(namebit, varname, 8);
	sprintf(errtext, "Variable name %s... is too long\n"
	       "(the max is 8 characters)", namebit);
	errbox(errtext);
	return 1;
    }
    if (!(isalpha(varname[0]))) {
	sprintf(errtext, "First char of name ('%c') is bad\n"
	       "(first must be alphabetical)", varname[0]);
	errbox(errtext);
	return 1;
    }
    for (i=1; i<n; i++) {
	if (!(isalpha(varname[i]))  
	    && !(isdigit(varname[i]))
	    && varname[i] != '_') {
	    sprintf(errtext, "Name contains an illegal char (in place %d)\n"
		    "Use only letters, digits and underscore", i + 1);
	    errbox(errtext);
	    return 1;
	}
    }
    return 0;
}	

/* .................................................................. */

void options_dialog (gpointer data) 
{
    GtkWidget *tempwid, *dialog, *notebook;

    dialog = gtk_dialog_new ();
    gtk_window_set_title (GTK_WINDOW (dialog), "gretl: options");
    gtk_container_border_width 
	(GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_border_width 
	(GTK_CONTAINER (GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 2);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX (GTK_DIALOG (dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect_object 
	(GTK_OBJECT (dialog), "delete_event", GTK_SIGNAL_FUNC 
	 (gtk_widget_destroy), GTK_OBJECT (dialog));
   
    notebook = gtk_notebook_new ();
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (dialog)->vbox), 
			notebook, TRUE, TRUE, 0);
    gtk_widget_show (notebook);

    make_prefs_tab (notebook, 1);
    make_prefs_tab (notebook, 2);
    make_prefs_tab (notebook, 3);
    make_prefs_tab (notebook, 4);
    make_prefs_tab (notebook, 5);
   
    tempwid = gtk_button_new_with_label ("OK");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (apply_changes), NULL);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label ("  Cancel  ");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (dialog));
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label ("Apply");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG 
				 (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (apply_changes), NULL);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    gtk_widget_show (dialog);
}

/* .................................................................. */

static void flip_sensitive (GtkWidget *w, gpointer data)
{
    GtkWidget *entry = GTK_WIDGET(data);
    
    gtk_widget_set_sensitive(entry, GTK_TOGGLE_BUTTON(w)->active);
}

/* .................................................................. */

static void make_prefs_tab (GtkWidget *notebook, int tab) 
{
    GtkWidget *box, *inttbl, *chartbl, *tempwid = NULL;
    int i, tbl_len, tbl_num, tbl_col;
    RCVARS *rc = NULL;
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show (box);

    if (tab == 1)
	tempwid = gtk_label_new ("General");
    else if (tab == 2)
	tempwid = gtk_label_new ("Databases");
    else if (tab == 3)
	tempwid = gtk_label_new ("Toolbar");
    else if (tab == 4)
	tempwid = gtk_label_new ("Open/Save path");
    else if (tab == 5)
	tempwid = gtk_label_new ("Data files");
    
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    chartbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (chartbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (chartbl), 5);
    gtk_box_pack_start (GTK_BOX (box), chartbl, FALSE, FALSE, 0);
    gtk_widget_show (chartbl);
   
    tbl_num = tbl_col = 0;
    inttbl = gtk_table_new (1, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (inttbl), 2);
    gtk_table_set_col_spacings (GTK_TABLE (inttbl), 5);
    gtk_box_pack_start (GTK_BOX (box), inttbl, FALSE, FALSE, 0);
    gtk_widget_show (inttbl);

    i = 0;
    while (rc_vars[i].key != NULL) {
	rc = &rc_vars[i];
	if (rc->tab == tab) {
	    if (rc->type == 'B' 
		&& rc->link == NULL) { /* simple boolean variable */
		tempwid = gtk_check_button_new_with_label 
		    (rc->description);
		gtk_table_attach_defaults 
		    (GTK_TABLE (inttbl), tempwid, tbl_col, tbl_col + 1, 
		     tbl_num, tbl_num + 1);
		if (*(int *)(rc->var))
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON (tempwid), TRUE);
		else
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON (tempwid), FALSE);
		/* special case: link between toggle and preceding entry */
		if (rc->len) {
		    gtk_widget_set_sensitive(rc_vars[i-1].widget,
					     GTK_TOGGLE_BUTTON(tempwid)->active);
		    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
				       GTK_SIGNAL_FUNC(flip_sensitive),
				       rc_vars[i-1].widget);
		}
		/* end link to entry */
		gtk_widget_show (tempwid);
		rc->widget = tempwid;
		tbl_col++;
		if (tbl_col == 2) {
		    tbl_col = 0;
		    tbl_num++;
		    gtk_table_resize (GTK_TABLE (inttbl), tbl_num + 1, 2);
		}
	    } 
	    else if (rc->type == 'B') { /* radio-button dichotomy */
		int val = *(int *)(rc->var);
		GSList *group;

		tbl_num += 2;
		gtk_table_resize (GTK_TABLE(inttbl), tbl_num + 1, 2);

		tempwid = gtk_radio_button_new_with_label(NULL, 
							  rc->description);
		gtk_table_attach_defaults 
		    (GTK_TABLE (inttbl), tempwid, tbl_col, tbl_col + 1, 
		     tbl_num - 2, tbl_num - 1);    
		if (val) 
		    gtk_toggle_button_set_active 
			(GTK_TOGGLE_BUTTON(tempwid), TRUE);
		gtk_widget_show (tempwid);
		rc->widget = tempwid;
		group = gtk_radio_button_group(GTK_RADIO_BUTTON(tempwid));
		tempwid = gtk_radio_button_new_with_label(group, rc->link);
		gtk_table_attach_defaults 
		    (GTK_TABLE (inttbl), tempwid, tbl_col, tbl_col + 1, 
		     tbl_num - 1, tbl_num);  
		if (!val)
		    gtk_toggle_button_set_active
			(GTK_TOGGLE_BUTTON(tempwid), TRUE);
		gtk_widget_show (tempwid);
	    } else { /* string variable */
		tbl_len++;
		gtk_table_resize (GTK_TABLE (chartbl), tbl_len, 2);
		tempwid = gtk_label_new (rc->description);
		gtk_misc_set_alignment (GTK_MISC (tempwid), 1, 0.5);
		gtk_table_attach_defaults (GTK_TABLE (chartbl), 
					   tempwid, 0, 1, tbl_len-1, tbl_len);
		gtk_widget_show (tempwid);

		tempwid = gtk_entry_new ();
		gtk_table_attach_defaults (GTK_TABLE (chartbl), 
					   tempwid, 1, 2, tbl_len-1, tbl_len);
		gtk_entry_set_text (GTK_ENTRY (tempwid), rc->var);
		gtk_widget_show (tempwid);
		rc->widget = tempwid;
	    } 
	}
	i++;
    }
}

/* .................................................................. */

static void apply_changes (GtkWidget *widget, gpointer data) 
{
    gchar *tempstr;
    extern void show_toolbar (void);
    int i = 0;

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].widget != NULL) {
	    if (rc_vars[i].type == 'B') {
		if (GTK_TOGGLE_BUTTON(rc_vars[i].widget)->active)
		    *(int *)(rc_vars[i].var) = TRUE;
		else *(int *)(rc_vars[i].var) = FALSE;
	    } 
	    if (rc_vars[i].type == 'U' || rc_vars[i].type == 'R') {
		tempstr = gtk_entry_get_text
		    (GTK_ENTRY(rc_vars[i].widget));
		if (tempstr != NULL && strlen(tempstr)) 
		    strncpy(rc_vars[i].var, tempstr, rc_vars[i].len - 1);
	    }
	}
	i++;
    }
    write_rc();
    if (toolbar_box == NULL && want_toolbar)
	show_toolbar();
    else if (toolbar_box != NULL && !want_toolbar) {
	gtk_widget_destroy(toolbar_box);
	toolbar_box = NULL;
    }
    proxy_init(dbproxy);
}

/* .................................................................. */

static void str_to_boolvar (char *s, void *b)
{
    if (strcmp(s, "true") == 0 || strcmp(s, "1") == 0)
	*(int *)b = TRUE;
    else
	*(int *)b = FALSE;	
}

/* .................................................................. */

static void boolvar_to_str (void *b, char *s)
{
    if (*(int *)b) strcpy(s, "true");
    else strcpy(s, "false");
}

/* .................................................................. */

#if defined(USE_GNOME)

void write_rc (void) 
{
    char gpath[MAXSTR];
    char val[6];
    int i = 0;

    while (rc_vars[i].key != NULL) {
	sprintf(gpath, "/gretl/%s/%s", rc_vars[i].description, rc_vars[i].key);
	if (rc_vars[i].type == 'B') {
	    boolvar_to_str(rc_vars[i].var, val);
	    gnome_config_set_string(gpath, val);
	} else
	    gnome_config_set_string(gpath, rc_vars[i].var);
	i++;
    }
    gnome_printfilelist(1); /* data files */
    gnome_printfilelist(2); /* session files */
    gnome_printfilelist(3); /* script files */    
    gnome_config_sync();
    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    int i = 0;
    gchar *value = NULL;
    char gpath[MAXSTR];

    while (rc_vars[i].key != NULL) {
	sprintf(gpath, "/gretl/%s/%s", 
		rc_vars[i].description, 
		rc_vars[i].key);
	if ((value = gnome_config_get_string(gpath)) != NULL) {
	    if (rc_vars[i].type == 'B')
		str_to_boolvar(value, rc_vars[i].var);
	    else
		strncpy(rc_vars[i].var, value, rc_vars[i].len - 1);
	    g_free(value);
	}
	i++;
    }

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
    /* get recent file lists */
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent data files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(datalist[i], value);
	    g_free(value);
	}
	else break;
    }    
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent session files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(sessionlist[i], value);
	    g_free(value);
	}
	else break;
    } 
    for (i=0; i<MAXRECENT; i++) {
	sprintf(gpath, "/gretl/recent script files/%d", i);
	if ((value = gnome_config_get_string(gpath)) != NULL) { 
	    strcpy(scriptlist[i], value);
	    g_free(value);
	}
	else break;
    }
    set_paths(&paths, 0, 1); /* 0 = not defaults, 1 = gui */
}

/* end of gnome versions, now win32 */
#elif defined(G_OS_WIN32)

void write_rc (void) 
{
    int i = 0;
    char val[6];

    while (rc_vars[i].key != NULL) {
	if (rc_vars[i].type == 'B') {
	    boolvar_to_str(rc_vars[i].var, val);
	    write_reg_val(HKEY_CURRENT_USER, rc_vars[i].key, val);
	} else
	    write_reg_val((rc_vars[i].type == 'R')? 
			  HKEY_CLASSES_ROOT : HKEY_CURRENT_USER, 
			  rc_vars[i].key, rc_vars[i].var);
	i++;
    }
    win_printfilelist(1); /* data files */
    win_printfilelist(2); /* session files */
    win_printfilelist(3); /* script files */
    set_paths(&paths, 0, 1);
}

void read_rc (void) 
{
    int i = 0;
    char rpath[MAXSTR], value[MAXSTR];

    while (rc_vars[i].key != NULL) {
	if (read_reg_val((rc_vars[i].type == 'R')? 
			 HKEY_CLASSES_ROOT : HKEY_CURRENT_USER, 
			 rc_vars[i].key, value) == 0) {
	    if (rc_vars[i].type == 'B') {
		str_to_boolvar(value, rc_vars[i].var);
	    } else
		strncpy(rc_vars[i].var, value, rc_vars[i].len - 1);
	}
	i++;
    }

    /* initialize lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
    /* get recent file lists */
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent data files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, rpath, value) == 0) 
	    strcpy(datalist[i], value);
	else break;
    }    
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent session files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, rpath, value) == 0) 
	    strcpy(sessionlist[i], value);
	else break;
    } 
    for (i=0; i<MAXRECENT; i++) {
	sprintf(rpath, "recent script files\\%d", i);
	if (read_reg_val(HKEY_CURRENT_USER, rpath, value) == 0) 
	    strcpy(scriptlist[i], value);
	else break;
    }
    set_paths(&paths, 0, 1);
}

#else /* end of win32 versions, now plain GTK */

void write_rc (void) 
{
    FILE *rc;
    int i;
    char val[6];

    rc = fopen(rcfile, "w");
    if (rc == NULL) {
	errbox("Couldn't open config file for writing");
	return;
    }
    fprintf(rc, "# config file written by gretl: do not edit\n");
    i = 0;
    while (rc_vars[i].var != NULL) {
	fprintf(rc, "# %s\n", rc_vars[i].description);
	if (rc_vars[i].type == 'B') {
	    boolvar_to_str(rc_vars[i].var, val);
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, val);
	} else
	    fprintf(rc, "%s = %s\n", rc_vars[i].key, rc_vars[i].var);
	i++;
    }
    printfilelist(1, rc); /* data files */
    printfilelist(2, rc); /* session files */
    printfilelist(3, rc); /* script files */
    fclose(rc);
    set_paths(&paths, 0, 1);
}

static void read_rc (void) 
{
    FILE *rc;
    int i, j;
    char line[MAXLEN], key[32], linevar[MAXLEN];
    int gotrecent = 0;

    if ((rc = fopen(rcfile, "r")) == NULL) return;

    i = 0;
    while (rc_vars[i].var != NULL) {
	if (fgets(line, MAXLEN, rc) == NULL) 
	    break;
	if (line[0] == '#') 
	    continue;
	if (!strncmp(line, "recent ", 7)) {
	    gotrecent = 1;
	    break;
	}
	if (sscanf(line, "%s", key) == 1) {
	    strcpy(linevar, line + strlen(key) + 3); 
	    chopstr(linevar); 
	    for (j=0; rc_vars[j].key != NULL; j++) {
		if (!strcmp(key, rc_vars[j].key)) {
		    if (rc_vars[j].type == 'B')
			str_to_boolvar(linevar, rc_vars[j].var);
		    else
			strcpy(rc_vars[j].var, linevar);
		}
	    }
	}
	i++;
    }

    /* get lists of recently opened files */
    for (i=0; i<MAXRECENT; i++) { 
	datalist[i][0] = 0;
	sessionlist[i][0] = 0;
	scriptlist[i][0] = 0;
    }
    if (gotrecent || (fgets(line, MAXLEN, rc) != NULL && 
	strncmp(line, "recent data files:", 18) == 0)) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    if (strncmp(line, "recent session files:", 21) == 0)
		break;
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(datalist[i++], line);
	}
    }
    if (strncmp(line, "recent session files:", 21) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    if (strncmp(line, "recent script files:", 20) == 0)
		break;
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(sessionlist[i++], line);
	}
    }
    if (strncmp(line, "recent script files:", 20) == 0) {
	i = 0;
	while (fgets(line, MAXLEN, rc) && i<MAXRECENT) {
	    chopstr(line);
	    if (strlen(line)) 
		strcpy(scriptlist[i++], line);
	}
    }
    fclose(rc);
    set_paths(&paths, 0, 1);
}

#endif /* end of "plain gtk" versions of read_rc, write_rc */

/* .................................................................. */

static void font_selection_ok (GtkWidget *w, GtkFontSelectionDialog *fs)
{
    gchar *fstring = gtk_font_selection_dialog_get_font_name(fs);

    if (strlen(fstring)) {
	strcpy(fontspec, fstring);
	gdk_font_unref(fixed_font);
	fixed_font = gdk_font_load(fontspec);
	write_rc();
    }
    g_free(fstring);
    gtk_widget_destroy(GTK_WIDGET (fs));
}

/* .................................................................. */

void font_selector (void)
{
    static GtkWidget *fontsel = NULL;
    gchar *spacings[] = { "c", "m", NULL };

    if (!fontsel) {
	fontsel = gtk_font_selection_dialog_new 
	    ("Font for gretl output windows");
	gtk_window_set_position (GTK_WINDOW (fontsel), GTK_WIN_POS_MOUSE);
	gtk_font_selection_dialog_set_filter 
	    (GTK_FONT_SELECTION_DIALOG (fontsel),
                                       GTK_FONT_FILTER_BASE, GTK_FONT_ALL,
                                       NULL, NULL, NULL, NULL, spacings, NULL);
	gtk_font_selection_dialog_set_font_name 
	    (GTK_FONT_SELECTION_DIALOG (fontsel), fontspec);

	gtk_signal_connect (GTK_OBJECT(fontsel), "destroy",
			    GTK_SIGNAL_FUNC(gtk_widget_destroyed),
			    &fontsel);

	gtk_signal_connect (GTK_OBJECT 
			    (GTK_FONT_SELECTION_DIALOG 
			     (fontsel)->ok_button),
			    "clicked", GTK_SIGNAL_FUNC(font_selection_ok),
			    GTK_FONT_SELECTION_DIALOG (fontsel));

	gtk_signal_connect_object (GTK_OBJECT 
				   (GTK_FONT_SELECTION_DIALOG 
				    (fontsel)->cancel_button),
				   "clicked", 
				   GTK_SIGNAL_FUNC(gtk_widget_destroy),
				   GTK_OBJECT (fontsel));
    }
    if (!GTK_WIDGET_VISIBLE (fontsel)) gtk_widget_show (fontsel);
    else gtk_widget_destroy (fontsel);
}

/* .................................................................. */

static void close_find_dialog (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy (widget);
    find_window = NULL;
}

/* .................................................................. */

static void find_in_help (GtkWidget *widget, gpointer data)
{
    int nIndex = 0, i, linecount = 0;
    char *haystack;   

    haystack = gtk_editable_get_chars(GTK_EDITABLE(help_view), 0,
	gtk_text_get_length(GTK_TEXT(help_view)));

    if (needle) g_free(needle);

    needle = gtk_editable_get_chars(GTK_EDITABLE (find_entry), 0, -1);
    nIndex = GTK_EDITABLE (help_view)->selection_end_pos;

    nIndex = look_for_string(haystack, needle, nIndex);

    if (nIndex >= 0) {
	gtk_text_freeze(GTK_TEXT(help_view));
        gtk_text_set_point (GTK_TEXT (help_view), nIndex);
        gtk_text_insert (GTK_TEXT (help_view), NULL, NULL, NULL, " ", 1);
        gtk_text_backward_delete (GTK_TEXT (help_view), 1);
	gtk_text_thaw(GTK_TEXT(help_view));
        gtk_editable_select_region (GTK_EDITABLE (help_view), 
				    nIndex, nIndex + strlen (needle));
	for (i=0; i<nIndex; i++) 
	    if (haystack[i] == '\n') linecount++;
	gtk_adjustment_set_value(GTK_TEXT(help_view)->vadj, 
				 (gfloat) (linecount - 2) *
				 GTK_TEXT(help_view)->vadj->upper / help_length);
	find_window = NULL;
    } else infobox("String was not found.");

    g_free(haystack);
}

/* .................................................................. */

static void find_in_clist (GtkWidget *w, gpointer data)
{
    int start, found = 0, n, i;
    gchar *tmp; 
    char haystack[MAXLEN];
    windata_t *dbdat;

    dbdat = (windata_t *) gtk_object_get_data(GTK_OBJECT(data), "windat");

    if (needle) g_free(needle);
    needle = gtk_editable_get_chars(GTK_EDITABLE (find_entry), 0, -1);
    lower(needle);

    start = dbdat->active_var + 1;
    n = GTK_CLIST(dbdat->listbox)->rows;

    for (i=start; i<n; i++) {  
	gtk_clist_get_text(GTK_CLIST(dbdat->listbox), i, 1, &tmp);
	strcpy(haystack, tmp);
	lower(haystack);
	found = look_for_string(haystack, needle, 0);
	if (found >= 0) break;
    }
    if (found >= 0) {
	gtk_clist_moveto(GTK_CLIST(dbdat->listbox), i, 0, 0, .1);
	gtk_clist_select_row(GTK_CLIST(dbdat->listbox), i, 0);
	dbdat->active_var = i;
	find_window = NULL;    
    } else {
	gtk_clist_select_row(GTK_CLIST(dbdat->listbox), 0, 0);
	dbdat->active_var = 0;
	infobox("String was not found.");
    }
}

/* .................................................................. */

static int look_for_string (char *haystack, char *needle, int start)
{
    int pos;
    int HaystackLength = strlen(haystack);
    int NeedleLength = strlen(needle);

    for (pos = start; pos < HaystackLength; pos++) {
        if (strncmp(&haystack[pos], needle, NeedleLength) == 0) 
             return pos;
    }
    return -1;
}

/* .................................................................. */
 
static void cancel_find (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
    find_window = NULL;
}

/* .................................................................. */

static void find_string_dialog (void (*YesFunc)(), void (*NoFunc)(),
				gpointer data)
{
    GtkWidget *label;
    GtkWidget *button;
    GtkWidget *hbox;
    windata_t *mydat = (windata_t *) data;

    if (find_window) {
	gtk_object_set_data(GTK_OBJECT(find_window), "windat", mydat); 
	return;
    }

    find_window = gtk_dialog_new();
    gtk_object_set_data(GTK_OBJECT(find_window), "windat", mydat);

    gtk_signal_connect (GTK_OBJECT (find_window), "destroy",
	                GTK_SIGNAL_FUNC (close_find_dialog),
	                find_window);
    gtk_window_set_title (GTK_WINDOW (find_window), "gretl: find");
    gtk_container_border_width (GTK_CONTAINER (find_window), 5);

    hbox = gtk_hbox_new(TRUE, TRUE);
    label = gtk_label_new(" Find what:");
    gtk_widget_show (label);
    find_entry = gtk_entry_new();

    if (needle) {
	gtk_entry_set_text(GTK_ENTRY (find_entry), needle);
	gtk_entry_select_region (GTK_ENTRY (find_entry), 0, 
				 strlen (needle));
    }
    gtk_signal_connect(GTK_OBJECT (find_entry), 
			"activate", 
			GTK_SIGNAL_FUNC (YesFunc),
	                find_window);
    gtk_widget_show (find_entry);

    gtk_box_pack_start (GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(hbox), find_entry, TRUE, TRUE, 0);
    gtk_widget_show (hbox);

    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->vbox), 
                        hbox, TRUE, TRUE, 0);

    gtk_box_set_spacing(GTK_BOX (GTK_DIALOG (find_window)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX 
			     (GTK_DIALOG (find_window)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW (find_window), GTK_WIN_POS_MOUSE);

    /* find button -- make this the default */
    button = gtk_button_new_with_label ("Find next");
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (YesFunc), find_window);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* cancel button */
    button = gtk_button_new_with_label ("Cancel");
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG (find_window)->action_area), 
		       button, TRUE, TRUE, FALSE);
    gtk_signal_connect(GTK_OBJECT (button), "clicked",
		       GTK_SIGNAL_FUNC (NoFunc), find_window);
    gtk_widget_show(button);

    gtk_widget_grab_focus(find_entry);
    gtk_widget_show (find_window);
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

void text_copy (gpointer data, guint how, GtkWidget *widget) 
{
    windata_t *mydata = (windata_t *) data;
    PRN *prn;

    /* descriptive statistics */
    if ((mydata->action == SUMMARY || mydata->action == VAR_SUMMARY)
	&& (how == COPY_LATEX || how == COPY_RTF)) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) mydata->data;
	
	if (bufopen(&prn)) return;
	if (how == COPY_LATEX) {
	    texprint_summary(summ, datainfo, prn);
	    prn_to_clipboard(prn);
	} else {
	    rtfprint_summary(summ, datainfo, prn);
#ifdef G_OS_WIN32
	    win_copy_rtf(prn);
#else
	    prn_to_clipboard(prn);
#endif
	}
	gretl_print_destroy(prn);
	return;
    }

    /* correlation matrix */
    if (mydata->action == CORR 
	&& (how == COPY_LATEX || how == COPY_RTF)) {
	CORRMAT *corr = (CORRMAT *) mydata->data;

	if (bufopen(&prn)) return;
	if (how == COPY_LATEX) {
	    texprint_corrmat(corr, datainfo, prn);
	    prn_to_clipboard(prn);
	} else {
	    rtfprint_corrmat(corr, datainfo, prn);
#ifdef G_OS_WIN32
	    win_copy_rtf(prn);
#else
	    prn_to_clipboard(prn);
#endif
	}
	gretl_print_destroy(prn);
	return;
    }

    /* or it's a model window we're copying from? */
    if (how == COPY_RTF) {
	MODEL *pmod = (MODEL *) mydata->data;

	model_to_rtf(pmod);	
	return;
    }
    else if (how == COPY_LATEX || how == COPY_HTML) {
	MODEL *pmod = (MODEL *) mydata->data;

	if (bufopen(&prn)) return;
	if (how == COPY_LATEX)
	    tex_print_model(pmod, datainfo, 0, prn);
	else
	    h_printmodel(pmod, datainfo, prn);
	prn_to_clipboard(prn);
	gretl_print_destroy(prn);
	return;
    }
    else if (how == COPY_LATEX_EQUATION) {
	MODEL *pmod = (MODEL *) mydata->data;

	if (bufopen(&prn)) return;
	tex_print_equation(pmod, datainfo, 0, prn);
	prn_to_clipboard(prn);
	gretl_print_destroy(prn);
	return;
    }

    /* otherwise just copying plain text from plain text window */
    else if (how == COPY_TEXT) {
	gtk_editable_select_region(GTK_EDITABLE(mydata->w), 0, -1);
	gtk_editable_copy_clipboard(GTK_EDITABLE(mydata->w));
    }
    else if (how == COPY_SELECTION) {
	gtk_editable_copy_clipboard(GTK_EDITABLE(mydata->w));
    }
}

/* .................................................................. */

#if defined(G_OS_WIN32) || defined (USE_GNOME)

void window_print (gpointer data, guint u, GtkWidget *widget) 
{
    windata_t *mydata = (windata_t *) data;
    char *buf, *selbuf = NULL;
    GtkEditable *gedit = GTK_EDITABLE(mydata->w);

    buf = gtk_editable_get_chars(gedit, 0, -1);
    if (gedit->has_selection)
	selbuf = gtk_editable_get_chars(gedit, 
					gedit->selection_start_pos,
					gedit->selection_end_pos);
    winprint(buf, selbuf);
}

#endif

/* .................................................................. */

void text_paste (windata_t *mydata, guint u, GtkWidget *widget)
{
    gtk_editable_paste_clipboard(GTK_EDITABLE(mydata->w));
}

/* .................................................................. */

void make_menu_item (gchar *label, GtkWidget *menu,
		     GtkSignalFunc func, gpointer data)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(label);
    gtk_menu_append(GTK_MENU(menu), item);
    gtk_signal_connect_object(GTK_OBJECT(item), "activate",
			      GTK_SIGNAL_FUNC(func), data);
    gtk_widget_show(item);
}

/* .................................................................. */

void allocate_fileptrs (void)
{
    int i;
    
    for (i=0; i<MAXRECENT; i++) {
	datap[i] = datalist[i];
	sessionp[i] = sessionlist[i];
	scriptp[i] = scriptlist[i];
    }
}

/* .................................................................. */

void mkfilelist (int filetype, const char *fname)
{
    char *tmp[MAXRECENT-1];
    char **filep;
    int i, match = -1;

    if (filetype == 1) filep = datap;
    else if (filetype == 2) filep = sessionp;
    else if (filetype == 3) filep = scriptp;
    else return;

    /* see if this file is already on the list */
    for (i=0; i<MAXRECENT; i++) {
        if (strcmp(filep[i], fname) == 0) {
            match = i;
            break;
        }
    }
    if (match == 0) return; /* no change in list */
    else  { /* clear menu files list before rebuilding */
	GtkWidget *w;
	char tmpname[MAXSTR];
	gchar itempath[80];
	gchar *pathstart[] = {"/File/Open data", 
			      "/Session/Open",
			      "/File/Open command file"};

	for (i=0; i<MAXRECENT; i++) {
	    sprintf(itempath, "%s/%d. %s", pathstart[filetype - 1],
		    i+1, endbit(tmpname, filep[i], -1));
	    w = gtk_item_factory_get_widget(mdata->ifac, itempath);
	    if (w != NULL) 
		gtk_item_factory_delete_item(mdata->ifac, itempath);
	}
    }
    
    /* save pointers to current order */
    for (i=0; i<MAXRECENT-1; i++) tmp[i] = filep[i];

    /* copy fname into array, if not already present */
    if (match == -1) {
        for (i=1; i<MAXRECENT; i++) {
            if (filep[i][0] == '\0') {
                strcpy(filep[i], fname);
                match = i;
                break;
	    }
	    if (match == -1) {
		match = MAXRECENT - 1;
		strcpy(filep[match], fname);
	    }
	}
    } 

    /* set first pointer to new file */
    filep[0] = filep[match];

    /* rearrange other pointers */
    for (i=1; i<=match; i++) filep[i] = tmp[i-1];

    add_files_to_menu(filetype);
}

/* .................................................................. */

char *endbit (char *dest, char *src, int addscore)
{
    /* take last part of src filename */
    if (strrchr(src, SLASH))
	strcpy(dest, strrchr(src, SLASH) + 1);
    else
	strcpy(dest, src);

    if (addscore != 0) {
	/* then either double (1) or delete (-1) any underscores */
	char mod[MAXSTR];
	size_t i, j, n;

	n = strlen(dest);
	j = 0;
	for (i=0; i<=n; i++) {
	    if (dest[i] != '_')
		mod[j++] = dest[i];
	    else {
		if (addscore == 1) {
		    mod[j++] = '_';
		    mod[j++] = dest[i];
		} 
	    }
	}
	strcpy(dest, mod);
    }
    return dest;
}

/* .................................................................. */

#if defined(USE_GNOME)

static void gnome_printfilelist (int filetype)
{
    int i;
    char **filep;
    char gpath[MAXLEN];
    static char *section[] = {"recent data files",
			      "recent session files",
			      "recent script files"};

    switch (filetype) {
    case 1: filep = datap; break;
    case 2: filep = sessionp; break;
    case 3: filep = scriptp; break;
    default: return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) { 
	    sprintf(gpath, "/gretl/%s/%d", section[filetype - 1], i);
	    gnome_config_set_string(gpath, filep[i]);
	} else break;
    }
}

#elif defined(G_OS_WIN32)

static void win_printfilelist (int filetype)
{
    int i;
    char **filep;
    char rpath[MAXLEN];
    static char *section[] = {"recent data files",
			      "recent session files",
			      "recent script files"};

    switch (filetype) {
    case 1: filep = datap; break;
    case 2: filep = sessionp; break;
    case 3: filep = scriptp; break;
    default: return;
    }

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) { 
	    sprintf(rpath, "%s\\%d", section[filetype - 1], i);
	    write_reg_val(HKEY_CURRENT_USER, rpath, filep[i]);
	} else break;
    }
}

#else /* "plain" version follows */

static void printfilelist (int filetype, FILE *fp)
{
    int i;
    char **filep;

    if (filetype == 1) {
	fprintf(fp, "recent data files:\n");
	filep = datap;
    } else if (filetype == 2) {
	fprintf(fp, "recent session files:\n");
	filep = sessionp;
    } else if (filetype == 3) {
	fprintf(fp, "recent script files:\n");
	filep = scriptp;
    } else 
	return;

    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) 
	    fprintf(fp, "%s\n", filep[i]);
	else break;
    }
}

#endif 

/* .................................................................. */

void add_files_to_menu (int filetype)
{
    int i;
    char **filep, tmp[MAXSTR];
    void (*callfunc)();
    GtkItemFactoryEntry fileitem;
    GtkWidget *w;
    gchar *msep[] = {"/File/Open data/sep",
		     "/Session/sep",
		     "/File/Open command file/sep"};
    gchar *mpath[] = {"/File/_Open data",
		     "/Session",
		     "/File/Open command file"};

    fileitem.path = NULL;

    if (filetype == 1) {
	callfunc = set_data_from_filelist;
	filep = datap;
    } else if (filetype == 2) {
	callfunc = set_session_from_filelist;
	filep = sessionp;
    } else if (filetype == 3) {
	callfunc = set_script_from_filelist;
	filep = scriptp;
    }
    else
	return;

    /* See if there are any files to add */
    if (filep[0][0] == '\0') return;
    else {
	gchar *itemtype = "<Separator>";
	GtkWidget *w;

	/* is a separator already in place? */
	w = gtk_item_factory_get_widget(mdata->ifac, msep[filetype - 1]);
	if (w == NULL) {
	    fileitem.path = mymalloc(80);
	    strcpy(fileitem.path, mpath[filetype - 1]);
	    strcat(fileitem.path, "/sep");
	    fileitem.accelerator = NULL;
	    fileitem.callback = NULL;
	    fileitem.callback_action = 0;
	    fileitem.item_type = itemtype;
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	}
    }

    /* put the files under the menu separator */
    for (i=0; i<MAXRECENT; i++) {
	if (filep[i][0]) {
	    if (fileitem.path == NULL)
		fileitem.path = mymalloc(80);
	    fileitem.accelerator = NULL;
	    fileitem.callback_action = i; 
	    fileitem.item_type = NULL;
	    sprintf(fileitem.path, "%s/%d. %s", mpath[filetype - 1],
		    i+1, endbit(tmp, filep[i], 1));
	    fileitem.callback = callfunc; 
	    gtk_item_factory_create_item(mdata->ifac, &fileitem, NULL, 1);
	    w = gtk_item_factory_get_widget_by_action(mdata->ifac, i);
	    if (w != NULL)
		gtk_tooltips_set_tip(gretl_tips, w, filep[i], NULL);
	} else break;
    }
    free(fileitem.path);
}

/* .................................................................. */

#ifndef G_OS_WIN32
# include <dlfcn.h>
#endif

int open_plugin (const char *plugin, void **handle)
{
    char pluginpath[MAXLEN];

#ifdef G_OS_WIN32
    sprintf(pluginpath, "%s\\%s.dll", paths.gretldir, plugin);
    *handle = LoadLibrary(pluginpath);
    if (*handle == NULL) {
	char buf[MAXLEN];

	sprintf(buf, "Couldn't load plugin %s", pluginpath);
	errbox(buf);
	return 1;
    }
#else
    sprintf(pluginpath, "%splugins/%s.so", paths.gretldir, plugin);
    *handle = dlopen(pluginpath, RTLD_LAZY);
    if (*handle == NULL) {
	errbox(dlerror());
	fprintf(stderr, "Failed to load plugin: %s\n", pluginpath);
	return 1;
    } 
#endif 
    return 0;
}

/* .................................................................. */

void *get_plugin_function (const char *funcname, void *handle)
{
    void *funp;
    char *error;

#ifdef G_OS_WIN32
    funp = GetProcAddress(handle, funcname);
    if (funp == NULL)  {
	errbox("Couldn't access plugin function");
	return NULL;
    }
#else
    funp = dlsym(handle, funcname);
    if ((error = dlerror()) != NULL)  {
	errbox(error);
	return NULL;
    }
#endif   
    return funp;
}

/* .................................................................. */

void close_plugin (void *handle)
{
#ifdef G_OS_WIN32
    FreeLibrary(handle);
#else
    dlclose(handle);
#endif
}

/* .................................................................. */

void get_default_dir (char *s)
{
    char *test = NULL;

    if (usecwd) {
	test = getcwd(s, MAXLEN);
	if (test == NULL) 
	    strcpy(s, paths.userdir);
	else
	    strcat(s, SLASHSTR);
    }
    else
	strcpy(s, paths.userdir);    
}
