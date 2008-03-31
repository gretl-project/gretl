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
#include "var.h"
#include "johansen.h"
#include "varprint.h"
#include "modelspec.h"
#include "forecast.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "system.h"
#include "matrix_extra.h"
#include "bootstrap.h"

#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

#include "model_table.h"
#include "series_view.h"
#include "console.h"
#include "session.h"
#include "textbuf.h"
#include "textutil.h"
#include "cmdstack.h"
#include "filelists.h"
#include "menustate.h"
#include "dlgutils.h"
#include "ssheet.h"
#include "datafiles.h"
#include "gpt_control.h"

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mini.tsplot.xpm"
#include "../pixmaps/mini.boxplot.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.pin.xpm"
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 8)
# include "../pixmaps/info_24.xpm"
#endif
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)
# include "../pixmaps/edit_24.xpm"
# include "../pixmaps/mini.edit.xpm"
#endif

/* for gretl toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.model.xpm"
#include "../pixmaps/mini.func.xpm"

char *storelist = NULL;

#define CONTENT_IS_CHANGED(w) (w->active_var == 1)

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static void view_window_save (GtkWidget *widget, windata_t *vwin);
static gint query_save_text (GtkWidget *w, GdkEvent *event, windata_t *vwin);
static void auto_save_script (windata_t *vwin);
static void auto_save_plot (windata_t *vwin);
static void add_model_dataset_items (windata_t *vwin);
static void add_model_tex_items (windata_t *vwin);
static void add_vars_to_plot_menu (windata_t *vwin);
static void add_dummies_to_plot_menu (windata_t *vwin);
static void add_system_menu_items (windata_t *vwin, int vecm);
static void add_x12_output_menu_item (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);
static void buf_edit_save (GtkWidget *widget, windata_t *vwin);
static void model_copy_callback (gpointer p, guint u, GtkWidget *w);
static void panel_heteroskedasticity_menu (windata_t *vwin);

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    GP_ITEM,
    PLOT_ITEM,
    RUN_ITEM,
    COPY_ITEM,
    TEX_ITEM,
    ADD_DATA_ITEM,
    ADD_MATRIX_ITEM,
    MAIL_ITEM,
    HELP_ITEM,
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    INDEX_ITEM,
    EDIT_SCRIPT_ITEM,
    STICKIFY_ITEM
} viewbar_flags;

static GtkWidget *get_toolbar_button_by_flag (GtkToolbar *tb, int flag)
{
    GList *kids;
    GtkToolbarChild *child;
    GtkWidget *w = NULL;
    int wflag;

    if (tb == NULL) {
	return NULL;
    }

    kids = tb->children;

    while (kids != NULL) {
	child = kids->data;
	if (child->type == GTK_TOOLBAR_CHILD_BUTTON) {
	    wflag = GPOINTER_TO_INT(g_object_get_data
				    (G_OBJECT(child->widget), "flag"));
	    if (wflag == flag) {
		w = child->widget;
		break;
	    }
	}
	kids = kids->next;
    }

    return w;
}

static void mark_content_changed (windata_t *vwin) 
{
    if (vwin->active_var == 0) {
	GtkWidget *w = get_toolbar_button_by_flag(GTK_TOOLBAR(vwin->mbar), 
						  SAVE_ITEM);

	if (w != NULL) {
	    gtk_widget_set_sensitive(w, TRUE);
	}
	vwin->active_var = 1;
    }
}

void mark_content_saved (windata_t *vwin) 
{
    GtkWidget *w = get_toolbar_button_by_flag(GTK_TOOLBAR(vwin->mbar), 
					      SAVE_ITEM);

    if (w != NULL) {
	gtk_widget_set_sensitive(w, FALSE);
    }
    vwin->active_var = 0;

    w = get_toolbar_button_by_flag(GTK_TOOLBAR(vwin->mbar), 
				   SAVE_AS_ITEM);
    if (w != NULL) {
	gtk_widget_set_sensitive(w, TRUE);
    }
}

static void close_model (gpointer data, guint close, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;

    if (window_is_busy(vwin)) {
	maybe_raise_dialog();
    } else {
	gtk_widget_destroy(vwin->dialog);
    }
}

static int arma_by_x12a (const MODEL *pmod)
{
    int ret = 0;

    if (pmod->ci == ARMA) {
	int acode = gretl_model_get_int(pmod, "arma_flags");

	if (acode & ARMA_X12A) {
	    ret = 1;
	}
    }

    return ret;
}

static int latex_is_ok (void)
{
    static int latex_ok = -1; 
  
    if (latex_ok == -1) {
	latex_ok = check_for_prog(latex);
    }

    return latex_ok;
}

static void model_output_save_callback (gpointer p, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;

    copy_format_dialog(vwin, W_SAVE);    
}

static gretlopt tex_eqn_opt;

static void eqn_set_show_stderrs (gpointer p, guint u, GtkWidget *w)
{
    if (u) {
	tex_eqn_opt = OPT_NONE;
    } else {
	tex_eqn_opt = OPT_T;
    }
}

gretlopt get_tex_eqn_opt (void)
{
    return tex_eqn_opt;
}

static GtkItemFactoryEntry model_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as..."), NULL, model_output_save_callback, 0, 
      "<StockItem>", GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as _icon"), NULL, model_add_as_icon, 
      0, NULL, GNULL },
    { N_("/File/Save as icon and cl_ose"), NULL, model_add_as_icon_and_close, 
      0, NULL, GNULL },
#ifdef NATIVE_PRINTING
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
#endif
    { N_("/File/_Close"), NULL, close_model, 0, "<StockItem>", GTK_STOCK_CLOSE },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 1, "<StockItem>", GTK_STOCK_COPY },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/Tests/_Omit variables"), NULL, selector_callback, OMIT, NULL, GNULL },
    { N_("/Tests/_Add variables"), NULL, selector_callback, ADD, NULL, GNULL },
    { N_("/Tests/_Sum of coefficients"), NULL, selector_callback, COEFFSUM, NULL, GNULL },
    { N_("/Tests/_Linear restrictions"), NULL, gretl_callback, RESTRICT, NULL, GNULL },
    { "/Tests/sep1", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/Non-linearity (s_quares)"), NULL, do_lmtest, LMTEST_SQUARES, NULL, GNULL },
    { N_("/Tests/Non-linearity (_logs)"), NULL, do_lmtest, LMTEST_LOGS, NULL, GNULL },
    { N_("/Tests/_Ramsey's RESET"), NULL, do_reset, RESET, NULL, GNULL },
    { "/Tests/sep2", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/_Heteroskedasticity"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Tests/Heteroskedasticity/White's test"), NULL, do_lmtest, LMTEST_WHITE, NULL, GNULL },
    { N_("/Tests/Heteroskedasticity/Breusch-Pagan"), NULL, do_lmtest, LMTEST_BP, NULL, GNULL },
    { N_("/Tests/Heteroskedasticity/Koenker"), NULL, do_lmtest, LMTEST_BPK, NULL, GNULL },
    { N_("/Tests/_Normality of residual"), NULL, do_resid_freq, TESTUHAT, NULL, GNULL },
    { N_("/Tests/_Influential observations"), NULL, do_leverage, LEVERAGE, NULL, GNULL },
    { N_("/Tests/_Collinearity"), NULL, do_vif, VIF, NULL, GNULL },
    { "/Tests/sep3", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/_Autocorrelation"), NULL, do_autocorr, LMTEST, NULL, GNULL },
    { N_("/Tests/A_RCH"), NULL, do_arch, ARCH, NULL, GNULL },
    { N_("/Tests/_Chow test"), NULL, do_chow_cusum, CHOW, NULL, GNULL },
    { N_("/Tests/_QLR test"), NULL, do_chow_cusum, QLRTEST, NULL, GNULL },
    { N_("/Tests/_CUSUM test"), NULL, do_chow_cusum, CUSUM, NULL, GNULL },
    { N_("/Tests/CUSUM_SQ test"), NULL, do_chow_cusum, CUSUMSQ, NULL, GNULL },
    { "/Tests/sep4", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/_Panel diagnostics"), NULL, do_panel_diagnostics, HAUSMAN, NULL, GNULL },
    { N_("/_Save"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/_Graphs"), NULL, NULL, 0, "<Branch>", GNULL }, 
    { N_("/Graphs/_Residual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Graphs/_Fitted, actual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/_Analysis"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Analysis/_Display actual, fitted, residual"), NULL, 
      display_fit_resid, 0, NULL, GNULL },
    { N_("/Analysis/_Forecasts..."), NULL, 
      gui_do_forecast, 0, NULL, GNULL },
    { N_("/Analysis/_Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL, GNULL },
    { N_("/Analysis/Confidence _ellipse..."), NULL, 
      selector_callback, ELLIPSE, NULL, GNULL },
    { N_("/Analysis/Coefficient covariance _matrix"), NULL, 
      do_outcovmx, 0, NULL, GNULL },
    { N_("/Analysis/ANOVA"), NULL, do_anova, 0, NULL, GNULL },
    { N_("/Analysis/Bootstrap..."), NULL, do_bootstrap, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL },
};

static GtkItemFactoryEntry model_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/_View"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/View/_Tabular"), NULL, model_tex_view, 
      GRETL_FORMAT_TEX, NULL, GNULL },
    { N_("/LaTeX/View/_Equation"), NULL, model_tex_view, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL, GNULL },
    { N_("/LaTeX/_Copy"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Copy/_Tabular"), NULL, window_copy, 
      GRETL_FORMAT_TEX, NULL, GNULL },
    { N_("/LaTeX/Copy/_Equation"), NULL, window_copy, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL, GNULL },
    { N_("/LaTeX/_Save"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Save/_Tabular"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX, NULL, GNULL },
    { N_("/LaTeX/Save/_Equation"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL, GNULL },
    { N_("/LaTeX/_Equation options"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/Equation options/Show _standard errors"), NULL, 
      eqn_set_show_stderrs, 1, "<RadioItem>", GNULL },
    { N_("/LaTeX/Equation options/Show _t-ratios"), NULL, 
      eqn_set_show_stderrs, 0, "/LaTeX/Equation options/Show standard errors", 
      GNULL },
    { N_("/LaTeX/_Tabular options..."), NULL, tex_format_dialog, 
      0, NULL, GNULL }
};

static GtkItemFactoryEntry sys_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/_View"), NULL, system_tex_callback, 0, NULL, GNULL },
    { N_("/LaTeX/_Copy"), NULL, system_tex_callback, 1, NULL, GNULL },
    { N_("/LaTeX/_Save"), NULL, system_tex_callback, 2, NULL, GNULL }
};

static GtkItemFactoryEntry system_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as..."), NULL, model_output_save_callback, 0, "<StockItem>", 
      GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as _icon"), NULL, model_add_as_icon, 
      0, NULL, GNULL },
    { N_("/File/Save as icon and cl_ose"), NULL, model_add_as_icon_and_close, 
      0, NULL, GNULL },
#ifdef NATIVE_PRINTING
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
#endif
    { N_("/File/_Close"), NULL, close_model, 0, "<StockItem>", GTK_STOCK_CLOSE },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, "<StockItem>", GTK_STOCK_COPY },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

static void model_copy_callback (gpointer p, guint u, GtkWidget *w)
{
    copy_format_dialog(p, W_COPY);
}

#ifdef ENABLE_NLS
gchar *menu_translate (const gchar *path, gpointer p)
{
    return (_(path));
}
#endif

int copyfile (const char *src, const char *dest) 
{
    FILE *srcfd, *destfd;
    char buf[GRETL_BUFSIZE];
    size_t n;

    if (!strcmp(src, dest)) return 1;
   
    if ((srcfd = gretl_fopen(src, "rb")) == NULL) {
	file_read_errbox(src);
	return 1; 
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
	file_write_errbox(dest);
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

int isdir (const char *path)
{
    struct stat buf;

    return (stat(path, &buf) == 0 && S_ISDIR(buf.st_mode)); 
}

FILE *gretl_tempfile_open (char *fname)
{
    FILE *fp = NULL;
    int fd;

    strcat(fname, ".XXXXXX");
#ifdef G_OS_WIN32
    fd = g_mkstemp(fname);
#else
    fd = mkstemp(fname);
#endif
    if (fd != -1) {
	fp = fdopen(fd, "w+");
	if (fp == NULL) {
	    file_write_errbox(fname);
	    close(fd);
	    remove(fname);
	}
    }

    return fp;
}

int gretl_tempname (char *fname)
{
    int fd, err = 0;

    strcat(fname, ".XXXXXX");
#ifdef G_OS_WIN32
    fd = g_mkstemp(fname);
#else
    fd = mkstemp(fname);
#endif
    if (fd == -1) {
	file_write_errbox(fname);
	err = 1;
    } else {
	close(fd);
	remove(fname);
    }

    return err;
}

int probably_native_datafile (const char *fname)
{
    char test[5];
    int len = strlen(fname);
    int ret = 0;

    if (len >= 5) {
	*test = 0;
	strncat(test, fname + len - 4, 4);
	lower(test);
	if (!strcmp(test, ".gdt")) {
	    ret = 1;
	} 
    }

    return ret;
}

int probably_script_file (const char *fname)
{
    char test[5];
    int len = strlen(fname);

    if (len < 5) return 0;

    *test = 0;
    strncat(test, fname + len - 4, 4);
    lower(test);

    return !strcmp(test, ".inp");
}

int probably_session_file (const char *fname)
{
    char test[7];
    int len = strlen(fname);

    if (len < 7) return 0;

    *test = 0;
    strncat(test, fname + len - 6, 6);
    lower(test);

    return !strcmp(test, ".gretl");
}

static int max_var_in_stacked_models (GtkWidget **wstack, int nwin)
{
    int i, role, mvm, vmax = 0;

    for (i=0; i<nwin; i++) {
	if (wstack[i] != NULL) {
	    role = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(wstack[i]), "role"));
	    if (role == VIEW_MODEL) {
		const MODEL *pmod;

		pmod = g_object_get_data(G_OBJECT(wstack[i]), "object");
		if (pmod != NULL) {
		    mvm = highest_numbered_var_in_model(pmod, datainfo);
		    if (mvm > vmax) {
			vmax = mvm;
		    }
		}
	    } else if (role == VAR || role == VECM) {
		const GRETL_VAR *var;

		var = g_object_get_data(G_OBJECT(wstack[i]), "object");
		if (var != NULL) {
		    mvm = gretl_VAR_get_highest_variable(var);
		    if (mvm > vmax) {
			vmax = mvm;
		    }		    
		}
	    } 
	}
    }    

    return vmax;
}

/* Below: Keep a record of (most) windows that are open, so they can
   be destroyed en masse when a new data file is opened, to prevent
   weirdness that could arise if (e.g.) a model window that pertains
   to a previously opened data file remains open after the data set
   has been changed.  Script windows are exempt, otherwise they are
   likely to disappear when their "run" control is activated, which we
   don't want.
*/

enum winstack_codes {
    STACK_INIT,
    STACK_ADD,
    STACK_REMOVE,
    STACK_DESTROY,
    STACK_QUERY,
    STACK_MATCH_FNAME,
    STACK_MAXVAR
};

static int 
winstack (int code, GtkWidget *w, gconstpointer ptest, GtkWidget **pw)
{
    static int n_windows;
    static GtkWidget **wstack;
    int i, ret = 0;

    switch (code) {

    case STACK_DESTROY:	
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] != NULL) {
		fprintf(stderr, "winstack: destroying widget at %p\n", 
			(void *) wstack[i]);
		gtk_widget_destroy(wstack[i]);
	    }
	}
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
	    GtkWidget **newstack;

	    newstack = myrealloc(wstack, (n_windows + 1) * sizeof *wstack);
	    if (newstack != NULL) { 
		wstack = newstack;
		wstack[n_windows] = w;
		n_windows++;
	    }
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
		    if (pw != NULL) {
			*pw = wstack[i];
		    }
		    ret = 1;
		    break;
		}
	    }
	}
	break;

    case STACK_MATCH_FNAME:
	if (wstack != NULL) {
	    const char *ctest = (const char *) ptest;

	    for (i=0; i<n_windows; i++) {
		if (wstack[i] != NULL) {
		    windata_t *vwin = 
			g_object_get_data(G_OBJECT(wstack[i]), "object");

		    if (vwin != NULL && strstr(ctest, vwin->fname)) {
			if (pw != NULL) {
			    *pw = wstack[i];
			}
			ret = 1;
			break;
		    }
		}
	    }
	}
	break;

    case STACK_MAXVAR:
	ret = max_var_in_stacked_models(wstack, n_windows);
	break;	

    default:
	break;
    }

    return ret;
}

void winstack_init (void)
{
    winstack(STACK_INIT, NULL, NULL, NULL);
}
    
void winstack_destroy (void)
{
    winstack(STACK_DESTROY, NULL, NULL, NULL);
}

int winstack_match_data (const gpointer p)
{
    return winstack(STACK_QUERY, NULL, p, NULL);
}

GtkWidget *match_window_by_data (const gpointer p)
{
    GtkWidget *w = NULL;

    winstack(STACK_QUERY, NULL, p, &w);
    return w;
}

GtkWidget *match_window_by_filename (const char *fname)
{
    GtkWidget *w = NULL;

    winstack(STACK_MATCH_FNAME, NULL, fname, &w);
    return w;
}

int highest_numbered_variable_in_winstack (void)
{
    return winstack(STACK_MAXVAR, NULL, NULL, NULL);
}

void winstack_add (GtkWidget *w)
{
    winstack(STACK_ADD, w, NULL, NULL);
}

void winstack_remove (GtkWidget *w)
{
    winstack(STACK_REMOVE, w, NULL, NULL);
}

/* ........................................................... */

static void delete_file (GtkWidget *widget, char *fname) 
{
    remove(fname);
    g_free(fname);
}

static void delete_file_viewer (GtkWidget *widget, windata_t *vwin) 
{
    gint resp = 0;

    if (window_is_busy(vwin)) {
	maybe_raise_dialog();
	return;
    }

    if ((vwin->role == EDIT_SCRIPT || vwin->role == EDIT_HEADER ||
	 vwin->role == EDIT_NOTES || vwin->role == GR_PLOT) && 
	CONTENT_IS_CHANGED(vwin)) {
	resp = query_save_text(NULL, NULL, vwin);
    }

    if (!resp) {
	gtk_widget_destroy(vwin->dialog); 
    }
}

void delete_widget (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
}

#if defined(HAVE_FLITE) || defined(G_OS_WIN32)

static int set_or_get_audio_stop (int set, int val)
{
    static int audio_quit;

    if (set) audio_quit = val;

    return audio_quit;
}

static int should_stop_talking (void)
{
    while (gtk_events_pending()) {
	gtk_main_iteration();
    }

    return set_or_get_audio_stop(0, 0);
}

void stop_talking (void)
{
    set_or_get_audio_stop(1, 1);
}

void audio_render_window (windata_t *vwin, int key)
{
    void *handle;
    int (*read_window_text) (windata_t *, const DATAINFO *, int (*)());

    if (vwin == NULL) {
	stop_talking();
	return;
    }

    read_window_text = gui_get_plugin_function("read_window_text", 
					       &handle);
    if (read_window_text == NULL) {
        return;
    }

    set_or_get_audio_stop(1, 0);

    if (key == AUDIO_LISTBOX) {
	(*read_window_text) (vwin, NULL, &should_stop_talking);
    } else {
	(*read_window_text) (vwin, datainfo, &should_stop_talking);
    }

    close_plugin(handle);
}

#endif

static gboolean Ctrl_C (windata_t *vwin)
{
#ifdef G_OS_WIN32 
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    if (gtk_text_buffer_get_selection_bounds(buf, NULL, NULL)) {
	window_copy(vwin, GRETL_FORMAT_SELECTION, NULL);
    } else if (MULTI_FORMAT_ENABLED(vwin->role)) {
	window_copy(vwin, GRETL_FORMAT_RTF, NULL);
    } else {
	window_copy(vwin, GRETL_FORMAT_TXT, NULL);
    }
    return TRUE;
#else
    return FALSE;
#endif
}

static gint catch_viewer_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods;
    int editing;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);
    editing = gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->w));

    if (mods & GDK_CONTROL_MASK) {
	if (key->keyval == GDK_f) {
	    /* Ctrl-S: find */
	    text_find_callback(NULL, vwin);
	    return TRUE;
	} else if (key->keyval == GDK_c) {
	    /* Ctrl-C: copy */
	    return Ctrl_C(vwin);
	} else if (editing) {
	    if (gdk_keyval_to_upper(key->keyval) == GDK_S) { 
		/* Ctrl-S: save */
		if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		    buf_edit_save(NULL, vwin);
		} else {
		    view_window_save(NULL, vwin);
		}
	    } else if (gdk_keyval_to_upper(key->keyval) == GDK_Q) {
		/* Ctrl-Q: quit */
		if (vwin->role == EDIT_SCRIPT && CONTENT_IS_CHANGED(vwin)) {
		    if (query_save_text(NULL, NULL, vwin) == FALSE) {
			gtk_widget_destroy(vwin->dialog);
		    }
		} else { 
		    gtk_widget_destroy(w);
		}
	    } 
	} 
    }

    if (editing) {
	/* we respond to plain keystrokes below: this won't do if we're
	   editing text */
	return FALSE;
    }

    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    } else if (key->keyval == GDK_s && Z != NULL && vwin->role == VIEW_MODEL) {
	model_add_as_icon_and_close(vwin, GRETL_OBJ_EQN, NULL);
    } else if (key->keyval == GDK_w) {
	if (mods & GDK_CONTROL_MASK) {
	    gtk_widget_destroy(w);
	    return TRUE;
	}	
    }
#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    else if (key->keyval == GDK_a) {
	audio_render_window(vwin, AUDIO_TEXT);
    } else if (key->keyval == GDK_x) {
	stop_talking();
    }
#endif

    return FALSE;
}

/* ........................................................... */

void nomem (void)
{
    errbox(_("Out of memory!"));
}

void *mymalloc (size_t size) 
{
    void *mem;
   
    if ((mem = malloc(size)) == NULL) {
	nomem();
    }

    return mem;
}

void *myrealloc (void *ptr, size_t size) 
{
    void *mem;

    if ((mem = realloc(ptr, size)) == NULL) {
	nomem();
    }

    return mem;
}

/* ........................................................... */

void mark_dataset_as_modified (void)
{
    data_status |= MODIFIED_DATA;
    set_sample_label(datainfo);
}

void register_data (char *fname, const char *user_fname,
		    int record)
{    
    /* basic accounting */
    data_status |= HAVE_DATA;
    orig_vars = datainfo->v;

    /* set appropriate data_status bits */
    if (fname == NULL) {
	data_status |= GUI_DATA;
	mark_dataset_as_modified();
    } else if (!(data_status & IMPORT_DATA)) {
	if (has_system_prefix(paths.datfile, &paths, DATA_SEARCH)) {
	    data_status |= BOOK_DATA;
	    data_status &= ~USER_DATA;
	} else {
	    data_status &= ~BOOK_DATA;
	    data_status |= USER_DATA; 
	}
	if (is_gzipped(paths.datfile)) {
	    data_status |= GZIPPED_DATA;
	} else {
	    data_status &= ~GZIPPED_DATA;
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
	gretl_command_sprintf("open %s", user_fname ? user_fname : fname);
	check_and_record_command();
    } 

    if (fname != NULL) {
	/* focus the data window */
	gtk_widget_grab_focus(mdata->listbox);
    }

    /* invalidate "remove extra obs" menu item */
    drop_obs_state(FALSE);
}

#define APPENDING(action) (action == APPEND_DATA || \
                           action == APPEND_CSV || \
                           action == APPEND_GNUMERIC || \
                           action == APPEND_EXCEL || \
                           action == APPEND_ODS || \
                           action == APPEND_ASCII || \
                           action == APPEND_WF1 || \
                           action == APPEND_DTA || \
                           action == APPEND_JMULTI)

int get_worksheet_data (char *fname, int datatype, int append)
{
    void *handle;
    PRN *errprn;
    const char *errbuf;
    FILE *fp;
    int (*sheet_get_data)(const char*, double ***, DATAINFO *, 
			  gretlopt, PRN *);
    int err = 0;
    
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	delete_from_filelist(FILE_LIST_DATA, fname);
	file_read_errbox(fname);
	return 1;
    } else {
	fclose(fp);
    }

    if (datatype == GRETL_GNUMERIC) {
	sheet_get_data = gui_get_plugin_function("gnumeric_get_data",
						 &handle);
    } else if (datatype == GRETL_EXCEL) {
	sheet_get_data = gui_get_plugin_function("xls_get_data",
						 &handle);
    } else if (datatype == GRETL_ODS) {
	sheet_get_data = gui_get_plugin_function("ods_get_data",
						 &handle);
    } else if (datatype == GRETL_WF1) {
	sheet_get_data = gui_get_plugin_function("wf1_get_data",
						 &handle);
    } else if (datatype == GRETL_DTA) {
	sheet_get_data = gui_get_plugin_function("dta_get_data",
						 &handle);
    } else if (datatype == GRETL_JMULTI) {
	sheet_get_data = gui_get_plugin_function("jmulti_get_data",
						 &handle);
    } else {
	errbox(_("Unrecognized data type"));
	return 1;
    }

    if (sheet_get_data == NULL) {
        return 1;
    }

    if (bufopen(&errprn)) {
	close_plugin(handle);
	return 1;
    }

    err = (*sheet_get_data)(fname, &Z, datainfo, OPT_G, errprn);
    close_plugin(handle);

    if (err == -1) {
	fprintf(stderr, "data import canceled\n");
	gretl_print_destroy(errprn);
	return E_CANCEL;
    }

    errbuf = gretl_print_get_buffer(errprn);

    if (err) {
	if (errbuf != NULL && *errbuf != '\0') {
	    errbox(errbuf);
	} else {
	    errbox(_("Failed to import data"));
	}
	gretl_print_destroy(errprn);
	delete_from_filelist(FILE_LIST_DATA, fname);
	return 1;
    } else {
	if (errbuf != NULL && *errbuf != '\0') {
	    infobox(errbuf);
	}
	if (datatype == GRETL_DTA) {
	    maybe_display_string_table();
	}
    }

    gretl_print_destroy(errprn);

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	data_status |= IMPORT_DATA;
	if (fname != paths.datfile) {
	    strcpy(paths.datfile, fname);
	}
	if (mdata != NULL) {
	    register_data(fname, NULL, 1);
	}
	if (!dataset_is_time_series(datainfo) && 
	    !dataset_is_panel(datainfo) && mdata != NULL) {
	    int resp;

	    resp = yes_no_dialog(_("gretl: open data"),
				 _("The imported data have been interpreted as undated\n"
				   "(cross-sectional).  Do you want to give the data a\n"
				   "time-series or panel interpretation?"),
				 0);
	    if (resp == GRETL_YES) {
		data_structure_wizard(NULL, 0, NULL);
	    }
	}
    }

    return err;
}

/* cases for do_open_data: 
   - called from dialog: user has said Yes to opening data file,
     although a data file is already open (or user wants to append
     data)
   - reached without dialog, when no datafile is open yet
*/

void do_open_data (GtkWidget *w, gpointer data, int code)
{
    gint datatype, err = 0;
    dialog_t *dlg = NULL;
    windata_t *fwin = NULL;
    int append = APPENDING(code);

    if (data != NULL) {    
	if (w == NULL) { /* not coming from edit_dialog */
	    fwin = (windata_t *) data;
	} else {
	    dlg = (dialog_t *) data;
	    fwin = (windata_t *) edit_dialog_get_data(dlg);
	}
    }

    if (code == OPEN_CSV || code == APPEND_CSV || code == OPEN_ASCII ||
	code == APPEND_ASCII) {
	datatype = GRETL_CSV_DATA;
    } else if (code == OPEN_GNUMERIC || code == APPEND_GNUMERIC) {
	datatype = GRETL_GNUMERIC;
    } else if (code == OPEN_ODS || code == APPEND_ODS) {
	datatype = GRETL_ODS;
    } else if (code == OPEN_EXCEL || code == APPEND_EXCEL) {
	datatype = GRETL_EXCEL;
    } else if (code == OPEN_OCTAVE || code == APPEND_OCTAVE) {
	datatype = GRETL_OCTAVE;
    } else if (code == OPEN_WF1 || code == APPEND_WF1) {
	datatype = GRETL_WF1;
    } else if (code == OPEN_DTA || code == APPEND_DTA) {
	datatype = GRETL_DTA;
    } else if (code == OPEN_JMULTI || code == APPEND_JMULTI) {
	datatype = GRETL_JMULTI;
    } else {
	/* no filetype specified: have to guess */
	PRN *prn;	

	if (bufopen(&prn)) return;
	datatype = detect_filetype(tryfile, &paths, prn);
	gretl_print_destroy(prn);
    }

    /* destroy the current data set, etc., unless we're explicitly appending */
    if (!append) {
	close_session(NULL, &Z, datainfo, OPT_NONE); /* FIXME opt */
    }

    if (datatype == GRETL_GNUMERIC || datatype == GRETL_EXCEL ||
	datatype == GRETL_WF1 || datatype == GRETL_DTA ||
	datatype == GRETL_JMULTI || datatype == GRETL_ODS) {
	get_worksheet_data(tryfile, datatype, append);
	return;
    } else if (datatype == GRETL_CSV_DATA) {
	do_open_csv_octave(tryfile, OPEN_CSV, append);
	return;
    } else if (datatype == GRETL_OCTAVE) {
	do_open_csv_octave(tryfile, OPEN_OCTAVE, append);
	return;
    } else { 
	/* native data */
	PRN *errprn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	if (datatype == GRETL_XML_DATA) {
	    err = gretl_read_gdt(&Z, datainfo, tryfile, &paths, 
				 OPT_P, errprn);
	} else {
	    err = gretl_get_data(&Z, datainfo, tryfile, &paths, 
				 OPT_NONE, errprn);
	}

	gretl_print_destroy(errprn);
    }

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_DATA, tryfile);
	return;
    }	

    /* trash the practice files window that launched the query? */
    if (fwin != NULL) {
	gtk_widget_destroy(fwin->w);
    }

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	strcpy(paths.datfile, tryfile);
	register_data(paths.datfile, NULL, 1);
    } 
}

/* give user choice of not opening selected datafile, if there's
   already a datafile open */

void verify_open_data (gpointer userdata, int code)
{
    if (dataset_locked()) {
	return;
    }

    if (data_status) {
	int resp = 
	    yes_no_dialog (_("gretl: open data"), 
			   _("Opening a new data file will automatically\n"
			     "close the current one.  Any unsaved work\n"
			     "will be lost.  Proceed to open data file?"), 0);

	if (resp != GRETL_YES) return;
    } 

    do_open_data(NULL, userdata, code);
}

/* give user choice of not opening session file, if there's already a
   datafile open */

void verify_open_session (void)
{
    if (!gretl_is_pkzip_file(tryfile)) {
	/* not a new-style zipped session file */
	do_open_script();
	return;
    }

    if (data_status) {
	int resp = 
	    yes_no_dialog (_("gretl: open session"), 
			   _("Opening a new session file will automatically\n"
			     "close the current session.  Any unsaved work\n"
			     "will be lost.  Proceed to open session file?"), 0);

	if (resp != GRETL_YES) return;
    }

    do_open_session();
}

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor(vwin->w, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
}

static void buf_edit_save (GtkWidget *widget, windata_t *vwin)
{
    char **pbuf = (char **) vwin->data;
    gchar *text;

    text = textview_get_text(vwin->w);

    if (text == NULL || *text == '\0') {
	errbox(_("Buffer is empty"));
	g_free(text);
	return;
    }

    /* swap the edited text into the buffer */
    free(*pbuf); 
    *pbuf = text;

    if (vwin->role == EDIT_HEADER) {
	mark_content_saved(vwin);
	mark_dataset_as_modified();
    } else if (vwin->role == EDIT_NOTES) {
	mark_content_saved(vwin);
	mark_session_changed();
    }
}

static int update_func_code (windata_t *vwin)
{
    int iface, err = 0;

    /* callback used when editing a function in the context of
       the "function package editor" */
	
    iface = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->w), "iface"));
    err = update_function_from_script(vwin->fname, iface);
    if (err) {
	gui_errmsg(err);
    }

    return err;
}

static void view_window_save (GtkWidget *widget, windata_t *vwin)
{
    if (strstr(vwin->fname, "script_tmp") || *vwin->fname == '\0') {
	/* special case: a newly created script */
	file_save(vwin, SAVE_SCRIPT, NULL);
	strcpy(vwin->fname, scriptfile);
	mark_content_saved(vwin);
    } else {
	FILE *fp;
	gchar *text;

	if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	    errbox(_("Can't open file for writing"));
	    return;
	} else {
	    text = textview_get_text(vwin->w);
	    system_print_buf(text, fp);
	    fclose(fp);
	    g_free(text);
	    mark_content_saved(vwin);
	    if (vwin->role == EDIT_FUNC_CODE) {
		update_func_code(vwin);
	    }
	}
    }
}

void windata_init (windata_t *vwin)
{
    vwin->dialog = NULL;
    vwin->vbox = NULL;
    vwin->listbox = NULL;
    vwin->mbar = NULL;
    vwin->w = NULL;
    vwin->status = NULL;
    vwin->popup = NULL;
    vwin->ifac = NULL;
    vwin->gretl_parent = NULL;
    vwin->gretl_children = NULL;
    vwin->data = NULL;
    vwin->active_var = 0;
    vwin->role = 0;
    vwin->n_model_tests = 0;
    vwin->n_gretl_children = 0;
    vwin->flags = 0;
    vwin->fname[0] = '\0';
    vwin->sbuf = NULL;
}

static int vwin_add_child (windata_t *parent, windata_t *child)
{
    int n = parent->n_gretl_children;
    int i, done = 0, err = 0;

    for (i=0; i<n; i++) {
	if (parent->gretl_children[i] == NULL) {
	    /* reuse a vacant slot */
	    parent->gretl_children[i] = child;
	    done = 1;
	    break;
	}
    }

    if (!done) {
	windata_t **children;

	children = realloc(parent->gretl_children, (n + 1) * sizeof *children);
	if (children == NULL) {
	    err = 1;
	} else {
	    parent->gretl_children = children;
	    parent->gretl_children[n] = child;
	    parent->n_gretl_children += 1;
	}
    }
    
    if (!err) {
	child->gretl_parent = parent;
    }

    return err;
}

static void vwin_nullify_child (windata_t *parent, windata_t *child)
{
    int i, n = parent->n_gretl_children;

    for (i=0; i<n; i++) {
	if (child == parent->gretl_children[i]) {
	    parent->gretl_children[i] = NULL;
	}
    }
}

windata_t *vwin_first_child (windata_t *vwin)
{
    int i, n = vwin->n_gretl_children;

    for (i=0; i<n; i++) {
	if (vwin->gretl_children[i] != NULL) {
	    return vwin->gretl_children[i];
	}
    }

    return NULL;
}

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin != NULL) {
	if (vwin->w != NULL) { 
	    gchar *undo = g_object_steal_data(G_OBJECT(vwin->w), "undo");
	    
	    if (undo != NULL) {
		g_free(undo);
	    }
	}

	/* notify parent, if any */
	if (vwin->gretl_parent != NULL) {
	    vwin_nullify_child(vwin->gretl_parent, vwin);
	}

	/* notify children, if any */
	if (vwin->n_gretl_children > 0) {
	    int i;

	    for (i=0; i<vwin->n_gretl_children; i++) {
		if (vwin->gretl_children[i] != NULL) {
		    vwin->gretl_children[i]->gretl_parent = NULL;
		}
	    }
	    free(vwin->gretl_children);
	}

	/* menu stuff */
	if (vwin->popup != NULL) {
	    gtk_widget_destroy(GTK_WIDGET(vwin->popup));
	}
	if (vwin->ifac != NULL) {
	    g_object_unref(G_OBJECT(vwin->ifac));
	}

	/* data specific to certain windows */
	if (vwin->role == SUMMARY) {
	    free_summary(vwin->data); 
	} else if (vwin->role == CORR || vwin->role == PCA || 
		   vwin->role == COVAR) {
	    free_vmatrix(vwin->data);
	} else if (vwin->role == FCAST || vwin->role == AFR) {
	    free_fit_resid(vwin->data);
	} else if (vwin->role == COEFFINT) {
	    free_coeff_intervals(vwin->data);
	} else if (vwin->role == VIEW_SERIES) {
	    free_series_view(vwin->data);
	} else if (vwin->role == VIEW_MODEL) {
	    gretl_object_unref(vwin->data, GRETL_OBJ_EQN);
	} else if (vwin->role == VAR || vwin->role == VECM) { 
	    gretl_object_unref(vwin->data, GRETL_OBJ_VAR);
	} else if (vwin->role == LEVERAGE) {
	    gretl_matrix_free(vwin->data);
	} else if (vwin->role == MAHAL) {
	    free_mahal_dist(vwin->data);
	} else if (vwin->role == XTAB) {
	    free_xtab(vwin->data);
	} else if (vwin->role == COINT2) {
	    gretl_VAR_free(vwin->data);
	} else if (vwin->role == SYSTEM) {
	    gretl_object_unref(vwin->data, GRETL_OBJ_SYS);
	} else if (vwin->role == PRINT && vwin->data != NULL) {
	    free_multi_series_view(vwin->data);
	} else if (vwin->role == GUI_HELP || vwin->role == GUI_HELP_EN) {
	    free(vwin->data); /* help file text */
	}

	if (window_delete_filename(vwin)) {
	    if (vwin->gretl_parent == NULL) {
		windata_t *child = vwin_first_child(vwin);

		if (child == NULL) {
		    remove(vwin->fname);
		}
	    }
	} else if (vwin->role == EDIT_FUNC_CODE) {
	    remove(vwin->fname);
	}

	if (vwin->role == NATIVE_SERIES) {
	    winstack_remove(vwin->w);
	} else if (vwin->dialog != NULL) {
	    winstack_remove(vwin->dialog);
	}

	free(vwin);
    }
}

static GtkIconFactory *gretl_stock_ifac;

void gretl_stock_icons_init (void)
{
    char **xpms[] = {
#if NO_INFO_ICON
	info_24_xpm,
#endif
#if NO_EDIT_ICON
	edit_24_xpm,
	mini_edit_xpm,
#endif
	mini_tex_xpm,
	mail_16_xpm,
	mini_tsplot_xpm,
	mini_boxplot_xpm,
	mini_pdf_xpm,
	mini_manual_xpm,
	mini_calc_xpm,
	mini_sh_xpm,
	mini_session_xpm,
	mini_plot_xpm,
	mini_model_xpm,
	mini_func_xpm,
	mini_pin_xpm
    };
    const char *stocks[] = {
#if NO_INFO_ICON
	GRETL_STOCK_INFO,
#endif
#if NO_EDIT_ICON
	GRETL_STOCK_EDIT,
	GRETL_STOCK_SCRIPT,
#endif
	GRETL_STOCK_TEX,
	GRETL_STOCK_MAIL,
	GRETL_STOCK_TS,
	GRETL_STOCK_BOX,
	GRETL_STOCK_PDF,
	GRETL_STOCK_BOOK,
	GRETL_STOCK_CALC,
	GRETL_STOCK_CONSOLE,
	GRETL_STOCK_ICONS,
	GRETL_STOCK_SCATTER,
	GRETL_STOCK_MODEL,
	GRETL_STOCK_FUNC,
	GRETL_STOCK_PIN
    };
    int n = sizeof stocks / sizeof stocks[0];

    if (gretl_stock_ifac == NULL) {
	GtkIconSource *source;
	GtkIconSet *set;
	GdkPixbuf *pbuf;
	int i;

	gretl_stock_ifac = gtk_icon_factory_new();

	for (i=0; i<n; i++) {
	    set = gtk_icon_set_new();
	    source = gtk_icon_source_new();
	    gtk_icon_source_set_size(source, GTK_ICON_SIZE_MENU);
	    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) xpms[i]);
	    gtk_icon_source_set_pixbuf(source, pbuf);
	    g_object_unref(pbuf);
	    gtk_icon_set_add_source(set, source);
	    gtk_icon_source_free(source);
	    gtk_icon_factory_add(gretl_stock_ifac, stocks[i], set);
	    gtk_icon_set_unref(set);
	}

	gtk_icon_factory_add_default(gretl_stock_ifac);
    }
}

#ifdef NATIVE_PRINTING
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(vwin, 0, w);
}
#endif

static void save_plot_commands_callback (GtkWidget *w, windata_t *vwin)
{
    auto_save_plot(vwin);
}

/* is any text selected? */

static int vwin_selection_present (gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GtkTextBuffer *buf;
    int ret = 0;

    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    if (gtk_text_buffer_get_selection_bounds(buf, NULL, NULL)) {
	ret = 1;
    }

    return ret;
}

#define editor_role(r) (r == EDIT_SCRIPT || \
                        r == EDIT_HEADER || \
                        r == EDIT_NOTES || \
                        r == EDIT_FUNC_CODE || \
                        r == GR_PLOT)

#define script_role(r) (r == VIEW_SCRIPT || r == VIEW_LOG)

static void text_copy_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin_selection_present(vwin)) {
	window_copy(vwin, GRETL_FORMAT_SELECTION, w);
    } else if (vwin->role == VIEW_SCALAR) {
	scalar_to_clipboard(vwin);
    } else if (!script_role(vwin->role) && !editor_role(vwin->role)) {
	copy_format_dialog(vwin, W_COPY);
    } else {
	window_copy(vwin, GRETL_FORMAT_TXT, w);
    }
}

static void text_paste_callback (GtkWidget *w, windata_t *vwin)
{
    text_paste(vwin, 0, w);
}

static void text_replace_callback (GtkWidget *w, windata_t *vwin)
{
    text_replace(vwin, 0, w);
}

static void text_undo_callback (GtkWidget *w, windata_t *vwin)
{
    text_undo(vwin, 0, w);
}

static void add_data_callback (GtkWidget *w, windata_t *vwin)
{
    int oldv = datainfo->v;

    if (vwin->role == PCA) {
	add_pca_data(vwin);
    } else if (vwin->role == LEVERAGE) {
	add_leverage_data(vwin);
    } else if (vwin->role == MAHAL) {
	add_mahalanobis_data(vwin);
    } else if (vwin->role == FCAST) {
	add_fcast_data(vwin);
    }

    if (datainfo->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }	
}

static void matrix_savename (GtkWidget *w, dialog_t *dlg)
{
    char *newname = (char *) edit_dialog_get_data(dlg);
    const gchar *buf = edit_dialog_get_text(dlg);

    if (buf == NULL || validate_varname(buf)) return;

    *newname = 0;
    strncat(newname, buf, VNAMELEN - 1);

    close_dialog(dlg);
}

static void add_matrix_callback (GtkWidget *w, windata_t *vwin)
{
    gretl_matrix *m = NULL;
    char mname[VNAMELEN];
    int err, cancel = 0;

    if (vwin->role == XTAB) {
	m = xtab_to_matrix(vwin->data);
	if (m == NULL) {
	    nomem();
	} else {
	    edit_dialog(_("gretl: save matrix"), 
			_("Enter a name"),
			NULL, matrix_savename, mname, 
			0, VARCLICK_NONE, &cancel);
	    if (cancel) {
		gretl_matrix_free(m);
	    } else {
		err = add_or_replace_user_matrix(m, mname);
		if (err) {
		    gui_errmsg(err);
		} else {
		    infobox(_("Saved matrix as %s"), mname);
		}
	    }
	}
    }
}

static void mail_script_callback (GtkWidget *w, windata_t *vwin)
{
    if (viewer_char_count(vwin) == 0) {
	infobox(_("Nothing to send"));
	return;
    }

    if (query_save_text(NULL, NULL, vwin)) {
	return;
    }
    
    send_file(vwin->fname);
}

static void window_help (GtkWidget *w, windata_t *vwin)
{
    context_help(NULL, GINT_TO_POINTER(vwin->role));
}

static void multi_save_as_callback (GtkWidget *w, windata_t *vwin)
{
    copy_format_dialog(vwin, W_SAVE);
}

static void script_index (GtkWidget *w, windata_t *vwin)
{
    display_files(NULL, PS_FILES, NULL);
}

static void set_output_sticky (GtkWidget *w, windata_t *vwin)
{
    const char *opts[] = {
	N_("Running a script replaces text"),
	N_("Running a script adds to text")
    };
    int resp, deflt;

    deflt = (vwin->flags & VWIN_STICKY)? 1 : 0;

    resp = radio_dialog(NULL, NULL, opts, 2, deflt, 0);

    if (resp == 0) {
	vwin->flags &= ~VWIN_STICKY;
    } else if (resp == 1) {
	vwin->flags |= VWIN_STICKY;
    }
}

typedef void (*toolfunc) (GtkWidget *w, windata_t *vwin);

struct viewbar_item {
    const char *str;
    const gchar *icon;
    toolfunc func;
    int flag;
};

static struct viewbar_item viewbar_items[] = {
    { N_("Save"), GTK_STOCK_SAVE, view_window_save, SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, file_save_callback, SAVE_AS_ITEM },
    { N_("Send to gnuplot"), GTK_STOCK_EXECUTE, gp_send_callback, GP_ITEM },
#ifdef NATIVE_PRINTING
    { N_("Print..."), GTK_STOCK_PRINT, window_print_callback, 0 },
#endif
    { N_("Run"), GTK_STOCK_EXECUTE, do_run_script, RUN_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, text_copy_callback, COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, text_paste_callback, EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, text_find_callback, 0 },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, text_replace_callback, EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, text_undo_callback, EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, series_view_sort, SORT_ITEM },    
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, series_view_sort_by, SORT_BY_ITEM },
    { N_("Configure..."), GTK_STOCK_PREFERENCES, script_tabs_dialog, EDIT_SCRIPT_ITEM },
    { N_("Send To..."), GRETL_STOCK_MAIL, mail_script_callback, MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, script_index, INDEX_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, activate_script_help, RUN_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, window_tex_callback, TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, series_view_graph, PLOT_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, series_view_format_dialog, FORMAT_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, add_data_callback, ADD_DATA_ITEM },
    { N_("Add as matrix..."), GTK_STOCK_ADD, add_matrix_callback, ADD_MATRIX_ITEM },
    { N_("Stickiness..."), GRETL_STOCK_PIN, set_output_sticky, STICKIFY_ITEM },
    { N_("Help"), GTK_STOCK_HELP, window_help, HELP_ITEM },
    { N_("Close"), GTK_STOCK_CLOSE, delete_file_viewer, 0 },
    { NULL, NULL, NULL, 0 }
};

static void set_plot_icon (struct viewbar_item *vitem)
{
    if (dataset_is_time_series(datainfo)) {
	vitem->icon = GRETL_STOCK_TS;
    } else {
	vitem->icon = GRETL_STOCK_BOX;
    }
}

#define run_ok(r) (r == EDIT_SCRIPT || \
                   r == VIEW_SCRIPT || \
                   r == VIEW_LOG)

#define edit_ok(r) (r == EDIT_SCRIPT || \
                    r == EDIT_HEADER || \
                    r == EDIT_NOTES || \
                    r == EDIT_FUNC_CODE || \
	            r == GR_PLOT || \
                    r == GR_BOX || \
                    r == SCRIPT_OUT)

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_FUNC_CODE && \
		       r != VIEW_SCALAR)

#define help_ok(r) (r == LEVERAGE || \
		    r == COINT2 || \
		    r == HURST || \
		    r == RMPLOT || \
		    r == MAHAL)

#define sort_ok(r)    (r == VIEW_SERIES)
#define format_ok(r)  (r == VIEW_SERIES || r == VIEW_SCALAR)
#define plot_ok(r)    (r == VIEW_SERIES)

#define add_data_ok(r) (r == PCA || r == LEVERAGE || \
                        r == MAHAL || r == FCAST)

/* Screen out unwanted menu items depending on the context; also
   adjust the callbacks associated with some items based on
   context.
*/

static toolfunc item_get_callback (struct viewbar_item *item, windata_t *vwin, 
				   int latex_ok, int sortby_ok)
{
    toolfunc func = item->func;
    int f = item->flag;
    int r = vwin->role;

    if (!edit_ok(r) && f == EDIT_ITEM) {
	return NULL;
    } else if (!run_ok(r) && f == RUN_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && f == MAIL_ITEM) {
	return NULL;
    } else if (!help_ok(r) && f == HELP_ITEM) {
	return NULL;
    } else if (r != GR_PLOT && f == GP_ITEM) {
	return NULL;
    } else if (r == VIEW_SCALAR && f == 0) {
	return NULL;
    } else if ((!latex_ok || !MULTI_FORMAT_ENABLED(r)) && f == TEX_ITEM) {
	return NULL;
    } else if (!add_data_ok(r) && f == ADD_DATA_ITEM) {
	return NULL;
    } else if (r != XTAB && f == ADD_MATRIX_ITEM) {
	return NULL;
    } else if (!sort_ok(r) && f == SORT_ITEM) {
	return NULL;
    } else if (!sortby_ok && f == SORT_BY_ITEM) {
	return NULL;
    } else if (!plot_ok(r) && f == PLOT_ITEM) {
	return NULL;
    } else if (!format_ok(r) && f == FORMAT_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && f == EDIT_SCRIPT_ITEM) {
	return NULL;
    } else if (r != VIEW_SCRIPT && f == INDEX_ITEM) {
	return NULL;
    } else if (r != SCRIPT_OUT && f == STICKIFY_ITEM) {
	return NULL;
    } else if (f == SAVE_ITEM) { 
	if (!edit_ok(r) || r == SCRIPT_OUT) {
	    /* script output doesn't already have a filename */
	    return NULL;
	}
	if (r == EDIT_HEADER || r == EDIT_NOTES) {
	    func = buf_edit_save;
	} else if (r == GR_PLOT) {
	    func = save_plot_commands_callback;
	}
    } else if (f == SAVE_AS_ITEM) {
	if (!save_as_ok(r)) {
	    return NULL;
	} else if (MULTI_FORMAT_ENABLED(r) || (r == PRINT && vwin->data != NULL)) {
	    func = multi_save_as_callback;
	}
    }

    return func;
}

static void make_viewbar (windata_t *vwin, int text_out)
{
    GtkWidget *hbox, *button, *w;
    toolfunc func;
    int sortby_ok = has_sortable_data(vwin);
    int latex_ok = latex_is_ok();
    int i;

    if (gretl_stock_ifac == NULL) {
	gretl_stock_icons_init();
    }

    if (text_out || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->dialog), "text_out", GINT_TO_POINTER(1));
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gtk_toolbar_new();
    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);

    for (i=0; viewbar_items[i].str != NULL; i++) {
	struct viewbar_item *vitem = &viewbar_items[i];

	func = item_get_callback(vitem, vwin, latex_ok, sortby_ok);
	if (func == NULL) {
	    continue;
	}

	if (vitem->flag == PLOT_ITEM) {
	    set_plot_icon(vitem);
	}

	button = gtk_image_new();
	gtk_image_set_from_stock(GTK_IMAGE(button), vitem->icon, 
				 GTK_ICON_SIZE_MENU);
        w = gtk_toolbar_append_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(vitem->str), NULL,
				    button, G_CALLBACK(func), vwin);
	g_object_set_data(G_OBJECT(w), "flag", 
			  GINT_TO_POINTER(vitem->flag));

	if (vitem->flag == SAVE_ITEM) { 
	    /* nothing to save just yet */
	    gtk_widget_set_sensitive(w, FALSE);
	} else if (vitem->flag == SAVE_AS_ITEM &&
	    strstr(vwin->fname, "script_tmp")) {
	    gtk_widget_set_sensitive(w, FALSE);
	}
    }

    gtk_widget_show(vwin->mbar);
    gtk_widget_show(hbox);
}

static void add_edit_items_to_viewbar (windata_t *vwin)
{
    struct viewbar_item *vitem;
    GtkWidget *w, *button;
    int i, pos = 0;

    for (i=0; viewbar_items[i].str != NULL; i++) {
	vitem = &viewbar_items[i];
	if (vitem->flag == SAVE_ITEM ||
	    vitem->flag == EDIT_ITEM ||
	    vitem->flag == EDIT_SCRIPT_ITEM) {
	    button = gtk_image_new();
	    gtk_image_set_from_stock(GTK_IMAGE(button), 
				     vitem->icon, 
				     GTK_ICON_SIZE_MENU);
	    w = gtk_toolbar_insert_item(GTK_TOOLBAR(vwin->mbar),
					NULL, _(vitem->str), NULL,
					button, G_CALLBACK(vitem->func), 
					vwin, pos);
	    g_object_set_data(G_OBJECT(w), "flag", 
			      GINT_TO_POINTER(vitem->flag));
	    if (vitem->flag == SAVE_ITEM) { 
		gtk_widget_set_sensitive(w, FALSE);
	    } 
	}
	if (vitem->flag != GP_ITEM) {
	    pos++;
	}
    }
}

static GtkWidget *build_text_popup (windata_t *vwin)
{
    struct viewbar_item *vitem;
    toolfunc func;
    GtkWidget *pmenu = gtk_menu_new();
    GtkWidget *w;
    int i;

    for (i=0; viewbar_items[i].str != NULL; i++) {
	vitem = &viewbar_items[i];
	func = item_get_callback(vitem, vwin, 0, 0);
	if (func != NULL) {
	    if (func == text_paste_callback) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);

		if (!gtk_clipboard_wait_is_text_available(cb)) {
		    continue;
		}
	    } else if (func == text_undo_callback &&
		       !text_can_undo(vwin)) {
		continue;
	    }
	    w = gtk_menu_item_new_with_label(_(vitem->str));
	    g_signal_connect(G_OBJECT(w), "activate",
			     G_CALLBACK(func),
			     vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	}
    }

    return pmenu;
}

static gboolean 
text_popup_handler (GtkWidget *w, GdkEventButton *event, gpointer p)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);

    if (mods & GDK_BUTTON3_MASK) {
	windata_t *vwin = (windata_t *) p;

	if (vwin->popup) {
	    gtk_widget_destroy(vwin->popup);
	    vwin->popup = NULL;
	}

	vwin->popup = build_text_popup(vwin);

	if (vwin->popup != NULL) {
	    gtk_menu_popup(GTK_MENU(vwin->popup), NULL, NULL, NULL, NULL,
			   event->button, event->time);
	    gtk_signal_connect(GTK_OBJECT(vwin->popup), "destroy",
			       GTK_SIGNAL_FUNC(gtk_widget_destroyed), 
			       &vwin->popup);
	}

	return TRUE;
    }

    return FALSE;
}

static gchar *title_from_filename (const char *fname)
{
    const char *p = strrchr(fname, SLASH);
    gchar *trfname, *title = NULL;

    if (p != NULL) {
	trfname = my_filename_to_utf8(p + 1);
    } else {
	trfname = my_filename_to_utf8(fname);
    }

    title = g_strdup_printf("gretl: %s", trfname);

    g_free(trfname);

    return title;
}

static gchar *make_viewer_title (int role, const char *fname)
{
    gchar *title = NULL;

    switch (role) {
    case GUI_HELP: 
	title = g_strdup(_("gretl: help")); break;
    case FUNCS_HELP:
	title = g_strdup(_("gretl: function reference")); break;
    case CLI_HELP:
	title = g_strdup(_("gretl: command reference")); break;
    case GUI_HELP_EN: 
	title = g_strdup("gretl: help"); break;
    case CLI_HELP_EN:
	title = g_strdup("gretl: command reference"); break;
    case VIEW_LOG:
	title = g_strdup(_("gretl: command log")); break;
    case CONSOLE:
	title = g_strdup(_("gretl console")); break;
    case EDIT_SCRIPT:
    case VIEW_SCRIPT:	
    case VIEW_FILE:
	if (strstr(fname, "script_tmp") || strstr(fname, "session.inp")) {
	    title = g_strdup(_("gretl: command script"));
	} else {
	    title = title_from_filename(fname);
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

static void content_changed (GtkWidget *w, windata_t *vwin)
{
    mark_content_changed(vwin);
}

static void attach_content_changed_signal (windata_t *vwin)
{
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    g_signal_connect(G_OBJECT(tbuf), "changed", 
		     G_CALLBACK(content_changed), vwin);
}

static windata_t *common_viewer_new (int role, const char *title, 
				     gpointer data, int record)
{
    windata_t *vwin;

    vwin = mymalloc(sizeof *vwin);
    if (vwin == NULL) return NULL;

    windata_init(vwin);

    vwin->role = role;
    vwin->data = data;
    vwin->dialog = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(vwin->dialog), title);

    if (record) {
	g_object_set_data(G_OBJECT(vwin->dialog), "object", data);
	g_object_set_data(G_OBJECT(vwin->dialog), "role", 
			  GINT_TO_POINTER(vwin->role));
	winstack_add(vwin->dialog);
    }

    return vwin;
}

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);

#ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(vwin->dialog), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif

    gtk_container_add(GTK_CONTAINER(vwin->dialog), vwin->vbox);
}

static void view_buffer_insert_text (windata_t *vwin, PRN *prn)
{
    if (prn != NULL) {
	const char *buf = gretl_print_get_buffer(prn);

	if (vwin->role == VIEW_FUNC_CODE || vwin->role == EDIT_FUNC_CODE) {
	    sourceview_insert_buffer(vwin, buf);
	} else if (vwin->role == SCRIPT_OUT) {
	    textview_set_text_colorized(vwin->w, buf);
	} else {
	    textview_set_text(vwin->w, buf);
	}
    }
}

static windata_t *reuse_script_out (windata_t *vwin, PRN *prn)
{
    int sticky = (vwin->flags & VWIN_STICKY);
    GtkTextBuffer *buf;
    const char *newtext;

    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    if (sticky) {
	/* append to previous content */
	GtkTextMark *mark;
	GtkTextIter iter;

	gtk_text_buffer_get_end_iter(buf, &iter);
	mark = gtk_text_buffer_create_mark(buf, NULL, &iter, TRUE);
	textview_append_text_colorized(vwin->w, newtext, 1);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->w), 
				     mark, 0.0, TRUE, 0, 0.05);
	gtk_text_buffer_delete_mark(buf, mark);
    } else {
	/* replace previous content */
	gtk_text_buffer_set_text(buf, "", -1);
	textview_set_text_colorized(vwin->w, newtext);
	cursor_to_top(vwin);
    }

    gretl_print_destroy(prn);

    gtk_window_present(GTK_WINDOW(vwin->dialog));

    return vwin;
}

static gboolean nullify_script_out (GtkWidget *w, windata_t **pvwin)
{
    *pvwin = NULL;
    return FALSE;
}

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			const char *title, int role, 
			gpointer data) 
{
    static windata_t *script_out;
    windata_t *vwin;
    int record = (role != SCRIPT_OUT);

    if (role == SCRIPT_OUT && script_out != NULL) {
	return reuse_script_out(script_out, prn);
    }

    if (title != NULL) {
	vwin = common_viewer_new(role, title, data, record);
    } else {
	gchar *tmp = make_viewer_title(role, NULL);

	vwin = common_viewer_new(role, tmp, data, record);
	g_free(tmp);
    }

    if (vwin == NULL) return NULL;

    viewer_box_config(vwin);

    if (role == VAR || role == VECM || role == SYSTEM) {
	/* special case: use a text-based menu bar */
	set_up_viewer_menu(vwin->dialog, vwin, system_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
	gretl_object_ref(data, (role == SYSTEM)? GRETL_OBJ_SYS : GRETL_OBJ_VAR);
	add_system_menu_items(vwin, role);
    } else if (role == VIEW_FUNC_CODE || role == EDIT_FUNC_CODE) {
	make_viewbar(vwin, 0);
    } else if (role != IMPORT) {
	make_viewbar(vwin, 1);
    }

    if (role == VIEW_FUNC_CODE) {
	create_source(vwin, hsize, vsize, FALSE);
    } else if (role == EDIT_FUNC_CODE) {
	create_source(vwin, hsize, vsize, TRUE);
    } else {
	create_text(vwin, hsize, vsize, FALSE);
	if (role == PRINT || role == SCRIPT_OUT ||
	    role == VIEW_MODELTABLE) {
	    text_set_word_wrap(vwin->w, 0);
	}
    }

    text_table_setup(vwin->vbox, vwin->w);

    if (role == SCRIPT_OUT && data != NULL) {
	/* partial output window for script */
	vwin_add_child((windata_t *) data, vwin);
    }

    /* arrange for clean-up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    /* register destruction of script output viewer */
    if (role == SCRIPT_OUT) {
	g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
			 G_CALLBACK(nullify_script_out), &script_out);
	script_out = vwin;
    }

    /* insert and then free the text buffer */
    view_buffer_insert_text(vwin, prn);
    gretl_print_destroy(prn);

    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    if (role == EDIT_FUNC_CODE) {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
	attach_content_changed_signal(vwin);
	g_signal_connect(G_OBJECT(vwin->dialog), "delete-event", 
			 G_CALLBACK(query_save_text), vwin);
    } 

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(text_popup_handler), vwin);
    cursor_to_top(vwin);

    return vwin;
}

#define view_file_use_sourceview(r) (r == EDIT_SCRIPT || \
                                     r == VIEW_SCRIPT || \
                                     r == VIEW_LOG || \
                                     r == GR_PLOT)

#define doing_script(r) (r == EDIT_SCRIPT || \
			 r == VIEW_SCRIPT || \
			 r == VIEW_LOG)

windata_t *view_file (const char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role)
{
    windata_t *vwin;
    FILE *fp;
    gchar *title = NULL;

    /* first check that we can open the specified file */
    fp = gretl_fopen(filename, "r");
    if (fp == NULL) {
	errbox(_("Can't open %s for reading"), filename);
	return NULL;
    } else {
	fclose(fp);
    }

    /* then start building the file viewer */
    title = make_viewer_title(role, filename);
    vwin = common_viewer_new(role, (title != NULL)? title : filename, 
			     NULL, !doing_script(role) && role != CONSOLE);
    g_free(title);

    if (vwin == NULL) {
	return NULL;
    }

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);
    make_viewbar(vwin, (role == VIEW_DATA || role == CONSOLE || role == VIEW_FILE));

    if (view_file_use_sourceview(role)) {
	create_source(vwin, hsize, vsize, editable);
    } else {
	create_text(vwin, hsize, vsize, editable);
    }

    text_table_setup(vwin->vbox, vwin->w);

    if (view_file_use_sourceview(role)) {
	sourceview_insert_file(vwin, filename);
    } else {
	textview_insert_file(vwin, filename);
    }

    /* grab the "changed" signal when editing a script or graph */
    if (role == EDIT_SCRIPT || role == GR_PLOT) {
	attach_content_changed_signal(vwin);
    }

    /* catch some special keystrokes */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (editable) {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
    }

    /* alert for unsaved changes on exit */
    if (role == EDIT_SCRIPT || role == GR_PLOT) {
	g_signal_connect(G_OBJECT(vwin->dialog), "delete-event", 
			 G_CALLBACK(query_save_text), vwin);
    }

    /* clean up when dialog is destroyed */
    if (del_file) {
	gchar *fname = g_strdup(filename);

	g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
			 G_CALLBACK(delete_file), (gpointer) fname);
    }

    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(text_popup_handler), vwin);

    cursor_to_top(vwin);
    gtk_widget_grab_focus(vwin->w);

    return vwin;
}

windata_t *
view_help_file (const char *filename, int role, GtkItemFactoryEntry *menu_items)
{
    windata_t *vwin;
    gchar *fbuf = NULL;
    gchar *title = NULL;
    GError *gerr = NULL;
    int hsize = 80, vsize = 400;

    /* grab content of the appropriate help file into a buffer */
    g_file_get_contents(filename, &fbuf, NULL, &gerr);
    if (gerr != NULL) {
	errbox(gerr->message);
	g_error_free(gerr);
	return NULL;
    }

    title = make_viewer_title(role, NULL);
    vwin = common_viewer_new(role, title, NULL, 0);
    g_free(title);

    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);
    vwin->data = fbuf;

    viewer_box_config(vwin);
    set_up_viewer_menu(vwin->dialog, vwin, menu_items);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    if (role == FUNCS_HELP) {
	vsize = 500;
    }

    create_text(vwin, hsize, vsize, FALSE);
    text_table_setup(vwin->vbox, vwin->w);

    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (vwin->role == CLI_HELP || vwin->role == CLI_HELP_EN ||
	vwin->role == FUNCS_HELP) {
	g_signal_connect(G_OBJECT(vwin->w), "button_press_event",
			 G_CALLBACK(help_popup_handler), 
			 vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
			 G_CALLBACK(text_popup_handler), vwin);
    }	

    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    /* make the helpfile variant discernible via vwin->w */
    g_object_set_data(G_OBJECT(vwin->w), "role", GINT_TO_POINTER(vwin->role));

    gtk_widget_grab_focus(vwin->w);

    return vwin;
}

void view_window_set_editable (windata_t *vwin)
{
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), TRUE);
    g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
    g_signal_connect(G_OBJECT(vwin->dialog), "delete-event", 
		     G_CALLBACK(query_save_text), vwin);
    vwin->role = EDIT_SCRIPT;
    add_edit_items_to_viewbar(vwin);
    attach_content_changed_signal(vwin);
}

static gint query_save_text (GtkWidget *w, GdkEvent *event, 
			     windata_t *vwin)
{
    if (CONTENT_IS_CHANGED(vwin)) {
	int resp = yes_no_dialog("gretl", _("Save changes?"), 1);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	}

	if (resp == GRETL_YES) {
	    if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		buf_edit_save(NULL, vwin);
	    } else if (vwin->role == EDIT_SCRIPT) {
		auto_save_script(vwin);
	    } else if (vwin->role == GR_PLOT) {
		auto_save_plot(vwin);
	    }
	}
    }

    return FALSE;
}

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role) 
{
    windata_t *vwin;

    vwin = common_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) {
	return NULL;
    }

    viewer_box_config(vwin); 

    /* add a menu bar */
    make_viewbar(vwin, 0);

    create_text(vwin, hsize, vsize, TRUE);
    text_table_setup(vwin->vbox, vwin->w);
    
    /* insert the buffer text */
    if (*pbuf) {
	GtkTextBuffer *tbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	gtk_text_buffer_set_text(tbuf, *pbuf, -1);
    }
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(text_popup_handler), vwin);
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    attach_content_changed_signal(vwin);

    /* alert for unsaved changes on exit */
    g_signal_connect(G_OBJECT(vwin->dialog), "delete-event",
		     G_CALLBACK(query_save_text), vwin);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    cursor_to_top(vwin);

    return vwin;
}

static gint 
check_delete_model_window (GtkWidget *w, GdkEvent *e, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gint ret = FALSE;

    if (window_is_busy(vwin)) {
	maybe_raise_dialog();
	ret = TRUE;
    }

    return ret;
}

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title) 
{
    windata_t *vwin;
    const char *buf;

    vwin = common_viewer_new(VIEW_MODEL, title, pmod, 1);
    if (vwin == NULL) {
	return 1;
    }

    /* Take responsibility for one reference to this model */
    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    viewer_box_config(vwin);

    set_up_viewer_menu(vwin->dialog, vwin, model_items);

    if (pmod->ci != MLE && pmod->ci != GMM) {
	add_vars_to_plot_menu(vwin);
	add_model_dataset_items(vwin);
    }

    if (latex_is_ok() && !pmod->errcode) {
	add_model_tex_items(vwin);
    }    

    if (pmod->ci != ARMA && pmod->ci != GARCH && 
	pmod->ci != NLS && pmod->ci != MLE && pmod->ci != GMM &&
	pmod->ci != PANEL && pmod->ci != ARBOND) {
	add_dummies_to_plot_menu(vwin);
    }

    g_signal_connect(G_OBJECT(vwin->mbar), "button_press_event", 
		     G_CALLBACK(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    create_text(vwin, hsize, vsize, FALSE);
    text_table_setup(vwin->vbox, vwin->w);

    /* insert and then free the model results buffer */
    buf = gretl_print_get_buffer(prn);
    textview_set_text(vwin->w, buf);
    gretl_print_destroy(prn);

    /* attach shortcuts */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(text_popup_handler), vwin);

    /* don't allow deletion of model window when a model
       test dialog is active */
    g_signal_connect(G_OBJECT(vwin->dialog), "delete-event", 
		     G_CALLBACK(check_delete_model_window), 
		     vwin);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), 
		     vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show_all(vwin->dialog);

    cursor_to_top(vwin);

    return 0;
}

static void auto_save_plot (windata_t *vwin)
{
    gchar *buf;

    buf = textview_get_text(vwin->w);
    if (buf == NULL) {
	return;
    }

    dump_plot_buffer(buf, vwin->fname, 1);
    g_free(buf);

    mark_content_saved(vwin);
}

static void auto_save_script (windata_t *vwin)
{
    FILE *fp;
    gchar *savestuff;
    int unsaved = 0;

    if (strstr(vwin->fname, "script_tmp") || *vwin->fname == '\0') {
	file_save(vwin, SAVE_SCRIPT, NULL);
	strcpy(vwin->fname, scriptfile);
	unsaved = 1;
    }

    if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	file_write_errbox(vwin->fname);
	return;
    }

    savestuff = textview_get_text(vwin->w);
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);

    mark_content_saved(vwin);
}

static void model_tex_equation_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/LaTeX/View/Equation", s);
    flip(ifac, "/LaTeX/Save/Equation", s);
    flip(ifac, "/LaTeX/Copy/Equation", s);
    flip(ifac, "/LaTeX/Equation options", s);
}

static void copy_no_underscore (char *targ, const char *src)
{
    while (*src) {
	if (*src != '_') {
	    *targ++ = *src;
	}
	src++;
    }
    *targ = 0;
}

static void set_tests_menu_state (GtkItemFactory *ifac, const MODEL *pmod)
{
    gretlopt opt = OPT_NONE;
    char path[128];
    int i, cmd_ci, ok;

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->ci == MPOLS) { 
	/* FIXME? */
	flip(ifac, "/Tests", FALSE);
	return;
    }

    for (i=0; model_items[i].path != NULL; i++) {
	if (model_items[i].item_type == NULL &&
	    strstr(model_items[i].path, "Tests")) {
	    cmd_ci = model_items[i].callback_action;

	    if (cmd_ci == LMTEST_SQUARES) {
		cmd_ci = LMTEST;
		opt = OPT_S;
	    } else if (cmd_ci == LMTEST_LOGS) {
		cmd_ci = LMTEST;
		opt = OPT_L;
	    } else if (cmd_ci == LMTEST_WHITE) {
		cmd_ci = LMTEST;
		opt = OPT_W;
	    } else if (cmd_ci == LMTEST_BP) {
		cmd_ci = LMTEST;
		opt = OPT_B;
	    } else if (cmd_ci == LMTEST_BPK) {
		cmd_ci = LMTEST;
		opt = OPT_B | OPT_R;
	    } else if (cmd_ci == ARCH) {
		cmd_ci = LMTEST;
		opt = OPT_H;
	    } else if (cmd_ci == LMTEST) { 
		/* unqualified: autocorrelation */
		opt = OPT_A;
	    } else if (cmd_ci == CUSUMSQ) {
		cmd_ci = CUSUM;
	    }
		
	    ok = model_test_ok(cmd_ci, opt, pmod, datainfo);
	    copy_no_underscore(path, model_items[i].path);
	    flip(ifac, path, ok);
	}
    }

    if (pmod->ncoeff == 1) {
	flip(ifac, "/Analysis/Confidence ellipse...", FALSE);
    }
}

static void model_save_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/File/Save to session as icon", s);
    flip(ifac, "/File/Save as icon and close", s);
}

static void arma_x12_menu_mod (windata_t *vwin)
{
    flip(vwin->ifac, "/Analysis/Coefficient covariance matrix", FALSE);
    add_x12_output_menu_item(vwin);
}

static void adjust_model_menu_state (windata_t *vwin, const MODEL *pmod)
{
    set_tests_menu_state(vwin->ifac, pmod);

    /* disallow saving an already-saved model */
    if (pmod->name != NULL) {
	model_save_state(vwin->ifac, FALSE);
    }

    if (pmod->ci == MLE || pmod->ci == GMM) {
	/* some of this could be relaxed later */
	flip(vwin->ifac, "/Analysis", FALSE);
	flip(vwin->ifac, "/Graphs", FALSE);
    } else if (pmod->ci == ARMA && arma_by_x12a(pmod)) {
	arma_x12_menu_mod(vwin);
    } 

    if (pmod->ci == GMM) {
	flip(vwin->ifac, "/Save", FALSE);
    }

    if (dataset_is_panel(datainfo) && pmod->ci == OLS) {
	panel_heteroskedasticity_menu(vwin);
    }

    if (pmod->ci == ARBOND) {
	flip(vwin->ifac, "/Analysis/Forecasts...", FALSE);
    }

    if (pmod->ci != OLS || !pmod->ifc || na(pmod->ess) ||
	na(pmod->tss)) {
	flip(vwin->ifac, "/Analysis/ANOVA", FALSE);
    }

    if (!bootstrap_ok(pmod->ci)) {
	flip(vwin->ifac, "/Analysis/Bootstrap...", FALSE);
    }
}

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[])
{
    gint n_items = 0;

    while (items[n_items].path != NULL) n_items++;

    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", NULL);

#ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(vwin->ifac, menu_translate, NULL, NULL);
#endif
    gtk_item_factory_create_items(vwin->ifac, n_items, items, vwin);
    vwin->mbar = gtk_item_factory_get_widget(vwin->ifac, "<main>");

    if (vwin->data == NULL) {
	return;
    }

    if (vwin->role == VIEW_MODEL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	adjust_model_menu_state(vwin, pmod);
    } else if (vwin->role == VAR || vwin->role == VECM || 
	       vwin->role == SYSTEM) {
	model_save_state(vwin->ifac, !is_session_model(vwin->data));
    }
}

static GtkItemFactoryEntry model_dataset_basic_items[] = {
    { N_("/Save/_Fitted values"), NULL, 
      fit_resid_callback, GENR_FITTED, NULL, GNULL },
    { N_("/Save/_Residuals"), NULL, 
      fit_resid_callback, GENR_RESID, NULL, GNULL },
    { N_("/Save/_Squared residuals"), NULL, 
      fit_resid_callback, GENR_RESID2, NULL, GNULL }
};

static GtkItemFactoryEntry ess_items[] = {
    { N_("/Save/_Error sum of squares"), NULL, 
      model_stat_callback, ESS, NULL, GNULL },
    { N_("/Save/_Standard error of residuals"), NULL, 
      model_stat_callback, SIGMA, NULL, GNULL }
}; 

static GtkItemFactoryEntry r_squared_items[] = {
    { N_("/Save/_R-squared"), NULL, 
      model_stat_callback, R2, NULL, GNULL },
    { N_("/Save/_T*R-squared"), NULL, 
      model_stat_callback, TR2, NULL, GNULL }
};   

static GtkItemFactoryEntry lnl_data_item = {
    N_("/Save/_Log likelihood"), NULL, 
    model_stat_callback, LNL, NULL, GNULL 
};

static GtkItemFactoryEntry criteria_items[] = {
    { N_("/Save/_Akaike Information Criterion"), NULL, 
      model_stat_callback, AIC, NULL, GNULL },
    { N_("/Save/_Bayesian Information Criterion"), NULL, 
      model_stat_callback, BIC, NULL, GNULL },
    { N_("/Save/_Hannan-Quinn Information Criterion"), NULL, 
      model_stat_callback, HQC, NULL, GNULL }
};

static GtkItemFactoryEntry garch_data_item = {
    N_("/Save/_Predicted error variance"), NULL, 
    fit_resid_callback, GENR_H, NULL, GNULL 
};

static GtkItemFactoryEntry fixed_effects_data_item = {
    N_("/Save/Per-unit _constants"), NULL, 
    fit_resid_callback, GENR_AHAT, NULL, GNULL 
};

static GtkItemFactoryEntry define_var_items[] = {
    { "/Save/sep1", NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Save/Define _new variable..."), NULL, model_genr_callback,
      MODEL_GENR, NULL, GNULL }
};

static int criteria_available (const MODEL *pmod)
{
    int i;

    for (i=0; i<C_MAX; i++) {
	if (na(pmod->criterion[i])) {
	    return 0;
	}
    }

    return 1;
}

static void add_model_dataset_items (windata_t *vwin)
{
    MODEL *pmod = vwin->data;
    int i, n;

    n = sizeof model_dataset_basic_items / 
	sizeof model_dataset_basic_items[0];

    for (i=0; i<n; i++) {
	gtk_item_factory_create_item(vwin->ifac, &model_dataset_basic_items[i], 
				     vwin, 1);
    }

    if (gretl_model_get_data(pmod, "ahat") != NULL) {
	gtk_item_factory_create_item(vwin->ifac, &fixed_effects_data_item, vwin, 1);
    }

    if (pmod->ci != GARCH) {
	n = sizeof ess_items / sizeof ess_items[0];
	for (i=0; i<n; i++) {
	    gtk_item_factory_create_item(vwin->ifac, &ess_items[i], vwin, 1);
	}
    }

    if (!ML_ESTIMATOR(pmod->ci) && pmod->ci != LAD && !na(pmod->rsq)) {
	n = sizeof r_squared_items / sizeof r_squared_items[0];
	for (i=0; i<n; i++) {
	    gtk_item_factory_create_item(vwin->ifac, &r_squared_items[i], 
					 vwin, 1);
	}
    }

    if (ML_ESTIMATOR(pmod->ci)) {
	gtk_item_factory_create_item(vwin->ifac, &lnl_data_item, vwin, 1);
    }

    if (criteria_available(pmod)) {
	n = sizeof criteria_items / sizeof criteria_items[0];
	for (i=0; i<n; i++) {
	    gtk_item_factory_create_item(vwin->ifac, &criteria_items[i], vwin, 1);
	}
    }

    if (pmod->ci == GARCH) {
	gtk_item_factory_create_item(vwin->ifac, &garch_data_item, vwin, 1);
    }

    for (i=0; i<2; i++) {
	gtk_item_factory_create_item(vwin->ifac, &define_var_items[i], 
				     vwin, 1);
    }
}

static void add_model_tex_items (windata_t *vwin)
{
    int i, n = sizeof model_tex_items / sizeof model_tex_items[0];
    MODEL *pmod = (MODEL *) vwin->data;
    int eqn_ok = command_ok_for_model(EQNPRINT, 0, pmod->ci);
    GtkWidget *w;

    for (i=0; i<n; i++) {
	gtk_item_factory_create_item(vwin->ifac, &model_tex_items[i], 
				     vwin, 1);
    }  

    model_tex_equation_state(vwin->ifac, !pmod->errcode && eqn_ok);
    w = gtk_item_factory_get_widget(vwin->ifac, 
				    "/LaTeX/Equation options/Show t-ratios");
    gtk_check_menu_item_set_active(GTK_CHECK_MENU_ITEM(w),
				   (gboolean) get_tex_eqn_opt());
}

#define VNAMELEN2 32

static void add_vars_to_plot_menu (windata_t *vwin)
{
    int i, j, varstart;
    GtkItemFactoryEntry varitem;
    const gchar *mpath[] = {
	N_("/Graphs/Residual plot"), 
	N_("/Graphs/Fitted, actual plot")
    };
    MODEL *pmod = vwin->data;
    char tmp[VNAMELEN2];

   varitem.accelerator = NULL; 
   varitem.item_type = NULL;
   varitem.callback_action = 0; 

    for (i=0; i<2; i++) {
	/* residual correlogram and spectrum */
	if (dataset_is_time_series(datainfo) && i == 0) {
	    varitem.path = g_strdup_printf(_("%s/_Correlogram"), mpath[i]);
	    varitem.callback = residual_correlogram;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	    varitem.path = g_strdup_printf(_("%s/_Spectrum"), mpath[i]);
	    varitem.callback = residual_periodogram;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}

	/* plot against time/obs number */
	if (dataset_is_time_series(datainfo)) {
	    varitem.path = g_strdup_printf(_("%s/_Against time"), mpath[i]);
	} else {
	    varitem.path = g_strdup_printf(_("%s/By _observation number"), mpath[i]);
	}
	varitem.callback = (i==0)? resid_plot : fit_actual_plot;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	if (pmod->ci == ARMA || pmod->ci == NLS || pmod->ci == GARCH ||
	    pmod->ci == PANEL || pmod->ci == ARBOND) { 
	    continue;
	}

	varstart = (i == 0)? 1 : 2;

	/* put the indep vars on the menu list */
	for (j=varstart; j<=pmod->list[0]; j++) {
	    if (pmod->list[j] == 0) continue;
	    if (pmod->list[j] == LISTSEP) break;
	    if (!strcmp(datainfo->varname[pmod->list[j]], "time")) {
		continue;
	    }

	    varitem.callback_action = pmod->list[j]; 
	    double_underscores(tmp, datainfo->varname[pmod->list[j]]);
	    varitem.path = 
		g_strdup_printf(_("%s/_Against %s"), mpath[i], tmp);
	    varitem.callback = (i == 0)? resid_plot : fit_actual_plot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}

	varitem.callback_action = 0;

	/* if the model has two independent vars, offer a 3-D fitted
	   versus actual plot */
	if (i == 1 && pmod->ifc && pmod->ncoeff == 3) {
	    char tmp2[VNAMELEN2];

	    double_underscores(tmp, datainfo->varname[pmod->list[3]]);
	    double_underscores(tmp2, datainfo->varname[pmod->list[4]]);
	    varitem.path =
		g_strdup_printf(_("%s/_Against %s and %s"),
				mpath[i], tmp, tmp2);
	    varitem.callback = fit_actual_splot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}	
    }
}

static void plot_dummy_call (gpointer data, guint v, GtkWidget *widget)
{
    GtkCheckMenuItem *item = GTK_CHECK_MENU_ITEM(widget);
    windata_t *vwin = (windata_t *) data;

    if (item->active) vwin->active_var = v; 
}

static void add_dummies_to_plot_menu (windata_t *vwin)
{
    GtkItemFactoryEntry dumitem;
    MODEL *pmod = vwin->data;
    const gchar *mpath[] = {
	"/Graphs/dumsep", 
	N_("/Graphs/Separation")
    };
    gchar *radiopath = NULL;
    char tmp[VNAMELEN2];
    int i, done_branch = 0;

    dumitem.path = NULL;
    dumitem.accelerator = NULL; 

    /* put the dummy independent vars on the menu list */
    for (i=2; i<=pmod->list[0]; i++) {

	if (pmod->list[i] == LISTSEP) {
	    break;
	}

	if (pmod->list[i] == 0 ||
	    !gretl_isdummy(datainfo->t1, datainfo->t2, Z[pmod->list[i]])) {
	    continue;
	}

	if (!done_branch) {
	    /* add separator */
	    dumitem.callback = NULL;
	    dumitem.callback_action = 0;
	    dumitem.item_type = "<Separator>";
	    dumitem.path = g_strdup_printf(_("%s"), mpath[0]);
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    g_free(dumitem.path);

	    /* add menu branch */
	    dumitem.item_type = "<Branch>";
	    dumitem.path = g_strdup_printf(_("%s"), mpath[1]);
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    g_free(dumitem.path);

	    /* add "none" option */
	    dumitem.callback = plot_dummy_call;
	    dumitem.item_type = "<RadioItem>";
	    dumitem.path = g_strdup_printf(_("%s/none"), mpath[1]);
	    radiopath = g_strdup(dumitem.path);
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    g_free(dumitem.path);

	    done_branch = 1;
	} 

	dumitem.callback_action = pmod->list[i]; 
	double_underscores(tmp, datainfo->varname[pmod->list[i]]);
	dumitem.callback = plot_dummy_call;	    
	dumitem.item_type = radiopath;
	dumitem.path = g_strdup_printf(_("%s/By %s"), mpath[1], tmp);
	gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	g_free(dumitem.path);
    }

    g_free(radiopath);
}

static void x12_output_callback (gpointer p, guint v, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = vwin->data;
    char *fname;

    if (pmod == NULL) return;

    fname = gretl_model_get_data(pmod, "x12a_output");
    if (fname != NULL) {
	char *p = strrchr(fname, '.');

	if (p != NULL && strlen(p) == 7) {
	    gchar *tmp = g_strdup(fname);

	    sprintf(p, ".%d", pmod->ID);
	    rename(tmp, fname);
	    g_free(tmp);
	}
	view_file(fname, 0, 0, 78, 350, VIEW_FILE);
    }
}

enum {
    SYS_DATA_RESIDS,
    SYS_DATA_FITTED,
    SYS_DATA_VCV
};

static void system_data_callback (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *M = NULL;
    gchar *wtitle = NULL;
    PRN *prn;
    int k = 0, err = 0;

    if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    } else {
	var = (GRETL_VAR *) vwin->data;
    } 

    if ((var == NULL && sys == NULL) || bufopen(&prn)) {
	return;
    }

    if (code == SYS_DATA_VCV) {
	if (var != NULL) {
	    wtitle = g_strdup(_("gretl: VAR covariance matrix"));
	    err = gretl_VAR_print_VCV(var, prn);
	} else {
	    wtitle = g_strdup(_("gretl: system covariance matrix"));
	    err = system_print_VCV(sys, prn);
	}
    } else if (code == SYS_DATA_RESIDS || code == SYS_DATA_FITTED) {
	const char *titles[] = {
	    N_("System residuals"),
	    N_("System fitted values")
	};
	const char *title;
	const char **heads = NULL;

	if (var != NULL) {
	    /* fitted values matrix not currently available */
	    M = (code == SYS_DATA_RESIDS)? gretl_VAR_get_residual_matrix(var) :
		NULL;
	} else {
	    M = (code == SYS_DATA_RESIDS)? sys->E : sys->yhat;
	}

	if (M == NULL) {
	    err = E_DATA;
	} else {
	    k = gretl_matrix_cols(M);
	    heads = malloc(k * sizeof *heads);
	    if (heads == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    int i, v;

	    for (i=0; i<k && !err; i++) {
		v = (var != NULL)? gretl_VAR_get_variable_number(var, i) :
		    sys->lists[i][1];
		if (v < 0 || v >= datainfo->v) {
		    err = E_DATA;
		} else {
		    heads[i] = datainfo->varname[v];
		}
	    }
	}

	if (!err) {
	    title = (code == SYS_DATA_RESIDS)? titles[0] : titles[1];
	    wtitle = g_strdup_printf("gretl: %s", _(title));
	    gretl_matrix_print_with_col_heads(M, _(title), heads, prn);
	}

	free(heads);
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	/* FIXME: add matrix as saveable data */
	view_buffer(prn, 80, 400, wtitle, PRINT, NULL);
    }

    g_free(wtitle);
}

static void VAR_model_data_callback (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = vwin->data;
    gchar *title;
    PRN *prn;
    int h = 0;
    int err;

    if (var == NULL) return;

    if (bufopen(&prn)) return;

    h = default_VAR_horizon(datainfo);
    title = g_strdup_printf("gretl: %s", 
			    (code == VAR_IRF)? _("impulse responses") :
			    _("variance decompositions"));
    err = checks_dialog(title, NULL, NULL, 0, NULL, 0, NULL,
			&h, _("forecast horizon (periods):"),
			2, datainfo->n / 2, 0);
    g_free(title);

    if (err < 0) {
	gretl_print_destroy(prn);
	return;
    } 

    if (code == VAR_IRF) {
	title = g_strdup(_("gretl: VAR impulse responses"));
	err = gretl_VAR_print_all_impulse_responses(var, datainfo, h, prn);
    } else if (code == VAR_DECOMP) {
	title = g_strdup(_("gretl: VAR variance decompositions"));
	err = gretl_VAR_print_all_fcast_decomps(var, datainfo, h, prn);
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	windata_t *viewer;

	viewer = view_buffer(prn, 80, 400, title, code, NULL);
	vwin_add_child(vwin, viewer);
	viewer->active_var = h;
    }

    g_free(title);
}

static void panel_heteroskedasticity_menu (windata_t *vwin)
{
    const gchar *tpath = "/Tests";
    GtkItemFactoryEntry hitem;

    gtk_item_factory_delete_item(vwin->ifac, "/Tests/Heteroskedasticity");

    hitem.accelerator = NULL;
    hitem.item_type = NULL;

    hitem.callback = do_lmtest;
    hitem.callback_action = LMTEST_WHITE;
    hitem.path = g_strdup_printf("%s/%s", tpath, _("Heteroskedasticity (_White's test)"));
    gtk_item_factory_create_item(vwin->ifac, &hitem, vwin, 1);
    g_free(hitem.path);

    hitem.callback = do_lmtest;
    hitem.callback_action =  LMTEST_GROUPWISE;
    hitem.path = g_strdup_printf("%s/%s", tpath, ("Heteroskedasticity (_groupwise)"));
    gtk_item_factory_create_item(vwin->ifac, &hitem, vwin, 1);
    g_free(hitem.path);
}

static void add_x12_output_menu_item (windata_t *vwin)
{
    GtkItemFactoryEntry item;
    const gchar *mpath = "/Analysis";

    item.accelerator = NULL; 
    item.callback_action = 0;

    /* separator */
    item.callback = NULL;
    item.item_type = "<Separator>";
    item.path = g_strdup_printf("%s/%s", mpath, _("x12sep"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    /* actual item */
    item.callback = x12_output_callback;
    item.item_type = NULL;
    item.path = g_strdup_printf("%s/%s", mpath, _("View X-12-ARIMA output"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);
}

static int 
impulse_response_setup (GRETL_VAR *var, int *horizon, int *bootstrap)
{
    gchar *title;
    int h = default_VAR_horizon(datainfo);
    const char *impulse_opts[] = {
	N_("include bootstrap confidence interval")
    };
    static int active[] = { 0 };
    int err;

    if (restricted_VECM(var)) {
	active[0] = -1;
    }

    title = g_strdup_printf("gretl: %s", _("impulse responses"));

    err = checks_dialog(title, NULL,
			impulse_opts, 
			1, 
			active,
			0, NULL,
			&h, _("forecast horizon (periods):"),
			2, datainfo->n / 2, IRF_BOOT);
    g_free(title);

    if (err < 0) {
	/* cancelled */
	*horizon = 0;
    } else {
	*horizon = h;
	*bootstrap = (active[0] > 0);
    }

    return err;
}

static void impulse_plot_call (gpointer p, guint shock, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    gint targ;
    const double **vZ = NULL;
    int err;

    targ = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "targ"));

    if (impulse_response_setup(var, &horizon, &bootstrap) < 0) {
	return;
    }

    if (bootstrap) {
	vZ = (const double **) Z;
    }

    err = gretl_VAR_plot_impulse_response(var, targ, shock, horizon,
					  vZ, datainfo);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static void multiple_irf_plot_call (gpointer p, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    const double **vZ = NULL;
    int err;

    if (impulse_response_setup(var, &horizon, &bootstrap) < 0) {
	return;
    }

    if (bootstrap) {
	vZ = (const double **) Z;
    }    

    err = gretl_VAR_plot_multiple_irf(var, horizon, vZ, datainfo);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static void system_forecast_callback (gpointer p, guint i, GtkWidget *w)
{
    static gretlopt gopt = OPT_P;
    windata_t *vwin = (windata_t *) p;
    int ci = vwin->role;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    FITRESID *fr;
    int t1, t2, t2est, resp;
    int premax, pre_n, dyn_ok;
    int static_model = 0;
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) vwin->data;
	t2est = gretl_VAR_get_t2(var);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) vwin->data;
	t2est = sys->t2;
	static_model = (sys->order == 0);
    } else {
	return;
    }

    t2 = datainfo->n - 1;

    /* if no out-of-sample obs are available, alert the user */
    if (t2 == t2est) {
	err = out_of_sample_info(1, &t2);
	if (err) {
	    return;
	}
	t2 = datainfo->n - 1;
    }

    /* max number of pre-forecast obs in "best case" */
    premax = datainfo->n - 1;

    /* if there are spare obs available, default to an
       out-of-sample forecast */
    if (t2 > t2est) {
	t1 = t2est + 1;
	pre_n = t2est / 2;
	if (pre_n > 100) {
	    pre_n = 100;
	}
	dyn_ok = !static_model;
    } else {
	if (var != NULL) {
	    t1 = effective_order(var);
	} else {
	    t1 = sys->order;
	}
	pre_n = 0;
	dyn_ok = 0;
    }

    /* FIXME pre_n with static fcast? */

    resp = forecast_dialog(t1, t1, &t1,
			   t1, t2, &t2,
			   0, premax, &pre_n,
			   dyn_ok, &gopt, NULL);
    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    }

    fr = get_system_forecast(vwin->data, ci, i, t1, t2, pre_n,
			     (const double **) Z, datainfo, 
			     opt, &err);

    if (err) {
	gui_errmsg(err);
    } else {
	int width = 78;
	PRN *prn;

	if (bufopen(&prn)) {
	    return;
	}

	err = text_print_forecast(fr, datainfo, gopt, prn);
	if (!err) {
	    register_graph();
	}
	if (fr->sderr == NULL) {
	    width = 50;
	}
	view_buffer(prn, width, 400, _("gretl: forecasts"), FCAST, fr);
    }
}

enum {
    SYS_AUTOCORR_TEST,
    SYS_ARCH_TEST,
    SYS_NORMALITY_TEST,
    SYS_RESTRICT
};

static void system_test_call (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    char title[72];
    PRN *prn;
    int order = 0;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    } else {
	var = (GRETL_VAR *) vwin->data;
    }

    if (code == SYS_AUTOCORR_TEST || code == SYS_ARCH_TEST) {
	order = default_lag_order(datainfo);
	set_window_busy(vwin);
	err = spin_dialog((code == SYS_AUTOCORR_TEST)?
			  _("gretl: autocorrelation") :
			  _("gretl: ARCH test"), NULL,
			  &order, _("Lag order for test:"),
			  1, datainfo->n / 2, LMTEST);
	unset_window_busy(vwin);
	if (err < 0) {
	    gretl_print_destroy(prn);
	    return;
	}
    }	

    if (code == SYS_AUTOCORR_TEST) {
	strcpy(title, _("gretl: autocorrelation"));
	if (var != NULL) {
	    err = gretl_VAR_autocorrelation_test(var, order, 
						 &Z, datainfo, 
						 prn);
	} else {
	    err = system_autocorrelation_test(sys, order, prn);
	}
    } else if (code == SYS_ARCH_TEST) {
	strcpy(title, _("gretl: ARCH test"));
	if (var != NULL) {
	    err = gretl_VAR_arch_test(var, order, datainfo, prn);
	} else {
	    err = system_arch_test(sys, order, prn);
	}
    } else if (code == SYS_NORMALITY_TEST) {
	sprintf(title, "gretl: %s", _("Test for normality of residual"));
	if (var != NULL) {
	    err = gretl_VAR_normality_test(var, prn);
	} else {
	    err = system_normality_test(sys, prn);
	}
    } else {
	err = 1;
    }

    if (err) {
	gui_errmsg(err);
	gretl_print_destroy(prn);
    } else {
	view_buffer(prn, 78, 400, title, PRINT, NULL); 
    }
}

static void VAR_roots_plot_call (gpointer p, guint u, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_roots_plot(var);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void system_resid_plot_call (gpointer p, guint ci, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    int err;

    err = gretl_system_residual_plot(vwin->data, ci, datainfo);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void system_resid_mplot_call (gpointer p, guint ci, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    int err;

    err = gretl_system_residual_mplot(vwin->data, ci, datainfo);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void add_system_menu_items (windata_t *vwin, int ci)
{
    GtkItemFactoryEntry item;
    const gchar *tpath = N_("/Tests");
    const gchar *gpath = N_("/Graphs");
    const gchar *mpath = N_("/Analysis");
    const gchar *fpath = N_("/Analysis/Forecasts");
    const gchar *dpath = N_("/Save");
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    int neqns, nfc, vtarg, vshock;
    char tmp[VNAMELEN2];
    int i, j;

    if (ci == SYSTEM) {
	sys = (equation_system *) vwin->data;
	neqns = sys->neqns;
	nfc = sys->neqns + sys->nidents;
    } else {
	var = (GRETL_VAR *) vwin->data;
	nfc = neqns = gretl_VAR_get_n_equations(var);
    }

    item.accelerator = NULL;
    item.callback = NULL;
    item.callback_action = 0;
    item.item_type = "<Branch>";

    item.path = g_strdup(_("/_Tests"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    item.path = g_strdup(_("/_Analysis"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    item.path = g_strdup(_("/_Graphs"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    item.path = g_strdup(_(fpath));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);
    
    item.path = g_strdup(_(dpath));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    /* FIXME: the following two tests should really be multivariate */

    /* univariate autocorrelation tests */
    item.path = g_strdup_printf("%s/%s", _(tpath), 
				_("Autocorrelation"));
    item.callback = system_test_call;
    item.callback_action = SYS_AUTOCORR_TEST;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    /* univariate ARCH tests */
    item.path = g_strdup_printf("%s/%s", _(tpath), 
				_("ARCH"));
    item.callback = system_test_call;
    item.callback_action = SYS_ARCH_TEST;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    /* multivariate normality test */
    item.path = g_strdup_printf("%s/%s", _(tpath), 
				_("Normality of residuals"));
    item.callback = system_test_call;
    item.callback_action = SYS_NORMALITY_TEST;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    if (ci == VECM || ci == SYSTEM) {
	/* linear restrictions (on cointegrating relations, for VECM) */
	item.path = g_strdup_printf("%s/%s", _(tpath), 
				    _("Linear restrictions"));
	item.callback = gretl_callback;
	item.callback_action = RESTRICT;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);
    } else if (ci == VAR) {
	/* regular VAR: omit exogenous variables test */
	const int *xlist;

	xlist = gretl_VAR_get_exo_list(var);
	if (xlist != NULL) {
	    item.path = g_strdup_printf("%s/%s", _(tpath), 
					_("Omit exogenous variables..."));
	    item.callback = selector_callback;
	    item.callback_action = VAROMIT;
	    item.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	    g_free(item.path);
	}	    
    }

    /* Display residual matrix */
    item.path = g_strdup_printf("%s/%s", _(mpath), 
				_("Display residuals, all equations"));
    item.callback = system_data_callback;
    item.callback_action = SYS_DATA_RESIDS;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    if (ci == SYSTEM) {
	/* Display fitted values matrix */
	item.path = g_strdup_printf("%s/%s", _(mpath), 
				    _("Display fitted values, all equations"));
	item.callback = system_data_callback;
	item.callback_action = SYS_DATA_FITTED;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);  
    }  

    /* Display VCV matrix */
    item.path = g_strdup_printf("%s/%s", _(mpath), 
				_("Cross-equation covariance matrix"));
    item.callback = system_data_callback;
    item.callback_action = SYS_DATA_VCV;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    if (ci == VAR || ci == VECM) {
	/* impulse response printout */
	item.path = g_strdup_printf("%s/%s", _(mpath), _("Impulse responses"));
	item.callback = VAR_model_data_callback;
	item.callback_action = VAR_IRF;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);    

	/* variance decomp printout */
	item.path = g_strdup_printf("%s/%s", _(mpath), 
				    _("Forecast variance decomposition"));
	item.callback = VAR_model_data_callback;
	item.callback_action = VAR_DECOMP;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path); 
    }

    if (neqns <= 6) {
	/* separate residual plot */
	item.path = g_strdup_printf("%s/%s", _(gpath), _("Residual plots"));
	item.callback = system_resid_mplot_call;
	item.callback_action = ci;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);
    }

    /* combined residual plot */
    item.path = g_strdup_printf("%s/%s", _(gpath), _("Combined residual plot"));
    item.callback = system_resid_plot_call;
    item.callback_action = ci;
    item.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);

    if (ci != SYSTEM) {
	/* VAR inverse roots */
	item.path = g_strdup_printf("%s/%s", _(gpath), _("VAR inverse roots"));
	item.callback = VAR_roots_plot_call;
	item.callback_action = 0;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);
    }

    if (ci != SYSTEM && neqns <= 4) {
	/* Multiple IRFs */
	item.path = g_strdup_printf("%s/%s", _(gpath), _("Impulse responses (combined)"));
	item.callback = multiple_irf_plot_call;
	item.callback_action = 0;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);
    }

    for (i=0; i<nfc; i++) {
	char maj[64], min[32];
	int dv;

	/* forecast items */
	if (var != NULL) {
	    dv = gretl_VAR_get_variable_number(var, i);
	} else {
	    dv = sys->ylist[i+1];
	}
	double_underscores(tmp, datainfo->varname[dv]);
	item.path = g_strdup_printf("%s/%s", _(fpath), tmp);
	item.callback = system_forecast_callback;
	item.callback_action = i;
	item.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);

	if (i < neqns) {
	    /* save resids items */
	    item.path = g_strdup_printf("%s/%s %d", _(dpath), 
					_("Residuals from equation"), i + 1);
	    item.callback = add_system_resid;
	    item.callback_action = i;
	    item.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	    g_free(item.path);
	}

	if (var == NULL) {
	    continue;
	}

	/* impulse response plots: make branch for target */
	vtarg = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, datainfo->varname[vtarg]);
	sprintf(maj, _("Response of %s"), tmp);

	item.path = g_strdup_printf("%s/%s", _(gpath), maj);
	item.callback = NULL;
	item.callback_action = 0;
	item.item_type = "<Branch>";
	gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	g_free(item.path);

	item.item_type = NULL;
	
	for (j=0; j<neqns; j++) {
	    GtkWidget *w;

	    /* impulse responses: subitems for shocks */
	    vshock = gretl_VAR_get_variable_number(var, j);
	    item.callback_action = j;
	    double_underscores(tmp, datainfo->varname[vshock]);
	    sprintf(min, _("to %s"), tmp);

	    item.path = g_strdup_printf("%s/%s/%s", _(gpath), maj, min);
	    item.callback = impulse_plot_call;
	    item.callback_action = j;
	    item.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	    g_free(item.path);
	    w = gtk_item_factory_get_widget_by_action(vwin->ifac, j);
	    g_object_set_data(G_OBJECT(w), "targ", GINT_TO_POINTER(i));
	}
    }

    if (ci == VECM) {
	/* save ECs items */
	for (i=0; i<jrank(var); i++) {
	    item.path = g_strdup_printf("%s/%s %d", _(dpath), 
					_("EC term"), i+1);
	    item.callback = VECM_add_EC_data;
	    item.callback_action = i;
	    item.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
	    g_free(item.path);
	}
    }

    if (latex_is_ok()) {
	int n = sizeof sys_tex_items / sizeof sys_tex_items[0];

	for (i=0; i<n; i++) {
	    gtk_item_factory_create_item(vwin->ifac, &sys_tex_items[i], 
					 vwin, 1);
	}
    }
}

static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data)
{
    windata_t *mwin = (windata_t *) data;
    MODEL *pmod = mwin->data;
    gboolean s;
    int ok = 1, graphs_ok = 1;

    if (Z == NULL) {
	flip(mwin->ifac, "/File/Save to session as icon", FALSE);
	flip(mwin->ifac, "/File/Save as icon and close", FALSE);
	flip(mwin->ifac, "/Edit/Copy all", FALSE);
	flip(mwin->ifac, "/Analysis", FALSE);
	flip(mwin->ifac, "/Tests", FALSE);
	flip(mwin->ifac, "/Graphs", FALSE);
	flip(mwin->ifac, "/Analysis", FALSE);
	flip(mwin->ifac, "/LaTeX", FALSE);

	return FALSE;
    }

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->ci == MPOLS) {
	return FALSE;
    }

    if (model_sample_problem(pmod, datainfo)) {
	ok = 0;
	graphs_ok = (pmod->dataset != NULL);
	if (!graphs_ok && add_dataset_to_model(pmod, datainfo) == 0) {
	    graphs_ok = 1;
	}
    }

    s = GTK_WIDGET_IS_SENSITIVE(gtk_item_factory_get_item(mwin->ifac, "/Tests"));
    if ((s && ok) || (!s && !ok)) {
	/* no need to flip state */
	return FALSE;
    }

    flip(mwin->ifac, "/Tests", ok);
    flip(mwin->ifac, "/Graphs", graphs_ok);
    flip(mwin->ifac, "/Analysis/Display actual, fitted, residual", ok);
    flip(mwin->ifac, "/Analysis/Forecasts...", ok);
    flip(mwin->ifac, "/Analysis/Confidence intervals for coefficients", ok);
    flip(mwin->ifac, "/Save/Fitted values", ok);
    flip(mwin->ifac, "/Save/Residuals", ok);
    flip(mwin->ifac, "/Save/Squared residuals", ok);
    flip(mwin->ifac, "/Save/Define new variable...", ok);

    if (!ok) {
	const char *msg = gretl_errmsg_get();

	if (msg != NULL && *msg != 0) {
	    infobox(msg);
	}
    } 

    return FALSE;
}

int validate_varname (const char *varname)
{
    int i, n = strlen(varname);
    char namebit[VNAMELEN];
    unsigned char c;
    int err = 0;

    *namebit = 0;
    
    if (n > VNAMELEN - 1) {
	strncat(namebit, varname, VNAMELEN - 1);
	errbox(_("Variable name %s... is too long\n"
		 "(the max is %d characters)"), namebit,
	       VNAMELEN - 1);
	err = 1;
    } else if (!(isalpha(*varname))) {
	errbox(_("First char of name ('%c') is bad\n"
		 "(first must be alphabetical)"), *varname);
	err = 1;
    } else {
	for (i=1; i<n && !err; i++) {
	    c = (unsigned char) varname[i];
	
	    if ((!(isalpha(c)) && !(isdigit(c)) && c != '_') || c > 127) {
		errbox(_("Name contains an illegal char (in place %d)\n"
			 "Use only unaccented letters, digits and underscore"), i + 1);
		err = 1;
	    }
	}
    }

    return err;
}

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

void add_popup_item (const gchar *label, GtkWidget *menu,
		     GCallback callback, 
		     gpointer data)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(label);
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    g_signal_connect (G_OBJECT(item), "activate",
		      G_CALLBACK(callback), data);
    gtk_widget_show(item);
}

void *gui_get_plugin_function (const char *funcname, 
			       void **phandle)
{
    void *func;

    func = get_plugin_function(funcname, phandle);
    if (func == NULL) {
	errbox(gretl_errmsg_get());
    }

    return func;
}

char *double_underscores (char *targ, const char *src)
{
    char *p = targ;

    while (*src) {
	if (*src == '_') {
	    *p++ = '_';
	    *p++ = '_';
	} else {
	    *p++ = *src;
	}
	src++;
    }
    *p = '\0';

    return targ;
}

#ifndef G_OS_WIN32

int browser_open (const char *url)
{
# if defined(USE_GNOME)
    gnome_url_show(url, NULL); 
# elif defined(OSX_BUILD)
    osx_open_url(url);
# else
    gchar *urlcmd;
    int err;
    
    urlcmd = g_strdup_printf("%s -remote \"openURLNewWindow(%s)\"", Browser, url);
    err = gretl_spawn(urlcmd);
    g_free(urlcmd);

    if (err) {
	gretl_fork("Browser", url);
    }
# endif /* !GNOME, !OSX */

    return 0;
}

#include <signal.h>

void startR (const char *Rcommand)
{
    char Rprofile[MAXLEN], Rdata[MAXLEN], Rline[MAXLEN];
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    int *list;
    FILE *fp;
    int enverr;
    int i;
    char *s0, *s1, *s2;
    pid_t pid;

    if (!data_status) {
	warnbox(_("Please open a data file first"));
	return;
    }

    build_path(Rprofile, paths.dotdir, "gretl.Rprofile", NULL);
    fp = fopen(Rprofile, "w");
    if (fp == NULL) {
	file_write_errbox(Rprofile);
	return;
    }

    enverr = setenv("R_PROFILE", Rprofile, 1);
    if (enverr) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 	

    build_path(Rdata, paths.dotdir, "Rdata.tmp", NULL);

    sprintf(Rline, "store \"%s\" -r", Rdata);
    list = command_list_from_string(Rline);

    if (list == NULL ||
	write_data(Rdata, list, (const double **) Z, datainfo, 
		   OPT_R, NULL)) {
	errbox(_("Write of R data file failed"));
	fclose(fp);
	return; 
    }

    free(list);

    if (dataset_is_time_series(datainfo)) {
	fputs("# load data from gretl\n", fp);
	fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fp);
	fputs("if (vnum > 1.89) library(stats) else library(ts)\n", fp);
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rdata);
    } else {
	char Rtmp[MAXLEN];
	FILE *fq;

	build_path(Rtmp, paths.dotdir, "Rtmp", NULL);
	fq = fopen(Rtmp, "w");
	if (fq != NULL) {
	    fputs("# load data from gretl\n", fq);
	    fputs("library(stats)\n", fq);
	    fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fq);
	    fputs("if (vnum > 2.41) library(utils)\n", fq);
	    fprintf(fq, "gretldata <- read.table(\"%s\", header=TRUE)\n", Rdata);
	    fprintf(fq, "attach(gretldata)\n");
	    fclose(fq);
	}
	fprintf(fp, "source(\"%s\", echo=TRUE)\n", Rtmp);
    }

    fclose(fp);

    s0 = mymalloc(64);
    s1 = mymalloc(32);
    s2 = mymalloc(32);
    if (s0 == NULL || s1 == NULL || s2 == NULL) return;

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
	if (i == 1) {
	    execlp(s0, s0, supp1, supp2, NULL);
	} else if (i == 2) {
	    execlp(s0, s0, s1, supp1, supp2, NULL);
	} else if (i == 3) {
	    execlp(s0, s0, s1, s2, supp1, supp2, NULL);
	}
	perror("execlp");
	_exit(EXIT_FAILURE);
    }

    free(s0); 
    free(s1); 
    free(s2);
}

#endif /* ! G_OS_WIN32 */
