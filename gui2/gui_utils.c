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

#include "gretl.h"
#include "var.h"
#include "modelspec.h"
#include "forecast.h"

#include <sys/stat.h>
#include <unistd.h>

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

#ifdef G_OS_WIN32
# include <windows.h>
#endif

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
#endif

char *storelist = NULL;

#ifdef OLD_GTK
#include "../pixmaps/stock_save_16.xpm"
#include "../pixmaps/stock_save_as_16.xpm"
#include "../pixmaps/stock_exec_16.xpm"
#include "../pixmaps/stock_copy_16.xpm"
#include "../pixmaps/stock_paste_16.xpm"
#include "../pixmaps/stock_search_16.xpm"
#include "../pixmaps/stock_search_replace_16.xpm"
#include "../pixmaps/stock_undo_16.xpm"
#include "../pixmaps/stock_help_16.xpm"
#include "../pixmaps/stock_add_16.xpm"
#include "../pixmaps/stock_close_16.xpm"
# if defined(USE_GNOME)
#  include "../pixmaps/stock_print_16.xpm"
# endif
#endif

#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mail_send.xpm"

#define CONTENT_IS_CHANGED(w) (w->active_var == 1)

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static void file_viewer_save (GtkWidget *widget, windata_t *vwin);
static gint query_save_text (GtkWidget *w, GdkEvent *event, windata_t *vwin);
static void auto_save_script (windata_t *vwin);
static void add_model_dataset_items (windata_t *vwin);
static void add_model_tex_items (windata_t *vwin);
static void add_vars_to_plot_menu (windata_t *vwin);
static void add_dummies_to_plot_menu (windata_t *vwin);
static void add_VAR_menu_items (windata_t *vwin);
static void add_x12_output_menu_item (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);
static void buf_edit_save (GtkWidget *widget, gpointer data);
static void model_copy_callback (gpointer p, guint u, GtkWidget *w);
static void panel_heteroskedasticity_menu (windata_t *vwin);

#ifndef OLD_GTK
static int maybe_recode_file (const char *fname);
#endif

/* in editable window toolbar, Save is first button */

static GtkWidget *toolbar_first_button (GtkToolbar *tb)
{
    GList *kids = tb->children;
    GtkToolbarChild *child;
    GtkWidget *w = NULL;

    while (kids != NULL) {
	child = kids->data;
	if (child->type == GTK_TOOLBAR_CHILD_BUTTON) {
	    w = child->widget;
	    break;
	}
	kids = kids->next;
    }

    return w;
}

static void mark_content_changed (windata_t *vwin) 
{
    if (vwin->active_var == 0) {
	GtkWidget *w = toolbar_first_button(GTK_TOOLBAR(vwin->mbar));

	if (w != NULL) {
	    gtk_widget_set_sensitive(w, TRUE);
	}
	vwin->active_var = 1;	
    }
}

static void mark_content_saved (windata_t *vwin) 
{
    GtkWidget *w = toolbar_first_button(GTK_TOOLBAR(vwin->mbar));

    if (w != NULL) {
	gtk_widget_set_sensitive(w, FALSE);
    }
    vwin->active_var = 0;
}

static void close_model (gpointer data, guint close, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) data;

    if (!window_is_busy(vwin)) {
	gtk_widget_destroy(vwin->dialog);
    }
}

static int arma_by_x12a (const MODEL *pmod)
{
    int ret = 0;

    if (pmod->ci == ARMA && gretl_model_get_int(pmod, "arma_by_x12a")) {
	ret = 1;
    }

    return ret;
}

static int latex_is_ok (void)
{
    static int latex_ok = -1; 
  
    if (latex_ok == -1) {
#ifdef G_OS_WIN32
	latex_ok = check_for_prog("latex.exe");
#else
	latex_ok = check_for_prog("latex");
#endif
    }

    return latex_ok;
}

#ifndef OLD_GTK

static GtkItemFactoryEntry model_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, 
      "<StockItem>", GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as icon"), NULL, remember_model, 0, NULL, GNULL },
    { N_("/File/Save as icon and close"), NULL, remember_model, 1, NULL, GNULL },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
# endif
    { N_("/File/Close"), NULL, close_model, 0, "<StockItem>", GTK_STOCK_CLOSE },
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, "<StockItem>", GTK_STOCK_COPY },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/Tests/omit variables"), NULL, selector_callback, OMIT, NULL, GNULL },
    { N_("/Tests/add variables"), NULL, selector_callback, ADD, NULL, GNULL },
    { N_("/Tests/sum of coefficients"), NULL, selector_callback, COEFFSUM, NULL, GNULL },
    { N_("/Tests/linear restrictions"), NULL, gretl_callback, RESTRICT, NULL, GNULL },
    { N_("/Tests/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/non-linearity (squares)"), NULL, do_lmtest, LMTEST_SQUARES, NULL, GNULL },
    { N_("/Tests/non-linearity (logs)"), NULL, do_lmtest, LMTEST_LOGS, NULL, GNULL },
    { N_("/Tests/Ramsey's RESET"), NULL, do_reset, RESET, NULL, GNULL },
    { N_("/Tests/sep2"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/heteroskedasticity"), NULL, do_lmtest, LMTEST_WHITE, NULL, GNULL },
    { N_("/Tests/normality of residual"), NULL, do_resid_freq, TESTUHAT, NULL, GNULL },
    { N_("/Tests/influential observations"), NULL, do_leverage, LEVERAGE, NULL, GNULL },
    { N_("/Tests/collinearity"), NULL, do_vif, VIF, NULL, GNULL },
    { N_("/Tests/Chow test"), NULL, do_chow_cusum, CHOW, NULL, GNULL },
    { N_("/Tests/sep3"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/autocorrelation"), NULL, do_autocorr, LMTEST, NULL, GNULL },
    { N_("/Tests/ARCH"), NULL, do_arch, ARCH, NULL, GNULL },
    { N_("/Tests/CUSUM test"), NULL, do_chow_cusum, CUSUM, NULL, GNULL },
    { N_("/Tests/sep4"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Tests/panel diagnostics"), NULL, do_panel_diagnostics, HAUSMAN, NULL, GNULL },
    { N_("/_Graphs"), NULL, NULL, 0, "<Branch>", GNULL }, 
    { N_("/Graphs/residual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Graphs/fitted, actual plot"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/_Model data"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Model data/Display actual, fitted, residual"), NULL, 
      display_fit_resid, 0, NULL, GNULL },
    { N_("/Model data/Forecasts..."), NULL, 
      do_forecast, FCASTERR, NULL, GNULL },
    { N_("/Model data/Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL, GNULL },
    { N_("/Model data/coefficient covariance matrix"), NULL, 
      do_outcovmx, 0, NULL, GNULL },
    { N_("/Model data/Add to data set"), NULL, NULL, 0, "<Branch>", GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

#else /* now old versions */

static GtkItemFactoryEntry model_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, NULL },
    { N_("/File/Save to session as icon"), NULL, remember_model, 0, NULL },
    { N_("/File/Save as icon and close"), NULL, remember_model, 1, NULL },
# if defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
# endif
    { N_("/File/Close"), NULL, close_model, 0, NULL },

    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, NULL },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>" },    
    { N_("/Tests/omit variables"), NULL, selector_callback, OMIT, NULL },
    { N_("/Tests/add variables"), NULL, selector_callback, ADD, NULL },
    { N_("/Tests/sum of coefficients"), NULL, selector_callback, COEFFSUM, NULL },
    { N_("/Tests/linear restrictions"), NULL, gretl_callback, RESTRICT, NULL },
    { N_("/Tests/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/non-linearity (squares)"), NULL, do_lmtest, LMTEST_SQUARES, NULL },
    { N_("/Tests/non-linearity (logs)"), NULL, do_lmtest, LMTEST_LOGS, NULL },
    { N_("/Tests/Ramsey's RESET"), NULL, do_reset, RESET, NULL },
    { N_("/Tests/sep2"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/heteroskedasticity"), NULL, do_lmtest, LMTEST_WHITE, NULL },
    { N_("/Tests/normality of residual"), NULL, do_resid_freq, TESTUHAT, NULL },
    { N_("/Tests/influential observations"), NULL, do_leverage, LEVERAGE, NULL },
    { N_("/Tests/collinearity"), NULL, do_vif, VIF, NULL },
    { N_("/Tests/Chow test"), NULL, do_chow_cusum, CHOW, NULL },
    { N_("/Tests/sep3"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/autocorrelation"), NULL, do_autocorr, LMTEST, NULL },
    { N_("/Tests/ARCH"), NULL, do_arch, ARCH, NULL },
    { N_("/Tests/CUSUM test"), NULL, do_chow_cusum, CUSUM, NULL },
    { N_("/Tests/sep4"), NULL, NULL, 0, "<Separator>" },
    { N_("/Tests/panel diagnostics"), NULL, do_panel_diagnostics, HAUSMAN, NULL },
    { N_("/_Graphs"), NULL, NULL, 0, "<Branch>" }, 
    { N_("/Graphs/residual plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/Graphs/fitted, actual plot"), NULL, NULL, 0, "<Branch>" },
    { N_("/_Model data"), NULL, NULL, 0, "<Branch>" },
    { N_("/Model data/Display actual, fitted, residual"), NULL, 
      display_fit_resid, 0, NULL },
    { N_("/Model data/Forecasts..."), NULL, 
      do_forecast, FCASTERR, NULL },
    { N_("/Model data/Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL },
    { N_("/Model data/coefficient covariance matrix"), NULL, 
      do_outcovmx, 0, NULL },
    { N_("/Model data/Add to data set"), NULL, NULL, 0, "<Branch>" },
    { NULL, NULL, NULL, 0, NULL}
};
#endif /* old versus new GTK */

#ifndef OLD_GTK

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
    { N_("/LaTeX/Save/Tabular"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX, NULL, GNULL },
    { N_("/LaTeX/Save/Equation"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL, GNULL }
};

static GtkItemFactoryEntry VAR_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/_View"), NULL, var_tex_callback, 0, NULL, GNULL },
    { N_("/LaTeX/_Copy"), NULL, var_tex_callback, 1, NULL, GNULL },
    { N_("/LaTeX/_Save"), NULL, var_tex_callback, 2, NULL, GNULL }
};

static GtkItemFactoryEntry VAR_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, "<StockItem>", 
      GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as icon"), NULL, remember_var, 0, NULL, GNULL },
    { N_("/File/Save as icon and close"), NULL, remember_var, 1, NULL, GNULL },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
# endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, "<StockItem>", GTK_STOCK_COPY },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

#else

static GtkItemFactoryEntry model_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/_View"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/View/_Tabular"), NULL, model_tex_view, 
      GRETL_FORMAT_TEX, NULL },
    { N_("/LaTeX/View/_Equation"), NULL, model_tex_view, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL },
    { N_("/LaTeX/_Copy"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Copy/_Tabular"), NULL, window_copy, 
      GRETL_FORMAT_TEX, NULL },
    { N_("/LaTeX/Copy/_Equation"), NULL, window_copy, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL },
    { N_("/LaTeX/_Save"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/Save/Tabular"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX, NULL },
    { N_("/LaTeX/Save/Equation"), NULL, model_tex_save, 
      GRETL_FORMAT_TEX | GRETL_FORMAT_EQN, NULL }
};

static GtkItemFactoryEntry VAR_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>" },
    { N_("/LaTeX/_View"), NULL, var_tex_callback, 0, NULL },
    { N_("/LaTeX/_Copy"), NULL, var_tex_callback, 1, NULL },
    { N_("/LaTeX/_Save"), NULL, var_tex_callback, 2, NULL }
};

static GtkItemFactoryEntry VAR_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },
    { N_("/File/_Save as text..."), NULL, file_save, SAVE_MODEL, NULL },
    { N_("/File/Save to session as icon"), NULL, remember_var, 0, NULL },
    { N_("/File/Save as icon and close"), NULL, remember_var, 1, NULL },
# if defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
# endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, NULL },
    { NULL, NULL, NULL, 0, NULL}
};

#endif /* old versus new GTK */

static void model_copy_callback (gpointer p, guint u, GtkWidget *w)
{
    copy_format_dialog((windata_t *) p, 1);
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
	sprintf(errtext, _("Couldn't open %s"), src);
	errbox(errtext);
	return 1; 
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
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

int isdir (const char *path)
{
    struct stat buf;

    return (stat(path, &buf) == 0 && S_ISDIR(buf.st_mode)); 
}

/* ........................................................... */

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
	    } else if (role == VAR) {
		const GRETL_VAR *var;

		var = g_object_get_data(G_OBJECT(wstack[i]), "object");
		if (var != NULL) {
		    mvm = gretl_VAR_get_highest_variable(var, datainfo);
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
    STACK_MAXVAR
};

static int winstack (int code, GtkWidget *w, gpointer ptest)
{
    static int n_windows;
    static GtkWidget **wstack;
    int i, ret = 0;

    switch (code) {

    case STACK_DESTROY:	
	for (i=0; i<n_windows; i++) {
	    if (wstack[i] != NULL) {
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
	    n_windows++;
	    wstack = myrealloc(wstack, n_windows * sizeof *wstack);
	    if (wstack != NULL) { 
		wstack[n_windows-1] = w;
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
		    ret = 1;
		    break;
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

int highest_numbered_variable_in_winstack (void)
{
    return winstack(STACK_MAXVAR, NULL, NULL);
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

static void delete_file_viewer (GtkWidget *widget, gpointer data) 
{
    windata_t *vwin = (windata_t *) data;
    gint resp = 0;

    if (window_is_busy(vwin)) {
	return;
    }

    if ((vwin->role == EDIT_SCRIPT || vwin->role == EDIT_HEADER ||
	 vwin->role == EDIT_NOTES) && CONTENT_IS_CHANGED(vwin)) {
	resp = query_save_text(NULL, NULL, vwin);
    }

    if (!resp) {
	gtk_widget_destroy(vwin->dialog); 
    }
}

static void delete_unnamed_model (GtkWidget *widget, gpointer data) 
{
    MODEL *pmod = (MODEL *) data;

    if (pmod->dataset != NULL) {
	free_model_dataset(pmod);
    }

    if (pmod->name == NULL) {
	free_model(pmod);
    }
}

void delete_widget (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
}

#ifndef OLD_GTK

static gint catch_button_3 (GtkWidget *w, GdkEventButton *event)
{
    GdkModifierType mods;

    gdk_window_get_pointer (w->window, NULL, NULL, &mods); 
    if (mods & GDK_BUTTON3_MASK) {
	return TRUE;
    }
    return FALSE;
}

#endif

#ifdef G_OS_WIN32

static void win_ctrl_c (windata_t *vwin)
{
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    if (gtk_text_buffer_get_selection_bounds(buf, NULL, NULL)) {
	window_copy(vwin, GRETL_FORMAT_SELECTION, NULL);
    } else if (MULTI_FORMAT_ENABLED(vwin->role)) {
	window_copy(vwin, GRETL_FORMAT_RTF, NULL);
    } else {
	window_copy(vwin, GRETL_FORMAT_TXT, NULL);
    }
}

#endif

static gint catch_edit_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods);

    if (key->keyval == GDK_F1 && vwin->role == EDIT_SCRIPT) { 
	set_window_help_active(vwin);
	edit_script_help(NULL, NULL, vwin);
    }

#if !defined(OLD_GTK) && !defined(USE_GTKSOURCEVIEW)
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
	} else if (gdk_keyval_to_upper(key->keyval) == GDK_Q) {
	    if (vwin->role == EDIT_SCRIPT && CONTENT_IS_CHANGED(vwin)) {
		gint resp;

		resp = query_save_text(NULL, NULL, vwin);
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

void audio_render_window (windata_t *vwin, int key)
{
    void *handle;
    int (*read_window_text) (windata_t *, const DATAINFO *, int (*)());

    if (vwin == NULL) {
	set_or_get_audio_stop(1, 1);
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

static gint catch_viewer_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
#ifndef OLD_GTK
    if (gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->w))) {
	return catch_edit_key(w, key, vwin);
    }
#else
    if (GTK_EDITABLE(vwin->w)->editable) {
	return catch_edit_key(w, key, vwin);
    }    
#endif

    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    } else if (key->keyval == GDK_s && Z != NULL && vwin->role == VIEW_MODEL) {
	remember_model(vwin, 1, NULL);
    } else if (key->keyval == GDK_w) {
	GdkModifierType mods;

	gdk_window_get_pointer(w->window, NULL, NULL, &mods); 
	if (mods & GDK_CONTROL_MASK) {
	    gtk_widget_destroy(w);
	    return TRUE;
	}	
    }
#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    else if (key->keyval == GDK_a) {
	audio_render_window(vwin, AUDIO_TEXT);
    } else if (key->keyval == GDK_x) {
	audio_render_window(NULL, AUDIO_TEXT);
    }
#endif
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

    /* Ctrl-F for find */
    if (key->keyval == GDK_f) {
	GdkModifierType mods;

	gdk_window_get_pointer(w->window, NULL, NULL, &mods); 
	if (mods & GDK_CONTROL_MASK) {
	    text_find_callback(NULL, vwin);
	    return TRUE;
	}
    }

    return FALSE;
}

gint catch_listbox_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    if (key->keyval == GDK_q) { 
	if (vwin != mdata) {
	    gtk_widget_destroy(vwin->w);
	}
	return TRUE;
    } else if (key->keyval == GDK_f) {
	GdkModifierType mods;

	if (vwin == mdata && !data_status) {
	    return TRUE;
	}

	gdk_window_get_pointer(w->window, NULL, NULL, &mods); 
	if (mods & GDK_CONTROL_MASK) {
	    menu_find(vwin, 1, NULL);
	    return TRUE;
	}	
    }

    return FALSE;
}

/* ........................................................... */

void *mymalloc (size_t size) 
{
    void *mem;
   
    if ((mem = malloc(size)) == NULL) 
	errbox(_("Out of memory!"));
    return mem;
}

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
	if (strstr(paths.datfile, paths.datadir) != NULL) {
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

    /* focus the data window */
    gtk_widget_grab_focus(mdata->listbox);

    /* invalidate "remove extra obs" menu item */
    drop_obs_state(FALSE);
}

#define APPENDING(action) (action == APPEND_DATA || \
                           action == APPEND_CSV || \
                           action == APPEND_GNUMERIC || \
                           action == APPEND_EXCEL || \
                           action == APPEND_ASCII)

int get_worksheet_data (char *fname, int datatype, int append,
			int *gui_get_data)
{
    void *handle;
    PRN *errprn;
    const char *errbuf;
    FILE *fp;
    int (*sheet_get_data)(const char*, double ***, DATAINFO *, PRN *);
    int err = 0;
    
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	sprintf(errtext, _("Couldn't open %s"), fname);
	errbox(errtext);
	return 1;
    } else {
	fclose(fp);
    }

    if (datatype == GRETL_GNUMERIC) {
	sheet_get_data = gui_get_plugin_function("wbook_get_data",
						 &handle);
    } else if (datatype == GRETL_EXCEL) {
	sheet_get_data = gui_get_plugin_function("excel_get_data",
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

    err = (*sheet_get_data)(fname, &Z, datainfo, errprn);
    close_plugin(handle);

    if (err == -1) {
	fprintf(stderr, "data import canceled\n");
	if (gui_get_data != NULL) {
	    *gui_get_data = 1;
	}
	gretl_print_destroy(errprn);
	return 0;
    }

    errbuf = gretl_print_get_buffer(errprn);

    if (err) {
	if (errbuf != NULL && *errbuf != '\0') {
	    errbox(errbuf);
	} else {
	    errbox(_("Failed to import spreadsheet data"));
	}
	gretl_print_destroy(errprn);
	return 1;
    } else {
	if (errbuf != NULL && *errbuf != '\0') {
	    infobox(errbuf);
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
    }

    return err;
}

static void copy_utf8_filename (char *targ, const char *src)
{
    strcpy(targ, src);
#if defined(ENABLE_NLS) && !defined(OLD_GTK)
    my_filename_to_utf8(targ);
#endif
}

/* cases for do_open_data: 
   - called from dialog: user has said Yes to opening data file,
     although a data file is already open (or user wants to append
     data)
   - reached without dialog, in expert mode or when no datafile
     is open yet
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
    } else if (code == OPEN_EXCEL || code == APPEND_EXCEL) {
	datatype = GRETL_EXCEL;
    } else if (code == OPEN_OCTAVE || code == APPEND_OCTAVE) {
	datatype = GRETL_OCTAVE;
    } else if (code == OPEN_BOX) {
	datatype = GRETL_BOX_DATA;
    } else {
	/* no filetype specified: have to guess */
	PRN *prn;	

	if (bufopen(&prn)) return;
	datatype = detect_filetype(trydatfile, &paths, prn);
	gretl_print_destroy(prn);
    }

    /* destroy the current data set, etc., unless we're explicitly appending */
    if (!append) {
	close_session();
    }

    if (datatype == GRETL_GNUMERIC || datatype == GRETL_EXCEL) {
	get_worksheet_data(trydatfile, datatype, append, NULL);
	return;
    } else if (datatype == GRETL_CSV_DATA) {
	do_open_csv_box(trydatfile, OPEN_CSV, append);
	return;
    } else if (datatype == GRETL_OCTAVE) {
	do_open_csv_box(trydatfile, OPEN_OCTAVE, append);
	return;
    } else if (datatype == GRETL_BOX_DATA) {
	do_open_csv_box(trydatfile, OPEN_BOX, 0);
	return;
    } else { /* native data */
	int clear_code = DATA_NONE;
	PRN *errprn;

	errprn = gretl_print_new(GRETL_PRINT_STDERR);

	if (append) {
	    clear_code = DATA_APPEND;
	} else if (data_status) {
	    clear_code = DATA_CLEAR;
	}

	if (datatype == GRETL_XML_DATA) {
	    err = get_xmldata(&Z, &datainfo, trydatfile, &paths, 
			      clear_code, errprn, 1);
	} else {
	    err = gretl_get_data(&Z, &datainfo, trydatfile, &paths, 
				 clear_code, errprn);
	}

	gretl_print_destroy(errprn);
    }

    if (err) {
	gui_errmsg(err);
	delete_from_filelist(FILE_LIST_DATA, trydatfile);
	return;
    }	

    /* trash the practice files window that launched the query? */
    if (fwin != NULL) {
	gtk_widget_destroy(fwin->w);
    }

    if (append) {
	register_data(NULL, NULL, 0);
    } else {
	copy_utf8_filename(paths.datfile, trydatfile);
	register_data(paths.datfile, NULL, 1);
    } 
}

/* give user choice of not opening selected datafile, if there's
   already a datafile open and we're not in "expert" mode */

void verify_open_data (gpointer userdata, int code)
{
    if (dataset_locked()) {
	return;
    }

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

/* give user choice of not opening session file, if there's already a
   datafile open and we're not in "expert" mode */

void verify_open_session (gpointer userdata)
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

void save_session (char *fname) 
{
    char msg[MAXLEN], savedir[MAXLEN], fname2[MAXLEN];
    char session_base[MAXLEN];
    int spos;
    FILE *fp;
    PRN *prn;

    *savedir = '\0';

    spos = slashpos(fname);
    if (spos) {
	safecpy(savedir, fname, spos);
    } 

#ifdef CMD_DEBUG
    dump_command_stack("stderr", 0);
#endif

    /* append ".gretl" to session filename? */
    if (haschar('.', fname) < 0) {
	strcat(fname, ".gretl");
    }

    /* save commands, by dumping the command stack */
    if (dump_command_stack(fname, 1)) {
	return;
    }

    get_base(session_base, fname, '.');

    /* get ready to save "session" */
    fp = gretl_fopen(fname, "a");
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
    prn = gretl_print_new_with_filename(fname2);
    if (prn == NULL) {
	errbox(_("Couldn't open output file for writing"));
	return;
    }

    /* preamble */
    gui_logo(prn);
    session_time(prn);
    pprintf(prn, _("Output from %s\n"), fname);

    /* actual commands output */
    execute_script(fname, NULL, prn, SAVE_SESSION_EXEC); 

    gretl_print_destroy(prn);

    /* output may need re-encoding, UTF-8 to locale? */
#ifndef OLD_GTK
    maybe_recode_file(fname2);
#endif
    
    sprintf(msg, _("session saved to %s -\n"), savedir);
    strcat(msg, _("commands: "));
    strcat(msg, (spos)? fname + spos + 1 : fname);
    strcat(msg, _("\noutput: "));
    spos = slashpos(fname2);
    strcat(msg, (spos)? fname2 + spos + 1 : fname2);
    infobox(msg);

    mkfilelist(FILE_LIST_SESSION, fname);
    set_session_saved(1);
    session_changed(0);

    return;
}

/* ........................................................... */

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
#ifndef OLD_GTK
    text_set_cursor(vwin->w, GDK_QUESTION_ARROW);
#else
    GdkCursor *cursor = gdk_cursor_new(GDK_QUESTION_ARROW);

    gdk_window_set_cursor(GTK_TEXT(vwin->w)->text_area, cursor);
    gdk_cursor_destroy(cursor);
#endif

    set_window_help_active(vwin);
}

/* ........................................................... */

static void buf_edit_save (GtkWidget *widget, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *text;
    char **pbuf = (char **) vwin->data;

#ifndef OLD_GTK
    text = textview_get_text(GTK_TEXT_VIEW(vwin->w));
#else
    text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
#endif

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
	session_changed(1);
    }
}

static void file_viewer_save (GtkWidget *widget, windata_t *vwin)
{
    if (strstr(vwin->fname, "script_tmp") || !strlen(vwin->fname)) {
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
#ifndef OLD_GTK
	    text = textview_get_text(GTK_TEXT_VIEW(vwin->w));
#else
	    text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
#endif
	    system_print_buf(text, fp);
	    fclose(fp);
	    g_free(text);
	    mark_content_saved(vwin);
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
#ifdef USE_GTKSOURCEVIEW
    vwin->sbuf = NULL;
#endif
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

void free_windata (GtkWidget *w, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (vwin != NULL) {
	if (vwin->w != NULL) { 
	    gchar *undo = g_object_get_data(G_OBJECT(vwin->w), "undo");
	    
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
#ifndef OLD_GTK 
	    gtk_widget_destroy(GTK_WIDGET(vwin->popup));
#else
	    gtk_object_unref(GTK_OBJECT(vwin->popup));
#endif
	}
	if (vwin->ifac != NULL) {
#ifndef OLD_GTK 
	    g_object_unref(G_OBJECT(vwin->ifac));
#else
	    gtk_object_unref(GTK_OBJECT(vwin->ifac));
#endif
	}

	/* data specific to certain windows */
	if (vwin->role == SUMMARY || vwin->role == VAR_SUMMARY) {
	    free_summary(vwin->data); 
	} else if (vwin->role == CORR || vwin->role == PCA || 
		   vwin->role == COVAR) {
	    free_vmatrix(vwin->data);
	} else if (vwin->role == FCASTERR || vwin->role == FCAST) {
	    free_fit_resid(vwin->data);
	} else if (vwin->role == COEFFINT) {
	    free_coeff_intervals(vwin->data);
	} else if (vwin->role == MPOLS) {
	    free_gretl_mp_results(vwin->data);
	} else if (vwin->role == VIEW_SERIES) {
	    free_series_view(vwin->data);
	} else if (vwin->role == VAR) { 
	    gretl_VAR_free_unnamed(vwin->data);
	} else if (vwin->role == LEVERAGE) {
	    gretl_matrix_free(vwin->data);
	} else if (vwin->role == MAHAL) {
	    free_mahal_dist(vwin->data);
	}

	if (vwin->dialog) {
	    winstack_remove(vwin->dialog);
	}

	free(vwin);
    }
}

#ifndef OLD_GTK
void gretl_stock_icons_init (void)
{
    static GtkIconFactory *ifac;

    if (ifac == NULL) {
	GtkIconSource *source;
	GtkIconSet *set;
	GdkPixbuf *pbuf;

	ifac = gtk_icon_factory_new();

	set = gtk_icon_set_new();
	source = gtk_icon_source_new();
	gtk_icon_source_set_size(source, GTK_ICON_SIZE_SMALL_TOOLBAR);
	pbuf = gdk_pixbuf_new_from_xpm_data((const char **) mini_tex_xpm);
	gtk_icon_source_set_pixbuf(source, pbuf);
	g_object_unref(pbuf);
	gtk_icon_set_add_source(set, source);
	gtk_icon_source_free(source);
	gtk_icon_factory_add(ifac, GRETL_STOCK_TEX, set);
	gtk_icon_set_unref(set);

	set = gtk_icon_set_new();

	source = gtk_icon_source_new();
	gtk_icon_source_set_size(source, GTK_ICON_SIZE_SMALL_TOOLBAR);
	pbuf = gdk_pixbuf_new_from_xpm_data((const char **) mini_mail_xpm);
	gtk_icon_source_set_pixbuf(source, pbuf);
	g_object_unref(pbuf);
	gtk_icon_set_add_source(set, source);
	gtk_icon_source_free(source);

	source = gtk_icon_source_new();
	gtk_icon_source_set_size(source, GTK_ICON_SIZE_MENU);
	pbuf = gdk_pixbuf_new_from_xpm_data((const char **) mail_send_xpm);
	gtk_icon_source_set_pixbuf(source, pbuf);
	g_object_unref(pbuf);
	gtk_icon_set_add_source(set, source);
	gtk_icon_source_free(source);

	gtk_icon_factory_add(ifac, GRETL_STOCK_MAIL, set);
	gtk_icon_set_unref(set);

	gtk_icon_factory_add_default(ifac);
    }
}
#endif

#if defined(G_OS_WIN32) || defined(USE_GNOME) 
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(vwin, 0, w);
}
#endif

static void choose_copy_format_callback (GtkWidget *w, windata_t *vwin)
{
    copy_format_dialog(vwin, MULTI_FORMAT_ENABLED(vwin->role));
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
    } else if (vwin->role == FCASTERR) {
	add_fcast_data(vwin);
    }

    if (datainfo->v > oldv) {
	populate_varlist();
	mark_dataset_as_modified();
    }	
}

static void window_help (GtkWidget *w, windata_t *vwin)
{
    context_help(NULL, GINT_TO_POINTER(vwin->role));
}

struct viewbar_item {
    const char *str;
#ifndef OLD_GTK
    const gchar *icon;
#else
    gchar **toolxpm;
#endif
    void (*toolfunc)();
    int flag;
};

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    GP_ITEM,
    RUN_ITEM,
    COPY_ITEM,
    TEX_ITEM,
    ADD_ITEM,
    HELP_ITEM
} viewbar_codes;

#ifndef OLD_GTK

static struct viewbar_item viewbar_items[] = {
    { N_("Save"), GTK_STOCK_SAVE, file_viewer_save, SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, file_save_callback, SAVE_AS_ITEM },
    { N_("Send to gnuplot"), GTK_STOCK_EXECUTE, gp_send_callback, GP_ITEM },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("Print..."), GTK_STOCK_PRINT, window_print_callback, 0 },
# endif
    { N_("Run"), GTK_STOCK_EXECUTE, run_script_callback, RUN_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, text_copy_callback, COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, text_paste_callback, EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, text_find_callback, 0 },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, text_replace_callback, EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, text_undo_callback, EDIT_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, activate_script_help, RUN_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, window_tex_callback, TEX_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, add_data_callback, ADD_ITEM },
    { N_("Help"), GTK_STOCK_HELP, window_help, HELP_ITEM },
    { N_("Close"), GTK_STOCK_CLOSE, delete_file_viewer, 0 },
    { NULL, NULL, NULL, 0 }};

#else

static struct viewbar_item viewbar_items[] = {
    { N_("Save"), stock_save_16_xpm, file_viewer_save, SAVE_ITEM },
    { N_("Save as..."), stock_save_as_16_xpm, file_save_callback, SAVE_AS_ITEM },
    { N_("Send to gnuplot"), stock_exec_16_xpm, gp_send_callback, GP_ITEM },
# ifdef USE_GNOME
    { N_("Print..."), stock_print_16_xpm, window_print_callback, 0 },
# endif
    { N_("Run"), stock_exec_16_xpm, run_script_callback, RUN_ITEM },
    { N_("Copy"), stock_copy_16_xpm, text_copy_callback, COPY_ITEM }, 
    { N_("Paste"), stock_paste_16_xpm, text_paste_callback, EDIT_ITEM },
    { N_("Find..."), stock_search_16_xpm, text_find_callback, 0 },
    { N_("Replace..."), stock_search_replace_16_xpm, text_replace_callback, EDIT_ITEM },
    { N_("Undo"), stock_undo_16_xpm, text_undo_callback, EDIT_ITEM },
    { N_("Help on command"), stock_help_16_xpm, activate_script_help, RUN_ITEM },
    { N_("LaTeX"), mini_tex_xpm, window_tex_callback, TEX_ITEM },
    { N_("Add to dataset..."), stock_add_16_xpm, add_data_callback, ADD_ITEM },
    { N_("Help"), stock_help_16_xpm, window_help, HELP_ITEM },
    { N_("Close"), stock_close_16_xpm, delete_file_viewer, 0 },
    { NULL, NULL, NULL, 0 }};

#endif /* old versus new GTK */

#define editor_role(r) (r == EDIT_SCRIPT || \
                        r == EDIT_HEADER || \
                        r == EDIT_NOTES || \
                        r == GR_PLOT)

static void make_viewbar (windata_t *vwin, int text_out)
{
    GtkWidget *hbox, *button;
#ifdef OLD_GTK
    GdkPixmap *icon;
    GdkBitmap *mask;
    GdkColormap *cmap;
#endif
    void (*toolfunc)() = NULL;
    int i;

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

    int help_ok = (vwin->role == LEVERAGE || 
		   vwin->role == COINT2 ||
		   vwin->role == MAHAL);

    int latex_ok = latex_is_ok();

#ifndef OLD_GTK
    if (MULTI_FORMAT_ENABLED(vwin->role) && latex_ok) {
	gretl_stock_icons_init();
    }
#endif

    if (text_out || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->dialog), "text_out", GINT_TO_POINTER(1));
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

#ifndef OLD_GTK
    vwin->mbar = gtk_toolbar_new();
#else
    vwin->mbar = gtk_toolbar_new(GTK_ORIENTATION_HORIZONTAL, GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_button_relief(GTK_TOOLBAR(vwin->mbar), GTK_RELIEF_NONE);
    gtk_toolbar_set_space_size(GTK_TOOLBAR(vwin->mbar), 3);
#endif
    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);

#ifdef OLD_GTK
    cmap = gdk_colormap_get_system();
    colorize_tooltips(GTK_TOOLBAR(vwin->mbar)->tooltips);
#endif

    for (i=0; viewbar_items[i].str != NULL; i++) {
	GtkWidget *w;

	toolfunc = viewbar_items[i].toolfunc;

	if (!edit_ok && viewbar_items[i].flag == EDIT_ITEM) {
	    continue;
	}

	if (!save_as_ok && viewbar_items[i].flag == SAVE_AS_ITEM) {
	    continue;
	}

	if (!run_ok && viewbar_items[i].flag == RUN_ITEM) {
	    continue;
	}

	if (!help_ok && viewbar_items[i].flag == HELP_ITEM) {
	    continue;
	}

	if (vwin->role != GR_PLOT && viewbar_items[i].flag == GP_ITEM) {
	    continue;
	}

	if ((!latex_ok || !MULTI_FORMAT_ENABLED(vwin->role)) && 
	    viewbar_items[i].flag == TEX_ITEM) {
	    continue;
	}

	if (vwin->role != PCA && vwin->role != LEVERAGE && 
	    vwin->role != MAHAL && vwin->role != FCASTERR &&
	    viewbar_items[i].flag == ADD_ITEM) {
	    continue;
	}

#ifndef OLD_GTK
	if (viewbar_items[i].flag == COPY_ITEM && 
	    !editor_role(vwin->role)) {
	    toolfunc = choose_copy_format_callback;
	}
#else
        if (viewbar_items[i].flag == COPY_ITEM && 
            MULTI_FORMAT_ENABLED(vwin->role)) {
            toolfunc = choose_copy_format_callback;
        }
#endif

	if (viewbar_items[i].flag == SAVE_ITEM) { 
	    if (!edit_ok || vwin->role == SCRIPT_OUT) {
		/* script output doesn't already have a filename */
		continue;
	    }
	    if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		toolfunc = buf_edit_save;
	    } else if (vwin->role == GR_PLOT) {
		toolfunc = save_plot_commands_callback;
	    }
	}

#ifndef OLD_GTK
	button = gtk_image_new();
	gtk_image_set_from_stock(GTK_IMAGE(button), viewbar_items[i].icon, 
				 GTK_ICON_SIZE_MENU);
        w = gtk_toolbar_append_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(viewbar_items[i].str), NULL,
				    button, toolfunc, vwin);
#else
	icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, NULL, 
						     viewbar_items[i].toolxpm);
	button = gtk_pixmap_new(icon, mask);
	w = gtk_toolbar_append_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(viewbar_items[i].str), NULL,
				    button, toolfunc, vwin);
	gtk_toolbar_append_space(GTK_TOOLBAR(vwin->mbar));
#endif

	if (viewbar_items[i].flag == SAVE_ITEM) { 
	    /* nothing to save just yet */
	    gtk_widget_set_sensitive(w, FALSE);
	}
    }

    gtk_widget_show(vwin->mbar);
    gtk_widget_show(hbox);
}

static void add_edit_items_to_viewbar (windata_t *vwin)
{
    GtkWidget *button;
#ifdef OLD_GTK
    GdkPixmap *icon;
    GdkBitmap *mask;
    GdkColormap *cmap = gdk_colormap_get_system();
#endif
    int i, pos = 0;

    for (i=0; viewbar_items[i].str != NULL; i++) {
	if (viewbar_items[i].flag == SAVE_ITEM ||
	    viewbar_items[i].flag == EDIT_ITEM) {

#ifndef OLD_GTK
	    button = gtk_image_new();
	    gtk_image_set_from_stock(GTK_IMAGE(button), 
				     viewbar_items[i].icon, 
				     GTK_ICON_SIZE_MENU);
	    gtk_toolbar_insert_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(viewbar_items[i].str), NULL,
				    button, viewbar_items[i].toolfunc, 
				    vwin, pos);
#else
	    icon = gdk_pixmap_colormap_create_from_xpm_d(NULL, cmap, &mask, NULL, 
							 viewbar_items[i].toolxpm);
	    button = gtk_pixmap_new(icon, mask);
	    gtk_toolbar_insert_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(viewbar_items[i].str), NULL,
				    button, viewbar_items[i].toolfunc, 
				    vwin, pos);
	    gtk_toolbar_insert_space(GTK_TOOLBAR(vwin->mbar), pos + 1);
#endif
	}
	if (viewbar_items[i].flag != GP_ITEM) {
#ifndef OLD_GTK	    
	    pos++;
#else
	    pos += 2;
#endif
	}
    }
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
	if (strstr(fname, "script_tmp") || strstr(fname, "session.inp")) {
	    title = g_strdup(_("gretl: command script"));
	} else {
	    const char *p = strrchr(fname, SLASH);

	    title = g_strdup_printf("gretl: %s", 
				    (p != NULL)? p + 1 : fname);
#if defined(ENABLE_NLS) && !defined(OLD_GTK)
	    my_filename_to_utf8(title);
#endif	    
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

static void content_changed (GtkWidget *w, windata_t *vwin)
{
    mark_content_changed(vwin);
}

/* ........................................................... */

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

/* ........................................................... */

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);

    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);

#ifndef OLD_GTK
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);
# ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(vwin->dialog), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
# endif
#else
    gtk_container_border_width(GTK_CONTAINER(vwin->vbox), 4);
    gtk_signal_connect_after(GTK_OBJECT(vwin->dialog), "realize", 
                             GTK_SIGNAL_FUNC(set_wm_icon), 
                             NULL);
#endif

    gtk_container_add(GTK_CONTAINER(vwin->dialog), vwin->vbox);
}

static void view_buffer_insert_text (windata_t *vwin, PRN *prn)
{
    const char *pbuf = gretl_print_get_buffer(prn);
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
    if (vwin->role == SCRIPT_OUT) {
	text_buffer_insert_colorized_buffer(tbuf, prn);
    } else {
	gtk_text_buffer_set_text(tbuf, pbuf, -1);
    }
#else
    if (vwin->role == SCRIPT_OUT) {
	text_buffer_insert_colorized_buffer(vwin->w, prn);
    } else {
	gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
			NULL, NULL, pbuf, strlen(pbuf));
    }
#endif
}

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			const char *title, int role, 
			gpointer data) 
{
    GtkWidget *close;
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;
#endif
    windata_t *vwin;

    vwin = common_viewer_new(role, 
			     (title != NULL)? title : make_viewer_title(role, NULL), 
			     data, 
			     1);
    if (vwin == NULL) return NULL;

#ifdef OLD_GTK
    create_text(vwin, hsize, vsize, FALSE);
#endif

    viewer_box_config(vwin);

    /* in a few special cases, add a text-based menu bar */
    if (role == VAR || role == VIEW_SERIES || role == VIEW_SCALAR) {
	GtkItemFactoryEntry *menu_items;

	if (role == VAR) {
	    menu_items = VAR_items;
	} else {
	    menu_items = get_series_view_menu_items(role);
	}
	set_up_viewer_menu(vwin->dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    } else if (role != IMPORT) {
	make_viewbar(vwin, 1);
    }

    if (role == VAR) {
	/* model-specific additions to menus */
	add_VAR_menu_items(vwin);
    }

#ifndef OLD_GTK
    create_text(vwin, &tbuf, hsize, vsize, FALSE);
#endif
    
    text_table_setup(vwin);

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
    view_buffer_insert_text(vwin, prn);
    gretl_print_destroy(prn);

    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

#ifndef OLD_GTK    
    g_signal_connect (G_OBJECT(vwin->w), "button_press_event", 
		      G_CALLBACK(catch_button_3), vwin->w);
#endif

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

#ifndef OLD_GTK    
    cursor_to_top(vwin);
#endif

    return vwin;
}

#ifdef OLD_GTK
static void set_file_view_style (GtkWidget *w)
{
    static GtkStyle *style;

    if (style == NULL) {
	style = gtk_style_new();
	gdk_font_unref(style->font);
	style->font = fixed_font;
    }
    gtk_widget_set_style(w, style);
}
#endif

/* ........................................................... */

windata_t *view_file (const char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role)
{
    GtkWidget *close;
#ifndef OLD_GTK
    GtkTextBuffer *tbuf = NULL;
# ifdef USE_GTKSOURCEVIEW
    GtkSourceBuffer *sbuf = NULL;
# endif
#endif
    FILE *fp;
    windata_t *vwin;
    gchar *title = NULL;
    int doing_script = (role == EDIT_SCRIPT ||
			role == VIEW_SCRIPT ||
			role == VIEW_LOG);

    /* first check that we can open the specified file */
    fp = gretl_fopen(filename, "r");
    if (fp == NULL) {
	sprintf(errtext, _("Can't open %s for reading"), filename);
	errbox(errtext);
	return NULL;
    } else {
	fclose(fp);
    }

    /* then start building the file viewer */
    title = make_viewer_title(role, filename);
    vwin = common_viewer_new(role, (title != NULL)? title : filename, 
			     NULL, !doing_script && role != CONSOLE);
    g_free(title);
    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);

#ifdef OLD_GTK
    create_text(vwin, hsize, vsize, editable);
#endif

    viewer_box_config(vwin);

    if (help_role(role)) {
	GtkItemFactoryEntry *menu_items;

	menu_items = get_help_menu_items(role);
	set_up_viewer_menu(vwin->dialog, vwin, menu_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), 
			   vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
    } else { /* was else if (role != VIEW_FILE) */
	make_viewbar(vwin, (role == VIEW_DATA || role == CONSOLE));
    }

#ifndef OLD_GTK
# ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	create_source(vwin, &sbuf, hsize, vsize, editable);
	tbuf = GTK_TEXT_BUFFER(sbuf);
	vwin->sbuf = sbuf;
    } else {
	create_text(vwin, &tbuf, hsize, vsize, editable);
    }
# else
    create_text(vwin, &tbuf, hsize, vsize, editable);
# endif
#endif

    text_table_setup(vwin);

#ifdef OLD_GTK
    set_file_view_style(GTK_WIDGET(vwin->w));
#endif

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

    /* make a Close button */
    close = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), 
		       close, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(close), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(close);

#ifndef OLD_GTK
# ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	source_buffer_insert_file(sbuf, filename, role);
    } else {
	text_buffer_insert_file(tbuf, filename, role);
    }
# else
    text_buffer_insert_file(tbuf, filename, role);
# endif
#else
    text_buffer_insert_file(vwin->w, filename, role);
#endif

#ifndef OLD_GTK
    g_object_set_data(G_OBJECT(vwin->w), "tbuf", tbuf);
#endif

    /* grab the "changed" signal when editing a script */
    if (role == EDIT_SCRIPT) {
#ifndef OLD_GTK
        g_signal_connect(G_OBJECT(tbuf), "changed", 
                         G_CALLBACK(content_changed), vwin);
#else
        gtk_signal_connect(GTK_OBJECT(vwin->w), "changed", 
                           GTK_SIGNAL_FUNC(content_changed), vwin);
#endif
    }

    /* catch some keystrokes */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (editable) {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
    }

#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
#endif

    /* offer chance to save script on exit */
    if (role == EDIT_SCRIPT) {
	g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
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

#ifndef OLD_GTK
    cursor_to_top(vwin);
#endif

    gtk_widget_grab_focus(vwin->w);

    return vwin;
}

void file_view_set_editable (windata_t *vwin)
{
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;

    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), TRUE);
    g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);

    tbuf = GTK_TEXT_BUFFER(g_object_get_data(G_OBJECT(vwin->w), "tbuf"));
    g_signal_connect(G_OBJECT(tbuf), "changed", 
		     G_CALLBACK(content_changed), vwin);
#else
    gtk_text_set_editable(GTK_TEXT(vwin->w), TRUE);
    gtk_object_set_data(GTK_OBJECT(vwin->dialog), "vwin", vwin);
    gtk_signal_connect(GTK_OBJECT(vwin->w), "changed", 
		       GTK_SIGNAL_FUNC(content_changed), vwin);
#endif

    vwin->role = EDIT_SCRIPT;
    add_edit_items_to_viewbar(vwin);
}

static gint query_save_text (GtkWidget *w, GdkEvent *event, 
			     windata_t *vwin)
{
    if (CONTENT_IS_CHANGED(vwin)) {
	int resp = yes_no_dialog("gretl", 
				 _("Save changes?"), 1);

	if (resp == GRETL_CANCEL) {
	    return TRUE;
	}

	if (resp == GRETL_YES) {
	    if (vwin->role == EDIT_HEADER || vwin->role == EDIT_NOTES) {
		buf_edit_save(NULL, vwin);
	    } else if (vwin->role == EDIT_SCRIPT) {
		auto_save_script(vwin);
	    }
	}
    }

    return FALSE;
}

windata_t *edit_buffer (char **pbuf, int hsize, int vsize, 
			char *title, int role) 
{
    GtkWidget *close;
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;
#endif
    windata_t *vwin;

    vwin = common_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) {
	return NULL;
    }

#ifdef OLD_GTK
    create_text(vwin, hsize, vsize, TRUE);
#endif

    viewer_box_config(vwin); 

    /* add a menu bar */
    make_viewbar(vwin, 0);

#ifndef OLD_GTK
    create_text(vwin, &tbuf, hsize, vsize, TRUE);
#endif

    text_table_setup(vwin);
    
    /* insert the buffer text */
#ifndef OLD_GTK
    if (*pbuf) {
	gtk_text_buffer_set_text(tbuf, *pbuf, -1);
    }

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);
#else
    if (*pbuf) {
	gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
			NULL, NULL, *pbuf, strlen(*pbuf));
    } else {
	gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
			NULL, NULL, "A", 1);
	gtk_editable_delete_text(GTK_EDITABLE(vwin->w), 0, -1);
    }

    gtk_signal_connect(GTK_OBJECT(vwin->dialog), "key_press_event", 
		       GTK_SIGNAL_FUNC(catch_edit_key), vwin);	
#endif	

    /* check on delete event? */
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tbuf), "changed", 
                     G_CALLBACK(content_changed), vwin);
#else
    gtk_signal_connect(GTK_OBJECT(vwin->w), "changed", 
                       GTK_SIGNAL_FUNC(content_changed), vwin);
#endif   
    g_signal_connect(G_OBJECT(vwin->dialog), "delete_event",
		     G_CALLBACK(query_save_text), vwin);

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

#ifndef OLD_GTK
    cursor_to_top(vwin);
#endif

    return vwin;
}

static gint 
check_delete_model_window (GtkWidget *w, GdkEvent *e, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    gint ret = FALSE;

    if (window_is_busy(vwin)) {
	ret = TRUE;
    }

    return ret;
}

int view_model (PRN *prn, MODEL *pmod, int hsize, int vsize, 
		char *title) 
{
    windata_t *vwin;
    GtkWidget *close;
    const char *pbuf;
#ifndef OLD_GTK
    GtkTextBuffer *tbuf;
#endif

    vwin = common_viewer_new(VIEW_MODEL, title, pmod, 1);
    if (vwin == NULL) {
	return 1;
    }

#ifdef OLD_GTK
    create_text(vwin, hsize, vsize, FALSE);
#endif

    viewer_box_config(vwin);

    set_up_viewer_menu(vwin->dialog, vwin, model_items);
    add_vars_to_plot_menu(vwin);
    add_model_dataset_items(vwin);
    if (latex_is_ok() && !pmod->errcode) {
	add_model_tex_items(vwin);
    }

    if (pmod->ci != ARMA && pmod->ci != GARCH && pmod->ci != NLS) {
	add_dummies_to_plot_menu(vwin);
    }

    g_signal_connect(G_OBJECT(vwin->mbar), "button_press_event", 
		     G_CALLBACK(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

#ifndef OLD_GTK
    create_text(vwin, &tbuf, hsize, vsize, FALSE);
#endif

    text_table_setup(vwin);

    /* close button */
    close = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), close, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(close), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(close);

    pbuf = gretl_print_get_buffer(prn);

    /* insert and then free the model buffer */
#ifndef OLD_GTK
    gtk_text_buffer_set_text(tbuf, pbuf, strlen(pbuf));
#else
    gtk_text_insert(GTK_TEXT(vwin->w), fixed_font, 
		    NULL, NULL, pbuf, strlen(pbuf));
#endif
    gretl_print_destroy(prn);

    /* attach shortcuts */
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
#else
    gtk_signal_connect(GTK_OBJECT(vwin->dialog), "key_press_event", 
		       GTK_SIGNAL_FUNC(catch_viewer_key), 
		       vwin);
#endif

    /* don't allow deletion of model window when a model
       test dialog is active */
    g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
		     G_CALLBACK(check_delete_model_window), 
		     vwin);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(delete_unnamed_model), 
		     vwin->data);
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), 
		     vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show_all(vwin->dialog);

#ifndef OLD_GTK
    cursor_to_top(vwin);
#endif

    return 0;
}

static void auto_save_script (windata_t *vwin)
{
    FILE *fp;
    char msg[MAXLEN];
    gchar *savestuff;
    int unsaved = 0;

    if (strstr(vwin->fname, "script_tmp") || *vwin->fname == '\0') {
	file_save(vwin, SAVE_SCRIPT, NULL);
	strcpy(vwin->fname, scriptfile);
	unsaved = 1;
    }

    if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	sprintf(msg, _("Couldn't write to %s"), vwin->fname);
	errbox(msg); 
	return;
    }

#ifndef OLD_GTK
    savestuff = textview_get_text(GTK_TEXT_VIEW(vwin->w));
#else
    savestuff = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
#endif
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);

    mark_content_saved(vwin);
}

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

static void model_tex_equation_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/LaTeX/View/Equation", s);
    flip(ifac, "/LaTeX/Save/Equation", s);
    flip(ifac, "/LaTeX/Copy/Equation", s);
}

static void minimal_model_check (GtkItemFactory *ifac, const MODEL *pmod)
{
    if (pmod->ncoeff == 1) {
	flip(ifac, "/Tests/omit variables", FALSE);
	flip(ifac, "/Tests/sum of coefficients", FALSE);

	if (pmod->ifc) {
	    flip(ifac, "/Tests/heteroskedasticity", FALSE);
	    flip(ifac, "/Tests/non-linearity (squares)", FALSE);
	    flip(ifac, "/Tests/non-linearity (logs)", FALSE);
	}
    }

    if (pmod->ncoeff - pmod->ifc <= 1) {
	flip(ifac, "/Tests/collinearity", FALSE);
    }

    if (pmod->missmask != NULL) {
	flip(ifac, "/Tests/autocorrelation", FALSE);
	flip(ifac, "/Tests/CUSUM test", FALSE);
	flip(ifac, "/Tests/ARCH", FALSE);
    }
}

static void set_tests_menu_state (GtkItemFactory *ifac, const MODEL *pmod)
{
    int i, cmd_ci, ok;

    for (i=0; model_items[i].path != NULL; i++) {
	if (model_items[i].item_type == NULL &&
	    strstr(model_items[i].path, "Tests")) {
	    cmd_ci = model_items[i].callback_action;
	    if (cmd_ci == LMTEST_SQUARES || 
		cmd_ci == LMTEST_LOGS ||
		cmd_ci == LMTEST_WHITE) {
		cmd_ci = LMTEST;
	    }
	    ok = command_ok_for_model(cmd_ci, pmod->ci);
	    flip(ifac, model_items[i].path, ok);
	}
    }

    /* cross-sectional data: disallow time-series tests */
    if (!dataset_is_time_series(datainfo)) {
	flip(ifac, "/Tests/CUSUM test", FALSE);
	flip(ifac, "/Tests/ARCH", FALSE);
    }

    if (!dataset_is_time_series(datainfo) && !dataset_is_panel(datainfo)) {
	flip(ifac, "/Tests/autocorrelation", FALSE);
    }

    minimal_model_check(ifac, pmod);
}

static void arch_menu_off (GtkItemFactory *ifac)
{
    flip(ifac, "/Tests/ARCH", FALSE);
}

static void model_save_state (GtkItemFactory *ifac, gboolean s)
{
    flip(ifac, "/File/Save to session as icon", s);
    flip(ifac, "/File/Save as icon and close", s);
}

static void arma_x12_menu_mod (windata_t *vwin)
{
    flip(vwin->ifac, "/Model data/coefficient covariance matrix", FALSE);
    add_x12_output_menu_item(vwin);
}

static void adjust_model_menu_state (windata_t *vwin, const MODEL *pmod)
{
    set_tests_menu_state(vwin->ifac, pmod);

    /* disallow saving an already-saved model */
    if (pmod->name != NULL) {
	model_save_state(vwin->ifac, FALSE);
    }

    if (pmod->ci == ARMA && arma_by_x12a(pmod)) {
	arma_x12_menu_mod(vwin);
    }	

    if (dataset_is_panel(datainfo)) {
	arch_menu_off(vwin->ifac);
	panel_heteroskedasticity_menu(vwin);
    }
}

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[])
{
#ifdef OLD_GTK
    GtkAccelGroup *accel = gtk_accel_group_new();
#endif
    gint n_items = 0;

    while (items[n_items].path != NULL) n_items++;

#ifdef OLD_GTK
    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", accel);
#else
    vwin->ifac = gtk_item_factory_new(GTK_TYPE_MENU_BAR, "<main>", NULL);
#endif

# ifdef ENABLE_NLS
    gtk_item_factory_set_translate_func(vwin->ifac, menu_translate, NULL, NULL);
# endif
    gtk_item_factory_create_items(vwin->ifac, n_items, items, vwin);
    vwin->mbar = gtk_item_factory_get_widget(vwin->ifac, "<main>");
#ifdef OLD_GTK
    gtk_accel_group_attach(accel, GTK_OBJECT(window));
#endif

    if (vwin->data == NULL) return;

    if (vwin->role == VIEW_MODEL) { 
	MODEL *pmod = (MODEL *) vwin->data;

	adjust_model_menu_state(vwin, pmod);
    } else if (vwin->role == VAR) {
	GRETL_VAR *var = (GRETL_VAR *) vwin->data;
	const char *name = gretl_VAR_get_name(var);

	if (name != NULL && *name != '\0') {
	    model_save_state(vwin->ifac, FALSE);
	}	
    }
}

#ifndef OLD_GTK

static GtkItemFactoryEntry model_dataset_basic_items[] = {
    { N_("/Model data/Add to data set/fitted values"), NULL, 
      fit_resid_callback, GENR_FITTED, NULL, GNULL },
    { N_("/Model data/Add to data set/residuals"), NULL, 
      fit_resid_callback, GENR_RESID, NULL, GNULL },
    { N_("/Model data/Add to data set/squared residuals"), NULL, 
      fit_resid_callback, GENR_RESID2, NULL, GNULL }
};

static GtkItemFactoryEntry ess_items[] = {
    { N_("/Model data/Add to data set/error sum of squares"), NULL, 
      model_stat_callback, ESS, NULL, GNULL },
    { N_("/Model data/Add to data set/standard error of residuals"), NULL, 
      model_stat_callback, SIGMA, NULL, GNULL }
}; 

static GtkItemFactoryEntry r_squared_items[] = {
    { N_("/Model data/Add to data set/R-squared"), NULL, 
      model_stat_callback, R2, NULL, GNULL },
    { N_("/Model data/Add to data set/T*R-squared"), NULL, 
      model_stat_callback, TR2, NULL, GNULL }
};   

static GtkItemFactoryEntry lnl_data_item = {
    N_("/Model data/Add to data set/log likelihood"), NULL, 
    model_stat_callback, LNL, NULL, GNULL 
};

static GtkItemFactoryEntry criteria_items[] = {
    { N_("/Model data/Add to data set/Akaike Information Criterion"), NULL, 
      model_stat_callback, AIC, NULL, GNULL },
    { N_("/Model data/Add to data set/Bayesian Information Criterion"), NULL, 
      model_stat_callback, BIC, NULL, GNULL }
};

static GtkItemFactoryEntry garch_data_item = {
    N_("/Model data/Add to data set/predicted error variance"), NULL, 
    fit_resid_callback, GENR_H, NULL, GNULL 
};

static GtkItemFactoryEntry define_var_items[] = {
    { N_("/Model data/sep1"), NULL, NULL, 0, "<Separator>", GNULL },
    { N_("/Model data/Define new variable..."), NULL, model_genr_callback,
      MODEL_GENR, NULL, GNULL }
};

#else /* old GTK versions */

static GtkItemFactoryEntry model_dataset_basic_items[] = {
    { N_("/Model data/Add to data set/fitted values"), NULL, 
      fit_resid_callback, GENR_FITTED, NULL },
    { N_("/Model data/Add to data set/residuals"), NULL, 
      fit_resid_callback, GENR_RESID, NULL },
    { N_("/Model data/Add to data set/squared residuals"), NULL, 
      fit_resid_callback, GENR_RESID2, NULL },
    { N_("/Model data/Add to data set/degrees of freedom"), NULL, 
      model_stat_callback, DF, NULL }
};

static GtkItemFactoryEntry ess_items[] = {
    { N_("/Model data/Add to data set/error sum of squares"), NULL, 
      model_stat_callback, ESS, NULL },
    { N_("/Model data/Add to data set/standard error of residuals"), NULL, 
      model_stat_callback, SIGMA, NULL }
}; 

static GtkItemFactoryEntry r_squared_items[] = {
    { N_("/Model data/Add to data set/R-squared"), NULL, 
      model_stat_callback, R2, NULL },
    { N_("/Model data/Add to data set/T*R-squared"), NULL, 
      model_stat_callback, TR2, NULL }
};  

static GtkItemFactoryEntry lnl_data_item = {
    N_("/Model data/Add to data set/log likelihood"), NULL, 
    model_stat_callback, LNL, NULL
};

static GtkItemFactoryEntry criteria_items[] = {
    { N_("/Model data/Add to data set/Akaike Information Criterion"), NULL, 
      model_stat_callback, AIC, NULL },
    { N_("/Model data/Add to data set/Bayesian Information Criterion"), NULL, 
      model_stat_callback, BIC, NULL }
};

static GtkItemFactoryEntry garch_data_item = {
    N_("/Model data/Add to data set/predicted error variance"), NULL, 
    fit_resid_callback, GENR_H, NULL
};

static GtkItemFactoryEntry define_var_items[] = {
    { N_("/Model data/sep1"), NULL, NULL, 0, "<Separator>" },
    { N_("/Model data/Define new variable..."), NULL, model_genr_callback,
      MODEL_GENR, NULL }
};

#endif /* GTK versions alternates */

static void add_model_dataset_items (windata_t *vwin)
{
    int i, n;
    MODEL *pmod = vwin->data;

    n = sizeof model_dataset_basic_items / 
	sizeof model_dataset_basic_items[0];

    for (i=0; i<n; i++) {
	gtk_item_factory_create_item(vwin->ifac, &model_dataset_basic_items[i], 
				     vwin, 1);
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

    n = sizeof criteria_items / sizeof criteria_items[0];
    for (i=0; i<n; i++) {
	gtk_item_factory_create_item(vwin->ifac, &criteria_items[i], vwin, 1);
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

    for (i=0; i<n; i++) {
	gtk_item_factory_create_item(vwin->ifac, &model_tex_items[i], 
				     vwin, 1);
    }  

    model_tex_equation_state(vwin->ifac, !pmod->errcode && pmod->ci != NLS);
}

static void add_vars_to_plot_menu (windata_t *vwin)
{
    int i, j, varstart;
    GtkItemFactoryEntry varitem;
    const gchar *mpath[] = {
	N_("/Graphs/residual plot"), 
	N_("/Graphs/fitted, actual plot")
    };
    MODEL *pmod = vwin->data;
    char tmp[16];

   varitem.accelerator = NULL; 
   varitem.item_type = NULL;
   varitem.callback_action = 0; 

    for (i=0; i<2; i++) {
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

	if (pmod->ci == ARMA || pmod->ci == NLS || pmod->ci == GARCH) 
	    continue;

	varstart = (i == 0)? 1 : 2;

	/* put the indep vars on the menu list */
	for (j=varstart; j<=pmod->list[0]; j++) {
	    if (pmod->list[j] == 0) continue;
	    if (pmod->list[j] == LISTSEP) break;
	    if (!strcmp(datainfo->varname[pmod->list[j]], "time")) 
		continue;

	    varitem.callback_action = pmod->list[j]; 
	    double_underscores(tmp, datainfo->varname[pmod->list[j]]);
	    varitem.path = 
		g_strdup_printf(_("%s/against %s"), mpath[i], tmp);
	    varitem.callback = (i==0)? resid_plot : fit_actual_plot;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	}

	varitem.callback_action = 0;

	/* if the model has two independent vars, offer a 3-D fitted
	   versus actual plot */
	if (i == 1 && pmod->ifc && pmod->ncoeff == 3) {
	    char tmp2[16];

	    double_underscores(tmp, datainfo->varname[pmod->list[3]]);
	    double_underscores(tmp2, datainfo->varname[pmod->list[4]]);
	    varitem.path =
		g_strdup_printf(_("%s/against %s and %s"),
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
	N_("/Graphs/dumsep"), 
	N_("/Graphs/Separation")
    };
    gchar *radiopath = NULL;
    char tmp[16];
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
	    dumitem.path = g_strdup(_(mpath[0]));
	    gtk_item_factory_create_item(vwin->ifac, &dumitem, vwin, 1);
	    g_free(dumitem.path);

	    /* add menu branch */
	    dumitem.item_type = "<Branch>";
	    dumitem.path = g_strdup(_(mpath[1]));
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
	dumitem.path = g_strdup_printf(_("%s/by %s"), mpath[1], tmp);
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

static void VAR_model_data_callback (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = vwin->data;
    int h, neqns;
    gchar *title;
    PRN *prn;
    int i, err;

    if (var == NULL) return;
    if (bufopen(&prn)) return;

    if (code != VAR_VCV) {
	neqns = gretl_VAR_get_n_equations(var);
	h = default_VAR_horizon(datainfo);
	title = g_strdup_printf("gretl: %s", 
				(code == VAR_IRF)? _("impulse responses") :
				_("variance decompositions"));
	err = checks_dialog(title, NULL, 0, NULL,
			    &h, _("forecast horizon (periods):"),
			    2, datainfo->n / 2, 0);
	g_free(title);
	if (err < 0) {
	    gretl_print_destroy(prn);
	    return;
	} 
    }

    if (code == VAR_VCV) {
	title = g_strdup(_("gretl: VAR covariance matrix"));
	err = gretl_VAR_print_VCV(var, prn);
    } else if (code == VAR_IRF) {
	title = g_strdup(_("gretl: VAR impulse responses"));
	for (i=0; i<neqns && !err; i++) {
	    err = gretl_VAR_print_impulse_response(var, i, h, datainfo, 0, prn);
	}
    } else if (code == VAR_DECOMP) {
	title = g_strdup(_("gretl: VAR variance decompositions"));
	for (i=0; i<neqns && !err; i++) {
	    err = gretl_VAR_print_fcast_decomp(var, i, h, datainfo, 0, prn);
	}
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
    }

    g_free(title);
}

static void panel_heteroskedasticity_menu (windata_t *vwin)
{
    GtkItemFactoryEntry hitem;

    gtk_item_factory_delete_item(vwin->ifac, "/Tests/heteroskedasticity");

    hitem.accelerator = NULL;
    hitem.item_type = NULL;

    hitem.callback = do_lmtest;
    hitem.callback_action = LMTEST_WHITE;
    hitem.path = g_strdup(_("/Tests/heteroskedasticity (White's test)"));
    gtk_item_factory_create_item(vwin->ifac, &hitem, vwin, 1);
    g_free(hitem.path);

    hitem.callback = do_lmtest;
    hitem.callback_action =  LMTEST_GROUPWISE;
    hitem.path = g_strdup(_("/Tests/heteroskedasticity (groupwise)"));
    gtk_item_factory_create_item(vwin->ifac, &hitem, vwin, 1);
    g_free(hitem.path);
}

static void add_x12_output_menu_item (windata_t *vwin)
{
    GtkItemFactoryEntry item;
    const gchar *mpath = "/Model data";

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
    item.path = g_strdup_printf("%s/%s", mpath, _("view X-12-ARIMA output"));
    gtk_item_factory_create_item(vwin->ifac, &item, vwin, 1);
    g_free(item.path);
}

static void impulse_plot_call (gpointer p, guint shock, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    gchar *title;
    int h = default_VAR_horizon(datainfo);
    gint targ;
    int err;
    const double **vZ = NULL;
    const char *impulse_opts[] = {
	N_("include bootstrap confidence interval")
    };
    static int active[] = { 0 };

    targ = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "targ"));

    title = g_strdup_printf("gretl: %s", _("impulse responses"));

    err = checks_dialog(title, impulse_opts, 1, active,
			&h, _("forecast horizon (periods):"),
			2, datainfo->n / 2, IRF_BOOT);
    g_free(title);

    if (err < 0) {
	return;
    } 

    if (active[0]) {
	vZ = (const double **) Z;
    }

    err = gretl_VAR_plot_impulse_response(var, targ, shock, h, vZ,
					  datainfo);

    if (!err) {
	register_graph();
    }
}

static void VAR_forecast_callback (gpointer p, guint i, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    const MODEL *pmod;
    FITRESID *fr;
    int t1, t2, resp;
    int premax, pre_n, dyn_ok;
    gretlopt opt = OPT_NONE;
    int err = 0;

    pmod = gretl_VAR_get_model(var, i);

    t2 = datainfo->n - 1;

    /* if no out-of-sample obs are available, alert the user */
    if (t2 == pmod->t2) {
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
    if (t2 > pmod->t2) {
	t1 = pmod->t2 + 1;
	pre_n = pmod->t2 / 2;
	if (pre_n > 100) {
	    pre_n = 100;
	}
	dyn_ok = 1;
    } else {
	t1 = 0;
	pre_n = 0;
	dyn_ok = 0;
    }

    resp = forecast_dialog(t1, t1, &t1,
			   t1, t2, &t2,
			   0, premax, &pre_n,
			   dyn_ok);
    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    }

    fr = get_VAR_forecast(var, i, t1, t2, pre_n, (const double **) Z, 
			  datainfo, opt);

    if (fr == NULL) {
	errbox("Forecast failed");
    } else {
	int width = 78;
	PRN *prn;

	if (bufopen(&prn)) {
	    return;
	}

	err = text_print_forecast(fr, &Z, datainfo, OPT_P, prn);
	if (!err) {
	    register_graph();
	}
	if (fr->sderr == NULL) {
	    width = 50;
	}
	view_buffer(prn, width, 400, _("gretl: forecasts"), FCASTERR, fr);
    }
}

enum {
    VAR_AUTOCORR_TEST,
    VAR_NORMALITY_TEST
};

static void VAR_test_call (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    char title[72];
    PRN *prn;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    if (code == VAR_AUTOCORR_TEST) {
	int order = default_lag_order(datainfo);

	set_window_busy(vwin);
	err = spin_dialog(_("gretl: autocorrelation"), 
			  &order, _("Lag order for test:"),
			  1, datainfo->n / 2, LMTEST);
	unset_window_busy(vwin);

	if (err < 0) {
	    gretl_print_destroy(prn);
	    return;
	}

	strcpy(title, _("gretl: LM test (autocorrelation)"));
	err = gretl_VAR_autocorrelation_test(var, order, 
					     &Z, datainfo, 
					     prn);
    } else if (code == VAR_NORMALITY_TEST) {
	sprintf(title, "gretl: %s", _("Test for normality of residual"));
	err = gretl_VAR_normality_test(var, prn);
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

static void VAR_resid_plot_call (gpointer p, guint i, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_residual_plot(var, &Z, datainfo);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void add_VAR_menu_items (windata_t *vwin)
{
    GtkItemFactoryEntry varitem;
    const gchar *tpath = N_("/Tests");
    const gchar *gpath = N_("/Graphs");
    const gchar *mpath = N_("/Model data");
    const gchar *fpath = N_("/Model data/Forecasts");
    const gchar *dpath = N_("/Model data/Add to data set");
    GRETL_VAR *var = vwin->data;
    int neqns = gretl_VAR_get_n_equations(var);
    int vtarg, vshock;
    char tmp[16];
    int i, j;

    varitem.accelerator = NULL;
    varitem.callback = NULL;
    varitem.callback_action = 0;
    varitem.item_type = "<Branch>";

    varitem.path = g_strdup(_("/_Tests"));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    varitem.path = g_strdup(_("/_Model data"));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    varitem.path = g_strdup(_("/_Graphs"));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    varitem.path = g_strdup(_(fpath));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);
    
    varitem.path = g_strdup(_(dpath));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* autocorrelation tests */
    varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				   _("autocorrelation"));
    varitem.callback = VAR_test_call;
    varitem.callback_action = VAR_AUTOCORR_TEST;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* multivariate normality test */
    varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				   _("normality of residuals"));
    varitem.callback = VAR_test_call;
    varitem.callback_action = VAR_NORMALITY_TEST;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* cross-equation VCV */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), 
				   _("Cross-equation covariance matrix"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_VCV;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* impulse response printout */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), _("impulse responses"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_IRF;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);    

    /* variance decomp printout */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), 
				   _("forecast variance decomposition"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_DECOMP;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path); 

    /* combined residual plot */
    varitem.path = g_strdup_printf("%s/%s", _(gpath), _("residual plot"));
    varitem.callback = VAR_resid_plot_call;
    varitem.callback_action = 0;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    for (i=0; i<neqns; i++) {
	char maj[32], min[16];
	int dv;

	/* forecast items */
	dv = gretl_VAR_get_variable_number(var, i);
	varitem.path = g_strdup_printf("%s/%s", _(fpath), 
				       datainfo->varname[dv]);
	varitem.callback = VAR_forecast_callback;
	varitem.callback_action = i;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	/* save resids items */
	varitem.path = g_strdup_printf("%s/%s %d", _(dpath), 
				       _("residuals from equation"), i + 1);
	varitem.callback = VAR_resid_callback;
	varitem.callback_action = i;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);	

	/* impulse response plots: make branch for target */
	vtarg = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, datainfo->varname[vtarg]);
	sprintf(maj, _("response of %s"), tmp);

	varitem.path = g_strdup_printf("%s/%s", _(gpath), maj);
	varitem.callback = NULL;
	varitem.callback_action = 0;
	varitem.item_type = "<Branch>";
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	varitem.item_type = NULL;
	
	for (j=0; j<neqns; j++) {
	    GtkWidget *w;

	    /* impulse responses: subitems for shocks */
	    vshock = gretl_VAR_get_variable_number(var, j);
	    varitem.callback_action = j;
	    double_underscores(tmp, datainfo->varname[vshock]);
	    sprintf(min, _("to %s"), tmp);

	    varitem.path = g_strdup_printf("%s/%s/%s", _(gpath), maj, min);
	    varitem.callback = impulse_plot_call;
	    varitem.callback_action = j;
	    varitem.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	    w = gtk_item_factory_get_widget_by_action(vwin->ifac, j);
#ifndef OLD_GTK
	    g_object_set_data(G_OBJECT(w), "targ", GINT_TO_POINTER(i));
#else
	    gtk_object_set_data(GTK_OBJECT(w), "targ", GINT_TO_POINTER(i));
#endif
	}
    }

    if (latex_is_ok()) {
	int n = sizeof VAR_tex_items / sizeof VAR_tex_items[0];

	for (i=0; i<n; i++) {
	    gtk_item_factory_create_item(vwin->ifac, &VAR_tex_items[i], 
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
	flip(mwin->ifac, "/Model data", FALSE);
	flip(mwin->ifac, "/Tests", FALSE);
	flip(mwin->ifac, "/Graphs", FALSE);
	flip(mwin->ifac, "/Model data", FALSE);
	flip(mwin->ifac, "/LaTeX", FALSE);

	return FALSE;
    }

    if (model_sample_issue(pmod, NULL, 0, datainfo)) {
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
    flip(mwin->ifac, "/Model data/Display actual, fitted, residual", ok);
    flip(mwin->ifac, "/Model data/Forecasts...", ok);
    flip(mwin->ifac, "/Model data/Confidence intervals for coefficients", ok);
    flip(mwin->ifac, "/Model data/Add to data set/fitted values", ok);
    flip(mwin->ifac, "/Model data/Add to data set/residuals", ok);
    flip(mwin->ifac, "/Model data/Add to data set/squared residuals", ok);
    flip(mwin->ifac, "/Model data/Define new variable...", ok);

    if (!ok) {
	const char *msg = get_gretl_errmsg();

	if (msg != NULL && *msg != 0) {
	    infobox(msg);
	}
    } 

    return FALSE;
}

int validate_varname (const char *varname)
{
    int i, n = strlen(varname);
    char namebit[USER_VLEN];
    unsigned char c;

    *namebit = 0;
    
    if (n > USER_VLEN - 1) {
	strncat(namebit, varname, USER_VLEN - 1);
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

#ifndef OLD_GTK
int my_utf_validate (char *s)
{
    if (!g_utf8_validate(s, -1, NULL)) {
	gchar *new = my_locale_to_utf8(s);

	if (new != NULL) {
	    strcpy(s, new);
	    g_free(new);
	} else {
	    *s = '\0';
	}
	return 1;
    }
    return 0;
}
#endif	

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
#ifndef OLD_GTK
		     GCallback callback, 
#else
		     GtkSignalFunc callback, 
#endif
		     gpointer data)
{
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(label);
#ifndef OLD_GTK
    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
#else
    gtk_menu_append(GTK_MENU(menu), item);
#endif
    g_signal_connect (G_OBJECT(item), "activate",
		      G_CALLBACK(callback), data);
    gtk_widget_show(item);
}

/* .................................................................. */

void *gui_get_plugin_function (const char *funcname, 
			       void **phandle)
{
    void *func;

    func = get_plugin_function(funcname, phandle);
    if (func == NULL) {
	errbox(get_gretl_errmsg());
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

#ifndef OLD_GTK

static int seven_bit_string (const unsigned char *s)
{
    while (*s) {
	if (*s > 127) return 0;
	s++;
    }
    return 1;
}

static int seven_bit_file (const char *fname)
{
    FILE *fp;
    char line[256];
    int ascii = 1;
    
    fp = gretl_fopen(fname, "r");
    if (fp == NULL) {
	return 1;
    }

    while (fgets(line, sizeof line, fp)) {
	if (!seven_bit_string(line)) {
	    ascii = 0;
	    break;
	}
    }

    fclose(fp);

    return ascii;
}

static int maybe_recode_file (const char *fname)
{
    const gchar *charset;

    if (g_get_charset(&charset)) {
	/* locale uses UTF-8 */
	return 0;
    }

    if (seven_bit_file(fname)) {
	return 0;
    } else {
	FILE *fin, *fout;
	char trname[MAXLEN];
	char line[128];
	gchar *trbuf;
	int err = 0;

	fin = gretl_fopen(fname, "r");
	if (fin == NULL) {
	    return 1;
	}

	sprintf(trname, "%s.tr", fname);

	fout = gretl_fopen(trname, "w");
	if (fout == NULL) {
	    fclose(fin);
	    return 1;
	}	

	while (fgets(line, sizeof line, fin) && !err) {
	    trbuf = my_locale_from_utf8(line);
	    if (trbuf != NULL) {
		fputs(trbuf, fout);
		g_free(trbuf);
	    } else {
		err = 1;
	    }
	}

	fclose(fin);
	fclose(fout);

	if (!err) {
	    err = copyfile(trname, fname);
	    remove(trname);
	}

	return err;
    }

    return 0;
}

gchar *my_filename_from_utf8 (char *fname)
{
    gchar *trfname;
    gsize bytes;
    GError *err = NULL;

    if (seven_bit_string(fname)) {
	return fname;
    }

    trfname = g_filename_from_utf8(fname, -1, NULL, &bytes, &err);

    if (err != NULL) {
	errbox("g_filename_from_utf8 failed");
	g_error_free(err);
    } else {
	strcpy(fname, trfname);
    }

    g_free(trfname);

    return fname;
}

#define TR_DEBUG 0

static gchar *
real_locale_from_utf8 (const gchar *src, int force)
{
    gchar *trstr;
    gsize bytes;
    GError *err = NULL;
    const gchar *cset = NULL;

#if TR_DEBUG
    gchar *msg;
    int u = g_get_charset(&cset);
    msg = g_strdup_printf("real_locale_from_utf8: force=%d, "
			  "g_get_charset returned %d and gave "
			  "'%s'\n", force, u, 
			  (cset != NULL)? cset : "NULL");
    infobox(msg);
    g_free(msg);
    if (!force && u) {
	return g_strdup(src);
    }
#else
    if (!force && g_get_charset(&cset)) {
	/* According to the glib manual, g_get_charset returns TRUE if 
	   the returned charset is UTF-8 */ 
	return g_strdup(src);
    }
#endif

    trstr = g_locale_from_utf8(src, -1, NULL, &bytes, &err);

    if (err != NULL) {
	if (cset != NULL) {
	    sprintf(errtext, "g_locale_from_utf8 failed for charset '%s'",
		    cset);
	} else {
	    strcpy(errtext, "g_locale_from_utf8 failed; "
		   "so did g_get_charset");
	}
	errbox(errtext);
	g_error_free(err);
    }

    return trstr;
}

gchar *my_locale_from_utf8 (const gchar *src)
{
    return real_locale_from_utf8(src, 0);
}

gchar *force_locale_from_utf8 (const gchar *src)
{
    return real_locale_from_utf8(src, 1);
}

gchar *my_filename_to_utf8 (char *fname)
{
    gchar *trfname;
    gsize bytes;
    GError *err = NULL;

    if (g_utf8_validate(fname, -1, NULL)) {
	return fname;
    }

    trfname = g_filename_to_utf8(fname, -1, NULL, &bytes, &err);

    if (err != NULL) {
	errbox("g_filename_to_utf8 failed");
	g_error_free(err);
    } else {
	strcpy(fname, trfname);
    }

    g_free(trfname);

    return fname;
}

gchar *my_locale_to_utf8 (const gchar *src)
{
    gchar *trstr;
    gsize bytes;
    GError *err = NULL;

    trstr = g_locale_to_utf8(src, -1, NULL, &bytes, &err);

    if (err != NULL) {
	const gchar *cset = NULL;

	g_get_charset(&cset);
	if (cset != NULL) {
	    sprintf(errtext, "g_locale_to_utf8 failed for charset '%s'",
		    cset);
	} else {
	    strcpy(errtext, "g_locale_to_utf8 failed; "
		   "so did g_get_charset");
	}
	errbox(errtext);
	g_error_free(err);
    }

    return trstr;
}

#endif /* ! OLD_GTK */

#ifndef G_OS_WIN32

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
	errbox(_("Please open a data file first"));
	return;
    }

    build_path(paths.userdir, "gretl.Rprofile", Rprofile, NULL);
    fp = fopen(Rprofile, "w");
    if (fp == NULL) {
	errbox(_("Couldn't write R startup file"));
	return;
    }

    enverr = setenv("R_PROFILE", Rprofile, 1);
    if (enverr) {
	errbox(_("Couldn't set R_PROFILE environment variable"));
	fclose(fp);
	return;
    } 	

    build_path(paths.userdir, "Rdata.tmp", Rdata, NULL);

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
	fputs("vnum <- as.double(R.version$major) + (as.double(R.version$minor) / 10.0)\n", fp);
	fputs("if (vnum > 1.89) library(stats) else library(ts)\n", fp);
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

#ifdef OLD_GTK /* for forwards compatibility */

static gint entry_activate (GtkWidget *w, GdkEventKey *key, gpointer p)
{
    GtkWidget *top = gtk_widget_get_toplevel(w);

    gtk_window_activate_default(GTK_WINDOW(top));

    return FALSE;
}

void gtk_entry_set_activates_default (GtkEntry *entry, gboolean setting)
{
    gtk_signal_connect(GTK_OBJECT(entry), "activate", 
		       GTK_SIGNAL_FUNC(entry_activate), NULL);
}

#endif /* old GTK */
