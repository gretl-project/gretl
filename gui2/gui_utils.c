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
#include "varprint.h"
#include "modelspec.h"
#include "forecast.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "system.h"

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

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

#ifdef USE_GTKSOURCEVIEW
# include <gtksourceview/gtksourceview.h>
#endif

char *storelist = NULL;

#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mini.tsplot.xpm"
#include "../pixmaps/mini.boxplot.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.manual.xpm"

#define CONTENT_IS_CHANGED(w) (w->active_var == 1)

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkItemFactoryEntry items[]);
static void view_window_save (GtkWidget *widget, windata_t *vwin);
static gint query_save_text (GtkWidget *w, GdkEvent *event, windata_t *vwin);
static void auto_save_script (windata_t *vwin);
static void add_model_dataset_items (windata_t *vwin);
static void add_model_tex_items (windata_t *vwin);
static void add_vars_to_plot_menu (windata_t *vwin);
static void add_dummies_to_plot_menu (windata_t *vwin);
static void add_VAR_menu_items (windata_t *vwin, int vecm);
static void add_SYS_menu_items (windata_t *vwin);
static void add_x12_output_menu_item (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);
static void buf_edit_save (GtkWidget *widget, gpointer data);
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
    ADD_ITEM,
    MAIL_ITEM,
    HELP_ITEM,
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    CODE_ITEM,
    INDEX_ITEM
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

    if (!window_is_busy(vwin)) {
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
#ifdef G_OS_WIN32
	latex_ok = check_for_prog("latex.exe");
#else
	latex_ok = check_for_prog("latex");
#endif
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
      GRETL_OBJ_EQN, NULL, GNULL },
    { N_("/File/Save as icon and cl_ose"), NULL, model_add_as_icon_and_close, 
      GRETL_OBJ_EQN, NULL, GNULL },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
# endif
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
    { N_("/Tests/_Heteroskedasticity"), NULL, do_lmtest, LMTEST_WHITE, NULL, GNULL },
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
      do_forecast, FCASTERR, NULL, GNULL },
    { N_("/Analysis/_Confidence intervals for coefficients"), NULL, 
      do_coeff_intervals, 0, NULL, GNULL },
    { N_("/Analysis/Confidence _ellipse..."), NULL, 
      selector_callback, ELLIPSE, NULL, GNULL },
    { N_("/Analysis/Coefficient covariance _matrix"), NULL, 
      do_outcovmx, 0, NULL, GNULL },
    { NULL, NULL, NULL, 0, NULL, GNULL }
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
      GNULL }
};

static GtkItemFactoryEntry VAR_tex_items[] = {
    { N_("/_LaTeX"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/LaTeX/_View"), NULL, var_tex_callback, 0, NULL, GNULL },
    { N_("/LaTeX/_Copy"), NULL, var_tex_callback, 1, NULL, GNULL },
    { N_("/LaTeX/_Save"), NULL, var_tex_callback, 2, NULL, GNULL }
};

static GtkItemFactoryEntry VAR_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/_Save as..."), NULL, model_output_save_callback, 0, "<StockItem>", 
      GTK_STOCK_SAVE_AS },
    { N_("/File/Save to session as _icon"), NULL, model_add_as_icon, 
      GRETL_OBJ_VAR, NULL, GNULL },
    { N_("/File/Save as icon and cl_ose"), NULL, model_add_as_icon_and_close, 
      GRETL_OBJ_VAR, NULL, GNULL },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
# endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 1, "<StockItem>", GTK_STOCK_COPY },
    { NULL, NULL, NULL, 0, NULL, GNULL }
};

static GtkItemFactoryEntry SYS_items[] = {
    { N_("/_File"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/File/Save to session as _icon"), NULL, model_add_as_icon, 
      GRETL_OBJ_SYS, NULL, GNULL },
    { N_("/File/Save as icon and cl_ose"), NULL, model_add_as_icon_and_close, 
      GRETL_OBJ_SYS, NULL, GNULL },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("/File/_Print..."), NULL, window_print, 0, "<StockItem>", GTK_STOCK_PRINT },
# endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>", GNULL },
    { N_("/Edit/_Copy"), "", model_copy_callback, 0, "<StockItem>", GTK_STOCK_COPY },
    { N_("/_Tests"), NULL, NULL, 0, "<Branch>", GNULL },    
    { N_("/Tests/_Linear restrictions"), NULL, gretl_callback, RESTRICT, NULL, GNULL },
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
	errbox(_("Couldn't open %s"), src);
	return 1; 
    }

    if ((destfd = gretl_fopen(dest, "wb")) == NULL) {
	errbox(_("Couldn't write to %s"), dest);
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

#ifndef G_OS_WIN32

int gretl_mkdir (const char *path)
{
    int err = 0;
    extern int errno;

    errno = 0;

    if (mkdir(path, 0755)) {
	if (errno != EEXIST) { 
	    fprintf(stderr, "%s: %s\n", path, strerror(errno));
	    err = 1;
	}
    }

    return err;
}

#endif

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
	    errbox(_("Couldn't open %s"), fname);
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
	errbox(_("Couldn't open %s"), fname);
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

static int winstack (int code, GtkWidget *w, gpointer ptest, GtkWidget **pw)
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
		    if (pw != NULL) {
			*pw = wstack[i];
		    }
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
    winstack(STACK_INIT, NULL, NULL, NULL);
}
    
void winstack_destroy (void)
{
    winstack(STACK_DESTROY, NULL, NULL, NULL);
}

int winstack_match_data (gpointer p)
{
    return winstack(STACK_QUERY, NULL, p, NULL);
}

GtkWidget *match_window_by_data (gpointer p)
{
    GtkWidget *w = NULL;

    winstack(STACK_QUERY, NULL, p, &w);
    return w;
}

int highest_numbered_variable_in_winstack (void)
{
    return winstack(STACK_MAXVAR, NULL, NULL, NULL);
}

static void winstack_add (GtkWidget *w)
{
    winstack(STACK_ADD, w, NULL, NULL);
}

static void winstack_remove (GtkWidget *w)
{
    winstack(STACK_REMOVE, w, NULL, NULL);
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

void delete_widget (GtkWidget *widget, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
}

static gint catch_button_3 (GtkWidget *w, GdkEventButton *event)
{
    GdkModifierType mods;

    gdk_window_get_pointer(w->window, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON3_MASK) {
	return TRUE;
    }

    return FALSE;
}

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
		view_window_save(NULL, vwin);
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

static gint catch_viewer_key (GtkWidget *w, GdkEventKey *key, windata_t *vwin)
{
    if (gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->w))) {
	return catch_edit_key(w, key, vwin);
    }

    if (key->keyval == GDK_q) { 
        gtk_widget_destroy(w);
    } else if (key->keyval == GDK_s && Z != NULL && vwin->role == VIEW_MODEL) {
	model_add_as_icon_and_close(vwin, GRETL_OBJ_EQN, NULL);
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
	stop_talking();
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
                           action == APPEND_ASCII || \
                           action == APPEND_WF1 || \
                           action == APPEND_DTA || \
                           action == APPEND_JMULTI)

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
	errbox(_("Couldn't open %s"), fname);
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
	    errbox(_("Failed to import data"));
	}
	gretl_print_destroy(errprn);
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

static void copy_utf8_filename (char *targ, const char *src)
{
    strcpy(targ, src);
#ifdef ENABLE_NLS
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
    } else if (code == OPEN_WF1 || code == APPEND_WF1) {
	datatype = GRETL_WF1;
    } else if (code == OPEN_DTA || code == APPEND_DTA) {
	datatype = GRETL_DTA;
    } else if (code == OPEN_JMULTI || code == APPEND_JMULTI) {
	datatype = GRETL_JMULTI;
    } else if (code == OPEN_BOX) {
	datatype = GRETL_BOX_DATA;
    } else {
	/* no filetype specified: have to guess */
	PRN *prn;	

	if (bufopen(&prn)) return;
	datatype = detect_filetype(tryfile, &paths, prn);
	gretl_print_destroy(prn);
    }

    /* destroy the current data set, etc., unless we're explicitly appending */
    if (!append) {
	close_session(NULL, &Z, &datainfo);
    }

    if (datatype == GRETL_GNUMERIC || datatype == GRETL_EXCEL ||
	datatype == GRETL_WF1 || datatype == GRETL_DTA ||
	datatype == GRETL_JMULTI) {
	get_worksheet_data(tryfile, datatype, append, NULL);
	return;
    } else if (datatype == GRETL_CSV_DATA) {
	do_open_csv_box(tryfile, OPEN_CSV, append);
	return;
    } else if (datatype == GRETL_OCTAVE) {
	do_open_csv_box(tryfile, OPEN_OCTAVE, append);
	return;
    } else if (datatype == GRETL_BOX_DATA) {
	do_open_csv_box(tryfile, OPEN_BOX, 0);
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
	    err = gretl_read_gdt(&Z, &datainfo, tryfile, &paths, 
				 clear_code, errprn, 1);
	} else {
	    err = gretl_get_data(&Z, &datainfo, tryfile, &paths, 
				 clear_code, errprn);
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
	copy_utf8_filename(paths.datfile, tryfile);
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

void verify_open_session (void)
{
    if (!gretl_is_pkzip_file(tryfile)) {
	/* not a new-style zipped session file */
	do_open_script();
	return;
    }

    if (data_status && !expert) {
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

static void buf_edit_save (GtkWidget *widget, gpointer data)
{
    windata_t *vwin = (windata_t *) data;
    gchar *text;
    char **pbuf = (char **) vwin->data;

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
	session_changed(1);
	set_replay_off();
    }
}

static void update_func_code (windata_t *vwin)
{
    int iface, err = 0;

    /* callback used when editing a function in the context of
       the "function package editor" */
	
    iface = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->w), "iface"));
    err = update_function_from_script(vwin->fname, iface);
    if (err) {
	gui_errmsg(err);
    }
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

static windata_t *vwin_first_child (windata_t *vwin)
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
	} else if (vwin->role == FCASTERR || vwin->role == FCAST) {
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

	if (vwin->dialog) {
	    winstack_remove(vwin->dialog);
	}

	free(vwin);
    }
}

void gretl_stock_icons_init (void)
{
    char **xpms[] = {
	mini_tex_xpm,
	mail_16_xpm,
	mini_tsplot_xpm,
	mini_boxplot_xpm,
	mini_pdf_xpm,
	mini_manual_xpm
    };
    const char *stocks[] = {
	GRETL_STOCK_TEX,
	GRETL_STOCK_MAIL,
	GRETL_STOCK_TS,
	GRETL_STOCK_BOX,
	GRETL_STOCK_PDF,
	GRETL_STOCK_BOOK
    };
    int n = sizeof stocks / sizeof stocks[0];

    static GtkIconFactory *ifac;

    if (ifac == NULL) {
	GtkIconSource *source;
	GtkIconSet *set;
	GdkPixbuf *pbuf;
	int i;

	ifac = gtk_icon_factory_new();

	for (i=0; i<n; i++) {
	    set = gtk_icon_set_new();
	    source = gtk_icon_source_new();
	    gtk_icon_source_set_size(source, GTK_ICON_SIZE_MENU);
	    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) xpms[i]);
	    gtk_icon_source_set_pixbuf(source, pbuf);
	    g_object_unref(pbuf);
	    gtk_icon_set_add_source(set, source);
	    gtk_icon_source_free(source);
	    gtk_icon_factory_add(ifac, stocks[i], set);
	    gtk_icon_set_unref(set);
	}

	gtk_icon_factory_add_default(ifac);
    }
}

#if defined(G_OS_WIN32) || defined(USE_GNOME) 
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(vwin, 0, w);
}
#endif

/* when copying from windows where a choice of copy format
   is appropriate */

static void choose_copy_format_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_SCALAR) {
	scalar_to_clipboard(vwin);
    } else {
	copy_format_dialog(vwin, W_COPY);
    }
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

/* copying when only plain text is appropriate */

static void text_copy_callback (GtkWidget *w, gpointer p)
{
    int fmt = GRETL_FORMAT_TXT;

    if (vwin_selection_present(p)) {
	fmt = GRETL_FORMAT_SELECTION;
    }

    window_copy(p, fmt, w);
}

static void text_paste_callback (GtkWidget *w, gpointer p)
{
    text_paste(p, 0, w);
}

static void text_replace_callback (GtkWidget *w, gpointer p)
{
    text_replace(p, 0, w);
}

static void text_undo_callback (GtkWidget *w, gpointer p)
{
    text_undo(p, 0, w);
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

static void view_code_callback (GtkWidget *w, windata_t *vwin)
{
    windata_t *child = vwin_first_child(vwin);

    if (child != NULL) {
	gdk_window_raise(child->dialog->window);
    } else {
	child = gui_show_function_info(vwin->fname, VIEW_FUNC_CODE);
	if (child != NULL) {
	    vwin_add_child(vwin, child);
	}
    }
}

struct viewbar_item {
    const char *str;
    const gchar *icon;
    void (*toolfunc)();
    int flag;
};

static struct viewbar_item viewbar_items[] = {
    { N_("Save"), GTK_STOCK_SAVE, view_window_save, SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, file_save_callback, SAVE_AS_ITEM },
    { N_("Send to gnuplot"), GTK_STOCK_EXECUTE, gp_send_callback, GP_ITEM },
# if defined(G_OS_WIN32) || defined(USE_GNOME)
    { N_("Print..."), GTK_STOCK_PRINT, window_print_callback, 0 },
# endif
    { N_("Run"), GTK_STOCK_EXECUTE, do_run_script, RUN_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, text_copy_callback, COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, text_paste_callback, EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, text_find_callback, 0 },
    { N_("View code"), GTK_STOCK_PROPERTIES, view_code_callback, CODE_ITEM },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, text_replace_callback, EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, text_undo_callback, EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, series_view_sort, SORT_ITEM },    
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, series_view_sort_by, SORT_BY_ITEM },    
    { N_("Send To..."), GRETL_STOCK_MAIL, mail_script_callback, MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, script_index, INDEX_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, activate_script_help, RUN_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, window_tex_callback, TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, series_view_graph, PLOT_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, series_view_format_dialog, FORMAT_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, add_data_callback, ADD_ITEM },
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

#define editor_role(r) (r == EDIT_SCRIPT || \
                        r == EDIT_HEADER || \
                        r == EDIT_NOTES || \
                        r == EDIT_FUNC_CODE || \
                        r == GR_PLOT)

static void make_viewbar (windata_t *vwin, int text_out)
{
    GtkWidget *hbox, *button;
    void (*toolfunc)() = NULL;
    int i;

    int run_ok = (vwin->role == EDIT_SCRIPT ||
		  vwin->role == VIEW_SCRIPT ||
		  vwin->role == VIEW_LOG);
    int edit_ok = (vwin->role == EDIT_SCRIPT ||
		   vwin->role == EDIT_HEADER ||
		   vwin->role == EDIT_NOTES ||
		   vwin->role == EDIT_FUNC_CODE ||
		   vwin->role == GR_PLOT || 
		   vwin->role == GR_BOX ||
		   vwin->role == SCRIPT_OUT);
    int save_as_ok = (vwin->role != EDIT_HEADER && 
		      vwin->role != EDIT_NOTES &&
		      vwin->role != EDIT_FUNC_CODE &&
                      vwin->role != VIEW_SCALAR);
    int help_ok = (vwin->role == LEVERAGE || 
		   vwin->role == COINT2 ||
		   vwin->role == HURST ||
		   vwin->role == RMPLOT ||
		   vwin->role == MAHAL);
    int sort_ok    = (vwin->role == VIEW_SERIES);
    int sort_by_ok = has_sortable_data(vwin);
    int format_ok  = (vwin->role == VIEW_SERIES || vwin->role == VIEW_SCALAR);
    int plot_ok    = (vwin->role == VIEW_SERIES);
    int latex_ok   = latex_is_ok();

    if (MULTI_FORMAT_ENABLED(vwin->role) && latex_ok) {
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
	GtkWidget *w;

	toolfunc = viewbar_items[i].toolfunc;

	if (!edit_ok && viewbar_items[i].flag == EDIT_ITEM) {
	    continue;
	}

	if (!run_ok && viewbar_items[i].flag == RUN_ITEM) {
	    continue;
	}

	if (vwin->role != EDIT_SCRIPT && viewbar_items[i].flag == MAIL_ITEM) {
	    continue;
	}	

	if (!help_ok && viewbar_items[i].flag == HELP_ITEM) {
	    continue;
	}

	if (vwin->role != GR_PLOT && viewbar_items[i].flag == GP_ITEM) {
	    continue;
	}

	if (vwin->role == VIEW_SCALAR && viewbar_items[i].flag == 0) {
	    continue;
	}

	if (vwin->role != VIEW_FUNC_INFO && viewbar_items[i].flag == CODE_ITEM) {
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

	if (!sort_ok && viewbar_items[i].flag == SORT_ITEM) {
	    continue;
	}

	if (!sort_by_ok && viewbar_items[i].flag == SORT_BY_ITEM) {
	    continue;
	}

	if (!plot_ok && viewbar_items[i].flag == PLOT_ITEM) {
	    continue;
	}

	if (!format_ok && viewbar_items[i].flag == FORMAT_ITEM) {
	    continue;
	}

	if (vwin->role != VIEW_SCRIPT && 
	    viewbar_items[i].flag == INDEX_ITEM) {
	    continue;
	}

	if (viewbar_items[i].flag == COPY_ITEM && 
	    !editor_role(vwin->role)) {
	    toolfunc = choose_copy_format_callback;
	}

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

	if (viewbar_items[i].flag == SAVE_AS_ITEM) {
	    if (!save_as_ok) {
		continue;
	    } else if (MULTI_FORMAT_ENABLED(vwin->role) ||
		       (vwin->role == PRINT && vwin->data != NULL)) {
		toolfunc = multi_save_as_callback;
	    }
	}

	if (viewbar_items[i].flag == PLOT_ITEM) {
	    set_plot_icon(&viewbar_items[i]);
	}

	button = gtk_image_new();
	gtk_image_set_from_stock(GTK_IMAGE(button), viewbar_items[i].icon, 
				 GTK_ICON_SIZE_MENU);
        w = gtk_toolbar_append_item(GTK_TOOLBAR(vwin->mbar),
				    NULL, _(viewbar_items[i].str), NULL,
				    button, toolfunc, vwin);

	g_object_set_data(G_OBJECT(w), "flag", 
			  GINT_TO_POINTER(viewbar_items[i].flag));

	if (viewbar_items[i].flag == SAVE_ITEM) { 
	    /* nothing to save just yet */
	    gtk_widget_set_sensitive(w, FALSE);
	} 
	if (viewbar_items[i].flag == SAVE_AS_ITEM &&
	    strstr(vwin->fname, "script_tmp")) {
	    gtk_widget_set_sensitive(w, FALSE);
	}
    }

    gtk_widget_show(vwin->mbar);
    gtk_widget_show(hbox);
}

static void add_edit_items_to_viewbar (windata_t *vwin)
{
    GtkWidget *button;
    int i, pos = 0;

    for (i=0; viewbar_items[i].str != NULL; i++) {
	if (viewbar_items[i].flag == SAVE_ITEM ||
	    viewbar_items[i].flag == EDIT_ITEM) {
	    GtkWidget *w;

	    button = gtk_image_new();
	    gtk_image_set_from_stock(GTK_IMAGE(button), 
				     viewbar_items[i].icon, 
				     GTK_ICON_SIZE_MENU);
	    w = gtk_toolbar_insert_item(GTK_TOOLBAR(vwin->mbar),
					NULL, _(viewbar_items[i].str), NULL,
					button, viewbar_items[i].toolfunc, 
					vwin, pos);
	    g_object_set_data(G_OBJECT(w), "flag", 
			      GINT_TO_POINTER(viewbar_items[i].flag));
	    if (viewbar_items[i].flag == SAVE_ITEM) { 
		gtk_widget_set_sensitive(w, FALSE);
	    } 
	}
	if (viewbar_items[i].flag != GP_ITEM) {
	    pos++;
	}
    }
}

static gchar *make_viewer_title (int role, const char *fname)
{
    gchar *title = NULL;

    switch (role) {
    case GUI_HELP: 
	title = g_strdup(_("gretl: help")); break;
    case CLI_HELP:
	title = g_strdup(_("gretl: command syntax")); break;
    case GUI_HELP_EN: 
	title = g_strdup("gretl: help"); break;
    case CLI_HELP_EN:
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
#ifdef ENABLE_NLS
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
    const char *buf = gretl_print_get_buffer(prn);

#ifdef USE_GTKSOURCEVIEW
    if (vwin->role == VIEW_FUNC_CODE || vwin->role == EDIT_FUNC_CODE) {
	sourceview_insert_buffer(vwin, buf);
	return;
    }
#endif

    if (vwin->role == SCRIPT_OUT) {
	textview_set_text_colorized(vwin->w, buf);
    } else {
	textview_set_text(vwin->w, buf);
    }
}

static void viewer_add_close_button (windata_t *vwin)
{
    GtkWidget *button;

    button = gtk_button_new_with_label(_("Close"));
    gtk_box_pack_start(GTK_BOX(vwin->vbox), button, FALSE, TRUE, 0);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_file_viewer), vwin);
    gtk_widget_show(button);
}

windata_t *view_buffer (PRN *prn, int hsize, int vsize, 
			const char *title, int role, 
			gpointer data) 
{
    windata_t *vwin;

    if (title != NULL) {
	vwin = common_viewer_new(role, title, data, 1);
    } else {
	gchar *tmp = make_viewer_title(role, NULL);

	vwin = common_viewer_new(role, tmp, data, 1);
	g_free(tmp);
    }

    if (vwin == NULL) return NULL;

    viewer_box_config(vwin);

    if (role == VAR || role == VECM || role == SYSTEM) {
	/* special case: use a text-based menu bar */
	set_up_viewer_menu(vwin->dialog, vwin, 
			   (role == SYSTEM)? SYS_items : VAR_items);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
	if (role == SYSTEM) {
	    gretl_object_ref(data, GRETL_OBJ_SYS);
	    add_SYS_menu_items(vwin);
	} else {
	    gretl_object_ref(data, GRETL_OBJ_VAR);
	    add_VAR_menu_items(vwin, role == VECM);
	}
    } else if (role == VIEW_FUNC_CODE || role == EDIT_FUNC_CODE) {
	make_viewbar(vwin, 0);
    } else if (role != IMPORT) {
	make_viewbar(vwin, 1);
    }

#ifdef USE_GTKSOURCEVIEW
    if (role == VIEW_FUNC_CODE) {
	create_source(vwin, hsize, vsize, FALSE);
    } else if (role == EDIT_FUNC_CODE) {
	create_source(vwin, hsize, vsize, TRUE);
    } else {
	vwin->w = create_text(vwin->dialog, hsize, vsize, FALSE);
    }
#else
    vwin->w = create_text(vwin->dialog, hsize, vsize, 
			  role == EDIT_FUNC_CODE);
#endif

    text_table_setup(vwin->vbox, vwin->w);

    /* arrange for clean-up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    /* "Close" button */
    viewer_add_close_button(vwin);

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
	g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
			 G_CALLBACK(query_save_text), vwin);
	/* FIXME add callback for updating fn code */
    }

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
    cursor_to_top(vwin);

    return vwin;
}

windata_t *view_file (const char *filename, int editable, int del_file, 
		      int hsize, int vsize, int role)
{
    windata_t *vwin;
    FILE *fp;
    gchar *title = NULL;
    int doing_script = (role == EDIT_SCRIPT ||
			role == VIEW_SCRIPT ||
			role == VIEW_LOG);

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
			     NULL, !doing_script && role != CONSOLE);
    g_free(title);

    if (vwin == NULL) {
	return NULL;
    }

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);
    make_viewbar(vwin, (role == VIEW_DATA || role == CONSOLE));

#ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	create_source(vwin, hsize, vsize, editable);
    } else {
	vwin->w = create_text(vwin->dialog, hsize, vsize, editable);
    }
#else
    vwin->w = create_text(vwin->dialog, hsize, vsize, editable);
#endif

    text_table_setup(vwin->vbox, vwin->w);

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

    /* "Close" button */
    viewer_add_close_button(vwin);

#ifdef USE_GTKSOURCEVIEW
    if (doing_script || role == GR_PLOT) {
	sourceview_insert_file(vwin, filename);
    } else {
	textview_insert_file(vwin, filename);
    }
#else
    textview_insert_file(vwin, filename);
#endif

    /* grab the "changed" signal when editing a script */
    if (role == EDIT_SCRIPT) {
	attach_content_changed_signal(vwin);
    }

    /* catch some keystrokes */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (editable) {
	g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
    }

    /* alert for unsaved changes on exit */
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

    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
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
    int hsize = 80, vsize = 400;

    title = make_viewer_title(role, NULL);
    vwin = common_viewer_new(role, title, NULL, 0);
    g_free(title);

    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);
    set_up_viewer_menu(vwin->dialog, vwin, menu_items);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    vwin->w = create_text(vwin->dialog, hsize, vsize, FALSE);
    text_table_setup(vwin->vbox, vwin->w);

    /* "Close" button */
    viewer_add_close_button(vwin);

    /* grab content of the appropriate help file into a buffer
       and stick it onto vwin as data */
    g_file_get_contents(filename, &fbuf, NULL, NULL);
    vwin->data = fbuf;

    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (vwin->role == CLI_HELP || vwin->role == CLI_HELP_EN) {
	g_signal_connect(G_OBJECT(vwin->w), "button_press_event",
			 G_CALLBACK(help_popup_handler), 
			 vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
			 G_CALLBACK(catch_button_3), vwin->w);
    }	

    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->dialog);

    gtk_widget_grab_focus(vwin->w);

    return vwin;
}

void view_window_set_editable (windata_t *vwin)
{
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->w), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->w), TRUE);
    g_object_set_data(G_OBJECT(vwin->dialog), "vwin", vwin);
    g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
		     G_CALLBACK(query_save_text), vwin);
    vwin->role = EDIT_SCRIPT;
    add_edit_items_to_viewbar(vwin);
    attach_content_changed_signal(vwin);
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
    windata_t *vwin;

    vwin = common_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) {
	return NULL;
    }

    viewer_box_config(vwin); 

    /* add a menu bar */
    make_viewbar(vwin, 0);

    vwin->w = create_text(vwin->dialog, hsize, vsize, TRUE);
    text_table_setup(vwin->vbox, vwin->w);
    
    /* insert the buffer text */
    if (*pbuf) {
	GtkTextBuffer *tbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

	gtk_text_buffer_set_text(tbuf, *pbuf, -1);
    }
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    attach_content_changed_signal(vwin);

    /* alert for unsaved changes on exit */
    g_signal_connect(G_OBJECT(vwin->dialog), "delete_event",
		     G_CALLBACK(query_save_text), vwin);

    /* clean up when dialog is destroyed */
    g_signal_connect(G_OBJECT(vwin->dialog), "destroy", 
		     G_CALLBACK(free_windata), vwin);

    /* "Close" button */
    viewer_add_close_button(vwin);

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

    if (pmod->ci != MLE) {
	add_vars_to_plot_menu(vwin);
	add_model_dataset_items(vwin);
	if (latex_is_ok() && !pmod->errcode) {
	    add_model_tex_items(vwin);
	}
    }	

    if (pmod->ci != ARMA && pmod->ci != GARCH && 
	pmod->ci != NLS && pmod->ci != MLE &&
	pmod->ci != PANEL && pmod->ci != ARBOND) {
	add_dummies_to_plot_menu(vwin);
    }

    g_signal_connect(G_OBJECT(vwin->mbar), "button_press_event", 
		     G_CALLBACK(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    vwin->w = create_text(vwin->dialog, hsize, vsize, FALSE);
    text_table_setup(vwin->vbox, vwin->w);

    /* "Close" button */
    viewer_add_close_button(vwin);

    /* insert and then free the model results buffer */
    buf = gretl_print_get_buffer(prn);
    textview_set_text(vwin->w, buf);
    gretl_print_destroy(prn);

    /* attach shortcuts */
    g_signal_connect(G_OBJECT(vwin->dialog), "key_press_event", 
		     G_CALLBACK(catch_viewer_key), vwin);
    g_signal_connect(G_OBJECT(vwin->w), "button_press_event", 
		     G_CALLBACK(catch_button_3), vwin->w);

    /* don't allow deletion of model window when a model
       test dialog is active */
    g_signal_connect(G_OBJECT(vwin->dialog), "delete_event", 
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

    if (pmod->ci == MLE || pmod->ci == MPOLS) { /* FIXME? */
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
    flip(vwin->ifac, "/Analysis/coefficient covariance matrix", FALSE);
    add_x12_output_menu_item(vwin);
}

static void adjust_model_menu_state (windata_t *vwin, const MODEL *pmod)
{
    set_tests_menu_state(vwin->ifac, pmod);

    /* disallow saving an already-saved model */
    if (pmod->name != NULL) {
	model_save_state(vwin->ifac, FALSE);
    }

    if (pmod->ci == MLE) {
	/* some of this could be relaxed later */
	flip(vwin->ifac, "/Analysis", FALSE);
	flip(vwin->ifac, "/Graphs", FALSE);
    } else if (pmod->ci == ARMA && arma_by_x12a(pmod)) {
	arma_x12_menu_mod(vwin);
    } 

    if (dataset_is_panel(datainfo) && pmod->ci == OLS) {
	panel_heteroskedasticity_menu(vwin);
    }

    if (pmod->ci == ARBOND) {
	flip(vwin->ifac, "/Graphs/Fitted, actual plot", FALSE);
	flip(vwin->ifac, "/Analysis/Forecasts...", FALSE);
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
    int eqn_ok = command_ok_for_model(EQNPRINT, pmod->ci);
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

static void add_vars_to_plot_menu (windata_t *vwin)
{
    int i, j, varstart;
    GtkItemFactoryEntry varitem;
    const gchar *mpath[] = {
	N_("/Graphs/Residual plot"), 
	N_("/Graphs/Fitted, actual plot")
    };
    MODEL *pmod = vwin->data;
    char tmp[16];

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

	if (pmod->ci == ARBOND) {
	    break;
	}

	if (pmod->ci == ARMA || pmod->ci == NLS || pmod->ci == GARCH ||
	    pmod->ci == PANEL) { 
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

    if (code != VAR_VCV) {
	h = default_VAR_horizon(datainfo);
	title = g_strdup_printf("gretl: %s", 
				(code == VAR_IRF)? _("impulse responses") :
				_("variance decompositions"));
	err = checks_dialog(title, NULL, 0, NULL, 0, NULL,
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

static int impulse_response_setup (int *horizon, int *bootstrap)
{
    gchar *title;
    int h = default_VAR_horizon(datainfo);
    const char *impulse_opts[] = {
	N_("include bootstrap confidence interval")
    };
    static int active[] = { 0 };
    int err;

    title = g_strdup_printf("gretl: %s", _("impulse responses"));

    err = checks_dialog(title, 
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
	*bootstrap = active[0];
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

    if (impulse_response_setup(&horizon, &bootstrap) < 0) {
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

static void multiple_irf_plot_call (gpointer p, guint vecm, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    const double **vZ = NULL;
    int err;

    if (impulse_response_setup(&horizon, &bootstrap) < 0) {
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
	t1 = var->order + var->ecm;
	pre_n = 0;
	dyn_ok = 0;
    }

    /* FIXME pre_n with static fcast? */

    resp = forecast_dialog(t1, t1, &t1,
			   t1, t2, &t2,
			   0, premax, &pre_n,
			   dyn_ok, NULL);
    if (resp < 0) {
	return;
    }

    if (resp == 1) {
	opt = OPT_D;
    } else if (resp == 2) {
	opt = OPT_S;
    }

    /* FIXME dating here */

    fr = get_VAR_forecast(var, i, t1 - pre_n, t1, t2, (const double **) Z, 
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
    VAR_ARCH_TEST,
    VAR_NORMALITY_TEST,
    VAR_RESTRICT
};

static void VAR_test_call (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    char title[72];
    PRN *prn;
    int order = 0;
    int err;

    if (bufopen(&prn)) {
	return;
    }

    if (code == VAR_AUTOCORR_TEST || code == VAR_ARCH_TEST) {
	order = default_lag_order(datainfo);
	set_window_busy(vwin);
	err = spin_dialog((code == VAR_AUTOCORR_TEST)?
			  _("gretl: autocorrelation") :
			  _("gretl: ARCH test"),
			  &order, _("Lag order for test:"),
			  1, datainfo->n / 2, LMTEST);
	unset_window_busy(vwin);
	if (err < 0) {
	    gretl_print_destroy(prn);
	    return;
	}
    }	

    if (code == VAR_AUTOCORR_TEST) {
	strcpy(title, _("gretl: LM test (autocorrelation)"));
	err = gretl_VAR_autocorrelation_test(var, order, 
					     &Z, datainfo, 
					     prn);
    } else if (code == VAR_ARCH_TEST) {
	strcpy(title, _("gretl: ARCH test"));
	err = gretl_VAR_arch_test(var, order, &Z, datainfo, prn);
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

static void VAR_roots_plot_call (gpointer p, guint vecm, GtkWidget *w)
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

static void VAR_resid_plot_call (gpointer p, guint vecm, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_residual_plot(var, datainfo);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void VAR_resid_mplot_call (gpointer p, guint vecm, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_residual_mplot(var, datainfo);
    
    if (err) {
	errbox(_("gnuplot command failed"));
    } else {
	register_graph();
    }
}

static void add_VAR_menu_items (windata_t *vwin, int vecm)
{
    GtkItemFactoryEntry varitem;

    const gchar *tpath = N_("/Tests");
    const gchar *gpath = N_("/Graphs");
    const gchar *mpath = N_("/Analysis");
    const gchar *fpath = N_("/Analysis/Forecasts");
    const gchar *dpath = N_("/Save");

    GRETL_VAR *var = NULL;
    int neqns, vtarg, vshock;
    char tmp[32];
    int i, j;

    var = (GRETL_VAR *) vwin->data;
    neqns = gretl_VAR_get_n_equations(var);

    varitem.accelerator = NULL;
    varitem.callback = NULL;
    varitem.callback_action = 0;
    varitem.item_type = "<Branch>";

    varitem.path = g_strdup(_("/_Tests"));
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    varitem.path = g_strdup(_("/_Analysis"));
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
				   _("Autocorrelation"));
    varitem.callback = VAR_test_call;
    varitem.callback_action = VAR_AUTOCORR_TEST;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* ARCH tests */
    varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				   _("ARCH"));
    varitem.callback = VAR_test_call;
    varitem.callback_action = VAR_ARCH_TEST;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* multivariate normality test */
    varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				   _("Normality of residuals"));
    varitem.callback = VAR_test_call;
    varitem.callback_action = VAR_NORMALITY_TEST;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    if (vecm) {
	/* linear restrictions on cointegrating relations */
	varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				       _("Linear restrictions"));
	varitem.callback = gretl_callback;
	varitem.callback_action = RESTRICT;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);
    } else {
	/* regular VAR: omit exogenous variables test */
	int err, *exolist;

	exolist = gretl_VAR_get_exo_list(var, &err);
	if (exolist != NULL) {
	    varitem.path = g_strdup_printf("%s/%s", _(tpath), 
				       _("Omit exogenous variables..."));
	    varitem.callback = selector_callback;
	    varitem.callback_action = VAROMIT;
	    varitem.item_type = NULL;
	    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	    g_free(varitem.path);
	    free(exolist);
	}	    
    }

    /* cross-equation VCV */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), 
				   _("Cross-equation covariance matrix"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_VCV;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* impulse response printout */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), _("Impulse responses"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_IRF;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);    

    /* variance decomp printout */
    varitem.path = g_strdup_printf("%s/%s", _(mpath), 
				   _("Forecast variance decomposition"));
    varitem.callback = VAR_model_data_callback;
    varitem.callback_action = VAR_DECOMP;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path); 

    if (neqns <= 6) {
	/* separate residual plot */
	varitem.path = g_strdup_printf("%s/%s", _(gpath), _("Residual plots"));
	varitem.callback = VAR_resid_mplot_call;
	varitem.callback_action = vecm;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);
    }

    /* combined residual plot */
    varitem.path = g_strdup_printf("%s/%s", _(gpath), _("Combined residual plot"));
    varitem.callback = VAR_resid_plot_call;
    varitem.callback_action = vecm;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    /* VAR inverse roots */
    varitem.path = g_strdup_printf("%s/%s", _(gpath), _("VAR inverse roots"));
    varitem.callback = VAR_roots_plot_call;
    varitem.callback_action = vecm;
    varitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
    g_free(varitem.path);

    if (neqns <= 4) {
	/* Multiple IRFs */
	varitem.path = g_strdup_printf("%s/%s", _(gpath), _("Impulse responses (combined)"));
	varitem.callback = multiple_irf_plot_call;
	varitem.callback_action = vecm;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);
    }

    for (i=0; i<neqns; i++) {
	char maj[64], min[32];
	int dv;

	/* forecast items */
	dv = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, datainfo->varname[dv]);
	varitem.path = g_strdup_printf("%s/%s", _(fpath), tmp);
	varitem.callback = VAR_forecast_callback;
	varitem.callback_action = i;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	/* save resids items */
	varitem.path = g_strdup_printf("%s/%s %d", _(dpath), 
				       _("Residuals from equation"), i + 1);
	varitem.callback = VAR_resid_callback;
	varitem.callback_action = i;
	varitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &varitem, vwin, 1);
	g_free(varitem.path);

	/* impulse response plots: make branch for target */
	vtarg = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, datainfo->varname[vtarg]);
	sprintf(maj, _("Response of %s"), tmp);

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
	    g_object_set_data(G_OBJECT(w), "targ", GINT_TO_POINTER(i));
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

enum {
    SYS_RESTRICT,
    SYS_NORMALITY
};

static void SYS_test_call (gpointer p, guint code, GtkWidget *w)
{
    windata_t *vwin = (windata_t *) p;
    gretl_equation_system *sys;
    char title[72];
    PRN *prn;
    int err;

    sys = (gretl_equation_system *) vwin->data;
    if (sys == NULL) {
	/* error message? */
	return;
    }    

    if (bufopen(&prn)) {
	return;
    }

    if (code == SYS_NORMALITY) {
	sprintf(title, "gretl: %s", _("Test for normality of residual"));
	err = system_normality_test(sys, prn);
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

static void add_SYS_menu_items (windata_t *vwin)
{
    GtkItemFactoryEntry sysitem;
    const gchar *tpath = N_("/Tests");
    const gchar *dpath = N_("/Save");
    gretl_equation_system *sys;
    int i, neqns;

    sys = (gretl_equation_system *) vwin->data;
    if (sys == NULL) {
	return;
    }

    neqns = sys->n_equations;

    sysitem.accelerator = NULL;
    sysitem.callback = NULL;
    sysitem.callback_action = 0;
    sysitem.item_type = "<Branch>";

#if 0
    /* model data menu path */
    sysitem.path = g_strdup(_("/_Analysis"));
    gtk_item_factory_create_item(vwin->ifac, &sysitem, vwin, 1);
    g_free(sysitem.path);
#endif

    /* add to dataset menu path */
    sysitem.path = g_strdup(_("/_Save"));
    gtk_item_factory_create_item(vwin->ifac, &sysitem, vwin, 1);
    g_free(sysitem.path);

    /* multivariate normality test */
    sysitem.path = g_strdup_printf("%s/%s", _(tpath), 
				   _("Normality of residuals"));
    sysitem.callback = SYS_test_call;
    sysitem.callback_action = SYS_NORMALITY;
    sysitem.item_type = NULL;
    gtk_item_factory_create_item(vwin->ifac, &sysitem, vwin, 1);
    g_free(sysitem.path);  

    for (i=0; i<neqns; i++) {
	/* save resids items */
	sysitem.path = g_strdup_printf("%s/%s %d", _(dpath), 
				       _("Residuals from equation"), i + 1);
	sysitem.callback = SYS_resid_callback;
	sysitem.callback_action = i;
	sysitem.item_type = NULL;
	gtk_item_factory_create_item(vwin->ifac, &sysitem, vwin, 1);
	g_free(sysitem.path);
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

    if (pmod->ci == MLE || pmod->ci == MPOLS) {
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

#ifndef G_OS_WIN32

int browser_open (const char *url)
{
# if defined(USE_GNOME)
    gnome_url_show(url, NULL); 
# elif defined(OSX_BUILD)
    osx_open_url(url);
# else
    int err;
    char ns_cmd[256];

    sprintf(ns_cmd, "%s -remote \"openURLNewWindow(%s)\"", Browser, url);
    err = gretl_spawn(ns_cmd);
    if (err) {
	gretl_fork(Browser, url);
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
	errbox(_("Please open a data file first"));
	return;
    }

    build_path(Rprofile, paths.userdir, "gretl.Rprofile", NULL);
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

    build_path(Rdata, paths.userdir, "Rdata.tmp", NULL);

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

	build_path(Rtmp, paths.userdir, "Rtmp", NULL);
	fq = fopen(Rtmp, "w");
	if (fq != NULL) {
	    fprintf(fq, "gretldata <- read.table(\"%s\")\n", Rdata);
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
