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
#include "forecast.h"
#include "objstack.h"
#include "gretl_xml.h"
#include "gretl_func.h"
#include "system.h"
#include "matrix_extra.h"
#include "bootstrap.h"
#include "gretl_foreign.h"

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
#include "fileselect.h"
#include "toolbar.h"
#include "winstack.h"
#include "fnsave.h"
#include "datawiz.h"
#include "selector.h"

#ifdef G_OS_WIN32
# include <windows.h>
# include "gretlwin32.h"
#endif

char *storelist = NULL;

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkActionEntry *items, const gchar *ui_info);
static void set_up_model_view_menu (GtkWidget *window, windata_t *vwin);
static void auto_save_script (windata_t *vwin);
static void add_system_menu_items (windata_t *vwin, int vecm);
static void add_x12_output_menu_item (windata_t *vwin);
static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data);
static void model_copy_callback (GtkAction *action, gpointer p);

static void close_model (GtkAction *action, gpointer data)
{
    windata_t *vwin = (windata_t *) data;

    if (window_is_busy(vwin)) {
	maybe_raise_dialog();
    } else {
	gtk_widget_destroy(vwin->main);
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

int latex_is_ok (void)
{
    static int latex_ok = -1; 
  
    if (latex_ok == -1) {
	latex_ok = check_for_prog(latex);
    }

    return latex_ok;
}

static void model_output_save (GtkAction *action, gpointer p)
{
    copy_format_dialog((windata_t *) p, W_SAVE);    
}

static gretlopt tex_eqn_opt;

static void set_tex_eqn_opt (GtkRadioAction *action)
{
    int v = gtk_radio_action_get_current_value(action);

    tex_eqn_opt = (v)? OPT_T : OPT_NONE;
}

gretlopt get_tex_eqn_opt (void)
{
    return tex_eqn_opt;
}

static void model_revise_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    
    selector_from_model(vwin->data, vwin->role);
}

static GtkActionEntry model_items[] = {
    { "File", NULL, N_("_File"), NULL, NULL, NULL },
    { "SaveAs", GTK_STOCK_SAVE_AS, N_("_Save as..."), NULL, NULL, G_CALLBACK(model_output_save) },
    { "SaveAsIcon", NULL, N_("Save to session as _icon"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
    { "SaveAndClose", NULL, N_("Save as icon and cl_ose"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
#ifdef NATIVE_PRINTING
    { "Print", GTK_STOCK_PRINT, N_("_Print..."), NULL, NULL, G_CALLBACK(window_print) },
#endif
    { "Close", GTK_STOCK_CLOSE, N_("_Close"), NULL, NULL, G_CALLBACK(close_model) },
    { "Edit", NULL, N_("_Edit"), NULL, NULL, NULL },    
    { "Copy", GTK_STOCK_COPY, N_("_Copy"), NULL, NULL, G_CALLBACK(model_copy_callback) },
    { "Revise", GTK_STOCK_EDIT, N_("_Modify model..."), NULL, NULL, 
      G_CALLBACK(model_revise_callback) },
    { "Tests", NULL, N_("_Tests"), NULL, NULL, NULL },    
    { "Save", NULL, N_("_Save"), NULL, NULL, NULL },    
    { "Graphs", NULL, N_("_Graphs"), NULL, NULL, NULL },    
    { "ResidPlot", NULL, N_("_Residual plot"), NULL, NULL, NULL },    
    { "FittedActualPlot", NULL, N_("_Fitted, actual plot"), NULL, NULL, NULL },    
    { "Analysis", NULL, N_("_Analysis"), NULL, NULL, NULL }, 
    { "DisplayAFR", NULL, N_("_Display actual, fitted, residual"), NULL, NULL, 
      G_CALLBACK(display_fit_resid) },    
    { "Forecasts", NULL, N_("_Forecasts..."), NULL, NULL, G_CALLBACK(gui_do_forecast) },    
    { "ConfIntervals", NULL, N_("_Confidence intervals for coefficients"), NULL, NULL, 
      G_CALLBACK(do_coeff_intervals) },    
    { "ConfEllipse", NULL, N_("Confidence _ellipse..."), NULL, NULL, G_CALLBACK(selector_callback) },    
    { "Covariance", NULL, N_("Coefficient covariance _matrix"), NULL, NULL, G_CALLBACK(do_outcovmx) },    
    { "ANOVA", NULL, N_("_ANOVA"), NULL, NULL, G_CALLBACK(do_anova) },    
    { "Bootstrap", NULL, N_("_Bootstrap..."), NULL, NULL, G_CALLBACK(do_bootstrap) }    
};

static GtkActionEntry model_test_items[] = {
    { "omit", NULL, N_("_Omit variables"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "add", NULL, N_("_Add variables"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "coeffsum", NULL, N_("_Sum of coefficients"), NULL, NULL, G_CALLBACK(selector_callback) },
    { "restrict", NULL, N_("_Linear restrictions"), NULL, NULL, G_CALLBACK(gretl_callback) },
    { "lmtest:s", NULL, N_("Non-linearity (s_quares)"), NULL, NULL, G_CALLBACK(do_lmtest) },
    { "lmtest:l", NULL, N_("Non-linearity (_logs)"), NULL, NULL, G_CALLBACK(do_lmtest) },
    { "reset", NULL, N_("_Ramsey's RESET"), NULL, NULL, G_CALLBACK(do_reset) },
    { "Hsk", NULL, N_("_Heteroskedasticity"), NULL, NULL, NULL },    
    { "normtest", NULL, N_("_Normality of residual"), NULL, NULL, G_CALLBACK(do_resid_freq) },
    { "leverage", NULL, N_("_Influential observations"), NULL, NULL, G_CALLBACK(do_leverage) },
    { "chow", NULL, N_("_Chow test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },    
    { "vif", NULL, N_("_Collinearity"), NULL, NULL, G_CALLBACK(do_vif) },
    { "lmtest:a", NULL, N_("_Autocorrelation"), NULL, NULL, G_CALLBACK(do_autocorr) },
    { "dwpval", NULL, N_("_Durbin-Watson p-value"), NULL, NULL, G_CALLBACK(do_dwpval) },
    { "lmtest:h", NULL, N_("A_RCH"), NULL, NULL, G_CALLBACK(do_arch) },
    { "qlrtest", NULL, N_("_QLR test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "cusum", NULL, N_("_CUSUM test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "cusum:r", NULL, N_("CUSUM_SQ test"), NULL, NULL, G_CALLBACK(do_chow_cusum) },
    { "hausman", NULL, N_("_Panel diagnostics"), NULL, NULL, G_CALLBACK(do_panel_tests) }
};

static GtkActionEntry base_hsk_items[] = {
    { "White", NULL, N_("White's test"), NULL, NULL, G_CALLBACK(do_lmtest) },
    { "WhiteSquares", NULL, N_("White's test (squares only)"), NULL, NULL, G_CALLBACK(do_lmtest) },
    { "BreuschPagan", NULL, "Breusch-Pagan", NULL, NULL, G_CALLBACK(do_lmtest) },
    { "Koenker", NULL, "Koenker", NULL, NULL, G_CALLBACK(do_lmtest) }
};

static GtkActionEntry panel_hsk_items[] = {
    { "White", NULL, N_("White's test"), NULL, NULL, G_CALLBACK(do_lmtest) },
    { "Groupwise", NULL, N_("_groupwise"), NULL, NULL, G_CALLBACK(do_lmtest) }
};

const gchar *model_tex_ui = 
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='LaTeX'>"
    "      <menu action='TeXView'>"
    "        <menuitem action='TabView'/>"
    "        <menuitem action='EqnView'/>"
    "      </menu>"
    "      <menu action='TeXCopy'>"
    "        <menuitem action='TabCopy'/>"
    "        <menuitem action='EqnCopy'/>"
    "      </menu>"
    "      <menu action='TeXSave'>"
    "        <menuitem action='TabSave'/>"
    "        <menuitem action='EqnSave'/>"
    "      </menu>"
    "      <menu action='EqnOpts'>"
    "        <menuitem action='TeXstderrs'/>"
    "        <menuitem action='TeXtratios'/>"
    "      </menu>"
    "      <menuitem action='TabOpts'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

const gchar *missing_tex_ui =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='LaTeX'>"
    "      <menuitem action='notex'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

static GtkActionEntry model_tex_items[] = {
    { "LaTeX",   NULL, N_("_LaTeX"), NULL, NULL, NULL },      
    { "TeXView", NULL, N_("_View"), NULL, NULL, NULL },      
    { "TabView", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_view) },      
    { "EqnView", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_view) },      
    { "TeXCopy", NULL, N_("_Copy"), NULL, NULL, NULL }, 
    { "TabCopy", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_copy) },      
    { "EqnCopy", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_copy) },      
    { "TeXSave", NULL, N_("_Save"), NULL, NULL, NULL }, 
    { "TabSave", NULL, N_("_Tabular"), NULL, NULL, G_CALLBACK(model_tex_save) },      
    { "EqnSave", NULL, N_("_Equation"), NULL, NULL, G_CALLBACK(model_tex_save) },      
    { "EqnOpts", NULL, N_("_Equation options"), NULL, NULL, NULL }, 
    { "TabOpts", NULL, N_("_Tabular options..."), NULL, NULL, G_CALLBACK(tex_format_dialog) }
};

static GtkRadioActionEntry tex_eqn_items[] = {
    { "TeXstderrs", NULL, N_("Show _standard errors"), NULL, NULL, 0 },
    { "TeXtratios", NULL, N_("Show _t-ratios"), NULL, NULL, 1 },
};

static GtkActionEntry missing_tex_items[] = {
    { "LaTeX", NULL, N_("_LaTeX"), NULL, NULL, NULL },    
    { "notex", NULL, "No TeX", NULL, NULL, G_CALLBACK(dummy_call) }
};

static GtkActionEntry system_items[] = {
    { "File", NULL, N_("_File"), NULL, NULL, NULL },      
    { "SaveAs", GTK_STOCK_SAVE_AS, N_("_Save as..."), NULL, NULL, G_CALLBACK(model_output_save) },      
    { "SaveAsIcon", NULL, N_("Save to session as _icon"), NULL, NULL, G_CALLBACK(model_add_as_icon) },      
    { "SaveAndClose", NULL, N_("Save as icon and cl_ose"), NULL, NULL, G_CALLBACK(model_add_as_icon) },
#ifdef NATIVE_PRINTING
    { "Print", GTK_STOCK_PRINT, N_("_Print..."), NULL, NULL, G_CALLBACK(window_print) },
#endif
    { "Close", GTK_STOCK_CLOSE, N_("_Close"), NULL, NULL, G_CALLBACK(close_model) },
    { "Edit", NULL, N_("_Edit"), NULL, NULL, NULL },      
    { "Copy", GTK_STOCK_COPY, N_("_Copy"), NULL, NULL, G_CALLBACK(model_copy_callback) }, 
    { "Revise", GTK_STOCK_EDIT, N_("_Revise specification..."), NULL, NULL, 
      G_CALLBACK(model_revise_callback) }, 
    { "Save", NULL, N_("_Save"), NULL, NULL, NULL },    
    { "Tests", NULL, N_("_Tests"), NULL, NULL, NULL },    
    { "Graphs", NULL, N_("_Graphs"), NULL, NULL, NULL },    
    { "Analysis", NULL, N_("_Analysis"), NULL, NULL, NULL },  
    { "Forecasts", NULL, N_("_Forecasts"), NULL, NULL, NULL },  
    { NULL, NULL, NULL, NULL, NULL, NULL }
};

static GtkActionEntry sys_tex_items[] = {
    { "LaTeX", NULL, N_("_LaTeX"), NULL, NULL, NULL },  
    { "TeXView", NULL, N_("_View"), NULL, NULL, G_CALLBACK(model_tex_view) },      
    { "TeXCopy", NULL, N_("_Copy"), NULL, NULL, G_CALLBACK(model_tex_copy) },      
    { "TeXSave", NULL, N_("_Save"), NULL, NULL, G_CALLBACK(model_tex_save) },
};

static const gchar *sys_ui =
    "<ui>"
    "  <menubar name='MenuBar'>"
    "    <menu action='File'>"
    "      <menuitem action='SaveAs'/>"
    "      <menuitem action='SaveAsIcon'/>"
    "      <menuitem action='SaveAndClose'/>"
#ifdef NATIVE_PRINTING
    "      <menuitem action='Print'/>"
#endif
    "      <menuitem action='Close'/>"
    "    </menu>"
    "    <menu action='Edit'>"
    "      <menuitem action='Copy'/>"
    "      <menuitem action='Revise'/>"
    "    </menu>"
    "    <menu action='Tests'/>"
    "    <menu action='Save'/>"
    "    <menu action='Graphs'/>"
    "    <menu action='Analysis'>"
    "      <menu action='Forecasts'/>"
    "    </menu>"
    "  </menubar>"
    "</ui>";

static void model_copy_callback (GtkAction *action, gpointer p)
{
    copy_format_dialog((windata_t *) p, W_COPY);
}

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

static void delete_file (GtkWidget *widget, char *fname) 
{
    remove(fname);
    g_free(fname);
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
    int (*read_window_text) (windata_t *, const DATAINFO *, int (*)());
    void *handle;

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
    GtkTextBuffer *buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (gtk_text_buffer_get_selection_bounds(buf, NULL, NULL)) {
	window_copy(vwin, GRETL_FORMAT_SELECTION);
    } else if (MULTI_FORMAT_ENABLED(vwin->role)) {
	window_copy(vwin, GRETL_FORMAT_RTF);
    } else {
	window_copy(vwin, GRETL_FORMAT_TXT);
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
    editing = gtk_text_view_get_editable(GTK_TEXT_VIEW(vwin->text));

    if (mods & GDK_CONTROL_MASK) {
	if (key->keyval == GDK_f) {
	    /* Ctrl-F: find */
	    text_find(NULL, vwin);
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
		if (vwin->role == EDIT_SCRIPT && 
		    (vwin->flags & VWIN_CONTENT_CHANGED)) {
		    if (query_save_text(NULL, NULL, vwin) == FALSE) {
			gtk_widget_destroy(vwin->main);
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
	model_add_as_icon(NULL, vwin);
    } else if (key->keyval == GDK_w) {
	if (mods & GDK_CONTROL_MASK) {
	    gtk_widget_destroy(w);
	    return TRUE;
	}	
    }

#if defined(HAVE_FLITE) || defined(G_OS_WIN32)
    /* respond to 'a' and 'x', but not if Alt-modified */
    else if (!(mods & GDK_MOD1_MASK)) {
	if (key->keyval == GDK_a) {
	    audio_render_window(vwin, AUDIO_TEXT);
	} else if (key->keyval == GDK_x) {
	    stop_talking();
	}
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

void mark_dataset_as_modified (void)
{
    data_status |= MODIFIED_DATA;
    set_sample_label(datainfo);
}

static void gui_record_data_opening (const char *fname, const int *list)
{
    const char *recname = (fname != NULL)? fname : paths.datfile;

    if (haschar(' ', recname) >= 0) {
	gretl_command_sprintf("open \"%s\"", recname);
    } else {
	gretl_command_sprintf("open %s", recname);
    }

    if (list != NULL && list[0] == 3) {
	/* record spreadsheet parameters */
	char parm[32];

	if (list[1] != 1) {
	    sprintf(parm, " --sheet=%d", list[1]);
	    gretl_command_strcat(parm);
	}
	if (list[2] != 0) {
	    sprintf(parm, " --coloffset=%d", list[2]);
	    gretl_command_strcat(parm);
	}
	if (list[3] != 0) {
	    sprintf(parm, " --rowoffset=%d", list[3]);
	    gretl_command_strcat(parm);
	}
    }

    check_and_record_command();

    if (*paths.datfile != '\0') {
	char tmp[FILENAME_MAX];

	strcpy(tmp, paths.datfile);
	mkfilelist(FILE_LIST_DATA, tmp);
    }
}

#define file_opened(f) (f == DATAFILE_OPENED || \
	                f == OPENED_VIA_CLI || \
                        f == OPENED_VIA_SESSION)

static void real_register_data (int flag, const char *user_fname, 
				const int *list)
{    
    /* basic accounting */
    data_status |= HAVE_DATA;
    orig_vars = datainfo->v;

    /* set appropriate data_status bits */
    if (file_opened(flag)) {
	if (!(data_status & IMPORT_DATA)) {
	    /* we opened a native data file */
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
    } else {
	/* we modified the current dataset somehow */
	data_status |= GUI_DATA;
	mark_dataset_as_modified();
    }

    /* sync main window with datafile */
    if (mdata != NULL) {
	populate_varlist();
	set_sample_label(datainfo);
	main_menubar_state(TRUE);
	session_menu_state(TRUE);
    }

    /* Record the opening of the data file in the GUI recent files
       list and command log; note that we don't do this if the file
       was opened via script or console, or if it was opened as a
       side effect of re-opening a saved session.
    */
    if (mdata != NULL && flag == DATAFILE_OPENED) {
	gui_record_data_opening(user_fname, list);
    } 

    if (mdata != NULL) {
	/* focus the data window */
	gtk_widget_grab_focus(mdata->listbox);
	/* invalidate "remove extra obs" menu item */
	drop_obs_state(FALSE);
    }
}

void register_data (int flag)
{
    real_register_data(flag, NULL, NULL);
}

void register_startup_data (const char *fname)
{
    real_register_data(DATAFILE_OPENED, fname, NULL);
}

static void finalize_data_open (const char *fname, int ftype,
				int import, int append, 
				const int *plist)
{
    if (import) {
	if (ftype == GRETL_CSV || ftype == GRETL_DTA || ftype == GRETL_SAV) {
	    maybe_display_string_table();
	}
	data_status |= IMPORT_DATA;
    }

    if (append) {
	register_data(DATA_APPENDED);
    } else {
	if (fname != paths.datfile) {
	    strcpy(paths.datfile, fname);
	}
	real_register_data(DATAFILE_OPENED, NULL, plist);
    }  

    if (import && !dataset_is_time_series(datainfo) && 
	!dataset_is_panel(datainfo) && mdata != NULL) {
	int resp;

	resp = yes_no_dialog(_("gretl: open data"),
			     _("The imported data have been interpreted as undated\n"
			       "(cross-sectional).  Do you want to give the data a\n"
			       "time-series or panel interpretation?"),
			     0);
	if (resp == GRETL_YES) {
	    data_structure_dialog();
	}
    }    
}

static int datafile_missing (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    int err = 0;

    if (fp == NULL) {
	delete_from_filelist(FILE_LIST_DATA, fname);
	file_read_errbox(fname);
	err = E_FOPEN;
    } else {
	fclose(fp);
    }

    return err;
}

/* below: get data of a sort that requires an import plugin */

int get_imported_data (char *fname, int ftype, int append)
{
    void *handle;
    PRN *errprn;
    const char *errbuf;
    int list[4] = {3, 1, 0, 0};
    int *plist = NULL;
    int (*ss_importer) (const char *, int *, char *, double ***, DATAINFO *, 
			gretlopt, PRN *);
    int (*misc_importer) (const char *, double ***, DATAINFO *, 
			  gretlopt, PRN *);
    int err = 0;

    if (datafile_missing(fname)) {
	return E_FOPEN;
    }

    ss_importer = NULL;
    misc_importer = NULL;

    if (ftype == GRETL_XLS) {
	ss_importer = gui_get_plugin_function("xls_get_data",
					      &handle);
	plist = list;
    } else if (ftype == GRETL_GNUMERIC) {
	ss_importer = gui_get_plugin_function("gnumeric_get_data",
					      &handle);
	plist = list;
    } else if (ftype == GRETL_ODS) {
	ss_importer = gui_get_plugin_function("ods_get_data",
					      &handle);
	plist = list;
    } else if (ftype == GRETL_DTA) {
	misc_importer = gui_get_plugin_function("dta_get_data",
						&handle);
    } else if (ftype == GRETL_SAV) {
	misc_importer = gui_get_plugin_function("sav_get_data",
						&handle);
    } else if (ftype == GRETL_JMULTI) {
	misc_importer = gui_get_plugin_function("jmulti_get_data",
						&handle);
    } else if (ftype == GRETL_WF1) {
	misc_importer = gui_get_plugin_function("wf1_get_data",
						&handle);
    } else {
	errbox(_("Unrecognized data type"));
	return 1;
    }

    if (ss_importer == NULL && misc_importer == NULL) {
        return 1;
    }

    if (bufopen(&errprn)) {
	close_plugin(handle);
	return 1;
    }

    if (SPREADSHEET_IMPORT(ftype)) {
	err = (*ss_importer)(fname, plist, NULL, &Z, datainfo, 
			     OPT_G, errprn);
    } else {
	err = (*misc_importer)(fname, &Z, datainfo, OPT_G, errprn);
    }
	
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
	    gui_errmsg(err);
	} 
	delete_from_filelist(FILE_LIST_DATA, fname);
    } else if (errbuf != NULL && *errbuf != '\0') {
	infobox(errbuf);
    }

    gretl_print_destroy(errprn);

    if (!err) {
	finalize_data_open(fname, ftype, 1, append, plist);
    }

    return err;
}

/* get "CSV" (or more generally, ASCII) data or GNU octave data:
   plugin is not required
*/

static int get_csv_data (char *fname, int ftype, int append)
{
    PRN *prn;
    gchar *title;
    int err = 0;

    if (datafile_missing(fname)) {
	return E_FOPEN;
    }

    if (bufopen(&prn)) {
	return 1;
    }

    if (ftype == GRETL_OCTAVE) {
	err = import_other(fname, ftype, &Z, datainfo, OPT_NONE, prn);
	title = g_strdup_printf(_("gretl: import %s data"), "Octave");
    } else {
	err = import_csv(fname, &Z, datainfo, OPT_NONE, prn);
	title = g_strdup_printf(_("gretl: import %s data"), "CSV");
    }

    /* show details regarding the import */
    view_buffer(prn, 78, 350, title, IMPORT, NULL); 
    g_free(title);

    if (err) {
	delete_from_filelist(FILE_LIST_DATA, fname);
    } else {
	finalize_data_open(fname, ftype, 1, append, NULL);
    }

    return err;
}

static int get_native_data (char *fname, int ftype, int append, 
			    windata_t *fwin)
{
    /* we'll send any output to stderr */
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);
    int err;

    if (ftype == GRETL_XML_DATA) {
	err = gretl_read_gdt(fname, &paths, &Z, datainfo, 
			     OPT_P, prn);
    } else {
	err = gretl_get_data(fname, &paths, &Z, datainfo, 
			     OPT_NONE, prn);
    }

    gretl_print_destroy(prn);

    if (fwin != NULL) {
	/* close the files browser window that launched the query */
	gtk_widget_destroy(fwin->main);
    }    

    if (err) {
	if (err == E_FOPEN) {
	    file_read_errbox(tryfile);
	} else {
	    gui_errmsg(err);
	}
	delete_from_filelist(FILE_LIST_DATA, fname);
    } else {
	finalize_data_open(fname, ftype, 0, append, NULL);
    } 

    return err;
}

#define APPENDING(action) (action == APPEND_DATA || \
                           action == APPEND_CSV || \
                           action == APPEND_GNUMERIC || \
                           action == APPEND_XLS || \
                           action == APPEND_ODS || \
                           action == APPEND_ASCII || \
                           action == APPEND_WF1 || \
                           action == APPEND_DTA || \
	                   action == APPEND_SAV || \
                           action == APPEND_JMULTI)

void do_open_data (windata_t *fwin, int code)
{
    int append = APPENDING(code);
    gint ftype;

    gretl_error_clear();

    if (code == OPEN_CSV || code == APPEND_CSV || code == OPEN_ASCII ||
	code == APPEND_ASCII) {
	ftype = GRETL_CSV;
    } else if (code == OPEN_GNUMERIC || code == APPEND_GNUMERIC) {
	ftype = GRETL_GNUMERIC;
    } else if (code == OPEN_ODS || code == APPEND_ODS) {
	ftype = GRETL_ODS;
    } else if (code == OPEN_XLS || code == APPEND_XLS) {
	ftype = GRETL_XLS;
    } else if (code == OPEN_OCTAVE || code == APPEND_OCTAVE) {
	ftype = GRETL_OCTAVE;
    } else if (code == OPEN_WF1 || code == APPEND_WF1) {
	ftype = GRETL_WF1;
    } else if (code == OPEN_DTA || code == APPEND_DTA) {
	ftype = GRETL_DTA;
    } else if (code == OPEN_SAV || code == APPEND_SAV) {
	ftype = GRETL_SAV;
    } else if (code == OPEN_JMULTI || code == APPEND_JMULTI) {
	ftype = GRETL_JMULTI;
    } else {
	/* no filetype specified: have to guess */
	ftype = detect_filetype(tryfile, &paths);
    }

    /* destroy the current data set, etc., unless we're explicitly appending */
    if (!append) {
	close_session(NULL, &Z, datainfo, OPT_NONE); /* FIXME opt? */
    }

    if (ftype == GRETL_CSV || ftype == GRETL_OCTAVE) {
	get_csv_data(tryfile, ftype, append);
    } else if (SPREADSHEET_IMPORT(ftype) || OTHER_IMPORT(ftype)) {
	get_imported_data(tryfile, ftype, append);
    } else { 
	get_native_data(tryfile, ftype, append, fwin);
    }
}

/* give user choice of not opening selected datafile, if there's
   already a datafile open */

void verify_open_data (windata_t *vwin, int code)
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

    do_open_data(vwin, code);
}

/* give user choice of not opening session file, if there's already a
   datafile open */

void verify_open_session (void)
{
    if (!gretl_is_pkzip_file(tryfile)) {
	/* not a new-style zipped session file */
	do_open_script(EDIT_SCRIPT);
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

void mark_vwin_content_changed (windata_t *vwin) 
{
    if (vwin->active_var == 0) {
	GtkWidget *w = g_object_get_data(G_OBJECT(vwin->mbar), "save_button");

	if (w != NULL) {
	    gtk_widget_set_sensitive(w, TRUE);
	}
	vwin->flags |= VWIN_CONTENT_CHANGED;
    }
}

void mark_vwin_content_saved (windata_t *vwin) 
{
    GtkWidget *w = g_object_get_data(G_OBJECT(vwin->mbar), "save_button");

    if (w != NULL) {
	gtk_widget_set_sensitive(w, FALSE);
    }

    vwin->flags &= ~VWIN_CONTENT_CHANGED;

    w = g_object_get_data(G_OBJECT(vwin->mbar), "save_as_button");
    if (w != NULL) {
	gtk_widget_set_sensitive(w, TRUE);
    }
}

void buf_edit_save (GtkWidget *widget, windata_t *vwin)
{
    char **pbuf = (char **) vwin->data;
    gchar *text;

    text = textview_get_text(vwin->text);

    if (text == NULL || *text == '\0') {
	errbox(_("Buffer is empty"));
	g_free(text);
	return;
    }

    /* swap the edited text into the buffer */
    free(*pbuf); 
    *pbuf = text;

    if (vwin->role == EDIT_HEADER) {
	mark_vwin_content_saved(vwin);
	mark_dataset_as_modified();
    } else if (vwin->role == EDIT_NOTES) {
	mark_vwin_content_saved(vwin);
	mark_session_changed();
    }
}

static int update_func_code (windata_t *vwin)
{
    int iface, err = 0;

    /* callback used when editing a function in the context of
       the "function package editor" */
	
    iface = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(vwin->text), "iface"));
    err = update_function_from_script(vwin->fname, iface);
    if (err) {
	gui_errmsg(err);
    }

    return err;
}

void view_window_save (GtkWidget *widget, windata_t *vwin)
{
    if (vwin->role == EDIT_FUNC_CODE && *vwin->fname == '\0') {
	/* function package, sample script window */
	update_sample_script(vwin);
	mark_vwin_content_saved(vwin);
    } else if (strstr(vwin->fname, "script_tmp") || *vwin->fname == '\0') {
	/* special case: a newly created script */
	if (vwin->role == EDIT_SCRIPT) {
	    file_selector(_("Save command script"), SAVE_SCRIPT, 
			  FSEL_DATA_VWIN, vwin);
	} else if (vwin->role == EDIT_GP) {
	    file_selector(_("Save gnuplot commands"), SAVE_GP_CMDS, 
			  FSEL_DATA_VWIN, vwin);
	} else if (vwin->role == EDIT_R) {
	    file_selector(_("Save R commands"), SAVE_R_CMDS, 
			  FSEL_DATA_VWIN, vwin);
	}
    } else {
	FILE *fp;
	gchar *text;

	if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	    errbox(_("Can't open file for writing"));
	    return;
	} else {
	    text = textview_get_text(vwin->text);
	    system_print_buf(text, fp);
	    fclose(fp);
	    g_free(text);
	    mark_vwin_content_saved(vwin);
	    if (vwin->role == EDIT_FUNC_CODE) {
		update_func_code(vwin);
	    }
	}
    }
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
	if (vwin->text != NULL) { 
	    gchar *undo = g_object_steal_data(G_OBJECT(vwin->text), "undo");
	    
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
	if (vwin->ui != NULL) {
	    g_object_unref(G_OBJECT(vwin->ui));
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

	winstack_remove(vwin->main);
	free(vwin);
    }
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
	    g_signal_connect(G_OBJECT(vwin->popup), "destroy",
			     G_CALLBACK(gtk_widget_destroyed), 
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
    case VIEW_CODEBOOK:
	if (strstr(fname, "script_tmp") || strstr(fname, "session.inp")) {
	    title = g_strdup(_("gretl: command script"));
	} else {
	    title = title_from_filename(fname);
	} 
	break;
    case EDIT_NOTES:
	title = g_strdup(_("gretl: session notes")); break;
    case EDIT_GP:
	title = g_strdup(_("gretl: edit plot commands")); break;
    case EDIT_R:
	title = g_strdup(_("gretl: edit R commands")); break;
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
    mark_vwin_content_changed(vwin);
}

static void attach_content_changed_signal (windata_t *vwin)
{
    GtkTextBuffer *tbuf;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    g_signal_connect(G_OBJECT(tbuf), "changed", 
		     G_CALLBACK(content_changed), vwin);
}

static void viewer_box_config (windata_t *vwin)
{
    vwin->vbox = gtk_vbox_new(FALSE, 1);
    gtk_box_set_spacing(GTK_BOX(vwin->vbox), 4);
    gtk_container_set_border_width(GTK_CONTAINER(vwin->vbox), 4);

#ifndef G_OS_WIN32
    g_signal_connect_after(G_OBJECT(vwin->main), "realize", 
			   G_CALLBACK(set_wm_icon), 
			   NULL);
#endif

    gtk_container_add(GTK_CONTAINER(vwin->main), vwin->vbox);
}

static void view_buffer_insert_text (windata_t *vwin, PRN *prn)
{
    if (prn != NULL) {
	const char *buf = gretl_print_get_trimmed_buffer(prn);

	if (vwin->role == VIEW_FUNC_CODE || vwin->role == EDIT_FUNC_CODE) {
	    sourceview_insert_buffer(vwin, buf);
	} else if (vwin->role == SCRIPT_OUT) {
	    textview_set_text_colorized(vwin->text, buf);
	} else {
	    textview_set_text(vwin->text, buf);
	}
    }
}

static windata_t *reuse_script_out (windata_t *vwin, PRN *prn)
{
    int sticky = (vwin->flags & VWIN_STICKY);
    GtkTextBuffer *buf;
    const char *newtext;

    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

    if (sticky) {
	/* append to previous content */
	GtkTextMark *mark;
	GtkTextIter iter;

	gtk_text_buffer_get_end_iter(buf, &iter);
	mark = gtk_text_buffer_create_mark(buf, NULL, &iter, TRUE);
	textview_append_text_colorized(vwin->text, newtext, 1);
	gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(vwin->text), 
				     mark, 0.0, TRUE, 0, 0.05);
	gtk_text_buffer_delete_mark(buf, mark);
    } else {
	/* replace previous content */
	gtk_text_buffer_set_text(buf, "", -1);
	textview_set_text_colorized(vwin->text, newtext);
	cursor_to_top(vwin);
    }

    gretl_print_destroy(prn);

    gtk_window_present(GTK_WINDOW(vwin->main));

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
    int w, nlines;

    if (role == SCRIPT_OUT && script_out != NULL) {
	return reuse_script_out(script_out, prn);
    }

    if (title != NULL) {
	vwin = gretl_viewer_new(role, title, data, record);
    } else {
	gchar *tmp = make_viewer_title(role, NULL);

	vwin = gretl_viewer_new(role, tmp, data, record);
	g_free(tmp);
    }

    if (vwin == NULL) return NULL;

    viewer_box_config(vwin);

    if (role == VAR || role == VECM || role == SYSTEM) {
	/* special case: use a text-based menu bar */
	set_up_viewer_menu(vwin->main, vwin, system_items, sys_ui);
	add_system_menu_items(vwin, role);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
	gtk_widget_show(vwin->mbar);
	gretl_object_ref(data, (role == SYSTEM)? GRETL_OBJ_SYS : GRETL_OBJ_VAR);
    } else if (role == VIEW_FUNC_CODE || 
	       role == EDIT_FUNC_CODE ||
	       role == VIEW_MODELTABLE) {
	vwin_add_viewbar(vwin, 0);
    } else if (role != IMPORT) {
	vwin_add_viewbar(vwin, 1);
    }

    gretl_print_get_size(prn, &w, &nlines);
#if 1
    if (role != SCRIPT_OUT && w > 0 && w + 2 < hsize) {
	hsize = w + 2;
    }
#endif

    if (role == VIEW_FUNC_CODE) {
	create_source(vwin, hsize, vsize, FALSE);
    } else if (role == EDIT_FUNC_CODE) {
	create_source(vwin, hsize, vsize, TRUE);
    } else {
	create_text(vwin, hsize, vsize, nlines, FALSE);
	if (role == PRINT || role == SCRIPT_OUT ||
	    role == VIEW_MODELTABLE) {
	    text_set_word_wrap(vwin->text, 0);
	}
    }

    text_table_setup(vwin->vbox, vwin->text);

    if (role == SCRIPT_OUT && data != NULL) {
	/* partial output window for script */
	vwin_add_child((windata_t *) data, vwin);
    }

    /* register destruction of script output viewer */
    if (role == SCRIPT_OUT) {
	g_signal_connect(G_OBJECT(vwin->main), "destroy", 
			 G_CALLBACK(nullify_script_out), &script_out);
	script_out = vwin;
    }

    /* insert and then free the text buffer */
    view_buffer_insert_text(vwin, prn);
    gretl_print_destroy(prn);

    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    if (role == EDIT_FUNC_CODE) {
	g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
	attach_content_changed_signal(vwin);
	g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
			 G_CALLBACK(query_save_text), vwin);
    } 

    g_signal_connect(G_OBJECT(vwin->text), "button-press-event", 
		     G_CALLBACK(text_popup_handler), vwin);
    cursor_to_top(vwin);

    return vwin;
}

#define view_file_use_sourceview(r) (r == EDIT_SCRIPT || \
                                     r == VIEW_SCRIPT || \
                                     r == VIEW_LOG || \
                                     r == EDIT_GP || \
				     r == EDIT_R)

#define record_on_winstack(r) (r != EDIT_SCRIPT && \
	                       r != EDIT_GP && \
                               r != EDIT_R && \
                               r != VIEW_LOG && \
                               r != VIEW_SCRIPT && \
                               r != CONSOLE)

#define editing_script(r) (r == EDIT_SCRIPT || \
	                   r == EDIT_GP || \
                           r == EDIT_R)

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
    vwin = gretl_viewer_new(role, (title != NULL)? title : filename, 
			    NULL, record_on_winstack(role));
    g_free(title);

    if (vwin == NULL) {
	return NULL;
    }

    strcpy(vwin->fname, filename);

    viewer_box_config(vwin);
    vwin_add_viewbar(vwin, (role == VIEW_DATA || role == CONSOLE || role == VIEW_FILE));

    if (view_file_use_sourceview(role)) {
	create_source(vwin, hsize, vsize, editable);
    } else {
	create_text(vwin, hsize, vsize, 0, editable);
    }

    text_table_setup(vwin->vbox, vwin->text);

    if (view_file_use_sourceview(role)) {
	sourceview_insert_file(vwin, filename);
    } else {
	textview_insert_file(vwin, filename);
    }

    /* catch some special keystrokes */
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (editable) {
	g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    }

    /* editing script or graph: grab the "changed" signal and
       set up alert for unsaved changes on exit */
    if (editing_script(role)) {
	attach_content_changed_signal(vwin);	
	g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
			 G_CALLBACK(query_save_text), vwin);
    }

    /* clean up when dialog is destroyed */
    if (del_file) {
	gchar *fname = g_strdup(filename);

	g_signal_connect(G_OBJECT(vwin->main), "destroy", 
			 G_CALLBACK(delete_file), (gpointer) fname);
    }

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    g_signal_connect(G_OBJECT(vwin->text), "button-press-event", 
		     G_CALLBACK(text_popup_handler), vwin);

    cursor_to_top(vwin);
    gtk_widget_grab_focus(vwin->text);

    return vwin;
}

windata_t *
view_help_file (const char *filename, int role, GtkActionEntry *menu_items,
		const gchar *ui_info)
{
    windata_t *vwin;
    gchar *fbuf = NULL;
    gchar *title = NULL;
    int hsize = 80, vsize = 400;

    /* grab content of the appropriate help file into a buffer */
    gretl_file_get_contents(filename, &fbuf);
    if (fbuf == NULL) {
	return NULL;
    }

    title = make_viewer_title(role, NULL);
    vwin = gretl_viewer_new(role, title, NULL, 0);
    g_free(title);

    if (vwin == NULL) return NULL;

    strcpy(vwin->fname, filename);
    vwin->data = fbuf;

    viewer_box_config(vwin);
    set_up_viewer_menu(vwin->main, vwin, menu_items, ui_info);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    if (role == FUNCS_HELP) {
	hsize = 82;
	vsize = 500;
    }

    create_text(vwin, hsize, vsize, 0, FALSE);
    text_table_setup(vwin->vbox, vwin->text);

    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    if (vwin->role == CLI_HELP || vwin->role == CLI_HELP_EN ||
	vwin->role == FUNCS_HELP) {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event",
			 G_CALLBACK(help_popup_handler), 
			 vwin);
    } else {
	g_signal_connect(G_OBJECT(vwin->text), "button-press-event", 
			 G_CALLBACK(text_popup_handler), vwin);
    }	

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

    /* make the helpfile variant discernible via vwin->text */
    g_object_set_data(G_OBJECT(vwin->text), "role", GINT_TO_POINTER(vwin->role));

    gtk_widget_grab_focus(vwin->text);

    return vwin;
}

void view_window_set_editable (windata_t *vwin)
{
    gtk_text_view_set_editable(GTK_TEXT_VIEW(vwin->text), TRUE);
    gtk_text_view_set_cursor_visible(GTK_TEXT_VIEW(vwin->text), TRUE);
    g_object_set_data(G_OBJECT(vwin->main), "vwin", vwin);
    g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
		     G_CALLBACK(query_save_text), vwin);
    vwin->role = EDIT_SCRIPT;
    viewbar_add_edit_items(vwin);
    attach_content_changed_signal(vwin);
}

gint query_save_text (GtkWidget *w, GdkEvent *event, 
		      windata_t *vwin)
{
    if (vwin->flags & VWIN_CONTENT_CHANGED) {
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

    vwin = gretl_viewer_new(role, title, pbuf, 1);
    if (vwin == NULL) {
	return NULL;
    }

    viewer_box_config(vwin); 

    /* add a tool bar */
    vwin_add_viewbar(vwin, 0);

    create_text(vwin, hsize, vsize, 0, TRUE);
    text_table_setup(vwin->vbox, vwin->text);
    
    /* insert the buffer text */
    if (*pbuf) {
	GtkTextBuffer *tbuf = 
	    gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));

	gtk_text_buffer_set_text(tbuf, *pbuf, -1);
    }
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event", 
		     G_CALLBACK(text_popup_handler), vwin);
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);

    attach_content_changed_signal(vwin);

    /* alert for unsaved changes on exit */
    g_signal_connect(G_OBJECT(vwin->main), "delete-event",
		     G_CALLBACK(query_save_text), vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show(vwin->main);

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
    int w, nlines;

    vwin = gretl_viewer_new(VIEW_MODEL, title, pmod, 1);
    if (vwin == NULL) {
	return 1;
    }

    /* Take responsibility for one reference to this model */
    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    viewer_box_config(vwin);

    set_up_model_view_menu(vwin->main, vwin);

    g_signal_connect(G_OBJECT(vwin->mbar), "button-press-event", 
		     G_CALLBACK(check_model_menu), vwin);

    gtk_box_pack_start(GTK_BOX(vwin->vbox), vwin->mbar, FALSE, TRUE, 0);
    gtk_widget_show(vwin->mbar);

    gretl_print_get_size(prn, &w, &nlines);
    if (w + 2 < hsize) {
	hsize = w + 2;
    }

    create_text(vwin, hsize, vsize, nlines, FALSE);
    text_table_setup(vwin->vbox, vwin->text);

    /* insert and then free the model results buffer */
    buf = gretl_print_get_trimmed_buffer(prn);
    textview_set_text(vwin->text, buf);
    gretl_print_destroy(prn);

    /* attach shortcuts */
    g_signal_connect(G_OBJECT(vwin->main), "key-press-event", 
		     G_CALLBACK(catch_viewer_key), vwin);
    g_signal_connect(G_OBJECT(vwin->text), "button-press-event", 
		     G_CALLBACK(text_popup_handler), vwin);

    /* don't allow deletion of model window when a model
       test dialog is active */
    g_signal_connect(G_OBJECT(vwin->main), "delete-event", 
		     G_CALLBACK(check_delete_model_window), 
		     vwin);

    gtk_widget_show(vwin->vbox);
    gtk_widget_show_all(vwin->main);

    cursor_to_top(vwin);

    return 0;
}

void auto_save_plot (windata_t *vwin)
{
    gchar *buf;

    buf = textview_get_text(vwin->text);
    if (buf == NULL) {
	return;
    }

    dump_plot_buffer(buf, vwin->fname, 1);
    g_free(buf);

    mark_vwin_content_saved(vwin);
}

static void auto_save_script (windata_t *vwin)
{
    FILE *fp;
    gchar *savestuff;
    int unsaved = 0;

    if (strstr(vwin->fname, "script_tmp") || *vwin->fname == '\0') {
	file_selector(_("Save command script"), SAVE_SCRIPT, 
		      FSEL_DATA_VWIN, vwin);
	strcpy(vwin->fname, scriptfile);
	unsaved = 1;
    }

    if ((fp = gretl_fopen(vwin->fname, "w")) == NULL) {
	file_write_errbox(vwin->fname);
	return;
    }

    savestuff = textview_get_text(vwin->text);
    fprintf(fp, "%s", savestuff);
    g_free(savestuff); 
    fclose(fp);

    mark_vwin_content_saved(vwin);
}

#define dw_pval_ok(m) ((m->ci == OLS || m->ci == PANEL) && !na(pmod->dw))

static void get_ci_and_opt (const gchar *s, int *ci, gretlopt *opt)
{
    char c, word[9];

    sscanf(s, "%8[^:]:%c", word, &c);
    *ci = gretl_command_number(word);
    *opt = opt_from_flag((unsigned char) c);
}

static void set_tests_menu_state (GtkUIManager *ui, const MODEL *pmod)
{
    gretlopt opt;
    char path[128];
    const gchar *s;
    int i, n, ci;

    if (pmod->ci == MPOLS) {
	/* can we relax this? */
	flip(ui, "/MenuBar/Tests", FALSE);
	return;
    }

    n = G_N_ELEMENTS(model_test_items);

    for (i=0; i<n; i++) {
	opt = OPT_NONE;
	s = model_test_items[i].name;
	if (strchr(s, ':')) {
	    get_ci_and_opt(s, &ci, &opt);
	} else if (!strcmp(s, "dwpval")) {
	    sprintf(path, "/MenuBar/Tests/%s", s);
	    flip(ui, path, dw_pval_ok(pmod));
	    continue;
	} else {
	    ci = gretl_command_number(s);
	}
	sprintf(path, "/MenuBar/Tests/%s", s);
	flip(ui, path, model_test_ok(ci, opt, pmod, datainfo));
    }
}

static void model_save_state (GtkUIManager *ui, gboolean s)
{
    flip(ui, "/MenuBar/File/SaveAsIcon", s);
    flip(ui, "/MenuBar/File/SaveAndClose", s);
}

static void arma_x12_menu_mod (windata_t *vwin)
{
    flip(vwin->ui, "/MenuBar/Analysis/Covariance", FALSE);
    add_x12_output_menu_item(vwin);
}

static void rq_coeff_intervals_mod (windata_t *vwin)
{
    flip(vwin->ui, "/MenuBar/Analysis/ConfIntervals", FALSE);
}

#define intervals_model(m) (m->ci == LAD && \
			    gretl_model_get_data(m, "coeff_intervals"))

static void adjust_model_menu_state (windata_t *vwin, const MODEL *pmod)
{
    set_tests_menu_state(vwin->ui, pmod);

    /* disallow saving an already-saved model */
    if (pmod->name != NULL) {
	model_save_state(vwin->ui, FALSE);
    }

    if (RQ_SPECIAL_MODEL(pmod)) {
	/* can we relax this later? */
	flip(vwin->ui, "/MenuBar/Tests", FALSE);
	flip(vwin->ui, "/MenuBar/Save", FALSE);
	flip(vwin->ui, "/MenuBar/Analysis", FALSE);
	return;
    }

    if (intervals_model(pmod)) {
	rq_coeff_intervals_mod(vwin);
    }

    if (pmod->ci == MLE || pmod->ci == GMM) {
	/* can we relax some of this later? */
	flip(vwin->ui, "/MenuBar/Analysis/DisplayAFR", FALSE);
	flip(vwin->ui, "/MenuBar/Analysis/Forecasts", FALSE);
	flip(vwin->ui, "/MenuBar/Graphs", FALSE);
    } else if (pmod->ci == ARMA && arma_by_x12a(pmod)) {
	arma_x12_menu_mod(vwin);
    } 

    if (pmod->ci == GMM) {
	/* FIXME */
	flip(vwin->ui, "/MenuBar/Save", FALSE);
    }

    if (pmod->ncoeff == 1) {
	flip(vwin->ui, "/MenuBar/Analysis/ConfEllipse", FALSE);
    }

    if (pmod->ci == ARBOND) {
	flip(vwin->ui, "/MenuBar/Analysis/Forecasts", FALSE);
    }

    if (pmod->ci != OLS || !pmod->ifc || na(pmod->ess) || na(pmod->tss)) {
	flip(vwin->ui, "/MenuBar/Analysis/ANOVA", FALSE);
    }

    if (!bootstrap_ok(pmod->ci)) {
	flip(vwin->ui, "/MenuBar/Analysis/Bootstrap", FALSE);
    }
}

static void set_up_viewer_menu (GtkWidget *window, windata_t *vwin, 
				GtkActionEntry *items,
				const gchar *ui_info)
{
    gint n = 0;

    while (items[n].name != NULL) {
	n++;
    }

    vwin_add_ui(vwin, items, n, ui_info);

    if (vwin->role == VAR || vwin->role == VECM || 
	vwin->role == SYSTEM) {
	model_save_state(vwin->ui, !is_session_model(vwin->data));
    }
}

static GtkActionEntry model_data_base_items[] = {
    { "yhat", NULL, N_("_Fitted values"), NULL, NULL, 
      G_CALLBACK(fit_resid_callback) },
    { "uhat", NULL, N_("_Residuals"), NULL, NULL, 
      G_CALLBACK(fit_resid_callback) },
    { "uhat2", NULL, N_("_Squared residuals"), NULL, NULL, 
      G_CALLBACK(fit_resid_callback) }
};

static GtkActionEntry ess_items[] = {
    { "ess", NULL, N_("_Error sum of squares"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) },
    { "se", NULL, N_("_Standard error of the regression"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) }
}; 

static GtkActionEntry r_squared_items[] = {
    { "rsq", NULL, N_("_R-squared"), NULL, NULL, G_CALLBACK(model_stat_callback) },
    { "trsq", NULL, N_("_T*R-squared"), NULL, NULL, G_CALLBACK(model_stat_callback) }
}; 

static GtkActionEntry lnl_data_items[] = {
    { "lnL", NULL, N_("_Log likelihood"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry criteria_items[] = {
    { "AIC", NULL, N_("_Akaike Information Criterion"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) },
    { "BIC", NULL, N_("_Bayesian Information Criterion"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) },
    { "HQC", NULL, N_("_Hannan-Quinn Information Criterion"), NULL, NULL, 
      G_CALLBACK(model_stat_callback) }
};

static GtkActionEntry garch_data_items[] = {
    { "h", NULL, N_("_Predicted error variance"), NULL, NULL, 
      G_CALLBACK(fit_resid_callback)
    }
};

static GtkActionEntry fixed_effects_data_items[] = {
    { "ahat", NULL, N_("Per-unit _constants"), NULL, NULL, 
      G_CALLBACK(fit_resid_callback)
    }
};

static GtkActionEntry define_var_items[] = {
    /* Under Save; Sep wanted */
    { "NewVar", NULL, N_("Define _new variable..."), NULL, NULL,
      G_CALLBACK(model_genr_callback) }
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
    const gchar *path = "/MenuBar/Save";
    MODEL *pmod = vwin->data;

    vwin_menu_add_items(vwin, path, model_data_base_items,
			G_N_ELEMENTS(model_data_base_items));
			
    if (gretl_model_get_data(pmod, "ahat") != NULL) {
	vwin_menu_add_items(vwin, path, fixed_effects_data_items,
			    G_N_ELEMENTS(fixed_effects_data_items));
    }

    if (pmod->ci != GARCH) {
	vwin_menu_add_items(vwin, path, ess_items,
			    G_N_ELEMENTS(ess_items));
    }

    if (!ML_ESTIMATOR(pmod->ci) && pmod->ci != LAD && !na(pmod->rsq)) {
	vwin_menu_add_items(vwin, path, r_squared_items,
			    G_N_ELEMENTS(r_squared_items));
    }

    if (!na(pmod->lnL)) {
	vwin_menu_add_items(vwin, path, lnl_data_items,
			    G_N_ELEMENTS(lnl_data_items));
    }

    if (criteria_available(pmod)) {
	vwin_menu_add_items(vwin, path, criteria_items,
			    G_N_ELEMENTS(criteria_items));
    }

    if (pmod->ci == GARCH) {
	vwin_menu_add_items(vwin, path, garch_data_items,
			    G_N_ELEMENTS(garch_data_items));
    }

    vwin_menu_add_separator(vwin, path);

    vwin_menu_add_items(vwin, path, define_var_items,
			G_N_ELEMENTS(define_var_items));
}

static void add_model_tex_items (windata_t *vwin)
{
    MODEL *pmod = (MODEL *) vwin->data;
    int eqn_ok = command_ok_for_model(EQNPRINT, 0, pmod->ci);
    GtkActionGroup *actions;
    GError *err = NULL;
    int imod = 0;

    gtk_ui_manager_add_ui_from_string(vwin->ui, model_tex_ui, -1, &err);

    if (err != NULL) {
	g_message("building LaTeX menu failed: %s", err->message);
	g_error_free(err);
	return;
    }	

    actions = gtk_action_group_new("ModelTeX");
#ifdef ENABLE_NLS
    gtk_action_group_set_translation_domain(actions, "gretl");
#endif
    gtk_action_group_add_actions(actions, model_tex_items, 
				 G_N_ELEMENTS(model_tex_items),
				 vwin);
    gtk_action_group_add_radio_actions(actions, tex_eqn_items, 
				       G_N_ELEMENTS(tex_eqn_items),
				       (get_tex_eqn_opt() == OPT_T),
				       G_CALLBACK(set_tex_eqn_opt),
				       vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    if (intervals_model(pmod)) {
	eqn_ok = 0;
	imod = 1;
    }

    if (!eqn_ok || pmod->errcode) {
	flip(vwin->ui, "/MenuBar/LaTeX/TeXView/EqnView", FALSE);
	flip(vwin->ui, "/MenuBar/LaTeX/TeXCopy/EqnCopy", FALSE);
	flip(vwin->ui, "/MenuBar/LaTeX/TeXSave/EqnSave", FALSE);
	flip(vwin->ui, "/MenuBar/LaTeX/EqnOpts", FALSE);
    }

    if (imod) {
	flip(vwin->ui, "/MenuBar/LaTeX/TabOpts", FALSE);
    }
}

/* dummy placeholder, for when TeX is not supported */

static void add_missing_tex_items (windata_t *vwin)
{
    GtkActionGroup *actions;
    GError *err = NULL;

    gtk_ui_manager_add_ui_from_string(vwin->ui, missing_tex_ui, -1, &err);
    if (err != NULL) {
	g_message("building menus failed: %s", err->message);
	g_error_free(err);
	return;
    }	

    actions = gtk_action_group_new("MissingTeX");
    gtk_action_group_add_actions(actions, missing_tex_items, 
				 G_N_ELEMENTS(missing_tex_items),
				 vwin);
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    flip(vwin->ui, "/MenuBar/LaTeX", FALSE);
}

#define VNAMELEN2 32

static void add_vars_to_plot_menu (windata_t *vwin)
{
    GtkActionEntry entry;
    const gchar *mpath[] = {
	"/MenuBar/Graphs/ResidPlot", 
	"/MenuBar/Graphs/FittedActualPlot"
    };
    MODEL *pmod = vwin->data;
    char tmp[VNAMELEN2], aname[16];
    gchar *alabel;
    int *xlist;
    int v1, v2;
    int i, j;

    action_entry_init(&entry);

    xlist = gretl_model_get_x_list(pmod);
    
    for (i=0; i<2; i++) {
	/* plot against time/obs number */
	if (dataset_is_time_series(datainfo)) {
	    entry.label = _("_Against time");
	} else {
	    entry.label = _("By _observation number");
	}
	entry.name = (i == 0)? "r:byobs" : "f:byobs";
	entry.callback = (i == 0)? G_CALLBACK(resid_plot) : 
	    G_CALLBACK(fit_actual_plot);
	vwin_menu_add_item(vwin, mpath[i], &entry);

	if (pmod->ci == NLS || 
	    pmod->ci == MLE || 
	    pmod->ci == GMM ||
	    pmod->ci == PANEL) {
	    continue;
	}

	/* if doing resid plot, put dependent var in menu */
	if (i == 0) {
	    v1 = gretl_model_get_depvar(pmod);
	    if (v1 > 0) {
		sprintf(aname, "r:xvar %d", v1); /* FIXME */
		double_underscores(tmp, datainfo->varname[v1]);
		alabel = g_strdup_printf(_("_Against %s"), tmp);
		entry.name = aname;
		entry.label = alabel;
		entry.callback = G_CALLBACK(resid_plot);
		vwin_menu_add_item(vwin, mpath[0], &entry);
		g_free(alabel);
	    }
	}

	if (xlist == NULL) {
	    continue;
	}

	/* put the independent vars on the menu list */
	for (j=1; j<=xlist[0]; j++) {
	    v1 = xlist[j];
	    if (v1 == 0) {
		continue;
	    }
	    if (!strcmp(datainfo->varname[v1], "time")) {
		continue;
	    }
	    sprintf(aname, "%c:xvar %d", (i == 0)? 'r' : 'f', v1);
	    double_underscores(tmp, datainfo->varname[v1]);
	    alabel = g_strdup_printf(_("_Against %s"), tmp);
	    entry.name = aname;
	    entry.label = alabel;
	    entry.callback = (i == 0)? G_CALLBACK(resid_plot) : 
		G_CALLBACK(fit_actual_plot);
	    vwin_menu_add_item(vwin, mpath[i], &entry);
	    g_free(alabel);
	}
    }

    /* time series models: residual correlogram, spectrum */
    if (dataset_is_time_series(datainfo)) {
	vwin_menu_add_separator(vwin, "/MenuBar/Graphs");
	entry.name = "Correlogram";
	entry.label = _("Residual _correlogram");
	entry.callback = G_CALLBACK(residual_correlogram);
	vwin_menu_add_item(vwin, "/MenuBar/Graphs", &entry);
	entry.name = "Spectrum";
	entry.label = _("Residual _spectrum");
	entry.callback = G_CALLBACK(residual_periodogram);
	vwin_menu_add_item(vwin, "/MenuBar/Graphs", &entry);
    }

    /* 3-D fitted versus actual plot? */
    if (xlist != NULL) {
	v1 = v2 = 0;
	if (pmod->ifc && xlist[0] == 3) {
	    v1 = xlist[2];
	    v2 = xlist[3];
	} else if (!pmod->ifc && xlist[0] == 2) {
	    v1 = xlist[1];
	    v2 = xlist[2];
	}
	if (v1 > 0 && v2 > 0) {
	    char tmp2[VNAMELEN2];

	    vwin_menu_add_separator(vwin, mpath[1]);
	    double_underscores(tmp, datainfo->varname[v1]);
	    double_underscores(tmp2, datainfo->varname[v2]);
	    alabel = g_strdup_printf(_("_Against %s and %s"), tmp, tmp2);	
	    entry.name = "splot";
	    entry.label = alabel;
	    entry.callback = G_CALLBACK(fit_actual_splot);
	    vwin_menu_add_item(vwin, mpath[1], &entry);
	    g_free(alabel);
	}
    }

    free(xlist);
}

static void plot_dummy_call (GtkRadioAction *action, 
			     GtkRadioAction *current,
			     windata_t *vwin)
{
    vwin->active_var = gtk_radio_action_get_current_value(action);
}

static void radio_action_init (GtkRadioActionEntry *a)
{
    a->stock_id = NULL;
    a->accelerator = NULL;
    a->tooltip = NULL;
}

static void add_dummies_to_plot_menu (windata_t *vwin)
{
    GtkActionEntry item;
    GtkRadioActionEntry *items;
    MODEL *pmod = vwin->data;
    const gchar *gpath = "/MenuBar/Graphs/ResidPlot";
    const gchar *spath = "/MenuBar/Graphs/ResidPlot/Separation";
    char tmp[VNAMELEN2];
    int *dlist = NULL;
    int i, vi, ndums;

    /* make a list of dummy independent variables */
    for (i=2; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	if (vi == LISTSEP) {
	    break;
	} else if (vi > 0 &&
	    gretl_isdummy(datainfo->t1, datainfo->t2, Z[vi])) {
	    gretl_list_append_term(&dlist, vi);
	}
    }

    if (dlist == NULL) {
	return;
    }

    ndums = dlist[0];
    items = malloc((ndums + 1) * sizeof *items);
    if (items == NULL) {
	free(dlist);
	return;
    }

    /* add separator */
    vwin_menu_add_separator(vwin, gpath);

    /* add menu branch */
    action_entry_init(&item);
    item.name = "Separation";
    item.label = _("Separation");
    vwin_menu_add_menu(vwin, gpath, &item);

    /* configure "none" radio option */
    radio_action_init(&items[0]);
    items[0].name = "none";
    items[0].label = _("none");
    items[0].value = 0;

    /* put the dummy independent vars on the menu list */
    for (i=1; i<=dlist[0]; i++) {
	vi = dlist[i];
	radio_action_init(&items[i]);
	double_underscores(tmp, datainfo->varname[vi]);
	items[i].name = g_strdup_printf("dum %d", vi);
	items[i].label = g_strdup_printf(_("By %s"), tmp);
	items[i].value = vi;
    }

    vwin_menu_add_radios(vwin, spath, items, ndums + 1, 0,
			 G_CALLBACK(plot_dummy_call));

    for (i=1; i<=dlist[0]; i++) {
	g_free((gchar *) items[i].name);
	g_free((gchar *) items[i].label);
    }

    free(items);
    free(dlist);
}

static void varnum_from_action (GtkAction *action, int *i)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "%*s %d", i);
}

static void tau_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    MODEL *pmod = (MODEL *) vwin->data;
    int v, err;

    varnum_from_action(action, &v);
    err = plot_tau_sequence(pmod, datainfo, v);

    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }    
}

static void add_tau_plot_menu (windata_t *vwin)
{
    GtkActionEntry item;
    MODEL *pmod = vwin->data;
    char tmp[VNAMELEN2], aname[16];
    int i;

    action_entry_init(&item);
    item.name = "TauMenu";
    item.label = _("tau sequence");
    vwin_menu_add_menu(vwin, "/MenuBar/Graphs", &item);
    
    item.callback = G_CALLBACK(tau_plot_call);

    /* put the independent vars on the menu list */
    for (i=2; i<=pmod->list[0]; i++) {
	sprintf(aname, "tauseq %d", i - 2);
	double_underscores(tmp, datainfo->varname[pmod->list[i]]);
	item.name = aname;
	item.label = tmp;
	vwin_menu_add_item(vwin, "/MenuBar/Graphs/TauMenu", &item);
    }
}

static void x12_output_callback (GtkAction *action, gpointer p)
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

static const gchar *model_ui =
    "<ui>"
    " <menubar name='MenuBar'>"
    "  <menu action='File'>"
    "   <menuitem action='SaveAs'/>"
    "   <menuitem action='SaveAsIcon'/>"
    "   <menuitem action='SaveAndClose'/>"
#ifdef NATIVE_PRINTING
    "   <menuitem action='Print'/>"
#endif
    "   <menuitem action='Close'/>"
    "  </menu>"    
    "  <menu action='Edit'>"
    "   <menuitem action='Copy'/>"
    "   <menuitem action='Revise'/>"
    "  </menu> "     
    "  <menu action='Tests'>"
    "   <menuitem action='omit'/>"
    "   <menuitem action='add'/>"
    "   <menuitem action='coeffsum'/>"
    "   <menuitem action='restrict'/>"
    "   <separator/>"
    "   <menuitem action='lmtest:s'/>"
    "   <menuitem action='lmtest:l'/>"
    "   <menuitem action='reset'/>"
    "   <separator/>"
    "   <menu action='Hsk'/>"
    "   <menuitem action='normtest'/>"
    "   <menuitem action='leverage'/>"
    "   <menuitem action='vif'/>"
    "   <menuitem action='chow'/>"
    "   <separator/>"
    "   <menuitem action='lmtest:a'/>"
    "   <menuitem action='dwpval'/>"
    "   <menuitem action='lmtest:h'/>"
    "   <menuitem action='qlrtest'/>"
    "   <menuitem action='cusum'/>"
    "   <menuitem action='cusum:r'/>"
    "   <separator/>"
    "   <menuitem action='hausman'/>"
    "  </menu>"      
    "  <menu action='Save'/>"
    "  <menu action='Graphs'>"
    "   <menu action='ResidPlot'/>"
    "   <menu action='FittedActualPlot'/>"
    "  </menu>"    
    "  <menu action='Analysis'>"
    "   <menuitem action='DisplayAFR'/>"
    "   <menuitem action='Forecasts'/>"
    "   <menuitem action='ConfIntervals'/>"
    "   <menuitem action='ConfEllipse'/>"
    "   <menuitem action='Covariance'/>"
    "   <menuitem action='ANOVA'/>"
    "   <menuitem action='Bootstrap'/>"
    "  </menu>"     
    "  <menu action='LaTeX'/>"
    " </menubar>"
    "</ui>";

static void 
set_up_model_view_menu (GtkWidget *window, windata_t *vwin) 
{
    MODEL *pmod = (MODEL *) vwin->data;
    GtkActionGroup *actions;
    GError *err = NULL;

    actions = gtk_action_group_new("ModelActions");
#ifdef ENABLE_NLS
    gtk_action_group_set_translation_domain(actions, "gretl");
#endif

    gtk_action_group_add_actions(actions, model_items, 
				 G_N_ELEMENTS(model_items), 
				 vwin);
    gtk_action_group_add_actions(actions, model_test_items, 
				 G_N_ELEMENTS(model_test_items),
				 vwin);

    vwin->ui = gtk_ui_manager_new();
    gtk_ui_manager_insert_action_group(vwin->ui, actions, 0);
    g_object_unref(actions);

    gtk_ui_manager_add_ui_from_string(vwin->ui, model_ui, -1, &err);
    if (err != NULL) {
	g_message("building menus failed: %s", err->message);
	g_error_free(err);
    }

    if (pmod->ci != MLE && pmod->ci != GMM) {
	if (RQ_SPECIAL_MODEL(pmod)) {
	    add_tau_plot_menu(vwin);
	} else {
	    add_vars_to_plot_menu(vwin);
	}
	add_model_dataset_items(vwin);
    }

    if (dataset_is_panel(datainfo) && pmod->ci == OLS) {
	vwin_menu_add_items(vwin, "/MenuBar/Tests/Hsk", 
			    panel_hsk_items, 
			    G_N_ELEMENTS(panel_hsk_items));
    } else if (model_test_ok(LMTEST, OPT_W, pmod, datainfo)) {
	vwin_menu_add_items(vwin, "/MenuBar/Tests/Hsk", 
			    base_hsk_items, 
			    G_N_ELEMENTS(base_hsk_items));
    }

    if (latex_is_ok() && !pmod->errcode && !RQ_SPECIAL_MODEL(pmod)) {
	add_model_tex_items(vwin);
    } else {
	add_missing_tex_items(vwin);
    }

    if (pmod->ci != ARMA && pmod->ci != GARCH && 
	pmod->ci != NLS && pmod->ci != MLE && pmod->ci != GMM &&
	pmod->ci != PANEL && pmod->ci != ARBOND) {
	add_dummies_to_plot_menu(vwin);
    }

    if (vwin->main != NULL) {
	gtk_window_add_accel_group(GTK_WINDOW(vwin->main), 
				   gtk_ui_manager_get_accel_group(vwin->ui));
    }

    vwin->mbar = gtk_ui_manager_get_widget(vwin->ui, "/MenuBar");

    /* disable some menu items if need be */
    adjust_model_menu_state(vwin, pmod);
}

enum {
    SYS_DATA_RESIDS,
    SYS_DATA_FITTED,
    SYS_DATA_SIGMA
};

static int sys_data_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "uhat")) {
	return SYS_DATA_RESIDS;
    } else if (!strcmp(s, "yhat")) {
	return SYS_DATA_FITTED;	
    } else if (!strcmp(s, "sigma")) {
	return SYS_DATA_SIGMA;
    } else {
	return SYS_DATA_RESIDS;
    }
}

static void system_data_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *M = NULL;
    gchar *wtitle = NULL;
    PRN *prn;
    int code, k = 0;
    int err = 0;

    if (vwin->role == SYSTEM) {
	sys = (equation_system *) vwin->data;
    } else {
	var = (GRETL_VAR *) vwin->data;
    } 

    if ((var == NULL && sys == NULL) || bufopen(&prn)) {
	return;
    }

    code = sys_data_code(action);

    if (code == SYS_DATA_SIGMA) {
	if (var != NULL) {
	    wtitle = g_strdup(_("gretl: VAR covariance matrix"));
	    err = gretl_VAR_print_sigma(var, prn);
	} else {
	    wtitle = g_strdup(_("gretl: system covariance matrix"));
	    err = system_print_sigma(sys, prn);
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

static int VAR_model_data_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "VarIrf")) {
	return VAR_IRF;
    } else if (!strcmp(s, "VarDecomp")) {
	return VAR_DECOMP;
    } else {
	return VAR_IRF;
    }
}

static void VAR_model_data_callback (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = vwin->data;
    gchar *title;
    PRN *prn;
    int code, h = 0;
    int err;

    if (var == NULL || bufopen(&prn)) return;

    code = VAR_model_data_code(action);

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

static void add_x12_output_menu_item (windata_t *vwin)
{
    const gchar *mpath = "/MenuBar/Analysis";
    GtkActionEntry entry;

    vwin_menu_add_separator(vwin, mpath);

    action_entry_init(&entry);
    entry.name = "x12aout";
    entry.label = _("View X-12-ARIMA output");
    entry.callback = G_CALLBACK(x12_output_callback);
    vwin_menu_add_item(vwin, mpath, &entry);
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

static void impulse_params_from_action (GtkAction *action, 
					int *targ,
					int *shock)
{
    const gchar *s = gtk_action_get_name(action);

    sscanf(s, "Imp:%d:%d", targ, shock);
}

static void impulse_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int horizon, bootstrap;
    gint shock, targ;
    const double **vZ = NULL;
    int err;

    impulse_params_from_action(action, &targ, &shock);

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

static void multiple_irf_plot_call (GtkAction *action, gpointer p)
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

static void system_forecast_callback (GtkAction *action, gpointer p)
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
    int i, err = 0;

    varnum_from_action(action, &i);

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
			   t1, t2, &t2, NULL,
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

static int sys_test_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "autocorr")) {
	return SYS_AUTOCORR_TEST;
    } else if (!strcmp(s, "ARCH")) {
	return SYS_ARCH_TEST;
    } else if (!strcmp(s, "normtest")) {
	return SYS_NORMALITY_TEST;
    } else if (!strcmp(s, "restrict")) {
	return SYS_RESTRICT;
    } else {
	return SYS_NORMALITY_TEST;
    }
}

static void system_test_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    gchar *title = NULL;
    gchar *cstr = NULL;
    PRN *prn;
    int code, order = 0;
    int err = 0;

    if (bufopen(&prn)) {
	return;
    }

    code = sys_test_code(action);

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
	title = g_strdup(_("gretl: autocorrelation"));
	cstr = g_strdup_printf("lmtest %d --autocorr", order);
	if (var != NULL) {
	    err = gretl_VAR_autocorrelation_test(var, order, 
						 &Z, datainfo, 
						 prn);
	} else {
	    err = system_autocorrelation_test(sys, order, prn);
	}
    } else if (code == SYS_ARCH_TEST) {
	title = g_strdup(_("gretl: ARCH test"));
	cstr = g_strdup_printf("lmtest %d --arch", order);
	if (var != NULL) {
	    err = gretl_VAR_arch_test(var, order, datainfo, prn);
	} else {
	    err = system_arch_test(sys, order, prn);
	}
    } else if (code == SYS_NORMALITY_TEST) {
	title = g_strdup_printf("gretl: %s", _("Test for normality of residual"));
	cstr = g_strdup("testuhat");
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
	add_command_to_stack(cstr);
	view_buffer(prn, 78, 400, title, PRINT, NULL); 
    }

    g_free(title);
    g_free(cstr);
}

static void VAR_roots_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = (GRETL_VAR *) vwin->data;
    int err;

    err = gretl_VAR_roots_plot(var);
    
    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static int sys_ci_from_action (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);
    char cmdword[9];

    sscanf(s, "%*s %8s", cmdword);
    return gretl_command_number(cmdword);
}

static void system_resid_plot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int ci = sys_ci_from_action(action);
    int err;

    err = gretl_system_residual_plot(vwin->data, ci, datainfo);
    
    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static void system_resid_mplot_call (GtkAction *action, gpointer p)
{
    windata_t *vwin = (windata_t *) p;
    int ci = sys_ci_from_action(action);
    int err;

    err = gretl_system_residual_mplot(vwin->data, ci, datainfo);
    
    if (err) {
	gui_errmsg(err);
    } else {
	register_graph();
    }
}

static void add_system_menu_items (windata_t *vwin, int ci)
{
    GtkActionEntry item;
    const gchar *top = "/MenuBar";
    const gchar *tests = "/MenuBar/Tests";
    const gchar *save = "/MenuBar/Save";
    const gchar *graphs = "/MenuBar/Graphs";
    const gchar *analysis = "/MenuBar/Analysis";
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    int neqns, nfc, vtarg, vshock;
    char tmp[VNAMELEN2], istr[16];
    char maj[64], min[32];
    const char *cmdword;
    int i, j;

    if (ci == SYSTEM) {
	sys = (equation_system *) vwin->data;
	neqns = sys->neqns;
	nfc = sys->neqns + sys->nidents;
    } else {
	var = (GRETL_VAR *) vwin->data;
	nfc = neqns = gretl_VAR_get_n_equations(var);
    }

    cmdword = gretl_command_word(ci);   
    action_entry_init(&item);

    /* FIXME: the following two tests should really be multivariate */

    if (dataset_is_time_series(datainfo)) {
	/* univariate autocorrelation tests */
	item.name = "autocorr";
	item.label = N_("_Autocorrelation");
	item.callback = G_CALLBACK(system_test_call);
	vwin_menu_add_item(vwin, tests, &item);

	/* univariate ARCH tests */
	item.name = "ARCH";
	item.label = N_("A_RCH");
	vwin_menu_add_item(vwin, tests, &item);
    }

    /* multivariate normality test */
    item.name = "normtest";
    item.label = N_("_Normality of residuals");
    vwin_menu_add_item(vwin, tests, &item);

    if (ci == VECM || ci == SYSTEM) {
	/* linear restrictions (on cointegrating relations, for VECM) */
	item.name = "restrict";
	item.label = N_("Linear restrictions");
	item.callback = G_CALLBACK(gretl_callback);
	vwin_menu_add_item(vwin, tests, &item);
    } else if (ci == VAR) {
	/* regular VAR: omit exogenous variables test */
	if (gretl_VAR_get_exo_list(var) != NULL) {
	    item.name = "VarOmit";
	    item.label = N_("Omit exogenous variables...");
	    item.callback = G_CALLBACK(selector_callback);
	    vwin_menu_add_item(vwin, tests, &item);
	}	    
    }

    /* Save residuals */
    for (i=0; i<neqns; i++) {
	sprintf(istr, "resid %d", i);
	sprintf(maj, "%s %d", _("Residuals from equation"), i + 1);
	item.name = istr;
	item.label = maj;
	item.callback = G_CALLBACK(add_system_resid);
	vwin_menu_add_item(vwin, save, &item);
    }

    /* Display residual matrix */
    item.name = "uhat";
    item.label = N_("Display residuals, all equations");
    item.callback = G_CALLBACK(system_data_callback);
    vwin_menu_add_item(vwin, analysis, &item);

    if (ci == SYSTEM) {
	/* Display fitted values matrix */
	item.name = "yhat";
	item.label = N_("Display fitted values, all equations");
	vwin_menu_add_item(vwin, analysis, &item);
    }  

    /* Display VCV matrix */
    item.name = "sigma";
    item.label = N_("Cross-equation covariance matrix");
    vwin_menu_add_item(vwin, analysis, &item);

    if (ci == VAR || ci == VECM) {
	/* impulse response printout */
	item.name = "VarIrf";
	item.label = N_("Impulse responses");
	item.callback = G_CALLBACK(VAR_model_data_callback);
	vwin_menu_add_item(vwin, analysis, &item);

	/* variance decomp printout */
	item.name = "VarDecomp";
	item.label = N_("Forecast variance decomposition");
	vwin_menu_add_item(vwin, analysis, &item);
    }

    if (neqns <= 6) {
	/* multiple residual plots */
	sprintf(min, "multiresid %s", cmdword);
	item.name = min;
	item.label = N_("Residual plots");
	item.callback = G_CALLBACK(system_resid_mplot_call);
	vwin_menu_add_item(vwin, graphs, &item);
    }

    /* combined residual plot */
    sprintf(min, "comboresid %s", cmdword);
    item.name = min;
    item.label = N_("Combined residual plot");
    item.callback = G_CALLBACK(system_resid_plot_call);
    vwin_menu_add_item(vwin, graphs, &item);

    if (ci != SYSTEM) {
	/* VAR inverse roots */
	item.name = "VarRoots";
	item.label = N_("VAR inverse roots");
	item.callback = G_CALLBACK(VAR_roots_plot_call);
	vwin_menu_add_item(vwin, graphs, &item);
    }

    if (ci != SYSTEM && neqns <= 4) {
	/* Multiple IRFs */
	item.name = "MultiIrf";
	item.label = N_("Impulse responses (combined)");
	item.callback = G_CALLBACK(multiple_irf_plot_call);
	vwin_menu_add_item(vwin, graphs, &item);
    }

    for (i=0; i<nfc; i++) {
	char newpath[64];
	int dv;

	/* forecast items */
	if (var != NULL) {
	    dv = gretl_VAR_get_variable_number(var, i);
	} else {
	    dv = sys->ylist[i+1];
	}
	double_underscores(tmp, datainfo->varname[dv]);
	sprintf(istr, "fcast %d", i);
	item.name = istr;
	item.label = tmp;
	item.callback = G_CALLBACK(system_forecast_callback);
	vwin_menu_add_item(vwin, "/MenuBar/Analysis/Forecasts", &item);

	if (var == NULL) {
	    continue;
	}

	/* impulse response plots: make menu for target */
	vtarg = gretl_VAR_get_variable_number(var, i);
	double_underscores(tmp, datainfo->varname[vtarg]);
	sprintf(istr, "targ_%d", i);
	sprintf(maj, _("Response of %s"), tmp);
	item.name = istr;
	item.label = maj;
	item.callback = NULL;
	vwin_menu_add_menu(vwin, graphs, &item);

	/* path under which to add shocks */
	sprintf(newpath, "/MenuBar/Graphs/targ_%d", i);
	
	for (j=0; j<neqns; j++) {
	    /* impulse responses: subitems for shocks */
	    vshock = gretl_VAR_get_variable_number(var, j);
	    double_underscores(tmp, datainfo->varname[vshock]);
	    sprintf(istr, "Imp:%d:%d", i, j);
	    sprintf(min, _("to %s"), tmp);
	    item.name = istr;
	    item.label = min;
	    item.callback = G_CALLBACK(impulse_plot_call);
	    vwin_menu_add_item(vwin, newpath, &item);
	}
    }

    if (ci == VECM) {
	/* save ECs items */
	for (i=0; i<jrank(var); i++) {
	    sprintf(istr, "EC %d", i);
	    sprintf(maj, "%s %d", _("EC term"), i+1);
	    item.name = istr;
	    item.label = maj;
	    item.callback = G_CALLBACK(VECM_add_EC_data);
	    vwin_menu_add_item(vwin, save, &item);
	}
    }

    if (latex_is_ok()) {
	int n = G_N_ELEMENTS(sys_tex_items);

	vwin_menu_add_menu(vwin, top, &sys_tex_items[0]);
	vwin_menu_add_items(vwin, "/MenuBar/LaTeX",
			    sys_tex_items + 1, n - 1);
    }
}

static gint check_model_menu (GtkWidget *w, GdkEventButton *eb, 
			      gpointer data)
{
    windata_t *mwin = (windata_t *) data;
    MODEL *pmod = mwin->data;
    GtkAction *action;
    gboolean s;
    int ok = 1, graphs_ok = 1;

    if (RQ_SPECIAL_MODEL(pmod)) {
	return FALSE;
    }

    if (Z == NULL) {
	flip(mwin->ui, "/MenuBar/File/SaveAsIcon", FALSE);
	flip(mwin->ui, "/MenuBar/File/SaveAndClose", FALSE);
	flip(mwin->ui, "/MenuBar/Edit/Copy", FALSE);
	flip(mwin->ui, "/MenuBar/Tests", FALSE);
	flip(mwin->ui, "/MenuBar/Graphs", FALSE);
	flip(mwin->ui, "/MenuBar/Analysis", FALSE);
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

    action = gtk_ui_manager_get_action(mwin->ui, "/MenuBar/Tests");
    s = gtk_action_is_sensitive(action);
    if ((s && ok) || (!s && !ok)) {
	/* no need to flip state */
	return FALSE;
    }

    flip(mwin->ui, "/MenuBar/Tests", ok);
    flip(mwin->ui, "/MenuBar/Graphs", graphs_ok);
    flip(mwin->ui, "/MenuBar/Analysis/DisplayAFR", ok);
    flip(mwin->ui, "/MenuBar/Analysis/Forecasts", ok);
    flip(mwin->ui, "/MenuBar/Analysis/ConfIntervals", ok);
    flip(mwin->ui, "/MenuBar/Save/yhat", ok);
    flip(mwin->ui, "/MenuBar/Save/uhat", ok);
    flip(mwin->ui, "/MenuBar/Save/uhat2", ok);
    flip(mwin->ui, "/MenuBar/Save/NewVar", ok);

    if (!ok) {
	const char *msg = gretl_errmsg_get();

	if (msg != NULL && *msg != 0) {
	    infobox(msg);
	}
    } 

    return FALSE;
}

static gchar *exists_string (const char *name, GretlType t)
{
    gchar *s = NULL;

    if (t == GRETL_TYPE_SERIES) {
	s = g_strdup_printf(_("A series named %s already exists"), name);
    } else if (t == GRETL_TYPE_MATRIX) {
	s = g_strdup_printf(_("A matrix named %s already exists"), name);
    } else if (t == GRETL_TYPE_DOUBLE) {
	s = g_strdup_printf(_("A scalar named %s already exists"), name);
    } else if (t == GRETL_TYPE_LIST) {
	s = g_strdup_printf(_("A list named %s already exists"), name);
    } else if (t == GRETL_TYPE_STRING) {
	s = g_strdup_printf(_("A string named %s already exists"), name);
    }

    return s;
}

static int object_overwrite_ok (const char *name, GretlType t)
{
    gchar *info = exists_string(name, t);
    gchar *msg = g_strdup_printf("%s\n%s", info, _("OK to overwrite it?"));
    int resp;

    resp = yes_no_dialog("gretl", msg, 0);
    g_free(info);
    g_free(msg);

    return (resp == GRETL_YES);
}	    

int gui_validate_varname (const char *name, GretlType t)
{
    int i, n = strlen(name);
    char namebit[VNAMELEN];
    unsigned char c;
    int err = 0;

    *namebit = 0;
    
    if (n > VNAMELEN - 1) {
	strncat(namebit, name, VNAMELEN - 1);
	errbox(_("Variable name %s... is too long\n"
		 "(the max is %d characters)"), namebit,
	       VNAMELEN - 1);
	err = 1;
    } else if (!(isalpha(*name))) {
	errbox(_("First char of name ('%c') is bad\n"
		 "(first must be alphabetical)"), *name);
	err = 1;
    } else {
	for (i=1; i<n && !err; i++) {
	    c = (unsigned char) name[i];
	
	    if ((!(isalpha(c)) && !(isdigit(c)) && c != '_') || c > 127) {
		errbox(_("Name contains an illegal char (in place %d)\n"
			 "Use only unaccented letters, digits and underscore"), i + 1);
		err = 1;
	    }
	}
    }

    if (!err && t != GRETL_TYPE_NONE) {
	/* check for collisions */
	GretlType t0 = gretl_type_from_name(name, datainfo);

	if (t0 != GRETL_TYPE_NONE) {
	    if (t == t0) {
		err = !object_overwrite_ok(name, t);
	    } else {
		/* won't work */
		gchar *msg = exists_string(name, t0);
		
		errbox(msg);
		g_free(msg);
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

#ifdef G_OS_WIN32

static void run_R_sync (void)
{
    char *Rterm;
    gchar *cmd;
    int err;

    Rterm = R_path_from_registry();
    if (Rterm == NULL) {
	gui_errmsg(E_EXTERNAL);
	return;
    }

    cmd = g_strdup_printf("\"%s\" --no-save --no-init-file --no-restore-data "
			  "--slave", Rterm);

    err = winfork(cmd, NULL, SW_SHOWMINIMIZED, CREATE_NEW_CONSOLE);

    if (err) {
	gui_errmsg(err);
    } else {
	gchar *Rout = g_strdup_printf("%sR.out", paths.dotdir);

	view_file(Rout, 0, 1, 78, 350, VIEW_FILE);
	g_free(Rout);
    }

    g_free(cmd);
    free(Rterm);
}

#else /* some non-Windows functions follow */

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

/* Start an R session in asynchronous (interactive) mode.
   Note that there's a separate win32 function for this
   in gretlwin32.c.
*/

static void start_R_async (void)
{
    const char *supp1 = "--no-init-file";
    const char *supp2 = "--no-restore-data";
    char *s0 = NULL, *s1 = NULL, *s2 = NULL;
    pid_t pid;
    int n;

    s0 = mymalloc(64);
    s1 = mymalloc(32);
    s2 = mymalloc(32);

    if (s0 == NULL || s1 == NULL || s2 == NULL) {
	goto bailout;
    }

    *s0 = *s1 = *s2 = '\0';

    n = sscanf(Rcommand, "%63s %31s %31s", s0, s1, s2);

    if (n == 0) {
	errbox(_("No command was supplied to start R"));
	goto bailout;
    }

    signal(SIGCHLD, SIG_IGN); 
    pid = fork();

    if (pid == -1) {
	errbox(_("Couldn't fork"));
	perror("fork");
	return;
    } else if (pid == 0) {  
	if (n == 1) {
	    execlp(s0, s0, supp1, supp2, NULL);
	} else if (n == 2) {
	    execlp(s0, s0, s1, supp1, supp2, NULL);
	} else if (n == 3) {
	    execlp(s0, s0, s1, s2, supp1, supp2, NULL);
	}
	perror("execlp");
	_exit(EXIT_FAILURE);
    }

 bailout:

    free(s0); 
    free(s1); 
    free(s2);
}

/* run R in synchronous (batch) mode and display the results
   in a gretl window
*/

static void run_R_sync (void)
{
    gchar *argv[6];
    gchar *out = NULL;
    gchar *errout = NULL;
    gint status = 0;
    GError *gerr = NULL;
    PRN *prn = NULL;

    argv[0] = "R";
    argv[1] = "--no-save";
    argv[2] = "--no-init-file";
    argv[3] = "--no-restore-data";
    argv[4] = "--slave";
    argv[5] = NULL;

    signal(SIGCHLD, SIG_DFL);

    g_spawn_sync(NULL, argv, NULL, G_SPAWN_SEARCH_PATH,
		 NULL, NULL, &out, &errout,
		 &status, &gerr);

    if (gerr != NULL) {
	errbox(gerr->message);
	g_error_free(gerr);
    } else if (status != 0) {
	if (errout != NULL) {
	    if (*errout == '\0') {
		errbox("R exited with status %d", status);
	    } else if (strlen(errout) < MAXLEN) {
		errbox(errout);
	    } else {
		prn = gretl_print_new_with_buffer(errout);
		errout = NULL;
	    }
	}
    } else if (out != NULL) {
	prn = gretl_print_new_with_buffer(out);
	out = NULL;
    } else {
	warnbox("Got no output");
    }

    if (prn != NULL) {
	view_buffer(prn, 78, 350, _("R output"), PRINT, NULL);
    }

    g_free(out);
    g_free(errout);
}

#endif /* !G_OS_WIN32 */

/* driver for starting R, either interactive or in batch mode */

void start_R (const char *buf, int send_data, int interactive)
{
    gretlopt Ropt = OPT_G;
    int err;

    if (send_data && !data_status) {
	warnbox(_("Please open a data file first"));
	return;
    }

    if (interactive) {
	Ropt |= OPT_I;
    }

    if (send_data) {
	Ropt |= OPT_D;
    }

     err = write_gretl_R_files(buf, (const double **) Z,
			       datainfo, Ropt);

    if (err) {
	gui_errmsg(err);
	delete_gretl_R_files();
    } else if (interactive) {
#ifdef G_OS_WIN32
	win32_start_R_async();
#else
	start_R_async();
#endif
    } else {
	run_R_sync();
    }
}

void verbose_gerror_report (GError *gerr, const char *src)
{
    fprintf(stderr, "GError details from %s\n"
	    " message: '%s'\n domain = %d, code = %d\n",
	    src, gerr->message, gerr->domain, gerr->code);
}

int gretl_file_get_contents (const gchar *fname, gchar **contents)
{
    GError *gerr = NULL;
    gboolean ok;

    ok = g_file_get_contents(fname, contents, NULL, &gerr);

    if (gerr != NULL) {
	verbose_gerror_report(gerr, "g_file_get_contents");
	if (g_error_matches(gerr, G_FILE_ERROR, G_FILE_ERROR_INVAL)) {
	    gchar *trfname = NULL;
	    gsize bytes;

	    verbose_gerror_report(gerr, "g_file_get_contents");
	    g_error_free(gerr);
	    gerr = NULL;

	    if (!g_utf8_validate(fname, -1, NULL)) {
		fprintf(stderr, "Trying g_locale_to_utf8 on filename\n");
		trfname = g_locale_to_utf8(fname, -1, NULL, &bytes, &gerr);
		if (trfname == NULL) {
		    verbose_gerror_report(gerr, "g_locale_to_utf8");
		}
	    } else {
		fprintf(stderr, "Trying g_locale_from_utf8 on filename\n");
		trfname = g_locale_from_utf8(fname, -1, NULL, &bytes, &gerr);
		if (trfname == NULL) {
		    verbose_gerror_report(gerr, "g_locale_from_utf8");
		}
	    }

	    if (trfname != NULL) {
		ok = g_file_get_contents(trfname, contents, NULL, &gerr);
		g_free(trfname);
		if (!ok) {
		    verbose_gerror_report(gerr, "g_file_get_contents");
		}
	    }
	}
	if (gerr != NULL) {
	    errbox(gerr->message);
	    g_error_free(gerr);
	}
    }

    return !ok;
}	
