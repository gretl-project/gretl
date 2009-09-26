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

/* toolbar.c: main-window toolbar, viewer window toolbars, etc. */

#include "gretl.h"
#include "console.h"
#include "session.h"
#include "datafiles.h"
#include "selector.h"
#include "textbuf.h"
#include "textutil.h"
#include "series_view.h"
#include "model_table.h"
#include "cmdstack.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "toolbar.h"

#include "usermat.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

/* for viewer window toolbars */
#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mini.tsplot.xpm"
#include "../pixmaps/mini.boxplot.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.pin.xpm"
#include "../pixmaps/mini.alpha.xpm"
#include "../pixmaps/mini.en.xpm"
#include "../pixmaps/mini.split.xpm"
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 8)
# include "../pixmaps/info_24.xpm"
#endif
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 6)
# include "../pixmaps/edit_24.xpm"
# include "../pixmaps/mini.edit.xpm"
#endif

/* for main-window toolbar */
#include "../pixmaps/mini.calc.xpm"
#include "../pixmaps/mini.sh.xpm"
#include "../pixmaps/mini.session.xpm"
#include "../pixmaps/mini.plot.xpm"
#include "../pixmaps/mini.model.xpm"
#include "../pixmaps/mini.func.xpm"

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    PLOT_ITEM,
    EXEC_ITEM,
    COPY_ITEM,
    TEX_ITEM,
    ADD_DATA_ITEM,
    ADD_MATRIX_ITEM,
    MAIL_ITEM,
    HELP_ITEM,
    CMD_HELP_ITEM,
    GP_HELP_ITEM,
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    INDEX_ITEM,
    EDIT_SCRIPT_ITEM,
    STICKIFY_ITEM,
    ALPHA_ITEM,
    REFRESH_ITEM,
    OPEN_ITEM,
    SPLIT_ITEM
} viewbar_flags;

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
	mini_pin_xpm,
	mini_alpha_xpm,
	mini_en_xpm,
	mini_split_xpm
    };
    const char *stocks[] = {
#if NO_INFO_ICON
	GRETL_STOCK_INFO,
#endif
#if NO_EDIT_ICON
	GRETL_STOCK_EDIT,    /* -> edit_24_xpm, 24 x 24 */
	GRETL_STOCK_SCRIPT,  /* -> mini_edit_xpm, 16 x 16 */  
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
	GRETL_STOCK_PIN,
	GRETL_STOCK_ALPHA,
	GRETL_STOCK_EN,
	GRETL_STOCK_SPLIT
    };
    int n = G_N_ELEMENTS(stocks);

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

/* callbacks for viewer window toolbar */

void save_as_callback (GtkWidget *w, windata_t *vwin)
{
    guint u = 0;

    if (g_object_get_data(G_OBJECT(vwin->main), "text_out")) {
	const char *opts[] = {
	    N_("Save to file"),
	    N_("Save to session as icon")
	};
	int resp;

	resp = radio_dialog(_("gretl: save text"), _("Save text"), 
			    opts, 2, 0, 0);
	if (resp < 0) {
	    return;
	} else if (resp == 1) {
	    save_output_as_text_icon(vwin);
	    return;
	} else {
	    u = SAVE_OUTPUT;
	}
    } else if (vwin->role == EDIT_SCRIPT ||
	       vwin->role == VIEW_SCRIPT ||
	       vwin->role == VIEW_LOG ||
	       vwin->role == VIEW_PKG_CODE) {
	u = SAVE_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	u = SAVE_GP_CMDS;
    } else if (vwin->role == EDIT_R) {
	u = SAVE_R_CMDS;
    } else if (vwin->role == EDIT_OX) {
	u = SAVE_OX_CMDS;
    } else if (vwin->role == VIEW_FILE) {
	u = SAVE_TEXT;
    } else {
	dummy_call();
	return;
    }

    file_save(vwin, u);
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

static void file_open_callback (GtkWidget *w, windata_t *vwin)
{
    if (query_save_text(NULL, NULL, vwin)) {
	return;
    }

    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);

    if (vwin->flags & VWIN_CONTENT_CHANGED) {
	mark_vwin_content_saved(vwin);
	vwin_set_filename(vwin, tryfile);
    }
}

static void toolbar_new_callback (GtkWidget *w, windata_t *vwin)
{
    do_new_script(vwin->role);
}

#ifdef NATIVE_PRINTING
static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    window_print(NULL, vwin);
}
#endif

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
    display_files(PS_FILES, NULL);
}

static void cmd_log_refresh (GtkWidget *w, windata_t *vwin)
{
    gchar *newtext;
    int err = 0;

    newtext = get_logfile_content(&err);

    if (err) {
	gui_errmsg(err);
    } else {
	GtkTextBuffer *buf;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	gtk_text_buffer_set_text(buf, "", -1);
	if (newtext != NULL) {
	    textview_set_text(vwin->text, newtext);
	    g_free(newtext);
	}
    }
}

static void matrix_savename (GtkWidget *w, dialog_t *dlg)
{
    char *newname = (char *) edit_dialog_get_data(dlg);
    const gchar *buf = edit_dialog_get_text(dlg);

    if (buf == NULL || gui_validate_varname(buf, GRETL_TYPE_MATRIX)) {
	return;
    }

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

static void real_coeffint_set_alpha (GtkWidget *w, GtkWidget *dialog)
{
    windata_t *vwin = g_object_get_data(G_OBJECT(dialog), "vwin");
    double *x = g_object_get_data(G_OBJECT(dialog), "xptr");
    CoeffIntervals *cf = vwin->data;
    GtkTextBuffer *buf;
    const char *newtext;
    PRN *prn;

    if (bufopen(&prn)) {
	return;
    }

    reset_coeff_intervals(cf, 1.0 - *x);
    text_print_model_confints(cf, prn);
    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    gtk_text_buffer_set_text(buf, "", -1);
    textview_set_text(vwin->text, newtext);
    gretl_print_destroy(prn); 

    gtk_widget_destroy(dialog);
}

static void alpha_button_callback (GtkToggleButton *b, double *x)
{
    if (gtk_toggle_button_get_active(b)) {
	int i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(b), "i"));

	if (i == 0) {
	    *x = 0.90;
	} else if (i == 1) {
	    *x = 0.95;
	} else if (i == 2) {
	    *x = 0.99;
	}
    } 
}

static void reformat_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_MODELTABLE) {
	format_model_table(vwin);
    } else {
	series_view_format_dialog(w, vwin);
    }
}

static void split_pane_callback (GtkWidget *w, windata_t *vwin)
{
    if (g_object_get_data(G_OBJECT(vwin->vbox), "sw") != NULL) {
	/* currently in single-view mode */
	viewer_split_pane(w, vwin);
    } else {
	viewer_close_pane(w, vwin);
    }
}

static void toggle_alpha_spin (GtkToggleButton *b, GtkWidget *w)
{
    gtk_widget_set_sensitive(w, gtk_toggle_button_get_active(b));
}

static void coeffint_set_alpha (GtkWidget *w, windata_t *vwin)
{
    CoeffIntervals *cf = vwin->data;
    GtkWidget *dialog, *tmp, *hbox;
    GtkWidget *vbox, *b, *hb2;
    GSList *group = NULL;
    GtkObject *adj;
    gchar txt[16];
    double x = 1.0 - cf->alpha;
    gboolean defset = FALSE;
    int i;

    if (maybe_raise_dialog()) {
	return;
    }

    dialog = gretl_dialog_new(_("gretl: coefficient confidence intervals"), 
			      vwin->main, GRETL_DLG_BLOCK);

    hbox = gtk_hbox_new(FALSE, 5);
    hb2 = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Confidence level"));
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 0);

    tmp = gtk_label_new("1 - α :");
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), hb2, TRUE, TRUE, 10);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    vbox = gtk_vbox_new(FALSE, 5);

    /* radio button for 90%, 95%, 99% confidence */

    for (i=0; i<3; i++) {
	double a = (i == 0)? 0.90 : (i == 1)? 0.95 : 0.99;

	sprintf(txt, "%.2f", a);
	b = gtk_radio_button_new_with_label(group, txt);
	gtk_box_pack_start(GTK_BOX(vbox), b, FALSE, FALSE, 0);
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
	if (a == x) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), TRUE);
	    defset = TRUE;
	} else {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), FALSE);
	}
	g_object_set_data(G_OBJECT(b), "i", GINT_TO_POINTER(i));
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(alpha_button_callback), 
			 &x);
    }

    /* radio button for "other" confidence level, plus spinner */

    hb2 = gtk_hbox_new(FALSE, 0);
    b = gtk_radio_button_new_with_label(group, _("Other"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), !defset);
    gtk_box_pack_start(GTK_BOX(hb2), b, FALSE, FALSE, 0);
    adj = gtk_adjustment_new(x, 0.60, 0.99, 0.01, 0, 0);
    tmp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.01, 2);
    g_signal_connect(tmp, "value-changed", 
		     G_CALLBACK(set_double_from_spinner), &x);
    gtk_widget_set_sensitive(tmp, !defset);
    g_signal_connect(G_OBJECT(b), "toggled", 
		     G_CALLBACK(toggle_alpha_spin),
		     tmp);
    gtk_entry_set_activates_default(GTK_ENTRY(tmp), TRUE);
    gtk_box_pack_start(GTK_BOX(hb2), tmp, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hb2, FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 10);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* Cancel button */
    cancel_options_button(hbox, dialog, NULL);

    g_object_set_data(G_OBJECT(dialog), "vwin", vwin);
    g_object_set_data(G_OBJECT(dialog), "xptr", &x);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(real_coeffint_set_alpha), dialog);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);
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

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor(vwin->text, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
}

static void delete_file_viewer (GtkWidget *widget, windata_t *vwin) 
{
    gint resp = 0;

    if (window_is_busy(vwin)) {
	maybe_raise_dialog();
	return;
    }

    if (vwin_is_editing(vwin) && vwin_content_changed(vwin)) {
	resp = query_save_text(NULL, NULL, vwin);
    }

    if (!resp) {
	gtk_widget_destroy(vwin->main); 
    }
}

static int edit_script_popup_item (GretlToolItem *item)
{
    return !strcmp(item->icon, GTK_STOCK_COPY) ||
	!strcmp(item->icon, GTK_STOCK_PASTE) ||
	!strcmp(item->icon, GTK_STOCK_FIND) ||
	!strcmp(item->icon, GTK_STOCK_UNDO) ||
	!strcmp(item->icon, GTK_STOCK_FIND_AND_REPLACE);
}

static void set_plot_icon (GretlToolItem *item)
{
    if (dataset_is_time_series(datainfo)) {
	item->icon = GRETL_STOCK_TS;
    } else {
	item->icon = GRETL_STOCK_BOX;
    }
}

static GretlToolItem viewbar_items[] = {
    { N_("New window"), GTK_STOCK_NEW, G_CALLBACK(toolbar_new_callback), OPEN_ITEM },
    { N_("Open..."), GTK_STOCK_OPEN, G_CALLBACK(file_open_callback), OPEN_ITEM },
    { N_("Save"), GTK_STOCK_SAVE, G_CALLBACK(vwin_save_callback), SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, G_CALLBACK(save_as_callback), SAVE_AS_ITEM },
#ifdef NATIVE_PRINTING
    { N_("Print..."), GTK_STOCK_PRINT, G_CALLBACK(window_print_callback), 0 },
#endif
    { N_("Run"), GTK_STOCK_EXECUTE, G_CALLBACK(do_run_script), EXEC_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, G_CALLBACK(vwin_copy_callback), COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, G_CALLBACK(text_paste), EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, G_CALLBACK(text_find), 0 },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, G_CALLBACK(text_replace), EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, G_CALLBACK(text_undo), EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(series_view_toggle_sort), SORT_ITEM },    
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(multi_series_view_sort_by), SORT_BY_ITEM },
    { N_("Configure tabs..."), GTK_STOCK_PREFERENCES, G_CALLBACK(script_tabs_dialog), EDIT_SCRIPT_ITEM },
    { N_("Send To..."), GRETL_STOCK_MAIL, G_CALLBACK(mail_script_callback), MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, G_CALLBACK(script_index), INDEX_ITEM },
    { N_("Confidence level..."), GRETL_STOCK_ALPHA, G_CALLBACK(coeffint_set_alpha), ALPHA_ITEM },
    { N_("Refresh"), GTK_STOCK_REFRESH, G_CALLBACK(cmd_log_refresh), REFRESH_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, G_CALLBACK(window_tex_callback), TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, G_CALLBACK(series_view_graph), PLOT_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, G_CALLBACK(reformat_callback), FORMAT_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, G_CALLBACK(add_data_callback), ADD_DATA_ITEM },
    { N_("Add as matrix..."), GTK_STOCK_ADD, G_CALLBACK(add_matrix_callback), ADD_MATRIX_ITEM },
    { N_("Stickiness..."), GRETL_STOCK_PIN, G_CALLBACK(set_output_sticky), STICKIFY_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT, G_CALLBACK(split_pane_callback), SPLIT_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, G_CALLBACK(activate_script_help), CMD_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(window_help), HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(display_gnuplot_help), GP_HELP_ITEM },
    { N_("Close"), GTK_STOCK_CLOSE, G_CALLBACK(delete_file_viewer), 0 }
};

static int n_viewbar_items = G_N_ELEMENTS(viewbar_items);

#define exec_ok(r) (vwin_editing_script(r) || \
		    r == VIEW_SCRIPT || \
	            r == EDIT_PKG_SAMPLE)

#define open_ok(r) (vwin_editing_script(r))

#define edit_ok(r) (vwin_editing_script(r) || \
		    vwin_editing_buffer(r) || \
                    r == EDIT_PKG_CODE || \
		    r == EDIT_PKG_SAMPLE)

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_PKG_CODE && \
		       r != EDIT_PKG_SAMPLE && \
	               r != CONSOLE)

#define help_ok(r) (r == LEVERAGE || \
		    r == COINT2 || \
		    r == HURST || \
		    r == RMPLOT || \
		    r == MAHAL)

#define cmd_help_ok(r) (r == EDIT_SCRIPT || \
	                r == EDIT_PKG_CODE || \
                        r == EDIT_PKG_SAMPLE || \
			r == VIEW_SCRIPT || \
			r == VIEW_LOG)

#define sort_ok(r)    (r == VIEW_SERIES)
#define plot_ok(r)    (r == VIEW_SERIES)

#define add_data_ok(r) (r == PCA || r == LEVERAGE || \
                        r == MAHAL || r == FCAST)

#define split_ok(r) (r == SCRIPT_OUT)

/* Screen out unwanted menu items depending on the context; also
   adjust the callbacks associated with some items based on
   context.
*/

static GCallback item_get_callback (GretlToolItem *item, windata_t *vwin, 
				    int latex_ok, int sortby_ok,
				    int format_ok, int save_ok)
{
    GCallback func = item->func;
    int f = item->flag;
    int r = vwin->role;

    if (!edit_ok(r) && f == EDIT_ITEM) {
	return NULL;
    } else if (!open_ok(r) && f == OPEN_ITEM) {
	return NULL;
    } else if (!exec_ok(r) && f == EXEC_ITEM) {
	return NULL;
    } else if (!cmd_help_ok(r) && f == CMD_HELP_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && f == MAIL_ITEM) {
	return NULL;
    } else if (!help_ok(r) && f == HELP_ITEM) {
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
    } else if (!split_ok(r) && f == SPLIT_ITEM) {
	return NULL;
    } else if (!format_ok && f == FORMAT_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && r != EDIT_PKG_CODE && 
	       r != EDIT_PKG_SAMPLE && f == EDIT_SCRIPT_ITEM) {
	return NULL;
    } else if (r != VIEW_SCRIPT && f == INDEX_ITEM) {
	return NULL;
    } else if (r != SCRIPT_OUT && f == STICKIFY_ITEM) {
	return NULL;
    } else if (r != COEFFINT && f == ALPHA_ITEM) {
	return NULL;
    } else if (r != VIEW_LOG && f == REFRESH_ITEM) {
	return NULL;
    } else if (r != EDIT_GP && f == GP_HELP_ITEM) {
	return NULL;
    } else if (f == SAVE_ITEM && !save_ok) {
	return NULL;
    } else if (f == SAVE_AS_ITEM) {
	if (!save_as_ok(r)) {
	    return NULL;
	} else if (MULTI_FORMAT_ENABLED(r) || (r == PRINT && vwin->data != NULL)) {
	    func = G_CALLBACK(multi_save_as_callback);
	}
    }

    return func;
}

GtkWidget *gretl_toolbar_new (void)
{
    GtkWidget *tb;

    tb = gtk_toolbar_new();
    gtk_toolbar_set_icon_size(GTK_TOOLBAR(tb), GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(tb), FALSE);

    return tb;
}

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 12)

static GtkTooltips *gretl_tips;

void gretl_tooltips_init (void)
{
    gretl_tips = gtk_tooltips_new();
    gtk_tooltips_enable(gretl_tips); /* redundant? */
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_tooltips_set_tip(gretl_tips, w, str, NULL);
}

#else /* new tooltips API */

void gretl_tooltips_init (void)
{
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_widget_set_tooltip_text(w, str);
}

#endif

GtkToolItem *gretl_toolbar_insert (GtkWidget *tbar,
				   GretlToolItem *item,
				   GCallback func,
				   gpointer data,
				   gint pos)
{
    GtkToolItem *button;

    button = gtk_tool_button_new_from_stock(item->icon);
#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 12)
    gtk_tool_item_set_tooltip(button, gretl_tips, _(item->tip), NULL);
#else
    gtk_widget_set_tooltip_text(GTK_WIDGET(button), _(item->tip));
#endif
    g_signal_connect(button, "clicked", func, data);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), button, pos);

    return button;
}

static void viewbar_add_items (windata_t *vwin, ViewbarFlags flags)
{
    int sortby_ok = has_sortable_data(vwin);
    int format_ok = can_format_data(vwin);
    int latex_ok = latex_is_ok();
    int save_ok = (flags & VIEWBAR_EDITABLE);
    GtkToolItem *button;
    GretlToolItem *item;
    GCallback func;
    int i;
 
    for (i=0; i<n_viewbar_items; i++) {
	item = &viewbar_items[i];

	func = item_get_callback(item, vwin, latex_ok, sortby_ok,
				 format_ok, save_ok);
	if (func == NULL) {
	    continue;
	}

	if (item->flag == PLOT_ITEM) {
	    set_plot_icon(item);
	}

	button = gretl_toolbar_insert(vwin->mbar, item, func, vwin, -1);

	if (item->flag == SAVE_ITEM) { 
	    if (vwin->role != CONSOLE) {
		/* nothing to save just yet */
		g_object_set_data(G_OBJECT(vwin->mbar), "save_button", button); 
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	    }
	} else if (item->flag == SAVE_AS_ITEM) {
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button", button);
	    if (strstr(vwin->fname, "script_tmp")) {
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	    }
	}
    }
}

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags)
{
    GtkWidget *hbox;

    if (gretl_stock_ifac == NULL) {
	gretl_stock_icons_init();
    }

    if ((flags & VIEWBAR_HAS_TEXT) || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->main), "text_out", GINT_TO_POINTER(1));
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);

    vwin->mbar = gretl_toolbar_new();

    viewbar_add_items(vwin, flags);

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

static void remove_child (GtkWidget *child, GtkWidget *cont)
{
    gtk_container_remove(GTK_CONTAINER(cont), child);
}

void viewbar_add_edit_items (windata_t *vwin)
{
    gtk_container_foreach(GTK_CONTAINER(vwin->mbar),
			  (GtkCallback) remove_child, vwin->mbar);
    viewbar_add_items(vwin, VIEWBAR_EDITABLE);
    gtk_widget_show_all(vwin->mbar);
}

GtkWidget *build_text_popup (windata_t *vwin)
{
    GtkWidget *pmenu = gtk_menu_new();
    GretlToolItem *item;
    GCallback func;
    GtkWidget *w;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	item = &viewbar_items[i];
	if (vwin->role == EDIT_SCRIPT) {
	    /* the script editor popup may have some special stuff
	       added: don't clutter it up */
	    if (edit_script_popup_item(item)) {
		func = item->func;
	    } else {
		func = NULL;
	    }
	} else {
	    func = item_get_callback(item, vwin, 0, 0, 0, 0);
	}
	if (func != G_CALLBACK(NULL)) {
	    if (func == G_CALLBACK(text_paste)) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);

		if (!gtk_clipboard_wait_is_text_available(cb)) {
		    continue;
		}
	    } else if (func == G_CALLBACK(text_undo) && !text_can_undo(vwin)) {
		continue;
	    }
	    w = gtk_menu_item_new_with_label(_(item->tip));
	    g_signal_connect(G_OBJECT(w), "activate", func, vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	}
    }

    if (vwin->role == SCRIPT_OUT) {
	if (g_object_get_data(G_OBJECT(vwin->vbox), "vpaned") == NULL) {
	    w = gtk_menu_item_new_with_label(_("Split pane"));
	    g_signal_connect(G_OBJECT(w), "activate", G_CALLBACK(viewer_split_pane), vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	} else {
	    w = gtk_menu_item_new_with_label(_("Close pane"));
	    g_signal_connect(G_OBJECT(w), "activate", G_CALLBACK(viewer_close_pane), vwin);
	    gtk_widget_show(w);
	    gtk_menu_shell_append(GTK_MENU_SHELL(pmenu), w);
	}
    }

    return pmenu;
}

/* callbacks for main-window toolbar icons */

static void tbar_calc (void)
{
#ifdef G_OS_WIN32
    create_child_process(calculator);
#else
    gretl_fork("calculator", NULL);
#endif 
}

static void tbar_open_data (void)
{
    display_files(TEXTBOOK_DATA, NULL);
}

static void tbar_users_guide (void)
{
    display_pdf_help(NULL);
}

static void tbar_command_ref (void)
{
    plain_text_cmdref(NULL);
}

static void tbar_xy_graph (void)
{
    if (data_status) {
	if (datainfo->v == 2) {
	    do_graph_var(mdata->active_var);
	} else if (mdata_selection_count() == 2) {
	    plot_from_selection(GR_XY);
	} else {
	    selection_dialog(_("gretl: define graph"), 
			     do_graph_from_selector,
			     GR_XY);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_model (void)
{
    if (data_status) {
	selection_dialog(_("gretl: specify model"), do_model, OLS);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_iconview (void)
{
    if (data_status) {
	view_session(NULL);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_new_script (void)
{
    do_new_script(EDIT_SCRIPT);
}

static void tbar_show_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, mdata);
}

/* end toolbar icon callbacks */

static GretlToolItem mainbar_items[] = {
    { N_("launch calculator"),  GRETL_STOCK_CALC,    G_CALLBACK(tbar_calc), 0 },
#if NO_EDIT_ICON
    { N_("new script"),         GRETL_STOCK_SCRIPT,  G_CALLBACK(tbar_new_script), 0 },
#else
    { N_("new script"),         GTK_STOCK_EDIT,      G_CALLBACK(tbar_new_script), 0 },
#endif
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, G_CALLBACK(gretl_console), 0 },
    { N_("session icon view"),  GRETL_STOCK_ICONS,   G_CALLBACK(tbar_iconview), 0 },
    { N_("function packages"),  GRETL_STOCK_FUNC,    G_CALLBACK(tbar_show_funcs), 0 },
    { N_("user's guide"),       GRETL_STOCK_PDF,     G_CALLBACK(tbar_users_guide), 0 },
    { N_("command reference"),  GTK_STOCK_HELP,      G_CALLBACK(tbar_command_ref), 0 },
    { N_("X-Y graph"),          GRETL_STOCK_SCATTER, G_CALLBACK(tbar_xy_graph), 0 },
    { N_("OLS model"),          GRETL_STOCK_MODEL,   G_CALLBACK(tbar_model), 0 },
    { N_("open dataset"),       GTK_STOCK_OPEN,      G_CALLBACK(tbar_open_data), 0 }
};

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *hbox, *toolbar;
    GretlToolItem *item;
    int i, n = G_N_ELEMENTS(mainbar_items);

    gretl_stock_icons_init();

    toolbar = gretl_toolbar_new();

    for (i=0; i<n; i++) {
	item = &mainbar_items[i];
	gretl_toolbar_insert(toolbar, item, item->func, mdata, -1);
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

void show_toolbar (void)
{
    GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->main), "vbox");

    make_toolbar(vbox);
}

