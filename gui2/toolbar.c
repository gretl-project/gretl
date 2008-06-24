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

/* toolbar.c: the gretl toolbar */

#include "gretl.h"
#include "console.h"
#include "session.h"
#include "datafiles.h"
#include "selector.h"
#include "textbuf.h"
#include "textutil.h"
#include "series_view.h"
#include "cmdstack.h"
#include "dlgutils.h"

#include "usermat.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include "../pixmaps/mini.tex.xpm"
#include "../pixmaps/mail_16.xpm"
#include "../pixmaps/mini.tsplot.xpm"
#include "../pixmaps/mini.boxplot.xpm"
#include "../pixmaps/mini.pdf.xpm"
#include "../pixmaps/mini.manual.xpm"
#include "../pixmaps/mini.pin.xpm"
#include "../pixmaps/mini.alpha.xpm"
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
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    INDEX_ITEM,
    EDIT_SCRIPT_ITEM,
    STICKIFY_ITEM,
    ALPHA_ITEM,
    REFRESH_ITEM
} viewbar_flags;

typedef void (*toolfunc) (GtkWidget *w, windata_t *vwin);

struct viewbar_item {
    const char *tip;
    const gchar *icon;
    toolfunc func;
    int flag;
};

struct toolbar_item {
    const char *tip;
    const gchar *icon;
    void (*toolfunc)();
};

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
	mini_alpha_xpm
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
	GRETL_STOCK_PIN,
	GRETL_STOCK_ALPHA
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

/* callbacks for viewer window toolbar */

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
	window_copy(vwin, GRETL_FORMAT_SELECTION);
    } else if (vwin->role == VIEW_SCALAR) {
	scalar_to_clipboard(vwin);
    } else if (!script_role(vwin->role) && !editor_role(vwin->role)) {
	copy_format_dialog(vwin, W_COPY);
    } else {
	window_copy(vwin, GRETL_FORMAT_TXT);
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

static void save_plot_commands_callback (GtkWidget *w, windata_t *vwin)
{
    auto_save_plot(vwin);
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
    gchar *logfile = NULL;
    gchar *newtext = NULL;
    int err;

    logfile = g_strdup_printf("%ssession.inp", paths.dotdir);
    err = dump_command_stack(logfile, 0);

    if (!err) {
	err = gretl_file_get_contents(logfile, &newtext);
    }

    if (!err) {
	GtkTextBuffer *buf;

	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	gtk_text_buffer_set_text(buf, "", -1);
	textview_set_text(vwin->w, newtext);
	g_free(newtext);
    }

    g_free(logfile);
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

#define xround(x) (((x-floor(x))>.5)? ceil(x) : floor(x))

static void coeffint_set_alpha (GtkWidget *w, windata_t *vwin)
{
    CoeffIntervals *cf = vwin->data;
    GtkTextBuffer *buf;
    const char *newtext;
    double alpha, x = 100 * (1 - cf->alpha);
    int cval = (int) xround(x);
    PRN *prn;
    int resp;

    resp = spin_dialog("gretl: alpha", NULL,
		       &cval, _("Confidence level, percent"),
		       60, 99, 0);

    if (resp < 0 || bufopen(&prn)) {
	return;
    }

    alpha = (100.0 - cval) / 100.0;
    reset_coeff_intervals(cf, alpha);
    text_print_model_confints(cf, prn);
    newtext = gretl_print_get_buffer(prn);
    buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));

    gtk_text_buffer_set_text(buf, "", -1);
    textview_set_text(vwin->w, newtext);
    gretl_print_destroy(prn);
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
    text_set_cursor(vwin->w, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
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
	(vwin->flags & VWIN_CONTENT_CHANGED)) {
	resp = query_save_text(NULL, NULL, vwin);
    }

    if (!resp) {
	gtk_widget_destroy(vwin->dialog); 
    }
}

static int edit_script_popup_item (struct viewbar_item *item)
{
    return !strcmp(item->icon, GTK_STOCK_COPY) ||
	!strcmp(item->icon, GTK_STOCK_PASTE) ||
	!strcmp(item->icon, GTK_STOCK_FIND) ||
	!strcmp(item->icon, GTK_STOCK_UNDO) ||
	!strcmp(item->icon, GTK_STOCK_FIND_AND_REPLACE);
}

static void set_plot_icon (struct viewbar_item *vitem)
{
    if (dataset_is_time_series(datainfo)) {
	vitem->icon = GRETL_STOCK_TS;
    } else {
	vitem->icon = GRETL_STOCK_BOX;
    }
}

static struct viewbar_item viewbar_items[] = {
    { N_("Save"), GTK_STOCK_SAVE, view_window_save, SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, file_save_callback, SAVE_AS_ITEM },
#ifdef NATIVE_PRINTING
    { N_("Print..."), GTK_STOCK_PRINT, window_print_callback, 0 },
#endif
    { N_("Run"), GTK_STOCK_EXECUTE, do_run_script, EXEC_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, text_copy_callback, COPY_ITEM }, 
    { N_("Paste"), GTK_STOCK_PASTE, text_paste, EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, (toolfunc) text_find, 0 },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, text_replace, EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, text_undo, EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, series_view_sort, SORT_ITEM },    
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, series_view_sort_by, SORT_BY_ITEM },
    { N_("Configure tabs..."), GTK_STOCK_PREFERENCES, script_tabs_dialog, EDIT_SCRIPT_ITEM },
    { N_("Send To..."), GRETL_STOCK_MAIL, mail_script_callback, MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, script_index, INDEX_ITEM },
    { N_("Confidence level..."), GRETL_STOCK_ALPHA, coeffint_set_alpha, ALPHA_ITEM },
    { N_("Refresh"), GTK_STOCK_REFRESH, cmd_log_refresh, REFRESH_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, window_tex_callback, TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, series_view_graph, PLOT_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, series_view_format_dialog, FORMAT_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, add_data_callback, ADD_DATA_ITEM },
    { N_("Add as matrix..."), GTK_STOCK_ADD, add_matrix_callback, ADD_MATRIX_ITEM },
    { N_("Stickiness..."), GRETL_STOCK_PIN, set_output_sticky, STICKIFY_ITEM },
    { N_("Help on command"), GTK_STOCK_HELP, activate_script_help, CMD_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, window_help, HELP_ITEM },
    { N_("Close"), GTK_STOCK_CLOSE, delete_file_viewer, 0 }
};

static int n_viewbar_items = G_N_ELEMENTS(viewbar_items);

#define exec_ok(r) (r == EDIT_SCRIPT || \
                    r == EDIT_GP || \
                    r == EDIT_R || \
	            r == VIEW_SCRIPT)

#define edit_ok(r) (r == EDIT_SCRIPT || \
                    r == EDIT_HEADER || \
                    r == EDIT_NOTES || \
                    r == EDIT_FUNC_CODE || \
	            r == EDIT_GP || \
                    r == EDIT_BOX || \
		    r == EDIT_R || \
                    r == SCRIPT_OUT)

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_FUNC_CODE && \
		       r != EDIT_BOX && \
		       r != VIEW_SCALAR)

#define help_ok(r) (r == LEVERAGE || \
		    r == COINT2 || \
		    r == HURST || \
		    r == RMPLOT || \
		    r == MAHAL)

#define cmd_help_ok(r) (r == EDIT_SCRIPT || \
			r == VIEW_SCRIPT || \
			r == VIEW_LOG)

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
    } else if (!exec_ok(r) && f == EXEC_ITEM) {
	return NULL;
    } else if (!cmd_help_ok(r) && f == CMD_HELP_ITEM) {
	return NULL;
    } else if (r != EDIT_SCRIPT && f == MAIL_ITEM) {
	return NULL;
    } else if (!help_ok(r) && f == HELP_ITEM) {
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
    } else if (r != COEFFINT && f == ALPHA_ITEM) {
	return NULL;
    } else if (r != VIEW_LOG && f == REFRESH_ITEM) {
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

void vwin_add_viewbar (windata_t *vwin, int text_out)
{
    GtkWidget *hbox;
    GtkToolItem *button;
    toolfunc func;
    const char *tooltip;
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
    gtk_toolbar_set_icon_size(GTK_TOOLBAR(vwin->mbar), GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(vwin->mbar), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(vwin->mbar), FALSE);

    for (i=0; i<n_viewbar_items; i++) {
	struct viewbar_item *vitem = &viewbar_items[i];

	func = item_get_callback(vitem, vwin, latex_ok, sortby_ok);
	if (func == NULL) {
	    continue;
	}

	if (vitem->flag == PLOT_ITEM) {
	    set_plot_icon(vitem);
	}

	tooltip = vitem->tip;

	if (vitem->flag == EXEC_ITEM) {
	    if (vwin->role == EDIT_GP) {
		tooltip = N_("Send to gnuplot");
	    } else if (vwin->role == EDIT_R) {
		tooltip = N_("Send to R");
	    }
	}

	button = gtk_tool_button_new_from_stock(vitem->icon);
	gtk_tool_item_set_tooltip(button, get_gretl_tips(), 
				  _(tooltip), NULL);
	g_signal_connect(button, "clicked", G_CALLBACK(func), vwin);
	gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), button, -1);

	if (vitem->flag == SAVE_ITEM) { 
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_button", button); 
	    /* nothing to save just yet */
	    gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	} else if (vitem->flag == SAVE_AS_ITEM) {
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button", button);
	    if (strstr(vwin->fname, "script_tmp")) {
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
	    }
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), vwin->mbar, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

void viewbar_add_edit_items (windata_t *vwin)
{
    struct viewbar_item *vitem;
    GtkToolItem *button;
    int i, pos = 0;

    for (i=0; i<n_viewbar_items; i++) {
	vitem = &viewbar_items[i];
	if (vitem->flag == SAVE_ITEM ||
	    vitem->flag == EDIT_ITEM ||
	    vitem->flag == EDIT_SCRIPT_ITEM) {
	    button = gtk_tool_button_new_from_stock(vitem->icon);
	    gtk_tool_item_set_tooltip(button, get_gretl_tips(), 
				      _(vitem->tip), NULL);
	    g_signal_connect(button, "clicked", G_CALLBACK(vitem->func), vwin);
	    gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), button, pos);
	    if (vitem->flag == SAVE_ITEM) {
		gtk_widget_set_sensitive(GTK_WIDGET(button), FALSE);
		g_object_set_data(G_OBJECT(vwin->mbar), "save_button",
				  button);
	    } else if (vitem->flag == SAVE_AS_ITEM) {
		g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button",
				  button);
	    }
	    gtk_widget_show(GTK_WIDGET(button));
	}
	if (vitem->flag != EDIT_GP) {
	    pos++;
	}
    }
}

GtkWidget *build_text_popup (windata_t *vwin)
{
    struct viewbar_item *vitem;
    toolfunc func;
    GtkWidget *pmenu = gtk_menu_new();
    GtkWidget *w;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	vitem = &viewbar_items[i];
	if (vwin->role == EDIT_SCRIPT) {
	    /* the script editor popup may have some special stuff
	       added: don't clutter it up */
	    if (edit_script_popup_item(vitem)) {
		func = vitem->func;
	    } else {
		func = NULL;
	    }
	} else {
	    func = item_get_callback(vitem, vwin, 0, 0);
	}
	if (func != NULL) {
	    if (func == text_paste) {
		GtkClipboard *cb = gtk_clipboard_get(GDK_NONE);

		if (!gtk_clipboard_wait_is_text_available(cb)) {
		    continue;
		}
	    } else if (func == text_undo && !text_can_undo(vwin)) {
		continue;
	    }
	    w = gtk_menu_item_new_with_label(_(vitem->tip));
	    g_signal_connect(G_OBJECT(w), "activate",
			     G_CALLBACK(func),
			     vwin);
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
			     GR_XY, 0);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_model (void)
{
    if (data_status) {
	selection_dialog(_("gretl: specify model"), do_model, OLS, 0);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_iconview (void)
{
    if (data_status) {
	view_session();
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
    display_files(FUNC_FILES, w);
}

/* end toolbar icon callbacks */



#if NO_EDIT_ICON
# define GTK_STOCK_EDIT "gretl-script"
#endif

static struct toolbar_item toolbar_items[] = {
    { N_("launch calculator"),  GRETL_STOCK_CALC,    tbar_calc },
    { N_("new script"),         GTK_STOCK_EDIT,      tbar_new_script },
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, show_gretl_console },
    { N_("session icon view"),  GRETL_STOCK_ICONS,   tbar_iconview },
    { N_("function packages"),  GRETL_STOCK_FUNC,    tbar_show_funcs },
    { N_("user's guide"),       GRETL_STOCK_PDF,     tbar_users_guide },
    { N_("command reference"),  GTK_STOCK_HELP,      tbar_command_ref },
    { N_("X-Y graph"),          GRETL_STOCK_SCATTER, tbar_xy_graph },
    { N_("OLS model"),          GRETL_STOCK_MODEL,   tbar_model },
    { N_("open dataset"),       GTK_STOCK_OPEN,      tbar_open_data }
};

static void make_toolbar (GtkWidget *vbox)
{
    GtkWidget *hbox, *toolbar;
    GtkToolItem *item;
    int i, n = G_N_ELEMENTS(toolbar_items);

    gretl_stock_icons_init();

    toolbar = gtk_toolbar_new();
    gtk_toolbar_set_icon_size(GTK_TOOLBAR(toolbar), GTK_ICON_SIZE_MENU);
    gtk_toolbar_set_style(GTK_TOOLBAR(toolbar), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(toolbar), FALSE);

    for (i=0; i<n; i++) {
	item = gtk_tool_button_new_from_stock(toolbar_items[i].icon);
	gtk_tool_item_set_tooltip(item, get_gretl_tips(), 
				  _(toolbar_items[i].tip),
				  NULL);
	g_signal_connect(item, "clicked", toolbar_items[i].toolfunc, mdata);
	gtk_toolbar_insert(GTK_TOOLBAR(toolbar), item, -1);
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), toolbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

void show_toolbar (void)
{
    GtkWidget *vbox = g_object_get_data(G_OBJECT(mdata->w), "vbox");

    make_toolbar(vbox);
}

