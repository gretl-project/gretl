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
#include "gui_utils.h"
#include "bundle_menus.h"
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
#include "winstack.h"
#include "tabwin.h"
#include "fncall.h"
#include "fnsave.h"
#include "database.h"
#include "toolbar.h"

#include "uservar.h"
#include "forecast.h"
#include "gretl_www.h"
#include "gretl_mdconv.h"

#ifdef G_OS_WIN32
# include "gretlwin32.h"
#endif

#include <gio/gio.h>

/* for viewer window toolbars */
#include "../pixmaps/mini.en.xpm"

/* for pop-up search entry */
#include "../pixmaps/close_16.xpm"

/* for window-finder menu */
#include "../pixmaps/mini.gretl.xpm"
#include "../pixmaps/mini.table.xpm"
#include "../pixmaps/mini.page.xpm"

/* for markdown editor */
#include "../pixmaps/eye.xpm"
#include "../pixmaps/eye-off.xpm"
#include "../pixmaps/markdown.xpm"

enum {
    SAVE_ITEM = 1,
    SAVE_AS_ITEM,
    EDIT_ITEM,
    PLOT_ITEM,
    EXEC_ITEM,
    COPY_ITEM,
    PRINT_ITEM,
    TEX_ITEM,
    ADD_DATA_ITEM,
    ADD_MATRIX_ITEM,
    MAIL_ITEM,
    HELP_ITEM,
    CMD_HELP_ITEM,
    GP_HELP_ITEM,
    X12A_HELP_ITEM,
    SORT_ITEM,
    SORT_BY_ITEM,
    FORMAT_ITEM,
    INDEX_ITEM,
    EDIT_HANSL_ITEM,
    STICKIFY_ITEM,
    ALPHA_ITEM,
    REFRESH_ITEM,
    OPEN_ITEM,
    SPLIT_H_ITEM,
    SPLIT_V_ITEM,
    EDITOR_ITEM,
    NOTES_ITEM,
    NEW_ITEM,
    BUNDLE_ITEM,
    FIND_ITEM,
    COPY_SCRIPT_ITEM,
    BUILD_ITEM,
    HMAP_ITEM,
    DIGITS_ITEM,
    DBN_ITEM,
    FCAST_ITEM,
    CLOSE_ITEM,
    CLEAR_ITEM
} viewbar_flags;

int toolbar_icon_size = GTK_ICON_SIZE_MENU;

enum {
    I_IN_MENU  = 1 << 0,
    I_HAS_DARK = 1 << 1
} icon_flags;

struct png_stock_maker {
    char *name;
    const char *id;
    gint8 flags;
};

struct png_stock_maker png_stocks[] = {
    { "calc",      GRETL_STOCK_CALC, 0 },
    { "database",  GRETL_STOCK_DB, 0 },
    { "fx",        GRETL_STOCK_FUNC, I_HAS_DARK },
    { "betahat",   GRETL_STOCK_MODEL, I_HAS_DARK },
    { "iconview",  GRETL_STOCK_ICONS, 0 },
    { "console",   GRETL_STOCK_CONSOLE, 0 },
    { "plot",      GRETL_STOCK_SCATTER, 0 },
    { "winlist",   GRETL_STOCK_WINLIST },
    { "bundle",    GRETL_STOCK_BUNDLE, I_IN_MENU },
    { "alpha",     GRETL_STOCK_ALPHA, 0 },
    { "tex",       GRETL_STOCK_TEX, 0 },
    { "bigger",    GRETL_STOCK_BIGGER, 0 },
    { "smaller",   GRETL_STOCK_SMALLER, 0 },
    { "menu",      GRETL_STOCK_MENU, 0 },
    { "tools",     GRETL_STOCK_TOOLS, 0 },
    { "dbnomics",  GRETL_STOCK_DBN, I_IN_MENU },
    { "fcast",     GRETL_STOCK_FCAST, 0 },
    { "heatmap",   GRETL_STOCK_HMAP, 0 },
    { "pushpin",   GRETL_STOCK_PIN, 0 },
    { "mail",      GRETL_STOCK_MAIL, I_IN_MENU },
    { "pdf",       GRETL_STOCK_PDF, I_IN_MENU },
    { "split_h",   GRETL_STOCK_SPLIT_H, 0 },
    { "join_h",    GRETL_STOCK_JOIN_H, 0 },
    { "split_v",   GRETL_STOCK_SPLIT_V, 0 },
    { "join_v",    GRETL_STOCK_JOIN_V, 0 },
    { "boxplot",   GRETL_STOCK_BOX, 0 },
    { "tsplot",    GRETL_STOCK_TS, 0 },
    { "coeffplot", GRETL_STOCK_CP, 0 },
    { "book",      GRETL_STOCK_BOOK, 0 },
    { "query",     GRETL_STOCK_QUERY, 0 },
};

struct xpm_stock_maker {
    char **xpm;
    const char *id;
};

#if GTK_MAJOR_VERSION == 3

static void try_auto_icon_sizing (int *bigger)
{
    GdkDisplay *display = gdk_display_get_default();
    GdkMonitor *monitor = gdk_display_get_primary_monitor(display);
    GdkScreen *screen = gdk_screen_get_default();

    if (monitor != NULL) {
	int mmw = gdk_monitor_get_width_mm(monitor);
	int pxw = gdk_screen_get_width(screen);
	double mm16 = 16 * mmw / (double) pxw;

	if (mm16 < 2.8) {
	    /* size of 16 pixels in millimeters */
	    fprintf(stderr, " auto-setting larger icons\n");
	    *bigger = 1;
	}
    }
}

#endif /* GTK3 */

void gretl_stock_icons_init (void)
{
    struct xpm_stock_maker xpm_stocks[] = {
	{mini_en_xpm, GRETL_STOCK_EN},
	{mini_gretl_xpm, GRETL_STOCK_GRETL},
	{mini_table_xpm, GRETL_STOCK_TABLE},
	{mini_page_xpm, GRETL_STOCK_PAGE},
	{close_16_xpm, GRETL_STOCK_CLOSE},
	{eye_xpm, GRETL_STOCK_EYE},
	{eye_off_xpm, GRETL_STOCK_EYE_OFF},
	{markdown_xpm, GRETL_STOCK_MD}
    };
    static GtkIconFactory *gretl_factory;
    int n1 = G_N_ELEMENTS(png_stocks);
    int n2 = G_N_ELEMENTS(xpm_stocks);

    if (gretl_factory == NULL) {
	int bigger = (get_icon_sizing() == ICON_SIZE_MEDIUM);
	char pngname[16], icon_path[48], menu_path[48];
	gchar *p, *pm, *respath;
	GResource *icons;
	GtkIconSource *isrc;
	GtkIconSet *iset;
	GdkPixbuf *pbuf;
	int dark = 0;
	int i;

#if GTK_MAJOR_VERSION == 3
	if (get_icon_sizing() == ICON_SIZE_AUTO) {
	    try_auto_icon_sizing(&bigger);
	}
#endif

	if (bigger) {
#if GTK_MAJOR_VERSION == 3
	    toolbar_icon_size = GTK_ICON_SIZE_LARGE_TOOLBAR;
#else
	    toolbar_icon_size = GTK_ICON_SIZE_SMALL_TOOLBAR;
#endif
	}

	gretl_factory = gtk_icon_factory_new();

	respath = g_strdup_printf("%sgretl-icons.gresource", gretl_home());
	icons = g_resource_load(respath, NULL);
	if (icons == NULL) {
	    fprintf(stderr, "g_resource_load: failed to load icons\n");
	    g_free(respath);
	    goto do_pixmaps;
	}

	g_resources_register(icons);

	if (bigger) {
	    strcpy(icon_path, "/gretl/icons/24x24/");
	    strcpy(menu_path, "/gretl/icons/16x16/");
	    pm = strrchr(menu_path, '/') + 1;
	} else {
	    strcpy(icon_path, "/gretl/icons/16x16/");
	}
	p = strrchr(icon_path, '/') + 1;

	dark = dark_theme_active();

	for (i=0; i<n1; i++) {
	    if (dark && (png_stocks[i].flags & I_HAS_DARK)) {
		sprintf(pngname, "%s_dk.png", png_stocks[i].name);
	    } else {
		sprintf(pngname, "%s.png", png_stocks[i].name);
	    }
	    strcat(icon_path, pngname);
	    pbuf = gdk_pixbuf_new_from_resource(icon_path, NULL);
	    if (pbuf == NULL) {
		fprintf(stderr, "Failed to load %s\n", icon_path);
		*p = '\0';
		continue;
	    }
	    if (bigger && (png_stocks[i].flags & I_IN_MENU)) {
		iset = gtk_icon_set_new();
		/* for toolbar use */
		isrc = gtk_icon_source_new();
		gtk_icon_source_set_pixbuf(isrc, pbuf);
		gtk_icon_source_set_size(isrc, toolbar_icon_size);
		gtk_icon_source_set_size_wildcarded(isrc, FALSE);
		gtk_icon_set_add_source(iset, isrc);
		g_object_unref(pbuf);
		/* for menu use */
		strcat(menu_path, pngname);
		pbuf = gdk_pixbuf_new_from_resource(menu_path, NULL);
		isrc = gtk_icon_source_new();
		gtk_icon_source_set_pixbuf(isrc, pbuf);
		gtk_icon_source_set_size(isrc, GTK_ICON_SIZE_MENU);
		gtk_icon_source_set_size_wildcarded(isrc, FALSE);
		gtk_icon_set_add_source(iset, isrc);
		g_object_unref(pbuf);
		*pm = '\0';
	    } else {
		/* we just need a single icon */
		iset = gtk_icon_set_new_from_pixbuf(pbuf);
	    }
	    gtk_icon_factory_add(gretl_factory, png_stocks[i].id, iset);
	    gtk_icon_set_unref(iset);
	    *p = '\0';
	}

	g_free(respath);
	g_resources_unregister(icons);
	g_resource_unref(icons);

    do_pixmaps:

	for (i=0; i<n2; i++) {
	    pbuf = gdk_pixbuf_new_from_xpm_data((const char **) xpm_stocks[i].xpm);
	    iset = gtk_icon_set_new_from_pixbuf(pbuf);
	    g_object_unref(pbuf);
	    gtk_icon_factory_add(gretl_factory, xpm_stocks[i].id, iset);
	    gtk_icon_set_unref(iset);
	}

	gtk_icon_factory_add_default(gretl_factory);
    }
}

/* callbacks for viewer window toolbar */

static void copy_to_editor (GtkWidget *w, windata_t *vwin)
{
    gchar *buf = textview_get_text(vwin->text);

    if (vwin->role == VIEW_LOG) {
	/* allow for the possibility that the buffer is empty */
	gchar *s = buf;
	int n = 0;

	while (s != NULL && n < 3) {
	    s = strchr(s, '\n');
	    if (s != NULL) {
		s++;
		n++;
	    }
	}

	if (s != NULL) {
	    gchar *modbuf;

	    modbuf = g_strdup_printf("# logged commands\n%s", s);
	    do_new_script(EDIT_HANSL, modbuf, NULL);
	    g_free(modbuf);
	} else {
	    do_new_script(EDIT_HANSL, buf, NULL);
	}
    } else {
	do_new_script(EDIT_HANSL, buf, NULL);
    }

    g_free(buf);
}

static void save_as_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *vmain = vwin_toplevel(vwin);
    guint u = 0;

    if (g_object_get_data(G_OBJECT(vmain), "text_out")) {
	const char *opts[] = {
	    N_("Save to file"),
	    N_("Save to session as icon")
	};
	int resp;

	resp = radio_dialog(_("gretl: save text"), _("Save text"),
			    opts, 2, 0, 0, vmain);
	if (resp < 0) {
	    return;
	} else if (resp == 1) {
	    save_output_as_text_icon(vwin);
	    return;
	} else {
	    u = SAVE_OUTPUT;
	}
    } else if (vwin->role == EDIT_HANSL) {
	u = SAVE_SCRIPT;
    } else if (vwin->role == EDIT_GP) {
	u = SAVE_GP_CMDS;
    } else if (vwin->role == EDIT_R) {
	u = SAVE_R_CMDS;
    } else if (vwin->role == EDIT_OX) {
	u = SAVE_OX_CMDS;
    } else if (vwin->role == EDIT_OCTAVE) {
	u = SAVE_OCTAVE_CMDS;
    } else if (vwin->role == EDIT_PYTHON) {
        u = SAVE_PYTHON_CMDS;
    } else if (vwin->role == EDIT_JULIA) {
	u = SAVE_JULIA_CODE;
    } else if (vwin->role == EDIT_DYNARE) {
	u = SAVE_DYNARE_CODE;
    } else if (vwin->role == EDIT_LPSOLVE) {
	u = SAVE_LPSOLVE_CODE;
    } else if (vwin->role == EDIT_STATA) {
	u = SAVE_STATA_CMDS;
    } else if (vwin->role == EDIT_SPEC) {
	u = SAVE_SPEC_FILE;
    } else if (vwin->role == VIEW_FILE) {
	u = SAVE_TEXT;
    } else if (vwin->role == EDIT_PKG_HELP ||
	       vwin->role == EDIT_PKG_GHLP) {
	u = SAVE_HELP_TEXT;
    } else if (vwin->role == EDIT_X12A) {
	u = SAVE_X13_SPC;
    } else {
	dummy_call();
	return;
    }

    file_save(vwin, u);
}

/* Adjust the number of significant figures used in printing
   coefficients and standard errors in a model viewer window,
   or descriptive statistics in a "summary" window. Also record
   the number of digits for use via Copy.
*/

static void display_digits_callback (GtkWidget *w, windata_t *vwin)
{
    int wdigits = widget_get_int(vwin->text, "digits");
    int save_digits = get_gretl_digits();
    int digits, resp;

    digits = wdigits > 0 ? wdigits : save_digits;

    if (vwin->role == SUMMARY) {
	/* desciptive statistics window */
	Summary *s = vwin->data;
	int dmax = (s->opt & OPT_S)? 4 : 5;

	digits = (digits > dmax)? dmax : digits;
	resp = spin_dialog(NULL, NULL, &digits,
			   _("Number of digits to show for statistics"),
			   3, dmax, 0, vwin_toplevel(vwin));
    } else {
	/* model window */
	resp = spin_dialog(NULL, NULL, &digits,
			   _("Number of digits to show for coefficients"),
			   3, 6, 0, vwin_toplevel(vwin));
    }

    if (resp != GRETL_CANCEL && digits != wdigits) {
	const char *buf;
	int save_digits;
	PRN *prn;

	if (bufopen(&prn)) return;
	save_digits = get_gretl_digits();
	set_gretl_digits(digits);
	if (vwin->role == SUMMARY) {
	    print_summary(vwin->data, dataset, prn);
	} else {
	    printmodel(vwin->data, dataset, OPT_NONE, prn);
	}
	buf = gretl_print_get_trimmed_buffer(prn);
	textview_set_text(vwin->text, buf);
	gretl_print_destroy(prn);
	set_gretl_digits(save_digits);
	widget_set_int(vwin->text, "digits", digits);
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

    send_attachment(vwin->fname);
}

/* callback for the "Open" icon in a script editing window,
   which enables the user to switch to a different script,
   or to open another tab if the editor is tab-enabled
*/

static void file_open_callback (GtkWidget *w, windata_t *vwin)
{
    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);
}

static void open_pkg_sample (GtkWidget *w, windata_t *vwin)
{
    if (viewer_char_count(vwin) > 0) {
	int resp;

	resp = yes_no_dialog(NULL, _("Really replace content with\n"
				     "a selected file?"),
			     vwin->main);
	if (resp != GRETL_YES) {
	    return;
	}
    }

    file_selector(OPEN_SCRIPT, FSEL_DATA_VWIN, vwin);
}

static void toolbar_new_callback (GtkWidget *w, windata_t *vwin)
{
    do_new_script(vwin->role, NULL, NULL);
}

static void window_print_callback (GtkWidget *w, windata_t *vwin)
{
    if (textview_use_highlighting(vwin->role)) {
	int resp = yes_no_cancel_dialog(NULL,
					_("Print with syntax highlighting?"),
					vwin_toplevel(vwin));

	if (resp == GRETL_YES) {
	    sourceview_print(vwin);
	} else if (resp == GRETL_NO) {
	    window_print(NULL, vwin);
	}
    } else {
	window_print(NULL, vwin);
    }
}

static void window_help (GtkWidget *w, windata_t *vwin)
{
    show_gui_help(vwin->role);
}

static void markdown_help (GtkWidget *w, windata_t *vwin)
{
    show_gui_help(MDHELP);
}

static void multi_save_as_callback (GtkWidget *w, windata_t *vwin)
{
    copy_format_dialog(vwin, W_SAVE);
}

static void script_index (GtkWidget *w, windata_t *vwin)
{
    display_files(PS_FILES, NULL);
}

static void toolbar_refresh (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_SERIES) {
	series_view_refresh(w, vwin);
    }
}

static void set_matrix_name (GtkWidget *widget, dialog_t *dlg)
{
    char *vname = (char *) edit_dialog_get_data(dlg);
    GtkWidget *parent =  edit_dialog_get_window(dlg);
    const gchar *s = edit_dialog_get_text(dlg);

    if (s == NULL || gui_validate_varname(s, GRETL_TYPE_MATRIX, parent)) {
	edit_dialog_reset(dlg);
    } else {
	strcpy(vname, s);
	edit_dialog_close(dlg);
    }
}

static void add_matrix_callback (GtkWidget *w, windata_t *vwin)
{
    char mname[VNAMELEN];
    gretl_matrix *m = NULL;
    int err, cancel = 0;

    if (vwin->role != XTAB && vwin->role != COEFFINT &&
	vwin->role != ALAGSEL && vwin->role != VLAGSEL) {
	dummy_call();
	return;
    }

    blocking_edit_dialog(0, _("gretl: save matrix"),
			 _("Enter a name"), NULL,
			 set_matrix_name, mname,
			 VARCLICK_NONE,
			 vwin_toplevel(vwin),
			 &cancel);
    if (!cancel) {
	if (vwin->role == XTAB) {
	    m = xtab_to_matrix(vwin->data);
	} else if (vwin->role == COEFFINT) {
	    m = conf_intervals_matrix(vwin->data);
	} else if (vwin->role == ALAGSEL || vwin->role == VLAGSEL) {
	    m = gretl_matrix_copy(vwin->data);
	}
	if (m == NULL) {
	    nomem();
	} else {
	    err = user_var_add_or_replace(mname,
					  GRETL_TYPE_MATRIX,
					  m);
	    if (err) {
		gretl_matrix_free(m);
		gui_errmsg(err);
	    } else {
		infobox_printf(_("Saved matrix as %s"), mname);
	    }
	}
    }
}

static void add_data_callback (GtkWidget *w, windata_t *vwin)
{
    int oldv = dataset->v;

    if (vwin->role == PCA) {
	add_pca_data(vwin);
    } else if (vwin->role == LEVERAGE) {
	add_leverage_data(vwin);
    } else if (vwin->role == MAHAL) {
	add_mahalanobis_data(vwin);
    } else if (vwin->role == FCAST) {
	add_fcast_data(vwin, M_FCAST);
    } else if (vwin->role == LOESS || vwin->role == NADARWAT) {
	add_nonparam_data(vwin);
    } else if (vwin->role == VIEW_DBNOMICS) {
	add_dbnomics_data(vwin);
    }

    if (dataset->v > oldv) {
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
    print_coeff_intervals(cf, prn);
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
	series_view_format_dialog(vwin);
    }
}

static void split_pane_callback (GtkWidget *w, windata_t *vwin)
{
    GtkWidget *hb = g_object_get_data(G_OBJECT(w), "hpane");
    GtkWidget *vb = g_object_get_data(G_OBJECT(w), "vpane");
    int vertical = 0;

    if (hb != NULL) {
	vb = w;
	vertical = 1;
    } else {
	hb = w;
    }

    /* Note: by "vertical" here we mean that the split runs vertically,
       dividing the pane into left- and right-hand sections; otherwise
       the split runs horizontally. In a "gnuplot commands" window we
       only offer a horizontal split, which means that @vb ("vertical
       button") may be NULL.
    */

    if (g_object_get_data(G_OBJECT(vwin->vbox), "sw") != NULL) {
	/* currently in single-view mode: so split */
	viewer_split_pane(vwin, vertical);
	if (vb != NULL) {
	    gtk_widget_set_sensitive(vb, vertical);
	}
	gtk_widget_set_sensitive(hb, !vertical);
	if (vertical) {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb),
					 GRETL_STOCK_JOIN_V);
	} else {
	    gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb),
					 GRETL_STOCK_JOIN_H);
	}
    } else {
	GtkWidget *paned;

	paned = g_object_get_data(G_OBJECT(vwin->vbox), "paned");

	if (paned != NULL) {
	    /* currently in split-view mode: so rejoin */
	    vertical = GTK_IS_HPANED(paned);
	    viewer_close_pane(vwin);
	    gtk_widget_set_sensitive(hb, TRUE);
	    if (vb != NULL) {
		gtk_widget_set_sensitive(vb, TRUE);
	    }
	    if (vertical) {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(vb),
					     GRETL_STOCK_SPLIT_V);
	    } else {
		gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(hb),
					     GRETL_STOCK_SPLIT_H);
	    }
	}
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
    GtkAdjustment *adj;
    gchar txt[16];
    double x = 1.0 - cf->alpha;
    gboolean defset = FALSE;
    int i;

    if (maybe_raise_dialog()) {
	return;
    }

    dialog = gretl_dialog_new(_("gretl: coefficient confidence intervals"),
			      vwin_toplevel(vwin), GRETL_DLG_BLOCK);

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
    adj = (GtkAdjustment *) gtk_adjustment_new(x, 0.60, 0.99, 0.01, 0, 0);
    tmp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.01, 2);
    g_signal_connect(G_OBJECT(tmp), "value-changed",
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
    cancel_delete_button(hbox, dialog);

    g_object_set_data(G_OBJECT(dialog), "vwin", vwin);
    g_object_set_data(G_OBJECT(dialog), "xptr", &x);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(real_coeffint_set_alpha), dialog);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);
}

static void stickiness_callback (GtkWidget *w, windata_t *vwin)
{
    output_policy_dialog(vwin, vwin, 1);
}

static void do_corr_plot (windata_t *vwin)
{
    VMatrix *corr = vwin->data;
    int err;

    err = plot_corrmat(corr, OPT_NONE);
    gui_graph_handler(err);
}

static void toolbar_plot_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_SERIES) {
	series_view_graph(w, vwin);
    } else if (vwin->role == VIEW_BUNDLE) {
	exec_bundle_special_function(vwin->data, BUNDLE_PLOT,
				     NULL, vwin->main);
    } else if (vwin->role == CORR) {
	do_corr_plot(vwin);
    } else if (vwin->role == VIEW_DBNOMICS) {
	show_dbnomics_data(vwin, 1);
    } else {
	do_nonparam_plot(vwin);
    }
}

static void toolbar_fcast_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == VIEW_BUNDLE) {
	exec_bundle_special_function(vwin->data, BUNDLE_FCAST,
				     NULL, vwin->main);
    }
}

static void dbnomics_show_series (GtkWidget *w, windata_t *vwin)
{
    show_dbnomics_data(vwin, 0);
}

static void editor_prefs_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin->role == CONSOLE) {
	console_prefs_dialog(vwin->main);
    } else {
	preferences_dialog(TAB_EDITOR, NULL, vwin_toplevel(vwin));
    }
}

static void build_pkg_callback (GtkWidget *w, windata_t *vwin)
{
    if (vwin_content_changed(vwin)) {
	int resp;

	resp = yes_no_cancel_dialog("gretl", _("Save changes?"),
				    vwin->main);
	if (resp == GRETL_CANCEL) {
	    return;
	}
	if (resp == GRETL_YES) {
	    vwin_save_callback(NULL, vwin);
	}
    }

    build_package_from_spec_file(vwin);
}

static void enable_markdown_editor (windata_t *vwin,
				    gboolean s)
{
    GtkTextView *view = GTK_TEXT_VIEW(vwin->text);
    GtkToolbar *tbar = GTK_TOOLBAR(vwin->mbar);
    int i, n = gtk_toolbar_get_n_items(tbar);
    GtkToolItem *item;

    /* item 0 is the Save button, which should be sensitive
       just if editable content is changed, and item n-1 is
       the preview toggle button, which should not be made
       insensitive (yet)
    */
    item = gtk_toolbar_get_nth_item(tbar, 0);
    gtk_widget_set_sensitive(GTK_WIDGET(item),
			     vwin->flags & VWIN_CONTENT_CHANGED);
    for (i=1; i<n-1; i++) {
	item = gtk_toolbar_get_nth_item(tbar, i);
	if (GTK_IS_TOOL_BUTTON(item)) {
	    gtk_widget_set_sensitive(GTK_WIDGET(item), s);
	}
    }

    gtk_text_view_set_editable(view, s);
    gtk_text_view_set_cursor_visible(view, s);
    widget_set_int(vwin->text, "preview_on", !s);
}

static void toggle_md_preview (GtkWidget *w, windata_t *vwin)
{
    GtkTextBuffer *tbuf;
    GtkToolItem *item;
    PRN *prn = NULL;
    gulong cid;

    tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
    cid = widget_get_int(tbuf, "changed_id");

    if (widget_get_int(vwin->text, "preview_on")) {
	/* return to editing mode */
	gchar *buf = g_object_get_data(G_OBJECT(vwin->mbar), "raw_text");

	textview_set_text(vwin->text, buf);
	item = g_object_get_data(G_OBJECT(vwin->mbar), "eye_button");
	gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(item), GRETL_STOCK_EYE);
	gtk_tool_item_set_tooltip_text(item, _("Preview as markdown"));
	enable_markdown_editor(vwin, TRUE);
	g_free(buf);
	g_object_set_data(G_OBJECT(vwin->mbar), "raw_text", NULL);
	g_signal_handler_unblock(G_OBJECT(tbuf), cid);
    } else if ((prn = gui_prn_new()) != NULL) {
	/* preview formatted markdown */
	gchar *buf = textview_get_text(vwin->text);
	char *pbuf;

	g_signal_handler_block(G_OBJECT(tbuf), cid);
	g_object_set_data(G_OBJECT(vwin->mbar), "raw_text", buf);
	md_to_gretl(buf, prn);
	pbuf = gretl_print_steal_buffer(prn);
	gretl_viewer_insert_formatted_buffer(vwin, pbuf);
	item = g_object_get_data(G_OBJECT(vwin->mbar), "eye_button");
	gtk_tool_button_set_stock_id(GTK_TOOL_BUTTON(item), GRETL_STOCK_EYE_OFF);
	gtk_tool_item_set_tooltip_text(item, _("Return to editing"));
	enable_markdown_editor(vwin, FALSE);
	gretl_print_destroy(prn);
    }
}

static void set_md_status (GtkWidget *w, windata_t *vwin)
{
    const char *opts[] = {
	N_("plain literal text"),
	N_("markdown")
    };
    int def = widget_get_int(vwin->text, "md");
    int resp;

    resp = radio_dialog(_("Text status"), NULL, opts, 2, def, 0, vwin->main);
    if (resp >= 0 && resp != def) {
	widget_set_int(vwin->text, "md", resp);
    }
}

static int bundle_plot_ok (windata_t *vwin)
{
    gretl_bundle *b = vwin->data;
    gchar *pf = get_bundle_special_function(b, BUNDLE_PLOT);
    int ret = 0;

    if (pf != NULL) {
	const char *s;

	ret = 1;
	g_free(pf);
    }

    return ret;
}

static int bundle_fcast_ok (windata_t *vwin)
{
    gretl_bundle *b = vwin->data;
    gchar *ff = get_bundle_special_function(b, BUNDLE_FCAST);
    int ret = 0;

    if (ff != NULL) {
	ret = 1;
	g_free(ff);
    }

    return ret;
}

static int suppress_hmap (VMatrix *corr)
{
    return corr->dim < 3;
}

static void activate_script_help (GtkWidget *widget, windata_t *vwin)
{
    text_set_cursor(vwin->text, GDK_QUESTION_ARROW);
    set_window_help_active(vwin);
}

static int edit_script_popup_item (GretlToolItem *item)
{
    if (item->icon == NULL) return 0;

    return !strcmp(item->icon, GTK_STOCK_COPY) ||
	!strcmp(item->icon, GTK_STOCK_PASTE) ||
	!strcmp(item->icon, GTK_STOCK_FIND) ||
	!strcmp(item->icon, GTK_STOCK_UNDO) ||
	!strcmp(item->icon, GTK_STOCK_FIND_AND_REPLACE);
}

static void set_plot_icon (GretlToolItem *item, windata_t *vwin)
{
    int r = vwin->role;

    if (r == LOESS || r == NADARWAT || r == VIEW_BUNDLE) {
	int done = 0;

	if (r == VIEW_BUNDLE) {
	    const char *s = gretl_bundle_get_creator(vwin->data);

	    if (s != NULL && !strcmp(s, "regls")) {
		item->icon = GRETL_STOCK_CP;
		done = 1;
	    }
	}
	if (!done) {
	    item->icon = GRETL_STOCK_SCATTER;
	}
    } else if (r == VIEW_SERIES && dataset_is_cross_section(dataset)) {
	item->icon = GRETL_STOCK_BOX;
    } /* else stay with the default, a time series plot icon */
}

static void vwin_cut_callback (GtkWidget *w, windata_t *vwin)
{
    gtk_text_buffer_cut_clipboard(gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text)),
				  gtk_clipboard_get(GDK_NONE),
				  TRUE);
}

static GretlToolItem viewbar_items[] = {
    { N_("New window"), GTK_STOCK_NEW, G_CALLBACK(toolbar_new_callback), NEW_ITEM },
    { N_("Open..."), GTK_STOCK_OPEN, G_CALLBACK(file_open_callback), OPEN_ITEM },
    { N_("Save"), GTK_STOCK_SAVE, G_CALLBACK(vwin_save_callback), SAVE_ITEM },
    { N_("Save as..."), GTK_STOCK_SAVE_AS, G_CALLBACK(save_as_callback), SAVE_AS_ITEM },
    { N_("Open in script editor"), GTK_STOCK_EDIT, G_CALLBACK(copy_to_editor), COPY_SCRIPT_ITEM },
    { N_("Save bundle content..."), GRETL_STOCK_BUNDLE, NULL, BUNDLE_ITEM },
    { N_("Print..."), GTK_STOCK_PRINT, G_CALLBACK(window_print_callback), PRINT_ITEM },
    { N_("Show/hide"), GRETL_STOCK_PIN, G_CALLBACK(session_notes_callback), NOTES_ITEM },
    { N_("Display values"), GTK_STOCK_MEDIA_PLAY, G_CALLBACK(dbnomics_show_series), DBN_ITEM },
    { N_("Run"), GTK_STOCK_EXECUTE, G_CALLBACK(do_run_script), EXEC_ITEM },
    { N_("Build package"), GRETL_STOCK_TOOLS, G_CALLBACK(build_pkg_callback), BUILD_ITEM },
    { N_("Cut"), GTK_STOCK_CUT, G_CALLBACK(vwin_cut_callback), EDIT_ITEM },
    { N_("Copy"), GTK_STOCK_COPY, G_CALLBACK(vwin_copy_callback), COPY_ITEM },
    { N_("Paste"), GTK_STOCK_PASTE, G_CALLBACK(text_paste), EDIT_ITEM },
    { N_("Find..."), GTK_STOCK_FIND, G_CALLBACK(text_find), FIND_ITEM },
    { N_("Replace..."), GTK_STOCK_FIND_AND_REPLACE, G_CALLBACK(text_replace), EDIT_ITEM },
    { N_("Undo"), GTK_STOCK_UNDO, G_CALLBACK(text_undo), EDIT_ITEM },
    { N_("Redo"), GTK_STOCK_REDO, G_CALLBACK(text_redo), EDIT_ITEM },
    { N_("Sort"), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(series_view_toggle_sort), SORT_ITEM },
    { N_("Sort by..."), GTK_STOCK_SORT_ASCENDING, G_CALLBACK(multi_series_view_sort_by), SORT_BY_ITEM },
    { N_("Preferences..."), GTK_STOCK_PREFERENCES, G_CALLBACK(editor_prefs_callback), EDIT_HANSL_ITEM },
    { N_("Auto-indent script"), GTK_STOCK_INDENT, G_CALLBACK(indent_hansl), EDIT_HANSL_ITEM },
    { N_("Send To..."), GRETL_STOCK_MAIL, G_CALLBACK(mail_script_callback), MAIL_ITEM },
    { N_("Scripts index"), GTK_STOCK_INDEX, G_CALLBACK(script_index), INDEX_ITEM },
    { N_("Confidence level..."), GRETL_STOCK_ALPHA, G_CALLBACK(coeffint_set_alpha), ALPHA_ITEM },
    { N_("LaTeX"), GRETL_STOCK_TEX, G_CALLBACK(window_tex_callback), TEX_ITEM },
    { N_("Graph"), GRETL_STOCK_TS, G_CALLBACK(toolbar_plot_callback), PLOT_ITEM },
    { N_("Forecast"), GRETL_STOCK_FCAST, G_CALLBACK(toolbar_fcast_callback), FCAST_ITEM },
    { N_("Heatmap"), GRETL_STOCK_HMAP, G_CALLBACK(toolbar_plot_callback), HMAP_ITEM },
    { N_("Reformat..."), GTK_STOCK_CONVERT, G_CALLBACK(reformat_callback), FORMAT_ITEM },
    { N_("Edit values..."), GTK_STOCK_EDIT, G_CALLBACK(series_view_edit), EDITOR_ITEM },
    { N_("Refresh"), GTK_STOCK_REFRESH, G_CALLBACK(toolbar_refresh), REFRESH_ITEM },
    { N_("Add to dataset..."), GTK_STOCK_ADD, G_CALLBACK(add_data_callback), ADD_DATA_ITEM },
    { N_("Add as matrix..."), GTK_STOCK_ADD, G_CALLBACK(add_matrix_callback), ADD_MATRIX_ITEM },
    { N_("Stickiness..."), GRETL_STOCK_PIN, G_CALLBACK(stickiness_callback), STICKIFY_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_H, G_CALLBACK(split_pane_callback), SPLIT_H_ITEM },
    { N_("Toggle split pane"), GRETL_STOCK_SPLIT_V, G_CALLBACK(split_pane_callback), SPLIT_V_ITEM },
    { N_("Clear"), GTK_STOCK_CLEAR, G_CALLBACK(clear_console), CLEAR_ITEM },
    { N_("Help on command"), GRETL_STOCK_QUERY, G_CALLBACK(activate_script_help), CMD_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(window_help), HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(display_gnuplot_help), GP_HELP_ITEM },
    { N_("Help"), GTK_STOCK_HELP, G_CALLBACK(display_x12a_help), X12A_HELP_ITEM },
    { N_("Digits..."), NULL, G_CALLBACK(NULL), DIGITS_ITEM }
};

static int n_viewbar_items = G_N_ELEMENTS(viewbar_items);

#define exec_ok(r) (vwin_editing_script(r) || \
		    r == VIEW_SCRIPT || \
		    r == VIEW_PKG_SAMPLE || \
	            r == EDIT_PKG_SAMPLE)

#define open_ok(r) (vwin_editing_script(r))

#define new_ok(r) (vwin_editing_script(r))

#define edit_ok(r) (vwin_editing_script(r) || \
		    vwin_editing_buffer(r) || \
                    r == EDIT_PKG_CODE || \
		    r == EDIT_PKG_SAMPLE || \
		    r == EDIT_PKG_HELP || \
		    r == EDIT_PKG_GHLP)

#define save_as_ok(r) (r != EDIT_HEADER && \
	               r != EDIT_NOTES && \
	               r != EDIT_PKG_CODE && \
		       r != EDIT_PKG_SAMPLE && \
		       r != EDIT_PKG_HELP && \
		       r != EDIT_PKG_GHLP && \
		       r != CONSOLE && \
		       r != VIEW_BUNDLE && \
		       r != VIEW_DBNOMICS)

#define help_ok(r) (r == LEVERAGE || \
		    r == COINT2 || \
		    r == HURST || \
		    r == RMPLOT || \
		    r == MAHAL)

#define cmd_help_ok(r) (r == EDIT_HANSL || \
	                r == EDIT_PKG_CODE || \
                        r == EDIT_PKG_SAMPLE || \
			r == VIEW_PKG_SAMPLE || \
			r == VIEW_PKG_CODE || \
			r == VIEW_SCRIPT || \
			r == VIEW_LOG || \
			r == CONSOLE)

/* for a non-editable script: can offer option to copy
   content into an editor window */
#define copy_script_ok(r) (r == VIEW_PKG_SAMPLE || \
			   r == VIEW_PKG_CODE || \
			   r == VIEW_SCRIPT || \
			   r == VIEW_LOG)

#define sort_ok(r) (r == VIEW_SERIES)

#define plot_ok(r) (r == VIEW_SERIES || \
		    r == LOESS || \
		    r == NADARWAT || \
		    r == VIEW_DBNOMICS)

#define add_data_ok(r) (r == PCA || r == LEVERAGE || \
                        r == MAHAL || r == FCAST || \
			r == LOESS || r == NADARWAT || \
			r == VIEW_DBNOMICS)

#define add_matrix_ok(r) (r == XTAB || r == COEFFINT || \
			  r == ALAGSEL || r == VLAGSEL)

#define split_h_ok(r) (r == SCRIPT_OUT || r == FNCALL_OUT || \
		       r == VIEW_LOG || r == VIEW_PKG_CODE || \
		       r == VIEW_BUNDLE || r == X12A || \
		       vwin_editing_script(r))

#define split_v_ok(r) (r == SCRIPT_OUT || r == FNCALL_OUT)

/* Screen out unwanted menu items depending on the context; also
   adjust the callbacks associated with some items based on
   context.
*/

static GCallback tool_item_get_callback (GretlToolItem *item, windata_t *vwin,
					 int latex_ok, int sortby_ok,
					 int format_ok, int save_ok)
{
    static int mail_ok = -1;
    GCallback func = item->func;
    int f = item->flag;
    int r = vwin->role;

    if (mail_ok < 0) {
	mail_ok = curl_does_smtp();
    }

    if (r == EDIT_SPEC) {
	/* This is a "special" that should maybe be regularized:
	   a bit like editing a script, but different...
	*/
	if (f == NEW_ITEM || f == OPEN_ITEM || f == EXEC_ITEM) {
	    return NULL;
	}
    } else if (f == BUILD_ITEM) {
	return NULL;
    }

    /* popup use only */
    if (f == DIGITS_ITEM) {
	return NULL;
    }

    if (use_toolbar_search_box(r) && f == FIND_ITEM) {
	/* using an "inline" search box: skip the
	   "Find" button */
	return NULL;
    }

    if (r == VIEW_DBNOMICS) {
	if (f == PRINT_ITEM || f == FIND_ITEM) {
	    return NULL;
	}
    } else if (f == DBN_ITEM) {
	return NULL;
    }

    if (copy_script_ok(r)) {
	if (f == SAVE_AS_ITEM) {
	    return NULL;
	}
    } else if (f == COPY_SCRIPT_ITEM) {
	return NULL;
    }

    if (f == CLEAR_ITEM && r != CONSOLE) {
	return NULL;
    }

    if (r == EDIT_PKG_SAMPLE && f == OPEN_ITEM) {
	return G_CALLBACK(open_pkg_sample);
    } else if (!edit_ok(r) && f == EDIT_ITEM) {
	return NULL;
    } else if (!open_ok(r) && f == OPEN_ITEM) {
	return NULL;
    } else if (!new_ok(r) && f == NEW_ITEM) {
	return NULL;
    } else if (!exec_ok(r) && f == EXEC_ITEM) {
	return NULL;
    } else if (!cmd_help_ok(r) && f == CMD_HELP_ITEM) {
	return NULL;
    } else if ((!mail_ok || r != EDIT_HANSL) && f == MAIL_ITEM) {
	return NULL;
    } else if (!help_ok(r) && f == HELP_ITEM) {
	return NULL;
    } else if ((!latex_ok || !multiple_formats_ok(vwin)) && f == TEX_ITEM) {
	return NULL;
    } else if (!add_data_ok(r) && f == ADD_DATA_ITEM) {
	return NULL;
    } else if (!add_matrix_ok(r) && f == ADD_MATRIX_ITEM) {
	return NULL;
    } else if (!sort_ok(r) && f == SORT_ITEM) {
	return NULL;
    } else if (!sortby_ok && f == SORT_BY_ITEM) {
	return NULL;
    } else if (!plot_ok(r) && f == PLOT_ITEM) {
	if (r == VIEW_BUNDLE && bundle_plot_ok(vwin)) {
	    ; /* alright then */
	} else {
	    return NULL;
	}
    } else if (f == FCAST_ITEM) {
	if (r == VIEW_BUNDLE && bundle_fcast_ok(vwin)) {
	    ; /* alright then */
	} else {
	    return NULL;
	}
    } else if (!split_h_ok(r) && f == SPLIT_H_ITEM) {
	return NULL;
    } else if (!split_v_ok(r) && f == SPLIT_V_ITEM) {
	return NULL;
    } else if (!format_ok && f == FORMAT_ITEM) {
	return NULL;
    } else if (r != VIEW_SERIES && f == EDITOR_ITEM) {
	return NULL;
    } else if (r == CONSOLE && func == G_CALLBACK(editor_prefs_callback)) {
	; /* alright then */
    } else if (r != EDIT_HANSL && r != EDIT_PKG_CODE &&
	       r != EDIT_PKG_SAMPLE && f == EDIT_HANSL_ITEM) {
	return NULL;
    } else if (r != VIEW_SCRIPT && f == INDEX_ITEM) {
	return NULL;
    } else if (r != SCRIPT_OUT && f == STICKIFY_ITEM) {
	return NULL;
    } else if (r != COEFFINT && f == ALPHA_ITEM) {
	return NULL;
    } else if (r != VIEW_SERIES && f == REFRESH_ITEM) {
	return NULL;
    } else if (r != EDIT_GP && f == GP_HELP_ITEM) {
	return NULL;
    } else if (r != EDIT_X12A && f == X12A_HELP_ITEM) {
	return NULL;
    } else if (f == SAVE_ITEM && !save_ok) {
	return NULL;
    } else if (r != EDIT_NOTES && f == NOTES_ITEM) {
	return NULL;
    } else if (r != VIEW_BUNDLE && r != VIEW_DBNOMICS && f == BUNDLE_ITEM) {
	return NULL;
    } else if (f == HMAP_ITEM) {
	if (r != CORR) {
	    return NULL;
	} else if (suppress_hmap(vwin->data)) {
	    return NULL;
	}
    } else if (f == SAVE_AS_ITEM) {
	if (!save_as_ok(r) || (vwin->flags & VWIN_NO_SAVE)) {
	    return NULL;
	} else if (multiple_formats_ok(vwin) ||
		   (vwin->flags & VWIN_MULTI_SERIES)) {
	    func = G_CALLBACK(multi_save_as_callback);
	}
    }

    return func;
}

static void fcast_save_call (GtkAction *action, gpointer p)
{
    const char *s = gtk_action_get_name(action);
    ModelDataIndex idx;

    idx = (strcmp(s, "SaveSderr") == 0)? M_FCSE : M_FCAST;

    add_fcast_data((windata_t *) p, idx);
}

static GtkWidget *make_fcast_save_menu (windata_t *vwin)
{
    GtkWidget *menu = gtk_menu_new();
    GtkAction *action;
    GtkWidget *item;
    int i;

    for (i=0; i<2; i++) {
	action = gtk_action_new(i == 0 ? "SaveFcast" : "SaveSderr",
				i == 0 ? _("_Save forecast...") :
				_("Save standard errors..."),
				NULL, NULL);
	g_signal_connect(G_OBJECT(action), "activate",
			 G_CALLBACK(fcast_save_call), vwin);
	item = gtk_action_create_menu_item(action);
	gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
    }

    return menu;
}

static GtkWidget *tool_item_get_menu (GretlToolItem *item,
				      windata_t *vwin,
				      int *insensitive)
{
    GtkWidget *menu = NULL;

    if (vwin->role == VIEW_BUNDLE) {
	if (item->flag == BUNDLE_ITEM) {
	    menu = make_bundle_content_menu(vwin);
	} else if (item->flag == PLOT_ITEM) {
	    menu = make_bundle_plot_menu(vwin, insensitive);
	} else if (item->flag == SAVE_ITEM) {
	    menu = make_bundle_save_menu(vwin);
	    if (menu != NULL) {
		item->tip = N_("Save...");
	    }
	}
    } else if (vwin->role == VIEW_DBNOMICS) {
	if (item->flag == BUNDLE_ITEM) {
	    menu = make_bundle_content_menu(vwin);
	} else if (item->flag == SAVE_ITEM) {
	    menu = make_bundle_save_menu(vwin);
	    if (menu != NULL) {
		item->tip = N_("Save...");
	    }
	}
    } else if (vwin->role == FCAST && item->flag == ADD_DATA_ITEM) {
	FITRESID *fr = vwin->data;

	if (fr->sderr != NULL) {
	    menu = make_fcast_save_menu(vwin);
	}
    }

    if (menu != NULL) {
	/* don't leak: record pointer to menu so it can
	   be destroyed when the window is closed */
	vwin_record_toolbar_popup(vwin, menu);
    }

    return menu;
}

static void gretl_toolbar_flat (GtkWidget *w)
{
    static int style_done;

    gtk_widget_set_name(w, "gretl_toolbar");

    if (!style_done) {
	gtk_rc_parse_string("style \"gretl-tb-style\"\n{\n"
			    "  GtkToolbar::shadow-type = GTK_SHADOW_NONE\n"
			    "}\n"
			    "widget \"*.gretl_toolbar\" style \"gretl-tb-style\"");
	style_done = 1;
    }
}

GtkWidget *gretl_toolbar_new (GtkWidget *sibling)
{
    GtkWidget *tb = gtk_toolbar_new();

    gtk_toolbar_set_icon_size(GTK_TOOLBAR(tb), toolbar_icon_size);
    gtk_toolbar_set_style(GTK_TOOLBAR(tb), GTK_TOOLBAR_ICONS);
    gtk_toolbar_set_show_arrow(GTK_TOOLBAR(tb), FALSE);

    if (sibling == NULL) {
	/* if we're not alongside a menu bar ("sibling"),
	   show the toolbar without a shadow
	*/
	gretl_toolbar_flat(tb);
    }

    return tb;
}

void gretl_tooltips_add (GtkWidget *w, const gchar *str)
{
    gtk_widget_set_tooltip_text(w, str);
}

static GtkToolItem *gretl_menu_button (const char *icon,
				       const char *tip,
				       GtkWidget **pw)
{
    GtkWidget *img, *button = gtk_button_new();
    GtkToolItem *item = gtk_tool_item_new();

    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _(tip));
    gtk_button_set_relief(GTK_BUTTON(button), GTK_RELIEF_NONE);
    img = gtk_image_new_from_stock(icon, toolbar_icon_size /* GTK_ICON_SIZE_MENU */);
    gtk_container_add(GTK_CONTAINER(button), img);
    gtk_container_add(GTK_CONTAINER(item), button);
    *pw = button;

    return item;
}

static void gretl_tool_item_set_tip (GtkToolItem *item,
				     GretlToolItem *tool)
{
    const char *accel = NULL;

    if (tool->flag == EXEC_ITEM) {
	accel = "Ctrl+R";
    } else if (tool->flag == COPY_ITEM) {
	accel = "Ctrl+C";
    } else if (tool->flag == SAVE_ITEM) {
	accel = "Ctrl+S";
    } else if (tool->flag == FIND_ITEM) {
	accel = "Ctrl+F";
    } else if (!strcmp(tool->icon, GTK_STOCK_FIND_AND_REPLACE)) {
	accel = "Ctrl+H";
    }

    if (accel != NULL) {
	gchar *s = g_strdup_printf("%s (%s)", _(tool->tip), accel);

	gtk_tool_item_set_tooltip_text(item, s);
	g_free(s);
    } else {
	gtk_tool_item_set_tooltip_text(item, _(tool->tip));
    }
}

GtkWidget *gretl_toolbar_insert (GtkWidget *tbar,
				 GretlToolItem *tool,
				 GCallback func,
				 gpointer data,
				 gint pos)
{
    GtkToolItem *item;

    item = gtk_tool_button_new_from_stock(tool->icon);
    gretl_tool_item_set_tip(item, tool);
    g_signal_connect(G_OBJECT(item), "clicked", func, data);
    gtk_widget_set_size_request(GTK_WIDGET(item), 30, -1);
    gtk_toolbar_insert(GTK_TOOLBAR(tbar), item, pos);

    return GTK_WIDGET(item);
}

static void button_menu_pos (GtkMenu *menu,
			     gint *x,
			     gint *y,
			     gboolean *push_in,
			     gpointer data)
{
    GtkWidget *button = data;
    gint wx, wy, tx, ty;

    gdk_window_get_origin(gtk_widget_get_window(button), &wx, &wy);
    gtk_widget_translate_coordinates(button, gtk_widget_get_toplevel(button),
				     0, 0, &tx, &ty);
    *x = wx + tx;
    *y = wy + ty + 26;
    *push_in = TRUE;
}

static void tool_item_popup (GtkWidget *button, GdkEvent *event,
			     GtkWidget *menu)
{
    gtk_menu_popup(GTK_MENU(menu), NULL, NULL,
		   button_menu_pos, button,
		   event->button.button, event->button.time);
}

/* right-click and middle-click actions for Run (exec) button */

static gint exec_press (GtkWidget *w, GdkEventButton *eb, windata_t *vwin)
{
    if (eb->button == 3) {
	run_script_silent(NULL, vwin);
	return TRUE;
    } else {
	return FALSE;
    }
}

GtkWidget *vwin_toolbar_insert (GretlToolItem *tool,
				GCallback func,
				GtkWidget *menu,
				windata_t *vwin,
				gint pos)
{
    GtkToolItem *item;

    if (menu != NULL) {
	/* make and insert a button that pops down a menu */
	GtkWidget *button;

	item = gretl_menu_button(tool->icon, tool->tip, &button);
	g_signal_connect(G_OBJECT(button), "button-press-event",
			 G_CALLBACK(tool_item_popup), menu);
    } else {
	/* make and insert a regular callback button */
	item = gtk_tool_button_new_from_stock(tool->icon);
	g_signal_connect(G_OBJECT(item), "clicked", func, vwin);
	if (tool->flag == NEW_ITEM && window_is_tab(vwin)) {
	    gtk_widget_set_tooltip_text(GTK_WIDGET(item), _("New tab"));
	} else {
	    gretl_tool_item_set_tip(item, tool);
	    if (tool->flag == EXEC_ITEM) {
		g_signal_connect(G_OBJECT(item), "button-press-event",
				 G_CALLBACK(exec_press), vwin);
	    }
	}
    }

    gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), item, pos);

    return GTK_WIDGET(item);
}

static void viewbar_add_items (windata_t *vwin, ViewbarFlags flags)
{
    int sortby_ok = has_sortable_data(vwin);
    int format_ok = can_format_data(vwin);
    int latex_ok = latex_is_ok();
    int save_ok = (flags & VIEWBAR_EDITABLE);
    GtkWidget *hpane = NULL, *vpane = NULL;
    GtkWidget *button;
    GtkWidget *menu;
    GretlToolItem *item;
    GCallback func;
    int insensitive;
    int i;

    for (i=0; i<n_viewbar_items; i++) {
	item = &viewbar_items[i];
	insensitive = 0;
	func = NULL;
	menu = NULL;

	/* Is there anything to hook up, in context? We
	   try first for a menu to attach to the toolbar
	   button; failing that we test for a "direct"
	   callback function.
	*/
	menu = tool_item_get_menu(item, vwin, &insensitive);
	if (menu == NULL && item->func != NULL) {
	    func = tool_item_get_callback(item, vwin, latex_ok, sortby_ok,
					  format_ok, save_ok);
	}
	if (func == NULL && menu == NULL) {
	    /* nothing to hook up */
	    continue;
	}

	if (item->flag == PLOT_ITEM) {
	    set_plot_icon(item, vwin);
	}

	button = vwin_toolbar_insert(item, func, menu, vwin, -1);
	if (insensitive) {
	    gtk_widget_set_sensitive(button, FALSE);
	}

	if (func == (GCallback) split_pane_callback) {
	    if (hpane == NULL) {
		hpane = button;
	    } else {
		vpane = button;
	    }
	}

	if (item->flag == SAVE_ITEM) {
	    if (vwin->role != CONSOLE &&
		vwin->role != VIEW_BUNDLE &&
		vwin->role != VIEW_DBNOMICS) {
		/* nothing to save just yet */
		g_object_set_data(G_OBJECT(vwin->mbar), "save_button", button);
		gtk_widget_set_sensitive(button, FALSE);
	    }
	} else if (item->flag == SAVE_AS_ITEM) {
	    g_object_set_data(G_OBJECT(vwin->mbar), "save_as_button", button);
	    if (strstr(vwin->fname, "script_tmp")) {
		gtk_widget_set_sensitive(button, FALSE);
	    }
	}
    }

    if (hpane != NULL) {
	g_object_set_data(G_OBJECT(hpane), "vpane", vpane);
    }
    if (vpane != NULL) {
	g_object_set_data(G_OBJECT(vpane), "hpane", hpane);
    }
}

static void destroy_md_text (GtkWidget *w, GtkWidget *mbar)
{
    gchar *buf = g_object_get_data(G_OBJECT(mbar), "raw_text");

    if (buf != NULL) {
	g_free(buf);
    }
}

static void add_markdown_items (windata_t *vwin)
{
    static GretlToolItem markdown_items[] = {
	{ N_("Help"), GTK_STOCK_HELP, G_CALLBACK(markdown_help), 0 },
	{ N_("Set markdown status"), GRETL_STOCK_MD, G_CALLBACK(set_md_status), 0 },
	{ N_("Preview as markdown"), GRETL_STOCK_EYE, G_CALLBACK(toggle_md_preview), 0 }
    };
    GretlToolItem *tool;
    GtkToolItem *item;
    int i;

    item = gtk_separator_tool_item_new();
    gtk_separator_tool_item_set_draw(GTK_SEPARATOR_TOOL_ITEM(item), TRUE);
    gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), item, -1);

    for (i=0; i<3; i++) {
	tool = &markdown_items[i];
	item = gtk_tool_button_new_from_stock(tool->icon);
	g_signal_connect(G_OBJECT(item), "clicked", tool->func, vwin);
	gretl_tool_item_set_tip(item, tool);
	gtk_toolbar_insert(GTK_TOOLBAR(vwin->mbar), item, -1);
    }

    g_object_set_data(G_OBJECT(vwin->mbar), "eye_button", item);
    /* try to ensure no memory leakage */
    g_signal_connect(G_OBJECT(vwin->main), "destroy",
		     G_CALLBACK(destroy_md_text), vwin->mbar);
}

void vwin_add_viewbar (windata_t *vwin, ViewbarFlags flags)
{
    if ((flags & VIEWBAR_HAS_TEXT) || vwin->role == SCRIPT_OUT) {
	g_object_set_data(G_OBJECT(vwin->main), "text_out",
			  GINT_TO_POINTER(1));
    }

    vwin->mbar = gretl_toolbar_new(NULL);
    viewbar_add_items(vwin, flags);
    if (vwin->role == EDIT_PKG_HELP ||
	vwin->role == EDIT_PKG_GHLP) {
	add_markdown_items(vwin);
    }
    vwin_pack_toolbar(vwin);
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
	func = G_CALLBACK(NULL);
	if (item->flag == SPLIT_H_ITEM || item->flag == SPLIT_V_ITEM) {
	    continue;
	} else if (item->flag == DIGITS_ITEM &&
		   (vwin->role == VIEW_MODEL || vwin->role == SUMMARY)) {
	    func = G_CALLBACK(display_digits_callback);
	} else if (vwin->role == EDIT_HANSL) {
	    /* the script editor popup may have some special stuff
	       added: don't clutter it up */
	    if (edit_script_popup_item(item)) {
		func = item->func;
	    } else {
		func = NULL;
	    }
	} else {
	    func = tool_item_get_callback(item, vwin, 0, 0, 0, 0);
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

    if (vwin->role != EDIT_HANSL) {
	if (window_is_undockable(vwin)) {
	    add_undock_popup_item(pmenu, vwin);
	} else if (window_is_dockable(vwin)) {
	    add_dock_popup_item(pmenu, vwin);
	}
    }

    return pmenu;
}

/* callbacks for main-window toolbar icons */

static void tbar_calc (void)
{
#ifdef G_OS_WIN32
    win32_run_async(calculator, NULL);
#else
    gretl_fork("calculator", NULL, NULL);
#endif
}

static void tbar_open_data (void)
{
    display_files(TEXTBOOK_DATA, NULL);
}

static void tbar_command_ref (void)
{
    display_text_help(NULL);
}

static void tbar_xy_graph (void)
{
    if (data_status) {
	if (dataset->v == 2) {
	    do_graph_var(mdata->active_var);
	} else if (mdata_selection_count() == 2) {
	    plot_from_selection(GR_XY);
	} else {
	    selection_dialog(GR_XY, _("gretl: define graph"),
			     NULL, do_graph_from_selector);
	}
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_model (void)
{
    if (data_status) {
	selection_dialog(OLS, _("gretl: specify model"), NULL, do_model);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

static void tbar_new_script (void)
{
    do_new_script(EDIT_HANSL, NULL, NULL);
}

static void tbar_show_funcs (GtkWidget *w, gpointer p)
{
    display_files(FUNC_FILES, NULL);
}

static void maybe_find_series (GtkWidget *w, gpointer p)
{
    if (data_status) {
	listbox_find(w, p);
    } else {
	warnbox(_("Please open a data file first"));
    }
}

/* end toolbar icon callbacks */

static GretlToolItem mainbar_items[] = {
    { N_("launch calculator"),  GRETL_STOCK_CALC,    G_CALLBACK(tbar_calc), 0 },
    { N_("new script"),         GTK_STOCK_EDIT,      G_CALLBACK(tbar_new_script), 0 },
    { N_("open gretl console"), GRETL_STOCK_CONSOLE, G_CALLBACK(gretl_console), 1 },
    { N_("session icon view"),  GRETL_STOCK_ICONS,   G_CALLBACK(view_session), 0 },
    { N_("function packages"),  GRETL_STOCK_FUNC,    G_CALLBACK(tbar_show_funcs), 0 },
    { N_("command reference"),  GTK_STOCK_HELP,      G_CALLBACK(tbar_command_ref), 0 },
    { N_("find series"),        GTK_STOCK_FIND,      G_CALLBACK(maybe_find_series), 0 },
    { N_("X-Y graph"),          GRETL_STOCK_SCATTER, G_CALLBACK(tbar_xy_graph), 0 },
    { N_("OLS model"),          GRETL_STOCK_MODEL,   G_CALLBACK(tbar_model), 0 },
    { N_("databases"),          GRETL_STOCK_DB,      G_CALLBACK(show_native_dbs), 0 },
    { N_("open dataset"),       GTK_STOCK_OPEN,      G_CALLBACK(tbar_open_data), 0 },
};

void add_mainwin_toolbar (GtkWidget *vbox)
{
    GretlToolItem *item;
    GtkWidget *hbox;
    int i, n = G_N_ELEMENTS(mainbar_items);

    mdata->mbar = gretl_toolbar_new(NULL);

    for (i=0; i<n; i++) {
	item = &mainbar_items[i];
	if (swallow && item->flag) {
	    continue;
	}
	gretl_toolbar_insert(mdata->mbar, item, item->func, mdata, -1);
    }

    hbox = gtk_hbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(hbox), mdata->mbar, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

/* Add a temporary menubar for use in a script output window, while
   we're waiting for the output. If the output window is being reused
   this is a bit more complicated; we have to "hide" the regular
   menubar before inserting the temporary one.
 */

void vwin_add_tmpbar (windata_t *vwin)
{
    GretlToolItem stop_item = {
	N_("Stop"),
	GTK_STOCK_STOP,
	G_CALLBACK(do_stop_script),
	0
    };
    GtkWidget *hbox, *tmp;

    hbox = g_object_get_data(G_OBJECT(vwin->main), "top-hbox");

    if (hbox != NULL) {
	/* We're replacing a "real" menubar temporarily: ref. the
	   widgets in @hbox before removing them so we can put
	   them back later.
	*/
	GtkWidget *winlist = g_object_get_data(G_OBJECT(hbox), "winlist");

	g_object_ref(G_OBJECT(vwin->mbar));
	gtk_container_remove(GTK_CONTAINER(hbox), vwin->mbar);
	if (vwin->finder != NULL) {
	    g_object_ref(G_OBJECT(vwin->finder));
	    gtk_container_remove(GTK_CONTAINER(hbox), vwin->finder);
	}
	if (winlist != NULL) {
	    g_object_ref(G_OBJECT(winlist));
	    gtk_container_remove(GTK_CONTAINER(hbox), winlist);
	}
    } else {
	/* starting from scratch */
	hbox = gtk_hbox_new(FALSE, 0);
	g_object_set_data(G_OBJECT(vwin->main), "top-hbox", hbox);
	gtk_box_pack_start(GTK_BOX(vwin->vbox), hbox, FALSE, FALSE, 0);
    }

    tmp = gretl_toolbar_new(NULL);
    gretl_toolbar_insert(tmp, &stop_item, stop_item.func, NULL, 0);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);

    start_wait_for_output(vwin, hbox);
    gtk_widget_show_all(hbox);
}
