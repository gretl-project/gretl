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

/* gpt_dialog.c for gretl -- GUI gnuplot controller dialog */

#include "gretl.h"
#include "plotspec.h"
#include "gpt_control.h"
#include "session.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "calculator.h"

#include "../pixmaps/mouse.xpm"

struct gpt_titles_t {
    char *description; /* How the field will show up in the options dialog */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
};

struct gpt_range_t {
    gint ID;
    GtkWidget *isauto;
    GtkWidget *min;
    GtkWidget *max;
    GtkWidget *lbase;
};

static GtkWidget **linetitle;
static GtkWidget **stylecombo;
static GtkWidget **yaxiscombo;
static GtkWidget **linescale;
static GtkWidget **linewidth;

static GtkWidget *labeltext[MAX_PLOT_LABELS];
static GtkWidget *labeljust[MAX_PLOT_LABELS];
static GtkWidget *labelpos[MAX_PLOT_LABELS];

static GtkWidget *gpt_control;
static GtkWidget *keycombo;
static GtkWidget *fitcombo;
static GtkWidget *termcombo;
static GtkWidget *border_check;
static GtkWidget *markers_check;
static GtkWidget *y2_check;
static GtkWidget *ttfcombo;
static GtkWidget *ttfspin;

#define MAX_AXES 3

struct gpt_range_t axis_range[MAX_AXES];

#define NTITLES 4
#define PLOT_LABEL_POS_LEN 32

struct gpt_titles_t gpt_titles[] = {
    { N_("Title of plot"),  0, NULL },
    { N_("Title for axis"), 1, NULL },
    { N_("Title for axis"), 2, NULL },
    { N_("Title for axis"), 3, NULL },
};

const gchar *fittype_strings[] = {
    N_("none"),
    N_("linear: y = a + b*x"),
    N_("quadratic: y = a + b*x + c*x^2"),
    N_("inverse: y = a + b*(1/x)"),
    N_("loess (locally weighted fit)"),
    NULL
};

static const char *get_font_filename (const char *showname);

static void close_plot_controller (GtkWidget *widget, gpointer data) 
{
    GPT_SPEC *spec = (GPT_SPEC *) data;
    png_plot *plot = (png_plot *) spec->ptr;

    gpt_control = NULL;

    if (plot != NULL) { 
	/* PNG plot window is open */
#ifdef G_OS_WIN32 /* z-order gets messed up */
	GtkWidget *shell = plot_get_shell(plot);

	gdk_window_raise(shell->window);
#endif
	plot_remove_controller(plot);
    } else {
	plotspec_destroy(spec); 
    }
} 

static void flip_manual_range (GtkWidget *widget, gpointer data)
{
    gint i = GPOINTER_TO_INT(data);
    gboolean s = GTK_TOGGLE_BUTTON(axis_range[i].isauto)->active;

    gtk_widget_set_sensitive(GTK_WIDGET(axis_range[i].min), !s);
    gtk_widget_set_sensitive(GTK_WIDGET(axis_range[i].max), !s);
}

static void disable_lbase (GtkWidget *b, GtkWidget *entry)
{
    gboolean s = GTK_TOGGLE_BUTTON(b)->active;

    gtk_widget_set_sensitive(entry, !s);
}

static void enable_lbase (GtkWidget *b, GtkWidget *entry)
{
    gboolean s = GTK_TOGGLE_BUTTON(b)->active;

    gtk_widget_set_sensitive(entry, s);
}

/* Take text from a gtkentry and write to gnuplot spec string */

static void entry_to_gp_string (GtkWidget *w, char *targ, size_t n)
{
    const gchar *str;

    *targ = '\0';

    g_return_if_fail(GTK_IS_ENTRY(w));

    str = gtk_entry_get_text(GTK_ENTRY(w));
    if (str != NULL || *str != '\0') {
	strncat(targ, str, n-1);
    }
}

static FitType fit_type_from_string (const char *s)
{
    FitType f = PLOT_FIT_NONE;
    int i;

    if (s != NULL && *s != '\0') {
	for (i=0; fittype_strings[i] != NULL; i++) {
	    if (!strcmp(s, _(fittype_strings[i]))) {
		f = i;
		break;
	    }
	}
    }

    return f;
}

static void fittype_from_combo (GtkComboBox *box, GPT_SPEC *spec)
{
    FitType f = PLOT_FIT_NONE;
    gchar *s = gtk_combo_box_get_active_text(box);

    f = fit_type_from_string(s);

    if (f == PLOT_FIT_OLS || f == PLOT_FIT_QUADRATIC || 
	f == PLOT_FIT_INVERSE || f == PLOT_FIT_LOESS) {
	plotspec_add_fit(spec, f);
	spec->flags &= ~GPT_FIT_HIDDEN;
    } else if (f == PLOT_FIT_NONE) {
	if (spec->n_lines == 2) {
	    spec->flags |= GPT_FIT_HIDDEN;
	}
	spec->fit = f;
    }

    g_free(s);
}

static gboolean fit_type_changed (GtkComboBox *box, GPT_SPEC *spec)
{
    const char *s1 = spec->yvarname;
    const char *s2 = spec->xvarname;
    gchar *s, *title = NULL;
    FitType f;

    if (*s1 == '\0' || *s2 == '\0') {
	return FALSE;
    }

    s = gtk_combo_box_get_active_text(box);
    f = fit_type_from_string(s);
    g_free(s);

    if (f == PLOT_FIT_OLS) {
	title = g_strdup_printf(_("%s versus %s (with least squares fit)"),
		s1, s2);
    } else if (f == PLOT_FIT_QUADRATIC) {
	title = g_strdup_printf(_("%s versus %s (with quadratic fit)"),
		s1, s2);
    } else if (f == PLOT_FIT_INVERSE) {
	title = g_strdup_printf(_("%s versus %s (with inverse fit)"),
		s1, s2);
    } else if (f == PLOT_FIT_LOESS) {
	title = g_strdup_printf(_("%s versus %s (with loess fit)"),
		s1, s2);
    }

    gtk_entry_set_text(GTK_ENTRY(gpt_titles[0].widget), title);
    g_free(title);
    
    return FALSE;
}

/* take a double (which might be NA) and format it for
   a gtkentry widget */

static void double_to_gp_entry (double x, GtkWidget *w)
{
    if (w != NULL && GTK_IS_ENTRY(w)) {
	gchar *numstr;

	if (na(x)) {
	    numstr = g_strdup("*");
	} else {
	    numstr = g_strdup_printf("%g", x);
	}
	gtk_entry_set_text(GTK_ENTRY(w), numstr);
	g_free(numstr);
    }
}

/* read a double from a gtkentry, with error checking */

static double entry_to_gp_double (GtkWidget *w)
{
    double ret = NADBL;

    if (w != NULL && GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));

#ifdef ENABLE_NLS
	if (s != NULL && *s != '\0') {
	    gchar *tmp = g_strdup(s);

	    charsub(tmp, ',', '.');
	    gretl_push_c_numeric_locale();	    
	    if (check_atof(tmp)) {
		errbox(gretl_errmsg_get());
	    } else {
		ret = atof(tmp);
	    }
	    gretl_pop_c_numeric_locale();
	    g_free(tmp);
	}
#else
	if (s != NULL && *s != '\0') {
	    if (check_atof(s)) {
		errbox(gretl_errmsg_get());
	    } else {
		ret = atof(s);
	    }
	}
#endif
    }

    return ret;
}

#define gp_string_to_entry(w,s) do { \
	gtk_entry_set_text(GTK_ENTRY(w),s); } while (0)

static int
get_label_pos_from_entry (GtkWidget *w, double *pos)
{
    int err = 0;

    if (GTK_IS_ENTRY(w)) {
	const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
	int chk;
    
	chk = sscanf(s, "%lf %lf", &pos[0], &pos[1]);
	if (chk != 2) {
	    errbox(_("Invalid label position, must be X Y"));
	    gtk_editable_select_region(GTK_EDITABLE(w), 0, strlen(s));
	    pos[0] = pos[1] = NADBL;
	    err = 1;
	} 
    } else {
	err = 1;
    }

    return err;
}

static int validate_range (double *r)
{
    int err = 0;

    if (na(r[0]) || na(r[1])) {
	r[0] = r[1] = NADBL;
	err = 1;
    } else if (r[1] <= r[0]) {
	r[0] = r[1] = NADBL;
	err = 1;
    }

    return err;
}

static int 
set_logscale_from_entry (GPT_SPEC *spec, int i, GtkWidget *entry)
{
    double base;
    int err = 0;

    base = entry_get_numeric_value(entry, C_POS_DBL);
    if (na(base)) {
	err = 1;
    } else if (!na(base) && base < 1.1) {
	err = 1;
	errbox("bad base");
    } else {
	spec->logbase[i] = base;
    }

    return err;
}

enum {
    ERRORBARS = 1,
    FILLEDCURVE
};

static void apply_gpt_changes (GtkWidget *widget, GPT_SPEC *spec) 
{
    gchar *s;
    int supress_y2 = 0;
    int i, k, err = 0;

    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].widget != NULL) {
	    entry_to_gp_string(gpt_titles[i].widget, spec->titles[i], 
			       sizeof spec->titles[0]);
	}
    }

    s = gtk_combo_box_get_active_text(GTK_COMBO_BOX(keycombo));
    strcpy(spec->keyspec, s);
    g_free(s);

    spec->flags &= ~GPT_Y2AXIS;

    if (y2_check != NULL) {
	if (GTK_TOGGLE_BUTTON(y2_check)->active) {
	    supress_y2 = 1;
	} 
    } 

    for (i=0; i<spec->n_lines; i++) {
	spec->lines[i].yaxis = 1;
	if (!supress_y2 && yaxiscombo[i] != NULL) {
	    fprintf(stderr, "line %d, yaxiscombo[i]=%p\n",
		    i, (void *) yaxiscombo[i]);
	    gchar *str = 
		gtk_combo_box_get_active_text(GTK_COMBO_BOX(yaxiscombo[i]));

	    if (!strcmp(str, "right")) {
		spec->lines[i].yaxis = 2;	
	    }
	    if (spec->lines[i].yaxis == 2) {
		spec->flags |= GPT_Y2AXIS;
	    }
	    g_free(str);
	}
    }

    if (spec->code == PLOT_REGULAR) {
	k = (spec->flags & GPT_Y2AXIS)? 3 : 2;
	for (i=0; i<k; i++) {
	    if (axis_range[i].isauto != NULL) {
		if (GTK_TOGGLE_BUTTON(axis_range[i].isauto)->active) {
		    spec->range[i][0] = NADBL;
		    spec->range[i][1] = NADBL;
		} else {
		    spec->range[i][0] = entry_to_gp_double(axis_range[i].min);
		    spec->range[i][1] = entry_to_gp_double(axis_range[i].max);
		    err = validate_range(spec->range[i]);
		}
	    }
	    if (axis_range[i].lbase != NULL) {
		if (GTK_WIDGET_SENSITIVE(axis_range[i].lbase)) {
		    err = set_logscale_from_entry(spec, i, axis_range[i].lbase);
		} else {
		    spec->logbase[i] = 0.0;
		}
	    }
	}
    }

    if (!err) {   
	for (i=0; i<spec->n_lines; i++) {
	    if (stylecombo[i] != NULL) {
		int oldalt = 0;

		if (!strncmp(spec->lines[i].style, "filled", 6)) {
		    oldalt = FILLEDCURVE;
		} else if (!strncmp(spec->lines[i].style, "error", 5)) {
		    oldalt = ERRORBARS;
		}
		entry_to_gp_string(GTK_COMBO(stylecombo[i])->entry, 
				   spec->lines[i].style, 
				   sizeof spec->lines[0].style);
		if (oldalt == FILLEDCURVE &&
		    !strncmp(spec->lines[i].style, "error", 5)) {
		    spec->flags |= GPT_ERR_SWITCH;
		    spec->flags &= ~GPT_FILL_SWITCH;
		} else if (oldalt == ERRORBARS &&
			   !strncmp(spec->lines[i].style, "filled", 6)) {
		    spec->flags |= GPT_FILL_SWITCH;
		    spec->flags &= ~GPT_ERR_SWITCH;
		}
	    }
	    if (linetitle[i] != NULL) {
		entry_to_gp_string(linetitle[i], 
				   spec->lines[i].title, 
				   sizeof spec->lines[0].title);
	    }
	    if (linescale[i] != NULL) {
		entry_to_gp_string(linescale[i], 
				   spec->lines[i].scale, 
				   sizeof spec->lines[0].scale);
	    }
	    if (linewidth[i] != NULL) {
		spec->lines[i].width = 
		    gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(linewidth[i]));
	    }
	}
    }

    for (i=0; i<MAX_PLOT_LABELS && !err; i++) {
	entry_to_gp_string(labeltext[i], spec->labels[i].text, 
			   sizeof spec->labels[0].text);
	if (string_is_blank(spec->labels[i].text)) {
	    continue;
	}
	err = get_label_pos_from_entry(labelpos[i], spec->labels[i].pos);
	if (err) {
	    break;
	}
	spec->labels[i].just = 
	    gtk_option_menu_get_history(GTK_OPTION_MENU(labeljust[i]));
    } 

    if (!err && border_check != NULL) {
	if (GTK_TOGGLE_BUTTON(border_check)->active) {
	    spec->flags &= ~GPT_MINIMAL_BORDER;
	} else {
	    spec->flags |= GPT_MINIMAL_BORDER;
	}
    } 

    if (!err && markers_check != NULL) {
	if (GTK_TOGGLE_BUTTON(markers_check)->active) {
	    free(spec->labeled);
	    spec->labeled = NULL;
	    spec->flags |= GPT_ALL_MARKERS;
	} else {
	    spec->flags &= ~GPT_ALL_MARKERS;
	}
    }

    if (!err && ttfcombo != NULL && ttfspin != NULL) {
	const gchar *tmp = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry));
	int ptsize = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(ttfspin));

	if (tmp != NULL && *tmp != '\0') {
	    const char *fname = get_font_filename(tmp);
	    char pngfont[128];

	    if (fname != NULL && ptsize > 5 && ptsize < 25) {
		sprintf(pngfont, "%s %d", fname, ptsize);
	    } else {
		*pngfont = '\0';
	    }
	    set_gretl_png_font(pngfont, &paths);
	}
    }

    if (!err && fitcombo != NULL) {
	fittype_from_combo(GTK_COMBO_BOX(fitcombo), spec);
    }

    if (!err) {
	png_plot *plot = (png_plot *) spec->ptr;

	set_plot_has_y2_axis(plot, spec->flags & GPT_Y2AXIS);
	redisplay_edited_plot(plot);
	mark_session_changed();
    }

    spec->flags &= ~(GPT_FILL_SWITCH | GPT_ERR_SWITCH);
}

static void set_keyspec_sensitivity (GPT_SPEC *spec)
{
    gboolean state = FALSE;
    const char *p;
    int i;

    for (i=0; i<spec->n_lines; i++) {
	p = gtk_entry_get_text(GTK_ENTRY(linetitle[i]));
	if (p != NULL && *p != 0) {
	    state = TRUE;
	    break;
	}
    }

    gtk_widget_set_sensitive(keycombo, state);
}

#define TAB_MAIN_COLS 3

struct font_info {
    const char *fname;
    const char *showname;
};

static struct font_info ttf_fonts[] = {
    { "arial", "Arial", },
    { "georgia", "Georgia", },
#ifndef G_OS_WIN32
    { "luxirr", "Luxi Serif" },
    { "luxisr", "Luxi Sans" },
    { "Vera", "Vera" },
#endif
    { "tahoma", "Tahoma" },
    { "trebuc", "Trebuchet" },
    { "verdana", "Verdana" }
};

static const char *get_font_filename (const char *showname)
{
    int i, nfonts = sizeof ttf_fonts / sizeof ttf_fonts[0];

    for (i=0; i<nfonts; i++) {
	if (!strcmp(ttf_fonts[i].showname, showname)) {
	    return ttf_fonts[i].fname;
	}
    }
    return NULL;
}

#ifndef G_OS_WIN32

static int font_is_ok (const char *fname)
{
    static int pngterm;
    char cmd[64];
    int err;

    if (pngterm == 0) {
	pngterm = gnuplot_png_terminal();
    }

    if (pngterm == GP_PNG_CAIRO) {
	sprintf(cmd, "set term pngcairo font \"%s,10\"", fname);
    } else {
	sprintf(cmd, "set term png font %s 10", fname);
    }

    err = gnuplot_test_command(cmd);

    return err == 0;
}

#endif

static struct font_info *get_gnuplot_ttf_list (int *nf)
{
    static struct font_info *retlist = NULL;
    static int goodfonts = -1;

    if (goodfonts >= 0) {
	*nf = goodfonts;
    } else {
	int i, j, nfonts = sizeof ttf_fonts / sizeof ttf_fonts[0];

	retlist = malloc(nfonts * sizeof *retlist);
	if (retlist == NULL) return NULL;

	j = 0;
	for (i=0; i<nfonts; i++) {
#ifdef G_OS_WIN32
	    retlist[j++] = ttf_fonts[i];
#else
	    if (font_is_ok(ttf_fonts[i].fname)) {
		retlist[j++] = ttf_fonts[i];
	    }
#endif
	}
	goodfonts = j;
	*nf = goodfonts;
    }

    return retlist;
}

static int font_match (const char *ttfname, const char *pngfont)
{
    return !strncmp(ttfname, pngfont, strlen(ttfname));
}

static int get_point_size (const char *font)
{
    int pts;

    if (sscanf(font, "%*s %d\n", &pts) == 1) {
	return pts;
    } else {
	return 10;
    }
}

static void strip_lr (gchar *txt)
{
    gchar test[16];
    gchar *p;

    sprintf(test, "(%s)", I_("left"));
    p = strstr(txt, test);
    if (p != NULL) {
	*p = '\0';
    } else {
	sprintf(test, "(%s)", I_("right"));
	p = strstr(txt, test);
	if (p != NULL) {
	   *p = '\0';
	}
    } 
}

static void toggle_axis_selection (GtkWidget *w, GPT_SPEC *spec)
{
    int no_y2 = GTK_TOGGLE_BUTTON(w)->active;
    int i;

    for (i=0; i<spec->n_lines; i++) {
	if (yaxiscombo[i] != NULL) {
	    gtk_widget_set_sensitive(yaxiscombo[i], !no_y2);
	}
    }
}

/* re-establish the default plot colors and reset the
   color selection buttons accordingly */

static void color_default_callback (GtkWidget *w, GtkWidget *book)
{
    GtkWidget *button;
    char id[32];
    int i;

    sprintf(id, "color-button%d", BOXCOLOR);
    button = g_object_get_data(G_OBJECT(book), id);

    if (button != NULL) {
	graph_palette_reset(BOXCOLOR);
	color_patch_button_reset(button, BOXCOLOR);
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    sprintf(id, "color-button%d", i);
	    button = g_object_get_data(G_OBJECT(book), id);
	    if (button != NULL) {
		if (i == 0) {
		    graph_palette_reset(i);
		}
		color_patch_button_reset(button, i);
	    }
	}
    }
}

static void table_add_row (GtkWidget *tbl, int *rows, int cols)
{
    *rows += 1;
    gtk_table_resize(GTK_TABLE(tbl), *rows, cols);    
}

static void add_color_selector (int i, GtkWidget *tbl, int *rows,
				GtkWidget *notebook)
{
    GtkWidget *button, *hbox;
    GtkWidget *label;
    char str[32];

    table_add_row(tbl, rows, TAB_MAIN_COLS);

    hbox = gtk_hbox_new(FALSE, 2);

    if (i == BOXCOLOR) {
	strcpy(str, _("Fill color"));
    } else {
	sprintf(str, _("Color %d"), i + 1);
    }

    label = gtk_label_new(str);
    gtk_container_add(GTK_CONTAINER(hbox), label);
    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 0, 1, 
			      *rows - 1, *rows);
    gtk_widget_show(label);
    gtk_widget_show(hbox);

    hbox = gtk_hbox_new(FALSE, 2);
    button = color_patch_button(i);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 1, 2, 
			      *rows - 1, *rows);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(graph_color_selector), 
		     GINT_TO_POINTER(i));

    gtk_widget_show_all(button);
    gtk_widget_show(hbox);

    sprintf(str, "color-button%d", i);
    g_object_set_data(G_OBJECT(notebook), str, button);

    if (i == BOXCOLOR || i == BOXCOLOR - 1) {
	table_add_row(tbl, rows, TAB_MAIN_COLS);
	hbox = gtk_hbox_new(FALSE, 2);
	button = gtk_button_new_with_label(_("Reset to default"));
	gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 10);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(color_default_callback), 
			 notebook);
	gtk_table_attach(GTK_TABLE(tbl), hbox, 0, 2, *rows - 1, *rows,
			 GTK_FILL, 0, 0, 5);
	gtk_widget_show(button);
	gtk_widget_show(hbox);
    }
}

/* PNG anti-aliasing switch */

static void set_aa_status (GtkWidget *w, int *ok)
{
    *ok = GTK_TOGGLE_BUTTON(w)->active;

    gnuplot_png_set_use_aa(*ok);
}

static GtkWidget *gp_dialog_table (int rows, int cols, 
				   GtkWidget *vbox)
{
    GtkWidget *tbl;

    tbl = gtk_table_new(rows, cols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(vbox), tbl, FALSE, FALSE, 0);

    return tbl;
}

static GtkWidget *gp_dialog_vbox (void)
{
    GtkWidget *vbox;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);

    return vbox;
}

static GtkWidget *gp_page_vbox (GtkWidget *notebook, char *str)
{
    GtkWidget *vbox;
    GtkWidget *label;

    vbox = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_widget_show(vbox);
    label = gtk_label_new(str);
    gtk_widget_show(label);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label); 

    return vbox;
}

static void gpt_tab_main (GtkWidget *notebook, GPT_SPEC *spec) 
{
    static int aa_ok = 1;
    static int show_aa_check = -1;
    GtkWidget *label, *vbox, *tbl;
    int i, rows = 1;
    int kactive = 0;
    gchar *keypos[] = {
	"left top",
	"right top",
	"left bottom",
	"right bottom",
	"outside",
	"none",
	NULL
    };

    vbox = gp_page_vbox(notebook, _("Main"));

    tbl = gp_dialog_table(rows, TAB_MAIN_COLS, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 0) {
	    GtkWidget *entry;

	    if (i > 0) {
		table_add_row(tbl, &rows, TAB_MAIN_COLS);
	    }

	    label = gtk_label_new(_(gpt_titles[i].description));
	    gtk_table_attach_defaults(GTK_TABLE (tbl), 
				      label, 0, 1, rows-1, rows);
	    gtk_widget_show(label);

	    entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      entry, 1, TAB_MAIN_COLS, 
				      rows-1, rows);
				      
            if (spec->titles[i] != NULL && *spec->titles[i] != '\0') {
		gp_string_to_entry(entry, spec->titles[i]);
            }		

	    g_signal_connect(G_OBJECT(entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     spec);

	    gtk_widget_show(entry);
	    gpt_titles[i].widget = entry;
	}
    }

    /* specify position of plot key or legend */
    table_add_row(tbl, &rows, TAB_MAIN_COLS);
    label = gtk_label_new(_("key position"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, rows-1, rows);
    gtk_widget_show(label);

    keycombo = gtk_combo_box_new_text();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      keycombo, 1, TAB_MAIN_COLS, rows-1, rows);
    for (i=0; keypos[i] != NULL; i++) {
	gtk_combo_box_append_text(GTK_COMBO_BOX(keycombo), keypos[i]);
	if (!strcmp(keypos[i], spec->keyspec)) {
	    kactive = i;
	}
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(keycombo), kactive);
    gtk_widget_show(keycombo);

    /* choice of fitted line type, if appropriate */
    if (spec->fit != PLOT_FIT_NA) {
	table_add_row(tbl, &rows, TAB_MAIN_COLS);

	label = gtk_label_new(_("fitted line"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, rows-1, rows);
	gtk_widget_show(label);

	fitcombo = gtk_combo_box_new_text();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  fitcombo, 1, TAB_MAIN_COLS, rows-1, rows);
	for (i=0; fittype_strings[i] != NULL; i++) {
	    gtk_combo_box_append_text(GTK_COMBO_BOX(fitcombo), _(fittype_strings[i]));
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(fitcombo), spec->fit);
	g_signal_connect(G_OBJECT(fitcombo), "changed",
			 G_CALLBACK(fit_type_changed), spec);
	gtk_widget_show(fitcombo);
    } else {
	fitcombo = NULL;
    }

    /* give option of removing top & right border */
    if (!(spec->flags & GPT_Y2AXIS)) { 
	y2_check = NULL;
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	border_check = gtk_check_button_new_with_label(_("Show full border"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  border_check, 0, TAB_MAIN_COLS, 
				  rows-1, rows);
	if (!(spec->flags & GPT_MINIMAL_BORDER)) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(border_check),
					 TRUE);
	}	
	gtk_widget_show(border_check);
    } else {
	border_check = NULL;
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	y2_check = gtk_check_button_new_with_label(_("Use only one y axis"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  y2_check, 0, TAB_MAIN_COLS, 
				  rows-1, rows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(y2_check),
				     FALSE);
	g_signal_connect(G_OBJECT(y2_check), "clicked", 
			 G_CALLBACK(toggle_axis_selection), spec);
	gtk_widget_show(y2_check);
    }

    /* give option of showing all case markers */
    if (spec->flags & GPT_ALL_MARKERS_OK) { 
	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	markers_check = gtk_check_button_new_with_label(_("Show all data labels"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  markers_check, 0, TAB_MAIN_COLS, 
				  rows-1, rows);
	if (spec->flags & GPT_ALL_MARKERS) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(markers_check),
					 TRUE);
	}	
	gtk_widget_show(markers_check);
    } else {
	markers_check = NULL;
    }

    if (show_aa_check < 0) {
	show_aa_check = (gnuplot_png_terminal() == GP_PNG_GD2);
    }

    /* give option of suppressing anti-aliasing for PNGs */
    if (show_aa_check) {
	GtkWidget *aa_check;

	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	aa_check = gtk_check_button_new_with_label(_("Allow anti-aliasing of lines"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  aa_check, 0, TAB_MAIN_COLS, 
				  rows-1, rows);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(aa_check),
				     aa_ok);
	g_signal_connect(G_OBJECT(aa_check), "clicked", 
			 G_CALLBACK(set_aa_status), &aa_ok);
	gtk_widget_show(aa_check);
    }	

    /* set TT font, if gnuplot uses freetype */
    if (gnuplot_has_ttf(0)) {
	GtkWidget *ebox, *hsep;
	GList *fontnames = NULL;
	struct font_info *ttflist;
	const char *default_font = NULL;
	int nfonts = 0;

	ttflist = get_gnuplot_ttf_list(&nfonts);

	for (i=0; i<nfonts; i++) {
	    fontnames = g_list_append(fontnames, (gpointer) ttflist[i].showname);
	    if (font_match(ttflist[i].fname, gretl_png_font())) {
		default_font = ttflist[i].showname;
	    }
	}

	fontnames = g_list_append(fontnames, _("None"));
	if (default_font == NULL) {
	    default_font = _("None");
	}

	/* first a separator */
	table_add_row(tbl, &rows, TAB_MAIN_COLS);	
	hsep = gtk_hseparator_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  rows-1, rows);  
	gtk_widget_show(hsep);

	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	ebox = gtk_event_box_new();
	label = gtk_label_new(_("TrueType font"));
	gtk_container_add(GTK_CONTAINER(ebox), label);
	gtk_table_attach_defaults(GTK_TABLE (tbl), ebox, 0, 1, 
				  rows-1, rows);
	gtk_widget_show(label);
	gtk_widget_show(ebox);

	/* FIXME max length of font name? */

	ttfcombo = gtk_combo_new();
	gtk_entry_set_max_length(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), 15);
	gtk_table_attach_defaults(GTK_TABLE(tbl), ttfcombo, 1, 2, 
				  rows-1, rows);
	gtk_combo_set_popdown_strings(GTK_COMBO(ttfcombo), fontnames); 
	gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), default_font);
	gtk_entry_set_width_chars(GTK_ENTRY(GTK_COMBO(ttfcombo)->entry), 15);
	g_signal_connect(G_OBJECT(GTK_COMBO(ttfcombo)->entry), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_widget_show(ttfcombo);
	g_list_free(fontnames);

	ttfspin = gtk_spin_button_new_with_range(6, 24, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(ttfspin), 
				  get_point_size(gretl_png_font()));
	gtk_table_attach_defaults(GTK_TABLE(tbl), ttfspin, 2, 3, 
				  rows - 1, rows);
	gtk_widget_show(ttfspin);
    } else {
	ttfcombo = NULL;
	ttfspin = NULL;
    }

    if (gnuplot_png_terminal() != GP_PNG_OLD) {
	GtkWidget *hsep = gtk_hseparator_new();

	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  rows - 1, rows);  
	gtk_widget_show(hsep);

	if (frequency_plot_code(spec->code)) {
	    add_color_selector(BOXCOLOR, tbl, &rows, notebook);
	} else {
	    for (i=0; i<BOXCOLOR; i++) {
		add_color_selector(i, tbl, &rows, notebook);
	    }
	}
    }
}

static void gpt_save_as_callback (GtkWidget *button, GPT_SPEC *spec) 
{
    GtkWidget *w = GTK_COMBO(termcombo)->entry;
    const gchar *str;

    str = gtk_entry_get_text(GTK_ENTRY(w));
    if (str == NULL || *str == '\0') {
	return;
    }

    *spec->termtype = '\0';
    strncat(spec->termtype, str, sizeof spec->termtype - 1);

    file_selector(_("Save gnuplot graph"), SAVE_GNUPLOT, 
		  FSEL_DATA_MISC, spec);
}

static void gpt_tab_output (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl;
    GtkWidget *button;
    int i, tbl_len;
    GList *termlist = NULL;
    gchar *termtypes[] = {
	"postscript",
	"postscript color",
	"PDF",
	"fig",
	"latex",
	"png",
	"plot commands",
	NULL
    }; 
    int pdf_ok = gnuplot_pdf_terminal();

    for (i=0; termtypes[i] != NULL; i++) {
	if (!pdf_ok && !strcmp(termtypes[i], "PDF")) {
	    continue;
	}
	termlist = g_list_append(termlist, termtypes[i]);
    }
   
    vbox = gp_page_vbox(notebook, _("Output to file"));

    tbl_len = 1;
    tbl = gp_dialog_table(tbl_len, 2, vbox);
    gtk_widget_show(tbl);
   
    tbl_len++;
    label = gtk_label_new(_("output format"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);

    termcombo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), termcombo, 1, 2, 
			      tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(termcombo), termlist); 
    gtk_widget_show(termcombo);
    g_list_free(termlist); 

    /* button to generate output to file */
    button = gtk_button_new_from_stock(GTK_STOCK_SAVE_AS);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    tbl_len++;
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      button, 1, 2, tbl_len-1, tbl_len);
    g_signal_connect (G_OBJECT(button), "clicked", 
		      G_CALLBACK(gpt_save_as_callback), 
		      spec);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);  
}

static void linetitle_callback (GtkWidget *w, GPT_SPEC *spec)
{
    set_keyspec_sensitivity(spec);
}

static void gpt_tab_lines (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *hbox;
    int i, tbl_len, tbl_num, tbl_col;
    char label_text[32];
    GList *stylist = NULL;
    int do_scale_axis = 0;

    if (spec->code == PLOT_REGULAR && (spec->flags & GPT_TS)) {
	do_scale_axis = 1;
    }

    if (frequency_plot_code(spec->code)) {
	stylist = g_list_append(stylist, "boxes");
    }

    if (spec->flags & GPT_TS) {
	stylist = g_list_append(stylist, "lines");
	stylist = g_list_append(stylist, "points");
    } else {
	stylist = g_list_append(stylist, "points");
	stylist = g_list_append(stylist, "lines");
    }

    stylist = g_list_append(stylist, "linespoints"); 
    stylist = g_list_append(stylist, "impulses");
    stylist = g_list_append(stylist, "dots");
    stylist = g_list_append(stylist, "steps");

    vbox = gp_dialog_vbox();
    gtk_widget_show(vbox);

    label = gtk_label_new(_("Lines"));
    gtk_widget_show(label);

    if (spec->n_lines > 4) {
	GtkWidget *scroller;
	
	scroller = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				       GTK_POLICY_AUTOMATIC, 
				       GTK_POLICY_AUTOMATIC);
	gtk_widget_show(scroller);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					      vbox);
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), scroller, label); 
    } else {
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), vbox, label); 
    }  

    tbl_len = 1;
    tbl = gp_dialog_table(tbl_len, 3, vbox);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    for (i=0; i<spec->n_lines; i++) {
	/* identifier and key or legend text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sprintf(label_text, _("line %d: "), i + 1);
	label = gtk_label_new(label_text);
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	label = gtk_label_new(_("legend"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	linetitle[i] = gtk_entry_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linetitle[i], 2, 3, tbl_len-1, tbl_len);

	strip_lr(spec->lines[i].title);
	gp_string_to_entry(linetitle[i], spec->lines[i].title);

	g_signal_connect(G_OBJECT(linetitle[i]), "changed", 
			 G_CALLBACK(linetitle_callback), 
			 spec);
	g_signal_connect(G_OBJECT(linetitle[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_widget_show(linetitle[i]);

	/* line type or style */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("type"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	stylecombo[i] = gtk_combo_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  stylecombo[i], 2, 3, tbl_len-1, tbl_len);

	/* the errorbars and filledcurves styles are not exchangeable
	   with the others */
	if (!strcmp(spec->lines[i].style, "errorbars") &&
	    !gnuplot_has_style_fill()) {
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(stylecombo[i])->entry), 
			       spec->lines[i].style); 
	    gtk_widget_set_sensitive(stylecombo[i], FALSE);
	} else if (!strcmp(spec->lines[i].style, "errorbars") ||
		   !strcmp(spec->lines[i].style, "filledcurve")) {
	    GList *altsty = NULL;

	    altsty = g_list_append(altsty, "errorbars"); 
	    altsty = g_list_append(altsty, "filledcurve");
	    gtk_combo_set_popdown_strings(GTK_COMBO(stylecombo[i]), altsty); 
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(stylecombo[i])->entry), 
			       spec->lines[i].style); 
	    g_list_free(altsty);
	} else {
	    gtk_combo_set_popdown_strings(GTK_COMBO(stylecombo[i]), stylist); 
	    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(stylecombo[i])->entry), 
			       spec->lines[i].style); 
	} 

	gtk_widget_show(stylecombo[i]);	

	if (!do_scale_axis) {
	    linescale[i] = NULL;
	    yaxiscombo[i] = NULL;
	} else {
	    /* scale factor for data? */
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	    label = gtk_label_new(_("scale"));
	    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      label, 1, 2, tbl_len-1, tbl_len);
	    gtk_widget_show(label);

	    linescale[i] = gtk_entry_new();
	    gtk_entry_set_max_length(GTK_ENTRY(linescale[i]), 6);
	    gtk_entry_set_text(GTK_ENTRY(linescale[i]), spec->lines[i].scale);
	    gtk_entry_set_width_chars(GTK_ENTRY(linescale[i]), 6);
	    g_signal_connect(G_OBJECT(linescale[i]), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     spec);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      linescale[i], 2, 3, tbl_len-1, tbl_len);
	    gtk_widget_show(linescale[i]);

	    /* use left or right y axis? */
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	    label = gtk_label_new(_("y axis"));
	    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      label, 1, 2, tbl_len-1, tbl_len);
	    gtk_widget_show(label);

	    yaxiscombo[i] = gtk_combo_box_new_text();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      yaxiscombo[i], 2, 3, tbl_len-1, tbl_len);
	    gtk_combo_box_append_text(GTK_COMBO_BOX(yaxiscombo[i]), "left");
	    gtk_combo_box_append_text(GTK_COMBO_BOX(yaxiscombo[i]), "right");
	    gtk_combo_box_set_active(GTK_COMBO_BOX(yaxiscombo[i]), 
				     (spec->lines[i].yaxis == 1)? 0 : 1);
	    gtk_widget_show(yaxiscombo[i]);
	}

	/* line-width adjustment */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("line width"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	hbox = gtk_hbox_new(FALSE, 5);
	linewidth[i] = gtk_spin_button_new_with_range(1, 6, 1);
	if (spec->lines[i].width > 1) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(linewidth[i]),
				      spec->lines[i].width);
	}
	g_signal_connect(G_OBJECT(linewidth[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_box_pack_start(GTK_BOX(hbox), linewidth[i], FALSE, FALSE, 0);
	gtk_widget_show(linewidth[i]);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, 3, 
				  tbl_len-1, tbl_len);
	gtk_widget_show(hbox);
    }

    g_list_free(stylist);
}

static void label_pos_to_entry (double *pos, GtkWidget *w)
{
    if (!na(pos[0]) && !na(pos[1])) {
	gchar *s = g_strdup_printf("%g %g", pos[0], pos[1]);

	gtk_entry_set_text(GTK_ENTRY(w), s);
	g_free(s);
    } else {
	gtk_entry_set_text(GTK_ENTRY(w), "");
    }
}

static void gpt_tab_labels (GtkWidget *notebook, GPT_SPEC *spec) 
{
    GtkWidget *label, *vbox, *tbl, *menu;
    int i, j, tbl_len, tbl_num, tbl_col;
    char label_text[32];
    png_plot *plot = (png_plot *) spec->ptr;

    vbox = gp_page_vbox(notebook, _("Labels"));
 
    tbl_len = 1;
    tbl = gp_dialog_table(tbl_len, 3, vbox);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    for (i=0; i<MAX_PLOT_LABELS; i++) {
	GtkWidget *hbox, *button, *image;
	GdkPixbuf *icon;

	/* label text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sprintf(label_text, _("label %d: "), i + 1);
	label = gtk_label_new(label_text);
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 0, 1, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	label = gtk_label_new(_("text"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	labeltext[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(labeltext[i]), PLOT_LABEL_TEXT_LEN);
	gp_string_to_entry(labeltext[i], spec->labels[i].text);
	gtk_entry_set_width_chars(GTK_ENTRY(labeltext[i]), PLOT_LABEL_TEXT_LEN);
	g_signal_connect (G_OBJECT(labeltext[i]), "activate", 
			  G_CALLBACK(apply_gpt_changes), 
			  spec);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  labeltext[i], 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(labeltext[i]);

	/* label placement */
	tbl_len++;

	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("position (X Y)"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	/* holder for entry and button */
	hbox = gtk_hbox_new(FALSE, 5);

	/* entry for coordinates */
	labelpos[i] = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(labelpos[i]), PLOT_LABEL_POS_LEN);
	label_pos_to_entry(spec->labels[i].pos, labelpos[i]);
	gtk_entry_set_width_chars(GTK_ENTRY(labelpos[i]), PLOT_LABEL_POS_LEN);
	g_signal_connect(G_OBJECT(labelpos[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_container_add(GTK_CONTAINER(hbox), labelpos[i]);
	gtk_widget_show(labelpos[i]);

	if (plot_is_mouseable(plot)) {
	    /* button to invoke mouse-assisted placement */
	    button = gtk_button_new();
	    g_object_set_data(G_OBJECT(button), "labelpos_entry", labelpos[i]);
	    g_signal_connect(G_OBJECT(button), "clicked",
			     G_CALLBACK(plot_label_position_click), spec->ptr);
	    icon = gdk_pixbuf_new_from_xpm_data((const char **) mini_mouse_xpm);
	    image = gtk_image_new_from_pixbuf(icon);
	    gtk_widget_set_size_request(button, 32, 26);
	    gtk_container_add (GTK_CONTAINER(button), image);
	    gtk_container_add(GTK_CONTAINER(hbox), button);
	    gtk_widget_show_all(button);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  hbox, 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(hbox);

	/* label justification */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("justification"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	labeljust[i] = gtk_option_menu_new();
	menu = gtk_menu_new();
	for (j=0; j<3; j++) {
	    GtkWidget *item;

	    item = gtk_menu_item_new_with_label(gp_justification_string(j));
	    gtk_menu_shell_append(GTK_MENU_SHELL(menu), item);
	}
	gtk_option_menu_set_menu(GTK_OPTION_MENU(labeljust[i]), menu);	
	gtk_option_menu_set_history(GTK_OPTION_MENU(labeljust[i]),
				    spec->labels[i].just);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  labeljust[i], 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show_all(labeljust[i]);	
    }
}

static void gpt_tab_XY (GtkWidget *notebook, GPT_SPEC *spec, gint axis) 
{
    png_plot *plot = (png_plot *) spec->ptr;
    GtkWidget *b1, *b2, *entry, *vbox, *tbl;
    GtkWidget *label = NULL;
    char *labelstr = NULL;
    int i, tbl_len;
   
    if (axis == 0) {
	labelstr = _("X-axis");
    } else if (axis == 1) {
	labelstr = _("Y-axis");
    } else if (axis == 2) {
	labelstr = _("Y2-axis");
    } else {
	return;
    }

    vbox = gp_page_vbox(notebook, labelstr);

    tbl_len = 1;
    tbl = gp_dialog_table(tbl_len, 2, vbox);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 1 + axis) {
	    GtkWidget *title_entry;

	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
            
	    label = gtk_label_new(_(gpt_titles[i].description));
	    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      label, 0, 1, tbl_len-1, tbl_len);
	    gtk_widget_show(label);

	    title_entry = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      title_entry, 1, 2, 
				      tbl_len-1, tbl_len);
	    gp_string_to_entry(title_entry, spec->titles[i]);

	    g_signal_connect(G_OBJECT(title_entry), "activate", 
			     G_CALLBACK(apply_gpt_changes), 
			     spec);

	    gtk_widget_show(title_entry);
	    gpt_titles[i].widget = title_entry;
	}
    } 

    if (spec->code != PLOT_REGULAR) return;

    axis_range[axis].ID = axis;

    /* axis range: "auto" button */
    tbl_len += 3;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);

    label = gtk_label_new("");
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      tbl_len-3, tbl_len-2);
    gtk_widget_show(label);
    b1 = gtk_radio_button_new_with_label(NULL, _("auto axis range"));
    g_signal_connect(G_OBJECT(b1), "clicked",
		     G_CALLBACK(flip_manual_range), 
		     GINT_TO_POINTER(axis));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 1, 
			      tbl_len-2, tbl_len-1);
    gtk_widget_show(b1);
    axis_range[axis].isauto = b1;

    /* axis range: manual range button */
    b2 = gtk_radio_button_new_with_label(gtk_radio_button_get_group 
					 (GTK_RADIO_BUTTON(b1)),
					 _("manual range:")); 
    g_signal_connect(G_OBJECT(b2), "clicked",
		     G_CALLBACK(flip_manual_range), 
		     GINT_TO_POINTER(axis));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 1, 
			      tbl_len-1, tbl_len);
    gtk_widget_show(b2);

    /* axis range min. entry */
    tbl_len++;
    label = gtk_label_new(_("minimum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, 
			      tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(entry), "");
    g_signal_connect(G_OBJECT(entry), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     spec);
    gtk_widget_show(entry);
    axis_range[axis].min = entry;

    /* axis range max. entry */
    tbl_len++;
    label = gtk_label_new(_("maximum"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
			      tbl_len-1, tbl_len);
    gtk_widget_show(label);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    entry = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, 
			      tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(entry), "");
    g_signal_connect(G_OBJECT(entry), "activate", 
		     G_CALLBACK(apply_gpt_changes), 
		     spec);
    gtk_widget_show(entry);
    axis_range[axis].max = entry;

    if (na(spec->range[axis][0])) {
	flip_manual_range(NULL, GINT_TO_POINTER(axis));
    } else {
	double_to_gp_entry(spec->range[axis][0], axis_range[axis].min);
	double_to_gp_entry(spec->range[axis][1], axis_range[axis].max);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(axis_range[axis].isauto), 
				     FALSE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    }

    /* axis scale: linear vs log? */

    if ((axis == 0 && plot_get_xmin(plot) >= 0.0) ||
	(axis == 1 && plot_get_ymin(plot) >= 0.0)) {
	GtkWidget *combo;
	GList *strs = NULL;

	strs = g_list_append(strs, "e");
	strs = g_list_append(strs, "2");
	strs = g_list_append(strs, "10");

	tbl_len += 3;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);

	label = gtk_label_new("");
	gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 
				  tbl_len-3, tbl_len-2);
	gtk_widget_show(label);
	b1 = gtk_radio_button_new_with_label(NULL, _("linear scale"));
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 1, 
				  tbl_len-2, tbl_len-1);
	gtk_widget_show(b1);

	b2 = gtk_radio_button_new_with_label(gtk_radio_button_get_group 
					     (GTK_RADIO_BUTTON(b1)),
					     _("logarithmic scale, base:")); 
	gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 1, 
				  tbl_len-1, tbl_len);
	gtk_widget_show(b2);    

	combo = gtk_combo_new();
	gtk_combo_set_popdown_strings(GTK_COMBO(combo), strs); 
	g_list_free(strs);
	gtk_table_attach_defaults(GTK_TABLE(tbl), combo, 1, 2, 
				  tbl_len-1, tbl_len);
	entry = GTK_COMBO(combo)->entry;
	g_signal_connect(G_OBJECT(entry), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_widget_show(combo);
	axis_range[axis].lbase = entry;

	if (spec->logbase[axis] > 1.1) {
	    gchar *txt = g_strdup_printf("%g", spec->logbase[axis]);

	    gtk_entry_set_text(GTK_ENTRY(entry), txt);
	    g_free(txt);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
	} else {
	    gtk_widget_set_sensitive(entry, FALSE);
	}

	g_signal_connect(G_OBJECT(b1), "clicked", G_CALLBACK(disable_lbase), 
			 entry);
	g_signal_connect(G_OBJECT(b2), "clicked", G_CALLBACK(enable_lbase), 
			 entry);
    }
}

void close_gnuplot_dialog (GtkWidget *w, gpointer p)
{
    free(linetitle);
    free(stylecombo);
    free(yaxiscombo);
    free(linescale);
    free(linewidth);

    linetitle = NULL;
    stylecombo = NULL;
    yaxiscombo = NULL;
    linescale = NULL;
    linewidth = NULL;

    gtk_widget_destroy(GTK_WIDGET(p));
}

static int gpt_allocate_widgets (int n)
{
    int i;

    if (n == 0) {
	return 0;
    }

    linetitle = malloc(n * sizeof *linetitle);
    stylecombo = malloc(n * sizeof *stylecombo);
    yaxiscombo = malloc(n * sizeof *yaxiscombo);
    linescale = malloc(n * sizeof *linescale);
    linewidth = malloc(n * sizeof *linewidth);
    
    if (linetitle == NULL || stylecombo == NULL ||
	yaxiscombo == NULL || linescale == NULL ||
	linewidth == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	linetitle[i] = NULL;
	stylecombo[i] = NULL;
	yaxiscombo[i] = NULL;
	linescale[i] = NULL;
	linewidth[i] = NULL;
    }

    return 0;
}

int show_gnuplot_dialog (GPT_SPEC *spec) 
{
    png_plot *plot = (png_plot *) spec->ptr;
    GtkWidget *button, *notebook;
    GtkWidget *hbox;
    int i;

    if (gpt_control != NULL) {
	errbox(_("You can only have one plot controller open\n"
		 "at any given time"));
	return 1;
    }

    if (gpt_allocate_widgets(spec->n_lines)) {
	nomem();
	return 1;
    }

    for (i=0; i<MAX_AXES; i++) {
	axis_range[i].isauto = NULL;
	axis_range[i].lbase = NULL;
    }

    for (i=0; i<NTITLES; i++) {
	gpt_titles[i].widget = NULL;
    }

    gpt_control = gretl_dialog_new(_("gretl plot controls"), NULL, 0);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 2);
    gtk_dialog_set_has_separator(GTK_DIALOG(gpt_control), FALSE);

    if (plot != NULL) {
	gtk_window_set_transient_for(GTK_WINDOW(gpt_control), 
				     GTK_WINDOW(plot_get_shell(plot)));
    }

    g_signal_connect(G_OBJECT(gpt_control), "destroy",
		     G_CALLBACK(close_plot_controller), 
		     (gpointer *) spec);
   
    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 
		       notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    gpt_tab_main(notebook, spec);
    gpt_tab_XY(notebook, spec, 0);
    gpt_tab_XY(notebook, spec, 1);

    if (spec->flags & GPT_Y2AXIS) {
	gpt_tab_XY(notebook, spec, 2);
    }

    if (spec->lines != NULL) {
	gpt_tab_lines(notebook, spec);
    }
    gpt_tab_labels(notebook, spec); 
    gpt_tab_output(notebook, spec);

    hbox = GTK_DIALOG(gpt_control)->action_area;

    /* "Apply" button */
    button = apply_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), spec);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    /* "OK" button (apply and close) */
    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(apply_gpt_changes), spec);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(close_gnuplot_dialog), gpt_control);
    gtk_widget_show(button);

    /* Close button (do not apply changes) */
    button = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    GTK_WIDGET_SET_FLAGS(button, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(close_gnuplot_dialog), gpt_control);
    gtk_widget_show(button);

    /* Help button */
    context_help_button(hbox, GR_PLOT);

    set_keyspec_sensitivity(spec);

    gtk_widget_show(gpt_control);

    return 0;
}

void raise_gpt_control_window (void)
{
    if (gpt_control != NULL) {
	gdk_window_raise(gpt_control->window);
    }
}

void destroy_gpt_control_window (void)
{
    gtk_widget_destroy(gpt_control);
}

