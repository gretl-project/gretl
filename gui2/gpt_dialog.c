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
#include "gppoints.h"

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
static GtkWidget **lineformula;
static GtkWidget **stylecombo;
static GtkWidget **yaxiscombo;
static GtkWidget **linescale;
static GtkWidget **linewidth;

static GtkWidget *fitformula;
static GtkWidget *fitlegend;

static GtkWidget *labeltext[MAX_PLOT_LABELS];
static GtkWidget *labeljust[MAX_PLOT_LABELS];
static GtkWidget *labelpos[MAX_PLOT_LABELS];

static GtkWidget *gpt_control;
static GtkWidget *keycombo;
static GtkWidget *fitcombo;
static GtkWidget *border_check;
static GtkWidget *markers_check;
static GtkWidget *y2_check;
static GtkWidget *ttfcombo;
static GtkWidget *ttfspin;

static int gui_nlines;

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

static void gpt_tab_lines (GtkWidget *notebook, GPT_SPEC *spec, int ins);

static int gpt_expand_widgets (int n);

static void widget_set_int (GtkWidget *w, const gchar *key, gint val)
{
    g_object_set_data(G_OBJECT(w), key, GINT_TO_POINTER(val));
}

static int widget_get_int (GtkWidget *w, const char *key)
{
    return GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), key));
}

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

/* graph color selection apparatus */

#define XPMROWS 19
#define XPMCOLS 17

#define scale_round(v) ((v) * 255.0 / 65535.0)

static GtkWidget *get_image_for_color (const gretlRGB *color)
{
    static char **xpm = NULL;
    GdkPixbuf *icon;
    GtkWidget *image;
    char colstr[8] = {0};
    int i;

    if (color == NULL) {
	return NULL;
    }

    if (xpm == NULL) {
	xpm = strings_array_new_with_length(XPMROWS, XPMCOLS);
	if (xpm == NULL) {
	    return NULL;
	}

	/* common set-up */
	strcpy(xpm[0], "16 16 2 1");
	strcpy(xpm[1], "X      c #000000");
	strcpy(xpm[2], ".      c #000000");
	strcpy(xpm[3], "................");

	for (i=4; i<XPMROWS-1; i++) {
	    strcpy(xpm[i], ".XXXXXXXXXXXXXX.");
	}

	strcpy(xpm[XPMROWS-1], "................");
    }

    /* write in the specific color we want */
    print_rgb_hash(colstr, color);
    for (i=0; i<6; i++) {
	xpm[1][10+i] = colstr[i+1];
    }    

    icon = gdk_pixbuf_new_from_xpm_data((const char **) xpm);
    image = gtk_image_new_from_pixbuf(icon);
    
    return image;
}

static void color_select_callback (GtkWidget *button, GtkWidget *w)
{
    GtkWidget *csel;
    GtkWidget *color_button, *image;
    GdkColor gcolor;
    gpointer data;
    gretlRGB rgb;
    gint i;

    color_button = g_object_get_data(G_OBJECT(w), "color_button");
    csel = GTK_COLOR_SELECTION_DIALOG(w)->colorsel;

    gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(csel), &gcolor);

    rgb.r = (unsigned char) (scale_round(gcolor.red));
    rgb.g = (unsigned char) (scale_round(gcolor.green));
    rgb.b = (unsigned char) (scale_round(gcolor.blue));

    i = widget_get_int(w, "colnum");
    data = g_object_get_data(G_OBJECT(color_button), "plotspec");
    
    if (data != NULL) {
	gretlRGB *prgb = malloc(sizeof *prgb);

	if (prgb != NULL) {
	    prgb->r = rgb.r;
	    prgb->g = rgb.g;
	    prgb->b = rgb.b;
	    g_object_set_data_full(G_OBJECT(color_button), "rgb",
				   prgb, free);
	}
    } else {
	set_graph_palette(i, rgb);
	update_persistent_graph_colors();
    } 

    /* update the "image" widget */
    image = g_object_get_data(G_OBJECT(color_button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(&rgb);
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(color_button), image);
    g_object_set_data(G_OBJECT(color_button), "image", image);
  
    gtk_widget_destroy(w);
}

static void color_cancel (GtkWidget *button, GtkWidget *w)
{
    gtk_widget_destroy(w);
}

/* reset color patch button after selecting the option to
   restore the default plot colors */

static void color_patch_button_reset (GtkWidget *button, int cnum)
{
    GtkWidget *image;

    image = g_object_get_data(G_OBJECT(button), "image");
    gtk_widget_destroy(image);
    image = get_image_for_color(get_graph_color(cnum));
    gtk_widget_show(image);
    gtk_container_add(GTK_CONTAINER(button), image);
    g_object_set_data(G_OBJECT(button), "image", image);

    if (cnum == BOXCOLOR || cnum == BOXCOLOR - 1) {
	update_persistent_graph_colors();
    }
}

static void graph_color_selector (GtkWidget *w, gpointer p)
{
    GPT_SPEC *spec;
    GtkWidget *cdlg;
    GtkWidget *button;
    gint i = GPOINTER_TO_INT(p);
    char colstr[8];
    GdkColor gcolor;

    spec = g_object_get_data(G_OBJECT(w), "plotspec");

    if (spec != NULL && spec->lines[i].rgb[0] != '\0') {
	strcpy(colstr, spec->lines[i].rgb);
    } else {
	const gretlRGB *rgb = get_graph_color(i);

	if (rgb == NULL) {
	    fprintf(stderr, "graph_color_selector: got NULL rgb\n");
	    return;
	}
	print_rgb_hash(colstr, rgb);
    }

    gdk_color_parse(colstr, &gcolor);

    cdlg = gtk_color_selection_dialog_new(_("gretl: graph color selection"));

    widget_set_int(cdlg, "colnum", i);
    g_object_set_data(G_OBJECT(cdlg), "color_button", w);

    gtk_color_selection_set_current_color(GTK_COLOR_SELECTION
					  (GTK_COLOR_SELECTION_DIALOG(cdlg)->colorsel),
					  &gcolor);					  

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->ok_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_select_callback), cdlg);

    button = GTK_COLOR_SELECTION_DIALOG(cdlg)->cancel_button;
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(color_cancel), cdlg);

    gtk_widget_show(cdlg);
    gtk_window_set_modal(GTK_WINDOW(cdlg), TRUE);
}

static GtkWidget *color_patch_button (int i)
{
    GtkWidget *image, *button;

    image = get_image_for_color(get_graph_color(i));

    if (image == NULL) {
	button = gtk_button_new_with_label(_("Select color"));
    } else {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
	g_object_set_data(G_OBJECT(button), "image", image);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(graph_color_selector), 
			 GINT_TO_POINTER(i));
    }	

    return button;
}

static int style_from_line_number (GPT_SPEC *spec, int i)
{
    int j, lt = spec->lines[i].type;

    if (lt > 0) {
	j = lt - 1;
    } else if (lt == LT_NONE) {
	j = i;
    } else {
	j = LT_NONE;
    }

    return j;
}

static GtkWidget *line_color_button (GPT_SPEC *spec, int i)
{
    GtkWidget *image, *button;
    int j = style_from_line_number(spec, i);

    if (j == LT_NONE) {
	return NULL;
    }

    if (spec->lines[j].rgb[0] != '\0') {
	gretlRGB color;

	gretl_rgb_get(&color, spec->lines[j].rgb); 
	image = get_image_for_color(&color);
    } else {
	image = get_image_for_color(get_graph_color(j));
    }

    if (image == NULL) {
	/* failsafe -- shouldn't happen */
	button = gtk_button_new_with_label(_("Select color"));
    } else {
	button = gtk_button_new();
	gtk_container_add(GTK_CONTAINER(button), image);
	g_object_set_data(G_OBJECT(button), "image", image);
	g_object_set_data(G_OBJECT(button), "plotspec", spec);
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(graph_color_selector), 
			 GINT_TO_POINTER(j));
    }	

    return button;
}

static void maybe_apply_line_color (GtkWidget *w, GPT_SPEC *spec,
				    int i)
{
    gretlRGB *rgb = NULL;
    GtkWidget *cb;
    int j;

    cb = g_object_get_data(G_OBJECT(w), "colorsel");

    if (cb != NULL) {
	rgb = g_object_get_data(G_OBJECT(cb), "rgb");
    }

    if (rgb != NULL) {
	j = style_from_line_number(spec, i);
	print_rgb_hash(spec->lines[j].rgb, rgb);
    }
}

/* end graph color selection apparatus */

static GdkPixbuf *get_pixbuf_for_line (int dots)
{
    static char **xpm = NULL;
    int i;

    if (xpm == NULL) {
	xpm = strings_array_new_with_length(14, 31);
	if (xpm == NULL) {
	    return NULL;
	}

	/* common set-up */
	strcpy(xpm[0], "30 11 2 1");
	strcpy(xpm[1], "  c None");
	strcpy(xpm[2], ". c #000000");

	for (i=3; i<14; i++) {
	    strcpy(xpm[i], "                              ");
	}
    }

    /* write in the specific pattern we want */
    if (dots) {
	strcpy(xpm[8], ".  .  .  .  .  .  .  .  .  .  ");
    } else {
	strcpy(xpm[8], "............................  ");
    }    

    return gdk_pixbuf_new_from_xpm_data((const char **) xpm);
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
    const gchar *s;

    *targ = '\0';

    g_return_if_fail(GTK_IS_ENTRY(w));

    s = gtk_entry_get_text(GTK_ENTRY(w));
    if (s != NULL && *s != '\0') {
	strncat(targ, s, n - 1);
    }
}

static void combo_to_gp_string (GtkWidget *w, char *targ, size_t n)
{
    gchar *s;

    *targ = '\0';

    g_return_if_fail(GTK_IS_COMBO_BOX(w));

    s = gtk_combo_box_get_active_text(GTK_COMBO_BOX(w));

    if (s != NULL && *s != '\0') {
	strncat(targ, s, n - 1);
    }

    g_free(s);
}

static void fittype_from_combo (GtkComboBox *box, GPT_SPEC *spec)
{
    int oldfit = widget_get_int(GTK_WIDGET(box), "oldfit");
    FitType f = gtk_combo_box_get_active(box);

    if (f == oldfit) {
	/* no change */
	return;
    }

    if (f == PLOT_FIT_OLS || f == PLOT_FIT_QUADRATIC || 
	f == PLOT_FIT_INVERSE || f == PLOT_FIT_LOESS) {
	plotspec_add_fit(spec, f);
	spec->flags &= ~GPT_FIT_HIDDEN;
    } else if (f == PLOT_FIT_NONE) {
	if (spec->n_lines >= 2) {
	    spec->flags |= GPT_FIT_HIDDEN;
	}
	spec->fit = f;
    }

    widget_set_int(GTK_WIDGET(box), "oldfit", f);
}

static gboolean fit_type_changed (GtkComboBox *box, GPT_SPEC *spec)
{
    const char *s1 = spec->yvarname;
    const char *s2 = spec->xvarname;
    gchar *title = NULL;
    FitType f = PLOT_FIT_NONE;

    if (*s1 == '\0' || *s2 == '\0') {
	return FALSE;
    }

    f = gtk_combo_box_get_active(box);

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

    if (title != NULL) {
	gtk_entry_set_text(GTK_ENTRY(gpt_titles[0].widget), title);
	g_free(title);
    }

    /* also re-jig the "Lines" tab entries for the fitted
       line */

    if (fitformula != NULL && fitlegend != NULL) {
	fittype_from_combo(box, spec);
	gtk_entry_set_text(GTK_ENTRY(fitformula), spec->lines[1].formula);
	gtk_entry_set_text(GTK_ENTRY(fitlegend), spec->lines[1].title);
    }
    
    return FALSE;
}

static void dot_callback (GtkComboBox *box, GPT_SPEC *spec)
{
    int i = widget_get_int(GTK_WIDGET(box), "linenum");
    int d = gtk_combo_box_get_active(box);
    GtkWidget *w;

    /* flip between colored solid line and black dotted */

    spec->lines[i].type = (d > 0)? 0 : LT_NONE;

    w = g_object_get_data(G_OBJECT(box), "colorsel");
    gtk_widget_set_sensitive(w, d == 0);
    w = g_object_get_data(G_OBJECT(box), "color-label");
    gtk_widget_set_sensitive(w, d == 0);
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

static void maybe_set_point_type (GPT_LINE *line, GtkWidget *w, int i)
{
    GtkWidget *ptsel = g_object_get_data(G_OBJECT(w), "pointsel");

    if (ptsel != NULL && GTK_WIDGET_SENSITIVE(ptsel)) {
	int pt = gtk_combo_box_get_active(GTK_COMBO_BOX(ptsel));
	int ptdef = (line->type == LT_NONE)? i : line->type - 1;

	if (pt != ptdef) {
	    /* point-type is not just the default */
	    line->ptype = pt + 1;
	}
    }
}

enum {
    ERRORBARS = 1,
    FILLEDCURVE
};

static void apply_gpt_changes (GtkWidget *w, GPT_SPEC *spec) 
{
    int suppress_y2 = 0;
    int i, k, err = 0;

    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].widget != NULL) {
	    entry_to_gp_string(gpt_titles[i].widget, spec->titles[i], 
			       sizeof spec->titles[0]);
	}
    }

    if (keycombo != NULL) {
	gchar *s = gtk_combo_box_get_active_text(GTK_COMBO_BOX(keycombo));

	strcpy(spec->keyspec, s);
	g_free(s);
    }

    spec->flags &= ~GPT_Y2AXIS;

    if (y2_check != NULL && GTK_TOGGLE_BUTTON(y2_check)->active) {
	suppress_y2 = 1;
    } 

    for (i=0; i<gui_nlines; i++) {
	GPT_LINE *line = &spec->lines[i];

	line->yaxis = 1;
	if (!suppress_y2 && yaxiscombo[i] != NULL) {
	    gchar *s = 
		gtk_combo_box_get_active_text(GTK_COMBO_BOX(yaxiscombo[i]));

	    if (!strcmp(s, "right")) {
		line->yaxis = 2;	
	    }
	    if (line->yaxis == 2) {
		spec->flags |= GPT_Y2AXIS;
	    }
	    g_free(s);
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
	for (i=0; i<gui_nlines; i++) {
	    GPT_LINE *line = &spec->lines[i];

	    if (stylecombo[i] != NULL) {
		int oldalt = 0;

		if (!strncmp(line->style, "filled", 6)) {
		    oldalt = FILLEDCURVE;
		} else if (!strncmp(line->style, "error", 5)) {
		    oldalt = ERRORBARS;
		}
		combo_to_gp_string(stylecombo[i], line->style, 
				   sizeof spec->lines[0].style);
		if (oldalt == FILLEDCURVE &&
		    !strncmp(line->style, "error", 5)) {
		    spec->flags |= GPT_ERR_SWITCH;
		    spec->flags &= ~GPT_FILL_SWITCH;
		} else if (oldalt == ERRORBARS &&
			   !strncmp(line->style, "filled", 6)) {
		    spec->flags |= GPT_FILL_SWITCH;
		    spec->flags &= ~GPT_ERR_SWITCH;
		}
		maybe_set_point_type(line, stylecombo[i], i);
	    }
	    if (linetitle[i] != NULL) {
		entry_to_gp_string(linetitle[i], 
				   line->title, 
				   sizeof spec->lines[0].title);
	    }
	    if (lineformula[i] != NULL && GTK_WIDGET_IS_SENSITIVE(lineformula[i])) {
		entry_to_gp_string(lineformula[i], 
				   line->formula, 
				   sizeof spec->lines[0].formula);
	    }
	    if (linescale[i] != NULL) {
		entry_to_gp_string(linescale[i], 
				   line->scale, 
				   sizeof spec->lines[0].scale);
	    }
	    if (linewidth[i] != NULL) {
		line->width = 
		    gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(linewidth[i]));
		maybe_apply_line_color(linewidth[i], spec, i);
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
	    gtk_combo_box_get_active(GTK_COMBO_BOX(labeljust[i]));
    } 

    if (!err && border_check != NULL) {
	if (GTK_TOGGLE_BUTTON(border_check)->active) {
	    /* full border */
	    spec->border = GP_BORDER_DEFAULT;
	} else {
	    /* left and bottom only */
	    spec->border = 3;
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
	gchar *tmp = gtk_combo_box_get_active_text(GTK_COMBO_BOX(ttfcombo));
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
	g_free(tmp);
    }

    if (!err && fitcombo != NULL && GTK_WIDGET_IS_SENSITIVE(fitcombo)) {
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

    for (i=0; i<gui_nlines; i++) {
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
    { "FreeSans", "Free Sans" },
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
    gchar *cmd;
    int err;

    if (pngterm == 0) {
	pngterm = gnuplot_png_terminal();
    }

    if (pngterm == GP_PNG_CAIRO) {
	cmd = g_strdup_printf("set term pngcairo font \"%s,10\"", fname);
    } else {
	cmd = g_strdup_printf("set term png font %s 10", fname);
    }

    err = gnuplot_test_command(cmd);
    g_free(cmd);

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

    for (i=0; i<gui_nlines; i++) {
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

static void add_color_selector (int i, GtkWidget *tbl, int cols,
				int *rows, GtkWidget *notebook)
{
    static int r0;
    int collen = (N_GP_COLORS - 1) / 2;
    GtkWidget *button, *hbox;
    GtkWidget *label;
    int row, cmin, cmax;
    char str[32];

    if (i == 0) {
	/* get baseline table row for color selectors */
	r0 = *rows;
    }

    if (i < collen || i == BOXCOLOR) {
	table_add_row(tbl, rows, cols);
	cmin = 0;
	cmax = 1;
	row = *rows;
    } else {
	/* place selector to the right */
	row = r0 + 1 + i - collen;
	cmin = 1;
	cmax = 2;
    }

    hbox = gtk_hbox_new(FALSE, 2);

    if (i == BOXCOLOR) {
	strcpy(str, _("Fill color"));
    } else {
	sprintf(str, _("Color %d"), i + 1);
    }

    label = gtk_label_new(str);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);

    button = color_patch_button(i);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);
    gtk_widget_show_all(button);
    gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, cmin, cmax, 
			      row - 1, row);
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
	widget_set_int(fitcombo, "oldfit", spec->fit);
	g_signal_connect(G_OBJECT(fitcombo), "changed",
			 G_CALLBACK(fit_type_changed), spec);
	gtk_widget_show(fitcombo);
    } else {
	fitcombo = NULL;
    }

    border_check = y2_check = NULL;

    /* give option of removing/adding top & right border? */
    if (!(spec->flags & GPT_Y2AXIS)) { 
	if (spec->border == GP_BORDER_DEFAULT || spec->border == 3) {
	    table_add_row(tbl, &rows, TAB_MAIN_COLS);
	    border_check = gtk_check_button_new_with_label(_("Show full border"));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      border_check, 0, TAB_MAIN_COLS, 
				      rows-1, rows);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(border_check),
					 spec->border == GP_BORDER_DEFAULT);
	    gtk_widget_show(border_check);
	}
    } else {
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

    markers_check = NULL;

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
	GtkWidget *ebox, *hsep, *entry;
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

	ttfcombo = gtk_combo_box_entry_new_text();
	entry = gtk_bin_get_child(GTK_BIN(ttfcombo));
	gtk_entry_set_max_length(GTK_ENTRY(entry), 15);
	gtk_table_attach_defaults(GTK_TABLE(tbl), ttfcombo, 1, 2, 
				  rows-1, rows);
	set_combo_box_strings_from_list(GTK_COMBO_BOX(ttfcombo), fontnames); 
	gtk_entry_set_text(GTK_ENTRY(entry), default_font);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), 15);
	g_signal_connect(G_OBJECT(entry), "activate", 
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

    if (gnuplot_png_terminal() != GP_PNG_OLD && 
	frequency_plot_code(spec->code)) {
	GtkWidget *hsep = gtk_hseparator_new();

	table_add_row(tbl, &rows, TAB_MAIN_COLS);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hsep, 0, TAB_MAIN_COLS, 
				  rows - 1, rows);  
	gtk_widget_show(hsep);

	add_color_selector(BOXCOLOR, tbl, TAB_MAIN_COLS, &rows, 
			   notebook);
    }
}

static void linetitle_callback (GtkWidget *w, GPT_SPEC *spec)
{
    set_keyspec_sensitivity(spec);
}

struct new_line_info_ {
    GtkWidget *dlg;
    GtkWidget *formula_entry;
    GPT_SPEC *spec;
};

typedef struct new_line_info_ new_line_info;

static void gpt_tab_new_line (new_line_info *nlinfo) 
{
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *hbox;
    int tbl_len, tbl_num, tbl_col;
    char label_text[32];

    vbox = GTK_DIALOG(nlinfo->dlg)->vbox;
    hbox = gtk_hbox_new(FALSE, 5);

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    /* identifier and formula text */
    tbl_len++;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
    sprintf(label_text, _("line %d: "), gui_nlines + 1);
    label = gtk_label_new(label_text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(label);

    label = gtk_label_new(_("formula"));
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 1, 2, tbl_len-1, tbl_len);
    gtk_widget_show(label);

    nlinfo->formula_entry = gtk_entry_new();
    gtk_entry_set_text(GTK_ENTRY(nlinfo->formula_entry), "");
    gtk_entry_set_width_chars(GTK_ENTRY(nlinfo->formula_entry), 32);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      nlinfo->formula_entry, 2, 3, tbl_len-1, tbl_len);
    gtk_entry_set_activates_default(GTK_ENTRY(nlinfo->formula_entry), TRUE);
    gtk_widget_show(nlinfo->formula_entry);
}

static void real_add_line (GtkWidget *w, new_line_info *nlinfo)
{
    GPT_SPEC *spec = nlinfo->spec;
    GPT_LINE *line;
    GtkWidget *notebook;
    const gchar *s;
    gint pgnum;
    int err = 0;

    s = gtk_entry_get_text(GTK_ENTRY(nlinfo->formula_entry));
    if (s == NULL || *s == '\0') {
	errbox(_("No formula was given"));
	return;
    }

    err = plotspec_add_line(spec);
    if (err) {
	nomem();
	gtk_widget_destroy(nlinfo->dlg);
	return;
    }

    line = &spec->lines[spec->n_lines - 1];

    entry_to_gp_string(nlinfo->formula_entry, line->formula, GP_MAXFORMULA);

    strcpy(line->style, "lines");
    line->type = spec->n_lines; /* assign next line style */
    strcpy(line->scale, "NA");  /* mark as a non-data line */
    line->flags = GP_LINE_USER;

    notebook = g_object_get_data(G_OBJECT(gpt_control), "notebook");
    pgnum = widget_get_int(notebook, "lines_page");

    err = gpt_expand_widgets(gui_nlines + 1);

    if (err) {
	nomem();
	spec->n_lines -= 1;
    } else {
	/* re-fill the "lines" notebook page */
	gtk_notebook_remove_page(GTK_NOTEBOOK(notebook), pgnum);
	gpt_tab_lines(notebook, spec, pgnum);
	gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), pgnum);
	if (spec->n_lines == 2) {
	    /* user-defined line has taken the place of a potential fitted line */
	    spec->fit = PLOT_FIT_NA;
	    if (fitcombo != NULL) {
		gtk_widget_set_sensitive(fitcombo, FALSE);
	    }
	}
    }

    gtk_widget_destroy(nlinfo->dlg);
}

static void add_line_callback (GtkWidget *w, GPT_SPEC *spec)
{
    new_line_info *nlinfo;
    GtkWidget *hbox;
    GtkWidget *button;

    nlinfo = mymalloc(sizeof *nlinfo);
    if (nlinfo == NULL) {
	return;
    }

    nlinfo->spec = spec;

    nlinfo->dlg = gretl_dialog_new(_("Add line"), gpt_control, GRETL_DLG_BLOCK);

    gpt_tab_new_line(nlinfo);

    hbox = GTK_DIALOG(nlinfo->dlg)->action_area;

    button = cancel_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(delete_widget), nlinfo->dlg);
    gtk_widget_show(button);

    button = ok_button(hbox);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(real_add_line), nlinfo);
    gtk_widget_grab_default(button);
    gtk_widget_show(button);

    context_help_button(hbox, GPT_ADDLINE);

    gtk_widget_show(nlinfo->dlg);

    free(nlinfo);
}

static void remove_line (GtkWidget *w, GPT_SPEC *spec)
{
    GtkWidget *notebook;
    int i, pgnum;

    i = widget_get_int(w, "linenum");
    plotspec_delete_line(spec, i);
    gui_nlines -= 1;

    notebook = g_object_get_data(G_OBJECT(gpt_control), "notebook");
    pgnum = widget_get_int(notebook, "lines_page");

    /* re-fill the "lines" notebook page */
    gtk_notebook_remove_page(GTK_NOTEBOOK(notebook), pgnum);
    gpt_tab_lines(notebook, spec, pgnum);
    gtk_notebook_set_current_page(GTK_NOTEBOOK(notebook), pgnum);
}

static void line_remove_button (GtkWidget *tbl, int row, 
				GPT_SPEC *spec, int i)
{
    GtkWidget *button; 

    button = gtk_button_new_with_label("Remove");
    widget_set_int(button, "linenum", i);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(remove_line), spec);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      button, 0, 1, row - 1, row);
    gtk_widget_show(button);
}

static void print_line_label (GtkWidget *tbl, int row, GPT_SPEC *spec, 
			      int i)
{
    char label_text[32];
    GtkWidget *label; 

    if (spec->code == PLOT_BOXPLOTS) {
	if (i == 0) {
	    sprintf(label_text, "%s: ", _("quartiles"));
	} else if (i == 1) {
	    sprintf(label_text, "%s: ", _("median"));
	} else if (i == 2 && spec->n_lines == 3) {
	    sprintf(label_text, "%s: ", _("mean"));
	} else if (i == 2) {
	    sprintf(label_text, "%s: ", _("lower bound"));
	} else if (i == 3) {
	    sprintf(label_text, "%s: ", _("upper bound"));
	}
    } else {
	sprintf(label_text, _("line %d: "), i+1);
    }

    label = gtk_label_new(label_text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 0, 1, row - 1, row);
    gtk_widget_show(label);
}

static void print_field_label (GtkWidget *tbl, int row,
			       const gchar *text)
{
    GtkWidget *label; 

    label = gtk_label_new(text);
    gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      label, 1, 2, row - 1, row);
    gtk_widget_show(label);
}

static GtkWidget *dash_types_combo (void)
{
    GtkWidget *dotsel;
    GtkCellRenderer *cell;
    GtkListStore *store;
    GtkTreeIter iter;
    GdkPixbuf *pbuf;
    int i;

    store = gtk_list_store_new(1, GDK_TYPE_PIXBUF);

    for (i=0; i<2; i++) {
	gtk_list_store_append(store, &iter);
	pbuf = get_pixbuf_for_line(i);
	gtk_list_store_set(store, &iter, 0, pbuf, -1);
	g_object_unref(pbuf);
    }

    dotsel = gtk_combo_box_new_with_model(GTK_TREE_MODEL(store));
    cell = gtk_cell_renderer_pixbuf_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(dotsel), cell, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(dotsel), cell,
				   "pixbuf", 0, NULL);
    
    return dotsel;
}

static GtkWidget *point_types_combo (void)
{
    GtkWidget *ptsel;
    GtkCellRenderer *cell;
    GtkListStore *store;
    GtkTreeIter iter;
    GdkPixbuf *pbuf;
    int i;

    store = gtk_list_store_new(1, GDK_TYPE_PIXBUF);

    for (i=0; i<13; i++) {
	gtk_list_store_append(store, &iter);
	pbuf = gdk_pixbuf_from_pixdata(gppoints[i], FALSE, NULL);
	gtk_list_store_set(store, &iter, 0, pbuf, -1);
	g_object_unref(pbuf);
    }

    ptsel = gtk_combo_box_new_with_model(GTK_TREE_MODEL(store));
    cell = gtk_cell_renderer_pixbuf_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(ptsel), cell, FALSE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(ptsel), cell,
				   "pixbuf", 0, NULL);
    
    return ptsel;
}

static int line_get_point_style (GPT_LINE *line, int i)
{
    if (line->ptype > 0) {
	/* a specific point-style has been selected: convert
	   to zero-based */
	return line->ptype - 1;
    } else if (line->type == LT_NONE) {
	/* line type is set by placement of line in plot */
	return i;
    } else {
	/* a specific line-type has been selected: give the 
	   associated point-style */
	return line->type - 1;
    }
}

static int gpt_style_as_int (const char *sty, GList *list)
{
    GList *mylist = list;
    int i = 0;

    while (mylist != NULL) {
	if (!strcmp(sty, (const char *) mylist->data)) {
	    return i;
	}
	i++;
	mylist = mylist->next;
    }

    return 0;
}

#define has_point(s) (!strcmp(s, "points") || !strcmp(s, "linespoints"))

static void flip_pointsel (GtkWidget *box, GtkWidget *targ)
{
    GtkWidget *label = g_object_get_data(G_OBJECT(targ), "label");
    gchar *s = gtk_combo_box_get_active_text(GTK_COMBO_BOX(box));
    int hp = has_point(s);

    gtk_widget_set_sensitive(targ, hp);
    gtk_widget_set_sensitive(label, hp);
}

static void gpt_tab_lines (GtkWidget *notebook, GPT_SPEC *spec, int ins)
{
    GtkWidget *label, *tbl;
    GtkWidget *vbox, *hbox, *sep;
    GtkWidget *button, *page;
    int i, tbl_len, tbl_num, tbl_col;
    GList *stylist = NULL;
    int do_scale_axis = 0;
    int pgnum = -1;

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

    if (gui_nlines > 4) {
	GtkWidget *scroller;
	
	scroller = gtk_scrolled_window_new(NULL, NULL);
	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				       GTK_POLICY_AUTOMATIC, 
				       GTK_POLICY_AUTOMATIC);
	gtk_widget_show(scroller);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					      vbox);
	page = scroller;
    } else {
	page = vbox;
    }

    if (ins > 0) {
	pgnum = gtk_notebook_insert_page(GTK_NOTEBOOK(notebook), page, label, ins); 
    } else {
	pgnum = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
    }  

    widget_set_int(notebook, "lines_page", pgnum);

    tbl_len = 1;
    tbl = gp_dialog_table(tbl_len, 3, vbox);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    fitformula = fitlegend = NULL;

    for (i=0; i<gui_nlines; i++) {
	GPT_LINE *line = &spec->lines[i];
	int label_done = 0;

	if (line->formula[0] != '\0' || (line->flags & GP_LINE_USER)) {
	    /* the line has a formula (or is user-defined) */
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);

	    print_line_label(tbl, tbl_len, spec, i);
	    label_done = 1;

	    print_field_label(tbl, tbl_len, _("formula"));

	    lineformula[i] = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  lineformula[i], 2, 3, tbl_len-1, tbl_len);
	    strip_lr(line->formula);
	    gp_string_to_entry(lineformula[i], line->formula);
	    if (i == 1 && (spec->flags & GPT_AUTO_FIT)) {
		/* fitted formula: not GUI-editable */
		gtk_widget_set_sensitive(lineformula[i], FALSE);
		fitformula = lineformula[i];
	    } else {
		g_signal_connect(G_OBJECT(lineformula[i]), "activate", 
				 G_CALLBACK(apply_gpt_changes), 
				 spec);
	    }
	    gtk_widget_show(lineformula[i]);

	    linescale[i] = NULL;
	    yaxiscombo[i] = NULL;
	}

	/* identifier and key or legend text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);

	if (!label_done) {
	    print_line_label(tbl, tbl_len, spec, i);
	}

	if (line->flags & GP_LINE_USER) {
	    line_remove_button(tbl, tbl_len, spec, i);
	}

	print_field_label(tbl, tbl_len, _("legend"));

	linetitle[i] = gtk_entry_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linetitle[i], 2, 3, tbl_len-1, tbl_len);
	strip_lr(line->title);
	gp_string_to_entry(linetitle[i], line->title);
	g_signal_connect(G_OBJECT(linetitle[i]), "changed", 
			 G_CALLBACK(linetitle_callback), 
			 spec);
	g_signal_connect(G_OBJECT(linetitle[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_widget_show(linetitle[i]);
	if (i == 1 && (spec->flags & GPT_AUTO_FIT)) {
	    fitlegend = linetitle[i];
	}

	if (line->formula[0] != '\0') {
	    goto line_width_adj;
	}

	if (!strcmp(line->style, "candlesticks")) {
	    goto line_width_adj;
	}

	/* line type (lines, points, etc.) */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	label = gtk_label_new(_("type"));
	gtk_misc_set_alignment(GTK_MISC(label), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  label, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(label);

	hbox = gtk_hbox_new(FALSE, 5);
	stylecombo[i] = gtk_combo_box_new_text();
	gtk_box_pack_start(GTK_BOX(hbox), stylecombo[i], FALSE, FALSE, 0);

	/* the errorbars and filledcurves styles are not exchangeable
	   with the others */
	if (!strcmp(line->style, "errorbars") && !gnuplot_has_style_fill()) {
	    set_combo_box_default_text(GTK_COMBO_BOX(stylecombo[i]), 
				       line->style); 
	    gtk_widget_set_sensitive(stylecombo[i], FALSE);
	} else if (!strcmp(line->style, "errorbars") ||
		   !strcmp(line->style, "filledcurve")) {
	    GList *altsty = NULL;

	    altsty = g_list_append(altsty, "errorbars"); 
	    altsty = g_list_append(altsty, "filledcurve");
	    set_combo_box_strings_from_list(GTK_COMBO_BOX(stylecombo[i]), altsty);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(stylecombo[i]), 
				     gpt_style_as_int(line->style, altsty));
	    g_list_free(altsty);
	} else {
	    GtkWidget *ptsel = point_types_combo();
	    int lt = gpt_style_as_int(line->style, stylist);
	    int hp, pt = line_get_point_style(line, i);

	    set_combo_box_strings_from_list(GTK_COMBO_BOX(stylecombo[i]), stylist);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(stylecombo[i]), lt);
	    hp = has_point(line->style);
	    label = gtk_label_new(_("point"));
	    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
	    gtk_box_pack_start(GTK_BOX(hbox), ptsel, FALSE, FALSE, 0);
	    g_object_set_data(G_OBJECT(stylecombo[i]), "pointsel", ptsel);
	    g_object_set_data(G_OBJECT(ptsel), "label", label);
	    gtk_combo_box_set_active(GTK_COMBO_BOX(ptsel), pt);
	    gtk_widget_set_sensitive(label, hp);
	    gtk_widget_set_sensitive(ptsel, hp);
	    g_signal_connect(G_OBJECT(stylecombo[i]), "changed",
			     G_CALLBACK(flip_pointsel), ptsel);
	} 

	gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show_all(hbox);	

	if (!do_scale_axis) {
	    linescale[i] = NULL;
	    yaxiscombo[i] = NULL;
	} else {
	    /* scale factor for data? */
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	    print_field_label(tbl, tbl_len, _("scale"));

	    linescale[i] = gtk_entry_new();
	    gtk_entry_set_max_length(GTK_ENTRY(linescale[i]), 6);
	    gtk_entry_set_text(GTK_ENTRY(linescale[i]), line->scale);
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
	    print_field_label(tbl, tbl_len, _("y axis"));

	    yaxiscombo[i] = gtk_combo_box_new_text();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      yaxiscombo[i], 2, 3, tbl_len-1, tbl_len);
	    gtk_combo_box_append_text(GTK_COMBO_BOX(yaxiscombo[i]), "left");
	    gtk_combo_box_append_text(GTK_COMBO_BOX(yaxiscombo[i]), "right");
	    gtk_combo_box_set_active(GTK_COMBO_BOX(yaxiscombo[i]), 
				     (line->yaxis == 1)? 0 : 1);
	    gtk_widget_show(yaxiscombo[i]);
	}

    line_width_adj:

	/* line-width adjustment */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	print_field_label(tbl, tbl_len, _("line width"));

	hbox = gtk_hbox_new(FALSE, 5);
	linewidth[i] = gtk_spin_button_new_with_range(1, 6, 1);
	if (line->width > 1) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(linewidth[i]),
				      line->width);
	}
	g_signal_connect(G_OBJECT(linewidth[i]), "activate", 
			 G_CALLBACK(apply_gpt_changes), 
			 spec);
	gtk_box_pack_start(GTK_BOX(hbox), linewidth[i], FALSE, FALSE, 0);
	gtk_widget_show(linewidth[i]);
	gtk_table_attach_defaults(GTK_TABLE(tbl), hbox, 2, 3, 
				  tbl_len-1, tbl_len);
	/* line color adjustment */
	if (i < 6 && !frequency_plot_code(spec->code)) {
	    button = line_color_button(spec, i);
	    if (button != NULL) {
		label = gtk_label_new(_("color"));
		gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
		gtk_widget_show(label);
		gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
		g_object_set_data(G_OBJECT(linewidth[i]), "colorsel",
				  button);
		gtk_widget_show_all(button);
	    }
	} else {
	    label = NULL;
	    button = NULL;
	}

	if (line->flags & GP_LINE_USER) {
	    /* dotted option */
	    GtkWidget *dotcombo = dash_types_combo();

	    gtk_combo_box_set_active(GTK_COMBO_BOX(dotcombo), 0); 
	    gtk_widget_show(dotcombo);
	    widget_set_int(dotcombo, "linenum", i);
	    g_object_set_data(G_OBJECT(dotcombo), "colorsel", button);
	    g_object_set_data(G_OBJECT(dotcombo), "color-label", label);
	    g_signal_connect(G_OBJECT(dotcombo), "changed", 
			     G_CALLBACK(dot_callback), spec);
	    gtk_box_pack_start(GTK_BOX(hbox), dotcombo, FALSE, FALSE, 5);
	}

	gtk_widget_show(hbox);

	/* separator */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sep = gtk_hseparator_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), sep, 0, 3, 
				  tbl_len-1, tbl_len);
	gtk_widget_show(sep);
    }

    if (spec->code == PLOT_REGULAR) {
	/* button for adding a line (formula) */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	button = gtk_button_new_with_label(_("Add line..."));
	g_signal_connect(G_OBJECT(button), "clicked", 
			 G_CALLBACK(add_line_callback), 
			 spec);
	gtk_widget_show(button);
	gtk_table_attach_defaults(GTK_TABLE(tbl), button, 0, 1, 
				  tbl_len-1, tbl_len);
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
    GtkWidget *label, *vbox, *tbl;
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

	print_field_label(tbl, tbl_len, _("text"));

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

	print_field_label(tbl, tbl_len, _("position (X Y)"));

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
	    gtk_container_add(GTK_CONTAINER(button), image);
	    gtk_container_add(GTK_CONTAINER(hbox), button);
	    gtk_widget_show_all(button);
	}

	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  hbox, 2, 3, tbl_len-1, tbl_len);
	gtk_widget_show(hbox);

	/* label justification */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);

	print_field_label(tbl, tbl_len, _("justification"));

	labeljust[i] = gtk_combo_box_new_text();
	for (j=0; j<3; j++) {
	    gtk_combo_box_append_text(GTK_COMBO_BOX(labeljust[i]),
				      gp_justification_string(j));
	}
	gtk_combo_box_set_active(GTK_COMBO_BOX(labeljust[i]), 
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

    if (spec->code != PLOT_REGULAR) {
	return;
    }

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

	combo = gtk_combo_box_entry_new_text();
	set_combo_box_strings_from_list(GTK_COMBO_BOX(combo), strs);
	g_list_free(strs);
	gtk_table_attach_defaults(GTK_TABLE(tbl), combo, 1, 2, 
				  tbl_len-1, tbl_len);
	entry = gtk_bin_get_child(GTK_BIN(combo));
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

static void gpt_tab_palette (GtkWidget *notebook) 
{
    GtkWidget *vbox, *label, *tbl;
    int i, rows = 1;

    vbox = gp_page_vbox(notebook, _("Palette"));
    tbl = gp_dialog_table(1, 2, vbox);

    label = gtk_label_new("These colors will be used unless overridden\n"
			  "by graph-specific choices\n");
    gtk_table_attach(GTK_TABLE(tbl), label, 0, 2, 0, 1, 0, 0, 0, 0);
    gtk_widget_show(label);

    for (i=0; i<BOXCOLOR; i++) {
	add_color_selector(i, tbl, 2, &rows, notebook);
    }

    gtk_widget_show(tbl);
}

void close_gnuplot_dialog (GtkWidget *w, gpointer p)
{
    free(linetitle);
    free(lineformula);
    free(stylecombo);
    free(yaxiscombo);
    free(linescale);
    free(linewidth);

    linetitle = NULL;
    lineformula = NULL;
    stylecombo = NULL;
    yaxiscombo = NULL;
    linescale = NULL;
    linewidth = NULL;

    gtk_widget_destroy(GTK_WIDGET(p));
}

static int gpt_expand_widgets (int n)
{
    int i;

    if (n == 0) {
	return 0;
    }

    linetitle = myrealloc(linetitle, n * sizeof *linetitle);
    lineformula = myrealloc(lineformula, n * sizeof *lineformula);
    stylecombo = myrealloc(stylecombo, n * sizeof *stylecombo);
    yaxiscombo = myrealloc(yaxiscombo, n * sizeof *yaxiscombo);
    linescale = myrealloc(linescale, n * sizeof *linescale);
    linewidth = myrealloc(linewidth, n * sizeof *linewidth);
    
    if (linetitle == NULL || lineformula == NULL ||
	stylecombo == NULL || yaxiscombo == NULL || 
	linescale == NULL || linewidth == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	linetitle[i] = NULL;
	lineformula[i] = NULL;
	stylecombo[i] = NULL;
	yaxiscombo[i] = NULL;
	linescale[i] = NULL;
	linewidth[i] = NULL;
    }

    gui_nlines = n;

    return 0;
}

static int gpt_allocate_widgets (int n)
{
    int i;

    if (n == 0) {
	return 0;
    }

    linetitle = malloc(n * sizeof *linetitle);
    lineformula = malloc(n * sizeof *lineformula);
    stylecombo = malloc(n * sizeof *stylecombo);
    yaxiscombo = malloc(n * sizeof *yaxiscombo);
    linescale = malloc(n * sizeof *linescale);
    linewidth = malloc(n * sizeof *linewidth);
    
    if (linetitle == NULL || lineformula == NULL ||
	stylecombo == NULL || yaxiscombo == NULL || 
	linescale == NULL || linewidth == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	linetitle[i] = NULL;
	lineformula[i] = NULL;
	stylecombo[i] = NULL;
	yaxiscombo[i] = NULL;
	linescale[i] = NULL;
	linewidth[i] = NULL;
    }

    gui_nlines = n;

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

    g_object_set_data(G_OBJECT(gpt_control), "notebook", notebook);

    gpt_tab_main(notebook, spec);
    gpt_tab_XY(notebook, spec, 0);
    gpt_tab_XY(notebook, spec, 1);

    if (spec->flags & GPT_Y2AXIS) {
	gpt_tab_XY(notebook, spec, 2);
    }

    if (spec->lines != NULL) {
	gpt_tab_lines(notebook, spec, 0);
    }

    gpt_tab_labels(notebook, spec);

    if (!frequency_plot_code(spec->code)) {
	gpt_tab_palette(notebook);
    }

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

