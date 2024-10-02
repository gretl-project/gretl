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
#include "plotspec.h"
#include "gpt_control.h"
#include "gpt_dialog.h"
#include "dlgutils.h"
#include "fileselect.h"
#include "graphics.h"

/* default values */
static double pwidth = 5.0;
static double pheight = 3.5;
static char font[64] = "Sans";
static int fontsize = 12;
static double lw_factor = 1.0;
static int mono;

#define CMFAC 2.54 /* centimeters per inch */

struct pdf_ps_saver {
    GtkWidget *dialog;
    GPT_SPEC *spec;
    int mono;
    int stdsize;
    double pwidth;
    double pheight;
    char font[64];
    int fontsize;
    double lw_factor;
    GtkWidget *w_in, *h_in;
    GtkWidget *w_cm, *h_cm;
    GtkWidget *combo;
};

static void map_pdf_ps_dims (png_plot *plot,
			     double *pw, double *ph)
{
    double hr;

    plot_get_pixel_dims(plot, pw, ph);
    hr = *ph / *pw;

    if (hr > 8.5 / 5.5) {
	/* tall and skinny */
	*ph = 8.5;
	*pw = *ph / hr;
    } else {
	*pw = 5.5;
	*ph = *pw * hr;
    }
}

static void set_pdf_ps_dims (struct pdf_ps_saver *s, GPT_SPEC *spec)
{
    png_plot *plot = (png_plot *) spec->ptr;
    PlotType ptype = spec->code;
    double w = pwidth, h = pheight;

    if (spec->code == PLOT_GEOMAP && plot != NULL) {
	map_pdf_ps_dims(plot, &w, &h);
	s->stdsize = 0;
    } else if (spec->flags & GPT_LETTERBOX) {
	/* for time series */
	w = (5.0 * GP_LB_WIDTH) / GP_WIDTH;
	h = (3.5 * GP_LB_HEIGHT) / GP_HEIGHT;
    } else if (spec->flags & GPT_XL) {
	/* large */
	w = (5.0 * GP_XL_WIDTH) / GP_WIDTH;
	h = (3.5 * GP_XL_HEIGHT) / GP_HEIGHT;
	s->stdsize = 0;
    } else if (spec->flags & GPT_XXL) {
	/* extra large */
	w = h = (5.0 * GP_XXL_WIDTH) / GP_WIDTH;
	s->stdsize = 0;
    } else if (ptype == PLOT_ROOTS || ptype == PLOT_QQ) {
	/* square plots */
	w = h;
    } 

    s->pwidth = w;
    s->pheight = h;
}

static void saver_init (struct pdf_ps_saver *s,
			GtkWidget *w,
			GPT_SPEC *spec)
{
    s->dialog = w;
    s->spec = spec;
    s->mono = mono;
    s->stdsize = 1;
    strcpy(s->font, font);
    s->fontsize = fontsize;
    s->lw_factor = lw_factor;

    set_pdf_ps_dims(s, spec);

    if (!s->stdsize) {
	s->fontsize *= 0.8;
    } 
}

static void saver_set_defaults (struct pdf_ps_saver *s)
{
    if (s->stdsize) {
	/* save user's preferred size */
	pwidth = s->pwidth;
	pheight = s->pheight;
    }

    strcpy(font, s->font);
    if (s->stdsize) {
	fontsize = s->fontsize;
    }

    mono = s->mono;
    lw_factor = s->lw_factor;    
}

static void set_dim_callback (GtkSpinButton *b, struct pdf_ps_saver *s)
{
    GtkWidget *w = GTK_WIDGET(b);

    if (w == s->w_in) {
	s->pwidth = gtk_spin_button_get_value(b);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->w_cm), s->pwidth * CMFAC);
    } else if (w == s->w_cm) {
	s->pwidth = gtk_spin_button_get_value(b) / CMFAC;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->w_in), s->pwidth);
    } else if (w == s->h_in) {
	s->pheight = gtk_spin_button_get_value(b);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->h_cm), s->pheight * CMFAC);
    } else if (w == s->h_cm) {
	s->pheight = gtk_spin_button_get_value(b) / CMFAC;
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->h_in), s->pheight);
    }
}

static void set_lw_callback (GtkSpinButton *b, struct pdf_ps_saver *s)
{
    s->lw_factor = gtk_spin_button_get_value(b);
}

static GtkWidget *pdf_ps_size_spin_buttons (struct pdf_ps_saver *s)
{
    GtkWidget *tbl, *label, *b;
    GtkWidget *vbox, *hbox;

    s->w_in = gtk_spin_button_new_with_range(1.5, 10, 0.01);
    s->h_in = gtk_spin_button_new_with_range(1.5, 10, 0.01);

    s->w_cm = gtk_spin_button_new_with_range(1.5 * CMFAC, 10 * CMFAC,
					     0.01);
    s->h_cm = gtk_spin_button_new_with_range(1.5 * CMFAC, 10 * CMFAC,
					     0.01);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->w_in), s->pwidth);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->h_in), s->pheight);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->w_cm), CMFAC * s->pwidth);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s->h_cm), CMFAC * s->pheight);

    g_signal_connect(G_OBJECT(s->w_in), "value-changed",
		     G_CALLBACK(set_dim_callback), s);
    g_signal_connect(G_OBJECT(s->h_in), "value-changed",
		     G_CALLBACK(set_dim_callback), s);

    g_signal_connect(G_OBJECT(s->w_cm), "value-changed",
		     G_CALLBACK(set_dim_callback), s);
    g_signal_connect(G_OBJECT(s->h_cm), "value-changed",
		     G_CALLBACK(set_dim_callback), s);

    tbl = gtk_table_new(3, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);

    label = gtk_label_new(_("width"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 0, 1, 0, 1);
    label = gtk_label_new(_("height"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 1, 2, 0, 1);

    gtk_table_attach_defaults(GTK_TABLE(tbl), s->w_in, 0, 1, 1, 2);
    gtk_table_attach_defaults(GTK_TABLE(tbl), s->h_in, 1, 2, 1, 2);
    label = gtk_label_new(_("inches"));
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 2, 3, 1, 2);

    gtk_table_attach_defaults(GTK_TABLE(tbl), s->w_cm, 0, 1, 2, 3);
    gtk_table_attach_defaults(GTK_TABLE(tbl), s->h_cm, 1, 2, 2, 3);
    label = gtk_label_new(_("cm"));
    gtk_misc_set_alignment(GTK_MISC(label), 0.0, 0.5);
    gtk_table_attach_defaults(GTK_TABLE(tbl), label, 2, 3, 2, 3);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), tbl, TRUE, FALSE, 0);

    vbox = gtk_vbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, FALSE, 0);

    vbox_add_hsep(vbox);

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("line width factor"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    b = gtk_spin_button_new_with_range(0.5, 12.0, 0.1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(b), s->lw_factor);
    g_signal_connect(G_OBJECT(b), "value-changed",
		     G_CALLBACK(set_lw_callback), s);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, FALSE, 0);

    return vbox;
}

static GtkWidget *label_in_hbox (const char *s, int center)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *label = gtk_label_new(s);

    if (center) {
	gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    } else {
	gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    }

    return hbox;
}

void set_color_mode (GtkToggleButton *b, struct pdf_ps_saver *s)
{
    s->mono = button_is_active(b);
}

static void color_mode_selector (struct pdf_ps_saver *s,
				 GtkWidget *vbox)
{
    GSList *group = NULL;
    GtkWidget *b;

    b = gtk_radio_button_new_with_label(NULL, _("color"));
    pack_in_hbox(b, vbox, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b));
    b = gtk_radio_button_new_with_label(group, _("monochrome"));
    g_signal_connect(b, "toggled", G_CALLBACK(set_color_mode), s);
    pack_in_hbox(b, vbox, 0);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), s->mono);
}

const char *pdf_ps_saver_current_font (gpointer p)
{
    static char fontname[68];
    struct pdf_ps_saver *s = p;

    sprintf(fontname, "%s %d", s->font, s->fontsize);
    return fontname;
}

void pdf_ps_saver_set_fontname (gpointer p, const char *fontname)
{
    struct pdf_ps_saver *s = p;
    char name[64];
    int psz = 0;

    *name = '\0';
    split_graph_fontspec(fontname, name, &psz);
    if (*name != '\0' && psz > 1) {
	strcpy(s->font, name);
	s->fontsize = psz;
    }
}

static void 
saver_make_term_string (struct pdf_ps_saver *s, char *termstr)
{
    char fontstr[84];
    char lwstr[32];
    const char *ttype;

    gretl_push_c_numeric_locale();

    if (s->lw_factor != 1.0) {
	sprintf(lwstr, " linewidth %g", s->lw_factor);
    } else {
	*lwstr = '\0';
    }

    if (s->spec->termtype == GP_TERM_PDF) {
	ttype = (s->mono)? "pdfcairo noenhanced dashed" : "pdfcairo noenhanced";
	sprintf(fontstr, "font \"%s,%d\"", s->font, s->fontsize);
    } else {
	ttype = (s->mono)? "epscairo noenhanced dashed" : "epscairo noenhanced";
	sprintf(fontstr, "font \"%s,%d\"", s->font, s->fontsize);
    }

    if (s->mono) {
	s->spec->flags |= GPT_MONO;
    } else {
	s->spec->flags &= ~GPT_MONO;
    }

    sprintf(termstr, "set term %s %s%s size %gin,%gin", ttype, fontstr, 
	    lwstr, s->pwidth, s->pheight);

    gretl_pop_c_numeric_locale();
}

static void preview_callback (GtkWidget *w, struct pdf_ps_saver *s)
{
    char termstr[256];

    saver_make_term_string(s, termstr);
    fprintf(stderr, "termstr: '%s'\n", termstr);
    saver_preview_graph(s->spec, termstr);
}

void save_graphic_to_file (gpointer data, const char *fname)
{
    struct pdf_ps_saver *s = data;
    char termstr[256];
    int err;

    saver_make_term_string(s, termstr);
    err = saver_save_graph(s->spec, termstr, fname);
    if (!err) {
	saver_set_defaults(s);
    }
}

static void 
pdf_ps_save_callback (GtkWidget *w, struct pdf_ps_saver *saver)
{
    file_selector_with_parent(SAVE_GRAPHIC, FSEL_DATA_MISC, 
			      saver, saver->dialog);
}

void pdf_ps_dialog (GPT_SPEC *spec, GtkWidget *parent)
{
    struct pdf_ps_saver saver;
    GtkWidget *dialog;
    GtkWidget *vbox, *hbox;
    GtkWidget *label, *b;
    gchar *title;
    int ps;

    ps = (spec->termtype == GP_TERM_EPS);

    title = g_strdup_printf("gretl: %s", _("save graph"));
    /* note: we need to block to keep 'saver' current */
    dialog = gretl_dialog_new(title, parent, GRETL_DLG_BLOCK);
    g_free(title);

    saver_init(&saver, dialog, spec);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    
    label = label_in_hbox(ps ? _("EPS file") : _("PDF file"), 1);
    gtk_box_pack_start(GTK_BOX(vbox), label, TRUE, TRUE, 5);

    label = label_in_hbox(_("Plot dimensions:"), 0);
    gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, FALSE, 0);
    
    gtk_container_add(GTK_CONTAINER(vbox), pdf_ps_size_spin_buttons(&saver));

    vbox_add_hsep(vbox);

    title = g_strdup_printf(_("font: %s"), 
			    pdf_ps_saver_current_font(&saver));
    hbox = gtk_hbox_new(FALSE, 5);
    b = gtk_button_new_with_label(title);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 5);
    gtk_container_add(GTK_CONTAINER(vbox), hbox);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(pdf_font_selector), 
		     &saver);
    g_free(title);

    vbox_add_hsep(vbox);
    color_mode_selector(&saver, vbox);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    b = gtk_button_new_with_label(_("Preview"));
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(preview_callback), &saver);

    b = gtk_button_new_from_stock(GTK_STOCK_SAVE);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect(G_OBJECT(b), "clicked", 
		     G_CALLBACK(pdf_ps_save_callback), &saver);

    b = gtk_button_new_from_stock(GTK_STOCK_CLOSE);
    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, 0);
    g_signal_connect_swapped(G_OBJECT(b), "clicked", 
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    
    gtk_widget_show_all(dialog);

    /* unset mono flag on exit */
    spec->flags &= ~GPT_MONO;
}

GPT_SPEC *graph_saver_get_plotspec (gpointer p)
{
    struct pdf_ps_saver *saver = p;

    return saver->spec;
}
