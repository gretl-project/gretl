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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* gpt_control.c for gretl -- GUI gnuplot controller */

#include "gretl.h"

#ifdef G_OS_WIN32
# include <windows.h>
#endif

struct gpt_titles_t {
    char *description; /* How the field will show up in the options dialog */
    short tab;         /* which tab (if any) does the item fall under? */
    GtkWidget *widget;
};

typedef struct {
    gint ID;
    GtkWidget *isauto;
    GtkWidget *min;
    GtkWidget *max;
} GPT_RANGE;

GtkWidget *linetitle[6];
GtkWidget *stylecombo[6];
GtkWidget *yaxiscombo[6];
GtkWidget *linescale[6];
static GtkWidget *gpt_control;
static GtkWidget *keycombo;
static GtkWidget *termcombo;
GtkWidget *filesavebutton;

GPT_RANGE axis_range[3];

#define NTITLES 4

struct gpt_titles_t gpt_titles[] = {
    {"Title of plot", 0, NULL},
    {"Title for axis", 1, NULL},
    {"Title for axis", 2, NULL},
    {"Title for axis", 3, NULL},
};

/* ........................................................... */

static void close_plot (GtkWidget *widget, gpointer data) 
{
    GPT_SPEC *plot = (GPT_SPEC *) data;

    gpt_control = NULL;
#ifdef G_OS_WIN32
    fclose(plot->fp);
#else
    pclose(plot->fp);
#endif
    free_plot(plot);
}    

/* ........................................................... */

static void flip_manual_range (GtkWidget *widget, gpointer data)
{
    gint axis = GPOINTER_TO_INT(data);

    if (GTK_TOGGLE_BUTTON (axis_range[axis].isauto)->active) {
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].min), FALSE);
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].max), FALSE);
    } else {
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].min), TRUE);
	gtk_widget_set_sensitive(GTK_WIDGET(axis_range[axis].max), TRUE);
    }
}

/* ........................................................... */

static void widget_to_str (GtkWidget *w, char *str, size_t n)
{
    char *tmp;

    str[0] = '\0';
    tmp = gtk_entry_get_text(GTK_ENTRY(w));
    if (tmp != NULL && strlen(tmp))
	safecpy(str, tmp, n-1);
}

/* ........................................................... */

static void apply_gpt_changes (GtkWidget *widget, gpointer data) 
{
    gchar *yaxis;
    int i, save = 0, k, numlines;
    GPT_SPEC *plot = (GPT_SPEC *) data;

#ifdef G_OS_WIN32
    plot->fp = fopen(paths.plotfile, "w");
    if (plot->fp == NULL) {
	errbox("Couldn't open plot file");
	return;
    }
#endif

    numlines = plot->list[0] - 1;
    if (widget == filesavebutton) {
	widget_to_str(GTK_COMBO(termcombo)->entry, plot->termtype, 
		      sizeof plot->termtype);
	if (strcmp(plot->termtype, "screen")) save = 1;
    }
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].widget != NULL) {
	    widget_to_str(gpt_titles[i].widget, plot->titles[i], 
			  sizeof plot->titles[0]);
	}
    }

    widget_to_str(GTK_COMBO(keycombo)->entry, plot->keyspec, 
		  sizeof plot->keyspec);

    plot->y2axis = 0;
    for (i=0; i<numlines; i++) {
	plot->lines[i].yaxis = 1;
	yaxis = 
	    gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(yaxiscombo[i])->entry));
	if (yaxis != NULL && strlen(yaxis) && !strcmp(yaxis, "right"))
	    plot->lines[i].yaxis = 2;	
	if (plot->lines[i].yaxis == 2) plot->y2axis = 1;
    }
    k = (plot->y2axis)? 3 : 2;
    for (i=0; i<k; i++) {
	if (axis_range[i].isauto != NULL) {
	    if (GTK_TOGGLE_BUTTON (axis_range[i].isauto)->active) {
		strcpy(plot->range[i][0], "*");
		strcpy(plot->range[i][1], "*");
	    } else {
		widget_to_str(axis_range[i].min, plot->range[i][0], 
			      sizeof plot->range[0][0]);
		widget_to_str(axis_range[i].max, plot->range[i][1], 
			      sizeof plot->range[0][1]);
	    }
	}
    }
    for (i=0; i<numlines; i++) {
	widget_to_str(GTK_COMBO(stylecombo[i])->entry, 
		      plot->lines[i].style, sizeof plot->lines[0].style);
	widget_to_str(linetitle[i], 
		      plot->lines[i].title, sizeof plot->lines[0].title);
	widget_to_str(linescale[i], plot->lines[i].scale, 
		      sizeof plot->lines[0].scale);
    }

    if (plot->edit == 2 || plot->edit == 3) {  /* silent update */
	plot->edit -= 2;
	return;
    }

    if (save) { /* do something other than a screen graph? */
#ifdef G_OS_WIN32
	fclose(plot->fp);
#endif
	file_selector("Save gnuplot graph", paths.userdir,
		      SAVE_GNUPLOT, plot);
    }
    else {
	go_gnuplot(plot, NULL, &paths);
    } 
}

/* ........................................................... */

static void save_session_graph (GtkWidget *w, gpointer data)
{
    GPT_SPEC *plot = (GPT_SPEC *) data;
    int err = 0;

    plot->edit += 2;
    apply_gpt_changes(NULL, plot);
    strcpy(plot->termtype, "plot commands");
    plot->edit += 2;
    err = go_gnuplot(plot, plot->fname, &paths);
    if (err == 1) errbox("Error saving graph");
    else infobox("graph saved");
}

/* ........................................................... */

static void gpt_tab_main (GtkWidget *notebook, GPT_SPEC *plot) 
{
    GtkWidget *tempwid, *box, *tbl;
    int i, tbl_len;
    GList *keypos = NULL;

    char *key_positions[6] = {
	"left top",
	"right top",
	"left bottom",
	"right bottom",
	"outside",
	"none"
    };

    for (i=0; i<6; i++)
	keypos = g_list_append(keypos, key_positions[i]);
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show(box);
    
    tempwid = gtk_label_new ("Main");
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (tbl), 5);
    gtk_box_pack_start (GTK_BOX (box), tbl, FALSE, FALSE, 0);
    gtk_widget_show (tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 0) {
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
	    tempwid = gtk_label_new(gpt_titles[i].description);
	    gtk_misc_set_alignment(GTK_MISC (tempwid), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE (tbl), 
				      tempwid, 0, 1, tbl_len-1, tbl_len);
	    gtk_widget_show(tempwid);
	    tempwid = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      tempwid, 1, 2, tbl_len-1, tbl_len);
	    gtk_entry_set_text(GTK_ENTRY(tempwid), plot->titles[i]);
	    gtk_signal_connect(GTK_OBJECT(tempwid), "activate", 
			       GTK_SIGNAL_FUNC(apply_gpt_changes), 
			       plot);
	    gtk_widget_show (tempwid);
	    gpt_titles[i].widget = tempwid;
	}
    }
    tbl_len++;
    tempwid = gtk_label_new("key position");
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(tempwid);

    keycombo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      keycombo, 1, 2, tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(keycombo), keypos); 
    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(keycombo)->entry), plot->keyspec);
    gtk_widget_show (keycombo);	
}

/* ........................................................... */

static void gpt_tab_output (GtkWidget *notebook, GPT_SPEC *plot) 
{
    GtkWidget *tempwid, *box, *tbl;
    int i, tbl_len;
    GList *termtype = NULL;

    char *terminal_types[5] = {
	"postscript",
	"fig",
	"latex",
	"png",
	"plot commands"
    };  

    for (i=0; i<5; i++)
	termtype = g_list_append(termtype, terminal_types[i]);
   
    box = gtk_vbox_new (FALSE, 0);
    gtk_container_border_width (GTK_CONTAINER (box), 10);
    gtk_widget_show(box);
    
    tempwid = gtk_label_new ("Output to file");
    gtk_widget_show (tempwid);
    gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new (tbl_len, 2, FALSE);
    gtk_table_set_row_spacings (GTK_TABLE (tbl), 5);
    gtk_table_set_col_spacings (GTK_TABLE (tbl), 5);
    gtk_box_pack_start (GTK_BOX (box), tbl, FALSE, FALSE, 0);
    gtk_widget_show (tbl);
   
    tbl_len++;
    tempwid = gtk_label_new("output format");
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(tempwid);

    termcombo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      termcombo, 1, 2, tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(termcombo), termtype);   
    gtk_widget_show(termcombo);

    /* button to generate output to file */
    filesavebutton = gtk_button_new_with_label ("Save to file...");
    GTK_WIDGET_SET_FLAGS(filesavebutton, GTK_CAN_DEFAULT);
    tbl_len++;
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      filesavebutton, 1, 2, tbl_len-1, tbl_len);
    gtk_signal_connect (GTK_OBJECT(filesavebutton), "clicked", 
                        GTK_SIGNAL_FUNC(apply_gpt_changes), 
			plot);
    gtk_widget_grab_default(filesavebutton);
    gtk_widget_show(filesavebutton);    
}

/* ........................................................... */

static void gpt_tab_lines (GtkWidget *notebook, GPT_SPEC *plot) 
{
    GtkWidget *tempwid, *box, *tbl;
    int i, tbl_len, tbl_num, tbl_col, numlines;
    char label_text[32];
    GList *plot_types = NULL;
    GList *yaxis_loc = NULL;

    numlines = plot->list[0] - 1;
    if (plot->ts) {
	plot_types = g_list_append(plot_types, "lines");
	plot_types = g_list_append(plot_types, "points");
    } else {
	plot_types = g_list_append(plot_types, "points");
	plot_types = g_list_append(plot_types, "lines");
    }
    plot_types = g_list_append(plot_types, "linespoints"); 
    plot_types = g_list_append(plot_types, "impulses");
    plot_types = g_list_append(plot_types, "dots");

    yaxis_loc = g_list_append(yaxis_loc, "left");
    yaxis_loc = g_list_append(yaxis_loc, "right");

    box = gtk_vbox_new(FALSE, 0);
    gtk_container_border_width(GTK_CONTAINER (box), 10);
    gtk_widget_show(box);

    tempwid = gtk_label_new("Lines");

    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 3, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    tbl_num = tbl_col = 0;

    for (i=0; i<numlines; i++) {
	/* identifier and key or legend text */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	sprintf(label_text, "line %d: ", i + 1);
	tempwid = gtk_label_new(label_text);
	gtk_misc_set_alignment(GTK_MISC(tempwid), 1, 0.5);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				   tempwid, 0, 1, tbl_len-1, tbl_len);
	gtk_widget_show(tempwid);

	tempwid = gtk_label_new("key");
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  tempwid, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(tempwid);

	linetitle[i] = gtk_entry_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linetitle[i], 2, 3, tbl_len-1, tbl_len);
	gtk_entry_set_text (GTK_ENTRY(linetitle[i]), plot->lines[i].title);
	gtk_signal_connect (GTK_OBJECT(linetitle[i]), "activate", 
			    GTK_SIGNAL_FUNC(apply_gpt_changes), 
			    plot);
	gtk_widget_show(linetitle[i]);
	/* line type or style */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	tempwid = gtk_label_new("type");
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  tempwid, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(tempwid);

	stylecombo[i] = gtk_combo_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  stylecombo[i], 2, 3, tbl_len-1, tbl_len);
	gtk_combo_set_popdown_strings(GTK_COMBO(stylecombo[i]), plot_types); 
	gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(stylecombo[i])->entry), 
			   plot->lines[i].style);  
	gtk_widget_show(stylecombo[i]);	
	/* scale factor for data */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	tempwid = gtk_label_new("scale");
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  tempwid, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show (tempwid);

	linescale[i] = gtk_entry_new_with_max_length (6);
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  linescale[i], 2, 3, tbl_len-1, tbl_len);
	gtk_entry_set_text(GTK_ENTRY(linescale[i]), plot->lines[i].scale);
	gtk_signal_connect(GTK_OBJECT(linescale[i]), "activate", 
			   GTK_SIGNAL_FUNC(apply_gpt_changes), 
			   plot);
	gtk_widget_show(linescale[i]);
	/* use left or right y axis? */
	tbl_len++;
	gtk_table_resize(GTK_TABLE(tbl), tbl_len, 3);
	tempwid = gtk_label_new("y axis");
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  tempwid, 1, 2, tbl_len-1, tbl_len);
	gtk_widget_show(tempwid);

	yaxiscombo[i] = gtk_combo_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), 
				  yaxiscombo[i], 2, 3, tbl_len-1, tbl_len);
	gtk_combo_set_popdown_strings(GTK_COMBO(yaxiscombo[i]), yaxis_loc); 
	gtk_entry_set_text (GTK_ENTRY(GTK_COMBO(yaxiscombo[i])->entry), 
			    (plot->lines[i].yaxis == 1)? "left": "right");  
	gtk_widget_show (yaxiscombo[i]);	
    }
}

/* ........................................................... */

static void gpt_tab_XY (GtkWidget *notebook, GPT_SPEC *plot, gint axis) 
{
    GtkWidget *box, *manual, *tbl, *tempwid = NULL;
    int i, tbl_len;
   
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show (box);

    if (axis == 0)
	tempwid = gtk_label_new("X-axis");
    else if (axis == 1)
	tempwid = gtk_label_new("Y-axis");
    else if (axis == 2)
	tempwid = gtk_label_new("Y2-axis");

    gtk_widget_show(tempwid);
    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tempwid);   

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);
   
    for (i=0; i<NTITLES; i++) {
	if (gpt_titles[i].tab == 1 + axis) {
	    tbl_len++;
	    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
            
	    tempwid = gtk_label_new(gpt_titles[i].description);
	    gtk_misc_set_alignment(GTK_MISC(tempwid), 1, 0.5);
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      tempwid, 0, 1, tbl_len-1, tbl_len);
	    gtk_widget_show(tempwid);

	    tempwid = gtk_entry_new();
	    gtk_table_attach_defaults(GTK_TABLE(tbl), 
				      tempwid, 1, 2, tbl_len-1, tbl_len);
	    gtk_entry_set_text(GTK_ENTRY(tempwid), plot->titles[i]);
	    gtk_signal_connect(GTK_OBJECT (tempwid), "activate", 
			       GTK_SIGNAL_FUNC(apply_gpt_changes), 
			       plot);
	    gtk_widget_show(tempwid);
	    gpt_titles[i].widget = tempwid;
	}
    }    

    /* axis range: auto versus manual buttons */
    axis_range[axis].ID = axis;
    tbl_len +=3;
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    tempwid = gtk_label_new("");
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 1, tbl_len-3, tbl_len-2);
    gtk_widget_show(tempwid);
    axis_range[axis].isauto = 
	gtk_radio_button_new_with_label(NULL, "auto axis range");
    gtk_signal_connect(GTK_OBJECT(axis_range[axis].isauto), "clicked",
		       GTK_SIGNAL_FUNC(flip_manual_range), 
		       GINT_TO_POINTER(axis_range[axis].ID));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON
				 (axis_range[axis].isauto), TRUE);
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      axis_range[axis].isauto, 
			      0, 1, tbl_len-2, tbl_len-1);
    gtk_widget_show(axis_range[axis].isauto);
    manual = 
	gtk_radio_button_new_with_label(gtk_radio_button_group 
					(GTK_RADIO_BUTTON 
					 (axis_range[axis].isauto)),
					"manual range:"); 
    gtk_signal_connect(GTK_OBJECT(manual), "clicked",
		       GTK_SIGNAL_FUNC(flip_manual_range), 
		       GINT_TO_POINTER(axis_range[axis].ID));
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      manual, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(manual);

    /* axis range min. entry */
    tbl_len++;
    tempwid = gtk_label_new("minimum");
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(tempwid);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    axis_range[axis].min = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      axis_range[axis].min, 1, 2, tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(axis_range[axis].min), "");
    gtk_signal_connect(GTK_OBJECT(axis_range[axis].min), "activate", 
		       GTK_SIGNAL_FUNC(apply_gpt_changes), 
		       plot);
    gtk_widget_show(axis_range[axis].min);

    /* axis range max. entry */
    tbl_len++;
    tempwid = gtk_label_new("maximum");
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      tempwid, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(tempwid);
    gtk_table_resize(GTK_TABLE(tbl), tbl_len, 2);
    axis_range[axis].max = gtk_entry_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      axis_range[axis].max, 1, 2, tbl_len-1, tbl_len);
    gtk_entry_set_text(GTK_ENTRY(axis_range[axis].max), "");
    gtk_signal_connect(GTK_OBJECT(axis_range[axis].max), "activate", 
		       GTK_SIGNAL_FUNC(apply_gpt_changes), 
		       plot);
    gtk_widget_show(axis_range[axis].max);
   
    if (strcmp(plot->range[axis][0], "*") == 0)
	flip_manual_range(NULL, GINT_TO_POINTER(axis_range[axis].ID));
    else {
	gtk_entry_set_text(GTK_ENTRY(axis_range[axis].min),
			   plot->range[axis][0]);
	gtk_entry_set_text(GTK_ENTRY(axis_range[axis].max),
			   plot->range[axis][1]);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON
				     (axis_range[axis].isauto), FALSE);
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON
				     (manual), TRUE);
    }
}

/* ........................................................... */

void gnuplot_dialog (GPT_SPEC *plot) 
{
    GtkWidget *tempwid, *notebook;
    int i;

    if (gpt_control != NULL) {
	errbox("You can only have one plot controller open\n"
	       "at any given time");
	return;
    }

    for (i=0; i<3; i++) axis_range[i].isauto = NULL;

    gpt_control = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(gpt_control), "gretl plot controls");
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(gpt_control)->vbox), 10);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(gpt_control)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(gpt_control), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (gpt_control), "destroy",
                        GTK_SIGNAL_FUNC (close_plot), (gpointer *) plot);
   
    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->vbox), 
		       notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    gpt_tab_main(notebook, plot);
    gpt_tab_XY(notebook, plot, 0);
    gpt_tab_XY(notebook, plot, 1);
    if (plot->y2axis) gpt_tab_XY(notebook, plot, 2);
    gpt_tab_lines(notebook, plot); 
    gpt_tab_output(notebook, plot);

    tempwid = gtk_button_new_with_label ("Redraw");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(apply_gpt_changes), plot);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label ("Save");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT(tempwid), "clicked", 
                        GTK_SIGNAL_FUNC(save_session_graph), plot);
    gtk_widget_show (tempwid);

    tempwid = gtk_button_new_with_label("Close");
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
		       tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
                       GTK_SIGNAL_FUNC(delete_widget), gpt_control);
    gtk_widget_show(tempwid);

    tempwid = gtk_button_new_with_label("Help");
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(gpt_control)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC(context_help), 
			GINT_TO_POINTER(GNUPLOT));
    gtk_widget_show (tempwid);

    gtk_widget_show (gpt_control);
}

/* Below: functions for saving last auto-generated graph */

void do_save_graph (const char *fname, char *savestr)
{
    FILE *fq;
    print_t *prn;
    char plottmp[MAXLEN], plotline[MAXLEN], plotcmd[MAXLEN];
    char termstr[MAXLEN];
    int cmds;

    if (!user_fopen("gptout.tmp", plottmp, &prn)) return;
    fq = fopen(paths.plotfile, "r");
    if (fq == NULL) {
	errbox("Couldn't access graph info");
	gretl_print_destroy(prn);
	return;
    } 
    cmds = termtype_to_termstr(savestr, termstr);  
    if (cmds) {
	if (copyfile(paths.plotfile, fname)) 
	    errbox("Failed to copy graph file");
	return;
    } else {
	pprintf(prn, "set term %s\n", termstr);
	pprintf(prn, "set output '%s'\n", fname);
	while (fgets(plotline, MAXLEN-1, fq))
	    pprintf(prn, "%s", plotline);
    }
    gretl_print_destroy(prn);
    fclose(fq);
    sprintf(plotcmd, "\"%s\" \"%s\"", paths.gnuplot, plottmp);
#ifdef G_OS_WIN32
    if (WinExec(plotcmd, SW_SHOWMINIMIZED) < 32)
	errbox("Gnuplot error creating graph");
#else
    if (system(plotcmd))
	errbox("Gnuplot error creating graph");
#endif
    remove(plottmp);
    infobox("Graph saved");
}

/* ........................................................... */

static void plot_save_filesel (GtkWidget *w, gpointer data)
{
    static char savestr[MAXLEN];
    GtkWidget *combo = (GtkWidget *) data;

    strcpy(savestr, 
	   gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(combo)->entry)));
    gtk_widget_destroy(GTK_WIDGET(combo->parent->parent->parent));
    file_selector("save graph", paths.userdir, SAVE_LAST_GRAPH, savestr);
}

/* ........................................................... */

static void kill_gpt_save_dialog (GtkWidget *w, gpointer data)
{
    GtkWidget **dialog = (GtkWidget **) data;

    gtk_widget_destroy(*dialog);
    *dialog = NULL;
}

/* ........................................................... */

void gpt_save_dialog (void)
{
    GtkWidget *tempwid, *tbl, *combo;
    static GtkWidget *dialog;
    gint tbl_len;
    GList *termtype = NULL;
    int i;
    char *ttypes[] = {"postscript","fig","latex","png",
		      "plot commands"};

    if (dialog != NULL) {
	gdk_window_raise(dialog->window);
	return;
    }

    for (i=0; i<5; i++)
	termtype = g_list_append(termtype, ttypes[i]);

    dialog = gtk_dialog_new();
    gtk_window_set_title(GTK_WINDOW(dialog), "gretl: save graph");
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dialog)->vbox), 10);
    gtk_container_border_width 
        (GTK_CONTAINER(GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 2);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->action_area), 15);
    gtk_box_set_homogeneous(GTK_BOX(GTK_DIALOG(dialog)->action_area), TRUE);
    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);
    gtk_signal_connect(GTK_OBJECT(dialog), "destroy",
		       GTK_SIGNAL_FUNC(kill_gpt_save_dialog), &dialog);

    tbl_len = 1;
    tbl = gtk_table_new(tbl_len, 2, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX (GTK_DIALOG(dialog)->vbox), 
		       tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl);

    tempwid = gtk_label_new("output type");
    gtk_table_attach_defaults(GTK_TABLE (tbl), 
			      tempwid, 0, 1, tbl_len-1, tbl_len);
    gtk_widget_show(tempwid);

    combo = gtk_combo_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), 
			      combo, 1, 2, tbl_len-1, tbl_len);
    gtk_combo_set_popdown_strings(GTK_COMBO(combo), termtype);   
    gtk_widget_show(combo);

    tempwid = gtk_button_new_with_label("Save");
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(plot_save_filesel), combo);
    gtk_widget_grab_default(tempwid);
    gtk_widget_show(tempwid);
   
    tempwid = gtk_button_new_with_label("Cancel");
    GTK_WIDGET_SET_FLAGS(tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->action_area), 
                        tempwid, TRUE, TRUE, 0);
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked", 
		       GTK_SIGNAL_FUNC(delete_widget), dialog);
    gtk_widget_show(tempwid);

    gtk_widget_show(dialog);
}

/* ........................................................... */

static void charsub (int c_old, int c_new, char *str)
{
    size_t i, n = strlen(str);

    for (i=0; i<n; i++) {
	if (str[i] == c_old) {
	    str[i] = c_new;
	    break;
	}
    }
}

/* ........................................................... */

static int chop_comma (char *str)
{
    size_t i, n = strlen(str);

    for (i=n-1; i>0; i--) {
	if (isspace((unsigned char) str[i])) continue;
	if (str[i] == ',') {
	    str[i] = 0;
	    return 1;
	}
    }		
    return 0;
}

/* ........................................................... */

static void get_gpt_data (char *line, double *x, double *y)
{
    if (x != NULL) {
	if (sscanf(line, "%lf %lf", x, y) == 2) return;
	if (sscanf(line, "%lf", x) == 1) {
	    *y = NADBL;
	    return;
	}
	if (sscanf(line, "%*s %lf", y) == 1) {
	    *x = NADBL;
	    return;
	}
	*x = NADBL; *y = NADBL;
	return;
    } else {
	if (sscanf(line, "%*f %lf", y) == 1) return;
	*y = NADBL;
    }
    return;
}

/* ........................................................... */

int read_plotfile (GPT_SPEC *plot, char *fname)
    /* read in plot struct from gnuplot command file.
       This is _not_ a general parser for gnuplot files; it is
       designed specifically for files auto-generated by gretl. */
{
    int i, j, t, n, done;
    char line[MAXLEN], *tmp;
    char set_thing[12], setting[MAXLEN], range[32];
    double *tmpy;
    size_t diff;
    FILE *fp;

    /* initialize the struct */
    if ((plot->lines = mymalloc(6 * sizeof(GPT_LINE))) == NULL)
	return 1;
    strcpy(plot->titles[0], "");
    strcpy(plot->titles[1], "");
    strcpy(plot->titles[2], "");
    strcpy(plot->titles[3], "");
    strcpy(plot->keyspec, "left top");
    strcpy(plot->xtics, "");
    strcpy(plot->mxtics, "");
    for (i=0; i<3; i++) {
	strcpy(plot->range[i][0], "*");
	strcpy(plot->range[i][1], "*");
    }
    plot->literal[0] = NULL;
    plot->code = GNUPLOT;

    /* open the file */
    fp = fopen(fname, "r");
    if (fp == NULL) {
	errbox("Couldn't open graph file");
	return 1;
    }
    
    /* first get the "set" lines */
    i = 0;
    while (fgets(line, MAXLEN - 1, fp)) {
	if (strncmp(line, "# mult", 6) == 0) {
	    errbox("Sorry, can't edit multiple scatterplots");
	    free(plot->lines);
	    return 1;
	}
	if (strncmp(line, "# CUSUM", 7) == 0) {
	    errbox("Sorry, can't edit CUSUM plots");
	    free(plot->lines);
	    return 1;
	}
	if (strncmp(line, "# sampl", 7) == 0) {
	    errbox("Sorry, can't edit sampling distribution plots");
	    free(plot->lines);
	    return 1;
	}
	if (strncmp(line, "# freq", 6) == 0 ||
	    strncmp(line, "# peri", 6) == 0) {
	    /* special cases */
	    if (line[2] == 'f') plot->code = FREQ;
	    else  plot->code = PERGM;
	    for (j=0; j<4; j++) {
		plot->literal[j] = mymalloc(MAXLABEL);
		if (plot->literal[j] == NULL) return 1;
		fgets(plot->literal[j], MAXLABEL - 1, fp);
		top_n_tail(plot->literal[j]);
	    }
	    continue;
	}
	if (strncmp(line, "# forecast", 10) == 0) {
	    plot->code = FCASTERR;
	    continue;
	}
	if (strncmp(line, "set ", 4)) 
	    break;
	sscanf(line + 4, "%s", set_thing);
	if (strcmp(set_thing, "y2tics") == 0) {
	    plot->y2axis = 1;
	    continue;
	}
	n = strlen(set_thing);
	strcpy(setting, line + 4 + n);
	top_n_tail(setting);
	if (strstr(set_thing, "range")) {
	    strcpy(range, setting + 1);
	    n = strlen(range);
	    for (j=0; j<n; j++) {
		if (range[j] == ':') {
		    range[j] = 0;
		    strcpy(plot->range[i][0], range);
		    break;
		}
	    }
	    strcpy(range, strchr(setting, ':') + 1);
	    delchar(']', range);
	    strcpy(plot->range[i][1], range);
	    i++;
	}	  
	else if (strcmp(set_thing, "title") == 0) 
	    strcpy(plot->titles[0], setting);
	else if (strcmp(set_thing, "xlabel") == 0)
	    strcpy(plot->titles[1], setting);
	else if (strcmp(set_thing, "ylabel") == 0)
	    strcpy(plot->titles[2], setting);
	else if (strcmp(set_thing, "y2label") == 0)
	    strcpy(plot->titles[3], setting);
	else if (strcmp(set_thing, "key") == 0)
	    strcpy(plot->keyspec, setting);
	else if (strcmp(set_thing, "xtics") == 0) 
	    safecpy(plot->xtics, setting, 15);
	else if (strcmp(set_thing, "mxtics") == 0) 
	    safecpy(plot->mxtics, setting, 3);
    } /* end of "set" lines */

    for (i=0; i<4; i++)
	delchar('\'', plot->titles[i]);
    if (strlen(plot->keyspec) == 0)
	strcpy(plot->keyspec, "none");

    /* then get the "plot" lines */
    if (strncmp(line, "plot ", 4) ||
	(strlen(line) < 10 && fgets(line, MAXLEN - 1, fp) == NULL)) {	
	errbox("Failed to parse gnuplot file");
	fprintf(stderr, "plotfile line: '%s'\n", line);
	fclose(fp);
	return 1;
    }
    i = 0;
    done = 0;
    while (1) {
	top_n_tail(line);
	if (!chop_comma(line)) done++;
	/* scale, [yaxis,] style */
	tmp = strstr(line, "using");
	if (tmp && tmp[11] == '*') {
            safecpy(plot->lines[i].scale, tmp + 12, 7);
	    charsub(')', '\0', plot->lines[i].scale);
	} else {
	    if (tmp) 
		strcpy(plot->lines[i].scale, "1.0");
	    else {
		strcpy(plot->lines[i].scale, "NA");
		tmp = strstr(line, "axes");
		if (tmp == NULL)
		    tmp = strstr(line, "title");
		diff = tmp - line;
		strncpy(plot->lines[i].formula, line, diff);
		plot->lines[i].formula[diff - 1] = 0;
		/*  printf("formula: '%s'\n", plot->lines[i].formula); */
	    }
	}
	tmp = strstr(line, "axes");
	if (tmp) {
	    if (tmp[8] == '2') plot->lines[i].yaxis = 2;
	    else plot->lines[i].yaxis = 1;
	} else plot->lines[i].yaxis = 1;
	tmp = strstr(line, "title"); 
	if (tmp) {
	    tmp += 7;
	    plot->lines[i].title[0] = '\'';
	    j = 0;
	    while (tmp[j] != '\'') {
		plot->lines[i].title[j] = tmp[j];
		j++;
	    }
	    plot->lines[i].title[j] = 0; 
	}
	tmp = strstr(line, " w ");
	if (tmp) {
	    strcpy(plot->lines[i].style, tmp + 3);
	    delchar(',', plot->lines[i].style);
	} else 
	    strcpy(plot->lines[i].style, "points");
	if (done) break;
	i++;
	fgets(line, MAXLEN - 1, fp);
	if (line == NULL) break;
    }

    /* free any unused lines */
    plot->lines = realloc(plot->lines, (i + 1) * sizeof(GPT_LINE));

    /* finally, get the plot data */
    plot->data = mymalloc(datainfo->n * (i + 2) * sizeof(double));
    tmpy = mymalloc(datainfo->n * sizeof *tmpy);
    if (plot->data == NULL || tmpy == NULL) return 1;
    j = 1;
    t = 0;
    n = 0;
    while (fgets(line, MAXLEN - 1, fp)) {
	if (line[0] == 'e') {
	    n = t;
	    t = 0;
	    j++;
	    continue;
	}
	if (j == 1)  
	    get_gpt_data(line, &(plot->data[t]), &(tmpy[t]));
	else
	    get_gpt_data(line, NULL, &(plot->data[j*n + t]));
	t++;
    }
    for (t=0; t<n; t++)
	plot->data[n + t] = tmpy[t];
    free(tmpy);
    plot->t1 = 0;
    plot->t2 = n - 1;
    
    fclose(fp);
    plot->list[0] = i+2;

    if (open_gnuplot_pipe(&paths, plot)) {
	errbox("gnuplot command failed");
	return 1;
    }
    else {
	strcpy(plot->fname, fname); 
	gnuplot_dialog(plot);
    } 
    
    return 0;
}

