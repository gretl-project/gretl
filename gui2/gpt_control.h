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

#ifndef GPT_CONTROL_H
#define GPT_CONTROL_H

typedef struct png_plot_t png_plot;

void mark_plot_as_saved (GPT_SPEC *spec);
int remove_png_term_from_plotfile (const char *fname, GPT_SPEC *spec);
void save_this_graph (GPT_SPEC *spec, const char *fname);
void display_session_graph_png (char *fname);
int gnuplot_show_png (const char *plotfile, GPT_SPEC *spec, int saved);

void plot_label_position_click (GtkWidget *w, GPT_SPEC *spec);

int redisplay_edited_png (png_plot *plot);
void plot_remove_controller (png_plot *plot);
void set_plot_has_y2_axis (png_plot *plot, gboolean s);
int plot_is_mouseable (const png_plot *plot);
GtkWidget *plot_get_shell (png_plot *plot);


#endif /* GPT_CONTROL_H */
