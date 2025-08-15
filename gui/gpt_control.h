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

#ifndef GPT_CONTROL_H
#define GPT_CONTROL_H

#include "plotspec.h"

typedef struct png_plot_t png_plot;

int remove_png_term_from_plot_by_name (const char *fname);

void display_session_graph (const char *path, const char *title,
			    void *ptr);

void gnuplot_view_session_graph (const char *fname);

void register_graph (void);

void adjust_plot_collection (const char *parm);

void reset_collection_count (void);

void saver_preview_graph (GPT_SPEC *spec, char *termstr);

int saver_save_graph (GPT_SPEC *spec, char *termstr, const char *fname);

void plot_position_click (GtkWidget *w, png_plot *plot);

int redisplay_edited_plot (png_plot *plot);

void start_editing_png_plot (png_plot *plot);

void set_plot_has_y2_axis (png_plot *plot, gboolean s);

int plot_is_saved (const png_plot *plot);

int plot_is_mouseable (const png_plot *plot);

GtkWidget *plot_get_shell (png_plot *plot);

double plot_get_xmin (png_plot *plot);

double plot_get_ymin (png_plot *plot);

void plot_get_pixel_dims (png_plot *plot,
			  double *pw,
			  double *ph);

int plot_get_coordinates (png_plot *plot,
			  double *xmin,
			  double *xmax,
			  double *ymin,
			  double *ymax);

GPT_SPEC *plot_get_spec (png_plot *plot);

void revise_distribution_plotspec (png_plot *plot, int d, int df1, int df2);

void activate_plot_font_choice (png_plot *plot, const char *fontname);

int gp_term_code (gpointer p, int action);

void save_graph_to_file (gpointer p, const char *fname);

void filter_gnuplot_file (int mono, const char *inpname,
			  FILE *fpin, FILE *fpout);

void run_gnuplot_script (gchar *buf, windata_t *vwin);

void launch_gnuplot_interactive (void);

void gnuplot_view_3d (const char *plotfile);

int gnuplot_show_map (gretl_bundle *b);

int dump_plot_buffer (const char *buf, const char *fname,
		      int addpause, const char *src);

int get_graph_scale (int i, double *s);

#ifndef G_OS_WIN32

int write_plot_for_copy (int target);

#endif

#endif /* GPT_CONTROL_H */
