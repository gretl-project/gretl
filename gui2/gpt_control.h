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

#ifdef GNUPLOT_PIPE
void start_editing_session_graph (const char *fname);
#endif

void mark_plot_as_saved (GPT_SPEC *spec);
int remove_png_term_from_plotfile (const char *fname, GPT_SPEC *spec);
void save_this_graph (GPT_SPEC *spec, const char *fname);
void display_session_graph_png (const char *fname);
int gnuplot_show_png (const char *plotfile, GPT_SPEC *spec, int saved);

#endif /* GPT_CONTROL_H */
