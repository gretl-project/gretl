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

int read_plotfile (const char *fname, GPT_SPEC *spec);

void do_save_graph (const char *fname, char *savestr);

void gpt_save_dialog (void);

#ifdef GNUPLOT_PNG
void save_this_graph (GPT_SPEC *plot, const char *fname);

int gnuplot_show_png (const char *plotfile, GPT_SPEC *spec, int saved);
#endif /* GNUPLOT_PNG */

#endif /* GPT_CONTROL_H */
