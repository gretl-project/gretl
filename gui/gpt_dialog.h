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

#ifndef GPT_DIALOG_H
#define GPT_DIALOG_H

GtkWidget *plot_add_editor (png_plot *plot);

void plot_show_font_selector (png_plot *plot, const char *currfont);

void png_font_selector (GtkButton *button, gpointer p);

void pdf_font_selector (GtkButton *button, gpointer p);

void set_plotbars_filename (const char *fname, gpointer data);

void show_color_tool (void);

#endif /* GPT_DIALOG_H */
