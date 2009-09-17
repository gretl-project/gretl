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

#ifndef GRAPHICS_H_
#define GRAPHICS_H_

void pdf_ps_dialog (GPT_SPEC *spec, GtkWidget *parent);

void save_graphic_to_file (gpointer data, const char *fname);

GPT_SPEC *graph_saver_get_plotspec (gpointer p);

void pdf_saver_set_fontname (gpointer p, const char *fontname);

const char *pdf_saver_current_font (gpointer p);

#endif /* GRAPHICS_H_ */
