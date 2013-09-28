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

#ifndef HELPFILES_H
#define HELPFILES_H

void helpfile_init (void);

void show_gui_help (int helpcode);

GtkWidget *context_help_button (GtkWidget *hbox, int helpcode);

void command_help_callback (int cmdnum, int en);

void function_help_callback (int fnum);

void plain_text_cmdref (GtkAction *action);

void genr_funcs_ref (GtkAction *action);

gint interactive_script_help (GtkWidget *widget, GdkEventButton *b,
			      windata_t *vwin);

void gretl_show_pdf (const char *fname);

void display_pdf_help (GtkAction *action);

void display_gnuplot_help (void);

void display_x12a_help (void);

void listbox_find (gpointer unused, gpointer data);

void text_find (gpointer unused, gpointer data);

void vwin_add_finder (windata_t *vwin);

int add_help_navigator (windata_t *vwin, GtkWidget *hp);

char *quoted_help_string (const char *s);

int function_help_index_from_word (const char *s);

int gui_console_help (const char *param);

void set_up_helpview_menu (windata_t *hwin);

void notify_string_not_found (GtkWidget *entry);

#endif /* HELPFILES_H */
