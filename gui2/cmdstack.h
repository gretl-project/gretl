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

#ifndef CMDSTACK_H
#define CMDSTACK_H

void free_command_stack (void);

int add_command_to_stack (const char *s);

int add_model_command_to_stack (const char *s, int model_ID);

gchar *get_logfile_content (int *err);

void view_command_log (void);

void set_session_log (const char *dirname, int code);

void maybe_suspend_session_log (void);

#endif /* CMDSTACK_H */
