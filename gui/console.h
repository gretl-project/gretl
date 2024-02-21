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

#ifndef CONSOLE_H
#define CONSOLE_H

windata_t *gretl_console (void);

void console_record_sample (const DATASET *pdinfo);

int console_sample_changed (const DATASET *pdinfo);

int console_is_busy (void);

void clear_console (GtkWidget *w, windata_t *vwin);

int emulate_console_command (const char *cmdline);

#endif /* CONSOLE_H */
