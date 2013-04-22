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

#ifndef GRETL_IPC_H
#define GRETL_IPC_H

#if defined(__linux) || defined(linux) || defined(WIN32)
# define GRETL_OPEN_HANDLER
#endif

int write_pid_to_file (void);

int delete_pid_from_file (void);

long gretl_prior_instance (void);

int gretl_sequence_number (void);

#ifdef GRETL_OPEN_HANDLER

int install_open_handler (void);

gboolean forward_open_request (long gpid, const char *fname);

#endif

#endif /* GRETL_IPC_H */

