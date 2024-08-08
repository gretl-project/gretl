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

#if defined(__linux) || defined(linux)
# define GRETL_PID_FILE
#elif defined(WIN32)
# define GRETL_PID_FILE
#elif defined(HAVE_LIBPROC_H) && defined(HAVE_SYS_PROC_INFO_H)
# define GRETL_PID_FILE
#endif

#if defined(__linux) || defined(linux) || defined(WIN32)
# define GRETL_OPEN_HANDLER
#endif

#define IPC_DEBUG 0

#if IPC_DEBUG
extern FILE *fipc;
#endif

#ifdef GRETL_PID_FILE

int write_pid_to_file (void);

void delete_pid_from_file (void);

int gretl_sequence_number (void);

#endif

#ifdef GRETL_OPEN_HANDLER

long gretl_prior_instance (void);

int install_open_handler (void);

gboolean forward_open_request (long gpid, const char *fname);

void record_gretl_binary_path (const char *argv0);

gchar *get_gretl_binary_path (void);

#endif

#endif /* GRETL_IPC_H */

