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

#ifndef GRETL_FOREIGN_H
#define GRETL_FOREIGN_H

typedef enum {
    LANG_R = 1,
    LANG_OX,
    LANG_OCTAVE,
    LANG_STATA,
    LANG_PYTHON,
    LANG_JULIA,
    LANG_MPI
} ForeignLang;

int foreign_start (int ci, const char *param, gretlopt opt, PRN *prn);

int foreign_append (const char *line, int context);

int foreign_execute (DATASET *dset, gretlopt opt, PRN *prn);

void foreign_destroy (void);

int check_set_R_home (void);

int execute_R_buffer (const char *buf,
		      const DATASET *dset,
		      gretlopt opt,
		      PRN *prn);

int check_R_depends (const char *pkgname,
		     const char *deps,
		     PRN *inprn);

int write_gretl_R_files (const char *buf,
			 const DATASET *dset,
			 gretlopt opt);

void delete_gretl_R_files (void);

int write_gretl_foreign_script (const char *buf, int lang,
				gretlopt opt, const DATASET *dset,
				const char **pfname);

int gretl_max_mpi_processes (void);

#ifdef HAVE_MPI

void set_mpi_variant (const char *pref);

int check_for_mpiexec (void);

int get_default_mpi_np (void);

#endif /* HAVE_MPI */

#ifdef USE_RLIB

int get_R_function_by_name (const char *name);

int gretl_R_get_call (const char *name, int argc);

int gretl_R_function_add_arg (void *ptr, GretlType type);

int gretl_R_function_add_series (double *x, const DATASET *dset, int v);

int gretl_R_function_exec (const char *name, int *rtype, void **ret);

void gretl_R_reset_error (void);

void gretl_R_cleanup (void);

#endif /* USE_RLIB */

#endif /* GRETL_FOREIGN_H */
