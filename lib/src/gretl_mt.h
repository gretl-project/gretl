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

#ifndef GRETL_MT_H
#define GRETL_MT_H

int gretl_n_processors (void);

int gretl_n_physical_cores (void);

void num_threads_init (int blas_type);

int get_omp_n_threads (void);

int set_omp_n_threads (int n);

#ifdef _OPENMP

int gretl_use_openmp (guint64 n);

int openmp_by_default (void);

int set_omp_mnk_min (int n);

#endif

#endif /* GRETL_MT_H */
