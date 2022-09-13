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

/* monte_carlo.h for gretl: handle command loops and conditionals */

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

typedef struct LOOPSET_ LOOPSET;

int gretl_compiling_loop (void);

void gretl_abort_compiling_loop (void);

int gretl_execute_loop (void);

int ok_in_loop (int ci);

int gretl_loop_append_line (ExecState *s, DATASET *dset);

int gretl_loop_append_line_full (ExecState *s, DATASET *dset,
				 LOOPSET **ploop);

int gretl_loop_exec (ExecState *s, DATASET *dset);

void gretl_loop_destroy (LOOPSET *loop);

int model_is_in_loop (const MODEL *pmod);

int scalar_is_read_only_index (const char *name);

void loop_reset_uvars (LOOPSET *loop);

int get_loop_renaming (void);

#endif /* MONTE_CARLO_H */
