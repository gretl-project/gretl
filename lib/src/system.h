/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#ifndef GRETL_EQUATION_SYSTEM_H
#define GRETL_EQUATION_SYSTEM_H

typedef struct _gretl_equation_system gretl_equation_system;

enum {
    SUR = 0,
    THREESLS
} gretl_system_types;

enum {
    GRETL_SYSTEM_SAVE_UHAT = 1 << 0,
    GRETL_SYSTEM_SAVE_YHAT = 1 << 1
};

gretl_equation_system *parse_system_start_line (const char *line);

int gretl_equation_system_append (gretl_equation_system *sys, 
				  int *list);


int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn);

void gretl_equation_system_destroy (gretl_equation_system *sys);

int system_save_uhat (const gretl_equation_system *sys);

int system_save_yhat (const gretl_equation_system *sys);

int system_n_equations (const gretl_equation_system *sys);

int system_n_indep_vars (const gretl_equation_system *sys);

int system_adjust_t1t2 (const gretl_equation_system *sys,
			int *t1, int *t2, const double **Z);

int *system_get_list (const gretl_equation_system *sys, int i);

int system_get_depvar (const gretl_equation_system *sys, int i);

const char *gretl_system_short_string (const MODEL *pmod);

int system_get_type (const gretl_equation_system *sys);

#endif /* GRETL_EQUATION_SYSTEM_H */
