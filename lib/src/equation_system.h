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

struct _gretl_equation_system {
    int type;
    int n_equations;
    int **lists;
};

gretl_equation_system *parse_system_start_line (const char *line);

int gretl_equation_system_expand (gretl_equation_system *sys, 
				  int *list);

int gretl_equation_system_print (gretl_equation_system *sys, PRN *prn);

void gretl_equation_system_destroy (gretl_equation_system *sys);


#endif /* GRETL_EQUATION_SYSTEM_H */
