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

#ifndef GRETL_LPSOLVE_H
#define GRETL_LPSOLVE_H

/* Definitions we need for the gretl lpsolve plugin */

typedef struct _lprec lprec;
typedef double REAL;

enum { LE = 1, GE, EQ };

/* reporting levels */
#define NEUTRAL    0
#define CRITICAL   1
#define SEVERE     2
#define IMPORTANT  3
#define NORMAL     4
#define DETAILED   5
#define FULL       6

/* solver status values */
#define UNKNOWNERROR  -5
#define DATAIGNORED   -4
#define NOBFP         -3
#define NOMEMORY      -2
#define NOTRUN        -1
#define OPTIMAL        0
#define SUBOPTIMAL     1
#define INFEASIBLE     2
#define UNBOUNDED      3
#define DEGENERATE     4
#define NUMFAILURE     5
#define USERABORT      6
#define TIMEOUT        7
#define RUNNING        8
#define PRESOLVED      9

#ifdef PRELINKED

/* Declarations we need for the gretl lpsolve plugin */

lprec *make_lp (int rows, int columns);
lprec *read_lp (FILE *fp, int verbose, char *lp_name);
unsigned char set_add_rowmode (lprec *lp, unsigned char s);
void delete_lp (lprec *lp);
void set_verbose (lprec *lp, int verbose);
void set_maxim (lprec *lp);
void set_minim (lprec *lp);
unsigned char set_lp_name (lprec *lp, char *s);
unsigned char set_obj_fn (lprec *lp, REAL *row);
unsigned char add_constraint (lprec *lp, REAL *row,
			      int constr_type, REAL rh);
unsigned char set_col_name (lprec *lp, int col, char *name);
unsigned char set_row_name (lprec *lp, int row, char *name);
char *get_col_name (lprec *lp, int col);
char *get_row_name (lprec *lp, int row);
unsigned char set_int (lprec *lp, int col, unsigned char s);
int get_Nrows (lprec *lp);
int get_Ncolumns (lprec *lp);
int solve (lprec *lp);
void print_objective (lprec *lp);
void print_solution (lprec *lp, int columns);
void print_constraints (lprec *lp, int columns);
void print_duals (lprec *lp);
REAL get_objective (lprec *lp);
REAL get_accuracy (lprec *lp);
unsigned char get_primal_solution (lprec *lp, REAL *pv);
unsigned char get_dual_solution (lprec *lp, REAL *duals);
unsigned char get_sensitivity_rhs (lprec *lp, REAL *duals,
				   REAL *from, REAL *till);
void set_outputstream (lprec *lp, FILE *fp);

#endif /* PRELINKED */

#endif /* GRETL_LPSOLVE_H */
