/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* monte_carlo.h for gretl: handle command loops and conditionals */

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

enum ifcodes {
    SET_FALSE,
    SET_TRUE,
    SET_ELSE,
    SET_ENDIF,
    IS_FALSE,
    RELAX
};

typedef struct LOOPSET_ LOOPSET;

/* functions follow */

int ok_in_loop (int ci, const LOOPSET *loop);

LOOPSET *parse_loopline (char *line, LOOPSET *ploop, int loopstack,
			 DATAINFO *pdinfo, const double **Z);

LOOPSET *gretl_loop_terminate (LOOPSET *loop);

int add_to_loop (LOOPSET *loop, char *line, int ci, 
		 gretlopt oflags);

void get_cmd_ci (const char *line, CMD *command);

int loop_exec (LOOPSET *loop, 
	       double ***pZ, DATAINFO **ppdinfo, 
	       MODEL **models, PATHS *paths, 
	       int echo_off, PRN *prn);

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo);

int ifstate (int code);

#endif /* MONTE_CARLO_H */
