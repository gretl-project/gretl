/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 2005 Allin Cottrell
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

/* gretl: private communication between monte_carlo.c (loop code),
   generate.c and interact.c (for boolean flow control)
*/

#ifndef LOOP_PRIVATE_H
#define LOOP_PRIVATE_H

enum {
    SET_FALSE,
    SET_TRUE,
    SET_ELSE,
    SET_ELIF,
    SET_ENDIF,
    IS_FALSE,
    DOINDENT,
    UNINDENT,
    GETINDENT,
    RELAX
};

int is_active_index_loop_char (int c);

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo);

int ifstate (int code);

#endif /* LOOP_PRIVATE_H */
