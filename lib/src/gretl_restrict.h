/*
 *  Copyright (c) 2004 by Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef GRETL_RESTRICT_H
#define GRETL_RESTRICT_H

typedef struct restriction_set_ gretl_restriction_set;

gretl_restriction_set *
restriction_set_start (const char *line, gretlopt opt, int *err);

gretl_restriction_set *
cross_restriction_set_start (const char *line, gretl_equation_system *sys);

gretl_restriction_set *
var_restriction_set_start (const char *line, GRETL_VAR *var);

gretl_restriction_set *
eqn_restriction_set_start (const char *line, MODEL *pmod);

int 
restriction_set_parse_line (gretl_restriction_set *rset, const char *line,
			    const DATAINFO *pdinfo);

int
gretl_restriction_set_finalize (gretl_restriction_set *rset, 
				const double **Z, const DATAINFO *pdinfo,
				PRN *prn);

int gretl_sum_test (const int *list, MODEL *pmod, const DATAINFO *pdinfo,
		    PRN *prn);

#endif /* GRETL_RESTRICT_H */


