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

/* functions follow */

double tprob (double x, int df);

double fdist (double x, int dfn, int dfd);

double chisq (double x, int df);

double normal (double x);

double rhocrit95 (int n);

double batch_pvalue (const char *str, 
		     double **Z, const DATAINFO *pdinfo, 
                     PRN *prn);

void interact_pvalue (void);

double f_crit_a (double a, int df1, int df2);

int print_critical (const char *line, PATHS *ppaths, PRN *prn);


