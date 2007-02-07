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

/* modelprint.h for gretl */

#ifndef MODELPRINT_H
#define MODELPRINT_H

int printmodel (MODEL *pmod, const DATAINFO *pdinfo, gretlopt opt,
		PRN *prn);

const char *estimator_string (int ci, PRN *prn);

void print_model_vcv_info (const MODEL *pmod, PRN *prn);

int ols_print_anova (const MODEL *pmod, PRN *prn);

#endif /* MODELPRINT_H */
