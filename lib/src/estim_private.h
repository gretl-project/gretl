/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) Allin Cottrell
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

#ifndef ESTIM_PRIVATE_H
#define ESTIM_PRIVATE_H

double dwstat (int order, MODEL *pmod, const double **Z);

double rhohat (int order, int t1, int t2, const double *uhat);

#endif /* ESTIM_PRIVATE_H */
