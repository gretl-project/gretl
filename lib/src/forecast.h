/*
 * Copyright (C) 1999-2005 Allin Cottrell
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

#ifndef FORECAST_H
#define FORECAST_H

/* estimators where a simple X*b does _not_ give the
   predicted value of the dependent variable */

#define FCAST_SPECIAL(c) (c == LOGIT || \
                          c == NLS || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == TOBIT)

FITRESID *fit_resid_new (int n, int errs);

void free_fit_resid (FITRESID *fr);


FITRESID *get_fit_resid (const MODEL *pmod, const double **Z, 
			 const DATAINFO *pdinfo);

FITRESID *get_fcast_with_errs (const char *str, MODEL *pmod, 
			       const double **Z, const DATAINFO *pdinfo); 

FITRESID *get_fcast_without_errs (const char *str, MODEL *pmod, 
				  double ***pZ, DATAINFO *pdinfo);

int fcast_with_errs (const char *str, MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, 
		     gretlopt opt, PRN *prn);

int fcast (const char *line, const MODEL *pmod, double ***pZ,
	   DATAINFO *pdinfo);

#endif /* FORECAST_H */


