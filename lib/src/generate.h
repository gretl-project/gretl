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

#ifndef GENERATE_H
#define GENERATE_H

#include <time.h>

#define _HIGHVALU 87     /* maximum value for exponent in genr formula */
#define _VSMALL 1.0e-14  /* adjustment for (int) conversion */

typedef enum {
    UHATNUM = 999,
    TNUM,
    INDEXNUM
} genr_numbers;

typedef struct {
    double *xvec;
    int varnum;
    char varname[9];
    char label[MAXLABEL];
    int errcode;
    char msg[MAXLABEL];
    int special;
    int scalar; 
} GENERATE;

typedef struct {
    int lag;
    int varnum;
    char varname[9];
} LAGVAR;

/* functions follow */
 
GENERATE generate (double ***pZ, DATAINFO *pdinfo, 
		   const char *line, int model_count, 
		   MODEL *pmod, int oflag);

int dummy (double ***pZ, DATAINFO *pdinfo);

int paneldum (double ***pZ, DATAINFO *pdinfo, 
	      int opt);

int plotvar (double ***pZ, DATAINFO *pdinfo, 
	     const char *period);

void varlist (const DATAINFO *pdinfo, PRN *prn);

int varindex (const DATAINFO *pdinfo, const char *varname);

int logs (const LIST list, 
	  double ***pZ, DATAINFO *pdinfo);

int lags (const LIST list, 
	  double ***pZ, DATAINFO *pdinfo);

int xpxgenr (const LIST list, 
	     double ***pZ, DATAINFO *pdinfo, 
	     const int opt, const int nodup);

int rhodiff (char *param, const LIST list, 
	     double ***pZ, DATAINFO *pdinfo);

int simulate (char *cmd, 
	      double ***pZ, DATAINFO *pdinfo); 

int genr_fit_resid (MODEL *pmod, 
		    double ***pZ, DATAINFO *pdinfo,
		    int code, int undo);

#endif /* GENERATE_H */

