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

#ifndef CSVDATA_H
#define CSVDATA_H

#define import_na_string(s) (!strcmp(s, "NA") || \
                             !strcmp(s, "N.A.") || \
                             !strcmp(s, "n.a.") || \
                             !strcmp(s, "na") || \
                             !strcmp(s, "N/A") || \
			     !strcmp(s, "#N/A") || \
                             !strcmp(s, "NaN") || \
                             !strcmp(s, ".NaN") || \
                             !strcmp(s, ".") || \
                             !strcmp(s, "..") || \
                             !strcmp(s, "-999") || \
                             !strcmp(s, "-9999"))

int import_obs_label (const char *s);

int test_markers_for_dates (double ***pZ, DATAINFO *pdinfo, 
			    int *reversed, char *skipstr, 
			    PRN *prn);

void reverse_data (double **Z, DATAINFO *pdinfo, PRN *prn);

#endif /* CSVDATA_H */
