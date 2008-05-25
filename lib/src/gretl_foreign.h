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

#ifndef GRETL_FOREIGN_H

int foreign_append_line(const char *line, gretlopt opt, PRN *prn);

int foreign_execute (const double **Z, const DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn);

int write_gretl_R_files (const char *buf,
			 const double **Z, const DATAINFO *pdinfo,
			 gretlopt opt);

void delete_gretl_R_files (void);

#endif /* GRETL_FOREIGN_H */
