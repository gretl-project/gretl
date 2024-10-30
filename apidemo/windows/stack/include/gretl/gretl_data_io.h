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

#ifndef GRETL_DATA_IO_H
#define GRETL_DATA_IO_H

int gretl_read_native_data (const char *fname, DATASET *dset);

int gretl_write_native_data (const char *fname, const int *list,
			     const DATASET *dset);

int gretl_read_foreign_data (const char *fname, GretlFileType file_type,
			     DATASET *dset, PRN *prn);

GretlFileType gretl_detect_filetype (const char *fname);


#endif /* GRETL_DATA_IO_H */

