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

#ifndef FOREIGN_DB_H
#define FOREIGN_DB_H

int get_rats_series_info (const char *series_name,
                          SERIESINFO *sinfo);

dbwrapper *read_rats_db (const char *fname, FILE *fp);

int get_rats_db_data (const char *fname,
                      SERIESINFO *sinfo,
		      double **Z);

int get_pcgive_series_info (const char *series,
                            SERIESINFO *sinfo);

int get_pcgive_db_data (const char *dbbase,
                        SERIESINFO *sinfo,
			double **Z);

dbwrapper *read_pcgive_db (const char *fname, FILE *fp);

#endif /* FOREIGN_DB_H */
