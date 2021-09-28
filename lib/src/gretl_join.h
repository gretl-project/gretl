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

#ifndef GRETL_JOIN_H
#define GRETL_JOIN_H

typedef enum {
    AGGR_NONE,
    AGGR_COUNT,
    AGGR_AVG,
    AGGR_SUM,
    AGGR_MIN,
    AGGR_MAX,
    AGGR_SEQ,
    AGGR_MIDAS
} AggrType;

int gretl_join_data (const char *fname,
                     const char **vnames,
                     int nvars,
                     DATASET *dset,
                     const int *ikeyvars,
                     const char *okey,
                     const char *filtstr,
                     const char *srcname,
                     AggrType aggr,
                     int seqval,
                     const char *auxname,
                     const char *tconvstr,
                     const char *tconvfmt,
                     int midas_pd,
                     gretlopt opt,
                     PRN *prn);

#endif /* GRETL_JOIN_H */
