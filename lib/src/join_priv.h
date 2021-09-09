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

#ifndef JOIN_PRIV_H
#define JOIN_PRIV_H

/* stuff that needs to be shared between csvdata.c and gretl_join.c */

typedef struct csvdata_ csvdata;

struct joinspec_ {
    int ncols;
    const char **colnames;
    const char *mdsbase;
    int *colnums;
    int *timecols;
    csvdata *c;
    DATASET *dset;
    int wildcard;
    int auto_midas;
    int midas_pd;
    char **wildnames;
    char **mdsnames;
    char **tmpnames;
    int n_tmp;
};

typedef struct joinspec_ joinspec;

DATASET *csvdata_get_dataset (csvdata *c);

void csvdata_free (csvdata *c);

int real_import_csv (const char *fname,
		     DATASET *dset,
		     const char *cols,
		     const char *rows,
		     joinspec *join,
		     void *probe,
		     gretl_matrix **pm,
		     gretlopt opt,
		     PRN *prn);

int timecol_get_format (const DATASET *dset, int v,
			char **pfmt, int *q);

#endif /* JOIN_PRIV_H */
