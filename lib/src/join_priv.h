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

enum {
    TCONV_FMT = 0,
    TKEY_FMT = 1
};

#define no_formats(map) (map.fmt == NULL)
#define no_tkey_format(map) (map.tname == NULL)
#define has_tconv_format(map) (map.fmt[TCONV_FMT] != NULL)
#define is_tkey_variable(name, map) (strcmp(name, map.tname) == 0)

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
typedef struct csvprobe_ csvprobe;

struct time_mapper {
    int ncols;         /* number of "timeconv" columns */
    char **colnames;   /* array of outer-dataset column names */
    char *tname;       /* the name of the "tkey", if among colnames, or NULL */
    char **fmt;        /* array of up to two time-format strings, or NULL */
    char m_means_q[2]; /* array of "monthly means quarterly" flags */
};

extern struct time_mapper tconv_map;

DATASET *csvdata_get_dataset (csvdata *c);

void csvdata_free (csvdata *c);

int real_import_csv (const char *fname,
		     DATASET *dset,
		     const char *cols,
		     const char *rows,
		     joinspec *join,
		     csvprobe *probe,
		     gretl_matrix **pm,
		     gretlopt opt,
		     PRN *prn);

#endif /* JOIN_PRIV_H */
