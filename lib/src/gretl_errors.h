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

#ifndef GRETL_ERRORS_H
#define GRETL_ERRORS_H

extern char gretl_errmsg[ERRLEN];

enum gretl_error_codes {
    E_DATA = 2,
    E_SINGULAR,     /* 3 */
    E_DF,           /* 4 */
    E_ZERO,         /* 5 */
    E_TSS,          /* 6 */
    E_ESS,          /* 7 */
    E_NOTIMP,       /* 8 */
    E_UNSPEC,       /* 9 */
    E_SYNTAX,      /* 10 */
    E_PDWRONG,     /* 11 */
    E_FOPEN,       /* 12 */
    E_ALLOC,       /* 13 */
    E_EQN,         /* 14 */
    E_UNKVAR,      /* 15 */
    E_ARGS,        /* 16 */
    E_OLSONLY,     /* 17 */
    E_INVARG,      /* 18 */
    E_PARSE,       /* 19 */
    E_NOVARS,      /* 20 */
    E_NOOMIT,      /* 21 */
    E_NOADD,       /* 22 */
    E_ADDDUP,      /* 23 */
    E_LOGS,        /* 24 */
    E_SQUARES,     /* 25 */
    E_LAGS,        /* 26 */
    E_SQRT,        /* 27 */
    E_HIGH,        /* 28 */
    E_OBS,         /* 29 */
    E_NOCONST,     /* 30 */
    E_BADSTAT,     /* 31 */
    E_NOMERGE,     /* 32 */
    E_NOCONV,      /* 33 */
    E_CANCEL,      /* 34 */
    E_MISSDATA,    /* 35 */
    E_NAN,         /* 36 */
    E_NONCONF,     /* 37 */
    E_TYPES,       /* 38 */
    E_DATATYPE,    /* 39 */
    E_BADOPT,      /* 40 */
    E_NOIDENT,     /* 41 */
    E_EXTERNAL,    /* 42 */
    E_TOOLONG,     /* 43 : command line too long */
    E_DB_DUP,      /* 44 : duplicate vars found when saving to database */
    E_OK,          /* 45 : not really an error */
    E_MAX          /* 46 */
}; 

void errmsg (int err, PRN *prn);

const char *errmsg_get_with_default (int err);

const char *gretl_errmsg_get (void);

void gretl_errmsg_set (const char *str);

void gretl_errmsg_sprintf (const char *fmt, ...);

void gretl_errmsg_set_from_errno (void);

void gretl_error_clear (void);

#endif /* GRETL_ERRORS_H */
