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

#ifndef GRETL_ERRORS_H
#define GRETL_ERRORS_H

extern char gretl_errmsg[ERRLEN];

enum gretl_error_codes {
    E_DATA = 2,
    E_SINGULAR,     /* 3 */
    E_DF,           /* 4 */
    E_YPY,          /* 5 */
    E_ZERO,         /* 6 */
    E_TSS,          /* 7 */
    E_ESS,          /* 8 */
    E_UNBAL,        /* 9 */
    E_NOTIMP,      /* 10 */
    E_UNSPEC,      /* 11 */
    E_SYNTAX,      /* 12 */
    E_PDWRONG,     /* 13 */
    E_FOPEN,       /* 14 */
    E_ALLOC,       /* 15 */
    E_EQN,         /* 16 */
    E_UNKVAR,      /* 17 */
    E_NODATA,      /* 18 */
    E_ARGS,        /* 19 */
    E_OLSONLY,     /* 20 */
    E_INVARG,      /* 21 */
    E_SPLIT,       /* 22 */
    E_PARSE,       /* 23 */
    E_NOVARS,      /* 24 */
    E_NOOMIT,      /* 25 */
    E_VARCHANGE,   /* 26 */
    E_NOADD,       /* 27 */
    E_ADDDUP,      /* 28 */
    E_LOGS,        /* 29 */
    E_SQUARES,     /* 30 */
    E_LAGS,        /* 31 */
    E_SQRT,        /* 32 */
    E_HIGH,        /* 33 */
    E_WTZERO,      /* 34 */
    E_WTNEG,       /* 35 */
    E_OBS,         /* 36 */
    E_NOCONST,     /* 37 */
    E_MISS,        /* 38 */
    E_BADSTAT,     /* 39 */
    E_NOMERGE,     /* 40 */
    E_NOCONV,      /* 41 */
    E_CANCEL,      /* 42 */
    E_MISSDATA,    /* 43 */
    E_NAN,         /* 44 */
    E_NONCONF,     /* 45 */
    E_TYPES,       /* 46 */
    E_DATATYPE,    /* 47 */
    E_BADOPT,      /* 48 */
    E_DB_DUP,      /* 49 : duplicate vars found when saving to database */
    E_OK,          /* 50 : not really an error */
    E_MAX          /* 51 */
}; 

void errmsg (int err, PRN *prn);

const char *errmsg_get_with_default (int err);

const char *gretl_errmsg_get (void);

void gretl_errmsg_set (const char *str);

void gretl_error_clear (void);

#endif /* GRETL_ERRORS_H */
