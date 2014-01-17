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

typedef enum {
    E_DATA = 2,
    E_SINGULAR,     /* 3 */
    E_DF,           /* 4 */
    E_ZERO,         /* 5 */
    E_TSS,          /* 6 */
    E_ESS,          /* 7 */
    E_NOTIMP,       /* 8 */
    E_UNSPEC,       /* 9 */
    E_PDWRONG,     /* 10 */
    E_FOPEN,       /* 11 */
    E_ALLOC,       /* 12 */
    E_EQN,         /* 13 */
    E_UNKVAR,      /* 14 */
    E_ARGS,        /* 15 */
    E_OLSONLY,     /* 16 */
    E_INVARG,      /* 17 */
    E_PARSE,       /* 18 */
    E_NOVARS,      /* 19 */
    E_NOOMIT,      /* 20 */
    E_NOADD,       /* 21 */
    E_ADDDUP,      /* 22 */
    E_LOGS,        /* 23 */
    E_SQUARES,     /* 24 */
    E_LAGS,        /* 25 */
    E_SQRT,        /* 26 */
    E_HIGH,        /* 27 */
    E_OBS,         /* 28 */
    E_NOCONST,     /* 29 */
    E_BADSTAT,     /* 30 */
    E_NOMERGE,     /* 31 */
    E_NOCONV,      /* 32 */
    E_CANCEL,      /* 33 */
    E_MISSDATA,    /* 34 */
    E_NAN,         /* 35 */
    E_NONCONF,     /* 36 */
    E_TYPES,       /* 37 */
    E_DATATYPE,    /* 38 */
    E_BADOPT,      /* 39 */
    E_NOIDENT,     /* 40 */
    E_EXTERNAL,    /* 41 */
    E_TOOLONG,     /* 42 */
    E_NODATA,      /* 43 */
    E_NOTPD,       /* 44 */
    E_JACOBIAN,    /* 45 */
    E_TOOFEW,      /* 46 */
    E_FUNCERR,     /* 47 : error set by function writer */
    E_STOP,        /* 48 : user aborted execution */
    E_DB_DUP,      /* 49 : duplicate vars found when saving to database */
    E_OK,          /* 50 : not really an error */
    E_MAX          /* 51 */
} GretlError; 

enum gretl_warning_codes {
    W_GRADIENT = 1,
    W_GENMISS,     /* 2 */
    W_GENNAN,      /* 3 */
    W_MAX          /* 4 */
};

void errmsg (int err, PRN *prn);

void warnmsg (PRN *prn);

const char *errmsg_get_with_default (int err);

const char *gretl_errmsg_get (void);

const char *gretl_warnmsg_get (void);

char *maybe_save_gretl_errmsg (int err);

void gretl_errmsg_set (const char *str);

void gretl_errmsg_ensure (const char *str);

void gretl_warnmsg_set (const char *str);

void gretl_errmsg_sprintf (const char *fmt, ...);

void gretl_warnmsg_sprintf (const char *fmt, ...);

void gretl_errmsg_set_from_errno (const char *s);

void gretl_error_clear (void);

void set_gretl_alarm (int val);

void set_gretl_errno (int err);

void set_gretl_warning (int w);

int get_gretl_errno (void);

int check_gretl_errno (void);

int check_gretl_warning (void);

int gretl_error_is_fatal (void);

int gretl_errmsg_is_set (void);

#endif /* GRETL_ERRORS_H */
