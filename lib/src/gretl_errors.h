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

/* errors.h -- prototypes of functions in errors.c */

#include <stdio.h>

extern int gretl_errno;
extern char gretl_errmsg[ERRLEN];
extern char gretl_msg[ERRLEN];

enum gretl_error_codes {
    E_DATA = 2,
    E_SINGULAR,     /* 3 */
    E_DF,           /* 4 */
    E_YPY,          /* 5 */
    E_ZERO,         /* 6 */
    E_TSS,          /* 7 */
    E_ESS,          /* 8 */
    E_NOTALPH,      /* 9 */
    E_CONST,       /* 10 */
    E_UNBAL,       /* 11 */
    E_NEST,        /* 12 */
    E_NOTINTG,     /* 13 */
    E_NOTIMP,      /* 14 */
    E_IGNONZERO,   /* 15 */
    E_CASEU,       /* 16 */
    E_UNSPEC,      /* 17 */
    E_SYNTAX,      /* 18 */
    E_ER,          /* 19 */
    E_BADOP,       /* 20 */
    E_PDWRONG,     /* 21 */
    E_OFLAG,       /* 22 */
    E_FOPEN,       /* 23 */
    E_ALLOC,       /* 24 */
    E_NOEQ,        /* 25 */
    E_EQN,         /* 26 */
    E_UNKVAR,      /* 27 */
    E_NODATA,      /* 28 */
    E_ARGS,        /* 29 */
    E_OLSONLY,     /* 30 */
    E_INVARG,      /* 31 */
    E_ADF,         /* 32 */
    E_SPLIT,       /* 33 */
    E_PARSE,       /* 34 */
    E_NOVARS,      /* 35 */
    E_NOOMIT,      /* 36 */
    E_VARCHANGE,   /* 37 */
    E_NOADD,       /* 38 */
    E_ADDDUP,      /* 39 */
    E_LOGS,        /* 40 */
    E_SQUARES,     /* 41 */
    E_LAGS,        /* 42 */
    E_RHO,         /* 43 */
    E_SQRT,        /* 44 */
    E_HIGH,        /* 45 */
    E_WTZERO,      /* 46 */
    E_OBS,         /* 47 */
    E_NOVAR,       /* 48 */
    E_NOCONST,     /* 49 */
    E_MISS,        /* 50 */
    E_BADSTAT,     /* 51 */
    E_NOMERGE,     /* 52 */
    E_NOCONV,      /* 53 */
    E_CANCEL,      /* 54 */
    E_NAN,         /* 55 */
    E_MAX          /* 56 */
}; 

/* functions follow */
 
char *get_errmsg (const int errcode, char *errtext, PRN *prn);

void errmsg (const int errcode, PRN *prn);

int get_gretl_errno (void);

const char *get_gretl_errmsg (void);

const char *get_gretl_msg (void);

int print_gretl_msg (PRN *prn);

int print_gretl_errmsg (PRN *prn);

void gretl_errmsg_set (const char *str);
