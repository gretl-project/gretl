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

typedef enum {
    E_DATA = 2,
    E_SINGULAR,
    E_DF,
    E_YPY,
    E_ZERO,
    E_TSS,
    E_ESS,
    E_NOTALPH,
    E_CONST,
    E_UNBAL,
    E_NEST,
    E_NOTINTG,
    E_NOTIMP,
    E_IGNONZERO,
    E_CASEU,
    E_UNSPEC,
    E_SYNTAX,
    E_ER,
    E_BADOP,
    E_PDWRONG,
    E_OFLAG,
    E_FOPEN,
    E_ALLOC,
    E_NOEQ,
    E_EQN,
    E_UNKVAR,
    E_NODATA,
    E_ARGS,
    E_OLSONLY,
    E_INVARG,
    E_ADF,
    E_SPLIT,
    E_PARSE,
    E_NOVARS,
    E_NOOMIT,
    E_VARCHANGE,
    E_NOADD,
    E_ADDDUP,
    E_LOGS,
    E_SQUARES,
    E_LAGS,
    E_RHO,
    E_SQRT,
    E_HIGH,
    E_WTZERO,
    E_OBS,
    E_NOVAR,
    E_NOCONST,
    E_MISS,
    E_BADSTAT,
    E_NOMERGE,
    E_NAN,
    E_MAX
} error_codes; 

/* functions follow */
 
char *get_errmsg (const int errcode, char *errtext, PRN *prn);

void errmsg (const int errcode, PRN *prn);

int get_gretl_errno (void);

char *get_gretl_errmsg (void);
