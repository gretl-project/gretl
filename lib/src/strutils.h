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

/* strutils.h for gretl */

#ifndef STRUTILS_H
#define STRUTILS_H

#include "generate.h"

#ifdef OS_WIN32
#define SLASH '\\'
#define SLASHSTR "\\"
#else
#define SLASH '/'
#define SLASHSTR "/"
#endif

/* functions follow */
 
int dotpos (const char *str);

int slashpos (const char *str);

void delchar (const int c, char *str);

int haschar (const char c, const char *str);

void lower (char *str);

void clear (char *str, const int len);

void chopstr (char *str);

char *switch_ext (char *targ, const char *src, char *ext);

int get_base (char *targ, const char *src, char c);

int top_n_tail (char *str);

void compress_spaces (char *str);

int pprintf (print_t *prn, const char *template, ...);

char *safecpy (char *targ, const char *src, int n);

#endif /* STRUTILS_H */
