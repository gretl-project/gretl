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

#include "libgretl.h"
#include <time.h>

#ifdef WIN32
#define SLASH '\\'
#define SLASHSTR "\\"
#else
#define SLASH '/'
#define SLASHSTR "/"
#endif

#define CTRLZ 26

extern char gretl_tmp_str[MAXLEN];

/* functions follow */

int string_is_blank (const char *s);

double dot_atof (const char *s);
 
size_t dotpos (const char *str);

int slashpos (const char *str);

void delchar (int c, char *str);

int haschar (char c, const char *s);

int lastchar (char c, const char *s);

int ends_with_backslash (const char *s);

char *gretl_strdup (const char *src);

char *charsub (char *str, char find, char repl);

void lower (char *str);

void clear (char *str, int len);

void chopstr (char *str);

char *switch_ext (char *targ, const char *src, char *ext);

int get_base (char *targ, const char *src, char c);

int top_n_tail (char *str);

void compress_spaces (char *str);

char *space_to_score (char *s);

char *safecpy (char *targ, const char *src, int n);

int doing_nls (void);

int reset_local_decpoint (void);

int get_local_decpoint (void);

char *get_obs_string (char *obs, int t, const DATAINFO *pdinfo);

char *get_full_obs_string (char *obs, int t, const DATAINFO *pdinfo);

double obs_str_to_double (const char *obs);

char *colonize_obs (char *obs);

#ifdef ENABLE_NLS
char *iso_gettext (const char *msgid);

char *maybe_iso_gettext (const char *msgid);

void set_gretl_charset (const char *s);

const char *get_gretl_charset (void);

const char *get_gnuplot_charset (void);

int get_utf_width (const char *str, int width);

# define UTF_WIDTH(s, w) get_utf_width(s, w) 
int get_utf_width (const char *str, int width);
#else
# define UTF_WIDTH(s, w)    w
#endif  /* ENABLE_NLS */

const char *print_time (const time_t *timep);

char *gretl_xml_encode (char *buf);

void unescape_url (char *url);

char *iso_to_ascii (char *s);

char *make_varname_unique (char *vname, int v, DATAINFO *pdinfo);

char *append_dir (char *fname, const char *dir);

int build_path (const char *dir, const char *fname, char *path, 
		const char *ext);

#ifndef USE_GTK2
int
utf8_to_iso_latin_1 (unsigned char* out, int outlen, 
		     unsigned char* in, int inlen);
#endif

#if defined(USE_GTK2) || defined (HAVE_FNMATCH_H)
int *varname_match_list (const DATAINFO *pdinfo, const char *pattern);
#endif

#endif /* STRUTILS_H */
