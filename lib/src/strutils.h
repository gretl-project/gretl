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

/* functions follow */

int string_is_blank (const char *s);

int has_suffix (const char *str, const char *sfx);

int numeric_string (const char *str);

int count_fields (const char *s);

double dot_atof (const char *s);
 
int dotpos (const char *str);

int slashpos (const char *str);

char *delchar (int c, char *str);

int haschar (char c, const char *s);

int lastchar (char c, const char *s);

int ends_with_backslash (const char *s);

int gretl_varchar_spn (const char *s);

char *gretl_trunc (char *str, size_t n);

char *gretl_delete (char *str, int idx, int count);

char *gretl_strdup (const char *src);

char *gretl_strndup (const char *src, size_t n);

char *gretl_word_strdup (const char *src, const char **ptr);

char *charsub (char *str, char find, char repl);

char *shift_string_left (char *str, size_t move);

char *lower (char *str);

void clear (char *str, int len);

char *chopstr (char *str);

char *switch_ext (char *targ, const char *src, char *ext);

int get_base (char *targ, const char *src, char c);

int equation_get_lhs_and_rhs (const char *s, char **plh, char **prh);

int top_n_tail (char *str);

char *tailstrip (char *str);

char *compress_spaces (char *s);

char *space_to_score (char *s);

char *safecpy (char *targ, const char *src, int n);

char **create_strings_array (int nstrs);

void free_strings_array (char **strs, int nstrs);

char *get_obs_string (char *obs, int t, const DATAINFO *pdinfo);

char *get_full_obs_string (char *obs, int t, const DATAINFO *pdinfo);

double obs_str_to_double (const char *obs);

char *colonize_obs (char *obs);

void modify_date_for_csv (char *s, int pd);

void csv_obs_to_prn (int t, const DATAINFO *pdinfo, PRN *prn);

const char *print_time (const time_t *timep);

char *gretl_xml_encode (char *buf);

void unescape_url (char *url);

char *make_varname_unique (char *vname, int v, DATAINFO *pdinfo);

int fix_varname_duplicates (DATAINFO *pdinfo);

char *append_dir (char *fname, const char *dir);

int build_path (const char *dir, const char *fname, char *path, 
		const char *ext);

#if defined(USE_GTK2) || defined (HAVE_FNMATCH_H)
int *varname_match_list (const DATAINFO *pdinfo, const char *pattern);
#endif

#endif /* STRUTILS_H */
