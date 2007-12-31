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

char *gretl_strdup_printf (const char *template, ...);

char *gretl_word_strdup (const char *src, const char **ptr);

char **gretl_string_split (const char *s, int *n);

char *gretl_str_expand (char **orig, const char *add, const char *sep);

double gretl_double_from_string (const char *s, const char **ptr);

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

char **strings_array_new (int nstrs);

int strings_array_add (char ***pS, int *n, const char *p);

char **strings_array_new_with_length (int nstrs, int len);

char **strings_array_dup (char **strs, int n);

void free_strings_array (char **strs, int nstrs);

char *get_obs_string (char *obs, int t, const DATAINFO *pdinfo);

char *get_full_obs_string (char *obs, int t, const DATAINFO *pdinfo);

double obs_str_to_double (const char *obs);

char *colonize_obs (char *obs);

void modify_date_for_csv (char *s, int pd);

void csv_obs_to_prn (int t, const DATAINFO *pdinfo, PRN *prn);

const char *print_time (const time_t *timep);

int gretl_xml_validate (const char *s);

char *gretl_xml_encode (const char *str);

int gretl_xml_encode_to_buf (char *targ, const char *src, int n);

void unescape_url (char *url);

char *make_varname_unique (char *vname, int v, DATAINFO *pdinfo);

int fix_varname_duplicates (DATAINFO *pdinfo);

char *append_dir (char *fname, const char *dir);

char *build_path (char *targ, const char *dirname, const char *fname, 
		  const char *ext);

const char *path_last_element (const char *path);

char *trim_slash (char *s);

int *varname_match_list (const DATAINFO *pdinfo, const char *pattern);

#endif /* STRUTILS_H */
