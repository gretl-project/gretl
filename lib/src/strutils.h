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

int has_native_data_suffix (const char *fname);

int numeric_string (const char *str);

int integer_string (const char *str);

int count_fields (const char *s, const char *sep);

int count_lines (const char *s);

double dot_atof (const char *s);

void set_atof_point (char c);
 
int gretl_dotpos (const char *str);

int gretl_slashpos (const char *str);

char *strrslash (const char *s);

char *gretl_delchar (int c, char *str);

int gretl_charpos (char c, const char *s);

int ends_with_backslash (const char *s);

int is_greek_letter (const char *s);

int gretl_namechar_spn (const char *s);

int double_quote_position (const char *s);

char *gretl_trunc (char *str, size_t n);

char *gretl_delete (char *str, int idx, int count);

char *gretl_unquote (char *str, int *err);

char *gretl_strdup (const char *src);

char *gretl_strndup (const char *src, size_t n);

char *gretl_strdup_printf (const char *format, ...);

char *gretl_word_strdup (const char *src, const char **ptr,
			 gretlopt opt, int *err);

char *gretl_quoted_string_strdup (const char *s, const char **ptr);

char **gretl_string_split (const char *s, int *n, const char *sep);

char **gretl_string_split_quoted (const char *s, int *n, 
				  const char *sep, int *err);

char **gretl_string_split_lines (const char *s, int *n);

char *gretl_str_expand (char **orig, const char *add, const char *sep);

char *gretl_charsub (char *str, char find, char repl);

char *gretl_substring (const char *str, int first, int last, int *err);

char *comma_separate_numbers (char *s);

char *shift_string_left (char *str, size_t move);

char *gretl_lower (char *str);

char *gretl_strstrip (char *str);

char *gretl_strstrip_copy (const char *str, int *err);

char *switch_ext (char *targ, const char *src, const char *ext);

char *switch_ext_in_place (char *fname, const char *ext);

char *switch_ext_new (const char *src, const char *ext);

int equation_get_lhs_and_rhs (const char *s, char **plh, char **prh);

int top_n_tail (char *str, size_t maxlen, int *err);

char *tailstrip (char *str);

char *compress_spaces (char *s);

char *space_to_score (char *s);

char **strings_array_new (int nstrs);

char **strings_array_realloc_with_length (char ***pS, 
					  int oldn, 
					  int newn,
					  int len);

int strings_array_add (char ***pS, int *n, const char *p);

int strings_array_donate (char ***pS, int *n, char *p);

int strings_array_add_uniq (char ***pS, int *n, const char *p,
			    int *pos);

int strings_array_prepend_uniq (char ***pS, int *n, const char *p);

char **strings_array_new_with_length (int nstrs, int len);

char **strings_array_dup (char **strs, int n);

char **strings_array_dup_selected (char **strs, int n,
				   const int *list);

int strings_array_sort (char ***pS, int *n, gretlopt opt);

int strings_array_cmp (char **strs1, char **strs2, int n);

int strings_array_position (char **strs, int n, const char *s);

int strings_array_diff (char **strs1, int n1,
			char **strs2, int n2,
			char ***extra, int *n_extra);

char **strings_array_reverse (char **strs, int nstrs);

void strings_array_free (char **strs, int nstrs);

char *get_obs_string (char *obs, int t, const DATASET *dset);

double obs_str_to_double (const char *obs);

char *colonize_obs (char *obs);

void modify_date_for_csv (char *s, int pd);

char *print_time (char *s);

int gretl_xml_validate (const char *s);

char *gretl_xml_encode (const char *str);

int gretl_xml_encode_to_buf (char *targ, const char *src, int n);

void unescape_url (char *url);

char *make_varname_unique (char *vname, int v, DATASET *dset);

int fix_varname_duplicates (DATASET *dset);

char *append_dir (char *fname, const char *dir);

const char *path_last_element (const char *path);

char *trim_slash (char *s);

int gretl_string_ends_with (const char *s, const char *test);

void get_column_widths (const char **strs, int *widths, int n);

char *gretl_utf8_strncat (char *dest, const char *src, size_t n);

char *gretl_utf8_strncat_trim (char *dest, const char *src, size_t n);

char *gretl_utf8_truncate (char *s, size_t nmax);

char *gretl_utf8_truncate_b (char *s, size_t bmax);

char *gretl_utf8_replace_char (const char *targ, const char *src,
			       int pos);

char *gretl_utf8_select (const char *s, const int *list);

int gretl_scan_varname (const char *src, char *targ);

int gretl_normalize_varname (char *targ, const char *src,
			     int underscore, int seq);

char *gretl_regexp_replace (const char *orig,
			    const char *match,
			    const char *repl,
			    int *err);

char *gretl_literal_replace (const char *orig,
			     const char *match,
			     const char *repl,
			     int *err);

#endif /* STRUTILS_H */
