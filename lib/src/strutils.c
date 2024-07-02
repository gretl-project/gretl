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

/* strutils.c for gretl */

#include "libgretl.h"

#include <errno.h>
#include <time.h>
#include <glib.h>

/**
 * SECTION:strutils
 * @short_description: miscellaneous string-handling utilities
 * @title: Strings
 * @include: libgretl.h
 *
 * Various functions for creating, testing and manipulating
 * strings and arrays of strings.
 */

/**
 * string_is_blank:
 * @s: the string to examine.
 *
 * Returns: 1 if the string is NULL, of length zero, or contains
 * nothing but space characters, otherwise returns 0.
 **/

int string_is_blank (const char *s)
{
    int ret = 1;

    if (s != NULL) {
	while (*s) {
	    if (!isspace((unsigned char) *s) &&
		*s != '\r' && *s != CTRLZ) {
		ret = 0;
		break;
	    }
	    s++;
	}
    }

    return ret;
}

static int atof_point;

void set_atof_point (char c)
{
    atof_point = c;
}

/**
 * dot_atof:
 * @s: the string to convert.
 *
 * Returns: the double-precision numeric interpretation of @s,
 * where the decimal point character is forced to be '.',
 * regardless of the current locale.
 **/

double dot_atof (const char *s)
{
#ifndef ENABLE_NLS
    return atof(s);
#else
    double x;

    if (atof_point == 0) {
	struct lconv *lc = localeconv();

	atof_point = *lc->decimal_point;
    }

    if (atof_point == '.') {
	x = atof(s);
    } else {
	gretl_push_c_numeric_locale();
	x = atof(s);
	gretl_pop_c_numeric_locale();
    }

    return x;
#endif
}

/**
 * gretl_dotpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last "." within @str,
 * or strlen(@str) in case a dot is not found, or the string
 * ends with a (backward or forward) slash.
 */

int gretl_dotpos (const char *str)
{
    int i, p = 0;

    if (str != NULL && *str != '\0') {
	p = strlen(str);
	for (i=p-1; i>0; i--) {
	    if (str[i] == '/' || str[i] == '\\') {
		break;
	    } else if (str[i] == '.') {
		p = i;
		break;
	    }
	}
    }

    return p;
}

#ifdef WIN32
# define is_slash(c) (c == '\\' || c == '/')
#else
# define is_slash(c) (c == SLASH)
#endif

/**
 * gretl_slashpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last #SLASH within @str,
 * or 0 in case a #SLASH is not found.
 */

int gretl_slashpos (const char *str)
{
    int i, p = 0;

    if (str != NULL && *str != '\0') {
	p = strlen(str);
	for (i=p-1; i>0; i--) {
	    if (is_slash(str[i])) {
		p = i;
		break;
	    }
	}
    }

    return p;
}

/**
 * strrslash:
 * @s: the string to examine.
 *
 * Returns: a pointer to the last occurrence of 'SLASH'
 * in the string @d, making allowance for the fact that
 * on MS Windows this may be either a backslash or a
 * forward slash, or NULL is no 'SLASH' is found.
 */

char *strrslash (const char *s)
{
#ifdef WIN32
    char *p1 = strrchr(s, '\\');
    char *p2 = strrchr(s, '/');

    if (p1 != NULL && p2 == NULL) {
	return p1;
    } else if (p2 != NULL && p1 == NULL) {
	return p2;
    } else if (p1 != NULL && p2 != NULL) {
	return p2 - p1 > 0 ? p2 : p1;
    } else {
	return NULL;
    }
#else
    return strrchr(s, '/');
#endif
}

/**
 * gretl_delchar:
 * @c: the character to delete.
 * @str: the string from which to delete @c.
 *
 * Deletes all instances of @c within @str.
 *
 * Returns: the possibly modified string.
 */

char *gretl_delchar (int c, char *str)
{
    int i, j;

    for (i=j=0; str[i] != '\0'; i++) {
	if (str[i] != c) {
	    str[j++] = str[i];
	}
    }

    str[j] = '\0';

    return str;
}

/**
 * gretl_delete:
 * @str: the string to process.
 * @idx: the starting point for deleting characters.
 * @count: the number of characters to delete.
 *
 * Deletes @count characters from @str, starting at position @idx.
 *
 * Returns: the modified string.
 */

char *gretl_delete (char *str, int idx, int count)
{
    size_t i, n = strlen(str);

    for (i=idx; i<=n-count; ++i) {
	str[i] = str[count+i];
    }

    return str;
}

/**
 * gretl_unquote:
 * @str: the string to process.
 * @err: location to receive error code.
 *
 * If @str begins with the ASCII double-quote character, checks
 * that the last character is also a double-quote, and in that
 * case trims the quotes from both ends. If the first character
 * is a double quote but the last is not, flags an error. If
 * the string is not quoted at all, returns the original
 * string.
 *
 * Returns: the input string, possibly modified in place.
 */

char *gretl_unquote (char *str, int *err)
{
    *err = 0;

    if (*str == '"') {
	int n = strlen(str);

	if (n > 1) {
	    if (str[n-1] == '"') {
		str[n-1] = '\0';
	    } else {
		*err = E_PARSE;
	    }
	} else {
	    *err = E_PARSE;
	}

	if (!*err) {
	    shift_string_left(str, 1);
	}
    }

    return str;
}

/**
 * gretl_charpos:
 * @c: the character to look for.
 * @s: the string to examine.
 *
 * Returns: the first position of @c in @s, or -1 if @c is not
 * found.
 */

int gretl_charpos (char c, const char *s)
{
    int i = 0;

    while (*s) {
	if (*s++ == c) {
	    return i;
	}
	i++;
    }

    return -1;
}

/**
 * gretl_charsub:
 * @str: the string to operate on.
 * @find: the character to replace.
 * @repl: the replacement character.
 *
 * Replaces all occurrences of @find with @repl in @str.
 *
 * Returns: the (possibly modified) string.
 */

char *gretl_charsub (char *str, char find, char repl)
{
    char *p = str;

    while (*p) {
	if (*p == find) {
	    *p = repl;
	}
	p++;
    }

    return str;
}

/**
 * comma_separate_numbers:
 * @s: the string to operate on.
 *
 * Given a string which contains two or more numbers
 * separated by spaces and/or commas, revise the
 * string to ensure that all the numbers are comma-separated.
 *
 * Returns: the (possibly modified) string.
 */

char *comma_separate_numbers (char *s)
{
    const char *numstart = "+-.0123456789";
    char *p = s;
    int i, n, done;

    while (*s) {
	n = strspn(s, " ,");
	if (n > 0 && s[n] != '\0' && strchr(numstart, s[n])) {
	    done = 0;
	    for (i=0; i<n && !done; i++) {
		if (s[i] == ',') {
		    done = 1;
		}
	    }
	    if (!done) {
		*s = ',';
	    }
	}
	s += (n > 0)? n : 1;
    }

    return p;
}

/**
 * has_suffix:
 * @str: the string to check.
 * @sfx: the suffix to check for, including the leading '.'
 *
 * Returns: 1 if @str ends with @sfx (on a case-insensitive
 * comparison), 0 otherwise.
 */

 int has_suffix (const char *str, const char *sfx)
 {
     int offset;

     if (str == NULL || sfx == NULL) {
	 return 0;
     }

     offset = strlen(str) - strlen(sfx);
     if (offset < 0) {
	 return 0;
     }

     return g_ascii_strcasecmp(str + offset, sfx) == 0;
}

/**
 * has_native_data_suffix:
 * @fname: the filename to check.
 *
 * Returns: 1 if @fname ends with a suffix indicating it is a
 * native gretl data file, 0 otherwise.
 */

int has_native_data_suffix (const char *fname)
{
    return has_suffix(fname, ".gdt") || has_suffix(fname, ".gdtb");
}

/**
 * numeric_string:
 * @str: the string to examine.
 *
 * Returns: 1 if the given @str is numeric, otherwise 0.
 */

int numeric_string (const char *str)
{
    char *test;
    int ret = 1;

    if (str == NULL || *str == '\0') {
	return 0;
    }

    if (strlen(str) == 3) {
	char screen[4];

	strcpy(screen, str);
	gretl_lower(screen);
	if (!strcmp(screen, "inf") || !strcmp(screen, "nan")) {
	    /* could be variable names: they are not valid numbers */
	    return 0;
	}
    }

    gretl_push_c_numeric_locale();
    errno = 0;
    strtod(str, &test);
    gretl_pop_c_numeric_locale();

    if (*test != '\0' || errno == ERANGE) {
	ret = 0;
    }

    return ret;
}

/**
 * integer_string:
 * @str: the string to examine.
 *
 * Returns: 1 if the given @str represents an integer, otherwise 0.
 */

int integer_string (const char *str)
{
    char *test;
    int ret = 1;

    if (str == NULL || *str == '\0') {
	return 0;
    }

    errno = 0;
    strtol(str, &test, 10);
    if (*test != '\0' || errno != 0) {
	ret = 0;
    }

    return ret;
}

/**
 * ends_with_backslash:
 * @s: the string to examine.
 *
 * Returns: 1 if the last non-space character in @s is a backslash,
 * otherwise 0.
 */

int ends_with_backslash (const char *s)
{
    int i, n = strlen(s);
    int bs = 0;

    for (i=n-1; i>=0; i--) {
	if (!isspace((unsigned char) s[i])) {
	    if (s[i] == '\\') {
		bs = 1;
	    }
	    break;
	}
    }

    return bs;
}

/**
 * gretl_lower:
 * @str: the string to transform.
 *
 * Converts any upper case characters in @str to lower case.
 *
 * Returns: the possibly modified string.
 */

char *gretl_lower (char *str)
{
    char *p = str;

    while (*p) {
        if (isupper((unsigned char) *p)) {
	    *p = tolower(*p);
	}
        p++;
    }

    return str;
}

/**
 * gretl_strdup:
 * @src: the string to duplicate.
 *
 * Returns: an allocated copy of @src, or NULL on error.
 */

char *gretl_strdup (const char *src)
{
    char *targ = NULL;

    if (src != NULL) {
	size_t n = strlen(src) + 1;

	targ = calloc(n, 1);
	if (targ != NULL) {
	    memcpy(targ, src, n);
	}
    }

    return targ;
}

/**
 * gretl_strdup_plus:
 * @src: the string to duplicate.
 * @plus: number of extra bytes to allocate.
 *
 * Returns: an allocated copy of @src, with extra trailing bytes
 * for possible appending, or NULL on error.
 */

char *gretl_strdup_plus (const char *src, int plus)
{
    char *targ = NULL;

    /* sanity guard */
    if (plus < 0) {
	plus = 0;
    }

    if (src != NULL) {
	size_t n = strlen(src) + 1 + plus;

	targ = calloc(n, 1);
	if (targ != NULL) {
	    memcpy(targ, src, n);
	}
    }

    return targ;
}

/**
 * gretl_strndup:
 * @src: the string to be copied.
 * @n: the maximum number of characters to copy.
 *
 * Returns: an allocated copy of at most @n characters from
 * @src, or NULL on error.
 */

char *gretl_strndup (const char *src, size_t n)
{
    char *targ = NULL;

    if (src != NULL) {
	size_t len = strlen(src);

	if (len > n) {
	    len = n;
	}

	targ = malloc(len + 1);
	if (targ != NULL) {
	    memcpy(targ, src, len);
	    targ[len] = '\0';
	}
    }

    return targ;
}

/**
 * gretl_strdup_printf:
 * @format: as in printf().
 * @Varargs: arguments to be printed.
 *
 * Print the arguments according to @format.
 *
 * Returns: allocated result of the printing, or NULL on failure.
 */

char *gretl_strdup_printf (const char *format, ...)
{
    va_list args;
    char *buf = NULL;
    int len;

#ifdef HAVE_VASPRINTF
    va_start(args, format);
    len = vasprintf(&buf, format, args);
    va_end(args);
    if (len < 0) {
	buf = NULL;
    }
#else
    int bsize = 2048;

    buf = malloc(bsize);
    if (buf == NULL) {
	return NULL;
    }

    memset(buf, 0, 1);

    va_start(args, format);
    len = vsnprintf(buf, bsize, format, args);
    va_end(args);

    if (len >= bsize) {
	fputs("gretl_strdup_printf warning: string was truncated\n",
	      stderr);
    }
#endif

    return buf;
}

/**
 * gretl_str_expand:
 * @orig: pointer to the base string.
 * @add: the string to be added.
 * @sep: string to be interpolated, or NULL.
 *
 * Creates a newly allocated string built by concatenating
 * @orig and @add, with @sep interpolated unless @sep is
 * NULL, and replaces the content of @orig with the new string.
 * As a special case, if @orig is NULL, or if the content of
 * @orig is NULL, we just duplicate @add.
 *
 * Returns: the reallocated string, or NULL on failure.  In case
 * of failure the content of @orig is freed, if @orig is not NULL,
 * to avoid memory leakage.
 */

char *gretl_str_expand (char **orig, const char *add, const char *sep)
{
    char *targ;
    int n;

    if (add == NULL) {
	return NULL;
    }

    if (orig == NULL || *orig == NULL) {
	return gretl_strdup(add);
    }

    n = strlen(*orig);
    if (sep != NULL) {
	n += strlen(sep);
    }
    n += strlen(add) + 1;

    targ = realloc(*orig, n);
    if (targ == NULL) {
	free(*orig);
	*orig = NULL;
	return NULL;
    }

    if (sep != NULL) {
	strcat(targ, sep);
    }
    strcat(targ, add);
    *orig = targ;

    return targ;
}

static int is_word_char (const char *s, int i, gretlopt opt)
{
    unsigned char c = *s;

    if (isalnum(c) || c == '_') {
	return 1;
    } else if ((opt & OPT_D) && c == '.') {
	return 1;
    } else if (opt & OPT_U) {
	/* allow Greeks */
	if ((c == 0xCE || c == 0xCF) && *(s+1)) {
	    return 1;
	} else if (i > 0) {
	    unsigned char prev = *(s-1);

	    return (prev == 0xCE && c >= 0x91 && c <= 0xBF) ||
		(prev == 0xCF && c >= 0x80 && c <= 0x89);
	}
    }

    return 0;
}

/**
 * gretl_word_strdup:
 * @src: the source string.
 * @ptr: location to receive end of word pointer, or NULL.
 * @opt: can include OPT_S for "strict" operation: in this
 * case an error is flagged if @src contains any characters
 * other than 'word' characters (see below), comma and space.
 * Also may include OPT_D to allow dot as a "word" character,
 * OPT_U to accept UTF-8 Greek letters.
 * @err: location to receive error code.
 *
 * Copies the first 'word' found in @src, where a word
 * is defined as consisting of alphanumeric characters
 * and the underscore.  If @ptr is not NULL, on exit it
 * points at the next position in @src after the copied
 * word.
 *
 * Returns: the allocated word or NULL in case no word is
 * found, or on error.
 */

char *gretl_word_strdup (const char *src, const char **ptr,
			 gretlopt opt, int *err)
{
    char *targ = NULL;

    if (src == NULL) {
	if (ptr != NULL) {
	    *ptr = NULL;
	}
    } else if (*src == '\0') {
	if (ptr != NULL) {
	    *ptr = src;
	}
    } else {
	const char *p;
	int len = 0;
	int i = 0;

	if (opt & OPT_S) {
	    /* strict: check for any junk */
	    while (*src && (*src == ' ' || *src == ',')) {
		src++; i++;
	    }
	    if (*src && !is_word_char(src, 0, opt)) {
		gretl_errmsg_sprintf(_("Unexpected symbol '%c'"), *src);
		*err = E_PARSE;
		return NULL;
	    }
	} else {
	    /* just skip to first 'word char' */
	    while (*src && !is_word_char(src, i, opt)) {
		src++; i++;
	    }
	}

	/* set to start of 'word' */
	p = src;

	while (is_word_char(src, i, opt)) {
	    len++;
	    src++;
	    i++;
	}

	if (opt & OPT_S) {
	    if (*src != '\0' && *src != ' ' && *src != ',') {
		gretl_errmsg_sprintf(_("Unexpected symbol '%c'"), *src);
		*err = E_PARSE;
		return NULL;
	    }
	}

	if (ptr != NULL) {
	    *ptr = src;
	}

	if (len > 0) {
	    targ = gretl_strndup(p, len);
	    if (targ == NULL) {
		*err = E_ALLOC;
	    }
	}
    }

    return targ;
}

/**
 * gretl_quoted_string_strdup:
 * @s: the source string.
 * @ptr: location to receive end pointer, or NULL.
 *
 * If @s starts with a quote (double or single), return a copy of
 * the portion of @s that is enclosed in quotes.  That is,
 * from @s + 1 up to but not including the next matching quote.
 * If @ptr is not NULL, on output it receives a pointer to
 * the next byte in @s after the closing quote.
 *
 * Returns: the allocated string or NULL on failure.
 */

char *gretl_quoted_string_strdup (const char *s, const char **ptr)
{
    char q, *ret = NULL;
    const char *p = NULL;

    if (s != NULL && (*s == '"' || *s == '\'')) {
	int gotit = 0;

	q = *s;
	s++;
	p = s;
	while (*p && !gotit) {
	    if (*p == q && *(p-1) != '\\') {
		/* found non-escaped matching quote */
		gotit = 1;
	    } else {
		p++;
	    }
	}
	if (!gotit) {
	    p = NULL;
	}
    }

    if (p == NULL) {
	if (ptr != NULL) {
	    *ptr = NULL;
	}
    } else {
	if (ptr != NULL) {
	    *ptr = p + 1;
	}
	ret = gretl_strndup(s, p - s);
    }

    return ret;
}

/* variant of gretl_string_split() that respects
   empty fields, including them in the output array
*/

static char **string_split_2 (const char *s, int *n,
			      const char *sep)
{
    char **S = NULL;
    gchar **tmp;
    gchar *mysep;
    int i, m = 0;

    *n = 0;

    if (sep[0] == '\t' && sep[1] == '\0') {
	mysep = g_strdup(sep);
    } else {
	mysep = g_strstrip(g_strdup(sep));
    }

    tmp = g_strsplit(s, mysep, -1);
    if (tmp != NULL) {
	for (i=0; tmp[i]; i++) {
	    m++;
	}
	if (m > 0) {
	    S = strings_array_new(m);
	    if (S != NULL) {
		for (i=0; i<m; i++) {
		    S[i] = gretl_strdup(g_strstrip(tmp[i]));
		}
	    }
	}
	g_strfreev(tmp);
    }

    g_free(mysep);
    *n = m;

    return S;
}

/* Re. the separator given to gretl_string_split():
   if it contains anything other than whitespace
   characters we'll respect empty fields
*/

static int respect_empty_fields (const char *s)
{
    if (s[0] == '\t' && s[1] == '\0') {
	/* treat single tab as "true" separator string */
	return 1;
    }
    while (*s) {
	if (!isspace(*s)) {
	    return 1;
	}
	s++;
    }

    return 0;
}

/**
 * gretl_string_split:
 * @s: the source string.
 * @n: location to receive the number of substrings.
 * @sep: string containing the character(s) to count as
 * field separators, or NULL. If @sep is NULL only the
 * space character counts.
 *
 * Parses @s into a set of zero or more substrings and
 * creates an array of those substrings. On sucessful exit
 * @n holds the number of substrings.
 *
 * Returns: the allocated array or NULL in case of failure.
 */

char **gretl_string_split (const char *s, int *n,
			   const char *sep)
{
    int i, k, m;
    char *word;
    char **S;

    *n = 0;
    if (s == NULL) {
	return NULL;
    }

    if (sep == NULL) {
	sep = " ";
    } else if (respect_empty_fields(sep)) {
	return string_split_2(s, n, sep);
    }

    m = count_fields(s, sep);
    if (m == 0) {
	return NULL;
    }

    S = strings_array_new(m);
    if (S == NULL) {
	return NULL;
    }

    for (i=0; i<m; i++) {
	s += strspn(s, sep);
	k = strcspn(s, sep);
	word = gretl_strndup(s, k);
	if (word == NULL) {
	    strings_array_free(S, m);
	    return NULL;
	}
	S[i] = word;
	s += k;
    }

    *n = m;

    return S;
}

/**
 * gretl_string_split_lines:
 * @s: the source string.
 * @n: location to receive the number of substrings.
 *
 * Parses @s into a set of zero or more substrings, one per
 * complete line of @s, and creates an array of those substrings.
 * On sucessful exit @n holds the number of substrings.
 *
 * Returns: the allocated array or NULL in case of failure.
 */

char **gretl_string_split_lines (const char *s, int *n)
{
    const char *p = s;
    int i, len, m = 0;
    int err = 0;
    char **S = NULL;

    *n = 0;

    while (*p) {
	if (*p == '\n') {
	    m++;
	}
	p++;
    }

    if (m == 0) {
	return NULL;
    }

    S = strings_array_new(m);
    if (S == NULL) {
	return NULL;
    }

    p = s;
    i = 0;

    while (*p && i < m) {
	len = strcspn(p, "\r\n");
	S[i] = gretl_strndup(p, len);
	if (S[i] == NULL) {
	    err = E_ALLOC;
	    break;
	}
	gretl_strstrip(S[i]);
	i++;
	p += len;
	if (*p == '\r') p++;
	if (*p == '\n') p++;
	if (*p == '\0') {
	    break;
	}
    }

    if (err) {
	strings_array_free(S, m);
	S = NULL;
    } else {
	*n = m;
    }

    return S;
}

/**
 * gretl_string_split_quoted:
 * @s: the source string.
 * @n: location to receive the number of substrings.
 * @sep: string containing the character(s) to count as
 * field separators, or NULL. If @sep is NULL only space,
 * tab and newline count.
 * @err: location to receive error code.
 *
 * Similar to gretl_string_split(), except that this variant
 * allows for the presence of double-quoted substrings
 * which may contain spaces. The quotes are removed in the
 * members of the returned array.
 *
 * Returns: allocated array of substrings or NULL in case of failure.
 */

char **gretl_string_split_quoted (const char *s, int *n,
				  const char *sep, int *err)
{
    const char *ignore;
    const char *q, *p = s;
    int i, len, m = 0;
    int grabit, quoted;
    char *substr;
    char **S;

    *err = 0;
    ignore = sep != NULL ? sep : " \t\n";

    *n = 0;

    while (*p) {
	p += strspn(p, ignore);
	if (*p == '"') {
	    /* quoted substring */
	    m++;
	    q = strchr(p + 1, '"');
	    if (q == NULL) {
		*err = E_PARSE;
		return NULL;
	    }
	    p = q;
	} else {
	    len = strcspn(p, ignore);
	    if (len > 0) {
		/* unquoted substring */
		m++;
		p += len - 1;
	    }
	}
	if (*p == '\0') {
	    break;
	}
	p++;
    }

    if (*err || m == 0) {
	return NULL;
    }

    S = strings_array_new(m);
    if (S == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    p = s;
    i = 0;

    while (*p && i < m) {
	grabit = quoted = 0;
	p += strspn(p, ignore);
	if (*p == '"') {
	    grabit = quoted = 1;
	    p++;
	    len = strcspn(p, "\"");
	} else {
	    len = strcspn(p, ignore);
	    grabit = (len > 0);
	}
	if (grabit) {
	    substr = gretl_strndup(p, len);
	    if (substr == NULL) {
		*err = E_ALLOC;
		strings_array_free(S, m);
		return NULL;
	    }
	    S[i++] = substr;
	    p += len + quoted;
	}
    }

    *n = m;

    return S;
}

/**
 * gretl_trunc:
 * @str: the string to truncate.
 * @n: the desired length of the truncated string.
 *
 * Truncates the given @str to the specified length.
 *
 * Returns: the possibly truncated string.
 */

char *gretl_trunc (char *str, size_t n)
{
    if (n < strlen(str)) {
	str[n] = '\0';
    }

    return str;
}

static const char *name_ok =
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "0123456789_";

/**
 * is_greek_letter:
 * @s: the string to test.
 *
 * Checks for the case where @s contains a greek letter,
 * followed by a nul byte or an ascii character that cannot
 * occur in a gretl identifier.

 * Returns: 1 if the case is met, otherwise 0.
 */

int is_greek_letter (const char *s)
{
    unsigned char u[2] = {0};
    int ret = 0;

    u[0] = s[0];
    if (s[0] != '\0') {
	u[1] = s[1];
    }

    if ((u[0] == 0xCE && u[1] >= 0x91 && u[1] <= 0xBF) ||
	(u[0] == 0xCF && u[1] >= 0x80 && u[1] <= 0x89)) {
	char nxt = s[2];

	if (nxt == '\0') {
	    ret = 1;
	} else if (nxt >= 32 && nxt <= 126 && !strchr(name_ok, nxt)) {
	    ret = 1;
	}
    }

    return ret;
}

/**
 * gretl_namechar_spn:
 * @s: the string to examine.
 *
 * Returns: the length of the initial segment of @s which
 * consists of characters that are valid in a gretl
 * variable or object name, namely a-z, A-Z, 0-9 and _,
 * starting with a letter, up to a maximum of %VNAMELEN - 1.
 */

int gretl_namechar_spn (const char *s)
{
    int ret = 0;

    if (isalpha(*s)) {
	ret = strspn(s, name_ok);
    }

    if (ret == 0 && is_greek_letter(s)) {
	ret = 2;
    }

    return ret;
}

/**
 * double_quote_position:
 * @s: the source string.
 *
 * Returns: the 0-based index of the position of the next
 * unescaped double-quote character in @s, or -1 if no
 * such character is found.
 */

int double_quote_position (const char *s)
{
    int nbs; /* number of preceding backslashes */
    int i, j;
    int ret = -1;

    for (i=0; s[i]; i++) {
	if (s[i] == '"') {
	    nbs = 0;
	    for (j=i-1; j>=0; j--) {
		if (s[j] == '\\') {
		    nbs++;
		} else {
		    break;
		}
	    }
	    if (nbs % 2 == 0) {
		/* got an unescaped double-quote */
		ret = i;
		break;
	    }
	}
    }

    return ret;
}

/**
 * count_fields:
 * @s: the string to process.
 * @sep: string containing the character(s) to count as
 * field separators, or NULL. If @sep is NULL only the
 * space character counts.
 *
 * Returns: the number of fields in @s.
 */

int count_fields (const char *s, const char *sep)
{
    int nf = 0;

    if (sep == NULL) {
	sep = " ";
    }

    if (s != NULL && *s != '\0') {
	const char *p;

	/* step past separator(s) */
	s += strspn(s, sep);

	if (*s != '\0') {
	    s++;
	    nf++;
	}

	while (*s) {
	    p = strpbrk(s, sep);
	    if (p != NULL) {
		s = p + strspn(p, sep);
		if (*s) {
		    nf++;
		}
	    } else {
		break;
	    }
	}
    }

    return nf;
}

/**
 * count_lines:
 * @s: the string to process.
 *
 * Returns: the number of complete lines (lines ending
 * with the newline character) in @s.
 */

int count_lines (const char *s)
{
    int nl = 0;

    if (s != NULL) {
	while (*s) {
	    if (*s == '\n') {
		nl++;
	    }
	    s++;
	}
    }

    return nl;
}

/**
 * shift_string_left:
 * @str: the string to process.
 * @move: the number of places to shift.
 *
 * Shifts the content of @str left by @move places, dropping
 * leading bytes as needed.
 *
 * Returns: the modified string.
 */

char *shift_string_left (char *str, size_t move)
{
    size_t n = strlen(str);

    if (move >= n) {
	*str = '\0';
    } else {
	memmove(str, str + move, n - move);
	str[n - move] = '\0';
    }

    return str;
}

/**
 * gretl_strstrip:
 * @str: the string to process.
 *
 * Removes leading and trailing white space from a string.
 *
 * Returns: the possibly modified string.
 */

char *gretl_strstrip (char *str)
{
    int i, n = strspn(str, " \t");

    if (n > 0) {
	shift_string_left(str, n);
    }

    n = strlen(str);

    for (i=n-1; i>=0; i--) {
	if (isspace(str[i]) || str[i] == '\r') {
	    str[i] = '\0';
	} else {
	    break;
	}
    }

    return str;
}

/**
 * gretl_strstrip_copy:
 * @str: the string to process.
 *
 * Returns: a copy of @str, from which both leading and
 * trailing white space have been removed.
 */

char *gretl_strstrip_copy (const char *str, int *err)
{
    char *ret = NULL;
    int i, n;

    if (str == NULL) {
	*err = E_INVARG;
	return NULL;
    }

    while (isspace(*str)) {
	str++;
    }

    n = strlen(str);

    for (i=n-1; i>=0; i--) {
	if (isspace(str[i]) || str[i] == '\r') {
	    n--;
	} else {
	    break;
	}
    }

    ret = gretl_strndup(str, n);
    if (ret == NULL) {
	*err = E_ALLOC;
    }

    return ret;
}

/**
 * switch_ext:
 * @targ: the target or output string (must be pre-allocated).
 * @src: the source or input string.
 * @ext: the extension or suffix to attach, without leading dot.
 *
 * For processing filenames: copies @src to @targ, minus any existing
 * filename extension, and adds to @targ the specified extension.
 *
 * Returns: the output string, @targ.
 */

char *switch_ext (char *targ, const char *src, const char *ext)
{
    int i = gretl_dotpos(src);

    if (targ != src) {
	*targ = '\0';
        strncat(targ, src, i);
    }

    targ[i] = '.';
    targ[i + 1] = '\0';
    strcat(targ, ext);

    return targ;
}

/**
 * switch_ext_in_place:
 * @fname: must have sufficient space to add the given extension.
 * @ext: the extension or suffix to attach, without dot.
 *
 * For processing filenames: removes any existing dot-extension on
 * @fname and appends a dot followed by @ext.
 *
 * Returns: the modified string, @fname.
 */

char *switch_ext_in_place (char *fname, const char *ext)
{
    int i = gretl_dotpos(fname);

    fname[i] = '.';
    fname[i + 1] = '\0';
    strcat(fname, ext);

    return fname;
}

/**
 * switch_ext_new:
 * @src: the original string.
 * @ext: the extension or suffix to attach (without leading '.').
 *
 * For processing filenames: creates a copy of @src in which
 * any existing dot-extension is removed and @ext is appended
 * (with a dot automatically inserted).
 *
 * Returns: the newly allocated string.
 */

char *switch_ext_new (const char *src, const char *ext)
{
    int len = strlen(src) + strlen(ext) + 2;
    const char *p = strrchr(src, '.');
    char *ret = NULL;

    if (p != NULL) {
	len -= strlen(p);
    }

    ret = calloc(len, 1);

    if (ret != NULL) {
	if (p != NULL) {
	    strncat(ret, src, p - src);
	} else {
	    strcat(ret, src);
	}
	strcat(ret, ".");
	strcat(ret, ext);
    }

    return ret;
}

static int ends_in_comment (const char *s, int n)
{
    int i, quoted = 0;

    /* the '#' character is inert (only) if it appears
       within a string literal */

    for (i=n; i>1; i--) {
	if (s[i] == '"') {
	    quoted = !quoted;
	} else if (!quoted && s[i] == '#') {
	    return 1;
	}
    }

    return 0;
}

#define LINE_CONT(c) (c == '\\' || c == ',' || c == '(')

/**
 * top_n_tail:
 * @str: the string to process.
 * @maxlen: maximum length of string, including NUL termination.
 * @err: location to receive error code, or NULL.
 *
 * Drop leading space and trailing space and newline from string,
 * then replace a trailing backslash (if any) with a space.
 * If @str does not end with a newline within the limit set by
 * @maxlen, and @err is not NULL, then E_TOOLONG is written
 * to @err.
 *
 * Returns: 1 if a trailing backslash, comma, semicolon, or left
 * parenthesis was found, otherwise 0.
 */

int top_n_tail (char *str, size_t maxlen, int *err)
{
    int i, n, cont = 0;

    if (str == NULL || *str == '\0' || *str == '\n' || *str == '\r') {
	return 0;
    }

    n = strlen(str) - 1;

    if (err != NULL && n > maxlen - 2 && str[n] != '\n') {
	*err = E_TOOLONG;
    }

    /* chop any trailing space */
    for (i=n; i>=0; i--) {
	if (isspace((unsigned char) str[i])) {
	    str[i] = '\0';
	    n--;
	} else {
	    break;
	}
    }

    if (*str != '\0') {
	/* Drop any leading spaces, also possible questionmark.  Try
	   to catch non-breaking spaces too -- ugh, Windows!
	   (NBSP is 0xA0 in Windows CP1252)
	*/
	i = 0;
	while (isspace((unsigned char) str[i]) ||
	       str[i] == '?' ||
	       str[i] == (char) 0xC2 ||
	       str[i] == (char) 0xA0) {
	    n--;
	    i++;
	}
	if (i > 0) {
	    shift_string_left(str, i);
	}

	if (*str == '#' || !strncmp(str, "/*", 2)) {
	    ; /* the line starts a comment: leave well alone */
	} else if (n >= 0 && LINE_CONT(str[n])) {
	    /* register line continuation characters at the end of
	       the line, but only if not preceded by the comment
	       character '#' (unquoted)
	    */
	    cont = !ends_in_comment(str, n - 1);
	    if (cont && str[n] == '\\') {
		/* replace backslash */
		str[n] = ' ';
	    }
	}
    }

    return cont;
}

/**
 * equation_get_lhs_and_rhs:
 * @s: equation in string form.
 * @plh: pointer to receive left-hand side expression.
 * @prh: pointer to receive right-hand side expression.
 *
 * Given a string @s, parse it into a left-hand side and a right-hand
 * side, separated by an equals sign.  Return in @plh and @prh
 * allocated copies of the respective sides, with any leading or trailing
 * white space trimmed.
 *
 * Returns: 0 on success, 1 on error.
 */

int equation_get_lhs_and_rhs (const char *s, char **plh, char **prh)
{
    const char *p;
    char *lh = NULL, *rh = NULL;
    int i, len, err = 0;

    if (s == NULL || plh == NULL || prh == NULL) {
	err = 1;
    }

    if (!err) {
	*plh = NULL;
	*prh = NULL;

	p = strchr(s, '=');
	if (p == NULL) {
	    err = 1;
	}
    }

    if (!err) {
	p = s;
	while (isspace(*p)) p++;
	len = strcspn(p, " =");
	if (len == 0) {
	    err = 1;
	} else {
	    lh = gretl_strndup(p, len);
	    if (lh == NULL) {
		err = 1;
	    }
	}
    }

    if (!err) {
	p = strchr(s, '=') + 1;
	while (isspace(*p)) p++;
	len = strlen(p);
	if (len == 0) {
	    err = 1;
	} else {
	    for (i=len-1; i>=0; i--) {
		if (isspace(p[i])) len--;
		else break;
	    }
	    rh = gretl_strndup(p, len);
	    if (rh == NULL) {
		err = 1;
	    }
	}
    }

    if (err) {
	free(lh);
	free(rh);
    } else {
	*plh = lh;
	*prh = rh;
    }

    return err;
}

/**
 * tailstrip:
 * @str: the string to process.
 *
 * Drop trailing space (and newline if any) from string.
 *
 * Returns: the modified string.
 */

char *tailstrip (char *str)
{
    int i, len;

    if (str == NULL || *str == '\0') {
	return str;
    }

    len = strlen(str);

    for (i=len-1; i>=0; i--) {
	if (isspace((unsigned char) str[i]) ||
	    str[i] == '\n' || str[i] == '\r') {
	    str[i] = '\0';
	} else {
	    break;
	}
    }

    return str;
}

/**
 * compress_spaces:
 * @s: the string to process.
 *
 * Reduce multiple contiguous space characters to single spaces
 * within @s.
 *
 * Returns: the compressed string.
 */

char *compress_spaces (char *s)
{
    int i = 0, inquote = 0;
    char *p, *q;

    if (s == NULL || *s == '\0') {
	return s;
    }

    p = q = s;

    while (*s) {
	if (*s == '"' && (i == 0 || *(s-1) != '\\')) {
	    inquote = !inquote;
	}
	if (!inquote) {
	    if (*s == '\t') {
		*s = ' '; /* trash tabs */
	    }
	    if (*s == ' ') {
		p = s + 1;
		if (*p == '\0') break;
		while (*p == ' ') p++;
		if (p - s > 1) {
		    memmove(s + 1, p, strlen(p) + 1);
		}
	    }
	}
	s++;
	i++;
    }

    return q;
}

/**
 * space_to_score:
 * @s: the string to process.
 *
 * Replace any spaces with underscores in @s.
 *
 * Returns: the (possibly) modified string.
 */

char *space_to_score (char *s)
{
    char *p = s;

    while (*p) {
	if (*p == ' ') *p = '_';
	p++;
    }

    return s;
}

/**
 * strings_array_new:
 * @nstrs: number of strings in array.
 *
 * Allocates storage for @nstrs strings and initializes all
 * to NULL.
 *
 * Returns: the allocated array, or NULL on failure.
 */

char **strings_array_new (int nstrs)
{
    char **s;
    int i;

    if (nstrs <= 0) {
	return NULL;
    }

    s = malloc(nstrs * sizeof *s);
    if (s != NULL) {
	for (i=0; i<nstrs; i++) {
	    s[i] = NULL;
	}
    }

    return s;
}

/**
 * strings_array_add:
 * @pS: pointer to strings array.
 * @n: location of present number of strings in array.
 * @p: string to add to array.
 *
 * Allocates storage for an extra member of @S and adds a
 * copy of string @p in the last position.  On success,
 * the content of @n is incremented by 1.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int strings_array_add (char ***pS, int *n, const char *p)
{
    char **Tmp;
    int m = *n;

    Tmp = realloc(*pS, (m + 1) * sizeof *Tmp);
    if (Tmp == NULL) {
	return E_ALLOC;
    }

    *pS = Tmp;

    if (p != NULL) {
	Tmp[m] = gretl_strdup(p);
	if (Tmp[m] == NULL) {
	    return E_ALLOC;
	}
    } else {
	Tmp[m] = NULL;
    }

    *n += 1;

    return 0;
}

/**
 * strings_array_prepend_uniq:
 * @pS: pointer to strings array.
 * @n: location of present number of strings in array.
 * @p: string to prepend to array.
 *
 * If @p is already present in the array at location @pS,
 * moves it into first position if it is not already there.
 * Otherwise allocates storage for an extra member of the
 * array and pushes @p into first position, in which case
 * the content of @n is incremented by 1.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int strings_array_prepend_uniq (char ***pS, int *n, const char *p)
{
    char **S = *pS;
    char *s0;
    int m = *n;
    int i, pos = -1;

    if (p == NULL) {
	return E_DATA;
    }

    for (i=0; i<m; i++) {
	if (!strcmp(S[i], p)) {
	    pos = i;
	    break;
	}
    }

    if (pos == 0) {
	/* already present in 1st position */
	return 0;
    } else if (pos > 0) {
	/* already present, not first */
	s0 = S[pos];
	for (i=pos; i>0; i--) {
	    S[i] = S[i-1];
	}
	S[0] = s0;
    } else if (pos < 0) {
	/* not present */
	S = realloc(*pS, (m + 1) * sizeof *S);
	if (S == NULL) {
	    return E_ALLOC;
	}
	s0 = gretl_strdup(p);
	if (s0 == NULL) {
	    return E_ALLOC;
	}
	*pS = S;
	*n = ++m;
	/* shuffle up */
	for (i=m; i>0; i--) {
	    S[i] = S[i-1];
	}
	/* and insert @p */
	S[0] = s0;
    }

    return 0;
}

/**
 * strings_array_donate:
 * @pS: pointer to strings array.
 * @n: location of present number of strings in array.
 * @p: string to append to array.
 *
 * Allocates storage for an extra member of @S and adds
 * string @p in the last position. Unlike strings_array_add(),
 * the array takes ownnership of @p rather than copying it.
 * On success, the content of @n is incremented by 1.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int strings_array_donate (char ***pS, int *n, char *p)
{
    char **Tmp;
    int m = *n;

    Tmp = realloc(*pS, (m + 1) * sizeof *Tmp);
    if (Tmp == NULL) {
	return E_ALLOC;
    }

    *pS = Tmp;

    if (p != NULL) {
	Tmp[m] = p;
	if (Tmp[m] == NULL) {
	    return E_ALLOC;
	}
    } else {
	Tmp[m] = NULL;
    }

    *n += 1;

    return 0;
}

/**
 * strings_array_add_uniq:
 * @pS: pointer to strings array.
 * @n: location of present number of strings in array.
 * @p: string to test for addition to array.
 * @pos: location to receive the position of @p in the
 * array (whether already matched or newly added), or NULL.
 *
 * If the array does not already include a copy of @p,
 * allocates storage for an extra member of @pS and adds a
 * copy of string @p in the last position. On successful
 * addition the content of @n is incremented by 1.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int strings_array_add_uniq (char ***pS, int *n, const char *p,
			    int *pos)
{
    char **Tmp, **S = *pS;
    int m = *n;
    int i;

    for (i=0; i<m; i++) {
	if (S[i] != NULL && strcmp(S[i], p) == 0) {
	    if (pos != NULL) {
		*pos = i;
	    }
	    return 0; /* no-op */
	}
    }

    Tmp = realloc(S, (m + 1) * sizeof *Tmp);
    if (Tmp == NULL) {
	return E_ALLOC;
    }

    if (p != NULL) {
	Tmp[m] = gretl_strdup(p);
	if (Tmp[m] == NULL) {
	    return E_ALLOC;
	}
    } else {
	Tmp[m] = NULL;
    }

    if (pos != NULL) {
	*pos = m;
    }

    *pS = Tmp;
    *n += 1;

    return 0;
}

/**
 * strings_array_new_with_length:
 * @nstrs: number of strings in array.
 * @len: number of bytes per string.
 *
 * Allocates storage for @nstrs strings, each of them
 * @len bytes long.  The first byte of each string is
 * initialized to 0.
 *
 * Returns: the allocated array, or NULL on failure.
 */

char **strings_array_new_with_length (int nstrs, int len)
{
    char **S;
    int i, j;

    if (nstrs <= 0) {
	return NULL;
    }

    S = malloc(nstrs * sizeof *S);
    if (S == NULL) return NULL;

    for (i=0; i<nstrs; i++) {
	S[i] = malloc(len);
	if (S[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(S[j]);
	    }
	    free(S);
	    return NULL;
	}
	S[i][0] = '\0';
    }

    return S;
}

/**
 * strings_array_realloc_with_length:
 * @pS: existing array to reallocate.
 * @oldn: original number of strings in the array.
 * @newn: new number of strings in array.
 * @len: number of bytes per string (or 0).
 *
 * Adjusts the storage in @pS to a size of @newn
 * strings, each of them @len bytes long.  If @len
 * is greater than zero the first byte of any
 * additional strings is initialized to 0.
 * This function may be used either to expand or to
 * shrink an existing array of strings.
 *
 * Returns: the new array, or NULL on failure.
 */

char **strings_array_realloc_with_length (char ***pS,
					  int oldn,
					  int newn,
					  int len)
{
    char **S;
    int i;

    if (pS == NULL) {
	return NULL;
    } else if (newn == oldn) {
	return *pS;
    } else if (newn <= 0) {
	strings_array_free(*pS, oldn);
	*pS = NULL;
	return NULL;
    }

    /* in case we're shrinking the array */
    for (i=newn; i<oldn; i++) {
	free((*pS)[i]);
	(*pS)[i] = NULL;
    }

    S = realloc(*pS, newn * sizeof *S);
    if (S == NULL) {
	strings_array_free(*pS, oldn);
	*pS = NULL;
	return NULL;
    }

    *pS = S;

    /* in case we're expanding the array */
    if (len == 0) {
	for (i=oldn; i<newn; i++) {
	    S[i] = NULL;
	}
    } else {
	int j;

	for (i=oldn; i<newn; i++) {
	    S[i] = malloc(len);
	    if (S[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(S[j]);
		}
		free(*pS);
		*pS = NULL;
		break;
	    } else {
		S[i][0] = '\0';
	    }
	}
    }

    return *pS;
}

/**
 * strings_array_dup:
 * @strs: array of strings to be copied.
 * @n: number of strings in array.
 *
 * Returns: an allocated copy of @strs, or NULL on failure.
 */

char **strings_array_dup (char **strs, int n)
{
    char **S = NULL;
    int i, err = 0;

    if (n <= 0 || strs == NULL) {
	return NULL;
    }

    S = strings_array_new(n);
    if (S == NULL) return NULL;

    for (i=0; i<n; i++) {
	if (strs[i] == NULL) {
	    S[i] = NULL;
	} else {
	    S[i] = gretl_strdup(strs[i]);
	    if (S[i] == NULL) {
		err = E_ALLOC;
		break;
	    }
	}
    }

    if (err) {
	strings_array_free(S, n);
	S = NULL;
    }

    return S;
}

/**
 * strings_array_dup_selected:
 * @strs: array of strings from which to copy.
 * @n: number of strings in array.
 * @list: list of indices of strings to copy.
 *
 * Returns: an allocated array holding the strings selected
 * from @strs, or NULL on failure.
 */

char **strings_array_dup_selected (char **strs, int n,
				   const int *list)
{
    char **S = NULL;
    int m, i, j, k;
    int err = 0;

    if (n <= 0 || strs == NULL) {
	return NULL;
    } else if (list == NULL || list[0] == 0) {
	return NULL;
    }

    m = list[0];
    S = strings_array_new(m);
    if (S == NULL) return NULL;

    for (i=1, k=0; i<=m; i++) {
	j = list[i];
	if (j < 0 || j >= n) {
	    err = E_DATA;
	    break;
	}
	if (strs[j] == NULL) {
	    S[k] = NULL;
	} else {
	    S[k] = gretl_strdup(strs[j]);
	    if (S[k] == NULL) {
		err = E_ALLOC;
		break;
	    }
	}
	k++;
    }

    if (err) {
	strings_array_free(S, m);
	S = NULL;
    }

    return S;
}

/**
 * strings_array_sort:
 * @pS: location of array of strings.
 * @n: location of the number of strings in the array.
 * @opt: may contain %OPT_U to trim the sorted array
 * so that it contains only unique entries.
 *
 * Sorts an array of strings in ascending lexicographical
 * order. If %OPT_U is given, @n holds the number of unique
 * strings on exit. It is assumed that storage for the
 * strings array was obtained via strings_array_new() or
 * a similar libgretl function.
 *
 * Returns: 0 on success, non-zero on error.
 */

int strings_array_sort (char ***pS, int *n, gretlopt opt)
{
    char **S;
    int ns;

    if (pS == NULL || n == NULL) {
	return E_DATA;
    }

    S = *pS;
    ns = *n;

    qsort(S, ns, sizeof *S, gretl_compare_strings);

    if (opt & OPT_U) {
	int i, j, m = ns;

	for (i=0; i<m-1; i++) {
	    if (!strcmp(S[i], S[i+1])) {
		free(S[i+1]);
		for (j=i+1; j<m-1; j++) {
		    S[j] = S[j+1];
		}
		S[m-1] = NULL;
		i--;
		m--;
	    }
	}
	if (m < ns) {
	    char **tmp = realloc(S, m * sizeof *S);

	    if (tmp != NULL) {
		*pS = tmp;
	    }
	    *n = m;
	}
    }

    return 0;
}

/**
 * strings_array_cmp:
 * @strs1: first array of strings.
 * @strs2: second array of strings.
 * @n: number of strings to examine.
 *
 * Compares for equality two arrays of strings, each of
 * which must contain at least @n elements.  Equality
 * of the arrays means that strcmp returns 0 for
 * each pair of strings @strs1[i], @strs2[i], for i
 * equals 0 to @n - 1.
 *
 * Returns: 0 if the arrays compare equal, non-zero
 * otherwise.
 */

int strings_array_cmp (char **strs1, char **strs2, int n)
{
    int i, ret = 0;

    for (i=0; i<n && !ret; i++) {
	ret = strcmp(strs1[i], strs2[i]);
    }

    return ret;
}

/**
 * strings_array_position:
 * @strs: array of strings.
 * @n: number of elements in @strs.
 * @s: string to test.
 *
 * Returns: the 0-based position of the first member of @strs
 * to compare equal to @s, or -1 if no match is found.
 */

int strings_array_position (char **strs, int n, const char *s)
{
    int i, ret = -1;

    if (s != NULL) {
	for (i=0; i<n && ret<0; i++) {
	    if (strs[i] != NULL && !strcmp(strs[i], s)) {
		ret = i;
	    }
	}
    }

    return ret;
}

/**
 * strings_array_diff:
 * @strs1: first array of strings.
 * @n1: number of strings in @strs1.
 * @strs2: second array of strings.
 * @n2: number of strings in @strs2.
 * @extra: location to receive array of strings present
 * in @strs2 but not in @strs1 (or NULL).
 * @n_extra: location to receive bumber of strings in
 * the @extra array (or NULL).
 *
 * Fills out @extra with any strings present in @strs2
 * but not in @strs1.
 *
 * Returns: 0 on success, non-zero otherwise.
 */

int strings_array_diff (char **strs1, int n1,
			char **strs2, int n2,
			char ***extra, int *n_extra)
{
    char **S = NULL;
    int i, j, n, found;
    int err = 0;

    if (extra != NULL && n_extra != NULL) {
	n = 0;
	for (j=0; j<n2 && !err; j++) {
	    found = 0;
	    for (i=0; i<n1 && !found; i++) {
		if (!strcmp(strs2[j], strs1[i])) {
		    found = 1;
		}
	    }
	    if (!found) {
		err = strings_array_add(&S, &n, strs2[j]);
	    }
	}
	if (!err) {
	    *extra = S;
	    *n_extra = n;
	}
    }

    return err;
}

/**
 * strings_array_reverse:
 * @strs: array of allocated strings.
 * @nstrs: number of strings in array.
 *
 * Returns: an array of strings with @nstrs in which the
 * strings are those in @strs, in reverse order -- or
 * NULL on failure.
 */

char **strings_array_reverse (char **strs, int nstrs)
{
    char **S = NULL;
    int i;

    S = strings_array_new(nstrs);

    if (S != NULL) {
	for (i=0; i<nstrs; i++) {
	    S[i] = gretl_strdup(strs[nstrs-i-1]);
	    if (S[i] == NULL) {
		strings_array_free(S, nstrs);
		break;
	    }
	}
    }

    return S;
}

/**
 * strings_array_free:
 * @strs: array of allocated strings.
 * @nstrs: number of strings in array.
 *
 * Frees each allocated string in @strs, then frees @strs itself.
 * Checks that @strs is not NULL before proceeding.
 */

void strings_array_free (char **strs, int nstrs)
{
    int i;

    if (strs != NULL) {
	for (i=0; i<nstrs; i++) {
	    free(strs[i]);
	}
	free(strs);
    }
}

/**
 * get_obs_string:
 * @obs: char array big enough to hold the observation (#OBSLEN).
 * @t: zero-based observation number.
 * @dset: pointer to dataset information.
 *
 * Returns: the observation string corresponding to @t.
 */

char *get_obs_string (char *obs, int t, const DATASET *dset)
{
    if (dataset_has_markers(dset)) {
	strcpy(obs, dset->S[t]);
    } else {
	ntolabel(obs, t, dset);
    }

    return obs;
}

/**
 * obs_str_to_double:
 * @obs: string representation of observation number.
 *
 * Returns: the floating-point counterpart of @obs,
 * or #NADBL on invalid input.
 */

double obs_str_to_double (const char *obs)
{
    char *p, *test, tmp[OBSLEN];
    double ret;

    strcpy(tmp, obs);
    p = tmp;

    while (*p) {
	if (*p == ':' || *p == ',') {
	    *p = '.';
	}
	p++;
    }

    errno = 0;

    gretl_push_c_numeric_locale();
    ret = strtod(tmp, &test);
    gretl_pop_c_numeric_locale();

    if (*test != '\0' || errno == ERANGE) {
	ret = NADBL;
    }

    return ret;
}

/**
 * colonize_obs:
 * @obs: string representation of observation number.
 *
 * Converts a decimal point in @obs to a colon.
 *
 * Returns: the (possibly) modified obs string.
 */

char *colonize_obs (char *obs)
{
    char *p = obs;

    while (*p) {
	if (*p == '.' || *p == ',') {
	    *p = ':';
	}
	p++;
    }

    return obs;
}

/**
 * modify_obs_for_csv:
 * @s: observation string (date).
 * @pd: data frequency.
 *
 * Modifies the observation string corresponding to obervation @t to
 * producing a form suitable for a CSV file.  This applies only to
 * time series data. The string @s should be obtained by calling
 * ntolabel();
 */

void modify_date_for_csv (char *s, int pd)
{
    if (pd == 4) {
	gretl_charsub(s, ':', 'Q');
    } else {
	gretl_charsub(s, ':', 'M');
    }
}

/**
 * print_time:
 * @s: string into which to print: must be at least 48 bytes.
 *
 * Returns: @s, which will contain a string representation of the
 * current date and time, in the format YYYY-mm-dd H:M.
 */

char *print_time (char *s)
{
    time_t now = time(NULL);
    struct tm *local;

    local = localtime(&now);
    strftime(s, 47, "%Y-%m-%d %H:%M", local);

    return s;
}

/**
 * gretl_xml_validate:
 * @s: string to be tested.
 *
 * Returns: 1 if @s is acceptable for insertion into an XML file
 * as is, 0 if it contains special characters that need to be
 * escaped.  See also gretl_xml_encode().
 */

int gretl_xml_validate (const char *s)
{
    while (*s) {
	if (*s == '&' || *s == '<' || *s == '>' || *s == '"') {
	    return 0;
	}
	s++;
    }

    return 1;
}

/**
 * gretl_xml_encode:
 * @str: NUL-terminated source string.
 *
 * Returns: an allocated re-write of @str, with characters that are
 * special in XML encoded as character entities.  See also
 * gretl_xml_validate().
 */

char *gretl_xml_encode (const char *str)
{
    char *targ, *p;
    const char *s = str;
    int len = strlen(s) + 1;

    while (*s) {
	if (*s == '&') len += 4;
	else if (*s == '<') len += 3;
	else if (*s == '>') len += 3;
	else if (*s == '"') len += 5;
	s++;
    }

    targ = malloc(len);
    if (targ == NULL) {
	gretl_errmsg_set(_("out of memory in XML encoding"));
	return NULL;
    }

    s = str;
    p = targ;

    while (*s) {
	if (*s == '&') {
	    strcpy(p, "&amp;");
	    p += 5;
	} else if (*s == '<') {
	    strcpy(p, "&lt;");
	    p += 4;
	} else if (*s == '>') {
	    strcpy(p, "&gt;");
	    p += 4;
	} else if (*s == '"') {
	    strcpy(p, "&quot;");
	    p += 6;
	} else {
	    *p++ = *s;
	}
	s++;
    }

    targ[len-1] = '\0';

#ifdef XML_DEBUG
    fprintf(stderr, "done gretl_xml_encode: targ='%s'\n", targ);
#endif

    return targ;
}

/**
 * gretl_xml_encode_to_buf:
 * @targ: target buffer.
 * @src: NUL-terminated source string.
 * @n: size of @targ in bytes.
 *
 * Writes into @targ a version of @src in which characters that are
 * special in XML are encoded as character entities.  See also
 * gretl_xml_encode() for the case where the encoding of @src is
 * of unknown size at compile time.
 *
 * Returns: 0 on success or 1 on error.  An error occurs if (a) the
 * encoded version of @src is longer than @n bytes (allowing for NUL
 * termination), or (b) @src does not validate as UTF-8.  On error
 * the conversion is not done.
 */

int gretl_xml_encode_to_buf (char *targ, const char *src, int n)
{
    const char *s = src;
    int len = strlen(s) + 1;

    if (!g_utf8_validate(src, -1, NULL)) {
	fprintf(stderr, "gretl_xml_encode_to_buf: source not UTF-8\n");
	return 1;
    }

    while (*s) {
	if (*s == '&') len += 4;
	else if (*s == '<') len += 3;
	else if (*s == '>') len += 3;
	else if (*s == '"') len += 5;
	s++;
    }

    *targ = '\0';

    if (len > n) {
	fprintf(stderr, "gretl_xml_encode_to_buf: buffer too small\n");
	return 1;
    }

    s = src;

    while (*s) {
	if (*s == '&') {
	    strcpy(targ, "&amp;");
	    targ += 5;
	} else if (*s == '<') {
	    strcpy(targ, "&lt;");
	    targ += 4;
	} else if (*s == '>') {
	    strcpy(targ, "&gt;");
	    targ += 4;
	} else if (*s == '"') {
	    strcpy(targ, "&quot;");
	    targ += 6;
	} else {
	    *targ++ = *s;
	}
	s++;
    }

    *targ = '\0';

    return 0;
}

static char x2c (char *s)
{
    register char digit;

    digit = (s[0] >= 'A' ? ((s[0] & 0xdf) - 'A') + 10 : (s[0] - '0'));
    digit *= 16;
    digit += (s[1] >= 'A' ? ((s[1] & 0xdf) - 'A') + 10 : (s[1] - '0'));
    return digit;
}

/**
 * unescape_url:
 * @url: string representing a URL.
 *
 */

void unescape_url (char *url)
{
    register int x, y;

    for (x=0, y=0; url[y]; ++x, ++y) {
        if ((url[x] = url[y]) == '%') {
            url[x] = x2c(&url[y+1]);
            y += 2;
        }
    }
    url[x] = '\0';
}

/**
 * make_varname_unique:
 * @vname: tentative name for variable.
 * @v: the ID number for the new variable.
 * @dset: dataset information.
 *
 * Given a tentative name for a new variable, check that it
 * is not a duplicate of an existing varname.  If it is,
 * modify the new name so that it becomes unique. The ID
 * number @v is required so that, if the variable has already
 * been added to the dataset, its name does not appear to
 * conflict with itself!  If the name to be tested is not
 * associated with an existing variable, pass 0 for @v.
 *
 * Returns: the (possibly modified) variable name.
 */

char *make_varname_unique (char *vname, int v, DATASET *dset)
{
    const char *sfx = "abcdefghijklmnopqrstuvwxzy"
	"ABCDEFGHIJKLMNOPQRSTUVWXZY";
    size_t n, nmax = VNAMELEN - 8;
    char tmp[5] = {0};
    int i, j, k, vi;
    int unique = 1;

    /* first off: see if the series name is already unique! */
    vi = current_series_index(dset, vname);
    if ((v == 0 && vi < 0) || vi == v) {
	return vname;
    }

    /* strategy: cut @vname down to a length that permits
       addition of a suffix (if necessary), then add a
       suffix composed of underscore and three (ASCII)
       letters. This allows for 52^3 = 140608 unique
       suffixes.
    */

    n = strlen(vname);
    if (n > nmax) {
	n = nmax;
    }

    tmp[0] = '_';
    unique = 0;

    for (i=0; i<52 && !unique; i++) {
	tmp[1] = sfx[i];
	for (j=0; j<52 && !unique; j++) {
	    tmp[2] = sfx[j];
	    for (k=0; k<52 && !unique; k++) {
		tmp[3] = sfx[k];
		vname[n] = '\0';
		strcat(vname, tmp);
		unique = 1;
		for (vi = 1; vi < dset->v; vi++) {
		    if (vi != v && !strcmp(vname, dset->varname[vi])) {
			unique = 0;
			break;
		    }
		}
	    }
	}
    }

    if (!unique) {
	fprintf(stderr, "make_varname_unique: unresolved conflict!\n");
    }

    return vname;
}

int fix_varname_duplicates (DATASET *dset)
{
    int msg_done, dups = 0;
    int i, j;

    for (i=1; i<dset->v; i++) {
	msg_done = 0;
	for (j=i+1; j<dset->v; j++) {
	    if (!strcmp(dset->varname[i], dset->varname[j])) {
		if (!msg_done) {
		    fprintf(stderr, "'%s' duplicated variable name\n",
			    dset->varname[i]);
		    msg_done = 1;
		}
		dups = 1;
		make_varname_unique(dset->varname[j], j, dset);
	    }
	}
    }

    return dups;
}

char *append_dir (char *fname, const char *dir)
{
    size_t len;

    if (dir == NULL) {
	return fname;
    }

    len = strlen(fname);

    if (fname[len - 1] == '/' || fname[len - 1] == '\\') {
        strcat(fname, dir);
    } else {
        strcat(fname, SLASHSTR);
        strcat(fname, dir);
    }

    strcat(fname, SLASHSTR);

    return fname;
}

/**
 * path_last_element:
 * @path: path to work on.
 *
 * Returns: a pointer to the last element of @path, that is,
 * the element following the last path separator character, if any.
 * If @path does not contain a separator, @path itself is returned.
 * Note that the return value may be the empty string, if @path
 * ends with a separator.
 */

const char *path_last_element (const char *path)
{
    const char *p = strrslash(path);

    if (p == NULL) {
	p = path;
    } else {
	p++;
    }

    return p;
}

/**
 * trim_slash:
 * @s: string to work on.
 *
 * If @s ends with #SLASH, remove this character.
 *
 * Returns: the (possibly) modified string.
 */

char *trim_slash (char *s)
{
    int n = strlen(s);

    if (n > 0 && (s[n-1] == SLASH)) {
	s[n-1] = '\0';
    }

    return s;
}

/**
 * gretl_string_ends_with:
 * @s: string to examine.
 * @test: string to test for.
 *
 * Returns: 1 if @s ends with @test, else 0.
 */

int gretl_string_ends_with (const char *s, const char *test)
{
    int nt = strlen(test);
    int n = strlen(s);
    int ret = 0;

    if (n >= nt) {
	const char *p = s + n - nt;

	ret = !strcmp(p, test);
    }

    return ret;
}

/**
 * get_column_widths:
 * @strs: array of @n strings.
 * @widths: array of @n default column widths.
 * @n: number of columns.
 *
 * If need be, increases the column widths in @widths to
 * accomodate the current translations of @strs.
 */

void get_column_widths (const char **strs, int *widths, int n)
{
    int i, len;

    for (i=0; i<n; i++) {
	len = g_utf8_strlen(_(strs[i]), -1);
	if (len > widths[i]) {
	    widths[i] = len;
	}
    }
}

/**
 * gretl_utf8_strncat:
 * @dest: destination string.
 * @src: source string.
 * @n: maximum number of bytes to append.
 *
 * Works just like strncat(), except that it ensures that we
 * don't end up with an incomplete UTF-8 character preceding
 * the terminating NUL byte.
 *
 * Returns: the destination string.
 */

char *gretl_utf8_strncat (char *dest, const char *src, size_t n)
{
    const char *p = src;
    size_t b, b0 = 0;

    while (p && *p) {
	p = g_utf8_next_char(p);
	if (p) {
	    b = p - src;
	    if (b > n) {
		break;
	    }
	    b0 = b;
	}
    }

    return strncat(dest, src, b0);
}

/**
 * gretl_utf8_strncat_trim:
 * @dest: destination string.
 * @src: source string.
 * @n: maximum number of bytes to append.
 *
 * The same as gretl_utf8_strncat(), except that any leading and/or
 * trailing white space is trimmed from @dest.
 *
 * Returns: the destination string.
 */

char *gretl_utf8_strncat_trim (char *dest, const char *src, size_t n)
{
    const char *p;
    size_t b, b0 = 0;
    int i;

    src += strspn(src, " \t\r\n");
    p = src;

    while (p && *p) {
	p = g_utf8_next_char(p);
	if (p) {
	    b = p - src;
	    if (b > n) {
		break;
	    }
	    b0 = b;
	}
    }

    strncat(dest, src, b0);

    n = strlen(dest);

    for (i=n-1; i>=0; i--) {
	if (isspace(dest[i]) || dest[i] == '\r') {
	    dest[i] = '\0';
	} else {
	    break;
	}
    }

    return dest;
}

/**
 * gretl_utf8_truncate:
 * @s: string to process.
 * @nmax: maximum number of characters to retain.
 *
 * Truncates @s to a maximum length of @nmax UTF-8 characters,
 * ensuring that we don't end up with an incomplete UTF-8
 * character preceding the terminating NUL byte.
 *
 * Returns: the (possibly truncated) string.
 */

char *gretl_utf8_truncate (char *s, size_t nmax)
{
    char *p = s;
    size_t n = 0;

    while (p && *p) {
	p = g_utf8_next_char(p);
	if (p && *p) {
	    if (++n == nmax) {
		*p = '\0';
		break;
	    }
	}
    }

    return s;
}

/**
 * gretl_utf8_truncate_b:
 * @s: string to process.
 * @bmax: maximum number of bytes to retain.
 *
 * Truncates @s to a maximum length of @bmax bytes,
 * ensuring that we don't end up with an incomplete UTF-8
 * character preceding the terminating NUL byte.
 *
 * Returns: the (possibly truncated) string.
 */

char *gretl_utf8_truncate_b (char *s, size_t bmax)
{
    char *p = s;
    size_t b = 0;

    while (p && *p) {
	p = g_utf8_next_char(p);
	b = p - s;
	if (b == bmax) {
	    *p = '\0';
	    break;
	} else if (b > bmax) {
	    p = g_utf8_prev_char(p);
	    *p = '\0';
	    break;
	}
    }

    return s;
}

/**
 * gretl_utf8_select:
 * @s: source string.
 * @list: list of (1-based) indices of characters to select.
 *
 * Constructs a string composed of the selected characters from @s,
 * which do not have to be contiguous or in order.
 *
 * Returns: the constructed string.
 */

char *gretl_utf8_select (const char *s, const int *list)
{
    char *ret = NULL;
    int i, k;

    if (gretl_is_ascii(s)) {
	int j = 0;

	ret = calloc(list[0]+1, 1);
	for (i=1; i<=list[0]; i++) {
	    k = list[i] - 1; /* 0-based */
	    ret[j++] = s[k];
	}
    } else {
	GString *gs = g_string_new(NULL);
	int len = g_utf8_strlen(s, -1);
	gchar const **bits;
	gchar *si;

	bits = calloc(len + 1, sizeof *bits);
	bits[0] = s;
	for (i=1; i<len; i++) {
	    s = bits[i] = g_utf8_find_next_char(s, NULL);
	}
	for (i=1; i<=list[0]; i++) {
	    k = list[i] - 1; /* 0-based */
	    si = g_utf8_substring(bits[k], 0, 1);
	    g_string_append(gs, si);
	    g_free(si);
	}

	ret = gretl_strdup(gs->str);
	g_string_free(gs, TRUE);
	free(bits);
    }

    return ret;
}

/**
 * gretl_utf8_replace_char:
 * @targ: the target for replacement.
 * @src: the UTF-8 character to insert (NUL terminated).
 * @pos: the position, in UTF-8 characters, at which to replace.
 *
 * Notes: @pos must be prechecked as within bounds, and
 * @src must be prechecked for containing a single character.
 *
 * Returns: newly allocated resulting string.
 */

char *gretl_utf8_replace_char (const char *targ, const char *src,
			       int pos)
{
    const char *s = g_utf8_offset_to_pointer(targ, pos);
    const char *p = g_utf8_find_next_char(s, NULL);
    char *ret;
    int len;

    len = strlen(targ) + strlen(src) + 1 - (p - s);
    ret = calloc(len, 1);
    strncat(ret, targ, s - targ);
    strcat(ret, src);
    strcat(ret, p);

    return ret;
}

/**
 * gretl_scan_varname:
 * @src: source string.
 * @targ: target string.
 *
 * Performs sscanf() on @src, using a conversion specifier
 * which allows for writing up to VNAMELEN-1 bytes into
 * @targ (stopping at white space); @targ therefore be at
 * least VNAMELEN bytes long. No checking is done for the
 * validity of the scanned string as a gretl identifier.
 *
 * Returns: the return value from sscanf().
 */

int gretl_scan_varname (const char *src, char *targ)
{
    char fmt[8];

    sprintf(fmt, "%%%ds", VNAMELEN-1);
    return sscanf(src, fmt, targ);
}

/**
 * gretl_normalize_varname:
 * @targ: target string.
 * @src: source string.
 * @underscore: flag to replace all illegal characters
 * with underscore.
 * @seq: sequence number in array of names, if applicable.
 *
 * Writes a valid gretl identifier to @targ, which must be
 * at least #VNAMELEN bytes in length, taking @src as basis
 * and replacing any illegal characters as described in the
 * documentation for the userland fixname function.
 *
 * Returns: 1 if any changes were required, 0 if not.
 */

int gretl_normalize_varname (char *targ, const char *src,
			     int underscore, int seq)
{
    const char *letters = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int n, ret = 0, i = 0;

    /* skip any leading non-letters */
    n = strcspn(src, letters);
    if (n > 0) {
	src += n;
	ret = 1;
    }

    while (*src && i < VNAMELEN - 1) {
	if (strspn(src, letters) > 0 || isdigit(*src) || *src == '_') {
	    /* transcribe valid characters */
	    targ[i++] = *src;
	} else {
	    if (*src == ' ' || underscore) {
		/* convert space to underscore */
		if (i > 0 && targ[i-1] == '_') {
		    ; /* skip */
		} else {
		    targ[i++] = '_';
		}
	    }
	    ret = 1;
	}
	src++;
    }

    if (i > 0) {
	targ[i] = '\0';
    } else {
	if (seq <= 0) {
	    strcpy(targ, "col[n]");
	} else {
	    sprintf(targ, "col%d", seq);
	}
	ret = 1;
    }

    return ret;
}

/**
 * gretl_regexp_replace:
 * @orig: the original string.
 * @match: the pattern to match.
 * @repl: the replacement expression for @match.
 * @err: location to receive error code.
 *
 * Builds a string based on @orig but in which all
 * occurrences of @match (which is interpreted as a
 * regular expression of the Perl type) are replaced
 * by means of @repl (also interpreted as a regular
 * expression).
 *
 * Returns: newly allocated string or NULL on failure.
 */

char *gretl_regexp_replace (const char *orig,
			    const char *match,
			    const char *repl,
			    int *err)
{
    GRegex *regex;
    GError *error = NULL;
    char *mod = NULL;

    regex = g_regex_new(match, 0, 0, &error);

    if (error == NULL) {
	mod = g_regex_replace(regex, orig, -1, 0, repl, 0, &error);
    }

    if (error != NULL) {
	*err = 1;
	gretl_errmsg_set(error->message);
	g_error_free(error);
    }

    if (regex != NULL) {
	g_regex_unref(regex);
    }

    return mod;
}

/**
 * gretl_literal_replace:
 * @orig: the original string.
 * @match: the substring to match.
 * @repl: the replacement string for @match.
 * @err: location to receive error code.
 *
 * Builds a string based on @orig but in which all
 * occurrences of @match (which is interpreted as a
 * straight string literal) are replaced by @repl (also
 * a straight string literal).
 *
 * Returns: newly allocated string or NULL on failure.
 */

char *gretl_literal_replace (const char *orig,
			     const char *match,
			     const char *repl,
			     int *err)
{
    char *mod = NULL;
    const char *q, *r;
    int mlen = strlen(match);
    int nrep = 0;

    if (mlen > 0) {
	/* count the occurrences of @match */
	q = orig;
	while ((r = strstr(q, match)) != NULL) {
	    nrep++;
	    q = r + mlen;
	}
    }

    if (nrep == 0) {
	/* no replacement needed */
	mod = gretl_strdup(orig);
    } else {
	int ldiff = nrep * (strlen(repl) - mlen);

	mod = malloc(strlen(orig) + ldiff + 1);
	if (mod != NULL) {
	    q = orig;
	    *mod = '\0';
	    while ((r = strstr(q, match)) != NULL) {
		strncat(mod, q, r - q);
		strcat(mod, repl);
		q = r + mlen;
	    }
	    if (*q) {
		strcat(mod, q);
	    }
	}
    }

    if (mod == NULL) {
	*err = E_ALLOC;
    }

    return mod;
}

/**
 * gretl_substring:
 * @str: the string to operate on.
 * @first: 1-based index of initial character.
 * @last: 1-based index of final character, or -1 to go to the end.
 * @err: location to receive error code.
 *
 * Returns: a substring of @str, from @first to @last.
 */

char *gretl_substring (const char *str, int first, int last, int *err)
{
    int len, ini, fin, sublen;
    char *ret;

    len = g_utf8_strlen(str, -1);
    if (last == -1) {
	last = len;
    }

    if (first <= 0 || last <= 0) {
	gretl_errmsg_sprintf("Index value %d is out of bounds",
			     first <= 0 ? first : last);
	*err = E_DATA;
    }

    ini = (first < 1) ? 1 : ((first > len) ? len : first);
    fin = (last < 1) ? 1 : ((last > len) ? len : last);
    sublen = (fin >= ini) ? fin - ini + 1 : 0;

    if (sublen == 0) {
	ret = calloc(1, 1);
    } else {
	const char *s1;
	int i;

	for (i=1; i<ini; i++) {
	    str = g_utf8_next_char(str);
	}
	s1 = str;
	for (i=ini; i<=last; i++) {
	    str = g_utf8_next_char(str);
	}
	len = str - s1;
	ret = calloc(len + 1, 1);
	if (ret != NULL) {
	    *ret = '\0';
	    gretl_utf8_strncat(ret, s1, len);
	}
    }

    if (ret == NULL) {
	*err = E_ALLOC;
    }

    return ret;
}

GretlCase gretl_case_from_string (const char *s)
{
    if (!strcmp(s, "lower")) {
        return CASE_LOWER;
    } else if (!strcmp(s, "upper")) {
        return CASE_UPPER;
    } else if (!strcmp(s, "camel")) {
        return CASE_CAMEL;
    } else if (!strcmp(s, "snake")) {
        return CASE_SNAKE;
    } else {
        return CASE_NONE;
    }
}

char *gretl_change_case (const char *s, GretlCase c, int *err)
{
    char *ret = NULL;
    int slen, clen;
    int i, j, delta = 0;

    if (s == NULL || c < CASE_LOWER || c > CASE_SNAKE) {
        *err = E_INVARG;
        return NULL;
    }

    clen = slen = strlen(s);

    for (i=0; i<slen; i++) {
        if (c == CASE_LOWER && isupper(s[i])) {
            delta++;
        } else if (c == CASE_UPPER && islower(s[i])) {
            delta++;
        } else if (c == CASE_CAMEL && s[i] == '_') {
            delta++;
            clen--;
        } else if (c == CASE_SNAKE && i > 0 && isupper(s[i])) {
            delta++;
            clen++;
        }
    }

    if (delta > 0) {
        ret = calloc(clen + 1, 1);
        if (ret == NULL) {
            *err = E_ALLOC;
            return NULL;
        }
        for (i=0, j=0; i<slen; i++) {
            if (c == CASE_LOWER) {
                ret[i] = isupper(s[i]) ? tolower(s[i]) : s[i];
            } else if (c == CASE_UPPER) {
                ret[i] = islower(s[i]) ? toupper(s[i]) : s[i];
            } else if (c == CASE_CAMEL) {
                if (s[i] == '_') {
                    if (islower(s[i+1])) {
                        ret[j++] = toupper(s[++i]);
                    }
                } else {
                    ret[j++] = tolower(s[i]);
                }
            } else if (c == CASE_SNAKE) {
                if (i > 0 && isupper(s[i])) {
                    ret[j++] = '_';
                    ret[j++] = tolower(s[i]);
                } else {
                    ret[j++] = s[i];
                }
            }
        }
    }

    return ret;
}
