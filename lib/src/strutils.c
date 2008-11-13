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
	    if (!isspace((unsigned char) *s) && *s != CTRLZ) {
		ret = 0;
		break;
	    }
	    s++;
	}
    }

    return ret;
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
    double x;

#ifdef ENABLE_NLS
    static int decpoint = 0;

    if (decpoint == 0) {
	struct lconv *lc;

	lc = localeconv();
	decpoint = *lc->decimal_point;
    }
    if (decpoint == '.') {
	x = atof(s);
    } else {
	gretl_push_c_numeric_locale();
	x = atof(s);
	gretl_pop_c_numeric_locale();
    }
#else
    x = atof(s);
#endif

    return x;
}

/**
 * dotpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last "." within @str,
 * or strlen(@str) in case a dot is not found, or the string
 * ends with a (backward or forward) slash.
 */

int dotpos (const char *str)
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

/**
 * slashpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last #SLASH within @str,
 * or 0 in case a #SLASH is not found.
 */

int slashpos (const char *str)
{ 
    int i, p = 0;

    if (str != NULL && *str != '\0') {
	p = strlen(str);
	for (i=p-1; i>0; i--) {
	    if (str[i] == SLASH) {
		p = i;
		break;
	    }
	}
    }

    return p;    
}

/**
 * delchar:
 * @c: the character to delete.
 * @str: the string from which to delete @c.
 *
 * Deletes all instances of @c within @str.
 *
 * Returns: the possibly modified string.
 */

char *delchar (int c, char *str)
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
 * haschar:
 * @c: the character to look for.
 * @s: the string to examine.
 *
 * Returns: the first position of @c in @s, or -1 if @c is not
 * found.
 */

int haschar (char c, const char *s)
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
 * lastchar:
 * @c: the character to look for.
 * @s: the string to examine.
 *
 * Returns: 1 if @c is the last character in @s, 0 otherwise
 */

int lastchar (char c, const char *s)
{
    int ret = 0;

    if (s != NULL && s[strlen(s) - 1] == c) {
	ret = 1;
    }

    return ret;
}

/**
 * charsub:
 * @str: the string to operate on.
 * @find: the character to replace.
 * @repl: the replacement character.
 *
 * Replaces all occurrences of @find with @repl in @str.
 *
 * Returns: the (possibly modified) string.
 */

char *charsub (char *str, char find, char repl)
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
 * string to ensure that all the numbers are comma-
 * separated.
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
 * @sfx: the suffix to check for.
 *
 * Returns: 1 if @str ends with @sfx (on a case-insensitive
 * comparison), 0 otherwise.
 */

int has_suffix (const char *str, const char *sfx)
{
    int diff, ret = 0;

    if (str != NULL && sfx != NULL) {
	diff = strlen(str) - strlen(sfx);
	if (diff >= 0) {
	    ret = 1;
	    str += diff;
	    while (*str) {
		if (*str != *sfx && *str != toupper(*sfx)) {
		    ret = 0;
		    break;
		}
		str++;
		sfx++;
	    }
	}
    }
    
    return ret;
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

    if (!strcmp(str, "inf") || !strcmp(str, "nan")) {
	/* could be variable names: they are not valid numbers */
	return 0;
    }

    gretl_push_c_numeric_locale();

    errno = 0;

    strtod(str, &test);
    if (*test != '\0' || errno == ERANGE) {
	ret = 0;
    }

    gretl_pop_c_numeric_locale();

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
 * lower:
 * @str: the string to transform.
 *
 * Converts any upper case characters in @str to lower case.
 *
 * Returns: the possibly modified string.
 */

char *lower (char *str)
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
 * Returns: an allocated copy of @src, or %NULL on error.
 */

char *gretl_strdup (const char *src)
{
    char *targ = NULL;

    if (src != NULL) {
	targ = malloc(strlen(src) + 1);
	if (targ != NULL) {
	    strcpy(targ, src);
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
 * @src, or %NULL on error.
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
	    *targ = '\0';
	    strncat(targ, src, len);
	}
    }

    return targ;
}

/**
 * gretl_strdup_printf:
 * @template: as in printf().
 * @Varargs: arguments to be printed.
 *
 * Print the arguments according to @format.
 * 
 * Returns: allocated result of the printing, or %NULL on failure.
 */

char *gretl_strdup_printf (const char *template, ...)
{
    va_list args;
    int plen, bsize = 2048;
    char *buf;

    buf = malloc(bsize);
    if (buf == NULL) {
	return NULL;
    }

    memset(buf, 0, 1);

    va_start(args, template);
    plen = vsnprintf(buf, bsize, template, args);
    va_end(args);

    if (plen >= bsize) {
	fputs("gretl_strdup_printf warning: string was truncated\n",
	      stderr);
    }

    return buf;
}

/**
 * gretl_str_expand:
 * @orig: pointer to the base string.
 * @add: the string to be added.
 * @sep: string to be interpolated, or %NULL.
 *
 * Creates a newly allocated string built by concatenating
 * @orig and @add, with @sep interpolated unless @sep is
 * %NULL, and replaces the content of @orig with the new string.
 * As a special case, if @orig is %NULL, or if the content of
 * @orig is %NULL, we just duplicate @add.
 *
 * Returns: the reallocated string, or %NULL on failure.  In case
 * of failure the content of @orig is freed, if @orig is not %NULL,
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

#define is_word_char(c) (isalnum((unsigned char) c) || c == '_')

/**
 * gretl_word_strdup:
 * @src: the source string.
 * @ptr: location to receive end of word pointer, or %NULL.
 *
 * Copies the first 'word' found in @src, where a word
 * is defined as consisting of alphanumeric characters
 * and the underscore.  If @ptr is not %NULL, on exit it
 * points at the next position in @src after the copied
 * word.
 *
 * Returns: the allocated word or %NULL in case no word is
 * found, or if allocation fails.
 */

char *gretl_word_strdup (const char *src, const char **ptr)
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

	while (*src && !is_word_char(*src)) {
	    src++;
	}

	p = src;

	while (is_word_char(*src)) {
	    len++;
	    src++;
	}

	if (ptr != NULL) {
	    *ptr = src;
	}

	if (len > 0) {
	    targ = gretl_strndup(p, len);
	}
    }

    return targ;
}

/**
 * gretl_quoted_string_strdup:
 * @s: the source string.
 * @ptr: location to receive end pointer, or %NULL.
 *
 * If @s starts with a quote (double or single), return a copy of  
 * the portion of @s that is enclosed in quotes.  That is, 
 * from @s + 1 up to but not including the next matching quote.
 * If @ptr is not %NULL, on output it receives a pointer to
 * the next byte in @s after the closing quote.
 *
 * Returns: the allocated string or %NULL on failure.
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

/**
 * gretl_string_split:
 * @s: the source string.
 * @n: location to receive the number of substrings.
 *
 * Parses @s into a set of zero or more substrings, separated
 * by one or more spaces, and creates an array of those substrings. 
 * On sucessful exit, @n holds the number of substrings. 
 *
 * Returns: the allocated array or %NULL in case of failure.
 */

char **gretl_string_split (const char *s, int *n)
{
    int i, k, m = count_fields(s);
    char *word;
    char **S;

    *n = 0;

    if (m == 0) {
	return NULL;
    }

    S = strings_array_new(m);
    if (S == NULL) {
	return NULL;
    }

    for (i=0; i<m; i++) {
	s += strspn(s, " ");
	k = strcspn(s, " ");
	word = gretl_strndup(s, k);
	if (word == NULL) {
	    free_strings_array(S, m);
	    return NULL;
	}
	S[i] = word;
	s += k;
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
	str[n] = 0;
    }

    return str;
}

/**
 * gretl_namechar_spn:
 * @s: the string to examine.
 *
 * Returns: the length of the intial segment of @s which
 * consists of characters that are valid in a gretl
 * variable or object name, namely a-z, A-Z, 0-9 and _,
 * starting with a letter.
 */

int gretl_namechar_spn (const char *s)
{
    const char *ok = "abcdefghijklmnopqrstuvwxyz"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
	"0123456789_";
    int ret = 0;

    if (isalpha(*s)) {
	ret = strspn(s, ok);
    }

    return ret;
}

/**
 * clear:
 * @str: the string to clear.
 * @len: the length of the string to be cleared.
 *
 * Sets all bytes in @str to 0.
 */

void clear (char *str, int len)
{
    memset(str, 0, len);
}

/**
 * count_fields:
 * @s: the string to process.
 *
 * Returns: the number of space-separated fields in @s.
 */

int count_fields (const char *s)
{
    int nf = 0;

    if (s != NULL && *s != '\0') {
	const char *p;

	/* step past any leading space */
	while (*s == ' ') {
	    s++;
	}

	if (*s != '\0' && *s != ' ') {
	    s++;
	    nf++;
	}

	while (*s) {
	    p = strpbrk(s, " ");
	    if (p != NULL) {
		s = p + strspn(p, " ");
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
 * chopstr:
 * @str: the string to process.
 *
 * Removes both leading and trailing space from a string.
 *
 * Returns: the possibly modified string.
 */

char *chopstr (char *str)
{
    int i, n = strspn(str, " \t");

    if (n > 0) {
	shift_string_left(str, n);
    }

    n = strlen(str);

    for (i = n - 1; i >= 0; i--) {
	if (isspace(str[i]) || str[i] == '\r') {
	    str[i] = '\0';
	} else {
	    break;
	}
    }

    return str;
}

/**
 * switch_ext:
 * @targ: the target or output string (must be pre-allocated).
 * @src: the source or input string.
 * @ext: the extension or suffix to attach.
 *
 * For processing filenames: copies @src to @targ, minus any existing
 * filename extension, and adds to @targ the specified extension.
 *
 * Returns: the output string, @targ.
 */

char *switch_ext (char *targ, const char *src, char *ext)
{
    int i = dotpos(src);

    if (targ != src) {
        strncpy(targ, src, i);
    }

    targ[i] = '.';
    targ[i + 1] = 0;
    strcat(targ, ext);

    return targ;
}

/**
 * get_base:
 * @targ: the target or output string (must be pre-allocated).
 * @src: the source or input string.
 * @c: the "base marker" character.
 *
 * If @c is found in @src, puts into @targ the portion of @src up to and 
 * including the last occurrence of @c within @src.
 *
 * Returns: 1 if @c is found in @str, otherwise 0.
 */

int get_base (char *targ, const char *src, char c)
{
    int ret = 0;

    if (src != NULL && *src != '\0') {
	int i, n = strlen(src);

	for (i=n-1; i>=0; i--) {
	    if (src[i] == c) {
		*targ = '\0';
		strncat(targ, src, i + 1);
		ret = 1;
		break;
	    }
	}
    }

    return ret;
}

/**
 * top_n_tail:
 * @str: the string to process.
 * @maxlen: maximum length of string, including NUL termination.
 * @err: location to receive error code, or %NULL.
 *
 * Drop leading space and trailing space and newline from string,
 * then replace a trailing backslash (if any) with a space.
 * If @str does not end with a newline within the limit set by
 * @maxlen, and @err is not %NULL, then %E_TOOLONG is written 
 * to @err.
 * 
 * Returns: 1 if a trailing backslash or comma was found, 
 * otherwise 0.
 */

int top_n_tail (char *str, size_t maxlen, int *err)
{
    int i, n, cont = 0;

    if (str == NULL || *str == 0 || *str == '\n' || *str == '\r') {
	return 0;
    }

    n = strlen(str) - 1;

    if (err != NULL && n > maxlen - 2 && str[n] != '\n') {
	*err = E_TOOLONG;
    }

    /* chop any trailing space */
    for (i=n; i>=0; i--) {
	if (isspace((unsigned char) str[i])) {
	    str[i] = 0;
	} else {
	    break;
	}
    }

    if (*str != 0) {
	/* Drop any leading spaces, also possible questionmark.  Try
	   to catch non-breaking spaces too -- ugh, Windows!
	*/
	i = 0;
	while (isspace((unsigned char) str[i]) || 
	       str[i] == '?' ||
	       str[i] == (char) 0xC2 ||
	       str[i] == (char) 0xA0) {
	    i++;
	}
	if (i > 0) {
	    shift_string_left(str, i);
	}

	/* does this line start a comment? */
	if (strchr(str, '#') || !strncmp(str, "/*", 2)) {
	    ; /* leave well alone */
	} else {
	    /* replace backslash, if present */
	    n = strlen(str) - 1;
	    if (n >= 0) {
		if (str[n] == '\\') {
		    str[n] = ' ';
		    cont = 1;
		} else if (str[n] == ',') {
		    cont = 1;
		}
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

    if (str == NULL || *str == 0) {
	return str;
    }

    len = strlen(str);

    for (i=len-1; i>=0; i--) {
	if (isspace((unsigned char) str[i]) ||
	    str[i] == '\n' || str[i] == '\r') {
	    str[i] = 0;
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

    if (s == NULL || *s == 0) {
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
		if (*p == 0) break;
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
 * safecpy:
 * @targ: target or output string (must be pre-allocated).
 * @src: source or input string.
 * @n: maximum length of target string.
 *
 * Copies at most @n characters from @src to @targ, and ensures that
 * @targ[@n] is a NUL byte.
 * 
 * Returns: the output string.
 */

char *safecpy (char *targ, const char *src, int n)
{
    *targ = 0;
    strncat(targ, src, n);
    return targ;
}

/**
 * strings_array_new:
 * @nstrs: number of strings in array.
 *
 * Allocates storage for @nstrs strings and initalizes all 
 * to %NULL.
 * 
 * Returns: the allocated array, or %NULL on failure.
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
 * strings_array_new_with_length:
 * @nstrs: number of strings in array.
 * @len: number of bytes per string.
 *
 * Allocates storage for @nstrs strings, each of them 
 * @len bytes long.  The first byte of each string is
 * initialized to 0.
 * 
 * Returns: the allocated array, or %NULL on failure.
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
 * strings_array_dup:
 * @strs: array of strings to be copied.
 * @n: number of strings in array.
 *
 * Returns: an allocated copy of @strs, or %NULL on failure.
 */

char **strings_array_dup (char **strs, int n)
{
    char **S;
    int i, j;

    if (n <= 0 || strs == NULL) {
	return NULL;
    }

    S = malloc(n * sizeof *S);
    if (S == NULL) return NULL;

    for (i=0; i<n; i++) {
	if (strs[i] == NULL) {
	    S[i] = NULL;
	} else {
	    S[i] = gretl_strdup(strs[i]);
	    if (S[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(S[j]);
		}
		free(S);
		return NULL;
	    }
	}
    }

    return S;
}

/**
 * free_strings_array:
 * @strs: array of allocated strings.
 * @nstrs: number of strings in array.
 *
 * Frees each allocated string in @strs, then frees @strs itself.
 * Checks that @strs is not %NULL before proceeding.
 */

void free_strings_array (char **strs, int nstrs)
{
    int i;

    if (strs != NULL) {
	for (i=0; i<nstrs; i++) {
	    free(strs[i]);
	}
	free(strs);
    }
}

static int maybe_year_plus (const char *s)
{
    int n = strlen(s);

    if (n > 8 && isdigit(s[0]) && isdigit(s[1]) &&
	isdigit(s[2]) && isdigit(s[3]) && !isdigit(s[4])) {
	return 1;
    }

    return 0;
}

static char *
real_get_obs_string (char *s, int t, const DATAINFO *pdinfo, int full)
{
    if (pdinfo->markers && pdinfo->S != NULL) { 
	/* data marker strings present */
	if (full) {
	    strcpy(s, pdinfo->S[t]);
	} else {
	    *s = '\0';
	    if (maybe_year_plus(pdinfo->S[t])) {
		/* daily date, year first? */
		strncat(s, pdinfo->S[t] + 2, 8);
	    } else {
		strncat(s, pdinfo->S[t], 8);
	    }
	}
    } else {
	if (full) {
	    ntodate_full(s, t, pdinfo);
	} else {
	    ntodate(s, t, pdinfo);
	}
    }

    return s;
}

/**
 * get_obs_string:
 * @obs: char array big enough to hold the observation (#OBSLEN).
 * @t: zero-based observation number.
 * @pdinfo: pointer to dataset information.
 *
 * Returns: the observation string corresponding to @t.
 */

char *get_obs_string (char *obs, int t, const DATAINFO *pdinfo)
{
    return real_get_obs_string(obs, t, pdinfo, 0);
}

/**
 * get_full_obs_string:
 * @obs: char array big enough to hold the observation (#OBSLEN).
 * @t: zero-based observation number.
 * @pdinfo: pointer to dataset information.
 *
 * Returns: the observation string corresponding to @t, using
 * a four-digit year in case of dated daily data.
 */

char *get_full_obs_string (char *obs, int t, const DATAINFO *pdinfo)
{
    return real_get_obs_string(obs, t, pdinfo, 1);
}

/**
 * obs_str_to_double:
 * @obs: string representation of observation number.
 *
 * Returns: the floating-point counterpart of @obs.
 */

double obs_str_to_double (const char *obs)
{
    char tmp[OBSLEN];
    char *p;

    strcpy(tmp, obs);
    p = tmp;

    while (*p) {
	if (*p == ':' || *p == ',') {
	    *p = '.';
	}
	p++;
    }

    return dot_atof(tmp);
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
 * ntodate();
 */

void modify_date_for_csv (char *s, int pd)
{
    if (pd == 4) {
	charsub(s, ':', 'Q');
    } else {
	charsub(s, ':', 'M');
    }
}

/**
 * csv_obs_string_to_prn:
 * @t: 0-based observation number.
 * @pdinfo: data set information struct.
 * @prn: printing struct.
 *
 * Prints the observation string corresponding to obervation @t to
 * @prn, in a format suitable for a CSV file.
 */

void csv_obs_to_prn (int t, const DATAINFO *pdinfo, PRN *prn)
{
    if (pdinfo->S != NULL) {
	pprintf(prn, "%s%c", pdinfo->S[t], pdinfo->delim);
    } else if (pdinfo->structure != CROSS_SECTION) {
	char tmp[OBSLEN];

	ntodate_full(tmp, t, pdinfo);
	if (quarterly_or_monthly(pdinfo)) {
	    modify_date_for_csv(tmp, pdinfo->pd);
	}
	pprintf(prn, "%s%c", tmp, pdinfo->delim);
    }
}
	
/**
 * print_time:
 * @timep: time to print.
 *
 * Returns: pointer to a static string containing a locale-dependent
 * representation of @timep.  In English, this will be in the format
 * Y/m/d H:M.
 */

const char *print_time (const time_t *timep)
{
    static char timestr[48];
    struct tm *local;

    local = localtime(timep);

    strftime(timestr, 47, "%Y/%m/%d %H:%M", local);

    return timestr;
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
	sprintf(gretl_errmsg, _("out of memory in XML encoding"));
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
 * @pdinfo: dataset information.
 *
 * Given a tenative name for a new variable, check that it
 * is not a duplicate of an existing varname.  If it is,
 * modify the new name so that it becomes unique. The ID
 * number @v is required so that, if the variable has already
 * been added to the dataset, its name does not appear to 
 * conflict with itself!  If the name to be tested is not
 * associated with an existing variable, pass 0 for @v.
 *
 * Returns: the (possibly modified) variable name.
 */

char *make_varname_unique (char *vname, int v, DATAINFO *pdinfo) 
{
    const char *add = "abcdefghijklmnopqrstuvwxyz";
    size_t n = strlen(vname);
    int i, j, conflict;

    if (n > 7) {
	n = 7;
    }

    for (j=0; j<26; j++) {
	conflict = 0;
	for (i=1; i<pdinfo->v && !conflict; i++) {
	    if (i != v && !strcmp(vname, pdinfo->varname[i])) {
		conflict = 1;
	    }
	}
	if (conflict) {
	    vname[n] = add[j];
	    vname[n+1] = '\0';
	} else {
	    break;
	}
    }

    return vname;
}

int fix_varname_duplicates (DATAINFO *pdinfo)
{
    int dups = 0;
    int i, j;

    for (i=1; i<pdinfo->v; i++) {
	for (j=i+1; j<pdinfo->v; j++) {
	    if (!strcmp(pdinfo->varname[i], pdinfo->varname[j])) {
		fprintf(stderr, "'%s' duplicated variable name\n", pdinfo->varname[i]);
		dups = 1;
		make_varname_unique(pdinfo->varname[j], j, pdinfo);
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
 * build_path:
 * @targ: target string to write to (must be pre-allocated).
 * @dirname: first part of path.
 * @fname: filename.
 * @ext: filename extension to be appended (or %NULL).
 *
 * Writes to @targ a full path composed of @dirname,
 * @fname and (optionally) @ext.  This function ensures
 * that an appropriate separator is inserted between 
 * @dirname and @fname, if @dirname is not already
 * terminated with such a separator.
 *
 * Returns: the target string, @targ.
 */

char *build_path (char *targ, const char *dirname, const char *fname, 
		  const char *ext)
{
    size_t len;

    if (dirname == NULL || fname == NULL || targ == NULL) {
	return NULL;
    }

    *targ = '\0';
    strcat(targ, dirname);
    len = strlen(targ);
    if (len == 0) {
	return NULL;
    }

    /* strip a trailing single dot */
    if (len > 1 && targ[len-1] == '.' && 
	(targ[len-2] == '/' || targ[len-2] == '\\')) {
	    targ[len-1] = '\0';
    }

    if (targ[len-1] == '/' || targ[len-1] == '\\') {
        /* dirname is already properly terminated */
        strcat(targ, fname);
    } else {
        /* otherwise put a separator in */
        strcat(targ, SLASHSTR);
        strcat(targ, fname);
    }

    if (ext != NULL) {
	strcat(targ, ext);
    }

    return targ;
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
    const char *p = strrchr(path, SLASH);

#ifdef WIN32
    if (p == NULL) {
	p = strrchr(path, '\\');
    }
#endif

    if (p == NULL) {
	p = path;
    } else {
	p++;
    }

    return p;
}

/**
 * ensure_slash:
 * @s: string to work on.
 *
 * If @s does not ends with #SLASH, append this character.
 *
 * Returns: the (possibly) modified string.
 */

char *ensure_slash (char *s)
{
    int n = strlen(s);

    if (n > 0 && s[n-1] != SLASH) {
	strcat(s, SLASHSTR);
    }

    return s;
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

int *varname_match_list (const DATAINFO *pdinfo, const char *pattern)
{
    GPatternSpec *pspec;
    int *list = NULL;
    int i, n = 0;

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<pdinfo->v; i++) { 
	if (g_pattern_match_string(pspec, pdinfo->varname[i])) {
	    n++;
	}
    }

    if (n > 0) {
	list = malloc((n + 1) * sizeof *list);
	if (list != NULL) {
	    int j = 1;

	    list[0] = n;
	    for (i=1; i<pdinfo->v; i++) { 
		if (g_pattern_match_string(pspec, pdinfo->varname[i])) {
		    list[j++] = i;
		}
	    }
	}
    }

    g_pattern_spec_free(pspec);

    return list;
}

int varname_match_any (const DATAINFO *pdinfo, const char *pattern)
{
    GPatternSpec *pspec;
    int i, ret = 0;

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<pdinfo->v; i++) { 
	if (g_pattern_match_string(pspec, pdinfo->varname[i])) {
	    ret = 1;
	    break;
	}
    }

    g_pattern_spec_free(pspec);

    return ret;
}
