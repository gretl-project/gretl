/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* strutils.c for gretl */

#include "libgretl.h"
#include <stdarg.h>

/**
 * dotpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last "." within @str,
 * or strlen(@str) in case a dot is not found, or the string
 * ends with a (backward or forward) slash.
 *
 */

int dotpos (const char *str)
{ 
    int i, n = strlen(str);

    for (i=n-1; i>0; i--) { 
	if (str[i] == '/' || str[i] == '\\') return n;
	if (str[i] == '.') return i;
    }
    return n;    
}

/**
 * slashpos:
 * @str: the string to examine.
 *
 * Returns: the integer position of the last #SLASH within @str,
 * or 0 in case a #SLASH is not found.
 *
 */

int slashpos (const char *str)
{ 
    size_t i, n;

    if (str == NULL) return 0;

    n = strlen(str);
    for (i=n-1; i>0; i--) 
	if (str[i] == SLASH) return i;
    return 0;    
}

/**
 * copy:
 * @str: the source string.
 * @indx: position in @str from which to start the copying.
 * @count: number of characters to copy.
 * @dest: destination string (must be pre-allocated).
 *
 * Copies @count characters from @indx in @str to @dest.
 *
 */

void copy (const char *str, int indx, int count, char *dest)
{
    int i;

    dest[0] = '\0';
    for (i=0; i<count; ++i) dest[i] = str[indx+i];
    dest[count] = '\0';
}

/**
 * delchar:
 * @c: the character to delete.
 * @str: the string from which to delete @c.
 *
 * Deletes all instances of @c within @str.
 *
 */

void delchar (int c, char *str)
{
    int i, j;

    for (i=j=0; str[i] != '\0'; i++)
	if (str[i] != c)
	    str[j++] = str[i];
    str[j] = '\0';

}

/**
 * _delete:
 * @str: the string to process.
 * @indx: the starting point for deleting characters.
 * @count: the number of characters to delete.
 *
 * Deletes @count characters from @str, starting at position @indx.
 *
 */

void _delete (char *str, int indx, int count)
{
    size_t i, n = strlen(str);

    for (i=indx; i<=n-count; ++i) 
	str[i] = str[count+i];
}

/**
 * haschar:
 * @c: the character to look for.
 * @str: the string to examine.
 *
 * Returns: the first position of @c in @str, or -1 if @c is not
 * found.
 *
 */

int haschar (char c, const char *str)
{
    size_t i = 0, n = strlen(str);

    do {
        if (str[i] == c) return i;
        i++;
    } while (i < n);
    return -1;
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
	if (*p == find) *p = repl;
	p++;
    }
    return str;
}

/**
 * _isnumber:
 * @str: the string to examine.
 *
 * Returns: 1 if the given @str is numeric, otherwise 0.
 *
 */

int _isnumber (const char *str)
{
    size_t i, n = strlen(str);
    int decimal = 0, efound = 0;
    char c;

    c = str[0];
    if (c != '+' && c !='-' && c != '.' && !isdigit((unsigned char) c))
        return 0;
    for (i=1; i<=n-1; i++) {
        c = str[i];
        if (c == '.') {
            if (decimal) return 0;
            else decimal = 1;
        }
        else if (c == 'e' || c == 'E')  {
            if (efound) return 0;
            i++;
            efound = 1;
        }
        else if (isdigit((unsigned char) c)) continue;
	else return 0;
    }
    return 1;
}

/**
 * lower:
 * @str: the string to transform.
 *
 * Converts any upper case characters in @str to lower case.
 *
 */

void lower (char *str)
{
    while (*str) {
        if (isupper((unsigned char) *str)) *str = tolower(*str);
        str++;
    }
}

/**
 * _esl_trunc:
 * @str: the string to truncate.
 * @n: the desired length of the truncated string.
 *
 * Truncates the given @str to the specified length.
 *
 */

void _esl_trunc (char *str, int n)
{
    size_t len = strlen(str);

    if (len > n) 
	_delete(str, n, len - n);  
}

/**
 * clear:
 * @str: the string to clear.
 * @len: the length of the string to be cleared.
 *
 * Sets all bytes in @str to 0.
 *
 */

void clear (char *str, int len)
{
    memset(str, 0, len);
}

/**
 * _count_fields:
 * @str: the string to process.
 *
 * Returns: the number of space-separated fields in @str.
 *
 */

int _count_fields (const char *str)
{
    int n = 0;
    char tmpstr[MAXLEN];

    strcpy(tmpstr, str);
    if (strtok(tmpstr, " ")) n++;
    while (strtok(NULL, " ")) n++;
    return n;
}

/**
 * _shiftleft:
 * @str: the string to process.
 * @move: the number of places to shift.
 *
 * Shifts the content of @str left by @move places, dropping
 * leading bytes as needed.
 *
 */

void _shiftleft (char *str, size_t move)
{
    size_t n = strlen(str);

    if (move >= n) {
	str[0] = '\0';
	return;
    }
    memmove(str, str + move, n - move);
    str[n - move] = '\0';
}

/**
 * chopstr:
 * @str: the string to process.
 *
 * Removes both leading and trailing space from a string.
 *
 */

void chopstr (char *str)
{
    int i = (int) strspn(str, " ");

    _shiftleft(str, i);
    for (i = strlen(str) - 1; i >= 0; i--)
	if (isspace((unsigned char) str[i])) str[i] = '\0';
	else break;
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
 * Returns: the output string.
 *
 */

char *switch_ext (char *targ, const char *src, char *ext)
{
    int i = dotpos(src);

    if (targ != src)
        strncpy(targ, src, i);
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
 * Returns: 1 if @c is found in @str, othewise 0.
 *
 */

int get_base (char *targ, const char *src, char c)
{
    int i, n = strlen(src);
	
    for (i=n-1; i>=0; i--) {
	if (src[i] == c) {
	    strncpy(targ, src, i+1);
	    targ[i+1] = '\0';
	    return 1;
	}
    }
    return 0;
}

/**
 * top_n_tail:
 * @str: the string to process.
 *
 * Drop leading space and trailing space and newline from string,
 * then replace a trailing backslash (if any) with a space.
 * 
 * Returns: 1 if a trailing backslash was found, otherwise 0.
 *
 */

int top_n_tail (char *str)
{
    int i;

    /* chop trailing space */
    i = strlen(str) - 1;
    while (isspace((unsigned char) str[i])) 
	str[i--] = '\0';
    /* drop leading spaces, also possible questionmark */
    i = 0;
    while (isspace((unsigned char) str[i]) || str[i] == '?') 
	i++;
    if (i) _shiftleft(str, i);
    /* then replace backslash, if present */
    if (str[strlen(str) - 1] == '\\') {
	str[strlen(str) - 1] = ' ';
	return 1;
    }
    return 0;
}  

/**
 * compress_spaces:
 * @str: the string to process.
 *
 * Reduce multiple contiguous space characters to single spaces
 * within @str.
 * 
 */

void compress_spaces (char *str)
{
    int i, j, n = strlen(str);

    for (i=0; i<n; i++) {
	if (str[i] == ' ' && str[i+1] == ' ') {
	    n--;
	    for (j=i+1; j<=n; j++) str[j] = str[j+1];
	    i--;
	}
    }
}   

/**
 * pprintf:
 * @prn: gretl printing struct.
 * @template: as in printf().
 * @Varargs: arguments to be printed.
 *
 * Multi-purpose printing function: can output to stream, to buffer
 * or to nowhere (silently discarding the output), depending on
 * how @prn was initialized.
 * 
 * Returns: 0 on successful completion, 1 on memory allocation
 * failure.
 * 
 */

int pprintf (PRN *prn, const char *template, ...)
{
    va_list args;
    size_t blen;

    if (prn == NULL) return 0;

    if (prn->fp != NULL) {
	va_start(args, template);
	vfprintf(prn->fp, template, args);
	va_end(args);
	return 0;
    }

    if (strncmp(template, "@init", 5) == 0) {
	prn->bufsize = 2048;
	prn->buf = malloc(prn->bufsize);
#ifdef PRN_DEBUG
  	fprintf(stderr, "pprintf: malloc'd %d bytes at %p\n", prn->bufsize,  
		(void *) prn->buf); 
#endif
	if (prn->buf == NULL) return 1;
	memset(prn->buf, 0, 1);
	return 0;
    }

    if (prn->buf == NULL) return 1;
    blen = strlen(prn->buf);
    if (prn->bufsize - blen < 1024) { 
	char *tmp;

#ifdef PRN_DEBUG
 	fprintf(stderr, "%d bytes left\ndoing realloc(%p, %d)\n",
 		prn->bufsize - blen, prn->buf, 2 * prn->bufsize);
#endif
	prn->bufsize *= 2; 
	tmp = realloc(prn->buf, prn->bufsize); 
	if (tmp == NULL) return 1;
	prn->buf = tmp;
#ifdef PRN_DEBUG
 	fprintf(stderr, "realloc: prn->buf is %d bytes at %p\n",
 		prn->bufsize, (void *) prn->buf);
#endif
	memset(prn->buf + blen, 0, 1);
    }
    va_start(args, template);
#ifdef PRN_DEBUG
     fprintf(stderr, "printing at %p\n", (void *) (prn->buf + blen));
#endif
    vsprintf(prn->buf + blen, template, args);
    va_end(args);
#ifdef PRN_DEBUG
     fprintf(stderr, "printed %d byte(s)\n", strlen(prn->buf) - blen);
#endif

    return 0;
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
 *
 */

char *safecpy (char *targ, const char *src, int n)
{
    strncpy(targ, src, n);
    targ[n] = 0;
    return targ;
}

/**
 * doing_nls:
 *
 * Returns: 1 if NLS translation is in effect, 0 otherwise.
 *
 */

int doing_nls (void)
{
#ifdef ENABLE_NLS
    static int called, nls;

    if (!called) {
	nls = (strcmp("/_File", _("/_File")) != 0);
	called = 1;
    }
    return nls;
#else
    return 0;
#endif
}

/**
 * get_local_decpoint:
 *
 * Returns: the character representing a decimal point in the current locale.
 *
 */

static int decpoint;

int reset_local_decpoint (void)
{
#ifdef ENABLE_NLS
    struct lconv *lc;

    lc = localeconv();
    decpoint = *lc->decimal_point;
    return decpoint;
#else
    return '.';
#endif
}


int get_local_decpoint (void)
{
#ifdef ENABLE_NLS
    if (decpoint == 0) decpoint = reset_local_decpoint();
    return decpoint;
#else
    return '.';
#endif
}

/**
 * obs_str_to_double:
 * @obs: string representation of observation number.
 *
 * Returns: the floating-point counterpart of @obs.
 */

double obs_str_to_double (const char *obs)
{
    char tmp[9];
    static int decpoint;

    if (decpoint == 0) decpoint = get_local_decpoint();

    strcpy(tmp, obs);
    charsub(tmp, ':', decpoint);
    return atof(tmp);
}

/**
 * colonize_obs:
 * @obs: string representation of observation number.
 *
 * Converts a decimal point in @obs to a colon.  Locale sensitive.
 *
 * Returns: the (possibly) modified obs string.
 */

char *colonize_obs (char *obs)
{
    static int decpoint;

    if (decpoint == 0) decpoint = get_local_decpoint();

    charsub(obs, decpoint, ':');
    if (decpoint != '.') charsub(obs, '.', ':');
    return obs;
}
