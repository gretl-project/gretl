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
#include <errno.h>

#if defined(ENABLE_NLS) && defined(USE_GTK2)
#include <glib.h>
#endif

char gretl_tmp_str[MAXLEN];

/**
 * string_is_blank:
 * @s: the string to examine.
 *
 * Returns: 1 if the string is NULL, of length zero, or contains
 * nothing but space characters, otherwise returns 0.
 */

int string_is_blank (const char *s)
{
    if (s == NULL) return 1;

    while (*s) {
        if (!isspace((unsigned char) *s) && *s != CTRLZ) return 0;
        s++;
    }

    return 1;
}

/**
 * dot_atof:
 * @s: the string to convert.
 *
 * Returns: the double-precision numeric interpretation of s,
 * where the decimal point character is forced to be '.', 
 * regardless of the current locale.
 */

double dot_atof (const char *s)
{
    double x;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif
    x = atof(s);
#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
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
 *
 */

size_t dotpos (const char *str)
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
    int i = 0;

    while (*str) {
	if (*str++ == c) return i;
	i++;
    }

    return -1;
}

/**
 * lastchar:
 * @c: the character to look for.
 * @str: the string to examine.
 *
 * Returns: 1 if @c is the last character in @s, 0 otherwise
 *
 */

int lastchar (char c, const char *s)
{
    if (s == NULL || *s == 0) return 0;
    if (s[strlen(s) - 1] == c) return 1;
    return 0;
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
    char *test;
    extern int errno;
    int ret = 1;

    if (!strcmp(str, "inf") || !strcmp(str, "nan")) {
	/* could be variable names: they are not valid numbers */
	return 0;
    }

    errno = 0;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    strtod(str, &test);
    if (*test != '\0' || !strcmp(str, test) || errno == ERANGE) {
	ret = 0;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return ret;
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

void _esl_trunc (char *str, size_t n)
{
    if (n < strlen(str)) str[n] = 0;
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
	*str = '\0';
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
    for (i = strlen(str) - 1; i >= 0; i--) {
	if (isspace((unsigned char) str[i])) str[i] = '\0';
	else break;
    }
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
 * Returns: 1 if @c is found in @str, otherwise 0.
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
    int i, len;

    if (str == NULL || *str == 0 || 
	*str == '\n' || *str == '\r') return 0;

    len = strlen(str);

    /* chop any trailing space */
    for (i=len-1; i>=0; i--) {
	if (isspace((unsigned char) str[i])) str[i] = 0;
	else break;
    }

    if (*str == 0) return 0;
	
    /* drop any leading spaces, also possible questionmark */
    i = 0;
    while (isspace((unsigned char) str[i]) || str[i] == '?') i++;
    if (i) _shiftleft(str, i);

    /* then replace backslash, if present */
    len = strlen(str);
    if (str[len - 1] == '\\') {
	str[len - 1] = ' ';
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

void compress_spaces (char *s)
{
    char *p;

    if (s == NULL || *s == 0) return;

    p = s;
    while (*s) {
	if (*s == '\t') *s = ' '; /* trash tabs */
	if (*s == ' ') {
	    p = s + 1;
	    if (*p == 0) break;
	    while (*p == ' ') p++;
	    if (p - s > 1) memmove(s + 1, p, strlen(p) + 1);
	}
	s++;
    }
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
 *
 */

char *safecpy (char *targ, const char *src, int n)
{
    *targ = 0;
    strncat(targ, src, n);
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

/* fudges for strings that should not be in utf8 under some 
   conditions */

#ifdef ENABLE_NLS

char *iso_gettext (const char *msgid)
{
   char *ret;
   static int cli;

   /* when running command-line client, ensure that translated
      messages will not appear in utf8 */
   if (!strcmp(msgid, "@CLI_INIT")) {
       cli = 1;
       return NULL;
   }

   if (cli) { /* command line program */
       return gettext(msgid);
   }

   bind_textdomain_codeset(PACKAGE, "ISO-8859-1");
   ret = gettext(msgid);
   bind_textdomain_codeset(PACKAGE, "UTF-8");
   return ret;
} 

/* library global */
int printing_to_console;

char *maybe_iso_gettext (const char *msgid)
{
   if (printing_to_console) return iso_gettext(msgid);
   else return gettext(msgid);
} 

#endif  

const char *print_time (const time_t *timep)
{
    static char timestr[48];
    struct tm *local;

    local = localtime(timep);

    strftime(timestr, 47, "%c", local);

    return timestr;
}

int get_utf_width (const char *str, int width)
{
#if defined(ENABLE_NLS) && defined(USE_GTK2)
    width += strlen(str) - g_utf8_strlen(str, -1);
#endif

    return width;
}

char *gretl_xml_encode (char *buf)
{
    char *xmlbuf, *p;
    size_t sz = strlen(buf) + 1;

#ifdef XML_DEBUG
    fprintf(stderr, "gretl_xml_encode: original buffer size=%d\n", sz);
#endif

    p = buf;
    while (*buf++) {
	if (*buf == '&') sz += 4;
	else if (*buf == '<') sz += 3;
	else if (*buf == '>') sz += 3;
    }
    buf = p;

    xmlbuf = malloc(sz);
    if (xmlbuf == NULL) {
	sprintf(gretl_errmsg, _("out of memory in XML encoding"));
	return NULL;
    }
#ifdef XML_DEBUG
    fprintf(stderr, "gretl_xml_encode: malloc'd xmlbuf at size %d\n", sz);
#endif
    p = xmlbuf;
    while (*buf) {
	if (*buf == '&') {
	    strcpy(p, "&amp;");
	    p += 5;
	} else if (*buf == '<') {
	    strcpy(p, "&lt;");
	    p += 4;
	} else if (*buf == '>') {
	    strcpy(p, "&gt;");
	    p += 4;
	} else {
	    *p++ = *buf;
	}
	buf++;
    }
    xmlbuf[sz-1] = '\0';

#ifdef XML_DEBUG
    fprintf(stderr, "done gretl_xml_encode: xmlbuf='%s'\n", xmlbuf);
#endif

    return xmlbuf;
}

static char x2c (char *what) 
{
    register char digit;

    digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A') + 10 : (what[0] - '0'));
    digit *= 16;
    digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A') + 10 : (what[1] - '0'));
    return digit;
}

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

char *make_varname_unique (char *vname, int v, DATAINFO *pdinfo) 
{
    int i, j, conflict;
    size_t n = strlen(vname);
    const char *add = "abcdefghijklmnopqrstuvwxyz";

    if (n > 7) n = 7;

    for (j=0; j<26; j++) {
	conflict = 0;
	for (i=1; i<pdinfo->v; i++) {
	    if (i != v && !strcmp(vname, pdinfo->varname[i])) {
		conflict = 1;
		break;
	    }
	}
	if (!conflict) break;
	vname[n] = add[j];
	vname[n+1] = '\0';
    }

    return vname;
}


