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
#include "gretl_private.h"

#include <errno.h>
#include <time.h>

#if defined(USE_GTK2)
# include <glib.h>
#else
# include <fnmatch.h>
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

    if (str == NULL || *str == '\0') return 0;

    n = strlen(str);
    for (i=n-1; i>0; i--) {
	if (str[i] == SLASH) return i;
    }

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

    for (i=j=0; str[i] != '\0'; i++) {
	if (str[i] != c) {
	    str[j++] = str[i];
	}
    }
    str[j] = '\0';
}

/**
 * gretl_delete:
 * @str: the string to process.
 * @indx: the starting point for deleting characters.
 * @count: the number of characters to delete.
 *
 * Deletes @count characters from @str, starting at position @indx.
 *
 */

void gretl_delete (char *str, int indx, int count)
{
    size_t i, n = strlen(str);

    for (i=indx; i<=n-count; ++i) {
	str[i] = str[count+i];
    }
}

/**
 * haschar:
 * @c: the character to look for.
 * @s: the string to examine.
 *
 * Returns: the first position of @c in @s, or -1 if @c is not
 * found.
 *
 */

int haschar (char c, const char *s)
{
    int i = 0;

    while (*s) {
	if (*s++ == c) return i;
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
 * numeric_string:
 * @str: the string to examine.
 *
 * Returns: 1 if the given @str is numeric, otherwise 0.
 *
 */

int numeric_string (const char *str)
{
    char *test;
    int ret = 1;

    if (!strcmp(str, "inf") || !strcmp(str, "nan")) {
	/* could be variable names: they are not valid numbers */
	return 0;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    errno = 0;

    strtod(str, &test);
    if (*test != '\0' || errno == ERANGE) {
	ret = 0;
    }

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    return ret;
}

/**
 * ends_with_backslash:
 * @s: the string to examine.
 *
 * Returns: 1 if the last non-space character in @s is a backslash,
 * otherwise 0.
 *
 */

int ends_with_backslash (const char *s)
{
    int i, n = strlen(s);
    int bs = 0;

    for (i=n-1; i>=0; i--) {
	if (!isspace((unsigned char) s[i])) {
	    if (s[i] == '\\') bs = 1;
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
 */

void lower (char *str)
{
    while (*str) {
        if (isupper((unsigned char) *str)) *str = tolower(*str);
        str++;
    }
}

/**
 * gretl_strdup:
 * @src: the string to duplicate.
 *
 * Returns: an allocated copy of @src, or NULL on error.
 */

char *gretl_strdup (const char *src)
{
    char *targ = malloc(strlen(src) + 1);

    if (src != NULL) {
        strcpy(targ, src);
    }

    return targ;
}


/**
 * gretl_trunc:
 * @str: the string to truncate.
 * @n: the desired length of the truncated string.
 *
 * Truncates the given @str to the specified length.
 *
 */

void gretl_trunc (char *str, size_t n)
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
 * count_fields:
 * @str: the string to process.
 *
 * Returns: the number of space-separated fields in @str.
 *
 */

int count_fields (const char *s)
{
    int nf = 0;
    const char *p;

    if (s == NULL || *s == '\0') return 0;

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
	    if (*s) nf++;
	} else {
	    break;
	}
    }
	    
    return nf;
}

/**
 * shift_left:
 * @str: the string to process.
 * @move: the number of places to shift.
 *
 * Shifts the content of @str left by @move places, dropping
 * leading bytes as needed.
 *
 */

void shift_left (char *str, size_t move)
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
    int i = strspn(str, " \t");

    shift_left(str, i);

    for (i = strlen(str) - 1; i >= 0; i--) {
	if (isspace((unsigned char) str[i]) || str[i] == '\r') {
	    str[i] = '\0';
	}
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
    int i, n;

    if (src == NULL || *src == '\0') return 0;

    n = strlen(src);
    for (i=n-1; i>=0; i--) {
	if (src[i] == c) {
	    *targ = '\0';
	    strncat(targ, src, i+1);
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
    if (i) shift_left(str, i);

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

    if (strchr(s, '"') != NULL) {
	/* don't mess with literals */
	return;
    }

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
	nls = (strcmp("/File/_Open data", _("/File/_Open data")) ||
	       strcmp("Test statistic", _("Test statistic")) ||
	       strcmp("annual", _("annual")));
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

char *real_get_obs_string (char *obs, int t, const DATAINFO *pdinfo, int full)
{
    if (pdinfo->markers && pdinfo->S != NULL) { 
	/* data marker strings present */
	strcpy(obs, pdinfo->S[t]);
    } else {
	if (full) {
	    ntodate_full(obs, t, pdinfo);
	} else {
	    ntodate(obs, t, pdinfo);
	}
    }

    if (!full) {
	if (strlen(obs) > 8) { 
	    char tmp[9];

	    if (isdigit(*obs) && isdigit(*(obs + 1))) {
		strcpy(tmp, obs + 2);
		strcpy(obs, tmp);
	    } else {
		*tmp = '\0';
		strncat(tmp, obs, 8);
		strcpy(obs, tmp);
	    }
	} 
    }

    return obs;
}

/**
 * get_obs_string:
 * @obs: char array big enough to hold the observation (OBSLEN).
 * @t: zero-based observation number.
 * @pdinfo: pointer to dataset information.
 *
 * Returns: the observation string corresponding to @t.
 */

char *get_obs_string (char *obs, int t, const DATAINFO *pdinfo)
{
    return real_get_obs_string(obs, t, pdinfo, 0);
}

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

/* fudges for strings that should not be in utf-8 under some 
   conditions: under gtk2, translations usually come out in
   utf-8 in the GUI, but when we're sending stuff to stderr,
   it should probably be in ISO-8859-N. */

#ifdef ENABLE_NLS

static int gretl_cset_maj;
static int gretl_cset_min;
# ifdef WIN32
static int gretl_cpage;
# endif

void set_gretl_charset (const char *s)
{
    char gretl_charset[32];
    const char *charset = NULL;
    int using_utf8 = 0;

# ifdef USE_GTK2
    using_utf8 = g_get_charset(&charset);
# else
    charset = s;
# endif

    *gretl_charset = '\0';

    if (!using_utf8 && charset != NULL && *charset != '\0') {
	char *p;

	strncat(gretl_charset, charset, 31);
	lower(gretl_charset);
	p = strstr(gretl_charset, "iso");
	if (p != NULL) {
	    char numstr[6];

	    while (*p && !isdigit((unsigned char) *p)) p++;
	    *numstr = '\0';
	    strncat(numstr, p, 4);
	    gretl_cset_maj = atoi(numstr);
	    if (strlen(p) > 4) {
		p += 4;
		while (*p && !isdigit((unsigned char) *p)) p++;
		gretl_cset_min = atoi(p);
	    }
	    
	    if (gretl_cset_maj < 0 || gretl_cset_maj > 9000) {
		gretl_cset_maj = gretl_cset_min = 0;
	    } else if (gretl_cset_min < 0 || gretl_cset_min > 30) {
		gretl_cset_maj = gretl_cset_min = 0;
	    }
	} 
# ifdef WIN32
	if (p == NULL) {
	    sscanf(gretl_charset, "cp%d", &gretl_cpage);
	}
# endif
    }
}

const char *get_gretl_charset (void)
{
    static char cset[12];

    if (gretl_cset_maj > 0 && gretl_cset_min > 0) {
	sprintf(cset, "ISO-%d-%d\n", gretl_cset_maj, gretl_cset_min);
	return cset;
    } 

# ifdef WIN32
    if (gretl_cpage > 0) {
	sprintf(cset, "CP%d\n", gretl_cpage);
	return cset;
    }
# endif

    return NULL;
}

const char *get_gnuplot_charset (void)
{
    static char gp_enc[12];

    if (gretl_cset_maj == 8859 && 
	(gretl_cset_min == 1 || 
	 gretl_cset_min == 2 ||
	 gretl_cset_min == 15)) {
	sprintf(gp_enc, "iso_%d_%d\n", gretl_cset_maj, gretl_cset_min);
	return gp_enc;
    } 

    return NULL;
}

int use_latin_2 (void)
{
    int l2 = gretl_cset_maj == 8859 && gretl_cset_min == 2;

# ifdef WIN32
    if (!l2) {
	l2 = gretl_cpage == 1250;
    }
# endif

    return l2;
}

char *iso_gettext (const char *msgid)
{
   char *ret;
   static int cli;
   static int iso_ok = -1;
   static const char *cset;

   /* the command-line program is "special": it doesn't emit
      utf-8 at all, so we omit the redundant switching of
      codesets */
   if (!strcmp(msgid, "@CLI_INIT")) {
       cli = 1;
       return NULL;
   }

   if (cli) { /* command line program: switch not required */
       return gettext(msgid);
   }

   if (iso_ok < 0) {
       cset = get_gretl_charset();
       if (cset == NULL) iso_ok = 0;
       else iso_ok = 1;
   }

   if (iso_ok) {
       bind_textdomain_codeset(PACKAGE, cset);
       ret = gettext(msgid);
       bind_textdomain_codeset(PACKAGE, "UTF-8");
   } else {
       ret = gettext(msgid);
   }

   return ret;
} 

static int printing_to_console;

char *maybe_iso_gettext (const char *msgid)
{
   if (printing_to_console) {
       return iso_gettext(msgid);
   } else {
       return gettext(msgid);
   }
} 

void check_for_console (PRN *prn)
{
    if (prn == NULL) return;

    if (prn->fp == stdout || prn->fp == stderr) {
	printing_to_console = 1;
    } else {
	printing_to_console = 0;
    }
}

void console_off (void)
{
    printing_to_console = 0;
}

#endif  

const char *print_time (const time_t *timep)
{
    static char timestr[48];
    struct tm *local;

    local = localtime(timep);

    strftime(timestr, 47, "%Y/%m/%d %H:%M", local);

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
    while (*buf) {
	if (*buf == '&') sz += 4;
	else if (*buf == '<') sz += 3;
	else if (*buf == '>') sz += 3;
	else if (*buf == '"') sz += 5;
	buf++;
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
	} else if (*buf == '"') {
	    strcpy(p, "&quot;");
	    p += 6;
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

static char x2c (char *s) 
{
    register char digit;

    digit = (s[0] >= 'A' ? ((s[0] & 0xdf) - 'A') + 10 : (s[0] - '0'));
    digit *= 16;
    digit += (s[1] >= 'A' ? ((s[1] & 0xdf) - 'A') + 10 : (s[1] - '0'));
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

#ifndef USE_GTK2

int
utf8_to_iso_latin_1 (unsigned char *out, int outlen, 
		     unsigned char *in, int inlen)
{
    unsigned char* outstart = out;
    unsigned char* outend = out + outlen;
    unsigned char* inend = in + inlen;
    unsigned char c;

    while (in < inend) {
        c = *in++;
        if (c < 0x80) {
            if (out >= outend) return -1;
            *out++ = c;
        }
        else if (((c & 0xFE) == 0xC2) && in < inend) {
            if (out >= outend) return -1;
            *out++ = ((c & 0x03) << 6) | (*in++ & 0x3F);
        }
        else return -2;
    }

    return out - outstart;
}

#endif

static char *real_iso_to_ascii (char *s, int latin)
{
    char *tmp, *p, *q;

    tmp = malloc(strlen(s) + 1);
    if (tmp == NULL) return NULL;

    p = tmp;
    q = s;

    if (latin != 1 && latin != 2) {
	/* fallback?? */
	latin = 1;
    }

    if (latin == 1) {
	while (*q) {
	    unsigned char c = *q;

	    if (c == '\t' || (c >= 32 && c <= 126)) {
		*p++ = c;
	    } else if (c >= 192 && c <= 198) {
		*p++ = 'A';
	    } else if (c == 199) {
		*p++ = 'C';
	    } else if (c >= 200 && c <= 203) {
		*p++ = 'E';
	    } else if (c >= 204 && c <= 207) {
		*p++ = 'I';
	    } else if (c == 208) {
		*p++ = 'D';
	    } else if (c == 209) {
		*p++ = 'N';
	    } else if (c >= 210 && c <= 214) {
		*p++ = 'O';
	    } else if (c == 216) {
		*p++ = 'O';
	    } else if (c >= 217 && c <= 220) {
		*p++ = 'U';
	    } else if (c == 221) {
		*p++ = 'Y';
	    } else if (c >= 224 && c <= 230) {
		*p++ = 'a';
	    } else if (c == 231) {
		*p++ = 'c';
	    } else if (c >= 232 && c <= 235) {
		*p++ = 'e';
	    } else if (c >= 236 && c <= 239) {
		*p++ = 'i';
	    } else if (c == 240) {
		*p++ = 'd';
	    } else if (c == 241) {
		*p++ = 'n';
	    } else if (c >= 242 && c <= 246) {
		*p++ = 'o';
	    } else if (c == 248) {
		*p++ = 'o';
	    } else if (c >= 249 && c <= 252) {
		*p++ = 'u';
	    } else if (c == 253) {
		*p++ = 'y';
	    }
	    q++;
	}
    } else if (latin == 2) {
	while (*q) {
	    unsigned char c = *q;

	    if (c == '\t' || (c >= 32 && c <= 126)) {
		*p++ = c;
	    }

#ifndef WIN32
	    if (c==161 || c==193 || c==194 || c==195 || c==196) {
		*p++ = 'A';
	    }
#else
	    if (c==165 || c==193 || c==194 || c==195 || c==196) {
		*p++ = 'A';
	    }
#endif
	    else if (c==198 || c==199 || c==200) {
		*p++ = 'C';
	    }
	    else if (c==207 || c==208) {
		*p++ = 'D';
	    }
	    else if (c==201 || c==202 || c==203 || c==204) {
		*p++ = 'E';
	    }
	    else if (c==205 || c==206) {
		*p++ = 'I';
	    }
#ifndef WIN32
	    else if (c==163 || c==165 || c==197) {
		*p++ = 'L';
	    }
#else
	    else if (c==163 || c==188 || c==197) {
		*p++ = 'L';
	    }
#endif
	    else if (c==209 || c==210) {
		*p++ = 'N';
	    }
	    else if (c==211 || c==212 || c==213 || c==214) {
		*p++ = 'O';
	    }
	    else if (c==192 || c==216) {
		*p++ = 'R';
	    }
#ifndef WIN32
	    else if (c==166 || c==169 || c==170) {
		*p++ = 'S';
	    }
#else
	    else if (c==138 || c==140 || c==170) {
		*p++ = 'S';
	    }
#endif
#ifndef WIN32
	    else if (c==171 || c==222) {
		*p++ = 'T';
	    }
#else
	    else if (c==141 || c==222) {
		*p++ = 'T';
	    }
#endif
	    else if (c==217 || c==218 || c==219 || c==220) {
		*p++ = 'U';
	    }
	    else if (c==221) {
		*p++ = 'Y';
	    }
#ifndef WIN32
	    else if (c==172 || c==174 || c==175) {
		*p++ = 'Z';
	    }
#else
	    else if (c==142 || c==143 || c==175) {
		*p++ = 'Z';
	    }
#endif
#ifndef WIN32
	    else if (c==177 || c==225 || c==226 || c==227 || c==228) {
		*p++ = 'a';
	    }
#else
	    else if (c==185 || c==225 || c==226 || c==227 || c==228) {
		*p++ = 'a';
	    }
#endif
	    else if (c==230 || c==231 || c==232) {
		*p++ = 'c';
	    }
	    else if (c==239 || c==240) {
		*p++ = 'd';
	    }
	    else if (c==233 || c==234 || c==235 || c==236) {
		*p++ = 'e';
	    }
	    else if (c==237 || c==238) {
		*p++ = 'i';
	    }
#ifndef WIN32
	    else if (c==179 || c==181 || c==229) {
		*p++ = 'l';
	    }
#else
	    else if (c==179 || c==190 || c==229) {
		*p++ = 'l';
	    }
#endif
	    else if (c==241 || c==242) {
		*p++ = 'n';
	    }
	    else if (c==243 || c==244 || c==245 || c==246) {
		*p++ = 'o';
	    }
	    else if (c==224 || c==248) {
		*p++ = 'r';
	    }
#ifndef WIN32
	    else if (c==182 || c==185 || c==186 || c==223) {
		*p++ = 's';
	    }
#else
	    else if (c==154 || c==156 || c==186 || c==223) {
		*p++ = 's';
	    }
#endif
#ifndef WIN32
	    else if (c==187 || c==254) {
		*p++ = 't';
	    }
#else
	    else if (c==157 || c==254) {
		*p++ = 't';
	    }
#endif
	    else if (c==249 || c==250 || c==251 || c==252) {
		*p++ = 'u';
	    }
	    else if (c==253) {
		*p++ = 'y';
	    }
#ifndef WIN32
	    else if (c==188 || c==190 || c==191) {
		*p++ = 'z';
	    }
#else
	    else if (c==158 || c==159 || c==191) {
		*p++ = 'z';
	    }
#endif
	    q++;
	}
    }

    *p = '\0';

    strcpy(s, tmp);
    free(tmp);

    return s;
}

char *iso_to_ascii (char *s) 
{
    return real_iso_to_ascii(s, 1);
}

char *get_month_name (char *mname, int m)
{
    struct tm mt;

    mt.tm_sec = 0;
    mt.tm_min = 0;
    mt.tm_hour = 0;
    mt.tm_mday = 1;
    mt.tm_mon = m - 1;
    mt.tm_year = 100;

    strftime(mname, 7, "%b", &mt);
    *mname = tolower(*mname);

#ifdef ENABLE_NLS
    real_iso_to_ascii(mname, gretl_cset_min);
#endif

    return mname;
}


#ifdef ENABLE_NLS

struct l2sym {
    int l2val;   /* iso-8859-2 (or CP1250) character code */
    int ucs2val; /* corresponding UCS-2 code */
};

#ifndef WIN32

/* iso-8859-2 */

static struct l2sym l2table[] = { 
    { 161, 260 }, /*  A; */
    { 162, 728 }, /*  '( */
    { 163, 321 }, /*  L/ */
    { 165, 317 }, /*  L< */
    { 166, 346 }, /*  S' */
    { 169, 352 }, /*  S< */
    { 170, 350 }, /*  S, */
    { 171, 356 }, /*  T< */
    { 172, 377 }, /*  Z' */
    { 174, 381 }, /*  Z< */
    { 175, 379 }, /*  Z. */
    { 177, 261 }, /*  a; */
    { 178, 731 }, /*  '; */
    { 179, 322 }, /*  l/ */
    { 181, 318 }, /*  l< */
    { 182, 347 }, /*  s' */
    { 183, 711 }, /*  '< */
    { 185, 353 }, /*  s< */
    { 186, 351 }, /*  s, */
    { 187, 357 }, /*  t< */
    { 188, 378 }, /*  z' */
    { 189, 733 }, /*  '" */
    { 190, 382 }, /*  z< */
    { 191, 380 }, /*  z. */
    { 192, 340 }, /*  R' */
    { 195, 258 }, /*  A( */
    { 197, 313 }, /*  L' */
    { 198, 262 }, /*  C' */
    { 200, 268 }, /*  C< */
    { 202, 280 }, /*  E; */
    { 204, 282 }, /*  E< */
    { 207, 270 }, /*  D< */
    { 208, 272 }, /*  D/ */
    { 209, 323 }, /*  N' */
    { 210, 327 }, /*  N< */
    { 213, 336 }, /*  O" */
    { 216, 344 }, /*  R< */
    { 217, 366 }, /*  U0 */
    { 219, 368 }, /*  U" */
    { 222, 354 }, /*  T, */
    { 224, 341 }, /*  r' */
    { 227, 259 }, /*  a( */
    { 229, 314 }, /*  l' */
    { 230, 263 }, /*  c' */
    { 232, 269 }, /*  c< */
    { 234, 281 }, /*  e; */
    { 236, 283 }, /*  e< */
    { 239, 271 }, /*  d< */
    { 240, 273 }, /*  d/ */
    { 241, 324 }, /*  n' */
    { 242, 328 }, /*  n< */
    { 245, 337 }, /*  o" */
    { 248, 345 }, /*  r< */
    { 249, 367 }, /*  u0 */
    { 251, 369 }, /*  u" */
    { 254, 355 }, /*  t, */
    { 255, 729 }  /*  '. */
};

#else

/* Windows codepage 1250 */

static struct l2sym l2table[] = { 
    { 128, 8364 }, /*  Eu */
    { 130, 8218 }, /*  .9 */
    { 132, 8222 }, /*  :9 */
    { 133, 8230 }, /*  .3 */
    { 134, 8224 }, /*  /- */
    { 135, 8225 }, /*  /= */
    { 137, 8240 }, /*  %0 */
    { 138,  352 }, /*  S< */
    { 139, 8249 }, /*  <1 */
    { 140,  346 }, /*  S' */
    { 141,  356 }, /*  T< */
    { 142,  381 }, /*  Z< */
    { 143,  377 }, /*  Z' */
    { 145, 8216 }, /*  '6 */
    { 146, 8217 }, /*  '9 */
    { 147, 8220 }, /*  "6 */
    { 148, 8221 }, /*  "9 */
    { 149, 8226 }, /*  sb */
    { 150, 8211 }, /*  -N */
    { 151, 8212 }, /*  -M */
    { 153, 8482 }, /*  TM */
    { 154,  353 }, /*  s< */
    { 155, 8250 }, /*  >1 */
    { 156,  347 }, /*  s' */
    { 157,  357 }, /*  t< */
    { 158,  382 }, /*  z< */
    { 159,  378 }, /*  z' */
    { 161,  711 }, /*  '< */
    { 162,  728 }, /*  '( */
    { 163,  321 }, /*  L/ */
    { 165,  260 }, /*  A; */
    { 170,  350 }, /*  S, */
    { 175,  379 }, /*  Z. */
    { 178,  731 }, /*  '; */
    { 179,  322 }, /*  l/ */
    { 185,  261 }, /*  a; */
    { 186,  351 }, /*  s, */
    { 188,  317 }, /*  L< */
    { 189,  733 }, /*  '" */
    { 190,  318 }, /*  l< */
    { 191,  380 }, /*  z. */
    { 192,  340 }, /*  R' */
    { 195,  258 }, /*  A( */
    { 197,  313 }, /*  L' */
    { 198,  262 }, /*  C' */
    { 200,  268 }, /*  C< */
    { 202,  280 }, /*  E; */
    { 204,  282 }, /*  E< */
    { 207,  270 }, /*  D< */
    { 208,  272 }, /*  D/ */
    { 209,  323 }, /*  N' */
    { 210,  327 }, /*  N< */
    { 213,  336 }, /*  O" */
    { 216,  344 }, /*  R< */
    { 217,  366 }, /*  U0 */
    { 219,  368 }, /*  U" */
    { 222,  354 }, /*  T, */
    { 224,  341 }, /*  r' */
    { 227,  259 }, /*  a( */
    { 229,  314 }, /*  l' */
    { 230,  263 }, /*  c' */
    { 232,  269 }, /*  c< */
    { 234,  281 }, /*  e; */
    { 236,  283 }, /*  e< */
    { 239,  271 }, /*  d< */
    { 240,  273 }, /*  d/ */
    { 241,  324 }, /*  n' */
    { 242,  328 }, /*  n< */
    { 245,  337 }, /*  o" */
    { 248,  345 }, /*  r< */
    { 249,  367 }, /*  u0 */
    { 251,  369 }, /*  u" */
    { 254,  355 }, /*  t, */
    { 255,  729 }  /*  '. */
};

#endif

static int l2_lookup (int c)
{
    int i, n = sizeof l2table / sizeof l2table[0];

    for (i=0; i<n; i++) {
	if (c == l2table[i].l2val) {
	    return l2table[i].ucs2val;
	}
    }

    return c;
}

static int ucs_lookup (int c)
{
    int i, n = sizeof l2table / sizeof l2table[0];

    for (i=0; i<n; i++) {
	if (c == l2table[i].ucs2val) {
	    return l2table[i].l2val;
	}
    }

    return c;
}

char *sprint_l2_to_html (char *targ, const char *s, size_t len)
{
    unsigned char c;
    char *p = targ;

    *p = '\0';

    while ((c = *s)) {
#ifndef WIN32
	if (c > 160) {
	    sprintf(p, "&#%d;", l2_lookup(c));
	    p = strchr(p, ';') + 1;
	} else if (c > 127) {
	    sprintf(p, "&#%d;", c);
	    p = strchr(p, ';') + 1;
	} else {
	    *p++ = c;
	}
#else
	if (c > 127) {
	    sprintf(p, "&#%d;", l2_lookup(c));
	    p = strchr(p, ';') + 1;
	} else {
	    *p++ = c;
	}
#endif
	s++;
	if (p - targ > len - 8) {
	    break;
	}
    }

    *p = '\0';

    return targ;
}

char *sprint_html_to_l2 (char *targ, const char *s)
{
    int u;
    char *p = targ;

    *p = '\0';

    while (*s) {
	if (sscanf(s, "&#%d;", &u)) {
	    *p++ = ucs_lookup(u);
	    s = strchr(s, ';') + 1;
	} else {
	    *p++ = *s++;
	}
    }

    *p = '\0';

    return targ;
}

int print_as_html (const char *s, FILE *fp)
{
    unsigned char c;
    int nread = 0;

    while ((c = *s)) {
        if (c > 160) {
            fprintf(fp, "&#%d;", l2_lookup(c));
        } else if (c > 127) {
            fprintf(fp, "&#%d;", c);
        } else {
            fputc(c, fp);
        }
        s++;
        nread++;
    }

    return nread;
}

int print_as_locale (const char *s, FILE *fp)
{
    int u, nwrote = 0;

    while (*s) {
	if (sscanf(s, "&#%d;", &u)) {
	    fputc(ucs_lookup(u), fp);
	    s = strchr(s, ';') + 1;
	} else {
	    fputc(*s++, fp);
	}
	nwrote++;
    }

    return nwrote;
}

#endif /* ENABLE_NLS */

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

char *append_dir (char *fname, const char *dir)
{
    size_t len;

    if (dir == NULL) return fname;

    len = strlen(fname);
    if (fname[len - 1] == '/' || fname[len - 1] == '\\')
        strcat(fname, dir);
    else {
        strcat(fname, SLASHSTR);
        strcat(fname, dir);
    }

    strcat(fname, SLASHSTR);

    return fname;
}

int build_path (const char *dir, const char *fname, char *path, 
		const char *ext)
{
    size_t len;

    if (dir == NULL || fname == NULL || path == NULL) return 1;

    *path = '\0';
    strcat(path, dir);
    len = strlen(path);
    if (len == 0) return 1;

    /* strip a trailing single dot */
    if (len > 1 && path[len-1] == '.' && 
	(path[len-2] == '/' || path[len-2] == '\\')) {
	    path[len-1] = '\0';
    }

    if (path[len-1] == '/' || path[len-1] == '\\') {
        /* dir is already properly terminated */
        strcat(path, fname);
    } else {
        /* otherwise put a separator in */
        strcat(path, SLASHSTR);
        strcat(path, fname);
    }

    if (ext != NULL) strcat(path, ext);

    return 0;
}

#if defined(USE_GTK2)

int *varname_match_list (const DATAINFO *pdinfo, const char *pattern)
{
    GPatternSpec *pspec;
    int *list = NULL;
    int i, n = 0;

    pspec = g_pattern_spec_new(pattern);
    for (i=1; i<pdinfo->v; i++) { 
	if (pdinfo->vector[i] &&
	    g_pattern_match_string(pspec, pdinfo->varname[i])) {
	    n++;
	}
    }

    if (n > 0) {
	list = malloc((n + 1) * sizeof *list);
	if (list != NULL) {
	    int j = 1;

	    list[0] = n;
	    for (i=1; i<pdinfo->v; i++) { 
		if (pdinfo->vector[i] &&
		    g_pattern_match_string(pspec, pdinfo->varname[i])) {
		    list[j++] = i;
		}
	    }
	}
    }

    g_pattern_spec_free(pspec);

    return list;
}

#elif defined(HAVE_FNMATCH_H) 

int *varname_match_list (const DATAINFO *pdinfo, const char *pattern)
{
    int *list = NULL;
    int i, n = 0;

    for (i=1; i<pdinfo->v; i++) { 
	if (pdinfo->vector[i] &&
	    fnmatch(pattern, pdinfo->varname[i], 0) == 0) {
	    n++;
	}
    }

    if (n > 0) {
	list = malloc((n + 1) * sizeof *list);
	if (list != NULL) {
	    int j = 1;

	    list[0] = n;
	    for (i=1; i<pdinfo->v; i++) { 
		if (pdinfo->vector[i] &&
		    fnmatch(pattern, pdinfo->varname[i], 0) == 0) {
		    list[j++] = i;
		}
	    }
	}
    }

    return list;
}

#endif
    



