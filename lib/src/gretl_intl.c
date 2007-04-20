/*
 *  Copyright (c) 2005 by Allin Cottrell
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

#include "libgretl.h"
#include "texprint.h"

#include <glib.h>

#ifdef ENABLE_NLS

static int numeric_c_locale_depth = 0;
static char *numeric_locale = NULL;

/**
 * gretl_push_c_numeric_locale:
 *
 * Description:  Saves the current %LC_NUMERIC locale and sets it to "C"
 * This way you can safely read write floating point numbers all in the
 * same format.  You should make sure that code between
 * gretl_push_c_numeric_locale() and gretl_pop_c_numeric_locale()
 * doesn't do any setlocale calls or locale may end up in a strange setting.
 * Also make sure to always pop the C numeric locale after you've pushed it.
 * The calls can be nested.
 **/

void gretl_push_c_numeric_locale (void)
{
    if (numeric_c_locale_depth == 0) {
	free(numeric_locale);
	numeric_locale = gretl_strdup(setlocale(LC_NUMERIC, NULL));
	setlocale(LC_NUMERIC, "C");
    }
    numeric_c_locale_depth++;
}

/**
 * gretl_pop_c_numeric_locale:
 *
 * Description:  Restores the LC_NUMERIC locale to what it was 
 * before the matching gretl_push_c_numeric_locale(). If these calls
 * were nested, then this is a no-op until we get to the most outermost
 * layer. Code in between these should not do any setlocale calls
 * to change the %LC_NUMERIC locale or things may come out very strange.
 **/

void gretl_pop_c_numeric_locale (void)
{
    if (numeric_c_locale_depth == 0) {
	return;
    }

    numeric_c_locale_depth--;

    if (numeric_c_locale_depth == 0) {
	setlocale(LC_NUMERIC, numeric_locale);
	free(numeric_locale);
	numeric_locale = NULL;
    }
}

#else

void gretl_push_c_numeric_locale (void)
{
    return;
}

void gretl_pop_c_numeric_locale (void)
{
    return;
}

#endif /* ENABLE_NLS */

/**
 * doing_nls:
 *
 * Returns: 1 if NLS translation is in effect, 0 otherwise.
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

static int decpoint;

/**
 * reset_local_decpoint:
 *
 * Returns: the character representing a decimal point in the current locale.
 */

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

/**
 * get_local_decpoint:
 *
 * Returns: the character representing a decimal point in the current locale.
 */

int get_local_decpoint (void)
{
#ifdef ENABLE_NLS
    if (decpoint == 0) {
	decpoint = reset_local_decpoint();
    }
    return decpoint;
#else
    return '.';
#endif
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
    const char *charset = NULL;
    char gretl_charset[32];
    int using_utf8 = 0;

    using_utf8 = g_get_charset(&charset);
    if (using_utf8) {
	set_tex_use_utf(1);
    }

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
	sprintf(cset, "ISO-%d-%d", gretl_cset_maj, gretl_cset_min);
	return cset;
    } 

# ifdef WIN32
    if (gretl_cpage > 0) {
	sprintf(cset, "CP%d", gretl_cpage);
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
       if (cset == NULL) {
	   iso_ok = 0;
       } else {
	   iso_ok = 1;
       }
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

#endif  /* ENABLE_NLS */

static void 
iso_to_ascii_translate (char *targ, const char *src, int latin)
{
    char *p;
    const char *q;

    p = targ;
    q = src;

    if (latin == 1) {
	while (*q) {
	    unsigned char c = *q;

	    if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
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

	    if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
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
}

static char *real_iso_to_ascii (char *s, int latin)
{
    char *tmp;

    tmp = malloc(strlen(s) + 1);
    if (tmp == NULL) {
	return NULL;
    }

    if (latin != 1 && latin != 2) {
	/* fallback?? */
	latin = 1;
    }

    iso_to_ascii_translate(tmp, s, latin);

    strcpy(s, tmp);
    free(tmp);

    return s;
}

char *iso_to_ascii (char *s) 
{
    return real_iso_to_ascii(s, 1);
}

#ifdef ENABLE_NLS

struct l2sym {
    int l2val;   /* iso-8859-2 (or CP1250) character code */
    int ucs2val; /* corresponding UCS-2 code */
};

# ifndef WIN32

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

# else

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

# endif

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
# ifndef WIN32
	if (c > 160) {
	    sprintf(p, "&#%d;", l2_lookup(c));
	    p = strchr(p, ';') + 1;
	} else if (c > 127) {
	    sprintf(p, "&#%d;", c);
	    p = strchr(p, ';') + 1;
	} else {
	    *p++ = c;
	}
# else
	if (c > 127) {
	    sprintf(p, "&#%d;", l2_lookup(c));
	    p = strchr(p, ';') + 1;
	} else {
	    *p++ = c;
	}
# endif
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

char *sprint_l2_to_ascii (char *targ, const char *s, size_t len)
{
    iso_to_ascii_translate(targ, s, 2);

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

    real_iso_to_ascii(mname, gretl_cset_min);

    return mname;
}

int get_utf_width (const char *str, int width)
{
    width += strlen(str) - g_utf8_strlen(str, -1);

    return width;
}

int get_translated_width (const char *str)
{
    int w = strlen(str);

    w += w - g_utf8_strlen(str, -1);

    return w;
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

    if (printing_to_standard_stream(prn)) {
	printing_to_console = 1;
    } else {
	printing_to_console = 0;
    }
}

void console_off (void)
{
    printing_to_console = 0;
}

#endif /* ENABLE_NLS */
