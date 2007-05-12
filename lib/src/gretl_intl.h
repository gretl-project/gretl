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

#ifndef GRETL_INTL_H
#define GRETL_INTL_H

void gretl_push_c_numeric_locale (void);

void gretl_pop_c_numeric_locale (void);

int doing_nls (void);

int reset_local_decpoint (void);

int get_local_decpoint (void);

char *iso_to_ascii (char *s);

char *get_month_name (char *mname, int m);

#ifdef ENABLE_NLS

char *iso_gettext (const char *msgid);

char *maybe_iso_gettext (const char *msgid);

void set_gretl_charset (const char *s);

const char *get_gretl_charset (void);

const char *get_gnuplot_charset (void);

int iso_latin_version (void);

char *sprint_l2_to_html (char *targ, const char *s, size_t len);

char *sprint_l2_to_ascii (char *targ, const char *s, size_t len);

char *sprint_html_to_l2 (char *targ, const char *s);

int print_as_html (const char *s, FILE *fp);

int print_as_locale (const char *s, FILE *fp);

char *utf8_to_latin (const char *s);

int get_utf_width (const char *str, int width);

int get_translated_width (const char *str);

void check_for_console (PRN *prn);

void console_off (void);

# define UTF_WIDTH(s, w) get_utf_width(s, w) 
# define TRANSLATED_WIDTH(s) get_translated_width(s)

#else

# define UTF_WIDTH(s, w) w
# define TRANSLATED_WIDTH(s) strlen(s)

#endif /* ENABLE_NLS */

#endif /* GRETL_INTL_H */
