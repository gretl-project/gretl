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

/* embryonic form of a "plot" block command */

#include "libgretl.h"
#include "gretl_plot.h"

static gretlopt plotopt;

static const char *
get_plot_word_and_advance (const char *s, char *word, size_t maxlen)
{
    size_t i = 0;

    while (isspace(*s)) s++;

    *word = '\0';

    while (*s && !isspace(*s)) {
	if (i < maxlen) {
	    word[i++] = *s;
	}
	s++;
    }

    word[i] = '\0';
    s += strspn(s, " ");

    return (*word != '\0')? s : NULL;
}

int gretl_plot_append_line (const char *s, const DATASET *dset)
{
    char word[32];
    int err = 0;

    fprintf(stderr, "gretl_plot_append_line: '%s'\n", s);
    s = get_plot_word_and_advance(s, word, 32);

    if (!strcmp(word, "option")) {
	fprintf(stderr, "got option '%s'\n", s);
    } else if (!strcmp(word, "data")) {
	fprintf(stderr, "got data source '%s'\n", s);
    } else {
	fprintf(stderr, "got word '%s'; s='%s'\n", word, s);
    }

    return err;
}

int gretl_plot_finalize (const char *s, gretlopt opt)
{
    int err = 0;

    fprintf(stderr, "gretl_plot_finalize: '%s'\n", s);
    return err;
}
