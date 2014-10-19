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

/* 2014-10-19: embryonic form of a "plot" block command. Right now
   it's assumed that a data source (series, list of matrix) is
   provided, but it could be generalized to handle the case of
   a fully user-defined plot.
*/
  
#include "libgretl.h"
#include "uservar.h"
#include "gretl_typemap.h"
#include "gretl_plot.h"

#define PDEBUG 1

typedef struct gretl_plot_ {
    int in_progress;
    gretlopt opt;
    char *datasource;
    GretlType datatype;
    char **lines;
    int nlines;
} gretl_plot;

gretl_plot plot;

static void clear_plot (void)
{
    plot.in_progress = 0;
    plot.opt = 0;
    free(plot.datasource);
    plot.datasource = NULL;
    plot.datatype = 0;
    strings_array_free(plot.lines, plot.nlines);
    plot.lines = NULL;
    plot.nlines = 0;
}

static int execute_plot (const DATASET *dset, gretlopt opt)
{
    int *list = NULL;
    char *literal = NULL;
    size_t litlen = 0;
    int free_list = 0;
    int i, err = 0;

    if (plot.datasource == NULL || plot.datatype == 0) {
	/* FIXME maybe this should be OK? */
	fprintf(stderr, "plot has no data source\n");
	return E_DATA;
    } else {
#if PDEBUG
	fprintf(stderr, "plot datasource = %s (type %s)\n", 
		plot.datasource, gretl_type_get_name(plot.datatype));
#endif
	if (plot.datatype == GRETL_TYPE_LIST) {
	    list = get_list_by_name(plot.datasource);
	} else if (plot.datatype == GRETL_TYPE_SERIES) {
	    list = gretl_list_new(1);
	    list[1] = current_series_index(dset, plot.datasource);
	    free_list = 1;
	}
    }

#if PDEBUG
    fprintf(stderr, "number of literal lines: %d\n", plot.nlines);
#endif

    for (i=0; i<plot.nlines; i++) {
#if PDEBUG
	fprintf(stderr, " %d: '%s'\n", i, plot.lines[i]);
#endif
	litlen += strlen(plot.lines[i]);
    }

    if (litlen > 0) {
	litlen += 2 * plot.nlines + 4;
	literal = malloc(litlen);
	if (literal == NULL) {
	    err = E_ALLOC;
	} else {
	    strcpy(literal, "{ ");
	    for (i=0; i<plot.nlines; i++) {
		strcat(literal, plot.lines[i]);
		strcat(literal, "; ");
	    }
	    strcat(literal, "}");
	}
#if PDEBUG
	fprintf(stderr, "composed literal: '%s'\n", literal);
#endif
    }	    

    if (!err) {
	err = gnuplot(list, literal, dset, plot.opt | opt);
    }

    if (free_list) {
	free(list);
    }

    return err;
}

static const char *
get_plot_field_and_advance (const char *s, char *field, size_t maxlen)
{
    const char *p = s;
    size_t i = 0;

    while (isspace(*s)) s++;

    *field = '\0';

    while (*s && !isspace(*s)) {
	if (i < maxlen) {
	    field[i++] = *s;
	} else {
	    fprintf(stderr, "plot field: overflow!\n");
	    *field = '\0';
	    return p;
	}
	s++;
    }

    field[i] = '\0';
    s += strspn(s, " \t\r\n");

    return s;
}

static int check_plot_data_source (const char *s, 
				   const DATASET *dset)
{
    int err = 0;

    /* we'll accept a single series, list or matrix as
       the plot's "datasource" */

    if (plot.datatype != 0) {
	/* duplicated "data" line */
	err = E_DATA;
    } else {
	int id = current_series_index(dset, s);

	if (id >= 0) {
	    plot.datatype = GRETL_TYPE_SERIES;
	} else {
	    GretlType type = 0;
    
	    user_var_get_value_and_type(s, &type);
	    if (type == GRETL_TYPE_MATRIX ||
		type == GRETL_TYPE_LIST) {
		plot.datatype = type;
	    }
	}
    }

    if (!err) {
	if (plot.datatype == 0) {
	    err = E_DATA;
	} else {
	    plot.datasource = gretl_strdup(s);
	}
    }

    return err;
}

static int check_plot_option (const char *s)
{
    const char *p = NULL;
    char *param = NULL;
    OptStatus status;
    gretlopt opt;
    int err = 0;

    if (!strncmp(s, "--", 2)) {
	s += 2;
    }

    if ((p = strchr(s, '=')) != NULL) {
	char *flag = gretl_strndup(s, p - s);

	param = gretl_strdup(p + 1);
	opt = valid_long_opt(GNUPLOT, flag, &status);
	free(flag);
    } else {
	opt = valid_long_opt(GNUPLOT, s, &status);
    }

    if (opt == OPT_NONE) {
	fprintf(stderr, "plot option: got OPT_NONE for '%s'\n", s);
	err = E_BADOPT;
    } else if (status == OPT_NO_PARM && param != NULL) {
	fprintf(stderr, "plot option: got spurious param: '%s'\n", s);
	err = E_BADOPT;
    } else if (status == OPT_NEEDS_PARM && param == NULL) {
	fprintf(stderr, "plot option: missing param: '%s'\n", s);
	err = E_BADOPT;
    } else {
	plot.opt |= opt;
	if (param != NULL) {
	    err = push_option_param(GNUPLOT, opt, param);
	    if (err) {
		fprintf(stderr, "plot option: error pushing param: '%s'\n", s);
	    } else {
		param = NULL;
	    }
	}
    }

    free(param);
    
    return err;
}

int gretl_plot_append_line (const char *s, const DATASET *dset)
{
    char field[64];
    int err = 0;

    s = get_plot_field_and_advance(s, field, 64);

    if (!strcmp(field, "plot")) {
	if (plot.in_progress) {
	    return E_PARSE;
	} else {
	    plot.in_progress = 1;
	    return 0;
	}
    } else if (!plot.in_progress) {
	return E_PARSE;
    }

    if (!strcmp(field, "option")) {
	s = get_plot_field_and_advance(s, field, 64);
	if (*field == '\0' || *s != '\0') {
	    err = E_PARSE;
	} else {
	    err = check_plot_option(field);
	}
	if (err) {
	    fprintf(stderr, "Invalid plot option '%s'\n", field);
	}
    } else if (!strcmp(field, "data")) {
	s = get_plot_field_and_advance(s, field, 64);
	if (*field == '\0' || *s != '\0') {
	    err = E_PARSE;
	} else {
	    err = check_plot_data_source(field, dset);
	}
	if (err) {
	    fprintf(stderr, "Invalid plot data source '%s'\n", field);
	}
    } else if (!strcmp(field, "literal")) {
	err = strings_array_add(&plot.lines, &plot.nlines, s);
    } else if (!strcmp(field, "printf")) {
	fprintf(stderr, "got printf '%s'\n", s);
	fprintf(stderr, " (not implemented yet)\n");
    } else {
	fprintf(stderr, "plot: invalid field '%s'\n", field);
	err = E_PARSE;
    }

    return err;
}

int gretl_plot_finalize (const char *s, const DATASET *dset, 
			 gretlopt opt)
{
    int err;

#if PDEBUG
    fprintf(stderr, "gretl_plot_finalize: '%s'\n", s);
#endif
    err = execute_plot(dset, opt);
    clear_plot();

    return err;
}
