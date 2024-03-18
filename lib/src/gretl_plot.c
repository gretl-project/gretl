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

/* Implementation of "plot" block command. Current set-up:

   The block must start with "plot" (possibly preceded by the GUI
   "savename <- " apparatus). This line also takes a required
   parameter: the name of the variable that supplies the data to
   be plotted. This can be a single series, a list or a matrix.

   Optional elements:

   * zero or more "option" lines, holding a single option each

   * zero or more "options" lines holding one or more options
   each

   * zero or more "literal" lines, to be passed literally to
   gnuplot

   * zero or more "printf" lines (which turn into "literal"
   lines once the printf is cashed out)

   Option lines take the form "option flag" or "option flag=val"
   depending on whether the option takes a paraneter or not.
   The usual double-dash before the option flag is not required.

   The block ends with "end plot". The --output=whatever option
   may be appended to the ending line.
*/

#include "libgretl.h"
#include "uservar.h"
#include "usermat.h"
#include "libset.h"
#include "gretl_plot.h"

#define PDEBUG 0

#if PDEBUG
# include "gretl_typemap.h"
#endif

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
    reset_effective_plot_ci();
}

static int no_data_plot (gretlopt opt)
{
    FILE *fp = NULL;
    int i, np = 0;
    int err = 0;

    /* check that there's something to plot */
    for (i=0; i<plot.nlines; i++) {
	if (!strncmp(plot.lines[i], "plot ", 5) ||
	    !strncmp(plot.lines[i], "splot ", 6)) {
	    np++;
	}
    }

    if (np == 0) {
	gretl_errmsg_set("plot: nothing to plot");
	return E_ARGS;
    }

    fp = open_plot_input_file(PLOT_USER, 0, &err);

    if (!err) {
	for (i=0; i<plot.nlines; i++) {
	    fputs(plot.lines[i], fp);
	    fputc('\n', fp);
	}
	err = finalize_plot_input_file(fp);
    }

    return err;
}

/* In the following function we turn the array of "literal"
   lines from "plot" into the form

   { literal 1; literal 2; ...}

   which is the form wanted by the gnuplot() function at
   present. For future reference, it would be more efficient to
   revise gnuplot() so that it accepts an array of strings as
   input.
*/

static char *construct_literal_arg (int *err)
{
    char *literal = NULL;
    size_t litlen = 0;
    int i;

    for (i=0; i<plot.nlines; i++) {
	litlen += strlen(plot.lines[i]);
    }

    if (litlen > 0) {
	litlen += 2 * plot.nlines + 4;
	literal = calloc(litlen, 1);
	if (literal == NULL) {
	    *err = E_ALLOC;
	} else {
	    strcpy(literal, "{ ");
	    for (i=0; i<plot.nlines; i++) {
		strcat(literal, plot.lines[i]);
		strcat(literal, "; ");
	    }
	    strcat(literal, "}");
	}
    }

    return literal;
}

static int execute_plot (const DATASET *dset, gretlopt opt)
{
    int *list = NULL;
    char *literal = NULL;
    int free_list = 0;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "plot datasource = %s (type %s)\n",
	    plot.datasource, gretl_type_get_name(plot.datatype));
    fprintf(stderr, "number of literal lines: %d\n", plot.nlines);
#endif

    plot.opt |= opt;

    if (plot.datatype == GRETL_TYPE_LIST) {
	list = get_list_by_name(plot.datasource);
    } else if (plot.datatype == GRETL_TYPE_SERIES) {
	/* convert to singleton list */
	list = gretl_list_new(1);
	list[1] = current_series_index(dset, plot.datasource);
	free_list = 1;
    }

    if (plot.datatype != 0 && plot.nlines > 0) {
	literal = construct_literal_arg(&err);
    }

    if (!err && list != NULL && list[0] == 1) {
	/* default to time-series plot */
	if (dataset_is_time_series(dset)) {
	    plot.opt |= OPT_T;
	}
    }

    if (!err) {
	if (plot.datatype == 0) {
	    err = no_data_plot(plot.opt);
	} else if (plot.datatype == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = get_matrix_by_name(plot.datasource);

	    err = matrix_plot(m, NULL, literal, plot.opt);
	} else {
	    err = gnuplot(list, literal, dset, plot.opt);
	}
#if PDEBUG
	fprintf(stderr, "actual exec returned %d\n", err);
#endif
    }

    free(literal);

    if (free_list) {
	free(list);
    }

    return err;
}

/* we'll accept a single series, list or matrix as the
   "datasource" for a plot block */

static int check_plot_data_source (const char *s,
				   const DATASET *dset)
{
    int err = 0;

    if (plot.datatype != 0) {
	/* duplicated data spec */
	err = E_DATA;
    } else if (s == NULL) {
	/* no param given */
	return 0;
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

static void maybe_unquote_param (char *s)
{
    if (*s == '"') {
	int n;

	shift_string_left(s, 1);
	n = strlen(s);
	if (s[n-1] == '"') {
	    s[n-1] = '\0';
	}
    }
}

static int check_plot_option (const char *s)
{
    const char *p = NULL;
    char *param = NULL;
    OptStatus status;
    gretlopt opt;
    int err = 0;

    if (!strncmp(s, "--", 2)) {
	/* tolerate (non-required) option dashes */
	s += 2;
    }

#if PDEBUG
    fprintf(stderr, "check_plot_option: '%s'\n", s);
#endif

    if ((p = strchr(s, '=')) != NULL) {
	char *flag = gretl_strndup(s, p - s);

	param = gretl_strdup(p + 1);
#if PDEBUG
	fprintf(stderr, " flag = '%s', param = '%s'\n", flag, param);
#endif
	opt = valid_long_opt(PLOT, flag, &status);
	free(flag);
    } else {
	opt = valid_long_opt(PLOT, s, &status);
    }

    if (opt == OPT_NONE) {
	fprintf(stderr, "plot option: got OPT_NONE for '%s'\n", s);
	err = E_BADOPT;
    } else if (status == OPT_NO_PARM && param != NULL) {
	fprintf(stderr, "plot option: got spurious param: '%s'\n", s);
	err = E_BADOPT;
    } else if (status == OPT_NEEDS_PARM && param == NULL) {
	fprintf(stderr, "plot option: missing param: '%s'\n", s);
	err = E_ARGS;
    } else {
	plot.opt |= opt;
	if (param != NULL) {
	    if (opt == OPT_W) {
		/* --font=whatever */
		maybe_unquote_param(param);
	    }
	    err = push_option_param(PLOT, opt, param);
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

static int plot_printf (const char *s, const DATASET *dset)
{
    char *genline;
    char *genout;
    int err = 0;

    genline = g_strdup_printf("sprintf(%s)", s);
    gretl_push_c_numeric_locale();
    genout = generate_string(genline, (DATASET *) dset, &err);
    gretl_pop_c_numeric_locale();

    if (err) {
	fprintf(stderr, "plot_printf error: genline='%s'\n",
		genline);
    } else {
	strings_array_add(&plot.lines, &plot.nlines, genout);
	free(genout);
    }

    g_free(genline);

    return err;
}

/* This could go into strutils.c? */

static int strings_array_splice (char ***pS1, int *pn1,
				 char **S2, int n2)
{
    char **S = NULL;
    int i, j, n;
    int err = 0;

    if (S2 == NULL || n2 == 0) {
	return 0;
    }

    n = *pn1 + n2;

    S = realloc(*pS1, n * sizeof *S);
    if (S == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n2; i++) {
	S[*pn1 + i] = NULL;
    }

    for (i=0; i<n2 && !err; i++) {
	j = *pn1 + i;
	if (S2[i] != NULL) {
	    S[j] = gretl_strdup(S2[i]);
	    if (S[j] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    *pS1 = S;
    *pn1 = n;

    return err;
}

static int plot_append_strings (char ***pS, int *pns, const char *s)
{
    gretl_array *A = get_array_by_name(s);
    int ns = 0, err = 0;

    if (A == NULL) {
	err = E_DATA;
    } else if (gretl_array_get_type(A) != GRETL_TYPE_STRINGS) {
	err = E_TYPES;
    } else {
	char **S = gretl_array_get_strings(A, &ns);

	if (S != NULL) {
	    err = strings_array_splice(pS, pns, S, ns);
	}
    }

    return err;
}

enum {
    LAST_FIELD = 1,
    EMPTY_OK
};

static const char *
get_plot_field_and_advance (const char *s, char *field,
			    size_t maxlen, int flag,
			    int *err)
{
    const char *p = s;
    int quoted = 0;
    size_t i = 0;

    while (isspace(*s)) s++;

    *field = '\0';

    while (*s) {
	if (*s == '"') {
	    quoted = !quoted;
	} else if (!quoted && isspace(*s)) {
	    break;
	}
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

    if (*field == '\0' && flag != EMPTY_OK) {
	*err = E_ARGS;
    } else if (flag == LAST_FIELD && *s != '\0') {
	gretl_errmsg_sprintf(_("Parse error at unexpected token '%s'"), s);
	*err = E_PARSE;
    }

    return s;
}

/**
 * gretl_plot_append_line:
 * @s: the line to append.
 * @dset: pointer to dataset.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_append_line (const char *s, const DATASET *dset)
{
    char field[256];
    int err = 0;

    if (!plot.in_progress) {
	return E_DATA;
    }

#if PDEBUG
    fprintf(stderr, "gretl_plot_append_line: '%s'\n", s);
#endif

    s = get_plot_field_and_advance(s, field, 16, 0, &err);

    if (!strcmp(field, "option")) {
	s = get_plot_field_and_advance(s, field, sizeof field, LAST_FIELD, &err);
	if (!err) {
	    err = check_plot_option(field);
	    if (err) {
		fprintf(stderr, "Invalid plot option '%s'\n", field);
	    }
	}
    } else if (!strcmp(field, "options")) {
	int flag = 0;

	while (1) {
	    s = get_plot_field_and_advance(s, field, sizeof field, flag, &err);
	    if (err || *field == '\0') {
		break;
	    } else {
		err = check_plot_option(field);
	    }
	    flag = EMPTY_OK;
	}
    } else if (!strcmp(field, "literal")) {
	err = strings_array_add(&plot.lines, &plot.nlines, s);
    } else if (!strcmp(field, "printf")) {
	err = plot_printf(s, dset);
    } else if (!strcmp(field, "strings")) {
	err = plot_append_strings(&plot.lines, &plot.nlines, s);
    } else {
	fprintf(stderr, "plot: invalid field '%s'\n", field);
	err = E_PARSE;
    }

    if (err) {
	clear_plot();
    }

    return err;
}

/**
 * gretl_plot_start:
 * @param: name of the variable supplying the data.
 * @dset: pointer to dataset.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_start (const char *param, const DATASET *dset)
{
    int err = 0;

#if PDEBUG
    fprintf(stderr, "gretl_plot_start: '%s'\n", param);
#endif

    if (plot.in_progress) {
	clear_plot();
	err = E_DATA;
    } else {
	err = check_plot_data_source(param, dset);
    }

    if (!err) {
	plot.in_progress = 1;
	/* for option-param handling in graphing.c */
	set_effective_plot_ci(PLOT);
    }

    return err;
}

/**
 * gretl_plot_finalize:
 * @s: unused.
 * @dset: pointer to dataset.
 * @opt: may contain OPT_U (output spec).
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_finalize (const char *s, const DATASET *dset,
			 gretlopt opt)
{
    int err = 0;

#if PDEBUG
    fprintf(stderr, "gretl_plot_finalize: '%s'\n", s);
    debug_print_option_flags("end plot", opt);
#endif

    if (!gretl_in_gui_mode() && getenv("CLI_NO_PLOTS")) {
	; /* skip it */
    } else {
	err = execute_plot(dset, opt);
    }
    clear_plot();

    return err;
}
