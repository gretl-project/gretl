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

/* Support for the commands "gpbuild" and "gridplot", which produce
   plots using gnuplot's "multiplot" apparatus. This facility was
   added in May/June 2023.
*/

#include "libgretl.h"
#include "uservar.h"
#include "usermat.h"
#include "gretl_gridplot.h"

#define GRID_DEBUG 0

typedef struct {
    gretlopt flag;
    int *target;
    int min;
    int max;
    int def;
} mplot_option;

static gretl_array *mp_array;
static gchar *array_name;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

/* called in graphing.c */

int gretl_gridplot_collecting (void)
{
    return mp_array != NULL && mp_collecting;
}

static const mplot_option mp_options[] = {
    { OPT_F, &mp_fontsize, 4, 24, 10 },
    { OPT_W, &mp_width,    200, 2048, 800 },
    { OPT_H, &mp_height,   200, 2048, 600 },
    { OPT_R, &mp_rows,     1, 12, 0 },
    { OPT_C, &mp_cols,     1, 12, 0 }
};

static int n_mp_options = G_N_ELEMENTS(mp_options);

static int set_multiplot_sizes (gretlopt opt)
{
    const mplot_option *mpo;
    int i, k, err = 0;

    for (i=0; i<n_mp_options && !err; i++) {
	mpo = &mp_options[i];
	if (opt & mpo->flag) {
	    k = get_optval_int(GRIDPLOT, mpo->flag, &err);
	    if (!err && (k < mpo->min || k > mpo->max)) {
		gretl_errmsg_set("gridplot: out-of-bounds option value");
		err = E_INVARG;
	    }
	    if (!err) {
		*mpo->target = k;
	    }
	}
    }

    return err;
}

static void set_multiplot_defaults (void)
{
    int i;

    for (i=0; i<n_mp_options; i++) {
	*(mp_options[i].target) = mp_options[i].def;
    }
}

/* Set up the strings array in which gpbuild will cumulate
   sub-plots, with checks that @param is not the name of
   an existing object of incompatible type, and is a
   valid gretl identifier.
*/

static int initialize_mp_array (const char *param, DATASET *dset)
{
    int err = 0;

    mp_array = get_strings_array_by_name(param);

    if (mp_array != NULL) {
	; /* OK, there's already a strings array of this name */
    } else if (gretl_is_user_var(param) ||
	       current_series_index(dset, param) >= 0) {
	/* there's already an object of this name, but of the wrong type */
	err = E_TYPES;
    } else {
	/* check validity of @param as identifier */
	err = check_identifier(param);
    }

    if (!err) {
	mp_array = gretl_array_new(GRETL_TYPE_STRINGS, 0, &err);
    }

    if (mp_array != NULL) {
	array_name = g_strdup(param);
	mp_collecting = 1;
#if GRID_DEBUG
	fprintf(stderr, "gpbuild: started strings array %p (%s)\n",
		(void *) mp_array, array_name);
#endif
    }

    return err;
}

void gretl_gridplot_clear (int err)
{
    if (mp_array != NULL) {
#if GRID_DEBUG
	fprintf(stderr, "gretl_gridplot_clear: err=%d, array length %d\n",
		err, gretl_array_get_length(mp_array));
#endif
	if (err) {
	    gretl_array_destroy(mp_array);
	} else {
	    err = user_var_add_or_replace(array_name,
					  GRETL_TYPE_STRINGS,
					  mp_array);
	}
        mp_array = NULL;
    }

    g_free(array_name);
    array_name = NULL;
    mp_collecting = 0;
}

/* called from interact.c: process_command_error() */

void gretl_gridplot_destroy (void)
{
    gretl_gridplot_clear(1);
}

/* This responds to the starting command for a "gpbuild" block */

int gretl_gridplot_start (const char *param, gretlopt opt,
			  DATASET *dset)
{
    int err = 0;

    if (mp_array == NULL) {
	err = initialize_mp_array(param, dset);
    } else {
        gretl_errmsg_set(_("gpbuild: cannot be nested"));
        err = E_DATA;
    }

    return err;
}

static int invalid_mp_error (int ci)
{
    gretl_errmsg_sprintf(_("%s: invalid (multiplot) plot specification"),
			 gretl_command_word(ci));
    return E_INVARG;
}

/* Append a plot specification to the @multiplot array. This is
   called in graphing.c, but only inside a gpbuild block.
*/

int gretl_gridplot_add_plot (gchar *buf)
{
    int err = 0;

    if (mp_array != NULL && buf != NULL) {
        if (strstr(buf, "set multiplot")) {
            err = invalid_mp_error(GPBUILD);
        } else {
	    gretl_array_append_string(mp_array, buf, 1);
        }
    } else {
	gretl_errmsg_set("gretl_gridplot_add_plot: failed");
        err = E_DATA;
    }

#if GRID_DEBUG
    fprintf(stderr, "gretl_gridplot_add_plot, err = %d\n", err);
#endif

    return err;
}

static int set_mp_grid (int n_plots)
{
    int err = 0;

    if (mp_rows > n_plots) {
	gretl_errmsg_sprintf("invalid --rows specification for %d plots", n_plots);
	err = E_INVARG;
    } else if (mp_cols > n_plots) {
	gretl_errmsg_sprintf("invalid --cols specification for %d plots", n_plots);
	err = E_INVARG;
    }

    if (!err) {
	if (mp_rows > 0) {
	    /* automatic cols value */
	    mp_cols = ceil(n_plots / (double) mp_rows);
	} else if (mp_cols > 0) {
	    /* automatic rows value */
	    mp_rows = ceil(n_plots / (double) mp_cols);
	} else {
	    /* fully automatic */
	    mp_rows = ceil(sqrt((double) n_plots));
	    mp_cols = ceil(n_plots / (double) mp_rows);
	}
    }

#if GRID_DEBUG
    fprintf(stderr, "set_mp_grid: %d x %d\n",  mp_rows, mp_cols);
#endif

    return err;
}

/* Write subplot buffer @buf to file, stripping out any
   "set term..." statements.
*/

static void filter_subplot (const char *buf, FILE *fp)
{
    size_t sz = 2048;
    char *line = calloc(sz, 1);

    bufgets_init(buf);
    while (safe_bufgets(&line, &sz, buf)) {
	if (strncmp(line, "set term", 8)) {
	    fputs(line, fp);
	}
    }
    bufgets_finalize(buf);
    free(line);
}

static int get_subplot_index (gretl_matrix *m, int i, int j)
{
    return (int) gretl_matrix_get(m, i, j) - 1;
}

/* Write a multiplot specification to file, drawing on the
   strings array @a.

   @np indicates the number of included plots and @maxp the
   total number of subplots available (the length of the
   relevant array).
*/

static int output_multiplot_script (gretl_array *a,
				    gretl_matrix *m,
				    int np, int maxp)
{
    const char *buf;
    int i, j, k, p;
    int err = 0;
    FILE *fp;

    fp = open_plot_input_file(PLOT_GRIDPLOT, 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    fprintf(fp, "# grid_params: plots=%d, fontsize=%d, width=%d, height=%d, ",
	    np, mp_fontsize, mp_width, mp_height);
    fprintf(fp, "rows=%d, cols=%d\n", mp_rows, mp_cols);
    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", mp_rows, mp_cols);
    gretl_push_c_numeric_locale();

    k = -1;
    p = 0;
    for (i=0; i<mp_rows; i++) {
	for (j=0; j<mp_cols; j++) {
	    if (m != NULL) {
		k = get_subplot_index(m, i, j);
	    } else {
		k++;
	    }
#if GRID_DEBUG
	    fprintf(stderr, "i=%d, j=%d, k=%d, maxp=%d\n", i, j, k, maxp);
#endif
	    if (k < 0) {
		fputs("set multiplot next\n", fp);
	    } else {
		buf = NULL;
		if (k < maxp) {
		    buf = gretl_array_get_data(a, k);
		}
		if (buf != NULL) {
		    if (i + j > 0) {
			fputs("reset\n", fp);
		    }
		    fprintf(fp, "# subplot %d\n", ++p);
		    if (strstr(buf, "set term")) {
			filter_subplot(buf, fp);
		    } else {
			fputs(buf, fp);
		    }
		}
	    }
	}
    }

    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);
    err = finalize_plot_input_file(fp);

    return err;
}

/* set multiplot layout using a matrix argument */

static int set_mp_layout (gretl_matrix **pm, int *np)
{
    const char *s = get_optval_string(GRIDPLOT, OPT_L);
    gretl_matrix *m = NULL;
    int err = 0;

    if (s != NULL) {
	m = get_matrix_by_name(s);
    }
    if (gretl_is_null_matrix(m)) {
	err = E_INVARG;
    } else {
	int i, n = m->rows * m->cols;
        int nonzero = 0;
	double mi;

	for (i=0; i<n; i++) {
	    mi = m->val[i];
	    if (na(mi) || mi != floor(mi) ||
		mi < 0 || mi > *np) {
		gretl_errmsg_set(_("Invalid layout specification"));
		err = E_INVARG;
		break;
	    } else if (mi != 0) {
                nonzero++;
            }
	}
	if (!err) {
	    *pm = m;
	    *np = nonzero;
	    mp_rows = m->rows;
	    mp_cols = m->cols;
	}
    }

    return err;
}

/* respond to "end gpbuild" */

int gretl_gridplot_finalize (gretlopt opt)
{
    int err = 0;

#if GRID_DEBUG
    fprintf(stderr, "gretl_gridplot_finalize\n");
#endif

    if (mp_array == NULL) {
	gretl_errmsg_set("end gpbuild: building not started");
	err = E_DATA;
    } else {
	gretl_gridplot_clear(0);
    }

    return err;
}

static int is_graphic_filename (const char *s)
{
    const char *exts[] = {
        ".png", ".pdf", ".svg", ".eps", ".fig", ".emf", ".html"
    };
    int i;

    for (i=0; i < G_N_ELEMENTS(exts); i++) {
        if (has_suffix(s, exts[i])) {
            return 1;
        }
    }

    return 0;
}

/* If we find that @s is not really a gnuplot commands
   buffer, try treating it as the name of file that
   may contain gnuplot commands. This requires ruling out
   the possibility that it's a graphic file (e.g. a PNG).
   If it looks plausible, return the contents of the file.
*/

static char *try_for_plot_filename (const char *s)
{
    char *ret = NULL;

    if (!is_graphic_filename(s)) {
        char fname[FILENAME_MAX];
        GretlFileType ft;

        strcpy(fname, s);
        ft = detect_filetype(fname, OPT_P);
        if (ft == GRETL_CSV) {
            /* generic text data: might be gnuplot commands */
            gchar *tmp = NULL;
            gsize sz = 0;

            g_file_get_contents(fname, &tmp, &sz, NULL);
            if (tmp != NULL) {
                ret = gretl_strdup(tmp);
                g_free(tmp);
            }
        }
    }

    return ret;
}

/* Find and check the strings array identified by @argname.
   Called when gridplot is used independently of gpbuild.
*/

static int retrieve_plots_array (const char *argname,
				 gretl_array **pa,
				 int *pnp)
{
    gretl_array *a = NULL;
    int msg_set = 0;
    int err = 0;

    a = get_array_by_name(argname);
    if (a == NULL) {
	err = E_INVARG;
    } else {
	int i, n = gretl_array_get_length(a);
        char *content = NULL;
	char *buf;

	for (i=0; i<n; i++) {
	    buf = gretl_array_get_data(a, i);
            if (buf != NULL && strchr(buf, '\n') == NULL) {
                /* @buf can't be an actual plot buffer */
                content = try_for_plot_filename(buf);
                if (content != NULL) {
                    free(buf);
                    gretl_array_set_data(a, i, content);
                }
                buf = content;
            }
            /* check for embedded multiplots */
            if (buf == NULL || strstr(buf, "set multiplot")) {
		err = invalid_mp_error(GRIDPLOT);
		msg_set = 1;
		break;
	    }
	}
	if (!err) {
	    *pa = a;
	    *pnp = n;
	}
    }

    if (err && !msg_set) {
	gretl_errmsg_set("Didn't get a valid array of plot strings");
    }

    return err;
}

/* This supports production of a multiplot from the array of
   plot-specification strings identified by @param. Such an
   array may be produced by "gpbuild", or may be assembled
   manually by the user.
*/

int gretl_gridplot_from_array (const char *param, gretlopt opt)
{
    gretl_array *a = NULL;
    gretl_matrix *m = NULL;
    int maxp, np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    err = retrieve_plots_array(param, &a, &np);
    if (err) {
        return err;
    }

    maxp = np;
    set_multiplot_defaults();

#if GRID_DEBUG
    fprintf(stderr, "gretl_gridplot_from_array: np = %d, opt = %d\n",
	    np, opt);
#endif

    if (opt) {
	err = set_multiplot_sizes(opt);
	if (!err) {
	    if (opt & OPT_L) {
		err = set_mp_layout(&m, &np);
	    } else {
		err = set_mp_grid(np);
	    }
	}
    } else {
	/* no options supplied, figure the default grid */
	err = set_mp_grid(np);
    }

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(a, m, np, maxp);
    }

    return err;
}

int check_gridplot_options (gretlopt opt)
{
    int err;

    /* can't have both --output and --outbuf */
    err = incompatible_options(opt, OPT_U | OPT_b);
    if (!err) {
	/* can't have more than one of --rows, --cols, --layout */
	err = incompatible_options(opt, OPT_R | OPT_C | OPT_L);
    }

    return err;
}
