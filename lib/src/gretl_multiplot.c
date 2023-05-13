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

/* Support for the built-in "gridplot" command for producing
   multiple plots.
*/

#include "libgretl.h"
#include "uservar.h"
#include "gretl_multiplot.h"

#define GRID_DEBUG 0

typedef struct {
    gretlopt flag;
    int *target;
    int min;
    int max;
    int def;
} mplot_option;

static GPtrArray *multiplot;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

int gretl_multiplot_active (void)
{
    return multiplot != NULL && mp_collecting;
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

static int update_multiplot_sizes (gretlopt opt, int *changes)
{
    const mplot_option *mpo;
    int new_rows = 0;
    int new_cols = 0;
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
                if (i < 3 && k != *mpo->target) {
                    *mpo->target = k;
                    *changes += 1;
                } else if (mpo->flag == OPT_R) {
                    new_rows = k;
                } else if (mpo->flag == OPT_C) {
                    new_cols = k;
                }
            }
	}
    }

    /* If an update involves a change to rows or cols (but not both)
       the prior complementary dimension may now be invalid. So we
       set it to zero, meaning that it will be set automatically in
       multiplot_set_grid().
    */

    if (!err) {
        if (new_rows > 0) {
            /* got a rows specification */
            if (new_rows != mp_rows) {
                mp_rows = new_rows;
                *changes += 1;
                if (new_cols == 0) {
                    /* make cols automatic */
                    mp_cols = 0;
                }
            }
        }
        if (new_cols > 0) {
            /* got a cols specification */
            if (new_cols != mp_cols) {
                mp_cols = new_cols;
                *changes += 1;
                if (new_rows == 0) {
                    /* make rows automatic */
                    mp_rows = 0;
                }
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

static void gretl_multiplot_destroy (void)
{
    if (multiplot != NULL) {
        g_ptr_array_unref(multiplot);
        multiplot = NULL;
    }
    set_multiplot_defaults();
}

int gretl_multiplot_start (gretlopt opt)
{
    int err = 0;

    if (multiplot == NULL) {
	if (opt == OPT_NONE) {
	    set_multiplot_defaults();
	} else {
	    err = set_multiplot_sizes(opt);
	}
	if (!err) {
	    multiplot = g_ptr_array_new_with_free_func(g_free);
	    mp_collecting = 1;
	}
    } else {
        gretl_errmsg_set("gridplot: cannot be nested");
        err = E_DATA;
    }

    return err;
}

/* Append a plot specification, in @buf, to the @multiplot array */

int gretl_multiplot_add_plot (gchar *buf)
{
    if (multiplot != NULL && buf != NULL) {
        g_ptr_array_add(multiplot, buf);
        return 0;
    } else {
	gretl_errmsg_set("gretl_multiplot_add_plot: failed");
        return E_DATA;
    }
}

static int multiplot_set_grid (int n)
{
    int err = 0;

    if (mp_rows == 0 && mp_cols == 0) {
	/* fully automatic grid */
	mp_rows = ceil(sqrt((double) n));
	mp_cols = ceil((double) n / mp_rows);
    } else if (mp_rows == 0) {
	/* automatic rows */
        if (mp_cols > n) {
            mp_cols = n;
            mp_rows = 1;
        } else {
            mp_rows = ceil((double) n / mp_cols);
        }
    } else if (mp_cols == 0) {
	/* automatic cols */
        if (mp_rows > n) {
            mp_rows = n;
            mp_cols = 1;
        } else {
            mp_cols = ceil((double) n / mp_rows);
        }
    } else if (mp_rows * mp_cols < n) {
	gretl_errmsg_sprintf("Specified grid (%d by %d) is too small "
			     "for %d sub-plots", mp_rows, mp_cols, n);
	err = E_INVARG;
    } else if (mp_rows * mp_cols > n) {
	int ar = ceil(sqrt((double) n));
	int ac = ceil((double) n / ar);

	if (mp_rows * mp_cols > ar * ac) {
	    gretl_errmsg_sprintf("Specified grid (%d by %d) is too big "
				 "for %d sub-plots", mp_rows, mp_cols, n);
	    err = E_INVARG;
	}
    }

#if GRID_DEBUG
    fprintf(stderr, "multiplot_set_grid: %d x %d for n=%d\n",
	    mp_rows, mp_cols, n);
#endif

    return err;
}

/* Write a multiplot specification to file, either using
   the @multiplot struct or an array of individual plot
   specification strings, @S.
*/

static int output_multiplot_script (const char **S, int np)
{
    gchar *buf;
    int i, err = 0;
    FILE *fp;

    /* insure against segfault */
    if (S == NULL && multiplot == NULL) {
        fprintf(stderr, "output_multiplot_script: internal error!\n");
        return E_DATA;
    }

    fp = open_plot_input_file(PLOT_USER_MULTI, 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    fprintf(fp, "# grid_params: fontsize=%d, width=%d, height=%d, "
            "rows=%d, cols=%d, plots=%d\n", mp_fontsize, mp_width,
            mp_height, mp_rows, mp_cols, np);
    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", mp_rows, mp_cols);
    gretl_push_c_numeric_locale();
    for (i=0; i<np; i++) {
        fprintf(fp, "# gridplot: subplot %d\n", i);
        if (i > 0) {
            fputs("reset\n", fp);
        }
	if (S != NULL) {
	    fputs(S[i], fp);
	} else {
	    buf = g_ptr_array_index(multiplot, i);
	    fputs(buf, fp);
	}
    }
    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);
    err = finalize_plot_input_file(fp);

    return err;
}

/* respond to "end gridplot" */

int gretl_multiplot_finalize (gretlopt opt)
{
    int np, err = 0;

    if (multiplot == NULL) {
	gretl_errmsg_set("end multiplot: multiplot not started");
	return E_DATA;
    }

    mp_collecting = 0;
    np = multiplot->len;
    if (np > 0) {
	err = multiplot_set_grid(np);
    }

    if (np > 0 && !err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(NULL, multiplot->len);
    }

    gretl_multiplot_destroy();

    return err;
}

static void filter_multiplot_buffer (const char *buf,
				     int np, FILE *fp)
{
    char line[1024];

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (!strncmp(line, "# grid_params: ", 15)) {
	    /* substitute revised parameters in comment */
	    fprintf(fp, "# grid_params: fontsize=%d, width=%d, height=%d, "
		    "rows=%d, cols=%d, plots=%d\n", mp_fontsize, mp_width,
		    mp_height, mp_rows, mp_cols, np);
	} else if (!strncmp(line, "set multiplot", 13)) {
	    /* and revise multiplot layout spec */
	    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n",
		    mp_rows, mp_cols);
	} else {
	    /* simply transcribe */
	    fputs(line, fp);
	}
    }

    bufgets_finalize(buf);
}

static int revise_multiplot_script (const char *buf,
				    int filter,
				    int np)
{
    int err = 0;
    FILE *fp;

    fp = open_plot_input_file(PLOT_USER_MULTI, 0, &err);
    if (err) {
        return err;
    }

    if (filter) {
	filter_multiplot_buffer(buf, np, fp);
    } else {
	fputs(buf, fp);
    }

    err = finalize_plot_input_file(fp);

    return err;
}

/* read the parameters from an existing gridplot buffer */

static int retrieve_grid_params (const char *buf, int *np)
{
    int parms[5] = {0};
    char line[256];
    int n = 0;
    int i = 0;
    int err = 0;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf) && 1 < 10) {
        if (!strncmp(line, "# grid_params: ", 15)) {
            n = sscanf(line+15, "fontsize=%d, width=%d, height=%d, "
		       "rows=%d, cols=%d, plots=%d", &parms[0], &parms[1],
		       &parms[2], &parms[3], &parms[4], np);
            break;
        }
        i++;
    }

    bufgets_finalize(buf);

#if GRID_DEBUG
    fprintf(stderr, "retrieve params: got %d x %d for n=%d\n",
	    parms[3], parms[4], *np);
#endif

    if (n == n_mp_options + 1) {
	/* transcribe to options array */
	for (i=0; i<n_mp_options; i++) {
	    *(mp_options[i].target) = parms[i];
	}
    } else {
	gretl_errmsg_set("Failed to retrieve gridplot specification");
	err = E_DATA;
    }

    return err;
}

static int get_prior_plot_spec (gretlopt opt,
				const char **pbuf,
				gchar **pgbuf)
{
    const char *s;
    int err = 0;

    err = incompatible_options(opt, OPT_i | OPT_I);
    if (err) {
	return err;
    }

    if (opt & OPT_i) {
	/* read a buffer */
	s = get_optval_string(GRIDPLOT, OPT_i);
	*pbuf = get_string_by_name(s);
    } else if (opt & OPT_I) {
	/* read a file */
	gboolean ok;

	s = get_optval_string(GRIDPLOT, OPT_I);
	ok = g_file_get_contents(s, pgbuf, NULL, NULL);
	if (ok) {
	    *pbuf = *pgbuf;
	}
    }

    if (pbuf == NULL) {
        gretl_errmsg_set("Couldn't find an input specification");
        err = E_DATA;
    }

    return err;
}

/* Revise an existing gridplot buffer or command file,
   presumably obtained via "end gridplot" with the --output
   or --outbuf option or perhaps via the "standalone"
   usage of gridplot. This may just be a matter of selecting
   an output format, or it may involve changes to options
   such as font size or layout.
*/

int gretl_multiplot_revise (gretlopt opt)
{
    const char *argname;
    const char *buf = NULL;
    gchar *gbuf = NULL;
    gretlopt myopt = opt;
    int filter = 0;
    int np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    /* we need an incoming gridplot specification */
    err = get_prior_plot_spec(opt, &buf, &gbuf);
    if (!err) {
	/* extract the dimensions recorded in @buf */
	err = retrieve_grid_params(buf, &np);
    }
    if (err) {
        return err;
    }

    /* what sort of output is wanted? */
    argname = get_optval_string(GRIDPLOT, OPT_U);
#if GRID_DEBUG
    if (argname != NULL) {
        fprintf(stderr, "gretl_multiplot_revise: output '%s'\n",
                argname);
    }
#endif

    /* let the current @opt override previous choices */
    myopt &= ~(OPT_i | OPT_I | OPT_U);
    if (myopt) {
        err = update_multiplot_sizes(myopt, &filter);
    }
    if (!err) {
	err = multiplot_set_grid(np);
    }

#if GRID_DEBUG > 1
    fprintf(stderr, "*** here's the prior input (filter = %d) ***\n", filter);
    fputs(buf, stderr);
#endif

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
	revise_multiplot_script(buf, filter, np);
    }

    if (gbuf != NULL) {
	g_free(gbuf);
    }

    return err;
}

/* This supports "standalone" usage of gridplot to process an
   array of individual plot-specification strings.
*/

int gretl_multiplot_from_array (gretlopt opt)
{
    const char *argname;
    gretl_array *a = NULL;
    const char **S = NULL;
    int np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    argname = get_optval_string(GRIDPLOT, OPT_S);
    a = get_array_by_name(argname);
    if (a == NULL) {
	err = E_DATA;
    } else {
	S = (const char **) gretl_array_get_strings(a, &np);
	if (S == NULL) {
	    err = E_DATA;
	}
    }

    if (!err) {
	/* pick up any options */
	gretlopt myopt = opt;

	set_multiplot_defaults();
	myopt &= ~OPT_S;
	if (myopt) {
	    err = set_multiplot_sizes(myopt);
	}
	if (!err) {
	    err = multiplot_set_grid(np);
	}
    }

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(S, np);
    }

    return err;
}
