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
static gretl_matrix *mp_layout;
static int mp_fontsize = 10;
static int mp_width = 800;
static int mp_height = 600;
static int mp_rows;
static int mp_cols;
static int mp_collecting;

int gretl_multiplot_collecting (void)
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

static int set_multiplot_sizes (int ci, gretlopt opt)
{
    const mplot_option *mpo;
    int i, k, err = 0;

    for (i=0; i<n_mp_options && !err; i++) {
	mpo = &mp_options[i];
	if (opt & mpo->flag) {
	    k = get_optval_int(ci, mpo->flag, &err);
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

/* respond to a validated case of the --layout=matrix option */

static void set_mp_layout (const gretl_matrix *m)
{
    if (mp_layout != NULL) {
	gretl_matrix_free(mp_layout);
    }
    mp_layout = gretl_matrix_copy(m);
    if (mp_layout != NULL) {
	if (mp_rows != m->rows || mp_cols != m->cols) {
	    mp_rows = m->rows;
	    mp_cols = m->cols;
	}
    }
}

static int update_multiplot_sizes (gretlopt opt)
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

    if (mp_layout != NULL) {
	gretl_matrix_free(mp_layout);
	mp_layout = NULL;
    }
    for (i=0; i<n_mp_options; i++) {
	*(mp_options[i].target) = mp_options[i].def;
    }
}

void gretl_multiplot_destroy (void)
{
    if (multiplot != NULL) {
        g_ptr_array_unref(multiplot);
        multiplot = NULL;
    }
    set_multiplot_defaults();
}

/* This responds to the starting command for a "gpbuild" block */

int gretl_multiplot_start (gretlopt opt)
{
    int err = 0;

    if (multiplot == NULL) {
	if (opt == OPT_NONE) {
	    set_multiplot_defaults();
	} else {
	    err = set_multiplot_sizes(GPBUILD, opt);
	}
	if (!err) {
	    multiplot = g_ptr_array_new_with_free_func(g_free);
	    mp_collecting = 1;
	}
    } else {
        gretl_errmsg_set(_("gridplot: cannot be nested"));
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

/* Append a plot specification to the @multiplot array */

int gretl_multiplot_add_plot (gchar *buf)
{
    int err = 0;

    if (multiplot != NULL && buf != NULL) {
        if (strstr(buf, "set multiplot")) {
            err = invalid_mp_error(GPBUILD);
        } else {
            g_ptr_array_add(multiplot, buf);
        }
    } else {
	gretl_errmsg_set("gretl_multiplot_add_plot: failed");
        err = E_DATA;
    }
#if GRID_DEBUG
    fprintf(stderr, "gretl_multiplot_add_plot, err = %d\n", err);
#endif
#if 0 /* let this be handled in interact.c ? */
    if (err) {
        gretl_multiplot_destroy();
    }
#endif

    return err;
}

static int multiplot_set_grid (int n)
{
    int err = 0;

#if GRID_DEBUG
    fprintf(stderr, "multiplot_set_grid: n=%d, prior size %d x %d\n",
	    n, mp_rows, mp_cols);
#endif

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
    fprintf(stderr, "multiplot_set_grid: set %d x %d\n",  mp_rows, mp_cols);
#endif

    return err;
}

static int get_subplot_index (int i, int j)
{
    return (int) gretl_matrix_get(mp_layout, i, j) - 1;
}

/* write a layout matrix into a plot file, in a form that can
   be easily reconstituted via generate_matrix()
*/

static void output_layout_matrix (gretl_matrix *m, FILE *fp)
{
    double mij;
    int i, j;

    fputs("layout={", fp);
    for (i=0; i<m->rows; i++) {
	for (j=0; j<m->cols; j++) {
	    mij = gretl_matrix_get(m, i, j);
	    fprintf(fp, "%d", (int) mij);
	    if (j < m->cols-1) {
		fputc(',', fp);
	    }
	}
	if (i < m->rows-1) {
	    fputc(';', fp);
	}
    }
    fputs("}\n", fp);
}

static void write_mp_spec_comment (int np, FILE *fp)
{
    fprintf(fp, "# grid_params: plots=%d, fontsize=%d, width=%d, height=%d, ",
	    np, mp_fontsize, mp_width, mp_height);
    if (mp_layout != NULL) {
	output_layout_matrix(mp_layout, fp);
    } else {
	fprintf(fp, "rows=%d, cols=%d\n", mp_rows, mp_cols);
    }
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

/* Write a multiplot specification to file, drawing on either
   @multiplot (g_ptr_array) or an array of individual plot
   specification strings, @S.

   @np indicates the number of included plots and @maxp the
   total number of subplots available (the length of the
   relevant array).
*/

static int output_multiplot_script (int ci, const char **S,
				    int np, int maxp)
{
    PlotType ptype;
    const char *buf;
    int i, j, k, p;
    int err = 0;
    FILE *fp;

    /* insure against segfault */
    if (S == NULL && multiplot == NULL) {
        fprintf(stderr, "output_multiplot_script: internal error!\n");
        return E_DATA;
    }

    ptype = ci == GRIDPLOT ? PLOT_GRIDPLOT : PLOT_GPBUILD;
    fp = open_plot_input_file(ptype, 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    write_mp_spec_comment(np, fp);

    fprintf(fp, "set multiplot layout %d,%d rowsfirst\n", mp_rows, mp_cols);
    gretl_push_c_numeric_locale();

    k = -1;
    p = 0;
    for (i=0; i<mp_rows; i++) {
	for (j=0; j<mp_cols; j++) {
	    if (mp_layout != NULL) {
		k = get_subplot_index(i, j);
	    } else {
		k++;
	    }
	    if (k < 0) {
		fputs("set multiplot next\n", fp);
	    } else {
		buf = NULL;
		if (S != NULL) {
		    if (k < maxp) {
			buf = S[k];
		    }
		} else if (k < multiplot->len) {
		    buf = g_ptr_array_index(multiplot, k);
		}
		if (buf != NULL) {
		    if (i+j > 0) {
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

static int maybe_set_mp_layout (int ci, int *np)
{
    const char *s = get_optval_string(ci, OPT_L);
    const gretl_matrix *m = NULL;
    int err = 0;

    if (s != NULL) {
	m = get_matrix_by_name(s);
    }
    if (m != NULL) {
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
	    *np = nonzero;
	}
	if (!err) {
#if GRID_DEBUG
	    gretl_matrix_print(m, "m, in maybe_set_mp_layout");
#endif
	    set_mp_layout(m);
	}
    }

    return err;
}

/* respond to "end gpbuild" */

int gretl_multiplot_finalize (gretlopt opt)
{
    int np, maxp, err = 0;

    if (multiplot == NULL) {
	gretl_errmsg_set("end gpbuild: building not started");
	return E_DATA;
    }

    mp_collecting = 0;
    np = maxp = multiplot->len;

    if (np > 0) {
	if (opt & OPT_L) {
	    err = maybe_set_mp_layout(GPBUILD, &np);
	} else {
	    err = multiplot_set_grid(np);
	}
    }

    if (np > 0 && !err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(GPBUILD, NULL, np, maxp);
    }

    gretl_multiplot_destroy();

    return err;
}

static int read_layout_matrix (const char *s, int *pr, int *pc, int *pn)
{
    gretl_matrix *m = NULL;
    int err = 0;

    s = strstr(s, "layout=") + 7;
    m = generate_matrix(s, NULL, &err);
    if (!err) {
	mp_layout = m;
	*pr = m->rows;
	*pc = m->cols;
	*pn += 2;
    }

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
	    if (strstr(line, "layout=") != NULL) {
		n = sscanf(line+15, "plots=%d, fontsize=%d, width=%d, height=%d",
			   np, &parms[0], &parms[1], &parms[2]);
		err = read_layout_matrix(line, &parms[3], &parms[4], &n);
	    } else {
		n = sscanf(line+15, "plots=%d, fontsize=%d, width=%d, height=%d, "
			   "rows=%d, cols=%d", np, &parms[0], &parms[1],
			   &parms[2], &parms[3], &parms[4]);
		break;
	    }
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

/* Given a gridplot buffer @buf, extract the individual
   subplots into a GLib pointer-array.
*/

static int disassemble_multiplot (const char *buf, int np)
{
    const char *e1, *e2, *q;
    GPtrArray *arr = NULL;
    gchar *subplot;
    const char *p = buf;
    int k, err = 0;

    arr = g_ptr_array_new_full(np, g_free);

    while ((p = strstr(p, "# subplot ")) != NULL) {
	k = atoi(p + 10);
#if GRID_DEBUG
	fprintf(stderr, "found subplot %d\n", k);
#endif
	e1 = strstr(p + 10, "set multiplot next");
	e2 = strstr(p + 10, "reset");
	if (e1 != NULL && e2 != NULL) {
	    q = e1 - e2 > 0 ? e2 : e1;
	} else if (e1 == NULL && e2 == NULL) {
	    q = strstr(p + 10, "unset multiplot");
	} else {
	    q = e1 != NULL ? e1 : e2;
	}
	if (q == NULL) {
	    err = 1;
	    break;
	}
        /* skip "# subplot ..." */
        p = strchr(p, '\n') + 1;
	subplot = g_strndup(p, q-p);
	g_ptr_array_insert(arr, k-1, subplot);
	p = q;
    }

    if (!err && arr->len != np) {
	err = 1;
    }
    if (err) {
	gretl_errmsg_set("Failed to recreate gridplot");
	g_ptr_array_unref(arr);
    } else {
	multiplot = arr;
    }

    return err;
}

/* Revise an existing gridplot buffer or command file,
   presumably obtained via "end gpbuild" with the --output
   or --outbuf option or perhaps via the "standalone"
   usage of gridplot. This may just be a matter of selecting
   an output format, or it may involve changes to options
   such as font size or layout.
*/

int gretl_multiplot_revise (gretlopt opt)
{
    const char *buf = NULL;
    gchar *gbuf = NULL;
    gretlopt myopt = opt;
    int maxp, np = 0;
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
    if (!err) {
	/* extract the subplot specs */
	err = disassemble_multiplot(buf, np);
    }
    if (err) {
        return err;
    }

#if GRID_DEBUG
    /* what sort of output is wanted? */
    const char *argname = get_optval_string(GRIDPLOT, OPT_U);
    if (argname != NULL) {
        fprintf(stderr, "gretl_multiplot_revise: output '%s'\n",
                argname);
    }
#endif

    maxp = np;

    /* let the current @opt override previous choices */
    myopt &= ~(OPT_i | OPT_I | OPT_U);
    if (myopt) {
        err = update_multiplot_sizes(myopt);
    }
    if (myopt & OPT_L) {
        err = maybe_set_mp_layout(GRIDPLOT, &np);
    }
    if (!err && mp_layout == NULL) {
	err = multiplot_set_grid(np);
    }

#if GRID_DEBUG > 1
    fprintf(stderr, "*** here's the prior input (filter = %d) ***\n", filter);
    fputs(buf, stderr);
#endif

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
	output_multiplot_script(GRIDPLOT, NULL, np, maxp);
    }

    if (gbuf != NULL) {
	g_free(gbuf);
    }

    gretl_multiplot_destroy();

    return err;
}

static int retrieve_plots_array (const char ***pS, int *pnp)
{
    const char *argname;
    gretl_array *a = NULL;
    int msg_set = 0;
    int err = 0;

    /* get the @name in --strings=name */
    argname = get_optval_string(GRIDPLOT, OPT_S);
    a = get_array_by_name(argname);
    if (a == NULL) {
	err = E_INVARG;
    } else {
	const char **S;
	int i;

	S = (const char **) gretl_array_get_strings(a, pnp);
	if (S == NULL) {
	    err = E_INVARG;
	} else {
	    /* check for embedded multiplots */
	    for (i=0; i<*pnp; i++) {
		if (strstr(S[i], "set multiplot")) {
		    err = invalid_mp_error(GRIDPLOT);
		    msg_set = 1;
		    break;
		}
	    }
	    if (!err) {
		*pS = S;
	    }
	}
    }

    if (err && !msg_set) {
	gretl_errmsg_set("Didn't get a valid array of plot strings");
    }

    return err;
}

/* This supports "standalone" usage of gridplot to process a
   "free form" array of individual plot-specification strings.
*/

int gretl_multiplot_from_array (gretlopt opt)
{
    const char **S = NULL;
    int maxp, np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return E_DATA;
    }

    err = retrieve_plots_array(&S, &np);
    maxp = np;

    if (!err) {
	/* pick up any options */
	gretlopt myopt = opt;

	set_multiplot_defaults();
	myopt &= ~OPT_S;
	if (myopt) {
	    err = set_multiplot_sizes(GRIDPLOT, myopt);
	}
	if (myopt & OPT_L) {
	    err = maybe_set_mp_layout(GRIDPLOT, &np);
	}
	if (!err && mp_layout == NULL) {
	    err = multiplot_set_grid(np);
	}
    }

    if (!err) {
	set_special_plot_size(mp_width, mp_height);
	set_special_font_size(mp_fontsize);
        err = output_multiplot_script(GRIDPLOT, S, np, maxp);
    }

    return err;
}

int check_multiplot_options (int ci, gretlopt opt)
{
    int err = 0;

    if (ci == GRIDPLOT) {
	/* no more than one of --input, --inbuf, --strings */
	err = incompatible_options(opt, OPT_I | OPT_i | OPT_S);
    }

    if (!err) {
	/* can't have both --output and --outbuf */
	err = incompatible_options(opt, OPT_U | OPT_b);
    }

    return err;
}
