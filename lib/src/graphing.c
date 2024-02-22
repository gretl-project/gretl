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

/* graphing.c for gretl */

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "libset.h"
#include "matrix_extra.h"
#include "forecast.h"
#include "plotspec.h"
#include "usermat.h"
#include "gretl_panel.h"
#include "missing_private.h"
#include "gretl_string_table.h"
#include "gretl_gridplot.h"
#include "uservar.h"
#include "gretl_midas.h"
#include "boxplots.h"
#include "mapinfo.h"
#include "gretl_func.h"
#include "plot_priv.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <unistd.h>
#include <errno.h>

#define GP_DEBUG 0

#ifdef WIN32
# include <windows.h>
#else
# include <signal.h>
# if HAVE_SYS_WAIT_H
#  include <sys/wait.h>
# endif
# ifndef WEXITSTATUS
#  define WEXITSTATUS(stat_val) ((unsigned)(stat_val) >> 8)
# endif
# ifndef WIFEXITED
#  define WIFEXITED(stat_val) (((stat_val) & 255) == 0)
# endif
#endif /* ! _WIN32 */

/* length of buffer for "set term ..." */
#define TERMLEN 256

static char gnuplot_path[MAXLEN];
static int gp_small_font_size;
static double default_png_scale = 1.0;
static int xwide = 0;

static char ad_hoc_font[64];
static char plot_buffer_name[32];
static char plot_buffer_idx[8];

enum {
    W_POINTS,
    W_LINES,
    W_IMPULSES,
    W_LP,
    W_BOXES,
    W_STEPS
};

#define MAX_LETTERBOX_LINES 8

#define ts_plot(g) ((g)->flags & GPT_TS)

#if GP_DEBUG
static void print_gnuplot_flags (int flags, int revised);
#endif

struct plot_type_info {
    PlotType ptype;
    const char *pstr;
};

struct plot_type_info ptinfo[] = {
    { PLOT_REGULAR,        NULL },
    { PLOT_CORRELOGRAM,    "correlogram" },
    { PLOT_CUSUM,          "CUSUM test" },
    { PLOT_FORECAST,       "forecasts with 95 pc conf. interval" },
    { PLOT_FREQ_SIMPLE,    "frequency plot (simple)" },
    { PLOT_FREQ_NORMAL,    "frequency plot (against normal)" },
    { PLOT_FREQ_GAMMA,     "frequency plot (against gamma)" },
    { PLOT_FREQ_DISCRETE,  "frequency plot (discrete)" },
    { PLOT_GARCH,          "GARCH residual plot" },
    { PLOT_HURST,          "rescaled range plot" },
    { PLOT_IRFBOOT,        "impulse response plot with quantiles" },
    { PLOT_KERNEL,         "kernel density plot" },
    { PLOT_LEVERAGE,       "leverage/influence plot" },
    { PLOT_MULTI_BASIC,    "multiple small plots" },
    { PLOT_PERIODOGRAM,    "periodogram" },
    { PLOT_RANGE_MEAN,     "range-mean plot" },
    { PLOT_H_TEST,         "sampling distribution" },
    { PLOT_PROB_DIST,      "probability distribution" },
    { PLOT_TRI_GRAPH,      "TRAMO / X12A tri-graph" },
    { PLOT_ROOTS,          "roots plot" },
    { PLOT_ELLIPSE,        "confidence ellipse plot" },
    { PLOT_MULTI_IRF,      "multiple impulse responses" },
    { PLOT_PANEL,          "multiple panel plots" },
    { PLOT_BI_GRAPH,       "double time-series plot" },
    { PLOT_MANY_TS,        "multiple timeseries" },
    { PLOT_RQ_TAU,         "tau sequence plot" },
    { PLOT_FACTORIZED,     "factorized scatter" },
    { PLOT_BOXPLOTS,       "boxplots" },
    { PLOT_CURVE,          "curve" },
    { PLOT_QQ,             "QQ plot" },
    { PLOT_USER,           "user-defined plot" },
    { PLOT_XCORRELOGRAM,   "cross-correlogram" },
    { PLOT_BAR,            "bars" },
    { PLOT_STACKED_BAR,    "stacked-bars" },
    { PLOT_3D,             "3-D plot" },
    { PLOT_BAND,           "band plot" },
    { PLOT_HEATMAP,        "heatmap" },
    { PLOT_GEOMAP,         "geoplot" },
    { PLOT_GRIDPLOT,       "user-multi" },
    { PLOT_TYPE_MAX,       NULL }
};

static void get_multiplot_layout (int n, int tseries,
				  int *rows, int *cols);
static char *gretl_emf_term_line (char *term_line,
				  PlotType ptype,
				  GptFlags flags);

#ifdef WIN32
static void win32_forwardize (gchar *fname);
#endif

#ifndef WIN32

#define SPAWN_DEBUG 0

/**
 * gnuplot_test_command:
 * @cmd: gnuplot command string.
 *
 * See if the installed version of gnuplot will accept a given
 * command.
 *
 * Returns: 0 if gnuplot successfully handles the given command,
 * 1 on error.
 */

int gnuplot_test_command (const char *cmd)
{
    int ok, ret = 1;
    int child_pid = 0, sinp = 0, serr = 0;
    GError *error = NULL;
    gchar *argv[] = {
	NULL,
	NULL
    };

    if (*gnuplot_path == '\0') {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    argv[0] = gnuplot_path;

    ok = g_spawn_async_with_pipes (NULL,
				   argv,
				   NULL,
				   G_SPAWN_SEARCH_PATH |
				   G_SPAWN_STDOUT_TO_DEV_NULL |
				   G_SPAWN_DO_NOT_REAP_CHILD,
				   NULL,
				   NULL,
				   &child_pid,
				   &sinp,
				   NULL,
				   &serr,
				   &error);

# if SPAWN_DEBUG
    fprintf(stderr, "Testing gnuplot command '%s'\n", cmd);
    fprintf(stderr, "ok=%d, child_pid=%d, sinp=%d\n",
	    ok, child_pid, sinp);
# endif

    if (ok) {
	char errbuf[128];
	int test, status;
	int errbytes = 0;

	errbytes += write(sinp, cmd, strlen(cmd));
	errbytes += write(sinp, "\n", 1);
	close(sinp);
	test = waitpid(child_pid, &status, 0);
# if SPAWN_DEBUG
	fprintf(stderr, "waitpid returned %d, WIFEXITED %d, "
		"WEXITSTATUS %d\n", test, WIFEXITED(status),
		WEXITSTATUS(status));
# endif
	if (test == child_pid && WIFEXITED(status)) {
	    ret = WEXITSTATUS(status);
	}
	errbytes = read(serr, errbuf, sizeof errbuf - 1);
	if (errbytes > 0) {
	    errbuf[errbytes] = '\0';
	    if (strstr(errbuf, "not find/open font")) {
# if SPAWN_DEBUG
		fprintf(stderr, "%s\n", errbuf);
# endif
		if (strstr(cmd, "font") != NULL) {
		    ret = 1;
		}
	    }
	}
	close(serr);
    } else {
	fprintf(stderr, "error: '%s'\n", error->message);
	g_error_free(error);
    }

# if SPAWN_DEBUG
    fprintf(stderr, "gnuplot test: ret = %d\n", ret);
# endif

    return ret;
}

#endif /* not WIN32 */

/* Apparatus for getting the gnuplot version, either
   in numerical form or as a string
*/

static char gpver_string[16];
static double gpver_number;

static void set_gpver_number (const char *s)
{
    const char *fmt = "%d.%d.%d";
    int maj, min, plev;

    if (sscanf(s, fmt, &maj, &min, &plev) == 3) {
	gpver_number = maj + min / 10.0 + plev / 100.0;
    }
}

static int get_gp_version_info (void)
{
    gchar *qname = NULL;
    gchar *fname = NULL;
    gchar *buf = NULL;
    FILE *fp;
    int err = 0;

    if (*gnuplot_path == '\0') {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    qname = gretl_make_dotpath("gpver_query");
    fname = gretl_make_dotpath("gpver.txt");
#ifdef WIN32
    win32_forwardize(fname);
#endif

    fp = gretl_fopen(qname, "w");
    if (fp == NULL) {
	err = E_FOPEN;
    } else {
	/* write a little query script */
	fputs("gpver = sprintf(\"%.1f.%s\", GPVAL_VERSION, GPVAL_PATCHLEVEL)\n", fp);
	fprintf(fp, "set print \"%s\"\n", fname);
	fputs("print gpver\n", fp);
	fclose(fp);
    }

    if (!err) {
	/* send script to gnuplot for execution */
	buf = g_strdup_printf("\"%s\" \"%s\"", gnuplot_path, qname);
	err = gretl_spawn(buf);
    }

    if (!err) {
	/* read version info from response */
	fp = gretl_fopen(fname, "r");
	if (fp != NULL) {
	    if (fgets(gpver_string, sizeof gpver_string, fp)) {
		g_strstrip(gpver_string);
		set_gpver_number(gpver_string);
	    }
	    fclose(fp);
	}
    }

    gretl_remove(qname);
    gretl_remove(fname);
    g_free(qname);
    g_free(fname);
    g_free(buf);

    return err;
}

double gnuplot_version (void)
{
    if (gpver_number == 0) {
	get_gp_version_info();
    }

    return gpver_number;
}

char *gnuplot_version_string (void)
{
    if (gpver_number == 0) {
	get_gp_version_info();
    }

    return gpver_number == 0 ? "unknown" : gpver_string;
}

/* end gnuplot version apparatus */

static int gp_list_pos (const char *s, const int *list,
			const DATASET *dset)
{
    int k;

    if (integer_string(s)) {
	k = atoi(s);
    } else {
	k = current_series_index(dset, s);
    }

    return in_gretl_list(list, k);
}

static int plot_ci = GNUPLOT;

void set_effective_plot_ci (int ci)
{
    plot_ci = ci;
}

int get_effective_plot_ci (void)
{
    return plot_ci;
}

/* When we get from the user something like

   --with-lines=foo,bar

   this indicates that the "with lines" format should be
   applied to selected y-axis variables, not all.
*/

static int gp_set_non_point_info (gnuplot_info *gi,
				  const int *list,
				  const DATASET *dset,
				  gretlopt opt)
{
    const char *s = get_optval_string(plot_ci, opt);
    int withval = W_POINTS;
    int i, imax = gi->withlist[0];

    if (opt == OPT_O) {
	withval = W_LINES;
    } else if (opt == OPT_M) {
	withval = W_IMPULSES;
    } else if (opt == OPT_P) {
	withval = W_LP;
    } else if (opt == OPT_B) {
	withval = W_BOXES;
    } else if (opt & OPT_Q) {
	withval = W_STEPS;
    }

    if (s == NULL) {
	/* spec applies to all members of list */
	for (i=1; i<=imax; i++) {
	    if (gi->withlist[i] == W_POINTS) {
		gi->withlist[i] = withval;
	    }
	}
    } else if (strchr(s, ',') != NULL) {
	/* spec has multiple components */
	gchar **strs = g_strsplit(s, ",", 0);
	int j;

	for (j=0; strs[j]!=NULL; j++) {
	    i = gp_list_pos(strs[j], list, dset);
	    if (i > 0 && i <= imax) {
		gi->withlist[i] = withval;
	    }
	}
	g_strfreev(strs);
    } else {
	/* just one component */
	i = gp_list_pos(s, list, dset);
	if (i > 0 && i <= imax) {
	    gi->withlist[i] = withval;
	}
    }

    return 0;
}

static int plain_lines_spec (gretlopt opt)
{
    if ((opt & OPT_O) && !(opt & (OPT_M | OPT_B | OPT_P | OPT_Q))) {
	return get_optval_string(plot_ci, OPT_O) == NULL;
    } else {
	return 0;
    }
}

static int plain_impulses_spec (gretlopt opt)
{
    if ((opt & OPT_M) && !(opt & (OPT_O | OPT_B | OPT_P | OPT_Q))) {
	return get_optval_string(plot_ci, OPT_M) == NULL;
    } else {
	return 0;
    }
}

static int plain_steps_spec (gretlopt opt)
{
    if ((opt & OPT_Q) && !(opt & (OPT_O | OPT_M | OPT_B | OPT_P))) {
	return get_optval_string(plot_ci, OPT_Q) == NULL;
    } else {
	return 0;
    }
}

static int get_fit_type (gnuplot_info *gi)
{
    const char *ftype = get_optval_string(plot_ci, OPT_F);
    int err = 0;

    if (ftype == NULL || *ftype == '\0') {
	err = E_DATA;
    } else if (!strcmp(ftype, "none")) {
	gi->flags |= GPT_FIT_OMIT;
    } else if (!strcmp(ftype, "linear")) {
	gi->fit = PLOT_FIT_OLS;
    } else if (!strcmp(ftype, "quadratic")) {
	gi->fit = PLOT_FIT_QUADRATIC;
    } else if (!strcmp(ftype, "cubic")) {
	gi->fit = PLOT_FIT_CUBIC;
    } else if (!strcmp(ftype, "inverse")) {
	gi->fit = PLOT_FIT_INVERSE;
    } else if (!strcmp(ftype, "loess")) {
	gi->fit = PLOT_FIT_LOESS;
    } else if (!strcmp(ftype, "semilog")) {
	gi->fit = PLOT_FIT_LOGLIN;
    } else if (!strcmp(ftype, "linlog")) {
	gi->fit = PLOT_FIT_LINLOG;
    } else {
	err = invalid_field_error(ftype);
    }

    return err;
}

static void maybe_record_font_choice (gretlopt opt)
{
    const char *s = get_optval_string(plot_ci, opt);

    if (s != NULL) {
	ad_hoc_font[0] = '\0';
	strcat(ad_hoc_font, s);
	gretl_charsub(ad_hoc_font, ',', ' ');
    }
}

/* We come here in response to the --buffer option in plot
   commands, and if applicable we set the static variables
   @plot_buffer_name and @plot_buffer_idx.
*/

static int set_plot_buffer_name (const char *bname)
{
    int err = 0;

    if (bname == NULL || *bname == '\0') {
	*plot_buffer_name = '\0';
        *plot_buffer_idx = '\0';
    } else if (is_strings_array_element(bname, plot_buffer_name, plot_buffer_idx)) {
	; /* handled */
    } else {
	*plot_buffer_idx = '\0';
	err = check_stringvar_name(bname, 1, NULL);
	if (!err) {
	    strcpy(plot_buffer_name, bname);
	}
    }

    return err;
}

int plot_output_to_buffer (void)
{
    return *plot_buffer_name != '\0' ||
	gretl_multiplot_collecting();
}

static int make_plot_commands_buffer (const char *fname)
{
    gchar *contents = NULL;
    GError *gerr = NULL;
    gsize size = 0;
    int err = 0;

    g_file_get_contents(fname, &contents, &size, &gerr);

    if (gerr != NULL) {
	gretl_errmsg_set(gerr->message);
	g_error_free(gerr);
	err = E_FOPEN;
    } else if (gretl_multiplot_collecting()) {
	err = gretl_multiplot_add_plot(contents);
    } else if (*plot_buffer_idx != '\0') {
	gretl_array *a = get_strings_array_by_name(plot_buffer_name);
        int i = generate_int(plot_buffer_idx, NULL, &err);

        if (a != NULL && !err) {
            err = gretl_array_set_string(a, i-1, contents, 1);
        }
    } else {
	char *buf = gretl_strdup(contents);

	err = user_var_add_or_replace(plot_buffer_name,
				      GRETL_TYPE_STRING,
				      buf);
    }

    g_free(contents);
    gretl_remove(fname);

    return err;
}

static int get_gp_flags (gnuplot_info *gi, gretlopt opt,
			 const int *list, const DATASET *dset)
{
    int n_yvars = list[0] - 1;
    int err = 0;

    gi->flags = 0;

    if (opt & (OPT_N | OPT_a)) {
	/* --band or --bands */
	if (opt & OPT_T) {
	    /* --time-series */
	    gi->flags |= (GPT_TS | GPT_IDX);
	    /* there's no xvar in @list */
	    n_yvars++;
	}
	gi->flags |= GPT_FIT_OMIT;
	gi->band = 1;
	goto linespec;
    }

    if (opt & OPT_W) {
	/* --font=<fontspec> */
	maybe_record_font_choice(OPT_W);
    }

    if (opt & OPT_L) {
	/* log y axis */
	const char *sbase = get_optval_string(GNUPLOT, OPT_L);

	gi->flags |= GPT_LOGY;
	gi->ybase = 10;
	if (sbase != NULL) {
	    gi->ybase = atof(sbase);
	    if (gi->ybase <= 0) {
		gi->ybase = 10;
	    }
	}
    }

    if (opt & OPT_S) {
	/* the old --suppress-fitted option may still be used
	   internally, for some plot types */
	gi->flags |= GPT_FIT_OMIT;
    }

    if (opt & OPT_R) {
	/* internal option for residual plot */
	gi->flags |= GPT_RESIDS;
    } else if (opt & OPT_A) {
	/* internal option for fitted-actual plot */
	gi->flags |= GPT_FA;
    }

    if (opt & OPT_Z) {
	/* --dummy */
	gi->flags |= GPT_DUMMY;
    } else if (opt & OPT_C) {
	/* --control */
	gi->flags |= GPT_XYZ;
    } else {
	if (opt & OPT_T) {
	    /* --time-series */
	    gi->flags |= GPT_IDX;
	    /* there's no xvar in @list */
	    n_yvars++;
	}
    }

 linespec:

    if (plain_lines_spec(opt)) {
	/* just using lines */
	gi->flags |= GPT_LINES;
    } else if (plain_impulses_spec(opt)) {
	/* just using impulses */
	gi->flags |= GPT_IMPULSES;
    } else if (plain_steps_spec(opt)) {
	/* just using steps */
	gi->flags |= GPT_STEPS;
    } else if (opt & (OPT_M | OPT_O | OPT_P | OPT_B)) {
	/* for handling per-variable "plot with" options */
	gi->withlist = gretl_list_new(n_yvars);
    }

    if (gi->withlist != NULL) {
	if (opt & OPT_M) {
	    /* --with-impulses */
	    gp_set_non_point_info(gi, list, dset, OPT_M);
	}
	if (opt & OPT_O) {
	    /* --with-lines */
	    gp_set_non_point_info(gi, list, dset, OPT_O);
	}
	if (opt & OPT_P) {
	    /* --with-lp */
	    gp_set_non_point_info(gi, list, dset, OPT_P);
	}
	if (opt & OPT_B) {
	    /* --with-boxes */
	    gp_set_non_point_info(gi, list, dset, OPT_B);
	}
    }

    if (opt & OPT_G) {
	/* internal option, saving as icon */
	gi->flags |= GPT_ICON;
    }

    gi->fit = PLOT_FIT_NONE;

    if (!(gi->flags & GPT_FIT_OMIT) && n_yvars == 1) {
	if (opt & OPT_F) {
	    /* the --fit=fitspec option */
	    err = get_fit_type(gi);
	}
    }

    if (xwide) {
	/* access file-scope global */
	gi->flags |= GPT_XW;
	xwide = 0;
    }

#if GP_DEBUG
    if (gi->flags) {
	print_gnuplot_flags(gi->flags, 0);
    }
#endif

    return err;
}

void write_gp_dataval (double x, FILE *fp, int final)
{
    if (final) {
	if (na(x)) {
	    fprintf(fp, "%s\n", GPNA);
	} else {
	    fprintf(fp, "%.10g\n", x);
	}
    } else {
	if (na(x)) {
	    fprintf(fp, "%s ", GPNA);
	} else {
	    fprintf(fp, "%.10g ", x);
	}
    }
}

static void printvars (FILE *fp, int t,
		       const int *list,
		       const DATASET *dset,
		       gnuplot_info *gi,
		       const char *label,
		       double offset)
{
    const double *x = (gi != NULL)? gi->x : NULL;
    double xt;
    int i;

    if (x != NULL) {
	xt = x[t] + offset;
	if (gi->flags & GPT_TIMEFMT) {
	    fprintf(fp, "%.0f ", xt);
	} else {
	    fprintf(fp, "%.10g ", xt);
	}
    }

    for (i=1; i<=list[0]; i++) {
	xt = dset->Z[list[i]][t];
	if (!na(xt) && x == NULL && i == 1) {
	    /* the x variable */
	    xt += offset;
	}
	write_gp_dataval(xt, fp, 0);
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);
}

static int factor_check (gnuplot_info *gi, const DATASET *dset)
{
    int err = 0;
    int v3 = 0;

    if (gi->list[0] != 3) {
	err = E_DATA;
    } else {
	v3 = gi->list[3];
	if (!accept_as_discrete(dset, v3, 0)) {
	    err = E_DATA;
	}
    }

    if (err) {
	gretl_errmsg_set(_("You must supply three variables, the last of "
			   "which is discrete"));
    } else {
	const double *d = dset->Z[v3] + gi->t1;
	int T = gi->t2 - gi->t1 + 1;

	gi->dvals = gretl_matrix_values(d, T, OPT_S, &err);
    }

    return err;
}

#ifndef WIN32

int gnuplot_has_wxt (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set term wxt");
    }

    return !err;
}

static int gnuplot_has_x11 (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set term x11");
    }

    return !err;
}

static int gnuplot_has_qt (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set term qt");
    }

    return !err;
}

static int gnuplot_has_tikz (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set term tikz");
    }

    return !err;
}

#else

int gnuplot_has_wxt (void)
{
    /* There's no WxWidgets support in the current
       Windows build of gnuplot 5
    */
    return 0;
}

static int gnuplot_has_tikz (void)
{
    /* There's no Lua/TikZ support in the current
       Windows build of gnuplot 5
    */
    return 0;
}

#endif /* !WIN32 or WIN32 */

static gretlRGB user_color[N_GP_LINETYPES] = {
    0xff0000,
    0x0000ff,
    0x00cc00,
    0xbf25b2,
    0x8faab3,
    0xffa500,
    0xe51e10,
    0x000000
};

/* apparatus for handling plot "extra" colors */

static gretlRGB extra_color[2] = {
    0x5f6b84,
    0xdddddd
};

static gretlRGB user_extra_color[2] = {
    0x5f6b84,
    0xdddddd
};

gretlRGB get_boxcolor (void)
{
    return user_extra_color[0];
}

void set_boxcolor (gretlRGB color) {
    user_extra_color[0] = color;
}

gretlRGB get_shadecolor (void)
{
    return user_extra_color[1];
}

void set_shadecolor (gretlRGB color) {
    user_extra_color[1] = color;
}

void print_rgb_hash (char *s, gretlRGB color)
{
    sprintf(s, "#%06X", color);
}

gretlRGB gretl_rgb_get (const char *s)
{
    gretlRGB x = 0;

    if (sscanf(s, "#%x", &x) != 1) {
	x = 0;
    }

    return x;
}

void print_palette_string (char *s)
{
    sprintf(s, "x%06X x%06X", user_extra_color[0],
	    user_extra_color[0]);
}

gretlRGB get_graph_color (int i)
{
    return (i >= 0 && i < N_GP_LINETYPES)? user_color[i] : 0;
}

void set_graph_color_from_string (int i, const char *s)
{
    int err = 0;

    if (i >= 0 && i < 2) {
	gretlRGB x;

	if (sscanf(s + 1, "%06x", &x) == 1) {
	    user_extra_color[i] = x;
	} else {
	    err = 1;
	}
    } else {
	err = 1;
    }

    if (err) {
	fprintf(stderr, "Error in set_graph_palette_from_string(%d, '%s')\n",
		i, s);
    }
}

void graph_palette_reset (int i)
{
    if (i >= 0 && i < 2) {
	user_extra_color[i] = extra_color[i];
    }
}

/* Given a string @s such as "Sans 8" or "Bodoni MT 12", write
   the name part into @name and the point-size part into @psz.
   Return 2 if we got both a name and a size, 1 if we just
   got a name, 0 if we got nothing.
*/

int split_graph_fontspec (const char *s, char *name, int *psz)
{
    int i, k = 0, n = strlen(s);
    int nf = 0;

    for (i=n-1; i>0; i--) {
	if (isdigit(s[i])) k++;
	else break;
    }

    if (k > 0) {
	/* got a size */
	char ptstr[8];

	*ptstr = *name = '\0';
	strncat(ptstr, s + n - k, k);
	*psz = atoi(ptstr);
	strncat(name, s, n - k - 1);
	nf = 2;
    } else if (*s != '\0') {
	nf = 1;
	strcpy(name, s);
    }

    return nf;
}

char *gretl_png_font_string (void)
{
    const char *s = gretl_png_font();
    char fstr[256];
    char name[128];
    int nf, ptsize = 0;

    fstr[0] = '\0';
    nf = split_graph_fontspec(s, name, &ptsize);
    if (nf == 2) {
	sprintf(fstr, " font \"%s,%d\"", name, ptsize);
    } else if (nf == 1) {
	sprintf(fstr, " font \"%s\"", name);
    }

    return gretl_strdup(fstr);
}

static void maybe_set_small_font (int nplots)
{
    gp_small_font_size = (nplots > 4)? 6 : 0;
}

/* apparatus for plots at custom sizes (e.g. maps) */

static float special_width;
static float special_height;
static int special_fontsize;

void set_special_plot_size (float width, float height)
{
    special_width = width;
    special_height = height;
}

void set_special_font_size (int fsize)
{
    special_fontsize = fsize;
}

static void clear_special_plot_size (void)
{
    special_width = special_height = 0;
}

static void clear_special_font_size (void)
{
    special_fontsize = 0;
}

static int special_plot_size_is_set (void)
{
    return special_width > 0 && special_height > 0;
}

static int special_font_size_is_set (void)
{
    return special_fontsize > 0;
}

/* end special size apparatus */

static void write_png_font_string (char *fstr,
				   char *ad_hoc_fontspec,
				   PlotType ptype,
				   const char *grfont,
				   double scale)
{
    int adhoc = 0;

    if (grfont == NULL) {
	if (ad_hoc_font[0] != '\0') {
	    adhoc = 1;
	    grfont = ad_hoc_font;
	} else {
	    grfont = gretl_png_font();
	}
    }

    if (*grfont == '\0') {
	grfont = getenv("GRETL_PNG_GRAPH_FONT");
    }

    if (*grfont == '\0') {
	*fstr = '\0';
	return;
    } else {
	char fname[128];
	int nf, fsize = 0;

	nf = split_graph_fontspec(grfont, fname, &fsize);

	if (special_font_size_is_set()) {
	    fsize = special_fontsize;
	} else if (nf == 2) {
	    if (maybe_big_multiplot(ptype) && gp_small_font_size > 0) {
		fsize = gp_small_font_size;
	    }
	    if (scale > 1.0) {
		fsize = round(scale * fsize);
	    }
	}
	if (fsize > 0 && fsize < 100) {
	    sprintf(fstr, " font \"%s,%d\"", fname, fsize);
	} else if (nf == 1) {
	    sprintf(fstr, " font \"%s\"", fname);
	}
	if (adhoc) {
	    strcpy(ad_hoc_fontspec, grfont);
	}
	/* ensure these settings don't outstay their welcome */
	ad_hoc_font[0] = '\0';
	clear_special_font_size();
    }
}

/* for gnuplot pdfcairo, epscairo output */

static gchar *write_other_font_string (int stdsize)
{
    gchar *fstr = NULL;

    if (ad_hoc_font[0] != '\0') {
	char fname[128];
	int nf, fsize = 0;

	nf = split_graph_fontspec(ad_hoc_font, fname, &fsize);
	if (nf == 2) {
	    fstr = g_strdup_printf("%s,%d", fname, fsize);
	} else if (nf == 1) {
	    fstr = g_strdup_printf("%s,%d", fname, stdsize);
	}
	ad_hoc_font[0] = '\0';
    } else if (special_font_size_is_set()) {
	fstr = g_strdup_printf("sans,%d", special_fontsize);
    } else {
	fstr = g_strdup_printf("sans,%d", stdsize);
    }

    return fstr;
}

/* gnuplot styles apparatus */

static char gp_style[32]; /* basename of style file */
static gchar *alt_sty;    /* content of style file */

static const char *classic_sty =
    "# gpstyle classic\n"
    "set linetype 1 pt 1 lc rgb \"#FF0000\"\n"  /* red */
    "set linetype 2 pt 2 lc rgb \"#0000FF\"\n"  /* blue */
    "set linetype 3 pt 3 lc rgb \"#00CC00\"\n"  /* non-standard green */
    "set linetype 4 pt 4 lc rgb \"#BF25B2\"\n"  /* purple */
    "set linetype 5 pt 5 lc rgb \"#8FAAB3\"\n"  /* gray-blue */
    "set linetype 6 pt 6 lc rgb \"#FFA500\"\n"  /* yellow-orange */
    "set linetype 7 pt 7 lc rgb \"#E51E10\"\n"  /* unnamed red */
    "set linetype 8 pt 8 lc rgb \"#000000\"\n"; /* black */

/* Read 8 line colors out of @s and transcribe to the
   @user_colors array. Fail if we can't get all 8.
*/

static int transcribe_style (const char *s)
{
    const char *p = s;
    gretlRGB rgb;
    int i, got;

    for (i=0; i<N_GP_LINETYPES; i++) {
	p = strstr(p, " lc rgb ");
	got = 0;
	if (p != NULL) {
	    p += 8;
	    p += strspn(p, " ");
	    if (sscanf(p+1, "#%x", &rgb) == 1) {
		user_color[i] = rgb;
		got = 1;
		p += 8;
	    }
	}
	if (!got) {
	    break;
	}
    }

    return i < N_GP_LINETYPES ? E_DATA : 0;
}

/* Given a style-name (provisionally copied to @gp_style),
   try to find its file and check it for conformance.
*/

static int try_set_alt_sty (void)
{
    GError *gerr = NULL;
    gsize len = 0;
    gchar *try = NULL;
    gchar *fname;
    int err = 0;

    fname = g_build_filename(gretl_home(), "data",
			     "gnuplot", gp_style, NULL);

    if (g_file_get_contents(fname, &try, &len, &gerr)) {
	err = transcribe_style(try);
	if (err) {
	    /* failed to conform to spec */
	    fprintf(stderr, "%s failed spec check\n", gp_style);
	    set_plotstyle("classic");
	} else {
	    /* OK, put the style in place */
	    g_free(alt_sty);
	    alt_sty = try;
	}
    } else {
	/* couldn't find the file */
	if (gerr != NULL) {
	    fprintf(stderr, "%s\n", gerr->message);
	    g_error_free(gerr);
	}
	err = E_FOPEN;
    }

    g_free(fname);

    return err;
}

/* callback from libset.c, in response to "set plot_style <name>" */

int set_plotstyle (const char *style)
{
    int to_classic = 0;

    if (!strcmp(style, "classic") || !strcmp(style, "default")) {
	/* "default" is just for backward compat */
	to_classic = 1;
    }

    if (!strcmp(style, gp_style)) {
	return 0; /* no-op */
    } else if (to_classic && *gp_style == '\0') {
	return 0; /* no-op */
    } else if (to_classic) {
	/* replace alt with classic */
	g_free(alt_sty);
	alt_sty = NULL;
	gp_style[0] = '\0';
	transcribe_style(classic_sty);
	return 0;
    } else {
	/* try replacing current with what's requested */
	sprintf(gp_style, "%s.gpsty", style);
	return try_set_alt_sty();
    }
}

const char *get_plotstyle (void)
{
    static char pstyle[32];
    char *p;

    strcpy(pstyle, gp_style);
    p = strrchr(pstyle, '.');
    if (p != NULL) {
	*p = '\0';
    }

    return pstyle;
}

/* Write the content of either the default, or an alternative,
   gnuplot style into @fp. The @offset argument allows for
   skipping one or more leading linetype definitions. The
   @linetypes argument allows for skipping the linetypes
   altogether.
*/

static void inject_gp_style (int offset, int linetypes, FILE *fp)
{
    const char *sty = alt_sty != NULL ? alt_sty : classic_sty;
    const char *sub = NULL;

    if (linetypes == 0) {
        /* just print what follows the linetypes, if anything */
        sub = strstr(sty, "set bord");
        if (sub != NULL) {
            fputs(sub, fp);
        }
        return;
    } else if (offset > 0) {
        /* start at a specified linetype */
	char targ[32];

	sprintf(targ, "set linetype %d", offset + 1);
	sub = strstr(sty, targ);
    }

    fputs(sub != NULL ? sub : sty, fp);
}

void write_plot_line_styles (int ptype, FILE *fp)
{
    char cstr[12];
    int i;

    if (ptype == PLOT_3D) {
	for (i=0; i<2; i++) {
	    print_rgb_hash(cstr, user_color[i]);
	    fprintf(fp, "set linetype %d lc rgb \"%s\"\n", i+1, cstr);
	}
    } else if (ptype == PLOT_BOXPLOTS) {
	for (i=0; i<2; i++) {
	    print_rgb_hash(cstr, user_color[i+1]);
	    fprintf(fp, "set linetype %d lc rgb \"%s\"\n", i+1, cstr);
	}
	inject_gp_style(0, 0, fp);
    } else if (frequency_plot_code(ptype)) {
	print_rgb_hash(cstr, get_boxcolor());
	fprintf(fp, "set linetype 1 lc rgb \"%s\"\n", cstr);
	fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
        inject_gp_style(0, 0, fp);
    } else if (ptype == PLOT_RQ_TAU) {
	fputs("set linetype 1 lc rgb \"#000000\"\n", fp);
	fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
	fputs("set linetype 3 lc rgb \"#0000FF\"\n", fp);
	fputs("set linetype 4 lc rgb \"#0000FF\"\n", fp);
	fputs("set linetype 5 lc rgb \"#0000FF\"\n", fp);
	fputs("set linetype 6 lc rgb \"#FFA500\"\n", fp);
	fputs("set linetype 7 lc rgb \"#E51E10\"\n", fp);
	fputs("set linetype 8 lc rgb \"#000000\"\n", fp);
        inject_gp_style(0, 0, fp);
    } else if (ptype == PLOT_HEATMAP || ptype == PLOT_GEOMAP) {
	; /* these are handled specially */
    } else {
	/* the primary default case */
	inject_gp_style(0, 1, fp);
    }
}

/* end gnuplot styles apparatus */

#ifdef WIN32

static void win32_forwardize (gchar *fname)
{
    gchar *p;

    while ((p = strchr(fname, '\\')) != NULL) {
	*p = '/';
	p++;
    }
}

#endif

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds (gnuplot >= 4.4.0).
*/

int write_plot_bounding_box_request (FILE *fp)
{
    gchar *fname = gretl_make_dotpath("gretltmp.png.bounds");

#ifdef WIN32
    win32_forwardize(fname);
#endif

    fprintf(fp, "set print \"%s\"\n", fname);
    g_free(fname);

    fputs("print \"pixel_bounds: \", GPVAL_TERM_XMIN, GPVAL_TERM_XMAX, "
	  "GPVAL_TERM_YMIN, GPVAL_TERM_YMAX\n", fp);
    fputs("print \"data_bounds: \", GPVAL_X_MIN, GPVAL_X_MAX, "
	  "GPVAL_Y_MIN, GPVAL_Y_MAX\n", fp);
    fputs("print \"term_size: \", GPVAL_TERM_XSIZE, "
	  "GPVAL_TERM_YSIZE, GPVAL_TERM_SCALE\n", fp);

    return 0;
}

static int do_plot_bounding_box (void)
{
    FILE *fp = gretl_fopen(gretl_plotfile(), "ab");
    int err = 0;

    if (fp != NULL) {
	err = write_plot_bounding_box_request(fp);
	fclose(fp);
    } else {
	err = E_FOPEN;
    }

    return err;
}

static void maybe_set_eps_pdf_dims (char *s, PlotType ptype, GptFlags flags)
{
    double w = 0, h = 0;

    if (special_plot_size_is_set()) {
	w = (5.0 * special_width) / GP_WIDTH;
	h = (3.5 * special_height) / GP_HEIGHT;
	clear_special_plot_size();
    } else if (flags & GPT_LETTERBOX) {
	/* for time series */
	w = (5.0 * GP_LB_WIDTH) / GP_WIDTH;
	h = (3.5 * GP_LB_HEIGHT) / GP_HEIGHT;
    } else if (flags & GPT_XL) {
	/* large */
	w = (5.0 * GP_XL_WIDTH) / GP_WIDTH;
	h = (3.5 * GP_XL_HEIGHT) / GP_HEIGHT;
    } else if (flags & GPT_XXL) {
	/* extra large */
	w = h = (5.0 * GP_XXL_WIDTH) / GP_WIDTH;
    } else if (flags & GPT_XW) {
	/* extra wide */
	w = (5.0 * GP_XW_WIDTH) / GP_WIDTH;
	h = 3.5;
    } else if (ptype == PLOT_ROOTS || ptype == PLOT_QQ) {
	/* square plots */
	w = h = 3.5;
    }

    if (w > 0 && h > 0) {
	char size_str[32];

	gretl_push_c_numeric_locale();
	sprintf(size_str, " size %.2f,%.2f", w, h);
	gretl_pop_c_numeric_locale();
	strcat(s, size_str);
    }
}

static void append_gp_encoding (char *s)
{
    strcat(s, "\nset encoding utf8");
}

/* In pdf and eps term lines: should "dashed" be appended when
   "mono" is specified? Try experimenting?
*/

static char *eps_pdf_term_line (char *term_line,
				PlotType ptype,
				GptFlags flags,
				TermType ttype)
{
    gchar *font_string;
    const char *tname;
    int ptsize;

    ptsize = (ptype == PLOT_MULTI_BASIC)? 6 : 12;
    tname = (ttype == GP_TERM_EPS)? "epscairo" : "pdfcairo";

    font_string = write_other_font_string(ptsize);
    sprintf(term_line, "set term %s noenhanced font \"%s\"",
	    tname, font_string);
    g_free(font_string);

    maybe_set_eps_pdf_dims(term_line, ptype, flags);
    append_gp_encoding(term_line);

    return term_line;
}

static char *tex_term_line (char *term_line,
			    PlotType ptype,
			    GptFlags flags)
{
    if (gnuplot_has_tikz()) {
	strcpy(term_line, "set term tikz");
    } else {
	strcpy(term_line, "set term cairolatex");
    }

    append_gp_encoding(term_line);

    return term_line;
}

void plot_get_scaled_dimensions (int *width, int *height, double scale)
{
    *width *= scale;
    *height *= scale;

    /* PNG: round up to an even number of pixels if need be */
    if (*width % 2) *width += 1;
    if (*height % 2) *height += 1;
}

static void write_png_size_string (char *s, PlotType ptype,
				   GptFlags flags, double scale)
{
    int w = GP_WIDTH, h = GP_HEIGHT;

    if (special_plot_size_is_set()) {
	w = (int) special_width;
	h = (int) special_height;
	clear_special_plot_size();
    } else if (flags & GPT_LETTERBOX) {
	/* time series plots */
	w = GP_LB_WIDTH;
	h = GP_LB_HEIGHT;
    } else if (flags & GPT_XL) {
	/* large */
	w = GP_XL_WIDTH;
	h = GP_XL_HEIGHT;
    } else if (flags & GPT_XXL) {
	/* extra large */
	w = GP_XXL_WIDTH;
	h = GP_XXL_HEIGHT;
    } else if (flags & GPT_XW) {
	/* extra wide */
	w = GP_XW_WIDTH;
	h = GP_HEIGHT;
    } else if (ptype == PLOT_ROOTS || ptype == PLOT_QQ) {
	/* square plots */
	w = h = GP_SQ_SIZE;
    }

    if (scale != 1.0 && ptype != PLOT_GRIDPLOT) {
	plot_get_scaled_dimensions(&w, &h, scale);
    }

    *s = '\0';
    sprintf(s, " size %d,%d", w, h);
}

/* platform-specific on-screen term string */

static char *var_term_line (char *term_line, int ptype, GptFlags flags)
{
    char font_string[140];
    char size_string[16];
    const char *varterm;

#ifdef WIN32
    varterm = "windows";
#else
    if (gnuplot_has_wxt()) {
	varterm = "wxt";
    } else if (gnuplot_has_qt()) {
	varterm = "qt";
    } else {
	varterm = "x11";
    }
#endif

    *font_string = *size_string = '\0';
    write_png_font_string(font_string, "", ptype, NULL, 1.0);
    write_png_size_string(size_string, ptype, flags, 1.0);

    sprintf(term_line, "set term %s%s%s noenhanced",
	    varterm, font_string, size_string);
    append_gp_encoding(term_line);

    return term_line;
}

static char *real_png_term_line (char *term_line,
				 PlotType ptype,
				 GptFlags flags,
				 const char *specfont,
				 double scale)
{
    char ad_hoc_fontspec[128];
    char font_string[140];
    char size_string[16];

    *font_string = *size_string = *ad_hoc_fontspec = '\0';

    write_png_font_string(font_string, ad_hoc_fontspec,
			  ptype, specfont, scale);
    write_png_size_string(size_string, ptype, flags, scale);

    sprintf(term_line, "set term pngcairo%s%s noenhanced",
	    font_string, size_string);
    append_gp_encoding(term_line);

    if (*ad_hoc_fontspec != '\0') {
	strcat(term_line, "\n# fontspec: ");
	strcat(term_line, ad_hoc_fontspec);
    }

#if GP_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", term_line);
#endif

    return term_line;
}

static char *gretl_png_term_line (char *term_line,
				  PlotType ptype,
				  GptFlags flags)
{
    if (ptype == PLOT_GEOMAP) {
	return real_png_term_line(term_line, ptype, flags, NULL, 1.0);
    } else {
	double s = default_png_scale;

	return real_png_term_line(term_line, ptype, flags, NULL, s);
    }
}

/**
 * gretl_gnuplot_term_line:
 * @ttype: code for the gnuplot "terminal" type.
 * @ptype: indication of the sort of plot to be made, which
 * may make a difference to the color palette chosen.
 * @flags: plot option flags.
 * @font: if non-NULL, try to respect a specified font.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the specified "terminal".
 *
 * Returns: a static char * pointer.
 */

const char *gretl_gnuplot_term_line (TermType ttype,
				     PlotType ptype,
				     GptFlags flags,
				     const char *font)
{
    static char term_line[TERMLEN];

    *term_line = '\0';

    if (font != NULL && ad_hoc_font[0] == '\0') {
	strcpy(ad_hoc_font, font);
    }

    if (ttype == GP_TERM_PNG) {
	double s = default_png_scale;

	real_png_term_line(term_line, ptype, flags, NULL, s);
    } else if (ttype == GP_TERM_EPS || ttype == GP_TERM_PDF) {
	eps_pdf_term_line(term_line, ptype, flags, ttype);
    } else if (ttype == GP_TERM_EMF) {
	gretl_emf_term_line(term_line, ptype, flags);
    } else if (ttype == GP_TERM_SVG) {
	strcpy(term_line, "set term svg noenhanced");
	append_gp_encoding(term_line);
    } else if (ttype == GP_TERM_HTM) {
	strcpy(term_line, "set term canvas noenhanced");
	append_gp_encoding(term_line);
    } else if (ttype == GP_TERM_FIG) {
	strcpy(term_line, "set term fig");
	append_gp_encoding(term_line);
    } else if (ttype == GP_TERM_TEX) {
	tex_term_line(term_line, ptype, flags);
    } else if (ttype == GP_TERM_VAR) {
	var_term_line(term_line, ptype, flags);
    }

    return term_line;
}

const char *get_png_line_for_plotspec (const GPT_SPEC *spec)
{
    static char term_line[TERMLEN];

    real_png_term_line(term_line, spec->code, spec->flags,
		       spec->fontstr, spec->scale);
    return term_line;
}

void set_default_png_scale (double s)
{
    if (s >= 0.5 && s <= 2.0) {
	default_png_scale = s;
    }
}

double get_default_png_scale (void)
{
    return default_png_scale;
}

static void write_emf_font_string (char *fstr)
{
    const char *src = NULL;
    int stdsize = 16;
    int adhoc = 0;

    if (ad_hoc_font[0] != '\0') {
	src = ad_hoc_font;
	adhoc = 1;
    } else {
	src = gretl_png_font();
    }

    if (src != NULL) {
	char fname[128];
	int nf, fsize = 0;

	nf = split_graph_fontspec(src, fname, &fsize);
	if (nf == 2) {
	    if (adhoc) {
		/* go with what the user specified */
		sprintf(fstr, "font \"%s,%d\"", fname, fsize);
	    } else {
		/* adjust size to avoid tiny text? */
		fsize = (fsize <= 8)? 12 : 16;
		sprintf(fstr, "font \"%s,%d\"", fname, fsize);
	    }
	} else if (nf == 1) {
	    sprintf(fstr, "font \"%s,%d\"", fname, stdsize);
	}
	ad_hoc_font[0] = '\0';
    } else {
	sprintf(fstr, "font \"sans,%d\"", stdsize);
    }
}

static char *gretl_emf_term_line (char *term_line,
				  PlotType ptype,
				  GptFlags flags)
{
    gchar *size_string = NULL;
    char font_string[140];

    *font_string = '\0';
    write_emf_font_string(font_string);

    if (special_plot_size_is_set()) {
	size_string = g_strdup_printf("size %d,%d ",
				      (int) special_width,
				      (int) special_height);
	clear_special_plot_size();
    }

    if (flags & GPT_MONO) {
	strcat(term_line, "set term emf dash noenhanced ");
    } else {
	strcat(term_line, "set term emf color noenhanced ");
    }

    if (size_string != NULL) {
	strcat(term_line, size_string);
	g_free(size_string);
    }

    if (*font_string != '\0') {
	strcat(term_line, font_string);
    }

    append_gp_encoding(term_line);

    return term_line;
}

/**
 * plot_type_from_string:
 * @str: initial comment line from plot file.
 *
 * Returns: the special plot code corresponding to the initial
 * comment string in the plot file, or %PLOT_REGULAR if no special
 * comment is recognized.
 */

PlotType plot_type_from_string (const char *str)
{
    int i, len, ret = PLOT_REGULAR;

    for (i=1; ptinfo[i].pstr != NULL; i++) {
	len = strlen(ptinfo[i].pstr);
	if (!strncmp(str + 2, ptinfo[i].pstr, len)) {
	    ret = ptinfo[i].ptype;
	    break;
	}
    }

    return ret;
}

int write_plot_type_string (PlotType ptype, GptFlags flags, FILE *fp)
{
    int i, ret = 0;

    if (ptype == PLOT_GEOMAP) {
	/* handled specially */
	return 0;
    }

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (ptype == ptinfo[i].ptype) {
	    if (flags & GPT_XL) {
		fprintf(fp, "# %s (large)\n", ptinfo[i].pstr);
	    } else if (flags & GPT_XXL) {
		fprintf(fp, "# %s (extra-large)\n", ptinfo[i].pstr);
	    } else if (flags & GPT_XW) {
		fprintf(fp, "# %s (extra-wide)\n", ptinfo[i].pstr);
	    } else {
		fprintf(fp, "# %s\n", ptinfo[i].pstr);
	    }
	    ret = 1;
	    break;
	}
    }

    if (ret == 0 && (flags & GPT_XW)) {
	fputs("# extra-wide\n", fp);
    }

    if (get_local_decpoint() == ',') {
	/* is this right? */
	fputs("set decimalsign ','\n", fp);
    }

    return ret;
}

static void print_term_string (int ttype, PlotType ptype,
			       GptFlags flags, FILE *fp)
{
    char term_line[TERMLEN];

    *term_line = '\0';

    if (ttype == GP_TERM_EPS || ttype == GP_TERM_PDF) {
	eps_pdf_term_line(term_line, ptype, flags, ttype);
    } else if (ttype == GP_TERM_PNG) {
	gretl_png_term_line(term_line, ptype, flags);
    } else if (ttype == GP_TERM_EMF) {
	gretl_emf_term_line(term_line, ptype, flags);
    } else if (ttype == GP_TERM_FIG) {
	strcpy(term_line, "set term fig\nset encoding utf8");
    } else if (ttype == GP_TERM_SVG) {
	strcpy(term_line, "set term svg noenhanced\nset encoding utf8");
    } else if (ttype == GP_TERM_HTM) {
	strcpy(term_line, "set term canvas noenhanced\nset encoding utf8");
    } else if (ttype == GP_TERM_TEX) {
	tex_term_line(term_line, ptype, flags);
    }

    if (*term_line != '\0') {
	fprintf(fp, "%s\n", term_line);
	if (flags & GPT_MONO) {
	    fputs("set mono\n", fp);
	} else if (ptype != PLOT_GRIDPLOT) {
	    write_plot_line_styles(ptype, fp);
	}
    }
}

static int gretl_plot_count;
static int this_term_type;

/* recorder for filename given via --output=foo */
static char gnuplot_outname[FILENAME_MAX];

static int set_term_type_from_fname (const char *fname)
{
    if (has_suffix(fname, ".eps")) {
	this_term_type = GP_TERM_EPS;
    } else if (has_suffix(fname, ".ps")) {
	this_term_type = GP_TERM_EPS;
    } else if (has_suffix(fname, ".pdf")) {
	this_term_type = GP_TERM_PDF;
    } else if (has_suffix(fname, ".png")) {
	this_term_type = GP_TERM_PNG;
    } else if (has_suffix(fname, ".fig")) {
	this_term_type = GP_TERM_FIG;
    } else if (has_suffix(fname, ".emf")) {
	this_term_type = GP_TERM_EMF;
    } else if (has_suffix(fname, ".svg")) {
	this_term_type = GP_TERM_SVG;
    } else if (has_suffix(fname, ".html") ||
	       has_suffix(fname, ".htm")) {
	this_term_type = GP_TERM_HTM;
    } else if (has_suffix(fname, ".tex")) {
	this_term_type = GP_TERM_TEX;
    } else if (!strcmp(fname, "gnuplot")) {
	this_term_type = GP_TERM_VAR;
    }

    return this_term_type;
}

int specified_gp_output_format (void)
{
    return this_term_type;
}

void reset_plot_count (void)
{
    gretl_plot_count = 0;
}

/* We're writing the 'set output...' line for a gnuplot script.
   If @path is non-NULL we use it, otherwise we make a path
   using dotdir and "gretltmp.png".
*/

int write_plot_output_line (const char *path, FILE *fp)
{
    const char *fname;
    gchar *tmp = NULL;

#ifdef WIN32
    if (path == NULL) {
	tmp = gretl_make_dotpath("gretltmp.png");
    } else {
	tmp = g_strdup(path);
    }
    win32_forwardize(tmp);
    fname = tmp;
#else
    if (path == NULL) {
	fname = tmp = gretl_make_dotpath("gretltmp.png");
    } else {
	fname = path;
    }
#endif

    fprintf(fp, "set output \"%s\"\n", fname);
    g_free(tmp);

    return 0;
}

static FILE *gp_set_up_batch (char *fname,
			      PlotType ptype,
			      GptFlags flags,
			      const char *outspec,
			      int *err)
{
    int fmt = GP_TERM_NONE;
    FILE *fp = NULL;

    if (plot_output_to_buffer()) {
	this_term_type = GP_TERM_PLT;
	sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
	fp = gretl_mktemp(fname, "w");
    } else if (outspec != NULL) {
	/* user gave --output=<filename> */
	fmt = set_term_type_from_fname(outspec);
	if (fmt) {
	    /* input needs processing */
	    strcpy(gnuplot_outname, outspec);
	    gretl_maybe_prepend_dir(gnuplot_outname);
	    sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
	    fp = gretl_mktemp(fname, "w");
	} else {
	    /* just passing gnuplot commands through */
	    this_term_type = GP_TERM_PLT;
	    strcpy(fname, outspec);
	    gretl_maybe_prepend_dir(fname);
	    fp = gretl_fopen(fname, "w");
	}
    } else {
	/* auto-constructed gnuplot commands filename */
	this_term_type = GP_TERM_PLT;
	sprintf(fname, "gpttmp%02d.plt", ++gretl_plot_count);
	gretl_maybe_prepend_dir(fname);
	fp = gretl_fopen(fname, "w");
    }

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	gretl_set_path_by_name("plotfile", fname);
	if (*gnuplot_outname != '\0') {
	    /* write terminal/style/output lines */
	    print_term_string(fmt, ptype, flags, fp);
	    write_plot_output_line(gnuplot_outname, fp);
	} else if (ptype != PLOT_GRIDPLOT) {
	    /* just write style lines */
	    write_plot_line_styles(ptype, fp);
	}
	if (get_local_decpoint() == ',') {
	    fputs("set decimalsign ','\n", fp);
	}
    }

    return fp;
}

static char *iact_gpfile;

#define IACT_SPECIFY_TERM 1 /* experiment 2022-05-30 */

/* Set-up for an "interactive" plot: we open a file in the user's
   dotdir into which gnuplot commands will be written.  If we're
   running the GUI program this command file will eventually be used
   to create a PNG file for display in a gretl window; otherwise
   (gretlcli) the commands will eventually be sent to gnuplot for
   "direct" display (e.g. using the wxt or windows terminal).

   In this function we just open the file for writing; if in GUI
   mode insert a suitable PNG terminal line and output spec line;
   and write some header-type material including our line style
   specifications.
*/

static FILE *gp_set_up_interactive (char *fname, PlotType ptype,
				    GptFlags flags, int *err)
{
    int gui = gretl_in_gui_mode();
    FILE *fp = NULL;

    if (iact_gpfile != NULL) {
	fname = iact_gpfile;
	fp = gretl_fopen(fname, "wb");
	iact_gpfile = NULL;
    } else if (gui) {
	/* the filename should be unique */
	sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
	fp = gretl_mktemp(fname, "wb");
    } else {
	/* gretlcli: no need for uniqueness */
	sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
	fp = gretl_fopen(fname, "wb");
    }

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	gretl_set_path_by_name("plotfile", fname);
	if (gui) {
	    /* set up for PNG output */
	    fprintf(fp, "%s\n", gretl_gnuplot_term_line(GP_TERM_PNG, ptype,
							flags, NULL));
	    if (default_png_scale != 1.0) {
		gretl_push_c_numeric_locale();
		fprintf(fp, "# scale = %.1f\n", default_png_scale);
		fputs("# auto linewidth\n", fp);
		gretl_pop_c_numeric_locale();
	    }
	    write_plot_output_line(NULL, fp);
	} else if (ptype == PLOT_GEOMAP) {
	    fprintf(fp, "%s\n", gretl_gnuplot_term_line(GP_TERM_VAR, ptype,
							flags, NULL));
	} else {
#if IACT_SPECIFY_TERM
	    fprintf(fp, "%s\n", gretl_gnuplot_term_line(GP_TERM_VAR, ptype,
							flags, NULL));
#else
	    fputs("set termoption noenhanced\n", fp);
#endif
	}
	write_plot_type_string(ptype, flags, fp);
	if (ptype != PLOT_GRIDPLOT) {
	    write_plot_line_styles(ptype, fp);
	}
    }

    return fp;
}

#ifndef WIN32

static int gnuplot_too_old (void)
{
    static double gpv;
    int msg_done = 0;
    int ret = 0;

    if (gpv == 0.0) {
	gpv = gnuplot_version();
    }

    if (gpv == 0.0) {
	if (!msg_done) {
	    gretl_errmsg_sprintf("'%s': gnuplot is broken or too old: must be >= version 5.0",
				 gretl_gnuplot_path());
	}
	ret = 1;
    } else if (gpv < 5.0) {
	gretl_errmsg_sprintf("Gnuplot %g is too old: must be >= version 5.0", gpv);
	ret = 1;
    }

    return ret;
}

#endif

static int got_display_option (const char *s)
{
    return s != NULL && !strcmp(s, "display");
}

static int got_none_option (const char *s)
{
    return s != NULL && !strcmp(s, "none");
}

static const char *plot_output_option (PlotType p, int *pci, int *err)
{
    int mp_mode = gretl_multiplot_collecting();
    int ci = plot_ci;
    const char *s;

    /* clear any previous setting */
    set_plot_buffer_name(NULL);

    /* set a more specific command index based on
       the plot type, if applicable */

    if (p == PLOT_MULTI_BASIC) {
	ci = SCATTERS;
    } else if (p == PLOT_BOXPLOTS) {
	ci = BXPLOT;
    } else if (p == PLOT_FORECAST) {
	ci = FCAST;
    } else if (p == PLOT_CORRELOGRAM) {
	ci = CORRGM;
    } else if (p == PLOT_XCORRELOGRAM) {
	ci = XCORRGM;
    } else if (p == PLOT_PERIODOGRAM) {
	ci = PERGM;
    } else if (p == PLOT_HURST) {
	ci = HURST;
    } else if (p == PLOT_QQ) {
	ci = QQPLOT;
    } else if (p == PLOT_RANGE_MEAN) {
	ci = RMPLOT;
    } else if (p == PLOT_LEVERAGE) {
	ci = LEVERAGE;
    } else if (p == PLOT_FREQ_SIMPLE ||
	       p == PLOT_FREQ_NORMAL ||
	       p == PLOT_FREQ_GAMMA ||
	       p == PLOT_FREQ_DISCRETE) {
	ci = FREQ;
    } else if (p == PLOT_HEATMAP) {
	ci = CORR;
    } else if (p == PLOT_CUSUM) {
	ci = CUSUM;
    } else if (p == PLOT_KERNEL) {
	ci = KDPLOT;
    } else if (p == PLOT_GRIDPLOT) {
	ci = GRIDPLOT;
    }

    s = get_optval_string(ci, OPT_U);

    if (mp_mode && s != NULL) {
	if (gretl_function_depth() > 0 && !strcmp(s, "display")) {
	    /* let this pass: hansl functions that do plots
	       generally seem to default to "display"
	    */
	    s = NULL;
	} else if (p == PLOT_GEOMAP) {
	    s = NULL;
	} else {
	    *err = E_BADOPT;
	    return NULL;
	}
    }

    if (s != NULL && *s == '\0') {
	s = NULL;
    } else if (s == NULL && !mp_mode) {
	/* try for --outbuf=<strname> */
	s = get_optval_string(ci, OPT_b);
	if (s != NULL && *s == '\0') {
	    s = NULL;
	} else {
	    *err = set_plot_buffer_name(s);
	}
    }

    if (pci != NULL) {
	*pci = ci;
    }

    return s;
}

/**
 * open_plot_input_file:
 * @ptype: indication of the sort of plot to be made.
 * @flags: may inflect some characteristics of plot.
 * @err: location to receive error code.
 *
 * Opens a file into which gnuplot commands will be written.
 * Depending on the prospective use of the stream, some
 * header-type material may be written into it (the primary
 * case being when we're going to produce PNG output
 * for display in the gretl GUI). The prospective use is
 * figured out based on the program state, @ptype and
 * @flags.
 *
 * Returns: writable stream on success, %NULL on failure.
 */

FILE *open_plot_input_file (PlotType ptype, GptFlags flags, int *err)
{
    char fname[FILENAME_MAX] = {0};
    const char *outspec = NULL;
    int ci, interactive = 0;
    FILE *fp = NULL;

    /* ensure we have 'gnuplot_path' in place (file-scope static var) */
    if (*gnuplot_path == '\0') {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

#ifndef WIN32
    if (gnuplot_too_old()) {
	*err = E_EXTERNAL;
	return NULL;
    }
#endif

    /* initialize */
    this_term_type = GP_TERM_NONE;
    *gnuplot_outname = '\0';

    /* check for --output=whatever or --outbuf=whatever */
    outspec = plot_output_option(ptype, &ci, err);
    if (*err) {
	return NULL;
    }

    if (gretl_multiplot_collecting()) {
	interactive = 0;
    } else if (got_display_option(outspec)) {
	/* --output=display specified */
	interactive = 1;
    } else if (outspec != NULL) {
	/* --output=filename or --buffer=starvar specified */
	interactive = 0;
    } else if (flags & GPT_ICON) {
	interactive = 1;
    } else {
	/* default */
	interactive = !gretl_in_batch_mode();
    }

#if GP_DEBUG
    fprintf(stderr, "outspec = '%s', interactive = %d\n",
	    outspec == NULL ? "null" : outspec, interactive);
#endif

    if (interactive) {
	fp = gp_set_up_interactive(fname, ptype, flags, err);
    } else {
	fp = gp_set_up_batch(fname, ptype, flags, outspec, err);
    }

#if GP_DEBUG
    fprintf(stderr, "open_plot_input_file: '%s'\n", gretl_plotfile());
#endif

    return fp;
}

/* For plot-types that are generated by commands that also produce
   textual output: based on the program state and command option,
   figure out if a plot is actually wanted or not. In interactive mode
   the answer will be Yes unless --plot=none has been given; in batch
   mode the answer will be No unless --plot=display or --plot=fname
   has been given.
*/

int gnuplot_graph_wanted (PlotType ptype, gretlopt opt)
{
    const char *optname = NULL;
    int ret = 0;
    int err = 0;

    if (opt & (OPT_U | OPT_b)) {
	/* check for --plot or --outbuf option */
	optname = plot_output_option(ptype, NULL, &err);
    }

    if (got_none_option(optname)) {
	/* --plot=none specified */
	ret = 0;
    } else if (optname != NULL) {
	/* --plot=display or --plot=fname specified */
	ret = 1;
    } else if (gretl_multiplot_collecting()) {
	ret = 1;
    } else {
	/* defaults */
	ret = !gretl_in_batch_mode();
    }

    return ret;
}

/**
 * gnuplot_cleanup:
 *
 * Removes any temporary gnuplot input file written in
 * the user's dot directory.
 */

void gnuplot_cleanup (void)
{
    const char *p, *fname = gretl_plotfile();

    p = strstr(fname, "gpttmp");

    if (p != NULL) {
	int pnum;

	if (sscanf(p, "gpttmp%d.plt", &pnum) == 0) {
	    gretl_remove(fname);
	}
    }

    if (alt_sty != NULL) {
	/* clean up alternative gnuplot stylesheet */
	g_free(alt_sty);
	alt_sty = NULL;
    }
}

static int graph_file_written;

int graph_written_to_file (void)
{
    return graph_file_written;
}

static int graph_file_shown;

int graph_displayed (void)
{
    return graph_file_shown;
}

static void remove_old_png (void)
{
    gchar *tmp = gretl_make_dotpath("gretltmp.png");

    gretl_remove(tmp);
    g_free(tmp);
}

/*
 * gnuplot_make_graph:
 *
 * Executes gnuplot to produce a graph: in the gretl GUI
 * in interactive mode this will be a PNG file. In the
 * CLI program in interactive mode there will be a direct
 * call to gnuplot to display the graph. In batch mode
 * the type of file written depends on the options selected
 * by the user.
 *
 * Returns: 0 on success, non-zero on error.
 */

static int gnuplot_make_graph (void)
{
    char buf[MAXLEN];
    const char *fname = gretl_plotfile();
    int gui = gretl_in_gui_mode();
    int fmt, err = 0;

    graph_file_written = 0;
    graph_file_shown = 0;
    fmt = specified_gp_output_format();

    if (fmt == GP_TERM_PLT) {
	/* just the gnuplot commands are wanted */
	if (plot_output_to_buffer()) {
	    err = make_plot_commands_buffer(fname);
	} else {
	    graph_file_written = 1;
	}
	return err;
    } else if (fmt == GP_TERM_NONE && gui) {
	do_plot_bounding_box();
	/* ensure we don't get stale output */
	remove_old_png();
    }

#ifdef WIN32
    if (gui || fmt) {
	sprintf(buf, "\"%s\" \"%s\"", gretl_gnuplot_path(), fname);
    } else {
	/* gretlcli, interactive */
	sprintf(buf, "\"%s\" -persist \"%s\"", gretl_gnuplot_path(), fname);
    }
    err = gretl_spawn(buf);
#else /* !WIN32 */
    if (gui || fmt) {
	sprintf(buf, "%s \"%s\"", gretl_gnuplot_path(), fname);
    } else {
	/* gretlcli, interactive */
	sprintf(buf, "%s -persist \"%s\"", gretl_gnuplot_path(), fname);
    }
    err = gretl_spawn(buf);
#endif

#if GP_DEBUG
    fprintf(stderr, "gnuplot_make_graph:\n"
	    " command='%s', fmt = %d, err = %d\n", buf, fmt, err);
#endif

    if (fmt) {
	/* got a user-specified output format */
	if (err) {
	    /* leave the bad file for diagnostic purposes */
	    fprintf(stderr, "err = %d: bad file is '%s'\n", err, fname);
	} else {
	    /* remove the temporary input file */
	    gretl_remove(fname);
	    if (fmt == GP_TERM_VAR) {
		graph_file_shown = 1;
	    } else {
		gretl_set_path_by_name("plotfile", gnuplot_outname);
		graph_file_written = 1;
	    }
	}
    }

    return err;
}

/**
 * finalize_plot_input_file:
 * @fp: stream to which gnuplot commands have been written.
 *
 * Closes @fp and attempts to "make" the graph that it specifies.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int finalize_plot_input_file (FILE *fp)
{
    int err;

    if (fp != NULL) {
	fclose(fp);
	err = gnuplot_make_graph();
	if (!err) {
	    /* for the benefit of gretl_cmd_exec() in interact.c,
	       indicate that we actually produced a plot
	    */
	    set_plot_produced();
	}
    } else {
	err = 1;
    }

    return err;
}

/**
 * finalize_3d_plot_input_file:
 * @fp: stream to which gnuplot commands have been written.
 *
 * Closes @fp and alerts libgretl to the fact that an interactive
 * 3-D plot is wanted.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int finalize_3d_plot_input_file (FILE *fp)
{
    int err = 0;

    if (fp != NULL) {
	fclose(fp);
	graph_file_written = 1;
    } else {
	err = 1;
    }

    return err;
}

enum {
    GTITLE_VLS,
    GTITLE_RESID,
    GTITLE_AF,
    GTITLE_AFV
} graph_titles;

static void make_gtitle (gnuplot_info *gi, int code,
			 const char *s1, const char *s2,
			 FILE *fp)
{
    char depvar[VNAMELEN];
    gchar *title = NULL;

    switch (code) {
    case GTITLE_VLS:
	if (gi->fit == PLOT_FIT_OLS) {
	    title = g_strdup_printf(_("%s versus %s (with least squares fit)"),
				    s1, s2);
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    title = g_strdup_printf(_("%s versus %s (with inverse fit)"),
				    s1, s2);
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    title = g_strdup_printf(_("%s versus %s (with quadratic fit)"),
				    s1, s2);
	} else if (gi->fit == PLOT_FIT_CUBIC) {
	    title = g_strdup_printf(_("%s versus %s (with cubic fit)"),
				    s1, s2);
	}
	break;
    case GTITLE_RESID:
	if (strncmp(s1, "residual for ", 13) == 0 &&
	    gretl_scan_varname(s1 + 13, depvar) == 1) {
	    title = g_strdup_printf(_("Regression residuals (= observed - fitted %s)"),
				    depvar);
	}
	break;
    case GTITLE_AF:
	title = g_strdup_printf(_("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || (gi->flags & GPT_TS)) {
	    title = g_strdup_printf(_("Actual and fitted %s"), s1);
	} else {
	    title = g_strdup_printf(_("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	break;
    }

    if (title != NULL) {
	fprintf(fp, "set title \"%s\"\n", title);
	g_free(title);
    }
}

static void print_axis_label (char axis, const char *s, FILE *fp)
{
    if (strchr(s, '\'')) {
	fprintf(fp, "set %clabel \"%s\"\n", axis, s);
    } else {
	fprintf(fp, "set %clabel '%s'\n", axis, s);
    }
}

static int literal_line_out (const char *s, int len,
			     FILE *fp)
{
    char *q, *p = malloc(len + 1);
    int n, warn = 0;

    if (p != NULL) {
	*p = '\0';
	strncat(p, s, len);
	q = p + strspn(p, " \t");
	n = strlen(q);
	if (n > 0) {
	    /* note: allow "set termoption ..." */
	    if (!strncmp(q, "set term ", 9)) {
		warn = 1;
	    } else {
		fputs(q, fp);
		if (q[n-1] != '\n') {
		    fputc('\n', fp);
		}
	    }
	}
	free(p);
    }

    return warn;
}

static int gnuplot_literal_from_string (const char *s,
					FILE *fp)
{
    const char *p;
    int braces = 1;
    int wi, warn = 0;

    s += strspn(s, " \t{");
    p = s;

    fputs("# start literal lines\n", fp);

    while (*s) {
	if (*s == '{') {
	    braces++;
	} else if (*s == '}') {
	    braces--;
	}
	if (braces == 0) {
	    break;
	}
	if (*s == ';') {
	    wi = literal_line_out(p, s - p, fp);
	    if (wi && !warn) {
		warn = 1;
	    }
	    p = s + 1;
	}
	s++;
    }

    fputs("# end literal lines\n", fp);

    return warn;
}

/* Alternative (undocumented!) means of supplying "literal"
   lines to the "gnuplot" command (as opposed to the time-
   honored "{...}" mechanism). Syntax is

   gnuplot <args> --tweaks=<name-of-array-of-strings>

   FIXME: document this or get rid of it! Although this
   mechanism has something going for it, maybe it's too
   late to substitute it for the old method.
*/

static char **literal_strings_from_opt (int ci, int *ns,
					int *real_ns)
{
    const char *aname = get_optval_string(ci, OPT_K);
    char **S = NULL;

    *ns = *real_ns = 0;

    if (aname != NULL) {
	GretlType type;
	gretl_array *A;
	int i;

	A = user_var_get_value_and_type(aname, &type);

	if (A != NULL && type == GRETL_TYPE_ARRAY) {
	    S = gretl_array_get_strings(A, ns);
	    if (*ns > 0) {
		for (i=0; i<*ns; i++) {
		    if (S[i] != NULL && S[i][0] != '\0') {
			*real_ns += 1;
		    }
		}
	    }
	}
    }

    return S;
}

static int gnuplot_literal_from_opt (int ci, FILE *fp)
{
    char *s, **S;
    int ns, real_ns;
    int warn = 0;

    S = literal_strings_from_opt(ci, &ns, &real_ns);

    if (real_ns > 0) {
	int i, n;

	fputs("# start literal lines\n", fp);

	for (i=0; i<ns; i++) {
	    if (S[i] != NULL) {
		s = S[i];
		s += strspn(s, " \t");
		n = strlen(s);
		if (n > 0) {
		    if (!strncmp(s, "set term", 8)) {
			warn = 1;
		    } else {
			fputs(s, fp);
			if (s[n-1] != '\n') {
			    fputc('\n', fp);
			}
		    }
		}
	    }
	}

	fputs("# end literal lines\n", fp);
    }

    return warn;
}

int print_gnuplot_literal_lines (const char *s, int ci,
				 gretlopt opt, FILE *fp)
{
    if (s != NULL && *s != '\0') {
	gnuplot_literal_from_string(s, fp);
    } else if (opt & OPT_K) {
	gnuplot_literal_from_opt(ci, fp);
    }

    return 0;
}

static void print_extra_literal_lines (char **S,
				       int ns,
				       FILE *fp)
{
    int i, n;

    for (i=0; i<ns; i++) {
	if (S[i] != NULL) {
	    n = strlen(S[i]);
	    if (n > 0) {
		fputs(S[i], fp);
		if (S[i][n-1] != '\n') {
		    fputc('\n', fp);
		}
	    }
	}
    }
}

static int loess_plot (gnuplot_info *gi, const char *literal,
		       const DATASET *dset)
{
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;
    gretl_matrix *yh = NULL;
    int xno, yno = gi->list[1];
    const double *yvar = dset->Z[yno];
    const double *xvar;
    FILE *fp = NULL;
    gchar *title = NULL;
    int t, T, d = 1;
    double q = 0.5;
    int err = 0;

    if (gi->x != NULL) {
	xno = 0;
	xvar = gi->x;
    } else {
	xno = gi->list[2];
	xvar = dset->Z[xno];
    }

    err = graph_list_adjust_sample(gi->list, gi, dset, 2);
    if (!err && gi->list[0] > 2) {
	err = E_DATA;
    }
    if (err) {
	return err;
    }

    fp = open_plot_input_file(PLOT_REGULAR, gi->flags, &err);
    if (err) {
	return E_FOPEN;
    }

    err = gretl_plotfit_matrices(yvar, xvar, PLOT_FIT_LOESS,
				 gi->t1, gi->t2, &y, &x);

    if (!err) {
	err = sort_pairs_by_x(x, y, NULL, NULL); /* markers! */
    }

    if (!err) {
	yh = loess_fit(x, y, d, q, OPT_R, &err);
    }

    if (err) {
	fclose(fp);
	goto bailout;
    }

    if (xno > 0) {
	const char *s1 = series_get_graph_name(dset, yno);
	const char *s2 = series_get_graph_name(dset, xno);

	title = g_strdup_printf(_("%s versus %s (with loess fit)"), s1, s2);
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	fprintf(fp, "set title \"%s\"\n", title);
	g_free(title);
	print_axis_label('y', s1, fp);
	print_axis_label('x', s2, fp);
    } else {
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	print_axis_label('y', series_get_graph_name(dset, yno), fp);
    }

    print_auto_fit_string(PLOT_FIT_LOESS, fp);

    print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    title = g_strdup_printf(_("loess fit, d = %d, q = %g"), d, q);
    fprintf(fp, " '-' using 1:2 title \"%s\" w lines\n", title);
    g_free(title);

    T = gretl_vector_get_length(yh);

    gretl_push_c_numeric_locale();

    for (t=0; t<T; t++) {
	fprintf(fp, "%.10g %.10g\n", x->val[t], y->val[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<T; t++) {
	fprintf(fp, "%.10g %.10g\n", x->val[t], yh->val[t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(x);
    gretl_matrix_free(yh);
    clear_gpinfo(gi);

    return err;
}

static int get_fitted_line (gnuplot_info *gi,
			    const DATASET *dset,
			    gchar **targ)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    const double *yvar, *xvar = NULL;
    double x0, s2 = 0, *ps2 = NULL;
    int allow_err = 0;
    int err = 0;

    if (gi->x != NULL && (dset->pd == 1 || dset->pd == 4 || dset->pd == 12)) {
	/* starting value of time index */
	x0 = gi->x[gi->t1];
    } else {
	xvar = dset->Z[gi->list[2]];
	x0 = NADBL;
    }

    yvar = dset->Z[gi->list[1]];

    if (gi->fit == PLOT_FIT_NONE) {
	/* Doing first-time automatic OLS: we want to check for
	   statistical significance of the slope coefficient
	   to see if it's worth drawing the fitted line, so
	   we have to allocate the variance matrix; otherwise
	   we only need the coefficients.
	*/
	V = gretl_matrix_alloc(2, 2);
	if (V == NULL) {
	    return E_ALLOC;
	}
	ps2 = &s2;
	/* if the fit attempt is automatic we'll allow it to
	   fail without aborting the plot */
	allow_err = 1;
    }

    err = gretl_plotfit_matrices(yvar, xvar, gi->fit,
				 dset->t1, dset->t2,
				 &y, &X);

    if (!err) {
	int  k = 2;

	if (gi->fit == PLOT_FIT_CUBIC) {
	    k = 4;
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    k = 3;
	}
	b = gretl_column_vector_alloc(k);
	if (b == NULL) {
	    err = E_ALLOC;
	} else if (gi->fit == PLOT_FIT_LOGLIN) {
	    ps2 = &s2;
	}
    }

    if (!err) {
	err = gretl_matrix_ols(y, X, b, V, NULL, ps2);
    }

    if (!err && gi->fit == PLOT_FIT_NONE) {
	/* the "automatic" case */
	double pv, v = gretl_matrix_get(V, 1, 1);
	int T = gretl_vector_get_length(y);

	pv = student_pvalue_2(T - 2, b->val[1] / sqrt(v));
	/* show the line if the two-tailed p-value for the slope coeff
	   is less than 0.1, otherwise discard it */
	if (pv < 0.10) {
	    gi->fit = PLOT_FIT_OLS;
	}
    }

    if (!err && gi->fit != PLOT_FIT_NONE) {
	GPT_LINE line = {0};
	double pd = dset->pd;

	if (gi->fit == PLOT_FIT_LOGLIN) {
	    b->val[0] += s2 / 2;
	}
	set_plotfit_line(&line, gi->fit, b->val, x0, pd);
	*targ = g_strdup_printf("%s title \"%s\" w lines\n",
				line.formula, line.title);
	g_free(line.formula);
	g_free(line.title);
	gi->flags |= GPT_AUTO_FIT;
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    if (err && allow_err) {
	err = 0;
    }

    return err;
}

/* support the "fit" options for a single time-series plot */

static int time_fit_plot (gnuplot_info *gi, const char *literal,
			  const DATASET *dset)
{
    int yno = gi->list[1];
    const double *yvar = dset->Z[yno];
    gchar *fitline = NULL;
    FILE *fp = NULL;
    PRN *prn;
    int t, err = 0;

    if (gi->x == NULL) {
	return E_DATA;
    }

    err = graph_list_adjust_sample(gi->list, gi, dset, 1);
    if (err) {
	return err;
    }

    err = get_fitted_line(gi, dset, &fitline);
    if (err) {
	return err;
    }

    gi->flags |= GPT_LETTERBOX;

    fp = open_plot_input_file(PLOT_REGULAR, gi->flags, &err);
    if (err) {
	g_free(fitline);
	return err;
    }

    prn = gretl_print_new_with_stream(fp);

    if (prn != NULL) {
	make_time_tics(gi, dset, 0, NULL, prn);
	gretl_print_detach_stream(prn);
	gretl_print_destroy(prn);
    }

    print_keypos_string(GP_KEY_LEFT_TOP, fp);
    print_axis_label('y', series_get_graph_name(dset, yno), fp);
    print_auto_fit_string(gi->fit, fp);
    print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w lines, \\\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, " %s", fitline);
    g_free(fitline);

    for (t=gi->t1; t<=gi->t2; t++) {
	if (gi->flags & GPT_TIMEFMT) {
	    fprintf(fp, "%.0f %.10g\n", gi->x[t], yvar[t]);
	} else {
	    fprintf(fp, "%.10g %.10g\n", gi->x[t], yvar[t]);
	}
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    clear_gpinfo(gi);

    return err;
}

static int check_tic_labels (double vmin, double vmax,
			     gnuplot_info *gi,
			     char axis)
{
    char s1[32], s2[32];
    int d, err = 0;

    for (d=6; d<12; d++) {
	sprintf(s1, "%.*g", d, vmin);
	sprintf(s2, "%.*g", d, vmax);
	if (strcmp(s1, s2)) {
	    break;
	}
    }

    if (d > 6) {
	if (axis == 'x') {
	    sprintf(gi->xfmt, "%% .%dg", d+1);
	    if (vmax > vmin) {
		sprintf(gi->xtics, "%.*g %#.6g", d+1, vmin,
			(vmax - vmin)/ 4.0);
	    }
	} else {
	    sprintf(gi->yfmt, "%% .%dg", d+1);
	}
    }

    return err;
}

static void check_y_tics (gnuplot_info *gi, const double **Z,
			  FILE *fp)
{
    double ymin, ymax;

    *gi->yfmt = '\0';

    gretl_minmax(gi->t1, gi->t2, Z[gi->list[1]], &ymin, &ymax);
    check_tic_labels(ymin, ymax, gi, 'y');

    if (*gi->yfmt != '\0') {
	fprintf(fp, "set format y \"%s\"\n", gi->yfmt);
    }
}

/* Find the minimum and maximum x-axis values and construct the gnuplot
   x-range.  We have to be a bit careful here to include only values
   that will actually display on the plot, i.e. x-values that are
   accompanied by at least one non-missing y-axis value.

   In the case of a "factorized" plot the factor variable must also
   be non-missing in order to include a given data point.

   We also have to avoid creating an "empty x range" that will choke
   gnuplot.
*/

static void print_x_range_from_list (gnuplot_info *gi,
				     const DATASET *dset,
				     const int *list,
				     FILE *fp)
{
    const double *x, *d = NULL;
    int k, l0 = list[0];

    if (gi->flags & GPT_DUMMY) {
	/* the factor variable comes last and the x variable
	   is in second-last place */
	d = dset->Z[list[l0]];
	k = l0 - 1;
    } else {
	/* the x variable comes last in the list */
	k = l0;
    }

    x = dset->Z[list[k]];

    if (gretl_isdummy(gi->t1, gi->t2, x)) {
	fputs("set xrange [-1:2]\n", fp);
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
	gi->xrange = 3;
    } else {
	double xmin, xmin0 = NADBL;
	double xmax, xmax0 = NADBL;
	int t, i, vy, obs_ok;

	for (t=gi->t1; t<=gi->t2; t++) {
	    obs_ok = 0;
	    if (!na(x[t]) && (d == NULL || !na(d[t]))) {
		for (i=1; i<k; i++) {
		    vy = list[i];
		    if (!na(dset->Z[vy][t])) {
			/* got x obs and at least one y obs */
			obs_ok = 1;
			break;
		    }
		}
	    }
	    if (obs_ok) {
		if (na(xmin0) || x[t] < xmin0) {
		    xmin0 = x[t];
		}
		if (na(xmax0) || x[t] > xmax0) {
		    xmax0 = x[t];
		}
	    }
	}

	gi->xrange = xmax0 - xmin0;

	if (gi->xrange == 0.0) {
	    /* construct a non-empty range */
	    xmin = xmin0 - 0.5;
	    xmax = xmin0 + 0.5;
	} else {
	    xmin = xmin0 - gi->xrange * .025;
	    if (xmin0 >= 0.0 && xmin < 0.0) {
		xmin = 0.0;
	    }
	    xmax = xmax0 + gi->xrange * .025;
	}

	fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
	gi->xrange = xmax - xmin;
	check_tic_labels(xmin0, xmax0, gi, 'x');
    }
}

void print_x_range (gnuplot_info *gi, FILE *fp)
{
    if (gretl_isdummy(gi->t1, gi->t2, gi->x)) {
	fputs("set xrange [-1:2]\n", fp);
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
	gi->xrange = 3;
    } else {
	double xmin0, xmin, xmax0, xmax;

	gretl_minmax(gi->t1, gi->t2, gi->x, &xmin0, &xmax0);
	gi->xrange = xmax0 - xmin0;
	xmin = xmin0 - gi->xrange * .025;
	if (xmin0 >= 0.0 && xmin < 0.0) {
	    xmin = 0.0;
	}
	xmax = xmax0 + gi->xrange * .025;
	fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
	gi->xrange = xmax - xmin;
    }
}

/* two or more y vars plotted against some x: test to see if we want
   to use two y axes */

void check_for_yscale (gnuplot_info *gi, const double **Z, int *oddman)
{
    double ymin[6], ymax[6];
    double ratio;
    int lmax = gi->list[0];
    int i, j, oddcount;

#if GP_DEBUG
    fprintf(stderr, "gnuplot: doing check_for_yscale, listlen %d\n",
	    gi->list[0]);
#endif

    if (gi->flags & GPT_IDX) {
	/* do this only if we haven't added a 0 at the end of
	   the list */
	if (gi->list[lmax] != 0) {
	    lmax++;
	}
    }

    /* find minima, maxima of the y-axis vars */
    for (i=1; i<lmax; i++) {
	gretl_minmax(gi->t1, gi->t2, Z[gi->list[i]],
		     &ymin[i-1], &ymax[i-1]);
    }

    gi->flags &= ~GPT_Y2AXIS;

    for (i=lmax-1; i>0; i--) {
	oddcount = 0;
	for (j=1; j<lmax; j++) {
	    if (j != i) {
		ratio = ymax[i-1] / ymax[j-1];
		if (ratio > 5.0 || ratio < 0.2) {
		    gi->flags |= GPT_Y2AXIS;
		    oddcount++;
		}
	    }
	}
	if (oddcount == lmax - 2) {
	    /* series at list position i differs considerably in scale
	       from all the others in the list */
	    *oddman = i;
	    break;
	}
    }

    if (*oddman == 0) {
	gi->flags &= ~GPT_Y2AXIS;
    }
}

static int print_gp_dummy_data (gnuplot_info *gi,
				const DATASET *dset,
				FILE *fp)
{
    const double *d = dset->Z[gi->list[3]];
    const double *y = dset->Z[gi->list[1]];
    double xt, yt;
    int i, t, n;

    n = gretl_vector_get_length(gi->dvals);

    for (i=0; i<n; i++) {
	for (t=gi->t1; t<=gi->t2; t++) {
	    if (gi->x != NULL) {
		xt = gi->x[t];
	    } else {
		xt = dset->Z[gi->list[2]][t];
		if (na(xt)) {
		    continue;
		}
	    }
	    if (na(d[t])) {
		continue;
	    }
	    yt = (d[t] == gi->dvals->val[i])? y[t] : NADBL;
	    if (na(yt)) {
		fprintf(fp, "%.10g %s\n", xt, GPNA);
	    } else {
		fprintf(fp, "%.10g %.10g", xt, yt);
		if (!(gi->flags & GPT_TS)) {
		    if (dset->markers) {
			fprintf(fp, " # %s", dset->S[t]);
		    } else if (dataset_is_time_series(dset)) {
			char obs[OBSLEN];

			ntolabel(obs, t, dset);
			fprintf(fp, " # %s", obs);
		    }
		}
		fputc('\n', fp);
	    }
	}
	fputs("e\n", fp);
    }

    return 0;
}

/* for printing panel time-series graph: insert a discontinuity
   between the panel units */

static void
maybe_print_panel_jot (int t, const DATASET *dset, FILE *fp)
{
    char obs[OBSLEN];
    int maj, min;

    ntolabel(obs, t, dset);
    sscanf(obs, "%d:%d", &maj, &min);
    if (maj > 1 && min == 1) {
	fprintf(fp, "%g %s\n", t + 0.5, GPNA);
    }
}

/* sanity check for totally empty graph */

static int
all_graph_data_missing (const int *list, int t, const double **Z)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (!na(Z[list[i]][t])) {
	    return 0;
	}
    }

    return 1;
}

static int use_impulses (gnuplot_info *gi)
{
    if (gi->withlist != NULL) {
	int i;

	for (i=1; i<=gi->withlist[0]; i++) {
	    if (gi->withlist[i] == W_IMPULSES) {
		return 1;
	    }
	}
    }

    return 0;
}

static int use_lines (gnuplot_info *gi)
{
    if (gi->withlist != NULL) {
	int i;

	for (i=1; i<=gi->withlist[0]; i++) {
	    if (gi->withlist[i] == W_LINES) {
		return 1;
	    }
	}
    }

    return 0;
}

/* list of series IDs for which we should skip observations
   with NAs when printing the plot data */
static int *na_skiplist;

static void print_gp_data (gnuplot_info *gi, const DATASET *dset,
			   FILE *fp)
{
    int n = gi->t2 - gi->t1 + 1;
    double offset = 0.0;
    int datlist[3];
    int lmax, ynum = 2;
    int nomarkers = 0;
    int i, t;

    /* multi impulse plot? calculate offset for lines */
    if (use_impulses(gi) && gi->list[0] > 2) {
	offset = 0.10 * gi->xrange / n;
    }

    if (gi->x != NULL) {
	lmax = gi->list[0] - 1;
	datlist[0] = 1;
	ynum = 1;
    } else {
	lmax = gi->list[0] - 1;
	datlist[0] = 2;
	datlist[1] = gi->list[gi->list[0]];
    }

    if (use_impulses(gi) || use_lines(gi)) {
	nomarkers = 1;
    }

    /* loop across the variables, printing x then y[i] for each i */

    for (i=1; i<=lmax; i++) {
	double xoff = offset * (i - 1);

	datlist[ynum] = gi->list[i];

	for (t=gi->t1; t<=gi->t2; t++) {
	    const char *label = NULL;
	    char obs[OBSLEN];

	    if (in_gretl_list(na_skiplist, datlist[ynum]) &&
		na(dset->Z[datlist[ynum]][t])) {
		continue;
	    } else if (gi->x == NULL &&
		       all_graph_data_missing(gi->list, t, (const double **) dset->Z)) {
		continue;
	    }
	    if (!(gi->flags & GPT_TS) && i == 1) {
		if (dset->markers) {
		    label = dset->S[t];
		} else if (!nomarkers && dataset_is_time_series(dset)) {
		    ntolabel(obs, t, dset);
		    label = obs;
		}
	    }
	    if ((gi->flags & GPT_TS) && dset->structure == STACKED_TIME_SERIES) {
		maybe_print_panel_jot(t, dset, fp);
	    }
	    printvars(fp, t, datlist, dset, gi, label, xoff);
	}

	fputs("e\n", fp);
    }
}

static int
gpinfo_init (gnuplot_info *gi, gretlopt opt, const int *list,
	     const char *literal, const DATASET *dset)
{
    int l0 = list[0];
    int err = 0;

    gi->withlist = NULL;
    gi->yformula = NULL;
    gi->x = NULL;
    gi->list = NULL;
    gi->dvals = NULL;
    gi->band = 0;

    err = get_gp_flags(gi, opt, list, dset);
    if (err) {
	return err;
    }

    if (gi->band) {
	/* force single y-axis in "band" case */
	opt |= OPT_Y;
    }

    if (gi->fit == PLOT_FIT_NONE) {
	gi->flags |= GPT_TS; /* may be renounced later */
    }

    if (dset->t2 - dset->t1 + 1 <= 0) {
	/* null sample range */
	return E_DATA;
    }

    gi->t1 = dset->t1;
    gi->t2 = dset->t2;
    gi->xrange = 0.0;
    gi->xtics[0] = '\0';
    gi->xfmt[0] = '\0';
    gi->yfmt[0] = '\0';

    if (l0 < 2 && !(gi->flags & GPT_IDX)) {
	return E_ARGS;
    }

    if ((gi->flags & GPT_DUMMY) && (gi->flags & GPT_IDX)) {
	return E_BADOPT;
    }

    gi->list = gretl_list_copy(list);
    if (gi->list == NULL) {
	return E_ALLOC;
    }

    if ((l0 > 2 || (l0 > 1 && (gi->flags & GPT_IDX))) &&
	l0 < 7 && !(gi->flags & GPT_RESIDS) && !(gi->flags & GPT_FA)
	&& !(gi->flags & GPT_DUMMY) && !(opt & OPT_Y)) {
	/* FIXME GPT_XYZ ? */
	/* allow probe for using two y axes */
#if GP_DEBUG
	fprintf(stderr, "l0 = %d, setting y2axis probe\n", l0);
#endif
	gi->flags |= GPT_Y2AXIS;
    }

    if ((gi->flags & GPT_FA) && literal != NULL &&
	!strncmp(literal, "yformula: ", 10)) {
	/* fitted vs actual plot with fitted given by formula */
	gi->yformula = literal + 10;
    }

    if (literal != NULL && strstr(literal, "set style data")) {
	gi->flags |= GPT_DATA_STYLE;
    }

#if GP_DEBUG
    if (gi->flags) {
	print_gnuplot_flags(gi->flags, 1);
    }
#endif

    return 0;
}

void clear_gpinfo (gnuplot_info *gi)
{
    free(gi->list);
    gretl_matrix_free(gi->dvals);
    free(gi->withlist);
}

#if GP_DEBUG

static void print_gnuplot_flags (int flags, int revised)
{
    if (revised) {
	fprintf(stderr, "*** gnuplot flags after initial revision:\n");
    } else {
	fprintf(stderr, "*** gnuplot() called with flags:\n");
    }

    if (flags & GPT_RESIDS) {
	fprintf(stderr, " GPT_RESIDS\n");
    }
    if (flags & GPT_FA) {
	fprintf(stderr, " GPT_FA\n");
    }
    if (flags & GPT_DUMMY) {
	fprintf(stderr, " GPT_DUMMY\n");
    }
    if (flags & GPT_XYZ) {
	fprintf(stderr, " GPT_XYZ\n");
    }
    if (flags & GPT_FIT_OMIT) {
	fprintf(stderr, " GPT_FIT_OMIT\n");
    }
    if (flags & GPT_DATA_STYLE) {
	fprintf(stderr, " GPT_DATA_STYLE\n");
    }
    if (flags & GPT_IDX) {
	fprintf(stderr, " GPT_IDX\n");
    }
    if (flags & GPT_TS) {
	fprintf(stderr, " GPT_TS\n");
    }
    if (flags & GPT_Y2AXIS) {
	fprintf(stderr, " GPT_Y2AXIS\n");
    }
    if (flags & GPT_AUTO_FIT) {
	fprintf(stderr, " GPT_AUTO_FIT\n");
    }
    if (flags & GPT_FIT_HIDDEN) {
	fprintf(stderr, " GPT_FIT_HIDDEN\n");
    }
}

#endif

static void set_lwstr (char *s)
{
    if (default_png_scale > 1.0) {
	strcpy(s, " lw 2");
    } else {
	*s = '\0';
    }
}

void set_plot_withstr (gnuplot_info *gi, int i, char *str)
{
    if (gi->flags & GPT_DATA_STYLE) {
	*str = '\0';
    } else if (gi->withlist != NULL) {
	int withval = W_POINTS;

	if (i > 0 && i <= gi->withlist[0]) {
	    withval = gi->withlist[i];
	}

	if (withval == W_LINES) {
	    strcpy(str, "w lines");
	} else if (withval == W_IMPULSES) {
	    strcpy(str, "w impulses");
	} else if (withval == W_LP) {
	    strcpy(str, "w linespoints");
	} else if (withval == W_BOXES) {
	    strcpy(str, "w boxes");
	} else if (withval == W_STEPS) {
	    strcpy(str, "w steps");
	} else {
	    strcpy(str, "w points");
	}
    } else if (gi->flags & GPT_LINES) {
        strcpy(str, "w lines");
    } else if (gi->flags & GPT_IMPULSES) {
	strcpy(str, "w impulses");
    } else if (gi->flags & GPT_STEPS) {
	strcpy(str, "w steps");
    } else {
	strcpy(str, "w points");
    }
}

int graph_list_adjust_sample (int *list,
			      gnuplot_info *ginfo,
			      const DATASET *dset,
			      int listmin)
{
    int t1min = ginfo->t1;
    int t2max = ginfo->t2;
    int i, t, vi, t_ok;
    int err = 0;

    for (t=t1min; t<=t2max; t++) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && !na(dset->Z[vi][t])) {
		t_ok = 1;
		break;
	    }
	}
	if (t_ok) {
	    break;
	}
	t1min++;
    }

    for (t=t2max; t>t1min; t--) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && !na(dset->Z[vi][t])) {
		t_ok = 1;
		break;
	    }
	}
	if (t_ok) {
	    break;
	}
	t2max--;
    }

    if (t2max > t1min) {
	for (i=1; i<=list[0]; i++) {
	    int all_missing = 1;

	    vi = list[i];
	    for (t=t1min; t<=t2max; t++) {
		if (!na(dset->Z[vi][t])) {
		    all_missing = 0;
		    break;
		}
	    }
	    if (all_missing) {
		gretl_list_delete_at_pos(list, i);
		i--;
	    }
	}
    }

    ginfo->t1 = t1min;
    ginfo->t2 = t2max;

    if (ginfo->t1 >= ginfo->t2 || list[0] < listmin) {
	err = E_MISSDATA;
    }

    return err;
}

static int timefmt_useable (const DATASET *dset)
{
    return dated_daily_data(dset) || dated_weekly_data(dset);
}

static int maybe_add_plotx (gnuplot_info *gi, int time_fit,
			    const DATASET *dset)
{
    gretlopt xopt = OPT_NONE;
    int k = gi->list[0];
    int add0 = 0;

    /* are we really doing a time-series plot? */
    if (k > 1 && !strcmp(dset->varname[gi->list[k]], "time")) {
	; /* yes */
    } else if (gi->flags & GPT_IDX) {
	add0 = 1; /* yes */
    } else {
	/* no: get out */
	gi->flags &= ~GPT_TS;
	return 0;
    }

    if (!time_fit && timefmt_useable(dset)) {
	/* experimental */
	gi->flags |= GPT_TIMEFMT;
	xopt = OPT_T;
    }

    gi->x = gretl_plotx(dset, xopt);
    if (gi->x == NULL) {
	return E_ALLOC;
    }

    /* a bit nasty, but add a dummy list entry for
       the 'virtual' plot variable
    */
    if (add0) {
	gretl_list_append_term(&gi->list, 0);
	if (gi->list == NULL) {
	    return E_ALLOC;
	}
    }

#if GP_DEBUG
    fprintf(stderr, "maybe_add_plotx: gi->x at %p\n",
	    (void *) gi->x);
    printlist(gi->list, "gi->list");
#endif

    return 0;
}

void gnuplot_missval_string (FILE *fp)
{
    fputs("set datafile missing \"?\"\n", fp);
}

static void graph_month_name (char *mname, int m)
{
    struct tm mt;

    mt.tm_sec = 0;
    mt.tm_min = 0;
    mt.tm_hour = 0;
    mt.tm_mday = 1;
    mt.tm_mon = m - 1;
    mt.tm_year = 100;

    strftime(mname, 7, "%b", &mt);

    if (!g_utf8_validate(mname, -1, NULL)) {
	/* we might be in a non-UTF-8 locale */
	gchar *tmp;
	gsize bytes;

	tmp = g_locale_to_utf8(mname, -1, NULL, &bytes, NULL);
	if (tmp != NULL) {
	    strcpy(mname, tmp);
	    g_free(tmp);
	}
    }
}

/* for short daily time-series plots: write month names
   into the xtics */

static void make_named_month_tics (const gnuplot_info *gi, double yrs,
				   PRN *prn)
{
    double t0 = gi->x[gi->t1];
    double t1 = gi->x[gi->t2];
    double x, tw = 1.0/12;
    int i, m, n = 0;
    char mname[16];
    int notfirst = 0;
    int scale = (int) (yrs * 1.5);

    if (scale == 0) {
	scale = 1;
    }

    t0 += (1.0 - (t0 - floor(t0)) * 12.0) / 12.0;
    for (x=t0; x<t1; x+=tw) n++;

    x = (t0 - floor(t0)) * 12;
    m = 1 + ((x - floor(x) > .8)? ceil(x) : floor(x));
    if (m > 12) m -= 12;

    pputs(prn, "set xtics (");
    x = t0;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	if (m == 1) {
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%4.0f\" %.10g", x, x);
	    notfirst = 1;
	} else if ((scale == 1) || (m % scale == 1)) {
	    graph_month_name(mname, m);
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%s\" %.10g", mname, x);
	    notfirst = 1;
	}
	m++;
	x += tw;
	if (m > 12) m -= 12;
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

static int continues_unit (const DATASET *dset, int t)
{
    return t / dset->pd == (t-1) / dset->pd;
}

/* Below: we're making a combined time series plot for panel data.
   That is, time series for unit 1, followed by time series for unit
   2, etc.  We'd like to show tic marks to represent the start of each
   unit's time series, but we have to watch out for the case where
   there are "too many" units -- we don't want a dense fudge of marks
   on the x-axis.  In that case we put a tic mark only for every k'th
   unit.
*/

static void make_panel_unit_tics (const DATASET *dset,
				  gnuplot_info *gi,
				  PRN *prn)
{
    int maxtics, ticskip;
    double ntics;
    int printed;
    int u, t, n;

    pputs(prn, "set xtics (");

    gretl_push_c_numeric_locale();

    /* how many panel units are included in the plot? */
    maxtics = gi->t2 / dset->pd - gi->t1 / dset->pd + 1;

    ntics = maxtics;
    while (ntics > 20) {
	ntics /= 1.5;
    }

    ticskip = maxtics / ceil(ntics);

    if (ticskip == 1 && ntics < maxtics) {
	/* otherwise we'll get an incomplete scale */
	ntics = maxtics;
    }

    n = printed = 0;
    u = gi->t1 / dset->pd;

    for (t=gi->t1; t<=gi->t2 && printed<ntics; t++) {
	if (t == gi->t1 || !continues_unit(dset, t)) {
	    u++;
	    if (n % ticskip == 0) {
		pprintf(prn, "\"%d\" %.10g", u, gi->x[t]);
		if (++printed < ntics) {
		    pputs(prn, ", ");
		}
	    }
	    n++;
	}
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

static void make_calendar_tics (const DATASET *dset,
				const gnuplot_info *gi,
				PRN *prn)
{
    int T = gi->t2 - gi->t1 + 1;
    double yrs;

    if (dset->pd == 52) {
	yrs = T / 52.0;
    } else {
	yrs = T / (dset->pd * 52.0);
    }

    if (yrs <= 3) {
	make_named_month_tics(gi, yrs, prn);
    } else if (yrs < 6) {
	/* don't show ugly "fractions of years" */
	pputs(prn, "set xtics 1\n");
	if (yrs < 3) {
	    /* put monthly minor tics */
	    pputs(prn, "set mxtics 12\n");
	} else if (yrs < 5) {
	    /* quarterly minor tics */
	    pputs(prn, "set mxtics 4\n");
	}
    }
}

static int multiple_groups (const DATASET *dset, int t1, int t2)
{
    int ret = 0;

    if (dataset_is_panel(dset)) {
	ret = (t2 / dset->pd > t1 / dset->pd);
    }

    return ret;
}

static int single_year_sample (const DATASET *dset,
			       int t1, int t2)
{
    char obs[OBSLEN];
    int y1, y2;

    ntolabel(obs, t1, dset);
    y1 = atoi(obs);
    ntolabel(obs, t2, dset);
    y2 = atoi(obs);

    return y2 == y1;
}

/* Use year-and-month major tics when plotting fewer than ??
   monthly observations.
*/

static void short_monthly_tics (gnuplot_info *gi, int T,
				const DATASET *dset,
				PRN *prn)
{
    GDateTime *dt0 = NULL;
    GDateTime *dt1 = NULL;
    gchar *tstr = NULL;
    int use_names = 1;
    int ticskip = 1;
    int y, m, t, j;
    double x;

    if (T > 12) {
	ticskip += (T > 12) + (T > 24);
    }

    date_maj_min(gi->t1, dset, &y, &m);
    if (use_names) {
	dt0 = g_date_time_new_local(y, m, 1, 0, 0, 0);
    }

    pprintf(prn, "# xtics lines = %d\n", T);
    pprintf(prn, "set xtics (");
    gretl_push_c_numeric_locale();

    for (t=gi->t1, j=0; t<=gi->t2; t++, j++) {
	x = y + (m - 1) / 12.0;
	if (j % ticskip == 0) {
	    /* major (labeled) tic */
	    if (use_names) {
		tstr = g_date_time_format(dt0, "%b %Y");
		pprintf(prn, "\"%s\" %g 0", tstr, x);
		dt1 = g_date_time_add_months(dt0, ticskip);
		g_date_time_unref(dt0);
		dt0 = dt1;
		g_free(tstr);
	    } else {
		pprintf(prn, "\"%d-%02d\" %g 0", y, m, x);
	    }
	} else {
	    /* minor (unlabeled) tic */
	    pprintf(prn, "\"\" %g 1", x);
	}
	pputs(prn, t < gi->t2 ? ", \\\n" : ")\n");
	if (m == 12) {
	    y++;
	    m = 1;
	} else {
	    m++;
	}
    }

    gretl_pop_c_numeric_locale();
    if (T > 8) {
	pputs(prn, "set xtics rotate by -45\n");
    }
    if (use_names) {
	g_date_time_unref(dt0);
    }
}

/* special tics for time series plots */

void make_time_tics (gnuplot_info *gi,
		     const DATASET *dset,
		     int many, char *xlabel,
		     PRN *prn)
{
    int T = gi->t2 - gi->t1 + 1;
    int few_years = 8;

    if (many) {
	pprintf(prn, "# multiple timeseries %d\n", dset->pd);
    } else {
	pprintf(prn, "# timeseries %d", dset->pd);
	gi->flags |= GPT_LETTERBOX;
	pputs(prn, " (letterbox)\n");
    }

    if (gi->flags & GPT_TIMEFMT) {
        double yr100 = -59011441438;

	pputs(prn, "set xdata time\n");
	pputs(prn, "set timefmt \"%s\"\n");
	if (single_year_sample(dset, gi->t1, gi->t2)) {
	    strcpy(gi->xfmt, "%m-%d");
	} else if (gi->x != NULL && gi->x[gi->t2] < yr100) {
            /* century not known */
            strcpy(gi->xfmt, "%y-%m-%d");
        } else {
	    strcpy(gi->xfmt, "%Y-%m-%d");
	}
	pprintf(prn, "set format x \"%s\"\n", gi->xfmt);
	pputs(prn, "set xtics rotate by -45\n");
	pputs(prn, "set xtics nomirror\n");
	return;
    }

    if (dset->pd == 1 && T < few_years) {
	pputs(prn, "set xtics nomirror 1\n");
    } else if (dset->pd == 4 && T / 4 < few_years) {
	pputs(prn, "set xtics nomirror 0,1\n");
	pputs(prn, "set mxtics 4\n");
    } else if (dset->pd == 12) {
	if (T < 36) {
	    short_monthly_tics(gi, T, dset, prn);
	} else if (T / 12 < few_years) {
	    pputs(prn, "set xtics nomirror 0,1\n");
	    pputs(prn, "set mxtics 12\n");
	}
    } else if (dated_daily_data(dset) || dated_weekly_data(dset)) {
	make_calendar_tics(dset, gi, prn);
    } else if (multiple_groups(dset, gi->t1, gi->t2)) {
	make_panel_unit_tics(dset, gi, prn);
	if (xlabel != NULL) {
	    strcpy(xlabel, _("time series by group"));
	}
    }
}

/* Handle the use of a matrix in the context of the "plot" command
   block, or create a plot directly from a matrix and a plot list.
*/

int matrix_plot (gretl_matrix *m, const int *list, const char *literal,
		 gretlopt opt)
{
    DATASET *dset = NULL;
    int *plotlist = NULL;
    int pmax, err = 0;

    if (gretl_is_null_matrix(m)) {
	err = E_DATA;
    } else if (list != NULL && list[0] == 0) {
	/* let an empty list be equivalent to NULL */
	dset = gretl_dataset_from_matrix(m, NULL, OPT_B, &err);
    } else {
	dset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
    }
    if (err) {
	return err;
    }

    pmax = dset->v - 1;
    if (pmax <= 0) {
	err = E_DATA;
    } else {
	plotlist = gretl_consecutive_list_new(1, pmax);
	if (plotlist == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	if (opt & (OPT_N | OPT_a)) {
	    /* --band or --bands */
	    gnuplot_info gi;

	    err = gpinfo_init(&gi, opt, plotlist, literal, dset);
	    if (!err) {
		err = maybe_add_plotx(&gi, 0, dset);
	    }
	    if (!err) {
		err = plot_with_band(BP_BLOCKMAT, &gi, literal,
				     dset, opt);
	    }
	} else {
	    err = gnuplot(plotlist, literal, dset, opt);
	}
    }

    destroy_dataset(dset);
    free(plotlist);

    return err;
}

static int plotlist_is_group_invariant (const int *list, const DATASET *dset)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (!series_is_group_invariant(dset, list[i])) {
	    return 0;
	}
    }

    return 1;
}

static int panel_group_invariant_plot (const int *plotlist,
				       const char *literal,
				       DATASET *dset,
				       gretlopt opt)
{
    DATASET orig = *dset;
    int err;

    /* limit sample to first group */
    dset->t1 = 0;
    dset->t2 = dset->pd - 1;

    /* and mark as time-series data */
    if (dset->panel_pd > 0) {
	dset->pd = dset->panel_pd;
	dset->sd0 = dset->panel_sd0;
	dset->structure = TIME_SERIES;
    } else {
	dset->structure = SPECIAL_TIME_SERIES;
	dset->pd = 1;
    }

    err = gnuplot(plotlist, literal, dset, opt);

    /* put everything back as it was */
    *dset = orig;

    return err;
}

static int time_fit_wanted (gretlopt *popt)
{
    if ((*popt & OPT_T) && (*popt & OPT_F)) {
        const char *s = get_optval_string(GNUPLOT, OPT_F);

	if (s == NULL || !strcmp(s, "none")) {
	    *popt &= ~OPT_F;
	} else {
	    return 1;
	}
    }

    return 0;
}

/**
 * gnuplot:
 * @plotlist: list of variables to plot, by ID number.
 * @literal: commands to be passed to gnuplot.
 * @dset: dataset struct.
 * @opt: option flags.
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gnuplot (const int *plotlist, const char *literal,
	     const DATASET *dset, gretlopt opt)
{
    PRN *prn = NULL;
    FILE *fp = NULL;
    int *list = NULL;
    char s1[MAXDISP] = {0};
    char s2[MAXDISP] = {0};
    char xlabel[MAXDISP] = {0};
    char withstr[16] = {0};
    char lwstr[8] = {0};
    char keystr[48] = {0};
    gchar *fitline = NULL;
    int time_fit = 0;
    int oddman = 0;
    int many = 0;
    int set_xrange = 1;
    PlotType ptype;
    gnuplot_info gi = {0};
    int i, err = 0;

    gretl_error_clear();

    if (time_fit_wanted(&opt)) {
	if (plotlist[0] > 1 || !dataset_is_time_series(dset)) {
	    return E_BADOPT;
	} else {
	    time_fit = 1;
	}
    }

    if (literal != NULL && strstr(literal, "set xdata time")) {
	set_xrange = 0;
    }

    if (dataset_is_panel(dset) &&
	plotlist_is_group_invariant(plotlist, dset)) {
#if GP_DEBUG
	fprintf(stderr, "doing panel_group_invariant_plot\n");
#endif
	return panel_group_invariant_plot(plotlist, literal,
					  (DATASET *) dset, opt);
    }

#if GP_DEBUG
    printlist(plotlist, "gnuplot: plotlist");
    fprintf(stderr, "incoming plot range: obs %d to %d\n", dset->t1, dset->t2);
#endif

    err = gpinfo_init(&gi, opt, plotlist, literal, dset);
    if (err) {
	goto bailout;
    }

#if GP_DEBUG
    fprintf(stderr, "after gpinfo_init: gi.fit = %d\n", gi.fit);
#endif

    err = maybe_add_plotx(&gi, time_fit, dset);
    if (err) {
	goto bailout;
    }

    /* hive off some special cases */
    if (gi.fit == PLOT_FIT_LOESS) {
	return loess_plot(&gi, literal, dset);
    } else if (time_fit) {
	return time_fit_plot(&gi, literal, dset);
    } else if (gi.band) {
	return plot_with_band(BP_REGULAR, &gi, literal,
			      (DATASET *) dset, opt);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
	many = 1;
    }

    /* convenience pointer */
    list = gi.list;

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GPT_TS)) {
	int v = (gi.flags & GPT_DUMMY)? list[2] : list[list[0]];

	strcpy(xlabel, series_get_graph_name(dset, v));
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	goto bailout;
    }

    /* adjust sample range, and reject if it's empty */
    err = graph_list_adjust_sample(list, &gi, dset, 2);
    if (err) {
	goto bailout;
    }

    /* add a regression line if appropriate */
    if (!use_impulses(&gi) && !(gi.flags & GPT_FIT_OMIT) && list[0] == 2 &&
	!(gi.flags & GPT_TS) && !(gi.flags & GPT_RESIDS)) {
	err = get_fitted_line(&gi, dset, &fitline);
	if (err) {
	    goto bailout;
	} else {
	    const char *xname = dset->varname[list[2]];
	    const char *yname = dset->varname[list[1]];

	    if (*xname != '\0' && *yname != '\0') {
		pprintf(prn, "# X = '%s' (%d)\n", xname, list[2]);
		pprintf(prn, "# Y = '%s' (%d)\n", yname, list[1]);
	    }
	}
    }

    ptype = PLOT_REGULAR;

    /* separation by dummy: create special vars */
    if (gi.flags & GPT_DUMMY) {
	err = factor_check(&gi, dset);
	if (err) {
	    goto bailout;
	}
	ptype = PLOT_FACTORIZED;
    }

    /* special tics for time series plots */
    if (gi.flags & GPT_TS) {
	make_time_tics(&gi, dset, many, xlabel, prn);
    }

    /* open file and, if that goes OK, dump the prn into it
       after writing the header
    */
    fp = open_plot_input_file(ptype, gi.flags, &err);
    if (err) {
	gretl_print_destroy(prn);
	goto bailout;
    }

    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    print_axis_label('x', xlabel, fp);
    fputs("set xzeroaxis\n", fp);
    gnuplot_missval_string(fp);

    if (gi.flags & GPT_LOGY) {
	fprintf(fp, "set logscale y %g\n", gi.ybase);
    }

    /* key: default to left top */
    strcpy(keystr, "set key left top\n");

    if (list[0] == 1) {
	/* only one variable (time series) */
	print_axis_label('y', series_get_graph_name(dset, list[1]), fp);
	strcpy(keystr, "set nokey\n");
    } else if (list[0] == 2) {
	/* plotting two variables */
	int no_key = 1;

	if (gi.flags & GPT_AUTO_FIT) {
	    print_auto_fit_string(gi.fit, fp);
	    if (gi.flags & GPT_FA) {
		make_gtitle(&gi, GTITLE_AFV, series_get_graph_name(dset, list[1]),
			    series_get_graph_name(dset, list[2]), fp);
	    } else {
		make_gtitle(&gi, GTITLE_VLS, series_get_graph_name(dset, list[1]),
			    xlabel, fp);
	    }
	    no_key = 0;
	}
	if (gi.flags & GPT_RESIDS) {
	    const char *vlabel = series_get_label(dset, list[1]);

	    make_gtitle(&gi, GTITLE_RESID, vlabel == NULL ? "residual" : vlabel,
			NULL, fp);
	    fprintf(fp, "set ylabel '%s'\n", _("residual"));
	} else {
	    print_axis_label('y', series_get_graph_name(dset, list[1]), fp);
	}
	if (no_key) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((gi.flags & GPT_RESIDS) && (gi.flags & GPT_DUMMY)) {
	const char *vlabel = series_get_label(dset, list[1]);

	make_gtitle(&gi, GTITLE_RESID, vlabel == NULL ? "residual" : vlabel,
		    NULL, fp);
	fprintf(fp, "set ylabel '%s'\n", _("residual"));
    } else if (gi.flags & GPT_FA) {
	if (list[3] == dset->v - 1) {
	    /* x var is just time or index: is this always right? */
	    make_gtitle(&gi, GTITLE_AF, series_get_graph_name(dset, list[2]),
			NULL, fp);
	} else {
	    make_gtitle(&gi, GTITLE_AFV, series_get_graph_name(dset, list[2]),
			series_get_graph_name(dset, list[3]), fp);
	}
	print_axis_label('y', series_get_graph_name(dset, list[2]), fp);
    } else if (gi.flags & GPT_DUMMY) {
	print_axis_label('y', series_get_graph_name(dset, list[1]), fp);
    }

    if (many) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    if (set_xrange) {
	if (gi.x != NULL) {
	    print_x_range(&gi, fp);
	} else {
	    print_x_range_from_list(&gi, dset, list, fp);
	}
    }

    if (!(gi.flags & GPT_TIMEFMT) && *gi.xfmt != '\0' && *gi.xtics != '\0') {
	/* remedial handling of broken tics */
	fprintf(fp, "set format x \"%s\"\n", gi.xfmt);
	fprintf(fp, "set xtics %s\n", gi.xtics);
    }

    if (gi.flags & GPT_Y2AXIS) {
	check_for_yscale(&gi, (const double **) dset->Z, &oddman);
	if (gi.flags & GPT_Y2AXIS) {
	    fputs("set ytics nomirror\n", fp);
	    fputs("set y2tics\n", fp);
	}
    } else if (gi.yformula == NULL && list[0] == 2) {
	check_y_tics(&gi, (const double **) dset->Z, fp);
    }

#if GP_DEBUG
    fprintf(stderr, "literal = '%s', yformula = '%s'\n", literal,
	    gi.yformula);
#endif

    if (gi.yformula != NULL) {
	/* cut out the "dummy" yvar that is in fact represented
	   by a formula rather than raw data */
	list[1] = list[2];
	list[2] = list[3];
	list[0] = 2;
    } else {
	print_gnuplot_literal_lines(literal, GNUPLOT, opt, fp);
    }

    /* now print the 'plot' lines */
    fputs("plot \\\n", fp);
    if (gi.flags & GPT_Y2AXIS) {
	/* using two y axes */
	int lmax = list[0];

	if ((gi.flags & GPT_IDX) && list[lmax] != 0) {
	    lmax++;
	}
	for (i=1; i<lmax; i++) {
	    set_lwstr(lwstr);
	    set_plot_withstr(&gi, i, withstr);
	    fprintf(fp, " '-' using 1:2 axes %s title \"%s (%s)\" %s%s%s",
		    (i == oddman)? "x1y2" : "x1y1",
		    series_get_graph_name(dset, list[i]),
		    (i == oddman)? _("right") : _("left"),
		    withstr,
		    lwstr,
		    (i == lmax - 1)? "\n" : ", \\\n");
	}
    } else if (gi.flags & GPT_DUMMY) {
	/* plot shows separation by discrete variable */
	int nd = gretl_vector_get_length(gi.dvals);
	int dv = list[3];
	series_table *st;

	strcpy(s2, series_get_graph_name(dset, dv));
	st = series_get_string_table(dset, dv);

	for (i=0; i<nd; i++) {
	    double di = gretl_vector_get(gi.dvals, i);

	    if (st != NULL) {
		fprintf(fp, " '-' using 1:2 title \"%s\" w points",
			series_table_get_string(st, di));
	    } else {
		fprintf(fp, " '-' using 1:2 title \"%s=%g\" w points",
			s2, di);
	    }
	    if (i < nd - 1) {
		fputs(", \\\n", fp);
	    } else {
		fputc('\n', fp);
	    }
	}
    } else if (gi.yformula != NULL) {
	/* we have a formula to plot, not just data */
	fprintf(fp, " '-' using 1:2 title \"%s\" w points, \\\n", _("actual"));
	fprintf(fp, "%s title '%s' w lines\n", gi.yformula, _("fitted"));
    } else if (gi.flags & GPT_FA) {
	/* this is a fitted vs actual plot */
	/* try reversing here: 2014-09-22 */
	int tmp = list[1];

	list[1] = list[2];
	list[2] = tmp;
	set_plot_withstr(&gi, 1, withstr);
	fprintf(fp, " '-' using 1:2 title \"%s\" %s, \\\n", _("actual"), withstr);
	fprintf(fp, " '-' using 1:2 title \"%s\" %s\n", _("fitted"), withstr);
    } else {
	/* all other cases */
	int lmax = list[0] - 1;

	for (i=1; i<=lmax; i++)  {
	    set_lwstr(lwstr);
	    if (list[0] == 2 && !(gi.flags & GPT_TIMEFMT)) {
		*s1 = '\0';
	    } else {
		strcpy(s1, series_get_graph_name(dset, list[i]));
	    }
	    set_plot_withstr(&gi, i, withstr);
	    fprintf(fp, " '-' using 1:2 title \"%s\" %s%s", s1, withstr, lwstr);
	    if (i < lmax || (gi.flags & GPT_AUTO_FIT)) {
	        fputs(", \\\n", fp);
	    } else {
	        fputc('\n', fp);
	    }
	}
    }

    if (fitline != NULL) {
        fputs(fitline, fp);
    }

    /* print the data to be graphed */
    if (gi.flags & GPT_DUMMY) {
	print_gp_dummy_data(&gi, dset, fp);
    } else {
	print_gp_data(&gi, dset, fp);
    }

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    g_free(fitline);
    clear_gpinfo(&gi);

    return err;
}

int theil_forecast_plot (const int *plotlist, const DATASET *dset,
			 gretlopt opt)
{
    FILE *fp = NULL;
    gnuplot_info gi;
    int vx, vy;
    int err = 0;

    gretl_error_clear();

    if (plotlist[0] != 2) {
	return E_DATA;
    }

    err = gpinfo_init(&gi, opt | OPT_S, plotlist, NULL, dset);
    if (err) {
	goto bailout;
    }

    /* ensure the time-series flag is unset */
    gi.flags &= ~GPT_TS;

    err = graph_list_adjust_sample(gi.list, &gi, dset, 1);
    if (err) {
	goto bailout;
    }

    fp = open_plot_input_file(PLOT_REGULAR, gi.flags, &err);
    if (err) {
	goto bailout;
    }

    vx = gi.list[2];
    vy = gi.list[1];

    print_axis_label('x', series_get_graph_name(dset, vx), fp);
    print_axis_label('y', series_get_graph_name(dset, vy), fp);

    fputs("set xzeroaxis\n", fp);
    gnuplot_missval_string(fp);
    fputs("set key left top\n", fp);

    gretl_push_c_numeric_locale();

    print_x_range_from_list(&gi, dset, gi.list, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    fprintf(fp, " x title \"%s\" w lines\n", _("actual = predicted"));

    print_gp_data(&gi, dset, fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    clear_gpinfo(&gi);

    return err;
}

static void multi_time_tics (const double *obs,
			     const DATASET *dset,
			     FILE *fp)
{
    double startdate = obs[dset->t1];
    double enddate = obs[dset->t2];
    double obsrange = enddate - startdate;
    int k1 = ceil(startdate);
    int k2 = floor(enddate);
    int nmaj = k2 - k1 + 1;

    fputs("set xtics nomirror\n", fp);

    if (obsrange > 8) {
	double incr = obsrange / 4;

	fprintf(fp, "set xrange [%g:%g]\n", floor(startdate), ceil(enddate));
	fprintf(fp, "set xtics %g,%g\n", ceil(startdate), floor(incr));
    } else if (nmaj == 0) {
	fputs("set format x ''\n", fp);
    } else {
	/* integer major tics plus minor */
	int T = dset->t2 - dset->t1 + 1;

	fprintf(fp, "set xrange [%g:%g]\n", startdate, enddate);
	fprintf(fp, "set xtics %g,1\n", floor(startdate));
	if (T < 55) {
	    fprintf(fp, "set mxtics %d\n", dset->pd);
	}
    }
}

static void multi_set_timefmt (const DATASET *dset,
			       const double *obs,
			       FILE *fp)
{
    double T = obs[dset->t2] - obs[dset->t1];
    int ntics = 6;

    fputs("set xdata time\n", fp);
    fputs("set timefmt \"%s\"\n", fp);
    if (single_year_sample(dset, dset->t1, dset->t2)) {
	fputs("set format x \"%m-%d\"\n", fp);
    } else {
	fputs("set format x \"%Y-%m-%d\"\n", fp);
    }
    fputs("set xtics rotate by -45\n", fp);
    fprintf(fp, "set xtics %g\n", round(T/ntics));
}

/**
 * multi_plots:
 * @list: list of variables to plot, by ID number.
 * @dset: dataset struct.
 * @opt: can include %OPT_O to use lines, %PT_T to
 * produce time-series plots, %OPT_U to direct output
 * to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the variables in @list, and calls gnuplot to make
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_plots (const int *list, const DATASET *dset,
		 gretlopt opt)
{
    GptFlags flags = 0;
    int xvar = 0, yvar = 0;
    const double *x = NULL;
    const double *y = NULL;
    const double *obs = NULL;
    int rows, cols, tseries = 0;
    int use_timefmt = 0;
    int tsdata = 0;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

    pos = gretl_list_separator_position(list);
    tsdata = dataset_is_time_series(dset);

    if (opt & OPT_O) {
	flags |= GPT_LINES;
    }

    if (opt & OPT_T) {
	/* the "tsplots" command */
	if (pos > 0) {
	    /* separator not accepted */
	    err = E_INVARG;
	} else if (!tsdata) {
	    err = E_PDWRONG;
	} else {
	    flags |= GPT_LINES;
	}
    }

    if (pos == 0) {
	/* plot against time or index */
	plotlist = gretl_list_copy(list);
	flags |= GPT_LINES;
	if (tsdata) {
	    tseries = 1;
	    if (calendar_data(dset)) {
		use_timefmt = 1;
	    }
	}
	obs = gretl_plotx(dset, use_timefmt ? OPT_T : OPT_S);
	if (obs == NULL) {
	    return E_ALLOC;
	}
    } else if (pos > 2) {
	/* plot several yvars against one xvar */
	plotlist = gretl_list_new(pos - 1);
	xvar = list[list[0]];
	x = dset->Z[xvar];
    } else {
	/* plot one yvar against several xvars */
	plotlist = gretl_list_new(list[0] - pos);
	yvar = list[1];
	y = dset->Z[yvar];
    }

    if (plotlist == NULL) {
	return E_ALLOC;
    }

    if (yvar) {
	for (i=1; i<=plotlist[0]; i++) {
	    plotlist[i] = list[i + pos];
	}
    } else if (xvar) {
	for (i=1; i<pos; i++) {
	    plotlist[i] = list[i];
	}
    }

    /* max 16 plots */
    if (plotlist[0] > 16) {
	plotlist[0] = 16;
    }

    nplots = plotlist[0];

    if (nplots > 1) {
	get_multiplot_layout(nplots, tseries, &rows, &cols);
	if (use_timefmt) {
	    gp_small_font_size = nplots > 2 ? 6 : 0;
	} else {
	    maybe_set_small_font(nplots);
	}
	if (nplots > 12) {
	    flags |= GPT_XXL;
	} else if (nplots > 9) {
	    flags |= GPT_XL;
	}
    }

    fp = open_plot_input_file(PLOT_MULTI_BASIC, flags, &err);
    if (err) {
	return err;
    }

    if (nplots > 1) {
	fprintf(fp, "set multiplot layout %d,%d\n", rows, cols);
    }
    fputs("set nokey\n", fp);

    if (opt & OPT_K) {
	/* --tweaks=foo */
	print_gnuplot_literal_lines(NULL, SCATTERS, opt, fp);
    }

    gretl_push_c_numeric_locale();

    if (use_timefmt) {
	fprintf(fp, "set xrange [%.0f:%.0f]\n", obs[dset->t1], obs[dset->t2]);
	multi_set_timefmt(dset, obs, fp);
    } else if (obs != NULL) {
	multi_time_tics(obs, dset, fp);
    } else {
	/* avoid having points sticking to the axes */
	fputs("set offsets graph 0.02, graph 0.02, graph 0.02, graph 0.02\n", fp);
	fputs("set noxtics\nset noytics\n", fp);
    }

    fputs("set xzeroaxis\n", fp);

    for (i=0; i<nplots; i++) {
	int j = plotlist[i+1];

	if (obs != NULL) {
	    fputs("set noxlabel\n", fp);
	    fputs("set noylabel\n", fp);
	    fprintf(fp, "set title '%s'\n", series_get_graph_name(dset, j));
	} else {
	    fprintf(fp, "set xlabel '%s'\n",
		    (yvar)? dset->varname[j] :
		    dset->varname[xvar]);
	    fprintf(fp, "set ylabel '%s'\n",
		    (yvar)? dset->varname[yvar] :
		    dset->varname[j]);
	}
	fputs("plot '-' using 1:2", fp);
	if (flags & GPT_LINES) {
	    fputs(" with lines", fp);
	}
	fputc('\n', fp);

	for (t=dset->t1; t<=dset->t2; t++) {
	    double xt = yvar ? dset->Z[j][t] : xvar ? x[t] : obs[t];
	    double yt = yvar ? y[t] : dset->Z[j][t];

	    write_gp_dataval(xt, fp, 0);
	    write_gp_dataval(yt, fp, 1);
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    if (nplots > 1) {
	fputs("unset multiplot\n", fp);
    }

    free(plotlist);

    return finalize_plot_input_file(fp);
}

static int matrix_plotx_ok (const gretl_matrix *m, const DATASET *dset,
			    int *pt1, int *pt2, int *ppd)
{
    if (dset == NULL) {
	return 0;
    } else if (m->rows == dset->n) {
	return 1;
    } else {
	int t1 = gretl_matrix_get_t1(m);
	int t2 = gretl_matrix_get_t2(m);

	if (t2 > t1 && t2 < dset->n) {
	    *pt1 = t1;
	    *pt2 = t2;
	    *ppd = dset->pd;
	    return 1;
	}
    }

    return 0;
}

static const double *matrix_col (const gretl_matrix *m, int j)
{
    const double *x = m->val;

    return x + (j-1) * m->rows;
}

static void plot_colname (char *s, const char **colnames, int j)
{
    if (colnames != NULL) {
	const char *name = colnames[j-1];

	*s = '\0';
	if (strlen(name) >= 16) {
	    strncat(s, name, 14);
	    strcat(s, "~");
	} else {
	    strncat(s, name, 15);
	}
    } else {
	sprintf(s, "col %d", j);
    }
}

static double get_obsx (const double *obs, int t, int s)
{
    return (obs != NULL)? obs[t] : s;
}

/**
 * matrix_multi_plots:
 * @m: matrix containing data to plot.
 * @list: list of columns to plot, or NULL.
 * @dset: dataset pointer, or NULL.
 * @opt: can include %OPT_O to use lines, %OPT_U to
 * direct output to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the data in @m, and calls gnuplot to make
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int matrix_multi_plots (const gretl_matrix *m, const int *list,
			const DATASET *dset, gretlopt opt)
{
    GptFlags flags = 0;
    const double *x = NULL;
    const double *y = NULL;
    const double *obs = NULL;
    const char **colnames = NULL;
    FILE *fp = NULL;
    int *plotlist = NULL;
    int need_list = 0;
    int rows, cols;
    int xcol = 0, ycol = 0;
    int t1 = 0, t2 = 0, pd = 1;
    int pos = 0, nplots = 0;
    int simple_obs = 0;
    int i, t, err = 0;

    if (gretl_is_null_matrix(m)) {
	return E_DATA;
    }

    if (opt & OPT_O) {
	flags |= GPT_LINES;
    }

    if (list != NULL) {
	for (i=1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) {
		pos = i;
	    } else if (list[i] < 1 || list[i] > m->cols) {
		err = E_INVARG;
		break;
	    }
	}
    }

    if (err) {
	return err;
    }

    t1 = 0;
    t2 = m->rows - 1;

    if (pos == 0) {
	/* plot against time or index */
	if (matrix_plotx_ok(m, dset, &t1, &t2, &pd)) {
	    obs = gretl_plotx(dset, OPT_NONE);
	    if (obs == NULL) {
		return E_ALLOC;
	    }
	} else {
	    simple_obs = 1;
	}
	if (list != NULL && list[0] > 0) {
	    need_list = 1;
	    plotlist = gretl_list_copy(list);
	}
	flags |= GPT_LINES;
    } else if (pos > 2) {
	/* plot several yvars against one xvar */
	need_list = 1;
	plotlist = gretl_list_new(pos - 1);
	xcol = list[list[0]];
	x = matrix_col(m, xcol);
    } else {
	/* plot one yvar against several xvars */
	need_list = 1;
	plotlist = gretl_list_new(list[0] - pos);
	ycol = list[1];
	y = matrix_col(m, ycol);
    }

    if (need_list && plotlist == NULL) {
	return E_ALLOC;
    }

    if (plotlist != NULL) {
	if (y != NULL) {
	    for (i=1; i<=plotlist[0]; i++) {
		plotlist[i] = list[i + pos];
	    }
	} else if (x != NULL) {
	    for (i=1; i<pos; i++) {
		plotlist[i] = list[i];
	    }
	}
	/* max 16 plots */
	if (plotlist[0] > 16) {
	    plotlist[0] = 16;
	}
	nplots = plotlist[0];
    } else {
	nplots = (m->cols > 16)? 16 : m->cols;
    }

    get_multiplot_layout(nplots, 0, &rows, &cols);
    maybe_set_small_font(nplots);

    if (nplots > 12) {
	flags |= GPT_XXL;
    } else if (nplots > 9) {
	flags |= GPT_XL;
    }

    fp = open_plot_input_file(PLOT_MULTI_BASIC, flags, &err);
    if (err) {
	return err;
    }

    colnames = gretl_matrix_get_colnames(m);

    fprintf(fp, "set multiplot layout %d,%d\n", rows, cols);
    fputs("set xzeroaxis\n", fp);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (obs != NULL) {
	double startdate = obs[t1];
	double enddate = obs[t2];
	int incr, T = t2 - t1 + 1;

	fprintf(fp, "set xrange [%g:%g]\n", floor(startdate), ceil(enddate));

	incr = (pd == 1)? (T / 6) : (T / (4 * pd));
	if (incr > 0) {
	    fprintf(fp, "set xtics %g, %d\n", ceil(startdate), incr);
	}
    } else if (simple_obs) {
	int incr = m->rows / 6;

	fprintf(fp, "set xrange [0:%d]\n", m->rows - 1);
	if (incr > 0) {
	    fprintf(fp, "set xtics 0, %d\n", incr);
	}
    } else {
	fputs("set noxtics\nset noytics\n", fp);
    }

    if (obs != NULL || simple_obs) {
	fputs("set noxlabel\n", fp);
	fputs("set noylabel\n", fp);
    }

    for (i=0; i<nplots; i++) {
	int j = (plotlist == NULL)? (i+1) : plotlist[i+1];
	const double *zj = matrix_col(m, j);
	char label[16];

	if (obs != NULL || simple_obs) {
	    plot_colname(label, colnames, j);
	    fprintf(fp, "set title '%s'\n", label);
	} else {
	    plot_colname(label, colnames, (y != NULL)? j : xcol);
	    fprintf(fp, "set xlabel '%s'\n", label);
	    plot_colname(label, colnames, (y != NULL)? ycol : j);
	    fprintf(fp, "set ylabel '%s'\n", label);
	}

	fputs("plot '-' using 1:2", fp);
	if (flags & GPT_LINES) {
	    fputs(" with lines", fp);
	}
	fputc('\n', fp);

	for (t=t1; t<=t2; t++) {
	    int s = t - t1;
	    double xt = ycol ? zj[s] : xcol ? x[s] : get_obsx(obs, t, s);
	    double yt = (y != NULL)? y[s] : zj[s];

	    write_gp_dataval(xt, fp, 0);
	    write_gp_dataval(yt, fp, 1);
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

    free(plotlist);

    return finalize_plot_input_file(fp);
}

static FILE *get_3d_input_file (int *err)
{
    FILE *fp = NULL;
    char fname[MAXLEN];

    sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
    fp = gretl_fopen(fname, "w");

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	gretl_set_path_by_name("plotfile", fname);
    }

    return fp;
}

static gchar *maybe_get_surface (const int *list,
				 DATASET *dset,
				 gretlopt opt)
{
    MODEL smod;
    double umin, umax, vmin, vmax;
    int olslist[5];
    gchar *ret = NULL;

    olslist[0] = 4;
    olslist[1] = list[3];
    olslist[2] = 0;
    olslist[3] = list[2];
    olslist[4] = list[1];

    gretl_minmax(dset->t1, dset->t2, dset->Z[list[2]], &umin, &umax);
    gretl_minmax(dset->t1, dset->t2, dset->Z[list[1]], &vmin, &vmax);

    smod = lsq(olslist, dset, OLS, OPT_A);

    if (!smod.errcode && !na(smod.fstt) &&
	(snedecor_cdf_comp(smod.dfn, smod.dfd, smod.fstt) < .10 || (opt & OPT_A))) {
	double uadj = (umax - umin) * 0.02;
	double vadj = (vmax - vmin) * 0.02;

	ret = g_strdup_printf("[u=%g:%g] [v=%g:%g] "
			      "%g+(%g)*u+(%g)*v notitle",
			      umin - uadj, umax + uadj,
			      vmin - vadj, vmax + vadj,
			      smod.coeff[0], smod.coeff[1],
			      smod.coeff[2]);
    }

    clear_model(&smod);

    return ret;
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @dset: pointer to dataset.
 * @opt: may include OPT_A to force display of fitted surface;
 * may include OPT_I to force an interactive (rotatable) plot.
 * Note that OPT_I may be removed on output if a suitable
 * gnuplot terminal is not present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		DATASET *dset, gretlopt *opt)
{
    FILE *fp = NULL;
    int t, t1 = dset->t1, t2 = dset->t2;
    int save_t1 = dset->t1, save_t2 = dset->t2;
    int lo = list[0];
    int datlist[4];
    int interactive = (*opt & OPT_I);
    const char *term = NULL;
    gchar *surface = NULL;
    int err = 0;

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return E_DATA;
    }

    list_adjust_sample(list, &t1, &t2, dset, NULL);

    /* if resulting sample range is empty, complain */
    if (t1 >= t2) {
	return E_MISSDATA;
    }

#ifndef WIN32
    if (interactive) {
	/* On Windows we let the gnuplot terminal default to
	   "win"; on other systems we need a suitable
	   terminal for interactive 3-D display.
	*/
	if (gnuplot_has_wxt()) {
	    term = "wxt size 640,420 noenhanced";
	} else if (gnuplot_has_x11()) {
	    term = "x11";
	} else if (gnuplot_has_qt()) {
	    term = "qt";
	} else {
	    *opt &= ~OPT_I;
	    interactive = 0;
	}
    }
#endif

    if (interactive) {
	fp = get_3d_input_file(&err);
    } else {
	fp = open_plot_input_file(PLOT_3D, 0, &err);
    }

    if (err) {
	return err;
    }

    dset->t1 = t1;
    dset->t2 = t2;

    if (interactive) {
	if (term != NULL) {
	    fprintf(fp, "set term %s\n", term);
	}
	write_plot_line_styles(PLOT_3D, fp);
    }

    gretl_push_c_numeric_locale();

    print_axis_label('x', series_get_graph_name(dset, list[2]), fp);
    print_axis_label('y', series_get_graph_name(dset, list[1]), fp);
    print_axis_label('z', series_get_graph_name(dset, list[3]), fp);

    gnuplot_missval_string(fp);

    print_gnuplot_literal_lines(literal, GNUPLOT, *opt, fp);

    surface = maybe_get_surface(list, dset, *opt);

    if (surface != NULL) {
	fprintf(fp, "splot %s, \\\n'-' notitle w p\n", surface);
	g_free(surface);
    } else {
	fputs("splot '-' notitle w p\n", fp);
    }

    datlist[0] = 3;
    datlist[1] = list[2];
    datlist[2] = list[1];
    datlist[3] = list[3];

    for (t=t1; t<=t2; t++) {
	const char *label = NULL;

	if (dset->markers) {
	    label = dset->S[t];
	}
	printvars(fp, t, datlist, dset, NULL, label, 0.0);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    if (interactive) {
	fputs("pause mouse close\n", fp);
	fclose(fp);
    } else {
	err = finalize_plot_input_file(fp);
    }

    return err;
}

/**
 * open_3d_plot_input_file:
 * @iact: on input, non-zero if an interactive plot is
 * preferred, 0 otherwise; on output, non-zero if interactive
 * status can be supported, 0 otherwise.
 *
 * Writes a gnuplot plot file to display a 3D plot
 * (interactive if requested and feasible).
 *
 * Returns: FILE pointer on success, NULL on error.
 */

FILE *open_3d_plot_input_file (int *iact)
{
    const char *term = NULL;
    FILE *fp = NULL;
    int err = 0;

    if (*iact != 0) {
#ifndef WIN32
	/* On Windows we let the gnuplot terminal default to
	   "win"; on other operating systems we need a suitable
	   terminal for interactive 3-D display.
	*/
	if (gnuplot_has_wxt()) {
	    term = "wxt size 640,420 noenhanced";
	} else if (gnuplot_has_x11()) {
	    term = "x11";
	} else if (gnuplot_has_qt()) {
	    term = "qt";
	} else {
	    /* can't do it? */
	    *iact = 0;
	}
#endif
    }

    if (*iact != 0) {
	fp = get_3d_input_file(&err);
    } else {
	fp = open_plot_input_file(PLOT_3D, 0, &err);
    }

    if (*iact) {
	if (term != NULL) {
	    fprintf(fp, "set term %s\n", term);
	}
	write_plot_line_styles(PLOT_3D, fp);
    }

    return fp;
}

static gchar *make_freq_test_label (int teststat, double v, double pv)
{
    gchar *s;

    gretl_pop_c_numeric_locale();
    if (teststat == GRETL_STAT_Z) {
	s = g_strdup_printf("z = %.3f [%.4f]", v, pv);
    } else if (teststat == GRETL_STAT_NORMAL_CHISQ) {
	s = g_strdup_printf("%s(2) = %.3f [%.4f]", _("Chi-square"), v, pv);
    }
    gretl_push_c_numeric_locale();

    return s;
}

static gchar *make_freq_dist_label (int dist, double x, double y)
{
    gchar *s;
    char c, test[10];

    gretl_pop_c_numeric_locale();
    sprintf(test, "%g", 0.5);
    c = strchr(test, ',') ? ' ' : ',';

    if (dist == D_NORMAL) {
	s = g_strdup_printf("N(%.5g%c%.5g)", x, c, y);
    } else if (dist == D_GAMMA) {
	s = g_strdup_printf("gamma(%.5g%c%.5g)", x, c, y);
    }
    gretl_push_c_numeric_locale();

    return s;
}

/* Below: a fix for the case where the y-range is by default
   degenerate, in which case gnuplot produces a graph OK, but
   issues a warning and returns non-zero.
*/

static void maybe_set_yrange (FreqDist *freq, double lambda, FILE *fp)
{
    double ymin = 1.0e+200;
    double ymax = -1.0e+200;
    int i;

    for (i=0; i<freq->numbins; i++) {
	if (freq->f[i] > ymax) {
	    ymax = freq->f[i];
	}
	if (freq->f[i] < ymin) {
	    ymin = freq->f[i];
	}
    }

    if (ymax == ymin) {
	fprintf(fp, "set yrange [%.10g:%.10g]\n", ymax * lambda * 0.99,
		ymax * lambda * 1.01);
    } else {
	fprintf(fp, "set yrange [0.0:%.10g]\n", ymax * lambda * 1.1);
    }
}

static double discrete_minskip (FreqDist *freq)
{
    double s, ms = freq->midpt[1] - freq->midpt[0];
    int i;

    for (i=2; i<freq->numbins; i++) {
	s = freq->midpt[i] - freq->midpt[i-1];
	if (s < ms) {
	    ms = s;
	}
    }

    return ms;
}

/**
 * plot_freq:
 * @freq: pointer to frequency distribution struct.
 * @dist: probability distribution code.
 *
 * Plot the actual frequency distribution for a variable versus a
 * theoretical distribution: Gaussian, gamma or none.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_freq (FreqDist *freq, DistCode dist, gretlopt opt)
{
    double alpha = 0.0, beta = 0.0, lambda = 1.0;
    FILE *fp = NULL;
    int i, K = freq->numbins;
    char withstr[32] = {0};
    gchar *label = NULL;
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth;
    const double *endpt;
    int plottype, use_boxes = 1;
    char **S = NULL;
    int ns = 0, real_ns = 0;
    int err = 0;

    if (K == 0) {
	return E_DATA;
    }

    if (K == 1) {
	gretl_errmsg_sprintf(_("'%s' is a constant"), freq->varname);
	return E_DATA;
    }

    if (freq->strvals) {
	dist = 0; /* just to be safe */
    }

    if (dist == D_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == D_GAMMA) {
	plottype = PLOT_FREQ_GAMMA;
    } else if (freq->discrete) {
	plottype = PLOT_FREQ_DISCRETE;
    } else {
	plottype = PLOT_FREQ_SIMPLE;
    }

    fp = open_plot_input_file(plottype, 0, &err);
    if (err) {
	return err;
    }

#if GP_DEBUG
    fprintf(stderr, "*** plot_freq called\n");
#endif

    if (freq->strvals) {
	endpt = NULL;
	barwidth = 1;
	use_boxes = 0;
    } else if (freq->discrete) {
	endpt = freq->midpt;
	barwidth = discrete_minskip(freq);
	use_boxes = 0;
    } else {
	/* equally sized bins, width to be determined */
	endpt = freq->endpt;
	barwidth = freq->endpt[K-1] - freq->endpt[K-2];
    }

    S = literal_strings_from_opt(FREQ, &ns, &real_ns);

    gretl_push_c_numeric_locale();

    if (dist) {
	int nlit = 2 + 2 * (!na(freq->test)) + real_ns;

	lambda = 1.0 / (freq->n * barwidth);

	if (dist == D_NORMAL) {
	    fprintf(fp, "# literal lines = %d\n", nlit);
	    fprintf(fp, "sigma = %g\n", freq->sdx);
	    fprintf(fp, "mu = %g\n", freq->xbar);

	    plotmin = endpt[0] - barwidth;
	    if (plotmin > freq->xbar - 3.3 * freq->sdx) {
		plotmin = freq->xbar - 3.3 * freq->sdx;
	    }
	    plotmax = endpt[K-1] + barwidth;
	    if (plotmax < freq->xbar + 3.3 * freq->sdx) {
		plotmax = freq->xbar + 3.3 * freq->sdx;
	    }
	    if (!na(freq->test)) {
		fprintf(fp, "set label \"%s:\" at graph .03, graph .97 front\n",
			_("Test statistic for normality"));
		label = make_freq_test_label(GRETL_STAT_NORMAL_CHISQ, freq->test,
					     chisq_cdf_comp(2, freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n",
			label);
		g_free(label);
	    }
	    if (real_ns > 0) {
		print_extra_literal_lines(S, ns, fp);
	    }
	} else if (dist == D_GAMMA) {
	    double var = freq->sdx * freq->sdx;

	    /* scale param = variance/mean */
	    beta = var / freq->xbar;
	    /* shape param = mean/scale */
	    alpha = freq->xbar / beta;

	    fprintf(fp, "# literal lines = %d\n", nlit);
	    fprintf(fp, "beta = %g\n", beta);
	    fprintf(fp, "alpha = %g\n", alpha);
	    plotmin = 0.0;
	    plotmax = freq->xbar + 4.0 * freq->sdx;

	    if (!na(freq->test)) {
		fprintf(fp, "set label '%s:' at graph .03, graph .97 front\n",
			_("Test statistic for gamma"));
		label = make_freq_test_label(GRETL_STAT_Z, freq->test,
					     normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n",
			label);
		g_free(label);
	    }
	    if (real_ns > 0) {
		print_extra_literal_lines(S, ns, fp);
	    }
	}

	/* adjust min, max if needed */
	if (freq->midpt[0] < plotmin) {
	    plotmin = freq->midpt[0];
	}
	if (freq->midpt[K-1] > plotmax) {
	    plotmax = freq->midpt[K-1];
	}

	fprintf(fp, "set xrange [%.10g:%.10g]\n", plotmin, plotmax);
	fputs("set key right top\n", fp);
    } else {
	/* plain frequency plot (no theoretical distribution shown) */
	lambda = 1.0 / freq->n;
	if (freq->strvals) {
	    plotmin = 0.5;
	    plotmax = K + 0.5;
	} else {
	    plotmin = freq->midpt[0] - barwidth;
	    plotmax = freq->midpt[K-1] + barwidth;
	}
	fprintf(fp, "set xrange [%.10g:%.10g]\n", plotmin, plotmax);
	maybe_set_yrange(freq, lambda, fp);
	fputs("set nokey\n", fp);

	if (real_ns > 0) {
	    fprintf(fp, "# literal lines = %d\n", real_ns);
	    print_extra_literal_lines(S, ns, fp);
	}
    }

    if (isnan(lambda)) {
	if (fp != NULL) {
	    fclose(fp);
	}
	return 1;
    }

    if (freq->strvals) {
	fprintf(fp, "set title \"%s: relative frequencies\"\n",
		freq->varname);
    } else {
	fprintf(fp, "set xlabel '%s'\n", freq->gname);
	if (dist) {
	    fprintf(fp, "set ylabel '%s'\n", _("Density"));
	} else {
	    fprintf(fp, "set ylabel '%s'\n", _("Relative frequency"));
	}
    }

    if (freq->strvals) {
	fputs("set xtics rotate by -45\n", fp);
	fputs("set xtics (", fp);
	for (i=0; i<K; i++) {
	    label = g_strdup(freq->S[i]);
	    gretl_utf8_truncate(label, 6);
	    fprintf(fp, "\"%s\" %d", label, i+1);
	    if (i < K-1) {
		fputs(", ", fp);
	    }
	    g_free(label);
	}
	fputs(")\n", fp);
    } else if (freq->discrete > 1 && K < 10 && fabs(freq->midpt[K-1]) < 1000) {
	/* few values, all integers: force integer tic marks */
	fprintf(fp, "set xtics %.0f, 1, %.0f\n", freq->midpt[0],
		freq->midpt[K-1]);
    }

    /* plot instructions */
    if (use_boxes) {
	fputs("set style fill solid 0.6\n", fp);
	strcpy(withstr, "w boxes");
    } else {
	strcpy(withstr, "w impulses lw 3");
    }

    if (!dist) {
	fprintf(fp, "plot '-' using 1:2 %s\n", withstr);
    } else if (dist == D_NORMAL) {
	label = make_freq_dist_label(dist, freq->xbar, freq->sdx);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title \"%s\" %s, \\\n"
		"1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) "
		"title \"%s\" w lines\n",
		_("relative frequency"), withstr, label);
	g_free(label);
    } else if (dist == D_GAMMA) {
	label = make_freq_dist_label(dist, alpha, beta);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s, \\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) "
		"title \"%s\" w lines\n",
		_("relative frequency"), withstr, label);
	g_free(label);
    }

    for (i=0; i<K; i++) {
	if (freq->midpt == NULL) {
	    fprintf(fp, "%d %.10g\n", i + 1, lambda * freq->f[i]);
	} else {
	    fprintf(fp, "%.10g %.10g\n", freq->midpt[i], lambda * freq->f[i]);
	}
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/**
 * plot_corrmat:
 * @corr: pointer to correlation matrix struct.
 * @opt: can use OPT_T for triangular representation.
 *
 * Produces a heatmap plot based on a correlation matrix.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_corrmat (VMatrix *corr, gretlopt opt)
{
    FILE *fp;
    double rcrit = 0.0;
    int i, j, df, n, idx;
    int allpos = 1;
    int err = 0;

    fp = open_plot_input_file(PLOT_HEATMAP, 0, &err);
    if (err) {
	return err;
    }

    n = corr->dim;

    /* are all the correlations non-negative? */
    for (i=0; i<n; i++) {
	for (j=i+1; j<n; j++) {
	    idx = ijton(i, j, n);
	    if (corr->vec[idx] < 0) {
		allpos = 0;
		break;
	    }
	}
    }

    df = corr->n - 2;
    if (df > 1) {
	/* determine 20% critical value */
	double tc = student_critval(df, 0.10);
	double t2 = tc * tc;

	rcrit = sqrt(t2 / (t2 + df));
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set title '%s'\n", _("Correlation matrix"));
    fputs("set nokey\n", fp);
    fputs("set tics nomirror\n", fp);

    if (allpos) {
	fputs("set cbrange [0:1]\n", fp);
	if (rcrit > 0) {
	    fprintf(fp, "set palette defined (0 'white', %.4f 'white', 1 'red')\n",
		    rcrit);
	} else {
	    fputs("set palette defined (0 'white', 1 'red')\n", fp);
	}
    } else {
	fputs("set cbrange [-1:1]\n", fp);
	if (rcrit > 0) {
	    fprintf(fp, "set palette defined (-1 'blue', %.4f 'white', %.4f 'white', 1 'red')\n",
		    -rcrit, rcrit);
	} else {
	    fputs("set palette defined (-1 'blue', 0 'white', 1 'red')\n", fp);
	}
    }

    if (opt & OPT_T) {
	fputs("set border 3\n", fp);
    }

    /* for grid lines */
    fputs("set x2tics 1 format '' scale 0,0.001\n", fp);
    fputs("set y2tics 1 format '' scale 0,0.001\n", fp);
    fputs("set mx2tics 2\n", fp);
    fputs("set my2tics 2\n", fp);

    /* y-axis tics */
    fputs("set ytics (", fp);
    for (i=0; i<n; i++) {
	fprintf(fp, "\"%s\" %d", corr->names[i], n-i-1);
	if (i < n - 1) {
	    fputs(", ", fp);
	}
    }
    fputs(") out\n", fp);

    /* x-axis tics */
    fputs("set xtics (", fp);
    for (i=0; i<n; i++) {
	fprintf(fp, "\"%s\" %d", corr->names[i], i);
	if (i < n - 1) {
	    fputs(", ", fp);
	}
    }
    fputs(") out\n", fp);
    fputs("set xtics rotate by 45 right\n", fp);

    /* note: "set link" requires gnuplot 5 */
    fputs("set autoscale fix\n", fp);
    fputs("set link x\n", fp);
    fputs("set link y\n", fp);
    fputs("set grid front mx2tics my2tics lw 2 lt -1 lc rgb 'white'\n", fp);

    gnuplot_missval_string(fp);
    fprintf(fp, "printcorr = %d\n", n <= 16 ? 1 : 0);

    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    if ((opt & OPT_T) && j > n-i-1) {
		write_gp_dataval(NADBL, fp, 0);
	    } else {
		idx = ijton(n-i-1, j, n);
		fprintf(fp, "%.4f ", corr->vec[idx]);
	    }
	}
	fputc('\n', fp);
    }
    fputs("EOD\n", fp);
    fputs("# end inline data\n", fp);
    fputs("if (printcorr) {\n", fp);
    fputs("plot $data matrix with image, $data matrix using 1:2:", fp);
    if (opt & OPT_T) {
	fputs("($3!=$3 ? \"\" : sprintf(\"%.1f\",$3)) with labels\n", fp);
    } else {
	fputs("(sprintf(\"%.1f\",$3)) with labels\n", fp);
    }
    fputs("} else {\n", fp);
    fputs("plot $data matrix with image\n", fp);
    fputs("}\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* print the y-axis data in the context of a forecast
   with errors plot
*/

static void fcast_print_y_data (const double *x,
				const double *y,
				int t0, int t1, int t2,
				FILE *fp)
{
    int i, t, n = t2 - t0 + 1;
    double yt;

    for (i=0; i<n; i++) {
	t = t0 + i;
	yt = t < t1 ? NADBL : y[t];
	fprintf(fp, "%.10g ", x[t]);
	write_gp_dataval(yt, fp, 1);
    }

    fputs("e\n", fp);
}

void print_user_y_data (const double *x,
			const double *y,
			int t1, int t2,
			FILE *fp)
{
    int t;

    for (t=t1; t<=t2; t++) {
	fprintf(fp, "%.10g ", x[t]);
	write_gp_dataval(y[t], fp, 1);
    }

    fputs("e\n", fp);
}

enum {
    CONF_BARS,
    CONF_FILL,
    CONF_LOW,
    CONF_HIGH
};

static void print_confband_data (const double *x,
				 const double *y,
				 const double *e,
				 int t0, int t1, int t2,
				 int mode, FILE *fp)
{
    int i, t, n = t2 - t0 + 1;
    double xt;

    for (i=0; i<n; i++) {
	t = t0 + i;
	xt = x[t];
	if (t < t1 || na(y[t]) || na(e[t])) {
	    if (mode == CONF_LOW || mode == CONF_HIGH) {
		fprintf(fp, "%.10g %s\n", xt, GPNA);
	    } else {
		fprintf(fp, "%.10g %s %s\n", xt, GPNA, GPNA);
	    }
	} else if (mode == CONF_FILL) {
	    fprintf(fp, "%.10g %.10g %.10g\n", xt, y[t] - e[t], y[t] + e[t]);
	} else if (mode == CONF_LOW) {
	    fprintf(fp, "%.10g %.10g\n", xt, y[t] - e[t]);
	} else if (mode == CONF_HIGH) {
	    fprintf(fp, "%.10g %.10g\n", xt, y[t] + e[t]);
	} else {
	    fprintf(fp, "%.10g %.10g %.10g\n", xt, y[t], e[t]);
	}
    }

    fputs("e\n", fp);
}

static void print_x_confband_data (const double *x, const double *y,
				   const double *se, const int *order,
				   double tval, int t1, int t2,
				   int mode, FILE *fp)
{
    int i, t, n = t2 - t1 + 1;
    double et;

    for (i=0; i<n; i++) {
	t = order[i];
	if (na(y[t]) || na(se[t])) {
	    if (!na(x[t])) {
		fprintf(fp, "%.10g %s\n", x[t], GPNA);
	    }
	} else {
	    et = tval * se[t];
	    if (mode == CONF_LOW) {
		fprintf(fp, "%.10g %.10g\n", x[t], y[t] - et);
	    } else if (mode == CONF_HIGH) {
		fprintf(fp, "%.10g %.10g\n", x[t], y[t] + et);
	    }
	}
    }

    fputs("e\n", fp);
}

struct fsorter {
    int obs;
    double y;
};

static int compare_fs (const void *a, const void *b)
{
    const struct fsorter *fa = a;
    const struct fsorter *fb = b;

    return (fa->y > fb->y) - (fa->y < fb->y);
}

static void print_filledcurve_line (const char *title,
				    const char *rgb,
				    FILE *fp)
{
    char cstr[10];

    if (rgb != NULL && *rgb != '\0') {
	*cstr = '\0';
	strncat(cstr, rgb, 9);
    } else {
	print_rgb_hash(cstr, get_shadecolor());
    }

    if (title == NULL) {
	fprintf(fp, "'-' using 1:2:3 notitle lc rgb \"%s\" w filledcurve, \\\n",
		cstr);
    } else {
	fprintf(fp, "'-' using 1:2:3 title '%s' lc rgb \"%s\" w filledcurve, \\\n",
		title, cstr);
    }
}

/* note: if @opt includes OPT_H, that says to show fitted
   values for the pre-forecast range
*/

int plot_fcast_errs (const FITRESID *fr, const double *maxerr,
		     const DATASET *dset, gretlopt opt)
{
    FILE *fp = NULL;
    const double *obs = NULL;
    GptFlags flags = 0;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int use_fill = 0, use_lines = 0;
    int do_errs = (maxerr != NULL);
    gchar *cistr = NULL;
    int t2 = fr->t2;
    int t1, yhmin;
    int t, n, err = 0;

    /* note: yhmin is the first obs at which to start plotting y-hat */
    if (do_errs) {
	t1 = fr->t0;
	yhmin = (opt & OPT_H)? fr->t0 : fr->t1;
    } else {
	t1 = (fr->t0 >= 0)? fr->t0 : 0;
	/* was: yhmin = t1; */
	yhmin = (opt & OPT_H)? t1 : fr->t1;
    }

    /* don't graph empty trailing portion of forecast */
    for (t=fr->t2; t>=t1; t--) {
	if (na(fr->actual[t]) && na(fr->fitted[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    n = t2 - t1 + 1;

    if (n < 3) {
	/* we won't draw a graph for 2 data points or less */
	return 1;
    }

    obs = gretl_plotx(dset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if (dataset_is_time_series(dset)) {
	flags |= GPT_LETTERBOX;
    }

    fp = open_plot_input_file(PLOT_FORECAST, flags, &err);
    if (err) {
	return err;
    }

    /* check that we have any values for the actual var */
    for (t=t1; t<=t2; t++) {
	if (!na(fr->actual[t])) {
	    depvar_present = 1;
	    break;
	}
    }

    if (do_errs) {
	if (opt & OPT_F) {
	    use_fill = 1;
	} else if (opt & OPT_L) {
	    use_lines = 1;
	}
    }

    gretl_minmax(t1, t2, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    if (dataset_is_time_series(dset)) {
	fprintf(fp, "# timeseries %d (letterbox)\n", dset->pd);
    }

    if (do_errs && !use_fill && !use_lines && n > 150) {
	use_fill = 1;
    }

    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("plot \\\n", fp);

    if (do_errs) {
	cistr = g_strdup_printf(_("%g percent interval"), 100 * (1 - fr->alpha));
    }

    if (use_fill) {
	/* plot the confidence band first so the other lines
	   come out on top */
	if (do_errs) {
	    print_filledcurve_line(cistr, NULL, fp);
	}
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1, \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines lt 2\n", _("forecast"));
    } else {
	/* plot confidence bands last */
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines", _("forecast"));
	if (do_errs) {
	    if (use_lines) {
		fprintf(fp, ", \\\n'-' using 1:2 title '%s' w lines, \\\n",
			cistr);
		fputs("'-' using 1:2 notitle '%s' w lines lt 3\n", fp);
	    } else {
		fprintf(fp, ", \\\n'-' using 1:2:3 title '%s' w errorbars\n",
			cistr);
	    }
	} else {
	    fputc('\n', fp);
	}
    }

    g_free(cistr);

    gretl_push_c_numeric_locale();

    /* write out the inline data, the order depending on whether
       or not we're using fill style for the confidence bands
    */

    if (use_fill) {
	if (do_errs) {
	    print_confband_data(obs, fr->fitted, maxerr,
				t1, yhmin, t2, CONF_FILL, fp);
	}
	if (depvar_present) {
	    fcast_print_y_data(obs, fr->actual, t1, t1, t2, fp);
	}
	fcast_print_y_data(obs, fr->fitted, t1, yhmin, t2, fp);
    } else {
	if (depvar_present) {
	    fcast_print_y_data(obs, fr->actual, t1, t1, t2, fp);
	}
	fcast_print_y_data(obs, fr->fitted, t1, yhmin, t2, fp);
	if (do_errs) {
	    if (use_lines) {
		print_confband_data(obs, fr->fitted, maxerr,
				    t1, yhmin, t2, CONF_LOW, fp);
		print_confband_data(obs, fr->fitted, maxerr,
				    t1, yhmin, t2, CONF_HIGH, fp);
	    } else {
		print_confband_data(obs, fr->fitted, maxerr,
				    t1, yhmin, t2, CONF_BARS, fp);
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

gretlRGB numeric_color_from_string (const char *s, int *err)
{
    gretlRGB ret = 0;
    char *test = NULL;
    int len = strlen(s);
    int offset = 0;

    if (s[0] == '#') {
	offset = 1;
    } else if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
	offset = 2;
    }

    if (offset > 0) {
	/* should be numeric [AA]RRGGBB */
	const char *p = s + offset;

	len -= offset;
	if (len != 6 && len != 8) {
	     *err = invalid_field_error(s);
	} else {
	    ret = (gretlRGB) strtoul(p, &test, 16);
	    if (*test != '\0') {
		*err = invalid_field_error(s);
	    }
	}
	return ret;
    }

    /* otherwise we should have a gnuplot color name */
    if (len < 3 || len > 17) {
	*err = invalid_field_error(s);
    } else {
	char fname[FILENAME_MAX];
	FILE *fp;

	sprintf(fname, "%sdata%cgnuplot%cgpcolors.txt",
		gretl_home(), SLASH, SLASH);
	fp = gretl_fopen(fname, "r");

	if (fp == NULL) {
	    *err = E_FOPEN;
	} else {
	    char line[32], cname[18];
	    int found = 0;
	    guint32 u;

	    while (fgets(line, sizeof line, fp)) {
		if (sscanf(line, "%s %x", cname, &u) == 2 &&
		    strcmp(s, cname) == 0) {
		    ret = u;
		    found = 1;
		    break;
		}
	    }
	    fclose(fp);
	    if (!found) {
		*err = invalid_field_error(s);
		ret = 0;
	    }
	}
    }

    return ret;
}

static int *get_x_sorted_order (const FITRESID *fr,
				const double *x,
				int t1, int *pt2)
{
    int *order = NULL;
    struct fsorter *fs;
    int nmiss = 0;
    int t2 = *pt2;
    int n = t2 - t1 + 1;
    int i, t;

    for (t=t1; t<=t2; t++) {
	if (na(fr->actual[t])) {
	    nmiss++;
	}
    }

    fs = malloc(n * sizeof *fs);
    if (fs == NULL) {
	return NULL;
    }

    order = malloc(n * sizeof *order);
    if (order == NULL) {
	free(fs);
	return NULL;
    }

    for (i=0, t=t1; t<=t2; t++, i++) {
	fs[i].obs = t;
	fs[i].y = x[t];
    }

    qsort(fs, n, sizeof *fs, compare_fs);

    for (i=0; i<n; i++) {
	order[i] = fs[i].obs;
    }

    free(fs);

    if (nmiss > 0) {
	/* chop off trailing NAs */
	*pt2 = (n - nmiss) + t1 - 1;
    }

    return order;
}

static void print_x_ordered_data (const double *x, const double *y,
				  const int *order, int t1, int t2,
				  FILE *fp)
{
    int i, t, n = t2 - t1 + 1;

    for (i=0; i<n; i++) {
	t = order[i];
	if (na(x[t])) {
	    continue;
	} else if (na(y[t])) {
	    fprintf(fp, "%.10g %s\n", x[t], GPNA);
	} else {
	    fprintf(fp, "%.10g %.10g\n", x[t], y[t]);
	}
    }

    fputs("e\n", fp);
}

/* Plotting routine for a simple regression where we want to show
   actual and forecast y, plus confidence bands, against x.
   We use lines to indicate the confidence bands.
*/

int plot_simple_fcast_bands (const MODEL *pmod,
			     const FITRESID *fr,
			     const DATASET *dset,
			     gretlopt opt)
{
    FILE *fp = NULL;
    const double *x = NULL;
    int *order = NULL;
    GptFlags flags = 0;
    double a, xmin, xmax, xrange, tval;
    int xv = pmod->list[3];
    gchar *cistr;
    int t2 = fr->t2;
    int t1, yhmin;
    int t, n, err = 0;

    /* note: yhmin is the first obs at which to start plotting y-hat */
    t1 = fr->t0;
    yhmin = (opt & OPT_H)? fr->t0 : fr->t1;

    /* don't graph empty trailing portion of forecast */
    for (t=fr->t2; t>=t1; t--) {
	if (na(fr->actual[t]) && na(fr->fitted[t])) {
	    t2--;
	} else {
	    break;
	}
    }

    n = t2 - t1 + 1;

    if (n < 3) {
	/* won't draw a graph for 2 data points or less */
	return 1;
    }

    x = dset->Z[xv];

    order = get_x_sorted_order(fr, x, t1, &t2);
    if (order == NULL) {
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_FORECAST, flags, &err);
    if (err) {
	return err;
    }

    gretl_minmax(t1, t2, x, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);
    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    fprintf(fp, "set xlabel \"%s\"\n", dset->varname[xv]);
    fprintf(fp, "set ylabel \"%s\"\n", fr->depvar);

    fputs("set key left top\n", fp);
    fputs("plot \\\n", fp);

    a = 100 * (1 - fr->alpha);
    tval = student_critval(fr->df, fr->alpha / 2);

    if (opt & OPT_M) {
	cistr = g_strdup_printf(_("%g%% interval for mean"), a);
    } else {
	cistr = g_strdup_printf(_("%g percent interval"), a);
    }

    fputs("'-' using 1:2 notitle w points, \\\n", fp);
    fputs("'-' using 1:2 notitle w lines, \\\n", fp);
    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n", cistr);
    fputs("'-' using 1:2 notitle '%s' w lines lt 3\n", fp);
    g_free(cistr);

    gretl_push_c_numeric_locale();

    print_x_ordered_data(x, fr->actual, order, t1, t2, fp);
    print_x_ordered_data(x, fr->fitted, order, yhmin, t2, fp);
    print_x_confband_data(x, fr->fitted, fr->sderr, order,
			  tval, yhmin, t2, CONF_LOW, fp);
    print_x_confband_data(x, fr->fitted, fr->sderr, order,
			  tval, yhmin, t2, CONF_HIGH, fp);

    gretl_pop_c_numeric_locale();

    free(order);

    return finalize_plot_input_file(fp);
}

#ifndef min
# define min(x,y) (((x)<(y))? (x):(y))
#endif

#ifndef max
# define max(x,y) (((x)>(y))? (x):(y))
#endif

int plot_tau_sequence (const MODEL *pmod, const DATASET *dset,
		       int k)
{
    FILE *fp;
    gretl_matrix *tau = gretl_model_get_data(pmod, "rq_tauvec");
    gretl_matrix *B = gretl_model_get_data(pmod, "rq_sequence");
    double tau_i, bi, se, blo, bhi;
    double alpha, cval, tcrit, olsband;
    double ymin[2], ymax[2];
    gchar *tmp;
    int ntau, bcols;
    int i, j, err = 0;

    if (tau == NULL || B == NULL) {
	return E_DATA;
    }

    ntau = gretl_vector_get_length(tau);
    if (ntau == 0) {
	return E_DATA;
    }

    fp = open_plot_input_file(PLOT_RQ_TAU, 0, &err);
    if (err) {
	return err;
    }

    bcols = gretl_matrix_cols(B);

    alpha = gretl_model_get_double(pmod, "rq_alpha");
    if (na(alpha)) {
	alpha = .05;
    }

    cval = 100 * (1 - alpha);
    tcrit = student_cdf_inverse(pmod->dfd, 1 - alpha/2);
    olsband = tcrit * pmod->sderr[k];

    /* Try to figure best placement of key */

    j = k * ntau;
    if (bcols == 3) {
	blo = gretl_matrix_get(B, j, 1);
	bhi = gretl_matrix_get(B, j, 2);
    } else {
	bi = gretl_matrix_get(B, j, 0);
	se = gretl_matrix_get(B, j, 1);
	blo = bi - tcrit * se;
	bhi = bi + tcrit * se;
    }
    ymin[0] = min(blo, pmod->coeff[k] - olsband);
    ymax[0] = max(bhi, pmod->coeff[k] + olsband);

    j += ntau - 1;
    if (bcols == 3) {
	blo = gretl_matrix_get(B, j, 1);
	bhi = gretl_matrix_get(B, j, 2);
    } else {
	bi = gretl_matrix_get(B, j, 0);
	se = gretl_matrix_get(B, j, 1);
	blo = bi - tcrit * se;
	bhi = bi + tcrit * se;
    }
    ymin[1] = min(blo, pmod->coeff[k] - olsband);
    ymax[1] = max(bhi, pmod->coeff[k] + olsband);

    fputs("set xrange [0.0:1.0]\n", fp);
    fputs("set xlabel 'tau'\n", fp);

    tmp = g_strdup_printf(_("Coefficient on %s"),
			  series_get_graph_name(dset, pmod->list[k+2]));
    fprintf(fp, "set title \"%s\"\n", tmp);
    g_free(tmp);

    fputs("set style fill solid 0.5\n", fp);

    if (ymax[0] < .88 * ymax[1]) {
	fputs("set key left top\n", fp);
    } else if (ymax[1] < .88 * ymax[0]) {
	fputs("set key right top\n", fp);
    } else if (ymin[0] < .88 * ymin[1]) {
	fputs("set key right bottom\n", fp);
    } else {
	fputs("set key left bottom\n", fp);
    }

    fputs("plot \\\n", fp);

    /* plot the rq confidence band first so the other lines
       come out on top */
    print_filledcurve_line(NULL, NULL, fp);

    /* rq estimates */
    tmp = g_strdup_printf(_("Quantile estimates with %g%% band"), cval);
    fprintf(fp, "'-' using 1:2 title '%s' w lp, \\\n", tmp);
    g_free(tmp);

    /* numeric output coming up! */
    gretl_push_c_numeric_locale();

    /* ols estimate plus (1 - alpha) band */
    tmp = g_strdup_printf(_("OLS estimate with %g%% band"), cval);
    fprintf(fp, "%g title '%s' w l, \\\n", pmod->coeff[k], tmp);
    g_free(tmp);
    fprintf(fp, "%g notitle w l dt 2, \\\n", pmod->coeff[k] + olsband);
    fprintf(fp, "%g notitle w l dt 2\n", pmod->coeff[k] - olsband);

    /* write out the interval values */

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
	tau_i = gretl_vector_get(tau, i);
	if (bcols == 3) {
	    blo = gretl_matrix_get(B, j, 1);
	    bhi = gretl_matrix_get(B, j, 2);
	} else {
	    bi = gretl_matrix_get(B, j, 0);
	    se = gretl_matrix_get(B, j, 1);
	    blo = bi - tcrit * se;
	    bhi = bi + tcrit * se;
	}
	fprintf(fp, "%.10g %.10g %.10g\n", tau_i, blo, bhi);
    }
    fputs("e\n", fp);

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
	tau_i = gretl_vector_get(tau, i);
	bi = gretl_matrix_get(B, j, 0);
	fprintf(fp, "%.10g %.10g\n", tau_i, bi);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int garch_resid_plot (const MODEL *pmod, const DATASET *dset)
{
    FILE *fp = NULL;
    const double *obs;
    const double *h;
    double sd2;
    int t, err = 0;

    h = gretl_model_get_data(pmod, "garch_h");
    if (h == NULL) {
	return E_DATA;
    }

    obs = gretl_plotx(dset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_GARCH, 0, &err);
    if (err) {
	return err;
    }

    fputs("set key left top\n", fp);

    fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines, \\\n"
	    "'-' using 1:2 title '%s' w lines lt 2, \\\n"
	    "'-' using 1:2 notitle w lines lt 2\n",
	    _("residual"), _("+- sqrt(h(t))"));

    gretl_push_c_numeric_locale();

    for (t=pmod->t1; t<=pmod->t2; t++) {
	fprintf(fp, "%.10g %.10g\n", obs[t], pmod->uhat[t]);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = -sqrt(h[t]);
	fprintf(fp, "%.10g %.10g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = sqrt(h[t]);
	fprintf(fp, "%.10g %.10g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int rmplot (const int *list, DATASET *dset,
	    gretlopt opt, PRN *prn)
{
    int (*range_mean_graph) (int, const DATASET *,
			     gretlopt, PRN *);

    range_mean_graph = get_plugin_function("range_mean_graph");
    if (range_mean_graph == NULL) {
	return 1;
    }

    return range_mean_graph(list[1], dset, opt, prn);
}

int hurstplot (const int *list, DATASET *dset, gretlopt opt, PRN *prn)
{
    int (*hurst_exponent) (int, const DATASET *, gretlopt, PRN *);

    hurst_exponent = get_plugin_function("hurst_exponent");
    if (hurst_exponent == NULL) {
	return 1;
    }

    return hurst_exponent(list[1], dset, opt, prn);
}

static void get_multiplot_layout (int n, int tseries,
				  int *rows, int *cols)
{
    if (n < 3) {
	if (tseries) {
	    *cols = 1;
	    *rows = 2;
	} else {
	    *cols = 2;
	    *rows = 1;
	}
    } else if (n < 5) {
	*cols = *rows = 2;
    } else if (n < 7) {
	*cols = 3;
	*rows = 2;
    } else if (n < 10) {
	*cols = *rows = 3;
    } else if (n < 13) {
	*cols = 4;
	*rows = 3;
    } else if (n < 17) {
	*cols = *rows = 4;
    } else {
	*cols = *rows = 0;
    }
}

static int panel_ytic_width (double ymin, double ymax)
{
    char s1[16], s2[16];
    int n1, n2;

    if (ymin < 0 && ymax > 0) {
	sprintf(s1, "% g", ymin);
	sprintf(s2, "% g", ymax);
    } else {
	sprintf(s1, "%g", ymin);
	sprintf(s2, "%g", ymax);
    }

    n1 = strlen(s1);
    n2 = strlen(s2);

    return (n1 > n2)? n1 : n2;
}

/* Panel: produce a time-series plot for the group mean of the
   series in question.
*/

static int panel_means_ts_plot (const int vnum,
				const char *literal,
				const DATASET *dset,
				gretlopt opt)
{
    DATASET *gset;
    int nunits, T = dset->pd;
    int list[2] = {1, 1};
    gchar *my_literal = NULL;
    gchar *title = NULL;
    int i, t, s, s0;
    int err = 0;

    nunits = panel_sample_size(dset);

    gset = create_auxiliary_dataset(2, T, 0);
    if (gset == NULL) {
	return E_ALLOC;
    }

    strcpy(gset->varname[1], dset->varname[vnum]);
    series_set_display_name(gset, 1, series_get_display_name(dset, vnum));

    time_series_from_panel(gset, dset);

    s0 = dset->t1;

    for (t=0; t<T; t++) {
	double xit, xsum = 0.0;
	int n = 0;

	for (i=0; i<nunits; i++) {
	    s = s0 + i * T + t;
	    xit = dset->Z[vnum][s];
	    if (!na(xit)) {
		xsum += xit;
		n++;
	    }
	}
	gset->Z[1][t] = (n == 0)? NADBL : xsum / n;
    }

    opt |= (OPT_O | OPT_T); /* use lines, time series */

    title = g_strdup_printf(_("mean %s"),
			    series_get_graph_name(dset, vnum));
    my_literal = g_strdup_printf("set ylabel \"%s\" ; set xlabel ;",
				 title);
    if (literal != NULL) {
	gchar *all_lit;

	all_lit = g_strdup_printf("%s %s", my_literal, literal);
	err = gnuplot(list, all_lit, gset, opt);
	g_free(all_lit);
    } else {
	err = gnuplot(list, my_literal, gset, opt);
    }

    g_free(title);
    g_free(my_literal);
    destroy_dataset(gset);

    return err;
}

int panel_means_XY_scatter (const int *list,
                            const char *literal,
                            const DATASET *dset,
			    gretlopt opt)
{
    DATASET *gset;
    int N, T = dset->pd;
    int glist[3] = {2, 1, 2};
    gchar *my_literal = NULL;
    int grpnames = 0;
    int yvar, xvar;
    int i, t, s;
    int err = 0;

    if (list == NULL || list[0] != 2) {
	return E_DATA;
    }

    N = panel_sample_size(dset);

    gset = create_auxiliary_dataset(3, N, 0);
    if (gset == NULL) {
	return E_ALLOC;
    }

    /* If we have valid panel group names, use them
       as obs markers here */
    grpnames = panel_group_names_ok(dset);
    if (grpnames) {
	dataset_allocate_obs_markers(gset);
    }

    yvar = list[1];
    xvar = list[2];

    strcpy(gset->varname[1], dset->varname[yvar]);
    series_set_display_name(gset, 1, series_get_display_name(dset, yvar));

    strcpy(gset->varname[2], dset->varname[xvar]);
    series_set_display_name(gset, 2, series_get_display_name(dset, xvar));

    s = dset->t1;

    for (i=0; i<N; i++) {
	double yit, ysum = 0.0;
	double xit, xsum = 0.0;
	int ny = 0, nx = 0;
	int s0 = s;

	for (t=0; t<T; t++) {
	    yit = dset->Z[yvar][s];
	    xit = dset->Z[xvar][s];
	    if (!na(yit)) {
		ysum += yit;
		ny++;
	    }
	    if (!na(xit)) {
		xsum += xit;
		nx++;
	    }
	    s++;
	}
	gset->Z[1][i] = ny == 0 ? NADBL : ysum / ny;
	gset->Z[2][i] = nx == 0 ? NADBL : xsum / nx;
	if (gset->S != NULL) {
	    strcpy(gset->S[i], get_panel_group_name(dset, s0));
	}
    }

    my_literal = g_strdup_printf("set title \"%s\";", _("Group means"));

    if (literal != NULL) {
        gchar *all_lit = g_strdup_printf("%s %s", my_literal, literal);

        err = gnuplot(glist, all_lit, gset, opt);
        g_free(all_lit);
    } else {
        err = gnuplot(glist, my_literal, gset, opt);
    }

    g_free(my_literal);
    destroy_dataset(gset);

    return err;
}

static void copy_string_stripped (char *targ, const char *src,
				  int strip)
{
    char *tmp = gretl_strdup(src);
    size_t len = strlen(tmp) - strip;

    strcpy(targ, gretl_utf8_truncate(tmp, len));
    free(tmp);
}

/* Here we're trying to find out if the observation labels
   for a panel dataset are such that they uniquely identify
   the units/individuals (e.g. country or city names,
   repeated for each time-series observation on the given
   unit).
*/

static int dataset_has_panel_labels (const DATASET *dset,
				     int maxlen, int *use,
				     int *strip)
{
    int t, u, ubak = -1;
    int len, lmax = 0, fail = 0;
    int ret = 0;

    if (dset->S == NULL) {
	return 0;
    }

    for (t=dset->t1; t<=dset->t2 && !fail; t++) {
	u = t / dset->pd;
	if (u == ubak && strcmp(dset->S[t], dset->S[t-1])) {
	    /* same unit, different label: no */
	    fail = 1;
	} else if (ubak >= 0 && u != ubak &&
		   !strcmp(dset->S[t], dset->S[t-1])) {
	    /* different unit, same label: no */
	    fail = 2;
	}
	if (!fail) {
	    len = strlen(dset->S[t]);
	    if (len > lmax) {
		lmax = len;
	    }
	}
	ubak = u;
    }

    if (!fail) {
	/* the full obs labels satisfy the criterion,
	   but are they perhaps too long? */
	if (maxlen > 0 && lmax > maxlen) {
	    ret = 0;
	} else {
	    ret = 1;
	}
    } else if (fail == 1) {
	/* Try for a leading portion of the obs labels: for
	   example we might have AUS1990, AUS1991, ...
	   followed by USA1990, USA1991, ... or some such.
	   We try to identify a trailing portion of the obs
	   string that varies by time, and which should be
	   omitted in forming "panel labels".
	*/
	const char *s;
	int i, n, len2t, len2 = 0;
	int obslen = 0;

	fail = 0;
	for (t=dset->t1; t<=dset->t2 && !fail; t++) {
	    s = dset->S[t];
	    n = strlen(s);
	    len2t = 0;
	    for (i=n-1; i>0; i--) {
		if (isdigit(s[i]) || s[i] == ':' ||
		    s[i] == '-' || s[i] == '_') {
		    len2t++;
		} else {
		    break;
		}
	    }
	    if (len2t == 0) {
		/* no "tail" string (e.g. year) */
		fail = 1;
	    } else if (len2 == 0) {
		/* starting */
		len2 = len2t;
		obslen = n;
	    } else if (len2t != len2) {
		/* the "tails" don't have a common length */
		fail = 1;
	    } else if (n != obslen) {
		/* the obs strings are of differing lengths */
		obslen = 0;
	    }
	}

	if (!fail) {
	    char s0[OBSLEN], s1[OBSLEN];

	    if (obslen > 0) {
		*use = obslen - len2;
	    } else {
		*strip = len2;
	    }
	    /* now check that the leading portion really
	       is in common for each unit/individual
	    */
	    *s0 = '\0';
	    ubak = -1;
	    for (t=dset->t1; t<=dset->t2 && !fail; t++) {
		u = t / dset->pd;
		*s1 = '\0';
		if (*use > 0) {
		    strncat(s1, dset->S[t], *use);
		} else {
		    copy_string_stripped(s1, dset->S[t], *strip);
		}
		if (u == ubak && strcmp(s1, s0)) {
		    /* same unit, different label: no */
		    fail = 1;
		} else if (ubak >= 0 && u != ubak && !strcmp(s1, s0)) {
		    /* different unit, same label: no */
		    fail = 2;
		}
		if (!fail && maxlen > 0 && strlen(s1) > maxlen) {
		    fail = 3;
		}
		strcpy(s0, s1);
		ubak = u;
	    }
	    if (fail) {
		*use = *strip = 0;
	    } else {
		ret = 1;
	    }
	}
    }

    /* There's a loophole above: unit m might have the same
       label as some other unit, although we've checked that
       it doesn't have the same label as unit m - 1. But
       that seems ike a corner case and I can't be bothered
       checking for it right now.
    */

    return ret;
}

/* Panel: plot one series using separate lines for each
   cross-sectional unit. The individuals' series are overlaid, in the
   same manner as a plot of several distinct time series. To do
   this we construct on the fly a notional time-series dataset.

   But note: if it turns out the series in question is invariant
   across groups, just show a single line.
*/

static int panel_overlay_ts_plot (const int vnum,
				  const char *literal,
				  const DATASET *dset,
				  gretlopt opt)
{
    DATASET *gset;
    int u0, nunits, T = dset->pd;
    int *list = NULL;
    gchar *my_literal = NULL;
    gchar *title = NULL;
    series_table *gst = NULL;
    const char *sval;
    int vg = 0;
    int panel_labels = 0;
    int maxlen = 0;
    int single_series;
    int use = 0, strip = 0;
    int i, s, s0;
    int err = 0;

    single_series = series_is_group_invariant(dset, vnum);

    if (single_series) {
	nunits = 1;
    } else {
	nunits = panel_sample_size(dset);
    }

    u0 = dset->t1 / T;

    gset = create_auxiliary_dataset(nunits + 1, T, OPT_B);
    if (gset == NULL) {
	return E_ALLOC;
    }

    /* add time series info to @gset */
    time_series_from_panel(gset, dset);

    list = gretl_consecutive_list_new(1, nunits);
    if (list == NULL) {
	destroy_dataset(gset);
	return E_ALLOC;
    }

    if (nunits > 80) {
	/* FIXME calibrate this properly */
	maxlen = 3;
    }

    if (!single_series) {
	gst = get_panel_group_table(dset, maxlen, &vg);
	if (gst == NULL && dset->S != NULL) {
	    /* maybe we have obs markers that are usable */
	    panel_labels =
		dataset_has_panel_labels(dset, maxlen, &use, &strip);
	}
    }

    s0 = dset->t1;

    /* compose varnames and transcribe data */
    for (i=0; i<nunits; i++) {
	s = s0 + i * T;
	if (single_series) {
	    strcpy(gset->varname[i+1], dset->varname[vnum]);
	} else if (gst != NULL) {
	    /* look up the string for this unit/group */
	    sval = series_table_get_string(gst, dset->Z[vg][s]);
	    if (sval != NULL) {
		strncat(gset->varname[i+1], sval, VNAMELEN-1);
	    } else {
		sprintf(gset->varname[i+1], "%d", u0+i+1);
	    }
	} else if (panel_labels) {
	    if (use > 0) {
		strncat(gset->varname[i+1], dset->S[s], use);
	    } else if (strip > 0) {
		copy_string_stripped(gset->varname[i+1], dset->S[s], strip);
	    } else {
		strcpy(gset->varname[i+1], dset->S[s]);
	    }
	} else {
	    sprintf(gset->varname[i+1], "%d", u0+i+1);
	}
	/* borrow the appropriate chunk of data */
	gset->Z[i+1] = dset->Z[vnum] + s;
    }

    if (nunits > 9 && T < 50) {
	opt |= OPT_P; /* lines/points */
    } else {
	opt |= OPT_O; /* use lines */
    }

    if (!single_series) {
	const char *gname = series_get_graph_name(dset, vnum);
	const char *vname = panel_group_names_varname(dset);

	if (vname != NULL) {
	    title = g_strdup_printf("%s by %s", gname, vname);
	} else {
	    title = g_strdup_printf("%s by group", gname);
	}
	my_literal = g_strdup_printf("set title \"%s\" ; set xlabel ;", title);
    }

    if (nunits > 80) {
	/* set file-scope global */
	xwide = 1;
    }

    if (my_literal != NULL && literal != NULL) {
	gchar *all_lit = g_strdup_printf("%s %s", my_literal, literal);

	err = gnuplot(list, all_lit, gset, opt | OPT_T);
	g_free(all_lit);
    } else if (literal != NULL) {
	err = gnuplot(list, literal, gset, opt | OPT_T);
    } else {
	err = gnuplot(list, my_literal, gset, opt | OPT_T);
    }

    xwide = 0;
    g_free(title);
    g_free(my_literal);
    destroy_dataset(gset);
    free(list);

    return err;
}

/* Panel: plot one variable as a time series, with separate plots for
   each cross-sectional unit.  By default we arrange the plots in a
   grid, but if OPT_V is given we make each plot full width and
   stack the plots vertically on the "page".
*/

static int panel_grid_ts_plot (int vnum, const DATASET *dset,
			       gretlopt opt)
{
    FILE *fp = NULL;
    int w, rows, cols;
    const double *y, *x = NULL;
    const char *vname;
    char uname[OBSLEN];
    double xt, yt, ymin, ymax, incr;
    int u0, nunits, T = dset->pd;
    int n_ok_units = 0;
    int panel_labels = 0;
    int use = 0, strip = 0;
    int *badlist = NULL;
    int i, s, t, t0;
    int err = 0;

    n_ok_units = nunits = panel_sample_size(dset);
    u0 = dset->t1 / dset->pd;
    y = dset->Z[vnum];

    /* check for "blank" units */
    t0 = dset->t1;
    for (i=0; i<nunits; i++) {
	int ok = 0;

	for (t=0; t<T; t++) {
	    if (!na(y[t+t0])) {
		ok = 1;
		break;
	    }
	}
	if (!ok) {
	    badlist = gretl_list_append_term(&badlist, i);
	    n_ok_units--;
	}
	t0 += T;
    }

    if (n_ok_units < 2) {
	free(badlist);
	return E_MISSDATA;
    }

    if (opt & OPT_V) {
	int xvar = plausible_panel_time_var(dset);

	if (xvar > 0) {
	    x = dset->Z[xvar];
	}
	cols = 1;
	rows = n_ok_units;
    } else {
	get_multiplot_layout(n_ok_units, 0, &rows, &cols);
    }

    if (rows == 0 || cols == 0) {
	return E_DATA;
    }

    maybe_set_small_font(nunits);

    fp = open_plot_input_file(PLOT_PANEL, 0, &err);
    if (err) {
	return err;
    }

    if (dset->S != NULL) {
	panel_labels = dataset_has_panel_labels(dset, 0, &use, &strip);
    }

    vname = dset->varname[vnum];
    gretl_minmax(dset->t1, dset->t2, y, &ymin, &ymax);
    w = panel_ytic_width(ymin, ymax);

    fputs("set key left top\n", fp);
    gnuplot_missval_string(fp);
    fputs("set xtics nomirror\n", fp);
    fputs("set ytics nomirror\n", fp);
    fprintf(fp, "set format y \"%%%dg\"\n", w);
    fprintf(fp, "set multiplot layout %d,%d\n", rows, cols);

    if (opt & OPT_V) {
	fputs("set noxlabel\n", fp);
    } else {
	fprintf(fp, "set xlabel '%s'\n", _("time"));
    }

    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    t0 = dset->t1;
    for (i=0; i<nunits; i++) {
	if (in_gretl_list(badlist, i)) {
	    t0 += T;
	    continue;
	}
	if (panel_labels) {
	    *uname = '\0';
	    s = (u0 + i) * dset->pd;
	    if (use > 0) {
		strncat(uname, dset->S[s], use);
	    } else if (strip > 0) {
		copy_string_stripped(uname, dset->S[s], strip);
	    } else {
		strcpy(uname, dset->S[s]);
	    }
	} else {
	    sprintf(uname, "%d", u0+i+1);
	}
	if (opt & OPT_V) {
	    gretl_minmax(t0, t0 + T - 1, y, &ymin, &ymax);
	    incr = (ymax - ymin) / 2.0;
	    fprintf(fp, "set ytics %g\n", incr);
	    fprintf(fp, "set ylabel '%s (%s)'\n", vname, uname);
	} else {
	    fprintf(fp, "set title '%s (%s)'\n", vname, uname);
	}

	fputs("plot \\\n'-' using 1:2 notitle w lines\n", fp);

	for (t=0; t<T; t++) {
	    if (x != NULL) {
		xt = x[t+t0];
	    } else {
		xt = t + 1;
	    }
	    yt = y[t+t0];
	    if (na(yt)) {
		fprintf(fp, "%g %s\n", xt, GPNA);
	    } else {
		fprintf(fp, "%g %.10g\n", xt, yt);
	    }
	}
	fputs("e\n", fp);
	t0 += T;
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

    free(badlist);

    return finalize_plot_input_file(fp);
}

int gretl_panel_ts_plot (int vnum, DATASET *dset, gretlopt opt)
{
    if (opt & OPT_S) {
	return panel_grid_ts_plot(vnum, dset, opt);
    } else if (opt & OPT_M) {
	/* group means */
	opt &= ~OPT_M;
	opt |= OPT_S;
	return panel_means_ts_plot(vnum, NULL, dset, opt);
    } else {
	return panel_overlay_ts_plot(vnum, NULL, dset, opt);
    }
}

/* The following implements the script command "panplot" */

int cli_panel_plot (const int *list, const char *literal,
		    const DATASET *dset, gretlopt opt)
{
    int N, vnum = list[1];
    int err;

    /* condition on multi_unit_panel_sample() ? */

    if (!dataset_is_panel(dset)) {
	gretl_errmsg_set(_("This command needs panel data"));
	err = E_DATA;
    } else {
	err = incompatible_options(opt, OPT_M | OPT_V | OPT_S |
				   OPT_D | OPT_A | OPT_B | OPT_C);
    }
    if (err) {
	return err;
    }

    N = panel_sample_size(dset);

    /* check for too many groups */
    if ((opt & (OPT_V | OPT_S)) && N > 130) {
	err = E_BADOPT;
    } else if ((opt & OPT_B) && N > 150) {
	err = E_BADOPT;
    } else if ((opt & OPT_D) && N > 16) {
	err = E_BADOPT;
    } else if ((opt & OPT_A) && N > 6) {
	err = E_BADOPT;
    }
    if (err) {
	gretl_errmsg_set("Too many groups for the specified plot");
	return err;
    }

    /* select a default if no panplot-specific option given */
    if (!(opt & (OPT_M | OPT_V | OPT_S | OPT_D |
		 OPT_A | OPT_B | OPT_C))) {
	if (N <= 130) {
	    opt |= OPT_V; /* --overlay */
	} else if (N <= 150) {
	    opt |= OPT_B; /* --boxplots */
	} else {
	    opt |= OPT_M; /* --means */
	}
    }

    if (opt & (OPT_U | OPT_b)) {
	/* handle output or outbuf spec */
	gretlopt outopt = (opt & OPT_U) ? OPT_U : OPT_b;
	const char *s = get_optval_string(PANPLOT, outopt);
	int pci = (opt & (OPT_B | OPT_C)) ? BXPLOT : GNUPLOT;

	if (s != NULL) {
	    set_optval_string(pci, outopt, s);
	}
    }

    if (opt & OPT_M) {
	/* --means */
	opt &= ~OPT_M;
	opt |= OPT_S;
	err = panel_means_ts_plot(vnum, literal, dset, opt);
    } else if (opt & OPT_V) {
	/* --overlay */
	opt &= ~OPT_V;
	err = panel_overlay_ts_plot(vnum, literal, dset, opt);
    } else if (opt & OPT_S) {
	/* --sequence */
	opt &= ~OPT_S;
	err = gnuplot(list, literal, dset, opt | OPT_O | OPT_T);
    } else if (opt & OPT_D) {
	/* --grid */
	opt &= ~OPT_D;
	err = panel_grid_ts_plot(vnum, dset, opt);
    } else if (opt & OPT_A) {
	/* --stack */
	opt &= ~OPT_A;
	err = panel_grid_ts_plot(vnum, dset, opt | OPT_S | OPT_V);
    } else if (opt & OPT_B) {
	/* --boxplots */
	opt &= ~OPT_B;
	err = boxplots(list, literal, dset, opt | OPT_P);
    } else if (opt & OPT_C) {
	/* --boxplot */
	opt &= ~OPT_C;
	err = boxplots(list, literal, dset, opt);
    }

    return err;
}

static int data_straddle_zero (const gretl_matrix *m)
{
    int t, lt0 = 0, gt0 = 0;

    for (t=0; t<m->rows; t++) {
	if (gretl_matrix_get(m, t, 1) < 0) {
	    lt0 = 1;
	}
	if (gretl_matrix_get(m, t, 2) > 0) {
	    gt0 = 1;
	}
	if (lt0 && gt0) {
	    return 1;
	}
    }

    return 0;
}

static void real_irf_print_plot (const gretl_matrix *resp,
				 const char *targname,
				 const char *shockname,
				 const char *perlabel,
				 double alpha,
				 int confint,
				 int use_fill,
				 FILE *fp)
{
    int periods = gretl_matrix_rows(resp);
    gchar *title = NULL;
    int t;

    if (!confint) {
	fputs("# impulse response plot\n", fp);
    }

    if (confint) {
	fputs("set key left top\n", fp);
	title = g_strdup_printf(_("response of %s to a shock in %s, "
				  "with bootstrap confidence interval"),
				targname, shockname);
    } else {
	fputs("set nokey\n", fp);
	title = g_strdup_printf(_("response of %s to a shock in %s"),
				targname, shockname);
    }

    fprintf(fp, "set xlabel '%s'\n", perlabel);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [-1:%d]\n", periods);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);

    if (confint) {
	double ql = alpha / 2;
	double qh = 1.0 - ql;

	fputs("plot \\\n", fp);
	if (use_fill) {
	    title = g_strdup_printf(_("%g percent confidence band"), 100 * (1 - alpha));
	    print_filledcurve_line(title, NULL, fp);
	    g_free(title);
	    if (data_straddle_zero(resp)) {
		fputs("0 notitle w lines lt 0, \\\n", fp);
	    }
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1\n", _("point estimate"));
	} else {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n",
		    _("point estimate"));
	    title = g_strdup_printf(_("%g and %g quantiles"), ql, qh);
	    fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n", title);
	    g_free(title);
	}
    } else {
	fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    if (confint && use_fill) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.10g %.10g\n", t,
		    gretl_matrix_get(resp, t, 1),
		    gretl_matrix_get(resp, t, 2));
	}
	fputs("e\n", fp);
    }

    for (t=0; t<periods; t++) {
	fprintf(fp, "%d %.10g\n", t, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint && !use_fill) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.10g %.10g %.10g\n", t,
		    gretl_matrix_get(resp, t, 0),
		    gretl_matrix_get(resp, t, 1),
		    gretl_matrix_get(resp, t, 2));
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();
}

int
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock,
				 int periods, double alpha,
				 const DATASET *dset,
				 gretlopt opt)
{
    int use_fill = !(opt & OPT_E);
    gretl_matrix *resp;
    int err = 0;

    if (alpha != 0 && (alpha < 0.01 || alpha > 0.5)) {
	return E_DATA;
    }

    resp = gretl_VAR_get_impulse_response(var, targ, shock, periods,
					  alpha, dset, &err);

    if (!err) {
	int vtarg = gretl_VAR_get_variable_number(var, targ);
	int vshock = gretl_VAR_get_variable_number(var, shock);
	int confint = (resp->cols > 1);
	FILE *fp;

	fp = open_plot_input_file((confint)? PLOT_IRFBOOT : PLOT_REGULAR, 0, &err);
	if (!err) {
	    real_irf_print_plot(resp, dset->varname[vtarg],
				dset->varname[vshock],
				dataset_period_label(dset),
				alpha, confint, use_fill,
				fp);
	    err = finalize_plot_input_file(fp);
	}
	gretl_matrix_free(resp);
    }

    return err;
}

int gretl_VAR_plot_FEVD (GRETL_VAR *var, int targ, int periods,
			 const DATASET *dset, gretlopt opt)
{
    FILE *fp = NULL;
    gretl_matrix *V;
    gchar *title;
    int i, t, v, histo;
    PlotType ptype;
    int err = 0;

    V = gretl_VAR_get_FEVD_matrix(var, targ, -1, periods, dset, &err);
    if (V == NULL) {
	return E_ALLOC;
    }

    histo = (opt & OPT_H)? 1 : 0;
    ptype = histo ? PLOT_STACKED_BAR : PLOT_REGULAR;

    fp = open_plot_input_file(ptype, 0, &err);
    if (err) {
	gretl_matrix_free(V);
	return err;
    }

    v = gretl_VAR_get_variable_number(var, targ);

    fprintf(fp, "set xlabel '%s'\n", dataset_period_label(dset));
    title = g_strdup_printf(_("forecast variance decomposition for %s"),
			    dset->varname[v]);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);

    if (histo) {
	fputs("set key outside\n", fp);
	fputs("# literal lines = 3\n", fp);
	fputs("set style fill solid 0.35\n", fp);
	fputs("set style histogram rowstacked\n", fp);
	fputs("set style data histogram\n", fp);
	fprintf(fp, "set xrange [-1:%d]\n", periods);
    } else {
	fputs("set key left top\n", fp);
	fputs("set xzeroaxis\n", fp);
    }

    fputs("set yrange [0:100]\n", fp);
    fputs("plot \\\n", fp);

    for (i=0; i<var->neqns; i++) {
	v = gretl_VAR_get_variable_number(var, i);
	if (histo) {
	    fprintf(fp, "'-' using 2 title \"%s\"", dset->varname[v]);
	} else {
	    fprintf(fp, "'-' using 1:2 title \"%s\" w lines", dset->varname[v]);
	}
	if (i < var->neqns - 1) {
	    fputs(", \\\n", fp);
	} else {
	    fputc('\n', fp);
	}
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<var->neqns; i++) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.4f\n", t, 100 * gretl_matrix_get(V, t, i));
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    gretl_matrix_free(V);

    return finalize_plot_input_file(fp);
}

#define NEW_IRF 1

int gretl_VAR_plot_multiple_irf (GRETL_VAR *var,
				 int periods, double alpha,
				 const DATASET *dset,
				 gretlopt opt)
{
    FILE *fp = NULL;
    GptFlags flags = 0;
    int confint = 0;
    int use_fill = !(opt & OPT_E);
    gchar *title = NULL;
    int n = var->neqns;
    int nplots = n * n;
    int vtarg, vshock;
#if NEW_IRF
    gretl_matrix *R = NULL;
    int Rcol, Rstep;
#endif
    int t, i, j;
    int err = 0;

    maybe_set_small_font(nplots);

    if (nplots > 12) {
	flags |= GPT_XXL;
    } else if (nplots > 9) {
	flags |= GPT_XL;
    }

    fp = open_plot_input_file(PLOT_MULTI_IRF, flags, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "set multiplot layout %d,%d\n", n, n);

    if (n < 4) {
	fprintf(fp, "set xlabel '%s'\n", dataset_period_label(dset));
    } else {
	fputs("set noxlabel\n", fp);
    }

    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [-1:%d]\n", periods);

    gretl_push_c_numeric_locale();

    /* Use facility to get all impulse responses
       via one call
    */
#if NEW_IRF
    R = gretl_VAR_get_impulse_response(var, -1, -1, periods,
				       alpha, dset, &err);
    if (!err && R->cols > nplots) {
	confint = 1;
    }
    Rcol = 0;
    Rstep = confint ? 3 : 1;

    for (i=0; i<n && !err; i++) {
	vtarg = gretl_VAR_get_variable_number(var, i);

	for (j=0; j<n; j++) {
	    vshock = gretl_VAR_get_variable_number(var, j);

	    if (i == 0 && j == 0) {
		/* the first plot */
		if (confint) {
		    fputs("set key left top\n", fp);
		} else {
		    fputs("set nokey\n", fp);
		}
	    }
	    title = g_strdup_printf("%s -> %s", dset->varname[vshock],
				    dset->varname[vtarg]);
	    fprintf(fp, "set title '%s'\n", title);
	    g_free(title);

	    fputs("plot \\\n", fp);
	    if (confint && use_fill) {
		print_filledcurve_line(NULL, NULL, fp);
		fputs("'-' using 1:2 notitle w lines lt 1\n", fp);
	    } else if (confint) {
		fputs("'-' using 1:2 notitle w lines, \\\n", fp);
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("'-' using 1:2 notitle w lines\n", fp);
	    }

	    if (confint && use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g\n", t,
			    gretl_matrix_get(R, t, Rcol+1),
			    gretl_matrix_get(R, t, Rcol+2));
		}
		fputs("e\n", fp);
	    }

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.10g\n", t, gretl_matrix_get(R, t, Rcol));
	    }
	    fputs("e\n", fp);

	    if (confint && !use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g %.10g\n", t,
			    gretl_matrix_get(R, t, Rcol),
			    gretl_matrix_get(R, t, Rcol+1),
			    gretl_matrix_get(R, t, Rcol+2));
		}
		fputs("e\n", fp);
	    }
	    Rcol += Rstep;
	}
    }
    gretl_matrix_free(R);
#else /* old IRF method */
    for (i=0; i<n && !err; i++) {
	vtarg = gretl_VAR_get_variable_number(var, i);

	for (j=0; j<n; j++) {
	    gretl_matrix *resp;

	    resp = gretl_VAR_get_impulse_response(var, i, j, periods,
						  alpha, dset, &err);
	    if (err) {
		break;
	    }

	    if (i == 0 && j == 0) {
		/* the first plot */
		if (gretl_matrix_cols(resp) > 1) {
		    confint = 1;
		    fputs("set key left top\n", fp);
		} else {
		    fputs("set nokey\n", fp);
		}
	    }

	    vshock = gretl_VAR_get_variable_number(var, j);
	    fprintf(fp, "set title '%s -> %s'\n", dset->varname[vshock],
		    dset->varname[vtarg]);

	    fputs("plot \\\n", fp);

	    if (confint && use_fill) {
		print_filledcurve_line(NULL, NULL, fp);
		fputs("'-' using 1:2 notitle w lines lt 1\n", fp);
	    } else if (confint) {
		fputs("'-' using 1:2 notitle w lines, \\\n", fp);
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("'-' using 1:2 notitle w lines\n", fp);
	    }

	    if (confint && use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g\n", t,
			    gretl_matrix_get(resp, t, 1),
			    gretl_matrix_get(resp, t, 2));
		}
		fputs("e\n", fp);
	    }

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.10g\n", t, gretl_matrix_get(resp, t, 0));
	    }
	    fputs("e\n", fp);

	    if (confint && !use_fill) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.10g %.10g %.10g\n", t,
			    gretl_matrix_get(resp, t, 0),
			    gretl_matrix_get(resp, t, 1),
			    gretl_matrix_get(resp, t, 2));
		}
		fputs("e\n", fp);
	    }

	    gretl_matrix_free(resp);
	}
    }
#endif /* NEW_IRF or not */

    gretl_pop_c_numeric_locale();

    if (err) {
	fclose(fp);
	return err;
    }

    fputs("unset multiplot\n", fp);

    return finalize_plot_input_file(fp);
}

int gretl_system_residual_plot (void *p, int ci, int eqn, const DATASET *dset)
{
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *E = NULL;
    FILE *fp = NULL;
    const double *obs;
    char lwstr[8];
    int single = 0;
    int nvars, nobs;
    int i, v, t, t1;
    int imin, imax;
    int err = 0;

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) p;
	E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	E = sys->E;
    }

    if (E == NULL) {
	return E_DATA;
    }

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);
    t1 = gretl_matrix_get_t1(E);

    if (eqn > 0 && eqn <= nvars) {
	imin = eqn - 1;
	imax = imin + 1;
	single = 1;
    } else {
	imin = 0;
	imax = nvars;
	single = (nvars == 1);
    }

    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	return err;
    }

    obs = gretl_plotx(dset, OPT_NONE);

    if (quarterly_or_monthly(dset)) {
	fprintf(fp, "# timeseries %d\n", dset->pd);
    }

    if (!single) {
	fputs("# system residual plot\n", fp);
    }

    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (ci == VAR) {
	fprintf(fp, "set title '%s'\n", _("VAR residuals"));
    } else {
	fprintf(fp, "set title '%s'\n", _("System residuals"));
    }

    set_lwstr(lwstr);

    if (single) {
	fputs("plot ", fp);
    } else {
	fputs("plot \\\n", fp);
    }

    for (i=imin; i<imax; i++) {
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines%s",
		dset->varname[v], lwstr);
	if (i == imax - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp);
	}
    }

    gretl_push_c_numeric_locale();

    for (i=imin; i<imax; i++) {
	for (t=0; t<nobs; t++) {
	    double eti = gretl_matrix_get(E, t, i);

	    if (obs != NULL) {
		fprintf(fp, "%g %.10g\n", obs[t+t1], eti);
	    } else {
		fprintf(fp, "%d %.10g\n", t+1, eti);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int gretl_VECM_combined_EC_plot (GRETL_VAR *var,
				 const DATASET *dset)
{
    const gretl_matrix *EC = NULL;
    FILE *fp = NULL;
    const double *obs;
    int nvars, nobs;
    int i, t, t1;
    int err = 0;

    EC = VECM_get_EC_matrix(var, dset, &err);
    if (err) {
	return err;
    }

    t1 = gretl_matrix_get_t1(EC);

    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	return err;
    }

    obs = gretl_plotx(dset, OPT_NONE);

    nvars = gretl_matrix_cols(EC);
    nobs = gretl_matrix_rows(EC);

    fputs("# VECM EC plot\n", fp);
    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (nvars > 1) {
	fprintf(fp, "set title '%s'\n", _("EC terms"));
    } else {
	fprintf(fp, "set title '%s'\n", _("EC term"));
    }

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
	if (nvars > 1) {
	    fprintf(fp, "'-' using 1:2 title 'EC %d' w lines", i + 1);
	} else {
	    fprintf(fp, "'-' using 1:2 notitle w lines");
	}
	if (i == nvars - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp);
	}
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<nvars; i++) {
	for (t=0; t<nobs; t++) {
	    double eti = gretl_matrix_get(EC, t, i);

	    if (obs != NULL) {
		fprintf(fp, "%g %.10g\n", obs[t+t1], eti);
	    } else {
		fprintf(fp, "%d %.10g\n", t+1, eti);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int gretl_system_residual_mplot (void *p, int ci, const DATASET *dset)
{
    const gretl_matrix *E = NULL;
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    FILE *fp = NULL;
    const double *obs;
    double startdate;
    double xmin, xmax, xrange;
    int nvars, nobs, incr;
    int rows, cols;
    int i, v, t, t1;
    int err = 0;

    if (ci == VAR || ci == VECM) {
	var = (GRETL_VAR *) p;
	E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == SYSTEM) {
	sys = (equation_system *) p;
	E = sys->E;
    }

    if (E == NULL) {
	return E_DATA;
    }

    nvars = gretl_matrix_cols(E);
    if (nvars > 6) {
	return 1;
    }

    obs = gretl_plotx(dset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    nobs = gretl_matrix_rows(E);
    t1 = gretl_matrix_get_t1(E);

    fp = open_plot_input_file(PLOT_MULTI_BASIC, 0, &err);
    if (err) {
	return err;
    }

    if (nvars <= 3) {
	rows = nvars;
	cols = 1;
    } else if (nvars == 4) {
	rows = cols = 2;
    } else {
	rows = 3;
	cols = ceil(nvars / 3);
    }

    fprintf(fp, "set multiplot layout %d,%d\n", rows, cols);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set noxlabel\n", fp);
    fputs("set noylabel\n", fp);

    gretl_push_c_numeric_locale();

    startdate = obs[t1];
    incr = nobs / (2 * dset->pd);
    if (incr > 0) {
	fprintf(fp, "set xtics %g, %d\n", ceil(startdate), incr);
    }

    gretl_minmax(t1, t1 + nobs - 1, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.10g:%.10g]\n", xmin, xmax);

    for (i=0; i<nvars; i++) {
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}

	fprintf(fp, "set title '%s'\n", dset->varname[v]);
	fputs("plot '-' using 1:2 with lines\n", fp);

	for (t=0; t<nobs; t++) {
	    double eti;

	    fprintf(fp, "%.10g\t", obs[t+t1]);
	    eti = gretl_matrix_get(E, t, i);
	    write_gp_dataval(eti, fp, 1);
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);

    return finalize_plot_input_file(fp);
}

int gretl_VAR_roots_plot (GRETL_VAR *var)
{
    const gretl_matrix *lam;
    FILE *fp = NULL;
    double x, y;
    double px, py;
    int i, n, err = 0;

    lam = gretl_VAR_get_roots(var, &err);
    if (err) {
	return err;
    }

    fp = open_plot_input_file(PLOT_ROOTS, 0, &err);
    if (err) {
	return err;
    }

    n = gretl_matrix_rows(lam);

    fprintf(fp, "set title '%s'\n",
	    _("VAR inverse roots in relation to the unit circle"));
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 w lines, \\\n'-' using 1:2 w points pt 7\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(lam, i, 0);
        y = gretl_matrix_get(lam, i, 1);
	/* in polar form */
	px = atan2(y, x);
	py = sqrt(x * x + y * y);
	fprintf(fp, "%.8f %.8f # %.4f,%.4f\n", px, py, x, y);
    }

    gretl_pop_c_numeric_locale();

    fputs("e\n", fp);

    return finalize_plot_input_file(fp);
}

/**
 * confidence_ellipse_plot:
 * @V: 2x2 covariance matrix.
 * @b: 2-vector containing point estimates
 * @tcrit: critical t-value for 1 - alpha confidence.
 * @Fcrit: critical F-value for 1 - alpha confidence.
 * @alpha: nominal non-coverage, as decimal.
 * @iname: name of first parameter.
 * @jname: name of second parameter.
 *
 * Plots a 95% confidence ellipse for the parameter estimates
 * in @b with covariance @V.
 *
 * Returns: 0 on success, non-zero on error.
 */

int confidence_ellipse_plot (gretl_matrix *V, double *b,
			     double tcrit, double Fcrit, double alpha,
			     const char *iname, const char *jname)
{
    FILE *fp = NULL;
    double maxerr[2];
    double xcoeff[2];
    double ycoeff[2];
    double cval = 100 * (1 - alpha);
    gretl_matrix *e = NULL;
    gchar *title;
    int i, err = 0;

    maxerr[0] = tcrit * sqrt(gretl_matrix_get(V, 0, 0));
    maxerr[1] = tcrit * sqrt(gretl_matrix_get(V, 1, 1));

    err = gretl_invert_symmetric_matrix(V);
    if (err) {
	return err;
    }

    e = gretl_symmetric_matrix_eigenvals(V, 1, &err);
    if (err) {
	return err;
    }

    for (i=0; i<2; i++) {
	e->val[i] = sqrt(1.0 / e->val[i] * Fcrit);
	xcoeff[i] = e->val[i] * gretl_matrix_get(V, 0, i);
	ycoeff[i] = e->val[i] * gretl_matrix_get(V, 1, i);
    }

    gretl_matrix_free(e);

    fp = open_plot_input_file(PLOT_ELLIPSE, 0, &err);
    if (err) {
	return err;
    }

    title = g_strdup_printf(_("%g%% confidence ellipse and %g%% marginal intervals"),
			    cval, cval);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);

    fputs("# literal lines = 9\n", fp);
    fputs("set parametric\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);

    fprintf(fp, "set xlabel '%s'\n", iname);
    fprintf(fp, "set ylabel '%s'\n", jname);
    fprintf(fp, "set label '%.3g, %.3g' at ", b[0], b[1]);

    gretl_push_c_numeric_locale();

    fprintf(fp, "%g,%g point lt 2 pt 1 offset 3,3\n", b[0], b[1]);

    fprintf(fp, "x(t) = %g*cos(t)%+g*sin(t)%+g\n", xcoeff[0], xcoeff[1], b[0]);
    fprintf(fp, "y(t) = %g*cos(t)%+g*sin(t)%+g\n", ycoeff[0], ycoeff[1], b[1]);

    fputs("plot x(t), y(t) notitle, \\\n", fp);
    fprintf(fp, "%g, y(t) notitle w lines lt 2, \\\n", b[0] - maxerr[0]);
    fprintf(fp, "%g, y(t) notitle w lines lt 2, \\\n", b[0] + maxerr[0]);
    fprintf(fp, "x(t), %g notitle w lines lt 2, \\\n", b[1] - maxerr[1]);
    fprintf(fp, "x(t), %g notitle w lines lt 2\n", b[1] + maxerr[1]);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static void corrgm_min_max (const double *acf, const double *pacf,
			    int m, double pm, double *ymin, double *ymax)
{
    int k;

    /* the range should include the plus/minus bands, but
       should not go outside (-1, 1) */
    *ymax = pm * 1.2;
    if (*ymax > 1) *ymax = 1;
    *ymin = -pm * 1.2;
    if (*ymin < -1) *ymin = -1;

    /* adjust based on min and max of ACF, PACF */
    for (k=0; k<m; k++) {
	if (acf[k] > *ymax) {
	    *ymax = acf[k];
	} else if (acf[k] < *ymin) {
	    *ymin = acf[k];
	}
	if (pacf[k] > *ymax) {
	    *ymax = pacf[k];
	} else if (pacf[k] < *ymin) {
	    *ymin = pacf[k];
	}
    }

    if (*ymax > 0.5) {
	*ymax = 1;
    } else {
	*ymax *= 1.2;
    }

    if (*ymin < -0.5) {
	*ymin = -1;
    } else {
	*ymin *= 1.2;
    }

    /* make the range symmetrical */
    if (fabs(*ymin) > *ymax) {
	*ymax = -*ymin;
    } else if (*ymax > fabs(*ymin)) {
	*ymin = -*ymax;
    }
}

static int real_correlogram_print_plot (const char *vname,
					const double *acf,
					const double *pacf,
					const gretl_matrix *PM,
					int m, double pm,
					gretlopt opt,
					FILE *fp)
{
    /* xgettext:no-c-format */
    const char *PM_title = N_("95% interval");
    char pm_title[16];
    double ymin, ymax;
    int k;

    sprintf(pm_title, "%.2f/T^%.1f", 1.96, 0.5);

    corrgm_min_max(acf, pacf, m, pm, &ymin, &ymax);

    gretl_push_c_numeric_locale();

    /* create two separate plots, if both are OK */
    if (pacf != NULL) {
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    }
    fputs("set xzeroaxis\n", fp);
    print_keypos_string(GP_KEY_RIGHT_TOP, fp);
    fprintf(fp, "set xlabel '%s'\n", _("lag"));

    fprintf(fp, "set yrange [%.2f:%.2f]\n", ymin, ymax);

    /* upper plot: Autocorrelation Function or ACF */
    if (pacf != NULL) {
	fputs("set origin 0.0,0.50\n", fp);
    }
    if (opt & OPT_R) {
	fprintf(fp, "set title '%s'\n", _("Residual ACF"));
    } else {
	fprintf(fp, "set title '%s %s'\n", _("ACF for"), vname);
    }
    fprintf(fp, "set xrange [0:%d]\n", m + 1);
    if (PM != NULL) {
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 notitle w impulses lw 5, \\\n"
		"'-' title '%s' w lines lt 2, \\\n"
		"'-' notitle w lines lt 2\n", _(PM_title));
    } else {
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 notitle w impulses lw 5, \\\n"
		"%g title '+- %s' lt 2, \\\n"
		"%g notitle lt 2\n", pm, pm_title, -pm);
    }
    for (k=0; k<m; k++) {
	fprintf(fp, "%d %g\n", k + 1, acf[k]);
    }
    fputs("e\n", fp);
    if (PM != NULL) {
	/* Bartlett-type confidence band data */
	for (k=0; k<m; k++) {
	    fprintf(fp, "%d %g\n", k + 1, gretl_matrix_get(PM, k, 1));
	}
	fputs("e\n", fp);
	for (k=0; k<m; k++) {
	    fprintf(fp, "%d -%g\n", k + 1, gretl_matrix_get(PM, k, 1));
	}
	fputs("e\n", fp);
    }

    if (pacf != NULL) {
	/* lower plot: Partial Autocorrelation Function or PACF */
	fputs("set origin 0.0,0.0\n", fp);
	if (opt & OPT_R) {
	    fprintf(fp, "set title '%s'\n", _("Residual PACF"));
	} else {
	    fprintf(fp, "set title '%s %s'\n", _("PACF for"), vname);
	}
	fprintf(fp, "set xrange [0:%d]\n", m + 1);
	fprintf(fp, "plot \\\n"
		"'-' using 1:2 notitle w impulses lw 5, \\\n"
		"%g title '+- %s' lt 2, \\\n"
		"%g notitle lt 2\n", pm, pm_title, -pm);
	for (k=0; k<m; k++) {
	    fprintf(fp, "%d %g\n", k + 1, pacf[k]);
	}
	fputs("e\n", fp);
    }

    if (pacf != NULL) {
	fputs("unset multiplot\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

int correlogram_plot (const char *vname,
		      const double *acf,
		      const double *pacf,
		      const gretl_matrix *PM,
		      int m, double pm,
		      gretlopt opt)
{
    FILE *fp;
    int err = 0;

    fp = open_plot_input_file(PLOT_CORRELOGRAM, 0, &err);

    if (!err) {
	real_correlogram_print_plot(vname, acf, pacf,
				    PM, m, pm, opt, fp);
	err = finalize_plot_input_file(fp);
    }

    return err;
}

static int roundup_mod (int i, double x)
{
    return (int) ceil((double) x * i);
}

/* options: OPT_R use radians as unit
   OPT_D use degrees as unit
   OPT_L use log scale
*/

static int real_pergm_plot (const char *vname,
			    int T, int L,
			    const double *dens,
			    gretlopt opt,
			    FILE *fp)
{
    char s[80];
    double ft;
    int T2 = T / 2;
    int k, t, err = 0;

    fputs("set xtics nomirror\n", fp);

    fprintf(fp, "set x2label '%s'\n", _("periods"));
    fprintf(fp, "set x2range [0:%d]\n", roundup_mod(T, 2.0));

    fputs("set x2tics (", fp);
    k = T2 / 6;
    for (t = 1; t <= T2; t += k) {
	fprintf(fp, "\"%.1f\" %d, ", (double) T / t, 4 * t);
    }
    fprintf(fp, "\"\" %d)\n", 2 * T);

    if (opt & OPT_R) {
	fprintf(fp, "set xlabel '%s'\n", _("radians"));
    } else if (opt & OPT_D) {
	fprintf(fp, "set xlabel '%s'\n", _("degrees"));
    } else {
	fprintf(fp, "set xlabel '%s'\n", _("scaled frequency"));
    }

    fputs("set xzeroaxis\n", fp);
    fputs("set nokey\n", fp);

    /* open gnuplot title string */
    fputs("set title '", fp);

    if (vname == NULL) {
	fputs(_("Residual spectrum"), fp);
    } else {
	sprintf(s, _("Spectrum of %s"), vname);
	fputs(s, fp);
    }

    if (opt & OPT_O) {
	fputs(" (", fp);
	fprintf(fp, _("Bartlett window, length %d"), L);
	fputc(')', fp);
    }

    if (opt & OPT_L) {
	fputs(" (", fp);
	fputs(_("log scale"), fp);
	fputc(')', fp);
    }

    /* close gnuplot title string */
    fputs("'\n", fp);

    gretl_push_c_numeric_locale();

    if (opt & OPT_R) {
	/* frequency scale in radians */
	fputs("set xrange [0:3.1416]\n", fp);
    } else if (opt & OPT_D) {
	/* frequency scale in degrees */
	fputs("set xrange [0:180]\n", fp);
    } else {
	/* data-scaled frequency */
	fprintf(fp, "set xrange [0:%d]\n", roundup_mod(T, 0.5));
    }

    if (!(opt & OPT_L)) {
	fprintf(fp, "set yrange [0:%g]\n", 1.2 * gretl_max(0, T/2 - 1, dens));
    }

    if (opt & OPT_R) {
	fputs("set xtics (\"0\" 0, \"π/4\" pi/4, \"π/2\" pi/2, "
	      "\"3π/4\" 3*pi/4, \"π\" pi)\n", fp);
    }

    fputs("plot '-' using 1:2 w lines\n", fp);

    for (t=1; t<=T2; t++) {
	if (opt & OPT_R) {
	    ft = M_PI * (double) t / T2;
	} else if (opt & OPT_D) {
	    ft = 180 * (double) t / T2;
	} else {
	    ft = t;
	}
	fprintf(fp, "%g %g\n", ft, (opt & OPT_L)? log(dens[t-1]) : dens[t-1]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return err;
}

int periodogram_plot (const char *vname,
		      int T, int L, const double *x,
		      gretlopt opt)
{
    FILE *fp;
    int err = 0;

    fp = open_plot_input_file(PLOT_PERIODOGRAM, 0, &err);

    if (!err) {
	real_pergm_plot(vname, T, L, x, opt, fp);
	err = finalize_plot_input_file(fp);
    }

    return err;
}

int arma_spectrum_plot (MODEL *pmod, const DATASET *dset,
			gretlopt opt)
{
    gretl_matrix *pdata = NULL;
    FILE *fp;
    int err = 0;

    pdata = arma_spectrum_plot_data(pmod, dset, &err);
    if (err) {
	return err;
    }

    fp = open_plot_input_file(PLOT_PERIODOGRAM, 0, &err);

    if (!err) {
	double px, pRe, pIm, scale = pmod->nobs * M_2PI;
	int i, grid = pdata->rows;

	gretl_push_c_numeric_locale();

	fprintf(fp, "set xrange [0:%g]\n", M_PI);
	switch (dset->pd) {
	case 12:
	    fputs("set xtics (\"0\" 0, \"π/6\" pi/6, "
		  "\"π/3\" pi/3, \"π/2\" pi/2, \"2π/3\" 2*pi/3, "
		  "\"5π/6\" 5*pi/6, \"π\" pi)\n", fp);
	    break;
	case 6:
	    fputs("set xtics (\"0\" 0, \"π/3\" pi/3, "
		  "\"2π/3\" 2*pi/3, \"π\" pi)\n", fp);
	    break;
	case 5:
	    fputs("set xtics (\"0\" 0, \"π/5\" pi/5, "
		  "\"2π/5\" 2*pi/5, \"3π/5\" 3*pi/5, "
		  "\"4π/5\" 4*pi/5, \"π\" pi)\n", fp);
	    break;
	default:
	    fputs("set xtics (\"0\" 0, \"π/4\" pi/4, \"π/2\" pi/2, "
		  "\"3π/4\" 3*pi/4, \"π\" pi)\n", fp);
	}
	fprintf(fp, "set title \"%s (%s)\"\n", _("Sample periodogram vs ARMA Spectrum"),
		_("log scale"));
	fprintf(fp, "plot '-' using 1:2 w lines title '%s' lw 2, \\\n", _("spectrum"));
	fprintf(fp, "'-' using 1:2 w lines title '%s' lw 0.5\n", _("periodogram"));

	for (i=0; i<grid; i++) {
	    fprintf(fp, "%7.5f %12.7f\n", gretl_matrix_get(pdata, i, 0),
		    log(gretl_matrix_get(pdata, i, 1)));
	}
	fputs("e\n", fp);

	for (i=0; i<grid; i++) {
	    pRe = gretl_matrix_get(pdata, i, 2);
	    pIm = gretl_matrix_get(pdata, i, 3);
	    px = (pRe * pRe + pIm * pIm) / scale;
	    fprintf(fp, "%7.5f %12.7f\n", gretl_matrix_get(pdata, i, 0), log(px));
	}
	fputs("e\n", fp);

	gretl_pop_c_numeric_locale();
	err = finalize_plot_input_file(fp);
    }

    return err;
}

#define MAKKONEN_POS 0

/* Probability of non-exceedance of the kth value in a set of n
   rank-ordered values.  See L. Makkonen, 'Bringing Closure to the
   Plotting Position Controversy', Communications in Statistics -
   Theory and Methods, vol 37, January 2008, for an argument in favor
   of using k / (n + 1); but also see many uses of (k - 1/2) / n in
   the literature.
*/

static double plotpos (int k, int n)
{
#if MAKKONEN_POS
    return k / (n + 1.0);
#else
    return (k - 0.5) / n;
#endif
}

static double quantile_interp (const double *y, int n,
			       double ftarg)
{
    double f, ret = NADBL;
    int i;

    for (i=0; i<n; i++) {
	f = plotpos(i+1, n);
	if (f >= ftarg) {
	    if (f > ftarg && i > 0) {
		double f0 = plotpos(i, n);
		double d = (ftarg - f0) / (f - f0);

		ret = (1-d) * y[i-1] + d * y[i];
	    } else {
		ret = y[i];
	    }
	    break;
	}
    }

    return ret;
}

static int qq_plot_two_series (const int *list,
			       const DATASET *dset)
{
    double *x = NULL;
    double *y = NULL;
    double f, qx, qy;
    FILE *fp = NULL;
    int vx = list[1];
    int vy = list[2];
    int nx = 10, ny = 10;
    int i, n, err = 0;

    x = gretl_sorted_series(vx, dset, OPT_NONE, &nx, &err);

    if (!err) {
	y = gretl_sorted_series(vy, dset, OPT_NONE, &ny, &err);
	if (err) {
	    free(x);
	    x = NULL;
	}
    }

    if (!err) {
	/* take the smaller sample as basis */
	n = (nx > ny)? ny : nx;
    }

    if (!err) {
	fp = open_plot_input_file(PLOT_QQ, 0, &err);
    }

    if (err) {
	free(x);
	free(y);
	return err;
    }

    fprintf(fp, "set title \"%s\"\n", _("Q-Q plot"));
    gnuplot_missval_string(fp);
    fputs("set key top left\n", fp);
    fprintf(fp, "set xlabel \"%s\"\n", series_get_graph_name(dset, vx));
    fprintf(fp, "set ylabel \"%s\"\n", series_get_graph_name(dset, vy));
    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    fputs(" x notitle w lines\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	f = plotpos(i+1, n);

	if (nx == ny) {
	    qx = x[i];
	    qy = y[i];
	} else if (nx == n) {
	    qx = x[i];
	    qy = quantile_interp(y, ny, f);
	} else {
	    qx = quantile_interp(x, nx, f);
	    qy = y[i];
	}

	if (!na(qx) && !na(qy)) {
	    fprintf(fp, "%.12g %.12g\n", qx, qy);
	}
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    free(x);
    free(y);

    return finalize_plot_input_file(fp);
}

static int normal_qq_plot (const int *list,
			   const DATASET *dset,
			   gretlopt opt)
{
    GptFlags flags = 0;
    gchar *title = NULL;
    int zscores = 0;
    double ym = 0, ys = 1;
    double p, qx, qy;
    double *y = NULL;
    FILE *fp = NULL;
    int v = list[1];
    int i, n = 20;
    int err = 0;

    y = gretl_sorted_series(v, dset, OPT_NONE, &n, &err);

    if (!err && y[0] == y[n-1]) {
	gretl_errmsg_sprintf(_("%s is a constant"), dset->varname[v]);
	err = E_DATA;
    }

    if (err) {
	return err;
    }

    if (opt & OPT_Z) {
	/* standardize the data */
	zscores = 1;
    }

    if (!(opt & OPT_R)) {
	ym = gretl_mean(0, n-1, y);
	ys = gretl_stddev(0, n-1, y);

	if (zscores) {
	    /* standardize y */
	    for (i=0; i<n; i++) {
		y[i] = (y[i] - ym) / ys;
	    }
	}
    }

    if (opt & OPT_G) {
	flags = GPT_ICON;
    }

    fp = open_plot_input_file(PLOT_QQ, flags, &err);
    if (err) {
	free(y);
	return err;
    }

    title = g_strdup_printf(_("Q-Q plot for %s"), series_get_graph_name(dset, v));
    fprintf(fp, "set title \"%s\"\n", title);
    g_free(title);
    gnuplot_missval_string(fp);
    fprintf(fp, "set xlabel \"%s\"\n", _("Normal quantiles"));

    if (opt & OPT_R) {
	fputs("set nokey\n", fp);
	fputs("plot \\\n", fp);
	fputs(" '-' using 1:2 notitle w points\n", fp);
    } else {
	fputs("set key top left\n", fp);
	fputs("plot \\\n", fp);
	fputs(" '-' using 1:2 notitle w points, \\\n", fp);
	fputs(" x title \"y = x\" w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	p = plotpos(i+1, n);
	/* empirical quantile */
	qy = y[i];
	/* normal quantile */
	qx = normal_critval(1 - p);
	if (!na(qx) && !zscores && !(opt & OPT_R)) {
	    qx = ys * qx + ym;
	}
	if (!na(qx) && !na(qy)) {
	    fprintf(fp, "%.12g %.12g\n", qx, qy);
	}
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    free(y);

    return finalize_plot_input_file(fp);
}

int qq_plot (const int *list, const DATASET *dset, gretlopt opt)
{
    int err;

    if (list[0] == 1) {
	/* one series against normal */
	err = normal_qq_plot(list, dset, opt);
    } else if (list[0] == 2) {
	/* two empirical series */
	err = qq_plot_two_series(list, dset);
    } else {
	err = E_DATA;
    }

    return err;
}

int kd_plot (const int *list, const DATASET *dset, gretlopt opt)
{
    int (*kdfunc) (const double *y, int, double,
		   const char *, gretlopt);
    int err = 0;

    kdfunc = get_plugin_function("kernel_density");

    if (kdfunc == NULL) {
	err = E_FOPEN;
    } else {
	int v = list[1];
	const double *y = dset->Z[v] + dset->t1;
	const char *label = dset->varname[v];
	int n = sample_size(dset);
	double bws = 1.0;

	if (opt & OPT_S) {
	    /* bandwidth scale */
	    bws = get_optval_double(KDPLOT, OPT_S, &err);
	}
	if (!err) {
	    err = (*kdfunc)(y, n, bws, label, opt);
	}
    }

    return err;
}

static int pd_from_compfac (const DATASET *dset,
			    int compfac,
			    char *stobs)
{
    int pd = -1;

    if (dset->pd == 1 && (compfac == 12 || compfac == 4)) {
	/* annual from monthly or quarterly */
	pd = compfac;
    } else if (dset->pd == 4 && compfac == 3) {
	/* quarterly from monthly */
	pd = 12;
    } else if (dset->pd == 4) {
	/* maybe quarterly from daily? */
	if (compfac >= 60 && compfac <= 69) {
	    pd = 5;
	} else if (compfac >= 71 && compfac <= 81) {
	    pd = 6;
	} else if (compfac >= 82 && compfac <= 93) {
	    return 7;
	}
    } else if (dset->pd == 12) {
	/* maybe monthly from daily? */
	if (compfac >= 20 && compfac <= 23) {
	    pd = 5;
	} else if (compfac >= 24 && compfac <= 27) {
	    pd = 6;
	} else if (compfac >= 28 && compfac <= 31) {
	    pd = 7;
	}
    }

    if (pd > 0) {
	char *p, tmp[OBSLEN];
	int y, q, m;

	ntolabel(tmp, dset->t1, dset);
	y = atoi(tmp);
	p = strchr(tmp, ':');

	if ((dset->pd == 4 || dset->pd == 12) && p == NULL) {
	    return -1;
	}

	if (dset->pd == 1) {
	    if (pd == 4) {
		sprintf(stobs, "%d:1", y);
	    } else if (pd == 12) {
		sprintf(stobs, "%d:01", y);
	    }
	} else if (dset->pd == 4) {
	    q = atoi(p + 1);
	    m = (q==1)? 1 : (q==2)? 4 : (q==3)? 7 : 10;
	    if (pd == 12) {
		sprintf(stobs, "%d:%02d", y, m);
	    } else {
		/* daily */
		sprintf(stobs, "%d-%02d-01", y, m);
	    }
	} else if (dset->pd == 12) {
	    /* daily */
	    m = atoi(p + 1);
	    sprintf(stobs, "%d-%02d-01", y, m);
	}
    }

    return pd;
}

static void transcribe_graph_name (DATASET *targ, int i,
				   const DATASET *src, int j)
{
    const char *s = series_get_display_name(src, j);

    if (s != NULL && *s != '\0') {
	series_record_display_name(targ, i, s);
    }
}

/* high-frequency plot for MIDAS */

int hf_plot (const int *list, const char *literal,
	     const DATASET *dset, gretlopt opt)
{
    DATASET *hset;
    double xit;
    char stobs[OBSLEN];
    int *gplist = NULL;
    int *hflist = NULL;
    int *lflist = NULL;
    gchar *mylit = NULL;
    char *p;
    gretlopt plotopt = OPT_T;
    int plotpd = 0;
    int nv, nlf = 0;
    int cfac;
    int i, s, t, T;
    int err;

    if (list == NULL || list[0] < 3) {
	return E_INVARG;
    } else if (!dataset_is_time_series(dset)) {
	return E_PDWRONG;
    }

    if (gretl_list_has_separator(list)) {
	err = gretl_list_split_on_separator(list, &hflist, &lflist);
	if (err) {
	    return err;
	} else {
	    cfac = hflist[0];
	    nlf = lflist[0];
	    nv = 2 + nlf;
	}
    } else {
	cfac = list[0];
	nv = 2;
    }

    T = sample_size(dset) * cfac;

    hset = create_auxiliary_dataset(nv, T, OPT_NONE);
    if (hset == NULL) {
	return E_ALLOC;
    }

    /* set the hf series name */
    strcpy(hset->varname[1], dset->varname[list[1]]);
    p = strrchr(hset->varname[1], '_');
    if (p != NULL) {
	*p = '\0';
    }
    transcribe_graph_name(hset, 1, dset, list[1]);

    s = 0;
    /* transcribe high-frequency data */
    for (t=dset->t1; t<=dset->t2; t++) {
	for (i=cfac; i>0; i--) {
	    xit = dset->Z[list[i]][t];
	    hset->Z[1][s++] = xit;
	}
    }

    if (lflist != NULL) {
	/* add low-frequency term(s), if any */
	for (i=1; i<=nlf; i++) {
	    int vi = lflist[i];

	    strcpy(hset->varname[i+1], dset->varname[vi]);
	    transcribe_graph_name(hset, i+1, dset, vi);
	    for (s=0; s<hset->n; s++) {
		hset->Z[i+1][s] = NADBL;
	    }
	    s = 0;
	    for (t=dset->t1; t<=dset->t2; t++) {
		xit = dset->Z[vi][t];
		hset->Z[i+1][s] = xit;
		s += cfac;
	    }
	}
    }

    gplist = gretl_consecutive_list_new(1, nv - 1);

    if (lflist != NULL) {
	free(lflist);
	lflist = gretl_consecutive_list_new(2, nv - 1);
    }

    /* try to set a suitable time-series interpretation
       on the data to be plotted
    */
    plotpd = pd_from_compfac(dset, cfac, stobs);
    if (plotpd > 0) {
	char numstr[12];

	sprintf(numstr, "%d", plotpd);
	set_obs(numstr, stobs, hset, OPT_T);
    }

    if (opt & OPT_O) {
	plotopt |= OPT_O;
    }
    if (opt & OPT_U) {
	plotopt |= OPT_U;
    }

    if (literal == NULL) {
	const char *pdstr = midas_pdstr(dset, cfac);
	gchar *title;

	title = g_strdup_printf("%s (%s)", hset->varname[1], _(pdstr));
	mylit = g_strdup_printf("{ set ylabel ''; set title '%s'; }", title);
	g_free(title);
    }

    set_effective_plot_ci(HFPLOT);
    na_skiplist = lflist; /* file-scope global */
    err = gnuplot(gplist, literal != NULL ? literal : mylit,
		  hset, plotopt);
    na_skiplist = NULL;
    set_effective_plot_ci(GNUPLOT);

    free(gplist);
    free(hflist);
    free(lflist);
    g_free(mylit);

    destroy_dataset(hset);

    return err;
}

/**
 * xy_plot_with_control:
 * @list: list of variables by ID number: Y, X, control.
 * @literal: extra gnuplot commands or %NULL.
 * @dset: dataset struct.
 * @opt: can add "gnuplot" options.
 *
 * Constructs a scatterplot of modified Y against modified X,
 * where the modification consists in taking the residuals from
 * OLS regression of the variable in question on the control variable,
 * a la Frisch-Waugh-Lovell.
 *
 * Returns: 0 on success, non-zero on error.
 */

int xy_plot_with_control (const int *list, const char *literal,
			  const DATASET *dset, gretlopt opt)
{
    int t1 = dset->t1, t2 = dset->t2;
    int mlist[4] = {3, 0, 0, 0};
    char dname[MAXDISP];
    MODEL mod;
    DATASET *gset = NULL;
    int vy, vx, vz;
    int s, t, T;
    int missvals = 0;
    int err = 0;

    if (list == NULL || list[0] != 3) {
	return E_DATA;
    }

    vy = list[1];
    vx = list[2];
    vz = list[3];

    list_adjust_sample(list, &t1, &t2, dset, &missvals);

    /* maximum usable observations */
    T = t2 - t1 + 1 - missvals;

    if (T < 3) {
	return E_DF;
    }

    /* create temporary dataset */

    gset = create_auxiliary_dataset(4, T, 0);
    if (gset == NULL) {
	return E_ALLOC;
    }

    sprintf(dname, _("adjusted %s"), dset->varname[vy]);
    series_set_display_name(gset, 1, dname);

    sprintf(dname, _("adjusted %s"), dset->varname[vx]);
    series_set_display_name(gset, 2, dname);

    s = 0;
    for (t=t1; t<=t2; t++) {
	if (!na(dset->Z[vy][t]) && !na(dset->Z[vx][t]) && !na(dset->Z[vz][t])) {
	    gset->Z[1][s] = dset->Z[vy][t];
	    gset->Z[2][s] = dset->Z[vx][t];
	    gset->Z[3][s] = dset->Z[vz][t];
	    s++;
	}
    }

    /* regress Y (1) on Z (3) and save the residuals in series 1 */

    mlist[1] = 1;
    mlist[3] = 3;
    mod = lsq(mlist, gset, OLS, OPT_A);
    err = mod.errcode;
    if (err) {
	clear_model(&mod);
	goto bailout;
    } else {
	for (t=0; t<mod.nobs; t++) {
	    gset->Z[1][t] = mod.uhat[t];
	}
	clear_model(&mod);
    }

    /* regress X (2) on Z and save the residuals in series 2 */

    mlist[1] = 2;
    mod = lsq(mlist, gset, OLS, OPT_A);
    err = mod.errcode;
    if (err) {
	clear_model(&mod);
	goto bailout;
    } else {
	for (t=0; t<mod.nobs; t++) {
	    gset->Z[2][t] = mod.uhat[t];
	}
	clear_model(&mod);
    }

    /* call for scatter of Y-residuals against X-residuals */

    mlist[0] = 2;
    mlist[1] = 1;
    mlist[2] = 2;
    err = gnuplot(mlist, literal, gset, opt | OPT_C);

 bailout:

    /* trash the temporary dataset */
    destroy_dataset(gset);

    return err;
}

int is_auto_fit_string (const char *s)
{
    /* FIXME? */
    if (strstr(s, "automatic fit")) return 1;
    if (strstr(s, _("with least squares fit"))) return 1;
    return 0;
}

static void inbuf_inject_literal (const char *buf,
				  const char *literal,
				  FILE *fp)
{
    char *p = strstr(buf, "\nplot ");

    if (p != NULL) {
	int n = p - buf;

	fwrite(buf, 1, n, fp);
	print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fp);
	fputs(p + 1, fp);
    } else {
	fputs(buf, fp);
    }
}

static void infile_inject_literal (FILE *fin,
				   const char *literal,
				   FILE *fout)
{
    char line[1024];

    while (fgets(line, sizeof line, fin)) {
	if (!strncmp(line, "plot ", 5)) {
	    print_gnuplot_literal_lines(literal, GNUPLOT, OPT_NONE, fout);
	}
	fputs(line, fout);
    }
}

/**
 * gnuplot_process_input:
 * @literal: literal "extra" gnuplot commands.
 * @opt: should be OPT_I (input file) or OPT_i (input buffer)
 * @prn: gretl printing struct.
 *
 * Respond to the "gnuplot" command with %OPT_I, to specify
 * that input should be taken from a user-created gnuplot
 * command file, or %OPT_i, to take input from a named
 * string variable. If @literal is non-NULL we'll inject
 * the literal gnuplot commands immediately preceding the
 * "plot" line in the file or buffer input.
 *
 * Returns: 0 on success, error code on error.
 */

int gnuplot_process_input (const char *literal, gretlopt opt, PRN *prn)
{
    const char *iname = NULL;
    const char *buf = NULL;
    FILE *fq, *fp = NULL;
    int err = 0;

    if (opt & OPT_I) {
        iname = get_optval_string(plot_ci, OPT_I);
    } else if (opt & OPT_i) {
        iname = get_optval_string(plot_ci, OPT_i);
    }

    if (iname != NULL && *iname != '\0') {
        if (opt & OPT_I) {
            /* open the input file for reading */
            fp = gretl_fopen(iname, "r");
        } else {
            /* find the input buffer */
            buf = get_string_by_name(iname);
        }
    }

    if (fp == NULL && buf == NULL) {
        gretl_errmsg_set("Couldn't find the specified input");
        return E_INVARG;
    }

    /* open our own file for writing */
    fq = open_plot_input_file(PLOT_USER, 0, &err);

    if (!err) {
	if (literal != NULL && *literal != '\0') {
	    if (fp != NULL) {
		infile_inject_literal(fp, literal, fq);
	    } else {
		inbuf_inject_literal(buf, literal, fq);
	    }
	} else if (fp != NULL) {
            char line[1024];

            while (fgets(line, sizeof line, fp)) {
                fputs(line, fq);
            }
        } else {
	    fputs(buf, fq);
	}
	err = finalize_plot_input_file(fq);
    }

    if (fp != NULL) {
        fclose(fp);
    }

    return err;
}

/* geoplot-specific functions */

static int inline_map_data (const char *datfile, FILE *fp)
{
    char buf[8192];
    FILE *fsrc;
    size_t n;

    fsrc = gretl_fopen(datfile, "rb");
    if (fsrc == NULL) {
	gretl_errmsg_sprintf(_("Couldn't open %s"), datfile);
	return E_FOPEN;
    }

    fputs("$MapData << EOD\n", fp);
    while ((n = fread(buf, 1, sizeof buf, fsrc)) > 0) {
	fwrite(buf, 1, n, fp);
    }
    fputs("EOD\n", fp);

    fclose(fsrc);
    gretl_remove(datfile); /* not needed any more */

    return 0;
}

static gretl_matrix *geoplot_dimensions (double *xlim,
					 double *ylim,
					 int height,
					 int have_payload,
					 int *non_standard)
{
    gretl_matrix *ret = gretl_matrix_alloc(1, 2);
    double xyr = (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]);
    int width;

    if (*non_standard == 0) {
	if (fabs(ylim[0]) > 180 || fabs(ylim[1]) > 180 ||
	    fabs(xlim[0]) > 360 || fabs(xlim[1]) > 360) {
	    /* Quick and dirty check in the absence of prior
	       information that the X, Y data are not in degrees
	       of latitude and longitude.
	    */
	    fprintf(stderr, "alt coordinates: X %g to %g, Y %g to %g\n",
		    xlim[0], xlim[1], ylim[0], ylim[1]);
	    *non_standard = 1;
	}
    }

    if (*non_standard == 0) {
	/* We'll calculate a width-to-height ratio which varies
	   inversely with the (mid-point of) latitude so that we don't
	   don't get severe stretching in the x dimension when showing
	   an area far from the equator. In effect this is a cheap and
	   cheerful version of Mercator.
	*/
	double ymid = (ylim[0] + ylim[1]) / 2;

	xyr *= cos(ymid * M_PI/180);
    }

    if (have_payload) {
	/* 1.05 is to compensate for the colorbox */
	width = floor(xyr * height * 1.05);
    } else {
	width = floor(xyr * height);
    }

    set_special_plot_size(width, height);

    ret->val[0] = width;
    ret->val[1] = height;

    return ret;
}

static void fputs_literal (const char *s, FILE *fp)
{
    gchar *tmp = g_strdup(s);

    fputs(g_strchomp(tmp), fp);
    fputc('\n', fp);
    g_free(tmp);
}

static const char *map_linecolor (const char *optlc,
				  int have_payload,
				  NaAction action)
{
    if (optlc != NULL) {
	/* respect the user's choice */
	return optlc;
    } else if (have_payload) {
	return (action == NA_OUTLINE)? "gray" : "white";
    } else {
	/* outlines only */
	return "black";
    }
}

static void output_map_plot_lines (mapinfo *mi,
				   const char *datasrc,
				   const char *optlc,
				   int have_payload,
                                   int do_key,
				   double linewidth,
				   FILE *fp)
{
    const char *lc = map_linecolor(optlc, have_payload, mi->na_action);
    gchar *bline = g_strdup_printf("lc '%s' lw %g", lc, linewidth);

    if (have_payload) {
	int do_lines = linewidth > 0;
	const char *with = "with filledcurves fc palette";
	const char *cont = ", \\\n";
	int i, nv = mi->n_codes;

        if (do_key) {
            const char *kp = gretl_bundle_get_string(mi->opts, "keypos", NULL);

            /* ensure colorbox is omitted and key boxes are filled */
            fputs("unset colorbox\n", fp); /* should be handled already? */
            fputs("set style fill solid\n", fp);
            if (kp != NULL && *kp != '\0') {
                fprintf(fp, "set key %s\n", kp);
            } else {
                /* this seems to work better than the default on average? */
                fputs("set key bottom right\n", fp);
            }
        }

	/* polygons */
	fprintf(fp, "plot for [i=0:*] %s index i %s notitle%s", datasrc, with,
		(do_lines || do_key)? cont : "\n");
	if (linewidth > 0) {
	    /* plus feature outlines */
	    fprintf(fp, "  %s using 1:2 with lines %s notitle%s", datasrc, bline,
		    do_key ? cont : "\n");
	}
	if (do_key) {
	    /* plus key for discrete payload */
	    int v, v0 = mi->zvals->val[0];

	    for (i=0; i<nv; i++) {
		v = i + v0;
		if (mi->zlabels != NULL) {
                    if (strcmp(mi->zlabels[i], "empty string")) {
                        fprintf(fp, "keyentry with boxes fc palette cb %d title \"%s\"%s",
                                v, mi->zlabels[i], (i < nv - 1)? cont : "\n");
                    } else if (i == nv - 1) {
                        fputc('\n', fp);
                    }
		} else {
		    fprintf(fp, "keyentry with boxes fc palette cb %d title \"%s=%d\"%s",
			    v, mi->zname, v, (i < nv - 1)? cont : "\n");
		}
	    }
	}
    } else {
	/* just show feature outlines */
	fprintf(fp, "plot %s using 1:2 with lines %s\n", datasrc, bline);
    }

    g_free(bline);
}

static int map_do_inlining (mapinfo *mi)
{
    if (gretl_bundle_get_int(mi->opts, "inlined", NULL)) {
        return 1;
    } else if (mi->flags & MAP_MULTI) {
        return 1;
    } else {
        return 0;
    }
}

/* called from the geoplot plugin to finalize a map */

int write_map_gp_file (void *ptr)
{
    mapinfo *mi = ptr;
    gretl_bundle *opts = mi->opts;
    double gpver = gnuplot_version();
    double xlim[2], ylim[2];
    gretl_matrix *dims = NULL;
    const char *optlc = NULL;
    const char *sval;
    FILE *fp = NULL;
    gchar *datasrc = NULL;
    double linewidth = 1.0;
    double margin = 0.02;
    int display, have_payload = 0;
    int non_standard;
    int use_arg0 = 0;
    int do_key = 0;
    int height = 600;
    int border = 1;
    int notics = 1;
    int err = 0;

    if (mi->zrange != NULL) {
        have_payload = 1;
        if (mi->n_codes > 0 && gpver >= 5.4) {
	    if (mi->zlabels != NULL || mi->zname != NULL) {
		do_key = 1;
	    }
        }
    } else if (opts == NULL) {
	/* the simple outlines case */
	border = 1;
	notics = 0;
    }

    display = (mi->flags & MAP_DISPLAY) != 0;
    non_standard = (mi->flags & MAP_NON_STD) != 0;

    set_map_plot_limits(mi, xlim, ylim, margin);

    if (gretl_bundle_has_key(opts, "height")) {
	height = gretl_bundle_get_scalar(opts, "height", &err);
	if (display && height <= 0) {
	    height = 600;
	}
    }

    gretl_push_c_numeric_locale();

    if (height > 0) {
	dims = geoplot_dimensions(xlim, ylim, height, have_payload,
				  &non_standard);
    }
    if (display) {
	set_optval_string(GNUPLOT, OPT_U, "display");
	if (mi->plotfile != NULL) {
	    iact_gpfile = (char *) mi->plotfile;
	}
    } else if (mi->flags & MAP_IS_IMAGE) {
	set_optval_string(GNUPLOT, OPT_U, mi->plotfile);
    }

    fp = open_plot_input_file(PLOT_GEOMAP, 0, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "# geoplot %g %g\n", dims->val[0], dims->val[1]);

    if (!do_key) {
        fputs("unset key\n", fp);
    }

    if (have_payload) {
	err = print_map_palette(mi, gpver, fp);
        if (err) {
            fclose(fp);
            gretl_remove(gretl_plotfile());
            return err;
        }
    }

    fprintf(fp, "set xrange [%g:%g]\n", xlim[0], xlim[1]);
    fprintf(fp, "set yrange [%g:%g]\n", ylim[0], ylim[1]);

    if (gretl_bundle_has_key(opts, "title")) {
	sval = gretl_bundle_get_string(opts, "title", NULL);
	if (sval != NULL) {
	    fprintf(fp, "set title \"%s\"\n", sval);
	}
    }

    if (gretl_bundle_get_int(opts, "tics", NULL)) {
	notics = 0;
    }
    if (notics) {
        fputs("set noxtics\n", fp);
        fputs("set noytics\n", fp);
    }

    if (gretl_bundle_get_int(opts, "logscale", NULL)) {
	fputs("set logscale cb\n", fp);
    }

    if (gretl_bundle_has_key(opts, "border")) {
	/* allow override of default */
	border = gretl_bundle_get_int(opts, "border", NULL);
    }

    if (border == 0) {
	fputs("unset border\n", fp);
    }

    if ((sval = gretl_bundle_get_string(opts, "literal", NULL))) {
	fputs_literal(sval, fp);
    }

    if (gretl_bundle_has_key(opts, "linewidth")) {
	double lw = gretl_bundle_get_scalar(opts, "linewidth", &err);

	if (!err && lw >= 0) {
	    if (have_payload) {
		linewidth = lw;
	    } else if (lw >= 0.1) {
		linewidth = lw;
	    }
	}
    }
    if (gretl_bundle_has_key(opts, "linecolor")) {
	sval = gretl_bundle_get_string(opts, "linecolor", &err);
	if (!err) {
	    optlc = sval;
	}
    }

    gnuplot_missval_string(fp);

    if (map_do_inlining(mi)) {
	err = inline_map_data(mi->datfile, fp);
	if (!err) {
	    datasrc = g_strdup("$MapData");
	}
    } else if (mi->flags & MAP_IS_IMAGE) {
	/* @plotfile and @datfile are both disposable, no need
	   to bother about name alignment */
	datasrc = g_strdup_printf("\"%s\"", mi->datfile);
    } else if (mi->plotfile != NULL) {
	/* the names of @plotfile and @datfile will already be
	   correctly aligned */
	use_arg0 = 1;
    } else {
	/* rename @datfile to match the auto-named plot file */
	gchar *tmp = g_strdup_printf("%s.dat", gretl_plotfile());

	gretl_copy_file(mi->datfile, tmp);
	gretl_remove(mi->datfile);
	g_free(tmp);
	use_arg0 = 1;
    }

    if (use_arg0) {
	fputs("datafile = sprintf(\"%s.dat\", ARG0)\n", fp);
	datasrc = g_strdup("datafile");
    }

    if (!err) {
	output_map_plot_lines(mi, datasrc, optlc, have_payload,
			      do_key, linewidth, fp);
    }

    err = finalize_plot_input_file(fp);
    if (!err) {
	if (display && gretl_in_gui_mode()) {
	    if (gretl_bundle_get_int(opts, "gui_auto", NULL)) {
		gretl_bundle_set_string(opts, "plotfile", gretl_plotfile());
		gretl_bundle_set_matrix(opts, "dims", dims);
	    } else {
		manufacture_gui_callback(GNUPLOT);
	    }
	}
    }

    gretl_pop_c_numeric_locale();

    gretl_matrix_free(dims);
    g_free(datasrc);

    return err;
}

/* Transcribe geoplot map file from @src to @dest, allowing for
   the possibility that data contained in a separate datafile
   have to be inlined. If @datname is NULL it is assumed that
   the datafile will be named as @src, with ".dat" appended.
*/

int transcribe_geoplot_file (const char *src,
			     const char *dest,
			     const char *datname)
{
    FILE *f1 = NULL, *f2 = NULL;
    const char *mapdata = "$MapData";
    char buf[8196];
    int integrate = -1;
    int n, err = 0;

    f1 = gretl_fopen(src, "rb");
    f2 = gretl_fopen(dest, "wb");

    if (f1 == NULL || f2 == NULL) {
	err = E_FOPEN;
	goto bailout;
    }

    while (integrate < 0 && fgets(buf, sizeof buf, f1)) {
	if (strstr(buf, mapdata)) {
	    integrate = 0;
	} else if (!strncmp(buf, "datafile =", 10)) {
	    integrate = 1;
	} else {
	    fputs(buf, f2);
	}
    }

    if (integrate == 1) {
	/* open the datafile and inject its content */
	gchar *s, *dattmp = NULL;
	FILE *fdat = NULL;
	int i;

	if (datname != NULL) {
	    fdat = gretl_fopen(datname, "rb");
	} else {
	    dattmp = g_strdup_printf("%s.dat", src);
	    fdat = gretl_fopen(dattmp, "rb");
	    g_free(dattmp);
	}

	if (fdat == NULL) {
	    err = E_FOPEN;
	} else {
	    /* inject data */
	    fprintf(f2, "%s << EOD\n", mapdata);
	    while ((n = fread(buf, 1, sizeof buf, fdat)) > 0) {
		fwrite(buf, 1, n, f2);
	    }
	    fputs("EOD\n", f2);
	    fclose(fdat);
	    buf[0] = '\0';

	    /* and transcribe the remainder of @src */
	    while (fgets(buf, sizeof buf, f1)) {
		if ((s = strstr(buf, "datafile")) != NULL) {
		    for (i=0; i<8; i++) {
			s[i] = mapdata[i];
		    }
		}
		fputs(buf, f2);
	    }
	}
    } else if (integrate == 0) {
	/* integration not needed */
	while ((n = fread(buf, 1, sizeof buf, f1)) > 0) {
	    fwrite(buf, 1, n, f2);
	}
    } else {
	/* ?? */
	err = E_DATA;
    }

 bailout:

    if (f1 != NULL) fclose(f1);
    if (f2 != NULL) fclose(f2);

    return err;
}

/* called from the interpolate plugin */

int write_tdisagg_plot (const gretl_matrix *YY, int mult,
			const char *title, DATASET *dset)
{
    const double *obs = NULL;
    char mstr[16] = {0};
    int t, T = YY->rows;
    double y0t;
    FILE *fp;
    int err = 0;

    set_optval_string(GNUPLOT, OPT_U, "display");
    fp = open_plot_input_file(PLOT_REGULAR, GPT_LETTERBOX, &err);
    if (err) {
	return err;
    }

    if (dset != NULL) {
	fprintf(fp, "# timeseries %d (letterbox)\n", dset->pd);
	obs = gretl_plotx(dset, OPT_NONE);
    } else {
	fputs("# timeseries 1 (letterbox)\n", fp);
    }
    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (title != NULL) {
	fprintf(fp, "set title \"%s\"\n", title);
    }

    gretl_push_c_numeric_locale();

    if (obs != NULL) {
	double d1 = obs[dset->t1];
	double d2 = obs[dset->t2];

	fprintf(fp, "set xrange [%g:%g]\n", floor(d1), ceil(d2));
    }

    gnuplot_missval_string(fp);
    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (t=0; t<T; t++) {
	if (obs != NULL) {
	    fprintf(fp, "%g ", obs[t+dset->t1]);
	} else {
	    fprintf(fp, "%d ", t + 1);
	}
	y0t = gretl_matrix_get(YY, t, 0);
	if (na(y0t)) {
	    fputs("? ", fp);
	} else {
	    fprintf(fp, "%.10g ", y0t);
	}
	fprintf(fp, "%.10g\n", gretl_matrix_get(YY, t, 1));
    }
    fputs("EOD\n", fp);

    if (mult > 1) {
	sprintf(mstr, " * %d", mult);
    }

    fprintf(fp, "plot $data using 1:2 title \"%s\" w steps, \\\n",
	    _("original data"));
    fprintf(fp, " $data using 1:3 title \"%s%s\" w lines\n",
	    _("final series"), mstr);

    err = finalize_plot_input_file(fp);

    if (!err && gretl_in_gui_mode()) {
	manufacture_gui_callback(GNUPLOT);
    }

    gretl_pop_c_numeric_locale();

    return err;
}
