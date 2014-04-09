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

#include <unistd.h>
#include <glib.h>

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

/* experimental, AC 2013-10-25 */
#define USE_TIMEFMT 1

static char gnuplot_path[MAXLEN];
static int gp_small_font_size;
static double default_png_scale = 1.0;

typedef struct gnuplot_info_ gnuplot_info;

struct gnuplot_info_ {
    GptFlags flags;
    FitType fit;
    int *list;
    int t1;
    int t2;
    double xrange;
    char timefmt[16];
    char xtics[64];
    char xfmt[16];
    char yfmt[16];
    const char *yformula;
    const double *x;
    gretl_matrix *dvals;
    int *withlist;
};

enum {
    W_POINTS,
    W_LINES,
    W_IMPULSES,
    W_LP
};

#define MAX_LETTERBOX_LINES 8

#define ts_plot(g)      ((g)->flags & GPT_TS)

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
    { PLOT_GARCH,          "GARCH residual plot" },
    { PLOT_HURST,          "rescaled range plot" },
    { PLOT_IRFBOOT,        "impulse response plot with quantiles" },
    { PLOT_KERNEL,         "kernel density plot" },
    { PLOT_LEVERAGE,       "leverage/influence plot" },
    { PLOT_MULTI_SCATTER,  "multiple scatterplots" },
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
    { PLOT_TYPE_MAX,       NULL }
};

static int graph_list_adjust_sample (int *list, 
				     gnuplot_info *ginfo,
				     const DATASET *dset,
				     int listmin);
static void clear_gpinfo (gnuplot_info *gi);
static void make_time_tics (gnuplot_info *gi,
			    const DATASET *dset,
			    int many, char *xlabel,
			    PRN *prn);
static void get_multiplot_layout (int n, int tseries,
				  int *rows, int *cols);
    
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

    // signal(SIGCHLD, SIG_DFL);

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

double gnuplot_version (void)
{
    static double vnum = 0.0;

    if (vnum == 0.0) {
	gboolean ok;
	gchar *sout = NULL;
	gchar *argv[] = {
	    NULL,
	    NULL,
	    NULL
	};

	if (*gnuplot_path == '\0') {
	    strcpy(gnuplot_path, gretl_gnuplot_path());
	}

	argv[0] = gnuplot_path;
	argv[1] = "--version";

	ok = g_spawn_sync (NULL,
			   argv,
			   NULL,
			   G_SPAWN_SEARCH_PATH |
			   G_SPAWN_STDERR_TO_DEV_NULL,
			   NULL,
			   NULL,
			   &sout,
			   NULL,
			   NULL,
			   NULL);

	if (ok && sout != NULL) {
	    if (!strncmp(sout, "gnuplot ", 8)) {
		/* e.g. "gnuplot 4.7 patchlevel 0" */
		vnum = dot_atof(sout + 8);
	    }
	    g_free(sout);
	}
    }

    return vnum;
}

#else /* MS Windows */

# ifdef WIN64
double gnuplot_version (void)
{
    /* As of the gretl 1.9.13 release, the package for 
       64-bit Windows includes gnuplot 4.7 */
    return 4.7;
}
# else
double gnuplot_version (void)
{
    /* As of the gretl 1.9.13 release, the package for 
       32-bit Windows includes gnuplot 4.6.3 */
    return 4.63;
}
# endif

#endif /* MS Windows or not */

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
    const char *s = get_optval_string(GNUPLOT, opt);
    int withval = W_POINTS;
    int i, imax = gi->withlist[0];

    if (opt == OPT_O) {
	withval = W_LINES;
    } else if (opt == OPT_M) {
	withval = W_IMPULSES;
    } else if (opt == OPT_P) {
	withval = W_LP;
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
    if ((opt & OPT_O) && !(opt & (OPT_M | OPT_P))) {
	return get_optval_string(GNUPLOT, OPT_O) == NULL;
    } else {
	return 0;
    }
}

static void get_gp_flags (gnuplot_info *gi, gretlopt opt, 
			  const int *list, const DATASET *dset)
{
    int n_yvars = list[0] - 1;

    gi->flags = 0;

    if (opt & OPT_R) {
	gi->flags |= GPT_RESIDS;
    } else if (opt & OPT_F) {
	gi->flags |= GPT_FA;
    }

    if (opt & OPT_Z) {
	gi->flags |= GPT_DUMMY;
    } else if (opt & OPT_C) {
	gi->flags |= GPT_XYZ;
    } else {
	if (opt & OPT_S) {
	    gi->flags |= GPT_FIT_OMIT;
	}
	if (opt & OPT_T) {
	    gi->flags |= GPT_IDX;
	    /* there's no xvar in @list */
	    n_yvars++;
	}
    }

    if (plain_lines_spec(opt)) {
	/* just using lines */
	gi->flags |= GPT_LINES;
    } else if (opt & (OPT_M | OPT_O | OPT_P)) {
	/* for handling per-variable "plot with" options */
	gi->withlist = gretl_list_new(n_yvars);
    }

    if (gi->withlist != NULL) {
	if (opt & OPT_M) {
	    gp_set_non_point_info(gi, list, dset, OPT_M);
	}
	if (opt & OPT_O) {
	    gp_set_non_point_info(gi, list, dset, OPT_O);
	}
	if (opt & OPT_P) {
	    gp_set_non_point_info(gi, list, dset, OPT_P);
	}
    }

    gi->fit = PLOT_FIT_NONE;

    if (!(gi->flags & GPT_FIT_OMIT) && n_yvars == 1) {
	if (opt & OPT_I) {
	    gi->fit = PLOT_FIT_INVERSE;
	} else if (opt & OPT_Q) {
	    gi->fit = PLOT_FIT_QUADRATIC;
	} else if (opt & OPT_B) {
	    gi->fit = PLOT_FIT_CUBIC;
	} else if (opt & OPT_L) {
	    gi->fit = PLOT_FIT_LOESS;
	} else if (opt & OPT_N) {
	    gi->fit = PLOT_FIT_OLS;
	} else if (opt & OPT_E) {
	    gi->fit = PLOT_FIT_LOGLIN;
	}
    }

#if GP_DEBUG
    print_gnuplot_flags(gi->flags, 0);
#endif
}

static void printvars (FILE *fp, int t, 
		       const int *list, 
		       const double **Z,
		       const double *x, 
		       const char *label, 
		       const char *date,
		       double offset)
{
    double xt;
    int i;

    if (date != NULL) {
	fprintf(fp, "%s ", date);
    } else if (x != NULL) {
	xt = x[t] + offset;
	fprintf(fp, "%.10g ", xt);
    }

    for (i=1; i<=list[0]; i++) {
	xt = Z[list[i]][t];
	if (na(xt)) {
	    fputs("? ", fp);
	} else {
	    if (x == NULL && i == 1) { 
		/* the x variable */
		xt += offset;
	    }
	    fprintf(fp, "%.10g ", xt);
	}
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);
}

static int factor_check (gnuplot_info *gi, const DATASET *dset)
{
    int err = 0;

    if (gi->list[0] != 3) {
	err = E_DATA;
    } else {
	int v3 = gi->list[3];

	if (!series_is_discrete(dset, v3) &&
	    !gretl_isdiscrete(gi->t1, gi->t2, dset->Z[v3])) {
	    err = E_DATA;
	}
    }

    if (err) {
	gretl_errmsg_set(_("You must supply three variables, the last of "
			   "which is discrete"));
    } else {
	const double *d = dset->Z[gi->list[3]] + gi->t1;
	int T = gi->t2 - gi->t1 + 1;

	gi->dvals = gretl_matrix_values(d, T, OPT_S, &err);
    }

    return err;
}

#ifdef WIN32

int gnuplot_has_ttf (int reset)
{
    /* we know the gnuplot supplied with gretl for win32
       does TrueType fonts */
    return 1;
}

int gnuplot_pdf_terminal (void)
{
    return GP_PDF_CAIRO;
}

int gnuplot_eps_terminal (void)
{
    return GP_EPS_CAIRO;
}

int gnuplot_png_terminal (void)
{
    return GP_PNG_CAIRO;
}
   
int gnuplot_has_cp950 (void)
{
    /* ... and that it supports CP950 */
    return 1;
}

#else /* !WIN32 */

int gnuplot_has_ttf (int reset)
{
    static int err = -1; 

    if (err == -1 || reset) {
	/* if we have cairo we know we (should be!) OK */
        err = gnuplot_test_command("set term pngcairo");
	if (err) {
	    /* otherwise (libgd) try some plausible ttf fonts */
	    err = gnuplot_test_command("set term png font Vera 8");
	    if (err) {
		err = gnuplot_test_command("set term png font luxisr 8");
	    }
	    if (err) {
		err = gnuplot_test_command("set term png font arial 8");
	    }
	}
    }

    return !err;
}

int gnuplot_has_cp950 (void)
{
    static int err = -1; 

    /* not OK in gnuplot 4.4.0 */

    if (err == -1) {
	err = gnuplot_test_command("set encoding cp950");
    }

    return !err;
}

int gnuplot_pdf_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
	int err = gnuplot_test_command("set term pdfcairo");

	if (!err) {
	    ret = GP_PDF_CAIRO;
	} else {
	    err = gnuplot_test_command("set term pdf");
	    if (!err) {
		ret = GP_PDF_PDFLIB;
	    } else {
		ret = GP_PDF_NONE;
	    }
	}
    }

    return ret;
}

int gnuplot_eps_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
	int err = gnuplot_test_command("set term epscairo");

	/* not OK in gnuplot 4.4.0 */

	if (!err) {
	    ret = GP_EPS_CAIRO;
	} else {
	    ret = GP_EPS_PS;
	}
    }

    return ret;
}

static int gnuplot_has_x11 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term x11");
    }

    return !err;
}

static int gnuplot_has_aqua (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term aqua");
    }

    return !err;
}

int gnuplot_has_wxt (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term wxt");
    }

    return !err;
}

int gnuplot_png_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
	int err = gnuplot_test_command("set term pngcairo");

	if (!err) {
	    fprintf(stderr, "gnuplot: using pngcairo driver\n");
	    ret = GP_PNG_CAIRO;
	} else {
	    fprintf(stderr, "gnuplot: using libgd png driver\n");
	    err = gnuplot_test_command("set term png truecolor");
	    ret = (err)? GP_PNG_GD1 : GP_PNG_GD2;
	}
    }

    return ret;
}

#endif /* !WIN32 */

static int gnuplot_png_use_aa = 1;

void gnuplot_png_set_use_aa (int s)
{
    gnuplot_png_use_aa = s;
}

/* apparatus for handling plot colors */

static const gretlRGB default_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 }, /* full-intensity green is not very legible */
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 },  /* box fill */
    { 0xdd, 0xdd, 0xdd },  /* shade fill */    
};

static gretlRGB user_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 },
    { 0xbf, 0x25, 0xb2 },
    { 0x8f, 0xaa, 0xb3 },
    { 0xff, 0xa5, 0x00 },
    { 0x5f, 0x6b, 0x84 },
    { 0xdd, 0xdd, 0xdd }    
};

void print_rgb_hash (char *s, const gretlRGB *color)
{
    sprintf(s, "#%02x%02x%02x", color->r, color->g, color->b);
}

void gretl_rgb_get (gretlRGB *color, const char *s)
{
    int n, r, g, b;

    n = sscanf(s, "#%2x%2x%2x", &r, &g, &b);

    if (n == 3) {
	color->r = r;
	color->g = g;
	color->b = b;
    } else {
	color->r = color->g = color->b = 0;
    }
}

void print_palette_string (char *s)
{
    char colstr[8];
    int i;

    *s = '\0';

    for (i=0; i<N_GP_COLORS; i++) {
	sprintf(colstr, "x%02x%02x%02x", user_color[i].r, user_color[i].g, 
		user_color[i].b);
	strcat(s, colstr);
	if (i < N_GP_COLORS - 1) {
	    strcat(s, " ");
	}
    }
}

const gretlRGB *get_graph_color (int i)
{
    return (i >= 0 && i < N_GP_COLORS)? &user_color[i] : NULL;
}

void set_graph_palette (int i, gretlRGB color)
{
    if (i >= 0 && i < N_GP_COLORS) {
	user_color[i] = color;
    } else {
	fprintf(stderr, "Out of bounds color index %d\n", i);
    }
}

void set_graph_palette_from_string (int i, const char *s)
{
    int err = 0;

    if (i >= 0 && i < N_GP_COLORS) {
	unsigned int x[3];

	if (sscanf(s + 1, "%02x%02x%02x", &x[0], &x[1], &x[2]) == 3) {
	    user_color[i].r = x[0];
	    user_color[i].g = x[1];
	    user_color[i].b = x[2];
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
    if (i == SHADECOLOR) {
	user_color[SHADECOLOR] = default_color[SHADECOLOR];
    } else if (i == BOXCOLOR) {
	user_color[BOXCOLOR] = default_color[BOXCOLOR];
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    user_color[i] = default_color[i];
	}
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

static void maybe_set_small_font (int nplots)
{
    gp_small_font_size = (nplots > 4)? 6 : 0;
}

static void 
write_gnuplot_font_string (char *fstr, PlotType ptype, int pngterm,
			   const char *grfont, double scale)
{
    if (grfont == NULL) {
	grfont = gretl_png_font();
    }

    if (*grfont == '\0') {
	grfont = getenv("GRETL_PNG_GRAPH_FONT");
    }

    if (grfont == NULL || *grfont == '\0') {
	*fstr = '\0';
	return;
    }

    if (pngterm == GP_PNG_CAIRO) {
	char fname[128];
	int nf, fsize = 0;

	nf = split_graph_fontspec(grfont, fname, &fsize);
	if (nf == 2) {
	    if (maybe_big_multiplot(ptype) && gp_small_font_size > 0) {
		fsize = gp_small_font_size;
	    }
	    if (scale > 1.0) {
		fsize = round(scale * fsize);
	    }
	    sprintf(fstr, " font \"%s,%d\"", fname, fsize);
	} else if (nf == 1) {
	    sprintf(fstr, " font \"%s\"", fname);
	}
    } else {
	int shrink = 0;

	if (maybe_big_multiplot(ptype) && gp_small_font_size > 0) {
	    char fname[64];
	    int fsize;

	    if (sscanf(grfont, "%s %d", fname, &fsize) == 2) {
		sprintf(fstr, " font %s %d", fname, gp_small_font_size);
		shrink = 1;
	    }
	}

	if (!shrink) {
	    sprintf(fstr, " font %s", grfont);
	}
    }
}

#ifndef WIN32

static void 
write_old_gnuplot_font_string (char *fstr, PlotType ptype)
{
    if (maybe_big_multiplot(ptype)) {
	strcpy(fstr, " tiny");
    } else {
	strcpy(fstr, " small");
    }
}

#endif

/* requires gnuplot 4.2 or higher */

void write_plot_line_styles (int ptype, FILE *fp)
{
    char cstr[8];
    int i;
    
    if (frequency_plot_code(ptype)) {
	print_rgb_hash(cstr, &user_color[BOXCOLOR]);
	fprintf(fp, "set style line 1 lc rgb \"%s\"\n", cstr);
	fputs("set style line 2 lc rgb \"#000000\"\n", fp);
    } else if (ptype == PLOT_RQ_TAU) {
	fputs("set style line 1 lc rgb \"#000000\"\n", fp);
	for (i=1; i<BOXCOLOR; i++) {
	    print_rgb_hash(cstr, &user_color[i]);
	    fprintf(fp, "set style line %d lc rgb \"%s\"\n", i+1, cstr);
	}
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    print_rgb_hash(cstr, &user_color[i]);
	    fprintf(fp, "set style line %d lc rgb \"%s\"\n", i+1, cstr);
	}
	print_rgb_hash(cstr, &user_color[SHADECOLOR]);
	fprintf(fp, "set style line %d lc rgb \"%s\"\n", 
		SHADECOLOR + 1, cstr);
    }

    fputs("set style increment user\n", fp);
}

#ifdef WIN32

static void reslash_filename (char *buf, const char *src)
{
    strcpy(buf, src);

    while (*buf) {
	if (*buf == '\\') *buf = '/';
	buf++;
    }
}

#endif

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds (gnuplot >= 4.4.0).
*/

void write_plot_bounding_box_request (FILE *fp)
{
#ifdef WIN32
    char buf[FILENAME_MAX];

    reslash_filename(buf, gretl_dotdir());
    fprintf(fp, "set print \"%sgretltmp.png.bounds\"\n", buf);
#else
    fprintf(fp, "set print \"%sgretltmp.png.bounds\"\n", gretl_dotdir());
#endif
    fputs("print \"pixel_bounds: \", GPVAL_TERM_XMIN, GPVAL_TERM_XMAX, "
	  "GPVAL_TERM_YMIN, GPVAL_TERM_YMAX\n", fp);
    fputs("print \"data_bounds: \", GPVAL_X_MIN, GPVAL_X_MAX, "
	  "GPVAL_Y_MIN, GPVAL_Y_MAX\n", fp);
}

static void do_plot_bounding_box (void)
{
    FILE *fp = fopen(gretl_plotfile(), "a");

    if (fp != NULL) {
	write_plot_bounding_box_request(fp);
	fclose(fp);
    }
}

static void maybe_set_eps_pdf_dims (char *s, PlotType ptype, GptFlags flags)
{
    double w = 0, h = 0;

    if (flags & GPT_LETTERBOX) {
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

const char *get_gretl_pdf_term_line (PlotType ptype, GptFlags flags)
{
    static char pdf_term_line[128];

    if (gnuplot_pdf_terminal() == GP_PDF_CAIRO) {
	int ptsize = 10;

	if (ptype == PLOT_MULTI_SCATTER) {
	    ptsize = 6;
	}
#ifndef WIN32
	if (gnuplot_version() <= 4.4) {
	    ptsize /= 2;
	}
#endif
	sprintf(pdf_term_line, "set term pdfcairo font \"sans,%d\"", 
		ptsize);
    } else {
	strcpy(pdf_term_line, "set term pdf");
    }

    maybe_set_eps_pdf_dims(pdf_term_line, ptype, flags);

    return pdf_term_line;
}

const char *get_gretl_eps_term_line (PlotType ptype, GptFlags flags)
{
    static char eps_term_line[128];

    if (flags & GPT_MONO) {
	strcpy(eps_term_line, "set term post eps enhanced mono");
    } else {
	strcpy(eps_term_line, "set term post eps enhanced color solid");
    }

    maybe_set_eps_pdf_dims(eps_term_line, ptype, flags);

    return eps_term_line;
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

    if (flags & GPT_LETTERBOX) {
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
    } else if (ptype == PLOT_ROOTS || ptype == PLOT_QQ) {
	/* square plots */
	w = h = GP_SQ_SIZE;
    }

    if (scale != 1.0) {
	plot_get_scaled_dimensions(&w, &h, scale);
    }

    *s = '\0';

    sprintf(s, " size %d,%d", w, h);
}

static const char *real_png_term_line (PlotType ptype, 
				       GptFlags flags,
				       const char *specfont,
				       double scale)
{
    static char png_term_line[256];
    char truecolor_string[12] = {0};
    char font_string[128];
    char size_string[16];
    int gpttf, pngterm = 0;

    *font_string = *size_string = '\0';

    pngterm = gnuplot_png_terminal();

#ifdef WIN32
    gpttf = 1;
#else
    gpttf = gnuplot_has_ttf(0);
#endif

    if (pngterm == GP_PNG_GD2 && gnuplot_png_use_aa) {
	strcpy(truecolor_string, " truecolor");
    }   

    if (gpttf) {
	write_gnuplot_font_string(font_string, ptype, pngterm,
				  specfont, scale);
    } 

#ifndef WIN32
    if (!gpttf) {
	write_old_gnuplot_font_string(font_string, ptype);
    }
#endif

    write_png_size_string(size_string, ptype, flags, scale);

    if (pngterm == GP_PNG_CAIRO) {
	sprintf(png_term_line, "set term pngcairo%s%s",
		font_string, size_string);
	strcat(png_term_line, "\nset encoding utf8");
    } else {
	sprintf(png_term_line, "set term png%s%s%s",
		truecolor_string, font_string, size_string); 
    }

#if GP_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", png_term_line);
#endif

    return png_term_line;
}

/**
 * get_gretl_png_term_line:
 * @ptype: indication of the sort of plot to be made, which
 * may made a difference to the color palette chosen.
 * @flags: plot option flags.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the PNG "terminal".  Checks the environment for setting of 
 * %GRETL_PNG_GRAPH_FONT.  Also appends a color-specification string 
 * if the gnuplot PNG driver supports this.
 *
 * Returns: the terminal string, "set term png ..."
 */

const char *get_gretl_png_term_line (PlotType ptype, GptFlags flags)
{
    double s = default_png_scale;

    return real_png_term_line(ptype, flags, NULL, s);
}

const char *get_png_line_for_plotspec (const GPT_SPEC *spec)
{
    return real_png_term_line(spec->code, spec->flags, 
			      spec->fontstr, spec->scale);
}

void gnuplot_png_set_default_scale (double s)
{
    if (s >= 0.5 && s <= 2.0) {
	default_png_scale = s;
    }
}

static void png_font_to_emf (const char *pngfont, char *emfline)
{
    char name[128];
    int pt;

    if (split_graph_fontspec(pngfont, name, &pt) == 2) {
	char ptstr[8];

	if (pt <= 8) {
	    pt = 12;
	} else {
	    pt = 16;
	}

	strcat(emfline, "'");
	strcat(emfline, name);
	strcat(emfline, "' ");
	sprintf(ptstr, "%d ", pt);
	strcat(emfline, ptstr);
    }
}

/**
 * get_gretl_emf_term_line:
 * @ptype: indication of the sort of plot to be made.
 * @color: 1 if graph is to be in color, else 0.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the EMF "terminal".  
 *
 * Returns: the term string, "set term emf ..."
 */

const char *get_gretl_emf_term_line (PlotType ptype, int color)
{
    static char tline[256];
    const char *grfont = NULL;
    
    strcpy(tline, "set term emf ");

    if (color) {
	strcat(tline, "color ");
    } else {
	strcat(tline, "mono dash ");
    }

    /* font spec */
    grfont = gretl_png_font();
    if (grfont != NULL && *grfont != '\0') {
	png_font_to_emf(grfont, tline);
    }

    return tline;
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

    for (i=1; i<PLOT_TYPE_MAX; i++) {
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

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (ptype == ptinfo[i].ptype) {
	    if (flags & GPT_XL) {
		fprintf(fp, "# %s (large)\n", ptinfo[i].pstr);
	    } else if (flags & GPT_XXL) {
		fprintf(fp, "# %s (extra-large)\n", ptinfo[i].pstr);
	    } else {
		fprintf(fp, "# %s\n", ptinfo[i].pstr);
	    }
	    ret = 1;
	    break;
	}
    }

    if (get_local_decpoint() == ',') {
	fputs("set decimalsign ','\n", fp);
    }

    return ret;
}

static void print_term_string (int ttype, PlotType ptype,
			       GptFlags flags, FILE *fp)
{
    const char *tstr = NULL;

    if (ttype == GP_TERM_EPS) {
	tstr = get_gretl_eps_term_line(ptype, flags);
    } else if (ttype == GP_TERM_PDF) {
	tstr = get_gretl_pdf_term_line(ptype, flags);
    } else if (ttype == GP_TERM_PNG) {
	tstr = get_gretl_png_term_line(ptype, 0);
    } else if (ttype == GP_TERM_EMF) {
	tstr = get_gretl_emf_term_line(ptype, 1);
    } else if (ttype == GP_TERM_FIG) {
	tstr = "set term fig";
    } else if (ttype == GP_TERM_SVG) {
	tstr = "set term svg";
    }

    if (tstr != NULL) {
	fprintf(fp, "%s\n", tstr);
	if (ttype != GP_TERM_EPS) {
	    write_plot_line_styles(PLOT_REGULAR, fp);
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

/* if @path is non-NULL we use it, otherwise we make a path
   using dotdir and gretltmp.png
*/

void write_plot_output_line (const char *path, FILE *fp)
{
#ifdef WIN32
    char buf[FILENAME_MAX];

    if (path == NULL) {
	reslash_filename(buf, gretl_dotdir());
	fprintf(fp, "set output \"%sgretltmp.png\"\n", buf);
    } else {
	reslash_filename(buf, path);
	fprintf(fp, "set output \"%s\"\n", buf);
    }
#else
    if (path == NULL) {
	fprintf(fp, "set output \"%sgretltmp.png\"\n", gretl_dotdir());
    } else {
	fprintf(fp, "set output \"%s\"\n", path);
    }
#endif
}

static FILE *gp_set_up_batch (char *fname, 
			      PlotType ptype, 
			      GptFlags flags,
			      const char *optname,
			      int *err)
{
    int fmt = GP_TERM_NONE;
    FILE *fp = NULL;

    if (optname != NULL) {
	/* user gave --output=<filename> */
	fmt = set_term_type_from_fname(optname);
	if (fmt) {
	    /* input needs processing */
	    strcpy(gnuplot_outname, optname);
	    gretl_maybe_prepend_dir(gnuplot_outname);
	    sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
	    fp = gretl_mktemp(fname, "w");
	} else {
	    /* just passing gnuplot commands through */
	    this_term_type = GP_TERM_PLT;
	    strcpy(fname, optname);
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
	set_gretl_plotfile(fname);
	if (*gnuplot_outname != '\0') {
	    /* write terminal/style/output lines */
	    print_term_string(fmt, ptype, flags, fp);
	    write_plot_output_line(gnuplot_outname, fp);
	} else {
	    /* just write style lines */
	    write_plot_line_styles(PLOT_REGULAR, fp);
	}
	if (get_local_decpoint() == ',') {
	    fputs("set decimalsign ','\n", fp);
	}
    }

    return fp;
}

/* Set-up for an "interactive" plot: we open a file in the user's
   dotdir into which gnuplot commands will be written.  If we're
   running the GUI program this command file will eventually be used
   to create a PNG file for display in a gretl window; otherwise
   (gretlcli) the commands will eventually be sent to gnuplot for
   "direct" display (e.g. using the x11 or windows terminal).

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

    if (gui) {
	/* the filename should be unique */
	sprintf(fname, "%sgpttmp.XXXXXX", gretl_dotdir());
	fp = gretl_mktemp(fname, "w");
    } else {
	/* gretlcli: no need for uniqueness */
	sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
	fp = gretl_fopen(fname, "w");
    }

    if (fp == NULL) {
	*err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
	if (gui) {
	    /* set up for PNG output */
	    fprintf(fp, "%s\n", get_gretl_png_term_line(ptype, flags));
	    if (default_png_scale != 1.0) {
		gretl_push_c_numeric_locale();
		fprintf(fp, "# scale = %.1f\n", default_png_scale);
		fputs("# auto linewidth\n", fp);
		gretl_pop_c_numeric_locale();
	    }
	    write_plot_output_line(NULL, fp);
	}
	write_plot_type_string(ptype, flags, fp);
	write_plot_line_styles(ptype, fp);
    }

    return fp;
}

#ifndef WIN32

static int gnuplot_too_old (void)
{
    static double gpv;

    if (gpv == 0.0) {
	gpv = gnuplot_version();
    }

    if (gpv < 4.4) {
	gretl_errmsg_set("Gnuplot is too old: must be >= version 4.4.0");
	return 1;
    } else {
	return 0;
    }
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

static const char *plot_output_option (PlotType p, int *pci)
{
    int ci = GNUPLOT;
    const char *s;

    /* set a more specific command index based on 
       the plot type, if applicable */

    if (p == PLOT_MULTI_SCATTER) {
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
	       p == PLOT_FREQ_GAMMA) {
	ci = FREQ;
    }

    s = get_optval_string(ci, OPT_U);
    if (s != NULL && *s == '\0') {
	s = NULL;
    }

    if (pci != NULL) {
	*pci = ci;
    }

    return s;
}

/* Open a file into which gnuplot commands will be written.

   Depending on the prospective use of the stream, we
   may write some header-type stuff into it, the primary
   case being when we're going to produce PNG output
   for display in the GUI.

   This internal function allows specification of non-zero
   GptFlags -- cf. open_plot_input_file(), which is a simpler
   wrapper for open_gp_stream().
*/

static FILE *open_gp_stream (PlotType ptype, GptFlags flags, 
			     int *err)
{
    char fname[FILENAME_MAX] = {0};
    const char *optname = NULL;
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

    /* check for --output=whatever option */
    optname = plot_output_option(ptype, &ci);

    if (ci == FREQ || got_display_option(optname)) {
	/* --output=display specified (or FREQ, a special) */
	interactive = 1;
    } else if (optname != NULL) {
	/* --output=filename specified */
	interactive = 0;
    } else {
	/* defaults */
	interactive = !gretl_in_batch_mode();
    }

#if GP_DEBUG
    fprintf(stderr, "optname = '%s', interactive = %d\n", 
	    optname, interactive);
#endif

    if (interactive) {
	fp = gp_set_up_interactive(fname, ptype, flags, err);
    } else {
	fp = gp_set_up_batch(fname, ptype, flags, optname, err);
    }

#if GP_DEBUG
    fprintf(stderr, "open_gp_stream_full: '%s'\n", gretl_plotfile());
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

    if (opt & OPT_U) {
	/* check for --plot=whatever option */
	optname = plot_output_option(ptype, NULL);
    }

    if (got_none_option(optname)) {
	/* --plot=none specified */
	ret = 0;
    } else if (optname != NULL) {
	/* --plot=display or --plot=fname specified */
	ret = 1;
    } else {
	/* defaults */
	ret = !gretl_in_batch_mode();
    }

    return ret;
}

/**
 * open_plot_input_file:
 * @ptype: indication of the sort of plot to be made.
 * @err: location to receive error code.
 *
 * Opens a file into which gnuplot commands will be written.
 * Depending on the prospective use of the stream, some
 * header-type material may be written into it (the primary
 * case being when we're going to produce PNG output
 * for display in the gretl GUI). The prospective use is
 * figured out based on the program state and @ptype.
 *
 * Returns: writable stream on success, %NULL on failure.
 */

FILE *open_plot_input_file (PlotType ptype, int *err)
{
    FILE *fp = open_gp_stream(ptype, 0, err);

    return fp;
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
}

static int graph_file_written;

int graph_written_to_file (void)
{
    return graph_file_written;
}

static void remove_old_png (char *buf)
{
    sprintf(buf, "%sgretltmp.png", gretl_dotdir());
    remove(buf);
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
    fmt = specified_gp_output_format();

    if (fmt == GP_TERM_PLT) {
	/* no-op: just the gnuplot commands are wanted */
	graph_file_written = 1;
	return 0;
    } else if (fmt == GP_TERM_PDF) {
	/* can we do this? */
	if (gnuplot_pdf_terminal() == GP_PDF_NONE) {
	    gretl_errmsg_set(_("Gnuplot does not support PDF output "
			       "on this system"));
	    return E_EXTERNAL;
	}
    } else if (fmt == GP_TERM_NONE && gui) {
	do_plot_bounding_box();
	/* ensure we don't get stale output */
	remove_old_png(buf);
    }

#ifdef WIN32
    sprintf(buf, "\"%s\" \"%s\"", gretl_gnuplot_path(), fname);
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
	    remove(fname);
	    set_gretl_plotfile(gnuplot_outname);
	    graph_file_written = 1;
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
	       indicate that we actually produced a plot */
	    set_plot_produced();
	}
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
    char title[128];

    switch (code) {
    case GTITLE_VLS:
	if (gi->fit == PLOT_FIT_OLS) {
	    sprintf(title, _("%s versus %s (with least squares fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    sprintf(title, _("%s versus %s (with inverse fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    sprintf(title, _("%s versus %s (with quadratic fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_CUBIC) {	 
	    sprintf(title, _("%s versus %s (with cubic fit)"),
		    s1, s2);
	}	    
	break;
    case GTITLE_RESID:
	if (strncmp(s1, "residual for ", 13) == 0 &&
	    gretl_scan_varname(s1 + 13, depvar) == 1) {
	    sprintf(title, _("Regression residuals (= observed - fitted %s)"), 
		    depvar);
	}
	break;
    case GTITLE_AF:
	sprintf(title, _("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || (gi->flags & GPT_TS)) {
	    sprintf(title, _("Actual and fitted %s"), s1);
	} else {
	    sprintf(title, _("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	*title = '\0';
	break;
    }

    if (*title != '\0') {
	fprintf(fp, "set title \"%s\"\n", title);
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

static const char *front_strip (const char *s)
{
    while (*s) {
	if (isspace(*s) || *s == '{') {
	    s++;
	} else {
	    break;
	}
    }
	
    return s;
}

static void line_out (const char *s, int len, FILE *fp)
{
    char *p = malloc(len + 1);

    if (p != NULL) {
	*p = '\0';
	strncat(p, s, len);
	fprintf(fp, "%s\n", front_strip(p));
	free(p);
    }
}

void print_gnuplot_literal_lines (const char *s, FILE *fp)
{
    const char *p;

    if (s == NULL || *s == '\0') {
	return;
    }

    p = s = front_strip(s);

    fputs("# start literal lines\n", fp);

    while (*s && *s != '}') {
	if (*s == ';') {
	    line_out(p, s - p, fp);
	    p = s + 1;
	}
	s++;
    }

    fputs("# end literal lines\n", fp);
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
    char title[96];
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

    fp = open_gp_stream(PLOT_REGULAR, gi->flags, &err);
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

	sprintf(title, _("%s versus %s (with loess fit)"), s1, s2);
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	fprintf(fp, "set title \"%s\"\n", title);
	print_axis_label('y', s1, fp);
	print_axis_label('x', s2, fp);
    } else {
	print_keypos_string(GP_KEY_LEFT_TOP, fp);
	print_axis_label('y', series_get_graph_name(dset, yno), fp);
    }

    print_auto_fit_string(PLOT_FIT_LOESS, fp);

    if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 title '' w points, \\\n", fp);
    sprintf(title, _("loess fit, d = %d, q = %g"), d, q);
    fprintf(fp, " '-' using 1:2 title \"%s\" w lines\n", title);

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
			    char *targ)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    const double *yvar, *xvar = NULL;
    double x0, s2 = 0, *ps2 = NULL;
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
	char title[MAXTITLE], formula[GP_MAXFORMULA];
	double pd = dset->pd;

	if (gi->fit == PLOT_FIT_LOGLIN) {
	    b->val[0] += s2 / 2;
	}

	set_plotfit_line(title, formula, gi->fit, b->val, x0, pd);
	sprintf(targ, "%s title \"%s\" w lines\n", formula, title);
	gi->flags |= GPT_AUTO_FIT;
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

/* support the "fit" options for a single time-series plot */

static int time_fit_plot (gnuplot_info *gi, const char *literal,
			  const DATASET *dset)
{
    int yno = gi->list[1];
    const double *yvar = dset->Z[yno];
    FILE *fp = NULL;
    char fitline[128] = {0};
    PRN *prn;
    int t, err = 0;

    if (gi->x == NULL) {
	return E_DATA;
    }

    err = graph_list_adjust_sample(gi->list, gi, dset, 1);
    if (err) {
	return err;
    }

    err = get_fitted_line(gi, dset, fitline);
    if (err) {
	return err;
    }

    gi->flags |= GPT_LETTERBOX;

    fp = open_gp_stream(PLOT_REGULAR, gi->flags, &err);
    if (err) {
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

    if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 title '' w lines, \\\n", fp);
    
    gretl_push_c_numeric_locale();

    fprintf(fp, " %s", fitline);

    for (t=gi->t1; t<=gi->t2; t++) {
	fprintf(fp, "%.10g %.10g\n", gi->x[t], yvar[t]);
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
	    sprintf(gi->xtics, "%.*g %#.6g", d+1, vmin, 
		    (vmax - vmin)/ 4.0);
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

static void 
print_x_range (gnuplot_info *gi, const double *x, FILE *fp)
{
    if (gretl_isdummy(gi->t1, gi->t2, x)) {
	fputs("set xrange [-1:2]\n", fp);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
	gi->xrange = 3;
    } else {
	double xmin0, xmin, xmax0, xmax;

	gretl_minmax(gi->t1, gi->t2, x, &xmin0, &xmax0);
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

static void 
print_x_range_from_dates (gnuplot_info *gi, const DATASET *dset,
			  FILE *fp)
{
    char obs[OBSLEN];
    double xmin, xmax;

    ntodate(obs, gi->t1, dset);
    xmin = gnuplot_time_from_date(obs, gi->timefmt);
    ntodate(obs, gi->t2, dset);
    xmax = gnuplot_time_from_date(obs, gi->timefmt);

    gi->xrange = xmax - xmin;
    xmin -= gi->xrange * .025;
    xmax += gi->xrange * .025;

    fprintf(fp, "set xrange [%.12g:%.12g]\n", xmin, xmax);
    gi->xrange = xmax - xmin;
}

/* two or more y vars plotted against some x: test to see if we want
   to use two y axes */

static void
check_for_yscale (gnuplot_info *gi, const double **Z, int *oddman)
{
    double ymin[6], ymax[6];
    double ratio;
    int i, j, oddcount;

#if GP_DEBUG
    fprintf(stderr, "gnuplot: doing check_for_yscale\n");
#endif

    /* find minima, maxima of the y-axis vars */
    for (i=1; i<gi->list[0]; i++) {
	gretl_minmax(gi->t1, gi->t2, Z[gi->list[i]], 
		     &ymin[i-1], &ymax[i-1]);
    }

    gi->flags &= ~GPT_Y2AXIS;

    for (i=1; i<gi->list[0]; i++) {
	oddcount = 0;
	for (j=1; j<gi->list[0]; j++) {
	    if (j == i) {
		continue;
	    }
	    ratio = ymax[i-1] / ymax[j-1];
	    if (ratio > 5.0 || ratio < 0.2) {
		gi->flags |= GPT_Y2AXIS;
		oddcount++;
	    }
	}
	if (oddcount == gi->list[0] - 2) {
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
		fprintf(fp, "%.10g ?\n", xt);
	    } else {
		fprintf(fp, "%.10g %.10g", xt, yt);
		if (!(gi->flags & GPT_TS)) {
		    if (dset->markers) {
			fprintf(fp, " # %s", dset->S[t]);
		    } else if (dataset_is_time_series(dset)) {
			char obs[OBSLEN];

			ntodate(obs, t, dset);
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

    ntodate(obs, t, dset);
    sscanf(obs, "%d:%d", &maj, &min);
    if (maj > 1 && min == 1) {
	fprintf(fp, "%g ?\n", t + 0.5);
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
	    if (gi->withlist[0] == W_IMPULSES) {
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
	    if (gi->withlist[0] == W_LINES) {
		return 1;
	    }
	}
    }

    return 0;
}

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

    if (gi->flags & GPT_TIMEFMT) {
	lmax = gi->list[0];
	datlist[0] = 1;
	ynum = 1;
    } else if (gi->x != NULL) {
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
	    const char *date = NULL;
	    char obs[OBSLEN];

	    if (gi->flags & GPT_TIMEFMT) {
		ntodate(obs, t, dset);
		date = obs;
	    } else if (gi->x == NULL && 
		all_graph_data_missing(gi->list, t, (const double **) dset->Z)) {
		continue;
	    }

	    if (!(gi->flags & GPT_TS) && i == 1) {
		if (dset->markers) {
		    label = dset->S[t];
		} else if (!nomarkers && dataset_is_time_series(dset)) {
		    ntodate(obs, t, dset);
		    label = obs;
		}
	    }

	    if ((gi->flags & GPT_TS) && dset->structure == STACKED_TIME_SERIES) {
		maybe_print_panel_jot(t, dset, fp);
	    }

	    printvars(fp, t, datlist, (const double **) dset->Z, 
		      gi->x, label, date, xoff);
	}

	fputs("e\n", fp);
    }
}

static int
gpinfo_init (gnuplot_info *gi, gretlopt opt, const int *list, 
	     const char *literal, const DATASET *dset)
{
    int l0 = list[0];

    gi->withlist = NULL;

    get_gp_flags(gi, opt, list, dset);

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
    gi->timefmt[0] = '\0';
    gi->xtics[0] = '\0';
    gi->xfmt[0] = '\0';
    gi->yfmt[0] = '\0';
    gi->yformula = NULL;

    gi->x = NULL;
    gi->list = NULL;
    gi->dvals = NULL;

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
	&& !(gi->flags & GPT_DUMMY)  & !(opt & OPT_Y)) { 
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
    print_gnuplot_flags(gi->flags, 1);
#endif

    return 0;
}

static void clear_gpinfo (gnuplot_info *gi)
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

static void set_lwstr (const DATASET *dset, int v, char *s)
{
    if (default_png_scale > 1.0) {
	strcpy(s, " lw 2");
    } else {
	*s = '\0';
    }
}

static void set_withstr (gnuplot_info *gi, int i, char *str)
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
	} else {
	    strcpy(str, "w points");
	}
    } else if (gi->flags & GPT_LINES) {
        strcpy(str, "w lines");
    } else {
	strcpy(str, "w points");
    }
}

static int graph_list_adjust_sample (int *list, 
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

static int maybe_add_plotx (gnuplot_info *gi, int time_fit,
			    const DATASET *dset)
{
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

#if USE_TIMEFMT
    if (!time_fit) {
	if (dated_daily_data(dset) || dated_weekly_data(dset)) {
	    if (sample_size(dset) < 300) {
		/* experimental */
		gi->flags |= GPT_TIMEFMT;
		return 0;
	    }
	}
    }
#endif

    gi->x = gretl_plotx(dset, OPT_NONE);
    if (gi->x == NULL) {
	return E_ALLOC;
    }

    /* a bit ugly, but add a dummy list entry for
       the 'virtual' plot variable */
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

/* special tics for time series plots */

static void make_time_tics (gnuplot_info *gi,
			    const DATASET *dset,
			    int many, char *xlabel,
			    PRN *prn)
{
    if (many) {
	pprintf(prn, "# multiple timeseries %d\n", dset->pd);
    } else {
	pprintf(prn, "# timeseries %d", dset->pd);
	gi->flags |= GPT_LETTERBOX;
	pputs(prn, " (letterbox)\n");
    }

    if (gi->flags & GPT_TIMEFMT) {
	if (gnuplot_version() < 4.7) {
	    pputs(prn, "set xdata time # ZERO_YEAR=2000\n");
	} else {
	    pputs(prn, "set xdata time\n");
	}
	strcpy(gi->timefmt, "%Y-%m-%d");
	pprintf(prn, "set timefmt x \"%s\"\n", gi->timefmt);
	strcpy(gi->xfmt, "%Y-%m-%d");
	pprintf(prn, "set format x \"%s\"\n", gi->xfmt);
	return;
    }

    if (dset->pd == 4 && (gi->t2 - gi->t1) / 4 < 8) {
	pputs(prn, "set xtics nomirror 0,1\n"); 
	pputs(prn, "set mxtics 4\n");
    } else if (dset->pd == 12 && (gi->t2 - gi->t1) / 12 < 8) {
	pputs(prn, "set xtics nomirror 0,1\n"); 
	pputs(prn, "set mxtics 12\n");
    } else if (dated_daily_data(dset) || dated_weekly_data(dset)) {
	make_calendar_tics(dset, gi, prn);
    } else if (multiple_groups(dset, gi->t1, gi->t2)) {
	make_panel_unit_tics(dset, gi, prn);
	if (xlabel != NULL) {
	    strcpy(xlabel, _("time series by group"));
	}
    }
}

/* Respond to use of the option --matrix=<matname> in the gnuplot
   command, or create a plot directly from a matrix and a plot list.
*/

int matrix_plot (gretl_matrix *m, const int *list, const char *literal, 
		 gretlopt opt)
{
    DATASET *dset = NULL;
    int *plotlist = NULL;
    int err = 0;

    if (gretl_is_null_matrix(m)) {
	return E_DATA;
    }

    if (list != NULL && list[0] == 0) {
	dset = gretl_dataset_from_matrix(m, NULL, OPT_B, &err);
    } else {
	dset = gretl_dataset_from_matrix(m, list, OPT_B, &err);
    }
 
    if (err) {
	return err;
    }

    plotlist = gretl_consecutive_list_new(1, dset->v - 1);
    if (plotlist == NULL) {
	err = E_ALLOC;
    } 

    if (!err) {
	opt &= ~OPT_X;
	err = gnuplot(plotlist, literal, dset, opt);
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

#define fit_opts (OPT_I | OPT_L | OPT_Q | OPT_N | OPT_E)

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
    char fit_line[128] = {0};
    int time_fit = 0;
    int oddman = 0;
    int lmin, many = 0;
    PlotType ptype;
    gnuplot_info gi;
    int i, err = 0;

    gretl_error_clear();

    err = incompatible_options(opt, fit_opts);
    if (err) {
	return err;
    }

    if ((opt & OPT_T) && (opt & fit_opts)) {
	if (plotlist[0] > 1 || !dataset_is_time_series(dset)) {
	    return E_BADOPT;
	} else {
	    time_fit = 1;
	}
    }

    if (dataset_is_panel(dset) && 
	plotlist_is_group_invariant(plotlist, dset)) {
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

    err = maybe_add_plotx(&gi, time_fit, dset);
    if (err) {
	goto bailout;
    }

    if (gi.fit == PLOT_FIT_LOESS) {
	return loess_plot(&gi, literal, dset);
    }

    if (time_fit) {
	return time_fit_plot(&gi, literal, dset);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
	many = 1;
    }

    /* convenience pointer */
    list = gi.list;

    /* below: did have "height 1 width 1 box" for win32,
       "width 1 box" otherwise */
    strcpy(keystr, "set key left top\n");

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GPT_TS)) {
	int v = (gi.flags & GPT_DUMMY)? list[2] : list[list[0]];

	strcpy(xlabel, series_get_graph_name(dset, v));
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) {
	goto bailout;
    }

    /* the minimum number of plotlist elements */
    lmin = (gi.flags & GPT_TIMEFMT) ? 1 : 2;

    /* adjust sample range, and reject if it's empty */
    err = graph_list_adjust_sample(list, &gi, dset, lmin);
    if (err) {
	goto bailout;
    }

    /* add a simple regression line if appropriate */
    if (!use_impulses(&gi) && !(gi.flags & GPT_FIT_OMIT) && list[0] == 2 && 
	!(gi.flags & GPT_TS) && !(gi.flags & GPT_RESIDS)) {
	const char *xname = dset->varname[list[2]];
	const char *yname = dset->varname[list[1]];

	get_fitted_line(&gi, dset, fit_line);
	if (*xname != '\0' && yname != '\0') {
	    pprintf(prn, "# X = '%s' (%d)\n", xname, list[2]);
	    pprintf(prn, "# Y = '%s' (%d)\n", yname, list[1]);
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

    /* open file and dump the prn into it: we delay writing
       the file header till we know a bit more about the plot
    */
    fp = open_gp_stream(ptype, gi.flags, &err);
    if (err) {
	gretl_print_destroy(prn);
	goto bailout;
    } 

    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    print_axis_label('x', xlabel, fp);
    fputs("set xzeroaxis\n", fp); 
    gnuplot_missval_string(fp);

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
	if (gi.flags & GPT_RESIDS && !(gi.flags & GPT_DUMMY)) { 
	    make_gtitle(&gi, GTITLE_RESID, series_get_label(dset, list[1]), 
			NULL, fp);
	    fprintf(fp, "set ylabel '%s'\n", _("residual"));
	} else if (gi.flags & GPT_TIMEFMT) {
	    no_key = 0;
	} else {
	    print_axis_label('y', series_get_graph_name(dset, list[1]), fp);
	}
	if (no_key) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((gi.flags & GPT_RESIDS) && (gi.flags & GPT_DUMMY)) { 
	make_gtitle(&gi, GTITLE_RESID, series_get_label(dset, list[1]), 
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
    } 

    if (many) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    if (gi.flags & GPT_TIMEFMT) {
	print_x_range_from_dates(&gi, dset, fp);
    } else if (gi.x != NULL) {
	print_x_range(&gi, gi.x, fp);
    } else {
	print_x_range_from_list(&gi, dset, list, fp);
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
    } else if (literal != NULL && *literal != '\0') {
	print_gnuplot_literal_lines(literal, fp);
    }

    /* now print the 'plot' lines */
    fputs("plot \\\n", fp);
    if (gi.flags & GPT_Y2AXIS) {
	/* using two y axes */
	for (i=1; i<list[0]; i++) {
	    set_lwstr(dset, list[i], lwstr);
	    set_withstr(&gi, i, withstr);
	    fprintf(fp, "'-' using 1:($2) axes %s title \"%s (%s)\" %s%s%s",
		    (i == oddman)? "x1y2" : "x1y1",
		    series_get_graph_name(dset, list[i]), 
		    (i == oddman)? _("right") : _("left"),
		    withstr,
		    lwstr,
		    (i == list[0] - 1)? "\n" : ", \\\n");
	}
    } else if (gi.flags & GPT_DUMMY) { 
	/* plot shows separation by discrete variable */
	int nd = gretl_vector_get_length(gi.dvals);

	strcpy(s1, (gi.flags & GPT_RESIDS)? _("residual") : 
	       series_get_graph_name(dset, list[1]));
	strcpy(s2, series_get_graph_name(dset, list[3]));
	for (i=0; i<nd; i++) {
	    fprintf(fp, " '-' using 1:($2) title \"%s (%s=%g)\" w points ", 
		    s1, s2, gretl_vector_get(gi.dvals, i));
	    if (i < nd - 1) {
		fputs(", \\\n", fp);
	    } else {
		fputc('\n', fp);
	    }
	}
    } else if (gi.yformula != NULL) {
	/* we have a formula to plot, not just data */
	fprintf(fp, " '-' using 1:($2) title \"%s\" w points , \\\n", _("actual"));	
	fprintf(fp, "%s title '%s' w lines\n", gi.yformula, _("fitted"));
    } else if (gi.flags & GPT_FA) {
	/* this is a fitted vs actual plot */
	set_withstr(&gi, 1, withstr);
	fprintf(fp, " '-' using 1:($2) title \"%s\" %s lt 2, \\\n", _("fitted"), withstr);
	fprintf(fp, " '-' using 1:($2) title \"%s\" %s lt 1\n", _("actual"), withstr);	
    } else {
	/* all other cases */
	int lmax = (gi.flags & GPT_TIMEFMT)? list[0] : list[0] - 1;

	for (i=1; i<=lmax; i++)  {
	    set_lwstr(dset, list[i], lwstr);
	    if (list[0] == 2 && !(gi.flags & GPT_TIMEFMT)) {
		*s1 = '\0';
	    } else {
		strcpy(s1, series_get_graph_name(dset, list[i]));
	    }
	    set_withstr(&gi, i, withstr);
	    fprintf(fp, " '-' using 1:($2) title \"%s\" %s%s", s1, withstr, lwstr);
	    if (i < lmax || (gi.flags & GPT_AUTO_FIT)) {
	        fputs(" , \\\n", fp); 
	    } else {
	        fputc('\n', fp);
	    }
	}
    } 

    if (*fit_line != '\0') {
        fputs(fit_line, fp);
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

    fp = open_gp_stream(PLOT_REGULAR, gi.flags, &err);
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
    fputs(" '-' using 1:($2) notitle w points , \\\n", fp);
    fprintf(fp, " x title \"%s\" w lines\n", _("actual = predicted"));

    print_gp_data(&gi, dset, fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    clear_gpinfo(&gi);

    return err;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @dset: dataset struct.
 * @opt: can include %OPT_L to use lines, %OPT_U to
 * direct output to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the variables in @list, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const int *list, const DATASET *dset, 
		    gretlopt opt)
{
    GptFlags flags = 0;
    int xvar = 0, yvar = 0;
    const double *x = NULL;
    const double *y = NULL;
    const double *obs = NULL;
    int rows, cols, tseries = 0;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

    if (opt & OPT_L) {
	flags |= GPT_LINES;
    }

    pos = gretl_list_separator_position(list);

    if (pos == 0) {
	/* plot against time or index */
	obs = gretl_plotx(dset, OPT_NONE);
	if (obs == NULL) {
	    return E_ALLOC;
	}
	plotlist = gretl_list_copy(list);
	flags |= GPT_LINES;
	if (dataset_is_time_series(dset)) {
	    tseries = 1;
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
    get_multiplot_layout(nplots, tseries, &rows, &cols);
    maybe_set_small_font(nplots);

    if (nplots > 12) {
	flags |= GPT_XXL;
    } else if (nplots > 9) {
	flags |= GPT_XL;
    }

    fp = open_gp_stream(PLOT_MULTI_SCATTER, flags, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "set multiplot layout %d,%d\n", rows, cols);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (obs != NULL) {
	double startdate = obs[dset->t1];
	double enddate = obs[dset->t2];
	int jump, T = dset->t2 - dset->t1 + 1;

	fprintf(fp, "set xrange [%g:%g]\n", floor(startdate), ceil(enddate));

	if (dset->pd == 1) {
	    jump = T / 6;
	} else {
	    jump = T / (4 * dset->pd);
	}

	fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);
    } else {
	fputs("set noxtics\nset noytics\n", fp);
    }

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
	    double xt, yt;

	    xt = yvar ? dset->Z[j][t] : xvar ? x[t] : obs[t];
	    yt = yvar ? y[t] : dset->Z[j][t];

	    if (na(xt)) {
		fputs("? ", fp);
	    } else {
		fprintf(fp, "%.10g ", xt);
	    }

	    if (na(yt)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.10g\n", yt);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

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
	    strncat(s, "~", 1);
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
 * matrix_scatters:
 * @m: matrix containing data to plot.
 * @list: list of columns to plot, or NULL.
 * @dset: dataset pointer, or NULL.
 * @opt: can include %OPT_L to use lines, %OPT_U to
 * direct output to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the data in @m, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int matrix_scatters (const gretl_matrix *m, const int *list, 
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

    if (opt & OPT_L) {
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

    fp = open_gp_stream(PLOT_MULTI_SCATTER, flags, &err);
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
	int jump, T = t2 - t1 + 1;

	fprintf(fp, "set xrange [%g:%g]\n", floor(startdate), ceil(enddate));

	jump = (pd == 1)? (T / 6) : (T / (4 * pd));
	fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);
    } else if (simple_obs) {
	fprintf(fp, "set xrange [0:%d]\n", m->rows - 1);
	fprintf(fp, "set xtics 0, %d\n", m->rows / 6);
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
	    double xt, yt;
	    int s = t - t1;

	    xt = ycol ? zj[s] : xcol ? x[s] : get_obsx(obs, t, s);
	    yt = (y != NULL)? y[s] : zj[s];

	    if (xna(xt)) {
		fputs("? ", fp);
	    } else {
		fprintf(fp, "%.10g ", xt);
	    }

	    if (xna(yt)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.10g\n", yt);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

    free(plotlist);

    return finalize_plot_input_file(fp);
}

static int get_3d_output_file (FILE **fpp)
{
    char fname[MAXLEN];
    int err = 0;

    sprintf(fname, "%sgpttmp.plt", gretl_dotdir());
    *fpp = gretl_fopen(fname, "w");

    if (*fpp == NULL) {
	err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
    }

    return err;
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
	(snedecor_cdf_comp(smod.dfn, smod.dfd, smod.fstt) < .10 || (opt & OPT_F))) {
	double uadj = (umax - umin) * 0.02;
	double vadj = (vmax - vmin) * 0.02;

	ret = g_strdup_printf("[u=%g:%g] [v=%g:%g] "
			      "%g+(%g)*u+(%g)*v title ''", 
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
 * @opt: may include OPT_F to force display of fitted surface.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		DATASET *dset, gretlopt opt)
{
    FILE *fq = NULL;
    int t, t1 = dset->t1, t2 = dset->t2;
    int orig_t1 = dset->t1, orig_t2 = dset->t2;
    int lo = list[0];
    int datlist[4];
    int addstyle = 0;
    gchar *surface = NULL;

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return E_DATA;
    }

    if (get_3d_output_file(&fq)) {
	return E_FOPEN;
    }

    list_adjust_sample(list, &t1, &t2, dset, NULL);

    /* if resulting sample range is empty, complain */
    if (t1 >= t2) {
	fclose(fq);
	return E_MISSDATA;
    }

    dset->t1 = t1;
    dset->t2 = t2;

#ifndef WIN32
    if (gnuplot_has_wxt()) {
	fputs("set term wxt\n", fq);
    } else if (gnuplot_has_x11()) {
	fputs("set term x11\n", fq);
    } else if (gnuplot_has_aqua()) {
	/* can't do rotation, but it's all we have */
	fputs("set term aqua\n", fq);
    } else {
	fclose(fq);
	return E_EXTERNAL;
    }
#endif

    gretl_push_c_numeric_locale();

    /* try to ensure we don't get "invisible" green datapoints */
    fprintf(fq, "set style line 2 lc rgb \"#0000ff\"\n");
    addstyle = 1;
    
    print_axis_label('x', series_get_graph_name(dset, list[2]), fq);
    print_axis_label('y', series_get_graph_name(dset, list[1]), fq);
    print_axis_label('z', series_get_graph_name(dset, list[3]), fq);

    gnuplot_missval_string(fq);

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

    surface = maybe_get_surface(list, dset, opt);

    if (surface != NULL) {
	if (addstyle) {
	    fprintf(fq, "splot %s, \\\n'-' title '' w p ls 2\n", surface);
	} else {
	    fprintf(fq, "splot %s, \\\n'-' title '' w p lt 3\n", surface);
	}
	g_free(surface);
    } else {
	if (addstyle) {
	    fputs("splot '-' title '' w p ls 2\n", fq);
	} else {
	    fputs("splot '-' title '' w p lt 3\n", fq);
	}
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
	printvars(fq, t, datlist, (const double **) dset->Z, 
		  NULL, label, NULL, 0.0);
    }	
    fputs("e\n", fq);

    gretl_pop_c_numeric_locale();

    dset->t1 = orig_t1;
    dset->t2 = orig_t2;

    fclose(fq);

    return 0;
}

static void print_freq_test_label (char *s, int teststat, 
				   double v, double pv)
{
    gretl_pop_c_numeric_locale();
    if (teststat == GRETL_STAT_Z) {
	sprintf(s, "z = %.3f [%.4f]", v, pv);
    } else if (teststat == GRETL_STAT_NORMAL_CHISQ) {
	sprintf(s, "%s(2) = %.3f [%.4f]", _("Chi-square"), v, pv);
    }
    gretl_push_c_numeric_locale();
}

static void print_freq_dist_label (char *s, int dist, double x, double y)
{
    int dcomma = 0;
    char test[8];

    gretl_pop_c_numeric_locale();

    sprintf(test, "%g", 0.5);
    if (strchr(test, ',')) {
	dcomma = 1;
    }

    if (dist == D_NORMAL) {
	sprintf(s, "N(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    } else if (dist == D_GAMMA) {
	sprintf(s, "gamma(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    }

    gretl_push_c_numeric_locale();
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

static double minskip (FreqDist *freq)
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

int plot_freq (FreqDist *freq, DistCode dist)
{
    double alpha = 0.0, beta = 0.0, lambda = 1.0;
    FILE *fp = NULL;
    int i, K = freq->numbins;
    char withstr[32] = {0};
    char label[80] = {0};
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth;
    const double *endpt;
    int plottype, use_boxes = 1;
    int err = 0;

    if (K == 0) {
	return E_DATA;
    }

    if (K == 1) {
	gretl_errmsg_sprintf(_("'%s' is a constant"), freq->varname);
	return E_DATA;
    }

    if (dist == D_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == D_GAMMA) {
	plottype = PLOT_FREQ_GAMMA;
    } else {
	plottype = PLOT_FREQ_SIMPLE;
    }

    fp = open_plot_input_file(plottype, &err);
    if (err) {
	return err;
    }  

#if GP_DEBUG
    fprintf(stderr, "*** plot_freq called\n");
#endif  

    if (freq->discrete) {
	endpt = freq->midpt;
	barwidth = minskip(freq); 
	use_boxes = 0;
    } else {
	/* equally sized bins, width to be determined */
	endpt = freq->endpt;
	barwidth = freq->endpt[K-1] - freq->endpt[K-2];
    }

    gretl_push_c_numeric_locale();

    if (dist) {
	lambda = 1.0 / (freq->n * barwidth);

	if (dist == D_NORMAL) {
	    fputs("# literal lines = 4\n", fp);
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
		print_freq_test_label(label, GRETL_STAT_NORMAL_CHISQ, freq->test, 
				      chisq_cdf_comp(2, freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n", 
			label);
	    }	
	} else if (dist == D_GAMMA) {
	    double var = freq->sdx * freq->sdx;

	    /* scale param = variance/mean */
	    beta = var / freq->xbar;
	    /* shape param = mean/scale */
	    alpha = freq->xbar / beta;

	    fputs("# literal lines = 4\n", fp);
	    fprintf(fp, "beta = %g\n", beta);
	    fprintf(fp, "alpha = %g\n", alpha);
	    plotmin = 0.0;
	    plotmax = freq->xbar + 4.0 * freq->sdx;

	    if (!na(freq->test)) {
		fprintf(fp, "set label '%s:' at graph .03, graph .97 front\n",
			_("Test statistic for gamma"));
		print_freq_test_label(label, GRETL_STAT_Z, freq->test, 
				      normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93 front\n", 
			label);
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
	plotmin = freq->midpt[0] - barwidth;
	plotmax = freq->midpt[K-1] + barwidth;
	fprintf(fp, "set xrange [%.10g:%.10g]\n", plotmin, plotmax);
	maybe_set_yrange(freq, lambda, fp);
	fputs("set nokey\n", fp);
    }

    fprintf(fp, "set xlabel '%s'\n", freq->varname);
    if (dist) {
	fprintf(fp, "set ylabel '%s'\n", _("Density"));
    } else {
	fprintf(fp, "set ylabel '%s'\n", _("Relative frequency"));
    }

    if (isnan(lambda)) {
	if (fp != NULL) {
	    fclose(fp);
	}
	return 1;
    }

    /* plot instructions */
    if (use_boxes) {
	fputs("set style fill solid 0.6\n", fp);
	strcpy(withstr, "w boxes");
    } else {
	strcpy(withstr, "w impulses linewidth 3");
    }

    if (!dist) {
	fprintf(fp, "plot '-' using 1:2 %s\n", withstr);
    } else if (dist == D_NORMAL) {
	print_freq_dist_label(label, dist, freq->xbar, freq->sdx);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title \"%s\" %s, \\\n"
		"1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) "
		"title \"%s\" w lines\n",
		freq->varname, withstr, label);
    } else if (dist == D_GAMMA) {
	print_freq_dist_label(label, dist, alpha, beta);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s, \\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) "
		"title \"%s\" w lines\n",
		freq->varname, withstr, label); 
    }

    for (i=0; i<K; i++) { 
	fprintf(fp, "%.10g %.10g\n", freq->midpt[i], lambda * freq->f[i]);
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static void print_y_data (const double *x, 
			  const double *y,
			  const int *order, 
			  int t0, int t1, int t2, 
			  FILE *fp)
{
    int i, t, n = t2 - t0 + 1;
    double xt;

    for (i=0; i<n; i++) {
	if (order != NULL) {
	    t = order[i];
	    xt = i + 1;
	} else {
	    t = t0 + i;
	    xt = x[t];
	}	
	if (t < t1 || na(y[t])) {
	    fprintf(fp, "%.10g ?\n", xt);
	} else {
	    fprintf(fp, "%.10g %.10g\n", xt, y[t]);
	}
    }

    fputs("e\n", fp);
}

enum {
    CONF_BARS,
    CONF_FILL,
    CONF_LOW,
    CONF_HIGH
};

static void print_confband_data (const double *x, const double *y,
				 const double *e, const int *order,
				 int t0, int t1, int t2, 
				 int mode, FILE *fp)
{
    int i, t, n = t2 - t0 + 1;
    double xt;

    for (i=0; i<n; i++) {
	if (order != NULL) {
	    t = order[i];
	    xt = i + 1;
	} else {
	    t = t0 + i;
	    xt = x[t];
	}
	if (t < t1 || na(y[t]) || na(e[t])) {
	    if (mode == CONF_LOW || mode == CONF_HIGH) {
		fprintf(fp, "%.10g ?\n", xt);
	    } else {
		fprintf(fp, "%.10g ? ?\n", xt);
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
		fprintf(fp, "%.10g ?\n", x[t]);
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

/* To get a somewhat more meangingful plot of actual and fitted
   with error bars for cross-sectional data, we consider
   ordering the observations by increasing value of the
   dependent variable. We first check to see if the data are
   already in sort-order; if not, we produce an array of
   ints that specifies the sorted order -- and if there are
   NAs for the dependent variable within the forecast range,
   we adjust the end-point so as to omit the useless
   observations.
*/

static int *get_sorted_fcast_order (const FITRESID *fr, 
				    int t1, int *pt2)
{
    int *order = NULL;
    struct fsorter *fs;
    int nmiss = 0, sorted = 1;
    int t2 = *pt2;
    int n = t2 - t1 + 1;
    int i, t;

    for (t=t1; t<=t2; t++) {
	if (na(fr->actual[t])) {
	    nmiss++;
	} else if (t < t2 && fr->actual[t] > fr->actual[t+1]) {
	    sorted = 0;
	}
    }

    if (sorted) {
	/* OK, leave well alone */
	return NULL;
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
	fs[i].y = fr->actual[t];
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

/* note: if @opt includes OPT_H, that says to show fitted 
   values for the pre-forecast range
*/

int plot_fcast_errs (const FITRESID *fr, const double *maxerr,
		     const DATASET *dset, gretlopt opt)
{
    FILE *fp = NULL;
    const double *obs = NULL;
    int *order = NULL;
    GptFlags flags = 0;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int use_fill = 0, use_lines = 0;
    int do_errs = (maxerr != NULL);
    char cistr[64];
    int t2 = fr->t2;
    int t1, yhmin;
    int t, n, err = 0;

    /* note: yhmin is the first obs at which to start plotting y-hat */
    if (do_errs) {
	t1 = fr->t0;
	yhmin = (opt & OPT_H)? fr->t0 : fr->t1;
    } else {
	t1 = (fr->t0 >= 0)? fr->t0 : 0;
	yhmin = t1;
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

    fp = open_gp_stream(PLOT_FORECAST, flags, &err);
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
	fprintf(fp, "# timeseries %d\n", dset->pd);
    } else if (dataset_is_cross_section(dset) && yhmin == t1) {
	order = get_sorted_fcast_order(fr, t1, &t2);
    } 

    if (!dataset_is_time_series(dset)) {
	if (order != NULL) {
	    gchar *text = g_strdup_printf(_("observations sorted by %s"),
					  fr->depvar);
	    fputs("unset xtics\n", fp);
	    fprintf(fp, "set xlabel \"%s\"\n", text);
	} else if (n < 33) {
	    fputs("set xtics 1\n", fp);
	}
    }

    if (do_errs && !use_fill && !use_lines && n > 150) {
	use_fill = 1;
    }

    if (use_fill) {
	fprintf(fp, "set style fill solid 0.4\n");
    }

    fputs("set key left top\n", fp);
    fputs("plot \\\n", fp);

    if (do_errs) {
	sprintf(cistr, _("%g percent interval"), 100 * (1 - fr->alpha));
    }

    if (use_fill) {
	/* plot the confidence bands first so the other lines
	   come out on top */
	if (do_errs) {
	    fprintf(fp, "'-' using 1:2:3 title '%s' w filledcurve lt 3 , \\\n",
		    cistr);
	} 
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1 , \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines lt 2\n", _("forecast"));
    } else {
	/* plot confidence bands last */
	if (depvar_present) {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines , \\\n",
		    fr->depvar);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines", _("forecast"));
	if (do_errs) {
	    if (use_lines) {
		fprintf(fp, " , \\\n'-' using 1:2 title '%s' w lines , \\\n",
			cistr);
		fputs("'-' using 1:2 notitle '%s' w lines lt 3\n", fp);
	    } else {
		fprintf(fp, " , \\\n'-' using 1:2:3 title '%s' w errorbars\n",
			cistr);
	    }
	} else {
	    fputc('\n', fp);
	}
    }

    gretl_push_c_numeric_locale();

    /* write out the inline data, the order depending on whether
       or not we're using fill style for the confidence
       bands
    */

    if (use_fill) {
	if (do_errs) {
	    print_confband_data(obs, fr->fitted, maxerr, order,
				t1, yhmin, t2, CONF_FILL, fp);
	}
	if (depvar_present) {
	    print_y_data(obs, fr->actual, order, t1, t1, t2, fp);
	}
	print_y_data(obs, fr->fitted, order, t1, yhmin, t2, fp);
    } else {
	if (depvar_present) {
	    print_y_data(obs, fr->actual, order, t1, t1, t2, fp);
	}
	print_y_data(obs, fr->fitted, order, t1, yhmin, t2, fp);
	if (do_errs) {
	    if (use_lines) {
		print_confband_data(obs, fr->fitted, maxerr, order,
				    t1, yhmin, t2, CONF_LOW, fp);
		print_confband_data(obs, fr->fitted, maxerr, order,
				    t1, yhmin, t2, CONF_HIGH, fp);
	    } else {
		print_confband_data(obs, fr->fitted, maxerr, order,
				    t1, yhmin, t2, CONF_BARS, fp);
	    }
	}	
    }

    gretl_pop_c_numeric_locale();

    if (order != NULL) {
	free(order);
    }

    return finalize_plot_input_file(fp);
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
	    fprintf(fp, "%.10g ?\n", x[t]);
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

    fp = open_gp_stream(PLOT_FORECAST, flags, &err);
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

    fputs("'-' using 1:2 notitle w points , \\\n", fp);
    fputs("'-' using 1:2 notitle w lines , \\\n", fp);
    fprintf(fp, "'-' using 1:2 title '%s' w lines , \\\n", cistr);
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

    fp = open_plot_input_file(PLOT_RQ_TAU, &err);
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

    fputs("set style fill solid 0.4\n", fp);

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
    fputs("'-' using 1:2:3 notitle w filledcurve lt 3 , \\\n", fp);

    /* rq estimates */
    tmp = g_strdup_printf(_("Quantile estimates with %g%% band"), cval);
    fprintf(fp, "'-' using 1:2 title '%s' w lp lt 1 , \\\n", tmp);
    g_free(tmp);

    /* numeric output coming up! */
    gretl_push_c_numeric_locale();

    /* ols estimate plus (1 - alpha) band */
    tmp = g_strdup_printf(_("OLS estimate with %g%% band"), cval);
    fprintf(fp, "%g title '%s' w lines lt 2 , \\\n", pmod->coeff[k], tmp);
    g_free(tmp);
    fprintf(fp, "%g notitle w dots lt 2 , \\\n", pmod->coeff[k] + olsband);
    fprintf(fp, "%g notitle w dots lt 2\n", pmod->coeff[k] - olsband);

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

    fp = open_plot_input_file(PLOT_GARCH, &err);
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
    void *handle = NULL;
    int err;

    range_mean_graph = get_plugin_function("range_mean_graph", &handle);
    if (range_mean_graph == NULL) {
	return 1;
    }

    err = range_mean_graph(list[1], dset, opt, prn);

    close_plugin(handle);

    return err;
}

int 
hurstplot (const int *list, DATASET *dset, gretlopt opt, PRN *prn)
{
    int (*hurst_exponent) (int, const DATASET *, gretlopt, PRN *);
    void *handle = NULL;
    int err;

    hurst_exponent = get_plugin_function("hurst_exponent", &handle);
    if (hurst_exponent == NULL) {
	return 1;
    }

    err = hurst_exponent(list[1], dset, opt, prn);

    close_plugin(handle);

    return err;
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
				const DATASET *dset,
				gretlopt opt)
{
    DATASET *gset;
    int nunits, T = dset->pd;
    int list[3] = {2, 1, 2};
    gchar *literal = NULL;
    gchar *title = NULL;
    const double *obs;
    int i, t, s, s0;
    int err = 0;

    nunits = panel_sample_size(dset);

    obs = gretl_plotx(dset, OPT_P);
    if (obs == NULL) {
	return E_ALLOC;
    }

    gset = create_auxiliary_dataset(3, T, 0);
    if (gset == NULL) {
	return E_ALLOC;
    }

    strcpy(gset->varname[1], dset->varname[vnum]);
    series_set_display_name(gset, 1, series_get_display_name(dset, vnum));

    s0 = dset->t1 * dset->pd;

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
	if (n == 0) {
	    gset->Z[1][t] = NADBL;
	} else {
	    gset->Z[1][t] = xsum / n;
	}
	gset->Z[2][t] = obs[t];
    }

    opt |= OPT_O; /* use lines */
    
    title = g_strdup_printf(_("mean %s"), 
			    series_get_graph_name(dset, vnum));
    literal = g_strdup_printf("set ylabel \"%s\" ; set xlabel ;", 
			      title);
    err = gnuplot(list, literal, gset, opt);

    g_free(title);
    g_free(literal);
    destroy_dataset(gset);

    return err;
}

/* Here we're trying to find out if the observation labels
   for a panel dataset are such that they uniquely identify
   the units/individuals (e.g. country or city names, 
   repeated for each time-series observation on the given
   unit).
*/

static int dataset_has_panel_labels (const DATASET *dset,
				     int *use, int *strip)
{
    int t, u, ubak = -1;
    int fail = 0;
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
	ubak = u;
    }

    if (!fail) {
	/* fine: the full obs labels satisfy the criterion */
	ret = 1;
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
		    strncat(s1, dset->S[t], strlen(dset->S[t]) - *strip);
		}
		if (u == ubak && strcmp(s1, s0)) {
		    /* same unit, different label: no */
		    fail = 1;
		} else if (ubak >= 0 && u != ubak && !strcmp(s1, s0)) {
		    /* different unit, same label: no */
		    fail = 2;
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
				  DATASET *dset,
				  gretlopt opt)
{
    DATASET *gset;
    int u0, nunits, T = dset->pd;
    int *list = NULL;
    gchar *literal = NULL;
    gchar *title = NULL;
    const double *obs = NULL;
    char const **grpnames = NULL;
    int nv, panel_labels = 0;
    int single_series;
    int use = 0, strip = 0;
    int i, t, s, s0;
    int err = 0;

    single_series = series_is_group_invariant(dset, vnum);

    if (single_series) {
	nunits = 1;
    } else {
	nunits = panel_sample_size(dset);
    }

    nv = nunits + 2;
    u0 = dset->t1 / T;

    obs = gretl_plotx(dset, OPT_P);
    if (obs == NULL) {
	return E_ALLOC;
    }

    gset = create_auxiliary_dataset(nv, T, 0);
    if (gset == NULL) {
	return E_ALLOC;
    }

    list = gretl_consecutive_list_new(1, nv - 1);
    if (list == NULL) {
	destroy_dataset(gset);
	return E_ALLOC;
    }

    if (!single_series) {
	grpnames = get_panel_group_names(dset);
	if (grpnames == NULL && dset->S != NULL) {
	    panel_labels = dataset_has_panel_labels(dset, &use, &strip);
	}
    }

    s0 = dset->t1 * dset->pd;

    for (i=0; i<nunits; i++) {
	s = s0 + i * T;
	if (single_series) {
	    strcpy(gset->varname[i+1], dset->varname[vnum]);
	} else if (grpnames != NULL) {
	    strncat(gset->varname[i+1], grpnames[u0+i], VNAMELEN-1);
	} else if (panel_labels) {
	    if (use > 0) {
		strncat(gset->varname[i+1], dset->S[s], use);
	    } else if (strip > 0) {
		strncat(gset->varname[i+1], dset->S[s],
			strlen(dset->S[s]) - strip);
	    } else {
		strcpy(gset->varname[i+1], dset->S[s]);
	    }
	} else {
	    sprintf(gset->varname[i+1], "%d", u0+i+1);
	}
	for (t=0; t<T; t++) {
	    gset->Z[i+1][t] = dset->Z[vnum][s++];
	}
    }

    for (t=0; t<T; t++) {
	gset->Z[nv-1][t] = obs[t];
    }	

    opt |= OPT_O; /* use lines */

    if (single_series) {
	opt |= OPT_S; /* suppress-fitted */
    } else {
	const char *gname = series_get_graph_name(dset, vnum);
	const char *vname = panel_group_names_varname(dset);

	if (vname != NULL) {
	    title = g_strdup_printf("%s by %s", gname, vname);
	} else {
	    title = g_strdup_printf("%s by group", gname);
	}
	literal = g_strdup_printf("set title \"%s\" ; set xlabel ;", title);
    }

    err = gnuplot(list, literal, gset, opt);

    g_free(title);
    g_free(literal);
    destroy_dataset(gset);
    free(list);

    return err;
}

/* Panel: plot one variable as a time series, with separate plots for
   each cross-sectional unit.  By default we arrange the plots in a
   grid, but if OPT_V is given we make each plot full width and
   stack the plots vertically on the "page".
*/

static int panel_grid_ts_plot (int vnum, DATASET *dset,
			       gretlopt opt) 
{
    FILE *fp = NULL;
    int w, rows, cols;
    const double *y, *x = NULL;
    const char *vname;
    char uname[OBSLEN];
    double xt, yt, ymin, ymax, incr;
    int u0, nunits, T = dset->pd;
    int panel_labels = 0;
    int use = 0, strip = 0;
    int i, s, t, t0;
    int err = 0;

    nunits = panel_sample_size(dset);
    u0 = dset->t1 / dset->pd;

    if (opt & OPT_V) {
	int xvar = plausible_panel_time_var(dset);

	if (xvar > 0) {
	    x = dset->Z[xvar];
	}
	cols = 1;
	rows = nunits;
    } else {
	get_multiplot_layout(nunits, 0, &rows, &cols);
    }

    if (rows == 0 || cols == 0) {
	return E_DATA;
    }

    maybe_set_small_font(nunits);

    fp = open_plot_input_file(PLOT_PANEL, &err);
    if (err) {
	return err;
    }

    if (dset->S != NULL) {
	panel_labels = dataset_has_panel_labels(dset, &use, &strip);
    }

    vname = dset->varname[vnum];
    y = dset->Z[vnum];
    gretl_minmax(dset->t1, dset->t2, y, &ymin, &ymax);
    w = panel_ytic_width(ymin, ymax);

    fputs("set key left top\n", fp);
    fputs("set datafile missing \"?\"\n", fp);
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
	if (panel_labels) {
	    *uname = '\0';
	    s = (u0 + i) * dset->pd;
	    if (use > 0) {
		strncat(uname, dset->S[s], use);
	    } else if (strip > 0) {
		strncat(uname, dset->S[s],
			strlen(dset->S[s]) - strip);
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

	fputs("plot \\\n'-' using 1:($2) notitle w lines\n", fp);

	for (t=0; t<T; t++) {
	    if (x != NULL) {
		xt = x[t+t0];
	    } else {
		xt = t + 1;
	    }
	    yt = y[t+t0];
	    if (na(yt)) {
		fprintf(fp, "%g ?\n", xt);
	    } else {
		fprintf(fp, "%g %.10g\n", xt, yt);
	    }
	}
	fputs("e\n", fp);

	t0 += T;
    }
    
    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

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
	return panel_means_ts_plot(vnum, dset, opt);
    } else {
	return panel_overlay_ts_plot(vnum, dset, opt);
    }
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
    char title[128];
    int t;

    if (!confint) {
	fputs("# impulse response plot\n", fp);
    }

    if (confint) {
	fputs("set key left top\n", fp);
	sprintf(title, _("response of %s to a shock in %s, "
			 "with bootstrap confidence interval"),
		targname, shockname);
    } else {
	fputs("set nokey\n", fp);
	sprintf(title, _("response of %s to a shock in %s"), 
		targname, shockname);
    }

    fprintf(fp, "set xlabel '%s'\n", perlabel);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [-1:%d]\n", periods);
    fprintf(fp, "set title '%s'\n", title);

    if (confint) {
	double ql = alpha / 2;
	double qh = 1.0 - ql;

	fputs("plot \\\n", fp);
	if (use_fill) {
	    sprintf(title, _("%g percent confidence band"), 100 * (1 - alpha));
	    fprintf(fp, "'-' using 1:2:3 title '%s' w filledcurve lt %d, \\\n", 
		    title, SHADECOLOR + 1);
	    if (data_straddle_zero(resp)) {
		fputs("0 notitle w lines lt 0, \\\n", fp);
	    }
	    fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1\n", _("point estimate"));
	} else {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n", 
		    _("point estimate"));
	    sprintf(title, _("%g and %g quantiles"), ql, qh);
	    fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n", title);
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

	fp = open_plot_input_file((confint)? PLOT_IRFBOOT : PLOT_REGULAR, &err);
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
    int err = 0;

    V = gretl_VAR_get_FEVD_matrix(var, targ, periods, dset, &err);
    if (V == NULL) {
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_REGULAR, &err);
    if (err) {
	gretl_matrix_free(V);
	return err;
    }

    histo = (opt & OPT_H)? 1 : 0;

    v = gretl_VAR_get_variable_number(var, targ);

    fprintf(fp, "set xlabel '%s'\n", dataset_period_label(dset));
    title = g_strdup_printf(_("forecast variance decomposition for %s"), 
			    dset->varname[v]);
    fprintf(fp, "set title '%s'\n", title);
    g_free(title);

    if (histo) {
	fputs("set key outside\n", fp);
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
	    fprintf(fp, "'-' using 2 title '%s'", dset->varname[v]);
	} else {
	    fprintf(fp, "'-' using 1:2 title '%s' w lines", dset->varname[v]);
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

int gretl_VAR_plot_multiple_irf (GRETL_VAR *var, 
				 int periods, double alpha,
				 const DATASET *dset,
				 gretlopt opt)
{
    FILE *fp = NULL;
    GptFlags flags = 0;
    int confint = 0;
    int use_fill = !(opt & OPT_E);
    char title[128];
    int n = var->neqns;
    int nplots = n * n;
    int t, i, j;
    int err = 0;

    maybe_set_small_font(nplots);

    if (nplots > 12) {
	flags |= GPT_XXL;
    } else if (nplots > 9) {
	flags |= GPT_XL;
    }

    fp = open_gp_stream(PLOT_MULTI_IRF, flags, &err);
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

    for (i=0; i<n && !err; i++) {
	int vtarg = gretl_VAR_get_variable_number(var, i);

	for (j=0; j<n; j++) {
	    gretl_matrix *resp;
	    int vshock;

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
	    sprintf(title, "%s -> %s", dset->varname[vshock], 
		    dset->varname[vtarg]);
	    fprintf(fp, "set title '%s'\n", title);

	    fputs("plot \\\n", fp);

	    if (confint && use_fill) {
		fprintf(fp, "'-' using 1:2:3 notitle w filledcurve lt %d, \\\n", 
			SHADECOLOR + 1);
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

    fp = open_plot_input_file(PLOT_REGULAR, &err);
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

    set_lwstr(NULL, 0, lwstr);

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
	fprintf(fp, "'-' using 1:2 title '%s' w lines%s", dset->varname[v], lwstr);
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

    fp = open_plot_input_file(PLOT_REGULAR, &err);
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
    int nvars, nobs, jump;
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

    fp = open_plot_input_file(PLOT_MULTI_SCATTER, &err);
    if (err) {
	return err;
    }

    fprintf(fp, "set multiplot layout %d,1\n", nvars);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set noxlabel\n", fp);
    fputs("set noylabel\n", fp);

    gretl_push_c_numeric_locale();

    startdate = obs[t1];
    jump = nobs / (2 * dset->pd);
    fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);

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
	    if (na(eti)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.10g\n", eti);
	    }
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

    fp = open_plot_input_file(PLOT_ROOTS, &err);
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
    fputs("plot 1 w lines, \\\n'-' w points pt 7\n", fp);

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

    fp = open_plot_input_file(PLOT_ELLIPSE, &err);
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

    fputs("plot x(t), y(t) title '', \\\n", fp);
    fprintf(fp, "%g, y(t) title '' w lines lt 2, \\\n", b[0] - maxerr[0]);
    fprintf(fp, "%g, y(t) title '' w lines lt 2, \\\n", b[0] + maxerr[0]);
    fprintf(fp, "x(t), %g title '' w lines lt 2, \\\n", b[1] - maxerr[1]);
    fprintf(fp, "x(t), %g title '' w lines lt 2\n", b[1] + maxerr[1]);

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
					int m, double pm, 
					gretlopt opt,
					FILE *fp)
{
    char crit_string[16];
    double ymin, ymax;
    int k;

    sprintf(crit_string, "%.2f/T^%.1f", 1.96, 0.5);

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
    fprintf(fp, "plot \\\n"
	    "'-' using 1:2 notitle w impulses lw 5, \\\n"
	    "%g title '+- %s' lt 2, \\\n"
	    "%g notitle lt 2\n", pm, crit_string, -pm);
    for (k=0; k<m; k++) {
	fprintf(fp, "%d %g\n", k + 1, acf[k]);
    }
    fputs("e\n", fp);

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
		"%g notitle lt 2\n", pm, crit_string, -pm);
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
		      int m, double pm, 
		      gretlopt opt)
{
    FILE *fp;
    int err = 0;

    fp = open_plot_input_file(PLOT_CORRELOGRAM, &err);

    if (!err) {
	real_correlogram_print_plot(vname, acf, pacf, 
				    m, pm, opt, fp);
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
			    const double *x,
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
	fprintf(fp, "set yrange [0:%g]\n", 1.2 * gretl_max(0, T/2, x));
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
	fprintf(fp, "%g %g\n", ft, (opt & OPT_L)? log(x[t]) : x[t]);
    }

    gretl_pop_c_numeric_locale();

    fputs("e\n", fp);

    return err;
}

int periodogram_plot (const char *vname,
		      int T, int L, const double *x,
		      gretlopt opt)
{
    FILE *fp;
    int err = 0;

    fp = open_plot_input_file(PLOT_PERIODOGRAM, &err);

    if (!err) {
	real_pergm_plot(vname, T, L, x, opt, fp);
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
	fp = open_plot_input_file(PLOT_QQ, &err);
    }

    if (err) {
	free(x);
	free(y);
	return err;
    }   

    fprintf(fp, "set title \"%s\"\n", _("Q-Q plot"));
    fputs("set datafile missing '?'\n", fp);
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
    char title[48];
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

    fp = open_plot_input_file(PLOT_QQ, &err);
    if (err) {
	free(y);
	return err;
    }

    sprintf(title, _("Q-Q plot for %s"), series_get_graph_name(dset, v));
    fprintf(fp, "set title \"%s\"\n", title);
    fputs("set datafile missing '?'\n", fp);
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

/**
 * xy_plot_with_control:
 * @list: list of variables by ID number: Y, X, control.
 * @literal: extra gnuplot commands or %NULL.
 * @dset: dataset struct.
 * @opt: use %OPT_G for GUI graph.
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

/** 
 * gnuplot_process_file:
 * @opt: may include %OPT_U for output to specified file.
 * @prn: gretl printing struct.
 *
 * Respond to the "gnuplot" command with %OPT_D, to specify
 * that input should be taken from a user-created gnuplot
 * command file.
 *
 * Returns: 0 on success, or if ignored; otherwise error code.
 */

int gnuplot_process_file (gretlopt opt, PRN *prn)
{
    const char *inname = get_optval_string(GNUPLOT, OPT_D);
    FILE *fp, *fq;
    int err = 0;

    if (inname == NULL && *inname == '\0') {
	return E_DATA;
    }

    fp = gretl_fopen(inname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    fq = open_plot_input_file(PLOT_USER, &err);

    if (err) {
	fclose(fp);
    } else {
	char line[1024];

	while (fgets(line, sizeof line, fp)) {
	    fputs(line, fq);
	}

	fclose(fp);
	err = finalize_plot_input_file(fq);
    }

    return err;
}

/* The complication below: up till version 4.6, gnuplot used
   a non-standard base of the year 2000 for its equivalent of
   time_t, so it's necessary to adjust by the number of seconds
   between the start of 1970 and the start of 2000 when
   calling *nix time functions. With version 4.7 (CVS) gnuplot
   switched to a base of 1970.

   The 64-bit Windows package of gretl includes gnuplot 4.7 so
   if WIN64 is defined we know we don't need the adjustment.
*/

#ifndef WIN64
# define GP_TIME_OFFSET 946684800
#endif

void date_from_gnuplot_time (char *targ, size_t tsize, 
			     const char *fmt, double x)
{
#if defined(WIN64)
    time_t etime = (time_t) x;

    strftime(targ, tsize, fmt, localtime(&etime));
#elif defined(WIN32)
    time_t etime;

    if (gnuplot_version() < 4.7) {
	x += GP_TIME_OFFSET;
    }

    etime = (time_t) x;
    strftime(targ, tsize, fmt, localtime(&etime));
#else
    struct tm t = {0};
    time_t etime;

    if (gnuplot_version() < 4.7) {
	x += GP_TIME_OFFSET;
    }   

    etime = (time_t) x;
    localtime_r(&etime, &t);
    strftime(targ, tsize, fmt, &t);
#endif
}

double gnuplot_time_from_date (const char *s, const char *fmt)
{
    double x = NADBL;

    if (fmt != NULL && *fmt != '\0') {
	struct tm t = {0};
	time_t etime;
	char *test;
	
	test = strptime(s, fmt, &t);
	if (test != NULL && *test == '\0') {
	    /* conversion went OK */
	    etime = mktime(&t);
	    x = (double) etime;
#ifndef WIN64
	    if (gnuplot_version() < 4.7) {
		x -= GP_TIME_OFFSET;
	    } 
#endif
	}
    }

    return x;
}
