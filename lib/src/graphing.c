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
#include "plotspec.h"

#include <unistd.h>

#define GP_DEBUG 0

#ifdef WIN32
# include <windows.h>
#else
# include <glib.h>
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

static char gnuplot_path[MAXLEN];
static int gp_small_font_size;

typedef struct gnuplot_info_ gnuplot_info;

struct gnuplot_info_ {
    GptFlags flags;
    FitType fit;
    int *list;
    int t1;
    int t2;
    double xrange;
    FILE *fp;
    const char *yformula;
    const double *x;
    double *yvar1;
    double *yvar2;
};

#define MAX_LETTERBOX_LINES 8
#define PREFER_CAIRO 1

#define ts_plot(g) ((g)->flags & GPT_TS)
#define use_impulses(g) ((g)->flags & GPT_IMPULSES)

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
    { PLOT_VAR_ROOTS,      "VAR inverse roots plot" },
    { PLOT_ELLIPSE,        "confidence ellipse plot" },
    { PLOT_MULTI_IRF,      "multiple impulse responses" },
    { PLOT_PANEL,          "multiple panel plots" },
    { PLOT_BI_GRAPH,       "double time-series plot" },
    { PLOT_MANY_TS,        "multiple timeseries" },
    { PLOT_TYPE_MAX,       NULL }
};

static void graph_list_adjust_sample (int *list, 
				      gnuplot_info *ginfo,
				      const double **Z);
static void clear_gpinfo (gnuplot_info *gi);
    
#ifndef WIN32

#define SPAWN_DEBUG 0  /* define != 0 for debugging statements */

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
    char errbuf[32];
    int child_pid = 0, sinp = 0, serr = 0;
    GError *error = NULL;
    gchar *argv[] = {
	NULL,
	NULL
    };

    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    argv[0] = gnuplot_path;

    signal(SIGCHLD, SIG_DFL);

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
	int test, status;
	int errbytes;

	write(sinp, cmd, strlen(cmd));
	write(sinp, "\n", 1);
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
		ret = 1;
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

#endif /* !WIN32 */

static GptFlags get_gp_flags (gretlopt opt, int k, FitType *f)
{
    GptFlags flags = 0;

    if (opt & OPT_B) {
	flags |= GPT_BATCH;
    }

    if (opt & OPT_G) {
	flags |= GPT_GUI;
    } 

    if (opt & OPT_R) {
	flags |= GPT_RESIDS;
    } else if (opt & OPT_F) {
	flags |= GPT_FA;
    }

    if (opt & OPT_M) {
	flags |= GPT_IMPULSES;
    } else if (opt & OPT_O) {
	flags |= GPT_LINES;
    }

    if (opt & OPT_Z) {
	flags |= GPT_DUMMY;
    } else {
	if (opt & OPT_S) {
	    flags |= GPT_FIT_OMIT;
	}
	if (opt & OPT_T) {
	    flags |= GPT_IDX;
	}
    }

    if (k == 2 && !(opt & OPT_S)) {
	/* OPT_S suppresses auto-fit */
	if (opt & OPT_I) {
	    *f = PLOT_FIT_INVERSE;
	} else if (opt & OPT_Q) {
	    *f = PLOT_FIT_QUADRATIC;
	} else if (opt & OPT_L) {
	    *f = PLOT_FIT_LOESS;
	} else if (opt & OPT_N) {
	    *f = PLOT_FIT_OLS;
	}
    }

#if GP_DEBUG
    print_gnuplot_flags(flags, 0);
#endif

    return flags;
}

static void printvars (FILE *fp, int t, const int *list, const double **Z,
		       const double *x, const char *label, double offset)
{
    double xt;
    int i;

    if (x == NULL && list[0] == 2) {
	/* skip missing obs for simple scatterplot */
	if (na(Z[list[1]][t]) || na(Z[list[2]][t])) {
	    return;
	}
    }

    if (x != NULL) {
	xt = x[t] + offset;
	fprintf(fp, "%.8g ", xt);
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
	    fprintf(fp, "%.8g ", xt);
	}
    }

    if (label != NULL) {
	fprintf(fp, "# %s", label);
    }

    fputc('\n', fp);
}

static int factorized_vars (gnuplot_info *gi, const double **Z)
{
    int ynum = gi->list[1];
    int dum = gi->list[3];
    int fn, t, i;

    fn = gi->t2 - gi->t1 + 1;

    gi->yvar1 = malloc(fn * sizeof *gi->yvar1);
    if (gi->yvar1 == NULL) {
	return 1;
    }

    gi->yvar2 = malloc(fn * sizeof *gi->yvar2);
    if (gi->yvar2 == NULL) {
	free(gi->yvar1);
	return 1;
    }

    i = 0;
    for (t=gi->t1; t<=gi->t2; t++) {
	if (na(Z[ynum][t])) {
	    gi->yvar1[i] = NADBL;
	    gi->yvar2[i] = NADBL;
	} else {
	    if (Z[dum][t] == 1.0) {
		gi->yvar1[i] = Z[ynum][t];
		gi->yvar2[i] = NADBL;
	    } else {
		gi->yvar1[i] = NADBL;
		gi->yvar2[i] = Z[ynum][t];
	    }
	}
	i++;
    }

    return 0;
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
#if PREFER_CAIRO
    return GP_PDF_CAIRO;
#else
    return GP_PDF_PDFLIB;
#endif
}

int gnuplot_png_terminal (void)
{
#if PREFER_CAIRO
    return GP_PNG_CAIRO;
#else
    return GP_PNG_GD2;
#endif
}
   
int gnuplot_has_style_fill (void)
{
    /* ... and that it does style fill */
    return 1;
}

static int gnuplot_uses_datafile_missing (void)
{
    /* yup */
    return 1;
}

const char *gnuplot_label_front_string (void)
{
    /* ... and that it handles "front" for labels */
    return " front";
}

int gnuplot_has_latin5 (void)
{
    /* ... and that it supports ISO-8859-9 */
    return 1;
}

int gnuplot_has_cp1250 (void)
{
    /* ... and that it supports CP1250 */
    return 1;
}

int gnuplot_has_cp1254 (void)
{
    /* ... and that it doesn't support CP1254 */
    return 0;
}

int gnuplot_has_rgb (void)
{
    /* ... and that it supports rgb line-color specs */
    return 1;
}

int gnuplot_has_bbox (void)
{
    /* ... and that it supports bounding box info */
    return 1;
}

int gnuplot_has_utf8 (void)
{
    /* ... and that it supports "set encoding utf8" */
    return 1;
}

#else

int gnuplot_has_ttf (int reset)
{
    static int err = -1; 

    /* try a range of ttf fonts that might plausibly be installed
       with X11 */

    if (err == -1 || reset) {
	err = gnuplot_test_command("set term png font luxisr 8");
	if (err) {
	    err = gnuplot_test_command("set term png font Vera 8");
	}
	if (err) {
	    err = gnuplot_test_command("set term png font arial 8");
	}
    }

    return !err;
}

static int gnuplot_has_size (void)
{
    static int err = -1; 
    
    if (err == -1) {
	err = gnuplot_test_command("set term png size 640,480");
    }

    return !err;
}

int gnuplot_has_latin5 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set encoding iso_8859_9");
    }

    return !err;
}

int gnuplot_has_cp1250 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set encoding cp1250");
    }

    return !err;
}

int gnuplot_has_cp1254 (void)
{
#if 1
    int err = 1;
#else /* not yet */
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set encoding cp1254");
    }
#endif
    return !err;
}

int gnuplot_pdf_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
#if PREFER_CAIRO
	int err = gnuplot_test_command("set term pdfcairo");
#else
	int err = 1;
#endif

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

static int gnuplot_has_x11 (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term x11");
    }

    return !err;
}

/* We should enable pngcairo as the default as soon as possible -- but
   for the present this poses a problem with regard to guessing the
   pixel bounds in the PNG file.  So we'll accept pngcairo only if
   we find a version that supports the TERM_XMIN printable
   variable.
*/

int gnuplot_png_terminal (void)
{
    static int ret = -1;

    if (ret == -1) {
#if PREFER_CAIRO
	int err = gnuplot_test_command("set term pngcairo");
#else
	int err = 1;
#endif

	if (!err) {
	    fprintf(stderr, "gnuplot: using pngcairo driver\n");
	    ret = GP_PNG_CAIRO;
	} else {
	    /* try the old-style command: if it fails, we have 
	       the libgd driver, we hope! */
	    err = gnuplot_test_command("set term png color");
	    if (!err) {
		fprintf(stderr, "gnuplot: got old png driver\n");
		ret = GP_PNG_OLD;
	    } else {
		fprintf(stderr, "gnuplot: using libgd png driver\n");
		err = gnuplot_test_command("set term png truecolor");
		ret = (err)? GP_PNG_GD1 : GP_PNG_GD2;
	    }
	}
    }

    return ret;
}

int gnuplot_has_bbox (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set term png ; "
				   "set output '/dev/null' ; "
				   "plot x ; print GPVAL_TERM_XMIN");
    }

    return !err;    
}

int gnuplot_has_utf8 (void)
{
    static int err = -1;

    if (err == -1) {
	err = gnuplot_test_command("set encoding utf8");
    }

    return !err;    
}

int gnuplot_has_style_fill (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style fill solid");
    }

    return !err;
}

static int gnuplot_uses_datafile_missing (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set datafile missing \"?\"");
    }

    return !err;
}

int gnuplot_has_rgb (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style line 2 lc rgb \"#0000ff\"");
    }

    return !err;
}

const char *gnuplot_label_front_string (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set label 'foo' at 0,0 front");
    }

    if (err) {
	return "";
    } else {
	return " front";
    }
}

#endif /* !WIN32 */

static int gnuplot_png_use_aa = 1;

void gnuplot_png_set_use_aa (int s)
{
    gnuplot_png_use_aa = s;
}

/* apparatus for handling plot colors */

enum {
    OLD_PNG_COLOR,  /* plain old "color" option */
    GD_PNG_COLOR,   /* pre-4.2.0 libgd-based PNG color spec */
    RGB_LINE_COLOR  /* per-line rgb settings */
};

static const gretlRGB default_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 }, /* full-intensity green is not very legible */
    { 0x9b, 0xa6, 0xbb }  /* color for box fill */
};

static gretlRGB user_color[N_GP_COLORS] = {
    { 0xff, 0x00, 0x00 },
    { 0x00, 0x00, 0xff },
    { 0x00, 0xcc, 0x00 },
    { 0x9b, 0xa6, 0xbb }
};

static void print_rgb_x (char *s, gretlRGB color)
{
    sprintf(s, "x%02x%02x%02x", color.r, color.g, color.b);
}

void print_rgb_hash (char *s, const gretlRGB *color)
{
    sprintf(s, "#%02x%02x%02x", color->r, color->g, color->b);
}

void print_palette_string (char *s)
{
    sprintf(s, "x%02x%02x%02x x%02x%02x%02x x%02x%02x%02x x%02x%02x%02x",
	    user_color[0].r, user_color[0].g, user_color[0].b,
	    user_color[1].r, user_color[1].g, user_color[1].b,
	    user_color[2].r, user_color[2].g, user_color[2].b,
	    user_color[3].r, user_color[3].g, user_color[3].b);
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
    if (i == BOXCOLOR) {
	user_color[BOXCOLOR] = default_color[BOXCOLOR];
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    user_color[i] = default_color[i];
	}
    }
}

static void 
write_gnuplot_font_string (char *fstr, const char *grfont, PlotType ptype,
			   int pngterm)
{
    if (pngterm == GP_PNG_CAIRO) {
	char fname[128];
	int fsize, nf;

	nf = sscanf(grfont, "%s %d", fname, &fsize);
	if (nf == 2) {
	    if ((ptype == PLOT_MULTI_IRF || ptype == PLOT_MULTI_SCATTER) 
		&& gp_small_font_size > 0) {
		sprintf(fstr, " font \"%s,%d\"", fname, gp_small_font_size);
	    } else {
		sprintf(fstr, " font \"%s,%d\"", fname, fsize);
	    }
	} else if (nf == 1) {
	    sprintf(fstr, " font \"%s\"", fname);
	}
    } else {
	int shrink = 0;

	if ((ptype == PLOT_MULTI_IRF || ptype == PLOT_MULTI_SCATTER) 
	    && gp_small_font_size > 0) {
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
    if (ptype == PLOT_MULTI_IRF || ptype == PLOT_MULTI_SCATTER) {
	strcpy(fstr, " tiny");
    } else {
	strcpy(fstr, " small");
    }
}

#endif

/* we need this only if we don't have per-line rgb
   settings, which are in gnuplot 4.2 and higher */

static char *make_png_colorspec (char *targ, int ptype)
{
    char cstr[8];
    int i;

    /* background etc. */
    strcpy(targ, " xffffff x000000 x202020");

    if (frequency_plot_code(ptype)) {
	strcat(targ, " ");
	print_rgb_x(cstr, user_color[BOXCOLOR]);
	strcat(targ, cstr);
	strcat(targ, " x000000");
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    strcat(targ, " ");
	    print_rgb_x(cstr, user_color[i]);
	    strcat(targ, cstr);
	}
    }

    return targ;
}

/* we use this mechanism with gnuplot 4.2 and higher */

void write_plot_line_styles (int ptype, FILE *fp)
{
    char cstr[8];
    int i;
    
    if (frequency_plot_code(ptype)) {
	print_rgb_hash(cstr, &user_color[BOXCOLOR]);
	fprintf(fp, "set style line 1 lc rgb \"%s\"\n", cstr);
	fputs("set style line 2 lc rgb \"#000000\"\n", fp);
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    print_rgb_hash(cstr, &user_color[i]);
	    fprintf(fp, "set style line %d lc rgb \"%s\"\n", i+1, cstr);
	}
    }

    fputs("set style increment user\n", fp);
}

/* end colors apparatus */

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds, if gnuplot supports this.
*/

void print_plot_bounding_box_request (FILE *fp)
{
    fprintf(fp, "set print '%sgretltmp.png.bounds'\n", gretl_dot_dir());
    fputs("print \"pixel_bounds: \", GPVAL_TERM_XMIN, GPVAL_TERM_XMAX, "
	  "GPVAL_TERM_YMIN, GPVAL_TERM_YMAX\n", fp);
    fputs("print \"data_bounds: \", GPVAL_X_MIN, GPVAL_X_MAX, "
	  "GPVAL_Y_MIN, GPVAL_Y_MAX\n", fp);
}

static void do_plot_bounding_box (void)
{
    FILE *fp = fopen(gretl_plotfile(), "a");

    if (fp != NULL) {
	print_plot_bounding_box_request(fp);
	fclose(fp);
    }
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
    static char png_term_line[256];
    char truecolor_string[12] = {0};
    char font_string[128];
    char size_string[16];
    char color_string[64];
    int gpcolors, gpttf = 1, gpsize = 1;
    int pngterm = 0;
    const char *grfont = NULL;

    *font_string = 0;
    *size_string = 0;
    *color_string = 0;

    pngterm = gnuplot_png_terminal();

#ifdef WIN32
    gpcolors = RGB_LINE_COLOR;
#else
    if (gnuplot_has_rgb()) {
	gpcolors = RGB_LINE_COLOR;
    } else if (pngterm == GP_PNG_OLD) {
	gpcolors = OLD_PNG_COLOR;
    } else {
	gpcolors = GD_PNG_COLOR;
    }

    gpttf = gnuplot_has_ttf(0);
    gpsize = gnuplot_has_size();
#endif

    if (pngterm == GP_PNG_GD2 && gnuplot_png_use_aa) {
	strcpy(truecolor_string, " truecolor");
    }    

    /* plot font setup */
    if (gpttf) {
	grfont = gretl_png_font();
	if (*grfont == 0) {
	    grfont = getenv("GRETL_PNG_GRAPH_FONT");
	}
	if (grfont != NULL && *grfont != 0) {
	    write_gnuplot_font_string(font_string, grfont, ptype, pngterm);
	}
    } 

#ifndef WIN32
    if (!gpttf) {
	write_old_gnuplot_font_string(font_string, ptype);
    }
#endif

    /* plot color setup */
    if (gpcolors == GD_PNG_COLOR) {
	make_png_colorspec(color_string, ptype);
    } else if (gpcolors == OLD_PNG_COLOR) {
	strcpy(color_string, " color");
    } else {
	/* handled via styles */
	*color_string = '\0';
    }

    if (gpsize) {
	if (flags & GPT_LETTERBOX) {
	    strcpy(size_string, " size 680,400");
	} else if (ptype == PLOT_VAR_ROOTS) {
	    strcpy(size_string, " size 480,480");
	}
    }

    if (pngterm == GP_PNG_CAIRO) {
	sprintf(png_term_line, "set term pngcairo%s%s",
		font_string, size_string);
	strcat(png_term_line, "\nset encoding utf8"); /* FIXME? */
    } else {
	sprintf(png_term_line, "set term png%s%s%s%s",
		truecolor_string, font_string, size_string, 
		color_string);
    }

#if GP_DEBUG
    fprintf(stderr, "png term line:\n'%s'\n", png_term_line);
#endif

    return png_term_line;
}

static void png_font_to_emf (const char *pngfont, char *emfline)
{
    char name[128];
    int pt;

    if (sscanf(pngfont, "%127s %d", name, &pt) == 2) {
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
    static char emf_term_line[256];
    const char *grfont = NULL;
    
    strcpy(emf_term_line, "set term emf ");

    if (color) {
	strcat(emf_term_line, "color ");
    } else {
	strcat(emf_term_line, "mono dash ");
    }

    /* font spec */
    grfont = gretl_png_font();
    if (grfont != NULL && *grfont != 0) {
	png_font_to_emf(grfont, emf_term_line);
    }

    return emf_term_line;
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

int write_plot_type_string (PlotType ptype, FILE *fp)
{
    int i, ret = 0;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (ptype == ptinfo[i].ptype) {
	    fprintf(fp, "# %s\n", ptinfo[i].pstr);
	    ret = 1;
	    break;
	}
    }

    return ret;
}

static int real_gnuplot_init (PlotType ptype, int flags, FILE **fpp)
{
    int gui = gretl_in_gui_mode();
    char plotfile[FILENAME_MAX] = {0};

    if (gretl_looping()) {
	return E_OK;
    }

    /* 'gnuplot_path' is file-scope static var */
    if (*gnuplot_path == 0) {
	strcpy(gnuplot_path, gretl_gnuplot_path());
    }

    if (gui) {
	sprintf(plotfile, "%sgpttmp.XXXXXX", gretl_dot_dir());
	if (mktemp(plotfile) == NULL) {
	    return E_FOPEN;
	}
    } else {
	sprintf(plotfile, "%sgpttmp.plt", gretl_dot_dir());
    }

    set_gretl_plotfile(plotfile);

    *fpp = gretl_fopen(plotfile, "w");
    if (*fpp == NULL) {
	fprintf(stderr, "gnuplot_init: couldn't write to %s\n", plotfile);
	return E_FOPEN;
    } 

    if (gui) {
	fprintf(*fpp, "%s\n", get_gretl_png_term_line(ptype, flags));
	fprintf(*fpp, "set output '%sgretltmp.png'\n", gretl_dot_dir());
    }

    write_plot_type_string(ptype, *fpp);

    if (gnuplot_has_rgb()) {
	write_plot_line_styles(ptype, *fpp);
    }

#if GP_DEBUG
    fprintf(stderr, "gnuplot_init: set plotfile = '%s'\n", 
	    plotfile);
#endif

    return 0;
}

/**
 * gnuplot_init:
 * @ptype: indication of the sort of plot to be made.
 * @fpp: pointer to stream to be opened.
 *
 * If we're in GUI mode: writes a unique temporary filename into
 * the internal variable #gretl_plotfile; opens plotfile for writing 
 * as @fpp; and writes initial lines into the output file to select 
 * the PNG terminal type, and direct gnuplot's ouput to a temporary
 * file in the gretl user directory.  
 *
 * If not in GUI mode, opens as @fpp the file %gpttmp.plt in the
 * gretl user directory.  
 *
 * This function is not used in batch mode.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gnuplot_init (PlotType ptype, FILE **fpp)
{
    return real_gnuplot_init(ptype, 0, fpp);
}

/**
 * gnuplot_make_graph:
 *
 * Executes gnuplot, passing as an argument the gretl plotfile.
 *
 * Returns: the return value from the system command.
 */

int gnuplot_make_graph (void)
{
    char plotcmd[MAXLEN];
    int err = 0;

    if (gretl_in_gui_mode() && gnuplot_has_bbox()) {
	do_plot_bounding_box();
    }

#ifdef WIN32
    sprintf(plotcmd, "\"%s\" \"%s\"", gretl_gnuplot_path(), gretl_plotfile());
    err = winfork(plotcmd, NULL, SW_SHOWMINIMIZED, 0);
#else
    sprintf(plotcmd, "%s%s \"%s\"", gretl_gnuplot_path(), 
	    (gretl_in_gui_mode())? "" : " -persist", gretl_plotfile());
    err = gretl_spawn(plotcmd);  
#endif

#if GP_DEBUG
    fprintf(stderr, "gnuplot_make_graph:\n"
	    " plotcmd='%s', err = %d\n", plotcmd, err);
#endif

    return err;
}

enum {
    GTITLE_VLS,
    GTITLE_RESID,
    GTITLE_AF,
    GTITLE_AFV
} graph_titles;

static void make_gtitle (gnuplot_info *gi, int code, 
			 const char *s1, const char *s2)
{
    char depvar[VNAMELEN];
    char title[128];

    switch (code) {
    case GTITLE_VLS:
	if (gi->fit == PLOT_FIT_OLS) {
	    sprintf(title, G_("%s versus %s (with least squares fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    sprintf(title, G_("%s versus %s (with inverse fit)"),
		    s1, s2);
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    sprintf(title, _("%s versus %s (with quadratic fit)"),
		    s1, s2);
	}	    
	break;
    case GTITLE_RESID:
	if (sscanf(s1, "residual for %15s", depvar) == 1) {
	    sprintf(title, G_("Regression residuals (= observed - fitted %s)"), 
		    depvar);
	}
	break;
    case GTITLE_AF:
	sprintf(title, G_("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || (gi->flags & GPT_TS)) {
	    sprintf(title, G_("Actual and fitted %s"), s1);
	} else {
	    sprintf(title, G_("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	*title = '\0';
	break;
    }

    if (*title != '\0') {
	fprintf(gi->fp, "set title \"%s\"\n", title);
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
	*p = 0;
	strncat(p, s, len);
	fprintf(fp, "%s\n", front_strip(p));
	free(p);
    }
}

static void print_gnuplot_literal_lines (const char *s, FILE *fp)
{
    const char *p;

    p = s = front_strip(s);

    while (*s && *s != '}') {
	if (*s == ';') {
	    line_out(p, s - p, fp);
	    p = s + 1;
	}
	s++;
    }
}

static int gretl_plot_count;

void reset_plot_count (void)
{
    gretl_plot_count = 0;
}

static int
get_gnuplot_output_file (FILE **fpp, GptFlags flags, int code)
{
    const char *plotfile = gretl_plotfile();
    int err = 0;

    *fpp = NULL;

    if ((flags & GPT_FILE) && *plotfile != '\0') {
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else if (flags & GPT_BATCH) {
	char fname[FILENAME_MAX];

	if (*plotfile == '\0' || strstr(plotfile, "gpttmp") != NULL) {
	    sprintf(fname, "%sgpttmp%02d.plt", gretl_work_dir(), 
		    ++gretl_plot_count);
	    set_gretl_plotfile(fname);
	} 
	plotfile = gretl_plotfile();
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else {
	/* note: gnuplot_init not used in batch mode */
	err = real_gnuplot_init(code, flags, fpp);
    }

    return err;
}

static int 
loess_plot (gnuplot_info *gi, const double **Z, const DATAINFO *pdinfo)
{
    gretl_matrix *y = NULL;
    gretl_matrix *x = NULL;
    gretl_matrix *yh = NULL;
    int yno = gi->list[1];
    int xno = gi->list[2];
    const char *s1, *s2;
    FILE *fp = NULL;
    char title[96];
    int t, T, d = 1;
    double q = 0.5;
    int err;

    graph_list_adjust_sample(gi->list, gi, Z);
    if (gi->t1 == gi->t2 || gi->list[0] != 2) {
	return GRAPH_NO_DATA;
    }

    if (get_gnuplot_output_file(&fp, gi->flags, PLOT_REGULAR)) {
	return E_FOPEN;
    } 

    err = gretl_plotfit_matrices(yno, xno, PLOT_FIT_LOESS, Z,
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

    s1 = var_get_graph_name(pdinfo, gi->list[1]);
    s2 = var_get_graph_name(pdinfo, gi->list[2]);

    sprintf(title, G_("%s versus %s (with loess fit)"), s1, s2);
    fputs("set key top left\n", fp);
    fprintf(fp, "set title \"%s\"\n", title);
    print_axis_label('x', s1, fp);
    print_axis_label('y', s2, fp);
    print_auto_fit_string(PLOT_FIT_LOESS, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 title '' w points, \\\n", fp);
    sprintf(title, G_("loess fit, d = %d, q = %g"), d, q);
    fprintf(fp, " '-' using 1:2 title \"%s\" w lines\n", title);

    T = gretl_vector_get_length(yh);

    gretl_push_c_numeric_locale();

    for (t=0; t<T; t++) {
	fprintf(fp, "%.8g %.8g\n", x->val[t], y->val[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<T; t++) {
	fprintf(fp, "%.8g %.8g\n", x->val[t], yh->val[t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (!(gi->flags & GPT_BATCH)) {
	err = gnuplot_make_graph();
    }    

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(x);
    gretl_matrix_free(yh);
    clear_gpinfo(gi);

    return err;
} 

static int get_fitted_line (gnuplot_info *gi, 
			    const double **Z, const DATAINFO *pdinfo, 
			    char *targ)
{
    gretl_matrix *y = NULL;
    gretl_matrix *X = NULL;
    gretl_matrix *b = NULL;
    gretl_matrix *V = NULL;
    int yno = gi->list[1];
    int xno = gi->list[2];
    double s2, *ps2 = NULL;
    FitType f = gi->fit;
    char title[72];
    int k, err;

    if (gi->fit == PLOT_FIT_NONE) {
	f = PLOT_FIT_OLS;
	ps2 = &s2;
    }

    k = (f == PLOT_FIT_QUADRATIC)? 3 : 2;

    err = gretl_plotfit_matrices(yno, xno, f, Z,
				 pdinfo->t1, pdinfo->t2,
				 &y, &X);

    if (!err) {
	b = gretl_column_vector_alloc(k);
	if (b == NULL) {
	    err = E_ALLOC;
	}
    }

    if (f == PLOT_FIT_OLS) {
	V = gretl_matrix_alloc(2, 2);
	if (V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_ols(y, X, b, V, NULL, ps2);
    }

    if (!err) {
	double *c = b->val;

	if (gi->fit == PLOT_FIT_NONE) {
	    double v = gretl_matrix_get(V, 1, 1);
	    int T = gretl_vector_get_length(y);
	    double pv = student_pvalue_2(c[1] / sqrt(v), T - k);

	    if (pv < .10) {
		sprintf(title, "Y = %#.3g %c %#.3gX", b->val[0],
			(c[1] > 0)? '+' : '-', fabs(c[1]));
		gretl_push_c_numeric_locale();
		sprintf(targ, "%g + %g*x title '%s' w lines\n", 
			c[0], c[1], title);
		gretl_pop_c_numeric_locale();
		gi->fit = PLOT_FIT_OLS;
	    }
	} else if (gi->fit == PLOT_FIT_OLS) {
	    sprintf(title, "Y = %#.3g %c %#.3gX", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]));
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%g + %g*x title '%s' w lines\n", 
		    c[0], c[1], title);
	    gretl_pop_c_numeric_locale();
	} else if (gi->fit == PLOT_FIT_INVERSE) {
	    sprintf(title, "Y = %#.3g %c %#.3g/X", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]));
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%g + %g/x title '%s' w lines\n", 
		    c[0], c[1], title);
	    gretl_pop_c_numeric_locale();
	} else if (gi->fit == PLOT_FIT_QUADRATIC) {
	    sprintf(title, "Y = %#.3g %c %#.3gX %c %#.3gX^2", c[0],
		    (c[1] > 0)? '+' : '-', fabs(c[1]),
		    (c[2] > 0)? '+' : '-', fabs(c[2])),
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%g + %g*x + %g*x**2 title '%s' w lines\n", 
		    c[0], c[1], c[2], title);
	    gretl_pop_c_numeric_locale();
	}

	if (gi->fit != PLOT_FIT_NONE) {
	    gi->flags |= GPT_AUTO_FIT;
	}
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

static void 
print_x_range (gnuplot_info *gi, const double *x)
{
    if (gretl_isdummy(gi->t1, gi->t2, x)) {
	fputs("set xrange [-1:2]\n", gi->fp);	
	fputs("set xtics (\"0\" 0, \"1\" 1)\n", gi->fp);
	gi->xrange = 3;
    } else {
	double xmin0, xmin, xmax;

	gretl_minmax(gi->t1, gi->t2, x, &xmin0, &xmax);
	gi->xrange = xmax - xmin0;
	xmin = xmin0 - gi->xrange * .025;
	if (xmin0 >= 0.0 && xmin < 0.0) {
	    xmin = 0.0;
	}
	xmax += gi->xrange * .025;
	fprintf(gi->fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);
	gi->xrange = xmax - xmin;
    }
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
				const double **Z, 
				const DATAINFO *pdinfo)
{
    double xx, yy;
    int i, s, t;

    for (i=0; i<2; i++) {
	for (t=gi->t1; t<=gi->t2; t++) {
	    s = t - gi->t1;
	    if (gi->x != NULL) {
		xx = gi->x[t];
	    } else {
		xx = Z[gi->list[2]][t];
		if (na(xx)) {
		    continue;
		}
	    }
	    yy = (i > 0)? gi->yvar2[s] : gi->yvar1[s];
	    if (na(yy)) {
		fprintf(gi->fp, "%.8g ?\n", xx);
	    } else {
		fprintf(gi->fp, "%.8g %.8g", xx, yy);
		if (!(gi->flags & GPT_TS)) {
		    if (pdinfo->markers) {
			fprintf(gi->fp, " # %s", pdinfo->S[t]);
		    } else if (dataset_is_time_series(pdinfo)) {
			char obs[OBSLEN];

			ntodate(obs, t, pdinfo);
			fprintf(gi->fp, " # %s", obs);
		    }
		}
		fputc('\n', gi->fp);
	    }
	}
	fputs("e\n", gi->fp);
    }

    return 0;
}

static void 
maybe_print_panel_jot (int t, const DATAINFO *pdinfo, FILE *fp)
{
    char obs[OBSLEN];
    int maj, min;

    ntodate(obs, t, pdinfo);
    sscanf(obs, "%d:%d", &maj, &min);
    if (maj > 1 && min == 1) {
	fprintf(fp, "%g ?\n", t + 0.5);
    }
}

static void
print_gp_data (gnuplot_info *gi, const double **Z, 
	       const DATAINFO *pdinfo)
{
    int n = gi->t2 - gi->t1 + 1;
    double offset = 0.0;
    int datlist[3];
    int ynum = 2;
    int i, t;

    /* multi impulse plot? calculate offset for lines */
    if (use_impulses(gi) && gi->list[0] > 2) { 
	offset = 0.10 * gi->xrange / n;
    }

    if (gi->x != NULL) {
	datlist[0] = 1;
	ynum = 1;
    } else {
	datlist[0] = 2;
	datlist[1] = gi->list[gi->list[0]];
    }

    /* loop across the variables, printing x then y[i] for each i */

    for (i=1; i<gi->list[0]; i++) {
	double xoff = offset * (i - 1);

	datlist[ynum] = gi->list[i];

	for (t=gi->t1; t<=gi->t2; t++) {
	    const char *label = NULL;
	    char obs[OBSLEN];

	    if (!(gi->flags & GPT_TS) && i == 1) {
		if (pdinfo->markers) {
		    label = pdinfo->S[t];
		} else if (dataset_is_time_series(pdinfo)) {
		    ntodate(obs, t, pdinfo);
		    label = obs;
		}
	    }

	    if ((gi->flags & GPT_TS) && pdinfo->structure == STACKED_TIME_SERIES) {
		maybe_print_panel_jot(t, pdinfo, gi->fp);
	    }

	    printvars(gi->fp, t, datlist, Z, gi->x, label, xoff);
	}

	fputs("e\n", gi->fp);
    }
}

static int
gpinfo_init (gnuplot_info *gi, gretlopt opt, const int *list, 
	     const char *literal, int t1, int t2)
{
    int l0 = list[0];

    gi->fit = PLOT_FIT_NONE;

    gi->flags = get_gp_flags(opt, l0, &gi->fit);
    gi->flags |= GPT_TS; /* may be renounced later */

    gi->t1 = t1;
    gi->t2 = t2;
    gi->xrange = 0.0;
    gi->yformula = NULL;
    gi->fp = NULL;

    gi->x = NULL;
    gi->yvar1 = NULL;
    gi->yvar2 = NULL;
    gi->list = NULL;

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
	&& !(gi->flags & GPT_DUMMY)) {
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
    free(gi->yvar1);
    free(gi->yvar2);
    free(gi->list);

    if (gi->fp != NULL) {
	fclose(gi->fp);
    }
}

#if GP_DEBUG
static void print_gnuplot_flags (int flags, int revised)
{
    if (revised) {
	fprintf(stderr, "*** gnuplot flags after initial revision:\n");
    } else {
	fprintf(stderr, "*** gnuplot() called with flags:\n");
    }

    if (flags & GPT_IMPULSES) {
	fprintf(stderr, " GPT_IMPULSES\n");
    }
    if (flags & GPT_LINES) {
	fprintf(stderr, " GPT_LINES\n");
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
    if (flags & GPT_BATCH) {
	fprintf(stderr, " GPT_BATCH\n");
    }
    if (flags & GPT_GUI) {
	fprintf(stderr, " GPT_GUI\n");
    }
    if (flags & GPT_FIT_OMIT) {
	fprintf(stderr, " GPT_FIT_OMIT\n");
    }
    if (flags & GPT_DATA_STYLE) {
	fprintf(stderr, " GPT_DATA_STYLE\n");
    }
    if (flags & GPT_FILE) {
	fprintf(stderr, " GPT_FILE\n");
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

static void set_lwstr (const DATAINFO *pdinfo, int v, char *s)
{
    int w = var_get_linewidth(pdinfo, v);

    if (w > 1) {
	sprintf(s, " lw %d", w);
    } else {
	*s = '\0';
    }
}

static void set_withstr (GptFlags flags, char *str)
{
    if (flags & GPT_DATA_STYLE) {
	*str = 0;
    } else if (flags & GPT_LINES) {
	strcpy(str, "w lines");
    } else {
	strcpy(str, "w points");
    }
}

static void graph_list_adjust_sample (int *list, 
				      gnuplot_info *ginfo,
				      const double **Z)
{
    int t1min = ginfo->t1;
    int t2max = ginfo->t2;
    int t_ok;
    int i, t, vi;

    for (t=t1min; t<=t2max; t++) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi > 0 && !na(Z[vi][t])) {
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
	    if (vi > 0 && !na(Z[vi][t])) {
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
		if (!na(Z[vi][t])) {
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
}

static int maybe_add_plotx (gnuplot_info *gi, 
			    const DATAINFO *pdinfo)
{
    int k = gi->list[0];
    int add0 = 0;

    /* are we really doing a time-series plot? */
    if (k > 1 && !strcmp(pdinfo->varname[gi->list[k]], "time")) {
	; /* yes */
    } else if (gi->flags & GPT_IDX) {
	add0 = 1; /* yes */
    } else {
	/* no: get out */
	gi->flags &= ~GPT_TS;
	return 0;
    }

    gi->x = gretl_plotx(pdinfo);
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
    if (gnuplot_uses_datafile_missing()) {
	fputs("set datafile missing \"?\"\n", fp);
    } else {
	fputs("set missing \"?\"\n", fp);
    }
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
    mname[4] = '\0';
}

/* for short daily time-series plots: write month names
   into the xtics */

static void make_named_month_tics (gnuplot_info *gi, double yrs, PRN *prn)
{
    double t0 = gi->x[gi->t1];
    double t1 = gi->x[gi->t2];
    double x, tw = 1.0/12;
    int i, m, n = 0;
    char mname[8];
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

    pputs(prn, "# literal lines = 1\n"); 
    pputs(prn, "set xtics ("); 
    x = t0;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
	if (m == 1) {
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%4.0f\" %.8g", x, x);
	    notfirst = 1;
	} else if ((scale == 1) || (m % scale == 1)) {
	    graph_month_name(mname, m);
	    if (notfirst) {
		pputs(prn, ", ");
	    }
	    pprintf(prn, "\"%s\" %.8g", mname, x);
	    notfirst = 1;
	}
	m++;
	x += tw;
	if (m > 12) m -= 12;
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

/**
 * gnuplot:
 * @plotlist: list of variables to plot, by ID number.
 * @literal: commands to be passed to gnuplot.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: option flags.
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gnuplot (const int *plotlist, const char *literal,
	     const double **Z, const DATAINFO *pdinfo, 
	     gretlopt opt)
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
    int oddman = 0;
    int many = 0;
    int i, err = 0;

    gnuplot_info gi;

    gretl_error_clear();

#if GP_DEBUG
    printlist(plotlist, "gnuplot: plotlist");
#endif

    /* below: did have "height 1 width 1 box" for win32,
       "width 1 box" otherwise */
    strcpy(keystr, "set key left top\n");

    err = gpinfo_init(&gi, opt, plotlist, literal, 
		      pdinfo->t1, pdinfo->t2);
    if (err) {
	goto bailout;
    }

    if (gi.fit == PLOT_FIT_LOESS) {
	return loess_plot(&gi, Z, pdinfo);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
	many = 1;
    }

    err = maybe_add_plotx(&gi, pdinfo);
    if (err) {
	goto bailout;
    }

    /* convenience pointer */
    list = gi.list;

    if (gi.flags & GPT_IMPULSES) {
	strcpy(withstr, "w i");
    }

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GPT_TS)) {
	int v = (gi.flags & GPT_DUMMY)? list[2] : list[list[0]];

	strcpy(xlabel, var_get_graph_name(pdinfo, v));
    }

    prn = gretl_print_new(GRETL_PRINT_BUFFER);
    if (prn == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    /* adjust sample range, and reject if it's empty */
    graph_list_adjust_sample(list, &gi, Z);
    if (gi.t1 == gi.t2 || list[0] < 2) {
	err = GRAPH_NO_DATA;
	goto bailout;
    }

    /* add a simple regression line if appropriate */
    if (!use_impulses(&gi) && !(gi.flags & GPT_FIT_OMIT) && list[0] == 2 && 
	!(gi.flags & GPT_TS) && !(gi.flags & GPT_RESIDS)) {
	get_fitted_line(&gi, Z, pdinfo, fit_line);
	pprintf(prn, "# X = '%s' (%d)\n", pdinfo->varname[list[2]], list[2]);
	pprintf(prn, "# Y = '%s' (%d)\n", pdinfo->varname[list[1]], list[1]);
    }

    /* separation by dummy: create special vars */
    if (gi.flags & GPT_DUMMY) { 
	if (list[0] != 3 || factorized_vars(&gi, Z)) {
	    err = E_DATA;
	    goto bailout;
	}
    } 

    /* special tics for time series plots */
    if (gi.flags & GPT_TS) {
	if (many) {
	    pprintf(prn, "# multiple timeseries %d\n", pdinfo->pd);
	} else {
	    int gpsize = 1;

#ifndef WIN32
	    gpsize = gnuplot_has_size();
#endif
	    pprintf(prn, "# timeseries %d", pdinfo->pd);
	    if (gpsize) {
		gi.flags |= GPT_LETTERBOX;
		pputs(prn, " (letterbox)\n");
	    } else {
		pputc(prn, '\n');
	    }
	} 
	if (pdinfo->pd == 4 && (gi.t2 - gi.t1) / 4 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 4\n");
	} else if (pdinfo->pd == 12 && (gi.t2 - gi.t1) / 12 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 12\n");
	} else if (dated_daily_data(pdinfo)) {
	    double yrs = (gi.t2 - gi.t1 + 1.0) / (pdinfo->pd * 52.0);

	    if (yrs <= 3) {
		make_named_month_tics(&gi, yrs, prn);
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
    } 

    /* open file and dump the prn into it: we delaying writing
       the file header till we know a bit more about the plot
    */
    if (get_gnuplot_output_file(&fp, gi.flags, PLOT_REGULAR)) {
	err = E_FOPEN;
	gretl_print_destroy(prn);
	goto bailout;
    } 

    gi.fp = fp;
    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    print_axis_label('x', xlabel, fp);
    fputs("set xzeroaxis\n", fp); 
    gnuplot_missval_string(fp);

    if (list[0] == 2) {
	/* only two variables */
	if (gi.flags & GPT_AUTO_FIT) {
	    print_auto_fit_string(gi.fit, fp);
	    if (gi.flags & GPT_FA) {
		make_gtitle(&gi, GTITLE_AFV, var_get_graph_name(pdinfo, list[1]), 
			    var_get_graph_name(pdinfo, list[2]));
	    } else {
		make_gtitle(&gi, GTITLE_VLS, var_get_graph_name(pdinfo, list[1]), 
			    xlabel);
	    }
	}
	if (gi.flags & GPT_RESIDS && !(gi.flags & GPT_DUMMY)) { 
	    make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	    fprintf(fp, "set ylabel '%s'\n", G_("residual"));
	} else {
	    print_axis_label('y', var_get_graph_name(pdinfo, list[1]), fp);
	}
	if (!(gi.flags & GPT_AUTO_FIT)) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((gi.flags & GPT_RESIDS) && (gi.flags & GPT_DUMMY)) { 
	make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	fprintf(fp, "set ylabel '%s'\n", G_("residual"));
    } else if (gi.flags & GPT_FA) {
	if (list[3] == pdinfo->v - 1) { 
	    /* x var is just time or index: is this always right? */
	    make_gtitle(&gi, GTITLE_AF, var_get_graph_name(pdinfo, list[2]), NULL);
	} else {
	    make_gtitle(&gi, GTITLE_AFV, var_get_graph_name(pdinfo, list[2]), 
			var_get_graph_name(pdinfo, list[3]));
	}
	print_axis_label('y', var_get_graph_name(pdinfo, list[2]), fp);
    } 

    if (many) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    if (gi.x != NULL) {
	print_x_range(&gi, gi.x);
    } else {
	int k = list[0];
	int v = (gi.flags & GPT_DUMMY)? list[k - 1] : list[k];

	print_x_range(&gi, Z[v]);
    }

    if (gi.flags & GPT_Y2AXIS) { 
	check_for_yscale(&gi, Z, &oddman);
	if (gi.flags & GPT_Y2AXIS) {
	    fputs("set ytics nomirror\n", fp);
	    fputs("set y2tics\n", fp);
	}
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
	for (i=1; i<list[0]; i++) {
	    set_lwstr(pdinfo, list[i], lwstr);
	    fprintf(fp, "'-' using 1:2 axes %s title \"%s (%s)\" %s%s%s",
		    (i == oddman)? "x1y2" : "x1y1",
		    var_get_graph_name(pdinfo, list[i]), 
		    (i == oddman)? G_("right") : G_("left"),
		    (use_impulses(&gi))? "w impulses" : 
		    (gi.flags & GPT_TS)? "w lines" : "w points",
		    lwstr,
		    (i == list[0] - 1)? "\n" : ", \\\n");
	}
    } else if (gi.flags & GPT_DUMMY) { 
	strcpy(s1, (gi.flags & GPT_RESIDS)? G_("residual") : 
	       var_get_graph_name(pdinfo, list[1]));
	strcpy(s2, var_get_graph_name(pdinfo, list[3]));
	fprintf(fp, " '-' using 1:2 title \"%s (%s=1)\", \\\n", s1, s2);
	fprintf(fp, " '-' using 1:2 title \"%s (%s=0)\"\n", s1, s2);
    } else if (gi.yformula != NULL) {
	fprintf(fp, " '-' using 1:2 title \"%s\" w points , \\\n", G_("actual"));	
	fprintf(fp, "%s title '%s' w lines\n", gi.yformula, G_("fitted"));
    } else if (gi.flags & GPT_FA) {
	set_withstr(gi.flags, withstr);
	fprintf(fp, " '-' using 1:2 title \"%s\" %s lt 2, \\\n", G_("fitted"), withstr);
	fprintf(fp, " '-' using 1:2 title \"%s\" %s lt 1\n", G_("actual"), withstr);	
    } else {
	for (i=1; i<list[0]; i++)  {
	    set_lwstr(pdinfo, list[i], lwstr);
	    if (list[0] == 2) {
		*s1 = '\0';
	    } else {
		strcpy(s1, var_get_graph_name(pdinfo, list[i]));
	    }
	    if (!use_impulses(&gi)) { 
		set_withstr(gi.flags, withstr);
	    }
	    fprintf(fp, " '-' using 1:2 title \"%s\" %s%s", s1, withstr, lwstr);
	    if (i < list[0] - 1 || (gi.flags & GPT_AUTO_FIT)) {
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
	print_gp_dummy_data(&gi, Z, pdinfo);
    } else {
	print_gp_data(&gi, Z, pdinfo);
    }

    /* flush stream */
    fclose(gi.fp);
    gi.fp = NULL;

    gretl_pop_c_numeric_locale();

    if (!(gi.flags & GPT_BATCH)) {
	err = gnuplot_make_graph();
    }

 bailout:

    clear_gpinfo(&gi);

    return err;
}

/**
 * multi_scatters:
 * @list: list of variables to plot, by ID number.
 * @Z: data array.
 * @pdinfo: data information struct.
 * @opt: can include %OPT_L to use lines.
 *
 * Writes a gnuplot plot file to display up to 6 small graphs
 * based on the variables in @list, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const int *list, const double **Z, 
		    const DATAINFO *pdinfo, gretlopt opt)
{
    GptFlags flags = 0;
    int xvar = 0, yvar = 0;
    const double *obs = NULL;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

    if (opt & OPT_B) {
	flags = GPT_BATCH;
    }

    if (opt & OPT_L) {
	flags |= GPT_LINES;
    }

    pos = gretl_list_separator_position(list);

    if (pos == 0) {
	/* plot against time or index */
	obs = gretl_plotx(pdinfo);
	if (obs == NULL) {
	    return E_ALLOC;
	}
	plotlist = gretl_list_copy(list);
	flags |= GPT_LINES;
    } else if (pos > 2) { 
	/* plot several yvars against one xvar */
	plotlist = gretl_list_new(pos - 1);
	xvar = list[list[0]];
    } else {       
	/* plot one yvar against several xvars */
	plotlist = gretl_list_new(list[0] - pos);
	yvar = list[1];
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

    /* max 6 plots */
    if (plotlist[0] > 6) {
	plotlist[0] = 6;
    }

    nplots = plotlist[0];
    gp_small_font_size = (nplots > 4)? 6 : 0;

    if (get_gnuplot_output_file(&fp, flags, PLOT_MULTI_SCATTER)) {
	return E_FOPEN;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (obs != NULL) {
	double startdate = obs[pdinfo->t1];
	int jump = (pdinfo->t2 - pdinfo->t1 + 1) / (2 * pdinfo->pd);

	fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);
    } else {
	fputs("set noxtics\nset noytics\n", fp);
    }

    for (i=0; i<nplots; i++) {  
	int pv = plotlist[i+1];

	if (nplots <= 4) {
	    fputs("set size 0.45,0.5\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.5,0.5\n", fp);
	    else if (i == 2) fputs("0.0,0.0\n", fp);
	    else if (i == 3) fputs("0.5,0.0\n", fp);
	} else {
	    fputs("set size 0.31,0.45\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.32,0.5\n", fp);
	    else if (i == 2) fputs("0.64,0.5\n", fp);
	    else if (i == 3) fputs("0.0,0.0\n", fp);
	    else if (i == 4) fputs("0.32,0.0\n", fp);
	    else if (i == 5) fputs("0.64,0.0\n", fp);
	}

	if (obs != NULL) {
	    fputs("set noxlabel\n", fp);
	    fputs("set noylabel\n", fp);
	    fprintf(fp, "set title '%s'\n", pdinfo->varname[pv]);
	} else {
	    fprintf(fp, "set xlabel '%s'\n",
		    (yvar)? pdinfo->varname[pv] :
		    pdinfo->varname[xvar]);
	    fprintf(fp, "set ylabel '%s'\n", 
		    (yvar)? pdinfo->varname[yvar] :
		    pdinfo->varname[pv]);
	}

	fputs("plot '-' using 1:2", fp);
	if (flags & GPT_LINES) {
	    fputs(" with lines", fp);
	}
	fputc('\n', fp);

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    double xx;

	    xx = (yvar)? Z[pv][t] : (xvar)? Z[xvar][t] : obs[t];

	    if (na(xx)) {
		fputs("? ", fp);
	    } else {
		fprintf(fp, "%.8g ", xx);
	    }

	    xx = (yvar)? Z[yvar][t] : Z[pv][t];

	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.8g\n", xx);
	    }
	}

	fputs("e\n", fp);
    } 

    gretl_pop_c_numeric_locale();

    fputs("set nomultiplot\n", fp);

    fclose(fp);

    if (!(flags & GPT_BATCH)) {
	err = gnuplot_make_graph();
    }

    free(plotlist);

    return err;
}

static int get_3d_output_file (FILE **fpp)
{
    char fname[MAXLEN];
    int err = 0;

    sprintf(fname, "%sgpttmp.plt", gretl_dot_dir());
    *fpp = gretl_fopen(fname, "w");

    if (*fpp == NULL) {
	err = E_FOPEN;
    } else {
	set_gretl_plotfile(fname);
    }

    return err;
}

static void 
maybe_add_surface (const int *list, double ***pZ, DATAINFO *pdinfo, 
		   gretlopt opt, char *surface)
{
    MODEL smod;
    double umin, umax, vmin, vmax;
    int olslist[5];

    olslist[0] = 4;
    olslist[1] = list[3];
    olslist[2] = 0;
    olslist[3] = list[2];
    olslist[4] = list[1];

    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[2]], &umin, &umax);
    gretl_minmax(pdinfo->t1, pdinfo->t2, (*pZ)[list[1]], &vmin, &vmax);

    smod = lsq(olslist, pZ, pdinfo, OLS, OPT_A);

    if (!smod.errcode && !na(smod.fstt) &&
	(snedecor_cdf_comp(smod.fstt, smod.dfn, smod.dfd) < .10 || (opt & OPT_F))) {
	double uadj = (umax - umin) * 0.02;
	double vadj = (vmax - vmin) * 0.02;

	sprintf(surface, "[u=%g:%g] [v=%g:%g] "
		"%g+(%g)*u+(%g)*v title '', ", 
		umin - uadj, umax + uadj, 
		vmin - vadj, vmax + vadj,
		smod.coeff[0], smod.coeff[1],
		smod.coeff[2]);
    } 

    clear_model(&smod);
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: unused at present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo,  
		gretlopt opt)
{
    FILE *fq = NULL;
    int t, t1 = pdinfo->t1, t2 = pdinfo->t2;
    int orig_t1 = pdinfo->t1, orig_t2 = pdinfo->t2;
    int lo = list[0];
    int datlist[4];
    char surface[128] = {0};

    if (lo != 3) {
	fprintf(stderr, "gnuplot_3d needs three variables (only)\n");
	return -1;
    }

    if (get_3d_output_file(&fq)) {
	return E_FOPEN;
    }

    varlist_adjust_sample(list, &t1, &t2, (const double **) *pZ);

    /* if resulting sample range is empty, complain */
    if (t2 == t1) {
	fclose(fq);
	return GRAPH_NO_DATA;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

#ifndef WIN32
    if (gnuplot_has_x11()) {
	/* wxt is too slow/jerky? */
	fputs("set term x11\n", fq);
    }
#endif

    gretl_push_c_numeric_locale();

    maybe_add_surface(list, pZ, pdinfo, opt, surface);
    
    print_axis_label('x', var_get_graph_name(pdinfo, list[2]), fq);
    print_axis_label('y', var_get_graph_name(pdinfo, list[1]), fq);
    print_axis_label('z', var_get_graph_name(pdinfo, list[3]), fq);

    gnuplot_missval_string(fq);

    if (literal != NULL && *literal != 0) {
	print_gnuplot_literal_lines(literal, fq);
    }

#ifdef WIN32
    fprintf(fq, "splot %s'-' title '' w p lt 3\n", surface);
#else
    fprintf(fq, "splot %s'-' title ''\n", surface);
#endif

    datlist[0] = 3;
    datlist[1] = list[2];
    datlist[2] = list[1];
    datlist[3] = list[3];

    for (t=t1; t<=t2; t++) {
	const char *label = NULL;

	if (pdinfo->markers) {
	    label = pdinfo->S[t];
	}
	printvars(fq, t, datlist, (const double **) *pZ, NULL, label, 0.0);
    }	
    fputs("e\n", fq);

    gretl_pop_c_numeric_locale();

    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    fclose(fq);

    return 0;
}

static void print_freq_test_label (char *s, const char *format, 
				   double v, double pv)
{
    gretl_pop_c_numeric_locale();
    sprintf(s, format, v, pv);
    gretl_push_c_numeric_locale();
}

static void print_freq_dist_label (char *s, int dist, double x, double y)
{
    int dcomma = 0;
#ifdef ENABLE_NLS
    char test[8];

    gretl_pop_c_numeric_locale();
    sprintf(test, "%g", 0.5);
    if (strchr(test, ',')) {
	dcomma = 1;
    }
#endif

    if (dist == D_NORMAL) {
	sprintf(s, "N(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    } else if (dist == D_GAMMA) {
	sprintf(s, "gamma(%.5g%c%.5g)", x, 
		((dcomma)? ' ' : ','), y);
    }

#ifdef ENABLE_NLS
    gretl_push_c_numeric_locale();
#endif
}

/* Below: a fix for the case where the y-range is by default
   degenerate, in which case gnuplot produces a graph OK, but
   issues a warning and returns non-zero.
*/

static void maybe_set_yrange (FreqDist *freq, double lambda, FILE *fp)
{
    double ymin = 1.0e+20;
    double ymax = -1.0e+20;
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
	fprintf(fp, "set yrange [%g:%g]\n", ymax * lambda * 0.99, 
		ymax * lambda * 1.01);
    } else {
	fprintf(fp, "set yrange [0.0:%g]\n", ymax * lambda * 1.1);
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
    int err;

    if (K == 0) {
	return E_DATA;
    }

    if (K == 1) {
	sprintf(gretl_errmsg, _("'%s' is a constant"), freq->varname);
	return E_DATA;
    }

    if (dist == D_NORMAL) {
	plottype = PLOT_FREQ_NORMAL;
    } else if (dist == D_GAMMA) {
	plottype = PLOT_FREQ_GAMMA;
    } else {
	plottype = PLOT_FREQ_SIMPLE;
    }

    if ((err = gnuplot_init(plottype, &fp))) {
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
		fprintf(fp, "set label \"%s:\" at graph .03, graph .97%s\n",
			G_("Test statistic for normality"),
			gnuplot_label_front_string());
		print_freq_test_label(label, G_("Chi-squared(2) = %.3f pvalue = %.5f"), 
				      freq->test, chisq_cdf_comp(freq->test, 2));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, gnuplot_label_front_string());
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
		fprintf(fp, "set label '%s:' at graph .03, graph .97%s\n",
			G_("Test statistic for gamma"),
			gnuplot_label_front_string());
		print_freq_test_label(label, G_("z = %.3f pvalue = %.5f"), 
				      freq->test, normal_pvalue_2(freq->test));
		fprintf(fp, "set label '%s' at graph .03, graph .93%s\n", 
			label, gnuplot_label_front_string());
	    }	
	}

	/* adjust min, max if needed */
	if (freq->midpt[0] < plotmin) {
	    plotmin = freq->midpt[0];
	}
	if (freq->midpt[K-1] > plotmax) {
	    plotmax = freq->midpt[K-1];
	}

	fprintf(fp, "set xrange [%.7g:%.7g]\n", plotmin, plotmax);
	fputs("set key right top\n", fp);
    } else { 
	/* plain frequency plot (no theoretical distribution shown) */
	lambda = 1.0 / freq->n;
	plotmin = freq->midpt[0] - barwidth;
	plotmax = freq->midpt[K-1] + barwidth;
	fprintf(fp, "set xrange [%.7g:%.7g]\n", plotmin, plotmax);
	maybe_set_yrange(freq, lambda, fp);
	fputs("set nokey\n", fp);
    }

    fprintf(fp, "set xlabel '%s'\n", freq->varname);
    if (dist) {
	fprintf(fp, "set ylabel '%s'\n", G_("Density"));
    } else {
	fprintf(fp, "set ylabel '%s'\n", G_("Relative frequency"));
    }

    if (isnan(lambda)) {
	if (fp != NULL) {
	    fclose(fp);
	}
	return 1;
    }

    /* plot instructions */
    if (use_boxes) {
	if (gnuplot_has_style_fill()) {
	    fputs("set style fill solid 0.6\n", fp);
	}
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
	fprintf(fp, "%.8g %.8g\n", freq->midpt[i], lambda * freq->f[i]);
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    if (fp != NULL) {
	fclose(fp);
    }

    return gnuplot_make_graph();
}

int plot_fcast_errs (int t1, int t2, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     int time_series)
{
    FILE *fp = NULL;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int do_errs = (maxerr != NULL);
    int t, n, err;

    /* don't graph empty portion of forecast */
    for (t=t2; t>=t1; t--) {
	if (na(depvar[t]) && na(yhat[t])) {
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

    if ((err = gnuplot_init(PLOT_FORECAST, &fp))) {
	return err;
    }    

    /* check that we have any values for the actual var */
    for (t=t1; t<=t2; t++) {
	if (!na(depvar[t])) {
	    depvar_present = 1;
	    break;
	}
    }

    fputs("# forecasts with 95 pc conf. interval\n", fp);

    gretl_minmax(t1, t2, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);

    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    if (time_series) {
	fprintf(fp, "# timeseries %d\n", time_series);
    } else if (n < 33) {
	fputs("set xtics 1\n", fp);
    }

    fputs("set key left top\nplot \\\n", fp);
    if (depvar_present) {
	fprintf(fp, "'-' using 1:2 title '%s' w lines , \\\n",
		varname);
    }
    fprintf(fp, "'-' using 1:2 title '%s' w lines", G_("forecast"));
    if (do_errs) {
	fprintf(fp, " , \\\n'-' using 1:2:3 title '%s' w errorbars\n",
		G_("95 percent confidence interval"));
    } else {
	fputc('\n', fp);
    }

    gretl_push_c_numeric_locale();

    if (depvar_present) {
	for (t=t1; t<=t2; t++) {
	    if (na(depvar[t])) {
		fprintf(fp, "%.8g ?\n", obs[t]);
	    } else {
		fprintf(fp, "%.8g %.8g\n", obs[t], depvar[t]);
	    }
	}
	fputs("e\n", fp);
    }

    for (t=t1; t<=t2; t++) {
	if (na(yhat[t])) {
	    fprintf(fp, "%.8g ?\n", obs[t]);
	} else {
	    fprintf(fp, "%.8g %.8g\n", obs[t], yhat[t]);
	}
    }
    fputs("e\n", fp);

    if (do_errs) {
	for (t=t1; t<=t2; t++) {
	    if (na(yhat[t]) || na(maxerr[t])) {
		fprintf(fp, "%.8g ? ?\n", obs[t]);
	    } else {
		fprintf(fp, "%.8g %.8g %.8g\n", obs[t], yhat[t], maxerr[t]);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int garch_resid_plot (const MODEL *pmod, const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    const double *obs;
    const double *h;
    double sd2;
    int t, err;

    h = gretl_model_get_data(pmod, "garch_h");
    if (h == NULL) {
	return E_DATA;
    }

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if ((err = gnuplot_init(PLOT_GARCH, &fp))) {
	return err;
    }

    fprintf(fp, "set key left top\n"
	    "plot \\\n'-' using 1:2 title '%s' w lines, \\\n"
	    "'-' using 1:2 title '%s' w lines lt 2, \\\n" 
	    "'-' using 1:2 notitle w lines lt 2\n", 
	    G_("residual"), G_("+- sqrt(h(t))"));

    gretl_push_c_numeric_locale();

    for (t=pmod->t1; t<=pmod->t2; t++) {
	fprintf(fp, "%.8g %.8g\n", obs[t], pmod->uhat[t]);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = -sqrt(h[t]);
	fprintf(fp, "%.8g %.8g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	sd2 = sqrt(h[t]);
	fprintf(fp, "%.8g %.8g\n", obs[t], sd2);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int 
rmplot (const int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int (*range_mean_graph) (int, const double **, const DATAINFO *, PRN *);
    void *handle = NULL;
    int err;

    range_mean_graph = get_plugin_function("range_mean_graph", &handle);
    if (range_mean_graph == NULL) {
        return 1;
    }

    err = range_mean_graph(list[1], Z, pdinfo, prn);

    close_plugin(handle);

    if (!err && !gretl_in_batch_mode() && !gretl_looping()) {
        err = gnuplot_make_graph();
    }

    return err;
}

int 
hurstplot (const int *list, const double **Z, DATAINFO *pdinfo, PRN *prn)
{
    int (*hurst_exponent) (int, const double **, const DATAINFO *, PRN *);
    void *handle = NULL;
    int err;

    hurst_exponent = get_plugin_function("hurst_exponent", &handle);
    if (hurst_exponent == NULL) {
        return 1;
    }

    err = hurst_exponent(list[1], Z, pdinfo, prn);

    close_plugin(handle);

    if (!err && !gretl_in_batch_mode() && !gretl_looping()) {
        err = gnuplot_make_graph();
    } 

    return err;
}

static void get_x_and_y_sizes (int n, int *x, int *y)
{
    if (n == 2) {
	*x = 2;
	*y = 1;
    } else if (n == 3 || n == 4) {
	*x = 2;
	*y = 2;
    } else if (n == 5 || n == 6) {
	*x = 3;
	*y = 2;
    } else if (n > 6 && n < 10) {
	*x = 3;
	*y = 3;
    } else {
	*x = 0;
	*y = 0;
    }
}

/* panel: plot one variable as a time series, with separate plots for
   each cross-sectional unit 
*/

int 
gretl_panel_ts_plot (const int *list, const double **Z, DATAINFO *pdinfo) 
{
    FILE *fp = NULL;
    int i, j, k;
    int t, t0;
    int xnum, ynum;
    float xfrac, yfrac;
    float xorig = 0.0;
    float yorig;
    int err = 0;

    int T = pdinfo->pd;
    int nunits = pdinfo->n / T;

    get_x_and_y_sizes(nunits, &xnum, &ynum);
    if (xnum == 0 || ynum == 0) {
	return E_DATA;
    }

    err = gnuplot_init(PLOT_PANEL, &fp);
    if (err) {
	return err;
    }

    xfrac = 1.0 / xnum;
    yfrac = 1.0 / ynum;

    fputs("set key top left\n", fp);
    fputs("set multiplot\n", fp);
    fprintf(fp, "set xlabel '%s'\n", _("time"));
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    if (yfrac > 1.4 * xfrac) {
	yfrac = 1.4 * xfrac;
    }
    fprintf(fp, "set size %g,%g\n", xfrac, yfrac);

    k = 0;
    t0 = 0;

    for (i=0; i<xnum; i++) {

	yorig = 1.0 - yfrac;

	for (j=0; j<ynum; j++) {
	    int vj = list[1];

	    if (k == nunits) {
		break;
	    }

	    fprintf(fp, "set origin %g,%g\n", xorig, yorig);
	    fprintf(fp, "set title '%s (%d)'\n", pdinfo->varname[vj], k+1);
	    fputs("plot \\\n'-' using 1:2 notitle w lines\n", fp);

	    for (t=0; t<T; t++) {
		fprintf(fp, "%d %.8g\n", t+1, Z[vj][t+t0]);
	    }
	    fputs("e\n", fp);

	    k++;
	    if (k == nunits) {
		break;
	    }

	    t0 += T;
	    yorig -= yfrac;
	}

	if (k == nunits) {
	    break;
	}	

	xorig += xfrac;
    }

    fputs("unset multiplot\n", fp);
    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 const double **Z,
				 const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    int confint = 0;
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err;

    resp = gretl_VAR_get_impulse_response(var, targ, shock, periods,
					  Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    err = gnuplot_init((confint)? PLOT_IRFBOOT : PLOT_REGULAR, &fp);
    if (err) {
	gretl_matrix_free(resp);
	return err;
    }

    vtarg = gretl_VAR_get_variable_number(var, targ);
    vshock = gretl_VAR_get_variable_number(var, shock);

    if (!confint) {
	fputs("# impulse response plot\n", fp);
    }

    if (confint) {
	fputs("set key top left\n", fp);
	sprintf(title, G_("response of %s to a shock in %s, "
			  "with bootstrap confidence interval"),
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    } else {
	fputs("set nokey\n", fp);
	sprintf(title, G_("response of %s to a shock in %s"), 
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    }

    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", title);

    if (confint) {
	fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines, \\\n", 
		G_("point estimate"));
	fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n",
		G_("0.025 and 0.975 quantiles"));
    } else {
	fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    for (t=0; t<periods; t++) {
	fprintf(fp, "%d %.8g\n", t+1, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint) {
	for (t=0; t<periods; t++) {
	    fprintf(fp, "%d %.8g %.8g %.8g\n", t+1, 
		    gretl_matrix_get(resp, t, 0),
		    gretl_matrix_get(resp, t, 1),
		    gretl_matrix_get(resp, t, 2));
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);
    gretl_matrix_free(resp);

    return gnuplot_make_graph();
}

int 
gretl_VAR_plot_multiple_irf (GRETL_VAR *var, int periods,
			     const double **Z,
			     const DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    int confint = 0;
    int vtarg, vshock;
    gretl_matrix *resp;
    char title[128];
    int t, err, i, j;

    int n = var->neqns;
    float plot_fraction = 1.0 / n;
    float xorig = 0.0;
    float yorig;

    gp_small_font_size = (n == 4)? 6 : 0;

    resp = gretl_VAR_get_impulse_response(var, 1, 1, periods, Z, pdinfo);
    if (resp == NULL) {
	return E_ALLOC;
    }

    if (gretl_matrix_cols(resp) > 1) {
	confint = 1;
    }

    err = gnuplot_init(PLOT_MULTI_IRF, &fp);
    if (err) {
	gretl_matrix_free(resp);
	return err;
    }

    if (!confint) {
	fputs("set nokey\n", fp);
    } else {
	fputs("set key top left\n", fp);
    }
    fputs("set multiplot\n", fp);
    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "set size %g,%g\n", plot_fraction, plot_fraction);

    for (i=0; i<n; i++) {

	yorig = 1.0 - plot_fraction;
	vtarg = gretl_VAR_get_variable_number(var, i);

	for (j=0; j<n; j++) {

	    fprintf(fp, "set origin %g,%g\n", xorig, yorig);
	    resp = gretl_VAR_get_impulse_response(var, i, j, periods, Z, pdinfo);
	    if (resp == NULL) {
		return E_ALLOC;
	    }

	    vshock = gretl_VAR_get_variable_number(var, j);
	    sprintf(title, "%s -> %s", pdinfo->varname[vshock], pdinfo->varname[vtarg]);
	    fprintf(fp, "set title '%s'\n", title);

	    if (confint) {
		fputs("plot \\\n'-' using 1:2 notitle w lines, \\\n", fp); 
		fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
	    } else {
		fputs("plot \\\n'-' using 1:2 w lines\n", fp);
	    }

	    for (t=0; t<periods; t++) {
		fprintf(fp, "%d %.8g\n", t+1, gretl_matrix_get(resp, t, 0));
	    }
	    fputs("e\n", fp);

	    if (confint) {
		for (t=0; t<periods; t++) {
		    fprintf(fp, "%d %.8g %.8g %.8g\n", t+1, 
			    gretl_matrix_get(resp, t, 0),
			    gretl_matrix_get(resp, t, 1),
			    gretl_matrix_get(resp, t, 2));
		}
		fputs("e\n", fp);
	    }
	    
	    yorig -= plot_fraction;
	}

	xorig += plot_fraction;
    }

    fputs("unset multiplot\n", fp);
    gretl_pop_c_numeric_locale();

    fclose(fp);
    gretl_matrix_free(resp);

    return gnuplot_make_graph();
}

int gretl_system_residual_plot (void *p, int ci, const DATAINFO *pdinfo)
{
    GRETL_VAR *var = NULL;
    equation_system *sys = NULL;
    const gretl_matrix *E = NULL;
    FILE *fp = NULL;
    const double *obs;
    int nvars, nobs;
    int i, v, t, t1, err;

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

    t1 = E->t1;

    err = gnuplot_init(PLOT_REGULAR, &fp);
    if (err) {
	return err;
    }

    obs = gretl_plotx(pdinfo);

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);

    fputs("# system residual plot\n", fp);
    fputs("set key top left\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (ci == VAR) {
	fprintf(fp, "set title '%s'\n", G_("VAR residuals"));
    } else {
	fprintf(fp, "set title '%s'\n", G_("System residuals"));
    }

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}
	fprintf(fp, "'-' using 1:2 title '%s' w lines", pdinfo->varname[v]);
	if (i == nvars - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(", \\\n", fp); 
	}
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<nvars; i++) {
	for (t=0; t<nobs; t++) {
	    double eti = gretl_matrix_get(E, t, i);

	    if (obs != NULL) {
		fprintf(fp, "%g %.8g\n", obs[t+t1], eti);
	    } else {
		fprintf(fp, "%d %.8g\n", t+1, eti);
	    }
	}
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return gnuplot_make_graph();
}

int gretl_system_residual_mplot (void *p, int ci, const DATAINFO *pdinfo) 
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

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    nobs = gretl_matrix_rows(E);
    t1 = E->t1;

    err = gnuplot_init(PLOT_MULTI_SCATTER, &fp);
    if (err) {
	return err;
    }

    fputs("set size 1.0,1.0\nset origin 0.0,0.0\n"
	  "set multiplot\n", fp);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    startdate = obs[t1];
    jump = nobs / (2 * pdinfo->pd);
    fprintf(fp, "set xtics %g, %d\n", ceil(startdate), jump);

    gretl_minmax(t1, t1 + nobs - 1, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);	

    for (i=0; i<nvars; i++) { 

	if (nvars <= 4) {
	    fputs("set size 0.45,0.5\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.5,0.5\n", fp);
	    else if (i == 2) fputs("0.0,0.0\n", fp);
	    else if (i == 3) fputs("0.5,0.0\n", fp);
	} else {
	    fputs("set size 0.31,0.45\n", fp);
	    fputs("set origin ", fp);
	    if (i == 0) fputs("0.0,0.5\n", fp);
	    else if (i == 1) fputs("0.32,0.5\n", fp);
	    else if (i == 2) fputs("0.64,0.5\n", fp);
	    else if (i == 3) fputs("0.0,0.0\n", fp);
	    else if (i == 4) fputs("0.32,0.0\n", fp);
	    else if (i == 5) fputs("0.64,0.0\n", fp);
	}

	fputs("set noxlabel\n", fp);
	fputs("set noylabel\n", fp);
	if (var != NULL) {
	    v = gretl_VAR_get_variable_number(var, i);
	} else {
	    v = system_get_depvar(sys, i);
	}
	fprintf(fp, "set title '%s'\n", pdinfo->varname[v]);

	fputs("plot '-' using 1:2 with lines\n", fp);

	for (t=0; t<nobs; t++) {
	    double xx;

	    fprintf(fp, "%.8g\t", obs[t+t1]);
	    xx = gretl_matrix_get(E, t, i);
	    if (na(xx)) {
		fputs("?\n", fp);
	    } else {
		fprintf(fp, "%.8g\n", xx);
	    }
	}

	fputs("e\n", fp);
    } 

    gretl_pop_c_numeric_locale();
    fputs("set nomultiplot\n", fp);
    fclose(fp);

    return gnuplot_make_graph();
}

int gretl_VAR_roots_plot (GRETL_VAR *var)
{
    const gretl_matrix *lam;
    FILE *fp = NULL;
    double x, y;
    double px, py;
    int i, n, err;

    lam = gretl_VAR_get_roots(var);
    if (lam == NULL) {
	return E_ALLOC;
    }

    err = gnuplot_init(PLOT_VAR_ROOTS, &fp);
    if (err) {
	return err;
    }

    n = gretl_matrix_rows(lam);

    fprintf(fp, "set title '%s'\n", 
	    G_("VAR inverse roots in relation to the unit circle"));
    fputs("# literal lines = 8\n", fp);
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 w lines, \\\n"
	  "'-' w points pt 7\n", fp);

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
    fclose(fp);

    return gnuplot_make_graph();
}

int confidence_ellipse_plot (gretl_matrix *V, double *b, double t, double c,
			     const char *iname, const char *jname)
{
    FILE *fp = NULL;
    double maxerr[2];
    double xcoeff[2];
    double ycoeff[2];
    gretl_matrix *e = NULL;
    int err = 0;

    maxerr[0] = t * sqrt(gretl_matrix_get(V, 0, 0));
    maxerr[1] = t * sqrt(gretl_matrix_get(V, 1, 1));

    err = gretl_invert_symmetric_matrix(V);
    if (err) {
	return err;
    }

    e = gretl_symmetric_matrix_eigenvals(V, 1, &err);
    if (err) {
	return err;
    }

    e->val[0] = sqrt(1.0 / e->val[0] * c);
    e->val[1] = sqrt(1.0 / e->val[1] * c);

    xcoeff[0] = e->val[0] * gretl_matrix_get(V, 0, 0);
    xcoeff[1] = e->val[1] * gretl_matrix_get(V, 0, 1);

    ycoeff[0] = e->val[0] * gretl_matrix_get(V, 1, 0);
    ycoeff[1] = e->val[1] * gretl_matrix_get(V, 1, 1);

    gretl_matrix_free(e);

    err = gnuplot_init(PLOT_ELLIPSE, &fp);
    if (err) {
	return err;
    }

    fprintf(fp, "set title '%s'\n",
	    /* xgettext:no-c-format */
	    G_("95% confidence ellipse and 95% marginal intervals"));
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

    fclose(fp);

    return gnuplot_make_graph();
}

int is_auto_fit_string (const char *s)
{
    /* FIXME? */
    if (strstr(s, "automatic fit")) return 1;
    if (strstr(s, G_("with least squares fit"))) return 1;
    return 0;
}

