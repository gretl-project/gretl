/*
 *  Copyright (c) by Allin Cottrell and Riccardo Lucchetti
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

/* graphing.c for gretl */

#include "libgretl.h"
#include "var.h"
#include "libset.h"
#include "matrix_extra.h"

#include <unistd.h>

#define GP_DEBUG 0

#ifdef _WIN32
# include <windows.h>
#else
# ifdef USE_GLIB2
#  define USE_GSPAWN
#  include <glib.h>
#  include <signal.h>
#  if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#  endif
#  ifndef WEXITSTATUS
#   define WEXITSTATUS(stat_val) ((unsigned)(stat_val) >> 8)
#  endif
#  ifndef WIFEXITED
#   define WIFEXITED(stat_val) (((stat_val) & 255) == 0)
#  endif
# endif /* GLIB2 */
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

#define use_impulses(g) ((g)->flags & GPT_IMPULSES)
#define ts_plot(g) ((g)->flags & GPT_TS)

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
    { PLOT_TYPE_MAX,       NULL }
};
    
#ifdef USE_GSPAWN

# define SPAWN_DEBUG 0  /* define != 0 for debugging statements */

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
	(*gnuplot_path == 0)? "gnuplot" : gnuplot_path,
	NULL
    };

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

#elif !defined(_WIN32)

#include <signal.h>

int gnuplot_test_command (const char *cmd)
{
    char fullcmd[512];
    int err;

    signal(SIGCHLD, SIG_DFL);

    sprintf(fullcmd, "echo \"%s\" | %s 2>/dev/null", cmd,
	    (*gnuplot_path == 0)? "gnuplot" : gnuplot_path);
    err =  system(fullcmd);

    return err;
}

#endif /* ! USE_GSPAWN, ! WIN32 */

GptFlags gp_flags (int batch, gretlopt opt)
{
    GptFlags flags = 0;

    if (batch) {
	flags |= GPT_BATCH;
    }

    if (opt & OPT_M) {
	flags |= GPT_IMPULSES;
    } else if (opt & OPT_L) {
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

int gnuplot_has_pdf (void)
{
    /* ... and we know it does PDF output */
    return 1;
}

int gnuplot_has_specified_colors (void)
{
    /* ... and we know it does specified colors */
    return 1;
}

int gnuplot_has_style_fill (void)
{
    /* ... and that it does style fill */
    return 1;
}

static int gnuplot_has_specified_emf_colors (void)
{
    /* ... and we know it does specified emf colors */
    return 1;
}

const char *gnuplot_label_front_string (void)
{
    /* ... and that it handles "front" for labels */
    return " front";
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

int gnuplot_has_pdf (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term pdf");
    }

    return !err;
}

int gnuplot_has_specified_colors (void)
{
    static int err = -1; 

    if (err == -1) {
	/* try the old-style command: 
	   if it fails, we have the new driver, we hope */
	err = gnuplot_test_command("set term png color");
    }

    return err;
}

int gnuplot_has_style_fill (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set style fill solid");
    }

    return !err;
}

static int gnuplot_has_specified_emf_colors (void)
{
    static int err = -1; 

    if (err == -1) {
	err = gnuplot_test_command("set term emf color xff0000");
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

#endif /* WIN32 */

/* apparatus for handling plot colors */

static const char default_palette[N_GP_COLORS][8] = {
    "xff0000", 
    "x0000ff", 
    "x00cc00",
    "x9ba6bb"
};

static char graph_palette[N_GP_COLORS][8] = {
    "xff0000", 
    "x0000ff", 
    "x00cc00",  /* full-intensity green is not very legible */
    "x9ba6bb"   /* color for box fill */
};

void print_palette_string (char *targ)
{
    sprintf(targ, "%s %s %s %s", graph_palette[0],
	    graph_palette[1], graph_palette[2],
	    graph_palette[3]);
}

const char *graph_color_string (int i)
{
#if GP_DEBUG
    fprintf(stderr, "get_graph_palette: i=%d\n", i);
#endif

    return (i >= 0 && i < N_GP_COLORS)? graph_palette[i] : "";
}

static int colstr_is_valid (const char *colstr)
{
    char *test = NULL;
    int cnum;

    if (*colstr != 'x' || strlen(colstr) != 7) {
	return 0;
    }

    cnum = strtol(colstr + 1, &test, 16);

    if (*test != '\0' || cnum < 0 || cnum > 0xffffff) {
	return 0;
    }

    return 1;
}

void set_graph_palette (int i, const char *colstr)
{
    if (i >= 0 && i < N_GP_COLORS && colstr_is_valid(colstr)) {
	*graph_palette[i] = '\0';
	strncat(graph_palette[i], colstr, 7);
    } else {
	fprintf(stderr, "Invalid color spec, '%s'\n", colstr);
    }
}

void graph_palette_reset (int i)
{
    if (i == BOXCOLOR) {
	strcpy(graph_palette[BOXCOLOR], default_palette[BOXCOLOR]);
    } else {
	for (i=0; i<BOXCOLOR; i++) {
	    strcpy(graph_palette[i], default_palette[i]);
	}
    }
}

static void 
write_gnuplot_font_string (char *fstr, const char *grfont, PlotType ptype)
{
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

/* end colors apparatus */

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
    char font_string[128];
    char size_string[16];
    char color_string[64];
    int gpcolors = 1, gpttf = 1;
    const char *grfont = NULL;

    *font_string = 0;
    *size_string = 0;
    *color_string = 0;

#ifndef WIN32
    gpcolors = gnuplot_has_specified_colors();
    gpttf = gnuplot_has_ttf(0);
#endif

    /* plot font setup */
    if (gpttf) {
	grfont = gretl_png_font();
	if (*grfont == 0) {
	    grfont = getenv("GRETL_PNG_GRAPH_FONT");
	}
	if (grfont != NULL && *grfont != 0) {
	    write_gnuplot_font_string(font_string, grfont, ptype);
	}
    }

    /* plot color setup */
    if (gpcolors) {
	int i;

	/* background etc. */
	strcpy(color_string, " xffffff x000000 x202020");

	if (frequency_plot_code(ptype)) {
	    strcat(color_string, " ");
	    strcat(color_string, graph_palette[BOXCOLOR]);
	    strcat(color_string, " x000000");
	} else {
	    for (i=0; i<BOXCOLOR; i++) {
		strcat(color_string, " ");
		strcat(color_string, graph_palette[i]);
	    }
	}
    } else {
	strcpy(color_string, " color"); /* old PNG driver */
    }

    if (flags & GPT_LETTERBOX) {
	strcpy(size_string, " size 680,400");
    } else if (ptype == PLOT_VAR_ROOTS) {
	strcpy(size_string, " size 480,480");
    }

    sprintf(png_term_line, "set term png%s%s%s",
	    font_string, size_string, color_string);

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

    if (color && gnuplot_has_specified_emf_colors()) {
	if (frequency_plot_code(ptype)) {
	    strcat(emf_term_line, graph_palette[BOXCOLOR]);
	    strcat(emf_term_line, " x000000");
	} else {
	    int i;

	    for (i=0; i<BOXCOLOR; i++) {
		strcat(emf_term_line, graph_palette[i]);
		strcat(emf_term_line, " ");
	    }
	}
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
    int i, ret = PLOT_REGULAR;

    for (i=1; i<PLOT_TYPE_MAX; i++) {
	if (!strcmp(str + 2, ptinfo[i].pstr)) {
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
	sprintf(plotfile, "%sgpttmp.XXXXXX", gretl_user_dir());
	if (mktemp(plotfile) == NULL) {
	    return E_FOPEN;
	}
    } else {
	sprintf(plotfile, "%sgpttmp.plt", gretl_user_dir());
    }

    set_gretl_plotfile(plotfile);

    *fpp = gretl_fopen(plotfile, "w");
    if (*fpp == NULL) {
	fprintf(stderr, "gnuplot_init: couldn't write to %s\n", plotfile);
	return E_FOPEN;
    } 

    if (gui) {
	fprintf(*fpp, "%s\n", get_gretl_png_term_line(ptype, flags));
	fprintf(*fpp, "set output '%sgretltmp.png'\n", gretl_user_dir());
    }

    write_plot_type_string(ptype, *fpp);

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

#ifdef ENABLE_NLS

#undef RECODE_DBG

static int recode_gnuplot_file (const char *fname)
{
    FILE *fp, *fq;
    const char *font;
    char oldline[512], newline[1024];
    char rname[FILENAME_MAX];
    int ttf = 0;

    fp = fopen(fname, "r");
    if (fp == NULL) {
	return 1;
    }

    strcpy(rname, fname);
    strcat(rname, "l2");

    fq = fopen(rname, "w");
    if (fq == NULL) {
	fclose(fp);
	return 1;
    }

    font = gretl_png_font();
    if (font != NULL && *font != '\0') {
	ttf = 1;
    }

    while (fgets(oldline, sizeof oldline, fp)) {
	if (isdigit((unsigned char) oldline[0])) {
	    fputs(oldline, fq);
	} else if (ttf) {
	    sprint_l2_to_html(newline, oldline, sizeof newline);
	    fputs(newline, fq);
#if RECODE_DBG
	    fprintf(stderr, "recode (sprint_l2_to_html):\n"
		    " original: '%s'\n modified: '%s'\n",
		    oldline, newline);
#endif
	} else {
	    fputs(oldline, fq);
	}
    }

    fclose(fp);
    fclose(fq);
	    
    remove(fname);
    rename(rname, fname);

    return 0;
}

#endif

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

#ifdef ENABLE_NLS  
    if (use_latin_2() && gnuplot_has_ttf(0)) {
# if GP_DEBUG
	fprintf(stderr, "gnuplot_make_graph: calling recode_gnuplot_file()\n");
# endif
	recode_gnuplot_file(gretl_plotfile());
    } 
#endif

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
    char title[128];
    char depvar[VNAMELEN];

    switch (code) {
    case GTITLE_VLS:
	sprintf(title, I_("%s versus %s (with least squares fit)"),
		s1, s2);
	break;
    case GTITLE_RESID:
	if (sscanf(s1, "residual for %15s", depvar) == 1) {
	    sprintf(title, I_("Regression residuals (= observed - fitted %s)"), 
		    depvar);
	}
	break;
    case GTITLE_AF:
	sprintf(title, I_("Actual and fitted %s"), s1);
	break;
    case GTITLE_AFV:
	if (s2 == NULL || (gi->flags & GPT_TS)) {
	    sprintf(title, I_("Actual and fitted %s"), s1);
	} else {
	    sprintf(title, I_("Actual and fitted %s versus %s"), s1, s2);
	}
	break;
    default:
	*title = '\0';
	break;
    }

    if (*title != '\0') {
	fprintf(gi->fp, "set title '%s'\n", title);
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

static const char *series_name (const DATAINFO *pdinfo, int v)
{
    const char *ret = pdinfo->varname[v];

    if (pdinfo->varinfo != NULL && pdinfo->varinfo[v] != NULL) {
	ret = DISPLAYNAME(pdinfo, v);
	if (ret[0] == '\0') {
	    ret = pdinfo->varname[v];
	}
    } 

    return ret;
}

static int
get_gnuplot_output_file (FILE **fpp, GptFlags flags, 
			 int *plot_count, int code)
{
    const char *plotfile = gretl_plotfile();
    int err = 0;

    *fpp = NULL;

    if ((flags & GPT_FILE) && *plotfile != '\0') {
	*fpp = gretl_fopen(plotfile, "w");
	if (*fpp == NULL) {
	    err = E_FOPEN;
	}
    } else if ((flags & GPT_BATCH) && plot_count != NULL) {  
	char fname[FILENAME_MAX];

	if (*plotfile == '\0' || strstr(plotfile, "gpttmp") != NULL) {
	    *plot_count += 1; 
	    sprintf(fname, "%sgpttmp%02d.plt", gretl_user_dir(), *plot_count);
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

static void print_auto_fit_string (gnuplot_info *gi)
{
    if (gi->fit == PLOT_FIT_OLS) {
	fputs("# plot includes automatic fit: OLS\n", gi->fp);
    } else if (gi->fit == PLOT_FIT_QUADRATIC) {
	fputs("# plot includes automatic fit: quadratic\n", gi->fp);
    } else if (gi->fit == PLOT_FIT_INVERSE) {
	fputs("# plot includes automatic fit: inverse\n", gi->fp);
    } else if (gi->fit == PLOT_FIT_LOESS) {
	fputs("# plot includes automatic fit: loess\n", gi->fp);
    }
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
    double s2, v, pv;
    int T, k = 2;
    int err;

    err = gretl_plotfit_matrices(yno, xno, PLOT_FIT_OLS, Z,
				 pdinfo->t1, pdinfo->t2,
				 &y, &X);

    if (!err) {
	b = gretl_vector_alloc(k);
	V = gretl_matrix_alloc(k, k);
	if (b == NULL || V == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	err = gretl_matrix_ols(y, X, b, V, NULL, &s2);
    }

    if (!err) {
	v = gretl_matrix_get(V, 1, 1);
	T = gretl_vector_get_length(y);
	pv = t_pvalue_2(b->val[1] / sqrt(v), T - k);

	if (pv < .10) {
	    char title[48];

	    sprintf(title, "Y = %#.3g %c %#.3gX", b->val[0],
		    (b->val[1] > 0)? '+' : '-', 
		    fabs(b->val[1]));
	    gretl_push_c_numeric_locale();
	    sprintf(targ, "%g + %g*x title '%s' w lines\n", 
		    b->val[0], b->val[1], title);
	    gretl_pop_c_numeric_locale();
	    gi->flags |= GPT_AUTO_FIT;
	    gi->fit = PLOT_FIT_OLS;
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
	double xmin, xmax;

	gretl_minmax(gi->t1, gi->t2, x, &xmin, &xmax);
	gi->xrange = xmax - xmin;
	xmin -= gi->xrange * .025;
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
		     &ymin[i], &ymax[i]);
    }

    gi->flags &= ~GPT_Y2AXIS;

    for (i=1; i<gi->list[0]; i++) {
	oddcount = 0;
	for (j=1; j<gi->list[0]; j++) {
	    if (j == i) {
		continue;
	    }
	    ratio = ymax[i] / ymax[j];
	    if (ratio > 5.0 || ratio < 0.2) {
		gi->flags |= GPT_Y2AXIS;
		oddcount++;
	    }
	}
	if (oddcount == gi->list[0] - 2) {
	    *oddman = i;
	    break;
	}
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
gp_info_init (gnuplot_info *gi, GptFlags flags, FILE *fp,
	      const int *list, const char *literal, int t1, int t2)
{
    int l0 = list[0];

    gi->flags = flags | GPT_TS;
    gi->t1 = t1;
    gi->t2 = t2;
    gi->xrange = 0.0;
    gi->yformula = NULL;
    gi->fp = fp;

    gi->x = NULL;
    gi->yvar1 = NULL;
    gi->yvar2 = NULL;
    gi->list = NULL;

    if (l0 < 2 && !(flags & GPT_IDX)) {
	return E_ARGS;
    }

    if ((flags & GPT_DUMMY) && (flags & GPT_IDX)) {
	return E_BADOPT;
    }

    gi->list = gretl_list_copy(list);
    if (gi->list == NULL) {
	return E_ALLOC;
    }

    if ((l0 > 2 || (l0 > 1 && (flags & GPT_IDX))) && 
	 l0 < 7 && !(flags & GPT_RESIDS) && !(flags & GPT_FA)
	&& !(flags & GPT_DUMMY)) {
	/* allow probe for using two y axes */
	gi->flags |= GPT_Y2AXIS;
    } 

    if ((flags & GPT_FA) && literal != NULL && 
	!strncmp(literal, "yformula: ", 10)) {
	/* fitted vs actual plot with fitted given by formula */
	gi->yformula = literal + 10;
    }

    if (literal != NULL && strstr(literal, "set style data")) {
	gi->flags |= GPT_DATA_STYLE;
    }

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
static void print_gnuplot_flags (GnuplotFlags flags)
{
    fprintf(stderr, "*** gnuplot() called with flags:\n");

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
}
#endif

static void set_lwstr (const DATAINFO *pdinfo, int v, char *s)
{
    int w = var_get_linewidth(pdinfo, v);

    if (w > 0) {
	sprintf(s, " lw %d", w);
    } else {
	strcpy(s, " lw 1");
    }
}

static void set_withstr (GptFlags flags, const int *lines, 
			 int i, char *str)
{
    int ltest = 0;

    if (flags & GPT_DATA_STYLE) {
	*str = 0;
	return;
    }

    if (lines != NULL) {
	ltest = lines[(flags & GPT_GUI)? i - 1 : 0];
    }

    if (ltest) {
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
    int i, t;

    for (t=t1min; t<=t2max; t++) {
	t_ok = 0;
	for (i=1; i<=list[0]; i++) {
	    if (!na(Z[list[i]][t])) {
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
	    if (!na(Z[list[i]][t])) {
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

	    for (t=t1min; t<=t2max; t++) {
		if (!na(Z[list[i]][t])) {
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
    if (!strcmp(pdinfo->varname[gi->list[k]], "time")) {
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

#define ts_letterbox 1

/**
 * gnuplot:
 * @plotlist: list of variables to plot, by ID number.
 * @lines: vector of 1s and 0s to indicate whether variables should
 * be represented by lines or not (or NULL).
 * @literal: commands to be passed to gnuplot.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @plot_count: pointer to count of graphs drawn so far.
 * @flags: bitwise OR of zero or more options from #GnuplotFlags
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gnuplot (const int *plotlist, const int *lines, const char *literal,
	     const double **Z, const DATAINFO *pdinfo, 
	     int *plot_count, GptFlags flags)
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
    int toomany = 0;
    int i, err = 0;

    gnuplot_info gi;

    gretl_error_clear();

#if GP_DEBUG
    print_gnuplot_flags(flags);
    printlist(plotlist, "gnuplot: plotlist");
#endif

    /* below: did have "height 1 width 1 box" for win32,
       "width 1 box" otherwise */
    strcpy(keystr, "set key left top\n");

    err = gp_info_init(&gi, flags, NULL, plotlist, literal, 
		       pdinfo->t1, pdinfo->t2);
    if (err) {
	goto bailout;
    }

    /* convenience pointer */
    list = gi.list;

    if (list[0] > MAX_PLOT_LINES + 1) {
	toomany = 1;
    }

    err = maybe_add_plotx(&gi, pdinfo);
    if (err) {
	goto bailout;
    }

    if (use_impulses(&gi)) {
	strcpy(withstr, "w i");
    }

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GPT_TS)) {
	if (flags & GPT_DUMMY) {
	    strcpy(xlabel, series_name(pdinfo, list[2])); 
	} else {
	    strcpy(xlabel, series_name(pdinfo, list[list[0]]));
	}
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
    if (!use_impulses(&gi) && !(flags & GPT_FIT_OMIT) && list[0] == 2 && 
	!(gi.flags & GPT_TS) && !(flags & GPT_RESIDS)) {
	get_fitted_line(&gi, Z, pdinfo, fit_line);
	if (gi.fit == PLOT_FIT_OLS) {
	    pprintf(prn, "# X = '%s' (%d)\n", pdinfo->varname[list[2]], list[2]);
	    pprintf(prn, "# Y = '%s' (%d)\n", pdinfo->varname[list[1]], list[1]);
	}
    }

    /* separation by dummy: create special vars */
    if (flags & GPT_DUMMY) { 
	if (list[0] != 3 || factorized_vars(&gi, Z)) {
	    err = E_DATA;
	    goto bailout;
	}
    } 

    /* special tics for short time series plots */
    if (gi.flags & GPT_TS) {
	if (toomany) {
	    pprintf(prn, "# multiple timeseries %d", pdinfo->pd);
	} else {
	    pprintf(prn, "# timeseries %d", pdinfo->pd);
	}
	if (ts_letterbox) {
	    flags |= GPT_LETTERBOX;
	    pputs(prn, " (letterbox)\n");
	} else {
	    pputc(prn, '\n');
	}
	if (pdinfo->pd == 4 && (gi.t2 - gi.t1) / 4 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 4\n");
	}
	if (pdinfo->pd == 12 && (gi.t2 - gi.t1) / 12 < 8) {
	    pputs(prn, "set xtics nomirror 0,1\n"); 
	    pputs(prn, "set mxtics 12\n");
	}
    } else if (toomany) {
	pputs(prn, "# multiple data series\n");
    }

    /* open file and dump the prn into it: we delaying writing
       the file header till we know a bit more about the plot
    */
    if (get_gnuplot_output_file(&fp, flags, plot_count, PLOT_REGULAR)) {
	err = E_FOPEN;
	gretl_print_destroy(prn);
	goto bailout;
    } 

    gi.fp = fp;
    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    fprintf(fp, "set xlabel '%s'\n", xlabel);
    fputs("set xzeroaxis\n", fp); 
    fputs("set missing \"?\"\n", fp);

    if (list[0] == 2) {
	/* only two variables */
	if (gi.flags & GPT_AUTO_FIT) {
	    print_auto_fit_string(&gi);
	    if (flags & GPT_FA) {
		make_gtitle(&gi, GTITLE_AFV, series_name(pdinfo, list[1]), 
			    series_name(pdinfo, list[2]));
	    } else {
		make_gtitle(&gi, GTITLE_VLS, series_name(pdinfo, list[1]), 
			    xlabel);
	    }
	}
	if (flags & GPT_RESIDS && !(flags & GPT_DUMMY)) { 
	    make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	    fprintf(fp, "set ylabel '%s'\n", I_("residual"));
	} else {
	    fprintf(fp, "set ylabel '%s'\n", series_name(pdinfo, list[1]));
	}
	if (!(gi.flags & GPT_AUTO_FIT)) {
	    strcpy(keystr, "set nokey\n");
	}
    } else if ((flags & GPT_RESIDS) && (flags & GPT_DUMMY)) { 
	make_gtitle(&gi, GTITLE_RESID, VARLABEL(pdinfo, list[1]), NULL);
	fprintf(fp, "set ylabel '%s'\n", I_("residual"));
    } else if (flags & GPT_FA) {
	if (list[3] == pdinfo->v - 1) { 
	    /* x var is just time or index: is this always right? */
	    make_gtitle(&gi, GTITLE_AF, series_name(pdinfo, list[2]), NULL);
	} else {
	    make_gtitle(&gi, GTITLE_AFV, series_name(pdinfo, list[2]), 
			series_name(pdinfo, list[3]));
	}
	fprintf(fp, "set ylabel '%s'\n", series_name(pdinfo, list[2]));
    } 

    if (toomany) {
	strcpy(keystr, "set key outside\n");
    }

    fputs(keystr, fp);

    gretl_push_c_numeric_locale();

    if (gi.x != NULL) {
	print_x_range(&gi, gi.x);
    } else {
	int k = list[0];
	int v = (flags & GPT_DUMMY)? list[k - 1] : list[k];

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
	    fprintf(fp, "'-' using 1:2 axes %s title '%s (%s)' %s%s%s",
		    (i == oddman)? "x1y2" : "x1y1",
		    series_name(pdinfo, list[i]), 
		    (i == oddman)? I_("right") : I_("left"),
		    (use_impulses(&gi))? "w impulses" : 
		    (gi.flags & GPT_TS)? "w lines" : "w points",
		    lwstr,
		    (i == list[0] - 1)? "\n" : " , \\\n");
	}
    } else if (flags & GPT_DUMMY) { 
	strcpy(s1, (flags & GPT_RESIDS)? I_("residual") : 
	       series_name(pdinfo, list[1]));
	strcpy(s2, series_name(pdinfo, list[3]));
	fprintf(fp, " '-' using 1:2 title '%s (%s=1)', \\\n", s1, s2);
	fprintf(fp, " '-' using 1:2 title '%s (%s=0)'\n", s1, s2);
    } else if (gi.yformula != NULL) {
	fprintf(fp, " '-' using 1:2 title '%s' w points , \\\n", I_("actual"));	
	fprintf(fp, "%s title '%s' w lines\n", gi.yformula, I_("fitted"));
    } else if (flags & GPT_FA) {
	set_withstr(flags, lines, 1, withstr);
	fprintf(fp, " '-' using 1:2 title '%s' %s lt 2, \\\n", I_("fitted"), withstr);
	set_withstr(flags, lines, 2, withstr);
	fprintf(fp, " '-' using 1:2 title '%s' %s lt 1\n", I_("actual"), withstr);	
    } else {
	for (i=1; i<list[0]; i++)  {
	    set_lwstr(pdinfo, list[i], lwstr);
	    if (list[0] == 2) {
		*s1 = '\0';
	    } else {
		strcpy(s1, series_name(pdinfo, list[i]));
	    }
	    if (!use_impulses(&gi)) { 
		set_withstr(gi.flags, lines, i, withstr);
	    }
	    fprintf(fp, " '-' using 1:2 title '%s' %s%s", s1, withstr, lwstr);
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
    if (flags & GPT_DUMMY) {
	print_gp_dummy_data(&gi, Z, pdinfo);
    } else {
	print_gp_data(&gi, Z, pdinfo);
    }

    /* flush stream */
    fclose(gi.fp);
    gi.fp = NULL;

    gretl_pop_c_numeric_locale();

    if (!(flags & GPT_BATCH)) {
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
 * @plot_count: count of graphs shown to date.
 * @flags: option flags.
 *
 * Writes a gnuplot plot file to display up to 6 small graphs
 * based on the variables in @list, and calls gnuplot to make 
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_scatters (const int *list, const double **Z, 
		    const DATAINFO *pdinfo, int *plot_count, 
		    GptFlags flags)
{
    int xvar = 0, yvar = 0;
    const double *obs = NULL;
    int *plotlist = NULL;
    int pos, nplots = 0;
    FILE *fp = NULL;
    int i, t, err = 0;

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

    if (get_gnuplot_output_file(&fp, flags, plot_count, PLOT_MULTI_SCATTER)) {
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

    sprintf(fname, "%sgpttmp.plt", gretl_user_dir());
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
		   unsigned int flags, char *surface)
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
	(f_cdf_comp(smod.fstt, smod.dfn, smod.dfd) < .10 || flags & GPT_FA)) {
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
 * @plot_count: pointer to counter variable for plots (unused)
 * @flags: unused at present.
 *
 * Writes a gnuplot plot file to display a 3D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo,  
		int *plot_count, GptFlags flags)
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

    gretl_push_c_numeric_locale();

    maybe_add_surface(list, pZ, pdinfo, flags, surface);

    fprintf(fq, "set xlabel '%s'\n", series_name(pdinfo, list[2]));
    fprintf(fq, "set ylabel '%s'\n", series_name(pdinfo, list[1]));
    fprintf(fq, "set zlabel '%s'\n", series_name(pdinfo, list[3]));

    fputs("set missing \"?\"\n", fq);

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
			I_("Test statistic for normality"),
			gnuplot_label_front_string());
		print_freq_test_label(label, I_("Chi-squared(2) = %.3f pvalue = %.5f"), 
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
			I_("Test statistic for gamma"),
			gnuplot_label_front_string());
		print_freq_test_label(label, I_("z = %.3f pvalue = %.5f"), 
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
	fprintf(fp, "set ylabel '%s'\n", I_("Density"));
    } else {
	fprintf(fp, "set ylabel '%s'\n", I_("Relative frequency"));
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
	    fputs("set style fill solid 0.8\n", fp);
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
	fprintf(fp, "'-' using 1:2 title '%s' %s , \\\n"
		"1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) "
		"title '%s' w lines\n",
		freq->varname, withstr, label);
    } else if (dist == D_GAMMA) {
	print_freq_dist_label(label, dist, alpha, beta);
	fputs("plot \\\n", fp);
	fprintf(fp, "'-' using 1:2 title '%s' %s ,\\\n"
		"x**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) "
		"title '%s' w lines\n",
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

    fputs("set missing \"?\"\n", fp);

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
    fprintf(fp, "'-' using 1:2 title '%s' w lines", I_("forecast"));
    if (do_errs) {
	fprintf(fp, " , \\\n'-' using 1:2:3 title '%s' w errorbars\n",
		I_("95 percent confidence interval"));
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
	    "plot \\\n'-' using 1:2 title '%s' w lines , \\\n"
	    "'-' using 1:2 title '%s' w lines lt 2, \\\n" 
	    "'-' using 1:2 notitle w lines lt 2\n", 
	    I_("residual"), I_("+- sqrt(h(t))"));

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
	sprintf(title, I_("response of %s to a shock in %s, "
			  "with bootstrap confidence interval"),
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    } else {
	fputs("set nokey\n", fp);
	sprintf(title, I_("response of %s to a shock in %s"), 
		pdinfo->varname[vtarg], pdinfo->varname[vshock]);
    }

    fprintf(fp, "set xlabel '%s'\n", _("periods"));
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", title);

    if (confint) {
	fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines,\\\n", 
		I_("point estimate"));
	fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n",
		I_("0.025 and 0.975 quantiles"));
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
		fputs("plot \\\n'-' using 1:2 notitle w lines,\\\n", fp); 
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

int gretl_VAR_residual_plot (const GRETL_VAR *var, const DATAINFO *pdinfo)
{
    const gretl_matrix *E;
    FILE *fp = NULL;
    const double *obs;
    int nvars, nobs;
    int i, v, t, t1, err;

    E = gretl_VAR_get_residual_matrix(var);
    if (E == NULL) {
	return E_DATA;
    }

    t1 = gretl_VAR_get_t1(var);

    err = gnuplot_init(PLOT_REGULAR, &fp);
    if (err) {
	return err;
    }

    obs = gretl_plotx(pdinfo);

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);

    fputs("# VAR residual plot\n", fp);
    fputs("set key top left\n", fp);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set title '%s'\n", I_("VAR residuals"));

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
	v = gretl_VAR_get_variable_number(var, i);
	fprintf(fp, "'-' using 1:2 title '%s' w lines", pdinfo->varname[v]);
	if (i == nvars - 1) {
	    fputc('\n', fp);
	} else {
	    fputs(",\\\n", fp); 
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

int gretl_VAR_residual_mplot (const GRETL_VAR *var, const DATAINFO *pdinfo) 
{
    const gretl_matrix *E;
    FILE *fp = NULL;
    const double *obs;
    double startdate;
    double xmin, xmax, xrange;
    int nvars, nobs, jump;
    int i, v, t, t1;
    int err = 0;

    E = gretl_VAR_get_residual_matrix(var);
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
    t1 = gretl_VAR_get_t1(var);

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
	v = gretl_VAR_get_variable_number(var, i);
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
	    I_("VAR inverse roots in relation to the unit circle"));
    fputs("# literal lines = 8\n", fp);
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 w lines , \\\n"
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
    double *e = NULL;
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

    e[0] = sqrt(1.0 / e[0] * c);
    e[1] = sqrt(1.0 / e[1] * c);

    xcoeff[0] = e[0] * gretl_matrix_get(V, 0, 0);
    xcoeff[1] = e[1] * gretl_matrix_get(V, 0, 1);

    ycoeff[0] = e[0] * gretl_matrix_get(V, 1, 0);
    ycoeff[1] = e[1] * gretl_matrix_get(V, 1, 1);

    free(e);

    err = gnuplot_init(PLOT_ELLIPSE, &fp);
    if (err) {
	return err;
    }

    fprintf(fp, "set title '%s'\n",
	    /* xgettext:no-c-format */
	    I_("95% confidence ellipse and 95% marginal intervals"));
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
    if (strstr(s, "automatic fit")) return 1;
    if (strstr(s, I_("with least squares fit"))) return 1;
    return 0;
}

#ifdef ENABLE_NLS

void pprint_gnuplot_encoding (const char *termstr, PRN *prn)
{
    if (strstr(termstr, "postscript")) {
	const char *enc = get_gnuplot_charset();

	if (enc != NULL) {
	    pprintf(prn, "set encoding %s\n", enc);
	}
    }
}

void fprint_gnuplot_encoding (const char *termstr, FILE *fp)
{
    if (strstr(termstr, "postscript")) {
	const char *enc = get_gnuplot_charset();

	if (enc != NULL) {
	    fprintf(fp, "set encoding %s\n", enc);
	}
    }
}

#endif 
