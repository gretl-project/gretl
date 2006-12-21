/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* graphing.h for gretl */

#ifndef GRAPHING_H
#define GRAPHING_H

#include <stdio.h>

#define GRAPH_NO_DATA -999

typedef enum {
    GP_IMPULSES   = 1 << 0,  /* use impulses for plotting */
    GP_LINES      = 1 << 1,  /* force use of lines for plotting */
    GP_RESIDS     = 1 << 2,  /* doing residual plot */
    GP_FA         = 1 << 3,  /* doing fitted/actual plot */
    GP_DUMMY      = 1 << 4,  /* using a dummy for separation */
    GP_BATCH      = 1 << 5,  /* working in batch mode */
    GP_GUI        = 1 << 6,  /* called from GUI context */
    GP_OLS_OMIT   = 1 << 7,  /* Don't draw fitted line on graph */
    GP_DATA_STYLE = 1 << 8,  /* data style is set by user */
    GP_FILE       = 1 << 9,  /* send output to named file */
    GP_IDX        = 1 << 10  /* plot against time or obs index */
} GnuplotFlags;

typedef enum {
    GPTSPEC_TS             = 1 << 0,
    GPTSPEC_Y2AXIS         = 1 << 1,
    GPTSPEC_AUTO_OLS       = 1 << 2,
    GPTSPEC_OLS_HIDDEN     = 1 << 3,
    GPTSPEC_MINIMAL_BORDER = 1 << 4,
    GPTSPEC_PNG_OUTPUT     = 1 << 5,
    GPTSPEC_ALL_MARKERS    = 1 << 6,
    GPTSPEC_ALL_MARKERS_OK = 1 << 7,
    GPTSPEC_NO_BORDER      = 1 << 8
} PlotSpecFlags; 

#define MAXTITLE 128
#define MAX_PLOT_LABELS 3
#define MAX_PLOT_LINES 8
#define N_GP_COLORS 4
#define BOXCOLOR (N_GP_COLORS - 1)

typedef struct {
    int varnum;            /* ID number of variable to plot */
    char title[MAXTITLE];  /* key or legend title */
    char formula[128];     /* expression to plot (rather than data) */
    char style[16];        /* lines, points, etc. */
    char scale[8];         /* string representation of scale factor */
    int yaxis;             /* 1 for left, 2 for right */
    int type;              /* 1, 2, ... (color) */
    int width;             /* default 1, could be bigger */
    int ncols;             /* number of data columns (0 for formula) */
} GPT_LINE;

#define PLOT_LABEL_TEXT_LEN 31

typedef enum {
    GP_JUST_LEFT,
    GP_JUST_CENTER,
    GP_JUST_RIGHT
} gp_just_codes;

typedef enum {
    PLOT_REGULAR = 0,
    PLOT_H_TEST,
    PLOT_PROB_DIST,
    PLOT_FORECAST,
    PLOT_GARCH,
    PLOT_FREQ_SIMPLE,
    PLOT_FREQ_NORMAL,
    PLOT_FREQ_GAMMA,
    PLOT_PERIODOGRAM,
    PLOT_CORRELOGRAM,
    PLOT_CUSUM,
    PLOT_MULTI_SCATTER,
    PLOT_TRI_GRAPH,
    PLOT_RANGE_MEAN,
    PLOT_HURST,
    PLOT_LEVERAGE,
    PLOT_IRFBOOT,
    PLOT_KERNEL,
    PLOT_VAR_ROOTS,
    PLOT_ELLIPSE,
    PLOT_MULTI_IRF,
    PLOT_PANEL,
    PLOT_BI_GRAPH,
    PLOT_TYPE_MAX
} PlotType;

typedef struct {
    char text[PLOT_LABEL_TEXT_LEN + 1]; 
    double pos[2];
    int just;
} GPT_LABEL;

typedef struct {
    FILE *fp;
    char fname[MAXLEN];        /* for gui purposes */
    PlotType code;             /* to deal with FREQ, FCASTERR... */
    PlotSpecFlags flags;        /* bitwise OR of options */
    int nobs;                  /* number of observations */
    char titles[4][MAXTITLE];  /* main, x, y, y2 */
    double range[3][2];        /* axis range specifiers */
    char keyspec[MAXTITLE];    /* position of key (or none) */
    char xtics[16];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    char termtype[MAXTITLE];   /* gnuplot "term" setting */
    int n_lines;               /* number of lines */
    int xzeroaxis;             /* show x == 0 (1) or not (0) */
    float boxwidth;            /* when using box style for frequency plots */
    GPT_LINE *lines;           /* details on individual lines */
    char **literal;            /* additional commands */
    int n_literal;             /* number of the above */
    double *data;              /* data to plot */
    char **markers;            /* data-point markers (not always present) */
    int n_markers;             /* number of such markers */
    GPT_LABEL labels[MAX_PLOT_LABELS];  /* textual labels written onto graph */
    int *reglist;              /* regression list for X-Y plot with fitted line */
    char *labeled;             /* for GUI use */
    void *ptr;                 /* for GUI use */
} GPT_SPEC;

#define frequency_plot_code(c) (c == PLOT_FREQ_SIMPLE || \
				c == PLOT_FREQ_NORMAL || \
				c == PLOT_FREQ_GAMMA)

#define set_png_output(p) (p->flags |= GPTSPEC_PNG_OUTPUT)
#define get_png_output(p) (p->flags & GPTSPEC_PNG_OUTPUT) 
    
/* functions follow */

const char *get_gretl_png_term_line (PlotType ptype);

const char *get_gretl_emf_term_line (PlotType ptype, int color);

const char *gp_justification_string (int j);

int gnuplot_init (PlotType ptype, FILE **fpp);

PlotType plot_type_from_string (const char *str);

int gnuplot_make_graph (void);

GnuplotFlags gp_flags (int batch, gretlopt opt);

int gnuplot (const int *plotlist, const int *lines, const char *literal,
	     double ***pZ, DATAINFO *pdinfo, 
	     int *plot_count, GnuplotFlags flags);

int multi_scatters (const int *list, const double **Z,
		    const DATAINFO *pdinfo, 
		    int *plot_count, GnuplotFlags flags);

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo, 
		int *plot_count, GnuplotFlags flags);

int plot_freq (FreqDist *freq, DistCode dist);

int garch_resid_plot (const MODEL *pmod, const DATAINFO *pdinfo); 

int print_plotspec_details (const GPT_SPEC *spec, FILE *fp);

int go_gnuplot (GPT_SPEC *spec, char *fname);

void free_plotspec (GPT_SPEC *spec);

int plotspec_add_line (GPT_SPEC *spec);

int get_termstr (const GPT_SPEC *spec, char *termstr);

int rmplot (const int *list, const double **Z, DATAINFO *pdinfo, 
	    PRN *prn);

int hurstplot (const int *list, const double **Z, DATAINFO *pdinfo, 
	       PRN *prn);

int 
gretl_panel_ts_plot (const int *list, const double **Z, DATAINFO *pdinfo);

int plot_fcast_errs (int t1, int t2, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     int time_series);

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 const double **Z,
				 const DATAINFO *pdinfo);

int 
gretl_VAR_plot_multiple_irf (GRETL_VAR *var, int periods,
			     const double **Z,
			     const DATAINFO *pdinfo);

int gretl_VAR_residual_plot (const GRETL_VAR *var, const DATAINFO *pdinfo);

int gretl_VAR_residual_mplot (const GRETL_VAR *var, const DATAINFO *pdinfo); 

int gretl_VAR_roots_plot (GRETL_VAR *var);

int confidence_ellipse_plot (gretl_matrix *V, double *b, double t, double c,
			     const char *iname, const char *jname);

int is_auto_ols_string (const char *s);

int gnuplot_has_ttf (int reset);

int gnuplot_has_pdf (void);

int gnuplot_has_specified_colors (void);

void set_graph_palette (int i, const char *colstr);

void graph_palette_reset (int i);

const char *graph_color_string (int i);

int gnuplot_test_command (const char *cmd);

#ifdef ENABLE_NLS
void pprint_gnuplot_encoding (const char *termstr, PRN *prn);
void fprint_gnuplot_encoding (const char *termstr, FILE *fp);
#endif

#endif /* GRAPHING_H */

