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

#include <stdio.h>

#define GRAPH_NO_DATA -999

enum gnuplot_flags {
    GP_IMPULSES = 1 << 0,  /* use impulses for plotting */
    GP_RESIDS   = 1 << 1,  /* doing residual plot */
    GP_FA       = 1 << 2,  /* doing fitted/actual plot */
    GP_DUMMY    = 1 << 3,  /* using a dummy for separation */
    GP_BATCH    = 1 << 4,  /* working in batch mode */
    GP_GUI      = 1 << 5,  /* called from GUI context */
    GP_OLS_OMIT = 1 << 6,  /* Don't draw fitted line on graph */
    GP_FILE     = 1 << 7   /* send output to named file */
};

enum gptspec_flags {
    GPTSPEC_TS            = 1 << 0,
    GPTSPEC_Y2AXIS        = 1 << 1,
    GPTSPEC_AUTO_OLS      = 1 << 2,
    GPTSPEC_OLS_HIDDEN    = 1 << 3,
    GPTSPEC_BORDER_HIDDEN = 1 << 4
}; 

#define MAXTITLE 128
#define MAX_PLOT_LABELS 3

typedef struct {
    int varnum;            /* ID number of variable to plot */
    char title[MAXTITLE];  /* key or legend title */
    char formula[128];     /* expression to plot (rather than data) */
    char style[16];        /* lines, points, etc. */
    char scale[8];         /* string repres. of scale factor */
    int yaxis;             /* 1 for left, 2 for right */
} GPT_LINE;

#define PLOT_LABEL_TEXT_LEN 31
#define PLOT_LABEL_JUST_LEN  7
#define PLOT_LABEL_POS_LEN  31

typedef struct {
    char text[PLOT_LABEL_TEXT_LEN + 1]; 
    char just[PLOT_LABEL_JUST_LEN + 1];
    char pos[PLOT_LABEL_POS_LEN + 1];
} GPT_LABEL;

typedef struct {
    FILE *fp;
    char fname[MAXLEN];        /* for gui purposes */
    int edit;                  /* 1 for editing existing plot */
    int code;                  /* to deal with FREQ, FCASTERR... */
    unsigned char flags;       /* bitwise OR of options (gptspec_flags) */
    int t1, t2;                /* starting and ending obs */
    char titles[4][MAXTITLE];  /* main, x, y, y2 */
    char range[3][2][12];      /* axis range specifiers */
    char keyspec[MAXTITLE];    /* position of key (or none) */
    char xtics[16];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    char termtype[MAXTITLE];   /* gnuplot "term" setting */
    int nlines;                /* number of lines */
    int n_y_series;            /* number of actual y-axis data series (not formulae) */
    GPT_LINE *lines;           /* details on individual lines */
    char *literal[4];          /* additional commands */
    double *data;              /* data to plot */
    char **labels;             /* data-point labels (not always present) */
    int nlabels;               /* number of labels */
    GPT_LABEL text_labels[MAX_PLOT_LABELS];  /* textual labels written onto graph */
    void *ptr;                 /* for GUI use */
} GPT_SPEC;

typedef enum {
    NORMAL = 1,
    GAMMA
} dist_codes;

typedef enum {
    PLOT_REGULAR = 0,
    PLOT_SAMPLING_DIST,
    PLOT_FORECAST,
    PLOT_FREQ_SIMPLE,
    PLOT_FREQ_NORMAL,
    PLOT_FREQ_GAMMA,
    PLOT_PERIODOGRAM,
    PLOT_CORRELOGRAM,
    PLOT_CUSUM,
    PLOT_MULTI_SCATTER,
    PLOT_TRI_GRAPH,
    PLOT_RANGE_MEAN,
    PLOT_LEVERAGE,
    PLOT_TYPE_MAX
} plot_type_codes;
    
#define GRETL_GUI(p) (p->binbase[0] && p->ratsbase[0] && p->dbhost_ip[0])

/* functions follow */
 
int plot (const LIST list, 
	  double **Z, const DATAINFO *pdinfo, 
	  unsigned long oflag, int pause, PRN *prn);

int graph (const LIST list, 
	   double **Z, const DATAINFO *pdinfo, 
	   unsigned long oflag, PRN *prn);

const char *get_gretl_png_term_line (const PATHS *ppaths, int plottype);

int gnuplot_init (PATHS *ppaths, int plottype, FILE **fpp);

int gnuplot_display (const PATHS *ppaths);

int gnuplot (LIST list, const int *lines, const char *literal,
	     double ***pZ, DATAINFO *pdinfo, PATHS *ppaths, 
	     int *plot_count, unsigned char flags);

int multi_scatters (const LIST list, int pos, 
		    double ***pZ, const DATAINFO *pdinfo, 
		    PATHS *ppaths, int *plot_count, 
		    unsigned char flags);

int gnuplot_3d (LIST list, const char *literal,
		double ***pZ, DATAINFO *pdinfo, PATHS *ppaths, 
		int *plot_count, unsigned char flags);

int plot_freq (FREQDIST *freq, PATHS *ppaths, int dist);

int print_plotspec_details (const GPT_SPEC *spec, FILE *fp);

int go_gnuplot (GPT_SPEC *spec, char *fname, PATHS *ppaths);

void free_plotspec (GPT_SPEC *spec);

int termtype_to_termstr (const char *termtype, char *termstr,
			 const PATHS *ppaths);

int rmplot (const LIST list, double **Z, DATAINFO *pdinfo, PRN *prn,
	    PATHS *ppaths);

int plot_fcast_errs (int n, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     int time_series, PATHS *ppaths);

int 
gretl_var_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 const DATAINFO *pdinfo,
				 PATHS *ppaths);

int is_auto_ols_string (const char *s);

int gnuplot_has_ttf (void);

int gnuplot_has_specified_colors (void);

void set_gnuplot_pallette (int i, const char *colstr);

const char *get_gnuplot_pallette (int i, int plottype);

int gnuplot_test_command (const char *cmd);

