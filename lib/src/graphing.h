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

#define MAXTITLE 80

typedef struct {
    int varnum;            /* ID number of variable to plot */
    char title[MAXTITLE];  /* key or legend title */
    char formula[128];     /* expression to plot (rather than data) */
    char style[16];        /* lines, points, etc. */
    char scale[8];         /* string repres. of scale factor */
    int yaxis;             /* 1 for left, 2 for right */
} GPT_LINE;

typedef struct {
    FILE *fp;
    char fname[MAXLEN];        /* for gui purposes */
    int edit;                  /* 1 for editing existing plot */
    int code;                  /* to deal with FREQ, FCASTERR... */
    int t1, t2;                /* starting and ending obs */
    char titles[4][MAXTITLE];  /* main, x, y, y2 */
    char range[3][2][12];      /* axis range specifiers */
    char keyspec[MAXTITLE];    /* position of key (or none) */
    char xtics[16];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    char termtype[MAXTITLE];   /* gnuplot "term" setting */
    int ts;                    /* time-series plot? (1 or 0) */ 
    int y2axis;                /* use second y axis? (1 or 0) */
    int list[8];               /* list of variables */
    GPT_LINE *lines;           /* details on individual lines */
    char *literal[4];          /* additional commands */
    double *data;              /* data to plot */
    char **labels;             /* data-point labels (not always present) */
    int nlabels;               /* number of labels */
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
    PLOT_RANGE_MEAN
} plot_type_codes;
    
#define GRETL_GUI(p) (p->binbase[0] && p->ratsbase[0] && p->dbhost_ip[0])

/* functions follow */
 
int plot (const LIST list, 
	  double **Z, const DATAINFO *pdinfo, 
	  int oflag, int pause, PRN *prn);

int graph (const LIST list, 
	   double **Z, const DATAINFO *pdinfo, 
	   int oflag, PRN *prn);

const char *get_gretl_png_term_line (void);

int gnuplot_init (PATHS *ppaths, FILE **fpp);

int gnuplot_display (const PATHS *ppaths);

int gnuplot (LIST list, const int *lines, const char *literal,
	     double ***pZ, DATAINFO *pdinfo, PATHS *ppaths, 
	     int *plot_count, int batch, int gui, int opt);

int multi_scatters (const LIST list, int pos, 
		    double ***pZ, const DATAINFO *pdinfo, 
		    PATHS *ppaths);

int plot_freq (FREQDIST *freq, PATHS *ppaths, int dist);

int print_plotspec_details (const GPT_SPEC *spec, FILE *fp);

int go_gnuplot (GPT_SPEC *spec, char *fname, PATHS *ppaths);

void free_plotspec (GPT_SPEC *spec);

int termtype_to_termstr (const char *termtype, char *termstr);

int rmplot (const LIST list, double **Z, DATAINFO *pdinfo, PRN *prn,
	    PATHS *ppaths);

int plot_fcast_errs (int n, const double *obs, 
		     const double *depvar, const double *yhat, 
		     const double *maxerr, const char *varname, 
		     PATHS *ppaths);






