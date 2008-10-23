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

#ifndef PLOTSPEC_H
#define PLOTSPEC_H

#define GP_MAXFORMULA 128
#define GP_MAXSTYLE    16
#define GP_MAXSCALE     8

typedef enum {
    GP_LINE_USER = 1 << 0
} gp_line_flags;

typedef struct {
    int varnum;                    /* ID number of variable to plot */
    char title[MAXTITLE];          /* key or legend title */
    char formula[GP_MAXFORMULA];   /* expression to plot (rather than data) */
    char style[GP_MAXSTYLE];       /* lines, points, etc. */
    char scale[GP_MAXSCALE];       /* string representation of scale factor */
    char yaxis;                    /* 1 for left, 2 for right */
    int type;                      /* 1, 2, ... (color) */
    int width;                     /* default 1, could be bigger */
    char ncols;                    /* number of data columns (0 for formula) */
    char flags;                    /* additional options */
} GPT_LINE;

#define PLOT_LABEL_TEXT_LEN 31

typedef enum {
    GP_JUST_LEFT,
    GP_JUST_CENTER,
    GP_JUST_RIGHT
} gp_just_codes;

typedef struct {
    char text[PLOT_LABEL_TEXT_LEN + 1]; 
    double pos[2];
    int just;
} GPT_LABEL;

typedef struct {
    FILE *fp;
    char fname[MAXLEN];        /* for gui purposes */
    PlotType code;             /* to deal with FREQ, FCASTERR... */
    GptFlags flags;            /* bitwise OR of options */
    FitType fit;               /* type of fitted line shown */
    int nobs;                  /* number of observations */
    int okobs;                 /* number of fully valid observations */
    int pd;                    /* frequency (time series data) */
    char xvarname[MAXDISP];    /* name of x variable */
    char yvarname[MAXDISP];    /* name of y variable */
    char titles[4][MAXTITLE];  /* main, x, y, y2 */
    double range[4][2];        /* axis range specifiers */
    double logbase[3];         /* axis log-scales base (0 for linear) */
    char keyspec[MAXTITLE];    /* position of key (or none) */
    char xtics[16];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    int termtype;              /* gnuplot "terminal" code */
    int n_lines;               /* number of lines */
    int samples;               /* number of samples for parametric plots */
    float boxwidth;            /* when using box style for frequency plots */
    GPT_LINE *lines;           /* details on individual lines */
    char **literal;            /* additional commands */
    int n_literal;             /* number of the above */
    double *data;              /* data to plot */
    char **markers;            /* data-point markers (not always present) */
    int n_markers;             /* number of such markers */
    GPT_LABEL labels[MAX_PLOT_LABELS];  /* textual labels written onto graph */
    int *reglist;              /* regression list for X-Y plot with fitted line */
    gretl_matrix *b_ols;       /* coeffs for linear fit */
    gretl_matrix *b_quad;      /* coeffs for quadratic fit */
    gretl_matrix *b_inv;       /* coeffs for inverse fit */
    char *labeled;             /* for GUI use */
    void *ptr;                 /* for GUI use */
} GPT_SPEC;

GPT_SPEC *plotspec_new (void);

void plotspec_destroy (GPT_SPEC *spec);

void plotspec_label_init (GPT_LABEL *lbl);

int plotspec_add_line (GPT_SPEC *spec);

int plotspec_delete_line (GPT_SPEC *spec, int i);

int plotspec_print (const GPT_SPEC *spec, FILE *fp);

int plotspec_add_fit (GPT_SPEC *spec, FitType f);

void print_auto_fit_string (FitType fit, FILE *fp);

#endif /* PLOTSPEC_H */
