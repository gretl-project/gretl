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

#define GP_BORDER_DEFAULT (-1)
#define LT_NONE (-2)

typedef enum {
    GP_STYLE_NONE,
    GP_STYLE_LINES,
    GP_STYLE_POINTS,
    GP_STYLE_LINESPOINTS,
    GP_STYLE_IMPULSES,
    GP_STYLE_DOTS,
    GP_STYLE_STEPS,
    GP_STYLE_BOXES,
    GP_STYLE_ERRORBARS,
    GP_STYLE_FILLEDCURVE,
    GP_STYLE_CANDLESTICKS
} GpLineStyle;

typedef enum {
    GP_KEY_LEFT_TOP,
    GP_KEY_RIGHT_TOP,
    GP_KEY_LEFT_BOTTOM,
    GP_KEY_RIGHT_BOTTOM,
    GP_KEY_OUTSIDE,
    GP_KEY_NONE
} GpKeyPos;

typedef struct gp_style_spec_ gp_style_spec;

struct gp_style_spec_ {
    int sty;
    const char *str;
};

typedef enum {
    GP_LINE_USER    = 1 << 0,
    GP_LINE_BOXDATA = 1 << 1
} gp_line_flags;

/* information about a line within a gnuplot graph */

typedef struct {
    int varnum;                    /* ID number of variable to plot */
    int style;                     /* lines, points, etc. */
    char title[MAXTITLE];          /* key or legend title */
    char formula[GP_MAXFORMULA];   /* expression to plot (rather than data) */
    double scale;                  /* scale factor for data */
    char rgb[8];                   /* rgb color specification */
    char yaxis;                    /* 1 for left, 2 for right */
    int type;                      /* 1, 2, ... (style) */
    int ptype;                     /* point type */
    int width;                     /* default 1, could be bigger */
    char ncols;                    /* number of data columns (0 for formula) */
    float whiskwidth;              /* whiskerbar width (boxplots) */
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

/* "global" information concerning a gnuplot graph specification */

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
    int keyspec;               /* position of key (or none) */
    char xtics[16];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    int termtype;              /* gnuplot "terminal" code */
    int n_lines;               /* number of lines */
    int samples;               /* number of samples for parametric plots */
    int border;                /* gnuplot border code */
    int bmargin;               /* bottom margin */
    float boxwidth;            /* when using box style for frequency plots */
    GPT_LINE *lines;           /* details on individual lines */
    char **literal;            /* additional commands */
    int n_literal;             /* number of the above */
    double *data;              /* data to plot */
    char **markers;            /* data-point markers (not always present) */
    int n_markers;             /* number of such markers */
    GPT_LABEL *labels;         /* textual labels written onto graph */
    int n_labels;              /* number of the above */
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

GPT_LINE *plotspec_clone_lines (GPT_SPEC *spec, int *err);

int plotspec_add_label (GPT_SPEC *spec);

int plotspec_delete_label (GPT_SPEC *spec, int i);

GPT_LABEL *plotspec_clone_labels (GPT_SPEC *spec, int *err);

int plotspec_print (const GPT_SPEC *spec, FILE *fp);

int plotspec_add_fit (GPT_SPEC *spec, FitType f);

void print_auto_fit_string (FitType fit, FILE *fp);

const char *gp_line_style_string (int t);

int gp_style_from_string (const char *s);

int gp_style_from_translation (const char *s);

gp_style_spec *get_style_spec (int t);

int gp_keypos_from_string (const char *s);

int gp_keypos_from_translation (const char *s);

gp_style_spec *get_keypos_spec (int t);

void print_keypos_string (int t, FILE *fp);

#endif /* PLOTSPEC_H */
