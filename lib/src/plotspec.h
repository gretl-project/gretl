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

#ifndef PLOTSPEC_H
#define PLOTSPEC_H

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
    char xvarname[MAXDISP];    /* name of x variable */
    char yvarname[MAXDISP];    /* name of y variable */
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

int plotspec_print (const GPT_SPEC *spec, FILE *fp);

int plotspec_get_term_string (const GPT_SPEC *spec, char *termstr);

int plotspec_ship_out (GPT_SPEC *spec, char *fname);

int plotspec_add_fit (GPT_SPEC *spec, FitType f);

#endif /* PLOTSPEC_H */
