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

#define GP_BORDER_DEFAULT (-1)
#define LT_AUTO (-2)

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
    GP_STYLE_CANDLESTICKS,
    GP_STYLE_AUTO
} GpLineStyle;

typedef enum {
    GP_KEY_LEFT_TOP,
    GP_KEY_RIGHT_TOP,
    GP_KEY_LEFT_BOTTOM,
    GP_KEY_RIGHT_BOTTOM,
    GP_KEY_OUTSIDE,
    GP_KEY_NONE
} GpKeyPos;

typedef enum {
    GP_X_RANGE,
    GP_Y_RANGE,
    GP_Y2_RANGE,
    GP_T_RANGE,
    GP_X2_RANGE
} GpRange;

typedef struct gp_style_spec_ gp_style_spec;

struct gp_style_spec_ {
    int id;
    const char *name;
    const char *trname;
};

typedef struct gp_key_spec_ gp_key_spec;

struct gp_key_spec_ {
    int id;
    const char *str;
};

typedef struct plotbars_ plotbars;

typedef enum {
    GP_LINE_USER    = 1 << 0,
    GP_LINE_BOXDATA = 1 << 1,
    GP_LINE_AUXDATA = 1 << 2
} gp_line_flags;

/* information about a line within a gnuplot graph */

typedef struct {
    int varnum;          /* ID number of variable to plot */
    GpLineStyle style;   /* lines, points, etc. */
    gchar *title;        /* key or legend title */
    gchar *formula;      /* expression to plot (rather than data) */
    char *ustr;          /* custom 'using' string */
    int *mcols;          /* data-matrix columns used */
    float pscale;        /* scale factor for points */
    char rgb[10];        /* rgb color specification */
    char yaxis;          /* 1 for left, 2 for right */
    int type;            /* 1, 2, ... (style reference) */
    int ptype;           /* point type */
    int dtype;           /* dash type */
    float width;         /* line width, default 1.0 */
    char ncols;          /* number of data columns (0 for formula) */
    float whiskwidth;    /* whiskerbar width (boxplots) */
    char flags;          /* additional options */
} GPT_LINE;

#define PLOT_LABEL_TEXT_LEN 31

typedef enum {
    GP_JUST_LEFT,
    GP_JUST_CENTER,
    GP_JUST_RIGHT
} gp_just_codes;

typedef struct {
    char text[PLOT_LABEL_TEXT_LEN + 1];
    double pos[2]; /* x,y position on plot */
    guint8 just;   /* justification of text, see gp_just_codes */
} GPT_LABEL;

typedef enum {
    GP_ARROW_HEAD = 1 << 0,
    GP_ARROW_DOTS = 1 << 1
} gp_arrow_flags;

typedef struct {
    double x0, y0, x1, y1;
    int flags;
} GPT_ARROW;

/* "global" information concerning a gnuplot graph specification */

struct GPT_SPEC_ {
    FILE *fp;
    char fname[MAXLEN];        /* for gui purposes */
    PlotType code;             /* to deal with FREQ, FCAST... */
    GptFlags flags;            /* bitwise OR of options */
    int heredata;              /* plot file uses "heredoc" for data */
    FitType fit;               /* type of fitted line shown */
    int nobs;                  /* number of observations */
    int okobs;                 /* number of fully valid observations */
    int *datacols;             /* list of data columns */
    int pd;                    /* frequency (time series data) */
    int nbars;                 /* number of time-series shaded bars */
    char xvarname[MAXDISP];    /* name of x variable */
    char yvarname[MAXDISP];    /* name of y variable */
    gchar *titles[5];          /* main, x, y, y2, x2 */
    double range[5][2];        /* axis range specifiers: x, y, y2, t, x2 */
    double logbase[3];         /* axis log-scales base (0 for linear) */
    int keyspec;               /* position of key (or none) */
    char xfmt[16];             /* x-axis tic format */
    char xtics[64];            /* x-axis tic marks */
    char mxtics[4];            /* minor tics */
    char yfmt[16];             /* y-axis tic format */
    char ytics[64];            /* y-axis tic marks */
    char *xticstr;             /* special xtics */
    char *x2ticstr;            /* special x2tics */
    int termtype;              /* gnuplot "terminal" code */
    int n_lines;               /* number of lines */
    int samples;               /* number of samples for parametric plots */
    int border;                /* gnuplot border code */
    char border_lc[10];        /* border line color */
    int bmargin;               /* bottom margin */
    float boxwidth;            /* when using box style for frequency plots */
    float fillfrac;            /* fraction for "set style fill solid" */
    GPT_LINE *lines;           /* details on individual lines */
    char **literal;            /* additional commands */
    int n_literal;             /* number of the above */
    char **multi_xtics;        /* array of xtics strings and positions */
    int n_xtics;               /* number of the above */
    gretl_matrix *data;        /* data to plot */
    char **markers;            /* data-point markers (not always present) */
    int n_markers;             /* number of such markers */
    GPT_LABEL *labels;         /* labels written onto graph */
    int n_labels;              /* number of the above */
    GPT_ARROW *arrows;         /* arrows drawn onto graph */
    int n_arrows;              /* number of the above */
    int *reglist;              /* regression list for X-Y plot with fitted line */
    gretl_matrix *b_ols;       /* coeffs for linear fit */
    gretl_matrix *b_quad;      /* coeffs for quadratic fit */
    gretl_matrix *b_cub;       /* coeffs for cubic fit */
    gretl_matrix *b_inv;       /* coeffs for inverse fit */
    gretl_matrix *b_log;       /* coeffs for semi-log fit */
    gretl_matrix *b_linlog;    /* coeffs for linear-log fit */
    double scale;              /* for GUI use */
    char *labeled;             /* for GUI use */
    void *ptr;                 /* for GUI use */
    plotbars *bars;            /* for GUI use */
    char *fontstr;             /* for GUI use */
    gretl_matrix *auxdata;     /* for GUI use */
};

GPT_SPEC *plotspec_new (void);

void plotspec_destroy (GPT_SPEC *spec);

void plotspec_label_init (GPT_LABEL *lbl);

int plotspec_add_line (GPT_SPEC *spec);

int plotspec_delete_line (GPT_SPEC *spec, int i);

int plotspec_max_line_width (GPT_SPEC *spec);

GPT_LINE *plotspec_clone_lines (GPT_SPEC *spec, int *err);

int plotspec_add_arrow (GPT_SPEC *spec);

int plotspec_delete_arrow (GPT_SPEC *spec, int i);

GPT_ARROW *plotspec_clone_arrows (GPT_SPEC *spec, int *err);

int plotspec_add_label (GPT_SPEC *spec);

int plotspec_delete_label (GPT_SPEC *spec, int i);

GPT_LABEL *plotspec_clone_labels (GPT_SPEC *spec, int *err);

int plotspec_print (GPT_SPEC *spec, FILE *fp);

int plotspec_add_fit (GPT_SPEC *spec, FitType f);

void print_auto_fit_string (FitType fit, FILE *fp);

int gp_line_data_columns (GPT_SPEC *spec, int i);

const char *gp_line_style_display_name (int t);

int gp_style_index_from_name (const char *s);

int gp_style_index_from_display_name (const char *s);

gp_style_spec *get_style_spec (int t);

int gp_keypos_from_name (const char *s);

int gp_keypos_from_display_name (const char *s);

gp_key_spec *get_keypos_spec (int t);

void print_keypos_string (int t, FILE *fp);

void set_plotfit_line (GPT_LINE *line,
		       FitType f, const double *b,
		       double x0, double pd);

int plotspec_add_bars_info (GPT_SPEC *spec,
			    double xmin, double xmax,
			    double ymin, double ymax,
			    const char *fname);

int plotspec_allocate_bars (GPT_SPEC *spec);

int plotspec_set_bar_info (GPT_SPEC *spec, int i,
			   double t1, double t2);

void plotspec_set_bars_limits (GPT_SPEC *spec,
			       double t1, double t2,
			       double ymin, double ymax);

void plotspec_remove_bars (GPT_SPEC *spec);

int plotspec_line_is_formula (const GPT_SPEC *spec, int i);

#endif /* PLOTSPEC_H */
