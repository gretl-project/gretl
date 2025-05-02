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

/* plot_priv.h : private header to share symbols between graphing.c
   and any other source files that work with plots at a low level.
   As of August 2023 that just means plotbands.c.
*/

#ifndef PLOT_PRIV_H
#define PLOT_PRIV_H

typedef enum {
    BP_REGULAR,
    BP_BLOCKMAT
} BPMode;

#define GPNA "NaN"

typedef struct gnuplot_info_ gnuplot_info;

struct gnuplot_info_ {
    GptFlags flags;
    FitType fit;
    int *list;
    int t1;
    int t2;
    double xrange;
    char xtics[64];
    char xfmt[16];
    char yfmt[16];
    const char *yformula;
    const double *x;
    gretl_matrix *fvals; /* factor values */
    int n_fvals;         /* number of factor values */
    int fid;             /* factor variable ID */
    int *withlist;
    int band;
    double ybase;
};

void clear_gpinfo (gnuplot_info *gi);

int get_effective_plot_ci (void);

int graph_list_adjust_sample (gnuplot_info *gi,
                              const DATASET *dset,
                              int listmin);

void make_time_tics (gnuplot_info *gi,
                     const DATASET *dset,
                     int many, char *xlabel,
                     PRN *prn);

void print_x_range (gnuplot_info *gi, FILE *fp);

void check_for_yscale (gnuplot_info *gi,
                       const double **Z,
                       int *oddman);

void set_plot_withstr (gnuplot_info *gi, int i, char *str);

void print_user_y_data (const double *x,
                        const double *y,
                        int t1, int t2,
                        FILE *fp);

int plot_with_band (BPMode mode, gnuplot_info *gi,
                    const char *literal,
                    DATASET *dset,
                    gretlopt opt);

#endif /* PLOT_PRIV_H */
