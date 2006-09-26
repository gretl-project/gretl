/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* Durbin-Watson table for gretl */

#include "libgretl.h"

#define NDVAL 12

typedef struct {
    int n;
    double dval[NDVAL];
} dw_t;

dw_t dw_vals[] = {
    { 15,{1.08,1.36, 0.95,1.54, 0.82,1.75, 0.69,1.97, 0.56,2.21, 0.00,0.00}},
    { 16,{1.10,1.37, 0.98,1.54, 0.86,1.73, 0.74,1.93, 0.62,2.15, 0.16,3.30}},
    { 17,{1.13,1.38, 1.02,1.54, 0.90,1.71, 0.78,1.90, 0.67,2.10, 0.20,3.18}},
    { 18,{1.16,1.39, 1.05,1.53, 0.93,1.69, 0.82,1.87, 0.71,2.06, 0.24,3.07}},
    { 19,{1.18,1.40, 1.08,1.53, 0.97,1.68, 0.86,1.85, 0.75,2.02, 0.29,2.97}},
    { 20,{1.20,1.41, 1.10,1.54, 1.00,1.68, 0.90,1.83, 0.79,1.99, 0.34,2.89}},
    { 21,{1.22,1.42, 1.13,1.54, 1.03,1.67, 0.93,1.81, 0.83,1.96, 0.38,2.81}},
    { 22,{1.24,1.43, 1.15,1.54, 1.05,1.66, 0.96,1.80, 0.86,1.94, 0.42,2.73}},
    { 23,{1.26,1.44, 1.17,1.54, 1.08,1.66, 0.99,1.79, 0.90,1.92, 0.47,2.67}},
    { 24,{1.27,1.45, 1.19,1.55, 1.10,1.66, 1.01,1.78, 0.93,1.90, 0.51,2.61}},
    { 25,{1.29,1.45, 1.21,1.55, 1.12,1.66, 1.04,1.77, 0.95,1.89, 0.54,2.57}},
    { 26,{1.30,1.46, 1.22,1.55, 1.14,1.65, 1.06,1.76, 0.98,1.88, 0.58,2.51}},
    { 27,{1.32,1.47, 1.24,1.56, 1.16,1.65, 1.08,1.76, 1.01,1.86, 0.62,2.47}},
    { 28,{1.33,1.48, 1.26,1.56, 1.18,1.65, 1.10,1.75, 1.03,1.85, 0.65,2.43}},
    { 29,{1.34,1.48, 1.27,1.56, 1.20,1.65, 1.12,1.74, 1.05,1.84, 0.68,2.40}},
    { 30,{1.35,1.49, 1.28,1.57, 1.21,1.65, 1.14,1.74, 1.07,1.83, 0.71,2.36}},
    { 31,{1.36,1.50, 1.30,1.57, 1.23,1.65, 1.16,1.74, 1.09,1.83, 0.74,2.33}},
    { 32,{1.37,1.50, 1.31,1.57, 1.24,1.65, 1.18,1.73, 1.11,1.82, 0.77,2.31}},
    { 33,{1.38,1.51, 1.32,1.58, 1.26,1.65, 1.19,1.73, 1.13,1.81, 0.80,2.28}},
    { 34,{1.39,1.51, 1.33,1.58, 1.27,1.65, 1.21,1.73, 1.15,1.81, 0.82,2.26}},
    { 35,{1.40,1.52, 1.34,1.53, 1.28,1.65, 1.22,1.73, 1.16,1.80, 0.85,2.24}},
    { 36,{1.41,1.52, 1.35,1.59, 1.29,1.65, 1.24,1.73, 1.18,1.80, 0.87,2.22}},
    { 37,{1.42,1.53, 1.36,1.59, 1.31,1.66, 1.25,1.72, 1.19,1.80, 0.89,2.20}},
    { 38,{1.43,1.54, 1.37,1.59, 1.32,1.66, 1.26,1.72, 1.21,1.79, 0.91,2.18}},
    { 39,{1.43,1.54, 1.38,1.60, 1.33,1.66, 1.27,1.72, 1.22,1.79, 0.93,2.16}},
    { 40,{1.44,1.54, 1.39,1.60, 1.34,1.66, 1.29,1.72, 1.23,1.79, 0.95,2.15}},
    { 45,{1.48,1.57, 1.43,1.62, 1.38,1.67, 1.34,1.72, 1.29,1.78, 1.04,2.09}},
    { 50,{1.50,1.59, 1.46,1.63, 1.42,1.67, 1.38,1.72, 1.34,1.77, 1.11,2.04}},
    { 55,{1.53,1.60, 1.49,1.64, 1.45,1.68, 1.41,1.72, 1.38,1.77, 1.17,2.01}},
    { 60,{1.55,1.62, 1.51,1.65, 1.48,1.69, 1.44,1.73, 1.41,1.77, 1.22,1.98}},
    { 65,{1.57,1.63, 1.54,1.66, 1.50,1.70, 1.47,1.73, 1.44,1.77, 1.27,1.96}},
    { 70,{1.58,1.64, 1.55,1.67, 1.52,1.70, 1.49,1.74, 1.46,1.77, 1.30,1.95}},
    { 75,{1.60,1.65, 1.57,1.68, 1.54,1.71, 1.51,1.74, 1.49,1.77, 1.34,1.94}},
    { 80,{1.61,1.66, 1.59,1.69, 1.56,1.72, 1.53,1.74, 1.51,1.77, 1.37,1.93}},
    { 85,{1.62,1.67, 1.60,1.70, 1.57,1.72, 1.55,1.75, 1.52,1.77, 1.40,1.92}},
    { 90,{1.63,1.68, 1.61,1.70, 1.59,1.73, 1.57,1.75, 1.54,1.78, 1.42,1.91}},
    { 95,{1.64,1.69, 1.62,1.71, 1.60,1.73, 1.58,1.75, 1.56,1.78, 1.44,1.90}},
    {100,{1.65,1.69, 1.63,1.72, 1.61,1.74, 1.59,1.76, 1.57,1.78, 1.46,1.90}}
};

static void other_tables (PRN *prn)
{
    pputs(prn, _("\nFor more comprehensive statistical tables, please consult "
		 "a statistics or\neconometrics text, e.g. Ramanathan's "
		 "Introductory Econometrics.\n"));
}

void dw_lookup (int n, PRN *prn)
{
    int ndw = sizeof dw_vals / sizeof dw_vals[0];
    int dist, mindist = 1000;
    int row = 0;
    int i, j;

    if (n < 15) {
	n = 15;
    } else if (n > 100) {
	n = 100;
    }

    for (i=0; i<ndw; i++) {
	dist = abs(dw_vals[i].n - n);
	if (dist == 0) {
	    row = i;
	    break;
	} else if (dist < mindist) {
	    mindist = dist;
	    row = i;
	} else if (dist >= mindist) {
	    break;
	}
    }

    pprintf(prn, "%s, n = %d\n\n",
	    /* xgettext:no-c-format */
	    _("5% critical values for Durbin-Watson statistic"), 
	    dw_vals[row].n);

    pprintf(prn, "%s:\n\n", 
	    _("       Number of explanatory variables (excluding the "
		 "constant)"));
    pputs(prn, "      1           2           3           4"
	  "           5          10\n");
    pputs(prn, "   dL   dU     dL   dU     dL   dU     dL   dU"
	  "     dL   dU     dL   dU\n\n");

    for (j=0; j<NDVAL; j++) {
	if (dw_vals[row].dval[j] == 0.0) {
	    break;
	}
	if (j % 2 == 0) {
	    pprintf(prn, "%6.2f ", dw_vals[row].dval[j]);
	} else {
	    pprintf(prn, "%4.2f ", dw_vals[row].dval[j]);
	}
    }

    pputc(prn, '\n');

    other_tables(prn);
}





