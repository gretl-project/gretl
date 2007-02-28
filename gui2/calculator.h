/*
 *  Copyright (c) by Allin Cottrell
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

#ifndef CALCULATOR_H
#define CALCULATOR_H

typedef enum {
    F_BINV,
    F_CHI,
    F_LOG2,
    F_BIGCHI,
    F_F
} FormulaCode;

void stats_calculator (gpointer p, guint code, GtkWidget *w);

double dist_xmax (int d, int df1, int df2);

gchar *dist_marker_line (int dist, int df1, int df2);

const char *dist_formula (FormulaCode c);

gchar *dist_graph_title (int dist, double x, int df1, int df2);

#endif /* CALCULATOR_H */
