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

#define MODEL_CODE(c) (c == OLS || c == CORC || c == HILU || c == WLS || \
                       c == POOLED || c == HCCM || c == HSK || \
                       c == TSLS || c == LOGIT || c == PROBIT || \
                       c == AR || c == VAR || c == COINT)

#define ADDVAR_CODE(c) (c == LOGS || c == LAGS || c == SQUARE || \
                        c == DIFF || c == LDIFF)

typedef struct {
    GtkWidget *dlg;
    GtkWidget *varlist;
    GtkWidget *depvar;
    GtkWidget *rightvars;
    GtkWidget *default_check;
    GtkWidget *extra;
    int code;
    char *cmdlist;
    gpointer data;
} selector;

void selection_dialog (const char *title, const char *oktxt, 
		       void (*okfunc)(), guint cmdcode);

void simple_selection (const char *title, const char *oktxt, 
		       void (*okfunc)(), guint cmdcode,
		       gpointer p);
