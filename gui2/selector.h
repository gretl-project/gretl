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

#ifndef SELECTOR_H
#define SELECTOR_H

#ifdef ENABLE_GMP
#define MODEL_CODE(c) (c == OLS || c == CORC || c == HILU || c == WLS || \
                       c == POOLED || c == HCCM || c == HSK || c == ARMA || \
                       c == TSLS || c == LOGIT || c == PROBIT || c == GARCH || \
                       c == AR || c == VAR || c == COINT || c == COINT2 || \
                       c == MPOLS || c == LAD || c == LOGISTIC || c == TOBIT || \
                       c == PWE)
#else
#define MODEL_CODE(c) (c == OLS || c == CORC || c == HILU || c == WLS || \
                       c == POOLED || c == HCCM || c == HSK || c == ARMA || \
                       c == TSLS || c == LOGIT || c == PROBIT || c == GARCH || \
                       c == AR || c == VAR || c == COINT || c == COINT2 || \
                       c == LAD || c == LOGISTIC || c == TOBIT || c == PWE)
#endif

#define ADDVAR_CODE(c) (c == LOGS || c == LAGS || c == SQUARE || \
                        c == DIFF || c == LDIFF)

#define GRAPH_CODE(c) (c == GR_PLOT || c == GR_XY || c == GR_IMP || GR_DUMMY)

typedef struct _selector selector;

void clear_selector (void);

void delete_selection_dialog (selector *sr);

void selection_dialog (const char *title, void (*okfunc)(), guint cmdcode);

void simple_selection (const char *title, void (*okfunc)(), guint cmdcode,
		       gpointer p);

char *mdata_selection_to_string (int n_required);

void data_save_selection_wrapper (int file_code);

int selector_code (const selector *sr);

const char *selector_list (const selector *sr);

gpointer selector_get_data (const selector *sr);

unsigned long selector_get_opts (const selector *sr);

int selector_error (const selector *sr);

#endif /* SELECTOR_H */
