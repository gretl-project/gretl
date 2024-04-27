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

#ifndef SELECTOR_H
#define SELECTOR_H

#define robust_conf(c) (c != LOGIT && c != PROBIT &&	\
                        c != OLOGIT && c != OPROBIT &&	\
                        c != QUANTREG && c != INTREG && \
                        c != MLOGIT && c != COUNTMOD && \
                        c != DURATION && c != HECKIT && \
			c != BIPROBIT && c != REPROBIT && \
			c != TOBIT && c != MIDASREG)

typedef struct iterinfo_t iterinfo;

struct iterinfo_t {
    int ci;
    int maxiters;
    double tol;
};

typedef struct gui_midas_spec_ gui_midas_spec;

struct gui_midas_spec_ {
    int nterms;
    char listname[VNAMELEN];
    int leadvar;
    int fratio;
    int ptype;
    int minlag;
    int maxlag;
    int nparm;
};

void clear_selector (void);

selector *selection_dialog (int ci, const char *title,
			    void *data, int (*callback)());

selector *simple_selection (int ci, const char *title, int (*callback)(), 
			    GtkWidget *parent);

selector *
simple_selection_for_viewer (int ci, const char *title, int (*callback)(), 
			     windata_t *vwin);

selector *
simple_selection_with_data (int ci, const char *title, int (*callback)(), 
			    GtkWidget *parent, gpointer data);

selector *
sublist_selection (int ci, const char *title, int (*callback)(),
		   GtkWidget *parent, const int *list,
		   const int *presel, void *data);

void modelspec_dialog (int ci);

void selector_set_varnum (int v);

void selector_from_model (windata_t *vwin);

void data_export_selection_wrapper (int file_code);

void functions_selection_wrapper (GtkWidget *parent);

void add_remove_functions_dialog (char **pubnames, int npub,
				  char **privnames, int npriv,
				  void *p1, void *p2);

int selector_code (const selector *sr);

const char *selector_list (const selector *sr);

int selector_list_hasconst (const selector *sr);

gpointer selector_get_data (const selector *sr);

gpointer selector_get_extra_data (const selector *sr);

gpointer selector_get_regls_bundle (void);

gretlopt selector_get_opts (const selector *sr);

int selector_get_depvar_number (const selector *sr);

int selector_get_VAR_order (const selector *sr);

const char *selector_entry_text (const selector *sr);

int selector_error (const selector *sr);

void maybe_clear_selector (const int *dlist);

GtkWidget *selector_get_window (const selector *sr);

gchar *get_selector_storelist (void);

void set_selector_storelist (const char *s);

void selector_register_genr (int newvars, gpointer p);

void selector_register_hc_choice (void);

void selector_cleanup (void);

#endif /* SELECTOR_H */
