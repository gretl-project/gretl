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

/* lib.h for gretl */

#ifndef LIB_H
#define LIB_H

/* general purpose functions */

char *user_fopen (const char *fname, char *fullname, print_t **pprn);

gint bufopen (print_t **pprn);

gint check_cmd (char *line);
 
gint cmd_init (char *line);

gint dump_cmd_stack (char *fname);

void do_menu_op (gpointer data, guint action, GtkWidget *widget);

void do_dialog_cmd (GtkWidget *widget, dialog_t *ddata);

void view_log (void);

void console (void);

void gui_errmsg (const int errcode, const char *msg);

/* sample-related functions */

void change_sample (GtkWidget *widget, dialog_t *ddata);

void do_sampledum (GtkWidget *widget, dialog_t *ddata);

void do_setobs (GtkWidget *widget, dialog_t *ddata);

void count_missing (void);

void do_add_markers (GtkWidget *widget, dialog_t *ddata);

/* model-related functions */

void do_forecast (GtkWidget *widget, dialog_t *ddata);

void do_add_omit (GtkWidget *widget, dialog_t *ddata);

void do_lmtest (gpointer data, guint aux_code, GtkWidget *widget);

void do_chow (GtkWidget *widget, dialog_t *ddata);

void do_cusum (gpointer data, guint u, GtkWidget *widget);

void do_arch (GtkWidget *widget, dialog_t *ddata);

void set_storelist (GtkWidget *widget, dialog_t *ddata);

void do_model (GtkWidget *widget, dialog_t *ddata);

/* variable-related functions */

void do_sim (GtkWidget *widget, dialog_t *ddata);

void do_simdata (GtkWidget *widget, dialog_t *ddata);

void do_genr (GtkWidget *widget, dialog_t *ddata);

void do_model_genr (GtkWidget *widget, dialog_t *ddata);

void do_random (GtkWidget *widget, dialog_t *ddata);

void do_seed (GtkWidget *widget, dialog_t *ddata);

void do_rename_var (GtkWidget *widget, dialog_t *ddata);

void delete_var (void);

void do_edit_label (GtkWidget *widget, dialog_t *ddata);

void do_resid_freq (gpointer data, guint action, GtkWidget *widget);

void do_freqplot (gpointer data, guint gamma, GtkWidget *widget);

void do_pergm (gpointer data, guint opt, GtkWidget *widget);

void do_outcovmx (gpointer data, guint action, GtkWidget *widget);

void add_dummies (gpointer data, guint action, GtkWidget *widget);

void add_time (gpointer data, guint index, GtkWidget *widget);

void add_logs_etc (GtkWidget *widget, dialog_t *ddata);

int add_fit_resid (MODEL *pmod, const int code, const int undo);

void add_model_stat (MODEL *pmod, const int which);

void resid_plot (gpointer data, guint xvar, GtkWidget *widget);

void fit_actual_plot (gpointer data, guint xvar, GtkWidget *widget);

void display_data (gpointer data, guint select, GtkWidget *widget);

void display_fit_resid (gpointer data, guint code, GtkWidget *widget);

void do_graph_var (void);

void do_boxplot_var (void);

void do_scatters (GtkWidget *widget, dialog_t *ddata);

void do_graph (GtkWidget *widget, dialog_t *ddata);

void do_box_graph (GtkWidget *widget, dialog_t *ddata);

void do_dummy_graph (GtkWidget *widget, dialog_t *ddata);

void display_selected (GtkWidget *widget, dialog_t *ddata);

void display_var (void);

/* script- and file-related functions */

void do_open_script (GtkWidget *w, GtkFileSelection *fs);

void view_script_default (void);

void do_new_script (gpointer data, guint action, GtkWidget *widget);

void do_open_csv_box (char *fname, int code);

void do_store (char *mydatfile, const int fmt);

void view_latex (gpointer data, guint prn_code, GtkWidget *widget);

void do_save_tex (char *fname, const int code, MODEL *pmod);

void do_save_html (char *fname, const int code, MODEL *pmod);

void do_save_text (char *fname, MODEL *pmod);

int execute_script (const char *runfile, 
		    SESSION *psession, session_t *rebuild,
		    print_t *prn, int exec_code);

#endif /* LIB_H */

