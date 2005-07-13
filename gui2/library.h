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

/* library.h for gretl */

#ifndef LIBRARY_H
#define LIBRARY_H

enum {
    LATEX_OK,
    LATEX_EXEC_FAILED,
    LATEX_ERROR
} tex_return_codes;

/* general purpose functions */

#ifdef G_OS_WIN32
void win_show_error (DWORD dw);

int winfork (char *cmdline, const char *dir, int wshow,
	     DWORD flags);
#endif

void library_command_init (void);

void library_command_free (void);

int replaying (void);

void set_replay_on (void);

void set_replay_off (void);

int gretl_command_sprintf (const char *template, ...);

int gretl_command_strcpy (const char *s);

int *command_list_from_string (char *s);

int user_fopen (const char *fname, char *fullname, PRN **pprn);

gint bufopen (PRN **pprn);

void do_menu_op (gpointer data, guint action, GtkWidget *widget);

void do_two_var_test (GtkWidget *widget, gpointer p);

void do_run_script (gpointer data, guint code, GtkWidget *w);

void gui_errmsg (const int errcode);

void register_graph (void);

void clear_data (void);

void exit_free_modelspec (void);

/* sample-related functions */

int bool_subsample (gretlopt opt);

void do_samplebool (GtkWidget *widget, dialog_t *dlg);

int do_set_sample (void);

void drop_all_missing (gpointer data, guint opt, GtkWidget *w);

void count_missing (void);

void do_add_markers (const char *fname);

void do_remove_markers (gpointer data, guint u, GtkWidget *w);

int dataset_is_restricted (void);

int maybe_restore_full_data (int action);

void gui_transpose_data (gpointer p, guint u, GtkWidget *w);

int dataset_is_subsampled (void);

void set_original_n (int n);

int get_original_n (void);

/* model-related functions */

void do_coint (GtkWidget *widget, gpointer p);

void do_coint2 (GtkWidget *widget, gpointer p);

void do_forecast (gpointer data, guint u, GtkWidget *w);

void do_coeff_sum (GtkWidget *widget, gpointer p);

void do_add_omit (GtkWidget *widget, gpointer p);

void do_lmtest (gpointer data, guint aux_code, GtkWidget *widget);

void do_autocorr (gpointer data, guint u, GtkWidget *widget);

void do_chow_cusum (gpointer data, guint action, GtkWidget *w);

void do_reset (gpointer data, guint u, GtkWidget *widget);

void unit_root_test (gpointer data, guint u, GtkWidget *widget);

void do_arch (gpointer data, guint u, GtkWidget *widget);

void do_restrict (GtkWidget *widget, dialog_t *dlg);

void do_nls_model (GtkWidget *widget, dialog_t *dlg);

void do_model (GtkWidget *widget, gpointer p);

void do_kernel (gpointer data, guint u, GtkWidget *w);

void do_vif (gpointer data, guint u, GtkWidget *w);

void do_leverage (gpointer data, guint u, GtkWidget *w);

void add_leverage_data (windata_t *vwin);

void do_coeff_intervals (gpointer data, guint i, GtkWidget *w);

void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w);

void do_spearman (GtkWidget *widget, gpointer p);

#ifdef ENABLE_GMP
void do_mp_ols (GtkWidget *widget, gpointer p);
#endif

int out_of_sample_info (int add_ok, int *t2);

/* variable-related functions */

void do_simdata (GtkWidget *widget, dialog_t *dlg);

void do_genr (GtkWidget *widget, dialog_t *dlg);

void do_model_genr (GtkWidget *widget, dialog_t *dlg);

void do_random (GtkWidget *widget, dialog_t *dlg);

void do_seed (GtkWidget *widget, dialog_t *dlg);

void do_global_setmiss (GtkWidget *widget, dialog_t *dlg);

void do_variable_setmiss (GtkWidget *widget, dialog_t *dlg);

void do_edit_label (GtkWidget *widget, dialog_t *dlg);

int do_rename_variable (int v, const char *newname, int full);

int record_varlabel_change (int v);

void do_resid_freq (gpointer data, guint action, GtkWidget *widget);

void do_freqplot (gpointer data, guint gamma, GtkWidget *widget);

void do_corrgm (gpointer data, guint u, GtkWidget *widget);

void do_pergm (gpointer data, guint opt, GtkWidget *widget);

#if defined (HAVE_TRAMO) || defined (HAVE_X12A)
void do_tramo_x12a (gpointer data, guint opt, GtkWidget *widget);
#endif

void do_range_mean (gpointer data, guint opt, GtkWidget *widget);

void do_hurst (gpointer data, guint opt, GtkWidget *widget);

void do_outcovmx (gpointer data, guint action, GtkWidget *widget);

void add_dummies (gpointer data, guint action, GtkWidget *widget);

void add_index (gpointer data, guint tm, GtkWidget *widget);

void do_add_obs (gpointer data, guint u, GtkWidget *widget);

void do_remove_obs (gpointer data, guint u, GtkWidget *widget);

void add_logs_etc (gpointer data, guint action, GtkWidget *widget);

int add_fit_resid (MODEL *pmod, int code, int undo);

int add_var_resid (GRETL_VAR *var, int eqnum);

void add_model_stat (MODEL *pmod, int which);

void resid_plot (gpointer data, guint xvar, GtkWidget *widget);

void fit_actual_plot (gpointer data, guint xvar, GtkWidget *widget);

void fit_actual_splot (gpointer data, guint u, GtkWidget *widget);

void display_data (gpointer data, guint select, GtkWidget *widget);

void display_fit_resid (gpointer data, guint code, GtkWidget *widget);

void do_graph_var (int varnum);

void do_boxplot_var (int varnum);

void ts_plot_var (gpointer data, guint opt, GtkWidget *widget);

void do_scatters (GtkWidget *widget, gpointer p);

void do_graph_from_selector (GtkWidget *widget, gpointer p);

void do_splot_from_selector (GtkWidget *widget, gpointer p);

void plot_from_selection (gpointer data, guint action, GtkWidget *widget);

void do_box_graph (GtkWidget *widget, dialog_t *dlg);

void do_dummy_graph (GtkWidget *widget, gpointer p);

void delete_selected_vars (int id);

void display_selected (gpointer data, guint action, GtkWidget *widget);

void display_var (void);

/* script- and file-related functions */

void do_open_script (void);

void open_info (gpointer data, guint edit, GtkWidget *widget);

void do_new_script (gpointer data, guint action, GtkWidget *widget);

void do_open_csv_box (char *fname, int code, int append);

int do_store (char *mydatfile, gretlopt oflag, int overwrite);

void view_latex (PRN *prn);

void save_latex (PRN *prn, const char *fname);

void do_save_text (char *fname, MODEL *pmod);

int execute_script (const char *runfile, const char *buf,
		    PRN *prn, int exec_code);

int gui_exec_line (char *line, 
		   LOOPSET **plp, int *plstack, int *plrun, 
		   PRN *prn, int exec_code, 
		   const char *myname); 

int check_and_record_command (void);

/* other */

int latex_compile (char *texshort);

void add_mahalanobis_data (windata_t *vwin);

void add_pca_data (windata_t *vwin);

void add_fcast_data (windata_t *vwin);

#endif /* LIBRARY_H */

