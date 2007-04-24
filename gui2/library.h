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

#ifdef OSX_BUILD
int osx_open_file (const char *path);

int osx_open_url (const char *url);
#endif

typedef struct _selector selector;

void library_command_init (void);

void library_command_free (void);

CMD *get_lib_cmd (void);

int gretl_command_sprintf (const char *template, ...);

int gretl_command_strcpy (const char *s);

int *command_list_from_string (char *s);

int user_fopen (const char *fname, char *fullname, PRN **pprn);

gint bufopen (PRN **pprn);

void do_menu_op (gpointer data, guint action, GtkWidget *w);

void do_run_script (GtkWidget *w, gpointer data);

void gui_errmsg (const int errcode);

void register_graph (void);

/* sample-related functions */

int bool_subsample (gretlopt opt);

void do_samplebool (GtkWidget *w, dialog_t *dlg);

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

int do_coint (selector *sr);

void do_forecast (gpointer p, guint u, GtkWidget *w);

void do_bootstrap (gpointer p, guint u, GtkWidget *w) ;

int do_coeff_sum (selector *sr);

int do_add_omit (selector *sr);

int do_VAR_omit (selector *sr);

int do_confidence_region (selector *sr);

void do_lmtest (gpointer data, guint aux_code, GtkWidget *w);

void do_autocorr (gpointer data, guint u, GtkWidget *w);

void do_chow_cusum (gpointer data, guint action, GtkWidget *w);

void do_reset (gpointer data, guint u, GtkWidget *w);

void unit_root_test (gpointer data, guint u, GtkWidget *w);

void do_arch (gpointer data, guint u, GtkWidget *w);

void do_restrict (GtkWidget *w, dialog_t *dlg);

void do_nls_model (GtkWidget *w, dialog_t *dlg);

void do_mle_model (GtkWidget *w, dialog_t *dlg);

void do_gmm_model (GtkWidget *w, dialog_t *dlg);

void do_eqn_system (GtkWidget *w, dialog_t *dlg);

void do_saved_eqn_system (GtkWidget *w, dialog_t *dlg);

int do_model (selector *sr);

int do_vector_model (selector *sr);

void do_graph_model (const int *list, int fit);

void do_gini (gpointer data, guint u, GtkWidget *w);

void do_kernel (gpointer data, guint u, GtkWidget *w);

void do_vif (gpointer data, guint u, GtkWidget *w);

void do_leverage (gpointer data, guint u, GtkWidget *w);

void add_leverage_data (windata_t *vwin);

void do_coeff_intervals (gpointer data, guint i, GtkWidget *w);

void do_panel_diagnostics (gpointer data, guint u, GtkWidget *w);

int do_spearman (selector *sr);

int out_of_sample_info (int add_ok, int *t2);

/* variable-related functions */

void do_minibuf (GtkWidget *w, dialog_t *dlg);

void do_genr (GtkWidget *w, dialog_t *dlg);

void do_model_genr (GtkWidget *w, dialog_t *dlg);

void add_rand_series (GtkWidget *w, dialog_t *dlg);

void do_global_setmiss (GtkWidget *w, dialog_t *dlg);

void do_variable_setmiss (GtkWidget *w, dialog_t *dlg);

void do_edit_label (GtkWidget *w, dialog_t *dlg);

int do_rename_variable (int v, const char *newname, int full);

int record_varlabel_change (int v);

void do_resid_freq (gpointer p, guint action, GtkWidget *w);

void do_freqplot (gpointer p, guint gamma, GtkWidget *w);

void do_corrgm (gpointer p, guint u, GtkWidget *w);

void residual_correlogram (gpointer p, guint u, GtkWidget *w);

void do_pergm (gpointer p, guint opt, GtkWidget *w);

void residual_periodogram (gpointer p, guint opt, GtkWidget *w);

#if defined (HAVE_TRAMO) || defined (HAVE_X12A)
void do_tramo_x12a (gpointer p, guint opt, GtkWidget *w);
#endif

void do_range_mean (gpointer p, guint opt, GtkWidget *w);

void do_hurst (gpointer p, guint opt, GtkWidget *w);

void do_outcovmx (gpointer p, guint u, GtkWidget *w);

void do_anova (gpointer p, guint u, GtkWidget *w);

void add_dummies (gpointer p, guint action, GtkWidget *w);

void add_index (gpointer p, guint tm, GtkWidget *w);

void do_add_obs (gpointer p, guint u, GtkWidget *w);

void do_remove_obs (gpointer p, guint u, GtkWidget *w);

void add_logs_etc (gpointer p, guint action, GtkWidget *w);

int add_system_resid (gpointer p, int eqnum, int ci);

int add_fit_resid (MODEL *pmod, int code, int undo);

int add_system_resid (gpointer p, int eqnum, int ci);

void add_model_stat (MODEL *pmod, int which);

void resid_plot (gpointer p, guint xvar, GtkWidget *w);

void fit_actual_plot (gpointer p, guint xvar, GtkWidget *w);

void fit_actual_splot (gpointer p, guint u, GtkWidget *w);

void display_fit_resid (gpointer p, guint code, GtkWidget *w);

void do_graph_var (int varnum);

void do_boxplot_var (int varnum);

void ts_plot_var (gpointer p, guint opt, GtkWidget *w);

int do_scatters (selector *sr);

int do_graph_from_selector (selector *sr);

int do_splot_from_selector (selector *sr);

void plot_from_selection (gpointer p, guint action, GtkWidget *w);

void do_box_graph (GtkWidget *w, dialog_t *dlg);

int do_dummy_graph (selector *sr);

int do_xcorrgm (selector *sr);

void delete_selected_vars (void);

void delete_single_var (int id);

void display_selected (gpointer p, guint action, GtkWidget *w);

void display_var (void);

/* script- and file-related functions */

void do_open_script (void);

void open_info (gpointer p, guint edit, GtkWidget *w);

void do_new_script (gpointer p, guint action, GtkWidget *w);

void do_open_csv_box (char *fname, int code, int append);

int do_store (char *mydatfile, gretlopt oflag);

void view_latex (PRN *prn);

void save_latex (PRN *prn, const char *fname);

void do_save_text (char *fname, MODEL *pmod);

int gui_exec_line (ExecState *s, double ***pZ, DATAINFO **ppdinfo);

int check_and_record_command (void);

/* other */

int latex_compile (char *texshort);

void add_mahalanobis_data (windata_t *vwin);

void add_pca_data (windata_t *vwin);

void add_fcast_data (windata_t *vwin);

void maybe_display_string_table (void);

#endif /* LIBRARY_H */

