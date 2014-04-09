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

/* library.h for gretl */

#ifndef LIBRARY_H
#define LIBRARY_H

/* general purpose functions */

#ifdef G_OS_WIN32
void win_show_error (DWORD dw);

int winfork (char *cmdline, const char *dir, int wshow,
	     DWORD flags);
#endif

#ifdef OS_OSX
int osx_open_file (const char *path);

int osx_open_url (const char *url);
#endif

typedef struct _selector selector;

void library_command_init (void);

void library_command_free (void);

CMD *get_lib_cmd (void);

int lib_command_sprintf (const char *template, ...);

int lib_command_strcpy (const char *s);

int lib_command_strcat (const char *s);

int record_command_verbatim (void);

int parse_lib_command (void);

int record_lib_command (void);

int execute_script (const char *runfile, const char *buf,
		    PRN *prn, int exec_code);

int user_fopen (const char *fname, char *fullname, PRN **pprn);

gint bufopen (PRN **pprn);

void menu_op_action (GtkAction *action, gpointer data);

void do_menu_op (int ci, const char *liststr, gretlopt opt);

void do_run_script (GtkWidget *w, windata_t *vwin);

void run_script_fragment (windata_t *vwin, gchar *buf);

void gui_errmsg (int errcode);

void gui_warnmsg (int errcode);

void errmsg_plus (int err, const char *plus);

void gui_graph_handler (int err);

/* sample-related functions */

int bool_subsample (gretlopt opt);

int do_set_sample (void);

void drop_all_missing (void);

void count_missing (void);

void do_add_markers (const char *fname);

int do_save_markers (const char *fname);

void markers_callback (void);

void do_add_labels (const char *fname);

int do_save_labels (const char *fname);

void labels_callback (void);

int dataset_is_restricted (void);

int maybe_restore_full_data (int action);

void gui_transpose_data (void);

void gui_sort_data (void);

void gui_resample_data (void);

int dataset_is_subsampled (void);

void set_original_n (int n);

int get_original_n (void);

/* model-related functions */

int do_coint (selector *sr);

void gui_do_forecast (GtkAction *action, gpointer p);

void do_bootstrap (GtkAction *action, gpointer p);

int do_coeff_sum (selector *sr);

int do_add_omit (selector *sr);

int do_VAR_omit (selector *sr);

void VAR_omit_auto (GtkAction *action, gpointer p);

int do_confidence_region (selector *sr);

void do_modtest (GtkAction *action, gpointer p);

void do_autocorr (GtkAction *action, gpointer p);

void do_dwpval (GtkAction *action, gpointer p);

void do_chow_cusum (GtkAction *action, gpointer p);

void do_reset (GtkAction *action, gpointer p);

void unit_root_test (int ci);

void ur_callback (GtkAction *action);

void do_arch (GtkAction *action, gpointer p);

void do_restrict (GtkWidget *w, dialog_t *dlg);

void do_nls_model (GtkWidget *w, dialog_t *dlg);

void do_mle_model (GtkWidget *w, dialog_t *dlg);

void do_gmm_model (GtkWidget *w, dialog_t *dlg);

void do_eqn_system (GtkWidget *w, dialog_t *dlg);

void do_saved_eqn_system (GtkWidget *w, dialog_t *dlg);

int do_model (selector *sr);

int do_nonparam_model (selector *sr);

int do_vector_model (selector *sr);

void do_graph_model (const int *list, int fit);

void do_nonparam_plot (windata_t *vwin);

void do_gini (void);

void do_qqplot (void);

void do_kernel (void);

void do_vif (GtkAction *action, gpointer p);

void do_leverage (GtkAction *action, gpointer p);

void add_leverage_data (windata_t *vwin);

void do_coeff_intervals (GtkAction *action, gpointer p);

void do_panel_tests (GtkAction *action, gpointer p);

int do_rankcorr (selector *sr);

int out_of_sample_info (int add_ok, int *t2);

/* variable-related functions */

void do_minibuf (GtkWidget *w, dialog_t *dlg);

void do_genr (GtkWidget *w, dialog_t *dlg);

void do_model_genr (GtkWidget *w, dialog_t *dlg);

void do_selector_genr (GtkWidget *w, dialog_t *dlg);

void do_fncall_genr (GtkWidget *w, dialog_t *dlg);

void do_range_dummy_genr (const gchar *buf);

void do_global_setmiss (GtkWidget *w, dialog_t *dlg);

void do_variable_setmiss (GtkWidget *w, dialog_t *dlg);

void do_edit_label (GtkWidget *w, dialog_t *dlg);

int do_rename_variable (int v, const char *newname);

int record_varlabel_change (int v);

void do_resid_freq (GtkAction *action, gpointer p);

void do_freq_dist (void);

void do_corrgm (void);

void residual_correlogram_callback (GtkAction *action, gpointer p);

void do_pergm (GtkAction *action);

void do_fractint (GtkAction *action);

void residual_periodogram_callback (GtkAction *action, gpointer p);

void residual_qq_plot (GtkAction *action, gpointer p);

#if defined (HAVE_TRAMO) || defined (HAVE_X12A)
void do_tramo_x12a (GtkAction *action, gpointer p);
#endif

void do_range_mean (void);

void do_hurst (void);

void do_outcovmx (GtkAction *action, gpointer p);

void do_anova (GtkAction *action, gpointer p);

void add_dummies (GtkAction *action);

void add_index (GtkAction *action);

void do_add_obs (void);

void do_remove_obs (void);

void add_logs_etc (int ci, int varnum);

void logs_etc_callback (GtkAction *action);

void add_system_resid (GtkAction *action, gpointer p);

int save_fit_resid (windata_t *vwin, int code);

int save_bundled_series (const double *x, const char *key,
			 const char *note, windata_t *vwin);

void add_model_stat (MODEL *pmod, int which, windata_t *vwin);

void resid_plot (GtkAction *action, gpointer p);

void fit_actual_plot (GtkAction *action, gpointer p);

void fit_actual_splot (GtkAction *action, gpointer p);

void display_fit_resid (GtkAction *action, gpointer p);

void do_graph_var (int varnum);

void do_boxplot_var (int varnum, gretlopt opt);

int do_factorized_boxplot (selector *sr);

void ts_plot_callback (void);

int do_scatters (selector *sr);

int do_graph_from_selector (selector *sr);

int do_splot_from_selector (selector *sr);

void plot_from_selection (int code);

void do_box_graph (GtkWidget *w, dialog_t *dlg);

int do_dummy_graph (selector *sr);

int do_xyz_graph (selector *sr);

int do_xcorrgm (selector *sr);

void delete_selected_vars (void);

void delete_single_var (int id);

void display_selected (void);

void display_var (void);

/* script- and file-related functions */

gboolean do_open_script (int action);

void dataset_info (void);

void do_new_script (int code);

void new_script_callback (GtkAction *action);

int do_store (char *filename, int action, gpointer data);

void set_csv_exclude_obs (gboolean s);

gboolean get_csv_exclude_obs (void);

void do_save_text (char *fname, MODEL *pmod);

int gui_exec_line (ExecState *s, DATASET *dset);

void start_wait_for_output (GtkWidget *w, gboolean big);

void stop_wait_for_output (GtkWidget *w);

int waiting_for_output (void);

int output_flush_in_progress (void);

/* other */

gchar *get_genr_string (GtkWidget *entry, dialog_t *dlg);

int menu_op_wrapper (selector *sr);

int max_untouchable_series_ID (void);

void add_mahalanobis_data (windata_t *vwin);

void add_pca_data (windata_t *vwin);

void add_fcast_data (windata_t *vwin);

void add_nonparam_data (windata_t *vwin);

void VECM_add_EC_data (GtkAction *action, gpointer p);

void maybe_display_string_table (void);

#endif /* LIBRARY_H */

