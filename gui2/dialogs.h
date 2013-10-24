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

#ifndef DIALOGS_H
#define DIALOGS_H

#include <stdarg.h>

#define GRETL_CANCEL (-1)
#define canceled(r) (r == -1)

enum {
    GRETL_YES,
    GRETL_NO
} buttons;

#ifndef BUILDING_PLUGIN
typedef struct dialog_t_ dialog_t;

void copy_format_dialog (windata_t *vwin, int action);
#endif

void errbox (const char *template, ...);

void infobox (const char *template, ...);

void warnbox (const char *template, ...);

void maybe_warn (void);

void file_read_errbox (const char *fname);

void file_write_errbox (const char *fname);

gint yes_no_dialog (const char *title, const char *msg, int cancel);

gint no_yes_dialog (const char *title, const char *msg);

gint yes_no_dialog_with_parent (const char *title, const char *msg, 
				int cancel, GtkWidget *parent);

void gretl_dialog_add_message (GtkWidget *dlg, const char *msg);

int make_default_storelist (void);

gboolean exit_check (void);

void menu_exit_check (void);

double gui_double_from_string (const char *str, int *err);

int csv_options_dialog (int ci, GretlObjType otype, GtkWidget *parent);

void rand_seed_dialog (void);

int select_list_dialog (char *listname);

void database_description_dialog (const char *binname);

int select_var_from_list_with_opt (const int *list, 
				   const char *query,
				   dialog_opts *opts,
				   int hcode,
				   GtkWidget *parent);

int select_var_from_list (const int *list, const char *query,
			  GtkWidget *parent);

void sample_range_dialog (GtkAction *action, gpointer p);

void range_dummy_dialog (GtkAction *action, gpointer p);

void sample_restrict_dialog (GtkAction *action, gpointer p);

int panel_graph_dialog (int *t1, int *t2);

void data_compact_dialog (int spd, int *target_pd, int *mon_start,
			  CompactMethod *method, int *repday, 
			  GtkWidget *parent);

void data_expand_dialog (int spd, int *interpol, GtkWidget *parent);

int pergm_dialog (gretlopt *opt, int *spinval, int spinmin, int spinmax,
		  GtkWidget *parent);

int density_dialog (int vnum, double *bw);

int radio_dialog (const char *title, const char *label, const char **opts, 
		  int nopts, int deflt, int helpcode, GtkWidget *parent);

int radio_dialog_with_spinner (const char *title, const char **opts, 
			       int nopts, int deflt, int helpcode,
			       int *spinvar, const char *spintxt,
			       int spinmin, int spinmax,
			       GtkWidget *parent);

int radio_dialog_with_check (const char *title, const char *label, 
			     const char **opts, int nopts, int deflt, int hcode,
			     int *checkvar, const char *checktxt,
			     GtkWidget *parent);

void set_checks_dialog_extra (int i, GtkWidget *extra);

GtkWidget *
build_checks_dialog (const char *title, const char *blurb,
		     const char **opts, 
		     int nopts, int *active, 
		     int check_min, int check_max,
		     int nradios, int *rvar, 
		     int *spinvar, const char *spintxt, 
		     int spinmin, int spinmax, 
		     int hcode, GtkWidget *parent, int *ret);

int checks_dialog (const char *title, const char *blurb,
		   const char **opts, 
		   int nopts, int *active, 
		   int check_min, int check_max,
		   int nradios, int *rvar,
		   int *spinvar, const char *spintxt, 
		   int spinmin, int spinmax, 
		   int hcode, GtkWidget *parent);

int checks_only_dialog (const char *title, const char *blurb,
			const char **opts, int nopts,
			int *active, int hcode,
			GtkWidget *parent);

int spin_dialog (const char *title, const char *blurb,
		 int *spinvar, const char *spintxt, 
		 int spinmin, int spinmax, int helpcode,
		 GtkWidget *parent);

int combo_selector_dialog (GList *list, const char *msg,
			   int deflt, GtkWidget *parent);

int yes_no_help_dialog (const char *msg, int hcode);

int add_obs_dialog (const char *blurb, int addmin, GtkWidget *parent);

int forecast_dialog (int t1min, int t1max, int *t1, 
		     int t2min, int t2max, int *t2,
		     int *k, int pmin, int pmax, int *p,
		     int dyn, gretlopt *optp,
		     double *conf, MODEL *pmod,
		     GtkWidget *parent);

void dialog_add_confidence_selector (GtkWidget *dlg, double *conf,
				     gretlopt *gopt);

int freq_dialog (const char *title, const char *blurb,
		 int *nbins, int nbmax, double *f0, double *fwid,
		 double xmin, double xmax, int *dist, int *plot);

int model_table_dialog (int *colhead_opt, int *se_opt, int *pv_opt,
			int *ast_opt, int *figs, char *fmt,
			GtkWidget *parent);

int bootstrap_dialog (windata_t *vwin, int *pp, int *pB,
		      gretlopt *popt);

int chow_dialog (int tmin, int tmax, int *t, int *dumv,
		 GtkWidget *parent);

int iter_control_dialog (int *optim, int *pmaxit, double *ptol, 
			 int *plmem, GtkWidget *parent);

void tex_format_dialog (GtkAction *action, gpointer data);

int paste_data_dialog (int *append);

int object_name_entry_dialog (char *name, GretlType type,
			      const char *labeltxt, int *show,
			      GtkWidget *parent);

int hc_config_dialog (char *vname, gretlopt opt, gboolean robust_conf, 
		      GtkWidget *parent);

#endif /* DIALOGS_H */
