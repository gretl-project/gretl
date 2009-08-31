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

enum {
    GRETL_YES,
    GRETL_NO,
    HELP_BUTTON
} buttons;

#ifndef BUILDING_PLUGIN
typedef struct dialog_t_ dialog_t;

void copy_format_dialog (windata_t *vwin, int action);
#endif

void errbox (const char *template, ...);

void infobox (const char *template, ...);

void warnbox (const char *template, ...);

void file_read_errbox (const char *fname);

void file_write_errbox (const char *fname);

gint yes_no_dialog (const char *title, const char *msg, int cancel);

int make_default_storelist (void);

gboolean exit_check (void);

void menu_exit_check (void);

double gui_double_from_string (const char *str, int *err);

int csv_options_dialog (gretlopt *optp);

void rand_seed_dialog (void);

void database_description_dialog (const char *binname);

int select_var_from_list_with_opt (const int *list, 
				   const char *query,
				   dialog_opts *opts,
				   int hcode);

int select_var_from_list (const int *list, const char *query);

void sample_range_dialog (GtkAction *action, gpointer p);

void panel_structure_dialog (DATAINFO *pdinfo);

void data_compact_dialog (GtkWidget *w, int spd, int *target_pd, 
			  int *mon_start, CompactMethod *method,
			  int *repday);

void data_expand_dialog (GtkWidget *w, int spd, int *target_pd);

int density_dialog (int vnum, double *bw);

int radio_dialog (const char *title, const char *label, const char **opts, 
		  int nopts, int deflt, int helpcode);

int radio_dialog_with_spinner (const char *title, const char **opts, 
			       int nopts, int deflt, int helpcode,
			       int *spinvar, const char *spintxt,
			       int spinmin, int spinmax);

int radio_dialog_with_check (const char *title, const char *label, 
			     const char **opts, int nopts, int deflt, int hcode,
			     int *checkvar, const char *checktxt);

GtkWidget *
build_checks_dialog (const char *title, const char *blurb,
		     const char **opts, int nopts,
		     int *active, int nradios, int *rvar, int *spinvar, 
		     const char *spintxt, int spinmin, int spinmax, 
		     int hcode, int *ret);

int checks_dialog (const char *title, const char *blurb,
		   const char **opts, int nopts, 
		   int *active, int nradios, int *rvar,
		   int *spinvar, const char *spintxt, 
		   int spinmin, int spinmax, int helpcode);

int checks_only_dialog (const char *title, const char *blurb,
			const char **opts, int nopts,
			int *active, int hcode);

int spin_dialog (const char *title, const char *blurb,
		 int *spinvar, const char *spintxt, 
		 int spinmin, int spinmax, int helpcode);

int get_obs_dialog (const char *title, const char *text,
		    const char *t1str, const char *t2str,
		    int t1min, int t1max, int *t1, 
		    int t2min, int t2max, int *t2);

int add_obs_dialog (const char *blurb, int addmin);

int forecast_dialog (int t1min, int t1max, int *t1, 
		     int t2min, int t2max, int *t2,
		     int *k, int pmin, int pmax, int *p,
		     int dyn, gretlopt *optp,
		     double *conf, MODEL *pmod);

void dialog_add_confidence_selector (GtkWidget *dlg, double *conf);

int freq_dialog (const char *title, const char *blurb,
		 int *nbins, int nbmax, double *f0, double *fwid,
		 double xmin, double xmax, int *dist, int plot);

void bootstrap_dialog (windata_t *vwin, int *pp, int *pB,
		       gretlopt *popt, int *cancelled);

void lmax_dialog (double *lmax, double ymax);

int chow_dialog (int tmin, int tmax, int *t, int *dumv);

void tex_format_dialog (void);

int gretl_all_done (void);

#endif /* DIALOGS_H */
