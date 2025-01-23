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

#ifndef GRETL_PANEL_H_
#define GRETL_PANEL_H_

int panel_diagnostics (MODEL *pmod, DATASET *dset,
		       gretlopt opt, PRN *prn);

MODEL real_panel_model (const int *list, DATASET *dset,
			gretlopt opt, PRN *prn);

MODEL panel_wls_by_unit (const int *list, DATASET *dset,
			 gretlopt opt, PRN *prn);

int panel_autocorr_test (MODEL *pmod, DATASET *dset,
			 gretlopt opt, PRN *prn);

int panel_xdepend_test (MODEL *pmod, DATASET *dset,
			gretlopt opt, PRN *prn);

int groupwise_hetero_test (MODEL *pmod, DATASET *dset,
			   gretlopt opt, PRN *prn);

int panel_DW_pval_ok (const MODEL *pmod);

double BFN_panel_DW_pvalue (MODEL *pmod, const DATASET *dset, int *err);

int panel_tsls_robust_vcv (MODEL *pmod, DATASET *dset);

int set_panel_structure_from_vars (int uv, int tv, DATASET *dset);

int set_panel_structure_from_varnames (const char *uname,
				       const char *tname,
				       DATASET *dset);

int set_panel_group_strings (const char *vname,
			     const char *grpnames,
			     DATASET *dset);

int switch_panel_orientation (DATASET *dset);

int balanced_panel (const DATASET *dset);

int undo_panel_padding (DATASET *dset);

int *panel_list_omit (const MODEL *orig, const int *drop, int *err);

int *panel_list_add (const MODEL *orig, const int *add, int *err);

int panel_variance_info (const double *x, const DATASET *dset,
			 double xbar, double *psw, double *psb);

int plausible_panel_time_var (const DATASET *dset);

int is_panel_time_var (const DATASET *dset, int v,
                       int tnax, int *minval,
                       int *incr);

int is_panel_unit_var (const DATASET *dset, int v, int tmax);

int panel_isconst (int t1, int t2, int pd, const double *x,
		   int bygroup);

int series_is_group_invariant (const DATASET *dset, int v);

int panel_padding_rows (const DATASET *dset);

int time_series_from_panel (DATASET *tset, const DATASET *pset);

int obs_index_from_aqm (const char *s, int pd, int t0, int n);

int usable_group_names_series_id (const DATASET *dset,
				  int maxlen);

#endif /* GRETL_PANEL_H_ */
