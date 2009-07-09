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

/* selector.c for gretl */

#include "gretl.h"
#include "selector.h"
#include "dlgutils.h"
#include "menustate.h"
#include "fileselect.h"
#include "fnsave.h"
#include "treeutils.h"
#include "lagpref.h"

#include "var.h"
#include "gretl_func.h"
#include "libset.h"
#include "johansen.h"

#if 0
# include "bootstrap.h"
#endif

#define LDEBUG 0
#define VLDEBUG 0

enum {
    SR_LVARS  = 1,
    SR_RVARS1,
    SR_RVARS2
};

#define N_EXTRA  8
#define N_RADIOS 3

struct _selector {
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *action_area;
    GtkWidget *lvars;
    GtkWidget *depvar;
    GtkWidget *rvars1;
    GtkWidget *rvars2;
    GtkWidget *default_check;
    GtkWidget *add_button;
    GtkWidget *remove_button;
    GtkWidget *lags_button;
    GtkWidget *hess_button;
    GtkWidget *x12a_button;
    GtkWidget *extra[N_EXTRA];
    GtkWidget *radios[N_RADIOS];
    int ci;
    int active_var;
    int error;
    int n_left;
    gretlopt opts;
    char *cmdlist;
    gpointer data;
    int (*callback)();
};

#define EXTRA_LAGS (N_EXTRA - 1)

/* single-equation estimation commands plus some GUI extensions */
#define MODEL_CODE(c) (MODEL_COMMAND(c) || c == CORC || c == HILU || \
                       c == PWE || c == PANEL_WLS || c == PANEL_B || \
                       c == OLOGIT || c == OPROBIT || c == MLOGIT || \
	               c == IV_LIML || c == IV_GMM)

#define IV_MODEL(c) (c == IVREG || c == IV_LIML || c == IV_GMM)

#define COINT_CODE(c) (c == COINT || c == COINT2)

#define VEC_CODE(c) (c == COINT || c == COINT2 || c == VAR || \
                     c == VECM || c == VLAGSEL)

#define VEC_MODEL_CODE(c) (c == VAR || c == VECM || c == VLAGSEL)

#define VECLAGS_CODE(c) (c == VAR || c == VECM)

#define ADDVAR_CODE(c) (c == LOGS || c == LAGS || c == SQUARE || \
                        c == DIFF || c == LDIFF)

#define GRAPH_CODE(c) (c == GR_PLOT || c == GR_XY || c == GR_IMP || \
		       c == GR_DUMMY || c == GR_XYZ || c == GR_3D)

#define TWO_VARS_CODE(c) (c == SPEARMAN || c == ELLIPSE || c == XCORRGM)

#define THREE_VARS_CODE(c) (c == GR_DUMMY || c == GR_XYZ || \
			    c == GR_3D || c == ANOVA)

#define WANT_TOGGLES(c) (c == ARBOND || \
                         c == ARMA || \
                         c == COINT || \
                         c == COINT2 || \
			 c == CORR || \
                         c == GARCH || \
                         c == HECKIT || \
                         c == HILU || \
                         c == INTREG || \
                         c == IVREG || \
                         c == LOGIT || \
                         c == OLOGIT || \
                         c == MLOGIT || \
                         c == MPOLS || \
                         c == OLS || \
                         c == PANEL || \
                         c == PANEL_WLS || \
                         c == PANEL_B || \
                         c == PROBIT || \
                         c == OPROBIT || \
	                 c == QUANTREG || \
			 c == SPEARMAN || \
                         c == TOBIT || \
                         c == VAR || \
                         c == VECM || \
                         c == VLAGSEL || \
                         c == WLS || \
                         c == XTAB)

#define USE_VECXLIST(c) (c == VAR || c == VLAGSEL || c == VECM || \
                         c == COINT2)

#define USE_RXLIST(c) (c == VECM || c == COINT2)

#define AUX_LAST(c) (c == IVREG || \
                     c == IV_LIML || \
                     c == IV_GMM || \
                     c == HECKIT || \
                     c == VAR || \
                     c == VLAGSEL || \
                     c == VECM || \
                     c == COINT2)

#define USE_ZLIST(c) (c == IVREG || c == IV_LIML || c == IV_GMM || c == HECKIT)

#define RHS_PREFILL(c) (c == CORR || \
	                c == MAHAL || \
			c == PCA || \
                	c == SUMMARY || \
			c == XTAB)

#define dataset_lags_ok(d) ((d)->structure == TIME_SERIES || \
			    (d)->structure == SPECIAL_TIME_SERIES || \
                            (d)->structure == STACKED_TIME_SERIES)

#define select_lags_primary(c) (MODEL_CODE(c))

#define select_lags_depvar(c) (MODEL_CODE(c) && c != ARMA && c != ARBOND) 

/* Should we have a lags button associated with auxiliary
   variable selector? */

#define select_lags_aux(c) (c == VAR || c == VLAGSEL || c == VECM || \
                            c == IVREG || c == IV_LIML || c == IV_GMM || \
                            c == HECKIT)

#define list_lag_special(i) (i < -1)

/* static state variables */

static int default_y = -1;
static int want_seasonals = 0;
static int default_order;
static int vartrend = 0;
static int varconst = 1;
static int lags_hidden;
static int arma_p = 1;
static int arma_P = 0;
static int arima_d = 0;
static int arima_D = 0;
static int arma_q = 1;
static int arma_Q = 0;
static int garch_p = 1;
static int garch_q = 1;
static int arma_const = 1;
static int arma_x12 = 0;
static int arma_hessian = 1;
static int selvar;
static int offvar;
static int jrank = 1;
static int jcase = J_UNREST_CONST;
static int verbose;
static int lovar;
static int hivar;
static int wtvar;

static int y_x_lags_enabled;
static int y_w_lags_enabled;

static int *xlist;
static int *instlist;
static int *veclist;
static int *vecxlist;

static char *arlags;
static char *malags;

static gretlopt model_opt;

static GtkWidget *multiplot_label;
static GtkWidget *multiplot_menu;

static selector *open_selector;

static gint listvar_special_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data);
static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event, 
			       selector *sr);
static gint listvar_flagcol_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data);
static int list_show_var (int v, int ci, int show_lags);
static int functions_list (selector *sr);
static void primary_rhs_varlist (selector *sr, GtkWidget *vbox);
static gboolean lags_dialog_driver (GtkWidget *w, selector *sr);

#define spinner_get_int(b) (gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(b)))

static int want_combo (selector *sr)
{
    return (sr->ci == ARMA ||
	    sr->ci == VECM || 
	    sr->ci == COINT ||
	    sr->ci == COINT2 ||
	    sr->ci == IV_GMM);
}

static int want_radios (selector *sr)
{
    int c = sr->ci;

    if (c == PANEL || c == SCATTERS || c == ARBOND || 
	c == LOGIT || c == PROBIT || c == HECKIT ||
	c == XTAB || c == SPEARMAN || c == PCA ||
	c == QUANTREG) {
	return 1;
    } else if (c == OMIT) {
	windata_t *vwin = (windata_t *) sr->data;
	MODEL *pmod = (MODEL *) vwin->data;

	if (pmod->ci == HECKIT) {
	    /* omit: Wald test only */
	    sr->opts |= OPT_W;
	} else {
	    return 1;
	}
    }

    return 0;
}

static int selection_at_max (selector *sr, int nsel)
{
    int ret = 0;

    if (TWO_VARS_CODE(sr->ci) && nsel == 2) {
	ret = 1;
    }

    return ret;
}

static int sr_get_lag_context (selector *sr, int locus)
{
    int c = 0;

    if (sr == NULL || !dataset_lags_ok(datainfo)) {
	return 0;
    }

    if (locus == SR_RVARS1 && select_lags_primary(sr->ci)) {
	c = LAG_X;
    } else if (locus == SR_RVARS2 && select_lags_aux(sr->ci)) {
	c = (USE_ZLIST(sr->ci))? LAG_W : LAG_X;
    }

    return c;
}

static int lag_context_from_widget (GtkWidget *w)
{
    gint locus = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "locus"));
    selector *sr = g_object_get_data(G_OBJECT(w), "sr");
    int context = 0;

    if (sr != NULL && locus) {
	context = sr_get_lag_context(sr, locus);
    }

    return context;
}

static int lags_button_relevant (selector *sr, int locus)
{
    if (sr->lags_button != NULL) {
	if (locus == SR_RVARS1 && select_lags_primary(sr->ci)) {
	    return 1;
	} else if (locus == SR_RVARS2 && select_lags_aux(sr->ci)) {
	    return 1;
	}
    }

    return 0;
}

void enable_lags_for_context (int context, gboolean s)
{
    if (context == LAG_Y_X) {
	y_x_lags_enabled = s;
    } else if (context == LAG_Y_W) {
	y_w_lags_enabled = s;
    }
}

gboolean lags_enabled_for_context (int context)
{
    if (context == LAG_Y_X) {
	return y_x_lags_enabled;
    } else if (context == LAG_Y_W) {
	return y_w_lags_enabled;
    } else {
	return TRUE;
    }
}

void clear_selector (void)
{
    default_y = -1;
    default_order = 0;
    selvar = 0;
    offvar = 0;
    wtvar = 0;
    vartrend = 0;
    varconst = 1;
    lovar = hivar = 0;

    arma_p = 1;
    arima_d = 0;
    arma_q = 1;
    arma_P = 0;
    arima_D = 0;
    arma_Q = 0;
    arma_const = 1;
    arma_hessian = 1;

    garch_p = 1;
    garch_q = 1;

    jrank = 1;
    jcase = J_UNREST_CONST;

    free(xlist);
    xlist = NULL;
    free(instlist);
    instlist = NULL;

    free(veclist);
    veclist = NULL;
    free(vecxlist);
    vecxlist = NULL;

    free(arlags);
    arlags = NULL;
    free(malags);
    malags = NULL;

    destroy_lag_preferences();
}

static int presel;

void selector_set_varnum (int v)
{
    presel = v;
}

void modelspec_dialog (int ci)
{
    selection_dialog(_("gretl: specify model"), do_model, ci);
}

static int varnum_from_keystring (MODEL *pmod, const char *key)
{
    char *s = (char *) gretl_model_get_data(pmod, key);
    int v = 0;

    if (s != NULL && *s != '\0') {
	v = series_index(datainfo, s);
	if (v == datainfo->v) {
	    v = 0;
	}
    }

    return v;
}

static int *lag_list_from_mask (char *mask, int k)
{
    int *list;
    int i, j;

    list = gretl_list_new(k);

    if (list != NULL) {
	j = 1;
	for (i=0; mask[i]; i++) {
	    if (mask[i] == '1') {
		list[j++] = i + 1;
	    } else {
		list[0] -= 1;
	    }
	}
    }

    return list;
}

static void retrieve_arma_info (MODEL *pmod)
{
    int acode = gretl_model_get_int(pmod, "arma_flags");
    int *laglist;
    char *mask;

    if (acode & ARMA_X12A) {
	model_opt |= OPT_X;
    }

    if (!(acode & ARMA_EXACT)) {
	model_opt |= OPT_C;
    } 

    if (pmod->opt & OPT_G) {
	/* use OPG for covariance matrix */
	arma_hessian = 0;
    }

    arma_p = arma_model_nonseasonal_AR_order(pmod);
    arma_q = arma_model_nonseasonal_MA_order(pmod);
    arima_d = gretl_model_get_int(pmod, "arima_d");
    arma_P = gretl_model_get_int(pmod, "arma_P");
    arma_Q = gretl_model_get_int(pmod, "arma_Q");
    arima_D = gretl_model_get_int(pmod, "arima_D");
    arma_const = pmod->ifc;

    free(arlags);
    arlags = NULL;

    if (arma_p > 0) {
	mask = (char *) gretl_model_get_data(pmod, "pmask");
	if (mask != NULL) {
	    laglist = lag_list_from_mask(mask, arma_p);
	    if (laglist != NULL) {
		arlags = gretl_list_to_string(laglist);
		free(laglist);
	    }
	}
    }

    free(malags);
    malags = NULL;

    if (arma_q > 0) {
	mask = (char *) gretl_model_get_data(pmod, "qmask");
	if (mask != NULL) {
	    laglist = lag_list_from_mask(mask, arma_q);
	    if (laglist != NULL) {
		free(malags);
		malags = gretl_list_to_string(laglist);
		free(laglist);
	    }
	} 
    }   
}

static void retrieve_AR_lags_info (MODEL *pmod)
{
    free(arlags);
    arlags = NULL;

    if (pmod->arinfo != NULL && pmod->arinfo->arlist != NULL) {
	arlags = gretl_list_to_string(pmod->arinfo->arlist);
    }
}

static void retrieve_heckit_info (MODEL *pmod, int *gotinst)
{
    int *zlist = gretl_model_get_secondary_list(pmod);

    if (zlist != NULL) {
	selvar = zlist[1];
	free(instlist);
	instlist = gretl_list_copy_from_pos(zlist, 2);
	if (instlist != NULL) {
	    *gotinst = 1;
	}
	free(zlist);
    }

    if (pmod->opt & OPT_T) {
	/* two-step estimation */
	model_opt |= OPT_T;
    }    
}

void selector_from_model (void *ptr, int ci)
{
    model_opt = OPT_NONE;

    if (ci == VIEW_MODEL) {
	/* single-equation model */
	MODEL *pmod = (MODEL *) ptr;
	int sel_ci = pmod->ci;
	int dv = -1, gotinst = 0;

	if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM) {
	    revise_nl_model(pmod);
	    return;
	}

	if (pmod->ci == INTREG) {
	    lovar = varnum_from_keystring(pmod, "lovar");
	    hivar = varnum_from_keystring(pmod, "hivar");
	} else {
	    dv = gretl_model_get_depvar(pmod);
	    if (dv >= 0 && dv < datainfo->v) {
		default_y = dv;
	    }
	}

	free(xlist);
	xlist = gretl_model_get_x_list(pmod);

	if (pmod->ci == WLS) {
	    wtvar = pmod->nwt;
	} else if (pmod->ci == ARMA) {
	    retrieve_arma_info(pmod);
	} else if (pmod->ci == GARCH) {
	    garch_p = pmod->list[1];
	    garch_q = pmod->list[2];
	} else if (pmod->ci == HECKIT) {
	    retrieve_heckit_info(pmod, &gotinst);
	} else if (pmod->ci == POISSON) {
	    offvar = gretl_model_get_int(pmod, "offset_var");
	} else if (pmod->ci == IVREG) {
	    free(instlist);
	    instlist = gretl_model_get_secondary_list(pmod);
	    if (instlist != NULL) {
		gotinst = 1;
	    }
	    if (pmod->opt & OPT_L) {
		sel_ci = IV_LIML;
	    } else if (pmod->opt & OPT_G) {
		sel_ci = IV_GMM;
	    }
	} else if (pmod->ci == ARCH) {
	    default_order = gretl_model_get_int(pmod, "arch_order");
	} else if (pmod->ci == AR) {
	    retrieve_AR_lags_info(pmod);
	} else if (pmod->ci == AR1) {
	    if (pmod->opt & OPT_P) {
		sel_ci = PWE;
	    } else if (pmod->opt & OPT_H) {
		sel_ci = HILU;
	    } else {
		sel_ci = CORC;
	    }
	} else if (pmod->ci == PANEL) {
	    if (pmod->opt & OPT_F) {
		model_opt |= OPT_F;
	    } else if (pmod->opt & OPT_U) {
		model_opt |= OPT_U;
	    } else if (pmod->opt & OPT_W) {
		sel_ci = PANEL_WLS;
	    } else if (pmod->opt & OPT_B) {
		sel_ci = PANEL_B;
	    }
	    if (pmod->opt & OPT_D) {
		model_opt |= OPT_D;
	    }
	} else if (pmod->ci == LAD) {
	    if (gretl_model_get_int(pmod, "rq")) {
		sel_ci = QUANTREG;
	    }
	    /* FIXME replicate rq_tauvec? */
	} else if (pmod->ci == LOGIT) {
	    if (gretl_model_get_int(pmod, "ordered")) {
		sel_ci = OLOGIT;
	    } else if (gretl_model_get_int(pmod, "multinom")) {
		sel_ci = MLOGIT;
	    }
	} else if (pmod->ci == PROBIT) {
	    if (gretl_model_get_int(pmod, "ordered")) {
		sel_ci = OPROBIT;
	    }
	}   

	if (pmod->opt & OPT_R) {
	    model_opt |= OPT_R;
	}

	y_x_lags_enabled = y_w_lags_enabled = 0;

	if ((dataset_is_time_series(datainfo) || 
	     dataset_is_panel(datainfo)) && 
	    (xlist != NULL || gotinst)) {
	    set_lag_prefs_from_model(dv, xlist, instlist);
	}

	modelspec_dialog(sel_ci);
    } else if (ci == VAR || ci == VECM) {
	GRETL_VAR *var = (GRETL_VAR *) ptr;

	varconst = (var->detflags & DET_CONST);
	vartrend = (var->detflags & DET_TREND);
	want_seasonals = (var->detflags & DET_SEAS);

	if (var->robust) {
	    model_opt |= OPT_R;
	}

	free(veclist);
	veclist = gretl_list_copy(var->ylist);

	free(vecxlist);
	if (var->rlist != NULL) {
	    vecxlist = gretl_lists_join_with_separator(var->xlist,
						       var->rlist);
	} else {
	    vecxlist = gretl_list_copy(var->xlist);
	}

	set_lag_prefs_from_VAR(var->lags, vecxlist);

	default_order = var->order;

	if (ci == VECM) {
	    jrank = gretl_VECM_rank(var);
	    jcase = jcode(var);
	    default_order += 1;
	}

	selection_dialog((ci == VAR)? _("gretl: VAR") : _("gretl: VECM"),
			 do_vector_model, ci);
    } else if (ci == SYSTEM) {
	revise_system_model(ptr);
    }

    model_opt = OPT_NONE;
}

#define UNRESTRICTED N_("U")
#define RESTRICTED   N_("R")

static char varflag[8];

static void set_varflag (const char *s)
{
    *varflag = '\0';
    strncat(varflag, s, 7);
}

int selector_get_depvar_number (const selector *sr)
{
    int ynum = -1;

    if (sr == NULL) {
	sr = open_selector;
    }

    if (sr != NULL && sr->depvar != NULL) {
	const char *s = gtk_entry_get_text(GTK_ENTRY(sr->depvar));

	if (s != NULL && *s != '\0') {
	    if (sr->ci == SAVE_FUNCTIONS) {
		ynum = user_function_index_by_name(s);
	    } else {
		ynum = series_index(datainfo, s);
		if (ynum == datainfo->v) {
		    ynum = -1;
		}
	    }
	}
    }

    return ynum;
}

static gint dblclick_lvars_row (GtkWidget *w, GdkEventButton *event, 
				selector *sr); 

/* when adding a lag to the exogenous vars box for a VECM,
   see if we can find a previously added lag (or non-lag) of
   this variable, and if so use that to set the Restricted/
   Unrestricted flag on the newly added lag
*/

static void set_varflag_from_parent (int v, int lag, GtkTreeModel *mod, 
				     GtkTreeIter *iter)
{
    GtkTreeIter piter;
    gchar *flag;
    gint pv, plag;
    int done = 0;

    if (gtk_tree_model_get_iter_first(mod, &piter)) {
	while (1) {
	    gtk_tree_model_get(mod, &piter, 0, &pv, 1, &plag, -1);
	    if (pv == v && plag != lag) {
		gtk_tree_model_get(mod, &piter, 3, &flag, -1);
		gtk_list_store_set(GTK_LIST_STORE(mod), iter, 3, flag, -1);
		done = 1;
		g_free(flag);
		break;
		
	    }
	    if (!gtk_tree_model_iter_next(mod, &piter)) {  
		break;
	    }
	}
    }

    if (!done) {
	gtk_list_store_set(GTK_LIST_STORE(mod), iter, 3, varflag, -1);
    }
}

static void 
real_varlist_set_var (int v, int lag, GtkTreeModel *mod, GtkTreeIter *iter)
{
    int ncols = gtk_tree_model_get_n_columns(mod);
    GtkListStore *store = GTK_LIST_STORE(mod);

    if (lag == 0) {
	gtk_list_store_set(store, iter, 0, v, 1, 0, 
			   2, datainfo->varname[v], 
			   -1);
    } else {
	char vstr[VNAMELEN + 8];

	sprintf(vstr, "%s(-%d)", datainfo->varname[v], lag);
	gtk_list_store_set(store, iter, 0, v, 1, lag, 
			   2, vstr, -1);
    }

    if (ncols == 4) {
	if (lag == 0) {
	    gtk_list_store_set(store, iter, 3, varflag, -1);
	} else {
	    set_varflag_from_parent(v, lag, mod, iter);
	}
    }
} 

static int 
varlist_remove_var_full (int v, GtkTreeModel *mod, GtkTreeIter *iter)
{
    GtkTreeIter *last = iter;
    int row = -1, ok = 1;
    int tv, i = 0;

#if VLDEBUG
    fprintf(stderr, "\nvarlist_remove_var_full: looking for var %d (%s)\n", v,
	    datainfo->varname[v]);
#endif

    if (gtk_tree_model_get_iter_first(mod, iter)) {
	last = iter;
	while (1) {
	    gtk_tree_model_get(mod, iter, 0, &tv, -1);
#if VLDEBUG
	    fprintf(stderr, "row %d: checking against %d\n", i, tv);
#endif
	    if (tv == v) {
		ok = gtk_list_store_remove(GTK_LIST_STORE(mod), iter);
#if VLDEBUG
		fprintf(stderr, "removed at row %d, now ok = %d\n", i, ok);
#endif
		if (row < 0) {
		    row = i;
		}
		if (ok) {
		    continue;
		} else {
		    break;
		}
	    }
	    if (!gtk_tree_model_iter_next(mod, iter)) {
		/* iter is now invalid */
		iter = last;
		break;
	    }
	    last = iter;
	    i++;
	} 
    } else {
	iter = last;
    }

    return row;
}   

static void
varlist_insert_var_full (int v, GtkTreeModel *mod, GtkTreeIter *iter, 
			 selector *sr, int locus)
{
    int lcontext = 0;

#if VLDEBUG
    fprintf(stderr, "varlist_insert_var_full: starting var %d\n", v);
#endif

    if (v > 0 && dataset_lags_ok(datainfo)) {
	lcontext = sr_get_lag_context(sr, locus);
    }

    if (lcontext) {
	int *laglist = get_lag_pref_as_list(v, lcontext);

	if (laglist != NULL) {
	    int i, row, append = 0;

	    row = varlist_remove_var_full(v, mod, iter);
	    if (row < 0) {
		append = 1;
	    }

#if VLDEBUG
	    fprintf(stderr, "got laglist, done prior removal, row=%d, append=%d\n",
		    row, append);
#endif
	    for (i=1; i<=laglist[0]; i++) {
		if (append) {
		    gtk_list_store_append(GTK_LIST_STORE(mod), iter);
		} else {
		    gtk_list_store_insert(GTK_LIST_STORE(mod), iter, row++);
		}
#if VLDEBUG
		fprintf(stderr, "adding var %d, lag %d\n", v, laglist[i]);
#endif
		real_varlist_set_var(v, laglist[i], mod, iter);
	    }
	    free(laglist);
	} else {
	    lcontext = 0;
	}
    }

    if (lcontext == 0) {
	gtk_list_store_append(GTK_LIST_STORE(mod), iter);
	real_varlist_set_var(v, 0, mod, iter); 
    }
}

static gboolean set_active_var (GtkWidget *widget, GdkEventButton *event,
				selector *sr)
{
    GtkTreeView *view = GTK_TREE_VIEW(widget);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreePath *path;

    if (gtk_tree_view_get_path_at_pos(view, event->x, event->y, &path, 
				      NULL, NULL, NULL)) { 
	GtkTreeIter iter;
	gint varnum, row;

	gtk_tree_model_get_iter(model, &iter, path);
	gtk_tree_model_get(model, &iter, 0, &varnum, -1);
	if (sr != NULL) {
	    sr->active_var = varnum;
	}
	row = tree_path_get_row_number(path);
	g_object_set_data(G_OBJECT(widget), "active_row",
			  GINT_TO_POINTER(row));
	gtk_tree_path_free(path);
    }

    return FALSE;
}

static void list_append_var_simple (GtkListStore *store, GtkTreeIter *iterp, int v)
{
    const char *vname = datainfo->varname[v];

    gtk_list_store_append(store, iterp);    
    gtk_list_store_set(store, iterp, 0, v, 1, 0, 2, vname, -1);
}

static void list_append_var (GtkTreeModel *mod, GtkTreeIter *iter,
			     int v, selector *sr, int locus)
{
    int i, lcontext = 0;

    if (v > 0 && dataset_lags_ok(datainfo)) {
	lcontext = sr_get_lag_context(sr, locus);
    }

    if (lcontext) {
	int *laglist = get_lag_pref_as_list(v, lcontext);

	if (laglist != NULL) {
	    for (i=1; i<=laglist[0]; i++) {
		gtk_list_store_append(GTK_LIST_STORE(mod), iter);
		real_varlist_set_var(v, laglist[i], mod, iter);
	    }
	    free(laglist);
	} else {
	    lcontext = 0;
	}
    }

    if (lcontext == 0) {
	gtk_list_store_append(GTK_LIST_STORE(mod), iter);
	real_varlist_set_var(v, 0, mod, iter);
    }
}

/* build a new liststore and associated tree view, and pack into the
   given box */

static GtkWidget *var_list_box_new (GtkBox *box, selector *sr, int locus) 
{
    GtkListStore *store; 
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer; 
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    int flagcol = 0;
    int width = 120;
    int height = -1;
    
    if (USE_RXLIST(sr->ci) && locus == SR_RVARS2) {
	store = gtk_list_store_new(4, G_TYPE_INT, G_TYPE_INT, 
				   G_TYPE_STRING, G_TYPE_STRING);
	flagcol = 1;
    } else {
	store = gtk_list_store_new(3, G_TYPE_INT, G_TYPE_INT, G_TYPE_STRING);
    }

    view = gtk_tree_view_new_with_model(GTK_TREE_MODEL(store));
    g_object_unref(G_OBJECT(store));

    g_object_set_data(G_OBJECT(view), "sr", sr);
    g_object_set_data(G_OBJECT(view), "locus", GINT_TO_POINTER(locus));

    renderer = gtk_cell_renderer_text_new();
    g_object_set(renderer, "ypad", 0, NULL);
    column = gtk_tree_view_column_new_with_attributes(NULL,
						      renderer,
						      "text", 2,
						      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);	

    if (flagcol) {
	column = gtk_tree_view_column_new_with_attributes(NULL,
							  renderer,
							  "text", 3,
							  NULL);
	gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
    }	

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(view), FALSE);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    gtk_tree_selection_set_mode(select, GTK_SELECTION_MULTIPLE);
    g_signal_connect(G_OBJECT(view), "motion-notify-event",
		     G_CALLBACK(listbox_drag), NULL);

    if (locus == SR_LVARS) { 
	/* left-hand box with the selectable vars */
	g_signal_connect(G_OBJECT(view), "button-press-event",
			 G_CALLBACK(lvars_right_click),
			 sr);
	g_signal_connect(G_OBJECT(view), "button-press-event",
			 G_CALLBACK(set_active_var),
			 sr);
	g_signal_connect(G_OBJECT(view), "button-press-event",
			 G_CALLBACK(dblclick_lvars_row),
			 sr);
    } else if (locus == SR_RVARS1 || locus == SR_RVARS2) { 
	/* lists of selected items */
	g_signal_connect(G_OBJECT(view), "button-press-event",
			 G_CALLBACK(set_active_var),
			 NULL);
	if (flagcol) {
	    g_signal_connect(G_OBJECT(view), "button-press-event",
			     G_CALLBACK(listvar_flagcol_click),
			     view);
	} else {
	    g_signal_connect(G_OBJECT(view), "button-press-event",
			     G_CALLBACK(listvar_special_click),
			     view);
	}
    } 

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW (scroller),
				   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW (scroller),
					GTK_SHADOW_IN);    
    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(box, scroller, TRUE, TRUE, 0);

    width *= gui_scale;
    gtk_widget_set_size_request(view, width, height);
    gtk_widget_show(view);
    gtk_widget_show(scroller);

    return view;
}

/* add to "extra" var slot the current selection from sr->lvars */

static void real_set_extra_var (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, selector *sr)
{
    gint vnum;
    gchar *vname;

    gtk_tree_model_get(model, iter, 0, &vnum, 2, &vname, -1);

    if (sr->ci == HECKIT) {
	if (!gretl_isdummy(datainfo->t1, datainfo->t2, Z[vnum])) {
	    errbox(_("The variable '%s' is not a 0/1 variable."), vname);
	    return;
	}
    }
    
    gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->extra[0]), "data",
		      GINT_TO_POINTER(vnum));
}

static void set_extra_var_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					real_set_extra_var,
					sr);
}

static void real_set_third_var (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, selector *sr)
{
    gint vnum;
    gchar *vname;
    
    gtk_tree_model_get(model, iter, 0, &vnum, 2, &vname, -1);
    gtk_entry_set_text(GTK_ENTRY(sr->rvars1), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->rvars1), "data",
		      GINT_TO_POINTER(vnum));
}

static void set_third_var_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					real_set_third_var,
					sr);
}

/* when adding a variable to the Endogenous or Exogenous listing
   for VAR, VECM, etc., check that the variable in question is not
   present in the other listing: if it is, then remove it.
*/

static void vecx_cleanup (selector *sr, int new_locus)
{
    GtkTreeModel *targ, *src;
    GtkTreeIter iter;
    gboolean ok;
    int *srclist = NULL;
    gint v;

    if (new_locus == SR_RVARS1) {
	src = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
	targ = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
    } else {
	src = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
	targ = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    }

    if (gtk_tree_model_get_iter_first(src, &iter)) {
	do {
	    ok = TRUE;
	    gtk_tree_model_get(src, &iter, 0, &v, -1);
	    gretl_list_append_term(&srclist, v);
        } while (ok && gtk_tree_model_iter_next(src, &iter));
    }

    if (srclist != NULL && gtk_tree_model_get_iter_first(targ, &iter)) {
	do {
	    ok = TRUE;
	    gtk_tree_model_get(targ, &iter, 0, &v, -1);
	    if (in_gretl_list(srclist, v)) {
		ok = gtk_list_store_remove(GTK_LIST_STORE(targ), &iter);
	    }
        } while (ok && gtk_tree_model_iter_next(targ, &iter));
    }

    free(srclist);
}

static void 
maybe_insert_or_revise_depvar_lags (selector *sr, int v, int lcontext,
				    int revise)
{
    int *laglist = NULL;
    GtkWidget *w;
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int locus, append = 1;
    int modv, row = 0;
    int jmin = 0, jmax = 1;
    int i, j;

    if (lcontext == LAG_Y_X) {
	jmin = 0;
	jmax = 1;
    } else if (lcontext == LAG_Y_W) {
	/* instrument, not indep var */
	jmin = 1;
	jmax = 2;
    }

    for (j=jmin; j<jmax; j++) {

	w = (j > 0)? sr->rvars2: sr->rvars1;
	if (w == NULL) {
	    return;
	}

	locus = (j > 0)? SR_RVARS2: SR_RVARS1;
	
	mod = gtk_tree_view_get_model(GTK_TREE_VIEW(w));

	if (lcontext == 0) {
	    lcontext = (j > 0)? LAG_Y_W : LAG_Y_X;
	}

	laglist = get_lag_pref_as_list(v, lcontext);
	if (laglist == NULL) {
	    if (revise) {
		varlist_remove_var_full(v, mod, &iter);
	    }
	    return;
	}

	varlist_remove_var_full(v, mod, &iter);

	if (gtk_tree_model_get_iter_first(mod, &iter)) {
	    do {
		gtk_tree_model_get(mod, &iter, 0, &modv, -1);
		if (modv > 0) {
		    append = 0;
		    break;
		}
		row++;
	    } while (gtk_tree_model_iter_next(mod, &iter));
	} 

	for (i=1; i<=laglist[0]; i++) {
	    if (append) {
		gtk_list_store_append(GTK_LIST_STORE(mod), &iter);
	    } else {
		gtk_list_store_insert(GTK_LIST_STORE(mod), &iter, row++);
	    }
#if VLDEBUG
	    fprintf(stderr, "depvar_lags: adding var %d, lag %d\n", v, laglist[i]);
#endif
	    real_varlist_set_var(v, laglist[i], mod, &iter); 
	}    

	free(laglist);
    }
}

static void maybe_insert_depvar_lags (selector *sr, int v, int lcontext)
{
    maybe_insert_or_revise_depvar_lags(sr, v, lcontext, 0);
}

static void maybe_revise_depvar_lags (selector *sr, int v, int lcontext)
{
    maybe_insert_or_revise_depvar_lags(sr, v, lcontext, 1);
}

static void remove_as_indep_var (selector *sr, gint v)
{
    GtkTreeView *view = GTK_TREE_VIEW(sr->rvars1);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeIter iter;
    gboolean ok;
    gint xnum;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
	do {
	    ok = TRUE;
	    gtk_tree_model_get(model, &iter, 0, &xnum, -1);
	    if (xnum == v) {
		ok = gtk_list_store_remove(GTK_LIST_STORE(model), &iter);
	    }
        } while (ok && gtk_tree_model_iter_next(model, &iter));
    }
}

static void dependent_var_cleanup (selector *sr, int newy)
{
    int oldy = selector_get_depvar_number(sr);

    if (GTK_IS_TREE_VIEW(sr->rvars1) && oldy != newy) {
	if (oldy > 0) {
	    remove_as_indep_var(sr, oldy); /* lags business */
	}
	y_x_lags_enabled = 0;
	y_w_lags_enabled = 0;
	remove_as_indep_var(sr, newy);
    }
}

static void set_public_from_active (selector *sr)
{
    gint v = sr->active_var;

    if (sr->depvar == NULL) {
	return;
    }

    dependent_var_cleanup(sr, v);
    gtk_entry_set_text(GTK_ENTRY(sr->depvar), 
		       user_function_name_by_index(v));
}

static void set_dependent_var_from_active (selector *sr)
{
    gint v = sr->active_var;

    if (sr->depvar == NULL) {
	return;
    }

    /* models: if we select foo as regressand, remove it from the list
       of regressors if need be; also remove lags associated with the
       previous dependent var, if any.
    */
    if (MODEL_CODE(sr->ci)) {
	dependent_var_cleanup(sr, v);
    }

    gtk_entry_set_text(GTK_ENTRY(sr->depvar), datainfo->varname[v]);

    if (select_lags_depvar(sr->ci)) {
	maybe_insert_depvar_lags(sr, v, 0);
    }
}

static void real_set_dependent_var (GtkTreeModel *model, GtkTreePath *path,
				    GtkTreeIter *iter, selector *sr)
{
    gchar *vname;
    gint v;

    gtk_tree_model_get(model, iter, 0, &v, 2, &vname, -1);

    if (MODEL_CODE(sr->ci)) {
	dependent_var_cleanup(sr, v);
    }

    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vname);

    if (select_lags_depvar(sr->ci)) {
	maybe_insert_depvar_lags(sr, v, 0);
    }
}

static void set_dependent_var_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (sr->depvar == NULL) return;

    selection = gtk_tree_view_get_selection (GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach (selection, 
					 (GtkTreeSelectionForeachFunc) 
					 real_set_dependent_var,
					 sr);
}

static void set_right_var_from_main (GtkTreeModel *model, GtkTreePath *path,
				     GtkTreeIter *iter, selector *sr)
{
    GtkTreeModel *rmod;
    GtkTreeIter r_iter;
    gchar *vnum = NULL;
    int v;

    gtk_tree_model_get(model, iter, 0, &vnum, -1);

    rmod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (rmod == NULL) {
	g_free(vnum);
	return;
    }

    v = atoi(vnum);

    if (gtk_tree_model_get_iter_first(rmod, &r_iter)) {
	while (gtk_tree_model_iter_next(rmod, &r_iter)) {
	    ;
	}
    }

    gtk_list_store_append(GTK_LIST_STORE(rmod), &r_iter);
    real_varlist_set_var(v, 0, rmod, &r_iter);   

    g_free(vnum);
}

static void set_vars_from_main (selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					set_right_var_from_main,
					sr);
}

/* Append a specified variable in the SR_RVARS1 locus: used when
   saving data and there's only one variable to save.
*/

static void select_singleton (selector *sr)
{
    GtkTreeModel *lmod, *rmod;
    GtkTreeIter iter;
    int v;

    lmod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars));
    if (lmod == NULL) {
	return;
    }    

    rmod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (rmod == NULL) {
	return;
    }

    gtk_tree_model_get_iter_first(lmod, &iter);
    gtk_tree_model_get(lmod, &iter, 0, &v, -1);

    gtk_tree_model_get_iter_first(rmod, &iter);
    gtk_list_store_append(GTK_LIST_STORE(rmod), &iter);
    gtk_list_store_set(GTK_LIST_STORE(rmod), &iter, 
		       0, v, 1, 0, 2, datainfo->varname[v], -1);
}

static int varflag_dialog (int v)
{
    const char *opts[] = {
	N_("Unrestricted"),
	N_("Restricted")
    };
    gchar *title, *label;
    int ret;

    title = g_strdup_printf("gretl: %s", _("add exogenous variable"));
    label = g_strdup_printf(_("Status of '%s' in VECM:"), datainfo->varname[v]);

    ret = radio_dialog(title, label, opts, 2, 0, 0);
    if (ret == 0) {
	set_varflag(UNRESTRICTED);
    } else if (ret == 1) {
	set_varflag(RESTRICTED);
    }

    g_free(title);
    g_free(label);

    return ret;
}

static void real_add_generic (GtkTreeModel *srcmodel, GtkTreeIter *srciter, 
			      selector *sr, int locus)
{
    GtkWidget *list;
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *vname = NULL;
    gint v, xnum;
    gint already_there = 0;
    gint at_max = 0;
    gint keep_names = 0;
    int nvars = 0;

    list = (locus == SR_RVARS2)? sr->rvars2 : sr->rvars1;

    if (!GTK_IS_TREE_VIEW(list)) {
	return;
    }

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(list));
    if (model == NULL) {
	return;
    }

    keep_names = 
	GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->lvars), 
					  "keep_names"));

    if (keep_names) {
	gtk_tree_model_get(srcmodel, srciter, 0, &v, 2, &vname, -1);
    } else {
	gtk_tree_model_get(srcmodel, srciter, 0, &v, -1);
    }

    nvars = 0;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
	do {
	    nvars++;
	    if (!at_max && selection_at_max(sr, nvars)) {
		at_max = 1;
	    }
	    if (!already_there) {
		gtk_tree_model_get(model, &iter, 0, &xnum, -1);
		if (xnum == v) {
		    already_there = 1; 
		}
	    }
	} while (gtk_tree_model_iter_next(model, &iter));
    }

    if (!already_there && !at_max) {
	if (vname != NULL) {
	    int ncols = gtk_tree_model_get_n_columns(model);

	    gtk_list_store_append(GTK_LIST_STORE(model), &iter);
	    gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
			       0, v, 1, 0, 2, vname, -1);
	    if (ncols == 4) {
		gtk_list_store_set(GTK_LIST_STORE(model), &iter, 
				   3, varflag, -1);
	    } 
	    g_free(vname);
	} else {
	    if (locus == SR_RVARS2 && USE_RXLIST(sr->ci)) {
		if (varflag_dialog(v) < 0) {
		    return;
		}
	    }
#if VLDEBUG
	    fprintf(stderr, "real_add_generic: calling varlist_insert_var_full\n");
#endif
	    varlist_insert_var_full(v, model, &iter, sr, locus);
	}
	nvars++;
    }

    if (sr->add_button != NULL && at_max) {
	gtk_widget_set_sensitive(sr->add_button, FALSE);
    }

    if (nvars > 0) {
	if (lags_button_relevant(sr, locus)) {
	    gtk_widget_set_sensitive(sr->lags_button, TRUE);
	} else if (VECLAGS_CODE(sr->ci) && locus == SR_RVARS1 &&
		   sr->extra[EXTRA_LAGS] != NULL) {
	    gtk_widget_set_sensitive(sr->extra[EXTRA_LAGS], TRUE);
	}
    }

    if (USE_VECXLIST(sr->ci)) {
	vecx_cleanup(sr, locus);
    }
}

static void add_to_rvars1 (GtkTreeModel *model, GtkTreePath *path,
			   GtkTreeIter *iter, selector *sr)
{
    /* models: don't add the regressand to the list of regressors;
       functions: keep public and private interfaces distinct */
    if (MODEL_CODE(sr->ci) || sr->ci == SAVE_FUNCTIONS) {
	gint xnum, ynum;
    
	gtk_tree_model_get(model, iter, 0, &xnum, -1);
	ynum = selector_get_depvar_number(sr);
	if (xnum == ynum) {
	    return;
	}
    }

    real_add_generic(model, iter, sr, SR_RVARS1);
}

static void add_to_rvars2 (GtkTreeModel *model, GtkTreePath *path,
			   GtkTreeIter *iter, selector *sr)
{
    real_add_generic(model, iter, sr, SR_RVARS2);
}

static void add_to_rvars2_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					add_to_rvars2,
					sr);
}

static void add_all_to_rvars1_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (!GTK_IS_TREE_VIEW(sr->lvars)) {
	return;
    }

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_select_all(selection);
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					add_to_rvars1,
					sr);
}

static void add_to_rvars1_callback (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    if (!GTK_IS_TREE_VIEW(sr->lvars)) {
	return;
    }

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection, 
					(GtkTreeSelectionForeachFunc) 
					add_to_rvars1,
					sr);
}

static void remove_from_right_callback (GtkWidget *w, gpointer data)
{
    GtkTreeView *view = GTK_TREE_VIEW(data);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(view);
    GtkTreePath *path;
    GtkTreeIter iter, last;
    selector *sr;
    int context = 0;
    int nsel = 0;

    if (model == NULL || selection == NULL) {
	return;
    }

    /* get to the last row */
    if (gtk_tree_model_get_iter_first(model, &iter)) {
	last = iter;
	nsel = 1;
	while (gtk_tree_model_iter_next(model, &iter)) {
	    last = iter;
	    nsel++;
	}
    } else {
	return;
    }
    
    context = lag_context_from_widget(GTK_WIDGET(view));
    
    /* work back up, deleting selected rows */
    path = gtk_tree_model_get_path(model, &last);
    while (1) {
	if (gtk_tree_model_get_iter(model, &last, path) &&
	    gtk_tree_selection_iter_is_selected(selection, &last)) {
	    if (context) {
		gint v, lag;

		gtk_tree_model_get(model, &last, 0, &v, 1, &lag, -1);
		remove_specific_lag(v, lag, context);
	    }
	    gtk_list_store_remove(GTK_LIST_STORE(model), &last);
	    nsel--;
	}
	if (!gtk_tree_path_prev(path)) {
	    break;
	}
    } 

    sr = g_object_get_data(G_OBJECT(data), "selector");

    if (sr != NULL && sr->add_button != NULL &&
	!GTK_WIDGET_SENSITIVE(sr->add_button) &&
	!selection_at_max(sr, nsel)) {
	gtk_widget_set_sensitive(sr->add_button, TRUE);
    }

    if (context && sr != NULL && sr->lags_button != NULL) {
	if (nsel == 0) {
	    gtk_widget_set_sensitive(sr->lags_button, FALSE);
	}
    }
}

/* callbacks from button presses in list boxes: double and right
   clicks do special stuff */

static gint 
dblclick_lvars_row (GtkWidget *w, GdkEventButton *event, selector *sr) 
{
    if (event != NULL && event->type == GDK_2BUTTON_PRESS) { 
	if (sr->ci == SAVE_FUNCTIONS) {
	    set_public_from_active(sr);
	} else {
	    set_dependent_var_from_active(sr);
	    if (sr->default_check != NULL) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
					     TRUE);
	    }
	}
    }

    return FALSE;
}

/* flip the flag that designates an exogenous veriable in a VECM
   as either Unrestricted or Restricted (to the cointegrating
   space) 
*/

static void maybe_flip_vecm_flag (GtkTreeModel *model, GtkTreePath *path,
				  GtkTreeIter *iter, GtkWidget *w)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gchar *orig = NULL, *repl = NULL;

    gtk_tree_model_get(model, iter, 3, &orig, -1);

    if (orig != NULL) {
	repl = (i == 0)? "U" : "R";
	if (strcmp(orig, repl)) {
	    gtk_list_store_set(GTK_LIST_STORE(model), iter, 3, repl, -1);
	}
	g_free(orig);
    }
}

static void flag_popup_activated (GtkWidget *w, gpointer data)
{
    GtkTreeView *view = GTK_TREE_VIEW(data);
    GtkTreeSelection *sel;

    sel = gtk_tree_view_get_selection(view);
    gtk_tree_selection_selected_foreach(sel, 
					(GtkTreeSelectionForeachFunc) 
					maybe_flip_vecm_flag,
					w);
    gtk_widget_destroy(w);
}

static void flag_popup_delete (GtkWidget *w, gpointer data)
{
    gtk_widget_destroy(GTK_WIDGET(data));
}

static void create_flag_item (GtkWidget *popup, int i, GtkWidget *view)
{
    static char *flag_strs[] = {
	N_("Unrestricted"),
	N_("Restricted")
    };
    GtkWidget *item;

    item = gtk_menu_item_new_with_label(_(flag_strs[i]));
    g_object_set_data(G_OBJECT(item), "opt", GINT_TO_POINTER(i));
    g_signal_connect(G_OBJECT(item), "activate",
		     G_CALLBACK(flag_popup_activated),
		     view);
    g_signal_connect(G_OBJECT(item), "destroy",
		     G_CALLBACK(flag_popup_delete),
		     popup);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);
}

static gint listvar_flagcol_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data)
{
    GtkWidget *view = GTK_WIDGET(data);
    GdkWindow *top = gtk_widget_get_parent_window(view);
    GdkModifierType mods;
    GtkWidget *flag_popup;
    int i;

    gdk_window_get_pointer(top, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON3_MASK) {
	flag_popup = gtk_menu_new();
	for (i=0; i<2; i++) {
	    create_flag_item(flag_popup, i, view);
	}
	gtk_menu_popup(GTK_MENU(flag_popup), NULL, NULL, NULL, NULL,
		       event->button, event->time);
	return TRUE;
    }

    return FALSE;
}

static gint listvar_special_click (GtkWidget *widget, GdkEventButton *event, 
				   gpointer data)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(GTK_WIDGET(data));
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    /* FIXME below: does this do anything useful?  I think
       it's dead code (AC, 2007-10-31). */

    if (mods & GDK_BUTTON2_MASK) {
	gtk_tree_view_set_reorderable(GTK_TREE_VIEW(data), TRUE);
    } else {
	gtk_tree_view_set_reorderable(GTK_TREE_VIEW(data), FALSE);
    }

    if (mods & GDK_BUTTON3_MASK) {
	remove_from_right_callback(NULL, data);
	return TRUE;
    } 

    return FALSE;
}

static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event, 
			       selector *sr)
{
    GdkWindow *topwin;
    GdkModifierType mods;

    topwin = gtk_widget_get_parent_window(sr->lvars);
    gdk_window_get_pointer(topwin, NULL, NULL, &mods); 

    if (mods & GDK_BUTTON3_MASK) {
	add_to_rvars1_callback(NULL, sr);
	return TRUE;
    }

    return FALSE;
}

/* end special click callbacks */

static void varlist_insert_const (GtkWidget *w)
{
    GtkTreeModel *mod = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
    GtkTreeIter iter;

    gtk_tree_model_get_iter_first(mod, &iter);
    gtk_list_store_append(GTK_LIST_STORE(mod), &iter);
    gtk_list_store_set(GTK_LIST_STORE(mod), &iter, 
		       0, 0, 1, 0, 2, "const", -1);
}

static void clear_vars (GtkWidget *w, selector *sr)
{
    GtkTreeSelection *selection;

    /* deselect all vars on left */
    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_unselect_all(selection);

    /* clear dependent var slot */
    if (sr->depvar != NULL) {
	gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
	if (sr->default_check != NULL) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
					 FALSE);
	}
	default_y = -1;
    }

    /* extra variable entry? */
    if (GTK_IS_ENTRY(sr->extra[0])) {
	gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), "");
    }

    if (THREE_VARS_CODE(sr->ci)) {
	/* clear special slot */
	gtk_entry_set_text(GTK_ENTRY(sr->rvars1), "");
    } else {
	/* empty lower right variable list */
	clear_varlist(sr->rvars1);
	if (sr->add_button != NULL) {
	    gtk_widget_set_sensitive(sr->add_button, TRUE);
	}
    }

    if (MODEL_CODE(sr->ci) && sr->ci != ARMA) {
	/* insert default const in regressors box */
	varlist_insert_const(sr->rvars1);
    }

    if (sr->rvars2 != NULL) {
	/* empty lower right variable list */
	clear_varlist(sr->rvars2);
	if (USE_ZLIST(sr->ci)) {
	    varlist_insert_const(sr->rvars2);
	}
    }

    if (sr->lags_button != NULL) {
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }

    if (VECLAGS_CODE(sr->ci) && sr->extra[EXTRA_LAGS] != NULL) {
	gtk_widget_set_sensitive(sr->extra[EXTRA_LAGS], FALSE);
	gtk_widget_set_sensitive(sr->extra[0], TRUE);
    }

    clear_selector();
}

static gint varlist_row_count (selector *sr, int locus, int *realrows)
{
    int lcontext = 0;
    GtkWidget *w;
    GtkTreeModel *mod;
    GtkTreeIter iter;
    gint v, lag, n = 0;

    w = (locus == SR_RVARS1)? sr->rvars1 : sr->rvars2;

    if (w == NULL || !GTK_WIDGET_IS_SENSITIVE(w)) {
	return 0;
    }

    if (realrows != NULL) {
	lcontext = sr_get_lag_context(sr, locus);
	*realrows = 0;
    }

    if (w != NULL) {
	mod = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
	if (GTK_IS_TREE_MODEL(mod) && 
	    gtk_tree_model_get_iter_first(mod, &iter)) {
	    do {
		n++;
		if (lcontext) {
		    gtk_tree_model_get(mod, &iter, 0, &v, 1, &lag, -1);
		    if (!is_lag_dummy(v, lag, lcontext)) {
			*realrows += 1;
		    }
		}
	    } while (gtk_tree_model_iter_next(mod, &iter));
	}	    
    }

    if (realrows != NULL && lcontext == 0) {
	*realrows = n;
    }

    return n;
}

static void topslot_empty (int ci)
{
    switch (ci) {
    case GR_XY:
    case GR_3D:
    case GR_IMP:
	warnbox(_("You must select an X-axis variable"));
	break;
    case SCATTERS:
	warnbox(_("You must select a Y-axis variable"));
	break;
    case SAVE_FUNCTIONS:
	warnbox(_("You must specify a public interface"));
	break;
    case INTREG:
	warnbox(_("You must select a lower bound variable"));
	break;
    default:
	warnbox(_("You must select a dependent variable"));
    }
}

static void reverse_list (char *list)
{
    char *tmp, *p;
    char istr[8];

    p = strchr(list, ';');
    if (p == NULL) return;

    tmp = malloc(strlen(list) + 4);
    if (tmp == NULL) return;

    sscanf(list, "%7s", istr);

    strcpy(tmp, p + 2);
    strcat(tmp, " ; ");
    strcat(tmp, istr);

    strcpy(list, tmp);

    free(tmp);
}

enum cmdlist_codes {
    ADD_NOW,
    ADD_AT_END
};

static int add_to_cmdlist (selector *sr, const char *add)
{
    int n = strlen(sr->cmdlist);
    char *cmdlist = NULL;
    int err = 0;

    if (n % MAXLEN > MAXLEN - 32) {
	int blocks = 2 + n / MAXLEN;

	cmdlist = realloc(sr->cmdlist, blocks * MAXLEN);
	if (cmdlist == NULL) {
	    err = 1;
	} else {
	    sr->cmdlist = cmdlist;
	}
    }

    if (!err) {
	strcat(sr->cmdlist, add);
    }

    return err;
}

static char *arma_lag_string (char *targ, const char *s)
{
    while (*s == ' ') s++;

    if (*s == '\0') {
	strcpy(targ, "0 ");
    } else if (isalpha(*s) || *s == '{') {
	sprintf(targ, "%s ", s);
    } else {
	sprintf(targ, "{%s} ", s);
    } 

    charsub(targ, ',', ' ');

    return targ;
}

static void arma_spec_to_cmdlist (selector *sr)
{
    const char *txt;
    char s[32];

    free(arlags);
    arlags = NULL;

    free(malags);
    malags = NULL;

    if (GTK_WIDGET_SENSITIVE(sr->extra[0])) {
	arma_p = spinner_get_int(sr->extra[0]);
	sprintf(s, "%d ", arma_p);
	add_to_cmdlist(sr, s);
    } else {
	txt = gtk_entry_get_text(GTK_ENTRY(sr->extra[1]));
	add_to_cmdlist(sr, arma_lag_string(s, txt)); 
	arlags = gretl_strdup(txt);
    }

    arima_d = spinner_get_int(sr->extra[2]);
    sprintf(s, "%d ", arima_d);
    add_to_cmdlist(sr, s);

    if (GTK_WIDGET_SENSITIVE(sr->extra[3])) {
	arma_q = spinner_get_int(sr->extra[3]);
	sprintf(s, "%d ; ", arma_q);
	add_to_cmdlist(sr, s);
    } else {
	txt = gtk_entry_get_text(GTK_ENTRY(sr->extra[4]));
	add_to_cmdlist(sr, arma_lag_string(s, txt));
	add_to_cmdlist(sr, "; ");
	malags = gretl_strdup(txt);
    }

    if (sr->extra[5] != NULL) {
	arma_P = spinner_get_int(sr->extra[5]);
	arima_D = spinner_get_int(sr->extra[6]);
	arma_Q = spinner_get_int(sr->extra[7]);

	if (arma_P > 0 || arima_D > 0 || arma_Q > 0) {
	    sprintf(s, "%d %d %d ; ", arma_P, arima_D, arma_Q);
	    add_to_cmdlist(sr, s);
	}
    }
}

static void add_pdq_vals_to_cmdlist (selector *sr)
{
    char s[32] = {0};

    if (sr->ci == GARCH) {
	garch_p = spinner_get_int(sr->extra[0]);
	garch_q = spinner_get_int(sr->extra[1]);
	sprintf(s, "%d %d ; ", garch_p, garch_q);
    } else if (sr->ci == ARBOND) {
	int p = spinner_get_int(sr->extra[0]);

	sprintf(s, "%d ; ", p);
    } else if (sr->ci == ARCH) {
	int p = spinner_get_int(sr->extra[0]);

	sprintf(s, "%d ", p);
    } 

    add_to_cmdlist(sr, s);
}

static void read_ellipse_alpha (selector *sr)
{
    if (sr->extra[0] != NULL) {
	char s[8];
	double cval;

	cval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[0]));
	sprintf(s, "%g ", 1 - cval);
	add_to_cmdlist(sr, s);
    }
}

/* Take the stored preferred laglist for a variable (if any) and
   convert to a string specification for adding to the regression
   command line.
*/

static char *
discrete_lags_string (const char *vname, const int *laglist,
		      char context)
{
    int len = strlen(vname) + 4;
    gchar *tmp;
    char *s;
    int i, li;
    int err = 0;

    s = malloc(len + laglist[0] * 6);

    if (s != NULL) {
	sprintf(s, " %s(", vname);
	for (i=1; i<=laglist[0]; i++) {
	    li = laglist[i];
	    if (li > 999) {
		err = 1;
		break;
	    }
	    tmp = g_strdup_printf("%s%d", (li > 0)? "-" : "", li);
	    strcat(s, tmp);
	    if (i < laglist[0]) {
		strcat(s, ", ");
	    } else {
		strcat(s, ")");
	    }
	    g_free(tmp);
	}

	if (err) {
	    free(s);
	    s = NULL;
	}
    }

    return s;
}  

int selector_get_VAR_order (const selector *sr)
{
    return gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(sr->extra[0]));
}

/* for use in constructing command list, possibly with
   embedded lags */

static char *get_lagpref_string (int v, char context,
				 selector *sr)
{
    const char *vname = datainfo->varname[v];
    const int *laglist;
    int lmin, lmax;
    char *s = NULL;

    if (v == 0) {
	if (context == LAG_Y_X || context == LAG_Y_W) {
	    /* const as dependent: empty lags string */
	    return g_strdup("");
	} else {
	    /* const as indep var: just return itself */
	    return g_strdup(" 0");
	}
    }

    get_lag_preference(v, &lmin, &lmax, &laglist, context, sr);

    if (laglist != NULL) {
	s = discrete_lags_string(vname, laglist, context);
    } else if (lmin != lmax) {
	s = g_strdup_printf(" %s(%s%d to -%d)", vname, (lmin > 0)? "-" : "",
			    lmin, lmax);
    } else if (lmin != 0) {
	s = g_strdup_printf(" %s(%s%d)", vname, (lmin > 0)? "-" : "",
			    lmin);
    } else if (context != LAG_Y_X && context != LAG_Y_W) {
	s = g_strdup_printf(" %d", v);
    }

#if LDEBUG
    if (s != NULL) {
	fprintf(stderr, "get_lagpref_string (v=%d, context=%d):\n"
		" constructed s = '%s'\n", v, (int) context, s);
    }
#endif

    return s;
} 

static int maybe_resize_recorder_lists (selector *sr, int n)
{
    int err = 0;

    if (MODEL_CODE(sr->ci) || VEC_CODE(sr->ci)) {
	int *newlist;

	if (MODEL_CODE(sr->ci)) {
	    newlist = gretl_list_resize(&xlist, n);
	} else {
	    newlist = gretl_list_resize(&veclist, n);
	}
	if (newlist == NULL) {
	    err = E_ALLOC;
	}
    }

    return err;
}

static int maybe_resize_exog_recorder_lists (selector *sr, int n)
{
    int err = 0;

    if (USE_ZLIST(sr->ci) || USE_VECXLIST(sr->ci)) {
	int *newlist;

	if (USE_ZLIST(sr->ci)) {
	    newlist = gretl_list_resize(&instlist, n);
	} else {
	    newlist = gretl_list_resize(&vecxlist, n);
	}
	if (newlist == NULL) {
	    err = E_ALLOC;
	}
    } 

    return err;
}

static void get_rvars1_data (selector *sr, int rows, int context)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint rvar, lag;
    gchar *rvstr;
    int added = 0;
    int gotconst = 0;
    int i, j = 1;

    if (SAVE_DATA_ACTION(sr->ci) && sr->ci != COPY_CSV && rows == sr->n_left) {
	/* saving/exporting all available series: leave the list blank
	   in case it overflows */
	return;
    }   

    if (sr->rvars1 == NULL) {
	return;
    }

    sr->n_left = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<rows; i++) {

	gtk_tree_model_get(model, &iter, 0, &rvar, 1, &lag, -1);
	if (rvar == 0) {
	    gotconst = 1;
	}

	if (is_lag_dummy(rvar, lag, context)) {
	    gtk_tree_model_iter_next(model, &iter);
	    continue;
	}

	if (context) {
	    rvstr = get_lagpref_string(rvar, context, sr);
	} else {
	    rvstr = g_strdup_printf(" %d", rvar);
	}

	if (rvstr == NULL) {
	    sr->error = E_ALLOC;
	    break;
	} else {
	    add_to_cmdlist(sr, rvstr);
	    g_free(rvstr);
	    added++;
	}

	/* save for future reference */
	if (MODEL_CODE(sr->ci) && xlist != NULL) {
	    xlist[j++] = rvar;
	} else if (VEC_CODE(sr->ci) && veclist != NULL) {
	    veclist[j++] = rvar;
	}

	gtk_tree_model_iter_next(model, &iter);
    }

    if (sr->ci == GARCH && !gotconst) {
	/* if the user has deliberately removed the constant,
	   don't put it back in */
	sr->opts |= OPT_N;
    } 

    if (sr->ci == ARMA && added && !(sr->opts & OPT_N)) {
	/* add const explicitly unless forbidden */
	add_to_cmdlist(sr, " 0");
    }
}

/* VECM: parse out the exogenous vars as either restricted
   or unrestricted.  Right now we don't attempt to combine
   this with auto lag selection ("Lags" button), but maybe
   we should */

static int get_vecm_exog_list (selector *sr, int rows, 
			       GtkTreeModel *mod)
{
    GtkTreeIter iter;
    int *xlist = NULL, *zlist = NULL;
    gchar *flag = NULL;
    gchar *tmp = NULL;
    int lag, i, j = 1;
    int v, err = 0;

    gtk_tree_model_get_iter_first(mod, &iter);

    for (i=0; i<rows && !err; i++) {
	flag = NULL;
	lag = 0;

	gtk_tree_model_get(mod, &iter, 0, &v, 1, &lag, 3, &flag, -1);

	if (flag == NULL) {
	    err = E_DATA;
	} else if (lag != 0) {
	    v = laggenr(v, lag, &Z, datainfo);
	    if (v < 0) {
		err = E_DATA;
	    }
	}

	if (!err) {
	    if (*flag == 'U') {
		gretl_list_append_term(&xlist, v);
		if (xlist == NULL) {
		    err = E_ALLOC;
		}
	    } else if (*flag == 'R') {
		gretl_list_append_term(&zlist, v);
		if (zlist == NULL) {
		    err = E_ALLOC;
		}		
	    } else {
		err = E_DATA;
	    }
	    g_free(flag);
	} 

	gtk_tree_model_iter_next(mod, &iter);
    }

    if (!err) {
	if (xlist != NULL) {
	    tmp = gretl_list_to_string(xlist);
	    add_to_cmdlist(sr, tmp);
	    g_free(tmp);
	    for (i=1; i<=xlist[0]; i++) {
		vecxlist[j++] = xlist[i];
	    }
	    free(xlist);
	} 
	if (zlist != NULL) {
	    int n = vecxlist[0];

	    gretl_list_resize(&vecxlist, n + 1);
	    if (vecxlist == NULL) {
		err = E_ALLOC;
	    } else {
		add_to_cmdlist(sr, " ;");
		vecxlist[j++] = LISTSEP;
		tmp = gretl_list_to_string(zlist);
		add_to_cmdlist(sr, tmp);
		g_free(tmp);
		for (i=1; i<=zlist[0]; i++) {
		    vecxlist[j++] = zlist[i];
		}
	    }
	    free(zlist);	    
	}
    }

    if (err) {
	gui_errmsg(err);
    }    

    return err;
}

/* get the component of the model specification from the secondary
   list box on the right, if applicable */

static int get_rvars2_data (selector *sr, int rows, int context)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint exog, lag, ynum;
    int *reclist = NULL;
    gchar *tmp;
    int ncols, i, j = 1;
    int err = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
    ncols = gtk_tree_model_get_n_columns(model);

    if (USE_RXLIST(sr->ci) && ncols == 4) {
	return get_vecm_exog_list(sr, rows, model);
    }

    gtk_tree_model_get_iter_first(model, &iter);

    if (USE_ZLIST(sr->ci)) {
	reclist = instlist;
    } else if (USE_VECXLIST(sr->ci)) {
	reclist = vecxlist;
    }

    ynum = selector_get_depvar_number(sr);

    for (i=0; i<rows; i++) {

	gtk_tree_model_get(model, &iter, 0, &exog, 1, &lag, -1);

	if (IV_MODEL(sr->ci) && exog == ynum && lag == 0) { /* HECKIT? */
	    errbox("You can't use the dependent variable as an instrument");
	    err = 1;
	    break;
	}
		
	if (is_lag_dummy(exog, lag, context)) {
	    gtk_tree_model_iter_next(model, &iter);
	    continue;
	}

	if (context) {
	    tmp = get_lagpref_string(exog, context, sr);
	} else {
	    tmp = g_strdup_printf(" %d", exog);
	}

	add_to_cmdlist(sr, tmp);
	g_free(tmp);

	if (reclist != NULL) {
	    reclist[j++] = exog;
	}

	gtk_tree_model_iter_next(model, &iter);
    }

    return err;
}

static void read_quantreg_extras (selector *sr)
{
    GtkWidget *e = GTK_BIN(sr->extra[0])->child;
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(e));
	
    if (s == NULL || *s == '\0') {
	warnbox(_("You must specify a quantile"));
	sr->error = 1;
    } else {
	/* convert the GUI string to what ought to be a valid
	   numerical matrix specification */
	gchar *tmp = g_strdup_printf("{%s} ", s);
	gretl_matrix *m;

	comma_separate_numbers(tmp);
	m = generate_matrix(tmp, &Z, datainfo, &sr->error);
	gretl_matrix_free(m);

	if (sr->error) {
	    warnbox(_("Invalid quantile specification"));
	} else {
	    add_to_cmdlist(sr, tmp);
	}
	g_free(tmp);
    }

    if (!sr->error && sr->extra[1] != NULL &&
	GTK_WIDGET_SENSITIVE(sr->extra[1])) {
	GtkAdjustment *adj;

	adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(sr->extra[1]));
	set_optval_double(QUANTREG, OPT_I, adj->value);
    }
} 

static void read_omit_cutoff (selector *sr)
{
    if (sr->extra[0] != NULL && GTK_WIDGET_IS_SENSITIVE(sr->extra[0])) {
	double val;

	val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[0]));
	if (val != 0.10) {
	    set_optval_double(OMIT, OPT_A, val);
	}
    }
}

#define extra_widget_get_int(c) (c == HECKIT || \
                                 c == INTREG || \
                                 c == POISSON || \
                                 c == WLS || \
                                 THREE_VARS_CODE(c))

static void parse_extra_widgets (selector *sr, char *endbit)
{
    const gchar *txt = NULL;
    char numstr[8];
    int k = 0;

    if (sr->ci == QUANTREG) {
	read_quantreg_extras(sr);
	return;
    }

    if (sr->ci == ELLIPSE) {
	read_ellipse_alpha(sr);
	return;
    }

    if (sr->ci == OMIT) {
	read_omit_cutoff(sr);
	return;
    }

    if (sr->ci == WLS || sr->ci == POISSON || 
	sr->ci == AR || sr->ci == HECKIT ||
	sr->ci == INTREG || THREE_VARS_CODE(sr->ci)) {
	txt = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));
	if (txt == NULL || *txt == '\0') {
	    if (sr->ci == WLS) {
		warnbox(_("You must select a weight variable"));
		sr->error = 1;
	    } else if (sr->ci == AR) {
		warnbox(_("You must specify a list of lags"));
		sr->error = 1;
	    } else if (sr->ci == HECKIT) {
		warnbox(_("You must specify a selection variable"));
		sr->error = 1;
	    } else if (sr->ci == INTREG) {
		warnbox(_("You must specify an upper bound variable"));
		sr->error = 1;
	    } else if (sr->ci == ANOVA) {
		warnbox(_("You must specify a treatment variable"));
		sr->error = 1;
	    } else if (THREE_VARS_CODE(sr->ci)) { 
		warnbox(("You must select a Y-axis variable"));
		sr->error = 1;
	    } else if (sr->ci == POISSON) {
		/* the 'extra' field is optional */
		return;
	    }
	}
    } 

    if (sr->error) {
	return;
    }

    if (extra_widget_get_int(sr->ci)) {
	k = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->extra[0]), "data"));
    }

    if (sr->ci == ANOVA) {
	sprintf(numstr, " %d", k);
	add_to_cmdlist(sr, numstr);
    } else if (sr->ci == WLS || THREE_VARS_CODE(sr->ci)) {
	sprintf(numstr, "%d ", k);
	add_to_cmdlist(sr, numstr);
	if (sr->ci == WLS) {
	    wtvar = k;
	}
    } else if (sr->ci == INTREG) {
	sprintf(numstr, " %d ", k);
	add_to_cmdlist(sr, numstr);
	hivar = k;
    } else if (sr->ci == POISSON) {
	sprintf(endbit, " ; %d", k);
	offvar = k;
    } else if (sr->ci == HECKIT) {
	sprintf(endbit, " %d", k);
	selvar = k;
    } else if (sr->ci == AR) {
	free(arlags);
	arlags = gretl_strdup(txt);
	add_to_cmdlist(sr, txt);
	add_to_cmdlist(sr, " ; ");
    } 
}

static char *VAR_dv_lags_string (const int *list, int *err)
{
    char *s = gretl_list_to_lags_string(list, err);
    char *ret = NULL;

    if (s != NULL) {
	ret = malloc(strlen(s) + 9);
	if (ret == NULL) {
	    *err = E_ALLOC;
	} else {
	    sprintf(ret, " --lags=%s", s);
	}
	free(s);
    }

    return ret;
}

static void vec_get_spinner_data (selector *sr, int *order,
				  char **dvlags)
{
    const int *llist;
    int lmax;
    char numstr[8];

    /* lag order from global spinner */
    lmax = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(sr->extra[0]));
    sprintf(numstr, "%d", lmax);
    add_to_cmdlist(sr, numstr);
    *order = lmax;

    /* possible list of specific lags */
    llist = get_VAR_lags_list();

    if (llist != NULL) {
	/* "gappy" lag specification */
	*dvlags = VAR_dv_lags_string(llist, &sr->error);
    } 

    if (sr->ci == VECM) {
	/* cointegration rank */
	jrank = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(sr->extra[1]));
	sprintf(numstr, " %d", jrank);
	add_to_cmdlist(sr, numstr);
    }
}

static void parse_depvar_widget (selector *sr, char *endbit, char **dvlags,
				 char **idvlags)
{
    int ynum = selector_get_depvar_number(sr);

    if (ynum < 0) {
	topslot_empty(sr->ci);
	sr->error = 1;
    } else {
	char numstr[8];

	if (sr->ci == GR_XY || sr->ci == GR_IMP) {
	    sprintf(endbit, " %d", ynum);
	} else {
	    sprintf(numstr, "%d", ynum);
	    add_to_cmdlist(sr, numstr);
	}
	if (select_lags_depvar(sr->ci) && dataset_lags_ok(datainfo)) {
	    *dvlags = get_lagpref_string(ynum, LAG_Y_X, NULL);
	    if (USE_ZLIST(sr->ci)) {
		*idvlags = get_lagpref_string(ynum, LAG_Y_W, NULL);
	    }
	}
	if (sr->default_check != NULL && 
	    gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(sr->default_check))) {
	    default_y = ynum;
	} else if (sr->ci == INTREG) {
	    lovar = ynum;
	}
    }
}

static void parse_third_var_slot (selector *sr)
{
    const gchar *txt = gtk_entry_get_text(GTK_ENTRY(sr->rvars1));
    char numstr[8];

    if (txt == NULL || *txt == '\0') {
	if (sr->ci == ANOVA) {
	    /* third var is optional */
	    return;
	} else if (sr->ci == GR_3D) {
	    warnbox(_("You must select a Z-axis variable"));
	} else if (sr->ci == GR_DUMMY) {
	    warnbox(_("You must select a factor variable"));
	} else {
	    warnbox(_("You must select a control variable"));
	}
	sr->error = 1;
    } else {
	int v = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->rvars1), 
						  "data"));
	sprintf(numstr, " %d", v);
	add_to_cmdlist(sr, numstr);
    }
}

static void selector_cancel_unavailable_options (selector *sr)
{
    if (sr->ci == ARMA) {
	if ((sr->opts & OPT_H) && !GTK_WIDGET_SENSITIVE(sr->hess_button)) {
	    sr->opts ^= OPT_H;
	}
    }
}

static void get_anova_list (selector *sr)
{
    /* get response var */
    parse_depvar_widget(sr, NULL, NULL, NULL);
    if (sr->error < 0) {
	return;
    } 

    /* get treatment var */
    parse_extra_widgets(sr, NULL);
    if (sr->error) {
	return;
    }

    /* get (optional) block var */
    parse_third_var_slot(sr);    
}

/* main function for building a command list from information stored
   in the various selector widgets */

static void compose_cmdlist (selector *sr)
{
    gint rows = 0, realrows = 0;
    char endbit[12] = {0};
    char *dvlags = NULL;
    char *idvlags = NULL;
    int context = 0;
    int order = 0;

    sr->error = 0;
    sr->cmdlist = mymalloc(MAXLEN); 
    if (sr->cmdlist == NULL) {
	return;
    }

    *sr->cmdlist = '\0';

    if (sr->ci == INTREG) {
	parse_depvar_widget(sr, endbit, &dvlags, &idvlags);
	if (!sr->error) {
	    parse_extra_widgets(sr, endbit);
	}
	goto int_next;
    } 

    if (sr->ci == ANOVA) {
	/* special: either 2 or 3 variables selected */
	get_anova_list(sr);
	return;
    }

    /* deal with content of "extra" widgets */
    if (sr->ci == ARMA) {
	arma_spec_to_cmdlist(sr);
    } else if (sr->ci == ARCH || 
	       sr->ci == GARCH || 
	       sr->ci == ARBOND) {
	add_pdq_vals_to_cmdlist(sr);
    } else if (VEC_CODE(sr->ci)) {
	vec_get_spinner_data(sr, &order, &dvlags);
    } else {
	parse_extra_widgets(sr, endbit);
    }

    /* deal with the "depvar" widget */
    if (!sr->error && sr->depvar != NULL) {
	parse_depvar_widget(sr, endbit, &dvlags, &idvlags);
    }

 int_next:

    /* bail out if things have gone wrong already */
    if (sr->error) {
	return;
    }

    if (THREE_VARS_CODE(sr->ci)) { 
	parse_third_var_slot(sr);
	return;
    } 

    /* count the rows (variables) in the "primary" right-hand selection
       list box */
    rows = varlist_row_count(sr, SR_RVARS1, &realrows);

    if (sr->ci == SCATTERS) {
	if (rows > 0) {
	    add_to_cmdlist(sr, " ;");
	} else {
	    sr->error = E_ARGS;
	    gui_errmsg(sr->error);
	    return;
	}
    }    

    if ((sr->ci == COINT || sr->ci == COINT2 || sr->ci == VECM) && rows < 2) {
	warnbox(_("You must select two or more endogenous variables"));
	sr->error = 1;
	return;
    } else if ((sr->ci == VAR || sr->ci == VLAGSEL) && rows < 1) {
	warnbox(_("You must select a dependent variable"));
	sr->error = 1;
	return;
    }	

    if (realrows > 0) {
	maybe_resize_recorder_lists(sr, realrows);
    }

    /* primary RHS varlist */
    context = sr_get_lag_context(sr, SR_RVARS1);
    get_rvars1_data(sr, rows, context);

    /* cases with a (possibly optional) secondary RHS list */
    if (USE_ZLIST(sr->ci) || USE_VECXLIST(sr->ci)) {
	rows = varlist_row_count(sr, SR_RVARS2, &realrows);
	if (rows > 0) {
	    if (realrows > 0) {
		maybe_resize_exog_recorder_lists(sr, realrows);
	    }
	    context = sr_get_lag_context(sr, SR_RVARS2);

	    if (USE_ZLIST(sr->ci) && dvlags != NULL) {
		add_to_cmdlist(sr, dvlags);
		free(dvlags);
		dvlags = NULL;
	    }

	    if (sr->ci == HECKIT) {
		add_to_cmdlist(sr, " ;");
		add_to_cmdlist(sr, endbit);
		*endbit = '\0';
	    } else if (*sr->cmdlist != '\0') {
		add_to_cmdlist(sr, " ;");
	    }

	    sr->error = get_rvars2_data(sr, rows, context);
	} else if (IV_MODEL(sr->ci)) {
	    warnbox(_("You must specify a set of instrumental variables"));
	    sr->error = 1;
	} else if (sr->ci == HECKIT) {
	    warnbox(_("You must specify regressors for the selection equation"));
	    sr->error = 1;
	} 
    }

    /* deal with any trailing strings */
    if (!sr->error) {
	if (endbit[0] != '\0') {
	    add_to_cmdlist(sr, endbit);
	} else if (dvlags != NULL) {
	    add_to_cmdlist(sr, dvlags);
	    free(dvlags);
	} else if (idvlags != NULL) {
	    add_to_cmdlist(sr, idvlags);
	    free(idvlags);
	}
    }

    if ((sr->ci == SCATTERS) && !sr->error) {
	int xstate;

	xstate = gtk_combo_box_get_active(GTK_COMBO_BOX(multiplot_menu));
	if (xstate) {
	    reverse_list(sr->cmdlist);
	}
    }

#if 0
    fprintf(stderr, "sr->cmdlist:\n'%s'\n", sr->cmdlist);
#endif

    if (!sr->error) {
	/* record some choices as defaults */
	if (sr->ci == VECM || sr->ci == VAR || sr->ci == VLAGSEL) {
	    want_seasonals = (sr->opts & OPT_D)? 1 : 0;
	}
	if (sr->ci == VECM || sr->ci == VAR || sr->ci == ARCH) {
	    default_order = (order < 0)? -order : order;
	}
	if (sr->ci == VAR || sr->ci == VLAGSEL) {
	    vartrend = (sr->opts & OPT_T)? 1 : 0;
	}
	if (sr->ci == ARMA) {
	    arma_const = (sr->opts & OPT_N)? 0 : 1;
	    arma_hessian = (sr->opts & OPT_G)? 0 : 1;
	    arma_x12 = (sr->opts & OPT_X)? 1 : 0;
	}
	if (sr->ci == GARCH) {
	    if (sr->opts & OPT_F) {
		sr->opts &= ~OPT_F;
		libset_set_bool(USE_FCP, 1);
	    } else {
		libset_set_bool(USE_FCP, 0);
	    }
	}
	verbose = (sr->opts & OPT_V)? 1 : 0;
    }

    selector_cancel_unavailable_options(sr);
}

static void cancel_selector (GtkWidget *widget, selector *sr)
{
    if (open_selector != NULL) {
	gtk_widget_destroy(sr->dlg);
    }
}

static void destroy_selector (GtkWidget *w, selector *sr) 
{
    if (SAVE_DATA_ACTION(sr->ci) || sr->ci == SAVE_FUNCTIONS ||
	sr->ci == DEFINE_MATRIX) {
	gtk_main_quit();
    }

    free(sr->cmdlist);
    free(sr);

    open_selector = NULL;
}

static char *est_str (int cmdnum)
{
    switch (cmdnum) {
    case OLS:
	return N_("OLS");
    case HSK:
	return N_("Heteroskedasticity corrected");
    case CORC:
	return N_("Cochrane-Orcutt");
    case HILU:
	return N_("Hildreth-Lu");
    case PWE:
	return N_("Prais-Winsten");
    case LOGIT:
	return N_("Logit");
    case OLOGIT:
	return N_("Ordered Logit");
    case MLOGIT:
	return N_("Multinomial Logit");
    case PROBIT:
	return N_("Probit");
    case OPROBIT:
	return N_("Ordered Probit");
    case TOBIT:
	return N_("Tobit");
    case HECKIT:
	return N_("Heckit");
    case LOGISTIC:
	return N_("Logistic");
    case POISSON:
	return N_("Poisson");
    case PANEL:
	return N_("Panel model");
    case PANEL_WLS:
	return N_("Groupwise WLS");
    case PANEL_B:
	return N_("Between-groups model");
    case ARBOND:
	return N_("Dynamic panel model");
    case WLS:
	return N_("Weighted least squares");
    case IVREG:
	return N_("Two-stage least squares");
    case IV_LIML:
	return N_("Limited information maximum likelihood");
    case IV_GMM:
	return N_("Generalized method of moments");
    case AR:
	return N_("Autoregressive model");
    case ARMA:
	return N_("ARIMA");
    case ARCH:
	return N_("ARCH");
    case GARCH:
	return N_("GARCH");
    case VAR:
	return N_("VAR");
    case VLAGSEL:
	return N_("VAR lag selection");
    case VECM:
	return N_("VECM");
    case LAD:
	return N_("LAD");
    case QUANTREG:
	return N_("Quantile regression");
    case INTREG:
	return N_("Interval regression");
    case COINT:
    case COINT2:
	return N_("Cointegration");
#ifdef ENABLE_GMP
    case MPOLS:
	return N_("Multiple precision OLS");
#endif
    default:
	return "";
    }
}

static char *extra_string (int ci)
{
    switch (ci) {
    case WLS:
	return N_("Weight variable");
    case POISSON:
	return N_("Offset variable");
    case HECKIT:
	return N_("Selection variable");
    case AR:
	return N_("List of AR lags");
    case QUANTREG:
	return N_("Desired quantile(s)");
    case INTREG:
	return N_("Upper bound variable");
    case GR_DUMMY:
    case GR_3D:
    case GR_XYZ:	
	return N_("Y-axis variable");
    case ANOVA:
	return N_("Treatment variable");
    default:
	return NULL;
    }
}

static gint flip_multiplot_axis (GtkComboBox *box, gpointer p)
{
    gint xstate = gtk_combo_box_get_active(box);

    if (xstate) {
	gtk_label_set_text(GTK_LABEL(multiplot_label), _("Y-axis variables"));
    } else {
	gtk_label_set_text(GTK_LABEL(multiplot_label), _("X-axis variables"));
    }

    return FALSE;
}

static GtkWidget *multiplot_popdown (int ci)
{
    GtkWidget *w = gtk_combo_box_new_text();

    gtk_combo_box_append_text(GTK_COMBO_BOX(w), _("Y-axis variable"));
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), _("X-axis variable"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
		     G_CALLBACK(flip_multiplot_axis), NULL);

    multiplot_menu = w;

    return w;
}

static gint set_gmm_est_option (GtkComboBox *box, selector *sr)
{
    gint i = gtk_combo_box_get_active(box);

    if (i == 0) {
	/* 1-step */
	sr->opts &= ~(OPT_I | OPT_T);
    } else if (i == 1) {
	/* 2-step */
	sr->opts |= OPT_T;
    } else {
	/* iterated */
	sr->opts |= OPT_I;
    }

    return FALSE;
}

static void build_gmm_popdown (selector *sr)
{
    GtkWidget *w = gtk_combo_box_new_text();
    GtkWidget *hbox;

    gtk_combo_box_append_text(GTK_COMBO_BOX(w), _("One-step estimation"));
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), _("Two-step estimation"));
    gtk_combo_box_append_text(GTK_COMBO_BOX(w), _("Iterated estimation"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
		     G_CALLBACK(set_gmm_est_option), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_widget_show(w);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);    
}

static GtkWidget *
entry_with_label_and_chooser (selector *sr, GtkWidget *vbox,
			      gchar *label_string,
			      int label_active,
			      void (*clickfunc)())
{
    GtkWidget *tmp, *x_hbox;
    GtkWidget *entry;

    if (label_active) {
	tmp = multiplot_popdown(sr->ci);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show_all(tmp);
    } else if (label_string != NULL) {
	tmp = gtk_label_new(label_string);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }

    x_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose->"));
    gtk_box_pack_start(GTK_BOX(x_hbox), tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(clickfunc), sr);
    gtk_widget_show(tmp); 

    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN - 1);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN + 3);

    gtk_box_pack_start(GTK_BOX(x_hbox), entry, FALSE, FALSE, 0);
    gtk_widget_show(entry); 

    gtk_box_pack_start(GTK_BOX(vbox), x_hbox, FALSE, FALSE, 0);
    gtk_widget_show(x_hbox); 

    if (label_active || label_string != NULL) {
	if (clickfunc != set_third_var_callback) {
	    vbox_add_hsep(vbox);
	}
    }

    return entry;
}

static void build_x_axis_section (selector *sr, GtkWidget *right_vbox)
{
    if (sr->ci == SCATTERS) {
	sr->depvar = entry_with_label_and_chooser(sr, right_vbox,
						  NULL, 1,
						  set_dependent_var_callback);
    } else {
	sr->depvar = entry_with_label_and_chooser(sr, right_vbox,
						  _("X-axis variable"), 0,
						  set_dependent_var_callback);
    }
}

static void maybe_activate_depvar_lags (GtkWidget *w, selector *sr)
{
    if (select_lags_depvar(sr->ci) && sr->lags_button != NULL) {
	const gchar *txt = gtk_entry_get_text(GTK_ENTRY(w));

	if (txt != NULL && *txt != 0) {
	    gtk_widget_set_sensitive(sr->lags_button, TRUE);
	}
    }
}

/* returns ID number of variable pre-inserted as dependent, or -1 */

static int build_depvar_section (selector *sr, GtkWidget *right_vbox,
				 int preselect)
{
    GtkWidget *tmp, *depvar_hbox;
    int defvar;

    if (sr->ci == INTREG) {
	defvar = (lovar > 0 && lovar < datainfo->v)? lovar : -1;
    } else {
	if (default_y >= datainfo->v) {
	    default_y = -1;
	}
	defvar = (preselect)? preselect : default_y;
    }

    if (sr->ci == INTREG) {
	tmp = gtk_label_new(_("Lower bound variable"));
    } else if (sr->ci == ANOVA) {
	tmp = gtk_label_new(_("Response variable"));
    } else {
	tmp = gtk_label_new(_("Dependent variable"));
    }

    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);

    depvar_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label (_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(depvar_hbox), tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
                      G_CALLBACK(set_dependent_var_callback), sr);
    gtk_widget_show(tmp); 

    sr->depvar = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(sr->depvar), VNAMELEN - 1);
    gtk_entry_set_width_chars(GTK_ENTRY(sr->depvar), VNAMELEN + 3);

    g_signal_connect(G_OBJECT(sr->depvar), "changed",
		     G_CALLBACK(maybe_activate_depvar_lags), sr);

    if (defvar >= 0) {
        gtk_entry_set_text(GTK_ENTRY(sr->depvar), datainfo->varname[defvar]);
    }

    gtk_box_pack_start(GTK_BOX(depvar_hbox), sr->depvar, FALSE, FALSE, 0);
    gtk_widget_show(sr->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), depvar_hbox, FALSE, FALSE, 0);
    gtk_widget_show(depvar_hbox); 

    if (sr->ci != INTREG) {
	sr->default_check = gtk_check_button_new_with_label(_("Set as default"));
	gtk_box_pack_start(GTK_BOX(right_vbox), sr->default_check, FALSE, FALSE, 0);
	gtk_widget_show(sr->default_check); 
    } 

    if (sr->ci != ARBOND) {
	vbox_add_hsep(right_vbox);
    }

    return defvar;
}

static void 
build_public_iface_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp, *pub_hbox;

    tmp = gtk_label_new(_("Public interface"));
    gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
    gtk_widget_show(tmp);

    pub_hbox = gtk_hbox_new(FALSE, 5); 

    tmp = gtk_button_new_with_label(_("Choose ->"));
    gtk_box_pack_start(GTK_BOX(pub_hbox), tmp, TRUE, TRUE, 0);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
                      G_CALLBACK(set_dependent_var_callback), sr);
    gtk_widget_show(tmp); 

    sr->depvar = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(sr->depvar), FN_NAMELEN - 1);
    gtk_entry_set_width_chars(GTK_ENTRY(sr->depvar), 18);

    gtk_box_pack_start(GTK_BOX(pub_hbox), sr->depvar, FALSE, FALSE, 0);
    gtk_widget_show(sr->depvar); 

    gtk_box_pack_start(GTK_BOX(right_vbox), pub_hbox, FALSE, FALSE, 0);
    gtk_widget_show(pub_hbox); 

    vbox_add_hsep(right_vbox);
}

/* In case we have a saved preference for the max lag of the
   endogenous vars in a VAR, set via the lags dialog, update this
   value from the global spinner 
*/

static void lag_order_sync (GtkSpinButton *b, selector *sr)
{
    if (GTK_WIDGET_IS_SENSITIVE(b) && sr->extra[EXTRA_LAGS] != NULL) {
	int lmax = gtk_spin_button_get_value_as_int(b);

	set_VAR_max_lag(lmax);
    }
}

enum {
    LAG_ONLY,
    LAG_AND_RANK
};

static void lag_order_spin (selector *sr, GtkWidget *vbox, int which)
{
    GtkWidget *tmp, *hbox;
    GtkObject *adj;
    gdouble lag; 
    gdouble minlag;
    gdouble maxlag;
    const char *labels[] = {
	N_("lag order:"),
	N_("cointegration rank:")
    };
    int i, nspin = (which == LAG_AND_RANK)? 2 : 1;

    maxlag = (datainfo->n < 72)? (datainfo->n / 2) : 36;
    minlag = (sr->ci == COINT)? 0 : 1;

    if (default_order > 0 && default_order <= maxlag) {
	lag = default_order;
    } else {
	lag = (datainfo->pd > 12)? 12 : datainfo->pd;
    }

    if (sr->ci == VLAGSEL) {
	minlag = 2;
	lag *= 2;
	if (lag > maxlag) {
	    lag = maxlag;
	}
    }

    for (i=0; i<nspin; i++) {
	hbox = gtk_hbox_new(FALSE, 5);

	if (sr->ci == VLAGSEL) {
	    tmp = gtk_label_new(_("maximum lag:"));
	} else {
	    tmp = gtk_label_new(_(labels[i]));
	}

	gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
	gtk_widget_show(tmp);
	gtk_misc_set_alignment(GTK_MISC(tmp), 0.0, 0.5);
	
	if (i == 0) {
	    /* lag order */
	    adj = gtk_adjustment_new(lag, minlag, maxlag, 1, 1, 0);
	} else {
	    /* rank */
	    adj = gtk_adjustment_new(jrank, 1, 10, 1, 1, 0);
	}

	sr->extra[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_box_pack_start(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
	gtk_widget_show(sr->extra[i]);

	if (i == 0) {
	    /* cross-connect with lag preferences dialog */
	    if (get_VAR_lags_list() != NULL) {
		gtk_widget_set_sensitive(sr->extra[i], FALSE);
	    }
	    g_signal_connect(G_OBJECT(sr->extra[i]), "value-changed",
			     G_CALLBACK(lag_order_sync), sr);
	}

	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox); 
    }
}

static void AR_order_spin (selector *sr, GtkWidget *vbox)
{
    GtkWidget *tmp, *hbox;
    GtkObject *adj;
    gdouble val, maxlag;

    hbox = gtk_hbox_new(FALSE, 5);

    if (sr->ci == ARCH) {
	tmp = gtk_label_new(_("ARCH order:"));
	val = datainfo->pd;
	maxlag = 2 * datainfo->pd;
    } else {
	/* arbond */
	tmp = gtk_label_new(_("AR order:"));
	val = 1;
	maxlag = 10;
	if (maxlag < datainfo->pd) {
	    maxlag = datainfo->pd;
	}
    }

    if (default_order > 0 && default_order <= maxlag) {
	val = default_order;
    }

    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);
    gtk_misc_set_alignment(GTK_MISC(tmp), 0.0, 0.5);
    adj = gtk_adjustment_new(val, 1, maxlag, 1, 1, 0);

    sr->extra[0] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_box_pack_start(GTK_BOX(hbox), sr->extra[0], FALSE, FALSE, 5);
    gtk_widget_show(sr->extra[0]);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox); 
}

static void third_var_box (selector *sr, GtkWidget *vbox)
{
    const gchar *label;

    if (sr->ci == GR_3D) {
	label = N_("Z-axis variable");
    } else if (sr->ci == GR_DUMMY) {
	label = _("Factor (dummy)");
    } else if (sr->ci == GR_XYZ) {
	label = N_("Control variable");
    } else if (sr->ci == ANOVA) {
	label = N_("Block variable (optional)");
    } else {
	return;
    }

    sr->rvars1 = entry_with_label_and_chooser(sr, vbox, _(label), 0,
					      set_third_var_callback);
}

static void extra_var_box (selector *sr, GtkWidget *vbox)
{
    int setvar = 0;

    sr->extra[0] = entry_with_label_and_chooser(sr, vbox,
						NULL, 0,
						set_extra_var_callback);

    if (sr->ci == WLS && wtvar > 0 && wtvar < datainfo->v) {
	setvar = wtvar;
    } else if (sr->ci == HECKIT && selvar > 0 && selvar < datainfo->v) {
	setvar = selvar;
    } else if (sr->ci == INTREG && hivar > 0 && hivar < datainfo->v) {
	setvar = hivar;
    } else if (sr->ci == POISSON && offvar > 0 && offvar < datainfo->v) {
	setvar = offvar;
    }

    if (setvar > 0) {
	gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), datainfo->varname[setvar]);
	g_object_set_data(G_OBJECT(sr->extra[0]), "data",
			  GINT_TO_POINTER(setvar));
    }	
}

static gboolean accept_right_arrow (GtkWidget *w, GdkEventKey *key,
				    gpointer p)
{
    if (key->keyval == GDK_Right) { 
        gtk_button_clicked(GTK_BUTTON(p));
	return TRUE;
    } else {
	return FALSE;
    }
}

static GtkWidget *add_button (selector *sr, GtkWidget *box)
{
    GtkWidget *w = gtk_button_new_with_label(_("Add ->"));

    g_signal_connect(G_OBJECT(sr->dlg), "key-press-event", 
		     G_CALLBACK(accept_right_arrow), w);
    gtk_box_pack_start(GTK_BOX(box), w, TRUE, FALSE, 0);
    return w;
}

static gboolean accept_left_arrow (GtkWidget *w, GdkEventKey *key,
				   gpointer p)
{
    if (key->keyval == GDK_Left) { 
        gtk_button_clicked(GTK_BUTTON(p));
	return TRUE;
    } else {
	return FALSE;
    }
}

static GtkWidget *remove_button (selector *sr, GtkWidget *box)
{
    GtkWidget *w = gtk_button_new_with_label(_("<- Remove"));

    g_signal_connect(G_OBJECT(sr->dlg), "key-press-event", 
		     G_CALLBACK(accept_left_arrow), w);
    gtk_box_pack_start(GTK_BOX(box), w, TRUE, FALSE, 0);
    return w;
}

static void auxiliary_rhs_varlist (selector *sr, GtkWidget *vbox)
{
    GtkTreeModel *mod;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *remove, *hbox, *bvbox;
    GtkWidget *tmp = NULL;
    int i;

    if (USE_VECXLIST(sr->ci)) {
	tmp = gtk_label_new(_("Exogenous variables"));
    } else if (IV_MODEL(sr->ci)) {
	tmp = gtk_label_new(_("Instruments"));
    } else if (sr->ci == HECKIT) {
	tmp = gtk_label_new(_("Selection equation regressors"));
    }

    if (tmp != NULL) {
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }

    hbox = gtk_hbox_new(FALSE, 5);
    bvbox = gtk_vbox_new(TRUE, 0);

    tmp = add_button(sr, bvbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(add_to_rvars2_callback), sr);
    
    remove = remove_button(sr, bvbox);

    if (sr->ci == VAR || sr->ci == VECM || sr->ci == VLAGSEL) {
	/* select lags of exogenous variables */
	sr->lags_button = gtk_button_new_with_label(_("lags..."));
	gtk_box_pack_start(GTK_BOX(bvbox), sr->lags_button, TRUE, FALSE, 0);
	g_signal_connect(G_OBJECT(sr->lags_button), "clicked", 
			 G_CALLBACK(lags_dialog_driver), sr);
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }

    gtk_box_pack_start(GTK_BOX(hbox), bvbox, TRUE, TRUE, 0);
    gtk_widget_show_all(bvbox);

    /* then the listing */
    sr->rvars2 = var_list_box_new(GTK_BOX(hbox), sr, SR_RVARS2);
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));

    store = GTK_LIST_STORE(mod);
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(mod, &iter);

    if (USE_ZLIST(sr->ci) && instlist != NULL) {
	for (i=1; i<=instlist[0]; i++) {
	    if (instlist[i] < datainfo->v) {
		list_append_var(mod, &iter, instlist[i], sr, SR_RVARS2);
	    }
	}
    } else if (USE_VECXLIST(sr->ci) && vecxlist != NULL) {
	set_varflag(UNRESTRICTED);
	for (i=1; i<=vecxlist[0]; i++) {
	    if (vecxlist[i] == LISTSEP) {
		set_varflag(RESTRICTED);
	    } else if (vecxlist[i] > 0) {
		list_append_var(mod, &iter, vecxlist[i], sr, SR_RVARS2);
		if (sr->lags_button != NULL) {
		    gtk_widget_set_sensitive(sr->lags_button, TRUE);
		}
	    }
	}
	set_varflag(UNRESTRICTED);
    } else if (!VEC_CODE(sr->ci)) {
	list_append_var(mod, &iter, 0, sr, SR_RVARS2);
    }

    /* hook up remove button to list box */
    g_signal_connect(G_OBJECT(remove), "clicked", 
		     G_CALLBACK(remove_from_right_callback), 
		     sr->rvars2);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    gtk_widget_show(hbox); 
}

static void make_tau_list (GtkWidget *w)
{
    GtkComboBox *box = GTK_COMBO_BOX(w);

    gtk_combo_box_append_text(box, "0.25 0.50 0.75");
    gtk_combo_box_append_text(box, ".05, .25 .50 .75, .95");
    gtk_combo_box_append_text(box, ".1 .2 .3 .4 .5 .6 .7 .8 .9");
} 

static int maybe_set_entry_text (GtkWidget *w, const char *s)
{
    if (s != NULL && *s != '\0') {
	gtk_entry_set_text(GTK_ENTRY(w), s + (*s == ' '));
	return 1;
    } else {
	return 0;
    }
}

static void build_mid_section (selector *sr, GtkWidget *right_vbox)
{
    GtkWidget *tmp;
    const char *str = _(extra_string(sr->ci));

    if (str != NULL) {
	tmp = gtk_label_new(str);
	gtk_box_pack_start(GTK_BOX(right_vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }	

    if (sr->ci == HECKIT) {
	extra_var_box(sr, right_vbox);
	vbox_add_hsep(right_vbox);
	primary_rhs_varlist(sr, right_vbox);
    } else if (sr->ci == WLS || sr->ci == INTREG || 
	       sr->ci == POISSON || THREE_VARS_CODE(sr->ci)) {
	extra_var_box(sr, right_vbox);
    } else if (USE_ZLIST(sr->ci)) {
	primary_rhs_varlist(sr, right_vbox);
    } else if (sr->ci == AR) {
	sr->extra[0] = gtk_entry_new();
	gtk_box_pack_start(GTK_BOX(right_vbox), sr->extra[0], 
			   FALSE, TRUE, 0);
	maybe_set_entry_text(sr->extra[0], arlags);
	gtk_widget_show(sr->extra[0]); 
    } else if (sr->ci == QUANTREG) {
	sr->extra[0] = gtk_combo_box_entry_new_text();
	make_tau_list(sr->extra[0]);
	gtk_entry_set_text(GTK_ENTRY(GTK_BIN(sr->extra[0])->child), "0.5");
	gtk_box_pack_start(GTK_BOX(right_vbox), sr->extra[0], 
			   FALSE, TRUE, 0);
	gtk_widget_show(sr->extra[0]); 
    } else if (sr->ci == VAR || sr->ci == VLAGSEL) {
	lag_order_spin(sr, right_vbox, LAG_ONLY);
	vbox_add_hsep(right_vbox);
	primary_rhs_varlist(sr, right_vbox);
    } else if (sr->ci == VECM) {
	lag_order_spin(sr, right_vbox, LAG_AND_RANK);
	vbox_add_hsep(right_vbox);
	primary_rhs_varlist(sr, right_vbox);
    } else if (sr->ci == COINT2) {
	lag_order_spin(sr, right_vbox, LAG_ONLY);
	vbox_add_hsep(right_vbox);
	primary_rhs_varlist(sr, right_vbox);
    } else if (VEC_CODE(sr->ci)) {
	lag_order_spin(sr, right_vbox, LAG_ONLY);
    } else if (sr->ci == ARBOND || sr->ci == ARCH) {
	AR_order_spin(sr, right_vbox);
    }
    
    vbox_add_hsep(right_vbox);
}

static void selector_init (selector *sr, guint ci, const char *title,
			   gpointer p, int (*callback)())
{
    GtkWidget *base;
    double x;
    int dlgx = -1, dlgy = 340;
    int i;

    sr->ci = ci;
    sr->opts = (ci == PANEL_WLS)? OPT_W : OPT_NONE;
    sr->data = p;
    
    if (MODEL_CODE(ci)) {
	if (datainfo->v > 9) {
	    dlgy += 80;
	} 
    } 

    if (ci == ARMA) {
	dlgy += 80;
    } else if (ci == WLS || ci == INTREG || ci == POISSON || ci == AR) {
	dlgy += 30;
    } else if (ci == HECKIT) {
	dlgy += 80;
    } else if (IV_MODEL(ci)) {
	dlgy += 60;
    } else if (ci == ANOVA) {
	dlgy -= 60;
    } else if (VEC_CODE(ci)) {
	dlgy = 450;
	if (ci == VAR || ci == VECM) {
	    dlgy += 50;
	} else if (ci == VLAGSEL) {
	    dlgy += 40;
	}
    } 

    if (WANT_TOGGLES(ci) && ci != COINT) {
	dlgy += 40;
    }

    if (want_radios(sr)) {
	dlgy += 60;
    }

    if (want_combo(sr)) {
	dlgy += 20;
    }    

    if (ci == ARMA && datainfo->pd > 1) {
	/* seasonal spinners */
	dlgy += 60;
    }

    if (ci == GARCH) {
	/* extra check boxes */
	dlgy += 50;
    }

    if (dataset_lags_ok(datainfo)) {
	if (MODEL_CODE(ci) && ci != ARMA) {
	    /* lag selector button at foot */
	    dlgy += 30;
	}
    }

    sr->lvars = NULL;
    sr->depvar = NULL;
    sr->rvars1 = NULL;
    sr->rvars2 = NULL;
    sr->default_check = NULL;
    sr->add_button = NULL;
    sr->remove_button = NULL;
    sr->lags_button = NULL;
    sr->hess_button = NULL;
    sr->x12a_button = NULL;

    for (i=0; i<N_EXTRA; i++) {
	sr->extra[i] = NULL;
    }

    for (i=0; i<N_RADIOS; i++) {
	sr->radios[i] = NULL;
    }    

    sr->cmdlist = NULL;
    sr->callback = callback;

    sr->active_var = 0;
    sr->error = 0;
    sr->n_left = 0;

    sr->dlg = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    open_selector = sr;

    gtk_window_set_title(GTK_WINDOW(sr->dlg), title);

    x = (double) dlgy * gui_scale;
    dlgy = x;

    if (ci == SAVE_FUNCTIONS) {
	x = (double) 460 * gui_scale;
	dlgx = x;
    }
    
    gtk_window_set_default_size(GTK_WINDOW(sr->dlg), dlgx, dlgy); 

    g_signal_connect(G_OBJECT(sr->dlg), "destroy", 
		     G_CALLBACK(destroy_selector), 
		     sr);
    g_signal_connect(G_OBJECT(sr->dlg), "key-press-event", 
		     G_CALLBACK(esc_kills_window), NULL);

    /* create equivalent of gtkdialog structure */
    base = gtk_vbox_new(FALSE, 5);
    gtk_container_add(GTK_CONTAINER(sr->dlg), base);
    gtk_widget_show(base);

    sr->vbox = gtk_vbox_new(FALSE, 0);
    gtk_widget_show(sr->vbox);

    /* make (upper) vbox expansible */
    gtk_box_pack_start(GTK_BOX(base), sr->vbox, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(sr->vbox), 5);
    gtk_box_set_spacing(GTK_BOX(sr->vbox), 5);

    vbox_add_hsep(base);

    sr->action_area = gtk_hbutton_box_new();
    gtk_button_box_set_layout(GTK_BUTTON_BOX(sr->action_area), 
			      GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(sr->action_area), 10);
    gtk_widget_show(sr->action_area);
    gtk_box_pack_start(GTK_BOX(base), sr->action_area,
		       FALSE, FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(sr->action_area), 5);
} 

static void option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	sr->opts |= opt;
    } else {
	sr->opts &= ~opt;
    }  

    if (sr->ci == PANEL) {
	GtkWidget *w = g_object_get_data(G_OBJECT(sr->dlg), "robust-button");
	
	if (w != NULL) {
	    gtk_widget_set_sensitive(w, !(sr->opts & OPT_U));
	}
    }
}

static void reverse_option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	sr->opts &= ~opt;
    } else {
	sr->opts |= opt;
    }    
}

static void garch_spin_check (GtkSpinButton *b, selector *sr)
{
    int p = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[0]));
    int q = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[1]));
    int i = (GTK_WIDGET(b) == sr->extra[1])? 1 : 0;

    if (p + q > 5) {
	/* limit p + q to 5 */
	if (i == 0) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[1]),
				      (gdouble) --q);
	} else {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[0]),
				      (gdouble) --p);
	}
    } else if (p > 0 && q == 0) {
	/* rule out pure AR in variance */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[1]), 1);
    }
}

static void build_garch_spinners (selector *sr)
{
    GtkWidget *tmp, *hbox;
    GtkObject *adj;
    gdouble val;
    const char *strs[] = {
	N_("GARCH p:"),
	N_("ARCH q:")
    };    
    int i;

    hbox = gtk_hbox_new(FALSE, 5);

    for (i=0; i<2; i++) {
	tmp = gtk_label_new(_(strs[i]));
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 5);
	val = (i==0)? garch_p : garch_q;
	adj = gtk_adjustment_new(val, 0, 4, 1, 1, 0);
	sr->extra[i] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	gtk_box_pack_start(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
	g_signal_connect(GTK_SPIN_BUTTON(sr->extra[i]), "value-changed",
			 G_CALLBACK(garch_spin_check), sr);
    }

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show_all(hbox);
}

static GtkWidget *arma_aux_label (int i)
{
    GtkWidget *hbox;
    GtkWidget *lbl;
    const char *strs[] = {
	N_("Non-seasonal"),
	N_("Seasonal")
    };

    hbox = gtk_hbox_new(FALSE, 5);
    lbl = gtk_label_new(_(strs[i]));
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);
    gtk_widget_show(lbl);

    return hbox;
}

static void toggle_p (GtkWidget *w, selector *sr)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	gtk_widget_set_sensitive(sr->extra[0], FALSE);
	gtk_widget_set_sensitive(sr->extra[1], TRUE);
    } else {
	gtk_widget_set_sensitive(sr->extra[0], TRUE);
	gtk_widget_set_sensitive(sr->extra[1], FALSE);
    }	
}

static void toggle_q (GtkWidget *w, selector *sr)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	gtk_widget_set_sensitive(sr->extra[3], FALSE);
	gtk_widget_set_sensitive(sr->extra[4], TRUE);
    } else {
	gtk_widget_set_sensitive(sr->extra[3], TRUE);
	gtk_widget_set_sensitive(sr->extra[4], FALSE);
    }	
}

static void build_arma_spinners (selector *sr)
{
    GtkWidget *lbl, *chk, *tab;
    GtkWidget *hbox;
    GtkObject *adj;
    gdouble vmax, val;
    gboolean freeform;
    const char *strs[] = {
	N_("AR order:"),
	N_("Difference:"),
	N_("MA order:")
    };
    int i, j;

    if (datainfo->pd > 1) {
	lbl = arma_aux_label(0);
	gtk_box_pack_start(GTK_BOX(sr->vbox), lbl, FALSE, FALSE, 0);
	gtk_widget_show(lbl);
    }

    tab = gtk_table_new(3, 4, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tab), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tab), 5);
    gtk_box_pack_start(GTK_BOX(sr->vbox), tab, FALSE, FALSE, 0);

    /* AR lags */
    lbl = gtk_label_new(_(strs[0]));
    gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 0, 1);
    adj = gtk_adjustment_new(arma_p, 0, 10, 1, 1, 0);
    sr->extra[0] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[0], 1, 2, 0, 1);
    chk = gtk_check_button_new_with_label(_("or specific lags"));
    g_signal_connect(G_OBJECT(chk), "clicked", G_CALLBACK(toggle_p), sr);
    gtk_table_attach_defaults(GTK_TABLE(tab), chk, 2, 3, 0, 1);
    /* or free-form lags */
    sr->extra[1] = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(sr->extra[1]), 16);
    freeform = maybe_set_entry_text(sr->extra[1], arlags);
    gtk_widget_set_sensitive(sr->extra[1], FALSE);
    gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[1], 3, 4, 0, 1);
    if (freeform) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), TRUE);
    }

    /* order for differencing */
    lbl = gtk_label_new(_(strs[1]));
    gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 1, 2);
    adj = gtk_adjustment_new(arima_d, 0, 2, 1, 1, 0);
    sr->extra[2] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[2], 1, 2, 1, 2);

    /* MA lags */
    lbl = gtk_label_new(_(strs[2]));
    gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 2, 3);
    adj = gtk_adjustment_new(arma_q, 0, 10, 1, 1, 0);
    sr->extra[3] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
    gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[3], 1, 2, 2, 3);
    chk = gtk_check_button_new_with_label(_("or specific lags"));
    g_signal_connect(G_OBJECT(chk), "clicked", G_CALLBACK(toggle_q), sr);
    gtk_table_attach_defaults(GTK_TABLE(tab), chk, 2, 3, 2, 3);
    /* or free-form lags */
    sr->extra[4] = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(sr->extra[4]), 16);
    freeform = maybe_set_entry_text(sr->extra[1], malags);
    gtk_widget_set_sensitive(sr->extra[4], FALSE);
    gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[4], 3, 4, 2, 3);
    if (freeform) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), TRUE);
    }

    gtk_widget_show_all(tab);

    j = 5;

    if (datainfo->pd > 1) {
	vbox_add_hsep(sr->vbox);

	lbl = arma_aux_label(1);
	gtk_box_pack_start(GTK_BOX(sr->vbox), lbl, FALSE, FALSE, 0);
	gtk_widget_show(lbl);

	hbox = gtk_hbox_new(FALSE, 5);
	for (i=0; i<3; i++) {
	    lbl = gtk_label_new(_(strs[i]));
	    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 0);
	    val = (i==0)? arma_P : (i==1)? arima_D : arma_Q;
	    vmax = (i == 1)? 2 : 4;
	    adj = gtk_adjustment_new(val, 0, vmax, 1, 1, 0);
	    sr->extra[j] = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 0);
	    gtk_box_pack_start(GTK_BOX(hbox), sr->extra[j++], FALSE, FALSE, 5);
	}

	gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
	gtk_widget_show_all(hbox);    
    }
}

static void hc_config (GtkWidget *w, selector *sr)
{
    options_dialog(TAB_VCV, NULL, sr->dlg);
}

static void pack_switch (GtkWidget *b, selector *sr,
			 gboolean checked, gboolean reversed, 
			 gretlopt opt, int child)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    gint offset = (child)? 15 : 0;
    gint i = opt;

    g_object_set_data(G_OBJECT(b), "opt", GINT_TO_POINTER(i));

    if (reversed) {
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(reverse_option_callback), sr);
	if (checked) {
	    sr->opts &= ~opt;
	} else {
	    sr->opts |= opt;
	}
    } else {
	g_signal_connect(G_OBJECT(b), "toggled", 
			 G_CALLBACK(option_callback), sr);
	if (checked) {
	    sr->opts |= opt;
	} else {
	    sr->opts &= ~opt;
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, offset);
    gtk_widget_show(b);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), checked);
}

#if 0 /* not ready yet */

static void call_iters_dialog (GtkWidget *w, selector *sr)
{
    iterinfo iinfo;
    int cancel = 0;

    iinfo.ci = sr->ci;

    /* FIXME need to figure out properly where to get the
       defaults from, given the command and possible options
    */

    if (1) {
	iinfo.maxiters = libset_get_int(BHHH_MAXITER);
	iinfo.tol = libset_get_double(BHHH_TOLER);
    } else if (sr->ci == NLS || sr->ci == MLE || sr->ci == GMM) {
	iinfo.maxiters = 400; /* FIXME: 100*(n+1) or 200*(n+1) */
	iinfo.tol = libset_get_double(NLS_TOLER);
    } 
    
    edit_dialog("Iterative estimation",
		"Tolerance for convergence",
		NULL,
		NULL,
		&iinfo,
		ITERATIONS,
		VARCLICK_NONE,
		&cancel);

    if (!cancel) {
	fprintf(stderr, "iters=%d, tol=%g\n", iinfo.maxiters,
		iinfo.tol);
	/* FIXME check and set values */
    }
}

#endif

#ifdef HAVE_X12A

static gboolean x12a_vs_hessian (GtkWidget *w, selector *sr)
{
    if (sr->hess_button != NULL) {
	gboolean s = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

	gtk_widget_set_sensitive(sr->hess_button, !s);
    }

    return FALSE;
}

#endif

static void gui_set_mp_bits (GtkComboBox *cb, gpointer p)
{
    gint i = gtk_combo_box_get_active(cb);
    int j, b = 256;

    for (j=0; j<i; j++) {
	b *= 2;
    }

    set_mp_bits(b);
}

static GtkWidget *mpols_bits_selector (void)
{
    GtkWidget *w, *hbox, *combo;
    char bstr[8];
    int bits = get_mp_bits();
    int b, i, deflt = 0;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Bits per floating-point value"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);

    combo = gtk_combo_box_new_text();

    i = 0;
    for (b=256; b<=4096; b*=2) {
	sprintf(bstr, "%d", b);
	gtk_combo_box_append_text(GTK_COMBO_BOX(combo), bstr);
	if (b == bits) {
	    deflt = i;
	}
	i++;
    }

    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), deflt);
    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(combo)), "changed",
		     G_CALLBACK(gui_set_mp_bits), NULL);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);

    return hbox;
}

#define robust_conf(c) (c != LOGIT && c != PROBIT && \
                        c != OLOGIT && c != OPROBIT && \
                        c != QUANTREG && c != INTREG && \
                        c != MLOGIT)

static void build_selector_switches (selector *sr) 
{
    GtkWidget *hbox, *tmp;

    if (sr->ci == OLS || sr->ci == WLS || sr->ci == INTREG ||
	sr->ci == GARCH || sr->ci == IVREG || sr->ci == VAR || 
	sr->ci == LOGIT || sr->ci == PROBIT || sr->ci == MLOGIT ||
	sr->ci == OLOGIT || sr->ci == OPROBIT ||
	sr->ci == PANEL || sr->ci == QUANTREG) {
	GtkWidget *b1;

	/* FIXME arma robust variant? */

	vbox_add_hsep(sr->vbox);

	if (sr->ci == QUANTREG) {
	    b1 = gtk_check_button_new_with_label(_("Robust standard errors/intervals"));
	} else {
	    b1 = gtk_check_button_new_with_label(_("Robust standard errors"));
	}

	g_object_set_data(G_OBJECT(b1), "opt", GINT_TO_POINTER(OPT_R));
	g_signal_connect(G_OBJECT(b1), "toggled",
			 G_CALLBACK(option_callback), sr);

	if (robust_conf(sr->ci) && using_hc_by_default()) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
	} else if (model_opt & OPT_R) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
	}

	if (sr->ci == PANEL) {
	    g_object_set_data(G_OBJECT(sr->dlg), "robust-button", b1);
	}

	hbox = gtk_hbox_new(FALSE, 5);
	gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 0);
	gtk_widget_show(b1);

	if (robust_conf(sr->ci)) {
	    GtkWidget *b2;

	    b2 = gtk_button_new_with_label(_("configure"));
	    g_signal_connect(G_OBJECT(b2), "clicked",
			     G_CALLBACK(hc_config), sr);
	    gtk_widget_set_sensitive(b2, using_hc_by_default());
	    sensitize_conditional_on(b2, b1);
	    gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 0);
	    gtk_widget_show(b2);
	}

	gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show(hbox);
    }

    if (sr->ci == TOBIT || sr->ci == ARMA || sr->ci == GARCH ||
	sr->ci == LOGIT || sr->ci == PROBIT || sr->ci == HECKIT ||
	sr->ci == OLOGIT || sr->ci == OPROBIT || sr->ci == MLOGIT) {
	if (sr->ci == ARMA) {
	    vbox_add_hsep(sr->vbox);
	    tmp = gtk_check_button_new_with_label(_("Include a constant"));
	    pack_switch(tmp, sr, arma_const, TRUE, OPT_N, 0);
	}
#if 0 /* testing */
	hbox = gtk_hbox_new(FALSE, 5);
	tmp = gtk_button_new_with_label(_("Configure"));
	g_signal_connect(G_OBJECT(tmp), "clicked",
			 G_CALLBACK(call_iters_dialog), sr);
	gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
	gtk_widget_show_all(hbox);
#else
	tmp = gtk_check_button_new_with_label(_("Show details of iterations"));
	pack_switch(tmp, sr, verbose, FALSE, OPT_V, 0);
#endif
    } else if (sr->ci == COINT2 || sr->ci == VECM || 
	       sr->ci == VAR || sr->ci == VLAGSEL) {
	if (sr->ci == VAR || sr->ci == VLAGSEL) {
	    tmp = gtk_check_button_new_with_label(_("Include a constant"));
	    pack_switch(tmp, sr, varconst, TRUE, OPT_N, 0);
	    tmp = gtk_check_button_new_with_label(_("Include a trend"));
	    pack_switch(tmp, sr, vartrend, FALSE, OPT_T, 0);
	} else {
	    tmp = gtk_check_button_new_with_label(_("Show details of regressions"));
	    pack_switch(tmp, sr, verbose, FALSE, OPT_V, 0);
	}
	tmp = gtk_check_button_new_with_label(_("Include seasonal dummies"));
	pack_switch(tmp, sr, 
		    want_seasonals && (datainfo->pd == 4 || datainfo->pd == 12),
		    FALSE, OPT_D, 0);
	if (datainfo->pd != 4 && datainfo->pd != 12) {
	    gtk_widget_set_sensitive(tmp, FALSE);
	}
    } else if (sr->ci == HILU) {
	tmp = gtk_check_button_new_with_label(_("Fine-tune using Cochrane-Orcutt"));
	pack_switch(tmp, sr, TRUE, TRUE, OPT_B, 0);
    } else if (sr->ci == COINT) {
	tmp = gtk_check_button_new_with_label(_("Test down from maximum lag order"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_E, 0);
	tmp = gtk_check_button_new_with_label(_("Skip initial DF tests"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_S, 0);
    } else if (sr->ci == PANEL_WLS) {
	tmp = gtk_check_button_new_with_label(_("Iterated weighted least squares"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_T, 0);
    } else if (sr->ci == PANEL || sr->ci == ARBOND) {
	tmp = gtk_check_button_new_with_label(_("Include time dummies"));
	pack_switch(tmp, sr, (model_opt & OPT_D), FALSE, OPT_D, 0);
    } else if (sr->ci == XTAB) {
	tmp = gtk_check_button_new_with_label(_("Show zeros explicitly"));
	pack_switch(tmp, sr, FALSE, FALSE, OPT_Z, 0);
    } else if (sr->ci == SPEARMAN) {
	tmp = gtk_check_button_new_with_label(_("Show rankings"));
	pack_switch(tmp, sr, verbose, FALSE, OPT_V, 0);
    } else if (sr->ci == CORR) {	
	tmp = gtk_check_button_new_with_label(_("Ensure uniform sample size"));
	pack_switch(tmp, sr, verbose, FALSE, OPT_U, 0);
    } 

    if (sr->ci == ARMA) {
	sr->hess_button = 
	    gtk_check_button_new_with_label(_("Parameter covariance matrix via Hessian"));
	pack_switch(sr->hess_button, sr, arma_hessian, TRUE, OPT_G, 0);
#ifdef HAVE_X12A   
	sr->x12a_button = gtk_check_button_new_with_label(_("Use X-12-ARIMA"));
	pack_switch(sr->x12a_button, sr, arma_x12, FALSE, OPT_X, 0);
	g_signal_connect(G_OBJECT(sr->x12a_button), "toggled",
			 G_CALLBACK(x12a_vs_hessian), sr);
	if (model_opt & OPT_X) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->x12a_button), TRUE);
	}
#endif
    } else if (sr->ci == GARCH) {
	tmp = gtk_check_button_new_with_label(_("Standardize the residuals"));
	pack_switch(tmp, sr, (model_opt & OPT_U), FALSE, OPT_U, 0);
	tmp = gtk_check_button_new_with_label(_("Use Fiorentini et al algorithm"));
	pack_switch(tmp, sr, libset_get_bool(USE_FCP), FALSE, OPT_F, 0);
    } else if (sr->ci == MPOLS) {
	hbox = mpols_bits_selector();
	gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
	gtk_widget_show_all(hbox);
    }
} 

static void unhide_lags_callback (GtkWidget *w, selector *sr)
{
    int i, show_lags;
    GtkListStore *store;
    GtkTreeIter iter;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    show_lags = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    for (i=1; i<datainfo->v; i++) {
	if (list_show_var(i, sr->ci, show_lags)) {
	    list_append_var_simple(store, &iter, i);
	}
    }
}

static void unhide_lags_switch (selector *sr) 
{
    GtkWidget *hbox;
    GtkWidget *b;

    vbox_add_hsep(sr->vbox);

    b = gtk_check_button_new_with_label(_("Show lagged variables"));
    g_signal_connect(G_OBJECT(b), "toggled",
		     G_CALLBACK(unhide_lags_callback), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

#if 0 /* not ready */ 

static void boot_switch_callback (GtkWidget *w, selector *sr)
{
    if (get_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)) {
	sr->opts |= OPT_P;
    } else {
	sr->opts &= ~OPT_P;
    }
}

static void test_boot_switch (selector *sr) 
{
    GtkWidget *hbox;
    GtkWidget *b;

    vbox_add_hsep(sr->vbox);

    b = gtk_check_button_new_with_label(_("Use bootstrap"));
    g_signal_connect(G_OBJECT(b), "toggled",
		     G_CALLBACK(boot_switch_callback), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show_all(hbox);
}

#endif 

static void pack_switch_with_extra (GtkWidget *b, selector *sr,
				    gboolean checked, gretlopt opt, 
				    int child, GtkWidget *extra,
				    const gchar *extra_text)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    gint offset = (child)? 15 : 0;
    gint i = opt;

    g_object_set_data(G_OBJECT(b), "opt", GINT_TO_POINTER(i));

    g_signal_connect(G_OBJECT(b), "toggled", 
		     G_CALLBACK(option_callback), sr);
    if (checked) {
	sr->opts |= opt;
    } else {
	sr->opts &= ~opt;
    }

    gtk_box_pack_start(GTK_BOX(hbox), b, TRUE, TRUE, offset);
    gtk_widget_show(b);

    if (extra_text != NULL) {
	GtkWidget *w = gtk_label_new(extra_text);

	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	gtk_widget_show(w);
    }

    gtk_box_pack_start(GTK_BOX(hbox), extra, TRUE, TRUE, offset);
    gtk_widget_show(extra);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), checked);
}

static GtkWidget *alpha_spinner (double deflt, double minval)
{
    GtkObject *adj;
    
    adj = gtk_adjustment_new(deflt, minval, 0.99, 0.01, 0.1, 0);
    return gtk_spin_button_new(GTK_ADJUSTMENT(adj), 1, 2);
}

static void build_quantreg_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Compute standard errors"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, 
					 _("Compute confidence intervals"));
    sr->extra[1] = alpha_spinner(0.90, 0.70);
    pack_switch_with_extra(b2, sr, FALSE, OPT_I, 0, sr->extra[1], "1 - α =");
    gtk_widget_set_sensitive(sr->extra[1], FALSE);
    sensitize_conditional_on(sr->extra[1], b2);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_ellipse_spinner (selector *sr)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *label;

    vbox_add_hsep(sr->vbox);

    label = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);
    gtk_widget_show(label);
    label = gtk_label_new(_("Confidence level"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);
    label = gtk_label_new("1 - α =");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_widget_show(label);
    sr->extra[0] = alpha_spinner(0.95, 0.70);
    gtk_box_pack_start(GTK_BOX(hbox), sr->extra[0], FALSE, FALSE, 0);
    gtk_widget_show(sr->extra[0]);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
}

static void build_rankcorr_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Spearman's rho"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Kendall's tau"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_K, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_pca_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Use correlation matrix"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Use covariance matrix"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_C, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_pvalues_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Show slopes at mean"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Show p-values"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_P, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_scatters_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Use points"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 1);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Use lines"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_L, 1);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void auto_omit_callback (GtkWidget *w, selector *sr)
{
    gboolean s = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    if (sr->extra[0] != NULL) {
	gtk_widget_set_sensitive(sr->extra[0], s);
    }
    if (sr->lvars != NULL) {
	gtk_widget_set_sensitive(sr->lvars, !s);
    }
    if (sr->rvars1 != NULL) {
	gtk_widget_set_sensitive(sr->rvars1, !s);
    }

    gtk_widget_set_sensitive(sr->add_button, !s);
    gtk_widget_set_sensitive(sr->remove_button, !s);
}

static void build_omit_test_radios (selector *sr)
{
    windata_t *vwin = (windata_t *) sr->data;
    MODEL *pmod = (MODEL *) vwin->data;
    GtkWidget *b1, *b2, *b3;
    GSList *group;

    vbox_add_hsep(sr->vbox);

    b1 = gtk_radio_button_new_with_label(NULL, _("Estimate reduced model"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);
    sr->radios[0] = b1;

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Wald test, based on covariance matrix"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_W, 0);
    sr->radios[1] = b2;

    if (pmod->ci != PANEL) {
	group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
	b3 = gtk_radio_button_new_with_label(group, _("Sequential elimination of variables\n"
						      "using two-sided p-value:"));
	g_signal_connect(G_OBJECT(b3), "toggled",
			 G_CALLBACK(auto_omit_callback), sr);

	sr->extra[0] = alpha_spinner(0.10, 0.01);
	pack_switch_with_extra(b3, sr, FALSE, OPT_A, 0, sr->extra[0], NULL);
	gtk_widget_set_sensitive(sr->extra[0], FALSE);

	sr->radios[2] = b3;
    }
}

static void build_panel_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;
    gboolean fe = TRUE;

    if (sr->opts & OPT_W) {
	/* panel weighted least squares */
	return;
    }

    if (model_opt & OPT_U) {
	fe = FALSE;
    }

    b1 = gtk_radio_button_new_with_label(NULL, _("Fixed effects"));
    pack_switch(b1, sr, fe, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Random effects"));
    pack_switch(b2, sr, !fe, FALSE, OPT_U, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_arbond_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("One-step estimation"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Two-step estimation"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_T, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_heckit_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;
    gboolean ml = TRUE;

    if (model_opt & OPT_T) {
	ml = FALSE;
    }

    b1 = gtk_radio_button_new_with_label(NULL, _("Maximum likelihood estimation"));
    pack_switch(b1, sr, ml, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("2-step estimation"));
    pack_switch(b2, sr, !ml, FALSE, OPT_T, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
}

static void build_xtab_radios (selector *sr)
{
    GtkWidget *b1, *b2, *b3;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Plain numerical values"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Show row percentages"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_R, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    b3 = gtk_radio_button_new_with_label(group, _("Show column percentages"));
    pack_switch(b3, sr, FALSE, FALSE, OPT_C, 0);

    sr->radios[0] = b1;
    sr->radios[1] = b2;
    sr->radios[2] = b3;
}

static gboolean arma_estimator_switch (GtkComboBox *box, selector *sr)
{
    if (sr->hess_button != NULL) {
	GtkWidget *xb = sr->x12a_button;

	if (xb == NULL || !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(xb))) {
	    gchar *s = gtk_combo_box_get_active_text(box);

	    gtk_widget_set_sensitive(sr->hess_button, 
				     !strcmp(s, _("Exact Maximum Likelihood")));
	    g_free(s);
	}
    }

    return FALSE;
}

static void build_arma_combo (selector *sr)
{
    GtkWidget *hbox, *combo;
    static const char *opt_strs[] = {
        N_("Exact Maximum Likelihood"),
        N_("Conditional Maximum Likelihood"),
	NULL
    };
    static gretlopt opts[] = { 
	OPT_NONE, 
	OPT_C, 
    };
    static combo_opts arma_opts;
    int deflt = 0;

    arma_opts.strs = opt_strs;
    arma_opts.vals = opts;
    arma_opts.optp = &sr->opts;

    if (model_opt & OPT_C) {
	deflt = 1;
    }

    combo = gretl_opts_combo_full(&arma_opts, deflt, 
				  G_CALLBACK(arma_estimator_switch),
				  sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    gtk_widget_show(combo);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);
}

static void build_coint_combo (selector *sr)
{
    GtkWidget *hbox, *combo;
    static const char *opt_strs[] = {
        N_("test without constant"),
        N_("test with constant"),
        N_("with constant and trend"),
        N_("with constant and quadratic trend"),
	NULL
    };
    static gretlopt opts[] = { 
	OPT_N, 
	OPT_NONE, 
	OPT_T, 
	OPT_R, 
    };
    static combo_opts coint_opts;
    int deflt = 1;

    coint_opts.strs = opt_strs;
    coint_opts.vals = opts;
    coint_opts.optp = &sr->opts;

    combo = gretl_opts_combo(&coint_opts, deflt);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    gtk_widget_show(combo);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);    
}

static void build_vecm_combo (selector *sr)
{
    GtkWidget *hbox, *combo;
    static const char *opt_strs[] = {
	N_("No constant"),
	N_("Restricted constant"),
	N_("Unrestricted constant"),
	N_("Restricted trend"),
	N_("Unrestricted trend"),
	NULL
    };
    static gretlopt opts[] = {
	OPT_N, 
	OPT_R, 
	OPT_NONE, 
	OPT_A, 
	OPT_T
    };    
    static combo_opts vecm_opts;

    vecm_opts.strs = opt_strs;
    vecm_opts.vals = opts;
    vecm_opts.optp = &sr->opts;

    combo = gretl_opts_combo(&vecm_opts, jcase);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    gtk_widget_show(combo);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);
}

static void build_selector_radios (selector *sr)
{
    if (sr->ci == PANEL) {
	build_panel_radios(sr);
    } else if (sr->ci == ARBOND) {
	build_arbond_radios(sr);
    } else if (sr->ci == SCATTERS) {
	build_scatters_radios(sr);
    } else if (sr->ci == OMIT) {
	build_omit_test_radios(sr);
    } else if (sr->ci == LOGIT || sr->ci == PROBIT) {
	build_pvalues_radios(sr);
    } else if (sr->ci == HECKIT) {
	build_heckit_radios(sr);
    } else if (sr->ci == XTAB) {
	build_xtab_radios(sr);
    } else if (sr->ci == SPEARMAN) {
	build_rankcorr_radios(sr);
    } else if (sr->ci == PCA) {
	build_pca_radios(sr);
    } else if (sr->ci == QUANTREG) {
	build_quantreg_radios(sr);
    }
}

static void build_selector_combo (selector *sr)
{
    if (sr->ci == ARMA) {
	build_arma_combo(sr);
    } else if (sr->ci == COINT) {
	build_coint_combo(sr);
    } else if (sr->ci == VECM || sr->ci == COINT2) {
	build_vecm_combo(sr);
    } else if (sr->ci == IV_GMM) {
	build_gmm_popdown(sr);
    }
}

static void lag_selector_button (selector *sr)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);

    sr->lags_button = gtk_button_new_with_label(_("lags..."));

    g_signal_connect(G_OBJECT(sr->lags_button), "clicked", 
		     G_CALLBACK(lags_dialog_driver), sr);
    if (varlist_row_count(sr, SR_RVARS1, NULL) < 2) {
	gtk_widget_set_sensitive(sr->lags_button, FALSE);
    }

    gtk_box_pack_start(GTK_BOX(hbox), sr->lags_button, FALSE, FALSE, 0);
    gtk_widget_show(sr->lags_button);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_show(hbox);
}    

static void selector_doit (GtkWidget *w, selector *sr)
{
    compose_cmdlist(sr);

    if (sr->error == 0) {
	int err = sr->callback(sr);

	if (!err && open_selector != NULL) {
	    gtk_widget_destroy(sr->dlg);
	}
    } 
}

static void build_selector_buttons (selector *sr)
{
    GtkWidget *tmp;

    if (sr->ci != PRINT && sr->ci != SAVE_FUNCTIONS &&
	sr->ci != DEFINE_LIST && sr->ci != DEFINE_MATRIX &&
	sr->ci != ELLIPSE && !SAVE_DATA_ACTION(sr->ci)) {
	int ci = sr->ci;

	if (sr->ci == OLOGIT || sr->ci == MLOGIT) {
	    ci = LOGIT;
	} else if (sr->ci == OPROBIT) {
	    ci = PROBIT;
	} else if (IV_MODEL(sr->ci)) {
	    ci = IVREG;
	}

	tmp = gtk_button_new_from_stock(GTK_STOCK_HELP);
	GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
	gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
	gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(sr->action_area),
					   tmp, TRUE);
	g_signal_connect(G_OBJECT(tmp), "clicked", 
			 G_CALLBACK(context_help), 
			 GINT_TO_POINTER(ci));
	gtk_widget_show(tmp);
    }

    tmp = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(clear_vars), sr);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(cancel_selector), sr);
    gtk_widget_show(tmp);

    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(selector_doit), sr);
    gtk_widget_show(tmp);
    gtk_widget_grab_default(tmp);
}

static int list_show_var (int v, int ci, int show_lags)
{
    int ret = 1;

    lags_hidden = 0;

    if (v == 0 && (ci == DEFINE_LIST || ci == DEFINE_MATRIX)) {
	;
    } else if (v == 0 && (!MODEL_CODE(ci) || ci == ARMA)) {
	ret = 0;
    } else if (var_is_hidden(datainfo, v)) {
	ret = 0;
    } else if (!show_lags && is_standard_lag(v, datainfo, NULL)) {
	lags_hidden = 1;
	ret = 0;
    } else if (ci == XTAB) {
	ret = 0;
	if (var_is_discrete(datainfo, v) || 
	    gretl_isdiscrete(datainfo->t1, datainfo->t2, Z[v])) {
	    ret = 1;
	}
    }

    return ret;
}

static GtkWidget *selection_dialog_top_label (int ci)
{
    gchar *s;

    if (MODEL_CODE(ci) || VEC_CODE(ci))
	s = _(est_str(ci));
    else if (ci == GR_XY)
	s = _("XY scatterplot");
    else if (ci == GR_IMP)
	s = _("plot with impulses");
    else if (ci == GR_3D)
	s = _("3D plot");
    else if (ci == SCATTERS)
	s = _("multiple scatterplots");
    else if (ci == GR_DUMMY)
	s = _("factorized plot");
    else if (ci == GR_XYZ)
	s = _("scatterplot with control");
    else if (ci == SAVE_FUNCTIONS) 
	s = _("functions to package");
    else if (ci == ANOVA) 
	s = _("ANOVA");
    else
	s = "fixme need string";

    return gtk_label_new(s);
}

static int has_0 (const int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    return 1;
	}
    }

    return 0;
}

static void primary_rhs_varlist (selector *sr, GtkWidget *vbox)
{
    GtkTreeModel *mod;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *remove;
    GtkWidget *hbox;
    GtkWidget *bvbox;
    GtkWidget *tmp = NULL;
    int i;

    if (COINT_CODE(sr->ci)) {
	tmp = gtk_label_new(_("Variables to test"));
    } else if (VEC_CODE(sr->ci)) {
	tmp = gtk_label_new(_("Endogenous variables"));
    } else if (MODEL_CODE(sr->ci)) {
	tmp = gtk_label_new(_("Independent variables"));
    } else if (sr->ci == GR_XY || sr->ci == GR_IMP) {
	tmp = gtk_label_new(_("Y-axis variables"));
    } else if (sr->ci == SCATTERS) {
	multiplot_label = tmp = gtk_label_new(_("X-axis variables"));
    } else if (sr->ci == SAVE_FUNCTIONS) {
	tmp = gtk_label_new(_("Helper functions"));
    }

    if (tmp != NULL) {
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    }

    hbox = gtk_hbox_new(FALSE, 5);

    /* push/pull buttons first, in their own little vbox */
    bvbox = gtk_vbox_new(TRUE, 0);

    tmp = add_button(sr, bvbox);
    g_signal_connect (G_OBJECT(tmp), "clicked", 
		      G_CALLBACK(add_to_rvars1_callback), sr);
    
    remove = remove_button(sr, bvbox);

    if (sr->ci == ARMA || sr->ci == VAR) { /* FIXME VECM */
	tmp = gtk_button_new_with_label(_("lags..."));
	gtk_box_pack_start(GTK_BOX(bvbox), tmp, TRUE, FALSE, 0);
	g_signal_connect(G_OBJECT(tmp), "clicked", 
			 G_CALLBACK(lags_dialog_driver), sr);
	gtk_widget_set_sensitive(tmp, FALSE);
	if (sr->ci == ARMA) {
	    sr->lags_button = tmp;
	} else {
	    sr->extra[EXTRA_LAGS] = tmp;
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), bvbox, TRUE, TRUE, 0);
    gtk_widget_show_all(bvbox);

    /* then the actual primary RHS listbox */
    sr->rvars1 = var_list_box_new(GTK_BOX(hbox), sr, SR_RVARS1);
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));

    store = GTK_LIST_STORE(mod);
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(mod, &iter);

    if (MODEL_CODE(sr->ci)) {
	if (sr->ci == ARMA) {
	    g_object_set_data(G_OBJECT(sr->rvars1), "selector", sr);
	} else if (xlist == NULL || has_0(xlist)) {
	    /* stick the constant in by default */
	    list_append_var(mod, &iter, 0, sr, SR_RVARS1);
	} 
	if (xlist != NULL) {
	    int nx = 0;

	    /* we have a saved list of regressors */
	    for (i=1; i<=xlist[0]; i++) {
		if (xlist[i] != 0) {
		    list_append_var(mod, &iter, xlist[i], sr, SR_RVARS1);
		    nx++;
		}
	    }
	    if (nx > 0 && sr->ci == ARMA) {
		gtk_widget_set_sensitive(sr->lags_button, TRUE);
	    }
	}
    } else if (VEC_CODE(sr->ci) && veclist != NULL && veclist[0] > 0) {
	for (i=1; i<=veclist[0]; i++) {
	    list_append_var(mod, &iter, veclist[i], sr, SR_RVARS1);
	}
	if (sr->extra[EXTRA_LAGS] != NULL) {
	    gtk_widget_set_sensitive(sr->extra[EXTRA_LAGS], TRUE);
	}
    }

    /* hook remove button to listing */
    g_signal_connect(G_OBJECT(remove), "clicked", 
		     G_CALLBACK(remove_from_right_callback), 
		     sr->rvars1);

    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
    gtk_widget_show(hbox);
}

void selection_dialog (const char *title, int (*callback)(), guint ci)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *right_vbox, *tmp;
    GtkWidget *big_hbox;
    selector *sr;
    int preselect;
    int yvar = 0;
    int i;

    preselect = presel;
    presel = 0;

    if (open_selector != NULL) {
	gdk_window_raise(open_selector->dlg->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) return;

    selector_init(sr, ci, title, NULL, callback);

    tmp = selection_dialog_top_label(ci);
    gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    /* the following encloses LHS lvars, depvar and indepvar stuff */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* LHS: list of elements to choose from */
    sr->lvars = var_list_box_new(GTK_BOX(big_hbox), sr, SR_LVARS);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (ci == SAVE_FUNCTIONS) {
	functions_list(sr);
    } else {
	for (i=0; i<datainfo->v; i++) {
	    if (list_show_var(i, ci, 0)) {
		list_append_var_simple(store, &iter, i);
	    }
	}
    }

    /* RHS: vertical holder */
    right_vbox = gtk_vbox_new(FALSE, 5);

    if (MODEL_CODE(ci) || ci == ANOVA) { 
	/* models: top right -> dependent variable */
	yvar = build_depvar_section(sr, right_vbox, preselect);
    } else if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY ||
	       ci == SCATTERS || ci == GR_3D || ci == GR_XYZ) {
	/* graphs: top right -> x-axis variable */
	build_x_axis_section(sr, right_vbox);
    } else if (ci == SAVE_FUNCTIONS) {
	build_public_iface_section(sr, right_vbox);
    }

    /* middle right: used for some estimators and factored plot */
    if (ci == WLS || ci == AR || ci == ARCH || USE_ZLIST(ci) ||
	VEC_CODE(ci) || ci == POISSON || ci == ARBOND ||
	ci == QUANTREG || ci == INTREG || THREE_VARS_CODE(ci)) {
	build_mid_section(sr, right_vbox);
    }
    
    if (THREE_VARS_CODE(ci)) {
	/* choose extra var for plot */
	third_var_box(sr, right_vbox);
    } else if (AUX_LAST(ci)) {
	auxiliary_rhs_varlist(sr, right_vbox);
    } else {
	/* all other uses: list of vars */
	primary_rhs_varlist(sr, right_vbox);
    }

    /* pack the whole RHS to the right of the LHS lvars */
    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, TRUE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    if (ci == ARMA) {
	/* AR, D, MA for ARIMA */
	build_arma_spinners(sr);
    } else if (ci == GARCH) { 
	/* P and Q for GARCH */
	build_garch_spinners(sr);
    }

    /* toggle switches for some cases */
    if (WANT_TOGGLES(ci)) {
	build_selector_switches(sr);
    }

    /* radio buttons for some */
    if (want_radios(sr)) {
	build_selector_radios(sr);
    }

    /* drop-down selector for some */
    if (want_combo(sr)) {
	build_selector_combo(sr);
    }    

    /* plus lag selection stuff, if relevant */
    if (dataset_lags_ok(datainfo)) {
	if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY || \
	    ci == SCATTERS || ci == GR_3D || ci == GR_XYZ) {
	    unhide_lags_switch(sr);
	}
	if (MODEL_CODE(ci) && ci != ARMA) {
	    lag_selector_button(sr);
	} 
	if (select_lags_depvar(ci) && yvar > 0) {
	    maybe_activate_depvar_lags(sr->depvar, sr);
	    maybe_insert_depvar_lags(sr, yvar, 0);
	}
    } 

    /* buttons: Help, Clear, Cancel, OK */
    build_selector_buttons(sr);

    gtk_widget_show(sr->dlg);
}

static char *get_topstr (int cmdnum)
{
    switch (cmdnum) {    
    case LOGS:
	return N_("Select variables for logging");
    case LAGS:
	return N_("Select variables for lagging");
    case SQUARE:
	return N_("Select variables to square");
    case DIFF:
	return N_("Select variables to difference");
    case LDIFF:
	return N_("Select variables to log-difference");
    case ADD:
	return N_("Select variables to add");
    case OMIT:
    case VAROMIT:
	return N_("Select variables to omit");
    case COEFFSUM:
	return N_("Select coefficients to sum");
    case SPEARMAN:
	return N_("Select two variables");
    case ELLIPSE:
	return N_("Confidence region: select two variables");
    case PRINT:
	return N_("Select variables to display");
    case GR_PLOT: 
    case GR_BOX: 
    case GR_NBOX:
	return N_("Select variables to plot");
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case EXPORT_CSV:
    case EXPORT_R:
    case EXPORT_OCTAVE:
    case EXPORT_JM:
    case EXPORT_DAT:
	return N_("Select variables to save");
    case COPY_CSV:
	return N_("Select variables to copy");
    case DEFINE_LIST:
	return N_("Define named list");
    default:
	return NULL;
    }
}

static int add_omit_list (gpointer p, selector *sr)
{
    windata_t *vwin = (windata_t *) p;
    GRETL_VAR *var = NULL;
    MODEL *pmod = NULL;
    GtkListStore *store;
    GtkTreeIter iter;
    int i, nvars = 0;

    if (sr->ci == VAROMIT) {
	var = (GRETL_VAR *) vwin->data;
    } else {
	pmod = (MODEL *) vwin->data;
    }

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (sr->ci == ELLIPSE) {
	char pname[VNAMELEN];
	int nc = gretl_model_get_int(pmod, "base-coeffs"); /* FIXME? */

	if (nc == 0) {
	    nc = pmod->ncoeff;
	}

	for (i=0; i<nc; i++) {
	    gretl_model_get_param_name(pmod, datainfo, i, pname);
	    gtk_list_store_append(store, &iter);
	    gtk_list_store_set(store, &iter, 
			       0, i, 1, 0,
			       2, pname,
			       -1);
	    nvars++;
	}
	g_object_set_data(G_OBJECT(sr->lvars), "keep_names", 
			  GINT_TO_POINTER(1));
    } else if (sr->ci == OMIT || sr->ci == ADD || sr->ci == COEFFSUM) {
	int *xlist = gretl_model_get_x_list(pmod);

	if (xlist == NULL) {
	    return 0;
	}

	if (sr->ci == ADD) {
	    int dv = gretl_model_get_depvar(pmod);

	    for (i=0; i<datainfo->v; i++) {
		if (!in_gretl_list(xlist, i) && i != dv &&
		    !var_is_hidden(datainfo, i)) {
		    gtk_list_store_append(store, &iter);
		    gtk_list_store_set(store, &iter, 
				       0, i, 1, 0, 
				       2, datainfo->varname[i],
				       -1);
		    nvars++;
		}
	    }
	} else {	    
	    for (i=1; i<=xlist[0]; i++) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, xlist[i], 1, 0,
				   2, datainfo->varname[xlist[i]],
				   -1);
		nvars++;
	    }
	} 
	free(xlist);
    } else if (sr->ci == VAROMIT) {
	const int *xlist;

	xlist = gretl_VAR_get_exo_list(var);
	if (xlist != NULL) {
	    for (i=1; i<=xlist[0]; i++) {
		gtk_list_store_append(store, &iter);
		gtk_list_store_set(store, &iter, 
				   0, xlist[i], 1, 0,
				   2, datainfo->varname[xlist[i]],
				   -1);
		nvars++;
	    }	    
	}
    } 

    return nvars;
}

static int functions_list (selector *sr)
{
    GtkListStore *store;
    GtkTreeIter iter;
    const char *fnname;
    int n = 0;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    function_names_init();
    
    while ((fnname = next_free_function_name()) != NULL) {
	gtk_list_store_append(store, &iter);
	gtk_list_store_set(store, &iter,
			   0, n++, 1, 0,
			   2, fnname, -1);
    }

    g_object_set_data(G_OBJECT(sr->lvars), "keep_names",
		      GINT_TO_POINTER(1));

    return n;
}

static GtkWidget *simple_selection_top_label (int ci)
{
    GtkWidget *label = NULL;
    const char *str = get_topstr(ci);

    if (str != NULL && *str != '\0') {
	label = gtk_label_new(_(str));
    } 

    return label;
}

static gboolean remove_busy_signal (GtkWidget *w, windata_t *vwin)
{
    if (vwin != NULL) {
	unset_window_busy(vwin);
    }
    return FALSE;
}

static void 
add_to_rvars1_from_named_list (selector *sr, const int *list)
{
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int i, v;

    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (mod == NULL) {
	return;
    }

    gtk_tree_model_get_iter_first(mod, &iter);
    for (i=1; i<=list[0]; i++) {
	v = list[i];
	gtk_list_store_append(GTK_LIST_STORE(mod), &iter);
	gtk_list_store_set(GTK_LIST_STORE(mod), &iter, 
			   0, v, 1, 0, 2, datainfo->varname[v], -1);
    }
}

static void maybe_set_listdef_vars (selector *sr)
{
    const char *lname = selector_entry_text(sr);

    if (lname != NULL && *lname != 0) {
	int *list = get_list_by_name(lname);

	if (list != NULL) {
	    add_to_rvars1_from_named_list(sr, list);
	}
    }
}

static void maybe_set_tsplot_vars (selector *sr)
{
    int mc = mdata_selection_count();

    if (mc > 1 && mc < 7) {
	set_vars_from_main(sr);
    }
}

static void selector_add_top_entry (selector *sr)
{
    GtkWidget *src;
    GtkWidget *hbox;
    GtkWidget *label;
    GtkWidget *entry;
    const char *lname;

    src = GTK_WIDGET(sr->data);
    lname = gtk_entry_get_text(GTK_ENTRY(src));

    hbox = gtk_hbox_new(FALSE, 0);
    label = gtk_label_new(_("Name for list:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), 31);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
    if (lname != NULL && *lname != 0 && strcmp(lname, "null")) {
	gtk_entry_set_text(GTK_ENTRY(entry), lname);
	gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    }
    gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);

    sr->extra[0] = entry;

    gtk_widget_show_all(hbox);
}

#if 0 /* not ready */
static int ols_omit_select (windata_t *vwin)
{
    if (vwin->role == VIEW_MODEL) {
	MODEL *pmod = vwin->data;

	return bootstrap_ok(pmod->ci);
    } else {
	return 0;
    }
}
#endif

static void maybe_prefill_RHS (selector *sr)
{
    int *list = main_window_selection_as_list();

    if (list != NULL && list[0] >= 2) {
	GtkListStore *store;
	GtkTreeIter iter;
	int i;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1)));
	gtk_list_store_clear(store);
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

	for (i=1; i<=list[0]; i++) {
	    list_append_var_simple(store, &iter, list[i]);
	}
    }

    free(list);
}

void simple_selection (const char *title, int (*callback)(), guint ci,
		       gpointer p) 
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *left_vbox, *button_vbox, *right_vbox, *tmp;
    GtkWidget *top_hbox, *big_hbox;
    selector *sr;
    int nleft = 0;
    int i;

    if (open_selector != NULL) {
	gdk_window_raise(open_selector->dlg->window);
	return;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) {
	return;
    }

    selector_init(sr, ci, title, p, callback);

    tmp = simple_selection_top_label(ci);
    if (tmp != NULL) {
	gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
	gtk_widget_show(tmp);
    } 

    /* entry field for some uses */
    if (ci == DEFINE_LIST) {
	selector_add_top_entry(sr);
    }

    /* for titles */
    top_hbox = gtk_hbox_new(FALSE, 0); 
    gtk_box_set_homogeneous(GTK_BOX(top_hbox), TRUE);

    tmp = gtk_label_new(_("Available vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(" ");
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    tmp = gtk_label_new(_("Selected vars"));
    gtk_box_pack_start(GTK_BOX(top_hbox), tmp, FALSE, FALSE, 5);
    gtk_widget_show(tmp);

    gtk_box_pack_start(GTK_BOX(sr->vbox), top_hbox, FALSE, FALSE, 5);
    gtk_widget_show(top_hbox);

    /* the following will enclose 3 or more vboxes */
    big_hbox = gtk_hbox_new(FALSE, 5); 

    /* holds list of elements available for selection */
    left_vbox = gtk_vbox_new(FALSE, 5);

    sr->lvars = var_list_box_new(GTK_BOX(left_vbox), sr, SR_LVARS);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (ci == OMIT || ci == ADD || ci == COEFFSUM ||
	ci == ELLIPSE || ci == VAROMIT) {
        nleft = add_omit_list(p, sr);
	g_signal_connect(G_OBJECT(sr->dlg), "destroy", 
			 G_CALLBACK(remove_busy_signal), 
			 p);
    } else {
	int start = (ci == DEFINE_LIST || ci == DEFINE_MATRIX)? 0 : 1;

	for (i=start; i<datainfo->v; i++) {
	    if (list_show_var(i, ci, 0)) {
		list_append_var_simple(store, &iter, i);
		nleft++;
	    }
	}
    }

    sr->n_left = nleft;

    gtk_box_pack_start(GTK_BOX(big_hbox), left_vbox, TRUE, TRUE, 0);
    gtk_widget_show(left_vbox);
    
    /* middle: vertical holder for push/pull buttons */
    button_vbox = gtk_vbox_new(TRUE, 0);

    sr->add_button = add_button(sr, button_vbox);
    g_signal_connect(G_OBJECT(sr->add_button), "clicked", 
		     G_CALLBACK(add_to_rvars1_callback), sr);

    if (SAVE_DATA_ACTION(sr->ci)) {
	tmp = gtk_button_new_with_label (_("All ->"));
	gtk_box_pack_start(GTK_BOX(button_vbox), tmp, TRUE, FALSE, 0);
	g_signal_connect (G_OBJECT(tmp), "clicked", 
			  G_CALLBACK(add_all_to_rvars1_callback), sr);
    }
    
    sr->remove_button = remove_button(sr, button_vbox);

    gtk_box_pack_start(GTK_BOX(big_hbox), button_vbox, TRUE, TRUE, 0);
    gtk_widget_show_all(button_vbox);

    /* RHS: vertical holder for selected vars */
    right_vbox = gtk_vbox_new(FALSE, 5);

    sr->rvars1 = var_list_box_new(GTK_BOX(right_vbox), sr, SR_RVARS1);
    g_object_set_data(G_OBJECT(sr->rvars1), "selector", sr);

    gtk_box_pack_start(GTK_BOX(big_hbox), right_vbox, TRUE, TRUE, 0);
    gtk_widget_show(right_vbox);

    /* pre-fill RHS box? Only if we have 2 or more vars selected in the
       main window and if the command is "suitable"
    */
    if (RHS_PREFILL(ci)) {
	maybe_prefill_RHS(sr);
    }

    /* connect removal from right signal */
    g_signal_connect(G_OBJECT(sr->remove_button), "clicked", 
		     G_CALLBACK(remove_from_right_callback), 
		     sr->rvars1);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), big_hbox, TRUE, TRUE, 0);
    gtk_widget_show(big_hbox);

    /* unhide lags check box? */
    if (SAVE_DATA_ACTION(sr->ci) && lags_hidden) {
	unhide_lags_switch(sr);
    }

    /* radio buttons? */
    if (want_radios(sr)) {
	build_selector_radios(sr);
    }

    /* toggle switches for some cases */
    if (WANT_TOGGLES(ci)) {
	build_selector_switches(sr);
    }

    if (ci == ELLIPSE) {
	build_ellipse_spinner(sr);
    }

#if 0 /* not ready */
    /* bootstrap check box? */
    if (ci == OMIT && ols_omit_select(p)) {
	test_boot_switch(sr);
    }
#endif

    /* buttons: Help, Clear, Cancel, OK */
    build_selector_buttons(sr);

    if (TWO_VARS_CODE(sr->ci) && sr->ci != ELLIPSE &&
	mdata_selection_count() == 2) {
	set_vars_from_main(sr);
    } else if (nleft == 1) {
	select_singleton(sr);
    } else if (sr->ci == DEFINE_LIST) {
	maybe_set_listdef_vars(sr);
    } else if (sr->ci == TSPLOTS) {
	maybe_set_tsplot_vars(sr);
    }

    if (nleft == 0) {
	gtk_widget_destroy(sr->dlg);
	warnbox(_("No suitable data are available"));
    } else if ((ci == COEFFSUM || ci == ELLIPSE) && nleft < 2) {
	gtk_widget_destroy(sr->dlg);
	warnbox(_("No suitable data are available"));
    } else {
	gtk_widget_show(sr->dlg);
    }

    if (SAVE_DATA_ACTION(sr->ci)) {
	gretl_set_window_modal(sr->dlg);
    } else if (sr->ci == DEFINE_MATRIX) {
	gtk_main();
    }
}

struct list_maker {
    char *liststr;
    int n_items;
    size_t len;
    int overflow;
};

static void selection_add_item (GtkTreeModel *model, GtkTreePath *path,
				GtkTreeIter *iter, struct list_maker *lmkr)
{
    gchar *varnum = NULL;

    if (lmkr->len > MAXLEN - 12) {
	lmkr->overflow = 1;
	return;
    }

    gtk_tree_model_get (model, iter, 0, &varnum, -1);
    strcat(lmkr->liststr, " ");
    strcat(lmkr->liststr, varnum);
    lmkr->len += strlen(varnum) + 1;
    g_free(varnum);
    lmkr->n_items += 1;
}

char *main_window_selection_as_string (void) 
{
    GtkTreeSelection *select;
    struct list_maker lmkr;

    lmkr.liststr = mymalloc(MAXLEN);
    if (lmkr.liststr == NULL) {
	return NULL;
    }

    lmkr.liststr[0] = 0;
    lmkr.n_items = lmkr.overflow = 0;
    lmkr.len = 0;

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(mdata->listbox));
    gtk_tree_selection_selected_foreach(select, 
					(GtkTreeSelectionForeachFunc) 
					selection_add_item,
					&lmkr); 

    if (lmkr.overflow) {
	errbox(_("Too many items were selected"));
	lmkr.liststr[0] = 0;
    }

    return lmkr.liststr;
}

static const char *data_save_title (int ci)
{
    switch (ci) {
    case EXPORT_CSV:
	return _("Save CSV data file");
    case EXPORT_R:
	return _("Save R data file");
    case EXPORT_OCTAVE:
	return _("Save octave data file");
    default:
	return _("Save data file");
    }
    return "";
}

static int data_save_selection_callback (selector *sr)
{
    gpointer data = sr->data;
    int ci = sr->ci;

    if ((sr->cmdlist == NULL || *sr->cmdlist == 0) && sr->n_left == 0) {
	/* nothing selected */
	return 0;
    }

    if (storelist != NULL) {
	free(storelist);
	storelist = NULL;
    }

    if (sr->cmdlist != NULL && *sr->cmdlist != '\0') {
	storelist = g_strdup(sr->cmdlist);
    }

    gtk_widget_destroy(sr->dlg);

    if (ci == SAVE_FUNCTIONS) {
	prepare_functions_save();
    } else if (ci != COPY_CSV) {
	file_selector(data_save_title(ci), ci, FSEL_DATA_MISC, data);
    }

    return 0;
}

void data_save_selection_wrapper (int file_ci, gpointer p)
{
    if (file_ci == SAVE_FUNCTIONS) {
	if (no_user_functions_check()) {
	    return;
	}
	selection_dialog(_("Save functions"), 
			 data_save_selection_callback, file_ci);
    } else {
	simple_selection((file_ci == COPY_CSV)? 
			 _("Copy data") : _("Save data"), 
			 data_save_selection_callback, file_ci, p);
    }

    gtk_main(); /* the corresponding gtk_main_quit() is in
		   the function destroy_selector() */
}

/* accessor functions */

int selector_code (const selector *sr)
{
    return (sr->ci == PANEL_WLS || sr->ci == PANEL_B)? 
	PANEL : sr->ci;
}

const char *selector_list (const selector *sr)
{
    const char *ret = NULL;

    if (sr->cmdlist != NULL && *sr->cmdlist != '\0') {
	ret = sr->cmdlist;
    }

    return ret;
}

int selector_list_hasconst (const selector *sr)
{
    int hc = sr->cmdlist != NULL && 
	strstr(sr->cmdlist, " 0") != NULL;

    return hc;
}

gpointer selector_get_data (const selector *sr)
{
    return sr->data;
}

gretlopt selector_get_opts (const selector *sr)
{
    gretlopt ret;

    if (sr->ci == PANEL_B) {
	ret = sr->opts | OPT_B;
    } else {
	ret = sr->opts;
    }

    if (sr->ci == VECM || sr->ci == COINT2) {
	/* record Johansen case */
	if (ret & OPT_A) {
	    /* --crt */
	    jcase = J_REST_TREND;
	} else if (ret & OPT_N) {
	    /* --nc */
	    jcase = J_NO_CONST;
	} else if (ret & OPT_R) {
	    /* --rc */
	    jcase = J_REST_CONST;
	} else if (ret & OPT_T) {
	    /* --ct */
	    jcase = J_UNREST_TREND;
	} else {
	    /* unrestricted constant */
	    jcase = J_UNREST_CONST;
	}
    }

    return ret;
}

const char *selector_entry_text (const selector *sr)
{
    g_return_val_if_fail(GTK_IS_ENTRY(sr->extra[0]), NULL);

    return gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));
}

int selector_error (const selector *sr)
{
    return sr->error;
}

void maybe_clear_selector (const int *dlist)
{
    int i, j;

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    for (j=1; j<=dlist[0]; j++) {
		if (xlist[i] >= dlist[j]) {
		    clear_selector();
		    return;
		}
	    }
	}
    }
}

/* ------------- lag selection apparatus -------------- */

#define NOT_LAG 66666

typedef struct var_lag_info_ var_lag_info;

struct var_lag_info_ {
    int v;              /* variable ID number or VDEFLT for default */
    int pos;            /* starting position in overall array of lag setters */
    int nvl;            /* number of siblings in array */
    int lmin;           /* minimum lag */
    int lmax;           /* maximum lag */
    char context;       /* LAG_X, LAG_W, ... */
    char *lspec;        /* string specification of particular lags */
    GtkWidget *spin1;   /* spinner for minimum lag */
    GtkWidget *spin2;   /* spinner for maximum lag */
    GtkWidget *entry;   /* text entry for specific lags */
    GtkWidget *toggle;  /* button to switch between spinners and entry */
    var_lag_info *vlp;  /* parent array */
};

#define depvar_row(c) (c == LAG_Y_X || c == LAG_Y_W)

/* handle action in the "specific lags" entry box */

static void lag_entry_callback (GtkWidget *w, gpointer p)
{
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));
    var_lag_info *vlinfo = (var_lag_info *) p;

    free(vlinfo->lspec);
    vlinfo->lspec = g_strdup(s);

    if (vlinfo->v == VDEFLT) {
	/* set the default for this context */
	var_lag_info *vlset = vlinfo->vlp;
	int i;

	for (i=vlinfo->pos+1; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context && vlset[i].entry != NULL) {
		gtk_entry_set_text(GTK_ENTRY(vlset[i].entry), s);
	    }
	}
    }
}

/* retrieve the min or max lag from a spinner and process
   the result */

static void lag_set_callback (GtkWidget *w, gpointer p)
{
    var_lag_info *vlinfo;
    int lag, *plag;

    vlinfo = (var_lag_info *) g_object_get_data(G_OBJECT(w), "vlinfo");

    plag = (w == vlinfo->spin1)? &vlinfo->lmin : &vlinfo->lmax;
    lag = *plag = spinner_get_int(w);

    /* force consistency if need be */

    if (w == vlinfo->spin1 && vlinfo->spin2 != NULL) {
	if (spinner_get_int(vlinfo->spin2) < lag) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo->spin2), lag);
	}
    } else if (w == vlinfo->spin2 && vlinfo->spin1 != NULL) {
	if (spinner_get_int(vlinfo->spin1) > lag) {
	    gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo->spin1), lag);
	}
    }	

    if (vlinfo->v == VDEFLT) {
	/* set the default value for this context */
	var_lag_info *vlset = vlinfo->vlp;
	GtkWidget *s;
	int i;

	for (i=vlinfo->pos; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context) {
		s = (w == vlinfo->spin1)? vlset[i].spin1 : vlset[i].spin2;
		if (s != NULL) {
		    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s), lag);
		}
	    }
	}
    }
}

/* switch from using min/max spinners to using the "specific lags"
   text entry, or vice versa */

static void activate_specific_lags (GtkWidget *w, var_lag_info *vlinfo)
{
    gboolean active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    if (active) {
	gtk_widget_set_sensitive(vlinfo->entry, TRUE);
	gtk_widget_set_sensitive(vlinfo->spin1, FALSE);
	gtk_widget_set_sensitive(vlinfo->spin2, FALSE);
    } else {
	gtk_widget_set_sensitive(vlinfo->entry, FALSE);
	gtk_widget_set_sensitive(vlinfo->spin1, TRUE);
	gtk_widget_set_sensitive(vlinfo->spin2, TRUE);
    } 

    if (vlinfo->v == VDEFLT) {
	/* set the default per context */
	var_lag_info *vlset = vlinfo->vlp;
	int i;

	for (i=vlinfo->pos+1; i<vlinfo->nvl; i++) {
	    if (vlset[i].context == vlinfo->context) {
		if (vlset[i].toggle != NULL) {
		    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vlset[i].toggle), 
						 active);
		}
	    }
	}
    } 
}

static void lag_toggle_register (GtkWidget *w, var_lag_info *vlinfo)
{
    var_lag_info *vlset = vlinfo->vlp;
    gboolean active;
    int i;

    if (vlset == NULL) return;

    for (i=0; i<vlinfo->nvl; i++) {
	if (depvar_row(vlset[i].context)) {
	    if (!GTK_WIDGET_SENSITIVE(vlset[i].spin1) &&
		!GTK_WIDGET_SENSITIVE(vlset[i].entry)) {
		/* dependent var lags are disabled */
		vlset[i].lmin = vlset[i].lmax = NOT_LAG; /* ?? */
		free(vlset[i].lspec);
		vlset[i].lspec = NULL;
		continue;
	    }
	}
	active = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(vlset[i].toggle));
	if (active) {
	    vlset[i].lmin = vlset[i].lmax = NOT_LAG; 
	} else {
	    free(vlset[i].lspec);
	    vlset[i].lspec = NULL;
	}
    }
}

static void activate_y_lags (GtkWidget *w, var_lag_info *vlinfo)
{
    gboolean active = 
	gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    gtk_widget_set_sensitive(vlinfo->spin1, active);
    gtk_widget_set_sensitive(vlinfo->spin2, active);
    gtk_widget_set_sensitive(vlinfo->toggle, active);

    if (active) {
	vlinfo->lmin = spinner_get_int(vlinfo->spin1);
	vlinfo->lmax = spinner_get_int(vlinfo->spin2);
    }

    if (vlinfo->context == LAG_Y_X) {
	y_x_lags_enabled = active;
    } else {
	y_w_lags_enabled = active;
    }
}

static void lagsel_spin_connect (GtkWidget *button)
{
    g_signal_connect(G_OBJECT(button), "value-changed", 
		     G_CALLBACK(lag_set_callback), NULL);
    gtk_entry_set_activates_default(GTK_ENTRY(button), TRUE);
}

static void lags_set_OK (GtkWidget *w, gpointer p)
{
    int *resp = (int *) p;

    *resp = GRETL_YES;
}

static void resensitive_selector (GtkWidget *w, gpointer p)
{
    if (open_selector != NULL) {
	gtk_widget_set_sensitive(open_selector->dlg, TRUE);
    }
}

/* The actual lag selection dialog: we provide spinners for a lag
   range and also a free-form entry field for non-contiguous lags.  In
   some circumstances we allow specification of lags for the dependent
   variable as well as the independent vars.
*/

static int 
lags_dialog (const int *list, var_lag_info *vlinfo, selector *sr) 
{
    GtkWidget *lbl, *dialog, *myvbox;
    GtkWidget *tbl, *tmp, *hbox;
    GtkWidget *y_check = NULL;
    gint tbl_len;
    double lmax, ldef;
    int VAR_special, insts;
    int i, j;
    int ret = GRETL_CANCEL;

    dialog = gretl_dialog_new(_("lag order"), sr->dlg, 
			      GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(resensitive_selector), NULL);

    myvbox = gtk_vbox_new(FALSE, 5);

    VAR_special = (vlinfo[0].v == VDEFLT && vlinfo[0].context == LAG_Y_V);
    insts = in_gretl_list(list, LAG_W);

    lmax = (datainfo->t2 - datainfo->t1) / list[0];
    ldef = datainfo->pd;

    if (VAR_special) {
	/* allow for gaps in VAR lag order */
	lbl = gtk_label_new(_("Lags of endogenous variables"));
	gtk_box_pack_start(GTK_BOX(myvbox), lbl, FALSE, FALSE, 0);
    } 

    /* allow for additional label row */
    tbl_len = list[0] + 1;

    tbl = gtk_table_new(tbl_len, 7, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(myvbox), tbl, FALSE, FALSE, 0);

    /* row 0 of table: heading(s) */

    if (insts) {
	/* there's an instruments section to follow */
	lbl = gtk_label_new(_("Regressors"));
	gtk_table_attach(GTK_TABLE(tbl), lbl, 0, 7, 0, 1, 0, 0, 0, 5);
    } else {
	if (!VAR_special) {
	    lbl = gtk_label_new(_("Variable"));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, 0, 1);
	    lbl = gtk_label_new(_("lags"));
	    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 1, 4, 0, 1);
	    lbl = gtk_label_new("  ");
	    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 4, 5, 0, 1);
	}

	lbl = gtk_label_new(_("or"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 5, 6, 0, 1);
	lbl = gtk_label_new(_("specific lags"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 6, 7, 0, 1);  
    }

    j = 0;
    for (i=1; i<=list[0]; i++) {
	var_lag_info *vlj;
	int li = list[i];
	int lmin = 0;

	if (list_lag_special(li)) {
	    if (li == LAG_W) {
		tmp = gtk_label_new(_("Instruments"));
		gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 7, i, i+1);
	    } else if (depvar_row(li)) { 
		y_check = gtk_check_button_new_with_label(_("Lags of dependent variable"));
		gtk_table_attach_defaults(GTK_TABLE(tbl), y_check, 0, 7, i, i+1);
	    } else {
		tmp = gtk_hseparator_new();
		gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 7, i, i+1);
	    }
	    continue;
	}

	if (!VAR_special) {
	    if (li == VDEFLT) {
		lbl = gtk_label_new(_("default"));
	    } else {
		lbl = gtk_label_new(datainfo->varname[li]);
	    }
	    gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, i, i+1);
	}

	vlj = &vlinfo[j++];

	if (depvar_row(vlj->context) || vlj->context == LAG_Y_V) {
	    lmin = 1;
	} 

	/* min. lag spinner */
	vlj->spin1 = gtk_spin_button_new_with_range(lmin, lmax, 1);
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlj->spin1, 1, 2, i, i+1);
	g_object_set_data(G_OBJECT(vlj->spin1), "vlinfo", vlj);
	lagsel_spin_connect(vlj->spin1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlj->spin1), vlj->lmin);

	lbl = gtk_label_new(_("to"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, i, i+1);

	/* max. lag spinner */
	vlj->spin2 = gtk_spin_button_new_with_range(lmin, lmax, 1);
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlj->spin2, 3, 4, i, i+1);
	g_object_set_data(G_OBJECT(vlj->spin2), "vlinfo", vlj);
	lagsel_spin_connect(vlj->spin2);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlj->spin2), vlj->lmax);

	/* spacer column */
	lbl = gtk_label_new("  ");
	gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 4, 5, i, i+1);

	/* toggle button for activating entry of specific lags */
	vlj->toggle = gtk_check_button_new();
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlj->toggle, 5, 6, i, i+1);
	g_signal_connect(G_OBJECT(vlj->toggle), "toggled", 
			 G_CALLBACK(activate_specific_lags), vlj);

	/* text entry widget for specific lags */
	vlj->entry = gtk_entry_new();
	gtk_entry_set_width_chars(GTK_ENTRY(vlj->entry), 16);
	gtk_table_attach_defaults(GTK_TABLE(tbl), vlj->entry, 6, 7, i, i+1);
	g_signal_connect(G_OBJECT(vlj->entry), "changed", 
			 G_CALLBACK(lag_entry_callback), vlj);
	gtk_entry_set_activates_default(GTK_ENTRY(vlj->entry), TRUE);
	if (vlj->lspec != NULL && *vlj->lspec != '\0') {
	    /* got a saved non-contiguous lag spec string: apply it */
	    gtk_entry_set_text(GTK_ENTRY(vlj->entry), vlj->lspec);
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(vlj->toggle), 
					 TRUE);
	} else {
	    /* disable this by default */
	    gtk_widget_set_sensitive(vlj->entry, FALSE);
	}

	/* set sensitivity of dependent variable lags apparatus */
	if (y_check != NULL && depvar_row(vlj->context)) {
	    g_signal_connect(G_OBJECT(y_check), "toggled", 
			     G_CALLBACK(activate_y_lags), vlj);
	    if ((vlj->context == LAG_Y_X && y_x_lags_enabled) ||
		(vlj->context == LAG_Y_W && y_w_lags_enabled)) {
		gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(y_check),
					     TRUE);
	    } else {
		gtk_widget_set_sensitive(vlj->spin1, FALSE);
		gtk_widget_set_sensitive(vlj->spin2, FALSE);
		gtk_widget_set_sensitive(vlj->toggle, FALSE);
		gtk_widget_set_sensitive(vlj->entry, FALSE);
	    }
	}
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), myvbox, TRUE, TRUE, 5);

    if (list[0] > 10) {
	GtkWidget *scroller = gtk_scrolled_window_new(NULL, NULL);

	gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
				       GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller), 
					      hbox);
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), scroller, TRUE, TRUE, 5);
	gtk_widget_show_all(scroller);
	gtk_widget_set_size_request(scroller, -1, 360);
    } else {
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
	gtk_widget_show_all(hbox);
    }
    
    hbox = GTK_DIALOG(dialog)->action_area;
	
    /* "Cancel" button */
    tmp = cancel_delete_button(hbox, dialog, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);

    /* "OK" button */
    tmp = ok_button(hbox);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(lags_set_OK), &ret);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lag_toggle_register), &vlinfo[0]);
    g_signal_connect(G_OBJECT(tmp), "clicked", 
		     G_CALLBACK(delete_widget), dialog);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    /* "Help" button */
    context_help_button(hbox, LAGS_DIALOG);

    gtk_widget_set_sensitive(sr->dlg, FALSE);
    gtk_widget_show(dialog);   

    return ret;
}

/* based on the lags dialog, revise the representation of a given
   variable within the selection boxes in the model specification
   dialog, using varlist_insert_var_full().
*/

static void 
revise_var_string (var_lag_info *vlinfo, selector *sr, int locus)
{
    GtkWidget *w;
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int modv;

    w = (locus == SR_RVARS2)? sr->rvars2 : sr->rvars1;
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
    gtk_tree_model_get_iter_first(mod, &iter);

    do {
	gtk_tree_model_get(mod, &iter, 0, &modv, -1);
	if (modv == vlinfo->v) {
#if VLDEBUG
	    fprintf(stderr, "revise_var_string: calling varlist_insert_var_full\n");
#endif
	    varlist_insert_var_full(vlinfo->v, mod, &iter, sr, locus);
	    break;
	}
    } while (gtk_tree_model_iter_next(mod, &iter));
}

/* go through a list box and find the vars that are candidates for
   lag selection: this excludes constant and trend, and also
   list entries that are themselves just dummy placeholders
   for already selected lags of a variable */

static int get_laggable_vars (GtkWidget *w, int context, int *list, int *i)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint v, lag;
    int n = 0;

    if (!GTK_IS_TREE_VIEW(w)) {
	return 0;
    }

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
    if (model == NULL) {
	return 0;
    }

    if (gtk_tree_model_get_iter_first(model, &iter)) {
	do {
	    int laggable = 1;

	    gtk_tree_model_get(model, &iter, 0, &v, 1, &lag, -1);

	    if (v == 0 || !strcmp(datainfo->varname[v], "time") ||
		!strcmp(datainfo->varname[v], "timesq")) {
		laggable = 0;
	    } else if (is_lag_dummy(v, lag, context)) {
		laggable = 0;
	    }

	    if (laggable) { 
		if (list != NULL) {
		    list[*i] = v;
		    *i += 1;
		}
		n++;
	    }
	} while (gtk_tree_model_iter_next(model, &iter));
    }

    return n;
}

/* Construct a special list which encodes the information needed
   for building the lag selection dialog
*/

static int *sr_get_stoch_list (selector *sr, int *pnset, int *pcontext)
{
    GtkWidget *listw[2] = { NULL, NULL };
    gint ynum = 0;
    int nv[2] = {0};
    int nset, nsep;
    int context;
    int i, j;
    int *slist = NULL;

    if (sr->ci != ARMA && sr->ci != VAR &&
	sr->ci != VECM && sr->ci != VLAGSEL &&
	sr->ci != ARBOND) { 
	ynum = selector_get_depvar_number(sr);
    }

    listw[0] = (select_lags_primary(sr->ci))? sr->rvars1 : sr->rvars2;

    if (USE_ZLIST(sr->ci)) {
	listw[1] = sr->rvars2;
    } 

    if (listw[0] != NULL) {
	nv[0] = get_laggable_vars(listw[0], LAG_X, NULL, NULL);
	if (nv[0] == 0) {
	    listw[0] = NULL;
	}
    }

    if (listw[1] != NULL) {
	nv[1] = get_laggable_vars(listw[1], LAG_W, NULL, NULL);
	if (nv[1] == 0) {
	    listw[1] = NULL;
	}	
    }    

    nset = nsep = 0;

#if LDEBUG
    fprintf(stderr, "sr_get_stoch_list: ynum = %d, nv[0] = %d, nv[1] = %d\n",
	    ynum, nv[0], nv[1]);
#endif

    if (ynum < 0 && nv[0] == 0 && nv[1] == 0) {
	/* no vars to deal with */
	errbox("Please add some variables to the model first");
	return NULL;
    }

    /* first pass: figure out how many elements the list should have */

    if (nv[0] > 0) {
	*pcontext = LAG_X;
	nset += nv[0] + 1; /* the Xs plus their defaults */
	if (ynum > 0) {
	    nsep++; /* depvar heading */
	    nset++; /* depvar row */
	}
    }

    if (nv[1] > 0) {
	if (nv[0] > 0) {
	    nsep++; /* separator from Xs */
	} else {
	    *pcontext = LAG_W;
	}
	nsep++;            /* "instruments" heading */
	nset += nv[1] + 1; /* instruments plus their defaults */
	if (ynum > 0) {
	    nsep++; /* depvar heading */
	    nset++; /* depvar row */
	}	    
    }

    if (nv[0] == 0 && nv[1] == 0) {
	/* only the dependent variable is present */
	*pcontext = LAG_Y_X;
	nsep++; /* depvar heading */
	nset++; /* depvar row */
	if (USE_ZLIST(sr->ci)) {
	    nsep++; /* list separator for insts */
	    nsep++; /* "instruments" heading */
	    nsep++; /* depvar heading */
	    nset++; /* depvar row */
	}
    }

    /* allocate the list */

    slist = gretl_list_new(nset + nsep);
    if (slist == NULL) {
	return NULL;
    }

    /* second pass: actually build the list, inserting special
       separators if needed
    */

    i = 1;
    for (j=0; j<2; j++) {
	if (listw[j] == NULL) {
	    continue;
	}
	context = (j == 1)? LAG_W : LAG_X;
	if (j == 1) {
	    if (nv[0] > 0) {
		slist[i++] = LISTSEP;
	    }
	    slist[i++] = LAG_W;
	}
	slist[i++] = VDEFLT;
	get_laggable_vars(listw[j], context, slist, &i);
	if (ynum > 0) {
	    slist[i++] = (j > 0)? LAG_Y_W : LAG_Y_X;
	    slist[i++] = ynum;
	}
    }

    /* special case where the dependent variable is the only laggable
       variable selected so far */

    if (nv[0] == 0 && nv[1] == 0) {
	slist[i++] = LAG_Y_X;
	slist[i++] = ynum;
	if (USE_ZLIST(sr->ci)) {
	    slist[i++] = LISTSEP;
	    slist[i++] = LAG_W;
	    slist[i++] = LAG_Y_W;
	    slist[i++] = ynum;
	}
    }

    *pnset = nset;

    return slist;
}

static void maybe_revise_var_string (var_lag_info *vlinfo, selector *sr)
{
    int locus = 0;

#if VLDEBUG
    fprintf(stderr, "maybe_revise_var_string: v = %d, context = %d\n",
	    vlinfo->v, vlinfo->context);
#endif

    if (vlinfo->context == LAG_X) {
	locus = (sr->ci == VAR || 
		 sr->ci == VECM ||
		 sr->ci == VLAGSEL)? SR_RVARS2 : SR_RVARS1;
    } else if (vlinfo->context == LAG_W) {
	locus = SR_RVARS2;
    } else if (depvar_row(vlinfo->context)) { 
#if VLDEBUG
	fprintf(stderr, " calling maybe_revise_depvar_lags, context = %d\n",
		vlinfo->context);
#endif
	maybe_revise_depvar_lags(sr, vlinfo->v, vlinfo->context);
    }

    if (locus > 0) {
#if VLDEBUG
	fprintf(stderr, " calling revise_var_string, locus = %d\n", locus);
#endif
	revise_var_string(vlinfo, sr, locus);
    } 
}

/* return 1 if lags changed, else 0 */

static int set_lags_for_var (var_lag_info *vlinfo, int yxlags, int ywlags)
{
    int *llist = NULL;
    int changed = 0;
    int err = 0;

#if LDEBUG
    fprintf(stderr, "set_lags_for_var: v=%d, lmin=%d, lmax=%d, lspec=%p\n",
	    vlinfo->v, vlinfo->lmin, vlinfo->lmax, (void *) vlinfo->lspec);
#endif

    if (vlinfo->lspec != NULL && *vlinfo->lspec != 0) {
	charsub(vlinfo->lspec, ',', ' ');
	llist = gretl_list_from_string(vlinfo->lspec, &err);
	if (!err) {
	    err = set_lag_prefs_from_list(vlinfo->v, llist, vlinfo->context,
					  &changed);
	    if (err) {
		free(llist);
	    }
	}
    } else if (vlinfo->lmin != NOT_LAG && vlinfo->lmax != NOT_LAG) {
	set_lag_prefs_from_minmax(vlinfo->v, vlinfo->lmin, vlinfo->lmax, 
				  vlinfo->context, &changed);
    } else {
	set_null_lagpref(vlinfo->v, vlinfo->context, &changed);
    }

    if (vlinfo->context == LAG_Y_X && yxlags != y_x_lags_enabled) {
	changed = 1;
    } else if (vlinfo->context == LAG_Y_W && ywlags != y_w_lags_enabled) {
	changed = 1;
    }

    return changed;
}

/* Keep the global "lag order" spinner in sync with lag specifications
   entered via the the (more complex) lags dialog.  We set the value
   shown in the spinner to the maximum lag; in addition we desensitize
   the global spinner if the user has specified gaps in the lag
   structure.
*/

static void sync_lag_order_spinner (var_lag_info *vlinfo,
				    selector *sr)
{
    const int *list = NULL;
    int lmin = 0, lmax = 0;

    get_lag_preference(vlinfo->v, &lmin, &lmax, &list,
		       vlinfo->context, sr);

    if (list != NULL) {
	lmin = list[1];
	lmax = list[list[0]];
    }

    if (lmax != 0 && lmax != spinner_get_int(sr->extra[0])) {
	/* update max lag as shown by spinner */
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[0]),
				  lmax);
    }

    if (lmin > 1 || list != NULL) {
	/* not 1-based consecutive lags: so disable spinner */
	gtk_widget_set_sensitive(sr->extra[0], FALSE);
    } else {
	gtk_widget_set_sensitive(sr->extra[0], TRUE);
    }
}

#if LDEBUG > 2
static void print_vlinfo (var_lag_info *vlj)
{
    fprintf(stderr, "\nCreated vlinfostruct:\n");
    fprintf(stderr, " v = %d\n", vlj->v);
    fprintf(stderr, " pos = %d\n", vlj->pos);
    fprintf(stderr, " nvl = %d\n", vlj->nvl);
    fprintf(stderr, " context = %d\n", (int) vlj->context);
    fprintf(stderr, " vlp = %p\n", (void *) vlj->vlp);
    
    fprintf(stderr, " lmin = %d\n", vlj->lmin);
    fprintf(stderr, " lmax = %d\n", vlj->lmax);

    if (vlj->lspec != NULL) {
	fprintf(stderr, " lspec = '%s'\n", vlj->lspec);
    } else {
	fprintf(stderr, " lspec = NULL\n");
    }
}
#endif

/* Respond to the user clicking the "Lags..." button: here we build
   the information needed to make the lag selection dialog, run
   the dialog, then implement any changes made via the GUI.
 */

static gboolean lags_dialog_driver (GtkWidget *w, selector *sr)
{
    var_lag_info *vlinfo;
    int yxlags = y_x_lags_enabled;
    int ywlags = y_w_lags_enabled;
    int context = 0;
    int i, j, resp, nvl;
    int *list;

    if (w == sr->extra[EXTRA_LAGS]) {
	/* VAR: use a dummy entry for all endogenous vars */
	list = gretl_list_new(1);
	if (list == NULL) {
	    return FALSE;
	}
	list[1] = VDEFLT;
	nvl = 1;
	context = LAG_Y_V;
    } else {
	/* get the list of stochastic (laggable) variables */
	list = sr_get_stoch_list(sr, &nvl, &context);
	if (list == NULL) {
	    return FALSE;
	}
    }

#if LDEBUG
    printlist(list, "stochastic vars list");
    fprintf(stderr, "number of setters = %d\n", nvl);
#endif

    /* allocate array of vlinfo, one per laggable var */
    vlinfo = mymalloc(nvl * sizeof *vlinfo);
    if (vlinfo == NULL) {
	free(list);
	return FALSE;
    }

    j = 0;
    for (i=1; i<=list[0]; i++) {
	const int *laglist = NULL;
	var_lag_info *vlj;
	int vi = list[i];

	if (list_lag_special(vi)) {
	    /* LISTSEP is used to switch "lag context" */
	    context = (vi == LISTSEP)? 0 : vi;
#if LDEBUG
	    fprintf(stderr, "list[%d] = %d: set context = %d\n",
		    i, vi, context);
#endif
	    continue;
	}

	vlj = &vlinfo[j++];

	/* pick up any saved preferences (including saved defaults) */
	get_lag_preference(vi, &vlj->lmin, &vlj->lmax, 
			   &laglist, context, sr);

	if (laglist != NULL) {
	    /* we got a list of specific lags for variable vi: convert
	       to a string for putting into a text entry box */
	    vlj->lspec = gretl_list_to_string(laglist);
	    vlj->lmin = laglist[1];
	    vlj->lmax = laglist[laglist[0]];
	} else {
	    vlj->lspec = NULL;
	}

	vlj->v = vi;
	vlj->pos = j;
	vlj->nvl = nvl;
	vlj->context = context;
	vlj->spin1 = NULL;
	vlj->spin2 = NULL;
	vlj->entry = NULL;
	vlj->toggle = NULL;
	vlj->vlp = vlinfo;
#if LDEBUG > 2
	print_vlinfo(vlj);
#endif
    }

    /* show the user the lags dialog box */
    resp = lags_dialog(list, vlinfo, sr);

    /* register any changes made by the user */
    if (resp != GRETL_CANCEL) {
	for (j=0; j<nvl; j++) {
	    int changed = set_lags_for_var(&vlinfo[j], yxlags, ywlags);

	    if (vlinfo[j].v != VDEFLT && changed) {
		maybe_revise_var_string(&vlinfo[j], sr);
	    }
	}
	if (context == LAG_Y_V) {
	    sync_lag_order_spinner(&vlinfo[0], sr);
	}
    }

    /* clean up any allocated lag spec strings */
    for (j=0; j<nvl; j++) {
	if (vlinfo[j].lspec != NULL) {
	    free(vlinfo[j].lspec);
	}
    }

    free(list);
    free(vlinfo);

    return FALSE;
}
