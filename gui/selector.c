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
#include "toolbar.h"
#include "winstack.h"
#include "tabwin.h"
#include "lagpref.h"
#include "fncall.h"
#include "textbuf.h"

#include "var.h"
#include "gretl_func.h"
#include "libset.h"
#include "uservar.h"
#include "johansen.h"
#include "gretl_bfgs.h"
#include "gretl_midas.h"

/* for graphical selector buttons */
#include "arrows.h"

#if 0
# include "bootstrap.h"
#endif

#define LDEBUG 0
#define VLDEBUG 0

/* codes for window components */
enum {
    SR_LVARS  = 1,
    SR_RVARS1,
    SR_RVARS2,
    SR_LVARS2
};

/* state flags */
enum {
    SR_BLOCKING     = 1 << 0,
    SR_STATE_PUSHED = 1 << 1,
    SR_LAGS_HIDDEN  = 1 << 2,
    SR_NO_SINGLE    = 1 << 3
};

#define blocking(s)     (s->flags & SR_BLOCKING)
#define state_pushed(s) (s->flags & SR_STATE_PUSHED)
#define lags_hidden(s)  (s->flags & SR_LAGS_HIDDEN)

#define N_EXTRA 8

struct _selector {
    GtkWidget *parent;
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *action_area;
    GtkWidget *lvars;
    GtkWidget *lvars2;
    GtkWidget *table;
    GtkWidget *depvar;
    GtkWidget *rvars1;
    GtkWidget *rvars2;
    GtkWidget *default_check;
    GtkWidget *add_button;
    GtkWidget *remove_button;
    GtkWidget *lags_button;
    GtkWidget *hess_button;
    GtkWidget *x12a_button;
    GtkWidget *xdiff_button;
    GtkWidget *hccme_button;
    GtkWidget *extra[N_EXTRA];
    int ci;
    int flags;
    int active_var;
    int error;
    int n_left;
    int row;
    int n_rows;
    gretlopt opts;
    char *cmdlist;
    size_t cmdsize;
    size_t cmdlen;
    gpointer data;
    gpointer extra_data;
    sr_callback callback;
};

enum {
    COL_ID = 0,
    COL_LAG,
    COL_NAME,
    COL_FLAG
};

enum {
    MCOL_M = 1,
    MCOL_NAME,
    MCOL_MINLAG,
    MCOL_MAXLAG,
    MCOL_TYPE,
    MCOL_K
};

enum {
    ARMA_p,
    ARIMA_d,
    ARMA_q,
    ARMA_plist,
    ARMA_qlist,
    ARMA_P,
    ARIMA_D,
    ARMA_Q
};

enum {
    REGLS_EST,
    REGLS_ALPHA,
    REGLS_LAMVAL,
    REGLS_NLAM,
    REGLS_NFOLDS,
    REGLS_FTYPE
};

enum {
    OMIT_B3,
    OMIT_ALPHA,
    OMIT_B4,
    OMIT_IC,
    OMIT_SEL
};

#define EXTRA_LAGS (N_EXTRA - 1)

#define VNAME_WIDTH 18
#define BUTTON_WIDTH 64

/* single-equation estimation commands plus some GUI extensions */
#define MODEL_CODE(c) (MODEL_COMMAND(c) || c == PANEL_WLS || c == PANEL_B || \
                       c == OLOGIT || c == OPROBIT || c == REPROBIT ||  \
                       c == MLOGIT || c == IV_LIML || c == IV_GMM ||    \
                       c == COUNTMOD || c == REGLS || c == FE_LOGISTIC || \
                       c == ALAGSEL)

#define ARMA_RELATED(c) (c == ARMA || c == ALAGSEL)

#define NO_X_OK(c) (c == ARMA || c == ALAGSEL || c == GARCH)

#define IV_MODEL(c) (c == IVREG || c == IV_LIML || c == IV_GMM)

#define NONPARAM_CODE(c) (c == LOESS || c == NADARWAT)

#define COINT_CODE(c) (c == COINT || c == COINT2)

#define VEC_CODE(c) (c == COINT || c == COINT2 || c == VAR ||   \
                     c == VECM || c == VLAGSEL)

#define VEC_MODEL_CODE(c) (c == VAR || c == VECM || c == VLAGSEL)

#define VECLAGS_CODE(c) (c == VAR || c == VECM)

#define ADDVAR_CODE(c) (c == LOGS || c == LAGS || c == SQUARE ||        \
                        c == DIFF || c == LDIFF)

#define TWO_VARS_CODE(c) (c == ELLIPSE || c == XCORRGM || c == QQPLOT)

#define THREE_VARS_CODE(c) (c == GR_DUMMY || c == GR_XYZ ||     \
                            c == GR_3D || c == ANOVA)

#define FNPKG_CODE(c) (c == SAVE_FUNCTIONS || c == EDIT_FUNCTIONS)

#define SHOW_LISTS_CODE(c) (c == SUMMARY || c == CORR || c == MAHAL || c == PCA)

#define LIST_USE_INTS(c) (c == ELLIPSE || c == SAVE_FUNCTIONS || c == REGLS_PLOTSEL)

#define WANT_TOGGLES(c) (c == DPANEL ||         \
                         c == ARMA ||           \
                         c == COINT ||          \
                         c == COINT2 ||         \
                         c == CORR ||           \
                         c == GARCH ||          \
                         c == HECKIT ||         \
                         c == HSK ||            \
                         c == BIPROBIT ||       \
                         c == INTREG ||         \
                         c == IVREG ||          \
                         c == LOGIT ||          \
                         c == OLOGIT ||         \
                         c == MLOGIT ||         \
                         c == LOGISTIC ||       \
                         c == MPOLS ||          \
                         c == OLS ||            \
                         c == PANEL ||          \
                         c == PANEL_WLS ||      \
                         c == PANEL_B ||        \
                         c == FE_LOGISTIC ||    \
                         c == COUNTMOD ||       \
                         c == DURATION ||       \
                         c == PROBIT ||         \
                         c == OPROBIT ||        \
                         c == REPROBIT ||       \
                         c == QUANTREG ||       \
                         c == MIDASREG ||       \
                         c == TOBIT ||          \
                         c == VAR ||            \
                         c == VECM ||           \
                         c == VLAGSEL ||        \
                         c == ALAGSEL ||        \
                         c == WLS ||            \
                         c == GR_BOX ||         \
			 c == SUMMARY ||        \
			 c == FSUMMARY ||       \
                         c == XTAB)

#define USE_VECXLIST(c) (c == VAR || c == VLAGSEL || c == VECM ||       \
                         c == COINT2)

#define USE_RXLIST(c) (c == VECM || c == COINT2)

#define AUX_LAST(c) (c == IVREG ||              \
                     c == IV_LIML ||            \
                     c == IV_GMM ||             \
                     c == HECKIT ||             \
                     c == BIPROBIT ||           \
                     c == VAR ||                \
                     c == VLAGSEL ||            \
                     c == VECM ||               \
                     c == COINT2 ||             \
                     c == MIDASREG ||           \
                     c == SAVE_FUNCTIONS ||     \
                     c == EDIT_FUNCTIONS)

#define USE_ZLIST(c) (c == IVREG || c == IV_LIML || c == IV_GMM ||      \
                      c == HECKIT || c == BIPROBIT)

#define RHS_PREFILL(c) (c == CORR ||            \
                        c == MAHAL ||           \
                        c == PCA ||             \
                        c == SUMMARY ||         \
                        c == XTAB)

#define dataset_lags_ok(d) ((d)->structure == TIME_SERIES ||            \
                            (d)->structure == SPECIAL_TIME_SERIES ||    \
                            (d)->structure == STACKED_TIME_SERIES)

#define select_lags_primary(c) (MODEL_CODE(c))

#define select_lags_depvar(c) (MODEL_CODE(c) && c != ARMA &&    \
                               c != DPANEL && c != MIDASREG)

/* Should we have a lags button associated with auxiliary
   variable selector? */

#define select_lags_aux(c) (c == VAR || c == VLAGSEL || c == VECM ||    \
                            c == IVREG || c == IV_LIML || c == IV_GMM || \
                            c == HECKIT || c == BIPROBIT)

#define list_lag_special(i) (i < -1)

typedef void (*click_func) (GtkWidget *, selector *);

/* static state variables */

static int default_y = -1;
static int want_seasonals = 0;
static int default_order;
static int vartrend = 0;
static int varconst = 1;
static int arma_p = 1;
static int arma_P = 0;
static int arima_d = 0;
static int arima_D = 0;
static int arma_q = 1;
static int arma_Q = 0;
static int garch_p = 1;
static int garch_q = 1;
static int arma_const = 1;
static int garch_const = 1;
static int arma_x12 = 0;
static int arma_hessian = 1;
static int arima_xdiff = 1;
static int selvar;
static int y2var;
static int offvar;
static int censvar;
static int jrank = 1;
static int jcase = J_UNREST_CONST;
static int verbose;
static int lovar;
static int hivar;
static int wtvar;
static int np_xvar;
static char lp_pvals;

static char dpd_2step;
static char dpd_asy;
static char dpd_dpd;
static char dpd_p;
static char dpd_coll;

static int y_x_lags_enabled;
static int y_w_lags_enabled;

static int *xlist;
static int *instlist;
static int *veclist;
static int *vecxlist;

static char *arlags;
static char *malags;
static char cluster_var[VNAMELEN];

static double tobit_lo = 0;
static double tobit_hi = NADBL;

static char mds_listname[VNAMELEN];
static int mds_quad[4] = {0, 0, 1, 2};
static int mds_order = 1;

static gretlopt model_opt;

static GtkWidget *multiplot_label;
static GtkWidget *multiplot_menu;

static gretl_bundle *regls_adv;

static selector *open_selector;

static gint listvar_reorder_click (GtkWidget *widget, GdkEventButton *event,
                                   gpointer data);
static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event,
                               selector *sr);
static gint listvar_flagcol_click (GtkWidget *widget, GdkEventButton *event,
                                   gpointer data);
static gint listvar_midas_click (GtkWidget *widget, GdkEventButton *event,
                                 selector *sr);
static int list_show_var (selector *sr, int v, int show_lags, gchar **exclude);
static void available_functions_list (selector *sr);
static void primary_rhs_varlist (selector *sr);
static gboolean lags_dialog_driver (GtkWidget *w, selector *sr);
static void call_iters_dialog (GtkWidget *w, GtkWidget *combo);
static void reset_arma_spins (selector *sr);
static void clear_midas_spec (void);
static int check_midas_rvars2 (GtkTreeModel *model, gboolean *have_beta1);

static int set_or_get_n_rvars1 (selector *sr, int n)
{
    static int nv;

    if (n >= 0) {
        nv = n;
        if (sr != NULL && sr->ci == CORR && sr->extra[0] != NULL) {
            gtk_widget_set_sensitive(sr->extra[0], n > 2);
        }
    }

    return nv;
}

static void set_n_rvars1 (selector *sr, int n)
{
    set_or_get_n_rvars1(sr, n);
}

static int get_n_rvars1 (void)
{
    return set_or_get_n_rvars1(NULL, -1);
}

static int want_combo (selector *sr)
{
    return (sr->ci == ARMA ||
            sr->ci == ALAGSEL ||
            sr->ci == VECM ||
            sr->ci == COINT ||
            sr->ci == COINT2 ||
            sr->ci == COUNTMOD ||
            sr->ci == DURATION ||
            sr->ci == IV_GMM ||
            sr->ci == EXPORT);
}

static int want_radios (selector *sr)
{
    int c = sr->ci;
    int ret = 0;

    if (c == PANEL || c == SCATTERS || c == AR1 ||
        c == LOGIT || c == PROBIT || c == HECKIT ||
        c == XTAB || c == PCA ||
        c == QUANTREG || c == DPANEL ||
        c == LOGISTIC || c == FE_LOGISTIC ||
        c == VAROMIT || c == REGLS) {
        ret = 1;
    } else if (c == ADD || c == OMIT) {
        windata_t *vwin = (windata_t *) sr->data;
        MODEL *pmod = (MODEL *) vwin->data;

        if (c == OMIT && pmod->ci == HECKIT) {
            /* omit: Wald test only */
            sr->opts |= OPT_W;
        } else if (c == ADD && pmod->ci != OLS) {
            /* add: LM and auto options are OLS-only */
            ret = 0;
        } else {
            ret = 1;
        }
    }

    return ret;
}

static void selector_set_blocking (selector *sr, int modal)
{
    if (modal) {
        gretl_set_window_modal(sr->dlg);
    }
    sr->flags |= SR_BLOCKING;
    gtk_main();
}

static int selection_at_max (selector *sr, GtkWidget *w, int nsel)
{
    if (TWO_VARS_CODE(sr->ci) && nsel == 2) {
	return 1;
    } else if (w != NULL) {
	int selmax = widget_get_int(w, "selmax");

	if (selmax > 0 && nsel == selmax) {
	    return 1;
	}
    }

    return 0;
}

static int sr_get_lag_context (selector *sr, int locus)
{
    int c = 0;

    if (sr == NULL || !dataset_lags_ok(dataset)) {
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

/* note: @nx is exclusive of const and (lags of) dependent variable */

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
    y2var = 0;
    offvar = 0;
    censvar = 0;
    wtvar = 0;
    np_xvar = 0;
    vartrend = 0;
    varconst = 1;
    lovar = hivar = 0;

    dpd_asy = dpd_2step = 0;
    dpd_dpd = dpd_coll = 0;
    dpd_p = 1;

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
    garch_const = 1;

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

    *cluster_var = '\0';

    tobit_lo = 0;
    tobit_hi = NADBL;

    lp_pvals = 0;

    clear_midas_spec();
    destroy_lag_preferences();
    call_iters_dialog(NULL, NULL);

    if (open_selector != NULL) {
        selector *sr = open_selector;

        if (sr->ci == ARMA) {
            reset_arma_spins(sr);
        }
    }
}

GtkWidget *selector_get_window (const selector *sr)
{
    return (sr != NULL)? sr->dlg : NULL;
}

static int presel;

void selector_set_varnum (int v)
{
    presel = v;
}

void modelspec_dialog (int ci)
{
    selection_dialog(ci, _("gretl: specify model"), NULL, do_model);
}

static int varnum_from_keystring (MODEL *pmod, const char *key)
{
    char *s = (char *) gretl_model_get_data(pmod, key);
    int v = 0;

    if (s != NULL && *s != '\0') {
        v = series_index(dataset, s);
        if (v == dataset->v) {
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

    if (!(acode & ARMA_EXACT)) {
        model_opt |= OPT_C;
    }

    if (acode & ARMA_X12A) {
        arma_x12 = 1;
    }

    if (pmod->opt & OPT_G) {
        /* use OPG for covariance matrix */
        arma_hessian = 0;
    }

    if (pmod->opt & OPT_Y) {
        /* don't difference ARIMAX regressors */
        arima_xdiff = 0;
    }

    if (pmod->opt & OPT_L) {
        model_opt |= OPT_L;
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
                arlags = gretl_list_to_numeric_string(laglist);
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
                malags = gretl_list_to_numeric_string(laglist);
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
        arlags = gretl_list_to_numeric_string(pmod->arinfo->arlist);
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

static void retrieve_biprobit_info (MODEL *pmod, int *gotinst)
{
    /* the biprobit list takes the form:
       y1 y2 x1 ; x2
    */
    if (pmod->list != NULL && pmod->list[0] > 2) {
        y2var = pmod->list[2];
    }
    free(instlist);
    instlist = gretl_model_get_secondary_list(pmod);
    *gotinst = instlist != NULL;
}

static void retrieve_tobit_info (MODEL *pmod)
{
    double x0 = gretl_model_get_double(pmod, "llimit");
    double x1 = gretl_model_get_double(pmod, "rlimit");

    if (na(x0) && na(x1)) {
        /* the standard set-up */
        tobit_lo = 0;
        tobit_hi = NADBL;
    } else {
        /* user-specified limits */
        tobit_lo = x0;
        tobit_hi = x1;
    }
}

static void retrieve_midas_info (MODEL *pmod)
{
    gretl_bundle *b = NULL;
    gretl_array *A;
    int err = 0;

    if (pmod->opt & OPT_L) {
        model_opt |= OPT_L;
    }

    A = gretl_model_get_data(pmod, "midas_info");

    if (A != NULL) {
        b = gretl_array_get_bundle(A, 0);
    }

    if (b != NULL) {
        strcpy(mds_listname, gretl_bundle_get_string(b, "lname", &err));
        mds_quad[0] = gretl_bundle_get_int(b, "minlag", &err);
        mds_quad[1] = gretl_bundle_get_int(b, "maxlag", &err);
        mds_quad[2] = gretl_bundle_get_int(b, "type", &err);
        mds_quad[3] = gretl_bundle_get_int(b, "nparm", &err);
    }
}

/* support for the "Modify model..." Edit menu item */

void selector_from_model (windata_t *vwin)
{
    void *ptr = vwin->data;
    int ci = vwin->role;

    model_opt = OPT_NONE;

    if (ci == VIEW_MODEL) {
        /* single-equation model (mostly) */
        MODEL *pmod = (MODEL *) ptr;
	const char *cname = NULL;
        int sel_ci = pmod->ci;
        int dv = -1, gotinst = 0;

        if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM) {
            revise_nl_model(pmod, vwin_toplevel(vwin));
            return;
        }

        if (pmod->ci == INTREG) {
            lovar = varnum_from_keystring(pmod, "lovar");
            hivar = varnum_from_keystring(pmod, "hivar");
        } else {
            dv = gretl_model_get_depvar(pmod);
            if (dv >= 0 && dv < dataset->v) {
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
            garch_const = pmod->ifc;
            if (pmod->opt & OPT_F) {
                model_opt |= OPT_F;
            }
        } else if (pmod->ci == HECKIT) {
            retrieve_heckit_info(pmod, &gotinst);
        } else if (COUNT_MODEL(pmod->ci)) {
            if (pmod->ci == NEGBIN) {
                if (pmod->opt & OPT_M) {
                    model_opt |= OPT_M;
                } else {
                    model_opt |= OPT_N;
                }
            }
            sel_ci = COUNTMOD;
            offvar = gretl_model_get_int(pmod, "offset_var");
        } else if (pmod->ci == DURATION) {
            if (pmod->opt & OPT_E) {
                model_opt |= OPT_E;
            } else if (pmod->opt & OPT_L) {
                model_opt |= OPT_L;
            } else if (pmod->opt & OPT_Z) {
                model_opt |= OPT_Z;
            }
            censvar = gretl_model_get_int(pmod, "cens_var");
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
                model_opt |= OPT_P;
            } else if (pmod->opt & OPT_H) {
                model_opt |= OPT_H;
            }
        } else if (pmod->ci == PANEL) {
            if (pmod->opt & OPT_F) {
                model_opt |= OPT_F;
            } else if (pmod->opt & OPT_U) {
                model_opt |= OPT_U;
            } else if (pmod->opt & OPT_H) {
                sel_ci = PANEL_WLS;
            } else if (pmod->opt & OPT_B) {
                sel_ci = PANEL_B;
            }
            if (pmod->opt & OPT_D) {
                model_opt |= OPT_D;
            }
        } else if (pmod->ci == DPANEL) {
            if (pmod->opt & OPT_A) {
                model_opt |= OPT_A;
            }
            if (pmod->opt & OPT_D) {
                model_opt |= OPT_D;
            }
            if (pmod->opt & OPT_T) {
                model_opt |= OPT_T;
            }
            if (pmod->opt & OPT_L) {
                model_opt |= OPT_L;
            }
            if (pmod->opt & OPT_X) {
                model_opt |= OPT_X;
            }
            if (pmod->opt & OPT_C) {
                model_opt |= OPT_C;
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
            } else if (pmod->opt & OPT_E) {
                sel_ci = REPROBIT;
            }
        } else if (pmod->ci == TOBIT) {
            retrieve_tobit_info(pmod);
        } else if (pmod->ci == BIPROBIT) {
            retrieve_biprobit_info(pmod, &gotinst);
        } else if (pmod->ci == MIDASREG) {
            retrieve_midas_info(pmod);
        } else if (pmod->ci == LOGISTIC) {
            if (pmod->opt & OPT_F) {
                sel_ci = FE_LOGISTIC;
            }
        }
        if (pmod->opt & OPT_R) {
            model_opt |= OPT_R;
        }

        *cluster_var = '\0';
        cname = gretl_model_get_cluster_vname(pmod);
        if (cname != NULL) {
	    strcpy(cluster_var, cname);
            model_opt |= OPT_C;
        }

        y_x_lags_enabled = y_w_lags_enabled = 0;

        if ((dataset_is_time_series(dataset) ||
             dataset_is_panel(dataset)) &&
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
        vecxlist = NULL;

        if (var->rlist != NULL) {
            vecxlist = gretl_lists_join_with_separator(var->xlist,
                                                       var->rlist);
        } else if (var->xlist != NULL) {
            vecxlist = gretl_list_copy(var->xlist);
        }

        set_lag_prefs_from_VAR(var->lags, vecxlist);
        default_order = var->order;

        if (ci == VECM) {
            jrank = gretl_VECM_rank(var);
            jcase = jcode(var);
            default_order += 1;
        }

        selection_dialog(ci, (ci == VAR)? _("gretl: VAR") : _("gretl: VECM"),
                         NULL, do_vector_model);
    } else if (ci == SYSTEM) {
        revise_system_model(ptr, vwin_toplevel(vwin));
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
            ynum = series_index(dataset, s);
            if (ynum == dataset->v) {
                ynum = -1;
            }
        }
    }

    return ynum;
}

static int depvar_selected (const selector *sr)
{
    return selector_get_depvar_number(sr) > 0;
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
            gtk_tree_model_get(mod, &piter, COL_ID, &pv, COL_LAG, &plag, -1);
            if (pv == v && plag != lag) {
                gtk_tree_model_get(mod, &piter, COL_FLAG, &flag, -1);
                gtk_list_store_set(GTK_LIST_STORE(mod), iter, COL_FLAG, flag, -1);
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
        gtk_list_store_set(GTK_LIST_STORE(mod), iter, COL_FLAG, varflag, -1);
    }
}

static void
real_varlist_set_var (int v, int lag, GtkTreeModel *mod, GtkTreeIter *iter)
{
    int ncols = gtk_tree_model_get_n_columns(mod);
    GtkListStore *store = GTK_LIST_STORE(mod);

    if (lag == 0) {
        gtk_list_store_set(store, iter, COL_ID, v, COL_LAG, 0,
                           COL_NAME, dataset->varname[v],
                           -1);
    } else {
        char vstr[VNAMELEN + 8];

        sprintf(vstr, "%s(-%d)", dataset->varname[v], lag);
        gtk_list_store_set(store, iter, COL_ID, v, COL_LAG, lag,
                           COL_NAME, vstr, -1);
    }

    if (ncols == 4) {
        if (lag == 0) {
            gtk_list_store_set(store, iter, COL_FLAG, varflag, -1);
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
            dataset->varname[v]);
#endif

    if (gtk_tree_model_get_iter_first(mod, iter)) {
        last = iter;
        while (1) {
            gtk_tree_model_get(mod, iter, COL_ID, &tv, -1);
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

    if (v > 0 && dataset_lags_ok(dataset)) {
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
        gtk_tree_model_get(model, &iter, COL_ID, &varnum, -1);
        if (sr != NULL) {
            sr->active_var = varnum;
        }
        row = tree_path_get_row_number(path);
        gtk_tree_path_free(path);
        /* note: the following is used in listbox_drag() */
        g_object_set_data(G_OBJECT(widget), "active_row",
                          GINT_TO_POINTER(row));
    }

    return FALSE;
}

static void list_append_var_simple (GtkListStore *store,
                                    GtkTreeIter *iterp,
                                    int v)
{
    const char *vname = dataset->varname[v];

    gtk_list_store_append(store, iterp);
    gtk_list_store_set(store, iterp,
                       COL_ID, v,
                       COL_LAG, 0,
                       COL_NAME, vname,
                       -1);
}

static void list_append_var (GtkTreeModel *mod,
                             GtkTreeIter *iter,
                             int v, selector *sr,
                             int locus)
{
    int i, lcontext = 0;

    if (v > 0 && dataset_lags_ok(dataset)) {
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

static void list_append_midas_var (GtkListStore *store,
                                   GtkTreeIter *iterp,
                                   int v, int m)
{
    char mname[VNAMELEN];
    const char *vname;
    int ID = v;

    vname = get_listname_by_consecutive_content(m, v);

    if (vname != NULL) {
        /* got a pre-existing MIDAS list */
        strcpy(mname, vname);
        ID = -1;
    } else {
        /* got the "anchor" of a potential list */
        char *p;

        vname = dataset->varname[v];
        strcpy(mname, vname);
        p = strrchr(mname, '_');
        if (p != NULL) {
            *p = '\0';
        }
    }

    gtk_list_store_append(store, iterp);
    gtk_list_store_set(store, iterp,
                       COL_ID, ID, MCOL_M, m,
                       MCOL_NAME, mname, -1);

    if (!strcmp(mname, mds_listname)) {
        /* FIXME */
        fprintf(stderr, "Should add %s on right?\n", mname);
    }
}

static void render_varname (GtkTreeViewColumn *column,
                            GtkCellRenderer *renderer,
                            GtkTreeModel *model,
                            GtkTreeIter *iter,
                            gpointer p)
{
    gint id;

    gtk_tree_model_get(model, iter, COL_ID, &id, -1);
    if (id < 0) {
        g_object_set(renderer, "weight", PANGO_WEIGHT_BOLD, NULL);
    } else {
        g_object_set(renderer, "weight", PANGO_WEIGHT_NORMAL, NULL);
    }
}

/* build a new liststore and associated tree view, and pack into the
   given @hbox */

static GtkWidget *var_list_box_new (GtkBox *hbox, selector *sr, int locus)
{
    GtkListStore *store;
    GtkWidget *view, *scroller;
    GtkCellRenderer *renderer;
    GtkTreeViewColumn *column;
    GtkTreeSelection *select;
    gboolean flagcol = FALSE;
    gboolean midascol = FALSE;
    int cw, width = 160;
    int height = -1;

    cw = get_char_width(sr->dlg);
    if (cw > 0) {
        width = 18 * cw;
    }

    if (USE_RXLIST(sr->ci) && locus == SR_RVARS2) {
        /* VECM special, with restricted/unrestricted flag column */
        store = gtk_list_store_new(4, G_TYPE_INT, G_TYPE_INT,
                                   G_TYPE_STRING, G_TYPE_STRING);
        flagcol = TRUE;
    } else if (sr->ci == MIDASREG && locus == SR_RVARS2) {
        /* MIDAS special with parameterization info */
        store = gtk_list_store_new(7, G_TYPE_INT, G_TYPE_INT,
                                   G_TYPE_STRING, G_TYPE_INT,
                                   G_TYPE_INT, G_TYPE_INT,
                                   G_TYPE_INT);
        midascol = TRUE;
    } else {
        /* ID number, lag or frequency ratio, varname */
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
                                                      "text",
                                                      COL_NAME,
                                                      NULL);
    gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);

    if (locus == SR_RVARS2 && sr->ci == MIDASREG) {
        gretl_tooltips_add(view, _("Select and right-click to edit specification"));
    }

    if (flagcol) {
        column = gtk_tree_view_column_new_with_attributes(NULL,
                                                          renderer,
                                                          "text",
                                                          COL_FLAG,
                                                          NULL);
        gtk_tree_view_append_column(GTK_TREE_VIEW(view), column);
    } else {
        gtk_tree_view_column_set_cell_data_func(column, renderer,
                                                render_varname,
                                                NULL, NULL);
    }

    gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(view), FALSE);
    gtk_tree_view_set_reorderable(GTK_TREE_VIEW(view), FALSE);

    select = gtk_tree_view_get_selection(GTK_TREE_VIEW(view));

    if (locus == SR_LVARS2 || (locus == SR_RVARS2 && sr->ci == MIDASREG)) {
        gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
    } else {
        gtk_tree_selection_set_mode(select, GTK_SELECTION_MULTIPLE);
        g_signal_connect(G_OBJECT(view), "motion-notify-event",
                         G_CALLBACK(listbox_drag), NULL);
    }

    /* enable interactive search on name */
    gtk_tree_view_set_search_column(GTK_TREE_VIEW(view), COL_NAME);
    gtk_tree_view_set_enable_search(GTK_TREE_VIEW(view), TRUE);

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
        } else if (midascol) {
            g_signal_connect(G_OBJECT(view), "button-press-event",
                             G_CALLBACK(listvar_midas_click),
                             sr);
        } else {
            g_signal_connect(G_OBJECT(view), "button-press-event",
                             G_CALLBACK(listvar_reorder_click),
                             view);
        }
    } else if (locus == SR_LVARS2) {
        g_signal_connect(G_OBJECT(view), "button-press-event",
                         G_CALLBACK(lvars_right_click),
                         sr);
    }

    scroller = gtk_scrolled_window_new(NULL, NULL);
    gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                   GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
    gtk_scrolled_window_set_shadow_type(GTK_SCROLLED_WINDOW(scroller),
                                        GTK_SHADOW_IN);
    gtk_container_add(GTK_CONTAINER(scroller), view);

    gtk_box_pack_start(hbox, scroller, TRUE, TRUE, 0);

    if (gui_scale > 1.0) {
        width *= 1.0 + 0.4 * (gui_scale - 1);
    }
    gtk_widget_set_size_request(view, width, height);
    gtk_widget_show(view);

#if GTK_MAJOR_VERSION >= 3
    gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(scroller),
                                              width);
#endif

    gtk_widget_show(scroller);

    return view;
}

static int binary_var_check (int v, const char *vname)
{
    if (!gretl_isdummy(dataset->t1, dataset->t2, dataset->Z[v])) {
        errbox_printf(_("The variable '%s' is not a 0/1 variable."), vname);
        return 1;
    }

    return 0;
}

/* add to "extra" var slot the current selection from sr->lvars */

static void real_set_extra_var (GtkTreeModel *model, GtkTreePath *path,
                                GtkTreeIter *iter, selector *sr)
{
    gchar *vname;
    gint v;

    gtk_tree_model_get(model, iter, COL_ID, &v, COL_NAME, &vname, -1);

    if (v < 0) {
        return;
    }

    if (sr->ci == HECKIT || sr->ci == BIPROBIT) {
        if (binary_var_check(v, vname)) {
            return;
        }
    } else if (NONPARAM_CODE(sr->ci)) {
        const gchar *test = gtk_entry_get_text(GTK_ENTRY(sr->depvar));

        if (!strcmp(vname, test)) {
            /* can't have the same var in both places */
            gtk_entry_set_text(GTK_ENTRY(sr->depvar), "");
        }
    }

    gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->extra[0]), "data",
                      GINT_TO_POINTER(v));
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
    gchar *vname = NULL;
    gint v;

    gtk_tree_model_get(model, iter, COL_ID, &v, COL_NAME, &vname, -1);

    if (v < 0) {
        g_free(vname);
        return;
    }

    gtk_entry_set_text(GTK_ENTRY(sr->rvars1), vname);
    g_free(vname);
    g_object_set_data(G_OBJECT(sr->rvars1), "data",
                      GINT_TO_POINTER(v));
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

/* When adding a variable to the Endogenous or Exogenous listing
   for VAR, VECM, etc., check that the variable in question is not
   present in the other listing: if it is, then remove it.

   And similarly, when constructing a list of function interfaces
   for a package, don't allow selection of a given function as
   both a public interface and a private ("Helper") function.
*/

static void dual_selection_fix_conflicts (selector *sr, int new_locus)
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

    /* construct list of vars in @src */

    if (gtk_tree_model_get_iter_first(src, &iter)) {
        do {
            gtk_tree_model_get(src, &iter, COL_ID, &v, -1);
            gretl_list_append_term(&srclist, v);
        } while (gtk_tree_model_iter_next(src, &iter));
    }

    /* check @srclist against @targ */

    if (srclist != NULL && gtk_tree_model_get_iter_first(targ, &iter)) {
        do {
            ok = TRUE;
            gtk_tree_model_get(targ, &iter, COL_ID, &v, -1);
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
    int append = 1;
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
                gtk_tree_model_get(mod, &iter, COL_ID, &modv, -1);
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
            gtk_tree_model_get(model, &iter, COL_ID, &xnum, -1);
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

static void np_xvar_cleanup (selector *sr, int newy)
{
    const gchar *test = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));

    if (test != NULL && *test != '\0') {
        int xv = current_series_index(dataset, test);

        if (xv == newy) {
            gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), "");
        }
    }
}

static int set_dependent_var_from_active (selector *sr)
{
    gint v = sr->active_var;
    const char *vname;

    if (v < 0 || v >= dataset->v || sr->depvar == NULL) {
        return 1;
    }

    vname = dataset->varname[v];

    if (sr->ci == PROBIT || sr->ci == BIPROBIT) {
        if (binary_var_check(v, vname)) {
            return 1;
        }
    }

    /* models: if we select foo as regressand, remove it from the list
       of regressors if need be; also remove lags associated with the
       previous dependent var, if any.
    */
    if (MODEL_CODE(sr->ci)) {
        dependent_var_cleanup(sr, v);
    } else if (NONPARAM_CODE(sr->ci)) {
        np_xvar_cleanup(sr, v);
    }

    gtk_entry_set_text(GTK_ENTRY(sr->depvar), vname);

    if (select_lags_depvar(sr->ci)) {
        maybe_insert_depvar_lags(sr, v, 0);
    }

    return 0;
}

static void real_set_dependent_var (GtkTreeModel *model, GtkTreePath *path,
                                    GtkTreeIter *iter, selector *sr)
{
    gchar *vname;
    gint v;

    gtk_tree_model_get(model, iter, COL_ID, &v, COL_NAME, &vname, -1);

    if (v < 0) {
        return;
    }

    if (sr->ci == PROBIT || sr->ci == BIPROBIT) {
        if (binary_var_check(v, vname)) {
            return;
        }
    }

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

    selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
    gtk_tree_selection_selected_foreach(selection,
                                        (GtkTreeSelectionForeachFunc)
                                        real_set_dependent_var,
                                        sr);
}

static void set_right_var_from_main (GtkTreeModel *model, GtkTreePath *path,
                                     GtkTreeIter *iter, selector *sr)
{
    GtkTreeModel *rmod;
    GtkTreeIter r_iter;
    gchar *idstr = NULL;
    int v;

    gtk_tree_model_get(model, iter, COL_ID, &idstr, -1);

    rmod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (rmod == NULL) {
        g_free(idstr);
        return;
    }

    v = atoi(idstr);

    if (gtk_tree_model_get_iter_first(rmod, &r_iter)) {
        while (gtk_tree_model_iter_next(rmod, &r_iter)) {
            ;
        }
    }

    gtk_list_store_append(GTK_LIST_STORE(rmod), &r_iter);
    real_varlist_set_var(v, 0, rmod, &r_iter);

    g_free(idstr);
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
    gtk_tree_model_get(lmod, &iter, COL_ID, &v, -1);

    gtk_tree_model_get_iter_first(rmod, &iter);
    gtk_list_store_append(GTK_LIST_STORE(rmod), &iter);
    gtk_list_store_set(GTK_LIST_STORE(rmod), &iter,
                       COL_ID, v, COL_LAG, 0,
                       COL_NAME, dataset->varname[v], -1);
}

static int varflag_dialog (int v, GtkWidget *parent)
{
    const char *opts[] = {
        N_("Unrestricted"),
        N_("Restricted")
    };
    gchar *title, *label;
    int ret;

    title = g_strdup_printf("gretl: %s", _("add exogenous variable"));
    label = g_strdup_printf(_("Status of '%s' in VECM:"), dataset->varname[v]);

    ret = radio_dialog(title, label, opts, 2, 0, 0, parent);
    if (ret == 0) {
        set_varflag(UNRESTRICTED);
    } else if (ret == 1) {
        set_varflag(RESTRICTED);
    }

    g_free(title);
    g_free(label);

    return ret;
}

static int arima_selected (selector *sr)
{
    int ret = 0;

    if (sr->extra[ARIMA_d] != NULL) {
        /* the arima_d spin button */
        ret = spin_get_int(sr->extra[ARIMA_d]);
    }

    if (!ret && sr->extra[ARIMA_D] != NULL) {
        /* the seasonal arima_D spin button */
        ret = spin_get_int(sr->extra[ARIMA_D]);
    }

    return ret;
}

static void xdiff_button_set_sensitive (selector *sr, gboolean s)
{
    if (sr->xdiff_button != NULL) {
        s = s && arima_selected(sr);
        gtk_widget_set_sensitive(sr->xdiff_button, s);
    }
}

static int rvars1_n_vars (selector *sr)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    int nv = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (model == NULL) {
        return 0;
    }

    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            nv++;
        } while (gtk_tree_model_iter_next(model, &iter));
    }

    return nv;
}

/* add a variable (or possibly a list of variables) to the listbox
   at @locus */

static void real_add_generic (GtkTreeModel *srcmodel,
                              GtkTreeIter *srciter,
                              selector *sr,
                              int locus)
{
    GtkWidget *w;
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *vname = NULL;
    const int *vlist = NULL;
    gint v, xnum;
    gint at_max = 0;
    gint keep_names = 0;
    int i, addvars = 1;
    int nvars = 0;
    int err = 0;

    w = (locus == SR_RVARS2)? sr->rvars2 : sr->rvars1;

    if (!GTK_IS_TREE_VIEW(w)) {
        return;
    }

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(w));
    if (model == NULL) {
        return;
    }

    /* get the 'source' info */
    gtk_tree_model_get(srcmodel, srciter, COL_ID, &v, COL_NAME, &vname, -1);

    if (v < 0) {
        /* we should have a gretl list, not a single variable */
        vlist = get_list_by_name(vname);
        if (vlist == NULL) {
            err = E_DATA;
        } else {
            addvars = vlist[0];
        }
    } else {
        keep_names =
            GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->lvars),
                                              "keep-names"));
    }

    if (err) {
        gui_errmsg(err);
        return;
    }

    if (!keep_names) {
        g_free(vname);
        vname = NULL;
    }

    for (i=0; i<addvars && !at_max; i++) {
        int already_there = 0;

        if (vlist != NULL) {
            v = vlist[i+1];
        }

        /* first check if we're maxed out, or if the variable to
           add is already present */

        if (gtk_tree_model_get_iter_first(model, &iter)) {
            do {
                if (i == 0 && selection_at_max(sr, w, ++nvars)) {
                    at_max = 1;
                }
                if (!at_max && !already_there) {
                    gtk_tree_model_get(model, &iter, COL_ID, &xnum, -1);
                    if (xnum == v) {
                        already_there = 1;
                    }
                }
            } while (gtk_tree_model_iter_next(model, &iter));
        }

        if (!already_there && !at_max) {
            /* OK to append */
            if (keep_names) {
                int ncols = gtk_tree_model_get_n_columns(model);

                gtk_list_store_append(GTK_LIST_STORE(model), &iter);
                gtk_list_store_set(GTK_LIST_STORE(model), &iter,
                                   COL_ID, v, COL_LAG, 0, COL_NAME, vname, -1);
                if (ncols == 4) {
                    gtk_list_store_set(GTK_LIST_STORE(model), &iter,
                                       COL_FLAG, varflag, -1);
                }
                g_free(vname);
            } else {
                if (locus == SR_RVARS2 && USE_RXLIST(sr->ci)) {
                    if (varflag_dialog(v, sr->dlg) < 0) {
                        return;
                    }
                }
#if VLDEBUG
                fprintf(stderr, "real_add_generic: calling varlist_insert_var_full\n");
#endif
                varlist_insert_var_full(v, model, &iter, sr, locus);
            }

            if (selection_at_max(sr, w, ++nvars)) {
                at_max = 1;
            }
        }
    }

    if (sr->add_button != NULL && at_max) {
        gtk_widget_set_sensitive(sr->add_button, FALSE);
    }

    if (locus == SR_RVARS1 && sr->ci == CORR) {
        set_n_rvars1(sr, nvars);
    }

    if (nvars > 0) {
        if (lags_button_relevant(sr, locus)) {
            gtk_widget_set_sensitive(sr->lags_button, TRUE);
        } else if (VECLAGS_CODE(sr->ci) && locus == SR_RVARS1 &&
                   sr->extra[EXTRA_LAGS] != NULL) {
            gtk_widget_set_sensitive(sr->extra[EXTRA_LAGS], TRUE);
        }
        if (ARMA_RELATED(sr->ci)) {
            xdiff_button_set_sensitive(sr, TRUE);
        }
        if (USE_VECXLIST(sr->ci) || FNPKG_CODE(sr->ci)) {
            dual_selection_fix_conflicts(sr, locus);
        }
    }
}

/* shift a MIDAS term from one location in the selection
   dialog to another */

static void move_midas_term (GtkTreeModel *src,
                             GtkTreeIter *srciter,
                             selector *sr)
{
    GtkTreeModel *targ;
    GtkTreeIter iter;
    gchar *vname = NULL;
    int v, m, src_cols;

    /* get the 'source' info */
    gtk_tree_model_get(src, srciter,
                       COL_ID, &v, MCOL_M, &m,
                       MCOL_NAME, &vname, -1);

    /* append to target */
    src_cols = gtk_tree_model_get_n_columns(src);
    if (src_cols == 3) {
        /* going left to right (selecting): defaults */
        gboolean have_beta1 = 0;
        gboolean no_beta1 = 0;
        int lmin = 1, lmax = 2*m;
        int ptype = mds_quad[2];
        int k = mds_quad[3];
        int nterms;

        if (mds_quad[0] < mds_quad[1]) {
            lmin = mds_quad[0];
            lmax = mds_quad[1];
        }

        targ = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));

        /* If we have a beta1 term in place on the right,
           don't allow adding anything else. And if there's
           anything in place already, don't allow adding a
           beta1 term.
        */
        nterms = check_midas_rvars2(targ, &have_beta1);
        if (have_beta1) {
            warnbox("A one-parameter beta term cannot be combined with others");
            return;
        } else if (nterms > 0) {
            no_beta1 = 1;
        }

        if (midas_term_dialog(vname, m, &lmin, &lmax, &ptype,
                              &k, no_beta1, sr->dlg) < 0) {
            return;
        }
        gtk_list_store_append(GTK_LIST_STORE(targ), &iter);
        gtk_list_store_set(GTK_LIST_STORE(targ), &iter,
                           COL_ID, v, MCOL_M, m, MCOL_NAME, vname,
                           MCOL_MINLAG, lmin, MCOL_MAXLAG, lmax,
                           MCOL_TYPE, ptype, MCOL_K, k,
                           -1);
    } else {
        /* going right to left (deselecting) */
        targ = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars2));
        gtk_list_store_append(GTK_LIST_STORE(targ), &iter);
        gtk_list_store_set(GTK_LIST_STORE(targ), &iter,
                           COL_ID, v, MCOL_M, m,
                           MCOL_NAME, vname, -1);
    }

    g_free(vname);

    /* and remove from source */
    gtk_list_store_remove(GTK_LIST_STORE(src), srciter);
}

static void add_to_rvars1 (GtkTreeModel *model, GtkTreePath *path,
                           GtkTreeIter *iter, selector *sr)
{
    if (MODEL_CODE(sr->ci)) {
        /* don't add the regressand to the list of regressors */
        gint xnum, ynum;

        gtk_tree_model_get(model, iter, COL_ID, &xnum, -1);
        ynum = selector_get_depvar_number(sr);
        if (ynum >= 0 && xnum == ynum) {
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
    GtkTreeSelection *sel;

    if (sr->ci == MIDASREG) {
        GtkTreeModel *model = NULL;
        GtkTreeIter iter;

        sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars2));
        if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
            move_midas_term(model, &iter, sr);
        }
    } else {
        sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
        gtk_tree_selection_selected_foreach(sel,
                                            (GtkTreeSelectionForeachFunc)
                                            add_to_rvars2,
                                            sr);
    }
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

static void remove_from_right (GtkWidget *w, selector *sr,
                               GtkWidget *listbox)
{
    GtkTreeView *view = GTK_TREE_VIEW(listbox);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeSelection *selection = gtk_tree_view_get_selection(view);
    GtkTreePath *path;
    GtkTreeIter iter, last;
    int context = 0;
    int v, lag;
    int *sellist = NULL;
    int ridx, nrows = 0;

    if (model == NULL || selection == NULL) {
        return;
    }

    /* determine the number of rows in the list box, create a
       list of selected row indices, and navigate to the last row
    */
    if (gtk_tree_model_get_iter_first(model, &iter)) {
        if (gtk_tree_selection_iter_is_selected(selection, &iter)) {
            sellist = gretl_list_append_term(&sellist, nrows);
        }
        last = iter;
        nrows++;
        while (gtk_tree_model_iter_next(model, &iter)) {
            if (gtk_tree_selection_iter_is_selected(selection, &iter)) {
                sellist = gretl_list_append_term(&sellist, nrows);
            }
            last = iter;
            nrows++;
        }
    }

    if (nrows == 0 || sellist == NULL) {
        /* "can't happen", but... */
        return;
    }

    context = lag_context_from_widget(GTK_WIDGET(view));

    /* work back upward, deleting the selected rows */
    path = gtk_tree_model_get_path(model, &last);

    while (1) {
        ridx = gtk_tree_path_get_indices(path)[0];
        if (gtk_tree_model_get_iter(model, &last, path) &&
            in_gretl_list(sellist, ridx)) {
            if (context) {
                gtk_tree_model_get(model, &last, COL_ID, &v, COL_LAG, &lag, -1);
                remove_specific_lag(v, lag, context);
            }
            gtk_list_store_remove(GTK_LIST_STORE(model), &last);
            nrows--;
        }
        if (!gtk_tree_path_prev(path)) {
            break;
        }
    }

    if (sr->add_button != NULL &&
        !gtk_widget_is_sensitive(sr->add_button) &&
        !selection_at_max(sr, w, nrows)) {
        gtk_widget_set_sensitive(sr->add_button, TRUE);
    }

    if (sr->ci == CORR) {
        set_n_rvars1(sr, nrows);
    }

    if (nrows == 0) {
        /* the listbox is now empty */
        if (context && sr->lags_button != NULL) {
            /* can't do independent var lags, but can we still
               do dependent var lags?
            */
            int y_lags_ok = select_lags_depvar(sr->ci) && depvar_selected(sr);

            gtk_widget_set_sensitive(sr->lags_button, y_lags_ok);
        }
        if (GTK_WIDGET(view) == sr->rvars1 && ARMA_RELATED(sr->ci)) {
            xdiff_button_set_sensitive(sr, FALSE);
        }
    } else if (sr->ci == MIDASREG && context && nrows == 1) {
        /* If only the constant remains, lag selection
           for X-vars should be turned off
        */
        gtk_tree_model_get_iter_first(model, &iter);
        gtk_tree_model_get(model, &iter, COL_ID, &v, -1);
        if (v == 0) {
            gtk_widget_set_sensitive(sr->lags_button, FALSE);
        }
    }

    free(sellist);
}

static void remove_from_rvars1_callback (GtkWidget *w, selector *sr)
{
    remove_from_right(w, sr, sr->rvars1);
}

static void remove_from_rvars2_callback (GtkWidget *w, selector *sr)
{
    if (sr->ci == MIDASREG) {
        GtkTreeSelection *sel;
        GtkTreeModel *model = NULL;
        GtkTreeIter iter;

        sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->rvars2));
        if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
            move_midas_term(model, &iter, sr);
        }
    } else {
        remove_from_right(w, sr, sr->rvars2);
    }
}

/* double-clicking sets the dependent variable and marks it as the
   default, if applicable */

static gint
dblclick_lvars_row (GtkWidget *w, GdkEventButton *event, selector *sr)
{
    if (sr->depvar != NULL && event != NULL &&
        event->type == GDK_2BUTTON_PRESS) {
        int err = set_dependent_var_from_active(sr);

        if (!err && sr->default_check != NULL) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
                                         TRUE);
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

    gtk_tree_model_get(model, iter, COL_FLAG, &orig, -1);

    if (orig != NULL) {
        repl = (i == 0)? "U" : "R";
        if (strcmp(orig, repl)) {
            gtk_list_store_set(GTK_LIST_STORE(model), iter, COL_FLAG, repl, -1);
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
    g_signal_connect_swapped(G_OBJECT(item), "destroy",
			     G_CALLBACK(gtk_widget_destroy),
			     popup);
    gtk_widget_show(item);
    gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);
}

static gint listvar_flagcol_click (GtkWidget *widget,
                                   GdkEventButton *event,
                                   gpointer data)
{
    GtkWidget *view = GTK_WIDGET(data);
    GtkWidget *flag_popup;
    int i;

    if (right_click(event)) {
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

/* Below: we're allowing the user to revise the MIDAS
   parameterization of a previously selected term
*/

static gint listvar_midas_click (GtkWidget *widget,
                                 GdkEventButton *event,
                                 selector *sr)
{
    if (right_click(event)) {
        GtkTreeSelection *sel;
        GtkTreeModel *model = NULL;
        GtkTreeIter iter;

        sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(widget));
        if (gtk_tree_selection_get_selected(sel, &model, &iter)) {
            int m, lmin, lmax, ptype, k, resp;
            gchar *vname = NULL;
            gboolean no_beta1;

            no_beta1 = (check_midas_rvars2(model, NULL) > 1);
            gtk_tree_model_get(model, &iter, MCOL_M, &m, MCOL_NAME, &vname,
                               MCOL_MINLAG, &lmin, MCOL_MAXLAG, &lmax,
                               MCOL_TYPE, &ptype, MCOL_K, &k, -1);
            resp = midas_term_dialog(vname, m, &lmin, &lmax, &ptype,
                                     &k, no_beta1, sr->dlg);
            if (resp != GRETL_CANCEL) {
                gtk_list_store_set(GTK_LIST_STORE(model), &iter,
                                   MCOL_MINLAG, lmin, MCOL_MAXLAG, lmax,
                                   MCOL_TYPE, ptype, MCOL_K, k,
                                   -1);
            }
        }

        return TRUE;
    }

    return FALSE;
}

static int selection_get_row_info (GtkTreeModel *model,
                                   GtkTreeSelection *sel,
                                   int *r0, int *r1, int *rmin)
{
    GtkTreeIter iter;
    int i, rmax = 0;

    *r0 = 10000;
    *r1 = -1;
    *rmin = 0;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
        gint id;

        gtk_tree_model_get(model, &iter, COL_ID, &id, -1);
        if (id == 0) {
            /* don't shift const from position 0 */
            *rmin = 1;
        }
        for (i=0; ; i++) {
            if (gtk_tree_selection_iter_is_selected(sel, &iter)) {
                if (i < *r0) {
                    *r0 = i;
                }
                if (i > *r1) {
                    *r1 = i;
                }
            }
            if (!gtk_tree_model_iter_next(model, &iter)) {
                break;
            }
        }
        rmax = i;
    }

    return rmax;
}

static void set_selection_from_list (GtkTreeView *view,
                                     GtkTreeModel *model,
                                     GtkTreeSelection *sel,
                                     const int *list)
{
    GtkTreeIter iter;
    int id, nsel = 0;

    gtk_tree_selection_unselect_all(sel);
    gtk_tree_model_get_iter_first(model, &iter);

    while (gtk_tree_model_iter_next(model, &iter) && nsel < list[0]) {
        gtk_tree_model_get(model, &iter, COL_ID, &id, -1);
        if (in_gretl_list(list, id)) {
            gtk_tree_selection_select_iter(sel, &iter);
            nsel++;
        }
    }
}

static int swap_row_content (GtkTreeModel *model, GtkTreeIter *iter0,
                             GtkTreeIter *iter1, int flags)
{
    GtkListStore *store = GTK_LIST_STORE(model);
    gint id0, id1, l0, l1;
    gchar *s0, *s1;

    gtk_tree_model_get(model, iter0, COL_ID, &id0, COL_LAG, &l0,
                       COL_NAME, &s0, -1);
    gtk_tree_model_get(model, iter1, COL_ID, &id1, COL_LAG, &l1,
                       COL_NAME, &s1, -1);
    gtk_list_store_set(store, iter0, COL_ID, id1, COL_LAG, l1,
                       COL_NAME, s1, -1);
    gtk_list_store_set(store, iter1, COL_ID, id0, COL_LAG, l0,
                       COL_NAME, s0, -1);

    g_free(s0);
    g_free(s1);

    if (flags) {
        gtk_tree_model_get(model, iter0, COL_FLAG, &s0);
        gtk_tree_model_get(model, iter1, COL_FLAG, &s1);
        gtk_list_store_set(store, iter0, COL_FLAG, s1);
        gtk_list_store_set(store, iter1, COL_FLAG, s0);

        g_free(s0);
        g_free(s1);
    }

    return id1;
}

static void move_selected_rows (GtkMenuItem *item, GtkTreeView *view)
{
    gint down = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item), "down"));
    GtkTreeSelection *sel = gtk_tree_view_get_selection(view);
    GtkTreeModel *model = gtk_tree_view_get_model(view);
    GtkTreeIter iter0, iter1;
    GtkTreePath *path;
    gboolean flags = FALSE;
    int id, *sellist = NULL;

    if (gtk_tree_view_get_column(view, 3) != NULL) {
        flags = TRUE;
    }

    if (down) {
        tree_model_get_iter_last(model, &iter1);
        while (tree_model_iter_prev(model, &iter1)) {
            if (gtk_tree_selection_iter_is_selected(sel, &iter1)) {
                path = gtk_tree_model_get_path(model, &iter1);
                gtk_tree_path_next(path);
                gtk_tree_model_get_iter(model, &iter0, path);
                id = swap_row_content(model, &iter0, &iter1, flags);
                gretl_list_append_term(&sellist, id);
                gtk_tree_path_free(path);
            }
        }
    } else {
        gtk_tree_model_get_iter_first(model, &iter1);
        while (gtk_tree_model_iter_next(model, &iter1)) {
            if (gtk_tree_selection_iter_is_selected(sel, &iter1)) {
                path = gtk_tree_model_get_path(model, &iter1);
                gtk_tree_path_prev(path);
                gtk_tree_model_get_iter(model, &iter0, path);
                id = swap_row_content(model, &iter0, &iter1, flags);
                gretl_list_append_term(&sellist, id);
                gtk_tree_path_free(path);
            }
        }
    }

    if (sellist != NULL) {
        set_selection_from_list(view, model, sel, sellist);
        free(sellist);
    }

    gtk_widget_destroy(GTK_WIDGET(item));
}

static gint listvar_reorder_click (GtkWidget *widget, GdkEventButton *event,
                                   gpointer data)
{
    if (right_click(event)) {
        GtkTreeView *view = GTK_TREE_VIEW(data);
        GtkTreeModel *model = gtk_tree_view_get_model(view);
        GtkTreeSelection *sel = gtk_tree_view_get_selection(view);
        int r0, r1, rmin, rmax;

        rmax = selection_get_row_info(model, sel, &r0, &r1, &rmin);

        if (r0 >= 0 && !(rmin == 1 && r0 == 0)) {
            const gchar *icons[] = {
                GTK_STOCK_GO_UP,
                GTK_STOCK_GO_DOWN
            };
            GtkWidget *popup = gtk_menu_new();
            GtkWidget *item;
            int i;

            for (i=0; i<2; i++) {
                if (i == 0 && r0 <= rmin) {
                    /* can't do move up */
                    continue;
                } else if (i == 1 && r1 >= rmax) {
                    /* can't do move down */
                    continue;
                }
                item = gtk_image_menu_item_new_from_stock(icons[i], NULL);
                g_object_set_data(G_OBJECT(item), "down",
                                  GINT_TO_POINTER(i));
                g_signal_connect(G_OBJECT(item), "activate",
                                 G_CALLBACK(move_selected_rows), view);
                g_signal_connect_swapped(G_OBJECT(item), "destroy",
					 G_CALLBACK(gtk_widget_destroy),
					 popup);
                gtk_widget_show_all(item);
                gtk_menu_shell_append(GTK_MENU_SHELL(popup), item);
            }

            gtk_menu_popup(GTK_MENU(popup), NULL, NULL, NULL, NULL,
                           event->button, event->time);
        }
        return TRUE;
    }

    return FALSE;
}

static gint lvars_right_click (GtkWidget *widget, GdkEventButton *event,
                               selector *sr)
{
    if (right_click(event)) {
        if (widget == sr->lvars2) {
            if (sr->ci == MIDASREG) {
                add_to_rvars2_callback(widget, sr);
            }
        } else {
            if (NONPARAM_CODE(sr->ci)) {
                set_extra_var_callback(NULL, sr);
            } else if (sr->ci == GR_FBOX || sr->ci == FSUMMARY) {
                set_third_var_callback(NULL, sr);
            } else {
                add_to_rvars1_callback(NULL, sr);
            }
        }
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
                       COL_ID, 0, COL_LAG, 0, COL_NAME, "const", -1);
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

    if (GTK_IS_ENTRY(sr->rvars1)) {
        /* clear special slot */
        gtk_entry_set_text(GTK_ENTRY(sr->rvars1), "");
    } else if (sr->rvars1 != NULL) {
        /* empty upper right variable list */
        clear_varlist(sr->rvars1);
        if (sr->add_button != NULL) {
            gtk_widget_set_sensitive(sr->add_button, TRUE);
        }
    }

    if (MODEL_CODE(sr->ci) && sr->ci != ARMA &&
	sr->ci != GARCH && sr->ci != REGLS) {
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

    if (w == NULL || !gtk_widget_is_sensitive(w)) {
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
                    gtk_tree_model_get(mod, &iter, COL_ID, &v, COL_LAG, &lag, -1);
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
    case INTREG:
        warnbox(_("You must select a lower bound variable"));
        break;
    default:
        warnbox(_("You must select a dependent variable"));
    }
}

static void reverse_list (char *list)
{
    char istr[VNAMELEN];
    char *tmp, *p;

    p = strchr(list, ';');
    if (p == NULL) return;

    tmp = malloc(strlen(list) + 4);
    if (tmp == NULL) return;

    sscanf(list, "%31s", istr);

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
    size_t addlen = strlen(add);
    size_t req = sr->cmdlen + addlen + 1;
    int err = 0;

    if (req > sr->cmdsize) {
        size_t newsize = sr->cmdsize + MAXLEN;
        char *tmp = NULL;

        if (newsize < req) {
            newsize = req;
        }

        tmp = realloc(sr->cmdlist, newsize);
        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            sr->cmdlist = tmp;
            sr->cmdsize = newsize;
        }
    }

    if (!err) {
        strcat(sr->cmdlist, add);
        sr->cmdlen += addlen;
    }

    return err;
}

/* append a space followed by either the ID number or
   the name of the series */

static void cmdlist_append_series (selector *sr,
                                   const char *s0,
                                   int id)
{
    char idstr[8];

    if (s0 != NULL) {
        add_to_cmdlist(sr, s0);
    }

    if (LIST_USE_INTS(sr->ci)) {
        sprintf(idstr, "%d", id);
        add_to_cmdlist(sr, idstr);
    } else if (id > 0 && id < dataset->v) {
        add_to_cmdlist(sr, dataset->varname[id]);
    } else {
        sprintf(idstr, "%d", id);
        add_to_cmdlist(sr, idstr);
    }
}

static void print_list_element (selector *sr,
                                char *targ,
                                const char *s0,
                                int id)
{
    if (LIST_USE_INTS(sr->ci)) {
        sprintf(targ, "%s%d", s0, id);
    } else if (id > 0 && id < dataset->v) {
        sprintf(targ, "%s%s", s0, dataset->varname[id]);
    } else {
        sprintf(targ, "%s%d", s0, id);
    }
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

    gretl_charsub(targ, ',', ' ');

    return targ;
}

static int widget_is_relevant (GtkWidget *w)
{
    if (w == NULL) {
        return 0;
    } else {
        return gtk_widget_is_sensitive(w);
    }
}

static void arma_spec_to_cmdlist (selector *sr)
{
    const char *txt;
    char s[32];

    free(arlags);
    arlags = NULL;

    free(malags);
    malags = NULL;

    if (widget_is_relevant(sr->extra[ARMA_plist])) {
        /* "gappy" AR lags activated */
        txt = gtk_entry_get_text(GTK_ENTRY(sr->extra[ARMA_plist]));
        add_to_cmdlist(sr, arma_lag_string(s, txt));
        arlags = gretl_strdup(txt);
    } else {
        /* regular max AR lag */
        arma_p = spin_get_int(sr->extra[ARMA_p]);
        sprintf(s, "%d ", arma_p);
        add_to_cmdlist(sr, s);
    }

    arima_d = spin_get_int(sr->extra[ARIMA_d]);
    sprintf(s, "%d ", arima_d);
    add_to_cmdlist(sr, s);

    if (widget_is_relevant(sr->extra[ARMA_qlist])) {
        /* "gappy" MA lags activated */
        txt = gtk_entry_get_text(GTK_ENTRY(sr->extra[ARMA_qlist]));
        add_to_cmdlist(sr, arma_lag_string(s, txt));
        add_to_cmdlist(sr, "; ");
        malags = gretl_strdup(txt);
    } else {
        /* regular max MA lag */
        arma_q = spin_get_int(sr->extra[ARMA_q]);
        sprintf(s, "%d ; ", arma_q);
        add_to_cmdlist(sr, s);
    }

    if (sr->extra[ARMA_P] != NULL) {
        arma_P = spin_get_int(sr->extra[ARMA_P]);
        arima_D = spin_get_int(sr->extra[ARIMA_D]);
        arma_Q = spin_get_int(sr->extra[ARMA_Q]);

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
        garch_p = spin_get_int(sr->extra[0]);
        garch_q = spin_get_int(sr->extra[1]);
        sprintf(s, "%d %d ; ", garch_p, garch_q);
    } else if (sr->ci == DPANEL) {
        dpd_p = spin_get_int(sr->extra[0]);
        sprintf(s, "%d ; ", dpd_p);
    } else if (sr->ci == ARCH) {
        int p = spin_get_int(sr->extra[0]);

        sprintf(s, "%d ", p);
    }

    add_to_cmdlist(sr, s);
}

static void read_ellipse_alpha (selector *sr)
{
    if (sr->extra[0] != NULL) {
        char s[16];
        double cval;

        cval = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[0]));
        sprintf(s, "%g", 1 - cval);
        add_to_cmdlist(sr, s);
    }
}

static void read_coint_opt_parm (selector *sr)
{
    GtkWidget *combo = sr->extra[EXTRA_LAGS];

    if (combo != NULL && gtk_widget_is_sensitive(combo)) {
        int method = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));

        if (method == 1) {
            sr->extra_data = gretl_strdup("BIC");
        } else if (method == 2) {
            sr->extra_data = gretl_strdup("tstat");
        }
    }
}

/* Take the stored preferred laglist for a variable (if any) and
   convert to a string specification for adding to the regression
   command line.
*/

static gchar *
discrete_lags_string (const char *vname, const int *laglist,
                      char context)
{
    gchar *ret = NULL;
    int nlags = laglist[0];
    int i, li, len;
    char tmp[64];

    /* allow 3 for space, '(' and ')' */
    len = 1 + nlags * (strlen(vname) + 3);

    for (i=1; i<=nlags; i++) {
        li = laglist[i];
        sprintf(tmp, "%d", li);
        len += strlen(tmp) + (li > 0);
    }

    ret = g_malloc0(len);

    if (ret != NULL) {
        for (i=1; i<=nlags; i++) {
            li = laglist[i];
            sprintf(tmp, " %s(%s%d)", vname, (li > 0)? "-" : "", li);
            strcat(ret, tmp);
        }
    }

    return ret;
}

int selector_get_VAR_order (const selector *sr)
{
    return spin_get_int(sr->extra[0]);
}

/* for use in constructing command list, possibly with
   embedded lags */

static char *get_lagpref_string (int v, char context,
                                 selector *sr)
{
    const char *vname = dataset->varname[v];
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
        s = g_strdup_printf(" %s", vname);
    }

#if LDEBUG
    if (s != NULL) {
        fprintf(stderr, "get_lagpref_string for v=%d (%s), context=%d:\n"
                " constructed s = '%s'\n", v, vname, (int) context, s);
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
            /* FIXME not needed for VECM? */
            newlist = gretl_list_resize(&vecxlist, n);
        }
        if (newlist == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static void nonparam_record_xvar (const char *s)
{
    int k;

    if (sscanf(s, "%d", &k) == 1 && k > 0) {
        np_xvar = k;
    }
}

static void get_rvars1_data (selector *sr, int rows, int context)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gint rvar, lag;
    gchar *rvstr;
    int added = 0;
    int i, j = 1;

    if ((SAVE_DATA_ACTION(sr->ci) || sr->ci == EXPORT) &&
        sr->ci != COPY_CSV && rows == sr->n_left) {
        /* saving/exporting all available series: leave the
           list blank in case it overflows
        */
        return;
    }

    sr->n_left = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<rows; i++) {
        gtk_tree_model_get(model, &iter, COL_ID, &rvar, COL_LAG, &lag, -1);

        if (is_lag_dummy(rvar, lag, context)) {
            gtk_tree_model_iter_next(model, &iter);
            continue;
        }

        if (context) {
            rvstr = get_lagpref_string(rvar, context, sr);
            if (rvstr == NULL) {
                sr->error = E_ALLOC;
                break;
            } else {
                add_to_cmdlist(sr, rvstr);
                g_free(rvstr);
                added++;
            }
        } else {
            cmdlist_append_series(sr, " ", rvar);
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

#if 0
    if (MODEL_CODE(sr->ci)) {
        printlist(xlist, "xlist");
    }
#endif

    if (sr->ci == ARMA && added && !(sr->opts & OPT_N)) {
        /* add const explicitly unless forbidden */
        add_to_cmdlist(sr, " 0");
    }
}

/* VECM: parse out the exogenous vars as either restricted
   or unrestricted */

static int get_vecm_exog_list (selector *sr, int rows,
                               GtkTreeModel *mod)
{
    GtkTreeIter iter;
    int *xlist = NULL;
    int *rlist = NULL;
    int i, v;
    int err = 0;

    gtk_tree_model_get_iter_first(mod, &iter);

    for (i=0; i<rows && !err; i++) {
        gchar *flag = NULL;
        int lag = 0;

        gtk_tree_model_get(mod, &iter, COL_ID, &v, COL_LAG, &lag,
                           COL_FLAG, &flag, -1);

        if (flag == NULL) {
            err = E_DATA;
        } else if (lag != 0) {
            v = laggenr(v, lag, dataset);
            if (v < 0) {
                err = E_DATA;
            }
        }

        if (!err) {
            if (*flag == 'U') {
                /* unrestricted terms */
                gretl_list_append_term(&xlist, v);
                if (xlist == NULL) {
                    err = E_ALLOC;
                }
            } else if (*flag == 'R') {
                /* restricted terms */
                gretl_list_append_term(&rlist, v);
                if (rlist == NULL) {
                    err = E_ALLOC;
                }
            } else {
                err = E_DATA;
            }
        }

        g_free(flag);
        gtk_tree_model_iter_next(mod, &iter);
    }

    if (!err) {
        if (xlist != NULL) {
            for (i=1; i<=xlist[0]; i++) {
                cmdlist_append_series(sr, " ", xlist[i]);
            }
        }
        if (rlist != NULL) {
            add_to_cmdlist(sr, " ;");
            for (i=1; i<=rlist[0]; i++) {
                cmdlist_append_series(sr, " ", rlist[i]);
            }
        }

        free(vecxlist);
        vecxlist = NULL;

        if (xlist != NULL && rlist != NULL) {
            vecxlist = gretl_lists_join_with_separator(xlist, rlist);
        } else if (xlist != NULL) {
            vecxlist = xlist;
            xlist = NULL;
        }
    }

    free(xlist);
    free(rlist);

    if (err) {
        gui_errmsg(err);
    }

    return err;
}

/* get the component of the model (or whatever) specification from the
   secondary list box on the right, if applicable */

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

        gtk_tree_model_get(model, &iter, COL_ID, &exog, COL_LAG, &lag, -1);

        if (IV_MODEL(sr->ci) && exog == ynum && lag == 0) {
            /* HECKIT? */
            errbox(_("You can't use the dependent variable as an instrument"));
            err = 1;
            break;
        }

        if (is_lag_dummy(exog, lag, context)) {
            gtk_tree_model_iter_next(model, &iter);
            continue;
        }

        if (context) {
            tmp = get_lagpref_string(exog, context, sr);
            add_to_cmdlist(sr, tmp);
            g_free(tmp);
        } else {
            cmdlist_append_series(sr, " ", exog);
        }

        if (reclist != NULL) {
            reclist[j++] = exog;
        }

        gtk_tree_model_iter_next(model, &iter);
    }

    return err;
}

static void clear_midas_spec (void)
{
    *mds_listname = '\0';
    mds_quad[0] = 0;
    mds_quad[1] = 0;
    mds_quad[2] = 1;
    mds_quad[3] = 2;
    mds_order = 1;
}

static void get_midas_specs (selector *sr)
{
    gui_midas_spec *specs;
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *vname;
    int i, rows = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            rows++;
        } while (gtk_tree_model_iter_next(model, &iter));
    }

    if (rows == 0) {
        warnbox("You must specify at least one MIDAS term");
        sr->error = 1;
        return;
    }

    specs = malloc(rows * sizeof *specs);
    if (specs == NULL) {
        nomem();
        sr->error = E_ALLOC;
        return;
    }

    gtk_tree_model_get_iter_first(model, &iter);

    for (i=0; i<rows; i++) {
        specs[i].nterms = rows;
        gtk_tree_model_get(model, &iter,
                           COL_ID, &specs[i].leadvar,
                           MCOL_M, &specs[i].fratio,
                           MCOL_NAME, &vname,
                           MCOL_MINLAG, &specs[i].minlag,
                           MCOL_MAXLAG, &specs[i].maxlag,
                           MCOL_TYPE, &specs[i].ptype,
                           MCOL_K, &specs[i].nparm, -1);
        if (i == 0) {
            /* remember some stuff */
            mds_quad[0] = specs[i].minlag;
            mds_quad[1] = specs[i].maxlag;
            mds_quad[2] = specs[i].ptype;
            mds_quad[3] = specs[i].nparm;
        }
        if (specs[i].leadvar > 0) {
            specs[i].listname[0] = '\0';
        } else {
            strcpy(specs[i].listname, vname);
            specs[i].leadvar = 0;
        }
        gtk_tree_model_iter_next(model, &iter);
    }

    if (sr->extra_data != NULL) {
        free(sr->extra_data);
    }

    sr->extra_data = specs;
}

static int check_midas_rvars2 (GtkTreeModel *model,
                               gboolean *have_beta1)
{
    GtkTreeIter iter;
    int ptype = -1;
    int nterms = 0;

    if (gtk_tree_model_get_iter_first(model, &iter)) {
        do {
            nterms++;
            if (have_beta1 != NULL) {
                gtk_tree_model_get(model, &iter, MCOL_TYPE,
                                   &ptype, -1);
                if (ptype == MIDAS_BETA1) {
                    *have_beta1 = 1;
                }
            }
        } while (gtk_tree_model_iter_next(model, &iter));
    }

    return nterms;
}

static gretl_bundle *regls_bundle;

void *selector_get_regls_bundle (void)
{
    return regls_bundle;
}

/* add "advanced" options from @src, if present */

static void regls_transcribe_advanced (gretl_bundle *rb,
                                       gretl_bundle *src,
                                       int xvalidate,
                                       int eid)
{
    int ccd = 0;

    if (gretl_bundle_get_int(src, "timer", NULL)) {
        gretl_bundle_set_int(rb, "timer", 1);
    }

    if ((eid == 0 && gretl_bundle_get_int(src, "lccd", NULL)) ||
        (eid == 1 && gretl_bundle_get_int(src, "rccd", NULL))) {
        ccd = 1;
    }
    gretl_bundle_set_int(rb, "ccd", ccd);

    if (xvalidate) {
        int use_1se = gretl_bundle_get_int(src, "use_1se", NULL);
        double s = gretl_bundle_get_scalar(src, "seed", NULL);

        gretl_bundle_set_int(rb, "use_1se", use_1se);
        if (gretl_bundle_get_int(src, "set_seed", NULL)) {
            gretl_bundle_set_scalar(rb, "seed", s);
        } else {
            gretl_bundle_delete_data(rb, "seed");
        }
    }

#ifdef HAVE_MPI
    if (xvalidate) {
        int no_mpi = gretl_bundle_get_int(src, "no_mpi", NULL);

        gretl_bundle_set_int(rb, "no_mpi", no_mpi);
    }
#endif
}

static void read_regls_extras (selector *sr)
{
    GtkWidget *w = sr->extra[REGLS_EST];
    gchar *estr = combo_box_get_active_text(w);
    gretl_bundle *rb = regls_bundle;
    int xvalidate = 0;
    int eid = 0;

    gretl_bundle_void_content(rb);
    gretl_bundle_set_int(rb, "gui", 1);

    if (!strcmp(estr, _("Elastic net"))) {
        GtkWidget *aspin = sr->extra[REGLS_ALPHA];
        double a = gtk_spin_button_get_value(GTK_SPIN_BUTTON(aspin));

        if (a == 1.0) {
            ; /* LASSO */
        } else if (a == 0) {
            gretl_bundle_set_int(rb, "ridge", 1);
            eid = 1;
        } else {
            gretl_bundle_set_scalar(rb, "alpha", a);
            eid = 2;
        }
    } else if (!strcmp(estr, _("Ridge"))) {
        gretl_bundle_set_int(rb, "ridge", 1);
        eid = 1;
    }

    if (gtk_widget_is_sensitive(sr->extra[REGLS_LAMVAL])) {
        GtkWidget *w = sr->extra[REGLS_LAMVAL];
        double lf = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));

        gretl_bundle_set_scalar(rb, "lfrac", lf);
    } else {
        int nlam = spin_get_int(sr->extra[REGLS_NLAM]);

        gretl_bundle_set_int(rb, "nlambda", nlam);
    }

    if (gtk_widget_is_sensitive(sr->extra[REGLS_NFOLDS])) {
        int nfolds = spin_get_int(sr->extra[REGLS_NFOLDS]);
        gchar *ft = combo_box_get_active_text(sr->extra[REGLS_FTYPE]);

        gretl_bundle_set_int(rb, "xvalidate", 1);
        gretl_bundle_set_int(rb, "nfolds", nfolds);
        if (!strcmp(ft, _("random"))) {
            gretl_bundle_set_int(rb, "randfolds", 1);
        }
        xvalidate = 1;
    }

    if (regls_adv != NULL) {
        regls_transcribe_advanced(rb, regls_adv, xvalidate, eid);
    }

    g_free(estr);
}

static void read_quantreg_extras (selector *sr)
{
    GtkWidget *w = gtk_bin_get_child(GTK_BIN(sr->extra[0]));
    const gchar *s = gtk_entry_get_text(GTK_ENTRY(w));

    if (s == NULL || *s == '\0') {
        warnbox(_("You must specify a quantile"));
        sr->error = 1;
    } else {
        /* convert the GUI string to what ought to be a valid
           numerical matrix specification */
        gchar *tmp = g_strdup_printf("{%s}", s);
        gretl_matrix *m;

        comma_separate_numbers(tmp);
        m = generate_matrix(tmp, dataset, &sr->error);
        gretl_matrix_free(m);

        if (sr->error) {
            warnbox(_("Invalid quantile specification"));
        } else {
            add_to_cmdlist(sr, tmp);
            add_to_cmdlist(sr, " ");
        }
        g_free(tmp);
    }

    w = sr->extra[1];

    if (!sr->error && w != NULL && gtk_widget_is_sensitive(w)) {
        GtkAdjustment *adj;

        adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(w));
        set_optval_double(QUANTREG, OPT_I, gtk_adjustment_get_value(adj));
    }
}

static void read_logistic_extras (selector *sr)
{
    GtkWidget *w = sr->extra[1];

    if (w != NULL && gtk_widget_is_sensitive(w)) {
        GtkAdjustment *adj;

        adj = gtk_spin_button_get_adjustment(GTK_SPIN_BUTTON(w));
        set_optval_double(sr->ci, OPT_M, gtk_adjustment_get_value(adj));
    }
}

static void read_omit_criterion (selector *sr)
{
    GtkWidget *w = sr->extra[OMIT_ALPHA];
    int err = 0;

    if (w != NULL && gtk_widget_is_sensitive(w)) {
        double val = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
        double orig = get_optval_double(OMIT, OPT_A, &err);

        if (val != orig) {
            set_optval_double(OMIT, OPT_A, val);
        }
        return;
    }

    w = sr->extra[OMIT_IC];
    if (w != NULL && gtk_widget_is_sensitive(w)) {
        gchar *val = combo_box_get_active_text(w);
        const char *orig = get_optval_string(OMIT, OPT_A);

        if (orig == NULL || strcmp(orig, val)) {
            set_optval_string(OMIT, OPT_A, val);
        }
        g_free(val);
    }
}

static void read_add_auto_param (selector *sr)
{
    if (sr->extra[0] != NULL && gtk_widget_is_sensitive(sr->extra[0])) {
        gchar *parm;

        if (gtk_widget_is_sensitive(sr->extra[1])) {
            double a = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[1]));
            gretl_push_c_numeric_locale();
            parm = g_strdup_printf("%g", a);
            gretl_pop_c_numeric_locale();
        } else {
            parm = combo_box_get_active_text(sr->extra[0]);
        }
        set_optval_string(ADD, OPT_A, parm);
        g_free(parm);
    }
}

static void read_reprobit_quadpoints (selector *sr)
{
    if (sr->extra[0] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[0])) {
        int qp = spin_get_int(sr->extra[0]);

        set_optval_int(PROBIT, OPT_G, qp);
    }
}

#define TOBIT_UNSET -1.0e300

static void read_tobit_limits (selector *sr)
{
    double lval = 0, rval = NADBL;
    const char *s;
    int err = 0;

    s = gtk_entry_get_text(GTK_ENTRY(sr->extra[0]));
    if (*s != '\0') {
        if (strcmp(s, "NA") == 0) {
            lval = TOBIT_UNSET;
        } else {
            err = check_atof(s);
            if (err) {
                warnbox(gretl_errmsg_get());
                sr->error = 1;
                return;
            } else {
                lval = atof(s);
            }
        }
    }

    s = gtk_entry_get_text(GTK_ENTRY(sr->extra[1]));
    if (*s != '\0' && strcmp(s, "NA")) {
        err = check_atof(s);
        if (err) {
            warnbox(gretl_errmsg_get());
            sr->error = 1;
            return;
        } else {
            rval = atof(s);
        }
    }

    /* record the user's choices */
    tobit_lo = (lval == TOBIT_UNSET)? NADBL : lval;
    tobit_hi = rval;

    if (lval == 0 && na(rval)) {
        ; /* the default, no-op */
    } else {
        if (lval != TOBIT_UNSET) {
            sr->opts |= OPT_L;
            set_optval_double(TOBIT, OPT_L, lval);
        }
        if (!na(rval)) {
            sr->opts |= OPT_M;
            set_optval_double(TOBIT, OPT_M, rval);
        }
    }
}

static void read_np_extras (selector *sr)
{
    char s[32];

    if (sr->ci == LOESS) {
        int d = spin_get_int(sr->extra[1]);
        double q;

        q = gtk_spin_button_get_value(GTK_SPIN_BUTTON(sr->extra[2]));
        sprintf(s, " d=%d q=%g", d, q);
        add_to_cmdlist(sr, s);
        if (button_is_active(sr->extra[3])) {
            sr->opts |= OPT_R;
        }
    } else if (sr->ci == NADARWAT) {
        GtkWidget *w = sr->extra[1];
        double h;

        if (w != NULL && gtk_widget_is_sensitive(w)) {
            h = gtk_spin_button_get_value(GTK_SPIN_BUTTON(w));
            sprintf(s, " h=%g", h);
            add_to_cmdlist(sr, s);
        }
        if (button_is_active(sr->extra[2])) {
            sr->opts |= OPT_O;
        }
    }
}

static int cluster_option_is_active (selector *sr)
{
    GtkWidget *hcb = sr->hccme_button;

    if (hcb != NULL && gtk_widget_is_sensitive(hcb) &&
        GTK_IS_BUTTON(hcb)) {
        const gchar *s = gtk_button_get_label(GTK_BUTTON(hcb));

        return s != NULL && strcmp(s, _("Cluster")) == 0;
    } else {
        return 0;
    }
}

static void maybe_read_cluster_var (selector *sr)
{
    if (cluster_option_is_active(sr)) {
        set_optval_string(sr->ci, OPT_C, cluster_var);
        sr->opts |= OPT_C;
    }
}

static void maybe_read_var_hac_option (selector *sr)
{
    GtkWidget *b = sr->hccme_button;

    if (b != NULL && gtk_widget_is_sensitive(b) &&
        GTK_IS_COMBO_BOX(b)) {
        gint i = gtk_combo_box_get_active(GTK_COMBO_BOX(b));

        if (i == 1) {
            sr->opts |= OPT_H;
        }
    }
}

#define extra_widget_get_int(c) (c == HECKIT ||         \
                                 c == BIPROBIT ||       \
                                 c == INTREG ||         \
                                 c == COUNTMOD ||       \
                                 c == DURATION ||       \
                                 c == WLS ||            \
                                 THREE_VARS_CODE(c) ||  \
                                 NONPARAM_CODE(c))

#define offer_cluster_option(c) (dataset_is_cross_section(dataset) &&   \
                                 cluster_option_ok(c))

static void parse_extra_widgets (selector *sr, char *endbit)
{
    const gchar *txt = NULL;
    int k = 0;

    if (offer_cluster_option(sr->ci)) {
        maybe_read_cluster_var(sr);
    }

    if (sr->ci == QUANTREG) {
        read_quantreg_extras(sr);
        return;
    } else if (sr->ci == REGLS) {
        read_regls_extras(sr);
        return;
    } else if (sr->ci == LOGISTIC || sr->ci == FE_LOGISTIC) {
        read_logistic_extras(sr);
        return;
    } else if (sr->ci == ELLIPSE) {
        read_ellipse_alpha(sr);
        return;
    } else if (sr->ci == OMIT) {
        read_omit_criterion(sr);
        return;
    } else if (sr->ci == ADD) {
        read_add_auto_param(sr);
        return;
    } else if (sr->ci == TOBIT) {
        read_tobit_limits(sr);
        return;
    } else if (sr->ci == REPROBIT) {
        read_reprobit_quadpoints(sr);
        return;
    }

    if (sr->ci == WLS || sr->ci == COUNTMOD || sr->ci == DURATION ||
        sr->ci == AR || sr->ci == HECKIT || sr->ci == BIPROBIT ||
        sr->ci == INTREG || THREE_VARS_CODE(sr->ci) ||
        NONPARAM_CODE(sr->ci)) {
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
            } else if (sr->ci == BIPROBIT) {
                warnbox(_("You must specify a second dependent variable"));
                sr->error = 1;
            } else if (sr->ci == INTREG) {
                warnbox(_("You must specify an upper bound variable"));
                sr->error = 1;
            } else if (sr->ci == ANOVA) {
                warnbox(_("You must specify a treatment variable"));
                sr->error = 1;
            } else if (NONPARAM_CODE(sr->ci)) {
                warnbox(_("You must specify an independent variable"));
                sr->error = 1;
            } else if (THREE_VARS_CODE(sr->ci)) {
                warnbox(("You must select a Y-axis variable"));
                sr->error = 1;
            } else if (sr->ci == COUNTMOD || sr->ci == DURATION) {
                /* empty 'extra' field is OK */
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
        /* treatment var */
        cmdlist_append_series(sr, " ", k);
    } else if (sr->ci == WLS || THREE_VARS_CODE(sr->ci)) {
        /* weight or y-axis var */
        cmdlist_append_series(sr, NULL, k);
        add_to_cmdlist(sr, " ");
        if (sr->ci == WLS) {
            wtvar = k;
        }
    } else if (sr->ci == INTREG) {
        /* upper-bound var */
        cmdlist_append_series(sr, " ", k);
        add_to_cmdlist(sr, " ");
        hivar = k;
    } else if (sr->ci == COUNTMOD) {
        /* offset variable */
        print_list_element(sr, endbit, " ; ", k);
        offvar = k;
    } else if (sr->ci == DURATION) {
        /* censoring variable */
        print_list_element(sr, endbit, " ; ", k);
        censvar = k;
    } else if (sr->ci == HECKIT) {
        /* selection variable */
        print_list_element(sr, endbit, " ", k);
        selvar = k;
    } else if (sr->ci == BIPROBIT) {
        /* depvar #2 */
        print_list_element(sr, endbit, " ", k);
    } else if (NONPARAM_CODE(sr->ci)) {
        /* independent var */
        print_list_element(sr, endbit, " ", k);
    } else if (sr->ci == AR) {
        /* lags */
        free(arlags);
        arlags = gretl_strdup(txt);
        add_to_cmdlist(sr, txt);
        add_to_cmdlist(sr, " ; ");
    }
}

static void vec_get_spin_data (selector *sr, int *order)
{
    const int *llist;
    char numstr[16];

    /* check for list of specific lags */
    llist = get_VAR_lags_list();

    if (llist != NULL) {
        /* "gappy" lag specification for VAR */
        char *dvlags;

        dvlags = gretl_list_to_lags_string(llist, &sr->error);
        if (dvlags != NULL) {
            add_to_cmdlist(sr, dvlags);
            free(dvlags);
        }
    } else {
        *order = spin_get_int(sr->extra[0]);
        sprintf(numstr, "%d", *order);
        add_to_cmdlist(sr, numstr);
    }

    if (sr->ci == VECM) {
        /* cointegration rank */
        jrank = spin_get_int(sr->extra[1]);
        sprintf(numstr, " %d", jrank);
        add_to_cmdlist(sr, numstr);
    }
}

static void parse_depvar_widget (selector *sr, char *endbit,
                                 char **dvlags,
                                 char **idvlags)
{
    int ynum = selector_get_depvar_number(sr);

    if (ynum < 0) {
        topslot_empty(sr->ci);
        sr->error = 1;
    } else {
        if (sr->ci == GR_XY || sr->ci == GR_IMP) {
            print_list_element(sr, endbit, " ", ynum);
        } else if (sr->ci == BIPROBIT) {
            cmdlist_append_series(sr, NULL, ynum);
            add_to_cmdlist(sr, endbit);
            *endbit = '\0';
        } else {
            cmdlist_append_series(sr, NULL, ynum);
        }

        if (select_lags_depvar(sr->ci) && dataset_lags_ok(dataset)) {
            *dvlags = get_lagpref_string(ynum, LAG_Y_X, NULL);
            if (USE_ZLIST(sr->ci)) {
                *idvlags = get_lagpref_string(ynum, LAG_Y_W, NULL);
            }
        }

        if (sr->default_check != NULL) {
            if (button_is_active(sr->default_check)) {
                default_y = ynum;
            } else {
                default_y = -1;
            }
        } else if (sr->ci == INTREG) {
            lovar = ynum;
        }
    }
}

static void parse_third_var_slot (selector *sr)
{
    const gchar *txt = gtk_entry_get_text(GTK_ENTRY(sr->rvars1));

    if (txt == NULL || *txt == '\0') {
        if (sr->ci == ANOVA) {
            /* third var is optional */
            return;
        } else if (sr->ci == GR_3D) {
            warnbox(_("You must select a Z-axis variable"));
        } else if (sr->ci == GR_DUMMY || sr->ci == GR_FBOX || sr->ci == FSUMMARY) {
            warnbox(_("You must select a factor variable"));
        } else {
            warnbox(_("You must select a control variable"));
        }
        sr->error = 1;
    } else {
        int v = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(sr->rvars1),
                                                  "data"));

        cmdlist_append_series(sr, " ", v);
    }
}

/* lag order for the dependent variable in midasreg */

static void midas_process_AR_spin (selector *sr)
{
    int yno = selector_get_depvar_number(sr);

    if (sr->extra[0] != NULL && yno > 0 && yno < dataset->v) {
        const char *yname = dataset->varname[yno];
        int p = spin_get_int(sr->extra[0]);
        gchar *bit = NULL;

        if (p == 1) {
            bit = g_strdup_printf(" %s(-1)", yname);
        } else if (p > 1) {
            bit = g_strdup_printf(" %s(-1 to -%d)", yname, p);
        }
        if (bit != NULL) {
            add_to_cmdlist(sr, bit);
            g_free(bit);
        }
        mds_order = p;
    }
}

static void selector_cancel_unavailable_options (selector *sr)
{
    if (sr->ci == ARMA) {
        if ((sr->opts & OPT_H) && !gtk_widget_is_sensitive(sr->hess_button)) {
            sr->opts ^= OPT_H;
        }
    } else if (sr->ci == CORR) {
        if (sr->opts & OPT_N && !gtk_widget_is_sensitive(sr->extra[0])) {
            sr->opts ^= OPT_N;
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
    char endbit[64] = {0};
    char *dvlags = NULL;
    char *idvlags = NULL;
    int context = 0;
    int order = 0;

    sr->error = 0;
    sr->cmdlist = mymalloc(MAXLEN);

    if (sr->cmdlist == NULL) {
        return;
    }

    sr->cmdsize = MAXLEN;
    sr->cmdlen = 0;
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
    if (ARMA_RELATED(sr->ci)) {
        arma_spec_to_cmdlist(sr);
    } else if (sr->ci == ARCH ||
               sr->ci == GARCH ||
               sr->ci == DPANEL) {
        add_pdq_vals_to_cmdlist(sr);
    } else if (VEC_CODE(sr->ci)) {
        vec_get_spin_data(sr, &order);
        if (!sr->error) {
            if (sr->ci == VAR) {
                maybe_read_var_hac_option(sr);
            } else if (sr->ci == COINT && !sr->error) {
                read_coint_opt_parm(sr);
            }
        }
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

    if (sr->ci == GR_FBOX || sr->ci == FSUMMARY || THREE_VARS_CODE(sr->ci)) {
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
    } else if (MODEL_CODE(sr->ci) && !NO_X_OK(sr->ci) && rows < 1) {
        warnbox(_("You must specify an independent variable"));
        sr->error = 1;
        return;
    }

    if (realrows > 0) {
        maybe_resize_recorder_lists(sr, realrows);
    }

    /* primary RHS varlist */
    if (sr->rvars1 != NULL) {
        context = sr_get_lag_context(sr, SR_RVARS1);
        get_rvars1_data(sr, rows, context);
    }

    if (sr->ci == MIDASREG) {
        /* FIXME placement of this? */
        midas_process_AR_spin(sr);
    }

    if (sr->ci == MIDASREG) {
        /* read special material from lower right */
        get_midas_specs(sr);
    } else if (USE_ZLIST(sr->ci) || USE_VECXLIST(sr->ci)) {
        /* cases with a (possibly optional) secondary RHS list */
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
        if (NONPARAM_CODE(sr->ci)) {
            read_np_extras(sr);
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
        } else if (sr->ci == ARMA) {
            arma_const = (sr->opts & OPT_N)? 0 : 1;
            arma_hessian = (sr->opts & OPT_G)? 0 : 1;
            arima_xdiff = (sr->opts & OPT_Y)? 0 : 1;
            arma_x12 = (sr->opts & OPT_X)? 1 : 0;
        } else  if (sr->ci == GARCH) {
            garch_const = (sr->opts & OPT_N)? 0 : 1;
        } else if (sr->ci == LOGIT || sr->ci == PROBIT) {
            lp_pvals = (sr->opts & OPT_P)? 1 : 0;
        } else if (sr->ci == DPANEL) {
            dpd_2step = (sr->opts & OPT_T)? 1 : 0;
            dpd_asy = (sr->opts & OPT_A)? 1 : 0;
            dpd_dpd = (sr->opts & OPT_X)? 1 : 0;
            dpd_coll = (sr->opts & OPT_C)? 1 : 0;
        } else if (NONPARAM_CODE(sr->ci)) {
            nonparam_record_xvar(endbit);
        }

        /* panel: scrub --nerlove if not doing random effects */
        if (sr->ci == PANEL) {
            if ((sr->opts & OPT_E) && !(sr->opts & OPT_U)) {
                sr->opts &= ~OPT_E;
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
    if (blocking(sr)) {
        gtk_main_quit();
    }

    if (state_pushed(sr)) {
        pop_program_state();
    }

    free(sr->cmdlist);
    free(sr->extra_data);
    free(sr);

    open_selector = NULL;
}

static char *estimator_label (int ci)
{
    switch (ci) {
    case OLS:
        return N_("OLS");
    case HSK:
        return N_("Heteroskedasticity corrected");
    case AR1:
        return N_("AR(1) errors");
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
    case REPROBIT:
        return N_("Random effects (binary) probit");
    case TOBIT:
        return N_("Tobit");
    case HECKIT:
        return N_("Heckit");
    case BIPROBIT:
        return N_("Bivariate probit");
    case LOGISTIC:
        return N_("Logistic");
    case FE_LOGISTIC:
        return N_("Fixed effects logistic model");
    case COUNTMOD:
        return N_("Count data model");
    case DURATION:
        return N_("Duration model");
    case PANEL:
        return N_("Panel model");
    case PANEL_WLS:
        return N_("Groupwise WLS");
    case PANEL_B:
        return N_("Between-groups model");
    case DPANEL:
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
        return N_("Autoregressive errors");
    case ARMA:
        return N_("ARIMA");
    case ALAGSEL:
        return N_("ARIMA lag selection");
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
        return N_("Cointegration (Engle-Granger)");
    case COINT2:
        return N_("Cointegration (Johansen)");
    case MPOLS:
        return N_("Multiple precision OLS");
    case LOESS:
        return N_("Loess");
    case NADARWAT:
        return N_("Nadaraya-Watson");
    case MIDASREG:
        return N_("MIDAS regression");
    case REGLS:
        return N_("Regularized least squares");
    default:
        return "";
    }
}

static char *extra_var_string (int ci)
{
    switch (ci) {
    case WLS:
        return N_("Weight variable");
    case COUNTMOD:
        return N_("Offset variable");
    case DURATION:
        return N_("Censoring variable");
    case HECKIT:
        return N_("Selection variable");
    case BIPROBIT:
        return N_("Dependent variable 2");
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
    case LOESS:
    case NADARWAT:
        return N_("Independent variable");
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
    GtkWidget *w = gtk_combo_box_text_new();

    combo_box_append_text(w, _("Y-axis variable"));
    combo_box_append_text(w, _("X-axis variable"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
                     G_CALLBACK(flip_multiplot_axis), NULL);

    multiplot_menu = w;

    return w;
}

static gint set_count_data_option (GtkComboBox *box, selector *sr)
{
    gint i = gtk_combo_box_get_active(box);

    if (i == 0) {
        /* Poisson */
        sr->opts &= ~(OPT_N | OPT_M);
    } else if (i == 1) {
        /* NegBin 2 */
        sr->opts &= ~OPT_M;
        sr->opts |= OPT_N;
    } else {
        /* NegBin 1 */
        sr->opts |= OPT_M;
    }

    return FALSE;
}

static void build_count_data_popdown (selector *sr)
{
    GtkWidget *w = gtk_combo_box_text_new();
    GtkWidget *hbox, *label;

    combo_box_append_text(w, _("Poisson"));
    combo_box_append_text(w, _("NegBin 2"));
    combo_box_append_text(w, _("NegBin 1"));

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
                     G_CALLBACK(set_count_data_option), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Distribution:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);

    if (model_opt & OPT_M) {
        sr->opts |= OPT_M;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 2);
    } else if (model_opt & OPT_N) {
        sr->opts |= OPT_N;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 1);
    } else {
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    }
}

static gint set_duration_option (GtkComboBox *box, selector *sr)
{
    gint i = gtk_combo_box_get_active(box);

    if (i == 0) {
        /* Weibull */
        sr->opts &= ~(OPT_E | OPT_L | OPT_Z);
    } else if (i == 1) {
        /* Exponential */
        sr->opts &= ~(OPT_L | OPT_Z);
        sr->opts |= OPT_E;
    } else if (i == 2) {
        /* Log-logistic */
        sr->opts &= ~(OPT_E | OPT_Z);
        sr->opts |= OPT_L;
    } else {
        /* Log-normal */
        sr->opts &= ~(OPT_E | OPT_L);
        sr->opts |= OPT_Z;
    }

    return FALSE;
}

static void build_duration_popdown (selector *sr)
{
    GtkWidget *w = gtk_combo_box_text_new();
    GtkWidget *hbox, *label;

    combo_box_append_text(w, _("Weibull"));
    combo_box_append_text(w, _("Exponential"));
    combo_box_append_text(w, _("Log-logistic"));
    combo_box_append_text(w, _("Log-normal"));

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
                     G_CALLBACK(set_duration_option), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Distribution:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);

    if (model_opt & OPT_E) {
        sr->opts |= OPT_E;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 1);
    } else if (model_opt & OPT_L) {
        sr->opts |= OPT_L;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 2);
    } else if (model_opt & OPT_Z) {
        sr->opts |= OPT_Z;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 3);
    } else {
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    }
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
    GtkWidget *w = gtk_combo_box_text_new();
    GtkWidget *hbox;

    combo_box_append_text(w, _("One-step estimation"));
    combo_box_append_text(w, _("Two-step estimation"));
    combo_box_append_text(w, _("Iterated estimation"));

    g_signal_connect(G_OBJECT(GTK_COMBO_BOX(w)), "changed",
                     G_CALLBACK(set_gmm_est_option), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);

    if (model_opt & OPT_I) {
        sr->opts |= OPT_I;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 2);
    } else if (model_opt & OPT_T) {
        sr->opts |= OPT_T;
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 1);
    } else {
        gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    }
}

static void table_add_left (selector *sr,
                            GtkWidget *child,
                            int startrow,
                            int endrow)
{
    guint xpad = 4, ypad = 2;

    gtk_table_attach(GTK_TABLE(sr->table), child,
                     0, 1, startrow, endrow,
                     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
                     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
                     xpad, ypad);
}

static void maybe_add_row (selector *sr)
{
    if (sr->row + 1 > sr->n_rows) {
        sr->n_rows += 1;
        gtk_table_resize(GTK_TABLE(sr->table),
                         sr->n_rows, 3);
    }
}

static void table_add_mid (selector *sr,
                           GtkWidget *child)
{
    guint xpad = 4, ypad = 2;

    maybe_add_row(sr);
    gtk_table_attach(GTK_TABLE(sr->table), child,
                     1, 2, sr->row, sr->row+1,
                     0, 0,
                     xpad, ypad);
}

static void table_add_right (selector *sr,
                             GtkWidget *child,
                             int fixed)
{
    guint xpad = 4, ypad = 2;

    maybe_add_row(sr);
    gtk_table_attach(GTK_TABLE(sr->table), child,
                     2, 3, sr->row, sr->row+1,
                     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
                     fixed ? 0 : (GTK_EXPAND | GTK_SHRINK | GTK_FILL),
                     xpad, ypad);
    /* finished a row, so advance */
    sr->row += 1;
}

static void alt_table_add_left (selector *sr,
                                GtkWidget *child,
                                int startrow,
                                int endrow)
{
    guint xpad = 4, ypad = 2;

    gtk_table_attach(GTK_TABLE(sr->table), child,
                     0, 1, startrow, endrow,
                     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
                     GTK_EXPAND | GTK_SHRINK | GTK_FILL,
                     xpad, ypad);
    sr->row += 1;
}

static void table_add_vwedge (selector *sr)
{
    GtkWidget *h = gtk_hbox_new(FALSE, 0);

    maybe_add_row(sr);
    gtk_table_attach(GTK_TABLE(sr->table), h,
                     2, 3, sr->row, sr->row+1,
                     0, 0, 0, 2);
    sr->row += 1;
}

static void vbox_add_vwedge (GtkWidget *vbox)
{
    GtkWidget *h = gtk_hbox_new(FALSE, 0);

    gtk_box_pack_start(GTK_BOX(vbox), h, FALSE, FALSE, 2);
    gtk_widget_show(h);
}

static GtkWidget *pix_button (const guint8 *src, gchar *tip)
{
    GtkWidget *img, *button;
    GdkPixbuf *pbuf;

    pbuf = gdk_pixbuf_new_from_inline(-1, src, FALSE, NULL);
    img = gtk_image_new_from_pixbuf(pbuf);
    button = gtk_button_new();
    gtk_widget_set_size_request(button, BUTTON_WIDTH, -1);
    // gtk_container_add(GTK_CONTAINER(button), img);
    gtk_button_set_image(GTK_BUTTON(button), img);
    gretl_tooltips_add(button, tip);
    g_object_unref(pbuf);

    return button;
}

static GtkWidget *name_entry_in_hbox (GtkWidget **pentry)
{
    GtkWidget *hbox, *entry;

    hbox = gtk_hbox_new(FALSE, 0);
    entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN - 1);
    gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAME_WIDTH);
    gtk_box_pack_start(GTK_BOX(hbox), entry, TRUE, TRUE, 0);

    if (pentry != NULL) {
        *pentry = entry;
    }

    return hbox;
}

static GtkWidget *
entry_with_label_and_chooser (selector *sr,
                              gchar *label_string,
                              int label_active,
                              click_func cf)
{
    GtkWidget *tmp, *hbox, *entry;

    if (label_active) {
        tmp = multiplot_popdown(sr->ci);
        table_add_right(sr, tmp, 1);
    } else if (label_string != NULL) {
        tmp = gtk_label_new(label_string);
        table_add_right(sr, tmp, 1);
    }

    tmp = pix_button(choose_inline, _("Choose"));
    table_add_mid(sr, tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(cf), sr);

    hbox = name_entry_in_hbox(&entry);
    table_add_right(sr, hbox, 1);

    if (label_active || label_string != NULL) {
        if (cf != set_third_var_callback) {
            table_add_vwedge(sr);
        }
    }

    return entry;
}

static void build_x_axis_section (selector *sr, int v)
{
    if (sr->ci == SCATTERS) {
        sr->depvar = entry_with_label_and_chooser(sr, NULL, 1,
                                                  set_dependent_var_callback);
    } else if (sr->ci == GR_FBOX) {
        sr->depvar = entry_with_label_and_chooser(sr, _("Variable to plot"), 0,
                                                  set_dependent_var_callback);
    } else if (sr->ci == FSUMMARY) {
        sr->depvar = entry_with_label_and_chooser(sr, _("Primary variable"), 0,
                                                  set_dependent_var_callback);
    } else {
        sr->depvar = entry_with_label_and_chooser(sr, _("X-axis variable"), 0,
                                                  set_dependent_var_callback);
    }

    if (v > 0 && v < dataset->v) {
        gtk_entry_set_text(GTK_ENTRY(sr->depvar), dataset->varname[v]);
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

static int build_depvar_section (selector *sr, int preselect)
{
    GtkWidget *tmp, *hbox;
    int defvar;

    if (sr->ci == INTREG) {
        defvar = (lovar > 0 && lovar < dataset->v)? lovar : -1;
    } else {
        if (default_y >= dataset->v) {
            default_y = -1;
        }
        defvar = (preselect)? preselect : default_y;
    }

    if (sr->ci == INTREG) {
        tmp = gtk_label_new(_("Lower bound variable"));
    } else if (sr->ci == ANOVA) {
        tmp = gtk_label_new(_("Response variable"));
    } else if (sr->ci == BIPROBIT) {
        tmp = gtk_label_new(_("Dependent variable 1"));
    } else {
        tmp = gtk_label_new(_("Dependent variable"));
    }

    table_add_right(sr, tmp, 1);

    tmp = pix_button(choose_inline, _("Choose"));
    table_add_mid(sr, tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(set_dependent_var_callback), sr);

    hbox = name_entry_in_hbox(&sr->depvar);
    g_signal_connect(G_OBJECT(sr->depvar), "changed",
                     G_CALLBACK(maybe_activate_depvar_lags), sr);
    if (defvar >= 0) {
        gtk_entry_set_text(GTK_ENTRY(sr->depvar), dataset->varname[defvar]);
    }
    table_add_right(sr, hbox, 1);

    if (sr->ci != INTREG && sr->ci != BIPROBIT) {
        sr->default_check = gtk_check_button_new_with_label(_("Set as default"));
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->default_check),
                                     default_y >= 0);
        table_add_right(sr, sr->default_check, 1);
    }

    if (sr->ci != DPANEL && sr->ci != BIPROBIT && sr->ci != TOBIT) {
        table_add_vwedge(sr);
    }

    return defvar;
}

/* In case we have a saved preference for the max lag of the
   endogenous vars in a VAR, set via the lags dialog, update this
   value from the global spin button
*/

static void lag_order_sync (GtkSpinButton *b, selector *sr)
{
    if (gtk_widget_is_sensitive(GTK_WIDGET(b)) &&
        sr->extra[EXTRA_LAGS] != NULL) {
        int lmax = gtk_spin_button_get_value_as_int(b);

        set_VAR_max_lag(lmax);
    }
}

enum {
    LAG_ONLY,
    LAG_AND_RANK
};

static void lag_order_spin (selector *sr, int which)
{
    gdouble lag, minlag, maxlag;
    GtkWidget *hbox, *label, *spin;
    GtkAdjustment *adj;
    const char *labels[] = {
        N_("lag order:"),
        N_("rank:")
    };
    int i, nspin = (which == LAG_AND_RANK)? 2 : 1;

    maxlag = (dataset->n < 72)? (dataset->n / 2) : 36;
    minlag = (sr->ci == COINT)? 0 : 1;

    if (default_order > 0 && default_order <= maxlag) {
        lag = default_order;
    } else {
        lag = (dataset->pd > 12)? 12 : dataset->pd;
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
            label = gtk_label_new(_("maximum lag:"));
        } else {
            label = gtk_label_new(_(labels[i]));
        }
        if (i == 0) {
            /* lag order */
            adj = (GtkAdjustment *) gtk_adjustment_new(lag, minlag, maxlag,
                                                       1, 1, 0);
        } else {
            /* rank */
            adj = (GtkAdjustment *) gtk_adjustment_new(jrank, 1, 10,
                                                       1, 1, 0);
        }
        spin = gtk_spin_button_new(adj, 1, 0);
        if (i == 1) {
            gretl_tooltips_add(label, _("Number of cointegrating vectors"));
        }
        gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, 5);
        gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 0);
        if (i == 0) {
            /* cross-connect with lag preferences dialog */
            if (get_VAR_lags_list() != NULL) {
                gtk_widget_set_sensitive(spin, FALSE);
            }
            g_signal_connect(G_OBJECT(spin), "value-changed",
                             G_CALLBACK(lag_order_sync), sr);
        }
        sr->extra[i] = spin;
        table_add_right(sr, hbox, 1);
    }
}

static void AR_order_spin (selector *sr)
{
    GtkWidget *tmp, *hbox;
    GtkAdjustment *adj;
    gdouble val, minlag, maxlag;

    hbox = gtk_hbox_new(FALSE, 5);

    if (sr->ci == ARCH) {
        tmp = gtk_label_new(_("ARCH order:"));
        val = dataset->pd;
        minlag = 1;
        maxlag = 2 * dataset->pd;
    } else if (sr->ci == MIDASREG) {
        tmp = gtk_label_new(_("AR order:"));
        val = mds_order;
        minlag = 0;
        maxlag = 100;
    } else {
        /* dpanel */
        tmp = gtk_label_new(_("AR order:"));
        val = dpd_p;
        minlag = 1;
        maxlag = 10;
        if (maxlag < dataset->pd - 2) {
            maxlag = dataset->pd - 2;
        }
    }

    if (default_order > 0 && default_order <= maxlag) {
        val = default_order;
    }

    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_misc_set_alignment(GTK_MISC(tmp), 0.0, 0.5);
    adj = (GtkAdjustment *) gtk_adjustment_new(val, minlag, maxlag, 1, 1, 0);
    sr->extra[0] = gtk_spin_button_new(adj, 1, 0);
    gtk_box_pack_start(GTK_BOX(hbox), sr->extra[0], FALSE, FALSE, 5);

    table_add_right(sr, hbox, 1);
}

static void extra_plotvar_box (selector *sr)
{
    const gchar *label;

    if (sr->ci == GR_3D) {
        label = N_("Z-axis variable");
    } else if (sr->ci == GR_DUMMY || sr->ci == GR_FBOX || sr->ci == FSUMMARY) {
        label = _("Factor (discrete)");
    } else if (sr->ci == GR_XYZ) {
        label = N_("Control variable");
    } else if (sr->ci == ANOVA) {
        label = N_("Block variable (optional)");
    } else {
        return;
    }

    sr->rvars1 = entry_with_label_and_chooser(sr, _(label), 0,
                                              set_third_var_callback);
}

static int get_nonparam_xvar (void)
{
    if (np_xvar > 0) {
        return np_xvar;
    } else if (xlist != NULL) {
        int i;

        for (i=1; i<=xlist[0]; i++) {
            if (xlist[i] != 0) {
                return xlist[i];
            }
        }
    }

    return 0;
}

static int get_setvar_value (int v)
{
    if (v > 0 && v < dataset->v) {
        return v;
    } else {
        return 0;
    }
}

/* selector for auxiliary series of some kind */

static void extra_var_box (selector *sr)
{
    int setvar = 0;

    sr->extra[0] = entry_with_label_and_chooser(sr, NULL, 0,
                                                set_extra_var_callback);

    if (sr->ci == WLS) {
        setvar = get_setvar_value(wtvar);
    } else if (sr->ci == HECKIT) {
        setvar = get_setvar_value(selvar);
    } else if (sr->ci == INTREG) {
        setvar = get_setvar_value(hivar);
    } else if (sr->ci == COUNTMOD) {
        setvar = get_setvar_value(offvar);
    } else if (sr->ci == DURATION) {
        setvar = get_setvar_value(censvar);
    } else if (sr->ci == BIPROBIT) {
        setvar = get_setvar_value(y2var);
    } else if (NONPARAM_CODE(sr->ci)) {
        setvar = get_nonparam_xvar();
    }

    if (setvar > 0) {
        gtk_entry_set_text(GTK_ENTRY(sr->extra[0]), dataset->varname[setvar]);
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

static void push_pull_buttons (selector *sr,
                               void (*addfunc)(GtkWidget *, selector *),
                               void (*remfunc)(GtkWidget *, selector *),
                               GtkWidget **lags_button,
                               gretlopt opt)
{
    GtkWidget *align, *vbox;
    GtkWidget *button;
    int spacing = 5;

    if ((opt & OPT_A) || lags_button == NULL) {
        spacing = 15;
    }

    align = gtk_alignment_new(0.5, 0.5, 0, 0);
    gtk_alignment_set_padding(GTK_ALIGNMENT(align),
                              0, 0, 0, 0);
    vbox = gtk_vbox_new(TRUE, 5);

    /* "Add" button */
    button = pix_button(add_inline, _("Add"));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(addfunc), sr);
    g_signal_connect(G_OBJECT(sr->dlg), "key-press-event",
                     G_CALLBACK(accept_right_arrow), button);
    gtk_box_pack_start(GTK_BOX(vbox), button, 0, 0, spacing);
    if (opt & OPT_A) {
        sr->add_button = button;
    }

    /* "All" button? */
    if (SAVE_DATA_ACTION(sr->ci) || sr->ci == EXPORT ||
	(sr->ci == REGLS_PLOTSEL && sr->n_left < 21)) {
        button = gtk_button_new_with_label(_("All ->"));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(add_all_to_rvars1_callback), sr);
        gtk_box_pack_start(GTK_BOX(vbox), button, 0, 0, spacing);
    }

    /* "Remove" button */
    button = pix_button(remove_inline, _("Remove"));
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(remfunc), sr);
    g_signal_connect(G_OBJECT(sr->dlg), "key-press-event",
                     G_CALLBACK(accept_left_arrow), button);
    gtk_box_pack_start(GTK_BOX(vbox), button, 0, 0, spacing);
    if (opt & OPT_R) {
        sr->remove_button = button;
    }

    /* "Lags" button? */
    if (lags_button != NULL) {
        button = gtk_button_new_with_label(_("lags..."));
        g_signal_connect(G_OBJECT(button), "clicked",
                         G_CALLBACK(lags_dialog_driver), sr);
        gtk_widget_set_sensitive(button, FALSE);
        *lags_button = button;
        gtk_box_pack_start(GTK_BOX(vbox), button, 0, 0, spacing);
    }

    gtk_container_add(GTK_CONTAINER(align), vbox);
    table_add_mid(sr, align);
}

static void secondary_rhs_varlist (selector *sr)
{
    GtkTreeModel *mod;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *hbox;
    GtkWidget *tmp = NULL;
    GtkWidget **lptr = NULL;
    int i;

    if (USE_VECXLIST(sr->ci)) {
        tmp = gtk_label_new(_("Exogenous variables"));
    } else if (IV_MODEL(sr->ci)) {
        tmp = gtk_label_new(_("Instruments"));
    } else if (sr->ci == HECKIT) {
        tmp = gtk_label_new(_("Selection regressors"));
    } else if (sr->ci == BIPROBIT) {
        tmp = gtk_label_new(_("Equation 2 regressors"));
    } else if (sr->ci == MIDASREG) {
        tmp = gtk_label_new(_("High-frequency"));
    } else if (FNPKG_CODE(sr->ci)) {
        tmp = gtk_label_new(_("Helper functions"));
    }

    if (tmp != NULL) {
        table_add_right(sr, tmp, 1);
    }

    hbox = gtk_hbox_new(FALSE, 5);

    sr->rvars2 = var_list_box_new(GTK_BOX(hbox), sr, SR_RVARS2);
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));

    if (sr->ci == VAR || sr->ci == VECM || sr->ci == VLAGSEL) {
        lptr = &sr->lags_button;
    }

    /* add push-pull buttons */
    push_pull_buttons(sr, add_to_rvars2_callback,
                      remove_from_rvars2_callback,
                      lptr, OPT_NONE);

    store = GTK_LIST_STORE(mod);
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(mod, &iter);

    if (sr->ci == MIDASREG) {
        ; /* FIXME do something here? */
    } else if (USE_ZLIST(sr->ci) && instlist != NULL) {
        for (i=1; i<=instlist[0]; i++) {
            if (instlist[i] < dataset->v) {
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
    } else if (!VEC_CODE(sr->ci) && !FNPKG_CODE(sr->ci)) {
        list_append_var(mod, &iter, 0, sr, SR_RVARS2);
    }

    table_add_right(sr, hbox, 0);
}

static void make_tau_list (GtkWidget *box)
{
    combo_box_append_text(box, "0.25 0.50 0.75");
    combo_box_append_text(box, ".05, .25 .50 .75, .95");
    combo_box_append_text(box, ".1 .2 .3 .4 .5 .6 .7 .8 .9");
}

static void tobit_limits_selector (selector *sr)
{
    const char *bstrs[] = {
        N_("left bound"),
        N_("right bound")
    };
    GtkWidget *hbox, *entry, *label;
    gchar *valstr;
    double val;
    int i;

    for (i=0; i<2; i++) {
        hbox = gtk_hbox_new(FALSE, 5);
        entry = gtk_entry_new();
        gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
        val = (i == 0)? tobit_lo : tobit_hi;
        valstr = na(val)? g_strdup("NA") : g_strdup_printf("%g", val);
        gtk_entry_set_text(GTK_ENTRY(entry), valstr);
        g_free(valstr);
        gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
        gtk_box_pack_end(GTK_BOX(hbox), entry, FALSE, FALSE, 5);
        label = gtk_label_new(_(bstrs[i]));
        gtk_box_pack_end(GTK_BOX(hbox), label, FALSE, FALSE, 5);
        table_add_right(sr, hbox, 1);
        sr->extra[i] = entry;
    }
}

static void add_np_controls (selector *sr)
{
    GtkAdjustment *adj;
    GtkWidget *hbox, *w;
    double bmin = 0.01;
    double bmax = 1.0;
    const char *optstr;
    int i = 1;

    if (sr->ci == LOESS) {
        /* polynomial order is specific to loess */
        hbox = gtk_hbox_new(FALSE, 5);
        w = gtk_label_new(_("Polynomial order"));
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        adj = (GtkAdjustment *) gtk_adjustment_new(1, 0, 2, 1, 1, 0);
        sr->extra[i] = gtk_spin_button_new(adj, 1, 0);
        gtk_box_pack_end(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
        table_add_right(sr, hbox, 1);
        i++;
    }

    /* bandwidth specification */
    if (sr->ci == LOESS) {
        hbox = gtk_hbox_new(FALSE, 5);
        w = gtk_label_new(_("Bandwidth"));
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        adj = (GtkAdjustment *) gtk_adjustment_new(0.5, bmin, bmax, 0.01, 0.1, 0);
        sr->extra[i] = gtk_spin_button_new(adj, 0.01, 2);
        gtk_box_pack_end(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
        table_add_right(sr, hbox, 1);
        i++;
    } else {
        double b0 = pow(sample_size(dataset), -0.2);
        GtkWidget *b1, *b2;
        GSList *group;

        hbox = gtk_hbox_new(FALSE, 5);
        b1 = gtk_radio_button_new_with_label(NULL, _("Automatic bandwidth"));
        gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 5);
        table_add_right(sr, hbox, 1);
        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
        hbox = gtk_hbox_new(FALSE, 5);
        b2 = gtk_radio_button_new_with_label(group,
                                             _("User-specified"));
        adj = (GtkAdjustment *) gtk_adjustment_new(b0, bmin, bmax, 0.01, 0.1, 0);
        w = sr->extra[i] = gtk_spin_button_new(adj, 0.01, 2);
        gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 5);
        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        table_add_right(sr, hbox, 1);
        gtk_widget_set_sensitive(w, FALSE);
        sensitize_conditional_on(w, b2);
        i++;
    }

    optstr = (sr->ci == LOESS)? N_("Use robust weights") :
        N_("Use \"leave one out\"");

    /* option checkbox */
    hbox = gtk_hbox_new(FALSE, 5);
    w = sr->extra[i] = gtk_check_button_new_with_label(_(optstr));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    table_add_right(sr, hbox, 1);
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

static void build_mid_section (selector *sr)
{
    const char *str = _(extra_var_string(sr->ci));
    GtkWidget *tmp;

    if (str != NULL) {
        tmp = gtk_label_new(str);
        table_add_right(sr, tmp, 1);
    }

    if (sr->ci == HECKIT || sr->ci == BIPROBIT) {
        extra_var_box(sr);
        table_add_vwedge(sr);
        primary_rhs_varlist(sr);
    } else if (sr->ci == WLS || sr->ci == INTREG ||
               sr->ci == COUNTMOD || sr->ci == DURATION ||
               THREE_VARS_CODE(sr->ci)) {
        extra_var_box(sr);
    } else if (NONPARAM_CODE(sr->ci)) {
        extra_var_box(sr);
        table_add_vwedge(sr);
        add_np_controls(sr);
    } else if (sr->ci == TOBIT) {
        tobit_limits_selector(sr);
        table_add_vwedge(sr);
    } else if (sr->ci == MIDASREG) {
        AR_order_spin(sr);
        primary_rhs_varlist(sr);
    } else if (USE_ZLIST(sr->ci)) {
        primary_rhs_varlist(sr);
    } else if (sr->ci == AR) {
        sr->extra[0] = gtk_entry_new();
        maybe_set_entry_text(sr->extra[0], arlags);
        table_add_right(sr, sr->extra[0], 1);
    } else if (sr->ci == QUANTREG) {
        GtkWidget *child;

        sr->extra[0] = combo_box_text_new_with_entry();
        make_tau_list(sr->extra[0]);
        child = gtk_bin_get_child(GTK_BIN(sr->extra[0]));
        gtk_entry_set_text(GTK_ENTRY(child), "0.5");
        table_add_right(sr, sr->extra[0], 1);
    } else if (sr->ci == VAR || sr->ci == VLAGSEL) {
        lag_order_spin(sr, LAG_ONLY);
        table_add_vwedge(sr);
        primary_rhs_varlist(sr);
    } else if (sr->ci == VECM) {
        lag_order_spin(sr, LAG_AND_RANK);
        table_add_vwedge(sr);
        primary_rhs_varlist(sr);
    } else if (sr->ci == COINT2) {
        lag_order_spin(sr, LAG_ONLY);
        table_add_vwedge(sr);
        primary_rhs_varlist(sr);
    } else if (VEC_CODE(sr->ci)) {
        lag_order_spin(sr, LAG_ONLY);
    } else if (sr->ci == DPANEL || sr->ci == ARCH) {
        AR_order_spin(sr);
    }

    table_add_vwedge(sr);
}

enum {
    SELECTOR_SIMPLE,
    SELECTOR_FULL
};

static GtkWidget *selector_dialog_new (selector *sr)
{
    GtkWidget *d = gretl_gtk_dialog();
    GtkWidget *base, *ca, *aa;

    g_signal_connect(G_OBJECT(d), "key-press-event",
                     G_CALLBACK(esc_kills_window), NULL);

    base = gtk_dialog_get_content_area(GTK_DIALOG(d));
    gtk_box_set_homogeneous(GTK_BOX(base), FALSE);
    gtk_box_set_spacing(GTK_BOX(base), 5);

    ca = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(base), ca, TRUE, TRUE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(ca), 5);
    gtk_box_set_spacing(GTK_BOX(ca), 5);

    aa = gtk_dialog_get_action_area(GTK_DIALOG(d));
    gtk_button_box_set_layout(GTK_BUTTON_BOX(aa), GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(aa), 10);
    gtk_container_set_border_width(GTK_CONTAINER(aa), 5);

    sr->vbox = ca;
    sr->action_area = aa;

    return d;
}

static int maybe_increase_vsize (selector *sr, float vsize)
{
    int ch = get_char_height(sr->dlg);
    float try = (ch / 18.0) * vsize;
    float adj = 0.5;
    int sh = get_screen_height();
    int ret = (int) vsize;

    if (try > vsize) {
        ret = (try <= adj * sh)? (int) try : (int) (adj * sh);
    }

    return ret;
}

static void selector_init (selector *sr, guint ci, const char *title,
                           sr_callback cb, GtkWidget *parent,
                           gpointer data, int selcode)
{
    int i, dlgx = -1, dlgy = 340;
    float x;

    sr->row = 0;
    sr->n_rows = 1;

    sr->ci = ci;
    sr->flags = 0;
    sr->opts = ci == TSPLOTS ? OPT_T : OPT_NONE;
    sr->parent = parent;
    sr->data = data;
    sr->extra_data = NULL;

    if (MODEL_CODE(ci) && dataset->v > 9) {
        dlgy += 80;
    }

    if (ci == ARMA) {
        dlgy += dataset->pd > 1 ? 140 : 80;
    } else if (ci == GARCH) {
        dlgy += 50;
    } else if (ci == WLS || ci == INTREG || ci == COUNTMOD ||
               ci == DURATION || ci == AR) {
        dlgy += 30;
    } else if (ci == HECKIT || ci == TOBIT) {
        dlgy += 80;
    } else if (ci == BIPROBIT) {
        dlgy += 110;
    } else if (IV_MODEL(ci)) {
        dlgy += 60;
    } else if (ci == ANOVA) {
        dlgy -= 60;
    } else if (ci == PANEL_WLS) {
        sr->opts |= OPT_H;
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
        if (ci == REGLS) {
            /* more stuff to show */
            dlgy += 240;
        } else {
            dlgy += 60;
        }
    }

    if (want_combo(sr)) {
        dlgy += 20;
    }

    if (dataset_lags_ok(dataset)) {
        if (MODEL_CODE(ci) && ci != ARMA) {
            /* lag selector button at foot */
            dlgy += 30;
        }
    }

    if (ci == DPANEL) {
        /* lots of option buttons */
        dlgy += 50;
    }

    sr->lvars = NULL;
    sr->lvars2 = NULL;
    sr->depvar = NULL;
    sr->rvars1 = NULL;
    sr->rvars2 = NULL;
    sr->default_check = NULL;
    sr->add_button = NULL;
    sr->remove_button = NULL;
    sr->lags_button = NULL;
    sr->hess_button = NULL;
    sr->x12a_button = NULL;
    sr->xdiff_button = NULL;
    sr->hccme_button = NULL;

    for (i=0; i<N_EXTRA; i++) {
        sr->extra[i] = NULL;
    }

    sr->cmdlist = NULL;
    sr->cmdlen = sr->cmdsize = 0;
    sr->callback = cb;

    sr->active_var = 0;
    sr->error = 0;
    sr->n_left = 0;

    if (ci == ARMA && push_program_state() == 0) {
        if (model_opt & OPT_L) {
            libset_set_bool(USE_LBFGS, 1);
        }
        sr->flags |= SR_STATE_PUSHED;
    }

    if (selcode == SELECTOR_SIMPLE) {
        sr->dlg = selector_dialog_new(sr);
        if (parent != NULL) {
            gtk_window_set_transient_for(GTK_WINDOW(sr->dlg),
                                         GTK_WINDOW(parent));
            gtk_window_set_destroy_with_parent(GTK_WINDOW(sr->dlg),
                                               TRUE);
        } else {
            gtk_window_set_transient_for(GTK_WINDOW(sr->dlg),
                                         GTK_WINDOW(mdata->main));
        }
    } else {
        sr->dlg = gretl_gtk_window();
        gretl_emulated_dialog_add_structure(sr->dlg,
                                            &sr->vbox,
                                            &sr->action_area);
    }

    gtk_window_set_title(GTK_WINDOW(sr->dlg), title);
    open_selector = sr;

    x = dlgy * gui_scale;
    dlgy = maybe_increase_vsize(sr, x);

    if (FNPKG_CODE(ci)) {
        x = 460 * gui_scale;
        dlgx = x;
    }

    gtk_window_set_default_size(GTK_WINDOW(sr->dlg), dlgx, dlgy);
#ifndef G_OS_WIN32
    set_wm_icon(sr->dlg);
#endif
    gtk_window_set_position(GTK_WINDOW(sr->dlg), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(sr->dlg), "destroy",
                     G_CALLBACK(destroy_selector),
                     sr);
}

static void option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (button_is_active(w)) {
        sr->opts |= opt;
    } else {
        sr->opts &= ~opt;
    }
}

static void reverse_option_callback (GtkWidget *w, selector *sr)
{
    gint i = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));
    gretlopt opt = i;

    if (button_is_active(w)) {
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

static void build_garch_spins (selector *sr)
{
    GtkWidget *tmp, *hbox;
    GtkAdjustment *adj;
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
        adj = (GtkAdjustment *) gtk_adjustment_new(val, 0, 4, 1, 1, 0);
        sr->extra[i] = gtk_spin_button_new(adj, 1, 0);
        gtk_box_pack_start(GTK_BOX(hbox), sr->extra[i], FALSE, FALSE, 5);
        g_signal_connect(GTK_SPIN_BUTTON(sr->extra[i]), "value-changed",
                         G_CALLBACK(garch_spin_check), sr);
    }

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
}

static GtkWidget *arma_aux_label (int ci, int i)
{
    GtkWidget *hbox;
    GtkWidget *lbl;
    const char *strs[] = {
        N_("Non-seasonal"),
        N_("Seasonal"),
        N_("Orders")
    };

    hbox = gtk_hbox_new(FALSE, 5);
    if (ci == ALAGSEL && i == 2) {
        lbl = gtk_label_new(_("Maxima"));
    } else {
        lbl = gtk_label_new(_(strs[i]));
    }
    gtk_box_pack_start(GTK_BOX(hbox), lbl, FALSE, FALSE, 5);

    return hbox;
}

static void toggle_p (GtkWidget *w, selector *sr)
{
    gboolean s = button_is_active(w);

    gtk_widget_set_sensitive(sr->extra[ARMA_p], !s);
    gtk_widget_set_sensitive(sr->extra[ARMA_plist], s);
}

static void toggle_q (GtkWidget *w, selector *sr)
{
    gboolean s = button_is_active(w);

    gtk_widget_set_sensitive(sr->extra[ARMA_q], !s);
    gtk_widget_set_sensitive(sr->extra[ARMA_qlist], s);
}

static void arima_callback (GtkWidget *w, selector *sr)
{
    gtk_widget_set_sensitive(sr->xdiff_button, arima_selected(sr) &&
                             rvars1_n_vars(sr) > 0);
}

static void build_arma_spins (selector *sr)
{
    GtkWidget *lbl, *chk, *tab;
    GtkAdjustment *adj;
    gdouble vmax, val;
    gboolean freeform;
    const char *strs[] = {
        N_("AR"),
        N_("I"),
        N_("MA")
    };
    int nrows, ncols = 7;
    int seasonals;
    int c, i, j;

    seasonals = dataset->pd > 1;
    nrows = seasonals ? 3 : 2;

    tab = gtk_table_new(nrows, ncols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tab), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tab), 5);
    gtk_box_pack_start(GTK_BOX(sr->vbox), tab, FALSE, FALSE, 0);

    lbl = arma_aux_label(sr->ci, seasonals? 0 : 2);
    gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 0, 1);
    c = 1;

    /* non-seasonal AR, I, MA spin-buttons */
    for (i=0; i<3; i++) {
        lbl = gtk_label_new(_(strs[i]));
        gtk_table_attach_defaults(GTK_TABLE(tab), lbl, c, c+1, 0, 1);
        c++;
        if (sr->ci == ALAGSEL) {
            val = (i==0)? 3 : (i==1)? 0 : 3;
        } else {
            val = (i==0)? arma_p : (i==1)? arima_d : arma_q;
        }
        vmax = (i == 1)? 2 : 10;
        adj = (GtkAdjustment *) gtk_adjustment_new(val, 0, vmax, 1, 1, 0);
        sr->extra[i] = gtk_spin_button_new(adj, 1, 0);
        if (i == 1) {
            g_signal_connect(G_OBJECT(sr->extra[i]), "value-changed",
                             G_CALLBACK(arima_callback), sr);
        }
        gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[i], c, c+1, 0, 1);
        c++;
    }

    if (sr->ci != ALAGSEL) {
        lbl = gtk_label_new(_("specific lags"));
        gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 1, 2);
        c = 1;
        j = 3;
        /* check buttons and entries for free-form lags */
        for (i=0; i<2; i++) {
            GCallback cb = (i == 0)? G_CALLBACK(toggle_p) : G_CALLBACK(toggle_q);
            char *fill = (i == 0)? arlags : malags;

            chk = gtk_check_button_new();
            g_signal_connect(G_OBJECT(chk), "clicked", cb, sr);
            gtk_table_attach_defaults(GTK_TABLE(tab), chk, c, c+1, 1, 2);
            c++;
            sr->extra[j] = gtk_entry_new();
            gtk_entry_set_max_length(GTK_ENTRY(sr->extra[j]), 16);
            gtk_entry_set_width_chars(GTK_ENTRY(sr->extra[j]), 8);
            gtk_widget_set_sensitive(sr->extra[j], FALSE);
            freeform = maybe_set_entry_text(sr->extra[j], fill);
            gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[j], c, c+1, 1, 2);
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk), freeform);
            j++;
            c += 3;
        }
    }

    if (seasonals) {
        gtk_table_set_row_spacing(GTK_TABLE(tab), 1, 10);
        lbl = arma_aux_label(sr->ci, 1);
        gtk_table_attach_defaults(GTK_TABLE(tab), lbl, 0, 1, 2, 3);
        c = 1;
        for (i=0; i<3; i++) {
            lbl = gtk_label_new(_(strs[i]));
            gtk_table_attach_defaults(GTK_TABLE(tab), lbl, c, c+1, 2, 3);
            c++;
            val = (i==0)? arma_P : (i==1)? arima_D : arma_Q;
            vmax = (i == 1)? 2 : 4;
            adj = (GtkAdjustment *) gtk_adjustment_new(val, 0, vmax, 1, 1, 0);
            sr->extra[j] = gtk_spin_button_new(adj, 1, 0);
            if (i == 1) {
                g_signal_connect(G_OBJECT(sr->extra[j]), "value-changed",
                                 G_CALLBACK(arima_callback), sr);
            }
            gtk_table_attach_defaults(GTK_TABLE(tab), sr->extra[j++], c, c+1, 2, 3);
            c++;
        }
    }
}

static void reset_arma_spins (selector *sr)
{
    if (sr->extra[ARMA_p] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARMA_p])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARMA_p]),
                                  (gdouble) arma_p);
    }
    if (sr->extra[ARIMA_d] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARIMA_d])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARIMA_d]),
                                  (gdouble) arima_d);
    }
    if (sr->extra[ARMA_q] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARMA_q])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARMA_q]),
                                  (gdouble) arma_q);
    }
    if (sr->extra[ARMA_plist] != NULL && GTK_IS_ENTRY(sr->extra[ARMA_plist])) {
        gtk_entry_set_text(GTK_ENTRY(sr->extra[ARMA_plist]), "");
    }
    if (sr->extra[ARMA_qlist] != NULL && GTK_IS_ENTRY(sr->extra[ARMA_qlist])) {
        gtk_entry_set_text(GTK_ENTRY(sr->extra[ARMA_qlist]), "");
    }
    if (sr->extra[ARMA_P] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARMA_P])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARMA_P]),
                                  (gdouble) arma_P);
    }
    if (sr->extra[ARIMA_D] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARIMA_D])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARIMA_D]),
                                  (gdouble) arima_D);
    }
    if (sr->extra[ARMA_Q] != NULL && GTK_IS_SPIN_BUTTON(sr->extra[ARMA_Q])) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[ARMA_Q]),
                                  (gdouble) arma_Q);
    }
}

static void hc_config (GtkWidget *w, selector *sr)
{
    if (offer_cluster_option(sr->ci)) {
        gboolean rconf = robust_conf(sr->ci);
        gretlopt c_opt = OPT_NONE;
        int resp;

        if (*cluster_var != '\0') {
            int v = current_series_index(dataset, cluster_var);

            if (v < 1) {
                *cluster_var = '\0';
            } else if (cluster_option_is_active(sr)) {
                c_opt = OPT_C;
            }
        }

        resp = hc_config_dialog(cluster_var, c_opt, rconf, sr->dlg);

        if (resp == 0) {
            /* regular HCCME */
            if (rconf) {
                preferences_dialog(TAB_VCV, NULL, sr->dlg);
            } else {
                gtk_button_set_label(GTK_BUTTON(sr->hccme_button),
                                     get_default_hc_string(sr->ci));
            }
        } else if (resp == 1) {
            /* cluster option selected */
            gtk_button_set_label(GTK_BUTTON(sr->hccme_button), _("Cluster"));
        }
    } else {
        preferences_dialog(TAB_VCV, NULL, sr->dlg);
    }
}

/* Callback for when the user selects something in the VCV
   tab of the preferences dialog: make sure that if there's
   a selection dialog open, and it's showing an HCCME
   button, the text on the button is synced with the
   user's selection.
*/

void selector_register_hc_choice (void)
{
    selector *sr = open_selector;

    if (sr != NULL && sr->hccme_button != NULL) {
        const char *txt = get_default_hc_string(sr->ci);
        GtkWidget *w = sr->hccme_button;

        if (GTK_IS_BUTTON(w)) {
            gtk_button_set_label(GTK_BUTTON(w), txt);
        } else if (GTK_IS_COMBO_BOX(w)) {
            gint i = gtk_combo_box_get_active(GTK_COMBO_BOX(w));

            combo_box_remove(w, 0);
            combo_box_prepend_text(w, txt);
            gtk_combo_box_set_active(GTK_COMBO_BOX(w), i);
        }
    }
}

static GtkWidget *pack_switch (GtkWidget *b, selector *sr,
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
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), checked);

    return hbox;
}

static void pack_switch_in (GtkWidget *hbox, GtkWidget *b, selector *sr,
                            gboolean checked, gboolean reversed,
                            gretlopt opt, int child)
{
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
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), checked);
}

static void call_iters_dialog (GtkWidget *w, GtkWidget *combo)
{
    static int optim = BFGS_MAX;
    selector *sr;
    int ci, active, maxit, lmem = 0;
    double tol;
    int resp;

    if (w == NULL) {
        /* clean-up signal */
        optim = libset_get_bool(USE_LBFGS)? LBFGS_MAX : BFGS_MAX;
        return;
    }

    sr = g_object_get_data(G_OBJECT(combo), "selector");
    active = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));

    ci = sr->ci == ALAGSEL ? ARMA : sr->ci;

    if (ci == ARMA && active == 1) {
        maxit = libset_get_int(BHHH_MAXITER);
        tol = libset_get_double(BHHH_TOLER);
        optim = BHHH_MAX;
    } else {
        optim = BFGS_MAX;
        BFGS_defaults(&maxit, &tol, ci);
        lmem = libset_get_int(LBFGS_MEM);
    }

    if (maxit <= 0) {
        maxit = 1000;
    }

    if (optim == BFGS_MAX && libset_get_bool(USE_LBFGS)) {
        optim = LBFGS_MAX;
    }

    resp = iter_control_dialog(&optim, &maxit, &tol,
                               &lmem, sr->dlg);

    if (!canceled(resp)) {
        int err = 0;

        if (optim == BHHH_MAX) {
            err = libset_set_int(BHHH_MAXITER, maxit);
            err += libset_set_double(BHHH_TOLER, tol);
        } else {
            err = libset_set_int(BFGS_MAXITER, maxit);
            err += libset_set_double(BFGS_TOLER, tol);
        }

        if (optim == LBFGS_MAX) {
            libset_set_int(LBFGS_MEM, lmem);
            sr->opts |= OPT_L;
        } else {
            sr->opts &= ~OPT_L;
        }

        if (err) {
            errbox("Error setting values");
        }
    }
}

#ifdef HAVE_X12A

static gboolean x12a_vs_native_opts (GtkWidget *w, selector *sr)
{
    gboolean s = button_is_active(w);

    if (sr->hess_button != NULL) {
        gtk_widget_set_sensitive(sr->hess_button, !s);
    }

    if (sr->xdiff_button != NULL) {
        if (s) {
            gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(sr->xdiff_button), s);
        }
        gtk_widget_set_sensitive(sr->xdiff_button, !s && arima_selected(sr));
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

    libset_set_int(GMP_BITS, b);
}

static GtkWidget *mpols_bits_selector (void)
{
    GtkWidget *w, *hbox, *combo;
    char bstr[16];
    int bits = get_mp_bits();
    int b, i, deflt = 0;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Bits per floating-point value"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);

    combo = gtk_combo_box_text_new();

    i = 0;
    for (b=256; b<=4096; b*=2) {
        sprintf(bstr, "%d", b);
        combo_box_append_text(combo, bstr);
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

static GtkWidget *pack_adf_test_down_sel (GtkWidget *hbox)
{
    GtkWidget *combo = gtk_combo_box_text_new();

    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 5);
    combo_box_append_text(combo, _("AIC"));
    combo_box_append_text(combo, _("BIC"));
    combo_box_append_text(combo, _("t-statistic"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    return combo;
}

void vbox_add_hwedge (GtkWidget *vbox)
{
    GtkWidget *h = gtk_hbox_new(FALSE, 0);

    gtk_box_pack_start(GTK_BOX(vbox), h, FALSE, FALSE, 2);
    gtk_widget_show(h);
}

static int seasonals_ok (void)
{
    if (dataset_is_time_series(dataset)) {
        if (dataset->pd == 4 || dataset->pd == 12) {
            return 1;
        } else if (dataset_is_daily(dataset)) {
            return 1;
        }
    }

    return 0;
}

static void build_selector_switches (selector *sr)
{
    int seas_ok = seasonals_ok();
    GtkWidget *hbox, *tmp;

    if (sr->ci == REPROBIT) {
        /* number of quadrature points for random effects probit */
        vbox_add_hwedge(sr->vbox);
        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_label_new(_("Quadrature points"));
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
        sr->extra[0] = tmp = gtk_spin_button_new_with_range(3, 64, 1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(tmp), 8);
        gtk_box_pack_start(GTK_BOX(hbox), tmp, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
        sr->opts |= OPT_G;
    }

    if (sr->ci == OLS || sr->ci == WLS || sr->ci == INTREG ||
        sr->ci == GARCH || sr->ci == IVREG || sr->ci == VAR ||
        sr->ci == LOGIT || sr->ci == PROBIT || sr->ci == MLOGIT ||
        sr->ci == OLOGIT || sr->ci == OPROBIT || sr->ci == COUNTMOD ||
        sr->ci == DURATION || sr->ci == PANEL || sr->ci == QUANTREG ||
        sr->ci == HECKIT || sr->ci == BIPROBIT || sr->ci == TOBIT ||
        sr->ci == MIDASREG || sr->ci == LOGISTIC ||
        sr->ci == FE_LOGISTIC) {
        GtkWidget *b1;

        /* FIXME arma robust variant? (and REPROBIT should be here?) */

        vbox_add_hwedge(sr->vbox);

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

        if (sr->ci == VAR) {
            /* special: we don't use HAC by default */
            GtkWidget *b2 = gtk_combo_box_text_new();

            sr->hccme_button = b2;
            combo_box_append_text(b2, get_default_hc_string(VAR));
            combo_box_append_text(b2, "HAC");
            gtk_combo_box_set_active(GTK_COMBO_BOX(b2), 0);
            gtk_widget_set_sensitive(b2, using_hc_by_default());
            sensitize_conditional_on(b2, b1);
            if (model_opt & OPT_R) {
                gtk_widget_set_sensitive(b2, TRUE);
            }
            gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 0);
        } else if (robust_conf(sr->ci) || offer_cluster_option(sr->ci)) {
            const char *deftxt;
            GtkWidget *b2;

            deftxt = get_default_hc_string(sr->ci);
            sr->hccme_button = b2 = gtk_button_new_with_label(deftxt);
            g_signal_connect(G_OBJECT(b2), "clicked",
                             G_CALLBACK(hc_config), sr);
            gtk_widget_set_sensitive(b2, using_hc_by_default());
            sensitize_conditional_on(b2, b1);
            if (model_opt & OPT_R) {
                gtk_widget_set_sensitive(b2, TRUE);
            }
            gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 0);
        }

        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    }

    if (ARMA_RELATED(sr->ci)) {
        hbox = gtk_hbox_new(FALSE, 5);
        vbox_add_vwedge(sr->vbox);
        tmp = gtk_check_button_new_with_label(_("Include a constant"));
        pack_switch_in(hbox, tmp, sr, arma_const, TRUE, OPT_N, 0);
        if (sr->ci == ARMA) {
            tmp = gtk_check_button_new_with_label(_("Show details of iterations"));
            pack_switch_in(hbox, tmp, sr, verbose, FALSE, OPT_V, 0);
        }
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    } else if (sr->ci == TOBIT || sr->ci == GARCH ||
               sr->ci == LOGIT || sr->ci == PROBIT || sr->ci == HECKIT ||
               sr->ci == OLOGIT || sr->ci == OPROBIT || sr->ci == MLOGIT ||
               sr->ci == BIPROBIT || sr->ci == REPROBIT) {
        if (sr->ci == GARCH) {
            tmp = gtk_check_button_new_with_label(_("Include a constant"));
            pack_switch(tmp, sr, garch_const, TRUE, OPT_N, 0);
        }
        tmp = gtk_check_button_new_with_label(_("Show details of iterations"));
        pack_switch(tmp, sr, verbose, FALSE, OPT_V, 0);
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
        pack_switch(tmp, sr, want_seasonals && seas_ok, FALSE, OPT_D, 0);
        gtk_widget_set_sensitive(tmp, seas_ok);
    } else if (sr->ci == COINT) {
        GtkWidget *combo;

        hbox = gtk_hbox_new(FALSE, 5);
        tmp = gtk_check_button_new_with_label(_("Test down from maximum lag order"));
        pack_switch_in(hbox, tmp, sr, FALSE, FALSE, OPT_E, 0);
        sr->extra[EXTRA_LAGS] = combo = pack_adf_test_down_sel(hbox);
        gtk_widget_set_sensitive(combo, button_is_active(tmp));
        sensitize_conditional_on(combo, tmp);
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
        tmp = gtk_check_button_new_with_label(_("Skip initial DF tests"));
        pack_switch(tmp, sr, FALSE, FALSE, OPT_S, 0);
        tmp = gtk_check_button_new_with_label(_("Show details of regressions"));
        pack_switch(tmp, sr, verbose, FALSE, OPT_V, 0);
        tmp = gtk_check_button_new_with_label(_("Include seasonal dummies"));
        pack_switch(tmp, sr, want_seasonals && seas_ok, FALSE, OPT_D, 0);
        gtk_widget_set_sensitive(tmp, seas_ok);
    } else if (sr->ci == PANEL_WLS) {
        tmp = gtk_check_button_new_with_label(_("Iterated weighted least squares"));
        pack_switch(tmp, sr, FALSE, FALSE, OPT_I, 0);
    } else if (sr->ci == PANEL || (sr->ci == OLS && dataset_is_panel(dataset))) {
        tmp = gtk_check_button_new_with_label(_("Include time dummies"));
        pack_switch(tmp, sr, (model_opt & OPT_D), FALSE, OPT_D, 0);
    } else if (sr->ci == DPANEL) {
        tmp = gtk_check_button_new_with_label(_("Include time dummies"));
        pack_switch(tmp, sr, (model_opt & OPT_D), FALSE, OPT_D, 0);
        tmp = gtk_check_button_new_with_label(_("Include levels equations (GMM-SYS)"));
        pack_switch(tmp, sr, (model_opt & OPT_L), FALSE, OPT_L, 0);
        tmp = gtk_check_button_new_with_label(_("Collapse instruments"));
        pack_switch(tmp, sr, (dpd_coll || (model_opt & OPT_C)), FALSE, OPT_C, 0);
        tmp = gtk_check_button_new_with_label(_("Asymptotic standard errors"));
        pack_switch(tmp, sr, (dpd_asy || (model_opt & OPT_A)), FALSE, OPT_A, 0);
        tmp = gtk_check_button_new_with_label(_("DPD-style initial covariance matrix"));
        pack_switch(tmp, sr, (dpd_dpd || (model_opt & OPT_X)), FALSE, OPT_X, 0);
    } else if (sr->ci == XTAB) {
        tmp = gtk_check_button_new_with_label(_("Show zeros explicitly"));
        pack_switch(tmp, sr, FALSE, FALSE, OPT_Z, 0);
    } else if (sr->ci == CORR) {
        sr->extra[0] = gtk_check_button_new_with_label(_("Ensure uniform sample size"));
        pack_switch(sr->extra[0], sr, verbose, FALSE, OPT_N, 0);
        gtk_widget_set_sensitive(sr->extra[0], get_n_rvars1() > 2);
    } else if (sr->ci == GR_3D) {
        tmp = gtk_check_button_new_with_label(_("Make plot interactive"));
        pack_switch(tmp, sr, TRUE, FALSE, OPT_I, 0);
    } else if (sr->ci == GR_BOX) {
        tmp = gtk_check_button_new_with_label(_("Show interval for median"));
        pack_switch(tmp, sr, FALSE, FALSE, OPT_O, 0);
    } else if (sr->ci == HSK) {
        tmp = gtk_check_button_new_with_label(_("Variance equation includes squares"));
        pack_switch(tmp, sr, TRUE, TRUE, OPT_N, 0);
    }

    if (sr->ci == ALAGSEL) {
        sr->xdiff_button = gtk_check_button_new_with_label(_("Difference the independent variables"));
        pack_switch(sr->xdiff_button, sr, arima_xdiff, TRUE, OPT_Y, 0);
        gtk_widget_set_sensitive(sr->xdiff_button, (arima_d > 0 || arima_D > 0) && !arma_x12);
    } else if (sr->ci == ARMA) {
        sr->hess_button =
            gtk_check_button_new_with_label(_("Parameter covariance matrix via Hessian"));
        pack_switch(sr->hess_button, sr, arma_hessian, TRUE, OPT_G, 0);
        sr->xdiff_button = gtk_check_button_new_with_label(_("Difference the independent variables"));
        pack_switch(sr->xdiff_button, sr, arima_xdiff, TRUE, OPT_Y, 0);
        gtk_widget_set_sensitive(sr->xdiff_button, (arima_d > 0 || arima_D > 0) && !arma_x12);
#ifdef HAVE_X12A
        sr->x12a_button = gtk_check_button_new_with_label(_("Use X-12-ARIMA"));
        pack_switch(sr->x12a_button, sr, arma_x12, FALSE, OPT_X, 0);
        g_signal_connect(G_OBJECT(sr->x12a_button), "toggled",
                         G_CALLBACK(x12a_vs_native_opts), sr);
#endif
    } else if (sr->ci == GARCH) {
        tmp = gtk_check_button_new_with_label(_("Standardize the residuals"));
        pack_switch(tmp, sr, (model_opt & OPT_Z), FALSE, OPT_Z, 0);
        tmp = gtk_check_button_new_with_label(_("Use Fiorentini et al algorithm"));
        pack_switch(tmp, sr, (model_opt & OPT_F), FALSE, OPT_F, 0);
    } else if (sr->ci == MPOLS) {
        hbox = mpols_bits_selector();
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
    } else if (sr->ci == MIDASREG) {
        tmp = gtk_check_button_new_with_label(_("Prefer NLS via Levenberg-Marquardt"));
        pack_switch(tmp, sr, (model_opt & OPT_L), FALSE, OPT_L, 0);
    } else if (sr->ci == SUMMARY || sr->ci == FSUMMARY) {
        tmp = gtk_check_button_new_with_label(_("Show all statistics"));
        pack_switch(tmp, sr, FALSE, TRUE, OPT_S, 0);
    }
}

static void unhide_lags_callback (GtkWidget *w, selector *sr)
{
    int i, imin, show_lags;
    GtkListStore *store;
    GtkTreeIter iter;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    imin = (sr->ci == DEFINE_LIST)? 0 : 1;
    show_lags = button_is_active(w);

    for (i=imin; i<dataset->v; i++) {
        if (list_show_var(sr, i, show_lags, NULL)) {
            list_append_var_simple(store, &iter, i);
        }
    }
}

static void unhide_lags_switch (selector *sr)
{
    GtkWidget *hbox;
    GtkWidget *b;

    vbox_add_vwedge(sr->vbox);

    b = gtk_check_button_new_with_label(_("Show lagged variables"));
    g_signal_connect(G_OBJECT(b), "toggled",
                     G_CALLBACK(unhide_lags_callback), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
}

#if 0 /* not ready */

static void boot_switch_callback (GtkWidget *w, selector *sr)
{
    if (button_is_active(w)) {
        sr->opts |= OPT_P;
    } else {
        sr->opts &= ~OPT_P;
    }
}

static void test_boot_switch (selector *sr)
{
    GtkWidget *hbox;
    GtkWidget *b;

    vbox_add_vwedge(sr->vbox);

    b = gtk_check_button_new_with_label(_("Use bootstrap"));
    g_signal_connect(G_OBJECT(b), "toggled",
                     G_CALLBACK(boot_switch_callback), sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), b, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
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

    if (extra_text != NULL) {
        GtkWidget *w = gtk_label_new(extra_text);

        gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    }

    gtk_box_pack_start(GTK_BOX(hbox), extra, TRUE, TRUE, offset);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b), checked);
}

static GtkWidget *alpha_spin (double deflt, double minval)
{
    GtkAdjustment *adj;

    adj = (GtkAdjustment *) gtk_adjustment_new(deflt, minval, 0.99,
                                               0.01, 0.1, 0);
    return gtk_spin_button_new(adj, 1, 2);
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
    sr->extra[1] = alpha_spin(0.90, 0.70);
    pack_switch_with_extra(b2, sr, FALSE, OPT_I, 0, sr->extra[1], "1 - α =");
    gtk_widget_set_sensitive(sr->extra[1], FALSE);
    sensitize_conditional_on(sr->extra[1], b2);
}

static GtkWidget *ymax_spin (void)
{
    GtkAdjustment *adj;

    adj = (GtkAdjustment *) gtk_adjustment_new(1.0, 0, 1e+10, 0.01, 0.1, 0);
    return gtk_spin_button_new(adj, 1, 1);
}

static GtkWidget *single_lambda_spin (double lam)
{
    GtkAdjustment *adj;

    adj = (GtkAdjustment *) gtk_adjustment_new(lam, 0, 1, 0.001, 0.1, 0);
    return gtk_spin_button_new(adj, 1, 3);
}

static GtkWidget *multi_lambda_spin (int nlam)
{
    GtkAdjustment *adj;

    adj = (GtkAdjustment *) gtk_adjustment_new(nlam, 4, 100, 1, 10, 0);
    return gtk_spin_button_new(adj, 1, 0);
}

static GtkWidget *regls_alpha_spin (double alpha)
{
    GtkAdjustment *adj;

    adj = (GtkAdjustment *) gtk_adjustment_new(alpha, 0, 1, 0.1, 0.1, 0);
    return gtk_spin_button_new(adj, 1, 1);
}

static void build_logistic_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    vbox_add_vwedge(sr->vbox);

    b1 = gtk_radio_button_new_with_label(NULL, _("Automatic maximum"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group,
                                         _("Specified maximum"));
    sr->extra[1] = ymax_spin();
    pack_switch_with_extra(b2, sr, FALSE, OPT_M, 0, sr->extra[1], NULL);
    gtk_widget_set_sensitive(sr->extra[1], FALSE);
    sensitize_conditional_on(sr->extra[1], b2);
}

static void regls_estim_switch (GtkComboBox *cb, selector *sr)
{
    GtkWidget *aspin = sr->extra[REGLS_ALPHA];
    gchar *estr = combo_box_get_active_text(cb);

    if (!strcmp(estr, _("LASSO"))) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(aspin), 1);
        gtk_widget_set_sensitive(aspin, FALSE);
    } else if (!strcmp(estr, _("Ridge"))) {
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(aspin), 0);
        gtk_widget_set_sensitive(aspin, FALSE);
    } else {
        /* elastic net */
        double a = gtk_spin_button_get_value(GTK_SPIN_BUTTON(aspin));

        if (a == 0 || a == 1) {
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(aspin), 0.5);
        }
        gtk_widget_set_sensitive(aspin, TRUE);
    }

    g_free(estr);
}

static void xvalidation_ok (GtkToggleButton *b, GtkWidget *targ)
{
    gtk_widget_set_sensitive(targ, button_is_active(b));
}

static void call_regls_advanced (GtkWidget *w, selector *sr)
{
    if (regls_adv == NULL) {
        double seed = (double) gretl_rand_get_seed();

        regls_adv = gretl_bundle_new();
        gretl_bundle_set_int(regls_adv, "lccd", 0);
        gretl_bundle_set_int(regls_adv, "rccd", 0);
        gretl_bundle_set_int(regls_adv, "use_1se", 0);
        gretl_bundle_set_int(regls_adv, "timer", 0);
        gretl_bundle_set_int(regls_adv, "set_seed", 0);
        gretl_bundle_set_scalar(regls_adv, "seed", seed);
#ifdef HAVE_MPI
        gretl_bundle_set_int(regls_adv, "no_mpi", 0);
#endif
    }

    regls_advanced_dialog(regls_adv, sr->dlg);
}

static int regls_int_default (const char *key)
{
    if (gretl_bundle_has_key(regls_bundle, key)) {
        return gretl_bundle_get_int(regls_bundle, key, NULL);
    } else if (!strcmp(key, "multi")) {
	return 1;
    } else if (!strcmp(key, "nlambda")) {
        return 25;
    } else if (!strcmp(key, "nfolds")) {
        return 10;
    } else if (!strcmp(key, "contiguous")) {
        return 1;
    } else if (!strcmp(key, "eid")) {
        if (gretl_bundle_get_int(regls_bundle, "ridge", NULL)) {
            return 1;
        } else if (gretl_bundle_has_key(regls_bundle, "alpha")) {
            return 2;
        }
    }

    return 0;
}

static double regls_scalar_default (const char *key)
{
    if (gretl_bundle_has_key(regls_bundle, key)) {
        return gretl_bundle_get_scalar(regls_bundle, key, NULL);
    } else if (!strcmp(key, "lfrac")) {
        return 0.5;
    } else if (!strcmp(key, "alpha")) {
        return 1.0;
    } else {
        return 0.0;
    }
}

static GtkWidget *cross_validation_options (selector *sr,
					    int nfolds,
					    int randfolds)
{
    GtkWidget *w, *hbox;
    int ids[2] = {
	REGLS_NFOLDS, REGLS_FTYPE
    };

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Folds:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    sr->extra[ids[0]] = w = gtk_spin_button_new_with_range(4, 20, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), nfolds);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    w = gtk_label_new(_("type:"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    sr->extra[ids[1]] = w = gtk_combo_box_text_new();
    combo_box_append_text(w, _("contiguous"));
    combo_box_append_text(w, _("random"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), randfolds? 1 : 0);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);

    return hbox;
}

static void build_regls_controls (selector *sr)
{
    GtkWidget *w, *hbox, *b1, *b2, *b3;
    int multi, nlambda, xvalidate, nfolds, randfolds;
    double lfrac, alpha;
    GSList *group;
    int eid;

    if (regls_bundle == NULL) {
        regls_bundle = gretl_bundle_new();
    }

    eid       = regls_int_default("eid");
    multi     = regls_int_default("multi");
    nlambda   = regls_int_default("nlambda");
    xvalidate = regls_int_default("xvalidate");
    nfolds    = regls_int_default("nfolds");
    randfolds = regls_int_default("randfolds");

    lfrac = regls_scalar_default("lfrac");
    alpha = regls_scalar_default("alpha");

    vbox_add_vwedge(sr->vbox);

    /* choice of estimator */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("Estimator"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    sr->extra[REGLS_EST] = w = gtk_combo_box_text_new();
    combo_box_append_text(w, _("LASSO"));
    combo_box_append_text(w, _("Ridge"));
    combo_box_append_text(w, _("Elastic net"));
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), eid);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    w = gtk_label_new("α =");
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    sr->extra[REGLS_ALPHA] = w = regls_alpha_spin(alpha);
    gtk_widget_set_sensitive(w, (eid == 2));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    g_signal_connect(G_OBJECT(sr->extra[REGLS_EST]), "changed",
                     G_CALLBACK(regls_estim_switch), sr);

    vbox_add_vwedge(sr->vbox);

    /* single lambda value */
    hbox = gtk_hbox_new(FALSE, 5);
    b1 = gtk_radio_button_new_with_label(NULL, _("Single λ-fraction"));
    gtk_box_pack_start(GTK_BOX(hbox), b1, FALSE, FALSE, 5);
    w = sr->extra[REGLS_LAMVAL] = single_lambda_spin(lfrac);
    sensitize_conditional_on(w, b1);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);

    /* multiple lambda values */
    hbox = gtk_hbox_new(FALSE, 5);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Multiple λ values"));
    gtk_box_pack_start(GTK_BOX(hbox), b2, FALSE, FALSE, 5);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), multi);
    w = sr->extra[REGLS_NLAM] = multi_lambda_spin(nlambda);
    gtk_widget_set_sensitive(w, multi);
    sensitize_conditional_on(w, b2);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);

    /* cross validation */
    hbox = gtk_hbox_new(FALSE, 5);
    b3 = gtk_check_button_new_with_label(_("Optimize via cross-validation"));
    gtk_widget_set_sensitive(b3, multi);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b3), xvalidate);
    gtk_box_pack_start(GTK_BOX(hbox), b3, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    hbox = cross_validation_options(sr, nfolds, randfolds);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    gtk_widget_set_sensitive(hbox, xvalidate);
    sensitize_conditional_on(hbox, b3);

    /* note: b2 = multiple lambdas, b3 = xvalidate */
    g_signal_connect(G_OBJECT(b2), "toggled",
                     G_CALLBACK(xvalidation_ok), b3);

    /* "advanced" controls */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_button_new_with_label(_("Advanced..."));
    g_signal_connect(G_OBJECT(w), "clicked",
                     G_CALLBACK(call_regls_advanced), sr);
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
}

static void build_ellipse_spin (selector *sr)
{
    GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
    GtkWidget *label;

    vbox_add_vwedge(sr->vbox);

    label = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, TRUE, 0);

    label = gtk_label_new(_("Confidence level"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    label = gtk_label_new("1 - α =");
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    sr->extra[0] = alpha_spin(0.95, 0.70);
    gtk_box_pack_start(GTK_BOX(hbox), sr->extra[0], FALSE, FALSE, 0);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
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

    if (lp_pvals) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    } else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    }
}

static void build_scatters_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Use points"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Use lines"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_O, 0);
}

static void auto_omit_restrict_callback (GtkWidget *w, selector *sr)
{
    gboolean r = button_is_active(w);

    if (sr->lvars != NULL) {
        gtk_widget_set_sensitive(sr->lvars, r);
    }

    if (sr->rvars1 != NULL) {
        gtk_widget_set_sensitive(sr->rvars1, r);
    }

    gtk_widget_set_sensitive(sr->add_button, r);
    gtk_widget_set_sensitive(sr->remove_button, r);
}

static void auto_omit_callback (GtkWidget *w, selector *sr)
{
    gboolean s = button_is_active(w);
    gboolean arrows = !s;

    if (w == sr->extra[OMIT_B3]) {
        /* the p-value criterion button */
        gtk_widget_set_sensitive(sr->extra[OMIT_ALPHA], s);
        if (sr->extra[OMIT_IC] != NULL) {
            gtk_widget_set_sensitive(sr->extra[OMIT_IC], !s);
        }
    } else {
        /* B4: the info criterion button */
        gtk_widget_set_sensitive(sr->extra[OMIT_IC], s);
        gtk_widget_set_sensitive(sr->extra[OMIT_ALPHA], !s);
    }

    if (sr->extra[OMIT_SEL] != NULL) {
        if (button_is_active(sr->extra[OMIT_SEL])) {
            arrows = TRUE;
        }
        gtk_widget_set_sensitive(sr->extra[OMIT_SEL], s);
    }

    if (sr->lvars != NULL) {
        gtk_widget_set_sensitive(sr->lvars, arrows);
    }

    if (sr->rvars1 != NULL) {
        gtk_widget_set_sensitive(sr->rvars1, arrows);
    }

    gtk_widget_set_sensitive(sr->add_button, arrows);
    gtk_widget_set_sensitive(sr->remove_button, arrows);
}

static void sensitize_alpha_spin (GtkComboBox *cb,
                                  GtkWidget *spin)
{
    gint i = gtk_combo_box_get_active(cb);

    gtk_widget_set_sensitive(spin, i == 3);
}

static void sensitize_stepwise (GtkToggleButton *tb,
                                selector *sr)
{
    if (gtk_toggle_button_get_active(tb)) {
        gint i = gtk_combo_box_get_active(GTK_COMBO_BOX(sr->extra[0]));

        gtk_widget_set_sensitive(sr->extra[0], TRUE);
        gtk_widget_set_sensitive(sr->extra[1], i == 3);
    } else {
        gtk_widget_set_sensitive(sr->extra[0], FALSE);
        gtk_widget_set_sensitive(sr->extra[1], FALSE);
    }
}

static void build_add_test_radios (selector *sr)
{
    const char *opts[] = {"AIC", "BIC", "HQC", "SSR"};
    GtkWidget *b1, *b2, *b3;
    GtkWidget *cb, *lbl, *spin;
    GtkWidget *hbox;
    GSList *group;
    int i;

    vbox_add_vwedge(sr->vbox);

    b1 = gtk_radio_button_new_with_label(NULL, _("Estimate augmented model"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("LM test using auxiliary regression"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_L, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    b3 = gtk_radio_button_new_with_label(group, _("Stepwise add, using"));
    hbox = pack_switch(b3, sr, FALSE, FALSE, OPT_A, 0);
    sr->extra[0] = cb = gtk_combo_box_text_new();
    for (i=0; i<4; i++) {
        combo_box_append_text(cb, opts[i]);
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(cb), 1);
    gtk_box_pack_start(GTK_BOX(hbox), cb, TRUE, TRUE, 0);
    lbl = gtk_label_new("α =");
    gtk_box_pack_start(GTK_BOX(hbox), lbl, TRUE, TRUE, 0);
    sr->extra[1] = spin = gtk_spin_button_new_with_range(0.01, 0.5, 0.01);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), 0.10);
    gtk_widget_set_sensitive(cb, FALSE);
    gtk_widget_set_sensitive(spin, FALSE);
    gtk_box_pack_start(GTK_BOX(hbox), spin, TRUE, TRUE, 0);
    g_signal_connect(GTK_COMBO_BOX(cb), "changed",
                     G_CALLBACK(sensitize_alpha_spin), spin);
    g_signal_connect(GTK_TOGGLE_BUTTON(b3), "toggled",
                     G_CALLBACK(sensitize_stepwise), sr);
}

static void build_omit_test_radios (selector *sr)
{
    windata_t *vwin = (windata_t *) sr->data;
    GtkWidget *b1, *b2, *b3, *b4;
    GSList *group;
    int auto_ok = 0;

    vbox_add_vwedge(sr->vbox);

    b1 = gtk_radio_button_new_with_label(NULL, _("Estimate reduced model"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Wald test, based on covariance matrix"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_W, 0);

    if (sr->ci != VAROMIT) {
        MODEL *pmod = (MODEL *) vwin->data;

        if (pmod != NULL) {
            auto_ok = (pmod->ci != PANEL);
        }
    }

    if (auto_ok) {
        GtkWidget *label = gtk_label_new(_("Or stepwise elimination"));
        GtkWidget *hbox = gtk_hbox_new(FALSE, 5);
        GtkWidget *aspin, *combo, *chk;

        vbox_add_vwedge(sr->vbox);
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);

        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
        b3 = gtk_radio_button_new_with_label(group, _("using two-sided p-value:"));
        sr->extra[OMIT_B3] = b3;
        g_signal_connect(G_OBJECT(b3), "toggled", G_CALLBACK(auto_omit_callback), sr);
        sr->extra[OMIT_ALPHA] = aspin = alpha_spin(0.10, 0.01);
        pack_switch_with_extra(b3, sr, FALSE, OPT_A, 0, aspin, NULL);
        gtk_widget_set_sensitive(aspin, FALSE);

        group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b3));
        b4 = gtk_radio_button_new_with_label(group, _("using information criterion: "));
        sr->extra[OMIT_B4] = b4;
        g_signal_connect(G_OBJECT(b4), "toggled", G_CALLBACK(auto_omit_callback), sr);
        sr->extra[OMIT_IC] = combo = gtk_combo_box_text_new();
        combo_box_append_text(combo, "AIC");
        combo_box_append_text(combo, "BIC");
        combo_box_append_text(combo, "HQC");
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
        pack_switch_with_extra(b4, sr, FALSE, OPT_A, 0, combo, NULL);
        gtk_widget_set_sensitive(combo, FALSE);

        chk = gtk_check_button_new_with_label(_("Test only selected variables"));
        sr->extra[OMIT_SEL] = chk;
        pack_switch(chk, sr, FALSE, FALSE, OPT_NONE, 0);
        gtk_widget_set_sensitive(chk, FALSE);
        g_signal_connect(G_OBJECT(chk), "toggled",
                         G_CALLBACK(auto_omit_restrict_callback), sr);
    }
}

static void select_re_method (GtkComboBox *box, selector *sr)
{
    int a = gtk_combo_box_get_active(box);

    if (a == 2) {
        sr->opts |= OPT_E;
        sr->opts &= ~OPT_X;
    } else {
        sr->opts &= ~OPT_E;
        if (a == 1) {
            sr->opts |= OPT_X;
        } else {
            sr->opts &= ~OPT_X;
        }
    }
}

static void build_panel_radios (selector *sr)
{
    GtkWidget *b1, *b2, *hbox, *w;
    GSList *group;
    gboolean fe = TRUE;

    if (sr->opts & OPT_H) {
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
    hbox = pack_switch(b2, sr, !fe, FALSE, OPT_U, 0);

    /* random effects transformation selector */
    w = gtk_combo_box_text_new();
    combo_box_append_text(w, _("Swamy-Arora"));
    combo_box_append_text(w, _("Swamy-Arora / Baltagi-Chang"));
    combo_box_append_text(w, _("Nerlove"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
    gtk_combo_box_set_active(GTK_COMBO_BOX(w), 0);
    gtk_widget_set_sensitive(w, !fe);
    g_signal_connect(G_OBJECT(w), "changed",
                     G_CALLBACK(select_re_method), sr);
    sensitize_conditional_on(w, b2);
}

static void build_dpanel_radios (selector *sr)
{
    GtkWidget *b1, *b2;
    GSList *group;
    gboolean twostep = FALSE;

    if (dpd_2step || (model_opt & OPT_T)) {
        twostep = TRUE;
    }

    b1 = gtk_radio_button_new_with_label(NULL, _("One-step estimation"));
    pack_switch(b1, sr, !twostep, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Two-step estimation"));
    pack_switch(b2, sr, twostep, FALSE, OPT_T, 0);
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
}

static void build_ar1_radios (selector *sr)
{
    GtkWidget *b1, *b2, *b3, *b4, *hbox;
    GSList *group;

    b1 = gtk_radio_button_new_with_label(NULL, _("Cochrane-Orcutt"));
    pack_switch(b1, sr, TRUE, FALSE, OPT_NONE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    b2 = gtk_radio_button_new_with_label(group, _("Prais-Winsten"));
    pack_switch(b2, sr, FALSE, FALSE, OPT_P, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    b3 = gtk_radio_button_new_with_label(group, _("Hildreth-Lu"));
    hbox = pack_switch(b3, sr, FALSE, FALSE, OPT_H, 0);

    b4 = gtk_check_button_new_with_label(_("Fine-tune using Cochrane-Orcutt"));
    pack_switch_in(hbox, b4, sr, TRUE, TRUE, OPT_B, 0);
    gtk_widget_set_sensitive(b4, FALSE);
    sensitize_conditional_on(b4, b3);

    if (model_opt & OPT_P) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    } else if (model_opt & OPT_H) {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b3), TRUE);
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b4),
                                     !(model_opt & OPT_B));
    } else {
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    }
}

static gboolean arma_estimator_switch (GtkComboBox *box, selector *sr)
{
    if (sr->hess_button != NULL) {
        GtkWidget *xb = sr->x12a_button;

        if (xb == NULL || !button_is_active(xb)) {
            gchar *s = combo_box_get_active_text(box);

            gtk_widget_set_sensitive(sr->hess_button,
                                     !strcmp(s, _("Exact Maximum Likelihood")));
            g_free(s);
        }
    }

    return FALSE;
}

static void build_arma_combo (selector *sr)
{
    GtkWidget *hbox, *combo, *button;
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

    combo = gretl_opts_combo_full(&arma_opts, deflt, NULL,
                                  G_CALLBACK(arma_estimator_switch),
                                  sr);
    g_object_set_data(G_OBJECT(combo), "selector", sr);

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);

    button = gtk_button_new_from_stock(GTK_STOCK_PREFERENCES);
    g_signal_connect(G_OBJECT(button), "clicked",
                     G_CALLBACK(call_iters_dialog), combo);
    gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 5);

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
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
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
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

    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
}

#define DSET_DB_OK(d) (d->pd == 1 || (d->structure == TIME_SERIES &&    \
                                      (d->pd == 4 || d->pd == 12)))

static void build_data_export_combo (selector *sr)
{
    GtkWidget *hbox, *label, *combo;
    static const char *opt_strs[] = {
        N_("CSV"),
        N_("GeoJSON"),
        N_("gretl datafile (.gdt)"),
        N_("gretl binary datafile (.gdtb)"),
        N_("gretl database (.bin)"),
        N_("space separated (R-friendly)"),
        N_("Octave"),
        N_("Stata"),
        N_("JMulTi"),
        N_("PcGive"),
        NULL
    };
    static gretlopt opts[] = {
        OPT_C,
        OPT_P,
        OPT_Z,
        OPT_B,
        OPT_D,
        OPT_R,
        OPT_M,
        OPT_S,
        OPT_J,
        OPT_G
    };
    static combo_opts export_opts;
    int deflt = 0;

    export_opts.strs = opt_strs;
    export_opts.vals = opts;
    export_opts.optp = &sr->opts;

    hbox = gtk_hbox_new(FALSE, 5);

    label = gtk_label_new(_("Select format"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);

    if (dataset->mapfile != NULL) {
        int masked[4] = {3, 4, 8, 9};

        deflt = 1;
        combo = gretl_opts_combo_masked(&export_opts, deflt, masked);
    } else if (DSET_DB_OK(dataset)) {
        int masked[2] = {1, 1};

        combo = gretl_opts_combo_masked(&export_opts, deflt, masked);
    } else {
        int masked[3] = {2, 1, 4};

        combo = gretl_opts_combo_masked(&export_opts, deflt, masked);
    }

    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
}

static void build_selector_radios (selector *sr)
{
    if (sr->ci == PANEL) {
        build_panel_radios(sr);
    } else if (sr->ci == DPANEL) {
        build_dpanel_radios(sr);
    } else if (sr->ci == SCATTERS) {
        build_scatters_radios(sr);
    } else if (sr->ci == ADD) {
        build_add_test_radios(sr);
    } else if (sr->ci == OMIT || sr->ci == VAROMIT) {
        build_omit_test_radios(sr);
    } else if (sr->ci == LOGIT || sr->ci == PROBIT) {
        build_pvalues_radios(sr);
    } else if (sr->ci == HECKIT) {
        build_heckit_radios(sr);
    } else if (sr->ci == XTAB) {
        build_xtab_radios(sr);
    } else if (sr->ci == PCA) {
        build_pca_radios(sr);
    } else if (sr->ci == QUANTREG) {
        build_quantreg_radios(sr);
    } else if (sr->ci == LOGISTIC || sr->ci == FE_LOGISTIC) {
        build_logistic_radios(sr);
    } else if (sr->ci == AR1) {
        build_ar1_radios(sr);
    } else if (sr->ci == REGLS) {
        build_regls_controls(sr);
    }
}

static void build_selector_combo (selector *sr)
{
    if (ARMA_RELATED(sr->ci)) {
        build_arma_combo(sr);
    } else if (sr->ci == COINT) {
        build_coint_combo(sr);
    } else if (sr->ci == VECM || sr->ci == COINT2) {
        build_vecm_combo(sr);
    } else if (sr->ci == IV_GMM) {
        build_gmm_popdown(sr);
    } else if (sr->ci == COUNTMOD) {
        build_count_data_popdown(sr);
    } else if (sr->ci == DURATION) {
        build_duration_popdown(sr);
    } else if (sr->ci == EXPORT) {
        build_data_export_combo(sr);
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
    gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
}

static int count_selected (GtkWidget *w)
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    int ns = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(w));

    if (gtk_tree_model_get_iter_first(model, &iter)) {
	do {
	    ns++;
	} while (gtk_tree_model_iter_next(model, &iter));
    }

    return ns;
}

static void selector_doit (GtkWidget *w, selector *sr)
{
    if (!FNPKG_CODE(sr->ci)) {
        /* we don't need this when gathering function names */
        compose_cmdlist(sr);
    }

    if (sr->ci == DEFINE_LIST && (sr->flags & SR_NO_SINGLE)) {
	/* respect a requirement coming from fncall.c */
	if (count_selected(sr->rvars1) < 2) {
	    warnbox(_("At least two series are needed"));
	    return;
	}
    }

    if (sr->error == 0) {
        int err = 0;

#ifdef __APPLE__
        /* Hiding the selector window prevents the "next" window
           (i.e. the one opened by sr->callback) from being hidden
           behind gretl's main window, on Apple's X11
        */
        if (open_selector != NULL) {
            gtk_widget_hide(sr->dlg);
        }
#endif

        err = sr->callback(sr);

        if (!err && open_selector != NULL) {
            gtk_widget_destroy(sr->dlg);
        }

#ifdef __APPLE__
        if (err && open_selector != NULL) {
            gtk_widget_show(sr->dlg);
        }
#endif
    }
}

static void build_selector_buttons (selector *sr)
{
    GtkWidget *tmp;

    if (sr->ci != PRINT && sr->ci != SUMMARY && !FNPKG_CODE(sr->ci) &&
        sr->ci != DEFINE_LIST && sr->ci != DEFINE_MATRIX &&
        sr->ci != ELLIPSE && sr->ci != CHOW && sr->ci != REGLS_PLOTSEL &&
        sr->ci != FSUMMARY && !SAVE_DATA_ACTION(sr->ci)) {
        /* add a Help button if appropriate */
        int ci = sr->ci;

        if (sr->ci == OLOGIT || sr->ci == MLOGIT) {
            ci = LOGIT;
        } else if (sr->ci == OPROBIT) {
            ci = PROBIT;
        } else if (sr->ci == FE_LOGISTIC) {
            ci = LOGISTIC;
        } else if (IV_MODEL(sr->ci)) {
            ci = IVREG;
        }
        context_help_button(sr->action_area, ci);
    }

    if (sr->ci != EDIT_FUNCTIONS) {
        tmp = gtk_button_new_from_stock(GTK_STOCK_CLEAR);
        gtk_widget_set_can_default(tmp, TRUE);
        gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
        g_signal_connect(G_OBJECT(tmp), "clicked",
                         G_CALLBACK(clear_vars), sr);
    }

    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(cancel_selector), sr);

    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_container_add(GTK_CONTAINER(sr->action_area), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(selector_doit), sr);
    gtk_widget_grab_default(tmp);
}

static int list_show_var (selector *sr, int v, int show_lags,
			  gchar **exclude)
{
    int ci = sr->ci;
    int i, ret = 1;

    if (exclude != NULL) {
	for (i=0; exclude[i] != NULL; i++) {
	    if (!strcmp(dataset->varname[v], exclude[i])) {
		return 0;
	    }
	}
    }

    if (ci == LOESS || ci == NADARWAT) {
        /* special: for nonparam models we should show
           lagged vars, since we don't display the full
           lag-selection mechanism */
        show_lags = 1;
    }

    if (v == 0 && (ci == DEFINE_LIST || ci == DEFINE_MATRIX)) {
        ;
    } else if (v == 0 && (!MODEL_CODE(ci) || ci == ARMA ||
			  ci == GARCH || ci == REGLS)) {
        ret = 0;
    } else if (v == 0 && (ci == LOESS || ci == NADARWAT)) {
        ret = 0;
    } else if (series_is_hidden(dataset, v)) {
        ret = 0;
    } else if (is_panel_group_names_series(dataset, v)) {
        ret = 0;
    } else if (!show_lags && series_get_lag(dataset, v)) {
        sr->flags |= SR_LAGS_HIDDEN;
        ret = 0;
    } else if (ci == XTAB) {
        ret = accept_as_discrete(dataset, v, 0);
    } else if (ci == MIDASREG && series_get_midas_period(dataset, v)) {
        ret = 0;
    }

    return ret;
}

/* callback from successful generation of a new variable or
   variables via click on the selector's "Add variable button"
*/

void selector_register_genr (int newvars, gpointer p)
{
    selector *sr = p;
    GtkListStore *store;
    GtkTreeIter iter;
    int i, v;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    tree_model_get_iter_last(GTK_TREE_MODEL(store), &iter);

    for (i=0; i<newvars; i++) {
        v = dataset->v - newvars + i;
        if (list_show_var(sr, v, 0, NULL)) {
            list_append_var_simple(store, &iter, v);
        }
    }
}

static void new_var_callback (GtkWidget *w, selector *sr)
{
    edit_dialog(GENR, _("gretl: add var"),
                _("Enter formula for new variable"),
                NULL, do_selector_genr, sr,
                VARCLICK_INSERT_NAME, sr->dlg);
}

static GtkWidget *add_var_button (selector *sr)
{
    GtkWidget *img = gtk_image_new_from_stock(GTK_STOCK_ADD,
                                              GTK_ICON_SIZE_MENU);
    GtkWidget *button = gtk_button_new();

    gtk_container_add(GTK_CONTAINER(button), img);
    gretl_tooltips_add(GTK_WIDGET(button), _("New variable"));
    g_signal_connect(button, "clicked", G_CALLBACK(new_var_callback), sr);

    return button;
}

static void selection_dialog_add_top_label (selector *sr)
{
    GtkWidget *label;
    gchar *s = NULL;
    int ci = sr->ci;

    if (MODEL_CODE(ci) || VEC_CODE(ci) || NONPARAM_CODE(ci)) {
        GtkWidget *hbox, *button;

        hbox = gtk_hbox_new(FALSE, 0);
        button = add_var_button(sr);
        gtk_box_pack_start(GTK_BOX(hbox), button, FALSE, FALSE, 0);
        s = estimator_label(ci);
        if (s != NULL) {
            label = gtk_label_new(_(s));
            gtk_box_pack_start(GTK_BOX(hbox), label, TRUE, FALSE, 0);
            label = gtk_label_new("");
            gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 16);
        }
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 0);
    } else if (ci == SAVE_FUNCTIONS) {
        GtkWidget *hbox, *entry;

        hbox = gtk_hbox_new(FALSE, 5);
        label = gtk_label_new(_("Name for package:"));
        gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
        entry = gtk_entry_new();
        gtk_entry_set_max_length(GTK_ENTRY(entry), 31);
        gtk_entry_set_width_chars(GTK_ENTRY(entry), 24);
        gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);
        gtk_box_pack_start(GTK_BOX(sr->vbox), hbox, FALSE, FALSE, 5);
        gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);
        sr->extra[0] = entry;
    } else {
        if (ci == GR_XY)
            s = N_("XY scatterplot");
        else if (ci == GR_IMP)
            s = N_("plot with impulses");
        else if (ci == GR_3D)
            s = N_("3D plot");
        else if (ci == SCATTERS)
            s = N_("multiple scatterplots");
        else if (ci == GR_DUMMY)
            s = N_("factorized plot");
        else if (ci == GR_FBOX)
            s = N_("factorized boxplot");
        else if (ci == FSUMMARY)
            s = N_("factorized statistics");
        else if (ci == GR_XYZ)
            s = N_("scatterplot with control");
        else if (ci == ANOVA)
            s = N_("ANOVA");

        if (s != NULL) {
            label = gtk_label_new(_(s));
            gtk_box_pack_start(GTK_BOX(sr->vbox), label, FALSE, FALSE, 5);
        }
    }
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

static void primary_rhs_varlist (selector *sr)
{
    GtkTreeModel *mod;
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *hbox;
    GtkWidget *tmp = NULL;
    GtkWidget **lptr = NULL;
    int i;

    if (COINT_CODE(sr->ci)) {
        tmp = gtk_label_new(_("Variables to test"));
    } else if (VEC_CODE(sr->ci)) {
        tmp = gtk_label_new(_("Endogenous variables"));
    } else if (sr->ci == BIPROBIT) {
        tmp = gtk_label_new(_("Equation 1 regressors"));
    } else if (MODEL_CODE(sr->ci)) {
        tmp = gtk_label_new(_("Regressors"));
    } else if (sr->ci == GR_XY || sr->ci == GR_IMP) {
        tmp = gtk_label_new(_("Y-axis variables"));
    } else if (sr->ci == SCATTERS) {
        multiplot_label = tmp = gtk_label_new(_("X-axis variables"));
    } else if (FNPKG_CODE(sr->ci)) {
        tmp = gtk_label_new(_("Public functions"));
    }

    if (tmp != NULL) {
        table_add_right(sr, tmp, 1);
    }

    hbox = gtk_hbox_new(FALSE, 5);

    /* add the actual primary RHS listbox */
    sr->rvars1 = var_list_box_new(GTK_BOX(hbox), sr, SR_RVARS1);
    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));

    if (ARMA_RELATED(sr->ci)) {
        lptr = &sr->lags_button;
    } else if (sr->ci == VAR) { /* FIXME VECM? */
        lptr = &sr->extra[EXTRA_LAGS];
    }

    /* add push/pull buttons */
    push_pull_buttons(sr, add_to_rvars1_callback,
                      remove_from_rvars1_callback,
                      lptr, OPT_NONE);

    store = GTK_LIST_STORE(mod);
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(mod, &iter);

    if (MODEL_CODE(sr->ci)) {
        if (ARMA_RELATED(sr->ci) || sr->ci == GARCH ||
	    sr->ci == REGLS || NONPARAM_CODE(sr->ci)) {
            ; /* skip */
        } else if (xlist == NULL || has_0(xlist)) {
            /* stick the constant in by default */
            list_append_var(mod, &iter, 0, sr, SR_RVARS1);
        }
        if (xlist != NULL) {
            /* we have a saved list of regressors */
            int nx = 0;

            for (i=1; i<=xlist[0]; i++) {
                if (xlist[i] != 0) {
                    list_append_var(mod, &iter, xlist[i], sr, SR_RVARS1);
                    nx++;
                }
            }
            if (nx > 0 && sr->lags_button != NULL && sr->ci == ARMA) {
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

    table_add_right(sr, hbox, 0);
}

/* On opening the selector dialog, select (if possible) the
   first relevant item in the listbox on the left.
*/

static void selector_set_focus (selector *sr)
{
    if (sr->lvars != NULL) {
        GtkTreeModel *mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars));
        GtkTreeSelection *sel;
        GtkTreeIter iter;
        gboolean do_sel;
        int v;

        gtk_widget_grab_focus(sr->lvars);
        do_sel = gtk_tree_model_get_iter_first(mod, &iter);
        if (do_sel) {
            if (!FNPKG_CODE(sr->ci) && list_show_var(sr, 0, 0, NULL)) {
                /* don't select the constant: skip a row */
                do_sel = gtk_tree_model_iter_next(mod, &iter);
            }
            while (do_sel && MODEL_CODE(sr->ci)) {
                /* skip named lists? */
                gtk_tree_model_get(mod, &iter, COL_ID, &v, -1);
                if (v > 0) {
                    break;
                }
                do_sel = gtk_tree_model_iter_next(mod, &iter);
            }
        }
        if (do_sel) {
            sel = gtk_tree_view_get_selection(GTK_TREE_VIEW(sr->lvars));
            gtk_tree_selection_select_iter(sel, &iter);
        }
    }
}

static void list_append_named_lists (GtkListStore *store,
                                     GtkTreeIter *iterp)
{
    GList *llist = user_var_names_for_type(GRETL_TYPE_LIST);
    GList *tail = llist;
    const int *list;

    while (tail != NULL) {
        list = get_list_by_name(tail->data);
        if (list != NULL && list[0] > 0 &&
            series_get_midas_period(dataset, list[1]) == 0) {
            gtk_list_store_append(store, iterp);
            gtk_list_store_set(store, iterp, COL_ID, -1, COL_LAG, 0,
                               COL_NAME, tail->data, -1);
        }
        tail = tail->next;
    }

    g_list_free(llist);
}

static int midas_special_left_panel (selector *sr,
                                     GtkWidget *left_box,
                                     int saverow)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *l2box, *lbl;
    int i, m, nmidas = 0;
    int err = 0;

    lbl = gtk_label_new("MIDAS vars");
    l2box = gtk_hbox_new(FALSE, 0);
    sr->lvars2 = var_list_box_new(GTK_BOX(l2box), sr, SR_LVARS2);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars2)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    alt_table_add_left(sr, left_box, 0, saverow);
    alt_table_add_left(sr, lbl, saverow, saverow + 1);
    alt_table_add_left(sr, l2box, saverow + 1, sr->n_rows);

    for (i=1; i<dataset->v; i++) {
        m = series_is_midas_anchor(dataset, i);
        if (m > 0 && i + m <= dataset->v) {
            int is_midas = 1;
            int j, p, p0 = m;

            for (j=i+1; j<i+m; j++) {
                p = series_get_midas_period(dataset, j);
                if (p != p0 - 1) {
                    is_midas = 0;
                    break;
                } else {
                    p0 = p;
                }
            }
            if (is_midas) {
                nmidas++;
                list_append_midas_var(store, &iter, i, m);
            }
        }
    }

    if (nmidas == 0) {
        /* "can't happen" */
        err = 1;
    }

    return err;
}

selector *selection_dialog (int ci, const char *title, void *data,
                            sr_callback callback)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *left_box;
    selector *sr;
    int preselect = 0;
    int saverow;
    int i, yvar = 0;

    if (presel > 0) {
	preselect = presel;
	presel = 0;
    } else if (dataset->v == 2 && NO_X_OK(ci)) {
	preselect = 1;
    }

    if (open_selector != NULL) {
        gtk_window_present(GTK_WINDOW(open_selector->dlg));
        return open_selector;
    }

    sr = mymalloc(sizeof *sr);

    if (sr == NULL) {
        return NULL;
    }

    selector_init(sr, ci, title, callback, NULL, data, SELECTOR_FULL);
    selection_dialog_add_top_label(sr);

    /* the following encloses LHS lvars, depvar and indepvar stuff */
    sr->table = gtk_table_new(sr->n_rows, 3, FALSE);

    /* LHS: list of elements to choose from */
    left_box = gtk_hbox_new(FALSE, 0);
    sr->lvars = var_list_box_new(GTK_BOX(left_box), sr, SR_LVARS);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    if (FNPKG_CODE(ci)) {
        available_functions_list(sr);
    } else {
        for (i=0; i<dataset->v; i++) {
            if (i == 1 && (MODEL_CODE(ci) || VEC_CODE(ci))) {
                list_append_named_lists(store, &iter);
            }
            if (list_show_var(sr, i, 0, NULL)) {
                list_append_var_simple(store, &iter, i);
            }
        }
    }

    if (MODEL_CODE(ci) || NONPARAM_CODE(ci) || ci == ANOVA) {
        /* models: top right -> dependent variable */
        yvar = build_depvar_section(sr, preselect);
    } else if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY ||
               ci == SCATTERS || ci == GR_3D || ci == GR_XYZ ||
               ci == GR_FBOX || ci == FSUMMARY) {
        /* top right -> x-axis variable or equivalent */
        build_x_axis_section(sr, preselect);
    } else if (FNPKG_CODE(ci)) {
        primary_rhs_varlist(sr);
    }

    /* middle right: used for some estimators and factored plot */
    if (ci == WLS || ci == AR || ci == ARCH || USE_ZLIST(ci) ||
        VEC_CODE(ci) || ci == COUNTMOD || ci == DURATION ||
        ci == QUANTREG || ci == INTREG || ci == TOBIT ||
        ci == DPANEL || ci == MIDASREG || THREE_VARS_CODE(ci) ||
        NONPARAM_CODE(ci)) {
        build_mid_section(sr);
    }

    saverow = sr->row;

    if (ci == GR_FBOX || ci == FSUMMARY || THREE_VARS_CODE(ci)) {
        /* choose extra var for plot or factorized stats */
        extra_plotvar_box(sr);
    } else if (AUX_LAST(ci)) {
        secondary_rhs_varlist(sr);
    } else if (!NONPARAM_CODE(ci)) {
        /* all other uses: list of vars */
        primary_rhs_varlist(sr);
    }

    if (ci == MIDASREG) {
        /* we need two left-hand list boxes, the second to
           hold high-frequency vars */
        midas_special_left_panel(sr, left_box, saverow);
    } else {
        /* add left-hand column now we know how many rows it should span */
        table_add_left(sr, left_box, 0, sr->n_rows);
    }

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), sr->table, TRUE, TRUE, 0);

    if (ARMA_RELATED(ci)) {
        /* AR, D, MA for ARIMA */
        build_arma_spins(sr);
    } else if (ci == GARCH) {
        /* P and Q for GARCH */
        build_garch_spins(sr);
    }

    /* toggle switches for some cases */
    if (WANT_TOGGLES(ci)) {
        build_selector_switches(sr);
    }

#ifdef GNUPLOT3D
    if (ci == GR_3D) {
        build_selector_switches(sr);
    }
#endif

    /* radio buttons for some */
    if (want_radios(sr)) {
        build_selector_radios(sr);
    }

    /* drop-down selector for some */
    if (want_combo(sr)) {
        build_selector_combo(sr);
    }

    /* plus lag selection stuff, if relevant */
    if (dataset_lags_ok(dataset)) {
        if (ci == GR_XY || ci == GR_IMP || ci == GR_DUMMY || \
            ci == SCATTERS || ci == GR_XYZ) {
            unhide_lags_switch(sr);
        }
        if (MODEL_CODE(ci) && !ARMA_RELATED(ci)) {
            lag_selector_button(sr);
        }
        if (select_lags_depvar(ci) && yvar > 0) {
            maybe_activate_depvar_lags(sr->depvar, sr);
            maybe_insert_depvar_lags(sr, yvar, 0);
        }
    }

    /* buttons: Help, Clear, Cancel, OK */
    build_selector_buttons(sr);

    gtk_widget_show_all(sr->dlg);
    selector_set_focus(sr);

    return sr;
}

static char *simple_sel_label (int ci)
{
    switch (ci) {
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
    case QQPLOT:
        return N_("Q-Q plot: select one or two variables");
    case ELLIPSE:
        return N_("Confidence region: select two variables");
    case PRINT:
        return N_("Select variables to display");
    case GR_BOX:
        return N_("Select variables for boxplot");
    case GR_PLOT:
    case TSPLOTS:
        return N_("Select variables to plot");
    case SAVE_DATA:
    case SAVE_DATA_AS:
    case EXPORT_CSV:
    case EXPORT_R:
    case EXPORT_OCTAVE:
    case EXPORT_JM:
    case EXPORT_DAT:
    case EXPORT:
        return N_("Select variables to save");
    case COPY_CSV:
        return N_("Select variables to copy");
    case CHOW:
        return N_("Select variables to test");
    case DEFINE_LIST:
        return N_("Define named list");
    case REGLS_PLOTSEL:
	return N_("Select coefficients to track (max 25)");
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
            /* note: special case, not using varnames as such */
            gretl_model_get_param_name(pmod, dataset, i, pname);
            gtk_list_store_append(store, &iter);
            gtk_list_store_set(store, &iter,
                               COL_ID, i, COL_LAG, 0,
                               COL_NAME, pname,
                               -1);
            nvars++;
        }
        g_object_set_data(G_OBJECT(sr->lvars), "keep-names",
                          GINT_TO_POINTER(1));
    } else if (sr->ci == OMIT || sr->ci == ADD || sr->ci == CHOW ||
               sr->ci == COEFFSUM) {
        int *xlist = gretl_model_get_x_list(pmod);

        if (xlist == NULL) {
            return 0;
        }

        if (sr->ci == ADD) {
            int dv = gretl_model_get_depvar(pmod);

            for (i=0; i<dataset->v; i++) {
                if (!in_gretl_list(xlist, i) && i != dv &&
                    !series_is_hidden(dataset, i)) {
                    gtk_list_store_append(store, &iter);
                    gtk_list_store_set(store, &iter,
                                       COL_ID, i, COL_LAG, 0,
                                       COL_NAME, dataset->varname[i],
                                       -1);
                    nvars++;
                }
            }
        } else {
            for (i=1; i<=xlist[0]; i++) {
                if (sr->ci == CHOW && xlist[i] == 0) {
                    continue;
                }
                gtk_list_store_append(store, &iter);
                gtk_list_store_set(store, &iter,
                                   COL_ID, xlist[i], COL_LAG, 0,
                                   COL_NAME, dataset->varname[xlist[i]],
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
                                   COL_ID, xlist[i], COL_LAG, 0,
                                   COL_NAME, dataset->varname[xlist[i]],
                                   -1);
                nvars++;
            }
        }
    }

    return nvars;
}

static void available_functions_list (selector *sr)
{
    GtkListStore *store;
    GtkTreeIter iter;
    const char *fnname;
    int idx;

    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    function_names_init();

    while ((fnname = next_available_function_name(sr->data, &idx)) != NULL) {
        gtk_list_store_append(store, &iter);
        gtk_list_store_set(store, &iter,
                           COL_ID, idx, COL_LAG, 0,
                           COL_NAME, fnname, -1);
    }

    g_object_set_data(G_OBJECT(sr->lvars), "keep-names",
                      GINT_TO_POINTER(1));
}

static GtkWidget *simple_selection_top_label (int ci, const char *title)
{
    const char *s = simple_sel_label(ci);
    GtkWidget *label = NULL;
    GtkWidget *hbox = NULL;

    if (s != NULL && *s != '\0') {
        label = gtk_label_new(_(s));
    } else if (title != NULL) {
        if (!strncmp(title, "gretl: ", 7)) {
            title += 7;
        }
        label = gtk_label_new(title);
    }

    if (label != NULL) {
        hbox = gtk_hbox_new(FALSE, 5);
        gtk_container_add(GTK_CONTAINER(hbox), label);
    }

    return hbox;
}

static void
add_to_rvars1_from_named_list (selector *sr, const int *list,
                               int clear)
{
    GtkTreeModel *mod;
    GtkTreeIter iter;
    int i, v;

    mod = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (mod == NULL) {
        return;
    }

    gtk_tree_model_get_iter_first(mod, &iter);

    if (clear) {
        gtk_list_store_clear(GTK_LIST_STORE(mod));
    }

    for (i=1; i<=list[0]; i++) {
        v = list[i];
        gtk_list_store_append(GTK_LIST_STORE(mod), &iter);
        gtk_list_store_set(GTK_LIST_STORE(mod), &iter,
                           COL_ID, v, COL_LAG, 0,
                           COL_NAME, dataset->varname[v], -1);
    }
}

static void maybe_set_plot_vars (selector *sr)
{
    int mc = mdata_selection_count();

    if (mc > 1 && mc < 7) {
        set_vars_from_main(sr);
    }
}

/* the following is called on start-up when defining a list */

static void maybe_set_listdef_vars (selector *sr)
{
    const char *lname = selector_entry_text(sr);

    if (lname != NULL && *lname != 0) {
        int *list = get_list_by_name(lname);

        if (list != NULL) {
            add_to_rvars1_from_named_list(sr, list, 0);
        }
    } else if (sr->data == NULL) {
        /* called from main window */
        int mc = mdata_selection_count();

        if (mc > 1) {
            set_vars_from_main(sr);
        }
    }
}

/* callback from change in combo selector, when defining or editing
   a list */

static void listdef_vars_callback (GtkComboBox *b, selector *sr)
{
    const char *lname = selector_entry_text(sr);

    if (lname != NULL && *lname != '\0') {
        int *list = get_list_by_name(lname);

        if (list != NULL) {
            add_to_rvars1_from_named_list(sr, list, 1);
        }
    }
}

static void set_name_from_fn_param (GtkWidget *w, selector *sr,
                                    GtkWidget *entry)
{
    gchar *tmp = NULL;

    get_fncall_param_info(sr->parent, NULL, &tmp);

    if (tmp == NULL) {
        int argnum = widget_get_int(w, "argnum");

        tmp = g_strdup_printf("arg%d", argnum + 1);
    }

    gtk_entry_set_text(GTK_ENTRY(entry), tmp);
    g_free(tmp);
}

static void selector_add_list_name_entry (selector *sr)
{
    const char *lname = NULL;
    GtkWidget *src = NULL;
    GList *lnames = NULL;
    GtkWidget *combo = NULL;
    GtkWidget *label;
    GtkWidget *entry;

    if (sr->data != NULL) {
        /* called from function call dialog */
        src = GTK_WIDGET(sr->data);
        lname = gtk_entry_get_text(GTK_ENTRY(src));
        if (current_series_index(dataset, lname) > 0) {
            /* @lname is actually the name of a series */
            lname = NULL;
        }
    } else {
        lnames = user_var_names_for_type(GRETL_TYPE_LIST);
    }

    label = gtk_label_new(_("Name of list"));
    table_add_right(sr, label, 1);

    if (lnames != NULL) {
        combo = combo_box_text_new_with_entry();
        set_combo_box_strings_from_list(combo, lnames);
        entry = gtk_bin_get_child(GTK_BIN(combo));
        g_signal_connect(G_OBJECT(combo), "changed",
                         G_CALLBACK(listdef_vars_callback), sr);
    } else {
        entry = gtk_entry_new();
    }

    sr->extra[0] = entry;

    gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN);
    gtk_entry_set_activates_default(GTK_ENTRY(entry), TRUE);

    if (src != NULL) {
        if (lname != NULL && *lname != '\0' && strcmp(lname, "null")) {
            gtk_entry_set_text(GTK_ENTRY(entry), lname);
        } else {
            set_name_from_fn_param(src, sr, entry);
        }
        gtk_editable_select_region(GTK_EDITABLE(entry), 0, -1);
    }

    if (combo != NULL) {
        table_add_right(sr, combo, 1);
    } else {
        table_add_right(sr, entry, 1);
    }

    table_add_vwedge(sr);
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

        if (sr->ci == CORR) {
            set_n_rvars1(sr, list[0]);
        }
    }

    free(list);
}

selector *
simple_selection_with_data (int ci, const char *title,
                            sr_callback callback,
                            GtkWidget *parent, gpointer data)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *left_box, *right_box;
    GtkWidget *tmp;
    selector *sr;
    gchar **exclude = NULL;
    int nleft = 0;
    int i, err = 0;

    if (open_selector != NULL) {
        gtk_window_present(GTK_WINDOW(open_selector->dlg));
        return open_selector;
    }

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) {
        return NULL;
    }

    selector_init(sr, ci, title, callback, parent, data, SELECTOR_SIMPLE);

    tmp = simple_selection_top_label(ci, title);
    if (tmp != NULL) {
        gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
    }

    /* tables to hold vboxes and buttons */
    sr->table = gtk_table_new(1, 3, FALSE);

    /* holds list of elements available for selection */
    left_box = gtk_vbox_new(FALSE, 5);

    sr->lvars = var_list_box_new(GTK_BOX(left_box), sr, SR_LVARS);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* check for any specific exclusions */
    if (ci == DEFINE_LIST && data != NULL) {
	int listmin = 1;

	exclude = get_listdef_exclude(data, &listmin);
	if (listmin > 1) {
	    sr->flags |= SR_NO_SINGLE;
	}
    }

    /* add selectable series */
    if (ci == OMIT || ci == ADD || ci == COEFFSUM ||
        ci == ELLIPSE || ci == VAROMIT || ci == CHOW) {
        nleft = add_omit_list(data, sr);
    } else {
        int start = (ci == DEFINE_LIST || ci == DEFINE_MATRIX)? 0 : 1;

        for (i=start; i<dataset->v; i++) {
            if (i == 1 && SHOW_LISTS_CODE(ci)) {
                list_append_named_lists(store, &iter);
            }
            if (list_show_var(sr, i, 0, exclude)) {
                list_append_var_simple(store, &iter, i);
                nleft++;
            }
        }
    }

    sr->n_left = nleft;

    right_box = gtk_vbox_new(FALSE, 5);
    sr->rvars1 = var_list_box_new(GTK_BOX(right_box), sr, SR_RVARS1);

    if (sr->ci == CORR) {
        /* ensure that RHS var count is set to zero */
        set_n_rvars1(sr, 0);
    }

    /* pre-fill RHS box? Only if we have 2 or more vars selected in the
       main window and if the command is "suitable"
    */
    if (RHS_PREFILL(ci)) {
        maybe_prefill_RHS(sr);
    }

    /* entry field for some uses */
    if (ci == DEFINE_LIST) {
        selector_add_list_name_entry(sr);
    }

    /* put buttons into mid-section */
    push_pull_buttons(sr, add_to_rvars1_callback,
                      remove_from_rvars1_callback,
                      NULL, OPT_A | OPT_R);

    /* pack RHS */
    table_add_right(sr, right_box, 0);

    /* pack left-hand stuff */
    table_add_left(sr, left_box, 0, sr->n_rows);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), sr->table, TRUE, TRUE, 0);

    /* "unhide lags" check box? */
    if ((sr->ci == DEFINE_LIST || sr->ci == EXPORT || SAVE_DATA_ACTION(sr->ci))
        && lags_hidden(sr)) {
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
    /* combo/dropdown list? */
    if (want_combo(sr)) {
        build_selector_combo(sr);
    }
    if (ci == ELLIPSE) {
        build_ellipse_spin(sr);
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
    } else if (sr->ci == TSPLOTS || sr->ci == GR_BOX) {
        maybe_set_plot_vars(sr);
    }

    if (nleft == 0) {
        err = E_DATA;
    } else if ((ci == COEFFSUM || ci == ELLIPSE) && nleft < 2) {
        err = E_DATA;
    } else {
        gtk_widget_show_all(sr->dlg);
    }

    if (err) {
        warnbox(_("No suitable data are available"));
        gtk_widget_destroy(sr->dlg);
        sr = NULL;
    } else if (sr->ci == DEFINE_MATRIX) {
        selector_set_blocking(sr, 1);
    }

    if (exclude != NULL) {
	g_strfreev(exclude);
    }

    return sr;
}

selector *
sublist_selection (int ci, const char *title,
                   sr_callback callback,
		   GtkWidget *parent, const int *list,
		   const int *presel, void *data)
{
    GtkListStore *store;
    GtkTreeIter iter;
    GtkWidget *left_box, *right_box;
    GtkWidget *tmp;
    selector *sr;
    int selmax = 0;
    int i;

    sr = mymalloc(sizeof *sr);
    if (sr == NULL) {
        return NULL;
    }

    selector_init(sr, ci, title, callback, parent, data, SELECTOR_SIMPLE);

    tmp = simple_selection_top_label(ci, title);
    if (tmp != NULL) {
        gtk_box_pack_start(GTK_BOX(sr->vbox), tmp, FALSE, FALSE, 0);
    }

    sr->table = gtk_table_new(1, 3, FALSE);
    left_box = gtk_vbox_new(FALSE, 5);
    sr->lvars = var_list_box_new(GTK_BOX(left_box), sr, SR_LVARS);
    store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->lvars)));
    gtk_list_store_clear(store);
    gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);

    /* add @list to left-hand box */
    for (i=1; i<=list[0]; i++) {
	list_append_var_simple(store, &iter, list[i]);
	sr->n_left += 1;
    }

    /* specify right-hand_box */
    right_box = gtk_vbox_new(FALSE, 5);
    sr->rvars1 = var_list_box_new(GTK_BOX(right_box), sr, SR_RVARS1);
    if (ci == REGLS_PLOTSEL) {
	selmax = 25;
	widget_set_int(sr->rvars1, "selmax", selmax);
    }
    if (presel != NULL) {
	/* put some content into the right-hand box */
	int n_right = 0;

	store = GTK_LIST_STORE(gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1)));
	gtk_list_store_clear(store);
	gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
	for (i=1; i<=presel[0]; i++) {
	    list_append_var_simple(store, &iter, presel[i]);
	    if (++n_right == selmax) {
		break;
	    }
	}
    }

    /* put buttons into mid-section */
    push_pull_buttons(sr, add_to_rvars1_callback,
                      remove_from_rvars1_callback,
                      NULL, OPT_A | OPT_R);

    /* pack RHS and LHS */
    table_add_right(sr, right_box, 0);
    table_add_left(sr, left_box, 0, sr->n_rows);

    /* pack the whole central section into the dialog's vbox */
    gtk_box_pack_start(GTK_BOX(sr->vbox), sr->table, TRUE, TRUE, 0);

    /* buttons: Help, Clear, Cancel, OK */
    build_selector_buttons(sr);

    gtk_widget_show_all(sr->dlg);

    return sr;
}

selector *simple_selection (int ci, const char *title,
                            sr_callback callback,
                            GtkWidget *parent)
{
    return simple_selection_with_data(ci, title, callback, parent, NULL);
}


static void restore_vwin_menu (GtkWidget *w, windata_t *vwin)
{
    if (vwin != NULL && vwin->mbar != NULL) {
        gtk_widget_set_sensitive(vwin->mbar, TRUE);
    }
}

selector *
simple_selection_for_viewer (int ci, const char *title,
                             sr_callback callback,
                             windata_t *vwin)
{
    selector *sr;

    sr = simple_selection_with_data(ci, title, callback,
                                    vwin_toplevel(vwin),
                                    vwin);

    if (sr != NULL && vwin->mbar != NULL) {
        if (window_is_tab(vwin)) {
            tabwin_register_dialog(sr->dlg, vwin_toplevel(vwin));
        } else {
            gtk_widget_set_sensitive(vwin->mbar, FALSE);
            g_signal_connect(sr->dlg, "destroy",
                             G_CALLBACK(restore_vwin_menu),
                             vwin);
        }
    }

    return sr;
}

static gchar *get_or_set_storelist (const char *s, int reset)
{
    static gchar *storelist;

    if (s != NULL) {
        /* set storelist */
        g_free(storelist);
        if (*s == '\0') {
            storelist = NULL;
        } else {
            storelist = g_strdup(s);
        }
        return NULL;
    } else if (reset) {
        g_free(storelist);
        storelist = NULL;
        return NULL;
    } else {
        /* retrieve storelist (and NULL-ify it) */
        gchar *ret = storelist;

        storelist = NULL;
        return ret;
    }
}

gchar *get_selector_storelist (void)
{
    return get_or_set_storelist(NULL, 0);
}

void set_selector_storelist (const char *s)
{
    if (s == NULL) {
        /* forces a reset */
        get_or_set_storelist(NULL, 1);
    } else {
        get_or_set_storelist(s, 0);
    }
}

static int get_selected_function_names (selector *sr,
                                        char ***S1, int *n1,
                                        char ***S2, int *n2)
{
    GtkTreeModel *src;
    GtkTreeIter iter;
    gchar *fname;

    /* public function names are in top right list */
    src = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
    if (gtk_tree_model_get_iter_first(src, &iter)) {
        do {
            gtk_tree_model_get(src, &iter, COL_NAME, &fname, -1);
            strings_array_add(S1, n1, gretl_strdup(fname));
            g_free(fname);
        } while (gtk_tree_model_iter_next(src, &iter));
    }

    if (*n1 == 0) {
        warnbox(_("You must specify a public interface"));
        return 1;
    }

    /* private function names are in lower right list */
    src = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
    gtk_tree_model_get_iter_first(src, &iter);
    if (gtk_tree_model_get_iter_first(src, &iter)) {
        do {
            gtk_tree_model_get(src, &iter, COL_NAME, &fname, -1);
            strings_array_add(S2, n2, gretl_strdup(fname));
            g_free(fname);
        } while (gtk_tree_model_iter_next(src, &iter));
    }

    return 0;
}

static int pkg_add_remove_callback (selector *sr)
{
    void *p = sr->data;
    char **pubnames = NULL;
    char **privnames = NULL;
    int npub = 0;
    int npriv = 0;
    int err;

    err = get_selected_function_names(sr, &pubnames, &npub,
                                      &privnames, &npriv);
    if (!err) {
        revise_function_package(p, pubnames, npub, privnames, npriv);
    }

    return err;
}

static int check_pkgname (const char *name,
                          GtkWidget *parent)
{
    int n = strlen(name);
    int err = 0;

    if (has_suffix(name, ".gfn")) {
        n -= 4;
    }

    if (n >= FN_NAMELEN) {
        /* too long */
        err = 1;
    } else if (gretl_namechar_spn(name) != n) {
        /* contains funny stuff */
        err = 1;
    }

    if (err) {
        msgbox(_("Invalid package name: the name must start with a letter,\n"
                 "must be less than 32 characters in length, and must include\n"
                 "only ASCII letters, numbers and '_'."),
               GTK_MESSAGE_ERROR, parent);
    }

    return err;
}

static int functions_selected_callback (selector *sr)
{
    GtkWidget *entry = sr->extra[0];
    char **pubnames = NULL;
    char **privnames = NULL;
    int npub = 0;
    int npriv = 0;
    const char *s;
    int err = 0;

    s = gtk_entry_get_text(GTK_ENTRY(entry));

    if (s == NULL || *s == '\0' || check_pkgname(s, sr->dlg)) {
        gtk_widget_grab_focus(entry);
        err = 1;
    } else {
        err = get_selected_function_names(sr, &pubnames, &npub,
                                          &privnames, &npriv);
    }

    if (!err) {
        gchar *pkgname = g_strdup(s);

        if (has_suffix(pkgname, ".gfn")) {
            char *p = strrchr(pkgname, '.');

            *p = '\0';
        }
        edit_new_function_package(pkgname, pubnames, npub,
                                  privnames, npriv);
    }

    return err;
}

static int data_export_selection_callback (selector *sr)
{
    int ci = sr->ci;

    if ((sr->cmdlist == NULL || *sr->cmdlist == '\0') && sr->n_left == 0) {
        warnbox(_("No variables are selected"));
        /* return non-zero to block closing of dialog */
        return 1;
    }

    if (ci == EXPORT) {
        /* set the specific export format based on the
           option from the series selector
        */
        if (sr->opts & OPT_C) {
            ci = EXPORT_CSV;
        } else if (sr->opts & OPT_R) {
            ci = EXPORT_R;
        } else if (sr->opts & OPT_M) {
            ci = EXPORT_OCTAVE;
        } else if (sr->opts & OPT_J) {
            ci = EXPORT_JM;
        } else if (sr->opts & OPT_G) {
            ci = EXPORT_DAT;
        } else if (sr->opts & OPT_S) {
            ci = EXPORT_DTA;
        } else if (sr->opts & OPT_D) {
            ci = EXPORT_DB;
        } else if (sr->opts & OPT_P) {
            ci = SAVE_MAP;
        } else if (sr->opts & OPT_B) {
            ci = EXPORT_GDTB;
        } else {
            ci = EXPORT_GDT;
        }
    }

    if (sr->cmdlist != NULL) {
        set_selector_storelist(sr->cmdlist);
    }

    gtk_widget_destroy(sr->dlg);

    if (ci == EXPORT_CSV) {
        int resp = csv_options_dialog(EXPORT_CSV, GRETL_OBJ_DSET,
                                      NULL);

        if (canceled(resp)) {
            set_selector_storelist(NULL);
            return 0;
        }
    }

    if (ci != COPY_CSV) {
        file_selector(ci, FSEL_DATA_NONE, NULL);
    }

    return 0;
}

void data_export_selection_wrapper (int file_ci)
{
    selector *sr;

    set_selector_storelist(NULL);

    sr = simple_selection(file_ci, (file_ci == COPY_CSV)?
                          _("Copy data") : _("Export data"),
                          data_export_selection_callback,
                          NULL);
    if (sr != NULL) {
        selector_set_blocking(sr, 0);
    }
}

void functions_selection_wrapper (GtkWidget *parent)
{
    int err;

    set_selector_storelist(NULL);
    err = no_user_functions_check(parent);

    if (!err) {
        selector *sr;

        sr = selection_dialog(SAVE_FUNCTIONS, _("Select functions"),
                              NULL, functions_selected_callback);
        if (sr != NULL) {
            selector_set_blocking(sr, 1);
        }
    }
}

void add_remove_functions_dialog (char **pubnames, int npub,
                                  char **privnames, int npriv,
                                  void *p1, void *p2)
{
    fnpkg *pkg = p1;
    void *finfo = p2;
    const char *title = NULL;
    selector *sr = NULL;

    if (pkg != NULL) {
        title = function_package_get_name(pkg);
    }

    if (title == NULL || *title == '\0') {
        title = _("Edit package list");
    }

    /* to start with, set @pkg as sr->data to enable the correct
       LHS listing of functions for the package */

    sr = selection_dialog(EDIT_FUNCTIONS, title,
                          pkg, pkg_add_remove_callback);

    if (sr != NULL) {
        GtkTreeModel *model;
        GtkListStore *store;
        GtkTreeIter iter;
        int i, idx;

        /* switch sr->data to point to the current 'editor' */
        sr->data = finfo;

        /* put current @pubnames into top right list box */
        if (npub > 0) {
            model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars1));
            store = GTK_LIST_STORE(model);
            gtk_list_store_clear(store);
            gtk_tree_model_get_iter_first(model, &iter);
            for (i=0; i<npub; i++) {
                idx = user_function_index_by_name(pubnames[i], pkg);
                if (idx >= 0) {
                    gtk_list_store_append(store, &iter);
                    gtk_list_store_set(store, &iter,
                                       COL_ID, idx, COL_LAG, 0,
                                       COL_NAME, pubnames[i], -1);
                }
            }
        }

        /* put current @privnames into lower right box */
        if (npriv > 0) {
            model = gtk_tree_view_get_model(GTK_TREE_VIEW(sr->rvars2));
            store = GTK_LIST_STORE(model);
            gtk_list_store_clear(store);
            gtk_tree_model_get_iter_first(model, &iter);
            for (i=0; i<npriv; i++) {
                idx = user_function_index_by_name(privnames[i], pkg);
                if (idx >= 0) {
                    gtk_list_store_append(store, &iter);
                    gtk_list_store_set(store, &iter,
                                       COL_ID, idx, COL_LAG, 0,
                                       COL_NAME, privnames[i], -1);
                }
            }
        }

        selector_set_blocking(sr, 1);
    }
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

gpointer selector_get_extra_data (const selector *sr)
{
    return sr->extra_data;
}

gretlopt selector_get_opts (const selector *sr)
{
    gretlopt ret;

    if (sr->ci == PANEL_B) {
        ret = sr->opts | OPT_B;
    } else if (sr->ci == FE_LOGISTIC) {
        ret = sr->opts | OPT_F;
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

void selector_cleanup (void)
{
    clear_selector();

    if (open_selector != NULL) {
        gtk_widget_destroy(open_selector->dlg);
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
    GtkWidget *spin1;   /* spin button for minimum lag */
    GtkWidget *spin2;   /* spin button for maximum lag */
    GtkWidget *entry;   /* text entry for specific lags */
    GtkWidget *toggle;  /* button to switch between spin buttons and entry */
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

/* retrieve the min or max lag from a spin button and process
   the result */

static void lag_set_callback (GtkWidget *w, gpointer p)
{
    var_lag_info *vlinfo;
    int lag, *plag;

    vlinfo = (var_lag_info *) g_object_get_data(G_OBJECT(w), "vlinfo");

    plag = (w == vlinfo->spin1)? &vlinfo->lmin : &vlinfo->lmax;
    lag = *plag = spin_get_int(w);

    /* force consistency if need be */

    if (w == vlinfo->spin1 && vlinfo->spin2 != NULL) {
        if (spin_get_int(vlinfo->spin2) < lag) {
            gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlinfo->spin2), lag);
        }
    } else if (w == vlinfo->spin2 && vlinfo->spin1 != NULL) {
        if (spin_get_int(vlinfo->spin1) > lag) {
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

/* switch from using min/max spins to using the "specific lags"
   text entry, or vice versa */

static void activate_specific_lags (GtkWidget *w, var_lag_info *vlinfo)
{
    gboolean active = button_is_active(w);

    if (active) {
        gtk_widget_set_sensitive(vlinfo->entry, TRUE);
        gtk_widget_set_sensitive(vlinfo->spin1, FALSE);
        gtk_widget_set_sensitive(vlinfo->spin2, FALSE);
        gtk_widget_grab_focus(vlinfo->entry);
    } else {
        gtk_widget_set_sensitive(vlinfo->entry, FALSE);
        gtk_widget_set_sensitive(vlinfo->spin1, TRUE);
        gtk_widget_set_sensitive(vlinfo->spin2, TRUE);
        gtk_widget_grab_focus(vlinfo->spin1);
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
            if (!gtk_widget_is_sensitive(vlset[i].spin1) &&
                !gtk_widget_is_sensitive(vlset[i].entry)) {
                /* dependent var lags are disabled */
                vlset[i].lmin = vlset[i].lmax = NOT_LAG; /* ?? */
                free(vlset[i].lspec);
                vlset[i].lspec = NULL;
                continue;
            }
        }
        active = button_is_active(vlset[i].toggle);
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
    gboolean active = button_is_active(w);

    gtk_widget_set_sensitive(vlinfo->spin1, active);
    gtk_widget_set_sensitive(vlinfo->spin2, active);
    gtk_widget_set_sensitive(vlinfo->toggle, active);

    if (active) {
        vlinfo->lmin = spin_get_int(vlinfo->spin1);
        vlinfo->lmax = spin_get_int(vlinfo->spin2);
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

static void resensitize_selector (GtkWidget *w, gpointer p)
{
    if (open_selector != NULL) {
        gtk_widget_set_sensitive(open_selector->dlg, TRUE);
    }
}

/* The actual lag selection dialog: we provide spins for a lag
   range and also a free-form entry field for non-contiguous lags.  In
   some circumstances we allow specification of lags for the dependent
   variable as well as the independent vars.
*/

static int
lags_dialog (const int *list, var_lag_info *vlinfo, selector *sr)
{
    GtkWidget *lbl, *dialog, *myvbox;
    GtkWidget *tbl, *tmp, *vbox, *hbox;
    GtkWidget *y_check = NULL;
    gint tbl_len;
    double lmax;
    int VAR_special, insts;
    int i, j;
    int ret = GRETL_CANCEL;

    dialog = gretl_dialog_new(_("lag order"), sr->dlg,
                              GRETL_DLG_BLOCK | GRETL_DLG_RESIZE);

    g_signal_connect(G_OBJECT(dialog), "destroy",
                     G_CALLBACK(resensitize_selector), NULL);

    myvbox = gtk_vbox_new(FALSE, 5);

    VAR_special = (vlinfo[0].v == VDEFLT && vlinfo[0].context == LAG_Y_V);
    insts = in_gretl_list(list, LAG_W);

    if (sr->ci == REGLS) {
        lmax = 32;
    } else {
        lmax = 2 * (dataset->t2 - dataset->t1) / list[0];
    }

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
                lbl = gtk_label_new(dataset->varname[li]);
            }
            gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 0, 1, i, i+1);
        }

        vlj = &vlinfo[j++];

        if (depvar_row(vlj->context) || vlj->context == LAG_Y_V) {
            lmin = 1;
        }

        /* min. lag spin button */
        vlj->spin1 = gtk_spin_button_new_with_range(lmin, lmax, 1);
        gtk_table_attach_defaults(GTK_TABLE(tbl), vlj->spin1, 1, 2, i, i+1);
        g_object_set_data(G_OBJECT(vlj->spin1), "vlinfo", vlj);
        lagsel_spin_connect(vlj->spin1);
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(vlj->spin1), vlj->lmin);

        lbl = gtk_label_new(_("to"));
        gtk_table_attach_defaults(GTK_TABLE(tbl), lbl, 2, 3, i, i+1);

        /* max. lag spin button */
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

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));

    if (list[0] > 10) {
        GtkWidget *scroller = gtk_scrolled_window_new(NULL, NULL);

        gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scroller),
                                       GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
        gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scroller),
                                              hbox);
        gtk_box_pack_start(GTK_BOX(vbox), scroller, TRUE, TRUE, 5);
        gtk_widget_set_size_request(scroller, -1, 360);
    } else {
        gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);
    }

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    /* "Cancel" button */
    tmp = cancel_delete_button(hbox, dialog);

    /* "OK" button */
    tmp = ok_validate_button(hbox, &ret, NULL);
    g_signal_connect(G_OBJECT(tmp), "clicked",
                     G_CALLBACK(lag_toggle_register), &vlinfo[0]);
    g_signal_connect_swapped(G_OBJECT(tmp), "clicked",
			     G_CALLBACK(gtk_widget_destroy),
			     dialog);
    gtk_widget_grab_default(tmp);

    /* "Help" button */
    context_help_button(hbox, LAGS_DIALOG);

    gtk_widget_set_sensitive(sr->dlg, FALSE);
    gtk_widget_show_all(dialog);

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
        gtk_tree_model_get(mod, &iter, COL_ID, &modv, -1);
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

            gtk_tree_model_get(model, &iter, COL_ID, &v, COL_LAG, &lag, -1);

            if (v == 0 || !strcmp(dataset->varname[v], "time") ||
                !strcmp(dataset->varname[v], "timesq")) {
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
        sr->ci != DPANEL && sr->ci != MIDASREG) {
        ynum = selector_get_depvar_number(sr);
    }

    listw[0] = (select_lags_primary(sr->ci))? sr->rvars1 : sr->rvars2;

    if (USE_ZLIST(sr->ci)) {
        listw[1] = sr->rvars2;
    }

    if (listw[0] != NULL) {
        /* number of laggable regressors */
        nv[0] = get_laggable_vars(listw[0], LAG_X, NULL, NULL);
        if (nv[0] == 0) {
            listw[0] = NULL;
        }
    }

    if (listw[1] != NULL) {
        /* number of laggable instrumentsm if applicable */
        nv[1] = get_laggable_vars(listw[1], LAG_W, NULL, NULL);
        if (nv[1] == 0) {
            listw[1] = NULL;
        }
    }

#if LDEBUG
    fprintf(stderr, "sr_get_stoch_list: ynum = %d, nv[0] = %d, nv[1] = %d\n",
            ynum, nv[0], nv[1]);
#endif

    if (ynum < 0 && nv[0] == 0 && nv[1] == 0) {
        /* nothing relevant was found */
        errbox(_("Please add some variables to the model first"));
        return NULL;
    }

    /* initialize number of setters and separators */
    nset = nsep = 0;

    /* first pass: figure out how many elements the list should have */

    if (nv[0] > 0) {
        /* regressors */
        *pcontext = LAG_X;
        nset += nv[0] + (nv[0] > 1); /* the Xs plus their defaults */
        if (ynum > 0) {
            nsep++; /* depvar heading */
            nset++; /* depvar row */
        }
    }

    if (nv[1] > 0) {
        /* instruments */
        if (nv[0] > 0) {
            nsep++; /* separator from Xs */
        } else {
            *pcontext = LAG_W;
        }
        nsep++;            /* "instruments" heading */
        nset += nv[1] + (nv[1] > 1); /* instruments plus their defaults */
        if (ynum > 0) {
            nsep++; /* depvar (as instrument) heading */
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

    fprintf(stderr, "nset=%d, nsep=%d\n", nset, nsep);

    slist = gretl_list_new(nset + nsep);
    if (slist == NULL) {
        return NULL;
    }

    /* second pass: actually build the list, inserting special
       list separators if needed
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
        if (nv[j] > 1) {
            /* insert default only if we have more than one entry */
            slist[i++] = VDEFLT;
        }
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
        gretl_charsub(vlinfo->lspec, ',', ' ');
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

/* Keep the global "lag order" spin button in sync with lag
   specifications entered via the the (more complex) lags dialog.  We
   set the value shown in the spin button to the maximum lag; in
   addition we desensitize the global spin button if the user has
   specified gaps in the lag structure.
*/

static void sync_lag_order_spin (var_lag_info *vlinfo,
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

    if (lmax != 0 && lmax != spin_get_int(sr->extra[0])) {
        /* update max lag as shown by spin */
        gtk_spin_button_set_value(GTK_SPIN_BUTTON(sr->extra[0]),
                                  lmax);
    }

    if (lmin > 1 || list != NULL) {
        /* not 1-based consecutive lags: so disable spin */
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
            vlj->lspec = gretl_list_to_numeric_string(laglist);
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

    if (canceled(resp)) {
        y_x_lags_enabled = yxlags;
        y_w_lags_enabled = ywlags;
    } else {
        /* register any changes made by the user */
        for (j=0; j<nvl; j++) {
            int changed = set_lags_for_var(&vlinfo[j], yxlags, ywlags);

            if (vlinfo[j].v != VDEFLT && changed) {
                maybe_revise_var_string(&vlinfo[j], sr);
            }
        }
        if (context == LAG_Y_V) {
            sync_lag_order_spin(&vlinfo[0], sr);
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
