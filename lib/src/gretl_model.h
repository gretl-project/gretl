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

#ifndef GRETL_MODEL_H
#define GRETL_MODEL_H

#include "objstack.h"

#ifdef  __cplusplus
extern "C" {
#endif

typedef enum {
    ARMA_X12A  = 1 << 0, /* using X-12-ARIMA (or X-13) to generate estimates */
    ARMA_EXACT = 1 << 1, /* using exact ML */
    ARMA_LS    = 1 << 2, /* using conditional ML, and O/NLS == CML */
    ARMA_OLS   = 1 << 3  /* OLS == MLE */
} ArmaFlags;

typedef struct CoeffIntervals_ CoeffIntervals;

struct CoeffIntervals_ {
    int asy;
    int ncoeff;
    double alpha;
    double t;
    char **names;
    double *coeff;
    double *maxerr;
    int df;
    gretlopt opt;
};

typedef enum {
    VCV_CLASSICAL,
    VCV_HC,
    VCV_HAC,
    VCV_ML,
    VCV_PANEL,
    VCV_RQ,
    VCV_CLUSTER
} VCVMajorType;

typedef enum {
    ML_UNSET,
    ML_HESSIAN,
    ML_IM,
    ML_OP,
    ML_QML,
    ML_BW,
    ML_HAC,
    ML_VCVMAX
} MLVCVType;

typedef enum {
    KERNEL_BARTLETT,
    KERNEL_PARZEN,
    KERNEL_QS,
    KERNEL_MAX
} HACKernel;

typedef enum {
    PANEL_HAC,  /* clustered by individual (Arellano) */
    PANEL_BK,   /* Beck-Katz PCSE */
    PANEL_TIME, /* clustered by period */
    PANEL_DK,   /* Driscoll-Kraay SCC */
    PANEL_BOTH  /* clustered by both unit and period */
} PanelVCVType;

typedef enum {
    RQ_ASY,   /* asymptotic */
    RQ_NID    /* sandwich */
} RQVCVType;

typedef enum {
    HAC_PREWHITEN = 1
} VCVFlags;

typedef struct VCVInfo_ VCVInfo;

struct VCVInfo_ {
    int vmaj;        /* general type of VCV (see VCVMajorType) */
    int vmin;        /* variant of general type */
    int order;       /* for use with HAC */
    VCVFlags flags;  /* includes prewhitening */
    double bw;       /* for use with QS HAC kernel */
    char *cv1;       /* name of (first) cluster variable */
    char *cv2;       /* name of second cluster var */
};

/* single-equation model commands */

#define MODEL_COMMAND(c) (c == AR || \
                          c == AR1 || \
                          c == ARCH || \
                          c == ARMA || \
			  c == DPANEL ||   \
			  c == DURATION || \
                          c == GARCH || \
                          c == GMM || \
		          c == HECKIT || \
                          c == HSK || \
                          c == INTREG || \
                          c == IVREG || \
                          c == LAD || \
                          c == LOGISTIC || \
                          c == LOGIT || \
			  c == MIDASREG || \
                          c == MLE || \
                          c == MPOLS || \
			  c == NEGBIN || \
                          c == NLS || \
                          c == OLS || \
                          c == PANEL || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == BIPROBIT || \
                          c == QUANTREG || \
                          c == TOBIT || \
                          c == WLS)

#define AR_MODEL(c) (c == AR || \
		     c == AR1 || \
                     c == ARMA || \
                     c == GARCH)

#define SIMPLE_AR_MODEL(c) (c == AR || c == AR1)

#define ML_ESTIMATOR(c) (c == ARMA || \
			 c == DURATION || \
                         c == GARCH || \
                         c == HECKIT || \
                         c == LOGIT || \
                         c == MLE || \
			 c == NEGBIN ||	\
                         c == POISSON || \
                         c == PROBIT || \
			 c == BIPROBIT || \
                         c == TOBIT)

#define LIMDEP(c) (c == LOGIT || \
                   c == PROBIT || \
                   c == TOBIT || \
                   c == INTREG)

#define COUNT_MODEL(c) (c == POISSON || c == NEGBIN)

#define LSQ_MODEL(c) (c == AR1 || \
                      c == HSK || \
                      c == OLS || \
                      c == WLS)

#define ASYMPTOTIC_MODEL(c) (c == ARMA || \
			     c == DPANEL ||   \
			     c == DURATION || \
                             c == GARCH || \
                             c == GMM || \
                             c == HECKIT || \
                             c == INTREG || \
                             c == IVREG || \
                             c == LOGIT || \
                             c == MLE || \
			     c == NEGBIN || \
                             c == POISSON || \
                             c == PROBIT || \
                             c == TOBIT || \
                             c == BIPROBIT)

#define EQN_SYSTEM_COMMAND(c) (c == VAR || c == VECM || c == SYSTEM)

/* model where the specification is not based on a list
   of variables */
#define NONLIST_MODEL(c) (c == NLS || c == MLE || c == GMM || c == MIDASREG)

#define is_model_ref_cmd(c) (c == ADD || \
	                     c == ARCH || \
	                     c == CHOW || \
	                     c == CUSUM || \
	                     c == FCAST || \
                             c == LEVERAGE || \
	                     c == MODTEST || \
                             c == OMIT || \
	                     c == RESTRICT || \
                             c == VIF)

#define RQ_SPECIAL_MODEL(m) ((m->ci == LAD || m->ci == QUANTREG) &&	\
                             NULL != gretl_model_get_data(m, "rq_tauvec"))

#define POOLED_MODEL(m) ((m->ci == OLS || m->ci == PANEL) && \
                         gretl_model_get_int(m, "pooled"))

typedef enum {
    GRETL_TEST_ADD,
    GRETL_TEST_ARCH,
    GRETL_TEST_AUTOCORR,
    GRETL_TEST_CHOW,
    GRETL_TEST_CUSUM,
    GRETL_TEST_QLR,
    GRETL_TEST_GROUPWISE,
    GRETL_TEST_LOGS,
    GRETL_TEST_NORMAL,
    GRETL_TEST_OMIT,
    GRETL_TEST_RESET,
    GRETL_TEST_SQUARES,
    GRETL_TEST_WHITES,
    GRETL_TEST_SARGAN,
    GRETL_TEST_IV_HAUSMAN,
    GRETL_TEST_PANEL_HAUSMAN,
    GRETL_TEST_PANEL_F,
    GRETL_TEST_PANEL_BP,
    GRETL_TEST_PANEL_TIMEDUM,
    GRETL_TEST_PANEL_AR,
    GRETL_TEST_HET_1,
    GRETL_TEST_BP,
    GRETL_TEST_CHOWDUM,
    GRETL_TEST_COMFAC,
    GRETL_TEST_INDEP,
    GRETL_TEST_RE,
    GRETL_TEST_WITHIN_F,
    GRETL_TEST_PANEL_WELCH,
    GRETL_TEST_RE_WALD,
    GRETL_TEST_XDEPEND,
    GRETL_TEST_MAX
} ModelTestType;

MODEL *gretl_model_new (void);

void gretl_model_init (MODEL *pmod, const DATASET *dset);

int gretl_model_allocate_storage (MODEL *pmod);

MODEL **gretl_model_array_new (int n);

MODEL *allocate_working_model (void);

void gretl_model_array_destroy (MODEL **models, int n);

void destroy_working_model (MODEL *models);

void gretl_model_smpl_init (MODEL *pmod, const DATASET *dset);

void impose_model_smpl (const MODEL *pmod, DATASET *dset);

void gretl_model_set_auxiliary (MODEL *pmod, ModelAuxCode aux);

void clear_model (MODEL *pmod);

void gretl_model_free (MODEL *pmod);

void clear_model_xpx (MODEL *pmod);

void gretl_model_free_on_exit (MODEL *pmod);

void display_model_data_items (const MODEL *pmod);

int bundlize_model_data_items (const MODEL *pmod, gretl_bundle *b);

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr,
			  GretlType type, size_t size);

int gretl_model_set_matrix_as_data (MODEL *pmod, const char *key,
				    gretl_matrix *m);

int gretl_model_set_list_as_data (MODEL *pmod, const char *key, int *list);

int gretl_model_set_string_as_data (MODEL *pmod, const char *key, char *str);

int gretl_model_set_array_as_data (MODEL *pmod, const char *key,
				   gretl_array *A);

int gretl_model_destroy_data_item (MODEL *pmod, const char *key);

int gretl_model_detach_data_item (MODEL *pmod, const char *key);

int gretl_model_set_int (MODEL *pmod, const char *key, int val);

int gretl_model_set_double (MODEL *pmod, const char *key, double val);

int gretl_model_set_vcv_info (MODEL *pmod, int vmaj, int vmin);

int gretl_model_set_hac_vcv_info (MODEL *pmod, int kern,
				  int order, int flags,
				  double bw);

int gretl_model_set_hac_order (MODEL *pmod, int order);

int gretl_model_set_cluster_vcv_info (MODEL *pmod,
				      const char *cv1,
				      const char *cv2);

int gretl_model_get_vcv_type (const MODEL *pmod);

int gretl_model_get_hc_version (const MODEL *pmod);

const char *gretl_model_get_cluster_vname (const MODEL *pmod);

const char *gretl_model_get_cluster_vname2 (const MODEL *pmod);

void *gretl_model_get_data (const MODEL *pmod, const char *key);

void *gretl_model_get_data_full (const MODEL *pmod, const char *key,
				 GretlType *type, int *copied,
				 size_t *sz);

void *gretl_model_steal_data (MODEL *pmod, const char *key);

int gretl_model_get_int (const MODEL *pmod, const char *key);

double gretl_model_get_double (const MODEL *pmod, const char *key);

double gretl_model_get_double_default (const MODEL *pmod,
				       const char *key,
				       double deflt);

int *gretl_model_get_list (const MODEL *pmod, const char *key);

char *gretl_model_get_param_name (const MODEL *pmod,
				  const DATASET *dset,
				  int i, char *targ);

gretl_array *gretl_model_get_param_names (const MODEL *pmod,
					  const DATASET *dset,
					  int *err);

int gretl_model_get_param_number (const MODEL *pmod,
				  const DATASET *dset,
				  const char *s);

void free_coeff_intervals (CoeffIntervals *cf);

CoeffIntervals *
gretl_model_get_coeff_intervals (const MODEL *pmod,
				 const DATASET *dset,
				 gretlopt opt);

int reset_coeff_intervals (CoeffIntervals *cf, double alpha);

gretl_matrix *conf_intervals_matrix (CoeffIntervals *cf);

int gretl_model_get_depvar (const MODEL *pmod);

const char *gretl_model_get_depvar_name (const MODEL *pmod,
					 const DATASET *dset);

int *gretl_model_get_x_list (const MODEL *pmod);

int *gretl_model_get_y_list (const MODEL *pmod);

int *gretl_model_get_secondary_list (const MODEL *pmod);

int arma_model_nonseasonal_AR_order (const MODEL *pmod);

int arma_model_nonseasonal_MA_order (const MODEL *pmod);

int arma_model_max_AR_lag (const MODEL *pmod);

int arma_model_max_MA_lag (const MODEL *pmod);

int arma_model_AR_MA_coeffs (const MODEL *pmod,
			     gretl_vector **phi_star,
			     gretl_vector **theta_star,
			     gretlopt opt);

int regarma_model_AR_coeffs (const MODEL *pmod,
			     double **phi0,
			     int *pp);

const double *arma_model_get_x_coeffs (const MODEL *pmod);

int arma_model_get_n_arma_coeffs (const MODEL *pmod);

int regarima_model_get_AR_coeffs (const MODEL *pmod,
				  double **phi0,
				  int *pp);

int *arima_delta_coeffs (int d, int D, int s);

gretl_matrix *arma_spectrum_plot_data (const MODEL *pmod,
				       const DATASET *dset,
				       int *err);

gretl_matrix *gretl_model_ahat_vec (const MODEL *pmod, int *err);

int gretl_model_set_coeff_separator (MODEL *pmod, const char *s, int pos);

int gretl_model_get_coeff_separator (const MODEL *pmod, char **ps, int *ppos);

int gretl_model_new_vcv (MODEL *pmod, int *nelem);

int gretl_model_write_vcv (MODEL *pmod, const gretl_matrix *V);

int gretl_model_add_QML_vcv (MODEL *pmod, int ci,
			     const gretl_matrix *H,
			     const gretl_matrix *G,
			     const DATASET *dset,
			     gretlopt opt,
			     gretl_matrix **pV);

int gretl_model_add_hessian_vcv (MODEL *pmod,
				 const gretl_matrix *H);

int gretl_model_add_OPG_vcv (MODEL *pmod,
			     const gretl_matrix *G,
			     gretl_matrix **pV);

VMatrix *gretl_model_get_vcv (MODEL *pmod, const DATASET *dset);

double gretl_model_get_vcv_element (const MODEL *pmod,
				    int i, int j,
				    int np);

int gretl_model_write_coeffs (MODEL *pmod, double *b, int k);

int gretl_model_add_arinfo (MODEL *pmod, int nterms);

MODEL *gretl_model_copy (MODEL *pmod);

void swap_models (MODEL *targ, MODEL *src);

int command_ok_for_model (int test_ci, gretlopt opt,
			  const MODEL *pmod);

int model_test_ok (int ci, gretlopt opt, const MODEL *pmod,
		   const DATASET *dset);

int gretl_is_simple_OLS (const MODEL *pmod);

int gretl_is_arima_model (const MODEL *pmod);

int gretl_is_between_model (const MODEL *pmod);

int gretl_is_regular_panel_model (const MODEL *pmod);

int get_first_model_stat (const char **word, const char **desc);

int get_next_model_stat (const char **word, const char **desc);

int get_model_count (void);

void set_model_count (int c);

int model_count_plus (void);

void model_count_minus (MODEL *pmod);

void set_model_id (MODEL *pmod, gretlopt opt);

ModelTest *model_test_new (ModelTestType ttype);

ModelTest *gretl_model_get_test (MODEL *pmod, ModelTestType ttype);

void model_test_free (ModelTest *test);

int maybe_add_test_to_model (MODEL *pmod, ModelTest *test);

void model_test_set_teststat (ModelTest *test, unsigned char ts);
void model_test_set_order (ModelTest *test, int order);
void model_test_set_dfn (ModelTest *test, int df);
void model_test_set_dfd (ModelTest *test, double df);
void model_test_set_value (ModelTest *test, double val);
void model_test_set_pvalue (ModelTest *test, double pval);
void model_test_set_param (ModelTest *test, const char *s);
void model_test_set_opt (ModelTest *test, gretlopt opt);
void model_test_set_allocated_param (ModelTest *test, char *s);
void model_test_set_crit_and_alpha (ModelTest *test,
				    double crit,
				    double alpha);

void gretl_model_test_print (const MODEL *pmod, int i, PRN *prn);
void gretl_model_print_last_test (const MODEL *pmod, PRN *prn);
void gretl_model_test_print_direct (const ModelTest *test, int heading,
				    PRN *prn);

const char *get_h0_string_for_test (ModelTestType ttype);

void gretl_model_destroy_tests (MODEL *pmod);

void model_list_to_string (int *list, char *buf);

int highest_numbered_var_in_model (const MODEL *pmod,
				   const DATASET *dset);

int mle_criteria (MODEL *pmod, int addk);

int model_use_zscore (const MODEL *pmod);

double coeff_pval (int ci, double x, int df);

double model_coeff_pval (const MODEL *pmod, double x);

int exact_fit_check (const MODEL *pmod, PRN *prn);

void maybe_suppress_time_dummies (MODEL *pmod, int ndum);

int gretl_model_allocate_param_names (MODEL *pmod, int k);

int gretl_model_set_param_name (MODEL *pmod, int i, const char *name);

int gretl_model_add_arma_varnames (MODEL *pmod, const DATASET *dset,
				   int yno, int p, int q,
				   const char *pmask, const char *qmask,
				   int P, int Q,
				   int r);

int gretl_model_add_panel_varnames (MODEL *pmod, const DATASET *dset,
				    const int *ulist);

void gretl_model_add_allocated_varnames (MODEL *pmod, char **vnames);

int gretl_model_add_y_median (MODEL *pmod, const double *y);

int gretl_model_add_normality_test (MODEL *pmod, double X2);

int gretl_model_get_normality_test (const MODEL *pmod, PRN *prn);

char *gretl_model_get_fitted_formula (const MODEL *pmod, int xvar,
				      const DATASET *dset);

void gretl_model_set_name (MODEL *pmod, const char *name);

const char *gretl_model_get_name (const MODEL *pmod);

double gretl_model_get_scalar (MODEL *pmod,
			       ModelDataIndex idx,
			       DATASET *dset,
			       int *err);

int gretl_model_get_series (double *x, MODEL *pmod,
			    const DATASET *dset,
			    ModelDataIndex idx);

gretl_matrix *gretl_model_get_matrix (MODEL *pmod,
				      ModelDataIndex idx,
				      int *err);

double
gretl_model_get_data_element (MODEL *pmod, int idx, const char *s,
			      const DATASET *dset, int *err);

int gretl_model_serialize (const MODEL *pmod, SavedObjectFlags flags,
			   PRN *prn);

#ifdef FULL_XML_HEADERS

int attach_model_tests_from_xml (MODEL *pmod, xmlNodePtr node);

MODEL *gretl_model_from_XML (xmlNodePtr node, xmlDocPtr doc,
			     const DATASET *dset,
			     int *err);
#endif

#ifdef  __cplusplus
}
#endif

#endif /* GRETL_MODEL_H */
