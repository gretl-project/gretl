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

typedef enum {
    ARMA_SEAS  = 1 << 0, /* includes seasonal component */
    ARMA_DSPEC = 1 << 1, /* input list includes differences */
    ARMA_X12A  = 1 << 2, /* using X-12-ARIMA to generate estimates */
    ARMA_EXACT = 1 << 3, /* using exact ML */
    ARMA_VECH  = 1 << 4, /* using vech representation when computing
			    variance matrix of state for Kalman filter
			 */
    ARMA_LS    = 1 << 5  /* using conditional ML, and O/NLS == CML */
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
    int ifc;
};

/* single-equation model commands */

#define MODEL_COMMAND(c) (c == AR || \
                          c == AR1 || \
                          c == ARBOND || \
                          c == ARCH || \
                          c == ARMA || \
                          c == GARCH || \
                          c == GMM || \
		          c == HECKIT || \
                          c == HSK || \
                          c == INTREG || \
                          c == IVREG || \
                          c == LAD || \
                          c == LOGISTIC || \
                          c == LOGIT || \
                          c == MLE || \
                          c == MPOLS || \
                          c == NLS || \
                          c == OLS || \
                          c == PANEL || \
                          c == POISSON || \
                          c == PROBIT || \
                          c == QUANTREG || \
                          c == TOBIT || \
                          c == WLS)

#define AR_MODEL(c) (c == AR || \
		     c == AR1 || \
                     c == ARMA || \
                     c == GARCH)

#define SIMPLE_AR_MODEL(c) (c == AR || c == AR1)

#define ML_ESTIMATOR(c) (c == ARMA || \
                         c == GARCH || \
                         c == HECKIT || \
                         c == LOGIT || \
                         c == MLE || \
                         c == POISSON || \
                         c == PROBIT || \
                         c == TOBIT)

#define LIMDEP(c) (c == LOGIT || \
                   c == PROBIT || \
                   c == TOBIT)

#define LSQ_MODEL(c) (c == AR1 || \
                      c == HSK || \
                      c == OLS || \
                      c == WLS)

#define ASYMPTOTIC_MODEL(c) (c == ARBOND || \
                             c == ARMA || \
                             c == GARCH || \
                             c == GMM || \
                             c == HECKIT || \
                             c == INTREG || \
                             c == IVREG || \
                             c == LOGIT || \
                             c == MLE || \
                             c == POISSON || \
                             c == PROBIT || \
                             c == TOBIT)

/* model where the specification is not based on a list
   of variables */
#define NONLIST_MODEL(c) (c == NLS || c == MLE || c == GMM)

#define is_model_ref_cmd(c) (c == ADD || \
	                     c == ARCH || \
	                     c == CHOW || \
	                     c == CUSUM || \
	                     c == FCAST || \
                             c == LEVERAGE || \
	                     c == LMTEST || \
                             c == OMIT || \
	                     c == RESTRICT || \
	                     c == TESTUHAT || \
                             c == VIF)

#define RQ_SPECIAL_MODEL(m) (m->ci == LAD && \
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
    GRETL_TEST_HET_1,
    GRETL_TEST_BP,
    GRETL_TEST_MAX
} ModelTestType;

#ifndef GRETLCLI
int attach_model_tests_from_xml (MODEL *pmod, xmlNodePtr node);
#endif

MODEL *gretl_model_new (void);

void gretl_model_init (MODEL *pmod);

int gretl_model_allocate_storage (MODEL *pmod);

MODEL **gretl_model_array_new (int n);

MODEL **allocate_working_models (int n);

void gretl_model_array_destroy (MODEL **models, int n);

void destroy_working_models (MODEL **models, int n);

void gretl_model_smpl_init (MODEL *pmod, const DATAINFO *pdinfo);

void impose_model_smpl (const MODEL *pmod, DATAINFO *pdinfo);

void gretl_model_set_auxiliary (MODEL *pmod, ModelAuxCode aux);

void clear_model (MODEL *pmod);

void gretl_model_free (MODEL *pmod);

void gretl_model_free_on_exit (MODEL *pmod);

void display_model_data_items (const MODEL *pmod);

int gretl_model_set_data_with_destructor (MODEL *pmod, const char *key, void *ptr, 
					  GretlType type, size_t size, 
					  void (*destructor) (void *));

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr, 
			  GretlType type, size_t size);

int gretl_model_set_matrix_as_data (MODEL *pmod, const char *key, 
				    gretl_matrix *m);

int gretl_model_set_list_as_data (MODEL *pmod, const char *key, int *list);

int gretl_model_set_string_as_data (MODEL *pmod, const char *key, char *str);

int gretl_model_destroy_data_item (MODEL *pmod, const char *key);

int gretl_model_detach_data_item (MODEL *pmod, const char *key);

int gretl_model_set_int (MODEL *pmod, const char *key, int val);

int gretl_model_set_double (MODEL *pmod, const char *key, double val);

void *gretl_model_get_data (const MODEL *pmod, const char *key);

void *gretl_model_get_data_and_size (const MODEL *pmod, const char *key,
				     size_t *sz);

int gretl_model_get_int (const MODEL *pmod, const char *key);

double gretl_model_get_double (const MODEL *pmod, const char *key);

int *gretl_model_get_list (const MODEL *pmod, const char *key);

char *gretl_model_get_param_name (const MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  int i, char *targ);

int gretl_model_get_param_number (const MODEL *pmod, 
				  const DATAINFO *pdinfo,
				  const char *s);

void free_coeff_intervals (CoeffIntervals *cf);

CoeffIntervals *
gretl_model_get_coeff_intervals (const MODEL *pmod, 
				 const DATAINFO *pdinfo);

int reset_coeff_intervals (CoeffIntervals *cf, double alpha);

int gretl_model_get_depvar (const MODEL *pmod);

const char *gretl_model_get_depvar_name (const MODEL *pmod,
					 const DATAINFO *pdinfo);

int *gretl_model_get_x_list (const MODEL *pmod);

int *gretl_model_get_secondary_list (const MODEL *pmod);

int arma_model_nonseasonal_AR_order (const MODEL *pmod);

int arma_model_nonseasonal_MA_order (const MODEL *pmod);

int arma_model_max_AR_lag (const MODEL *pmod);

int arma_model_max_MA_lag (const MODEL *pmod);

int arma_model_integrated_AR_MA_coeffs (const MODEL *pmod,
					double **phi_star,
					double **theta_star);

int regarma_model_AR_coeffs (const MODEL *pmod,
			     double **phi0,
			     int *pp);

const double *arma_model_get_x_coeffs (const MODEL *pmod);

int regarima_model_get_AR_coeffs (const MODEL *pmod,
				  double **phi0,
				  int *pp);

int gretl_model_set_coeff_separator (MODEL *pmod, const char *s, int pos);

int gretl_model_get_coeff_separator (const MODEL *pmod, const char **ps, int *ppos);

int gretl_model_new_vcv (MODEL *pmod, int *nelem);

int gretl_model_write_vcv (MODEL *pmod, const gretl_matrix *V);

VMatrix *gretl_model_get_vcv (MODEL *pmod, const DATAINFO *pdinfo);

int gretl_model_add_arinfo (MODEL *pmod, int nterms);

MODEL *gretl_model_copy (const MODEL *pmod);

void swap_models (MODEL *targ, MODEL *src);

int command_ok_for_model (int test_ci, gretlopt opt, int mci);

int model_test_ok (int ci, gretlopt opt, const MODEL *pmod, 
		   const DATAINFO *pdinfo);

int gretl_is_arima_model (const MODEL *pmod);

int get_first_model_stat (const char **word, const char **desc);

int get_next_model_stat (const char **word, const char **desc);

int get_model_count (void);

void reset_model_count (void);

int model_count_plus (void);

void model_count_minus (void);

void set_model_id (MODEL *pmod);

ModelTest *model_test_new (ModelTestType ttype);
void model_test_free (ModelTest *test);

int maybe_add_test_to_model (MODEL *pmod, ModelTest *test);

void model_test_set_teststat (ModelTest *test, unsigned char ts);
void model_test_set_order (ModelTest *test, int order);
void model_test_set_dfn (ModelTest *test, int df);
void model_test_set_dfd (ModelTest *test, int df);
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

void gretl_model_destroy_tests (MODEL *pmod);

void model_list_to_string (int *list, char *buf);

int highest_numbered_var_in_model (const MODEL *pmod, 
				   const DATAINFO *pdinfo);

int mle_criteria (MODEL *pmod, int addk);

double coeff_pval (int ci, double x, int df);

int exact_fit_check (const MODEL *pmod, PRN *prn);

int gretl_model_allocate_params (MODEL *pmod, int k);

int gretl_model_add_arma_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				   int yno, int p, int q, 
				   const char *pmask, const char *qmask,
				   int P, int Q, 
				   int r);

int gretl_model_add_panel_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				    const int *ulist);

void gretl_model_add_allocated_varnames (MODEL *pmod, char **vnames);

int gretl_model_add_y_median (MODEL *pmod, const double *y);

char *gretl_model_get_fitted_formula (const MODEL *pmod, int xvar,
				      const double **Z,
				      const DATAINFO *pdinfo);

void gretl_model_set_name (MODEL *pmod, const char *name);

const char *gretl_model_get_name (const MODEL *pmod);

double gretl_model_get_scalar (const MODEL *pmod, ModelDataIndex idx, 
			       double ***pZ, DATAINFO *pdinfo,
			       int *err);

double *
gretl_model_get_series (const MODEL *pmod, const DATAINFO *pdinfo, 
			ModelDataIndex idx, int *err);

gretl_matrix *gretl_model_get_matrix (MODEL *pmod, ModelDataIndex idx, 
				      int *err);

double 
gretl_model_get_data_element (MODEL *pmod, int idx, const char *s,
			      const DATAINFO *pdinfo, int *err);

int gretl_model_serialize (const MODEL *pmod, SavedObjectFlags flags,
			   FILE *fp);

#ifndef GRETLCLI

MODEL *gretl_model_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err);

#endif

#endif /* GRETL_MODEL_H */
