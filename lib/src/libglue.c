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

/* Miscellaneous functions to bridge between gretl commands and
   the corresponding libgretl functions. This glue allows a cleaner
   interface for the latter, hiding some parsing of strings
   coming from the command line.
*/

#include "libgretl.h"
#include "var.h"
#include "system.h"
#include "gretl_panel.h"
#include "libglue.h"

/*
 * model_test_driver:
 * @param: auxiliary parameter for some uses.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: controls which test(s) will be performed; %OPT_Q
 * gives less verbose results.
 * @prn: gretl printing struct.
 * 
 * Performs some subset of gretl's "modtest" tests on the
 * model last estimated, and prints the results to @prn.
 * 
 * Returns: 0 on successful completion, error code on error.
 */

int model_test_driver (const char *param, 
		       double ***pZ, DATAINFO *pdinfo, 
		       gretlopt opt, PRN *prn)
{
    GretlObjType type;
    gretlopt testopt;
    void *ptr;
    int k = 0;
    int err = 0;

    if (opt == OPT_NONE || opt == OPT_Q) {
	pprintf(prn, "modtest: no options selected\n");
	return 0;
    }

    err = incompatible_options(opt, OPT_A | OPT_H | OPT_L | OPT_S |
			       OPT_N | OPT_P | OPT_W | OPT_X);
    if (err) {
	return err;
    }

    ptr = get_last_model(&type);  
    if (ptr == NULL) {
	return E_DATA;
    }

    if (type == GRETL_OBJ_EQN && exact_fit_check(ptr, prn)) {
	return 0;
    }

    if (opt & (OPT_A | OPT_H)) {
	/* autocorrelation and arch: lag order */
	k = atoi(param);
	if (k == 0) {
	    k = pdinfo->pd;
	}
    }

    testopt = (opt & OPT_Q)? OPT_Q : OPT_NONE;

    /* non-linearity (squares) */
    if (!err && (opt & OPT_S)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, pZ, pdinfo, 
				    AUX_SQ, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* non-linearity (logs) */
    if (!err && (opt & OPT_L)) {
	if (type == GRETL_OBJ_EQN) {
	    err = nonlinearity_test(ptr, pZ, pdinfo, 
				    AUX_LOG, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* heteroskedasticity (White or Breusch-Pagan) */
    if (!err && (opt & (OPT_W | OPT_X | OPT_B))) {
	if (type == GRETL_OBJ_EQN) {
	    transcribe_option_flags(&testopt, opt, OPT_B | OPT_X);
	    if ((opt & OPT_B) && (opt & OPT_R)) {
		testopt |= OPT_R;
	    }
	    err = whites_test(ptr, pZ, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* autocorrelation */
    if (!err && (opt & OPT_A)) {
	if (type == GRETL_OBJ_EQN) {
	    err = autocorr_test(ptr, k, pZ, pdinfo, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_autocorrelation_test(ptr, k, pZ, pdinfo, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_autocorrelation_test(ptr, k, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* ARCH */
    if (!err && (opt & OPT_H)) {
	if (type == GRETL_OBJ_EQN) {
	    err = arch_test(ptr, k, pdinfo, testopt, prn);
	} else if (type == GRETL_OBJ_VAR) {
	    err = gretl_VAR_arch_test(ptr, k, pdinfo, prn);
	} else if (type == GRETL_OBJ_SYS) {
	    err = system_arch_test(ptr, k, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    /* normality of residual */
    if (!err && (opt & OPT_N)) {
	err = last_model_test_uhat(pZ, pdinfo, testopt, prn);
    }

    /* groupwise heteroskedasticity */
    if (!err && (opt & OPT_P)) {
	if (type == GRETL_OBJ_EQN) {
	    err = groupwise_hetero_test(ptr, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }

    /* common factor restriction */
    if (!err && (opt & OPT_C)) {
	if (type == GRETL_OBJ_EQN) {
	    err = comfac_test(ptr, pZ, pdinfo, testopt, prn);
	} else {
	    err = E_NOTIMP;
	}
    }    

    return err;
}

static int get_chow_dummy (const char *s, const double **Z,
			   const DATAINFO *pdinfo, int *err)
{
    int v = current_series_index(pdinfo, s);

    if (v < 0) {
	*err = E_UNKVAR;
    } else if (!gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[v])) {
	*err = E_DATA;
    }

    return v;
}

/*
 * chow_test_driver:
 * @line: command line for parsing.
 * @pmod: pointer to model to be tested.
 * @pZ: pointer to data array.
 * @pdinfo: dataset information.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_D included, do the Chow test based on a given dummy
 * variable.
 * @prn: gretl printing struct.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int chow_test_driver (const char *line, MODEL *pmod, double ***pZ,
		      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    char chowstr[VNAMELEN];
    int chowparm = 0;
    int err = 0;

    /* chowstr should hold either an observation at which to
       split the sample, or the name of a dummy variable
       to be used to divide the sample (if given OPT_D)
    */

    if (sscanf(line, "%*s %15s", chowstr) != 1) {
	err = E_PARSE;
    } else if (opt & OPT_D) {
	chowparm = get_chow_dummy(chowstr, (const double **) *pZ, pdinfo, &err);
    } else {
	chowparm = dateton(chowstr, pdinfo);
    }

    if (!err) {
	if (opt & OPT_D) {
	    err = chow_test_from_dummy(chowparm, pmod, pZ, pdinfo, opt, prn);
	} else {
	    err = chow_test(chowparm, pmod, pZ, pdinfo, opt, prn);
	}
    }

    return err;
}

/* parse the tau vector out of @parm before calling the
   "real" quantreg function
*/

MODEL quantreg_driver (const char *parm, const int *list, 
		       double ***pZ, DATAINFO *pdinfo,
		       gretlopt opt, PRN *prn)
{
    gretl_vector *tau;
    MODEL mod;
    int err = 0;

    tau = generate_matrix(parm, pZ, pdinfo, &err);

    if (!err && gretl_vector_get_length(tau) == 0) {
	err = E_DATA;
    }

    if (err) {
	gretl_model_init(&mod);
	mod.errcode = err;
    } else {
	mod = quantreg(tau, list, pZ, pdinfo, opt, prn);
    }

    gretl_matrix_free(tau);

    return mod;
}

/* wrapper for the various sorts of logit and probit models
   that gretl supports
*/

MODEL logit_probit (int *list, double ***pZ, DATAINFO *pdinfo, 
		    int ci, gretlopt opt, PRN *prn)
{
    int yv = list[1];

    if (ci == LOGIT && (opt & OPT_M)) {
	return multinomial_logit(list, pZ, pdinfo, opt, prn);
    } else if (gretl_isdummy(pdinfo->t1, pdinfo->t2, (*pZ)[yv])) {
	if (ci == LOGIT) {
	    return binary_logit(list, pZ, pdinfo, opt, prn);
	} else {
	    return binary_probit(list, pZ, pdinfo, opt, prn);
	}
    } else {
	if (ci == LOGIT) {
	    return ordered_logit(list, pZ, pdinfo, opt, prn);
	} else {
	    return ordered_probit(list, pZ, pdinfo, opt, prn);
	}
    } 
}

/* parse out optional "ymax=..." parameter before called the real
   logistic model function 
*/

MODEL logistic_driver (const int *list, double ***pZ, DATAINFO *pdinfo,
		       const char *param) 
{
    double lmax;

    if (param == NULL || sscanf(param, "ymax=%lf", &lmax) != 1) {
	lmax = NADBL;
    }

    return logistic_model(list, lmax, pZ, pdinfo);
}

    

