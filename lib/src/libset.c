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
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

/* libset.c for gretl */

#include "libgretl.h"
#include "gretl_private.h"
#include "libset.h"

static int use_qr = -1;
static int halt_on_error = -1;
static double hp_lambda;

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE
};

static struct {
    int auto_lag;
    int user_lag;
    int hc_version;
} robust_opts;

int get_hc_version (void)
{
    return robust_opts.hc_version;
}

double get_hp_lambda (void)
{
    return hp_lambda;
}

static int get_or_set_force_hc (int f)
{
    static int force;

    if (f >= 0) force = f;
    return force;
}

static int get_or_set_garch_vcv (int v)
{
    static int variant;

    if (v >= 0) variant = v;
    return variant;
}

static int get_or_set_garch_robust_vcv (int v)
{
    static int variant;

    if (v >= 0) variant = v;
    return variant;
}

static int set_garch_vcv_variant (const char *s)
{
    int vopt = VCV_UNSET;

    if (!strcmp(s, "hessian"))  vopt = VCV_HESSIAN;
    else if (!strcmp(s, "im"))  vopt = VCV_IM;
    else if (!strcmp(s, "op"))  vopt = VCV_OP;
    else if (!strcmp(s, "qml")) vopt = VCV_QML;
    else if (!strcmp(s, "bw"))  vopt = VCV_BW;

    get_or_set_garch_vcv(vopt);
    
    return (vopt == VCV_UNSET);
}

int get_garch_vcv_version (void)
{
    return get_or_set_garch_vcv(-1);
}

int get_garch_robust_vcv_version (void)
{
    return get_or_set_garch_robust_vcv(-1);
}

static void set_force_hc (int f)
{
    get_or_set_force_hc(f);
}

int get_force_hc (void)
{
    return get_or_set_force_hc(-1);
}

int get_hac_lag (int m)
{
    /* Variants of Newey-West */

    if (robust_opts.user_lag != 0 && robust_opts.user_lag < m - 2) {
	/* FIXME upper limit? */
	return robust_opts.user_lag;
    }

    if (robust_opts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return 0.75 * pow(m, 1.0 / 3.0);
    } else if (robust_opts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	return 4.0 * pow(m / 100.0, 2.0 / 9.0);
    }

    /* fallback -- should not be reached */
    return 0.75 * pow(m, 1.0 / 3.0);
}

static int parse_hc_variant (const char *s)
{
    int err = 1;

    if (!strcmp(s, "0") || !strcmp(s, "1") ||
	!strcmp(s, "2") || !strcmp(s, "3")) {
	robust_opts.hc_version = atoi(s);
	err = 0;
    } else if (!strcmp(s, "3a")) {
	robust_opts.hc_version = 4;
	err = 0;
    }

    if (err) {
	int hcv;

	if (!strcmp(s, "hc3a")) {
	    robust_opts.hc_version = 4;
	    err = 0;
	} else if (sscanf(s, "hc%d", &hcv)) {
	    if (hcv >= 0 && hcv <= 4) {
		robust_opts.hc_version = hcv;
		err = 0;
	    }
	}
    }

    return err;
}

void set_xsect_hccme (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    parse_hc_variant(scpy);
    free(scpy);
}

void set_tseries_hccme (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    if (parse_hc_variant(scpy) == 0) {
	set_force_hc(1);
    } else {
	set_force_hc(0);
    }
    free(scpy);
}

void set_garch_robust_vcv (const char *s)
{
    char *scpy = gretl_strdup(s);

    if (scpy == NULL) return;

    lower(scpy);
    if (!strcmp(s, "qml")) {
	get_or_set_garch_robust_vcv(VCV_QML);
    } else if (!strcmp(s, "bw")) {
	get_or_set_garch_robust_vcv(VCV_BW);
    }
    free(scpy);
}

int parse_set_line (const char *line, int *echo_off, PRN *prn)
{
    char setobj[16], setarg[16];
    int nw, err = E_PARSE;

    *setobj = *setarg = '\0';

    nw = sscanf(line, "%*s %15s %15s", setobj, setarg);
    
    if (nw == 1) {
	if (!strcmp(setobj, "echo")) {
	    *echo_off = 0;
	    err = 0;
	}
    } else if (nw == 2) {
	lower(setarg);

	/* set echo on/off */
	if (!strcmp(setobj, "echo")) {
	    if (echo_off != NULL) {
		if (!strcmp(setarg, "off")) {
		    *echo_off = 1;
		    err = 0;
		} else if (!strcmp(setarg, "on")) {
		    *echo_off = 0;
		    err = 0;
		}
	    } 
	} else if (!strcmp(setobj, "hac_lag")) {
	    /* set max lag for HAC estimation */
	    if (!strcmp(setarg, "nw1")) {
		robust_opts.auto_lag = AUTO_LAG_STOCK_WATSON;
		robust_opts.user_lag = 0;
		err = 0;
	    } else if (!strcmp(setarg, "nw2")) {
		robust_opts.auto_lag = AUTO_LAG_WOOLDRIDGE;
		robust_opts.user_lag = 0;
		err = 0;
	    } else if (isdigit(*setarg)) {
		robust_opts.user_lag = atoi(setarg);
		err = 0;
	    }
	} else if (!strcmp(setobj, "hc_version")) {
	    /* set HCCM variant */
	    err = parse_hc_variant(setarg);
	} else if (!strcmp(setobj, "force_hc")) {
	    /* use HCCM, not HAC, even for time series */
	    if (!strcmp(setarg, "on")) { 
		set_force_hc(1);
		err = 0;
	    } else if (!strcmp(setarg, "off")) { 
		set_force_hc(0);
		err = 0;
	    }
	} else if (!strcmp(setobj, "garch_vcv")) {
	    /* set GARCH VCV variant */
	    err = set_garch_vcv_variant(setarg);
	} else if (!strcmp(setobj, "qr")) {
	    /* switch QR vs Cholesky decomposition */
	    if (!strcmp(setarg, "on")) {
		use_qr = 1;
		err = 0;
	    }
	    else if (!strcmp(setarg, "off")) {
		use_qr = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "halt_on_error")) {
	    if (!strcmp(setarg, "on")) {
		halt_on_error = 1;
		err = 0;
	    }
	    else if (!strcmp(setarg, "off")) {
		halt_on_error = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "seed")) {
	    /* seed for PRNG */
	    if (isdigit(*setarg)) {
		int k = atoi(setarg);

		gretl_rand_set_seed(k);
		pprintf(prn, 
			_("Pseudo-random number generator seeded with %d\n"), k);
		err = 0;

	    }
	} else if (!strcmp(setobj, "hp_lambda")) {
	    /* Hodrick-Prescott filter parameter */
	    if (!strcmp(setarg, "auto")) {
		hp_lambda = 0.0;
		err = 0;
	    } else if (check_atof(setarg) == 0) {
		hp_lambda = atof(setarg);
		err = 0;
	    }
	}	
    }
		    
    return err;
}

/* use Cholesky or QR for regression? */

void set_use_qr (int set)
{
    use_qr = set;
}

int get_use_qr (void)
{
    /* if use_qr has not been set explicitly, try env */
    if (use_qr == -1) {
	char *s = getenv("GRETL_USE_QR");

	if (s != NULL && *s != '\0' && *s != '0') {
	    use_qr = 1;
	} else {
	    use_qr = 0;
	}
    } 

    return use_qr;
}

int get_halt_on_error (void)
{
    /* if halt_on_error has not been set explicitly, try env */
    if (halt_on_error == -1) {
	char *s = getenv("GRETL_KEEP_GOING");

	if (s != NULL && *s != '\0' && *s != '0') {
	    halt_on_error = 0;
	} else {
	    halt_on_error = 1;
	}
    } 

    return halt_on_error;
}

/* pause between screens of output? (cli operation, not in
   batch mode) */

static int gretl_text_pause;

void gretl_set_text_pause (int p)
{
    gretl_text_pause = p;
}

int gretl_get_text_pause (void)
{
    return gretl_text_pause;
}

