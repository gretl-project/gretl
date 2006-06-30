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
#include "libset.h"

#define PDEBUG 0

enum {
    AUTO_LAG_STOCK_WATSON,
    AUTO_LAG_WOOLDRIDGE
};

/* for values that really want a non-negative integer */
#define UNSET_INT -1
#define is_unset(i) (i == UNSET_INT)

#define INITVALS_MAX 20

typedef struct set_vars_ set_vars;

struct robust_opts {
    int auto_lag;
    int user_lag;
    int hc_version;
    int force_hc;
};

struct garch_opts {
    int vcv_variant;
    int robust_vcv_variant;
};

struct bkbp_opts {
    int k;
    int periods[2];
};

struct sample_info {
    int t1;
    int t2;
    char submode;
    char *submask;
};
    
struct bhhh_opts {
    double tolerance;
    int maxiter;
};

struct set_vars_ {
    int use_qr;                 /* use QR decomposition? */
    unsigned int seed;          /* for PRNG */
    int halt_on_error;          /* halt cli program on script error? */
    int shell_ok;               /* shell commands permitted? */
    double hp_lambda;           /* for Hodrick-Prescott filter */
    int horizon;                /* for VAR impulse responses */ 
    double nls_toler;           /* NLS convergence criterion */
    int gretl_echo;             /* echoing commands or not */
    int gretl_msgs;             /* emitting non-error messages or not */
    char delim;                 /* delimiter for CSV data export */
    int longdigits;             /* digits for printing data in long form */
    double *initvals;           /* for parameter initialization */
    struct robust_opts ropts;   /* robust standard error options */
    struct garch_opts gopts;    /* GARCH covariance matrix */
    struct bkbp_opts bkopts;    /* Baxter-King filter */
    struct sample_info sinfo;   /* record of dataset sample state */
    struct bhhh_opts maxopts;   /* options for BHHH maximisation */
};

/* global state */
set_vars *state;

static void robust_opts_init (struct robust_opts *opts)
{
    opts->auto_lag = AUTO_LAG_STOCK_WATSON;
    opts->user_lag = 0;
    opts->hc_version = 0;
    opts->force_hc = 0; 
}

static void robust_opts_copy (struct robust_opts *opts)
{
    opts->auto_lag = state->ropts.auto_lag;
    opts->user_lag = state->ropts.user_lag;
    opts->hc_version = state->ropts.hc_version;
    opts->force_hc = state->ropts.force_hc; 
}

static void garch_opts_init (struct garch_opts *opts)
{
    opts->vcv_variant = VCV_UNSET;
    opts->robust_vcv_variant = VCV_UNSET;
}

static void garch_opts_copy (struct garch_opts *opts)
{
    opts->vcv_variant = state->gopts.vcv_variant;
    opts->robust_vcv_variant = state->gopts.robust_vcv_variant;
}

static void bkbp_opts_init (struct bkbp_opts *opts)
{
    opts->k = 8;
    opts->periods[0] = 8;
    opts->periods[1] = 32;
}

static void bkbp_opts_copy (struct bkbp_opts *opts)
{
    opts->k = state->bkopts.k;
    opts->periods[0] = state->bkopts.periods[0];
    opts->periods[1] = state->bkopts.periods[1];
}

static void bhhh_opts_init (struct bhhh_opts *opts)
{
    opts->tolerance = NADBL;
    opts->maxiter = 500;
}

static void bhhh_opts_copy (struct bhhh_opts *opts)
{
    opts->tolerance = state->maxopts.tolerance;
    opts->maxiter = state->maxopts.maxiter;
}

static double *initvals_copy (double *ivals)
{
    int n;

    if (ivals == NULL) {
	return NULL;
    }

    for (n=0; ivals[n] != NADBL; n++) {
	if (n > INITVALS_MAX) {
	    return NULL;
	}
    }

    return copyvec(ivals, n);
}

static void print_initvals (double *ivals, PRN *prn)
{
    int i, n = 0;

    if (ivals != NULL) {
	for (n=0; ivals[n] != NADBL; n++) {
	    if (n > INITVALS_MAX) {
		n = 0;
		break;
	    }
	}
    }

    if (n == 0) {
	pputs(prn, " initvals = auto\n");
    } else {
	pputs(prn, " initvals =");
	for (i=0; i<n - 1; i++) {
	    pprintf(prn, " %#.8g", ivals[i]);
	}
	pputc(prn, '\n');
    }
}

static void sample_info_init (struct sample_info *sinfo)
{
    sinfo->t1 = UNSET_INT;
    sinfo->t2 = UNSET_INT;
    sinfo->submask = NULL;
    sinfo->submode = 0;
}

#define sinfo_is_set(s) (s.t1 != UNSET_INT && s.t2 != UNSET_INT)

static void check_for_state (void) 
{
    if (state != NULL) {
	libset_init();
    }
}

static void state_vars_copy (set_vars *sv, const DATAINFO *pdinfo)
{
#if PDEBUG
    fprintf(stderr, "state_vars_copy called: pdinfo=%p\n", (void *) pdinfo);
#endif
    sv->use_qr = state->use_qr;
    sv->seed = state->seed;
    sv->halt_on_error = state->halt_on_error;
    sv->shell_ok = state->shell_ok;
    sv->hp_lambda = state->hp_lambda;
    sv->horizon = state->horizon;
    sv->nls_toler = state->nls_toler;
    sv->gretl_echo = state->gretl_echo; 
    sv->gretl_msgs = state->gretl_msgs; 
    sv->delim = state->delim; 
    sv->longdigits = state->longdigits; 
    sv->initvals = initvals_copy(state->initvals);

    robust_opts_copy(&sv->ropts);
    garch_opts_copy(&sv->gopts);
    bkbp_opts_copy(&sv->bkopts);
    bhhh_opts_copy(&sv->maxopts);

    if (pdinfo != NULL) {
	sv->sinfo.t1 = pdinfo->t1;
	sv->sinfo.t2 = pdinfo->t2;
	sv->sinfo.submode = pdinfo->submode;
	sv->sinfo.submask = copy_datainfo_submask(pdinfo);
#if PDEBUG
	fprintf(stderr, " sinfo.t1=%d, sinfo.t2=%d, "
		"sinfo.submode=%d, sinfo.submask=%p\n", sv->sinfo.t1,
		sv->sinfo.t2, sv->sinfo.submode,
		(void *) sv->sinfo.submask);
#endif
    } else {
	sample_info_init(&sv->sinfo);
    }
}

static void state_vars_init (set_vars *sv)
{
#if PDEBUG
    fprintf(stderr, "state_vars_init called\n");
#endif
    sv->use_qr = UNSET_INT; 
    sv->seed = 0;
    sv->halt_on_error = UNSET_INT;
    sv->shell_ok = 0;
    sv->hp_lambda = NADBL;
    sv->horizon = UNSET_INT;
    sv->nls_toler = NADBL;
    sv->gretl_echo = 1; 
    sv->gretl_msgs = 1; 
    sv->delim = UNSET_INT;
    sv->longdigits = 10;
    sv->initvals = NULL;

    robust_opts_init(&sv->ropts);
    garch_opts_init(&sv->gopts);
    bkbp_opts_init(&sv->bkopts);
    bhhh_opts_init(&sv->maxopts);
    sample_info_init(&sv->sinfo);
}

int get_hc_version (void)
{
    check_for_state();
    return state->ropts.hc_version;
}

double get_hp_lambda (void)
{
    check_for_state();
    return state->hp_lambda;
}

int get_bkbp_k (void)
{
    check_for_state();
    return state->bkopts.k;
}

void get_bkbp_periods (int *periods)
{
    check_for_state();
    periods[0] = state->bkopts.periods[0];
    periods[1] = state->bkopts.periods[1];
}

double get_bhhh_toler (void)
{
    check_for_state();
    return state->maxopts.tolerance;
}

int get_bhhh_maxiter (void)
{
    check_for_state();
    return state->maxopts.maxiter;
}

int set_bhhh_toler (double tol)
{
    int err = 0;

    check_for_state();

    if (tol <= 0.0) {
	err = 1;
    } else {
	state->maxopts.tolerance = tol;
    }

    return err;
}

int set_bhhh_maxiter (int n)
{
    int err = 0;

    check_for_state();

    if (n < 1) {
	err = 1;
    } else {
	state->maxopts.maxiter = n;
    }

    return err;
}

int set_long_digits (int n)
{
    int err = 0;

    check_for_state();

    if (n < 1 || n > 20) {
	err = 1;
    } else {
	state->longdigits = n;
    }

    return err;
}


int get_VAR_horizon (void)
{
    check_for_state();
    return state->horizon;
}

double get_nls_toler (void)
{
    check_for_state();

    if (na(state->nls_toler)) {
	state->nls_toler = get_default_nls_toler();
    }

    return state->nls_toler;
}

int set_nls_toler (double tol)
{
    int err = 0;

    check_for_state();

    if (tol <= 0.0) {
	err = 1;
    } else {
	state->nls_toler = tol;
    }

    return err;
}

static int get_or_set_force_hc (int f)
{
    check_for_state();

    if (f >= 0) {
	state->ropts.force_hc = f;
    }

    return state->ropts.force_hc;
}

static int get_or_set_garch_vcv (int v)
{
    check_for_state();

    if (v >= 0) {
	state->gopts.vcv_variant = v;
    }

    return state->gopts.vcv_variant;
}

static int get_or_set_garch_robust_vcv (int v)
{
    check_for_state();

    if (v >= 0) {
	state->gopts.robust_vcv_variant = v;
    }

    return state->gopts.robust_vcv_variant;
}

struct g_vcv {
    int opt;
    char *str;
};

static struct g_vcv g_vcv_types[] = {
    { VCV_UNSET,   "unset" },
    { VCV_HESSIAN, "hessian" },
    { VCV_IM,      "im" },
    { VCV_OP,      "op" },
    { VCV_QML,     "qml" },
    { VCV_BW,      "bw" }
};

static const char *garch_vcv_string (void)
{
    int i, n = sizeof g_vcv_types / sizeof g_vcv_types[0];
    const char *ret = "";

    for (i=0; i<n; i++ ) {
	if (state->gopts.vcv_variant == g_vcv_types[i].opt) {
	    ret = g_vcv_types[i].str;
	    break;
	}
    }

    return ret;
}

static int set_garch_vcv_variant (const char *s)
{
    int i, vopt = VCV_UNSET;
    int n = sizeof g_vcv_types / sizeof g_vcv_types[0];

    for (i=0; i<n; i++ ) {
	if (!strcmp(s, g_vcv_types[i].str)) {
	    vopt = g_vcv_types[i].opt;
	    break;
	}
    }	

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

void set_gretl_echo (int e)
{
    check_for_state();
    state->gretl_echo = e;
}

int gretl_echo_on (void)
{
    check_for_state();
    return state->gretl_echo;
}

void set_gretl_messages (int e)
{
    check_for_state();
    state->gretl_msgs = e;
}

int gretl_messages_on (void)
{
    check_for_state();
    return state->gretl_msgs;
}

char get_csv_delim (const DATAINFO *pdinfo)
{
    check_for_state();
    if (state->delim > 0) {
	return state->delim;
    } else {
	return pdinfo->delim;
    }
}

int get_long_digits (void)
{
    check_for_state();
    return state->longdigits;
}

const double *get_init_vals (int *pn)
{
    int n;

    if (state->initvals == NULL) {
	*pn = 0;
	return NULL;
    }

    for (n=0; state->initvals[n] != NADBL; n++) {
	if (n > INITVALS_MAX) {
	    *pn = 0;
	    return NULL;
	}
    }

    *pn = n;

    return state->initvals;
}    

int get_hac_lag (int m)
{
    check_for_state();

    /* Variants of Newey-West */

    if (state->ropts.user_lag != 0 && state->ropts.user_lag < m - 2) {
	/* FIXME upper limit? */
	return state->ropts.user_lag;
    }

    if (state->ropts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return 0.75 * pow(m, 1.0 / 3.0);
    } else if (state->ropts.auto_lag == AUTO_LAG_WOOLDRIDGE) {
	return 4.0 * pow(m / 100.0, 2.0 / 9.0);
    }

    /* fallback -- should not be reached */
    return 0.75 * pow(m, 1.0 / 3.0);
}

static char *get_hac_lag_string (void)
{
    check_for_state();

    if (state->ropts.user_lag > 0 && state->ropts.user_lag < 1000) {
	static char lagstr[6];

	sprintf(lagstr, "%d", state->ropts.user_lag);
	return lagstr;
    } else if (state->ropts.auto_lag == AUTO_LAG_STOCK_WATSON) {
	return "nw1";
    } else {
	return "nm2";
    }
}

static int parse_hc_variant (const char *s)
{
    int err = 1;

    check_for_state();

    if (!strcmp(s, "0") || !strcmp(s, "1") ||
	!strcmp(s, "2") || !strcmp(s, "3")) {
	state->ropts.hc_version = atoi(s);
	err = 0;
    } else if (!strcmp(s, "3a")) {
	state->ropts.hc_version = 4;
	err = 0;
    }

    if (err) {
	int hcv;

	if (!strcmp(s, "hc3a")) {
	    state->ropts.hc_version = 4;
	    err = 0;
	} else if (sscanf(s, "hc%d", &hcv)) {
	    if (hcv >= 0 && hcv <= 4) {
		state->ropts.hc_version = hcv;
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

static int set_bkbp_periods (const char *s0, const char *s1,
			     PRN *prn)
{
    int p0, p1;

    if (!isdigit((unsigned char) *s0) || !isdigit((unsigned char) *s1)) {
	return 1;
    }

    p0 = atoi(s0);
    p1 = atoi(s1);

    if (p1 < p0) {
	/* 2nd entry should be bigger than 1st one */
	int tmp = p1;

	p1 = p0;
	p0 = tmp;
    }

    pprintf(prn, _("Baxter-King band = %d-%d periods\n"), p0, p1);

    state->bkopts.periods[0] = p0;
    state->bkopts.periods[1] = p1;

    return 0;
}

static int set_initvals (const char *s)
{
    int i, n;
    double *x = NULL;
    int err = 0;

    /* skip past "set initvals" */
    s += 12;

    n = count_fields(s);

    if (n == 0) {
	free(state->initvals);
	state->initvals = NULL;
	return 0;
    } else if (n > INITVALS_MAX) {
	return E_DATA;
    }

    x = malloc((n + 1) * sizeof *x);
    if (x == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<n; i++) {
	x[i] = gretl_double_from_string(s, &s);
	if (na(x[i])) {
	    err = E_DATA;
	    break;
	}
    }

    if (!err) {
	free(state->initvals);
	x[n] = NADBL; /* sentinel */
	state->initvals = x;
    }

    return err;
}

static int parse_set_plotfile (const char *s)
{
    char *fname;
    int err = 0;

    while (isspace((unsigned char) *s)) {
	s++;
    }

    /* now skip two words, "set" and "plotfile" */
    s += strcspn(s, " ");
    s += strspn(s, " ");
    s += strcspn(s, " ");
    s += strspn(s, " ");

    fname = gretl_strdup(s);
    if (fname != NULL) {
	tailstrip(fname);
	set_gretl_plotfile(fname);
	free(fname);
    } else {
	err = E_ALLOC;
    }

    return err;
}

static char delim_from_arg (const char *s)
{
    char ret = 0;

    if (!strcmp(s, "comma")) {
	ret = ',';
    } else if (!strcmp(s, "space")) {
	ret = ' ';
    } else if (!strcmp(s, "tab")) {
	ret = '\t';
    } else if (!strcmp(s, "semicolon")) {
	ret = ';';
    }

    return ret;
}

static const char *arg_from_delim (char c)
{
    const char *ret = "unset";

    if (c == ',') {
	ret = "comma";
    } else if (c == ' ') {
	ret = "space";
    } else if (c == '\t') {
	ret = "tab";
    } else if (c == ';') {
	ret = "semicolon";
    }

    return ret;
}

static int display_settings (PRN *prn)
{
    double dval;
    unsigned int uval;
    int ival;

    pputs(prn, "Variables that can be set using \"set\" (do \"help set\""
	  " for details):\n");

    pprintf(prn, " echo = %d\n", state->gretl_echo);

    ival = get_use_qr(); /* checks env */
    pprintf(prn, " qr = %d\n", state->use_qr);

    uval = get_gretl_random_seed();
    pprintf(prn, " seed = %u\n", uval);

    pprintf(prn, " hac_lag = %s\n", get_hac_lag_string());
    pprintf(prn, " hc_version = %d\n", state->ropts.hc_version);
    pprintf(prn, " force_hc = %d\n", state->ropts.force_hc);

    pprintf(prn, " garch_vcv = %s\n", garch_vcv_string());

    if (na(state->hp_lambda)) {
	pputs(prn, " hp_lambda: auto\n");
    } else {
	pprintf(prn, " hp_lambda = %g\n", state->hp_lambda);
    }

    pprintf(prn, " bkbp_limits = (%d, %d)\n", state->bkopts.periods[0], 
	    state->bkopts.periods[1]);
    pprintf(prn, " bkbp_k = %d\n", state->bkopts.k);

    if (is_unset(state->horizon)) {
	pputs(prn, " horizon: auto\n");
    } else {
	pprintf(prn, " horizon = %d\n", state->horizon);
    }

    pprintf(prn, " nls_toler = %g\n", get_nls_toler());

    dval = get_bhhh_toler();
    if (na(dval)) {
	pputs(prn, " bhhh_toler = default\n");
    } else {
	pprintf(prn, " bhhh_toler = %g\n", dval);
    }

    pprintf(prn, " bhhh_maxiter = %d\n", get_bhhh_maxiter());
    pprintf(prn, " messages = %d\n", state->gretl_msgs);

    ival =  get_halt_on_error(); /* checks env */
    pprintf(prn, " halt_on_error = %d\n", state->halt_on_error);

    pprintf(prn, " shell_ok = %d\n", state->shell_ok);
    pprintf(prn, " csv_delim = %s\n", arg_from_delim(state->delim));
    pprintf(prn, " longdigits = %s\n", state->longdigits);
    print_initvals(state->initvals, prn);

    return 0;
}

#define boolean_on(s) (!strcmp(s, "on") || !strcmp(s, "1") || \
                       !strcmp(s, "true"))

#define boolean_off(s) (!strcmp(s, "off") || !strcmp(s, "0") || \
                        !strcmp(s, "false"))

int execute_set_line (const char *line, PRN *prn)
{
    char setobj[16], setarg[16], setarg2[16];
    int nw, err = E_PARSE;

    check_for_state();

    *setobj = *setarg = *setarg2 = '\0';

    nw = sscanf(line, "%*s %15s %15s %15s", setobj, setarg, setarg2);

    if (nw <= 0) {
	return display_settings(prn);
    }

    /* specials which need the whole line */
    if (nw > 1) {
	if (!strcmp(setobj, "plotfile")) {
	    return parse_set_plotfile(line);
	} else if (!strcmp(setobj, "initvals")) {
	    return set_initvals(line);
	}
    }

    if (nw == 1) {
	if (!strcmp(setobj, "echo")) {
	    state->gretl_echo = 1;
	    err = 0;
	}
    } else if (nw == 2) {
	lower(setarg);

	/* set command echo on/off */
	if (!strcmp(setobj, "echo")) {
	    if (boolean_off(setarg)) {
		state->gretl_echo = 0;
		err = 0;
	    } else if (boolean_on(setarg)) {
		state->gretl_echo = 1;
		err = 0;
	    }
	} else if (!strcmp(setobj, "messages")) {
	    if (boolean_off(setarg)) {
		state->gretl_msgs = 0;
		err = 0;
	    } else if (boolean_on(setarg)) {
		state->gretl_msgs = 1;
		err = 0;
	    }
	} else if (!strcmp(setobj, "csv_delim")) {
	    char c = delim_from_arg(setarg);

	    if (c > 0) {
		state->delim = c;
		err = 0;
	    }
	} else if (!strcmp(setobj, "hac_lag")) {
	    /* set max lag for HAC estimation */
	    if (!strcmp(setarg, "nw1")) {
		state->ropts.auto_lag = AUTO_LAG_STOCK_WATSON;
		state->ropts.user_lag = 0;
		err = 0;
	    } else if (!strcmp(setarg, "nw2")) {
		state->ropts.auto_lag = AUTO_LAG_WOOLDRIDGE;
		state->ropts.user_lag = 0;
		err = 0;
	    } else if (isdigit(*setarg)) {
		state->ropts.user_lag = atoi(setarg);
		err = 0;
	    }
	} else if (!strcmp(setobj, "hc_version")) {
	    /* set HCCM variant */
	    err = parse_hc_variant(setarg);
	} else if (!strcmp(setobj, "force_hc")) {
	    /* use HCCM, not HAC, even for time series */
	    if (boolean_on(setarg)) { 
		set_force_hc(1);
		err = 0;
	    } else if (boolean_off(setarg)) { 
		set_force_hc(0);
		err = 0;
	    }
	} else if (!strcmp(setobj, "garch_vcv")) {
	    /* set GARCH VCV variant */
	    err = set_garch_vcv_variant(setarg);
	} else if (!strcmp(setobj, "qr")) {
	    /* switch QR vs Cholesky decomposition */
	    if (boolean_on(setarg)) {
		state->use_qr = 1;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->use_qr = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "halt_on_error")) {
	    if (boolean_on(setarg)) {
		state->halt_on_error = 1;
		err = 0;
	    } else if (boolean_off(setarg)) {
		state->halt_on_error = 0;
		err = 0;
	    }
	} else if (!strcmp(setobj, "shell_ok")) {
	    pprintf(prn, "You can only set this variable via the gretl GUI\n");
	} else if (!strcmp(setobj, "seed")) {
	    /* seed for PRNG */
	    if (isdigit(*setarg)) {
		int k = atoi(setarg);

		gretl_rand_set_seed((unsigned int) k);
		pprintf(prn, 
			_("Pseudo-random number generator seeded with %d\n"), k);
		state->seed = k;
		err = 0;

	    }
	} else if (!strcmp(setobj, "hp_lambda")) {
	    /* Hodrick-Prescott filter parameter */
	    if (!strcmp(setarg, "auto")) {
		state->hp_lambda = NADBL;
		err = 0;
	    } else if (check_atof(setarg) == 0) {
		state->hp_lambda = atof(setarg);
		err = 0;
	    }
	} else if (!strcmp(setobj, "bkbp_k")) {
	    /* Baxter-King approximation order */
	    if (isdigit(*setarg)) {
		state->bkopts.k = atoi(setarg);
		pprintf(prn, 
			_("Baxter-King approximation = %d\n"), state->bkopts.k);
		err = 0;
	    }
	} else if (!strcmp(setobj, "horizon")) {
	    /* horizon for VAR impulse responses */
	    if (!strcmp(setarg, "auto")) {
		state->horizon = UNSET_INT;
		err = 0;
	    } else {
		state->horizon = atoi(setarg);
		if (state->horizon > 0) {
		    err = 0;
		} else {
		    state->horizon = UNSET_INT;
		}
	    }	    
	} else if (!strcmp(setobj, "nls_toler")) {
	    double tol;

	    if (sscanf(setarg, "%lf", &tol)) {
		err = set_nls_toler(tol);
	    }
	} else if (!strcmp(setobj, "bhhh_toler")) {
	    /* Tolerance for BHHH (ARMA, Tobit) */
	    double tol;

	    if (!strcmp(setarg, "default")) {
		set_bhhh_toler(NADBL);
		err = 0;
	    } else if (sscanf(setarg, "%lf", &tol)) {
		err = set_bhhh_toler(tol);
	    }
	} else if (!strcmp(setobj, "bhhh_maxiter")) {
	    /* Maximum iterations for BHHH (ARMA, Tobit) */
	    if (isdigit(*setarg)) {
		err = set_bhhh_maxiter(atoi(setarg));
	    }
	} else if (!strcmp(setobj, "longdigits")) {
	    if (isdigit(*setarg)) {
		err = set_long_digits(atoi(setarg));
	    }
	} 
    } else if (nw == 3) {
	if (!strcmp(setobj, "bkbp_limits")) {
	    err = set_bkbp_periods(setarg, setarg2, prn);
	}
    }
		    
    return err;
}

/* use Cholesky or QR for regression? */

void set_use_qr (int set)
{
    check_for_state();
    state->use_qr = set;
}

int get_use_qr (void)
{
    check_for_state();

    if (is_unset(state->use_qr)) {
	char *s = getenv("GRETL_USE_QR");

	if (s != NULL && *s != '\0' && *s != '0') {
	    state->use_qr = 1;
	} else {
	    state->use_qr = 0;
	}
    } 

    return state->use_qr;
}

int get_halt_on_error (void)
{
    check_for_state();

    if (is_unset(state->halt_on_error)) {
	char *s = getenv("GRETL_KEEP_GOING");

	if (s != NULL && *s != '\0' && *s != '0') {
	    state->halt_on_error = 0;
	} else {
	    state->halt_on_error = 1;
	}
    } 

    return state->halt_on_error;
}

#ifndef WIN32
static int read_cli_shell_status (void)
{
    char shellstamp[FILENAME_MAX];
    FILE *fp;
    int ok = 0;

    sprintf(shellstamp, "%s.gretl_shell_stamp", gretl_user_dir());
    fp = fopen(shellstamp, "r");
    if (fp != NULL) {
	ok = 1;
	fclose(fp);
    }

    return ok;
}
#endif

void set_shell_ok (int set)
{
    check_for_state();
    state->shell_ok = set;
}

int get_shell_ok (void)
{
    check_for_state();

#ifndef WIN32
    if (!gretl_in_gui_mode()) {
	state->shell_ok = read_cli_shell_status();
    }
#endif

    return state->shell_ok;
}

static void update_sample_info (set_vars *sv, const DATAINFO *pdinfo)
{
    if (pdinfo != NULL) {
	sv->sinfo.t1 = pdinfo->t1;
	sv->sinfo.t2 = pdinfo->t2;
	sv->sinfo.submode = pdinfo->submode;
	free(sv->sinfo.submask);
	sv->sinfo.submask = copy_datainfo_submask(pdinfo);
    }	 
}

/* Mechanism for pushing and popping program state for user-defined
   functions. push_program_state() is used when a function starts
   execution: the function gets a copy of the current program state,
   while that state is pushed onto the stack for restoration when the
   function exits.
*/

int n_states;
static set_vars **state_stack;

int push_program_state (const DATAINFO *pdinfo)
{
    set_vars **sstack;
    set_vars *newstate;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "push_program_state: ns=%d\n", ns);
    if (pdinfo != NULL) {
	fprintf(stderr, " pdinfo=%p, t1=%d, t2=%d\n",
		(void *) pdinfo, pdinfo->t1, pdinfo->t2);
    } else {
	fprintf(stderr, " pdinfo is NULL\n");
    }
#endif

    newstate = malloc(sizeof *newstate);

    if (newstate == NULL) {
	err = 1;
    } else {
	sstack = realloc(state_stack, (ns + 1) * sizeof *sstack);
	if (sstack == NULL) {
	    free(newstate);
	    err = 1;
	}
    }

    if (!err) {
	if (ns == 0) {
	    /* set all defaults */
	    state_vars_init(newstate);
	} else {
	    /* set sample on existing state */
	    update_sample_info(state, pdinfo);
	    /* copy existing state */
	    state_vars_copy(newstate, pdinfo);
	}
	state_stack = sstack;
	state = state_stack[ns] = newstate;
	n_states++;
    }

#if PDEBUG
    if (!err) {
	fprintf(stderr, " state is now state_stack[%d]\n", ns);
    }
#endif

    return err;
}

static int 
subsampled_at_this_level (const DATAINFO *pdinfo, set_vars *sv)
{
    int ret = 0;

    if (!complex_subsampled()) {
	return 0;
    }

    if (pdinfo->submask != NULL) {
	if (sv->sinfo.submask == NULL) {
	    ret = 1;
	} else if (submask_cmp(pdinfo->submask, sv->sinfo.submask)) {
	    ret = 1;
	}
    } 

    return ret;
}

static void free_state (set_vars *sv)
{
    if (sv->sinfo.submask != NULL) {
	free(sv->sinfo.submask);
    }

    if (sv->initvals != NULL) {
	free(sv->initvals);
    }

    free(sv);
}

/* Called when a new-style function exits: restores the program state
   that was in force when the function started executing.  This may
   involve putting the sample back the way it was.  But if the dataset
   was complex subsampled either before or after function execution,
   we will not (currently) mess with the sample state.  That case
   requires more careful thought.
*/

int pop_program_state (double ***pZ, DATAINFO **ppdinfo)
{
    set_vars **sstack;
    int ns = n_states;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "pop_program_state called: ns=%d, pdinfo=%p\n",
	    ns, (void *) *ppdinfo);
#endif

    if (ns < 2) {
	err = 1;
    } else {
	free_state(state_stack[ns - 1]);
	state_stack[ns - 1] = NULL;
	sstack = realloc(state_stack, (ns - 1) * sizeof *sstack);
	if (sstack == NULL) {
	    err = 1;
	}	
    }

    if (!err) {
	state_stack = sstack;
	state = state_stack[ns - 2];
	n_states--;
    }

#if PDEBUG
    fprintf(stderr, " state is now state_stack[%d]\n", ns - 2);
#endif

    if (!err && *ppdinfo != NULL && sinfo_is_set(state->sinfo)) {
	if (subsampled_at_this_level(*ppdinfo, state)) {
#if PDEBUG
	    fprintf(stderr, " undoing current sub-sampling\n");
#endif
	    restore_full_sample(pZ, ppdinfo);
	    if (state->sinfo.submask != NULL) {
#if PDEBUG
		fprintf(stderr, " reimposing outer sample restriction\n");
#endif
		restrict_sample_from_mask(state->sinfo.submask, state->sinfo.submode, 
					  pZ, ppdinfo);
	    }
	}
	(*ppdinfo)->t1 = state->sinfo.t1;
	(*ppdinfo)->t2 = state->sinfo.t2;
    }

#if PDEBUG
    if (*ppdinfo != NULL) {
	fprintf(stderr, " err=%d, pdinfo->t1=%d, pdinfo->t2=%d\n",
		err, (*ppdinfo)->t1, (*ppdinfo)->t2);
    } 
#endif

    return err;
}

/* initialization of all user-settable settings */

int libset_init (void)
{
    static int done;
    int err = 0;

#if PDEBUG
    fprintf(stderr, "libset_init called, done=%d\n", done);
#endif

    if (!done) {
	err = push_program_state(NULL);
	done = 1;
    }

    return err;
}

void libset_cleanup (void)
{
    int i;

    for (i=0; i<n_states; i++) {
	free_state(state_stack[i]);
    }

    free(state_stack);
    state_stack = NULL;
    n_states = 0;
}

int libset_restore_state_zero (double ***pZ, DATAINFO **ppdinfo)
{
    int i, ns = n_states;
    int err = 0;

    for (i=ns; i>=2 && !err; i--) {
	err = pop_program_state(pZ, ppdinfo);
    }

    return err;
}

/* switches for looping, batch mode, and pausing between screens of
   output: these depend on the state of the program calling libgretl,
   but are not user-settable
*/

static int gretl_text_pause;
static int loop_on;
static int batch_mode;
static int gui_mode;

int gretl_get_text_pause (void)
{
    return gretl_text_pause;
}

void set_loop_on (void)
{
    loop_on = 1;
    gretl_text_pause = 0;
}

void set_loop_off (void)
{
    loop_on = 0;
    if (!batch_mode && !gui_mode) {
	gretl_text_pause = 1;
    }
}

int gretl_looping (void)
{
    return loop_on;
}

void gretl_set_batch_mode (int b)
{
    batch_mode = b;
    if (batch_mode) {
	gretl_text_pause = 0;
    }	
}

int gretl_in_batch_mode (void)
{
    return batch_mode;
}

void gretl_set_gui_mode (int g)
{
    gui_mode = g;
}

int gretl_in_gui_mode (void)
{
    return gui_mode;
}


