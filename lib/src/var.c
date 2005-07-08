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
 *   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111, USA.
 *
 */

/* var.c - vector autoregressions */  

#include "libgretl.h" 
#include "var.h"  
#include "libset.h"

/* in transforms.c */
extern int 
real_list_laggenr (const int *list, double ***pZ, DATAINFO *pdinfo,
		   int maxlag, int **lagnums);

static gretl_matrix *
gretl_VAR_get_fcast_decomp (GRETL_VAR *var, int targ, int periods);

static gretl_matrix *irf_bootstrap (const GRETL_VAR *var, 
				    int targ, int shock, int periods,
				    const double **Z, 
				    const DATAINFO *pdinfo);

#define VAR_DEBUG 0

struct GRETL_VAR_ {
    int neqns;         /* number of equations in system */
    int order;         /* lag order */
    int n;             /* number of observations */
    int ifc;           /* equations include a constant (1) or not (0) */
    gretl_matrix *A;   /* augmented coefficient matrix */
    gretl_matrix *E;   /* residuals matrix */
    gretl_matrix *C;   /* augmented Cholesky-decomposed error matrix */
    gretl_matrix *F;   /* optional forecast matrix */
    MODEL **models;    /* pointers to individual equation estimates */
    double *Fvals;     /* hold results of F-tests */
    char *name;        /* for use in session management */
};

struct var_resids {
    int *levels_list;
    double **uhat;
    int m;
    int t1, t2;
};

struct var_lists {
    int *detvars;
    int *stochvars;
    int *reglist;
    int *testlist;
    int **lagvlist;
};

enum {
    ADF_EG_TEST   = 1 << 0,
    ADF_PRINT_ACK = 1 << 1
} adf_flags;

#define TREND_FAILED 9999

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags, 
			  PRN *prn);

static void pad_var_coeff_matrix (gretl_matrix *A, int neqns, int order)
{
    int i, j;
    int rowmax = neqns * order;

    for (i=neqns; i<rowmax; i++) {
	for (j=0; j<rowmax; j++) {
	    gretl_matrix_set(A, i, j, (j == i - neqns)? 1.0 : 0.0);
	}
    }
}

static GRETL_VAR *gretl_VAR_new (int neqns, int order, const DATAINFO *pdinfo)
{
    GRETL_VAR *var;
    int rows;
    int i, j;
    int err = 0;

    if (neqns == 0 || order == 0) {
	return NULL;
    }

    var = malloc(sizeof *var);
    if (var == NULL) {
	return NULL;
    }

    var->neqns = neqns;
    var->order = order;

    var->A = NULL;
    var->E = NULL;
    var->C = NULL;
    var->F = NULL;

    var->models = NULL;
    var->Fvals = NULL;
    var->name = NULL;

    rows = neqns * order;

    var->A = gretl_matrix_alloc(rows, rows);
    if (var->A == NULL) {
	err = 1;
    } else {
	pad_var_coeff_matrix(var->A, var->neqns, var->order);
    }

    if (!err) {
	var->C = gretl_matrix_alloc(rows, neqns);
	if (var->C == NULL) {
	    err = 1;
	    gretl_matrix_free(var->A);
	    var->A = NULL;
	} else {
	    gretl_matrix_zero(var->C);
	}
    }

    if (!err) {
	var->models = malloc(neqns * sizeof *var->models);
	if (var->models == NULL) {
	    err = 1;
	} else {
	    for (i=0; i<neqns; i++) {
		var->models[i] = gretl_model_new();
		if (var->models[i] == NULL) {
		    err = 1;
		    for (j=0; j<i; j++) {
			free(var->models[i]);
		    }
		    free(var->models);
		    var->models = NULL;
		}
	    }
	}
    } 

    if (!err) {
	int m = neqns * neqns + neqns;
	
	var->Fvals = malloc(m  * sizeof *var->Fvals);
	if (var->Fvals == NULL) {
	    err = 1;
	}
    }

    if (err) {
	gretl_VAR_free(var);
	var = NULL;
    }

    return var;
}

void gretl_VAR_free (GRETL_VAR *var)
{
    int i;

    if (var == NULL) return;

    gretl_matrix_free(var->A);
    gretl_matrix_free(var->E);
    gretl_matrix_free(var->C);
    gretl_matrix_free(var->F);

    free(var->Fvals);
    free(var->name);

    if (var->models != NULL) {
	for (i=0; i<var->neqns; i++) {
	    clear_model(var->models[i]);
	    free(var->models[i]);
	}
	free(var->models);
    }

    free(var);
}

void gretl_VAR_free_unnamed (GRETL_VAR *var)
{
    if (var == NULL) return;

    if (var->name == NULL || *var->name == '\0') {
	gretl_VAR_free(var);
    }
}

static int
gretl_VAR_add_forecast (GRETL_VAR *var, int t1, int t2, const double **Z, 
			const DATAINFO *pdinfo, gretlopt opt)
{
    const MODEL *pmod;
    gretl_matrix *F;
    double fti, xti;
    int i, j, k, s, t;
    int nf, ns, lag, vj, m;
    int staticfc, fcols;

    pmod = var->models[0];

    nf = t2 - t1 + 1;

    staticfc = (opt & OPT_S);
    if (staticfc) {
	fcols = var->neqns;
    } else {
	fcols = 2 * var->neqns;
    }

    /* rows = number of forecast periods; cols = 1 to hold forecast
       for each variable, plus 1 to hold variance for each variable
       if forecast is dynamic.
    */
    F = gretl_matrix_alloc(nf, fcols);
    if (F == NULL) {
	return E_ALLOC;
    }

    gretl_matrix_zero(F);

    ns = var->order * var->neqns;

    for (t=t1; t<=t2; t++) {
	int miss = 0;

	s = t - t1;
	for (i=0; i<var->neqns; i++) {
	    pmod = var->models[i];
	    fti = 0.0;
	    lag = 1;
	    k = 0;
	    for (j=0; j<pmod->ncoeff; j++) {
		vj = pmod->list[j + 2];
		if (j < ns + pmod->ifc && vj > 0) {
		    /* stochastic var */
		    if (staticfc || s - lag < 0) {
			/* pre-forecast value */
			m = (j - pmod->ifc) / var->order;
			vj = var->models[m]->list[1];
			if (t - lag < 0) {
			    xti = NADBL;
			} else {
			    xti = Z[vj][t-lag];
			}
			if (na(xti)) {
			    miss = 1;
			}
		    } else {
			/* prior forecast value */
			xti = gretl_matrix_get(F, s - lag, k);
		    }
		    lag++;
		    if (lag > var->order) {
			lag = 1;
			k++;
		    }
		} else {
		    /* deterministic var: value from dataset */
		    xti = Z[vj][t];
		    if (na(xti)) {
			miss = 1;
		    }
		}
		if (miss) {
		    fti = NADBL;
		} else {
		    fti += pmod->coeff[j] * xti;
		}
	    }
	    gretl_matrix_set(F, s, i, fti);
	}
    }

    /* now get variances, if not static */
    if (!staticfc) {
	double vti;
	int totcol;

	for (i=0; i<var->neqns; i++) {
	    gretl_matrix *vd;
	    vd = gretl_VAR_get_fcast_decomp(var, i, nf);
	    if (vd != NULL) {
		totcol = gretl_matrix_cols(vd) - 1;
		for (s=0; s<nf; s++) {
		    vti = gretl_matrix_get(vd, s, totcol);
		    gretl_matrix_set(F, s, var->neqns + i, vti);
		}
		gretl_matrix_free(vd);
	    } else {
		for (s=0; s<nf; s++) {
		    gretl_matrix_set(F, s, var->neqns + i, NADBL);
		}
	    }
	}
    }

#if 0
    gretl_matrix_print(F, "var->F", NULL);
#endif

    gretl_matrix_set_int(F, t1);

    var->F = F;

    return 0;
}

const gretl_matrix *
gretl_VAR_get_forecast_matrix (GRETL_VAR *var, int t1, int t2, const double **Z, 
			       const DATAINFO *pdinfo, gretlopt opt)
{
    if (var->F != NULL) {
	int ncols, nf = t2 - t1 + 1;
	int ft1 = gretl_matrix_get_int(var->F);

	ncols = (opt & OPT_S)? var->neqns: 2 * var->neqns;

	if (nf == gretl_matrix_rows(var->F) && t1 == ft1 && 
	    ncols == gretl_matrix_cols(var->F)) {
	    ; /* already done, fine */
	} else {
	    gretl_matrix_free(var->F);
	    var->F = NULL;
	}
    }

    if (var->F == NULL) {
	gretl_VAR_add_forecast(var, t1, t2, Z, pdinfo, opt);
    }

    return var->F;
}

int gretl_VAR_print_VCV (const GRETL_VAR *var, PRN *prn)
{
    gretl_matrix *V;
    double ldet;
    int err = 0;

    if (var->E == NULL) {
	err = 1;
    } else {
	V = gretl_matrix_vcv(var->E);
	if (V == NULL) {
	    err = 1;
	} else {
	    ldet = print_contemp_covariance_matrix(V, prn);
	    if (na(ldet)) {
		err = 1;
	    }
	    gretl_matrix_free(V);
	}
    }

    return err;
}

static int gretl_VAR_do_error_decomp (int n, int neqns,
				      const gretl_matrix *E,
				      gretl_matrix *C)
{
    gretl_matrix *tmp = NULL;
    int i, j, err = 0;

    tmp = gretl_matrix_alloc(neqns, neqns);
    if (tmp == NULL) {
	err = E_ALLOC;
    }

    /* form e'e */
    if (!err && gretl_matrix_multiply_mod (E, GRETL_MOD_TRANSPOSE,
					   E, GRETL_MOD_NONE,
					   tmp)) {
	err = 1;
    }

    /* divide by T (or use df correction?) to get sigma-hat.
       Note: RATS 4 uses straight T.
    */
    if (!err) {
	gretl_matrix_divide_by_scalar(tmp, (double) n);
    }

#if VAR_DEBUG
    if (!err) {
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	gretl_matrix_print(tmp, "Sigma-hat from VAR system", prn);
	gretl_print_destroy(prn);
    }
#endif

    /* lower-triangularize and decompose */
    if (!err) {
	for (i=0; i<neqns-1; i++) {
	    for (j=i+1; j<neqns; j++) {
		gretl_matrix_set(tmp, i, j, 0.0);
	    }
	}
	err = gretl_matrix_cholesky_decomp(tmp);
    }

    /* write the decomposition into the C matrix */
    if (!err) {
	for (i=0; i<neqns; i++) {
	    for (j=0; j<neqns; j++) {
		double x = gretl_matrix_get(tmp, i, j);

		gretl_matrix_set(C, i, j, x);
	    }
	}
    }

    if (tmp != NULL) {
	gretl_matrix_free(tmp);
    }

    return err;
}

int gretl_VAR_get_variable_number (const GRETL_VAR *var, int k)
{
    return (var->models[k])->list[1];
}

int gretl_VAR_get_n_equations (const GRETL_VAR *var)
{
    return var->neqns;
}

int gretl_VAR_get_t1 (const GRETL_VAR *var)
{
    return var->models[0]->t1;
}

int gretl_VAR_get_t2 (const GRETL_VAR *var)
{
    return var->models[0]->t2;
}

const MODEL *gretl_VAR_get_model (const GRETL_VAR *var, int i)
{
    if (i < var->neqns) {
	return var->models[i];
    } else {
	return NULL;
    }
}

#define VARS_IN_ROW 4

static void tex_print_double (double x, PRN *prn)
{
    char number[16];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (x < 0.) {
	pprintf(prn, "$-$%s", number + 1);
    } else {
	pputs(prn, number);
    }
}

static int periods_from_pd (int pd)
{
    int periods = 10;

    if (pd == 4) {
	/* quarterly: try 5 years */
	periods = 20;
    } else if (pd == 12) {
	/* monthly: two years */
	periods = 24;
    } else if (pd == 7 || pd == 6 || pd == 5) {
	/* daily: three weeks */
	periods = 3 * pd;
    } 

    return periods;
}

int default_VAR_horizon (const DATAINFO *pdinfo)
{
    int h = get_VAR_horizon();

    if (h <= 0) {
	h = periods_from_pd(pdinfo->pd);
    }

    return h;
}

static int 
gretl_VAR_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn)
{
    int i, t;
    int vsrc;
    int rows = var->neqns * var->order;
    gretl_matrix *rtmp, *ctmp;
    int block, blockmax;
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return 1;
    }  

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	return E_ALLOC;
    }

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	gretl_matrix_free(rtmp);
	return E_ALLOC;
    }

    vsrc = (var->models[shock])->list[1];

    blockmax = var->neqns / VARS_IN_ROW;
    if (var->neqns % VARS_IN_ROW) {
	blockmax++;
    }

    for (block=0; block<blockmax && !err; block++) {
	int vtarg, k;
	char vname[16];
	double r;

	if (tex_format(prn)) {
	    pputs(prn, "\\vspace{1em}\n\n");
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    tex_escape(vname, pdinfo->varname[vsrc]));

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", I_("continued"));
	    }
	    pputs(prn, "\\vspace{1em}\n\n"
		  "\\begin{longtable}{rcccc}\n");
	} else {
	    pprintf(prn, _("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[vsrc]);

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", _("continued"));
	    }
	}

	if (tex_format(prn)) {
	    pprintf(prn, "%s & ", I_("period"));
	} else {
	    pprintf(prn, "%s ", _("period"));
	}

	for (i=0; i<VARS_IN_ROW; i++) {
	    k = VARS_IN_ROW * block + i;
	    if (k >= var->neqns) {
		break;
	    }
	    vtarg = (var->models[k])->list[1];
	    if (tex_format(prn)) {
		pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[vtarg]));
		if (i < VARS_IN_ROW - 1 && k < var->neqns - 1) {
		    pputs(prn, "& ");
		} else {
		    pputs(prn, "\\\\");
		}
	    } else {
		pprintf(prn, "  %8s  ", pdinfo->varname[vtarg]);
	    }
	}

	pputs(prn, "\n\n");

	for (t=0; t<periods && !err; t++) {
	    pprintf(prn, " %3d  ", t + 1);
	    if (tex_format(prn)) {
		pputs(prn, "& ");
	    }
	    if (t == 0) {
		/* calculate initial estimated responses */
		err = gretl_matrix_copy_values(rtmp, var->C);
	    } else {
		/* calculate further estimated responses */
		err = gretl_matrix_multiply(var->A, rtmp, ctmp);
		gretl_matrix_copy_values(rtmp, ctmp);
	    }

	    if (err) break;

	    /* matrix rtmp holds the responses */

	    for (i=0; i<VARS_IN_ROW; i++) {
		k = VARS_IN_ROW * block + i;
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(rtmp, k, shock);
		if (tex_format(prn)) {
		    tex_print_double(r, prn);
		    if (i < VARS_IN_ROW - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else {
		    pprintf(prn, "%#12.5g ", r);
		}
	    }
	    if (tex_format(prn)) {
		pputs(prn, "\\\\\n");
	    } else {
		pputc(prn, '\n');
	    }
	}

	if (tex_format(prn)) {
	    pputs(prn, "\\end{longtable}\n\n");
	} else {
	    pputc(prn, '\n');
	}

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    return err;
}

static gretl_matrix *
gretl_VAR_get_point_responses (GRETL_VAR *var, int targ, int shock,
			       int periods) 
{
    int rows = var->neqns * var->order;
    gretl_matrix *rtmp = NULL;
    gretl_matrix *ctmp = NULL;
    gretl_matrix *resp = NULL;
    double rt;
    int t, err = 0;

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return NULL;
    }  

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return NULL;
    } 

    if (periods <= 0) {
	fprintf(stderr, "Invalid number of periods\n");
	return NULL;
    }

    resp = gretl_matrix_alloc(periods, 1);
    if (resp == NULL) {
	return NULL;
    }

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	gretl_matrix_free(resp);
	return NULL;
    }

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	free(resp);
	gretl_matrix_free(rtmp);
	return NULL;
    }

    for (t=0; t<periods && !err; t++) {
	if (t == 0) {
	    /* calculate initial estimated responses */
	    err = gretl_matrix_copy_values(rtmp, var->C);
	} else {
	    /* calculate further estimated responses */
	    err = gretl_matrix_multiply(var->A, rtmp, ctmp);
	    gretl_matrix_copy_values(rtmp, ctmp);
	}

	if (!err) {
	    rt = gretl_matrix_get(rtmp, targ, shock);
	    gretl_matrix_set(resp, t, 0, rt);
	}
    }

    gretl_matrix_free(rtmp);
    gretl_matrix_free(ctmp);

    return resp;    
}

gretl_matrix *
gretl_VAR_get_impulse_response (GRETL_VAR *var, 
				int targ, int shock, int periods,
				const double **Z,
				const DATAINFO *pdinfo)
{
    gretl_matrix *point = NULL;
    gretl_matrix *full = NULL;
    gretl_matrix *ret = NULL;
    int i;

    point = gretl_VAR_get_point_responses(var, targ, shock, periods);

    if (Z == NULL) {
	/* no data matrix given: just return point estimate */
	ret = point;
    } else if (point != NULL) {
	full = irf_bootstrap(var, targ, shock, periods, Z, pdinfo);
	if (full != NULL) {
	    double p;

	    for (i=0; i<periods; i++) {
		p = gretl_matrix_get(point, i, 0);
		gretl_matrix_set(full, i, 0, p);
	    }
	}
	gretl_matrix_free(point);
	ret = full;
    }

    return ret;
}

static gretl_matrix *
gretl_VAR_get_fcast_decomp (GRETL_VAR *var, int targ, int periods) 
{
    int i, t;
    int rows = var->neqns * var->order;
    gretl_matrix *ctmp = NULL, *idx = NULL, *vtmp = NULL;
    gretl_matrix *cic = NULL, *vt = NULL;
    gretl_matrix *vd = NULL;
    int err = 0;

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return NULL;
    } 

    if (periods <= 0) {
	fprintf(stderr, "Invalid number of periods\n");
	return NULL;
    }

    vd = gretl_matrix_alloc(periods, var->neqns + 1);
    ctmp = gretl_matrix_alloc(var->neqns, rows);
    idx = gretl_matrix_alloc(var->neqns, var->neqns); 
    cic = gretl_matrix_alloc(rows, rows);
    vt = gretl_matrix_alloc(rows, rows);
    vtmp = gretl_matrix_alloc(rows, rows);

    if (vd == NULL || ctmp == NULL || idx == NULL ||
	cic == NULL || vt == NULL || vtmp == NULL) {
	gretl_matrix_free(vd);
	gretl_matrix_free(ctmp);
	gretl_matrix_free(idx);
	gretl_matrix_free(cic);
	gretl_matrix_free(vt);
	gretl_matrix_free(vtmp);
	return NULL;
    }

    for (i=0; i<var->neqns; i++) {
	double vti;

	/* make appropriate index matrix */
	gretl_matrix_zero(idx);
	gretl_matrix_set(idx, i, i, 1.0);

	for (t=0; t<periods && !err; t++) {

	    if (t == 0) {
		/* calculate initial variances */
		err = gretl_matrix_multiply_mod(idx, GRETL_MOD_NONE,
						var->C, GRETL_MOD_TRANSPOSE,
						ctmp);
		err = gretl_matrix_multiply(var->C, ctmp, cic);
		gretl_matrix_copy_values(vt, cic);
	    } else {
		/* calculate further variances */
		err = gretl_matrix_multiply_mod(vt, GRETL_MOD_NONE,
						var->A, GRETL_MOD_TRANSPOSE,
						vtmp);
		err = gretl_matrix_multiply(var->A, vtmp, vt);
		gretl_matrix_add_to(vt, cic);
	    }

	    if (err) break;

	    vti = gretl_matrix_get(vt, targ, targ);
	    gretl_matrix_set(vd, t, i, vti);
	}
    }

    /* normalize variance contributions as percentage shares */
    for (t=0; t<periods && !err; t++) {
	double vtot = 0.0;
	double vi;

	for (i=0; i<var->neqns; i++) {
	    vtot += gretl_matrix_get(vd, t, i);
	}

	for (i=0; i<var->neqns; i++) {
	    vi = gretl_matrix_get(vd, t, i);
	    gretl_matrix_set(vd, t, i, 100.0 * vi / vtot);
	}

	gretl_matrix_set(vd, t, var->neqns, sqrt(vtot));
    }

    gretl_matrix_free(ctmp);
    gretl_matrix_free(idx);
    gretl_matrix_free(cic);
    gretl_matrix_free(vt);
    gretl_matrix_free(vtmp);

    return vd;
}

#define VDROWMAX 5

static int 
gretl_VAR_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn)
{
    int i, t;
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int err = 0;

    if (prn == NULL) {
	return 0;
    }

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return 1;
    } 

    vd = gretl_VAR_get_fcast_decomp(var, targ, periods);
    if (vd == NULL) {
	return E_ALLOC;
    }

    vtarg = (var->models[targ])->list[1];

    blockmax = (var->neqns + 1) / VDROWMAX;
    if ((var->neqns + 1) % VDROWMAX) {
	blockmax++;
    }

    for (block=0; block<blockmax; block++) {
	int k, vsrc;
	char vname[16];
	double r;

	/* print block header */
	if (tex_format(prn)) {
	    pputs(prn, "\\vspace{1em}\n\n");
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    tex_escape(vname, pdinfo->varname[vtarg]));

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", I_("continued"));
	    }
	    pputs(prn, "\\vspace{1em}\n\n"
		  "\\begin{longtable}{rccccc}\n");
	} else {
	    pprintf(prn, _("Decomposition of variance for %s"), 
		    pdinfo->varname[vtarg]);

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", _("continued"));
	    }
	}

	/* first column: print period/step label */
	if (tex_format(prn)) {
	    pprintf(prn, "%s & ", I_("period"));
	} else {
	    pprintf(prn, "%s ", _("period"));
	}

	/* print variable names row */
	for (i=0; i<VDROWMAX; i++) {
	    k = VDROWMAX * block + i - 1;
	    if (k < 0) {
		if (tex_format(prn)) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else {
		    pprintf(prn, " %12s ", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) {
		break;
	    }
	    vsrc = (var->models[k])->list[1];
	    if (tex_format(prn)) {
		pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[vsrc]));
		if (i < VDROWMAX - 1 && k < var->neqns - 1) pputs(prn, "& ");
		else pputs(prn, "\\\\");
	    } else {
		pprintf(prn, "  %8s ", pdinfo->varname[vsrc]);
	    }
	}

	pputs(prn, "\n\n");

	/* print block of numbers */
	for (t=0; t<periods && !err; t++) {
	    pprintf(prn, " %3d  ", t + 1);
	    if (tex_format(prn)) pputs(prn, "& ");

	    for (i=0; i<VDROWMAX; i++) {
		k = VDROWMAX * block + i - 1;
		if (k < 0) {
		    r = gretl_matrix_get(vd, t, var->neqns);
		    if (tex_format(prn)) {
			pprintf(prn, "%g & ", r);
		    } else {
			pprintf(prn, " %14g ", r);
		    }
		    continue;
		}
		if (k >= var->neqns) {
		    break;
		}
		r = gretl_matrix_get(vd, t, k);
		if (tex_format(prn)) {
		    pprintf(prn, "$%.4f$", r);
		    if (i < VDROWMAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else {
		    pprintf(prn, "%10.4f ", r);
		}
	    }
	    if (tex_format(prn)) {
		pputs(prn, "\\\\\n");
	    } else {
		pputc(prn, '\n');
	    }
	}

	if (tex_format(prn)) {
	    pputs(prn, "\\end{longtable}\n\n");
	} else {
	    pputc(prn, '\n');
	}

	if (pause && block < blockmax - 1) {
	    scroll_pause();
	}
    }

    if (vd != NULL) {
	gretl_matrix_free(vd);
    }

    return err;
}

static int gettrend (double ***pZ, DATAINFO *pdinfo, int square)
{
    int index;
    int t, n = pdinfo->n, v = pdinfo->v;
    double x;

    if (square) {
	index = varindex(pdinfo, "timesq");
    } else {
	index = varindex(pdinfo, "time");
    }

    if (index < v) {
	return index;
    }
    
    if (dataset_add_series(1, pZ, pdinfo)) {
	return TREND_FAILED;
    }

    for (t=0; t<n; t++) {
	x = (double) t + 1;
	(*pZ)[v][t] = (square)? x * x : x;
    }

    if (square) {
	strcpy(pdinfo->varname[v], "timesq");
	strcpy(VARLABEL(pdinfo, v), _("squared time trend variable"));
    } else {
	strcpy(pdinfo->varname[v], "time");
	strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
    }
	    
    return index;
}

static void var_lists_free (struct var_lists *vl)
{
    if (vl->lagvlist != NULL && vl->stochvars != NULL) {
	int i, ns = vl->stochvars[0];

	for (i=0; i<ns; i++) {
	    free(vl->lagvlist[i]);
	}
	free(vl->lagvlist);
    }

    free(vl->detvars);
    free(vl->stochvars);
    free(vl->reglist);
    free(vl->testlist);
}

static int **lagvlist_construct (int nstoch, int order)
{
    int **lvlist;
    int i, j;

    lvlist = malloc(nstoch * sizeof *lvlist);
    if (lvlist == NULL) {
	return NULL;
    }

    for (i=0; i<nstoch; i++) {
	lvlist[i] = malloc((order + 1) * sizeof **lvlist);
	if (lvlist[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(lvlist[j]);
	    }
	    free(lvlist);
	    lvlist = NULL;
	}
	lvlist[i][0] = order;
    }

    return lvlist;
}

static int var_lists_init (struct var_lists *vl,
			   int ndet, int nstoch, 
			   int order)
{
    int nreg = 1 + ndet + nstoch * order;
    int ntest = nstoch * 2 + ndet;

    vl->detvars = NULL;
    vl->stochvars = NULL;
    vl->reglist = NULL;
    vl->testlist = NULL;
    vl->lagvlist = NULL;

    vl->detvars = malloc((ndet + 1) * sizeof *vl->detvars);
    vl->stochvars = malloc((nstoch + 1) * sizeof *vl->stochvars);
    vl->reglist = malloc((nreg + 1) * sizeof *vl->reglist);
    vl->testlist = malloc((ntest + 1) * sizeof *vl->testlist);

    if (vl->detvars == NULL || vl->stochvars == NULL ||
	vl->reglist == NULL || vl->testlist == NULL) {
	goto bailout;
    }

    vl->detvars[0] = ndet;
    vl->stochvars[0] = nstoch;
    vl->reglist[0] = nreg;
    vl->testlist[0] = ntest;

    vl->lagvlist = lagvlist_construct(nstoch, order);
    if (vl->lagvlist == NULL) {
	goto bailout;
    }

    return 0;
    
 bailout:

    var_lists_free(vl);

    return E_ALLOC;
}

int var_max_order (const int *list, const DATAINFO *pdinfo)
{
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int nstoch = 0, ndet = 0;
    int gotsep = 0;
    int order = 1;
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (!gotsep) {
	    nstoch++;
	} else {
	    ndet++;
	}
    }

    order = (T - ndet) / nstoch;

    while (order > 0) {
	int t1 = (order > pdinfo->t1)? order : pdinfo->t1;

	T = pdinfo->t2 - t1 + 1;
	if (nstoch * order + ndet > T) {
	    order--;
	} else {
	    break;
	}
    }

    return order - 1;
}

/* Given an incoming regression list, separate it into deterministic
   components (constant, trend, dummy variables) and stochastic
   components, and construct a list for each sort of variable.  Also
   allocate a VAR regression list that is long enough to hold the
   deterministic vars plus order lags of each stochastic var.
*/

static int organize_var_lists (const int *list, const double **Z,
			       const DATAINFO *pdinfo, int order,
			       struct var_lists *vlists)
{
    int ndet = 0, nstoch = 0;
    int gotsep = 0;
    char *d;
    int i, j, k, li;
    
    d = calloc(list[0] + 1, 1);
    if (d == NULL) {
	return E_ALLOC;
    }

    /* figure out the lengths of the lists */
    for (i=1; i<=list[0]; i++) {
	li = list[i];
	if (li == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (gotsep || 
	    !strcmp(pdinfo->varname[li], "const") ||	   
	    !strcmp(pdinfo->varname[li], "time") ||
	    gretl_isdummy(pdinfo->t1, pdinfo->t2, Z[li])) {
	    d[i] = 1;
	    ndet++;
	} else {
	    nstoch++;
	}
    }

    /* check for degrees of freedom */
    if (nstoch * order + ndet > pdinfo->t2 - pdinfo->t1 + 1) {
	free(d);
	return E_DF;
    }

    /* allocate the lists */
    if (var_lists_init(vlists, ndet, nstoch, order)) {
	free(d);
	return E_ALLOC;
    }

    /* fill out the detvars and stochvars lists */
    j = k = 1;
    for (i=1; i<=list[0]; i++) {
	if (list[i] != LISTSEP) {
	    if (d[i]) {
		vlists->detvars[j++] = list[i];
	    } else {
		vlists->stochvars[k++] = list[i];
	    }
	}
    }

    free(d);

#if VAR_DEBUG
    printlist(vlists->detvars, "deterministic vars");
    printlist(vlists->stochvars, "stochastic vars");
#endif

    return 0;
}

/* compose a VAR regression list: it may be complete, or one variable
   may be omitted (to run an F-test), or the order may be one less
   than the full VAR order (again, for an F-test)
*/

static int
compose_varlist (struct var_lists *vl, int depvar, int order, int omit, 
		 const DATAINFO *pdinfo)
{
    int l0 = 1 + vl->detvars[0] + order * vl->stochvars[0];
    int i, j, pos;
    int err = 0;

    if (omit) {
	l0 -= order;
    } 

    vl->reglist[0] = l0;
    vl->reglist[1] = depvar;

    pos = 2;
    for (i=1; i<=vl->stochvars[0]; i++) {
	if (i != omit) {
	    /* insert order lags of the given var */
	    for (j=1; j<=order; j++) {
		vl->reglist[pos++] = vl->lagvlist[i-1][j];
	    }
	}
    }

    /* append the deterministic vars */
    for (i=1; i<=vl->detvars[0]; i++) {
	vl->reglist[pos++] = vl->detvars[i];
    }

    /* now build the test list (to screen missing values) */
    pos = 1;
    for (i=1; i<=vl->stochvars[0]; i++) {
	vl->testlist[pos++] = vl->stochvars[i];
	vl->testlist[pos++] = vl->lagvlist[i-1][order];
    }
    for (i=1; i<=vl->detvars[0]; i++) {
	vl->testlist[pos++] = vl->detvars[i];
    }    

#if VAR_DEBUG
    printlist(vl->reglist, "composed VAR list");
    printlist(vl->testlist, "composed test list");
#endif

    return err;
}

static int add_model_data_to_var (GRETL_VAR *var, const MODEL *pmod, int k)
{
    int i, j;
    int v = 0, lag = 0;
    int start = pmod->ifc;
    int rowmax = var->neqns * var->order + start;
    int err = 0;

    if (k == 0) {
	/* first equation: set up storage for residuals */
	var->n = pmod->t2 - pmod->t1 + 1;
	var->E = gretl_matrix_alloc(var->n, var->neqns);
	if (var->E == NULL) {
	    err = 1;
	} else {
	    var->ifc = pmod->ifc;
	}
    }

    /* save residuals */
    if (!err) {
	for (i=0; i<var->n; i++) {
	    gretl_matrix_set(var->E, i, k, pmod->uhat[pmod->t1 + i]);
	}
    }	

    /* save coefficients */
    if (!err) {
	for (i=start; i<rowmax; i++) {
	    if ((i - start) % var->order == 0) {
		v++;
		lag = 1;
	    } else {
		lag++;
	    }
	    j = (lag - 1) * var->neqns + v - 1;
	    gretl_matrix_set(var->A, k, j, pmod->coeff[i]);
	}
    }

    return err;
}

static int var_compute_F_tests (MODEL *varmod, GRETL_VAR *var,
				struct var_lists *vl,
				double ***pZ, DATAINFO *pdinfo,
				int i, int *k) 
{
    MODEL testmod;
    double F = NADBL;
    int robust = gretl_model_get_int(varmod, "robust");
    int depvar = vl->stochvars[i + 1];
    int *outlist = NULL;
    int j, err = 0;

    if (robust) {
	outlist = malloc(varmod->list[0] * sizeof *outlist);
	if (outlist == NULL) {
	    return E_ALLOC;
	}
    }

    /* restrictions for all lags of specific variables */
    for (j=0; j<var->neqns && !err; j++) {

	compose_varlist(vl, depvar, var->order, j + 1, pdinfo);	

	if (robust) {
	    gretl_list_diff(outlist, varmod->list, vl->reglist);
	    F = robust_omit_F(outlist, varmod);
	    if (na(F)) {
		err = 1;
	    }
	} else {
	    testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    err = testmod.errcode;
	    if (!err) {
		F = ((testmod.ess - varmod->ess) / var->order) / 
		    (varmod->ess / varmod->dfd);
	    }
	    clear_model(&testmod);
	}

	if (!err) {
	    var->Fvals[*k] = F;
	    *k += 1;
	}
    }
    
    /* restrictions for last lag, all variables */
    if (!err) {
	compose_varlist(vl, depvar, var->order - 1, 0, pdinfo);	

	if (robust) {
	    gretl_list_diff(outlist, varmod->list, vl->reglist);
	    F = robust_omit_F(outlist, varmod);
	    if (na(F)) {
		err = 1;
	    }
	} else {
	    testmod = lsq(vl->reglist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    err = testmod.errcode;
	    if (!err) {
		F = ((testmod.ess - varmod->ess) / var->neqns) / 
		    (varmod->ess / varmod->dfd);
	    }
	    clear_model(&testmod);
	}

	if (!err) {
	    var->Fvals[*k] = F;
	    *k += 1;
	}
    }

    if (outlist != NULL) {
	free(outlist);
    }

    return err;
}

/* construct the respective VAR lists by adding the appropriate
   number of lags ("order") to the variables in list 

   Say the list is "x_1 const time x_2 x_3", and the order is 2.
   Then the first list should be

   x_1 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

   the second:

   x_2 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

   and so on.

   Run the regressions and print the results.
*/

static GRETL_VAR *real_var (int order, const int *inlist, 
			    double ***pZ, DATAINFO *pdinfo,
			    gretlopt opt, int *err)
{
    GRETL_VAR *var = NULL;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    struct var_lists vlists;
    gretlopt lsqopt = OPT_A;
    int i, k, neqns;

    if (order < 1) {
	fprintf(stderr, I_("Not much point in a zero-order \"VAR\" surely?\n"));
	*err = 1;
	return NULL;
    }

    if (opt & OPT_R) {
	lsqopt |= OPT_R;
    }

    *err = organize_var_lists(inlist, (const double **) *pZ, pdinfo, 
			      order, &vlists);
    if (*err) {
	return NULL;
    }

    /* generate the required lags */
    if (real_list_laggenr(vlists.stochvars, pZ, pdinfo, 
			  order, vlists.lagvlist)) {
	*err = E_ALLOC;
	goto var_bailout;
    }

    neqns = vlists.stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for t1 and t2 while we're at it */
    *err = compose_varlist(&vlists, vlists.stochvars[1], 
			   order, 0, pdinfo);
    if (*err) {
	*err = E_DATA;
	goto var_bailout;
    }

    /* sort out sample range */
    if (check_for_missing_obs(vlists.testlist, &pdinfo->t1, &pdinfo->t2,
			      (const double **) *pZ, NULL)) {
	*err = E_MISSDATA;
	goto var_bailout;
    }
    
    var = gretl_VAR_new(neqns, order, pdinfo);
    if (var == NULL) {
	*err = E_ALLOC;
	goto var_bailout;
    }

    k = 0;

    for (i=0; i<neqns && !*err; i++) {
	MODEL *pmod = var->models[i];

	compose_varlist(&vlists, vlists.stochvars[i + 1], 
			order, 0, pdinfo);

	*pmod = lsq(vlists.reglist, pZ, pdinfo, VAR, lsqopt, 0.0);

	if (pmod->errcode) {
	    *err = pmod->errcode;
	} else {
	    pmod->aux = AUX_VAR;
	    pmod->ID = i + 1;
	}

	if (!*err) {
	    *err = add_model_data_to_var(var, pmod, i);
	}

	if (!*err) {
	    *err = var_compute_F_tests(pmod, var, &vlists, pZ, pdinfo, i, &k);
	}
    }

 var_bailout:

    var_lists_free(&vlists);

    /* reset sample range */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    if (!*err) {
	*err = gretl_VAR_do_error_decomp(var->n, var->neqns, var->E, var->C);
    }

    if (*err) {
	gretl_VAR_free(var);
	var = NULL;
    }

    return var;
}

/**
 * simple_VAR:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opts: if OPT_R, use robust VCV,
 *        if OPT_V, print impulse responses.
 * @prn: gretl printing struct.
 *
 * Estimate a vector auto-regression (VAR) and print the results.
 *
 * Returns: 0 on successful completion, 1 on error.
 */

int simple_VAR (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		gretlopt opt, PRN *prn)
{
    GRETL_VAR *var;
    int err = 0;

    var = real_var(order, list, pZ, pdinfo, opt, &err);

    if (var != NULL) {
	gretl_VAR_print(var, pdinfo, opt, prn);
	gretl_VAR_free(var);
    }

    return err;
}

/* "full" version returns pointer to VAR struct -- invoked by gui */

GRETL_VAR *full_VAR (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opt, PRN *prn)
{
    GRETL_VAR *var;
    int err = 0;

    var = real_var(order, list, pZ, pdinfo, opt, &err);

    if (var != NULL) {
	gretl_VAR_print(var, pdinfo, opt, prn);
    }

    return var;
}

static double df_pvalue_from_plugin (double tau, int n, int niv, int itv)
{
    char datapath[FILENAME_MAX];
    void *handle;
    double (*mackinnon_pvalue)(double, int, int, int, char *);
    double pval = NADBL;
    static int nodata;
    
    if (nodata) {
	return pval;
    }

    mackinnon_pvalue = get_plugin_function("mackinnon_pvalue", &handle);
    if (mackinnon_pvalue == NULL) {
	nodata = 1;
        return pval;
    }

    strcpy(datapath, gretl_lib_path());
#ifdef WIN32
    append_dir(datapath, "plugins");
#endif

    pval = (*mackinnon_pvalue)(tau, n, niv, itv, datapath);

#if 0
    fprintf(stderr, "getting pval: tau=%g, n=%d, niv=%d, itv=%d: pval=%g\n",
	    tau, n, niv, itv, pval);
#endif

    close_plugin(handle);

    if (*datapath == '\0') {
	nodata = 1;
    } 

    return pval;
}

/**
 * coint:
 * @order: lag order for the test.
 * @list: specifies the variables to use.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: if OPT_N, do not an include a constant in the
 *       cointegrating regression.
 * @prn: gretl printing struct.
 *
 * Test for cointegration.  
 *
 * Returns: 0 on successful completion.
 *
 */

int coint (int order, const int *list, double ***pZ, 
	   DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    int i, t, n, nv, l0 = list[0];
    int hasconst = gretl_list_has_const(list);
    MODEL cmod;
    int *cointlist = NULL;

    if (order <= 0 || list[0] - hasconst < 2) {
	strcpy(gretl_errmsg, "coint: needs a positive lag order "
	       "and at least two variables");
	return 1;
    }

    gretl_model_init(&cmod);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (list[i] == 0) {
	    continue;
	}
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	real_adf_test(list[i], order, 1, pZ, pdinfo, OPT_NONE, 
		      ADF_EG_TEST, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (!hasconst && !(opt & OPT_N)) {
	/* add const to coint regression list */
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<=l0; i++) {
	    cointlist[i] = list[i];
	}
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else {
	cointlist = gretl_list_copy(list);
	if (cointlist == NULL) {
	    return E_ALLOC;
	}
    }

    pprintf(prn, _("Step %d: cointegrating regression\n"), l0 + 1);
    
    cmod = lsq(cointlist, pZ, pdinfo, OLS, OPT_NONE, 0.0); 
    cmod.aux = AUX_COINT;
    printmodel(&cmod, pdinfo, OPT_NONE, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (dataset_add_series(1, pZ, pdinfo)) {
	return E_ALLOC;
    }
    nv = pdinfo->v - 1;

    for (t=0; t<cmod.t1; t++) {
	(*pZ)[nv][t] = NADBL;
    }
    for (t=cmod.t1; t<=cmod.t2; t++) {
	(*pZ)[nv][t] = cmod.uhat[t];
    }
    for (t=cmod.t2+1; t<n; t++) {
	(*pZ)[nv][t] = NADBL;
    }

    strcpy(pdinfo->varname[nv], "uhat");

    pputc(prn, '\n');
    pprintf(prn, _("Step %d: Dickey-Fuller test on residuals\n"), l0 + 2);

    /* Run (A)DF test on the residuals */
    real_adf_test(pdinfo->v - 1, order, 1 + cmod.ncoeff - cmod.ifc, 
		  pZ, pdinfo, OPT_N, ADF_EG_TEST | ADF_PRINT_ACK, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));

    /* clean up and get out */
    clear_model(&cmod);
    free(cointlist);
    dataset_drop_last_variables(1, pZ, pdinfo);

    return 0;
}

static int *adf_prepare_vars (int order, int varno,
			      double ***pZ, DATAINFO *pdinfo)
{
    int i, orig_t1 = pdinfo->t1;
    int *list;
    int err = 0;

    if (varno == 0) {
	return NULL;
    }

    list = malloc((6 + order) * sizeof *list);
    if (list == NULL) {
	return NULL;
    }

    /* temporararily reset sample */
    pdinfo->t1 = 0;

    /* generate first difference of the given variable */
    list[1] = diffgenr(varno, pZ, pdinfo, 0);
    if (list[1] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }	

    /* generate lag of given var */
    list[2] = laggenr(varno, 1, pZ, pdinfo); 
    if (list[2] < 0) {
	pdinfo->t1 = orig_t1;
	free(list);
	return NULL;
    }

    /* undo reset sample */
    pdinfo->t1 = orig_t1;

    /* generate lags of difference for augmented test */
    for (i=1; i<=order && !err; i++) {
	int lnum = laggenr(list[1], i, pZ, pdinfo);

	if (lnum < 0) {
	    fprintf(stderr, "Error generating lag variable\n");
	    err = 1;
	} else {
	    list[2 + i] = lnum;
	} 
    } 

    return list;
}

#define ADF_DEBUG 0

static int auto_adjust_order (int *list, int order_max,
			      double ***pZ, DATAINFO *pdinfo,
			      PRN *prn)
{
    MODEL kmod;
    double tstat, pval = 1.0;
    int i, k = order_max;

    for (k=order_max; k>0; k--) {
	int j = k;

	if (list[list[0]] == 0) j++;

	kmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);

	if (kmod.errcode) {
	    clear_model(&kmod);
	    fprintf(stderr, "adf: model failed in auto_adjust_order()\n");
	    k = -1;
	    break;
	}

#if ADF_DEBUG
	printmodel(&kmod, pdinfo, OPT_NONE, prn);
#endif

	tstat = kmod.coeff[j] / kmod.sderr[j];
	clear_model(&kmod);
	pval = normal_pvalue_2(tstat);

	if (pval > 0.10) {
#if ADF_DEBUG
	    fprintf(stderr, "auto_adjust_order: lagged difference not "
		    "significant at order %d (t = %g)\n", k, tstat);
#endif
	    if (k == 1) {
		k = 0;
		break;
	    } else {
		for (i=k+2; i<list[0]; i++) {
		    list[i] = list[i+1];
		}
		list[0] -= 1;
	    }
	} else {
#if ADF_DEBUG
	    fprintf(stderr, "auto_adjust_order: lagged difference is "
		    "significant at order %d (t = %g)\n", k, tstat);
#endif
	    break;
	}
    }

    return k;
}

static void copy_list_values (int *targ, const int *src)
{
    int i;

    for (i=0; i<=src[0]; i++) {
	targ[i] = src[i];
    }
}

static void 
print_adf_results (int order, double DFt, double pv, const MODEL *dfmod,
		   int dfnum, const char *vname, int *blurb_done,
		   unsigned char flags, int i, PRN *prn)
{
    const char *models[] = {
	"(1 - L)y = (a-1)*y(-1) + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + e"
    };
    const char *aug_models[] = {
	"(1 - L)y = (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + (a-1)*y(-1) + ... + e",
	"(1 - L)y = b0 + b1*t + b2*t^2 + (a-1)*y(-1) + ... + e"
    };
    const char *teststrs[] = {
	N_("test without constant"),
	N_("test with constant"),
	N_("with constant and trend"),
	N_("with constant and quadratic trend")
    };

    char pvstr[48];

    if (prn == NULL) return;

    if (na(pv)) {
	sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
    } else {
	sprintf(pvstr, "%s %.4g", 
		(order > 0)? _("asymptotic p-value") : _("p-value"), 
		pv);
    } 

    if (*blurb_done == 0) {
	if (order > 0) {
	    pprintf(prn, _("\nAugmented Dickey-Fuller tests, order %d, for %s\n"),
		    order, vname);
	} else {
	    pprintf(prn, _("\nDickey-Fuller tests for %s\n"), vname);
	}
	pprintf(prn, _("sample size %d\n"), dfmod->nobs);
	pputs(prn, _("unit-root null hypothesis: a = 1"));
	pputs(prn, "\n\n");
	*blurb_done = 1;
    }

    pprintf(prn, "   %s\n", _(teststrs[i]));

    if (!(flags & ADF_EG_TEST)) {
	pprintf(prn, "   %s: %s\n", _("model"), 
		(order > 0)? aug_models[i] : models[i]);
    }

    pprintf(prn, "   %s: %g\n"
	    "   %s: t = %g\n"
	    "   %s\n",
	    _("estimated value of (a - 1)"), dfmod->coeff[dfnum],
	    _("test statistic"), DFt,
	    pvstr);	
}

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, unsigned char flags,
			  PRN *prn)
{
    MODEL dfmod;

    int orig_nvars = pdinfo->v;
    int blurb_done = 0;
    int auto_order = 0;
    int order_max = 0;
    int *list;
    int *biglist = NULL;
    double DFt = NADBL;
    double pv = NADBL;
    char mask[4] = {0};
    int i, itv;
    int err = 0;

#if ADF_DEBUG
    fprintf(stderr, "real_adf_test: got order = %d\n", order);
#endif

    if (order < 0) {
	auto_order = 1;
	order = -order;
    }

    order_max = order;

    list = adf_prepare_vars(order, varno, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    if (auto_order) {
	int tmp = list[0];

	list[0] = order + 5;
	biglist = gretl_list_copy(list);
	if (biglist == NULL) {
	    free(list);
	    return E_ALLOC;
	}
	list[0] = tmp;
    }

    gretl_model_init(&dfmod);

    if (opt == OPT_NONE || opt == OPT_V) {
	/* default display */
	mask[1] = mask[2] = mask[3] = 1;
    } else {
	if (opt & OPT_N) {
	    /* nc model */
	    mask[0] = 1;
	}
	if (opt & OPT_C) {
	    /* c */
	    mask[1] = 1;
	}
	if (opt & OPT_T) {
	    /* ct */
	    mask[2] = 1;
	}
	if (opt & OPT_R) {
	    /* ctt */
	    mask[3] = 1;
	}
    }

    for (i=0; i<4; i++) {
	int dfnum = (i > 0);

	if (mask[i] == 0) {
	    continue;
	}

	if (auto_order) {
	    order = order_max;
	    copy_list_values(list, biglist);
	}

	list[0] = 2 + order + i;

	if (i > 0) {
	    list[list[0]] = 0;
	} 

	if (i >= 2) {
	    list[3 + order] = gettrend(pZ, pdinfo, 0);
	    if (list[3 + order] == TREND_FAILED) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (i > 2) {
	    list[4 + order] = gettrend(pZ, pdinfo, 1);
	    if (list[4 + order] == TREND_FAILED) {
		err = E_ALLOC;
		goto bailout;
	    }
	}

	if (auto_order) {
	    order = auto_adjust_order(list, order_max, pZ, pdinfo, prn);
	    if (order < 0) {
		err = 1;
		clear_model(&dfmod);
		goto bailout;
	    }
	}

	dfmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (dfmod.errcode) {
	    fprintf(stderr, "adf_test: dfmod.errcode = %d\n", 
		    dfmod.errcode);
	    err = dfmod.errcode;
	    clear_model(&dfmod);
	    goto bailout;
	}

	DFt = dfmod.coeff[dfnum] / dfmod.sderr[dfnum];

	itv = (i == 0)? UR_NO_CONST :
	    (i == 1)? UR_CONST : 
	    (i == 2)? UR_TREND :
	    UR_TREND_SQUARED;

	pv = df_pvalue_from_plugin(DFt, 
				   /* use asymptotic p-value for augmented case */
				   (order > 0)? 0 : dfmod.nobs, 
				   niv, itv);

	if (!(opt & OPT_Q)) {
	    print_adf_results(order, DFt, pv, &dfmod, dfnum, pdinfo->varname[varno],
			      &blurb_done, flags, i, prn);
	}

	if (opt & OPT_V) {
	    /* verbose */
	    dfmod.aux = (order > 0)? AUX_ADF : AUX_DF;
	    if (!na(pv)) {
		gretl_model_set_int(&dfmod, "dfnum", dfnum + 2);
		gretl_model_set_double(&dfmod, "dfpval", pv);
	    }
	    printmodel(&dfmod, pdinfo, OPT_NONE, prn);
	} else if (!(opt & OPT_Q)) {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (!err) {
	if (!(flags & ADF_EG_TEST)) {
	    record_test_result(DFt, pv, "Dickey-Fuller");
	}
	if ((flags & ADF_PRINT_ACK) && !(opt & OPT_Q)) {
	    pputs(prn, _("P-values based on MacKinnon (JAE, 1996)\n"));
	}	
    }

 bailout:

    free(list);

    if (biglist != NULL) {
	free(biglist);
    }

    dataset_drop_last_variables(pdinfo->v - orig_nvars, pZ, pdinfo);

    return err;
}

/**
 * adf_test:
 * @order: lag order for the test.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller test for 
 * a unit root.
 *
 * Returns: 0 on successful completion, non-zero on error.
 *
 */

int adf_test (int order, int varno, double ***pZ,
	      DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    return real_adf_test(varno, order, 1, pZ, pdinfo, opt, 
			 ADF_PRINT_ACK, prn);
}

/**
 * kpss_test:
 * @order: window size for Bartlett smoothing.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @opt: option flag.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the KPSS test for 
 * stationarity.
 *
 * Returns: 0 on successful completion, non-zero on error.
 *
 */

int kpss_test (int order, int varno, double ***pZ,
	       DATAINFO *pdinfo, gretlopt opt, PRN *prn)
{
    MODEL KPSSmod;
    int list[4];
    int hastrend = 0;
    double s2 = 0.0;
    double cumsum = 0.0, cumsum2 = 0.0;
    double teststat;
    double *autocov;
    double et;

    int i, t;
    int t1, t2, T;

    /* sanity check */
    if (order < 0 || varno <= 0 || varno >= pdinfo->v) {
	return 1;
    }

    if (opt & OPT_T) {
	hastrend = 1;
    }

    list[0] = (2 + hastrend);
    list[1] = varno;
    list[2] = 0;
    if (hastrend) {
	list[3] = gettrend(pZ, pdinfo, 0);
    }

    /* OPT_M: reject missing values within sample range */
    KPSSmod = lsq(list, pZ, pdinfo, OLS, OPT_A | OPT_M, 0.0);
    if (KPSSmod.errcode) {
	clear_model(&KPSSmod);
	return KPSSmod.errcode;
    }

    t1 = KPSSmod.t1;
    t2 = KPSSmod.t2;
    T = KPSSmod.nobs;

    if (opt & OPT_V) {
	KPSSmod.aux = AUX_KPSS;
	printmodel(&KPSSmod, pdinfo, OPT_NONE, prn);
    }
  
    autocov = malloc(order * sizeof *autocov);
    if (autocov == NULL) {
	return E_ALLOC;
    }
  
    for (i=0; i<order; i++) {
	autocov[i] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	et = KPSSmod.uhat[t];
	if (na(et)) {
	    continue;
	}
	cumsum += et;
	cumsum2 += cumsum * cumsum;
	s2 += et * et;
	for (i=0; i<order; i++) {
	    int s = i + 1;

	    if (t - s >= t1) {
		autocov[i] += et * KPSSmod.uhat[t - s];
	    }
	}
#ifdef KPSS_DEBUG
	fprintf(stderr, "%d: %#12.4g %#12.4g %#12.4g %#12.4g \n", 
		t, et, KPSSmod.uhat[t-1], s2, cumsum2);
#endif
    }

    for (i=0; i<order; i++) {
	double wt = 1.0 - ((double) (i + 1)) / (order + 1);

	s2 += 2.0 * wt * autocov[i];
    }

    s2 /= T;
    teststat = cumsum2 / (s2 * T * T);

    if (opt & OPT_V) {
	pprintf(prn, "  %s: %g\n", _("Robust estimate of variance"), s2);
	pprintf(prn, "  %s: %g\n", _("Sum of squares of cumulated residuals"), 
		cumsum2);
    }

    if (!(opt & OPT_Q)) {
	pprintf(prn, _("\nKPSS test for %s %s\n\n"), pdinfo->varname[varno],
		(hastrend)? _("(including trend)") : _("(without trend)"));

	pprintf(prn, _("Lag truncation parameter = %d\n"), order);

	pprintf(prn, "%s = %g\n\n", _("Test statistic"), teststat);

	pprintf(prn, "		    10%%\t   5%%\t 2.5%%\t   1%%\n");

	if (hastrend) {
	    pprintf(prn, "%s: 0.119\t0.146\t0.176\t0.216\n\n", _("Critical values"));
	} else {
	    pprintf(prn, "%s: 0.347\t0.463\t0.574\t0.739\n\n", _("Critical values"));
	}
    }

    record_test_result(teststat, NADBL, "KPSS");

    clear_model(&KPSSmod);

    free(autocov);

    return 0;
}

static int 
has_time_trend (const int *varlist, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    int origv = pdinfo->v;
    int tlist[4];
    int trends = 0;
    MODEL tmod;

    gretl_model_init(&tmod);

    tlist[0] = 3;
    tlist[2] = 0;

    for (i=1; i<=varlist[0]; i++) {
	int vi = varlist[i];
	double tstat;

	if (vi == 0) {
	    continue;
	}

	tlist[1] = vi;

	tlist[3] = laggenr(vi, 1, pZ, pdinfo);
	if (tlist[3] < 0) {
	    trends = -1;
	    break;
	}

	tmod = lsq(tlist, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (tmod.errcode) {
	    trends = -1;
	    clear_model(&tmod);
	    break;
	}

	tstat = tmod.coeff[0] / tmod.sderr[0];
	if (t_pvalue_2(tstat, tmod.dfd) < 0.05) {
	    trends = 1;
	}

	clear_model(&tmod);

	if (trends) break;
    }

    dataset_drop_last_variables(pdinfo->v - origv, pZ, pdinfo);

    return trends;
}

static int allocate_sigmas (double ***X, double ***Y, double ***Z, int k)
{
    int i, j;
    double **Suu, **Svv, **Suv;

    Suu = malloc(k * sizeof *Suu);
    Svv = malloc(k * sizeof *Svv);
    Suv = malloc(k * sizeof *Suv);

    if (Suu == NULL || Svv == NULL || Suv == NULL) return 1;

    for (i=0; i<k; i++) {
	Suu[i] = malloc(k * sizeof **Suu);
	Svv[i] = malloc(k * sizeof **Svv);
	Suv[i] = malloc(k * sizeof **Suv);
	if (Suu[i] == NULL || Svv[i] == NULL || Suv[i] == NULL) {
	    free(Suu);
	    free(Svv);
	    free(Suv);
	    return 1;
	}
	for (j=0; j<k; j++) {
	    Suu[i][j] = 0.0;
	    Svv[i][j] = 0.0;
	    Suv[i][j] = 0.0;
	}
    }

    *X = Suu;
    *Y = Svv;
    *Z = Suv;

    return 0;
}

static void free_sigmas (double **X, double **Y, double **Z, int k)
{
    int i;

    for (i=0; i<k; i++) {
	if (X != NULL) free(X[i]);
	if (Y != NULL) free(Y[i]);
	if (Z != NULL) free(Z[i]);
    }

    free(X);
    free(Y);
    free(Z);
}

static void scatter_product (const double **u, const double **v, 
			     double **X, int T, int k)
{
    int i, j, t;

    for (t=0; t<T; t++) {
	for (i=0; i<k; i++) {
	    for (j=0; j<k; j++) {
		X[i][j] += u[i][t] * v[j][t];
	    }
	}
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    X[i][j] /= (double) T;
	}
    }
}

static void 
print_sigmas (const double **X, const double **Y, const double **Z, 
	      int k, PRN *prn)
{
    int i, j, l;
    const double **P = NULL;

    pprintf(prn, "\n%s\n\n", _("Sample variance-covariance matrices for residuals"));

    for (l=0; l<3; l++) {
	if (l == 0) {
	    P = X;
	    pprintf(prn, " %s\n\n", _("VAR system in first differences"));
	}
	else if (l == 1) {
	    P = Y;
	    pprintf(prn, " %s\n\n", _("System with levels as dependent variable"));
	}
	else {
	    P = Z;
	    pprintf(prn, " %s\n\n", _("Cross-products"));
	}
	for (i=0; i<k; i++) {
	    for (j=0; j<k; j++) {
		pprintf(prn, "%#12.6g", P[i][j]);
	    }
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');
    }
}

static int 
johansen_complete (const double **X, const double **Y, const double **Z,
		   int k, int T, int trends, PRN *prn)
{
    void *handle = NULL;
    int (*johansen) (const double **, const double **, const double **,
		     int, int, int, PRN *);
    int err = 0;

    *gretl_errmsg = 0;
    
    johansen = get_plugin_function("johansen_eigenvals", &handle);

    if (johansen == NULL) {
	err = 1;
    } else {
	err = (* johansen) (X, Y, Z, k, T, trends, prn);
	close_plugin(handle);
    }
    
    return err;
}

#undef JOHANSEN_DEBUG

static int johansen_VAR (int order, const int *inlist, 
			 double ***pZ, DATAINFO *pdinfo,
			 struct var_resids *resids, 
			 gretlopt opt, PRN *prn)
{
    int i, k, neqns;
    int oldt1 = pdinfo->t1;
    int oldt2 = pdinfo->t2;
    struct var_lists vlists;
    MODEL var_model;
    MODEL jmod;
    int err = 0;

    err = organize_var_lists(inlist, (const double **) *pZ, pdinfo, 
			     order, &vlists);
    if (err) {
	return err;
    }

    /* generate the required lags */
    if (real_list_laggenr(vlists.stochvars, pZ, pdinfo, 
			  order, vlists.lagvlist)) {
	err = E_ALLOC;
	goto var_bailout;
    }

    neqns = vlists.stochvars[0];    

    /* compose base VAR list (entry 1 will vary across equations);
       assemble test list for t1 and t2 while we're at it */
    err = compose_varlist(&vlists, vlists.stochvars[1], 
			  order, 0, pdinfo);
    if (err) {
	err = E_DATA;
	goto var_bailout;
    }

    if (check_for_missing_obs(vlists.testlist, &pdinfo->t1, &pdinfo->t2,
			      (const double **) *pZ, NULL)) {
	err = E_MISSDATA;
	goto var_bailout;
    }
    
    gretl_model_init(&var_model);

    if (opt & OPT_V) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), order);
    }

    /* apparatus for saving the residuals */
    resids->m = 2 * neqns;
    resids->uhat = malloc(2 * neqns * sizeof *resids->uhat);
    if (resids->uhat == NULL) {
	err = E_ALLOC;
	goto var_bailout;
    }

    k = 0;
    for (i=0; i<neqns; i++) {

	compose_varlist(&vlists, vlists.stochvars[i + 1], 
			order, 0, pdinfo);

	var_model = lsq(vlists.reglist, pZ, pdinfo, VAR, OPT_A, 0.0);

	if ((err = var_model.errcode)) {
	    goto var_bailout;
	}

	/* save the residuals */
	resids->t1 = var_model.t1;
	resids->t2 = var_model.t2;
	resids->uhat[i] = var_model.uhat;
	var_model.uhat = NULL;

	if (opt & OPT_V) {
	    var_model.aux = AUX_VAR;
	    var_model.ID = i + 1;
	    printmodel(&var_model, pdinfo, OPT_NONE, prn);
	}

	/* estimate equations for Johansen test */
	vlists.reglist[1] = resids->levels_list[i + 1]; 
	jmod = lsq(vlists.reglist, pZ, pdinfo, VAR, OPT_A, 0.0);

	if (opt & OPT_V) {
	    jmod.aux = AUX_JOHANSEN;
	    jmod.ID = -1;
	    printmodel(&jmod, pdinfo, OPT_NONE, prn);
	}

	resids->uhat[i + neqns] = jmod.uhat;
	jmod.uhat = NULL;
	clear_model(&jmod);

	clear_model(&var_model);
    }

    pputc(prn, '\n');

 var_bailout:

    var_lists_free(&vlists);

    /* reset sample range */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    return err;
}

int johansen_test (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn)
{
    PRN *varprn = NULL;
    struct var_resids resids;
    int err = 0;
    int i, j;
    int orig_t1 = pdinfo->t1;
    int orig_v = pdinfo->v;
    int *varlist;
    int hasconst = 0;
    int l0 = list[0];
    int trends = 0;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    hasconst = 1;
	    break;
	}
    }

    if (order <= 0 || list[0] - hasconst < 2) {
	strcpy(gretl_errmsg, "coint2: needs a positive lag order "
	       "and at least two variables");
	return 1;
    }

    /* we're assuming that the list we are fed is in levels */
    resids.levels_list = gretl_list_new(l0);
    if (resids.levels_list == NULL) {
	return E_ALLOC;
    }

    varlist = gretl_list_new(l0 + 1);
    if (varlist == NULL) {
	free(resids.levels_list);
	return E_ALLOC;
    }

    varlist[0] = resids.levels_list[0] = l0 - hasconst;

    j = 1;
    for (i=1; i<=list[0]; i++) {
	int lnum;

	if (list[i] != 0) {
	    lnum = laggenr(list[i], 1, pZ, pdinfo);
	    if (lnum < 0) {
		free(varlist);
		free(resids.levels_list);
		return E_DATA;
	    }
	    resids.levels_list[j++] = lnum;
	}
    }

    /* now get first differences and put them into list */
    j = 1;
    for (i=1; i<=list[0]; i++) {
	if (list[i] != 0) {
	    varlist[j] = diffgenr(list[i], pZ, pdinfo, 0);
	    if (varlist[j] < 0) {
		free(varlist);
		free(resids.levels_list);
		return E_DATA;
	    } 
	    j++;
	}
    }

    /* add the constant to the VAR list */
    varlist[0] += 1;
    varlist[varlist[0]] = 0;

    if (opt & OPT_V) {
	varprn = prn;
    } else {
	varprn = gretl_print_new(GRETL_PRINT_NULL);
    }

    /* FIXME? */
    pdinfo->t1 += (order + 1);
    /* Check Hamilton: what if order for test = 1? */
    err = johansen_VAR(order - 1, varlist, pZ, pdinfo, 
		       &resids, opt, varprn); 
    
    if (varprn != NULL) {
	gretl_print_destroy(varprn);
    }

    if (!err) {
	int k = resids.m / 2;
	int T = resids.t2 - resids.t1 + 1;
	double **Suu, **Svv, **Suv;
	double **u = NULL, **v = NULL;
	char stobs[OBSLEN], endobs[OBSLEN];

	if (allocate_sigmas(&Suu, &Svv, &Suv, k)) {
	    err = E_ALLOC;
	    goto johansen_bailout;
	}

	u = malloc(k * sizeof *u);
	v = malloc(k * sizeof *v);

	if (u == NULL || v == NULL) {
	    err = E_ALLOC;
	    goto johansen_bailout;
	}

	for (i=0; i<k; i++) {
	    u[i] = &(resids.uhat[i][resids.t1]);
	    v[i] = &(resids.uhat[i + k][resids.t1]);
	}

	scatter_product((const double **) u, (const double **) u, Suu, T, k);
	scatter_product((const double **) v, (const double **) v, Svv, T, k);
	scatter_product((const double **) u, (const double **) v, Suv, T, k);

	pprintf(prn, "%s:\n", _("Johansen test"));
	pprintf(prn, "%s = %d\n", _("Number of equations"), k);
	pprintf(prn, "%s: %s - %s (T = %d)\n", _("Estimation period"),
		ntodate(stobs, resids.t1, pdinfo), 
		ntodate(endobs, resids.t2, pdinfo), T);

	if (opt & OPT_V) {
	    print_sigmas((const double **) Suu, 
			 (const double **) Svv, 
			 (const double **) Suv, k, prn);
	}

#ifdef JOHANSEN_DEBUG
	for (i=0; i<resids.m; i++) {
	    char datestr[OBSLEN];
	    int t;

	    pprintf(prn, "Residuals from VAR model %d\n", i);
	    for (t=resids.t1; t<=resids.t2; t++) {
		ntodate(datestr, t, pdinfo);
		pprintf(prn, "%8s %#.*g\n", datestr, 
			GRETL_DIGITS, resids.uhat[i][t]);
	    }
	}
#endif

	trends = has_time_trend(list, pZ, pdinfo);
	if (trends == -1) {
	    pprintf(prn, "%s\n", _("Error checking for time trends"));
	    goto johansen_bailout;
	}

	/* now get johansen plugin to finish the job */
	err = johansen_complete((const double **) Suu, 
				(const double **) Svv, 
				(const double **) Suv, k, T, trends, prn);

    johansen_bailout:
	
	for (i=0; i<resids.m; i++) {
	    free(resids.uhat[i]);
	}
	free(resids.uhat);

	free_sigmas(Suu, Svv, Suv, k);
	free(u);
	free(v);
    } 

    free(resids.levels_list);
    free(varlist);

    pdinfo->t1 = orig_t1;

    dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);

    return err;
}

/**
 * gretl_VAR_print:
 * @var: pointer to VAR struct.
 * @pdinfo: dataset information.
 * @opt: if %OPT_V, include impulse responses and forecast
 * variance decomposition.
 * @prn: pointer to printing struct.
 *
 * Prints the models in @var along with relevant F-tests.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_VAR_print (GRETL_VAR *var, const DATAINFO *pdinfo, gretlopt opt, 
		     PRN *prn)
{
    int i, j, k, v;
    int dfd = (var->models[0])->dfd;
    int tex = tex_format(prn);
    int pause = 0;

    if (prn == NULL) {
	return 0;
    }

    if (tex) {
	pputs(prn, "\\noindent");
	pprintf(prn, I_("\nVAR system, lag order %d\n\n"), var->order);
    } else {
	pause = gretl_get_text_pause();
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), var->order);
    }

    k = 0;

    for (i=0; i<var->neqns; i++) {

	printmodel(var->models[i], pdinfo, OPT_NONE, prn);

	if (pause) {
	    scroll_pause();
	}

	if (tex) {
	    pputs(prn, "\n\\begin{center}\n");
	    pprintf(prn, "%s\\\\[1em]\n", I_("F-tests of zero restrictions"));
	    pputs(prn, "\\begin{tabular}{lll}\n");
	} else {
	    pputs(prn, _("\nF-tests of zero restrictions:\n\n"));
	}

	for (j=0; j<var->neqns; j++) {
	    v = (var->models[j])->list[1];
	    if (tex) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, I_("p-value %f"), fdist(var->Fvals[k], var->order, dfd));
		pputs(prn, "\\\\\n");
	    } else {
		pprintf(prn, _("All lags of %-8s "), pdinfo->varname[v]);
		pprintf(prn, "F(%d, %d) = %10g, ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->order, dfd));
	    }
	    k++;
	}

	if (var->order > 1) {
	    if (tex) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->neqns, dfd));
	    } else {
		pprintf(prn, _("All vars, lag %-6d "), var->order);
		pprintf(prn, "F(%d, %d) = %10g, ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->neqns, dfd));
	    } 
	    k++;
	}

	if (tex) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else {
	    pputc(prn, '\n');
	    if (pause) {
		scroll_pause();
	    }
	}
    }

    pputc(prn, '\n');

    if (opt & OPT_V) {
	int horizon = default_VAR_horizon(pdinfo);

	for (i=0; i<var->neqns; i++) {
	    gretl_VAR_print_impulse_response(var, i, horizon, pdinfo, pause, prn);
	    gretl_VAR_print_fcast_decomp(var, i, horizon, pdinfo, pause, prn);
	}
    }

    return 0;
}

void gretl_VAR_assign_name (GRETL_VAR *var)
{
    static int n = 0;

    if (var->name != NULL) {
	free(var->name);
    }

    var->name = malloc(8);
    if (var->name != NULL) {
	sprintf(var->name, "%s %d", _("VAR"), ++n);
    }
}

void gretl_VAR_assign_specific_name (GRETL_VAR *var, const char *name)
{
    if (var->name != NULL) {
	free(var->name);
    }

    var->name = gretl_strdup(name);
}

const char *gretl_VAR_get_name (const GRETL_VAR *var)
{
    return var->name;
}

int gretl_VAR_add_resids_to_dataset (GRETL_VAR *var, int eqnum,
				     double ***pZ, DATAINFO *pdinfo)
{
    MODEL *pmod = var->models[eqnum];
    int i, t;

    if (dataset_add_series(1, pZ, pdinfo)) return E_ALLOC;

    i = pdinfo->v - 1;

    for (t=0; t<pdinfo->n; t++) {
	if (t < pmod->t1 || t > pmod->t2) {
	    (*pZ)[i][t] = NADBL;
	} else {
	    (*pZ)[i][t] = pmod->uhat[t];
	}
    }

    sprintf(pdinfo->varname[i], "uhat%d", eqnum + 1);
    sprintf(VARLABEL(pdinfo, i), _("residual from VAR system, equation %d"), 
	    eqnum + 1);

    return 0;
}

int gretl_VAR_get_highest_variable (const GRETL_VAR *var,
				    const DATAINFO *pdinfo)
{
    int vmax = 0;

    if (var->models != NULL && var->neqns >= 1) {
	vmax = highest_numbered_var_in_model(var->models[0], pdinfo);
    }

    return vmax;
}

#include "irfboot.c"
