/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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
#include "gretl_private.h"

#undef VAR_DEBUG 

struct GRETL_VAR_ {
    int neqns;         /* number of equations in system */
    int order;         /* lag order */
    int n;             /* number of observations */
    int ifc;           /* equations include a constant (1) or not (0) */
    gretl_matrix *A;   /* augmented coefficient matrix */
    gretl_matrix *E;   /* residuals matrix */
    gretl_matrix *C;   /* augmented Cholesky-decomposed error matrix */
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

enum {
    VAR_PRINT_MODELS =      1 << 0,
    VAR_DO_FTESTS    =      1 << 1,
    VAR_PRINT_PAUSE  =      1 << 2,
    VAR_IMPULSE_RESPONSES = 1 << 3,
    VAR_SAVE =              1 << 4
} var_flags;

#define TREND_FAILED 9999

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, int cointcode, PRN *prn);

/* ...................................................................  */

static void pad_var_coeff_matrix (GRETL_VAR *var)
{
    int i, j;
    int rowmax = var->neqns * var->order;

    for (i=var->neqns; i<rowmax; i++) {
	for (j=0; j<rowmax; j++) {
	    gretl_matrix_set(var->A, i, j, (j == i - var->neqns)? 1.0 : 0.0);
	}
    }
}

static int 
gretl_var_init (GRETL_VAR *var, int neqns, int order, const DATAINFO *pdinfo,
		char flags)
{
    int i, j, err = 0;
    int rows = neqns * order;

    var->neqns = neqns;
    var->order = order;
    var->A = NULL;
    var->E = NULL;
    var->C = NULL;
    var->models = NULL;
    var->Fvals = NULL;
    var->name = NULL;

    if (neqns > 0) {
	var->A = gretl_matrix_alloc(rows, rows);
	if (var->A == NULL) {
	    err = 1;
	} else {
	    pad_var_coeff_matrix(var);
	}
    }

    if (!err && neqns > 0) {
	var->C = gretl_matrix_alloc(rows, neqns);
	if (var->C == NULL) {
	    err = 1;
	    gretl_matrix_free(var->A);
	    var->A = NULL;
	} else {
	    gretl_matrix_zero(var->C);
	}
    }

    if (!err && neqns > 0) {
	var->models = malloc(neqns * sizeof *var->models);
	if (var->models == NULL) err = 1;
    } else {
	var->models = NULL;
    }

    if (var->models != NULL) {
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

    if (!err && (flags & VAR_SAVE)) {
	int m = neqns * neqns;
	
	if (order > 1) m += neqns;
#ifdef VAR_DEBUG
	fprintf(stderr, "var->Fvals: allocating %d terms\n", m);
#endif
	var->Fvals = malloc(m  * sizeof *var->Fvals);
    }

    return err;
}

void gretl_var_free (GRETL_VAR *var)
{
    int i;

    if (var == NULL) return;

    gretl_matrix_free(var->A);
    gretl_matrix_free(var->E);
    gretl_matrix_free(var->C);

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

void gretl_var_free_unnamed (GRETL_VAR *var)
{
    if (var == NULL) return;

    if (var->name == NULL || *var->name == '\0') {
	gretl_var_free(var);
    }
}

GRETL_VAR *gretl_var_new (int neqns, int order, const DATAINFO *pdinfo,
			  char flags)
{
    GRETL_VAR *var;

    var = malloc(sizeof *var);
    if (var == NULL) return NULL;

    if (gretl_var_init(var, neqns, order, pdinfo, flags)) {
	free(var);
	return NULL;
    } 

    return var;
}

/* ......................................................  */

static int gretl_var_do_error_decomp (GRETL_VAR *var)
{
    gretl_matrix *tmp = NULL;
    int i, j, err = 0;

    tmp = gretl_matrix_alloc(var->neqns, var->neqns);
    if (tmp == NULL) err = E_ALLOC;

    /* form e'e */
    if (!err && gretl_matrix_multiply_mod (var->E, GRETL_MOD_TRANSPOSE,
					   var->E, GRETL_MOD_NONE,
					   tmp)) {
	err = 1;
    }

    /* divide by T (or use df correction?) to get sigma-hat.
       Note: RATS 4 uses straight T.
    */
    gretl_matrix_divide_by_scalar(tmp, (double) var->n);

#ifdef VAR_DEBUG
    if (1) {
	PRN *prn = gretl_print_new(GRETL_PRINT_STDERR, NULL);

	gretl_matrix_print(tmp, "Sigma-hat from VAR system", prn);
	gretl_print_destroy(prn);
    }
#endif

    /* lower-triangularize and decompose */
    if (!err) {
	for (i=0; i<var->neqns-1; i++) {
	    for (j=i+1; j<var->neqns; j++) {
		gretl_matrix_set(tmp, i, j, 0.0);
	    }
	}
	err = gretl_matrix_cholesky_decomp(tmp);
    }

    /* write the decomposition into the C matrix */
    if (!err) {
	for (i=0; i<var->neqns; i++) {
	    for (j=0; j<var->neqns; j++) {
		double x = gretl_matrix_get(tmp, i, j);

		gretl_matrix_set(var->C, i, j, x);
	    }
	}
    }

    if (tmp != NULL) gretl_matrix_free(tmp);

    /* we're done with the full-length residual series */
    gretl_matrix_free(var->E);
    var->E = NULL;

    return err;
}

/* ......................................................  */

int gretl_var_get_variable_number (const GRETL_VAR *var, int k)
{
    return (var->models[k])->list[1];
}

int gretl_var_get_n_equations (const GRETL_VAR *var)
{
    return var->neqns;
}

/* ......................................................  */

#define VARS_IN_ROW 4

#define PLAIN_PRN(p) (p->format == GRETL_PRINT_FORMAT_PLAIN)
#define RTF_PRN(p)   (p->format == GRETL_PRINT_FORMAT_RTF)
#define TEX_PRN(p)   (p->format == GRETL_PRINT_FORMAT_TEX || \
                      p->format == GRETL_PRINT_FORMAT_TEX_DOC)

static void tex_print_double (double x, PRN *prn)
{
    char number[16];

    x = screen_zero(x);

    sprintf(number, "%#.*g", GRETL_DIGITS, x);

    if (x < 0.) pprintf(prn, "$-$%s", number + 1);
    else pputs(prn, number);
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

int 
gretl_var_print_impulse_response (GRETL_VAR *var, int shock,
				  int periods, const DATAINFO *pdinfo, 
				  int pause, PRN *prn)
{
    int i, t;
    int vsrc;
    int rows = var->neqns * var->order;
    gretl_matrix *rtmp, *ctmp;
    int block, blockmax;
    int err = 0;

    if (prn == NULL) return 0;

    if (shock >= var->neqns) {
	fprintf(stderr, "Shock variable out of bounds\n");
	return 1;
    }  

    if (periods == 0) {
	periods = periods_from_pd(pdinfo->pd);
    }

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) return E_ALLOC;

    ctmp = gretl_matrix_alloc(rows, var->neqns);
    if (ctmp == NULL) {
	gretl_matrix_free(rtmp);
	return E_ALLOC;
    }

    vsrc = (var->models[shock])->list[1];

    blockmax = var->neqns / VARS_IN_ROW;
    if (var->neqns % VARS_IN_ROW) blockmax++;

    for (block=0; block<blockmax && !err; block++) {
	int vtarg, k;
	char vname[16];
	double r;

	if (TEX_PRN(prn)) {
	    pputs(prn, "\\vspace{1em}\n\n");
	    pprintf(prn, I_("Responses to a one-standard error shock in %s"), 
		    tex_escape(vname, pdinfo->varname[vsrc]));

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", I_("continued"));
	    }
	    pputs(prn, "\\vspace{1em}\n\n"
		  "\\begin{tabular}{rcccc}\n");
	} else {
	    pprintf(prn, _("Responses to a one-standard error shock in %s"), 
		    pdinfo->varname[vsrc]);

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", _("continued"));
	    }
	}

	if (TEX_PRN(prn)) {
	    pprintf(prn, "%s & ", I_("period"));
	} else {
	    pprintf(prn, "%s ", _("period"));
	}

	for (i=0; i<VARS_IN_ROW; i++) {
	    k = VARS_IN_ROW * block + i;
	    if (k >= var->neqns) break;
	    vtarg = (var->models[k])->list[1];
	    if (TEX_PRN(prn)) {
		pprintf(prn, " %s ", tex_escape(vname, pdinfo->varname[vtarg]));
		if (i < VARS_IN_ROW - 1 && k < var->neqns - 1) pputs(prn, "& ");
		else pputs(prn, "\\\\");
	    } else {
		pprintf(prn, "  %8s  ", pdinfo->varname[vtarg]);
	    }
	}

	pputs(prn, "\n\n");

	for (t=0; t<periods && !err; t++) {
	    pprintf(prn, " %3d  ", t + 1);
	    if (TEX_PRN(prn)) pputs(prn, "& ");
	    if (t == 0) {
		/* calculate initial estimated responses */
		err = gretl_matrix_copy_values(rtmp, var->C);
	    } else {
		/* calculate further estimated responses */
		err = gretl_matrix_multiply(var->A, rtmp, ctmp);
		gretl_matrix_copy_values(rtmp, ctmp);
	    }

	    if (err) break;

	    for (i=0; i<VARS_IN_ROW; i++) {
		k = VARS_IN_ROW * block + i;
		if (k >= var->neqns) break;
		r = gretl_matrix_get(rtmp, k, shock);
		if (TEX_PRN(prn)) {
		    tex_print_double(r, prn);
		    if (i < VARS_IN_ROW - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else {
		    pprintf(prn, "%#12.5g ", r);
		}
	    }
	    if (TEX_PRN(prn)) pputs(prn, "\\\\\n");
	    else pputc(prn, '\n');
	}

	if (TEX_PRN(prn)) {
	    pputs(prn, "\\end{tabular}\n\n");
	} else {
	    pputc(prn, '\n');
	}

	if (pause && block < blockmax - 1) {
	    takenotes(0);
	}
    }

    if (rtmp != NULL) gretl_matrix_free(rtmp);
    if (ctmp != NULL) gretl_matrix_free(ctmp);

    return err;
}

double *
gretl_var_get_impulse_responses (GRETL_VAR *var, int targ, int shock,
				 int periods) 
{
    int t;
    int rows = var->neqns * var->order;
    gretl_matrix *rtmp, *ctmp;
    double *resp;
    int err = 0;

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

    resp = malloc(periods * sizeof *resp);
    if (resp == NULL) return NULL;

    rtmp = gretl_matrix_alloc(rows, var->neqns);
    if (rtmp == NULL) {
	free(resp);
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

	if (err) break;

	resp[t] = gretl_matrix_get(rtmp, targ, shock);
    }

    gretl_matrix_free(rtmp);
    gretl_matrix_free(ctmp);

    return resp;
}

static gretl_matrix *
gretl_var_get_fcast_decomp (GRETL_VAR *var, int targ, int periods) 
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

int 
gretl_var_print_fcast_decomp (GRETL_VAR *var, int targ,
			      int periods, const DATAINFO *pdinfo, 
			      int pause, PRN *prn)
{
    int i, t;
    int vtarg;
    gretl_matrix *vd = NULL;
    int block, blockmax;
    int err = 0;

    if (prn == NULL) return 0;

    if (targ >= var->neqns) {
	fprintf(stderr, "Target variable out of bounds\n");
	return 1;
    } 

    if (periods == 0) {
	periods = periods_from_pd(pdinfo->pd);
    } 

    vd = gretl_var_get_fcast_decomp(var, targ, periods);
    if (vd == NULL) return E_ALLOC;

    vtarg = (var->models[targ])->list[1];

    blockmax = (var->neqns + 1) / VDROWMAX;
    if ((var->neqns + 1) % VDROWMAX) blockmax++;

    for (block=0; block<blockmax; block++) {
	int k, vsrc;
	char vname[16];
	double r;

	/* print block header */
	if (TEX_PRN(prn)) {
	    pputs(prn, "\\vspace{1em}\n\n");
	    pprintf(prn, I_("Decomposition of variance for %s"), 
		    tex_escape(vname, pdinfo->varname[vtarg]));

	    if (block == 0) {
		pputs(prn, "\n\n");
	    } else {
		pprintf(prn, " (%s)\n\n", I_("continued"));
	    }
	    pputs(prn, "\\vspace{1em}\n\n"
		  "\\begin{tabular}{rccccc}\n");
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
	if (TEX_PRN(prn)) {
	    pprintf(prn, "%s & ", I_("period"));
	} else {
	    pprintf(prn, "%s ", _("period"));
	}

	/* print variable names row */
	for (i=0; i<VDROWMAX; i++) {
	    k = VDROWMAX * block + i - 1;
	    if (k < 0) {
		if (TEX_PRN(prn)) {
		    pprintf(prn, " %s & ", I_("std. error"));
		} else {
		    pprintf(prn, " %12s ", _("std. error"));
		}
		continue;
	    }
	    if (k >= var->neqns) break;
	    vsrc = (var->models[k])->list[1];
	    if (TEX_PRN(prn)) {
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
	    if (TEX_PRN(prn)) pputs(prn, "& ");

	    for (i=0; i<VDROWMAX; i++) {
		k = VDROWMAX * block + i - 1;
		if (k < 0) {
		    r = gretl_matrix_get(vd, t, var->neqns);
		    if (TEX_PRN(prn)) {
			pprintf(prn, "%g & ", r);
		    } else {
			pprintf(prn, " %14g ", r);
		    }
		    continue;
		}
		if (k >= var->neqns) break;
		r = gretl_matrix_get(vd, t, k);
		if (TEX_PRN(prn)) {
		    pprintf(prn, "$%.4f$", r);
		    if (i < VDROWMAX - 1 && k < var->neqns - 1) {
			pputs(prn, " & ");
		    }
		} else {
		    pprintf(prn, "%10.4f ", r);
		}
	    }
	    if (TEX_PRN(prn)) pputs(prn, "\\\\\n");
	    else pputc(prn, '\n');
	}

	if (TEX_PRN(prn)) {
	    pputs(prn, "\\end{tabular}\n\n");
	} else {
	    pputc(prn, '\n');
	}

	if (pause && block < blockmax - 1) {
	    takenotes(0);
	}
    }

    if (vd != NULL) gretl_matrix_free(vd);

    return err;
}

/* ......................................................  */

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
    
    if (dataset_add_vars(1, pZ, pdinfo)) {
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

/* ...................................................................  */

static void reset_list (int *list1, int *list2)
{
    int i;
    
    for (i=2; i<=list1[0]; i++) {
	list1[i] = list2[i];
    }
}

/* parse varlist (for a VAR) and determine how long the augmented list
   will be, once all the appropriate lag terms are inserted
*/

static int 
get_var_listlen (int *varlist, char *detlist, 
		 int order, int *depvars,
		 double **Z, const DATAINFO *pdinfo)
{
    int i, j = 1, v = 1;
    int gotsep = 0;

    depvars[0] = 0;

    for (i=1; i<=varlist[0]; i++) {
	if (varlist[i] == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (gotsep || 
	    strcmp(pdinfo->varname[varlist[i]], "time") == 0 ||
	    strcmp(pdinfo->varname[varlist[i]], "const") == 0 ||
	    isdummy(Z[varlist[i]], pdinfo->t1, pdinfo->t2)) {
	    v++;
	    detlist[j] = 1;
	} else {
	    v += order;
	    detlist[j] = 0;
	    depvars[0] += 1;
	    depvars[depvars[0]] = varlist[i];
	}
	varlist[j++] = varlist[i];
    }

    if (gotsep) {
	varlist[0] -= 1;
    }

    detlist[0] = varlist[0];

    return v;
}

/* ...................................................................  */

static int add_model_data_to_var (GRETL_VAR *var, const MODEL *pmod, int k)
{
    int i, j;
    int v = 0, lag = 0;
    int start = pmod->ifc;
    int rowmax = var->neqns * var->order + start;

    if (k == 0) {
	/* first equation: set up storage for residuals */
	var->n = pmod->t2 - pmod->t1 + 1;
	var->E = gretl_matrix_alloc(var->n, var->neqns);
	if (var->E == NULL) return 1;
	var->ifc = pmod->ifc;
    }

    /* save residuals */
    for (i=0; i<var->n; i++) {
	gretl_matrix_set(var->E, i, k, pmod->uhat[pmod->t1 + i]);
    }	

    /* save coefficients */
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

    return 0;
}

static int var_F_tests (MODEL *varmod, GRETL_VAR *var,
			int *varlist, int *depvars,
			double ***pZ, DATAINFO *pdinfo,
			int order, int neqns, int idx, int *k, 
			int pause, PRN *prn)
{
    MODEL testmod;
    double F = NADBL;
    int robust = gretl_model_get_int(varmod, "robust");
    int *shortlist = NULL, *outlist = NULL;
    int end, err = 0;
    int j, l;

    if (robust) {
	outlist = malloc(varmod->list[0] * sizeof *outlist);
	if (outlist == NULL) return E_ALLOC;
    }

    shortlist = malloc((varmod->list[0] + 1) * sizeof *shortlist);
    if (shortlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }
    
    shortlist[0] = varlist[0] - order;
    shortlist[1] = varlist[1];

    pputs(prn, _("\nF-tests of zero restrictions:\n\n"));

    for (j=0; j<neqns; j++) {
	reset_list(shortlist, varlist);
	for (l=1; l<=order; l++) {
	    idx = l + 1 + j * order;
	    if (idx > shortlist[0]) break;
	    shortlist[idx] = varlist[idx + order];
	}
	end = 0;
	for (l=shortlist[0]; l>idx; l--) {
	    shortlist[l] = varlist[varlist[0] - end];
	    end++;
	}
	pprintf(prn, _("All lags of %-8s "), 
		pdinfo->varname[depvars[j + 1]]);

	if (robust) {
	    gretl_list_diff(outlist, varmod->list, shortlist);
	    F = robust_omit_F(outlist, varmod);
	    if (na(F)) err = 1;
	} else {
	    testmod = lsq(shortlist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    err = testmod.errcode;
	    if (!err) {
		F = ((testmod.ess - varmod->ess) / order) / 
		    (varmod->ess / varmod->dfd);
	    }
	    clear_model(&testmod);
	}
	if (err) break;

	pprintf(prn, "F(%d, %d) = %f, ", order, varmod->dfd, F);
	pprintf(prn, _("p-value %f\n"), fdist(F, order, varmod->dfd));
	if (var != NULL) {
	    var->Fvals[*k] = F;
	    *k += 1;
	}
    }

    if (order > 1 && !err) {
	pprintf(prn, _("All vars, lag %-6d "), order);
	reset_list(shortlist, varlist);
	idx = 2;
	for (j=1; j<=neqns*(order); j++) {
	    if (j % order) {
		shortlist[idx] = varlist[j+1];
		idx++;
	    }
	}
	end = 0;
	for (l=shortlist[0]; l>=idx; l--) {
	    shortlist[l] = varlist[varlist[0]-end];
	    end++;
	}

	if (robust) {
	    gretl_list_diff(outlist, varmod->list, shortlist);
	    F = robust_omit_F(outlist, varmod);
	    if (na(F)) err = 1;
	} else {
	    testmod = lsq(shortlist, pZ, pdinfo, VAR, OPT_A, 0.0);
	    err = testmod.errcode;
	    if (!err) {
		F = ((testmod.ess - varmod->ess) / neqns) / 
		    (varmod->ess / varmod->dfd);
	    }
	    clear_model(&testmod);
	}

	if (!err) {
	    pprintf(prn, "F(%d, %d) = %f, ", neqns, varmod->dfd, F);
	    pprintf(prn, _("p-value %f\n"), fdist(F, neqns, varmod->dfd)); 
	    if (var != NULL) {
		var->Fvals[*k] = F;
		*k += 1;
	    }
	}
    }

    pputc(prn, '\n');
    if (!err && pause) takenotes(0);

 bailout:
    free(shortlist);
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

static int real_var (int order, const int *inlist, 
		     double ***pZ, DATAINFO *pdinfo,
		     GRETL_VAR **pvar, struct var_resids *resids, 
		     PRN *prn, gretlopt opts, char flags)
{
    int i, idx, k, l, listlen, end, neqns;
    int *list = NULL, *varlist = NULL, *depvars = NULL;
    char *detlist = NULL;
    int t1, t2, oldt1, oldt2;
    int missv = 0, misst = 0;
    MODEL var_model;
    GRETL_VAR *var = NULL;
    int save = (flags & VAR_SAVE);
    int pause = (flags & VAR_PRINT_PAUSE);
    int err = 0;

    oldt1 = pdinfo->t1;
    oldt2 = pdinfo->t2;

    if (resids == NULL && order < 1) {
	fprintf(stderr, I_("Not much point in a zero-order \"VAR\" surely?\n"));
	return 1;
    }

    listlen = inlist[0] + 1;
    detlist = malloc(listlen * sizeof *detlist);
    if (detlist == NULL) {
	return E_ALLOC;
    }

    list = copylist(inlist);
    if (list == NULL) {
	free(detlist);
	return E_ALLOC;
    }

    depvars = malloc(listlen * sizeof *depvars);
    if (depvars == NULL) {
	free(detlist);
	free(list);
	return E_ALLOC;
    }

    /* how long will our list have to be? */
    listlen = get_var_listlen(list, detlist, order, depvars,
			      *pZ, pdinfo);

    varlist = gretl_list_new(listlen);
    if (varlist == NULL) {
	err = E_ALLOC;
	goto var_bailout;
    }

    neqns = depvars[0];

    /* generate the required lags */
    if (real_list_laggenr(depvars, pZ, pdinfo, order)) {
	err = E_ALLOC;
	goto var_bailout;
    }

    /* initialize */
    idx = 2; /* skip beyond the counter and the dep var */
    end = listlen;

    /* and fill out the list */
    for (i=1; i<=list[0]; i++) {
	if (detlist[i]) {
	    /* deterministic var: put at end of list */
	    varlist[end] = list[i];
	    end--;
	    continue;	    
	}
	/* otherwise it's a "real" variable and we replace it with
	   <order> lags of itself */
	if (varindex(pdinfo, pdinfo->varname[list[i]]) < pdinfo->v) {
	    for (l=1; l<=order; l++) {
		int lnum = lagvarnum(list[i], l, pdinfo);

		if (lnum > 0) {
		    varlist[idx++] = lnum; 
		} else {
		    err = E_ALLOC;
		    goto var_bailout;
		}
	    }
	}
    }

    /* sort out sample range */
    t1 = pdinfo->t1;
    t2 = pdinfo->t2;
    varlist[1] = depvars[1];

    if ((missv = adjust_t1t2(NULL, varlist, &t1, &t2, 
			     (const double **) *pZ, &misst))) {
	err = 1;
	goto var_bailout;
    }

    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    if (flags & VAR_IMPULSE_RESPONSES) {
	var = gretl_var_new(neqns, order, pdinfo, flags);
	if (var == NULL) {
	    err = E_ALLOC;
	    goto var_bailout;
	}
	if (save && var->Fvals == NULL) {
	    save = 0;
	}
    } 

    /* even in case of VAR_IMPULSE_RESPONSES this may be used */
    gretl_model_init(&var_model);

    if (flags & VAR_PRINT_MODELS) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), order);
    }

    /* apparatus for saving the residuals */
    if (resids != NULL) {
	resids->m = 2 * neqns;
	resids->uhat = malloc(2 * neqns * sizeof *resids->uhat);
	if (resids->uhat == NULL) {
	    err = E_ALLOC;
	    goto var_bailout;
	}
    }

    /* run the several regressions */
    k = 0;
    for (i=0; i<neqns; i++) {
	MODEL *pmod;

	if (flags & VAR_IMPULSE_RESPONSES) {
	    pmod = var->models[i];
	} else {
	    pmod = &var_model;
	}

	varlist[1] = depvars[i + 1];

	/* run an OLS regression for the current dependent var */
	*pmod = lsq(varlist, pZ, pdinfo, VAR, (opts | OPT_A), 0.0);
	pmod->aux = VAR;
	pmod->ID = i + 1;

	/* save the residuals if required */
	if (resids != NULL) {
	    resids->t1 = pmod->t1;
	    resids->t2 = pmod->t2;
	    resids->uhat[i] = pmod->uhat;
	    pmod->uhat = NULL;
	} 

	if (flags & VAR_PRINT_MODELS) {
	    printmodel(pmod, pdinfo, OPT_NONE, prn);
	}

	if (flags & VAR_IMPULSE_RESPONSES) {
	    /* store info in var structure */
	    err = add_model_data_to_var(var, pmod, i);
	    if (err) goto var_bailout;
	} 

	if (resids != NULL) {
	    /* estimate equations for Johansen test */
	    MODEL jmod;

	    varlist[1] = resids->levels_list[i + 1]; 
	    jmod = lsq(varlist, pZ, pdinfo, VAR, (opts | OPT_A), 0.0);
	    if (flags & VAR_PRINT_MODELS) {
		jmod.aux = VAR;
		printmodel(&jmod, pdinfo, OPT_NONE, prn);
	    }
	    resids->uhat[i + neqns] = jmod.uhat;
	    jmod.uhat = NULL;
	    clear_model(&jmod);
	}

	if (flags & VAR_DO_FTESTS) {
	    var_F_tests(pmod, (save)? var : NULL,
			varlist, depvars,
			pZ, pdinfo, order, neqns, 
			idx, &k, pause, prn);
	}

	if (1) { /* FIXME? */
	    clear_model(&var_model);
	}

    }
    pputc(prn, '\n');

 var_bailout:

    free(list);
    free(varlist);
    free(depvars);
    free(detlist);

    /* reset sample range to what it was before */
    pdinfo->t1 = oldt1;
    pdinfo->t2 = oldt2;

    if (flags & VAR_IMPULSE_RESPONSES) {
	if (!err) {
#ifdef VAR_DEBUG
	    gretl_matrix_print(var->A, "var->A", prn);
#endif
	    err = gretl_var_do_error_decomp(var);
	    if (!err) {
#ifdef VAR_DEBUG
		gretl_matrix_print(var->C, "var->C", prn);
#endif
		for (i=0; i<var->neqns; i++) {
		    /* FIXME: make horizon configurable */
		    gretl_var_print_impulse_response(var, i, 0, pdinfo, 
						     pause, prn);
		    gretl_var_print_fcast_decomp(var, i, 0, pdinfo, 
						 pause, prn);
		}
	    } else {
		fprintf(stderr, "failed: gretl_var_do_error_decomp\n");
	    }
	}
	if ((flags & VAR_SAVE) && pvar != NULL) {
	    *pvar = var;
	} else {
	    gretl_var_free(var);
	}
    }

    return err;
}

/**
 * simple_var:
 * @order: lag order for the VAR
 * @list: specification for the first model in the set.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @pause: if = 1, pause after showing each model.
 * @prn: gretl printing struct.
 *
 * Estimate a vector auto-regression (VAR) and print the results.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int simple_var (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		int pause, gretlopt opts, PRN *prn)
{
    char flags = VAR_PRINT_MODELS | VAR_DO_FTESTS;
    PRN *myprn;

#if 1
    flags |= VAR_IMPULSE_RESPONSES;
#endif

    if (pause) {
	flags |= VAR_PRINT_PAUSE;
    }

    if (opts & OPT_Q) myprn = NULL;
    else myprn = prn;

    return real_var(order, list, pZ, pdinfo, NULL, NULL, myprn, opts, flags);
}

/* "full" version returns pointer to VAR struct */

GRETL_VAR *full_var (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		     gretlopt opts, PRN *prn)
{
    GRETL_VAR *var = NULL;
    int err;

    err = real_var(order, list, pZ, pdinfo, &var, NULL, prn, opts, 
		   VAR_PRINT_MODELS | VAR_DO_FTESTS | 
		   VAR_IMPULSE_RESPONSES | VAR_SAVE);
    if (err) {
	gretl_var_free(var);
	return NULL;
    } else {
	return var;
    }
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
 * @prn: gretl printing struct.
 *
 * Test for cointegration.  
 *
 * Returns: 0 on successful completion.
 *
 */

int coint (int order, const int *list, double ***pZ, 
	   DATAINFO *pdinfo, PRN *prn)
{
    int i, t, n, nv, l0 = list[0];
    MODEL cmod;
    int *cointlist = NULL;
    int adfcode = UR_CONST;

    gretl_model_init(&cmod);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (list[i] == 0) {
	    continue;
	}
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	real_adf_test(list[i], order, 1, pZ, pdinfo, 0L, -1, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (gretl_hasconst(list) == 0) {
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
	cointlist = copylist(list);
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
    if (dataset_add_vars(1, pZ, pdinfo)) {
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
		  pZ, pdinfo, OPT_N, adfcode, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
		 "(a) The unit-root hypothesis is not rejected for the individual"
		 " variables.\n(b) The unit-root hypothesis is rejected for the "
		 "residuals (uhat) from the \n    cointegrating regression.\n"));

    /* clean up and get out */
    clear_model(&cmod);
    free(cointlist);
    dataset_drop_vars(1, pZ, pdinfo);

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

static int real_adf_test (int varno, int order, int niv,
			  double ***pZ, DATAINFO *pdinfo, 
			  gretlopt opt, int cointcode,
			  PRN *prn)
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
    int orig_nvars = pdinfo->v;
    int blurb_done = 0;
    MODEL dfmod;
    int *list;
    char pvstr[48];
    double DFt, pv;
    char mask[4] = {0};
    int i, itv;

    list = adf_prepare_vars(order, varno, pZ, pdinfo);
    if (list == NULL) {
	return E_ALLOC;
    }

    gretl_model_init(&dfmod);

    if (opt == 0L || opt == OPT_V) {
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

	list[0] = 2 + order + i;

	if (i > 0) {
	    list[list[0]] = 0;
	} 

	if (i >= 2) {
	    list[3 + order] = gettrend(pZ, pdinfo, 0);
	    if (list[3 + order] == TREND_FAILED) {
		return 1;
	    }
	}
	if (i > 2) {
	    list[4 + order] = gettrend(pZ, pdinfo, 1);
	    if (list[4 + order] == TREND_FAILED) {
		return 1;
	    }
	}	    

	dfmod = lsq(list, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (dfmod.errcode) {
	    return dfmod.errcode;
	}

	DFt = dfmod.coeff[dfnum] / dfmod.sderr[dfnum];

	if (cointcode > 0) {
	    itv = cointcode;
	} else {
	    itv = (i == 0)? UR_NO_CONST :
		(i == 1)? UR_CONST : 
		(i == 2)? UR_TREND :
		UR_TREND_SQUARED;
	}

	pv = df_pvalue_from_plugin(DFt, 
				   /* use asymptotic p-value for augmented case */
				   (order > 0)? 0 : dfmod.nobs, 
				   niv, itv);

	if (na(pv)) {
	    sprintf(pvstr, "%s %s", _("p-value"), _("unknown"));
	} else {
	    sprintf(pvstr, "%s %.4g", 
		    (order > 0)? _("asymptotic p-value") : _("p-value"), 
		    pv);
	} 

	if (!blurb_done) {
	    if (order > 0) {
		pprintf(prn, _("\nAugmented Dickey-Fuller tests, order %d, for %s\n"),
			order, pdinfo->varname[varno]);
	    } else {
		pprintf(prn, _("\nDickey-Fuller tests for %s\n"),
			pdinfo->varname[varno]);
	    }
	    pprintf(prn, _("sample size %d\n"), dfmod.nobs);
	    pputs(prn, _("unit-root null hypothesis: a = 1"));
	    pputs(prn, "\n\n");
	    blurb_done = 1;
	}

	pprintf(prn, "   %s\n", _(teststrs[i]));
	if (cointcode == 0) {
	    pprintf(prn, "   %s: %s\n", _("model"), 
		    (order > 0)? aug_models[i] : models[i]);
	}
	pprintf(prn, "   %s: %g\n"
		"   %s: t = %g\n"
		"   %s\n",
		_("estimated value of (a - 1)"), dfmod.coeff[dfnum],
		_("test statistic"), DFt,
		pvstr);	

	if (opt & OPT_V) {
	    /* verbose */
	    dfmod.aux = (order > 0)? AUX_ADF : AUX_DF;
	    if (!na(pv)) {
		gretl_model_set_int(&dfmod, "dfnum", dfnum + 2);
		gretl_model_set_double(&dfmod, "dfpval", pv);
	    }
	    printmodel(&dfmod, pdinfo, OPT_NONE, prn);
	} else {
	    pputc(prn, '\n');
	}

	clear_model(&dfmod);
    }

    if (cointcode >= 0) {
	pputs(prn, _("P-values based on MacKinnon (JAE, 1996)\n"));
    }

    free(list);
    dataset_drop_vars(pdinfo->v - orig_nvars, pZ, pdinfo);

    return 0;
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
    return real_adf_test(varno, order, 1, pZ, pdinfo, opt, 0, prn);
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
	if (tprob(tstat, tmod.dfd) < 0.05) {
	    trends = 1;
	}

	clear_model(&tmod);

	if (trends) break;
    }

    dataset_drop_vars(pdinfo->v - origv, pZ, pdinfo);

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

int johansen_test (int order, const int *list, double ***pZ, DATAINFO *pdinfo,
		   gretlopt opt, PRN *prn)
{
    PRN *varprn = NULL;
    struct var_resids resids;
    char flags;
    int err = 0;
    int i, j;
    int orig_t1 = pdinfo->t1;
    int orig_v = pdinfo->v;
    int *varlist;
    int verbose = (opt & OPT_O);
    int hasconst = 0;
    int trends = 0;

    /* we're assuming that the list we are fed is in levels */
    resids.levels_list = malloc((1 + list[0]) * sizeof *list);
    if (resids.levels_list == NULL) {
	return E_ALLOC;
    }

    varlist = malloc((2 + list[0]) * sizeof *list);
    if (varlist == NULL) {
	free(resids.levels_list);
	return E_ALLOC;
    }

    resids.levels_list[0] = varlist[0] = list[0];

    j = 1;
    for (i=1; i<=list[0]; i++) {
	int lnum;

	if (list[i] == 0) {
	    resids.levels_list[0] -= 1;
	    hasconst = 1;
	    continue;
	}
	lnum = laggenr(list[i], 1, pZ, pdinfo);
	if (lnum < 0) {
	    free(varlist);
	    free(resids.levels_list);
	    return E_DATA;
	}
	resids.levels_list[j++] = lnum;
    }

    /* now get differences and put them into list */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) {
	    continue;
	}
	varlist[i] = diffgenr(list[i], pZ, pdinfo, 0);
	if (varlist[i] < 0) {
	    free(varlist);
	    free(resids.levels_list);
	    return E_DATA;
	}
    }

    if (!hasconst) {
	varlist[0] += 1;
	varlist[varlist[0]] = 0;
    }

    if (verbose) {
	flags = VAR_PRINT_MODELS;
	varprn = prn;
    } else {
	flags = 0;
	varprn = gretl_print_new(GRETL_PRINT_NULL, NULL);
    }

    /* FIXME? */
    pdinfo->t1 += (order + 1);
    /* Check Hamilton: what if order for test = 1? */
    err = real_var(order - 1, varlist, pZ, pdinfo, 
		   NULL, &resids, varprn, OPT_NONE, flags); 
    
    if (!verbose) {
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

	if (verbose) {
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

    dataset_drop_vars(pdinfo->v - orig_v, pZ, pdinfo);

    return err;
}

int gretl_var_print (GRETL_VAR *var, const DATAINFO *pdinfo, PRN *prn)
{
    int i, j, k, v;
    int dfd = (var->models[0])->dfd;

    if (prn == NULL) return 0;

    if (TEX_PRN(prn)) {
	pputs(prn, "\\noindent");
	pprintf(prn, I_("\nVAR system, lag order %d\n\n"), var->order);
    } else {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), var->order);
    }

    k = 0;
    for (i=0; i<var->neqns; i++) {

	printmodel(var->models[i], pdinfo, OPT_NONE, prn);

	if (var->Fvals == NULL) continue;

	if (TEX_PRN(prn)) {
	    pputs(prn, "\n\\begin{center}\n");
	    pprintf(prn, "%s\\\\[1em]\n", I_("F-tests of zero restrictions"));
	    pputs(prn, "\\begin{tabular}{lll}\n");
	} else {
	    pputs(prn, _("\nF-tests of zero restrictions:\n\n"));
	}

	for (j=0; j<var->neqns; j++) {
	    v = (var->models[j])->list[1];
	    if (TEX_PRN(prn)) {
		pprintf(prn, I_("All lags of %-8s "), pdinfo->varname[v]);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, I_("p-value %f"), fdist(var->Fvals[k], var->order, dfd));
		pputs(prn, "\\\\\n");
	    } else {
		pprintf(prn, _("All lags of %-8s "), pdinfo->varname[v]);
		pprintf(prn, "F(%d, %d) = %g, ", var->order, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->order, dfd));
	    }
	    k++;
	}

	if (var->order > 1) {
	    if (TEX_PRN(prn)) {
		pprintf(prn, I_("All vars, lag %-6d "), var->order);
		pputs(prn, "& ");
		pprintf(prn, "$F(%d, %d) = %g$ & ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->neqns, dfd));
	    } else {
		pprintf(prn, _("All vars, lag %-6d "), var->order);
		pprintf(prn, "F(%d, %d) = %g, ", var->neqns, dfd, var->Fvals[k]);
		pprintf(prn, _("p-value %f\n"), fdist(var->Fvals[k], var->neqns, dfd));
	    } 
	    k++;
	}
	if (TEX_PRN(prn)) {
	    pputs(prn, "\\end{tabular}\n"
		  "\\end{center}\n\n"
		  "\\clearpage\n\n");
	} else {
	    pputc(prn, '\n');
	}
    }

    pputc(prn, '\n');

    for (i=0; i<var->neqns; i++) {
	gretl_var_print_impulse_response(var, i, 0, pdinfo, 0, prn);
	gretl_var_print_fcast_decomp(var, i, 0, pdinfo, 0, prn);
    }

    return 0;
}

void gretl_var_assign_name (GRETL_VAR *var)
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

void gretl_var_assign_specific_name (GRETL_VAR *var, const char *name)
{
    if (var->name != NULL) {
	free(var->name);
    }
    var->name = malloc(strlen(name) + 1);
    if (var->name != NULL) {
	strcpy(var->name, name);
    }    
}

const char *gretl_var_get_name (const GRETL_VAR *var)
{
    return var->name;
}

int gretl_var_add_resids_to_dataset (GRETL_VAR *var, int eqnum,
				     double ***pZ, DATAINFO *pdinfo)
{
    MODEL *pmod = var->models[eqnum];
    int i, t;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

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

int gretl_var_get_highest_variable (const GRETL_VAR *var,
				    const DATAINFO *pdinfo)
{
    int vmax = 0;

    if (var->models != NULL && var->neqns >= 1) {
	vmax = highest_numbered_var_in_model(var->models[0], pdinfo);
    }

    return vmax;
}

