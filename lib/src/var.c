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
#include "internal.h"

/* #define VAR_DEBUG */

struct _GRETL_VAR {
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
	    var->models[i] = gretl_model_new(pdinfo);
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
	    clear_model(var->models[i], NULL);
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
	if (pdinfo->pd == 4) periods = 20;
	else if (pdinfo->pd == 12) periods = 24;
	else periods = 10;
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
	    page_break(0, NULL, 0);
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
	if (pdinfo->pd == 4) periods = 20;
	else if (pdinfo->pd == 12) periods = 24;
	else periods = 10;
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
	    page_break(0, NULL, 0);
	}
    }

    if (vd != NULL) gretl_matrix_free(vd);

    return err;
}

/* ......................................................  */

static int gettrend (double ***pZ, DATAINFO *pdinfo)
{
    int index;
    int t, n = pdinfo->n, v = pdinfo->v;

    index = varindex(pdinfo, "time");
    if (index < v) return index;
    
    if (dataset_add_vars(1, pZ, pdinfo)) return 999;

    for (t=0; t<n; t++) (*pZ)[v][t] = (double) t+1;
    strcpy(pdinfo->varname[v], "time");
    strcpy(VARLABEL(pdinfo, v), _("time trend variable"));
	    
    return index;
}

/* ...................................................................  */

static int diffvarnum (int index, const DATAINFO *pdinfo)
     /* Given an "ordinary" variable name, construct the name of the
	corresponding first difference and find its ID number */
{
    char diffname[16], s[16];
    
    strcpy(s, pdinfo->varname[index]);
    _esl_trunc(s, 6);
    strcpy(diffname, "d_");
    strcat(diffname, s);
    return varindex(pdinfo, diffname);
}

/* ......................................................  */

static int diffgenr (int iv, double ***pZ, DATAINFO *pdinfo)
{
    char word[32];
    char s[32];
    int t, t1, n = pdinfo->n, v = pdinfo->v;
    double x0, x1;

    strcpy(word, pdinfo->varname[iv]);
    _esl_trunc(word, 6);
    strcpy(s, "d_");
    strcat(s, word);

    /* "s" should now contain the new variable name --
     check whether it already exists: if so, get out */
    if (varindex(pdinfo, s) < v) return 0;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->time_series == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}	
	x0 = (*pZ)[iv][t];
	x1 = (*pZ)[iv][t-1];
	if (na(x0) || na(x1)) {
	    (*pZ)[v][t] = NADBL;
	} else {				      
	    (*pZ)[v][t] = x0 - x1;
	}
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(VARLABEL(pdinfo, v), _("%s = first difference of %s"),
	    pdinfo->varname[v], pdinfo->varname[iv]);
	    
    return 0;
}

/* ......................................................  */

static int ldiffgenr (int iv, double ***pZ, DATAINFO *pdinfo)
{
    char word[32];
    char s[32];
    int t, t1, n = pdinfo->n, v = pdinfo->v;
    double x0, x1;

    strcpy(word, pdinfo->varname[iv]);
    _esl_trunc(word, 5);
    strcpy(s, "ld_");
    strcat(s, word);

    /* "s" should now contain the new variable name --
     check whether it already exists: if so, get out */
    if (varindex(pdinfo, s) < v) return 0;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    for (t=0; t<n; t++) (*pZ)[v][t] = NADBL;
    t1 = (pdinfo->t1 > 1)? pdinfo->t1 : 1;
    for (t=t1; t<=pdinfo->t2; t++) {
	if (pdinfo->time_series == STACKED_TIME_SERIES &&
	    panel_unit_first_obs(t, pdinfo)) {
	    continue;
	}
	x0 = (*pZ)[iv][t];
	x1 = (*pZ)[iv][t-1];
	if (na(x0) || na(x1) || x0 / x1 <= 0.) {
	    (*pZ)[v][t] = NADBL;
	} else {			      
	    (*pZ)[v][t] = log(x0 / x1);
	}
    }

    strcpy(pdinfo->varname[v], s);
    sprintf(VARLABEL(pdinfo, v), _("%s = log difference of %s"),
	    pdinfo->varname[v], pdinfo->varname[iv]);
	    
    return 0;
}

/**
 * list_diffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate first-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_diffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++) {
	if (diffgenr(list[i], pZ, pdinfo)) return 1;
    }
    return 0;
}

/**
 * list_ldiffgenr:
 * @list: list of variables to process.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 *
 * Generate log-differences of the variables in @list, and add them
 * to the data set.
 *
 * Returns: 0 on successful completion, 1 on error.
 *
 */

int list_ldiffgenr (const LIST list, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    
    for (i=1; i<=list[0]; i++) {
	if (ldiffgenr(list[i], pZ, pdinfo)) return 1;
    }
    return 0;
}

/**
 * lagvarnum:
 * @iv: ID number of the variable.
 * @lag: Desired lag length.
 * @pdinfo: data information struct.
 *
 * Given an "ordinary" variable, construct the name of the
 * corresponding lagged variable and find its ID number.
 *
 * Returns: the ID number of the lagged variable.
 *
 */

static int lagvarnum (int iv, int lag, const DATAINFO *pdinfo)
{
    char lagname[16], ext[6];

    strcpy(lagname, pdinfo->varname[iv]);

    if (pdinfo->pd >=10) _esl_trunc(lagname, 5);
    else _esl_trunc(lagname, 6);

    sprintf(ext, "_%d", lag);
    strcat(lagname, ext);

    return varindex(pdinfo, lagname);
}

/* ...................................................................  */

static void reset_list (int *list1, int *list2)
{
    int i;
    
    for (i=2; i<=list1[0]; i++) list1[i] = list2[i];
}

/* ...................................................................  */

static int get_listlen (int *varlist, char *detlist, int order, 
			double **Z, const DATAINFO *pdinfo)
     /* parse varlist (for a VAR) and determine how long the augmented 
	list will be, once all the appropriate lag terms are inserted */
{
    int i, j = 1, v = 1;
    int gotsep = 0;

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
	}
	varlist[j++] = varlist[i];
    }

    if (gotsep) varlist[0] -= 1;
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

/* ...................................................................  */

static int real_var (int order, const LIST inlist, 
		     double ***pZ, DATAINFO *pdinfo,
		     GRETL_VAR **pvar, struct var_resids *resids, 
		     PRN *prn, char flags)
{
    /* construct the respective lists by adding the appropriate
       number of lags ("order") to the variables in list 

       Say the list is "x_1 const time x_2 x_3", and the order is 2.
       Then the first list should be

       x_1 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

       the second:

       x_2 const time x_1(-1) x_1(-2) x_2(-1) x_2(-2) x_3(-1) x_3(-2)

       and so on.

       Run the regressions and print the results.
    */

    int i, j, index, k, l, listlen, end, neqns = 0;
    int *list = NULL, *varlist = NULL, *depvars = NULL, *shortlist = NULL;
    char *detlist = NULL;
    int t1, t2, oldt1, oldt2, dfd = 0;
    int missv = 0, misst = 0;
    double essu = 0.0, F = 0.0;
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

    detlist = malloc(inlist[0] + 1);
    if (detlist == NULL) return E_ALLOC;

    if (copylist(&list, inlist)) {
	free(detlist);
	return E_ALLOC;
    }

    /* how long will our list have to be? */
    listlen = get_listlen(list, detlist, order, *pZ, pdinfo);

    varlist = malloc((listlen + 1) * sizeof *varlist);
    depvars = malloc((listlen + 1) * sizeof *depvars);
    shortlist = malloc(listlen * sizeof *shortlist);

    if (varlist == NULL || depvars == NULL || shortlist == NULL) {
	err = E_ALLOC;
	goto var_bailout;
    }

    varlist[0] = listlen;
    index = 2; /* skip beyond the counter and the dep var */
    end = listlen;

    /* now fill out the list */
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
	    depvars[neqns] = list[i];
	    neqns++;
	    for (l=1; l<=order; l++) {
		int lnum = laggenr(list[i], l, 1, pZ, pdinfo);

		if (lnum > 0) {
		    varlist[index] = lnum; 
		    index++;
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
    varlist[1] = depvars[0];

    if ((missv = _adjust_t1t2(NULL, varlist, &t1, &t2, 
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
    gretl_model_init(&var_model, pdinfo);

    if (flags & VAR_PRINT_MODELS) {
	pprintf(prn, _("\nVAR system, lag order %d\n\n"), order);
    }
    shortlist[0] = listlen - order;

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

	varlist[1] = depvars[i];

	/* run an OLS regression for the current dependent var */
	*pmod = lsq(varlist, pZ, pdinfo, VAR, 1, 0.0);
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
	    printmodel(pmod, pdinfo, prn);
	}

	if (flags & VAR_DO_FTESTS) {
	    /* keep some results for hypothesis testing */
	    essu = pmod->ess;
	    dfd = pmod->dfd;	    
	}

	if (flags & VAR_IMPULSE_RESPONSES) {
	    /* store info in var structure */
	    err = add_model_data_to_var(var, pmod, i);
	    if (err) goto var_bailout;
	} else {
	    clear_model(&var_model, pdinfo);
	}

	if (resids != NULL) {
	    /* estimate equations for Johansen test too (use var_model) */
	    varlist[1] = resids->levels_list[i + 1]; 
	    var_model = lsq(varlist, pZ, pdinfo, VAR, 0, 0.0);
	    if (flags & VAR_PRINT_MODELS) {
		var_model.aux = VAR;
		printmodel(&var_model, pdinfo, prn);
	    }
	    resids->uhat[i + neqns] = var_model.uhat;
	    var_model.uhat = NULL;
	    clear_model(&var_model, pdinfo);
	}

	if (flags & VAR_DO_FTESTS) {
	    /* now build truncated lists for hypothesis tests */
	    shortlist[1] = varlist[1];
	    pputs(prn, _("\nF-tests of zero restrictions:\n\n"));

	    for (j=0; j<neqns; j++) {
		reset_list(shortlist, varlist);
		for (l=1; l<=order; l++) {
		    index = l + 1 + j * order;
		    if (index > shortlist[0]) break;
		    shortlist[index] = varlist[index + order];
		}
		end = 0;
		for (l=shortlist[0]; l>index; l--) {
		    shortlist[l] = varlist[varlist[0] - end];
		    end++;
		}
		pprintf(prn, _("All lags of %-8s "), 
			pdinfo->varname[depvars[j]]);
		var_model = lsq(shortlist, pZ, pdinfo, VAR, 0, 0.0);
		F = ((var_model.ess - essu) / order) / (essu / dfd);
		clear_model(&var_model, pdinfo);
		pprintf(prn, "F(%d, %d) = %f, ", order, dfd, F);
		pprintf(prn, _("p-value %f\n"), fdist(F, order, dfd));
		if (save) var->Fvals[k++] = F;
	    }

	    if (order > 1) {
		pprintf(prn, _("All vars, lag %-6d "), order);
		reset_list(shortlist, varlist);
		index = 2;
		for (j=1; j<=neqns*(order); j++) {
		    if (j % order) {
			shortlist[index] = varlist[j+1];
			index++;
		    }
		}
		end = 0;
		for (l=shortlist[0]; l>=index; l--) {
		    shortlist[l] = varlist[varlist[0]-end];
		    end++;
		}
		var_model = lsq(shortlist, pZ, pdinfo, VAR, 0, 0.0);
		F = ((var_model.ess - essu) / neqns) / (essu / dfd);
		clear_model(&var_model, pdinfo);
		pprintf(prn, "F(%d, %d) = %f, ", neqns, dfd, F);
		pprintf(prn, _("p-value %f\n"), fdist(F, neqns, dfd)); 
		if (save) var->Fvals[k++] = F;
	    }

	    pputc(prn, '\n');
	    if (pause) page_break(0, NULL, 0);
	}
    }
    pputc(prn, '\n');

 var_bailout:

    free(list);
    free(varlist);
    free(shortlist);
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

int simple_var (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
		int pause, PRN *prn)
{
    char flags = VAR_PRINT_MODELS | VAR_DO_FTESTS;

#if 1
    flags |= VAR_IMPULSE_RESPONSES;
#endif

    if (pause) {
	flags |= VAR_PRINT_PAUSE;
    }

    return real_var(order, list, pZ, pdinfo, NULL, NULL, prn, flags);
}

/* "full" version returns pointer to VAR struct */

GRETL_VAR *full_var (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
		     PRN *prn)
{
    GRETL_VAR *var = NULL;
    int err;

    err = real_var(order, list, pZ, pdinfo, &var, NULL, prn,
		   VAR_PRINT_MODELS | VAR_DO_FTESTS | 
		   VAR_IMPULSE_RESPONSES | VAR_SAVE);
    if (err) {
	gretl_var_free(var);
	return NULL;
    } else {
	return var;
    }
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

int coint (int order, const LIST list, double ***pZ, 
	   DATAINFO *pdinfo, PRN *prn)
     /* FIXME - need proper error checking here */
{
    int i, t, n, nv, l0 = list[0];
    MODEL coint_model;
    int *cointlist = NULL;

    gretl_model_init(&coint_model, pdinfo);

    /* step 1: test all the vars for unit root */
    for (i=1; i<=l0; i++) {
	if (i > 1) pputc(prn, '\n');
	pprintf(prn, _("Step %d: testing for a unit root in %s\n"),
		i, pdinfo->varname[list[i]]);
	adf_test(order, list[i], pZ, pdinfo, prn);
    }

    /* step 2: carry out the cointegrating regression */
    if (_hasconst(list) == 0) {
	cointlist = malloc((l0 + 2) * sizeof *cointlist);
	if (cointlist == NULL) return E_ALLOC;
	for (i=0; i<=l0; i++) cointlist[i] = list[i];
	cointlist[l0 + 1] = 0;
	cointlist[0] += 1;
    } else {
	copylist(&cointlist, list);
    }

    pputc(prn, '\n');
    pprintf(prn, _("Step %d: cointegration\n"), l0 + 1);
    
    coint_model = lsq(cointlist, pZ, pdinfo, OLS, 1, 0.0); 
    coint_model.aux = AUX_COINT;
    printmodel(&coint_model, pdinfo, prn);

    /* add residuals from cointegrating regression to data set */
    n = pdinfo->n;
    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
    nv = pdinfo->v - 1;

    for (t=0; t<coint_model.t1; t++) {
	(*pZ)[nv][t] = NADBL;
    }
    for (t = coint_model.t1; t<=coint_model.t2; t++) {
	(*pZ)[nv][t] = coint_model.uhat[t];
    }
    for (t=coint_model.t2 + 1; t<n; t++) {
	(*pZ)[nv][t] = NADBL;
    }

    strcpy(pdinfo->varname[nv], "uhat");

    /* Run ADF test on these residuals */
    pputc(prn, '\n');
    adf_test(order, pdinfo->v - 1, pZ, pdinfo, prn);

    pputs(prn, _("\nThere is evidence for a cointegrating relationship if:\n"
	    "(a) The unit-root hypothesis is not rejected for the individual"
	    " variables.\n(b) The unit-root hypothesis is rejected for the "
	    "residuals (uhat) from the \n    cointegrating regression.\n"
	    "\n(Note that significance levels for the D-W and F statistics here "
	    "cannot be \nread from the usual statistical tables.)\n"));

    /* clean up and get out */
    clear_model(&coint_model, pdinfo);
    free(cointlist);
    dataset_drop_vars(1, pZ, pdinfo);

    return 0;
}

/**
 * adf_test:
 * @order: lag order for the test.
 * @varno: ID number of the variable to test.
 * @pZ: pointer to data matrix.
 * @pdinfo: data information struct.
 * @prn: gretl printing struct.
 *
 * Carries out and prints the results of the Augmented Dickey-Fuller test for 
 * a unit root.
 *
 * Returns: 0 on successful completion, non-zero on error.
 *
 */

int adf_test (int order, int varno, double ***pZ,
	      DATAINFO *pdinfo, PRN *prn)
{
    int i, l, T, k, row, orig_nvars = pdinfo->v;
    int *adflist = NULL;
    int *shortlist = NULL;
    MODEL adf_model;
    double essu, F, DFt;
    char pval[40];

                                 /* 1%    2.5%    5%    10% */
    double t_crit_vals[6][4] = {{ -3.75, -3.33, -3.00, -2.62 },  /* T=25 */
				{ -3.58, -3.22, -2.93, -2.60 },  /* T=50 */
				{ -3.51, -3.17, -2.89, -2.58 },  /* T=100 */
				{ -3.46, -3.14, -2.88, -2.57 },  /* T=250 */
				{ -3.44, -3.13, -2.87, -2.57 },  /* T=500 */
				{ -3.43, -3.12, -2.86, -2.57 }}; /* inf */

                                /* 1%   2.5%   5%    10% */
    double f_crit_vals[6][4] = {{ 5.91, 7.24, 8.65, 10.61 },  /* T = 25 */
 			        { 5.61, 6.73, 7.81,  9.31 },  /* T = 50 */
			        { 5.47, 6.49, 7.44,  8.73 },  /* T = 100 */
			        { 5.39, 6.34, 7.25,  8.43 },  /* T = 250 */
			        { 5.36, 6.30, 7.20,  8.34 },  /* T = 500 */
			        { 5.34, 6.25, 7.16,  8.27 }}; /* inf */
    

    if (varno == 0) return E_DATA;

    gretl_model_init(&adf_model, pdinfo);
    k = 3 + order;

    adflist = malloc((5 + order) * sizeof *adflist);
    shortlist = malloc(k * sizeof *shortlist);

    if (adflist == NULL || shortlist == NULL) {
	free(adflist);
	free(shortlist);
	return E_ALLOC;
    }

    i = pdinfo->t1;
    pdinfo->t1 = 0;
    diffgenr(varno, pZ, pdinfo);
    if (laggenr(varno, 1, 1, pZ, pdinfo) < 0) {
	free(adflist);
	free(shortlist);
	return E_DATA;
    }
    pdinfo->t1 = i;

    adflist[1] = diffvarnum(varno, pdinfo);

    /* do the more familiar Dickey-Fuller t-test first */
    adflist[0] = 3;
    adflist[2] = 0;
    adflist[3] = lagvarnum(varno, 1, pdinfo);

    adf_model = lsq(adflist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode) {
	return adf_model.errcode;
    }

    DFt = adf_model.coeff[1] / adf_model.sderr[1];
    T = adf_model.nobs;

    row = (T > 500)? 5 : 
	(T > 450)? 4 : 
	(T > 240)? 3 : 
	(T > 90)? 2 : 
	(T > 40)? 1 : 
	(T > 24)? 0 : -1;

    if (row < 0) {
	sprintf(pval, _("significance level unknown"));
    } else {
	if (DFt < t_crit_vals[row][0])
	    sprintf(pval, _("significant at the 1 percent level"));
	else if (DFt < t_crit_vals[row][1])
	    sprintf(pval, _("significant at the 2.5 percent level"));
	else if (DFt < t_crit_vals[row][2])
	    sprintf(pval, _("significant at the 5 percent level"));
	else if (DFt < t_crit_vals[row][3])
	    sprintf(pval, _("significant at the 10 percent level"));
	else
	    sprintf(pval, _("not significant at the 10 percent level"));
    }
    
    pprintf(prn, _("\nDickey-Fuller test with constant\n\n"
	    "   model: (1 - L)%s = m + g * %s(-1) + e\n"
	    "   unit-root null hypothesis: g = 0\n"
	    "   estimated value of g: %f\n"
	    "   test statistic: t = %f, with sample size %d\n"
	    "   %s\n"),
	    pdinfo->varname[varno], pdinfo->varname[varno],
	    adf_model.coeff[1], DFt, adf_model.nobs, pval);

    clear_model(&adf_model, pdinfo);

    /* then do ADF test using F-statistic */
    adflist[0] = 4 + order;
    adflist[3] = lagvarnum(varno, 1, pdinfo);

    for (l=1; l<=order; l++) {
	int lnum = laggenr(adflist[1], l, 1, pZ, pdinfo);

	/* FIXME: handle laggenr error */
	if (lnum > 0) {
	    adflist[l+3] = lnum;
	} 
    }

    adflist[adflist[0]] = 0;
    if ((adflist[2] = gettrend(pZ, pdinfo)) == 999) {
	free(adflist);
	free(shortlist);
	return E_ALLOC;
    }

    adf_model = lsq(adflist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode)
	return adf_model.errcode;
    adf_model.aux = AUX_ADF;
    printmodel(&adf_model, pdinfo, prn);
    essu = adf_model.ess;
    T = adf_model.nobs;
    clear_model(&adf_model, pdinfo);

    shortlist[0] = adflist[0] - 2;
    shortlist[1] = adflist[1];
    for (i=0; i<=order; i++) {
	shortlist[2+i] = adflist[4+i];
    }

    adf_model = lsq(shortlist, pZ, pdinfo, OLS, 0, 0.0);
    if (adf_model.errcode) {
	/* FIXME: clean up */
	return adf_model.errcode;
    }	

    F = (adf_model.ess - essu) * (T - k)/(2 * essu);
    clear_model(&adf_model, pdinfo);

    row = -1;
    if (T > 500) row = 5;
    else if (T > 250) row = 4;
    else if (T > 100) row = 3;
    else if (T > 50) row = 2;
    else if (T > 25) row = 1;
    else if (T > 23) row = 0;

    if (row == -1) strcpy(pval, _("unknown pvalue"));
    else {
	if (F > f_crit_vals[row][3]) strcpy(pval, _("pvalue < .01"));
	else if (F > f_crit_vals[row][2]) strcpy(pval, _(".025 > pvalue > .01"));
	else if (F > f_crit_vals[row][1]) strcpy(pval, _(".05 > pvalue > .025"));
	else if (F > f_crit_vals[row][0]) strcpy(pval, _(".10 > pvalue > .05"));
	else strcpy(pval, _("pvalue > .10"));
    }

    pprintf(prn, _("Augmented Dickey-Fuller test on %s:\n   F(2, %d) = %f, "
	   "with %s\n"), pdinfo->varname[varno], T - k, F, pval);
    pprintf(prn, _("The null hypothesis is that %s has a unit root, i.e. "
	    "the parameters on\nthe time trend and %s are both zero.\n"),
	    pdinfo->varname[varno], pdinfo->varname[adflist[3]]);

    free(adflist);
    free(shortlist);
    dataset_drop_vars(pdinfo->v - orig_nvars, pZ, pdinfo);

    return 0;
}

/* ....................................................... */

#ifdef notyet

int ma_model (LIST list, double ***pZ, DATAINFO *pdinfo, PRN *prn)
{
    int t, v = pdinfo->v, err = 0;
    int malist[4], iv = list[2];
    double a, aopt, essmin, diff;
    int step, t0 = pdinfo->t1, T = pdinfo->t2;
    MODEL mamod;

    if (list[0] != 2) {
	pputs(prn, "mvavg: takes a list of two variables\n");
	return 1;
    }
    
    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;
    strcpy(pdinfo->varname[v], "Z_t");

    malist[0] = 3;
    malist[1] = list[1]; /* original dependent variable */
    malist[2] = v;       /* new var: moving average of indep var */
    malist[3] = 0;

    gretl_model_init(&mamod, pdinfo);

    a = aopt = 0.0;
    essmin = 0.0;
    diff = 0.01;
    for (step=1; step<=100; step++) {
	a += diff;
	if (a > 0.995) break;
	(*pZ)[v][t0] = (*pZ)[iv][t0] / (1 - a);
	for (t=t0+1; t<T; t++) { 
	    (*pZ)[v][t] = (*pZ)[iv][t] + a * (*pZ)[v][t-1];
	    /*  printf("newvars[%d] %g %g\n", t, 
		(*pZ)[v*n + t], (*pZ)[(v+1)*n + t]); */
	}
	clear_model(&mamod, pdinfo);
	mamod = lsq(malist, pZ, pdinfo, OLS, 0, 0.0);
	if ((err = mamod.errcode)) {
	    clear_model(&mamod, pdinfo);
	    return err;
	}	
	if (step == 1) {
	    pputs(prn, "\n ADJ       ESS      ADJ       ESS      "
		    "ADJ       ESS      ADJ       ESS     \n");
	}
	pprintf(prn, "%5.2f %10.4g", a, mamod.ess);
	if (step%4 == 0) pputc(prn, '\n');
	else _bufspace(3, prn);
	if (step == 1 || mamod.ess < essmin) {
	    essmin = mamod.ess;
	    aopt = a;
	}
    }

    pprintf(prn, "\n\nESS is minimum for adj = %.2f\n\n", aopt);
    a = aopt;
    (*pZ)[v][t0] = (*pZ)[iv][t0] / (1 - a);
    for (t=t0+1; t<T; t++) { 
	(*pZ)[v][t] = (*pZ)[iv][t] + a * (*pZ)[v][t-1];
    }
    mamod = lsq(malist, pZ, pdinfo, OLS, 1, 0.0);
    printmodel(&mamod, pdinfo, prn);

    pputs(prn, "\nEstimates of original parameters:\n");
    pprintf(prn, "constant: %.4g\n", mamod.coeff[0]);
    pprintf(prn, "slope:    %.4g\n", mamod.coeff[1] / (1 - a));
    pprintf(prn, "adaptive coefficient: %.2f\n", a);
	   
    clear_model(&mamod, pdinfo);

    dataset_drop_vars(1, pZ, pdinfo);

    return 0;
}

#endif /* notyet */

static int 
has_time_trend (LIST varlist, double ***pZ, DATAINFO *pdinfo)
{
    int i;
    int origv = pdinfo->v;
    int tlist[4];
    int trends = 0;
    MODEL tmod;

    gretl_model_init(&tmod, pdinfo);

    tlist[0] = 3;
    tlist[2] = 0;

    for (i=1; i<=varlist[0]; i++) {
	double tstat;
	int v, vl;

	v = varlist[i];
	if (v == 0) continue;
	if (diffgenr(v, pZ, pdinfo)) {
	    trends = -1;
	    break;
	}
	vl = lagvarnum(v, 1, pdinfo);
	tlist[1] = v;
	tlist[3] = vl;
	tmod = lsq(tlist, pZ, pdinfo, OLS, 0, 0.0);
	if (tmod.errcode) {
	    trends = -1;
	    clear_model(&tmod, pdinfo);
	    break;
	}
	tstat = tmod.coeff[0] / tmod.sderr[0];
	if (tprob(tstat, tmod.dfd) < 0.05) {
	    trends = 1;
	}
	clear_model(&tmod, pdinfo);
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

int johansen_test (int order, const LIST list, double ***pZ, DATAINFO *pdinfo,
		   unsigned long opt, PRN *prn)
{
    PRN *varprn = NULL;
    struct var_resids resids;
    char flags;
    int err = 0;
    int i, j;
    int orig_t1 = pdinfo->t1;
    int orig_v = pdinfo->v;
    int *varlist;
    int verbose = (opt & OPT_V);
    int hasconst = 0;
    int trends = 0;

    /* we're assuming that the list we are fed is in levels */
    resids.levels_list = malloc((1 + list[0]) * sizeof *list);
    if (resids.levels_list == NULL) return E_ALLOC;
    resids.levels_list[0] = list[0];

    varlist = malloc((2 + list[0]) * sizeof *list);
    if (varlist == NULL) return E_ALLOC;
    varlist[0] = list[0];

    j = 1;
    for (i=1; i<=list[0]; i++) {
	int lnum;

	if (list[i] == 0) {
	    resids.levels_list[0] -= 1;
	    hasconst = 1;
	    continue;
	}
	lnum = laggenr(list[i], 1, 1, pZ, pdinfo);
	/* FIXME: handle laggenr error */
	if (lnum > 0) {
	    resids.levels_list[j++] = lnum;
	}
    }

    /* now get differences and put them into list */
    for (i=1; i<=list[0]; i++) {
	if (list[i] == 0) continue;
	diffgenr(list[i], pZ, pdinfo);
	varlist[i] = diffvarnum(list[i], pdinfo);
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
		   NULL, &resids, varprn, flags); 
    
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
	    int t;

	    pprintf(prn, "Residuals from VAR model %d\n", i);
	    for (t=resids.t1; t<=resids.t2; t++) {
		ntodate(date, t, pdinfo);
		pprintf(prn, "%8s %#.*g\n", ntodate(stobs, t, pdinfo), 
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

	printmodel(var->models[i], pdinfo, prn);

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

    if (var->name != NULL) free(var->name);
    var->name = malloc(8);
    if (var->name != NULL) {
	sprintf(var->name, "%s %d", _("VAR"), ++n);
    }
}

void gretl_var_assign_specific_name (GRETL_VAR *var, const char *name)
{
    if (var->name != NULL) free(var->name);
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
    char vname[VNAMELEN], vlabel[MAXLABEL];
    MODEL *pmod = var->models[eqnum];
    int i, n, t, t1 = pmod->t1, t2 = pmod->t2;

    if (dataset_add_vars(1, pZ, pdinfo)) return E_ALLOC;

    i = pdinfo->v - 1;
    n = pdinfo->n;

    if (pmod->data != NULL) t2 += get_misscount(pmod);

    for (t=0; t<t1; t++) (*pZ)[i][t] = NADBL;
    for (t=t2+1; t<n; t++) (*pZ)[i][t] = NADBL;

    sprintf(vname, "uhat%d", eqnum + 1);
    sprintf(vlabel, _("residual from VAR system, equation %d"), eqnum + 1);

    for (t=t1; t<=t2; t++) {
	(*pZ)[i][t] = pmod->uhat[t];
    }

    strcpy(pdinfo->varname[i], vname);
    strcpy(VARLABEL(pdinfo, i), vlabel);

    return 0;
}
