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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "libgretl.h"
#include "gretl_private.h"

struct _gretl_equation_system {
    int type;
    int n_equations;
    char flags;
    int **lists;
    int *endog_vars;
    int *exog_vars;
};

const char *gretl_system_type_strings[] = {
    "sur",
    "3sls",
    "fiml",
    NULL
};

const char *gretl_system_short_strings[] = {
    N_("SUR"),
    N_("3SLS"),
    N_("FIML"),
    NULL
};

const char *gretl_system_long_strings[] = {
    N_("Seemingly Unrelated Regressions"),
    N_("Three-Stage Least Squares"),
    N_("Full Information Maximum Likelihood"),
    NULL
};

const char *nosystem = N_("No system of equations has been defined");
const char *badsystem = N_("Unrecognized equation system type");
const char *toofew = N_("An equation system must have at least two equations");

static int gretl_system_type_from_string (const char *str)
{
    int i = 0;

    while (gretl_system_type_strings[i] != NULL) {
	if (!strcmp(str, gretl_system_type_strings[i]))
	    return i;
	i++;
    }

    return -1;
}

static int on_exo_list (const int *exo_list, int k)
{
    int i;

    for (i=1; i<=exo_list[0]; i++) {
	if (exo_list[i] < 0) break;
	if (k == exo_list[i]) return 1;
    }

    return 0;
}

static int is_endogenous (const int *endo_list, int k)
{
    int i;

    for (i=1; i<=endo_list[0]; i++) {
	if (k == endo_list[i]) return 1;
    }

    return 0;
}

static int make_exo_list (gretl_equation_system *sys)
{
    int i, j, k, nexo, maxnexo = 0;

    for (i=0; i<sys->n_equations; i++) {
	for (j=2; j<=sys->lists[i][0]; j++) {
	    k = sys->lists[i][j];
	    if (!is_endogenous(sys->endog_vars, k)) {
		maxnexo++;
	    }
	}
    }

    sys->exog_vars = malloc((maxnexo + 1) * sizeof *sys->exog_vars);
    if (sys->exog_vars == NULL) {
	return 1;
    }

    sys->exog_vars[0] = maxnexo;
    for (i=1; i<=maxnexo; i++) {
	sys->exog_vars[i] = -1;
    }

    /* form list of exogenous variables */
    nexo = 0;
    for (i=0; i<sys->n_equations; i++) {
	for (j=2; j<=sys->lists[i][0]; j++) {
	    k = sys->lists[i][j];
	    if (!is_endogenous(sys->endog_vars, k) && 
		!on_exo_list(sys->exog_vars, k)) {
		nexo++;
		sys->exog_vars[nexo] = k;
	    }
	}
    } 

    if (nexo < maxnexo) {
	sys->exog_vars[0] = nexo;
    }

    return 0;
}

static int fiml_categorize_vars (gretl_equation_system *sys)
{
    int nendog = sys->n_equations;
    int i, j;

    sys->endog_vars = malloc((nendog + 1) * sizeof *sys->endog_vars);
    if (sys->endog_vars == NULL) {
	return 1;
    }

    sys->endog_vars[0] = nendog;
    for (i=0; i<sys->n_equations; i++) {
	sys->endog_vars[i+1] = sys->lists[i][1];
    }

    /* check for duplication */
    for (i=1; i<=nendog; i++) {
	for (j=i+1; j<=nendog; j++) {
	    if (sys->endog_vars[i] == sys->endog_vars[j]) {
		fprintf(stderr, "Duplicated endogenous variable, number %d\n", 
			sys->endog_vars[i]);
		return 1;
	    }
	}
    }

    if (make_exo_list(sys)) {
	return 1;
    }

#if 1
    fprintf(stderr, "FIML system\n"
	    " Number of endogenous vars = %d\n"
	    " Number of exogenous vars = %d\n",
	    sys->endog_vars[0], sys->exog_vars[0]);
    printlist(sys->endog_vars, "endog");
    printlist(sys->exog_vars, "exog");
#endif

    return 0;
}

static gretl_equation_system *gretl_equation_system_new (int type)
{
    gretl_equation_system *sys;

    if (type < 0) return NULL;

    sys = malloc(sizeof *sys);
    if (sys == NULL) return NULL;

    sys->type = type;
    sys->n_equations = 0;
    sys->flags = 0;
    sys->lists = NULL;
    sys->endog_vars = NULL;
    sys->exog_vars = NULL;

    return sys;
}

void gretl_equation_system_destroy (gretl_equation_system *sys)
{
    int i;

    if (sys == NULL || sys->lists == NULL) return;

    for (i=0; i<sys->n_equations; i++) {
	free(sys->lists[i]);
    }
    free(sys->lists);
    sys->lists = NULL;

    free(sys->endog_vars);
    free(sys->exog_vars);

    free(sys);
}

int gretl_equation_system_append (gretl_equation_system *sys, 
				  int *list)
{
    int i, neq;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    neq = sys->n_equations;

    sys->lists = realloc(sys->lists, (neq + 1) * sizeof *sys->lists);
    if (sys->lists == NULL) return E_ALLOC;

    sys->lists[neq] = malloc((list[0] + 1) * sizeof *list);
    if (sys->lists[neq] == NULL) {
	for (i=0; i<neq; i++) {
	    free(sys->lists[i]);
	}
	free(sys->lists);
	sys->lists = NULL;
	return E_ALLOC;
    }

    for (i=0; i<=list[0]; i++) {
	sys->lists[neq][i] = list[i];
    }

    if (sys->type == SUR) {
	rearrange_list(sys->lists[neq]);
    }

    sys->n_equations += 1;

    return 0;
}

gretl_equation_system *parse_system_start_line (const char *line)
{
    char sysstr[9];
    gretl_equation_system *sys = NULL;
    int systype = -1;

    if (sscanf(line, "system type=%8s\n", sysstr) == 1) {
	lower(sysstr);
	systype = gretl_system_type_from_string(sysstr);
    } 

    if (systype >= 0) {
	sys = gretl_equation_system_new(systype);
    } else {
	strcpy(gretl_errmsg, _(badsystem));
    }

    if (strstr(line, "save=")) {
	if (strstr(line, "resids") || strstr(line, "uhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_UHAT;
	}
	if (strstr(line, "fitted") || strstr(line, "yhat")) {
	    sys->flags |= GRETL_SYSTEM_SAVE_YHAT;
	}
    }

    return sys;
}

int gretl_equation_system_finalize (gretl_equation_system *sys, 
				    double ***pZ, DATAINFO *pdinfo,
				    PRN *prn)
{
    int err = 0;
    void *handle = NULL;
    int (*system_est) (gretl_equation_system *, 
		       double ***, DATAINFO *, PRN *);

    *gretl_errmsg = 0;

    if (sys == NULL) {
	strcpy(gretl_errmsg, _(nosystem));
	return 1;
    }

    if (sys->type != SUR && sys->type != THREESLS && sys->type != FIML) {
	err = 1;
	strcpy(gretl_errmsg, _(badsystem));
	goto system_bailout;
    }

    if (sys->n_equations < 2) {
	err = 1;
	strcpy(gretl_errmsg, _(toofew));
	goto system_bailout;
    }

    if (sys->type == FIML) {
	err = fiml_categorize_vars(sys);
	if (err) {
	    goto system_bailout; 
	}
	goto system_bailout; 
    }

    system_est = get_plugin_function("system_estimate", &handle);

    if (system_est == NULL) {
	err = 1;
        goto system_bailout;
    }

    pputc(prn, '\n');
    pprintf(prn, _("Equation system, %s\n\n"),
	    gretl_system_long_strings[sys->type]);

    err = (* system_est) (sys, pZ, pdinfo, prn);
    
 system_bailout:
    if (handle != NULL) {
	close_plugin(handle);
    }

    /* for now, we'll free the system after printing */
    gretl_equation_system_destroy(sys);

    return err;
}

static int get_real_list_length (const int *list)
{
    int i, len = list[0];

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    len = i - 1;
	    break;
	}
    }

    return len;
}

int system_max_indep_vars (const gretl_equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->n_equations; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	if (nvi > nv) nv = nvi;
    }

    return nv;
}

int system_n_indep_vars (const gretl_equation_system *sys)
{
    int i, nvi, nv = 0;

    for (i=0; i<sys->n_equations; i++) {
	nvi = get_real_list_length(sys->lists[i]) - 1;
	nv += nvi;
    }

    return nv;
}

const char *gretl_system_short_string (const MODEL *pmod)
{
    int i = gretl_model_get_int(pmod, "systype");

    return gretl_system_short_strings[i];
}

int system_adjust_t1t2 (const gretl_equation_system *sys,
			int *t1, int *t2, const double **Z)
{
    int i, misst, err = 0;

    for (i=0; i<sys->n_equations && !err; i++) {
	err = adjust_t1t2(NULL, sys->lists[i], t1, t2, Z, &misst);
    }

    return err;
}

/* simple accessor functions */

int system_save_uhat (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYSTEM_SAVE_UHAT;
}

int system_save_yhat (const gretl_equation_system *sys)
{
    return sys->flags & GRETL_SYSTEM_SAVE_YHAT;
}

int system_n_equations (const gretl_equation_system *sys)
{
    return sys->n_equations;
}

int *system_get_list (const gretl_equation_system *sys, int i)
{
    if (i >= sys->n_equations) return NULL;

    return sys->lists[i];
}

int system_get_depvar (const gretl_equation_system *sys, int i)
{
    if (i >= sys->n_equations) return 0;

    return sys->lists[i][1];
}

int system_get_type (const gretl_equation_system *sys)
{
    return sys->type;
}

int *system_get_endog_vars (const gretl_equation_system *sys)
{
    return sys->endog_vars;
}

int *system_get_exog_vars (const gretl_equation_system *sys)
{
    return sys->endog_vars;
}
