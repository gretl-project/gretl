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

typedef struct id_atom_ id_atom;
typedef struct identity_ identity;

enum {
    OP_PLUS,
    OP_MINUS
} identity_ops;

enum {
    ENDOG_LIST,
    INSTR_LIST
} aux_list_types;

struct id_atom_ {
    int op;
    int varnum;
};

struct identity_ {
    int n_atoms;
    int depvar;
    id_atom *atoms;
};

struct _gretl_equation_system {
    int type;
    int n_equations;
    int n_identities;
    int n_obs;
    char flags;
    int **lists;
    int *endog_vars;
    int *instr_vars;
    identity **idents;
    gretl_matrix *uhat;
    MODEL **models;
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

static void destroy_ident (identity *pident);
static int make_instrument_list (gretl_equation_system *sys);

/* FIML checking stuff */

static void 
print_ident (const identity *pident, const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    pprintf(prn, "Identity: %s = %s ", 
	    pdinfo->varname[pident->depvar],
	    pdinfo->varname[pident->atoms[0].varnum]);
    for (i=1; i<pident->n_atoms; i++) {
	pprintf(prn, "%c %s ", (pident->atoms[i].op == OP_PLUS)? '+' : '-',
		pdinfo->varname[pident->atoms[i].varnum]);
    }
    pputc(prn, '\n');
}

void print_fiml_sys_info (const gretl_equation_system *sys, 
			  const DATAINFO *pdinfo, PRN *prn)
{
    int i;

    for (i=0; i<sys->n_identities; i++) {
	print_ident(sys->idents[i], pdinfo, prn);
    }

    if (sys->endog_vars != NULL) {
	pputs(prn, "Endogenous variables:");
	for (i=1; i<=sys->endog_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->endog_vars[i]]);
	}
	pputc(prn, '\n');
    }

    if (sys->instr_vars != NULL) {
	pputs(prn, "Instruments:");
	for (i=1; i<=sys->instr_vars[0]; i++) {
	    pprintf(prn, " %s", pdinfo->varname[sys->instr_vars[i]]);
	}
	pputc(prn, '\n');
    }

}

/* end FIML checking stuff */

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

static gretl_equation_system *gretl_equation_system_new (int type)
{
    gretl_equation_system *sys;

    if (type < 0) return NULL;

    sys = malloc(sizeof *sys);
    if (sys == NULL) return NULL;

    sys->type = type;
    sys->n_equations = 0;
    sys->n_identities = 0;
    sys->n_obs = 0;
    sys->flags = 0;
    sys->lists = NULL;
    sys->endog_vars = NULL;
    sys->instr_vars = NULL;
    sys->idents = NULL;
    sys->uhat = NULL;
    sys->models = NULL;

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

    for (i=0; i<sys->n_identities; i++) {
	destroy_ident(sys->idents[i]);
    }
    free(sys->idents);

    free(sys->endog_vars);
    free(sys->instr_vars);

    gretl_matrix_free(sys->uhat);

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

gretl_equation_system *system_start (const char *line)
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
	err = make_instrument_list(sys);
	if (err) goto system_bailout;
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

int *compose_tsls_list (const gretl_equation_system *sys, int i)
{
    int *list;
    int j, k1, k2;

    if (i >= sys->n_equations || sys->instr_vars == NULL) {
	return NULL;
    }

    k1 = sys->lists[i][0];
    k2 = sys->instr_vars[0];

    list = malloc((k1 + k2 + 2) * sizeof *list);
    if (list == NULL) return NULL;

    list[0] = k1 + k2 + 1;
    for (j=1; j<=list[0]; j++) {
	if (j <= k1) {
	    list[j] = sys->lists[i][j];
	} else if (j == k1 + 1) {
	    list[j] = LISTSEP;
	} else {
	    list[j] = sys->instr_vars[j - (k1 + 1)];
	}
    }

    return list;
}

void system_set_n_obs (gretl_equation_system *sys, int n)
{
    sys->n_obs = n;
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

int system_n_indentities (const gretl_equation_system *sys)
{
    return sys->n_identities;
}

int system_n_obs (const gretl_equation_system *sys)
{
    return sys->n_obs;
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

int *system_get_instr_vars (const gretl_equation_system *sys)
{
    return sys->instr_vars;
}

void system_attach_uhat (gretl_equation_system *sys, gretl_matrix *u)
{
    if (sys->uhat != NULL) {
	gretl_matrix_free(sys->uhat);
    }
    sys->uhat = u;
}

const gretl_matrix *system_get_uhat (const gretl_equation_system *sys)
{
    return sys->uhat;
}

void system_attach_models (gretl_equation_system *sys, MODEL **models)
{
    sys->models = models;
}

MODEL *system_get_model (const gretl_equation_system *sys, int i)
{
    if (sys->models == NULL || i >= sys->n_equations) {
	return NULL; 
    }   

    return sys->models[i];
}

/* dealing with identities (FIML) */

int rhs_var_in_identity (const gretl_equation_system *sys, int lhsvar,
			 int rhsvar)
{
    const identity *ident;
    int i, j;

    for (i=0; i<sys->n_identities; i++) {
	ident = sys->idents[i];
	if (ident->depvar == lhsvar) {
	    for (j=0; j<ident->n_atoms; j++) {
		if (ident->atoms[j].varnum == rhsvar) {
		    return (ident->atoms[j].op == OP_PLUS)? 1 : -1;
		}
	    }
	}
    }

    return 0;
}

int eval_identity (double *targ, identity *ident,
		   const double **Z, int t1, int t2)
{
    int i, k, t;

    for (t=t1; t<=t2; t++) {
	targ[t] = 0.0;
	for (i=0; i<ident->n_atoms; i++) {
	    k = ident->atoms[i].varnum;
	    if (ident->atoms[i].op == OP_PLUS) {
		targ[t] += Z[k][t];
	    } else {
		targ[t] -= Z[k][t];
	    }
	}
    }

    return 0;
}

static void destroy_ident (identity *pident)
{
    free(pident->atoms);
    free(pident);
}

static identity *ident_new (int nv)
{
    identity *pident;

    pident = malloc(sizeof *pident);
    if (pident == NULL) return NULL;

    pident->n_atoms = nv;
    pident->atoms = malloc(nv * sizeof *pident->atoms);
    if (pident->atoms == NULL) {
	free(pident);
	pident = NULL;
    }

    return pident;
}

static identity *
parse_identity (const char *str, const DATAINFO *pdinfo, int *err)
{
    identity *pident;
    const char *p;
    char f1[24], f2[16];
    char op, vname1[VNAMELEN], vname2[VNAMELEN];
    int i, nv;

    sprintf(f1, "%%%ds = %%%d[^+ -]", VNAMELEN - 1, VNAMELEN - 1);
    sprintf(f2, "%%c %%%d[^+ -]", VNAMELEN - 1);

    if (sscanf(str, f1, vname1, vname2) != 2) {
	*err = E_PARSE;
	return NULL;
    }

    p = str;
    nv = 1;
    while (*p) {
	if (*p == '+' || *p == '-') nv++;
	p++;
    }

    pident = ident_new(nv);
    if (pident == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    pident->depvar = varindex(pdinfo, vname1);
    if (pident->depvar == pdinfo->v) {
	destroy_ident(pident);
	*err = E_UNKVAR;
	return NULL;
    }

    pident->atoms[0].op = OP_PLUS;
    pident->atoms[0].varnum = varindex(pdinfo, vname2);
    if (pident->atoms[0].varnum == pdinfo->v) {
	destroy_ident(pident);
	*err = E_UNKVAR;
	return NULL;
    }

    p = str;
    for (i=1; i<nv && !*err; i++) {
	p += strcspn(p, "+-");
	sscanf(p, f2, &op, vname1);
	if (op == '+') op = OP_PLUS;
	else if (op == '-') op = OP_MINUS;
	else *err = E_PARSE;
	if (!*err) {
	    pident->atoms[i].op = op;
	    pident->atoms[i].varnum = varindex(pdinfo, vname1);
	    if (pident->atoms[i].varnum == pdinfo->v) {
		*err = E_UNKVAR;
	    }
	}
	p++;
    }

    if (*err) {
	destroy_ident(pident);
	pident = NULL;
    }
       
    return pident;
}

static int 
add_identity_to_sys (gretl_equation_system *sys, const char *line,
		     const DATAINFO *pdinfo)
{
    identity **ppident;
    identity *pident;
    int ni = sys->n_identities;
    int err = 0;

    pident = parse_identity(line, pdinfo, &err);
    if (pident == NULL) return err;

    /* connect the identity to the equation system */
    ppident = realloc(sys->idents, (ni + 1) * sizeof *sys->idents);
    if (ppident == NULL) {
	destroy_ident(pident);
	return E_ALLOC;
    }

    sys->idents = ppident;
    sys->idents[ni] = pident;
    sys->n_identities += 1;

    return 0;
}

static int
add_aux_list_to_sys (gretl_equation_system *sys, const char *line,
		     const DATAINFO *pdinfo, int which)
{
    const char *p;
    char vname[VNAMELEN];
    int *list;
    int i, v, nf, len, cplen;
    int err = 0;

    if (which == ENDOG_LIST) {
	if (sys->endog_vars != NULL) {
	    strcpy(gretl_errmsg, "Only one list of endogenous variables may be given");
	    return 1;
	}
    } else if (which == INSTR_LIST) {
	if (sys->instr_vars != NULL) {
	    strcpy(gretl_errmsg, "Only one list of instruments may be given");
	    return 1;
	}
	if (sys->type != THREESLS) {
	    strcpy(gretl_errmsg, "Instruments may only be specified for 3SLS");
	    return 1;
	}
    } else {
	return 1;
    }

    nf = count_fields(line);
    if (nf < 1) return 1;

    list = malloc((nf + 1) * sizeof *list);
    if (list == NULL) return E_ALLOC;

    list[0] = nf;
    
    p = line;
    for (i=1; i<=nf && !err; i++) {
	while (isspace(*p)) p++;
	*vname = '\0';
	cplen = len = strcspn(p, " \t\n");
	if (cplen > VNAMELEN - 1) {
	    cplen = VNAMELEN - 1;
	}
	strncat(vname, p, cplen);
	if (isdigit(*vname)) {
	    v = atoi(vname);
	} else {
	    v = varindex(pdinfo, vname);
	}
	if (v < 0 || v >= pdinfo->v) {
	    sprintf(gretl_errmsg, "Undefined variable '%s'.", vname);
	    err = 1;
	} else {
	    list[i] = v;
	}
	p += len;
    }
	
    if (err) {
	free(list);
	return err;
    }

    if (which == ENDOG_LIST) {
	sys->endog_vars = list;
    } else {
	sys->instr_vars = list;
    }

    return 0;
}

int 
system_parse_line (gretl_equation_system *sys, const char *line,
		   const DATAINFO *pdinfo)
{
    int err = 1;

    *gretl_errmsg = '\0';

    if (strncmp(line, "identity", 8) == 0) {
	err = add_identity_to_sys(sys, line + 8, pdinfo);
    } 
    else if (strncmp(line, "endog", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, ENDOG_LIST);
    }
    else if (strncmp(line, "instr", 5) == 0) {
	err = add_aux_list_to_sys(sys, line + 5, pdinfo, INSTR_LIST);
    }

    if (err) {
	gretl_equation_system_destroy(sys);
    }

    return err;
}

/* More FIML-related functionality */

static int in_list (const int *list, int k)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] < 0) break;
	if (k == list[i]) return 1;
    }

    return 0;
}

static int make_instrument_list (gretl_equation_system *sys)
{
    int *ilist, *elist = sys->endog_vars;
    int i, j, k, nexo, maxnexo = 0;

    if (elist == NULL) return 1;

    /* First pass: get a count of the max possible number of
       exogenous variables (probably an over-estimate due to
       double-counting).
    */

    for (i=0; i<sys->n_equations; i++) {
	const int *slist = sys->lists[i];

	for (j=2; j<=slist[0]; j++) {
	    if (!in_list(elist, slist[j])) maxnexo++;
	}
    }

    for (i=0; i<sys->n_identities; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    if (!in_list(elist, ident->atoms[j].varnum)) maxnexo++;
	}
    }

    ilist = malloc((maxnexo + 1) * sizeof *ilist);
    if (ilist == NULL) {
	return E_ALLOC;
    }

    ilist[0] = maxnexo;
    for (i=1; i<=maxnexo; i++) {
	ilist[i] = -1;
    }

    /* Form list of exogenous variables, drawing on both the
       stochastic equations and the identities, if any. 
    */

    nexo = 0;

    for (i=0; i<sys->n_equations; i++) {
	const int *slist = sys->lists[i];

	for (j=2; j<=slist[0]; j++) {
	    k = slist[j];
	    if (!in_list(elist, k) && 
		!in_list(ilist, k)) {
		ilist[++nexo] = k;
	    }
	}
    } 

    for (i=0; i<sys->n_identities; i++) {
	const identity *ident = sys->idents[i];

	for (j=0; j<ident->n_atoms; j++) {
	    k = ident->atoms[j].varnum;
	    if (!in_list(elist, k) &&
		!in_list(ilist, k)) {
		ilist[++nexo] = k;
	    }
	}
    }

    /* Trim the list (remove effect of double-counting) */

    if (nexo < maxnexo) {
	int *plist = realloc(ilist, (nexo + 1) * sizeof *plist);

	if (plist != NULL) {
	    ilist = plist;
	}
	ilist[0] = nexo;
    }

    sys->instr_vars = ilist;

    return 0;
}
