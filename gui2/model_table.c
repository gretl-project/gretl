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

/* model_table.c for gretl */

#include "gretl.h"
#include "model_table.h"

static const MODEL **model_list;
static int model_list_len;
static int *grand_list;

static void tex_print_model_table (void);

GtkItemFactoryEntry model_table_items[] = {
#ifdef USE_GNOME
    { N_("/_File"), NULL, NULL, 0, "<Branch>" },     
    { N_("/File/_Print..."), NULL, window_print, 0, NULL },
#endif
    { N_("/_Edit"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/_Copy selection"), NULL, text_copy, COPY_SELECTION, NULL },
    { N_("/Edit/Copy _all"), NULL, NULL, 0, "<Branch>" },
    { N_("/Edit/Copy all/as plain _text"), NULL, text_copy, COPY_TEXT, NULL },
    { N_("/Edit/Copy all/as _LaTeX"), NULL, tex_print_model_table, 0, NULL },
#if 0
    { N_("/Edit/Copy all/as _RTF"), NULL, text_copy, COPY_RTF, NULL },
#endif
    { NULL, NULL, NULL, 0, NULL }
};


static int model_already_listed (const MODEL *pmod)
{
    int i;

    for (i=0; i<model_list_len; i++) {
	if (pmod == model_list[i]) return 1;
    }

    return 0;
}

int start_model_list (const MODEL *pmod, int add_mode)
{
    model_list = mymalloc(sizeof *model_list);
    if (model_list == NULL) return 1;

    model_list_len = 1;
    model_list[0] = pmod;

    if (add_mode == MODEL_ADD_FROM_MENU) {
	infobox(_("Model added to table")); 
    }   

    return 0;
}

void remove_from_model_list (const MODEL *pmod)
{
    int i;

    if (model_list_len == 0 || model_list == NULL) 
	return;

    for (i=0; i<model_list_len; i++) {
	if (pmod == model_list[i]) {
	    model_list[i] = NULL;
	}
    }
}

int add_to_model_list (const MODEL *pmod, int add_mode)
{
    const MODEL **tmp;

    /* check that list is really started */
    if (model_list_len == 0) {
	return start_model_list(pmod, add_mode);
    }

    /* check that the dependent variable is in common */
    if (pmod->list[1] != (model_list[0])->list[1]) {
	errbox(_("Can't add model to table -- this model has a "
		 "different dependent variable"));
	return 1;
    }

    /* check that model is not already on the list */
    if (model_already_listed(pmod)) {
	errbox(_("Model is already included in the table"));
	return 0;
    }

    model_list_len++;
    tmp = myrealloc(model_list, model_list_len * sizeof *model_list);
    if (tmp == NULL) {
	free(model_list);
	return 1;
    }

    model_list = tmp;
    model_list[model_list_len - 1] = pmod;

    if (add_mode == MODEL_ADD_FROM_MENU) {
	infobox(_("Model added to table"));
    }

    return 0;
}

void free_model_list (void)
{
    free(model_list);
    model_list = NULL;
    free(grand_list);
    grand_list = NULL;
    model_list_len = 0;
#if 0
    infobox(_("Model table cleared"));
#endif
}

static int var_is_in_model (int v, const MODEL *pmod)
{
    int i;

    for (i=2; i<=pmod->list[0]; i++) {
	if (v == pmod->list[i]) return i;
    }

    return 0;    
}

static int on_grand_list (int v)
{
    int i;

    for (i=2; i<=grand_list[0]; i++) {
	if (v == grand_list[i]) return 1;
    }

    return 0;
}

static void add_to_grand_list (const int *list)
{
    int i, j = grand_list[0] + 1;

    for (i=2; i<=list[0]; i++) {
	if (!on_grand_list(list[i])) {
	    grand_list[0] += 1;
	    grand_list[j++] = list[i];
	}
    }
}

static int make_grand_varlist (void)
{
    int i, j;
    int l0 = 0;
    const MODEL *pmod;

    free(grand_list);

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] == NULL) continue;
	l0 += (model_list[i])->list[0];
    }

    grand_list = mymalloc((l0 + 1) * sizeof *grand_list);
    if (grand_list == NULL) return 1;

    for (i=0; i<model_list_len; i++) {
	pmod = model_list[i];
	if (pmod == NULL) continue;
	if (i == 0) {
	    for (j=0; j<=pmod->list[0]; j++) {
		grand_list[j] = pmod->list[j];
	    }
	} else {
	    add_to_grand_list(pmod->list);
	}
    }

    return 0;
}

static int model_list_empty (void)
{
    int i, real_n_models = 0;

    if (model_list_len == 0 || model_list == NULL) 
	return 1;

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] != NULL) 
	    real_n_models++;
    }

    return (real_n_models == 0);
}

static int common_estimator (void)
{
    int i, ci0;

    ci0 = (model_list[0])->ci;

    if (model_list_len == 1) return ci0;

    for (i=1; i<model_list_len; i++) {
	if ((model_list[i])->ci != ci0) return 0;
    }  

    return ci0;
}

static int common_df (void)
{
    int i, dfn0, dfd0;

    if (model_list_len == 1) return 1;

    dfn0 = (model_list[0])->dfn;
    dfd0 = (model_list[0])->dfd;

    for (i=1; i<model_list_len; i++) {
	if ((model_list[i])->dfn != dfn0) return 0;
	if ((model_list[i])->dfd != dfd0) return 0;
    }  

    return 1;
}

static void center_in_field (const char *s, int width, PRN *prn)
{
    int rem = width - strlen(s);

    if (rem <= 1) {
	pprintf(prn, "%s", s);
    }
    else {
	int i, off = rem / 2;

	for (i=0; i<off; i++) {
	    pputs(prn, " ");
	}
	pprintf(prn, "%-*s", width - off, s);
    }
}

static const char *short_estimator_string (int ci, int format)
{
    if (ci == HSK) return N_("HSK");
    else if (ci == CORC) return N_("CORC");
    else if (ci == HILU) return N_("HILU");
    else if (ci == ARCH) return N_("ARCH");
    else if (ci == POOLED) return N_("OLS");
    else return estimator_string (ci, format);
}

int display_model_table (void)
{
    int i, j, gl0, ci;
    int same_df;
    const MODEL *pmod;
    PRN *prn;
    char se[16];

    if (model_list_empty()) {
	errbox(_("The model table is empty"));
	return 1;
    }

    if (make_grand_varlist()) return 1;

    if (bufopen(&prn)) {
	free_model_list();
	return 1;
    }

    ci = common_estimator();

    gl0 = grand_list[0];

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, _("%s estimates"), 
		_(estimator_string(ci, prn->format)));
	pputs(prn, "\n");
    }

    pprintf(prn, _("Dependent variable: %s\n"),
	    datainfo->varname[grand_list[1]]);

    pputs(prn, _("Standard errors in parentheses\n\n"));

    pputs(prn, "            ");
    for (j=0; j<model_list_len; j++) {
	char modhd[16];

	if (model_list[j] == NULL) continue;
	sprintf(modhd, _("Model %d"), (model_list[j])->ID);
	center_in_field(modhd, 12, prn);
    }
    pputs(prn, "\n");
    
    if (ci == 0) {
	char est[12];	

	pputs(prn, "            ");
	for (j=0; j<model_list_len; j++) {
	    if (model_list[j] == NULL) continue;
	    strcpy(est, 
		   _(short_estimator_string((model_list[j])->ci,
					    prn->format)));
	    center_in_field(est, 12, prn);
	}
	pputs(prn, "\n");
    }

    pputs(prn, "\n");    

    /* print coefficients, standard errors */
    for (i=2; i<=gl0; i++) {
	int k, v = grand_list[i];

	pprintf(prn, "%8s ", datainfo->varname[v]);
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		pprintf(prn, "%#12.5g", pmod->coeff[k-1]);
	    } else {
		pputs(prn, "            ");
	    }
	}
	pputs(prn, "\n          ");
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		sprintf(se, "(%#.5g)", pmod->sderr[k-1]);
		pprintf(prn, "%12s", se);
	    } else {
		pputs(prn, "            ");
	    }
	}
	pputs(prn, "\n\n");
    }

    /* print sample sizes, R-squared */
    pprintf(prn, "%8s ", _("n"));
    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	pprintf(prn, "%12d", pmod->nobs);
    }
    pputs(prn, "\n");

    same_df = common_df();
    pprintf(prn, "%9s", (same_df)? _("R-squared") : _("Adj. R**2"));

    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    /* McFadden */
	    pprintf(prn, "%#12.4g", pmod->rsq);
	} else {
	    pprintf(prn, "%#12.4g", (same_df)? pmod->rsq : pmod->adjrsq);
	}
    }
    pputs(prn, "\n");

    view_buffer(prn, 78, 450, _("gretl: model table"), PRINT, 
		model_table_items);

    return 0;
}

static void tex_print_model_table (void)
{
    int i, j, gl0, ci;
    int same_df;
    const MODEL *pmod;
    PRN *prn;

    if (model_list_empty()) {
	errbox(_("The model table is empty"));
	return;
    }

    if (make_grand_varlist()) return;

    if (bufopen(&prn)) return;

    ci = common_estimator();

    gl0 = grand_list[0];

    pputs(prn, "\\begin{center}\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, I_("%s estimates"), 
		I_(estimator_string(ci, prn->format)));
	pputs(prn, "\\\\\n");
    }

    pprintf(prn, "%s: %s \\\\\n", I_("Dependent variable"),
	    datainfo->varname[grand_list[1]]);

    pputs(prn, I_("Standard errors in parentheses\n\n"));

    pputs(prn, "\\vspace{1em}\n\n");
    pputs(prn, "\\begin{tabular}{l");
    for (j=0; j<model_list_len; j++) {
	pputs(prn, "c");
    }
    pputs(prn, "}\n");

    for (j=0; j<model_list_len; j++) {
	char modhd[16];

	if (model_list[j] == NULL) continue;
	sprintf(modhd, I_("Model %d"), (model_list[j])->ID);
	pprintf(prn, " & %s ", modhd);
    }
    pputs(prn, "\\\\ ");
    
    if (ci == 0) {
	char est[12];

	pputs(prn, "\n");

	for (j=0; j<model_list_len; j++) {
	    if (model_list[j] == NULL) continue;
	    strcpy(est, 
		   I_(short_estimator_string((model_list[j])->ci,
					    prn->format)));
	    pprintf(prn, " & %s ", est);
	}
	pputs(prn, "\\\\ ");
    }

    pputs(prn, " [6pt] \n");    

    /* print coefficients, standard errors */
    for (i=2; i<=gl0; i++) {
	int k, v = grand_list[i];

	pprintf(prn, "%s ", datainfo->varname[v]);
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		double x = pmod->coeff[k-1];

		if (x < 0) {
		    pprintf(prn, "& $-$%#.5g ", fabs(x));
		} else {
		    pprintf(prn, "& %#.5g ", x);
		}
	    } else {
		pputs(prn, "& ");
	    }
	}
	pputs(prn, "\\\\\n");
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		pprintf(prn, "& (%#.5g) ", pmod->sderr[k-1]);
	    } else {
		pputs(prn, "& ");
	    }
	}
	pputs(prn, "\\\\ [4pt] \n");
    }

    /* print sample sizes, R-squared */
    pprintf(prn, "$%s$ ", _("n"));
    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	pprintf(prn, "& %d ", pmod->nobs);
    }
    pputs(prn, "\\\\\n");

    same_df = common_df();
    pputs(prn, (same_df)? "$R^2$" : "$\\bar R^2$ ");

    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	    /* McFadden */
	    pprintf(prn, "& %.4f ", pmod->rsq);
	} else {
	    pprintf(prn, "& %.4f ", (same_df)? pmod->rsq : pmod->adjrsq);
	}
    }
    pputs(prn, "\n");

    pputs(prn, "\\end{tabular}\n\\end{center}");

    prn_to_clipboard(prn, COPY_LATEX);
}
