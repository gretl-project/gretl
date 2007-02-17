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
#include "session.h"
#include "textutil.h"
#include "texprint.h"

static MODEL **table_models;
static int n_models;
static int *grand_list;
static int use_tstats;

static void print_rtf_row_spec (PRN *prn, int tall);

#define MAX_TABLE_MODELS 6

static void mtable_errmsg (char *msg, int gui)
{
    if (gui) {
	errbox(msg);
    } else {
	gretl_errmsg_set(msg);
    }
}

static int real_table_n_models (void)
{
    int i, len = 0;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    len++;
	}
    }

    return len;    
}

static int model_table_too_many (int gui)
{
    if (real_table_n_models() == MAX_TABLE_MODELS) {
	mtable_errmsg(_("Model table is full"), gui);
	return 1;
    }

    return 0;
}

int in_model_table (const MODEL *pmod)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (pmod == table_models[i]) {
	    return 1;
	}
    }

    return 0;
}

int model_table_n_models (void)
{
    return n_models;
}

MODEL *model_table_model_by_index (int i)
{
    if (i >= 0 && i < n_models) {
	return table_models[i];
    } else {
	return NULL;
    }
}

void clear_model_table (PRN *prn)
{
    int i;

    for (i=0; i<n_models; i++) {
	/* reduce refcount on the model pointer */
	if (table_models[i] != NULL) {
	    gretl_object_unref(table_models[i], GRETL_OBJ_EQN);
	}
    }

    free(table_models);
    table_models = NULL;

    free(grand_list);
    grand_list = NULL;
    n_models = 0;
    
    if (prn != NULL) {
	pputs(prn, _("Model table cleared"));
	pputc(prn, '\n');
    }
}

static int model_table_depvar (void)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    return table_models[i]->list[1];
	}
    }

    return -1;
}

int add_to_model_table (MODEL *pmod, int add_mode, PRN *prn)
{
    int gui = (add_mode != MODEL_ADD_BY_CMD);

    if (pmod == NULL) {
	return 1;
    }

    /* NLS, MLE and GMM models won't work */
    if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM ||
	pmod->ci == ARBOND) {
	mtable_errmsg(_("Sorry, this model can't be put in the model table"),
		      gui);
	return 1;
    }

    /* nor will ARMA */
    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	mtable_errmsg(_("Sorry, ARMA models can't be put in the model table"),
		      gui);
	return 1;
    }

    /* nor TSLS */
    if (pmod->ci == TSLS) {
	mtable_errmsg(_("Sorry, TSLS models can't be put in the model table"),
		      gui);
	return 1;
    }    

    /* is the list started or not? */
    if (n_models == 0) {
	table_models = mymalloc(sizeof *table_models);
	if (table_models == NULL) {
	    return 1;
	}
	n_models = 1;
    } else {
	int dv = model_table_depvar();
	MODEL **mods;

	/* check that the dependent variable is in common */
	if (pmod->list[1] != dv) {
	    mtable_errmsg(_("Can't add model to table -- this model has a "
			    "different dependent variable"), gui);
	    return 1;
	}

	/* check that model is not already on the list */
	if (in_model_table(pmod)) {
	    mtable_errmsg(_("Model is already included in the table"), 0);
	    return 1;
	}

	/* check that the model table is not already full */
	if (model_table_too_many(gui)) {
	    return 1;
	}

	n_models++;
	mods = myrealloc(table_models, n_models * sizeof *mods);
	if (mods == NULL) {
	    clear_model_table(NULL);
	    return 1;
	}

	table_models = mods;
    }

    table_models[n_models - 1] = pmod;

    /* augment refcount so model won't get deleted */
    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    if (add_mode == MODEL_ADD_FROM_MENU) {
	infobox(_("Model added to table"));
    } else if (add_mode == MODEL_ADD_BY_CMD) {
	pputs(prn, _("Model added to table"));
	pputc(prn, '\n');
    }

    return 0;
}

static int var_is_in_model (int v, const MODEL *pmod)
{
    int i;

    for (i=2; i<=pmod->list[0]; i++) {
	if (pmod->list[i] == LISTSEP) {
	    break;
	}
	if (v == pmod->list[i]) {
	    return i;
	}
    }

    return 0;    
}

static int on_grand_list (int v)
{
    int i;

    for (i=2; i<=grand_list[0]; i++) {
	if (v == grand_list[i]) {
	    return 1;
	}
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

static int real_list_length (int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    return i - 1;
	}
    }

    return list[0];
}

static int make_grand_varlist (void)
{
    int i, j, first = 1;
    int l0 = 0;
    const MODEL *pmod;

    free(grand_list);

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    l0 += real_list_length(table_models[i]->list);
	}
    }

    grand_list = gretl_list_new(l0);
    if (grand_list == NULL) {
	return 1;
    }

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    pmod = table_models[i];
	    if (first) {
		for (j=0; j<=pmod->list[0]; j++) {
		    if (pmod->list[j] == LISTSEP) {
			break;
		    }
		    grand_list[j] = pmod->list[j];
		}
		first = 0;
	    } else {
		add_to_grand_list(pmod->list);
	    }
	}
    }

    return 0;
}

static int model_table_is_empty (void)
{
    int i, n = 0;

    if (n_models == 0 || table_models == NULL) { 
	return 1;
    }

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    n++;
	}
    }

    return (n == 0);
}

static int common_estimator (void)
{
    int i, ci0 = -1;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    if (ci0 == -1) {
		ci0 = table_models[i]->ci;
	    } else if (table_models[i]->ci != ci0) {
		return 0;
	    }
	}
    }  

    return ci0;
}

static int common_df (void)
{
    int i, dfn0 = -1, dfd0 = -1;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    if (dfn0 == -1) {
		dfn0 = table_models[i]->dfn;
		dfd0 = table_models[i]->dfd;
	    } else {
		if (table_models[i]->dfn != dfn0 ||
		    table_models[i]->dfd != dfd0) {
		    return 0;
		}
	    }
	}
    }  

    return 1;
}

static const char *short_estimator_string (int ci, PRN *prn)
{
    if (ci == HSK) return N_("HSK");
    else if (ci == CORC) return N_("CORC");
    else if (ci == HILU) return N_("HILU");
    else if (ci == PWE) return N_("PWE");
    else if (ci == ARCH) return N_("ARCH");
    else return estimator_string(ci, prn);
}

static const char *get_asts (double pval)
{
    return (pval >= 0.1)? "  " : (pval >= 0.05)? "* " : "**";
}

static const char *tex_get_asts (double pval)
{
    return (pval >= 0.1)? "" : (pval >= 0.05)? "$^{*}$" : "$^{**}$";
}

static const char *get_pre_asts (double pval)
{
    return (pval >= 0.1)? "" : (pval >= 0.05)? "$\\,$" : "$\\,\\,$";
}

static void print_model_table_coeffs (int nwidth, PRN *prn)
{
    int i, j, k;
    const MODEL *pmod;
    char tmp[32];
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);

    /* loop across all variables that appear in any model */
    for (i=2; i<=grand_list[0]; i++) {
	int v = grand_list[i];
	int f = 1;

	if (tex) {
	    tex_escape(tmp, datainfo->varname[v]);
	    pprintf(prn, "%s ", tmp);
	} else if (rtf) {
	    print_rtf_row_spec(prn, 0);
	    pprintf(prn, "\\intbl \\qc %s\\cell ", datainfo->varname[v]);
	} else {
	    pprintf(prn, "%-*s ", nwidth, datainfo->varname[v]);
	}

	/* print the coefficient estimates across a row */
	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) {
		continue;
	    }
	    if ((k = var_is_in_model(v, pmod))) {
		double x = screen_zero(pmod->coeff[k-2]);
		double s = screen_zero(pmod->sderr[k-2]);
		double pval;
		char numstr[32];

		if (floateq(s, 0.0)) {
		    if (floateq(x, 0.0)) {
			pval = 1.0;
		    } else {
			pval = 0.0001;
		    }
		} else {
		    pval = coeff_pval(pmod->ci, x / s, pmod->dfd);
		}

		sprintf(numstr, "%#.4g", x);
		gretl_fix_exponent(numstr);

		if (tex) {
		    if (x < 0) {
			pprintf(prn, "& %s$-$%s%s ", get_pre_asts(pval),
				numstr + 1, tex_get_asts(pval));
		    } else {
			pprintf(prn, "& %s%s%s ", get_pre_asts(pval), 
				numstr, tex_get_asts(pval));
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %s%s\\cell ", numstr, get_asts(pval));
		} else {
		    /* note: strlen(asts) = 2 */
		    pprintf(prn, "%*s%s", (f == 1)? 12 : 10,
			    numstr, get_asts(pval));
		}
		f = 0;
	    } else {
		/* variable not present in this column */
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    pputs(prn, "            "); /* 12 */
		}
	    }
	}

	/* terminate the coefficient row and start the next one,
	   which holds standard errors */
	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (rtf) {
	    pputs(prn, "\\intbl \\row\n");
	    print_rtf_row_spec(prn, 1);
	    pputs(prn, "\\intbl ");
	} else {
	    pputc(prn, '\n');
	    bufspace(nwidth + 2, prn);
	}

	/* print the t-stats or standard errors across a row */
	f = 1;
	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) {
		continue;
	    }
	    if ((k = var_is_in_model(v, pmod))) {
		char numstr[32];
		double val;

		if (use_tstats) {
		    val = pmod->coeff[k-2] / pmod->sderr[k-2];
		} else {
		    val = pmod->sderr[k-2];
		}

		sprintf(numstr, "%#.4g", val);
		gretl_fix_exponent(numstr);

		if (tex) {
		    if (val < 0) {
			pprintf(prn, "& \\footnotesize{($-$%s)} ", numstr + 1);
		    } else {
			pprintf(prn, "& \\footnotesize{(%s)} ", numstr);
		    }
		} else if (rtf) {
		    if (f == 1) {
			pputs(prn, "\\qc \\cell ");
		    }
		    pprintf(prn, "\\qc (%s)\\cell ", numstr);
		    f = 0;
		} else {
		    sprintf(tmp, "(%s)", numstr);
		    pprintf(prn, "%12s", tmp);
		}
	    } else {
		/* variable not present in this column */
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    pputs(prn, "            "); /* 12 */
		}
	    }
	}

	if (tex) {
	    pputs(prn, "\\\\ [4pt] \n");
	} else if (rtf) {
	    pputs(prn, "\\intbl \\row\n");
	} else {
	    pputs(prn, "\n\n");
	}
    }
}

static int any_log_lik (void)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (table_models[i] == NULL) {
	    continue;
	}
	if (!na(table_models[i]->lnL)) {
	    return 1;
	}
    }

    return 0;
}

static int any_r_squared (void)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (table_models[i] == NULL) {
	    continue;
	}
	if (!na(table_models[i]->rsq)) {
	    return 1;
	}
    }

    return 0;
}

static void print_n_r_squared (int wid, PRN *prn, int *binary)
{
    int same_df, any_R2, any_ll;
    const MODEL *pmod;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int j;

    if (rtf) {
	print_rtf_row_spec(prn, 0);
    }

    if (tex) {
	pprintf(prn, "$%s$ ", _("n"));
    } else if (rtf) {
	pprintf(prn, "\\intbl \\qc %s\\cell ", _("n"));
    } else {
	pprintf(prn, "%*s", wid, _("n"));
    }

    for (j=0; j<n_models; j++) {
	pmod = table_models[j];
	if (pmod != NULL) {
	    if (tex) {
		pprintf(prn, "& %d ", pmod->nobs);
	    } else if (rtf) {
		pprintf(prn, "\\qc %d\\cell ", pmod->nobs);
	    } else {
		pprintf(prn, "%12d", pmod->nobs);
	    }
	}
    }

    if (tex) {
	pputs(prn, "\\\\\n");
    } else if (rtf) {
	pputs(prn, "\\intbl \\row\n\\intbl ");
    } else {
	pputc(prn, '\n');
    }

    same_df = common_df();
    any_R2 = any_r_squared();
    any_ll = any_log_lik();

    if (any_R2) {
	/* print R^2 values */
	if (tex) {
	    pputs(prn, (same_df)? "$R^2$" : "$\\bar R^2$ ");
	} else if (rtf) {
	    pprintf(prn, "\\qc %s\\cell ", 
		    (same_df)? "R{\\super 2}" : _("Adj. R{\\super 2}"));
	} else {
	    pprintf(prn, "%*s", wid, (same_df)? _("R-squared") : _("Adj. R**2"));
	}

	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) continue;
	    if (na(pmod->rsq)) {
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    pputs(prn, "            ");
		}		
	    } else if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
		*binary = 1;
		/* McFadden */
		if (tex) {
		    pprintf(prn, "& %.4f ", pmod->rsq);
		} else if (rtf) {
		    pprintf(prn, "\\qc %.4f\\cell ", pmod->rsq);
		} else {
		    pprintf(prn, "%#12.4g", pmod->rsq);
		}
	    } else {
		double rsq = (same_df)? pmod->rsq : pmod->adjrsq;

		if (tex) {
		    pprintf(prn, "& %.4f ", rsq);
		} else if (rtf) {
		    pprintf(prn, "\\qc %.4f\\cell ", rsq);
		} else {
		    pprintf(prn, "%#12.4g", rsq);
		}
	    }
	}

	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (rtf) {
	    pputs(prn, "\\intbl \\row\n");
	} else {
	    pputc(prn, '\n');
	    if (!any_ll) {
		pputc(prn, '\n');
	    }
	}
    }

    if (any_ll) {
	/* print log-likelihoods */

	if (tex) {
	    pputs(prn, "$\\ell$");
	} else if (rtf) {
	    pputs(prn, "\\qc lnL\\cell ");
	} else {
	    pprintf(prn, "%*s", wid, "lnL");
	}

	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) continue;
	    if (na(pmod->lnL)) {
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    pputs(prn, "            ");
		}		
	    } else {
		if (tex) {
		    if (pmod->lnL > 0) {
			pprintf(prn, "& %.2f ", pmod->lnL);
		    } else {
			pprintf(prn, "& $-$%.2f ", -pmod->lnL);
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %.3f\\cell ", pmod->lnL);
		} else {
		    pprintf(prn, "%#12.6g", pmod->lnL);
		}
	    }
	}

	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (rtf) {
	    pputs(prn, "\\intbl \\row\n");
	} else {
	    pputs(prn, "\n\n");
	}
    }
}

static int grand_list_namelen (void)
{
    int i, len, maxlen = 8;

    for (i=2; i<=grand_list[0]; i++) {
	len = strlen(datainfo->varname[grand_list[i]]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    return maxlen;
}

int display_model_table (int gui)
{
    int j, ci;
    int binary = 0;
    int winwidth = 78;
    int namelen;
    PRN *prn;

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), gui);
	return 1;
    }

    if (make_grand_varlist()) {
	return 1;
    }

    namelen = grand_list_namelen();

    if (bufopen(&prn)) {
	clear_model_table(NULL);
	return 1;
    }

    ci = common_estimator();

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, _("%s estimates"), 
		_(estimator_string(ci, prn)));
	pputc(prn, '\n');
    }

    pprintf(prn, _("Dependent variable: %s\n"),
	    datainfo->varname[grand_list[1]]);

    pputc(prn, '\n');
    bufspace(namelen + 4, prn);

    for (j=0; j<n_models; j++) {
	char modhd[32];

	if (table_models[j] == NULL) {
	    continue;
	}
	if (table_models[j]->name != NULL) {
	    *modhd = '\0';
	    strncat(modhd, table_models[j]->name, 31);
	} else {
	    sprintf(modhd, _("Model %d"), table_models[j]->ID);
	}
	print_centered(modhd, 12, prn);
    }
    pputc(prn, '\n');
    
    if (ci == 0) {
	char est[32];	

	bufspace(namelen + 4, prn);
	for (j=0; j<n_models; j++) {
	    if (table_models[j] != NULL) {
		strcpy(est, 
		       _(short_estimator_string(table_models[j]->ci,
						prn)));
		print_centered(est, 12, prn);
	    }
	}
	pputc(prn, '\n');
    }

    pputc(prn, '\n'); 

    print_model_table_coeffs(namelen, prn);
    print_n_r_squared(namelen + 1, prn, &binary);

    if (use_tstats) {
	pprintf(prn, "%s\n", _("t-statistics in parentheses"));
    } else {
	pprintf(prn, "%s\n", _("Standard errors in parentheses"));
    }

    pprintf(prn, "%s\n", _("* indicates significance at the 10 percent level"));
    pprintf(prn, "%s\n", _("** indicates significance at the 5 percent level"));
   
    if (binary) {
	pprintf(prn, "%s\n", _("For logit and probit, R-squared is "
			       "McFadden's pseudo-R-squared"));
    }

    if (real_table_n_models() > 5) {
	winwidth = 90;
    }

    view_buffer(prn, winwidth, 450, _("gretl: model table"), VIEW_MODELTABLE, 
		NULL);

    return 0;
}

static int tex_print_model_table (PRN *prn)
{
    int j, ci;
    int binary = 0;
    char tmp[32];

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), 1);
	return 1;
    }

    if (make_grand_varlist()) {
	return 1;
    }

    gretl_print_set_format(prn, GRETL_FORMAT_TEX);

    ci = common_estimator();

    pputs(prn, "\\begin{center}\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, I_("%s estimates"), 
		I_(estimator_string(ci, prn)));
	pputs(prn, "\\\\\n");
    }

    tex_escape(tmp, datainfo->varname[grand_list[1]]);
    pprintf(prn, "%s: %s \\\\\n", I_("Dependent variable"), tmp);

    pputs(prn, "\\vspace{1em}\n\n");
    pputs(prn, "\\begin{tabular}{l");
    for (j=0; j<n_models; j++) {
	pputs(prn, "c");
    }
    pputs(prn, "}\n");

    for (j=0; j<n_models; j++) {
	char modhd[48];

	if (table_models[j] == NULL) {
	    continue;
	}
	if (table_models[j]->name != NULL) {
	    *modhd = '\0';
	    strncat(modhd, table_models[j]->name, 16);
	    tex_escape(tmp, modhd);
	    strcpy(modhd, tmp);
	} else {
	    sprintf(modhd, I_("Model %d"), table_models[j]->ID);
	}
	pprintf(prn, " & %s ", modhd);
    }
    pputs(prn, "\\\\ ");
    
    if (ci == 0) {
	char est[32];

	pputc(prn, '\n');

	for (j=0; j<n_models; j++) {
	    if (table_models[j] == NULL) {
		continue;
	    }
	    strcpy(est, 
		   I_(short_estimator_string(table_models[j]->ci,
					     prn)));
	    pprintf(prn, " & %s ", est);
	}
	pputs(prn, "\\\\ ");
    }

    pputs(prn, " [6pt] \n");   

    print_model_table_coeffs(0, prn);
    print_n_r_squared(0, prn, &binary);

    pputs(prn, "\\end{tabular}\n\n");
    pputs(prn, "\\vspace{1em}\n");

    if (use_tstats) {
	pprintf(prn, "%s\\\\\n", I_("$t$-statistics in parentheses"));
    } else {
	pprintf(prn, "%s\\\\\n", I_("Standard errors in parentheses"));
    }

    pprintf(prn, "{}%s\\\\\n", 
	    I_("* indicates significance at the 10 percent level"));
    pprintf(prn, "{}%s\\\\\n", 
	    I_("** indicates significance at the 5 percent level"));

    if (binary) {
	pprintf(prn, "%s\\\\\n", I_("For logit and probit, $R^2$ is "
				    "McFadden's pseudo-$R^2$"));
    }

    pputs(prn, "\\end{center}\n");

    return 0;
}

static void print_rtf_row_spec (PRN *prn, int tall)
{
    int i, cols = 1 + real_table_n_models();
    int col1 = 1000;
    int ht = (tall)? 362 : 262;

    pprintf(prn, "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh%d", ht);
    for (i=0; i<cols; i++) {
	pprintf(prn, "\\cellx%d", col1 +  i * 1400);
    }
    pputc(prn, '\n');
}

static int rtf_print_model_table (PRN *prn)
{
    int j, ci;
    int binary = 0;

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), 1);
	return 1;
    }

    if (make_grand_varlist()) return 1;

    gretl_print_set_format(prn, GRETL_FORMAT_RTF);

    ci = common_estimator();

    pputs(prn, "{\\rtf1\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pputs(prn, "\\par \\qc ");
	pprintf(prn, I_("%s estimates"), 
		I_(estimator_string(ci, prn)));
	pputc(prn, '\n');
    }

    pprintf(prn, "\\par \\qc %s: %s\n\\par\n\\par\n{", 
	    I_("Dependent variable"),
	    datainfo->varname[grand_list[1]]);

    /* RTF row stuff */
    print_rtf_row_spec(prn, 1);

    pputs(prn, "\\intbl \\qc \\cell ");
    for (j=0; j<n_models; j++) {
	char modhd[32];

	if (table_models[j] == NULL) {
	    continue;
	}
	if (table_models[j]->name != NULL) {
	    *modhd = '\0';
	    strncat(modhd, table_models[j]->name, 31);
	} else {
	    sprintf(modhd, I_("Model %d"), table_models[j]->ID);
	}
	pprintf(prn, "\\qc %s\\cell ", modhd);
    }
    pputs(prn, "\\intbl \\row\n");
    
    if (ci == 0) {
	char est[32];

	pputs(prn, "\\intbl \\qc \\cell ");

	for (j=0; j<n_models; j++) {
	    if (table_models[j] == NULL) continue;
	    strcpy(est, 
		   I_(short_estimator_string(table_models[j]->ci, prn)));
	    pprintf(prn, "\\qc %s\\cell ", est);
	}
	pputs(prn, "\\intbl \\row\n");
    }

    print_model_table_coeffs(0, prn);
    print_n_r_squared(0, prn, &binary);

    pputs(prn, "}\n\n");

    pprintf(prn, "\\par \\qc %s\n", I_("Standard errors in parentheses"));
    pprintf(prn, "\\par \\qc %s\n", 
	    I_("* indicates significance at the 10 percent level"));
    pprintf(prn, "\\par \\qc %s\n", 
	    I_("** indicates significance at the 5 percent level"));

    if (binary) {
	pprintf(prn, "\\par \\qc %s\n", I_("For logit and probit, "
					   "R{\\super 2} is "
					   "McFadden's pseudo-R{\\super 2}"));
    }

    pputs(prn, "\\par\n}\n");

    return 0;
}

int special_print_model_table (PRN *prn)
{
    if (tex_format(prn)) {
	return tex_print_model_table(prn);
    } else if (rtf_format(prn)) {
	return rtf_print_model_table(prn);
    } else {
	return 1;
    }
}

int modeltab_parse_line (const char *line, MODEL *pmod, PRN *prn)
{
    char cmdword[9];
    int err = 0;

    if (sscanf(line, "%*s %8s", cmdword) != 1) {
	return E_PARSE;
    }

    if (!strcmp(cmdword, "add")) {
	if (pmod == NULL || pmod->ID == 0) {
	    gretl_errmsg_set(_("No model is available"));
	    err = 1;
	} else {
	    err = add_to_model_table(pmod, MODEL_ADD_BY_CMD, prn);
	}
    } else if (!strcmp(cmdword, "show")) {
	err = display_model_table(0);
    } else if (!strcmp(cmdword, "free")) {
	if (model_table_is_empty()) {
	    mtable_errmsg(_("The model table is empty"), 0);
	    err = 1;
	} else {
	    clear_model_table(prn);
	}
    }

    return err;
}

void model_table_dialog (void)
{
    const char *opts[] = {
	N_("standard errors in parentheses"),
	N_("t-statistics in parentheses")
    };
    int opt;

    opt = radio_dialog(_("model table options"), NULL, opts, 2, use_tstats, 0);

    if (opt >= 0) {
	use_tstats = opt;
    }
}

