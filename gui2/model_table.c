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
static void rtf_print_model_table (void);
static void print_rtf_row_spec (PRN *prn, int tall);

#define MAX_TABLE_MODELS 6

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
    { N_("/Edit/Copy all/as _RTF"), NULL, rtf_print_model_table, 0, NULL },
    { NULL, NULL, NULL, 0, NULL }
};

static int real_model_table_list_length (void)
{
    int i, len = 0;

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] != NULL) len++;
    }

    return len;    
}

static int model_table_too_many (void)
{
    if (real_model_table_list_length() == MAX_TABLE_MODELS) {
	errbox(_("Model table is full"));
	return 1;
    }
    return 0;
}

static int model_already_in_table (const MODEL *pmod)
{
    int i;

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] == NULL) continue;
	if (pmod == model_list[i]) return 1;
    }

    return 0;
}

int start_model_table_list (const MODEL *pmod, int add_mode)
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

void remove_from_model_table_list (const MODEL *pmod)
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

int add_to_model_table_list (const MODEL *pmod, int add_mode)
{
    const MODEL **tmp;

    /* NLS models won't work */
    if (pmod->ci == NLS) {
	errbox(_("Sorry, NLS models can't be put in the model table"));
	return 1;
    }	

    /* check that list is really started */
    if (model_list_len == 0) {
	return start_model_table_list(pmod, add_mode);
    }

    /* check that the dependent variable is in common */
    if (pmod->list[1] != (model_list[0])->list[1]) {
	errbox(_("Can't add model to table -- this model has a "
		 "different dependent variable"));
	return 1;
    }

    /* check that model is not already on the list */
    if (model_already_in_table(pmod)) {
	errbox(_("Model is already included in the table"));
	return 0;
    }

    /* check that the model table is not already full */
    if (model_table_too_many()) return 1;

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

void free_model_table_list (void)
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
    int i, j, f = 1;
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
	if (model_list[i] == NULL) continue;
	pmod = model_list[i];
	if (f == 1) {
	    for (j=0; j<=pmod->list[0]; j++) {
		grand_list[j] = pmod->list[j];
	    }
	    f = 0;
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
    int i, ci0 = -1;

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] == NULL) continue;
	if (ci0 == -1) {
	    ci0 = (model_list[i])->ci;
	} else {
	    if ((model_list[i])->ci != ci0) return 0;
	}
    }  

    return ci0;
}

static int common_df (void)
{
    int i, dfn0 = -1, dfd0 = -1;

    for (i=0; i<model_list_len; i++) {
	if (model_list[i] == NULL) continue;
	if (dfn0 == -1) {
	    dfn0 = (model_list[i])->dfn;
	    dfd0 = (model_list[i])->dfd;
	} else {
	    if ((model_list[i])->dfn != dfn0) return 0;
	    if ((model_list[i])->dfd != dfd0) return 0;
	}
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

static const char *get_asts (double pval)
{
    return (pval >= 0.1)? "  " : (pval >= 0.05)? "* " : "**";
}

static void print_model_table_coeffs (PRN *prn)
{
    int i, j, k;
    const MODEL *pmod;
    char tmp[16];
    int tex = (prn->format == GRETL_PRINT_FORMAT_TEX);
    int rtf = (prn->format == GRETL_PRINT_FORMAT_RTF);

    /* loop across all variables that appear in any model */
    for (i=2; i<=grand_list[0]; i++) {
	int v = grand_list[i];
	int f = 1;

	if (tex) {
	    tex_escape(tmp, datainfo->varname[v]);
	    pprintf(prn, "%s ", tmp);
	}
	else if (rtf) {
	    print_rtf_row_spec(prn, 0);
	    pprintf(prn, "\\intbl \\qc %s\\cell ", datainfo->varname[v]);
	} else {
	    pprintf(prn, "%8s ", datainfo->varname[v]);
	}

	/* print the coefficient estimates across a row */
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		double x = screen_zero(pmod->coeff[k-1]);
		double s = screen_zero(pmod->sderr[k-1]);
		double pval;
		char numstr[32];

		if (floateq(s, 0.0)) {
		    if (floateq(x, 0.0)) pval = 1.0;
		    else pval = 0.0001;
		} else {
		    pval = tprob(x / s, pmod->dfd);
		}

		if (!tex) {
		    sprintf(numstr, "%#.4g", x);
		    gretl_fix_exponent(numstr);
		}

		if (tex) {
		    if (x < 0) {
			pprintf(prn, "& $-$%#.4g%s ", fabs(x), get_asts(pval));
		    } else {
			pprintf(prn, "& %#.4g%s ", x, get_asts(pval));
		    }
		} else if (rtf) {
		    pprintf(prn, "\\qc %s%s\\cell ", numstr, get_asts(pval));
		} else {
		    pprintf(prn, "%*s%s", (f == 1)? 12 : 10,
			    numstr, get_asts(pval));
		}
		f = 0;
	    } else {
		if (tex) pputs(prn, "& ");
		else if (rtf) pputs(prn, "\\qc \\cell ");
		else pputs(prn, "            ");
	    }
	}

	/* terminate the coefficient row and start the next one */
	if (tex) {
	    pputs(prn, "\\\\\n");
	} else if (rtf) {
	    pputs(prn, "\\intbl \\row\n");
	    print_rtf_row_spec(prn, 1);
	    pputs(prn, "\\intbl ");
	} else {
	    pputs(prn, "\n          ");
	}

	/* print the estimated standard errors across a row */
	f = 1;
	for (j=0; j<model_list_len; j++) {
	    pmod = model_list[j];
	    if (pmod == NULL) continue;
	    if ((k = var_is_in_model(v, pmod))) {
		if (tex) {
		    pprintf(prn, "& (%#.4g) ", pmod->sderr[k-1]);
		} else {
		    char numstr[32];

		    sprintf(numstr, "%#.4g", pmod->sderr[k-1]);
		    gretl_fix_exponent(numstr);
		    if (rtf) {
			if (f == 1) pputs(prn, "\\qc \\cell ");
			pprintf(prn, "\\qc (%s)\\cell ", numstr);
			f = 0;
		    } else {
			sprintf(tmp, "(%s)", numstr);
			pprintf(prn, "%12s", tmp);
		    }
		}
	    } else {
		if (tex) pputs(prn, "& ");
		else if (rtf) pputs(prn, "\\qc \\cell ");
		else pputs(prn, "            ");
	    }
	}
	if (tex) pputs(prn, "\\\\ [4pt] \n");
	else if (rtf) pputs(prn, "\\intbl \\row\n");
	else pputs(prn, "\n\n");
    }
}

static void print_n_r_squared (PRN *prn, int *binary)
{
    int j;
    int same_df;
    const MODEL *pmod;
    int tex = (prn->format == GRETL_PRINT_FORMAT_TEX);
    int rtf = (prn->format == GRETL_PRINT_FORMAT_RTF);

    if (rtf) print_rtf_row_spec(prn, 0);

    if (tex) pprintf(prn, "$%s$ ", _("n"));
    else if (rtf) pprintf(prn, "\\intbl \\qc %s\\cell ", _("n"));
    else pprintf(prn, "%8s ", _("n"));

    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	if (tex) pprintf(prn, "& %d ", pmod->nobs);
	else if (rtf) pprintf(prn, "\\qc %d\\cell ", pmod->nobs);
	else pprintf(prn, "%12d", pmod->nobs);
    }

    if (tex) pputs(prn, "\\\\\n");
    else if (rtf) pputs(prn, "\\intbl \\row\n\\intbl ");
    else pputs(prn, "\n");

    same_df = common_df();
    if (tex) {
	pputs(prn, (same_df)? "$R^2$" : "$\\bar R^2$ ");
    } else if (rtf) {
	pprintf(prn, "\\qc %s\\cell ", 
		(same_df)? "R{\\super 2}" : _("Adj. R{\\super 2}"));
    } else {
	pprintf(prn, "%9s", (same_df)? _("R-squared") : _("Adj. R**2"));
    }

    for (j=0; j<model_list_len; j++) {
	pmod = model_list[j];
	if (pmod == NULL) continue;
	if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
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
    if (rtf) pputs(prn, "\\intbl \\row\n");
    else pputs(prn, "\n\n");
}

int display_model_table (void)
{
    int j, ci;
    int binary = 0;
    int winwidth = 78;
    PRN *prn;

    if (model_list_empty()) {
	errbox(_("The model table is empty"));
	return 1;
    }

    if (make_grand_varlist()) return 1;

    if (bufopen(&prn)) {
	free_model_table_list();
	return 1;
    }

    ci = common_estimator();

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, _("%s estimates"), 
		_(estimator_string(ci, prn->format)));
	pputs(prn, "\n");
    }

    pprintf(prn, _("Dependent variable: %s\n"),
	    datainfo->varname[grand_list[1]]);

    pputs(prn, "\n            ");
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

    print_model_table_coeffs(prn);
    print_n_r_squared(prn, &binary);

    pprintf(prn, "%s\n", _("Standard errors in parentheses"));
    pprintf(prn, "%s\n", _("* indicates significance at the 10 percent level"));
    pprintf(prn, "%s\n", _("** indicates significance at the 5 percent level"));
   
    if (binary) {
	pprintf(prn, "%s\n", _("For logit and probit, R-squared is "
				 "McFadden's pseudo-R-squared"));
    }

    if (real_model_table_list_length() > 5) winwidth = 90;

    view_buffer(prn, winwidth, 450, _("gretl: model table"), PRINT, 
		model_table_items);

    return 0;
}

static void tex_print_model_table (void)
{
    int j, ci;
    int binary = 0;
    char tmp[16];
    PRN *prn;

    if (model_list_empty()) {
	errbox(_("The model table is empty"));
	return;
    }

    if (make_grand_varlist()) return;

    if (bufopen(&prn)) return;
    prn->format = GRETL_PRINT_FORMAT_TEX;

    ci = common_estimator();

    pputs(prn, "\\begin{center}\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, I_("%s estimates"), 
		I_(estimator_string(ci, prn->format)));
	pputs(prn, "\\\\\n");
    }

    tex_escape(tmp, datainfo->varname[grand_list[1]]);
    pprintf(prn, "%s: %s \\\\\n", I_("Dependent variable"), tmp);

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

    print_model_table_coeffs(prn);
    print_n_r_squared(prn, &binary);

    pputs(prn, "\\end{tabular}\n\n");
    pputs(prn, "\\vspace{1em}\n");

    pprintf(prn, "%s\\\\\n", I_("Standard errors in parentheses"));
    pprintf(prn, "{}%s\\\\\n", 
	    I_("* indicates significance at the 10 percent level"));
    pprintf(prn, "{}%s\\\\\n", 
	    I_("** indicates significance at the 5 percent level"));

    if (binary) {
	pprintf(prn, "%s\\\\\n", I_("For logit and probit, $R^2$ is "
				    "McFadden's pseudo-$R^2$"));
    }

    pputs(prn, "\\end{center}\n");

    prn_to_clipboard(prn, COPY_LATEX);
}

static void print_rtf_row_spec (PRN *prn, int tall)
{
    int i, cols = 1 + real_model_table_list_length();
    int col1 = 1000;
    int ht = (tall)? 362 : 262;

    pprintf(prn, "\\trowd \\trqc \\trgaph30\\trleft-30\\trrh%d", ht);
    for (i=0; i<cols; i++) {
	pprintf(prn, "\\cellx%d", col1 +  i * 1400);
    }
    pputs(prn, "\n");
}

static void rtf_print_model_table (void)
{
    int j, ci;
    int binary = 0;
    PRN *prn;

    if (model_list_empty()) {
	errbox(_("The model table is empty"));
	return;
    }

    if (make_grand_varlist()) return;

    if (bufopen(&prn)) return;
    prn->format = GRETL_PRINT_FORMAT_RTF;

    ci = common_estimator();

    pputs(prn, "{\\rtf1\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pputs(prn, "\\par \\qc ");
	pprintf(prn, I_("%s estimates"), 
		I_(estimator_string(ci, prn->format)));
	pputs(prn, "\n");
    }

    pprintf(prn, "\\par \\qc %s: %s\n\\par\n\\par\n{", 
	    I_("Dependent variable"),
	    datainfo->varname[grand_list[1]]);

    /* RTF row stuff */
    print_rtf_row_spec(prn, 1);

    pputs(prn, "\\intbl \\qc \\cell ");
    for (j=0; j<model_list_len; j++) {
	char modhd[16];

	if (model_list[j] == NULL) continue;
	sprintf(modhd, I_("Model %d"), (model_list[j])->ID);
	pprintf(prn, "\\qc %s\\cell ", modhd);
    }
    pputs(prn, "\\intbl \\row\n");
    
    if (ci == 0) {
	char est[12];

	pputs(prn, "\\intbl \\qc \\cell ");

	for (j=0; j<model_list_len; j++) {
	    if (model_list[j] == NULL) continue;
	    strcpy(est, 
		   I_(short_estimator_string((model_list[j])->ci,
					    prn->format)));
	    pprintf(prn, "\\qc %s\\cell ", est);
	}
	pputs(prn, "\\intbl \\row\n");
    }

    print_model_table_coeffs(prn);
    print_n_r_squared(prn, &binary);

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

    prn_to_clipboard(prn, COPY_RTF);
}
