/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* model_table.c for gretl */

#include "gretl.h"
#include "model_table.h"
#include "session.h"
#include "textutil.h"
#include "textbuf.h"
#include "texprint.h"

static MODEL **table_models;
static int n_models;
static char **pnames;
static int n_params;
static int depvarnum;

/* may be governed by, and saved in, rc file */
static int colheads;
static int use_tstats;
static int do_pvals;
static int do_asts = 1;
static int mt_figs = 4;   /* figures for printing */
static char mt_fmt = 'g'; /* floating-point format ('g' or 'f') */

static void print_rtf_row_spec (PRN *prn, int tall);

#define MAX_PORTRAIT_MODELS 6
#define MAX_TABLE_MODELS 12

#define MT_DEBUG 0

enum {
    COLHEAD_ARABIC,
    COLHEAD_ROMAN,
    COLHEAD_ALPHA,
    COLHEAD_NAMES
};

/* transcribe file-scope values to settings.c */

static void record_mtable_settings (void)
{
    set_model_table_prefs(colheads,
			  use_tstats,
			  do_pvals,
			  do_asts,
			  mt_figs,
			  mt_fmt);
}

/* load file-scope values from settings.c */

static void load_mtable_settings (void)
{
    get_model_table_prefs(&colheads,
			  &use_tstats,
			  &do_pvals,
			  &do_asts,
			  &mt_figs,
			  &mt_fmt);
}

static int mtable_opt_from_bundle (gretl_bundle *b,
				   const char *key, int j)
{
    int limits[6][2] = {
	{1,4}, {0,1}, {0,1}, {0,1}, {2,6}, {2,6}
    };
    int *targets[] = {
	&colheads, &use_tstats, &do_pvals, &do_asts
    };
    int val, err = 0;

    val = gretl_bundle_get_int(b, key, &err);
    if (!err) {
	if (val < limits[j][0] || val > limits[j][1]) {
	    err = E_INVARG;
	}
    }
    if (!err) {
	if (j == 0) {
	    /* convert to 0-based */
	    *targets[j] = val - 1;
	} else if (j < 4) {
	    *targets[j] = val;
	} else {
	    mt_figs = val;
	    mt_fmt = (j == 4)? 'g' : 'f';
	}
    }

    return err;
}

static void set_mtable_defaults (void)
{
    colheads = 0;
    use_tstats = 0;
    do_pvals = 0;
    do_asts = 1;
    mt_figs = 4;
    mt_fmt = 'g';
    record_mtable_settings();
}

static int mtable_read_bundle (gretl_bundle *b)
{
    const char *ok_keys[] = {
	"colheads", "tstats", "pvalues", "asterisks",
	"digits", "decplaces", NULL
    };
    gretl_array *a;
    int err = 0;

    a = gretl_bundle_get_keys(b, &err);

    if (!err) {
	const char *key;
	int i, j, n, found;
	int got_digits = 0;

	n = gretl_array_get_length(a);
	for (i=0; i<n && !err; i++) {
	    key = gretl_array_get_data(a, i);
	    found = 0;
	    for (j=0; ok_keys[j] != NULL; j++) {
		if (!strcmp(key, ok_keys[j])) {
		    found = 1;
		    if (j == 4) {
			got_digits = 1;
		    } else if (got_digits && j == 5) {
			err = E_INVARG;
		    }
		    if (!err) {
			err = mtable_opt_from_bundle(b, key, j);
		    }
		}
	    }
	    if (!err && !found) {
		gretl_errmsg_sprintf("modeltab: invalid option-key '%s'", key);
		err = E_INVARG;
	    }
	}
    }

    if (!err) {
	record_mtable_settings();
    }

    gretl_array_destroy(a);

    return err;
}

static void mtable_errmsg (char *msg, int gui)
{
    if (gui) {
	errbox(msg);
    } else {
	gretl_errmsg_set(msg);
    }

    fprintf(stderr, "%s\n", msg);
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
    int i, ID;

    if (pmod == NULL || n_models == 0) {
	return 0;
    }

    for (i=0; i<n_models; i++) {
	if (table_models[i] == NULL) {
	    continue;
	}
	ID = table_models[i]->ID;
	if (pmod == table_models[i] || (pmod->ID > 0 && pmod->ID == ID)) {
	    return 1;
	}
    }

    return 0;
}

GList *table_model_list (void)
{
    GList *list = NULL;

    if (n_models > 0) {
	int i;

	for (i=0; i<n_models; i++) {
	    if (table_models[i] != NULL) {
		list = g_list_append(list, table_models[i]);
	    }
	}
    }

    return list;
}

int model_table_n_models (void)
{
    return n_models;
}

int model_table_landscape (void)
{
    return n_models > MAX_PORTRAIT_MODELS;
}

MODEL *model_table_model_by_index (int i)
{
    if (i >= 0 && i < n_models) {
	return table_models[i];
    } else {
	return NULL;
    }
}

int model_table_position (const MODEL *pmod)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (pmod == table_models[i]) {
	    return i + 1;
	}
    }

    return 0;
}

void clear_model_table (int on_exit, PRN *prn)
{
    int i;

    if (!on_exit && n_models > 0) {
	mark_session_changed();
    }

    for (i=0; i<n_models; i++) {
	/* reduce refcount on the model pointer */
	if (table_models[i] != NULL) {
	    gretl_object_unref(table_models[i], GRETL_OBJ_EQN);
	}
    }

    free(table_models);
    table_models = NULL;

    strings_array_free(pnames, n_params);
    pnames = NULL;
    n_params = 0;

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
	if (table_models[i] != NULL &&
	    table_models[i]->list != NULL &&
	    table_models[i]->list[0] > 0) {
	    return table_models[i]->list[1];
	}
    }

    return -1;
}

static int model_table_precheck (MODEL *pmod, int add_mode)
{
    int gui = (add_mode != MODEL_ADD_BY_CMD);

    if (pmod == NULL) {
	return 1;
    }

    /* various sorts of models that will not work */
    if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM ||
	pmod->ci == DPANEL || pmod->ci == INTREG || pmod->ci == BIPROBIT) {
	mtable_errmsg(_("Sorry, this model can't be put in the model table"),
		      gui);
	return 1;
    }

    /* nor will ARMA, GARCH */
    if (pmod->ci == ARMA || pmod->ci == GARCH) {
	mtable_errmsg(_("Sorry, this model can't be put in the model table"),
		      gui);
	return 1;
    }

    if (n_models > 0) {
	int dv = model_table_depvar();

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
    }

    return 0;
}

/* @pos will usually be 0, which means: add to end of model table array.
   But when reconstituting a session, @pos may be a 1-based index of
   the position within the array that this model should occupy.
*/

static int real_add_to_model_table (MODEL *pmod, int add_mode,
				    int pos, PRN *prn)
{
    int i, n = (pos == 0)? n_models + 1 : pos;

    /* is the list started or not? */
    if (n_models == 0) {
	table_models = mymalloc(n * sizeof *table_models);
	if (table_models == NULL) {
	    return 1;
	}
	for (i=0; i<n; i++) {
	    table_models[i] = NULL;
	}
	n_models = n;
    } else if (pos == 0 || pos > n_models) {
	MODEL **mods;

	mods = myrealloc(table_models, n * sizeof *mods);
	if (mods == NULL) {
	    clear_model_table(0, NULL);
	    return 1;
	}

	for (i=n_models; i<n; i++) {
	    mods[i] = NULL;
	}

	table_models = mods;
	n_models = n;
    }

    table_models[n-1] = pmod;

    /* augment refcount so model won't get deleted */
    gretl_object_ref(pmod, GRETL_OBJ_EQN);

    if (add_mode == MODEL_ADD_FROM_MENU) {
	infobox(_("Model added to table"));
    } else if (add_mode == MODEL_ADD_BY_CMD) {
	pputs(prn, _("Model added to table"));
	pputc(prn, '\n');
    }

    mark_session_changed();

    return 0;
}

int add_to_model_table (MODEL *pmod, int add_mode, int pos, PRN *prn)
{
    if (model_table_precheck(pmod, add_mode)) {
	fprintf(stderr, "add_to_model_table: precheck failed\n");
	return 1;
    }

    return real_add_to_model_table(pmod, add_mode, pos, prn);
}

void remove_from_model_table (MODEL *pmod)
{
    if (n_models > 0) {
	int pos = model_table_position(pmod);

	if (pos == 0) {
	    /* not present */
	    return;
	} else if (n_models == 1) {
	    /* it's the only model in the table */
	    clear_model_table(0, NULL);
	} else {
	    /* remove and reshuffle */
	    int i;

	    gretl_object_unref(pmod, GRETL_OBJ_EQN);

	    for (i=pos-1; i<n_models-1; i++) {
		table_models[i] = table_models[i+1];
	    }

	    table_models[n_models-1] = NULL;
	    n_models--;
	}
    }
}

static int on_param_list (const char *pname)
{
    int i;

    for (i=0; i<n_params; i++) {
	if (!strcmp(pname, pnames[i])) {
	    return 1;
	}
    }

    return 0;
}

static int add_to_param_list (const MODEL *pmod)
{
    char pname[VNAMELEN];
    int i, err = 0;

    for (i=0; i<pmod->ncoeff && !err; i++) {
	gretl_model_get_param_name(pmod, dataset, i, pname);
	if (!on_param_list(pname)) {
	    err = strings_array_add(&pnames, &n_params, pname);
	}
    }

    return err;
}

static int make_full_param_list (void)
{
    const MODEL *pmod;
    int first = 1;
    int i, err = 0;

    strings_array_free(pnames, n_params);
    pnames = NULL;
    n_params = 0;

    for (i=0; i<n_models && !err; i++) {
	if (table_models[i] != NULL) {
	    pmod = table_models[i];
	    if (first) {
		depvarnum = gretl_model_get_depvar(pmod);
		first = 0;
	    }
	    err = add_to_param_list(pmod);
	}
    }

    return err;
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

static const char *short_estimator_string (const MODEL *pmod, PRN *prn)
{
    if (pmod->ci == HSK) {
	return N_("HSK");
    } else if (pmod->ci == ARCH) {
	return N_("ARCH");
    } else if (pmod->ci == WLS) {
	if (gretl_model_get_int(pmod, "iters")) {
	    return N_("MLE");
	} else {
	    return N_("WLS");
	}
    } else if (pmod->ci == PANEL) {
	if (pmod->opt & OPT_F) {
	    return N_("Within");
	} else if (pmod->opt & OPT_U) {
	    return N_("GLS");
	} else {
	    return N_("Between");
	}
    } else if (pmod->ci == AR1) {
	if (pmod->opt & OPT_H) {
	    return N_("HILU");
	} else if (pmod->opt & OPT_P) {
	    return N_("PWE");
	} else {
	    return N_("CORC");
	}
    } else {
	return estimator_string(pmod, prn);
    }
}

static const char *get_asts (double pval, int tex)
{
    if (pval < 0.01) {
	return tex ? "$^{***}$" : "***";
    } else if (pval < 0.05) {
	return tex ? "$^{**}$" : "**";
    } else if (pval < 0.10) {
	return tex ? "$^{*}$" : "*";
    } else {
	return tex ? "" : "   ";
    }
}

static const char *get_pre_asts (double pval)
{
    if (pval < 0.01) {
	return "$\\,\\,\\,$";
    } else if (pval < 0.05) {
	return "$\\,\\,$";
    } else if (pval < 0.10) {
	return "$\\,$";
    } else {
	return "";
    }
}

static void terminate_coeff_row (int namewidth, PRN *prn)
{
    if (tex_format(prn)) {
	pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
	pputs(prn, "\\intbl \\row\n");
	print_rtf_row_spec(prn, 1);
	pputs(prn, "\\intbl ");
    } else {
	pputc(prn, '\n');
	bufspace(namewidth + 2, prn);
    }
}

static double modtab_get_pval (const MODEL *pmod, int k)
{
    double x = pmod->coeff[k];
    double s = pmod->sderr[k];
    double pval = NADBL;

    if (!na(x) && !na(s)) {
	pval = coeff_pval(pmod->ci, x / s, pmod->dfd);
    }

    return pval;
}

static void mt_print_value (char *s, double x)
{
    if (mt_fmt == 'f') {
	sprintf(s, "%.*f", mt_figs, x);
    } else {
	sprintf(s, "%#.*g", mt_figs, x);
	gretl_fix_exponent(s);
    }
}

static void print_model_table_coeffs (int namewidth, int colwidth, PRN *prn)
{
    const MODEL *pmod;
    char numstr[32], tmp[64];
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    int i, j, k;

    /* loop across all variables that appear in any model */

    for (i=0; i<n_params; i++) {
	char *pname = pnames[i];
	int first_coeff = 1;
	int first_se = 1;
	int first_pval = 1;

	if (tex) {
	    tex_escape(tmp, pname);
	    pprintf(prn, "%s ", tmp);
	} else if (rtf) {
	    print_rtf_row_spec(prn, 0);
	    pprintf(prn, "\\intbl \\qc %s\\cell ", pname);
	} else if (strlen(pnames[i]) > namewidth) {
	    sprintf(tmp, "%.*s...", namewidth - 3, pname);
	    pprintf(prn, "%-*s ", namewidth, tmp);
	} else {
	    pprintf(prn, "%-*s ", namewidth, pname);
	}

	/* print the coefficient estimates across a row */

	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) {
		continue;
	    }
	    if ((k = gretl_model_get_param_number(pmod, dataset, pname)) >= 0) {
		double x = screen_zero(pmod->coeff[k]);

		mt_print_value(numstr, x);

		if (do_asts) {
		    double pval = modtab_get_pval(pmod, k);

		    if (tex) {
			if (x < 0) {
			    pprintf(prn, "& %s$-$%s%s ", get_pre_asts(pval),
				    numstr + 1, get_asts(pval, 1));
			} else {
			    pprintf(prn, "& %s%s%s ", get_pre_asts(pval),
				    numstr, get_asts(pval, 1));
			}
		    } else if (rtf) {
			pprintf(prn, "\\qc %s%s\\cell ", numstr, get_asts(pval, 0));
		    } else {
			/* note: strlen(asts) = 3 */
			pprintf(prn, "%*s%s", (first_coeff)? colwidth : colwidth - 3,
				numstr, get_asts(pval, 0));
		    }
		} else {
		    /* not showing asterisks */
		    if (tex) {
			if (x < 0) {
			    pprintf(prn, "& $-$%s ", numstr + 1);
			} else {
			    pprintf(prn, "& %s ", numstr);
			}
		    } else if (rtf) {
			pprintf(prn, "\\qc %s\\cell ", numstr);
		    } else {
			pprintf(prn, "%*s", colwidth, numstr);
		    }
		}
		first_coeff = 0;
	    } else {
		/* variable not present in this column */
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    bufspace(colwidth, prn);
		}
	    }
	}

	terminate_coeff_row(namewidth, prn);

	/* print the t-stats or standard errors across a row */

	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) {
		continue;
	    }
	    if ((k = gretl_model_get_param_number(pmod, dataset, pname)) >= 0) {
		double val;

		if (use_tstats) {
		    val = pmod->coeff[k] / pmod->sderr[k];
		} else {
		    val = pmod->sderr[k];
		}

		mt_print_value(numstr, val);

		if (tex) {
		    if (val < 0) {
			pprintf(prn, "& \\subsize{($-$%s)} ", numstr + 1);
		    } else {
			pprintf(prn, "& \\subsize{(%s)} ", numstr);
		    }
		} else if (rtf) {
		    if (first_se) {
			pputs(prn, "\\qc \\cell ");
		    }
		    pprintf(prn, "\\qc (%s)\\cell ", numstr);
		} else {
		    sprintf(tmp, "(%s)", numstr);
		    pprintf(prn, "%*s", colwidth, tmp);
		}
		first_se = 0;
	    } else {
		/* variable not present in this column */
		if (tex) {
		    pputs(prn, "& ");
		} else if (rtf) {
		    pputs(prn, "\\qc \\cell ");
		} else {
		    bufspace(colwidth, prn);
		}
	    }
	}

	if (do_pvals) {
	    terminate_coeff_row(namewidth, prn);
	    for (j=0; j<n_models; j++) {
		pmod = table_models[j];
		if (pmod == NULL) {
		    continue;
		}
		if ((k = gretl_model_get_param_number(pmod, dataset, pname)) >= 0) {
		    double pval = modtab_get_pval(pmod, k);

		    if (na(pval)) {
			strcpy(numstr, "NA");
		    } else {
			sprintf(numstr, "%.*f", mt_figs, pval);
		    }

		    if (tex) {
			pprintf(prn, "& \\subsize{[%s]} ", numstr);
		    } else if (rtf) {
			if (first_pval) {
			    pputs(prn, "\\qc \\cell ");
			}
			pprintf(prn, "\\qc [%s]\\cell ", numstr);
		    } else {
			sprintf(tmp, "[%s]", numstr);
			pprintf(prn, "%*s", colwidth, tmp);
		    }
		    first_pval = 0;
		} else {
		    /* variable not present in this column */
		    if (tex) {
			pputs(prn, "& ");
		    } else if (rtf) {
			pputs(prn, "\\qc \\cell ");
		    } else {
			bufspace(colwidth, prn);
		    }
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

enum {
    MT_LNL,
    MT_RSQ,
    MT_FSTAT
};

static int any_stat (int s)
{
    int i;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    if (s == MT_LNL && !na(table_models[i]->lnL)) {
		return 1;
	    } else if (s == MT_RSQ && !na(table_models[i]->rsq)) {
		return 1;
	    } else if (s == MT_FSTAT && !na(table_models[i]->fstt)) {
		return 1;
	    }
	}
    }

    return 0;
}

static int catch_bad_point (char *s)
{
    int len = strlen(s);
    int c = s[len-1];

    if (c == '.' || c == ',') {
	s[len-1] = '\0';
	return 1;
    } else {
	return 0;
    }
}

static void print_equation_stats (int width0, int colwidth, PRN *prn,
				  int *binary)
{
    const MODEL *pmod;
    int same_df, any_R2, any_ll;
    int tex = tex_format(prn);
    int rtf = rtf_format(prn);
    double rsq;
    int j;

    if (rtf) {
	print_rtf_row_spec(prn, 0);
    }

    if (tex) {
	pprintf(prn, "$%s$ ", _("n"));
    } else if (rtf) {
	pprintf(prn, "\\intbl \\qc %s\\cell ", _("n"));
    } else {
	pprintf(prn, "%*s", width0, _("n"));
    }

    for (j=0; j<n_models; j++) {
	pmod = table_models[j];
	if (pmod != NULL) {
	    if (tex) {
		pprintf(prn, "& %d ", pmod->nobs);
	    } else if (rtf) {
		pprintf(prn, "\\qc %d\\cell ", pmod->nobs);
	    } else {
		pprintf(prn, "%*d", colwidth, pmod->nobs);
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

    any_R2 = any_stat(MT_RSQ);
    any_ll = any_stat(MT_LNL);

    if (any_R2) {
	/* print R^2 values */
	if (tex) {
	    pputs(prn, same_df ? "$R^2$" : "$\\bar R^2$ ");
	} else if (rtf) {
	    pprintf(prn, "\\qc %s\\cell ",
		    (same_df)? "R{\\super 2}" : _("Adj. R{\\super 2}"));
	} else {
	    pprintf(prn, "%*s", width0, same_df ? _("R-squared") :
		    _("Adj. R^2"));
	}

	for (j=0; j<n_models; j++) {
	    pmod = table_models[j];
	    if (pmod == NULL) continue;
	    if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
		rsq = pmod->rsq;
	    } else {
		rsq = (same_df)? pmod->rsq : pmod->adjrsq;
	    }
	    if (na(rsq)) {
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
		    pprintf(prn, "& %.*f ", mt_figs, pmod->rsq);
		} else if (rtf) {
		    pprintf(prn, "\\qc %.*f\\cell ", mt_figs, pmod->rsq);
		} else {
		    pprintf(prn, "%*.*f", colwidth, mt_figs, pmod->rsq);
		}
	    } else {
		if (tex) {
		    pprintf(prn, "& %.*f ", mt_figs, rsq);
		} else if (rtf) {
		    pprintf(prn, "\\qc %.*f\\cell ", mt_figs, rsq);
		} else {
		    pprintf(prn, "%*.*f", colwidth, mt_figs, rsq);
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
	    pprintf(prn, "%*s", width0, "lnL");
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
		gchar *numstr = NULL;

		if (tex) {
		    if (pmod->lnL > 0) {
			numstr = g_strdup_printf("%#.*g", mt_figs, pmod->lnL);
			catch_bad_point(numstr);
			pprintf(prn, "& %s ", numstr);
		    } else {
			numstr = g_strdup_printf("%#.*g", mt_figs, -pmod->lnL);
			catch_bad_point(numstr);
			pprintf(prn, "& $-$%s ", numstr);
		    }
		} else if (rtf) {
		    numstr = g_strdup_printf("%#.*g", mt_figs, pmod->lnL);
		    catch_bad_point(numstr);
		    pprintf(prn, "\\qc %s\\cell ", numstr);
		} else {
		    /* plain text */
		    numstr = g_strdup_printf("%#*.*g", colwidth, mt_figs, pmod->lnL);
		    if (catch_bad_point(numstr)) {
			g_free(numstr);
			numstr = g_strdup_printf("%#*.*g", colwidth + 1, mt_figs, pmod->lnL);
			catch_bad_point(numstr);
		    }
		    pputs(prn, numstr);
		}
		g_free(numstr);
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

static int mtab_max_namelen (void)
{
    int i, len, maxlen = 8;

    for (i=0; i<n_params; i++) {
	len = strlen(pnames[i]);
	if (len > maxlen) {
	    maxlen = len;
	}
    }

    if (maxlen > 15) {
	maxlen = 15;
    }

    return maxlen;
}

static int mtab_get_colwidth (void)
{
    int maxlen = 0;
    int cw = 13;

    if (colheads == COLHEAD_NAMES) {
	const MODEL *pmod;
	int i, len;

	for (i=0; i<n_models; i++) {
	    pmod = table_models[i];
	    if (pmod != NULL) {
		if (pmod->name != NULL) {
		    len = strlen(pmod->name);
		} else {
		    gchar *tmp = g_strdup_printf(_("Model %d"), pmod->ID);

		    len = strlen(tmp);
		    g_free(tmp);
		}
		if (len > maxlen) {
		    maxlen = (len > 31)? 31 : len;
		}
	    }
	}
    }

    if (maxlen < 11 && mt_figs > 0 && mt_figs < 4) {
	cw -= 4 - mt_figs;
    }

    cw = (cw < maxlen + 2)? maxlen + 2 : cw;

    return cw;
}

static void print_estimator_strings (int colwidth, PRN *prn)
{
    const char *s;
    char est[32];
    int i;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    s = short_estimator_string(table_models[i], prn);
	    if (tex_format(prn)) {
		strcpy(est, _(s));
		pprintf(prn, " & %s ", est);
	    } else if (rtf_format(prn)) {
		strcpy(est, _(s));
		pprintf(prn, "\\qc %s\\cell ", est);
	    } else {
		strcpy(est, _(s));
		print_centered(est, colwidth, prn);
	    }
	}
    }
}

static void print_model_head (const MODEL *pmod, int j, int colwidth,
			      PRN *prn)
{
    gchar *targ = NULL;

    if (colheads == COLHEAD_ARABIC) {
	targ = g_strdup_printf("(%d)", j + 1);
    } else if (colheads == COLHEAD_ROMAN) {
	const char *R[] = {
	    "I", "II", "III", "IV", "V", "VI",
	    "VII", "VIII", "IX", "X", "XI", "XII"
	};

	targ = g_strdup_printf("%s", R[j]);
    } else if (colheads == COLHEAD_ALPHA) {
	targ = g_strdup_printf("%c", 'A' + j);
    } else if (tex_format(prn)) {
	if (pmod->name != NULL) {
	    targ = tex_escape_new(pmod->name);
	} else {
	    targ = g_strdup_printf(_("Model %d"), pmod->ID);
	}
    } else if (rtf_format(prn)) {
	if (pmod->name != NULL) {
	    targ = g_strdup(pmod->name);
	    gretl_utf8_truncate(targ, 31);
	} else {
	    targ = g_strdup_printf(_("Model %d"), pmod->ID);
	}
    } else {
	if (pmod->name != NULL) {
	    targ = g_strdup(pmod->name);
	    gretl_utf8_truncate(targ, 31);
	} else {
	    targ = g_strdup_printf(_("Model %d"), pmod->ID);
	}
    }

    if (tex_format(prn)) {
	pprintf(prn, " & %s ", targ);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\qc %s\\cell ", targ);
    } else {
	print_centered(targ, colwidth, prn);
    }

    g_free(targ);
}

static void print_column_heads (int colwidth, PRN *prn)
{
    int i, j = 0;

    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    print_model_head(table_models[i], j++, colwidth, prn);
	}
    }
}

static const char *sigstrs[] = {
    N_("significant at the 10 percent level"),
    N_("significant at the 5 percent level"),
    N_("significant at the 1 percent level")
};

static void plain_print_model_table (PRN *prn)
{
    int namelen = mtab_max_namelen();
    int colwidth = mtab_get_colwidth();
    int ci = common_estimator();
    int binary = 0;

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, _("%s estimates"),
		_(estimator_string(table_models[0], prn)));
	pputc(prn, '\n');
    }

    pprintf(prn, _("Dependent variable: %s\n"), dataset->varname[depvarnum]);

    pputc(prn, '\n');
    bufspace(namelen + 4, prn);
    print_column_heads(colwidth, prn);
    pputc(prn, '\n');

    if (ci == 0) {
	bufspace(namelen + 4, prn);
	print_estimator_strings(colwidth, prn);
	pputc(prn, '\n');
    }

    pputc(prn, '\n');

    print_model_table_coeffs(namelen, colwidth, prn);
    print_equation_stats(namelen + 1, colwidth, prn, &binary);

    if (use_tstats) {
	pprintf(prn, "%s\n", _("t-statistics in parentheses"));
    } else {
	pprintf(prn, "%s\n", _("Standard errors in parentheses"));
    }

    if (do_pvals) {
	pprintf(prn, "%s\n", _("p-values in brackets"));
    }

    if (do_asts) {
	pprintf(prn, "*   %s\n", _(sigstrs[0]));
	pprintf(prn, "**  %s\n", _(sigstrs[1]));
	pprintf(prn, "*** %s\n", _(sigstrs[2]));
    }

    if (binary) {
	pprintf(prn, "%s\n", _("For logit and probit, R-squared is "
			       "McFadden's pseudo-R-squared"));
    }
}

int display_model_table (int gui)
{
    int winwidth;
    PRN *prn;

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), gui);
	return 1;
    } else if (make_full_param_list()) {
	return 1;
    } else if (bufopen(&prn)) {
	clear_model_table(0, NULL);
	return 1;
    }

    load_mtable_settings();
    plain_print_model_table(prn);
    winwidth = real_table_n_models() > 5 ? 90 : 78;
    view_buffer(prn, winwidth, 450, _("gretl: model table"),
		VIEW_MODELTABLE, NULL);

    return 0;
}

static int tex_print_model_table (PRN *prn)
{
    int binary = 0;
    char tmp[32];
    int i, ci;

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), 1);
	return 1;
    }

    if (make_full_param_list()) {
	return 1;
    }

    if (tex_doc_format(prn)) {
	gretl_tex_preamble(prn, GRETL_FORMAT_MODELTAB);
    }

    pputs(prn, "\n\\newcommand{\\subsize}[1]{\\footnotesize{#1}}\n\n");
    pputs(prn, "\\begin{center}\n");

    ci = common_estimator();

    if (ci > 0) {
	/* all models use same estimation procedure */
	pprintf(prn, _("%s estimates"),
		_(estimator_string(table_models[0], prn)));
	pputs(prn, "\\\\\n");
    }

    tex_escape(tmp, dataset->varname[depvarnum]);
    pprintf(prn, "%s: %s \\\\\n", _("Dependent variable"), tmp);

    pputs(prn, "\\vspace{1em}\n\n");
    pputs(prn, "\\begin{longtable}{l");
    for (i=0; i<n_models; i++) {
	if (table_models[i] != NULL) {
	    pputc(prn, 'c');
	}
    }
    pputs(prn, "}\n");

    print_column_heads(0, prn);
    pputs(prn, "\\\\ ");

    if (ci == 0) {
	pputc(prn, '\n');
	print_estimator_strings(0, prn);
	pputs(prn, "\\\\ ");
    }

    pputs(prn, " [6pt] \n");

    print_model_table_coeffs(0, 0, prn);
    print_equation_stats(0, 0, prn, &binary);

    pputs(prn, "\\end{longtable}\n\n");
    pputs(prn, "\\vspace{1em}\n");

    if (use_tstats) {
	pprintf(prn, "%s\\\\\n", _("$t$-statistics in parentheses"));
    } else {
	pprintf(prn, "%s\\\\\n", _("Standard errors in parentheses"));
    }

    if (do_pvals) {
	pprintf(prn, "%s\\\\\n", _("$p$-values in brackets"));
    }

    if (do_asts) {
	pprintf(prn, "{}* %s\\\\\n", _(sigstrs[0]));
	pprintf(prn, "{}** %s\\\\\n", _(sigstrs[1]));
	pprintf(prn, "{}*** %s\\\\\n", _(sigstrs[2]));
    }

    if (binary) {
	pprintf(prn, "%s\\\\\n", _("For logit and probit, $R^2$ is "
				   "McFadden's pseudo-$R^2$"));
    }

    pputs(prn, "\\end{center}\n");

    if (tex_doc_format(prn)) {
	pputs(prn, "\n\\end{document}\n");
    }

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
    int ci, binary = 0;

    if (model_table_is_empty()) {
	mtable_errmsg(_("The model table is empty"), 1);
	return 1;
    }

    if (make_full_param_list()) {
	return 1;
    }

    ci = common_estimator();

    pputs(prn, "{\\rtf1\n");

    if (ci > 0) {
	/* all models use same estimation procedure */
	pputs(prn, "\\par \\qc ");
	pprintf(prn, _("%s estimates"),
		_(estimator_string(table_models[0], prn)));
	pputc(prn, '\n');
    }

    pprintf(prn, "\\par \\qc %s: %s\n\\par\n\\par\n{",
	    _("Dependent variable"), dataset->varname[depvarnum]);

    print_rtf_row_spec(prn, 1);
    pputs(prn, "\\intbl \\qc \\cell ");
    print_column_heads(0, prn);
    pputs(prn, "\\intbl \\row\n");

    if (ci == 0) {
	pputs(prn, "\\intbl \\qc \\cell ");
	print_estimator_strings(0, prn);
	pputs(prn, "\\intbl \\row\n");
    }

    print_model_table_coeffs(0, 0, prn);
    print_equation_stats(0, 0, prn, &binary);

    pputs(prn, "}\n\n");

    if (use_tstats) {
	pprintf(prn, "\\par \\qc %s\n", _("t-statistics in parentheses"));
    } else {
	pprintf(prn, "\\par \\qc %s\n", _("Standard errors in parentheses"));
    }

    if (do_pvals) {
	pprintf(prn, "\\par \\qc %s\n", _("p-values in brackets"));
    }

    if (do_asts) {
	pprintf(prn, "\\par \\qc * %s\n", _(sigstrs[0]));
	pprintf(prn, "\\par \\qc ** %s\n", _(sigstrs[1]));
	pprintf(prn, "\\par \\qc *** %s\n", _(sigstrs[2]));
    }

    if (binary) {
	pprintf(prn, "\\par \\qc %s\n", _("For logit and probit, "
					  "R{\\super 2} is "
					  "McFadden's pseudo-R{\\super 2}"));
    }

    pputs(prn, "\\par\n}\n");

    return 0;
}

int special_print_model_table (PRN *prn)
{
    load_mtable_settings();

    if (tex_format(prn)) {
	return tex_print_model_table(prn);
    } else if (rtf_format(prn)) {
	return rtf_print_model_table(prn);
    } else {
	return 1;
    }
}

static int cli_modeltab_add (PRN *prn)
{
    GretlObjType type;
    void *ptr = get_last_model(&type);
    int err = 0;

#if MT_DEBUG
    fprintf(stderr, "cli_modeltab_add, ptr = %p\n", ptr);
#endif

    if (type != GRETL_OBJ_EQN) {
	gretl_errmsg_set(_("No model is available"));
	err = 1;
    } else {
	MODEL *pmod = (MODEL *) ptr;
	MODEL *cpy = NULL;
	int freeit = 0;

	err = model_table_precheck(pmod, MODEL_ADD_BY_CMD);
	if (err) {
	    fprintf(stderr, "cli_modeltab_add: precheck error\n");
	    return err;
	}

	cpy = get_model_by_ID(pmod->ID);
	if (cpy == NULL) {
	    cpy = gretl_model_copy(pmod);
	    if (cpy == NULL) {
		err = E_ALLOC;
	    } else {
		freeit = 1;
	    }
	}

	if (!err) {
	    err = real_add_to_model_table(cpy, MODEL_ADD_BY_CMD, 0, prn);
#if MT_DEBUG
	    fprintf(stderr, "real_add: cpy = %p, err = %d\n", (void *) cpy, err);
#endif
	}

	if (err && freeit) {
	    gretl_model_free(cpy);
	}
    }

    return err;
}

static int print_model_table_direct (const char *fname,
				     gretlopt opt,
				     PRN *msgprn)
{
    char outfile[MAXLEN];
    PRN *prn;
    int err = 0;

    if (model_table_is_empty()) {
	gretl_errmsg_set(_("The model table is empty"));
	return E_DATA;
    }

    if (make_full_param_list()) {
	return E_DATA;
    }

    strcpy(outfile, fname);
    gretl_maybe_switch_dir(fname);

    prn = gretl_print_new_with_filename(outfile, &err);
    if (err) {
	return err;
    }

    if (has_suffix(fname, ".tex")) {
	gretl_print_set_format(prn, GRETL_FORMAT_TEX);
	if (opt & OPT_C) {
	    gretl_print_toggle_doc_flag(prn);
	}
    } else if (has_suffix(fname, ".rtf")) {
	gretl_print_set_format(prn, GRETL_FORMAT_RTF);
    }

    if (tex_format(prn)) {
	err = tex_print_model_table(prn);
    } else if (rtf_format(prn)) {
	err = rtf_print_model_table(prn);
    } else {
	plain_print_model_table(prn);
    }

    gretl_print_destroy(prn);

    if (!err) {
	pprintf(msgprn, _("wrote %s\n"), outfile);
    }

    return err;
}

int modeltab_exec (const char *param, gretlopt opt, PRN *prn)
{
    int err = 0;

#if MT_DEBUG
    fprintf(stderr, "*** modeltab_exec ***\n");
#endif

    if (param != NULL && (opt & OPT_O)) {
	/* the --output option is incompatible with the
	   various "action" parameters
	*/
	return E_BADOPT;
    }

    if (opt & OPT_B) {
	/* --options */
	const char *optstr;
	gretl_bundle *b = NULL;

	optstr = get_optval_string(MODELTAB, OPT_B);
	if (optstr != NULL && !strcmp(optstr, "defaults")) {
	    set_mtable_defaults();
	} else if ((b = get_bundle_by_name(optstr)) == NULL) {
	    err = E_DATA;
	} else {
	    err = mtable_read_bundle(b);
	}
    }

    if (!err && (opt & OPT_O)) {
	/* --output */
	const char *fname = get_optval_string(MODELTAB, OPT_O);

	if (fname == NULL) {
	    err = E_PARSE;
	} else {
	    err = print_model_table_direct(fname, opt, prn);
	}
	return err;
    }

    if (!err && param != NULL && *param != '\0') {
	if (!strcmp(param, "add")) {
	    err = cli_modeltab_add(prn);
	} else if (!strcmp(param, "show")) {
	    err = display_model_table(0);
	} else if (!strcmp(param, "free")) {
	    if (!model_table_is_empty()) {
		clear_model_table(0, prn);
	    }
	} else {
	    err = E_PARSE;
	}
    }

    if (err) {
	gui_errmsg(err);
    }

    return err;
}

/* Implement on-the-fly reformatting of the model table, in
   response to clicking on the "Reformat" button in a model
   table viewer window.
*/

void format_model_table (windata_t *vwin)
{
    int colhead_opt;
    int se_opt, pv_opt;
    int ast_opt, figs;
    char fmt;
    int resp;

    /* load the stored values from settings.c */
    load_mtable_settings();

    /* and record these to temporary variables */
    colhead_opt = colheads;
    se_opt = use_tstats;
    pv_opt = do_pvals;
    ast_opt = do_asts;
    figs = mt_figs;
    fmt = mt_fmt;

    /* put up a dialog giving formatting choices */
    resp = model_table_dialog(&colhead_opt, &se_opt, &pv_opt, &ast_opt,
			      &figs, &fmt, vwin->main);

    if (resp == GRETL_CANCEL) {
	return;
    }

    if (colhead_opt == colheads && se_opt == use_tstats &&
	pv_opt == do_pvals && ast_opt == do_asts && figs == mt_figs) {
	/* no-op: nothing was changed */
	return;
    } else {
	GtkTextBuffer *buf;
	const char *newtext;
	PRN *prn;

	/* update the active variables */
	colheads = colhead_opt;
	use_tstats = se_opt;
	do_pvals = pv_opt;
	do_asts = ast_opt;
	mt_figs = figs;
	mt_fmt = fmt;

	/* and push their values back to settings.c */
	record_mtable_settings();

	if (bufopen(&prn)) {
	    return;
	}
	/* reprint the table and replace the prior window content */
	plain_print_model_table(prn);
	newtext = gretl_print_get_buffer(prn);
	buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->text));
	gtk_text_buffer_set_text(buf, "", -1);
	textview_set_text(vwin->text, newtext);
	gretl_print_destroy(prn);
    }
}
